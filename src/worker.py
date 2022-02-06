import os
import re
import time
import queue
import pprint
import numpy as np
import torch
from tempfile import TemporaryDirectory
from threading import Thread
import multiprocessing.managers
import multiprocessing as mp
from pyfaidx import Fasta
from cyvcf2 import VCF

from . models import VariantTokenizer, VariantFilterModel, ModelInputStruct
from . predict import VariantToVector
from . cbuf import CircularBuffer
from . base import Site, Genotype

class Worker(mp.Process):
    def __init__(self, manager=None, **kw):
        super().__init__(**kw)
        assert manager
        self.manager = manager
        self._running = mp.Event()

    @property
    def running(self):
        return self._running.is_set()
    
    def shutdown(self):
        self._running.clear()

    def _run(self):
        pass

    def run(self):
        pwd = os.getcwd()
        with TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            self._running.set()
            msg = f"Starting {self.__class__.__name__} worker run={self.running}, tmpdir={tmpdir}"
            print(msg)
            time.sleep(1)
            self._run()
        os.chdir(pwd)
        msg = f"{self.__class__.__name__} shutting down"
        print(msg)

class GatherWorker(Worker):
    def __init__(self, vcf_in_path=None, vcf_idx_path=None, vcf_out_path=None, **kw):
        super().__init__(**kw)
        self.vcf_in_path = vcf_in_path
        self.vcf_idx_path = vcf_idx_path
        self.vcf_out_path = vcf_out_path
        self.in_q = self.manager.to_gather
        self.site_cache = {}
        self.results_cache = {}

    def process_item(self, item=None):
        site_id = item.site_id
        if isinstance(item, Site):
            if site_id not in self.site_cache:
                self.site_cache[site_id] = item
            else:
                site = self.site_cache[site_id]
                site.update(item)
            return
        if isinstance(item, Genotype):
            # XXX: place holder for place holder
            site = self.site_cache.get(site_id)
            site.update_genotype(item)
            return
        raise ValueError(type(item))

    def wait_on(self, site_id=None):
        while True:
            site = self.site_cache.get(site_id)
            if site and not site.is_pending:
                return self.site_cache.pop(site_id)

            try:
                info = self.in_q.get(timeout=1)
            except queue.Empty:
                print(f"{self.__class__.__name__} wait_in(#{site_id}): in_q empty. outstanding={len(self.site_cache)}")
                pending = [x.site_id for x in self.site_cache.values() if x.is_pending]
                print(pending)
                continue
            self.in_q.task_done()
            if type(info) == list:
                for item in info:
                    self.process_item(item)
            else:
                self.process_item(info)

    def _run(self):
        vcf_in = VCF(self.vcf_in_path)
        if self.vcf_idx_path:
            vcf_in.set_index(self.vcf_idx_path)

        for (site_id, site) in enumerate(vcf_in):
            site_info = self.wait_on(site_id=site_id)
            assert site_id == site_info.site_id
            assert site.POS == site_info.pos
            assert site.CHROM == site_info.chrom
            print(site_info)
            if site_id == 5000:
                break
        self.manager.flush_model.set()
        self.manager.flush_vectorizer.set()

class ModelWorker(Worker):
    def __init__(self, model_path=None, batch_size=None, klen=None, **kw):
        super().__init__(**kw)
        self.model_path = model_path
        self.batch_size = batch_size
        self.klen = klen
        self.in_q = self.manager.to_model
        self.out_q = self.manager.to_gather
        self.processed_count = 0

    def get_batch(self, timeout=5):
        # XXX: need to flush somehow
        batch = []
        while (len(batch) < self.batch_size):
            try:
                item = self.in_q.pop(timeout=timeout)
            except TimeoutError:
                print(f"{self.__class__.__name__} get_batch: buf={self.in_q.qsize.value}, processed={self.processed_count}")
                break
            batch.append(item)
        if not batch:
            return (None, None)
        keys = [ds.get_key() for ds in batch]
        batch = np.array([ds.as_numpy() for ds in batch])
        batch = batch.transpose([1, 0, 2])
        batch = torch.tensor(batch)
        headers = ('input_ids', 'attention_mask', 'token_type_ids')
        batch = dict(zip(headers, batch))
        return (keys, batch)

    def _run(self):
        model = VariantFilterModel(model_path=self.model_path, klen=self.klen)

        while self.running:
            (batch_keys, batch) = self.get_batch()
            if batch:
                outp = model.predict(batch)
                result = []
                assert len(outp['log_odds']) == len(batch_keys)
                for (log_odds, keys) in zip(outp['log_odds'], batch_keys):
                    ns = keys.copy()
                    ns['log_odds'] = log_odds.tolist()
                    ns['status'] = 'called'
                    gt = Genotype(**ns)
                    result.append(gt)
                    self.processed_count += 1
                self.out_q.put(result)
            #
            if self.in_q.qsize.value == 0 and self.manager.flush_model.is_set():
                break

class VariantToVectorWorker(Worker):
    def __init__(self, ref_path=None, tokenizer_config=None, window=96, **kw):
        super().__init__(**kw)
        self.ref_path = ref_path
        self.tokenizer_config = tokenizer_config
        self.window = window
        self.in_q = self.manager.to_vectorizer
        self.to_model = self.manager.to_model
        self.to_gather = self.manager.to_gather
        self.sent_to_model = 0

    def dispatch_site(self, site_info=None, gt_inps=None):
        site_id = site_info.site_id
        model_inputs = []
        for (genotype_id, gt_inp) in gt_inps.items():
            if gt_inp is None:
                site_info.genotypes[genotype_id].status = 'skipped'
                continue
            gt_inp["site_id"] = site_info.site_id
            gt_inp["genotype_id"] = genotype_id
            inp = ModelInputStruct(**gt_inp)
            model_inputs.append(inp)
        assert model_inputs or site_info.is_pending == False

        # dispatch
        self.to_gather.put(site_info)
        for inp in model_inputs:
            self.sent_to_model += 1
            self.to_model.push(inp)

    def _run(self):
        ref = Fasta(self.ref_path)
        tokenizer = VariantTokenizer(**self.tokenizer_config)
        vectorizer = VariantToVector(ref=ref, tokenizer=tokenizer, window=self.window)

        self._running.set()
        while self.running:
            try:
                site_info = self.in_q.get(timeout=1)
                self.in_q.task_done()
            except queue.Empty:
                if self.manager.flush_vectorizer.is_set():
                    break
                print(f"{self.__class__.__name__} loop: in_q empty. buf={self.to_model.qsize.value}" + f" sent={self.sent_to_model}")
                continue

            gt_inps = vectorizer.process_site(site_info)
            self.dispatch_site(site_info=site_info, gt_inps=gt_inps)
            
class ScatterWorker(Worker):
    def __init__(self, vcf_in_path=None, vcf_idx_path=None, **kw):
        super().__init__(**kw)
        self.vcf_in_path = vcf_in_path
        self.vcf_idx_path = vcf_idx_path
        self.to_vectorizer_que = self.manager.to_vectorizer
        self.to_gather_que = self.manager.to_gather

    def push_site(self, site=None, site_id=None):
        site_info = Site.load_from_site(site=site, site_id=site_id)
        if not (site.is_snp or site.is_indel):
            site_info.set_site_status("skipped")
            self.to_gather_que.put(site_info)
            return
        if site.FILTER:
            site_info.set_site_status("skipped")
            self.to_gather_que.put(site_info)
            return
        self.to_vectorizer_que.put(site_info)
    
    def _run(self):
        vcf_in = VCF(self.vcf_in_path)
        if self.vcf_idx_path:
            vcf_in.set_index(self.vcf_idx_path)

        for (site_id, site) in enumerate(vcf_in):
            self.push_site(site=site, site_id=site_id)
            if site_id == 5000:
                break
        self.to_vectorizer_que.join()
        self.to_gather_que.join()

def init_manager(batch_size=64, factor=4):
    manager = mp.managers.SyncManager()
    manager.start()
    maxsize = int(round(batch_size * factor))
    manager.to_vectorizer = manager.JoinableQueue(maxsize=maxsize)
    manager.flush_vectorizer = manager.Event()
    manager.to_model = CircularBuffer(ctype=ModelInputStruct, size=maxsize)
    manager.flush_model = manager.Event()
    manager.to_gather = manager.JoinableQueue()
    return manager

def main(ref_path=None, vcf_in_path=None, vcf_idx_path=None, vcf_out_path=None, model_path=None, batch_size=None, klen=None, window=96, n_workers=1):
    manager = init_manager(batch_size=batch_size)
    tokenizer_config = dict(klen=klen)
    vtv_init = lambda: VariantToVectorWorker(
        manager=manager, 
        ref_path=ref_path, 
        tokenizer_config=tokenizer_config, 
        window=window
    )
    vtv_workers = [vtv_init() for x in range(n_workers)]
    model_worker = ModelWorker(manager=manager, model_path=model_path, batch_size=batch_size, klen=klen)
    scatter_worker = ScatterWorker(manager=manager, vcf_in_path=vcf_in_path, vcf_idx_path=vcf_idx_path)
    gather_worker = GatherWorker(manager=manager, vcf_in_path=vcf_in_path, vcf_idx_path=vcf_idx_path, vcf_out_path=vcf_out_path)
    workers = vtv_workers + [model_worker, scatter_worker, gather_worker]
    for worker in workers[::-1]:
        worker.start()
    #
    for worker in workers:
        worker.join()
