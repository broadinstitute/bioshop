import os
import re
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
from . predict import VariantToVector, Batcher
from . cbuf import CircularBuffer

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
            self._run()
        os.chdir(pwd)

class VariantToVectorWorker(Worker):
    def __init__(self, ref_path=None, tokenizer_config=None, window=96, **kw):
        super().__init__(**kw)
        self.ref_path = ref_path
        self.tokenizer_config = tokenizer_config
        self.window = window
        self.in_q = self.manager.to_vectorizer
        self.to_model = self.manager.to_model
        self.to_gather = self.manager.to_gather

    def build_model_inputs(self, call_info=None, gt_inps=None):
        call_id = call_info["call_id"]
        results = []
        model_inputs = []
        for (genotype_id, gt_inp) in enumerate(gt_inps):
            if gt_inp is None:
                results.append(dict(skipped=True))
                continue
            results.append(dict(pending=True))
            gt_inp["call_id"] = call_info["call_id"]
            gt_inp["genotype_id"] = genotype_id
            inp = ModelInputStruct(**gt_inp)
            model_inputs.append(inp)
        call_info["gt_results"] = results
        return (call_info, model_inputs)

    def _run(self):
        ref = Fasta(self.ref_path)
        tokenizer = VariantTokenizer(**self.tokenizer_config)
        vectorizer = VariantToVector(ref=ref, tokenizer=tokenizer, window=self.window)

        self._running.set()
        while self.running:
            try:
                call_info = self.in_q.get(timeout=1)
            except queue.Empty:
                print(f"{self.__class__.__name__} loop: in_q empty. run={self.running}")
                continue

            gt_inps = vectorizer.process_call(call_info)
            (call_info, model_inputs) = self.build_model_inputs(call_info=call_info, gt_inps=gt_inps)
            self.to_gather.put(call_info)
            for inp in model_inputs:
                self.to_model.push(inp)
            
class ModelWorker(Worker):
    def __init__(self, model_path=None, batch_size=None, klen=None, **kw):
        super().__init__(**kw)
        self.model_path = model_path
        self.batch_size = batch_size
        self.klen = klen
        self.in_q = self.manager.to_model
        self.out_q = self.manager.to_gather

    def get_batch(self):
        # XXX: need to flush somehow
        batch = [self.in_q.pop() for _ in range(self.batch_size)]
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
            outp = model.predict(batch)
            result = []
            for (log_odds, keys) in zip(outp['log_odds'], batch_keys):
                ns = keys.copy()
                ns['log_odds'] = log_odds.tolist()
                result.append(ns)
            self.out_q.put(result)

class GatherWorker(Worker):
    def __init__(self, vcf_in_path=None, vcf_idx_path=None, vcf_out_path=None, **kw):
        super().__init__(**kw)
        self.vcf_in_path = vcf_in_path
        self.vcf_idx_path = vcf_idx_path
        self.vcf_out_path = vcf_out_path
        self.in_q = self.manager.to_gather
        self.ready_calls = {}

    def post_process_call(self, call_info=None):
        if 'skipped' not in call_info:
            var_type = call_info["var_type"]
            var_idx = 0 if var_type == "snp" else 1
            for al in call_info['gt_scores']:
                if call_info['gt_scores'][al] == None:
                    continue
                call_info['gt_scores'][al]['score'] = call_info['gt_scores'][al]['log_odds'][var_idx]
        call_id = call_info['call_id']
        self.ready_calls[call_id] = call_info

    def wait_on(self, call_id=None):
        while True:
            try:
                call_info = self.in_q.get(timeout=1)
                print(call_info)
            except queue.Empty:
                continue
            #self.post_process_call(call_info=call_info)
            #self.in_q.task_done()
            if call_id in self.ready_calls:
                return

    def _run(self):
        vcf_in = VCF(self.vcf_in_path)
        if self.vcf_idx_path:
            vcf_in.set_index(self.vcf_idx_path)

        for (call_id, call) in enumerate(vcf_in):
            if call_id not in self.ready_calls:
                self.wait_on(call_id=call_id)
            call_info = self.ready_calls.pop(call_id)
            assert call_id == call_info['call_id']
            assert call.POS == call_info['pos']
            assert call.CHROM == call_info['chrom']
            #print(call_id, call.CHROM, call.POS)

class ScatterWorker(Worker):
    def __init__(self, vcf_in_path=None, vcf_idx_path=None, **kw):
        super().__init__(**kw)
        self.vcf_in_path = vcf_in_path
        self.vcf_idx_path = vcf_idx_path
        self.to_vectorizer_que = self.manager.to_vectorizer
        self.to_gather_que = self.manager.to_gather

    def get_genotypes(self, call=None):
        # XXX: it might be possible to encode phase
        gt_bases = [tuple(re.split('[/|]', bases)) for bases in call.gt_bases]
        gt_list = list(set(gt_bases))
        return gt_list

    def process_call(self, call=None, call_id=None):
        gt_list = self.get_genotypes(call)
        gt_inps = {}
        ref_key = (call.REF, call.REF)
        for gt in gt_list:
            if (gt == ('.', '.')) or (gt == ref_key):
                gt_inps[gt] = None
                continue
        call_info = {
            "call_id": call_id,
            "chrom": call.CHROM,
            "pos": call.POS,
            "locus": f"{call.CHROM}:{call.POS}",
            "gt_list": gt_list,
            "var_type": call.var_type,
            "ref": call.REF,
        }
        return call_info

    def push_call(self, call=None, call_id=None):
        call_info = {"chrom": call.CHROM, "pos": call.POS, "call_id": call_id}
        if not (call.is_snp or call.is_indel):
            call_info["skipped"] = "var_type"
            self.to_gather_que.put(call_info)
            return
        if call.FILTER:
            call_info["skipped"] = "filtered"
            self.to_gather_que.put(call_info)
            return
        call_info = self.process_call(call=call, call_id=call_id)
        self.to_vectorizer_que.put(call_info)
    
    def _run(self):
        vcf_in = VCF(self.vcf_in_path)
        if self.vcf_idx_path:
            vcf_in.set_index(self.vcf_idx_path)

        for (call_id, call) in enumerate(vcf_in):
            self.push_call(call=call, call_id=call_id)

def init_manager(batch_size=64, factor=64):
    manager = mp.managers.SyncManager()
    manager.start()
    maxsize = int(round(batch_size * factor))
    manager.to_vectorizer = manager.Queue(maxsize=maxsize)
    manager.to_model = CircularBuffer(ctype=ModelInputStruct, size=maxsize)
    manager.to_gather = manager.Queue()
    return manager

def main(ref_path=None, vcf_in_path=None, vcf_idx_path=None, vcf_out_path=None, model_path=None, batch_size=None, klen=None, window=96, n_workers=2):
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
    for worker in workers:
        worker.start()
    #
    for worker in workers:
        worker.join()
