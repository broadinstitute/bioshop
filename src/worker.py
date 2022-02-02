import queue
import pprint
from threading import Thread
import multiprocessing.managers
import multiprocessing as mp
from . models import VariantTokenizer, VariantFilterModel
from pyfaidx import Fasta
from . predict import VariantToVector, Batcher
from cyvcf2 import VCF

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
        self._running.set()
        msg = f"Starting {self.__class__.__name__} worker run={self.running}"
        print(msg)
        self._run()

class VariantToVectorWorker(Worker):
    def __init__(self, ref_path=None, vcf_path=None, tokenizer_config=None, window=96, **kw):
        super().__init__(**kw)
        self.ref_path = ref_path
        self.vcf_path = vcf_path
        self.tokenizer_config = tokenizer_config
        self.window = window
        self.in_q = self.manager.to_scatter
        self.out_q = self.manager.to_model

    def _run(self):
        ref = Fasta(self.ref_path)
        vcf = VCF(self.vcf_path)
        tokenizer = VariantTokenizer(**self.tokenizer_config)
        vectorizer = VariantToVector(ref=ref, tokenizer=tokenizer, window=self.window)

        self._running.set()
        while self.running:
            try:
                ref_key = self.in_q.get(timeout=1)
            except queue.Empty:
                print(f"{self.__class__.__name__} loop: in_q empty. run={self.running}")
                continue
            (chrom, pos) = ref_key
            site = str.join(':', (chrom, str(pos)))
            region = vcf(site)
            call = next(region)
            while call.POS < pos:
                call = next(region)
            assert call.POS == pos, f"{call.POS} != {pos}"
            vec = vectorizer.process_call(call)
            self.out_q.put(vec)
            
class ModelRunner(Worker):
    def __init__(self, model_path=None, batch_size=None, klen=None, **kw):
        super().__init__(**kw)
        self.model_path = model_path
        self.batch_size = batch_size
        self.klen = klen
        self.in_q = self.manager.to_model
        self.out_q = self.manager.to_gather

    def batch_thread(self, batcher=None):
        while self.running:
            try:
                call = self.in_q.get(timeout=1)
            except queue.Empty:
                print(f"{self.__class__.__name__} thread: empty queue, run={self.running}")
                continue
            batcher.add_call(call)

    def _run(self):
        model = VariantFilterModel(model_path=self.model_path, klen=self.klen)
        batcher = Batcher(model=model, batch_size=self.batch_size)
        batch_thread = Thread(target=self.batch_thread, args=(batcher, ))
        batch_thread.start()
        while self.running:
            batcher.do_batch()
            calls = batcher.get_completed_calls()
            for call in calls:
                self.out_q.put(call)
        batch_thread.join()
        batcher.flush()
        calls = batcher.get_completed_calls()
        for call in calls:
            self.out_q.put(call)

class MainRunner(Worker):
    def __init__(self, vcf_in_path=None, vcf_out_path=None, **kw):
        super().__init__(**kw)
        self.vcf_in_path = vcf_in_path
        self.vcf_out_path = vcf_out_path
        self.out_q = self.manager.to_scatter
        self.in_q = self.manager.to_gather
        self.call_order = []
        self.pending_calls = {}
        self.ready_calls = {}

    def scatter_calls(self, vcf_in=None):
        for call in vcf_in:
            key = (call.CHROM, call.POS)
            self.call_order.append(key)
            if not (call.is_snp or call.is_indel):
                self.ready_calls[key] = call
                continue
            if call.FILTER:
                self.ready_calls[key] = call
                continue
            self.pending_calls[key] = call
            yield key
    
    def post_process_call(self, call=None):
        var_type = call["var_type"]
        var_idx = 0 if var_type == "snp" else 1
        print(call['locus'], call["var_type"], "REF", call["ref"])
        return
        for al in call['gt_scores']:
            if call['gt_scores'][al] == None:
                continue
            call['gt_scores'][al]['score'] = call['gt_scores'][al]['logodds'][var_idx]
        pprint.pprint(call['gt_scores'])
        print()

    def gather_calls(self):
        while True:
            try:
                call_info = self.in_q.get(block=False)
            except queue.Empty:
                return
            self.post_process_call(call_info)

    def _run(self):
        vcf_in = VCF(self.vcf_in_path)
        itr = self.scatter_calls(vcf_in=vcf_in)
        for key in itr:
            self.out_q.put(key)
            self.gather_calls()

def init_manager(batch_size=16, factor=8):
    manager = mp.managers.SyncManager()
    manager.start()
    maxsize = int(round(batch_size * factor))
    manager.to_scatter = manager.Queue(maxsize=maxsize)
    manager.to_model = manager.Queue(maxsize=maxsize)
    manager.to_gather = manager.Queue(maxsize=maxsize)
    return manager

def main(ref_path=None, vcf_in_path=None, vcf_out_path=None, model_path=None, batch_size=None, klen=None, window=96, n_workers=4):
    manager = init_manager(batch_size=batch_size)
    tokenizer_config = dict(klen=klen)
    vtv_init = lambda: VariantToVectorWorker(
        manager=manager, 
        ref_path=ref_path, 
        vcf_path=vcf_in_path, 
        tokenizer_config=tokenizer_config, 
        window=window
    )
    vtv_workers = [vtv_init() for x in range(n_workers)]
    model_worker = ModelRunner(manager=manager, model_path=model_path, batch_size=batch_size, klen=klen)
    main_worker = MainRunner(manager=manager, vcf_in_path=vcf_in_path, vcf_out_path=vcf_out_path)
    workers = vtv_workers + [model_worker, main_worker]
    for worker in workers:
        worker.start()
    #
    for worker in workers:
        worker.join()
