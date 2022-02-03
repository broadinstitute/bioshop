import os
import queue
import pprint
from tempfile import TemporaryDirectory
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
        pwd = os.getcwd()
        with TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            self._running.set()
            msg = f"Starting {self.__class__.__name__} worker run={self.running}, tmpdir={tmpdir}"
            print(msg)
            self._run()
        os.chdir(pwd)

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
                work = self.in_q.get(timeout=1)
            except queue.Empty:
                #print(f"{self.__class__.__name__} loop: in_q empty. run={self.running}")
                continue
            (chrom, pos) = (work['chrom'], work['pos'])
            site = f"{chrom}:{pos}"
            region = vcf(site)
            call = next(region)
            while call.POS < pos:
                call = next(region)
            assert call.POS == pos, f"{call.POS} != {pos}"
            vec = vectorizer.process_call(call)
            vec['call_id'] = work['call_id']
            self.out_q.put(vec)
            self.in_q.task_done()
            
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
                #print(f"{self.__class__.__name__} thread: empty queue, run={self.running}")
                continue
            batcher.add_call(call)
            self.in_q.task_done()

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
    def __init__(self, vcf_in_path=None, vcf_idx_path=None, vcf_out_path=None, **kw):
        super().__init__(**kw)
        self.vcf_in_path = vcf_in_path
        self.vcf_idx_path = vcf_idx_path
        self.vcf_out_path = vcf_out_path
        self.out_q = self.manager.to_scatter
        self.in_q = self.manager.to_gather
        self.next_call_id = 0
        self.pending_calls = {}
        self.ready_calls = {}

    def post_process_call(self, call=None):
        var_type = call["var_type"]
        var_idx = 0 if var_type == "snp" else 1
        for al in call['gt_scores']:
            if call['gt_scores'][al] == None:
                continue
            call['gt_scores'][al]['score'] = call['gt_scores'][al]['logodds'][var_idx]
        call_id = call['call_id']
        pending_call = self.pending_calls.pop(call_id)
        self.ready_calls[call_id] = (pending_call, call)
        #print(call['locus'], call["var_type"], "REF", call["ref"])
        #pprint.pprint(call['gt_scores'])
        #print()

    def push_call(self, call=None, call_id=None):
        call_info = dict(call_id=call_id)
        if not (call.is_snp or call.is_indel):
            call_info["skipped"] = "var_type"
            self.ready_calls[call_id] = (call, call_info)
            return
        if call.FILTER:
            call_info["skipped"] = "filtered"
            self.ready_calls[call_id] = (call, call_info)
            return
        self.pending_calls[call_id] = call
        work = {"chrom": call.CHROM, "pos": call.POS, "call_id": call_id}
        self.out_q.put(work)
    
    def gather_calls(self):
        ret = []
        print("gather")
        while True:
            try:
                call_info = self.in_q.get(timeout=1)
            except queue.Empty:
                break
            self.post_process_call(call_info)
            self.in_q.task_done()
        while self.next_call_id in self.ready_calls:
            call_info = self.ready_calls.pop(self.next_call_id)
            ret.append(call_info)
            self.next_call_id += 1
        return ret

    def scatter_thread(self):
        print("Starting scatter thread")
        vcf_in = VCF(self.vcf_in_path)
        if self.vcf_idx_path:
            vcf_in.set_index(self.vcf_idx_path)

        for (call_id, call) in enumerate(vcf_in):
            self.push_call(call=call, call_id=call_id)

    def _run(self):
        scatter_thread = Thread(target=self.scatter_thread)
        scatter_thread.start()
        self.next_call_id = 0
        while True:
            calls_out = self.gather_calls()
            for (vcf_call, predict_call) in calls_out:
                print(predict_call['call_id'], vcf_call.CHROM, vcf_call.POS)

def init_manager(batch_size=32, factor=8):
    manager = mp.managers.SyncManager()
    manager.start()
    maxsize = None
    maxsize = int(round(batch_size * factor))
    manager.to_scatter = manager.Queue(maxsize=maxsize)
    manager.to_model = manager.Queue(maxsize=maxsize)
    manager.to_gather = manager.Queue(maxsize=maxsize)
    return manager

def main(ref_path=None, vcf_in_path=None, vcf_idx_path=None, vcf_out_path=None, model_path=None, batch_size=None, klen=None, window=96, n_workers=4):
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
    main_worker = MainRunner(manager=manager, vcf_in_path=vcf_in_path, vcf_idx_path=vcf_idx_path, vcf_out_path=vcf_out_path)
    #workers = vtv_workers + [model_worker, main_worker]
    workers = vtv_workers + [model_worker]
    for worker in workers:
        worker.start()
    main_worker.run()
    #
    for worker in workers:
        worker.join()
