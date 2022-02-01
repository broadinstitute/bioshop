import queue
import threading is th
import multiprocessing as mp
from pyfaidx import Fasta
from . predict import VariantToVector, Batcher

class Worker(mp.Process):
    def __init__(self, manager=None, **kw):
        assert manager
        self.manager = manager
        self._running = mp.Event()

    @property
    def running(self):
        return self._running.is_set()
    
    def shutdown(self):
        self._running.clear()

class VariantToVectorWorker(Worker):
    def __init__(self, ref_path=None, vcf_path=None, tokenizer_config=None, window=96, **kw):
        super().__init__(**kw)
        self.ref_path = ref_path
        self.vcf_path = vcf_path
        self.tokenizer_config = tokenizer_config
        self.window = window
        self.in_q = self.manager.to_vectorizer
        self.out_q = self.manager.to_model

    def run(self):
        ref = Fasta(self.ref_path)
        tokenizer = VariantTokenizer(**self.tokenizer_config)
        vectorizer = VariantToVector(ref=ref, tokenizer=tokenizer, window=self.window)

        self._running.set()
        while self.running:
            try:
                ref_key = self.in_q.get(timeout=1)
            except queue.Empty:
                continue
            (chrom, pos) = ref_key
            region = vcf(pos)
            call = next(region)
            vec = vectorizer.process_call(call)
            self.out_q.put(vec)
            

class ModelRunner(Worker):
    def __init__(self, model_path=None, batch_size=None, klen=None, **kw):
        super().__init__(**kw)
        self.model_path = model_path
        self.batch_size = batch_szie
        self.klen = klen
        self.in_q = self.manager.to_model
        self.out_q = self.manager.to_gather

    def batch_thread(self, batcher=None):
        model = models.VariantFilerModel(model_path=self.model_path, klen=self.klen)
        while self.running:
            try:
                call = self.in_q.get(timeout=1)
            except queue.Empty:
                continue
            batcher.add_call(call)

    def run(self):
        batcher = predict.Batcher(model=model, batch_size=batch_size)
        batch_thread = th.Thread(target=self.batch_thread, args=(batcher, ))
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
    def __init__(self, vcf_in_path=None, vcf_out_path=None):
        super().__init__(**kw)
        self.in_q = self.manager.to_gather
        self.call_order = []
        self.pending_calls = {}
        self.ready_calls = {}

    def scatter_calls(self, vcf_in=None):
        for call in vcf_in:
            key = (call.CHR, call.POS)
            self.call_order.append(key)
            if not (call.is_snp or call.is_indel):
                self.ready_calls[key] = call
                continue
            if call.FILTER:
                self.ready_calls[key] = call
                continue
            self.pending_calls[key] = call
            yield key
    
    def post_process_calls(self, call=None):
        var_type = call["var_type"]
        var_idx = 0 if var_type == "snp" else 1
        print(call['locus'], call["var_type"], "REF", call["ref"])
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

    def run(self):
        vcf_in = VCF(self.vcf_in_path)
        itr = self.scatter_calls(vcf_in=vcf_in)
        for key in itr:
            self.out_q.put(key)
            self.gather_calls()

def main(vcf_in_path=None, vcf_out_path=None, model_path=None, batch_size=None, klen=None, window=None, n_workers=4):
    lambda vtv_init: VariantToVectorWorker(
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

