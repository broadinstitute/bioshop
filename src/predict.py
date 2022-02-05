import re
import torch
import numpy as np
from threading import Lock, Condition
from . models import VariantTokenizer

class VariantToVector(object):
    def __init__(self, ref=None, tokenizer=None, window=96):
        self.ref = ref
        self.tokenizer = tokenizer
        self.window = window
    
    def get_reference_context(self, call_info=None):
        pos = call_info['pos']
        chrom = call_info['chrom']
        ref = call_info['ref']

        start = pos - 1
        end = start + len(ref)
        pre = self.ref[chrom][start - self.window:start]
        post = self.ref[chrom][end:end + self.window]
        return (pre, post)
    
    def process_call(self, call_info=None):
        (pre, post) = self.get_reference_context(call_info=call_info)
        gt_list = call_info['gt_list']
        gt_inps = []
        ref_key = (call_info['ref'], call_info['ref'])
        for gt in gt_list:
            if (gt == ('.', '.')) or (gt == ref_key):
                gt_inps.append(None)
                continue
            gt_ref = [f"{pre}{var if var != '.' else call.REF}{post}" for var in gt]
            if set(gt_ref[0] + gt_ref[1]) - set('AGTC'):
                gt_ref = None
            else:
                gt_ref = self.tokenizer.tokenize(gt_ref)
            gt_inps.append(gt_ref)
        return gt_inps

class Batcher(object):
    def __init__(self, model=None, batch_size=None):
        self.model = model
        self.batch_size = batch_size
        self.pending_calls = dict()
        self.batch_pool = list()
        self.completed_calls = list()
        self._lock = Lock()
        self._batch_ready_cv = Condition(lock=self._lock)

    @property
    def num_pending_batches(self):
        return len(self.batch_pool) // self.batch_size

    @property
    def num_pending_calls(self):
        return len(self.pending_calls)

    @property
    def num_completed_calls(self):
        return len(self.completed_calls)

    @property
    def batch_ready(self):
        return bool(self.num_pending_batches)

    @property
    def completed_ready(self):
        return bool(self.num_completed_calls)
    
    def add_call(self, call):
        call['gt_scores'] = {}
        call_id = call['call_id']
        pool = []
        for (gt, toks) in call['gt_inps'].items():
            if toks is None:
                call['gt_scores'][gt] = None
                continue
            key = (call_id, gt)
            pair = (key, toks)
            pool.append(pair)
        if not pool:
            self.completed_calls.append(call)
            return
        with self._batch_ready_cv:
            self.batch_pool += pool
            assert call_id not in self.pending_calls, f"{locus} in self.pending_calls"
            self.pending_calls[call_id] = call
            if self.batch_ready:
                self._batch_ready_cv.notify()

    def get_completed_calls(self):
        with self._lock:
            ret = self.completed_calls
            self.completed_calls = []
        return ret

    def pivot_batch(self, batch):
        t_batch = {key: list() for key in batch[0]}
        for item in batch:
            for key in item:
                t_batch[key].append(item[key])
        for key in t_batch:
            t_batch[key] = torch.tensor(t_batch[key])
        return t_batch

    def get_batch(self, force=False):
        with self._batch_ready_cv:
            if force or not self.batch_ready:
                # XXX: timeout?
                print("starving")
                ok = self._batch_ready_cv.wait(timeout=1)
                if not ok:
                    return None
            batch = self.batch_pool[:self.batch_size]
            self.batch_pool = self.batch_pool[self.batch_size:]
        (keys, batch) = list(zip(*batch))
        batch = self.pivot_batch(batch)
        return (keys, batch)

    def postprocess_batch(self, keys=None, results=None):
        with self._lock:
            for (idx, key) in enumerate(keys):
                (call_id, gt) = key
                call = self.pending_calls[call_id]
                outp = {k: np.squeeze(v[idx, :]) for (k, v) in results.items()}
                call['gt_scores'][gt] = outp
                if len(call['gt_scores']) == len(call['gt_inps']):
                    del self.pending_calls[call_id]
                    self.completed_calls.append(call)

    def do_batch(self, force=False):
        batch = self.get_batch(force=force)
        if batch is None:
            return
        (keys, batch) = batch
        results = self.model.predict(batch)
        self.postprocess_batch(keys=keys, results=results)

    def flush(self):
        self.do_batch(force=True)

def build_worker(ref_path=None, klen=None):
    from pyfaidx import Fasta
    ref = Fasta(ref_path)
    tokenizer = VariantTokenizer(klen=klen)
    worker = VariantWorker(ref=ref, tokenizer=tokenizer)
    return worker
