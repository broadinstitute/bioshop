import re
import torch
import numpy as np
from . models import VariantTokenizer

class VariantWorker(object):
    def __init__(self, ref=None, tokenizer=None, window=96):
        self.ref = ref
        self.tokenizer = tokenizer
        self.window = window
    
    def get_reference_context(self, call=None):
        start = call.POS - 1
        end = start + len(call.REF)
        pre = self.ref[call.CHROM][start - self.window:start]
        post = self.ref[call.CHROM][end:end + self.window]
        return (pre, post)
    
    def get_genotypes(self, call=None):
        # XXX: it might be possible to encode phase
        gt_bases = [tuple(re.split('[/|]', bases)) for bases in call.gt_bases]
        gt_uniqs = set(gt_bases)
        gt_counts = {gt: gt_bases.count(gt) for gt in gt_uniqs}
        return gt_counts

    def process_call(self, call=None):
        (pre, post) = self.get_reference_context(call)
        gt_counts = self.get_genotypes(call)
        gt_inps = {}
        for gt in gt_counts:
            gt_ref = [f"{pre}{var if var != '.' else call.REF}{post}" for var in gt]
            if set(gt_ref[0] + gt_ref[1]) - set('AGTC'):
                gt_ref = None
            else:
                gt_ref = self.tokenizer.tokenize(gt_ref)
            gt_inps[gt] = gt_ref
        ret = {
            "locus": f"{call.CHROM}:{call.POS}",
            "gt_inps": gt_inps,
            "gt_counts": gt_counts,
            "var_type": call.var_type
        }
        return ret

class Batcher(object):
    def __init__(self, batch_size=120):
        self.batch_size = batch_size
        self.pending_calls = {}
        self.batch_pool = []
        self.completed_calls = []

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
    
    def flush(self):
        self.do_batch(force=True)

    def add_call(self, call):
        call['gt_scores'] = {}
        key = call['locus']
        self.pending_calls[key] = call
        for (gt, toks) in call['gt_inps'].items():
            if toks is None:
                call['gt_scores'][gt] = None
                continue
            key = (call['locus'], gt)
            pair = (key, toks)
            self.batch_pool.append(pair)

    def get_completed_calls(self):
        ret = self.completed_calls
        self.completed_calls = []
        return ret

    def process_batch(self, batch):
        return [[1,1,1,1] for it in batch]
    
    def do_batch(self, force=False):
        if force or not self.batch_ready:
            return
        batch = self.batch_pool[:self.batch_size]
        self.batch_pool = self.batch_pool[self.batch_size:]
        (keys, batch) = list(zip(*batch))
        batch = torch.tensor(batch)
        batch_results = self.process_batch(batch)
        for (idx, key) in enumerate(keys):
            (locus, gt) = key
            call = self.pending_calls[locus]
            outp = {k: np.squeeze(v[idx, :]) for (k, v) in batch_results.items()}
            call['gt_scores'][gt] = outp
            if len(call['gt_scores']) == len(call['gt_inps']):
                del self.pending_calls[locus]
                self.completed_calls.append(call)

def build_worker(ref_path=None, klen=None):
    from pyfaidx import Fasta
    ref = Fasta(ref_path)
    tokenizer = VariantTokenizer(klen=klen)
    worker = VariantWorker(ref=ref, tokenizer=tokenizer)
    return worker
