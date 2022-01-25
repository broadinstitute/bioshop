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
        gt_bases = [re.split('[/|]', bases) for bases in call.gt_bases]
        gt_uniqs = set(gq_bases)
        gt_counts = {gt: gt_bases.count(gt) for gt in gt_uniqs}
        return gt_counts

    def process_call(self, call=None):
        (pre, post) = self.get_reference_context(call)
        gt_counts = self.get_genotypes(call)
        gt_inps = {}
        for gt in gt_counts:
            gt_ref = [f"{pre}{var}{post}" for var in gt]
            if set(gt_ref[0] + gt_ref[1]) - set('AGTC'):
                gt_ref = None
            else:
                gt_ref = self.tokenizer.tokenize(gt_ref)
            gt_inps[gt] = gt_ref
        ret = {
            "locus": f"{call.CHROM}:{call.POS}",
            "gt_inps": gt_inps,
            "gt_counts": gt_counts,
        }
        return ret

def build_worker(ref_path=None, klen=None):
    from pyfaidx import Fasta
    ref = Fasta(ref_path)
    tokenizer = VariantTokenizer(klen=klen)
    worker = VariantWorker(ref=ref, tokenizer=tokenizer)
    return worker
