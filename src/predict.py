from . models import VariantTokenizer

class VariantToVector(object):
    def __init__(self, ref=None, tokenizer=None, window=96):
        self.ref = ref
        self.tokenizer = tokenizer
        self.window = window
    
    def get_reference_context(self, site_info=None):
        pos = site_info.pos
        chrom = site_info.chrom
        ref = site_info.ref

        start = pos - 1
        end = start + len(ref)
        pre = self.ref[chrom][start - self.window:start]
        post = self.ref[chrom][end:end + self.window]
        return (pre, post)
    
    def process_site(self, site_info=None):
        (pre, post) = self.get_reference_context(site_info=site_info)
        gt_toks = {}
        ref_bases = (site_info.ref, site_info.ref)
        for (genotype_id, gt) in site_info.genotypes.items():
            if (gt.bases == ('.', '.')) or (gt.bases == ref_bases):
                gt_ref = None
            else:
                gt_ref = [f"{pre}{var if var != '.' else site_info.ref}{post}" for var in gt.bases]
                if set(gt_ref[0] + gt_ref[1]) - set('AGTC'):
                    gt_ref = None
                else:
                    gt_ref = self.tokenizer.tokenize(gt_ref)
            gt_toks[genotype_id] = gt_ref
        return gt_toks
