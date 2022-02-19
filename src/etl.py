import re
import glob

class TrainingVcfLoader(object):
    Header = ('chrom', 'pos', 'ref', 'alt', 'var_type', 'label')

    def iter_vcf(self, vcf=None, region=None):
        pass

class TrainingVcfEvalLoader(TrainingVcfLoader):

    def __init__(self, vcfeval_call_annot='CALL'):
        self.vcfeval_call_annot = vcfeval_call_annot
    
    def get_label(self, site=None):
        if site.is_snp:
            var_type = "SNP"
        elif site.is_indel:
            var_type = "INDEL"
        else:
            return None
        truth = site.INFO.get(self.vcfeval_call_annot)
        if truth not in ('FP', 'TP'):
            return None
        label_name = f"{var_type}_{truth}"
        label_info = {
            "var_type": var_type,
            "label": label_name,
        }
        return label_info

    def process_site(self, site=None):
        label = self.get_label(site=site)
        if not label:
            return None
        ret = {
            "chrom": site.CHROM,
            "pos": site.POS,
            "ref": site.REF,
            "alt": str.join(',', site.ALT),
        }
        ret.update(label)
        return ret

    def iter_vcf(self, vcf=None, region=None):
        #sample_idx = vcf.samples.index(self.vcfeval_sample_name)
        itr = vcf(region) if region is not None else iter(vcf)
        for site in itr:
            ret = self.process_site(site=site)
            if ret is None:
                continue
            yield ret
