import re
import random
from cyvcf2 import VCF
from sklearn.preprocessing import QuantileTransformer
import pandas as pd

class VariantSampler(object):
    def random_region(self, vcf, seqmap):
        (seqname, seqlen) = random.choice(seqmap)
        pos = random.randint(0, seqlen)
        coord = f"{seqname}:{pos}"
        region = vcf(coord)
        return region

    def sample(self, vcf=None, n_sites=100, n_samples=100):
        re_chr = re.compile("^chr[\dX]{1,2}$")
        seqnames = filter(re_chr.match, vcf.seqnames)
        seqmap = list(zip(seqnames, vcf.seqlens))

        for site_idx in range(n_sites):
            sample_cnt = 0
            region = self.random_region(vcf, seqmap)
            while sample_cnt < n_samples:
                try:
                    rec = next(region)
                except StopIteration:
                    region = self.random_region(vcf, seqmap)
                    continue
                yield rec
                sample_cnt += 1

class AnnotationTransformer(object):
    DefaultAnnotations = ("DP", "QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR")

    def __init__(self, method="quantile", ann_list=None):
        if method == "quantile":
            self.transformer = QuantileTransformer()
        else:
            raise ValueError(method)
        self.method = method
        self.ann_list = ann_list or self.DefaultAnnotations
        self._trained = False

    @property
    def is_trained(self):
        return self._trained

    def transform(self, *args, **kw):
        return self.transformer.transform(*args, **kw)
    
    def __getstate__(self):
        return {
            'method': self.method,
            'ann_list': self.ann_list,
            'transformer': self.transformer,
            '_trained': self._trained,
        }
    
    def __setstate__(self, state):
        self.method = state['method']
        self.ann_list = state['ann_list']
        self.transformer = state['transformer']
        self._trained = state['_trained']

    def get_annotations(self, vcf_site=None):
        info = vcf_site.INFO
        info = {key: info.get(key, None) for key in self.ann_list}
        return info

    def sample_annotations(self, vcf=None, n_sites=100, n_samples=100):
        def get_info(itr=None):
            for vcf_site in itr:
                info = self.get_annotations(vcf_site=vcf_site)
                yield info

        vs = VariantSampler()
        itr = vs.sample(vcf=vcf, n_sites=n_sites, n_samples=n_samples)
        itr = get_info(itr=itr)
        samples = list(itr)
        return pd.DataFrame(samples).fillna(0)

    def train(self, vcf_path=None, n_sites=100, n_samples=100):
        assert not self.is_trained
        msg = f"Training numerical transformer with {n_samples * n_sites} variants from {n_sites} random sites"
        print(msg)
        vcf = VCF(vcf_path)
        df = self.sample_annotations(vcf=vcf)
        samples = df[list(self.ann_list)]
        self.transformer.fit(samples)
        self._trained = True
