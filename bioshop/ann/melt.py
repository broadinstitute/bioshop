from .. rep.region import Region

import random
from Bio.SeqUtils.MeltingTemp import Tm_NN
import numpy as np

def degen_resolver(seq=None, method='to_any'):
    # XXX: set a random seed?
    # XXX: case sensitive?
    # XXX: speed?
    canonical = tuple('CGAT')
    if method.lower() == 'to_any':
        acceptable = set(canonical + ('N',))
        resolver = lambda bs: bs if bs in acceptable else 'N'
    elif method.lower() == 'to_random':
        acceptable = set(canonical)
        resolver = lambda bs: bs if bs in acceptable else random.choice(canonical)
    else:
        raise TypeError(method)
    seq = map(resolver, seq.upper())
    return str.join('', seq)

# made from four different cheeses
class OligoMelt(object):
    def __init__(self, reference=None, read_len=256, window_len=32, step_len=16, degen_method='to_random'):
        self.reference = reference
        self.read_len = read_len
        self.window_len = window_len
        self.step_len = step_len
        self.degen_method = degen_method

    def melt(self, chrom=None, pos=None):
        seq = self.reference[chrom]
        start = max(0, pos - self.read_len)
        stop = min(len(seq) - 1, pos + self.read_len)
        region = Region(f'{chrom}:{start}-{stop}')
        samples = region.window(width=self.window_len, step=self.step_len)
        melts = []
        for sample in samples:
            subseq = seq[sample.start:sample.stop]
            subseq = degen_resolver(subseq, method=self.degen_method)
            tm = Tm_NN(subseq)
            melts.append(tm)
        melts = np.array(melts)
        return (np.mean(melts), np.std(melts))

    # at some point, this can move to a chr / pos lookup, 
    # but by using the site object from the VCF, we can 
    # assume (site.pos - 1) is the correct thing to do
    def melt_site(self, site=None):
        return self.melt(chrom=site.chrom, pos=site.pos - 1)

