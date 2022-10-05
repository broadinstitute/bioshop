import random

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

class VariantFlanks(object):
    def __init__(self, reference=None, flank_len=50, degen_method='to_any'):
        self.reference = reference
        self.flank_len = flank_len
        self.degen_method = degen_method
        self._cache = {}

    # at some point, this can move to a chr / pos lookup, 
    # but by using the site object from the VCF, we can do to things:
    #
    #  1) assume (site.pos - 1) is the correct thing to do
    #  2) assert (flank_up + ref + flank_down) == genome[x:y]
    #
    def get_flanks(self, site=None):
        chrom = site.chrom
        if chrom not in self._cache:
            seq = self.reference[chrom].upper()
            if self.degen_method:
                seq = degen_resolver(seq, method=self.degen_method)
            self._cache[chrom] = seq
        seq = self._cache[chrom]
        var_start = site.pos - 1
        var_end = var_start + len(site.ref)

        up = seq[var_start - self.flank_len:var_start]
        down = seq[var_end:var_end + self.flank_len]

        #assert seq[var_start - self.flank_len:var_end + self.flank_len] == (up + site.ref + down)
        if seq[var_start - self.flank_len:var_end + self.flank_len] != (up + site.ref + down):
            msg = f'Reference and VCF mismatch ({chrom}:{var_start}-{var_end}) != {site.ref})'
            raise ValueError(msg)
        #return {'chrom': chrom, 'flank_up': up, 'flank_down': down}
        return (up, down)
