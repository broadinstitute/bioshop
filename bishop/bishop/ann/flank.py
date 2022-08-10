class VariantFlanks(object):
    def __init__(self, assembly=None, flank_len=50, as_scheme=None):
        self.assembly = assembly
        self.as_scheme = as_scheme
        self.flank_len = flank_len
        self._cache = {}

    def get_flanks(self, site=None):
        chrom = site.chrom
        if self.as_scheme:
            chrom = self.assembly.as_scheme(chrom, as_scheme=self.as_scheme)
        if chrom not in self._cache:
            self._cache[chrom] = self.assembly.genome[chrom].upper()
        seq = self._cache[chrom]
        var_start = site.pos - 1
        var_end = var_start + len(site.ref)

        up = seq[var_start - self.flank_len:var_start]
        down = seq[var_end:var_end + self.flank_len]

        assert seq[var_start - self.flank_len:var_end + self.flank_len] == (up + site.ref + down)
        return {'chrom': chrom, 'flank_up': up, 'flank_down': down}
