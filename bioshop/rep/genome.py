class Genome(object):
    def __init__(self, name=None, reference=None):
        self.name = name
        self.reference = reference

class GenomeReference(object):
    def __init__(self, path_ref=None):
        self.path_ref = path_ref
        self.fa_refs = Fasta(self.path_ref)
        self._cache = {}

    def lookup(self, site, chrom=None):
        chrom = chrom or site.chrom
        chrom = norm_chrom(chrom, self.fa_refs.keys())
        if chrom not in self._cache:
            self._cache[chrom] = str(self.fa_refs[chrom]).upper()
        seq = self._cache[chrom]
        var_start = site.pos - 1
        var_end = var_start + len(site.ref)
        up = seq[var_start - self.window_size:var_start]
        down = seq[var_end:var_end + self.window_size]
        assert seq[var_start - self.window_size:var_end + self.window_size] == (up + site.ref + down)
        return (up, down)

    __call__ = lookup
