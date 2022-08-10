from . alignment import RefAltAlignment, AltAltAlignment, LiteralFingerprint, CigarFingerprint

class AlleleFingerprint(object):
    def __init__(self, chrom=None, pos=None, ref=None, alt=None, flanks=None):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.flanks = flanks
        self.alignment = RefAltAlignment(allele=self)
    
    def __str__(self):
        return str(self.literal_fingerprint)
    
    def __repr__(self):
        return str(self)
    
    def match(self, other, debug=False):
        al = AltAltAlignment(self, other)
        if debug:
            al.debug_match()
        return al.is_match

    @property
    def alt_span(self):
        up_len = len(self.flanks[0])
        alt_len = len(self.alt)
        coords = (up_len, up_len + alt_len)
        assert self.altseq[slice(*coords)] == self.alt
        return coords

    @property
    def refseq(self):
        return  self.flanks[0] + self.ref + self.flanks[1] 
    
    @property
    def altseq(self):
        return self.flanks[0] + self.alt + self.flanks[1]

    @property
    def cigar_fingerprint(self):
        offset = self.pos - len(self.flanks[0])
        inner = self.alignment.cigar
        if inner.parts[0][1] == '=':
            offset += inner.parts[0][0]
            inner = inner[1:]
        if inner.parts[-1][1] == '=':
            inner = inner[:-1]
        return CigarFingerprint(chrom=self.chrom, pos=offset, cigar=inner)

    @property
    def literal_fingerprint(self):
        return LiteralFingerprint(chrom=self.chrom, pos=self.pos, ref=self.ref, alt=self.alt)

    @classmethod
    def from_site(cls, site=None, alt=None, flanks=None, chrom=None):
        chrom = chrom or site.chrom
        return cls(
            chrom=chrom,
            pos=site.pos,
            ref=site.ref,
            alt=alt,
            flanks=flanks,
        )
