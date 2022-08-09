class AlleleIndex(object):
    def __init__(self, alleles=None):
        self.alleles = tuple(alleles)
        self.index_alleles()

    def index_alleles(self):
        self.fingerprints = {}
        for allele in self.alleles:
            # note: it is impossible for the two types of fingerprints to overlap
            self.fingerprints[allele.literal_fingerprint] = allele
            self.fingerprints[allele.cigar_fingerprint] = allele
        
    def match(self, allele=None, debug=False):
        if allele.literal_fingerprint in self.fingerprints:
            # this is only an assert, but it is also critical to the conditional branch
            assert allele.cigar_fingerprint == self.fingerprints[allele.literal_fingerprint].cigar_fingerprint
            return True
        if allele.cigar_fingerprint in self.fingerprints:
            matched_allele = self.fingerprints[allele.cigar_fingerprint]
            # If the literal fingerprint shape matches, but the subsequent key does not
            # then, by definition, the Cigar fingerprint will not match either.  This mainly
            # is a relief valve for SNPs that do not match the literal fingerprint and trickle down
            # to Cigar check 
            if allele.literal_fingerprint.chrom == matched_allele.literal_fingerprint.chrom and \
                allele.literal_fingerprint.pos == matched_allele.literal_fingerprint.pos and \
                len(allele.literal_fingerprint.ref) == len(matched_allele.literal_fingerprint.ref) and \
                len(allele.literal_fingerprint.alt) == len(matched_allele.literal_fingerprint.alt):
                    return False
            if allele.match(matched_allele, debug=debug):
                return True
        # no match
        return False
