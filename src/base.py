import re

class BaseObject(object):
    pass

class Site(BaseObject):
    def __init__(self,
        site_id=None,
        chrom=None,
        pos=None,
        ref=None,
        filter=None,
        var_type=None,
        genotypes=None,
    ):
        self.site_id = site_id
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.var_type = var_type
        self.genotypes = genotypes
        self.filter = filter

    @classmethod
    def load_from_site(cls, site=None, site_id=None, **kw):
        # load genotypes first
        # XXX: it might be possible to encode phase
        seen = set()
        genotypes = {}
        for (bases, var_type) in zip(site.gt_bases, site.gt_types):
            bases = tuple(re.split('[/|]', bases))
            if bases in seen:
                continue
            seen.add(bases)
            genotype_id = len(genotypes)
            gt = Genotype(
                site_id=site_id, 
                genotype_id=genotype_id,
                bases=bases,
                var_type=var_type,
            )
            genotypes[genotype_id] = gt

        site_obj = cls(
            site_id=site_id,
            chrom=site.CHROM,
            pos=site.POS,
            filter=site.FILTER,
            genotypes=genotypes,
            var_type=site.var_type,
            ref=site.REF,
            **kw
        )
        return site_obj

    @property
    def is_placeholder(self):
        return (self.chrom is None) or \
            (self.pos is None) or \
            (self.ref is None)

    def set_site_status(self, status):
        for gt in self.genotypes.values():
            gt.status = status

    @property
    def is_pending(self):
        if self.is_placeholder:
            return True
        gt_all = [gt.is_pending for gt in self.genotypes.values()]
        return any(gt_all)
    
    def update(self, other):
        assert self.site_id == other.site_id
        if self.placeholder:
            self.chrom = other.chrom
            self.pos = other.pos
            self.ref = other.ref
            self.var_type = other.var_type
        if other.filter is not None:
            self.filter = other.filter
        for genotype in other.genotypes.values():
            self.update_genotype(genotype=genotype)

    def update_genotype(self, genotype=None, genotype_id=None):
        if genotype_id is None:
            genotype_id = genotype.genotype_id
        if genotype_id not in self.genotypes:
            self.genotypes[genotype_id] = genotype
        else:
            self.genotypes[genotype_id].update(genotype)

    def __repr__(self):
        gts = str.join('\n', [f"  {gt}" for gt in self.genotypes.values()])
        msg = f"{self.__class__.__name__}<#{self.site_id}> {self.chrom}:{self.pos}, pending={self.is_pending}\n{gts}"
        return msg

class Genotype(BaseObject):
    def __init__(self,
        site_id=None,
        genotype_id=None,
        var_type=None,
        ref=None,
        bases=None,
        status='pending',
        log_odds=None
    ):
        self.site_id = site_id
        self.genotype_id = genotype_id
        self.var_type = var_type
        self.ref = ref
        self.bases = bases
        self.status = status
        self.log_odds = log_odds

    def __repr__(self):
        msg = f"{self.__class__.__name__}<#{self.site_id},{self.genotype_id}> pending={self.is_pending}, alleles={str.join(', ', self.bases)}, log_odds={self.log_odds}"
        return msg

    @property
    def is_pending(self):
        return self.status == "pending"

    @property
    def is_homref(self):
        # XXX: hardwired
        return self.var_type == 0

    @property
    def is_symbolic(self):
        assert self.ref is not None
        assert self.bases is not None
        symbols = self.ref + str.join('', self.bases)
        symbols = set(symbols.upper())
        concrete = set('AGTC')
        return bool(symbols - concrete)

    @property
    def is_snp(self):
        assert self.ref is not None
        assert self.bases is not None
        return (not self.is_symbolic) and \
            all([len(self.ref) == len(bases) for bases in self.bases])

    @property
    def is_indel(self):
        assert self.ref is not None
        assert self.bases is not None
        return (not self.is_symbolic) and \
            any([len(self.ref) != len(bases) for bases in self.bases])

    def update(self, other):
        assert self.site_id == other.site_id
        assert self.genotype_id == other.genotype_id
        if other.status is not None:
            self.status = other.status
        if other.log_odds is not None:
            self.log_odds = other.log_odds
