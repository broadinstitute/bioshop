import re
import numpy as np

class BaseObject(object):
    pass

class Site(BaseObject):
    def __init__(self,
        site_id=None,
        chrom=None,
        pos=None,
        ref=None,
        info=None,
        filter=None,
        var_type=None,
        genotypes=None,
    ):
        self.site_id = site_id
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.var_type = var_type
        self.info = info
        self.genotypes = genotypes
        self.filter = filter

    @classmethod
    def load_from_site(cls, site=None, site_id=None, **kw):
        # load genotypes first
        # XXX: it might be possible to encode phase
        seen = set()
        genotypes = {}
        # XXX: this is code with possible intention for
        # genotype specific information
        #format={key: site.format(key)[sample_idx] for key in site.FORMAT},
        for bases in site.gt_bases:
            if '/' in bases:
                bases = tuple(bases.split('/'))
                phased = False
            elif '|' in bases:
                bases = tuple(bases.split('|'))
                phased = True
            gt_key = bases + (phased,)
            if gt_key in seen:
                continue
            seen.add(gt_key)
            genotype_id = len(genotypes)
            gt = Genotype(
                site_id=site_id, 
                genotype_id=genotype_id,
                bases=bases,
                phased=phased,
                ref=site.REF,
            )
            genotypes[genotype_id] = gt

        site_obj = cls(
            site_id=site_id,
            chrom=site.CHROM,
            pos=site.POS,
            ref=site.REF,
            filter=site.FILTER,
            genotypes=genotypes,
            var_type=site.var_type,
            info=dict(site.INFO),
            **kw
        )
        return site_obj

    def call_site(self, site=None):
        gt_map = {tuple(gt.bases) + (gt.phased,): gt for gt in self.genotypes.values()}
        scores = []
        nc = 0
        for (idx, bases) in enumerate(site.gt_bases):
            if '/' in bases:
                bases = tuple(bases.split('/'))
                phased = False
            elif '|' in bases:
                bases = tuple(bases.split('|'))
                phased = True
            gt_key = bases + (phased,)
            gt = gt_map[gt_key]
            if not gt.log_odds:
                nc += 1
                scores.append(0)
                continue
            if gt.is_snp:
                score = gt.log_odds[0]
            else:
                assert gt.is_indel
                score = gt.log_odds[1]
            # sample level filtering (disabled)
            #if score < 0:
                #site.genotypes[idx] = ([-1] * site.ploidy) + [False]
            scores.append(score)
        site.set_format("BT", np.array(scores))
        if len(scores) != nc:
            all_fail = all([score < 0 for score in scores if score != 0])
            if all_fail:
                site.FILTER = "BERT"
            else:
                site.FILTER = "PASS"
            best_score = max([score for score in scores if score != 0])
            site.INFO["BLOD"] = best_score
        site.genotypes = site.genotypes
        return site
                
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
        ref=None,
        bases=None,
        phased=False,
        status='pending',
        log_odds=None
    ):
        self.site_id = site_id
        self.genotype_id = genotype_id
        self.ref = ref
        self.bases = bases
        self.phased = phased
        self.status = status
        self.log_odds = log_odds

    def __repr__(self):
        pchar = '|' if self.phased else '/'
        bases = str.join(pchar, self.bases)
        msg = f"{self.__class__.__name__}<#{self.site_id},{self.genotype_id}> pending={self.is_pending}, alleles={bases}, log_odds={self.log_odds}"
        return msg

    @property
    def is_pending(self):
        return self.status == "pending"

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
