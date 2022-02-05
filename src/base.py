
class BaseObject(object):
    pass

class Site(BaseObject):
    def __init__(self,
        site_id=None,
        chrom=None,
        pos=None,
        ref=None,
        ignore=False,
        genotypes=None,
        status=None
    ):
        self.site_id = site_id
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.ignore = ignore
        self.genotypes = genotypes
        self.status = status

    def load_from_site(cls, site=None, site_id=None):
        # load genotypes first
        # XXX: it might be possible to encode phase
        seen = set()
        genotypes = []
        for (bases, var_type) in zip(site.gt_bases, site.gt_types):
            bases = [tuple(re.split('[/|]', bases)) for bases in bases]
            if gt_bases in seen:
                continue
            seen.add(gt_bases)
            genotype_id = len(genotypes)
            gt = Genotype(
                site_id=site_id, 
                genotype_id=genotype_id,
                bases=gt_bases,
                var_type=var_type,
            )
            genotypes.append(gt)

        site_obj = cls(
            site_id=site_id,
            chrom=site.CHROM,
            pos=site.POS,
            genotypes=genotypes,
            var_type=site.var_type,
            ref=site.REF,
        }
        return site_obj

class Genotype(BaseObject):
    def __init__(self,
        site_id=None,
        genotype_id=None,
        var_type=None,
        ref=None,
        bases=None,
        ignore=False,
        log_odds=None
    ):
        self.site_id = site_id
        self.genotype_id = genotype_id
        self.var_type = var_type
        self.ref = ref
        self.bases = bases
        self.ignore = ignore
        self.log_odds = log_odds

    def is_homref(self):
        return self.var_type == 0

