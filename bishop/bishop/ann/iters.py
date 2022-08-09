from .. rep.region import Region
from .. rep.fingerprint import AlleleFingerprint
from .. utils import is_concrete_nucleotides

def iter_sites(vcf=None, with_index=False, region=None):
    if isinstance(region, Region):
        region = str(region)
    sites = vcf.fetch(region=region)
    for (site_idx, site) in enumerate(sites):
        row = dict(site=site)
        if with_index:
            row.update(dict(site_idx=site_idx))
        yield row

def flank_site(itr=None, flanker=None):
    for row in itr:
        site = row['site']
        flanks = flanker.get_flanks(site=site)
        row['flanks'] = flanks
        yield row

def skip_site(itr=None, skip_filtered=True, skip_ambiguous_bases=True):
    for row in itr:
        site = row['site']
        if skip_filtered:
            filters = set(site.filter)
            if filters and ('PASS' not in filters):
                row['skip_site'] = True
                row['skip'] = True
        if skip_ambiguous_bases and 'flanks' in row:
            (up, down) = (row['flanks']['flank_up'], row['flanks']['flank_down'])
            if not is_concrete_nucleotides(up + down):
                row['skip_flanks'] = True
                row['skip'] = True
        yield row

def iter_alleles(itr=None, with_index=False, with_ref_allele=False):
    for row in itr:
        site = row['site']
        alleles = list(site.alts)
        if with_ref_allele:
            alleles.append(site.ref)
        for (allele_idx, allele) in enumerate(alleles):
            allele_row = row.copy()
            allele_row.update(dict(allele=allele))
            if with_index:
                allele_row.update(dict(allele_idx=allele_idx))
            yield allele_row

def skip_allele(itr=None, skip_ambiguous_bases=True):
    concrete_bases = set('AGTC')
    for row in itr:
        allele = row['allele']
        if skip_ambiguous_bases and \
            set(str(allele).upper()) - concrete_bases:
                row['skip_allele'] = True
                row['skip'] = True
        yield row