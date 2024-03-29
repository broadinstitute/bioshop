import pandas as pd

from . precis import Precis

from .. rep.region import Region
from .. rep.fingerprint import AlleleFingerprint
from .. utils import is_concrete_nucleotides

def iter_sites(vcf=None, region=None, assembly=None, as_scheme=None):
    if not isinstance(region, Region):
        region = Region(region)
    sites = vcf.fetch(region=str(region))
    ff_offset = 0
    for (site_idx, site) in enumerate(sites):
        if site.pos < region.start:
            ff_offset += 1
            continue
        if site.pos > region.stop:
            break
        row = Precis()
        row.cache.site = site
        row.meta.site_idx = (site_idx - ff_offset)
        row.meta.pos = site.pos
        row.meta.ref = site.ref
        if assembly and as_scheme:
            row.meta.chrom = assembly.as_scheme(site.chrom, as_scheme=as_scheme)
        else:
            row.meta.chrom = site.chrom
        yield row

def pos_monitor(itr=None, remote=None, name='nucs'):
    last_pos = None
    for row in itr:
        if last_pos is not None:
            delta = row.meta.pos - last_pos
            remote.add(name=name, value=delta)
        last_pos = row.meta.pos
        yield row

def iter_monitor(itr=None, remote=None, name=None):
    for row in itr:
        remote.add(name=name, value=1)
        yield row

def flank_site(itr=None, flanker=None):
    for row in itr:
        flanks = flanker.get_flanks(site=row.cache.site)
        row.cache.flanks = flanks
        yield row

def overlaps_with_site(itr=None, overlaps=None, slop=5):
    for row in itr:
        if not row.filter:
            # XXX: index match? does slop take care of this?
            start = row.meta.pos - slop
            stop = row.meta.pos + len(row.meta.ref) + slop
            region = Region(row.meta.chrom, start, stop)
            hits = set(overlaps.overlaps_with(region))
            ovmap = {f'overlap_{nm}': int(nm in hits) for nm in overlaps.names}
            row.feature.update(ovmap)
        yield row

def filter_by_site(itr=None, skip_filtered=False, skip_ambiguous_bases=True):
    for row in itr:
        site = row.cache.site
        if skip_filtered:
            filters = set(site.filter)
            if filters and ('PASS' not in filters):
                row.filter.set_filter('site filtered by VCF')
        if skip_ambiguous_bases and 'flanks' in row.cache:
            flanks = row.cache.flanks[0] + row.cache.flanks[1]
            if not is_concrete_nucleotides(flanks):
                row.filter.set_filter('ambiguous base in genomic flanks')
        yield row

def iter_alleles(itr=None, with_ref_allele=False):
    for row in itr:
        site = row.cache.site
        if site.alts is None:
            alleles = list()
        else:
            alleles = list(site.alts)
        if with_ref_allele:
            alleles.append(site.ref)
        if alleles:
            for (allele_idx, allele) in enumerate(alleles):
                allele_row = row.copy()
                allele_row.meta.allele = allele
                allele_row.meta.allele_idx = allele_idx
                yield allele_row
        else:
            yield row

def filter_by_allele(itr=None, skip_ambiguous_bases=True):
    concrete_bases = set('AGTC')
    for row in itr:
        if not row.filter:
            if 'allele' not in row.meta:
                row.filter.set_filter('missing allele')
            #
            elif skip_ambiguous_bases and \
                set(str(row.meta.allele).upper()) - concrete_bases:
                    row.filter.set_filter('allele is symbolic')
        yield row

def custom_itr(itr=None, custom_func=None):
    for row in itr:
        row = custom_func(row)
        yield row

def to_dataframe(itr, include_domains=('meta', 'feature', 'label', 'filter')):
    rows = []
    for row in itr:
        if row.filter:
            continue
        row = row.flatten(include_domains=include_domains)
        rows.append(row)
    if not rows:
        return
    return pd.DataFrame(rows)

def annotate_alleles_from_dataframe(itr=None, df=None, columns=None):
    # XXX
    if columns is None:
        columns = [('score', 'AS_BLOD')]
    # site-level
    for row in itr:
        site = row.cache.site
        if not row.filter:
            #hits = df[df.meta_site_idx == row.meta.site_idx]
            hits = df[df.meta_pos == row.meta.pos]
            hits = hits.sort_values(['meta_allele_idx'], ascending=True)
            assert len(hits) == len(site.alts)
            for col in columns:
                if len(col) == 1:
                    (from_col, to_col) = (col, col)
                elif len(col) == 2:
                    (from_col, to_col) = col
                else:
                    raise ValueError(col)
                site.info[to_col] = tuple(hits[from_col].to_list())
        yield row
