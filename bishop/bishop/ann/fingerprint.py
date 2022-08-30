import multiprocessing as mp
import time

from . iters import *
from .. rep.region import Region
from .. rep.fingerprint import AlleleFingerprint
from .. utils import region_progress_bar

def build_allele_index(itr):
    fingerprints = {}
    for row in itr:
        if row.filter:
            continue
        allele_fp = row.cache.allele_fingerprint
        fingerprints[allele_fp.literal_fingerprint] = allele_fp
        fingerprints[allele_fp.cigar_fingerprint] = allele_fp
    return fingerprints

def fingerprint_allele(itr):
    for row in itr:
        if not row.filter:
            row.cache.allele_fingerprint = AlleleFingerprint.from_site(
                site=row.cache.site, 
                alt=row.meta.allele, 
                flanks=row.cache.flanks, 
                chrom=row.meta.chrom
            )
        yield row

def fingerprint_vcf(vcf=None, region=None, flanker=None, overlaps=None, slop=50, assembly=None, as_scheme=None):
    if slop > 0:
        region.start = max(0, region.start - slop)
        region.stop = region.stop + slop
    itr = iter_sites(vcf=vcf, region=region, assembly=assembly, as_scheme=as_scheme)
    itr = flank_site(itr=itr, flanker=flanker)
    itr = filter_by_site(itr=itr)
    if overlaps is not None:
        itr = overlaps_with_site(itr, overlaps=overlaps)
    itr = iter_alleles(itr=itr) 
    itr = filter_by_allele(itr=itr)
    return fingerprint_allele(itr=itr)

def fingerprint_and_index_vcf(vcf=None, region=None, flanker=None, assembly=None, as_scheme=None):
    itr = fingerprint_vcf(vcf=vcf, region=region, flanker=flanker, assembly=assembly, as_scheme=as_scheme)
    fingerprints = build_allele_index(itr)
    return AlleleIndex(region=region, fingerprints=fingerprints)

class AlleleIndex(object):
    def __init__(self, region=None, fingerprints=None):
        self.region = region
        self.fingerprints = fingerprints

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

class ComparisonTask:
    def __init__(self, query_vcf=None, target_vcf=None, flanker=None, overlaps=None, annotate=None, slop=50, assembly=None, as_scheme=None, progress_bar=True):
        self.query_vcf = query_vcf
        self.target_vcf = target_vcf
        self.flanker = flanker
        self.overlaps = overlaps
        self.annotate = annotate
        self.slop = slop
        self.assembly = assembly
        self.as_scheme = as_scheme
        self.progress_bar = progress_bar
    
    def __call__(self, region=None):
        target_prints = fingerprint_and_index_vcf(vcf=self.target_vcf, region=region, flanker=self.flanker, assembly=self.assembly, as_scheme=self.as_scheme)
        query_prints = fingerprint_vcf(vcf=self.query_vcf, region=region, flanker=self.flanker, overlaps=self.overlaps, slop=self.slop, assembly=self.assembly, as_scheme=self.as_scheme)
        if self.annotate is not None:
            query_prints = custom_itr(query_prints, self.annotate)
        for row in query_prints:
            if not row.filter:
                row.label.fingerprint_match = \
                    int(target_prints.match(row.cache.allele_fingerprint))
            # not pickle-able
            del row.cache
            yield row

    def batch_call(self, region=None, **kw):
        batch = self(region=region, **kw)
        return (region, list(batch))

    def compare_region(self, region=None, chunk_size=100_000):
        if not isinstance(region, Region):
            region = Region(region)
        if self.progress_bar:
            pbar = region_progress_bar(region=region)
        else:
            pbar = None
        regions = region.split(chunk_size)
        all_rows = []
        with mp.Pool() as pool:
            itr = pool.imap(self.batch_call, regions)
            for (cur_region, rows) in itr:
                if pbar:
                    pbar(pos=cur_region.stop)
                all_rows.extend(rows)

        return all_rows

    def compare_regions_simple(self, region=None, chunk_size=10_000, **kw):
        if not isinstance(region, Region):
            region = Region(region)
        if self.progress_bar:
            pbar = region_progress_bar(region=region)
        else:
            pbar = None
        regions = region.split(chunk_size)
        all_rows = []
        for reg in regions:
            all_rows.extend(self.batch_call(reg))
        return all_rows
