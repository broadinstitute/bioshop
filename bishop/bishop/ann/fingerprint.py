import multiprocessing as mp

from . iters import *
from .. rep.region import Region
from .. rep.fingerprint import AlleleFingerprint

def build_allele_index(itr):
    fingerprints = {}
    for row in itr:
        if 'skip' in row:
            continue
        allele_fp = row['allele_fingerprint']
        fingerprints[allele_fp.literal_fingerprint] = allele_fp
        fingerprints[allele_fp.cigar_fingerprint] = allele_fp
    return fingerprints

def fingerprint_allele(itr):
    for row in itr:
        if 'skip' not in row:
            site = row['site']
            flanks = tuple(row['flanks'].values())
            chrom = row.get('chrom')
            alt = row['allele']
            row['allele_fingerprint'] = AlleleFingerprint.from_site(
                site=site, alt=alt, flanks=flanks, chrom=chrom
            )
        yield row

def fingerprint_vcf(vcf=None, region=None, flanker=None):
    itr = iter_sites(vcf=vcf, with_index=True, region=region)
    itr = flank_site(itr=itr, flanker=flanker)
    itr = skip_site(itr=itr)
    itr = iter_alleles(itr=itr, with_index=True)
    itr = skip_allele(itr=itr)
    return fingerprint_allele(itr=itr)

def fingerprint_and_index_vcf(vcf=None, region=None, flanker=None):
    itr = fingerprint_vcf(vcf=vcf, region=region, flanker=flanker)
    fingerprints = build_allele_index(itr)
    return AlleleIndex(region=region, fingerprints=fingerprints)

class IndexTask:
    def __init__(self, vcf=None, flanker=None):
        self.vcf = vcf
        self.flanker = flanker

    def __call__(self, region=None):
        return fingerprint_and_index_vcf(vcf=self.vcf, region=region, flanker=self.flanker)

    def fingerprint(self, region=None, chunk_size=10_000):
        if not isinstance(region, Region):
            region = Region(region)
        regions = region.split(chunk_size)
        pool = mp.Pool()
        yield from pool.imap(self, regions)

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
    def __init__(self, query_vcf=None, target_vcf=None, flanker=None):
        self.query_vcf = query_vcf
        self.target_vcf = target_vcf
        self.flanker = flanker
        self.index_task = IndexTask(self.target_vcf, flanker=self.flanker)

    def __call__(self, region=None, slop=50):
        cache = []
        target_prints = self.index_task.fingerprint(region=region)
        query_prints = fingerprint_vcf(vcf=self.query_vcf, region=region, flanker=self.flanker)
        last_pos = 0
        for row in query_prints:
            if 'skip' in row:
                continue
            site = row['site']
            while site.pos > (last_pos - slop):
                if target_prints is None:
                    break
                try:
                    prints = next(target_prints)
                except StopIteration:
                    target_prints = None
                    break
                last_pos = prints.region.interval.upper
                cache = cache[-2:] + [prints]
            alfp = row['allele_fingerprint']
            for index in cache:
                if index.match(alfp):
                    row['matched'] = True
                    break
            yield row


