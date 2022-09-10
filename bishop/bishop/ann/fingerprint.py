import os
from shutil import rmtree
import multiprocessing as mp

from . iters import *
from .. rep.region import Region
from .. rep.fingerprint import AlleleFingerprint
from .. io.monitor import get_remote_monitor

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

def fingerprint_vcf(vcf=None, region=None, flanker=None, overlaps=None, slop=0, assembly=None, as_scheme=None, remote=None):
    if slop > 0:
        region.start = max(1, region.start - slop)
        region.stop = region.stop + slop
    itr = iter_sites(vcf=vcf, region=region, assembly=assembly, as_scheme=as_scheme)
    if remote:
        itr = pos_monitor(itr, remote)
        itr = iter_monitor(itr, remote, 'sites')
    itr = flank_site(itr=itr, flanker=flanker)
    itr = filter_by_site(itr=itr)
    if overlaps is not None:
        itr = overlaps_with_site(itr, overlaps=overlaps)
    itr = iter_alleles(itr=itr) 
    itr = filter_by_allele(itr=itr)
    itr = fingerprint_allele(itr=itr)
    if remote:
        itr = iter_monitor(itr, remote, 'alleles')
    return itr

def fingerprint_and_index_vcf(vcf=None, region=None, flanker=None, assembly=None, as_scheme=None, remote=None, slop=None):
    itr = fingerprint_vcf(
        vcf=vcf, region=region, flanker=flanker, 
        assembly=assembly, as_scheme=as_scheme, 
        slop=slop, remote=remote
    )
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
    def __init__(self, query_vcf=None, target_vcf=None, flanker=None, overlaps=None, annotate=None, slop=50, assembly=None, as_scheme=None):
        self.query_vcf = query_vcf
        self.target_vcf = target_vcf
        self.flanker = flanker
        self.overlaps = overlaps
        self.annotate = annotate
        self.slop = slop
        self.assembly = assembly
        self.as_scheme = as_scheme
        self.index_remote = get_remote_monitor(domain='IDX')
        self.fingerprint_remote = get_remote_monitor(domain='FP')
    
    def __call__(self, region=None):
        target_prints = fingerprint_and_index_vcf(
                vcf=self.target_vcf, region=region,
                flanker=self.flanker, slop=self.slop,
                assembly=self.assembly, as_scheme=self.as_scheme,
                remote=self.index_remote
        )
        query_prints = fingerprint_vcf(
            vcf=self.query_vcf, region=region, flanker=self.flanker, 
            overlaps=self.overlaps, assembly=self.assembly, 
            as_scheme=self.as_scheme, remote=self.fingerprint_remote
        )
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
        df = to_dataframe(list(batch))
        return (region, df)

    def compare_region(self, region=None, chunk_size=1_000_000):
        if not isinstance(region, Region):
            region = Region(region)
        if len(region) == 0:
            vcf_contig = self.query_vcf.header.contigs[region.contig]
            region = region.clone(start=1, stop=vcf_contig.length)
        # XXX: set at construction
        #n_threads = mp.cpu_count()
        #regions = region.shard(n_threads)
        regions = region.split(chunk_size)
        regions = list(regions)
        df_list = []
        with mp.Pool() as pool:
            itr = pool.imap(self.batch_call, regions)
            for (cur_region, df) in itr:
                df_list.append(df)
        # XXX: single threaded debug
        #itr = map(self.batch_call, regions)
        #df_list = [df for (_, df) in itr]
        df = pd.concat(df_list)
        return df
