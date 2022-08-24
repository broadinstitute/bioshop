import time
import random
from pprint import pprint

from .. rep.assembly import GenomeAssemblyMetadata
from .. rep.vcf import VCF
from .. rep.region import Region
from .. ann.flank import VariantFlanks
from .. ann.iters import *
from .. ann.fingerprint import ComparisonTask
from .. io.intervals import load_interval_lists

from pysam import VariantFile

import sys

field_names = ['AS_BaseQRankSum', 'AS_FS', 'AS_InbreedingCoeff', 'AS_MQ', 'AS_MQRankSum', 'AS_QD', 'AS_ReadPosRankSum', 'AS_SOR']

def annotate_func(row):
    if 'skip' not in row:
        site = row['site']
        al_idx = row['allele_idx']
        allele = row['allele']
        if len(site.ref) == len(allele):
            row['variant_type'] = 'SNP'
        else:
            row['variant_type'] = 'INDEL'
        info = {fn: site.info[fn][al_idx] for fn in field_names}
        row.update(info)
    return row

def call_allele(itr=None, model=None):
    pass

def build_iter(vcf=None, region=None, overlaps=None, model=None):
    itr = iter_sites(vcf=vcf, with_index=True, region=region)
    itr = flank_site(itr=itr, flanker=flanker)
    if overlaps is not None:
        itr = overlaps_with_site(itr, overlaps=overlaps)
    itr = skip_site(itr=itr)
    itr = iter_alleles(itr=itr, with_index=True)
    itr = skip_allele(itr=itr)
    return fingerprint_allele(itr=itr)

def ack(query_path, region):
    region = Region(region)
    query_vcf = VCF(query_path, metadata=ga, ignore_missing=True)
    overlaps = load_interval_lists(interval_files, astype='dataframe')
    cmp = ComparisonTask(
        query_vcf=query_vcf,
        target_vcf=target_vcf,
        flanker=flanker,
        overlaps=overlaps,
        annotate=annotate_func
    )
    print('running compare')
    itr = cmp.compare_region(region=region)
    #itr = custom_itr(itr, annotate_func)
    df = to_dataframe(itr)
    df.to_pickle(f'df-{region}.pickle')
    return df

if __name__ == '__main__':
    call_vcf = sys.argv[1]
    region = sys.argv[2]
    ack(call_vcf, target_path, region)
