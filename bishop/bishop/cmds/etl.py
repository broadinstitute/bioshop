import time
import random
from pprint import pprint

from bishop.rep.assembly import GenomeAssemblyMetadata
from bishop.rep.vcf import VCF
from bishop.rep.region import Region
from bishop.ann.flank import VariantFlanks
from bishop.ann.iters import *
from bishop.ann.fingerprint import ComparisonTask
from bishop.io.intervals import load_interval_lists

from pysam import VariantFile

import sys


def random_regions(cnt):
    for x in range(cnt):
        x = random.randint(0, 200_000_000)
        y = x + random.randint(0, 100)
        yield Region('chr1', x, y)

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

def ack(query_path, target_path, region):
    region = Region(region)
    ga = GenomeAssemblyMetadata.load('GRCh38.p14')
    overlaps = load_interval_lists(interval_files, astype='dataframe')
    target_vcf = VCF(target_path, metadata=ga, ignore_missing=True)
    query_vcf = VCF(query_path, metadata=ga, ignore_missing=True)
    flanker = VariantFlanks(assembly=ga, as_scheme='ucsc')
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
    query_path = sys.argv[1]
    target_path = sys.argv[2]
    region = 'chr1:50000000-100000000'
    region = 'chr2:50000000-100000000'
    region = 'chr3:50000000-100000000'
    ack(query_path, target_path, region)
