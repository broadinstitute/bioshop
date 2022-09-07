import time
import random
from functools import partial
from pprint import pprint

from .. rep.assembly import GenomeAssemblyMetadata
from .. rep.vcf import VCF
from .. rep.region import Region
from .. ann.flank import VariantFlanks
from .. ann.iters import *
from .. ann.fingerprint import ComparisonTask
from .. ann.classify import AnnotateCozy
from .. io.intervals import load_interval_lists
from .. io.monitor import Monitor

from pysam import VariantFile

import sys
import argparse

parser = argparse.ArgumentParser(description='Train allele specific classification')
parser.add_argument(
    '--query_vcf',
    dest='query_vcf_path',
    required=True,
    help='Path to VCF to call'
)
parser.add_argument(
    '--target_vcf',
    dest='target_vcf_path',
    required=True, 
    help='Path to VCF with valid calls from population'
)
parser.add_argument(
    '--assembly',
    dest='assembly_name',
    default='GRCh38.p14',
    help='Name of the geome assembly to use'
)
parser.add_argument(
    '--output',
    dest='output_path',
    default='dataframe.pickle',
    help='Path for generated Pandas dataframe'
)
# XXX: add support for interval lists
# XXX: do entire genome if not provided
parser.add_argument(
    '-R', '--region',
    required=True,
    type=str,
    help='Region to generate results from'
)

#parser.add_argument('--skip_filtered', action='store_true', default=False, help='While building training set, skip filtered sites')

parser.add_argument(
    '-S', '--stratification',
    dest='strat_intervals',
    action='append', 
    type=str,
    help='Interval file for labeling lookup'
)

def etl(
    query_vcf_path=None, 
    target_vcf_path=None,
    assembly_name=None,
    strat_intervals=None,
    region=None,
    as_scheme='ucsc'
):
    ga = GenomeAssemblyMetadata.load(assembly_name)
    overlaps = load_interval_lists(strat_intervals, astype='dataframe')
    target_vcf = VCF(target_vcf_path, metadata=ga, ignore_missing=True)
    query_vcf = VCF(query_vcf_path, metadata=ga, ignore_missing=True)
    flanker = VariantFlanks(assembly=ga, as_scheme=as_scheme)
    annotate_func = AnnotateCozy()
    cmp = ComparisonTask(
        query_vcf=query_vcf,
        target_vcf=target_vcf,
        flanker=flanker,
        overlaps=overlaps,
        annotate=annotate_func,
        assembly=ga,
        as_scheme=as_scheme,
    )
    mon = Monitor()
    mon.enable_reporting(2)
    with mon:
        df = cmp.compare_region(region=region)
    return df

def main(args):
    df = etl(
        query_vcf_path=args.query_vcf_path,
        target_vcf_path=args.target_vcf_path,
        assembly_name=args.assembly_name,
        strat_intervals=args.strat_intervals,
        region=args.region,
    )
    df.to_pickle(args.output_path)

def validate_args(args):
    pass

def main_cli():
    intn = lambda it: it.split('=') if '=' in it else (it.split('/')[-1], it)
    hdr = ('name', 'path')
    args = parser.parse_args()
    args.strat_intervals = [dict(zip(hdr, intn(it))) for it in args.strat_intervals]
    args.region = Region(args.region)
    validate_args(args)
    main(args)

if __name__ == "__main__":
    main_cli()
