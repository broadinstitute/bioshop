import time
import random
from functools import partial
from pprint import pprint

from .. rep.reference import Reference
from .. rep.vcf import VCF
from .. rep.region import Region
from .. ann.flank import VariantFlanks
from .. ann.melt import OligoMelt
from .. ann.iters import *
from .. ann.fingerprint import ComparisonTask
from .. ann.classify import AnnotateCozy
from .. io.intervals import load_interval_lists
from .. io.monitor import Monitor

from pysam import VariantFile

import sys
import argparse

CLI_NAME = 'etl'
CLI_ALIASES = ()
CLI_DESC = 'Train allele specific classification'

def get_cli_parser(parent=None):
    if parent:
        parser = parent.add_parser(
            CLI_NAME, 
            aliases=CLI_ALIASES,
            description=CLI_DESC
        )
    else:
        parser = argparse.ArgumentParser(description=CLI_DESC)
    #
    parser.add_argument(
        '--query_vcf',
        dest='query_vcf_path',
        required=True,
        help='Path to VCF to call'
    )
    parser.add_argument(
        '--query_vcf_index',
        dest='query_vcf_index_path',
        required=False,
        help='Path to VCF index to call'
    )
    parser.add_argument(
        '--target_vcf',
        dest='target_vcf_path',
        required=True, 
        help='Path to VCF with valid calls from population'
    )
    parser.add_argument(
        '--target_vcf_index',
        dest='target_vcf_index_path',
        required=False, 
        help='Path to VCF index with valid calls from population'
    )
    parser.add_argument(
        '--reference',
        dest='reference_path',
        required=True,
        help='Path to genome reference'
    )
    parser.add_argument(
        '--reference_index',
        dest='reference_index_path',
        required=False,
        help='Path to index for genome reference'
    )
    parser.add_argument(
        '-o', '--output',
        dest='output_path',
        default='dataframe.pickle',
        help='Path for generated Pandas dataframe'
    )
    # XXX: do entire genome if not provided
    parser.add_argument(
        '-R', '--region',
        required=False,
        action='append', 
        type=str,
        help='Region(s) to generate results from'
    )
    parser.add_argument(
        '-I', '--intervals',
        required=False,
        action='append', 
        type=str,
        help='Interval(s) to generate results from'
    )

    #parser.add_argument('--skip_filtered', action='store_true', default=False, help='While building training set, skip filtered sites')

    parser.add_argument(
        '-S', '--stratification',
        dest='strat_intervals',
        action='append', 
        type=str,
        help='Interval file for labeling lookup'
    )
    return parser

def etl(
    query_vcf_path=None, 
    query_vcf_index_path=None, 
    target_vcf_path=None,
    target_vcf_index_path=None,
    reference_path=None,
    reference_index_path=None,
    strat_intervals=None,
    intervals=None
):

    if strat_intervals:
        overlaps = load_interval_lists(strat_intervals, astype='dataframe')
    else:
        overlaps = None

    reference = Reference(
        reference_path=reference_path, 
        reference_index_path=reference_index_path
    )
    target_vcf = VCF(target_vcf_path, index_filename=target_vcf_index_path, drop_samples=True)
    query_vcf = VCF(query_vcf_path, index_filename=query_vcf_index_path, drop_samples=True)
    flanker = VariantFlanks(reference=reference)
    melter = OligoMelt(reference=reference)
    annotate_func = AnnotateCozy()
    cmp = ComparisonTask(
        query_vcf=query_vcf,
        target_vcf=target_vcf,
        flanker=flanker,
        melter=melter,
        overlaps=overlaps,
        annotate=annotate_func,
    )
    mon = Monitor()
    mon.enable_reporting()
    with mon:
        df = cmp.compare_regions(region_list=intervals)
    return df

def main(args):
    df = etl(
        reference_path=args.reference_path,
        reference_index_path=args.reference_index_path,
        query_vcf_path=args.query_vcf_path,
        query_vcf_index_path=args.query_vcf_index_path,
        target_vcf_path=args.target_vcf_path,
        target_vcf_index_path=args.target_vcf_index_path,
        strat_intervals=args.strat_intervals,
        intervals=args.intervals,
    )
    df.to_pickle(args.output_path)

def validate_args(args):
    if len(args.intervals) < 1:
        msg = f'Currently, you must provide either a region or interval'
        raise TypeError(msg)

def main_cli(args=None):
    if args is None:
        parser = get_cli_parser()
        args = parser.parse_args()
    intn = lambda it: it.split('=') if '=' in it else (it.split('/')[-1], it)
    hdr = ('name', 'path')
    args.strat_intervals = args.strat_intervals or ()
    args.strat_intervals = [dict(zip(hdr, intn(it))) for it in args.strat_intervals]
    intervals = []
    if args.intervals:
        for interval_list in load_interval_lists(args.intervals):
            intervals += interval_list
    if args.region:
        intervals += [Region(reg) for reg in args.region]
    # XXX: flatten intervals?
    args.intervals = intervals
    validate_args(args)
    main(args)

if __name__ == "__main__":
    main_cli()
