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
from .. io.intervals import load_interval_lists
from .. ann.classify import Classifier, AnnotateCozy, ClassifyTask
from .. io.monitor import Monitor

from pysam import VariantFile

import sys
import argparse

CLI_NAME = 'call'
CLI_ALIASES = ()
CLI_DESC = 'Call VCF sites based on allelic training data'

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
        '--classifier',
        dest='classifier_path',
        default='classifier.pickle',
        help='Path to pickled classifier'
    )
    parser.add_argument(
        '--assembly',
        dest='assembly_name',
        default='GRCh38.p14',
        help='Name of the geome assembly to use'
    )
    parser.add_argument(
        '-o',
        '--output',
        dest='output_vcf_path',
        default='called.vcf',
        help='Path to generated VCF',
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


def call(
    query_vcf_path=None, 
    output_vcf_path=None, 
    classifier_path=None,
    assembly_name=None,
    strat_intervals=None,
    intervals=None,
    as_scheme='ucsc'
):
    ga = GenomeAssemblyMetadata.load(assembly_name)
    overlaps = load_interval_lists(strat_intervals, astype='dataframe')
    # XXX: drop samples?
    query_vcf = VCF(query_vcf_path, metadata=ga, ignore_missing=True)
    specs = [
        {'ID': 'BLOD', 'Description': 'Bishop LOD', 'Type': 'Float', 'Number': 1},
        {'ID': 'AS_BLOD', 'Description': 'Allele Specific Bishop LOD', 'Type': 'Float', 'Number': 'A'},
    ]
    order = ('ID', 'Number', 'Type', 'Description')
    for spec in specs:
        items = [(key, spec[key]) for key in order if key in spec]
        query_vcf.header.add_meta(key='INFO', items=items)
    output_vcf = query_vcf.to_writer(output_vcf_path)
    annotate_func = AnnotateCozy()
    # XXX: support pandas dataframe saving as well?
    cls = ClassifyTask(
        query_vcf=query_vcf,
        classifier_path=classifier_path,
        overlaps=overlaps,
        annotate=annotate_func,
        assembly=ga,
        as_scheme=as_scheme,
    )
    for region in intervals:
        cls.call_vcf_sites(output_vcf=output_vcf, region=region)

def main(args):
    mon = Monitor()
    mon.enable_reporting()
    with mon:
        call(
            query_vcf_path=args.query_vcf_path,
            output_vcf_path=args.output_vcf_path,
            classifier_path=args.classifier_path,
            assembly_name=args.assembly_name,
            strat_intervals=args.strat_intervals,
            intervals=args.intervals,
        )

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
    args.strat_intervals = [dict(zip(hdr, intn(it))) for it in args.strat_intervals]
    intervals = []
    if args.intervals:
        for interval_list in load_interval_lists(args.intervals):
            intervals += interval_list
    if args.region:
        intervals += [Region(reg) for reg in region_list]
    # XXX: flatten intervals?
    args.intervals = intervals
    validate_args(args)
    main(args)

if __name__ == "__main__":
    main_cli()
