import argparse
from . io import IntervalLists

parser = argparse.ArgumentParser(description='Train allele specific classification')
parser.add_argument('--query_vcf', help='Path to VCF to call')
parser.add_argument('--target_vcf', help='Path to VCF with valid calls from population')
parser.add_argument('--reference', help='Path to FastA reference')
parser.add_argument('--model_out', default='model.out', help='Path for generated model')
parser.add_argument('--skip_filtered', action='store_true', default=False, help='While building training set, skip filtered sites')
parser.add_argument('-I', '--interval', dest='intervals', action='append', help='Interval file for training examples')
parser.add_argument('--unbalanced', action='store_true', default=False, help='Do not balance training set')
parser.add_argument('--eval-ratio', type=float, default=0.2, help='Ratio of dataset to reserve for evaluation, between 0 and 1')

def main(args):
    intervals = IntervalLists.load(args.intervals)

def validate_args(args):
    pass

def main_cli():
    args = parser.parse_args()
    validate_args(args)
    main(args)

if __name__ == __main__:
    main_cli()
