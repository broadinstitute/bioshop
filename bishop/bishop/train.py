import argparse
from . io import IntervalLists

DefaultFields = ['AS_BaseQRankSum', 'AS_FS', 'AS_InbreedingCoeff', 'AS_MQ', 'AS_MQRankSum', 'AS_QD', 'AS_ReadPosRankSum', 'AS_SOR']

parser = argparse.ArgumentParser(description='Train allele specific classification')
parser.add_argument('--query_vcf', help='Path to VCF to call')
parser.add_argument('--target_vcf', help='Path to VCF with valid calls from population')
parser.add_argument('--reference', help='Path to FastA reference')
parser.add_argument('--assembly', default='GRCh38.p14', help='Name of the geome assembly to use')
parser.add_argument('--model_out', default='model.out', help='Path for generated model')
parser.add_argument('--skip_filtered', action='store_true', default=False, help='While building training set, skip filtered sites')
parser.add_argument('-S', '--stratification', dest='strat_intervals', action='append', help='Interval file for labeling lookup')
parser.add_argument('--unbalanced', action='store_true', default=False, help='Do not balance training set')
parser.add_argument('--eval-ratio', type=float, default=0.2, help='Ratio of dataset to reserve for evaluation, between 0 and 1')


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

def ack(query_path):
    ga = GenomeAssemblyMetadata.load('GRCh38.p14')
    region = Region('chr1', start=1_000_000, stop=2_000_000)
    overlaps = load_interval_lists(interval_files, astype='dataframe')
    target_vcf = VCF(dbsnp_url, metadata=ga, ignore_missing=True)
    query_vcf = VCF(query_path, metadata=ga, ignore_missing=True)
    flanker = VariantFlanks(assembly=ga, as_scheme='ucsc')
    cmp = ComparisonTask(
        query_vcf=query_vcf,
        target_vcf=target_vcf,
        flanker=flanker,
        overlaps=overlaps,
    )
    itr = cmp(region=region)
    itr = custom_itr(itr, annotate_func)
    df = to_dataframe(itr)
    return df

if __name__ == '__main__':
    query_path = sys.argv[1]
    ack(query_path)
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
