import os
from cyvcf2 import VCF
from . utils import vcf_progress_bar

import pandas as pd
import numpy as np

DefaultAnnotations = ["DP", "QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR"]

class TrainingVcfLoader(object):
    Header = ('chrom', 'pos', 'ref', 'alt', 'var_type', 'label')

    def iter_vcf(self, vcf=None, region=None):
        pass

class TrainingVcfEvalLoader(TrainingVcfLoader):
    def __init__(self, vcfeval_call_annot='CALL'):
        self.vcfeval_call_annot = vcfeval_call_annot
    
    def get_label(self, site=None):
        if site.is_snp:
            var_type = "SNP"
        elif site.is_indel:
            var_type = "INDEL"
        else:
            return None
        truth = site.INFO.get(self.vcfeval_call_annot)
        if truth not in ('FP', 'TP'):
            return None
        label_name = f"{var_type}_{truth}"
        label_info = {
            "var_type": var_type,
            "label": label_name,
        }
        return label_info

    def process_site(self, site=None):
        label = self.get_label(site=site)
        if not label:
            return None
        ret = {
            "chrom": site.CHROM,
            "pos": site.POS,
            "ref": site.REF,
            "alt": str.join(',', site.ALT),
        }
        ret.update(label)
        return ret

    def iter_vcf(self, vcf=None, region=None):
        #sample_idx = vcf.samples.index(self.vcfeval_sample_name)
        itr = vcf(region) if region is not None else iter(vcf)
        for site in itr:
            ret = self.process_site(site=site)
            if ret is None:
                continue
            yield ret

def balance_dataset(all_examples=None, seed=None):
    rng = np.random.default_rng(seed=seed)
    label_counts = all_examples['label'].value_counts()
    n_choose = label_counts.min()
    print(f"n_choose={n_choose}")
    print(label_counts)
    balanced_examples = []

    for label_name in all_examples['label'].unique():
        example_set = all_examples \
                    [all_examples.label == label_name] \
                    .sample(n_choose, random_state=rng.bit_generator)
        balanced_examples.append(example_set)

    balanced_examples = pd.concat(balanced_examples)
    print(f"{len(balanced_examples)} balanced examples")
    return balanced_examples

def annotate_dataset(balanced_examples=None, vcf_src=None):
    def rename_chrom(val):
        if val.startswith('chr'):
            val = val[3:]
        try: 
            #val = f'{val:02d}'
            val = int(val)
        except ValueError:
            pass
        return val
    balanced_examples['chrom_'] = balanced_examples['chrom'].map(rename_chrom)
    sorted_examples = balanced_examples.sort_values(['chrom_', 'pos', 'sample_name'], ascending=(True, True, True))
    pbar = vcf_progress_bar(vcf=vcf_src)
    src_itr = iter(vcf_src)
    src_site = None
    rows = []
    sample_map = {vcf_src.samples[idx]: idx for idx in range(len(vcf_src.samples))}
    for row in sorted_examples.to_dict(orient="records"):
        del row['chrom_']
        pos = row['pos']
        chrom = row['chrom']
        if (
            (src_site is None) or \
            (chrom != src_site.CHROM) or \
            ((pos - src_site.POS) > 100000)
        ):
            coord = f'{chrom}:{pos}'
            src_itr = vcf_src(coord)
            src_site = next(src_itr)
        if src_site is None or chrom != src_site.CHROM:
            src_itr = vcf_src(chrom)
            src_site = next(src_itr)
        assert src_site.CHROM == chrom
        while src_site.POS < pos:
            src_site = next(src_itr)
        assert src_site.CHROM == chrom
        assert src_site.POS == pos
        info = dict(src_site.INFO)
        info = {f"INFO_{key}": info[key] for key in info}
        smp_idx = sample_map[row['sample_name']]
        bases = src_site.gt_bases[smp_idx]
        if '/' in bases:
            row['gt_bases'] = bases.split('/')
            row['phased'] = False
        elif '|' in bases:
            row['gt_bases'] = bases.split('|')
            row['phased'] = True
        row.update(info)
        rows.append(row)
        pbar(chrom, pos)
    pbar()
    return pd.DataFrame(rows)

def split_dataset(balanced_examples=None, train_frac=0.9, seed=None):
    rng = np.random.default_rng(seed=seed)
    balanced_examples = balanced_examples.sample(frac=1, random_state=rng.bit_generator)

    train_cnt = int(round(len(balanced_examples) * train_frac))
    train_mask = np.arange(len(balanced_examples)) < train_cnt
    train_ds = balanced_examples[train_mask]
    test_ds = balanced_examples[~train_mask]
    print(f"{len(train_ds)} train examples, {len(test_ds)} test examples")
    return (train_ds, test_ds)

def load_vcfeval(vcf_path=None, sample_name=None, region=None):
    vcf = VCF(vcf_path)
    pbar = vcf_progress_bar(vcf=vcf)
    vcfeval_etl = TrainingVcfEvalLoader()
    itr = vcfeval_etl.iter_vcf(vcf=vcf, region=region)

    print(f"Scanning {os.path.split(vcf_path)[-1]} for training examples")
    for site in itr:
        pbar(site['chrom'], site['pos'])
        site['sample_name'] = sample_name
        yield site
    pbar()
    print()

def pile_vcfeval(vcf_sample_list=None, pile_fn=None):
    header = TrainingVcfEvalLoader.Header + ('sample_name', )
    with open(pile_fn, 'w') as fh:
        hdr = str.join('\t', header) + '\n'
        fh.write(hdr)
        for vcf_info in vcf_sample_list:
            itr = load_vcfeval(**vcf_info)
            for row in itr:
                row = [str(row[key]) for key in header]
                row = str.join('\t', row) + '\n'
                fh.write(row)


def build_training_dataset(
        vcf_sample_list=None, 
        vcf_src_path=None,
        output_dir='mostly-training-data', 
        method='vcfeval', 
        seed=None,
        train_frac=0.9,
    ):

    assert method == 'vcfeval'
    os.makedirs(output_dir, exist_ok=True)
    pile_fn = os.path.join(output_dir, "all-unbalanced.tsv")
    balanced_fn = os.path.join(output_dir, "balanced.tsv")
    train_fn = os.path.join(output_dir, "_train.tsv")
    test_fn = os.path.join(output_dir, "_eval.tsv")

    if not os.path.exists(pile_fn):
        pile_vcfeval(vcf_sample_list=vcf_sample_list, pile_fn=pile_fn)
    else:
        print("Using previously built example pile.")

    if not os.path.exists(balanced_fn):
        print("Balancing dataset")
        df = pd.read_csv(pile_fn, sep='\t')
        df = balance_dataset(df, seed=seed)
        if vcf_src_path:
            print("Annotating balanced dataset")
            vcf_src = VCF(vcf_src_path)
            df = annotate_dataset(df, vcf_src=vcf_src)
        df.to_csv(balanced_fn, sep='\t', index=False)
    else:
        print("Using previously built balanced dataset.")
        df = pd.read_csv(balanced_fn, sep='\t')

    print("Splitting dataset")
    (train_df, test_df) = split_dataset(df, train_frac=train_frac, seed=seed)

    train_df.to_csv(train_fn, sep='\t', index=False)
    test_df.to_csv(test_fn, sep='\t', index=False)
