import os
from cyvcf2 import VCF
from . utils import vcf_progress_bar

import pandas as pd
import numpy as np

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

def split_dataset(all_examples=None, train_frac=0.9, seed=None):
    rng = np.random.default_rng(seed=seed)
    label_counts = all_examples['label'].value_counts()
    n_choose = label_counts.min()
    print(f"n_choose={n_choose}")
    print(label_counts)
    sampled_examples = []

    for label_name in all_examples['label'].unique():
        example_set = all_examples \
                    [all_examples.label == label_name] \
                    .sample(n_choose, random_state=rng.bit_generator)
        sampled_examples.append(example_set)

    sampled_examples = pd.concat(sampled_examples)
    sampled_examples = sampled_examples.sample(frac=1, random_state=rng.bit_generator)
    print(f"{len(sampled_examples)} total examples")

    train_cnt = int(round(len(sampled_examples) * train_frac))
    train_mask = np.arange(len(sampled_examples)) > train_cnt
    train_ds = sampled_examples[train_mask]
    test_ds = sampled_examples[~train_mask]
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

def build_training_dataset(
        vcf_sample_list=None, 
        output_dir='mostly-training-data', 
        method='vcfeval', 
        seed=None,
        train_frac=0.9,
    ):

    assert method == 'vcfeval'
    os.makedirs(output_dir, exist_ok=True)
    pile_fn = os.path.join(output_dir, "all-unbalanced.tsv")
    train_fn = os.path.join(output_dir, "train.tsv")
    test_fn = os.path.join(output_dir, "eval.tsv")

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

    print("Balancing dataset")
    df = pd.read_csv(pile_fn, sep='\t')
    (train_df, test_df) = split_dataset(df, train_frac=train_frac, seed=seed)

    print("Writing balanced+split datasets to disk")
    train_df.to_csv(train_fn, sep='\t', index=False)
    test_df.to_csv(test_fn, sep='\t', index=False)
