import sys
import argparse

from sklearn.gaussian_process import *
from sklearn.preprocessing import *
from sklearn.neighbors import *
from sklearn.ensemble import *
from sklearn.neural_network import *

import pickle
import numpy as np
import pandas as pd
import sys

from .. ann.classify import Classifier

Classifiers = {
    'rf': RandomForestClassifier,
    'gb': GradientBoostingClassifier,
    'mlp': MLPClassifier,
    'kn': KNeighborsClassifier,
    'ada': AdaBoostClassifier, 
}

parser = argparse.ArgumentParser(description='Fit classifier to data')

parser.add_argument(
    '-i', '--input',
    dest='input_list',
    action='append',
    required=True,
    help='Path to one or more panda dataframes'
)
parser.add_argument(
    '-o', '--output',
    dest='classifier_path',
    type=str,
    default='classifier.pickle',
    help='Path to write classifier state',
)
parser.add_argument(
    '--classifier',
    dest='classifier',
    choices=tuple(Classifiers.keys()),
    default='rf',
    help='Classifier to use',
)
parser.add_argument(
    '--test-frac',
    dest='test_frac',
    default=.3,
    type=float,
    help='Fraction to hold in reserve for testing'
)
parser.add_argument(
    '--random-seed',
    dest='random_seed',
    default=None,
    type=int,
    help='Random seed to use (randomly assigned if left unset)'
)

def prepare_dataframe(df_list=None, random_seed=None):
    frames = (pd.read_pickle(fn) for fn in df_list)
    df = pd.concat(list(frames))

    df['label'] = df['fingerprint_match'].astype(int)
    for colname in df.columns:
        if colname.startswith('overlaps_with_'):
            df[colname] = df[colname].astype(int)

    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.fillna(0)

    df = df.sample(frac=1, random_state=random_seed)
    n_classes = df['label'].nunique()
    class_counts = df['label'].value_counts().to_list()
    min_class = np.argmin(class_counts)
    n_examples = np.min(class_counts)
    #print(len(df), n_classes, class_counts, min_class, n_examples)

    balanced_ds = []
    for class_idx in range(n_classes):
        subset = df[df['label'] == class_idx].sample(n=n_examples, random_state=random_seed)
        balanced_ds.append(subset)

    df = pd.concat(balanced_ds)
    df = df.sample(frac=1, random_state=random_seed)
    df = df.drop('site_idx', axis=1)
    df['allele_len'] = df.allele.str.len()
    df = df.drop('allele', axis=1)
    df = df.drop('allele_idx', axis=1)
    df = df.drop('chrom', axis=1)
    df = df.drop('fingerprint_match', axis=1)
    df = df.replace(dict(variant_type=dict(SNP=0, INDEL=1)))
    return df

def fit_classifier(clf_class=None, df=None, test_frac=None, random_seed=None):
    clf_name = clf_class.__name__
    kw = {'random_state': random_seed}
    if clf_name == 'RandomForestClassifier':
        kw['n_jobs'] = -1
    clf = clf_class(**kw)
    clf = Classifier(classifier=clf)
    acc = clf.fit_and_score(df=df, test_frac=test_frac)
    rpt = f"{clf_name} accuracy with test labels {acc * 100:.02f}%"
    print(rpt)
    return clf

def main(args):
    df = prepare_dataframe(
        df_list=args.input_list,
        random_seed=args.random_seed
    )
    clf_class = Classifiers[args.classifier]
    clf = fit_classifier(
        clf_class=clf_class, 
        df=df, 
        test_frac=args.test_frac, 
        random_seed=args.random_seed
    )
    clf.save_classifier(args.classifier_path)

def validate_args(args):
    pass

def main_cli():
    args = parser.parse_args()
    main(args)

if __name__ == "__main__":
    main_cli()
