import os
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

from .. ann.classify import Classifier, numlint, balance_dataframe
from .. utils import concat_saved_dataframes

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
parser.add_argument(
    '--combine',
    dest='combine_models',
    default=False,
    type=bool,
    help='Generate two different models for snp/non-snp (default) or a single, combined model'
)

def fit_classifier(clf_class=None, df=None, test_frac=None, random_seed=None):
    clf_name = clf_class.__name__
    kw = {'random_state': random_seed}
    if clf_name == 'RandomForestClassifier':
        kw['n_jobs'] = -1
    clf = clf_class(**kw)
    clf = Classifier(classifier=clf)
    acc = clf.fit_and_score(df=df, test_frac=test_frac)
    return clf

def create_combined_models(args):
    df = concat_saved_dataframes(args.input_list)
    df = numlint(df)
    df = balance_dataframe(df=df, random_seed=args.random_seed)
    clf_class = Classifiers[args.classifier]
    clf = fit_classifier(
        clf_class=clf_class, 
        df=df, 
        test_frac=args.test_frac, 
        random_seed=args.random_seed
    )
    clf.save_classifier(args.classifier_path)

def create_seperate_models(args):
    (root, cfn) = os.path.split(args.classifier_path)
    df = concat_saved_dataframes(args.input_list)
    df = numlint(df)
    # SNP
    df_snp = df[df.feature_is_snp == True].drop(columns='feature_is_snp')
    df_snp = balance_dataframe(df=df_snp, random_seed=args.random_seed)
    clf_class = Classifiers[args.classifier]
    print('Fitting SNP model')
    snp_clf = fit_classifier(
        clf_class=clf_class, 
        df=df_snp, 
        test_frac=args.test_frac, 
        random_seed=args.random_seed
    )
    snp_cfn = os.path.join(root, f'SNP_{cfn}')
    snp_clf.save_classifier(snp_cfn)
    print('=' * 50)
    # INDEL
    df_indel = df[df.feature_is_snp == False].drop(columns='feature_is_snp')
    df_indel = balance_dataframe(df=df_indel, random_seed=args.random_seed)
    clf_class = Classifiers[args.classifier]
    print('Fitting INDEL model')
    indel_clf = fit_classifier(
        clf_class=clf_class, 
        df=df_indel, 
        test_frac=args.test_frac, 
        random_seed=args.random_seed
    )
    indel_cfn = os.path.join(root, f'INDEL_{cfn}')
    indel_clf.save_classifier(indel_cfn)
    print('=' * 50)

def main(args):
    if args.combine_models:
        create_combined_models(args)
    else:
        create_seperate_models(args)

def validate_args(args):
    pass

def main_cli():
    args = parser.parse_args()
    main(args)

if __name__ == "__main__":
    main_cli()
