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

from .. ann.classify import Classifier, prepare_dataframe, balance_dataframe
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

def fit_classifier(clf_class=None, df=None, test_frac=None, random_seed=None):
    clf_name = clf_class.__name__
    kw = {'random_state': random_seed}
    if clf_name == 'RandomForestClassifier':
        kw['n_jobs'] = -1
    clf = clf_class(**kw)
    clf = Classifier(classifier=clf)
    acc = clf.fit_and_score(df=df, test_frac=test_frac)
    return clf

def main(args):
    df = concat_saved_dataframes(args.input_list)
    df = prepare_dataframe(df=df)
    df = balance_dataframe(df=df, random_seed=args.random_seed)
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
