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
    '--split',
    dest='split_ratio',
    default=.3,
    type=float,
    help='Split ratio between train and eval',
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

def split_dataset(df=None, split_ratio=0.3, random_seed=None):
    training_cols = [col for col in df.columns if col != 'label']
    train_ds = df[training_cols].to_numpy()
    #XXX: options for pre-procesing?
    #train_ds = StandardScaler().fit_transform(train_ds)
    labels_ds = df['label'].to_numpy()
    test_count = int(round(len(train_ds) * split_ratio))
    (train_input, train_labels) = (train_ds[:test_count], labels_ds[:test_count])
    (test_input, test_labels) = (train_ds[test_count:], labels_ds[test_count:])
    rpt = [
        f'train set N={np.sum(train_labels)}, {np.sum(train_labels) / len(train_labels) * 100:.02f}% positive',
        f'test set N={np.sum(test_labels)}, {np.sum(test_labels) / len(test_labels) * 100:.02f}% positive'
    ]
    print(str.join('\n', rpt))
    return (train_input, train_labels, test_input, test_labels)

def fit_classifier(clf_class=None, df=None, split_ratio=None, random_seed=None):
    (train_input, train_labels, test_input, test_labels) = \
        split_dataset(df, random_seed=random_seed)
    clf = clf_class(random_state=random_seed)
    clf_name = clf_class.__name__
    clf.fit(train_input, train_labels)
    acc = clf.score(test_input, test_labels)
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
        split_ratio=args.split_ratio, 
        random_seed=args.random_seed
    )
    with open(args.classifier_path, 'wb') as fh:
        pickle.dump(clf, fh)

def validate_args(args):
    pass

def main_cli():
    args = parser.parse_args()
    main(args)

if __name__ == "__main__":
    main_cli()
