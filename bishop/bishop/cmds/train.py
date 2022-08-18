from sklearn.gaussian_process import *
from sklearn.preprocessing import *
from sklearn.neighbors import *
from sklearn.ensemble import *
from sklearn.neural_network import *
import pickle
import numpy as np
import pandas as pd
import sys

RANDOM_SEED = 0
N_CLASSES = 2
split_ratio = .3

dffn = sys.argv[1]
others = sys.argv[2:]
frames = (pd.read_pickle(fn) for fn in sys.argv[1:])
df = pd.concat(list(frames))

df['label'] = df['fingerprint_match'].astype(int)
for colname in df.columns:
    if colname.startswith('overlaps_with_'):
        df[colname] = df[colname].astype(int)

df = df.replace([np.inf, -np.inf], np.nan)
df = df.fillna(0)

df = df.sample(frac=1, random_state=RANDOM_SEED)
class_counts = df['label'].value_counts().to_list()
min_class = np.argmin(class_counts)
n_examples = np.min(class_counts)

balanced_ds = []
for class_idx in range(N_CLASSES):
    subset = df[df['label'] == class_idx].sample(n=n_examples, random_state=RANDOM_SEED)
    balanced_ds.append(subset)

df = pd.concat(balanced_ds)
df = df.sample(frac=1, random_state=RANDOM_SEED)
df = df.drop('site_idx', axis=1)
df['allele_len'] = df.allele.str.len()
df = df.drop('allele', axis=1)
df = df.drop('allele_idx', axis=1)
df = df.drop('chrom', axis=1)
df = df.drop('fingerprint_match', axis=1)
df = df.replace(dict(variant_type=dict(SNP=0, INDEL=1)))

cat_cols = ['variant_type', 'overlaps_with_homo', 'overlaps_with_dup', 'overlaps_with_high_gc', 'overlaps_with_low_gc', 'overlaps_with_low_complex']

#train_ds = df[forest_fields].to_numpy()
training_cols = [col for col in df.columns if col != 'label']
#training_cols = [col for col in df.columns if not (col == 'label' or col.startswith('AS'))]
#training_cols = [col for col in df.columns if not (col == 'label' or col.startswith('overlaps'))]
cat_mask = [col in cat_cols for col in training_cols]
train_ds = df[training_cols].to_numpy()
#print(train_ds.shape)
#print(list(np.mean(train_ds, axis=0)))
#print(list(np.std(train_ds, axis=0)))
#train_ds = StandardScaler().fit_transform(train_ds)
print(sum(list(np.mean(train_ds, axis=0))))
print(sum(list(np.std(train_ds, axis=0))))
labels_ds = df['label'].to_numpy()
test_count = int(round(len(train_ds) * split_ratio))
(train_input, train_labels) = (train_ds[:test_count], labels_ds[:test_count])
(test_input, test_labels) = (train_ds[test_count:], labels_ds[test_count:])
print(f'train set N={np.sum(train_labels)}, {np.sum(train_labels) / len(train_labels) * 100:.02f}% positive')
print(f'test set N={np.sum(test_labels)}, {np.sum(test_labels) / len(test_labels) * 100:.02f}% positive')

cat_kw = {
    'categorical_features': cat_mask
}
kw={
    'random_state': RANDOM_SEED,
}

clfs = [
    RandomForestClassifier(n_jobs=-1, **kw),
    GradientBoostingClassifier(**kw),
    #HistGradientBoostingRegressor(**kw),
    #AdaBoostClassifier(**kw),
    #KNeighborsClassifier(n_neighbors=2),
    #MLPClassifier()
]

for clf in clfs:
    clf_name = clf.__class__.__name__.split('.')[-1]
    clf.fit(train_input, train_labels)
    msg = f"{clf_name} accuracy with test labels {clf.score(test_input, test_labels) * 100:.02f}%"
    print(msg)
    #pfn = f'models/{clf_name}.pickle'
    #with open(pfn, 'wb') as fh:
        #pickle.dump(clf, fh)
