import time
import pickle
import multiprocessing as mp
from textwrap import wrap

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from . iters import *
from .. rep.region import Region
from .. rep.fingerprint import AlleleFingerprint
from .. utils import region_progress_bar
from .. io.monitor import get_remote_monitor

def get_columns_by_prefix(df=None, prefix=None):
    return [col for col in df.columns if col.startswith(prefix)]

def get_feature_columns(df=None, prefix='feature_'):
    return get_columns_by_prefix(df=df, prefix=prefix)

def get_label_columns(df=None, prefix='label_'):
    return get_columns_by_prefix(df=df, prefix=prefix)

class Classifier:
    def __init__(self, classifier=None, label_cols=None, feature_cols=None, is_trained=False):
        self.classifier = classifier
        if hasattr(self.classifier, 'n_jobs'):
            self.classifier.n_jobs = -1
        self._label_cols = label_cols
        self._feature_cols = feature_cols
        self.is_trained = is_trained

    def __getstate__(self):
        return (self.classifier, self._label_cols, self._feature_cols, self.is_trained)

    def __setstate__(self, state):
        (self.classifier, self._label_cols, self._feature_cols, self.is_trained) = state

    def get_label_cols(self):
        return self._label_cols
    def set_label_cols(self, val):
        if self._label_cols is not None:
            msg = 'label column is already configured'
            raise ValueError(msg)
        self._label_cols = val
    label_cols = property(get_label_cols, set_label_cols)

    def get_feature_cols(self):
        return self._feature_cols
    def set_feature_cols(self, val):
        if self._feature_cols is not None:
            msg = 'features column is already configured'
            raise ValueError(msg)
        self._feature_cols = val
    feature_cols = property(get_feature_cols, set_feature_cols)

    @property
    def classifier_class(self):
        return self.classifier.__class__

    @property
    def classifier_name(self):
        return self.classifier_class.__name__

    def infer_columns(self, df=None, label_cols=None, feature_cols=None):
        self.label_cols = label_cols or self.label_cols or get_label_columns(df=df)
        self.feature_cols = feature_cols or self.feature_cols or get_feature_columns(df=df)
        return (self.label_cols, self.feature_cols)

    def fit(self, df=None, label_cols=None, feature_cols=None):
        if self.is_trained:
            msg = 'classifier is already trained'
            raise ValueError(msg)

        self.infer_columns(df=df, label_cols=label_cols, feature_cols=feature_cols)
        rpt = f'Fitting {self.classifier_name} to {len(df)} samples with {len(self.feature_cols)} features:\n'
        rpt += str.join('\n', ['  ' + ln for ln in wrap(str.join(' ', self.feature_cols))]) + '\n'
        print(rpt)
        #
        inp = df[self.feature_cols].to_numpy()
        gt = np.squeeze(df[self.label_cols].to_numpy())
        self.classifier.fit(inp, gt)
        self.is_trained = True

    def score(self, df=None):
        if not self.is_trained:
            msg = 'classifier is not trained'
            raise ValueError(msg)
        inp = df[self.feature_cols].to_numpy()
        gt = df[self.label_cols].to_numpy()
        return self.classifier.score(inp, gt)

    def fit_and_score(self, df=None, test_df=None, test_frac=None, **kw):
        if test_df is None:
            test_size = int(round(len(df) * test_frac))
            (df, test_df) = train_test_split(df, test_size=test_size)
        self.fit(df=df, **kw)
        acc = self.score(df=test_df)
        rpt = f'Final {self.classifier_name} accuracy is {acc * 100:.02f}% using {test_size} test samples'
        print(rpt)
        return acc

    def predict(self, df=None, mode=None, score_col='score', epsilon=1e-5):
        if not self.is_trained:
            msg = 'classifier is not trained'
            raise ValueError(msg)
        if mode == 'log_proba':
            func = self.classifier.predict_log_proba
        elif mode in ('proba', 'logodds'):
            func = self.classifier.predict_proba
        elif mode is None:
            func = self.classifier.predict
        else:
            raise ValueError(mode)
        features = df[self.feature_cols]
        features = numlint(features)
        with np.errstate(divide='ignore'):
            res = func(features.to_numpy())
            if mode == 'logodds':
                res[res == np.inf] = 1
                res[res == -np.inf] = 0
                res[res == 0] = epsilon
                assert res.shape[-1] == 2
                res = np.squeeze(np.log(res[:, 1] / res[:, 0]))
        if len(res.shape) == 2 and res.shape[-1] == 2:
            if res.shape[-1] != 2:
                raise NotImplementedError('n_classes != 2')
            res = np.squeeze(res[:, 1])
        df[score_col] = res
        return df
    
    def save_classifier(self, path=None):
        with open(path, 'wb') as fh:
            pickle.dump(self, fh)

    @classmethod
    def load_classifier(cls, path=None):
        with open(path, 'rb') as fh:
            return pickle.load(fh)

# XXX: task specific
class AnnotateCozy:
    DefaultFields = [
        'AS_BaseQRankSum', 'AS_FS', 'AS_InbreedingCoeff', 'AS_MQ',
        'AS_MQRankSum', 'AS_QD', 'AS_ReadPosRankSum', 'AS_SOR'
    ]

    def __init__(self, field_names=None, default_value=0):
        self.field_names = tuple(field_names or self.DefaultFields)
        self.default_value = 0

    def __call__(self, row):
        if not row.filter:
            #row.feature.delta_length = abs(len(row.meta.ref) - len(row.meta.allele))
            row.feature.delta_length = len(row.meta.ref) - len(row.meta.allele)
            row.feature.is_snp = (len(row.meta.ref) == len(row.meta.allele))
            site = row.cache.site
            for fn in self.field_names:
                if fn.startswith('AS_'):
                    try:
                        row.feature[fn] = site.info[fn][row.meta.allele_idx]
                    except KeyError:
                        row.feature[fn] = self.default_value
                else:
                    try:
                        row.feature[fn] = site.info[fn]
                    except KeyError:
                        row.feature[fn] = self.default_value
        return row

def balance_dataframe(df=None, label_cols=None, random_seed=None):
    if label_cols == None:
        label_cols = get_label_columns(df)
    n_classes = len(label_cols)
    class_counts = []
    for label in label_cols:
        class_counts += df[label].value_counts().to_list()
    min_class = np.argmin(class_counts)
    n_examples = np.min(class_counts)
    balanced_ds = []
    for label_col in label_cols:
        for value in df[label_col].unique():
            subset = df[df[label_col] == value].sample(n=n_examples, random_state=random_seed)
            balanced_ds.append(subset)
    df = pd.concat(balanced_ds)
    df = df.sample(frac=1, random_state=random_seed)
    return df

def classify_vcf(vcf=None, region=None, classifier=None, overlaps=None, annotate=None, batch_size=10_000, assembly=None, as_scheme=None, remote=None):
    itr = iter_sites(vcf=vcf, region=region, assembly=assembly, as_scheme=as_scheme)
    if remote:
        itr = pos_monitor(itr, remote)
        itr = iter_monitor(itr, remote, 'sites')
    #itr = filter_by_site(itr=itr)
    if overlaps is not None:
        itr = overlaps_with_site(itr, overlaps=overlaps)
    itr = iter_alleles(itr=itr)
    if remote:
        itr = iter_monitor(itr, remote, 'alleles')
    #itr = filter_by_allele(itr=itr)
    if annotate:
        itr = custom_itr(itr, annotate)
    df = to_dataframe(itr)
    return df

def numlint(thing, posinf=0, neginf=0, nan=0, **extra):
    rules = {np.inf: posinf, -np.inf: neginf, np.nan: nan}
    rules.update(extra)
    if isinstance(thing, pd.DataFrame):
        return thing.replace(rules)
    if isinstance(thing, np.ndarray):
        for (old, new) in rules.items():
            thing[thing == old] = new
        return thing
    raise TypeError(type(thing))

class ClassifyTask:
    def __init__(self, query_vcf=None, classifier=None, overlaps=None, annotate=None, assembly=None, as_scheme=None, progress_bar=True):
        self.query_vcf = query_vcf
        self.classifier = classifier
        self.overlaps = overlaps
        self.annotate = annotate
        self.assembly = assembly
        self.as_scheme = as_scheme
        self.remote = get_remote_monitor(domain='CLS')
    
    def __call__(self, region=None):
        df = classify_vcf(
            vcf=self.query_vcf, 
            region=region, 
            classifier=self.classifier,
            overlaps=self.overlaps, 
            annotate=self.annotate, 
            assembly=self.assembly, 
            as_scheme=self.as_scheme,
            remote=self.remote
        )
        return (region, df)

    def classify_region(self, region=None, mode='logodds', chunk_size=1_000_000):
        if not isinstance(region, Region):
            region = Region(region)
        if len(region) == 0:
            vcf_contig = self.query_vcf.header.contigs[region.contig]
            region = region.clone(start=1, stop=vcf_contig.length)
        regions = region.split(chunk_size)
        for region in regions:
            (cur_region, df)  = self(region)
            if df is None:
                continue
            if len(df) > 0:
                df = self.classifier.predict(df, mode=mode)
            yield (cur_region, df)

    def call_vcf_sites(self, output_vcf=None, columns=None, region=None, **kw):
        df_itr = self.classify_region(region=region, **kw)
        last_region = None
        for (cur_region, df) in df_itr: 
            if not len(df):
                continue
            assert last_region is None or last_region.stop < cur_region.start
            last_region = cur_region
            itr = iter_sites(
                vcf=self.query_vcf, region=cur_region, 
                assembly=self.assembly, as_scheme=self.as_scheme
            )
            itr = annotate_alleles_from_dataframe(itr=itr, df=df, columns=columns)
            for row in itr:
                # XXX not here, hardwired
                #row.cache.site.info['BLOD'] = min(row.cache.site.info['AS_BLOD'])
                row.cache.site.info['BLOD'] = min(row.cache.site.info['AS_BLOD'])
                output_vcf.write(row.cache.site)
