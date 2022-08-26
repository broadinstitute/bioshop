import time
import pickle
import multiprocessing as mp
from textwrap import wrap

from sklearn.model_selection import train_test_split

from . iters import *
from .. rep.region import Region
from .. rep.fingerprint import AlleleFingerprint
from .. utils import region_progress_bar

class Classifier:
    def __init__(self, classifier=None, label=None, features=None, is_trained=False):
        self.classifier = classifier
        self._label = label
        self._features = features
        self.is_trained = is_trained

    def get_label(self):
        return self._label
    def set_label(self, val):
        if self._label is not None:
            msg = 'label column is already configured'
            raise ValueError(msg)
        self._label = val
    label = property(get_label, set_label)

    def get_features(self):
        return self._features
    def set_features(self, val):
        if self._features is not None:
            msg = 'features column is already configured'
            raise ValueError(msg)
        self._features = val
    features = property(get_features, set_features)

    @property
    def classifier_class(self):
        return self.classifier.__class__

    @property
    def classifier_name(self):
        return self.classifier_class.__name__

    def fit(self, df=None, label=None, features=None):
        if self.is_trained:
            msg = 'classifier is already trained'
            raise ValueError(msg)
        label = label or self.label or 'label'
        assert label is not None
        features = features or self.features
        if features is None:
            features = sorted(set(df.columns) - set([label]))
        self.label = label
        self.features = features
        rpt = f'Fitting {self.classifier_name} to {len(df)} samples with {len(self.features)} features:\n'
        rpt += str.join('\n', ['  ' + ln for ln in wrap(str.join(' ', self.features))]) + '\n'
        print(rpt)
        #
        inp = df[self.features].to_numpy()
        gt = df[self.label].to_numpy()
        self.classifier.fit(inp, gt)
        self.is_trained = True

    def score(self, df=None):
        if not self.is_trained:
            msg = 'classifier is not trained'
            raise ValueError(msg)
        inp = df[self.features].to_numpy()
        gt = df[self.label].to_numpy()
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

    def predict(self, df=None, mode=None):
        if not self.is_trained:
            msg = 'classifier is not trained'
            raise ValueError(msg)
        if mode == 'log':
            func = self.classifier.predict_log_proba
        elif mode == 'proba':
            func = self.classifier.predict_proba
        elif mode is None:
            func = self.classifier.predict
        else:
            raise ValueError(mode)
        features = df[self.features].numpy()
        return func(features)
    
    def save_classifier(self, path=None):
        with open(path, 'wb') as fh:
            pickle.dump(self, fh)

    @classmethod
    def load_classifier(cls, path=None):
        with open(path, 'rb') as fh:
            return pickle.load(fh)

def classify_vcf(vcf=None, region=None, overlaps=None):
    itr = iter_sites(vcf=vcf, with_index=True, region=region)
    if overlaps is not None:
        itr = overlaps_with_site(itr, overlaps=overlaps)
    itr = skip_site(itr=itr)
    itr = iter_alleles(itr=itr, with_index=True)
    itr = skip_allele(itr=itr)
    return fingerprint_allele(itr=itr)

class ClassifyTask:
    def __init__(self, vcf=None, classifier=None):
        self.vcf = vcf
        self.classifier = classifier

    def __call__(self, region=None):
        return classify_vcf(vcf=self.vcf, region=region, classifier=self.classifier)

    def classify(self, region=None, chunk_size=50_000):
        if not isinstance(region, Region):
            region = Region(region)
        regions = region.split(chunk_size)
        pool = mp.Pool()
        yield from pool.imap(self, regions)

