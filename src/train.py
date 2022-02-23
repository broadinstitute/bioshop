import os
import ctypes
import torch
import pandas as pd
import numpy as np
from pyfaidx import Fasta
from torch.utils.data import IterableDataset
from transformers import TrainingArguments, Trainer
from sklearn.preprocessing import QuantileTransformer

from . models import VariantTokenizer, ModelInputStruct, VariantToVector, TabularVariantFilterModel
from . work import worker
from . work import cbuf

class TrainDumb(IterableDataset):
    def __init__(self, klen=3, size=1000):
        self.size = size
        self.klen = klen

    def __len__(self):
        return self.size

    def __iter__(self):
        for idx in range(self.size):
            size = 512
            n_cats = n_nums = 5
            inp = dict(
                input_ids=([0] * size),
                attention_mask=([0] * size),
                token_type_ids=([0] * size),
                cat_feats=[float(x) for x in range(n_cats)],
                numerical_feats=[float(x) for x in range(n_nums)],
                labels=[1],
            )
            inp = {key: torch.tensor(val) for (key, val) in inp.items()}
            yield inp

class TrainingDataIter(object):
    #DefaultColNums = ["INFO_DP", "INFO_QD", "INFO_MQRankSum", "INFO_ReadPosRankSum", "INFO_FS", "INFO_SOR"]
    DefaultColNums = ("DP", "QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR")
    DefaultColCats = ("phased", )
    DefaultLabels = ("SNP_TP", "SNP_FP", "INDEL_TP", "INDEL_FP")

    def __init__(self, dataset=None, ref_path=None, tokenizer_config=None, num_col_list=None, cat_col_list=None, labels=None, **kw):
        super().__init__(**kw)
        self.dataset = dataset.replace(np.nan, 0)
        self.ref_path = ref_path
        self.tokenizer_config = tokenizer_config
        # categorical
        self.cat_col_list = list(cat_col_list or self.DefaultColCats)
        # numerical
        self.num_col_list = list(num_col_list or self.DefaultColNums)
        # labels
        self.labels = labels or self.DefaultLabels
        self.transform_numerical_values()
        self.relabel_dataset()

    def transform_numerical_values(self):
        # XXX: all?
        self.numerical_transformer = QuantileTransformer()
        col_list = [f'INFO_{key}' for key in self.num_col_list]
        num_vals = self.dataset[col_list]
        num_vals = self.numerical_transformer.fit_transform(num_vals)
        self.dataset[col_list] = num_vals

    def relabel_dataset(self):
        #
        label_map = {lb: self.labels.index(lb) for lb in self.labels}
        label_f = lambda label: label_map[label]
        self.dataset['label_value'] = self.dataset['label'].apply(label_f)
        #
        self.dataset['phased_value'] = self.dataset['phased'].apply(int)

    @classmethod
    def load(cls, dataset_path=None, **kw):
        print(f"Loadining dataset from {dataset_path}")
        dataset = pd.read_csv(dataset_path, sep='\t')
        return cls(dataset=dataset, **kw)

    def __len__(self):
        return len(self.dataset)

    def __iter__(self):
        ref = Fasta(self.ref_path)
        tokenizer = VariantTokenizer(**self.tokenizer_config)
        vectorizer = VariantToVector(ref=ref, tokenizer=tokenizer)
        num_col_list = [f'INFO_{key}' for key in self.num_col_list]
        cat_col_list = ['phased_value']

        for (idx, row) in self.dataset.iterrows():
            # XXX: bug, kill eval()
            row.gt_bases = eval(row.gt_bases)
            vecs = vectorizer.process_training_site(row)
            if vecs is None:
                continue
            feats = {}
            feats['numerical_feats'] = row[num_col_list].to_numpy(dtype=np.float32)
            feats['cat_feats'] = row[cat_col_list].to_numpy(dtype=np.float32)
            #feats = {key: torch.tensor(val).float() for (key, val) in feats.items()}
            #vecs = {key: torch.tensor(val) for (key, val) in vecs.items()}
            feats.update(vecs)
            #feats['labels'] = torch.tensor(row['label_value'])
            feats['labels'] = np.array([row['label_value']], dtype=np.int32)
            yield feats

class ModelTrainingStruct(ctypes.Structure):
    N_NUM_FEATS = 6
    N_CAT_FEATS = 1

    _fields_ = [
        ('input_ids', ctypes.c_int * 512), 
        ('attention_mask', ctypes.c_int * 512),
        ('token_type_ids', ctypes.c_int * 512),
        ('numerical_feats', ctypes.c_float * N_NUM_FEATS),
        ('cat_feats', ctypes.c_float * N_CAT_FEATS),
        ('labels', ctypes.c_int * 1),
    ]

    def __init__(self, **kw):
        c_cast = lambda it: np.ctypeslib.as_ctypes(it) if isinstance(it, np.ndarray) else it
        kw = {key: c_cast(val) for (key, val) in kw.items()}
        super().__init__(**kw)

    def as_numpy(self):
        arrays = [
            np.ctypeslib.as_array(self.input_ids),
            np.ctypeslib.as_array(self.attention_mask),
            np.ctypeslib.as_array(self.token_type_ids),
            np.ctypeslib.as_array(self.numerical_feats),
            np.ctypeslib.as_array(self.cat_feats),
            np.ctypeslib.as_array(self.labels),
        ]
        return np.array(arrays)

    def as_dict(self):
        arrays = self.as_numpy()
        keys = [field[0] for field in self._fields_]
        return dict(zip(keys, arrays))

class TrainingWorker(worker.Worker):
    def __init__(self, worker_config=None, maxsize=1024, **kw):
        super().__init__(**kw)
        self.worker_config = worker_config
        self.cbuf = cbuf.CircularBuffer(ctype=ModelTrainingStruct, size=maxsize)

    def _run(self):
        itr = TrainingDataIter.load(**self.worker_config)
        for item in itr:
            item = ModelTrainingStruct(**item)
            self.cbuf.push(item)

    def pop(self, timeout=None):
        return self.cbuf.pop(timeout=timeout)

class TrainingDataset(IterableDataset):
    def __init__(self, dataset_path=None, ref_path=None, tokenizer_config=None, num_col_list=None, cat_col_list=None, labels=None, **kw):
        self.worker_config = dict(
            dataset_path=dataset_path,
            ref_path=ref_path,
            tokenizer_config=tokenizer_config,
            num_col_list=num_col_list,
            cat_col_list=cat_col_list,
            labels=labels,
        )
        self.worker_config.update(kw)
        self.worker = None

    def __len__(self):
        # quick and dirty `wc -l`
        ds_path = self.worker_config['dataset_path']
        count = 0
        with open(ds_path) as fh:
            fh.readline()
            for line in fh:
                count += 1
        return count

    def __iter__(self):
        worker = TrainingWorker(worker_config=self.worker_config)
        worker.start()
        worker.wait_until_running()

        while worker.running:
            try:
                item = worker.pop(timeout=.1)
                item = item.as_dict()
                yield item
            except TimeoutError:
                print("starving")

def train(ref_path=None, klen=None, checkpoint_dir=None, train_fn=None, eval_fn=None):
    tokenizer_config = {'klen':3}
    checkpoint_dir = f"checkpoints/null"

    vfm = TabularVariantFilterModel(klen=klen)
    model = vfm.model

    train_fn = "mostly-training-data/train.tsv"
    eval_fn = "mostly-training-data/eval.tsv"
    short_fn = "mostly-training-data/short.tsv"
    short_fn = os.path.abspath(short_fn)
    eval_fn = train_fn = short_fn
    ref_path = os.path.abspath(ref_path)


    train_ds = TrainingDataset(dataset_path=train_fn, ref_path=ref_path, tokenizer_config=tokenizer_config)
    eval_ds = TrainingDataset(dataset_path=eval_fn, ref_path=ref_path, tokenizer_config=tokenizer_config)

    training_args = TrainingArguments(
        output_dir=checkpoint_dir,
        num_train_epochs=1,
        #max_steps=100,
        per_device_train_batch_size=8,
        per_device_eval_batch_size=2,
        logging_first_step=True,
        logging_strategy="steps",
        logging_steps=50,
    )

    ds_train = TrainDumb()
    ds_eval = TrainDumb()
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_ds,
        eval_dataset=eval_ds,
    )

    trainer.train()

if __name__ == "__main__":
    train()
