import torch
import pandas as pd
import numpy as np
from pyfaidx import Fasta
from torch.utils.data import IterableDataset
from transformers import TrainingArguments, Trainer
from sklearn.preprocessing import QuantileTransformer

from mostly.models import VariantTokenizer, ModelInputStruct, VariantToVector, TabularVariantFilterModel

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

class TrainingDataIter(IterableDataset):
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
    def load(cls, tsv_fn=None, **kw):
        print(f"Loadining dataset from {tsv_fn}")
        dataset = pd.read_csv(tsv_fn, sep='\t')
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
            feats['numerical_feats'] = row[num_col_list].to_numpy(dtype=np.float)
            feats['cat_feats'] = row[cat_col_list].to_numpy(dtype=np.float)
            feats = {key: torch.tensor(val).float() for (key, val) in feats.items()}
            vecs = {key: torch.tensor(val) for (key, val) in vecs.items()}
            feats.update(vecs)
            feats['labels'] = torch.tensor(row['label_value'])
            yield feats

def train(ref_path=None, klen=None, checkpoint_dir=None, train_fn=None, eval_fn=None):
    tokenizer_config = {'klen':3}
    checkpoint_dir = f"checkpoints/null"

    vfm = TabularVariantFilterModel(klen=klen)
    model = vfm.model

    train_fn = "mostly-training-data/train.tsv"
    eval_fn = "mostly-training-data/eval.tsv"

    train_ds = TrainingDataIter.load(train_fn, ref_path=ref_path, tokenizer_config=tokenizer_config)
    eval_ds = TrainingDataIter.load(eval_fn, ref_path=ref_path, tokenizer_config=tokenizer_config)

    training_args = TrainingArguments(
        output_dir=checkpoint_dir,
        num_train_epochs=1,
        #max_steps=100,
        per_device_train_batch_size=8,
        per_device_eval_batch_size=2,
        logging_first_step=True,
        logging_strategy="steps",
        logging_steps=1,
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
