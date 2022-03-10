import numpy as np
import torch
from .. models import VariantFilterModel, TabularVariantFilterModel
from .. base import Genotype
from . worker import Worker

class ModelWorker(Worker):
    def __init__(self, model_path=None, batch_size=None, klen=None, **kw):
        super().__init__(**kw)
        self.model_path = model_path
        self.batch_size = batch_size
        self.klen = klen
        self.in_q = self.manager.to_model
        self.out_q = self.manager.to_gather

    def get_batch_data(self, timeout=1):
        batch = []
        while (len(batch) < self.batch_size):
            try:
                item = self.in_q.pop(timeout=timeout)
            except TimeoutError:
                #print(f"{self.__class__.__name__} get_batch: buf={self.in_q.qsize.value}")
                break
            batch.append(item)
        return batch

    def build_batch(self, batch_rows=None):
        batch_fields = (
            'input_ids', 'attention_mask', 
            'token_type_ids', 'numerical_feats', 
            'cat_feats', 
        )
        (batch_keys, batch_data) = zip(*[(row.get_key(), row.as_dict()) for row in batch_rows])
        to_torch = lambda field: torch.tensor([row[field] for row in batch_data])
        batch = {field: to_torch(field) for field in batch_fields}
        return (batch_keys, batch)

    def get_batch(self, timeout=1):
        batch_rows = self.get_batch_data(timeout=timeout)
        if not batch_rows:
            return (None, None)
        (batch_keys, batch) = self.build_batch(batch_rows=batch_rows)
        return (batch_keys, batch)

    def _run(self):
        #model = VariantFilterModel(model_path=self.model_path, klen=self.klen)
        model = TabularVariantFilterModel(model_path=self.model_path, klen=self.klen)

        while self.running:
            (batch_keys, batch) = self.get_batch()
            if batch:
                outp = model.predict(batch)
                result = []
                assert len(outp['log_odds']) == len(batch_keys)
                for (log_odds, keys) in zip(outp['log_odds'], batch_keys):
                    ns = keys.copy()
                    ns['log_odds'] = log_odds.tolist()
                    ns['status'] = 'called'
                    gt = Genotype(**ns)
                    result.append(gt)
                self.out_q.put(result)
            #
            if (
                (self.manager.flush_model.is_set()) and \
                (self.in_q.qsize.value == 0)
            ):
                break
