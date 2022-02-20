import numpy as np
import torch
from .. models import VariantFilterModel
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

    def get_batch(self, timeout=1):
        batch = []
        while (len(batch) < self.batch_size):
            try:
                item = self.in_q.pop(timeout=timeout)
            except TimeoutError:
                #print(f"{self.__class__.__name__} get_batch: buf={self.in_q.qsize.value}")
                break
            batch.append(item)
        if not batch:
            return (None, None)
        keys = [ds.get_key() for ds in batch]
        batch = np.array([ds.as_numpy() for ds in batch])
        batch = batch.transpose([1, 0, 2])
        batch = torch.tensor(batch)
        headers = ('input_ids', 'attention_mask', 'token_type_ids')
        batch = dict(zip(headers, batch))
        return (keys, batch)

    def _run(self):
        model = VariantFilterModel(model_path=self.model_path, klen=self.klen)

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
