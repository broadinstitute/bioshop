import ctypes
import torch
import numpy as np
from transformers import BertForSequenceClassification, BertTokenizer

def get_base_model_path(klen=3):
    assert klen in (3, 4, 5, 6)
    model_path = f"armheb/DNA_bert_{klen}"
    return model_path

class ModelInputStruct(ctypes.Structure):
    _fields_ = [
        ('input_ids', ctypes.c_int * 512), 
        ('attention_mask', ctypes.c_int * 512),
        ('token_type_ids', ctypes.c_int * 512),
        ('call_id', ctypes.c_int),
        ('genotype_id', ctypes.c_int),
    ]

    def __init__(self, **kw):
        c_cast = lambda it: np.ctypeslib.as_ctypes(it) if isinstance(it, np.ndarray) else it
        kw = {key: c_cast(val) for (key, val) in kw.items()}
        super().__init__(**kw)

    def as_numpy(self):
        arrays = [
            np.ctypeslib.as_array(self.input_ids),
            np.ctypeslib.as_array(self.attention_mask),
            np.ctypeslib.as_array(self.token_type_ids)
        ]
        return np.array(arrays)

    def get_key(self):
        return dict(call_id=self.call_id, genotype_id=self.genotype_id)

class VariantFilterModel(object):
    DefaultLabels = ("SNP_TP", "SNP_FP", "INDEL_TP", "INDEL_FP")
    DefaultDevice = "cuda:0"

    def __init__(self, klen=None, labels=None, model_path=None, device=None):
        assert klen in (3, 4, 5, 6)
        self.klen = klen
        if model_path is None:
            model_path = get_base_model_path(klen=self.klen)
        self.labels = labels or self.DefaultLabels
        self.num_labels = len(self.labels)
        self.model_path = model_path
        self.device = device or self.DefaultDevice
        self.model = BertForSequenceClassification.from_pretrained(self.model_path, num_labels=self.num_labels).to(device=self.device)
        self.model.to(self.device)
        self.model.eval()
        self.softmax_op = torch.nn.Softmax(dim=-1)

    def predict(self, inp):
        if type(inp) != dict:
            inp = {"input_ids": inp}
        inp = {key: val.to(device=self.device) for (key, val) in inp.items()}
        with torch.no_grad():
            outp = self.model(**inp)
        logits = outp["logits"]
        # XXX: hard wired for 2x2 classes
        logits = logits.reshape(logits.shape[0], 2, 2)
        softmax = self.softmax_op.forward(logits)
        argmax = torch.argmax(softmax, dim=1)[:, None]
        log_odds = torch.log(softmax[:, :, 0] / softmax[:, :, 1])
        ret = dict(
            logits=logits.cpu().numpy(),
            softmax=softmax.cpu().numpy(),
            argmax=argmax.cpu().numpy(),
            log_odds=log_odds.cpu().numpy(),
        )
        return ret

class VariantTokenizer(object):
    def __init__(self, klen=None, model_path=None, max_length=512):
        # XXX: maxlength from config?
        self.max_length = max_length
        self.klen = klen
        model_path = model_path or get_base_model_path(klen=klen)
        self._tokenizer = BertTokenizer.from_pretrained(model_path, do_lower_case=False)
        self.vocab = self._tokenizer.vocab
    
    def kmerize(self, seq):
        assert type(seq) is str
        k_or_mask = lambda kmer: kmer if 'M' not in kmer else '[MASK]'
        kmers = [k_or_mask(seq[i:i + self.klen]) for i in range(len(seq) - self.klen + 1)]
        return kmers

    def tokenize_only(self, gt=None):
        (a1, a2) = [self.kmerize(al) for al in gt]
        while (len(a1) + len(a2)) > (self.max_length - 3):
            a1 = a1[1:-1]
            a2 = a2[1:-1]
        inp = ["[CLS]"] + a1 + ["[SEP]"] + a2 + ["[SEP]"]
        toks = list(map(self.vocab.__getitem__, inp))
        toks += [0] * (self.max_length - len(toks))
        ret = {
            "input_ids": np.array(toks, dtype=ctypes.c_int),
        }
        return ret

    def tokenize(self, gt=None):
        ret = self.tokenize_only(gt=gt)
        toks = ret["input_ids"]
        ret["attention_mask"] = np.array(toks != 0, dtype=ctypes.c_int)
        ret["token_type_ids"] = np.zeros_like(toks, dtype=ctypes.c_int)
        return ret

if __name__ == "__main__":
    vfm = VariantFilterModel(klen=6)
    inp = [[0, 1, 2, 3], [3, 4, 3, 4]]
    inp = torch.tensor(inp)
    outp = vfm.predict(inp)
    print(outp)
