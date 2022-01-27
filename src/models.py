import numpy as np
import torch
from transformers import BertForSequenceClassification, BertTokenizer

def get_base_model_path(klen=3):
    assert klen in (3, 4, 5, 6)
    model_path = f"armheb/DNA_bert_{klen}"
    return model_path

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
        self.softmax_op = torch.nn.Softmax(dim=-1)

    def predict(self, inp):
        inp = inp.to(device=self.device)
        outp = self.model(input_ids=inp)
        logits = outp["logits"].detach()
        softmax = self.softmax_op.forward(logits)
        argmax = torch.argmax(softmax, dim=1)[:, None]

        ret = dict(
            logits=logits.cpu().numpy(),
            softmax=softmax.cpu().numpy(),
            argmax=argmax.cpu().numpy()
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

    def tokenize(self, gt=None):
        (a1, a2) = [self.kmerize(al) for al in gt]
        while (len(a1) + len(a2)) > (self.max_length - 3):
            a1 = a1[1:-1]
            a2 = a2[1:-1]
        inp = ["[CLS]"] + a1 + ["[SEP]"] + a2 + ["[SEP]"]
        toks = list(map(self.vocab.__getitem__, inp))
        toks += [0] * (self.max_length - len(toks))
        return toks

    def tokenize_for_training(self, gt=None):
        toks = self.tokenize(gt=gt)
        mask = list(map(bool, toks))
        types = [0] * self.max_length
        ret = {
            "input_ids": toks,
            "attention_mask": list(map(lambda val: int(bool(val)), toks)),
            "token_type_ids": ([0] * self.max_length)
        }
        return ret

if __name__ == "__main__":
    vfm = VariantFilterModel(klen=6)
    inp = [[0, 1, 2, 3], [3, 4, 3, 4]]
    inp = torch.tensor(inp)
    outp = vfm.predict(inp)
    print(outp)