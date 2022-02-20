import ctypes
import torch
import numpy as np
from transformers import BertForSequenceClassification, BertTokenizer, BertConfig, BertModel
from multimodal_transformers.model import BertWithTabular, TabularConfig

def get_base_model_path(klen=3):
    assert klen in (3, 4, 5, 6)
    model_path = f"armheb/DNA_bert_{klen}"
    return model_path

class ModelInputStruct(ctypes.Structure):
    _fields_ = [
        ('input_ids', ctypes.c_int * 512), 
        ('attention_mask', ctypes.c_int * 512),
        ('token_type_ids', ctypes.c_int * 512),
        ('site_id', ctypes.c_int),
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
        return dict(site_id=self.site_id, genotype_id=self.genotype_id)

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
        self.model = BertForSequenceClassification.from_pretrained(self.model_path, num_labels=self.num_labels).eval().to(device=self.device)
        self.softmax_op = torch.nn.Softmax(dim=-1)

    def predict(self, inp):
        if type(inp) != dict:
            inp = {"input_ids": inp}
        inp = {key: val.to(device=self.device) for (key, val) in inp.items()}
        with torch.no_grad():
            # return loss, logits, classifier_layer_outputs
            logits = self.model(**inp)["logits"]
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


class TabularVariantFilterModel(VariantFilterModel):
    def __init__(self, klen=None, labels=None, model_path=None, device=None, bert_config=None):
        assert klen in (3, 4, 5, 6)
        self.klen = klen
        if model_path is None:
            model_path = get_base_model_path(klen=self.klen)
        self.labels = labels or self.DefaultLabels
        self.num_labels = len(self.labels)
        self.model_path = model_path
        self.device = device or self.DefaultDevice
        if bert_config is None:
            bert_config = BertConfig(
                vocab_size=4 ** self.klen + 5
            )
        self.bert_config = bert_config
        self.tabular_config = TabularConfig(
            combine_feat_method="mlp_on_concatenated_cat_and_numerical_feats_then_concat",
            cat_feat_dim=5,
            numerical_feat_dim=5,
            numerical_bn=True,
            num_labels=4,
        )
        self.bert_config.tabular_config = self.tabular_config
        self.model = BertWithTabular.from_pretrained(self.model_path, config=self.bert_config).eval().to(device=self.device)
        self.softmax_op = torch.nn.Softmax(dim=-1)

    def predict(self, inp):
        if type(inp) != dict:
            inp = {"input_ids": inp}
        inp = {key: val.to(device=self.device) for (key, val) in inp.items()}
        with torch.no_grad():
            # return loss, logits, classifier_layer_outputs
            logits = self.model(**inp)[1]
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

class VariantToVector(object):
    def __init__(self, ref=None, tokenizer=None, window=96):
        self.ref = ref
        self.tokenizer = tokenizer
        self.window = window
    
    def get_reference_context(self, site_info=None):
        pos = site_info.pos
        chrom = site_info.chrom
        ref = site_info.ref

        start = pos - 1
        end = start + len(ref)
        pre = self.ref[chrom][start - self.window:start]
        post = self.ref[chrom][end:end + self.window]
        return (pre, post)
    
    def process_site(self, site_info=None):
        (pre, post) = self.get_reference_context(site_info=site_info)
        gt_toks = {}
        ref_bases = (site_info.ref, site_info.ref)
        for (genotype_id, gt) in site_info.genotypes.items():
            if (gt.bases == ('.', '.')) or (gt.bases == ref_bases):
                gt_ref = None
            else:
                gt_ref = [f"{pre}{var if var != '.' else site_info.ref}{post}" for var in gt.bases]
                if set(gt_ref[0] + gt_ref[1]) - set('AGTC'):
                    gt_ref = None
                else:
                    gt_ref = self.tokenizer.tokenize(gt_ref)
            gt_toks[genotype_id] = gt_ref
        return gt_toks

if __name__ == "__main__":
    vfm = TabularVariantFilterModel(klen=6)
    #vfm = VariantFilterModel(klen=6)
    size = 512
    n_cats = n_nums = 5
    inp = dict(
        input_ids=[list(range(size))],
        attention_mask=[[0] * size],
        token_type_ids=[[0] * size],
        cat_feats=[[float(x) for x in range(n_cats)]],
        numerical_feats=[[float(x) for x in range(n_nums)]],
    )
    inp = {key: torch.tensor(val) for (key, val) in inp.items()}
    outp = vfm.predict(inp)
    print(outp)
