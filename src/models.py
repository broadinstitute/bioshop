import ctypes
import torch
import numpy as np
from transformers import BertForSequenceClassification, BertTokenizer, BertConfig, BertModel
from multimodal_transformers.model import BertWithTabular, TabularConfig

def get_base_model_path(klen=3):
    assert klen in (3, 4, 5, 6)
    model_path = f"armheb/DNA_bert_{klen}"
    return model_path

class MM_ModelInputStruct(ctypes.Structure):
    N_NUM_FEATS = 6
    N_CAT_FEATS = 1

    _fields_ = [
        ('input_ids', ctypes.c_int * 512), 
        ('attention_mask', ctypes.c_int * 512),
        ('token_type_ids', ctypes.c_int * 512),
        ('label', ctypes.c_int),
        ('numerical_feats', ctypes.c_float * N_NUM_FEATS),
        ('cat_feats', ctypes.c_int * N_CAT_FEATS),
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
            np.ctypeslib.as_array(self.label),
            np.ctypeslib.as_array(self.numerical_feats),
        ]
        return np.array(arrays)

    def get_key(self):
        return dict(site_id=self.site_id, genotype_id=self.genotype_id)

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
            cat_feat_dim=1,
            numerical_feat_dim=6,
            numerical_bn=True,
            num_labels=4,
        )
        self.bert_config.tabular_config = self.tabular_config
        #self.model = BertWithTabular.from_pretrained(self.model_path, config=self.bert_config).eval().to(device=self.device)
        self.model = BertWithTabular.from_pretrained(self.model_path, config=self.bert_config)
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
    Default_MinTemplateLength = 32
    Token_CLS = "[CLS]"
    Token_SEP = "[SEP]"

    def __init__(self, klen=None, model_path=None, max_length=512, min_template_length=Default_MinTemplateLength):
        # XXX: maxlength from config?
        self.max_length = max_length
        self.klen = klen
        model_path = model_path or get_base_model_path(klen=klen)
        self._tokenizer = BertTokenizer.from_pretrained(model_path, do_lower_case=False)
        self.min_template_length = min_template_length
        self.vocab = self._tokenizer.vocab

    def kmerize(self, seq):
        assert type(seq) is str
        k_or_mask = lambda kmer: kmer if 'M' not in kmer else '[MASK]'
        kmers = [k_or_mask(seq[i:i + self.klen]) for i in range(len(seq) - self.klen + 1)]
        return kmers

    def tokenize(self, gt=None):
        ret = self.tokenize_only(gt=gt)
        toks = ret["input_ids"]
        ret["attention_mask"] = np.array(toks != 0, dtype=ctypes.c_int)
        ret["token_type_ids"] = np.zeros_like(toks, dtype=ctypes.c_int)
        return ret

    def tokenize_only(self, up=None, down=None, varlist=None):
        klen_bridge = lambda vs: up[-(self.klen - 1):] + vs + down[:self.klen - 1]
        middle_out = []
        for varstr in varlist:
            varstr = klen_bridge(varstr)
            k_varstr = self.kmerize(varstr)
            middle_out += [self.Token_SEP] + k_varstr
        middle_out.append(self.Token_SEP)
        if (self.max_length - len(middle_out)) < self.min_template_length:
            return None

        up_k = self.kmerize(up)
        down_k = self.kmerize(down)
        prompt_k = up_k + middle_out + down_k

        while len(prompt_k) > (self.max_length - 1):
            prompt_k = prompt_k[1:-1]
        inp = [self.Token_CLS] + prompt_k
        toks = list(map(self.vocab.__getitem__, inp))
        toks += [0] * (self.max_length - len(toks))
        ret = {
            "input_ids": np.array(toks, dtype=ctypes.c_int),
        }
        return ret

    def tokenize(self, up=None, down=None, varlist=None):
        ret = self.tokenize_only(up=up, down=down, varlist=varlist)
        if ret is None:
            return None
        toks = ret["input_ids"]
        ret["attention_mask"] = np.array(toks != 0, dtype=ctypes.c_int)
        ret["token_type_ids"] = np.zeros_like(toks, dtype=ctypes.c_int)
        return ret

class VariantToVector(object):
    def __init__(self, ref=None, tokenizer=None):
        self.ref = ref
        self.tokenizer = tokenizer
        self.window = self.tokenizer.max_length // 2
    
    def get_ref_up_down(self, site_info=None):
        pos = site_info.pos
        chrom = site_info.chrom
        ref = site_info.ref

        start = pos - 1
        end = start + len(ref)
        up = self.ref[chrom][start - self.window:start]
        down = self.ref[chrom][end:end + self.window]
        return (str(up), str(down))
    
    def process_training_site(self, site_info=None):
        (up, down) = self.get_ref_up_down(site_info=site_info)
        varlist = [site_info.ref] + site_info.gt_bases
        try:
            gt_tok = self.tokenizer.tokenize(up=up, down=down, varlist=varlist)
        except KeyError:
            return None
        return gt_tok

    def process_site(self, site_info=None):
        (up, down) = self.get_ref_up_down(site_info=site_info)
        gt_toks = {}
        ref_bases = (site_info.ref, site_info.ref)
        for (genotype_id, gt) in site_info.genotypes.items():
            if (gt.bases == ('.', '.')) or (gt.bases == ref_bases):
                gt_ref = None
            else:
                all_chars = str.join('', gt.bases) + site_info.ref + up + down
                if set(all_chars.upper()) - set('AGTC'):
                    gt_ref = None
                else:
                    varlist = [site_info.ref] + list(gt.bases)
                    gt_ref = self.tokenizer.tokenize(up=up, down=down, varlist=varlist)
            gt_toks[genotype_id] = gt_ref
        return gt_toks

def test():
    return
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

if __name__ == "__main__":
    test()
