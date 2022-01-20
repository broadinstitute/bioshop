import re
import glob

class CacheFile(object):
    def __init__(self, filename=None, poslist=None):
        self.filename = filename
        self.poslist = list() if poslist is None else poslist
        if self.poslist:
            assert os.path.exists(self.filename)
            self.fh = open(self.filename, 'a')
        else:
            self.fh = open(self.filename, 'w')
    
    def save(self):
        fn = self.filename + "-config.json"
        cfg = {
            "filename": self.filename,
            "poslist": self.poslist,
        }
        with open(fn, 'w') as fh:
            json.dump(cfg, fh)
    
    @classmethod
    def load(cls, filename=None):
        fn = filename + "-config.json"
        with open(fn) as fh:
            cfg = json.load(fh)
        return cls(**cfg)

    def add(self, call=None):
        self.poslist.append(self.fh.tell())
        outp = json.dumps(call) + '\n'
        self.fh.write(outp)
    
    def __len__(self):
        return len(self.poslist)
    
    def select(self, index_list=None):
        if not self.fh.closed:
            self.fh.close()
        endpos = os.path.getsize(self.filename)
        with open(self.filename) as fh:
            for idx in index_list:
                pos = self.poslist[idx]
                fh.seek(pos)
                if idx == (len(self.poslist) - 1):
                    rlen = endpos - pos
                else:
                    rlen = self.poslist[idx+1] - pos
                js = fh.read(rlen)
                yield json.loads(js)

class CallsetCache(object):
    def __init__(self, path=None, labels=None, files=None):
        self.path = os.path.abspath(path or "callset_cache")
        os.makedirs(self.path, exist_ok=True)
        self.labels = labels
        self.files = files
        if self.files is None:
            self.files = {lb: CacheFile(f"{self.path}/{lb}.json")  for lb in self.labels}
    
    def save(self):
        fn = self.path + "/cache-config.json"
        cfg = {
            "path": self.path,
            "labels": self.labels,
        }
        with open(fn, 'w') as fh:
            json.dump(cfg, fh)
        for fo in self.files:
            self.files[fo].save()

    @classmethod
    def load(cls, path=None):
        fn = path + "/cache-config.json"
        with open(fn) as fh:
            cfg = json.load(fh)
        files = {}
        for label in cfg["labels"]:
            filename = f"{cfg['path']}/{label}.json"
            files[label] = CacheFile.load(filename=filename)
        return cls(files=files, **cfg)
        
    def label_sizes(self):
        sizes = {lb: len(self.files[lb]) for lb in self.labels}
        return sizes
    
    def select(self, label=None, index_list=None):
        return self.files[label].select(index_list=index_list)
    
    def add(self, call=None, label=None):
        self.files[label].add(call=call)
    
class Callset(object):
    Labels = ("SNP_TP", "SNP_FP", "INDEL_TP", "INDEL_FP")
    Cache = None

    def __init__(self, vcf=None, ref=None, labels=None, cache_path=None):
        self.vcf = vcf
        self.ref = ref
        self.labels = labels or self.Labels
        if self.__class__.Cache is None:
            self.__class__.Cache = CallsetCache(labels=self.labels, path=cache_path)
        self.cache = self.__class__.Cache
    
    def get_reference_context(self, call=None, window=None):
        start = call.POS - 1
        end = start + len(call.REF)
        pre = self.ref[call.CHROM][start - window:start]
        post = self.ref[call.CHROM][end:end + window]
        return (pre, post)
    
    def get_alleles(self, call=None, sample_idx=None):
        bases = call.gt_bases[sample_idx]
        if '/' in bases:
            bases = bases.split('/')
        elif '|' in bases:
            bases = bases.split('|')
        else:
            assert 0, "bug"
        assert len(bases) == 2
        return bases

    def get_label(self, call=None):
        if call.is_snp:
            var_type = "SNP"
        elif call.is_indel:
            var_type = "INDEL"
        else:
            return None
        truth = call.INFO.get("CALL")
        if truth not in ('FP', 'TP'):
            return None
        label_name = f"{var_type}_{truth}"
        label_info = {
            "var_type": var_type,
            "truth": truth,
            "label_name": label_name,
            "label": self.labels.index(label_name)
        }
        return label_info

    def process_call(self, call=None, window=None, sample_idx=None):
        (pre, post) = self.get_reference_context(call, window=window)
        bases = self.get_alleles(call, sample_idx=sample_idx)
        alleles = [f"{pre}{var}{post}" for var in bases]
        if set(alleles[0] + alleles[1]) - set('AGTC'):
            return None
        ret = {
            "allele_a": alleles[0],
            "allele_b": alleles[1],
            "chrom": call.CHROM,
            "pos": call.POS,
        }
        return ret

    def load_region(self, region=None, window=96, sample_idx=-1):
        vcf = self.vcf(region)
        counts = {lb: 0 for lb in self.labels}
        for call in vcf:
            label = self.get_label(call)
            if not label:
                continue
            call = self.process_call(call, window=window, sample_idx=sample_idx)
            if not call:
                continue
            call.update(label)
            counts[call["label_name"]] += 1
            self.cache.add(call=call, label=call["label_name"])
        
        n_snp_tp = counts["SNP_TP"]
        n_snp_fp = counts["SNP_FP"]
        n_indel_tp = counts["INDEL_TP"]
        n_indel_fp = counts["INDEL_FP"]
        msg = f"[[ {region} ]] #SNP TP:{n_snp_tp}, FP:{n_snp_fp} #INDEL TP:{n_indel_tp}, FP:{n_indel_fp}"
        print(msg)

    def load_vcf(self, window=96, sample_name="CALLS"):
        re_chr = re.compile("^chr\d+$")
        sample_idx = self.vcf.samples.index(sample_name)
        for region in self.vcf.seqnames:
            if not re_chr.match(region):
                continue
            self.load_region(region=region, window=window)
        print(f"Saving cache")
        self.cache.save()
        return self.cache
    
class GenotypeTokenizer(object):
    def __init__(self, klen=None, max_length=512):
        self.max_length = max_length
        self.klen = klen
        model_path = f"armheb/DNA_bert_{self.klen}"
        self._tokenizer = BertTokenizer.from_pretrained(model_path, do_lower_case=False)
        self.vocab = self._tokenizer.vocab
    
    def kmerize(self, seq):
        assert type(seq) is str
        k_or_mask = lambda kmer: kmer if 'M' not in kmer else '[MASK]'
        kmers = [k_or_mask(seq[i:i + klen]) for i in range(len(seq) - klen + 1)]
        return kmers

    def tokenize(self, row):
        alleles = [row["allele_a"], row["allele_b"]]
        (a1, a2) = [self.kmerize(al) for al in alleles]
        while (len(a1) + len(a2)) > (self.max_length - 3):
            a1 = a1[1:-1]
            a2 = a2[1:-1]
        inp = ["[CLS]"] + a1 + ["[SEP]"] + a2 + ["[SEP]"]
        toks = list(map(self.vocab.__getitem__, inp))
        toks += [0] * (self.max_length - len(toks))
        mask = list(map(bool, toks))
        types = [0] * self.max_length
        ret = {
            "label": row["label"],
            "input_ids": toks,
            "attention_mask": list(map(lambda val: int(bool(val)), toks)),
            "token_type_ids": ([0] * self.max_length)
        }
        return ret

def build_dataset(vcf=None, ref=None, klen=None, window=None, dataset_path=None):
    ds_cache_path = f"{dataset_path}/cache"
    ds_baked_dir = f"{dataset_path}/baked"
    fn_train = f"{ds_baked_dir}/K{klen}-train.json"
    fn_eval = f"{ds_baked_dir}/K{klen}-eval.json"
    if os.path.exists(fn_train) and os.path.exists(fn_eval):
        print("Loading previously baked splits")
        ds_train = Dataset.from_json(fn_train)
        ds_eval = Dataset.from_json(fn_eval)
        return (ds_train, ds_eval)

    try:
        cache = CallsetCache.load(path=ds_cache_path)
        print("Cached callset loaded")
    except:
        cache = None
    if cache is None:
        cs = Callset(vcf=vcf, ref=ref, cache_path=ds_cache_path)
        cache = cs.load_vcf(window=window)
    
    tokenizer = GenotypeTokenizer(klen=klen)
    sizes = cache.label_sizes()
    n_samples = min(sizes.values())
    baked_calls = []
    msg = f"Baking {n_samples} per label ({n_samples * 4} total)"
    print(msg)
    for label in cache.labels:
        msg = f"  - Baking '{label}' at K{klen}"
        print(msg)
        index_list = list(range(sizes[label]))
        random.shuffle(index_list)
        index_list = index_list[:n_samples]
        calls = cache.select(label=label, index_list=index_list)
        calls = list(map(tokenizer.tokenize, calls))
        baked_calls += calls
    
    random.shuffle(baked_calls)
    n_train = int(len(baked_calls) * .9)
    n_eval = (len(baked_calls) - n_train)
    train_calls = baked_calls[:n_train]
    eval_calls = baked_calls[n_train:]
    
    with open(fn_train, 'w') as fh:
        for call in train_calls:
            js = json.dumps(call) + '\n'
            fh.write(js)
    
    with open(fn_eval, 'w') as fh:
        for call in eval_calls:
            js = json.dumps(call) + '\n'
            fh.write(js)
    
    ds_train = Dataset.from_json(fn_train)
    ds_eval = Dataset.from_json(fn_eval)
    return (ds_trian, ds_eval)
