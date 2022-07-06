#!/usr/bin/env python

import ctypes
import time
from pprint import pprint
from queue import Full, Empty
from math import ceil
from collections import defaultdict
from pyfaidx import Fasta
from transformers import BertTokenizer, BertForMaskedLM
from transformers.utils.logging import set_verbosity_error, set_verbosity_info
import multiprocessing as mp
import numpy as np
import torch
import torch.nn.functional as F
from tqdm import tqdm
from cyvcf2 import VCF, Writer
import click

Status_PENDING = 1
Status_SKIPPED = 2
Status_CALLED = 3
Status_PAD = 4

def default_model_path(klen=None):
    assert klen >= 3 and klen <= 6
    return f'armheb/DNA_bert_{klen}'

def is_dna(seq=None):
    return not bool(set(seq.upper()) - set('AGTC'))

def load_model(model_path=None, klen=None, head_class=BertForMaskedLM):
    model_path = model_path or default_model_path(klen=klen)
    set_verbosity_error()
    model = head_class.from_pretrained(model_path)
    set_verbosity_info()
    return model

def load_model_tokenizer(model_path=None, klen=None):
    model_path = model_path or default_model_path(klen=klen)
    set_verbosity_error()
    tokenizer = BertTokenizer.from_pretrained(model_path, do_lower_case=False)
    set_verbosity_info()
    return tokenizer

class BaseStruct(ctypes.Structure):
    @classmethod
    def from_numpy(cls, **kw):
        c_cast = lambda it: np.ctypeslib.as_ctypes(np.squeeze(it))
        kw = {key: c_cast(val) for (key, val) in kw.items()}
        #pprint({key: type(val) for (key, val) in kw.items()})
        return cls(**kw)

    def as_numpy(self):
        to_np = lambda attr: np.ctypeslib.as_array(getattr(self, attr))
        arrays = (to_np(name) for (name, _) in self._fields_)
        return tuple(arrays)

    def as_dict(self):
        arrays = self.as_numpy()
        keys = [field[0] for field in self._fields_]
        ret = dict(zip(keys, arrays))
        return ret

def get_batch_struct(batch_size=1, sequence_length=512):
    _fields_ = [
        ('input_ids', ctypes.c_int32 * sequence_length * batch_size), 
        ('attention_mask', ctypes.c_int32 * sequence_length * batch_size),
        ('token_type_ids', ctypes.c_int32 * sequence_length * batch_size),
        ('label_ids', ctypes.c_int32 * sequence_length * batch_size),
        ('status', ctypes.c_int64 * batch_size),
        ('allele_hash', ctypes.c_int64 * batch_size),
    ]
    return type("BatchStruct", (BaseStruct,), {'_fields_': _fields_})

class CircularBuffer(object):
    def __init__(self, ctype=None, size=None):
        self.ctype = ctype
        self.size = size
        self.buffer = mp.RawArray(self.ctype, self.size)
        self._read_ptr = mp.RawValue(ctypes.c_uint, 0)
        self._write_ptr = mp.RawValue(ctypes.c_uint, 0)
        self._qsize = mp.Value(ctypes.c_uint, 0)
        self.buffer_cv = mp.Condition()

    def get_qsize(self):
        return self._qsize.value
    def set_qsize(self, val):
        self._qsize.value = val
    qsize = property(get_qsize, set_qsize)

    def get_read_ptr(self):
        return self._read_ptr.value
    def set_read_ptr(self, val):
        self._read_ptr.value = val
    read_ptr = property(get_read_ptr, set_read_ptr)

    def get_write_ptr(self):
        return self._write_ptr.value
    def set_write_ptr(self, val):
        self._write_ptr.value = val
    write_ptr = property(get_write_ptr, set_write_ptr)

    @property
    def full(self):
        return self.qsize == self.size
    
    @property
    def empty(self):
        return self.qsize == 0

    def put(self, item, timeout=None):
        _timeout = timeout or 1
        with self.buffer_cv:
            while self.full:
                timeout_flag = not self.buffer_cv.wait(timeout=_timeout)
                if timeout and timeout_flag:
                    raise Full
            self.buffer[self.write_ptr] = item
            self.write_ptr = (self.write_ptr + 1) % self.size
            self.qsize += 1
            self.buffer_cv.notify_all()

    def get(self, timeout=None):
        _timeout = timeout or 1
        with self.buffer_cv:
            while self.empty:
                timeout_flag = not self.buffer_cv.wait(timeout=_timeout)
                if timeout and timeout_flag:
                    raise Empty
            src = self.buffer[self.read_ptr]
            # critical: copy the data before passing the pointer
            item = type(src)()
            ctypes.pointer(item)[0] = src
            self.read_ptr = (self.read_ptr + 1) % self.size
            self.qsize -= 1
            self.buffer_cv.notify_all()
        return item

class Worker(mp.Process):
    def __init__(self):
        super().__init__()
        self.daemon = True
        self._flush = mp.Value(ctypes.c_bool, 0)
        self._running = mp.Value(ctypes.c_bool, 0)

    def get_flush(self):
        return self._flush.value
    def set_flush(self, value):
        self._flush.value = int(value)
    flush = property(get_flush, set_flush)

    def get_running(self):
        return self._running.value
    def set_running(self, value):
        self._running.value = int(value)
    running = property(get_running, set_running)

def vcf_progress_bar(vcf=None):
    seqlen_map = dict(zip(vcf.seqnames, vcf.seqlens))
    ns = dict(
        chrom=None,
        pbar=None,
        last_pos=0
    )
    def update(chrom=None, pos=None):
        pbar = ns['pbar']
        if pbar and chrom is None:
            pbar.close()
            return
        last_chrom = ns['chrom']
        if last_chrom != chrom:
            seqlen = seqlen_map[chrom]
            if pbar:
                pbar.close()
            pbar = tqdm(
                total=seqlen,
                desc=chrom,
                unit="base",
                unit_scale=True,
                colour='green'
            )
            ns['chrom'] = chrom
            ns['pbar'] = pbar
            ns['last_pos'] = 0
        last_pos = ns['last_pos']
        pbar.update(pos - last_pos)
        ns['last_pos'] = pos
    return update

class Gather(Worker):
    def __init__(self, path_input_vcf=None, region=None, path_output_vcf=None, in_q=None):
        super().__init__()
        self.path_input_vcf = path_input_vcf
        self.path_output_vcf = path_output_vcf
        self.region = region
        self.in_q = in_q
        self.results = {}

    def transpose_batch(self, batch):
        keys = list(batch.keys())
        ret = [dict(zip(keys, it)) for it in zip(*[batch[key].tolist() for key in keys])]
        return ret

    def wait_on(self, allele_hash=None, timeout=.1):
        while self.running:
            try:
                results = self.in_q.get(timeout=timeout)
                if 'loss' in results:
                    results = self.transpose_batch(results)
                else:
                    results = [results]
                for item in results:
                    if item['status'] == Status_PAD:
                        continue
                    self.results[item['allele_hash']] = item
                if (allele_hash is None) or (allele_hash in self.results):
                    return True
            except Empty:
                pass
        return False

    def score_alleles(self, site=None):
        alleles = [site.REF] + list(site.ALT)
        results = {}
        for allele in alleles:
            allele_hash = get_allele_hash(site, allele)
            if allele_hash not in self.results:
                ok = self.wait_on(allele_hash=allele_hash)
                if not ok:
                    msg = 'Failed to wait on allele hash'
                    raise ValueError(msg)
            assert allele_hash in self.results
            res = self.results.pop(allele_hash)
            results[allele] = res
        #
        if results[site.REF]['status'] == Status_SKIPPED:
            scores = {alt: -1 for alt in alleles}
        else:
            scores = {}
            ref_loss = results[site.REF]['loss']
            for allele in alleles:
                if results[allele]['status'] == Status_SKIPPED:
                    scores[allele] = -1
                else:
                    #scores[allele] = results[allele]['loss'] / ref_loss
                    scores[allele] = results[allele]['loss']
        return tuple([scores[al] for al in alleles])

    def run(self):
        self.running = True

        vcf_in = VCF(self.path_input_vcf)
        pbar = vcf_progress_bar(vcf_in)
        vcf_in.add_info_to_header({'ID': 'BLOSS', 'Description': 'BERT Loss', 'Type': 'Float', 'Number': 'R'})
        vcf_out = Writer(self.path_output_vcf, vcf_in)

        if self.region:
            vcf_in = vcf_in(self.region)

        try:
            for site in vcf_in:
                if not self.running:
                    break
                pbar(site.CHROM, site.POS)
                site.INFO['BLOSS'] = self.score_alleles(site)
                vcf_out.write_record(site)
        finally:
            vcf_out.close()
            vcf_in.close()

class ModelRunner(Worker):
    def __init__(self, model_path=None, klen=None, batch_size=None, qsize=64):
        super().__init__()
        self.klen = klen
        self.model_path = model_path
        self.batch_size = batch_size
        mt = load_model_tokenizer(klen=self.klen)
        self.mask_token_id = mt.mask_token_id
        self.batch_struct = get_batch_struct(batch_size=batch_size)
        self.in_q = CircularBuffer(ctype=self.batch_struct, size=qsize)
        self.out_q = mp.Queue()

    def run(self):
        self.running = True
        model = load_model(model_path=self.model_path, klen=self.klen).cuda().eval()
        while self.running:
            try:
                batch = self.in_q.get(timeout=0.1)
                batch = batch.as_dict()
            except Empty:
                if self.flush:
                    break
                continue
            input_ids = torch.from_numpy(batch.pop('input_ids')).cuda()
            attention_mask = torch.from_numpy(batch.pop('attention_mask')).cuda()
            token_type_ids = torch.from_numpy(batch.pop('token_type_ids')).cuda()
            label_ids = torch.from_numpy(batch.pop('label_ids')).cuda()
            with torch.no_grad():
                outp = model(
                    input_ids=input_ids,
                    attention_mask=attention_mask,
                    token_type_ids=token_type_ids,
                )

                label_mask = (input_ids == self.mask_token_id)
                labels = (~label_mask * -100) + (label_mask * label_ids)
                labels = labels.to(outp.logits.device)
                loss = F.cross_entropy(outp.logits.permute(0, 2, 1), labels, reduction='none')
                loss = torch.sum(loss, axis=-1) / torch.sum(label_mask, axis=-1)
                batch['loss'] = loss.cpu()
            #pprint({key: (batch[key].shape, batch[key].dtype) for key in batch})
            self.out_q.put(batch)

class ReferenceLookup(object):
    def __init__(self, path_ref=None, window_size=300):
        self.path_ref = path_ref
        self.ref = Fasta(self.path_ref)
        self.window_size = window_size
        self._cache = {}

    def lookup(self, site):
        if site.CHROM not in self._cache:
            self._cache[site.CHROM] = str(self.ref[site.CHROM]).upper()
        seq = self._cache[site.CHROM]
        var_start = site.POS - 1
        var_end = var_start + len(site.REF)
        up = seq[var_start - self.window_size:var_start]
        down = seq[var_end:var_end + self.window_size]
        assert seq[var_start - self.window_size:var_end + self.window_size] == (up + site.REF + down)
        return (up, down)

    __call__ = lookup

class Masker(object):
    def __init__(self, mask_length=None, mask_pos=None, mask_value='M'):
        self.mask_pos = mask_pos
        self.mask_length = mask_length
        self.mask_value = mask_value
        self.mask_str = self.mask_value * self.mask_length

    def mask_seq(self, up=None, down=None):
        if self.mask_pos in ('up', 'both'):
            up = up[:-self.mask_length] + self.mask_str
        if self.mask_pos in ('down', 'both'):
            down = self.mask_str + down[self.mask_length:]
        return (up, down)

    __call__ = mask_seq

class Batcher(object):
    def __init__(self, batch_size=None):
        self.batch_size = batch_size

    def transpose_batch(self, batch=None):
        batch_t = defaultdict(list)
        for row in batch:
            for key in row:
                batch_t[key].append(row[key])

        for key in ('input_ids', 'attention_mask', 'token_type_ids'):
            batch_t[key] = np.vstack(batch_t[key])

        for key in ('status', 'allele_hash'):
            batch_t[key] = np.array(batch_t[key])
        return dict(batch_t)

    def batch_itr(self, itr=None):
        batch = []
        for item in itr:
            batch.append(item)
            if len(batch) == self.batch_size:
                yield self.transpose_batch(batch=batch)
                batch = []
        if batch:
            pad_len = self.batch_size - len(batch)
            tmpl = batch[0].copy()
            tmpl['status'] = Status_PAD
            tmpl['allele_hash'] = -1
            batch += [tmpl.copy() for _ in range(self.batch_size - len(batch))]
            yield self.transpose_batch(batch=batch)

    __call__ = batch_itr

class Tokenizer(object):
    Token_CLS = '[CLS]'
    Token_SEP = '[SEP]'
    Token_MASK = '[MASK]'

    def __init__(self, klen=None):
        self.klen = klen
        mt = load_model_tokenizer(klen=klen)
        self.vocab = mt.vocab.copy()
        self.max_length = mt.model_max_length

    def kmerize(self, seq):
        assert type(seq) is str
        k_or_mask = lambda kmer: kmer if 'M' not in kmer else self.Token_MASK
        kmers = [k_or_mask(seq[i:i + self.klen]) for i in range(len(seq) - self.klen + 1)]
        return kmers

    def tokenize(self, seq=None):
        kmers = self.kmerize(seq=seq)
        if len(kmers) > (self.max_length - 2):
            rem = ceil((len(kmers) - (self.max_length - 2)) / 2)
            kmers = kmers[rem:-rem]
        inp = [self.Token_CLS] + kmers + [self.Token_SEP]
        input_ids = list(map(self.vocab.__getitem__, inp))
        input_ids += [0] * (self.max_length - len(input_ids))
        assert len(input_ids) == self.max_length
        input_ids = np.array(input_ids, dtype=ctypes.c_int32)
        ret = {
            "input_ids": input_ids,
            "attention_mask": np.array(input_ids != 0, dtype=ctypes.c_int32),
            "token_type_ids": np.zeros_like(input_ids, dtype=ctypes.c_int32),
        }
        return ret

    __call__ = tokenize

def get_allele_hash(site=None, allele=None):
    key = f'{site.CHROM}:{site.POS}#{allele}'
    return hash(key)

def gather_vcf(path_vcf_out=None, region=None, query_func=None):
    vcf_in = VCF(path_input_vcf)
    # XXX: add interval file support
    if self.region:
        vcf_in = vcf_in(self.region)

    for (site_id, site) in enumerate(vcf_in):
        result = query_func(site_id=site_id)

def iter_vcf(path_input_vcf=None, region=None):
    vcf_in = VCF(path_input_vcf)
    # XXX: add interval file support
    if region:
        vcf_in = vcf_in(region)

    for (site_id, site) in enumerate(vcf_in):
        yield (site_id, site)

def build_inputs(vcf_itr=None, ref_lookup=None, masker=None, tokenizer=None, out_q=None):
    for (site_id, site) in vcf_itr:
        tmpl = {'status': Status_PENDING}
        (ref_up, ref_down) = ref_lookup(site)
        (mask_up, mask_down) = masker(ref_up, ref_down)
        if not (is_dna(ref_up) and is_dna(ref_down)):
            tmpl['status'] = Status_SKIPPED
        assert is_dna(site.REF)
        alleles = list(site.ALT) + [site.REF]
        for allele in alleles:
            inp = tmpl.copy()
            inp['allele_hash'] = get_allele_hash(site, allele)
            if not is_dna(allele):
                inp['status'] = Status_SKIPPED
            if inp['status'] == Status_SKIPPED:
                out_q.put(inp)
            elif inp['status'] == Status_PENDING:
                allele_seq = mask_up + allele + mask_down
                label_seq = ref_up + allele + ref_down
                inp.update(tokenizer(allele_seq))
                inp['label_ids'] = tokenizer(label_seq)['input_ids']
                yield inp
            else:
                raise ValueError(inp['status'])

def run_newt(
    path_input_vcf=None,
    path_input_ref=None,
    path_output_vcf=None,
    region=None,
    klen=None,
    mask_length=None,
    mask_pos=None,
    batch_size=None
):
    vcf_itr = iter_vcf(path_input_vcf=path_input_vcf, region=region)
    ref_lookup = ReferenceLookup(path_ref=path_input_ref)
    masker = Masker(mask_length=mask_length, mask_pos=mask_pos)
    tokenizer = Tokenizer(klen=klen)
    batcher = Batcher(batch_size=batch_size)
    model = ModelRunner(klen=klen, batch_size=batch_size)
    gather = Gather(
        path_input_vcf=path_input_vcf, 
        path_output_vcf=path_output_vcf, 
        in_q=model.out_q
    )
    inputs = build_inputs(
        vcf_itr=vcf_itr, 
        ref_lookup=ref_lookup, 
        masker=masker, 
        tokenizer=tokenizer, 
        out_q=gather.in_q
    )
    batches = batcher(inputs)

    try:
        model.start()
        gather.start()
        for batch in batches:
            batch = model.in_q.ctype.from_numpy(**batch)
            model.in_q.put(batch)
    except:
        # crash exit
        print("CRASH!")
        model.running = False
        gather.running = False
        model.join()
        gather.join()
        raise
    # normal exit
    model.flush = True
    model.join()
    gather.join()

@click.command()
@click.option('--vcf', required=True, help='Path to input VCF')
@click.option('--ref', required=True, help='Path to VCF Reference File')
@click.option('--out', help='Path to output VCF')
@click.option('--region', help='Region to call')
@click.option('--klen', default=6, help='kmer length')
@click.option('--mask_length', default=13, help='mask length in bases')
@click.option('--mask_pos', default='up', help='mask position (up, down, both)')
@click.option('--klen', default=6, help='kmer length')
@click.option('--batch_size', default=16, help='Batch size')
def main(**kw):
    if kw['out'] is None:
        fn_inp = kw['vcf'].split('/')[-1]
        parts = fn_inp.split('.')
        parts[0] += '_mostly'
        fn_out = str.join('.', parts)
        kw['out'] = fn_out
    kw['path_input_vcf'] = kw.pop('vcf')
    kw['path_output_vcf'] = kw.pop('out')
    kw['path_input_ref'] = kw.pop('ref')
    run_newt(**kw)

if __name__ == '__main__':
    main()
