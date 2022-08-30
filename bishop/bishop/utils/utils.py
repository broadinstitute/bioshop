import os
import zlib
import requests
import tempfile

import pandas as pd
from tqdm import tqdm
from cachier import cachier

from .. rep.region import Region

__all__ = [
    'download_genome',
    'get_genome_path',
    'get_cache_dir',
    'is_concrete_nucleotides',
    'seq_progress_bar',
    'vcf_progress_bar',
    'region_progress_bar',
    'cache_func',
    'softhash',
    'concat_saved_dataframes',
]

def concat_saved_dataframes(df_list=None):
    frames = [pd.read_pickle(fn) for fn in df_list]
    return pd.concat(frames)

def region_progress_bar(region=None):
    if not isinstance(region, Region):
        region = Region(region)

    pbar = tqdm(
        total=len(region),
        desc=str(region), 
        unit="base",
        unit_scale=True,
        colour='green'
    )
    ns = dict(pbar=pbar, last_pos=region.start)
    def update(pos=None):
        ns['pbar'].update(pos - ns['last_pos'])
        ns['last_pos'] = pos
    return update

def seq_progress_bar(seqlen_map=None):
    ns = dict(chrom=None, pbar=None, last_pos=0)
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

def vcf_progress_bar(vcf=None):
    seq_len_f = lambda it: (it.name, it.length)
    seqlen_map = dict(map(seq_len_f, vcf.header.contigs.values()))
    return seq_progress_bar(seqlen_map)

def get_cache_dir():
    cache_dir = os.environ.get('CACHE_DIR')
    if cache_dir is None:
        user_homedir = os.path.expanduser('~')
        cache_dir = f'{user_homedir}/.cache'
    return cache_dir

def cache_func(*args, **kw):
    cache_dir = get_cache_dir()
    return cachier(*args, cache_dir=cache_dir, **kw)

def download_genome(url, output_fn, chunk_size=None):
    chunk_size = chunk_size or (1 << 20)
    root_dir = os.path.split(os.path.abspath(output_fn))[0]
    (tmp_fd, tmp_path) = tempfile.mkstemp(dir=root_dir)
    tmp_fh = os.fdopen(tmp_fd, mode='wb')
    try:
        with requests.get(url, stream=True) as resp:
            decom = zlib.decompressobj(16 + zlib.MAX_WBITS)
            for chunk in resp.iter_content(chunk_size):
                tmp_fh.write(decom.decompress(chunk))
            tmp_fh.write(decom.flush())
        tmp_fh.close()
        os.rename(tmp_path, output_fn)
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
    return output_fn

def get_genome_path(name=None, url=None, cache_path=None):
    cache_path = cache_path or get_cache_path()
    path_genome = os.path.join(cache_path, name)
    os.makedirs(path_genome, exist_ok=True)
    cache_fn = os.path.join(path_genome, f'{name}.fa')
    if not os.path.exists(cache_fn):
        print(f'Downloading {name}')
        assert download_genome(url=url, output_fn=cache_fn)
    return cache_fn

CONCRETE_NUCLEOTIDES = set('AGTC')

def is_concrete_nucleotides(nucstr):
    nucs = set(str(nucstr).upper())
    return not bool(nucs - CONCRETE_NUCLEOTIDES)

def softhash(content: bytes, salt: int=0, mask: int=0xFFFFFFFF):
    return zlib.adler32(content, salt) & mask

