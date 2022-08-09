import os
import zlib
import requests
import tempfile

__all__ = [
    'download_genome',
    'get_genome_path',
    'get_cache_dir',
    'is_concrete_nucleotides'
]

def get_cache_dir():
    cache_dir = os.environ.get('CACHE_DIR')
    if cache_dir is None:
        user_homedir = os.path.expanduser('~')
        cache_dir = f'{user_homedir}/.cache'
    return cache_dir

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

