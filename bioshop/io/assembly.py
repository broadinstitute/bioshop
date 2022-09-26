import os
import re
import zlib
import json
import tempfile
from uuid import uuid4

import requests
import pandas as pd

from .. utils import get_cache_dir
from . entrez import esummary, esearch

def search_for_assembly(query=None, cache=None):
    if cache and query in cache['query']:
        build_key = cache['query'][query]
        return cache['summary'][build_key]

    re_acc = re.compile('^(G..)_(\d\d\d)(\d\d\d)(\d\d\d).*$')
    url_base = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'
    
    idlist = None
    for query_type in ('Organism', 'Accession', 'Assembly'):
        term = f'{query}[{query_type}]'
        res = esearch(db='genome', term=term)
        if res['Count'] > 0:
            idlist = res['IdList']
            break
    if idlist is None:
        return None
    summary = esummary(db='genome', id=idlist[-1])[0]
    asm_acc = summary.get('Assembly_Accession')
    asm_name = summary.get('Assembly_Name')
    m = re_acc.match(asm_acc)
    build_key = f'{asm_acc}_{asm_name}'
    url_parts = list(m.groups()) + [build_key, f'{build_key}_genomic.fna.gz']
    url_asm_path = str.join('/', url_parts)
    url = url_base + url_asm_path
    summary['build_key'] = build_key
    summary['genomic_fna'] = url
    summary['query'] = query
    if cache:
        cache['query'][query] = build_key
        cache['summary'][build_key] = summary
    return summary

def parse_header_from_report(url):
    resp = requests.get(url)
    last_line = None
    for line in resp.iter_lines():
        line = line.decode()
        if line[0] == '#':
            last_line = line
            continue
        break

    assert last_line is not None
    header = last_line[1:].strip().split('\t')
    return header

def download_genome_assembly_metadata(asm_name=None):
    genome_info = search_for_assembly(asm_name)
    assert asm_name == genome_info['Assembly_Name']
    asm_id = genome_info['AssemblyID']
    asc_id = genome_info['Assembly_Accession']
    resp = esummary(db='assembly', id=asm_id, report='full')
    doc = resp['DocumentSummarySet']['DocumentSummary'][0]
    rpt_url = doc['FtpPath_Assembly_rpt'].replace('ftp://', 'https://')
    names = parse_header_from_report(rpt_url)
    df = pd.read_csv(rpt_url, sep='\t', comment='#', names=names)
    units = df.to_dict(orient='records')
    genome_info['Units'] = units
    return genome_info

def download_genome_assembly(url, output_fn, chunk_size=None):
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

def load_assembly(asm_name=None, cache_dir=None):
    cache_dir = cache_dir or get_cache_dir()
    asm_path = os.path.join(cache_dir, 'genome_assemblies', asm_name)
    asm_jsfn = f'{asm_name}.json'
    asm_jsfn = os.path.join(asm_path, asm_jsfn)
    # XXX: save compressed genome!
    asm_fafn = os.path.join(asm_path, f'{asm_name}.fna')
    if not os.path.isfile(asm_jsfn):
        os.makedirs(asm_path, exist_ok=True)
        asm_metadata = download_genome_assembly_metadata(asm_name)
        asm_url = asm_metadata['genomic_fna']
        download_genome_assembly(asm_url, asm_fafn)
        with open(asm_jsfn, 'w') as fh:
            json.dump(asm_metadata, fh, indent=2)
    else:
        with open(asm_jsfn) as fh:
            asm_metadata = json.load(fh)
    # update to point locally
    asm_metadata['local_genomic_fna'] = asm_fafn
    return asm_metadata
