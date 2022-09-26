from uuid import uuid4
from Bio import Entrez

Entrez.email = f'{str(uuid4())}@example.com'

def esearch(block_size=1000, limit=None, **kw):
    req = Entrez.esearch(**kw, retmax=1)
    res = Entrez.read(req)
    total = res['Count'] = int(res['Count'])
    if limit:
        total = min(total, limit)
    block_size = min(block_size, total)
    while len(res['IdList']) < total:
        pos = len(res['IdList'])
        assert len(res['IdList']) == len(set(res['IdList']))
        req = Entrez.esearch(**kw, retstart=pos, retmax=block_size)
        res['IdList'] += Entrez.read(req)['IdList']
    return res

def esummary(**kw):
    req = Entrez.esummary(**kw)
    res = Entrez.read(req)
    return res

def efetch(**kw):
    req = Entrez.efetch(**kw)
    res = Entrez.read(req)
    return res
