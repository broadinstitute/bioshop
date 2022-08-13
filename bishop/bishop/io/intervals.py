import fsspec as FS
import pandas as pd
from .. rep.region import Region, RegionList

def load_interval_file_core(path=None, comment=None, index_offset=0, sep='\t', header=('chrom', 'start', 'stop'), astype='regionlist'):
    if astype not in ('region', 'regionlist', 'dataframe'):
        raise TypeError(f'illegal value for astype {astype}')
    intervals = pd.read_csv(path, sep=sep, comment=comment, header=None, names=header)
    intervals.start += index_offset
    intervals.stop += index_offset
    if astype == 'dataframe':
        return intervals
    to_region = lambda row: Region(row.chrom, row.start, row.stop)
    regions = intervals.apply(to_region, axis=1).to_list()
    if astype == 'regions':
        return regions
    region_list = RegionList()
    region_list.add_regions(regions)
    return region_list

def load_picard_interval_file(path=None, **kw):
    header = ('chrom', 'start', 'stop', 'strand', 'target')
    new_kw = dict(comment='@', index_offset=-1, header=header)
    new_kw.update(kw)
    return load_interval_file_core(path=path, **new_kw)

def load_gatk_interval_file(path=None, **kw):
    header = ('chrom', 'start', 'stop')
    new_kw = dict(header=header)
    new_kw.update(kw)
    return load_interval_file_core(path=path, **new_kw)

def load_bed_interval_file(path=None, **kw):
    raise NotImplementedError

def detect_interval_filetype(path=None):
    hints = {'picard': 0, 'bed': 0, 'gatk': 0}
    ext = path.lower().split('.')[-1]
    if ext == 'interval_list':
        hints['picard'] += 1
    elif ext in ('list', 'intervals'):
        hints['gatk'] += 1
    elif ext in ('bed', ):
        hints['bed'] += 1
    fs = FS.open(path, 'r')
    with fs as fh:
        line = fh.readline()
    if line[0] == '@':
        hints['picard'] += 1
    if '\t' in line:
        parts = line.split('\t')
        if len(parts) == 3:
            hints['bed'] += 1
        elif len(parts) == 5:
            hints['picard'] += 1
    else:
        try:
            Region(line)
            hints['gatk'] += 1
        except TypeError:
            pass
    hints = sorted(hints.items(), key=lambda it: it[1])
    return hints[-1][0]

def load_interval_file(path=None, filetype=None, **kw):
    if not filetype:
        filetype = detect_interval_filetype(path=path)
    filetype = filetype.lower()
    if filetype == 'picard':
        return load_picard_interval_file(path=path, **kw)
    elif filetype == 'gatk':
        return load_gatk_interval_file(path=path, **kw)
    elif filetype == 'bed':
        return load_bed_interval_file(path=path, **kw)
    else:
        raise TypeError(f'unknown interval file type')
