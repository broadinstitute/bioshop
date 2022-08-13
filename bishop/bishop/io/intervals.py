import fsspec as FS
import pandas as pd
from .. rep.region import Region, RegionList, RegionMap, PandasRegionMap

def load_interval_file_core(path=None, name=None, comment=None, index_offset=0, sep='\t', header=('chrom', 'start', 'stop'), astype='regionlist'):
    if astype not in ('region', 'regionlist', 'dataframe'):
        raise TypeError(f'illegal value for astype {astype}')
    if name is None:
        # grab the filename without the extension
        name = path.split('/')[-1].split('.')[0]
    intervals = pd.read_csv(path, sep=sep, comment=comment, header=None, names=header)
    intervals.start += index_offset
    intervals.stop += index_offset
    intervals['name'] = name
    if astype == 'dataframe':
        to_interval = lambda row: pd.Interval(row.start, row.stop, closed='both')
        intervals['interval'] = intervals.apply(to_interval, axis=1)
        return intervals
    to_region = lambda row: Region(row.chrom, row.start, row.stop)
    regions = intervals.apply(to_region, axis=1).to_list()
    if astype == 'regions':
        return regions
    region_list = RegionList(name=name)
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
    if sum(hints.values()) == 0:
        msg = f'unknown interval filetype'
        raise TypeError(msg)
    hints = sorted(hints.items(), key=lambda it: it[1])
    return hints[-1][0]

def load_interval_list(path=None, filetype=None, **kw):
    if not filetype:
        filetype = detect_interval_filetype(path=path)
    filetype = filetype.lower()
    stem = path.split('/')[-1]
    msg = f'Loading intervals from {stem} ({filetype} format)'
    print(msg)
    if filetype == 'picard':
        return load_picard_interval_file(path=path, **kw)
    elif filetype == 'gatk':
        return load_gatk_interval_file(path=path, **kw)
    elif filetype == 'bed':
        return load_bed_interval_file(path=path, **kw)
    else:
        raise TypeError(f'unknown interval file type')

def load_interval_lists(interval_files, astype=None):
    interval_lists = []
    astype = astype or 'regionlist'
    for load_spec in interval_files:
        load_spec = load_spec.copy()
        load_spec['astype'] = astype
        interval_lists.append(
            load_interval_list(**load_spec)
        )
    if astype == 'regionlist':
        interval_map = RegionMap()
        for interval_list in interval_lists:
            interval_map.add_region_list(interval_list)
        return interval_map
    # pandas
    if astype == 'dataframe':
        df = pd.concat(interval_lists)
        df.index = pd.IntervalIndex(df.interval)
        gb = df.groupby('chrom')
        by_chrom = {ch: gb.get_group(ch) for ch in gb.groups}
        return PandasRegionMap(by_chrom=by_chrom)
    # default
    return interval_list
