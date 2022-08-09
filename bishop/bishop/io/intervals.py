from .. rep import Region

def load_interval_file(path=None, comment=None, index_offset=0, sep='\t', header=('chrom', 'start', 'stop'), as_regions=True):
    intervals = pd.read_csv(path, sep=sep, comment=comment, header=None, names=header)
    intervals.start += index_offset
    intervals.stop += index_offset
    if as_regions:
        to_region = lambda row: Region(row.chrom, row.start, row.stop)
        intervals = intervals.apply(to_region, axis=1).to_list()
    return intervals

def load_picard_interval_file(path=None, **kw):
    header = ('chrom', 'start', 'stop', 'strand', 'target')
    new_kw = dict(comment='@', index_offset=-1, header=header)
    new_kw.update(kw)
    return load_interval_file(path=path, **new_kw)
