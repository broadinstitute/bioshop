import re
import math
import portion as P

class Region(object):
    re_region = re.compile('(\w+):?(\d+)?-?(\d+)?')

    def __init__(self, chrom=None, start=None, stop=None, contig=None):
        assert (bool(chrom) ^ bool(contig))
        chrom = chrom or contig
        (self.chrom, start, stop) = \
            self._parse_region(chrom=chrom, start=start, stop=stop)
        if start is None:
            self.interval = P.empty()
        elif stop is None:
            self.interval = P.singleton(start)
        else:
            self.interval = P.closed(start, stop)
    
    def _parse_region(self, chrom=None, start=None, stop=None):
        m = self.re_region.match(chrom)
        if not m:
            raise ValueError(f'unknown region format')
        vals = [m.groups()[0]] + [int(val) if val is not None else None for val in m.groups()[1:]]
        if start is not None and vals[1] is not None:
            raise ValueError('conflicting values for start')
        start = start if start is not None else vals[1]
        if stop is not None and vals[2] is not None:
            raise ValueError('conflicting values for stop')
        stop = stop if stop is not None else vals[2]
        return (
            vals[0],
            start if start is not None else vals[1],
            stop if stop is not None else vals[2],
        )

    def get_contig(self):
        return self.chrom
    def set_contig(self, val):
        self.chrom = val
    contig = property(get_contig, set_contig)

    def get_stop(self):
        return self.interval.upper
    def set_stop(self, val):
        #val = max(self.interval.lower, val)
        self.interval = P.closed(self.interval.lower, val)
    stop = property(get_stop, set_stop)

    def get_start(self):
        return self.interval.lower
    def set_start(self, val):
        #val = min(val, self.interval.upper)
        self.interval = P.closed(val, self.interval.upper)
    start = property(get_start, set_start)

    def __repr__(self):
        return f"{self.__class__.__name__}('{str(self)}')"
    
    def __str__(self):
        return f'{self.chrom}:{self.interval.lower}-{self.interval.upper}'

    def __len__(self):
        return len(self.interval)

    def split(self, step=None):
        for pos in P.iterate(self.interval, step=step):
            yield self.__class__(self.chrom, pos, pos + step)

    def shard(self, n_shards=None):
        step = int(math.ceil(len(self) / n_shards))
        for pos in P.iterate(self.interval, step=step):
            yield self.__class__(self.chrom, pos, pos + step)

