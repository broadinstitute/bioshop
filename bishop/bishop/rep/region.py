import re
import math
import portion as P

def interval_cmp(func):
    def cmpfunc(self, other):
        if type(other) is str:
            # try to convert to region
            other = Region(other)
        if isinstance(other, Region):
            if isinstance(self, Region) and self.chrom != other.chrom:
                msg = f'regions have mismatched chrom ({self.chrom} != {other.chrom})'
                raise ValueError(msg)
            other = other.interval
        elif type(other) in (int, float):
            other = P.singleton(other)
        if not isinstance(other, P.Interval):
            raise TypeError(other)
        return func(self, other)
    return cmpfunc

class Region(object):
    re_region = re.compile('(\w+):?(\d+)?-?(\d+)?')

    def __init__(self, chrom=None, start=None, stop=None, contig=None):
        if not (bool(chrom) ^ bool(contig)):
            raise TypeError(f'chrom must be set')
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
            raise TypeError(f'unknown region format')
        vals = [m.groups()[0]] + [int(val) if val is not None else None for val in m.groups()[1:]]
        if start is not None and vals[1] is not None:
            raise TypeError('conflicting values for start')
        start = start if start is not None else vals[1]
        if stop is not None and vals[2] is not None:
            raise TypeError('conflicting values for stop')
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
        self.interval = P.closed(self.interval.lower, val)
    stop = property(get_stop, set_stop)

    def get_start(self):
        return self.interval.lower
    def set_start(self, val):
        self.interval = P.closed(val, self.interval.upper)
    start = property(get_start, set_start)

    def __repr__(self):
        return f"{self.__class__.__name__}('{str(self)}')"
    
    def __str__(self):
        return f'{self.chrom}:{self.interval.lower}-{self.interval.upper}'

    def split(self, step=None):
        for pos in P.iterate(self.interval, step=step):
            yield self.__class__(self.chrom, pos, pos + step)

    def shard(self, n_shards=None):
        step = int(math.ceil(len(self) / n_shards))
        for pos in P.iterate(self.interval, step=step):
            yield self.__class__(self.chrom, pos, pos + step)

    def __len__(self):
        return abs(self.interval.upper - self.interval.lower)

    @interval_cmp
    def __lt__(self, other):
        return self.interval.__lt__(other)

    @interval_cmp
    def __gt__(self, other):
        return self.interval.__gt__(other)

    @interval_cmp
    def __le__(self, other):
        return self.interval.__le__(other)

    @interval_cmp
    def __ge__(self, other):
        return self.interval.__ge__(other)

    @interval_cmp
    def __eq__(self, other):
        return self.interval.__eq__(other)

    @interval_cmp
    def __ne__(self, other):
        return self.interval.__ne__(other)

    @interval_cmp
    def overlaps(self, other):
        return self.interval.overlaps(other)

    @interval_cmp
    def contains(self, other):
        return self.interval.contains(other)

    @interval_cmp
    def __contains__(self, other):
        return self.interval.contains(other)

class RegionList:
    def __init__(self, by_chrom=None):
        if by_chrom is None:
            by_chrom = dict()
        self.by_chrom = by_chrom.copy()

    def add_regions(self, regions=None):
        by_chrom_update = {}
        for region in regions:
            if region.chrom not in by_chrom_update:
                by_chrom_update[region.chrom] = list()
            by_chrom_update[region.chrom].append(region.interval)
        by_chrom_update = {key: P.Interval(*val) for (key, val) in by_chrom_update.items()}
        for chrom in by_chrom_update:
            if chrom not in self.by_chrom:
                self.by_chrom[chrom] = by_chrom_update[chrom]
            else:
                self.by_chrom[chrom] |= by_chrom_update[chrom]

    def contains(self, other):
        if not isinstance(other, Region):
            other = Region(other)
        return other.interval in self.by_chrom.get(other.chrom, P.empty())

    def __contains__(self, other):
        return self.contains(other)

def run_tests():
    r1 = Region('chr1:100-200')
    r2 = Region('chr1:100-200')
    assert r1 == r2
    r2 = Region('chr1:200-300')
    assert r2 >= r1
    assert r1.overlaps(r2)
    r2 = Region('chr1:201-300')
    assert r2 > r1
    assert not r1.overlaps(r2)
    assert 10 < r1
    assert 210 > r1
    assert 110 in r1
    assert r1.contains(110)
    assert not r1.contains(r2)
    r2 = Region('chr1:110-150')
    assert r2 in r1
    assert r1 == 'chr1:100-200'
    assert r1 == Region('chr1:100-200')
    r1 = Region('chr1', 1_000, 5_000)
    rl = RegionList()
    shards = [x for (i,x) in enumerate(r1.shard(6)) if i % 2]
    rl.add_regions(shards)
    assert 'chr1:1700' in rl
    assert 'chr2:1700' not in rl

if __name__ == '__main__':
    run_tests()
