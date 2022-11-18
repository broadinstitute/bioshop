from bioshop.rep.region import Region

def test_shard():
    reg = Region('chr3:50000000-55000000')
    rl = list(reg.shard(8))
    assert len(rl) == 8
    assert set(map(len, rl)) == set([625000])

def test_window():
    windows = 1024
    width = 4096
    overlap = 512
    step = width - overlap
    length = windows * step
    chrom = 'chr3'
    reg = Region(f'{chrom}:0-{length}')
    rl = list(reg.window(width=width, step=step))
    assert len(rl) == windows
    assert rl[0].start == reg.start
    assert rl[-1].stop == reg.stop
    assert set([chrom]) == set([r.chrom for r in rl])
