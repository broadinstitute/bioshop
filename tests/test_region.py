from bioshop.rep.region import Region

def test_shard():
    reg = Region('chr3:50000000-55000000')
    rl = list(reg.shard(8))
    assert len(rl) == 8
    assert set(map(len, rl)) == set([625000])

if __name__ == '__main__':
    test_shard()
