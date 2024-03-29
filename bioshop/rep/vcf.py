from pysam import VariantFile

from . assembly import CuratedAssembly
from . region import Region

class VCF(CuratedAssembly, VariantFile):
    def get_contig_names(self):
        return tuple(self.header.contigs.keys())

    def fetch(self, contig=None, region=None, **kw):
        if contig is None and region is None:
            return super().fetch(**kw)
        #
        assert bool(contig) ^ bool(region)
        if contig:
            region = Region(contig=contig)
        elif type(region) is str:
            region = Region(region)
        region.contig = self.translate_contig_name(region.contig)
        region = str(region)
        return super().fetch(region=region, **kw)
    
    def __reduce__(self):
        kw = self._reduce_hook()
        args = ((self.filename,), kw)
        return (self.__expand__, args)

    @classmethod
    def __expand__(cls, args, kw):
        return cls(*args, **kw)

    # XXX: support append?
    def to_writer(self, path=None):
        return self.__class__(path, 'w', header=self.header)
