import re

class Cigar(object):
    re_cigar = re.compile('(\d+)([=XID])')

    def __init__(self, cigar=None):
        self.cigar = cigar
        self._parts = None

    @property
    def parts(self):
        if self._parts == None:
            _parts = self.re_cigar.findall(self.cigar)
            self._parts = tuple(map(lambda it: (int(it[0]), it[1]), _parts))
        return self._parts
    
    def __getitem__(self, idx):
        cigar = str.join('', (f'{it[0]}{it[1]}' for it in self.parts[idx]))
        return self.__class__(cigar=cigar)
    
    def __str__(self):
        return self.cigar
    
    def __len__(self):
        return len(self.parts)
    
    def __eq__(self, other):
        return str(self) == str(other)


