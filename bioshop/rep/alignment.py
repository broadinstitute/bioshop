from dataclasses import dataclass
import edlib
from . cigar import Cigar

@dataclass
class LiteralFingerprint:
    chrom: str
    pos: int
    ref: str
    alt: str

    def __str__(self):
        return f'{self.chrom}:{self.pos}#{self.ref}>{self.alt}'

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(str(self))
    
    def __eq__(self, other):
        return str(self) == str(other)


@dataclass
class CigarFingerprint:
    chrom: str
    pos: int
    cigar: Cigar

    def __str__(self):
        return f'{self.chrom}:{self.pos}#{self.cigar}'

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(str(self))
    
    def __eq__(self, other):
        return str(self) == str(other)


class BaseAlignment(object):
    ModeMap = {'global': 'NW', 'infix': 'HW', 'prefix': 'SHW'}

    def __init__(self, mode='global'):
        super().__init__()
        self._al = None
        self.mode = mode

    def __getstate__(self):
        return (self._al, self.mode)
    
    def __setstate__(self, state):
        (self._al, self.mode) = state
    
    @property
    def alignment(self):
        if self._al is None:
            # infix
            mode_name = self.ModeMap[self.mode]
            self._al = edlib.align(self.target, self.query, task='path', mode=mode_name)
        return self._al

    @property
    def locations(self):
        return self.alignment['locations']

    @property
    def pretty_alignment(self):
        al = edlib.getNiceAlignment(self.alignment, self.target, self.query)
        al_pretty = '{query_aligned}\n{matched_aligned}\n{target_aligned}'
        return al_pretty.format(**al)

    @property
    def cigar(self):
        return Cigar(self.alignment['cigar'])
    
class RefAltAlignment(BaseAlignment):
    def __init__(self, allele=None):
        super().__init__()
        self.allele = allele

    def __getstate__(self):
        state = super().__getstate__()
        return (state, self.allele)

    def __setstate__(self, state):
        (state, self.allele) = state
        super().__setstate__(state)
    
    @property
    def target(self):
        return self.allele.refseq

    @property
    def query(self):
        return self.allele.altseq

class AltAltAlignment(BaseAlignment):
    def __init__(self, query_allele=None, target_allele=None):
        super().__init__(mode='infix')
        self.query_allele = query_allele
        self.target_allele = target_allele
    
    def __getstate__(self):
        state = super().__getstate__()
        return (state, self.query_allele, self.target_allele)

    def __setstate__(self, state):
        (state, self.query_allele, self.target_allele) = state
        super().__setstate__(state)
    
    @property
    def target(self):
        return self.target_allele.altseq

    @property
    def query(self):
        return self.query_allele.altseq
    
    @property
    def is_match(self):
        if self.locations[0][0] != 0:
            return False
        min_match = min(len(self.target), len(self.query))
        return self.cigar.parts[0] == (min_match, '=')

    def debug_match(self):
        (alt_start, alt_end) = self.target_allele.alt_span
        alt_len = len(self.target_allele.alt)
        target_alt_str = (' ' * alt_start) + ('*' * alt_len)
        #
        (alt_start, alt_end) = self.query_allele.alt_span
        alt_len = len(self.query_allele.alt)
        query_alt_str = (' ' * alt_start) + ('*' * alt_len)

        print(f'target: {self.target_allele.literal_fingerprint}')
        print(f'query: {self.query_allele.literal_fingerprint}')
        print(f'cigar: {self.cigar} match: {self.is_match}')
        print(target_alt_str)
        print(self.pretty_alignment)
        print(query_alt_str)
