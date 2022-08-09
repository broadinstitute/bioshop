from pprint import pprint
import json

from .. io.assembly import load_assembly

from pyfaidx import Fasta

class CuratedAssembly:
    def __init__(self, *args, contig_names=None, metadata=None, scheme=None, ignore_missing=False, **kw):
        super().__init__(*args, **kw)
        self.metadata = metadata
        if scheme is None:
            if contig_names is None:
                contig_names = self.get_contig_names()
            scheme = self.metadata.detect_scheme(contig_names, ignore_missing=ignore_missing)
        self.scheme = scheme

    def get_contig_names(self):
        raise NotImplementedError

    def translate_contig_name(self, name=None):
        return self.translate_contig_name_as(name=name, as_scheme=self.scheme)

    def translate_contig_name_as(self, name=None, as_scheme=None):
        if self.metadata:
            as_scheme = as_scheme or self.scheme
            return self.metadata.as_scheme(name, as_scheme=as_scheme)
        return name

class GenomeAssembly(CuratedAssembly):
    def __init__(self, fnfa=None, **kw):
        self.fnfa = fnfa
        self.seqs = Fasta(self.fnfa)
        names = tuple([seq.name for seq in self.seqs])
        super().__init__(contig_names=names, **kw)

    def get_sequence(self, name=None, as_type=str):
        name = self.translate_contig_name(name=name)
        seq = self.seqs[name]
        if as_type:
            seq = as_type(seq)
        return seq
    
    def __getitem__(self, name):
        return self.get_sequence(name=name)

class GenomeAssemblyMetadata:
    def __init__(self, name=None, organism=None, accession=None, genomic_fna=None, units=None):
        self.name = name
        self.organism = organism
        self.accession = accession
        self.genomic_fna = genomic_fna
        self.units = units
        self._genome = None

    @property
    def genome(self):
        if self._genome is None:
            assert self.genomic_fna is not None
            self._genome = GenomeAssembly(self.genomic_fna, metadata=self)
        return self._genome

    def as_scheme(self, name=None, as_scheme='ucsc'):
        if name is None:
            return None
        for unit in self.units:
            if name in unit:
                return unit.as_scheme(as_scheme)
        raise KeyError(name)

    def detect_scheme(self, seq_names=None, ignore_missing=False):
        find_scheme = lambda nm: [unit.detect_scheme(nm) for unit in self.units if nm in unit]
        schemes = [find_scheme(nm) for nm in seq_names]
        schemes = [it[0] for it in schemes if it]
        if len(set(schemes)) != 1:
            msg = f'Ambiguous sequence name'
            raise ValueError(msg)
        if (not ignore_missing) and (len(schemes) != len(seq_names)):
            n_missing = len(seq_names) - len(schemes)
            msg = f'Could not identify naming scheme for {n_missing} sequence names'
            raise ValueError(msg)
        return schemes[0]

    def contig_map(self, **kw):
        return ContigMapper(metadata=self, **kw)

    @classmethod
    def load_from_data(cls, data=None):
        acc = data.get('Assembly_Accession')
        name = data.get('Assembly_Name')
        org = data.get('Organism_name')
        genomic_fna = data.get('local_genomic_fna')
        units = data.get('Units', list())
        units = [GenomeAssemblyUnit.load_from_data(data=unit) for unit in units]
        return cls(name=name, organism=org, accession=acc, genomic_fna=genomic_fna, units=units)

    @classmethod
    def load(cls, name=None):
        asm_info = load_assembly(asm_name=name)
        return cls.load_from_data(data=asm_info)

class GenomeAssemblyUnit:
    def __init__(self, name=None, length=None, role=None, type=None, assigned=None, aliases=None):
        self.name = name
        self.length = length
        self.role = role
        self.type = type
        self.aliases = aliases
        self.alias_to_scheme = {val:key for (key, val) in self.aliases.items()}
        self.assigned = assigned

    def __contains__(self, name):
        return name in self.alias_to_scheme

    def as_scheme(self, scheme_name=None):
        return self.aliases[scheme_name]

    def detect_scheme(self, name=None):
        if name not in self:
            return KeyError(name)
        return self.alias_to_scheme[name]

    @classmethod
    def load_from_data(cls, data=None):
        assigned = data.get('Assigned-Molecule')
        name = data.get('Sequence-Name')
        role = data.get('Sequence-Role')
        _type = data.get('Assigned-Molecule-Location/Type')
        length = int(data.get('Sequence-Length', -1))
        ns_map = (
            ('ncbi', 'Sequence-Name'), 
            ('genbank', 'GenBank-Accn'),
            ('refseq', 'RefSeq-Accn'),
            ('ucsc', 'UCSC-style-name'),
        )
        aliases = {}
        for (ns_name, col_name) in ns_map:
            alias = data.get(col_name)
            alias = None if alias == 'na' else alias
            aliases[ns_name] = alias
        return cls(
            name=name,
            length=length,
            role=role,
            type=_type,
            assigned=assigned,
            aliases=aliases,
        )

def load_genome_assembly(jsfn=None):
    with open(jsfn) as fh:
        data = json.load(fh)
    asm = GenomeAssemblyMetadata.load_from_data(data=data)
    return asm
