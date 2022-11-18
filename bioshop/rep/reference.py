from pyfaidx import Fasta

class Reference:
    def __init__(self, reference_path=None, reference_index_path=None):
        self.reference_path = reference_path
        self.reference_index_path = reference_index_path
        self._reference = None
        self._seq_cache = {}

    @property
    def reference(self):
        if self._reference == None:
            self._reference = Fasta(
                self.reference_path, 
                indexname=self.reference_index_path
            )
        return self._reference

    def get_sequence(self, name=None):
        if name not in self._seq_cache:
            seq = self.reference[name]
            # sync from disk
            seq = str(seq)
            self._seq_cache[name] = seq
        return self._seq_cache[name]
    
    def __getitem__(self, name):
        return self.get_sequence(name=name)
