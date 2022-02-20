import queue
from pyfaidx import Fasta
from . worker import Worker
from .. models import VariantTokenizer, ModelInputStruct, VariantToVector

class VariantToVectorWorker(Worker):
    def __init__(self, ref_path=None, tokenizer_config=None, **kw):
        super().__init__(**kw)
        self.ref_path = ref_path
        self.tokenizer_config = tokenizer_config
        self.in_q = self.manager.to_vectorizer
        self.to_model = self.manager.to_model
        self.to_gather = self.manager.to_gather

    def dispatch_site(self, site_info=None, gt_inps=None):
        site_id = site_info.site_id
        model_inputs = []
        for (genotype_id, gt_inp) in gt_inps.items():
            if gt_inp is None:
                site_info.genotypes[genotype_id].status = 'skipped'
                continue
            gt_inp["site_id"] = site_info.site_id
            gt_inp["genotype_id"] = genotype_id
            inp = ModelInputStruct(**gt_inp)
            model_inputs.append(inp)
        assert model_inputs or site_info.is_pending == False

        # dispatch
        self.to_gather.put(site_info)
        for inp in model_inputs:
            self.to_model.push(inp)

    def _run(self):
        ref = Fasta(self.ref_path)
        tokenizer = VariantTokenizer(**self.tokenizer_config)
        vectorizer = VariantToVector(ref=ref, tokenizer=tokenizer)

        self._running.set()
        while self.running:
            try:
                site_info = self.in_q.get(timeout=1)
                self.in_q.task_done()
            except queue.Empty:
                if self.manager.flush_vectorizer.is_set():
                    break
                #print(f"{self.__class__.__name__} loop: in_q empty. buf={self.to_model.qsize.value}")
                continue

            gt_inps = vectorizer.process_site(site_info)
            self.dispatch_site(site_info=site_info, gt_inps=gt_inps)
            
