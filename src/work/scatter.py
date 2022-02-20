from .. base import Site
from . worker import Worker
from cyvcf2 import VCF

class ScatterWorker(Worker):
    def __init__(self, vcf_in_path=None, vcf_idx_path=None, region=None, **kw):
        super().__init__(**kw)
        self.vcf_in_path = vcf_in_path
        self.vcf_idx_path = vcf_idx_path
        self.to_vectorizer_que = self.manager.to_vectorizer
        self.to_gather_que = self.manager.to_gather
        self.region = region

    def push_site(self, site=None, site_id=None):
        site_info = Site.load_from_site(site=site, site_id=site_id)
        if not (site.is_snp or site.is_indel):
            site_info.set_site_status("skipped")
            self.to_gather_que.put(site_info)
            return
        """
        if site.FILTER:
            site_info.set_site_status("skipped")
            self.to_gather_que.put(site_info)
            return
        """
        self.to_vectorizer_que.put(site_info)
    
    def _run(self):
        vcf_in = VCF(self.vcf_in_path)
        if self.vcf_idx_path:
            vcf_in.set_index(self.vcf_idx_path)
        if self.region:
            vcf_in = vcf_in(self.region)

        for (site_id, site) in enumerate(vcf_in):
            self.push_site(site=site, site_id=site_id)

        self.to_vectorizer_que.join()
        self.to_gather_que.join()
