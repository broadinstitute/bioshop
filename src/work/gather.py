import os
import re
import time
import queue
from cyvcf2 import VCF, Writer
from . worker import Worker
from .. utils import vcf_progress_bar
from .. base import Site, Genotype

class GatherWorker(Worker):
    def __init__(self, vcf_in_path=None, vcf_out_path=None, region=None, **kw):
        super().__init__(**kw)
        self.vcf_in_path = vcf_in_path
        self.vcf_out_path = vcf_out_path
        self.in_q = self.manager.to_gather
        self.region = region
        self.site_cache = {}
        self.results_cache = {}

    def process_item(self, item=None):
        site_id = item.site_id
        if isinstance(item, Site):
            if site_id not in self.site_cache:
                self.site_cache[site_id] = item
            else:
                site = self.site_cache[site_id]
                site.update(item)
            return
        if isinstance(item, Genotype):
            # XXX: place holder for place holder
            site = self.site_cache.get(site_id)
            site.update_genotype(item)
            return
        raise ValueError(type(item))

    def wait_on(self, site_id=None):
        while True:
            site = self.site_cache.get(site_id)
            if site and not site.is_pending:
                return self.site_cache.pop(site_id)

            try:
                info = self.in_q.get(timeout=1)
            except queue.Empty:
                #print(f"{self.__class__.__name__} wait_in(#{site_id}): in_q empty. outstanding={len(self.site_cache)}")
                #pending = [x.site_id for x in self.site_cache.values() if x.is_pending]
                #print(pending)
                continue
            self.in_q.task_done()
            if type(info) == list:
                for item in info:
                    self.process_item(item)
            else:
                self.process_item(info)

    def _run(self):
        vcf_in = VCF(self.vcf_in_path)
        pbar = vcf_progress_bar(vcf_in)

        vcf_in.add_info_to_header({'ID': 'BLOD', 'Description': 'BERT LOD', 'Type':'Float', 'Number': '1'})
        vcf_in.add_filter_to_header({'ID': 'BERT', 'Description': 'mostly dnabert'})
        vcf_in.add_format_to_header({'ID': 'BT', 'Description': 'BERT LOD', 'Type':'Float', 'Number': '1'})
        vcf_out = Writer(self.vcf_out_path, vcf_in)

        if self.region:
            vcf_in = vcf_in(self.region)

        try:
            for (site_id, site) in enumerate(vcf_in):
                site_info = self.wait_on(site_id=site_id)
                assert site_id == site_info.site_id
                assert site.POS == site_info.pos
                assert site.CHROM == site_info.chrom
                site = site_info.call_site(site)
                vcf_out.write_record(site)
                pbar(site.CHROM, site.POS)
                #print(site_info)
        finally:
            vcf_out.close()
            vcf_in.close()

        self.manager.flush_model.set()
        self.manager.flush_vectorizer.set()
