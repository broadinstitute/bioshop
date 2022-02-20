import multiprocessing.managers
import multiprocessing as mp
from . model import ModelWorker
from . scatter import ScatterWorker
from . gather import GatherWorker
from . vector import VariantToVectorWorker
from . cbuf import CircularBuffer

def init_manager(batch_size=None, factor=4):
    manager = mp.managers.SyncManager()
    manager.start()
    maxsize = int(round(batch_size * factor))
    manager.to_vectorizer = manager.JoinableQueue(maxsize=maxsize)
    manager.flush_vectorizer = manager.Event()
    manager.to_model = CircularBuffer(ctype=ModelInputStruct, size=maxsize)
    manager.flush_model = manager.Event()
    manager.to_gather = manager.JoinableQueue()
    return manager

def main(ref_path=None, vcf_in_path=None, vcf_idx_path=None, vcf_out_path=None, model_path=None, batch_size=None, klen=None, window=96, region=None, n_workers=1):
    manager = init_manager(batch_size=batch_size)
    tokenizer_config = dict(klen=klen)
    vtv_init = lambda: VariantToVectorWorker(
        manager=manager, 
        ref_path=ref_path, 
        tokenizer_config=tokenizer_config, 
        window=window
    )
    vtv_workers = [vtv_init() for x in range(n_workers)]
    model_worker = ModelWorker(manager=manager, model_path=model_path, batch_size=batch_size, klen=klen)
    scatter_worker = ScatterWorker(manager=manager, vcf_in_path=vcf_in_path, vcf_idx_path=vcf_idx_path, region=region)
    gather_worker = GatherWorker(manager=manager, vcf_in_path=vcf_in_path, vcf_idx_path=vcf_idx_path, vcf_out_path=vcf_out_path, region=region)
    workers = vtv_workers + [model_worker, scatter_worker, gather_worker]
    for worker in workers[::-1]:
        worker.start()
    #
    for worker in workers:
        worker.join()
