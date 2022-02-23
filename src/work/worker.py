import os
import time
from tempfile import TemporaryDirectory
import multiprocessing as mp

class Worker(mp.Process):
    def __init__(self, manager=None, **kw):
        super().__init__(**kw)
        #assert manager
        self.manager = manager
        self._running = mp.Event()

    @property
    def running(self):
        return self._running.is_set()
    
    def shutdown(self):
        self._running.clear()

    def _run(self):
        pass
    
    def wait_until_running(self, timeout=None):
        return self._running.wait(timeout=timeout)

    def run(self):
        pwd = os.getcwd()
        with TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            msg = f"Starting {self.__class__.__name__} worker run={self.running}, tmpdir={tmpdir}"
            print(msg)
            self._running.set()
            self._run()
        os.chdir(pwd)
        msg = f"{self.__class__.__name__} shutting down"
        print(msg)
        self.shutdown()
