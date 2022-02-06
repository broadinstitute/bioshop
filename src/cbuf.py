import ctypes
import multiprocessing as mp

class CircularBuffer(object):
    def __init__(self, ctype=None, size=None):
        self.ctype = ctype
        self.size = size
        self.buffer = mp.RawArray(self.ctype, self.size)
        self.read_ptr = mp.RawValue(ctypes.c_uint, 0)
        self.write_ptr = mp.RawValue(ctypes.c_uint, 0)
        self.qsize = mp.RawValue(ctypes.c_uint, 0)
        self.buffer_cv = mp.Condition()

    def push(self, item):
        with self.buffer_cv:
            is_full = lambda: self.qsize.value == self.size
            while is_full():
                self.buffer_cv.wait(timeout=1)
            self.buffer[self.write_ptr.value] = item
            self.write_ptr.value = (self.write_ptr.value + 1) % self.size
            self.qsize.value += 1
            self.buffer_cv.notify_all()

    def pop(self, timeout=None):
        with self.buffer_cv:
            is_empty = lambda: self.qsize.value == 0
            while is_empty():
                if not self.buffer_cv.wait(timeout=timeout):
                    raise TimeoutError(timeout)
            src = self.buffer[self.read_ptr.value]
            # critical, don't forget to copy the data out!
            item = type(src)()
            ctypes.pointer(item)[0] = src
            self.read_ptr.value = (self.read_ptr.value + 1) % self.size
            self.qsize.value -= 1
            self.buffer_cv.notify_all()
        return item

