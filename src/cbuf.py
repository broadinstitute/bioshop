import ctypes
import multiprocessing as mp

class CircularBuffer(object):
    def __init__(self, ctype=None, size=None):
        self.ctype = ctype
        self.size = size
        self.buffer = mp.RawArray(self.ctype, self.size)
        self.read_ptr = mp.RawValue(ctypes.c_uint, 0)
        self.write_ptr = mp.RawValue(ctypes.c_uint, 0)
        self.buffer_cv = mp.Condition()

    def qsize(self):
        with self.buffer_cv:
            r_val = self.read_ptr.value
            w_val = self.write_ptr.value
        if w_val == r_val:
            return 0
        if w_val > r_val:
            return w_val - r_val
        return (self.size - r_val) + w_val

    def push(self, item):
        is_full = lambda: ((self.write_ptr.value + 1) % self.size) == self.read_ptr.value
        with self.buffer_cv:
            while is_full():
                self.buffer_cv.wait(timeout=1)
            self.buffer[self.write_ptr.value] = item
            self.write_ptr.value = (self.write_ptr.value + 1) % self.size
            self.buffer_cv.notify_all()

    def pop(self):
        is_empty = lambda: self.read_ptr.value == self.write_ptr.value
        with self.buffer_cv:
            while is_empty():
                self.buffer_cv.wait(timeout=1)
            item = self.buffer[self.read_ptr.value]
            self.read_ptr.value = (self.read_ptr.value + 1) % self.size
            self.buffer_cv.notify_all()
        return item
