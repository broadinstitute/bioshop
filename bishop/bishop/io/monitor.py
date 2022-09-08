import time
import threading as th
import multiprocessing as mp

from humanize import metric

def elapsed_time_banner(elapsed):
    hours = int(elapsed / (60 * 60))
    if hours:
        elapsed -= hours * (60 * 60)
    minutes = int(elapsed / 60)
    if minutes:
        elapsed -= minutes * 60
    seconds = int(elapsed)
    if hours:
        return f'{hours:02d}h{minutes:02d}m'
    else:
        return f'{minutes:02d}m{seconds:02d}s'

def timestamp():
    return time.time()

class Singleton(object):
    _instance = None
    
    def __new__(cls, *args, **kw):
        if cls._instance is None:
            cls._instance = super().__new__(cls, *args, **kw)
            cls._instance._init(*args, **kw)
        return cls._instance
    
    def _init(self, *args, **kw):
        pass

class Monitor(Singleton):
    def _init(self):
        (self._pipe_recv, self._pipe_send) = mp.Pipe()
        self.running = False
        self.thread = None
        self.tput = Throughput()
        self.seconds_per_report = 0
        self.last_report = None
        self.enter_ts = self.exit_ts = None

    def enable_reporting(self, seconds_per_report=30):
        self.seconds_per_report = seconds_per_report

    def disable_reporting(self):
        self.seconds_per_report = 0

    def print_report(self, domain=None):
        rpt = self.tput.get_report(domain=domain)
        #if not rpt:
            #return
        if self.enter_ts is not None:
            elapsed = time.time() - self.enter_ts
            elapsed = elapsed_time_banner(elapsed)
            ts = f'[{elapsed}] '
            tslen = len(ts)
            spc = (' ' * tslen)
            rpt = rpt.split('\n')
            rpt = [(ts if idx == 0 else spc) + ln for (idx, ln) in enumerate(rpt)]
            rpt = str.join('\n', rpt)
        print(rpt)

    def do_report(self, domain=None):
        ts = timestamp()
        if self.seconds_per_report <= 0:
            return
        if self.last_report is None:
            self.last_report = ts
            return
        elapsed = (ts - self.last_report)
        if elapsed < self.seconds_per_report:
            return
        self.print_report(domain=domain)
        self.last_report = ts

    def report_runtime(self):
        elapsed = self.exit_ts - self.enter_ts 
        rt_banner = elapsed_time_banner(elapsed)
        rpt = f'Runtime [{rt_banner}]'
        print(rpt)

    def get_remote(self, domain=None):
        return TelemetryRemote(cnx=self._pipe_send, domain=domain)
    
    def main_loop(self, timeout=.5, domain=None):
        self.running = True
        while self.running:
            self.do_report(domain=domain)
            if self._pipe_recv.poll(timeout):
                item = self._pipe_recv.recv()
                op = item.pop('op')
                if op == 'update':
                    self.tput.update(**item)
                elif op == 'add':
                    self.tput.add(**item)
                elif op == 'noop':
                    continue
                else:
                    raise ValueError(op)

    def stop(self):
        if self.thread:
            self.running = False
            pckt = dict(op='noop')
            self._pipe_send.send(pckt)
            self.thread.join()
            self.thread = None
            self.exit_ts = time.time()
            self.report_runtime()

    def start(self, domain=None):
        kw = {'domain': domain}
        assert self.running is False
        self.enter_ts = time.time()
        self.thread = th.Thread(target=self.main_loop, kwargs=kw, daemon=True)
        self.thread.start()
        return self

    def __enter__(self):
        if not self.running:
            self.start()

    def __exit__(self, *args):
        self.stop()

def get_remote_monitor(*args, **kw):
    mon = Monitor()
    return mon.get_remote(*args, **kw)

class Counter(object):
    def __init__(self, name=None, value=0, _type=None):
        self.name = name
        self.type = _type or type(value)
        self.value = self.type(value)
        self.start_ts = time.time()

    @property
    def rate(self):
        delta = (timestamp() - self.start_ts)
        return self.value / delta

    def __str__(self):
        if self.value > 9999:
            value = metric(self.value)
        else:
            value = str(self.value)
        if self.rate > 9999:
            rate = metric(self.rate)
        else:
            rate = f'{self.rate:.02f}'
        return f'{self.name}: {value} ({rate}/s)'

    def __add__(self, other: float):
        self.value += self.type(other)
        return self

    def update(self, other: float):
        self.value = self.type(other)
        return self

    def __getstate__(self):
        return (self.start_ts, self.name, self.value)
    
    def __setstate__(self, state):
        (self.start_ts, self.name, self.value) = state
        self.type = type(self.value)

class PerformanceLog(Singleton):
    def _init(self, *args, **kw):
        self.tput = Throughput
        self.pipe

class TelemetryRemote(object):
    def __init__(self, cnx=None, rps=1, domain=None):
        self.cnx = cnx
        self.rps = rps
        self.domain = domain
        self.last_report = 0
        self._cache = {}

    def flush_report(self):
        for (name, value) in self._cache.items():
            packet = dict(name=name, domain=self.domain, value=value, op='add')
            self.cnx.send(packet)
        self._cache = {}
        self.last_report = timestamp()

    def update(self, name=None, value=None):
        packet = dict(name=name, domain=self.domain, value=value, op='update')
        self.cnx.send(packet)

    def add(self, name=None, value=None):
        if name not in self._cache:
            self._cache[name] = 0
        self._cache[name] += value
        if (timestamp() + self.last_report) > self.rps:
            self.flush_report()

class Throughput(object):
    def __init__(self):
        self.domains = {}

    def get_domain_counters(self, domain=None):
        domain = domain or 'default'
        if domain not in self.domains:
            self.domains[domain] = {}
        return self.domains[domain]

    def add(self, name=None, domain=None, value=None):
        counters = self.get_domain_counters(domain=domain)
        if name not in counters:
            counters[name] = Counter(name=name, value=value)
        else:
            counters[name] += value

    def update(self, name=None, domain=None, value=None):
        counters = self.get_domain_counters(domain=domain)
        if name not in counters:
            counters[name] = Counter(name=name, value=value)
        else:
            counters[name].update(value)

    def get_report(self, domain=None):
        if domain is None:
            domains = self.domains.keys()
        else:
            domains = [domain]
        rpt = []
        for domain in domains:
            counters = self.get_domain_counters(domain=domain)
            line = str.join(', ', map(str, counters.values()))
            line = f'({domain}): {line}'
            rpt.append(line)
        return str.join('\n', rpt)

def test():
    import random
    mon = Monitor()
    rt = mon.get_remote()
    names = ['bob', 'harry', 'jill', 'steve']
    mon.enable_reporting(5)
    with mon:
        for x in range(100):
            slp = random.random() * .5
            val = int(random.random() * 100)
            #slp = 1
            #val = 100
            time.sleep(slp)
            nm = random.choice(names)
            rt.add(nm, val)
    #tr.stop()

if __name__ == '__main__':
    test()
