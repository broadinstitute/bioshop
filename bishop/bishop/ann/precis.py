# précis - a summary or abstract of a text or speech.

class CasualNamespace:
    _ns = None

    def __init__(self, *args, **ns):
        context = dict()
        if args:
            ns = args[0]
        self._init_attr('_ns', dict(ns))

    def copy(self):
        ns = self._ns.copy()
        return self.__class__(**ns)

    def _init_attr(self, key, val):
        super().__setattr__(key, val)

    def __getitem__(self, key):
        return self._ns[key]
    
    def __setitem__(self, key, val):
        self._ns[key] = val

    def __delitem__(self, key):
        del self._ns[key]

    def __contains__(self, key):
        return key in self._ns

    def __getattr__(self, key):
        if key in self:
            return self[key]
        raise AttributeError(key)

    def __delattr__(self, key):
        if key in self:
            del self[key]
        else:
            super().__delattr__(key)

    def __setattr__(self, key, val):
        self[key] = val

    def __iter__(self):
        return iter(self._ns.items())
    
    def as_dict(self):
        return dict(iter(self))

    def __repr__(self):
        return f'{self.__class__.__name__}({repr(self._ns)})'

    def keys(self):
        return self._ns.keys()

    def items(self):
        return self._ns.items()

    def values(self):
        return self._ns.values()
    
    def update(self, ns):
        self._ns.update(ns)

    def __getstate__(self):
        return (self._ns.copy(),)

    def __setstate__(self, state):
        self._init_attr('_ns', dict(state[0]))

class Domain(CasualNamespace):
    def __init__(self, *args, domain_name=None, **kw):
        super().__init__(*args, **kw)
        self._init_attr('domain_name', domain_name)

    def copy(self):
        ns = self._ns.copy()
        return self.__class__(domain_name=self.domain_name, **ns)

    def __getstate__(self):
        state = super().__getstate__()
        return (self.domain_name, state)

    def __setstate__(self, state):
        self._init_attr('domain_name', state[0])
        super().__setstate__(state[-1])

class FilterDomain(CasualNamespace):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        if 'count' not in self:
            self.count = 0
        if 'log' not in self:
            self.log = tuple()

    def __bool__(self):
        return self.count > 0

    def set_filter(self, log=None):
        self.count += 1
        log = log or 'no log provided'
        self.log = tuple(list(self.log) + [log])
    
class Precis(CasualNamespace):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        if 'filter' not in self:
            self['filter'] = FilterDomain()

    def copy(self, include_domains=None, exclude_domains=None):
        if (include_domains and exclude_domains):
            raise TypeError('include / exclude mutually exclusive')
        ns = {}
        for (name, domain) in self._ns.items():
            if include_domains and name not in include_domains:
                continue
            if exclude_domains and name in exclude_domains:
                continue
            ns[name] = domain.copy()
        return self.__class__(**ns)

    def flatten(self, include_domains=None, exclude_domains=None):
        if (include_domains and exclude_domains):
            raise TypeError('include / exclude mutually exclusive')
        flat_ns = {}
        for (name, domain) in self._ns.items():
            if include_domains and name not in include_domains:
                continue
            if exclude_domains and name in exclude_domains:
                continue
            flat_ns.update(domain._ns)
        return self.__class__(**flat_ns)

    def __getitem__(self, key):
        return self._ns[key].as_dict()
    
    def __getattr__(self, key):
        if key not in self:
            self[key] = Domain(domain_name=key)
        return self._ns[key]

    def __setattr__(self, key, val):
        raise NotImplementedError()

    def __iter__(self):
        for (key, val) in self._ns.items():
            val = val.as_dict()
            yield (key, val)

def run_tests():
    cn = CasualNamespace()
    cn.one = '1'
    cn.two = 2
    assert(cn.one == '1')
    assert(cn.two == 2)
    assert dict(one='1', two=2) == dict(cn)
    assert 'one' in cn
    assert cn['one'] == '1'
    cn['three'] = 3
    assert cn.three == 3
    del cn['three']
    assert 'three' not in cn
    cn.three = 3
    assert 'three' in cn
    del cn.three
    assert 'three' not in cn
    dm = Domain(domain_name='bob')
    assert dm.domain_name == 'bob'
    pr = Precis()
    pr.abc.xzy = '123456'
    assert pr.abc.xzy == '123456' 
    to_match = {'filter': {'count': 0, 'log': ()}, 'abc': {'xzy': '123456'}}
    assert dict(pr) == to_match
    pr2 = pr.copy()
    assert dict(pr) == dict(pr2)
    pr2.filter.set_filter('rude')
    to_match = {'filter': {'count': 1, 'log': ('rude', )}, 'abc': {'xzy': '123456'}}
    assert dict(pr2) == to_match
    assert dict(pr) != dict(pr2)
    import pickle
    pr_pickle = pickle.dumps(pr)
    pr3 = pickle.loads(pr_pickle)
    assert dict(pr) == dict(pr3)

if __name__ == '__main__':
    run_tests()
