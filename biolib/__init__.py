import os
import math
#import progressbar
import operator
import itertools

import platform

# given an iterable of pairs return the key corresponding to the greatest value
def argmax(pairs,f=max):
    return f(pairs, key=operator.itemgetter(1))[0]

def argmin(pairs):
    return argmax(pairs,f=min)

# given an iterable of values return the index of the greatest value
def argmax_index(values,f=max):
    return argmax(enumerate(values),f)
def argmin_index(values):
    return argmax_index(values,f=min)

# given an iterable of keys and a function f, return the key with largest f(key)
def argmax_f(keys, f):
    return max(keys, key=f)

def argmin_f(keys, f):
    return min(keys, key=f)

try:
    if platform.python_implementation() == "PyPy":
        import numpypy as np
        pypy = True
    else:
        import numpy as np
        pypy = False
except ImportError as err:
        print err, "(problem importing NumPy or NumPyPy)"

def log2_type(val):
    return float(np.log2(float(val)))

def pseudocount_log2_type(val):
    return float(np.log2(float(val)+1))

class Cseq(object):
    def __init__(self,name,pls,mns,header='',filename=None,default=0,refs=None):
        self.name = name
        if isinstance(pls,dict):
            self._pls = pls
        else:
            self._pls = { self.name:pls }
        if isinstance(mns,dict):
            self._mns = mns
        else:
            self._mns = { self.name:mns }
        if filename is not None:
            print 'warning: filename property is deprecated and is now name'
        self.default = default
        self.cur = self._pls.keys()[0]
        self.header = header
    @property
    def pls(self):
        return self._pls[self.cur]
    @property
    def mns(self):
        return self._mns[self.cur]
    @property
    def filename(self):
        return self.name
    @property
    def refs(self):
        return self._pls.keys()
    def ref(self,refname):
        if refname not in self._pls:
            raise ValueError, 'Invalid refname %s' % refname
        self.cur = refname
        return self
    @classmethod
    def from_file(cls,f,vtype=int,*args,**kwargs):
        header = ''
        last_pos = 0
        for line in f:
            if line[:2] == '##':
                header += line
                last_pos += len(line)
            else:
                break
        f.seek(last_pos)
        name = os.path.basename(f.name)
        cur = f.next().lstrip('#').strip()
        pls = {cur:[]}
        mns = {cur:[]}
        for line in f:
            if line[0] == '#':
                cur = line.lstrip('#').strip()
                pls[cur] = []
                mns[cur] = []
                continue
            plsv,mnsv = line.strip().split()
            pls[cur].append(vtype(plsv))
            mns[cur].append(vtype(mnsv))
        return cls(name,pls,mns,header,*args,**kwargs)
    def write(self,file,fmt='%d %d'):
        if not hasattr(file,'write'):
            file = open(file,'w')
        file.write(self.header)
        if len(self.header) and not self.header.endswith('\n'):
            file.write('\n')
        for ref in self.refs:
            file.write('#%s\n' % ref)
            for p,m in iter(self.ref(ref)):
                file.write(fmt%(p,m) + '\n')
        file.close()
    def get_strand(self,hint):
        if hint == '+' or hint == 1 or hint == 'plus' or hint == 'pls':
            return self.pls
        elif hint == '-' or hint == -1 or hint == 'minus' or hint == 'mns':
            return self.mns
        else:
            raise ValueError, 'Unknown strand key %s' % str(hint)
    def apply(self,func):
        for ref in self.refs:
            for i in range(len(self.ref(ref))):
                self.pls[i] = func(self.pls[i])
                self.mns[i] = func(self.mns[i])
    def __iter__(self):
        return itertools.izip_longest(self.pls,self.mns,fillvalue=self.default)
    def __getitem__(self,k):
        try:
            p = self.pls.__getitem__(k)
        except IndexError: 
            p = self.default
        try:
            m = self.mns.__getitem__(k)
        except IndexError:
            m = self.default
        return p,m
    def __len__(self):
        return max(len(self.pls),len(self.mns))

class CseqFileIter(object):
    def __init__(self,file,vtype=int):
        if not hasattr(file,'read'):
            self.file = open(file)
        else:
            self.file = file
        #self.name = self.file.next().lstrip('#').strip()
        self.name = self.file.name
        self.vtype = vtype
    def __iter__(self):
        return self
    def next(self):
        return [self.vtype(v) for v in self.file.next().split()]
     
class Eseq(object):
    def __init__(self,name,values,filename=None):
        self.name = name
        self.values = values
        self.filename = filename
    @classmethod
    def from_file(cls,f,vtype=int):
        name = f.next().lstrip('>').strip()
        values = []
        for line in f:
            for x in line.strip().split():
                values.append(vtype(x))
        return cls(name,values,filename=f.name)
    def log_transform(self,base=2):
        #old = np.seterr(all='ignore') # avoid log2(0) error warnings
        self.values = [float(np.log2(v)) for v in self.values]
        #np.seterr(**old)
    def __iter__(self):
        return iter(self.values)
    def __getitem__(self,k):
        return self.values.__getitem__(k)
    def __len__(self):
        return len(self.values)
        
def parse_cseq(filename,vtype=int):
    return Cseq.from_file(open(filename),vtype)
def parse_float_cseq(filename):
    return Cseq.from_file(open(filename),float)
def parse_cseq_log2_transform(filename):
    if not pypy:
        old = np.seterr(all='ignore') # avoid log2(0) error warnings
    eseq = Cseq.from_file(open(filename),log2_type)
    if not pypy:
        np.seterr(**old)
    return eseq
def parse_cseq_pseudocount_log2_transform(filename):
    if not pypy:
        old = np.seterr(all='ignore') # avoid log2(0) error warnings
    eseq = Cseq.from_file(open(filename),pseudocount_log2_type)
    if not pypy:
        np.seterr(**old)
    return eseq

def cseq_file_iter(filename,vtype=int):
    return CseqFileIter(filename,vtype)

def parse_eseq(filename,vtype=int):
    return Eseq.from_file(open(filename),vtype)
def parse_float_eseq(filename):
    return Eseq.from_file(open(filename),float)
def parse_eseq_log2_transform(filename):
    if not pypy:
        old = np.seterr(all='ignore') # avoid log2(0) error warnings
    eseq = Eseq.from_file(open(filename),log2_type)
    if not pypy:
        np.seterr(**old)
    return eseq

def write_eseq(stream,iterable,name='eseq',width=60,fmt='%d '):
    cur_line = ''
    stream.write('>%s\n' % name)
    for n in iterable:
        string = fmt % n
        if len(cur_line) + len(string) > width:
            stream.write('%s\n' % cur_line)
            cur_line = ''
        cur_line += string
def write_cseq(stream,pls,mns,name='eseq',fmt='%d %d'):
    Cseq(name,pls,mns).write(stream,fmt=fmt)

'''
def standard_progressbar(widgets=[progressbar.Percentage(),
                                      progressbar.Bar(),
                                      progressbar.ETA()], **kwargs):
    return progressbar.ProgressBar(widgets=widgets,**kwargs)
'''
