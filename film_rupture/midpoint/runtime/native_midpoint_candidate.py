"""C FV kernels under the unchanged Python Newton/GMRES/midpoint algorithm.

Prototype only. A Film instance is used by one trajectory. Cache reuse is enabled
only inside step(), where the reference creates immutable midpoint arrays.
"""
import ctypes
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import numpy as np
from midpoint_candidate import Film as ReferenceFilm

ROOT=Path(__file__).resolve().parent


def build(folder):
    folder=Path(folder);folder.mkdir(exist_ok=False)
    source=ROOT.parent/'sle_smid_fv_c99.c'
    frozen=folder/source.name;frozen.write_bytes(source.read_bytes())
    library=folder/('kernel.dll' if sys.platform=='win32' else 'kernel.so')
    command=['gcc','-std=c99','-O3','-fno-fast-math','-ffp-contract=off','-shared','-static-libgcc']
    if sys.platform!='win32':command+=['-fPIC']
    command+=[str(frozen),'-o',str(library),'-lm']
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    (folder/'compile.log').write_text(result.stdout,encoding='utf-8')
    if result.returncode:raise RuntimeError(result.stdout)
    manifest=dict(schema=1,command=command,compiler=shutil.which('gcc'),
        compiler_version=subprocess.check_output(['gcc','--version'],universal_newlines=True),
        system=platform.system(),machine=platform.machine(),library=library.name,
        library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(frozen.read_bytes()).hexdigest())
    (folder/'build.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return library


class Film(ReferenceFilm):
    def __init__(self,*args,library,**kwargs):
        super().__init__(*args,**kwargs)
        self.library=ctypes.CDLL(str(Path(library).resolve()))
        integer=ctypes.c_int;real=ctypes.c_double;pointer=ctypes.c_void_p
        self.library.film_prepare.argtypes=[integer,real,real,real,pointer,pointer]
        self.library.film_force.argtypes=[integer,real,pointer,pointer]
        self.library.film_noise.argtypes=[integer,real,real,pointer,pointer,pointer]
        self.library.film_derivative.argtypes=[integer,real,real,real]+[pointer]*6
        for name in ('film_prepare','film_force','film_noise','film_derivative'):
            getattr(self.library,name).restype=None
        self.buffer=np.empty((8,self.n,self.n));self.work=np.empty((5,self.n,self.n))
        self._cached_height=None;self._inside_step=False

    def _prepare(self,h):
        if not self._inside_step or self._cached_height is not h:
            field=np.ascontiguousarray(h,dtype=np.float64)
            if field.shape!=(self.n,self.n) or not np.isfinite(field).all() or field.min()<=0:
                raise ValueError('Nonpositive/nonfinite native midpoint; no repair')
            self.library.film_prepare(self.n,self.dx,self.a,self.slip,field.ctypes.data,self.buffer.ctypes.data)
            self._cached_height=h

    def force(self,h):
        if self.pressure!='energy':return super().force(h)
        self._prepare(h);out=np.empty((self.n,self.n),dtype=np.float64)
        self.library.film_force(self.n,self.dx,self.buffer.ctypes.data,out.ctypes.data)
        return out

    def stochastic(self,h,dw):
        if self.pressure!='energy':return super().stochastic(h,dw)
        self._prepare(h);noise=np.ascontiguousarray(dw,dtype=np.float64)
        assert noise.shape==(2,self.n,self.n)
        out=np.empty((self.n,self.n))
        self.library.film_noise(self.n,self.dx,self.phi,self.buffer.ctypes.data,noise.ctypes.data,out.ctypes.data)
        return out

    def derivatives(self,h,v,dw):
        if self.pressure!='energy':return super().derivatives(h,v,dw)
        self._prepare(h);direction=np.ascontiguousarray(v,dtype=np.float64)
        noise=np.ascontiguousarray(dw,dtype=np.float64)
        assert direction.shape==(self.n,self.n) and noise.shape==(2,self.n,self.n)
        df=np.empty((self.n,self.n));dn=np.empty_like(df)
        self.library.film_derivative(self.n,self.dx,self.a,self.phi,self.buffer.ctypes.data,
            direction.ctypes.data,noise.ctypes.data,self.work.ctypes.data,df.ctypes.data,dn.ctypes.data)
        return df,dn

    def step(self,*args,**kwargs):
        self._inside_step=True;self._cached_height=None
        try:return super().step(*args,**kwargs)
        finally:self._inside_step=False;self._cached_height=None
