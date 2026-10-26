"""Thin Python boundary for the qcc/Basilisk C midpoint solve."""
import ctypes
import hashlib
import json
from pathlib import Path
import numpy as np
from midpoint_candidate import Film as ReferenceFilm
from adaptive_runner import Trajectory as ReferenceTrajectory, BrownianTape, runtime

ROOT=Path(__file__).resolve().parent


class Film(ReferenceFilm):
    def __init__(self,height,n,library):
        super().__init__(height,n)
        self.library=ctypes.CDLL(str(Path(library).resolve()),mode=ctypes.RTLD_GLOBAL)
        p=ctypes.c_void_p;d=ctypes.c_double;i=ctypes.c_int
        self.library.bm_init.argtypes=[i,d];self.library.bm_init.restype=i
        self.library.bm_eval.argtypes=[p,p,p,d,i,p];self.library.bm_eval.restype=i
        self.library.bm_step.argtypes=[p,d,p,d,i,p,p];self.library.bm_step.restype=i
        assert self.library.bm_init(n,height)==0,'Basilisk MPI grid coverage mismatch'
        self.rank=self.library.bm_rank();self.ranks=self.library.bm_ranks()

    def evaluate(self,h,dw,v,step,kind):
        field=np.ascontiguousarray(h,dtype=np.float64)
        noise=np.ascontiguousarray(dw,dtype=np.float64)
        direction=np.ascontiguousarray(v,dtype=np.float64)
        out=np.empty_like(field)
        rc=self.library.bm_eval(field.ctypes.data,noise.ctypes.data,direction.ctypes.data,step,kind,out.ctypes.data)
        if rc:raise ValueError('Basilisk nonpositive/nonfinite midpoint')
        return out

    def force(self,h):return self.evaluate(h,np.zeros((2,self.n,self.n)),h,1.,0)
    def stochastic(self,h,dw):return self.evaluate(h,dw,h,1.,1)
    def derivatives(self,h,v,dw):
        zero=np.zeros_like(dw)
        df=self.evaluate(h,zero,v,1.,3)
        return df,self.evaluate(h,dw,v,0.,3)

    def step(self,old,step,dw,tolerance=1e-10,maxiter=15):
        initial=np.ascontiguousarray(old,dtype=np.float64)
        noise=np.ascontiguousarray(dw,dtype=np.float64)
        result=np.empty_like(initial);diagnostic=np.zeros(3)
        rc=self.library.bm_step(initial.ctypes.data,step,noise.ctypes.data,tolerance,maxiter,result.ctypes.data,diagnostic.ctypes.data)
        if rc:raise RuntimeError('Basilisk midpoint Newton/GMRES failure code=%d'%rc)
        return result,dict(newton=int(diagnostic[0]),residual=float(diagnostic[1]),evaluations=int(diagnostic[2]))


def fingerprint(library):
    names=['basilisk_midpoint.c','basilisk_midpoint_candidate.py','adaptive_runner.py','midpoint_candidate.py']
    return dict(sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
        library_SHA256=hashlib.sha256(Path(library).read_bytes()).hexdigest(),runtime=runtime())


class Trajectory(ReferenceTrajectory):
    def __init__(self,height,n,seed,dt_cap=.3,h=None,tiered=True,event_resolution=1e-5,
                 base_dt=.064,depth=18,library=None):
        self.library_path=str(Path(library or ROOT/'compiled/basilisk_midpoint.so').resolve())
        super().__init__(height,n,seed,dt_cap=dt_cap,h=h,tiered=tiered,event_resolution=event_resolution,depth=depth)
        self.film=Film(height,n,self.library_path)
        self.noise=BrownianTape(seed,(2,n,n),base_dt=base_dt,depth=depth)
        self.cap_ticks=self.noise.ticks(dt_cap);self.proposal=self.cap_ticks
        self.event_resolution=max(1,self.noise.ticks(event_resolution))
        assert abs(self.cap_ticks*self.noise.quantum-dt_cap)<1e-12
        self._fingerprint=fingerprint(self.library_path)

    def metadata(self):
        assert fingerprint(self.library_path)==self._fingerprint,'Basilisk source/library/runtime changed'
        result=super().metadata()
        result.update(schema=3,kernel='basilisk_midpoint',library_path=self.library_path,
            fingerprint=self._fingerprint,sources=self._fingerprint['sources'],runtime=runtime(),mpi_ranks=self.film.ranks)
        return result

    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as z:
            h=z['h'].copy();values=z['noise_values'].copy();payload=str(z['metadata'].item());digest=str(z['sha256'].item())
        assert hashlib.sha256(payload.encode()+h.tobytes()+values.tobytes()).hexdigest()==digest
        d=json.loads(payload);assert d['schema']==3
        assert fingerprint(d['library_path'])==d['fingerprint'],'Frozen checkpoint fingerprint mismatch'
        obj=cls(d['height'],d['n'],d['seed'],dt_cap=d['cap_ticks']*d['base_dt']/(1<<d['depth']),
            h=h if h.min()>.01 else np.ones_like(h),tiered=d['tiered'],
            event_resolution=d['event_resolution']*d['base_dt']/(1<<d['depth']),base_dt=d['base_dt'],depth=d['depth'],library=d['library_path'])
        assert obj.film.ranks==d['mpi_ranks'],'Restart MPI configuration changed without validation'
        obj.h=h
        for name in ('ticks','proposal','streak','step_count','event_resolution','status','event','initial_mean','ledger','last_interval'):
            setattr(obj,name,d[name])
        obj.noise.cursor=d['noise_cursor']
        obj.noise.cache={tuple(k):v.copy() for k,v in zip(d['noise_keys'],values)}
        return obj
