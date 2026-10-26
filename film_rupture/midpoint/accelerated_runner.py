"""Full fingerprinted native/reference trajectories with explicit Brownian root blocks.

Inherits the verified same-path adaptive step. Reference campaigns are untouched.
Schema 2 is deliberately distinct from the old schema-1 reference checkpoints.
"""
import copy
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import numpy as np
from adaptive_runner import Trajectory as ReferenceTrajectory,BrownianTape,runtime
from native_midpoint_candidate import Film as NativeFilm

ROOT=Path(__file__).resolve().parent


def fingerprint(kernel,library=None):
    names=['accelerated_runner.py','adaptive_runner.py','midpoint_candidate.py']
    details=None
    if kernel=='native':
        names+=['native_midpoint_candidate.py','native_midpoint_kernel.c']
        library=Path(library).resolve();manifest=library.parent/'build.json'
        details=json.loads(manifest.read_text())
        assert details['schema']==1 and details['library']==library.name
        assert details['library_sha256']==hashlib.sha256(library.read_bytes()).hexdigest(),'Native binary changed'
        assert details['source_sha256']==hashlib.sha256((ROOT/'native_midpoint_kernel.c').read_bytes()).hexdigest(),'Native C source/build mismatch'
        details=dict(build_manifest=details,build_manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest())
    return dict(kernel=kernel,sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
        runtime=dict(runtime(),system=platform.system(),machine=platform.machine(),byteorder=sys.byteorder),native=details)


class Trajectory(ReferenceTrajectory):
    def __init__(self,height,n,seed,dt_cap=.03,h=None,tiered=True,event_resolution=1e-5,
                 base_dt=.001,depth=12,kernel='native',library=None):
        assert kernel in ('native','reference') and math.isfinite(base_dt) and base_dt>0
        assert kernel=='native' or library is None
        self.kernel=kernel;self.library_path=str(Path(library).resolve()) if library is not None else None
        self._fingerprint=fingerprint(kernel,self.library_path)
        super().__init__(height,n,seed,dt_cap=dt_cap,h=h,tiered=tiered,event_resolution=event_resolution,depth=depth)
        self.noise=BrownianTape(seed,(2,n,n),base_dt=base_dt,depth=depth)
        self.cap_ticks=self.noise.ticks(dt_cap)
        assert self.cap_ticks>0 and abs(self.cap_ticks*self.noise.quantum-dt_cap)<1e-12
        self.proposal=self.cap_ticks
        self.event_resolution=max(1,self.noise.ticks(event_resolution))
        if kernel=='native':self.film=NativeFilm(height,n,library=self.library_path)

    def metadata(self):
        assert self._fingerprint==fingerprint(self.kernel,self.library_path),'Trajectory source/runtime/native build changed'
        result=super().metadata()
        result.update(schema=2,kernel=self.kernel,library_path=self.library_path,fingerprint=self._fingerprint,
            sources=self._fingerprint['sources'],runtime=self._fingerprint['runtime'])
        return result

    def fork(self,seed):
        self.metadata()  # Refuse to fork a live object after changing its code/library.
        child=type(self)(self.film.height,self.film.n,seed,dt_cap=self.cap_ticks*self.noise.quantum,
            h=self.h if self.h.min()>.01 else np.ones_like(self.h),tiered=self.tiered,
            event_resolution=self.event_resolution*self.noise.quantum,base_dt=self.noise.base_dt,depth=self.noise.depth,
            kernel=self.kernel,library=self.library_path)
        child.h=self.h.copy();child.noise=copy.deepcopy(self.noise);child.noise.seed=int(seed)
        for name in ('cap_ticks','proposal','streak','step_count','ticks','event_resolution','status','event',
                     'initial_mean','ledger','last_interval'):
            setattr(child,name,copy.deepcopy(getattr(self,name)))
        return child

    @classmethod
    def load(cls,path,library=None):
        with np.load(path,allow_pickle=False) as z:
            h=z['h'].copy();values=z['noise_values'].copy();payload=str(z['metadata'].item());expected=str(z['sha256'].item())
        assert hashlib.sha256(payload.encode()+h.tobytes()+values.tobytes()).hexdigest()==expected,'Checkpoint SHA mismatch'
        data=json.loads(payload)
        assert data['schema']==2,'Use the frozen original reader for old reference checkpoints'
        actual_library=library if library is not None else data['library_path']
        current=fingerprint(data['kernel'],actual_library)
        assert data['fingerprint']==current,'Checkpoint source/runtime/native fingerprint mismatch'
        assert data['sources']==current['sources'] and data['runtime']==current['runtime'],'Inconsistent checkpoint provenance'
        obj=cls(data['height'],data['n'],data['seed'],dt_cap=data['cap_ticks']*data['base_dt']/(1<<data['depth']),
            h=h if h.min()>.01 else np.ones_like(h),tiered=data['tiered'],
            event_resolution=data['event_resolution']*data['base_dt']/(1<<data['depth']),
            base_dt=data['base_dt'],depth=data['depth'],kernel=data['kernel'],library=actual_library)
        assert obj.cap_ticks==data['cap_ticks'] and obj.noise.base_dt==data['base_dt'] and obj.noise.depth==data['depth']
        assert obj.film.a==data['A'] and obj.film.phi==data['phi'] and obj.film.slip==data['slip'] and obj.film.pressure==data['pressure']
        obj.h=h
        assert values.shape==(len(data['noise_keys']),2,data['n'],data['n'])
        for name in ('ticks','proposal','streak','step_count','event_resolution','status','event','initial_mean','ledger','last_interval'):
            setattr(obj,name,data[name])
        obj.noise.cursor=data['noise_cursor']
        obj.noise.cache={tuple(key):value.copy() for key,value in zip(data['noise_keys'],values)}
        assert obj.ticks==obj.noise.cursor
        return obj
