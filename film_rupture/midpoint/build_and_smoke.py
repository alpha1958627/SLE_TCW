"""Compile Basilisk shared kernel in the allocated node, then operator/step gate."""
import hashlib,json,os,subprocess,time
from pathlib import Path
import numpy as np
from midpoint_candidate import Film as ReferenceFilm

ROOT=Path(__file__).resolve().parent
if int(os.environ.get('SLURM_PROCID','0'))==0:
    out=ROOT/'compiled';out.mkdir(exist_ok=True)
    command=[os.environ['BASILISK']+'/qcc','-D_GNU_SOURCE','-D_MPI=1','-O3','-fno-fast-math','-ffp-contract=off','-shared','-fPIC',
        'basilisk_midpoint.c','-o','compiled/basilisk_midpoint.so','-lm']
    p=subprocess.run(command,cwd=str(ROOT),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    (out/'compile.log').write_text(p.stdout)
    if p.returncode:raise RuntimeError(p.stdout)
    manifest=dict(command=command,source_SHA256=hashlib.sha256((ROOT/'basilisk_midpoint.c').read_bytes()).hexdigest(),
        library_SHA256=hashlib.sha256((out/'basilisk_midpoint.so').read_bytes()).hexdigest())
    (out/'basilisk_build.json').write_text(json.dumps(manifest,indent=2))
from basilisk_midpoint_candidate import Film
film=Film(1.05,64,ROOT/'compiled/basilisk_midpoint.so');ref=ReferenceFilm(1.05,64)
rng=np.random.default_rng(2026102501);h=1+.01*rng.normal(size=(64,64));v=rng.normal(size=(64,64));dw=.03*rng.normal(size=(2,64,64))
errors={}
errors['force']=float(np.max(abs(film.force(h)-ref.force(h))))
errors['noise']=float(np.max(abs(film.stochastic(h,dw)-ref.stochastic(h,dw))))
df,dn=film.derivatives(h,v,dw);rf,rn=ref.derivatives(h,v,dw)
errors['dforce']=float(np.max(abs(df-rf)));errors['dnoise']=float(np.max(abs(dn-rn)))
assert max(errors.values())<1e-11,errors
new,d=film.step(h,.3,dw);expected,rd=ref.step(h,.3,dw)
errors['step']=float(np.max(abs(new-expected)))
assert errors['step']<1e-9 and d['residual']<=1e-10
assert abs(new.mean()-h.mean())<1e-12
if film.rank==0:
    result=dict(PASS=True,operators=errors,diagnostic=d,mpi_ranks=film.ranks,mean_error=float(new.mean()-h.mean()))
    (ROOT/'SMOKE_PASS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
