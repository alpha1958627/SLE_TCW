"""All-height fixed-time TCW gate; five seeds per height, paired dt caps."""
import argparse,json,time
from pathlib import Path
import numpy as np
from basilisk_midpoint_candidate import Trajectory

ROOT=Path(__file__).resolve().parent
HEIGHTS=[.6,.7,.8,.9,1.,1.1,1.2,1.3]
parser=argparse.ArgumentParser();parser.add_argument('task',type=int);args=parser.parse_args()
height=HEIGHTS[args.task//5];seed=2026102201+args.task%5
rows=[];fields={}
for cap in [.03,.3]:
    s=Trajectory(height,64,seed,dt_cap=cap)
    start=time.perf_counter();end=s.noise.ticks(300.)
    while s.ticks<end and s.status=='RUNNING':s.advance(stop_ticks=end)
    assert s.status=='RUNNING' and s.time==300. and s.h.min()>.01
    power=abs(np.fft.fft2(s.h-s.h.mean())/64**2)**2
    film=s.film;omega=film.omega;d=film.phi*(1+film.slip)*film.lam/film.length**2
    expected=np.zeros_like(omega);mask=abs(omega)>1e-18
    expected[mask]=-d[mask]/omega[mask]*np.expm1(-2*omega[mask]*300.)
    expected[~mask]=2*d[~mask]*300.
    k=np.fft.fftfreq(64)*64;r=np.hypot(k[:,None],k[None,:]);high=(r>24)&(r<=46)
    row=dict(height=height,seed=seed,cap=cap,time=300.,wall=time.perf_counter()-start,
        steps=s.step_count,high_power=float(power[high].sum()),FD_ratio=float(power[high].sum()/expected[high].sum()),
        hmin=float(s.h.min()),mean=float(s.h.mean()),max_residual=s.ledger['max_residual'],ledger=s.ledger)
    rows.append(row);fields[str(cap)]=s.h.copy()
    if s.film.rank==0:
        folder=ROOT/'tcw_probe'/('task_%02d'%args.task);folder.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(folder/('cap%g.npz'%cap),h=s.h,power=power,expected=expected,time=300.)
        print('TCW',json.dumps(row),flush=True)
if s.film.rank==0:
    (folder/'summary.json').write_text(json.dumps(rows,indent=2))
