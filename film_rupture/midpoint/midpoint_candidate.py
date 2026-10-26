"""Conservative FV stochastic midpoint prototype; not a validated rupture solver.

Native square grid, arithmetic face heights, no clipping, no noise resampling.
The face noise's Ito-to-Stratonovich drift vanishes because one face transfers
height between its two adjacent cells without changing their face average.
Nonlinear invariant-measure and first-passage errors still require testing.
"""
import argparse
import csv
import hashlib
import inspect
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
GMRES_RELATIVE_KEY = 'rtol' if 'rtol' in inspect.signature(gmres).parameters else 'tol'


class Film:
    def __init__(self, height, n, unit_box=False, pressure='energy'):
        self.height, self.n = height, n
        self.unit_box=unit_box
        self.pressure=pressure
        self.mean_height = height if unit_box else 1.
        self.length = 1. if unit_box else 100/height
        self.dx = self.length/n
        self.a = 13.159 if unit_box else .0013159/height**2
        self.phi = .01 if unit_box else .01/height**2
        self.slip = .2 if unit_box else .2/height
        k = np.fft.fftfreq(n)*n
        self.lam = 4/self.dx**2*(np.sin(np.pi*k[:,None]/n)**2 + np.sin(np.pi*k[None,:]/n)**2)
        mobility0 = self.mean_height**3+self.slip*self.mean_height**2
        self.omega = mobility0*self.lam*(self.lam-3*self.a/self.mean_height**4)

    def lap(self, h):
        return sum(np.roll(h,1,axis)+np.roll(h,-1,axis)-2*h for axis in (0,1))/self.dx**2

    def chemical(self, h):
        if not np.all(np.isfinite(h)) or np.min(h)<=0:
            raise ValueError('Nonpositive/nonfinite midpoint; no floor repair is allowed')
        return -self.lap(h)+self.a/h**3

    def energy(self, h):
        gradient = sum(((np.roll(h,-1,axis)-h)/self.dx)**2 for axis in (0,1))
        return float(np.sum(.5*gradient-self.a/(2*h**2))*self.dx**2)

    def coefficients(self, h, axis):
        face = .5*(h+np.roll(h,-1,axis))
        return face, face**3+self.slip*face**2

    def div(self, flux, axis):
        return (flux-np.roll(flux,1,axis))/self.dx

    def mobility(self, h, v):
        return -sum(self.div(self.coefficients(h,axis)[1]*(np.roll(v,-1,axis)-v)/self.dx,axis) for axis in (0,1))

    def force(self, h):
        if self.pressure=='legacy_chain':
            return self.legacy_chain_force(h)
        return -self.mobility(h,self.chemical(h))

    def legacy_chain_force(self, h):
        """Diagnostic for old source's face-midpoint attraction chain rule.

        This differs nonlinearly from the energy-gradient potential difference.
        It is not used by the new candidate or its matched EK calculation.
        """
        pressure=-self.lap(h)
        result=np.zeros_like(h)
        for axis in (0,1):
            face,mobility=self.coefficients(h,axis)
            flux=mobility*(np.roll(pressure,-1,axis)-pressure)/self.dx
            flux-=3*self.a*mobility/face**4*(np.roll(h,-1,axis)-h)/self.dx
            result+=self.div(flux,axis)
        return result

    def stochastic(self, h, dw):
        return math.sqrt(2*self.phi)/self.dx*sum(
            self.div(np.sqrt(self.coefficients(h,axis)[1])*dw[axis],axis) for axis in (0,1))

    def derivatives(self, h, v, dw):
        mu = self.chemical(h)
        dmu = -self.lap(v)-3*self.a*v/h**4
        dforce = np.zeros_like(h)
        dnoise = np.zeros_like(h)
        for axis in (0,1):
            face,mob = self.coefficients(h,axis)
            df = .5*(v+np.roll(v,-1,axis))
            dm = (3*face**2+2*self.slip*face)*df
            if self.pressure=='energy':
                flux = dm*(np.roll(mu,-1,axis)-mu)/self.dx + mob*(np.roll(dmu,-1,axis)-dmu)/self.dx
            else:
                pressure=-self.lap(h);dpressure=-self.lap(v)
                beta=-3*self.a*mob/face**4
                dbeta=-3*self.a*(dm/face**4-4*mob*df/face**5)
                flux=dm*(np.roll(pressure,-1,axis)-pressure)/self.dx+mob*(np.roll(dpressure,-1,axis)-dpressure)/self.dx
                flux+=dbeta*(np.roll(h,-1,axis)-h)/self.dx+beta*(np.roll(v,-1,axis)-v)/self.dx
            dforce += self.div(flux,axis)
            dnoise += self.div(.5*dm/np.sqrt(mob)*dw[axis],axis)
        return dforce, math.sqrt(2*self.phi)/self.dx*dnoise

    def step(self, old, step, dw, tolerance=1e-10, maxiter=15):
        denominator = 1+.5*step*self.omega
        def precondition(v):
            return np.fft.ifft2(np.fft.fft2(v.reshape(old.shape))/denominator).real.ravel()
        pre = LinearOperator((old.size,old.size),matvec=precondition,dtype=float)
        y = old + precondition(step*self.force(old)+self.stochastic(old,dw)).reshape(old.shape)
        evaluations = 0
        def residual(new):
            nonlocal evaluations
            evaluations += 1
            mid = .5*(old+new)
            return new-old-step*self.force(mid)-self.stochastic(mid,dw)
        for iteration in range(maxiter+1):
            res = residual(y)
            norm = float(np.max(np.abs(res)))
            if norm<=tolerance:
                return y,dict(newton=iteration,residual=norm,evaluations=evaluations)
            if iteration==maxiter:
                raise RuntimeError('Newton iteration limit: residual=%g'%norm)
            mid = .5*(old+y)
            def jacobian(v):
                v=v.reshape(old.shape)
                df,dn=self.derivatives(mid,v,dw)
                return (v-.5*step*df-.5*dn).ravel()
            jac=LinearOperator((old.size,old.size),matvec=jacobian,dtype=float)
            correction,info=gmres(jac,-res.ravel(),M=pre,atol=1e-13,restart=25,maxiter=8,
                                  callback_type='pr_norm',**{GMRES_RELATIVE_KEY:min(.05,max(1e-9,.1*norm))})
            if info:
                raise RuntimeError('GMRES did not converge: %s'%info)
            alpha=1.
            for _ in range(15):
                trial=y+alpha*correction.reshape(old.shape)
                if np.min(.5*(old+trial))>0:
                    trial_norm=float(np.max(np.abs(residual(trial))))
                    if trial_norm<(1-1e-4*alpha)*norm:
                        y=trial
                        break
                alpha*=.5
            else:
                raise RuntimeError('Newton line search failed')


def self_test():
    rng=np.random.default_rng(20261021)
    film=Film(1.1,8)
    h=1+.03*rng.standard_normal((8,8))
    v=rng.standard_normal(h.shape)
    dw=.1*rng.standard_normal((2,8,8))
    eps=1e-6
    analytic,noise_analytic=film.derivatives(h,v,dw)
    finite=(film.force(h+eps*v)-film.force(h-eps*v))/(2*eps)
    noise_finite=(film.stochastic(h+eps*v,dw)-film.stochastic(h-eps*v,dw))/(2*eps)
    assert np.max(abs(analytic-finite))<1e-9
    assert np.max(abs(noise_analytic-noise_finite))<1e-9
    old=Film(1.1,8,pressure='legacy_chain')
    old_jacobian,_=old.derivatives(h,v,dw)
    old_finite=(old.force(h+eps*v)-old.force(h-eps*v))/(2*eps)
    assert np.max(abs(old_jacobian-old_finite))<1e-9
    de=(film.energy(h+eps*v)-film.energy(h-eps*v))/(2*eps)
    exact=film.dx**2*np.sum(film.chemical(h)*v)
    assert abs(de-exact)<1e-8
    # Independent full covariance construction on a small grid.
    small=Film(1.1,4)
    hs=1+.02*rng.standard_normal((4,4))
    columns=[]
    drift=np.zeros(hs.size)
    for axis in (0,1):
        for j in range(16):
            basis=np.zeros((2,4,4)); basis[axis].ravel()[j]=1.
            b=small.stochastic(hs,basis)
            columns.append(b.ravel())
            # Jacobian of each diffusion column along that same column.
            _,db=small.derivatives(hs,b,basis)
            drift+=db.ravel()
    b=np.column_stack(columns)
    k=np.column_stack([small.mobility(hs,np.eye(16)[:,i].reshape(4,4)).ravel() for i in range(16)])
    expected=2*small.phi/small.dx**2*k
    covariance_error=float(np.max(abs(b@b.T-expected)))
    assert covariance_error<1e-15 and np.max(abs(drift))<1e-15
    assert np.max(abs(k-k.T))<1e-15 and np.max(abs(k@np.ones(16)))<1e-15
    # Actual nonlinear step and its independently recomputed physical residual.
    new,diagnostic=film.step(h,.05,dw)
    residual=new-h-.05*film.force(.5*(h+new))-film.stochastic(.5*(h+new),dw)
    assert np.max(abs(residual))<1e-10
    assert abs(np.mean(new)-np.mean(h))<1e-12
    # The converged nonlinear deterministic scheme linearizes to CN, not BE.
    film64=Film(1.2,64)
    mode=np.cos(2*np.pi*17*np.arange(64)[:,None]/64)*np.ones((1,64))
    amplitude=1e-6
    initial=1+amplitude*mode
    propagated,check=film64.step(initial,.03,np.zeros((2,64,64)),tolerance=1e-13)
    ratio=np.sum((propagated-1)*mode)/(amplitude*np.sum(mode**2))
    x=film64.omega[17,0]*.03
    expected_ratio=(1-x/2)/(1+x/2)
    assert abs(ratio-expected_ratio)<1e-7
    unit=Film(1.1,8,unit_box=True)
    assert abs(unit.energy(1.1*h)-1.1**2*film.energy(h))<1e-12
    assert np.max(abs(unit.force(1.1*h)-1e8*film.force(h)))<1e-10
    assert np.max(abs(unit.legacy_chain_force(1.1*h)-1e8*film.legacy_chain_force(h)))<1e-10
    return dict(analytic_Jacobian='PASS',energy_gradient='PASS',covariance_max_abs_error=covariance_error,
                Ito_Stratonovich_drift_max_abs=float(np.max(abs(drift))),nonlinear_step=diagnostic,
                volume_mean_difference=float(np.mean(new)-np.mean(h)),linear_CN_ratio=float(ratio),
                linear_CN_expected=float(expected_ratio),unit_box_energy_and_time_scaling='PASS',nonlinear_rupture_validated=False)


def probe(args):
    assert abs(args.time/.15-round(args.time/.15))<1e-10, 'End time must be a multiple of .15'
    args.out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).read_bytes()
    (args.out/'source_snapshot.py').write_bytes(source)
    rows=[]
    steps=(.003,.01,.03,.05)
    start=time.perf_counter()
    for height in args.heights:
        film=Film(height,args.n,pressure=args.pressure)
        for seed in args.seeds:
            rng=np.random.default_rng(seed)
            fields={step:np.ones((args.n,args.n)) for step in steps}
            mode_indices=np.array([(1,0),(args.n//8,0),(args.n//4,0),(args.n//2,0),(args.n//2,args.n//2)])
            traces={step:[np.zeros(5,dtype=complex)] for step in steps}
            records={step:dict(state='RUNNING',t=0.,max_residual=0.,max_newton=0,volume_error=0.,wall_seconds=0.,steps=0) for step in steps}
            for block in range(round(args.time/.15)):
                increments=math.sqrt(.001)*rng.standard_normal((150,2,args.n,args.n))
                for step in steps:
                    record=records[step]
                    if record['state']!='RUNNING':
                        continue
                    window=round(step/.001)
                    for j in range(0,150,window):
                        dw=increments[j:j+window].sum(axis=0)
                        clock=time.perf_counter()
                        try:
                            if args.linearized:
                                old=fields[step]
                                stochastic=film.stochastic(np.ones_like(old),dw)
                                numerator=(1-.5*step*film.omega)*np.fft.fft2(old-1)+np.fft.fft2(stochastic)
                                new=1+np.fft.ifft2(numerator/(1+.5*step*film.omega)).real
                                average=.5*(old+new)-1
                                linear_force=np.fft.ifft2(-film.omega*np.fft.fft2(average)).real
                                d=dict(newton=0,residual=float(np.max(abs(new-old-step*linear_force-stochastic))))
                            else:
                                new,d=film.step(fields[step],step,dw)
                        except (ValueError,RuntimeError) as exc:
                            record.update(state='NUMERICAL_FAILURE',error=str(exc))
                            break
                        record['wall_seconds']+=time.perf_counter()-clock
                        record['max_residual']=max(record['max_residual'],d['residual'])
                        record['max_newton']=max(record['max_newton'],d['newton'])
                        record['volume_error']=max(record['volume_error'],abs(float(np.mean(new))-1))
                        record['t']+=step
                        record['steps']+=1
                        fields[step]=new
                        if args.record_modes:
                            transformed=np.fft.fft2(new)/args.n**2
                            traces[step].append(transformed[mode_indices[:,0],mode_indices[:,1]])
                        if new.min()<=0:
                            record.update(state='NUMERICAL_OVERSHOOT',error='Nonpositive endpoint; no rupture claim')
                            break
                        if new.min()<=.01:
                            record.update(state='RUPTURE_BRACKET',event_lower=record['t']-step,event_upper=record['t'])
                            break
            for step in steps:
                h=fields[step]; record=records[step]
                if record['state']=='RUNNING':
                    record['state']='CENSORED_AT_PROBE_END'
                row=dict(height=height,n=args.n,seed=seed,dt=step,**record,
                         hmin=float(h.min()),hmax=float(h.max()),rms=float(np.std(h)),energy=film.energy(h),
                         rms_difference_to_dt003=float(np.sqrt(np.mean((h-fields[.003])**2))))
                rows.append(row)
                np.savez_compressed(args.out/('h%g_seed%d_dt%g.npz'%(height,seed,step)),h=h,
                                    rng_state=json.dumps(rng.bit_generator.state),t=record['t'],restartable=False,
                                    mode_indices=mode_indices,mode_series=np.array(traces[step]),
                                    mode_times=step*np.arange(len(traces[step])))
            print('Finished h=%g seed=%d elapsed=%.2fs'%(height,seed,time.perf_counter()-start),flush=True)
    (args.out/'probe.json').write_text(json.dumps(dict(rows=rows,common_Brownian_base_dt=.001,
        macro_block=.15,clipping=False,noise_resampling=False,wall_seconds=time.perf_counter()-start,
        source_sha256=hashlib.sha256(source).hexdigest(),
        args=dict(n=args.n,time=args.time,heights=args.heights,seeds=args.seeds,pressure=args.pressure,linearized=args.linearized)),indent=2),encoding='utf-8')
    keys=sorted(set().union(*(row.keys() for row in rows)))
    with (args.out/'probe.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=keys); writer.writeheader(); writer.writerows(rows)
    if any(row['state'].startswith('NUMERICAL_') for row in rows):
        raise RuntimeError('Numerical failures preserved in probe.json; no successful-gate claim')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--n',type=int,default=64)
    parser.add_argument('--time',type=float,default=.3)
    parser.add_argument('--heights',type=float,nargs='+',default=[1.,1.05,1.1,1.15,1.2])
    parser.add_argument('--seeds',type=int,nargs='+',default=[2026102101,2026102102])
    parser.add_argument('--out',type=Path)
    parser.add_argument('--record-modes',action='store_true')
    parser.add_argument('--pressure',choices=['energy','legacy_chain'],default='energy')
    parser.add_argument('--linearized',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        result=self_test()
        print(json.dumps(result,indent=2))
        if args.out:
            args.out.write_text(json.dumps(result,indent=2),encoding='utf-8')
    else:
        if not args.out:
            parser.error('--out is required for an independent probe folder')
        probe(args)
