"""Resident, independent-lane FP64 Newton/GMRES on one CUDA GPU.

No CuPy GMRES call: Arnoldi, Givens and back-substitution stay on device.
The space/time equations equal the existing reference; the inner solver differs.
"""
from pathlib import Path
import numpy as np
import cupy as cp

ROOT=Path(__file__).resolve().parent

class BatchedFilm:
    def __init__(self,heights,n,restart=8):
        self.B=len(heights);self.n=n;self.restart=restart;self.shape=(self.B,n,n)
        h=cp.asarray(heights,dtype=cp.float64)
        self.dx=100/h/n;self.a=.0013159/h**2;self.phi=.01/h**2;self.slip=.2/h
        k=cp.fft.fftfreq(n)*n
        lam=4/self.dx[:,None,None]**2*(cp.sin(cp.pi*k[None,:,None]/n)**2+cp.sin(cp.pi*k[None,None,:]/n)**2)
        self.omega=(1+self.slip[:,None,None])*lam*(lam-3*self.a[:,None,None])
        self.buf=cp.empty((self.B,8,n,n),dtype=cp.float64)
        self.dmu=cp.empty(self.shape,dtype=cp.float64)
        self.module=cp.RawModule(code=(ROOT.parent/'sle_smid_fv_cuda.cu').read_text(),options=('--std=c++11','--fmad=false'),
                                 name_expressions=['prepare','drift','dchemical','jacobian','arnoldi','gmres_update'])
        self.kernels={name:self.module.get_function(name) for name in ['prepare','drift','dchemical','jacobian','arnoldi','gmres_update']}
        self.synchronizations=0;self.diagnostic_bytes=0

    def launch(self,name,args,total=None):
        total=self.B*self.n*self.n if total is None else total
        self.kernels[name](((total+255)//256,),(256,),args)

    def all(self,mask):
        self.synchronizations+=1
        return bool(cp.all(mask))

    def prepare(self,h):
        self.launch('prepare',(h,self.dx,self.a,self.slip,self.buf,np.int32(self.n),np.int32(h.size)))

    def drift(self,dw,dt):
        out=cp.empty(self.shape,dtype=cp.float64)
        self.launch('drift',(self.buf,dw,self.dx,self.phi,dt,out,np.int32(self.n),np.int32(out.size)))
        return out

    def jacobian(self,v,dw,dt):
        v=cp.ascontiguousarray(v);out=cp.empty_like(v)
        self.launch('dchemical',(v,self.buf,self.dx,self.a,self.dmu,np.int32(self.n),np.int32(v.size)))
        self.launch('jacobian',(v,self.buf,dw,self.dmu,self.dx,self.phi,dt,out,np.int32(self.n),np.int32(v.size)))
        return out

    @staticmethod
    def norm(v):return cp.sqrt(cp.sum(v*v,axis=(-2,-1)))

    def gmres(self,rhs,dw,dt,denominator,relative,enabled):
        def pre(v):return cp.fft.ifft2(cp.fft.fft2(v,axes=(-2,-1))/denominator,axes=(-2,-1)).real
        def operator(v):return pre(self.jacobian(v,dw,dt))
        b=pre(rhs);threshold=cp.maximum(1e-13,relative*self.norm(b))
        x=cp.zeros_like(b);done=~enabled;iterations=cp.zeros(self.B,dtype=cp.int32)
        m=self.restart
        # 25*8=200 is the reference's 8*25 Krylov-vector budget.
        for cycle in range((200+m-1)//m):
            residual=b-operator(x) if cycle else b.copy()
            beta=self.norm(residual);done|=beta<=threshold
            if self.all(done):break
            V=cp.zeros((m+1,self.B,self.n,self.n),dtype=cp.float64)
            H=cp.zeros((self.B,m+1,m),dtype=cp.float64)
            cs=cp.zeros((self.B,m),dtype=cp.float64);sn=cp.zeros_like(cs)
            g=cp.zeros((self.B,m+1),dtype=cp.float64);g[:,0]=cp.where(done,0,beta)
            active=(~done).astype(cp.uint8)
            V[0]=cp.where(done[:,None,None],0,residual/cp.maximum(beta,1e-300)[:,None,None])
            for j in range(m):
                w=cp.ascontiguousarray(operator(V[j]))
                self.kernels['arnoldi']((self.B,),(256,),
                    (w,V,H,cs,sn,g,active,threshold,iterations,np.int32(j),np.int32(m),np.int32(self.B),np.int32(self.n*self.n)),
                    shared_mem=8*256*(m+1))
            self.kernels['gmres_update']((self.B,),(256,),
                (H,g,V,x,np.int32(m),np.int32(self.B),np.int32(self.n*self.n)),shared_mem=8*m)
            done|=self.norm(b-operator(x))<=threshold
            if self.all(done):break
        return x,done,iterations

    def step(self,old,dt,dw,enabled=None,tolerance=1e-10,maxiter=15):
        dt=cp.asarray(dt,dtype=cp.float64);dw=cp.ascontiguousarray(dw,dtype=cp.float64)
        # Keep device references, with no copies, for independent last-step
        # residual reconstruction after a timed run.
        self.last_old=old;self.last_dw=dw;self.last_dt=dt
        enabled=cp.ones(self.B,dtype=cp.bool_) if enabled is None else cp.asarray(enabled,dtype=cp.bool_)
        denominator=1+.5*dt[:,None,None]*self.omega
        self.prepare(old);increment=self.drift(dw,dt)
        y=old+cp.fft.ifft2(cp.fft.fft2(increment,axes=(-2,-1))/denominator,axes=(-2,-1)).real
        done=~enabled;failed=cp.zeros(self.B,dtype=cp.bool_)
        counts=cp.zeros(self.B,dtype=cp.int32);linear_iters=cp.zeros_like(counts)
        norm=cp.zeros(self.B,dtype=cp.float64)
        for iteration in range(maxiter+1):
            mid=.5*(old+y)
            valid=cp.all(cp.isfinite(mid)&(mid>0),axis=(-2,-1));failed|=~valid&~done
            safe_mid=cp.where((done|failed)[:,None,None],old,mid)
            self.prepare(safe_mid);res=y-old-self.drift(dw,dt)
            norm=cp.max(cp.abs(res),axis=(-2,-1))
            done|=(norm<=tolerance)&~failed
            if self.all(done|failed):break
            if iteration==maxiter:failed|=~done;break
            relative=cp.minimum(.05,cp.maximum(1e-9,.1*norm))
            correction,linear_ok,krylov=self.gmres(cp.where((done|failed)[:,None,None],0,-res),dw,dt,denominator,relative,~(done|failed))
            linear_iters+=krylov;failed|=~linear_ok&~done
            accepted=done|failed;alpha=cp.ones(self.B,dtype=cp.float64)
            for attempt in range(15):
                trial=y+alpha[:,None,None]*correction;tmid=.5*(old+trial)
                good=cp.all(cp.isfinite(tmid)&(tmid>0),axis=(-2,-1))
                self.prepare(cp.where(good[:,None,None],tmid,old))
                trialnorm=cp.max(cp.abs(trial-old-self.drift(dw,dt)),axis=(-2,-1))
                take=good&(trialnorm<(1-1e-4*alpha)*norm)&~accepted
                y=cp.where(take[:,None,None],trial,y);counts+=take.astype(cp.int32)
                accepted|=take
                if self.all(accepted):break
                alpha=cp.where(accepted,alpha,.5*alpha)
            failed|=~accepted
        # Do not replace a converged lane's midpoint with old in its final
        # physical residual. Each returned candidate is checked independently.
        mid=.5*(old+y);valid=cp.all(cp.isfinite(mid)&(mid>0),axis=(-2,-1))
        self.prepare(cp.where(valid[:,None,None],mid,old))
        norm=cp.max(cp.abs(y-old-self.drift(dw,dt)),axis=(-2,-1))
        failed|=(norm>tolerance)|~cp.isfinite(norm)|~valid
        finite=cp.all(cp.isfinite(y),axis=(-2,-1))
        diagnostic=cp.stack([norm,counts.astype(cp.float64),failed.astype(cp.float64),
            cp.min(y,axis=(-2,-1)),cp.max(y,axis=(-2,-1)),cp.mean(y,axis=(-2,-1)),
            finite.astype(cp.float64),linear_iters.astype(cp.float64)],axis=1).get()
        self.diagnostic_bytes+=diagnostic.nbytes;self.synchronizations+=1
        return y,diagnostic
