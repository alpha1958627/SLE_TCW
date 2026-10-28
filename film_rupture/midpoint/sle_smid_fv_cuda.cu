// FP64 periodic FV kernels; no RNG, clipping, fast math or volume repair.
extern "C" {
__global__ void prepare(const double* h,const double* dx,const double* a,
                        const double* slip,double* buf,int n,int total) {
 int z=blockDim.x*blockIdx.x+threadIdx.x;if(z>=total)return;
 int p=z%(n*n),b=z/(n*n),i=p/n,j=p%n,P=n*n,base=b*P;
 int xp=((i+1)%n)*n+j,xm=((i+n-1)%n)*n+j,yp=i*n+(j+1)%n,ym=i*n+(j+n-1)%n;
 const double* v=h+base;double* q=buf+b*8*P;
 double value=v[p],fx=.5*(value+v[xp]),fy=.5*(value+v[yp]);
 double lap=((v[xp]+v[xm]-2*value)+(v[yp]+v[ym]-2*value))/(dx[b]*dx[b]);
 q[p]=-lap+a[b]/(value*value*value);q[7*P+p]=value*value*value*value;
 q[P+p]=fx*fx*fx+slip[b]*fx*fx;q[2*P+p]=fy*fy*fy+slip[b]*fy*fy;
 q[3*P+p]=3*fx*fx+2*slip[b]*fx;q[4*P+p]=3*fy*fy+2*slip[b]*fy;
 q[5*P+p]=sqrt(q[P+p]);q[6*P+p]=sqrt(q[2*P+p]);
}
__global__ void drift(const double* q0,const double* dw0,const double* dx,
                      const double* phi,const double* dt,double* out,int n,int total) {
 int z=blockDim.x*blockIdx.x+threadIdx.x;if(z>=total)return;
 int P=n*n,p=z%P,b=z/P,i=p/n,j=p%n;
 int xp=((i+1)%n)*n+j,xm=((i+n-1)%n)*n+j,yp=i*n+(j+1)%n,ym=i*n+(j+n-1)%n;
 const double* q=q0+b*8*P;const double* dw=dw0+b*2*P;double d=dx[b];
 double force=(q[P+p]*(q[xp]-q[p])/d-q[P+xm]*(q[p]-q[xm])/d)/d
             +(q[2*P+p]*(q[yp]-q[p])/d-q[2*P+ym]*(q[p]-q[ym])/d)/d;
 double noise=sqrt(2*phi[b])/d*((q[5*P+p]*dw[p]-q[5*P+xm]*dw[xm])/d
             +(q[6*P+p]*dw[P+p]-q[6*P+ym]*dw[P+ym])/d);
 out[z]=dt[b]*force+noise;
}
__global__ void dchemical(const double* v,const double* q,const double* dx,
                          const double* a,double* dm,int n,int total) {
 int z=blockDim.x*blockIdx.x+threadIdx.x;if(z>=total)return;
 int P=n*n,p=z%P,b=z/P,i=p/n,j=p%n,base=b*P;
 int xp=((i+1)%n)*n+j,xm=((i+n-1)%n)*n+j,yp=i*n+(j+1)%n,ym=i*n+(j+n-1)%n;
 double lap=((v[base+xp]+v[base+xm]-2*v[z])+(v[base+yp]+v[base+ym]-2*v[z]))/(dx[b]*dx[b]);
 dm[z]=-lap-3*a[b]*v[z]/q[b*8*P+7*P+p];
}
__device__ double df(int p,int next,int axis,const double* q,const double* v,
                    const double* dm,int P,double dx) {
 double dface=q[(3+axis)*P+p]*.5*(v[p]+v[next]);
 return dface*(q[next]-q[p])/dx+q[(1+axis)*P+p]*(dm[next]-dm[p])/dx;
}
__device__ double dn(int p,int next,int axis,const double* q,const double* v,
                    const double* dw,int P) {
 double dface=q[(3+axis)*P+p]*.5*(v[p]+v[next]);
 return .5*dface/q[(5+axis)*P+p]*dw[axis*P+p];
}
__global__ void jacobian(const double* v0,const double* q0,const double* dw0,
                         const double* dm0,const double* dx,const double* phi,
                         const double* dt,double* out,int n,int total) {
 int z=blockDim.x*blockIdx.x+threadIdx.x;if(z>=total)return;
 int P=n*n,p=z%P,b=z/P,i=p/n,j=p%n;
 int xp=((i+1)%n)*n+j,xm=((i+n-1)%n)*n+j,yp=i*n+(j+1)%n,ym=i*n+(j+n-1)%n;
 const double* q=q0+b*8*P;const double* v=v0+b*P;const double* dw=dw0+b*2*P;
 const double* dm=dm0+b*P;double d=dx[b];
 double f=(df(p,xp,0,q,v,dm,P,d)-df(xm,p,0,q,v,dm,P,d))/d
         +(df(p,yp,1,q,v,dm,P,d)-df(ym,p,1,q,v,dm,P,d))/d;
 double s=sqrt(2*phi[b])/d*((dn(p,xp,0,q,v,dw,P)-dn(xm,p,0,q,v,dw,P))/d
         +(dn(p,yp,1,q,v,dw,P)-dn(ym,p,1,q,v,dw,P))/d);
 out[z]=v[p]-.5*dt[b]*f-.5*s;
}
__global__ void rotate(double* H,double* cs,double* sn,double* g,
                       unsigned char* active,const double* threshold,int* iters,
                       int j,int m,int B) {
 int b=blockDim.x*blockIdx.x+threadIdx.x;if(b>=B)return;
 double* h=H+b*(m+1)*m;double* c=cs+b*m;double* s=sn+b*m;double* e=g+b*(m+1);
 if(!active[b]){h[j*m+j]=1;h[(j+1)*m+j]=0;return;}
 for(int i=0;i<j;i++){
   double x=h[i*m+j],y=h[(i+1)*m+j];
   h[i*m+j]=c[i]*x+s[i]*y;h[(i+1)*m+j]=-s[i]*x+c[i]*y;
 }
 double x=h[j*m+j],y=h[(j+1)*m+j],r=hypot(x,y);
 c[j]=r>0?x/r:1;s[j]=r>0?y/r:0;
 h[j*m+j]=r>0?r:1;h[(j+1)*m+j]=0;
 double rhs=e[j];e[j]=c[j]*rhs;e[j+1]=-s[j]*rhs;iters[b]++;
 if(fabs(e[j+1])<=threshold[b] || r==0)active[b]=0;
}
__global__ void backsolve(const double* H,const double* g,double* coeff,int m,int B) {
 int b=blockDim.x*blockIdx.x+threadIdx.x;if(b>=B)return;
 const double* h=H+b*(m+1)*m;const double* e=g+b*(m+1);double* y=coeff+b*m;
 for(int i=m-1;i>=0;i--){double rhs=e[i];for(int j=i+1;j<m;j++)rhs-=h[i*m+j]*y[j];y[i]=rhs/h[i*m+i];}
}
// One block per independent lane. Both Gram-Schmidt passes, column norm and
// all Givens rotations are fused; no host Hessenberg matrix or scalar reads.
__global__ void arnoldi(double* w,double* V,double* H,double* cs,double* sn,
                        double* g,unsigned char* active,const double* threshold,
                        int* iters,int j,int m,int B,int P) {
 int b=blockIdx.x,t=threadIdx.x,T=blockDim.x;
 extern __shared__ double sh[];
 double* h=H+b*(m+1)*m;double* ww=w+b*P;
 for(int pass=0;pass<2;pass++){
   for(int k=0;k<=j;k++){
     double value=0;const double* v=V+(k*B+b)*P;
     for(int p=t;p<P;p+=T)value+=ww[p]*v[p];
     sh[k*T+t]=value;
   }
   __syncthreads();
   for(int stride=T/2;stride>0;stride/=2){
     if(t<stride)for(int k=0;k<=j;k++)sh[k*T+t]+=sh[k*T+t+stride];
     __syncthreads();
   }
   if(t==0)for(int k=0;k<=j;k++)h[k*m+j]+=sh[k*T];
   for(int p=t;p<P;p+=T){double projection=0;
     for(int k=0;k<=j;k++)projection+=sh[k*T]*V[(k*B+b)*P+p];
     ww[p]-=projection;
   }
   __syncthreads();
 }
 double square=0;for(int p=t;p<P;p+=T)square+=ww[p]*ww[p];
 __syncthreads();sh[t]=square;__syncthreads();
 for(int stride=T/2;stride>0;stride/=2){if(t<stride)sh[t]+=sh[t+stride];__syncthreads();}
 if(t==0){
   double beta=sqrt(sh[0]);h[(j+1)*m+j]=beta;sh[T]=beta;
   double* c=cs+b*m;double* s=sn+b*m;double* e=g+b*(m+1);
   if(!active[b]){h[j*m+j]=1;h[(j+1)*m+j]=0;}
   else{
     for(int i=0;i<j;i++){double x=h[i*m+j],y=h[(i+1)*m+j];h[i*m+j]=c[i]*x+s[i]*y;h[(i+1)*m+j]=-s[i]*x+c[i]*y;}
     double x=h[j*m+j],y=h[(j+1)*m+j],r=hypot(x,y);
     c[j]=r>0?x/r:1;s[j]=r>0?y/r:0;h[j*m+j]=r>0?r:1;h[(j+1)*m+j]=0;
     double rhs=e[j];e[j]=c[j]*rhs;e[j+1]=-s[j]*rhs;iters[b]++;
     if(fabs(e[j+1])<=threshold[b] || r==0)active[b]=0;
   }
 }
 __syncthreads();double beta=sh[T];
 for(int p=t;p<P;p+=T)V[((j+1)*B+b)*P+p]=(active[b] && beta>1e-300)?ww[p]/beta:0;
}
__global__ void gmres_update(const double* H,const double* g,const double* V,
                             double* x,int m,int B,int P) {
 int b=blockIdx.x,t=threadIdx.x,T=blockDim.x;extern __shared__ double coeff[];
 if(t==0){
   const double* h=H+b*(m+1)*m;const double* e=g+b*(m+1);
   for(int i=m-1;i>=0;i--){double rhs=e[i];for(int j=i+1;j<m;j++)rhs-=h[i*m+j]*coeff[j];coeff[i]=rhs/h[i*m+i];}
 }
 __syncthreads();
 for(int p=t;p<P;p+=T){double correction=0;for(int j=0;j<m;j++)correction+=coeff[j]*V[(j*B+b)*P+p];x[b*P+p]+=correction;}
}
}
