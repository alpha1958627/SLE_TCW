/* Conservative stochastic implicit midpoint on Basilisk uniform periodic cells.
   Python supplies the frozen Brownian increments; the nonlinear solve is C.
   No floor, post-step mass projection or random-number generation here. */
#include "grid/multigrid.h"
#include "utils.h"

scalar bm_h[], bm_mu[], bm_v[], bm_dmu[];
face vector bm_m[], bm_dm[], bm_s[], bm_flux[];
static int bm_n = 0, bm_started = 0;
static double bm_a, bm_phi, bm_slip, bm_dx;
static int bm_spectral = 1;
static double * bm_cap = NULL, * bm_dimag = NULL;
static void bm_fft2(double * re,double * im,int inverse);
int bm_set_spectral(int enabled) { bm_spectral = enabled != 0; return bm_spectral; }
/* Positive spectral curvature: variational derivative of the Fourier
   Dirichlet energy. Mobility and noise retain the same FV face operator. */
static void bm_curvature(const double * input,double * output) {
  int n=bm_n,s=n*n;
  memcpy(output,input,s*sizeof(double));memset(bm_dimag,0,s*sizeof(double));
  bm_fft2(output,bm_dimag,0);
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
    int ki=i<=n/2?i:i-n,kj=j<=n/2?j:j-n;
    double q2=4.*sq(pi)/sq(L0)*(ki*ki+kj*kj);int k=i*n+j;
    output[k]*=q2;bm_dimag[k]*=q2;
  }
  bm_fft2(output,bm_dimag,1);
}

int bm_rank(void) { return pid(); }
int bm_ranks(void) { return npe(); }
void bm_sync(void) {
#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}
int bm_any(int value) {
#if _MPI
  int result; MPI_Allreduce(&value,&result,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);return result;
#else
  return value;
#endif
}

int bm_init(int n, double height) {
  if (!bm_started) { init_solver(); _init_solver(); bm_started = 1; }
  if (grid) free_grid();
  if (bm_cap) free(bm_cap); if (bm_dimag) free(bm_dimag);
  bm_cap=calloc(n*n,sizeof(double));bm_dimag=calloc(n*n,sizeof(double));
  bm_n = n; bm_a = .0013159/(height*height);
  bm_phi = .01/(height*height); bm_slip = .2/height;
  size(100./height); origin(0.,0.); bm_dx = L0/n;
  periodic(right); periodic(top); init_grid(n);
  long count = 0;
  foreach(reduction(+:count)) count++;
  return count == (long)n*n ? 0 : 1;
}

static int bm_index(double x, double y) {
  int i = ((int)floor(x/bm_dx + 1e-8) + bm_n)%bm_n;
  int j = ((int)floor(y/bm_dx + 1e-8) + bm_n)%bm_n;
  return i*bm_n+j;
}

static void bm_export(scalar f, double * out) {
  memset(out,0,bm_n*bm_n*sizeof(double));
  foreach(serial) out[bm_index(x,y)] = f[];
#if _MPI
  if (npe() > 1) MPI_Allreduce(MPI_IN_PLACE,out,bm_n*bm_n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
}

static int bm_prepare(const double * h) {
  for (int k=0;k<bm_n*bm_n;k++) if (!(h[k]>0) || !isfinite(h[k])) return 1;
  foreach() bm_h[] = h[bm_index(x,y)];
  boundary({bm_h});
  if (bm_spectral) bm_curvature(h,bm_cap);
  foreach() {
    double lap = (bm_h[1]+bm_h[-1]+bm_h[0,1]+bm_h[0,-1]-4*bm_h[])/(Delta*Delta);
    bm_mu[] = (bm_spectral ? bm_cap[bm_index(x,y)] : -lap) + bm_a/pow(bm_h[],3.);
  }
  boundary({bm_mu});
  foreach_face(x) {
    double f = .5*(bm_h[]+bm_h[-1]);
    bm_m.x[] = f*f*f+bm_slip*f*f;
    bm_dm.x[] = 3*f*f+2*bm_slip*f; bm_s.x[] = sqrt(bm_m.x[]);
  }
  foreach_face(y) {
    double f = .5*(bm_h[]+bm_h[0,-1]);
    bm_m.y[] = f*f*f+bm_slip*f*f;
    bm_dm.y[] = 3*f*f+2*bm_slip*f; bm_s.y[] = sqrt(bm_m.y[]);
  }
  boundary((scalar *){bm_m,bm_dm,bm_s});
  return 0;
}

/* kind 0 force, 1 noise, 2 dt*force+noise, 3 derivative of that,
   kind 4 Jacobian I - .5*(dt*dforce+dnoise), with prepared midpoint. */
static void bm_action(const double * dw, const double * v, double step, int kind, double * out) {
  if (kind>=3) {
    foreach() bm_v[] = v[bm_index(x,y)];
    boundary({bm_v});
    if (bm_spectral) bm_curvature(v,bm_cap);
    foreach() bm_dmu[] = (bm_spectral ? bm_cap[bm_index(x,y)] :
                         -(bm_v[1]+bm_v[-1]+bm_v[0,1]+bm_v[0,-1]-4*bm_v[])/(Delta*Delta))
                         -3*bm_a*bm_v[]/pow(bm_h[],4.);
    boundary({bm_dmu});
  }
  double factor = sqrt(2*bm_phi)/bm_dx;
  foreach_face(x) {
    int p = bm_index(x-Delta/2.,y);
    double force = bm_m.x[]*(bm_mu[]-bm_mu[-1])/Delta;
    double noise = factor*bm_s.x[]*dw[p];
    if (kind>=3) {
      double dm = bm_dm.x[]*.5*(bm_v[]+bm_v[-1]);
      force = dm*(bm_mu[]-bm_mu[-1])/Delta + bm_m.x[]*(bm_dmu[]-bm_dmu[-1])/Delta;
      noise = factor*.5*dm/bm_s.x[]*dw[p];
    }
    bm_flux.x[] = kind==0 ? force : kind==1 ? noise : step*force+noise;
  }
  foreach_face(y) {
    int p = bm_index(x,y-Delta/2.);
    double force = bm_m.y[]*(bm_mu[]-bm_mu[0,-1])/Delta;
    double noise = factor*bm_s.y[]*dw[bm_n*bm_n+p];
    if (kind>=3) {
      double dm = bm_dm.y[]*.5*(bm_v[]+bm_v[0,-1]);
      force = dm*(bm_mu[]-bm_mu[0,-1])/Delta + bm_m.y[]*(bm_dmu[]-bm_dmu[0,-1])/Delta;
      noise = factor*.5*dm/bm_s.y[]*dw[bm_n*bm_n+p];
    }
    bm_flux.y[] = kind==0 ? force : kind==1 ? noise : step*force+noise;
  }
  boundary((scalar *){bm_flux});
  foreach() bm_dmu[] = (bm_flux.x[1]-bm_flux.x[]+bm_flux.y[0,1]-bm_flux.y[])/Delta;
  bm_export(bm_dmu,out);
  if (kind==4) for(int k=0;k<bm_n*bm_n;k++) out[k] = v[k]-.5*out[k];
}

int bm_eval(const double * h, const double * dw, const double * v, double step, int kind, double * out) {
  if (bm_prepare(h)) return 1;
  bm_action(dw,v,step,kind,out); return 0;
}

/* Radix-two FFT preconditioner, same constant-coefficient FV symbol as the
   reference NumPy FFT preconditioner. Full small-grid vectors are replicated
   across MPI ranks; MPI scaling is measured, never assumed. */
static void bm_fft1(double * re, double * im, int n, int inverse) {
  for(int i=1,j=0;i<n;i++) {
    int bit=n>>1; for(;j&bit;bit>>=1) j^=bit; j^=bit;
    if(i<j) { double t=re[i]; re[i]=re[j]; re[j]=t; t=im[i];im[i]=im[j];im[j]=t; }
  }
  for(int len=2;len<=n;len<<=1) {
    double angle=(inverse?2.:-2.)*3.14159265358979323846/len;
    double ar=cos(angle),ai=sin(angle);
    for(int i=0;i<n;i+=len) {
      double wr=1.,wi=0.;
      for(int j=0;j<len/2;j++) {
        int u=i+j,v=i+j+len/2;
        double vr=re[v]*wr-im[v]*wi,vi=re[v]*wi+im[v]*wr;
        re[v]=re[u]-vr;im[v]=im[u]-vi;re[u]+=vr;im[u]+=vi;
        double nr=wr*ar-wi*ai;wi=wr*ai+wi*ar;wr=nr;
      }
    }
  }
  if(inverse) for(int i=0;i<n;i++) { re[i]/=n;im[i]/=n; }
}

static void bm_fft2(double * re,double * im,int inverse) {
  int n=bm_n; double rr[n],ii[n];
  for(int i=0;i<n;i++) bm_fft1(re+i*n,im+i*n,n,inverse);
  for(int j=0;j<n;j++) {
    for(int i=0;i<n;i++) { rr[i]=re[i*n+j];ii[i]=im[i*n+j]; }
    bm_fft1(rr,ii,n,inverse);
    for(int i=0;i<n;i++) { re[i*n+j]=rr[i];im[i*n+j]=ii[i]; }
  }
}

static void bm_precondition(const double * input,double * output,double step,double * imag) {
  int n=bm_n,s=n*n;memcpy(output,input,s*sizeof(double));memset(imag,0,s*sizeof(double));
  bm_fft2(output,imag,0);
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
    double sx=sin(3.14159265358979323846*i/n),sy=sin(3.14159265358979323846*j/n);
    double lam=4/(bm_dx*bm_dx)*(sx*sx+sy*sy);
    int ki=i<=n/2?i:i-n,kj=j<=n/2?j:j-n;
    double curvature=bm_spectral ? 4.*sq(pi)/sq(L0)*(ki*ki+kj*kj) : lam;
    double denom=1+.5*step*(1+bm_slip)*lam*(curvature-3*bm_a);
    int k=i*n+j;output[k]/=denom;imag[k]/=denom;
  }
  bm_fft2(output,imag,1);
}

static double bm_dot(const double * x,const double * y,int n) {
  double sum=0.;for(int k=0;k<n;k++) sum+=x[k]*y[k];return sum;
}
static double bm_max(const double * x,int n) {
  double m=0.;for(int k=0;k<n;k++) if(fabs(x[k])>m) m=fabs(x[k]);return m;
}

/* Left-preconditioned restarted GMRES. Newton always checks the original
   physical residual, regardless of this inner Krylov stopping criterion. */
static int bm_gmres(const double * rhs,const double * dw,double step,double tol,double * solution) {
  int n=bm_n*bm_n,m=25;
  double *V=calloc((m+1)*n,sizeof(double)),*w=calloc(n,sizeof(double)),*z=calloc(n,sizeof(double));
  double *r=calloc(n,sizeof(double)),*imag=calloc(n,sizeof(double));
  double H[26][25],cs[25],sn[25],g[26],y[25];
  memset(solution,0,n*sizeof(double));
  bm_precondition(rhs,r,step,imag);
  double initial=sqrt(bm_dot(r,r,n)),target=tol*initial;int success=0;
  if(initial==0) { success=1;goto finish; }
  for(int cycle=0;cycle<8;cycle++) {
    double beta=sqrt(bm_dot(r,r,n));
    if(beta<=target) { success=1;break; }
    for(int k=0;k<n;k++) V[k]=r[k]/beta;
    memset(H,0,sizeof(H));memset(g,0,sizeof(g));g[0]=beta;int used=0;
    for(int j=0;j<m;j++) {
      bm_action(dw,V+j*n,step,4,z);bm_precondition(z,w,step,imag);
      for(int i=0;i<=j;i++) {
        H[i][j]=bm_dot(w,V+i*n,n);
        for(int k=0;k<n;k++) w[k]-=H[i][j]*V[i*n+k];
      }
      H[j+1][j]=sqrt(bm_dot(w,w,n));double sub=H[j+1][j];
      if(sub>1e-30) for(int k=0;k<n;k++) V[(j+1)*n+k]=w[k]/sub;
      for(int i=0;i<j;i++) {
        double t=cs[i]*H[i][j]+sn[i]*H[i+1][j];
        H[i+1][j]=-sn[i]*H[i][j]+cs[i]*H[i+1][j];H[i][j]=t;
      }
      double den=hypot(H[j][j],H[j+1][j]);
      if(den==0) { used=j;break; }
      cs[j]=H[j][j]/den;sn[j]=H[j+1][j]/den;
      H[j][j]=den;H[j+1][j]=0;
      g[j+1]=-sn[j]*g[j];g[j]*=cs[j];used=j+1;
      if(fabs(g[j+1])<=target || sub<=1e-30) break;
    }
    if(used==0) break;
    for(int i=used-1;i>=0;i--) {
      double value=g[i];for(int j=i+1;j<used;j++) value-=H[i][j]*y[j];
      y[i]=value/H[i][i];
    }
    for(int i=0;i<used;i++) for(int k=0;k<n;k++) solution[k]+=y[i]*V[i*n+k];
    bm_action(dw,solution,step,4,z);
    for(int k=0;k<n;k++) w[k]=rhs[k]-z[k];
    bm_precondition(w,r,step,imag);
    if(sqrt(bm_dot(r,r,n))<=target) { success=1;break; }
  }
finish:
  free(V);free(w);free(z);free(r);free(imag);return success?0:1;
}

static int bm_residual(const double * old,const double * next,const double * dw,double step,
                       double * mid,double * residual,double * value) {
  int n=bm_n*bm_n;
  for(int k=0;k<n;k++) mid[k]=.5*(old[k]+next[k]);
  if(bm_prepare(mid)) return 1;
  bm_action(dw,NULL,step,2,value);
  for(int k=0;k<n;k++) residual[k]=next[k]-old[k]-value[k];
  return 0;
}

int bm_step(const double * old,double step,const double * dw,double tolerance,int maxiter,
            double * next,double * diagnostic) {
  int n=bm_n*bm_n,code=1,evaluations=0;
  double *mid=calloc(n,sizeof(double)),*res=calloc(n,sizeof(double)),*value=calloc(n,sizeof(double));
  double *corr=calloc(n,sizeof(double)),*trial=calloc(n,sizeof(double)),*rhs=calloc(n,sizeof(double));
  double *imag=calloc(n,sizeof(double));
  if(bm_prepare(old)) goto cleanup;
  bm_action(dw,NULL,step,2,value);bm_precondition(value,corr,step,imag);
  for(int k=0;k<n;k++) next[k]=old[k]+corr[k];
  for(int iteration=0;iteration<=maxiter;iteration++) {
    evaluations++;
    if(bm_residual(old,next,dw,step,mid,res,value)) goto cleanup;
    double norm=bm_max(res,n);
    diagnostic[0]=iteration;diagnostic[1]=norm;diagnostic[2]=evaluations;
    if(!isfinite(norm)) goto cleanup;
    if(norm<=tolerance) { code=0;break; }
    if(iteration==maxiter) break;
    for(int k=0;k<n;k++) rhs[k]=-res[k];
    double rtol=fmin(.05,fmax(1e-9,.1*norm));
    if(bm_gmres(rhs,dw,step,rtol,corr)) { code=2;break; }
    double alpha=1.;int accepted=0;
    for(int line=0;line<15;line++) {
      for(int k=0;k<n;k++) trial[k]=next[k]+alpha*corr[k];
      evaluations++;
      if(!bm_residual(old,trial,dw,step,mid,res,value) && bm_max(res,n)<(1-1e-4*alpha)*norm) {
        memcpy(next,trial,n*sizeof(double));accepted=1;break;
      }
      alpha*=.5;
    }
    if(!accepted) { code=3;break; }
  }
cleanup:
  free(mid);free(res);free(value);free(corr);free(trial);free(rhs);free(imag);
  return code;
}

/* qcc uses main to generate global scalar attributes; the Python driver calls
   bm_init rather than this executable entry point. */
int main(void) { return 0; }
