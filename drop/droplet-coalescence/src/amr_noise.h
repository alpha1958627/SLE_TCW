#ifndef AMR_NOISE_H
#define AMR_NOISE_H
#define AMR_FACE_STRIP_PROJECTION 1
#include "ou_kernel.h"
/* Each hierarchy entry averages fine-face flux ALONG the coarse face.
   The normal control width stays fixed at the finest fluid width.
   Therefore restriction commutes with finite-volume divergence.
   This is a rectangular strip average, not a point interpolation.
   A background at fluid_max+1 aligns all dyadic face-control-volume boxes.
   Shared coarse/fine faces are restricted from the actual fine fluxes by
   the tree's face restriction; the unused coarse-box value is not imposed.
   All levels retain one global normalization fixed at fluid_max. */
static int amr_noise_min=0,amr_noise_max=0,amr_latent_n=0;
static size_t amr_offsets[20],amr_pack_size=0;
static double *amr_pack[2]={NULL,NULL},*amr_prefix=NULL,*amr_c1=NULL;
static double amr_noise_sd2=1.,amr_czero=0.,amr_norm2=1.;
static double amr_gen_seconds=0.,amr_project_seconds=0.,amr_bcast_seconds=0.;
static double amr_clock (void) {
#if _MPI
  return MPI_Wtime();
#else
  struct timespec z;clock_gettime(CLOCK_MONOTONIC,&z);return z.tv_sec+1e-9*z.tv_nsec;
#endif
}
static void amr_noise_clear (void) {
  for(int c=0;c<2;c++){free(tcw_se_field[c]);free(tcw_se_endpoint[c]);tcw_se_field[c]=tcw_se_endpoint[c]=NULL;free(amr_pack[c]);amr_pack[c]=NULL;}
  free(tcw_se_first);free(tcw_se_innovation);free(tcw_se_filter);free(tcw_se_fft);free(tcw_se_col);
  tcw_se_first=tcw_se_innovation=tcw_se_filter=NULL;tcw_se_fft=tcw_se_col=NULL;
  free(amr_prefix);free(amr_c1);amr_prefix=amr_c1=NULL;tcw_noise_n=tcw_noise_initialized=0;
}
static void amr_noise_configure (int minimum,int maximum,int mode,
                                 double sigma,double tau,uint64_t seed) {
  amr_noise_clear();amr_noise_min=minimum;amr_noise_max=maximum;
  LEVEL_MIN=tcw_ouse_level=maximum+1;amr_latent_n=1<<LEVEL_MIN;
  tcw_noise_mode=mode;tcw_sigma=sigma;tcw_tau=tau;tcw_noise_seed=seed;
  amr_pack_size=0;
  if(mode<0)return;
  for(int l=minimum;l<=maximum;l++){amr_offsets[l]=amr_pack_size;amr_pack_size+=(size_t)(1<<l)*(1<<l);}
  for(int c=0;c<2;c++)amr_pack[c]=calloc(amr_pack_size,sizeof(double));
  amr_c1=calloc(amr_latent_n,sizeof(double));
  if(!amr_pack[0]||!amr_pack[1]||!amr_c1)exit(70);
  double pars[3]={.5,0.,1.};
  if(pid()==0) {
    amr_prefix=calloc((size_t)(amr_latent_n+1)*(amr_latent_n+1),sizeof(double));if(!amr_prefix)exit(70);
    if(mode==2) {
      tcw_se_ensure();
      double complex *q=malloc(amr_latent_n*sizeof(double complex));if(!q)exit(70);
      for(int k=0;k<amr_latent_n;k++)q[k]=sq(tcw_se_filter[k]);
      tcw_fft1d(q,amr_latent_n,1);
      for(int k=0;k<amr_latent_n;k++)amr_c1[k]=creal(q[k]);
      free(q);
      pars[1]=sq(sq(tcw_se_filter[0])/amr_latent_n);
      pars[2]=sq(amr_c1[0])-pars[1];
      pars[0]=sqrt((sq((amr_c1[0]+amr_c1[1])*.5)-pars[1])/pars[2]);
    } else amr_c1[0]=1.;
  }
#if _MPI
  MPI_Bcast(pars,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(amr_c1,amr_latent_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
  amr_noise_sd2=pars[0];amr_czero=pars[1];amr_norm2=pars[2];
  if(!(amr_noise_sd2>0.)||!isfinite(amr_noise_sd2))exit(71);
}
static double amr_rect (int x0,int x1,int y0,int y1) {
  int n=amr_latent_n+1;
  return amr_prefix[(size_t)y1*n+x1]-amr_prefix[(size_t)y0*n+x1]
        -amr_prefix[(size_t)y1*n+x0]+amr_prefix[(size_t)y0*n+x0];
}
static double amr_box (int sx,int sy,int mx,int my) {
  int n=amr_latent_n;
  sx=tcw_wrap_n(sx,n);sy=tcw_wrap_n(sy,n);
  int nx=min(mx,n-sx),ny=min(my,n-sy);
  double sum=amr_rect(sx,sx+nx,sy,sy+ny);
  if(nx<mx)sum+=amr_rect(0,mx-nx,sy,sy+ny);
  if(ny<my)sum+=amr_rect(sx,sx+nx,0,my-ny);
  if(nx<mx&&ny<my)sum+=amr_rect(0,mx-nx,0,my-ny);
  return sum/(mx*(double)my*amr_noise_sd2);
}
static void amr_noise_prepare (int step,double dt) {
  double start=amr_clock();
  if(pid()==0) {
    if(tcw_noise_mode==2)tcw_prepare_noise_local(step,dt);
    amr_gen_seconds+=amr_clock()-start;start=amr_clock();
    int n=amr_latent_n,stride=n+1;
    for(int c=0;c<2;c++) {
      memset(amr_prefix,0,(size_t)stride*stride*sizeof(double));
      for(int j=0;j<n;j++) {
        double row=0.;
        for(int i=0;i<n;i++) {
          double z=tcw_noise_mode==2?tcw_se_field[c][(size_t)j*n+i]:tcw_hash_gaussian(i,j,c,step,n,0x243f6a8885a308d3ULL);
          row+=z;amr_prefix[(size_t)(j+1)*stride+i+1]=row+amr_prefix[(size_t)j*stride+i+1];
        }
      }
      for(int l=amr_noise_min;l<=amr_noise_max;l++) {
        int nn=1<<l,m=n/nn;
        for(int j=0;j<nn;j++)for(int i=0;i<nn;i++) {
          int sx=i*m-(c==0?1:0),sy=j*m-(c==1?1:0);
          amr_pack[c][amr_offsets[l]+(size_t)j*nn+i]=amr_box(sx,sy,c==0?2:m,c==1?2:m);
        }
      }
    }
    amr_project_seconds+=amr_clock()-start;
  }
  start=amr_clock();
#if _MPI
  for(int c=0;c<2;c++)if(MPI_Bcast(amr_pack[c],amr_pack_size,MPI_DOUBLE,0,MPI_COMM_WORLD)!=MPI_SUCCESS)MPI_Abort(MPI_COMM_WORLD,72);
#endif
  amr_bcast_seconds+=amr_clock()-start;
}
static double amr_noise_at (double xx,double yy,int component,int l) {
  if(l<amr_noise_min||l>amr_noise_max)exit(73);
  int n=1<<l;double d=L0/n;
  int i=(int)llround((xx-X0)/d-(component?.5:0.));
  int j=(int)llround((yy-Y0)/d-(component?0.:.5));
  return amr_pack[component][amr_offsets[l]+(size_t)tcw_wrap_n(j,n)*n+tcw_wrap_n(i,n)];
}
/* Exact covariance of two projected boxes of the represented Gaussian law.
   c1 is independently obtained from the actual nonnegative spectrum. */
static double amr_noise_cov (int la,int ia,int ja,int lb,int ib,int jb,int c) {
  int ma=amr_latent_n/(1<<la),mb=amr_latent_n/(1<<lb);
  int ax=ia*ma-(c==0?1:0),ay=ja*ma-(c==1?1:0);
  int bx=ib*mb-(c==0?1:0),by=jb*mb-(c==1?1:0);
  int anx=c==0?2:ma,any=c==1?2:ma,bnx=c==0?2:mb,bny=c==1?2:mb;
  double cx=0.,cy=0.;
  for(int i=0;i<anx;i++)for(int j=0;j<bnx;j++)cx+=amr_c1[tcw_wrap_n(ax+i-bx-j,amr_latent_n)];
  for(int i=0;i<any;i++)for(int j=0;j<bny;j++)cy+=amr_c1[tcw_wrap_n(ay+i-by-j,amr_latent_n)];
  cx/=anx*(double)bnx;cy/=any*(double)bny;
  return (cx*cy-amr_czero)/(amr_norm2*sq(amr_noise_sd2));
}
#endif
