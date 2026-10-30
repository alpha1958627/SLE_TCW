#ifndef NONLINEAR_PAIR_H
#define NONLINEAR_PAIR_H
#include <float.h>
#include "persistent_activation.h"
/* Lagged variable mobility, energy-inner-product CG, and adjoint-stencil limiter.
   No MPI claim. Geometry and OU kernel are shared with AMR013. */
static double *nl_work=NULL,*nl_faces=NULL,*nl_latent=NULL;
static double *nl_hf,*nl_mob,*nl_amp,*nl_wa,*nl_wb,*nl_flux,*nl_noise;
static int nl_n=0,nl_ne=0;
struct NlStats {long first,hits;double raw2,removed2;};
static void nl_clear(void){free(nl_work);free(nl_faces);free(nl_latent);nl_work=nl_faces=nl_latent=NULL;nl_n=nl_ne=0;}
static void nl_geometry(void){
 nl_clear();nl_n=pm.n;nl_ne=pm.ne;nl_work=calloc(14*pm.n,sizeof(double));nl_faces=calloc(7*pm.ne,sizeof(double));nl_latent=calloc((size_t)amr_latent_n*amr_latent_n,sizeof(double));
 if(!nl_work||!nl_faces||!nl_latent)exit(80);
 nl_hf=nl_faces;nl_mob=nl_hf+pm.ne;nl_amp=nl_mob+pm.ne;nl_wa=nl_amp+pm.ne;nl_wb=nl_wa+pm.ne;nl_flux=nl_wb+pm.ne;nl_noise=nl_flux+pm.ne;
 for(int f=0;f<pm.ne;f++){
  struct PairEdge *e=&pm.edge[f];double ga=0.,gb=0.;for(int k=0;k<e->ng;k++){if(e->id[k]==e->a)ga=e->g[k]/pm.w[e->a];if(e->id[k]==e->b)gb=e->g[k]/pm.w[e->b];}
  if(!(ga<0.&&gb>0.)){fprintf(stderr,"MOBILITY_GEOMETRY_SIGN_FAIL face=%d a=%.17g b=%.17g\n",f,ga,gb);exit(81);}
  nl_wa[f]=gb/(gb-ga);nl_wb[f]=-ga/(gb-ga);
 }
}
static void nl_prepare(const double *height,int mode,double phi,double dt){
 double fine=L0/pm.N;
 double white=sqrt(2.*phi/dt)/fine;
 double ou=sqrt(phi*pow(tcw_sigma,5)*tcw_tau*sq(tcw_Upsilon(fine/tcw_sigma))*tcw_Theta(dt/tcw_tau)/(pow(pi,1.5)*pow(fine,4)*dt*dt));
 for(int f=0;f<pm.ne;f++){
  struct PairEdge *e=&pm.edge[f];double h=nl_wa[f]*height[e->a]+nl_wb[f]*height[e->b];if(!(h>0.)||!isfinite(h))exit(82);
  nl_hf[f]=h;nl_mob[f]=h*h*h;nl_amp[f]=mode<0?0.:(mode==0?white*sqrt(nl_mob[f]):ou*sqrt(fmax(tcw_Lambda(h/tcw_sigma),0.)));
 }
}
static void nl_lap(const double *u,double *v,const double *mob){
 memset(v,0,pm.n*sizeof(double));for(int f=0;f<pm.ne;f++){
  struct PairEdge *e=&pm.edge[f];double a=0.;for(int k=0;k<e->ng;k++)a+=e->g[k]*u[e->id[k]];a*=e->t*mob[f];
  for(int k=0;k<e->ng;k++)v[e->id[k]]-=e->g[k]*a/pm.w[e->id[k]];
 }
}
static void nl_apply(const double *u,double *v,double *tmp,double dt){pm_lap(u,tmp);nl_lap(tmp,v,nl_mob);for(int i=0;i<pm.n;i++)v[i]=u[i]+dt*v[i];}
static double nl_mean(const double *u){double s=0.,a=0.;for(int i=0;i<pm.n;i++){s+=pm.w[i]*u[i];a+=pm.w[i];}return s/a;}
static int nl_solve(double *u,const double *rhs,double dt,double *residual){
 int n=pm.n;double *x=nl_work,*b=x+n,*r=b+n,*v=r+n,*Av=v+n,*Lp=Av+n,*Kr=Lp+n,*Au=Kr+n,*tmp=Au+n;
 double mean=nl_mean(rhs),oldmean=nl_mean(u);for(int i=0;i<n;i++){x[i]=u[i]-oldmean;b[i]=rhs[i]-mean;}
 nl_apply(x,Au,tmp,dt);for(int i=0;i<n;i++)r[i]=b[i]-Au[i];double drift=nl_mean(r);for(int i=0;i<n;i++)r[i]-=drift;
 pm_lap(r,Kr);for(int i=0;i<n;i++){Kr[i]=-Kr[i];v[i]=r[i];}
 double rr=pm_dot(r,Kr),err=0.;for(int i=0;i<n;i++)err=fmax(err,fabs(r[i]));int iter=0;
 while(err>1e-10&&iter++<2000){
  pm_lap(v,Lp);nl_lap(Lp,Av,nl_mob);for(int i=0;i<n;i++)Av[i]=v[i]+dt*Av[i];double denom=-pm_dot(Lp,Av);
  if(!(denom>0.)||!(rr>0.)||!isfinite(denom)){fprintf(stderr,"ENERGY_CG_BREAKDOWN it=%d rr=%.17g denom=%.17g\n",iter,rr,denom);exit(83);}
  double alpha=rr/denom;for(int i=0;i<n;i++){x[i]+=alpha*v[i];r[i]-=alpha*Av[i];}
  drift=nl_mean(r);err=0.;for(int i=0;i<n;i++){r[i]-=drift;err=fmax(err,fabs(r[i]));}
  if(err<=1e-10)break;
  pm_lap(r,Kr);for(int i=0;i<n;i++)Kr[i]=-Kr[i];double next=pm_dot(r,Kr),beta=next/rr;for(int i=0;i<n;i++)v[i]=r[i]+beta*v[i];rr=next;
 }
 for(int i=0;i<n;i++)u[i]=x[i]+mean;nl_apply(u,Au,tmp,dt);*residual=0.;for(int i=0;i<n;i++)*residual=fmax(*residual,fabs(Au[i]-rhs[i]));
 if(iter>=2000||!isfinite(*residual)||*residual>1e-8){fprintf(stderr,"ENERGY_CG_FAIL it=%d residual=%.17g\n",iter,*residual);exit(84);}return iter;
}
static void nl_div(const double *flux,double *out){
 memset(out,0,pm.n*sizeof(double));for(int f=0;f<pm.ne;f++){
  struct PairEdge *e=&pm.edge[f];double z=e->length*flux[f];for(int k=0;k<e->ng;k++)out[e->id[k]]-=e->g[k]*z/pm.w[e->id[k]];
 }
}
static int nl_active_faces=0;
static void nl_random_flux(int step,double dt,int mode,const int *first_active){
 nl_active_faces=0;memset(nl_noise,0,pm.ne*sizeof(double));if(mode<0)return;if(mode==2)tcw_prepare_noise_local(step,dt);int n=amr_latent_n;
 for(int c=0;c<2;c++){
  for(int j=0;j<n;j++)for(int i=0;i<n;i++)nl_latent[j*n+i]=mode==2?tcw_se_field[c][j*n+i]:tcw_hash_gaussian(i,j,c,step,n,0x243f6a8885a308d3ULL);
  pm_prefix(nl_latent);for(int f=0;f<pm.ne;f++){struct PairEdge *e=&pm.edge[f];if(e->c==c && activation_face_allowed(e->ng,e->id,e->g,first_active)){nl_noise[f]=nl_amp[f]*amr_box(e->sx,e->sy,e->mx,e->my);nl_active_faces++;}}
 }
}
static void nl_total_flux(const double *u,double *pressure){
 pm_lap(u,pressure);for(int i=0;i<pm.n;i++)pressure[i]=-pressure[i];
 for(int f=0;f<pm.ne;f++){struct PairEdge *e=&pm.edge[f];double g=0.;for(int k=0;k<e->ng;k++)g+=e->g[k]*pressure[e->id[k]];nl_flux[f]=nl_mob[f]*g/e->distance+nl_noise[f];}
}
static void nl_limit(const double *old,double dt,double lower,int step,struct NlStats *stats){
 double *loss=nl_work+10*pm.n,*factor=loss+pm.n;memset(loss,0,pm.n*sizeof(double));
 for(int f=0;f<pm.ne;f++){struct PairEdge *e=&pm.edge[f];double z=e->length*nl_flux[f];for(int k=0;k<e->ng;k++)loss[e->id[k]]+=fmax(e->g[k]*z/pm.w[e->id[k]],0.);}
 for(int i=0;i<pm.n;i++)factor[i]=loss[i]>0.?fmin(1.,fmax(old[i]-lower,0.)/(dt*loss[i])):1.;
 for(int f=0;f<pm.ne;f++){
  struct PairEdge *e=&pm.edge[f];double before=nl_flux[f],a=1.;for(int k=0;k<e->ng;k++)if(e->g[k]*before>0.)a=fmin(a,factor[e->id[k]]);
  if(a<1.)a*=1.-64.*DBL_EPSILON;nl_flux[f]*=a;stats->raw2+=before*before*e->length;stats->removed2+=sq(nl_flux[f]-before)*e->length;
  if(a<1.-1e-12&&fabs(before)>1e-20){stats->hits++;if(stats->first<0)stats->first=step;}
 }
}
#endif
