#ifndef PAIRED_GRAPH_H
#define PAIRED_GRAPH_H
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
/* One shared face graph for curvature, dissipation and noise.
   Serial verification implementation; no claim of MPI scalability. */
struct PairEdge {int a,b,c,k,s,count,sx,sy,mx,my,ng,id[32];double length,distance,t,g[32];};
struct CellGrad {int n,id[12];double dx[12],dy[12],gx[12],gy[12];};
struct PairMesh {int n,ne,N,latent;double *x,*y,*d,*w,*diag;struct PairEdge *edge;};
static struct PairMesh pm={0};
static double pm_periodic_distance(double x){return x-L0*round(x/L0);}
static void pm_add_coefficient(struct PairEdge *e,int id,double value){
 if(value==0.)return;
 for(int k=0;k<e->ng;k++)if(e->id[k]==id){e->g[k]+=value;return;}
 if(e->ng>=32)exit(48);e->id[e->ng]=id;e->g[e->ng++]=value;
}
static void pm_add_neighbor(struct CellGrad *g,int i,int j){
 for(int k=0;k<g[i].n;k++)if(g[i].id[k]==j)return;
 int k=g[i].n++;if(k>=12)exit(49);g[i].id[k]=j;g[i].dx[k]=pm_periodic_distance(pm.x[j]-pm.x[i]);g[i].dy[k]=pm_periodic_distance(pm.y[j]-pm.y[i]);
}
static void pm_clear(void){free(pm.x);free(pm.y);free(pm.d);free(pm.w);free(pm.diag);free(pm.edge);memset(&pm,0,sizeof(pm));}
static int edge_compare(const void *aa,const void *bb){
 const struct PairEdge *a=aa,*b=bb;
 if(a->c!=b->c)return a->c-b->c;if(a->a!=b->a)return a->a-b->a;
 if(a->b!=b->b)return a->b-b->b;if(a->k!=b->k)return a->k-b->k;return a->s-b->s;
}
static void pm_build(int maximum){
 pm_clear();pm.N=1<<maximum;pm.latent=2*pm.N;
 foreach(serial)pm.n++;
 pm.x=malloc(pm.n*sizeof(double));pm.y=malloc(pm.n*sizeof(double));pm.d=malloc(pm.n*sizeof(double));pm.w=malloc(pm.n*sizeof(double));pm.diag=calloc(pm.n,sizeof(double));
 int *labels=malloc(pm.N*pm.N*sizeof(int));for(int i=0;i<pm.N*pm.N;i++)labels[i]=-1;
 double fine=L0/pm.N;int at=0;
 foreach(serial){
  pm.x[at]=x;pm.y[at]=y;pm.d[at]=Delta;pm.w[at]=sq(Delta);
  int ix=llround((x-.5*Delta-X0)/fine),iy=llround((y-.5*Delta-Y0)/fine),m=llround(Delta/fine);
  for(int j=iy;j<iy+m;j++)for(int i=ix;i<ix+m;i++){if(labels[j*pm.N+i]!=-1)exit(41);labels[j*pm.N+i]=at;}
  at++;
 }
 struct PairEdge *raw=calloc(2*pm.N*pm.N,sizeof(struct PairEdge));int nr=0;
 for(int c=0;c<2;c++)for(int k=0;k<pm.N;k++)for(int s=0;s<pm.N;s++){
  int a=c==0?labels[s*pm.N+(k+pm.N-1)%pm.N]:labels[((k+pm.N-1)%pm.N)*pm.N+s];
  int b=c==0?labels[s*pm.N+k]:labels[k*pm.N+s];
  if(a<0||b<0)exit(42);if(a==b)continue;
  raw[nr++]=(struct PairEdge){.a=a,.b=b,.c=c,.k=k,.s=s,.count=1};
 }
 free(labels);qsort(raw,nr,sizeof(struct PairEdge),edge_compare);
 pm.edge=calloc(nr,sizeof(struct PairEdge));
 for(int i=0;i<nr;i++){
  struct PairEdge *e=pm.ne?&pm.edge[pm.ne-1]:NULL;
  if(e&&e->c==raw[i].c&&e->a==raw[i].a&&e->b==raw[i].b&&e->k==raw[i].k&&e->s+e->count==raw[i].s)e->count++;
  else pm.edge[pm.ne++]=raw[i];
 }
 free(raw);
 double latent=fine/2.;int *cover=calloc(2*pm.latent*pm.latent,sizeof(int));
 for(int f=0;f<pm.ne;f++){
  struct PairEdge *e=&pm.edge[f];e->length=e->count*fine;e->distance=.5*(pm.d[e->a]+pm.d[e->b]);e->t=e->length/e->distance;
  int normal=llround((e->k*fine-.5*pm.d[e->a])/latent),tangent=2*e->s;
  int width=llround(e->distance/latent),along=2*e->count;
  e->sx=e->c==0?normal:tangent;e->sy=e->c==0?tangent:normal;
  e->mx=e->c==0?width:along;e->my=e->c==0?along:width;
  for(int j=0;j<e->my;j++)for(int i=0;i<e->mx;i++)cover[(e->c*pm.latent+tcw_wrap_n(e->sy+j,pm.latent))*pm.latent+tcw_wrap_n(e->sx+i,pm.latent)]++;
 }
 for(int i=0;i<2*pm.latent*pm.latent;i++)if(cover[i]!=1){fprintf(stderr,"DUAL_BOX_TILING_FAIL pixel=%d count=%d\n",i,cover[i]);exit(43);}free(cover);
 /* Tangential least-squares reconstruction, followed by its exact adjoint.
    Plain two-point gradients are inconsistent at non-orthogonal centres. */
 struct CellGrad *grad=calloc(pm.n,sizeof(struct CellGrad));
 for(int f=0;f<pm.ne;f++){struct PairEdge *e=&pm.edge[f];pm_add_neighbor(grad,e->a,e->b);pm_add_neighbor(grad,e->b,e->a);}
 for(int i=0;i<pm.n;i++){
  double xx=0.,xy=0.,yy=0.;for(int k=0;k<grad[i].n;k++){double x=grad[i].dx[k],y=grad[i].dy[k],w=1./(x*x+y*y);xx+=w*x*x;xy+=w*x*y;yy+=w*y*y;}
  double det=xx*yy-xy*xy;if(det<=0.)exit(50);
  for(int k=0;k<grad[i].n;k++){double x=grad[i].dx[k],y=grad[i].dy[k],w=1./(x*x+y*y);grad[i].gx[k]=w*(yy*x-xy*y)/det;grad[i].gy[k]=w*(xx*y-xy*x)/det;}
 }
 for(int f=0;f<pm.ne;f++){
  struct PairEdge *e=&pm.edge[f];pm_add_coefficient(e,e->a,-1.);pm_add_coefficient(e,e->b,1.);
  for(int side=0;side<2;side++){
   int i=side?e->b:e->a;int centre=llround((e->c==0?pm.y[i]-Y0:pm.x[i]-X0)/latent),diff=2*e->s+e->count-centre;
   if(diff>pm.latent/2)diff-=pm.latent;if(diff<-pm.latent/2)diff+=pm.latent;
   double offset=diff*latent,sign=side?1.:-1.;
   for(int k=0;k<grad[i].n;k++){
    double a=sign*offset*(e->c==0?grad[i].gy[k]:grad[i].gx[k]);pm_add_coefficient(e,grad[i].id[k],a);pm_add_coefficient(e,i,-a);
   }
  }
  int keep=0;for(int k=0;k<e->ng;k++)if(e->g[k]!=0.){e->id[keep]=e->id[k];e->g[keep++]=e->g[k];}e->ng=keep;
  for(int k=0;k<e->ng;k++)pm.diag[e->id[k]]+=e->t*sq(e->g[k])/pm.w[e->id[k]];
 }
 free(grad);
}
static void pm_read(scalar u,double *a){int i=0;foreach(serial)a[i++]=u[];}
static void pm_write(const double *a,scalar u){int i=0;foreach(serial)u[]=a[i++];boundary({u});}
static void pm_lap(const double *a,double *b){
 memset(b,0,pm.n*sizeof(double));
 for(int f=0;f<pm.ne;f++){struct PairEdge *e=&pm.edge[f];double v=0.;for(int k=0;k<e->ng;k++)v+=e->g[k]*a[e->id[k]];v*=e->t;for(int k=0;k<e->ng;k++)b[e->id[k]]-=e->g[k]*v/pm.w[e->id[k]];}
}
static double pm_dot(const double *a,const double *b){double v=0.;for(int i=0;i<pm.n;i++)v+=pm.w[i]*a[i]*b[i];return v;}
static double pm_energy(const double *a){double v=0.;for(int f=0;f<pm.ne;f++){struct PairEdge *e=&pm.edge[f];double g=0.;for(int k=0;k<e->ng;k++)g+=e->g[k]*a[e->id[k]];v+=.5*e->t*g*g;}return v;}
static void pm_apply(const double *a,double *b,double *tmp,double dt){pm_lap(a,tmp);pm_lap(tmp,b);for(int i=0;i<pm.n;i++)b[i]=a[i]+dt*b[i];}
/* Weighted PCG: I+dt K^2 is self-adjoint positive in W for frozen M=1. */
static double pm_solve(double *u,const double *rhs,double dt){
 int n=pm.n;double *buf=calloc(8*n,sizeof(double));double *r=buf,*z=r+n,*v=z+n,*Av=v+n,*tmp=Av+n,*diag=tmp+n,*Au=diag+n;
 for(int i=0;i<n;i++)diag[i]=1.+dt*sq(pm.diag[i]);
 pm_apply(u,Au,tmp,dt);for(int i=0;i<n;i++){r[i]=rhs[i]-Au[i];z[i]=r[i]/diag[i];v[i]=z[i];}
 double rz=pm_dot(r,z),initial=sqrt(pm_dot(rhs,rhs)),res=sqrt(pm_dot(r,r));int it=0;
 double tolerance=1e-12*fmax(initial,1.);
 while(res>tolerance&&it++<2000){
  pm_apply(v,Av,tmp,dt);double denom=pm_dot(v,Av);if(!(denom>0.)||!isfinite(denom))exit(44);double alpha=rz/denom;
  for(int i=0;i<n;i++){u[i]+=alpha*v[i];r[i]-=alpha*Av[i];z[i]=r[i]/diag[i];}
  double next=pm_dot(r,z),beta=rz>0.?next/rz:0.;for(int i=0;i<n;i++)v[i]=z[i]+beta*v[i];rz=next;res=sqrt(pm_dot(r,r));
 }
 pm_apply(u,Au,tmp,dt);double exact=0.;for(int i=0;i<n;i++)exact=fmax(exact,fabs(Au[i]-rhs[i]));
 if(it>=2000||!isfinite(exact)||exact>1e-9){fprintf(stderr,"PCG_FAIL it=%d residual=%.17g\n",it,exact);exit(45);}free(buf);return exact;
}
static void pm_prefix(const double *z){
 int n=amr_latent_n,st=n+1;memset(amr_prefix,0,(size_t)st*st*sizeof(double));
 for(int j=0;j<n;j++){double row=0.;for(int i=0;i<n;i++){row+=z[j*n+i];amr_prefix[(j+1)*st+i+1]=row+amr_prefix[j*st+i+1];}}
}
static void pm_component_div(int component,double scale,double *out){
 for(int f=0;f<pm.ne;f++){struct PairEdge *e=&pm.edge[f];if(e->c!=component)continue;
  double v=scale*e->length*amr_box(e->sx,e->sy,e->mx,e->my);for(int k=0;k<e->ng;k++)out[e->id[k]]-=e->g[k]*v/pm.w[e->id[k]];}
}
static void pm_noise(int step,double dt,double amplitude,double *out){
 memset(out,0,pm.n*sizeof(double));if(tcw_noise_mode==2)tcw_prepare_noise_local(step,dt);
 int n=amr_latent_n;double *z=calloc(n*n,sizeof(double));
 for(int c=0;c<2;c++){
  for(int j=0;j<n;j++)for(int i=0;i<n;i++)z[j*n+i]=tcw_noise_mode==2?tcw_se_field[c][j*n+i]:tcw_hash_gaussian(i,j,c,step,n,0x243f6a8885a308d3ULL);
  pm_prefix(z);pm_component_div(c,amplitude,out);
 }
 free(z);
}
static uint64_t pm_hidden_hash(void){
 uint64_t v=1469598103934665603ULL;int n=tcw_noise_n;
 for(int c=0;c<2;c++)if(tcw_se_endpoint[c]){unsigned char *p=(unsigned char*)tcw_se_endpoint[c];for(size_t i=0;i<(size_t)n*n*sizeof(double);i++){v^=p[i];v*=1099511628211ULL;}}
 return v;
}
#endif
