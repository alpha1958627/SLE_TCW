// v4: y-index-corrected wall MG with consistent multiplicative noise flux.
#include "poisson.h"
#include "run.h"

#include "grid/tree.h"



#include <sys/stat.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>

face vector prescribed_noise[];
double noise_phi=0.;
int outer_count, mg_cycles;
double outer_error, mg_res;
struct Layer {
    scalar f;
    double dt;
    double m;
    double n;
    double dd;
    face vector D;
    face vector Noise1;
    face vector beta;
    face vector gamma;
    scalar lambda;
    double tolerance;
  };
  
  static double residual_layer (scalar * al, scalar * bl,
                    scalar * resl, void * data)
  {
    scalar h = al[0], b = bl[0], res = resl[0];
    struct Layer * p = (struct Layer *) data;
    (const) face vector alpha = p->D;
    (const) face vector Noise1 = p->Noise1;
    (const) face vector beta = p->beta;
    (const) face vector gamma = p->gamma;
    (const) scalar lambda = p->lambda;
    double mytolerance = p->tolerance;
    double maxres = 0.;
    
    face vector g[], hf[] ,No[];
  
    foreach_face()
      hf.x[] = face_gradient_x (h, 0);
    boundary_flux ({hf});
  
    scalar pres[];
    foreach() {
      
      pres[] = 0.;
      foreach_dimension()
        pres[] -= (hf.x[1]- hf.x[0])/Delta;
    }
    boundary({pres});
    
    foreach_face() {
      
      g.x[] = alpha.x[]*face_gradient_x (pres, 0);
      hf.x[] = beta.x[]*face_gradient_x (h, 0) +
        gamma.x[]*(h[]+h[-1])/2.;
      No.x[]=Noise1.x[]*(h[]+h[-1])/2.; // v4: same face-height factor as relax_layer
    }
    boundary_flux ({g, hf,No});
    foreach (reduction(max:maxres)) {
      res[] = b[] - lambda[]*h[];
      foreach_dimension()
          res[] -= (g.x[1] - g.x[])/Delta + (hf.x[1] - hf.x[])/Delta +(No.x[1] - No.x[])/Delta;
      if (fabs (res[]) > maxres)
        maxres = fabs (res[]);
    }
    boundary (resl);
    return maxres;
  }
  static void relax_layer (scalar * al, scalar * bl, int l, void * data)
  {
    scalar h = al[0], b = bl[0];
    struct Layer * p = (struct Layer *) data;
    (const) face vector alpha = p->D;
    (const) face vector Noise1 = p->Noise1;
    (const) face vector beta = p->beta;
    (const) face vector gamma = p->gamma;
    (const) scalar lambda = p->lambda;
  
    
  #if JACOBI
    scalar c[];
  #else
    scalar c = h;
  #endif
  
    foreach_level_or_leaf (l) {
      double n =  sq(sq(Delta))*b[], d =  lambda[]*sq(sq(Delta));
      d -= 5*(alpha.x[] + alpha.x[1] + alpha.y[] + alpha.y[0,1]);
      d += 0.5*Delta*sq(Delta)*(gamma.x[1]-gamma.x[] + gamma.y[0,1]-gamma.y[]);
      d += 0.5*Delta*sq(Delta)*(Noise1.x[1]-Noise1.x[] + Noise1.y[0,1]-Noise1.y[]);
      foreach_dimension() {
        n -= (beta.x[1]*h[1] + beta.x[]*h[-1])*sq(Delta);
        d -= (beta.x[1] + beta.x[])*sq(Delta);
      }
  
      n -= alpha.x[]*(-h[-2] - h[-1,-1] + 5*h[-1] - h[-1,1] + h[0,-1] + h[0,1] + h[1]) +
        alpha.x[1]*(h[-1] + h[0,-1] + h[0,1] - h[1,-1] + 5*h[1] - h[1,1] - h[2]) +
        alpha.y[]*(-h[-1,-1] + h[-1] - h[0,-2] + 5*h[0,-1] + h[0,1] - h[1,-1] + h[1]) +
        alpha.y[0,1]*(h[-1] - h[-1,1] + h[0,-1] + 5*h[0,1] - h[0,2] + h[1] - h[1,1]);
  
      n -= 0.5*Delta*sq(Delta)*(h[1]*gamma.x[1] - h[-1]*gamma.x[] +
                    h[0,1]*gamma.y[0,1] - h[0,-1]*gamma.y[]); 
      n -= 0.5*Delta*sq(Delta)*(h[1]*Noise1.x[1] - h[-1]*Noise1.x[] +
                    h[0,1]*Noise1.y[0,1] - h[0,-1]*Noise1.y[]); 
      c[] = n/d;
    }
    
  #if JACOBI
    foreach_level_or_leaf (l)
      h[] = (h[] + 2.*c[])/3.;
  #endif
    
  #if TRASH
    scalar a1[];
    foreach_level_or_leaf (l)
      a1[] = h[];
    trash ({h});
    foreach_level_or_leaf (l)
      h[] = a1[];
  #endif
  }
  trace
  mgstats layer (scalar myf,double mydt,double mym,double myn,double mydd, int timestep){
  struct Layer p = {0};
  p.tolerance = 1e-3;
  mg_cycles = 0;
  p.dt=mydt;
  p.f=myf;
  p.m=mym;
  p.n=myn;
  p.dd=mydd;
  
    if (p.dt == 0.) {
      mgstats s = {0};
      return s;
    }
  
    scalar f = p.f;
    scalar b[], fold[], fa[];
    const scalar lambda[] = - 1./p.dt;
    double m = p.m, n = p.n ,dd=p.dd;
    face vector D[], beta[], gamma[] ,Noise1[];
      foreach() {
      fold[] = f[];
      }
  
    double maxres;
    mgstats s;
    int niter = 0, nitermax = 10;
  
    double defaultol = TOLERANCE;
    if (p.tolerance>1e-10)
      TOLERANCE = p.tolerance;
  
    // 在循环外部计算Noise场
  
//   int NN =256;
//   double Noise8[NN][NN];
//   double dl1=100.0/NN;
//   for (int iy = 0; iy < NN; iy++) {
//       // double yy = (iy + 0.5) * dl1;
//       for (int ix = 0; ix < NN; ix++) {
//           // double xx = (ix + 0.5) * dl1;
//           double del=100.0/NN;
//           Noise8[ix][iy] = noise()/(sqrt(del*p.dt*del));
//       }
//   }
  
    do {
        foreach_face()
            beta.x[] = face_gradient_x(f, 0);
        boundary_flux({beta});

        scalar pres[];
        foreach() {
            pres[] = 0.;
            foreach_dimension()
                pres[] -= (beta.x[1] - beta.x[0]) / Delta;
        }
        boundary({pres});
        
        face vector g[];
        face vector NRight[];
        
        foreach_face() {
            double aa = 0.2;
            double varphi = noise_phi;
            double  hstar = pow(0.1, dd - 3.);
            double ff = (f[] + f[-1]) / 2.;
          
            double Noiseff = prescribed_noise.x[];

            D.x[] = pow(ff, n);
            Noise1.x[] = 1.5 * Noiseff * sqrt(2 * varphi) * pow(ff, 0.5);
            beta.x[] = -3 * aa * pow(ff, n - m - 1.) + dd * aa * hstar * pow(ff, n - dd - 1.);
            gamma.x[] = n * pow(ff, n - 1.) * (pres[] - pres[-1]) / Delta -
                        3 * aa * (n - m - 1.) * pow(ff, n - m - 2.) * (f[] - f[-1]) / Delta +
                        dd * aa * hstar * (n - dd - 1.) * pow(ff, n - dd - 2.) * (f[] - f[-1]) / Delta;
            g.x[] = gamma.x[] * ff;
            NRight.x[] = 0.5 * Noiseff * sqrt(2 * varphi) * pow(ff, 1.5);
        }
        boundary_flux({g});

        foreach() {
            b[] = -fold[] / p.dt;
            fa[] = f[];
            foreach_dimension()
                b[] += (g.x[1] - g.x[]) / Delta + 1 * (NRight.x[1] - NRight.x[]) / Delta;
        }

        restriction({beta, D, lambda, gamma, Noise1});
        p.lambda = lambda;
        p.D = D;
        p.Noise1 = Noise1;
        p.beta = beta;
        p.gamma = gamma;
        s = mg_solve({f}, {b}, residual_layer, relax_layer, &p);
        mg_cycles += s.i;

        maxres = 0;
        foreach(reduction(max:maxres)) {
            double res = fabs(f[] - fa[]);
            if(res > maxres)
                maxres = res;
        }
        boundary({fa});
        niter++;
    } while(maxres > 1e-8 && niter < nitermax);

    if(p.tolerance > 1e-10)
        TOLERANCE = defaultol;

    outer_count=niter; outer_error=maxres; mg_res=s.resa;
    return s;
}



#include <stdint.h>
#include <time.h>
#include <omp.h>
static double coarse[2][64][64];
static uint64_t mix(uint64_t x) {
  x+=UINT64_C(0x9e3779b97f4a7c15); x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);
  x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb); return x^(x>>31);
}
static double uniform(uint64_t key) {return ((mix(key)>>11)+.5)*0x1p-53;}
static double gaussian(uint64_t key) {return sqrt(-2*log(uniform(2*key)))*cos(2*pi*uniform(2*key+1));}
static double interp(int axis,double xx,double yy) {
  double gx=xx/(100./64)-(axis==0 ? 0. : .5),gy=yy/(100./64)-(axis==1 ? 0. : .5);
  int ix=(int)floor(gx),iy=(int)floor(gy); double fx=gx-ix,fy=gy-iy;
  int a=(ix%64+64)%64,b=(iy%64+64)%64,c=(a+1)%64,d=(b+1)%64;
  return (1-fx)*(1-fy)*coarse[axis][a][b]+fx*(1-fy)*coarse[axis][c][b]+(1-fx)*fy*coarse[axis][a][d]+fx*fy*coarse[axis][c][d];
}
scalar height[], previous_height[];
static double equation_error(double dt) {
  scalar pressure[]; face vector slope[],flux[];
  foreach_face() slope.x[]=face_gradient_x(height,0);
  boundary_flux({slope});
  foreach(){pressure[]=0.;foreach_dimension() pressure[]-=(slope.x[1]-slope.x[])/Delta;}
  boundary({pressure});
  foreach_face(){double hf=(height[]+height[-1])/2.;
    flux.x[]=pow(hf,3)*face_gradient_x(pressure,0)+(-.6/hf+.01/pow(hf,3))*face_gradient_x(height,0)+sqrt(2*noise_phi)*pow(hf,1.5)*prescribed_noise.x[];
  }
  boundary_flux({flux});double maximum=0.;
  foreach(reduction(max:maximum)){double r=height[]-previous_height[];foreach_dimension()r-=dt*(flux.x[1]-flux.x[])/Delta;if(fabs(r)>maximum)maximum=fabs(r);}
  return maximum;
}
int main(int argc,char **argv) {
  if(argc!=5 && argc!=6 && argc!=8) return 80;
  noise_phi=atof(argv[1]);int steps=atoi(argv[2]),warm=atoi(argv[3]);char *out=argv[4];
  double dt_input=argc>=6 ? atof(argv[5]) : .001;
  int base_offset=argc==8 ? atoi(argv[7]) : 0;
  if(base_offset<0)return 80;
  int factor=(int)lround(dt_input/.001);
  if(noise_phi<0 || steps<1 || warm<0 || factor<1 || factor>300 || fabs(dt_input-factor*.001)>1e-12)return 80;
  periodic(right);periodic(top);origin(0,0);size(100);init_grid(1<<7);
  foreach() {
    int ix=(int)floor(x/Delta),iy=(int)floor(y/Delta);
    height[]=1.+.01*(2*uniform(7781+ix*128+iy)-1);
  }
  if(argc==8) {
    FILE *in=fopen(argv[6],"rb");if(!in)return 83;
    foreach(serial) {
      double row[3];if(fread(row,sizeof(double),3,in)!=3)exit(84);
      if(fabs(row[0]-x)>1e-12 || fabs(row[1]-y)>1e-12 || !isfinite(row[2]) || row[2]<=0)exit(85);
      height[]=row[2];
    }
    if(fgetc(in)!=EOF)return 86;fclose(in);
  }
  boundary({height}); double v0=statsf(height).sum;
  FILE *f=fopen("history.csv","w");if(!f)return 81;
  fprintf(f,"step,time,dt,wall_seconds,hmin,hmax,relative_mass,outer_iterations,outer_update,mg_cycles,mg_residual,equation_residual\n");
  double total=0.,noise_sum=0.,lastmin=0.,lastmax=0.,drift=0.;int totalcycles=0,totalouter=0;
  for(int step=0;step<warm+steps;step++) {
    double start=omp_get_wtime(),dt=dt_input;
    for(int a=0;a<2;a++)for(int ix=0;ix<64;ix++)for(int iy=0;iy<64;iy++) {
      double sum=0.;
      for(int sub=0;sub<factor;sub++) sum+=gaussian(UINT64_C(61101)+((uint64_t)base_offset+(uint64_t)step*factor+sub)*8192+a*4096+ix*64+iy);
      coarse[a][ix][iy]=sum/sqrt((double)factor)/((100./64)*sqrt(dt));
    }
    noise_sum+=(coarse[0][0][0]+coarse[1][63][63])*dt; // shared-path integrated checksum
    foreach_face(x) prescribed_noise.x[]=interp(0,x,y);
    foreach_face(y) prescribed_noise.y[]=interp(1,x,y);
    boundary((scalar*){prescribed_noise});
    foreach() previous_height[]=height[];
    layer(height,dt,3.,3.,5.,step);
    double elapsed=omp_get_wtime()-start;
    double physical_residual=equation_error(dt);
    if(step>=warm){total+=elapsed;totalcycles+=mg_cycles;totalouter+=outer_count;}
    stats st=statsf(height);lastmin=st.min;lastmax=st.max;drift=(st.sum-v0)/v0;
    fprintf(f,"%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%d,%.17g,%.17g\n",step+1,base_offset*.001+(step+1)*dt,dt,elapsed,lastmin,lastmax,drift,outer_count,outer_error,mg_cycles,mg_res,physical_residual);
    if(!isfinite(lastmin)||!isfinite(lastmax)||lastmin<=0||fabs(drift)>.01){fclose(f);return 82;}
  }
  fclose(f);
  FILE *field=fopen("final.bin","wb");if(!field)return 81;
  foreach(serial){double row[]={x,y,height[]};fwrite(row,sizeof(double),3,field);}fclose(field);
  FILE *fp=fopen(out,"w");if(!fp)return 81;
  fprintf(fp,"{\"phi\":%.17g,\"threads\":%d,\"grid\":128,\"noise_grid\":64,\"steps\":%d,\"warmup\":%d,\"dt\":%.17g,\"seconds\":%.17g,\"seconds_per_step\":%.17g,\"t_star_per_hour\":%.17g,\"hmin\":%.17g,\"hmax\":%.17g,\"relative_mass\":%.17g,\"mg_cycles\":%d,\"outer_iterations\":%d,\"noise_fingerprint\":%.17g}\n",noise_phi,omp_get_max_threads(),steps,warm,dt_input,total,total/steps,steps*dt_input*3600/total,lastmin,lastmax,drift,totalcycles,totalouter,noise_sum);
  fclose(fp);return 0;
}
