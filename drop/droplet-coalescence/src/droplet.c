#include "grid/quadtree.h"
#include "poisson.h"
#include <time.h>
#include <float.h>
#include "amr_noise.h"
#include "paired_graph.h"
#include "nonlinear_pair.h"
scalar h[],p[];
static void amr_grid_sync(void){long n=0;int hi=0;foreach(serial){n++;hi=max(hi,level);}grid->tn=n;grid->maxdepth=hi;}
#include "amr_height.h"
static int maximum=8,minimum=6;
static double phi=1e-4,stepdt=.003125,lower_bound=1e-6,activation_threshold=.03;
static void native_output(const char *name){
 FILE *f=fopen(name,"wb");if(!f)exit(90);int cols=6;uint64_t n=pm.n;fwrite(&cols,4,1,f);fwrite(&n,8,1,f);fwrite(&L0,8,1,f);
 foreach(serial){double a[]={x,y,Delta,h[],p[],level};fwrite(a,8,6,f);}if(fclose(f))exit(90);
}
static double mass_vector(const double *v){double a=0.;for(int i=0;i<pm.n;i++)a+=v[i]*pm.w[i];return a;}
static int save_new_activations(FILE *f,int *first,const double *height,int step){
 int added=activation_update(pm.n,first,height,activation_threshold,step);
 if(added)for(int i=0;i<pm.n;i++)if(first[i]==step)
  fprintf(f,"%d,%d,%.17g,%.17g,%.17g,%.17g\n",step,i,pm.x[i],pm.y[i],pm.d[i],height[i]);
 return added;
}
int main(int argc,char **argv){
 int mode=argc>1?atoi(argv[1]):2;uint64_t seed=argc>2?strtoull(argv[2],NULL,10):102901;int uniform=argc>3?atoi(argv[3]):0;
 if(argc>4)stepdt=atof(argv[4]);int steps=argc>5?atoi(argv[5]):3200;
 if(argc>6)maximum=atoi(argv[6]);if(argc>7)phi=atof(argv[7]);if(argc>8)lower_bound=atof(argv[8]);if(argc>9)activation_threshold=atof(argv[9]);if(!isfinite(activation_threshold)||activation_threshold<0.)return 2;if(!isfinite(phi)||phi<0.||!isfinite(lower_bound)||lower_bound<=0.||lower_bound>=.01)return 2;if(maximum<6||maximum>9||steps<1||stepdt<=0.||npe()!=1)return 2;
 if(mode<0)phi=0.;L0=64.;origin(0.,0.);periodic(right);periodic(top);init_grid(1<<(uniform?maximum:minimum));
 h.refine=amr_refine_positive;h.prolongation=refine_linear;h.restriction=restriction_volume_average;
 if(!uniform){
  refine(level<maximum-1&&hypot(x-32.,y-32.)-.5*sqrt(2.)*Delta<=18.);
  refine(level<maximum&&hypot(x-32.,y-32.)-.5*sqrt(2.)*Delta<=14.);
 }
 amr_grid_sync();pm_build(maximum);foreach(){h[]=amr_initial_average(x,y,Delta,maximum);p[]=0.;}boundary({h,p});
 amr_noise_configure(minimum,maximum,mode,.05,.01,seed);nl_geometry();
 FILE *mask_cfg=fopen("activation_config.json","w");if(!mask_cfg)exit(97);
 fprintf(mask_cfg,"{\"threshold\":%.17g,\"comparison\":\">\",\"permanent\":true,\"face_rule\":\"all_nonzero_stencil_cells_active\",\"mesh\":\"fixed\",\"restart_supported\":false}\n",activation_threshold);fclose(mask_cfg);
 FILE *cfg=fopen("runtime_config.json","w");fprintf(cfg,"{\"phi\":%.17g,\"lower_bound\":%.17g,\"sigma\":0.05,\"tau\":0.01,\"dt\":%.17g,\"steps\":%d,\"mode\":%d,\"seed\":%llu,\"hard_floor\":0,\"white_mask\":false}\n",phi,lower_bound,stepdt,steps,mode,(unsigned long long)seed);fclose(cfg);
 double *buffer=calloc(6*pm.n,sizeof(double)),*old=buffer,*u=old+pm.n,*rhs=u+pm.n,*div=rhs+pm.n,*pressure=div+pm.n;pm_read(h,old);
 int *first_active=malloc(pm.n*sizeof(int));if(!first_active)exit(97);
 for(int i=0;i<pm.n;i++)first_active[i]=-1;
 FILE *activation_events=fopen("activation_events.csv","w");if(!activation_events)exit(97);
 fprintf(activation_events,"step,cell,x,y,delta,height\n");
 int active_total=save_new_activations(activation_events,first_active,old,0);
 FILE *activation_metrics=fopen("activation_metrics.csv","w");if(!activation_metrics)exit(97);
 fprintf(activation_metrics,"step,t,threshold,active_before,active_after,allowed_faces,inactive_source_max,weighted_noise_sum\n");
 double initial=mass_vector(old),start=amr_clock(),maxmass=0.,minall=HUGE_VAL,maxsolve=0.,maxrecon=0.,maxguard=0.;int maxiter=0;long firststrict=-1;struct NlStats stats={.first=-1};
 int protected_fail=0;for(int i=0;i<pm.n;i++)if(hypot(pm.x[i]-32.,pm.y[i]-32.)-.5*sqrt(2.)*pm.d[i]<=14.&&pm.d[i]>L0/(1<<maximum)+1e-12)protected_fail++;
 if(protected_fail)exit(91);pm_lap(old,pressure);for(int i=0;i<pm.n;i++)pressure[i]=-pressure[i];pm_write(pressure,p);native_output("initial.native");
 int countlevels[12]={0},actualmin=20,actualmax=0;for(int i=0;i<pm.n;i++){int lev=llround(log2(L0/pm.d[i]));countlevels[lev]++;actualmin=min(actualmin,lev);actualmax=max(actualmax,lev);}
 FILE *mesh=fopen("mesh_summary.json","w");fprintf(mesh,"{\"uniform\":%d,\"cells\":%d,\"uniform_fine_reference_cells\":%d,\"actual_finest_cells\":%d,\"cells_level6\":%d,\"cells_level7\":%d,\"cells_level8\":%d,\"protected_radius\":14,\"buffer_radius\":18,\"actual_min_level\":%d,\"actual_max_level\":%d,\"protected_coverage_pass\":true,\"initial_volume\":%.17g}\n",uniform,pm.n,(1<<maximum)*(1<<maximum),countlevels[maximum],countlevels[6],countlevels[7],countlevels[8],actualmin,actualmax,initial);fclose(mesh);
 FILE *f=fopen("metrics.csv","w");fprintf(f,"step,t,dt,cells,hmin,hmax,volume,mass_relative,energy,solver_iter,solver_residual,first_limiter_step,limiter_hits,first_strict_step,removed_total_flux_RMS,strict_area,near_area,core_hmin,far_hmin,guard_ratio,wall_seconds\n");
 double energy_prev=pm_energy(old);int stride=max(1,(int)llround(.1/stepdt)),native_stride=max(1,(int)llround(.1/stepdt));
 for(int step=0;step<=steps;step++){
  int iterations=0;double residual=0.,inactive_source_max=0.,noise_sum=0.,noise_abs=0.;int active_before=active_total;nl_active_faces=0;
  if(step){
   pm_read(h,old);nl_prepare(old,mode,phi,stepdt);nl_random_flux(step-1,stepdt,mode,first_active);nl_div(nl_noise,div);
   for(int i=0;i<pm.n;i++){
    if(first_active[i]<0)inactive_source_max=fmax(inactive_source_max,fabs(div[i]));
    noise_sum+=pm.w[i]*div[i];noise_abs+=fabs(pm.w[i]*div[i]);
   }
   if(inactive_source_max!=0.||fabs(noise_sum)>1e-10*fmax(noise_abs,1.))exit(98);
   for(int i=0;i<pm.n;i++){rhs[i]=old[i]+stepdt*div[i];u[i]=old[i];}
   iterations=nl_solve(u,rhs,stepdt,&residual);maxiter=max(maxiter,iterations);maxsolve=fmax(maxsolve,residual);nl_total_flux(u,pressure);nl_div(nl_flux,div);
   for(int i=0;i<pm.n;i++)maxrecon=fmax(maxrecon,fabs(old[i]+stepdt*div[i]-u[i]));
   if(maxrecon>1e-8)exit(92);
   nl_limit(old,stepdt,lower_bound,step,&stats);nl_div(nl_flux,div);for(int i=0;i<pm.n;i++)u[i]=old[i]+stepdt*div[i];pm_write(u,h);
   pm_lap(u,pressure);for(int i=0;i<pm.n;i++)pressure[i]=-pressure[i];pm_write(pressure,p);
  }else pm_read(h,u);
  if(step)active_total+=save_new_activations(activation_events,first_active,u,step);
  if(step%stride==0||step==steps){
   fprintf(activation_metrics,"%d,%.17g,%.17g,%d,%d,%d,%.17g,%.17g\n",step,step*stepdt,activation_threshold,active_before,active_total,nl_active_faces,inactive_source_max,noise_sum);
   fflush(activation_metrics);fflush(activation_events);
  }
  double lo=HUGE_VAL,hi=-HUGE_VAL,core=HUGE_VAL,far=HUGE_VAL,area=0.,near=0.,guard=0.,guardarea=0.;
  for(int i=0;i<pm.n;i++){
   if(!isfinite(u[i])||u[i]<=0.||u[i]<lower_bound-64*DBL_EPSILON*.01)exit(93);
   lo=fmin(lo,u[i]);hi=fmax(hi,u[i]);double r=hypot(pm.x[i]-32.,pm.y[i]-32.);
   if(r<=14.)core=fmin(core,u[i]);else far=fmin(far,u[i]);
   if(u[i]<=lower_bound){area+=pm.w[i];if(firststrict<0)firststrict=step;}if(u[i]<=lower_bound+64*DBL_EPSILON*.01)near+=pm.w[i];
   if(r>=12.&&r<13.){guard+=(u[i]-.01)*pm.w[i];guardarea+=pm.w[i];}
  }
  double vol=mass_vector(u),mass=fabs((vol-initial)/initial),energy=pm_energy(u),ratio=guardarea>0.?guard/guardarea/fmax(hi-.01,1e-12):0.;maxguard=fmax(maxguard,ratio);maxmass=fmax(maxmass,mass);minall=fmin(minall,lo);
  if(mass>1e-9||ratio>.05)exit(94);
  if(mode<0&&stats.first<0&&energy>energy_prev*(1.+1e-8)+1e-10)exit(95);energy_prev=energy;
  if(step%stride==0||step==steps){fprintf(f,"%d,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%ld,%ld,%ld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",step,step*stepdt,stepdt,pm.n,lo,hi,vol,(vol-initial)/initial,energy,iterations,residual,stats.first,stats.hits,firststrict,stats.raw2>0.?sqrt(stats.removed2/stats.raw2):0.,area,near,core,far,ratio,amr_clock()-start);fflush(f);}
  if(step&&(step==1||step==steps||step==stats.first||step==firststrict||step%native_stride==0||fabs(step*stepdt-.1)<1e-12||fabs(step*stepdt-.5)<1e-12)){
   char name[80];sprintf(name,"field_step%09d.native",step);native_output(name);
  }
 }
 fclose(f);fclose(activation_events);fclose(activation_metrics);
 f=fopen("activation_final.bin","wb");if(!f)exit(97);int active_n=pm.n;
 if(fwrite(&active_n,sizeof(int),1,f)!=1||fwrite(first_active,sizeof(int),pm.n,f)!=(size_t)pm.n||fclose(f))exit(97);
 native_output("final.native");f=fopen("metadata.json","w");
 fprintf(f,"{\"mode\":%d,\"seed\":%llu,\"uniform\":%d,\"L\":64,\"R0\":5,\"hp\":0.01,\"excess_peak\":1,\"A\":0,\"phi\":%.17g,\"lower_bound\":%.17g,\"sigma\":0.05,\"tau\":0.01,\"dt\":%.17g,\"steps\":%d,\"T\":%.17g,\"cells\":%d,\"maximum\":%d,\"initial_volume\":%.17g,\"minimum_h\":%.17g,\"max_mass_relative\":%.17g,\"first_limiter_step\":%ld,\"limiter_hits\":%ld,\"first_strict_step\":%ld,\"removed_total_flux_RMS\":%.17g,\"max_solver_residual\":%.17g,\"max_solver_iterations\":%d,\"max_reconstruction_error\":%.17g,\"guard_ratio_max\":%.17g,\"hard_floor_added\":0,\"white_mask\":false,\"ranks\":1,\"MPI\":false,\"wall_seconds\":%.17g}\n",mode,(unsigned long long)seed,uniform,phi,lower_bound,stepdt,steps,steps*stepdt,pm.n,maximum,initial,minall,maxmass,stats.first,stats.hits,firststrict,stats.raw2>0.?sqrt(stats.removed2/stats.raw2):0.,maxsolve,maxiter,maxrecon,maxguard,amr_clock()-start);fclose(f);
 printf("DROPLET_COMPLETE mode=%d T=%.17g min=%.17g mass=%.17g limiter=%ld\n",mode,steps*stepdt,minall,maxmass,stats.first);free(first_active);free(buffer);nl_clear();amr_noise_clear();pm_clear();return 0;
}
