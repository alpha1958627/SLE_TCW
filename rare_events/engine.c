/* Shared fixed-grid, checkpointable SLE engine. 2026-10-19.
 * Deterministic solver and conservative periodic noise come from physics.h.
 * A checkpoint contains a source signature and every physical parameter.
 */
#include "grid/multigrid.h"
#include "poisson.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <limits.h>

#include "build/signature.h"
#define LEVEL_MIN 6
#define NN 64
static double param_A, param_phi, tcw_slip;
static uint64_t tcw_noise_seed;
#define COEFF_A param_A
#define COEFF_VARPHI param_phi
#include "physics.h"

scalar h[], h_old[];
static double ratio = 1.01, slip_ref = .2, step_dt = .003;
static double event_h = .01, floor_h = .005, mass0, best_q = 1., clamp_mass;
static uint64_t completed = 0, end_step = 2000;
static int save_records = 0, observe = 1000;

static void fail (const char *message) {
  if (pid() == 0) fprintf(stderr, "ERROR %s\n", message);
#if _MPI
  MPI_Abort(MPI_COMM_WORLD, 2);
#endif
  exit(2);
}

static double *full_field (void) {
  double *buf = calloc(NN*NN, sizeof(double));
  if (!buf) fail("field allocation");
  foreach() {
    int ix = (int)floor((x-X0)/Delta), iy = (int)floor((y-Y0)/Delta);
    if (ix < 0 || ix >= NN || iy < 0 || iy >= NN) fail("field index");
    buf[iy*NN+ix] = h[];
  }
#if _MPI
  MPI_Allreduce(MPI_IN_PLACE, buf, NN*NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return buf;
}

static void checkpoint_write (const char *name) {
  double *field = full_field();
  if (pid() == 0) {
    char temp[512]; snprintf(temp, sizeof(temp), "%s.tmp", name);
    FILE *fp = fopen(temp, "wb"); if (!fp) fail("checkpoint open");
    const char magic[8] = {'R','E','S','T','F','E','2','6'};
    char signature[64]; memset(signature, 0, sizeof(signature));
    strncpy(signature, ENGINE_SIGNATURE, sizeof(signature));
    uint32_t n = NN, version = 1;
    double p[11] = {step_dt,ratio,slip_ref,param_phi,param_A,L0,event_h,floor_h,mass0,best_q,clamp_mass};
    if (fwrite(magic, 8, 1, fp) != 1 || fwrite(signature, 64, 1, fp) != 1 ||
        fwrite(&n, 4, 1, fp) != 1 || fwrite(&version, 4, 1, fp) != 1 ||
        fwrite(&completed, 8, 1, fp) != 1 || fwrite(&tcw_noise_seed, 8, 1, fp) != 1 ||
        fwrite(p, sizeof(double), 11, fp) != 11 || fwrite(field, sizeof(double), NN*NN, fp) != NN*NN)
      fail("checkpoint write");
    if (fclose(fp) || rename(temp, name)) fail("checkpoint commit");
  }
  free(field);
}

static void checkpoint_read (const char *name) {
  FILE *fp = fopen(name, "rb"); if (!fp) fail("restore open");
  char magic[8], signature[64]; uint32_t n, version; uint64_t old_seed;
  double p[11], *field = malloc(NN*NN*sizeof(double));
  if (!field || fread(magic, 8, 1, fp) != 1 || fread(signature, 64, 1, fp) != 1 ||
      fread(&n, 4, 1, fp) != 1 || fread(&version, 4, 1, fp) != 1 ||
      fread(&completed, 8, 1, fp) != 1 || fread(&old_seed, 8, 1, fp) != 1 ||
      fread(p, sizeof(double), 11, fp) != 11 || fread(field, sizeof(double), NN*NN, fp) != NN*NN || fgetc(fp) != EOF)
    fail("restore length");
  fclose(fp);
  if (memcmp(magic, "RESTFE26", 8) || memcmp(signature, ENGINE_SIGNATURE, 64) || n != NN || version != 1)
    fail("restore signature/version mismatch");
  double expected[8] = {step_dt,ratio,slip_ref,param_phi,param_A,L0,event_h,floor_h};
  for (int j = 0; j < 8; j++) if (p[j] != expected[j]) fail("restore physical parameters mismatch");
  if (completed > end_step) fail("restore past horizon");
  foreach() {
    int ix = (int)floor((x-X0)/Delta), iy = (int)floor((y-Y0)/Delta);
    h[] = field[iy*NN+ix];
    if (!isfinite(h[]) || h[] < floor_h) fail("restore field invalid");
  }
  free(field); mass0=p[8]; best_q=p[9]; clamp_mass=p[10];
  boundary({h});
}

static double measure (double *minimum, double *maximum, double *mass) {
  double lo=HUGE_VAL, hi=-HUGE_VAL, vol=0., box=HUGE_VAL;
  foreach(reduction(min:lo) reduction(max:hi) reduction(+:vol) reduction(min:box)) {
    lo = min(lo, h[]); hi = max(hi, h[]); vol += h[]*sq(Delta);
    double average=0.;
    for (int dy=-1; dy<=1; dy++)
      for (int dx=-1; dx<=1; dx++) average += h[dx,dy]/9.;
    box=min(box,average);
  }
  *minimum=lo; *maximum=hi; *mass=vol;
  return lo + max(0.,min(1.,(lo-event_h)/(.2-event_h)))*(box-lo);
}

int main (int argc, char **argv) {
  const char *restore=NULL;
  for (int j=1; j<argc; j++) {
    if (!strcmp(argv[j],"--records")) save_records=1;
    else if (j+1>=argc) fail("option requires value");
    else if (!strcmp(argv[j],"--seed")) tcw_noise_seed=strtoull(argv[++j],NULL,10);
    else if (!strcmp(argv[j],"--end-step")) end_step=strtoull(argv[++j],NULL,10);
    else if (!strcmp(argv[j],"--dt")) step_dt=strtod(argv[++j],NULL);
    else if (!strcmp(argv[j],"--ratio")) ratio=strtod(argv[++j],NULL);
    else if (!strcmp(argv[j],"--slip")) slip_ref=strtod(argv[++j],NULL);
    else if (!strcmp(argv[j],"--observe")) observe=atoi(argv[++j]);
    else if (!strcmp(argv[j],"--restore")) restore=argv[++j];
    else fail("unknown option");
  }
  if (!tcw_noise_seed || !(step_dt>0.) || !(ratio>0.) || slip_ref<0. || end_step>INT_MAX || observe<0)
    fail("invalid configuration");
  if (strlen(ENGINE_SIGNATURE)!=64) fail("build signature missing");
  param_A=.0013159/sq(ratio); param_phi=.01/sq(ratio); tcw_slip=slip_ref/ratio;
  periodic(right); periodic(top); origin(0.,0.); L0=100./ratio;
  init_grid(NN);
  if (restore) checkpoint_read(restore);
  else {
    foreach() h[]=1.;
    boundary({h}); mass0=sq(L0);
  }
  double lo,hi,mass,q=measure(&lo,&hi,&mass),raw_min=lo;
  FILE *trace_fp=NULL,*records=NULL;
  if (pid()==0) {
    trace_fp=fopen("trajectory.tsv","w"); if (!trace_fp) fail("trace open");
    fprintf(trace_fp,"step\tt\thmin\tq3\tmass_relative_drift\traw_min\tclamp_mass\n");
    if (save_records) {records=fopen("records.tsv","w"); if (!records) fail("records open");}
  }
  int hit = lo<=event_h, invalid=0;
  while (completed<end_step && !hit && !invalid) {
    foreach() h_old[]=h[];
    boundary({h_old});
    if (hi-lo>1e-14) solve_thin_film(h,step_dt,3.,3.,5.,(int)completed);
    apply_split_stochastic_flux(h,h_old,step_dt,(int)completed);
    raw_min=HUGE_VAL; double added=0.; int bad=0;
    foreach(reduction(min:raw_min) reduction(+:added) reduction(max:bad)) {
      if (!isfinite(h[])) bad=1;
      raw_min=min(raw_min,h[]);
      if (h[]<floor_h) {added+=(floor_h-h[])*sq(Delta); h[]=floor_h;}
    }
    if (bad) {completed++; checkpoint_write("failed.chk"); fail("nonfinite physical field retained");}
    boundary({h}); completed++; clamp_mass+=added;
    q=measure(&lo,&hi,&mass); hit=lo<=event_h;
    int record=q<best_q;
    if (record) {
      best_q=q;
      if (save_records) {
        char name[80];snprintf(name,sizeof(name),"record_%012llu.chk",(unsigned long long)completed);
        checkpoint_write(name);
        if (pid()==0) {fprintf(records,"%llu\t%.17g\t%.17g\t%s\n",(unsigned long long)completed,completed*step_dt,q,name);fflush(records);}
      }
    }
    if (hi>100. || fabs((mass-mass0)/mass0)>1e-5) invalid=1;
    if (pid()==0 && ((observe && completed%observe==0) || hit || completed==end_step)) {
      fprintf(trace_fp,"%llu\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
              (unsigned long long)completed,completed*step_dt,lo,q,(mass-mass0)/mass0,raw_min,clamp_mass);
      fflush(trace_fp);
    }
  }
  checkpoint_write("final.chk");
  if (pid()==0) {
    if (records) fclose(records);fclose(trace_fp);
    FILE *out=fopen("summary.json","w");if (!out) fail("summary open");
    fprintf(out,"{\"reason\":\"%s\",\"step\":%llu,\"time\":%.17g,\"hmin\":%.17g,\"raw_min\":%.17g,\"q3\":%.17g,\"best_q3\":%.17g,\"mass_relative_drift\":%.17g,\"clamp_mass\":%.17g,\"seed\":%llu}\n",
      invalid?"invalid":hit?"rupture":"horizon",(unsigned long long)completed,completed*step_dt,lo,raw_min,q,best_q,
      (mass-mass0)/mass0,clamp_mass,(unsigned long long)tcw_noise_seed);
    if (fclose(out)) fail("summary close");
  }
  return invalid?3:0;
}
