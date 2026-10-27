/**
 * @file droplet_spreading.c
 * @brief Cos^4 thin-film droplet spreading with deterministic, white-noise,
 *        and OU-SE stochastic-flux modes.
 * @version 1021.13-test (same coupled update; Brownian-coupled temporal control)
 *
 * Default mode: deterministic. Optional mode 0 is Gaussian white noise with a
 * precursor-film cutoff; mode 2 is periodic OU-SE. OU-SE uses the Eq. (52)
 * amplitude correction unless legacy scaling is explicitly requested.
 * No hard height floor is applied by default.
 */

 #include "poisson.h"
 #include "run.h"
 #include "view.h"
 #include "grid/tree.h"
 #include <sys/stat.h>
 #include <string.h>
 #include <stdbool.h>
 #include <stdlib.h>
 #include <limits.h>
 #include <stdio.h>
 #include <math.h>
 #include <time.h>
#include <stdint.h>
#include <complex.h>
#if _OPENMP
#include <omp.h>
#endif

 /* Global parameters and compile-time settings. */

 #define T2OR 0          // 0: first-order, 1: second-order time integration
#ifndef TCW_MG_MINLEVEL
#define TCW_MG_MINLEVEL 3 /* fourth-order +/-2 stencil needs a nondegenerate periodic coarse grid */
#endif
#define TCW_FROZEN_MOBILITY 1 /* first-order lagged h^3; same pure-capillary PDE, linear implicit step */
#define TCW_POSITIVE_STOCHASTIC 1 /* conservative flux limit; white floor hp, OU floor 1e-6 */
#ifndef TCW_COUPLED_STOCHASTIC
#define TCW_COUPLED_STOCHASTIC 1 /* one raw-noise + implicit-capillary conservative step */
#endif
#define TCW_OUSE_BIN_AVERAGES 1 /* exact OU time-bin integrals; lateral cell-average SE kernel */
static int tcw_bridge_depth=0; /* 0: original counters; 1/2: same Brownian path at dt/2 or dt/4 */

 // Physical coefficients.
 #define COEFF_A 0.0 /* 2026-10-16 pure-capillary verification and free spreading */
 #ifndef COEFF_VARPHI
 #define COEFF_VARPHI 0.03
 #endif
#define TCW_DIAGNOSTIC_NO_DETERMINISTIC 0

// Optional modes: -1 = deterministic, 0 = Gaussian white, 2 = periodic OU-SE.
// Default to no fluctuations; select a stochastic mode explicitly at run time.
#define TCW_NOISE_MODE_DEFAULT -1
#define TCW_NOISE_LEVEL 6

/* 2026-09-03 OU-SE (mode 2) switches */
static int tcw_ouse_legacy_scale = 0; /* 0: Eq.(52) ratio; 1: 20260902 legacy */
static int tcw_ouse_level = 0;       /* 0: noise grid = LEVEL_MIN; else 1<<level */
static double tcw_ouse_amp_multiplier = 1.; /* 2D face-flux calibration only */
static double tcw_varphi_scale = 1.;        /* Fig.8 noise-strength sweep */
static double tcw_r_rms = 0.0;       /* rms(h - <h>) interface roughness */
static double tcw_mg_tolerance = 1e-8;
static double tcw_picard_tolerance = 1e-8;
static int tcw_picard_max = 30;
/* argv[19]: cap for Basilisk's global NITERMAX (poisson.h default 100).
   vs=10 sweeps can need more MG cycles per Picard solve; raising the cap
   keeps mg_tol unchanged instead of loosening the tolerance. */
static int tcw_mg_nitermax = 100;
static double tcw_height_floor = 0.;
static long tcw_solver_calls = 0, tcw_solver_failures = 0;
static int tcw_solver_max_picard = 0, tcw_solver_max_mg_cycles = 0;
static double tcw_solver_max_update = 0., tcw_solver_max_residual = 0.;
static int tcw_run_failed = 0;
static double tcw_Lambda (double xi);
static double tcw_Upsilon (double xi);
static double tcw_Theta (double xi);
static double tcw_grid_delta (void);
static int tcw_ouse_grid_level (void);

 // Output directories.
 char dirFlowFields[] = "dat/";
 char dirFigures[]    = "png/";

 // Global physical parameters.

 double n_exp = 3.0;     // Mobility exponent n.
 double m_exp = 3.0;     // Attractive exponent m.
 double dd_exp = 5.0;    // Repulsive exponent d.

 // Time-step parameters.
 double DTMAX = 5e-5;
 double DTMIN = 5e-5;
double tcw_dt = 1e-3;
int tcw_test_steps = 3;
int tcw_stop_index = 2;

/* 2026-09-15: some MPI runtimes (hpcx on xh) unmask SSE/x87 FP exceptions, so
   an inf/NaN produced mid-step kills the process with SIGFPE before the
   numeric_guard event can stop the run cleanly.  Mask all FP exceptions so
   non-finite values propagate to the guard exactly as on SIPF.  Done both as
   a constructor and once inside the first timestep event: the MPI runtime may
   (re)unmask exceptions during MPI_Init, after process constructors.
   NB: do NOT #include <fenv.h> here - qcc re-emits glibc inline decls and
   chokes; declare fedisableexcept by hand (glibc x86-64, FE_ALL_EXCEPT=0x3f). */
extern int fedisableexcept (int __excepts);
#define TCW_FP_ALL_EXCEPT 0x3f /* FE_INVALID|DENORM|DIVBYZERO|OVER|UNDER|INEXACT */

static void tcw_mask_fp_exceptions (void) {
    fedisableexcept (TCW_FP_ALL_EXCEPT);
}

__attribute__((constructor)) static void tcw_mask_fp_exceptions_ctor (void) {
    tcw_mask_fp_exceptions ();
}

event tcw_fp_mask (i = 0) {
    tcw_mask_fp_exceptions ();
    return 0;
}
 double dtold;
 double errorstep = 0.;

 // Volume statistics.
double V_total  = 0.0;   // ∫ h dA
double V_excess = 0.0;   // Integral of max(h - HOUT, 0) over the domain.
double HOUT = 0.1;
static uint64_t tcw_noise_seed = 1;       // Counter-based noise seed.


// 0: flat TCW; 1: cos^4 droplet, matching the referenced geometry.
static int jcp_case = 1;
static double jcp_initial_volume = 0.;
static double campaign_initial_excess = 0., campaign_domain = 32.;
static double campaign_floor_volume = 0.;
static double campaign_pre_floor_minimum = HUGE_VAL;
static unsigned long long campaign_floor_cells = 0;
static struct timespec campaign_started;
static double opt_noise_seconds=0., opt_solve_seconds=0.;
static double opt_clock (void) {
    struct timespec now; clock_gettime(CLOCK_MONOTONIC,&now);
    return now.tv_sec+1e-9*now.tv_nsec;
}
static int scan_step=-1, scan_first[3]={-1,-1,-1};
static double scan_removed2[3]={0.,0.,0.}, scan_minweight[3]={1.,1.,1.};
static int dt_coupling_ratio=1;
static double *dt_sum[2]={NULL,NULL};
static int campaign_log_stride = 1;
static double campaign_flux_raw2[3]={0.,0.,0.}, campaign_flux_limited2[3]={0.,0.,0.};
static double campaign_flux_active[3]={0.,0.,0.}, campaign_flux_hits[3]={0.,0.,0.};
static double campaign_joint_raw_minimum=HUGE_VAL;
static void campaign_limit_flux (scalar base, face vector flux, double dt_step, int kind);
static const double jcp_hp = .01, jcp_radius = 5., jcp_peak = 1.;
static int tcw_noise_mode = TCW_NOISE_MODE_DEFAULT;
static double tcw_tau = 0.;
static double tcw_sigma = 0.;
static int tcw_fixed_dt = 1;
static int tcw_sample_every = 10000;
static int tcw_noise_snapshot_every = 10000;
static int tcw_noise_probe_every = 10;
static int tcw_noise_snapshot_index = 0;
static int tcw_noise_initialized = 0;
static int tcw_noise_n = 0;
static double *tcw_se_field[2] = {NULL, NULL};
static double *tcw_se_endpoint[2] = {NULL, NULL};
static double *tcw_se_first = NULL;
static double *tcw_se_innovation = NULL;
static double *tcw_se_filter = NULL;
static double tcw_se_norm = 0.;
static double complex *tcw_se_fft = NULL;
static double complex *tcw_se_col = NULL;
static int tcw_ouse_bin_selftest (double sigma);

 // Grid and solver state.
 scalar h[], hi[], l[], h_noise_old[], tcw_stochastic_rhs[];
 face vector tcw_unlimited_stochastic_flux[];
 static int tcw_stochastic_rhs_active = 0;
 #if T2OR
 scalar h2[], h1[], err[];
 #endif

 double hmin, hmin1;
 int lm = 8;
 int levelm;
 int LEVEL_MAX = 8;
 int LEVEL_MIN = 8;
 /* Output management. */

 /**
  * @brief Create a directory if it does not already exist.
  */

void write_hfield_binary(const char *path) {
    int n = 1 << LEVEL_MIN;
    double *buf = calloc((size_t)n * (size_t)n, sizeof(double));
    foreach() {
        int ix = (int)floor((x - X0)/Delta);
        int iy = (int)floor((y - Y0)/Delta);
        if (ix >= 0 && ix < n && iy >= 0 && iy < n)
            buf[(size_t)iy * (size_t)n + (size_t)ix] = h[];
    }
#if _MPI
    double *sum = NULL;
    if (pid() == 0) sum = calloc((size_t)n * (size_t)n, sizeof(double));
    MPI_Reduce(buf, sum, n * n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    free(buf); buf = sum;
#endif
    if (pid() == 0) {
        FILE *fp = fopen(path, "wb");
        fwrite(&n, sizeof(int), 1, fp);
        fwrite(buf, sizeof(double), (size_t)n * (size_t)n, fp);
        fclose(fp);
    }
    if (buf) free(buf);
}

void create_directory(const char *path) {
    if (pid() == 0) {
        char temp_path[256];
        strncpy(temp_path, path, 255);
        temp_path[255] = '\0';

     // Ensure the path ends with a slash for consistent handling.
        size_t len = strlen(temp_path);
        if (len > 0 && temp_path[len - 1] != '/')
            strncat(temp_path, "/", 255 - len);

        struct stat st = {0};
        if (stat(temp_path, &st) == -1 && mkdir(temp_path, 0700) != 0)
            fprintf(stderr, "Warning: Could not create directory %s\n", temp_path);
    }
#if _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

 /**
  * @brief Save scalar fields in binary form.
  *
  * @param list Scalars to save (for example, {h}).
  * @param fname Output path.
  */
 void save_binary_data(scalar *list, const char *fname) {
     FILE *fh = fopen(fname, "w");
     if (!fh) {
         fprintf(stderr, "Error: Cannot open file %s for writing.\n", fname);
         return;
     }

     for (scalar f in list) {
         boundary({f});

 #if _MPI
         // MPI-parallel output.
         scalar index = {-1};
         index = new scalar;
         z_indexing(index, true);

         // Compute the record size for each cell.
         #if dimension == 2
             const int len = 5 * sizeof(double); // level, x, y, Delta, val
         #elif dimension == 3
             const int len = 6 * sizeof(double); // level, x, y, z, Delta, val
         #endif

         foreach() {
             if (is_local(cell)) {
                 long offset = index[] * len;
                 fseek(fh, offset, SEEK_SET);

                 double lvl = (double)level;
                 double val = f[];

                 fwrite(&lvl, sizeof(double), 1, fh);
                 fwrite(&x, sizeof(double), 1, fh);
                 fwrite(&y, sizeof(double), 1, fh);
                 #if dimension == 3
                 fwrite(&z, sizeof(double), 1, fh);
                 #endif
                 fwrite(&Delta, sizeof(double), 1, fh);
                 fwrite(&val, sizeof(double), 1, fh);
             }
         }
         delete({index});
 #else
         // Serial output.
         foreach() {
             double lvl = (double)level;
             double val = f[];

             fwrite(&lvl, sizeof(double), 1, fh);
             fwrite(&x, sizeof(double), 1, fh);
             fwrite(&y, sizeof(double), 1, fh);
             #if dimension == 3
             fwrite(&z, sizeof(double), 1, fh);
             #endif
             fwrite(&Delta, sizeof(double), 1, fh);
             fwrite(&val, sizeof(double), 1, fh);
         }
 #endif
     }
     fclose(fh);
 }

 /* Solver structures and core routines. */

 struct EquationParams {
     scalar f;
     double dt;
     double m;
     double n;
     double dd;
     face vector D;
     face vector beta;
     face vector gamma;
     scalar lambda;
     double tolerance;
 };

 /**
  * @brief Compute the coupled residual.
  */
 static double residual_thin_film(scalar *al, scalar *bl, scalar *resl, void *data) {
     scalar h_sol = al[0], b_rhs = bl[0], res = resl[0];
     struct EquationParams *p = (struct EquationParams *)data;

     (const) face vector alpha = p->D;
     (const) face vector beta = p->beta;
     (const) face vector gamma = p->gamma;
     (const) scalar lambda = p->lambda;

     double maxres = 0.;
     int invalid = 0;
     face vector g[], hf[];

     foreach_face()
         hf.x[] = face_gradient_x(h_sol, 0);

     boundary((scalar *){hf});

     scalar pres[];
     foreach() {
         pres[] = 0.;
         foreach_dimension()
             pres[] -= (hf.x[1] - hf.x[0]) / Delta;
     }
     boundary({pres});

     foreach_face() {
         g.x[] = alpha.x[] * face_gradient_x(pres, 0);
         hf.x[] = beta.x[] * face_gradient_x(h_sol, 0) +
                  gamma.x[] * (h_sol[] + h_sol[-1]) / 2.;
     }
     // Apply boundary conditions to the face flux.
     boundary((scalar *){g, hf});

     foreach(reduction(max:maxres) reduction(max:invalid)) {
         res[] = b_rhs[] - lambda[] * h_sol[];
         foreach_dimension()
             res[] -= (g.x[1] - g.x[]) / Delta +
                      (hf.x[1] - hf.x[]) / Delta ;

         if (!isfinite(res[])) invalid = 1;
         else if (fabs(res[]) > maxres)
             maxres = fabs(res[]);
     }
     boundary(resl);
     return invalid ? HUGE_VAL : maxres;
 }

 /**
  * @brief Relaxation/smoothing function.
  */
 static void relax_thin_film(scalar *al, scalar *bl, int l, void *data) {
     scalar h_sol = al[0], b_rhs = bl[0];
     struct EquationParams *p = (struct EquationParams *)data;

     (const) face vector alpha = p->D;
     (const) face vector beta = p->beta;
     (const) face vector gamma = p->gamma;
     (const) scalar lambda = p->lambda;

 #if JACOBI
     scalar c[];
 #else
     scalar c = h_sol;
 #endif

     foreach_level_or_leaf(l) {
         double num = sq(sq(Delta)) * b_rhs[];
         double den = lambda[] * sq(sq(Delta));

         den -= 5 * (alpha.x[] + alpha.x[1] + alpha.y[] + alpha.y[0,1]);
         den += 0.5 * Delta * sq(Delta) * (gamma.x[1] - gamma.x[] + gamma.y[0,1] - gamma.y[]);

         foreach_dimension() {
             num -= (beta.x[1] * h_sol[1] + beta.x[] * h_sol[-1]) * sq(Delta);
             den -= (beta.x[1] + beta.x[]) * sq(Delta);
         }

         num -= alpha.x[] * (-h_sol[-2] - h_sol[-1,-1] + 5*h_sol[-1] - h_sol[-1,1] + h_sol[0,-1] + h_sol[0,1] + h_sol[1]) +
                alpha.x[1] * (h_sol[-1] + h_sol[0,-1] + h_sol[0,1] - h_sol[1,-1] + 5*h_sol[1] - h_sol[1,1] - h_sol[2]) +
                alpha.y[] * (-h_sol[-1,-1] + h_sol[-1] - h_sol[0,-2] + 5*h_sol[0,-1] + h_sol[0,1] - h_sol[1,-1] + h_sol[1]) +
                alpha.y[0,1] * (h_sol[-1] - h_sol[-1,1] + h_sol[0,-1] + 5*h_sol[0,1] - h_sol[0,2] + h_sol[1] - h_sol[1,1]);

         num -= 0.5 * Delta * sq(Delta) * (h_sol[1]*gamma.x[1] - h_sol[-1]*gamma.x[] +
                                           h_sol[0,1]*gamma.y[0,1] - h_sol[0,-1]*gamma.y[]);
         c[] = num / den;
     }

 #if JACOBI
     foreach_level_or_leaf(l)
         h_sol[] = (h_sol[] + 2. * c[]) / 3.;
 #endif
 }

 /**
  * @brief Main thin-film solver.
  */

/* Deterministic implicit midpoint (including midpoint h^n mobility).
   The stochastic face flux was prepared ONCE from h_old and is unchanged.
   In particular white multiplicative noise remains left-endpoint/Ito.
   This is not a full stochastic-midpoint change of interpretation. */
trace
mgstats solve_thin_film (scalar f, double dt_step, double m_val,
                        double n_val, double dd_val, int timestep) {
    struct EquationParams p;
    p.dt=dt_step;p.f=f;p.m=m_val;p.n=n_val;p.dd=dd_val;
    p.tolerance=tcw_mg_tolerance;
    mgstats s={0};
    if (dt_step==0.) return s;
    scalar old[],previous[],b[],pold[],pnew[];
    const scalar lambda[]=-1./dt_step;
    face vector D[],beta[],gamma[],mobility[],oldflux[];
    foreach() old[]=f[];
    boundary({old});
    foreach() {
        pold[]=0.;
        foreach_dimension() pold[]-=(old[1]+old[-1]-2.*old[])/sq(Delta);
    }
    boundary({pold});
    double update=HUGE_VAL;int iter=0;
    do {
        foreach() previous[]=f[];
        foreach_face() {
            double hmid=fmax(.25*(old[]+old[-1]+f[]+f[-1]),1e-9);
            mobility.x[]=pow(hmid,n_val);
            D.x[]=.5*mobility.x[];
            beta.x[]=gamma.x[]=0.;
            oldflux.x[]=mobility.x[]*face_gradient_x(pold,0);
        }
        boundary((scalar *){D,beta,gamma,mobility,oldflux});
        foreach() {
            b[]=-old[]/dt_step;
            foreach_dimension() b[]-=.5*(oldflux.x[1]-oldflux.x[])/Delta;
            if(tcw_stochastic_rhs_active) b[]-=tcw_stochastic_rhs[];
        }
        boundary({b});restriction({D,beta,gamma,lambda});
        p.lambda=lambda;p.D=D;p.beta=beta;p.gamma=gamma;
        s=mg_solve({f},{b},residual_thin_film,relax_thin_film,&p,
                   4,NULL,TCW_MG_MINLEVEL,tcw_mg_tolerance);
        update=0.;int invalid=0;
        foreach(reduction(max:update) reduction(max:invalid)) {
            double v=fabs(f[]-previous[]);
            if(!isfinite(v))invalid=1;else update=fmax(update,v);
        }
        if(invalid)update=HUGE_VAL;
        boundary({f});iter++;
    } while(update>tcw_picard_tolerance && iter<tcw_picard_max);
    tcw_solver_calls++;
    tcw_solver_max_picard=max(tcw_solver_max_picard,iter);
    tcw_solver_max_mg_cycles=max(tcw_solver_max_mg_cycles,s.i);
    tcw_solver_max_update=fmax(tcw_solver_max_update,update);
    tcw_solver_max_residual=fmax(tcw_solver_max_residual,s.resa);
    if(!isfinite(update)||!isfinite(s.resa)||update>tcw_picard_tolerance||s.resa>tcw_mg_tolerance) {
        tcw_solver_failures++;
        if(pid()==0)fprintf(stderr,"TCW_SOLVER_FAIL midpoint step=%d picard=%d update=%.17g mg_residual=%.17g\n",timestep,iter,update,s.resa);
        return s;
    }
    foreach() {
        pnew[]=0.;
        foreach_dimension()pnew[]-=(f[1]+f[-1]-2.*f[])/sq(Delta);
    }
    boundary({pnew});
    face vector flux[];
    foreach_face()flux.x[]=.5*mobility.x[]*(face_gradient_x(pold,0)+face_gradient_x(pnew,0));
    int joint=TCW_COUPLED_STOCHASTIC&&tcw_stochastic_rhs_active&&tcw_noise_mode>=0;
    if(joint) {
        foreach_face()flux.x[]+=tcw_unlimited_stochastic_flux.x[];
        campaign_joint_raw_minimum=fmin(campaign_joint_raw_minimum,statsf(f).min);
    }
    boundary((scalar *){flux});
    if(TCW_POSITIVE_STOCHASTIC&&tcw_noise_mode>=0)campaign_limit_flux(old,flux,dt_step,joint?2:1);
    foreach() {
        f[]=old[];
        foreach_dimension()f[]+=dt_step*(flux.x[1]-flux.x[])/Delta;
    }
    boundary({f});return s;
}

 /**
  * @brief Adapt the time step.
  */
 double adjust_timestep(double h1, double h2, double dt_curr, int lm_curr, double dtmax_curr, int step_idx) {
     // Check for numerical instability.
     int error = (step_idx > 100 && h2 > 100);

 #if _MPI
     MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
 #endif

     if (error) {
         if (pid() == 0)
             fprintf(stderr, "Error: Instability detected. Step=%d, h_max=%g\n", step_idx, h2);
 #if _MPI
         MPI_Finalize();
 #endif
         exit(EXIT_FAILURE);
     }

     // Find the global minimum height.
     double h_min_local = h1;
     double h_min_global = h1;

 #if _MPI
     MPI_Allreduce(&h_min_local, &h_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
 #endif

     // Adjust dtmax using the stability estimate.
     if (h_min_global < 0.5) {
         dtmax_curr = 0.000134 * pow(55.901699, h_min_global);
     }

     if (pid() == 0) {
         printf("Status: h_min=%g dt=%g dtmax=%g lm=%d t=%g step=%08d\n",
                h_min_global, dt_curr, dtmax_curr, lm_curr, t, step_idx);
     }

     return dtmax_curr;
 }

#include "white_restart.h"
 /* Main program. */
static int jcp_selftest (int kind);
static int jcp_gaussian_selftest (void);
static int coupled_update_selftest (void);
static int tcw_brownian_selftest (void);
int main(int argc, char *argv[]) {
    if (argc == 2 && !strcmp(argv[1], "--brownian-selftest")) return tcw_brownian_selftest();
    if (argc == 2 && !strcmp(argv[1], "--coupled-selftest")) return coupled_update_selftest();
    if (argc >= 2 && !strcmp(argv[1], "--ouse-bin-selftest")) return tcw_ouse_bin_selftest(argc > 2 ? strtod(argv[2],NULL) : .5);
    if (argc == 2 && !strcmp(argv[1], "--gaussian-selftest")) return jcp_gaussian_selftest();
    if (argc == 2 && !strcmp(argv[1], "--operator-selftest")) return jcp_selftest(0);
    if (argc == 2 && !strcmp(argv[1], "--nonfinite-selftest")) return jcp_selftest(1);
    if (argc == 2 && !strcmp(argv[1], "--face-index-selftest")) return jcp_selftest(2);
    if (argc == 2 && !strcmp(argv[1], "--forcing-selftest")) return jcp_selftest(3);
    if (argc == 2 && !strcmp(argv[1], "--implicit-selftest")) return jcp_selftest(4);
    if (argc == 2 && !strcmp(argv[1], "--positivity-selftest")) return jcp_selftest(5);
    unsigned int seed = argc > 1 ? (unsigned int)strtoul(argv[1], NULL, 10) : 1u;
    tcw_noise_seed = seed;
    if (argc > 2) tcw_dt = strtod(argv[2], NULL);
    if (argc > 3) tcw_test_steps = (int)strtol(argv[3], NULL, 10);
    if (argc > 4) tcw_tau = strtod(argv[4], NULL);
    if (argc > 5) tcw_noise_mode = (int)strtol(argv[5], NULL, 10);
    if (argc > 6) tcw_sample_every = (int)strtol(argv[6], NULL, 10);
    if (argc > 7) tcw_noise_snapshot_every = (int)strtol(argv[7], NULL, 10);
    if (argc > 8) tcw_noise_probe_every = (int)strtol(argv[8], NULL, 10);
    double tcw_sigma_cli = argc > 9 ? strtod(argv[9], NULL) : 0.;
    if (argc > 10) tcw_ouse_level = (int)strtol(argv[10], NULL, 10);
    if (argc > 11) tcw_ouse_legacy_scale = (int)strtol(argv[11], NULL, 10);
    if (argc > 12) tcw_ouse_amp_multiplier = strtod(argv[12], NULL);
    if (argc > 13) LEVEL_MIN = LEVEL_MAX = (int)strtol(argv[13], NULL, 10);
    if (argc > 14) tcw_mg_tolerance = strtod(argv[14], NULL);
    if (argc > 15) tcw_picard_tolerance = strtod(argv[15], NULL);
    if (argc > 16) tcw_picard_max = (int)strtol(argv[16], NULL, 10);
    if (argc > 17) tcw_height_floor = strtod(argv[17], NULL);
    if (argc > 18) tcw_varphi_scale = strtod(argv[18], NULL);
    if (argc > 19) tcw_mg_nitermax = (int)strtol(argv[19], NULL, 10);
    if (!(tcw_dt > 0.) || tcw_test_steps < 2) return 2;
    if (tcw_noise_mode == 2 && !(tcw_tau > 0.)) return 3;
    if (tcw_noise_mode == 2 && !(tcw_ouse_amp_multiplier > 0.)) return 5;
    if (!(tcw_varphi_scale >= 0.)) return 8;
    if (LEVEL_MIN < 6 || LEVEL_MIN > 10 || !(tcw_mg_tolerance > 0.) ||
        !(tcw_picard_tolerance > 0.) || tcw_picard_max < 1 ||
        tcw_height_floor < 0.) return 6;
    if (argc > 20) jcp_case = (int)strtol(argv[20], NULL, 10);
    if (argc > 21) campaign_domain = strtod(argv[21], NULL);
    if (argc > 22) tcw_varphi_scale = strtod(argv[22], NULL)/COEFF_VARPHI;
    if (argc > 23) tcw_bridge_depth=(int)strtol(argv[23],NULL,10);
    if (argc > 24) dt_coupling_ratio=(int)strtol(argv[24],NULL,10);
    if (argc > 25) white_lower=strtod(argv[25],NULL);
    if (argc > 26) restart_stride=(int)strtol(argv[26],NULL,10);
    if (argc > 27 && strcmp(argv[27],"-")) restart_path=argv[27];
    restart_total_steps=tcw_test_steps;
    if (argc>25 && (tcw_noise_mode!=0 || white_lower<=0. || white_lower>jcp_hp || restart_stride<0)) return 40;
    if (dt_coupling_ratio!=1 && dt_coupling_ratio!=2 && dt_coupling_ratio!=4) return 32;
    if (tcw_bridge_depth<0 || tcw_bridge_depth>4 || (tcw_bridge_depth && tcw_noise_mode!=0)) return 18;
    if (!(campaign_domain > 2.*jcp_radius) || !isfinite(campaign_domain) ||
        !(tcw_varphi_scale >= 0.) || !isfinite(tcw_varphi_scale)) return 12;
    if (jcp_case < 0 || jcp_case > 5 ||
        (tcw_noise_mode != -1 && tcw_noise_mode != 0 && tcw_noise_mode != 2)) return 9;
    tcw_stop_index = tcw_test_steps - 1;
    periodic(right);
    periodic(top);
    origin(0., 0.);
    L0 = jcp_case == 1 ? campaign_domain : 12.8;
    HOUT = jcp_case == 1 ? jcp_hp : 1.;
    init_grid(1 << LEVEL_MIN);
    if (restart_path) {
        WhiteCheckpoint saved; restart_read_header(&saved);
        restart_offset=saved.completed;
        tcw_test_steps=restart_total_steps-restart_offset;
        tcw_stop_index=tcw_test_steps-1;
    }
    { extern int NITERMAX; NITERMAX = tcw_mg_nitermax; }
    tcw_sigma = tcw_sigma_cli > 0. ? tcw_sigma_cli : 2.*L0/(1 << LEVEL_MIN);
    if (tcw_noise_mode == 2 && !(tcw_sigma > 0.)) return 4;
    if (pid() == 0)
        fprintf(stderr, "TCW_CONFIG level=%d n=%d L0=%.13g dt=%.13g tau=%.13g sigma=%.13g mode=%d steps=%d ouse_grid=%d scale=%s amp_multiplier=%.13g mg_tol=%.13g picard_tol=%.13g picard_max=%d height_floor=%.13g varphi_scale=%.13g mg_nitermax=%d\n",
                LEVEL_MIN, 1 << LEVEL_MIN, L0, tcw_dt, tcw_tau, tcw_sigma,
                tcw_noise_mode, tcw_test_steps, 1 << tcw_ouse_grid_level(),
                tcw_ouse_legacy_scale ? "legacy20260902" : "eq52",
                tcw_ouse_amp_multiplier, tcw_mg_tolerance,
                tcw_picard_tolerance, tcw_picard_max, tcw_height_floor,
                tcw_varphi_scale, tcw_mg_nitermax);
    if (pid() == 0 && tcw_noise_mode == 2 && !tcw_ouse_legacy_scale) {
        /* one-shot Eq.(52) amplitude-ratio diagnostic at h0 = 1. */
        double dxc = tcw_grid_delta();
        double lam = tcw_Lambda(1./tcw_sigma);
        double ups = tcw_Upsilon(dxc/tcw_sigma);
        double th  = tcw_Theta(tcw_dt/tcw_tau);
        double r2  = sq(sq(tcw_sigma))*tcw_tau*lam*ups*th/
                     (2.*pi*dxc*tcw_dt);
        r2 *= tcw_sigma*tcw_Upsilon(dxc/tcw_sigma)/(sqrt(pi)*dxc);
        fprintf(stderr, "TCW_OUSE_RATIO eq52_amp_factor=%.8g base_factor=%.8g multiplier=%.8g lambda=%.6g upsilon=%.6g theta=%.6g dx=%.6g\n",
                (r2 > 0. ? sqrt(r2) : 0.)*tcw_ouse_amp_multiplier,
                r2 > 0. ? sqrt(r2) : 0., tcw_ouse_amp_multiplier,
                lam, ups, th, dxc);
    }
    campaign_log_stride = (int)ceil(.1/tcw_dt - 1e-10);
    if (campaign_log_stride < 1) campaign_log_stride = 1;
    clock_gettime(CLOCK_MONOTONIC, &campaign_started);
    if (pid() == 0)
        fprintf(stderr, "CAMPAIGN_1021 ranks=%d omp_max=%d varphi=%.17g L=%.17g dx=%.17g floor=%.17g sigma=%.17g tau=%.17g\n",
                npe(),
#if _OPENMP
                omp_get_max_threads(),
#else
                1,
#endif
                COEFF_VARPHI*tcw_varphi_scale, L0, tcw_grid_delta(), tcw_height_floor, tcw_sigma, tcw_tau);
    if (pid()==0) fprintf(stderr,"CAPILLARY_MIDPOINT deterministic_midpoint=1 noise_mobility=left noise_generated_once=1\n");
    if (pid() == 0) fprintf(stderr,"CAMPAIGN_SCHEME frozen_mobility=%d mg_minlevel=%d conservative_stochastic_limit=%d positive_floor=%.17g\n",TCW_FROZEN_MOBILITY,TCW_MG_MINLEVEL,TCW_POSITIVE_STOCHASTIC,tcw_noise_mode == 0 ? white_lower : 1e-6);
    if (pid() == 0) fprintf(stderr,"CAMPAIGN_COUPLED_RAW_NOISE enabled=%d limiter_kind=total_flux\n",TCW_COUPLED_STOCHASTIC);
    if (pid()==0 && tcw_noise_mode==0) fprintf(stderr,"WHITE_BROWNIAN_COUPLING enabled=%d depth=%d ratio=%d coarse_dt=%.17g fine_dt=%.17g\n",tcw_bridge_depth>0,tcw_bridge_depth,1<<tcw_bridge_depth,(1<<tcw_bridge_depth)*tcw_dt,tcw_dt);
    if (pid() == 0 && tcw_noise_mode == 2)
        fprintf(stderr,"OUSE_DISCRETE_COV exact_time_bins=%d lateral_cell_bins=%d\n",
                TCW_OUSE_BIN_AVERAGES && !tcw_ouse_legacy_scale,
                TCW_OUSE_BIN_AVERAGES && !tcw_ouse_legacy_scale);
    run();
    return tcw_run_failed ? 7 : 0;
}

event init(i = 0) {
    create_directory(dirFlowFields);
    create_directory(dirFigures);
    create_directory("noise");
    create_directory("contour");
    create_directory("observations");
    create_directory("checkpoints");
    foreach() {
        double r = hypot(x - L0/2., y - L0/2.);
        h[] = jcp_case == 1 ? jcp_hp + (r < jcp_radius ?
            jcp_peak*pow(cos(pi*r/(2.*jcp_radius)), 4.) : 0.) : 1.;
        if (jcp_case == 2 || jcp_case == 3)
            h[] = 1. + 1e-4*cos(2.*pi*(jcp_case == 2 ? x : y)/L0);
        if (jcp_case == 4 || jcp_case == 5) {
            double u = jcp_case == 4 ? x : y, v = jcp_case == 4 ? y : x;
            h[] = 1. + .2*cos(2.*pi*u/L0) + .1*sin(4.*pi*v/L0);
        }
    }
    boundary({h});
    stats hs = statsf(h);
    hmin = hmin1 = hs.min;
    jcp_initial_volume = hs.sum;
    campaign_initial_excess = jcp_initial_volume - HOUT*sq(L0);
    if (restart_path) restart_load();
    write_hfield_binary(restart_path ? "restart_initial.hfield" : "initial.hfield");
    if (pid()==0) {
        FILE *idx=fopen("observations/index.csv","w");
        if (!idx) exit(49);
        fprintf(idx,"step,time,file,n\n"); fclose(idx);
    }
    if (pid() == 0)
        fprintf(stderr, "JCP_INITIAL case=%d volume=%.17g hmin=%.17g L0=%.17g A=0 noise_mode=%d white_precursor_cutoff=%s floor=%.17g amr=OFF\n",
                jcp_case, jcp_initial_volume, hmin, L0, tcw_noise_mode,
                tcw_noise_mode == 0 ? "ON" : "OFF", tcw_height_floor);
}

static uint64_t tcw_mix64 (uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30))*0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27))*0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static double tcw_u01 (uint64_t key) {
    return ((double)(tcw_mix64(key) >> 11) + .5)*(1./9007199254740992.);
}

static int tcw_wrap_n (int i, int n) {
    return (i % n + n) % n;
}

/* --- 2026-09-03: variance functions, Lockerby et al. Eqs. (20)/(49)/(52) --- */
static double tcw_Lambda (double xi) {
    double x2 = xi*xi;
    if (fabs(xi) < .2)
        return x2*x2*(1.5 + x2*(-1./6. + x2*(1./40. + x2*(-1./280. +
               x2*(1./2160. + x2*(-1./18480. + x2/174720.))))));
    return 1. - 3.*x2 - exp(-x2)*(1. - 2.*x2)
         + 2.*sqrt(pi)*x2*xi*erf(xi);
}
static double tcw_Upsilon (double xi) {
    return expm1(-xi*xi) + sqrt(pi)*xi*erf(xi);
}
static double tcw_Theta (double xi) {
    if (fabs(xi) < 1e-3)
        return xi*xi*(.5 + xi*(-1./6. + xi*(1./24. + xi*(-1./120. + xi/720.))));
    return xi + expm1(-xi);
}
/* OU-SE noise-grid level: 0 = follow LEVEL_MIN (20260902 behavior). */
static int tcw_ouse_grid_level (void) {
    if (tcw_ouse_level > 0 && tcw_ouse_level < 17) return tcw_ouse_level;
    return LEVEL_MIN;
}

static double tcw_hash_gaussian (int ix, int iy, unsigned component,
                                 int timestep, int n, uint64_t tag) {
    uint64_t key = tcw_noise_seed ^ tag;
    key ^= (uint64_t)tcw_wrap_n(ix, n)*0x9e3779b97f4a7c15ULL;
    key ^= (uint64_t)tcw_wrap_n(iy, n)*0xbf58476d1ce4e5b9ULL;
    key ^= (uint64_t)component*0x94d049bb133111ebULL;
    key ^= (uint64_t)timestep*0xd6e8feb86659fd93ULL;
    double u1 = tcw_u01(key), u2 = tcw_u01(key ^ 0xa0761d6478bd642fULL);
    return sqrt(-2.*log(u1))*cos(2.*pi*u2);
}

static double tcw_grid_delta (void);

static double tcw_face_gaussian_native (double xx, double yy,
                                        unsigned component, int timestep) {
    int n = 1 << LEVEL_MIN;
    double d = tcw_grid_delta();
    int ix = (int)llround((xx - X0)/d - (component ? .5 : 0.));
    int iy = (int)llround((yy - Y0)/d - (component ? 0. : .5));
    return tcw_hash_gaussian(ix, iy, component, timestep, n,
                             0x243f6a8885a308d3ULL);
}

static double tcw_grid_delta (void) {
    return L0/(1 << LEVEL_MIN);
}

static double tcw_noise_delta (void) {
    return L0/(1 << TCW_NOISE_LEVEL);
}

/* Haar refinement: sum(fine Gaussian)/sqrt(ratio) is exactly the original
   coarse Gaussian. Distinct children retain independent N(0,1) marginals.
   Spatial sampling, amplitude, mask and physical white-noise law are unchanged. */
static double tcw_face_gaussian_halfstep (double xx, double yy,
                                          unsigned component, int fine_step) {
    int n=1 << LEVEL_MIN, ratio=1 << tcw_bridge_depth, coarse_step=fine_step/ratio;
    int sub=fine_step%ratio;
    double d=tcw_grid_delta();
    int ix=(int)llround((xx-X0)/d-(component ? .5 : 0.));
    int iy=(int)llround((yy-Y0)/d-(component ? 0. : .5));
    double value=tcw_face_gaussian_native(xx,yy,component,coarse_step);
    for (int level=0;level<tcw_bridge_depth;level++) {
        int node=sub >> (tcw_bridge_depth-level);
        int bit=(sub >> (tcw_bridge_depth-level-1)) & 1;
        uint64_t salt=0x9b05688c2b3e6c1fULL ^ (0x94d049bb133111ebULL*(uint64_t)(level+1));
        double bridge=tcw_hash_gaussian(ix,iy,component,coarse_step*(1 << level)+node,n,salt);
        value=(value+(bit ? -bridge : bridge))/sqrt(2.);
    }
    return value;
}

static int tcw_brownian_selftest (void) {
    LEVEL_MIN=LEVEL_MAX=6;L0=16.;tcw_noise_seed=102125;
    periodic(right);periodic(top);init_grid(1 << LEVEL_MIN);
    int n=1 << LEVEL_MIN;double d=tcw_grid_delta();
    for (int depth=1;depth<=4;depth*=2) {
        tcw_bridge_depth=depth;
        int ratio=1 << depth;
        double means[16]={0.}, squares[16]={0.}, cross[16][16]={{0.}};
        double count=0.,max_error=0.;
        for (int it=0;it<16;it++) for (unsigned c=0;c<2;c++)
            for (int iy=0;iy<n;iy++) for (int ix=0;ix<n;ix++) {
                double xx=(ix+(c ? .5 : 0.))*d,yy=(iy+(c ? 0. : .5))*d;
                double z[16],sum=0.;
                for (int sub=0;sub<ratio;sub++) {
                    z[sub]=tcw_face_gaussian_halfstep(xx,yy,c,ratio*it+sub);
                    sum+=z[sub];means[sub]+=z[sub];squares[sub]+=sq(z[sub]);
                }
                for (int j=0;j<ratio;j++) for (int k=j+1;k<ratio;k++) cross[j][k]+=z[j]*z[k];
                double g=tcw_face_gaussian_native(xx,yy,c,it);
                max_error=fmax(max_error,fabs(sum/sqrt((double)ratio)-g));
                max_error=fmax(max_error,fabs(z[0]-tcw_face_gaussian_halfstep(xx+L0,yy,c,ratio*it)));
                max_error=fmax(max_error,fabs(z[0]-tcw_face_gaussian_halfstep(xx,yy+L0,c,ratio*it)));
                count+=1.;
            }
        double max_mean=0.,max_var_error=0.,max_cov=0.;
        for (int j=0;j<ratio;j++) {
            means[j]/=count;
            double variance=squares[j]/count-sq(means[j]);
            max_mean=fmax(max_mean,fabs(means[j]));
            max_var_error=fmax(max_var_error,fabs(variance-1.));
        }
        for (int j=0;j<ratio;j++) for (int k=j+1;k<ratio;k++)
            max_cov=fmax(max_cov,fabs(cross[j][k]/count-means[j]*means[k]));
        int pass=max_error<1e-12 && max_mean<.02 && max_var_error<.02 && max_cov<.02;
        if (pid()==0) fprintf(stdout,"BROWNIAN_TREE_SELFTEST depth=%d ratio=%d aggregate_periodic_error=%.17g max_mean=%.12g max_var_error=%.12g max_offdiag_cov=%.12g result=%s\n",
            depth,ratio,max_error,max_mean,max_var_error,max_cov,pass ? "PASS" : "FAIL");
        if (!pass) return 17;
    }
    return 0;
}

static double tcw_face_gaussian_l6 (double xx, double yy, unsigned component,
                                    int timestep) {
    double d = tcw_noise_delta();
    double gx = (xx - X0)/d, gy = (yy - Y0)/d;
    int ix = (int)floor(gx), iy = (int)floor(gy);
    double tx = gx - ix, ty = gy - iy;
    double z = (1.-tx)*(1.-ty)*tcw_hash_gaussian(ix, iy, component, timestep,
                                                  1 << TCW_NOISE_LEVEL,
                                                  0x13198a2e03707344ULL)
             + tx*(1.-ty)*tcw_hash_gaussian(ix+1, iy, component, timestep,
                                             1 << TCW_NOISE_LEVEL,
                                             0x13198a2e03707344ULL)
             + (1.-tx)*ty*tcw_hash_gaussian(ix, iy+1, component, timestep,
                                            1 << TCW_NOISE_LEVEL,
                                            0x13198a2e03707344ULL)
             + tx*ty*tcw_hash_gaussian(ix+1, iy+1, component, timestep,
                                       1 << TCW_NOISE_LEVEL,
                                       0x13198a2e03707344ULL);
    double var = ((1.-tx)*(1.-tx) + tx*tx)*
                 ((1.-ty)*(1.-ty) + ty*ty);
    return z/sqrt(var);
}

static void tcw_fft1d (double complex *a, int n, int inverse) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double complex tmp = a[i]; a[i] = a[j]; a[j] = tmp;
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        double angle = (inverse ? 2. : -2.)*pi/len;
        double complex wlen = cos(angle) + I*sin(angle);
        for (int i = 0; i < n; i += len) {
            double complex w = 1.;
            for (int j = 0; j < len/2; j++) {
                double complex u = a[i+j];
                double complex v = a[i+j+len/2]*w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse)
        for (int i = 0; i < n; i++) a[i] /= n;
}

static void tcw_fft2d (double complex *a, int n, int inverse) {
    for (int y = 0; y < n; y++) tcw_fft1d(a + (size_t)y*n, n, inverse);
    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) tcw_se_col[y] = a[(size_t)y*n+x];
        tcw_fft1d(tcw_se_col, n, inverse);
        for (int y = 0; y < n; y++) a[(size_t)y*n+x] = tcw_se_col[y];
    }
}

static void tcw_se_ensure (void) {
    int want = 1 << tcw_ouse_grid_level();
    if (tcw_noise_n == want && tcw_se_field[0]) return;
    tcw_noise_n = want;
    size_t nn = (size_t)tcw_noise_n*tcw_noise_n;
    for (int c = 0; c < 2; c++)
        tcw_se_field[c] = calloc(nn, sizeof(double));
    for (int c = 0; c < 2; c++)
        tcw_se_endpoint[c] = calloc(nn, sizeof(double));
    tcw_se_first = calloc(nn, sizeof(double));
    tcw_se_innovation = calloc(nn, sizeof(double));
    tcw_se_fft = calloc(nn, sizeof(double complex));
    tcw_se_col = calloc(tcw_noise_n, sizeof(double complex));
    tcw_se_filter = calloc(tcw_noise_n, sizeof(double));
    if (!tcw_se_field[0] || !tcw_se_field[1] || !tcw_se_endpoint[0] ||
        !tcw_se_endpoint[1] || !tcw_se_first || !tcw_se_innovation ||
        !tcw_se_fft || !tcw_se_col || !tcw_se_filter) {
        fprintf(stderr, "OU-SE allocation failed for %d x %d field\n",
                tcw_noise_n, tcw_noise_n);
        exit(EXIT_FAILURE);
    }
    if (TCW_OUSE_BIN_AVERAGES && !tcw_ouse_legacy_scale) {
        /* Exact periodic covariance of SE cell averages. Taking the DFT of
           covariance includes all sampling aliases; truncating exp*sinc at
           Nyquist is not accurate when sigma is smaller than a cell.
           The unchanged Eq.(52) supplies the physical marginal variance. */
        double dc=tcw_grid_delta()/tcw_sigma, cellvar=tcw_Upsilon(dc)/sq(dc);
        double *cov=calloc(tcw_noise_n,sizeof(double));
        if (!cov || !(cellvar>0.)) exit(11);
        for (int j=0;j<tcw_noise_n;j++) {
            double value=0.;
            for (int image=-1;image<=1;image++) {
                double distance=fabs((j+image*tcw_noise_n)*tcw_grid_delta());
                if (distance>tcw_grid_delta()+8.*tcw_sigma) continue;
                double u=distance/tcw_sigma;
                value+=(tcw_Upsilon(u+dc)+tcw_Upsilon(u-dc)-2.*tcw_Upsilon(u))/(2.*sq(dc));
            }
            cov[j]=value/cellvar;
        }
        for (int k=0;k<tcw_noise_n;k++) {
            double eigen=0.;
            for (int j=0;j<tcw_noise_n;j++)
                eigen+=cov[j]*cos(2.*pi*k*j/tcw_noise_n);
            if (eigen < -1e-10) {
                fprintf(stderr,"OUSE_CELL_COV_NOT_PSD k=%d eigen=%.17g\\n",k,eigen);
                exit(11);
            }
            tcw_se_filter[k]=sqrt(fmax(0.,eigen));
        }
        free(cov);
    } else {
        for (int j=0;j<tcw_noise_n;j++) {
            int k=j<=tcw_noise_n/2 ? j : j-tcw_noise_n;
            double q=2.*pi*k/L0;
            tcw_se_filter[j]=exp(-.125*sq(tcw_sigma*q));
        }
    }
    double sumw2 = 0.;
    for (int j = 0; j < tcw_noise_n; j++) sumw2 += sq(tcw_se_filter[j]);
    /* Unit marginal variance after deleting the actual filtered DC mode. */
    tcw_se_norm = sqrt(sumw2*sumw2 - sq(sq(tcw_se_filter[0])))/tcw_noise_n;
    if (pid()==0) fprintf(stderr,"COVARIANCE_NORM_AUDIT filter0=%.17g sum_eigen=%.17g norm=%.17g legacy_norm=%.17g\n",tcw_se_filter[0],sumw2,tcw_se_norm,sqrt(sumw2*sumw2-1.)/tcw_noise_n);
    if (!(tcw_se_norm > 0.) || !isfinite(tcw_se_norm)) {
        fprintf(stderr, "Unresolved zero-mean SE covariance\n"); exit(11);
    }
}

static void tcw_se_innovation_field (unsigned component, int timestep) {
    int n = tcw_noise_n;
    size_t nn = (size_t)n*n;
    for (int iy = 0; iy < n; iy++)
        for (int ix = 0; ix < n; ix++) {
            size_t p = (size_t)iy*n + ix;
            tcw_se_fft[p] = tcw_hash_gaussian(ix, iy, component, timestep, n,
                                               0x3c6ef372fe94f82bULL);
        }
    tcw_fft2d(tcw_se_fft, n, 0);
    double dx = L0/n;
    for (int iy = 0; iy < n; iy++) {
        for (int ix = 0; ix < n; ix++) {
            double amp = tcw_se_filter[ix]*tcw_se_filter[iy];
            size_t p = (size_t)iy*n + ix;
            tcw_se_fft[p] = p ? tcw_se_fft[p]*amp : 0.;
        }
    }
    tcw_fft2d(tcw_se_fft, n, 1);
    /* Fixed discrete covariance normalization; no random RMS conditioning.
       Zero Fourier mode is omitted as in the original mean-removal convention. */
    for (size_t p = 0; p < nn; p++)
        tcw_se_innovation[p] = creal(tcw_se_fft[p])/tcw_se_norm;
    (void)dx;
}

/* Real SE filtering is linear and identical for the two independent
   fields. Pack their original RNG keys into real/imaginary parts. */
static void tcw_se_innovation_pair (unsigned c1, unsigned c2, int timestep) {
    int n=tcw_noise_n;
    for (int iy=0;iy<n;iy++) for (int ix=0;ix<n;ix++) {
        size_t p=(size_t)iy*n+ix;
        double a=tcw_hash_gaussian(ix,iy,c1,timestep,n,0x3c6ef372fe94f82bULL);
        double b=tcw_hash_gaussian(ix,iy,c2,timestep,n,0x3c6ef372fe94f82bULL);
        tcw_se_fft[p]=a+I*b;
    }
    tcw_fft2d(tcw_se_fft,n,0);
    for (int iy=0;iy<n;iy++) for (int ix=0;ix<n;ix++) {
        size_t p=(size_t)iy*n+ix;
        tcw_se_fft[p]=p ? tcw_se_fft[p]*tcw_se_filter[ix]*tcw_se_filter[iy] : 0.;
    }
    tcw_fft2d(tcw_se_fft,n,1);
    for (size_t p=0;p<(size_t)n*n;p++) {
        tcw_se_first[p]=creal(tcw_se_fft[p])/tcw_se_norm;
        tcw_se_innovation[p]=cimag(tcw_se_fft[p])/tcw_se_norm;
    }
}

static void tcw_prepare_noise_local (int timestep, double dt_step) {
    if (tcw_noise_mode != 2) return;
    tcw_se_ensure();
    double alpha = exp(-dt_step/tcw_tau);
    double beta = sqrt(fmax(0., 1.-alpha*alpha));
    size_t nn = (size_t)tcw_noise_n*tcw_noise_n;
    if (TCW_OUSE_BIN_AVERAGES && !tcw_ouse_legacy_scale) {
        /* Joint exact law of the endpoint and integral of a stationary OU
           process. The normalized interval field has Eq.(52)'s variance AND
           the covariance of adjacent time-bin averages, not exp(-dt/tau). */
        double x = dt_step/tcw_tau, one_minus = -expm1(-x);
        double endpoint_sd = sqrt(-expm1(-2.*x));
        double var_eta = fabs(x)<.01 ?
            sq(tcw_tau)*x*x*x*(2./3.+x*(-.5+x*(7./30.+x*(-1./12.+x*31./1260.)))) :
            sq(tcw_tau)*(2.*x-4.*one_minus-expm1(-2.*x));
        double shared = tcw_tau*sq(one_minus)/endpoint_sd;
        double independent2 = var_eta-sq(shared);
        double integral_sd = tcw_tau*sqrt(2.*tcw_Theta(x));
        if (!(independent2>0.) || !(integral_sd>0.) || !isfinite(independent2)) {
            fprintf(stderr,"OUSE_BIN_BAD_COEFFICIENT dt=%.17g tau=%.17g\n",dt_step,tcw_tau);
            exit(14);
        }
        double independent = sqrt(independent2);
        for (int c = 0; c < 2; c++) {
            if (!tcw_noise_initialized) {
                tcw_se_innovation_field(c+4,0);
                memcpy(tcw_se_endpoint[c],tcw_se_innovation,nn*sizeof(double));
            }
            tcw_se_innovation_pair(c,c+2,timestep);
            for (size_t p = 0; p < nn; p++) {
                double old = tcw_se_endpoint[c][p], z1 = tcw_se_first[p];
                tcw_se_field[c][p] = (tcw_tau*one_minus*old+shared*z1+
                                      independent*tcw_se_innovation[p])/integral_sd;
                tcw_se_endpoint[c][p] = alpha*old+endpoint_sd*z1;
            }
        }
        tcw_noise_initialized = 1;
        return;
    }
    for (int c = 0; c < 2; c++) {
        tcw_se_innovation_field(c, timestep);
        if (!tcw_noise_initialized)
            for (size_t p = 0; p < nn; p++) tcw_se_field[c][p] = tcw_se_innovation[p];
        else
            for (size_t p = 0; p < nn; p++)
                tcw_se_field[c][p] = alpha*tcw_se_field[c][p] + beta*tcw_se_innovation[p];
    }
    tcw_noise_initialized = 1;
}

/* Generate the same global counter-based fields once. Only the two
   interval flux fields are distributed; endpoint state stays on rank zero.
   This does not change the OU endpoint/integral law or random-number keys. */
static void tcw_prepare_noise (int timestep, double dt_step) {
    if (tcw_noise_mode != 2) return;
#if _MPI
    tcw_se_ensure();
    if (pid()==0) {
        if (dt_coupling_ratio==1) tcw_prepare_noise_local(timestep,dt_step);
        else {
            size_t nn=(size_t)tcw_noise_n*tcw_noise_n;
            for (int c=0;c<2;c++) {
                if (!dt_sum[c]) dt_sum[c]=calloc(nn,sizeof(double));
                if (!dt_sum[c]) exit(33);
                for (size_t p=0;p<nn;p++) dt_sum[c][p]=0.;
            }
            double fine=dt_step/dt_coupling_ratio;
            for (int j=0;j<dt_coupling_ratio;j++) {
                tcw_prepare_noise_local(timestep*dt_coupling_ratio+j,fine);
                for (int c=0;c<2;c++) for (size_t p=0;p<nn;p++)
                    dt_sum[c][p]+=tcw_se_field[c][p];
            }
            /* Sum exact OU integrals, then normalize as a single coarse bin. */
            double scale=sqrt(tcw_Theta(fine/tcw_tau)/tcw_Theta(dt_step/tcw_tau));
            for (int c=0;c<2;c++) for (size_t p=0;p<nn;p++)
                tcw_se_field[c][p]=scale*dt_sum[c][p];
        }
    }
    for (int c=0;c<2;c++)
        if (MPI_Bcast(tcw_se_field[c],tcw_noise_n*tcw_noise_n,MPI_DOUBLE,
                      0,MPI_COMM_WORLD)!=MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD,31);
    tcw_noise_initialized=1;
#else
    tcw_prepare_noise_local(timestep,dt_step);
#endif
}

static int tcw_ouse_bin_selftest (double sigma) {
    LEVEL_MIN=LEVEL_MAX=8;tcw_ouse_level=8;L0=64.;
    tcw_noise_mode=2;tcw_tau=.01;tcw_sigma=sigma;tcw_noise_seed=1021;
    periodic(right);periodic(top);init_grid(1 << LEVEL_MIN);
    int n=1 << LEVEL_MIN, count=256, block=32;
    size_t nn=(size_t)n*n;
    double *previous=calloc(nn,sizeof(double)), *sum=calloc(nn,sizeof(double));
    if (!previous || !sum) return 15;
    double variance=0.,lag=0.,spatial=0.,cross_components=0.,block_variance=0.,block_count=0.;
    for (int it=0;it<count;it++) {
        tcw_prepare_noise(it,tcw_tau);
        for (size_t p=0;p<nn;p++) {
            double z=tcw_se_field[0][p];variance+=z*z;
            size_t neighbor=(p/(size_t)n)*(size_t)n+(p+1)%(size_t)n;
            spatial+=z*tcw_se_field[0][neighbor];
            cross_components+=z*tcw_se_field[1][p];
            if (it) lag+=z*previous[p];
            previous[p]=z;sum[p]+=z;
        }
        if ((it+1)%block==0) {
            for (size_t p=0;p<nn;p++) {block_variance+=sq(sum[p]);sum[p]=0.;}
            block_count+=1.;
        }
    }
    variance/=count*(double)nn;
    double lag1=lag/((count-1.)*nn)/variance;
    double expected_lag=sq(-expm1(-1.))/(2.*tcw_Theta(1.));
    double dcell=tcw_grid_delta()/tcw_sigma, ups=tcw_Upsilon(dcell);
    double bin_variance=sq(ups/sq(dcell)), zero_mode=pi*sq(tcw_sigma/L0);
    double rho_space=(tcw_Upsilon(2.*dcell)-2.*ups)/(2.*ups);
    double expected_space=(rho_space*bin_variance-zero_mode)/(bin_variance-zero_mode);
    double spatial1=spatial/(count*(double)nn)/variance;
    double cross_xy=cross_components/(count*(double)nn);
    double expected_block=tcw_Theta((double)block)/tcw_Theta(1.);
    double block_ratio=block_variance/(block_count*nn*expected_block);
    int pass=fabs(variance-1.)<.025 && fabs(lag1-expected_lag)<.025 && fabs(block_ratio-1.)<.04 &&
             fabs(spatial1-expected_space)<.01 && fabs(cross_xy)<.02;
    if (pid()==0) fprintf(stdout,"OUSE_BIN_SELFTEST sigma=%.12g variance=%.12g lag1=%.12g expected_lag1=%.12g block32_ratio=%.12g spatial1=%.12g expected_spatial1=%.12g cross_xy=%.12g result=%s\n",
                          tcw_sigma,variance,lag1,expected_lag,block_ratio,spatial1,expected_space,cross_xy,pass ? "PASS" : "FAIL");
    free(previous);free(sum);
    return pass ? 0 : 16;
}

static double tcw_se_sample (double xx, double yy, unsigned component) {
    int n = tcw_noise_n;
    double d = L0/n;
    /* component 0: x faces (i,j+1/2); component 1: y faces (i+1/2,j).
       A stationary periodic FFT field is translation invariant, so assigning the
       two independent arrays to these shifted lattices preserves the SE spectrum. */
    double gx = (xx - X0)/d - (component ? .5 : 0.);
    double gy = (yy - Y0)/d - (component ? 0. : .5);
    int ix = (int)floor(gx), iy = (int)floor(gy);
    double tx = gx - ix, ty = gy - iy;
    ix = tcw_wrap_n(ix, n); iy = tcw_wrap_n(iy, n);
    int ix1 = tcw_wrap_n(ix+1, n), iy1 = tcw_wrap_n(iy+1, n);
    double *a = tcw_se_field[component];
    return (1.-tx)*(1.-ty)*a[(size_t)iy*n+ix]
         + tx*(1.-ty)*a[(size_t)iy*n+ix1]
         + (1.-tx)*ty*a[(size_t)iy1*n+ix]
         + tx*ty*a[(size_t)iy1*n+ix1];
}

static double tcw_noise_value (double xx, double yy, unsigned component,
                               int timestep) {
    if (tcw_noise_mode == -1) return 0.;
    if (tcw_noise_mode == 0)
        return tcw_bridge_depth ? tcw_face_gaussian_halfstep(xx, yy, component, timestep) : tcw_face_gaussian_native(xx, yy, component, timestep);
    if (tcw_noise_mode == 2)
        return tcw_se_sample(xx, yy, component);
    return tcw_face_gaussian_l6(xx, yy, component, timestep);
}

/* White-noise-only precursor-film mask. It is zero at h <= h_p and rises
   smoothly over 0.005; OU-SE retains the full-field pilot convention. */
static double tcw_white_precursor_weight (double height) {
    double q = (height - jcp_hp)/0.005;
    q = q < 0. ? 0. : (q > 1. ? 1. : q);
    return q*q*(3. - 2.*q);
}

static double tcw_time_scale (double dt_step) {
    if (tcw_noise_mode != 2 || !tcw_ouse_legacy_scale)
        return 1./sqrt(dt_step); /* Eq.(52) carries all tau corrections via tcw_Theta */
    double x = dt_step/tcw_tau;
    double theta = x + expm1(-x);
    return sqrt(tcw_tau*theta)/dt_step;
}

static void tcw_write_noise_snapshot (int timestep) {
    if (tcw_noise_snapshot_every <= 0 || timestep <= 0 ||
        timestep % tcw_noise_snapshot_every) return;
    if (pid() != 0) return;
    int n = tcw_noise_mode == 2 ? tcw_noise_n : 1 << LEVEL_MIN;
    double dn = L0/n;
    char name[128];
    for (int c = 0; c < 2; c++) {
        snprintf(name, sizeof(name), "noise/noise_%c_%06d.bin", c ? 'y' : 'x', timestep);
        FILE *fp = fopen(name, "wb");
        if (!fp) continue;
        fwrite(&n, sizeof(int), 1, fp);
        for (int iy = 0; iy < n; iy++)
            for (int ix = 0; ix < n; ix++) {
                double z;
                if (tcw_noise_mode == 2) z = tcw_se_field[c][(size_t)iy*n+ix];
                else z = tcw_noise_value(X0+(ix+.5)*dn, Y0+(iy+.5)*dn,
                                         c, timestep);
                fwrite(&z, sizeof(double), 1, fp);
            }
        fclose(fp);
    }
    /* Save the dimensionless noise at the exact x/y face locations used below.
       The raw OU-grid files above remain useful for checking the FFT generator. */
    int nf = 1 << LEVEL_MIN;
    double df = L0/nf;
    for (int c = 0; c < 2; c++) {
        snprintf(name, sizeof(name), "noise/noise_face_%c_%06d.bin", c ? 'y' : 'x', timestep);
        FILE *fp = fopen(name, "wb");
        if (!fp) continue;
        fwrite(&nf, sizeof(int), 1, fp);
        for (int iy = 0; iy < nf; iy++)
            for (int ix = 0; ix < nf; ix++) {
                double xx = X0 + (ix + (c ? .5 : 0.))*df;
                double yy = Y0 + (iy + (c ? 0. : .5))*df;
                double z = tcw_noise_value(xx, yy, c, timestep);
                fwrite(&z, sizeof(double), 1, fp);
            }
        fclose(fp);
    }
    tcw_noise_snapshot_index++;
}

static void tcw_write_noise_probe (int timestep, double dt_step) {
    static int header_written = 0;
    if (tcw_noise_probe_every <= 0 || timestep <= 0 ||
        timestep % tcw_noise_probe_every) return;
    if (pid() != 0) return;
    FILE *fp = fopen("noise_probe.csv", header_written ? "a" : "w");
    if (!fp) return;
    if (!header_written) {
        fprintf(fp, "step,t,dt,mode,xface0,yface0,xface_mid,yface_mid\n");
        header_written = 1;
    }
    double d = tcw_grid_delta();
    double xm = X0 + .5*L0, ym = Y0 + .5*L0;
    fprintf(fp, "%d,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g\n",
            timestep, t, dt_step, tcw_noise_mode,
            tcw_noise_value(X0,Y0+.5*d,0,timestep),
            tcw_noise_value(X0+.5*d,Y0,1,timestep),
            tcw_noise_value(xm,ym+.5*d,0,timestep),
            tcw_noise_value(xm+.5*d,ym,1,timestep));
    fclose(fp);
}

static void campaign_limit_flux (scalar base, face vector flux, double dt_step, int kind) {
    scalar fraction[];
    foreach() {
        double outgoing=0.;
        foreach_dimension() outgoing+=fmax(-flux.x[1],0.)+fmax(flux.x[],0.);
        double lower=tcw_noise_mode == 0 ? white_lower : 1e-6;
        double available=fmax(base[]-lower,0.);
        fraction[]=outgoing>0. ? fmin(1.,available*Delta/(dt_step*outgoing)) : 1.;
    }
    boundary({fraction});
    double raw2=0., limited2=0., active=0., hits=0., removed2=0., minweight=1.;
    foreach_face(reduction(+:raw2) reduction(+:limited2) reduction(+:active) reduction(+:hits) reduction(+:removed2) reduction(min:minweight)) {
        double before=flux.x[];
        double weight=before>0. ? fraction[] : fraction[-1];
        flux.x[]*=weight;
        /* Count each physical face once across MPI partitions. */
        if (is_local(cell)) {
            raw2+=sq(before); limited2+=sq(flux.x[]);
            removed2+=sq(flux.x[]-before);
            if (fabs(before)>1e-20) minweight=fmin(minweight,weight);
            if (fabs(before)>1e-20) {active+=1.; if (weight<1.-1e-12) hits+=1.;}
        }
    }
    boundary((scalar *){flux});
    campaign_flux_raw2[kind]+=raw2; campaign_flux_limited2[kind]+=limited2;
    campaign_flux_active[kind]+=active; campaign_flux_hits[kind]+=hits;
    scan_removed2[kind]+=removed2;
    scan_minweight[kind]=fmin(scan_minweight[kind],minweight);
    if (hits>0. && scan_first[kind]<0) scan_first[kind]=scan_step+1;
}

void prepare_stochastic_rhs (scalar f_old, double dt_step, int timestep) {
    scan_step=timestep;
    double opt_started=opt_clock();
    tcw_prepare_noise(timestep, dt_step);
    opt_noise_seconds+=opt_clock()-opt_started;
    tcw_write_noise_snapshot(timestep);
    tcw_write_noise_probe(timestep, dt_step);
    double dnoise = tcw_noise_mode == 1 ? tcw_noise_delta() : tcw_grid_delta();
    double noise_scale = tcw_time_scale(dt_step);
    /* 2026-09-03: Eq.(52)/Eq.(51) amplitude ratio,
       sqrt(sigma^4 tau Lambda(h/sigma) Upsilon(Dx/sigma) Theta(dt/tau)/(2 pi h^3 Dx dt));
       -> 1 as Dx/sigma, dt/tau -> inf (exact white-noise limit). Dx is the solver
       cell size Delta: the element the scheme averages over (paper Sec. 4.2). */
    int use_eq52 = tcw_noise_mode == 2 && !tcw_ouse_legacy_scale;
    face vector flux[];
    foreach_face(x) {
        double ff = (f_old[] + f_old[-1])/2.;
        if (ff <= 0.) ff = 1e-9;
        double precursor_weight = tcw_noise_mode == 0 ? tcw_white_precursor_weight(ff) : 1.;
        double amp = precursor_weight*sqrt(2.*COEFF_VARPHI*tcw_varphi_scale)*pow(ff, 1.5)*noise_scale/dnoise;
        if (use_eq52)
            amp *= sqrt(fmax(sq(sq(tcw_sigma))*tcw_tau*
                             tcw_Lambda(ff/tcw_sigma)*
                             tcw_Upsilon(Delta/tcw_sigma)*
                             tcw_Theta(dt_step/tcw_tau), 0.)/
                        (2.*pi*ff*ff*ff*Delta*dt_step))*tcw_ouse_amp_multiplier*
                   sqrt(tcw_sigma*tcw_Upsilon(Delta/tcw_sigma)/(sqrt(pi)*Delta));
        /* Second lateral average: Eq.(49) separability, absent in the 1-D Eq.(52). */
        flux.x[] = amp*tcw_noise_value(x, y, 0, timestep);
    }
    foreach_face(y) {
        double ff = (f_old[] + f_old[0,-1])/2.;
        if (ff <= 0.) ff = 1e-9;
        double precursor_weight = tcw_noise_mode == 0 ? tcw_white_precursor_weight(ff) : 1.;
        double amp = precursor_weight*sqrt(2.*COEFF_VARPHI*tcw_varphi_scale)*pow(ff, 1.5)*noise_scale/dnoise;
        if (use_eq52)
            amp *= sqrt(fmax(sq(sq(tcw_sigma))*tcw_tau*
                             tcw_Lambda(ff/tcw_sigma)*
                             tcw_Upsilon(Delta/tcw_sigma)*
                             tcw_Theta(dt_step/tcw_tau), 0.)/
                        (2.*pi*ff*ff*ff*Delta*dt_step))*tcw_ouse_amp_multiplier*
                   sqrt(tcw_sigma*tcw_Upsilon(Delta/tcw_sigma)/(sqrt(pi)*Delta));
        /* Second lateral average: Eq.(49) separability, absent in the 1-D Eq.(52). */
        flux.y[] = amp*tcw_noise_value(x, y, 1, timestep);
    }
    if (timestep==0) {
        double sum2=0.,faces=0.,thin2=0.,thinfaces=0.;
        foreach_face(reduction(+:sum2) reduction(+:faces) reduction(+:thin2) reduction(+:thinfaces)) {
            if (is_local(cell)) {
                sum2+=sq(flux.x[]);faces+=1.;
                double hf=(f_old[]+f_old[-1])/2.;
                if (hf<=.015) {thin2+=sq(flux.x[]);thinfaces+=1.;}
            }
        }
        if (pid()==0) fprintf(stderr,"FLUX_AMPLITUDE_AUDIT mode=%d phi=%.17g raw_rms=%.17g precursor_rms=%.17g faces=%.17g precursor_faces=%.17g\n",tcw_noise_mode,COEFF_VARPHI*tcw_varphi_scale,sqrt(sum2/faces),thinfaces>0.?sqrt(thin2/thinfaces):0.,faces,thinfaces);
    }
    boundary((scalar *){flux});
    if (TCW_COUPLED_STOCHASTIC && tcw_noise_mode>=0) {
        foreach_face() tcw_unlimited_stochastic_flux.x[]=flux.x[];
        boundary((scalar *){tcw_unlimited_stochastic_flux});
    }
    else if (TCW_POSITIVE_STOCHASTIC && tcw_noise_mode >= 0)
        campaign_limit_flux(f_old,flux,dt_step,0);
    foreach() {
        tcw_stochastic_rhs[] = 0.;
        foreach_dimension()
            tcw_stochastic_rhs[] += (flux.x[1] - flux.x[])/Delta;
    }
    boundary({tcw_stochastic_rhs});
}

/* 2026-09-10: per-step Fig.7 observable, same definition as the validated 1-D
   solver (space mean of |h_new-h_old|/dt per step). hnew is the post-limiter
   film, hold the pre-update film; foreach reductions handle MPI themselves. */
static void tcw_write_film_speed (scalar hnew, scalar hold,
                                  double dt_step, int timestep,
                                  int limiter_hits) {
    double spd = 0., area = 0.;
    foreach(reduction(+:spd) reduction(+:area)) {
        spd += fabs(hnew[] - hold[])/dt_step*sq(Delta);
        area += sq(Delta);
    }
    if (pid() == 0) {
        static FILE *fp = NULL;
        if (!fp) {
            fp = fopen("film_speed.csv", "w");
            if (!fp) {
                perror("film_speed.csv");
                exit(EXIT_FAILURE);
            }
            fprintf(fp, "step,t_start,t_end,dt,mean_abs_hdot,limiter_hits\n");
        }
        fprintf(fp, "%d,%.17g,%.17g,%.17g,%.17g,%d\n",
                timestep + 1, restart_now()-dt_step, restart_now(), dt_step,
                area > 0. ? spd/area : NAN, limiter_hits);
        fflush(fp);
    }
}
event stability(i++) {
    if (tcw_fixed_dt) {
        dt = dtnext(tcw_dt);
        return 0.;
    }
    double cap = hmin < 0.1 ? 1e-3 : (hmin < 0.5 ? 5e-3 : 1e-2);
    if (i > 0 && i % 20 == 0 && tcw_dt < cap) tcw_dt *= 1.25;
    if (tcw_dt > cap) tcw_dt = cap;
    if (tcw_dt < 1e-3) tcw_dt = 1e-3;
    dt = dtnext(tcw_dt);
}



event integration(i++) {
 #if T2OR
     // Second-order time integration.
     dt = dtnext(ckstp(DTMIN, DTMAX, dtold, errorstep));
     dtold = dt;
     foreach() {
         h2[] = h[];
         h1[] = h[];
     }
     boundary({h2, h1});

     solve_thin_film(h1, dt, m_exp, n_exp, dd_exp, i);
     solve_thin_film(h2, dt / 2., m_exp, n_exp, dd_exp, i);
     solve_thin_film(h2, dt / 2., m_exp, n_exp, dd_exp, i);

     double errormax = 0;
     foreach() {
         h[] = 2. * h2[] - h1[];
         err[] = fabs(h2[] - h1[]);
         if (err[] > errormax)
             errormax = err[];
     }
     boundary({h2, h1});
     stats s = statsf(err);
     errorstep = sqrt(s.sum);
#else
     // First-order IMEX: old-state multiplicative noise enters the same
     // implicit solve as a fixed right-hand-side flux divergence.
     dtold = dt;
     foreach() h_noise_old[] = h[];
     boundary({h_noise_old});
    prepare_stochastic_rhs(h_noise_old, dt, i+restart_offset);
    if (!TCW_COUPLED_STOCHASTIC && TCW_POSITIVE_STOCHASTIC && tcw_noise_mode >= 0) {
        foreach() h[]=h_noise_old[]+dt*tcw_stochastic_rhs[];
        boundary({h});
    }
#if !TCW_DIAGNOSTIC_NO_DETERMINISTIC
    tcw_stochastic_rhs_active = TCW_COUPLED_STOCHASTIC || !(TCW_POSITIVE_STOCHASTIC && tcw_noise_mode >= 0);
    double opt_solve_started=opt_clock();
    solve_thin_film(h, dt, m_exp, n_exp, dd_exp, i+restart_offset);
    opt_solve_seconds+=opt_clock()-opt_solve_started;
    tcw_stochastic_rhs_active = 0;
#else
    foreach() h[] += dt*tcw_stochastic_rhs[];
    boundary({h});
#endif
#endif

     // Update diagnostics and apply the stochastic flux.
     double minres = 1000;
     int levelmax = 1;
     levelm = -1;
     hmin1 = hmin;
     double limhits = 0.;
     double floor_added = 0.;
     double rawmin = HUGE_VAL;
     int invalid_height = 0;
     foreach(reduction(min:minres) reduction(min:rawmin) reduction(max:levelmax) reduction(+:limhits) reduction(+:floor_added) reduction(max:invalid_height)) {
         if (!isfinite(h[])) invalid_height = 1;
         rawmin=fmin(rawmin,h[]);
         if (tcw_height_floor > 0. && h[] < tcw_height_floor) {
             if (isfinite(h[])) floor_added += (tcw_height_floor-h[])*sq(Delta);
             h[] = tcw_height_floor;
             limhits += 1.;
         }

         hi[] = h[] != 0. ? 1.0 / h[] : HUGE_VAL;
         double res = h[];
         int ll = level;

         if (res < minres) minres = res;
         if (ll > levelmax) levelmax = ll;
     }
     hmin = invalid_height ? NAN : minres;
     lm = levelmax;
     campaign_floor_volume += floor_added;
     campaign_pre_floor_minimum=fmin(campaign_pre_floor_minimum,rawmin);
     campaign_floor_cells += (unsigned long long)limhits;
     boundary({h,hi});
#if !T2OR
     if (i % campaign_log_stride == 0 || i == tcw_stop_index)
         tcw_write_film_speed(h, h_noise_old, dt, i+restart_offset, (int)(limhits + .5));
#endif
 }
//  #if TREE
// event adapt(i++) {
//     static int refined_once = 0;

//     // 1. Initial full-domain refinement.
//     if (!refined_once && hmin < 0.15) {
//         refine (level < 10);
//         refined_once = 1;
//         if (pid() == 0)
//             fprintf(stderr, "Refined whole domain to level 11 at i=%d (hmin=%g, t=%g)\n", i, hmin, t);
//     }

//     // 2. Mark core cells.
//     scalar mark[];
//     foreach() {
//         mark[] = (h[] < 0.) ? 1.0 : 0.0;
//     }
//     boundary({mark});

//     // 3. Dilate marked cells to cover a radius of 1.0.
//     // Estimate passes from the finest grid spacing.
//     // L0 = 100, Level_Max = 11 -> Delta_min ≈ 100 / 2048 ≈ 0.0488
//     // Use a few extra passes beyond radius/Delta_min.

//     // Estimate dilation passes from the maximum grid level.
//     double delta_min = L0 / (1 << LEVEL_MAX);
//     int n_expand = (int)(1.0 / delta_min) + 2;

//     // Bound the pass count to avoid pathological runtimes.
//     if (n_expand > 50) n_expand = 50;

//     // Dilate through a temporary field to avoid directional bias.
//     for (int k = 0; k < n_expand; k++) {
//         scalar mark_new[];
//         foreach() {
//             mark_new[] = mark[]; // Preserve existing marks.
//         }

//         // Mark cells adjacent to an already marked cell.
//         foreach() {
//             if (mark[] < 0.5) { // Check only unmarked cells.
//                 // Include diagonal neighbors for more isotropic dilation.
//                 if (mark[1,0] > 0.5 || mark[-1,0] > 0.5 ||
//                     mark[0,1] > 0.5 || mark[0,-1] > 0.5 ||
//                     mark[1,1] > 0.5 || mark[1,-1] > 0.5 ||
//                     mark[-1,1] > 0.5 || mark[-1,-1] > 0.5) {
//                     mark_new[] = 1.0;
//                 }
//             }
//         }

//         // Replace the old mark field with the dilated field.
//         foreach() {
//             mark[] = mark_new[];
//         }
//         boundary({mark});
//     }

//     // 4. Refine the marked region.
//     refine(mark[] > 0.5 && level < 11);
// }
// #endif

 event volume_stats (i++) {
    if (i % campaign_log_stride && i != tcw_stop_index) return 0;
    double vtot = 0.0;
    double vex  = 0.0;

    double hsum = 0.0, hsum2 = 0.0, ncell = 0.0;
    foreach (reduction(+:vtot) reduction(+:vex)
             reduction(+:hsum) reduction(+:hsum2) reduction(+:ncell)) {
      double dA = sq(Delta);   // Cell area in two dimensions.
      double hh = h[];
      vtot += hh * dA;
      hsum += hh; hsum2 += hh*hh; ncell += 1.;

      double ex = hh - HOUT;
      if (ex > 0.) vex += ex * dA;
    }

  /* 2026-09-03: no explicit MPI_Allreduce here - the foreach(reduction(+:...))
     above already sums across ranks; an extra Allreduce multiplied V_total by
     the rank count (x8 on 8 ranks). */

    V_total  = vtot;
    V_excess = vex;
    { /* 2026-09-03: rms roughness for Fig.7-style dt-convergence checks */
      double hmean = ncell > 0. ? hsum/ncell : 0.;
      double hvar  = ncell > 0. ? hsum2/ncell - hmean*hmean : 0.;
      tcw_r_rms = hvar > 0. ? sqrt(hvar) : 0.;
    }
  }

 /* Logging and output events. */

 event logWriting (i++) {
    if (i % campaign_log_stride && i != tcw_stop_index) return 0;
    if (pid() != 0) return 0.;

    const char *fname = "log_restart";

    // Write the header only if the file is missing or empty.
    static int header_checked = 0;
    if (!header_checked) {
      struct stat st;
      int need_header = (stat(fname, &st) != 0) || (st.st_size == 0);

      FILE *fp0 = fopen(fname, need_header ? "w" : "a");
      if (fp0) {
        if (need_header) {
            fprintf(fp0, "i h1 h2 t dt lm V_total V_excess r_rms\n");
        }
        fclose(fp0);
      }
      header_checked = 1;
    }

    // Append a timestamp every 100 steps.
    if (i % 100 == 0) {
      time_t current_time = time(NULL);

      // Set the timezone for consistent timestamps.
      setenv("TZ", "Asia/Shanghai", 1);
      tzset();

      struct tm *local_time = localtime(&current_time);
      char time_str[100];
      strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", local_time);

      printf("Step %d: Time %s, LevelMax = %d\n", i, time_str, lm);

      FILE *time_log = fopen(fname, "a");
      if (time_log) {
        fprintf(time_log, "# Time: Step %d: %s, lm = %d\n", i, time_str, lm);
        fclose(time_log);
      }
    }

    // Append the current data row.
    FILE *fp = fopen(fname, "a");
    if (fp) {
        fprintf(fp, "%d %g %g %g %g %d %.15g %.15g %.15g\n",
                i + 1 + restart_offset, hmin1, hmin, restart_now(), dt, lm, V_total, V_excess, tcw_r_rms);
      fclose(fp);
    }
    fp = fopen("solver_convergence.log", i == 0 ? "w" : "a");
    if (fp) {
        if (i == 0)
            fprintf(fp, "i calls failures max_picard max_mg_cycles max_picard_update max_mg_residual\n");
        fprintf(fp, "%d %ld %ld %d %d %.15g %.15g\n", i,
                tcw_solver_calls, tcw_solver_failures,
                tcw_solver_max_picard, tcw_solver_max_mg_cycles,
                tcw_solver_max_update, tcw_solver_max_residual);
        fclose(fp);
        fp=fopen("flux_limiter.log",i == 0 ? "w" : "a");
        if (!fp) {perror("flux_limiter.log");exit(10);}
        if (i == 0) fprintf(fp,"step time noise_hits noise_active noise_retained_rms capillary_hits capillary_active capillary_retained_rms\n");
        fprintf(fp,"%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",i+1+restart_offset,restart_now(),
                campaign_flux_hits[0],campaign_flux_active[0],campaign_flux_raw2[0]>0.?sqrt(campaign_flux_limited2[0]/campaign_flux_raw2[0]):1.,
                campaign_flux_hits[1],campaign_flux_active[1],campaign_flux_raw2[1]>0.?sqrt(campaign_flux_limited2[1]/campaign_flux_raw2[1]):1.);
        fclose(fp);
        if (TCW_COUPLED_STOCHASTIC && tcw_noise_mode>=0) {
            fp=fopen("joint_flux_limiter.log",i == 0 ? "w" : "a");
            if (!fp) {perror("joint_flux_limiter.log");exit(10);}
            if (i==0) fprintf(fp,"step time joint_hits joint_active joint_retained_rms raw_implicit_hmin_min\n");
            fprintf(fp,"%d %.17g %.17g %.17g %.17g %.17g\n",i+1+restart_offset,restart_now(),
                    campaign_flux_hits[2],campaign_flux_active[2],
                    campaign_flux_raw2[2]>0.?sqrt(campaign_flux_limited2[2]/campaign_flux_raw2[2]):1.,
                    campaign_joint_raw_minimum);
            fclose(fp);
        }
    }
  }


 // Legacy safe_figure time event removed: it forces t=28 and t=10000.
//  event outputH_byStep(i += 100) {
//      char filename[256];
//      snprintf(filename, sizeof(filename), "%s/h_ada_%08d.bin", dirFlowFields, i);

//      // Use the binary output helper above.
//      save_binary_data({h}, filename);
//  }
// event outputsolutionadapt(i++)
// {

// }
event outputsolutionadapt(i++ ) {
    if (i%10<0) {

      char name[80];
      sprintf (name, "dat/adasnapshot-%08d.dat", i);

      // Count cells across all threads.
      int n_cells = 0;
      foreach(reduction(+:n_cells))
        n_cells++;

        FILE * fpx = fopen(name, "w");

        // Write the file header.
        fprintf(fpx, "VARIABLES= X, Y, h \n");

        // Write cell vertices and connectivity.
        fprintf(fpx, "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n",
                n_cells*4, n_cells);

        // Serialize output so only one thread writes at a time.
        foreach(serial) {
          // The original debug print was removed.

          // Write the four cell vertices.
          fprintf(fpx, "%g\t%g\t%g\n", x - Delta/2., y - Delta/2., interpolate(h, x - Delta/2., y - Delta/2.));
          fprintf(fpx, "%g\t%g\t%g\n", x + Delta/2., y - Delta/2., interpolate(h, x + Delta/2., y - Delta/2.));
          fprintf(fpx, "%g\t%g\t%g\n", x + Delta/2., y + Delta/2., interpolate(h, x + Delta/2., y + Delta/2.));
          fprintf(fpx, "%g\t%g\t%g\n", x - Delta/2., y + Delta/2., interpolate(h, x - Delta/2., y + Delta/2.));
        }

        // Write element connectivity.
        int node = 1;
        for (int i = 0; i < n_cells; i++) {
          fprintf(fpx, "%d %d %d %d\n", node, node+1, node+2, node+3);
          node += 4;
        }

        fclose(fpx);

      }
    }

    event outputsolution(i++ ) {

        if (0)
       {
         char dump_name[80];
         sprintf(dump_name, "Lay/Layer-%08d", i);
         dump (file = dump_name);
       //   if(i==2) {
       //  dump (file = "hinit");}
           int NN= pow(2,LEVEL_MIN);
           double dl1=L0/NN;
         char name[80];
         sprintf (name, "dat/snapshot-%08d.dat", i);
         FILE *  fpx = fopen(name, "w");
         fprintf(fpx,"VARIABLES= X, Y, h \n");
         fprintf(fpx,"ZONE I=%d,J=%d, F=POINT \n",NN,NN);
         for (double y = 0.5*dl1 ; y <=L0-0.5*dl1 ; y += dl1)
         {
           for (double x = 0.5*dl1 ; x <=L0-0.5*dl1 ; x += dl1)
           {
           fprintf(fpx,"%g\t %g\t %g\n ",x,y,interpolate(h,x,y));
           }
         }
         fclose(fpx);

       }
       }


event jcp_metrics(i++) {
    if (i % campaign_log_stride && i != tcw_stop_index) return 0;
    double volume=0., excess=0., minimum=HUGE_VAL, maximum=-HUGE_VAL;
    double xc=0., xs=0., yc=0., ys=0.;
    foreach(reduction(+:volume) reduction(+:excess) reduction(min:minimum) reduction(max:maximum)
            reduction(+:xc) reduction(+:xs) reduction(+:yc) reduction(+:ys)) {
        double ex = fmax(h[] - HOUT, 0.);
        volume += h[]*sq(Delta);
        excess += ex*sq(Delta);
        minimum = fmin(minimum, h[]);
        maximum = fmax(maximum, h[]);
        double core = h[] > .05 ? ex*sq(Delta) : 0.;
        xc += core*cos(2.*pi*x/L0); xs += core*sin(2.*pi*x/L0);
        yc += core*cos(2.*pi*y/L0); ys += core*sin(2.*pi*y/L0);
    }
    double cx=atan2(xs,xc)*L0/(2.*pi), cy=atan2(ys,yc)*L0/(2.*pi);
    if (cx < 0.) cx += L0;
    if (cy < 0.) cy += L0;
    double moment=0., origin_moment=0., wet_area=0., rmax=0., edge=0.;
    foreach(reduction(+:moment) reduction(+:origin_moment) reduction(+:wet_area)
            reduction(max:rmax) reduction(+:edge)) {
        double ex = fmax(h[]-HOUT,0.)*sq(Delta);
        double rx=x-cx, ry=y-cy;
        rx -= L0*floor(rx/L0+.5); ry -= L0*floor(ry/L0+.5);
        double rr=sqrt(sq(rx)+sq(ry));
        moment += ex*sq(rr);
        origin_moment += ex*(sq(x-L0/2.)+sq(y-L0/2.));
        if (rr > .4*L0) edge += ex;
        if (h[] > .03) { wet_area += sq(Delta); rmax=fmax(rmax,rr); }
    }
    if (pid() == 0) {
        struct timespec now; clock_gettime(CLOCK_MONOTONIC,&now);
        double wall=(now.tv_sec-campaign_started.tv_sec)+1e-9*(now.tv_nsec-campaign_started.tv_nsec);
        FILE *fp = fopen("metrics.csv", i == 0 ? "w" : "a");
        if (!fp) { perror("metrics.csv"); exit(10); }
        if (i == 0) fprintf(fp, "step,t_end,dt,hmin,hmax,volume,mass_relative,mass_balance_relative,floor_added,floor_over_excess0,ell_centroid,ell_origin,cx,cy,R_area_h03,R_max_h03,edge_excess_fraction,wall_seconds,solver_failures\n");
        fprintf(fp, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%ld\n",
                i+1+restart_offset,restart_now(),dt,minimum,maximum,volume,(volume-jcp_initial_volume)/jcp_initial_volume,
                (volume-jcp_initial_volume-campaign_floor_volume)/jcp_initial_volume,
                campaign_floor_volume,campaign_floor_volume/campaign_initial_excess,
                excess > 0. ? sqrt(moment/(2.*excess)) : 0.,
                excess > 0. ? sqrt(origin_moment/excess) : 0.,cx,cy,
                sqrt(wet_area/pi),rmax,excess > 0. ? edge/excess : 0.,wall,tcw_solver_failures);
        fclose(fp);
        fp=fopen("hard_floor.log",i == 0 ? "w" : "a");
        if (!fp) {perror("hard_floor.log"); exit(10);}
        if (i == 0) fprintf(fp,"step time cumulative_added_volume cumulative_cell_hits fraction_of_initial_excess minimum_before_floor\n");
        fprintf(fp,"%d %.17g %.17g %llu %.17g %.17g\n",i+1+restart_offset,restart_now(),campaign_floor_volume,campaign_floor_cells,
                campaign_floor_volume/campaign_initial_excess,campaign_pre_floor_minimum);
        fclose(fp);
    }
}

static void campaign_save_contour (int target) {
    char name[128];
    snprintf(name,sizeof(name),"contour/hfield_t%03d.bin",target);
    write_hfield_binary(name);
    if (pid() != 0) return;
    FILE *in=fopen(name,"rb");
    if (!in) {perror(name); exit(10);}
    int n=0;
    if (fread(&n,sizeof(int),1,in) != 1 || n != (1 << LEVEL_MIN)) exit(10);
    double *a=malloc((size_t)n*n*sizeof(double));
    if (!a || fread(a,sizeof(double),(size_t)n*n,in) != (size_t)n*n) exit(10);
    fclose(in);
    double dx=L0/n;
    snprintf(name,sizeof(name),"contour/hfield_t%03d.plt",target);
    FILE *out=fopen(name,"w");
    if (!out) {perror(name);exit(10);}
    fprintf(out,"TITLE=\"1021 droplet, actual t=%.17g\"\nVARIABLES=\"X\",\"Y\",\"h\"\nZONE I=%d,J=%d,F=POINT\n",restart_now(),n,n);
    for (int iy=0;iy<n;iy++) for(int ix=0;ix<n;ix++)
        fprintf(out,"%.12g %.12g %.17g\n",(ix+.5)*dx,(iy+.5)*dx,a[(size_t)iy*n+ix]);
    fclose(out);
    snprintf(name,sizeof(name),"contour/profile_t%03d.dat",target);
    out=fopen(name,"w");
    if (!out) {perror(name);exit(10);}
    fprintf(out,"# actual_time=%.17g; s hx(y=L/2) hy(x=L/2), neighbouring cell-centre average\n",restart_now());
    for (int k=0;k<n;k++) fprintf(out,"%.12g %.17g %.17g\n",(k+.5)*dx-L0/2.,
        .5*(a[(size_t)(n/2)*n+k]+a[(size_t)(n/2-1)*n+k]),
        .5*(a[(size_t)k*n+n/2]+a[(size_t)k*n+n/2-1]));
    fclose(out); free(a);
}

event tcw_samples(i++) {
    int completed = i + 1 + restart_offset;
    if (i == 0 || completed == (int)llround(1./tcw_dt) || completed == restart_total_steps ||
        (tcw_sample_every > 0 && completed % tcw_sample_every == 0)) {
        char name[80];
        snprintf(name, sizeof(name), "tcw_t%06d.hfield", completed);
        write_hfield_binary(name);
        if (pid() == 0)
            fprintf(stderr, "TCW_SAMPLE step=%d t=%.17g file=%s\n", completed, restart_now(), name);
    }
    static int next_target=0;
    int targets[]={1,10,50,100,200,300,400,500,1000,2000};
    while (next_target < 10 && restart_offset*tcw_dt >= targets[next_target]-1e-9) next_target++;
    while (next_target < 10 && restart_now() >= targets[next_target]-1e-9) {
        campaign_save_contour(targets[next_target]); next_target++;
    }
    return 0;
}

event campaign_observations(i++) {
    int completed=i+1+restart_offset, stride=(int)floor(10./tcw_dt+.5);
    if (stride<1) stride=1;
    if (completed % stride == 0) {
        char name[128];
        snprintf(name,sizeof(name),"observations/hfield_step%09d.bin",completed);
        write_hfield_binary(name);
        if (pid() == 0) {
            FILE *fp=fopen("observations/index.csv","a");
            if (!fp) {perror("observations/index.csv");exit(10);}
            fprintf(fp,"%d,%.17g,%s,%d\n",completed,restart_now(),name,1 << LEVEL_MIN);
            fclose(fp);
        }
    }
    return 0;
}


/* Observations only. First donor-limit activation is a numerical intervention,
   not a claim that h has reached the lower bound or that rupture occurred. */
static int lower_first_h6=-1,lower_first_h5=-1,lower_recorded_first=-1;
event lower_region_diagnostics(i++) {
    int completed=i+1+restart_offset;
    int new6=lower_first_h6<0 && hmin<=1e-6;
    int new5=lower_first_h5<0 && hmin<=1e-5;
    int newlim=scan_first[2]>=0 && lower_recorded_first<0;
    if(new6) lower_first_h6=completed;
    if(new5) lower_first_h5=completed;
    if(newlim) lower_recorded_first=scan_first[2];
    int stride=(int)llround(1./tcw_dt);
    if(i!=0 && completed!=restart_total_steps && completed%stride && !new6 && !new5 && !newlim) return 0;
    double a6=0.,a5=0.,near6=0.,a5inner=0.;
    double tol=64.*2.2204460492503131e-16*jcp_hp;
    foreach(reduction(+:a6) reduction(+:a5) reduction(+:near6) reduction(+:a5inner)) {
        if(h[]<=1e-6) a6+=sq(Delta);
        if(h[]<=1e-6+tol) near6+=sq(Delta);
        if(h[]<=1e-5) {
            a5+=sq(Delta);
            if(sq(x-L0/2.)+sq(y-L0/2.)<=sq(jcp_radius)) a5inner+=sq(Delta);
        }
    }
    if(pid()==0) {
        FILE *fp=fopen("threshold_metrics.csv",i==0 ? "w" : "a");
        if(!fp) exit(61);
        if(i==0) fprintf(fp,"step,t,hmin,area_le_1e6,area_near_1e6,area_le_1e5,area_le_1e5_inside_initial_R5,domain_fraction_le_1e5,first_limiter_step,first_h_le_1e6_step,first_h_le_1e5_step,limited_faces_cumulative,removed_total_flux_rms_fraction,floor_added\n");
        double rms=campaign_flux_raw2[2]>0. ? sqrt(scan_removed2[2]/campaign_flux_raw2[2]) : 0.;
        fprintf(fp,"%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%.17g,%.17g,%.17g\n",completed,completed*tcw_dt,hmin,a6,near6,a5,a5inner,a5/sq(L0),scan_first[2],lower_first_h6,lower_first_h5,campaign_flux_hits[2],rms,campaign_floor_volume);
        fclose(fp);
        fp=fopen("threshold_latest.json.tmp","w"); if(!fp) exit(61);
        fprintf(fp,"{\"step\":%d,\"t\":%.17g,\"hmin\":%.17g,\"area_le_1e6\":%.17g,\"area_near_1e6\":%.17g,\"near_tolerance\":%.17g,\"area_le_1e5\":%.17g,\"area_le_1e5_inside_initial_R5\":%.17g,\"first_limiter_step\":%d,\"first_limiter_time\":%.17g,\"first_h_le_1e6_step\":%d,\"first_h_le_1e5_step\":%d,\"removed_total_flux_rms_fraction\":%.17g,\"floor_added\":%.17g}\n",completed,completed*tcw_dt,hmin,a6,near6,tol,a5,a5inner,scan_first[2],scan_first[2]<0 ? -1. : scan_first[2]*tcw_dt,lower_first_h6,lower_first_h5,rms,campaign_floor_volume);
        if(fclose(fp) || rename("threshold_latest.json.tmp","threshold_latest.json")) exit(61);
    }
    if(new6 || new5 || newlim) {
        char name[128];
        snprintf(name,sizeof(name),"observations/first_threshold_step%09d.bin",completed);
        write_hfield_binary(name);
    }
    return 0;
}

event numeric_guard(i++) {
    if (!isfinite(hmin) || hmin <= 0. || tcw_solver_failures > 0) {
        tcw_run_failed = 1;
        if (pid() == 0)
            fprintf(stderr,
                    "NUMERIC_GUARD hmin=%.12g t=%.12g i=%d solver_failures=%ld\n",
                    hmin, restart_now(), i+1+restart_offset, tcw_solver_failures);
        return 1;
    }
    return 0;
}

event white_checkpoints(i++) {
    int completed=i+1+restart_offset;
    if (!tcw_run_failed && restart_stride>0 &&
        (completed%restart_stride==0 || completed==restart_total_steps ||
         completed==(int)llround(1./tcw_dt) || completed==(int)llround(10./tcw_dt)))
        restart_save(completed);
    return 0;
}

event tcw_stop(i = tcw_stop_index) {
    if (pid()==0) {
        FILE *fp=fopen("limiter_diagnostics.json","w");
        if (!fp) {perror("limiter_diagnostics.json");exit(10);}
        fprintf(fp,"{\"completed_steps\":%d,\"actual_time\":%.17g,\"kinds\":[",i+1+restart_offset,restart_now());
        for (int k=0;k<3;k++) {
            fprintf(fp,"%s{\"kind\":%d,\"first_step\":%d,\"first_time\":%.17g,\"hits\":%.17g,\"active\":%.17g,\"min_donor_weight\":%.17g,\"removed_flux_rms_fraction\":%.17g}",
                    k ? "," : "",k,scan_first[k],scan_first[k]>0 ? scan_first[k]*tcw_dt : -1.,
                    campaign_flux_hits[k],campaign_flux_active[k],scan_minweight[k],
                    campaign_flux_raw2[k]>0. ? sqrt(scan_removed2[k]/campaign_flux_raw2[k]) : 0.);
        }
        fprintf(fp,"]}\n");fclose(fp);
    }

    if (pid()==0) fprintf(stderr,"OPT_PROFILE noise_s=%.17g solve_s=%.17g\n",opt_noise_seconds,opt_solve_seconds);
    if (pid() == 0)
        fprintf(stderr,
                "TCW_SOLVER_SUMMARY calls=%ld failures=%ld max_picard=%d max_mg_cycles=%d max_picard_update=%.12g max_mg_residual=%.12g\n",
                tcw_solver_calls, tcw_solver_failures,
                tcw_solver_max_picard, tcw_solver_max_mg_cycles,
                tcw_solver_max_update, tcw_solver_max_residual);
    return 1;
}

static int jcp_selftest (int kind) {
    periodic(right); periodic(top); origin(0., 0.); L0 = 12.8;
    LEVEL_MIN = LEVEL_MAX = 6; init_grid(1 << LEVEL_MIN);
    if (kind == 5) {
        tcw_noise_mode=0;
        face vector flux[];
        foreach() h[]=jcp_hp+.02+.01*cos(2.*pi*x/L0)*sin(2.*pi*y/L0);
        boundary({h});
        foreach_face() flux.x[]=100.*sin(2.*pi*x/L0)+70.*cos(2.*pi*y/L0);
        boundary((scalar *){flux});
        double initial=statsf(h).sum;
        campaign_limit_flux(h,flux,.01,0);
        foreach() {foreach_dimension() h[]+=.01*(flux.x[1]-flux.x[])/Delta;}
        stats hs=statsf(h);
        if (pid() == 0) printf("POSITIVE_FLUX_SELFTEST hmin=%.17g volume_error=%.17g hits=%.0f\n",hs.min,hs.sum-initial,campaign_flux_hits[0]);
        return hs.min<jcp_hp-1e-12 || fabs(hs.sum-initial)>1e-10 ? 28 : 0;
    }
    if (kind == 3) {
        tcw_noise_mode=0; tcw_varphi_scale=1.;
        foreach() {
            double r=hypot(x-L0/2.,y-L0/2.);
            h[]=jcp_hp+(r<jcp_radius ? pow(cos(pi*r/(2.*jcp_radius)),4.) : 0.);
        }
        boundary({h});
        prepare_stochastic_rhs(h,.001,0);
        double integral=0., maxvalue=0.; int bad=0;
        foreach(reduction(+:integral) reduction(max:maxvalue) reduction(max:bad)) {
            if (!isfinite(tcw_stochastic_rhs[])) bad=1;
            integral += tcw_stochastic_rhs[]*sq(Delta);
            maxvalue=fmax(maxvalue,fabs(tcw_stochastic_rhs[]));
        }
        if (pid() == 0) printf("FORCING_SELFTEST bad=%d integral=%.17g max=%.17g\n",bad,integral,maxvalue);
        return bad || fabs(integral)>1e-8 ? 26 : 0;
    }
    if (kind == 4) {
        scalar expected[];
        double dx=L0/(1 << LEVEL_MIN), dtime=.01;
        double kx=2.*sin(pi*8./(1 << LEVEL_MIN))/dx;
        double ky=2.*sin(pi*5./(1 << LEVEL_MIN))/dx;
        foreach() {
            h[]=1.+.1*cos(2.*pi*8.*x/L0)+.05*sin(2.*pi*5.*y/L0);
            expected[]=1.+.1*cos(2.*pi*8.*x/L0)*(1.-.5*dtime*sq(sq(kx)))/(1.+.5*dtime*sq(sq(kx)))+
                .05*sin(2.*pi*5.*y/L0)*(1.-.5*dtime*sq(sq(ky)))/(1.+.5*dtime*sq(sq(ky)));
        }
        boundary({h,expected});
        solve_thin_film(h,dtime,3.,0.,5.,0);
        double error=0.;
        foreach(reduction(max:error)) error=fmax(error,fabs(h[]-expected[]));
        if (pid() == 0) printf("IMPLICIT_FOURIER_SELFTEST error=%.17g failures=%ld\n",error,tcw_solver_failures);
        return !isfinite(error) || error>1e-8 || tcw_solver_failures ? 27 : 0;
    }
    if (kind == 2) {
        int bad = 0; double d = L0/(1 << LEVEL_MIN);
        for (int c = 0; c < 2; c++)
            for (int j = 0; j < (1 << LEVEL_MIN); j++)
                for (int k = 0; k < (1 << LEVEL_MIN); k++) {
                    double xx = (k + (c ? .5 : 0.))*d;
                    double yy = (j + (c ? 0. : .5))*d;
                    double a = tcw_face_gaussian_native(xx, yy, c, 7);
                    double b = tcw_hash_gaussian(k, j, c, 7, 1 << LEVEL_MIN, 0x243f6a8885a308d3ULL);
                    if (a != b) bad++;
                }
        if (pid() == 0) printf("FACE_INDEX_SELFTEST mismatches=%d\n", bad);
        return bad ? 22 : 0;
    }
    scalar b[], res[], reference[], lambda[];
    face vector alpha[], beta[], gamma[];
    foreach() {
        h[] = 1. + .2*cos(2.*pi*x/L0) + .1*sin(4.*pi*y/L0);
        reference[] = h[]; b[] = 0.; lambda[] = -1000.;
        if (kind == 1) h[] = NAN;
    }
    foreach_face(x) {
        alpha.x[] = 1. + .2*sin(2.*pi*y/L0);
        beta.x[] = .2 + .05*cos(2.*pi*x/L0);
        gamma.x[] = .1*sin(2.*pi*(x+y)/L0);
    }
    foreach_face(y) {
        alpha.y[] = .7 + .15*cos(2.*pi*y/L0) + .1*sin(4.*pi*x/L0);
        beta.y[] = .3 + .05*sin(2.*pi*y/L0);
        gamma.y[] = .15*cos(2.*pi*(x-y)/L0);
    }
    boundary({h, b, reference, lambda}); boundary((scalar *){alpha,beta,gamma});
    struct EquationParams p = {.D=alpha, .beta=beta, .gamma=gamma, .lambda=lambda};
    double residual = residual_thin_film({h},{b},{res},&p);
    if (kind == 1) {
        if (pid() == 0) printf("NONFINITE_SELFTEST rejected=%d\n", !isfinite(residual));
        return isfinite(residual) ? 23 : 0;
    }
    foreach() b[] = -res[];
    boundary({b});
    relax_thin_film({h},{b},LEVEL_MIN,&p);
    double error = 0.;
    foreach(reduction(max:error)) error = fmax(error, fabs(h[]-reference[]));
    if (pid() == 0) printf("OPERATOR_FIXED_POINT_SELFTEST error=%.17g\n",error);
    return error > 1e-12 ? 21 : 0;
}

static int jcp_gaussian_selftest (void) {
    L0 = 12.8; tcw_sigma = .1; tcw_ouse_level = 6;
    tcw_se_ensure();
    double emin = HUGE_VAL, emax = 0., esum = 0., fourth = 0.;
    size_t nn = (size_t)tcw_noise_n*tcw_noise_n;
    for (int r = 0; r < 64; r++) {
        tcw_se_innovation_field(0, r);
        double e = 0., mean = 0.;
        for (size_t j = 0; j < nn; j++) {
            double v = tcw_se_innovation[j]; mean += v; e += v*v; fourth += v*v*v*v;
        }
        if (fabs(mean/nn) > 1e-12) return 25;
        e /= nn; emin = fmin(emin,e); emax = fmax(emax,e); esum += e;
    }
    printf("GAUSSIAN_FIXED_NORMALIZATION_SELFTEST energy_min=%.17g energy_max=%.17g mean_energy=%.17g fourth=%.17g\n",
           emin,emax,esum/64,fourth/(64*nn));
    return (emax-emin > 1e-5 && fabs(esum/64-1.) < .02 && fabs(fourth/(64*nn)-3.) < .15) ? 0 : 24;
}

/* Exact mixed-mode stochastic-RHS response, then a strong positive update.
   This exercises the actual combined reconstruction, not only the limiter. */
static int coupled_update_selftest (void) {
    if (!TCW_COUPLED_STOCHASTIC) return 31;
    periodic(right);periodic(top);origin(0.,0.);L0=12.8;
    LEVEL_MIN=LEVEL_MAX=6;init_grid(1 << LEVEL_MIN);
    tcw_noise_mode=0;tcw_stochastic_rhs_active=1;
    double dx=L0/(1 << LEVEL_MIN),dt_test=.01;
    double kx=2.*sin(pi*8./(1 << LEVEL_MIN))/dx;
    double ky=2.*sin(pi*5./(1 << LEVEL_MIN))/dx;
    scalar expected[];
    for (int strong=0;strong<2;strong++) {
        double amplitude=strong ? 1e4 : .1;
        foreach() {
            h[]=1.;
            expected[]=1.+dt_test*amplitude*cos(2.*pi*8.*x/L0)*cos(2.*pi*5.*y/L0)/
                (1.+.5*dt_test*sq(sq(kx)+sq(ky)));
        }
        boundary({h,expected});
        foreach_face(x) tcw_unlimited_stochastic_flux.x[]=
            .5*amplitude/kx*sin(2.*pi*8.*x/L0)*cos(2.*pi*5.*y/L0);
        foreach_face(y) tcw_unlimited_stochastic_flux.y[]=
            .5*amplitude/ky*cos(2.*pi*8.*x/L0)*sin(2.*pi*5.*y/L0);
        boundary((scalar *){tcw_unlimited_stochastic_flux});
        foreach() {
            tcw_stochastic_rhs[]=0.;
            foreach_dimension() tcw_stochastic_rhs[]+=
                (tcw_unlimited_stochastic_flux.x[1]-tcw_unlimited_stochastic_flux.x[])/Delta;
        }
        boundary({tcw_stochastic_rhs});
        solve_thin_film(h,dt_test,3.,0.,5.,strong);
        double error=0.;
        foreach(reduction(max:error)) error=fmax(error,fabs(h[]-expected[]));
        stats hs=statsf(h);
        if (pid()==0) printf("COUPLED_UPDATE_SELFTEST strong=%d error=%.17g hmin=%.17g volume_error=%.17g joint_hits=%.0f failures=%ld\n",
            strong,error,hs.min,hs.sum-L0*L0,campaign_flux_hits[2],tcw_solver_failures);
        if (!isfinite(hs.min) || hs.min<white_lower-1e-12 || fabs(hs.sum-L0*L0)>1e-10 ||
            tcw_solver_failures || (!strong && error>1e-8) ||
            (strong && campaign_flux_hits[2]<=0.)) return 32;
    }
    tcw_stochastic_rhs_active=0;
    return 0;
}
