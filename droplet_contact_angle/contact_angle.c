/**
 * @file    thinfilm_stochastic.c
 * @brief   Stochastic Thin-Film Equation (STFE) solver .
 * @date    2026-01-04
 * @author  <Liang jiahua/Department of Modern Mechanics,USTC><liangjh@mail.ustc.edu.cn>
 * @version 1.0.0
 */

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
 #include <time.h>
#include <stdint.h>
 
 /* --- 全局参数与宏定义 (原 parameters.h & layer.c) --- */
 
 #define T2OR 0          // 0: 一阶时间积分, 1: 二阶时间积分
 
 #ifndef NOISE_METHOD
#define NOISE_METHOD 3
#endif
#ifndef NOISE_VARPHI
#define NOISE_VARPHI 0.
#endif
#ifndef NOISE_CORR_LEVEL
#define NOISE_CORR_LEVEL 6
#endif
#ifndef NOISE_ALL_FIELD
#define NOISE_ALL_FIELD 0
#endif
#ifndef RUPTURE_HEIGHT
#define RUPTURE_HEIGHT 1e-4
#endif
#ifndef LOG_EVERY
#define LOG_EVERY 1000
#endif
#ifndef RADIUS_EVERY
#define RADIUS_EVERY 1000
#endif
#ifndef HP001_ADAPT_EVERY
#define HP001_ADAPT_EVERY 20
#endif

// 物理系数
#ifndef COEFF_A
#define COEFF_A 0.0
#endif
#ifndef HP_VALUE
#define HP_VALUE 0.01
#endif
#ifndef HP001_TEST_KIND
#define HP001_TEST_KIND 0
#endif
#ifndef HP001_GHOST_MODE
#define HP001_GHOST_MODE 2
#endif
#ifndef HP001_TEST_GRID_LEVEL
#define HP001_TEST_GRID_LEVEL 7
#endif
#ifndef HP001_TEST_AMPLITUDE
#define HP001_TEST_AMPLITUDE (HP_VALUE*1e-3)
#endif
#ifndef HP001_TEST_TRANSPOSE
#define HP001_TEST_TRANSPOSE 0
#endif
#ifndef HP001_FORCE_COUPLED_SOLVER
#define HP001_FORCE_COUPLED_SOLVER (HP001_TEST_KIND > 0)
#endif
#ifndef HP001_OUTPUT_LEVEL
#define HP001_OUTPUT_LEVEL 10
#endif
 double coeff_varphi = NOISE_VARPHI;
#define TCW_DIAGNOSTIC_NO_DETERMINISTIC 0
#define DYNAMIC_CONTACT_LINE_AMR 0

 
 // 目录路径
 char dirFlowFields[] = "dat/";
 char dirFigures[]    = "png/";
 
 // 全局物理参数

 double n_exp = 3.0;     // 指数 n
 double m_exp = 3.0;     // 指数 m
 double dd_exp = 5.0;    // 指数 dd (分离压)
 
 // 时间步长参数
 double DTMAX = 1e-3;
 double DTMIN = 1e-4;
double tcw_dt = 1e-4;
double tcw_end_time = 500.;
int last_picard_iter = 0;
double last_picard_res = 0.;
int last_mg_i = 0;
double last_mg_resa = 0., last_mg_resb = 0.;
static double V0 = NAN, max_rel_volume_drift = 0.;

 double dtold;
 double errorstep = 0.;

 // --- 体积统计 ---
double V_total  = 0.0;   // ∫ h dA
double V_excess = 0.0;   // ∫ max(h - HOUT, 0) dA（可选）
double HOUT = HP_VALUE;
const double HP = HP_VALUE, R0 = 5., H0 = 1.,
             NOISE_ACTIVATE_HEIGHT = HP_VALUE + 0.2;
const double H_CONTACT = 1.1*HP_VALUE;
int contour_every = 1000;
/* Disable legacy contour snapshots; checkpoint_output owns the five-time schedule. */
static const double contour_targets[] = {1e99, 1e99, 1e99, 1e99, 1e99, 1e99, 1e99, 1e99};
static const int n_contour_targets = 8;
static int contour_target_next = 0;
double R_contact = 0.;
static uint64_t tcw_noise_seed = 1;       // 若你有基底厚度（hout），填对应值；没有就保持 0
 // 网格与求解器状态变量
 scalar h[], hi[], l[], h_noise_old[], noise_active[], hlog[];
 scalar h_prev[], pchem[], delta_h[], delta_p[], pi_slope[], rhs_h[], rhs_p[];
 face vector hp_mobility[];
static volatile int hp001_solver_invalid = 0;
  static double hp001_resume_base_time = NAN;
  /*
   * A restored checkpoint is written immediately before the integration event
   * at the checkpoint time.  Basilisk's restore-side event pass can consume
   * that event without applying the coupled update.  The resume branch below
   * calls this helper once to fill that one missing physical step.  The fresh
   * initialization path and the normal integration event remain unchanged.
   */
  static void hp001_coupled_step (double dt_step, int timestep);
 #if T2OR
 scalar h2[], h1[], err[];
 #endif
 
 double hmin, hmin1;
 int lm = 4;
 int levelm;
int INIT_LEVEL = 10;
int LEVEL_MAX = 13;
 int LEVEL_MIN = 4;
 
 /* --- 输出管理功能 (原 outputManager.h & myoutput.h) --- */
 
 /**
  * @brief 创建目录 (如果不存在)
  */
 
void write_hfield_binary(const char *path) {
    /* Keep the level-10 export resolution independent of the native max level. */
    volatile int export_level = HP001_OUTPUT_LEVEL;
    int n = 1 << export_level;
    double *buf = calloc((size_t)n * (size_t)n, sizeof(double));
    double dx = L0/n, pixel_area = sq(dx);
    if (!buf) {
        if (pid() == 0) fprintf(stderr, "FATAL: cannot allocate level-10 hfield (%d^2)\n", n);
        exit(EXIT_FAILURE);
    }
    /*
     * Accumulate exact leaf-cell overlap into the level-10 pixels.  The
     * grids are dyadic, but the overlap form also handles a leaf coarser
     * than the requested export level without point-sampling or refinement.
     */
    foreach(serial) {
        double x0 = x - Delta/2., x1 = x + Delta/2.;
        double y0 = y - Delta/2., y1 = y + Delta/2.;
        int ix0 = (int)floor(x0/dx), ix1 = (int)ceil(x1/dx) - 1;
        int iy0 = (int)floor(y0/dx), iy1 = (int)ceil(y1/dx) - 1;
        if (ix0 < 0 || iy0 < 0 || ix1 >= n || iy1 >= n) {
            if (pid() == 0) fprintf(stderr, "FATAL: hfield leaf outside export domain\n");
            free(buf);
            exit(EXIT_FAILURE);
        }
        for (int iy = iy0; iy <= iy1; iy++) {
            double oy = min(y1, (iy + 1)*dx) - max(y0, iy*dx);
            for (int ix = ix0; ix <= ix1; ix++) {
                double ox = min(x1, (ix + 1)*dx) - max(x0, ix*dx);
                if (ox > 0. && oy > 0.)
                    buf[(size_t)iy*n + ix] += h[]*ox*oy/pixel_area;
            }
        }
    }
    if (pid() == 0) {
        FILE *fp = fopen(path, "wb");
        if (!fp) {
            free(buf);
            fprintf(stderr, "FATAL: cannot open hfield output %s\n", path);
            exit(EXIT_FAILURE);
        }
        fwrite(&n, sizeof(int), 1, fp);
        fwrite(buf, sizeof(double), (size_t)n * (size_t)n, fp);
        fclose(fp);
    }
    free(buf);
}

void write_centerline_profile(const char *path) {
    int n = 1 << LEVEL_MAX;
    double dx = L0/n;
    if (pid() == 0) {
        FILE *fp = fopen(path, "w");
        fprintf(fp, "# s h_horizontal(y=L0/2) h_vertical(x=L0/2)\n");
        for (int j = 0; j < n; j++) {
            double q = (j + .5)*dx;
            double hx = interpolate(h, q, L0/2.);
            double hy = interpolate(h, L0/2., q);
            fprintf(fp, "%.12g %.12g %.12g\n", q - L0/2., hx, hy);
        }
        fclose(fp);
    }
}

void write_native_pressure_stats(const char *tag, int step, double state_time) {
    face vector hp_grad[];
    scalar pressure_native[];
    double body_sum = 0., body_sq = 0., body_area = 0.;
    double film_sum = 0., film_sq = 0., film_area = 0., leaves = 0.;
    double film_h_sum = 0., film_h_min = HUGE_VAL, film_h_max = -HUGE_VAL;
    double peak_excess = 0.;
    boundary({h});
    foreach(reduction(max:peak_excess))
        peak_excess = max(peak_excess, max(h[] - HP, 0.));
    foreach_face()
        hp_grad.x[] = face_gradient_x(h, 0);
    boundary((scalar *){hp_grad});
    foreach(reduction(+:body_sum) reduction(+:body_sq) reduction(+:body_area)
            reduction(+:film_sum) reduction(+:film_sq) reduction(+:film_area)
            reduction(+:leaves) reduction(+:film_h_sum)
            reduction(min:film_h_min) reduction(max:film_h_max)) {
        double pres = 0.;
        foreach_dimension()
            pres -= (hp_grad.x[1] - hp_grad.x[])/Delta;
        double hh = max(h[], 1e-12);
        double pi_h = COEFF_A*(pow(HP_VALUE, dd_exp - m_exp)*pow(hh, -dd_exp) -
                               pow(hh, -m_exp));
        double ptotal = pres - pi_h;
        pressure_native[] = ptotal;
        double area = sq(Delta), excess = max(hh - HP, 0.);
        double r = hypot(x - L0/2., y - L0/2.);
        leaves += 1.;
        if (peak_excess > 0. && excess >= 0.4*peak_excess &&
            excess <= 0.8*peak_excess) {
            body_sum += ptotal*area;
            body_sq += sq(ptotal)*area;
            body_area += area;
        }
        else if (r > 5. && r < 7. && hh < H_CONTACT) {
            film_sum += ptotal*area;
            film_sq += sq(ptotal)*area;
            film_area += area;
            film_h_sum += hh*area;
            film_h_min = min(film_h_min, hh);
            film_h_max = max(film_h_max, hh);
        }
    }
    boundary({pressure_native});
    face vector physical_flux[];
    foreach_face() {
        double hf = (h[] + h[-1])/2.;
        physical_flux.x[] = hf*hf*hf*face_gradient_x(pressure_native, 0);
    }
    boundary((scalar *){physical_flux});
    double steady_rhs_max = 0., steady_rhs_body_max = 0.;
    foreach(reduction(max:steady_rhs_max) reduction(max:steady_rhs_body_max)) {
        double rhs = 0.;
        foreach_dimension()
            rhs += (physical_flux.x[1] - physical_flux.x[])/Delta;
        steady_rhs_max = max(steady_rhs_max, fabs(rhs));
        double excess = h[] - HP;
        if (peak_excess > 0. && excess >= .4*peak_excess && excess <= .8*peak_excess)
            steady_rhs_body_max = max(steady_rhs_body_max, fabs(rhs));
    }
    if (pid() == 0) {
        static int header_written = 0;
        FILE *fp = fopen("pressure_native.csv", header_written ? "a" : "w");
        if (!fp) {
            fprintf(stderr, "FATAL: cannot open pressure_native.csv\n");
            exit(EXIT_FAILURE);
        }
        if (!header_written)
            fprintf(fp, "step,time,tag,leaf_count,peak_excess,body_area,body_mean,body_std,film_area,film_mean,film_std,film_h_mean,film_h_min,film_h_max,steady_rhs_max,steady_rhs_body_max\n");
        double body_mean = body_area > 0. ? body_sum/body_area : NAN;
        double film_mean = film_area > 0. ? film_sum/film_area : NAN;
        double body_var = body_area > 0. ? max(0., body_sq/body_area - sq(body_mean)) : NAN;
        double film_var = film_area > 0. ? max(0., film_sq/film_area - sq(film_mean)) : NAN;
        fprintf(fp, "%d,%.17g,%s,%.0f,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                step, state_time, tag, leaves, peak_excess, body_area, body_mean, sqrt(body_var),
                film_area, film_mean, sqrt(film_var),
                film_area > 0. ? film_h_sum/film_area : NAN,
                film_area > 0. ? film_h_min : NAN, film_area > 0. ? film_h_max : NAN,
                steady_rhs_max, steady_rhs_body_max);
        fclose(fp);
        header_written = 1;
    }
}

void write_snapshot(const char *tag, int step, double state_time) {
    char field[128], profile[128], meta[128];
    snprintf(field, sizeof(field), "hfield_%s.bin", tag);
    snprintf(profile, sizeof(profile), "profile_%s.dat", tag);
    snprintf(meta, sizeof(meta), "snapshot_%s.meta", tag);
    write_hfield_binary(field);
    write_centerline_profile(profile);
    write_native_pressure_stats(tag, step, state_time);
    if (pid() == 0) {
        FILE *fp = fopen(meta, "w");
        if (fp) {
            double ctheta = fmin(1., fmax(-1.,
                                          1. - COEFF_A/(4.*HP_VALUE*HP_VALUE)));
            fprintf(fp, "step=%d\ntime=%.17g\nL0=%.17g\nN=%d\nlevel=%d\n"
                        "hfield_sampling=area_weighted_native_leaf_average\n"
                        "profile_sampling=interpolated_native_maxlevel\n"
                        "hstar=%.17g\nhout=%.17g\nh_contact=%.17g\n"
                        "n_repulsive=5\nm_attractive=3\nrepulsive_coefficient=%.17g\n"
                        "A=%.17g\nk=%.17g\ngamma=1\nnoise=%.17g\n"
                        "theta_exact_degrees=%.12g\ntheta_lub_degrees=%.12g\n",
                    step, state_time, L0,
                    1 << HP001_OUTPUT_LEVEL, LEVEL_MAX,
                    (double)HP_VALUE, (double)HP_VALUE, (double)H_CONTACT,
                    pow((double)HP_VALUE, dd_exp - m_exp),
                    (double)COEFF_A, (double)COEFF_A/pow(HP_VALUE, 3.), coeff_varphi,
                    acos(ctheta)*180./pi, sqrt(COEFF_A/(2.*HP_VALUE*HP_VALUE))*180./pi);
            fclose(fp);
        }
    }
}

void update_volume_stats(void) {
    double vtot = 0., vex = 0.;
    foreach (reduction(+:vtot) reduction(+:vex)) {
        double dA = sq(Delta), excess = max(h[] - HP, 0.);
        vtot += h[]*dA;
        vex += excess*dA;
    }
    V_total = vtot;
    V_excess = vex;
    if (!isfinite(V0)) V0 = vtot;
    double rel = V0 > 0. ? fabs(vtot - V0)/V0 : NAN;
    if (isfinite(rel) && rel > max_rel_volume_drift) max_rel_volume_drift = rel;
}

int hp001_height_invalid(void) {
    int invalid = 0;
    foreach(reduction(max:invalid))
        if (!isfinite(h[]) || h[] <= 0.)
            invalid = 1;
    return invalid;
}

void update_hp001_indicator(void) {
    foreach() {
        hlog[] = log(h[]/HP_VALUE);
    }
    boundary({hlog});
}

void write_native_mesh_stats(const char *label, int step, double state_time,
                             int nf, int nc, double v_before, double v_after) {
    double leaves = 0., min_delta = HUGE_VAL, max_level = 0.;
    foreach(reduction(+:leaves) reduction(min:min_delta) reduction(max:max_level)) {
        leaves += 1.;
        min_delta = min(min_delta, Delta);
        max_level = max(max_level, (double)level);
    }
    if (pid() == 0) {
        static int header_written = 0;
        FILE *fp = fopen("amr_native.csv", header_written ? "a" : "w");
        if (!fp) {
            fprintf(stderr, "FATAL: cannot open amr_native.csv\n");
            exit(EXIT_FAILURE);
        }
        if (!header_written)
            fprintf(fp, "step,time,label,leaf_count,min_delta,max_level,refined,coarsened,V_before,V_after,rel_transfer,V0\n");
        double rel_transfer = v_before > 0. ? fabs(v_after - v_before)/v_before : NAN;
        fprintf(fp, "%d,%.17g,%s,%.0f,%.17g,%.0f,%d,%d,%.17g,%.17g,%.17g,%.17g\n",
                step, state_time, label, leaves, min_delta, max_level, nf, nc,
                v_before, v_after, rel_transfer, V0);
        fclose(fp);
        header_written = 1;
        fprintf(stderr, "HP001_MESH step=%d time=%.12g label=%s leaves=%.0f minDelta=%.12g maxLevel=%.0f nf=%d nc=%d Vbefore=%.17g Vafter=%.17g relTransfer=%.6g V0=%.17g\n",
                step, state_time, label, leaves, min_delta, max_level, nf, nc,
                v_before, v_after, rel_transfer, V0);
    }
}

void create_directory(const char *path) {
     char temp_path[256];
     strncpy(temp_path, path, 255);
     
     // 确保路径以 '/' 结尾 (为了兼容性，虽然 mkdir 不需要)
     size_t len = strlen(temp_path);
     if (len > 0 && temp_path[len - 1] != '/') {
         strncat(temp_path, "/", 255 - len);
     }
 
     struct stat st = {0};
     if (stat(temp_path, &st) == -1) {
         // 权限设置为 0700 (rwx------)
         if (mkdir(temp_path, 0700) != 0) {
             if (pid() == 0) fprintf(stderr, "Warning: Could not create directory %s\n", temp_path);
         }
     }
 }
 
 /**
  * @brief 保存二进制数据 (原 myoutput.h 中的 output_data)
  * 
  * @param list 要保存的标量列表 (如 {h})
  * @param fname 完整文件路径
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
         // MPI 并行模式下的输出逻辑
         scalar index = {-1};
         index = new scalar;
         z_indexing(index, true);
 
         // 计算每个单元的数据大小
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
         // 串行模式下的输出逻辑
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
 
 /* --- 求解器结构体与核心函数 (原 layer.c) --- */
 
struct CoupledParams {
    double dt;
    face vector mobility;
    scalar pi_prime;
};

static double hp001_pi (double height) {
    const double hstar = pow(HP_VALUE, dd_exp - m_exp);
    return COEFF_A*(hstar*pow(height, -dd_exp) - pow(height, -m_exp));
}
static double hp001_pi_prime (double height) {
    const double hstar = pow(HP_VALUE, dd_exp - m_exp);
    return COEFF_A*(-dd_exp*hstar*pow(height, -dd_exp - 1.) +
                    m_exp*pow(height, -m_exp - 1.));
}

/* Correction operator A1=dh/dt-div(M grad(dp)); A2=dp+Pi'(hi)dh+lap(dh). */
static void hp001_mark_solver_invalid (void) {
    __sync_lock_test_and_set(&hp001_solver_invalid, 1);
}
static void hp001_clear_solver_invalid (void) {
    __sync_lock_test_and_set(&hp001_solver_invalid, 0);
}
static int hp001_any_solver_invalid (void) {
    int local = hp001_solver_invalid, global = local;
#if _MPI
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    return global;
}

/* Correction operator A1=dh/dt-div(M grad(dp)); A2=dp+Pi'(hi)dh+lap(dh). */
static double residual_hp001_coupled (scalar *al, scalar *bl,
                                      scalar *resl, void *data) {
    scalar dh = al[0], dp = al[1];
    scalar rh = bl[0], rp = bl[1];
    scalar r_h = resl[0], r_p = resl[1];
    struct CoupledParams *cp = (struct CoupledParams *) data;
    (const) face vector M = cp->mobility;
    (const) scalar Pi1 = cp->pi_prime;
    face vector gp[], flux[], gh[];
    foreach_face() {
        gp.x[] = face_gradient_x(dp, 0);
        flux.x[] = M.x[]*gp.x[];
        gh.x[] = face_gradient_x(dh, 0);
    }
    boundary((scalar *){gp, flux, gh});
    double maxres = 0.;
    int bad = 0;
    foreach(reduction(max:maxres) reduction(max:bad)) {
        double div_flux = 0., lap_dh = 0.;
        foreach_dimension() {
            div_flux += (flux.x[1] - flux.x[])/Delta;
            lap_dh += (gh.x[1] - gh.x[])/Delta;
        }
        double rhv = rh[] - dh[]/cp->dt + div_flux;
        double rpv = rp[] - dp[] - Pi1[]*dh[] - lap_dh;
        if (!isfinite(rhv) || !isfinite(rpv)) {
            bad = 1;
            r_h[] = r_p[] = 0.;
        }
        else {
            r_h[] = rhv;
            r_p[] = rpv;
            maxres = max(maxres, max(fabs(rhv), fabs(rpv)));
        }
    }
    boundary(resl);
    if (bad) hp001_mark_solver_invalid();
    return bad ? HUGE_VAL : maxres;
}

/* Pointwise 2x2 block Gauss-Seidel relaxation on each multigrid level. */
static void relax_hp001_coupled (scalar *al, scalar *bl, int lev,
                                 void *data) {
    scalar dh = al[0], dp = al[1];
    scalar rh = bl[0], rp = bl[1];
    struct CoupledParams *cp = (struct CoupledParams *) data;
    (const) face vector M = cp->mobility;
    (const) scalar Pi1 = cp->pi_prime;
    foreach_level_or_leaf(lev) {
        double invdx2 = 1./sq(Delta);
        double sumM = M.x[] + M.x[1] + M.y[] + M.y[0,1];
        double b = sumM*invdx2;
        double a = 1./cp->dt;
        double c = Pi1[] - 4.*invdx2;
        double d = 1.;
        double rhs1 = rh[] +
          (M.x[1]*dp[1] + M.x[]*dp[-1] +
           M.y[0,1]*dp[0,1] + M.y[]*dp[0,-1])*invdx2;
        double rhs2 = rp[] -
          (dh[1] + dh[-1] + dh[0,1] + dh[0,-1])*invdx2;
        double det = a*d - b*c;
        if (!isfinite(det) || det == 0. || !isfinite(rhs1) || !isfinite(rhs2)) {
            hp001_mark_solver_invalid();
            continue;
        }
        double new_dh = (d*rhs1 - b*rhs2)/det;
        double new_dp = (a*rhs2 - c*rhs1)/det;
        if (!isfinite(new_dh) || !isfinite(new_dp)) {
            hp001_mark_solver_invalid();
            continue;
        }
        dh[] = new_dh;
        dp[] = new_dp;
    }
}
static void hp001_fail_with_dump (const char *reason, int step, int picard) {
    dump(file = "native_invalid");
    if (pid() == 0)
        fprintf(stderr, "FATAL: coupled HP001 solve failed: %s step=%d picard=%d; native_invalid saved\n",
                reason, step, picard);
    exit(EXIT_FAILURE);
}

static void hp001_actual_residuals (scalar f, scalar hold, double dt,
                                    double *mass_max, double *pressure_max) {
    face vector M[], gp[], flux[], gh[];
    foreach_face() {
        double hf = (f[] + f[-1])/2.;
        if (!isfinite(hf) || hf <= 0.) {
            hp001_mark_solver_invalid();
            M.x[] = gp.x[] = flux.x[] = gh.x[] = 0.;
        }
        else {
            M.x[] = hf*hf*hf;
            gp.x[] = face_gradient_x(pchem, 0);
            flux.x[] = M.x[]*gp.x[];
            gh.x[] = face_gradient_x(f, 0);
        }
    }
    boundary((scalar *){M, gp, flux, gh});
    if (hp001_any_solver_invalid())
        hp001_fail_with_dump("invalid face state in residual audit", iter, -1);
    double maxm = 0., maxp = 0.;
    int bad = 0;
    foreach(reduction(max:maxm) reduction(max:maxp) reduction(max:bad)) {
        double div_flux = 0., lap_h = 0.;
        foreach_dimension() {
            div_flux += (flux.x[1] - flux.x[])/Delta;
            lap_h += (gh.x[1] - gh.x[])/Delta;
        }
        double rm = hold[]/dt - f[]/dt + div_flux;
        double rp = -hp001_pi(f[]) - pchem[] - lap_h;
        if (!isfinite(rm) || !isfinite(rp)) bad = 1;
        else {
            maxm = max(maxm, fabs(rm));
            maxp = max(maxp, fabs(rp));
        }
    }
    if (bad || !isfinite(maxm) || !isfinite(maxp))
        hp001_fail_with_dump("nonfinite physical residual", iter, -1);
    *mass_max = maxm;
    *pressure_max = maxp;
}

trace
mgstats solve_thin_film (scalar f, double dt_step, double m_val,
                         double n_val, double dd_val, int timestep) {
    (void) m_val; (void) n_val; (void) dd_val; (void) timestep;
    mgstats s = {0};
    if (!(dt_step > 0.) || !isfinite(dt_step))
        hp001_fail_with_dump("invalid timestep", iter, 0);
    int invalid = 0;
    foreach(reduction(max:invalid))
        if (!isfinite(f[]) || f[] <= 0.) invalid = 1;
    if (invalid) hp001_fail_with_dump("invalid input height", iter, 0);
    foreach() {
        h_prev[] = f[];
        hi[] = f[];
    }
    boundary({h_prev, hi});

    /* Initialize p on the nonlinear chemical-potential manifold. */
    face vector gh0[];
    foreach_face()
        gh0.x[] = face_gradient_x(f, 0);
    boundary((scalar *){gh0});
    foreach() {
        double lap_h = 0.;
        foreach_dimension()
            lap_h += (gh0.x[1] - gh0.x[])/Delta;
        pchem[] = -lap_h - hp001_pi(f[]);
    }
    boundary({pchem});

    const int picard_limit = 30;
    const double nonlinear_tol = 1e-8;
    double max_dh = HUGE_VAL, mass_res = HUGE_VAL, pressure_res = HUGE_VAL;
    int picard = 0;
    struct CoupledParams cp;
    cp.dt = dt_step;
    cp.mobility = hp_mobility;
    face vector mobility = cp.mobility;
    cp.pi_prime = pi_slope;

    for (picard = 1; picard <= picard_limit; picard++) {
        invalid = 0;
        foreach(reduction(max:invalid)) {
            if (!isfinite(hi[]) || hi[] <= 0.) {
                invalid = 1;
                pi_slope[] = 0.;
            }
            else {
                pi_slope[] = hp001_pi_prime(hi[]);
                if (!isfinite(pi_slope[])) {
                    invalid = 1;
                    pi_slope[] = 0.;
                }
            }
        }
        if (invalid) hp001_fail_with_dump("invalid Picard height", iter, picard);
        hp001_clear_solver_invalid();
        foreach_face() {
            double hf = (hi[] + hi[-1])/2.;
            if (!isfinite(hf) || hf <= 0.) {
                hp001_mark_solver_invalid();
                mobility.x[] = 0.;
            }
            else mobility.x[] = hf*hf*hf;
        }
        if (hp001_any_solver_invalid()) hp001_fail_with_dump("invalid Picard face height", iter, picard);
        boundary({pi_slope});
        boundary((scalar *){mobility});
        restriction({pi_slope, mobility});

        face vector gp[], flux[], gh[];
        foreach_face() {
            gp.x[] = face_gradient_x(pchem, 0);
            flux.x[] = mobility.x[]*gp.x[];
            gh.x[] = face_gradient_x(f, 0);
        }
        boundary((scalar *){gp, flux, gh});
        foreach() {
            double div_flux = 0., lap_h = 0.;
            foreach_dimension() {
                div_flux += (flux.x[1] - flux.x[])/Delta;
                lap_h += (gh.x[1] - gh.x[])/Delta;
            }
            rhs_h[] = h_prev[]/dt_step - f[]/dt_step + div_flux;
            rhs_p[] = -hp001_pi(hi[]) + pi_slope[]*hi[] -
                  pchem[] - pi_slope[]*f[] - lap_h;
        }
        boundary({rhs_h, rhs_p});
        restriction({rhs_h, rhs_p});
        delta_h.refine = refine_linear;
        delta_h.gradient = minmod2;
        set_prolongation(delta_h, HP001_GHOST_MODE ? refine_bilinear : refine_linear);
        set_restriction(delta_h, restriction_volume_average);
        delta_p.refine = refine_linear;
        delta_p.gradient = minmod2;
        set_prolongation(delta_p, HP001_GHOST_MODE ? refine_bilinear : refine_linear);
        set_restriction(delta_p, restriction_volume_average);
        hp001_clear_solver_invalid();
        foreach() {
            delta_h[] = 0.;
            delta_p[] = 0.;
        }
        boundary({delta_h, delta_p});
        s = mg_solve({delta_h, delta_p}, {rhs_h, rhs_p},
                     residual_hp001_coupled, relax_hp001_coupled, &cp);
        if (hp001_any_solver_invalid() || !isfinite(s.resa) || !isfinite(s.resb))
            hp001_fail_with_dump("nonfinite multigrid residual", iter, picard);

        max_dh = 0.;
        invalid = 0;
        foreach(reduction(max:max_dh) reduction(max:invalid)) {
            double candidate_h = f[] + delta_h[];
            double candidate_p = pchem[] + delta_p[];
            if (!isfinite(candidate_h) || candidate_h <= 0. || !isfinite(candidate_p) || !isfinite(delta_h[])) {
                invalid = 1;
                continue;
            }
            f[] = candidate_h;
            pchem[] = candidate_p;
            max_dh = max(max_dh, fabs(delta_h[]));
        }
        if (invalid) hp001_fail_with_dump("invalid coupled update", iter, picard);
        boundary({f, pchem});
        foreach() hi[] = f[];
        boundary({hi});

        hp001_actual_residuals(f, h_prev, dt_step, &mass_res, &pressure_res);
        if (pid() == 0 && timestep < 2)
            fprintf(stderr, "HP001_PICARD step=%d iteration=%d dh=%.17g mass=%.17g pressure=%.17g mg_after=%.17g mg_before=%.17g\n",
                    timestep, picard, max_dh, mass_res, pressure_res, s.resa, s.resb);
        if (!isfinite(max_dh) || !isfinite(mass_res) || !isfinite(pressure_res))
            hp001_fail_with_dump("nonfinite nonlinear residual", iter, picard);
        if (max_dh <= nonlinear_tol && mass_res <= nonlinear_tol &&
            pressure_res <= nonlinear_tol)
            break;
    }
    last_picard_iter = min(picard, picard_limit);
    last_picard_res = max(max_dh, max(mass_res, pressure_res));
    last_mg_i = s.i;
    last_mg_resa = s.resa;
    last_mg_resb = s.resb;
    if (picard > picard_limit)
        hp001_fail_with_dump("Picard/nonlinear criteria not met in 30 iterations", iter, picard_limit);
    if (pid() == 0)
        fprintf(stderr, "HP001_COUPLED step=%d picard=%d dh=%.3e mass=%.3e pressure=%.3e mg=(%.3e,%.3e)\n",
                timestep, last_picard_iter, max_dh, mass_res, pressure_res, s.resa, s.resb);
    return s;
}

/* h/p variables remain on the same adaptive tree; output and checkpoint events
 * from the base source are retained unchanged. */ double adjust_timestep(double h1, double h2, double dt_curr, int lm_curr, double dtmax_curr, int step_idx) {
     // 检查数值不稳定性
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
 
     // 寻找全局最小高度
     double h_min_local = h1;
     double h_min_global = h1;
 
 #if _MPI
     MPI_Allreduce(&h_min_local, &h_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
 #endif
 
     // 经验公式调整 dtmax
     if (h_min_global < 0.5) {
         dtmax_curr = 0.000134 * pow(55.901699, h_min_global);
     }
 
     if (pid() == 0) {
         printf("Status: h_min=%g dt=%g dtmax=%g lm=%d t=%g step=%08d\n", 
                h_min_global, dt_curr, dtmax_curr, lm_curr, t, step_idx);
     }
 
     return dtmax_curr;
 }
 
 /* --- 主程序 --- */
 
 int main(int argc, char *argv[]) {
     unsigned int seed = argc > 1 ? (unsigned int)strtoul(argv[1], NULL, 10) : 1u;
     tcw_noise_seed = seed;
     if (argc > 2) tcw_dt = strtod(argv[2], NULL);
     if (argc > 3) tcw_end_time = strtod(argv[3], NULL);
     if (argc > 4) coeff_varphi = strtod(argv[4], NULL);
     if (argc > 5) contour_every = (int)strtol(argv[5], NULL, 10);
    if (argc > 6) LEVEL_MAX = (int)strtol(argv[6], NULL, 10);
    if (!(tcw_dt > 0.) || !(tcw_end_time > 0.) || contour_every < 1 ||
        LEVEL_MAX < INIT_LEVEL || LEVEL_MAX > 14 || coeff_varphi != 0.) return 2;
    TOLERANCE = 1e-8;
    contour_every = INT_MAX;
     periodic(right);
     periodic(top);
     origin(0., 0.);
    if (HP001_TEST_KIND > 0) {
        INIT_LEVEL = HP001_TEST_GRID_LEVEL;
        LEVEL_MIN = LEVEL_MAX = HP001_TEST_GRID_LEVEL;
    }
    init_grid(1 << INIT_LEVEL);
     L0 = 32.;
     run();
 }
 
static int restored_global_step(void) { return iter; }
static double restored_global_time(void) { return t; }

 event init(i = 0) {
     create_directory(dirFlowFields);
     create_directory(dirFigures);
     create_directory("contour");

    /* Keep h conservative under every tree refine/restrict operation. */
    h.refine = refine_linear;
    h.gradient = minmod2;
    /* Separate conservative native refinement from correction halos. */
    set_prolongation(h, HP001_GHOST_MODE >= 2 ? refine_bilinear : refine_linear);
    if (HP001_GHOST_MODE) {
        hi.gradient = minmod2;
        set_prolongation(hi, HP001_GHOST_MODE >= 2 ? refine_bilinear : refine_linear);
        set_restriction(hi, restriction_volume_average);
    }
    set_prolongation(pchem, refine_bilinear);
    set_restriction(h, restriction_volume_average);
     
     const char *restart_path = getenv("RESTART_FILE");
    int has_restart = restart_path && *restart_path;
    int restored = has_restart ? restore(file = restart_path) : 0;
    if (has_restart && !restored) {
        if (pid() == 0) fprintf(stderr, "FATAL: checkpoint restore failed: %s\n", restart_path);
        exit(EXIT_FAILURE);
    }
    if (!has_restart) {
        if (pid() == 0) printf("Initializing simulation...\n");
        hmin = HP;
        hmin1 = HP + 1e-9;
        foreach() {
            double r = hypot(x - L0/2., y - L0/2.);
#if HP001_TEST_KIND == 1
            h[] = HP_VALUE;
#elif HP001_TEST_KIND == 2
            h[] = HP_VALUE + HP001_TEST_AMPLITUDE*cos(16.*pi*(x - L0/2.)/L0);
#elif HP001_TEST_KIND == 3
            /* Unequal wavelengths expose x/y mistakes with variable mobility. */
            double u = HP001_TEST_TRANSPOSE ? y : x;
            double v = HP001_TEST_TRANSPOSE ? x : y;
            h[] = HP_VALUE + HP001_TEST_AMPLITUDE*
              cos(16.*pi*(u - L0/2.)/L0)*cos(8.*pi*(v - L0/2.)/L0);
#else
            h[] = HP + (r < R0 ? H0*pow(cos(pi*r/(2.*R0)), 4.) : 0.);
#endif
            noise_active[] = h[] > NOISE_ACTIVATE_HEIGHT;
        }
        boundary({h, noise_active});
        if (hp001_height_invalid()) {
            dump(file = "native_invalid");
            if (pid() == 0)
                fprintf(stderr, "FATAL: non-positive or nonfinite h; native_invalid saved\n");
            exit(EXIT_FAILURE);
        }
        /* Start from the common level-10 analytic field, then let the dual
         * indicator coarsen the bulk and refine the contact region.  hlog is
         * mesh-only and is never fed back into the physics. */
        update_volume_stats();
#if HP001_TEST_KIND == 0
        for (int pass = 0; pass < 8; pass++) {
            update_hp001_indicator();
            double v_before = V_total;
            astats ast = adapt_wavelet({h, hlog}, {1e-4, 1e-3}, LEVEL_MAX, LEVEL_MIN,
                                       list = {h, hlog, noise_active});
            boundary({h, noise_active});
            update_volume_stats();
            write_native_mesh_stats("init", pass, 0., ast.nf, ast.nc,
                                    v_before, V_total);
            if (!ast.nf && !ast.nc)
                break;
        }
#else
        write_native_mesh_stats("test_uniform", 0, 0., 0, 0, V_total, V_total);
#endif
        boundary({h, noise_active});
        update_volume_stats();
        update_hp001_indicator();
        write_snapshot("t0000", 0, 0.);
        dump(file = "native_initial");
        if (pid() == 0) {
            FILE *fp = fopen("volume.csv", "w");
            FILE *state = fopen("volume_state.txt", "w");
            if (fp) {
                fprintf(fp, "step,time,V_total,V_excess,rel_drift,max_rel_drift,picard_iter,picard_res,mg_i,mg_resa_after,mg_resb_before\n");
                fprintf(fp, "0,%.12g,%.15g,%.15g,0,0,%d,%.8g,%d,%.8g,%.8g\n",
                        0., V_total, V_excess, last_picard_iter, last_picard_res,
                        last_mg_i, last_mg_resa, last_mg_resb);
                fclose(fp);
            }
            if (state) {
                fprintf(state, "V0=%.17g\nmax_rel_volume_drift=%.17g\n", V0, max_rel_volume_drift);
                fclose(state);
            }
        }
    } else {
        FILE *state = fopen("volume_state.txt", "r");
        int state_ok = state && fscanf(state, "V0=%lf\nmax_rel_volume_drift=%lf", &V0, &max_rel_volume_drift) == 2 &&
                       isfinite(V0) && V0 > 0. && isfinite(max_rel_volume_drift) && max_rel_volume_drift >= 0.;
        if (state) fclose(state);
        if (!state_ok) {
            if (pid() == 0) fprintf(stderr, "FATAL: invalid or missing volume_state.txt for checkpoint %s\n", restart_path);
            exit(EXIT_FAILURE);
        }
        update_volume_stats();
        if (hp001_height_invalid()) {
            dump(file = "native_invalid");
            if (pid() == 0)
                fprintf(stderr, "FATAL: non-positive or nonfinite h; native_invalid saved\n");
            exit(EXIT_FAILURE);
        }
        int restored_step = restored_global_step();
        double restored_time = restored_global_time();
        double restored_hmin = HUGE_VAL;
        int restored_lm = 0;
        foreach(reduction(min:restored_hmin) reduction(max:restored_lm)) {
            if (h[] < restored_hmin) restored_hmin = h[];
            if (level > restored_lm) restored_lm = level;
        }
        hmin = hmin1 = restored_hmin;
        lm = restored_lm;
        while (contour_target_next < n_contour_targets &&
               contour_targets[contour_target_next] <= restored_time + 1e-10)
            contour_target_next++;
        update_hp001_indicator();
        write_snapshot("resume", restored_step, restored_time);
        if (pid() == 0) fprintf(stderr, "RESUME step=%d time=%.17g file=%s V0=%.17g hmin=%.17g lm=%d\n",
                                restored_step, restored_time, restart_path, V0, hmin, lm);
        /*
         * restore() runs the event scheduler once.  For the native HP001
         * checkpoints used here that pass advances the integration clock from
         * (restored_step, restored_time) without applying the coupled update.
         * Fill exactly that missing first-order step before normal scheduling
         * resumes.  The state is not renormalized and no clock/index is
         * advanced here; the event scheduler owns both after init returns.
         */
        double resume_dt = min(tcw_dt, tcw_end_time - restored_time);
        if (restored_time < tcw_end_time && resume_dt > 0. && isfinite(resume_dt)) {
            hp001_resume_base_time = restored_time;
            if (pid() == 0)
                fprintf(stderr, "HP001_RESTORE_COUPLED_BEGIN step=%d time=%.17g dt=%.17g\n",
                        restored_step, restored_time, resume_dt);
            hp001_coupled_step(resume_dt, restored_step);
        }
     }
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
static int64_t tcw_periodic_noise_index (int64_t i, int64_t n) {
    i %= n;
    return i < 0 ? i + n : i;
}
static double tcw_noise_corr_delta (void) {
    return L0/(double)(1 << NOISE_CORR_LEVEL);
}
static double tcw_coarse_gaussian_node (int64_t ix, int64_t iy,
                                         unsigned component, int timestep) {
    const int64_t n = (int64_t)1 << NOISE_CORR_LEVEL;
    ix = tcw_periodic_noise_index(ix, n);
    iy = tcw_periodic_noise_index(iy, n);
    uint64_t key = tcw_noise_seed ^ ((uint64_t)ix*0x9e3779b97f4a7c15ULL) ^
                   ((uint64_t)iy*0xbf58476d1ce4e5b9ULL) ^
                   ((uint64_t)component*0x94d049bb133111ebULL) ^
                   ((uint64_t)timestep*0xd6e8feb86659fd93ULL);
    double u1 = tcw_u01(key), u2 = tcw_u01(key ^ 0xa0761d6478bd642fULL);
    return sqrt(-2.*log(u1))*cos(2.*pi*u2);
}
/* L6 Gaussian nodes, bilinearly interpolated to the adaptive face grid. */
static double tcw_correlated_face_gaussian (double xx, double yy,
                                             unsigned component, int timestep) {
    double ell = tcw_noise_corr_delta();
    double gx = xx/ell, gy = yy/ell;
    int64_t ix = (int64_t)floor(gx), iy = (int64_t)floor(gy);
    double tx = gx - ix, ty = gy - iy;
    double g00 = tcw_coarse_gaussian_node(ix,   iy,   component, timestep);
    double g10 = tcw_coarse_gaussian_node(ix+1, iy,   component, timestep);
    double g01 = tcw_coarse_gaussian_node(ix,   iy+1, component, timestep);
    double g11 = tcw_coarse_gaussian_node(ix+1, iy+1, component, timestep);
    return (1. - tx)*(1. - ty)*g00 + tx*(1. - ty)*g10 +
           (1. - tx)*ty*g01 + tx*ty*g11;
}
static double smooth01 (double q) {
    q = q < 0. ? 0. : (q > 1. ? 1. : q);
    return q*q*(3. - 2.*q);
}
static double tcw_noise_weight (double hh) {
#if NOISE_ALL_FIELD
    return 1.;
#else
#if NOISE_METHOD == 1
    return smooth01((hh - HP - 0.06)/0.08);
#elif NOISE_METHOD == 2
    return smooth01((hh - HP - 0.05)/0.06)*
           (1. - smooth01((hh - HP - 0.70)/0.20));
#elif NOISE_METHOD == 3
    /* Whole-drop noise; only the precursor-film interface is tapered. */
    return smooth01((hh - HP)/0.005);
#else
#error "NOISE_METHOD must be 1 (core), 2 (rim), or 3 (whole drop)"
#endif
#endif
}
void apply_split_stochastic_flux (scalar f, scalar f_old, double dt_step,
                                  int timestep) {
    face vector flux[];
    foreach_face(x) {
        double ff = max((f_old[] + f_old[-1])/2., HP);
        flux.x[] = sqrt(2.*coeff_varphi)*pow(ff, 1.5)*tcw_noise_weight(ff)*
          tcw_correlated_face_gaussian(x, y, 0, timestep)/(tcw_noise_corr_delta()*sqrt(dt_step));
    }
    foreach_face(y) {
        double ff = max((f_old[] + f_old[0,-1])/2., HP);
        flux.y[] = sqrt(2.*coeff_varphi)*pow(ff, 1.5)*tcw_noise_weight(ff)*
          tcw_correlated_face_gaussian(x, y, 1, timestep)/(tcw_noise_corr_delta()*sqrt(dt_step));
    }
    boundary((scalar *){flux});
    scalar alpha[];
    foreach() {
        /* Sum only faces that remove liquid from this cell. */
        double out = 0.;
        if (flux.x[] > 0.) out += flux.x[];
        if (flux.x[1] < 0.) out -= flux.x[1];
        if (flux.y[] > 0.) out += flux.y[];
        if (flux.y[0,1] < 0.) out -= flux.y[0,1];
        double cap = max(f[] - HP, 0.)*Delta/dt_step;
        alpha[] = out > 0. ? min(1., cap/out) : 1.;
    }
    boundary({alpha});
    foreach_face(x)
        flux.x[] *= flux.x[] > 0. ? alpha[] : alpha[-1];
    foreach_face(y)
        flux.y[] *= flux.y[] > 0. ? alpha[] : alpha[0,-1];
    boundary((scalar *){flux});
    foreach() foreach_dimension() f[] += dt_step*(flux.x[1] - flux.x[])/Delta;
    boundary({f});
}

/*
 * Apply the one first-order HP001 update that can be skipped by the
 * restore-side events(false) pass.  This is deliberately a separate helper:
 * the fresh-run integration event below is left byte-for-byte unchanged, so
 * this repair cannot alter fresh initialization or normal time stepping.
 */
static void hp001_resume_fail_with_dump (const char *reason, int step,
                                          double state_time, double dt_step) {
    dump(file = "native_invalid");
    if (pid() == 0) {
        FILE *failure = fopen("RESUME_SOLVER_FAILURE.txt", "w");
        if (failure) {
            fprintf(failure,
                    "reason=%s\nstep=%d\nstate_time=%.17g\ndt=%.17g\n"
                    "picard_iter=%d\npicard_res=%.17g\n"
                    "mg_i=%d\nmg_resa_after=%.17g\nmg_resb_before=%.17g\n",
                    reason, step, state_time, dt_step, last_picard_iter,
                    last_picard_res, last_mg_i, last_mg_resa, last_mg_resb);
            fclose(failure);
        }
        fprintf(stderr,
                "FATAL: resume HP001 coupled step failed: %s step=%d "
                "state_time=%.12g dt=%.12g picard_iter=%d picard_res=%.8g "
                "mg_i=%d mg_resa_after=%.8g mg_resb_before=%.8g; "
                "native_invalid saved\n",
                reason, step, state_time, dt_step, last_picard_iter,
                last_picard_res, last_mg_i, last_mg_resa, last_mg_resb);
    }
    exit(EXIT_FAILURE);
}

static void hp001_coupled_step (double dt_step, int timestep) {
    double resume_state_time = isfinite(hp001_resume_base_time) ?
                               hp001_resume_base_time + dt_step : t + dt_step;
    if (!(dt_step > 0.) || !isfinite(dt_step))
        hp001_resume_fail_with_dump("invalid resume timestep", timestep,
                                    resume_state_time, dt_step);

#if T2OR
    /* The source is fixed to the first-order production branch (T2OR=0). */
    hp001_resume_fail_with_dump("resume fix requires first-order branch", timestep,
                                resume_state_time, dt_step);
#else
    dtold = dt_step;
    foreach() h_noise_old[] = h[];
    boundary({h_noise_old});
#if !TCW_DIAGNOSTIC_NO_DETERMINISTIC
    double h_lo = HUGE_VAL, h_hi = -HUGE_VAL;
    foreach(reduction(min:h_lo) reduction(max:h_hi)) {
        if (h[] < h_lo) h_lo = h[];
        if (h[] > h_hi) h_hi = h[];
    }
    last_picard_iter = 0;
    last_picard_res = 0.;
    last_mg_i = 0;
    last_mg_resa = 0.;
    last_mg_resb = 0.;
    if (h_hi - h_lo > 1e-14 || HP001_FORCE_COUPLED_SOLVER)
        solve_thin_film(h, dt_step, m_exp, n_exp, dd_exp, timestep);
#endif
    apply_split_stochastic_flux(h, h_noise_old, dt_step, timestep);
#endif

    double minres = 1000.;
    int levelmax = 1, rupture = 0;
    levelm = -1;
    hmin1 = hmin;
    foreach(reduction(min:minres) reduction(max:levelmax) reduction(max:rupture)) {
        /* No hard floor: clipping would inject mass into the restored state. */
        if (!isfinite(h[]) || h[] < RUPTURE_HEIGHT) rupture = 1;
        hi[] = 1.0/max(h[], 1e-12);
        double res = fabs(h[]);
        int ll = level;
        if (res < minres) minres = res;
        if (ll > levelmax) levelmax = ll;
    }
    hmin = minres;
    lm = levelmax;
    int solver_bad = !isfinite(last_picard_res) || last_picard_res > 1e-8 ||
                     !isfinite(last_mg_resa) || last_mg_resa > TOLERANCE;
    if (solver_bad)
        hp001_resume_fail_with_dump("coupled residual guard", timestep,
                                    resume_state_time, dt_step);
    if (rupture)
        hp001_resume_fail_with_dump("rupture threshold", timestep,
                                    resume_state_time, dt_step);
    boundary({h, hi});
    if (hp001_height_invalid())
        hp001_resume_fail_with_dump("non-positive or nonfinite restored update",
                                    timestep, resume_state_time, dt_step);
    if (pid() == 0)
        fprintf(stderr,
                "HP001_RESTORE_COUPLED_DONE step=%d state_time=%.17g dt=%.17g "
                "hmin=%.17g level=%d picard_iter=%d picard_res=%.8g "
                "mg_i=%d mg_resa_after=%.8g mg_resb_before=%.8g\n",
                timestep, resume_state_time, dt_step, hmin, lm, last_picard_iter,
                last_picard_res, last_mg_i, last_mg_resa, last_mg_resb);
}

event stability(i++) {
    dt = dtnext(tcw_dt);
}

#if DYNAMIC_CONTACT_LINE_AMR
event contact_line_refine(i++) {
    for (int pass = LEVEL_MIN; pass < LEVEL_MAX; pass++) {
        refine (h[] > HP + 1e-3 && h[] < HP + .15 && level < LEVEL_MAX);
        boundary({h, noise_active});
    }
}
#endif

event hp001_adapt (i += HP001_ADAPT_EVERY) {
#if HP001_TEST_KIND > 0
    return 0;
#endif
    if (i == 0)
        return 0;
    update_hp001_indicator();
    update_volume_stats();
    double v_before = V_total;
    astats ast = adapt_wavelet({h, hlog}, {1e-4, 1e-3}, LEVEL_MAX, LEVEL_MIN,
                               list = {h, hlog, noise_active});
    boundary({h, noise_active});
    update_volume_stats();
    write_native_mesh_stats("dynamic", i, t, ast.nf, ast.nc,
                            v_before, V_total);
    return 0;
}

//  #if TREE
// event adapt(i++) {
//     static int refined_once = 0;
    
//     // 1. 初始阶段全域加密逻辑
//     if (!refined_once && hmin < 0.15) {
//         refine (level < 10);
//         refined_once = 1;
//         if (pid() == 0)
//             fprintf(stderr, "Refined whole domain to level 11 at i=%d (hmin=%g, t=%g)\n", i, hmin, t);
//     }

//     // 2. 标记核心坏点
//     scalar mark[];
//     foreach() {
//         mark[] = (h[] < 0.) ? 1.0 : 0.0;
//     }
//     boundary({mark});

//     // 3. 扩散标记以覆盖半径 1.0
//     // 计算当前最细网格尺寸对应的步数
//     // 半径 R = 1.0
//     // 需要步数 N ≈ R / Delta_min
//     // L0 = 100, Level_Max = 11 -> Delta_min ≈ 100 / 2048 ≈ 0.0488
//     // N ≈ 1.0 / 0.0488 ≈ 20.5 -> 取 22 次
    
//     // 动态计算需要的扩散次数 (基于当前最大层级)
//     // 注意：这里用全局最大层级估算，保证覆盖范围足够大
//     double delta_min = L0 / (1 << LEVEL_MAX);
//     int n_expand = (int)(1.0 / delta_min) + 2; 
    
//     // 限制一下上限，避免极端情况卡死
//     if (n_expand > 50) n_expand = 50;

//     // 迭代扩散
//     // 使用临时场来避免迭代过程中的方向性偏差
//     for (int k = 0; k < n_expand; k++) {
//         scalar mark_new[];
//         foreach() {
//             mark_new[] = mark[]; // 继承当前状态
//         }
        
//         // 遍历更新：如果邻居被标记，自己也被标记
//         foreach() {
//             if (mark[] < 0.5) { // 只有未被标记的点才需要检查
//                 // 检查周围 8 个邻居 (包括对角线，这样扩散形状更接近圆/八边形)
//                 if (mark[1,0] > 0.5 || mark[-1,0] > 0.5 || 
//                     mark[0,1] > 0.5 || mark[0,-1] > 0.5 ||
//                     mark[1,1] > 0.5 || mark[1,-1] > 0.5 ||
//                     mark[-1,1] > 0.5 || mark[-1,-1] > 0.5) {
//                     mark_new[] = 1.0;
//                 }
//             }
//         }
        
//         //以此更新 mark 场
//         foreach() {
//             mark[] = mark_new[];
//         }
//         boundary({mark});
//     }

//     // 4. 执行加密
//     refine(mark[] > 0.5 && level < 11);
// }
// #endif

 event volume_stats (i++) {
    if (i == 0) return 0;
    update_volume_stats();
    if (i % LOG_EVERY == 0 && pid() == 0) {
        FILE *fp = fopen("volume.csv", "a");
        if (fp) {
            double rel = V0 > 0. ? fabs(V_total - V0)/V0 : NAN;
        fprintf(fp, "%d,%.12g,%.15g,%.15g,%.9g,%.9g,%d,%.8g,%d,%.8g,%.8g\n",
                    i, t, V_total, V_excess, rel, max_rel_volume_drift,
                    last_picard_iter, last_picard_res, last_mg_i, last_mg_resa, last_mg_resb);
            fclose(fp);
        }
        FILE *state = fopen("volume_state.txt", "w");
        if (state) {
            fprintf(state, "V0=%.17g\nmax_rel_volume_drift=%.17g\n", V0, max_rel_volume_drift);
            fclose(state);
        }
    }
  }

event droplet_radius(i += RADIUS_EVERY) {
    double wet_area = 0., volume = 0., xsum = 0., ysum = 0.;
    foreach(reduction(+:wet_area) reduction(+:volume) reduction(+:xsum) reduction(+:ysum)) {
        double g = max(h[] - HP, 0.);
        if (h[] > H_CONTACT) wet_area += sq(Delta);
        volume += g*sq(Delta);
        xsum += x*g*sq(Delta);
        ysum += y*g*sq(Delta);
    }
    R_contact = sqrt(wet_area/pi);
    double xc = volume > 0. ? xsum/volume : L0/2.;
    double yc = volume > 0. ? ysum/volume : L0/2.;
    double m2x = 0., m2y = 0.;
    foreach(reduction(+:m2x) reduction(+:m2y)) {
        double g = max(h[] - HP, 0.);
        m2x += sq(x - xc)*g*sq(Delta);
        m2y += sq(y - yc)*g*sq(Delta);
    }
    double ell = volume > 0. ? sqrt((m2x + m2y)/(2.*volume)) : 0.;
    if (pid() == 0) {
        static int header_written = 0;
        FILE *fp = fopen("radius.dat", header_written ? "a" : "w");
        if (fp) {
            if (!header_written) fprintf(fp, "# i t ell_second_moment R_area wet_area V_excess xc yc dt\n");
            fprintf(fp, "%d %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",
                    i, t, ell, R_contact, wet_area, V_excess, xc, yc, dt);
            fclose(fp);
            header_written = 1;
        }
    }
}

event contour_output(i++) {
    if (contour_target_next >= n_contour_targets)
        return 0;
        while (contour_target_next < n_contour_targets &&
           t >= contour_targets[contour_target_next]) {
        double target_t = contour_targets[contour_target_next];
        int tag = (int)(target_t + 0.5);
        char snapshot_tag[32];
        snprintf(snapshot_tag, sizeof(snapshot_tag), "t%03d", tag);
        write_snapshot(snapshot_tag, i, t);
        if (pid() == 0)
            fprintf(stderr, "CONTOUR target_t=%.12g actual_t=%.12g step=%d\n", target_t, t, i);
        contour_target_next++;
    }
    return 0;
}

/* --- 日志与输出事件 --- */
 
 event logWriting (i++) {
    if (i % LOG_EVERY) return 0.;
    if (pid() != 0) return 0.;
  
    const char *fname = "log_restart";
  
    // --- 确保表头只在“文件不存在或为空”时写一次 ---
    static int header_checked = 0;
    if (!header_checked) {
      struct stat st;
      int need_header = (stat(fname, &st) != 0) || (st.st_size == 0);
  
      FILE *fp0 = fopen(fname, need_header ? "w" : "a");
      if (fp0) {
        if (need_header) {
            fprintf(fp0, "i h1 h2 t dt lm V_total V_excess picard_iter picard_res mg_i mg_resa_after mg_resb_before\n");
        }
        fclose(fp0);
      }
      header_checked = 1;
    }
  
    // --- 每 100 步写时间戳（追加注释行）---
    if (i % 100 == 0) {
      time_t current_time = time(NULL);
  
      // 更推荐把 TZ 设置移到 main() 一次性做；这里先保留你的写法
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
  
    // --- 正常数据行追加 ---
    FILE *fp = fopen(fname, "a");
    if (fp) {
        fprintf(fp, "%d %g %g %.12g %g %d %.15g %.15g %d %.8g %d %.8g %.8g\n", i, hmin1, hmin, t, tcw_dt, lm, V_total, V_excess, last_picard_iter, last_picard_res, last_mg_i, last_mg_resa, last_mg_resb);
      fclose(fp);
    }
  }
  
 
 // Legacy safe_figure time event removed: it forces t=28 and t=10000.
//  event outputH_byStep(i += 100) {
//      char filename[256];
//      snprintf(filename, sizeof(filename), "%s/h_ada_%08d.bin", dirFlowFields, i);
     
//      // 使用合并进来的二进制保存函数
//      save_binary_data({h}, filename);
//  }
// event outputsolutionadapt(i++)
// {

// }
event outputsolutionadapt(i++ ) {
    if (i % contour_every == 0) {
      
      char name[80];
      sprintf (name, "dat/adasnapshot-%08d.dat", i);
      
      // 计算总单元数 - 所有线程参与
      int n_cells = 0;
      foreach(reduction(+:n_cells))
        n_cells++;
    
        FILE * fpx = fopen(name, "w");
        
        // 输出文件头
        fprintf(fpx, "VARIABLES= X, Y, h \n");
        
        // 输出单元和节点信息
        fprintf(fpx, "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n", 
                n_cells*4, n_cells);
        
        // 使用foreach_serial禁用并行，确保只有线程1执行
        foreach(serial) {
          // 输出当前进程编号
          // printf("输出节点: 进程 %d, 位置 (%g, %g)\n", omp_get_thread_num(), x, y);
          
          // 输出四个顶点
          fprintf(fpx, "%g\t%g\t%g\n", x - Delta/2., y - Delta/2., interpolate(h, x - Delta/2., y - Delta/2.));
          fprintf(fpx, "%g\t%g\t%g\n", x + Delta/2., y - Delta/2., interpolate(h, x + Delta/2., y - Delta/2.));
          fprintf(fpx, "%g\t%g\t%g\n", x + Delta/2., y + Delta/2., interpolate(h, x + Delta/2., y + Delta/2.));
          fprintf(fpx, "%g\t%g\t%g\n", x - Delta/2., y + Delta/2., interpolate(h, x - Delta/2., y + Delta/2.));
        }
        
        // 输出单元连接关系
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

event checkpoint_output(t = 5.; t += 5.) {
    char path[96], tag[32];
    int checkpoint_time = (int)(t + 0.5);
    update_volume_stats();
    update_hp001_indicator();
    snprintf(path, sizeof(path), "checkpoint_t%04d", checkpoint_time);
    dump(file = path);
    snprintf(tag, sizeof(tag), "t%04d", checkpoint_time);
    write_snapshot(tag, i, t);
    if (pid() == 0) {
        FILE *state = fopen("volume_state.txt", "w");
        if (state) {
            fprintf(state, "V0=%.17g\nmax_rel_volume_drift=%.17g\n", V0, max_rel_volume_drift);
            fclose(state);
        }
        fprintf(stderr, "CHECKPOINT step=%d time=%.12g file=%s\n", i, t, path);
    }
}

event tcw_stop(t = tcw_end_time) {
    update_volume_stats();
    update_hp001_indicator();
    dump(file = "native_final");
    write_snapshot("final", i, tcw_end_time);
    if (pid() == 0) {
        FILE *fp = fopen("volume.csv", "a");
        if (fp) {
            double rel = V0 > 0. ? fabs(V_total - V0)/V0 : NAN;
            fprintf(fp, "%d,%.12g,%.15g,%.15g,%.9g,%.9g,%d,%.8g,%d,%.8g,%.8g\n",
                    i, t, V_total, V_excess, rel, max_rel_volume_drift,
                    last_picard_iter, last_picard_res, last_mg_i, last_mg_resa, last_mg_resb);
            fclose(fp);
        }
        FILE *state = fopen("volume_state.txt", "w");
        if (state) {
            fprintf(state, "V0=%.17g\nmax_rel_volume_drift=%.17g\n", V0, max_rel_volume_drift);
            fclose(state);
        }
    }
    return 1;
}

// Event registration follows source order; keep integration after diagnostics and the stop event.
event integration(i++) {
     if (t >= tcw_end_time) return 0;
     if (dt > tcw_end_time - t) dt = tcw_end_time - t;
 #if T2OR
     // 二阶时间积分
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
     // 一阶时间积分
     dtold = dt;
     // dt is scheduled in stability() before this event
     foreach() h_noise_old[] = h[];
     boundary({h_noise_old});
#if !TCW_DIAGNOSTIC_NO_DETERMINISTIC
    // A uniform film has zero deterministic RHS; avoid the legacy implicit
    // solve on that first step, then apply the conservative stochastic flux.
    double h_lo = HUGE_VAL, h_hi = -HUGE_VAL;
    foreach(reduction(min:h_lo) reduction(max:h_hi)) {
        if (h[] < h_lo) h_lo = h[];
        if (h[] > h_hi) h_hi = h[];
    }
    last_picard_iter = 0;
    last_picard_res = 0.;
    last_mg_i = 0;
    last_mg_resa = 0.;
    last_mg_resb = 0.;
    if (h_hi - h_lo > 1e-14 || HP001_FORCE_COUPLED_SOLVER)
        solve_thin_film(h, dt, m_exp, n_exp, dd_exp, i);
#endif
    apply_split_stochastic_flux(h, h_noise_old, dt, i);
 #endif
 
     // 更新统计信息与限制器
     double minres = 1000;
     int levelmax = 1, rupture = 0;
     levelm = -1;
     hmin1 = hmin;
     
     foreach(reduction(min:minres) reduction(max:levelmax) reduction(max:rupture)) {
         // No hard floor: clipping injects mass.
         if (!isfinite(h[]) || h[] < RUPTURE_HEIGHT) rupture = 1;
         hi[] = 1.0/max(h[], 1e-12);
         double res = fabs(h[]);
         int ll = level;
         
         if (res < minres) minres = res;
         if (ll > levelmax) levelmax = ll;
     }
     hmin = minres;
     lm = levelmax;
     int solver_bad = !isfinite(last_picard_res) || last_picard_res > 1e-8 ||
                      !isfinite(last_mg_resa) || last_mg_resa > TOLERANCE;
     if (solver_bad) {
         write_snapshot("solver_failure", i, t + dt);
         if (pid() == 0) {
             FILE *failure = fopen("SOLVER_FAILURE.txt", "w");
             if (failure) {
                 fprintf(failure, "step=%d\nstate_time=%.17g\ndt=%.17g\npicard_iter=%d\npicard_res=%.17g\npicard_tolerance=1e-8\nmg_i=%d\nmg_resa_after=%.17g\nmg_resb_before=%.17g\nmg_tolerance=%.17g\n",
                         i, t + dt, dt, last_picard_iter, last_picard_res,
                         last_mg_i, last_mg_resa, last_mg_resb, TOLERANCE);
                 fclose(failure);
             }
             fprintf(stderr, "SOLVER_FAILURE step=%d state_time=%.12g dt=%.12g picard_iter=%d picard_res=%.8g picard_tol=1e-8 mg_i=%d mg_resa_after=%.8g mg_resb_before=%.8g mg_tol=%.8g\n",
                     i, t + dt, dt, last_picard_iter, last_picard_res,
                     last_mg_i, last_mg_resa, last_mg_resb, TOLERANCE);
         }
         exit(EXIT_FAILURE);
     }
     if (rupture) {
         dump(file = "native_invalid");
         if (pid() == 0) fprintf(stderr, "FATAL: rupture threshold step=%d t=%.12g hmin=%.16g; native_invalid saved\n", i, t, hmin);
         exit(EXIT_FAILURE);
     }
}

event hp001_invalid_guard (i++) {
    if (hp001_height_invalid()) {
        dump(file = "native_invalid");
        if (pid() == 0)
            fprintf(stderr, "FATAL: non-positive or nonfinite h; native_invalid saved\n");
        exit(EXIT_FAILURE);
    }
    return 0;
}

/* Required for deterministic fresh/restart event timing on this solver. */
event restart_step_clock (t = tcw_dt; t += tcw_dt) {}
