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
#include <assert.h>
 
 /* --- 全局参数与宏定义 (原 parameters.h & layer.c) --- */
 
 #define T2OR 0          // 0: 一阶时间积分, 1: 二阶时间积分
 
 // 物理系数
 #define COEFF_A 0.0
 #define COEFF_VARPHI (SLE_PHI*noise_factor)
#define SLIP_LENGTH_X 0.16382252559726962
#define SLIP_LENGTH_Y 0.16382252559726962
#define TCW_DIAGNOSTIC_NO_DETERMINISTIC 0
 

#define SLE_REQUIRE(test) do { if (!(test)) { fprintf(stderr,"SLE_REQUIRE failed: %s at %d\n",#test,__LINE__); exit(36); } } while (0)
#define SLE_PHI 0.0089933890196434232
#define TIME_UNIT_NS 0.16596907894736843
static double noise_factor=1., tcw_end=11.875122248143917;
static int init_mode=0;
static double max_solver_residual=0.,max_nonlinear_update=0.;

 // 目录路径
 char dirFlowFields[] = "dat/";
 char dirFigures[]    = "png/";
 
 // 全局物理参数

 double n_exp = 3.0;     // 指数 n
 double m_exp = 3.0;     // 指数 m
 double dd_exp = 5.0;    // 指数 dd (分离压)
 
 // 时间步长参数
 double DTMAX = 5e-5;
 double DTMIN = 5e-5;
double tcw_dt = 5e-5;
int tcw_test_steps = 3;
 double dtold;
 double errorstep = 0.;

 // --- 体积统计 ---
double V_total  = 0.0;   // ∫ h dA
double V_excess = 0.0;   // ∫ max(h - HOUT, 0) dA（可选）
double HOUT = 0.1;
static uint64_t tcw_noise_seed = 1;       // 若你有基底厚度（hout），填对应值；没有就保持 0
 // 网格与求解器状态变量
 scalar h[], hi[], l[], h_noise_old[];
 #if T2OR
 scalar h2[], h1[], err[];
 #endif
 
 double hmin, hmin1;
 int lm = 10;
 int levelm;
 int LEVEL_MAX = 6;
 int LEVEL_MIN = 6;
 
 /* --- 输出管理功能 (原 outputManager.h & myoutput.h) --- */
 
 /**
  * @brief 创建目录 (如果不存在)
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
 
 /* Anisotropic Navier slip:
 * Mx(h)=h^3+3*lx*h^2, My(h)=h^3+3*ly*h^2.
 * Keep the two directions explicit so deterministic, linearized, and
 * stochastic face fluxes cannot silently collapse to one isotropic M(h).
 */
static inline double mobility_x (double hh) {
    return hh*hh*(hh + 3.*SLIP_LENGTH_X);
}

static inline double mobility_y (double hh) {
    return hh*hh*(hh + 3.*SLIP_LENGTH_Y);
}

static inline double mobility_prime_x (double hh) {
    return 3.*hh*(hh + 2.*SLIP_LENGTH_X);
}

static inline double mobility_prime_y (double hh) {
    return 3.*hh*(hh + 2.*SLIP_LENGTH_Y);
}



static inline double disjoining_prime (double hh, const struct EquationParams *p,
                                       double hstar) {
    return -p->m*COEFF_A*pow(hh, -p->m - 1.) +
            p->dd*COEFF_A*hstar*pow(hh, -p->dd - 1.);
}

static inline double disjoining_second (double hh, const struct EquationParams *p,
                                        double hstar) {
    return p->m*(p->m + 1.)*COEFF_A*pow(hh, -p->m - 2.) -
           p->dd*(p->dd + 1.)*COEFF_A*hstar*pow(hh, -p->dd - 2.);
}

/**
 * @brief 计算残差
 */
static double residual_thin_film(scalar *al, scalar *bl, scalar *resl, void *data) {
     scalar h_sol = al[0], b_rhs = bl[0], res = resl[0];
     struct EquationParams *p = (struct EquationParams *)data;
     
     (const) face vector alpha = p->D;
     (const) face vector beta = p->beta;
     (const) face vector gamma = p->gamma;
     (const) scalar lambda = p->lambda;
     
     double maxres = 0.;
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
     // 替换 boundary_flux 为 boundary
     boundary((scalar *){g, hf});
 
     foreach(reduction(max:maxres)) {
         res[] = b_rhs[] - lambda[] * h_sol[];
         foreach_dimension()
             res[] -= (g.x[1] - g.x[]) / Delta + 
                      (hf.x[1] - hf.x[]) / Delta ;
         
         if (fabs(res[]) > maxres)
             maxres = fabs(res[]);
     }
     boundary(resl);
     return maxres;
 }
 
 /**
  * @brief 松弛/平滑函数
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
  * @brief 求解器主函数
  */
 trace
 mgstats solve_thin_film(scalar f_scalar, double dt_step, double m_val, double n_val, double dd_val, int timestep) {
     struct EquationParams p;
     p.dt = dt_step;
     p.f = f_scalar;
     p.m = m_val;
     p.n = n_val;
     p.dd = dd_val;
     p.tolerance = TOLERANCE;
 
     if (p.dt == 0.) {
         mgstats s = {0};
         return s;
     }
 
     scalar f = p.f;
     scalar b[], fold[], fa[];
     const scalar lambda[] = -1. / p.dt;
     
     face vector D[], beta[], gamma[];
     
     foreach() {
         fold[] = f[];
     }
 
     double maxres;
     mgstats s;
     int niter = 0, nitermax = 10;
     double default_tol = TOLERANCE;
 
     // 初始化噪声场
 
     do {
         foreach_face()
             beta.x[] = face_gradient_x(f, 0);
         
         boundary((scalar *){beta});
 
         scalar pres[];
         foreach() {
             pres[] = 0.;
             foreach_dimension()
                 pres[] -= (beta.x[1] - beta.x[0]) / Delta;
         }
         boundary({pres});
 
         face vector g[];
 
         foreach_face(x) {
             double hstar = pow(0.1, p.dd - 3.);
             double ff = (f[] + f[-1]) / 2.;
             SLE_REQUIRE(isfinite(ff) && ff > 0.);
             double mob = mobility_x(ff);
             double mob_prime = mobility_prime_x(ff);
             double pi_prime = disjoining_prime(ff, &p, hstar);
             double beta_prime = mob_prime*pi_prime +
                                 mob*disjoining_second(ff, &p, hstar);

             D.x[] = mob;
             beta.x[] = mob*pi_prime;
             gamma.x[] = mob_prime*(pres[] - pres[-1])/Delta +
                         beta_prime*(f[] - f[-1])/Delta;
             g.x[] = gamma.x[] * ff;
         }
         foreach_face(y) {
             double hstar = pow(0.1, p.dd - 3.);
             double ff = (f[] + f[0,-1]) / 2.;
             SLE_REQUIRE(isfinite(ff) && ff > 0.);
             double mob = mobility_y(ff);
             double mob_prime = mobility_prime_y(ff);
             double pi_prime = disjoining_prime(ff, &p, hstar);
             double beta_prime = mob_prime*pi_prime +
                                 mob*disjoining_second(ff, &p, hstar);

             D.y[] = mob;
             beta.y[] = mob*pi_prime;
             gamma.y[] = mob_prime*(pres[] - pres[0,-1])/Delta +
                         beta_prime*(f[] - f[0,-1])/Delta;
             g.y[] = gamma.y[] * ff;
         }
         boundary((scalar *){g});
 
         foreach() {
             b[] = -fold[] / p.dt;
             fa[] = f[];
             foreach_dimension()
                 b[] += (g.x[1] - g.x[]) / Delta;
         }
 
         restriction({beta, D, lambda, gamma});
         
         p.lambda = lambda;
         p.D = D;
         p.beta = beta;
         p.gamma = gamma;
         
         s = mg_solve({f}, {b}, residual_thin_film, relax_thin_film, &p);
 
         maxres = 0;
         foreach(reduction(max:maxres)) {
             double res = fabs(f[] - fa[]);
             if (res > maxres)
                 maxres = res;
         }
         boundary({fa});
         niter++;
 
     } while (maxres > 1e-8 && niter < nitermax);
 
     if (p.tolerance > 1e-10)
         TOLERANCE = default_tol;
 

     if (!isfinite(maxres) || maxres > 1e-8 || !isfinite(s.resa) || s.resa > 1e-6) {
         fprintf(stderr,"SOLVER_FAIL t=%.17g dt=%.17g update=%.17g mgres=%.17g iterations=%d\n",t,dt_step,maxres,s.resa,niter);
         exit(34);
     }
     if (s.resa > max_solver_residual) max_solver_residual=s.resa;
     if (maxres > max_nonlinear_update) max_nonlinear_update=maxres;
     return s;
 }

 /**
  * @brief 自适应时间步长调整 (原 lininterp)
  */
 double adjust_timestep(double h1, double h2, double dt_curr, int lm_curr, double dtmax_curr, int step_idx) {
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
 

 static uint64_t tcw_mix64 (uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30))*0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27))*0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}


static double tcw_uniform_open01 (uint64_t key) {
    return ((double)(tcw_mix64(key) >> 11) + .5)*(1./9007199254740992.);
}

static double tcw_face_gaussian_noise (double xx, double yy, double delta,
                                       unsigned component, int timestep) {
    int64_t ix = llround(xx/delta);
    int64_t iy = llround(yy/delta);
    int64_t n=(int64_t)llround(L0/delta);
    ix=(ix%n+n)%n; iy=(iy%n+n)%n;
    uint64_t key = tcw_noise_seed ^ ((uint64_t)ix*0x9e3779b97f4a7c15ULL) ^
                   ((uint64_t)iy*0xbf58476d1ce4e5b9ULL) ^
                   ((uint64_t)component*0x94d049bb133111ebULL) ^
                   ((uint64_t)timestep*0xd6e8feb86659fd93ULL);
    double u1 = tcw_uniform_open01(key);
    double u2 = tcw_uniform_open01(key ^ 0xa0761d6478bd642fULL);
    return sqrt(-2.*log(u1))*cos(2.*acos(-1.)*u2);
}

void apply_split_stochastic_flux (scalar f, scalar f_old, double dt_step,
                                  int timestep) {
    face vector flux[];
    foreach_face(x) {
        double ff = (f_old[] + f_old[-1])/2.;
        SLE_REQUIRE(isfinite(ff) && ff > 0.);
        double eta_x = tcw_face_gaussian_noise(x, y, Delta, 0, timestep);
        flux.x[] = sqrt(2.*COEFF_VARPHI*mobility_x(ff))*eta_x/
                   (Delta*sqrt(dt_step));
    }
    foreach_face(y) {
        double ff = (f_old[] + f_old[0,-1])/2.;
        SLE_REQUIRE(isfinite(ff) && ff > 0.);
        double eta_y = tcw_face_gaussian_noise(x, y, Delta, 1, timestep);
        flux.y[] = sqrt(2.*COEFF_VARPHI*mobility_y(ff))*eta_y/
                   (Delta*sqrt(dt_step));
    }
    boundary((scalar *){flux});
    foreach() {
        foreach_dimension()
            f[] += dt_step*(flux.x[1] - flux.x[])/Delta;
    }
    boundary({f});
}


static void field_stats(double *mean,double *lo,double *hi) {
    double sum=0.,minimum=HUGE_VAL,maximum=-HUGE_VAL;
    foreach(reduction(+:sum) reduction(min:minimum) reduction(max:maximum)) {
        SLE_REQUIRE(isfinite(h[]) && h[]>0.);
        sum+=h[]*sq(Delta); minimum=fmin(minimum,h[]); maximum=fmax(maximum,h[]);
    }
    *mean=sum/sq(L0); *lo=minimum; *hi=maximum;
    if (fabs(*mean-1.)>1e-8) {
        fprintf(stderr,"MASS_FAIL t=%.17g mean=%.17g\n",t,*mean);exit(35);
    }
}

int main(int argc,char **argv) {
    if (argc>1) tcw_noise_seed=strtoull(argv[1],NULL,10);
    if (argc>2) tcw_dt=strtod(argv[2],NULL);
    if (argc>3) tcw_end=strtod(argv[3],NULL);
    if (argc>4) LEVEL_MIN=LEVEL_MAX=atoi(argv[4]);
    if (argc>5) noise_factor=strtod(argv[5],NULL);
    if (argc>6) init_mode=atoi(argv[6]);
    SLE_REQUIRE(tcw_dt>0. && tcw_end>0. && LEVEL_MIN>=5 && LEVEL_MIN<=7);
    SLE_REQUIRE((noise_factor==0. || noise_factor==1.) && init_mode>=0 && init_mode<=2);
    SLE_REQUIRE(fabs(mobility_x(1.)-(1.+3.*SLIP_LENGTH_X))<1e-14);
    SLE_REQUIRE(fabs(mobility_y(1.)-(1.+3.*SLIP_LENGTH_Y))<1e-14);
    SLE_REQUIRE(fabs(mobility_prime_x(1.)-(3.+6.*SLIP_LENGTH_X))<1e-14);
    SLE_REQUIRE(fabs(mobility_prime_y(1.)-(3.+6.*SLIP_LENGTH_Y))<1e-14);
    TOLERANCE=1e-7;
    periodic(right);periodic(top);
    size(10.716382252559727);origin(0.,0.);init_grid(1<<LEVEL_MIN);
    fprintf(stdout,"CONF seed=%llu dt=%.17g tend=%.17g level=%d slip_over_H=%.17g phi=%.17g binding=0 time_unit_ns=%.17g\n",(unsigned long long)tcw_noise_seed,tcw_dt,tcw_end,LEVEL_MIN,SLIP_LENGTH_X,COEFF_VARPHI,TIME_UNIT_NS);
    run();
}

event init(i=0) {
    foreach() {
        double phase=init_mode==1 ? x : y;
        h[]=1.+(init_mode ? .005*cos(4.*pi*phase/L0) : 0.);
    }
    boundary({h});
}

event stability(i++) {dt=dtnext(tcw_dt);}

/* Sampling precedes integration: labels refer to the current physical state. */
event sle_samples(t={0,0.25815483148138951,0.51630966296277903,0.77446449444416854,1.0326193259255581,1.2907741574069476,1.5489289888883371,1.8070838203697264,2.0652386518511161,2.3233934833325054,2.5815483148138951,2.8397031462952844,3.0978579777766742,3.3560128092580634,3.6141676407394527,3.8723224722208425,4.1304773037022322,4.3886321351836211,4.6467869666650108,4.9049417981464005,5.1630966296277903,5.4212514611091791,5.6794062925905688,5.9375611240719586,6.1957159555533483,6.4538707870347372,6.7120256185161269,6.9701804499975157,7.2283352814789055,7.4864901129602961,7.7446449444416849,8.0027997759230747,8.2609546074044644,8.5191094388858541,8.7772642703672421,9.0354191018486318,9.2935739333300216,9.5517287648114113,9.809883596292801,10.068038427774191,10.326193259255581,10.584348090736968,10.842502922218358,11.100657753699748,11.358812585181138,11.616967416662527,11.875122248143917}) {
    static int index=0;
    const int mdstep[]={0,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,55000,60000,65000,70000,75000,80000,85000,90000,95000,100000,105000,110000,115000,120000,125000,130000,135000,140000,145000,150000,155000,160000,165000,170000,175000,180000,185000,190000,195000,200000,205000,210000,215000,220000,225000,230000};
    char name[80];snprintf(name,sizeof(name),"sle_mdstep_%d.hfield",mdstep[index]);
    double mean,lo,hi;field_stats(&mean,&lo,&hi);
    write_hfield_binary(name);
    FILE *fp=fopen("samples.csv",index ? "a":"w");SLE_REQUIRE(fp);
    if (!index) fprintf(fp,"index,mdstep,t_dimensionless,time_ns,mean,hmin,hmax\n");
    fprintf(fp,"%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",index,mdstep[index],t,t*TIME_UNIT_NS,mean,lo,hi);fclose(fp);
    fprintf(stdout,"SLE_SAMPLE mdstep=%d t=%.17g time_ns=%.17g mean=%.17g\n",mdstep[index],t,t*TIME_UNIT_NS,mean);fflush(stdout);
    index++;
}

event first_kick(i=1) {
    write_hfield_binary("first_kick.hfield");
    FILE *fp=fopen("first_kick_time.txt","w");SLE_REQUIRE(fp);fprintf(fp,"%.17g\n",t);fclose(fp);
}

event sle_stop(t=tcw_end) {
    double mean,lo,hi;field_stats(&mean,&lo,&hi);
    write_hfield_binary("final.hfield");
    FILE *fp=fopen("DONE.json","w");SLE_REQUIRE(fp);
    fprintf(fp,"{\"steps\":%d,\"t_dimensionless\":%.17g,\"time_ns\":%.17g,\"mean\":%.17g,\"hmin\":%.17g,\"hmax\":%.17g,\"max_solver_residual\":%.17g,\"max_nonlinear_update\":%.17g}\n",i,t,t*TIME_UNIT_NS,mean,lo,hi,max_solver_residual,max_nonlinear_update);fclose(fp);
    fprintf(stdout,"SLE_COMPLETE steps=%d t=%.17g time_ns=%.17g mean=%.17g\n",i,t,t*TIME_UNIT_NS,mean);return 1;
}

event integration(i++) {
    foreach() h_noise_old[]=h[];
    boundary({h_noise_old});
    double lo=HUGE_VAL,hi=-HUGE_VAL;
    foreach(reduction(min:lo) reduction(max:hi)) {lo=fmin(lo,h[]);hi=fmax(hi,h[]);}
    if (hi-lo>1e-14) solve_thin_film(h,dt,m_exp,n_exp,dd_exp,i);
    apply_split_stochastic_flux(h,h_noise_old,dt,i);
}

event progress(i+=100) {
    double mean,lo,hi;field_stats(&mean,&lo,&hi);
    FILE *fp=fopen("progress.csv",i ? "a":"w");SLE_REQUIRE(fp);
    if (!i) fprintf(fp,"i,t_state,time_ns,mean,hmin,hmax\n");
    fprintf(fp,"%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",i,t+dt,(t+dt)*TIME_UNIT_NS,mean,lo,hi);fclose(fp);
}
