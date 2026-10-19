/**
 * @file    thinfilm_stochastic.c
 * @brief   Stochastic Thin-Film Equation (STFE) solver .
 * @date    2026-01-04
 * @author  <Liang jiahua/Department of Modern Mechanics,USTC><liangjh@mail.ustc.edu.cn>
 * @version 1.0.0
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
 
 /* --- 全局参数与宏定义 (原 parameters.h & layer.c) --- */
 
 #define T2OR 0          // 0: 一阶时间积分, 1: 二阶时间积分
 
 // 物理系数
double tcw_A = 0.0010875206611570248; double tcw_slip = 0.0;
double tcw_varphi = 0.0082644628099173554; double tcw_h0 = 0.5;
#define COEFF_A (tcw_A)
#define COEFF_VARPHI (tcw_varphi)
#define TCW_DIAGNOSTIC_NO_DETERMINISTIC 0
 
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
double tcw_dt = 1e-3;
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
 int lm = 7;
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
 
         foreach_face() {
             double hstar = 0.0; /* h^-5 repulsion disabled; h^-3 attraction retained. */
             double ff = (f[] + f[-1]) / 2.;
 
            if (ff <= 0) ff = 1e-9;

            const double mob = 1.0 + tcw_slip / ff;   /* Navier slip: m = h^3 (1+ell/h) */

            D.x[] = pow(ff, p.n) * mob;

            beta.x[] = (-p.n * COEFF_A * pow(ff, p.n - p.m - 1.) +
                        p.dd * COEFF_A * hstar * pow(ff, p.n - p.dd - 1.)) * mob;

            gamma.x[] = (p.n * pow(ff, p.n - 1.) * (pres[] - pres[-1]) / Delta -
                         p.n * COEFF_A * (p.n - p.m - 1.) * pow(ff, p.n - p.m - 2.) * (f[] - f[-1]) / Delta +
                         p.dd * COEFF_A * hstar * (p.n - p.dd - 1.) * pow(ff, p.n - p.dd - 2.) * (f[] - f[-1]) / Delta) * mob;
             
             g.x[] = gamma.x[] * ff;
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
 
 int main(int argc, char *argv[]) {
     unsigned int seed = argc > 1 ? (unsigned int)strtoul(argv[1], NULL, 10) : 1u;
     tcw_noise_seed = seed;
    if (argc > 2) tcw_dt = strtod(argv[2], NULL);
    if (argc > 3) tcw_test_steps = (int)strtol(argv[3], NULL, 10);
    if (!(tcw_dt > 0.) || tcw_test_steps < 2) return 2;
    if (argc > 4) tcw_slip = strtod(argv[4], NULL);
    if (argc > 5) tcw_h0   = strtod(argv[5], NULL);
    if (pid() == 0)
        fprintf(stdout, "CONF slip=%g h0=%g A=%g varphi=%g LEVEL_MIN=%d dt=%g steps=%d seed_arg=%u\n",
                tcw_slip, tcw_h0, tcw_A, tcw_varphi, LEVEL_MIN, tcw_dt, tcw_test_steps, seed);
     periodic(right);
     periodic(top);
     origin(0., 0.);
     init_grid(1 << LEVEL_MIN);
     L0 = 90.909090909090909;
     run();
 }
 
 event init(i = 0) {
     create_directory(dirFlowFields);
     create_directory(dirFigures);
     
     // 尝试从 checkpoint 恢复
     if (!restore(file = "Layer-00029100")) {
         if (pid() == 0) {
             printf("Initializing simulation...\n");
         }
          hmin = tcw_h0;
          hmin1 = tcw_h0 + 1e-9;

          foreach() {
              h[] = tcw_h0;
          }
          boundary({h});
     }
 }
 
 static uint64_t tcw_mix64 (uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30))*0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27))*0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static double tcw_u01 (uint64_t key) { return ((double)(tcw_mix64(key) >> 11)+.5)*(1./9007199254740992.); }
static double tcw_face_gaussian_noise (double xx,double yy,double delta,unsigned component,int timestep) {
 int64_t ix=llround(xx/delta),iy=llround(yy/delta);
 int64_t n = 1 << LEVEL_MIN;
 ix = (ix % n + n) % n;
 iy = (iy % n + n) % n;
 uint64_t key=tcw_noise_seed^((uint64_t)ix*0x9e3779b97f4a7c15ULL)^((uint64_t)iy*0xbf58476d1ce4e5b9ULL)^((uint64_t)component*0x94d049bb133111ebULL)^((uint64_t)timestep*0xd6e8feb86659fd93ULL); double u1=tcw_u01(key),u2=tcw_u01(key^0xa0761d6478bd642fULL); return sqrt(-2.*log(u1))*cos(2.*pi*u2);
}

void apply_split_stochastic_flux (scalar f, scalar f_old, double dt_step,
                                  int timestep) {
    face vector flux[];
    foreach_face(x) {
        double ff = (f_old[] + f_old[-1])/2.;
        if (ff <= 0.) ff = 1e-9;
        flux.x[] = sqrt(2.*COEFF_VARPHI)*pow(ff, 1.5)*sqrt(1.0 + tcw_slip/ff)*
                   tcw_face_gaussian_noise(x, y, Delta, 0, timestep)/
                   (Delta*sqrt(dt_step));
    }
    foreach_face(y) {
        double ff = (f_old[] + f_old[0,-1])/2.;
        if (ff <= 0.) ff = 1e-9;
        flux.y[] = sqrt(2.*COEFF_VARPHI)*pow(ff, 1.5)*sqrt(1.0 + tcw_slip/ff)*
                   tcw_face_gaussian_noise(x, y, Delta, 1, timestep)/
                   (Delta*sqrt(dt_step));
    }
    boundary((scalar *){flux});
    foreach() {
        foreach_dimension()
            f[] += dt_step*(flux.x[1] - flux.x[])/Delta;
    }
    boundary({f});
}
event stability(i++) {
double cap = 3e-3;   /* dt 封顶 3e-3: 抑制长亚稳期质量漂移 */
    if (i > 0 && i % 20 == 0 && tcw_dt < cap) tcw_dt *= 1.25;
    if (tcw_dt > cap) tcw_dt = cap;
    if (tcw_dt < 1e-3) tcw_dt = 1e-3;
    dt = dtnext(tcw_dt);
}

event integration(i++) {
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
    if (h_hi - h_lo > 1e-14)
        solve_thin_film(h, dt, m_exp, n_exp, dd_exp, i);
#endif
    apply_split_stochastic_flux(h, h_noise_old, dt, i);
 #endif
 
     // 更新统计信息与限制器
     double minres = 1000;
     int levelmax = 1;
     levelm = -1;
     hmin1 = hmin;
     
     foreach(reduction(min:minres) reduction(max:levelmax)) {
         // 限制最小厚度
         if (h[] < 0.005) h[] = 0.005;
         
         hi[] = 1.0 / h[];
         double res = fabs(h[]);
         int ll = level;
         
         if (res < minres) minres = res;
         if (ll > levelmax) levelmax = ll;
     }
     hmin = minres;
     lm = levelmax;
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
    double vtot = 0.0;
    double vex  = 0.0;
  
    foreach (reduction(+:vtot) reduction(+:vex)) {
      double dA = sq(Delta);   // 2D 单元面积
      double hh = h[];
      vtot += hh * dA;
  
      double ex = hh - HOUT;
      if (ex > 0.) vex += ex * dA;
    }
  
  /* foreach(reduction) already returns the global integrals. */
  
    V_total  = vtot;
    V_excess = vex;
  }
  
 /* --- 日志与输出事件 --- */
 
event logWriting (i += 1000) {
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
            fprintf(fp0, "i h1 h2 t dt lm V_total V_excess\n");
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
        fprintf(fp, "%d %g %g %g %g %d %.15g %.15g\n", i, hmin1, hmin, t, dt, lm, V_total, V_excess);
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
    if (i%10<0) {
      
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

event tcw_samples(i++) {
    int index = -1;
    if (i == 50) index = 0;
    else if (i == 500) index = 1;
    else if (i == 5000) index = 2;
    if (index >= 0) {
        char name[64];
        snprintf(name, sizeof(name), "tcw_t%04d.hfield", index);
        write_hfield_binary(name);
        if (pid() == 0) fprintf(stderr, "TCW_SAMPLE index=%d step=%d t=%.12g file=%s\\n", index, i, t, name);
    }
    return 0;
}

double tcw_nextout = 0.05;   /* 液膜输出时刻表: t 几何增长 ×1.5 */
event film_output(i++) {
    if (t >= tcw_nextout) {
        char nm[96];
        snprintf(nm, sizeof nm, "film_t%013.6f.hfield", t);
        write_hfield_binary(nm);
        if (pid() == 0) {
            FILE *fl = fopen("film_times.txt", "a");
            if (fl) { fprintf(fl, "%d %.6g %s\n", i, t, nm); fclose(fl); }
        }
        tcw_nextout = t * 1.5 + 0.05;
    }
}

event hfield_periodic (t = 0.; t += 100.) {
    char name[128];
    snprintf(name, sizeof(name), "hfield_t%09.1f.bin", t);
    write_hfield_binary(name);
    if (pid() == 0) fprintf(stderr, "HFIELD t=%.12g file=%s\n", t, name);
}

event rupture_stop(i++) { if (hmin <= 0.01) {
    if (pid() == 0) fprintf(stderr,"RUPTURE hmin=%.12g t=%.12g i=%d\n",hmin,t,i);
    /* 破碎瞬间快照: 只在停车步写, 不消耗随机数, 不改变轨迹; MPI 按进程分文件 */
    write_hfield_binary("rupture.hfield");
    char nmf[64], nms[64];
    snprintf(nmf, sizeof nmf, "rupture_field_r%d.dat", pid());
    snprintf(nms, sizeof nms, "rupture_section_r%d.dat", pid());
    FILE *f2 = fopen(nmf, "w");
    foreach()
        fprintf(f2, "%.6g %.6g %.6g\n", x, y, h[]);
    fclose(f2);
    FILE *fs = fopen(nms, "w");
    double ymid = L0 / 2.;
    foreach()
        if (fabs(y - ymid) < Delta)
            fprintf(fs, "%.6g %.6g\n", x, h[]);
    fclose(fs);
    return 1;
  }
  return 0;
}

event tcw_stop(i = tcw_test_steps) {
    write_hfield_binary("stop.hfield");
    return 1;
}
