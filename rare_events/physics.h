/* Generated from film_rupture at fed7bb3bf963afe44466017e61199c61613b417d.
 * Only beta/gamma repulsion removal changes the physical force.
 * Original source remains unchanged; see provenance.json. */
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
             /* User model: no h^-5 repulsive contribution. */
             double ff = (f[] + f[-1]) / 2.;
 
            if (ff <= 0) ff = 1e-9;

            const double mob = 1.0 + tcw_slip / ff;   /* Navier slip: m = h^3 (1+ell/h) */

            D.x[] = pow(ff, p.n) * mob;

            beta.x[] = -p.n * COEFF_A * pow(ff, p.n - p.m - 1.) * mob;

            gamma.x[] = (p.n * pow(ff, p.n - 1.) * (pres[] - pres[-1]) / Delta -
                         p.n * COEFF_A * (p.n - p.m - 1.) * pow(ff, p.n - p.m - 2.) * (f[] - f[-1]) / Delta) * mob;
             
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
