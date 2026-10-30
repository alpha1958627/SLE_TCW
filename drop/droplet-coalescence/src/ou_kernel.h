/* Exact Gaussian SE cell-average kernel and exact OU endpoint/bin law.
   Extracted verbatim from the tested baseline; root-only use in AMR wrapper. */
#ifndef AMR_OU_KERNEL_H
#define AMR_OU_KERNEL_H
#include <stdint.h>
#include <complex.h>
#include <string.h>
#define TCW_OUSE_BIN_AVERAGES 1
static int LEVEL_MIN=7,tcw_ouse_level=7,tcw_noise_mode=2;
static int tcw_ouse_legacy_scale=0,tcw_noise_initialized=0,tcw_noise_n=0;
static double tcw_sigma=.05,tcw_tau=.01;
static uint64_t tcw_noise_seed=1;
static double *tcw_se_field[2]={NULL,NULL},*tcw_se_endpoint[2]={NULL,NULL};
static double *tcw_se_first=NULL,*tcw_se_innovation=NULL,*tcw_se_filter=NULL;
static double complex *tcw_se_fft=NULL,*tcw_se_col=NULL;
static double tcw_se_norm=0.;
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

static double tcw_grid_delta (void) {
    return L0/(1 << LEVEL_MIN);
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
#endif
