/* Periodic bilinear interpolation of the original L6 staggered Gaussian field.
 * L0 and tcw_noise_seed are supplied by the unchanged thin-film solver.
 * The flux denominator MUST use L0/64, with no variance renormalization.
 */
#include <math.h>
#include <stdint.h>
#ifndef pi
/* Match Basilisk common.h exactly, including its original precision. */
#define pi 3.14159265358979
#endif

static uint64_t tcw_mix64 (uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30))*0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27))*0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static double tcw_u01 (uint64_t key) {
    return ((double)(tcw_mix64(key) >> 11)+.5)*(1./9007199254740992.);
}
static double tcw_native_face_gaussian_noise (double xx, double yy, double delta,
                                              unsigned component, int timestep) {
    int64_t ix=llround(xx/delta), iy=llround(yy/delta);
    int64_t nidx=(int64_t)llround(L0/delta);
    ix=(ix%nidx+nidx)%nidx; iy=(iy%nidx+nidx)%nidx;
    uint64_t key=tcw_noise_seed^((uint64_t)ix*0x9e3779b97f4a7c15ULL)^
        ((uint64_t)iy*0xbf58476d1ce4e5b9ULL)^
        ((uint64_t)component*0x94d049bb133111ebULL)^
        ((uint64_t)timestep*0xd6e8feb86659fd93ULL);
    double u1=tcw_u01(key), u2=tcw_u01(key^0xa0761d6478bd642fULL);
    return sqrt(-2.*log(u1))*cos(2.*pi*u2);
}
static double tcw_coarse_noise (int ix, int iy, unsigned component, int step) {
    /* Lazy per-step cache: deterministic on every MPI rank, no RNG state. */
    static int stamps[2][64][64];
    static double values[2][64][64];
    ix=(ix%64+64)%64; iy=(iy%64+64)%64;
    if (stamps[component][iy][ix] != step+1) {
        double dc=L0/64.;
        double ox=component == 0 ? 0. : .5;
        double oy=component == 0 ? .5 : 0.;
        values[component][iy][ix]=tcw_native_face_gaussian_noise(
            (ix+ox)*dc, (iy+oy)*dc, dc, component, step);
        stamps[component][iy][ix]=step+1;
    }
    return values[component][iy][ix];
}
static double tcw_face_gaussian_noise (double xx, double yy, double delta,
                                      unsigned component, int timestep) {
    double dc=L0/64.;
    /* Exact original call at L6, including its floating-point coordinates. */
    if (fabs(delta-dc) < 1e-12*dc)
        return tcw_native_face_gaussian_noise(xx, yy, delta, component, timestep);
    double ox=component == 0 ? 0. : .5;
    double oy=component == 0 ? .5 : 0.;
    double sx=xx/dc-ox, sy=yy/dc-oy;
    int ix=(int)floor(sx), iy=(int)floor(sy);
    double ax=sx-ix, ay=sy-iy;
    return (1.-ay)*((1.-ax)*tcw_coarse_noise(ix,iy,component,timestep)+
                   ax*tcw_coarse_noise(ix+1,iy,component,timestep))+
           ay*((1.-ax)*tcw_coarse_noise(ix,iy+1,component,timestep)+
               ax*tcw_coarse_noise(ix+1,iy+1,component,timestep));
}
