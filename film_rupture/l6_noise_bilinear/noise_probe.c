#include <assert.h>
#include <stdio.h>
#include <stdint.h>
double L0;
static uint64_t tcw_noise_seed;
#include "noise_l6_bilinear.h"

int main(void) {
    for (int domain=0; domain<2; domain++) {
        L0=domain ? 100./.6 : 200.;
        double dc=L0/64.;
        /* Unique step per domain keeps the per-step cache contract explicit. */
        int step=domain+100;
        tcw_noise_seed=2026112101;
        for (unsigned c=0; c<2; c++) {
            double ox=c==0 ? 0. : .5, oy=c==0 ? .5 : 0.;
            for (int y=0; y<64; y++) for (int x=0; x<64; x++) {
                double xx=(x+ox)*dc, yy=(y+oy)*dc;
                double native=tcw_native_face_gaussian_noise(xx,yy,dc,c,step);
                assert(tcw_face_gaussian_noise(xx,yy,dc,c,step)==native);
                assert(tcw_coarse_noise(x,y,c,step)==native);
                double node=tcw_face_gaussian_noise(xx,yy,dc/4.,c,step);
                assert(fabs(node-native)<1e-12);
            }
            double a=tcw_coarse_noise(63,63,c,step);
            double b=tcw_coarse_noise(0,63,c,step);
            double d=tcw_coarse_noise(63,0,c,step);
            double e=tcw_coarse_noise(0,0,c,step);
            double xx=(63.25+ox)*dc, yy=(63.75+oy)*dc;
            double expected=.1875*a+.0625*b+.5625*d+.1875*e;
            double actual=tcw_face_gaussian_noise(xx,yy,dc/4.,c,step);
            assert(fabs(actual-expected)<1e-12);
            assert(fabs(actual-tcw_face_gaussian_noise(xx-L0,yy+L0,dc/4.,c,step))<1e-12);
        }
    }
    puts("PASS: L6 exact nodes, staggered bilinear weights, periodic seams, both domains");
    return 0;
}
