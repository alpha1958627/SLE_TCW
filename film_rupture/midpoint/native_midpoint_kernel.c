/* Fixed square FV operators only. No timestep, RNG, floor or volume repair. */
#include <math.h>
#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __attribute__((visibility("default")))
#endif

EXPORT void film_prepare(int n, double dx, double a, double slip,
                         const double *h, double *b) {
    int size=n*n;
    double *mu=b, *mx=b+size, *my=b+2*size, *dmx=b+3*size, *dmy=b+4*size;
    double *sx=b+5*size, *sy=b+6*size, *h4=b+7*size;
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
        int p=i*n+j, xp=(i+1<n?i+1:0)*n+j, xm=(i?i-1:n-1)*n+j;
        int yp=i*n+(j+1<n?j+1:0), ym=i*n+(j?j-1:n-1);
        double lap=((h[xp]+h[xm]-2*h[p])+(h[yp]+h[ym]-2*h[p]))/(dx*dx);
        double fx=.5*(h[p]+h[xp]), fy=.5*(h[p]+h[yp]);
        mu[p]=-lap+a/pow(h[p],3.0); h4[p]=pow(h[p],4.0);
        mx[p]=pow(fx,3.0)+slip*pow(fx,2.0);
        my[p]=pow(fy,3.0)+slip*pow(fy,2.0);
        dmx[p]=3*pow(fx,2.0)+2*slip*fx; dmy[p]=3*pow(fy,2.0)+2*slip*fy;
        sx[p]=sqrt(mx[p]); sy[p]=sqrt(my[p]);
    }
}

EXPORT void film_force(int n, double dx, const double *b, double *out) {
    int size=n*n;
    const double *mu=b, *mx=b+size, *my=b+2*size;
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
        int p=i*n+j, xp=(i+1<n?i+1:0)*n+j, xm=(i?i-1:n-1)*n+j;
        int yp=i*n+(j+1<n?j+1:0), ym=i*n+(j?j-1:n-1);
        double fx=mx[p]*(mu[xp]-mu[p])/dx, bx=mx[xm]*(mu[p]-mu[xm])/dx;
        double fy=my[p]*(mu[yp]-mu[p])/dx, by=my[ym]*(mu[p]-mu[ym])/dx;
        out[p]=(fx-bx)/dx+(fy-by)/dx;
    }
}

EXPORT void film_noise(int n, double dx, double phi, const double *b,
                       const double *dw, double *out) {
    int size=n*n;
    const double *sx=b+5*size, *sy=b+6*size;
    double factor=sqrt(2*phi)/dx;
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
        int p=i*n+j, xm=(i?i-1:n-1)*n+j, ym=i*n+(j?j-1:n-1);
        out[p]=factor*((sx[p]*dw[p]-sx[xm]*dw[xm])/dx+
                      (sy[p]*dw[size+p]-sy[ym]*dw[size+ym])/dx);
    }
}

EXPORT void film_derivative(int n, double dx, double a, double phi,
                            const double *b, const double *v, const double *dw,
                            double *work, double *df, double *dn) {
    int size=n*n;
    const double *mu=b, *mx=b+size, *my=b+2*size, *dmx=b+3*size, *dmy=b+4*size;
    const double *sx=b+5*size, *sy=b+6*size, *h4=b+7*size;
    double *dmu=work, *fx=work+size, *fy=work+2*size, *nx=work+3*size, *ny=work+4*size;
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
        int p=i*n+j, xp=(i+1<n?i+1:0)*n+j, xm=(i?i-1:n-1)*n+j;
        int yp=i*n+(j+1<n?j+1:0), ym=i*n+(j?j-1:n-1);
        double lap=((v[xp]+v[xm]-2*v[p])+(v[yp]+v[ym]-2*v[p]))/(dx*dx);
        dmu[p]=-lap-3*a*v[p]/h4[p];
    }
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
        int p=i*n+j, xp=(i+1<n?i+1:0)*n+j, yp=i*n+(j+1<n?j+1:0);
        double vx=dmx[p]*(.5*(v[p]+v[xp])), vy=dmy[p]*(.5*(v[p]+v[yp]));
        fx[p]=vx*(mu[xp]-mu[p])/dx+mx[p]*(dmu[xp]-dmu[p])/dx;
        fy[p]=vy*(mu[yp]-mu[p])/dx+my[p]*(dmu[yp]-dmu[p])/dx;
        nx[p]=.5*vx/sx[p]*dw[p]; ny[p]=.5*vy/sy[p]*dw[size+p];
    }
    double factor=sqrt(2*phi)/dx;
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
        int p=i*n+j, xm=(i?i-1:n-1)*n+j, ym=i*n+(j?j-1:n-1);
        df[p]=(fx[p]-fx[xm])/dx+(fy[p]-fy[ym])/dx;
        dn[p]=factor*((nx[p]-nx[xm])/dx+(ny[p]-ny[ym])/dx);
    }
}
