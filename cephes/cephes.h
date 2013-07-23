#ifndef CEPHES_H
#define CEPHES_H

extern double MACHEP;
extern double MINLOG;
extern double MAXLOG;
extern double MAXNUM;

extern double PI, PIO4;
extern double SQRTH;
 
double cephes_exp (double);              /* unity.c */
double expx2 (double, int);              /* expx2.c */
double cephes_gamma (double);            /* gamma.c */
double igam (double, double);            /* igam.c */
double igamc (double, double);           /* igam.c */
double igami (double, double);           /* igami.c */
double incbet (double, double, double);  /* incbet.c */
double incbi (double, double, double);   /* incbi.c */
double lgam (double);                    /* gamma.c */
double cephes_log (double);              /* unity.c */
double cephes_exp (double);              /* unity.c */
double ndtri (double);                   /* ndtri.c */
double p1evl (double, double *, int);    /* polevl.c */
double polevl (double, double *, int);   /* polevl.c */
double chbevl (double, double *, int);   /* chbevl.c */
double cbrt (double);                    /* cbrt.c: cube root */
double psi (double);                     /* psi.c: digamma fn */

int airy (double x, double *ai, double *aip, double *bi, double *bip);
double hyp2f0 (double a, double b, double c, int t, double *err);
double hyp2f1 (double a, double b, double c, double x);
double hyperg (double a, double b, double x);
double cephes_hankel (double n, double x);
double cephes_bessel_Jn (int n, double x);
double cephes_bessel_Yn (int n, double x);
double cephes_bessel_Jv (double n, double x);
double cephes_bessel_Yv (double n, double x);
double cephes_bessel_Iv (double v, double x);
double cephes_bessel_Kn (int nn, double x);
double cephes_bessel_I0 (double x);
double cephes_bessel_I1 (double x);
double cephes_bessel_K0 (double x);
double cephes_bessel_K1 (double x);

double cephes_j0 (double x);
double cephes_j1 (double x);
double cephes_y0 (double x);
double cephes_y1 (double x);

/* interloper from elsewhere in netlib repository */
double netlib_bessel_K (double v, double x, int ize);

#ifndef INFINITY
extern double INFINITY;
#endif

#ifndef NAN
extern double NAN;
#endif

#endif /* CEPHES_H */

