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
double ndtri (double);                   /* ndtri.c */
double p1evl (double, double *, int);    /* polevl.c */
double polevl (double, double *, int);   /* polevl.c */

#ifndef INFINITY
extern double INFINITY;
#endif

#ifdef INFINITIES
int isfinite (double);
#endif

#ifndef NAN
extern double NAN;
#endif

#ifdef NANS
int isnan (double);
#endif
 
#endif /* CEPHES_H */

