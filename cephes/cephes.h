#ifndef CEPHES_H
#define CEPHES_H

extern double MACHEP;
extern double MINLOG;
extern double MAXLOG;
extern double MAXNUM;

extern double PI, PIO4;
extern double SQRTH;
 
double expm1 (double);                   /* unity.c */
double expx2 (double, int);              /* expx2.c */
double gamma (double);                   /* gamma.c */
double igam (double, double);            /* igam.c */
double igamc (double, double);           /* igam.c */
double igami (double, double);           /* igami.c */
double incbet (double, double, double);  /* incbet.c */
double incbi (double, double, double);   /* incbi.c */
double lgam (double);                    /* gamma.c */
double log1p (double);                   /* unity.c */
double ndtri (double);                   /* ndtri.c */
double p1evl (double, double *, int);      /* polevl.c */
double polevl (double, double *, int);     /* polevl.c */

#ifdef INFINITIES
extern double INFINITY;
int isfinite (double);
#endif

#ifdef NANS
extern double NAN;
int isnan (double);
#endif
 
#endif /* CEPHES_H */

