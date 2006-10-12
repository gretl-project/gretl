#ifndef MOD_GARCH_H
#define MOD_GARCH_H

enum {
    INIT_VAR_THEO,
    INIT_VAR_OLS,
    INIT_VAR_RESID
};

enum {
    DIST_NORM,
    DIST_T
};

typedef struct garch_container_ garch_container;
struct garch_container_ {
    double *y;
    double **X;
    int t1;
    int t2;
    int nobs;
    int ncm;
    int p;
    int q;
    int initmeth;
    int distrib;
    double *e;
    double *e2;
    double *h;
    int ascore;
    double **score_e;
    double **score_h;
    double **blockglue;
    double **G;
    double *tot_score;
};


int garch_estimate_mod (int t1, int t2, int nobs, 
			double **X, int nx, double *coeff, int nc, 
			double *vcv, double *res2, double *res, double *h,
			double *y, double *amax, double *b, 
			double scale, int *iters, PRN *prn, int vopt);

#endif /* MOD_GARCH_H */