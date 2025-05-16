/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#ifndef ARMA_PRIV_H
#define ARMA_PRIV_H

typedef enum {
    ARMA_SEAS   = 1 << 0, /* includes seasonal component */
    ARMA_DSPEC  = 1 << 1, /* input list includes differences */
    ARMA_XDIFF  = 1 << 2, /* ARIMA: exogenous regressors are differenced */
    ARMA_LBFGS  = 1 << 3, /* using L-BFGS-B with native exact ML */
    ARMA_VECH   = 1 << 4, /* using vech representation when computing
			     variance matrix of state for Kalman filter */
    ARMA_NAOK   = 1 << 5, /* allow missing observations */
    ARMA_NAS    = 1 << 6, /* sample contains NAs */
    ARMA_LEV    = 1 << 7, /* doing ARIMA via levels formulation */
    ARMA_YDIFF  = 1 << 8, /* ainfo->y contains differenced y */
    ARMA_AVGLL  = 1 << 9, /* passing average likelihood option to Kalman */
    ARMA_STDX   = 1 << 10 /* exogenous regressors standardized */
} PrivFlags;

typedef enum {
    INI_USER = 1,
    INI_HR,
    INI_SMALL,
    INI_NLS,
    INI_OLS
} IniMethod;

typedef struct arma_info_ arma_info;

struct arma_info_ {
    int yno;            /* ID of dependent variable */
    ArmaFlags flags;    /* specification flags */
    PrivFlags pflags;   /* "private" flags for estimation */
    int *alist;         /* copy of incoming list */
    const int *pqspec;  /* auxiliary list with specific AR, MA lags */
    char *pmask;        /* specific AR lags included */
    char *qmask;        /* specific MA lags included */
    double ll;          /* log-likelihood */
    IniMethod init;     /* initialization method */
    int ifc;            /* specification includes a constant? */
    int p;              /* max non-seasonal AR order */
    int d;              /* non-seasonal difference */
    int q;              /* max non-seasonal MA order */
    int P;              /* seasonal AR order */
    int D;              /* seasonal difference */
    int Q;              /* seasonal MA order */
    int np;             /* total non-seasonal AR lags */
    int nq;             /* total non-seasonal MA lags */
    int maxlag;         /* longest lag in model */
    int nexo;           /* number of other regressors (ARMAX) */
    int nc;             /* total number of coefficients */
    int t1;             /* starting observation */
    int t2;             /* ending observation */
    int pd;             /* periodicity of data */
    int T;              /* number of valid observations in sample */
    int fullT;          /* total obs (possibly including interior NAs) */
    int r0;             /* size of ARMA state vector */
    int fncount;        /* count of function evaluations */
    int grcount;        /* count of gradient evaluations */
    double *y;          /* dependent variable (possibly differenced) */
    double *e;          /* forecast errors */
    const double **Z;   /* virtual dataset */
    double yscale;      /* scale factor for y */
    double yshift;      /* shift factor for y */
    int *xlist;         /* list of regressors (ARMAX) */
    int *misslist;      /* list of missing observations */
    gretl_matrix *xstats; /* mean and std dev of regressors */
    gretl_matrix *dX;   /* differenced regressors (ARIMAX) */
    gretl_matrix *G;    /* score matrix */
    gretl_matrix *V;    /* covariance matrix */
    int n_aux;          /* number of auxiliary arrays */
    double **aux;       /* auxiliary arrays */
    PRN *prn;           /* verbose printer */
};

#define arma_by_x12a(a)        ((a)->flags & ARMA_X12A)
#define arma_exact_ml(a)       ((a)->flags & ARMA_EXACT)
#define arma_least_squares(a)  ((a)->flags & ARMA_LS)

#define set_arma_least_squares(a) ((a)->flags |= ARMA_LS)

#define arma_has_seasonal(a)   ((a)->pflags & ARMA_SEAS)
#define arma_is_arima(a)       ((a)->pflags & ARMA_DSPEC)
#define arma_xdiff(a)          ((a)->pflags & ARMA_XDIFF)
#define arma_lbfgs(a)          ((a)->pflags & ARMA_LBFGS)
#define arma_using_vech(a)     ((a)->pflags & ARMA_VECH)
#define arma_na_ok(a)          ((a)->pflags & ARMA_NAOK)
#define arma_missvals(a)       ((a)->pflags & ARMA_NAS)
#define arima_levels(a)        ((a)->pflags & ARMA_LEV)
#define arima_ydiff(a)         ((a)->pflags & ARMA_YDIFF)
#define arma_avg_ll(a)         ((a)->pflags & ARMA_AVGLL)
#define arma_stdx(a)           ((a)->pflags & ARMA_STDX)    

#define set_arma_has_seasonal(a)  ((a)->pflags |= ARMA_SEAS)
#define set_arma_is_arima(a)      ((a)->pflags |= ARMA_DSPEC)
#define unset_arma_is_arima(a)    ((a)->pflags &= ~ARMA_DSPEC)
#define set_arma_use_vech(a)      ((a)->pflags |= ARMA_VECH)
#define set_arma_na_ok(a)         ((a)->pflags |= ARMA_NAOK)
#define set_arma_missvals(a)      ((a)->pflags |= ARMA_NAS)
#define set_arima_levels(a)       ((a)->pflags |= ARMA_LEV)
#define set_arma_avg_ll(a)        ((a)->pflags |= ARMA_AVGLL)
#define set_arima_ydiff(a)        ((a)->pflags |= ARMA_YDIFF)
#define unset_arima_ydiff(a)      ((a)->pflags &= ~ARMA_YDIFF)
#define set_arma_stdx(a)          ((a)->pflags |= ARMA_STDX)

#define AR_included(a,i) (a->pmask == NULL || a->pmask[i] == '1')
#define MA_included(a,i) (a->qmask == NULL || a->qmask[i] == '1')

int flip_poly (double *theta, arma_info *ainfo,
	       int ar, int seasonal);

int maybe_correct_MA (arma_info *ainfo,
		      double *theta,
		      double *Theta);

void maybe_set_yscale (arma_info *ainfo);

int hr_arma_init (double *coeff, const DATASET *dset,
		  arma_info *ainfo);

int ar_arma_init (double *coeff, const DATASET *dset,
		  arma_info *ainfo, MODEL *pmod,
		  gretlopt opt);

int arma_by_simple_ols (const double *coeff, const DATASET *dset,
                        arma_info *ainfo, MODEL *pmod);

int arma_by_ls (const double *coeff, const DATASET *dset,
		arma_info *ainfo, MODEL *pmod);

int bhhh_arma (double *theta, const DATASET *dset,
	       arma_info *ainfo, MODEL *pmod,
	       gretlopt opt);

void transform_arma_const (double *b, arma_info *ainfo);

int arma_list_y_position (arma_info *ainfo);

int arma_model_add_roots (MODEL *pmod, arma_info *ainfo,
			  const double *coeff);

void write_arma_model_stats (MODEL *pmod, arma_info *ainfo,
			     const DATASET *dset);

int arima_difference (arma_info *ainfo, const DATASET *dset, 
		      int fullX);

void arima_difference_undo (arma_info *ainfo, const DATASET *dset);

#endif /* ARMA_PRIV_H */
