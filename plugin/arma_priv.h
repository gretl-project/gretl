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

typedef struct arma_info_ arma_info;

struct arma_info_ {
    int yno;            /* ID of dependent variable */
    ArmaFlags flags;    /* specification flags */
    int *alist;         /* copy of incoming list */
    const char *pqspec; /* input string with specific AR, MA lags */
    char *pmask;        /* specific AR lags included */
    char *qmask;        /* specific MA lags included */
    double ll;          /* log-likelihood */
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
    int T;              /* sample size for estimation */
    double *y;          /* dependent variable (possibly differenced) */
    double *e;          /* forecast errors */
    const double **X;   /* exogenous vars */
    double yscale;      /* scale factor for y */
    int *xlist;         /* list of regressors (ARMAX) */
    gretl_matrix *dX;   /* differenced regressors (ARIMAX) */
    gretl_matrix *G;    /* score matrix */
    gretl_matrix *V;    /* covariance matrix */
    int n_aux;          /* number of auxiliary arrays */
    double **aux;       /* auxiliary arrays */
    PRN *prn;           /* verbose printer */
};

#define arma_has_seasonal(a)   ((a)->flags & ARMA_SEAS)
#define arma_is_arima(a)       ((a)->flags & ARMA_DSPEC)
#define arma_by_x12a(a)        ((a)->flags & ARMA_X12A)
#define arma_exact_ml(a)       ((a)->flags & ARMA_EXACT)
#define arma_using_vech(a)     ((a)->flags & ARMA_VECH)
#define arma_least_squares(a)  ((a)->flags & ARMA_LS)

#define set_arma_has_seasonal(a)  ((a)->flags |= ARMA_SEAS)
#define set_arma_is_arima(a)      ((a)->flags |= ARMA_DSPEC)
#define unset_arma_is_arima(a)    ((a)->flags &= ~ARMA_DSPEC)
#define set_arma_use_vech(a)      ((a)->flags |= ARMA_VECH)
#define set_arma_least_squares(a) ((a)->flags |= ARMA_LS)

#define AR_included(a,i) (a->pmask == NULL || a->pmask[i] == '1')
#define MA_included(a,i) (a->qmask == NULL || a->qmask[i] == '1')

int ma_out_of_bounds (arma_info *ainfo, const double *theta,
		      const double *Theta);

void bounds_checker_cleanup (void);

int hr_arma_init (double *coeff, const double **Z, 
		  const DATAINFO *pdinfo,
		  arma_info *ainfo, int *done);

int ar_arma_init (double *coeff, const double **Z, 
		  const DATAINFO *pdinfo,
		  arma_info *ainfo, MODEL *pmod);

int arma_by_ls (const double *coeff, const double **Z, 
		const DATAINFO *pdinfo,
		arma_info *ainfo, MODEL *pmod);

#endif /* ARMA_PRIV_H */
