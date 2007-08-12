/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 2005 Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef JOHANSEN_H_
#define JOHANSEN_H_

#include "gretl_matrix.h"

typedef enum {
    J_NO_CONST = 0,
    J_REST_CONST,
    J_UNREST_CONST,
    J_REST_TREND,
    J_UNREST_TREND
} JohansenCode;

typedef struct JohansenInfo_ JohansenInfo;

struct JohansenInfo_ {
    int ID;               /* for identifying saved vars */
    JohansenCode code;    /* see enumeration above */
    int rank;             /* if specified, chosen cointegration rank, else 0 */
    int seasonals;        /* number of seasonal dummies included */
    gretl_matrix *R0;     /* residuals, VAR in differences */
    gretl_matrix *R1;     /* residuals, second regressions */
    gretl_matrix *S00;    /* cross-products of residuals */
    gretl_matrix *S11;    /* cross-products of residuals */
    gretl_matrix *S01;    /* cross-products of residuals */
    gretl_matrix *Beta;   /* matrix of eigenvectors */
    gretl_matrix *Alpha;  /* matrix of adjustments */
    gretl_matrix *Bvar;   /* variance matrix of beta */
    gretl_matrix *Bse;    /* standard errors of beta */
    gretl_matrix *Ase;    /* standard errors of alpha */
    gretl_matrix *R;      /* beta-restriction LHS matrix */
    gretl_matrix *q;      /* beta-restrictions RHS matrix */
    gretl_matrix *Ra;     /* alpha-restriction LHS matrix */
    gretl_matrix *qa;     /* alpha-restrictions RHS matrix */
    double ll0;           /* unrestricted log-likelihood */
    int lrdf;             /* df for likelihood ratio test */
};

struct GRETL_VAR_ {
    int ci;              /* command index (VAR or VECM) */
    int refcount;        /* for saving/deleting */
    int err;             /* error code */
    int neqns;           /* number of equations in system */
    int order;           /* lag order */
    int t1;              /* starting observation */
    int t2;              /* ending observation */
    int T;               /* number of observations */
    int ifc;             /* equations include a constant (1) or not (0) */
    int ncoeff;          /* total coefficients per equation */
    int *ylist;          /* list of stochastic vars */
    int *xlist;          /* list of exogenous variables */
    int detflags;        /* record of automatic deterministic vars added */
    int robust;          /* computing robust std errors? */
    int qr;              /* using QR decomposition? */
    gretl_matrix *Y;     /* matrix of dependent variables */
    gretl_matrix *X;     /* matrix of independent variables */
    gretl_matrix *B;     /* basic coefficient matrix */
    gretl_matrix *XTX;   /* X'X */
    gretl_matrix *A;     /* augmented coefficient matrix (companion form) */
    gretl_matrix *L;     /* lambda: inverse roots of A(L) polynomial */
    gretl_matrix *E;     /* residuals matrix */
    gretl_matrix *C;     /* augmented Cholesky-decomposed error matrix */
    gretl_matrix *S;     /* cross-equation variance matrix */
    gretl_matrix *F;     /* optional forecast matrix */
    MODEL **models;      /* pointers to individual equation estimates */
    double *Fvals;       /* hold results of F-tests */
    double *Ivals;       /* hold results of info criteria comparisons */
    double ldet;         /* log-determinant of S */
    double ll;           /* log-likelihood */
    double AIC;          /* Akaike criterion */
    double BIC;          /* Bayesian criterion */
    double HQC;          /* Hannan-Quinn criterion */
    double LR;           /* for likelihood-ratio testing */
    double LB;           /* Ljung-Box (Portmanteau) test statistic */
    int LBs;             /* order for for Portmanteau test */
    JohansenInfo *jinfo; /* extra information for VECMs */
    char *name;          /* for use in session management */
};
    
#define restricted(v) (v->jinfo != NULL && \
                       (v->jinfo->code == J_REST_CONST || \
                        v->jinfo->code == J_REST_TREND))

#define jcode(v) (v->jinfo->code)
#define jrank(v) (v->jinfo->rank)

#define effective_order(v) (v->order+(v->ci==VECM))

GRETL_VAR *johansen_test (int order, const int *list, 
			  const double **Z, const DATAINFO *pdinfo,
			  gretlopt opt, PRN *prn);

int johansen_test_simple (int order, const int *list, 
			  const double **Z, const DATAINFO *pdinfo,
			  gretlopt opt, PRN *prn);

void print_Johansen_test_case (JohansenCode jcode, PRN *prn);

int gretl_VECM_id (GRETL_VAR *vecm);

int 
gretl_VAR_do_error_decomp (const gretl_matrix *S, gretl_matrix *C);


#endif /* JOHANSEN_H_ */

