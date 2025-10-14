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
    gretl_matrix *evals;  /* vector of eigenvalues */
    gretl_matrix *Beta;   /* matrix of eigenvectors */
    gretl_matrix *Alpha;  /* matrix of adjustments */
    gretl_matrix *Gamma;  /* I - sum of \Gamma_i matrices */
    gretl_matrix *JC;     /* Johansen's 'C' (Lütkepohl's $\Xi$) */
    gretl_matrix *Bvar;   /* variance matrix of beta */
    gretl_matrix *Bse;    /* standard errors of beta */
    gretl_matrix *Ase;    /* standard errors of alpha */
    gretl_matrix *R;      /* beta-restriction LHS matrix */
    gretl_matrix *q;      /* beta-restrictions RHS matrix */
    gretl_matrix *Ra;     /* alpha-restriction LHS matrix */
    gretl_matrix *qa;     /* alpha-restrictions RHS matrix */
    gretl_matrix *YY;     /* double-size Y matrix */
    gretl_matrix *RR;     /* double-size residuals matrix */
    gretl_matrix *BB;     /* double-size coefficient matrix */
    gretl_matrix **G;     /* array of \Gamma_i matrices */
    double ll0;           /* unrestricted log-likelihood */
    int lrdf;             /* df for likelihood ratio test */
    double prior_ll;      /* ll for prior model in restriction sequence */
    int prior_df;         /* df for prior model in restriction sequence */
};

#define jcode(v) ((v->jinfo == NULL)? 0 : v->jinfo->code)

#define jrank(v) ((v->jinfo == NULL)? 0 : v->jinfo->rank)

/* vecm contains an "automatic" restricted term */
#define auto_restr(v) (v->jinfo != NULL && \
                       (v->jinfo->code == J_REST_CONST || \
                        v->jinfo->code == J_REST_TREND))

/* number of extra terms confined to the cointegrating space */
int n_restricted_terms (const GRETL_VAR *v);

void print_Johansen_test_case (JohansenCode jcode, PRN *prn);

int gretl_VECM_id (GRETL_VAR *vecm);

int *VAR_list_composite (const int *ylist, const int *xlist,
			 const int *rlist);

#endif /* JOHANSEN_H_ */

