/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef GRETL_SYSTEM_PRIVATE_H
#define GRETL_SYSTEM_PRIVATE_H

enum {
    GRETL_SYSTEM_SAVE_UHAT = 1 << 0,
    GRETL_SYSTEM_SAVE_YHAT = 1 << 1,
    GRETL_SYSTEM_DFCORR    = 1 << 2,
    GRETL_SYS_VCV_GEOMEAN  = 1 << 3,
    GRETL_SYS_SAVE_VCV     = 1 << 4,
    GRETL_SYS_RESTRICT     = 1 << 5,
    GRETL_SYS_ITERATE      = 1 << 6
};

typedef struct id_atom_ id_atom;
typedef struct identity_ identity;

struct _gretl_equation_system {
    char *name;                 /* user-specified name for system, or NULL */
    int t1;                     /* starting observation number */
    int t2;                     /* ending observation number */
    int method;                 /* estimation method */
    int n_equations;            /* number of stochastic equations */
    int n_identities;           /* number of identities */
    int n_obs;                  /* number of observations per equation */
    int iters;                  /* number of iterations taken */
    char flags;                 /* to record options (e.g. save residuals) */
    double ll;                  /* log-likelihood (restricted) */
    double llu;                 /* unrestricted log-likelihood */
    double X2;                  /* chi-square test value */
    double ess;                 /* total error sum of squares */
    double diag;                /* test stat for diagonal covariance matrix */
    double bdiff;               /* summary stat for change in coefficients */
    int **lists;                /* regression lists for stochastic equations */
    int *endog_vars;            /* list of endogenous variables */
    int *instr_vars;            /* list of instruments (exogenous vars) */
    identity **idents;          /* set of identities */
    gretl_matrix *b;            /* coefficient estimates */
    gretl_matrix *vcv;          /* covariance matrix of coefficients */
    gretl_matrix *sigma;        /* cross-equation covariance matrix */
    gretl_matrix *R;            /* LHS of any linear restrictions */
    gretl_matrix *q;            /* RHS of any linear restrictions */  
    gretl_matrix *uhat;         /* residuals, all equations */
    MODEL **models;             /* set of pointers to per-equation models: just
				   convenience pointers -- these should NOT be
				   freed as part of sys cleanup
				*/
};

void make_system_data_info (gretl_equation_system *sys, int eqn, 
			    DATAINFO *pdinfo, int v, int code);

#endif /* GRETL_SYSTEM_PRIVATE_H */
