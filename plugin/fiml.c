/*
 *  Copyright (c) 2004 by Allin Cottrell
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

#include "libgretl.h"
#include "gretl_matrix.h"
#include "gretl_matrix_private.h"

#define FDEBUG

typedef struct fiml_system_ fiml_system;

struct fiml_system_ {
    int n;                  /* number of observations per equation */
    int g;                  /* number of (stochastic) equations */
    int gn;                 /* convenience: g * n = number of obs in stacked vectors */
    int totk;               /* total right-hand side vars */
    int nendo;              /* total number of endogenous vars */
    int nexo;               /* total number of exogenous vars */

    gretl_matrix *uhat;     /* structural-form residuals, all equations */
    gretl_matrix *sigma;    /* cross-equation covariance matrix */
    gretl_matrix *psi;      /* Cholesky decomp of sigma-inverse */

    gretl_matrix *G;        /* Gamma matrix: coeffs for endogenous vars */
    gretl_matrix *B;        /* coeffs for exogenous or predetermined vars */

    gretl_vector *arty;     /* stacked gn-vector: LHS of artificial regression */
    gretl_matrix *artx;     /* RHS: transformed stacked matrix of indep vars */
    gretl_matrix *artb;     /* coefficient vector from artificial regression */

    gretl_matrix *WB1;      /* exog vars times coeffs */
    gretl_matrix *WB2;      /* exog vars times coeffs, times Gamma-inverse */

    gretl_equation_system *sys; /* pointer to "parent" equation system */
};

#if 0

static int 
fill_fiml_udot (fiml_system *fsys, const double **Z, int t1)
{
    int i, t;

    /* stack the g dependent variables vertically */
    for (i=0; i<fsys->g; i++) {
	int k = system_get_depvar(fsys->sys, i);

	for (t=0; t<fsys->n; t++) {
	    int gt = i * fsys->n + t;

	    gretl_vector_set(fsys->udot, gt, Z[k][t + t1]);
	}
    }

    /* subtract ML fitted values to get ML residuals */
    

    return 0;
}

#endif

static void fiml_system_destroy (fiml_system *fsys)
{
    gretl_matrix_free(fsys->uhat);
    gretl_matrix_free(fsys->sigma);
    gretl_matrix_free(fsys->psi);

    gretl_matrix_free(fsys->G);
    gretl_matrix_free(fsys->B);

    gretl_vector_free(fsys->arty);
    gretl_matrix_free(fsys->artx);
    gretl_vector_free(fsys->artb);

    gretl_matrix_free(fsys->WB1);
    gretl_matrix_free(fsys->WB2);

    free(fsys);
}

static fiml_system *fiml_system_new (gretl_equation_system *sys)
{
    fiml_system *fsys;
    int *endog_vars;
    int *exog_vars;

    fsys = malloc(sizeof *fsys);
    if (fsys == NULL) return NULL;

    fsys->sys = sys;

    fsys->g = system_n_equations(sys);
    fsys->n = system_n_obs(sys);
    fsys->gn = fsys->g * fsys->n;
    fsys->totk = system_n_indep_vars(sys);

    endog_vars = system_get_endog_vars(sys);
    exog_vars = system_get_instr_vars(sys);

    fsys->nendo = endog_vars[0];
    fsys->nexo = exog_vars[0];

    fsys->uhat = NULL;
    fsys->sigma = NULL;
    fsys->psi = NULL;

    fsys->G = NULL;
    fsys->B = NULL;

    fsys->arty = NULL;
    fsys->artx = NULL;
    fsys->artb = NULL;

    fsys->WB1 = NULL;
    fsys->WB2 = NULL;

    fsys->uhat = gretl_matrix_alloc(fsys->n, fsys->g);
    fsys->sigma = gretl_matrix_alloc(fsys->g, fsys->g);
    fsys->psi = gretl_matrix_alloc(fsys->g, fsys->g);

    fsys->G = gretl_matrix_alloc(fsys->nendo, fsys->nendo);
    fsys->B = gretl_matrix_alloc(fsys->nexo, fsys->nendo);

    fsys->arty = gretl_column_vector_alloc(fsys->gn);
    fsys->artx = gretl_matrix_alloc(fsys->gn, fsys->totk);
    fsys->artb = gretl_column_vector_alloc(fsys->totk);

    fsys->WB1 = gretl_matrix_alloc(fsys->n, fsys->nendo);
    fsys->WB2 = gretl_matrix_alloc(fsys->n, fsys->nendo);

    if (fsys->uhat == NULL || fsys->sigma == NULL || fsys->psi == NULL ||
	fsys->G == NULL || fsys->B == NULL ||
	fsys->arty == NULL || fsys->artx == NULL || fsys->artb == NULL ||
	fsys->WB1 == NULL || fsys->WB2 == NULL) {
	fiml_system_destroy(fsys);
	fsys = NULL;
    }	

    return fsys;
}

/* use the full residuals matrix to form the cross-equation covariance
   matrix; then invert this and do a Cholesky decomposition
*/ 

static int form_fiml_sigma_and_psi (fiml_system *fsys)
{
    int err;

    err = gretl_matrix_multiply_mod(fsys->uhat, GRETL_MOD_TRANSPOSE,
				    fsys->uhat, GRETL_MOD_NONE,
				    fsys->sigma);

    gretl_matrix_divide_by_scalar(fsys->sigma, fsys->n);

#ifdef FDEBUG
    gretl_matrix_print(fsys->sigma, "fiml sigma", NULL);
#endif

    if (!err) {
	gretl_matrix_copy_values(fsys->psi, fsys->sigma);
	err = gretl_invert_symmetric_matrix(fsys->psi);
    }

    if (!err) {
	err = gretl_matrix_cholesky_decomp(fsys->psi);
	gretl_square_matrix_transpose(fsys->psi);
	gretl_matrix_zero_lower(fsys->psi);
    }

#ifdef FDEBUG
    gretl_matrix_print(fsys->psi, "fiml psi", NULL);
#endif

    return err;
}

static int fiml_transcribe_results (fiml_system *fsys)
{
    int err = 0;

    /* to be written */

    return err;
}

/* form the LHS stacked vector for the artificial regression */

static void fiml_form_depvar (fiml_system *fsys)
{
    double u, p, x;
    int i, j, k, t;

    k = 0;
    for (i=0; i<fsys->g; i++) { /* loop across equations */
	for (t=0; t<fsys->n; t++) { /* loop across obs */
	    x = 0.0;
	    for (j=0; j<fsys->g; j++) {
		p = gretl_matrix_get(fsys->psi, i, j);
		u = gretl_matrix_get(fsys->uhat, t, j);
		x += p * u;
	    }
	    gretl_vector_set(fsys->arty, k++, x);
	}
    }

#ifdef FDEBUG
    gretl_matrix_print(fsys->arty, "fiml artificial Y", NULL);
#endif
}

static int on_exo_list (const int *exlist, int v)
{
    int i;

    for (i=1; i<=exlist[0]; i++) {
	if (exlist[i] == v) return 1;
    }

    return 0;
}

static int endo_var_number (const int *enlist, int v)
{
    int i;

    for (i=1; i<=enlist[0]; i++) {
	if (enlist[i] == v) return i - 1;
    }

    return -1;
}

/* form the RHS matrix for the artificial regression */

static void 
fiml_form_indepvars (fiml_system *fsys, const double **Z, int t1)
{
    const int *enlist = system_get_endog_vars(fsys->sys);
    const int *exlist = system_get_instr_vars(fsys->sys);
    int i, j, k, t;
    int bigrow, bigcol = 0;
    double p;

    gretl_matrix_zero(fsys->artx);

    for (i=0; i<fsys->g; i++) { /* loop across equations */
	const int *list = system_get_list(fsys->sys, i);

	for (j=2; j<=list[0]; j++) { /* loop across RHS vars */
	    const double *xj = NULL;
	    double xjt;
	    int vj = 0;

	    if (on_exo_list(exlist, list[j])) {
		xj = Z[list[j]] + t1;
	    } else {
		vj = endo_var_number(enlist, list[j]);
	    }

	    for (t=0; t<fsys->n; t++) { /* loop across obs */
		for (k=0; k<fsys->g; k++) { /* loop across vertical blocks */
		    bigrow = k * fsys->n + t;
		    p = gretl_matrix_get(fsys->psi, k, i);
		    if (p == 0.0) {
			gretl_matrix_set(fsys->artx, bigrow, bigcol, 0.0);
		    } else {
			if (xj != NULL) {
			    xjt = xj[t];
			} else {
			    xjt = gretl_matrix_get(fsys->WB2, t, vj);
			}
			gretl_matrix_set(fsys->artx, bigrow, bigcol, xjt * p);
		    }
		}
	    }
	    bigcol++;
	}
    }

#ifdef FDEBUG
    gretl_matrix_print(fsys->artx, "fiml artificial X", NULL);
#endif
}

/* initialize residuals matrix based on 3SLS */

static int fiml_uhat_init (fiml_system *fsys)
{
    const gretl_matrix *uhat = system_get_uhat(fsys->sys);
    double x;
    int i, t;

    for (i=0; i<fsys->g; i++) {
	for (t=0; t<fsys->n; t++) {
	    x = gretl_matrix_get(uhat, i, t);
	    gretl_matrix_set(fsys->uhat, t, i, x);
	}
    }

    return 0;
}

static int 
rhs_var_in_eqn (const gretl_equation_system *sys, int eq, int v)
{
    const int *list = system_get_list(sys, eq);

    if (list != NULL) {
	int i;

	for (i=2; i<=list[0]; i++) {
	    if (list[i] == v) return i;
	}
    }

    return 0;
}

static int fiml_G_init (fiml_system *fsys, const DATAINFO *pdinfo)
{
    const int *enlist = system_get_endog_vars(fsys->sys);
    MODEL *pmod;
    int lv, rv;
    int i, j, vi;

    for (j=0; j<fsys->nendo; j++) {
	/* outer loop across columns (equations) */
	lv = enlist[j + 1];
#ifdef FDEBUG
	printf("G: working on equation %d (%s, %s)\n", j+1, 
	       pdinfo->varname[lv], (j < fsys->g)? "stochastic" : "identity");
#endif
	for (i=0; i<fsys->nendo; i++) {
	    rv = enlist[i + 1];
	    if (j == i) {
		gretl_matrix_set(fsys->G, i, j, 1.0);
	    } else if (j < fsys->g) {
		/* column pertains to stochastic equation */
		vi = rhs_var_in_eqn(fsys->sys, j, rv);
		if (vi > 0) {
#ifdef FDEBUG
		    printf(" endog var %s is included\n", pdinfo->varname[rv]);
#endif
		    pmod = system_get_model(fsys->sys, j);
		    gretl_matrix_set(fsys->G, i, j, -pmod->coeff[vi-2]);
		} else {
#ifdef FDEBUG
		    printf(" endog var %s is excluded\n", pdinfo->varname[rv]);
#endif
		    gretl_matrix_set(fsys->G, i, j, 0.0);
		}
	    } else {
		/* column pertains to identity */
		vi = -1 * rhs_var_in_identity(fsys->sys, lv, rv);
#ifdef FDEBUG
		if (vi) {
		    printf(" var %s has coeff %d\n", pdinfo->varname[rv], vi);
		}
#endif
		gretl_matrix_set(fsys->G, i, j, vi);
	    }
	}
    }

#ifdef FDEBUG
    gretl_matrix_print(fsys->G, "fiml G", NULL);
#endif

    return 0;
}

static int fiml_B_init (fiml_system *fsys, const DATAINFO *pdinfo)
{
    const int *enlist = system_get_endog_vars(fsys->sys);
    const int *exlist = system_get_instr_vars(fsys->sys);
    MODEL *pmod;
    int lv, rv;
    int i, j, vi;

    for (j=0; j<fsys->nendo; j++) {
	lv = enlist[j + 1];
	/* outer loop across columns (equations) */
#ifdef FDEBUG
	printf("B: working on equation %d (%s, %s)\n", j+1, 
	       pdinfo->varname[lv], (j < fsys->g)? "stochastic" : "identity");
#endif
	for (i=0; i<fsys->nexo; i++) {
	    rv = exlist[i + 1];
	    if (j < fsys->g) {
		/* column pertains to stochastic equation */
		vi = rhs_var_in_eqn(fsys->sys, j, rv);
		if (vi > 0) {
#ifdef FDEBUG
		    printf(" exog var %s is included\n", pdinfo->varname[rv]);
#endif
		    pmod = system_get_model(fsys->sys, j);
		    gretl_matrix_set(fsys->B, i, j, pmod->coeff[vi-2]);
		} else {
#ifdef FDEBUG
		    printf(" exog var %s is excluded\n", pdinfo->varname[rv]);
#endif
		    gretl_matrix_set(fsys->B, i, j, 0.0);
		}
	    } else {
		vi = rhs_var_in_identity(fsys->sys, lv, rv);
#ifdef FDEBUG
		if (vi) {
		    printf(" var %s has coeff %d\n", pdinfo->varname[rv], vi);
		}
#endif
		gretl_matrix_set(fsys->B, i, j, vi);
	    }
	}
    }

#ifdef FDEBUG
    gretl_matrix_print(fsys->B, "fiml B", NULL);
#endif

    return 0;
}

#define LN_2_PI 1.837877066409345

static double fiml_ll (fiml_system *fsys)
{
    double ll = 0.0;

    /* should copy these matrices first? */
    double detG = gretl_LU_determinant(fsys->G);
    double detS = gretl_LU_determinant(fsys->sigma);

    ll -= (fsys->gn / 2.0) * LN_2_PI;
    ll += fsys->n * log(fabs(detG));
    ll -= (fsys->n / 2.0) * log(detS);
    
    /* now for the other bit, involving the FIML residuals... */

    return ll;
}

static int fiml_endog_rhs (fiml_system *fsys, const double **Z, int t1)
{
    double x;
    int i, j, t;
    int err;

    for (j=0; j<fsys->nendo; j++) {
	for (t=0; t<fsys->n; t++) {
	    x = 0.0;
	    for (i=0; i<fsys->nexo; i++) {
		x += Z[i][t + t1] * gretl_matrix_get(fsys->B, i, j);
	    } 
	    gretl_matrix_set(fsys->WB1, t, j, x);
	}
    }

    err = gretl_invert_general_matrix(fsys->G);

#ifdef FDEBUG
    if (!err) {
	gretl_matrix_print(fsys->G, "G-inverse", NULL);
    }
#endif

    if (!err) {
	gretl_matrix_multiply(fsys->WB1, fsys->G, fsys->WB2);
    }

    if (err) {
	fprintf(stderr, "fiml_endog_rhs failed\n");
    } 

    return err;
}

/* Below: this is meant to grow into a driver function for FIML as
   described in Davidson and MacKinnon, ETM, chap 12, sect 5.
*/

int fiml_driver (gretl_equation_system *sys, const double **Z, const DATAINFO *pdinfo)
{
    fiml_system *fsys;
    int i, err = 0;

    fsys = fiml_system_new(sys);
    if (fsys == NULL) {
	return E_ALLOC;
    }

    /* initialize uhat based on 3SLS estimates */
    fiml_uhat_init(fsys);

    /* intialize Gamma coefficient matrix */
    fiml_G_init(fsys, pdinfo);

    /* intialize B coefficient matrix */
    fiml_B_init(fsys, pdinfo);

    /* FIXME: need to finish writing the following loop, with
       appropriate convergence/fail conditions */

    while (1) {

	/* form \hat{Sigma} (ETM, equation 12.81), invert, 
	   and Cholesky-decompose to get \Psi */
	form_fiml_sigma_and_psi(fsys);

	/* form artificial dependent variable */
	fiml_form_depvar(fsys);

	/* instrument the RHS endog vars */
	fiml_endog_rhs(fsys, Z, pdinfo->t1);

	/* form artificial matrix of independent variables */
	fiml_form_indepvars(fsys, Z, pdinfo->t1);

	/* run artificial regression (ETM, equation 12.86) */
	gretl_matrix_ols(fsys->arty, fsys->artx, fsys->artb, NULL, NULL);
	for (i=0; i<fsys->totk; i++) {
	    printf("fsys->artb[%d] = %g\n", i, gretl_vector_get(fsys->artb, i));
	}

	/* now adjust parameter estimates based on fsys->artb,
	   recompute fiml residuals, and go back */

	break;

    }

#ifdef notyet
    /* write the results into the parent system */
    fiml_transcribe_results(fsys);
#endif
    
    /* clean up */
    fiml_system_destroy(fsys);

    return err;
}



    
	

    
    
    

    
    
