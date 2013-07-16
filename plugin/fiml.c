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

#include "libgretl.h"
#include "gretl_matrix.h"
#include "system.h"
#include "sysml.h"

#define FDEBUG 0

typedef struct fiml_system_ fiml_system;

struct fiml_system_ {
    int n;                /* number of observations per equation */
    int g;                /* number of (stochastic) equations */
    int gn;               /* g * n = number of obs in stacked vectors */
    int totk;             /* total right-hand side vars */
    int nendo;            /* total number of endogenous vars */
    int nexo;             /* total number of exogenous vars */

    double ll;            /* log-likelihood */
    double llu;           /* unrestricted log-likelihood */

    gretl_matrix *uhat;   /* structural-form residuals, all equations */
    gretl_matrix *sigma;  /* cross-equation covariance matrix */
    gretl_matrix *psi;    /* Cholesky decomp of sigma-inverse */
    gretl_matrix *Stmp;   /* workspace */

    gretl_matrix *G;      /* Gamma matrix: coeffs for endogenous vars */
    gretl_matrix *B;      /* coeffs for exogenous and predetermined vars */
    gretl_matrix *Gtmp;   /* workspace */

    gretl_vector *arty;   /* stacked gn-vector: LHS of artificial regression */
    gretl_matrix *artx;   /* stacked matrix of transformed indep vars: RHS */
    gretl_matrix *artb;   /* coefficient vector from artificial regression */
    gretl_matrix *btmp;   /* workspace */

    gretl_matrix *WB1;    /* exog vars times coeffs */
    gretl_matrix *WB2;    /* exog vars times coeffs, times Gamma-inverse */

    equation_system *sys; /* pointer to "parent" equation system */
};

static void fiml_system_destroy (fiml_system *fsys)
{
    gretl_matrix_free(fsys->uhat);
    gretl_matrix_free(fsys->sigma);
    gretl_matrix_free(fsys->psi);
    gretl_matrix_free(fsys->Stmp);

    gretl_matrix_free(fsys->G);
    gretl_matrix_free(fsys->B);
    gretl_matrix_free(fsys->Gtmp);

    gretl_vector_free(fsys->arty);
    gretl_matrix_free(fsys->artx);
    gretl_vector_free(fsys->artb);
    gretl_vector_free(fsys->btmp);

    gretl_matrix_free(fsys->WB1);
    gretl_matrix_free(fsys->WB2);

    free(fsys);
}

static fiml_system *fiml_system_new (equation_system *sys, int *err)
{
    fiml_system *fsys;
    int *endog_vars;
    int *exog_vars;

    endog_vars = system_get_endog_vars(sys);
    exog_vars = system_get_instr_vars(sys);

    if (endog_vars == NULL || exog_vars == NULL) {
	gretl_errmsg_set(_("No list of endogenous variables was given"));
	*err = E_DATA;
	return NULL;
    }

    fsys = malloc(sizeof *fsys);
    if (fsys == NULL) return NULL;

    fsys->sys = sys;

    fsys->g = sys->neqns;
    fsys->n = sys->T;
    fsys->gn = fsys->g * fsys->n;
    fsys->totk = system_n_indep_vars(sys);

    fsys->nendo = endog_vars[0];
    fsys->nexo = exog_vars[0];

    fsys->ll = 0.0;
    fsys->llu = 0.0;

    fsys->uhat = NULL;
    fsys->sigma = NULL;
    fsys->psi = NULL;
    fsys->Stmp = NULL;

    fsys->G = NULL;
    fsys->B = NULL;
    fsys->Gtmp = NULL;

    fsys->arty = NULL;
    fsys->artx = NULL;
    fsys->artb = NULL;
    fsys->btmp = NULL;

    fsys->WB1 = NULL;
    fsys->WB2 = NULL;

    clear_gretl_matrix_err();

    fsys->uhat = gretl_matrix_alloc(fsys->n, fsys->g);
    fsys->sigma = gretl_matrix_alloc(fsys->g, fsys->g);
    fsys->psi = gretl_matrix_alloc(fsys->g, fsys->g);
    fsys->Stmp = gretl_matrix_alloc(fsys->g, fsys->g);

    fsys->G = gretl_matrix_alloc(fsys->nendo, fsys->nendo);
    fsys->B = gretl_matrix_alloc(fsys->nexo, fsys->nendo);
    fsys->Gtmp = gretl_matrix_alloc(fsys->nendo, fsys->nendo);

    fsys->arty = gretl_column_vector_alloc(fsys->gn);
    fsys->artx = gretl_matrix_alloc(fsys->gn, fsys->totk);
    fsys->artb = gretl_column_vector_alloc(fsys->totk);
    fsys->btmp = gretl_column_vector_alloc(fsys->totk);

    fsys->WB1 = gretl_matrix_alloc(fsys->n, fsys->nendo);
    fsys->WB2 = gretl_matrix_alloc(fsys->n, fsys->nendo);

    if (get_gretl_matrix_err()) {
	fiml_system_destroy(fsys);
	fsys = NULL;
    }	

    return fsys;
}

/* estimate the unrestricted reduced-form equations to get the
   unrestricted log-likelihood for the system 
*/

static int fiml_overid_test (fiml_system *fsys, DATASET *dset)
{
    const int *enlist = system_get_endog_vars(fsys->sys);
    const int *exlist = system_get_instr_vars(fsys->sys);
    int t1 = dset->t1;
    gretl_matrix *uru = NULL;
    gretl_matrix *urv = NULL;
    MODEL umod;
    double ldetS;
    int *list;
    int i, t;
    int err = 0;

    if (system_get_overid_df(fsys->sys) <= 0) {
	return 1;
    }

    list = malloc((fsys->nexo + 2) * sizeof *list);
    if (list == NULL) {
	return E_ALLOC;
    }

    uru = gretl_matrix_alloc(fsys->n, fsys->g);
    if (uru == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    urv = gretl_matrix_alloc(fsys->g, fsys->g);
    if (urv == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    list[0] = fsys->nexo + 1;
    for (i=2; i<=list[0]; i++) {
	list[i] = exlist[i - 1];
    }

    for (i=0; i<fsys->g; i++) {
	list[1] = enlist[i + 1];
	umod = lsq(list, dset, OLS, OPT_A);
	if (umod.errcode) {
	    err = umod.errcode;
	    goto bailout;
	}
	for (t=0; t<fsys->n; t++) {
	    gretl_matrix_set(uru, t, i, umod.uhat[t + t1]);
	}
	clear_model(&umod);
    }

    err = gretl_matrix_multiply_mod(uru, GRETL_MOD_TRANSPOSE,
				    uru, GRETL_MOD_NONE,
				    urv, GRETL_MOD_NONE);

    if (err) {
	goto bailout;
    }

    gretl_matrix_divide_by_scalar(urv, fsys->n);
    ldetS = gretl_matrix_log_determinant(urv, &err);
    if (na(ldetS)) {
	goto bailout;
    }

    fsys->llu = - (fsys->gn / 2.0) * (LN_2_PI + 1.0);
    fsys->llu -= (fsys->n / 2.0) * ldetS;

 bailout:

    gretl_matrix_free(uru);
    gretl_matrix_free(urv);
    free(list);
    
    return err;
}

/* calculate FIML residuals as YG - WB */

static void fiml_form_uhat (fiml_system *fsys, 
			    const DATASET *dset, 
			    int t1)
{
    const int *enlist = system_get_endog_vars(fsys->sys);
    const int *exlist = system_get_instr_vars(fsys->sys);
    double bij, gij;
    double y, x;
    int i, j, t;

    for (j=0; j<fsys->nendo; j++) {
	for (t=0; t<fsys->n; t++) {
	    y = 0.0;
	    for (i=0; i<fsys->nendo; i++) {
		gij = gretl_matrix_get(fsys->G, i, j);
		y += dset->Z[enlist[i + 1]][t + t1] * gij;
	    }
	    x = 0.0;
	    for (i=0; i<fsys->nexo; i++) {
		bij = gretl_matrix_get(fsys->B, i, j);
		x += dset->Z[exlist[i + 1]][t + t1] * bij;
	    } 
	    gretl_matrix_set(fsys->WB1, t, j, x);
	    if (j < fsys->g) {
		gretl_matrix_set(fsys->uhat, t, j, y - x);
	    }	
	}
    }

#if FDEBUG
    gretl_matrix_print(fsys->uhat, "fiml uhat");
#endif
}

/* use the full residuals matrix to form the cross-equation covariance
   matrix; then invert this and do a Cholesky decomposition to find
   psi-transpose
*/ 

static int 
fiml_form_sigma_and_psi (fiml_system *fsys, const DATASET *dset, 
			 int t1)
{
    int err;

    /* YG - WB */
    fiml_form_uhat(fsys, dset, t1);

    /* Davidson and MacKinnon, ETM, equation (12.81) */

    err = gretl_matrix_multiply_mod(fsys->uhat, GRETL_MOD_TRANSPOSE,
				    fsys->uhat, GRETL_MOD_NONE,
				    fsys->sigma, GRETL_MOD_NONE);

    gretl_matrix_divide_by_scalar(fsys->sigma, fsys->n);

#if FDEBUG
    gretl_matrix_print(fsys->sigma, "fiml Sigma");
#endif

    if (!err) {
	gretl_matrix_copy_values(fsys->psi, fsys->sigma);
	err = gretl_invert_symmetric_matrix(fsys->psi);
    }

#if FDEBUG
    gretl_matrix_print(fsys->psi, "Sigma-inverse");
#endif

    if (!err) {
	err = gretl_matrix_cholesky_decomp(fsys->psi);
	/* we actually want the transpose of psi (ETM, under eq (12.86) */
	gretl_square_matrix_transpose(fsys->psi);
    }

#if FDEBUG
    gretl_matrix_print(fsys->psi, "fiml Psi-transpose");
#endif

    return err;
}

static void 
fiml_transcribe_results (fiml_system *fsys, const DATASET *dset, 
			 int t1, int iters)
{
    MODEL *pmod;
    const double *y;
    double u;
    int i, j, k, t;

    /* correct uhat and yhat; correct ESS/SSR and standard error,
       per equation; update vectorized coeffs, b */

    k = 0;
    for (i=0; i<fsys->g; i++) {
	pmod = system_get_model(fsys->sys, i);
	y = dset->Z[pmod->list[1]];
	pmod->ess = 0.0;
	for (t=0; t<fsys->n; t++) {
	    u = gretl_matrix_get(fsys->uhat, t, i);
	    pmod->uhat[t + t1] = u;
	    pmod->yhat[t + t1] = y[t + t1] - u;
	    pmod->ess += u * u;
	}
	pmod->sigma = sqrt(pmod->ess / pmod->nobs);
	for (j=0; j<pmod->ncoeff; j++) {
	    fsys->sys->b->val[k++] = pmod->coeff[j];
	}
    }

    /* not using df correction for pmod->sigma or sigma matrix */
    system_attach_sigma(fsys->sys, fsys->sigma);
    fsys->sigma = NULL;

    system_attach_uhat(fsys->sys, fsys->uhat);
    fsys->uhat = NULL;

    /* record restricted and unrestricted log-likelihood */
    fsys->sys->ll = fsys->ll;
    fsys->sys->llu = fsys->llu;

    /* record number of iterations taken */
    fsys->sys->iters = iters;
}

/* form the LHS stacked vector for the artificial regression */

static void fiml_form_depvar (fiml_system *fsys)
{
    double u, p, x;
    int i, j, k, t;

    k = 0;
    for (i=0; i<fsys->g; i++) { 
	/* loop across equations */
	for (t=0; t<fsys->n; t++) { 
	    /* loop across observations */
	    x = 0.0;
	    for (j=0; j<fsys->g; j++) {
		p = gretl_matrix_get(fsys->psi, i, j);
		u = gretl_matrix_get(fsys->uhat, t, j);
		x += p * u;
	    }
	    gretl_vector_set(fsys->arty, k++, x);
	}
    }

#if FDEBUG > 1
    gretl_matrix_print(fsys->arty, "fiml artificial Y");
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
fiml_form_indepvars (fiml_system *fsys, const DATASET *dset, 
		     int t1)
{
    const int *enlist = system_get_endog_vars(fsys->sys);
    const int *exlist = system_get_instr_vars(fsys->sys);
    int i, j, k, t;
    int xrow, xcol = 0;
    double p, xjt;

    gretl_matrix_zero(fsys->artx);

    for (i=0; i<fsys->g; i++) { 
	/* loop across equations */
	const int *list = system_get_list(fsys->sys, i);

	for (j=2; j<=list[0]; j++) { 
	    /* loop across RHS vars */
	    const double *xj = NULL;
	    int vj = 0;

	    if (on_exo_list(exlist, list[j])) {
		/* the variable is exogenous or predetermined */
		xj = dset->Z[list[j]] + t1;
	    } else {
		/* RHS endogenous variable */
		vj = endo_var_number(enlist, list[j]);
	    }

	    for (t=0; t<fsys->n; t++) { 
		/* loop across obs */
		for (k=0; k<fsys->g; k++) { 
		    /* loop across vertical blocks */
		    xrow = k * fsys->n + t;
		    p = gretl_matrix_get(fsys->psi, k, i);
		    if (p != 0.0) {
			if (xj != NULL) {
			    xjt = xj[t];
			} else {
			    xjt = gretl_matrix_get(fsys->WB2, t, vj);
			}
			gretl_matrix_set(fsys->artx, xrow, xcol, xjt * p);
		    }
		}
	    }
	    xcol++;
	}
    }

#if FDEBUG > 1
    gretl_matrix_print(fsys->artx, "fiml artificial X");
#endif
}

#if FDEBUG

/* check: set initial residual matrix based on 3SLS */

static void fiml_uhat_init (fiml_system *fsys)
{
    gretl_matrix *uhat = fsys->sys->uhat;
    double x;
    int i, t;

    for (i=0; i<fsys->g; i++) {
	for (t=0; t<fsys->n; t++) {
	    x = gretl_matrix_get(uhat, t, i);
	    gretl_matrix_set(fsys->uhat, t, i, x);
	}
    }

    gretl_matrix_print(fsys->uhat, "uhat from 3SLS");
}

#endif

static int 
rhs_var_in_eqn (const int *list, int v)
{
    if (list != NULL) {
	int i;

	for (i=2; i<=list[0]; i++) {
	    if (list[i] == v) {
		return i;
	    }
	}
    }

    return 0;
}

/* initialize Gamma matrix based on 3SLS estimates plus identities */

static void 
fiml_G_init (fiml_system *fsys, const DATASET *dset)
{
    const int *enlist = system_get_endog_vars(fsys->sys);
    const int *slist; 
    const MODEL *pmod;
    int lv, rv;
    int i, j, vi;

    for (j=0; j<fsys->nendo; j++) {
	/* outer loop across columns (equations) */
	if (j < fsys->g) {
	    slist = system_get_list(fsys->sys, j);
	} else {
	    slist = NULL;
	}
	    
	lv = enlist[j + 1];

	/* inner loop across variables in equation */
	for (i=0; i<fsys->nendo; i++) {
	    rv = enlist[i + 1];

	    if (slist != NULL) {
		/* column pertains to stochastic equation */
		if (rv == slist[1]) {
		    gretl_matrix_set(fsys->G, i, j, 1.0);
		} else {
		    vi = rhs_var_in_eqn(slist, rv);
		    if (vi > 0) {
			pmod = system_get_model(fsys->sys, j);
			gretl_matrix_set(fsys->G, i, j, -pmod->coeff[vi-2]);
		    } else {
			gretl_matrix_set(fsys->G, i, j, 0.0);
		    }
		}
	    } else {
		/* column pertains to identity */
		if (lv == rv) {
		    vi = 1.0;
		} else {
		    vi = -1 * rhs_var_in_identity(fsys->sys, lv, rv);
		}
		gretl_matrix_set(fsys->G, i, j, vi);
	    }
	}
    }

#if FDEBUG
    printf("Order of columns (and rows):");
    for (i=1; i<=enlist[0]; i++) {
	printf(" %s", dset->varname[enlist[i]]);
    }
    putchar('\n');
    gretl_matrix_print(fsys->G, "fiml Gamma");
#endif
}

/* update Gamma matrix with revised parameter estimates */

static void fiml_G_update (fiml_system *fsys)
{
    const int *enlist = system_get_endog_vars(fsys->sys);   
    const int *slist;
    const MODEL *pmod;
    int i, j, rv, vi;

    for (j=0; j<fsys->g; j++) {
	slist = system_get_list(fsys->sys, j);
	for (i=0; i<fsys->nendo; i++) {
	    rv = enlist[i + 1];
	    if (rv != slist[1]) {
		vi = rhs_var_in_eqn(slist, rv);
		if (vi > 0) {
		    pmod = system_get_model(fsys->sys, j);
		    gretl_matrix_set(fsys->G, i, j, -pmod->coeff[vi-2]);
		} 
	    }
	}
    }

#if FDEBUG
    gretl_matrix_print(fsys->G, "fiml Gamma");
#endif
}

/* initialize B matrix based on 3SLS estimates and identities */

static void fiml_B_init (fiml_system *fsys, const DATASET *dset)
{
    const int *enlist = system_get_endog_vars(fsys->sys);
    const int *exlist = system_get_instr_vars(fsys->sys);
    const int *slist;
    const MODEL *pmod;
    int lv, rv;
    int i, j, vi;

    for (j=0; j<fsys->nendo; j++) {
	slist = system_get_list(fsys->sys, j);
	lv = enlist[j + 1];
	/* outer loop across columns (equations) */
	for (i=0; i<fsys->nexo; i++) {
	    rv = exlist[i + 1];
	    if (j < fsys->g) {
		/* column pertains to stochastic equation */
		vi = rhs_var_in_eqn(slist, rv);
		if (vi > 0) {
		    pmod = system_get_model(fsys->sys, j);
		    gretl_matrix_set(fsys->B, i, j, pmod->coeff[vi-2]);
		} else {
		    gretl_matrix_set(fsys->B, i, j, 0.0);
		}
	    } else {
		vi = rhs_var_in_identity(fsys->sys, lv, rv);
		gretl_matrix_set(fsys->B, i, j, vi);
	    }
	}
    }

#if FDEBUG
    printf("Order of columns:");
    for (i=1; i<=enlist[0]; i++) {
	printf(" %s", dset->varname[enlist[i]]);
    }
    putchar('\n');
    printf("Order of rows:");
    for (i=1; i<=exlist[0]; i++) {
	printf(" %s", dset->varname[exlist[i]]);
    }
    putchar('\n');
    gretl_matrix_print(fsys->B, "fiml B");
#endif
}

/* update B matrix with revised parameter estimates */

static void fiml_B_update (fiml_system *fsys)
{
    const int *exlist = system_get_instr_vars(fsys->sys);
    const int *slist;
    const MODEL *pmod;
    int i, j, vi;

    for (j=0; j<fsys->g; j++) {
	slist = system_get_list(fsys->sys, j);
	for (i=0; i<fsys->nexo; i++) {
	    vi = rhs_var_in_eqn(slist, exlist[i + 1]);
	    if (vi > 0) {
		pmod = system_get_model(fsys->sys, j);
		gretl_matrix_set(fsys->B, i, j, pmod->coeff[vi-2]);
	    } 
	}
    }

#if FDEBUG
    gretl_matrix_print(fsys->B, "fiml B");
#endif
}

/* calculate log-likelihood for FIML system */

static int fiml_ll (fiml_system *fsys, const DATASET *dset,
		    int t1)
{
    double tr;
    double ldetG;
    double ldetS;
    int i, j, t;
    int err = 0;

    fsys->ll = 0.0;

    /* form \hat{\Sigma} (ETM, equation 12.81); invert and
       Cholesky-decompose to get \Psi while we're at it 
    */
    err = fiml_form_sigma_and_psi(fsys, dset, t1);
    if (err) {
	fputs("fiml_form_sigma_and_psi: failed\n", stderr);
	return err;
    }

    /* note: make copy because this determinant calculation
       destroys the original matrix */
    gretl_matrix_copy_values(fsys->Gtmp, fsys->G);
    ldetG = gretl_matrix_log_abs_determinant(fsys->Gtmp, &err);
    if (err) {
	return err;
    }

    /* vcv_log_determinant() doesn't overwrite */
    ldetS = gretl_vcv_log_determinant(fsys->sigma, &err);
    if (err) {
	return err;
    }

    /* Davidson and MacKinnon, ETM, equation (12.80) */

    fsys->ll -= (fsys->gn / 2.0) * LN_2_PI;
    fsys->ll -= (fsys->n / 2.0) * ldetS;
    fsys->ll += fsys->n * ldetG;

    gretl_matrix_copy_values(fsys->Stmp, fsys->sigma);
    err = gretl_invert_symmetric_matrix(fsys->Stmp);   
    if (err) {
	return err;
    }

    tr = 0.0;
    for (i=0; i<fsys->g; i++) {
	double epe, eti, etj, sij;

	for (j=0; j<fsys->g; j++) {
	    epe = 0.0;
	    for (t=0; t<fsys->n; t++) {
		eti = gretl_matrix_get(fsys->uhat, t, i);
		etj = gretl_matrix_get(fsys->uhat, t, j);
		epe += eti * etj;
	    }
	    sij = gretl_matrix_get(fsys->Stmp, i, j);
	    tr += sij * epe;
	}
    }

    fsys->ll -= 0.5 * tr;

    return 0;
}

/* calculate instrumented version of endogenous variables, using
   the "restricted reduced form": WB\Gamma^{-1}.  Davidson and
   MacKinnon, ETM, equation (12.70)
*/

static int fiml_endog_rhs (fiml_system *fsys, const DATASET *dset,
			   int t1)
{
    int err;

    gretl_matrix_copy_values(fsys->Gtmp, fsys->G);
    err = gretl_invert_general_matrix(fsys->Gtmp);

    if (err) {
	fputs("inversion of G failed\n", stderr);
    } else {
#if FDEBUG
	gretl_matrix_print(fsys->Gtmp, "G-inverse");
#endif
	gretl_matrix_multiply(fsys->WB1, fsys->Gtmp, fsys->WB2);
    } 

    return err;
}

static void copy_estimates_to_btmp (fiml_system *fsys)
{
    const MODEL *pmod;
    int i, j, k = 0;

    for (i=0; i<fsys->g; i++) {
	pmod = system_get_model(fsys->sys, i);
	for (j=0; j<pmod->ncoeff; j++) {
	    gretl_vector_set(fsys->btmp, k++, pmod->coeff[j]);
	}
    }
}

/* adjust parameter estimates based on results of the artificial
   regression 
*/

static int
fiml_adjust_estimates (fiml_system *fsys, const DATASET *dset,
		       int t1, double *instep)
{
    MODEL *pmod;
    double llbak = fsys->ll;
    double minstep = 1.0e-06;
    double step = 4.0;
    int improved = 0;
    int err = 0;

    /* make a backup copy of the current parameter estimates */
    copy_estimates_to_btmp(fsys);

#if FDEBUG
    gretl_matrix_print(fsys->btmp, "parameter estimates");
    gretl_matrix_print(fsys->artb, "estimated gradients");
#endif

    while (!improved && !err && step > minstep) {
	double bk, delta;
	int i, j, k = 0;

	/* new coeff = old + gradient * step */
	for (i=0; i<fsys->g; i++) {
	    pmod = system_get_model(fsys->sys, i);
	    for (j=0; j<pmod->ncoeff; j++) {
		bk = gretl_vector_get(fsys->btmp, k);
		delta = gretl_vector_get(fsys->artb, k) * step;
		pmod->coeff[j] = bk + delta;
		k++;
	    }
	}

	/* write the new estimates into the G and B matrices */
	fiml_G_update(fsys);
	fiml_B_update(fsys);

	/* has the likelihood improved? */
	err = fiml_ll(fsys, dset, t1);
	if (!err) {
	    if (fsys->ll > llbak) {
		improved = 1;
	    } else {
		step /= 2.0;
	    } 
	}
    }

    *instep = step;

    return err;
}

/* get standard errors for FIML estimates from the covariance
   matrix of the artificial OLS regression
*/

static int 
fiml_get_std_errs (fiml_system *fsys, const gretl_matrix *R)
{
    gretl_matrix *vcv;
    int ldv = fsys->totk;
    int err;

    if (R != NULL) {
	ldv += R->rows;
    }

    vcv = gretl_matrix_alloc(ldv, ldv);
    if (vcv == NULL) {
	return E_ALLOC;
    }

    /* These are "GLS-type" standard errors: see Calzolari
       and Panattoni */
    
    if (R != NULL) {
	err = gretl_matrix_restricted_ols(fsys->arty, fsys->artx, R, NULL,
					  fsys->artb, vcv, NULL, NULL);
    } else {
	err = gretl_matrix_SVD_ols(fsys->arty, fsys->artx, fsys->artb, 
				   vcv, NULL, NULL);
    }

    if (!err) {
	MODEL *pmod;
	int i, j, k = 0;

	for (i=0; i<fsys->g; i++) {
	    pmod = system_get_model(fsys->sys, i);
	    for (j=0; j<pmod->ncoeff; j++) {
		pmod->sderr[j] = sqrt(gretl_matrix_get(vcv, k, k));
		k++;
	    }
	}
    }

    if (!err) {
	gretl_matrix_replace(&fsys->sys->vcv, vcv);
    } else {
	gretl_matrix_free(vcv);
    }

    return err;
}

static void fiml_print_gradients (const gretl_matrix *b, PRN *prn)
{
    int i;

    pprintf(prn, "\n%s:\n\n", _("Gradients at last iteration"));

    for (i=0; i<b->rows; i++) {
	pprintf(prn, " %14e ", b->val[i]);
	if ((i + 1) % 4 == 0) {
	    pputc(prn, '\n');
	}
    }
    pputc(prn, '\n');
}

/* Driver function for FIML as described in Davidson and MacKinnon,
   ETM, chap 12, section 5.
*/

#define FIML_ITER_MAX 250

int fiml_driver (equation_system *sys, DATASET *dset, 
		 gretlopt opt, PRN *prn)
{
    const gretl_matrix *R = NULL;
    fiml_system *fsys;
    int t1 = dset->t1;
    double llbak;
    double crit = 1.0;
    double tol = 1.0e-12; /* over-ambitious? */
    double bigtol = 1.0e-9;
    int verbose = 0;
    int iters = 0;
    int err = 0;

    fsys = fiml_system_new(sys, &err);
    if (err) {
	return err;
    }

    if ((opt & OPT_V) && !(opt & OPT_Q)) {
	verbose = 1;
    }

#if FDEBUG
    /* check uhat calculation: set intial uhat based on 3SLS */
    fiml_uhat_init(fsys);
#endif

    /* intialize Gamma coefficient matrix */
    fiml_G_init(fsys, dset);

    /* intialize B coefficient matrix */
    fiml_B_init(fsys, dset);

    /* initial loglikelihood */
    err = fiml_ll(fsys, dset, t1);
    if (err) {
	fputs("fiml_ll: failed\n", stderr);
	goto bailout;
    } else {
	llbak = fsys->ll;
	if (verbose) {
	    pprintf(prn, "*** initial ll = %.8g\n", fsys->ll);
	}
    } 

    if ((sys->flags & SYSTEM_RESTRICT) && sys->R != NULL) {
	R = sys->R;
    }

    while (crit > tol && iters < FIML_ITER_MAX) {
	double step;

	/* form LHS vector for artificial regression */
	fiml_form_depvar(fsys);

	/* instrument the RHS endog vars */
	err = fiml_endog_rhs(fsys, dset, t1);
	if (err) {
	    fputs("fiml_endog_rhs: failed\n", stderr);
	    break;
	}	

	/* form RHS matrix for artificial regression */
	fiml_form_indepvars(fsys, dset, t1);

	/* run artificial regression (ETM, equation 12.86) */
	if (R != NULL) {
	    err = gretl_matrix_restricted_ols(fsys->arty, fsys->artx, 
					      R, NULL, fsys->artb, 
					      NULL, NULL, NULL);
	} else {
#if 0
	    err = gretl_matrix_svd_ols(fsys->arty, fsys->artx, 
				       fsys->artb, NULL, NULL, NULL);
#else
	    err = gretl_matrix_ols(fsys->arty, fsys->artx, 
				   fsys->artb, NULL, NULL, NULL);
#endif
	}

	if (err) {
	    fputs("gretl_matrix_ols: failed\n", stderr);
	    break;
	}

	/* adjust param estimates based on gradients in fsys->artb */
	err = fiml_adjust_estimates(fsys, dset, t1, &step);
	if (err) {
	    break;
	}

	if (verbose) {
	    pprintf(prn, "*** iteration %3d: step = %g, ll = %.8g\n", 
		    iters + 1, step, fsys->ll);
	}

	crit = fsys->ll - llbak;
	llbak = fsys->ll;

	iters++;
    }

    if (verbose) {
	if (crit < tol) {
	    pprintf(prn, "\nTolerance %g, criterion %g\n", tol, crit);
	} else if (crit < bigtol) {
	    pprintf(prn, "\nTolerance %g, criterion %g\n", bigtol, crit);
	} else {
	    pputc(prn, '\n');
	    pprintf(prn, "Tolerance of %g was not met\n", bigtol);
	    err = 1;
	}
    }

    if (!err) {
	if (verbose) {
	    fiml_print_gradients(fsys->artb, prn);
	}
	err = fiml_get_std_errs(fsys, R);
    }

    if (R != NULL && verbose) {
	fiml_overid_test(fsys, dset);
    }

    /* write the results into the parent system */
    fiml_transcribe_results(fsys, dset, t1, iters);

 bailout:
    
    fiml_system_destroy(fsys);

    return err;
}



    
	

    
    
    

    
    
