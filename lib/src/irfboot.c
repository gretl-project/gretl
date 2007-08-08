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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

/* bootstrapped confidence intervals for impulse response functions */

#include "libgretl.h" 
#include "var.h"  
#include "vartest.h"
#include "matrix_extra.h"

#define BDEBUG 0

#if BDEBUG
# define BOOT_ITERS 5
#else
# define BOOT_ITERS 999
#endif

typedef struct irfboot_ irfboot;

struct irfboot_ {
    gretlopt opt;       /* VECM options flags */
    int ncoeff;         /* number of coefficients per equation */
    int horizon;        /* horizon for impulse responses */
    gretl_matrix *A;    /* augmented coefficient matrix */
    gretl_matrix *C;    /* error covariance matrix */
    gretl_matrix *S;    /* covariance matrix of residuals */
    gretl_matrix *rE;   /* matrix of resampled original residuals */
    gretl_matrix *Xt;   /* row t of X matrix */
    gretl_matrix *Yt;   /* yhat at t */
    gretl_matrix *Et;   /* residuals at t */
    gretl_matrix *rtmp; /* temporary storage */
    gretl_matrix *ctmp; /* temporary storage */
    gretl_matrix *resp; /* impulse response matrix */
    gretl_matrix *C0;   /* initial coefficient estimates (VECM only) */
    gretl_matrix *Y;    /* levels of dependent vars (VECM only) */
    int *sample;        /* resampling array */
};

static void irf_boot_free (irfboot *b)
{
    if (b == NULL) {
	return;
    }

    gretl_matrix_free(b->A);
    gretl_matrix_free(b->C);
    gretl_matrix_free(b->S);
    gretl_matrix_free(b->rE);
    gretl_matrix_free(b->Xt);
    gretl_matrix_free(b->Yt);
    gretl_matrix_free(b->Et);
    gretl_matrix_free(b->rtmp);
    gretl_matrix_free(b->ctmp);
    gretl_matrix_free(b->resp);
    gretl_matrix_free(b->C0);
    gretl_matrix_free(b->Y);

    free(b->sample);
    free(b);
}

/* Is this redundant? */

static int boot_A_C_init (irfboot *b, const GRETL_VAR *v)
{
    int n = v->neqns * v->order;
    int i, j, err = 0;
    
    b->A = gretl_matrix_alloc(n, n);

    if (b->A == NULL) {
	err = E_ALLOC;
    } else {
	for (i=v->neqns; i<n; i++) {
	    for (j=0; j<n; j++) {
		gretl_matrix_set(b->A, i, j, (j == i - v->neqns)? 1.0 : 0.0);
	    }
	}
    }

    if (!err) {
	b->C = gretl_zero_matrix_new(n, v->neqns);
	if (b->C == NULL) {
	    err = E_ALLOC;
	    gretl_matrix_free(b->A);
	    b->A = NULL;
	} 
    }

    return err;
}

static int boot_tmp_init (irfboot *b, const GRETL_VAR *v)
{
    int n = v->neqns * (v->order + v->ecm);

    b->rtmp = gretl_matrix_alloc(n, v->neqns);
    b->ctmp = gretl_matrix_alloc(n, v->neqns);

    if (b->rtmp == NULL || b->ctmp == NULL) {
	return E_ALLOC;
    }

    return 0;
}

static gretlopt boot_get_opt (const GRETL_VAR *v)
{
    gretlopt opt = OPT_NONE;

    if (v->jinfo != NULL) {
	if (v->jinfo->code == J_NO_CONST) {
	    opt = OPT_N;
	} else if (v->jinfo->code == J_REST_CONST) {
	    opt = OPT_R;
	} else if (v->jinfo->code == J_REST_TREND) {
	    opt = OPT_A;
	} else if (v->jinfo->code == J_UNREST_TREND) {
	    opt = OPT_T;
	}

	if (v->jinfo->seasonals) {
	    opt |= OPT_D;
	}

	/* OPT_S: save (don't delete) the variables added
	   to the dataset in the Johansen process.
	   OPT_B: flags the fact that we're doing impulse
	   response bootstrap.
	*/
	opt |= (OPT_S | OPT_B);
    }

    return opt;
}

static irfboot *irf_boot_new (const GRETL_VAR *var, int periods)
{
    irfboot *b;
    int err = 0;

    b = malloc(sizeof *b);
    if (b == NULL) {
	return NULL;
    }

    b->C = NULL;
    b->A = NULL;
    b->S = NULL;
    b->rE = NULL;
    b->Xt = NULL;
    b->Yt = NULL;
    b->Et = NULL;
    b->rtmp = NULL;
    b->ctmp = NULL;
    b->resp = NULL;
    b->C0 = NULL;

    b->sample = NULL;
    b->horizon = periods;

    if (var->jinfo != NULL) {
	b->ncoeff = var->ncoeff + (var->neqns - jrank(var)) + 
	    restricted(var); /* FIXME? */
    } else {
	b->ncoeff = var->ncoeff;
    }

    b->opt = boot_get_opt(var);

#if BDEBUG
    fprintf(stderr, "boot: t1=%d, t2=%d, ncoeff=%d, horizon=%d\n",
	    var->t1, var->t2, b->ncoeff, b->horizon);
    fprintf(stderr, " T=%d, neqns=%d, order=%d, ecm=%d, ifc=%d\n",
	    var->T, var->neqns, var->order, var->ecm, var->ifc);
#endif

    err = boot_tmp_init(b, var);

    if (!err && !var->ecm) {
	/* we don't need these matrices for VECMs */
	err = boot_A_C_init(b, var);
	if (!err) {
	    b->S = gretl_matrix_alloc(var->neqns, var->neqns);
	    if (b->S == NULL) {
		err = E_ALLOC;
	    }
	}	    
    }

    if (!err) {
	b->rE = gretl_matrix_alloc(var->T, var->neqns);
	b->resp = gretl_matrix_alloc(b->horizon, BOOT_ITERS);
	b->Xt = gretl_matrix_alloc(1, var->X->cols);
	b->Yt = gretl_matrix_alloc(1, var->neqns);
	b->Et = gretl_matrix_alloc(1, var->neqns);
	b->sample = malloc(var->T * sizeof *b->sample);

	if (b->rE == NULL || b->resp == NULL || b->Xt == NULL ||
	    b->Yt == NULL || b->Et == NULL || b->sample == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	irf_boot_free(b);
	b = NULL;
    }

    return b;
}

static void
recalculate_impulse_responses (irfboot *b, const gretl_matrix *A, 
			       const gretl_matrix *C,
			       int targ, int shock, int iter)
{
    double rt;
    int t;

    for (t=0; t<b->horizon; t++) {
	if (t == 0) {
	    /* initial estimated responses */
	    gretl_matrix_copy_values(b->rtmp, C);
	} else {
	    /* calculate further estimated responses */
	    gretl_matrix_multiply(A, b->rtmp, b->ctmp);
	    gretl_matrix_copy_values(b->rtmp, b->ctmp);
	}

	rt = gretl_matrix_get(b->rtmp, targ, shock);
	gretl_matrix_set(b->resp, t, iter, rt);
    }
}

/* transcribe coefficients from re-estimated VAR model to the 
   boot->A matrix */

static int write_boot_A_matrix (irfboot *b, GRETL_VAR *var)
{
    int i, j, v, lag;
    int dim = var->neqns * var->order;
    double bij;

    for (j=0; j<var->neqns; j++) {
	v = lag = 0;
	for (i=0; i<dim; i++) {
	    bij = gretl_matrix_get(var->B, i+var->ifc, j);
	    gretl_matrix_set(b->A, j, var->neqns * lag + v, bij);
	    if (lag < var->order - 1) {
		lag++;
	    } else {
		lag = 0;
		v++;
	    }
	}
    }

    return 0;
}

/* In re-estimation of VAR or VECM we'll tolerate at most 4 cases of
   near-perfect collinearity (which can arise by chance): maybe this
   should be more flexible? 
*/

#define MAXSING 5
#define VAR_FATAL(e,i,s) (e && (e != E_SINGULAR || i == 0 || s >= MAXSING))

#if 0

/* FIXME !! */

static int 
re_estimate_VECM (irfboot *b, GRETL_VAR *jvar, int targ, int shock, 
		  int iter, int scount)
{
    static int (*jbr) (GRETL_VAR *, double **, DATAINFO *, int) = NULL;
    static void *handle = NULL;
    int err = 0;

    gretl_error_clear();

    if (iter == 0) {
	/* first round: open the Johansen plugin */
	jbr = get_plugin_function("johansen_boot_round", &handle);
	if (jbr == NULL) {
	    err = 1;
	}
    }

    if (!err) {
	err = johansen_stage_1(jvar, Z, pdinfo, b->opt, NULL);
    }

    if (!err) {
	/* call the plugin function */
	err = jbr(jvar, (const double **) b->Z, boot->dinfo, iter);
    }

    if (!err) {   
#if BDEBUG
	gretl_matrix_print(jvar->S, "jvar->S (Omega)");
#endif
	err = gretl_VAR_do_error_decomp(jvar->S, jvar->C);
    }

    if (!err) {
	recalculate_impulse_responses(b, jvar->A, jvar->C, targ, shock, iter);
    }

    if (iter == BOOT_ITERS - 1 || VAR_FATAL(err, iter, scount)) {
	if (handle != NULL) {
	    close_plugin(handle);
	    handle = NULL;
	    jbr = NULL;
	}
    }

    return err;
}

#endif

static int re_estimate_VAR (irfboot *b, GRETL_VAR *v, int targ, int shock, 
			    int iter)
{
    int err;

    if (v->qr) {
	err = gretl_matrix_QR_ols(v->Y, v->X, v->B, v->E,
				  NULL, NULL);
    } else {
	err = gretl_matrix_multi_ols(v->Y, v->X, v->B, v->E,
				     NULL);
    }

    if (!err) {
	write_boot_A_matrix(b, v);
    }
    
    if (!err) {    
	gretl_matrix_multiply_mod(v->E, GRETL_MOD_TRANSPOSE,
				  v->E, GRETL_MOD_NONE,
				  b->S, GRETL_MOD_NONE);
	gretl_matrix_divide_by_scalar(b->S, v->T);
	err = gretl_VAR_do_error_decomp(b->S, b->C);
    }

    if (!err) {
	recalculate_impulse_responses(b, b->A, b->C, targ, shock, iter);
    }

    return err;
}

/* VECM: make a vector, rbeta, containing the coefficients on the
   "restricted" constant or trend, if this is required.
*/

static gretl_matrix *make_restricted_coeff_vector (const GRETL_VAR *v)
{
    gretl_matrix *b = NULL;
    gretl_matrix *rbeta = NULL;
    int rank = v->jinfo->rank;
    double x;
    int j, err = 0;

    b = gretl_column_vector_alloc(rank);
    if (b == NULL) {
	return NULL;
    }

    for (j=0; j<rank; j++) {
	/* extract the last row of \beta, transposed */
	x = gretl_matrix_get(v->jinfo->Beta, v->neqns, j);
	gretl_vector_set(b, j, x);
    }

    rbeta = gretl_matrix_multiply_new(v->jinfo->Alpha, b, &err);

#if BDEBUG > 1
    gretl_matrix_print(b, "restricted var row of beta'");
    gretl_matrix_print(rbeta, "coeffs for restricted term");
#endif

    gretl_matrix_free(b);

    return rbeta;
}

/* VECM: make a record of all the original coefficients in the VAR
   representation of the VECM, so we have these on hand in convenient
   form for recomputing the dataset after resampling the residuals,
   or for forecasting.  

   We get these coefficients from three sources: (1) the 
   companion matrix, A; (2) the VECM models (unrestricted constant
   and/or trend, seasonals if applicable); and (3) the implied
   coefficient on a restricted constant or trend.

   Note that the construction here depends on the order in which
   variables are added to the dataset; that order can't be changed 
   without breaking stuff here.
*/

gretl_matrix *VAR_coeff_matrix_from_VECM (const GRETL_VAR *var)
{
    gretl_matrix *C0 = NULL;
    gretl_matrix *rbeta = NULL;
    int order = var->order + 1;
    int nexo = (var->xlist != NULL)? var->xlist[0] : 0;
    int ndelta = var->order * var->neqns;
    int nseas = var->jinfo->seasonals;
    int ncoeff;
    int X0, S0, T0;
    double aij;
    int i, j, k;

    /* total coeffs in VAR representation */
    ncoeff = var->ncoeff + (var->neqns - var->jinfo->rank) + 
	restricted(var);

    if (restricted(var)) {
	rbeta = make_restricted_coeff_vector(var);
	if (rbeta == NULL) {
	    return NULL;
	}
    }

    C0 = gretl_matrix_alloc(var->neqns, ncoeff);
    if (C0 == NULL) {
	gretl_matrix_free(rbeta);
	return NULL;
    }

    /* position of first exog var coeff in VECM model */
    X0 = var->ifc + ndelta;
    /* position of first seasonal coeff in VECM model */
    S0 = X0 + nexo;
    /* position of trend coeff in VECM model */
    T0 = S0 + nseas;

    for (i=0; i<var->neqns; i++) {
	const MODEL *pmod = var->models[i];
	int col = 0;

	/* constant, if present */
	if (var->ifc) {
	    gretl_matrix_set(C0, i, col++, pmod->coeff[0]);
	}

	/* endogenous vars: use companion matrix */
	for (j=0; j<var->neqns; j++) {
	    for (k=0; k<order; k++) {
		aij = gretl_matrix_get(var->A, i, k * var->neqns + j);
		gretl_matrix_set(C0, i, col++, aij);
	    }
	}

	/* exogenous vars */
	for (j=0; j<nexo; j++) {
	    gretl_matrix_set(C0, i, col++, pmod->coeff[X0+j]);
	}

	/* seasonals, if present */
	for (j=0; j<nseas; j++) {
	    gretl_matrix_set(C0, i, col++, pmod->coeff[S0+j]);
	}

	/* unrestricted trend, if present */
	if (jcode(var) == J_UNREST_TREND) {
	    gretl_matrix_set(C0, i, col++, pmod->coeff[T0]);
	}

	/* restricted term (const or trend), if present */
	if (rbeta != NULL) {
	    aij = gretl_vector_get(rbeta, i);
	    gretl_matrix_set(C0, i, col, aij);
	} 	
    }

    if (rbeta != NULL) {
	gretl_matrix_free(rbeta);
    }

#if BDEBUG > 1
    gretl_matrix_print(var->A, "var->A");
    gretl_matrix_print(C0, "C0");
#endif

    return C0;
}

#if 0

/* VECM: copy levels of Y vars from external dataset into
   b->Y matrix
*/

static int init_VECM_dataset (irfboot *b, const GRETL_VAR *var,
			      const double **Z)
{
    int i, vi, s, t;

    b->Y = gretl_matrix_alloc(var->T, var->neqns);
    if (b->Y == NULL) {
	return E_ALLOC;
    }
    
    for (i=0; i<var->neqns; i++) {
	vi = var->ylist[i+1];
	s = 0;
	for (t=var->t1; t<=var->t2; t++) {
	    gretl_matrix_set(b->Y, s++, i, Z[vi][t];
	}
    }

    return 0;
}

/* broken for now (and in fact, was broken before now?).  
   has to be reworked.  notice we're using the VAR
   representation of the vecm here!
*/

static void 
compute_VECM_dataset (irfboot *boot, const GRETL_VAR *var, int iter)
{
    int order = var->order + 1;
    int nexo = (var->xlist != NULL)? var->xlist[0] : 0;
    int nseas = var->jinfo->seasonals;
    int i, j, k, vj, t;

#if BDEBUG
    fprintf(stderr, "compute_VECM_dataset: order=%d, nexo=%d, nseas=%d, t1=%d\n",
	    order, nexo, nseas, var->t1);
#endif

    for (t=0; t<var->T; t++) {
	for (i=0; i<var->neqns; i++) {
	    double xti, cij, eti, bti = 0.0;
	    int xcol = var->ifc + var->neqns * var->order;
	    int col = 0;

	    /* unrestricted constant, if present */
	    if (var->ifc) {
		cij = gretl_matrix_get(boot->C0, i, col++);
		bti += cij;
	    }

	    /* lags of endogenous vars */
	    for (j=0; j<boot->neqns; j++) {
		vj = var->ylist[j+1];
		for (k=1; k<=order; k++) {
		    cij = gretl_matrix_get(boot->C0, i, col++);
		    bti += cij * Z[vj][t-k]; /* wrong: use b->Y */
		}
	    }

	    /* exogenous vars, if present */
	    for (j=0; j<nexo; j++) {
		cij = gretl_matrix_get(boot->C0, i, col++);
		xti = gretl_matrix_get(var->X, t, xcol++);
		bti += cij * xti;
	    }

	    /* seasonals, if present */
	    for (j=0; j<nseas; j++) {
		cij = gretl_matrix_get(boot->C0, i, col++);
		xti = gretl_matrix_get(var->X, t, xcol++);
		bti += cij * xti;
	    }

	    if (jcode(var) == J_UNREST_TREND) {
		/* unrestricted trend */
		cij = gretl_matrix_get(boot->C0, i, col);
		bti += cij * (t + var->t1 + 1);
	    } else if (jcode(var) == J_REST_CONST) {
		/* restricted constant */
		bti += gretl_matrix_get(boot->C0, i, col);
	    } else if (jcode(var) == J_REST_TREND) {
		/* restricted trend */
		cij = gretl_matrix_get(boot->C0, i, col);
		bti += cij * (t + var->t1);
	    }

	    /* set value of dependent var to fitted + re-sampled error */
	    eti = gretl_matrix_get(boot->rE, t - boot->t1, i);
	    gretl_matrix_set(b->Y, t, i, bti + ati);
	}
    }

    if (iter > 0) {
	int vl = 1 + var->neqns + nexo + nseas;
	int vd = vl + var->neqns;

	/* recompute first lags and first differences */
	for (i=0; i<boot->neqns; i++) {
	    for (t=1; t<boot->dinfo->n; t++) {
		boot->Z[vl][t] = boot->Z[i+1][t-1];
		boot->Z[vd][t] = boot->Z[i+1][t] - boot->Z[i+1][t-1];
	    }
	    vl++;
	    vd++;
	}
    }

#if BDEBUG > 1
    gretl_matrix_print(b->Y, "b->Y (vecm, levels)");
#endif
}

#endif

/* (Re-)fill the bootstrap dataset with artificial data, based on the
   re-sampled residuals from the original VAR (simple VAR, not VECM).
*/

static void compute_VAR_dataset (irfboot *b, GRETL_VAR *v,
				 const GRETL_VAR *vbak)
{
    double x;
    int i, j, k, t;

    for (t=0; t<v->T; t++) {
	/* extract "X" at t */
	gretl_matrix_extract_matrix(b->Xt, v->X, t, 0, GRETL_MOD_NONE);

	/* multiply into original coeff matrix */
	gretl_matrix_multiply(b->Xt, vbak->B, b->Yt);

	/* get resampled residuals at t */
	gretl_matrix_extract_matrix(b->Et, b->rE, t, 0, GRETL_MOD_NONE);

	/* add resampled residual */
	gretl_matrix_add_to(b->Yt, b->Et);

	/* write into Y */
	gretl_matrix_inscribe_matrix(v->Y, b->Yt, t, 0, GRETL_MOD_NONE);

	/* revise lagged Y columns in X */
	k = v->ifc;
	for (i=0; i<v->neqns; i++) {
	    x = b->Yt->val[i];
	    for (j=1; j<=v->order && t+j < v->T; j++) {
		gretl_matrix_set(v->X, t+j, k++, x);
	    }
	}
    }

#if BDEBUG > 1
    gretl_matrix_print(v->Y, "var->Y after resampling");
    gretl_matrix_print(v->X, "var->X after resampling");
#endif
}

/* Resample the original VAR or VECM residuals.  Note the option to
   "resample" _without_ actually changing the order, if BDEBUG > 1.
   This is useful for checking that the IRF bootstrap rounds are
   idempotent: we should then get exactly the same set of responses as
   in the original estimation of the VAR/VECM.
*/

static void resample_resids (irfboot *b, const GRETL_VAR *vbak)
{
    double eti;
    int i, t;

    /* construct sampling array */
    for (t=0; t<vbak->T; t++) {
#if BDEBUG > 1
	b->sample[t] = t;
#else
	b->sample[t] = gretl_rand_int_max(vbak->T);
#endif
#if BDEBUG
	fprintf(stderr, "boot->sample[%d] = %d\n", t, b->sample[t]);
#endif
    }

    /* draw sample from the original residuals */
    for (t=0; t<vbak->T; t++) {
	for (i=0; i<vbak->neqns; i++) {
	    eti = gretl_matrix_get(vbak->E, b->sample[t], i);
	    gretl_matrix_set(b->rE, t, i, eti);
	}
    }
}

static int irf_boot_quantiles (irfboot *b, gretl_matrix *R)
{
    double alpha = 0.05;
    double *rk;
    int k, ilo, ihi;

    rk = malloc(BOOT_ITERS * sizeof *rk);
    if (rk == NULL) {
	return E_ALLOC;
    }

#if BDEBUG
    ilo = 1;
    ihi = BOOT_ITERS;
    fprintf(stderr, "IRF bootstrap (%d iters), min and max values\n", BOOT_ITERS);
#else
    ilo = (BOOT_ITERS + 1) * alpha / 2.0;
    ihi = (BOOT_ITERS + 1) * (1.0 - alpha / 2.0);
#endif

    for (k=0; k<b->horizon; k++) {
	gretl_matrix_row_to_array(b->resp, k, rk);
	qsort(rk, BOOT_ITERS, sizeof *rk, gretl_compare_doubles);
	gretl_matrix_set(R, k, 1, rk[ilo-1]);
	gretl_matrix_set(R, k, 2, rk[ihi-1]);
    }

    free(rk);

    return 0;
}

/* make a back-up copy of the VAR data that we'll want to restore on
   exit from the IRF bootstrap operation */

static GRETL_VAR *back_up_VAR (const GRETL_VAR *v)
{
    GRETL_VAR *vbak;
    int err = 0;

    vbak = malloc(sizeof *vbak);
    if (vbak == NULL) {
	return NULL;
    }

    gretl_VAR_clear(vbak);

    clear_gretl_matrix_err();

    vbak->Y = gretl_matrix_copy(v->Y);
    vbak->X = gretl_matrix_copy(v->X);
    vbak->B = gretl_matrix_copy(v->B);
    vbak->E = gretl_matrix_copy(v->E);
    vbak->A = gretl_matrix_copy(v->A);
    vbak->C = gretl_matrix_copy(v->C);
    vbak->S = gretl_matrix_copy(v->S);

    err = get_gretl_matrix_err();

    if (!err && v->jinfo != NULL) {
	vbak->jinfo = malloc(sizeof *vbak->jinfo);
	if (vbak->jinfo == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && vbak->jinfo != NULL) {
	vbak->jinfo->R0 = gretl_matrix_copy(v->jinfo->R0);
	vbak->jinfo->R1 = gretl_matrix_copy(v->jinfo->R1);
	vbak->jinfo->S00 = gretl_matrix_copy(v->jinfo->S00);
	vbak->jinfo->S11 = gretl_matrix_copy(v->jinfo->S11);
	vbak->jinfo->S01 = gretl_matrix_copy(v->jinfo->S01);
	vbak->jinfo->Beta = gretl_matrix_copy(v->jinfo->Beta);
	vbak->jinfo->Alpha = gretl_matrix_copy(v->jinfo->Alpha);
	err = get_gretl_matrix_err();
    }

    if (err) {
	gretl_VAR_free(vbak);
	vbak = NULL;
    } else {
	vbak->T = v->T;
	vbak->neqns = v->neqns;
	vbak->order = v->order;
    }

    return vbak;
}

static void restore_VAR_data (GRETL_VAR *v, GRETL_VAR *vbak)
{
    gretl_matrix_free(v->Y);
    gretl_matrix_free(v->X);
    gretl_matrix_free(v->B);
    gretl_matrix_free(v->E);
    gretl_matrix_free(v->A);
    gretl_matrix_free(v->C);
    gretl_matrix_free(v->S);

    v->Y = vbak->Y;
    v->X = vbak->X;
    v->B = vbak->B;
    v->E = vbak->E;
    v->A = vbak->A;
    v->C = vbak->C;
    v->S = vbak->S;

    if (v->jinfo != NULL) {
	gretl_matrix_free(v->jinfo->R0);
	gretl_matrix_free(v->jinfo->R1);
	gretl_matrix_free(v->jinfo->S00);
	gretl_matrix_free(v->jinfo->S11);
	gretl_matrix_free(v->jinfo->S01);
	gretl_matrix_free(v->jinfo->Beta);
	gretl_matrix_free(v->jinfo->Alpha);

	v->jinfo->R0 = vbak->jinfo->R0;
	v->jinfo->R1 = vbak->jinfo->R1;
	v->jinfo->S00 = vbak->jinfo->S00;
	v->jinfo->S11 = vbak->jinfo->S11;
	v->jinfo->S01 = vbak->jinfo->S01;
	v->jinfo->Beta = vbak->jinfo->Beta;
	v->jinfo->Alpha = vbak->jinfo->Alpha;

	free(vbak->jinfo);
    }

    free(vbak);
}

/* "locally public" function, called from var.c */

gretl_matrix *irf_bootstrap (GRETL_VAR *var, 
			     int targ, int shock, int periods,
			     const double **Z, 
			     const DATAINFO *pdinfo)
{
    gretl_matrix *R = NULL; /* the return value */
    GRETL_VAR *vbak = NULL;
    irfboot *boot = NULL;
    int scount = 0;
    int iter, err = 0;

#if BDEBUG
    fprintf(stderr, "\n*** irf_bootstrap() called\n");
#endif

    R = gretl_zero_matrix_new(periods, 3);
    if (R == NULL) {
	return NULL;
    }    

    vbak = back_up_VAR(var);
    if (vbak == NULL) {
	gretl_matrix_free(R);
	return NULL;
    }

    boot = irf_boot_new(var, periods);
    if (boot == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (var->ecm) {
	boot->C0 = VAR_coeff_matrix_from_VECM(var);
	if (boot->C0 == NULL) {
	    err = E_ALLOC;
	}
    }

    for (iter=0; iter<BOOT_ITERS && !err; iter++) {
#if BDEBUG
	fprintf(stderr, "starting iteration %d\n", iter);
#endif
	resample_resids(boot, vbak);
	if (var->ecm) {
#if 0
	    /* broken :-( */
	    compute_VECM_dataset(boot, var, iter);
	    err = re_estimate_VECM(boot, var, targ, shock, iter, scount);
#endif
	} else {
	    compute_VAR_dataset(boot, var, vbak);
	    err = re_estimate_VAR(boot, var, targ, shock, iter);
	}
	if (err && !VAR_FATAL(err, iter, scount)) {
	    /* excessive collinearity: try again, unless this is becoming a habit */
	    scount++;
	    iter--;
	    err = 0;
	}
    }

    if (err && scount == MAXSING) {
	strcpy(gretl_errmsg, "Excessive collinearity in resampled datasets");
    }

    if (!err) {
	err = irf_boot_quantiles(boot, R);
    }

    irf_boot_free(boot);

 bailout:

    restore_VAR_data(var, vbak);

    if (err) {
	gretl_matrix_free(R);
	R = NULL;
    }

    return R;
}


