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

/* bootstrapped confidence intervals for impulse response functions */

#include "libgretl.h"
#include "var.h"
#include "johansen.h"
#include "vartest.h"
#include "matrix_extra.h"
#include "libset.h"

#define BDEBUG 0

#if BDEBUG
# define DBG_BOOT_ITERS 5
#endif

typedef struct irfboot_ irfboot;

struct irfboot_ {
    int ncoeff;         /* number of coefficients per equation */
    int horizon;        /* horizon for impulse responses */
    int iters;          /* number of iterations */
    gretl_matrix_block *MB; /* wrapper for some of the following */
    gretl_matrix *rE;   /* matrix of resampled original residuals */
    gretl_matrix *Xt;   /* row t of X matrix */
    gretl_matrix *Yt;   /* yhat at t */
    gretl_matrix *Et;   /* residuals at t */
    gretl_matrix *rtmp; /* temporary storage */
    gretl_matrix *ctmp; /* temporary storage */
    gretl_matrix *resp; /* impulse response matrix */
    gretl_matrix *C0;   /* initial coefficient estimates (VECM only) */
    int *sample;        /* resampling array */
    DATASET *dset;      /* dummy dataset for levels (VECM only) */
};

static void irf_boot_free (irfboot *b)
{
    if (b == NULL) {
	return;
    }

    gretl_matrix_block_destroy(b->MB);

    gretl_matrix_free(b->Xt);
    gretl_matrix_free(b->Yt);
    gretl_matrix_free(b->Et);
    gretl_matrix_free(b->C0);

    if (b->dset != NULL) {
	destroy_dataset(b->dset);
    }

    free(b->sample);
    free(b);
}

static int boot_allocate (irfboot *b, const GRETL_VAR *v)
{
    int n = v->neqns * effective_order(v);

    b->MB = gretl_matrix_block_new(&b->rtmp, n, v->neqns,
				   &b->ctmp, n, v->neqns,
				   &b->rE, v->T, v->neqns,
				   &b->resp, b->horizon, b->iters,
				   NULL);
    if (b->MB == NULL) {
	return E_ALLOC;
    }

    b->sample = malloc(v->T * sizeof *b->sample);
    if (b->sample == NULL) {
	return E_ALLOC;
    }

    if (v->ci == VAR) {
	/* not needed for VECM */
	b->Xt = gretl_matrix_alloc(1, v->X->cols);
	b->Yt = gretl_matrix_alloc(1, v->neqns);
	b->Et = gretl_matrix_alloc(1, v->neqns);

	if (b->Xt == NULL || b->Yt == NULL || b->Et == NULL) {
	    return E_ALLOC;
	}
    }

    return 0;
}

static irfboot *irf_boot_new (const GRETL_VAR *var, int periods)
{
    irfboot *b;
    int err = 0;

    b = malloc(sizeof *b);
    if (b == NULL) {
	return NULL;
    }

    b->MB = NULL;
    b->Xt = NULL;
    b->Yt = NULL;
    b->Et = NULL;
    b->C0 = NULL;
    b->sample = NULL;
    b->dset = NULL;

    b->horizon = periods;
#if BDEBUG
    b->iters = DBG_BOOT_ITERS;
#else
    b->iters = libset_get_int(BOOT_ITERS);
#endif

    if (var->jinfo != NULL) {
	b->ncoeff = var->ncoeff + (var->neqns - jrank(var)) +
	    n_restricted_terms(var); /* FIXME? */
    } else {
	b->ncoeff = var->ncoeff;
    }

#if BDEBUG
    fprintf(stderr, "boot: t1=%d, t2=%d, ncoeff=%d, horizon=%d\n",
	    var->t1, var->t2, b->ncoeff, b->horizon);
    fprintf(stderr, " T=%d, neqns=%d, order=%d, rank=%d, ifc=%d\n",
	    var->T, var->neqns, var->order, jrank(var), var->ifc);
#endif

    err = boot_allocate(b, var);

    if (err) {
	irf_boot_free(b);
	b = NULL;
    }

    return b;
}

static int
recalculate_impulse_responses (irfboot *b, GRETL_VAR *var,
			       int targ, int shock, int iter)
{
    gretl_matrix *C = var->C;
    double x;
    int t, err = 0;

    if (var->ord != NULL) {
	C = reorder_responses(var, &err);
	if (err) {
	    return err;
	}
    }

    for (t=0; t<b->horizon; t++) {
	if (t == 0) {
	    /* initial estimated responses */
	    copy_north_west(b->rtmp, C, 0);
	} else {
	    /* calculate further estimated responses */
	    gretl_matrix_multiply(var->A, b->rtmp, b->ctmp);
	    gretl_matrix_copy_values(b->rtmp, b->ctmp);
	}
	x = gretl_matrix_get(b->rtmp, targ, shock);
	gretl_matrix_set(b->resp, t, iter, x);
    }

    if (C != var->C) {
	gretl_matrix_free(C);
    }

    return err;
}

static void maybe_resize_vecm_matrices (GRETL_VAR *v)
{
    int nc0 = (v->order * v->neqns) + v->ifc + v->jinfo->seasonals;

    if (v->xlist != NULL) {
	nc0 += v->xlist[0];
    }

    if (v->detflags & DET_TREND) {
	nc0++;
    }

    if (v->X != NULL && v->X->cols > nc0) {
	/* for stage 1: skip the extra EC-term columns in X */
	gretl_matrix_reuse(v->X, -1, nc0);
	gretl_matrix_reuse(v->B, nc0, -1);
    }
}

/* In re-estimation of VAR or VECM we'll tolerate only so many cases
   of near-perfect collinearity. These can surely arise by chance but
   if their frequency is too high this would seem to indicate an
   inherent problem.
*/

#define MAXSING 0.10 /* no more than 10 percent such cases */

static int irf_fatal (int err, irfboot *b, int iter, int scount)
{
    if (err != E_SINGULAR || iter == 0) {
	return 1;
    } else if ((scount / (double) b->iters) > MAXSING) {
	return 1;
    } else {
	return 0;
    }
}

static int
re_estimate_VECM (irfboot *b, GRETL_VAR *v, int targ, int shock,
		  int iter, int scount)
{
    static int (*jbr) (GRETL_VAR *, const DATASET *) = NULL;
    int err = 0;

    gretl_error_clear();

    if (jbr == NULL) {
	/* open the Johansen plugin */
	jbr = get_plugin_function("johansen_boot_round");
	if (jbr == NULL) {
	    return E_FOPEN;
	}
    }

    /* The various VECM matrices may need to be re-set to the sizes
       expected by johansen_stage_1() */
    maybe_resize_vecm_matrices(v);

    err = johansen_stage_1(v, b->dset, OPT_NONE, NULL);
#if BDEBUG
    fprintf(stderr, "johansen_stage_1: err = %d\n", err);
#endif

    if (!err) {
	/* call the plugin function */
	err = jbr(v, b->dset);
#if BDEBUG
	fprintf(stderr, "johansen_boot_round: err = %d\n", err);
#endif
    }

#if BDEBUG
    if (!err) {
	gretl_matrix_print(v->jinfo->Beta, "var->jinfo->Beta");
	gretl_matrix_print(v->jinfo->Alpha, "var->jinfo->Alpha");
	gretl_matrix_print(v->S, "var->S (Omega)");
	if (v->B != NULL) {
	    gretl_matrix_print(v->B, "var->B");
	}
    }
#endif

    if (!err) {
	err = gretl_VAR_do_error_decomp(v->S, v->C, v->ord);
    }

    if (!err) {
	err = recalculate_impulse_responses(b, v, targ, shock, iter);
    }

    return err;
}

static int re_estimate_VAR (irfboot *b, GRETL_VAR *v, int targ, int shock,
			    int iter)
{
    int err;

    err = gretl_matrix_multi_ols(v->Y, v->X, v->B, v->E, NULL);

    if (!err) {
	VAR_write_A_matrix(v);
    }

    if (!err) {
	gretl_matrix_multiply_mod(v->E, GRETL_MOD_TRANSPOSE,
				  v->E, GRETL_MOD_NONE,
				  v->S, GRETL_MOD_NONE);
	gretl_matrix_divide_by_scalar(v->S, v->df); /* was v->T in denom. */
	err = gretl_VAR_do_error_decomp(v->S, v->C, v->ord);
    }

    if (!err) {
	err = recalculate_impulse_responses(b, v, targ, shock, iter);
    }

    return err;
}

/* VECM: make a matrix, rbeta, containing the coefficients on the
   "restricted" constant or trend and/or restricted exogenous
   variables, if this is required.
*/

static gretl_matrix *make_restricted_coeff_matrix (const GRETL_VAR *v)
{
    gretl_matrix *b = NULL;
    gretl_matrix *rbeta = NULL;
    int rank = v->jinfo->rank;
    int nr = n_restricted_terms(v);
    double x;
    int j, k, err = 0;

    b = gretl_matrix_alloc(rank, nr);
    if (b == NULL) {
	return NULL;
    }

    for (j=0; j<rank; j++) {
	/* extract the last row(s) of \beta, transposed */
	for (k=0; k<nr; k++) {
	    x = gretl_matrix_get(v->jinfo->Beta, v->neqns + k, j);
	    gretl_matrix_set(b, j, k, x);
	}
    }

    rbeta = gretl_matrix_multiply_new(v->jinfo->Alpha, b, &err);

#if BDEBUG > 1
    gretl_matrix_print(b, "restricted var row(s) of beta'");
    gretl_matrix_print(rbeta, "coeffs for restricted term(s)");
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

gretl_matrix *VAR_coeff_matrix_from_VECM (GRETL_VAR *var)
{
    gretl_matrix *C0 = NULL;
    gretl_matrix *rbeta = NULL;
    gretl_matrix *B = var->B;
    int order = var->order + 1;
    int nexo = (var->xlist != NULL)? var->xlist[0] : 0;
    int utrend = jcode(var) == J_UNREST_TREND;
    int ndelta = var->order * var->neqns;
    int nseas = var->jinfo->seasonals;
    int nr = n_restricted_terms(var);
    int ncoeff, nb;
    int X0, S0, T0;
    double aij;
    int i, j, k;

    /* total coeffs in VAR representation */
    ncoeff = var->ncoeff + (var->neqns - var->jinfo->rank) + nr;

    /* Is var->B required? And if so how many B rows do we need? */
    nb = var->ifc + nexo + nseas + utrend;

    if (nb > 0) {
	if (B == NULL) {
	    gretl_errmsg_set("VAR coefficient matrix B is missing!");
	    return NULL;
	} else if (B->rows != nb && B->rows != var->ncoeff) {
	    gretl_errmsg_set("VAR coefficient matrix B is of wrong size!");
	    fprintf(stderr, "B should have %d or %d rows, but has %d\n",
		    nb, var->ncoeff, B->rows);
	    gretl_matrix_print(B, "var->B");
	    return NULL;
	} else {
	    nb = B->rows;
	}
    }

    if (nr > 0) {
	rbeta = make_restricted_coeff_matrix(var);
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
	const double *coeff = NULL;
	int col = 0;

	if (nb > 0) {
	    coeff = B->val + i * nb;
	}

	/* constant, if present */
	if (var->ifc) {
	    gretl_matrix_set(C0, i, col++, coeff[0]);
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
	    gretl_matrix_set(C0, i, col++, coeff[X0+j]);
	}

	/* seasonals, if present */
	for (j=0; j<nseas; j++) {
	    gretl_matrix_set(C0, i, col++, coeff[S0+j]);
	}

	/* unrestricted trend, if present */
	if (utrend) {
	    gretl_matrix_set(C0, i, col++, coeff[T0]);
	}

	/* additional restricted term(s), if present */
	if (rbeta != NULL) {
	    for (j=0; j<nr; j++) {
		aij = gretl_matrix_get(rbeta, i, j);
		gretl_matrix_set(C0, i, col++, aij);
	    }
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

/* VECM: copy levels of Y vars from main Z into temporary dataset (for
   passing to Johansen stage 1).  If there are restricted exogenous
   vars they must be included in the temp dataset too.
*/

static int init_VECM_dataset (irfboot *b, GRETL_VAR *var,
			      const DATASET *dset)
{
    int nv = var->neqns + 1;
    int i, j, vi, t;

    if (var->rlist != NULL) {
	nv += var->rlist[0];
    }

    b->dset = create_auxiliary_dataset(nv, dset->n, 0);
    if (b->dset == NULL) {
	return E_ALLOC;
    }

    copy_dataset_obs_info(b->dset, dset);
    b->dset->t1 = dset->t1;
    b->dset->t2 = dset->t2;

    /* copy levels of Y into boot->Z and adjust ylist */
    for (i=0, j=1; i<var->neqns; i++, j++) {
	vi = var->ylist[i+1];
	for (t=0; t<dset->n; t++) {
	    b->dset->Z[j][t] = dset->Z[vi][t];
	}
	var->ylist[j] = j;
    }

    if (var->rlist != NULL) {
	/* copy restricted X into boot->Z and adjust rlist */
	for (i=1; i<=var->rlist[0]; i++, j++) {
	    vi = var->rlist[i];
	    for (t=0; t<dset->n; t++) {
		b->dset->Z[j][t] = dset->Z[vi][t];
	    }
	    var->rlist[i] = j;
	}
    }

    return 0;
}

/* Compute VECM Y (in levels) using the VAR representation. */

static void
compute_VECM_dataset (irfboot *b, GRETL_VAR *var, int iter)
{
    int order = var->order + 1;
    int nexo = (var->xlist != NULL)? var->xlist[0] : 0;
    int nseas = var->jinfo->seasonals;
    double cij, eti, xti;
    int i, j, k, s, t;

#if BDEBUG
    fprintf(stderr, "compute_VECM_dataset: order=%d, nexo=%d, nseas=%d, t1=%d\n",
	    order, nexo, nseas, var->t1);
    if (var->X != NULL) {
	gretl_matrix_print(var->X, "var->X, before resampling");
    }
#endif

    for (t=var->t1, s=0; t<=var->t2; t++, s++) {
	for (i=0; i<var->neqns; i++) {
	    double bti = 0.0;
	    int xcol = var->ifc + var->neqns * var->order;
	    int col = 0;

	    /* unrestricted constant, if present */
	    if (var->ifc) {
		cij = gretl_matrix_get(b->C0, i, col++);
		bti += cij;
	    }

	    /* lags of endogenous vars */
	    for (j=1; j<=var->neqns; j++) {
		for (k=1; k<=order; k++) {
		    cij = gretl_matrix_get(b->C0, i, col++);
		    bti += cij * b->dset->Z[j][t-k];
		}
	    }

	    /* exogenous vars, if present */
	    for (j=0; j<nexo; j++) {
		cij = gretl_matrix_get(b->C0, i, col++);
		xti = gretl_matrix_get(var->X, s, xcol++);
		bti += cij * xti;
	    }

	    /* seasonals, if present */
	    for (j=0; j<nseas; j++) {
		cij = gretl_matrix_get(b->C0, i, col++);
		xti = gretl_matrix_get(var->X, s, xcol++);
		bti += cij * xti;
	    }

	    if (jcode(var) == J_UNREST_TREND) {
		/* unrestricted trend */
		cij = gretl_matrix_get(b->C0, i, col++);
		bti += cij * (t + 1);
	    } else if (jcode(var) == J_REST_CONST) {
		/* restricted constant */
		bti += gretl_matrix_get(b->C0, i, col++);
	    } else if (jcode(var) == J_REST_TREND) {
		/* restricted trend */
		cij = gretl_matrix_get(b->C0, i, col++);
		bti += cij * t;
	    }

	    /* restricted exogenous vars, if present */
	    if (var->rlist != NULL) {
		for (j=1; j<=var->rlist[0]; j++) {
		    cij = gretl_matrix_get(b->C0, i, col++);
		    bti += cij * b->dset->Z[j+var->neqns][t-1];
		}
	    }

	    /* set level of Y(t, i) to fitted value + re-sampled error */
	    eti = gretl_matrix_get(b->rE, s, i);
	    b->dset->Z[i+1][t] = bti + eti;
	}
    }

#if BDEBUG > 1
    fprintf(stderr, "VECM: recomputed levels\n\n");
    for (t=0; t<b->dset->n; t++) {
	for (i=1; i<=var->neqns; i++) {
	    fprintf(stderr, "%12.5g", b->dset->Z[i][t]);
	}
	fputc('\n', stderr);
    }
#endif

    if (var->X != NULL) {
	/* now rewrite lagged differences into X matrix */
	k = var->ifc;
	for (i=1; i<=var->neqns; i++) {
	    for (j=1; j<=var->order; j++) {
		s = 0;
		for (t=var->t1; t<=var->t2; t++) {
		    xti = b->dset->Z[i][t-j] - b->dset->Z[i][t-j-1];
		    gretl_matrix_set(var->X, s++, k, xti);
		}
		k++;
	    }
	}
#if BDEBUG > 1
	gretl_matrix_print(var->X, "var->X (vecm, resampled)");
#endif
    }
}

/* (Re-)fill the bootstrap dataset with artificial data, based on the
   re-sampled residuals from the original VAR (case of simple VAR, not
   VECM).
*/

static void compute_VAR_dataset (irfboot *b, GRETL_VAR *var,
				 const GRETL_VAR *vbak)
{
    double x;
    int i, j, k, t;
    int nl = var_n_lags(var);

#if BDEBUG
    gretl_matrix_print(var->Y, "var->Y before resampling");
    gretl_matrix_print(var->X, "var->X before resampling");
#endif

    for (t=0; t<var->T; t++) {
	/* extract row of var->X at t */
	gretl_matrix_extract_matrix(b->Xt, var->X, t, 0, GRETL_MOD_NONE);

	/* multiply Xt into original coeff matrix, forming Yt */
	gretl_matrix_multiply(b->Xt, vbak->B, b->Yt);

	/* extract resampled residuals at t */
	gretl_matrix_extract_matrix(b->Et, b->rE, t, 0, GRETL_MOD_NONE);

	/* add resampled residual to Yt */
	gretl_matrix_add_to(b->Yt, b->Et);

	/* write into big Y matrix */
	gretl_matrix_inscribe_matrix(var->Y, b->Yt, t, 0, GRETL_MOD_NONE);

	/* revise lagged Y columns in X */
	k = var->ifc;
	for (i=0; i<var->neqns; i++) {
	    x = b->Yt->val[i];
	    for (j=1; j<=nl && t+j < var->T; j++) {
		gretl_matrix_set(var->X, t+j, k++, x);
	    }
	}
    }

#if BDEBUG > 1
    gretl_matrix_print(var->Y, "var->Y after resampling");
    gretl_matrix_print(var->X, "var->X after resampling");
#endif
}

/* Resample the original VAR or VECM residuals, stored in
   vbak->E, writing the new sample into b->rE.

   Note the option to "resample" _without_ actually changing the
   order, if BDEBUG > 1.  This is useful for checking that the IRF
   bootstrap rounds are idempotent: we should then get exactly the
   same set of responses as in the original estimation of the
   VAR/VECM.
*/

static void irf_resample_resids (irfboot *b, const GRETL_VAR *vbak)
{
    double eti;
    int i, t;

    /* construct sampling array */

    for (t=0; t<vbak->T; t++) {
#if BDEBUG > 1
	b->sample[t] = t; /* fake it */
#else
	b->sample[t] = gretl_rand_int_max(vbak->T);
#endif
#if 0
	fprintf(stderr, "boot->sample[%d] = %d\n", t, b->sample[t]);
#endif
    }

    /* draw from the original residuals */

    for (t=0; t<vbak->T; t++) {
	for (i=0; i<vbak->neqns; i++) {
	    eti = gretl_matrix_get(vbak->E, b->sample[t], i);
	    gretl_matrix_set(b->rE, t, i, eti);
	}
    }
}

static int irf_boot_quantiles (irfboot *b, gretl_matrix *R, double alpha)
{
    double *rk;
    int k, ilo, ihi;

    rk = malloc(b->iters * sizeof *rk);
    if (rk == NULL) {
	return E_ALLOC;
    }

#if BDEBUG
    ilo = 1;
    ihi = b->iters;
    fprintf(stderr, "IRF bootstrap (%d iters), min and max values\n", b->iters);
#else
    ilo = (b->iters + 1) * alpha / 2.0;
    ihi = (b->iters + 1) * (1.0 - alpha / 2.0);
#endif

    for (k=0; k<b->horizon; k++) {
	gretl_matrix_row_to_array(b->resp, k, rk);
	qsort(rk, b->iters, sizeof *rk, gretl_compare_doubles);
	gretl_matrix_set(R, k, 1, rk[ilo-1]);
	gretl_matrix_set(R, k, 2, rk[ihi-1]);
    }

#if BDEBUG
    gretl_matrix_print(R, "Resp");
#endif

    free(rk);

    return 0;
}

static void restore_VAR_data (GRETL_VAR *v, GRETL_VAR *vbak)
{
    gretl_matrix_replace(&v->Y, vbak->Y);
    gretl_matrix_replace(&v->X, vbak->X);
    gretl_matrix_replace(&v->B, vbak->B);
    gretl_matrix_replace(&v->XTX, vbak->XTX);
    gretl_matrix_replace(&v->A, vbak->A);
    gretl_matrix_replace(&v->E, vbak->E);
    gretl_matrix_replace(&v->C, vbak->C);
    gretl_matrix_replace(&v->S, vbak->S);

    if (vbak->jinfo != NULL) {
	gretl_matrix_replace(&v->jinfo->R0, vbak->jinfo->R0);
	gretl_matrix_replace(&v->jinfo->R1, vbak->jinfo->R1);
	gretl_matrix_replace(&v->jinfo->S00, vbak->jinfo->S00);
	gretl_matrix_replace(&v->jinfo->S11, vbak->jinfo->S11);
	gretl_matrix_replace(&v->jinfo->S01, vbak->jinfo->S01);
	gretl_matrix_replace(&v->jinfo->Beta, vbak->jinfo->Beta);
	gretl_matrix_replace(&v->jinfo->Alpha, vbak->jinfo->Alpha);

	if (vbak->ylist != NULL) {
	    free(v->ylist);
	    v->ylist = vbak->ylist;
	}

	if (vbak->rlist != NULL) {
	    free(v->rlist);
	    v->rlist = vbak->rlist;
	}

	free(vbak->jinfo);
    }

    free(vbak);
}

/* make a back-up copy of the VAR data that we'll want to restore on
   exit from the IRF bootstrap operation */

static void matrix_swap_copy (gretl_matrix **targ,
			      gretl_matrix **src)
{
    *targ = *src;
    *src = gretl_matrix_copy(*targ);
}

static int list_swap_copy (int **targ, int **src)
{
    *targ = *src;
    *src = gretl_list_copy(*targ);
    return *src == NULL ? E_ALLOC : 0;
}

static GRETL_VAR *back_up_VAR (GRETL_VAR *v)
{
    GRETL_VAR *vbak;
    int err = 0;

    vbak = malloc(sizeof *vbak);
    if (vbak == NULL) {
	return NULL;
    }

    gretl_VAR_clear(vbak);
    clear_gretl_matrix_err();

    if (v->Y != NULL) {
	matrix_swap_copy(&vbak->Y, &v->Y);
    }
    if (v->B != NULL) {
	matrix_swap_copy(&vbak->B, &v->B);
    }

    if (v->X != NULL) {
	if (v->xcols > v->X->cols) {
	    int save_k = v->X->cols;

	    vbak->X = v->X;
	    gretl_matrix_reuse(vbak->X, -1, v->xcols);
	    v->X = gretl_matrix_copy(vbak->X);
	    gretl_matrix_reuse(vbak->X, -1, save_k);
	    gretl_matrix_reuse(v->X, -1, save_k);
	} else {
	    matrix_swap_copy(&vbak->X, &v->X);
	}
    }

    if (v->XTX != NULL) {
	matrix_swap_copy(&vbak->XTX, &v->XTX);
    }
    matrix_swap_copy(&vbak->A, &v->A);
    matrix_swap_copy(&vbak->E, &v->E);
    matrix_swap_copy(&vbak->C, &v->C);
    matrix_swap_copy(&vbak->S, &v->S);

    err = get_gretl_matrix_err();

    if (!err && v->jinfo != NULL) {
	vbak->jinfo = calloc(1, sizeof *vbak->jinfo);
	if (vbak->jinfo == NULL) {
	    err = E_ALLOC;
	} else {
	    err = list_swap_copy(&vbak->ylist, &v->ylist);
	}
    }

    if (!err && v->rlist != NULL) {
	err = list_swap_copy(&vbak->rlist, &v->rlist);
    }

    if (!err && vbak->jinfo != NULL) {
	matrix_swap_copy(&vbak->jinfo->R0, &v->jinfo->R0);
	matrix_swap_copy(&vbak->jinfo->R1, &v->jinfo->R1);
	matrix_swap_copy(&vbak->jinfo->S00, &v->jinfo->S00);
	matrix_swap_copy(&vbak->jinfo->S11, &v->jinfo->S11);
	matrix_swap_copy(&vbak->jinfo->S01, &v->jinfo->S01);
	matrix_swap_copy(&vbak->jinfo->Beta, &v->jinfo->Beta);
	matrix_swap_copy(&vbak->jinfo->Alpha, &v->jinfo->Alpha);
	err = get_gretl_matrix_err();
    }

    if (err) {
	restore_VAR_data(v, vbak);
	gretl_VAR_free(vbak);
	vbak = NULL;
    } else {
	vbak->T = v->T;
	vbak->neqns = v->neqns;
	vbak->order = v->order;
    }

    return vbak;
}

/* public bootstrapping function, called from var.c */

gretl_matrix *irf_bootstrap (GRETL_VAR *var,
			     int targ, int shock,
			     int periods, double alpha,
			     const DATASET *dset,
			     int *err)
{
    gretl_matrix *R = NULL; /* the return value */
    GRETL_VAR *vbak = NULL;
    irfboot *boot = NULL;
    int scount = 0;
    int iter;

    if (0 && (var->X == NULL || var->Y == NULL)) {
	gretl_errmsg_set("X and/or Y matrix missing, can't do this");
	*err = E_DATA;
	return NULL;
    }

#if BDEBUG
    fprintf(stderr, "\n*** irf_bootstrap() called: var->Y = %p\n",
	    (void *) var->Y);
    gretl_matrix_print(var->A, "initial var->A");
#endif

    R = gretl_zero_matrix_new(periods, 3);
    if (R == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    vbak = back_up_VAR(var);
    if (vbak == NULL) {
	*err = E_ALLOC;
	gretl_matrix_free(R);
	return NULL;
    }

    boot = irf_boot_new(var, periods);
    if (boot == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    if (var->ci == VECM) {
	boot->C0 = VAR_coeff_matrix_from_VECM(var);
	if (boot->C0 == NULL) {
	    *err = E_ALLOC;
	} else {
	    *err = init_VECM_dataset(boot, var, dset);
	}
    }

#if 0 /* just checking, for mild debugging */
    fprintf(stderr, "boot->iters = %d\n", boot->iters);
#endif

    for (iter=0; iter<boot->iters && !*err; iter++) {
#if BDEBUG
	fprintf(stderr, "starting iteration %d\n", iter);
#endif
	irf_resample_resids(boot, vbak);
	if (var->ci == VECM) {
	    compute_VECM_dataset(boot, var, iter);
	    *err = re_estimate_VECM(boot, var, targ, shock, iter, scount);
#if BDEBUG
	    if (*err) {
		fprintf(stderr, " got err = %d from re_estimate_VECM\n", *err);
	    }
#endif
	} else {
	    compute_VAR_dataset(boot, var, vbak);
	    *err = re_estimate_VAR(boot, var, targ, shock, iter);
	}
	if (*err && !irf_fatal(*err, boot, iter, scount)) {
	    /* excessive collinearity: try again, unless this is
	       becoming a serious habit
	    */
	    scount++;
	    iter--;
	    *err = 0;
	}
    }

    if (*err && scount / (double) boot->iters >= MAXSING) {
	gretl_errmsg_set("Excessive collinearity in resampled datasets");
    }

    if (!*err) {
	*err = irf_boot_quantiles(boot, R, alpha);
    }

    irf_boot_free(boot);

 bailout:

    restore_VAR_data(var, vbak);

    if (*err) {
	gretl_matrix_free(R);
	R = NULL;
    }

    return R;
}
