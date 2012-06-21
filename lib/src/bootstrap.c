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
#include "libset.h"
#include "gretl_restrict.h"
#include "gretl_xml.h"
#include "bootstrap.h"

#define BDEBUG 0

enum {
    BOOT_CI          = 1 << 0,  /* compute confidence interval */
    BOOT_PVAL        = 1 << 1,  /* compute p-value */
    BOOT_RESAMPLE_U  = 1 << 2,  /* resample the empirical residuals */
    BOOT_NORMAL_U    = 1 << 3,  /* simulate normal residuals */
    BOOT_STUDENTIZE  = 1 << 4,  /* studentize, when doing confidence interval */
    BOOT_GRAPH       = 1 << 5,  /* graph the distribution */
    BOOT_RESTRICT    = 1 << 6,  /* called via "restrict" command */
    BOOT_F_FORM      = 1 << 7,  /* compute F-statistics */
    BOOT_FREE_RQ     = 1 << 8,  /* free restriction matrices */
    BOOT_SAVE        = 1 << 9,  /* save results vector */
    BOOT_VERBOSE     = 1 << 10, /* verbose output */
    BOOT_SILENT      = 1 << 11  /* suppress printed output */
};

#define resampling(b) (b->flags & BOOT_RESAMPLE_U)
#define verbose(b) (b->flags & BOOT_VERBOSE)
#define boot_Ftest(b) (b->flags & BOOT_F_FORM)

typedef struct boot_ boot;
typedef struct ldvinfo_ ldvinfo;

struct boot_ {
    int flags;          /* option flags */
    int B;              /* number of replications */
    int k;              /* number of coefficients in model */
    int T;              /* number of observations used */
    int p;              /* index number of coeff to examine */
    int g;              /* number of restrictions */
    int mci;            /* model command index */
    gretl_matrix *y;    /* holds original, then artificial, dep. var. */
    gretl_matrix *X;    /* independent variables */
    gretl_matrix *b0;   /* coefficients used to generate dep var */
    gretl_matrix *u0;   /* original residuals for resampling */
    gretl_matrix *R;    /* LHS restriction matrix */
    gretl_matrix *q;    /* RHS restriction matrix */
    gretl_matrix *w;    /* weights for WLS */
    double SE;          /* original std. error of residuals */
    double point;       /* point estimate of coeff */
    double se0;         /* original std. error for coeff of interest */
    double test0;       /* original test statistic */
    double b_p;         /* test value for coeff */
    double a;           /* alpha, for confidence interval */
    double pval;        /* p-value for test */
    char vname[VNAMELEN]; /* name of variable analysed */
    ldvinfo *ldv;       /* lagged depvar info, if applicable */
};

struct ldvinfo_ {
    int n;      /* number of lagged dependent variable terms */
    int *xcol;  /* location of such terms (0-based column of X) */
    int *lag;   /* lag order of each such term */
};

static gretl_vector *bs_data;
static char bs_vname[VNAMELEN];

static void free_ldvinfo (ldvinfo *ldv);

static void boot_destroy (boot *bs)
{
    if (bs == NULL) {
	return;
    }

    gretl_matrix_free(bs->y);
    gretl_matrix_free(bs->X);
    gretl_matrix_free(bs->b0);
    gretl_matrix_free(bs->u0);
    gretl_matrix_free(bs->w);

    if (bs->flags & BOOT_FREE_RQ) {
	gretl_matrix_free(bs->R);
	gretl_matrix_free(bs->q);
    }

    free_ldvinfo(bs->ldv);

    free(bs);
}

static int 
make_model_matrices (boot *bs, const MODEL *pmod, 
		     const DATASET *dset)
{
    double xti;
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    int needw = 0;
    int i, s, t;

    bs->y = gretl_column_vector_alloc(T);
    bs->X = gretl_matrix_alloc(T, k);
    bs->b0 = gretl_column_vector_alloc(k);
    bs->u0 = gretl_column_vector_alloc(T);

    if (pmod->ci == WLS || pmod->ci == HSK) {
	needw = 1;
	bs->w = gretl_column_vector_alloc(T);
    }

    if (bs->y == NULL || bs->X == NULL || bs->b0 == NULL || 
	bs->u0 == NULL || (needw && bs->w == NULL)) {
	return E_ALLOC;
    }

    for (i=0; i<k; i++) {
	gretl_vector_set(bs->b0, i, pmod->coeff[i]);
    }    

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    bs->y->val[s] = dset->Z[pmod->list[1]][t];
	    if (pmod->ci == WLS) {
		bs->w->val[s] = sqrt(dset->Z[pmod->nwt][t]);
		bs->y->val[s] *= bs->w->val[s];
		bs->u0->val[s] = bs->y->val[s];
	    } else {
		gretl_vector_set(bs->u0, s, pmod->uhat[t]);
	    }
	    for (i=2; i<=pmod->list[0]; i++) {
		xti = dset->Z[pmod->list[i]][t];
		if (pmod->ci == WLS) {
		    xti *= bs->w->val[s];
		    bs->u0->val[s] -= bs->b0->val[i-2] * xti;
		}
		gretl_matrix_set(bs->X, s, i-2, xti);
	    }
	    s++;
	}
    }

    return 0;
}

static int make_bs_flags (gretlopt opt)
{
    int flags = 0;

    if (opt & OPT_P) {
	flags |= BOOT_PVAL;
    } else {
	flags |= BOOT_CI;
    }

    if (opt & OPT_N) {
	flags |= BOOT_NORMAL_U;
    } else {
	flags |= BOOT_RESAMPLE_U;
    }

    if (opt & OPT_G) {
	flags |= BOOT_GRAPH;
    }

    if (opt & OPT_T) {
	flags |= BOOT_STUDENTIZE;
    }

    if (opt & OPT_R) {
	flags |= BOOT_RESTRICT;
    }

    if (opt & OPT_F) {
	flags |= BOOT_F_FORM;
    }

    if (opt & OPT_A) {
	flags |= BOOT_SAVE;
    }

    if (opt & OPT_S) {
	flags |= BOOT_SILENT;
    } else if (opt & OPT_V) {
	flags |= BOOT_VERBOSE;
    }

    return flags;
}

/* check in with libset for current default number of replications if
   need be; in addition, alpha * (B + 1) should be an integer, when
   constructing confidence intervals
*/

int maybe_adjust_B (int B, double a, int flags)
{
    if (B <= 0) {
	B = libset_get_int(BOOTREP);
    }

    if (flags & BOOT_CI) {
	double x;

	if (B % 10 == 0) {
	    B--;
	}

	x = a * (B + 1);
	while (x - floor(x) > 1e-13) {
	    x = a * (++B + 1);
	}
    }

    return B;
}

static void free_ldvinfo (ldvinfo *ldv)
{
    if (ldv != NULL) {
	free(ldv->xcol);
	free(ldv->lag);
	free(ldv);
    }
}

/* add information on lagged dependent variables among the
   regressors in a model to be subject to bootstrap
   analysis
*/

static int add_ldvinfo (boot *b, int n)
{
    ldvinfo *ldv;

    ldv = malloc(sizeof *ldv);
    if (ldv == NULL) {
	return E_ALLOC;
    }

    ldv->xcol = malloc(n * sizeof *ldv->xcol);
    ldv->lag = malloc(n * sizeof *ldv->lag);

    if (ldv->xcol == NULL || ldv->lag == NULL) {
	free_ldvinfo(ldv);
	return E_ALLOC;
    }

    ldv->n = n;
    b->ldv = ldv;

    return 0;
}

static int make_boot_ldvinfo (boot *b, const MODEL *pmod,
			      const DATASET *dset)
{
    int xnum, ynum = pmod->list[1];
    int i, p, nly = 0;
    int leads = 0;
    int err = 0;

    for (i=0; i<pmod->ncoeff && !err; i++) {
	xnum = pmod->list[i+2];
	p = standard_lag_of(xnum, ynum, dset);
	if (p > 0) {
	    nly++;
	} else if (p < 0) {
	    leads++;
	}
    }

    if (nly > 0 && leads > 0) {
	gretl_errmsg_set("model contains leads of the dependent "
			 "variable: not supported");
	err = E_DATA;
    }	

    if (!err && nly > 0) {
	int k = 0;

	err = add_ldvinfo(b, nly);

	if (!err) {
	    for (i=0; i<pmod->ncoeff; i++) {
		xnum = pmod->list[i+2];
		p = standard_lag_of(xnum, ynum, dset);
		if (p > 0) {
		    b->ldv->xcol[k] = i;
		    b->ldv->lag[k] = p;
		    k++;
		}
	    }	    
	}
    }

    return err;
}

static boot *boot_new (const MODEL *pmod,
		       const DATASET *dset,
		       int B, gretlopt opt,
		       int *err)
{
    boot *bs;

    bs = malloc(sizeof *bs);
    if (bs == NULL) {
	return NULL;
    }

    bs->ldv = NULL;

    *err = make_boot_ldvinfo(bs, pmod, dset);
    if (*err) {
	free(bs);
	return NULL;
    }

    bs->y = NULL;
    bs->X = NULL;
    bs->b0 = NULL;
    bs->u0 = NULL;
    bs->w = NULL;

    bs->R = NULL;
    bs->q = NULL;

    if (make_model_matrices(bs, pmod, dset)) {
	boot_destroy(bs);
	return NULL;
    }

    bs->flags = make_bs_flags(opt);

    bs->mci = pmod->ci;
    bs->a = 0.05; /* make configurable? */
    bs->B = maybe_adjust_B(B, bs->a, bs->flags);

    bs->p = 0;
    bs->g = 0;

    *bs->vname = '\0';

    bs->SE = NADBL;
    bs->point = NADBL;
    bs->se0 = NADBL;
    bs->test0 = NADBL;
    bs->b_p = NADBL;
    bs->pval = NADBL;

    bs->k = bs->X->cols;
    bs->T = bs->X->rows;

    return bs;
}

static int ldv_lag (boot *bs, int k)
{
    if (bs->ldv != NULL) {
	int i;

	for (i=0; i<bs->ldv->n; i++) {
	    if (bs->ldv->xcol[i] == k) {
		return bs->ldv->lag[i];
	    }
	}
    }

    return 0;
}

static void make_normal_y (boot *bs)
{
    double xti;
    int i, t, p;

    /* generate scaled normal errors */
    gretl_matrix_random_fill(bs->y, D_NORMAL);
    gretl_matrix_multiply_by_scalar(bs->y, bs->SE);

    /* construct y recursively */
    for (t=0; t<bs->X->rows; t++) {
	for (i=0; i<bs->X->cols; i++) {
	    p = ldv_lag(bs, i);
	    if (p > 0 && t >= p) {
		gretl_matrix_set(bs->X, t, i, bs->y->val[t-p]);
	    } 
	    xti = gretl_matrix_get(bs->X, t, i);
	    bs->y->val[t] += bs->b0->val[i] * xti;
	}
    }  	
}

static void 
resample_vector (const gretl_matrix *u0, gretl_matrix *u, int *z)
{
    int t, T = u->rows;

    /* generate T uniform drawings from [0 .. T-1] */
    gretl_rand_int_minmax(z, T, 0, T-1);

    /* sample from source vector based on indices */
    for (t=0; t<T; t++) {
	u->val[t] = u0->val[z[t]];
    }
}

static void make_resampled_y (boot *bs, int *z)
{
    double xti;
    int i, t, p;

    /* resample the residuals, into y */
    resample_vector(bs->u0, bs->y, z);

    /* construct y recursively */
    for (t=0; t<bs->X->rows; t++) {
	for (i=0; i<bs->X->cols; i++) {
	    p = ldv_lag(bs, i);
	    if (p > 0 && t >= p) {
		gretl_matrix_set(bs->X, t, i, bs->y->val[t-p]);
	    }
	    xti = gretl_matrix_get(bs->X, t, i);
	    bs->y->val[t] += bs->b0->val[i] * xti;
	}
    }
}

/* when doing a bootstrap p-value: run the restricted regression; save
   the coefficient vector in bs->b0 and residuals in bs->u0
*/

static int do_restricted_ols (boot *bs)
{
    double s2 = 0.0;
    int err = 0;
    
    err = gretl_matrix_restricted_ols(bs->y, bs->X, bs->R, bs->q, bs->b0, 
				      NULL, bs->u0, &s2);

#if BDEBUG
    if (1) {
	int i;

	fprintf(stderr, "Restricted estimates (err = %d):\n", err);
	for (i=0; i<bs->b0->rows; i++) {
	    fprintf(stderr, "b[%d] = %g\n", i, bs->b0->val[i]);
	}
	fprintf(stderr, "s2 = %g\n", s2);
    }
#endif

    if (!err) {
	bs->SE = sqrt(s2);
    }

    return err;
}

static int bs_store_result (boot *bs, double *xi)
{
    int i;

    if (bs_data != NULL) {
	gretl_matrix_free(bs_data);
	bs_data = NULL;
    }

    bs_data = gretl_column_vector_alloc(bs->B);
    if (bs_data == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<bs->B; i++) {
	gretl_vector_set(bs_data, i, xi[i]);
    }

    *bs_vname = '\0';

    if (bs->flags & BOOT_F_FORM) {
	strcpy(bs_vname, "F_test");
    } else {
	if (bs->flags & (BOOT_PVAL | BOOT_STUDENTIZE)) {
	    strcpy(bs_vname, "t_");
	} else {
	    strcpy(bs_vname, "b_");
	}
	strncat(bs_vname, bs->vname, VNAMELEN - 3);
    }

    return 0;
}

static void bs_print_result (boot *bs, double *xi, int tail, PRN *prn)
{
    if (bs->flags & BOOT_F_FORM) {
	pprintf(prn, "\n%s: ", _("bootstrap F-test"));
    } else if (bs->flags & BOOT_RESTRICT) {
	pputs(prn, "\n  ");
    } else {
	pprintf(prn, _("For the coefficient on %s (point estimate %g)"), 
		bs->vname, bs->point);
	pputs(prn, ":\n\n  ");
    }

    if (bs->flags & (BOOT_CI | BOOT_GRAPH)) {
	qsort(xi, bs->B, sizeof *xi, gretl_compare_doubles);
    }

    if (bs->flags & BOOT_PVAL) {
	pprintf(prn, "%s = %d / %d = %g", _("p-value"), tail, bs->B, bs->pval);
    } else {
	double ql, qu;
	int i, j;

	i = bs->a * (bs->B + 1) / 2.0;
	ql = xi[i-1];

	j = bs->B - i + 1;
	qu = xi[j-1];

	if (bs->flags & BOOT_STUDENTIZE) {
	    double cl = ql;

	    ql = bs->point - bs->se0 * qu;
	    qu = bs->point - bs->se0 * cl;
	    pprintf(prn, _("Studentized %g%% confidence interval = %g to %g"), 
		    100 * (1 - bs->a), ql, qu);
	} else {
	    pprintf(prn, _("%g%% confidence interval = %g to %g"), 
		    100 * (1 - bs->a), ql, qu);
	}
    }
    
    pputs(prn, "\n\n");
    pprintf(prn, _("Based on %d replications"), bs->B);
    pputs(prn, ", ");
    if (bs->flags & BOOT_RESAMPLE_U) {
	pputs(prn, _("using resampled residuals"));
    } else {
	pputs(prn, _("with simulated normal errors"));
    }
    pputc(prn, '\n');

    if (bs->ldv != NULL) {
	pprintf(prn, "(%s)",  _("recognized lagged dependent variable"));
	pputc(prn, '\n');
    }

    if (bs->flags & BOOT_GRAPH) {
	int (*kdfunc) (const double *, int, const char *);
	void *handle;
	char label[48];

	kdfunc = get_plugin_function("array_kernel_density", &handle);
	if (kdfunc == NULL) {
	    return;
	}

	if (bs->flags & BOOT_F_FORM) {
	    strcpy(label, _("bootstrap F-test"));
	} else if (bs->flags & (BOOT_PVAL | BOOT_STUDENTIZE)) {
	    strcpy(label, _("bootstrap t-ratio"));
	} else {
	    strcpy(label, _("bootstrap coefficient"));
	} 

	(*kdfunc)(xi, bs->B, label);
	close_plugin(handle);
    }
}

/* Davidson and MacKinnon, ETM, p. 163 */

static void rescale_residuals (boot *bs)
{
    double s;
    int k = bs->k;

    if (bs->flags & BOOT_PVAL) {
	k--;
    }

    s = sqrt((double) bs->T / (bs->T - k));

    gretl_matrix_multiply_by_scalar(bs->u0, s);
}

/* we do the following on each replication if we're bootstrapping
   an F-test (not for just a single-coefficient p-value) */

static double bs_F_test (const gretl_matrix *b,
			 const gretl_matrix *V, 
			 boot *bs, int *errp)
{
    gretl_vector *br = NULL;
    gretl_matrix *RVR = NULL;
    double test = 0.0;
    int err = 0;

    br = gretl_column_vector_alloc(bs->g);
    RVR = gretl_matrix_alloc(bs->R->rows, bs->R->rows);

    if (br == NULL || RVR == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = gretl_matrix_multiply(bs->R, b, br);
    if (err) {
	fprintf(stderr, "Failed: gretl_matrix_multiply(R, b, br)\n");
	goto bailout;
    }

    if (!gretl_is_zero_matrix(bs->q)) {
	err = gretl_matrix_subtract_from(br, bs->q);
	if (err) {
	    fprintf(stderr, "Failed: gretl_matrix_subtract_from(br, q)\n");
	    goto bailout;
	}
    }

    gretl_matrix_qform(bs->R, GRETL_MOD_NONE, V,
		       RVR, GRETL_MOD_NONE);

    err = gretl_invert_symmetric_matrix(RVR);

    if (!err) {
	test = gretl_scalar_qform(br, RVR, &err);
    }

    if (!err) {
	test /= bs->g;
    }

 bailout:

    gretl_vector_free(br);
    gretl_matrix_free(RVR);

    *errp = err;

    return test;
}

#if 0

static char *squarable_cols_mask (const gretl_matrix *X, int *n)
{
    const double *xi = X->val;
    char *mask;
    int i, t;

    mask = calloc(X->cols, 1);
    if (mask == NULL) {
	return NULL;
    }

    for (i=0; i<X->cols; i++) {
	for (t=0; t<X->rows; t++) {
	    if (xi[t] != 1.0 && xi[t] != 0.0) {
		mask[i] = 1;
		*n += 1;
		break;
	    }
	}
	xi += X->rows;
    }

    return mask;
}

static int hsk_transform_data (boot *bs, gretl_matrix *b,
			       gretl_matrix *yh)
{
    gretl_matrix *X = NULL;
    gretl_matrix *g = NULL;
    int Xcols = bs->k;
    double *xi, *xj;
    double ut;
    char *mask = NULL;
    int bigX = 0;
    int i, t, err = 0;

    gretl_matrix_multiply(bs->X, b, yh);

    /* calculate log of uhat squared */
    for (t=0; t<bs->T; t++) {
	ut = bs->y->val[t] - yh->val[t];
	bs->w->val[t] = log(ut * ut);
    } 

    mask = squarable_cols_mask(bs->X, &Xcols);
    if (mask == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    X = gretl_matrix_alloc(bs->T, Xcols);
    if (X == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (Xcols > bs->k) {
	bigX = 1;
    }
	    
    g = gretl_column_vector_alloc(Xcols);
    if (g == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    xi = bs->X->val;
    xj = X->val;

    for (i=0; i<bs->k; i++) {
	for (t=0; t<bs->T; t++) {
	    xj[t] = xi[t];
	}
	xi += bs->T;
	xj += bs->T;
    }

    if (bigX) {
	xi = bs->X->val;
	for (i=0; i<bs->k; i++) {
	    if (mask[i]) {
		for (t=0; t<bs->T; t++) {
		    xj[t] = xi[t] * xi[t];
		}
		xj += bs->T;
	    }
	    xi += bs->T;
	}
    }

    err = gretl_matrix_ols(bs->w, X, g, NULL, NULL, NULL);
    if (err) {
	goto bailout;
    }

    /* form weighted dataset */

    gretl_matrix_multiply(X, g, yh);
    for (t=0; t<bs->T; t++) {
	gretl_vector_set(bs->w, t, sqrt(1.0 / exp(yh->val[t])));
    }

    gretl_matrix_reuse(X, bs->T, bs->k);
    xi = X->val;
    for (i=0; i<bs->k; i++) {
	for (t=0; t<bs->T; t++) {
	    xi[t] *= gretl_vector_get(bs->w, t);
	}
	xi += bs->T;
    }
    for (t=0; t<bs->T; t++) {
	bs->w->val[t] *= bs->y->val[t];
    }

    /* do WLS */
    err = gretl_matrix_ols(bs->w, X, b, NULL, NULL, NULL);

 bailout:

    gretl_matrix_free(g);
    gretl_matrix_free(X);
    free(mask);

    return err;
}

#endif

static void print_test_round (boot *bs, int i, double test, PRN *prn)
{
    pprintf(prn, "round %d: test = %g%s\n", i+1, test,
	    test > bs->test0 ? " *" : "");
}

/* do the actual bootstrap analysis: the objective is either to form a
   confidence interval or to compute a p-value; the methodology is
   either to resample the original residuals or to simulate normal
   errors with the empirically given variance.
*/

static int real_bootstrap (boot *bs, PRN *prn)
{
    gretl_matrix *XTX = NULL;   /* X'X */
    gretl_matrix *XTXI = NULL;  /* X'X^{-1} */
    gretl_matrix *b = NULL;     /* re-estimated coeffs */
    gretl_matrix *yh = NULL;    /* fitted values */
    gretl_matrix *V = NULL;     /* covariance matrix */
    int *z = NULL;
    double *xi = NULL;
    int k = bs->k;
    int p = bs->p;
    int tail = 0;
    int i, t, err = 0;

    if (bs->flags & BOOT_PVAL) {
	err = do_restricted_ols(bs);
	if (err) {
	    return err;
	}
    }

    b = gretl_column_vector_alloc(k);
    XTX = gretl_matrix_alloc(k, k);
    yh = gretl_column_vector_alloc(bs->T);

    if (b == NULL || XTX == NULL || yh == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (resampling(bs)) {
	/* resampling index array */
	z = malloc(bs->T * sizeof *z);
	if (z == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
	rescale_residuals(bs);
    }

    if (bs->flags & (BOOT_CI | BOOT_GRAPH | BOOT_SAVE)) {
	/* array for storing results */
	xi = malloc(bs->B * sizeof *xi);
	if (xi == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }	

    if (bs->flags & BOOT_F_FORM) {
	/* covariance matrix for F-test */
	V = gretl_matrix_alloc(XTX->rows, XTX->cols);
	if (V == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }	    

    gretl_matrix_multiply_mod(bs->X, GRETL_MOD_TRANSPOSE,
			      bs->X, GRETL_MOD_NONE,
			      XTX, GRETL_MOD_NONE);

    XTXI = gretl_matrix_alloc(XTX->rows, XTX->cols);
    if (XTXI == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = gretl_matrix_cholesky_decomp(XTX);
    if (!err) {
	err = gretl_inverse_from_cholesky_decomp(XTXI, XTX);
    }

    if (verbose(bs)) {
	if (boot_Ftest(bs)) {
	    pputc(prn, '\n');
	} else {
	    pprintf(prn, "%13s %13s\n", "b", "tval");
	}
    }

    /* carry out B replications */

    for (i=0; i<bs->B && !err; i++) {
	double vpp, se, SSR, s2, ut, test;

#if BDEBUG > 1
	fprintf(stderr, "real_bootstrap: round %d\n", i);
#endif

	if (bs->flags & BOOT_NORMAL_U) {
	    make_normal_y(bs);
	} else {
	    make_resampled_y(bs, z); 
	}

	if (bs->ldv != NULL) {
	    /* X matrix includes lagged dependent variable(s), so it has
	       to be modified */
	    int j, p;

	    for (j=0; j<bs->X->cols; j++) {
		p = ldv_lag(bs, j);
		if (p > 0) {
		    for (t=p; t<bs->T; t++) {
			gretl_matrix_set(bs->X, t, j, bs->y->val[t-p]);
		    }
		}
	    }
	    gretl_matrix_multiply_mod(bs->X, GRETL_MOD_TRANSPOSE,
				      bs->X, GRETL_MOD_NONE,
				      XTX, GRETL_MOD_NONE);
	    err = gretl_matrix_cholesky_decomp(XTX);
	    if (!err) {
		err = gretl_inverse_from_cholesky_decomp(XTXI, XTX);
	    }
	}

	if (!err) {
	    /* form X'y */
	    gretl_matrix_multiply_mod(bs->X, GRETL_MOD_TRANSPOSE,
				      bs->y, GRETL_MOD_NONE,
				      b, GRETL_MOD_NONE);
	}

	if (!err) {
	    /* solve for current parameter estimates */
	    err = gretl_cholesky_solve(XTX, b);
	}

#if 0
	if (!err && bs->mci == HSK) {
	    err = hsk_transform_data(bs, b, yh);
	}
#endif

	if (err) {
	    break;
	}	

	/* form fitted values */
	gretl_matrix_multiply(bs->X, b, yh);

	/* residuals, etc */
	SSR = 0.0;
	for (t=0; t<bs->T; t++) {
	    ut = bs->y->val[t] - yh->val[t];
	    SSR += ut * ut;
	} 
	s2 = SSR / (bs->T - k);

	/* F-test, if wanted */
	if (bs->flags & BOOT_F_FORM) {
	    gretl_matrix_copy_values(V, XTXI);
	    gretl_matrix_multiply_by_scalar(V, s2);
	    test = bs_F_test(b, V, bs, &err);
	    if (verbose(bs)) {
		print_test_round(bs, i, test, prn);
	    }
	    if (test > bs->test0) {
		tail++;
	    }
	    if (bs->flags & (BOOT_GRAPH | BOOT_SAVE)) {
		xi[i] = test;
	    }
	    continue;
	}

	/* coeff, standard error and t-stat */
	vpp = gretl_matrix_get(XTXI, p, p);
	se = sqrt(s2 * vpp);
	test = (b->val[p] - bs->b_p) / se;

	if (verbose(bs)) {
	    pprintf(prn, "%13g %13g\n", b->val[p], test);
	}

	if (bs->flags & BOOT_CI) {
	    if (bs->flags & BOOT_STUDENTIZE) {
		xi[i] = test;
	    } else {
		xi[i] = b->val[p];
	    }
	} else {
	    /* doing p-value */
	    if (bs->flags & (BOOT_GRAPH | BOOT_SAVE)) {
		xi[i] = test;
	    }
	    if (fabs(test) > fabs(bs->test0)) {
		tail++;
	    }
	}
    }

    if (!err) {
	if (bs->flags & BOOT_PVAL) {
	    bs->pval = (double) tail / bs->B;
	    record_test_result(bs->test0, bs->pval, _("bootstrap test"));
	}
	if (bs->flags & BOOT_SAVE) {
	    bs_store_result(bs, xi);
	}
	if (!(bs->flags & BOOT_SILENT)) {
	    bs_print_result(bs, xi, tail, prn);
	}
    }

 bailout:

    gretl_matrix_free(b);
    gretl_matrix_free(XTX);
    gretl_matrix_free(XTXI);
    gretl_matrix_free(yh);
    gretl_matrix_free(V);

    free(z);
    free(xi);
    
    return err;
}

/* add basic restriction matrices R and q when doing a p-value
   calculation for a single variable */

static int bs_add_restriction (boot *bs, int p)
{
    bs->R = gretl_zero_matrix_new(1, bs->b0->rows);
    bs->q = gretl_zero_matrix_new(1, 1);

    if (bs->R == NULL || bs->q == NULL) {
	return E_ALLOC;
    }

    gretl_vector_set(bs->R, p, 1.0);

    return 0;
}

/**
 * bootstrap_analysis:
 * @pmod: model to be examined.
 * @p: 0-based index number of the coefficient to analyse.
 * @B: number of replications.
 * @dset: dataset struct.
 * @opt: option flags -- may contain %OPT_P to compute p-value
 * (default is to calculate confidence interval), %OPT_N
 * to use simulated normal errors (default is to resample the
 * empirical residuals), %OPT_G to display graph, %OPT_S for
 * silent operation.
 * @prn: printing struct.
 *
 * Calculates a bootstrap confidence interval or p-value for
 * a given coefficient in a given model, estimated via OLS.
 * If the first lag of the dependent variable is present as a
 * regressor it is handled correctly but more complex lag
 * autoregressive schemes are not (yet) handled.
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int bootstrap_analysis (MODEL *pmod, int p, int B, 
			const DATASET *dset, gretlopt opt,
			PRN *prn)
{
    boot *bs = NULL;
    int err = 0;

    if (!bootstrap_ok(pmod->ci)) {
	return E_NOTIMP;
    }

    if (p < 0 || p >= pmod->ncoeff) {
	return E_DATA;
    }

    bs = boot_new(pmod, dset, B, opt, &err);

    if (!err && (bs->flags & BOOT_PVAL)) {
	err = bs_add_restriction(bs, p);
    }

    if (!err) {
	int v = pmod->list[p+2];

	bs->p = p;  /* coeff to examine */
	if (pmod->ci == HSK) {
	    bs->SE = gretl_model_get_double(pmod, "sigma_orig");
	} else {
	    bs->SE = pmod->sigma;
	}
	strcpy(bs->vname, dset->varname[v]);
	bs->point = pmod->coeff[p];
	bs->se0 = pmod->sderr[p];
	bs->test0 = pmod->coeff[p] / pmod->sderr[p];
	if (bs->flags & BOOT_PVAL) {
	    bs->b_p = 0.0;
	} else {
	    bs->b_p = bs->point;
	}
	err = real_bootstrap(bs, prn);
    }

    boot_destroy(bs);

    return err;
}

/**
 * bootstrap_test_restriction:
 * @pmod: model to be examined.
 * @R: left-hand restriction matrix, as in Rb = q.
 * @q: right-hand restriction matrix.
 * @test: initial test statistic.
 * @g: number of restrictions.
 * @dset: pointer to dataset.
 * @opt: options passed to the restrict command.
 * @prn: printing struct.
 *
 * Calculates a bootstrap p-value for the restriction on the
 * coefficients of @pmod represented by the matrices @R and @q.
 * If lags of the dependent variable are present as 
 * regressors they should be handled correctly so long as
 * the lagged terms were generated in the standard gretl manner.
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int bootstrap_test_restriction (MODEL *pmod, gretl_matrix *R, 
				gretl_matrix *q, double test, int g,
				const DATASET *dset, 
				gretlopt opt, PRN *prn)
{
    gretlopt bopt = OPT_P | OPT_R | OPT_F;
    boot *bs = NULL;
    int B = 0;
    int err = 0;

#if BDEBUG
    fprintf(stderr, "bootstrap_test_restriction: on input test = %g, g = %d\n",
	    test, g);
    gretl_matrix_print(R, "R");
    gretl_matrix_print(q, "q");
#endif    

    if (opt & OPT_S) {
	/* silent */
	bopt |= OPT_S;
    } else if (opt & OPT_V) {
	/* verbose */
	bopt |= OPT_V;
    }

    gretl_restriction_get_boot_params(&B, &bopt);

    bs = boot_new(pmod, dset, B, bopt, &err);

    if (!err) {
	bs->R = R;
	bs->q = q;
	bs->g = g;
	bs->test0 = test;
	strcpy(bs->vname, "F-test");
	err = real_bootstrap(bs, prn);
    }

    boot_destroy(bs);

    return err;
}

int bootstrap_ok (int ci)
{
    return (ci == OLS || ci == WLS); /* HSK?? */
}

int bootstrap_save_data (const char *fname)
{
    char **S = NULL;
    int err, ns = 0;

    if (bs_data == NULL) {
	return E_DATA;
    }

    err = strings_array_add(&S, &ns, bs_vname);
    if (err) {
	return err;
    }

    err = gretl_write_matrix_as_gdt(fname, bs_data, (const char **) S, 
				    NULL);

    gretl_matrix_free(bs_data);
    bs_data = NULL;
    strings_array_free(S, ns);
    *bs_vname = '\0';

    return err;
}



