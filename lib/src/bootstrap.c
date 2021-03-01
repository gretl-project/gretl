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
#include "qr_estimate.h"
#include "bootstrap.h"

#define BDEBUG 0

enum {
    BOOT_CI          = 1 << 0,  /* compute confidence interval */
    BOOT_PVAL        = 1 << 1,  /* compute p-value */
    BOOT_RESAMPLE_U  = 1 << 2,  /* resample the empirical residuals */
    BOOT_NORMAL_U    = 1 << 3,  /* simulate normal residuals */
    BOOT_PAIRS       = 1 << 4,  /* resample y and X "pairs" */
    BOOT_WILD        = 1 << 5,  /* "wild" bootstrap */
    BOOT_STUDENTIZE  = 1 << 6,  /* studentize, when doing confidence interval */
    BOOT_GRAPH       = 1 << 7,  /* graph the distribution */
    BOOT_RESTRICT    = 1 << 8,  /* called via "restrict" command */
    BOOT_F_FORM      = 1 << 9,  /* compute F-statistics */
    BOOT_FREE_RQ     = 1 << 10, /* free restriction matrices */
    BOOT_SAVE        = 1 << 11, /* save results vector */
    BOOT_VERBOSE     = 1 << 12, /* verbose output */
    BOOT_SILENT      = 1 << 13, /* suppress printed output */
    BOOT_WILD_M      = 1 << 14, /* using Mammen form for wild bootstrap */
    BOOT_HAC         = 1 << 15  /* use HAC estimator in bootstrap */
};

#define resampling_u(b)     (b->flags & BOOT_RESAMPLE_U)
#define resampling_pairs(b) (b->flags & BOOT_PAIRS)
#define resampling(b)       (b->flags & (BOOT_RESAMPLE_U | BOOT_PAIRS))
#define wild_boot(b)        (b->flags & BOOT_WILD)
#define verbose(b)          (b->flags & BOOT_VERBOSE)
#define doing_Ftest(b)      (b->flags & BOOT_F_FORM)
#define tau_wanted(b)       (b->flags & (BOOT_PVAL | BOOT_STUDENTIZE))
#define boot_use_hac(b)     (b->flags & BOOT_HAC)
#define studentizing(b)     (b->flags & BOOT_STUDENTIZE)

typedef struct boot_ boot;
typedef struct ldvinfo_ ldvinfo;

struct boot_ {
    int flags;          /* option flags */
    int B;              /* number of replications */
    int k;              /* number of coefficients in model */
    int T;              /* number of observations used */
    int p;              /* index number of specific coeff to examine */
    int g;              /* number of restrictions */
    int mci;            /* model command index */
    gretl_matrix *y;    /* holds original, then artificial, dep. var. */
    gretl_matrix *X;    /* independent variables */
    gretl_matrix *y0;   /* holds original dep var (for pairs) */
    gretl_matrix *X0;   /* original independent variables (for pairs) */
    gretl_matrix *b0;   /* coefficients used to generate dep var */
    gretl_matrix *u0;   /* original residuals for resampling */
    gretl_matrix *R;    /* LHS restriction matrix */
    gretl_matrix *q;    /* RHS restriction matrix */
    gretl_matrix *w;    /* weights for WLS */
    int hc_version;     /* HCCME variant, or -1 for none */
    int blocklen;       /* block-length, for resampling by blocks */
    double SER0;        /* original std. error of regression */
    double point;       /* point estimate of coeff */
    double bp0;         /* reference value of coefficient */
    double sep0;        /* original std. error for coeff p */
    double test0;       /* original test statistic */
    double a;           /* alpha, for confidence interval */
    double pval;        /* p-value for bootstrap test */
    char vname[VNAMELEN]; /* name of variable analysed */
    VCVInfo *vi;        /* covariance matrix info from model */
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
    gretl_matrix_free(bs->y0);
    gretl_matrix_free(bs->X0);
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
    int need0 = 0;
    int i, s, t, yno;

    bs->y = gretl_column_vector_alloc(T);
    bs->X = gretl_matrix_alloc(T, k);
    bs->b0 = gretl_column_vector_alloc(k);
    bs->u0 = gretl_column_vector_alloc(T);

    if (pmod->ci == WLS) {
	needw = 1;
	bs->w = gretl_column_vector_alloc(T);
    }

    if (resampling_pairs(bs)) {
	need0 = 1;
	bs->y0 = gretl_column_vector_alloc(T);
	bs->X0 = gretl_matrix_alloc(T, k);
    }    

    if (bs->y == NULL || bs->X == NULL ||
	bs->b0 == NULL || bs->u0 == NULL) {
	return E_ALLOC;
    }

    if ((needw && bs->w == NULL) ||
	(need0 && (bs->y0 == NULL || bs->X0 == NULL))) {
	return E_ALLOC;
    }

    for (i=0; i<k; i++) {
	gretl_vector_set(bs->b0, i, pmod->coeff[i]);
    }

    yno = pmod->list[1];

    /* transcribe the data into matrix form, skipping
       any missing values */

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    bs->y->val[s] = dset->Z[yno][t];
	    if (bs->y0 != NULL) {
		bs->y0->val[s] = dset->Z[yno][t];
	    }
	    if (pmod->ci == WLS) {
		bs->w->val[s] = sqrt(dset->Z[pmod->nwt][t]);
		bs->y->val[s] *= bs->w->val[s];
		if (bs->y0 != NULL) {
		    bs->y0->val[s] = bs->y->val[s];
		}
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
		if (bs->X0 != NULL) {
		    gretl_matrix_set(bs->X0, s, i-2, xti);
		}
	    }
	    s++;
	}
    }

    return 0;
}

static int make_bs_flags (gretlopt opt)
{
    int flags = 0;

    /* p-value versus confidence interval */
    if (opt & OPT_P) {
	flags |= BOOT_PVAL;
    } else {
	flags |= BOOT_CI;
    }

    /* boostrap method */
    if (opt & OPT_N) {
	flags |= BOOT_NORMAL_U;
    } else if (opt & OPT_X) {
	flags |= BOOT_PAIRS;
    } else if (opt & OPT_W) {
	flags |= BOOT_WILD;
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

static void maybe_set_bs_blocklen (boot *bs)
{
    int blen = 0;
    
    if (resampling_u(bs)) {
	if (bs->vi->vmin == KERNEL_QS) {
	    blen = (int) ceil(bs->vi->bw);
	} else {
	    blen = bs->vi->order;
	}
    }

    if (blen > 1) {
	fprintf(stderr, "HAC: setting blocklen = %d\n", blen);
	bs->blocklen = blen;
    }
}

static boot *boot_new (const MODEL *pmod,
		       const DATASET *dset,
		       int B, double alpha,
		       gretlopt opt, int *err)
{
    boot *bs;

    bs = malloc(sizeof *bs);
    if (bs == NULL) {
	return NULL;
    }

    bs->ldv = NULL;

    if (!(opt & OPT_X)) {
	/* we don't need to do this for the pairs bootstrap */
	*err = make_boot_ldvinfo(bs, pmod, dset);
	if (*err) {
	    free(bs);
	    return NULL;
	}
    }

    bs->y = NULL;
    bs->X = NULL;
    bs->y0 = NULL;
    bs->X0 = NULL;
    bs->b0 = NULL;
    bs->u0 = NULL;
    bs->w = NULL;

    bs->R = NULL;
    bs->q = NULL;

    bs->flags = make_bs_flags(opt);

    if (make_model_matrices(bs, pmod, dset)) {
	boot_destroy(bs);
	return NULL;
    }

    bs->mci = pmod->ci;
    bs->hc_version = gretl_model_get_hc_version(pmod);
    bs->blocklen = 0;
    bs->a = alpha;
    bs->B = maybe_adjust_B(B, bs->a, bs->flags);

    bs->vi = gretl_model_get_data(pmod, "vcv_info");

    if (bs->vi != NULL && bs->vi->vmaj == VCV_HAC) {
	bs->flags |= BOOT_HAC;
	maybe_set_bs_blocklen(bs);
    }

    bs->p = 0;
    bs->g = 0;

    *bs->vname = '\0';

    bs->SER0 = NADBL;
    bs->point = NADBL;
    bs->sep0 = NADBL;
    bs->test0 = NADBL;
    bs->bp0 = NADBL;
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
    gretl_matrix_multiply_by_scalar(bs->y, bs->SER0);

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

#define HAC_DEBUG 0

static void make_resampled_y (boot *bs, int *z)
{
    double xti;
    int i, t, p;

#if HAC_DEBUG
    /* check replication on identical data */
    return;
#endif    

    /* resample the residuals, into y */
    if (bs->blocklen > 1) {
	gretl_matrix_block_resample2(bs->y, bs->u0, bs->blocklen, z);
    } else {
	resample_vector(bs->u0, bs->y, z);
    }

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

/* Davidson-Flachaire: the bootstrap value of the 
   dependent variable is

   y^*_t = X_t\hat{\beta} + f(\hat{u}_t) v^*_t

   where f(\hat{u}_t) = \hat{u}_t / (1 - h_t)^{1/2},
   for h_t the diagonal elements of the "hat" matrix.

   The v^*_t values are either:

   1 with probability 0.5 and -1 with probability 0.5
   (Rademacher variables), or

   -(sqrt(5)-1)/2 with probability (sqrt(5)+1)/(2*sqrt(5))
   and (sqrt(5)+1)/2 with the complementary probability
   (Mammen, 1993). The latter case is, approximately:
   -0.62 with probability 0.73 and 1.62 with probability
   0.28. (Mammen, 1993)
*/

static void make_wild_y (boot *bs, int *z, double *xz)
{
    static double pminus, mminus, mplus;
    double xti;
    int i, t, p;

    if (bs->flags & BOOT_WILD_M) {
	/* Mammen */
	if (pminus == 0.0) {
	    double r5 = sqrt(5.0);

	    pminus = (r5 + 1)/(2*r5);
	    mminus = -(r5 - 1)/2.0;
	    mplus = (r5 + 1)/2.0;
	}
	gretl_rand_uniform(xz, 0, bs->T - 1);
    } else {
	/* Rademacher */
	gretl_rand_int_minmax(z, bs->T, 0, 1);
    }

    /* construct y recursively */
    for (t=0; t<bs->X->rows; t++) {
	bs->y->val[t] = bs->u0->val[t];
	if (bs->flags & BOOT_WILD_M) {
	    bs->y->val[t] *= (xz[t] < pminus)? mminus : mplus;
	} else {
	    bs->y->val[t] *= z[t] ? 1 : -1;
	}
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

static void make_resampled_pairs (boot *bs, int *z)
{
    double xti;
    int i, s, t;

    /* fill the resampling array */
    gretl_rand_int_minmax(z, bs->T, 0, bs->T-1);

    /* fill y and X with resampled "pairs" */
    for (t=0; t<bs->T; t++) {
	s = z[t];
	bs->y->val[t] = bs->y0->val[s];
	for (i=0; i<bs->X->cols; i++) {
	    xti = gretl_matrix_get(bs->X0, s, i);
	    gretl_matrix_set(bs->X, t, i, xti);
	}
    }
}

/* When computing a bootstrap p-value: the coefficients used in the
   bootstrap DGP should be in agreement with the null hypothesis; so
   here we run the restricted regression, saving the coefficient
   vector in bs->b0 and the residuals in bs->u0.
*/

static int do_restricted_ols (boot *bs)
{
    double s2 = 0.0;
    int err = 0;
    
    err = gretl_matrix_restricted_ols(bs->y, bs->X, bs->R, bs->q,
				      bs->b0, NULL, bs->u0, &s2);

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
	bs->SER0 = sqrt(s2);
    }

    return err;
}

static void bs_store_result (boot *bs, gretl_matrix **rp)
{
    if (bs_data != NULL) {
	gretl_matrix_free(bs_data);
	bs_data = NULL;
    }

    bs_data = *rp;
    *rp = NULL;

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
}

static void bs_print_result (boot *bs, gretl_matrix *xi, int tail, PRN *prn)
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
	qsort(xi->val, bs->B, sizeof *xi->val, gretl_compare_doubles);
    }

    if (bs->flags & BOOT_PVAL) {
	pprintf(prn, "%s = %d / %d = %g", _("p-value"), tail, bs->B, bs->pval);
    } else {
	double ql, qu;
	int i, j;

	/* find the alpha/2 quantile */
	i = bs->a * (bs->B + 1) / 2.0;
	ql = xi->val[i-1];

	/* find the 1 - alpha/2 quantile */
	j = bs->B - i + 1;
	qu = xi->val[j-1];

#if BS_DEBUG
	fprintf(stderr, "i=%d, ql=xi[%d] = %g\n", i, i-1, ql);
	fprintf(stderr, "j=%d, qu=xi[%d] = %g\n", j, j-1, qu);
#endif	

	if (studentizing(bs)) {
	    /* the "percentile t" method */
	    double cl = ql;

	    ql = bs->point - bs->sep0 * qu;
	    qu = bs->point - bs->sep0 * cl;
	    pprintf(prn, _("Studentized %g%% confidence interval = %g to %g"), 
		    100 * (1 - bs->a), ql, qu);
	} else {
	    /* the "naive" percentile method */
	    pprintf(prn, _("%g%% confidence interval = %g to %g"), 
		    100 * (1 - bs->a), ql, qu);
	}
    }
    
    pputs(prn, "\n\n");
    pprintf(prn, _("Based on %d replications"), bs->B);
    pputs(prn, ", ");
    if (resampling_u(bs)) {
	pputs(prn, _("using resampled residuals"));
    } else if (resampling_pairs(bs)) {
	pputs(prn, _("using resampled y,X \"pairs\""));
    } else if (wild_boot(bs)) {
	pputs(prn, _("using \"wild\" bootstrap"));
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
	char label[48];

	kdfunc = get_plugin_function("array_kernel_density");
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

	(*kdfunc)(xi->val, bs->B, label);
    }
}

static void bs_calc_ci (boot *bs, gretl_matrix *xi, gretl_matrix *ci)
{
    double ql, qu;
    int i, j;
    
    qsort(xi->val, bs->B, sizeof *xi->val, gretl_compare_doubles);

    /* find the alpha/2 quantile */
    i = bs->a * (bs->B + 1) / 2.0;
    ql = xi->val[i-1];

    /* find the 1 - alpha/2 quantile */
    j = bs->B - i + 1;
    qu = xi->val[j-1];

#if BDEBUG
    fprintf(stderr, "B = %d, ci indices = %d, %d\n", bs->B, i, j);
    if (!studentizing(bs)) {
	fprintf(stderr, "\"basic\" interval: %g to %g\n",
		2*bs->point - qu, 2*bs->point - ql);
    }
#endif    

    if (studentizing(bs)) {
	/* the "percentile t" method */
	ci->val[0] = bs->point - bs->sep0 * qu;
	ci->val[1] = bs->point - bs->sep0 * ql;
    } else if (ci->rows == 2) {
	/* "naive" percentile plus "basic" method */
	gretl_matrix_set(ci, 0, 0, ql);
	gretl_matrix_set(ci, 0, 1, qu);
	gretl_matrix_set(ci, 1, 0, 2*bs->point - qu);
	gretl_matrix_set(ci, 1, 1, 2*bs->point - ql);
    } else {
	/* the "naive" percentile method */
	ci->val[0] = ql;
	ci->val[1] = qu;
    }
}

/* Davidson and MacKinnon, ETM, p. 163; see also
   Davidson and Flachaire, 2001 */

static void rescale_residuals (boot *bs, gretl_matrix *h)
{
    int done = 0;

    if (bs->hc_version == 0) {
	/* no scaling */
	done = 1;
    } else if (h != NULL) {
	int t;

	if (bs->hc_version == 2) {
	    for (t=0; t<bs->T; t++) {
		bs->u0->val[t] /= sqrt(1.0 - h->val[t]);
	    }
	    done = 1;
	} else if (bs->hc_version == 3) {
	    for (t=0; t<bs->T; t++) {
		bs->u0->val[t] /= 1.0 - h->val[t];
	    }
	    done = 1;
	}
    }

    if (!done) {
	/* HC1 or no HC applied */
	int k = bs->k;
	double s;

	if (bs->flags & BOOT_PVAL) {
	    k--;
	}

	s = sqrt((double) bs->T / (bs->T - k));
	gretl_matrix_multiply_by_scalar(bs->u0, s);
    }
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

static void recreate_ldv_X (boot *bs)
{
    int j, p, t;

    for (j=0; j<bs->X->cols; j++) {
	p = ldv_lag(bs, j);
	if (p > 0) {
	    for (t=p; t<bs->T; t++) {
		gretl_matrix_set(bs->X, t, j, bs->y->val[t-p]);
	    }
	}
    }
}

static void print_test_round (boot *bs, int i, double test, PRN *prn)
{
    pprintf(prn, "round %d: test = %g%s\n", i+1, test,
	    test > bs->test0 ? " *" : "");
}

static void fill_hat_vec (gretl_matrix *Q, gretl_matrix *h)
{
    int t, i;
    
    for (t=0; t<Q->rows; t++) {
	double q, ht = 0.0;

	for (i=0; i<Q->cols; i++) {
	    q = gretl_matrix_get(Q, t, i);
	    ht += q * q;
	}
	h->val[t] = ht;
    }
}

static int boot_calc_1 (boot *bs,
			gretl_matrix *XTX,
			gretl_matrix *XTXI,
			gretl_matrix *Q,
			gretl_matrix *R,
			gretl_matrix *h)
{
    int err = 0;

    if (Q != NULL) {
	/* using QR */
	gretl_matrix_copy_values(Q, bs->X);
	err = gretl_matrix_QR_decomp(Q, R);
	if (!err) {
	    err = gretl_invert_triangular_matrix(R, 'U');
	}
	if (!err) {
	    gretl_matrix_multiply_mod(R, GRETL_MOD_NONE,
				      R, GRETL_MOD_TRANSPOSE,
				      XTXI, GRETL_MOD_NONE);
	    if (h != NULL) {
		fill_hat_vec(Q, h);
	    }
	}
    } else {
	/* using Cholesky */
	gretl_matrix_multiply_mod(bs->X, GRETL_MOD_TRANSPOSE,
				  bs->X, GRETL_MOD_NONE,
				  XTX, GRETL_MOD_NONE);
	err = gretl_matrix_cholesky_decomp(XTX);
	if (!err) {
	    err = gretl_inverse_from_cholesky_decomp(XTXI, XTX);
	}	
    }

    return err;
}

static int boot_calc_2 (boot *bs,
			gretl_matrix *XTX,
			gretl_matrix *Q,
			gretl_matrix *R,
			gretl_matrix *g,
			gretl_matrix *b,
			gretl_matrix *yh,
			double *ps2)
{
    int err = 0;

    if (Q != NULL) {
	/* using QR */
	gretl_matrix_multiply_mod(Q, GRETL_MOD_TRANSPOSE,
				  bs->y, GRETL_MOD_NONE, 
				  g, GRETL_MOD_NONE);
	gretl_matrix_multiply(R, g, b);
	gretl_matrix_multiply(Q, g, yh); 
    } else {
	/* using Cholesky */
	gretl_matrix_multiply_mod(bs->X, GRETL_MOD_TRANSPOSE,
				  bs->y, GRETL_MOD_NONE,
				  b, GRETL_MOD_NONE);
	err = gretl_cholesky_solve(XTX, b);
	if (!err) {
	    gretl_matrix_multiply(bs->X, b, yh);
	}
    }

    if (!err) {
	/* residual-related statistics */
	double ut, SSR = 0.0;
	int t;
	
	for (t=0; t<bs->T; t++) {
	    ut = bs->y->val[t] - yh->val[t];
	    if (bs->hc_version >= 0) {
		/* re-use to hold squared residuals */
		yh->val[t] = ut * ut;
	    } else if (boot_use_hac(bs)) {
		/* re-use to hold plain residuals */
		yh->val[t] = ut;
	    } else {
		SSR += ut * ut;
	    }
	}
	if (bs->hc_version < 0 && !boot_use_hac(bs)) {
	    *ps2 = SSR / (bs->T - bs->k);
	}
    }

    return err;
}

static int boot_hac_vcv (boot *bs,
			 const gretl_matrix *XTXI,
			 gretl_matrix *d,
			 gretl_matrix *V)
{
    gretl_matrix *XOX;
    int err = 0;

    XOX = HAC_XOX(bs->X, d, bs->vi, 1, &err);

    if (!err) {
	err = gretl_matrix_qform(XTXI, GRETL_MOD_TRANSPOSE, XOX,
				 V, GRETL_MOD_NONE);
	if (err) {
	    fprintf(stderr, "qform error in boot_hac_vcv\n");
	    abort();
	}
	gretl_matrix_free(XOX);
    }

    return err;
}

static double boot_hac_tau (boot *bs,
			   const gretl_matrix *XTXI,
			   const gretl_matrix *b,
			   gretl_matrix *d,
			   gretl_matrix *V,
			   int *err)
{
    gretl_matrix *XOX;
    double se, tau = NADBL;
    int j = bs->p;

    XOX = HAC_XOX(bs->X, d, bs->vi, 1, err);

    if (!*err) {
	gretl_matrix_qform(XTXI, GRETL_MOD_TRANSPOSE, XOX,
			   V, GRETL_MOD_NONE);
	se = sqrt(gretl_matrix_get(V, j, j));
	tau = (b->val[j] - bs->bp0) / se;
	gretl_matrix_free(XOX);
#if HAC_DEBUG
	fprintf(stderr, "b, se(b), tau = %g, %g, %g\n",
		b->val[j], se, tau);
#endif    
	
    }

    return tau;
}

static double boot_hc_tau (boot *bs,
			   const gretl_matrix *XTXI,
			   const gretl_matrix *b,
			   const gretl_matrix *h,
			   gretl_matrix *d,
			   gretl_matrix *V,
			   int *err)
{
    double se, tau = NADBL;
    int j = bs->p;

    *err = qr_matrix_hccme(bs->X, h, XTXI, d,
			   V, bs->hc_version);
    if (!*err) {
	se = sqrt(gretl_matrix_get(V, j, j));
	tau = (b->val[j] - bs->bp0) / se;
    }
    
    return tau;
}

static double boot_tau (boot *bs,
			const gretl_matrix *XTXI,
			const gretl_matrix *b,
			double s2)
{
    double vpp, se;
    int j = bs->p;

    vpp = gretl_matrix_get(XTXI, j, j);
    se = sqrt(s2 * vpp);
    
    return (b->val[j] - bs->bp0) / se;
}

/* Do the actual bootstrap analysis: the objective is either to form a
   confidence interval or to compute a p-value; the methodology is
   one of

   - resampling the original (scaled) residuals
   - resampling the y, X pairs
   - wild bootstrap (Davidson-Flachaire)
   - simulate normal errors with the empirically given variance
*/

static int real_bootstrap (boot *bs, gretl_matrix *ci, PRN *prn)
{
    gretl_matrix *XTX = NULL;   /* X'X */
    gretl_matrix *XTXI = NULL;  /* X'X^{-1} */
    gretl_matrix *Q = NULL;     /* for use with QR decomp */
    gretl_matrix *R = NULL;     /* for use with QR decomp */
    gretl_matrix *g = NULL;     /* workspace, QR decomp */
    gretl_matrix *h = NULL;     /* "hat" vector (QR) */
    gretl_matrix *d = NULL;     /* workspace */
    gretl_matrix *b = NULL;     /* re-estimated coeffs */
    gretl_matrix *V = NULL;     /* covariance matrix */
    gretl_matrix *r = NULL;     /* recorder for results */
    int *z = NULL;              /* integer resampling array */
    double *xz = NULL;          /* random doubles */
    int k = bs->k;
    int p = bs->p;
    int tail = 0;
    int use_qr = 0;
    int use_h = 0;
    int j, err = 0;

    if ((bs->flags & BOOT_PVAL) && !resampling_pairs(bs)) {
	/* no point in doing this if we're resampling
	   data pairs, since we can't impose H0 
	*/
	err = do_restricted_ols(bs);
	if (err) {
	    return err;
	}
    }

    if (bs->hc_version >= 0 || (bs->flags & BOOT_WILD)) {
	use_qr = use_h = 1;
    }

    b = gretl_column_vector_alloc(k);
    d = gretl_column_vector_alloc(bs->T);
    XTXI = gretl_matrix_alloc(k, k);

    if (b == NULL || d == NULL || XTXI == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (use_qr) {
	Q = gretl_matrix_alloc(bs->T, k);
	R = gretl_matrix_alloc(k, k);
	g = gretl_matrix_alloc(k, 1);
	if (Q == NULL || R == NULL || g == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
	if (use_h) {
	    h = gretl_matrix_alloc(bs->T, 1);
	    if (h == NULL) {
		err = E_ALLOC;
		goto bailout;
	    }
	}
    } else {
	/* Cholesky */
	XTX = gretl_matrix_alloc(k, k);
	if (XTX == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    if (bs->flags & BOOT_WILD_M) {
	/* wild bootstrap with Mammen distribution */
	xz = malloc(bs->T * sizeof *xz);
	if (xz == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    } else if (resampling(bs) || wild_boot(bs)) {
	/* random integer array */
	int nz = bs->T;

	if (bs->blocklen > 1) {
	    nz = bs->T / bs->blocklen + (bs->T % bs->blocklen > 0);
	}
	z = malloc(nz * sizeof *z);
	if (z == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    if (bs->flags & (BOOT_CI | BOOT_GRAPH | BOOT_SAVE)) {
	/* storage for results */
	r = gretl_matrix_alloc(bs->B, 1);
	if (r == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }	

    if (bs->hc_version >= 0 || boot_use_hac(bs) || doing_Ftest(bs)) {
	/* covariance matrix needed */
	V = gretl_matrix_alloc(k, k);
	if (V == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    err = boot_calc_1(bs, XTX, XTXI, Q, R, h);

    if (resampling_u(bs) || wild_boot(bs)) {
	rescale_residuals(bs, h);
    }

    if (!err && verbose(bs)) {
	if (doing_Ftest(bs)) {
	    pputc(prn, '\n');
	} else {
	    pprintf(prn, "%13s %13s\n", "b", "tval");
	}
    }

    /* carry out B replications */

    for (j=0; j<bs->B && !err; j++) {
	double s2 = 0, tau = 0;

#if BDEBUG > 1
	fprintf(stderr, "real_bootstrap: round %d\n", j);
#endif

	if (resampling_u(bs)) {
	    make_resampled_y(bs, z);
	} else if (resampling_pairs(bs)) {
	    make_resampled_pairs(bs, z);
	} else if (wild_boot(bs)) {
	    make_wild_y(bs, z, xz);
	} else {
	    make_normal_y(bs);
	}

	if (bs->ldv != NULL || resampling_pairs(bs)) {
	    /* If the X matrix includes lags of the dependent variable,
	       it has to be rewritten, and X'X-inverse (or Q and R) 
	       recalculated. If we're doing the pairs bootstrap, X will
	       have been revised already but again X'X-inverse or Q, R
	       need redoing.
	    */
	    if (bs->ldv != NULL) {
		recreate_ldv_X(bs);
	    }
	    err = boot_calc_1(bs, XTX, XTXI, Q, R, NULL);
	}

	if (!err) {
	    err = boot_calc_2(bs, XTX, Q, R, g, b, d, &s2);
	}

	if (err) {
	    break;
	}	

	if (doing_Ftest(bs)) {
	    double test = 0;
	    
	    if (bs->hc_version >= 0) {
		err = qr_matrix_hccme(bs->X, h, XTXI, d,
				      V, bs->hc_version);
	    } else if (boot_use_hac(bs)) {
		err = boot_hac_vcv(bs, XTXI, d, V);
	    } else {
		gretl_matrix_copy_values(V, XTXI);
		gretl_matrix_multiply_by_scalar(V, s2);
	    }
	    if (!err) {
		test = bs_F_test(b, V, bs, &err);
		if (verbose(bs)) {
		    print_test_round(bs, j, test, prn);
		}
	    }
	    if (test > bs->test0) {
		tail++;
	    }
	    if (bs->flags & (BOOT_GRAPH | BOOT_SAVE)) {
		r->val[j] = test;
	    }
	    continue;
	}

	if (tau_wanted(bs)) {
	    /* bootstrap t-statistic */
	    if (bs->hc_version >= 0) {
		tau = boot_hc_tau(bs, XTXI, b, h, d, V, &err);
	    } else if (boot_use_hac(bs)) {
		tau = boot_hac_tau(bs, XTXI, b, d, V, &err);
	    } else {
		tau = boot_tau(bs, XTXI, b, s2);
	    }
	    if (verbose(bs)) {
		pprintf(prn, "%13g %13g\n", b->val[p], tau);
	    }
	}

	if (bs->flags & BOOT_CI) {
	    /* doing a confidence interval */
	    if (studentizing(bs)) {
		/* record bootstrap t-stat */
		r->val[j] = tau;
	    } else {
		/* record bootstrap coeff */
		r->val[j] = b->val[p];
	    }
	} else {
	    /* doing p-value */
	    if (bs->flags & (BOOT_GRAPH | BOOT_SAVE)) {
		r->val[j] = tau;
	    }
	    if (fabs(tau) > fabs(bs->test0)) {
		tail++;
	    }
	}
    }

    if (!err) {
	if (ci != NULL) {
	    bs_calc_ci(bs, r, ci);
	} else if (bs->flags & BOOT_PVAL) {
	    bs->pval = (double) tail / bs->B;
	    record_test_result(bs->test0, bs->pval);
	}
	if (!(bs->flags & BOOT_SILENT)) {
	    bs_print_result(bs, r, tail, prn);
	}	
	if (bs->flags & BOOT_SAVE) {
	    bs_store_result(bs, &r);
	}
    }

 bailout:

    gretl_matrix_free(b);
    gretl_matrix_free(XTX);
    gretl_matrix_free(XTXI);
    gretl_matrix_free(Q);
    gretl_matrix_free(R);
    gretl_matrix_free(g);
    gretl_matrix_free(h);
    gretl_matrix_free(d);
    gretl_matrix_free(V);
    gretl_matrix_free(r);

    free(z);
    free(xz);
    
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

    bs->flags |= BOOT_FREE_RQ;
    gretl_vector_set(bs->R, p, 1.0);

    return 0;
}

/**
 * bootstrap_analysis:
 * @pmod: model to be examined.
 * @p: 0-based index number of the coefficient to analyse.
 * @B: number of replications.
 * @alpha: for use when calculating confidence interval.
 * @dset: dataset struct.
 * @opt: option flags: may contain %OPT_P to compute p-value
 * (the default is to calculate confidence interval); %OPT_N
 * to use simulated normal errors, %OPT_X to resample "pairs",
 * or %OPT_W to do the wild bootstrap (the default method 
 * being to resample the empirical residuals); %OPT_G to display 
 * graph, %OPT_S for silent operation.
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
			double alpha, const DATASET *dset,
			gretlopt opt, PRN *prn)
{
    boot *bs = NULL;
    int err = 0;

    if (!bootstrap_ok(pmod->ci)) {
	return E_NOTIMP;
    }

    if (model_sample_problem(pmod, dset)) {
	return E_DATA;
    }      

    if (p < 0 || p >= pmod->ncoeff) {
	return E_DATA;
    }

    bs = boot_new(pmod, dset, B, alpha, opt, &err);

    if (!err && (bs->flags & BOOT_PVAL)) {
	err = bs_add_restriction(bs, p);
    }

    if (!err) {
	int v = pmod->list[p+2];

	bs->p = p;  /* coeff to examine */
	bs->SER0 = pmod->sigma;
	strcpy(bs->vname, dset->varname[v]);
	bs->point = pmod->coeff[p];
	bs->sep0 = pmod->sderr[p];
	bs->test0 = pmod->coeff[p] / pmod->sderr[p];
	if (bs->flags & BOOT_PVAL) {
	    if (opt & OPT_X) {
		bs->bp0 = pmod->coeff[p];
	    } else {
		bs->bp0 = 0.0;
	    }
	} else {
	    bs->bp0 = bs->point;
	}
	err = real_bootstrap(bs, NULL, prn);
    }

    boot_destroy(bs);

    return err;
}

static int opt_from_method (gretlopt *opt, int method)
{
    int err = 0;

    if (method == BOOT_METHOD_PAIRS) {
	/* resample pairs */
	*opt |= OPT_X;
    } else if (method == BOOT_METHOD_WILD) {
	/* wild */
	*opt |= OPT_W;
    } else if (method == BOOT_METHOD_PARAMETRIC) {
	/* normal errors */
	*opt |= OPT_N;
    } else if (method != BOOT_METHOD_RESIDUALS) {
	err = E_DATA;
    }

    return err;
}

#define TRY_BASIC 1

gretl_matrix *bootstrap_ci_matrix (const MODEL *pmod,
				   const DATASET *dset,
				   int p, int B,
				   double alpha,
				   int method,
				   int studentize,
				   int *err)
{
    gretl_matrix *ci = NULL;
    gretlopt opt = OPT_S; /* silent */
    boot *bs = NULL;

    if (!bootstrap_ok(pmod->ci)) {
	*err = E_NOTIMP;
	return NULL;
    }

    if (model_sample_problem(pmod, dset)) {
	*err = E_DATA;
	return NULL;
    }    

    /* convert coefficient index @p to zero-based */
    p -= 1;

    if (p < 0 || p >= pmod->ncoeff) {
	*err = E_DATA;
	return NULL;
    }

    if (!na(alpha) && (alpha < .001 || alpha > .999)) {
	*err = E_DATA;
	return NULL;
    }

    *err = opt_from_method(&opt, method);
    if (*err) {
	return NULL;
    }

    if (studentize) {
	opt |= OPT_T;
	ci = gretl_zero_matrix_new(1, 2);
    } else {
#if TRY_BASIC
	/* calculate c.i. two ways */
	ci = gretl_zero_matrix_new(2, 2);
#else
	ci = gretl_zero_matrix_new(1, 2);
#endif
    }
    
    if (ci == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (na(alpha)) {
	alpha = 0.05;
    }

    bs = boot_new(pmod, dset, B, alpha, opt, err);

    if (method == BOOT_METHOD_WILD) {
	if (libset_get_int(WILDBOOT_DIST) > 0) {
	    bs->flags |= BOOT_WILD_M;
	}
    }

    if (!*err) {
	bs->p = p;  /* coeff to examine */
	bs->SER0 = pmod->sigma;
	bs->point = pmod->coeff[p];
	bs->sep0 = pmod->sderr[p];
	bs->test0 = pmod->coeff[p] / pmod->sderr[p];
	bs->bp0 = bs->point;
	*err = real_bootstrap(bs, ci, NULL);
    }

    if (*err) {
	gretl_matrix_free(ci);
	ci = NULL;
    }

    boot_destroy(bs);

    return ci;
}

double bootstrap_pvalue (const MODEL *pmod,
			 const DATASET *dset,
			 int p, int B,
			 int method,
			 int *err)
{
    double pval = NADBL;
    gretlopt opt = OPT_P | OPT_S; /* p-value, silent */
    boot *bs = NULL;

    if (!bootstrap_ok(pmod->ci)) {
	*err = E_NOTIMP;
	return NADBL;
    }

    if (model_sample_problem(pmod, dset)) {
	*err = E_DATA;
	return NADBL;
    }    

    /* convert coefficient index @p to zero-based */
    p -= 1;

    if (p < 0 || p >= pmod->ncoeff) {
	*err = E_DATA;
	return NADBL;
    }

    *err = opt_from_method(&opt, method);
    if (*err) {
	return NADBL;
    }   

    bs = boot_new(pmod, dset, B, 0.0, opt, err);

    if (!*err) {
	*err = bs_add_restriction(bs, p);
    }

    if (method == BOOT_METHOD_WILD) {
	if (libset_get_int(WILDBOOT_DIST) > 0) {
	    bs->flags |= BOOT_WILD_M;
	}
    }    

    if (!*err) {
	bs->p = p;  /* index of coeff to examine */
	bs->SER0 = pmod->sigma;
	bs->point = pmod->coeff[p];
	bs->sep0 = pmod->sderr[p];
	bs->test0 = pmod->coeff[p] / pmod->sderr[p];
	bs->bp0 = bs->point;
	if (opt & OPT_X) {
	    bs->bp0 = pmod->coeff[p];
	} else {
	    bs->bp0 = 0.0;
	}
	*err = real_bootstrap(bs, NULL, NULL);
    }

    if (!*err) {
	pval = bs->pval;
    }

    boot_destroy(bs);

    return pval;
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
 * @method: bootstrap method.
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
				gretlopt opt, int method,
				PRN *prn)
{
    gretlopt bopt = OPT_P | OPT_R | OPT_F;
    boot *bs = NULL;
    int B = 0;
    int err = 0;

    if (method == BOOT_METHOD_PAIRS) {
	/* This has to be disallowed unless we can come up with
	   a way of testing a relevant H0 that is true under the 
	   bootstrap DGP (easy enough for a single zero
	   restriction, not so easy in general).
	*/
	gretl_errmsg_set("The pairs method cannot be used for this test");
	return E_DATA;
    }

    if (model_sample_problem(pmod, dset)) {
	return E_DATA;
    }

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

    err = opt_from_method(&bopt, method);

    if (!err) {
	gretl_restriction_get_boot_params(&B, &bopt);
	bs = boot_new(pmod, dset, B, 0, bopt, &err);
    }

    if (!err && method == BOOT_METHOD_WILD) {
	if (libset_get_int(WILDBOOT_DIST) > 0) {
	    bs->flags |= BOOT_WILD_M;
	}
    }      

    if (!err) {
	bs->R = R;
	bs->q = q;
	bs->g = g;
	bs->test0 = test;
	strcpy(bs->vname, "F-test");
	err = real_bootstrap(bs, NULL, prn);
    }

    boot_destroy(bs);

    return err;
}

int bootstrap_ok (int ci)
{
    return (ci == OLS || ci == WLS);
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

    err = gretl_matrix_write_as_gdt(fname, bs_data, (const char **) S, 
				    NULL);

    gretl_matrix_free(bs_data);
    bs_data = NULL;
    strings_array_free(S, ns);
    *bs_vname = '\0';

    return err;
}
