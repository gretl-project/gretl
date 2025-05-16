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

/* dpanel.c: implementation of dunamic panel data models of the sort
   developed by Arellano, Bond and Blundell
*/

#include "libgretl.h"
#include "version.h"
#include "matrix_extra.h"
#include "uservar.h"

#define ADEBUG 0
#define WRITE_MATRICES 0
#define DPDEBUG 0
#define IVDEBUG 0

enum {
    DPD_TWOSTEP  = 1 << 0,
    DPD_TIMEDUM  = 1 << 1,
    DPD_WINCORR  = 1 << 2,
    DPD_SYSTEM   = 1 << 3,
    DPD_DPDSTYLE = 1 << 4,
    DPD_REDO     = 1 << 5,
    DPD_COLLAPSE = 1 << 6
};

#define gmm_sys(d)   (d->flags & DPD_SYSTEM)
#define dpd_style(d) (d->flags & DPD_DPDSTYLE)
#define collapse(d)  (d->flags & DPD_COLLAPSE)

#define LEVEL_ONLY 2

typedef struct dpmod_ dpmod;
typedef struct unit_info_ unit_info;
typedef struct diag_info_ diag_info;

struct unit_info_ {
    int t1;      /* first usable obs in differences for unit */
    int t2;      /* last usable obs */
    int nobs;    /* number of usable observations (in the system case
		    this is the sum of the differenced and level
		    observations) */
    int nlev;    /* number of obs in levels (0 in non-system case) */
};

struct diag_info_ {
    int v;        /* ID number of variable */
    int depvar;   /* is the target variable the dependent variable (1/0) */
    int minlag;   /* minimum lag order */
    int maxlag;   /* maximum lag order */
    int level;    /* instrument spec is for levels */
    int rows;     /* max rows occupied in Zi */
    int tbase;    /* first obs with potentially available instruments */
    int collapse; /* "collapse" the instruments? */
};

struct dpmod_ {
    int flags;            /* option flags */
    int step;             /* what step are we on? (1 or 2) */
    int yno;              /* ID number of dependent var */
    int p;                /* lag order for dependent variable */
    int nx;               /* number of independent variables */
    int ifc;              /* model includes a constant */
    int nzr;              /* number of regular instruments */
    int nzb;              /* number of block-diagonal instruments (other than y) */
    int nz;               /* total columns in instrument matrix */
    int pc0;              /* column in Z where predet vars start */
    int xc0;              /* column in Z where exog vars start */
    int N;                /* total number of units in sample */
    int effN;             /* number of units with usable observations */
    int T;                /* total number of observations per unit */
    int minTi;            /* minimum equations (> 0) for any given unit */
    int maxTi;            /* maximum diff. equations for any given unit */
    int max_ni;           /* max number of (possibly stacked) obs per unit
			     (equals maxTi except in system case) */
    int k;                /* number of parameters estimated */
    int nobs;             /* total observations actually used (diffs or levels) */
    int totobs;           /* total observations (diffs + levels) */
    int t1;               /* initial offset into dataset */
    int t1min;            /* first usable observation, any unit */
    int t2max;            /* last usable observation, any unit */
    int ntdum;            /* number of time dummies to use */
    double SSR;           /* sum of squared residuals */
    double s2;            /* residual variance */
    double AR1;           /* z statistic for AR(1) errors */
    double AR2;           /* z statistic for AR(2) errors */
    double sargan;        /* overidentification test statistic (1-step) */
    double hansen;        /* overidentification test statistic (FEGMM) */
    double wald[2];       /* Wald test statistic(s) */
    int wdf[2];           /* degrees of freedom for Wald test(s) */
    int *xlist;           /* list of independent variables */
    int *ilist;           /* list of regular instruments */
    gretl_matrix_block *B1; /* matrix holder */
    gretl_matrix_block *B2; /* matrix holder */
    gretl_matrix *beta;   /* parameter estimates */
    gretl_matrix *vbeta;  /* parameter variance matrix */
    gretl_matrix *uhat;   /* residuals, differenced version */
    gretl_matrix *H;      /* step 1 error covariance matrix */
    gretl_matrix *A;      /* \sum Z'_i H Z_i */
    gretl_matrix *Acpy;   /* back-up of A matrix */
    gretl_matrix *V;      /* covariance matrix */
    gretl_matrix *ZT;     /* transpose of full instrument matrix */
    gretl_matrix *Zi;     /* per-unit instrument matrix */
    gretl_matrix *Y;      /* transformed dependent var */
    gretl_matrix *X;      /* lagged differences of y, indep vars, etc. */
    gretl_matrix *kmtmp;  /* workspace */
    gretl_matrix *kktmp;
    gretl_matrix *M;
    gretl_matrix *L1;
    gretl_matrix *XZA;
    gretl_matrix *ZY;
    gretl_matrix *XZ;
    diag_info *d;          /* info on block-diagonal instruments */
    unit_info *ui;         /* info on panel units */
    char *used;            /* global record of observations used */
    int ndiff;             /* total differenced observations */
    int nlev;              /* total levels observations */
    int nzb2;              /* number of block-diagonal specs, levels eqns */
    int nzdiff;            /* number of insts specific to eqns in differences */
    int nzlev;             /* number of insts specific to eqns in levels */
    int *laglist;          /* (possibly discontinuous) list of lags */
    diag_info *d2;         /* info on block-diagonal instruments, levels eqns
			      (note: not independently allocated) */
    int dcols;             /* number of columns, differenced data */
    int dcolskip;          /* initial number of skipped obs, differences */
    int lcol0;             /* column adjustment for data in levels */
    int lcolskip;          /* initial number of skipped obs, levels */
};

#define data_index(dpd,i) (i * dpd->T + dpd->t1)

static void dpanel_residuals (dpmod *dpd);
static int dpd_process_list (dpmod *dpd, int *list, const int *ylags);

static void dpmod_free (dpmod *dpd)
{
    if (dpd == NULL) {
	return;
    }

    gretl_matrix_block_destroy(dpd->B1);
    gretl_matrix_block_destroy(dpd->B2);

    gretl_matrix_free(dpd->V);

    free(dpd->xlist);
    free(dpd->ilist);
    free(dpd->laglist);

    free(dpd->d);
    free(dpd->used);
    free(dpd->ui);

    free(dpd);
}

static int dpd_allocate_matrices (dpmod *dpd)
{
    int T = dpd->max_ni;

    dpd->B1 = gretl_matrix_block_new(&dpd->beta,  dpd->k, 1,
				     &dpd->vbeta, dpd->k, dpd->k,
				     &dpd->uhat,  dpd->totobs, 1,
				     &dpd->ZT,    dpd->nz, dpd->totobs,
				     &dpd->H,     T, T,
				     &dpd->A,     dpd->nz, dpd->nz,
				     &dpd->Acpy,  dpd->nz, dpd->nz,
				     &dpd->Zi,    T, dpd->nz,
				     &dpd->Y,     dpd->totobs, 1,
				     &dpd->X,     dpd->totobs, dpd->k,
				     NULL);
    if (dpd->B1 == NULL) {
	return E_ALLOC;
    }

    dpd->B2 = gretl_matrix_block_new(&dpd->kmtmp, dpd->k, dpd->nz,
				     &dpd->kktmp, dpd->k, dpd->k,
				     &dpd->M,     dpd->k, dpd->k,
				     &dpd->L1,    1, dpd->nz,
				     &dpd->XZA,   dpd->k, dpd->nz,
				     &dpd->ZY,    dpd->nz, 1,
				     &dpd->XZ,    dpd->k, dpd->nz,
				     NULL);

    if (dpd->B2 == NULL) {
	return E_ALLOC;
    }

    return 0;
}

static int dpd_add_unit_info (dpmod *dpd)
{
    int i, err = 0;

    dpd->ui = malloc(dpd->N * sizeof *dpd->ui);

    if (dpd->ui == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<dpd->N; i++) {
	    dpd->ui[i].t1 = 0;
	    dpd->ui[i].t2 = 0;
	    dpd->ui[i].nobs = 0;
	    dpd->ui[i].nlev = 0;
	}
    }

    return err;
}

static int dpd_flags_from_opt (gretlopt opt)
{
    /* apply Windmeijer correction unless OPT_A is given */
    int f = (opt & OPT_A)? 0 : DPD_WINCORR;

    if (opt & OPT_D) {
	/* include time dummies */
	f |= DPD_TIMEDUM;
    }
    if (opt & OPT_T) {
	/* two-step estimation */
	f |= DPD_TWOSTEP;
    }
    if (opt & OPT_L) {
	/* system GMM: include levels equations */
	f |= DPD_SYSTEM;
    }
    if (opt & OPT_X) {
	/* compute H as per Ox/DPD */
	f |= DPD_DPDSTYLE;
    }
    if (opt & OPT_C) {
	/* "collapse" all block-diagonal instruments as per Roodman */
	f |= DPD_COLLAPSE;
    }

    return f;
}

static dpmod *dpmod_new (int *list, const int *ylags,
			 const DATASET *dset, gretlopt opt,
			 diag_info *d, int nzb, int *err)
{
    dpmod *dpd = NULL;
    int NT;

    if (list[0] < 3) {
	*err = E_PARSE;
	return NULL;
    }

    dpd = malloc(sizeof *dpd);
    if (dpd == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* set pointer members to NULL just in case */
    dpd->B1 = dpd->B2 = NULL;
    dpd->V = NULL;
    dpd->ui = NULL;
    dpd->used = NULL;
    dpd->xlist = NULL;
    dpd->ilist = NULL;
    dpd->laglist = NULL;

    dpd->flags = dpd_flags_from_opt(opt);

    dpd->d = d;
    dpd->nzb = nzb;
    dpd->step = 1;
    dpd->nx = dpd->nzr = 0;
    dpd->nzdiff = 0;
    dpd->nzlev = 0;
    dpd->t1min = dpd->t2max = 0;
    dpd->ntdum = 0;
    dpd->ifc = 0;

    /* "system"-specific */
    dpd->d2 = NULL;
    dpd->nzb2 = 0;

    *err = dpd_process_list(dpd, list, ylags);
    if (*err) {
	goto bailout;
    }

    NT = sample_size(dset);

    dpd->t1 = dset->t1;             /* start of sample range */
    dpd->T = dset->pd;              /* max obs. per individual */
    dpd->effN = dpd->N = NT / dpd->T; /* included individuals */
    dpd->minTi = dpd->maxTi = 0;
    dpd->max_ni = 0;
    dpd->nz = 0;
    dpd->pc0 = 0;
    dpd->xc0 = 0;
    dpd->nobs = dpd->totobs = 0;
    dpd->ndiff = dpd->nlev = 0;
    dpd->SSR = NADBL;
    dpd->s2 = NADBL;
    dpd->AR1 = NADBL;
    dpd->AR2 = NADBL;
    dpd->sargan = NADBL;
    dpd->hansen = NADBL;
    dpd->wald[0] = dpd->wald[1] = NADBL;
    dpd->wdf[0] = dpd->wdf[1] = 0;

    /* # of coeffs on lagged dep var, indep vars */
    if (dpd->laglist != NULL) {
	dpd->k = dpd->laglist[0] + dpd->nx;
    } else {
	dpd->k = dpd->p + dpd->nx;
    }

#if ADEBUG
    fprintf(stderr, "yno = %d, p = %d, nx = %d, k = %d\n",
	    dpd->yno, dpd->p, dpd->nx, dpd->k);
    fprintf(stderr, "t1 = %d, T = %d, N = %d\n", dpd->t1, dpd->T, dpd->N);
#endif

    *err = dpd_add_unit_info(dpd);

    if (!*err) {
	dpd->used = calloc(dset->n, 1);
	if (dpd->used == NULL) {
	    *err = E_ALLOC;
	}
    }

 bailout:

    if (*err) {
	dpmod_free(dpd);
	dpd = NULL;
    }

    return dpd;
}

/* if the const has been included among the regressors but not
   the instruments, add it to the instruments */

static int maybe_add_const_to_ilist (dpmod *dpd)
{
    int i, addc = dpd->ifc;
    int err = 0;

    if (dpd->xlist == NULL || (dpd->nzr == 0 && dpd->nzb == 0)) {
	/* no x's, or all x's treated as exogenous already */
	return 0;
    }

    if (addc && dpd->ilist != NULL) {
	for (i=1; i<=dpd->ilist[0]; i++) {
	    if (dpd->ilist[i] == 0) {
		addc = 0;
		break;
	    }
	}
    }

    if (addc) {
	dpd->ilist = gretl_list_append_term(&dpd->ilist, 0);
	if (dpd->ilist == NULL) {
	    err = E_ALLOC;
	} else {
	    dpd->nzr += 1;
	}
    }

    return err;
}

static int dpd_make_lists (dpmod *dpd, const int *list, int xpos)
{
    int i, nz = 0, spos = 0;
    int err = 0;

    /* do we have a separator, between exog vars and instruments? */
    for (i=xpos; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    spos = i;
	    break;
	}
    }

#if ADEBUG
    printlist(list, "incoming list in dpd_make_lists");
    fprintf(stderr, "separator pos = %d\n", spos);
#endif

    if (spos > 0) {
	dpd->nx = spos - xpos;
	nz = list[0] - spos;
    } else {
	dpd->nx = list[0] - (xpos - 1);
    }

#if ADEBUG
    fprintf(stderr, "got intial dpd->nx = %d, nz = %d\n", dpd->nx, nz);
#endif

    if (dpd->nx > 0) {
	/* compose indep vars list */
	dpd->xlist = gretl_list_new(dpd->nx);
	if (dpd->xlist == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<dpd->nx; i++) {
		dpd->xlist[i+1] = list[xpos + i];
		if (dpd->xlist[i+1] == 0) {
		    dpd->ifc = 1;
		}
	    }
#if ADEBUG
	    printlist(dpd->xlist, "dpd->xlist");
#endif
	    if (nz == 0 && dpd->nzb == 0) {
		/* implicitly all x vars are exogenous */
		dpd->ilist = gretl_list_copy(dpd->xlist);
		if (dpd->ilist == NULL) {
		    err = E_ALLOC;
		} else {
		    dpd->nzr = dpd->ilist[0];
		}
	    }
	}
    }

    if (!err && nz > 0) {
	/* compose regular instruments list */
	dpd->ilist = gretl_list_new(nz);
	if (dpd->ilist == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<nz; i++) {
		dpd->ilist[i+1] = list[spos + i + 1];
	    }
	    dpd->nzr = nz;
	}
    }

    if (!err) {
	err = maybe_add_const_to_ilist(dpd);
    }

#if ADEBUG
    printlist(dpd->ilist, "dpd->ilist");
#endif

    return err;
}

static int dpanel_make_laglist (dpmod *dpd, const int *list,
				int seppos, const int *ylags)
{
    int nlags = seppos - 1;
    int i, err = 0;

    if (nlags != 1) {
	err = E_INVARG;
    } else if (ylags != NULL) {
	dpd->laglist = gretl_list_copy(ylags);
	if (dpd->laglist == NULL) {
	    err = E_ALLOC;
	} else {
	    nlags = dpd->laglist[0];
	}
    } else {
	/* only got p; make "fake" list 1, ..., p */
	dpd->p = list[1];
	if (dpd->p < 1) {
	    err = E_INVARG;
	} else {
	    dpd->laglist = gretl_consecutive_list_new(1, dpd->p);
	    if (dpd->laglist == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (ylags != NULL && !err) {
	/* sort lags and check for duplicates */
	gretl_list_sort(dpd->laglist);
	dpd->p = dpd->laglist[nlags];
	for (i=1; i<=nlags; i++) {
	    if (dpd->laglist[i] <= 0) {
		err = E_INVARG;
	    } else if (i > 1 && dpd->laglist[i] == dpd->laglist[i-1]) {
		err = E_INVARG;
	    }
	}
    }

#if ADEBUG
    printlist(dpd->laglist, "dpd->laglist");
#endif

    return err;
}

/* Process the incoming list and figure out various parameters
   including the maximum lag of the dependent variable, dpd->p,
   and the number of exogenous regressors, dpd->nx.

   The first sublist (before the first LISTSEP) contains either
   a single integer value for p or a set of specific lags of y to
   use as regressors (given in the form of integers in braces or
   as the name of a matrix). In the former case the auxiliary
   @ylags list will be NULL; in the latter @ylags records the
   specific lags requested. At the user level, therefore,

   dpanel 2 ; y ...   means use lags 1 and 2
   dpanel {2} ; y ... means use lag 2 only

   Note that placing a limit on the max lag of y _as instrument_ is
   done via "GMM(y,min,max)".
*/

static int dpd_process_list (dpmod *dpd, int *list,
			     const int *ylags)
{
    int seppos = gretl_list_separator_position(list);
    int xpos, err = 0;

    if (seppos < 2 || list[0] == seppos) {
	/* there must be at least one element before (p) and
	   one element after (y) the first list separator
	*/
	return E_PARSE;
    }

    dpd->yno = list[seppos + 1];
    xpos = seppos + 2;

    if (!dpd_style(dpd) && !gmm_sys(dpd)) {
	/* maybe drop the constant (differenced out) */
	int i;

	for (i=list[0]; i>=xpos; i--) {
	    if (list[i] == 0) {
		gretl_list_delete_at_pos(list, i);
	    }
	}
    }

    /* the auxiliary 'ylags' list may contain specific y lags,
       otherwise there should be just one field */
    err = dpanel_make_laglist(dpd, list, seppos, ylags);

    if (!err) {
	if (dpd->p < 1) {
	    fprintf(stderr, "dpanel lag order = %d < 1\n", dpd->p);
	    err = E_INVARG;
	}
    }

    if (!err && list[0] >= xpos) {
	err = dpd_make_lists(dpd, list, xpos);
    }

    return err;
}

/* Figure the 0-based index of the position of the constant
   in the coefficient vector (or return -1 if the constant is
   not included).
*/

static int dpd_const_pos (dpmod *dpd)
{
    int cpos = -1;

    if (dpd->xlist != NULL) {
	int i;

	for (i=1; i<=dpd->xlist[0]; i++) {
	    if (dpd->xlist[i] == 0) {
		cpos = i - 1;
		break;
	    }
	}

	if (cpos >= 0) {
	    /* we have the position in the array of coeffs on
	       exogenous vars, but these follow the lagged
	       values of the dependent variable.
	    */
	    if (dpd->laglist != NULL) {
		cpos += dpd->laglist[0];
	    } else {
		cpos += dpd->p;
	    }
	}
    }

    return cpos;
}

/* We do either one or two tests here: first we test for the joint
   significance of all regular regressors (lags of y plus any
   exogenous variables, except for the constant, if present).
   Then, if time dummies are present, we test for their joint
   significance.
*/

static int dpd_wald_test (dpmod *dpd)
{
    gretl_matrix *b, *V;
    double x = 0.0;
    int cpos = dpd_const_pos(dpd);
    int k1, knd;
    int i, j, ri, rj;
    int err = 0;

    /* total number of coefficients excluding any time dummies */
    knd = dpd->k - dpd->ntdum;

    /* the number of coefficients to be tested at the first step
       (we exclude the constant)
    */
    k1 = (cpos > 0)? (knd - 1) : knd;

    b = gretl_matrix_reuse(dpd->kmtmp, k1, 1);
    V = gretl_matrix_reuse(dpd->kktmp, k1, k1);

    ri = 0;
    for (i=0; i<knd; i++) {
	if (i != cpos) {
	    b->val[ri++] = dpd->beta->val[i];
	}
    }

    ri = 0;
    for (i=0; i<knd; i++) {
	if (i != cpos) {
	    rj = 0;
	    for (j=0; j<knd; j++) {
		if (j != cpos) {
		    x = gretl_matrix_get(dpd->vbeta, i, j);
		    gretl_matrix_set(V, ri, rj++, x);
		}
	    }
	    ri++;
	}
    }

    err = gretl_invert_symmetric_matrix(V);

    if (!err) {
	x = gretl_scalar_qform(b, V, &err);
    }

    if (!err) {
	dpd->wald[0] = x;
	dpd->wdf[0] = k1;
    }

    if (!err && dpd->ntdum > 0) {
	/* time dummies: these are always at the end of the
	   coeff vector
	*/
	int kt = dpd->ntdum;

	b = gretl_matrix_reuse(dpd->kmtmp, kt, 1);
	V = gretl_matrix_reuse(dpd->kktmp, kt, kt);
	gretl_matrix_extract_matrix(b, dpd->beta, knd, 0,
				    GRETL_MOD_NONE);
	gretl_matrix_extract_matrix(V, dpd->vbeta, knd, knd,
				    GRETL_MOD_NONE);
	err = gretl_invert_symmetric_matrix(V);
	if (!err) {
	    x = gretl_scalar_qform(b, V, &err);
	}
	if (!err) {
	    dpd->wald[1] = x;
	    dpd->wdf[1] = kt;
	}
    }

    gretl_matrix_reuse(dpd->kmtmp, dpd->k, dpd->nz);
    gretl_matrix_reuse(dpd->kktmp, dpd->k, dpd->k);

    if (err) {
	fprintf(stderr, "dpd_wald_test failed: %s\n",
		errmsg_get_with_default(err));
    }

    return err;
}

static int dpd_overid_test (dpmod *dpd)
{
    int L1rows = dpd->L1->rows;
    int L1cols = dpd->L1->cols;
    gretl_matrix *ZTE;
    double test;
    int err = 0;

    ZTE = gretl_matrix_reuse(dpd->L1, dpd->nz, 1);
    gretl_matrix_multiply(dpd->ZT, dpd->uhat, ZTE);

    gretl_matrix_divide_by_scalar(dpd->A, dpd->effN);
    test = gretl_scalar_qform(ZTE, dpd->A, &err);

    gretl_matrix_reuse(dpd->L1, L1rows, L1cols);

    if (!err && test < 0) {
	test = NADBL;
	err = E_NOTPD;
    }

    if (!err && dpd->step == 1) {
	/* allow for scale factor in H matrix */
	double adj = dpd_style(dpd)? 2.0 : 1.0;

	test *= adj / dpd->s2;
    }

    if (err) {
	fprintf(stderr, "dpd_overid_test failed: %s\n",
		errmsg_get_with_default(err));
    } else if (dpd->step == 1 || dpd_style(dpd)) {
	dpd->sargan = test;
    } else {
	dpd->hansen = test;
    }

    return err;
}

/* \sigma^2 H_1, sliced and diced for unit i */

static void make_asy_Hi (dpmod *dpd, int i, gretl_matrix *H,
			 char *mask)
{
    int T = dpd->T;
    int t, s = data_index(dpd, i);

    for (t=0; t<T; t++) {
	mask[t] = (dpd->used[t+s] == 1)? 0 : 1;
    }

    gretl_matrix_reuse(H, T, T);
    gretl_matrix_zero(H);
    gretl_matrix_set(H, 0, 0, 1);

    for (t=1; t<T; t++) {
	gretl_matrix_set(H, t, t, 1);
	gretl_matrix_set(H, t-1, t, -0.5);
	gretl_matrix_set(H, t, t-1, -0.5);
    }

    gretl_matrix_multiply_by_scalar(H, dpd->s2);
    gretl_matrix_cut_rows_cols(H, mask);
}

/*
  AR(1) and AR(2) tests for residuals, see
  Doornik, Arellano and Bond, DPD manual (2006), pp. 28-29.

  AR(m) = \frac{d_0}{\sqrt{d_1 + d_2 + d_3}}

  d_0 = \sum_i w_i' u_i

  d_1 = \sum_i w_i' H_i w_i

  d_2 = -2(\sum_i w_i' X_i) M^{-1}
        (\sum_i X_i' Z_i) A_N (\sum Z_i' H_i w_i)

  d_3 = (\sum_i w_i' X_i) var(\hat{\beta}) (\sum_i X_i' w_i)

  The matrices with subscript i in the above equations are all
  in differences. In the "system" case the last term in d2 is
  modified as

  (\sum Z_i^f' u_i^f u_i' w_i)

  where the f superscript indicates that all observations are
  used, differences and levels.

  one step:         H_i = \sigma^2 H_{1,i}
  one step, robust: H_i = H_{2,i}; M = M_1
  two step:         H_i = H_{2,i}; M = M_2

  H_{2,i} = v^*_{1,i} v^*_{1,i}'

*/

static int dpd_ar_test (dpmod *dpd)
{
    double x, d0, d1, d2, d3;
    gretl_matrix_block *B = NULL;
    gretl_matrix *ui, *wi;
    gretl_matrix *Xi, *Zi;
    gretl_matrix *Hi, *ZU;
    gretl_matrix *wX, *ZHw;
    gretl_matrix *Tmp;
    char *hmask = NULL;
    int HT, T = dpd->maxTi;
    int save_rows, save_cols;
    int nz = dpd->nz;
    int asy, ZU_cols;
    int i, j, k, s, t;
    int nlags, m = 1;
    int err = 0;

    save_rows = gretl_matrix_rows(dpd->Zi);
    save_cols = gretl_matrix_cols(dpd->Zi);

    /* if non-robust and on first step, Hi is computed differently */
    if ((dpd->flags & DPD_TWOSTEP) || (dpd->flags & DPD_WINCORR)) {
	asy = 0;
	HT = T;
    } else {
	asy = 1;
	HT = dpd->T;
	hmask = malloc(HT);
	if (hmask == NULL) {
	    err = E_ALLOC;
	    goto finish;
	}
    }

    /* The required size of the workspace matrix ZU depends on
       whether or not we're in the system case */
    ZU_cols = (gmm_sys(dpd))? 1 : T;

    B = gretl_matrix_block_new(&ui,  T, 1,
			       &wi,  T, 1,
			       &Xi,  T, dpd->k,
			       &Hi,  HT, HT,
			       &ZU,  nz, ZU_cols,
			       &wX,  1, dpd->k,
			       &ZHw, nz, 1,
			       &Tmp, dpd->k, nz,
			       NULL);

    if (B == NULL) {
	err = E_ALLOC;
	goto finish;
    }

    Zi = gretl_matrix_reuse(dpd->Zi, T, nz);

 restart:

    /* initialize cumulators */
    d0 = d1 = 0.0;
    gretl_matrix_zero(wX);
    gretl_matrix_zero(ZHw);
    nlags = 0;

    s = 0;

    for (i=0; i<dpd->N; i++) {
	unit_info *unit = &dpd->ui[i];
	int s_, t0, ni = unit->nobs;
	int Ti = ni - unit->nlev;
	int nlags_i = 0;
	double uw;

	if (Ti == 0) {
	    continue;
	}

	gretl_matrix_reuse(ui, Ti, -1);
	gretl_matrix_reuse(wi, Ti, -1);
	gretl_matrix_reuse(Xi, Ti, -1);
	gretl_matrix_reuse(Zi, Ti, -1);

	if (!asy) {
	    gretl_matrix_reuse(Hi, Ti, Ti);
	}

	if (!gmm_sys(dpd)) {
	    gretl_matrix_reuse(ZU, -1, Ti);
	}

	gretl_matrix_zero(wi);
	t0 = data_index(dpd, i);

	/* Construct the per-unit matrices ui, wi, Xi and Zi. We have
	   to be careful with the dpd->used values for the lags here:
	   1 indicates an observation in (both levels and) differences
	   while observations that are levels-only have a "used" value
	   of LEVEL_ONLY.
	*/

	k = 0;
	for (t=unit->t1; t<=unit->t2; t++) {
	    if (dpd->used[t0+t]) {
		ui->val[k] = dpd->uhat->val[s];
		if (m == 1) {
		    if (dpd->used[t0+t-1] == 1) {
			wi->val[k] = dpd->uhat->val[s-1];
			nlags_i++;
		    }
		} else if (dpd->used[t0+t-2] == 1) {
		    /* m = 2: due to missing values the lag-2
		       residual might be at slot s-1, not s-2
		    */
		    s_ = (dpd->used[t0+t-1] == 1)? (s-2) : (s-1);
		    wi->val[k] = dpd->uhat->val[s_];
		    nlags_i++;
		}
		for (j=0; j<dpd->k; j++) {
		    x = gretl_matrix_get(dpd->X, s, j);
		    gretl_matrix_set(Xi, k, j, x);
		}
		for (j=0; j<nz; j++) {
		    x = gretl_matrix_get(dpd->ZT, j, s);
		    gretl_matrix_set(Zi, k, j, x);
		}
		k++;
		s++;
	    }
	}

	if (nlags_i == 0) {
	    /* no lagged residuals available for this unit */
	    s += unit->nlev;
	    continue;
	}

	nlags += nlags_i;

	/* d0 += w_i' * u_i */
	uw = gretl_matrix_dot_product(wi, GRETL_MOD_TRANSPOSE,
				      ui, GRETL_MOD_NONE, &err);
	d0 += uw;

	if ((dpd->flags & DPD_TWOSTEP) || (dpd->flags & DPD_WINCORR)) {
	    /* form H_i (robust or step 2 variant) */
	    gretl_matrix_multiply_mod(ui, GRETL_MOD_NONE,
				      ui, GRETL_MOD_TRANSPOSE,
				      Hi, GRETL_MOD_NONE);
	} else {
	    make_asy_Hi(dpd, i, Hi, hmask);
	}

	/* d1 += w_i' H_i w_i */
	d1 += gretl_scalar_qform(wi, Hi, &err);

	/* wX += w_i' X_i */
	gretl_matrix_multiply_mod(wi, GRETL_MOD_TRANSPOSE,
				  Xi, GRETL_MOD_NONE,
				  wX, GRETL_MOD_CUMULATE);

	if (gmm_sys(dpd)) {
	    /* Here we need (Z_i^f' u_i^f) * (u_i' w_i), which
	       mixes full series and differences.
	    */
	    gretl_matrix_multiply_mod(Zi, GRETL_MOD_TRANSPOSE,
				      ui, GRETL_MOD_NONE,
				      ZU, GRETL_MOD_NONE);
	    for (t=0; t<unit->nlev; t++) {
		/* catch the levels terms */
		for (j=0; j<nz; j++) {
		    x = gretl_matrix_get(dpd->ZT, j, s);
		    ZU->val[j] += x * dpd->uhat->val[s];
		}
		s++;
	    }
	    gretl_matrix_multiply_by_scalar(ZU, uw);
	    gretl_matrix_add_to(ZHw, ZU);
	} else {
	    /* differences only: ZHw += Z_i' H_i w_i */
	    gretl_matrix_multiply_mod(Zi, GRETL_MOD_TRANSPOSE,
				      Hi, GRETL_MOD_NONE,
				      ZU, GRETL_MOD_NONE);
	    gretl_matrix_multiply_mod(ZU, GRETL_MOD_NONE,
				      wi, GRETL_MOD_NONE,
				      ZHw, GRETL_MOD_CUMULATE);
	}
    }

    if (m == 1) {
	/* form M^{-1} * X'Z * A -- this doesn't have to
	   be repeated for m = 2 */
	gretl_matrix_multiply(dpd->M, dpd->XZ, dpd->kmtmp);
	gretl_matrix_multiply(dpd->kmtmp, dpd->A, Tmp);
    }

    if (nlags == 0) {
	/* no data for AR(m) test at current m */
	if (m == 1) {
	    m = 2;
	    goto restart;
	} else {
	    goto finish;
	}
    }

    /* d2 = -2 * w'X * (M^{-1} * XZ * A) * Z'Hw */
    gretl_matrix_multiply(wX, Tmp, dpd->L1);
    d2 = gretl_matrix_dot_product(dpd->L1, GRETL_MOD_NONE,
				  ZHw, GRETL_MOD_NONE, &err);

    if (!err) {
	d2 *= -2.0;
	/* form w'X * var(\hat{\beta}) * X'w */
	d3 = gretl_scalar_qform(wX, dpd->vbeta, &err);
    }

    if (!err) {
	double z, den = d1 + d2 + d3;

#if ADEBUG
	fprintf(stderr, "ar_test: m=%d, d0=%g, d1=%g, d2=%g, d3=%g, den=%g\n",
		m, d0, d1, d2, d3, den);
#endif
	if (den <= 0) {
	    err = E_NAN;
	} else {
	    z = d0 / sqrt(den);
	    if (m == 1) {
		dpd->AR1 = z;
	    } else {
		dpd->AR2 = z;
	    }
	}
    }

    if (!err && m == 1) {
	m = 2;
	goto restart;
    }

 finish:

    gretl_matrix_block_destroy(B);
    free(hmask);

    /* restore original dimensions */
    gretl_matrix_reuse(dpd->Zi, save_rows, save_cols);

    if (err) {
	fprintf(stderr, "dpd_ar_test failed: %s\n",
		errmsg_get_with_default(err));
    } else if (na(dpd->AR1) && na(dpd->AR1)) {
	fprintf(stderr, "dpd_ar_test: no data\n");
    }

    return err;
}

/* Windmeijer, Journal of Econometrics, 126 (2005), page 33:
   finite-sample correction for step-2 variance matrix.  Note: from a
   computational point of view this calculation is expressed more
   clearly in the working paper by Bond and Windmeijer, "Finite Sample
   Inference for GMM Estimators in Linear Panel Data Models" (2002),
   in particular the fact that the matrix named "D" below must be
   built column by column, taking the derivative of the "W" (or "A")
   matrix with respect to the successive independent variables.
*/

static int windmeijer_correct (dpmod *dpd, const gretl_matrix *uhat1,
			       const gretl_matrix *varb1)
{
    gretl_matrix_block *B;
    gretl_matrix *aV;  /* standard asymptotic variance */
    gretl_matrix *D;   /* finite-sample factor */
    gretl_matrix *dWj; /* one component of the above */
    gretl_matrix *ui;  /* per-unit residuals */
    gretl_matrix *xij; /* per-unit X_j values */
    gretl_matrix *mT;  /* workspace follows */
    gretl_matrix *km;
    gretl_matrix *k1;
    gretl_matrix *R1;
    gretl_matrix *Zui;
    gretl_matrix *Zxi;
    int i, j, t;
    int err = 0;

    aV = gretl_matrix_copy(dpd->vbeta);
    if (aV == NULL) {
	return E_ALLOC;
    }

    B = gretl_matrix_block_new(&D,   dpd->k, dpd->k,
			       &dWj, dpd->nz, dpd->nz,
			       &ui,  dpd->max_ni, 1,
			       &xij, dpd->max_ni, 1,
			       &mT,  dpd->nz, dpd->totobs,
			       &km,  dpd->k, dpd->nz,
			       &k1,  dpd->k, 1,
			       &Zui, dpd->nz, 1,
			       &Zxi, 1, dpd->nz,
			       NULL);
    if (B == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    R1 = dpd->ZY; /* borrowed */

    /* form -(1/N) * asyV * XZW^{-1} */
    gretl_matrix_multiply(aV, dpd->XZA, dpd->kmtmp);
    gretl_matrix_multiply_by_scalar(dpd->kmtmp, -1.0 / dpd->effN);

    /* form W^{-1}Z'v_2 */
    gretl_matrix_multiply(dpd->A, dpd->ZT, mT);
    gretl_matrix_multiply(mT, dpd->uhat, R1);

    for (j=0; j<dpd->k; j++) { /* loop across the X's */
	int s = 0;

	gretl_matrix_zero(dWj);

	/* form dWj = -(1/N) \sum Z_i'(X_j u' + u X_j')Z_i */

	for (i=0; i<dpd->N; i++) {
	    int ni = dpd->ui[i].nobs;

	    if (ni == 0) {
		continue;
	    }

	    gretl_matrix_reuse(ui, ni, 1);
	    gretl_matrix_reuse(xij, ni, 1);
	    gretl_matrix_reuse(dpd->Zi, ni, dpd->nz);

	    /* extract ui (first-step residuals) */
	    for (t=0; t<ni; t++) {
		ui->val[t] = uhat1->val[s++];
	    }

	    /* extract xij */
	    gretl_matrix_extract_matrix(xij, dpd->X, s - ni, j,
					GRETL_MOD_NONE);

	    /* extract Zi */
	    gretl_matrix_extract_matrix(dpd->Zi, dpd->ZT, 0, s - ni,
					GRETL_MOD_TRANSPOSE);

	    gretl_matrix_multiply_mod(dpd->Zi, GRETL_MOD_TRANSPOSE,
				      ui, GRETL_MOD_NONE,
				      Zui, GRETL_MOD_NONE);

	    gretl_matrix_multiply_mod(xij, GRETL_MOD_TRANSPOSE,
				      dpd->Zi, GRETL_MOD_NONE,
				      Zxi, GRETL_MOD_NONE);

	    gretl_matrix_multiply_mod(Zui, GRETL_MOD_NONE,
				      Zxi, GRETL_MOD_NONE,
				      dWj, GRETL_MOD_CUMULATE);
	}

	gretl_matrix_add_self_transpose(dWj);
        gretl_matrix_multiply_by_scalar(dWj, -1.0 / dpd->effN);

        /* D[.,j] = -aV * XZW^{-1} * dWj * W^{-1}Z'v_2 */
        gretl_matrix_multiply(dpd->kmtmp, dWj, km);
        gretl_matrix_multiply(km, R1, k1);

        /* write into appropriate column of D (k x k) */
        for (i=0; i<dpd->k; i++) {
            gretl_matrix_set(D, i, j, k1->val[i]);
	}
    }

    /* add to AsyV: D * AsyV */
    gretl_matrix_multiply_mod(D, GRETL_MOD_NONE,
			      aV, GRETL_MOD_NONE,
			      dpd->vbeta, GRETL_MOD_CUMULATE);

    /* add to AsyV: AsyV * D' */
    gretl_matrix_multiply_mod(aV, GRETL_MOD_NONE,
			      D, GRETL_MOD_TRANSPOSE,
			      dpd->vbeta, GRETL_MOD_CUMULATE);

    /* add to AsyV: D * var(\hat{\beta}_1) * D' */
    gretl_matrix_qform(D, GRETL_MOD_NONE, varb1,
		       dpd->vbeta, GRETL_MOD_CUMULATE);

 bailout:

    gretl_matrix_block_destroy(B);
    gretl_matrix_free(aV);

    return err;
}

/* second-step (asymptotic) variance */

static int dpd_variance_2 (dpmod *dpd,
			   gretl_matrix *u1,
			   gretl_matrix *V1)
{
    int err;

    err = gretl_matrix_qform(dpd->XZ, GRETL_MOD_NONE, dpd->V,
			     dpd->vbeta, GRETL_MOD_NONE);

    if (!err) {
	err = gretl_invert_symmetric_matrix(dpd->vbeta);
    }

    if (!err) {
	gretl_matrix_multiply_by_scalar(dpd->vbeta, dpd->effN);
    }

    if (!err && u1 != NULL && V1 != NULL) {
	err = windmeijer_correct(dpd, u1, V1);
    }

    return err;
}

/*
   Compute the robust step-1 robust variance matrix:

   N * M^{-1} * (X'*Z*A_N*\hat{V}_N*A_N*Z'*X) * M^{-1}

   where

   M = X'*Z*A_N*Z'*X

   and

   \hat{V}_N = N^{-1} \sum Z_i'*v_i*v_i'*Z_i,

   (v_i being the step-1 residuals).

   (But if we're doing 1-step and get the asymptotic
   flag, just do the simple thing, \sigma^2 M^{-1})
*/

static int dpd_variance_1 (dpmod *dpd)
{
    gretl_matrix *V, *ui;
    int i, t, k, c;
    int err = 0;

    if (!(dpd->flags & DPD_TWOSTEP) && !(dpd->flags & DPD_WINCORR)) {
	/* one step asymptotic variance */
	err = gretl_matrix_copy_values(dpd->vbeta, dpd->M);
	gretl_matrix_multiply_by_scalar(dpd->vbeta, dpd->s2 * dpd->effN/2.0);
	return 0;
    }

    V = gretl_zero_matrix_new(dpd->nz, dpd->nz);
    ui = gretl_column_vector_alloc(dpd->max_ni);

    if (V == NULL || ui == NULL) {
	gretl_matrix_free(V);
	gretl_matrix_free(ui);
	return E_ALLOC;
    }

    c = k = 0;

    for (i=0; i<dpd->N; i++) {
	int ni = dpd->ui[i].nobs;

	if (ni == 0) {
	    continue;
	}

	/* get per-unit instruments matrix, Zi */
	gretl_matrix_reuse(dpd->Zi, ni, dpd->nz);
	gretl_matrix_reuse(ui, ni, 1);
	gretl_matrix_extract_matrix(dpd->Zi, dpd->ZT, 0, c,
				    GRETL_MOD_TRANSPOSE);
	c += ni;

	/* load residuals into the ui vector */
	for (t=0; t<ni; t++) {
	    ui->val[t] = dpd->uhat->val[k++];
	}

	gretl_matrix_multiply_mod(ui, GRETL_MOD_TRANSPOSE,
				  dpd->Zi, GRETL_MOD_NONE,
				  dpd->L1, GRETL_MOD_NONE);
	gretl_matrix_multiply_mod(dpd->L1, GRETL_MOD_TRANSPOSE,
				  dpd->L1, GRETL_MOD_NONE,
				  V, GRETL_MOD_CUMULATE);
    }

    gretl_matrix_divide_by_scalar(V, dpd->effN);

    /* form X'Z A_N Z'X  */
    gretl_matrix_multiply(dpd->XZ, dpd->A, dpd->kmtmp);
    gretl_matrix_qform(dpd->kmtmp, GRETL_MOD_NONE, V,
		       dpd->kktmp, GRETL_MOD_NONE);

    if (!err) {
	/* pre- and post-multiply by M^{-1} */
	gretl_matrix_qform(dpd->M, GRETL_MOD_NONE, dpd->kktmp,
			   dpd->vbeta, GRETL_MOD_NONE);
	gretl_matrix_multiply_by_scalar(dpd->vbeta, dpd->effN);
    }

    if (!err && (dpd->flags & DPD_TWOSTEP)) {
	/* preserve V for the second stage */
	dpd->V = V;
    } else {
	gretl_matrix_free(V);
    }

    gretl_matrix_free(ui);

    return err;
}

/* In the "system" case dpd->uhat contains both differenced
   residuals and levels residuals (stacked per unit). In that case
   trim the vector down so that it contains only one set of
   residuals prior to transcribing in series form. Which set
   we preserve is governed by @save_levels.
*/

static int dpanel_adjust_uhat (dpmod *dpd,
			       const DATASET *dset,
			       int save_levels)
{
    double *tmp;
    int ntmp, nd;
    int i, k, s, t;

    ntmp = (save_levels)? dpd->nlev : dpd->ndiff;

    tmp = malloc(ntmp * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    k = s = 0;
    for (i=0; i<dpd->N; i++) {
	nd = dpd->ui[i].nobs - dpd->ui[i].nlev;
	if (save_levels) {
	    /* skip residuals in differences */
	    k += nd;
	    for (t=0; t<dpd->ui[i].nlev; t++) {
		tmp[s++] = dpd->uhat->val[k++];
	    }
	} else {
	    /* use only residuals in differences */
	    for (t=0; t<nd; t++) {
		tmp[s++] = dpd->uhat->val[k++];
	    }
	    k += dpd->ui[i].nlev;
	}
    }

    for (t=0; t<ntmp; t++) {
	dpd->uhat->val[t] = tmp[t];
    }

    gretl_matrix_reuse(dpd->uhat, ntmp, 1);
    free(tmp);

    if (!save_levels) {
	for (t=0; t<dset->n; t++) {
	    if (dpd->used[t] == LEVEL_ONLY) {
		dpd->used[t] = 0;
	    }
	}
    }

    return 0;
}

static void dpd_add_param_names (MODEL *pmod, dpmod *dpd,
				 const DATASET *dset,
				 int full)
{
    char tmp[32];
    int i, j = 0;

    if (dpd->laglist != NULL) {
	for (i=1; i<=dpd->laglist[0]; i++) {
	    sprintf(tmp, "%.10s(-%d)", dset->varname[dpd->yno],
		    dpd->laglist[i]);
	    gretl_model_set_param_name(pmod, j++, tmp);
	}
    } else {
	for (i=0; i<dpd->p; i++) {
	    sprintf(tmp, "%.10s(-%d)", dset->varname[dpd->yno], i+1);
	    gretl_model_set_param_name(pmod, j++, tmp);
	}
    }

    for (i=0; i<dpd->nx; i++) {
	gretl_model_set_param_name(pmod, j++, dset->varname[dpd->xlist[i+1]]);
    }

    if (full) {
	for (i=0; i<dpd->ntdum; i++) {
	    if (dpd->ifc) {
		sprintf(tmp, "T%d", dpd->t1min + i + 2);
		gretl_model_set_param_name(pmod, j++, tmp);
	    } else {
		sprintf(tmp, "T%d", dpd->t1min + i + 1);
		gretl_model_set_param_name(pmod, j++, tmp);
	    }
	}
    }
}

static int dpd_finalize_model (MODEL *pmod,
			       dpmod *dpd,
			       int *list,
			       const int *ylags,
			       const char *istr,
			       const DATASET *dset,
			       gretlopt opt)
{
    const double *y = dset->Z[dpd->yno];
    int keep_extra = opt & OPT_K;
    int i, err = 0;

    if (dpd->flags & DPD_REDO) {
	goto restart;
    }

    pmod->ci = DPANEL;
    pmod->t1 = dset->t1;
    pmod->t2 = dset->t2;
    pmod->dfn = dpd->k;
    pmod->dfd = dpd->nobs - dpd->k;

    pmod->list = list; /* donated, don't free */

    gretl_model_set_int(pmod, "yno", dpd->yno);
    gretl_model_set_int(pmod, "n_included_units", dpd->effN);
    gretl_model_set_int(pmod, "t1min", dpd->t1min + 1);
    gretl_model_set_int(pmod, "t2max", dpd->t2max + 1);
    gretl_model_set_int(pmod, "ntdum", dpd->ntdum);
    gretl_model_set_int(pmod, "ifc", dpd->ifc);

    pmod->ncoeff = dpd->k;
    pmod->nobs = dpd->nobs;
    pmod->full_n = dset->n;

    pmod->rsq = pmod->adjrsq = NADBL;
    pmod->fstt = pmod->chisq = NADBL;
    pmod->lnL = NADBL;

    if (pmod->params == NULL) {
	/* this step not already done */
	gretl_model_allocate_param_names(pmod, dpd->k);
	if (pmod->errcode) {
	    return pmod->errcode;
	}
	dpd_add_param_names(pmod, dpd, dset, 1);
    }

    err = gretl_model_allocate_storage(pmod);

 restart:

    if (!err) {
	pmod->ess = dpd->SSR;
	if (dpd->s2 >= 0) {
	    pmod->sigma = sqrt(dpd->s2);
	}
	for (i=0; i<dpd->k; i++) {
	    pmod->coeff[i] = dpd->beta->val[i];
	}
	err = gretl_model_write_vcv(pmod, dpd->vbeta);
    }

    /* add uhat, yhat */

    if (!err && dpd->nlev > 0) {
	err = dpanel_adjust_uhat(dpd, dset, 1); /* or 0 */
    }

    if (!err) {
	int t, s = 0;

	for (t=0; t<dset->n; t++) {
	    if (dpd->used[t]) {
		pmod->uhat[t] = dpd->uhat->val[s++];
		pmod->yhat[t] = y[t] - pmod->uhat[t];
	    }
	}
    }

    /* additional dpd-specific data */

    if (!err) {
	gretl_model_set_int(pmod, "step", dpd->step);
	if (dpd->step == 2) {
	    pmod->opt |= OPT_T;
	}
	if (!na(dpd->AR1)) {
	    gretl_model_set_double(pmod, "AR1", dpd->AR1);
	}
	if (!na(dpd->AR2)) {
	    gretl_model_set_double(pmod, "AR2", dpd->AR2);
	}
	if (!na(dpd->sargan)) {
	    gretl_model_set_int(pmod, "sargan_df", dpd->nz - dpd->k);
	    gretl_model_set_double(pmod, "sargan", dpd->sargan);
	}
	if (!na(dpd->hansen)) {
	    gretl_model_set_int(pmod, "hansen_df", dpd->nz - dpd->k);
	    gretl_model_set_double(pmod, "hansen", dpd->hansen);
	}
	gretl_model_set_int(pmod, "wald_df", dpd->wdf[0]);
	if (!na(dpd->wald[0])) {
	    gretl_model_set_double(pmod, "wald", dpd->wald[0]);
	}
	if (!na(dpd->wald[1])) {
	    gretl_model_set_int(pmod, "wald_time_df", dpd->wdf[1]);
	    gretl_model_set_double(pmod, "wald_time", dpd->wald[1]);
	}
	if (!(dpd->flags & DPD_WINCORR)) {
	    gretl_model_set_int(pmod, "asy", 1);
	    pmod->opt |= OPT_A;
	}
	if (dpd->flags & DPD_COLLAPSE) {
	    gretl_model_set_int(pmod, "collapse", 1);
	    pmod->opt |= OPT_C;
	}
	if (dpd->A != NULL) {
	    gretl_model_set_int(pmod, "ninst", dpd->A->rows);
	    if (keep_extra) {
		gretl_matrix *A = gretl_matrix_copy(dpd->A);

		gretl_model_set_matrix_as_data(pmod, "wgtmat", A);
	    }
	}
	if (keep_extra && dpd->ZT != NULL) {
	    gretl_matrix *Z = gretl_matrix_copy(dpd->ZT);

	    gretl_model_set_matrix_as_data(pmod, "GMMinst", Z);
	}
	if (opt & OPT_D) {
	    maybe_suppress_time_dummies(pmod, dpd->ntdum);
	}
    }

    if (!err && !(dpd->flags & DPD_REDO)) {
	/* things that only need to be done once */
	if (istr != NULL && *istr != '\0') {
	    gretl_model_set_string_as_data(pmod, "istr", gretl_strdup(istr));
	}
	if (ylags != NULL) {
	    gretl_model_set_list_as_data(pmod, "ylags", gretl_list_copy(ylags));
	}
	if (dpd->flags & DPD_TIMEDUM) {
	    pmod->opt |= OPT_D;
	}
	if (dpd->flags & DPD_SYSTEM) {
	    pmod->opt |= OPT_L;
	}
	if (dpd_style(dpd)) {
	    pmod->opt |= OPT_X;
	}
	if (dpd->maxTi > 0 && dpd->minTi > 0 && dpd->maxTi > dpd->minTi) {
	    gretl_model_set_int(pmod, "Tmin", dpd->minTi);
	    gretl_model_set_int(pmod, "Tmax", dpd->maxTi);
	}
    }

    return err;
}

/* Based on reduction of the A matrix, trim ZT to match and
   adjust the sizes of all workspace matrices that have a
   dimension involving dpd->nz. Note that this is called
   on the first step only, and it is called before computing
   XZ' and Z'Y. The only matrices we need to actually "cut"
   are A (already done when we get here) and ZT; XZ' and Z'Y
   should just have their nz dimension changed.
*/

static void dpd_shrink_matrices (dpmod *dpd, const char *mask)
{
#if IVDEBUG
    fprintf(stderr, "dpanel: dpd_shrink_matrices: cut nz from %d to %d\n",
	    dpd->nz, dpd->A->rows);
#endif

    gretl_matrix_cut_rows(dpd->ZT, mask);
    dpd->nz = dpd->A->rows;

    gretl_matrix_reuse(dpd->Acpy,  dpd->nz, dpd->nz);
    gretl_matrix_reuse(dpd->kmtmp, -1, dpd->nz);
    gretl_matrix_reuse(dpd->L1,    -1, dpd->nz);
    gretl_matrix_reuse(dpd->XZA,   -1, dpd->nz);
    gretl_matrix_reuse(dpd->XZ,    -1, dpd->nz);
    gretl_matrix_reuse(dpd->ZY,    dpd->nz, -1);
}

static int dpd_step_2_A (dpmod *dpd)
{
    int err = 0;

#if ADEBUG
    gretl_matrix_print(dpd->V, "V, in dpd_step_2_A");
#endif

    if (gretl_matrix_rows(dpd->V) > dpd->effN) {
	/* we know this case requires special treatment */
	err = gretl_SVD_invert_matrix(dpd->V);
	if (!err) {
	    gretl_matrix_xtr_symmetric(dpd->V);
	}
    } else {
	gretl_matrix_copy_values(dpd->Acpy, dpd->V);
 	err = gretl_invert_symmetric_matrix(dpd->V);
	if (err) {
	    /* revert the data in dpd->V and try again */
	    gretl_matrix_copy_values(dpd->V, dpd->Acpy);
	    err = gretl_SVD_invert_matrix(dpd->V);
	    if (!err) {
		gretl_matrix_xtr_symmetric(dpd->V);
	    }
	}
    }

    if (!err) {
	/* A <- V^{-1} */
	gretl_matrix_copy_values(dpd->A, dpd->V);
	dpd->step = 2;
    }

    return err;
}

/* compute \hat{\beta} from the moment matrices */

static int dpd_beta_hat (dpmod *dpd)
{
    int err = 0;

    if (dpd->step == 2) {
	/* form the new A matrix */
	err = dpd_step_2_A(dpd);
    }

    if (!err) {
	/* M <- XZ * A * XZ' */
	err = gretl_matrix_qform(dpd->XZ, GRETL_MOD_NONE,
				 dpd->A, dpd->M, GRETL_MOD_NONE);
#if 0
	gretl_matrix_print(dpd->XZ, "XZ");
	gretl_matrix_print(dpd->A, "A");
	gretl_matrix_print(dpd->M, "M");
#endif
    }

    if (!err) {
	/* dpd->XZA <- XZ * A */
	gretl_matrix_multiply(dpd->XZ, dpd->A, dpd->XZA);
	/* using beta as workspace */
	gretl_matrix_multiply(dpd->XZA, dpd->ZY, dpd->beta);
	gretl_matrix_copy_values(dpd->kktmp, dpd->M);
	err = gretl_cholesky_decomp_solve(dpd->kktmp, dpd->beta);
    }

    if (!err) {
	/* prepare M^{-1} for variance calculation */
	err = gretl_inverse_from_cholesky_decomp(dpd->M, dpd->kktmp);
    }

#if ADEBUG > 1
    gretl_matrix_print(dpd->M, "M");
#endif

    return err;
}

/* recompute \hat{\beta} and its variance matrix */

static int dpd_step_2 (dpmod *dpd)
{
    gretl_matrix *u1 = NULL;
    gretl_matrix *V1 = NULL;
    int err;

    dpd->step = 2;

    err = dpd_beta_hat(dpd);

    if (!err && (dpd->flags & DPD_WINCORR)) {
	/* we'll need a copy of the step-1 residuals, and also
	   the step-1 var(\hat{\beta})
	*/
	u1 = gretl_matrix_copy(dpd->uhat);
	V1 = gretl_matrix_copy(dpd->vbeta);
	if (u1 == NULL || V1 == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	dpanel_residuals(dpd);
	err = dpd_variance_2(dpd, u1, V1);
    }

    gretl_matrix_free(u1);
    gretl_matrix_free(V1);

    if (!err) {
	dpd_ar_test(dpd);
	dpd_overid_test(dpd);
	dpd_wald_test(dpd);
    }

    return err;
}

/* This function is used on the first step only */

static int dpd_invert_A_N (dpmod *dpd)
{
    int err = 0;

    /* just in case */
    gretl_matrix_xtr_symmetric(dpd->A);

    /* make a backup in case the first attempt fails */
    gretl_matrix_copy_values(dpd->Acpy, dpd->A);

    /* first try straight inversion */
    err = gretl_invert_symmetric_matrix(dpd->A);

    if (err) {
	/* try again, first reducing A based on its rank */
	char *mask = NULL;

	fprintf(stderr, "inverting dpd->A failed on first pass\n");

	gretl_matrix_copy_values(dpd->A, dpd->Acpy);
	mask = gretl_matrix_rank_mask(dpd->A, &err);

	if (!err) {
	    err = gretl_matrix_cut_rows_cols(dpd->A, mask);
	}

	if (!err) {
	    err = gretl_invert_symmetric_matrix(dpd->A);
	    if (!err) {
		/* OK, now register effects of reducing nz */
		dpd_shrink_matrices(dpd, mask);
	    } else {
		fprintf(stderr, "inverting dpd->A failed on second pass\n");
	    }
	}

	free(mask);
    }

#if ADEBUG
    gretl_matrix_print(dpd->A, "A_N");
#endif

    return err;
}

static int dpd_step_1 (dpmod *dpd, gretlopt opt)
{
    int err;

    /* invert A_N: we allow two attempts */
    err = dpd_invert_A_N(dpd);

    if (!err) {
	/* construct additional moment matrices: we waited
	   until we knew what size these should really be
	*/
	gretl_matrix_multiply(dpd->ZT, dpd->Y, dpd->ZY);
	gretl_matrix_multiply_mod(dpd->X, GRETL_MOD_TRANSPOSE,
				  dpd->ZT, GRETL_MOD_TRANSPOSE,
				  dpd->XZ, GRETL_MOD_NONE);
    }

#if ADEBUG > 1
    gretl_matrix_print(dpd->XZ, "XZ'");
    gretl_matrix_print(dpd->ZY, "Z'Y");
#endif

    if (!err) {
	err = dpd_beta_hat(dpd);
    }

    if (!err) {
	dpanel_residuals(dpd);
	err = dpd_variance_1(dpd);
    }

    if (!err) {
	if (!(dpd->flags & DPD_TWOSTEP)) {
	    /* step 1 is the final destination */
	    if (dpd->nzb2 > 0 && (opt & OPT_A)) {
		/* GMMlev + asymptotic + 1step: as per Ox/DPD,
		   skip the AR tests (which won't work)
		*/
		;
	    } else {
		dpd_ar_test(dpd);
	    }
	    dpd_overid_test(dpd);
	    dpd_wald_test(dpd);
	} else if (!dpd_style(dpd) || (opt & OPT_V)) {
	    /* we're going on to step 2, but record the
	       result of the first-step Sargan test
	    */
	    dpd_overid_test(dpd);
	}
    }

    return err;
}

static int diag_try_list (const char *vname, int *vp, int *nd,
			  diag_info **dp, int **listp)
{
    int *list = get_list_by_name(vname);

    if (list == NULL) {
	return E_UNKVAR;
    } else {
	int nv = list[0];

	if (nv == 1) {
	    *vp = list[1];
	} else if (nv < 1) {
	    return E_DATA;
	} else {
	    /* nv > 1 */
	    int newlen = *nd + nv - 1;
	    diag_info *dnew;

	    dnew = realloc(*dp, newlen * sizeof *dnew);
	    if (dnew == NULL) {
		return E_ALLOC;
	    } else {
		*listp = list;
		*dp = dnew;
		*nd = newlen;
	    }
	}
    }

    return 0;
}

/* Parse a particular entry in the (optional) incoming array of
   "GMM(foo,m1,m2[,collapse])" specifications.  We check that foo
   exists and that m1 and m2 have sane values.
*/

static int parse_diag_info (const char *s, diag_info **dp,
			    int *ip, int *nd, const DATASET *dset)
{
    char vname[VNAMELEN];
    char copt[9];
    char fmt[32];
    int level = 0;
    int i, m1, m2, c;
    int err = 0;

    if (!strncmp(s, "GMM(", 4)) {
	s += 4;
    } else if (!strncmp(s, "GMMlevel(", 9)) {
	level = 1;
	s += 9;
    }
    s += strspn(s, " ");

    /* count the comma separators */
    c = 0;
    for (i=0; s[i] != '\0'; i++) {
	if (s[i] == ',') c++;
    }

    if (c == 3) {
	/* three commas; we need four fields */
	c = 0;
	sprintf(fmt, "%%%d[^, ] , %%d , %%d , %%%d[^ )])", VNAMELEN-1, 8);
	if (sscanf(s, fmt, vname, &m1, &m2, copt) != 4) {
	    err = E_PARSE;
	} else if (!strcmp(copt, "collapse")) {
	    c = 1;
	} else {
	    err = E_PARSE;
	}
    } else if (c == 2) {
	/* we must have three fields */
	c = 0;
	sprintf(fmt, "%%%d[^, ] , %%d , %%d)", VNAMELEN-1);
	if (sscanf(s, fmt, vname, &m1, &m2) != 3) {
	    err = E_PARSE;
	}
    } else {
	err = E_PARSE;
    }

    if (!err) {
	int v = current_series_index(dset, vname);
	int *vlist = NULL;

	if (m1 < 0 || m2 < m1) {
	    err = E_DATA;
	} else if (v < 0) {
	    err = diag_try_list(vname, &v, nd, dp, &vlist);
	}

	if (!err) {
	    diag_info *d, *darray = *dp;
	    int nv = (vlist == NULL)? 1 : vlist[0];
	    int j, i = *ip;

	    for (j=0; j<nv; j++) {
		if (vlist != NULL) {
		    v = vlist[j+1];
		}
		d = &darray[i+j];
		d->depvar = 0;
		d->level = level;
		d->v = v;
		d->minlag = m1;
		d->maxlag = m2;
		d->rows = 0;
		d->collapse = c;
	    }

	    *ip = i + nv;
	}
    }

    return err;
}

/* parse requests of the form

      GMM(xvar, minlag, maxlag[, collapse])

   which call for inclusion of lags of xvar in block-diagonal
   (or perhaps collapsed) fashion
*/

static int
parse_GMM_instrument_spec (const char *spec,
			   const DATASET *dset,
			   diag_info **pd,
			   int *pnspec,
			   int *pnlevel,
			   gretlopt opt)
{
    diag_info *d = NULL;
    const char *s;
    int nspec = 0;
    int err = 0;

    /* first rough check: count closing parentheses */
    s = spec;
    while (*s) {
	if (*s == ')') {
	    nspec++;
	}
	s++;
    }

    if (nspec == 0) {
	/* spec is junk */
	err = E_PARSE;
    } else {
	/* allocate info structs */
	d = calloc(nspec, sizeof *d);
	if (d == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* parse and record individual GMM instrument specs */
	char test[48];
	const char *p;
	int len, i = 0;

	s = spec;
	while (*s && !err) {
	    while (*s == ' ') s++;
	    p = s;
	    while (*p && *p != ')') p++;
	    len = p - s + 1;
	    if (len > 47) {
		err = E_PARSE;
	    } else {
		*test = '\0';
		strncat(test, s, len);
		err = parse_diag_info(test, &d, &i, &nspec, dset);
		s = p + 1;
	    }
	}
    }

    if (!err) {
	int i, j;

	for (i=0; i<nspec && !err; i++) {
	    if (pnlevel != NULL && d[i].level) {
		*pnlevel += 1;
	    }
	    if (opt & OPT_C) {
		/* we got the "global" --collapse option */
		d[i].collapse = 1;
	    }
	    for (j=0; j<i && !err; j++) {
		if (d[i].v == d[j].v && d[i].level == d[j].level) {
		    gretl_errmsg_sprintf(_("variable %d duplicated "
					   "in the command list."),
					 d[i].v);
		    err = E_DATA;
		}
	    }
	}
    }

    if (err) {
	free(d);
	*pnspec = 0;
    } else {
	*pd = d;
	*pnspec = nspec;
    }

    return err;
}

/* Populate the residual vector, dpd->uhat. In the system case
   we stack the residuals in levels under the residuals in
   differences, per unit. Calculate SSR and \sigma^2 while
   we're at it.
*/

static void dpanel_residuals (dpmod *dpd)
{
    const double *b = dpd->beta->val;
    double SSRd = 0.0, SSRl = 0.0;
    double x, ut;
    int i, j, k, t;

    k = 0;

    for (i=0; i<dpd->N; i++) {
	unit_info *unit = &dpd->ui[i];
	int ndiff = unit->nobs - unit->nlev;

	for (t=0; t<ndiff; t++) {
	    /* differences */
	    ut = dpd->Y->val[k];
	    for (j=0; j<dpd->k; j++) {
		x = gretl_matrix_get(dpd->X, k, j);
		ut -= b[j] * x;
	    }
	    SSRd += ut * ut;
	    dpd->uhat->val[k++] = ut;
	}
	for (t=0; t<unit->nlev; t++) {
	    /* levels, if applicable */
	    ut = dpd->Y->val[k];
	    for (j=0; j<dpd->k; j++) {
		x = gretl_matrix_get(dpd->X, k, j);
		ut -= b[j] * x;
	    }
	    SSRl += ut * ut;
	    dpd->uhat->val[k++] = ut;
	}
    }

    if (gmm_sys(dpd)) {
	dpd->nobs = dpd->nlev;
	dpd->SSR = SSRl;
    } else {
	dpd->nobs = dpd->ndiff;
	dpd->SSR = SSRd;
    }

    if (dpd_style(dpd)) {
	dpd->s2 = dpd->SSR / (dpd->nobs - dpd->k);
    } else {
	/* The following produces agreement with xtabond2's Sargan
	   test, though xtabond's reported 'sig2' value is half of
	   the dpd->s2 we calculate here.
	*/
	dpd->SSR = SSRd;
	dpd->s2 = dpd->SSR / (2 * dpd->nobs);
    }
}

/*
   This is the first thing we do for each panel unit: construct a list
   holding the usable observations. A usable observation is defined
   as one for which we have all the data for the equation _in levels_.

   The list of good observations takes this form: the first element
   gives the number of good observations and the following elements
   give their positions, i.e. the 0-based places relative to the start
   of the time-series for the unit. While we're at it we record the
   positions of the good observations relative to the full dataset,
   in the big array dpd->used, so that we can place the residuals
   correctly later on.

   (Note that @t0 gives the offset in the full dataset of the start of
   the unit's data.)

   We return the number of "good observations" minus one, to allow
   for differencing; the return value (if positive) will be the
   max number of observations that are actually usable.
*/

static int check_unit_obs (dpmod *dpd, int *goodobs,
			   const DATASET *dset,
			   int t0)
{
    const double *y = dset->Z[dpd->yno];
    int i, s, t, ok;

    goodobs[0] = 0;

    for (t=0; t<dpd->T; t++) {
	int big_t = t + t0;

	/* do we have the dependent variable? */
	ok = !na(y[big_t]);

	/* and lags of dependent variable? */
	for (i=1; i<=dpd->laglist[0] && ok; i++) {
	    s = t - dpd->laglist[i];
	    if (s < 0) {
		ok = 0;
	    } else {
		ok &= !na(y[s+t0]);
	    }
	}

	if (ok && dpd->xlist != NULL) {
	    /* and regressors? */
	    for (i=1; i<=dpd->xlist[0] && ok; i++) {
		ok &= !na(dset->Z[dpd->xlist[i]][big_t]);
	    }
	}

	if (ok) {
	    goodobs[0] += 1;
	    goodobs[goodobs[0]] = t;
	    if (goodobs[0] > 1) {
		dpd->used[big_t] = 1;
	    } else if (gmm_sys(dpd)) {
		dpd->used[big_t] = LEVEL_ONLY;
	    }
	}
    }

    ok = goodobs[0];

    /* allow for differencing (but don't set ok < 0) */
    if (ok > 0) ok--;

    return ok;
}

static void copy_diag_info (diag_info *targ, diag_info *src)
{
    *targ = *src;
}

/* diff_iv_accounts:

   On input tmin should be the first available obs in levels;
   tmax should be the second-last available obs in levels
   (in each case, for any unit). These indices are based at 0
   for the first period in the unit's data, and they
   represent the subtractive terms in the first and last
   feasible observations in differences, respectively.

   minlag and maxlag represent the minimum and maximum
   lags that have been specified for the given instrument.
   These lags are relative to the "base" of the differenced
   observation, that is, the x_t from which x_{t-k} is
   subtracted to form a difference. With non-gappy data
   k = 1 (the differences are x_t - x_{t-1}) but with
   gappy data we may have k > 1 for some observations;
   nonetheless, we "impute" a difference-base index of
   (the subtractive term's index + 1). This ensures
   that we don't use level instruments that are entangled
   in the difference they are supposed to be instrumenting.
*/

int diff_iv_accounts (dpmod *dpd, int tmin, int tmax)
{
    int t, tbot, ttop;
    int k, i, nrows = 0;

    tbot = tmin + 1;
    ttop = tmax + 1;

#if IVDEBUG
    fprintf(stderr, "*** diff_iv_accounts: tbot = %d, ttop = %d\n", tbot, ttop);
#endif

    for (i=0; i<dpd->nzb; i++) {
	int minlag = dpd->d[i].minlag;
	int maxlag = dpd->d[i].maxlag;
	int usable_maxlag = 0;
	int tbase = tmax + 2;
	int ii, imax = 0;

	dpd->d[i].rows = 0;

#if IVDEBUG
	fprintf(stderr, "GMM spec %d, incoming: minlag = %d, maxlag = %d\n",
		i, minlag, maxlag);
#endif

	/* find tbase = the 'base' of the first differenced observation
	   for which there can be any usable instruments */

	for (t=tbot; t<=ttop; t++) {
	    if (t - minlag >= 0) {
		tbase = t;
		break;
	    }
	}

	if (tbase > ttop) {
	    fprintf(stderr, " no usable instruments for this spec\n");
	    dpd->nzb -= 1;
	    for (k=i; k<dpd->nzb; k++) {
		copy_diag_info(&dpd->d[k], &dpd->d[k+1]);
	    }
	    i--;
	    continue;
	}

#if IVDEBUG
	fprintf(stderr, " tbase = %d\n", tbase);
#endif

	/* step forward, cumulating in-principle usable instruments */
	for (t=tbase; t<=ttop; t++) {
	    ii = 0;
	    for (k=minlag; k<=maxlag && t-k >= 0; k++) {
		ii++;
		if (k > usable_maxlag) {
		    usable_maxlag = k;
		}
	    }
#if IVDEBUG
	    fprintf(stderr, "  max insts at t=%d = %d\n", t, ii);
#endif
	    if (dpd->d[i].collapse) {
		if (ii > imax) imax = ii;
	    } else {
		imax += ii;
	    }
	}

#if IVDEBUG
	fprintf(stderr, " total insts = %d\n", imax);
	fprintf(stderr, " usable maxlag = %d\n", usable_maxlag);
#endif

	dpd->d[i].tbase = tbase;
	dpd->d[i].rows = imax;
	dpd->d[i].maxlag = usable_maxlag;
	nrows += imax;
    }

    return nrows;
}

/* lev_iv_accounts:

   On input tbot should be the first available obs in levels
   and ttop should be the last available obs in levels
   (in each case, for any unit). These indices are based at 0
   for the first period in the unit's data.

   minlag and maxlag represent the minimum and maximum
   lags that have been specified for the given instrument.
*/

int lev_iv_accounts (dpmod *dpd, int tbot, int ttop)
{
    int i, t, k, nrows = 0;

#if IVDEBUG
    fprintf(stderr, "*** lev_iv_accounts: tbot = %d, ttop = %d\n", tbot, ttop);
#endif

    for (i=0; i<dpd->nzb2; i++) {
	int minlag = dpd->d2[i].minlag;
	int maxlag = dpd->d2[i].maxlag;
	int usable_maxlag = 0;
	int tbase = ttop + 1;
	int ii, imax = 0;

	dpd->d2[i].rows = 0;

#if IVDEBUG
	fprintf(stderr, "spec %d: minlag = %d, maxlag = %d\n",
		i, minlag, maxlag);
#endif

	/* find tbase = the first obs in levels for which there can
	   be any usable instruments; since these instruments are
	   differences we need to go back one step beyond minlag
	*/

	for (t=tbot; t<=ttop; t++) {
	    if (t - minlag - 1 >= 0) {
		tbase = t;
		break;
	    }
	}

	if (tbase > ttop) {
	    fprintf(stderr, " no usable instruments for this spec\n");
	    dpd->nzb2 -= 1;
	    for (k=i; k<dpd->nzb2; k++) {
		copy_diag_info(&dpd->d2[k], &dpd->d2[k+1]);
	    }
	    i--;
	    continue;
	}

#if IVDEBUG
	fprintf(stderr, " tbase = %d\n", tbase);
#endif

	/* step forward, cumulating in-principle usable instruments */
	for (t=tbase; t<=ttop; t++) {
	    ii = 0;
	    for (k=minlag; k<=maxlag && t-k-1 >= 0; k++) {
		ii++;
		if (k > usable_maxlag) {
		    usable_maxlag = k;
		}
	    }
#if IVDEBUG
	    fprintf(stderr, "  max insts at t=%d = %d\n", t, ii);
#endif
	    if (dpd->d[i].collapse) {
		if (ii > imax) imax = ii;
	    } else {
		imax += ii;
	    }
	}

#if IVDEBUG
	fprintf(stderr, " total insts = %d\n", imax);
	fprintf(stderr, " usable maxlag = %d\n", usable_maxlag);
#endif

	dpd->d2[i].tbase = tbase;
	dpd->d2[i].rows = imax;
	dpd->d2[i].maxlag = usable_maxlag;
	nrows += imax;
    }

    return nrows;
}

/* Work through the array of block-diagonal instrument
   specifications. Discard any useless ones, trim the
   maxlag values to what is supported on the data,
   and for each spec, count and record the implied number
   of instrument rows that will appear in the Z matrix.
*/

static int block_instrument_count (dpmod *dpd, int t1lev, int t2pen)
{
    int nrows;

    nrows = diff_iv_accounts(dpd, t1lev, t2pen);
    dpd->nzdiff = nrows;
    dpd->nz += dpd->nzdiff;

#if IVDEBUG
    fprintf(stderr, "block_instrument_count, diffs: got %d rows (total = %d)\n",
	    dpd->nzdiff, dpd->nz);
#endif

    if (gmm_sys(dpd)) {
	nrows = lev_iv_accounts(dpd, dpd->t1min, dpd->t2max);
	dpd->nzlev = nrows;
	dpd->nz += dpd->nzlev;
    } else {
	nrows = 0;
    }

#if IVDEBUG
    fprintf(stderr, "block_instrument_count, levels: got %d rows (total = %d)\n",
	    dpd->nzlev, dpd->nz);
#endif

    return 0;
}

/* We do this accounting first so that we know the sizes of the
   various (and numerous) matrices that we'll need for the whole
   analysis at the outset; we can then allocate memory en bloc.
*/

static void do_unit_accounting (dpmod *dpd, const DATASET *dset,
				int **Goodobs)
{
    /* t1lev = index of first good obs in levels,
       t1dif = index of first good obs in differences,
       t2pen = index of penultimate good obs in levels
    */
    int t1lev = dpd->T, t1dif = dpd->T, t2pen = 0;
    int i, t;

    /* just make sure these are zeroed */
    dpd->nzdiff = dpd->nzlev = 0;

    /* total instruments, so far */
    dpd->nz = dpd->nzr;

    /* initialize observation counts */
    dpd->effN = dpd->ndiff = dpd->nlev = 0;
    dpd->minTi = dpd->T;

    /* initialize other accounts */
    dpd->t2max = 0;

    for (i=0, t=dpd->t1; i<dpd->N; i++, t+=dpd->T) {
	int *goodobs = Goodobs[i];
	int Ti = check_unit_obs(dpd, goodobs, dset, t);
	int gmax = goodobs[0];

#if DPDEBUG > 1
	fprintf(stderr, "unit %d: Ti = %d\n", i+1, Ti);
#endif
	if (Ti > 0) {
	    dpd->effN += 1;
	    dpd->ndiff += Ti;
	    if (gmm_sys(dpd)) {
		dpd->nlev += Ti + 1;
	    }
	    if (Ti > dpd->maxTi) {
		dpd->maxTi = Ti;
	    }
	    if (Ti < dpd->minTi) {
		dpd->minTi = Ti;
	    }
	    if (goodobs[1] < t1lev) {
		t1lev = goodobs[1];
	    }
	    if (goodobs[2] < t1dif) {
		t1dif = goodobs[2];
	    }
	    if (goodobs[gmax] > dpd->t2max) {
		dpd->t2max = goodobs[gmax];
	    }
	    if (goodobs[gmax-1] > t2pen) {
		t2pen = goodobs[gmax-1];
	    }
	}
    }

    dpd->t1min = (gmm_sys(dpd))? t1lev : t1dif;

    /* figure number of time dummies, if wanted */
    if (dpd->flags & DPD_TIMEDUM) {
	dpd->ntdum = dpd->t2max - dpd->t1min;
	if (dpd->ifc == 0) {
	    dpd->ntdum += 1;
	}
	dpd->k += dpd->ntdum;
	dpd->nz += dpd->ntdum;
    }

    if (dpd->nzb > 0) {
	/* figure number of extra block-diagonal instruments */
	block_instrument_count(dpd, t1lev, t2pen);
    }

    /* figure the required number of columns for the Yi and Xi data
       matrices: this must be great enough to span the data range
       for all units taken together
    */
    dpd->max_ni = dpd->t2max - dpd->t1min + 1;
    if (dpd->flags & DPD_SYSTEM) {
	dpd->max_ni += dpd->max_ni - 1;
    }

    dpd->dcols = dpd->t2max - t1dif + 1;
    dpd->dcolskip = dpd->p + 1;
    if (t1dif > dpd->dcolskip) {
	dpd->dcolskip = t1dif;
    }
    dpd->lcolskip = (t1lev > dpd->p)? t1lev : dpd->p;
    dpd->lcol0 = dpd->dcols - dpd->lcolskip;

    /* sum the total observations overall */
    dpd->totobs = dpd->ndiff + dpd->nlev;

#if DPDEBUG
    fprintf(stderr, "*** after dpanel accounting:\n"
	    " effN=%d, max_ni=%d, k=%d, ntdum=%d, nz=%d\n",
	    dpd->effN, dpd->max_ni, dpd->k, dpd->ntdum, dpd->nz);
    fprintf(stderr, " maxTi=%d, minTi=%d\n", dpd->maxTi, dpd->minTi);
    fprintf(stderr, " t1min=%d, t2max=%d\n", dpd->t1min, dpd->t2max);
    fprintf(stderr, " t1lev=%d, t1dif=%d\n", t1lev, t1dif);
#endif
}

/* Based on the accounting of good observations for a unit recorded
   in the @goodobs list, fill matrix D (which is used to construct
   H unless we're doing things "dpdstyle").
*/

static void build_unit_D_matrix (dpmod *dpd, int *goodobs, gretl_matrix *D)
{
    int usable = goodobs[0] - 1;
    int i, j, i0, i1;

    gretl_matrix_zero(D);

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	j = i1 - dpd->dcolskip;
	gretl_matrix_set(D, i0, j, -1);
	gretl_matrix_set(D, i1, j,  1);
    }

    /* levels */
    if (gmm_sys(dpd)) {
	for (i=1; i<=goodobs[0]; i++) {
	    i1 = goodobs[i];
	    j = i1 + dpd->lcol0;
	    gretl_matrix_set(D, i1, j, 1);
	}
    }

#if DPDEBUG > 2
    gretl_matrix_print(D, "D");
#endif
}

static void build_unit_H_matrix (dpmod *dpd, int *goodobs,
				 gretl_matrix *D)
{
    build_unit_D_matrix(dpd, goodobs, D);
    gretl_matrix_multiply_mod(D, GRETL_MOD_TRANSPOSE,
			      D, GRETL_MOD_NONE,
			      dpd->H, GRETL_MOD_NONE);
}

static void make_dpdstyle_H (gretl_matrix *H, int nd)
{
    int i;

    gretl_matrix_zero(H);
    gretl_matrix_set(H, 0, 0, 2);

    for (i=1; i<H->rows; i++) {
	if (i < nd) {
	    /* the differences portion */
	    gretl_matrix_set(H, i, i, 2);
	    gretl_matrix_set(H, i-1, i, -1);
	    gretl_matrix_set(H, i, i-1, -1);
	} else {
	    /* the levels portion */
	    gretl_matrix_set(H, i, i, 1);
	}
    }

#if DPDEBUG > 1
    gretl_matrix_print(H, "dpdstyle H");
#endif
}

static int timedum_level (dpmod *dpd, int j, int t)
{
    if (dpd->ifc) {
	return (t == j + 1 + dpd->t1min)? 1 : 0;
    } else {
	return (t == j + dpd->t1min)? 1 : 0;
    }
}

static double timedum_diff (dpmod *dpd, int j, int t)
{
    int d0 = timedum_level(dpd, j, t);
    int d1 = timedum_level(dpd, j, t-1);

    return d0 - d1;
}

/* Build row vector of dependent variable values in @Yi using
   differences, followed by levels if wanted.
*/

static int build_Y (dpmod *dpd, int *goodobs,
		    const DATASET *dset,
		    int t, gretl_matrix *Yi)
{
    const double *y = dset->Z[dpd->yno];
    int i, usable = goodobs[0] - 1;
    int t0, t1, i0, i1, ycol;
    double dy;

    gretl_matrix_zero(Yi);

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	t0 = t + i0;
	t1 = t + i1;
	dy = y[t1] - y[t0];
	ycol = i1 - dpd->dcolskip;
	if (ycol >= Yi->cols) {
	    fprintf(stderr, "Bzzt! scribbling off the end of Yi (diffs)\n"
		    " Yi->cols = %d; i1 - dcolskip = %d - %d = %d\n",
		    Yi->cols, i1, dpd->dcolskip, ycol);
	    return E_DATA;
	} else {
	    gretl_vector_set(Yi, ycol, dy);
	}
    }

    if (gmm_sys(dpd)) {
	/* levels */
	for (i=0; i<=usable; i++) {
	    i1 = goodobs[i+1];
	    t1 = t + i1;
	    ycol = i1 + dpd->lcol0;
	    if (ycol >= Yi->cols) {
		fprintf(stderr, "Bzzt! scribbling off the end of Yi (levels)\n"
			" Yi->cols = %d; i1 + lcol0 = %d + %d = %d\n",
			Yi->cols, i1, dpd->lcol0, ycol);
		fprintf(stderr, " note: dpd->t1min = %d, 1 + dpd->p = %d\n",
			dpd->t1min, 1 + dpd->p);
		return E_DATA;
	    } else {
		gretl_vector_set(Yi, ycol, y[t1]);
	    }
	}
    }

    return 0;
}

/* Build matrix of right-hand side variable values in @Xi using
   differences, followed by levels if wanted.
*/

static void build_X (dpmod *dpd, int *goodobs,
		     const DATASET *dset,
		     int t, gretl_matrix *Xi)
{
    const double *y = dset->Z[dpd->yno];
    int usable = goodobs[0] - 1;
    int nlags = dpd->laglist[0];
    const double *xj;
    int t0, t1, i0, i1;
    int i, j, lj;
    int row, col;
    double dx;

    gretl_matrix_zero(Xi);

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	t0 = t + i0;
	t1 = t + i1;

	row = 0;
	col = i1 - dpd->dcolskip;

	for (j=1; j<=nlags; j++) {
	    lj = dpd->laglist[j];
	    dx = y[t1-lj] - y[t0-lj];
	    gretl_matrix_set(Xi, row++, col, dx);
	}

	for (j=1; j<=dpd->nx; j++) {
	    /* Note: we don't difference away the constant
	       here, but if dpdstyle is not in force the
	       constant will have been removed already.
	    */
	    if (!gmm_sys(dpd) && dpd->xlist[j] == 0) {
		dx = 1.0;
	    } else {
		xj = dset->Z[dpd->xlist[j]];
		dx = xj[t1] - xj[t0];
	    }
	    gretl_matrix_set(Xi, row++, col, dx);
	}

	if (dpd->ntdum > 0) {
	    for (j=0; j<dpd->ntdum; j++) {
		if (gmm_sys(dpd)) {
		    dx = timedum_diff(dpd, j, i1);
		} else if (dpd_style(dpd)) {
		    /* as per DPD: leave dummies in levels */
		    dx = timedum_level(dpd, j, i1);
		} else {
		    /* as per xtabond2: difference the dummies */
		    dx = timedum_diff(dpd, j, i1);
		}
		gretl_matrix_set(Xi, row++, col, dx);
	    }
	}
    }

    if (gmm_sys(dpd)) {
	for (i=0; i<=usable; i++) {
	    i1 = goodobs[i+1];
	    t1 = t + i1;
	    row = 0;
	    col = i1 + dpd->lcol0;

	    for (j=1; j<=nlags; j++) {
		lj = dpd->laglist[j];
		gretl_matrix_set(Xi, row++, col, y[t1-lj]);
	    }

	    for (j=1; j<=dpd->nx; j++) {
		xj = dset->Z[dpd->xlist[j]];
		gretl_matrix_set(Xi, row++, col, xj[t1]);
	    }

	    for (j=0; j<dpd->ntdum; j++) {
		dx = timedum_level(dpd, j, i1);
		gretl_matrix_set(Xi, row++, col, dx);
	    }
	}
    }
}

/* note: this amounts to t1*(t1-1)/2 in the
   straightforward case */

static int row_increment (diag_info *d, int t1)
{
    int k1 = d->level ? 1 : 0;
    int t, k, r = 0;

    for (t=d->tbase; t<t1; t++) {
	for (k=d->minlag; k<=d->maxlag && t-k-k1 >= 0; k++) {
	    r++;
	}
    }

    return r;
}

#if IVDEBUG

/* verify that we're not writing instrument values to rows of
   the Z matrix that are incompatible with what was figured
   out by block_instrument_count (see above).
*/

static int bad_write_check (dpmod *dpd, int row, int lev)
{
    if (!lev && row >= dpd->nzdiff) {
	fprintf(stderr, "*** ERROR in gmm_inst_diff: writing to "
		"bad row %d (max is %d)\n", row,
		dpd->nzdiff - 1);
	return 1;
    } else if (lev && (row < dpd->nzdiff || row >= dpd->nzdiff + dpd->nzlev)) {
	fprintf(stderr, "*** ERROR in gmm_inst_lev: writing to "
		"bad row %d (min is %d, max is %d)\n", row,
		dpd->nzdiff, dpd->nzdiff + dpd->nzlev - 1);
	return 1;
    }

    return 0;
}

#endif

/* GMM-style instruments in levels for the eqns in differences */

static int gmm_inst_diff (dpmod *dpd, int bnum, const double *x,
			  int s, int *goodobs, int row0, int col0,
			  gretl_matrix *Zi)
{
    int maxlag = dpd->d[bnum].maxlag;
    int minlag = dpd->d[bnum].minlag;
    int tmax = goodobs[goodobs[0]];
    int i, t, t1, t2;
    int col, row;
    double xt;

    for (i=1; i<goodobs[0]; i++) {
	t1 = goodobs[i];
	t2 = goodobs[i+1];
	col = col0 + t2 - dpd->dcolskip;
	if (dpd->d[bnum].collapse) {
	    row = row0;
	} else {
	    row = row0 + row_increment(&dpd->d[bnum], t1+1);
	}
	for (t=tmax; t>=0; t--) {
	    /* 2021-11-16: the iteration was: t=0; t<=tmax; t++ */
	    if (t1 - t >= minlag - 1 && t1 - t < maxlag) {
		/* the criterion here needs care! */
		xt = x[s+t];
		if (!na(xt)) {
#if IVDEBUG
		    bad_write_check(dpd, row, 0);
#endif
		    gretl_matrix_set(Zi, row, col, xt);
		}
		row++;
	    }
	}
    }

    return row0 + dpd->d[bnum].rows;
}

/* GMM-style instruments in differences for the eqns in levels */

static int gmm_inst_lev (dpmod *dpd, int bnum, const double *x,
			 int s, int *goodobs, int row0, int col0,
			 gretl_matrix *Zi)
{
    int maxlag = dpd->d2[bnum].maxlag;
    int minlag = dpd->d2[bnum].minlag;
    int tmax = goodobs[goodobs[0]];
    int i, k, t, t1;
    int col, row;
    double x0, x1;

    for (i=1; i<=goodobs[0]; i++) {
	t1 = goodobs[i];
	col = col0 + t1 - dpd->lcolskip;
	if (dpd->d[bnum].collapse) {
	    row = row0;
	} else {
	    row = row0 + row_increment(&dpd->d2[bnum], t1);
	}
	for (t=tmax; t>=1; t--) {
	    /* 2021-11-16: the iteration was: t=1; t<=tmax; t++ */
	    k = t1 - t;
	    if (k <= maxlag && k >= minlag) {
		x0 = x[s+t-1];
		x1 = x[s+t];
		if (!na(x1) && !na(x0)) {
#if IVDEBUG
		    bad_write_check(dpd, row, 1);
#endif
		    gretl_matrix_set(Zi, row, col, x1 - x0);
		}
		row++;
	    }
	}
    }

    return row0 + dpd->d2[bnum].rows;
}

/* Build the matrix of per-unit instrument values in @Zi, which
   has the instruments in rows and the observations in columns.

   Note that each unit's Zi is the same size, padded with zero columns
   for missing observations as needed. The number of columns in Zi
   equals the maximal span of the data for all units taken together,
   counting both observations in differences and observations in
   levels, if applicable.

   We pack the instruments in the following order:

   1) G1: GMM-style instruments in levels for equations in
      differences

   3) G2: GMM-style instruments in differences for equations in
      levels, if present

   5) I1: "Regular" instruments, differenced exog vars for eqns
      in differences

   6) I2: "Regular" instruments, levels of exog vars for eqns
      in levels, if any

   7) D1: Time dummies for eqns in differences, if specified and if
      "system" estimation is not being done

   8) D2: Time dummies for eqns in levels, if specified

   The pattern for the non-system case is

        Z' = | G1 : I1 : D1 |

   and for the full system case it is

        Z' = | G1 :  0 : I1 :  0 |
             |  0 : G2 : I2 : D2 |
*/

static void build_Z (dpmod *dpd, int *goodobs,
		     const DATASET *dset,
		     int t, gretl_matrix *Zi,
		     int unit)
{
    const int usable = goodobs[0] - 1;
    const double *x;
    double dx;
    /* k2 is the starting row for "regular" instruments */
    int k2 = dpd->nzdiff + dpd->nzlev;
    /* k3 marks the starting row for time dummies */
    int k3 = k2 + dpd->nzr;
    int t0, t1, i0, i1;
    int i, j, col, row = 0;

    gretl_matrix_zero(Zi);

#if IVDEBUG
    if (unit == 0) {
	fprintf(stderr, "Z is %d x %d\n", Zi->rows, Zi->cols);
	fprintf(stderr, "  nzb (levels for diffs)  = %d\n", dpd->nzb);
	fprintf(stderr, "  nzb2 (diffs for levels) = %d\n", dpd->nzb2);
	fprintf(stderr, "  nzr (differenced exog vars) = %d\n", dpd->nzr);
	fprintf(stderr, "  ntdum (time dummies for diffs) = %d\n",
		gmm_sys(dpd) ? 0 : dpd->ntdum);
	fprintf(stderr, "  nzr (levels of exog vars) = %d\n",
		!gmm_sys(dpd) ? 0 : dpd->nzr);
	fprintf(stderr, "  ntdum (time dummies for levels) = %d\n",
		!gmm_sys(dpd) ? 0 : dpd->ntdum);
    }
#endif

    /* GMM-style instruments in levels for diffs equations */
    for (i=0; i<dpd->nzb; i++) {
	x = dset->Z[dpd->d[i].v];
	row = gmm_inst_diff(dpd, i, x, t, goodobs, row, 0, Zi);
    }

    col = dpd->t2max - dpd->t1min;

    /* GMM-style instruments in diffs for levels equations */
    for (i=0; i<dpd->nzb2; i++) {
	x = dset->Z[dpd->d2[i].v];
	row = gmm_inst_lev(dpd, i, x, t, goodobs, row, col, Zi);
    }

    /* equations in differences: differenced exog vars */
    if (dpd->nzr > 0) {
	for (i=0; i<usable; i++) {
	    i0 = goodobs[i+1];
	    i1 = goodobs[i+2];
	    col = i1 - dpd->dcolskip;
	    t0 = t + i0;
	    t1 = t + i1;
	    for (j=0; j<dpd->nzr; j++) {
		/* we don't difference the constant, but unless
		   dpdstyle is in force it will have been
		   dropped by this point
		*/
		if (!gmm_sys(dpd) && dpd->ilist[j+1] == 0) {
		    dx = 1.0;
		} else {
		    x = dset->Z[dpd->ilist[j+1]];
		    dx = x[t1] - x[t0];
		}
		gretl_matrix_set(Zi, k2 + j, col, dx);
	    }
	}
    }

    /* equations in differences: time dummies */
    if (dpd->ntdum > 0 && !gmm_sys(dpd)) {
	for (i=0; i<usable; i++) {
	    i1 = goodobs[i+2];
	    col = i1 - dpd->dcolskip;
	    for (j=0; j<dpd->ntdum; j++) {
		dx = timedum_level(dpd, j, i1);
		gretl_matrix_set(Zi, k3 + j, col, dx);
	    }
	}
    }

    if (gmm_sys(dpd)) {
	/* equations in levels: levels of exog vars */
	if (dpd->nzr > 0) {
	    for (i=0; i<=usable; i++) {
		i1 = goodobs[i+1];
		t1 = t + i1;
		col = i1 + dpd->lcol0;
		for (j=0; j<dpd->nzr; j++) {
		    x = dset->Z[dpd->ilist[j+1]];
		    gretl_matrix_set(Zi, k2 + j, col, x[t1]);
		}
	    }
	}

	/* equations in levels: time dummies */
	if (dpd->ntdum > 0) {
	    for (i=0; i<=usable; i++) {
		i1 = goodobs[i+1];
		col = i1 + dpd->lcol0;
		for (j=0; j<dpd->ntdum; j++) {
		    dx = timedum_level(dpd, j, i1);
		    gretl_matrix_set(Zi, k3 + j, col, dx);
		}
	    }
	}
    }
}

static int trim_zero_inst (dpmod *dpd, PRN *prn)
{
    char *mask;
    int err = 0;

#if DPDEBUG
    fprintf(stderr, "before trimming, order of A = %d\n", dpd->A->rows);
#endif

#if WRITE_MATRICES
    gretl_matrix_write_to_file(dpd->A, "dpd-bigA.bin", 0);
#endif

    mask = gretl_matrix_zero_diag_mask(dpd->A, &err);

    if (mask != NULL) {
	err = gretl_matrix_cut_rows_cols(dpd->A, mask);
	if (!err) {
	    if (prn != NULL) {
		pprintf(prn, _("%d redundant instruments dropped, leaving %d\n"),
			dpd->nz - dpd->A->rows, dpd->A->rows);
	    }
	    dpd_shrink_matrices(dpd, mask);
	}
	free(mask);
    }

#if DPDEBUG
    gretl_matrix_print(dpd->A, "dpd->A, after trim_zero_inst");
#endif

    if (!err) {
	gretl_matrix_divide_by_scalar(dpd->A, dpd->effN);
    }

    return err;
}

/* allocate temporary storage needed by do_units() */

static int make_units_workspace (dpmod *dpd, gretl_matrix **D,
				 gretl_matrix **Yi, gretl_matrix **Xi)
{
    int err = 0;

    if (dpd_style(dpd)) {
	/* Ox/DPD-style H matrix: D matrix is not needed */
	*D = NULL;
    } else {
	*D = gretl_matrix_alloc(dpd->T, dpd->max_ni);
	if (*D == NULL) {
	    return E_ALLOC;
	}
    }

    *Yi = gretl_matrix_alloc(1, dpd->max_ni);
    *Xi = gretl_matrix_alloc(dpd->k, dpd->max_ni);

    if (*Yi == NULL || *Xi == NULL) {
	gretl_matrix_free(*D);
	gretl_matrix_free(*Yi);
	gretl_matrix_free(*Xi);
	err = E_ALLOC;
    }

    return err;
}

/* Stack the per-unit data matrices from unit @unum for future use,
   skipping unused observations and recording the numbers of
   observations in differences and in levels.
*/

static void stack_unit_data (dpmod *dpd,
			     const gretl_matrix *Yi,
			     const gretl_matrix *Xi,
			     const gretl_matrix *Zi,
			     int *goodobs, int unum,
			     int *row)
{
    unit_info *unit = &dpd->ui[unum];
    double x;
    int i, j, k, s = *row;

    for (i=2; i<=goodobs[0]; i++) {
	k = goodobs[i] - dpd->dcolskip;
	gretl_vector_set(dpd->Y, s, Yi->val[k]);
	for (j=0; j<Xi->rows; j++) {
	    x = gretl_matrix_get(Xi, j, k);
	    gretl_matrix_set(dpd->X, s, j, x);
	}
	for (j=0; j<dpd->nz; j++) {
	    x = gretl_matrix_get(Zi, j, k);
	    gretl_matrix_set(dpd->ZT, j, s, x);
	}
	s++;
    }

    /* record the indices of the first and last
       differenced observations */
    unit->t1 = goodobs[2];
    unit->t2 = goodobs[goodobs[0]];

    /* record the number of differenced obs */
    unit->nobs = (goodobs[0] > 0)? (goodobs[0] - 1) : 0;

    if (gmm_sys(dpd)) {
	for (i=1; i<=goodobs[0]; i++) {
	    k = goodobs[i] + dpd->lcol0;
	    if (k >= Yi->cols) {
		fprintf(stderr, "*** stack_unit_data: reading off "
			"end of Yi (k=%d, Yi->cols=%d)\n", k, Yi->cols);
		fprintf(stderr, " at goodobs[%d] = %d\n", i, goodobs[i]);
		continue;
	    }
	    gretl_vector_set(dpd->Y, s, Yi->val[k]);
	    for (j=0; j<Xi->rows; j++) {
		x = gretl_matrix_get(Xi, j, k);
		gretl_matrix_set(dpd->X, s, j, x);
	    }
	    for (j=0; j<dpd->nz; j++) {
		x = gretl_matrix_get(Zi, j, k);
		gretl_matrix_set(dpd->ZT, j, s, x);
	    }
	    s++;
	}

	/* record the number of levels obs and augment total */
	unit->nlev = goodobs[0];
	unit->nobs += unit->nlev;
    }

#if WRITE_MATRICES
    gretl_matrix_write_to_file(dpd->ZT, "dpdZT.bin", 0);
#endif

    *row = s;
}

/* Main driver for system GMM: the core is a loop across
   the panel units to build the data and instrument matrices
   and cumulate A = \sum_i Z_i H_i Z_i'.

   At this point we have already done the observations
   accounts, which are recorded in the Goodobs lists.
*/

static int do_units (dpmod *dpd, const DATASET *dset,
		     int **Goodobs)
{
#if DPDEBUG
    char istr[32];
#endif
    gretl_matrix *D = NULL;
    gretl_matrix *Yi = NULL;
    gretl_matrix *Xi = NULL;
    gretl_matrix *Zi = NULL;
    int i, t, Yrow;
    int err = 0;

    err = make_units_workspace(dpd, &D, &Yi, &Xi);
    if (err) {
	return err;
    }

    Zi = dpd->Zi;
    gretl_matrix_reuse(Zi, dpd->nz, dpd->max_ni);

    if (D == NULL) {
	/* the H matrix will not vary by unit */
	int tau = dpd->t2max - dpd->t1min + 1;

	if (gmm_sys(dpd)) {
	    /* t1min is actually "levels-only" */
	    tau--;
	}
	make_dpdstyle_H(dpd->H, tau);
    }

    /* initialize cumulators */
    gretl_matrix_zero(dpd->XZ);
    gretl_matrix_zero(dpd->A);
    gretl_matrix_zero(dpd->ZY);

    /* initialize data stacker */
    Yrow = 0;

#if DPDEBUG
    /* this should not be necessary if stack_unit_data() is
       working correctly */
    gretl_matrix_zero(dpd->Y);
    gretl_matrix_zero(dpd->X);
    gretl_matrix_zero(dpd->ZT);
#endif

    for (i=0; i<dpd->N; i++) {
	int *goodobs = Goodobs[i];
	int Ti = goodobs[0] - 1;

	if (Ti == 0) {
	    continue;
	}

	t = data_index(dpd, i);
	err = build_Y(dpd, goodobs, dset, t, Yi);
	if (err) {
	    break;
	}
	build_X(dpd, goodobs, dset, t, Xi);
	build_Z(dpd, goodobs, dset, t, Zi, i);
#if DPDEBUG
	sprintf(istr, "do_units: Y[%d]", i);
	gretl_matrix_print(Yi, istr);
	sprintf(istr, "do_units: X[%d]", i);
	gretl_matrix_print(Xi, istr);
	sprintf(istr, "do_units: Z[%d]", i);
	gretl_matrix_print(Zi, istr);
#endif
	if (D != NULL) {
	    build_unit_H_matrix(dpd, goodobs, D);
	}
	gretl_matrix_qform(Zi, GRETL_MOD_NONE,
			   dpd->H, dpd->A, GRETL_MOD_CUMULATE);
	/* stack the individual data matrices for future use */
	stack_unit_data(dpd, Yi, Xi, Zi, goodobs, i, &Yrow);
    }

#if DPDEBUG
    gretl_matrix_print(dpd->Y, "dpd->Y");
    gretl_matrix_print(dpd->X, "dpd->X");
    gretl_matrix_print(dpd->H, "dpd->H");
#endif

#if WRITE_MATRICES
    gretl_matrix_write_to_file(dpd->Y, "dpdY.bin", 0);
    gretl_matrix_write_to_file(dpd->X, "dpdX.bin", 0);
#endif

    gretl_matrix_free(D);
    gretl_matrix_free(Yi);
    gretl_matrix_free(Xi);

    return err;
}

/* the user hasn't supplied a block-diagonal spec for
   y in the differences equations: here we set up the
   default version, with unlimited lags
*/

static int add_default_ydiff_spec (dpmod *dpd)
{
    diag_info *d;

    d = realloc(dpd->d, (dpd->nzb + 1) * sizeof *d);

    if (d == NULL) {
	return E_ALLOC;
    } else {
	/* insert the y spec in first place, moving
	   any other specs up */
	int i;

	dpd->d = d;

	for (i=dpd->nzb; i>0; i--) {
	    copy_diag_info(&dpd->d[i], &dpd->d[i-1]);
	}

	d = &dpd->d[0];
	d->v = dpd->yno;
	d->depvar = 1;
	d->minlag = 2;
	d->maxlag = 99;
	d->level = 0;
	d->rows = 0;
	d->collapse = collapse(dpd);

	dpd->nzb += 1;
    }

    return 0;
}

/* the user has specified "system" but hasn't supplied a
   block-diagonal spec for y in the levels equations: here
   we set up the default version, with 1 lag
*/

static int add_default_ylev_spec (dpmod *dpd)
{
    diag_info *d;

    d = realloc(dpd->d, (dpd->nzb + 1) * sizeof *d);

    if (d == NULL) {
	return E_ALLOC;
    } else {
	dpd->d = d;

	d = &dpd->d[dpd->nzb];
	d->v = dpd->yno;
	d->depvar = 1;
	d->minlag = 1;
	d->maxlag = 1;
	d->level = 1;
	d->rows = 0;
	d->collapse = collapse(dpd);

	dpd->nzb += 1;
	dpd->nzb2 += 1;
    }

    return 0;
}

static int compare_gmm_specs (const void *a, const void *b)
{
    const diag_info *da = a;
    const diag_info *db = b;
    int ret = da->level - db->level;

    if (ret == 0) {
	ret = db->depvar - da->depvar;
    }

    return ret;
}

/* Given the info on instrument specification returned by the
   parser, make any adjustments that may be needed in the
   system case.
*/

static int dpanel_adjust_GMM_spec (dpmod *dpd)
{
    int have_ydiff_spec = 0;
    int have_ylev_spec = 0;
    int i, err = 0;

    /* check whether we have:
       - a GMM-style spec for the dep var, differenced equations
       - any block-diagonal specs for levels eqns
       - a GMM-style spec for dep var, levels equations
    */

    for (i=0; i<dpd->nzb; i++) {
	if (dpd->d[i].level == 0) {
	    if (dpd->d[i].v == dpd->yno) {
		dpd->d[i].depvar = 1;
		have_ydiff_spec = 1;
	    }
	} else {
	    dpd->nzb2 += 1;
	    if (dpd->d[i].v == dpd->yno) {
		dpd->d[i].depvar = 1;
		have_ylev_spec = 1;
	    }
	}
    }

    if (!have_ydiff_spec) {
	err = add_default_ydiff_spec(dpd);
	if (err) {
	    return err;
	}
    }

    if (gmm_sys(dpd) && !have_ylev_spec) {
	err = add_default_ylev_spec(dpd);
	if (err) {
	    return err;
	}
    }

    if (dpd->nzb2 > 0 && dpd->nzb2 < dpd->nzb) {
	/* ensure the levels-equations specs come last */
	qsort(dpd->d, dpd->nzb, sizeof *dpd->d, compare_gmm_specs);
    }

    if (dpd->nzb2 > 0) {
	/* henceforth dpd->nzb refers to the number of specs in
	   differences, not the total number */
	dpd->nzb -= dpd->nzb2;
	dpd->d2 = dpd->d + dpd->nzb;
	dpd->flags |= DPD_SYSTEM; /* in case it's not present */
    }

    return err;
}

static char *maxlag_string (char *targ, diag_info *d)
{
    if (d->maxlag == 99) {
	strcpy(targ, "maximum");
    } else {
	sprintf(targ, "%d", d->maxlag);
    }
    return targ;
}

static void print_instrument_specs (dpmod *dpd, const char *ispec,
				    const DATASET *dset, PRN *prn)
{
    char lmax[16];
    int i;

#if IVDEBUG
    if (ispec != NULL) {
	pputs(prn, "Incoming instrument specification:\n");
	pprintf(prn, "  '%s'\n", ispec);
    }
#endif

    pputc(prn, '\n');

    pputs(prn, _("GMM-style instruments, differences equation:\n"));
    for (i=0; i<dpd->nzb; i++) {
	pprintf(prn, _("  %s: %s %d %s %s\n"), dset->varname[dpd->d[i].v],
		_("lags"), dpd->d[i].minlag, _("to"), maxlag_string(lmax, &dpd->d[i]));
    }

    if (dpd->nzb2 > 0) {
	pputs(prn, _("GMM-style instruments, levels equation:\n"));
	for (i=0; i<dpd->nzb2; i++) {
	    pprintf(prn, _("  %s: %s %d %s %s\n"), dset->varname[dpd->d2[i].v],
		    _("lags"), dpd->d2[i].minlag, _("to"), 
		    maxlag_string(lmax, &dpd->d2[i]));
	}
    }
}

/* we're doing two-step, but print the one-step results for
   reference */

static int print_step_1 (dpmod *dpd, MODEL *pmod,
			 const DATASET *dset,
			 PRN *prn)
{
    gretl_array *pnames = NULL;
    gretl_matrix *cse;
    double sei;
    int i, err = 0;

    cse = gretl_matrix_alloc(dpd->k, 2);
    if (cse == NULL) {
	return E_ALLOC;
    }

    gretl_model_allocate_param_names(pmod, dpd->k);
    if (pmod->errcode) {
	return pmod->errcode;
    }
    dpd_add_param_names(pmod, dpd, dset, 1);

    pnames = gretl_array_from_strings(pmod->params, dpd->k,
				      0, &err);

    if (!err) {
	pputc(prn, '\n');
	pprintf(prn, _("Step 1 parameter estimates, using %d observations"),
		dpd->nobs);
	pputc(prn, '\n');
	for (i=0; i<dpd->k; i++) {
	    gretl_matrix_set(cse, i, 0, dpd->beta->val[i]);
	    sei = sqrt(gretl_matrix_get(dpd->vbeta, i, i));
	    gretl_matrix_set(cse, i, 1, sei);
	}
	err = print_model_from_matrices(cse, NULL, pnames, 0,
					OPT_NONE, prn);
	if (!na(dpd->sargan)) {
	    int df = dpd->nz - dpd->k;

	    pputs(prn, "  ");
	    pprintf(prn, _("Sargan test: Chi-square(%d) = %g [%.4f]\n"),
		    df, dpd->sargan, chisq_cdf_comp(df, dpd->sargan));
	}
    }

    gretl_matrix_free(cse);
    if (pnames != NULL) {
	gretl_array_nullify_content(pnames);
	gretl_array_destroy(pnames);
    }

    return err;
}

/* Public interface for the dpanel command */

MODEL dpd_estimate (const int *list, const int *laglist,
		    const char *ispec, const DATASET *dset,
		    gretlopt opt, PRN *prn)
{
    diag_info *d = NULL;
    dpmod *dpd = NULL;
    PRN *vprn = NULL;
    int **Goodobs = NULL;
    int *dlist = NULL;
    int verbose = 0;
    int nlevel = 0;
    int nzb = 0;
    MODEL mod;
    int err = 0;

    gretl_model_init(&mod, dset);

    if (opt & OPT_V) {
	verbose = 1;
	vprn = prn;
    } else if (opt & OPT_Q) {
	prn = NULL;
    }

    /* parse GMM instrument info, if present */
    if (ispec != NULL && *ispec != '\0') {
	mod.errcode = parse_GMM_instrument_spec(ispec, dset, &d,
						&nzb, &nlevel, opt);
	if (mod.errcode) {
	    return mod;
	}
    }

    if (nlevel > 0 && !(opt & OPT_L)) {
	/* make --system (levels) explicit */
	opt |= OPT_L;
    }

    dlist = gretl_list_copy(list);
    if (dlist == NULL) {
	mod.errcode = E_ALLOC;
	return mod;
    }

    dpd = dpmod_new(dlist, laglist, dset, opt, d, nzb, &mod.errcode);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in dpd_init\n", mod.errcode);
	return mod;
    }

    dpanel_adjust_GMM_spec(dpd);

    Goodobs = gretl_list_array_new(dpd->N, dpd->T);
    if (Goodobs == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	do_unit_accounting(dpd, dset, Goodobs);
	err = dpd_allocate_matrices(dpd);
#if ADEBUG
	if (!err) {
	    fprintf(stderr, "allocate_matrices: dpd->Zi is %d x %d\n",
		    dpd->Zi->rows, dpd->Zi->cols);
	}
#endif
    }

    if (!err) {
	/* build the moment matrices */
	err = do_units(dpd, dset, Goodobs);
    }

    gretl_list_array_free(Goodobs, dpd->N);

    if (!err) {
	err = trim_zero_inst(dpd, vprn);
    }

    if (verbose && (dpd->nzb > 0 || dpd->nzb2 > 0)) {
	print_instrument_specs(dpd, ispec, dset, prn);
    }

    if (!err) {
	err = dpd_step_1(dpd, opt);
    }

    if (!err && (opt & OPT_T)) {
	/* second step, if wanted */
	if (verbose) {
	    err = print_step_1(dpd, &mod, dset, prn);
	}
	if (!err) {
	    err = dpd_step_2(dpd);
	}
    }

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	/* write estimation info into model struct */
	mod.errcode = dpd_finalize_model(&mod, dpd, dlist, laglist, ispec,
					 dset, opt);
    } else {
	free(dlist);
    }

    dpmod_free(dpd);

    return mod;
}
