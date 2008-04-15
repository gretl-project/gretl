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

#define ADEBUG 0

typedef struct arbond_ arbond;

struct unit_info {
    int t1;      /* first usable obs for unit */
    int t2;      /* last usable obs */
    char *skip;  /* mask for obs to be skipped, if any */ 
};

struct diag_info {
    int v;       /* ID number of variable */
    int minlag;  /* minimum lag order */
    int maxlag;  /* maximum lag order */
};

struct arbond_ {
    gretlopt opt;         /* option flags */
    int doZX;             /* calculation flag */
    int doZy;             /* calculation flag */
    int step;             /* what step are we on? */
    int yno;              /* ID number of dependent var */
    int p;                /* lag order for dependent variable */
    int qmax;             /* longest lag of y used as instrument */
    int qmin;             /* shortest lag of y used as instrument */
    int nx;               /* number of independent variables */
    int nz;               /* number of regular instruments */
    int nzb;              /* number of block-diagonal instruments (other than y) */
    int m;                /* number of columns in instrument matrix, Z */
    int pc0;              /* column in Z where predet vars start */
    int xc0;              /* column in Z where exog vars start */
    int N;                /* total number of units in sample */
    int effN;             /* number of units with usable observations */
    int T;                /* total number of observations per unit */
    int maxTi;            /* maximum equations for any given unit */
    int k;                /* number of parameters estimated */
    int nobs;             /* total observations used */
    int t1;               /* initial offset into dataset */
    int t1min;            /* first usable observation, any unit */
    int ndum;             /* number of time dummies to use */
    double SSR;           /* sum of squared residuals */
    double s2;            /* residual variance */
    double AR1;           /* z statistic for AR(1) errors */
    double AR2;           /* z statistic for AR(2) errors */
    double sargan;        /* overidentification test statistic */
    double wald;          /* Wald test statistic */
    int wdf;              /* degrees of freedom for Wald test */
    int *xlist;           /* list of independent variables */
    int *ilist;           /* list of instruments */
    gretl_matrix *beta;   /* parameter estimates */
    gretl_matrix *vbeta;  /* parameter variance matrix */
    gretl_matrix *uhat;   /* residuals, differenced version */
    gretl_matrix *H;      /* step 1 error covariance matrix */
    gretl_matrix *A;
    gretl_matrix *Acpy;
    gretl_matrix *V;
    gretl_matrix *ZT;     /* transpose of full instrument matrix */
    gretl_matrix *Zi;     /* per-unit instrument matrix */
    gretl_matrix *dy;     /* first difference of dependent var */
    gretl_matrix *dX;     /* lagged differences of y and indep vars */
    gretl_matrix *tmp1;
    gretl_matrix *kmtmp;
    gretl_matrix *kktmp;
    gretl_matrix *den;
    gretl_matrix *L1;
    gretl_matrix *XZA;
    gretl_matrix *R1;
    gretl_matrix *ZX;
    struct diag_info *d;  /* info on block-diagonal instruments */
    struct unit_info *ui; /* info on panel units */
};

#define data_index(ab,i) (i * ab->T + ab->t1)

static void arbond_free (arbond *ab)
{
    int i;

    gretl_matrix_free(ab->beta);
    gretl_matrix_free(ab->vbeta);
    gretl_matrix_free(ab->uhat);
    gretl_matrix_free(ab->H);
    gretl_matrix_free(ab->A);
    gretl_matrix_free(ab->Acpy);
    gretl_matrix_free(ab->V);
    gretl_matrix_free(ab->Zi);
    gretl_matrix_free(ab->ZT);
    gretl_matrix_free(ab->dy);
    gretl_matrix_free(ab->dX);
    gretl_matrix_free(ab->tmp1);
    gretl_matrix_free(ab->kmtmp);
    gretl_matrix_free(ab->kktmp);
    gretl_matrix_free(ab->den);
    gretl_matrix_free(ab->L1);
    gretl_matrix_free(ab->XZA);
    gretl_matrix_free(ab->R1);
    gretl_matrix_free(ab->ZX);

    free(ab->xlist);
    free(ab->ilist);

    free(ab->d);

    for (i=0; i<ab->N; i++) {
	free(ab->ui[i].skip);
    }
    free(ab->ui);
}

static int arbond_allocate (arbond *ab)
{
    int T2 = ab->maxTi;

    ab->beta = gretl_matrix_alloc(ab->k, 1);
    ab->vbeta = gretl_matrix_alloc(ab->k, ab->k);
    ab->uhat = gretl_matrix_alloc(ab->nobs, 1);
    ab->ZT = gretl_matrix_alloc(ab->m, ab->nobs);
    ab->H = gretl_matrix_alloc(T2, T2);
    ab->A = gretl_matrix_alloc(ab->m, ab->m);
    ab->Acpy = gretl_matrix_alloc(ab->m, ab->m);
    ab->Zi = gretl_matrix_alloc(T2, ab->m);
    ab->dy = gretl_column_vector_alloc(ab->nobs);
    ab->dX = gretl_matrix_alloc(ab->nobs, ab->k);

    ab->tmp1 = gretl_matrix_alloc(ab->m, ab->m);
    ab->kmtmp = gretl_matrix_alloc(ab->k, ab->m);
    ab->kktmp = gretl_matrix_alloc(ab->k, ab->k);
    ab->den = gretl_matrix_alloc(ab->k, ab->k);
    ab->L1 = gretl_matrix_alloc(1, ab->m);
    ab->XZA = gretl_matrix_alloc(ab->k, ab->m);
    ab->R1 = gretl_matrix_alloc(ab->m, 1);
    ab->ZX = gretl_matrix_alloc(ab->m, ab->k);

    if (ab->ZT == NULL || ab->H == NULL || ab->A == NULL ||
	ab->Zi == NULL || ab->dy == NULL || ab->dX == NULL ||
	ab->tmp1 == NULL || ab->kmtmp == NULL ||
	ab->kktmp == NULL || ab->L1 == NULL || ab->XZA == NULL ||
	ab->R1 == NULL || ab->ZX == NULL || ab->den == NULL ||
	ab->beta == NULL || ab->vbeta == NULL ||
	ab->uhat == NULL || ab->Acpy == NULL) {
	return E_ALLOC;
    } else {
	return 0;
    }
}

/* if the const has been included among the regressors but not the
   instruments, add it to the instruments */

static int maybe_add_const_to_ilist (arbond *ab)
{
    int i, addc = 0;
    int err = 0;

    if (ab->xlist == NULL || (ab->nz == 0 && ab->nzb == 0)) {
	/* no x's, or all x's treated as exogenous already */
	return 0;
    }

    for (i=1; i<=ab->xlist[0]; i++) {
	if (ab->xlist[i] == 0) {
	    addc = 1;
	    break;
	}
    }

    if (addc && ab->ilist != NULL) {
	for (i=1; i<=ab->ilist[0]; i++) {
	    if (ab->ilist[i] == 0) {
		addc = 0;
		break;
	    }
	}
    }

    if (addc) {
	ab->ilist = gretl_list_append_term(&ab->ilist, 0);
	if (ab->ilist == NULL) {
	    err = E_ALLOC;
	} else {
	    ab->nz += 1;
	}
    }

    return err;
}

static int arbond_make_lists (arbond *ab, const int *list, int xpos)
{
    int i, nz = 0, spos = 0;
    int err = 0;

    for (i=xpos; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    spos = i;
	    break;
	}
    }

#if ADEBUG
    printlist(list, "incoming list in arbond_make_lists");
    fprintf(stderr, "separator pos = %d\n", spos);
#endif

    if (spos > 0) {
	ab->nx = spos - xpos;
	nz = list[0] - spos;
    } else {
	ab->nx = list[0] - (xpos - 1);
    }

#if ADEBUG
    fprintf(stderr, "got intial ab->nx = %d, nz = %d\n", ab->nx, nz);
#endif

    if (ab->nx > 0) {
	/* compose indep vars list */
	ab->xlist = gretl_list_new(ab->nx);
	if (ab->xlist == NULL) {
	    return E_ALLOC;
	} 
	for (i=0; i<ab->nx; i++) {
	    ab->xlist[i+1] = list[xpos + i];
	}
#if ADEBUG
	printlist(ab->xlist, "ab->xlist");
#endif
	if (nz == 0 && ab->nzb == 0) {
	    /* implicitly all x vars are exogenous */
	    ab->ilist = gretl_list_copy(ab->xlist);
	    if (ab->ilist == NULL) {
		return E_ALLOC;
	    }
	    ab->nz = ab->ilist[0];
	}
    }

    if (nz > 0) {
	/* compose regular instruments list */
	ab->ilist = gretl_list_new(nz);
	if (ab->ilist == NULL) {
	    return E_ALLOC;
	} 
	for (i=0; i<nz; i++) {
	    ab->ilist[i+1] = list[spos + i + 1];
	}
	ab->nz = nz;
    } 

    err = maybe_add_const_to_ilist(ab);

#if ADEBUG
    printlist(ab->ilist, "ab->ilist");
#endif

    return err;
}

static void arbond_zero_matrices (arbond *ab)
{
    ab->beta = NULL;
    ab->vbeta = NULL;
    ab->uhat = NULL;
    ab->ZT = NULL;
    ab->H = NULL;
    ab->A = NULL;
    ab->Acpy = NULL;
    ab->V = NULL;
    ab->Zi = NULL;
    ab->dy = NULL;
    ab->dX = NULL;
    ab->tmp1 = NULL;
    ab->kmtmp = NULL;
    ab->kktmp = NULL;
    ab->den = NULL;
    ab->L1 = NULL;
    ab->XZA = NULL;
    ab->R1 = NULL;
    ab->ZX = NULL;
}

static void arbond_free_lists (arbond *ab)
{
    free(ab->xlist);
    free(ab->ilist);
}

static int 
arbond_init (arbond *ab, const int *list, const DATAINFO *pdinfo,
	     gretlopt opt, struct diag_info *d, int nzb)
{
    int i, xpos = 0;
    int err = 0;

    ab->d = d;
    ab->nzb = nzb;

    if (list[0] < 3) {
	return E_PARSE;
    }

    ab->p = list[1];
    ab->qmin = 2;

    if (list[2] == LISTSEP) {
	ab->qmax = 0;
	ab->yno = list[3];
	xpos = 4;
    } else if (list[3] == LISTSEP && list[0] >= 4) {
	ab->qmax = list[2];
	ab->yno = list[4];
	xpos = 5;
    } else if (list[4] == LISTSEP && list[0] >= 5) {
	ab->qmax = list[2];
	ab->qmin = list[3];
	ab->yno = list[5];
	xpos = 6;
    } else {
	return E_PARSE;
    }

    if (ab->p < 1 || (ab->qmax != 0 && ab->qmax < ab->p + 1) ||
	(ab->qmax != 0 && ab->qmin > ab->qmax)) {
	/* is this all right? */
	return E_INVARG;
    }  

    /* FIXME: qmin doesn't actually do anything below, yet */

    ab->xlist = NULL;
    ab->ilist = NULL;

    ab->opt = opt;
    ab->step = 1;
    ab->doZX = 1;
    ab->doZy = 1;
    ab->nx = 0;
    ab->nz = 0;
    ab->t1min = 0;
    ab->ndum = 0;

    if (list[0] >= xpos) {
	err = arbond_make_lists(ab, list, xpos);
	if (err) {
	    return err;
	}
    }

    ab->t1 = pdinfo->t1;
    ab->T = pdinfo->pd;
    ab->effN = ab->N = (pdinfo->t2 - ab->t1 + 1) / ab->T;
    ab->k = ab->p + ab->nx;
    ab->maxTi = 0;
    ab->m = 0;
    ab->pc0 = 0;
    ab->xc0 = 0;
    ab->nobs = 0;
    ab->SSR = NADBL;
    ab->s2 = NADBL;
    ab->AR1 = NADBL;
    ab->AR2 = NADBL;
    ab->sargan = NADBL;
    ab->wald = NADBL;
    ab->wdf = 0;

#if ADEBUG
    fprintf(stderr, "yno = %d, p = %d, qmax = %d, qmin = %d, nx = %d, k = %d\n",
	    ab->yno, ab->p, ab->qmax, ab->qmin, ab->nx, ab->k);
    fprintf(stderr, "t1 = %d, T = %d, N = %d\n", ab->t1, ab->T, ab->N);
#endif

    arbond_zero_matrices(ab);

    ab->ui = malloc(ab->N * sizeof *ab->ui);
    if (ab->ui == NULL) {
	arbond_free_lists(ab);
	err = E_ALLOC;
    } else {
	for (i=0; i<ab->N; i++) {
	    ab->ui[i].skip = NULL;
	}
    }

    return err;
}

/* See if we have valid values for the dependent variable (in
   differenced or deviations form) plus p lags of same, and all of
   the independent variables, at the given observation, s.
*/

static int obs_is_usable (arbond *ab, const double **Z, int s)
{
    int imax = ab->p + 1;
    int i;

    for (i=0; i<=imax; i++) {
	if (na(Z[ab->yno][s-i])) {
	    return 0;
	}
    }

    if (ab->xlist != NULL) {
	/* check the independent vars */
	for (i=1; i<=ab->xlist[0]; i++) {
	    if (na(Z[ab->xlist[i]][s])) {
		return 0;
	    }
	}
    }

    return 1;
}

static int bzcols (arbond *ab, int i)
{
    int j, k, nc = 0;
    int t = i + ab->p + 1; /* ?? */

    for (j=0; j<ab->nzb; j++) {
	for (k=ab->d[j].minlag; k<=ab->d[j].maxlag; k++) {
	    if (t - k >= 0) {
		nc++;
	    }
	}
    }

    return nc;
}

static int arbond_sample_check (arbond *ab, const double **Z)
{
    const double *y = Z[ab->yno];
    char *mask = NULL;
    int tmin, tmax;
    int t1min = ab->T - 1;
    int t1imin = ab->T - 1;
    int t2max = 0;
    int i, s, t;
    int err = 0;

    for (i=0; i<ab->N; i++) {
	/* find the first y observation, all units */
	if (t1min > 0) {
	    s = data_index(ab, i);
	    for (t=0; t<ab->T; t++) {
		if (!na(y[s+t])) {
		    if (t < t1min) {
			t1min = t;
		    }
		    break;
		}
	    }
	}
    }

#if ADEBUG
    fprintf(stderr, "arbond_sample_check, initial scan: ab->T = %d, t1min = %d\n", 
	    ab->T, t1min);
#endif

    mask = malloc(ab->T);
    if (mask == NULL) {
	return E_ALLOC;
    }

    if (ab->qmax == 0) {
	ab->qmax = ab->T;
    }

    tmin = ab->p + 1;
    tmax = ab->T;

    for (i=0; i<ab->N; i++) {
	int t1i = ab->T - 1, t2i = 0; 
	int Ti = 0, usable = 0;

#if ADEBUG
	fprintf(stderr, "Checking unit %d\n", i);
#endif

	for (t=0; t<ab->T; t++) {
	    mask[t] = 0;
	}

	/* Identify the observations at which we can form the required
	   Delta y terms, have the requisite independent variables,
	   and can construct at least one orthogonality condition
	   using a lagged level of y.
	*/

	s = data_index(ab, i);
	for (t=tmin; t<tmax; t++) {
	    if (obs_is_usable(ab, Z, s + t)) {
		usable++;
		if (t < t1i) t1i = t;
		if (t > t2i) t2i = t;
	    } else {
		mask[t] = 1;
	    }
	}

	if (usable == 0) {
	    ab->effN -= 1;
	    ab->ui[i].t1 = -1;
	    ab->ui[i].t2 = -1;
#if ADEBUG
	    fprintf(stderr, "unit not usable\n");
#endif
	    continue;
	}

	Ti = t2i - t1i + 1;

	if (usable < Ti) {
	    /* there were gaps: steal and replace the skip mask */
	    ab->ui[i].skip = mask;
	    mask = NULL;
	    if (i < ab->N - 1) {
		mask = malloc(ab->T);
		if (mask == NULL) {
		    return E_ALLOC;
		}
	    }
	}

	if (usable > ab->maxTi) {
	    ab->maxTi = usable;
	}
	ab->nobs += usable;
	if (t1i < t1imin) {
	    t1imin = t1i;
	}	
	if (t2i > t2max) {
	    t2max = t2i;
	}
	ab->ui[i].t1 = t1i;
	ab->ui[i].t2 = t2i;
#if ADEBUG
	fprintf(stderr, "t1 = %d, t2 = %d, Ti = %d, usable obs = %d\n", 
		t1i, t2i, Ti, usable);
#endif
    }

    /* record first usable obs, any unit */
    ab->t1min = t1imin;

    /* figure number of time dummies, if wanted */
    if (ab->opt & OPT_D) {
	ab->ndum = t2max - t1imin;
	ab->k += ab->ndum;
    }

    /* is t1min actually "reachable"? */
    if (t1min < t1imin - ab->qmax) {
	t1min = t1imin - ab->qmax;
    }

#if ADEBUG
    fprintf(stderr, "Number of units with usable observations = %d\n", ab->effN);
    fprintf(stderr, "Total usable observations = %d\n", ab->nobs);
    fprintf(stderr, "Max equations for any given unit = %d\n", ab->maxTi);
    fprintf(stderr, "Maximal relevant time-series span: %d to %d = %d\n", 
	    t1min, t2max, t2max - t1min + 1);
#endif

    if (ab->effN == 0) {
	err = E_MISSDATA;
    } else {
	/* compute the number of columns in Zi */
	int tau = t2max - t1min + 1;
	int nblocks = tau - ab->p - 1;
	int bcols, cols = 0;

#if ADEBUG
	fprintf(stderr, "\ntau = %d (ab->p = %d)\n", tau, ab->p);
#endif
	for (i=0; i<nblocks; i++) {
	    /* lagged y values */
	    cols = (ab->p + i > ab->qmax - 1)? ab->qmax - 1 : ab->p + i;
	    ab->m += cols;
#if ADEBUG
	    fprintf(stderr, "i=%d: adding %d cols for y-lags\n", i, cols);
#endif
	    if (ab->nzb > 0) {
		/* other block-diagonal instruments */
		bcols = bzcols(ab, i);
		ab->m += bcols;
#if ADEBUG
		fprintf(stderr, " plus %d cols for z-lags\n", bcols);
#endif
	    }
	}
#if ADEBUG
	fprintf(stderr, "'basic' m = %d\n", ab->m);
#endif
	ab->qmax = cols + 1;
	/* record the column where the exogenous vars start */
	ab->xc0 = ab->m;
	ab->m += ab->nz;
	ab->m += ab->ndum;
#if ADEBUG
	fprintf(stderr, "total m = %d (dummies=%d, exog=%d)\n", 
		ab->m, ab->ndum, ab->nz);
#endif
    }

    free(mask);

    return err;
}

static int unit_nobs (arbond *ab, int i)
{
    int n;

    if (ab->ui[i].t1 < 0) return 0;

    n = ab->ui[i].t2 - ab->ui[i].t1 + 1;

    if (ab->ui[i].skip != NULL) {
	int t;

	for (t=ab->ui[i].t1; t<=ab->ui[i].t2; t++) {
	    if (ab->ui[i].skip[t]) {
		n--;
	    }
	}
    }

    return n;
}

static int skip_obs (arbond *ab, int i, int t)
{
    return (ab->ui[i].skip == NULL)? 0 : ab->ui[i].skip[t];
}

static int skip_unit (arbond *ab, int i)
{
    return ab->ui[i].t1 < 0;
}

static int ar_skip_unit (arbond *ab, int i, int k)
{
    int t;

    for (t=ab->ui[i].t1+k; t<=ab->ui[i].t2; t++) {
	if (!skip_obs(ab, i, t) && !skip_obs(ab, i, t-k)) {
	    return 0;
	}
    }

    return 1;
}

/* see if we have sufficient data to calculate A & B's z statistic for
   AR(k) errors: return the length of the needed arrays, or 0 if we
   can't do it
 */

static int ar_data_check (arbond *ab, int k)
{
    int i, t, T = 0;

    for (i=0; i<ab->N; i++) {
	if (skip_unit(ab, i)) {
	    continue;
	}
	for (t=ab->ui[i].t1+k; t<=ab->ui[i].t2; t++) {
	    if (!skip_obs(ab, i, t) && !skip_obs(ab, i, t-k)) {
		T++;
	    }
	}
    }

    return T;
}

static int arbond_const_pos (arbond *ab)
{
    int i;

    if (ab->xlist == NULL) {
	return 0;
    }

    for (i=1; i<=ab->xlist[0]; i++) {
	if (ab->xlist[i] == 0) {
	    return i;
	}
    }

    return 0;
}

/* FIXME try to detect and omit periodic dummies? */

static int arbond_wald_test (arbond *ab)
{
    gretl_matrix *vcv = NULL;
    gretl_vector *b = NULL;
    double x;
    int cpos = arbond_const_pos(ab);
    int i, j, k, kc = ab->p + ab->nx;
    int ri, rj;
    int err;

    /* position of const in coeff vector? */
    if (cpos != 0) {
	k = kc - 1;
	cpos += ab->p - 1;
    } else {
	k = kc;
	cpos = -1;
    }

    b = gretl_column_vector_alloc(k);
    vcv = gretl_matrix_alloc(k, k);
    if (b == NULL || vcv == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    ri = 0;
    for (i=0; i<kc; i++) {
	if (i != cpos) {
	    b->val[ri++] = ab->beta->val[i];
	}
    }

    ri = 0;
    for (i=0; i<kc; i++) {
	if (i != cpos) {
	    rj = 0;
	    for (j=0; j<kc; j++) {
		if (j != cpos) {
		    x = gretl_matrix_get(ab->vbeta, i, j);
		    gretl_matrix_set(vcv, ri, rj++, x);
		}
	    }
	    ri++;
	}
    } 

    err = gretl_invert_symmetric_matrix(vcv);
    if (err) {
	fprintf(stderr, "arbond_wald_test, error inverting vcv\n");
	goto bailout;
    }
    
    x = gretl_scalar_qform(b, vcv, &err);
    if (err) {
	fprintf(stderr, _("Failed to compute test statistic\n"));
	goto bailout;
    }

    if (!err) {
	ab->wald = x;
	ab->wdf = k;
    }

#if ADEBUG
    fprintf(stderr, "Wald chi^2(%d) = %g\n", k, x);
#endif

 bailout:

    gretl_matrix_free(vcv);
    gretl_vector_free(b);
    
    return err;
}

static int sargan_test (arbond *ab)
{
    gretl_matrix *Zu = NULL;
    gretl_matrix *m1 = NULL;
    int err = 0;

    Zu = gretl_matrix_alloc(ab->m, 1);
    m1 = gretl_matrix_alloc(ab->m, 1);

    if (Zu == NULL || m1 == NULL) {
	gretl_matrix_free(Zu);
	gretl_matrix_free(m1);
	return E_ALLOC;
    }

    gretl_matrix_multiply(ab->ZT, ab->uhat, Zu);
    gretl_matrix_divide_by_scalar(ab->A, ab->effN);
    gretl_matrix_multiply(ab->A, Zu, m1);

    ab->sargan = gretl_matrix_dot_product(Zu, GRETL_MOD_TRANSPOSE,
					  m1, GRETL_MOD_NONE,
					  &err);

    if (ab->step == 1) {
	/* allow for scale factor in H matrix */
	if (ab->opt & OPT_H) {
	    ab->sargan /= ab->s2;
	} else {
	    ab->sargan *= 2.0 / ab->s2; 
	}
    }

#if ADEBUG
    fprintf(stderr, "Sargan df = m - k = %d - %d\n",
	    ab->m, ab->k);
    fprintf(stderr, "Sargan test: Chi-square(%d) = %g\n",
	    ab->m - ab->k, ab->sargan);
#endif

    gretl_matrix_free(Zu);
    gretl_matrix_free(m1);

    return err;
}

/* Compute the z test for AR errors, if possible.  This should
   perhaps be rolled into arbond_variance() for the sake of
   efficiency 
*/

static int ar_test (arbond *ab, const gretl_matrix *C)
{
    gretl_matrix *v = NULL;
    gretl_matrix *vk = NULL;
    gretl_matrix *X = NULL; /* this is the trimmed "X_{*}" */
    gretl_matrix *vkX = NULL; 
    gretl_matrix *tmpk = NULL;
    gretl_matrix *ui = NULL; 
    gretl_matrix *m1 = NULL; 
    gretl_matrix *SZv = NULL;
    
    double x, num, den;
    double den2, den3;
    int s, t, q, Q;
    int sbig, k = 1;
    int i, j, err = 0;

 restart:

    Q = ar_data_check(ab, k);
    if (Q == 0) {
	if (k == 1) {
	    return 0;
	} else {
	    /* don't leak */
	    goto bailout;
	}
    }

#if ADEBUG 
    fprintf(stderr, "AR(%d) test: number of usable obs = %d\n", k, Q);
#endif

    if (k == 1) {
	v = gretl_column_vector_alloc(Q);
	vk = gretl_column_vector_alloc(Q);
	X = gretl_matrix_alloc(Q, ab->k);
	vkX = gretl_matrix_alloc(1, ab->k);
	tmpk = gretl_matrix_alloc(1, ab->k);
	ui = gretl_column_vector_alloc(ab->maxTi);
	m1 = gretl_matrix_alloc(ab->m, 1);
	SZv = gretl_zero_matrix_new(ab->m, 1);

	if (v == NULL || vk == NULL || X == NULL || 
	    vkX == NULL || tmpk == NULL || ui == NULL ||
	    m1 == NULL || SZv == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    } else {
	gretl_matrix_reuse(v, Q, 1);
	gretl_matrix_reuse(vk, Q, 1);
	gretl_matrix_reuse(X, Q, -1);
	gretl_matrix_reuse(m1, ab->m, 1);
	gretl_matrix_zero(SZv);
    }

    den = 0.0;
    q = s = sbig = 0;

    for (i=0; i<ab->N; i++) {
	int Ti = unit_nobs(ab, i);
	double den1i = 0.0;
	int si = 0;

	if (Ti == 0) {
	    continue;
	}

	if (ar_skip_unit(ab, i, k)) {
	    sbig += Ti;
	    s += Ti;
	    continue;
	}

	gretl_matrix_reuse(ui, Ti, -1);
	gretl_matrix_reuse(ab->Zi, ab->m, Ti);

	/* extract full-length Z'_i and u_i */

	for (t=ab->ui[i].t1; t<=ab->ui[i].t2; t++) {
	    if (!skip_obs(ab, i, t)) {
		for (j=0; j<ab->m; j++) {
		    x = gretl_matrix_get(ab->ZT, j, sbig);
		    gretl_matrix_set(ab->Zi, j, si, x);
		}
		x = ab->uhat->val[sbig];
		gretl_vector_set(ui, si, x);
		sbig++;
		si++;
	    }
	}

	for (t=ab->ui[i].t1; t<ab->ui[i].t1+k; t++) {
	    /* skip any obs prior to t1 + k */
	    if (!skip_obs(ab, i, t)) s++;
	}

	/* extract lagged residuals vk along with v_{*} and X_{*} */

	for (t=ab->ui[i].t1+k; t<=ab->ui[i].t2; t++) {
	    if (!skip_obs(ab, i, t)) {
		if (!skip_obs(ab, i, t-k)) {
		    v->val[q] = ab->uhat->val[s];
		    vk->val[q] = ab->uhat->val[s-k];
		    den1i += v->val[q] * vk->val[q];
		    for (j=0; j<ab->k; j++) {
			x = gretl_matrix_get(ab->dX, s, j);
			gretl_matrix_set(X, q, j, x);
		    }
		    q++;
		}
		s++;
	    }
	}

	/* cumulate Z'_i u_i u_i'_* u_{i,-k} */
	gretl_matrix_multiply_by_scalar(ui, den1i);
	gretl_matrix_multiply_mod(ab->Zi, GRETL_MOD_NONE,
				  ui, GRETL_MOD_NONE,
				  SZv, GRETL_MOD_CUMULATE);

	den += den1i * den1i;
    }

#if ADEBUG > 1
    gretl_matrix_print(vk, "lagged residuals");
    gretl_matrix_print(v, "current residuals");
#endif

    /* the numerator */
    num = 0.0;
    for (i=0; i<Q; i++) {
	num += vk->val[i] * v->val[i];
    }

    /* the denominator has three components, the first of which
       is already handled */

    /* required component "vkX" = \hat{v}'_{-k} X_{*} */
    gretl_matrix_multiply_mod(vk, GRETL_MOD_TRANSPOSE,
			      X, GRETL_MOD_NONE,
			      vkX, GRETL_MOD_NONE);
    
    /* vk' X_* (X'ZAZ'X)^{-1} X'ZA(sum stuff) */
    gretl_matrix_multiply(vkX, C, tmpk);
    gretl_matrix_reuse(m1, 1, ab->m);
    gretl_matrix_multiply(tmpk, ab->XZA, m1);
    den2 = gretl_matrix_dot_product(m1, GRETL_MOD_NONE,
				    SZv, GRETL_MOD_NONE,
				    &err);
    den -= 2.0 * den2;

    /* additive term vk' X_* vbeta X_*' vk */
    gretl_matrix_multiply(vkX, ab->vbeta, tmpk);
    den3 = gretl_matrix_dot_product(tmpk, GRETL_MOD_NONE,
				    vkX, GRETL_MOD_TRANSPOSE,
				    &err);
    den += den3;

    if (den < 0) {
	err = E_NAN;
	goto bailout;
    }

    /* now "x" = m1 or m2 */
    x = num / sqrt(den);

#if ADEBUG
    fprintf(stderr, "AR(%d) test: z = %.4g [%.3f]\n", k, x, 
	    normal_pvalue_2(x));
#endif

    if (k == 1) {
	ab->AR1 = x;
	k = 2;
	goto restart;
    } else {
	ab->AR2 = x;
    }

 bailout:

    gretl_matrix_free(v);
    gretl_matrix_free(vk);
    gretl_matrix_free(X);
    gretl_matrix_free(vkX);
    gretl_matrix_free(tmpk);
    gretl_matrix_free(ui);
    gretl_matrix_free(m1);
    gretl_matrix_free(SZv);

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

static int windmeijer_correct (arbond *ab, const gretl_matrix *uhat1,
			       const gretl_matrix *varb1)
{
    gretl_matrix *aV = NULL;  /* standard asymptotic variance */
    gretl_matrix *D = NULL;   /* finite-sample factor */
    gretl_matrix *dWj = NULL; /* one component of the above */
    gretl_matrix *ui = NULL;  /* per-unit residuals */
    gretl_matrix *xij = NULL; /* per-unit X_j values */
    gretl_matrix *TT = NULL;  /* workspace follows */
    gretl_matrix *mT = NULL;  
    gretl_matrix *km = NULL;  
    gretl_matrix *k1 = NULL;  
    int i, j, t;
    int err = 0;

    aV = gretl_matrix_copy(ab->vbeta);

    D = gretl_matrix_alloc(ab->k, ab->k);
    dWj = gretl_matrix_alloc(ab->m, ab->m);
    ui = gretl_column_vector_alloc(ab->maxTi);
    xij = gretl_column_vector_alloc(ab->maxTi);
    TT = gretl_matrix_alloc(ab->maxTi, ab->maxTi);
    mT = gretl_matrix_alloc(ab->m, ab->nobs);
    km = gretl_matrix_alloc(ab->k, ab->m);
    k1 = gretl_matrix_alloc(ab->k, 1);

    if (aV == NULL || D == NULL || dWj == NULL || 
	TT == NULL || mT == NULL || km == NULL || 
	k1 == NULL || ui == NULL || xij == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* form -(1/N) * asyV * XZW^{-1} */
    gretl_matrix_multiply(aV, ab->XZA, ab->kmtmp);
    gretl_matrix_multiply_by_scalar(ab->kmtmp, -1.0 / ab->effN);

    /* form W^{-1}Z'v_2 */
    gretl_matrix_multiply(ab->A, ab->ZT, mT);
    gretl_matrix_multiply(mT, ab->uhat, ab->R1);

    for (j=0; j<ab->k; j++) { /* loop across the X's */
	int s = 0;

	gretl_matrix_zero(dWj);

	/* form dWj = -(1/N) \sum Z_i'(X_j u' + u X_j')Z_i */

	for (i=0; i<ab->N; i++) {
	    int Ti = unit_nobs(ab, i);

	    if (Ti == 0) {
		continue;
	    }

	    gretl_matrix_reuse(ui, Ti, 1);
	    gretl_matrix_reuse(xij, Ti, 1);
	    gretl_matrix_reuse(ab->Zi, Ti, ab->m);
	    gretl_matrix_reuse(TT, Ti, Ti);

	    /* extract ui */
	    for (t=0; t<Ti; t++) {
		ui->val[t] = uhat1->val[s++];
	    }

	    /* extract xij */
	    gretl_matrix_extract_matrix(xij, ab->dX, s - Ti, j,
					GRETL_MOD_NONE);
	    gretl_matrix_multiply_mod(ui, GRETL_MOD_NONE,
				      xij, GRETL_MOD_TRANSPOSE,
				      TT, GRETL_MOD_NONE);
	    gretl_matrix_add_self_transpose(TT);

	    /* extract Zi */
	    gretl_matrix_extract_matrix(ab->Zi, ab->ZT, 0, s - Ti,
					GRETL_MOD_TRANSPOSE);

	    gretl_matrix_qform(ab->Zi, GRETL_MOD_TRANSPOSE,
			       TT, dWj, GRETL_MOD_CUMULATE);
	}

	gretl_matrix_multiply_by_scalar(dWj, -1.0 / ab->effN);

	/* D[.,j] = -aV * XZW^{-1} * dWj * W^{-1}Z'v_2 */
	gretl_matrix_multiply(ab->kmtmp, dWj, km);
	gretl_matrix_multiply(km, ab->R1, k1);

	/* write into appropriate column of D (k x k) */
	for (i=0; i<ab->k; i++) {
	    gretl_matrix_set(D, i, j, k1->val[i]);
	}
    }

    /* add to AsyV: D * AsyV */
    gretl_matrix_multiply_mod(D, GRETL_MOD_NONE,
			      aV, GRETL_MOD_NONE,
			      ab->vbeta, GRETL_MOD_CUMULATE);

    /* add to AsyV: AsyV * D' */
    gretl_matrix_multiply_mod(aV, GRETL_MOD_NONE,
			      D, GRETL_MOD_TRANSPOSE,
			      ab->vbeta, GRETL_MOD_CUMULATE);

    /* add to AsyV: D * var(\hat{\beta}_1) * D' */
    gretl_matrix_qform(D, GRETL_MOD_NONE, varb1,
		       ab->vbeta, GRETL_MOD_CUMULATE);

 bailout:

    gretl_matrix_free(D);
    gretl_matrix_free(dWj);
    gretl_matrix_free(ui);
    gretl_matrix_free(xij);
    gretl_matrix_free(TT);
    gretl_matrix_free(mT);
    gretl_matrix_free(km);
    gretl_matrix_free(k1);
    gretl_matrix_free(aV);

    return err;
}

/* 
   Compute the variance matrix:

   N * C^{-1} * (X'*Z*A_N*\hat{V}_N*A_N*Z'*X) * C^{-1} 

   where C = X'*Z*A_N*Z'*X and
    \hat{V}_N = N^{-1} \sum Z_i'*v_i*v_i'*Z_i,
    (v_i being the step-1 residuals)

   we compute the residuals while we're at it
*/

static int arbond_variance (arbond *ab)
{
    gretl_matrix *kk = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *ui = NULL;
    gretl_matrix *u1 = NULL;

    gretl_matrix *C = ab->den;

    double x, ut, SSR;
    int i, j, t, k, c;
    int err = 0;

    kk = gretl_matrix_alloc(ab->k, ab->k);
    if (kk == NULL) {
	return E_ALLOC;
    }

    if (ab->step == 1) {
	V = gretl_zero_matrix_new(ab->m, ab->m);
	ui = gretl_column_vector_alloc(ab->maxTi);
	if (V == NULL || ui == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    } else {
	u1 = gretl_matrix_copy(ab->uhat);
	if (u1 == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
	/* copy 1-step coefficient variance */
	gretl_matrix_copy_values(kk, ab->vbeta);
    }

    SSR = 0.0;
    c = k = 0;

    for (i=0; i<ab->N; i++) {
	int Ti = unit_nobs(ab, i);

	if (Ti == 0) {
	    continue;
	}

	if (ab->step == 1) {
	    gretl_matrix_reuse(ab->Zi, Ti, ab->m);
	    gretl_matrix_reuse(ui, Ti, 1);
	    gretl_matrix_extract_matrix(ab->Zi, ab->ZT, 0, c,
					GRETL_MOD_TRANSPOSE);
	    c += Ti;
	}

	for (t=0; t<Ti; t++) {
	    ut = ab->dy->val[k]; 
	    for (j=0; j<ab->k; j++) {
		x = gretl_matrix_get(ab->dX, k, j);
		ut -= ab->beta->val[j] * x;
	    }
	    SSR += ut * ut;
	    ab->uhat->val[k] = ut;
	    if (ui != NULL) {
		ui->val[t] = ut;
	    }	
	    k++;
	}

	if (ab->step == 1) {
	    gretl_matrix_multiply_mod(ui, GRETL_MOD_TRANSPOSE,
				      ab->Zi, GRETL_MOD_NONE,
				      ab->L1, GRETL_MOD_NONE);
	    gretl_matrix_multiply_mod(ab->L1, GRETL_MOD_TRANSPOSE,
				      ab->L1, GRETL_MOD_NONE,
				      V, GRETL_MOD_CUMULATE);
	}
    }

    if (ab->step == 1) {
	gretl_matrix_divide_by_scalar(V, ab->effN);

	/* form X'Z A_N Z'X  */
	gretl_matrix_multiply_mod(ab->ZX, GRETL_MOD_TRANSPOSE,
				  ab->A, GRETL_MOD_NONE,
				  ab->kmtmp, GRETL_MOD_NONE);
	gretl_matrix_qform(ab->kmtmp, GRETL_MOD_NONE, V,
			   kk, GRETL_MOD_NONE);

	/* pre- and post-multiply by C^{-1} */
	err = gretl_invert_symmetric_matrix(C);
	gretl_matrix_qform(C, GRETL_MOD_NONE, kk, ab->vbeta,
			   GRETL_MOD_NONE);

	if (ab->opt & OPT_T) {
	    /* preserve V for second stage */
	    ab->V = V;
	    V = NULL;
	}
    } else {
	/* second step */
	gretl_matrix_qform(ab->ZX, GRETL_MOD_TRANSPOSE, ab->V,
			   ab->vbeta, GRETL_MOD_NONE);
	err = gretl_invert_symmetric_matrix(ab->vbeta);
	if (!err) {
	    err = gretl_invert_symmetric_matrix(C); /* for AR test */
	}
    } 

    if (err) {
	goto bailout;
    }

    gretl_matrix_multiply_by_scalar(ab->vbeta, ab->effN);

    ab->SSR = SSR;
    ab->s2 = SSR / (ab->nobs - ab->k);

    if (ab->step == 2 && !(ab->opt & OPT_A)) {
	windmeijer_correct(ab, u1, kk);
    }   

#if ADEBUG
    gretl_matrix_print(ab->vbeta, "Var(beta)");
    for (i=0; i<ab->k; i++) {
	x = gretl_matrix_get(ab->vbeta, i, i);
	fprintf(stderr, "se(beta[%d]) = %g\n", i, sqrt(x));
    }
    fprintf(stderr, "\nSSR = %.11g\n", ab->SSR);
    fprintf(stderr, "sigma^2 = %.7g\n", ab->s2);
    fprintf(stderr, "sigma = %.7g\n", sqrt(ab->s2));
#endif

    /* while we're at it... */
    if (!(ab->opt & OPT_H)) {
	/* FIXME */
	ar_test(ab, C);
    }
    sargan_test(ab);
    arbond_wald_test(ab);

 bailout:

    gretl_matrix_free(V);
    gretl_matrix_free(ui);
    gretl_matrix_free(u1);
    gretl_matrix_free(kk);

    return err;
}

static int next_obs (arbond *ab, int i, int j0, int n)
{
    int j;

    for (j=j0; j<n; j++) {
	if (ab->ui[i].skip[j+ab->ui[i].t1] == 0) {
	    return j;
	}
    }

    return 0;
}

/* construct the H matrix for first-differencing
   as applied to unit i */

static int make_first_diff_matrix (arbond *ab, int i)
{
    static int *rc;

    int n, m;
    double x;
    int k, j, adjacent, skip = 0;

    if (ab == NULL) {
	/* clean-up signal */
	free(rc);
	rc = NULL;
	return 0;
    }

    if (rc == NULL) {
	rc = malloc((ab->T) * sizeof *rc);
	if (rc == NULL) {
	    return E_ALLOC;
	}
    }

    n = ab->ui[i].t2 - ab->ui[i].t1 + 1;
    m = unit_nobs(ab, i);

    if (m < n) {
	skip = 1;
	j = next_obs(ab, i, 0, n);
	for (k=0; k<m; k++) {
	    rc[k] = j;
	    j = next_obs(ab, i, j+1, n);
	}
    }

    gretl_matrix_reuse(ab->H, m, m);

    for (j=0; j<m; j++) {
	for (k=j; k<m; k++) {
	    if (skip) {
		adjacent = (abs(rc[k] - rc[j]) == 1); 
	    } else {
		adjacent = (abs(k-j) == 1);
	    }
	    x = (k==j)? 2 : (adjacent)? -1 : 0;
	    gretl_matrix_set(ab->H, j, k, x);
	    gretl_matrix_set(ab->H, k, j, x);
	}
    }

    return 0;
}

static int arbond_prepare_model (MODEL *pmod, arbond *ab,
				 const int *list, const char *istr,
				 const double **X, const DATAINFO *pdinfo)
{
    const double *y = X[ab->yno];
    char prefix;
    double x;
    int i, j;
    int err = 0;

    pmod->t1 = pdinfo->t1;
    pmod->t2 = pdinfo->t2;
    pmod->dfn = ab->k;
    pmod->dfd = ab->nobs - ab->k;

    pmod->list = gretl_list_copy(list);
    if (pmod->list == NULL) {
	pmod->errcode = E_ALLOC;
	return pmod->errcode;
    }

    gretl_model_set_int(pmod, "yno", ab->yno);
    gretl_model_set_int(pmod, "n_included_units", ab->effN);

    pmod->ci = ARBOND;
    pmod->ncoeff = ab->k;
    pmod->nobs = ab->nobs;
    pmod->full_n = pdinfo->n;
    pmod->ess = ab->SSR;
    if (ab->s2 >= 0) {
	pmod->sigma = sqrt(ab->s2);
    }

    pmod->rsq = pmod->adjrsq = NADBL;
    pmod->fstt = NADBL;
    pmod->lnL = NADBL;
  
    gretl_model_allocate_params(pmod, ab->k);
    if (pmod->errcode) {
	return pmod->errcode;
    }

    prefix = (ab->opt & OPT_H)? 'O' : 'D';

    j = 0;
    for (i=0; i<ab->p; i++) {
	/* FIXME possible varname overflow */
	sprintf(pmod->params[j++], "%c%.10s(-%d)", prefix, pdinfo->varname[ab->yno], i+1);
    }

    for (i=0; i<ab->nx; i++) {
	strcpy(pmod->params[j++], pdinfo->varname[ab->xlist[i+1]]);
    }

    for (i=0; i<ab->ndum; i++) {
	sprintf(pmod->params[j++], "T%d", i + 2);
    }

    err = gretl_model_allocate_storage(pmod);

    if (!err) {
	gretl_model_new_vcv(pmod, NULL);
    }

    if (!err) {
	for (i=0; i<ab->k; i++) {
	    pmod->coeff[i] = ab->beta->val[i];
	    for (j=0; j<=i; j++) {
		x = gretl_matrix_get(ab->vbeta, i, j);
		pmod->vcv[ijton(i, j, ab->k)] = x;
		if (i == j) {
		    pmod->sderr[i] = sqrt(x);
		}
	    }
	}
    }

    if (!err) {
	int s, t, k = 0;

	/* add uhat, yhat */

	for (i=0; i<ab->N; i++) {
	    if (skip_unit(ab, i)) {
		continue;
	    }
	    for (t=0; t<ab->T; t++) {
		if (t >= ab->ui[i].t1 && t <= ab->ui[i].t2) {
		    if (!skip_obs(ab, i, t)) {
			s = data_index(ab, i) + t;
			pmod->uhat[s] = ab->uhat->val[k];
			pmod->yhat[s] = y[s] - pmod->uhat[s];
			k++;
		    }
		}
	    }
	}
    }

    /* additional arbond-specific data */

    if (!err) {
	gretl_model_set_int(pmod, "step", ab->step);
	if (!na(ab->AR1)) {
	    gretl_model_set_double(pmod, "AR1", ab->AR1);
	}
	if (!na(ab->AR2)) {
	    gretl_model_set_double(pmod, "AR2", ab->AR2);
	}
	if (!na(ab->sargan)) {
	    gretl_model_set_int(pmod, "sargan_df", ab->m - ab->k);
	    gretl_model_set_double(pmod, "sargan", ab->sargan);
	}
	if (!na(ab->wald)) {
	    gretl_model_set_int(pmod, "wald_df", ab->wdf);
	    gretl_model_set_double(pmod, "wald", ab->wald);
	}
	if (istr != NULL && *istr != '\0') {
	    gretl_model_set_string_as_data(pmod, "istr", gretl_strdup(istr));
	}
	if (ab->opt & OPT_D) {
	    gretl_model_set_int(pmod, "time-dummies", 1);
	}
	if ((ab->opt & OPT_T) && (ab->opt & OPT_A)) {
	    gretl_model_set_int(pmod, "asy", 1);
	}
    }

    return err;
}

/* exclude any independent variables that are zero at all
   relevant observations */

static int arbond_zero_check (arbond *ab, const double **Z)
{
    const double *x;
    int drop = 0;
    int i, j, k, t;

    for (j=0; j<ab->nx; j++) {
	int all0 = 1;

	x = Z[ab->xlist[j+1]];
	for (i=0; i<ab->N && all0; i++) {
	    if (skip_unit(ab, i)) {
		continue;
	    }
	    for (t=ab->ui[i].t1; t<=ab->ui[i].t2 && all0; t++) {
		k = data_index(ab, i) + t;
		if (x[k] != 0.0 && !na(x[k])) {
		    all0 = 0;
		}
	    }
	}
	if (all0) {
	    ab->xlist[j+1] = -1;
	    drop++;
	}
    }

    if (drop > 0) {
	for (i=1; i<=ab->xlist[0]; i++) {
	    if (ab->xlist[i] < 0) {
		gretl_list_delete_at_pos(ab->xlist, i);
		i--;
	    }
	}
	ab->nx = ab->xlist[0];
	ab->k -= drop;
    }

    return 0;
}

static int count_mask (const char *mask, int n)
{
    int i, c = 0;

    for (i=0; i<n; i++) {
	c += mask[i] == 0;
    }

    return c;
}

static gretl_matrix *
matrix_copy_masked (const gretl_matrix *m, const char *mask)
{
    int n = count_mask(mask, m->rows);
    gretl_matrix *a;
    double x;
    int i, j, k, l;

    a = gretl_matrix_alloc(n, n);
    if (a == NULL) {
	return NULL;
    }

    k = 0;
    for (i=0; i<m->rows; i++) {
	if (!mask[i]) {
	    l = 0;
	    for (j=0; j<m->cols; j++) {
		if (!mask[j]) {
		    x = gretl_matrix_get(m, i, j);
		    gretl_matrix_set(a, k, l++, x);
		}
	    }
	    k++;
	}
    }

    return a;
}

static char *zero_row_mask (const gretl_matrix *m, int *err)
{
    char *mask = NULL;
    int any0 = 0, row0, i, j;

    mask = calloc(m->rows, 1);
    if (mask == NULL) {
	*err = E_ALLOC;
	return NULL;
    }
    
    for (i=0; i<m->rows; i++) {
	row0 = 1;
	for (j=0; j<m->cols; j++) {
	    if (gretl_matrix_get(m, i, j) != 0.0) {
		row0 = 0;
		break;
	    }
	}
	if (row0) {
	    mask[i] = 1;
	    any0 = 1;
	}
    }

    if (!any0) {
	free(mask);
	mask = NULL;
    }

    return mask;
}

static void cut_rows (gretl_matrix *a, const char *mask)
{
    int n = count_mask(mask, a->rows);
    int i, j, k;
    double x;

    for (j=0; j<a->cols; j++) {
	k = 0;
	for (i=0; i<a->rows; i++) {
	    if (!mask[i]) {
		x = gretl_matrix_get(a, i, j);
		a->val[j * n + k] = x;
		k++;
	    }
	}
    }

    a->rows = n;
}

static char *
make_rank_mask (const gretl_matrix *A, int *err)
{
    gretl_matrix *Q = NULL;
    gretl_matrix *R = NULL;
    char *mask = NULL;
    double test;
    int i, n = A->cols;

    Q = gretl_matrix_copy(A);
    if (Q == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    R = gretl_matrix_alloc(n, n);
    if (R == NULL) {
	gretl_matrix_free(Q);
	*err = E_ALLOC;
	return NULL;
    }

    *err = gretl_matrix_QR_decomp(Q, R);

    if (!*err) {
	mask = calloc(n, 1);
	if (mask == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	for (i=0; i<n; i++) {
	    test = gretl_matrix_get(R, i, i);
	    if (fabs(test) < R_DIAG_MIN) {
		mask[i] = 1;
	    }
	}
    }

    if (*err) {
	free(mask);
	mask = NULL;
    }

    gretl_matrix_free(Q);
    gretl_matrix_free(R);
    
    return mask;
}

/* Based on reduction of the A matrix, trim ZT to match and
   adjust the sizes of all workspace matrices that have a 
   dimension involving ab->m.
*/

static void real_shrink_matrices (arbond *ab, const char *mask)
{
    fprintf(stderr, "A matrix: shrinking m from %d to %d\n", 
	    ab->m, ab->A->rows);

    cut_rows(ab->ZT, mask);

    ab->m = ab->A->rows;
    gretl_matrix_reuse(ab->Acpy, ab->m, ab->m);
    gretl_matrix_reuse(ab->tmp1, ab->m, ab->m);
    gretl_matrix_reuse(ab->kmtmp, -1, ab->m);
    gretl_matrix_reuse(ab->L1, -1, ab->m);
    gretl_matrix_reuse(ab->XZA, -1, ab->m);
    gretl_matrix_reuse(ab->R1, ab->m, -1);
    gretl_matrix_reuse(ab->ZX, ab->m, -1);
}

/* Remove zero rows/cols from A, as indicated by mask, and delete the
   corresponding columns from Z.  At this point we leave open the
   question of whether the reduced A matrix is positive definite.
*/

static int 
reduce_Z_and_A (arbond *ab, const char *mask)
{
    gretl_matrix *A;
    int err = 0;

    A = matrix_copy_masked(ab->A, mask);

    if (A == NULL) {
	err = E_ALLOC;
    } else {
	gretl_matrix_free(ab->A);
	ab->A = A;
	real_shrink_matrices(ab, mask);
    }

    return err;
}

/* We already removed any zero rows/columns from A, but it couldn't
   be inverted.  Now we try reducing the dimension of A to its rank.
   We need to read from the backup, Acpy, since ab->A will have been
   mangled in the failed inversion attempt.  We proceed to check that
   the reduced A is invertible, and if not we flag an error.
*/

static int try_alt_inverse (arbond *ab)
{
    char *mask = NULL;
    int err = 0;

    mask = make_rank_mask(ab->Acpy, &err);
    if (err) {
	return err;
    }

    gretl_matrix_free(ab->A);

    ab->A = matrix_copy_masked(ab->Acpy, mask);
    if (ab->A == NULL) {
	free(mask);
	return E_ALLOC;
    }

    err = gretl_invert_symmetric_matrix(ab->A);
    if (!err) {
	real_shrink_matrices(ab, mask);
    } else {
	fprintf(stderr, "try_alt_inverse: error inverting\n");
    }

    free(mask);

    return err;
}

static int arbond_calculate (arbond *ab)
{
    int err = 0;

    /* find Z'X, if this job is not already done */
    if (ab->doZX) {
	gretl_matrix_multiply(ab->ZT, ab->dX, ab->ZX);
	ab->doZX = 0;
    }

#if ADEBUG
    gretl_matrix_print(ab->ZX, "Z'X (arbond_calculate)");
#endif

    /* find X'Z A */
    gretl_matrix_multiply_mod(ab->ZX, GRETL_MOD_TRANSPOSE,
			      ab->A, GRETL_MOD_NONE,
			      ab->XZA, GRETL_MOD_NONE);

    /* calculate "numerator", X'ZAZ'y */
    if (ab->doZy) {
	gretl_matrix_multiply(ab->ZT, ab->dy, ab->R1);
    } else {
	ab->doZy = 1;
    }
#if ADEBUG
    gretl_matrix_print(ab->R1, "Z'y (arbond_calculate)");
#endif
    gretl_matrix_multiply(ab->XZA, ab->R1, ab->beta);

    /* calculate "denominator", X'ZAZ'X */
    gretl_matrix_multiply(ab->XZA, ab->ZX, ab->den);
    gretl_matrix_xtr_symmetric(ab->den);

    gretl_matrix_copy_values(ab->kktmp, ab->den);
    err = gretl_cholesky_decomp_solve(ab->kktmp, ab->beta);

#if ADEBUG
    if (!err) {
	int i;

	fprintf(stderr, "%d-step estimates:\n\n", ab->step);
	for (i=0; i<ab->k; i++) {
	    fprintf(stderr, "beta[%d] = %g\n", i, ab->beta->val[i]);
	}
	fputc('\n', stderr);
    }
#endif

    if (!err) {
	err = arbond_variance(ab);
    }

    return err;
}

static int arbond_step_2 (arbond *ab, PRN *prn)
{
    int err;

#if ADEBUG
    gretl_matrix_print(ab->V, "V, in arbond_step_2");
#endif

    if (gretl_matrix_rows(ab->V) > ab->effN) {
	err = 1; /* we know this case won't work */
    } else {
	gretl_matrix_copy_values(ab->Acpy, ab->V);
	err = gretl_invert_symmetric_matrix(ab->V);
	if (err) {
	    gretl_matrix_copy_values(ab->V, ab->Acpy);
	}
    }

    if (err) {
	err = gretl_SVD_invert_matrix(ab->V);
	if (err) {
	    return err;
	} 
	gretl_matrix_xtr_symmetric(ab->V);
    }

    gretl_matrix_copy_values(ab->A, ab->V);

    ab->step = 2;
    err = arbond_calculate(ab);
    if (err) {
	fprintf(stderr, "step 2: arbond_calculate returned %d\n", err);
    }

    return err;
}

static double odev_at_lag (const double *x, int t, int lag, int pd)
{
    double ret, xbar = 0.0;
    int s, Tt, n = 0;

    t -= lag + 1;

    if (t < 0 || na(x[t])) {
	return NADBL;
    }

    Tt = pd - (t % pd) - (lag + 1);

    for (s=1; s<=Tt; s++) {
	if (!na(x[t+s]) && !na(x[t+s+lag])) {
	    xbar += x[t+s];
	    n++;
	}
    }

    if (n > 0) {
	xbar /= n;
	ret = sqrt(n / (n + 1.0)) * (x[t] - xbar);
    } else {
	ret = NADBL;
    }

    return ret;
}

static int arbond_make_dy_and_X (arbond *ab, const double *y,
				 const double **Z, 
				 const DATAINFO *pdinfo)
{
    int odev = (ab->opt & OPT_H);
    int i, j, s, t, k = 0;
    double x;

    for (i=0; i<ab->N; i++) {
	if (skip_unit(ab, i)) {
	    continue;
	}
	for (t=ab->ui[i].t1; t<=ab->ui[i].t2; t++) {
	    if (skip_obs(ab, i, t)) {
		continue;
	    }
	    s = data_index(ab, i) + t;
	    /* current difference of dependent var */
	    if (odev) {
		x = odev_at_lag(y, s, 0, pdinfo->pd);
		gretl_vector_set(ab->dy, k, x);
	    } else {
		gretl_vector_set(ab->dy, k, y[s] - y[s-1]);
	    }
	    for (j=0; j<ab->p; j++) {
		/* lagged difference of dependent var */
		if (odev) {
		    x = odev_at_lag(y, s, j+1, pdinfo->pd);
		    gretl_matrix_set(ab->dX, k, j, x);
		} else {
		    gretl_matrix_set(ab->dX, k, j, y[s-j-1] - y[s-j-2]);
		}
	    }
	    for (j=0; j<ab->nx; j++) {
		/* independent vars */
		x = Z[ab->xlist[j+1]][s];
		gretl_matrix_set(ab->dX, k, j + ab->p, x);
	    }
	    for (j=0; j<ab->ndum; j++) {
		/* time dummies */
		x = (t - ab->t1min - 1 == j)? 1 : 0;
		gretl_matrix_set(ab->dX, k, j + ab->p + ab->nx, x);
	    }	    
	    k++;
	}
    }

#if ADEBUG
    gretl_matrix_print(ab->dy, "dy");
    gretl_matrix_print(ab->dX, "X");
#endif

    return 0;
}

#define TRY_CUM 0

static int arbond_make_Z_and_A (arbond *ab, const double *y,
				const double **Z)
{
#if TRY_CUM
    gretl_matrix *dXi = NULL, *dyi = NULL;
    int cum_row = 0;
#endif
    int i, j, k, s, t, c = 0;
    int zi, zj, zk;
    double x;
    char *zmask;
#if ADEBUG
    char zstr[8];
#endif
    int err = 0;

#if TRY_CUM
    dXi = gretl_matrix_alloc(ab->Zi->rows, ab->dX->cols);
    dyi = gretl_matrix_alloc(ab->Zi->rows, 1);
#endif

    gretl_matrix_zero(ab->A);
    gretl_matrix_zero(ab->ZX);
    gretl_matrix_zero(ab->R1);

    for (i=0; i<ab->N && !err; i++) {
	int ycols = ab->p; /* intial y block width */
	int offj = 0;      /* initialize column offset */
	int Ti = unit_nobs(ab, i);

	if (Ti == 0) {
	    continue;
	}

	gretl_matrix_reuse(ab->Zi, Ti, ab->m);
	gretl_matrix_zero(ab->Zi);

	for (t=ab->p+1; t<ab->ui[i].t1; t++) {
	    /* unbalanced case: compute initial column offset */
	    offj += ycols;
	    if (ycols < ab->qmax - 1) {
		ycols++;
	    }
	    zj = 0;
	    for (zi=0; zi<ab->nzb; zi++) {
		for (zk=ab->d[zi].minlag; zk<=ab->d[zi].maxlag; zk++) {
		    if (t - zk >= 0) {
			zj++;
		    }
		}
	    }
	    offj += zj;
	}	    

	k = 0;
	for (t=ab->ui[i].t1; t<=ab->ui[i].t2; t++) {
	    int skip = skip_obs(ab, i, t);
	    int offincr = ycols;

	    if (!skip) {
		/* lagged y (GMM instr) columns */
		for (j=0; j<ycols; j++) {
		    s = data_index(ab, i) + t - (ycols + 1) + j;
		    if (!na(y[s])) {
			gretl_matrix_set(ab->Zi, k, j + offj, y[s]);
		    }
		}
	    }
	    
	    /* additional block-diagonal columns, if required --
	       needs checking for the unbalanced case 
	    */
	    zj = 0;
	    for (zi=0; zi<ab->nzb; zi++) {
		for (zk=ab->d[zi].minlag; zk<=ab->d[zi].maxlag; zk++) {
		    if (t - zk >= 0) { /* ?? */
			if (!skip) {
			    s = data_index(ab, i) + t - zk;
			    x = Z[ab->d[zi].v][s];
			    if (!na(x)) {
				gretl_matrix_set(ab->Zi, k, zj + offj + ycols, x);
			    }
			}
			zj++;
		    }
		}
	    }
	    offincr += zj;

	    if (!skip) {
		/* additional full-length instrument columns */
		s = data_index(ab, i) + t;
		for (j=0; j<ab->nz; j++) {
		    x = Z[ab->ilist[j+1]][s];
		    gretl_matrix_set(ab->Zi, k, ab->xc0 + j, x);
		}
		/* plus time dummies, if wanted */
		for (j=0; j<ab->ndum; j++) {
		    x = (t - ab->t1min - 1 == j)? 1 : 0;
		    gretl_matrix_set(ab->Zi, k, ab->xc0 + ab->nz + j, x);
		}
		k++; /* increment target row */	
	    }

	    offj += offincr; /* starting column for next block */
	    if (ycols < ab->qmax - 1) {
		/* increment y block width, if we're not already at max */
		ycols++;
	    }
	}

#if ADEBUG
	sprintf(zstr, "Z_%d", i + 1);
	gretl_matrix_print(ab->Zi, zstr);
#endif

#if TRY_CUM
	/* Try cumulating Z'X and Z'y here: we need to extract the
	   relevant Ti-related slices of dX and dy.
	*/
	int merr = 0;

	gretl_matrix_reuse(dXi, Ti, -1);
	gretl_matrix_reuse(dyi, Ti, -1);
	
	merr += gretl_matrix_extract_matrix(dXi, ab->dX, cum_row, 0,
					    GRETL_MOD_NONE);
	gretl_matrix_print(dXi, "dXi");

	merr += gretl_matrix_extract_matrix(dyi, ab->dy, cum_row, 0,
					    GRETL_MOD_NONE);
	gretl_matrix_print(dyi, "dyi");

	merr += gretl_matrix_multiply_mod(ab->Zi, GRETL_MOD_TRANSPOSE,
					  dXi, GRETL_MOD_NONE,
					  ab->ZX, GRETL_MOD_CUMULATE);

	merr += gretl_matrix_multiply_mod(ab->Zi, GRETL_MOD_TRANSPOSE,
					  dyi, GRETL_MOD_NONE,
					  ab->R1, GRETL_MOD_CUMULATE);

	fprintf(stderr, "merr = %d\n", merr);
	cum_row += Ti;
#if ADEBUG
	gretl_matrix_print(ab->ZX, "Z'X");
	gretl_matrix_print(ab->R1, "Z'y");
#endif
#endif

	/* Cumulate Z_i' H Z_i into A_N */
	if (ab->opt & OPT_H) {
	    /* orthogonal deviations: "H" is identity matrix */
	    gretl_matrix_multiply_mod(ab->Zi, GRETL_MOD_TRANSPOSE,
				      ab->Zi, GRETL_MOD_NONE,
				      ab->A, GRETL_MOD_CUMULATE);
	} else {
	    err = make_first_diff_matrix(ab, i);
	    gretl_matrix_qform(ab->Zi, GRETL_MOD_TRANSPOSE,
			       ab->H, ab->A, GRETL_MOD_CUMULATE);
	}

	/* Write Zi into ZT at offset 0, c */
	gretl_matrix_inscribe_matrix(ab->ZT, ab->Zi, 0, c, GRETL_MOD_TRANSPOSE);
	c += Ti;
    }

    if (!(ab->opt & OPT_H)) {
	/* clean up */
	make_first_diff_matrix(NULL, 0);
    }

#if TRY_CUM
    ab->doZX = 0;
    ab->doZy = 0;
    gretl_matrix_free(dXi);
    gretl_matrix_free(dyi);
#else
    if (!err) {
	/* mask zero rows of ZT, if required */
	zmask = zero_row_mask(ab->ZT, &err);
	if (zmask != NULL) {
	    err = reduce_Z_and_A(ab, zmask);
	    free(zmask);
	}
    } 
#endif

#if ADEBUG
    gretl_matrix_print(ab->A, "\\sum Z_i' H Z_i");
#endif

    if (!err) {
	gretl_matrix_divide_by_scalar(ab->A, ab->N);
    }

#if ADEBUG
    gretl_matrix_print(ab->ZT, "ZT");
    gretl_matrix_print(ab->A, "N^{-1} * \\sum Z_i' H Z_i");
#endif

    return err;
}

static int parse_diag_info (const char *s, struct diag_info *d,
			    const DATAINFO *pdinfo)
{
    char vname[VNAMELEN];
    int v, m1, m2;
    int err = 0;

    if (s == NULL) {
	err = E_ALLOC;
    } else if (sscanf(s, "GMM(%15[^, ],%d,%d)", vname, &m1, &m2) != 3) {
	err = E_PARSE;
    } else {
	if (m2 == 0) {
	    m2 = 99;
	}
	v = varindex(pdinfo, vname);
	if (v == pdinfo->v) {
	    err = E_UNKVAR;
	} else if (m1 < 0 || m2 < m1) {
	    err = E_DATA;
	} else {
	    d->v = v;
	    d->minlag = m1;
	    d->maxlag = m2;
	}
    }

    return err;
}

static int arbond_parse_istr (const char *istr, const DATAINFO *pdinfo,
			      struct diag_info **pd, int *ns)
{
    struct diag_info *d = NULL;
    char *s0 = gretl_strdup(istr);
    char *s, *p, *spec;
    int i, err = 0;

    if (s0 == NULL) {
	return E_ALLOC;
    }

    /* count ')'-terminated fields */
    s = s0;
    while (*s) {
	if (*s == ')') {
	    *ns += 1;
	}
	s++;
    }

    if (*ns == 0) {
	err = E_PARSE;
    } else {
	d = malloc(*ns * sizeof *d);
	if (d == NULL) {
	    err = E_ALLOC;
	}
    }

    /* parse and record individual GMM instrument specs */
    s = s0;
    i = 0;
    while (*s && !err) {
	while (*s == ' ') s++;
	p = s;
	while (*p && *p != ')') p++;
	spec = gretl_strndup(s, p - s + 1);
	err = parse_diag_info(spec, &d[i++], pdinfo);
	free(spec);
	s = p + 1;
    }

    free(s0);

    if (err) {
	free(d);
	*ns = 0;
    } else {
	*pd = d;
    }

    return err;
}

MODEL
arbond_estimate (const int *list, const char *istr, const double **X, 
		 const DATAINFO *pdinfo, gretlopt opt,
		 PRN *prn)
{
    struct diag_info *d = NULL;
    int nzb = 0;
    MODEL mod;
    arbond ab;
    int err = 0;

    gretl_model_init(&mod);
    gretl_model_smpl_init(&mod, pdinfo);

    if (istr != NULL && *istr != 0) {
	mod.errcode = arbond_parse_istr(istr, pdinfo, &d, &nzb);
	if (mod.errcode) {
	    fprintf(stderr, "Error %d in arbond_parse_istr\n", mod.errcode);
	    return mod;
	}
    }

    mod.errcode = arbond_init(&ab, list, pdinfo, opt, d, nzb);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in arbond_init\n", mod.errcode);
	free(d);
	return mod;
    }

    mod.errcode = arbond_sample_check(&ab, X);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in arbond_sample_check\n", mod.errcode);
	goto bailout;
    }

    if (ab.nx > 0) {
	arbond_zero_check(&ab, X);
    }

    mod.errcode = arbond_allocate(&ab);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in arbond_allocate\n", mod.errcode);
	goto bailout;
    }

    /* build the differenced vector \bar{y} plus the
       X matrix
    */
    err = arbond_make_dy_and_X(&ab, X[ab.yno], X, pdinfo);

    /* build instrument matrix blocks, Z_i, and insert into
       big Z' matrix; cumulate first-stage A_N as we go */
    if (!err) {
	err = arbond_make_Z_and_A(&ab, X[ab.yno], X);
    }

    /* invert A_N: we back up ab.A in case inversion fails */
    if (!err) {
	gretl_matrix_copy_values(ab.Acpy, ab.A);
	err = gretl_invert_symmetric_matrix(ab.A);
	if (err) {
	    fprintf(stderr, "inverting ab.A failed on first pass\n");
	    /* failed: try again, reducing A based on its rank */
	    err = try_alt_inverse(&ab);
	}
    }

#if ADEBUG
    gretl_matrix_print(ab.A, "A_N");
#endif

    if (!err) {
	err = arbond_calculate(&ab);
    }

    if (!err && (opt & OPT_T)) {
	err = arbond_step_2(&ab, prn);
    }

 bailout:

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	mod.errcode = arbond_prepare_model(&mod, &ab, list, istr, 
					   X, pdinfo);
    }

    arbond_free(&ab);

    return mod;
}
