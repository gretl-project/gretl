/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2006 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "libgretl.h"

#define ADEBUG 1

typedef struct arbond_ arbond;

struct unit_info {
    int t1;      /* first usable obs for unit */
    int t2;      /* last usable obs */
    char *skip;  /* mask for obs to be skipped, if any */ 
};

struct arbond_ {
    gretlopt opt;         /* option flags */
    int step;             /* what step are we on? */
    int yno;              /* ID number of dependent var */
    int p;                /* lag order for dependent variable */
    int q;                /* longest lag of y used as instrument */
    int nx;               /* number of independent variables */
    int nz;               /* number of instruments (other than y lags) */
    int m;                /* number of columns in instrument matrix, Z */
    int pc0;              /* column in Z where predet vars start */
    int xc0;              /* column in Z where exog vars start */
    int maxc;             /* number of Zi columns of maximal width */
    int N;                /* number of units */
    int T;                /* total number of observations per unit */
    int maxTi;            /* maximum equations for any given unit */
    int k;                /* number of parameters estimated */
    int nobs;             /* total observations used */
    double SSR;           /* sum of squared residuals */
    double AR1;           /* z statistic for AR(1) errors */
    double AR2;           /* z statistic for AR(2) errors */
    int *xlist;           /* list of independent variables */
    int *ilist;           /* list of instruments */
    int *plist;           /* list of predetermined indep vars */
    gretl_matrix *beta;   /* parameter estimates */
    gretl_matrix *vbeta;  /* parameter variance matrix */
    gretl_matrix *uhat;   /* residuals, differenced version */
    gretl_matrix *H;
    gretl_matrix *A;
    gretl_matrix *V;
    gretl_matrix *ZT;     /* transpose of full instrument matrix */
    gretl_matrix *Zi;     /* per-unit instrument matrix */
    gretl_matrix *dy;     /* first difference of dependent var */
    gretl_matrix *dX;     /* lagged differences of y and indep vars */
    gretl_matrix *tmp0;   /* workspace */
    gretl_matrix *tmp1;
    gretl_matrix *tmp2;
    gretl_matrix *L1;
    gretl_matrix *Lk;
    gretl_matrix *R1;
    gretl_matrix *Rk;
    struct unit_info *ui;
};

static void arbond_free (arbond *ab)
{
    int i;

    gretl_matrix_free(ab->beta);
    gretl_matrix_free(ab->vbeta);
    gretl_matrix_free(ab->uhat);
    gretl_matrix_free(ab->H);
    gretl_matrix_free(ab->A);
    gretl_matrix_free(ab->V);
    gretl_matrix_free(ab->Zi);
    gretl_matrix_free(ab->ZT);
    gretl_matrix_free(ab->dy);
    gretl_matrix_free(ab->dX);
    gretl_matrix_free(ab->tmp0);
    gretl_matrix_free(ab->tmp1);
    gretl_matrix_free(ab->tmp2);
    gretl_matrix_free(ab->L1);
    gretl_matrix_free(ab->Lk);
    gretl_matrix_free(ab->R1);
    gretl_matrix_free(ab->Rk);

    free(ab->xlist);
    free(ab->ilist);
    free(ab->plist);

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
    ab->Zi = gretl_matrix_alloc(T2, ab->m);
    ab->dy = gretl_column_vector_alloc(ab->nobs);
    ab->dX = gretl_matrix_alloc(ab->nobs, ab->k);

    ab->tmp0 = gretl_matrix_alloc(T2, ab->m);
    ab->tmp1 = gretl_matrix_alloc(ab->m, ab->m);
    ab->tmp2 = gretl_matrix_alloc(ab->k, ab->m);
    ab->L1 = gretl_matrix_alloc(1, ab->m);
    ab->Lk = gretl_matrix_alloc(ab->k, ab->m);
    ab->R1 = gretl_matrix_alloc(ab->m, 1);
    ab->Rk = gretl_matrix_alloc(ab->m, ab->k);

    if (ab->ZT == NULL || ab->H == NULL || ab->A == NULL ||
	ab->Zi == NULL || ab->dy == NULL || ab->dX == NULL ||
	ab->tmp0 == NULL || ab->tmp1 == NULL || ab->tmp2 == NULL ||
	ab->L1 == NULL || ab->Lk == NULL ||
	ab->R1 == NULL || ab->Rk == NULL ||
	ab->beta == NULL || ab->vbeta == NULL ||
	ab->uhat == NULL) {
	return E_ALLOC;
    } else {
	return 0;
    }
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

    if (spos > 0) {
	ab->nx = spos - xpos;
	nz = list[0] - spos;
    } else {
	ab->nx = list[0] - (xpos - 1);
    }

    if (ab->nx > 0) {
	/* compose indep vars list */
	ab->xlist = gretl_list_new(ab->nx);
	if (ab->xlist == NULL) {
	    return E_ALLOC;
	} 
	for (i=0; i<ab->nx; i++) {
	    ab->xlist[i+1] = list[xpos + i];
	}
	printlist(ab->xlist, "ab->xlist");
	if (nz == 0) {
	    /* implicit: all x vars are exogenous */
	    ab->ilist = gretl_list_copy(ab->xlist);
	    if (ab->ilist == NULL) {
		return E_ALLOC;
	    }
	    ab->nz = ab->nx;
	}
    }

    if (nz > 0) {
	/* compose instruments list */
	ab->ilist = gretl_list_new(nz);
	if (ab->ilist == NULL) {
	    return E_ALLOC;
	} 
	for (i=0; i<nz; i++) {
	    ab->ilist[i+1] = list[spos + i + 1];
	}
	printlist(ab->ilist, "ab->ilist");
	ab->nz = nz;
    } 

    if (ab->nx > 0 && nz > 0) {
	/* compose predetermined vars list */
	int j, np = 0;

	for (i=1; i<=ab->xlist[0]; i++) {
	    if (gretl_list_position(ab->xlist[i], ab->ilist) == 0) {
		np++;
	    }
	}

	if (np > 0) {
	    ab->plist = gretl_list_new(np);
	    if (ab->plist == NULL) {
		return E_ALLOC;
	    }
	    j = 1;
	    for (i=1; i<=ab->xlist[0]; i++) {
		if (gretl_list_position(ab->xlist[i], ab->ilist) == 0) {
		    ab->plist[j++] = ab->xlist[i];
		}
	    }
	    printlist(ab->plist, "ab->plist");
	}
    }

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
    ab->V = NULL;
    ab->Zi = NULL;
    ab->dy = NULL;
    ab->dX = NULL;
    ab->tmp0 = NULL;
    ab->tmp1 = NULL;
    ab->tmp2 = NULL;
    ab->L1 = NULL;
    ab->Lk = NULL;
    ab->R1 = NULL;
    ab->Rk = NULL;
}

static void arbond_free_lists (arbond *ab)
{
    free(ab->xlist);
    free(ab->ilist);
    free(ab->plist);
}

static int 
arbond_init (arbond *ab, const int *list, const DATAINFO *pdinfo,
	     gretlopt opt)
{
    int i, xpos = 0;
    int err = 0;

    if (list[0] < 3) {
	return E_PARSE;
    }

    ab->p = list[1];

    /* allow the 'q' field to be missing (implicit 0) */
    if (list[2] == LISTSEP) {
	ab->q = 0;
	ab->yno = list[3];
	xpos = 4;
    } else if (list[3] == LISTSEP && list[0] >= 4) {
	ab->q = list[2];
	ab->yno = list[4];
	xpos = 5;
    } else {
	return E_PARSE;
    }

    if (ab->p < 1 || (ab->q != 0 && ab->q < ab->p + 1)) {
	/* is this right? */
	return E_DATA;
    }    

    ab->xlist = NULL;
    ab->ilist = NULL;
    ab->plist = NULL;

    ab->opt = opt;
    ab->step = 1;
    ab->nx = 0;
    ab->nz = 0;

    if (list[0] >= xpos) {
	err = arbond_make_lists(ab, list, xpos);
	if (err) {
	    return err;
	}
    } 

    ab->N = pdinfo->paninfo->nunits;
    ab->T = pdinfo->n / ab->N;
    ab->k = ab->p + ab->nx;
    ab->maxTi = 0;
    ab->m = 0;
    ab->pc0 = 0;
    ab->xc0 = 0;
    ab->maxc = 0;
    ab->nobs = 0;
    ab->SSR = NADBL;
    ab->AR1 = NADBL;
    ab->AR2 = NADBL;

#if ADEBUG
    fprintf(stderr, "yno = %d, p = %d, q = %d, nx = %d, k = %d\n",
	    ab->yno, ab->p, ab->q, ab->nx, ab->k);
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

/* see if there are NAs for the dependent variable, or
   for any of the independent variables, at a given
   observation */

static int anymiss (arbond *ab, const double **Z, int s)
{
    int i;

    for (i=0; i<=ab->p+1; i++) {
	if (na(Z[ab->yno][s-i])) {
	    return 1;
	}
    }

    if (ab->xlist != NULL) {
	for (i=1; i<=ab->xlist[0]; i++) {
	    if (na(Z[ab->xlist[i]][s])) {
		return 1;
	    }
	}
    }

    return 0;
}

static int 
arbond_sample_check (arbond *ab, const int *list, 
		     const double **Z, const DATAINFO *pdinfo)
{
    const double *y = Z[ab->yno];
    char *mask = NULL;
    int t1min = ab->T - 1;
    int t1imin = ab->T - 1;
    int t2max = 0;
    int N = ab->N;
    int i, s, t;
    int err = 0;

    fprintf(stderr, "arbond_sample_check: ab->T = %d\n", ab->T);

    for (i=0; i<ab->N; i++) {
	/* find the first y observation, all units */
	if (t1min > 0) {
	    for (t=0; t<ab->T; t++) {
		s = i * ab->T + t;
		if (!na(y[s])) {
		    if (t < t1min) {
			t1min = t;
		    }
		    break;
		}
	    }
	}
    }

    fprintf(stderr, "initial scan: t1min = %d\n", t1min);

    mask = malloc(ab->T);
    if (mask == NULL) {
	return E_ALLOC;
    }

    if (ab->q == 0) {
	ab->q = ab->T;
    }

    for (i=0; i<ab->N; i++) {
	int t1i = ab->T - 1, t2i = 0; 
	int Ti = 0, usable = 0;

	fprintf(stderr, "Checking unit %d\n", i);

	for (t=0; t<ab->T; t++) {
	    mask[t] = 0;
	}

	/* identify the observations at which we can form the required
	   Delta y terms, have the requisite independent variables,
	   and can construct at least one orthogonality condition
	   using a lagged level of y
	*/

	for (t=ab->p+1; t<ab->T; t++) {
	    s = i * ab->T + t;
	    if (anymiss(ab, Z, s)) {
		mask[t] = 1;
	    } else {
		usable++;
		if (t < t1i) t1i = t;
		if (t > t2i) t2i = t;
	    } 
	}

	if (usable == 0) {
	    N--;
	    ab->ui[i].t1 = -1;
	    ab->ui[i].t2 = -1;
	    fprintf(stderr, "not usable\n");
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
	fprintf(stderr, "t1 = %d, t2 = %d, Ti = %d, usable obs = %d\n", 
		t1i, t2i, Ti, usable);
    }

    /* is t1min actually "reachable"? */
    if (t1min < t1imin - ab->q) {
	t1min = t1imin - ab->q;
    }

    fprintf(stderr, "Number of units with usable observations = %d\n", N);
    fprintf(stderr, "Total usable observations = %d\n", ab->nobs);
    fprintf(stderr, "Max equations for any given unit = %d\n", ab->maxTi);
    fprintf(stderr, "Maximal relevant time-series span: %d to %d = %d\n", 
	    t1min, t2max, t2max - t1min + 1);

    if (N == 0) {
	err = E_MISSDATA;
    } else {
	/* compute the number of lagged-y columns in Zi */
	int tau = t2max - t1min + 1;
	int cbak = ab->p, cols = 0;

	fprintf(stderr, "tau = %d (ab->p = %d)\n", tau, ab->p);
	ab->m = ab->p;
	ab->maxc = 1;
	for (i=1; i<tau-2; i++) {
	    cols = (ab->p + i > ab->q - 1)? ab->q - 1 : ab->p + i;
	    ab->m += cols;
	    if (cols == cbak) {
		ab->maxc += 1;
	    }
	    cbak = cols;
	}
	fprintf(stderr, "'basic' m = %d\n", ab->m);
	ab->q = cols + 1;
#if 0
	/* instruments for predetermined vars */
	if (ab->plist != NULL) {
	    int j, np = ab->plist[0];

	    ab->pc0 = ab->m;
	    for (i=0; i<np; i++) {
		for (j=2; j<tau; j++) {
		    ab->m += j; /* FIXME q */
		}
	    }
	    fprintf(stderr, "m, including predet cols, = %d\n", ab->m);
	}
#endif
	/* record the column where the exogenous vars start */
	ab->xc0 = ab->m;
	ab->m += ab->nx;
	fprintf(stderr, "total m = %d\n", ab->m);
    }

    free(mask);

    return err;
}

static int unit_nobs (arbond *ab, int i)
{
    int t1 = ab->ui[i].t1;
    int t2 = ab->ui[i].t2;
    int n = t2 - t1 + 1;

    if (ab->ui[i].skip != NULL) {
	int t;

	for (t=t1; t<=t2; t++) {
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

/* see if we have sufficient data to calculate A & B's z statistic for
   AR(k) errors: return the length of the needed arrays, or 0 if we
   can't do it
 */

static int ar_data_check (arbond *ab, int k)
{
    int i, t, q, T = 0;

    for (i=0; i<ab->N; i++) {
	q = 0;
	for (t=ab->ui[i].t1+k; t<=ab->ui[i].t2; t++) {
	    if (!skip_obs(ab, i, t) && !skip_obs(ab, i, t-k)) {
		q++;
	    }
	}
	if (q + ab->p + 1 < 5) {
	    fprintf(stderr, "unit %d has only %d ARk rows: m2 is undefined\n", 
		    i, q);
	    return 0;
	} else {
	    T += q;
	}
    }

    return T;
}

/* compute the z test for AR errors, if possible.  FIXME this should
   probably be rolled into arbond_variance() for the sake of
   efficiency */

static int ar_test (arbond *ab, const gretl_matrix *C, int lag,
		    PRN *prn)
{
    gretl_matrix *v = NULL;
    gretl_matrix *v2 = NULL;
    gretl_matrix *X = NULL; /* this is trimmed "X_{*}" */
    gretl_matrix *v2X = NULL; 
    gretl_matrix *tmpk = NULL;
    gretl_matrix *ui = NULL; 
    gretl_matrix *Zi = NULL;
    gretl_matrix *m1 = NULL; 
    gretl_matrix *SZv = NULL;
    gretl_matrix *XZA = ab->Lk; /* already calculated */
    
    double x, num, den;
    double den2, den3;
    int s, t, q, Q = 0;
    int sbig = 0;
    int i, j, err = 0;

    Q = ar_data_check(ab, lag);
    if (Q == 0) {
	return 0;
    }

#if ADEBUG 
    fprintf(stderr, "AR(%d) test: number of usable obs = %d\n", lag, Q);
#endif

    v = gretl_column_vector_alloc(Q);
    v2 = gretl_column_vector_alloc(Q);
    X = gretl_matrix_alloc(Q, ab->k);
    v2X = gretl_matrix_alloc(1, ab->k);
    tmpk = gretl_matrix_alloc(1, ab->k);
    ui = gretl_column_vector_alloc(ab->maxTi);
    Zi = gretl_matrix_alloc(ab->m, ab->maxTi); /* transpose */
    m1 = gretl_matrix_alloc(ab->m, 1);
    SZv = gretl_matrix_alloc(ab->m, 1);

    if (v == NULL || v2 == NULL || X == NULL || 
	v2X == NULL || tmpk == NULL || ui == NULL ||
	Zi == NULL || m1 == NULL || SZv == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_zero(SZv);

    den = 0.0;
    q = s = 0;

    for (i=0; i<ab->N; i++) {
	int Ti = unit_nobs(ab, i);
	double den1i = 0.0;
	int si = 0;

	gretl_matrix_reuse(ui, Ti, -1);
	gretl_matrix_reuse(Zi, -1, Ti);

	for (t=ab->ui[i].t1; t<=ab->ui[i].t2; t++) {
	    if (!skip_obs(ab, i, t)) {
		/* extract full-length Z'_i and u_i */
		for (j=0; j<ab->m; j++) {
		    x = gretl_matrix_get(ab->ZT, j, sbig);
		    gretl_matrix_set(Zi, j, si, x);
		}
		x = ab->uhat->val[sbig];
		gretl_vector_set(ui, si, x);
		sbig++;
		si++;
	    }
	}

	s += lag; /* skip first lag obs per unit */

	for (t=ab->ui[i].t1+lag; t<=ab->ui[i].t2; t++) {
	    if (!skip_obs(ab, i, t)) {
		if (!skip_obs(ab, i, t-lag)) {
		    v->val[q] = ab->uhat->val[s];
		    v2->val[q] = ab->uhat->val[s-lag];
		    den1i += v->val[q] * v2->val[q];
		    for (j=0; j<ab->k; j++) {
			x = gretl_matrix_get(ab->dX, s, j);
			gretl_matrix_set(X, q, j, x);
		    }
		    q++;
		}
		s++;
	    }
	}

	/* cumulate Z'_i u_i u_i'_* u_{i,-2} */
	err = gretl_matrix_multiply(Zi, ui, m1);
	gretl_matrix_multiply_by_scalar(m1, den1i);
	gretl_matrix_add_to(SZv, m1);

	den += den1i * den1i;
    }

#if ADEBUG > 1
    gretl_matrix_print(v2, "lagged residuals");
    gretl_matrix_print(v, "current residuals");
#endif

    /* the numerator */
    num = 0.0;
    for (i=0; i<Q; i++) {
	num += v2->val[i] * v->val[i];
    }

    /* the denominator has three components, the first of which
       is already handled */

    /* required component "v2X" = \hat{v}'_{-2} X_{*} */
    gretl_matrix_multiply_mod(v2, GRETL_MOD_TRANSPOSE,
			      X, GRETL_MOD_NONE,
			      v2X);
    
    /* v2' X_* (X'ZAZ'X)^{-1} X'ZA(sum stuff) */
    gretl_matrix_multiply(v2X, C, tmpk);
    gretl_matrix_reuse(m1, 1, ab->m);
    gretl_matrix_multiply(tmpk, XZA, m1);
    den2 = gretl_matrix_dot_product(m1, GRETL_MOD_NONE,
				    SZv, GRETL_MOD_NONE,
				    &err);
    den -= 2.0 * den2;

    /* additive term v2' X_* vbeta X_*' v2 */
    gretl_matrix_multiply(v2X, ab->vbeta, tmpk);
    den3 = gretl_matrix_dot_product(tmpk, GRETL_MOD_NONE,
				    v2X, GRETL_MOD_TRANSPOSE,
				    &err);
    den += den3;

    /* now "x" = m2 */
    x = num / sqrt(den);

    if (lag == 1) {
	ab->AR1 = x;
    } else if (lag == 2) {
	ab->AR2 = x;
    }

#if ADEBUG
    pprintf(prn, "AR(%d) test: z = %.4g [%.3f]\n", lag, x, 
	    normal_pvalue_2(x));
#endif

    gretl_matrix_free(v);
    gretl_matrix_free(v2);
    gretl_matrix_free(X);
    gretl_matrix_free(v2X);
    gretl_matrix_free(tmpk);
    gretl_matrix_free(ui);
    gretl_matrix_free(Zi);
    gretl_matrix_free(m1);
    gretl_matrix_free(SZv);

    return err;
}

/* 
   Compute the variance matrix:

   N * C^{-1} * (X'*Z*A_N*\hat{V}_N*A_N*Z'*X) * C^{-1} 

   where C = X'*Z*A_N*Z'*X and
    \hat{V}_N = N^{-1} \sum Z_i'*v_i*v_i'*Z_i,
    (v_i being the step-1 residuals)
*/

static int arbond_variance (arbond *ab, gretl_matrix *C, PRN *prn)
{
    static gretl_matrix *Vcpy;
    gretl_matrix *tmp = NULL;
    gretl_matrix *num = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *u = NULL;
    double x, s2;
    int T2 = ab->maxTi;
    int i, j, t, k, c;
    int err = 0;

    tmp = gretl_matrix_alloc(ab->k, ab->k);
    num = gretl_matrix_alloc(ab->k, ab->k);
    V = gretl_zero_matrix_new(ab->m, ab->m);
    u = gretl_column_vector_alloc(T2);

    if (tmp == NULL || num == NULL || V == NULL || u == NULL) {
	return E_ALLOC;
    }

    s2 = 0.0;
    c = 0;
    k = 0;

    if (Vcpy != NULL) {
	gretl_matrix_copy_values(V, Vcpy);
	gretl_matrix_free(Vcpy);
	Vcpy = NULL;
	goto have_v;
    }

    for (i=0; i<ab->N; i++) {
	int Ti;

	if (ab->ui[i].t1 < 0) {
	    continue;
	}

	Ti = unit_nobs(ab, i);

	gretl_matrix_reuse(ab->Zi, Ti, ab->m);
	gretl_matrix_reuse(u, Ti, 1);

	gretl_matrix_extract_matrix(ab->Zi, ab->ZT, 0, c,
				    GRETL_MOD_TRANSPOSE);

	for (t=0; t<Ti; t++) {
	    u->val[t] = ab->dy->val[k];
	    for (j=0; j<ab->k; j++) {
		x = gretl_matrix_get(ab->dX, k, j);
		u->val[t] -= ab->beta->val[j] * x;
	    }
	    s2 += u->val[t] * u->val[t];
	    ab->uhat->val[k] = u->val[t];
	    k++;
	}

#if ADEBUG > 1
	gretl_matrix_print(u, "ui");
#endif

	gretl_matrix_multiply_mod(u, GRETL_MOD_TRANSPOSE,
				  ab->Zi, GRETL_MOD_NONE,
				  ab->L1);

	gretl_matrix_multiply_mod(ab->L1, GRETL_MOD_TRANSPOSE,
				  ab->L1, GRETL_MOD_NONE,
				  ab->tmp1);

#if ADEBUG > 1
	gretl_matrix_print(ab->tmp1, "Z_i'*v_i*v_i'*Z_i");
#endif

	gretl_matrix_add_to(V, ab->tmp1);
	c += Ti;
    }

    gretl_matrix_divide_by_scalar(V, ab->N);

#if ADEBUG > 1
    gretl_matrix_print(V, "V");
#endif

    if (ab->step == 1 && (ab->opt & OPT_T)) {
	Vcpy = gretl_matrix_copy(V);
	ab->V = gretl_matrix_copy(V);
    }

 have_v:

    /* find the central term, A_N * \hat{V}_N * A_N :
       re-use V for this result */
    err += gretl_matrix_multiply(V, ab->A, ab->tmp1);
    err += gretl_matrix_multiply(ab->A, ab->tmp1, V);

    /* find the enclosing term, Z' * X */
    err += gretl_matrix_multiply(ab->ZT, ab->dX, ab->Rk);

    /* complete the large "middle bit" */
    gretl_matrix_reuse(ab->tmp2, ab->m, ab->k);
    err += gretl_matrix_multiply(V, ab->Rk, ab->tmp2);
    err += gretl_matrix_multiply_mod(ab->Rk, GRETL_MOD_TRANSPOSE,
				     ab->tmp2, GRETL_MOD_NONE,
				     num);

    /* pre- and post-multiply by C^{-1} */
    gretl_invert_symmetric_matrix(C);
    gretl_matrix_multiply(num, C, tmp);
    gretl_matrix_multiply(C, tmp, ab->vbeta);
    gretl_matrix_multiply_by_scalar(ab->vbeta, ab->N);

    gretl_matrix_print_to_prn(ab->vbeta, "Var(beta)", prn);

    /* just in case, restore original dim of tmp2 */
    gretl_matrix_reuse(ab->tmp2, ab->k, ab->m);

#if ADEBUG
    for (i=0; i<ab->k; i++) {
	x = gretl_matrix_get(ab->vbeta, i, i);
	pprintf(prn, "se(beta[%d]) = %g\n", i, sqrt(x));
    }
#endif

    ab->SSR = s2;

#if ADEBUG
    pprintf(prn, "\nRSS = %.11g\n", s2);
    s2 /= ab->nobs - ab->k; /* df correction wanted? (Ox uses it) */
    pprintf(prn, "sigma^2 = %.7g\n", s2);
    pprintf(prn, "sigma = %.7g\n", sqrt(s2));
#endif

    /* while we're at it... */
    ar_test(ab, C, 1, prn);
    ar_test(ab, C, 2, prn);

    gretl_matrix_free(V);
    gretl_matrix_free(u);
    gretl_matrix_free(tmp);
    gretl_matrix_free(num);

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

static void make_first_diff_matrix (arbond *ab, int i)
{
    int n, m = unit_nobs(ab, i);
    double x;
    int k, j;

    if (getenv("ARBOND_ADJUST_H") != NULL) {
	n = ab->ui[i].t2 - ab->ui[i].t1 + 1;
    } else {
	n = m;
    }

    gretl_matrix_reuse(ab->H, n, n);

    for (k=0; k<n; k++) {
	for (j=0; j<n; j++) {
	    x = (k == j)? 2 : (k == j-1 || k == j+1)? -1 : 0;
	    gretl_matrix_set(ab->H, k, j, x);
	}
    }

    if (m < n) {
	gretl_matrix *P = gretl_zero_matrix_new(m, n);
	gretl_matrix *L = gretl_matrix_alloc(m, n);
	gretl_matrix *R = gretl_matrix_alloc(m, m);

	j = next_obs(ab, i, 0, n);
	for (k=0; k<m; k++) {
	    gretl_matrix_set(P, k, j, 1.0);
	    j = next_obs(ab, i, j+1, n);
	}

	gretl_matrix_multiply(P, ab->H, L);
	gretl_matrix_multiply_mod(L, GRETL_MOD_NONE,
				  P, GRETL_MOD_TRANSPOSE,
				  R);
	gretl_matrix_reuse(ab->H, m, m);
	gretl_matrix_copy_values(ab->H, R);

	gretl_matrix_free(P);
	gretl_matrix_free(L);
	gretl_matrix_free(R);
    }
}

static int arbond_prepare_model (MODEL *pmod, arbond *ab,
				 const int *list,
				 const DATAINFO *pdinfo)
{
    double x;
    int i, j;
    int err = 0;

    pmod->t1 = pdinfo->t1;
    pmod->t2 = pdinfo->t2;
    pmod->dfd = ab->k;
    pmod->dfd = ab->nobs - ab->k;

    pmod->list = gretl_list_copy(list);
    if (pmod->list == NULL) {
	pmod->errcode = E_ALLOC;
	return pmod->errcode;
    }

    gretl_model_set_int(pmod, "yno", ab->yno);

    pmod->ci = ARBOND;
    pmod->ncoeff = ab->k;
    pmod->nobs = ab->nobs;
    pmod->full_n = pdinfo->n;
    pmod->ess = ab->SSR;
    
    x = ab->SSR / pmod->dfd;
    if (x >= 0.0) {
	pmod->sigma = sqrt(x);
    }

    gretl_model_allocate_params(pmod, ab->k);
    if (pmod->errcode) {
	return pmod->errcode;
    }

    j = 0;
    for (i=0; i<ab->p; i++) {
	/* FIXME possible varname overflow */
	sprintf(pmod->params[j++], "D%.10s(-%d)", pdinfo->varname[ab->yno], i+1);
    }

    for (i=0; i<ab->nx; i++) {
	strcpy(pmod->params[j++], pdinfo->varname[ab->xlist[i+1]]);
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

	for (i=0; i<ab->N; i++) {
	    for (t=0; t<ab->T; t++) {
		if (t >= ab->ui[i].t1 && t <= ab->ui[i].t2) {
		    if (!skip_obs(ab, i, t)) {
			s = i * ab->T + t;
			pmod->uhat[s] = ab->uhat->val[k];
			pmod->yhat[s] = ab->dy->val[k] - ab->uhat->val[k];
			k++;
		    }
		}
	    }
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
	    for (t=ab->ui[i].t1; t<=ab->ui[i].t2 && all0; t++) {
		k = i * ab->T + t;
		if (x[k] != 0.0) {
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

static void 
eliminate_zero_rows (gretl_matrix *a, const char *mask)
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
rank_mask (arbond *ab, const gretl_matrix *A, int *err)
{
    char *mask;
    int pos, cut, n;
    int i, k, r;

    r = gretl_matrix_rank(A, err);
    if (!*err && r == 0) {
	*err = E_SINGULAR;
    }

    if (*err) {
	return NULL;
    }

    cut = A->rows - r;

    mask = calloc(A->rows, 1);
    if (mask == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    pos = ab->xc0; 
    n  = 0;

    /* work down from max lag */
    for (k=ab->q; k>=2 && n<cut; k--) {
	int psub = ab->q - 1;

	for (i=0; i<ab->maxc; i++) {
	    pos -= psub;
	    mask[pos] = 1;
	    n++;
	}
	while (psub >= k) {
	    pos -= psub--;
	    mask[pos] = 1;
	    n++;
	}
	pos = ab->xc0 + (ab->q - k + 1);
    }

    return mask;
}

static int try_alt_inverse (arbond *ab, gretl_matrix *Acpy)
{
    char *mask = NULL;
    gretl_matrix *A = NULL;
    int err = 0;

    mask = zero_row_mask(ab->ZT, &err);

    if (mask == NULL && !err) {
	mask = rank_mask(ab, Acpy, &err);
    } 

    if (err) {
	return err;
    }

    A = matrix_copy_masked(Acpy, mask);
    if (A == NULL) {
	free(mask);
	return E_ALLOC;
    }

    err = gretl_invert_symmetric_matrix(A);

    if (!err && A->rows < ab->m) {
	printf("try_alt_inverse: shrinking matrices\n");
	eliminate_zero_rows(ab->ZT, mask);
	ab->m = A->rows;
	gretl_matrix_free(ab->A);
	ab->A = A;
	A = NULL;
	gretl_matrix_reuse(ab->tmp0, -1, ab->m);
	gretl_matrix_reuse(ab->tmp1, ab->m, ab->m);
	gretl_matrix_reuse(ab->tmp2, -1, ab->m);
	gretl_matrix_reuse(ab->L1, -1, ab->m);
	gretl_matrix_reuse(ab->Lk, -1, ab->m);
	gretl_matrix_reuse(ab->R1, ab->m, -1);
	gretl_matrix_reuse(ab->Rk, ab->m, -1);
    }

    if (A != NULL && A != Acpy) {
	gretl_matrix_free(A);
    }
    free(mask);

    return err;
}

static int arbond_calculate (arbond *ab, gretl_matrix *num,
			     gretl_matrix *den,
			     gretl_matrix *tmp,
			     PRN *prn)
{
    int err = 0;

    /* find "Lk" = X' Z A_N */
    gretl_matrix_multiply_mod(ab->dX, GRETL_MOD_TRANSPOSE,
			      ab->ZT, GRETL_MOD_TRANSPOSE,
			      ab->tmp2);
    gretl_matrix_multiply(ab->tmp2, ab->A, ab->Lk);

    /* calculate "numerator", X' Z A_N Z' y */

    gretl_matrix_multiply(ab->ZT, ab->dy, ab->R1);
    gretl_matrix_multiply(ab->Lk, ab->R1, num);

    /* calculate "denominator", X' Z A_N Z' X */
    
    gretl_matrix_multiply(ab->ZT, ab->dX, ab->Rk);
    gretl_matrix_multiply(ab->Lk, ab->Rk, den);

    gretl_matrix_copy_values(ab->beta, num);
    gretl_matrix_copy_values(tmp, den);

    err = gretl_LU_solve(tmp, ab->beta);

#if ADEBUG
    if (!err) {
	int i;

	pprintf(prn, "%d-step estimates:\n\n", ab->step);
	for (i=0; i<ab->k; i++) {
	    pprintf(prn, "beta[%d] = %g\n", i, ab->beta->val[i]);
	}
	pputc(prn, '\n');
    }
#endif

    if (!err) {
	err = arbond_variance(ab, den, prn);
    }

    return err;
}

static int arbond_step_2 (arbond *ab, gretl_matrix *num,
			  gretl_matrix *den, 
			  gretl_matrix *tmp,
			  PRN *prn)
{
    int err;

    if (ab->V == NULL) {
	return E_ALLOC;
    }

#if ADEBUG
    gretl_matrix_print(ab->V, "V, in arbond_step_2");
#endif

    err = gretl_invert_symmetric_indef_matrix(ab->V);
    if (err) {
	fprintf(stderr, "arbond_step_2: inversion failed\n");
	return err;
    }

    gretl_matrix_copy_values(ab->A, ab->V);

    ab->step = 2;
    err = arbond_calculate(ab, num, den, tmp, prn);
    if (err) {
	fprintf(stderr, "step 2: arbond_calculate returned %d\n", err);
    }

    return err;
}

static int skip_unit (arbond *ab, int i)
{
    return ab->ui[i].t1 < 0;
}

MODEL
arbond_estimate (const int *list, const double **X, 
		 const DATAINFO *pdinfo, gretlopt opt,
		 PRN *prn)
{
    MODEL mod;
    arbond ab;
    gretl_matrix *num = NULL;
    gretl_matrix *den = NULL;
    gretl_matrix *tmp = NULL;
    gretl_matrix *Acpy = NULL;
    const double *y;
    double x;
    char zstr[8];
    int i, j, k;
    int c, s, t;
    int err = 0;

    gretl_model_init(&mod);

    mod.errcode = arbond_init(&ab, list, pdinfo, opt);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in arbond_init\n", mod.errcode);
	return mod;
    }

    mod.errcode = arbond_sample_check(&ab, list, X, pdinfo);
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

    y = X[ab.yno];

    num = gretl_zero_matrix_new(ab.k, 1);
    den = gretl_matrix_alloc(ab.k, ab.k);
    tmp = gretl_matrix_alloc(ab.k, ab.k);
    Acpy = gretl_matrix_alloc(ab.m, ab.m);

    if (num == NULL || den == NULL || tmp == NULL || Acpy == NULL) {
	mod.errcode = E_ALLOC;
	goto bailout;
    }

    /* build the differenced vector \bar{y} plus the
       X matrix
    */

    k = 0;
    for (i=0; i<ab.N; i++) {
	int t1 = ab.ui[i].t1;
	int t2 = ab.ui[i].t2;

	if (skip_unit(&ab, i)) {
	    continue;
	}

	for (t=t1; t<=t2; t++) {
	    if (skip_obs(&ab, i, t)) {
		continue;
	    }
	    s = i * ab.T + t;
	    /* current difference of dependent var */
	    gretl_vector_set(ab.dy, k, y[s] - y[s-1]);
	    for (j=0; j<ab.p; j++) {
		/* lagged difference of dependent var */
		gretl_matrix_set(ab.dX, k, j, y[s-j-1] - y[s-j-2]);
	    }
	    for (j=0; j<ab.nx; j++) {
		/* independent vars */
		x = X[ab.xlist[j+1]][s];
		gretl_matrix_set(ab.dX, k, j + ab.p, x);
	    }
	    k++;
	}
    }

#if ADEBUG
    gretl_matrix_print(ab.dy, "dy");
    gretl_matrix_print(ab.dX, "X");
#endif

    /* build instrument matrix blocks, Z_i, and insert into
       big Z' matrix; cumulate first-stage A_N as we go */

    gretl_matrix_zero(ab.A);

    c = 0;
    for (i=0; i<ab.N; i++) {
	int t1 = ab.ui[i].t1;
	int t2 = ab.ui[i].t2;
	int Ti, ncols, offj;
	int npcols, offpj;

	if (skip_unit(&ab, i)) {
	    continue;
	}

	Ti = unit_nobs(&ab, i);
	
	gretl_matrix_reuse(ab.Zi, Ti, ab.m);
	gretl_matrix_zero(ab.Zi);

	offpj = offj = 0;
	ncols = ab.p;
	npcols = 2;

	for (t=ab.p+1; t<t1; t++) {
	    /* pre-shift fields right as needed */
	    offj += ncols;
	    if (ncols < ab.q - 1) {
		ncols++;
	    }
#if 0
	    if (ab.plist != NULL) {
		offpj += npcols;
		if (npcols < ab.q) { /* FIXME q */
		    npcols++;
		}
	    }
#endif
	}	    

	k = 0;
	for (t=t1; t<=t2; t++) {
	    int skip = skip_obs(&ab, i, t);

	    if (!skip) {
		/* lagged y (GMM instr) columns */
		for (j=0; j<ncols; j++) {
		    s = i * ab.T + t - (ncols + 1) + j;
		    if (t - s <= ab.q && !na(y[s])) {
			gretl_matrix_set(ab.Zi, k, j + offj, y[s]);
		    }
		}
	    }
	    offj += ncols;
	    if (ncols < ab.q - 1) {
		ncols++;
	    }
#if 0
	    if (!skip && ab.plist != NULL) {
		/* predetermined independent vars (FIXME q) */
		int ip;

		for (ip=1; ip<=ab.plist[0]; ip++) {
		    for (j=0; j<npcols; j++) {
			s = i * ab.T + t - (npcols + 1) + j;
			if (t - s <= ab.q && !na(X[ab.plist[ip]][s])) {
			    gretl_matrix_set(ab.Zi, k, ab.pc0 + j + offpj, 
					     X[ab.plist[ip]][s]);
			}
		    }
		}
	    }
	    offpj += npcols;
	    if (npcols < ab.q) {
		npcols++;
	    }
#endif
	    if (!skip) {
		/* additional full-length instrument columns */
		s = i * ab.T + t;
		for (j=0; j<ab.nz; j++) {
		    x = X[ab.ilist[j+1]][s];
		    gretl_matrix_set(ab.Zi, k, ab.xc0 + j, x);
		}
		k++; /* increment target row */	
	    }
	}

#if ADEBUG
	sprintf(zstr, "Z_%d", i + 1);
	gretl_matrix_print(ab.Zi, zstr);
#endif
	
	make_first_diff_matrix(&ab, i);
	gretl_matrix_reuse(ab.tmp0, Ti, ab.m);

	/* Add Z_i' H Z_i to A_N */
	gretl_matrix_multiply(ab.H, ab.Zi, ab.tmp0); 
	gretl_matrix_multiply_mod(ab.Zi, GRETL_MOD_TRANSPOSE,
				  ab.tmp0, GRETL_MOD_NONE,
				  ab.tmp1); 
	gretl_matrix_add_to(ab.A, ab.tmp1);

#if ADEBUG > 1
	gretl_matrix_print(ab.tmp1, "Z'_i H Z_i");
#endif

	/* Write Zi into ZT at offset 0, c */
	gretl_matrix_inscribe_matrix(ab.ZT, ab.Zi, 0, c, GRETL_MOD_TRANSPOSE);
	c += Ti;
    }

    gretl_matrix_divide_by_scalar(ab.A, ab.N);
    gretl_matrix_copy_values(Acpy, ab.A); /* in case inversion fails */

#if ADEBUG
    gretl_matrix_print(ab.ZT, "ZT");
    gretl_matrix_print(ab.A, "N^{-1} * \\sum Z_i' H Z_i");
#endif

    err = gretl_invert_symmetric_matrix(ab.A);
    if (err) {
	err = try_alt_inverse(&ab, Acpy);
    }

#if ADEBUG
    gretl_matrix_print(ab.A, "A_N");
#endif

    if (!err) {
	err = arbond_calculate(&ab, num, den, tmp, prn);
    }

    if (!err && (opt & OPT_T)) {
	err = arbond_step_2(&ab, num, den, tmp, prn);
    }

 bailout:

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	mod.errcode = arbond_prepare_model(&mod, &ab, list, pdinfo);
    }

    arbond_free(&ab);
    gretl_matrix_free(num);
    gretl_matrix_free(den);
    gretl_matrix_free(tmp);
    gretl_matrix_free(Acpy);

    return mod;
}
