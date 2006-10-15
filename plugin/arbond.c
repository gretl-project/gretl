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

#define XPOS 5

typedef struct arbond_ arbond;

struct unit_info {
    int t1;
    int t2;
};

struct arbond_ {
    int yno;              /* ID number of dependent var */
    int p;                /* lag order for dependent variable */
    int q;                /* max lags of y used as instruments */
    int nx;               /* number of exogenous variables */
    int m;                /* number of columns in instrument matrix */
    int N;                /* number of units */
    int T;                /* number of observations per unit */
    int k;                /* number of parameters estimated */
    gretl_matrix *beta;   /* parameter estimates */
    gretl_matrix *vbeta;  /* parameter variance matrix */
    gretl_matrix *H;
    gretl_matrix *A;
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
    gretl_matrix_free(ab->beta);
    gretl_matrix_free(ab->vbeta);
    gretl_matrix_free(ab->H);
    gretl_matrix_free(ab->A);
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

    free(ab->ui);
}

static int 
arbond_init (arbond *ab, const int *list, const DATAINFO *pdinfo)
{
    int T2, NT2;

    ab->p = list[1];
    ab->q = list[2];
    /* FIXME list[3] should be sep: check this */
    ab->yno = list[4];
    if (list[0] > 4) {
	ab->nx = list[0] - XPOS + 1; /* FIXME instr list */
    } else {
	ab->nx = 0;
    }
    ab->N = pdinfo->paninfo->nunits;
    ab->T = pdinfo->n / ab->N;
    ab->m = (ab->T - 2) * (ab->T - 1) / 2 + ab->nx; /* FIXME? */
    ab->k = ab->p + ab->nx;

#if ADEBUG
    fprintf(stderr, "yno = %d, p = %d, q = %d, nx = %d, m = %d, k = %d\n",
	    ab->yno, ab->p, ab->q, ab->nx, ab->m, ab->k);
#endif

    T2 = ab->T - 2;
    NT2 = ab->N * T2;

    ab->beta = NULL;
    ab->vbeta = NULL;
    ab->ZT = NULL;
    ab->H = NULL;
    ab->A = NULL;
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

    ab->beta = gretl_matrix_alloc(ab->k, 1);
    ab->vbeta = gretl_matrix_alloc(ab->k, ab->k);
    ab->ZT = gretl_matrix_alloc(ab->m, NT2);
    ab->H = gretl_matrix_alloc(T2, T2);
    ab->A = gretl_matrix_alloc(ab->m, ab->m);
    ab->Zi = gretl_matrix_alloc(T2, ab->m);
    ab->dy = gretl_column_vector_alloc(NT2);
    ab->dX = gretl_matrix_alloc(NT2, ab->k);

    ab->tmp0 = gretl_matrix_alloc(T2, ab->m);
    ab->tmp1 = gretl_matrix_alloc(ab->m, ab->m);
    ab->tmp2 = gretl_matrix_alloc(ab->k, ab->m);
    ab->L1 = gretl_matrix_alloc(1, ab->m);
    ab->Lk = gretl_matrix_alloc(ab->k, ab->m);
    ab->R1 = gretl_matrix_alloc(ab->m, 1);
    ab->Rk = gretl_matrix_alloc(ab->m, ab->k);

    ab->ui = malloc(ab->N * sizeof *ab->ui);

    if (ab->ZT == NULL || ab->H == NULL || ab->A == NULL ||
	ab->Zi == NULL || ab->dy == NULL || ab->dX == NULL ||
	ab->tmp0 == NULL || ab->tmp1 == NULL || ab->tmp2 == NULL ||
	ab->L1 == NULL || ab->Lk == NULL ||
	ab->R1 == NULL || ab->Rk == NULL ||
	ab->beta == NULL || ab->vbeta == NULL || ab->ui == NULL) {
	arbond_free(ab);
	return E_ALLOC;
    } else {
	return 0;
    }
}

static int 
arbond_sample_check (arbond *ab, const int *list, 
		     const double **Z, const DATAINFO *pdinfo)
{
    int *vlist;
    int i, j, s, t;

    vlist = gretl_list_new(1 + ab->nx);
    vlist[1] = ab->yno;
    for (i=0; i<ab->nx; i++) {
	vlist[i+2] = list[XPOS + i];
    }

    for (i=0; i<ab->N; i++) {
	int t1 = 0, t2 = ab->T - 1;
	int miss = 0;

	s = i * ab->T;
	for (t=t1; t<=t2; t++) {
	    for (j=1; j<=vlist[0]; j++) {
		if (na(Z[vlist[j]][s])) {
		    t1++;
		} else {
		    break;
		}
	    }
	    s++;
	}

	s = (i+1) * ab->T - 1;
	for (t=t2; t>=t1; t--) {
	    for (j=1; j<=vlist[0]; j++) {
		if (na(Z[vlist[j]][s])) {
		    t2--;
		} else {
		    break;
		}
	    }
	    s--;
	}

	s = i * ab->T;
	for (t=t1; t<=t2; t++) {
	    for (j=1; j<=vlist[0]; j++) {
		if (na(Z[vlist[j]][s])) {
		    miss = 1;
		    break;
		}
	    }
	    s++;
	}
	
	/* FIXME figure minimum number of usable observations */

	if (!miss) {
	    ab->ui[i].t1 = t1;
	    ab->ui[i].t2 = t2;
	} else {
	    ab->ui[i].t1 = -1;
	    ab->ui[i].t2 = -1;
	}
	fprintf(stderr, "unit %d: t1=%d, t2=%d, miss=%d\n", 
		i, t1, t2, miss);
    }

    free(vlist);

    return 0;
}

/* 
   Compute the variance matrix:

   N * C^{-1} * (\bar{X}'*Z*A_N*\hat{V}_N*A_N*Z'*\bar{X}) * C^{-1} 

   where C = \bar{X}'*Z*A_N*Z'*\bar{X} and
    \hat{V}_N = N^{-1} \sum Z_i'*v_i*v_i'*Z_i,
    (v_i being the step-1 residuals)
*/

static int arbond_variance (arbond *ab, gretl_matrix *den, PRN *prn)
{
    gretl_matrix *tmp = NULL;
    gretl_matrix *num = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *u = NULL;
    double x, s2;
    int T2 = ab->T - 2;
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

    for (i=0; i<ab->N; i++) {
	gretl_matrix_extract_matrix(ab->Zi, ab->ZT, 0, c,
				    GRETL_MOD_TRANSPOSE);

	for (t=0; t<T2; t++) {
	    k = i * T2 + t;
	    u->val[t] = ab->dy->val[k];
	    for (j=0; j<ab->k; j++) {
		x = gretl_matrix_get(ab->dX, k, j);
		u->val[t] -= ab->beta->val[j] * x;
	    }
	    s2 += u->val[t] * u->val[t];
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
	c += T2;
    }

    gretl_matrix_divide_by_scalar(V, ab->N);

#if ADEBUG > 1
    gretl_matrix_print(V, "V");
#endif

    /* find the central term, A_N * \hat{V}_N * A_N :
       re-use V for this result */
    err += gretl_matrix_multiply(V, ab->A, ab->tmp1);
    err += gretl_matrix_multiply(ab->A, ab->tmp1, V);

    /* find the enclosing term, Z' * \bar{X} */
    err += gretl_matrix_multiply(ab->ZT, ab->dX, ab->Rk);

    /* complete the large "middle bit" */
    gretl_matrix_reuse(ab->tmp2, ab->m, ab->k);
    err += gretl_matrix_multiply(V, ab->Rk, ab->tmp2);
    err += gretl_matrix_multiply_mod(ab->Rk, GRETL_MOD_TRANSPOSE,
				     ab->tmp2, GRETL_MOD_NONE,
				     num);

    /* pre- and post-multiply by C^{-1} */
    gretl_invert_symmetric_matrix(den);
    gretl_matrix_multiply(num, den, tmp);
    gretl_matrix_multiply(den, tmp, ab->vbeta);
    gretl_matrix_multiply_by_scalar(ab->vbeta, ab->N);

    gretl_matrix_print_to_prn(ab->vbeta, "Var(beta)", prn);

    /* just in case, restore original dim of tmp2 */
    gretl_matrix_reuse(ab->tmp2, ab->k, ab->m);

    for (i=0; i<ab->k; i++) {
	x = gretl_matrix_get(ab->vbeta, i, i);
	pprintf(prn, "se(beta[%d]) = %g\n", i, sqrt(x));
    }

    pprintf(prn, "\nRSS = %.11g\n", s2);
    s2 /= ab->N * T2 - 2; /* df correction wanted? (Ox uses it) */
    pprintf(prn, "sigma^2 = %.7g\n", s2);
    pprintf(prn, "sigma = %.7g\n", sqrt(s2));

    gretl_matrix_free(V);
    gretl_matrix_free(u);
    gretl_matrix_free(tmp);
    gretl_matrix_free(num);

    return err;
}

static void make_first_diff_matrix (gretl_matrix *m)
{
    int n = m->rows;
    double x;
    int i, j;

    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    x = (i == j)? 2 : (i == j-1 || i == j+1)? -1 : 0;
	    gretl_matrix_set(m, i, j, x);
	}
    }

#if ADEBUG
    gretl_matrix_print(m, "H");
#endif
}

/* not really ready yet, just for testing.  list should be
   p q ; y xvars */

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
    const double *y;
    double x;
    char zstr[8];
    int i, j, j0, k;
    int s, t, v;
    int c, xc, T2;
    int err = 0;

    gretl_model_init(&mod);

    err = arbond_init(&ab, list, pdinfo);
    if (err) {
	fprintf(stderr, "do_arbond: error %d in arbond_init\n", err);
	mod.errcode = err;
	return mod;
    }

    arbond_sample_check(&ab, list, X, pdinfo);

    T2 = ab.T - 2;
    y = X[ab.yno];

    num = gretl_zero_matrix_new(ab.k, 1);
    den = gretl_matrix_alloc(ab.k, ab.k);
    tmp = gretl_matrix_alloc(ab.k, ab.k);

    if (num == NULL || den == NULL || tmp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* build the differenced vector \bar{y} plus the differenced
       differenced matrix \bar{X}
    */

    s = 0;
    for (i=0; i<ab.N; i++) {
	for (t=2; t<ab.T; t++) {
	    k = i * ab.T + t;
	    gretl_vector_set(ab.dy, s, y[k] - y[k-1]);
	    for (j=0; j<ab.k; j++) {
		if (j < ab.p) {
		    /* lagged difference of dependent var */
		    x = y[k-j-1] - y[k-j-2];
		} else {
		    /* differenced (?) independent var */
		    v = list[XPOS + j - ab.p];
		    x = X[v][k] - X[v][k-1];
		}
		gretl_matrix_set(ab.dX, s, j, x);
	    }
	    s++;
	}
    }

#if ADEBUG
    gretl_matrix_print(ab.dy, "dy");
    gretl_matrix_print(ab.dX, "dX");
#endif

    make_first_diff_matrix(ab.H);

    /* build instrument matrix blocks, Z_i, and insert into
       big Z' matrix; cumulate first-stage A_N as we go */

    gretl_matrix_zero(ab.A);

    c = 0;
    for (i=0; i<ab.N; i++) {
	gretl_matrix_zero(ab.Zi);

	/* lagged dependent var (level) columns */
	j0 = 0;
	for (t=0; t<T2; t++) {
	    k = i * ab.T;
	    for (j=j0; j<j0+t+1; j++) {
		x = y[k++]; 
		gretl_matrix_set(ab.Zi, t, j, x);
	    }
	    j0 = j;
	}

	/* exogenous var (instr) columns */
	xc = (ab.T - 2) * (ab.T - 1) / 2;
	for (j=0; j<ab.nx; j++) {
	    k = i * ab.T + 2;
	    for (t=0; t<T2; t++) {
		/* should we be doing automatic differencing here? */
#if 0
		x = X[list[XPOS + j]][k];
#else
		x = X[list[XPOS + j]][k] - X[list[XPOS + j]][k-1];
#endif
		gretl_matrix_set(ab.Zi, t, xc, x);
		k++;
	    }
	    xc++;
	}

#if ADEBUG
	sprintf(zstr, "Z_%d", i + 1);
	gretl_matrix_print(ab.Zi, zstr);
#endif

	/* Add Z_i' H Z_i to A_N */
	gretl_matrix_multiply(ab.H, ab.Zi, ab.tmp0); 
	gretl_matrix_multiply_mod(ab.Zi, GRETL_MOD_TRANSPOSE,
				  ab.tmp0, GRETL_MOD_NONE,
				  ab.tmp1); 
	gretl_matrix_add_to(ab.A, ab.tmp1);

	/* Write Zi into ZT at offset 0, c */
	gretl_matrix_inscribe_matrix(ab.ZT, ab.Zi, 0, c, GRETL_MOD_TRANSPOSE);
	c += T2;
    }

    gretl_matrix_print(ab.ZT, "ZT");

    gretl_matrix_divide_by_scalar(ab.A, ab.N);
    gretl_matrix_print(ab.A, "N^{-1} * \\sum Z_i' H Z_i");

    gretl_invert_symmetric_matrix(ab.A);
    gretl_matrix_print(ab.A, "A_N");

    /* find "Lk" = \bar{X}' Z A_N */
    gretl_matrix_multiply_mod(ab.dX, GRETL_MOD_TRANSPOSE,
			      ab.ZT, GRETL_MOD_TRANSPOSE,
			      ab.tmp2);
    gretl_matrix_multiply(ab.tmp2, ab.A, ab.Lk);

    /* calculate "numerator", \bar{X}' Z A_N Z' \bar{y} */

    gretl_matrix_multiply(ab.ZT, ab.dy, ab.R1);
    gretl_matrix_multiply(ab.Lk, ab.R1, num);

    /* calculate "denominator", \bar{X}' Z A_N Z' \bar{X} */
    
    gretl_matrix_multiply(ab.ZT, ab.dX, ab.Rk);
    gretl_matrix_multiply(ab.Lk, ab.Rk, den);

    gretl_matrix_copy_values(ab.beta, num);
    gretl_matrix_copy_values(tmp, den);

    err = gretl_LU_solve(tmp, ab.beta);

    pputs(prn, "one-step estimates:\n\n");
    for (i=0; i<ab.k; i++) {
	pprintf(prn, "beta[%d] = %g\n", i, ab.beta->val[i]);
    }
    pputc(prn, '\n');

    err = arbond_variance(&ab, den, prn);

 bailout:

    arbond_free(&ab);
    gretl_matrix_free(num);
    gretl_matrix_free(den);
    gretl_matrix_free(tmp);

    /* We need to package the results into a MODEL struct,
       but for now we just flag an error */
    mod.errcode = 1;

    return mod;
}

