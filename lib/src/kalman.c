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
#include "uservar.h"
#include "gretl_func.h"
#include "matrix_extra.h"
#include "libset.h"
#include "kalman.h"

/**
 * SECTION:kalman
 * @short_description: The Kalman filter
 * @title: Kalman
 * @include: gretl/libgretl.h, gretl/kalman.h
 *
 */

#define KDEBUG 0

/* notation (quasi-Hamilton):

   State:   S_{t+1} = F*S_t + v_t,         E(v_t*v_t') = Q

   Obs:     y_t = A'*x_t + H'*S_t + w_t,   E(w_t*w_t') = R

*/

typedef struct crossinfo_ crossinfo;

struct crossinfo_ {
    gretl_matrix *B;  /* r x p */
    gretl_matrix *C;  /* n x p */
    gretl_matrix *BB; /* BB': r x r */
    gretl_matrix *CC; /* CC': n x n */
    gretl_matrix *BC; /* BC': r x n */
};

typedef struct stepinfo_ stepinfo;

struct stepinfo_ {
    gretl_matrix *F;  /* T x (r * r) */
    gretl_matrix *H;  /* T x (r * n) */
};

struct kalman_ {
    int flags;   /* for recording any options */
    int fnlevel; /* level of function execution */

    int r;  /* rows of S = number of elements in state */
    int n;  /* columns of y = number of observables */
    int k;  /* columns of A = number of exogenous vars in obs eqn */
    int p;  /* length of combined disturbance vector */
    int T;  /* rows of y = number of observations */
    int okT; /* T - number of missing observations */
    int t;   /* current time step, when filtering */

    int ifc; /* boolean: obs equation includes an implicit constant? */

    double SSRw;    /* \sum_{t=1}^T e_t^{\prime} V_t^{-1} e_t */
    double sumldet; /* \sum_{t=1}^T ln |V_t| */
    double loglik;  /* log-likelihood */
    double s2;      /* = SSRw / k */

    int nonshift; /* When F is a companion matrix (e.g. in arma), the
		     number of rows of F that do something other than
		     shifting down the elements of the state vector
		  */

    /* continuously updated matrices */
    gretl_matrix *S0; /* r x 1: state vector, before updating */
    gretl_matrix *S1; /* r x 1: state vector, after updating */
    gretl_matrix *P0; /* r x r: MSE matrix, before updating */
    gretl_matrix *P1; /* r x r: MSE matrix, after updating */
    gretl_matrix *e;  /* n x 1: one-step forecast error(s), time t */

    /* input data matrices: note that the order matters for various 
       functions, including user_kalman_recheck_matrices() and the
       matrix_is_varying() macro
    */
    const gretl_matrix *F; /* r x r: state transition matrix */
    const gretl_matrix *A; /* k x n: coeffs on exogenous vars, obs eqn */
    const gretl_matrix *H; /* r x n: coeffs on state variables, obs eqn */
    const gretl_matrix *Q; /* r x r: contemp covariance matrix, state eqn */
    const gretl_matrix *R; /* n x n: contemp covariance matrix, obs eqn */
    const gretl_matrix *mu; /* r x 1: constant term in state transition */
    const gretl_matrix *y; /* T x n: dependent variable vector (or matrix) */
    const gretl_matrix *x; /* T x k: independent variables matrix */
    const gretl_matrix *Sini; /* r x 1: S_{1|0} */
    const gretl_matrix *Pini; /* r x r: P_{1|0} */

    /* optional array of names of input matrices */
    char **mnames;

    /* optional array of function-call strings */
    char **matcalls;

    /* optional matrices for recording extra info */
    gretl_matrix *LL;  /* T x 1: loglikelihood, all time-steps */

    /* optional run-time export matrices */
    gretl_matrix *E;   /* T x n: forecast errors, all time-steps */
    gretl_matrix *V;   /* T x nn: MSE for observables, all time-steps */
    gretl_matrix *S;   /* T x r: state vector, all time-steps */
    gretl_matrix *P;   /* T x nr: MSE for state, all time-steps */
    gretl_matrix *K;   /* T x rn: gain matrix, all time-steps */

    /* structure needed only for cross-correlated case */
    crossinfo *cross;

    /* structure need only when smoothing in the time-varying case */
    stepinfo *step;
    
    /* workspace matrices */
    gretl_matrix_block *Blk; /* holder for the following */
    gretl_matrix *PH;
    gretl_matrix *HPH;
    gretl_matrix *FPH;
    gretl_matrix *Vt;
    gretl_matrix *Ve;
    gretl_matrix *PHV;
    gretl_matrix *Ax;
    gretl_matrix *Kt;
    gretl_matrix *Tmpnn;
    gretl_matrix *Tmprr;
    gretl_matrix *Tmprr_2a;
    gretl_matrix *Tmprr_2b;
    gretl_matrix *Tmpr1;

    void *data; /* handle for attching additional info */
    PRN *prn;   /* verbose printer */
};

/* max number of time-varying matrices: F, A, H, Q, R, mu */
#define NMATCALLS 6

#define arma_ll(K) (K->flags & KALMAN_ARMA_LL)

#define set_kalman_running(K) (K->flags |= KALMAN_FORWARD)
#define set_kalman_stopped(K) (K->flags &= ~KALMAN_FORWARD)
#define kalman_is_running(K)  (K->flags & KALMAN_FORWARD)
#define kalman_simulating(K)  (K->flags & KALMAN_SIM)
#define kalman_checking(K)    (K->flags & KALMAN_CHECK)

/* the matrix in question is not an external named user-matrix,
   and is "owned" by the Kalman struct */
#define kalman_owns_matrix(K,i) (K->mnames != NULL && K->mnames[i][0] == '$')

/* the matrix in question is time-varying (it has an associated 
   function call) */
#define matrix_is_varying(K,i) (K->matcalls != NULL && K->matcalls[i] != NULL)

#define filter_is_varying(K) (K->matcalls != NULL)

#define GENERIC_MATNAME "$Kmat"

static const char *kalman_matrix_name (int sym);

/* symbolic identifiers for input matrices: note that potentially
   time-varying matrices must appear first in the enumeration, and
   the order must match the order of the "input data matrices" in
   the Kalman struct (above).
*/

enum {
    K_F = 0,
    K_A,
    K_H,
    K_Q,
    K_R,
    K_m,
    K_y,
    K_x,
    K_S,
    K_P,
    K_MMAX /* sentinel */
};

void free_stepinfo (kalman *K)
{
    if (K->step != NULL) {
	gretl_matrix_free(K->step->F);
	gretl_matrix_free(K->step->H);
	free(K->step);
	K->step = NULL;
    }
}

void free_crossinfo (crossinfo *c)
{
    gretl_matrix_free(c->B);
    gretl_matrix_free(c->C);
    gretl_matrix_free(c->BB);
    gretl_matrix_free(c->CC);
    gretl_matrix_free(c->BC);

    free(c);
}

#define Q_is_cross_pointer(K) (K->cross != NULL && K->cross->BB != NULL \
                               && K->Q == K->cross->BB)

#define R_is_cross_pointer(K) (K->cross != NULL && K->cross->CC != NULL \
                               && K->R == K->cross->CC)

void kalman_free (kalman *K)
{
    if (K == NULL) {
	return;
    }

    gretl_matrix_free(K->S0);
    gretl_matrix_free(K->P0);
    gretl_matrix_free(K->S1);
    gretl_matrix_free(K->P1);
    gretl_matrix_free(K->e);
    gretl_matrix_free(K->LL);

    if (Q_is_cross_pointer(K)) {
	/* avoid double freeing */
	K->Q = NULL;
    }

    if (R_is_cross_pointer(K)) {
	/* avoid double freeing */
	K->R = NULL;
    }    

    gretl_matrix_block_destroy(K->Blk);

    if (K->mnames != NULL) {
	const gretl_matrix **mptr[] = {
	    &K->F, &K->A, &K->H, &K->Q, &K->R,
	    &K->mu, &K->y, &K->x, &K->Sini, &K->Pini
	};
	int i;

	for (i=0; i<K_MMAX; i++) {
	    if (kalman_owns_matrix(K, i)) {
		gretl_matrix_free((gretl_matrix *) *mptr[i]);
	    }
	}

	strings_array_free(K->mnames, K_MMAX);
    }    

    if (K->matcalls != NULL) {
	strings_array_free(K->matcalls, NMATCALLS);
    }

    if (K->cross != NULL) {
	free_crossinfo(K->cross);
    }

    if (K->step != NULL) {
	free_stepinfo(K);
    }    

    free(K);
}

static kalman *kalman_new_empty (int flags)
{
    kalman *K = malloc(sizeof *K);

    if (K != NULL) {
	K->Sini = K->Pini = NULL;
	K->S0 = K->S1 = NULL;
	K->P0 = K->P1 = NULL;
	K->LL = NULL;
	K->e = NULL;
	K->Blk = NULL;
	K->F = K->A = K->H = NULL;
	K->Q = K->R = NULL;
	K->E = K->V = K->S = K->P = K->K = NULL;
	K->y = K->x = NULL;
	K->mu = NULL;
	K->mnames = NULL;
	K->matcalls = NULL;
	K->cross = NULL;
	K->step = NULL;
	K->flags = flags;
	K->fnlevel = 0;
	K->t = 0;
	K->prn = NULL;
	K->data = NULL;
    }

    return K;
}

#define kappa 1.0e7

static void diffuse_Pini (kalman *K)
{
    int i;

    gretl_matrix_zero(K->P0);

    for (i=0; i<K->r; i++) {
	gretl_matrix_set(K->P0, i, i, kappa);
    }
}

static int F_out_of_bounds (kalman *K)
{
    gretl_matrix *Fcpy;
    gretl_matrix *evals;
    double r, c, x;
    int i, err = 0;

    Fcpy = gretl_matrix_copy(K->F);
    if (Fcpy == NULL) {
	return E_ALLOC;
    }

    evals = gretl_general_matrix_eigenvals(Fcpy, 0, &err);
    gretl_matrix_free(Fcpy);

    for (i=0; i<evals->rows && !err; i++) {
	r = gretl_matrix_get(evals, i, 0);
	c = gretl_matrix_get(evals, i, 1);
	x = sqrt(r*r + c*c);
	if (x >= 1.0) {
	    fprintf(stderr, "F: modulus of eigenvalue %d = %g\n", i+1, x);
	    err = E_SINGULAR;
	}
    }

    gretl_matrix_free(evals);

    return err;
}

/* If the user has not given an initial value for P_{1|0}, compute
   this automatically as per Hamilton, ch 13, p. 378.  This works only
   if the eigenvalues of K->F lie inside the unit circle.  Failing
   that, or if the --diffuse option is given for the user Kalman
   filter, we apply a diffuse initialization.
*/

static int construct_Pini (kalman *K)
{
    gretl_matrix *Svar;
    gretl_matrix *vQ; 
    int r2, err = 0;

    if (K->flags & KALMAN_DIFFUSE) {
	diffuse_Pini(K);
	return 0;
    }

    r2 = K->r * K->r;

    Svar = gretl_matrix_alloc(r2, r2);
    vQ = gretl_column_vector_alloc(r2);

    if (Svar == NULL || vQ == NULL) {
	gretl_matrix_free(Svar);
	gretl_matrix_free(vQ);
	return E_ALLOC;
    }

    gretl_matrix_kronecker_product(K->F, K->F, Svar);
    gretl_matrix_I_minus(Svar);
    gretl_matrix_vectorize(vQ, K->Q);

    err = gretl_LU_solve(Svar, vQ);
    if (err) {
	/* failed: are some of the eigenvalues out of bounds? */
	err = F_out_of_bounds(K);
	if (err == E_SINGULAR) {
	    err = 0;
	    diffuse_Pini(K);
	    K->flags |= KALMAN_DIFFUSE;
	}
    } else {
	gretl_matrix_unvectorize(K->P0, vQ);
    }

    gretl_matrix_free(Svar);
    gretl_matrix_free(vQ);

    return err;
}

static int check_matrix_dims (kalman *K, const gretl_matrix *m, int i)
{
    int r = 0, c = 0, symm = (i == K_Q || i == K_R);
    int err = 0;

    if (i == K_F || i == K_Q || i == K_P) {
	r = c = K->r;
    } else if (i == K_A) {
	r = K->k;
	c = K->n;
    } else if (i == K_H) {
	r = K->r;
	c = K->n;
    } else if (i == K_R)  {
	r = c = K->n;
    } else if (i == K_S || i == K_m) {
	r = K->r;
	c = 1;
    } 

    if (m->rows != r || m->cols != c) {
	gretl_errmsg_sprintf("kalman: %s is %d x %d, should be %d x %d\n", 
			     kalman_matrix_name(i), m->rows, m->cols, r, c);
	err = E_NONCONF;
    } else if (symm && !gretl_matrix_is_symmetric(m)) {
	gretl_errmsg_sprintf("kalman: %s is not symmetric\n", 
		kalman_matrix_name(i));
	err = E_NONCONF;
    }	

    return err;
}

enum {
    K_E,
    K_V,
    K_BIG_S,
    K_BIG_P,
    K_K,
    K_LL,
};

static int maybe_resize_export_matrix (kalman *K, gretl_matrix *m, int i)
{
    int r = K->T, c = 0;
    int err = 0;

    if (i == K_E) {
	c = K->n;
    } else if (i == K_V) {
	c = (K->n * K->n + K->n) / 2;
    } else if (i == K_BIG_S) {
	c = K->r;
    } else if (i == K_BIG_P) {
	c = (K->r * K->r + K->r) / 2;
    } else if (i == K_LL) {
	c = 1;
    } else if (i == K_K) {
	c = K->r * K->n;
    } else {
	err = E_DATA;
    }

    if (!err && (m->rows != r || m->cols != c)) {
	err = gretl_matrix_realloc(m, r, c);
    }

    return err;
}

static int missing_matrix_error (const char *name)
{
    if (name == NULL) {
	gretl_errmsg_set(_("kalman: a required matrix is missing"));
    } else {
	gretl_errmsg_sprintf(_("kalman: required matrix %s is missing"), 
			     name);
    }
    return E_DATA;
}

static int missing_kalman_error (void)
{
    gretl_errmsg_set(_("No Kalman filter is defined"));
    return E_DATA;
}

static int kalman_check_dimensions (kalman *K)
{
    int err = 0;

    if (K->r < 1 || K->n < 1 || K->T < 2) {
	/* the state and observation vectors must have at least one
	   element, and there must be at least two observations */
	err = E_DATA;
    }

    /* F is mandatory, should be r x r */
    if (!err) {
	err = check_matrix_dims(K, K->F, K_F);
    }

    /* H is mandatory, should be r x n */
    if (!err) {
	err = check_matrix_dims(K, K->H, K_H);
    }

    /* Q is mandatory, should be r x r and symmetric */
    if (!err) {
	err = check_matrix_dims(K, K->Q, K_Q);
    }

    /* R should be n x n and symmetric, if present */
    if (!err && K->R != NULL) {
	err = check_matrix_dims(K, K->R, K_R);
    }

    /* initial S should be r x 1, if present */
    if (!err && K->Sini != NULL) {
	err = check_matrix_dims(K, K->Sini, K_S);
    }

    /* initial P should be r x r, if present */
    if (!err && K->Pini != NULL) {
	err = check_matrix_dims(K, K->Pini, K_P);
    }

    /* A should be k x n, if present (A->rows defines k; n is defined by y) */
    if (!err) {
	if (K->A != NULL) {
	    err = check_matrix_dims(K, K->A, K_A);
	} else if (K->x != NULL) {
	    /* A is NULL => can't have a non-NULL x */
	    err = E_NONCONF;
	}
    }

    /* mu should be r x 1, if present */
    if (!err && K->mu != NULL) {
	err = check_matrix_dims(K, K->mu, K_m);
    }    

    if (err) {
	return err;
    }

    K->ifc = 0;

    /* x should have T rows to match y; and it should have either k or k - 1
       columns (the latter case indicating an implicit const) */
    if (K->x != NULL) {
	if (K->x->rows != K->T) {
	    fprintf(stderr, "kalman: %s has %d rows, should have %d\n", 
		    kalman_matrix_name(K_x), K->x->rows, K->T);
	    return E_NONCONF;
	} else if (K->x->cols != K->k && K->x->cols != K->k - 1) {
	    fprintf(stderr, "kalman: %s has %d columns, should have %d or %d\n", 
		    kalman_matrix_name(K_x), K->x->cols, K->k, K->k - 1);
	    return E_NONCONF;
	} else if (K->x->cols == K->k - 1) {
	    /* register the implicit const */
	    K->ifc = 1;
	}
    } else if (K->k == 1) {
	/* A has one row but there's no x => implicit const */
	K->ifc = 1;
    } else if (K->k > 1) {
	/* A has more than one row but there's no x => error */
	return missing_matrix_error("xmat");
    }

    /* Below we have the optional "export" matrices for shipping out 
       results.  If these are present but not sized correctly we'll
       try to fix them up.
    */

    /* E should be T x n */
    if (K->E != NULL) {
	err = maybe_resize_export_matrix(K, K->E, K_E);
    }

    /* V should be T x nn */
    if (!err && K->V != NULL) {
	err = maybe_resize_export_matrix(K, K->V, K_V);
    }

    /* big S should be T x r */
    if (!err && K->S != NULL) {
	err = maybe_resize_export_matrix(K, K->S, K_BIG_S);
    }

    /* big P should be T x nr */
    if (!err && K->P != NULL) {
	err = maybe_resize_export_matrix(K, K->P, K_BIG_P);
    }

    /* LL should be T x 1 */
    if (!err && K->LL != NULL) {
	err = maybe_resize_export_matrix(K, K->LL, K_LL);
    } 

    /* K (gain) should be T x (r * n) */
    if (!err && K->K != NULL) {
	err = maybe_resize_export_matrix(K, K->K, K_K);
    }    

    return 0;
}

/* Write the vech of @src into row @t of @targ */

static void load_to_vech (gretl_matrix *targ, const gretl_matrix *src,
			  int n, int t)
{
    int i, j, m = 0;
    double x;

    for (i=0; i<n; i++) {
	for (j=i; j<n; j++) {
	    x = gretl_matrix_get(src, i, j);
	    gretl_matrix_set(targ, t, m++, x);
	}
    }
}

/* Write the vec of @src into row @t of @targ */

static void load_to_vec (gretl_matrix *targ, const gretl_matrix *src,
			 int t)
{
    int i;

    for (i=0; i<targ->cols; i++) {
	gretl_matrix_set(targ, t, i, src->val[i]);
    }
}

/* copy from vector @src into row @t of @targ */

static void load_to_row (gretl_matrix *targ, const gretl_vector *src, 
			 int t)
{
    double x;
    int i;

    for (i=0; i<targ->cols; i++) {
	x = gretl_vector_get(src, i);
	gretl_matrix_set(targ, t, i, x);
    }
}

static void set_row (gretl_matrix *targ, int t, double x)
{
    int i;

    for (i=0; i<targ->cols; i++) {
	gretl_matrix_set(targ, t, i, x);
    }
}

/* checks if row has a 1 in position n and 0s elsewhere */

static int ok_companion_row (const gretl_matrix *F, 
			     int cols, int row, 
			     int n) 
{
    double x;
    int j;

    for (j=0; j<cols; j++) {
	x = gretl_matrix_get(F, row, j);
	if ((j == n && x != 1.0) || (j != n && x != 0.0)) {
	    return 0;
	}
    }

    return 1;
}

static int count_nonshifts (const gretl_matrix *F)
{
    int rmax = gretl_matrix_rows(F);
    double x;
    int i, j, n = 0;
    int ret = rmax;

    /* examine bottom row */
    for (j=0; j<rmax; j++) {
	x = gretl_matrix_get(F, rmax - 1, j);
	if (x != 0.0) {
	    if (n == 0 && j > 0 && x == 1.0) {
		n = j;
	    } else {
		return rmax;
	    }
	}
    }

    /* at this point n holds the candidate; now check n rows from 
       the (r - n - 1)-th */
    if (n > 0) {
	ret = n;
	j = 0;
	for (i=rmax-n; i<rmax; i++) {
	    if (!ok_companion_row(F, rmax, i, j++)) {
		return rmax;
	    }
	}
    }

    return ret;
}

static int kalman_init (kalman *K)
{
    int err = 0;

    K->SSRw = NADBL;
    K->sumldet = NADBL;
    K->loglik = NADBL;
    K->s2 = NADBL;

    clear_gretl_matrix_err();

    if (K->Sini != NULL) {
	K->S0 = gretl_matrix_copy(K->Sini);
	K->S1 = gretl_matrix_copy(K->Sini);
    } else {
	K->S0 = gretl_zero_matrix_new(K->r, 1);
	K->S1 = gretl_zero_matrix_new(K->r, 1);
    }	

    if (K->Pini != NULL) {
	K->P0 = gretl_matrix_copy(K->Pini);
	K->P1 = gretl_matrix_copy(K->Pini);
    } else {
	K->P0 = gretl_zero_matrix_new(K->r, K->r);
	K->P1 = gretl_zero_matrix_new(K->r, K->r);
    }	

    /* forecast error vector, per observation */
    K->e = gretl_matrix_alloc(K->n, 1);

    err = get_gretl_matrix_err();
    if (err) {
	return err;
    }

    K->nonshift = -1;

    K->Blk = gretl_matrix_block_new(&K->PH,  K->r, K->n, /* P*H */
				    &K->FPH, K->r, K->n, /* F*P*H */
				    &K->HPH, K->n, K->n, /* H'*P*H */
				    &K->Vt,  K->n, K->n, /* (H'*P*H + R)^{-1} */
				    &K->Ve,  K->n, 1,    /* (H'*P*H + R)^{-1} * e */
				    &K->PHV, K->r, K->n, /* P*H*V */
				    &K->Ax,  K->n, 1,    /* A'*x at obs t */
				    &K->Kt,  K->r, K->n, /* gain at t */
				    &K->Tmpnn, K->n, K->n,
				    &K->Tmprr, K->r, K->r,
				    &K->Tmprr_2a, K->r, K->r,
				    &K->Tmprr_2b, K->r, K->r,
				    &K->Tmpr1, K->r, 1,
				    NULL);

    if (K->Blk == NULL) {
	err = E_ALLOC;
    }

    if (!err && K->Pini == NULL && !(K->flags & KALMAN_USER)) {
	/* in the "user" case we do this later */
	err = construct_Pini(K);
    }

    return err;
}

static void kalman_set_dimensions (kalman *K, gretlopt opt)
{
    K->r = gretl_matrix_rows(K->F); /* F->rows defines r */
    K->k = gretl_matrix_rows(K->A); /* A->rows defines k */
    K->T = gretl_matrix_rows(K->y); /* y->rows defines T */
    K->n = gretl_matrix_cols(K->y); /* y->cols defines n */

    K->okT = K->T;

    /* K->p is non-zero only under cross-correlation; in that case the
       matrix given as 'Q' in Kalman set-up in fact represents B (as
       in v_t = B \varepsilon_t) and it must be r x p, where p is the
       number of elements in the "combined" disturbance vector
       \varepsilon_t.
    */
    K->p = (opt & OPT_C)? gretl_matrix_cols(K->Q): 0;
}

/**
 * kalman_new:
 * @S: r x 1 initial state vector.
 * @P: r x r initial precision matrix.
 * @F: r x r state transition matrix.
 * @A: n x k matrix of coefficients on exogenous variables in the
 * observation equation.
 * @H: n x r matrix of coefficients on the state variables in the
 * observation equation.
 * @Q: r x r contemporaneous covariance matrix for the errors in the
 * state equation.
 * @R: n x n contemporaneous covariance matrix for the errors in the 
 * observation equation (or NULL if this is not applicable).
 * @y: T x n matrix of dependent variable(s).
 * @x: T x k matrix of exogenous variable(s).  May be NULL if there
 * are no exogenous variables, or if there's only a constant.
 * @m: r x 1 vector of constants in the state transition, or NULL.
 * @E: T x n matrix in which to record forecast errors (or NULL if
 * this is not required).
 * @err: location to receive error code.
 *
 * Allocates and initializes a Kalman struct, which can subsequently
 * be used for forecasting with kalman_forecast().  The nomenclature
 * for the various required matrices is that in Hamilton's Time
 * Series Analysis (1994, chapter 13), except that "S" is used in
 * place of Hamilton's \xi for the state vector.
 *
 * Returns: pointer to allocated struct, or NULL on failure, in
 * which case @err will receive a non-zero code.
 */

kalman *kalman_new (const gretl_matrix *S, const gretl_matrix *P,
		    const gretl_matrix *F, const gretl_matrix *A,
		    const gretl_matrix *H, const gretl_matrix *Q,
		    const gretl_matrix *R, const gretl_matrix *y,
		    const gretl_matrix *x, const gretl_matrix *m,
		    gretl_matrix *E, int *err)
{
    kalman *K;

    *err = 0;

    if (F == NULL || H == NULL || Q == NULL || y == NULL) {
	*err = missing_matrix_error(NULL);
	return NULL;
    }

    K = kalman_new_empty(0);
    if (K == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* use const pointers for const matrices, don't copy */
    K->F = F;
    K->A = A;
    K->H = H;
    K->Q = Q;
    K->R = R;
    K->y = y;
    K->x = x;
    K->Sini = S;
    K->Pini = P;
    K->mu = m;

    /* non-const, but again use external pointer */
    K->E = E;

    kalman_set_dimensions(K, OPT_NONE);

    *err = kalman_check_dimensions(K);
    if (*err) {
	fprintf(stderr, "failed on kalman_check_dimensions\n");
	free(K);
	return NULL;
    }

    *err = kalman_init(K);

    if (*err) {
	kalman_free(K);
	K = NULL;
    } else {
	gretl_matrix_zero(K->e);
    }

    return K;
}

enum {
    UPDATE_INIT, /* initialization of matrices */
    UPDATE_STEP  /* refreshing matrices per time-step */
};

/* After reading 'Q' = B and 'R' = C from the user, either at
   (re-)initializaton or at a given time-step in the case where either
   of these matrices is time-varying, record the user input in K->B
   and K->C and form the 'real' Q and R.  But note that in the
   time-step case it may be that only one of Q, R needs to be treated
   in this way (if only one is time-varying, only one will have 
   been redefined via a function call).
*/

static int kalman_update_crossinfo (kalman *K, int mode)
{
    crossinfo *c = K->cross;
    int err = 0;

    /* Note that B and C may be needed for simulation */

    if (mode == UPDATE_INIT || matrix_is_varying(K, K_Q)) {
	err = gretl_matrix_copy_values(c->B, K->Q);
	if (!err) {
	    err = gretl_matrix_multiply_mod(c->B, GRETL_MOD_NONE,
					    c->B, GRETL_MOD_TRANSPOSE,
					    c->BB, GRETL_MOD_NONE);
	}
    }

    if (!err && (mode == UPDATE_INIT || matrix_is_varying(K, K_R))) {
	err = gretl_matrix_copy_values(c->C, K->R);
	if (!err) {
	    err = gretl_matrix_multiply_mod(c->C, GRETL_MOD_NONE,
					    c->C, GRETL_MOD_TRANSPOSE,
					    c->CC, GRETL_MOD_NONE);
	}
    }

    if (!err && (mode == UPDATE_INIT || matrix_is_varying(K, K_Q) ||
		 matrix_is_varying(K, K_R))) {
	err = gretl_matrix_multiply_mod(c->B, GRETL_MOD_NONE,
					c->C, GRETL_MOD_TRANSPOSE,
					c->BC, GRETL_MOD_NONE);
    }

    if (mode == UPDATE_STEP) {
	if (matrix_is_varying(K, K_Q)) {
	    K->Q = c->BB;
	}
	if (matrix_is_varying(K, K_R)) {
	    K->R = c->CC;
	}
    }

    return err;
}

static int kalman_add_crossinfo (kalman *K)
{
    crossinfo *c = malloc(sizeof *c);
    
    if (c == NULL) {
	return E_ALLOC;
    } 

    c->B = gretl_matrix_alloc(K->r, K->p);
    c->C = gretl_matrix_alloc(K->n, K->p);

    c->BB = gretl_matrix_alloc(K->r, K->r);
    c->CC = gretl_matrix_alloc(K->n, K->n);
    c->BC = gretl_matrix_alloc(K->r, K->n);

    if (c->B == NULL || c->C == NULL || c->BB == NULL || 
	c->CC == NULL || c->BC == NULL) {
	free_crossinfo(c);
	return E_ALLOC;
    }

    K->cross = c;

    return 0;
}

/* kalman_revise_variance: the user has actually given B in place of Q
   and C in place of R; so we have to form Q = BB', R = CC' and BC'.

   This function is called in the course of initial set-up of a
   filter, and also when the Kalman matrices are being re-checked at
   the start of filtering, smoothing or simulation.
*/

static int kalman_revise_variance (kalman *K)
{
    int err = 0;

    if (K->Q == NULL || K->R == NULL) {
	return missing_matrix_error("'statevar' or 'obsvar'");
    }

    if (K->cross == NULL) {
	/* not allocated yet: this should be the case only 
	   on initial set-up */
	err = kalman_add_crossinfo(K);
    }

    if (!err) {
	err = kalman_update_crossinfo(K, UPDATE_INIT);
    }

    if (err) {
	return err;
    }

    /* K->Q and K->R might not be user-supplied matrices: they could
       be matrices constructed from scalars and "owned" by the filter.
       In that case we should free them at this point in order to
       avoid leaking memory.
    */

    if (kalman_owns_matrix(K, K_Q)) {
	gretl_matrix_free((gretl_matrix *) K->Q);
    }

    if (kalman_owns_matrix(K, K_R)) {
	gretl_matrix_free((gretl_matrix *) K->R);
    }    

    K->Q = K->cross->BB;
    K->R = K->cross->CC;

#if 0
    gretl_matrix_print(K->cross->BB, "BB");
    gretl_matrix_print(K->cross->CC, "CC");
    gretl_matrix_print(K->cross->BC, "BC");
#endif

    return 0;
}

static int user_kalman_setup (kalman *K, gretlopt opt)
{
    int err = 0;

    if (K->F == NULL || K->H == NULL || K->Q == NULL || K->y == NULL) {
	return missing_matrix_error(NULL);
    }

    kalman_set_dimensions(K, opt);

    if (K->p > 0) {
	err = kalman_revise_variance(K);
	if (err) {
	    fprintf(stderr, "failed in kalman_revise_variance\n");
	}
    }

    if (!err) {
	err = kalman_check_dimensions(K);
	if (err) {
	    fprintf(stderr, "failed in kalman_check_dimensions\n");
	}
    }

    if (err) {
	return err;
    }

    err = kalman_init(K);

    if (!err) {
	gretl_matrix_zero(K->e);
	if (opt & OPT_C) {
	    K->flags |= KALMAN_CROSS;
	}
	if (opt & OPT_D) {
	    K->flags |= KALMAN_DIFFUSE;
	}
    }

    return err;
}

static int matrix_diff (const gretl_matrix *a, const gretl_matrix *b)
{
    int i, n = a->rows * a->cols;

    /* note: we ignore the 0,0 element here */

    for (i=1; i<n; i++) {
	if (b->val[i] != a->val[i]) {
	    return 1;
	}
    }

    return 0;
}

/* below: if postmult is non-zero, we're post-multiplying by the
   transpose of F */

static int multiply_by_F (kalman *K, const gretl_matrix *A, 
			  gretl_matrix *B, int postmult)
{
    int ret = 0;

#if 0 /* uninitialized stuff here? */
    gretl_matrix_print(A, "A, in multiply_by_F");
#endif

    if (gretl_is_zero_matrix(A)) {
	gretl_matrix_zero(B);
    } else if (K->nonshift == K->r) {
	if (postmult) {
	    gretl_matrix_multiply_mod(A, GRETL_MOD_NONE,
				      K->F, GRETL_MOD_TRANSPOSE,
				      B, GRETL_MOD_NONE);
	} else {
	    gretl_matrix_multiply(K->F, A, B);
	}
    } else { 
	gretl_matrix *topF = K->Tmprr_2a;
	gretl_matrix *top = K->Tmprr_2b;
	int r1 = K->nonshift;
	int r2 = K->r - r1;
	int i, j, c;
	double x;

	if (postmult) {
	    /* if post-multiplying by F', "top" actually means "left" */

	    topF->rows = K->r;
	    topF->cols = r1;
	    c = A->rows;
	    top->rows = c;
	    top->cols = r1;

	    /* copy from F to topF */
	    for (i=0; i<r1; i++) {
		for (j=0; j<K->r; j++) {
		    x = gretl_matrix_get(K->F, i, j);
		    gretl_matrix_set(topF, j, i, x);
		}
	    }

	    gretl_matrix_multiply(A, topF, top);

	    for (i=0; i<c; i++) {
		for (j=0; j<r1; j++) {
		    x = gretl_matrix_get(top, i, j);
		    gretl_matrix_set(B, i, j, x);
		}
	    }

	    for (i=0; i<c; i++) {
		for (j=0; j<r2; j++) {
		    x = gretl_matrix_get(A, i, j);
		    gretl_matrix_set(B, i, j + r1, x);
		}
	    }

	} else {
	    /* pre-multiplying by F */

	    topF->rows = r1;
	    topF->cols = K->r;
	    c = A->cols;
	    top->rows = r1;
	    top->cols = c;

	    for (i=0; i<r1; i++) {
		for (j=0; j<K->r; j++) {
		    x = gretl_matrix_get(K->F, i, j);
		    gretl_matrix_set(topF, i, j, x);
		}
	    }

	    gretl_matrix_multiply(topF, A, top);

	    for (i=0; i<r1; i++) {
		for (j=0; j<c; j++) {
		    x = gretl_matrix_get(top, i, j);
		    gretl_matrix_set(B, i, j, x);
		}
	    }

	    for (i=0; i<r2; i++) {
		for (j=0; j<c; j++) {
		    x = gretl_matrix_get(A, i, j);
		    gretl_matrix_set(B, i + r1, j, x);
		}
	    }
	}
    }

    return ret;
}

/* Simplified version of the function below, for ARMA:
   in this case we have H = (r x 1) and S = (r x 1),
   although F and P are (r x r).
*/

static int kalman_arma_iter_1 (kalman *K, int missobs)
{
    double Ve;
    int i, err = 0;

    /* write F*S into S+ */
    err += multiply_by_F(K, K->S0, K->S1, 0);

    if (missobs) {
	return err;
    }

    /* form e = y - A'x - H'S */
    K->e->val[0] -= K->Ax->val[0];
    for (i=0; i<K->r; i++) {
	K->e->val[0] -= K->H->val[i] * K->S0->val[i];
    }

    /* form FPH */
    err += multiply_by_F(K, K->PH, K->FPH, 0);
   
    /* form (H'PH + R)^{-1} * (y - Ax - H'S) = "Ve" */
    Ve = K->Vt->val[0] * K->e->val[0];

    /* form FPH * (H'PH + R)^{-1} * (y - A'x - H'S),
       and complete calculation of S+ */
    for (i=0; i<K->r; i++) {
	K->S1->val[i] += K->FPH->val[i] * Ve;
    }

    if (!err) {
	/* scalar quadratic form */
	K->SSRw += Ve * K->e->val[0];
    }

    if (K->flags & KALMAN_ETT) {
	K->e->val[0] = Ve;
    }

    return err;
}

/* Hamilton (1994) equation [13.2.23] page 381, in simplified notation:

   S+ = FS + FPH(H'PH + R)^{-1} * (y - A'x - H'S) 

   "S" is Hamilton's \xi (state vector)
*/

static int kalman_iter_1 (kalman *K, int missobs, double *llt)
{
    int err = 0;

    /* write F*S into S+ */
    err += multiply_by_F(K, K->S0, K->S1, 0);

    /* add \mu if present */
    if (K->mu != NULL) {
	gretl_matrix_add_to(K->S1, K->mu);
    }

    if (missobs) {
	/* the observable is not in fact observed at t */
	*llt = 0.0;
	if (K->K != NULL) {
	    set_row(K->K, K->t, 0.0);
	}
	return err;
    }   

    /* form e = y - A'x - H'S (e is already initialized to y) */
    err += gretl_matrix_subtract_from(K->e, K->Ax);
    err += gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
				     K->S0, GRETL_MOD_NONE,
				     K->e, GRETL_MOD_DECREMENT);

    if (!err) {
	/* contribution to log-likelihood -- see Hamilton
	   (1994) equation [13.4.1] page 385.
	*/
	double x = gretl_scalar_qform(K->e, K->Vt, &err);

	if (err) {
	    *llt = NADBL;
	} else {
	    *llt -= 0.5 * x;
	    K->SSRw += x;
	}
    }

    /* form the gain, Kt = (FPH + BC') * (H'PH + R)^{-1} */
    err += multiply_by_F(K, K->PH, K->FPH, 0);
    if (K->p > 0) {
	/* cross-correlated case */
	gretl_matrix_add_to(K->FPH, K->cross->BC);
    }
    err += gretl_matrix_multiply(K->FPH, K->Vt, K->Kt);

    /* form K_t * e_t and add to S+ */
    err += gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
				     K->e,  GRETL_MOD_NONE,
				     K->S1, GRETL_MOD_CUMULATE);

    if (!err && K->K != NULL) {
	/* record the gain */
	load_to_vec(K->K, K->Kt, K->t);
    }

    return err;
}

/* Hamilton (1994) equation [13.2.22] page 380, in simplified notation:

      P+ = F[P - PH(H'PH + R)^{-1}H'P]F' + Q 

   Or use the alternative formulation:

      P+ = FPF' - KVK + Q  where K = gain and V = error variance
*/

static int kalman_iter_2 (kalman *K, int missobs)
{
    int err = 0;

    if (K->p == 0 && !missobs) {
	/* revise "P0" as P_{t|t} = P - PH(H'PH + R)^{-1}H'P */
	err = gretl_matrix_qform(K->PH, GRETL_MOD_NONE, K->Vt,
				 K->P0, GRETL_MOD_DECREMENT);
    }

    /* pre-multiply by F, post-multiply by F' */
    if (K->nonshift == K->r) {
	err += gretl_matrix_qform(K->F, GRETL_MOD_NONE, K->P0,
				  K->P1, GRETL_MOD_NONE);
    } else {
	err += multiply_by_F(K, K->P0, K->Tmprr, 0);
	err += multiply_by_F(K, K->Tmprr, K->P1, 1);
	gretl_matrix_xtr_symmetric(K->P1);
    }

    if (K->p > 0 && !missobs) {
	/* subtract K(H'PH + R)K' */
	gretl_matrix_qform(K->Kt, GRETL_MOD_NONE, K->HPH, 
			   K->P1, GRETL_MOD_DECREMENT);
    }

    /* add Q */
    err += gretl_matrix_add_to(K->P1, K->Q);

    return err;
}

#if KDEBUG > 1
static void kalman_print_state (kalman *K)
{
    int j;

    /* if (t > 5) return; */

    fprintf(stderr, "t = %d:\n", K->t);

    for (j=0; j<K->n; j++) {
	fprintf(stderr, "y[%d] = %.10g, err[%d] = %.10g\n", j, 
		gretl_matrix_get(K->y, K->t, j), 
		j, gretl_vector_get(K->e, j));
    }

    gretl_matrix_print(K->S0, "K->S0");
    gretl_matrix_print(K->P0, "K->P0");
}
#endif

static void kalman_record_state (kalman *K)
{
    if (K->S != NULL) {
	load_to_row(K->S, K->S0, K->t);
    }

    if (K->P != NULL) {
	load_to_vech(K->P, K->P0, K->r, K->t);
    }
}

void kalman_set_nonshift (kalman *K, int n)
{
    K->nonshift = n;
}

void kalman_set_options (kalman *K, int opts)
{
    K->flags |= opts;
}

int kalman_get_options (kalman *K)
{
    return K->flags;
}

/* Read from the appropriate row of x (T x k) and multiply by A' to
   form A'x_t.  Note: the flag K->ifc is used to indicate that the
   observation equation has an implicit constant, with an entry in
   the A matrix (the first) but no explicit entry in the x matrix.
*/

static void kalman_set_Ax (kalman *K, int *missobs)
{
    double xjt, axi;
    int i, j;

    for (i=0; i<K->n; i++) {
	axi = 0.0;
	for (j=0; j<K->k; j++) {
	    if (K->ifc) {
		xjt = (j == 0)? 1.0 : gretl_matrix_get(K->x, K->t, j - 1);
	    } else {
		xjt = gretl_matrix_get(K->x, K->t, j);
	    }
	    if (isnan(xjt)) {
#if KDEBUG
		fprintf(stderr, "x: got nan at obs %d\n", K->t);
#endif
		*missobs += 1;
	    }
	    axi += xjt * gretl_matrix_get(K->A, j, i);
	}
	gretl_vector_set(K->Ax, i, axi);
    }
}

/* read from the appropriate row of y (T x n) and transcribe to
   the current e (n x 1)
*/

static void kalman_initialize_error (kalman *K, int *missobs)
{
    double yti;
    int i;

    for (i=0; i<K->n; i++) {
	yti = gretl_matrix_get(K->y, K->t, i);
	K->e->val[i] = yti;
	if (isnan(yti)) {
#if KDEBUG
	    fprintf(stderr, "y: got nan at obs %d\n", K->t);
#endif
	    *missobs += 1;
	}
    }
}

/* The user gave a function call that should be used for updating
   a given Kalman coefficient matrix: here we call it.
*/

static gretl_matrix *kalman_update_matrix (kalman *K, int i, 
					   PRN *prn, int *err)
{
    gretl_matrix *m = NULL;

    *err = generate(K->matcalls[i], NULL, OPT_O, prn);

    if (*err) {
	fprintf(stderr, "kalman_update_matrix: call='%s', err=%d\n", 
		K->matcalls[i], *err);
    } else {
	m = get_matrix_by_name(K->mnames[i]);
	if (m == NULL) {
	    *err = E_DATA;
	}
    }

    return m;
}

/* If we have any time-varying coefficient matrices, refresh these for
   the current time step.  This is called on a forward filtering pass.
*/

static int kalman_refresh_matrices (kalman *K, PRN *prn)
{
    const gretl_matrix **cptr[] = {
	&K->F, &K->A, &K->H, &K->Q, &K->R, &K->mu
    };  
    int cross_update = 0;
    int i, err = 0;

    for (i=0; i<NMATCALLS && !err; i++) {
	if (matrix_is_varying(K, i)) {
	    *cptr[i] = kalman_update_matrix(K, i, prn, &err);
	    if (!err) {
		if (K->p > 0 && i >= K_Q) {
		    cross_update = 1;
		} else {
		    err = check_matrix_dims(K, *cptr[i], i);
		}
	    }	    
	    if (err) {
		fprintf(stderr, "kalman_refresh_matrices: err = %d at t = %d\n", 
			err, K->t);
	    } 
	}
    }

    if (!err && K->step != NULL) {
	/* keep a record of F and/or H at the given time step */
	if (K->step->F != NULL) {
	    load_to_vec(K->step->F, K->F, K->t); 
	} 
	if (K->step->H != NULL) {
	    load_to_vec(K->step->H, K->H, K->t); 
	} 	
    }

    if (!err && cross_update) {
	/* cross-correlated case */
	err = kalman_update_crossinfo(K, UPDATE_STEP);
    }

    return err;
}

/**
 * kalman_forecast:
 * @K: pointer to Kalman struct: see kalman_new().
 * @prn: printing apparatus (or NULL).
 *
 * Generates a series of one-step ahead forecasts for y, based on
 * information entered initially using kalman_new(), and possibly
 * modified using kalman_set_initial_state_vector() and/or
 * kalman_set_initial_MSE_matrix().  The log-likelihood is
 * calculated for the sequence of forecast errors on the assumption
 * of normality: this can be accessed using kalman_get_loglik().
 *
 * Returns: 0 on success, non-zero on error.
 */

int kalman_forecast (kalman *K, PRN *prn)
{
    double ldet;
    int smoothing, update_P = 1;
    int Tmiss = 0;
    int i, err = 0;

#if KDEBUG
    fprintf(stderr, "kalman_forecast: T = %d\n", K->T);
#endif 

    smoothing = (K->flags & KALMAN_SMOOTH)? 1 : 0;

    if (K->nonshift < 0) {
	K->nonshift = count_nonshifts(K->F);
    }

    K->SSRw = K->sumldet = K->loglik = 0.0;
    K->s2 = NADBL;
    K->okT = K->T;

    if (K->x == NULL) {
	/* no exogenous vars */
	if (K->A != NULL) {
	    /* implicit const case: A is 1 x n and A'x is n x 1 */
	    gretl_vector_copy_values(K->Ax, K->A);
	} else {
	    gretl_matrix_zero(K->Ax);
	}
    } 

    set_kalman_running(K);

    for (K->t = 0; K->t < K->T && !err; K->t += 1) {
	int missobs = 0;
	double llt = 0.0;

#if KDEBUG > 1
	kalman_print_state(K);
#endif

	if (K->S != NULL || K->P != NULL) {
	    kalman_record_state(K);
	}

	if (filter_is_varying(K)) {
	    /* we have time-varying coefficients */
	    err = kalman_refresh_matrices(K, prn);
	    if (err) {
		K->loglik = NADBL;
		break;
	    }
	}

	/* read slice from y */
	kalman_initialize_error(K, &missobs);

	if (K->x != NULL) {
	    /* read from x if applicable 
	       note 2010-09-18: this was conditional on missobs < K->n
	       FIXME?
	    */
	    kalman_set_Ax(K, &missobs);
	}

	if (missobs) {
	    /* 2010-09-18: this was right after kalman_initialize_error 
	       FIXME?
	     */
	    Tmiss++;
	}	

	/* initial matrix calculations: form PH and H'PH 
	   (note that we need PH later) */
	gretl_matrix_multiply(K->P0, K->H, K->PH);
	if (K->n == 1) {
	    /* slight speed-up for univariate observable */
	    double x = (K->R == NULL)? 0.0 : K->R->val[0];

	    for (i=0; i<K->r; i++) {
		x += K->H->val[i] * K->PH->val[i];
	    }
	    if (x <= 0.0) {
		err = E_NAN;
	    } else {
		K->HPH->val[0] = x;
		ldet = log(x);
		K->Vt->val[0] = 1.0 / x;
	    }
	} else {
	    gretl_matrix_qform(K->H, GRETL_MOD_TRANSPOSE,
			       K->P0, K->HPH, GRETL_MOD_NONE);
	    if (K->R != NULL) {
		gretl_matrix_add_to(K->HPH, K->R);
	    }
	    gretl_matrix_copy_values(K->Vt, K->HPH);
	    err = gretl_invert_symmetric_matrix2(K->Vt, &ldet);
	    if (err) {
		fprintf(stderr, "kalman_forecast: failed to invert V\n");
	    }
	}

	/* likelihood bookkeeping */
	if (err) {
	    K->loglik = NADBL;
	    break;
	} else if (!missobs) {
	    K->sumldet += ldet;
	} 

	if (K->V != NULL) {
	    /* record MSE for estimate of observables */
	    if (missobs) {
		if (smoothing) {
		    set_row(K->V, K->t, 0.0); 
		} else {
		    set_row(K->V, K->t, 1.0/0.0); 
		}
	    } else if (smoothing) {
		/* record inverse */
		load_to_vech(K->V, K->Vt, K->n, K->t);
	    } else {
		/* record uninverted matrix */
		load_to_vech(K->V, K->HPH, K->n, K->t);
	    }
	}

	/* first stage of dual iteration */
	if (arma_ll(K) && !smoothing) {
	    err = kalman_arma_iter_1(K, missobs);
	} else {
	    err = kalman_iter_1(K, missobs, &llt);
	    if (K->LL != NULL) {
		if (na(llt) || missobs) {
		    llt = M_NA;
		} else {
		    llt -= 0.5 * (K->n * LN_2_PI + ldet);
		}
		gretl_vector_set(K->LL, K->t, llt);
	    }
	}
	
	/* record forecast errors if wanted */
	if (!err && K->E != NULL) {
	    if (missobs && smoothing) {
		set_row(K->E, K->t, 0.0);
	    } else {
		load_to_row(K->E, K->e, K->t);
	    }
	}

	/* update state vector */
	if (!err) {
	    gretl_matrix_copy_values(K->S0, K->S1);
	}

	if (!err && update_P) {
	    /* second stage of dual iteration */
	    err = kalman_iter_2(K, missobs);
	}

	if (!err) {
	    /* update MSE matrix, if needed */
	    if (arma_ll(K) && !smoothing && update_P && K->t > 20) {
		if (!matrix_diff(K->P1, K->P0)) {
		    K->P0->val[0] += 1.0;
		    update_P = 0;
		}
	    }
	    if (update_P) {
		gretl_matrix_copy_values(K->P0, K->P1);
	    } 
	}
    }

    set_kalman_stopped(K);

    if (isnan(K->loglik) || isinf(K->loglik)) {
	K->loglik = NADBL;
    }

    if (na(K->loglik)) {
	goto bailout;
    }

    K->okT = K->T - Tmiss;

    if (arma_ll(K)) {
	double ll1 = 1.0 + LN_2_PI + log(K->SSRw / K->okT);

	if (K->flags & KALMAN_AVG_LL) {
	    K->loglik = -0.5 * (ll1 + K->sumldet / K->okT);
	} else {
	    K->loglik = -0.5 * (K->okT * ll1 + K->sumldet);
	}
    } else {
	/* For K->s2 see Koopman, Shephard and Doornik, "Statistical
	   algorithms for models in state space using SsfPack 2.2",
	   Econometrics Journal, 1999, vol. 2, pp. 113-166; also
	   available at http://www.ssfpack.com/ .  For the role of
	   'd', see in addition Koopman's 1997 JASA article.
	*/
	int nT = K->n * K->okT;
	int d = (K->flags & KALMAN_DIFFUSE)? K->r : 0;

	if (d > 0) {
	    K->loglik = -0.5 * ((nT - d) * LN_2_PI + K->sumldet + K->SSRw)
		+ (d / 2.0) * log(kappa);
	} else {
	    K->loglik = -0.5 * (nT * LN_2_PI + K->sumldet + K->SSRw);
	}
	K->s2 = K->SSRw / (nT - d);
    }

    if (isnan(K->loglik) || isinf(K->loglik)) {
	K->loglik = NADBL;
	err = E_NAN;
    }     

    bailout:

#if KDEBUG
    fprintf(stderr, "kalman_forecast: err=%d, ll=%#.12g\n", err, 
	    K->loglik);
#endif

    return err;
}

/**
 * kalman_get_loglik:
 * @K: pointer to Kalman struct.
 * 
 * Retrieves the log-likelhood calculated via a run of 
 * kalman_forecast().
 * 
 * Returns: ll value, or #NADBL on failure.
 */

double kalman_get_loglik (const kalman *K)
{
    return K->loglik;
}

/**
 * kalman_get_arma_variance:
 * @K: pointer to Kalman struct.
 * 
 * Retrieves the estimated variance for an ARMA model
 * estimated using the Kalman filter.
 * 
 * Returns: sigma-squared value, or #NADBL on failure.
 */

double kalman_get_arma_variance (const kalman *K)
{
    if (na(K->SSRw)) {
	return NADBL;
    } else {
	return K->SSRw / K->okT;
    }
}

/**
 * kalman_set_initial_state_vector:
 * @K: pointer to Kalman struct.
 * @S: matrix of values to set.
 * 
 * Resets the initial value of the state vector in a Kalman
 * struct, using the values from @S.  See also kalman_new().
 * 
 * Returns: 0 on success, non-zero on error.
 */

int kalman_set_initial_state_vector (kalman *K, const gretl_matrix *S)
{
    return gretl_matrix_copy_values(K->S0, S);
}

/**
 * kalman_set_initial_MSE_matrix:
 * @K: pointer to Kalman struct.
 * @P: matrix of values to set.
 * 
 * Resets the initial value of the MSE matrix in a Kalman
 * struct, using the values from @P.  See also kalman_new().
 * 
 * Returns: 0 on success, non-zero on error.
 */

int kalman_set_initial_MSE_matrix (kalman *K, const gretl_matrix *P)
{
    return gretl_matrix_copy_values(K->P0, P);
}

/* apparatus for user-defined Kalman filters */

typedef struct user_kalman_ user_kalman;

struct user_kalman_ {
    kalman *K;
    int fnlevel;
};

/* At present there can be at most one user_kalman at each
   depth of function execution */

static user_kalman **uK; /* user-defined Kalman structs */
static int n_uK;         /* number of same */

struct K_input_mat {
    int sym;
    const char *name;
};

/* mapping to names used in "kalman ... end kalman" block */

struct K_input_mat K_input_mats[] = {
    { K_y, "obsy" },
    { K_H, "obsymat" },
    { K_x, "obsx" },
    { K_A, "obsxmat" },
    { K_R, "obsvar" },
    { K_F, "statemat" },
    { K_Q, "statevar" },
    { K_m, "stconst" },
    { K_S, "inistate" },
    { K_P, "inivar" }
};

/* Add storage to record function calls for updating matrices.
   We do this if at least one coefficient matrix has a time-varying
   specification.  We size the array to the max number of time-varying
   matrices, but some slots may be left NULL if the corresponding
   matrix is not in fact time-varying.
*/

static int kalman_add_matcalls (kalman *K)
{
    K->matcalls = strings_array_new(NMATCALLS);

    if (K->matcalls == NULL) {
	return E_ALLOC;
    } else {
	return 0;
    }
}

/* We found a function call given by way of specification of one
   of the coefficient matrices in a Kalman filter.
*/

static int add_matrix_fncall (kalman *K, const char *s, int i,
			      char *mname, gretl_matrix **pm)
{
    char *p, *q, *fncall;
    int n, err = 0;

    s += strspn(s, " ");

    /* We need the name of a matrix followed by a function call, the
       two elements separated by a semicolon: the separating ';' must
       come before the opening parenthesis of the function call.
    */

    p = strchr(s, ';');
    q = strchr(s, '(');

    if (p != NULL && q - p > 0) {
	/* OK, we have the two separated elements */
	n = strcspn(s, " ;");
	if (n > VNAMELEN - 1) {
	    n = VNAMELEN - 1;
	}
	strncat(mname, s, n);
	s = p + 1;
	s += strspn(s, " ");
    } else {
	err = E_PARSE;
    }

    if (!err && K->matcalls == NULL) {
	err = kalman_add_matcalls(K);
    }

    if (!err) {
	fncall = gretl_strdup(s);
	if (fncall == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	return err;
    }

    tailstrip(fncall); 

    *pm = get_matrix_by_name(mname);
    if (*pm == NULL) {
	err = E_UNKVAR;
    }

    if (err) {
	free(fncall);
    } else {
	/* record the function call */
	K->matcalls[i] = fncall;
    }

    return err;
}

/* We didn't find a matrix of the name @s: here we try for a series
   or a named list, and if found, make a matrix out of it.  This
   matrix will be "hard-wired" -- we don't do any further look-up
   of values -- and it will be owned by the Kalman struct.
*/

static gretl_matrix *k_matrix_from_dataset (char *s,
					    const DATASET *dset,
					    int *err)
{
    int v = current_series_index(dset, s);
    gretl_matrix *m = NULL;
	
    if (v >= 0) {
	int L[2] = {1, v};

	m = gretl_matrix_data_subset(L, dset, dset->t1, dset->t2,
				     M_MISSING_OK, err);
    } else {
	int *list = get_list_by_name(s);

	if (list != NULL) {
	    m = gretl_matrix_data_subset(list, dset,
					 dset->t1, dset->t2,
					 M_MISSING_OK, err);
	} 
    }

    if (m != NULL) {
	strcpy(s, GENERIC_MATNAME);
    }

    return m;
}

/* We didn't find a matrix of the name @s: here we try for a scalar
   (either a named scalar variable or a numeric constant), and if
   found, make a matrix out of it.  If we get a named scalar, we'll
   record its name so we're able to update its value later if need be.
*/

static gretl_matrix *k_matrix_from_scalar (char *s, int *errp)
{
    gretl_matrix *m = NULL;
    double x;
    int err = 0;

    x = gretl_double_from_string(s, &err);

    if (!err && na(x)) {
	*errp = err = E_MISSDATA;
    }

    if (!err) {
	m = gretl_matrix_from_scalar(x);
	if (m == NULL) {
	    *errp = err = E_ALLOC;
	}
    }

    if (!err) {
	if (gretl_is_scalar(s)) {
	    /* got a scalar variable: record its name */
	    char tmp[VNAMELEN];

	    strcpy(tmp, s);
	    sprintf(s, "$%s", tmp);
	} else {
	    strcpy(s, GENERIC_MATNAME);
	}
    }

    return m;
}

/* attach_input_matrix: the string @s should contain:

   (a) the name of a user-defined variable, alone (either a matrix or,
   if the dimensions are OK, a series, named list or scalar), or

   (b) a function call that returns a matrix, or 

   (c) the name of a matrix plus a void function that updates that
   matrix.

   Options (b) and (c) are designed to support time-varying 
   coefficient matrices.

   The integer @i is an ID number that identifies the role of the
   matrix in question within the Kalman struct.  Given this info
   we check for errors and, if all is OK, hook things up.
*/

static int 
attach_input_matrix (kalman *K, const char *s, int i,
		     const DATASET *dset)
{
    char mname[VNAMELEN+1] = {0};
    gretl_matrix *m = NULL;
    int err = 0;

    if (i < 0 || i >= K_MMAX) {
	/* "can't happen" */
	return E_DATA;
    }

    if (K->mnames == NULL) {
	K->mnames = strings_array_new_with_length(K_MMAX, VNAMELEN+1);
	if (K->mnames == NULL) {
	    return E_ALLOC;
	}
    }

    if (i <= K_R && strchr(s, '(')) {
	/* we have a function call? */
	err = add_matrix_fncall(K, s, i, mname, &m);
    } else if (gretl_scan_varname(s, mname) == 1) {
	m = get_matrix_by_name(mname);
    } else {
	err = E_PARSE;
    }

    if (!err && m == NULL) {
	/* didn't find a matrix */
	if (i == K_y || i == K_x) {
	    /* if we're looking at osby or obsx, try for a series or
	       list; this should fail gracefully if there's no dataset
	       in place
	    */
	    m = k_matrix_from_dataset(mname, dset, &err);
	} else {
	    /* parameter matrix: try a scalar */
	    m = k_matrix_from_scalar(mname, &err);
	}
    }

    if (!err && m == NULL) {
	/* out of options now */
	gretl_errmsg_sprintf(_("'%s': no such matrix"), mname);
	err = E_UNKVAR;
    }

    if (!err) {
	if (i == K_F) {
	    K->F = m;
	} else if (i == K_A) {
	    K->A = m;
	} else if (i == K_H) {
	    K->H = m;
	} else if (i == K_Q) {
	    K->Q = m;
	} else if (i == K_R) {
	    K->R = m;
	} else if (i == K_m) {
	    K->mu = m;
	} else if (i == K_y) {
	    K->y = m;
	} else if (i == K_x) {
	    K->x = m;
	} else if (i == K_S) {
	    K->Sini = m;
	} else if (i == K_P) {
	    K->Pini = m;
	} 

	/* record name of matrix */
	strcpy(K->mnames[i], mname);
    }

    return err;
}

static gretl_matrix *kalman_retrieve_matrix (const char *name,
					     int level, int cfd)
{
    gretl_matrix *m;

    if (*name == '\0') {
	return NULL;
    }

    /* We allow for the possibility that one or more of the matrices
       that are attached to a kalman struct have been temporarily
       "promoted", via matrix-pointer arguments, to the level of a
       user-defined function.  This can happen if a user-defined
       function is using a Kalman filter defined by the caller.
    */

    m = get_matrix_by_name_at_level(name, level);
    if (m == NULL && level < cfd) {
	m = get_matrix_by_name_at_level(name, cfd);
    }

    return m;
}

static int obsy_error (kalman *K)
{
    if (K->y == NULL) {
	return missing_matrix_error("obsy");
    } else if (K->y->rows != K->T || K->y->cols != K->n) {
	return E_NONCONF;
    } else {
	return 0;
    }
}

/* Did the user give the name of a scalar variable in place of a
   1 x 1 Kalman input matrix?  If so, re-read the value of that 
   scalar */

static int update_scalar_matrix (gretl_matrix *m, const char *name)
{
    int err = 0;

    if (gretl_matrix_is_scalar(m) && strcmp(name, GENERIC_MATNAME)) {
	double x = gretl_scalar_get_value(name + 1, NULL);

	if (na(x)) {
	    /* or should we just ignore x in this case? */
	    return E_MISSDATA;
	} else {
	    m->val[0] = x;
	}
    }

    return err;
}

/* When (re-)running a user-defined filter, check that no relevant
   matrices have disappeared or been resized.  In addition,
   re-initialize the state and variance.
*/

static int user_kalman_recheck_matrices (user_kalman *u, PRN *prn)
{
    int cfd = gretl_function_depth();
    kalman *K = u->K;
    const gretl_matrix **cptr[] = {
	&K->F, &K->A, &K->H, &K->Q, &K->R,
	&K->mu, &K->y, &K->x, &K->Sini, &K->Pini
    };
    int i, err = 0;

    for (i=0; i<K_MMAX; i++) {
	if (!kalman_owns_matrix(K, i)) {
	    /* pointer may be invalid, needs refreshing */
	    *cptr[i] = NULL;
	}
    }

    K->flags |= KALMAN_CHECK;

    for (i=0; i<K_MMAX && !err; i++) {
	if (i <= K_m && matrix_is_varying(K, i)) {
	    *cptr[i] = kalman_update_matrix(K, i, prn, &err);
	} else if (kalman_owns_matrix(K, i)) {
	    err = update_scalar_matrix(*(gretl_matrix **) cptr[i], K->mnames[i]);
	} else {
	    *cptr[i] = kalman_retrieve_matrix(K->mnames[i], u->fnlevel, cfd);
	}
    }

    K->flags ^= KALMAN_CHECK;

    if (err) {
	return err;
    }

    if (K->F == NULL || K->H == NULL || K->Q == NULL) {
	err = missing_matrix_error(NULL);
    } else if (gretl_matrix_rows(K->F) != K->r ||
	       gretl_matrix_rows(K->A) != K->k) {
	err = E_NONCONF;
    } else if (!kalman_simulating(K)) {
	err = obsy_error(K);
    }

    if (!err && K->p > 0) {
	err = kalman_revise_variance(K);
    }

    if (!err) {
	err = kalman_check_dimensions(K);
    }

    if (!err) {
	if (K->Sini != NULL) {
	    gretl_matrix_copy_values(K->S0, K->Sini);
	} else {
	    gretl_matrix_zero(K->S0);
	}
	if (K->Pini != NULL) {
	    gretl_matrix_copy_values(K->P0, K->Pini);
	} else {
	    err = construct_Pini(K);
	}
    }
    
    return err;
}

/* Try to find a user-defined Kalman filter at the current level of
   function execution.  Failing that, if @level is given as -1, we
   also try looking for a filter available at parent levels of 
   execution. 
*/

static user_kalman *get_user_kalman (int level)
{
    int i, go_up = 0;

    if (level == -1) {
	level = gretl_function_depth();
	if (level > 0) {
	    go_up = 1;
	}
    }

    for (i=0; i<n_uK; i++) {
	if (uK[i] != NULL && uK[i]->fnlevel == level) {
	    return uK[i];
	} 
    }

    if (go_up) {
	while (--level >= 0) {
	    for (i=0; i<n_uK; i++) {
		if (uK[i] != NULL && uK[i]->fnlevel == level) {
		    return uK[i];
		} 
	    }
	}
    }

    return NULL;
}

static int real_destroy_user_kalman (PRN *prn)
{
    int i, lev = gretl_function_depth();

    for (i=0; i<n_uK; i++) {
	if (uK[i] != NULL && uK[i]->fnlevel == lev) {
	    kalman_free(uK[i]->K);
	    free(uK[i]);
	    uK[i] = NULL;
	    if (prn != NULL && gretl_messages_on()) {
		pputs(prn, "Deleted kalman filter\n");
	    }
	    return 0;
	}
    }

    return E_UNKVAR;
}

static void destroy_user_kalman (void)
{
    real_destroy_user_kalman(NULL);
}

void kalman_cleanup (void)
{
    int i;

    for (i=0; i<n_uK; i++) {
	if (uK[i] != NULL) {
	    kalman_free(uK[i]->K);
	    free(uK[i]);
	}
    }

    free(uK);
    n_uK = 0;
}

int delete_kalman (PRN *prn)
{
    return real_destroy_user_kalman(prn);
}

void kalman_attach_data (kalman *K, void *data)
{
    if (K != NULL) {
	K->data = data;
    }
}

void *kalman_get_data (const kalman *K)
{
    return (K != NULL)? K->data : NULL;
}

void kalman_attach_printer (kalman *K, PRN *prn)
{
    if (K != NULL) {
	K->prn = prn;
    }
}

PRN *kalman_get_printer (const kalman *K)
{
    return (K != NULL)? K->prn : NULL;
}

static int add_user_kalman (void)
{
    user_kalman *u = NULL;
    int i, n = -1;
    int err = 0;

    /* clean out any existing entry */
    destroy_user_kalman();

    for (i=0; i<n_uK; i++) {
	/* do we have a blank slot? */
	if (uK[i] == NULL) {
	    /* yes: record the slot number */
	    n = i;
	    break;
	}
    }

    if (n < 0) {
	/* no: have to allocate space */
	user_kalman **tmp;

	tmp = realloc(uK, (n_uK + 1) * sizeof *uK);
	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    uK = tmp;
	    n = n_uK;
	    uK[n] = NULL;
	    n_uK++;
	}
    }

    if (!err) {
	u = malloc(sizeof *u);
	if (u == NULL) {
	    err = E_ALLOC;
	} else {
	    u->K = kalman_new_empty(KALMAN_USER);
	    if (u->K == NULL) {
		free(u);
		err = E_ALLOC;
	    } else {
		/* finalize and place on stack */
		u->fnlevel = gretl_function_depth();
		u->K->fnlevel = u->fnlevel;
		uK[n] = u;
	    }
	}
    }
	
    return err;
}

/* check the content of @ps for a recognized matrix code-name, and if
   successful, move the char pointer beyond the code-name
*/

static int get_kalman_matrix_id (const char **ps)
{
    const char *test;
    int i, n;

    for (i=0; i<K_MMAX; i++) {
	test = K_input_mats[i].name;
	n = strlen(test);
	if (!strncmp(*ps, test, n) && (*ps)[n] == ' ') {
	    *ps += n + 1;
	    return K_input_mats[i].sym;
	}
    }

    /* failed */
    gretl_errmsg_sprintf("kalman: invalid input '%s'", *ps);

    return -1;
}

static const char *kalman_matrix_name (int sym)
{
    int i;

    for (i=0; i<K_MMAX; i++) {
	if (K_input_mats[i].sym == sym) {
	    return K_input_mats[i].name;
	}
    }

    /* failed */
    return "matrix";
}

/*
 * kalman_parse_line:
 * @line: "kalman" to start, "end kalman" to end; otherwise
 * this string should contain a matrix specification on the 
 * pattern "key value".
 * @dset: dataset struct (may be NULL).
 * @opt: may contain %OPT_D for diffuse initialization of the
 * Kalman filter, %OPT_C to specify that the disturbances are 
 * correlated across the two equations.
 *
 * Parses @line and either (a) starts a filter definition or
 * (b) adds a matrix specification to the filter or (c)
 * completes the filter set-up.
 *
 * Returns: 0 on successful completion, non-zero error code
 * otherwise.
 */

int kalman_parse_line (const char *line, const DATASET *dset, 
		       gretlopt opt)
{
    user_kalman *u;
    int err = 0;

    if (opt == OPT_NONE && !strncmp(line, "kalman", 6)) {
	/* starting: allocate and return */
	return add_user_kalman();
    }

    u = get_user_kalman(gretl_function_depth());

    if (u == NULL) {
	return missing_kalman_error();
    }    

    if (!strncmp(line, "end ", 4)) {
	/* "end kalman": complete the set-up */
	err = user_kalman_setup(u->K, opt);
    } else {
	/* we're supposed to find a matrix spec */ 
	const char *s = line;
	int m = get_kalman_matrix_id(&s);

	if (m < 0) {
	    err = E_PARSE;
	} else {
	    err = attach_input_matrix(u->K, s, m, dset);
	}
    }

    if (err) {
	fprintf(stderr, "kalman_parse_line: '%s', err = %d\n",
		line, err);
	destroy_user_kalman();
    }

    return err;
}

#if 0
static void kalman_timestamp_matrix (gretl_matrix *m,
				     const DATASET *dset)
{
    if (K->y != NULL) {
	int t1 = K->y->t1;
	int t2 = K->y->t2;

	if (t1 > 0 && t1 < dset->n &&
	    t2 > 0 && t2 <= dset->n &&
	    t2 > t1) {
	    m->t1 = t1;
	    m->t2 = t2;
	}
    }
}
#endif

static gretl_matrix *series_export_matrix (kalman *K,
					   const char *name,
					   const DATASET *dset,
					   int i, int *v,
					   int *err)
{
    gretl_matrix *m = NULL;
    int T = sample_size(dset);

    if (T == K->T) {
	if ((i == K_E || i == K_V) && K->n == 1) {
	    /* just one observable */
	    *v = current_series_index(dset, name);
	} else if ((i == K_BIG_S || i == K_BIG_P) && K->r == 1) {
	    /* just one state variable */
	    *v = current_series_index(dset, name);
	} else if (i == K_K && K->n == 1 && K->r == 1) {
	    /* observable and state are both 1-vectors */
	    *v = current_series_index(dset, name);
	}
    }

    if (*v > 0) {
	/* make temporary matrix for series export */
	m = gretl_column_vector_alloc(T);
	if (m == NULL) {
	    *err = E_ALLOC;
	} 
    }

    return m;
}

static int attach_export_matrix (kalman *K,
				 const char *name, 
				 const DATASET *dset,
				 int *smat, int i)
{
    gretl_matrix *m = NULL;
    int err = 0;

    if (name != NULL && strcmp(name, "null")) {
	m = get_matrix_by_name(name);
	if (m == NULL) {
	    int v = 0;

	    m = series_export_matrix(K, name, dset, i, &v, &err);
	    if (m != NULL) {
		smat[i] = v;
	    } else if (!err) {
		gretl_errmsg_sprintf(_("'%s': no such matrix"), name);
		err = E_UNKVAR;
	    }
	}
    }

    if (!err) {
	if (i == K_E) {
	    K->E = m;
	} else if (i == K_V) {
	    K->V = m;
	} else if (i == K_BIG_S) {
	    K->S = m;
	} else if (i == K_BIG_P) {
	    K->P = m;
	} else if (i == K_K) {
	    K->K = m;
	}
    }

    return err;
}

/* transcribe from temporary matrix to export series (if
   things went OK), then free the temporary matrix
*/

static void transcribe_and_free (double *x, gretl_matrix *m,
				 int err)
{
    if (!err) {
	/* transcribe only if there was no error */
	int t;

	for (t=0; t<m->rows; t++) {
	    x[t] = m->val[t];
	}
    }

    gretl_matrix_free(m);
}

/* Called on behalf of the kfilter() function: run a user-defined Kalman
   filter in forecasting mode. The doc for kfilter says that it returns
   0 on successful completion or 1 if numerical problems are encountered.
*/

int user_kalman_run (const char *E, const char *V, 
		     const char *S, const char *P,
		     const char *G, const DATASET *dset,
		     PRN *prn, int *errp)
{
    user_kalman *u = get_user_kalman(-1);
    int smat[K_K+1] = {0};
    kalman *K;
    int err = 0;

    if (u == NULL) {
	err = missing_kalman_error();
	return 1;
    }

    K = u->K;

    if (!err) {
	/* forecast errors */
	err = attach_export_matrix(K, E, dset, smat, K_E);
    }

    if (!err) {
	/* MSE for observables */
	err = attach_export_matrix(K, V, dset, smat, K_V);
    } 

    if (!err) {
	/* estimate of state */
	err = attach_export_matrix(K, S, dset, smat, K_BIG_S);
    }

    if (!err) {
	/* MSE of estimate of state */
	err = attach_export_matrix(K, P, dset, smat, K_BIG_P);
    } 

    if (!err) {
	/* Kalman gain */
	err = attach_export_matrix(K, G, dset, smat, K_K);
    } 

    if (!err && K->LL == NULL) {
	/* log-likelihood: available via accessor */
	K->LL = gretl_matrix_alloc(K->T, 1);
	if (K->LL == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = user_kalman_recheck_matrices(u, prn);
    }

    if (!err) {
	err = kalman_forecast(K, prn);
    }

    if (dset != NULL) {
	int t1 = dset->t1;

	if (smat[K_E] > 0) 
	    transcribe_and_free(dset->Z[smat[K_E]] + t1, K->E, err);
	if (smat[K_V > 0]) 
	    transcribe_and_free(dset->Z[smat[K_V]] + t1, K->V, err);
	if (smat[K_BIG_S] > 0) 
	    transcribe_and_free(dset->Z[smat[K_BIG_S]] + t1, K->S, err);
	if (smat[K_BIG_P] > 0) 
	    transcribe_and_free(dset->Z[smat[K_BIG_P]] + t1, K->P, err);
	if (smat[K_K] > 0) 
	    transcribe_and_free(dset->Z[smat[K_K]] + t1, K->K, err);
    }

    if (err != E_NAN) {
	*errp = err;
    } else {
	/* we'll flag E_NAN with a return value of 1 but 
	   won't count it as a 'true' error */
	*errp = 0;
    }

#if KDEBUG
    fprintf(stderr, "user_kalman_run: returning %d, *errp = %d\n",
	    err, *errp);
#endif

    /* detach matrices */
    K->E = K->V = K->S = NULL;
    K->P = K->K = NULL;

    return err;    
}

/* Copy row @t from @src into @targ; or add row @t of @src to
   @targ; or subtract row @t of @src from @targ.  We allow the
   possibility that the length of vector @targ is less than
   the number of columns in @src, but not the converse.
*/

static int load_from_row (gretl_vector *targ, const gretl_matrix *src, 
			  int t, int mod)
{
    int i, n = gretl_vector_get_length(targ);
    double x;

    if (n > src->cols) {
	fprintf(stderr, "load_from_row: targ length = %d, but src "
		"has %d columns\n", n, src->cols);
	return 1;
    }

    for (i=0; i<n; i++) {
	x = gretl_matrix_get(src, t, i);
	if (mod == GRETL_MOD_CUMULATE) {
	    targ->val[i] += x;
	} else if (mod == GRETL_MOD_DECREMENT) {
	    targ->val[i] -= x;
	} else {
	    targ->val[i] = x;
	}
    }

    return 0;
}

/* Row @t of @src represents the vech of an n x n matrix: extract the
   row and apply the inverse operation of vech to reconstitute the
   matrix in @targ -- or subtract the newly reconstituted matrix
   from @targ.
*/

static void load_from_vech (gretl_matrix *targ, const gretl_matrix *src,
			    int n, int t, int mod)
{
    int i, j, m = 0;
    double x;

    for (i=0; i<n; i++) {
	for (j=i; j<n; j++) {
	    x = gretl_matrix_get(src, t, m++);
	    if (mod == GRETL_MOD_DECREMENT) {
		x = gretl_matrix_get(targ, i, j) - x;
	    }
	    gretl_matrix_set(targ, i, j, x);
	    if (i != j) {
		gretl_matrix_set(targ, j, i, x);
	    }
	}
    }
}

/* Row @t of @src represents the vec of a certain matrix: extract the
   row and reconstitute the matrix in @targ.
*/

static int load_from_vec (gretl_matrix *targ, const gretl_matrix *src,
			   int t)
{
    int i, k = targ->rows * targ->cols;

    for (i=0; i<k; i++) {
	targ->val[i] = gretl_matrix_get(src, t, i);
    }

    return 0;
}

/* Partial implementation of Koopman's "disturbance smoother".
   See Koopman, Shephard and Doornik, section 4.4.  Needs
   more work.
*/

static int koopman_smooth (kalman *K, gretl_matrix *U)
{
    gretl_matrix_block *B;
    gretl_matrix *u, *D, *L, *R;
    gretl_matrix *r0, *r1, *r2, *N1, *N2;
    gretl_matrix *n1 = NULL;
    gretl_matrix *Mp1 = NULL;
    double x;
    int i, t, err = 0;

    if (K->p > 0) {
	Mp1 = gretl_matrix_alloc(K->p, 1);
	if (Mp1 == NULL) {
	    return E_ALLOC;
	}
    }

    B = gretl_matrix_block_new(&u,  K->n, 1,
			       &D,  K->n, K->n,
			       &L,  K->r, K->r,
			       &R,  K->T, K->r,
			       &r0, K->r, 1,
			       &r1, K->r, 1,
			       &r2, K->r, 1,
			       &N1, K->r, K->r,
			       &N2, K->r, K->r,
			       NULL);

    if (B == NULL) {
	gretl_matrix_free(Mp1);
	return E_ALLOC;
    }

    gretl_matrix_zero(r1);
    gretl_matrix_zero(N1);

    for (t=K->T-1; t>=0 && !err; t--) {
	/* load et, Vt and Kt for time t */
	load_from_row(K->e, K->E, t, GRETL_MOD_NONE);
	load_from_vech(K->Vt, K->V, K->n, t, GRETL_MOD_NONE);
	load_from_vec(K->Kt, K->K, t);

	/* u_t = V_t^{-1}e_t - K_t'r_t */
	gretl_matrix_multiply(K->Vt, K->e, u);
	if (t < K->T - 1) {
	    gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_TRANSPOSE,
				      r1, GRETL_MOD_NONE,
				      u, GRETL_MOD_DECREMENT);
	}
	/* save u_t values in E */
	load_to_row(K->E, u, t);

	/* D_t = V_t^{-1} + K_t'N_tK_t */
	gretl_matrix_copy_values(D, K->Vt);
	if (t < K->T - 1) {
	    gretl_matrix_qform(K->Kt, GRETL_MOD_TRANSPOSE,
			       N1, D, GRETL_MOD_CUMULATE);
	}

	/* L_t = F - KH' */
	gretl_matrix_copy_values(L, K->F);
	gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
				  K->H, GRETL_MOD_TRANSPOSE,
				  L, GRETL_MOD_DECREMENT);

	/* save r_t values in R */
	load_to_row(R, r1, t);

	/* r_{t-1} = HV_t^{-1}e_t + L_t'r_t */
	gretl_matrix_multiply(K->H, K->Vt, K->Tmpr1);
	gretl_matrix_multiply(K->Tmpr1, K->e, r2);
	if (t < K->T - 1) {
	    gretl_matrix_multiply_mod(L, GRETL_MOD_TRANSPOSE, 
				      r1, GRETL_MOD_NONE,
				      r2, GRETL_MOD_CUMULATE);
	}
	gretl_matrix_copy_values(r1, r2);

	/* preserve r_0 for smoothing of state */
	if (t == 0) {
	    gretl_matrix_copy_values(r0, r2);
	}

	/* N_{t-1} = HV_t^{-1}H' + L'NL */
	gretl_matrix_qform(K->H, GRETL_MOD_NONE,
			   K->Vt, N2, GRETL_MOD_NONE);
	if (t < K->T - 1) {
	    gretl_matrix_qform(L, GRETL_MOD_TRANSPOSE,
			       N1, N2, GRETL_MOD_CUMULATE);
	}
	gretl_matrix_copy_values(N1, N2);
    }

#if 0
    gretl_matrix_print(K->E, "u_t, all t");
    gretl_matrix_print(R, "r_t, all t");
#endif

    if (K->R != NULL && U != NULL) {
	/* we need an n x 1 temp matrix */
	n1 = gretl_matrix_reuse(K->Tmpnn, K->n, 1);
    }

    /* smoothed disturbances, all time steps */
    for (t=0; t<K->T; t++) {
#if 0
	if (K->p > 0) {
	    /* cross-correlated disturbances */
	    load_from_row(K->e, K->E, t, GRETL_MOD_NONE);
	    load_from_row(r1, R, t, GRETL_MOD_NONE);
	    gretl_matrix_multiply(K->C, GRETL_MOD_TRANSPOSE,
				  K->e, GRETL_MOD_NONE,
				  Mp1, GRETL_MOD_NONE);
	    gretl_matrix_multiply(K->B, GRETL_MOD_TRANSPOSE,
				  r1, GRETL_MOD_NONE,
				  Mp1, GRETL_MOD_CUMULATE);
	    gretl_matrix_multiply(K->B, Mp1, r2);
	    load_to_row(R, r2, t);
	} 
#else
	load_from_row(r1, R, t, GRETL_MOD_NONE);
	gretl_matrix_multiply(K->Q, r1, r2);
	load_to_row(R, r2, t);
	if (U != NULL) {
	    for (i=0; i<K->r; i++) {
		x = gretl_vector_get(r2, i);
		gretl_matrix_set(U, t, i, x);
	    }
	}
	if (n1 != NULL) {
	    load_from_row(K->e, K->E, t, GRETL_MOD_NONE);
	    gretl_matrix_multiply(K->R, K->e, n1);
	    for (i=0; i<K->n; i++) {
		x = gretl_vector_get(n1, i);
		gretl_matrix_set(U, t, K->r + i, x);
	    }
	}	
#endif
    }

    if (K->R != NULL && U != NULL) {
	gretl_matrix_reuse(K->Tmpnn, K->n, K->n);
    }

    /* write initial smoothed state into first row of S */
    if (K->Pini != NULL) {
	gretl_matrix_multiply(K->Pini, r0, K->S0);
    } else {
	construct_Pini(K);
	gretl_matrix_multiply(K->P0, r0, K->S0);
    }
    if (K->Sini != NULL) {
	gretl_matrix_add_to(K->S0, K->Sini);
    }
    load_to_row(K->S, K->S0, 0);

    /* smoothed state, remaining time steps */
    for (t=1; t<K->T; t++) {
	/* S_{t+1} = FS_t + v_t */
	load_from_row(K->S0, K->S, t-1, GRETL_MOD_NONE);
	gretl_matrix_multiply(K->F, K->S0, K->S1);
	load_from_row(K->S1, R, t-1, GRETL_MOD_CUMULATE);
	load_to_row(K->S, K->S1, t);
    }

    gretl_matrix_block_destroy(B);

    return err;
}

static int retrieve_Ft (kalman *K, int t)
{
    if (K->step == NULL || K->step->F == NULL) {
	return E_DATA;
    } else {
	return load_from_vec((gretl_matrix *) K->F, K->step->F, t);
    }
}

static int retrieve_Ht (kalman *K, int t)
{
    if (K->step == NULL || K->step->H == NULL) {
	return E_DATA;
    } else {
	return load_from_vec((gretl_matrix *) K->H, K->step->H, t);
    }
}

/* Anderson-Moore Kalman smoothing: see Iskander Karibzhanov's
   exposition at http://www.econ.umn.edu/~karib003/help/kalcvs.htm
   This is much the clearest account I have seen (AC 2009-04-14).

   This method uses S_{t|t-1} and P_{t|t-1} for all t, but we can
   overwrite these with the smoothed values as we go. We also need
   stored values for the prediction error, its MSE, and the gain at
   each time step.  Note that u_t and U_t are set to zero for 
   t = T - 1.
*/

static int anderson_moore_smooth (kalman *K)
{
    gretl_matrix_block *B;
    gretl_matrix *L = K->Tmprr;
    gretl_matrix *u, *u1, *U, *U1;
    gretl_matrix *StT, *PtT;
    int t, err = 0;

    B = gretl_matrix_block_new(&StT, K->r, 1,
			       &PtT, K->r, K->r,
			       &u,   K->r, 1,
			       &u1,  K->r, 1,
			       &U,   K->r, K->r,
			       &U1,  K->r, K->r,
			       NULL);
    if (B == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_zero(u);
    gretl_matrix_zero(U);

    for (t=K->T-1; t>=0 && !err; t--) {
	/* get F_t and/or H_t if need be */
	if (matrix_is_varying(K, K_F)) {
	    err = retrieve_Ft(K, t);
	}
	if (!err && matrix_is_varying(K, K_H)) {
	    err = retrieve_Ht(K, t);
	}
	if (err) {
	    break;
	}

	/* L_t = F_t - K_t H_t' */
	gretl_matrix_copy_values(L, K->F);
	load_from_vec(K->Kt, K->K, t);
	gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
				  K->H, GRETL_MOD_TRANSPOSE,
				  L, GRETL_MOD_DECREMENT);

	/* u_{t-1} = H_t V_t e_t + L_t' u_t */
	load_from_vech(K->Vt, K->V, K->n, t, GRETL_MOD_NONE);
	load_from_row(K->e, K->E, t, GRETL_MOD_NONE);
	gretl_matrix_multiply(K->Vt, K->e, K->Ve);
	gretl_matrix_multiply(K->H, K->Ve, u1);
	if (t == K->T - 1) {
	    gretl_matrix_multiply(K->H, K->Ve, u);
	} else {
	    gretl_matrix_multiply_mod(L, GRETL_MOD_TRANSPOSE,
				      u, GRETL_MOD_NONE,
				      u1, GRETL_MOD_CUMULATE);
	    gretl_matrix_copy_values(u, u1);
	}

	/* U_{t-1} = H_t V_t H_t' + L_t' U_t L_t */
	if (t == K->T - 1) {
	    gretl_matrix_qform(K->H, GRETL_MOD_NONE,
			       K->Vt, U, GRETL_MOD_NONE);
	} else {
	    gretl_matrix_qform(K->H, GRETL_MOD_NONE,
			       K->Vt, U1, GRETL_MOD_NONE);
	    gretl_matrix_qform(L, GRETL_MOD_TRANSPOSE,
			       U, U1, GRETL_MOD_CUMULATE);
	    gretl_matrix_copy_values(U, U1);
	}

	/* S_{t|T} = S_{t|t-1} + P_{t|t-1} u_{t-1} */
	load_from_row(StT, K->S, t, GRETL_MOD_NONE);
	load_from_vech(K->P0, K->P, K->r, t, GRETL_MOD_NONE);
	gretl_matrix_multiply_mod(K->P0, GRETL_MOD_NONE,
				  u, GRETL_MOD_NONE,
				  StT, GRETL_MOD_CUMULATE);
	load_to_row(K->S, StT, t);

	/* P_{t|T} = P_{t|t-1} - P_{t|t-1} U_{t-1} P_{t|t-1} */
	gretl_matrix_copy_values(PtT, K->P0);
	gretl_matrix_qform(K->P0, GRETL_MOD_NONE,
			   U, PtT, GRETL_MOD_DECREMENT);
	load_to_vech(K->P, PtT, K->r, t);
    }

    gretl_matrix_block_destroy(B);

    return err;
}

/* If we're doing smoothing for a system that has time-varying
   coefficients in K->F or K->H we'll record the vec of the
   coefficient matrices for each time-step on the forward pass.  Here
   we allocate the required storage.
*/

static int kalman_add_stepinfo (kalman *K)
{
    int err = 0;

    K->step = malloc(sizeof *K->step);

    if (K->step == NULL) {
	return E_ALLOC;
    }

    K->step->F = K->step->H = NULL;

    if (matrix_is_varying(K, K_F)) {
	K->step->F = gretl_matrix_alloc(K->T, K->r * K->r);
	if (K->step->F == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && matrix_is_varying(K, K_H)) {
	K->step->H = gretl_matrix_alloc(K->T, K->r * K->n);
	if (K->step->H == NULL) {
	    err = E_ALLOC;
	}
    }
    
    if (err) {
	free_stepinfo(K);
    }

    return err;
}

/**
 * kalman_arma_smooth:
 * @K: pointer to Kalman struct.
 * @err: location to receive error code.
 * 
 * Runs a filtering pass followed by a smoothing pass.
 *
 * Returns: matrix containing the smoothed estimate of the
 * dependent variable, or NULL on error.
 */

gretl_matrix *kalman_arma_smooth (kalman *K, int *err)
{
    int nr = (K->r * K->r + K->r) / 2;
    int nn = (K->n * K->n + K->n) / 2;
    gretl_vector *ys = NULL;

    /* Set up the matrices we need to store computed results from all
       time steps on the forward pass: state S_{t|t-1},(inverse) error 
       variance, gain, and MSE of state, P_{t|t-1}.
    */

    K->S = gretl_matrix_alloc(K->T, K->r);
    K->V = gretl_matrix_alloc(K->T, nn);
    K->K = gretl_matrix_alloc(K->T, K->r * K->n);
    K->P = gretl_matrix_alloc(K->T, nr);

    if (K->S == NULL || K->V == NULL || K->K == NULL || K->P == NULL) {
	*err = E_ALLOC;
    } else {
	double yst, Sti;
	int i, t, miss = 0;

	K->flags |= KALMAN_SMOOTH;
	*err = kalman_forecast(K, NULL);
	K->flags &= ~KALMAN_SMOOTH;
	K->t = 0;

	if (!*err) {
	    *err = anderson_moore_smooth(K);
	}

	if (!*err) {
	    ys = gretl_column_vector_alloc(K->T);
	    if (ys == NULL) {
		*err = E_ALLOC;
	    } else {
		for (t=0; t<K->T; t++) {
		    yst = 0.0;
		    for (i=0; i<K->r; i++) {
			Sti = gretl_matrix_get(K->S, t, i);
			yst += K->H->val[i] * Sti;
		    }
		    if (K->Ax != NULL) {
			K->t = t;
			kalman_set_Ax(K, &miss);
			for (i=0; i<K->n; i++) {
			    yst += gretl_vector_get(K->Ax, i);
			}
		    } 
		    gretl_vector_set(ys, t, yst);
		}
		K->t = 0;
	    }
	}

	gretl_matrix_replace(&K->S, NULL);
	gretl_matrix_replace(&K->V, NULL);
	gretl_matrix_replace(&K->K, NULL);
	gretl_matrix_replace(&K->P, NULL);
    }

    if (*err && ys != NULL) {
	gretl_matrix_free(ys);
	ys = NULL;
    }

    return ys;
}

/**
 * kalman_smooth:
 * @K: pointer to kalman struct.
 * @pP: pointer to matrix in which to retrieve the MSE of the
 * smoothed state (or NULL if this is not required).
 * @pU: pointer to matrix in which to retrieve the smoothed
 * disturbances (or NULL if this is not required).
 * @err: location to receive error code.
 * 
 * Runs a filtering pass followed by a backward, smoothing pass. 
 * At present the @pU argument is experimental and a bodge: it will 
 * not actually do anything unless @pP is left NULL.
 *
 * Returns: matrix containing the smoothed estimate of the
 * state, or NULL on error.
 */

gretl_matrix *kalman_smooth (kalman *K,
			     gretl_matrix **pP,
			     gretl_matrix **pU,
			     int *err)
{
    gretl_matrix *E, *S, *P = NULL;
    gretl_matrix *G, *V, *U = NULL;
    int nr, nn;

    if (pP == NULL && pU != NULL) {
	/* optional accessor for smoothed disturbances a la Koopman:
	   experimental, and avilable only if @pP is not given
	*/
	int Ucols = K->r;

	if (K->R != NULL) {
	    Ucols += K->n;
	}

	U = gretl_matrix_alloc(K->T, Ucols);
	if (U == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err) {
	return NULL;
    }

    nr = (K->r * K->r + K->r) / 2;
    nn = (K->n * K->n + K->n) / 2;

    /* Set up the matrices we need to store computed results from all
       time steps on the forward pass: prediction error, (inverse)
       error variance, gain, state S_{t|t-1} and MSE of state,
       P_{t|t-1}.
    */

    E = gretl_matrix_alloc(K->T, K->n);
    V = gretl_matrix_alloc(K->T, nn);
    G = gretl_matrix_alloc(K->T, K->r * K->n);
    S = gretl_matrix_alloc(K->T, K->r);
    P = gretl_matrix_alloc(K->T, nr);

    if (E == NULL || V == NULL || G == NULL || 
	S == NULL || P == NULL) {
	*err = E_ALLOC;
	goto bailout;
    } 

    if (matrix_is_varying(K, K_F) || matrix_is_varying(K, K_H)) {
	/* add recorder for F_t and/or H_t */
	*err = kalman_add_stepinfo(K);
	if (*err) {
	    goto bailout;
	}
    }

    /* attach all export matrices to Kalman */
    K->E = E;
    K->V = V;
    K->K = G;
    K->S = S;
    K->P = P;

#if 0
    /* and recheck dimensions */
    *err = user_kalman_recheck_matrices(K, prn);
#endif

    if (!*err) {
	/* forward pass */
	K->flags |= KALMAN_SMOOTH;
	*err = kalman_forecast(K, NULL);
	K->flags &= ~KALMAN_SMOOTH;
    }

    K->t = 0;

    if (!*err) {
	/* bodge */
	if (U != NULL) {
	    *err = koopman_smooth(K, U);
	} else {
	    *err = anderson_moore_smooth(K);
	}
    }

    /* detach matrices */
    K->E = NULL;
    K->V = NULL;
    K->K = NULL;
    K->S = NULL;
    K->P = NULL; 

    /* and trash the "stepinfo" storage */
    free_stepinfo(K);

 bailout:

    gretl_matrix_free(E);
    gretl_matrix_free(V);
    gretl_matrix_free(G);

    if (!*err && pP != NULL) {
	*pP = P;
    } else {
	gretl_matrix_free(P);
    }

    if (!*err && pU != NULL) {
	*pU = U;
    } else {
	gretl_matrix_free(U);
    }    

    if (*err) {
	gretl_matrix_free(S);
	S = NULL;
    }

    return S;
}

/*
 * user_kalman_smooth:
 * @Pname: name of matrix in which to retrieve the MSE of the
 * smoothed state (or NULL if this is not required).
 * @Uname: name of matrix in which to retrieve the smoothed
 * disturbances (or NULL if this is not required).
 * @err: location to receive error code.
 * 
 * If a user-defined Kalman filter is found, runs a filtering
 * pass followed by a backward, smoothing pass.  At present
 * the @Uname argument is experimental and a bodge: it will 
 * not actually do anything unless @Pname is left null.
 *
 * Returns: matrix containing the smoothed estimate of the
 * state, or NULL on error.
 */

gretl_matrix *user_kalman_smooth (const char *Pname, 
				  const char *Uname,
				  int *err)
{
    user_kalman *u = get_user_kalman(-1);
    gretl_matrix *S = NULL;

    if (u == NULL) {
	*err = missing_kalman_error();
    } else {
	gretl_matrix **pP = NULL, **pU = NULL;
	gretl_matrix *UP = NULL, *UU = NULL;
	gretl_matrix *P = NULL, *U = NULL;

	if (Pname != NULL && strcmp(Pname, "null")) {
	    UP = get_matrix_by_name(Pname);
	    if (UP == NULL) {
		*err = E_UNKVAR;
	    } else {
		pP = &P;
	    }
	} else if (Uname != NULL && strcmp(Uname, "null")) {
	    UU = get_matrix_by_name(Uname);
	    if (UU == NULL) {
		*err = E_UNKVAR;
	    } else {
		pU = &U;
	    }	
	}

	if (!*err) {
	    *err = user_kalman_recheck_matrices(u, NULL);
	}

	if (!*err) {
	    S = kalman_smooth(u->K, pP, pU, err);
	}

	if (!*err && P != NULL) {
	    user_matrix_replace_matrix_by_name(Pname, P);
	}

	if (!*err && U != NULL) {
	    user_matrix_replace_matrix_by_name(Uname, U);
	}	
    }

    return S;
}

/* See the account in Koopman, Shephard and Doornik, Econometrics
   Journal, 1999 (volume 2, pp. 113-166), section 4.2, regarding
   the initialization of the state under simulation.
*/

static int sim_state_0 (kalman *K, const gretl_matrix *V)
{
    gretl_matrix *Q, *v0 = NULL, *s0 = NULL;
    int err = 0;
    
    Q = gretl_matrix_copy(K->P0);

    if (Q == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_matrix_psd_root(Q);
    }

    if (!err) {
	v0 = gretl_matrix_alloc(K->r, 1);
	s0 = gretl_matrix_alloc(K->r, 1);
	if (v0 == NULL || s0 == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* FIXME cross-correlated? */
	load_from_row(v0, V, 0, GRETL_MOD_NONE);
	err = gretl_matrix_multiply(Q, v0, s0);
    }

    if (!err) {
	/* mess with K->S0 only on success */
	gretl_matrix_add_to(K->S0, s0);
    }

    gretl_matrix_free(Q);
    gretl_matrix_free(v0);
    gretl_matrix_free(s0);

    return err;
}

static int kalman_simulate (kalman *K, 
			    const gretl_matrix *V,
			    const gretl_matrix *W,
			    gretl_matrix *Y, 
			    gretl_matrix *S,
			    PRN *prn)
{
    gretl_matrix *yt, *et = NULL;
    int err = 0;

    yt = gretl_matrix_alloc(K->n, 1);
    if (yt == NULL) {
	return E_ALLOC;
    }

    if (K->p > 0) {
	et = gretl_matrix_alloc(K->p, 1);
	if (et == NULL) {
	    gretl_matrix_free(yt);
	    return E_ALLOC;
	}
    }	

    sim_state_0(K, V);

    if (K->x == NULL) {
	/* no exogenous vars */
	if (K->A != NULL) {
	    /* implicit const case: A is 1 x n and A'x is n x 1 */
	    gretl_vector_copy_values(K->Ax, K->A);
	} else {
	    gretl_matrix_zero(K->Ax);
	}
    } 

    for (K->t = 0; K->t < K->T; K->t += 1) {
	int missobs = 0;

	if (filter_is_varying(K)) {
	    err = kalman_refresh_matrices(K, prn);
	    if (err) {
		break;
	    }
	}

	/* y_t = A'*x_t + H'*S_t + w_t */
	gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
				  K->S0, GRETL_MOD_NONE,
				  yt, GRETL_MOD_NONE);
	if (K->x != NULL) {
	    kalman_set_Ax(K, &missobs);
	}
	if (K->A != NULL) {
	    gretl_matrix_add_to(yt, K->Ax);
	}
	if (K->p > 0) {
	    /* C \varepsilon_t */
	    load_from_row(et, V, K->t, GRETL_MOD_NONE);
	    gretl_matrix_multiply(K->cross->C, et, K->e);
	} else if (W != NULL) {
	    load_from_row(K->e, W, K->t, GRETL_MOD_NONE);
	}
	gretl_matrix_add_to(yt, K->e);

	/* record the observables */
	load_to_row(Y, yt, K->t);
	if (S != NULL) {
	    /* and the state, if wanted */
	    load_to_row(S, K->S0, K->t);
	}
	
	/* S_{t+1} = F*S_t + v_t */
	gretl_matrix_multiply(K->F, K->S0, K->S1);
	if (K->p > 0) {
	    /* B \varepsilon_t */
	    gretl_matrix_multiply_mod(K->cross->B, GRETL_MOD_NONE,
				      et, GRETL_MOD_NONE,
				      K->S1, GRETL_MOD_CUMULATE);
	} else {	    
	    load_from_row(K->S1, V, K->t, GRETL_MOD_CUMULATE);
	} 

	if (K->mu != NULL) {
	    gretl_matrix_add_to(K->S1, K->mu);
	}

	gretl_matrix_copy_values(K->S0, K->S1);
    }

    gretl_matrix_free(yt);
    gretl_matrix_free(et);

    return err;
}

/*
 * user_kalman_simulate:
 * @V: artificial disturbance to state.
 * @W: artificial disturbance to observation.
 * @Sname: name of matrix to retrieve simulated state (or NULL).
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 * 
 * If a user-defined Kalman filter is found, use it to 
 * construct a simulation based on the artificial disturbances
 * @V and (possibly) @W.  If the disturbances are not correlated 
 * across the two equations, then @V should contain the 
 * disturbances in the state equation and @W should contain 
 * those in the observation equation (if any).  But if the 
 * disturbances are correlated, then @V should contain the
 * "combined" disturbance vector at each time step, and
 * @W should be left NULL.
 *
 * Returns: matrix containing the simulated values of the
 * observables, or NULL on failure.
 */

gretl_matrix *user_kalman_simulate (const gretl_matrix *V, 
				    const gretl_matrix *W,
				    const char *Sname, 
				    PRN *prn, int *err)
{
    user_kalman *u = get_user_kalman(-1);
    gretl_matrix *Y = NULL, *S = NULL;
    kalman *K;
    int T, saveT;

    if (V == NULL) {
	fprintf(stderr, "ksimul: V is NULL\n");
	*err = missing_matrix_error(NULL);
    } else if (u == NULL) {
	*err = missing_kalman_error();
    } 
    
    if (*err) {
	return NULL;
    }

    K = u->K;

    if (K->p > 0) {
	/* If K->p > 0 we're in the cross-correlated case: we'll interpret
	   @V as the underlying disturbance matrix (and ignore @W).
	*/
	if (V->cols != K->p) {
	    *err = E_NONCONF;
	}
    } else {	
	/* Otherwise @V must have r columns, and @W must be given if
	   K->R is non-null.
	*/
	if (V->cols != K->r) {
	    *err = E_NONCONF;
	} else if (K->R != NULL) {
	    if (W == NULL) {
		fprintf(stderr, "ksimul: W is NULL\n");
		*err = missing_matrix_error("W");
	    } else if (W->rows != V->rows || W->cols != K->n) {
		*err = E_NONCONF;
	    }
	}	    
    }

    if (*err) {
	return NULL;
    }

    /* we let V provisionally define the sample length */
    saveT = K->T;
    T = V->rows;

    /* now, are the other needed matrices in place? */
    K->flags |= KALMAN_SIM;
    K->T = T;
    *err = user_kalman_recheck_matrices(u, prn);

    /* optional accessor for simulated state */
    if (!*err && Sname != NULL && strcmp(Sname, "null")) {
	S = get_matrix_by_name(Sname);
	if (S == NULL) {
	    *err = E_UNKVAR;
	} else if (S->rows != K->T || S->cols != K->r) {
	    *err = gretl_matrix_realloc(S, K->T, K->r);
	}
    }

    /* return matrix to hold simulated observables */
    if (!*err) {
	Y = gretl_matrix_alloc(K->T, K->n);
	if (Y == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	*err = kalman_simulate(K, V, W, Y, S, prn);
    }

    if (*err) {
	gretl_matrix_free(Y);
	Y = NULL;
    }

    /* restore state */
    K->flags &= ~KALMAN_SIM;
    K->T = saveT;

    return Y;
}

/*
 * user_kalman_get_loglik:
 * 
 * Retrieves the log-likelhood calculated via the last run of 
 * a kalman forecast, if applicable.
 * 
 * Returns: ll value, or #NADBL on failure.
 */

double user_kalman_get_loglik (void)
{
    user_kalman *u = get_user_kalman(-1);

    if (u == NULL || u->K == NULL) {
	return NADBL;
    } else {
	return u->K->loglik;
    }
}

/*
 * user_kalman_get_matrix:
 * @idx: identifier for matrix.
 * @err: location to receive error code.
 * 
 * Retrieves a matrix, specified by @idx, from the last
 * run of a kalman forecast, if applicable.
 * 
 * Returns: allocated matrix, or NULL on failure.
 */

gretl_matrix *user_kalman_get_matrix (int idx, int *err)
{
    user_kalman *u = get_user_kalman(-1);
    gretl_matrix *m = NULL;

    if (u == NULL || u->K == NULL) {
	*err = E_BADSTAT;
    } else {
	const gretl_matrix *src = NULL;

	if (idx == M_KLLT) {
	    src = u->K->LL;
	} else if (idx == M_KUHAT) {
	    src = u->K->e;
	}

	if (src == NULL) {
	    *err = E_BADSTAT;
	} else {
	    m = gretl_matrix_copy(src);
	    if (m == NULL) {
		*err = E_ALLOC;
	    }
	}
    } 

    return m;
}

/*
 * user_kalman_get_s2:
 * 
 * Retrieves the scale factor, \hat{\sigma}^2, calculated 
 * via the last run of a kalman forecast, if applicable.
 * 
 * Returns: scale value, or #NADBL on failure.
 */

double user_kalman_get_s2 (void)
{
    user_kalman *u = get_user_kalman(-1);

    if (u == NULL || u->K == NULL) {
	return NADBL;
    } else {
	return u->K->s2;
    }
}

/*
 * user_kalman_get_time_step:
 * 
 * Retrieves the time step, t, from the current run of a
 * kalman forecast, if applicable.
 * 
 * Returns: time-step value (>= 1), or 0 on failure.
 */

int user_kalman_get_time_step (void)
{
    user_kalman *u = get_user_kalman(-1);

    if (u == NULL || u->K == NULL) {
	return 0;
    } else if (!kalman_is_running(u->K)) {
	return kalman_checking(u->K) ? 1 : 0;
    } else {
	return (int) u->K->t + 1;
    }
}
