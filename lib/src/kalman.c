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
#include "usermat.h"
#include "gretl_func.h"
#include "matrix_extra.h"
#include "libset.h"
#include "gretl_scalar.h"
#include "kalman.h"

#define KDEBUG 0

typedef struct ksmoother_ ksmoother;

struct ksmoother_ {
    gretl_matrix *Stt;   /* T x r: S_{t|t} */
    gretl_matrix *Ptt;   /* T x rr: P_{t|t} */
};

/* notation (cf. James Hamilton):

   State:   S_{t+1} = F*S_t + v{t+1},      E(v_t*v_t') = Q

   Obs:     y_t = A'*x_t + H'*S_t + w_t,   E(w_t*w_t') = R

*/

struct kalman_ {
    int flags;   /* for recording any options */
    int fnlevel; /* level of function execution */

    int r;  /* rows of S = number of elements in state */
    int n;  /* columns of y = number of observables */
    int k;  /* columns of A = number of exogenous vars in obs eqn */
    int p;  /* length of combined disturbance vector */
    int T;  /* rows of y = number of observations */
    int t;  /* current time step, when filtering */

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

    /* input data matrices: the order matters for various functions */
    const gretl_matrix *F; /* r x r: state transition matrix */
    const gretl_matrix *A; /* k x n: coeffs on exogenous vars, obs eqn */
    const gretl_matrix *H; /* r x n: coeffs on state variables, obs eqn */
    const gretl_matrix *Q; /* r x r: contemp covariance matrix, state eqn */
    const gretl_matrix *R; /* n x n: contemp covariance matrix, obs eqn */
    const gretl_matrix *y; /* T x n: dependent variable vector (or matrix) */
    const gretl_matrix *x; /* T x k: independent variables matrix */
    const gretl_matrix *Sini; /* r x 1: S_{1|0} */
    const gretl_matrix *Pini; /* r x r: P_{1|0} */

    /* optional array of names of input matrices */
    char **mnames;

    /* optional array of function-call strings */
    char **matcalls;

    /* optional array of time-varying matrices */
    gretl_matrix **Mt;

    /* optional smoothing-related info */
    ksmoother *smoother;

    /* optional matrices for recording extra info */
    gretl_matrix *LL;    /* T x 1: loglikelihood, all time-steps */

    /* optional run-time export matrices */
    gretl_matrix *E;   /* T x n: forecast errors, all time-steps */
    gretl_matrix *V;   /* T x nn: MSE for observables, all time-steps */
    gretl_matrix *S;   /* T x r: state vector, all time-steps */
    gretl_matrix *P;   /* T x nr: MSE for state, all time-steps */
    gretl_matrix *K;   /* T x rn: gain matrix, all time-steps */

    /* matrices needed only for cross-correlated case */
    gretl_matrix *BB;
    gretl_matrix *CC;
    gretl_matrix *BC;
    
    /* workspace matrices */
    gretl_matrix_block *B; /* holder for the following */
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
};

#define NMATCALLS 5  /* max number of time-varying matrices */

#define arma_ll(K) (K->flags & KALMAN_ARMA_LL)

#define set_kalman_running(K) (K->flags |= KALMAN_FORWARD)
#define set_kalman_stopped(K) (K->flags &= ~KALMAN_FORWARD)
#define kalman_is_running(K)  (K->flags & KALMAN_FORWARD)
#define kalman_simulating(K)  (K->flags & KALMAN_SIM)

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
    K_y,
    K_x,
    K_S,
    K_P,
    K_MMAX /* sentinel */
};

/* the matrix in question is not an external named user-matrix,
   and is "owned" by the Kalman struct */

#define kalman_owns_matrix(K,i) (K->mnames != NULL && K->mnames[i][0] == '$')

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
    gretl_matrix_free(K->BC);

    if (K->BB != NULL) {
	if (K->Q == K->BB) {
	    K->Q = NULL;
	}
	gretl_matrix_free(K->BB);
    }

    if (K->CC != NULL) {
	if (K->R == K->CC) {
	    K->R = NULL;
	}
	gretl_matrix_free(K->CC);
    }    

    gretl_matrix_block_destroy(K->B);

    if (K->mnames != NULL) {
	const gretl_matrix **mptr[] = {
	    &K->F, &K->A, &K->H, &K->Q, &K->R,
	    &K->y, &K->x, &K->Sini, &K->Pini
	};
	int i;

	for (i=0; i<K_MMAX; i++) {
	    if (kalman_owns_matrix(K, i)) {
		gretl_matrix_free((gretl_matrix *) *mptr[i]);
	    }
	}

	free_strings_array(K->mnames, K_MMAX);
    }    

    if (K->matcalls != NULL) {
	free_strings_array(K->matcalls, NMATCALLS);
    }

    if (K->Mt != NULL) {
	gretl_matrix_array_free(K->Mt, NMATCALLS);
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
	K->B = NULL;
	K->F = K->A = K->H = NULL;
	K->Q = K->R = NULL;
	K->E = K->V = K->S = K->P = K->K = NULL;
	K->y = K->x = NULL;
	K->BB = K->CC = K->BC = NULL;
	K->smoother = NULL;
	K->mnames = NULL;
	K->matcalls = NULL;
	K->Mt = NULL;
	K->flags = flags;
	K->fnlevel = 0;
	K->t = 0;
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
    } else if (i == K_S) {
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

static int missing_matrix_error (void)
{
    gretl_errmsg_set(_("kalman: a required matrix is missing"));
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
	return missing_matrix_error();
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

    K->B = gretl_matrix_block_new(&K->PH,  K->r, K->n, /* P*H */
				  &K->FPH, K->r, K->n, /* F*P*H */
				  &K->HPH, K->n, K->n, /* H'*P*H */
				  &K->Vt,  K->n, K->n, /* (H'*P*H + R)^{-1} */
				  &K->Ve,  K->n, 1,    /* (H'*P*H + R)^{-1} * E */
				  &K->PHV, K->r, K->n, /* P*H*V */
				  &K->Ax,  K->n, 1,    /* A'*x at obs t */
				  &K->Kt,  K->r, K->n, /* gain at t */
				  &K->Tmpnn, K->n, K->n,
				  &K->Tmprr, K->r, K->r,
				  &K->Tmprr_2a, K->r, K->r,
				  &K->Tmprr_2b, K->r, K->r,
				  &K->Tmpr1, K->r, 1,
				  NULL);

    if (K->B == NULL) {
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

    /* K->p is non-zero only under cross-correlation */
    K->p = (opt & OPT_C)? gretl_matrix_rows(K->Q): 0;
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
 * observation equation (or %NULL if this is not applicable).
 * @y: T x n matrix of dependent variable(s).
 * @x: T x k matrix of exogenous variable(s).  May be %NULL if there
 * are no exogenous variables, or if there's only a constant.
 * @E: T x n matrix in which to record forecast errors (or %NULL if
 * this is not required).
 * @err: location to receive error code.
 *
 * Allocates and initializes a Kalman struct, which can subsequently
 * be used for forecasting with kalman_forecast().  The nomenclature
 * for the various required matrices is that in Hamilton's Time
 * Series Analysis (1994, chapter 13), except that "S" is used in
 * place of Hamilton's \xi for the state vector.
 *
 * Returns: pointer to allocated struct, or %NULL on failure, in
 * which case @err will receive a non-zero code.
 */

kalman *kalman_new (const gretl_matrix *S, const gretl_matrix *P,
		    const gretl_matrix *F, const gretl_matrix *A,
		    const gretl_matrix *H, const gretl_matrix *Q,
		    const gretl_matrix *R, const gretl_matrix *y,
		    const gretl_matrix *x, gretl_matrix *E,
		    int *err)
{
    kalman *K;

    *err = 0;

    if (F == NULL || H == NULL || Q == NULL || y == NULL) {
	*err = missing_matrix_error();
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

/* kalman_revise_variance: the user has actually given B in place
   of Q and C in place of R; so we have to form K->Q = BB' and
   K->R = CC', and also K->BC = BC'.  
*/

static int kalman_revise_variance (kalman *K)
{
    if (K->Q == NULL || K->R == NULL) {
	return missing_matrix_error();
    }

    if (K->BB == NULL) {
	K->BB = gretl_matrix_alloc(K->r, K->p);
    }

    if (K->CC == NULL) {
	K->CC = gretl_matrix_alloc(K->n, K->p);
    }

    if (K->BC == NULL) {
	K->BC = gretl_matrix_alloc(K->r, K->n);
    }   

    if (K->BB == NULL || K->CC == NULL || K->BC == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_multiply_mod(K->Q, GRETL_MOD_NONE,
			      K->Q, GRETL_MOD_TRANSPOSE,
			      K->BB, GRETL_MOD_NONE);

    gretl_matrix_multiply_mod(K->R, GRETL_MOD_NONE,
			      K->R, GRETL_MOD_TRANSPOSE,
			      K->CC, GRETL_MOD_NONE);

    gretl_matrix_multiply_mod(K->Q, GRETL_MOD_NONE,
			      K->R, GRETL_MOD_TRANSPOSE,
			      K->BC, GRETL_MOD_NONE);

    if (kalman_owns_matrix(K, K_Q)) {
	gretl_matrix_free((gretl_matrix *) K->Q);
    }
    K->Q = K->BB;

    if (kalman_owns_matrix(K, K_R)) {
	gretl_matrix_free((gretl_matrix *) K->R);
    }    
    K->R = K->CC;

#if 0
    gretl_matrix_print(K->BB, "BB");
    gretl_matrix_print(K->CC, "CC");
    gretl_matrix_print(K->BC, "BC");
#endif

    return 0;
}

static int user_kalman_setup (kalman *K, gretlopt opt)
{
    int err = 0;

    if (K->F == NULL || K->H == NULL || K->Q == NULL || K->y == NULL) {
	return missing_matrix_error();
    }

    kalman_set_dimensions(K, opt);

    if (K->p > 0) {
	err = kalman_revise_variance(K);
    }

    if (!err) {
	err = kalman_check_dimensions(K);
    }

    if (err) {
	fprintf(stderr, "failed in user_kalman_setup\n");
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

static int kalman_arma_iter_1 (kalman *K)
{
    double Ve;
    int i, err = 0;

    /* write F*S into S+ */
    err += multiply_by_F(K, K->S0, K->S1, 0);

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

    return err;
}

static void smoother_record_state_as_is (kalman *K)
{
    load_to_row(K->smoother->Stt, K->S0, K->t);
}

/* record S_{t|t} for smoothing: FIXME cross-correlated */

static void smoother_record_state (kalman *K)
{
    gretl_matrix_multiply(K->Vt, K->e, K->Ve);
    gretl_matrix_multiply(K->PH, K->Ve, K->Tmpr1);
    gretl_matrix_add_to(K->Tmpr1, K->S0);
    load_to_row(K->smoother->Stt, K->Tmpr1, K->t);
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

    if (missobs) {
	/* the observable is not in fact observed at t */
	*llt = 0.0;
	if (K->smoother != NULL) {
	    smoother_record_state_as_is(K);
	}
	if (K->K != NULL) {
	    set_row(K->K, K->t, 0.0);
	}
	return 0;
    }   

    /* form e = y - A'x - H'S */
    err += gretl_matrix_subtract_from(K->e, K->Ax);
    gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
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
    if (K->BC != NULL) {
	/* cross-correlated case */
	err += gretl_matrix_add_to(K->FPH, K->BC);
    }
    err += gretl_matrix_multiply(K->FPH, K->Vt, K->Kt);

    /* form K_t * e_t and add to S+ */
    err += gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
				     K->e,  GRETL_MOD_NONE,
				     K->S1, GRETL_MOD_CUMULATE);

    if (!err && K->smoother != NULL) {
	smoother_record_state(K);
    }

    if (!err && K->K != NULL) {
	/* record the gain */
	load_to_vec(K->K, K->Kt, K->t);
    }

    return err;
}

#define P0MIN 1.0e-15

static int record_Ptt (kalman *K)
{
    double x;
    int i;

    for (i=0; i<K->r; i++) {
	x = gretl_matrix_get(K->P0, i, i);
	if (x < 0) {
	    /* problem! */
	    if (fabs(x) < P0MIN) {
		gretl_matrix_set(K->P0, i, i, 0.0);
	    } else {
		fprintf(stderr, "P(%d,%d) = %g\n", i, i, x);
		return E_DATA;
	    }
	}
    }

    /* record vech of P_{t|t} for smoothing */
    load_to_vech(K->smoother->Ptt, K->P0, K->r, K->t);

    return 0;
}

/* Hamilton (1994) equation [13.2.22] page 380, in simplified notation:

      P+ = F[P - PH(H'PH + R)^{-1}H'P]F' + Q 

   Or use the alternative calculation

      P+ = FPF' - KVK + Q  where K = gain and V = error variance
*/

static int kalman_iter_2 (kalman *K, int missobs)
{
    int err = 0;

    if (K->p > 0) {
	/* record P for applying alternative update using gain */
	gretl_matrix_copy_values(K->Tmprr, K->P0);
    }

    if (!missobs && (K->p == 0 || K->smoother != NULL)) {
	/* form P_{t|t} = P - PH(H'PH + R)^{-1}H'P */
	err = gretl_matrix_qform(K->PH, GRETL_MOD_NONE, K->Vt,
				 K->P0, GRETL_MOD_DECREMENT);
    }

    if (!err && K->smoother != NULL) {
	/* FIXME cross-correlated */
	record_Ptt(K);
    }

    if (K->p > 0) {
	/* FPF' */
	gretl_matrix_qform(K->F, GRETL_MOD_NONE, K->Tmprr,
			   K->P1, GRETL_MOD_NONE);
	if (!missobs) {
	    /* - K(H'PH + R)K' */
	    gretl_matrix_qform(K->Kt, GRETL_MOD_NONE, K->HPH, 
			       K->P1, GRETL_MOD_DECREMENT);
	}
    } else {
	/* pre-multiply by F, post-multiply by F' */
	if (K->nonshift == K->r) {
	    err += gretl_matrix_qform(K->F, GRETL_MOD_NONE, K->P0,
				      K->P1, GRETL_MOD_NONE);
	} else {
	    err += multiply_by_F(K, K->P0, K->Tmprr, 0);
	    err += multiply_by_F(K, K->Tmprr, K->P1, 1);
	    gretl_matrix_xtr_symmetric(K->P1);
	}
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

static void kalman_set_Ax (kalman *K)
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
	    *missobs += 1;
	}
    }
}

/* the user gave a function call that should be used for updating
   a given Kalman input matrix
*/

static gretl_matrix *kalman_update_matrix (kalman *K, int i, int *err)
{
    const char *fncall = K->matcalls[i];
    gretl_matrix *m = NULL;

    if (K->mnames[i][0] != '\0') {
	/* named matrix updated via void function call */
	*err = generate(fncall, NULL, NULL, OPT_U, NULL);
	if (!*err) {
	    m = get_matrix_by_name(K->mnames[i]);
	    if (m == NULL) {
		*err = E_DATA;
	    }
	}
    } else {
	/* anonymous matrix generated by function call */
	m = generate_matrix(fncall, NULL, NULL, err);
    }

    return m;
}

/* in case we have any time-varying coefficient matrices, refresh
   these for the current time step
*/

static int kalman_refresh_matrices (kalman *K)
{
    const gretl_matrix **cptr[] = {
	&K->F, &K->A, &K->H, &K->Q, &K->R
    };    
    int i, err = 0;

    for (i=0; i<NMATCALLS && !err; i++) {
	if (K->matcalls[i] != NULL) {
	    *cptr[i] = kalman_update_matrix(K, i, &err);
	    if (K->Mt[i] != NULL) {
		gretl_matrix_free(K->Mt[i]);
		K->Mt[i] = (gretl_matrix *) *cptr[i];
	    }
	    if (!err) {
		err = check_matrix_dims(K, *cptr[i], i);
	    }	    
	    if (err) {
		fprintf(stderr, "kalman_refresh_matrices: err = %d at t = %d\n", 
			err, K->t);
	    } 
	}
    }

    if (K->p > 0) {
	/* cross-correlated case */
	if (K->BB != NULL) {
	    err = gretl_matrix_multiply_mod(K->Q, GRETL_MOD_NONE,
					    K->Q, GRETL_MOD_TRANSPOSE,
					    K->BB, GRETL_MOD_NONE);
	} 
	if (K->CC != NULL && !err) {
	    err = gretl_matrix_multiply_mod(K->R, GRETL_MOD_NONE,
					    K->R, GRETL_MOD_TRANSPOSE,
					    K->CC, GRETL_MOD_NONE);
	} 
	if (K->BC != NULL && !err) {
	    err = gretl_matrix_multiply_mod(K->Q, GRETL_MOD_NONE,
					    K->R, GRETL_MOD_TRANSPOSE,
					    K->BC, GRETL_MOD_NONE);
	} 
    }

    return err;
}

/**
 * kalman_forecast:
 * @K: pointer to Kalman struct: see kalman_new().
 *
 * Generates a series of one-step ahead forecasts for y, based on
 * information entered initially using kalman_new(), and possibly
 * modified using kalman_set_initial_state_vector() and/or
 * kalman_set_initial_MSE_matrix().  The log-likelihood is
 * calculated for the sequence of forecast errors on the assumption
 * of normality: this can be accessed using kalman_get_loglik().
 * If @E is non-%NULL, the forecast errors are recorded in this
 * matrix.
 *
 * Returns: 0 on success, non-zero on error.
 */

int kalman_forecast (kalman *K)
{
    double ldet;
    int update_P = 1;
    int Tmiss = 0;
    int i, err = 0;

#if KDEBUG
    fprintf(stderr, "kalman_forecast: T = %d\n", K->T);
#endif  

    if (K->nonshift < 0) {
	K->nonshift = count_nonshifts(K->F);
    }

    K->SSRw = K->sumldet = K->loglik = 0.0;
    K->s2 = NADBL;

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

	if (K->matcalls != NULL) {
	    /* we have time-varying coefficients */
	    err = kalman_refresh_matrices(K);
	    if (err) {
		K->loglik = NADBL;
		break;
	    }
	}

	/* read slice from y */
	kalman_initialize_error(K, &missobs);

	if (missobs) {
	    Tmiss++;
	}

	if (K->x != NULL && missobs < K->n) {
	    /* read from x if applicable */
	    kalman_set_Ax(K);
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
	    /* record variance for est. of observables */
	    if (missobs) {
		set_row(K->V, K->t, 1.0/0.0);
	    } else {
		load_to_vech(K->V, K->HPH, K->n, K->t);
	    }
	}

	/* first stage of dual iteration */
	if (arma_ll(K)) {
	    err = kalman_arma_iter_1(K);
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
	    load_to_row(K->E, K->e, K->t);
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
	    if (arma_ll(K) && update_P && K->t > 20) {
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

    if (arma_ll(K)) {
	double ll1 = 1.0 + LN_2_PI + log(K->SSRw / K->T);

	if (K->flags & KALMAN_AVG_LL) {
	    K->loglik = -0.5 * (ll1 + K->sumldet / K->T);
	} else {
	    K->loglik = -0.5 * (K->T * ll1 + K->sumldet);
	}
    } else {
	/* For K->s2 see Koopman, Shephard and Doornik, "Statistical
	   algorithms for models in state space using SsfPack 2.2",
	   Econometrics Journal, 1999, vol. 2, pp. 113-166; also
	   available at http://www.ssfpack.com/ .  For the role of
	   'd', see in addition Koopman's 1997 JASA article.
	*/
	int nT = K->n * (K->T - Tmiss);
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

    if (K->matcalls != NULL) {
	gretl_matrix_array_free(K->Mt, NMATCALLS);
	K->Mt = NULL;
    }

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
	return K->SSRw / K->T;
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
    { K_S, "inistate" },
    { K_P, "inivar" }
};

static int add_matrix_fncall (kalman *K, const char *s, int i,
			      char *mname, gretl_matrix **pm)
{
    char *p, *q, *fncall;
    int n, err = 0;

    if (K->matcalls == NULL) {
	K->matcalls = strings_array_new(NMATCALLS);
	if (K->matcalls == NULL) {
	    return E_ALLOC;
	}
	K->Mt = gretl_matrix_array_new(NMATCALLS);
	if (K->Mt == NULL) {
	    return E_ALLOC;
	}
    }

    s += strspn(s, " ");

    p = strchr(s, ';');
    q = strchr(s, '(');

    /* We may have a straight function call, or the name of a matrix
       followed by a function call, the two elements separated by
       a semicolon: in that case the separating ';' must come before
       the opening parenthesis of the function call.
    */

    if (p != NULL && q - p > 0) {
	/* we have two elements */
	n = strcspn(s, " ;");
	if (n > VNAMELEN - 1) {
	    n = VNAMELEN - 1;
	}
	strncat(mname, s, n);
	s = p + 1;
	s += strspn(s, " ");
    } 

    fncall = gretl_strdup(s);
    if (fncall == NULL) {
	return E_ALLOC;
    }

    tailstrip(fncall); 

    if (*mname == '\0') {
	/* straight function call: try it out */
	*pm = generate_matrix(fncall, NULL, NULL, &err);
	if (err) {
	    fprintf(stderr, "add_matrix_fncall: err = %d from '%s'\n", err, fncall);
	}
    } else {
	/* got a named matrix */
	*pm = get_matrix_by_name(mname);
	if (*pm == NULL) {
	    err = E_UNKVAR;
	}
    }

    if (err) {
	free(fncall);
    } else {
	K->matcalls[i] = fncall;
    }

    return err;
}

/* We didn't find a matrix of the name @s: now we try for a series
   or a named list, and if found, make a matrix out of it.  This
   matrix will be "hard-wired" -- we don't do any further look-up
   of values -- and will be owned by the Kalman struct.
*/

static gretl_matrix *k_matrix_from_dataset (char *s,
					    const double **Z, 
					    const DATAINFO *pdinfo,
					    int *err)
{
    int v = current_series_index(pdinfo, s);
    gretl_matrix *m = NULL;
	
    if (v >= 0) {
	int L[2] = {1, v};

	m = gretl_matrix_data_subset(L, Z, pdinfo->t1, pdinfo->t2,
				     M_MISSING_OK, err);
    } else {
	int *list = get_list_by_name(s);

	if (list != NULL) {
	    m = gretl_matrix_data_subset(list, Z, pdinfo->t1, pdinfo->t2,
					 M_MISSING_OK, err);
	} 
    }

    if (m != NULL) {
	strcpy(s, GENERIC_MATNAME);
    }	

    return m;
}

/* We didn't find a matrix of the name @s: now we try for a scalar
   (either a named scalar variable or a numeric constant), and if
   found, make a matrix out of it.  If we get a named scalar, we'll
   record its name so we're able to update its value later if need be.
*/

static gretl_matrix *k_matrix_from_scalar (char *s, int *err)
{
    gretl_matrix *m = NULL;
    double x;

    x = gretl_double_from_string(s, err);

    if (!*err && na(x)) {
	*err = E_MISSDATA;
    }

    if (!*err) {
	m = gretl_matrix_from_scalar(x);
	if (m == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
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
		     const double **Z, const DATAINFO *pdinfo)
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
    } else if (sscanf(s, "%15s", mname) == 1) {
	m = get_matrix_by_name(mname);
    } else {
	err = E_PARSE;
    }

    if (!err && m == NULL) {
	/* didn't find a matrix */
	if (i == K_y || i == K_x) {
	    /* if osby or obsx, try series / list */
	    m = k_matrix_from_dataset(mname, Z, pdinfo, &err);
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
	} else if (i == K_y) {
	    K->y = m;
	} else if (i == K_x) {
	    K->x = m;
	} else if (i == K_S) {
	    K->Sini = m;
	} else if (i == K_P) {
	    K->Pini = m;
	} 

	if (*mname != '\0') {
	    /* record name of matrix */
	    strcpy(K->mnames[i], mname);
	} else {
	    /* take ownership of generated matrix, in the case
	       where we got a matrix-returning function call 
	    */
	    K->Mt[i] = m;
	} 
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
	return missing_matrix_error();
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
	double x = gretl_scalar_get_value(name + 1);

	if (na(x)) {
	    /* or should we just ignore x in this case? */
	    return E_MISSDATA;
	} else {
	    m->val[0] = x;
	}
    }

    return err;
}

#define use_fncall(K,i) (K->matcalls != NULL && K->matcalls[i] != NULL)

/* When (re-)running a user-defined filter, check that no relevant
   matrices have disappeared or been resized.  In addition,
   re-initialize the state and precision.
*/

static int user_kalman_recheck_matrices (user_kalman *u)
{
    int cfd = gretl_function_depth();
    kalman *K = u->K;
    const gretl_matrix **cptr[] = {
	&K->F, &K->A, &K->H, &K->Q, &K->R,
	&K->y, &K->x, &K->Sini, &K->Pini
    };
    int i, err = 0;

    for (i=0; i<=K_P; i++) {
	if (!kalman_owns_matrix(K, i)) {
	    /* pointer may be invalid, needs refreshing */
	    *cptr[i] = NULL;
	}
    }

    for (i=0; i<=K_P && !err; i++) {
	if (i <= K_R && use_fncall(K, i)) {
	    *cptr[i] = kalman_update_matrix(K, i, &err);
	    if (K->Mt[i] != NULL) {
		gretl_matrix_free(K->Mt[i]);
		K->Mt[i] = (gretl_matrix *) *cptr[i];
	    }	
	} else if (kalman_owns_matrix(K, i)) {
	    err = update_scalar_matrix(*(gretl_matrix **) cptr[i], K->mnames[i]);
	} else {
	    *cptr[i] = kalman_retrieve_matrix(K->mnames[i], u->fnlevel, cfd);
	}
    }

    if (err) {
	return err;
    }

    if (K->F == NULL || K->H == NULL || K->Q == NULL) {
	err = missing_matrix_error();
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
   also try looking at the parent level.
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
	level--;
	for (i=0; i<n_uK; i++) {
	    if (uK[i] != NULL && uK[i]->fnlevel == level) {
		return uK[i];
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
	    if (gretl_messages_on()) {
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

/* check @s for a recognized matrix code-name */

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

/* Given a line that is part of a "kalman ... end kalman" build block,
   parse and respond appropriately.
*/

int kalman_parse_line (const char *line, const double **Z,
		       const DATAINFO *pdinfo, gretlopt opt)
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
	    err = attach_input_matrix(u->K, s, m, Z, pdinfo);
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
				     const DATAINFO *pdinfo)
{
    if (K->y != NULL) {
	int t1 = K->y->t1;
	int t2 = K->y->t2;

	if (t1 > 0 && t1 < pdinfo->n &&
	    t2 > 0 && t2 <= pdinfo->n &&
	    t2 > t1) {
	    m->t1 = t1;
	    m->t2 = t2;
	}
    }
}
#endif

static gretl_matrix *series_export_matrix (kalman *K,
					   const char *name,
					   const DATAINFO *pdinfo,
					   int i, int *v,
					   int *err)
{
    gretl_matrix *m = NULL;
    int T = sample_size(pdinfo);

    if (T == K->T) {
	if ((i == K_E || i == K_V) && K->n == 1) {
	    /* just one observable */
	    *v = current_series_index(pdinfo, name);
	} else if ((i == K_BIG_S || i == K_BIG_P) && K->r == 1) {
	    /* just one state variable */
	    *v = current_series_index(pdinfo, name);
	} else if (i == K_K && K->n == 1 && K->r == 1) {
	    /* observable and state are both 1-vectors */
	    *v = current_series_index(pdinfo, name);
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
				 const DATAINFO *pdinfo,
				 int *smat, int i)
{
    gretl_matrix *m = NULL;
    int err = 0;

    if (name != NULL && strcmp(name, "null")) {
	m = get_matrix_by_name(name);
	if (m == NULL) {
	    int v = 0;

	    m = series_export_matrix(K, name, pdinfo, i, &v, &err);
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

/* called by the kfilter() function: run a user-defined Kalman
   filter in forecasting mode */

int user_kalman_run (const char *E, const char *V, 
		     const char *S, const char *P,
		     const char *G, double **Z,
		     const DATAINFO *pdinfo,
		     int *err)
{
    user_kalman *u = get_user_kalman(-1);
    int smat[K_K+1] = {0};
    kalman *K;
    int ret = 0;

    if (u == NULL) {
	*err = missing_kalman_error();
	return 1;
    }

    K = u->K;

    if (!*err) {
	/* forecast errors */
	*err = attach_export_matrix(K, E, pdinfo, smat, K_E);
    }

    if (!*err) {
	/* MSE for observables */
	*err = attach_export_matrix(K, V, pdinfo, smat, K_V);
    } 

    if (!*err) {
	/* estimate of state */
	*err = attach_export_matrix(K, S, pdinfo, smat, K_BIG_S);
    }

    if (!*err) {
	/* MSE of estimate of state */
	*err = attach_export_matrix(K, P, pdinfo, smat, K_BIG_P);
    } 

    if (!*err) {
	/* Kalman gain */
	*err = attach_export_matrix(K, G, pdinfo, smat, K_K);
    } 

    if (!*err && K->LL == NULL) {
	/* log-likelihood: available via accessor */
	K->LL = gretl_matrix_alloc(K->T, 1);
	if (K->LL == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	*err = user_kalman_recheck_matrices(u);
    }

    if (!*err) {
	*err = kalman_forecast(K);
    }

    if (pdinfo != NULL) {
	int t1 = pdinfo->t1;

	if (smat[K_E] > 0) 
	    transcribe_and_free(Z[smat[K_E]] + t1, K->E, *err);
	if (smat[K_V > 0]) 
	    transcribe_and_free(Z[smat[K_V]] + t1, K->V, *err);
	if (smat[K_BIG_S] > 0) 
	    transcribe_and_free(Z[smat[K_BIG_S]] + t1, K->S, *err);
	if (smat[K_BIG_P] > 0) 
	    transcribe_and_free(Z[smat[K_BIG_P]] + t1, K->P, *err);
	if (smat[K_K] > 0) 
	    transcribe_and_free(Z[smat[K_K]] + t1, K->K, *err);
    }

    /* reverted (??) */
    if (*err != 0 && *err != E_NAN) {
	ret = 1;
    }

    /* detach matrices */
    K->E = K->V = K->S = NULL;
    K->P = K->K = NULL;

    return ret;    
}

/* Copy row @t from @src into @targ; or add row @t of @src to
   @targ; or subtract row @t of @src from @targ.
*/

static void load_from_row (gretl_vector *targ, const gretl_matrix *src, 
			   int t, int mod)
{
    double x;
    int i;

    for (i=0; i<src->cols; i++) {
	x = gretl_matrix_get(src, t, i);
	if (mod == GRETL_MOD_CUMULATE) {
	    targ->val[i] += x;
	} else if (mod == GRETL_MOD_DECREMENT) {
	    targ->val[i] -= x;
	} else {
	    targ->val[i] = x;
	}
    }
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

/* Below: implementation of Kalman smoothing as described in Hamilton,
   Time Series Analysis, pp. 394-397.  Needs more checking.
 */

static int kalman_smooth (kalman *K, int nr)
{
    gretl_matrix_block *B;
    gretl_matrix *J, *Pl, *Pr, *M1, *M2;
    gretl_matrix *SS, *SP, *S_, *P_;
    gretl_matrix *Stt, *Ptt;
    double x;
    int i, t;
    int err = 0;

    B = gretl_matrix_block_new(&J, K->r, K->r,
			       &Pl, K->r, K->r,
			       &Pr, K->r, K->r,
			       &M1, K->r, 1,
			       &M2, K->r, 1,
			       &S_, 1, K->r,
			       &P_, 1, nr,
			       NULL);
    if (B == NULL) {
	return E_ALLOC;
    }

    /* the t|t values */
    Stt = K->smoother->Stt;
    Ptt = K->smoother->Ptt;

    /* these will hold the smoothed values */
    SS = K->S;
    SP = K->P;

    /* In the backward recursion below we need to have available the
       "previous" (that is, t + 1) values of the original, unsmoothed
       state estimate and its MSE.  We therefore keep copies of these
       values in the vectors S_ and P_, as we overwrite the original
       matrices with the smoothed values.
    */

    /* the final value of S_{t|t} is already S_{t|T} */
    for (i=0; i<K->r; i++) {
	x = gretl_matrix_get(K->S, K->T-1, i);
	gretl_vector_set(S_, i, x); /* backup */
	x = gretl_matrix_get(Stt, K->T-1, i);
	gretl_matrix_set(SS, K->T-1, i, x);
    }

    /* also final MSE is OK */
    for (i=0; i<nr; i++) {
	x = gretl_matrix_get(K->P, K->T-1, i);
	gretl_vector_set(P_, i, x); /* backup */
	x = gretl_matrix_get(Ptt, K->T-1, i);
	gretl_matrix_set(SP, K->T-1, i, x);
    }	

    for (t=K->T-2; t>=0 && !err; t--) {
	/* J_t = P_{t|t} F' P^{-1}_{t+1|t} */
	load_from_vech(Pl, Ptt, K->r, t, GRETL_MOD_NONE);
	load_from_vech(Pr, P_, K->r, 0, GRETL_MOD_NONE);
	err = gretl_maybe_invpd(Pr);
	if (err) {
	    /* OK, try generalized inverse */
	    load_from_vech(Pr, P_, K->r, 0, GRETL_MOD_NONE);
	    err = gretl_matrix_moore_penrose(Pr);
	} 
	if (err) {
	    fprintf(stderr, "kalman_smooth: failed to invert P_{%d|%d}\n", 
		    t + 2, t + 1);
	    break;
	}
	gretl_matrix_multiply_mod(K->F, GRETL_MOD_TRANSPOSE,
				  Pr, GRETL_MOD_NONE,
				  K->Tmprr, GRETL_MOD_NONE);
	gretl_matrix_multiply(Pl, K->Tmprr, J);

	/* "M2" = J_t * (S_{t+1|T} - S{t+1|t}) */
	load_from_row(M1, SS, t+1, GRETL_MOD_NONE);
	load_from_row(M1, S_, 0, GRETL_MOD_DECREMENT);
	gretl_matrix_multiply(J, M1, M2);

	/* S_{t|T} = S_{t|t} + M2 */
	load_from_row(M1, Stt, t, GRETL_MOD_NONE);
	gretl_matrix_add_to(M1, M2);
	/* record unsmoothed state estimate */
	load_from_row(S_, K->S, t, GRETL_MOD_NONE);
	/* and set smoothed value for t */
	load_to_row(SS, M1, t);

	/* P_{t|T} = P_{t|t} + J_t(P_{t+1|T} - P_{t+1|t})J'_t */
	load_from_vech(Pr, SP, K->r, t+1, GRETL_MOD_NONE);
	load_from_vech(Pr, P_, K->r, 0, GRETL_MOD_DECREMENT);
	gretl_matrix_qform(J, GRETL_MOD_NONE,
			   Pr, K->Tmprr, GRETL_MOD_NONE);
	gretl_matrix_add_to(Pl, K->Tmprr);
	/* record vech of unsmoothed MSE */
	load_from_row(P_, K->P, t, GRETL_MOD_NONE);
	/* and set smoothed value for t */
	load_to_vech(SP, Pl, K->r, t);
    }

    gretl_matrix_block_destroy(B);

    return err;
}

static void kalman_destroy_smoother (kalman *K)
{
    gretl_matrix_free(K->smoother->Stt);
    gretl_matrix_free(K->smoother->Ptt);
    free(K->smoother);
    K->smoother = NULL;
}

static int kalman_add_smoother (kalman *K, int nr)
{
    int err = 0;

    K->smoother = malloc(sizeof *K->smoother);
    
    if (K->smoother == NULL) {
	return E_ALLOC;
    }

    K->smoother->Stt = gretl_matrix_alloc(K->T, K->r);
    K->smoother->Ptt = gretl_matrix_alloc(K->T, nr);

    if (K->smoother->Stt == NULL || K->smoother->Ptt == NULL) {
	err = E_ALLOC;
	kalman_destroy_smoother(K);
    }

    return err;
}

/* called by the ksmooth() function: run a user-defined Kalman
   filter in forecasting mode, then recurse backwards to 
   produce the smoothed estimate of the state */

gretl_matrix *user_kalman_smooth (const char *Pname, int *err)
{
    user_kalman *u = get_user_kalman(-1);
    gretl_matrix *E, *S, *P = NULL;
    kalman *K;
    int nr, user_P = 0;

    if (u == NULL) {
	*err = missing_kalman_error();
	return NULL;
    }

    K = u->K;

    /* optional accessor for MSE of smoothed state estimate */
    if (Pname != NULL && strcmp(Pname, "null")) {
	P = get_matrix_by_name(Pname);
	if (P == NULL) {
	    *err = E_UNKVAR;
	    return NULL;
	} else {
	    user_P = 1;
	}
    }

    nr = (K->r * K->r + K->r) / 2;

    E = gretl_matrix_alloc(K->T, K->n);
    S = gretl_matrix_alloc(K->T, K->r);
    if (P == NULL) {
	P = gretl_matrix_alloc(K->T, nr);
    }

    if (E == NULL || S == NULL || P == NULL) {
	gretl_matrix_free(E);
	gretl_matrix_free(S);
	if (!user_P) {
	    gretl_matrix_free(P);
	}
	*err = E_ALLOC;
	return NULL;
    } 

    *err = kalman_add_smoother(K, nr);

    if (!*err) {
	/* attach all matrices to Kalman */
	K->E = E;
	K->S = S;
	K->P = P;
	/* and recheck dimensions */
	*err = user_kalman_recheck_matrices(u);
    }

    if (!*err) {
	*err = kalman_forecast(K);
    }

    K->t = 0;

    if (!*err) {
	*err = kalman_smooth(K, nr);
    }

    kalman_destroy_smoother(K);

    /* detach matrices */
    K->E = NULL;
    K->S = NULL;
    K->P = NULL; 

    gretl_matrix_free(E);

    if (*err) {
	gretl_matrix_free(S);
	S = NULL;
    }

    if (!user_P) {
	gretl_matrix_free(P);
    }

    return S;
}

/* See the account in Koopman, Shephard and Doornik, Econometrics
   Journal, 1999, section 4.2, regarding the initialization of the
   state under simulation.
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
			    gretl_matrix *S)
{
    gretl_matrix *yt;
    int err = 0;

    yt = gretl_matrix_alloc(K->n, 1);
    if (yt == NULL) {
	return E_ALLOC;
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
	/* y_t = A'*x_t + H'*S_t + w_t */
	gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
				  K->S0, GRETL_MOD_NONE,
				  yt, GRETL_MOD_NONE);
	if (W != NULL) {
	    load_from_row(yt, W, K->t, GRETL_MOD_CUMULATE);
	} 
	if (K->x != NULL) {
	    kalman_set_Ax(K);
	    gretl_matrix_add_to(yt, K->Ax);
	}
	/* record the observables */
	load_to_row(Y, yt, K->t);
	if (S != NULL) {
	    /* and the state, if wanted */
	    load_to_row(S, K->S0, K->t);
	}	
	/* S_{t+1} = F*S_t + v{t+1} */
	gretl_matrix_multiply(K->F, K->S0, K->S1);
	load_from_row(K->S1, V, K->t, GRETL_MOD_CUMULATE);
	gretl_matrix_copy_values(K->S0, K->S1);
    }

    gretl_matrix_free(yt);

    return err;
}

gretl_matrix *user_kalman_simulate (const gretl_matrix *V, 
				    const gretl_matrix *W,
				    const char *Sname, 
				    int *err)
{
    user_kalman *u = get_user_kalman(-1);
    gretl_matrix *Y = NULL, *S = NULL;
    kalman *K;
    int T, saveT;

    if (V == NULL) {
	fprintf(stderr, "ksimul: V is NULL\n");
	*err = missing_matrix_error();
    } else if (u == NULL) {
	*err = missing_kalman_error();
    } 
    
    if (*err) {
	return NULL;
    }

    K = u->K;

    /* the state disurbance matrix, V, must have r columns */
    if (V->cols != K->r) {
	*err = E_NONCONF;
	return NULL;
    }

    /* we let V provisionally define the sample length */
    saveT = K->T;
    T = V->rows;

    /* if R is non-null we need a disturbance for the obs equation */
    if (K->R != NULL && W == NULL) {
	fprintf(stderr, "ksimul: W is NULL\n");
	*err = missing_matrix_error();
	return NULL;
    }

    /* the matrix W must be T x n, if present */
    if (W != NULL && (W->rows != T || W->cols != K->n)) {
	*err = E_NONCONF;
	return NULL;
    }

    /* now, are the other needed matrices in place? */
    K->flags |= KALMAN_SIM;
    K->T = T;
    *err = user_kalman_recheck_matrices(u);

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
	*err = kalman_simulate(K, V, W, Y, S);
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

/**
 * user_kalman_get_loglik:
 * @K: pointer to Kalman struct.
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

/**
 * user_kalman_get_matrix:
 * @K: pointer to Kalman struct.
 * @idx: identifier for matrix.
 * @err: location to receive error code.
 * 
 * Retrieves a matrix, specified by @idx, from the last
 * run of a kalman forecast, if applicable.
 * 
 * Returns: allocated matrix, or %NULL on failure.
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

/**
 * user_kalman_get_s2:
 * @K: pointer to Kalman struct.
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

/**
 * user_kalman_get_time_step:
 * @K: pointer to Kalman struct.
 * 
 * Retrieves the time step, t, from the current run of a
 * kalman forecast, if applicable.
 * 
 * Returns: scale value, or #NADBL on failure.
 */

int user_kalman_get_time_step (void)
{
    user_kalman *u = get_user_kalman(-1);

    if (u == NULL || u->K == NULL || !kalman_is_running(u->K)) {
	return 0;
    } else {
	return (int) u->K->t + 1;
    }
}
