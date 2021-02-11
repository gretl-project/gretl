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

#define FULL_XML_HEADERS 1

#include "libgretl.h"
#include "uservar.h"
#include "gretl_func.h"
#include "gretl_xml.h"
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

    int r;   /* rows of S = number of elements in state */
    int n;   /* columns of y = number of observables */
    int k;   /* columns of A = number of exogenous vars in obs eqn */
    int p;   /* length of combined disturbance vector */
    int T;   /* rows of y = number of observations */
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
       functions, including matrix_is_varying()
    */
    gretl_matrix *F; /* r x r: state transition matrix */
    gretl_matrix *A; /* k x n: coeffs on exogenous vars, obs eqn */
    gretl_matrix *H; /* r x n: coeffs on state variables, obs eqn */
    gretl_matrix *Q; /* r x r: contemp covariance matrix, state eqn */
    gretl_matrix *R; /* n x n: contemp covariance matrix, obs eqn */
    gretl_matrix *mu; /* r x 1: constant term in state transition */
    gretl_matrix *y;  /* T x n: dependent variable vector (or matrix) */
    gretl_matrix *x;  /* T x k: independent variables matrix */
    gretl_matrix *Sini; /* r x 1: S_{1|0} */
    gretl_matrix *Pini; /* r x r: P_{1|0} */

    /* user inputs for cross-correlated disturbances */
    gretl_matrix *B; /* r x p: BB' = Q */
    gretl_matrix *C; /* n x p: CC' = R */

    /* apparatus for registering time-variation of matrices */
    char *matcall;
    char *varying;

    /* optional matrices for recording extra info */
    gretl_matrix *LL;  /* T x 1: loglikelihood, all time-steps */

    /* optional run-time export matrices */
    gretl_matrix *E;   /* T x n: forecast errors, all time-steps */
    gretl_matrix *V;   /* T x nn: MSE for observables, all time-steps */
    gretl_matrix *S;   /* T x r: state vector, all time-steps */
    gretl_matrix *P;   /* T x nr: MSE for state, all time-steps */
    gretl_matrix *K;   /* T x rn: gain matrix, all time-steps */
    gretl_matrix *U;   /* T x ??: smoothed disturbances */
    gretl_matrix *Vsd; /* Variance of smoothed disturbance */

    /* structure needed only for cross-correlated case */
    crossinfo *cross;

    /* structure needed only when smoothing in the time-varying case */
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

    gretl_bundle *b; /* the bundle of which this struct is a member */
    void *data;      /* handle for attaching additional info */
    PRN *prn;        /* verbose printer */
};

/* max number of time-varying matrices: F, A, H, Q, R, mu */
#define K_N_MATCALLS 6

#define arma_ll(K) (K->flags & KALMAN_ARMA_LL)

#define set_kalman_running(K) (K->flags |= KALMAN_FORWARD)
#define set_kalman_stopped(K) (K->flags &= ~KALMAN_FORWARD)
#define kalman_is_running(K)  (K->flags & KALMAN_FORWARD)
#define kalman_simulating(K)  (K->flags & KALMAN_SIM)
#define kalman_checking(K)    (K->flags & KALMAN_CHECK)
#define kalman_xcorr(K)       (K->flags & KALMAN_CROSS)
#define kalman_ssfsim(K)      (K->flags & KALMAN_SSFSIM)

#define filter_is_varying(K) (K->matcall != NULL)

static const char *kalman_matrix_name (int sym);
static int kalman_revise_variance (kalman *K);
static int check_for_matrix_updates (kalman *K, ufunc *uf);

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
    gretl_matrix_free(c->BB);
    gretl_matrix_free(c->CC);
    gretl_matrix_free(c->BC);

    free(c);
}

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

    gretl_matrix_block_destroy(K->Blk);

    if (K->flags & KALMAN_BUNDLE) {
	gretl_matrix **mptr[] = {
	    &K->F, &K->A, &K->H, &K->Q, &K->R,
	    &K->mu, &K->y, &K->x, &K->Sini, &K->Pini
	};
	int i;

	if (kalman_xcorr(K)) {
	    mptr[3] = &K->B;
	    mptr[4] = &K->C;
	}

	for (i=0; i<K_MMAX; i++) {
	    gretl_matrix_free(*mptr[i]);
	}

	/* we also own these "export" matrices */
	gretl_matrix_free(K->E);
	gretl_matrix_free(K->V);
	gretl_matrix_free(K->S);
	gretl_matrix_free(K->P);
	gretl_matrix_free(K->K);
	gretl_matrix_free(K->U);
	gretl_matrix_free(K->Vsd);
    }

    /* time-variation info */
    free(K->matcall);
    free(K->varying);

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
	K->B = K->C = NULL;
	K->E = K->V = K->S = K->P = K->K = NULL;
	K->y = K->x = NULL;
	K->mu = NULL;
	K->U = NULL;
	K->Vsd = NULL;
	K->matcall = NULL;
	K->varying = NULL;
	K->cross = NULL;
	K->step = NULL;
	K->flags = flags;
	K->fnlevel = 0;
	K->t = 0;
	K->prn = NULL;
	K->data = NULL;
	K->b = NULL;
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
    gretl_matrix *evals;
    double r, c, x;
    int i, err = 0;

    evals = gretl_general_matrix_eigenvals(K->F, &err);

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
    int rows = K->T, cols = 0;
    int err = 0;

    if (i == K_E) {
	cols = K->n;
    } else if (i == K_V) {
	cols = (K->n * K->n + K->n) / 2;
    } else if (i == K_BIG_S) {
	cols = K->r;
    } else if (i == K_BIG_P) {
	cols = (K->r * K->r + K->r) / 2;
    } else if (i == K_LL) {
	cols = 1;
    } else if (i == K_K) {
	cols = K->r * K->n;
    } else {
	err = E_DATA;
    }

    if (!err && (m->rows != rows || m->cols != cols)) {
	err = gretl_matrix_realloc(m, rows, cols);
	if (!err) {
	    gretl_matrix_zero(m);
	}
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

    if (K->B != NULL) {
	/* cross-correlated disturbances */
	if (K->C == NULL) {
	    err = missing_matrix_error("obsymat");
	} else if (K->B->rows != K->r || K->B->cols != K->p ||
		   K->C->rows != K->n || K->C->cols != K->p) {
	    err = E_NONCONF;
	}
    } else {
	/* Q is mandatory, should be r x r and symmetric */
	if (!err) {
	    err = check_matrix_dims(K, K->Q, K_Q);
	}

	/* R should be n x n and symmetric, if present */
	if (!err && K->R != NULL) {
	    err = check_matrix_dims(K, K->R, K_R);
	}
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
	goto bailout;
    }

    K->ifc = 0;

    /* x should have T rows to match y; and it should have either k or k - 1
       columns (the latter case indicating an implicit const) */
    if (K->x != NULL) {
	if (K->x->rows < K->T) {
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
	return missing_matrix_error("obsxmat");
    }

    /* Below we have the optional "export" matrices for shipping out 
       results. If these are present but not sized correctly we'll
       try to fix them up -- but note that they are not used in a
       simulation run.
    */

    if (kalman_simulating(K)) {
	return err;
    }

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

 bailout:

    if (err) {
	fprintf(stderr, "kalman_check_dimensions: err = %d\n", err);
    }

    return err;
}

/* Write the vech of @src into row @t of @targ */

static void load_to_vech (gretl_matrix *targ,
			  const gretl_matrix *src,
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

static void load_to_vec (gretl_matrix *targ,
			 const gretl_matrix *src,
			 int t)
{
    int i;

    for (i=0; i<targ->cols; i++) {
	gretl_matrix_set(targ, t, i, src->val[i]);
    }
}

/* Write the square root of diagonal of square matrix @src
   into row @t of @targ, starting at column offset @j
*/

static void load_to_diag (gretl_matrix *targ,
			  const gretl_matrix *src,
			  int t, int j)
{
    int i, n = gretl_vector_get_length(src);
    double x;

    for (i=0; i<n; i++) {
	x = gretl_matrix_get(src, i, i);
	if (x <= 0.0) {
	    gretl_matrix_set(targ, t, i+j, 0.0);
	} else {
	    gretl_matrix_set(targ, t, i+j, sqrt(x));
	}
    }
}

/* copy from vector @src into row @t of @targ */

static void load_to_row (gretl_matrix *targ,
			 const gretl_vector *src, 
			 int t)
{
    double x;
    int i;

    for (i=0; i<targ->cols; i++) {
	x = gretl_vector_get(src, i);
	gretl_matrix_set(targ, t, i, x);
    }
}

/* copy from vector @src into row @t of @targ,
   starting at column offset @j in @targ */

static void load_to_row_offset (gretl_matrix *targ,
				const gretl_vector *src, 
				int t, int j)
{
    int i, n = gretl_vector_get_length(src);

    for (i=0; i<n; i++) {
	gretl_matrix_set(targ, t, i+j, src->val[i]);
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

static void kalman_set_dimensions (kalman *K)
{
    K->r = gretl_matrix_rows(K->F); /* F->rows defines r */
    K->k = gretl_matrix_rows(K->A); /* A->rows defines k */
    K->n = gretl_matrix_cols(K->y); /* y->cols defines n */

    if (!kalman_simulating(K)) {
	/* y->rows defines T, except when simulating */
	K->T = gretl_matrix_rows(K->y);
    }

    K->okT = K->T;

    /* K->p is non-zero only under cross-correlation; in that case the
       matrix given as 'Q' in Kalman set-up in fact represents B (as
       in v_t = B \varepsilon_t) and it must be r x p, where p is the
       number of elements in the "combined" disturbance vector
       \varepsilon_t.
    */
    K->p = (K->B != NULL)? gretl_matrix_cols(K->B): 0;
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

kalman *kalman_new (gretl_matrix *S, gretl_matrix *P,
		    gretl_matrix *F, gretl_matrix *A,
		    gretl_matrix *H, gretl_matrix *Q,
		    gretl_matrix *R, gretl_matrix *y,
		    gretl_matrix *x, gretl_matrix *m,
		    gretl_matrix *E, int *err)
{
    kalman *K;

    *err = 0;

    if (y == NULL || F == NULL || H == NULL || Q == NULL) {
	fprintf(stderr, "kalman_new: y=%p, F=%p, H=%p, Q=%p\n",
		(void *) y, (void *) F, (void *) H, (void *) Q);
	*err = missing_matrix_error(NULL);
	return NULL;
    }

    K = kalman_new_empty(0);
    if (K == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* use pointers for input matrices, don't copy */
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

    /* output, but again use external pointer */
    K->E = E;

    kalman_set_dimensions(K);

    *err = kalman_check_dimensions(K);
    if (*err) {
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

/* supports hansl function for creating a named Kalman bundle */

kalman *kalman_new_minimal (gretl_matrix *M[], int copy[],
			    int nmat, int *err)
{
    gretl_matrix **targ[5];
    kalman *K;
    int i;

    *err = 0;

    if (M[0] == NULL || M[1] == NULL || M[2] == NULL || M[3] == NULL) {
	fprintf(stderr, "kalman_new_minimal: nmat=%d, y=%p, H=%p, F=%p, Q=%p\n",
		nmat, (void *) M[0], (void *) M[1], (void *) M[2], (void *) M[3]);
	*err = missing_matrix_error(NULL);
	return NULL;
    }

    K = kalman_new_empty(KALMAN_USER | KALMAN_BUNDLE);
    if (K == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    targ[0] = &K->y;
    targ[1] = &K->H;
    targ[2] = &K->F;

    if (nmat == 5) {
	K->flags |= KALMAN_CROSS;
	targ[3] = &K->B;
	targ[4] = &K->C;
    } else {
	targ[3] = &K->Q;
    }

    for (i=0; i<nmat; i++) {
	if (copy[i]) {
	    *targ[i] = gretl_matrix_copy(M[i]);
	} else {
	    *targ[i] = M[i];
	}
    }

#if 0
    gretl_matrix_print(K->y, "K->y");
    gretl_matrix_print(K->H, "K->H");
    gretl_matrix_print(K->F, "K->F");
#endif

    kalman_set_dimensions(K);

    if (K->p > 0) {
	*err = kalman_revise_variance(K);
    }

    if (!*err) {
	*err = kalman_check_dimensions(K);
    }

    if (!*err) {
	*err = kalman_init(K);
    }
    
    if (*err) {
	kalman_free(K);
	K = NULL;
    } else {
	gretl_matrix_zero(K->e);
    }

    return K;
}

static int matrix_is_varying (kalman *K, int i)
{
    if (K->matcall != NULL) {
	if (K->varying == NULL) {
	    check_for_matrix_updates(K, NULL);
	}
	if (K->varying != NULL) {
	    return K->varying[i];
	}
    }

    return 0;
}

enum {
    UPDATE_INIT, /* initialization of matrices */
    UPDATE_STEP  /* refreshing matrices per time-step */
};

/* After reading 'Q' = B and 'R' = C from the user, either at
   (re-)initialization or at a given time-step in the case where either
   of these matrices is time-varying, record the user input in K->B
   and K->C and form the 'real' Q and R.  But note that in the
   time-step case it may be that only one of Q, R needs to be treated
   in this way (if only one is time-varying, only one will have 
   been redefined via a function call).
*/

static int kalman_update_crossinfo (kalman *K, int mode)
{
    int err = 0;

    /* Note that B and C may be needed as such for simulation */

    if (mode == UPDATE_INIT || matrix_is_varying(K, K_Q)) {
	/* recreate BB using modified B */
	err = gretl_matrix_multiply_mod(K->B, GRETL_MOD_NONE,
					K->B, GRETL_MOD_TRANSPOSE,
					K->cross->BB, GRETL_MOD_NONE);
    }

    if (!err && (mode == UPDATE_INIT || matrix_is_varying(K, K_R))) {
	/* recreate CC using modified C */
	err = gretl_matrix_multiply_mod(K->C, GRETL_MOD_NONE,
					K->C, GRETL_MOD_TRANSPOSE,
					K->cross->CC, GRETL_MOD_NONE);
    }

    if (!err && (mode == UPDATE_INIT || matrix_is_varying(K, K_Q) ||
		 matrix_is_varying(K, K_R))) {
	/* recreate BC using modified B and/or C */
	err = gretl_matrix_multiply_mod(K->B, GRETL_MOD_NONE,
					K->C, GRETL_MOD_TRANSPOSE,
					K->cross->BC, GRETL_MOD_NONE);
    }

    /* redundant? */
    if (!err && mode == UPDATE_STEP) {
	if (matrix_is_varying(K, K_Q)) {
	    K->Q = K->cross->BB;
	}
	if (matrix_is_varying(K, K_R)) {
	    K->R = K->cross->CC;
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

    c->BB = gretl_matrix_alloc(K->r, K->r);
    c->CC = gretl_matrix_alloc(K->n, K->n);
    c->BC = gretl_matrix_alloc(K->r, K->n);

    if (c->BB == NULL || c->CC == NULL || c->BC == NULL) {
	free_crossinfo(c);
	return E_ALLOC;
    }

    K->cross = c;

    return 0;
}

/* kalman_revise_variance: the user has actually given B in place of Q
   and C in place of R; so we have to form Q = BB', R = CC', and BC'
   (cross-correlated disturbances).

   This function is called in the course of initial set-up of a
   filter, and also when the Kalman matrices are being re-checked at
   the start of filtering, smoothing or simulation.
*/

static int kalman_revise_variance (kalman *K)
{
    int err = 0;

    if (K->B == NULL || K->C == NULL) {
	return missing_matrix_error("'statevar' or 'obsvar'");
    }

#if 0
    fprintf(stderr, "kalman_revise_variance\n");
    gretl_matrix_print(K->B, "B");
    gretl_matrix_print(K->C, "C");
#endif    

    if (K->cross == NULL) {
	/* not allocated yet: this should be the case only 
	   on initial set-up */
	err = kalman_add_crossinfo(K);
    }

    if (!err) {
	err = kalman_update_crossinfo(K, UPDATE_INIT);
    }

    if (!err) {
	/* establish convenience pointers */
	K->Q = K->cross->BB;
	K->R = K->cross->CC;
#if 0
	gretl_matrix_print(K->cross->BB, "BB");
	gretl_matrix_print(K->cross->CC, "CC");
	gretl_matrix_print(K->cross->BC, "BC");
#endif
    }

    if (err) {
	fprintf(stderr, "kalman_revise_variance: err = %d\n", err);
    }    

    return err;
}

static int matrix_diff (const gretl_matrix *a,
			const gretl_matrix *b,
			double tol)
{
    int i, n = a->rows * a->cols;

    /* note: we ignore the 0,0 element here */

    for (i=1; i<n; i++) {
	if (fabs(b->val[i] - a->val[i]) > tol) {
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

    /* 2018-04-07: the following (which is now standard)
       was confined to the case of KALMAN_ETT being set,
       and I believe the formula was wrong, namely
       K->e->val[0] = Ve, which amounts dividing e->val[0]
       by its estimated variance rather than standard
       deviation. AC. */
    K->e->val[0] = sqrt(K->Vt->val[0]) * K->e->val[0];

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

/* Given a unified function to update one or more of the
   potentially time-varying matrices, try to figure out
   which matrix or matrices are actually modified by
   this function.
*/

static int check_for_matrix_updates (kalman *K, ufunc *uf)
{
    char **lines;
    int nlines = 0;

    if (K->varying != NULL) {
	free(K->varying);
	K->varying = NULL;
    }

    if (uf == NULL) {
	uf = get_user_function_by_name(K->matcall);
	if (uf == NULL) {
	    gretl_errmsg_sprintf("Couldn't find function '%s'", K->matcall);
	    return E_DATA;
	}
    }
    
    K->varying = calloc(K_N_MATCALLS, 1);

    lines = gretl_function_retrieve_code(uf, &nlines);
    
    if (lines != NULL) {
	const char *bname = fn_param_name(uf, 0);
	char test[VNAMELEN+1];
	const char *s;
	int n = strlen(bname) + 1;
	int i, j;
	
	sprintf(test, "%s.", bname);
	for (i=0; i<nlines; i++) {
	    if (!strncmp(lines[i], test, n)) {
		for (j=K_F; j<=K_m; j++) {
		    s = kalman_matrix_name(j);
		    if (!strncmp(lines[i] + n, s, strlen(s))) {
			fprintf(stderr, "matrix %s is varying\n", s);
			K->varying[j] = 1;
			break;
		    }
		}
	    }
	}
	free(lines);
    }

    return 0;
}

/* Function to update any time-varying matrices, for use
   with a kalman bundle. Bypasses the regular "genr" apparatus,
   passing the attached bundle directly to the given user
   function after is has been found by name.
*/

static int kalman_update_matrices (kalman *K, PRN *prn)
{
    ufunc *uf;
    fncall *fc;
    int err = 0;

    uf = get_user_function_by_name(K->matcall);
    
    if (uf == NULL) {
	gretl_errmsg_sprintf("Couldn't find function '%s'", K->matcall);
	return E_DATA;
    }

    if (K->varying == NULL) {
	check_for_matrix_updates(K, uf);
    }

    fc = fncall_new(uf, 0);
    err = push_anon_function_arg(fc, GRETL_TYPE_BUNDLE_REF, K->b);

    if (!err) {
	err = gretl_function_exec(fc, GRETL_TYPE_NONE, NULL, NULL,
				  NULL, prn);
    }
    
    if (err) {
	fprintf(stderr, "kalman_update_matrices: call='%s', err=%d\n", 
		K->matcall, err);
    }

    return err;
}

/* If we have any time-varying coefficient matrices, refresh these for
   the current time step. This is called on a forward filtering pass.
*/

static int kalman_refresh_matrices (kalman *K, PRN *prn)
{
    gretl_matrix **mptr[] = {
	&K->F, &K->A, &K->H, &K->Q, &K->R, &K->mu
    };  
    int cross_update = 0;
    int i, err = 0;

    if (kalman_xcorr(K)) {
	mptr[3] = &K->B;
	mptr[4] = &K->C;
    }

    if (K->matcall != NULL) {
	err = kalman_update_matrices(K, prn);
    }

    for (i=0; i<K_N_MATCALLS && !err; i++) {
	if (matrix_is_varying(K, i)) {
	    if (kalman_xcorr(K) && (i == K_Q || i == K_R)) {
		/* handle revised B and/or C */
		cross_update = 1;
	    } else {
		err = check_matrix_dims(K, *mptr[i], i);
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

/* Variant of the above for use when Koopman-smoothing */

static int ksmooth_refresh_matrices (kalman *K, PRN *prn)
{
    gretl_matrix **mptr[] = {
	&K->Q, &K->R
    };
    int idx[] = {
	K_Q, K_R
    };
    int cross_update = 0;
    int i, ii, err = 0;

    if (kalman_xcorr(K)) {
	mptr[0] = &K->B;
	mptr[1] = &K->C;
    }

    if (K->matcall != NULL) {
	err = kalman_update_matrices(K, prn);
    }

    for (i=0; i<2 && !err; i++) {
	ii = idx[i];
	if (matrix_is_varying(K, ii)) {
	    if (kalman_xcorr(K) && (ii == K_Q || ii == K_R)) {
		/* handle revised B and/or C */
		cross_update = 1;
	    } else {
		err = check_matrix_dims(K, *mptr[i], ii);
	    }
	    if (err) {
		fprintf(stderr, "ksmooth_refresh_matrices: err = %d at t = %d\n", 
			err, K->t);
	    }
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
		gretl_matrix_print(K->Vt, "V");
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
		    llt = NADBL;
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
		if (!matrix_diff(K->P1, K->P0, 1.0e-20)) {
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

static int obsy_check (kalman *K)
{
    if (K->y == NULL) {
	return missing_matrix_error("obsy");
    } else if (K->y->rows != K->T || K->y->cols != K->n) {
	fprintf(stderr, "obsy_check: K->y should be %d x %d, is %d x %d\n",
		K->T, K->n, gretl_matrix_rows(K->y),
		gretl_matrix_cols(K->y));
	return E_NONCONF;
    } else {
	return 0;
    }
}

static int kalman_bundle_recheck_matrices (kalman *K, PRN *prn)
{
    int err = 0;

    K->flags |= KALMAN_CHECK;

    if (filter_is_varying(K)) {
	err = kalman_update_matrices(K, prn);
    }

    K->flags ^= KALMAN_CHECK;

    if (!err && (K->H == NULL || K->F == NULL || K->Q == NULL)) {
	fprintf(stderr, "kalman_bundle_kalman_recheck_matrices: H=%p, F=%p, Q=%p\n",
		K->H, K->F, K->Q);
	err = missing_matrix_error(NULL);
    }

    if (err) {
	return err;
    }

    /* redundant? */
    kalman_set_dimensions(K);

    if (gretl_matrix_rows(K->F) != K->r ||
	gretl_matrix_rows(K->A) != K->k) {
	err = E_NONCONF;
    } else if (!kalman_simulating(K)) {
	err = obsy_check(K);
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

static int kalman_ensure_output_matrices (kalman *K)
{
    int err = 0;
    
    if (K->E == NULL) {
	K->E = gretl_null_matrix_new();
    }
    if (K->V == NULL) {
	K->V = gretl_null_matrix_new();
    }   
    if (K->S == NULL) {
	K->S = gretl_null_matrix_new();
    }
    if (K->P == NULL) {
	K->P = gretl_null_matrix_new();
    }
    if (K->K == NULL) {
	K->K = gretl_null_matrix_new();
    }

    if (K->E == NULL || K->V == NULL || K->S == NULL ||
	K->P == NULL || K->K == NULL) {
	err = E_ALLOC;
    }

    return err;
}

int kalman_bundle_run (gretl_bundle *b, PRN *prn, int *errp)
{
    kalman *K = gretl_bundle_get_private_data(b);
    int err;

    K->b = b; /* attach bundle pointer */
    err = kalman_ensure_output_matrices(K);

    if (!err) {
	gretl_matrix_zero(K->e);
	err = kalman_bundle_recheck_matrices(K, prn);
    }

    if (!err && K->LL == NULL) {
	K->LL = gretl_matrix_alloc(K->T, 1);
	if (K->LL == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = kalman_forecast(K, prn);
    }

    if (err != E_NAN) {
	*errp = err;
    } else {
	/* we'll flag E_NAN with a return value of 1 but 
	   won't count it as a 'true' error */
	*errp = 0;
    }

    return err;    
}

/* Copy row @t from @src into @targ; or add row @t of @src to
   @targ; or subtract row @t of @src from @targ.  We allow the
   possibility that the length of vector @targ is less than
   the number of columns in @src, but not the converse.
*/

static int load_from_row (gretl_vector *targ,
			  const gretl_matrix *src, 
			  int t, GretlMatrixMod mod)
{
    int i, n = gretl_vector_get_length(targ);
    double x;

    if (n > src->cols) {
	fprintf(stderr, "load_from_row: targ length = %d, but src "
		"has %d columns\n", n, src->cols);
	return E_NONCONF;
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

/* As load_from_row(), except that a column offset, @j, is 
   supported for the reading of a row from @src, and we 
   don't support the @mod option.
*/

static int load_from_row_offset (gretl_vector *targ,
				 const gretl_matrix *src, 
				 int t, int j)
{
    int i, n = gretl_vector_get_length(targ);
    double x;

    if (n > src->cols - j) {
	return E_NONCONF;
    }

    for (i=0; i<n; i++) {
	x = gretl_matrix_get(src, t, i + j);
	targ->val[i] = x;
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

/* For disturbance smoothing: ensure we have on hand
   matrices that are correctly sized to hold estimates
   of the variance of the disturbance(s) in the state
   and (if applicable) observation equations.
*/

static int maybe_resize_dist_mse (kalman *K, 
				  gretl_matrix **vvt,
				  gretl_matrix **vwt)
{
    int n = K->R == NULL ? 0 : K->n;
    int k, err = 0;

    /* combined results: how many columns do we need? */
    k = K->r + n;

    if (K->Vsd == NULL) {
	K->Vsd = gretl_matrix_alloc(K->T, k);
	if (K->Vsd == NULL) {
	    err = E_ALLOC;
	}
    } else if (K->Vsd->rows != K->T || K->Vsd->cols != k) {
	err = gretl_matrix_realloc(K->Vsd, K->T, k);
    }

    if (!err) {
	/* step-t square state matrix */
	*vvt = gretl_matrix_alloc(K->r, K->r);
	if (*vvt == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && n > 0) {
	/* step-t square obs matrix */
	*vwt = gretl_matrix_alloc(K->n, K->n);
	if (*vwt == NULL) {
	    err = E_ALLOC;
	}
    }    

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

/* Calculate the variance of the smoothed disturbances
   for the cross-correlated case. See Koopman, Shephard
   and Doornik (1998), page 19, var(\varepsilon_t|Y_n).
*/

static int combined_dist_variance (kalman *K,
				   gretl_matrix *D,
				   gretl_matrix *N,
				   gretl_matrix *vv,
				   gretl_matrix *vw,
				   gretl_matrix_block *BX,
				   int dkstyle)
{
    gretl_matrix *DC, *KN, *Veps, *NB, *NK;

    DC   = gretl_matrix_block_get_matrix(BX, 0);
    KN   = gretl_matrix_block_get_matrix(BX, 1);
    Veps = gretl_matrix_block_get_matrix(BX, 2);
    NB   = gretl_matrix_block_get_matrix(BX, 3);

    KN = gretl_matrix_reuse(KN, K->n, K->r);

    /* First chunk of Veps:
       Koopman's notation: G_t'(D_t*G_t - K_t'*N_t*H_t)
       In gretl:           C_t'(D_t*C_t - K_t'*N_t*B_t)
    */

    gretl_matrix_multiply(D, K->C, DC);
    gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_TRANSPOSE,
			      N, GRETL_MOD_NONE,
			      KN, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(KN, GRETL_MOD_NONE,
			      K->B, GRETL_MOD_NONE,
			      DC, GRETL_MOD_DECREMENT);
    gretl_matrix_multiply_mod(K->C, GRETL_MOD_TRANSPOSE,
			      DC, GRETL_MOD_NONE,
			      Veps, GRETL_MOD_NONE);

    NK = gretl_matrix_reuse(KN, K->r, K->n);

    /* Second chunk of Veps:
       Koopman's notation: H_t'(N_t*H_t - N_t*K_t*G_t)
       In gretl:           B_t'(N_t*B_t - N_t*K_t*C_t), 
       and add to the first chunk calculated above.
    */    

    gretl_matrix_multiply(N, K->B, NB);
    gretl_matrix_multiply_mod(N, GRETL_MOD_TRANSPOSE,
			      K->Kt, GRETL_MOD_NONE,
			      NK, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(NK, GRETL_MOD_NONE,
			      K->C, GRETL_MOD_NONE,
			      NB, GRETL_MOD_DECREMENT);
    gretl_matrix_multiply_mod(K->B, GRETL_MOD_TRANSPOSE,
			      NB, GRETL_MOD_NONE,
			      Veps, GRETL_MOD_CUMULATE);

    if (dkstyle) {
	/* Veps = I_p - Veps */
	double vii;
	int i;

	gretl_matrix_multiply_by_scalar(Veps, -1.0);

	for (i=0; i<K->p; i++) {
	    vii = gretl_matrix_get(Veps, i, i);
	    gretl_matrix_set(Veps, i, i, 1.0 + vii);
	}
    }

    /* Veps (p x p) holds the variance of \epsilon_t 
       conditional on Y_n: now form the per-equation
       disturbance variance matrices, @vv and @vw, for
       this time-step.
    */
    
    gretl_matrix_qform(K->B, GRETL_MOD_NONE, Veps,
		       vv, GRETL_MOD_NONE);
    gretl_matrix_qform(K->C, GRETL_MOD_NONE, Veps,
		       vw, GRETL_MOD_NONE);

    return 0;
}

/* See Koopman, Shephard and Doornik, section 4.4 */

static int koopman_smooth (kalman *K, int dkstyle)
{
    gretl_matrix_block *B, *BX = NULL;
    gretl_matrix *u, *D, *L, *R;
    gretl_matrix *r0, *r1, *r2, *N1, *N2;
    gretl_matrix *n1 = NULL;
    gretl_matrix *Vvt = NULL;
    gretl_matrix *Vwt = NULL;
    gretl_matrix *DC = NULL;
    gretl_matrix *KN = NULL;
    gretl_matrix *RHS = NULL;
    gretl_matrix *NB = NULL;
    gretl_matrix *Ut = NULL; 
    double x;
    int i, t, err = 0;

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
	return E_ALLOC;
    }

    if (K->b != NULL) {
	/* for variance of smoothed disturbances */
	err = maybe_resize_dist_mse(K, &Vvt, &Vwt);
    }

    if (K->b != NULL && K->p > 0) {
	BX = gretl_matrix_block_new(&DC,  K->n, K->p,
				    &KN,  K->n, K->r,
				    &RHS, K->p, K->p,
				    &NB,  K->r, K->p,
				    &Ut,  K->p, 1,
				    NULL);
	if (BX == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	gretl_matrix_block_destroy(B);
	gretl_matrix_block_destroy(BX);
	gretl_matrix_free(Vvt);
	gretl_matrix_free(Vwt);
	return err;
    }

    gretl_matrix_zero(r1);
    gretl_matrix_zero(N1);

    /* The backward recursion */

    for (t=K->T-1; t>=0 && !err; t--) {
	/* load et, Vt and Kt for time t */
	load_from_row(K->e, K->E, t, GRETL_MOD_NONE);
	load_from_vech(K->Vt, K->V, K->n, t, GRETL_MOD_NONE);
	load_from_vec(K->Kt, K->K, t);

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

	if (filter_is_varying(K)) {
	    /* K->Q and/or K->R may be time-varying */
	    K->t = t;
	    ksmooth_refresh_matrices(K, NULL);
	}

	/* u_t = V_t^{-1} e_t - K_t' r_t */
	gretl_matrix_multiply(K->Vt, K->e, u);
	if (t < K->T - 1) {
	    gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_TRANSPOSE,
				      r1, GRETL_MOD_NONE,
				      u, GRETL_MOD_DECREMENT);
	}
	
	/* save u_t values in E */
	load_to_row(K->E, u, t);

	if (K->Vsd != NULL && K->p == 0) {
	    /* variance of state disturbance */
	    if (dkstyle) {
		/* Q_t - Q_t N_t Q_t */
		gretl_matrix_copy_values(Vvt, K->Q);
		gretl_matrix_qform(K->Q, GRETL_MOD_TRANSPOSE,
				   N1, Vvt, GRETL_MOD_DECREMENT);

	    } else {
		/* Q_t N_t Q_t */
		gretl_matrix_qform(K->Q, GRETL_MOD_TRANSPOSE,
				   N1, Vvt, GRETL_MOD_NONE);
	    }
	    load_to_diag(K->Vsd, Vvt, t, 0);
	}

	/* D_t = V_t^{-1} + K_t' N_t K_t */
	gretl_matrix_copy_values(D, K->Vt);
	if (t < K->T - 1) {
	    gretl_matrix_qform(K->Kt, GRETL_MOD_TRANSPOSE,
			       N1, D, GRETL_MOD_CUMULATE);
	}

	if (K->R != NULL && K->Vsd != NULL && K->p == 0) {
	    /* variance of obs disturbance */
	    if (dkstyle) {
		/* R_t - R_t D_t R_t */
		gretl_matrix_copy_values(Vwt, K->R);
		gretl_matrix_qform(K->R, GRETL_MOD_TRANSPOSE,
				   D, Vwt, GRETL_MOD_DECREMENT);

	    } else {
		/* R_t D_t R_t */
		gretl_matrix_qform(K->R, GRETL_MOD_TRANSPOSE,
				   D, Vwt, GRETL_MOD_NONE);
	    }
	    load_to_diag(K->Vsd, Vwt, t, K->r);
	}

	if (K->Vsd != NULL && K->p > 0) {
	    /* variance of combined disturbance */
	    err = combined_dist_variance(K, D, N1, Vvt, Vwt, BX,
					 dkstyle);
	    if (err) {
		break;
	    } else {
		load_to_diag(K->Vsd, Vvt, t, 0);
		load_to_diag(K->Vsd, Vwt, t, K->r);
	    }
	}

	/* L_t = F - KH' */
	gretl_matrix_copy_values(L, K->F);
	gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
				  K->H, GRETL_MOD_TRANSPOSE,
				  L, GRETL_MOD_DECREMENT);

	/* save r_t values in R */
	load_to_row(R, r1, t);

	/* r_{t-1} = H V_t^{-1}e_t + L_t' r_t */
	gretl_matrix_multiply(K->H, K->Vt, K->Tmpr1);
	gretl_matrix_multiply(K->Tmpr1, K->e, r2);
	if (t < K->T - 1) {
	    gretl_matrix_multiply_mod(L, GRETL_MOD_TRANSPOSE, 
				      r1, GRETL_MOD_NONE,
				      r2, GRETL_MOD_CUMULATE);
	}
	/* transcribe for next step */
	gretl_matrix_copy_values(r1, r2);

	/* preserve r_0 for smoothing of state */
	if (t == 0) {
	    gretl_matrix_copy_values(r0, r2);
	}

	/* N_{t-1} = H V_t^{-1}H' + L' N L */
	gretl_matrix_qform(K->H, GRETL_MOD_NONE,
			   K->Vt, N2, GRETL_MOD_NONE);
	if (t < K->T - 1) {
	    gretl_matrix_qform(L, GRETL_MOD_TRANSPOSE,
			       N1, N2, GRETL_MOD_CUMULATE);
	}
	/* transcribe for next step */
	gretl_matrix_copy_values(N1, N2);
    }

#if 0
    gretl_matrix_print(K->E, "u_t, all t");
    gretl_matrix_print(R, "r_t, all t");
#endif

    if (K->R != NULL && K->U != NULL) {
	/* we need an n x 1 temp matrix */
	n1 = gretl_matrix_reuse(K->Tmpnn, K->n, 1);
    }

    /* Smoothed disturbances, all time steps */

    if (K->p > 0) {
	/* eps = B' r_t + C' e_t */
	for (t=0; t<K->T; t++) {
	    if (filter_is_varying(K)) {
		K->t = t;
		ksmooth_refresh_matrices(K, NULL);
	    }	    
	    load_from_row(r1, R, t, GRETL_MOD_NONE);
	    gretl_matrix_multiply_mod(K->B, GRETL_MOD_TRANSPOSE,
				      r1, GRETL_MOD_NONE,
				      Ut, GRETL_MOD_NONE);
	    load_from_row(K->e, K->E, t, GRETL_MOD_NONE);
	    gretl_matrix_multiply_mod(K->C, GRETL_MOD_TRANSPOSE,
				      K->e, GRETL_MOD_NONE,
				      Ut, GRETL_MOD_CUMULATE);
	    gretl_matrix_multiply(K->B, Ut, r2);
	    load_to_row(R, r2, t);
	    gretl_matrix_multiply(K->C, Ut, n1);
	    for (i=0; i<K->r; i++) {
		x = gretl_vector_get(r2, i);
		gretl_matrix_set(K->U, t, i, x);
	    }
	    for (i=0; i<K->n; i++) {
		x = gretl_vector_get(n1, i);
		gretl_matrix_set(K->U, t, K->r + i, x);
	    }	    
	}
    } else {
	/* the independent case */
	for (t=0; t<K->T; t++) {
	    if (filter_is_varying(K)) {
		K->t = t;
		ksmooth_refresh_matrices(K, NULL);
	    }	    
	    load_from_row(r1, R, t, GRETL_MOD_NONE);
	    gretl_matrix_multiply(K->Q, r1, r2);
	    load_to_row(R, r2, t);
	    if (K->U != NULL) {
		for (i=0; i<K->r; i++) {
		    x = gretl_vector_get(r2, i);
		    gretl_matrix_set(K->U, t, i, x);
		}
	    }
	    if (n1 != NULL) {
		load_from_row(K->e, K->E, t, GRETL_MOD_NONE);
		gretl_matrix_multiply(K->R, K->e, n1);
		for (i=0; i<K->n; i++) {
		    x = gretl_vector_get(n1, i);
		    gretl_matrix_set(K->U, t, K->r + i, x);
		}
	    }	
	}
    }

    if (K->R != NULL && K->U != NULL) {
	/* reset to full size */
	gretl_matrix_reuse(K->Tmpnn, K->n, K->n);
    }

    /* Write initial smoothed state into first row of S */
    
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

    /* Smoothed state, remaining time steps */

    for (t=1; t<K->T; t++) {
	/* S_{t+1} = FS_t + v_t (or + B*eps_t) */
	load_from_row(K->S0, K->S, t-1, GRETL_MOD_NONE);
	gretl_matrix_multiply(K->F, K->S0, K->S1);
	load_from_row(K->S1, R, t-1, GRETL_MOD_CUMULATE);
	load_to_row(K->S, K->S1, t);
    }

    gretl_matrix_block_destroy(B);
    gretl_matrix_block_destroy(BX);
    gretl_matrix_free(Vvt);
    gretl_matrix_free(Vwt);

    return err;
}

/* Anderson-Moore Kalman smoothing: see Iskander Karibzhanov's
   exposition at http://karibzhanov.com/help/kalcvs.htm
   This is much the clearest account I have seen (AC 2009-04-14,
   URL updated 2016-03-24).

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
   coefficient matrices for each time-step on the forward pass.
   Here we allocate the required storage.
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

/* optional matrix to hold smoothed disturbances a la Koopman */

static int ensure_U_matrix (kalman *K)
{
    int Ucols = K->r;
    int Urows = K->T;
    int err = 0;

    if (K->R != NULL) {
	Ucols += K->n;
    }

    if (K->U == NULL) {
	K->U = gretl_matrix_alloc(Urows, Ucols);
	if (K->U == NULL) {
	    err = E_ALLOC;
	}	
    } else if (K->U->rows != Urows || K->U->cols != Ucols) {
	err = gretl_matrix_realloc(K->U, Urows, Ucols);
    }

    return err;
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
    gretl_matrix *G, *V;
    int nr, nn;

    if (pP == NULL && pU != NULL) {
	/* optional accessor for smoothed disturbances a la Koopman:
	   experimental, and avilable only if @pP is not given
	*/
	*err = ensure_U_matrix(K);
	if (*err) {
	    return NULL;
	}	
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
	if (K->U != NULL) {
	    *err = koopman_smooth(K, 0);
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
	*pU = K->U;
    } else {
	gretl_matrix_free(K->U);
    }

    K->U = NULL;

    if (*err) {
	gretl_matrix_free(S);
	S = NULL;
    }

    return S;
}

int kalman_bundle_smooth (gretl_bundle *b, int dist, PRN *prn)
{
    kalman *K = gretl_bundle_get_private_data(b);    
    int err;

    if (K == NULL) {
	fprintf(stderr, "kalman_bundle_smooth: K is NULL\n");
	return E_DATA;
    }

    K->b = b; /* attach bundle pointer */
    
    err = kalman_ensure_output_matrices(K);

    if (!err && dist) {
	err = ensure_U_matrix(K);
    }

    if (err) {
	return err;
    }    

    if (matrix_is_varying(K, K_F) || matrix_is_varying(K, K_H)) {
	/* add recorder for F_t and/or H_t */
	err = kalman_add_stepinfo(K);
	if (err) {
	    goto bailout;
	}
    }

    if (!err) {
	err = kalman_bundle_recheck_matrices(K, prn);
    }

    if (!err) {
	/* forward pass */
	K->flags |= KALMAN_SMOOTH;
	err = kalman_forecast(K, NULL);
	K->flags &= ~KALMAN_SMOOTH;	
    }

    K->t = 0;

    if (!err) {
	if (dist > 1) {
	    err = koopman_smooth(K, 1);
	} else if (dist == 1) {
	    err = koopman_smooth(K, 0);
	} else {
	    err = anderson_moore_smooth(K);
	}
    }

 bailout:    

    /* trash the "stepinfo" storage */
    free_stepinfo(K);

    return err;
}

static gretl_matrix *extract_Q (kalman *K,
				const gretl_matrix *Sim0)
{
    gretl_matrix *Q;
    double x;
    int i, j;

    Q = gretl_matrix_alloc(K->r, K->r);
    
    if (Q != NULL) {
	for (i=0; i<K->r; i++) {
	    for (j=0; j<K->r; j++) {
		x = gretl_matrix_get(Sim0, i, j);
		gretl_matrix_set(Q, i, j, x);
	    }
	}
    }

    return Q;
}

/* See the account in Koopman, Shephard and Doornik, Econometrics
   Journal, 1999 (volume 2, pp. 113-166), section 4.2, regarding
   the initialization of the state under simulation.
*/

static int sim_state_0 (kalman *K, const gretl_matrix *U,
			const gretl_matrix *Sim0)
{
    gretl_matrix *Q, *v0 = NULL, *bv = NULL;
    int getroot = 1;
    int i, err = 0;

    if (!kalman_ssfsim(K)) {
	if (Sim0 != NULL) {
	    /* Sim0 contains the state for t = 1 */
	    err = gretl_matrix_copy_values(K->S0, Sim0);
	}
	/* error or not, we're done */
	return err;
    }

    /* now we're in the "ssfsim" case, emulating ssfpack */

    if (Sim0 != NULL) {
	/* Sim0 contains state variance factor
	   plus the state for t = 0 
	*/
	Q = extract_Q(K, Sim0);
	getroot = 0;
    } else {	
	Q = gretl_matrix_copy(K->P0);
    }

    if (Q == NULL) {
	err = E_ALLOC;
    } else if (getroot) {
	err = gretl_matrix_psd_root(Q, 0);
    }

    if (!err) {
	int vlen = K->p > 0 ? K->p : K->r;
	
	v0 = gretl_matrix_alloc(vlen, 1);
	if (v0 == NULL) {
	    err = E_ALLOC;
	}
    }

    if (K->p > 0) {
	bv = gretl_matrix_alloc(K->r, 1);
	if (bv == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && Sim0 != NULL) {
	/* set S0 from last row of Sim0 */
	for (i=0; i<K->r; i++) {
	    K->S0->val[i] = gretl_matrix_get(Sim0, K->r, i);
	}
    }

    if (!err) {
	/* handle the t = 0 disturbance */
	load_from_row(v0, U, 0, GRETL_MOD_NONE);
	if (K->p > 0) {
	    /* cross-correlated */
	    gretl_matrix_multiply(K->B, v0, bv);
	    gretl_matrix_multiply_mod(Q, GRETL_MOD_NONE, 
				      bv, GRETL_MOD_NONE,
				      K->S0, GRETL_MOD_CUMULATE);
	} else {
	    gretl_matrix_multiply_mod(Q, GRETL_MOD_NONE, 
				      v0, GRETL_MOD_NONE,
				      K->S0, GRETL_MOD_CUMULATE);
	}
    }

    gretl_matrix_free(Q);
    gretl_matrix_free(v0);
    gretl_matrix_free(bv);

    return err;
}

/* note: it's OK for @S to be NULL (if the simulated
   state is not wanted), so watch out for that!
*/

static int kalman_simulate (kalman *K, 
			    const gretl_matrix *U,
			    const gretl_matrix *Sim0,
			    gretl_matrix *Y, 
			    gretl_matrix *S,
			    PRN *prn)
{
    gretl_matrix *yt, *et = NULL;
    int obs_offset = 0;
    int obsdist = 0;
    int tmin = 0;
    int err = 0;

    yt = gretl_zero_matrix_new(K->n, 1);
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

    if (Y->cols == K->r + K->n) {
	/* combined (state, obs) in @Y */
	S = Y;
	obs_offset = K->r;
    }

    err = sim_state_0(K, U, Sim0);

    if (!err && kalman_ssfsim(K)) {
	if (S != NULL) {
	    load_to_row_offset(S, K->S0, 0, 0);
	}
	load_to_row_offset(Y, yt, 0, obs_offset);
	/* the first row of output is handled */
	tmin = 1;
    }

    if (!err && K->x == NULL) {
	/* no exogenous vars */
	if (K->A != NULL) {
	    /* implicit const case: A is 1 x n and A'x is n x 1 */
	    gretl_vector_copy_values(K->Ax, K->A);
	} else {
	    gretl_matrix_zero(K->Ax);
	}
    }

    if (K->p == 0 && K->R != NULL) {
	/* we want to read observation disturbances */
	obsdist = 1;
    }

    for (K->t = tmin; K->t < K->T && !err; K->t += 1) {
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
	    load_from_row(et, U, K->t, GRETL_MOD_NONE);
	    gretl_matrix_multiply(K->C, et, K->e);
	} else if (obsdist) {
	    load_from_row_offset(K->e, U, K->t, K->r);
	}
	gretl_matrix_add_to(yt, K->e);

	/* record the t-dated observables */
	load_to_row_offset(Y, yt, K->t, obs_offset);

	/* record the t-dated state? */
	if (S != NULL && tmin == 0) {
	    load_to_row_offset(S, K->S0, K->t, 0);
	}
	
	/* S_{t+1} = F*S_t + v_t */
	gretl_matrix_multiply(K->F, K->S0, K->S1);
	if (K->p > 0) {
	    /* B \varepsilon_t */
	    gretl_matrix_multiply_mod(K->B, GRETL_MOD_NONE,
				      et, GRETL_MOD_NONE,
				      K->S1, GRETL_MOD_CUMULATE);
	} else {	    
	    load_from_row(K->S1, U, K->t, GRETL_MOD_CUMULATE);
	} 

	if (K->mu != NULL) {
	    gretl_matrix_add_to(K->S1, K->mu);
	}

	/* record the (t+1)-dated state? */
	if (S != NULL && tmin == 1) {
	    load_to_row_offset(S, K->S1, K->t, 0);
	}	

	gretl_matrix_copy_values(K->S0, K->S1);
    }

    gretl_matrix_free(yt);
    gretl_matrix_free(et);

    return err;
}

static int check_simul_inputs (kalman *K,
			       const gretl_matrix *U,
			       const gretl_matrix *Sim0,
			       const gretl_matrix *SimX,
			       int ssfsim,
			       PRN *prn)
{
    int err = 0;

    if (U == NULL) {
	err = missing_matrix_error("U");
    } else {
	int ncols;

	if (K->p > 0) {
	    /* cross-correlated */
	    ncols = K->p;
	} else {
	    ncols = K->R == NULL ? K->r : K->r + K->n;
	}

	if (U->cols != ncols) {
	    pprintf(prn, "U should have %d columns but has %d\n",
		    ncols, U->cols);
	    err = E_NONCONF;
	}
    }

    if (!err && Sim0 != NULL) {
	int r = ssfsim ? K->r + 1 : K->r;
	int c = ssfsim ? K->r : 1;

	if (Sim0->rows != r || Sim0->cols != c) {
	    pprintf(prn, "simstart should be %d x %d, is %d x %d\n",
		    r, c, Sim0->rows, Sim0->cols);
	}
    }

    if (!err && K->x != NULL) {
	/* do we have enough "obsx" data? */
	const gretl_matrix *X = SimX != NULL ? SimX : K->x;

	if (X->rows < U->rows) {
	    pprintf(prn, "obsx should have %d rows but has %d\n",
		    U->rows, X->rows);
	    err = E_NONCONF;
	}
    }

    return err;
}

gretl_matrix *kalman_bundle_simulate (gretl_bundle *b,
				      const gretl_matrix *U, 
				      int get_state,
				      PRN *prn, int *err)
{
    kalman *K = gretl_bundle_get_private_data(b);
    const gretl_matrix *Sim0 = NULL;
    const gretl_matrix *SimX = NULL;
    gretl_matrix *Ret = NULL;
    gretl_matrix *savex = NULL;
    double ssfx;
    int ssfsim = 0;
    int saveT;

    if (K == NULL) {
	*err = E_DATA;
	return NULL;
    }

    /* try accessing auxiliary info from the bundle */
    Sim0 = gretl_bundle_get_matrix(b, "simstart", NULL);
    ssfx = gretl_bundle_get_scalar(b, "ssfsim", NULL);
    if (K->x != NULL) {
	SimX = gretl_bundle_get_matrix(b, "simx", NULL);
    }

    ssfsim = !na(ssfx) && ssfx != 0;

    *err = check_simul_inputs(K, U, Sim0, SimX, ssfsim, prn);
    if (*err) {
	return NULL;
    }

    K->b = b; /* attach bundle pointer */

    saveT = K->T;
    savex = K->x;

    /* we let U temporarily define the sample length */
    K->T = U->rows;

    /* and we allow temporary replacement of K->x */
    if (SimX != NULL) {
	K->x = (gretl_matrix *) SimX;
    }

    /* set state */
    if (ssfsim) {
	K->flags |= (KALMAN_SIM | KALMAN_SSFSIM);
    } else {
	K->flags |= KALMAN_SIM;
    }

    /* now, are the other needed matrices in place? */
    *err = kalman_bundle_recheck_matrices(K, prn);

    /* matrix to hold simulated observables, and state
       if wanted */
    if (!*err) {
	int ncols = get_state ? (K->r + K->n) : K->n;
	
	Ret = gretl_matrix_alloc(K->T, ncols);
	if (Ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	*err = kalman_simulate(K, U, Sim0, Ret, NULL, prn);
    }

    if (*err) {
	gretl_matrix_free(Ret);
	Ret = NULL;
    }

    /* restore state */
    K->flags &= ~KALMAN_SIM;
    K->flags &= ~KALMAN_SSFSIM;
    K->T = saveT;
    K->x = savex;

    return Ret;
}

static int matrix_is_diagonal (const gretl_matrix *m)
{
    double x;
    int i, j;

    for (j=0; j<m->cols; j++) {
	for (i=0; i<m->rows; i++) {
	    x = gretl_matrix_get(m, i, j);
	    if (i != j && x != 0.0) return 0;
	}
    }

    return 1;
}

static int simdata_refresh_QR (kalman *K, PRN *prn)
{
    int err = 0;

    if (matrix_is_varying(K, K_Q) || matrix_is_varying(K, K_R)) {
	err = kalman_update_matrices(K, prn);
    }

    return err;
}
 
/* Return a matrix in which the standard normal variates
   in @U are scaled according to K->Q, and K->R if present.
*/

gretl_matrix *kalman_bundle_simdata (gretl_bundle *b,
				     const gretl_matrix *U,
				     PRN *prn, int *err)
{
    kalman *K = gretl_bundle_get_private_data(b);
    gretl_matrix *E = NULL;

    if (K == NULL || U == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (K->p > 0) {
	if (U->cols != K->p) {
	    *err = E_DATA;
	    return NULL;
	} else {
	    /* nothing to be done, really */
	    E = gretl_matrix_copy(U);
	}
    } else {
	int n = K->R == NULL ? 0 : K->n;
	int t, j, rn = K->r + n;
	int T = U->rows;
	int varying = 0;
	double vjj, utj;

	if (U->cols != rn) {
	    *err = E_DATA;
	    return NULL;
	}

	E = gretl_matrix_alloc(U->rows, rn);
	if (E == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	}

	if (matrix_is_varying(K, K_Q) || matrix_is_varying(K, K_R)) {
	    varying = 1;
	}

	K->b = b;
	set_kalman_running(K);

	if (matrix_is_diagonal(K->Q) &&
	    (K->R == NULL || matrix_is_diagonal(K->R))) {
	    for (t=0; t<T && !*err; t++) {
		if (varying) {
		    K->t = t;
		    *err = simdata_refresh_QR(K, prn);
		}
		for (j=0; j<rn && !*err; j++) {
		    if (j < K->r) {
			vjj = gretl_matrix_get(K->Q, j, j);
		    } else {
			vjj = gretl_matrix_get(K->R, j-K->r, j-K->r);
		    }
		    utj = gretl_matrix_get(U, t, j);
		    gretl_matrix_set(E, t, j, sqrt(vjj) * utj);
		}
	    }		
	} else {
	    gretl_matrix *V = gretl_zero_matrix_new(rn, rn);
	    gretl_matrix *Ut = NULL;
	    gretl_matrix *Et = NULL;

	    if (V == NULL) {
		*err = E_ALLOC;
		goto bailout;
	    }

	    if (varying) {
		Ut = gretl_matrix_alloc(1, rn);
		Et = gretl_matrix_alloc(1, rn);
		if (Ut == NULL || Et == NULL) {
		    gretl_matrix_free(V);
		    *err = E_ALLOC;
		    goto bailout;
		}
	    }

	    if (varying) {
		for (t=0; t<T && !*err; t++) {
		    K->t = t;
		    *err = simdata_refresh_QR(K, prn);
		    if (!*err) {
			gretl_matrix_inscribe_matrix(V, K->Q, 0, 0,
						     GRETL_MOD_NONE);
			if (n > 0) {
			    gretl_matrix_inscribe_matrix(V, K->R, K->r, K->r,
							 GRETL_MOD_NONE);
			}
			*err = gretl_matrix_psd_root(V, 0);
			if (*err) {
			    gretl_errmsg_set("Failed to compute factor of Omega_t");
			} else {
			    load_from_row(Ut, U, t, GRETL_MOD_NONE);
			    gretl_matrix_multiply_mod(Ut, GRETL_MOD_NONE,
						      V, GRETL_MOD_TRANSPOSE,
						      Et, GRETL_MOD_NONE);
			    load_to_row(E, Et, t);
			}			
		    }
		}
		
		gretl_matrix_free(Ut);
		gretl_matrix_free(Et);
	    } else {
		gretl_matrix_inscribe_matrix(V, K->Q, 0, 0,
					     GRETL_MOD_NONE);
		if (n > 0) {
		    gretl_matrix_inscribe_matrix(V, K->R, K->r, K->r,
						 GRETL_MOD_NONE);
		}
		*err = gretl_matrix_psd_root(V, 0);
		if (*err) {
		    gretl_errmsg_set("Failed to compute factor of Omega");
		} else {
		    gretl_matrix_multiply_mod(U, GRETL_MOD_NONE,
					      V, GRETL_MOD_TRANSPOSE,
					      E, GRETL_MOD_NONE);
		}
	    }
	    
	    gretl_matrix_free(V);
	}
    }

 bailout:

    K->t = 0;
    set_kalman_stopped(K);
	    
    if (E == NULL && !*err) {
	*err = E_ALLOC;
    } else if (E != NULL && *err) {
	gretl_matrix_free(E);
	E = NULL;
    }

    return E;
}

/* below: functions to support the "wrapping" of a Kalman
   struct in a gretl bundle */

static int check_replacement_dims (const gretl_matrix *orig,
				   const gretl_matrix *repl,
				   int id)
{
    if (id == K_H || id == K_A || id == K_R || id == K_F ||
	id == K_Q || id == K_m || id == K_S || id == K_P) {
	if (repl->rows != orig->rows ||
	    repl->cols != orig->cols) {
	    gretl_errmsg_set("You cannot resize a state-space "
			     "system matrix");
	    return E_DATA;
	}
    } else if (id == K_y || id == K_x) {
	if (repl->cols != orig->cols) {
	    gretl_errmsg_set("You cannot resize a state-space "
			     "system matrix");
	    return E_DATA;
	}
    }

    return 0;
}

static int add_or_replace_k_matrix (gretl_matrix **targ,
				    gretl_matrix *src,
				    int id, int copy)
{
    int err = 0;
    
    if (*targ != src) {
	if (*targ != NULL) {
	    err = check_replacement_dims(*targ, src, id);
	    if (err) {
		return err;
	    }
	    /* destroy old Kalman-owned matrix */
	    gretl_matrix_free(*targ);
	}
	if (copy) { 
	    *targ = gretl_matrix_copy(src);
	    if (*targ == NULL) {
		err = E_ALLOC;
	    }
	} else {
	    *targ = src;
	}
    }

    return err;
}

static gretl_matrix **
get_input_matrix_target_by_id (kalman *K, int i)
{
    gretl_matrix **targ = NULL;
    
    if (i == K_F) {
	targ = &K->F;
    } else if (i == K_A) {
	targ = &K->A;
    } else if (i == K_H) {
	targ = &K->H;
    } else if (i == K_Q) {
	if (kalman_xcorr(K)) {
	    targ = &K->B;
	} else {
	    targ = &K->Q;
	}	
    } else if (i == K_R) {
	if (kalman_xcorr(K)) {
	    targ = &K->C;
	} else {
	    targ = &K->R;
	}
    } else if (i == K_m) {
	targ = &K->mu;
    } else if (i == K_y) {
	targ = &K->y;
    } else if (i == K_x) {
	targ = &K->x;
    } else if (i == K_S) {
	targ = &K->Sini;
    } else if (i == K_P) {
	targ = &K->Pini;
    }

    return targ;
}

/* try attaching a matrix to a Kalman bundle */

static int 
kalman_bundle_try_set_matrix (kalman *K, void *data,
			      GretlType vtype, int i,
			      int copy)
{
    gretl_matrix **targ;
    int err = 0;

    targ = get_input_matrix_target_by_id(K, i);

    if (targ == NULL) {
	err = E_DATA;
    } else {
	gretl_matrix *m;
	
	if (vtype == GRETL_TYPE_MATRIX) {
	    m = data;
	    err = add_or_replace_k_matrix(targ, m, i, copy);
	} else if (vtype == GRETL_TYPE_DOUBLE) {
	    m = gretl_matrix_from_scalar(*(double *) data);
	    if (m == NULL) {
		err = E_ALLOC;
	    } else {
		err = add_or_replace_k_matrix(targ, m, i, 0);
	    }
	} else {
	    err = E_TYPES;
	}
    }

    return err;
}

static int 
kalman_bundle_set_matrix (kalman *K, gretl_matrix *m, int i, int copy)
{
    gretl_matrix **targ;
    int err = 0;

    targ = get_input_matrix_target_by_id(K, i);

    if (targ == NULL) {
	err = E_DATA;
    } else {
	err = add_or_replace_k_matrix(targ, m, i, copy);
    }

    return err;
}

static gretl_matrix **kalman_output_matrix (kalman *K,
					    const char *key)
{
    gretl_matrix **pm = NULL;

    if (!strcmp(key, "prederr")) {
	pm = &K->E;
    } else if (!strcmp(key, "pevar")) {
	pm = &K->V;
    } else if (!strcmp(key, "state")) {
	pm = &K->S;
    } else if (!strcmp(key, "stvar")) {
	pm = &K->P;
    } else if (!strcmp(key, "gain")) {
	pm = &K->K;
    } else if (!strcmp(key, "llt")) {
	pm = &K->LL;
    } else if (!strcmp(key, "smdist")) {
	pm = &K->U;
    } else if (!strcmp(key, "smdisterr")) {
	pm = &K->Vsd;
    } else if (!strcmp(key, "uhat")) {
	pm = &K->e;
    }

    return pm;
}

#define K_N_OUTPUTS 9

static const char *kalman_output_matrix_names[K_N_OUTPUTS] = {
    "prederr",
    "pevar",
    "state",
    "stvar",
    "gain",
    "llt",
    "smdist",
    "smdisterr",
    "uhat"
};

static int output_matrix_slot (const char *s)
{
    int i;

    for (i=0; i<K_N_OUTPUTS; i++) {
	if (!strcmp(s, kalman_output_matrix_names[i])) {
	    return i;
	}
    }

    return -1;
}

#define K_N_SCALARS 9

enum {
    Ks_t = 0,
    Ks_DIFFUSE,
    Ks_CROSS,
    Ks_S2,
    Ks_LNL,
    Ks_r,
    Ks_n,
    Ks_T,
    Ks_p,     
};

static const char *kalman_output_scalar_names[K_N_SCALARS] = {
    "t",
    "diffuse",
    "cross",
    "s2",
    "lnl",
    "r",
    "n",
    "T",
    "p"     
};

static double *kalman_output_scalar (kalman *K,
				     const char *key)
{
    /* static storage for on-the-fly scalars */
    static double retval[K_N_SCALARS];
    int i, idx = -1;

    for (i=0; i<K_N_SCALARS; i++) {
	if (!strcmp(key, kalman_output_scalar_names[i])) {
	    idx = i;
	    break;
	}
    }

    if (idx < 0) {
	return NULL;
    }

    switch (idx) {
    case Ks_t:
	if (kalman_is_running(K)) {
	    retval[idx] = K->t + 1;
	} else {
	    retval[idx] = kalman_checking(K) ? 1 : 0;
	}
	break;
    case Ks_DIFFUSE:
	retval[idx] = (K->flags & KALMAN_DIFFUSE)? 1 : 0;
	break;
    case Ks_CROSS:
	retval[idx] = (K->flags & KALMAN_CROSS)? 1 : 0;
	break;
    case Ks_S2:
	retval[idx] = K->s2;
	break;
    case Ks_LNL:
	retval[idx] = K->loglik;
	break;
    case Ks_r:
	retval[idx] = K->r;
	break;
    case Ks_n:
	retval[idx] = K->n;
	break;
    case Ks_T:
	retval[idx] = K->T;
	break;
    case Ks_p:
	retval[idx] = K->p;
	break;
    default:
	break;
    }

    return &retval[idx];
}

static const gretl_matrix *k_input_matrix_by_id (kalman *K, int i)
{
    const gretl_matrix *m = NULL;
    
    if (i == K_F) {
	m = K->F;
    } else if (i == K_A) {
	m = K->A;
    } else if (i == K_H) {
	m = K->H;
    } else if (i == K_Q) {
	if (kalman_xcorr(K)) {
	    m = K->B;
	} else {
	    m = K->Q;
	}	
    } else if (i == K_R) {
	if (kalman_xcorr(K)) {
	    m = K->C;
	} else {
	    m = K->R;
	}
    } else if (i == K_m) {
	m = K->mu;
    } else if (i == K_y) {
	m = K->y;
    } else if (i == K_x) {
	m = K->x;
    } else if (i == K_S) {
	m = K->Sini;
    } else if (i == K_P) {
	m = K->Pini;
    }

    return m;
}

static int input_matrix_slot (const char *s)
{
    int i;

    for (i=0; i<K_MMAX; i++) {
	if (!strcmp(s, K_input_mats[i].name)) {
	    return K_input_mats[i].sym;
	}
    }

    return -1;
};

static GretlType kalman_extra_type (const char *key)
{
    if (!strcmp(key, "ssfsim")) {
	return GRETL_TYPE_DOUBLE;
    } else if (!strcmp(key, "simstart") ||
	       !strcmp(key, "simx")) {
	return GRETL_TYPE_MATRIX;
    } else {
	return GRETL_TYPE_NONE;
    }
}

int maybe_set_kalman_element (void *kptr,
			      const char *key,
			      void *vptr,
			      GretlType vtype,
			      int copy,
			      int *err)
{
    GretlType targ;
    kalman *K = kptr;
    int fncall = 0;
    int i, id = -1;
    int Kflag = 0;
    int done = 0;

    if (K == NULL) {
	*err = E_DATA;
	return 0;
    }

    /* check for optional "extra" kalman items that
       live outside of the kalman struct itself
    */
    targ = kalman_extra_type(key);
    if (targ != GRETL_TYPE_NONE) {
	if (vtype != targ) {
	    *err = E_TYPES;
	}
	return 0;
    }

    if (!strcmp(key, "diffuse")) {
	Kflag = KALMAN_DIFFUSE;
    }

    if (Kflag) {
	if (vtype == GRETL_TYPE_DOUBLE) {
	    double x = *(double *) vptr;

	    if (na(x)) {
		*err = E_DATA;
		return 0;
	    } else if (x == 0) {
		K->flags &= ~Kflag;
	    } else {
		K->flags |= Kflag;
	    }
	    return 1;
	} else {
	    *err = E_TYPES;
	    return 0;
	}
    }

    if (!strcmp(key, "timevar_call")) {
	/* try for a function call specifier (string) */
	if (vtype == GRETL_TYPE_STRING) {
	    fncall = 1;
	} else {
	    *err = E_TYPES;
	}
    } else {
	/* try for a matrix specifier */
	i = input_matrix_slot(key);
	if (i >= 0) {
	    if (vtype == GRETL_TYPE_MATRIX ||
		vtype == GRETL_TYPE_DOUBLE) {
		id = i;
	    } else {
		*err = E_TYPES;
	    }
	}
    }

    if (*err) {
	return 0;
    } else if (fncall) {
	if (copy) {
	    K->matcall = gretl_strdup((char *) vptr);
	} else {
	    K->matcall = (char *) vptr;
	}
	/* re-evaluate what's actually varying */
	*err = check_for_matrix_updates(K, NULL);
	if (!*err) {
	    done = 1;
	}
    } else if (id < 0) {
	if (kalman_output_matrix(K, key) != NULL ||
	    kalman_output_scalar(K, key) != NULL) {
	    *err = E_DATA;
	    gretl_errmsg_sprintf("The member %s is read-only", key);
	}
    } else {
	*err = kalman_bundle_try_set_matrix(K, vptr, vtype, id, copy);
	done = (*err == 0);
    }

    return done;
}

int maybe_delete_kalman_element (void *kptr,
				 const char *key,
				 int *err)
{
    kalman *K = kptr;
    gretl_matrix **pm;
    int done = 0;

    if (K == NULL) {
	return 0;
    }   

    if (kalman_output_scalar(K, key) != NULL ||
	input_matrix_slot(key) >= 0 || !strcmp(key, "uhat")) {
	/* note: the matrix under the key "uhat" is part of
	   the internal kalman apparatus */
	gretl_errmsg_sprintf("%s: cannot be deleted", key);
	*err = E_DATA;
    } else if ((pm = kalman_output_matrix(K, key)) != NULL) {
	/* OK to delete a user-output matrix */
	gretl_matrix_free(*pm);
	*pm = NULL;
    } else if (!strcmp(key, "timevar_call")) {
	/* OK to delete time-variation call */
	if (K->matcall != NULL) {
	    free(K->matcall);
	    K->matcall = NULL;
	    free(K->varying);
	    K->varying = NULL;
	    done = 1;
	} else {
	    *err = E_DATA;
	}
    }

    return done;
}

void *maybe_retrieve_kalman_element (void *kptr,
				     const char *key,
				     GretlType *type,
				     int *reserved,
				     int *err)
{
    kalman *K = kptr;
    void *ret = NULL;
    int i, id = -1;

    *type = GRETL_TYPE_NONE;

    if (K == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (!strcmp(key, "timevar_call")) {
	/* function call specifier? */
	*reserved = 1;
	if (K->matcall != NULL) {
	    ret = K->matcall;
	    *type = GRETL_TYPE_STRING;
	}
    } else {
	/* try for an input matrix specifier */
	for (i=0; i<K_MMAX; i++) {
	    if (!strcmp(key, K_input_mats[i].name)) {
		id = K_input_mats[i].sym;
		ret = (gretl_matrix *) k_input_matrix_by_id(K, id);
		if (ret != NULL) {
		    *type = GRETL_TYPE_MATRIX;
		}
		break;
	    }
	}
	if (id < 0) {
	    /* try for an output matrix */
	    gretl_matrix **pm = kalman_output_matrix(K, key);

	    if (pm != NULL) {
		*reserved = 1;
		if (*pm != NULL) {
		    ret = *pm;
		    *type = GRETL_TYPE_MATRIX;
		}
	    }
	}
	if (id < 0 && *reserved == 0) {
	    /* nothing matched yet: try scalar member */
	    ret = kalman_output_scalar(K, key);
	    if (ret != NULL) {
		*type = GRETL_TYPE_DOUBLE;
	    }
	}
    }

    if (id >= 0 && *reserved == 0) {
	/* flag the fact that @key was a kalman-reserved
	   identifier */
	*reserved = 1;
    }

    if (*reserved && ret == NULL) {
	gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
	*err = E_DATA;
    }

    return ret;
}

static int output_matrix_count (kalman *K)
{
    gretl_matrix **pm;
    int i, n = 0;
    
    for (i=0; i<K_N_OUTPUTS; i++) {
	pm = kalman_output_matrix(K, kalman_output_matrix_names[i]);
	n += (pm != NULL && *pm != NULL);
    }

    return n;
}

int print_kalman_bundle_info (void *kptr, PRN *prn)
{
    kalman *K = kptr;
    int err = 0;

    if (K == NULL) {
	pputs(prn, "Kalman struct: empty\n");
	err = E_DATA;
    } else {
	const gretl_matrix *m;
	gretl_matrix **pm;
	double *px;
	const char *name;
	int i, id;

	pputs(prn, "\nKalman input matrices\n");
	
	for (i=0; i<K_MMAX; i++) {
	    id = K_input_mats[i].sym;
	    m = k_input_matrix_by_id(K, id);
	    if (m != NULL) {
		pprintf(prn, " %s: ", K_input_mats[i].name);
		pprintf(prn, "%d x %d\n", m->rows, m->cols);
	    }
	}

	if (output_matrix_count(K) > 0) {
	    pputs(prn, "\nKalman output matrices\n");
	    for (i=0; i<K_N_OUTPUTS; i++) {
		name = kalman_output_matrix_names[i];
		pm = kalman_output_matrix(K, name);
		if (pm != NULL && *pm != NULL) {
		    m = *pm;
		    pprintf(prn, " %s: ", name);
		    pprintf(prn, "%d x %d\n", m->rows, m->cols);
		}
	    }
	}

	pputs(prn, "\nKalman scalars\n");

	for (i=0; i<K_N_SCALARS; i++) {
	    name = kalman_output_scalar_names[i];
	    pprintf(prn, " %s: ", name);
	    px = kalman_output_scalar(K, name);
	    if (px == NULL || na(*px)) {
		pputs(prn, "NA\n");
	    } else {
		pprintf(prn, "%g\n", *px);
	    }
	}

	if (K->matcall != NULL) {
	    pputs(prn, "\nKalman strings\n");
	    pprintf(prn, " timevar_call: %s\n", K->matcall);
	}
    }

    return err;
}

/* for use in context of a kalman bundle: serialize the
   information in the kalman struct to XML
*/

int kalman_serialize (void *kptr, PRN *prn)
{
    kalman *K = kptr;
    const gretl_matrix *m;
    gretl_matrix **pm;
    double *px;
    const char *name;
    int i, err = 0;

    if (K == NULL) {
	fputs("kalman_serialize: got NULL\n", stderr);
	return E_DATA;
    }

    pputs(prn, "<gretl-kalman>\n");

    for (i=0; i<K_MMAX; i++) {
	m = k_input_matrix_by_id(K, K_input_mats[i].sym);
	if (m != NULL) {
	    gretl_matrix_serialize(m, K_input_mats[i].name, prn);
	}
    }

    for (i=0; i<K_N_OUTPUTS; i++) {
	name = kalman_output_matrix_names[i];
	pm = kalman_output_matrix(K, name);
	if (pm != NULL && *pm != NULL) {
	    gretl_matrix_serialize(*pm, name, prn);
	}
    }

    for (i=0; i<=Ks_LNL; i++) {
	name = kalman_output_scalar_names[i];
	px = kalman_output_scalar(K, name);
	if (px != NULL && !na(*px)) {
	    gretl_finite_scalar_serialize(*px, name, prn);
	}
    }

    if (K->matcall != NULL) {
	gretl_string_serialize(K->matcall, "timevar_call", prn);
    }

    pputs(prn, "</gretl-kalman>\n");

    return err;
}

static int required_matrix_slot (const char *s)
{
    if (!strcmp(s, "obsy"))     return 0;
    if (!strcmp(s, "obsymat"))  return 1;
    if (!strcmp(s, "statemat")) return 2;
    if (!strcmp(s, "statevar")) return 3;
    if (!strcmp(s, "obsvar"))   return 4;
    return -1;
};

/* for use in context of a kalman bundle: deserialize the
   kalman struct from XML
*/

gretl_bundle *kalman_deserialize (void *p1, void *p2, int *err)
{
    xmlNodePtr cur, node = p1;
    xmlDocPtr doc = p2;
    gretl_matrix *Mreq[5] = {NULL};
    gretl_matrix *Mopt[K_MMAX] = {NULL};
    gretl_matrix *Mout[K_N_OUTPUTS] = {NULL};
    char *tvcall = NULL;
    double s2 = NADBL;
    double lnl = NADBL;
    int copy[5] = {0};
    int i, nmats = 0;
    int Kflags = 0;
    gretl_matrix *m;
    double x;
    char *key, *strv;
    gretl_bundle *b = NULL;

    while (node != NULL && !*err) {
        if (!xmlStrcmp(node->name, (XUC) "gretl-kalman")) {
	    cur = node->xmlChildrenNode;
	    while (cur != NULL && !*err) {
		key = (char *) xmlGetProp(cur, (XUC) "name");
		if (!xmlStrcmp(cur->name, (XUC) "gretl-matrix")) {
		    /* pick up kalman matrices */
		    m = gretl_xml_get_matrix(cur, doc, err);
		    if ((i = required_matrix_slot(key)) >= 0) {
			nmats++;
			Mreq[i] = m;
		    } else if ((i = input_matrix_slot(key)) >= 0) {
			Mopt[i] = m;
		    } else if ((i = output_matrix_slot(key)) >= 0) {
			Mout[i] = m;
		    }
		} else if (!xmlStrcmp(cur->name, (XUC) "scalar")) {
		    /* pick up kalman scalars */
		    if (gretl_xml_get_prop_as_double(cur, "value", &x)) {
			if (!strcmp(key, "diffuse") && x > 0) {
			    Kflags |= KALMAN_DIFFUSE;
			} else if (!strcmp(key, "cross") && x > 0) {
			    Kflags |= KALMAN_CROSS;
			} else if (!strcmp(key, "s2")) {
			    s2 = x;
			} else if (!strcmp(key, "lnl")) {
			    lnl = x;
			}
		    }
		} else if (!xmlStrcmp(cur->name, (XUC) "string")) {
		    /* pick up kalman strings */
		    if (!strcmp(key, "timevar_call") && 
			gretl_xml_get_prop_as_string(cur, "value", &strv)) {
			    tvcall = strv;
		    }
		}
		free(key);
		cur = cur->next;
	    }
	    break;
	}
	node = node->next;
    }

    if (nmats == 5 && !(Kflags & KALMAN_CROSS)) {
	/* drop obsvar from initialization */
	Mopt[K_R] = Mreq[4];
	Mreq[4] = NULL;
	nmats--;
    }

    if (((Kflags & KALMAN_CROSS) && nmats != 5) ||
	(!(Kflags & KALMAN_CROSS) && nmats != 4)) {
	*err = E_DATA;
    } else {
	b = kalman_bundle_new(Mreq, copy, nmats, err);
	if (b != NULL) {
	    kalman *K = gretl_bundle_get_private_data(b);
	    gretl_matrix **pm;
	    const char *name;

	    K->flags = Kflags;
	    K->s2 = s2;
	    K->loglik = lnl;

	    for (i=0; i<K_MMAX; i++) {
		if (Mopt[i] != NULL) {
		    kalman_bundle_set_matrix(K, Mopt[i], i, 0);
		}
	    }
	    for (i=0; i<K_N_OUTPUTS; i++) {
		if (Mout[i] != NULL) {
		    name = kalman_output_matrix_names[i];
		    pm = kalman_output_matrix(K, name);
		    *pm = Mout[i];
		}
	    }
	    K->matcall = tvcall;
	}
    }

    if (*err) {
	/* clean up */
	for (i=0; i<5; i++) {
	    gretl_matrix_free(Mreq[i]);
	}
	for (i=0; i<K_MMAX; i++) {
	    gretl_matrix_free(Mopt[i]);
	}
	for (i=0; i<K_N_OUTPUTS; i++) {
	    gretl_matrix_free(Mout[i]);
	}
	free(tvcall);
    }
    
    return b;
}

/* Called from gretl_bundle.c to meet the case where the user calls
   for a kalman bundle to be copied: here we create a new kalman
   struct and copy across the required elements (since they are
   not regular bundle members).
*/

gretl_bundle *kalman_bundle_copy (const gretl_bundle *src, int *err)
{
    kalman *K, *Knew;
    gretl_bundle *b = NULL;
    gretl_matrix *M[5] = {NULL};
    gretl_matrix *m, **pm, **pm1;
    const char *name;
    int copy[5] = {1, 1, 1, 1, 1};
    int i, id, k = 4;

    K = gretl_bundle_get_private_data((gretl_bundle *) src);

    if (K == NULL) {
	*err = E_DATA;
	return NULL;
    }

    M[0] = K->y;
    M[1] = K->H;
    M[2] = K->F;

    if (kalman_xcorr(K)) {
	M[3] = K->B;
	M[4] = K->C;
	k = 5;
    } else {
	M[3] = K->Q;
    }

    b = kalman_bundle_new(M, copy, k, err);

    if (*err) {
	return b;
    }

    Knew = gretl_bundle_get_private_data(b);
    Knew->flags = K->flags;

    for (i=k; i<K_MMAX && !*err; i++) {
	id = K_input_mats[i].sym;
	m = (gretl_matrix *) k_input_matrix_by_id(K, id);
	if (m != NULL) {
	    *err = kalman_bundle_set_matrix(Knew, m, i, 1);
	}
    }

    for (i=0; i<K_N_OUTPUTS && !*err; i++) {
	name = kalman_output_matrix_names[i];
	pm = kalman_output_matrix(K, name);
	if (pm != NULL && *pm != NULL) {
	    pm1 = kalman_output_matrix(Knew, name);
	    *pm1 = gretl_matrix_copy(*pm);
	}
    }

    Knew->ifc = K->ifc;
    Knew->s2 = K->s2;
    Knew->loglik = K->loglik;

    if (K->flags & KALMAN_CROSS) {
	Knew->flags |= KALMAN_CROSS;
    }

    if (K->flags & KALMAN_DIFFUSE) {
	Knew->flags |= KALMAN_DIFFUSE;
    }

    if (K->matcall != NULL) {
	Knew->matcall = gretl_strdup(K->matcall);
    }

    return b;
}

/* for use in constructing GUI bundle save menu */

char **kalman_bundle_get_matrix_names (kalman *K, int *ns)
{
    char **S = NULL;
    gretl_matrix **pm;
    const char *name;
    int i, id, err = 0;

    *ns = 0;

    for (i=0; i<K_MMAX && !err; i++) {
	id = K_input_mats[i].sym;
	if (k_input_matrix_by_id(K, id) != NULL) {
	    err = strings_array_add(&S, ns, K_input_mats[i].name);
	}
    }

    for (i=0; i<K_N_OUTPUTS && !err; i++) {
	name = kalman_output_matrix_names[i];
	pm = kalman_output_matrix(K, name);
	if (pm != NULL && *pm != NULL) {
	    err = strings_array_add(&S, ns, name);
	}
    }

    return S;
}

/* also for use in constructing GUI bundle save menu */

char **kalman_bundle_get_scalar_names (kalman *K, int *ns)
{
    char **S;

    *ns = K_N_SCALARS -1 - na(K->s2) - na(K->loglik);
     S = strings_array_new(*ns);

    if (S != NULL) {
	int i = 0;

	/* flags */
	S[i++] = gretl_strdup("cross");
	S[i++] = gretl_strdup("diffuse");

	/* actual numerical outputs */
	if (!na(K->s2)) {
	    S[i++] = gretl_strdup("s2");
	}
	if (!na(K->loglik)) {
	    S[i++] = gretl_strdup("lnl");
	}

	/* system dimensions */
	S[i++] = gretl_strdup("r");
	S[i++] = gretl_strdup("n");
	S[i++] = gretl_strdup("T");
	S[i++] = gretl_strdup("p");
    }	

    return S;
}

/* to support the nelem() function for kalman bundles */

int kalman_bundle_n_members (gretl_bundle *b)
{
    kalman *K = gretl_bundle_get_private_data(b);
    int n = 0;

    if (K != NULL) {
	int i, id;
	
	n = K_N_SCALARS;

	for (i=0; i<K_MMAX; i++) {
	    id = K_input_mats[i].sym;
	    n += (k_input_matrix_by_id(K, id) != NULL);
	}

	n += output_matrix_count(K);
    }

    return n;
}
