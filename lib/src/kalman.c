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
#include "kalman.h"

#define KDEBUG 0

/* notation (cf. James Hamilton):

   State:        S_{t+1} = F*S_t + v{t+1},     var(v_t) = Q

   Observation:  y_t = A'*x_t + H*S_t + w_t,   var(w_t) = R

 */

struct kalman_ {
    int flags;   /* for recording any options */
    int fnlevel; /* level of function execution */

    int r;  /* rows of S = number of elements in state */
    int n;  /* columns of y = number of observables */
    int k;  /* columns of A = number of exogenous vars in obs eqn */
    int T;  /* rows of y = number of observations */

    int ifc; /* boolean: obs equation includes an implicit constant? */

    double SSRw;   /* \sigma_{t=1}^T e_t^{\prime} V_t^{-1} e_t */
    double sumVt;  /* \sigma_{t=1}^T ln |V_t| */
    double loglik; /* log-likelihood */

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

    /* constant data matrices */
    const gretl_matrix *F; /* r x r: state transition matrix */
    const gretl_matrix *A; /* k x n: coeffs on exogenous vars, obs eqn */
    const gretl_matrix *H; /* r x n: coeffs on state variables, obs eqn */
    const gretl_matrix *Q; /* r x r: contemp covariance matrix, state eqn */
    const gretl_matrix *R; /* n x n: contemp covariance matrix, obs eqn */
    const gretl_matrix *y; /* T x n: dependent variable vector (or matrix) */
    const gretl_matrix *x; /* T x k: independent variables matrix */
    const gretl_matrix *S; /* r x 1: initial state vector */
    const gretl_matrix *P; /* r x r: initial precision matrix */

    /* optional run-time export matrices */
    gretl_matrix *E;       /* T x n: forecast errors, all time-steps */
    gretl_matrix *Sigma;   /* T x np: variance, all time-steps */
    gretl_matrix *State;   /* T x r: state vector, all time-steps */
    gretl_matrix *LL;      /* T x 1: loglikelihood, all time-steps */
    
    /* workspace matrices */
    gretl_matrix_block *B; /* holder for the following */
    gretl_matrix *PH;
    gretl_matrix *HPH;
    gretl_matrix *FPH;
    gretl_matrix *V;
    gretl_matrix *VE;
    gretl_matrix *PHV;
    gretl_matrix *Ax;
    gretl_matrix *Tmpnn;
    gretl_matrix *Tmprr;
    gretl_matrix *Tmprr_2a;
    gretl_matrix *Tmprr_2b;
};

#define arma_ll(K) (K->flags & KALMAN_ARMA_LL)

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

    gretl_matrix_block_destroy(K->B);

    free(K);
}

static kalman *kalman_new_empty (void)
{
    kalman *K = malloc(sizeof *K);

    if (K != NULL) {
	K->S = K->P = NULL;
	K->S0 = K->S1 = NULL;
	K->P0 = K->P1 = NULL;
	K->e = NULL;
	K->B = NULL;
	K->F = K->A = K->H = NULL;
	K->Q = K->R = NULL;
	K->E = K->Sigma = K->State = K->LL = NULL;
	K->y = K->x = NULL;
	K->flags = 0;
	K->fnlevel = 0;
    }

    return K;
}

static int kalman_check_dimensions (kalman *K)
{
    if (K->r < 1 || K->n < 1 || K->T < 2) {
	/* the state and observation vectors must have at least one
	   element, and there must be at least two observations */
	return E_DATA;
    }

    /* S should be r x 1 (S->rows defines r) */
    if (K->S->cols != 1) {
	fprintf(stderr, "kalman: matrix S is %d x %d, should have 1 column\n",
		K->S->rows, K->S->cols);
	return E_NONCONF;
    }

    /* P should be r x r */
    if (K->P->rows != K->r || K->P->cols != K->r) {
	fprintf(stderr, "kalman: matrix P is %d x %d, should be %d x %d\n", 
		K->P->rows, K->P->cols, K->r, K->r);
	return E_NONCONF;
    }

    /* F should be r x r */
    if (K->F->rows != K->r || K->F->cols != K->r) {
	fprintf(stderr, "kalman: matrix F is %d x %d, should be %d x %d\n", 
		K->F->rows, K->F->cols, K->r, K->r);
	return E_NONCONF;
    }

    /* H should be r x n (where y defines n) */
    if (K->H->rows != K->r || K->H->cols != K->n) {
	fprintf(stderr, "kalman: matrix H is %d x %d, should be %d x %d\n", 
		K->H->rows, K->H->cols, K->r, K->n);
	return E_NONCONF;
    }    

    /* Q should be r x r */
    if (K->Q->rows != K->r || K->Q->cols != K->r) {
	fprintf(stderr, "kalman: matrix Q is %d x %d, should be %d x %d\n", 
		K->Q->rows, K->Q->cols, K->r, K->r);
	return E_NONCONF;
    }

    /* A should be k x n, if present (A->rows defines k; n is defined by y) */
    if (K->A != NULL) {
	if (K->A->cols != K->n) {
	    fprintf(stderr, "kalman: matrix A has %d columns, should have %d\n", 
		    K->A->cols, K->n);
	    return E_NONCONF;
	}
    } else if (K->x != NULL) {
	/* A is NULL => can't have a non-NULL x */
	return E_NONCONF;
    }

    K->ifc = 0;

    /* x should have T rows to match y; and it should have either k or k - 1
       columns (the latter case indicating an implicit const) */
    if (K->x != NULL) {
	if (K->x->rows != K->T) {
	    fprintf(stderr, "kalman: matrix x has %d rows, should have %d\n", 
		    K->x->rows, K->T);
	    return E_NONCONF;
	} else if (K->x->cols != K->k && K->x->cols != K->k - 1) {
	    fprintf(stderr, "kalman: matrix x has %d columns, should have %d or %d\n", 
		    K->x->cols, K->k, K->k - 1);
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
	gretl_errmsg_set(_("kalman: a required matrix is missing"));
	return E_DATA;
    }

    /* R should be n x n, if present */
    if (K->R != NULL) {
	if (K->R->rows != K->n || K->R->cols != K->n) {
	    fprintf(stderr, "kalman: matrix R is %d x %d, should be %d x %d\n", 
		    K->R->rows, K->R->cols, K->n, K->n);
	    return E_NONCONF;
	}
    }

    /* E should be T x n, if present */
    if (K->E != NULL) {
	if (K->E->rows != K->T || K->E->cols != K->n) {
	    fprintf(stderr, "kalman: matrix E is %d x %d, should be %d x %d\n", 
		    K->E->rows, K->E->cols, K->T, K->n);
	    return E_NONCONF;
	}
    }

    /* Sigma should be T x np, if present */
    if (K->Sigma != NULL) {
	int np = (K->r * K->r + K->r) / 2;;

	if (K->Sigma->rows != K->T || K->Sigma->cols != np) {
	    fprintf(stderr, "kalman: matrix Sigma is %d x %d, should be %d x %d\n", 
		    K->Sigma->rows, K->Sigma->cols, K->T, np);
	    return E_NONCONF;
	}
    }

    /* State should be T x r, if present */
    if (K->State != NULL) {
	if (K->State->rows != K->T || K->State->cols != K->r) {
	    fprintf(stderr, "kalman: matrix State is %d x %d, should be %d x %d\n", 
		    K->State->rows, K->State->cols, K->T, K->r);
	    return E_NONCONF;
	}
    }

    /* LL should be T x 1, if present */
    if (K->LL != NULL) {
	if (K->LL->rows != K->T || K->LL->cols != 1) {
	    fprintf(stderr, "kalman: matrix LL is %d x %d, should be %d x %d\n", 
		    K->LL->rows, K->LL->cols, K->T, 1);
	    return E_NONCONF;
	}
    }    

    return 0;
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
    K->sumVt = NADBL;
    K->loglik = NADBL;

    clear_gretl_matrix_err();

    K->S0 = gretl_matrix_copy(K->S);
    K->S1 = gretl_matrix_copy(K->S);

    K->P0 = gretl_matrix_copy(K->P);
    K->P1 = gretl_matrix_copy(K->P);

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
				  &K->V,   K->n, K->n, /* (H'*P*H + R)^{-1} */
				  &K->VE,  K->n, 1,    /* (H'*P*H + R)^{-1} * E */
				  &K->PHV, K->r, K->n, /* P*H*V */
				  &K->Ax,  K->n, 1,    /* A'*x at obs t */
				  &K->Tmpnn, K->n, K->n,
				  &K->Tmprr, K->r, K->r,
				  &K->Tmprr_2a, K->r, K->r,
				  &K->Tmprr_2b, K->r, K->r,
				  NULL);

    if (K->B == NULL) {
	err = E_ALLOC;
    }

    return err;
}

static void kalman_set_dimensions (kalman *K)
{
    K->r = gretl_matrix_rows(K->S); /* S->rows defines r */
    K->k = gretl_matrix_rows(K->A); /* A->rows defines k */
    K->T = gretl_matrix_rows(K->y); /* y->rows defines T */
    K->n = gretl_matrix_cols(K->y); /* y->cols defines n */
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

    if (S == NULL || P == NULL || F == NULL || 
	H == NULL || Q == NULL || y == NULL) {
	gretl_errmsg_set(_("kalman: a required matrix is missing"));
	*err = E_DATA;
	return NULL;
    }

    K = kalman_new_empty();
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
    K->S = S;
    K->P = P;

    K->E = E;

    kalman_set_dimensions(K);

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

static int user_kalman_setup (kalman *K)
{
    int err = 0;

    if (K->S == NULL || K->P == NULL || K->F == NULL || 
	K->H == NULL || K->Q == NULL || K->y == NULL) {
	gretl_errmsg_set(_("kalman: a required matrix is missing"));
	return E_DATA;
    }

    kalman_set_dimensions(K);

    err = kalman_check_dimensions(K);
    if (err) {
	fprintf(stderr, "failed on kalman_check_dimensions\n");
	return err;
    }

    err = kalman_init(K);

    if (!err) {
	K->flags = KALMAN_USER;
	gretl_matrix_zero(K->e);
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

#if 0

/* no-brainer, just let the compiler handle things */

static inline int multiply_by_F (kalman *K, const gretl_matrix *A, 
				 gretl_matrix *B, int postmult)
{
    if (postmult) {
	gretl_matrix_multiply_mod(A, GRETL_MOD_NONE,
				  K->F, GRETL_MOD_TRANSPOSE,
				  B, GRETL_MOD_NONE);
    } else {
	gretl_matrix_multiply(K->F, A, B);
    }

    return 0;
}

#else

/* carefully optimized */

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
	int i, j, c;
	int r1 = K->nonshift;
	int r2 = K->r - r1;
	gretl_matrix *topF = K->Tmprr_2a;
	gretl_matrix *top = K->Tmprr_2b;
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

#endif

/* Simplified version of the function below, for ARMA:
   in this case we have H = (r x 1) and S = (r x 1),
   although F and P are (r x r).
*/

static int kalman_arma_iter_1 (kalman *K)
{
    double ve;
    int i, err = 0;

    /* write F*S into S+ */
    err += multiply_by_F(K, K->S0, K->S1, 0);

    /* form E = y - A'x - H'S */
    K->e->val[0] -= K->Ax->val[0];
    for (i=0; i<K->r; i++) {
	K->e->val[0] -= K->H->val[i] * K->S0->val[i];
    }
   
    /* form (H'PH + R)^{-1} * (y - Ax - H'S) = "ve" */
    ve = K->e->val[0] * K->V->val[0];

    /* form FPH */
    err += multiply_by_F(K, K->PH, K->FPH, 0);

    /* form FPH * (H'PH + R)^{-1} * (y - A'x - H'S),
       and complete calculation of S+ */
    for (i=0; i<K->r; i++) {
	K->S1->val[i] += K->FPH->val[i] * ve;
    }

    return err;
}

/* Hamilton (1994) equation [13.2.23] page 381, in simplified notation:

   S+ = FS + FPH(H'PH + R)^{-1} * (y - A'x - H'S) 

   "S" is Hamilton's \xi (state vector)
*/

static int kalman_iter_1 (kalman *K, double *llt)
{
    int err = 0;

    /* write F*S into S+ */
    err += multiply_by_F(K, K->S0, K->S1, 0);

    /* form E = y - A'x - H'S */
    err += gretl_matrix_subtract_from(K->e, K->Ax);
    gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
			      K->S0, GRETL_MOD_NONE,
			      K->e, GRETL_MOD_DECUMULATE);
    
    /* form (H'PH + R)^{-1} * (y - A'x - H'S) = "VE" */
    gretl_matrix_multiply(K->V, K->e, K->VE);

    /* form (y - A'x - H'S)' * (H'PH + R)^{-1} * (y - A'x - H'S) */
    gretl_matrix_multiply_mod(K->e, GRETL_MOD_TRANSPOSE,
			      K->VE, GRETL_MOD_NONE,
			      K->Tmpnn, GRETL_MOD_NONE);

    /* contribution to log-likelihood of the above -- see Hamilton
       (1994) equation [13.4.1] page 385.
    */
    *llt -= .5 * K->Tmpnn->val[0];

    /* form FPH */
    err += multiply_by_F(K, K->PH, K->FPH, 0);

    /* form FPH * (H'PH + R)^{-1} * (y - A'x - H'S) 
       and add to S+ */
    gretl_matrix_multiply_mod(K->FPH, GRETL_MOD_NONE,
			      K->VE, GRETL_MOD_NONE,
			      K->S1, GRETL_MOD_CUMULATE);

    return err;
}

/* Hamilton (1994) equation [13.2.22] page 380, in simplified notation:

   P+ = F[P - PH(H'PH + R)^{-1}H'P]F' + Q 
*/

static int kalman_iter_2 (kalman *K)
{
    int err = 0;

    /* form P - PH(H'PH + R)^{-1}H'P */
    err += gretl_matrix_qform(K->PH, GRETL_MOD_NONE, K->V,
			      K->P0, GRETL_MOD_DECUMULATE);

    /* pre-multiply by F, post-multiply by F' */
    if (K->nonshift == K->r) {
	gretl_matrix_qform(K->F, GRETL_MOD_NONE, K->P0,
			   K->P1, GRETL_MOD_NONE);
    } else {
	err += multiply_by_F(K, K->P0, K->Tmprr, 0);
	err += multiply_by_F(K, K->Tmprr, K->P1, 1);
	gretl_matrix_xtr_symmetric(K->P1);
    }

    /* add Q */
    err += gretl_matrix_add_to(K->P1, K->Q);

    return err;
}

#if KDEBUG > 1
static void kalman_print_state (kalman *K, int t)
{
    int j;

    /* if (t > 5) return; */

    fprintf(stderr, "Iteration %d:\n", t);

    for (j=0; j<K->n; j++) {
	fprintf(stderr, "y[%d] = %.10g, err[%d] = %.10g\n", j, 
		gretl_matrix_get(K->y, t, j), 
		j, gretl_vector_get(K->e, j));
    }

    gretl_matrix_print(K->S0, "K->S0");
    gretl_matrix_print(K->P0, "K->P0");
}
#endif

static void kalman_record_state (kalman *K, int t)
{
    int i, j, m = 0;
    double x;

    if (K->State != NULL) {
	for (i=0; i<K->r; i++) {
	    x = gretl_vector_get(K->S0, i);
	    gretl_matrix_set(K->State, t, i, x);
	}
    } 

    if (K->Sigma != NULL) {
	for (i=0; i<K->r; i++) {
	    for (j=i; j<K->r; j++) {
		x = gretl_matrix_get(K->P0, i, j);
		gretl_matrix_set(K->Sigma, t, m++, x);
	    }
	}
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

#define kalman_ll_average(K) (K->flags & KALMAN_AVG_LL)

/* Read from the appropriate row of x (T x k) and multiply by A' to
   form A'x_t.  Note: the flag K->ifc is used to indicate that the
   observation equation has an implicit constant, with an entry in
   the A matrix (the first) but no explicit entry in the x matrix.
*/

static void kalman_set_Ax (kalman *K, int t)
{
    double xjt, axi;
    int i, j;

    for (i=0; i<K->n; i++) {
	axi = 0.0;
	for (j=0; j<K->k; j++) {
	    if (K->ifc) {
		xjt = (j == 0)? 1.0 : gretl_matrix_get(K->x, t, j - 1);
	    } else {
		xjt = gretl_matrix_get(K->x, t, j);
	    }
	    axi += xjt * gretl_matrix_get(K->A, j, i);
	}
	gretl_vector_set(K->Ax, i, axi);
    }
}

/* read from the appropriate row of y (T x n) and transcribe to
   the current e (n x 1)
*/

static void
kalman_initialize_error (kalman *K, int t)
{
    int i;

    for (i=0; i<K->n; i++) {
	K->e->val[i] = gretl_matrix_get(K->y, t, i);
    }
}

/* read the current forecast error and write it into the appropriate
   row of the recorder matrix, K->E
*/

static void kalman_record_error (kalman *K, int t)
{
    double eti;
    int i;

    for (i=0; i<K->n; i++) {
	eti = gretl_vector_get(K->e, i);
	gretl_matrix_set(K->E, t, i, eti);
    }
}

/* increment the "weighted SSR": add e'_t V_t^{-1} e_t */

static int kalman_incr_S (kalman *K)
{
    double x;
    int err = 0;

    x = gretl_scalar_qform(K->e, K->V, &err);
    if (!err) {
	K->SSRw += x;
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
    int i, t, err = 0;

#if KDEBUG
    fprintf(stderr, "kalman_forecast: T = %d\n", K->T);
#endif  

    if (K->nonshift < 0) {
	K->nonshift = count_nonshifts(K->F);
    }

    if (arma_ll(K)) {
	K->SSRw = 0.0;
	K->sumVt = 0.0;
    }

    K->loglik = 0.0;

    if (K->x == NULL) {
	/* no exogenous vars */
	if (K->A != NULL) {
	    /* implicit const case: A is 1 x n and A'x is n x 1 */
	    gretl_vector_copy_values(K->Ax, K->A);
	} else {
	    gretl_matrix_zero(K->Ax);
	}
    } 

    for (t=0; t<K->T && !err; t++) {
	double llt = 0.0;

#if KDEBUG > 1
	kalman_print_state(K, t);
#endif

	if (K->State != NULL || K->Sigma != NULL) {
	    kalman_record_state(K, t);
	}

	/* read slice from y */
	kalman_initialize_error(K, t);

	/* and from x if applicable */
	if (K->x != NULL) {
	    kalman_set_Ax(K, t);
	}	

	/* initial matrix calculations: form PH and H'PH 
	   (note that we need PH later) */
	gretl_matrix_multiply(K->P0, K->H, K->PH);
	if (arma_ll(K)) {
	    K->HPH->val[0] = 0.0;
	    for (i=0; i<K->r; i++) {
		K->HPH->val[0] += K->H->val[i] * K->PH->val[i];
	    }
	    ldet = log(K->HPH->val[0]);
	    K->V->val[0] = 1.0 / K->HPH->val[0];
	} else {
	    gretl_matrix_qform(K->H, GRETL_MOD_TRANSPOSE,
			       K->P0, K->HPH, GRETL_MOD_NONE);
	    if (K->R != NULL) {
		gretl_matrix_add_to(K->HPH, K->R);
	    }
	    gretl_matrix_copy_values(K->V, K->HPH);
	    err = gretl_invert_symmetric_matrix2(K->V, &ldet);
	    if (err) {
		fprintf(stderr, "kalman_forecast: failed to invert V\n");
	    }
	}	    

	/* likelihood bookkeeping */
	if (err) {
	    K->loglik = NADBL;
	    break;
	} else if (arma_ll(K)) {
	    K->sumVt += ldet;
	} else {
	    llt = -(K->n / 2.0) * LN_2_PI - .5 * ldet;
	}

	/* first stage of dual iteration */
	if (arma_ll(K)) {
	    err = kalman_arma_iter_1(K);
	    if (!err) {
		err = kalman_incr_S(K);
	    }
	} else {
	    err = kalman_iter_1(K, &llt);
	    if (na(llt)) {
		if (K->LL != NULL) {
		    gretl_vector_set(K->LL, t, 0.0 / 0.0);
		}
		K->loglik = NADBL;
		err = E_NAN;
	    } else {
		if (K->LL != NULL) {
		    gretl_vector_set(K->LL, t, llt);
		}
		K->loglik += llt;
		if (isnan(K->loglik)) {
		    err = E_NAN;
		}
	    }
	}
	
	/* record forecast errors if wanted */
	if (!err && K->E != NULL) {
	    kalman_record_error(K, t);
	}

	/* update state vector */
	if (!err) {
	    gretl_matrix_copy_values(K->S0, K->S1);
	}

	if (!err && update_P) {
	    /* second stage of dual iteration */
	    err = kalman_iter_2(K);
	}

	if (!err) {
	    /* update MSE matrix, if needed */
	    if (arma_ll(K) && update_P && t > 20) {
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

    if (isnan(K->loglik) || isinf(K->loglik)) {
	K->loglik = NADBL;
    } 

    if (arma_ll(K) && !na(K->loglik)) {
	if (kalman_ll_average(K)) {
	    K->loglik = -0.5 * 
		(LN_2_PI + 1 + log(K->SSRw / K->T) + K->sumVt / K->T);
	} else {
	    double k = -(K->T / 2.0) * LN_2_PI;

	    K->loglik = k - (K->T / 2.0) * (1.0 + log(K->SSRw / K->T))
		- 0.5 * K->sumVt;
	}
	if (isnan(K->loglik) || isinf(K->loglik)) {
	    K->loglik = NADBL;
	    err = E_NAN;
	} 
    }

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

/* user-defined Kalman apparatus below */

#define KNMAT 10 /* the (max) number of user-attached matrices */

/* symbolic identifiers for attached matrix positions */

enum {
    K_y = 0, /* obsy */
    K_H,     /* obsymat */
    K_x,     /* obsex */
    K_A,     /* obsxmat */
    K_R,     /* obsvar */
    K_F,     /* strmat */
    K_Q,     /* strvar */
    K_S,     /* inistate */
    K_P,     /* iniprec */
    K_E      /* fcast errs */
};

typedef struct user_kalman_ user_kalman;

struct user_kalman_ {
    kalman *K;
    int fnlevel;
    char **mnames;
};

/* At present there can be at most one user_kalman at each
   depth of function execution */

static user_kalman **uK; /* user-defined Kalman structs */
static int n_uK;         /* number of same */

/* given the name of a user-defined matrix, attach it to the
   user Kalman struct and record its name */

static gretl_matrix *attach_matrix (user_kalman *u, const char *s, 
				    int i, int *err)
{
    char mname[VNAMELEN] = {0};
    gretl_matrix *m = NULL;

    if (u->mnames == NULL) {
	u->mnames = strings_array_new_with_length(KNMAT, VNAMELEN);
	if (u->mnames == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	}
    }

    if (sscanf(s, "%15s", mname) != 1) {
	*err = E_PARSE;
    } else {
	m = get_matrix_by_name(mname);
	if (m == NULL) {
	    gretl_errmsg_sprintf(_("'%s': no such matrix\n"), mname);
	    *err = E_UNKVAR;
	} else {
	    strcpy(u->mnames[i], mname);
	}
    }

    return m;
}

static gretl_matrix *kalman_retrieve_matrix (const char *name,
					     int level, int cfd)
{
    gretl_matrix *m;

    /* We allow for the possibility that one or more of the matrices
       that are attached to a kalman struct have been temporarily
       "promoted", via matrix-pointer arguments, to the level of a
       user-defined function.
    */

    m = get_matrix_by_name_at_level(name, level);
    if (m == NULL && level < cfd) {
	m = get_matrix_by_name_at_level(name, cfd);
    }

    return m;
}

/* When (re-)running a user-defined filter, check that no relevant
   matrices have disappeared or been resized.  In addition,
   re-initalize the state and precision.
*/

static int user_kalman_recheck_matrices (user_kalman *u)
{
    int cfd = gretl_function_depth();
    kalman *K = u->K;
    int err = 0;

    K->y = kalman_retrieve_matrix(u->mnames[K_y], u->fnlevel, cfd);
    K->H = kalman_retrieve_matrix(u->mnames[K_H], u->fnlevel, cfd);
    K->x = kalman_retrieve_matrix(u->mnames[K_x], u->fnlevel, cfd);
    K->A = kalman_retrieve_matrix(u->mnames[K_A], u->fnlevel, cfd);
    K->R = kalman_retrieve_matrix(u->mnames[K_R], u->fnlevel, cfd);
    K->F = kalman_retrieve_matrix(u->mnames[K_F], u->fnlevel, cfd);
    K->Q = kalman_retrieve_matrix(u->mnames[K_Q], u->fnlevel, cfd);
    K->S = kalman_retrieve_matrix(u->mnames[K_S], u->fnlevel, cfd);
    K->P = kalman_retrieve_matrix(u->mnames[K_P], u->fnlevel, cfd);

    if (K->S == NULL || K->P == NULL || K->F == NULL || 
	K->H == NULL || K->Q == NULL || K->y == NULL) {
	gretl_errmsg_set(_("kalman: a required matrix is missing"));
	err = E_DATA;
    } else if (gretl_matrix_rows(K->S) != K->r ||
	       gretl_matrix_rows(K->A) != K->k ||
	       gretl_matrix_rows(K->y) != K->T ||
	       gretl_matrix_cols(K->y) != K->n) {
	err = E_NONCONF;
    } else {
	err = kalman_check_dimensions(K);
	if (!err) {
	    gretl_matrix_copy_values(K->S0, K->S);
	    gretl_matrix_copy_values(K->P0, K->P);
	}
    }
    
    return err;
}

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

void destroy_user_kalman (void)
{
    int i, lev = gretl_function_depth();

    for (i=0; i<n_uK; i++) {
	if (uK[i] != NULL && uK[i]->fnlevel == lev) {
	    kalman_free(uK[i]->K);
	    free_strings_array(uK[i]->mnames, KNMAT);
	    free(uK[i]);
	    uK[i] = NULL;
	    return;
	}
    }
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
	    u->K = kalman_new_empty();
	    if (u->K == NULL) {
		free(u);
		err = E_ALLOC;
	    } else {
		/* finalize and place on stack */
		u->mnames = NULL;
		u->fnlevel = gretl_function_depth();
		u->K->fnlevel = u->fnlevel;
		uK[n] = u;
	    }
	}
    }
	
    return err;
}

/* Given a line that is part of a "kalman ... end kalman" build process,
   or either of the stand-alone lines "kalman --filter" or "kalman
   --delete", parse and respond appropriately.
*/

int kalman_parse_line (const char *line, gretlopt opt)
{
    user_kalman *u;
    int err = 0;

    if (opt == OPT_NONE && !strncmp(line, "kalman", 6)) {
	/* starting: allocate */
	return add_user_kalman();
    }

    u = get_user_kalman(gretl_function_depth());

    if (u == NULL) {
	/* there's no Kalman in progress */
	return E_DATA;
    }    

    if (!strncmp(line, "kalman", 6)) {
	if (opt & OPT_F) {
	    /* --filter option */
	    err = user_kalman_recheck_matrices(u);
	    if (!err) {
		err = kalman_forecast(u->K);
	    }
	    return err;
	} else if (opt & OPT_D) {
	    /* --delete option */
	    destroy_user_kalman();
	    return 0;
	} else {
	    /* huh? */
	    return E_DATA;
	}
    }
	    
    if (!strncmp(line, "end ", 4)) {
	/* "end kalman": complete the set-up */
	err = user_kalman_setup(u->K);
    } else if (!strncmp(line, "obsy ", 5)) {
	u->K->y = attach_matrix(u, line + 5, K_y, &err);
    } else if (!strncmp(line, "obsymat ", 8)) {
	u->K->H = attach_matrix(u, line + 8, K_H, &err);
    } else if (!strncmp(line, "obsex ", 6)) {
	u->K->x = attach_matrix(u, line + 6, K_x, &err);
    } else if (!strncmp(line, "obsxmat ", 8)) {
	u->K->A = attach_matrix(u, line + 8, K_A, &err);
    } else if (!strncmp(line, "obsvar ", 7)) {
	u->K->R = attach_matrix(u, line + 7, K_R, &err);
    } else if (!strncmp(line, "strmat ", 7)) {
	u->K->F = attach_matrix(u, line + 7, K_F, &err);
    } else if (!strncmp(line, "strvar ", 7)) {
	u->K->Q = attach_matrix(u, line + 7, K_Q, &err);
    } else if (!strncmp(line, "inistate ", 9)) {
	u->K->S = attach_matrix(u, line + 9, K_S, &err);
    } else if (!strncmp(line, "iniprec ", 8)) {
	u->K->P = attach_matrix(u, line + 8, K_P, &err);
    } else if (!strncmp(line, "errors ", 7)) {
	u->K->E = attach_matrix(u, line + 7, K_E, &err);
    } else {
	err = E_PARSE;
    }

#if 0
    fprintf(stderr, "kalman_parse_line: '%s' -> err = %d\n", line, err);
#endif

    if (err) {
	destroy_user_kalman();
    }

    return err;
}

static gretl_matrix *attach_export_matrix (const char *mname, int *err)
{
    gretl_matrix *m;

    if (mname == NULL || !strcmp(mname, "null")) {
	return NULL;
    }

    m = get_matrix_by_name(mname);
    
    if (m == NULL) {
	gretl_errmsg_sprintf(_("'%s': no such matrix\n"), mname);
	*err = E_UNKVAR;
    }

    return m;
}

/* FIXME something is broken with this apparatus */

int user_kalman_run (const char *E, const char *S, const char *P,
		     const char *L, int *err)
{
    user_kalman *u = get_user_kalman(-1);
    int ret = 0;

    if (u == NULL) {
	*err = E_DATA;
    }

    if (*err) {
	u->K->E = attach_export_matrix(E, err);
    }

    if (!*err) {
	u->K->State = attach_export_matrix(S, err);
    }

    if (!*err) {
	u->K->Sigma = attach_export_matrix(P, err);
    } 

    if (!*err) {
	u->K->LL = attach_export_matrix(L, err);
    } 

    if (!*err) {
	*err = user_kalman_recheck_matrices(u);
	if (!*err) {
	    *err = kalman_forecast(u->K);
	}
    }

    if (*err != 0 && *err != E_NAN) {
	ret = 1;
    }

    return ret;    
}

/**
 * user_kalman_get_loglik:
 * @K: pointer to Kalman struct.
 * 
 * Retrieves the log-likelhood calculated via the last run of 
 * kalman --forecast, if applicable.
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
