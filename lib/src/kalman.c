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
#include "kalman.h"

#define KDEBUG 0

struct kalman_ {
    int flags;  /* for recording any options */

    int r; /* rows of S */
    int n; /* rows of y_t */
    int k; /* rows of x_t */
    int T; /* number of observations */

    int ncoeff;    /* number of adjustable coefficients */
    int ifc;       /* include a constant (1) or not (0) */

    double SSRw;   /* \sigma_{t=1}^T e_t^{\prime} V_t^{-1} e_t */
    double sumVt;  /* \sigma_{t=1}^T ln |V_t| */
    double loglik; /* log-likelihood */

    int nonshift; /* When F is a companion matrix (e.g. in arma), the
		     number of rows of F that do something other than
		     shifting down the elements of the state vector
		  */

    /* continuously updated matrices */
    gretl_matrix *S0; /* state vector, before updating */
    gretl_matrix *S1; /* state vector, after updating */
    gretl_matrix *P0; /* MSE matrix, before updating */
    gretl_matrix *P1; /* MSE matrix, after updating */
    gretl_matrix *e;  /* one-step forecast error(s), time t */

    /* constant data matrices */
    const gretl_matrix *F; /* state transition matrix */
    const gretl_matrix *A; /* coeffs on exogenous vars, observation eqn */
    const gretl_matrix *H; /* coeffs on state variables, observation eqn */
    const gretl_matrix *Q; /* contemp covariance matrix, state eqn */
    const gretl_matrix *R; /* contemp covariance matrix, obs eqn */

    const gretl_matrix *y; /* dependent variable vector (or matrix) */
    const gretl_matrix *x; /* independent variables matrix */
    
    gretl_matrix *E;       /* forecast errors, all time-steps */

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
    gretl_matrix_free(K->S1);
    gretl_matrix_free(K->P0);
    gretl_matrix_free(K->P1);
    gretl_matrix_free(K->e);

    gretl_matrix_block_destroy(K->B);

    free(K);
}

static kalman *kalman_new_empty (void)
{
    kalman *K = malloc(sizeof *K);

    if (K != NULL) {
	K->S0 = K->S1 = NULL;
	K->P0 = K->P1 = NULL;
	K->e = NULL;
	K->B = NULL;
	K->F = K->A = K->H = NULL;
	K->Q = K->R = K->E = NULL;
	K->y = K->x = NULL;
    }

    return K;
}

static int kalman_check_dimensions (kalman *K, const gretl_matrix *S,
				    const gretl_matrix *P)
{
    K->r = gretl_matrix_rows(S);
    K->k = gretl_matrix_rows(K->A);
    K->n = gretl_matrix_cols(K->H);
    K->T = gretl_matrix_rows(K->y);

    /* S should be r x 1 */
    if (S->cols != 1) {
	fprintf(stderr, "kalman: matrix S is %d x %d, should have 1 column\n",
		S->rows, S->cols);
	return 1;
    }

    /* P should be r x r */
    if (P->rows != K->r || P->cols != K->r) {
	fprintf(stderr, "kalman: matrix P is %d x %d, should be %d x %d\n", 
		P->rows, P->cols, K->r, K->r);
	return 1;
    }

    /* F should be r x r */
    if (K->F->rows != K->r || K->F->cols != K->r) {
	fprintf(stderr, "kalman: matrix F is %d x %d, should be %d x %d\n", 
		K->F->rows, K->F->cols, K->r, K->r);
	return 1;
    }

    /* A should be k x n */
    if (K->A->cols != K->n) {
	fprintf(stderr, "kalman: matrix A has %d columns, should have %d\n", 
		K->A->cols, K->n);
	return 1;
    }    

    /* H should be r x n */
    if (K->H->rows != K->r) {
	fprintf(stderr, "kalman: matrix H has %d rows, should have %d\n", 
		K->H->rows, K->r);
	return 1;
    }    

    /* Q should be r x r */
    if (K->Q->rows != K->r || K->Q->cols != K->r) {
	fprintf(stderr, "kalman: matrix Q is %d x %d, should be %d x %d\n", 
		K->Q->rows, K->Q->cols, K->r, K->r);
	return 1;
    }

    /* R should be n x n, if present */
    if (K->R != NULL) {
	if (K->R->rows != K->n || K->R->cols != K->n) {
	    fprintf(stderr, "kalman: matrix R is %d x %d, should be %d x %d\n", 
		    K->R->rows, K->R->cols, K->n, K->n);
	    return 1;
	}
    }

    /* x should be T x (k - 1), if present (const is implicit) */
    if (K->x != NULL) {
	if (K->x->rows != K->T || K->x->cols != K->k - 1) {
	    fprintf(stderr, "kalman: matrix x is %d x %d, should be %d x %d\n", 
		    K->x->rows, K->x->cols, K->T, K->k - 1);
	    return 1;
	}
    }

    /* E should be T x n, if present */
    if (K->E != NULL) {
	if (K->E->rows != K->T || K->E->cols != K->n) {
	    fprintf(stderr, "kalman: matrix E is %d x %d, should be %d x %d\n", 
		    K->E->rows, K->E->cols, K->T, K->n);
	    return 1;
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

static int kalman_init (kalman *K, 
			const gretl_matrix *S,
			const gretl_matrix *P,
			int ncoeff, int ifc)
{
    int err = 0;

    K->B = NULL;
    K->e = NULL;

    K->flags = 0;
    K->ncoeff = ncoeff;
    K->ifc = ifc;

    K->SSRw = NADBL;
    K->sumVt = NADBL;
    K->loglik = NADBL;

    clear_gretl_matrix_err();

    if (K->S0 == NULL) {
	K->S0 = gretl_matrix_copy(S);
    }
    K->S1 = gretl_matrix_copy(S);

    if (K->P0 == NULL) {
	K->P0 = gretl_matrix_copy(P);
    }
    K->P1 = gretl_matrix_copy(P);

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

/**
 * kalman_new:
 * @S: initial state vector.
 * @P: initial MSE matrix.
 * @F: state transition matrix.
 * @A: matrix of coefficients on exogenous variables in the
 * observation equation.
 * @H: matrix of coefficients on the state variables in the
 * observation equation.
 * @Q: contemporaneous covariance matrix for the errors in the state
 * equation.
 * @R: contemporaneous covariance matrix for the errors in the 
 * observation equation (or %NULL if this is not applicable).
 * @y: T x n matrix of dependent variable(s).
 * @x: T x k matrix of exogenous variable(s).  May be %NULL if there
 * are no exogenous variables, or if there's only a constant.
 * @E: T x n matrix in which to record forecast errors (or %NULL).
 * @ncoeff: number of adjustable coefficients (used when the
 * Kalman filter is employed for estimation purposes).
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
		    int ncoeff, int ifc, int *err)
{
    kalman *K;

    *err = 0;

    if (S == NULL || P == NULL || F == NULL || 
	A == NULL || H == NULL || Q == NULL ||
	y == NULL) {
	fprintf(stderr, "A required matrix is NULL\n");
	*err = E_DATA;
	return NULL;
    }

    K = malloc(sizeof *K);
    if (K == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* just use const pointers for const matrices, don't copy */
    K->F = F;
    K->A = A;
    K->H = H;
    K->Q = Q;
    K->R = R;
    K->y = y;
    K->x = x;
    K->E = E;

    K->S0 = K->S1 = NULL;
    K->P0 = K->P1 = NULL;

    if (kalman_check_dimensions(K, S, P)) {
	fprintf(stderr, "failed on kalman_check_dimensions\n");
	*err = E_NONCONF;
	free(K);
	return NULL;
    }

    *err = kalman_init(K, S, P, ncoeff, ifc);

    if (*err) {
	kalman_free(K);
	K = NULL;
    } else {
	gretl_matrix_zero(K->e);
    }

    return K;
}

static int user_kalman_setup (kalman *K, int ncoeff, int ifc)
{
    int err = 0;

    if (K->S0 == NULL || K->P0 == NULL || K->F == NULL || 
	K->A == NULL || K->H == NULL || K->Q == NULL ||
	K->y == NULL) {
	fprintf(stderr, "A required matrix is NULL\n");
	return E_DATA;
    }

    if (kalman_check_dimensions(K, K->S0, K->P0)) {
	fprintf(stderr, "failed on kalman_check_dimensions\n");
	return E_NONCONF;
    }

    err = kalman_init(K, K->S0, K->P0, ncoeff, ifc);

    if (!err) {
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
    
    /* form (H'PH + R)^{-1} * (y - Ax - H'S) = "VE" */
    gretl_matrix_multiply(K->V, K->e, K->VE);

    if (llt != NULL) {
	/* form (y - Ax - H'S)' * (H'PH + R)^{-1} * (y - Ax - H'S) */
	gretl_matrix_multiply_mod(K->e, GRETL_MOD_TRANSPOSE,
				  K->VE, GRETL_MOD_NONE,
				  K->Tmpnn, GRETL_MOD_NONE);

	/* contribution to log-likelihood of the above -- see Hamilton
	   (1994) equation [13.4.1] page 385.
	*/
	*llt -= .5 * K->Tmpnn->val[0];
    }

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
   form A'x_t.  Note this complication: if there's a constant as well
   as other exogenous vars present, the constant does _not_ have a
   column of 1s in the x matrix (the column of 1s is implicit), though
   it _does_ have a coefficient entry in the A matrix.
*/

static void kalman_set_Ax (kalman *K, int t)
{
    double aji, xjt, axi;
    int i, j;

    /* note: in all existing applications, K->n = 1 */

    for (i=0; i<K->n; i++) {
	axi = 0.0;
	/* case j == 0 */
	if (K->ifc) {
	    axi += gretl_matrix_get(K->A, 0, i); /* \times 1.0, implicitly */
	}
	for (j=1; j<K->k; j++) {
	    aji = gretl_matrix_get(K->A, j, i);
	    xjt = gretl_matrix_get(K->x, t, j - 1);
	    axi += aji * xjt;
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
	eti = gretl_vector_get(K->e, i); /* K->e is a row vector */
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
    int Sdim = K->r * K->n;
    int Pdim = K->r * K->r;
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
	gretl_matrix_copy_values(K->Ax, K->A);
    } 

    for (t=0; t<K->T && !err; t++) {
	double llt = 0.0;

#if KDEBUG > 1
	kalman_print_state(K, t);
#endif

	/* read slice from y */
	kalman_initialize_error(K, t);

	/* and from x if applicable */
	if (K->x != NULL) {
	    kalman_set_Ax(K, t);
	}	

	/* initial matrix calculations */
	gretl_matrix_multiply(K->P0, K->H, K->PH);
	if (arma_ll(K)) {
	    K->HPH->val[0] = 0.0;
	    for (i=0; i<K->r; i++) {
		K->HPH->val[0] += K->H->val[i] * K->PH->val[i];
	    }
	    ldet = log(K->HPH->val[0]);
	    K->V->val[0] = 1.0 / K->HPH->val[0];
	} else {
	    gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
				      K->PH, GRETL_MOD_NONE,
				      K->HPH, GRETL_MOD_NONE);
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
		K->loglik = NADBL;
		err = E_NAN;
	    } else {
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
	    for (i=0; i<Sdim; i++) {
		K->S0->val[i] = K->S1->val[i];
	    }
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
		for (i=0; i<Pdim; i++) {
		    K->P0->val[i] = K->P1->val[i];
		}		
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
 * kalman_get_ncoeff:
 * @K: pointer to Kalman struct.
 * 
 * Returns: the number of adjustable coefficients
 * associated with the Kalman filter. 
 */

int kalman_get_ncoeff (const kalman *K)
{
    return K->ncoeff;
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

static gretl_matrix *attach_matrix (const char *s, int *err)
{
    char mname[VNAMELEN] = {0};
    gretl_matrix *m = NULL;

    if (sscanf(s, "%15s", mname) != 1) {
	*err = E_PARSE;
    } else {
	m = get_matrix_by_name(mname);
	if (m == NULL) {
	    *err = E_UNKVAR;
	}
    }

    return m;
}

static kalman *uK; /* user-defined Kalman struct */

/*

  kalman
    obsy y
    obsymat H
    obsex x
    obsxmat A
    obsvar R
    strmat F
    strex c ??
    strvar Q
    inistate S
    iniprec P
  end kalman

 */

/* The following is at present just a "proof of concept" thing,
   showing how we can assemble a "user" Kalman struct.  It's not
   yet hooked up in any useful way.
*/

int kalman_parse_line (const char *line, gretlopt opt)
{
    int err = 0;

    if (opt == OPT_NONE && !strncmp(line, "kalman", 6)) {
	/* starting: allocate */
	if (uK != NULL) {
	    kalman_free(uK);
	}
	uK = kalman_new_empty();
	if (uK == NULL) {
	    return E_ALLOC;
	} else {
	    return 0;
	}
    }

    if (uK == NULL) {
	/* there's no Kalman in progress */
	return E_DATA;
    }    

    if (!strncmp(line, "kalman", 6)) {
	if (opt & OPT_F) {
	    /* FIXME */
	    err = kalman_forecast(uK);
	    fprintf(stderr, "kalman_forecast: err = %d\n", err);
	    return err;
	} else {
	    return E_DATA;
	}
    }
	    
    if (!strncmp(line, "end ", 4)) {
	/* "end kalman": complete the set-up */
	err = user_kalman_setup(uK, 0, 1); /* FIXME */
    } else if (!strncmp(line, "obsy ", 5)) {
	uK->y = attach_matrix(line + 5, &err);
    } else if (!strncmp(line, "obsymat ", 8)) {
	uK->H = attach_matrix(line + 8, &err);
    } else if (!strncmp(line, "obsex ", 6)) {
	uK->x = attach_matrix(line + 6, &err);
    } else if (!strncmp(line, "obsxmat ", 8)) {
	uK->A = attach_matrix(line + 8, &err);
    } else if (!strncmp(line, "obsvar ", 7)) {
	uK->R = attach_matrix(line + 7, &err);
    } else if (!strncmp(line, "strmat ", 7)) {
	uK->F = attach_matrix(line + 7, &err);
    } else if (!strncmp(line, "strvar ", 7)) {
	uK->Q = attach_matrix(line + 7, &err);
    } else if (!strncmp(line, "inistate ", 9)) {
	uK->S0 = attach_matrix(line + 9, &err);
    } else if (!strncmp(line, "iniprec ", 8)) {
	uK->P0 = attach_matrix(line + 8, &err);
    } else {
	err = E_PARSE;
    }

    fprintf(stderr, "kalman_parse_line: '%s' -> err = %d\n", line, err);

    if (err && uK != NULL) {
	kalman_free(uK);
	uK = NULL;
    }

    return err;
}



