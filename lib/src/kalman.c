/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "libgretl.h"
#include "kalman.h"

#define KDEBUG 0

struct kalman_ {
    int r; /* rows of S */
    int n; /* rows of y_t */
    int k; /* rows of x_t */

    double loglik; /* log-likelihood */

    /* continuously updated matrices */
    gretl_matrix *S0; /* state vector, before updating */
    gretl_matrix *S1; /* state vector, after updating */
    gretl_matrix *P0; /* MSE matrix, before updating */
    gretl_matrix *P1; /* MSE matrix, after updating */
    gretl_matrix *E;  /* one-step forecast error(s), time t */

    /* constant data matrices */
    const gretl_matrix *F; /* state transition matrix */
    const gretl_matrix *A; /* coeffs on exogenous vars, observation eqn */
    const gretl_matrix *H; /* coeffs on state variables, observation eqn */
    const gretl_matrix *Q; /* contemp covariance matrix, state eqn */
    const gretl_matrix *R; /* contemp covariance matrix, obs eqn */

    /* workspace matrices (may be able to economize on these?) */
    gretl_matrix *PH;
    gretl_matrix *HPH;
    gretl_matrix *FPH;
    gretl_matrix *V;
    gretl_matrix *VE;
    gretl_matrix *PHV;
    gretl_matrix *Ax;
    gretl_matrix *Tmprn;
    gretl_matrix *Tmpnn;
    gretl_matrix *Tmprr;
};

void kalman_free (kalman *K)
{
    gretl_matrix_free(K->S0);
    gretl_matrix_free(K->S1);
    gretl_matrix_free(K->P0);
    gretl_matrix_free(K->P1);
    gretl_matrix_free(K->E);

    gretl_matrix_free(K->PH);
    gretl_matrix_free(K->HPH);
    gretl_matrix_free(K->FPH);
    gretl_matrix_free(K->V);
    gretl_matrix_free(K->VE);
    gretl_matrix_free(K->PHV);
    gretl_matrix_free(K->Ax);

    gretl_matrix_free(K->Tmprn);
    gretl_matrix_free(K->Tmpnn);
    gretl_matrix_free(K->Tmprr);

    free(K);
}

static int 
kalman_check_dimensions (kalman *K, 
			 const gretl_matrix *S, const gretl_matrix *P,
			 const gretl_matrix *F, const gretl_matrix *A,
			 const gretl_matrix *H, const gretl_matrix *Q,
			 const gretl_matrix *R)
{
    K->r = gretl_matrix_rows(S);
    K->k = gretl_matrix_rows(A);
    K->n = gretl_matrix_cols(H);

    /* S should be r x 1 */
    if (gretl_matrix_cols(S) != 1) {
	return 1;
    }

    /* P should be r x r */
    if (gretl_matrix_rows(P) != K->r ||
	gretl_matrix_cols(P) != K->r) {
	return 1;
    }

    /* F should be r x r */
    if (gretl_matrix_rows(F) != K->r ||
	gretl_matrix_cols(F) != K->r) {
	return 1;
    }

    /* A should be k x n */
    if (gretl_matrix_cols(A) != K->n) {
	return 1;
    }    

    /* H should be r x n */
    if (gretl_matrix_rows(H) != K->r) {
	return 1;
    }    

    /* Q should be r x r */
    if (gretl_matrix_rows(Q) != K->r ||
	gretl_matrix_cols(Q) != K->r) {
	return 1;
    }

    /* R should be n x n, if present */
    if (R != NULL) {
	if (gretl_matrix_rows(R) != K->n ||
	    gretl_matrix_cols(R) != K->n) {
	    return 1;
	}
    }

    return 0;
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
		    const gretl_matrix *R, int *err)
{
    kalman *K;

    *err = 0;

    K = malloc(sizeof *K);
    if (K == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (kalman_check_dimensions(K, S, P, F, A, H, Q, R)) {
	*err = E_NONCONF;
	free(K);
	return NULL;
    }

    K->loglik = NADBL;

    K->S0 = NULL;
    K->S1 = NULL;
    K->P0 = NULL;
    K->P1 = NULL;
    K->E = NULL;

    K->PH = NULL;
    K->HPH = NULL;
    K->FPH = NULL;
    K->V = NULL;
    K->VE = NULL;
    K->PHV = NULL;
    K->Ax = NULL;

    K->Tmprn = NULL;
    K->Tmpnn = NULL;
    K->Tmprr = NULL;

    K->S0 = gretl_matrix_copy(S);
    K->S1 = gretl_matrix_copy(S);

    K->P0 = gretl_matrix_copy(P);
    K->P1 = gretl_matrix_copy(P);

    /* forecast error vector, per observation */
    K->E = gretl_matrix_alloc(K->n, 1);

    /* just use const pointers for const matrices, don't copy */
    K->F = F;
    K->A = A;
    K->H = H;
    K->Q = Q;
    K->R = R;

    /* will hold P*H */
    K->PH = gretl_matrix_alloc(K->r, K->n);

    /* will hold F*P*H */
    K->FPH = gretl_matrix_alloc(K->r, K->n);

    /* will hold H'*P*H */
    K->HPH = gretl_matrix_alloc(K->n, K->n);

    /* will hold (H'*P*H + R)^{-1} */
    K->V = gretl_matrix_alloc(K->n, K->n);

    /* will hold (H'*P*H + R)^{-1} * E */
    K->VE = gretl_matrix_alloc(K->n, 1);

    /* will hold P*H*V */
    K->PHV = gretl_matrix_alloc(K->r, K->n);

    /* will hold A'*x at obs t */
    K->Ax = gretl_matrix_alloc(K->n, 1);

    /* will hold various intermediate products */
    K->Tmprn = gretl_matrix_alloc(K->r, K->n);
    K->Tmpnn = gretl_matrix_alloc(K->n, K->n);
    K->Tmprr = gretl_matrix_alloc(K->r, K->r);

    if (K->S0 == NULL || K->S1 == NULL || K->P0 == NULL || K->P1 == NULL ||
	K->E == NULL || K->F == NULL || K->A == NULL ||
	K->H == NULL || K->Q == NULL || 
	K->PH == NULL || K->HPH == NULL ||
	K->FPH == NULL || K->V == NULL || K->VE == NULL ||
	K->PHV == NULL || K->Ax == NULL ||
	K->Tmprn == NULL || K->Tmpnn == NULL ||
	K->Tmprr == NULL) {
	*err = E_ALLOC;
	kalman_free(K);
	K = NULL;
    } else {
	gretl_matrix_zero(K->E);
    }

    return K;
}

/* Hamilton (1994) equation [13.2.23] page 381, in simplified notation:

   S+ = FS + FPH(H'PH + R)^{-1} * (y - A'x - H'S) 

   "S" is Hamilton's \xi (state vector)
*/

static int kalman_iter_1 (kalman *K, double *llt)
{
    int err = 0;

    /* write F*S into S+ */
    err += gretl_matrix_multiply(K->F, K->S0, K->S1);

    /* form E = y - A'x - H'S */
    err += gretl_matrix_subtract_from(K->E, K->Ax);
    err += gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
				     K->S0, GRETL_MOD_NONE,
				     K->Tmpnn);
    err += gretl_matrix_subtract_from(K->E, K->Tmpnn);

    /* form (H'PH + R)^{-1} * (y - Ax - H'S) = "VE" */
    err += gretl_matrix_multiply(K->V, K->E, K->VE);

    /* form (y - Ax - H'S)' * (H'PH + R)^{-1} * (y - Ax - H'S) */
    err += gretl_matrix_multiply_mod(K->E, GRETL_MOD_TRANSPOSE,
				     K->VE, GRETL_MOD_NONE,
				     K->Tmpnn);

    /* contribution to log-likelihood of the above -- see Hamilton
       (1994) equation [13.4.1] page 385.
    */
    *llt -= .5 * gretl_matrix_get(K->Tmpnn, 0, 0);

    /* form FPH */
    err += gretl_matrix_multiply(K->F, K->PH, K->FPH);

    /* form FPH * (H'PH + R)^{-1} * (y - A'x - H'S) */
    err += gretl_matrix_multiply(K->FPH, K->VE, K->Tmprn);

    /* complete calculation of S+ */
    err += gretl_matrix_add_to(K->S1, K->Tmprn);

    return err;
}

/* Hamilton (1994) equation [13.2.22] page 380, in simplified notation:

   P+ = F[P - PH(H'PH + R)^{-1}H'P]F' + Q 
*/

static int kalman_iter_2 (kalman *K)
{
    gretl_matrix *PHV = K->Tmprn;
    gretl_matrix *HP;
    int err = 0;

    /* form P - PH(H'PH + R)^{-1}H'P */
    err += gretl_matrix_multiply(K->PH, K->V, PHV);
    HP = gretl_matrix_reuse(K->PH, K->n, K->r);
    err += gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
				     K->P0, GRETL_MOD_NONE,
				     HP);
    err += gretl_matrix_multiply(PHV, HP, K->Tmprr);
    gretl_matrix_subtract_from(K->P0, K->Tmprr);

    /* pre-multiply by F, post-multiply by F' */
    err += gretl_matrix_multiply(K->F, K->P0, K->Tmprr);
    err += gretl_matrix_multiply_mod(K->Tmprr, GRETL_MOD_NONE,
				     K->F, GRETL_MOD_TRANSPOSE,
				     K->P1);
    
    /* add Q */
    err += gretl_matrix_add_to(K->P1, K->Q);

    /* put K->PH back the way we found it */
    gretl_matrix_reuse(K->PH, K->r, K->n);
    
    return err;
}

#if KDEBUG
static void kalman_print_state (kalman *K, int i, double y)
{
    int j;

    fprintf(stderr, "Iteration %d:\n", i);

    for (j=0; j<K->n; j++) {
	fprintf(stderr, "y[%d] = %.8g, err[%d] = %.8g\n", j, 
		gretl_matrix_get(y, t, j), 
		j, gretl_vector_get(K->E, j));
    }

    gretl_matrix_print(K->S0, "K->S0");
    gretl_matrix_print(K->P0, "K->P0");
}
#endif

/* read from the appropriate row of x (T x k) and multiply by
   A' to form A'x_t.
*/

static void kalman_set_Ax (kalman *K, const gretl_matrix *x, int t)
{
    double aji, xj, axi;
    int i, j;

    for (i=0; i<K->n; i++) {
	axi = 0.0;
	for (j=0; j<K->k; j++) {
	    aji = gretl_matrix_get(K->A, j, i);
	    xj = gretl_matrix_get(x, t, j);
	    axi += aji * xj;
	}
	gretl_vector_set(K->Ax, i, axi);
    }

}

/* read from the appropriate row of y (T x n) and transcribe to
   the current E (n x 1)
*/

static void
kalman_initialize_error (kalman *K, const gretl_matrix *y, int t)
{
    double yti;
    int i;

    for (i=0; i<K->n; i++) {
	yti = gretl_matrix_get(y, t, i);
	gretl_vector_set(K->E, i, yti);    
    }
}

/* read the current forecast error and write it into the appropriate
   row of the recorder matrix, E
*/

static void
kalman_record_error (gretl_matrix *E, kalman *K, int t)
{
    double eti;
    int i;

    for (i=0; i<K->n; i++) {
	eti = gretl_vector_get(K->E, i);
	gretl_matrix_set(E, t, i, eti);    
    }
}

/**
 * kalman_forecast:
 * @K: pointer to Kalman struct: see kalman_new().
 * @y: T x n matrix of dependent variable(s).
 * @x: T x k matrix of exogenous variable(s).  May be %NULL if there
 * are no exogenous variables, or if there's only a constant.
 * @E: T x n matrix to hold one-step ahead forecast errors (or %NULL
 * if these do not have to be recorded).
 *
 * Generates a series of one-step ahead forecasts for @y, based on
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

int kalman_forecast (kalman *K, const gretl_matrix *y, const gretl_matrix *x,
		     gretl_matrix *E)
{
    int T = gretl_matrix_rows(y);
    double ldet, llt = 0.0;
    int t, err = 0;

#if KDEBUG
    fprintf(stderr, "kalman_forecast: T = %d\n", T);
#endif  

    K->loglik = 0.0;

    if (x == NULL) {
	/* no exogenous vars */
	gretl_matrix_copy_values(K->Ax, K->A);
    } else if (gretl_matrix_rows(x) != T ||
	       gretl_matrix_cols(x) != K->k) {
	return E_NONCONF;
    }

    if (E != NULL) {
	if (gretl_matrix_rows(E) != T || 
	    gretl_matrix_cols(E) != K->n) {
	    return E_NONCONF;
	}
    }

    for (t=0; t<T && !err; t++) {
#if KDEBUG
	kalman_print_state(K, y, t);
#endif
	/* intial matrix calculations */
	gretl_matrix_multiply(K->P0, K->H, K->PH);
	gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
				  K->PH, GRETL_MOD_NONE,
				  K->HPH);
	if (K->R != NULL) {
	    gretl_matrix_add_to(K->HPH, K->R);
	}
	gretl_matrix_copy_values(K->V, K->HPH);

	ldet = gretl_matrix_log_determinant(K->HPH, &err);
	if (err) {
	    K->loglik = llt = NADBL;
	    break;
	} else {
	    llt = -(K->n / 2.0) * LN_2_PI - .5 * ldet;
	}

	err = gretl_invert_symmetric_matrix(K->V);
	if (err) {
	    break;
	}

	/* read slice from y */
	kalman_initialize_error(K, y, t);
	/* and from x if applicable */
	if (x != NULL) {
	    kalman_set_Ax(K, x, t);
	}

	/* first stage of dual interation */
	err = kalman_iter_1(K, &llt);
	if (na(llt)) {
	    K->loglik = NADBL;
	    err = 1;
	} else {
	    K->loglik += llt;
	    if (isnan(K->loglik)) {
		err = 1;
	    }
	}

	/* record forecast errors if wanted */
	if (E != NULL) {
	    kalman_record_error(E, K, t);
	}

	if (!err) {
	    /* second stage of dual interation */
	    err = kalman_iter_2(K);
	}

	if (!err) {
	    /* update state vector and MSE matrix */
	    gretl_matrix_copy_values(K->S0, K->S1);
	    gretl_matrix_copy_values(K->P0, K->P1);
	}
    }

    if (isnan(K->loglik)) {
	K->loglik = NADBL;
    }    

#if KDEBUG
    fprintf(stderr, "kalman_forecast: err = %d, ll = %.10g\n", err, 
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
