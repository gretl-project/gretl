/*
 *  Copyright (c) by Allin Cottrell and Riccardo "Jack" Lucchetti
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "gretl_matrix.h"
#include "var.h"
#include "jprivate.h"

/* Computations concerned with testing restrictions on the alpha
   (adjustments) matrix in a VECM.  At present only homogeneous,
   common restrictions are handled.  Notation based on S. Johansen,
   "Likelihood-Based Inference in Cointegrated Vector Autoregressive
   Models" (Oxford, 1995).  See section 8.2.1, "The same restriction
   on all \alpha", p. 124ff.
*/

#define ADEBUG 1

/* \bar{A} = A(A'A)^{-1}, where A is the orthogonal complement of the
   restriction R imposed on \alpha in the form R\alpha = 0
*/

static gretl_matrix *make_A_bar (const gretl_matrix *R,
				 gretl_matrix **pA,
				 int *err)
{
    gretl_matrix *A = NULL;
    gretl_matrix *Abar = NULL;
    gretl_matrix *Tmp = NULL;

    A = gretl_matrix_right_nullspace(R, err);
    if (*err) {
	return NULL;
    }

    Tmp = gretl_matrix_alloc(A->cols, A->cols);
    if (Tmp == NULL) {
	*err = E_ALLOC;
	gretl_matrix_free(A);
	return NULL;
    }

    gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
			      A, GRETL_MOD_NONE,
			      Tmp, GRETL_MOD_NONE);

    *err = gretl_invert_symmetric_matrix(Tmp);

    if (!*err) {
	Abar = gretl_matrix_multiply_new(A, Tmp, err);
    }

#if ADEBUG
    gretl_matrix_print(A, "A");
    gretl_matrix_print(Tmp, "Tmp");
    gretl_matrix_print(Abar, "Abar");
#endif

    if (!*err && pA != NULL) {
	*pA = A;
    } else {
	gretl_matrix_free(A);
    }

    gretl_matrix_free(Tmp);

    return Abar;
}

/* "AS00" is the counterpart to S_{00}^{-1} in the unconstrained
   case.  It is given by

   \bar{A} (\bar{A}'S_{00.A_{\perp}}\bar{A})^{-1} \bar{A}'

*/

static gretl_matrix *make_AS00 (const gretl_matrix *R,
				const gretl_matrix *S00a,
				gretl_matrix **pA,
				int *err)
{
    gretl_matrix *AS00 = NULL;
    gretl_matrix *Abar = NULL;
    gretl_matrix *Tmp = NULL;

    Abar = make_A_bar(R, pA, err);
    if (*err) {
	return NULL;
    }

    Tmp = gretl_matrix_alloc(Abar->cols, Abar->cols);
    if (Tmp == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    *err = gretl_matrix_qform(Abar, GRETL_MOD_TRANSPOSE,
			      S00a, Tmp, GRETL_MOD_NONE);

#if ADEBUG
    gretl_matrix_print(Tmp, "Abar'*S00a*Abar");
#endif

    if (!*err) {
	*err = gretl_invert_symmetric_matrix(Tmp);
    }

    if (!*err) {
	AS00 = gretl_matrix_alloc(Abar->rows, Abar->rows);
	if (AS00 == NULL) {
	    *err = E_ALLOC;
	    goto bailout;
	}
    }

    if (!*err) {
	*err = gretl_matrix_qform(Abar, GRETL_MOD_NONE, Tmp,
				  AS00, GRETL_MOD_NONE);
    }

#if ADEBUG
    gretl_matrix_print(AS00, "AS00");
#endif

 bailout:

    gretl_matrix_free(Tmp);
    gretl_matrix_free(Abar);

    return AS00;
}

/* Solve the eigenvalue problem

   |\tilde{\lambda}S_{11.A_{\perp}} - S_{10.A_{\perp}}
      \bar{A}(\bar{A}' S_{00.A_{\perp}} \bar{A})^{-1}
        \bar{A}'S_{01.A_{\perp}}| = 0

   (Johansen, p. 126)
*/

int alt_get_eigenvalues (gretl_matrix *AS00,
			 const gretl_matrix *S01a,
			 const gretl_matrix *S11a,
			 gretl_matrix *M,
			 gretl_matrix **evals,
			 gretl_matrix *Tmp,
			 int rank)
{
    int err = 0;

    gretl_matrix_qform(S01a, GRETL_MOD_TRANSPOSE, 
		       AS00, Tmp, GRETL_MOD_NONE);

    *evals = gretl_gensymm_eigenvals(Tmp, S11a, M, &err);

    if (!err) {
	err = gretl_symmetric_eigen_sort(*evals, M, rank);
    }

    return err;
}

static int alpha_compute_alpha (JohansenInfo *jv,
				const gretl_matrix *S11a,
				const gretl_matrix *S01a)
{
    const gretl_matrix *B = jv->Beta;
    gretl_matrix *a = NULL;
    gretl_matrix *BSB = NULL;
    gretl_matrix *Tmp = NULL;
    int err = 0;

    BSB = gretl_matrix_alloc(B->cols, B->cols);
    Tmp = gretl_matrix_alloc(B->rows, B->cols);
    a = gretl_matrix_alloc(jv->S01->rows, B->cols);

    if (BSB == NULL || Tmp == NULL || a == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	err = gretl_matrix_qform(B, GRETL_MOD_TRANSPOSE, S11a,
				 BSB, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(BSB);
    }

    if (!err) {
	gretl_matrix_multiply(B, BSB, Tmp);
	gretl_matrix_multiply(S01a, Tmp, a);
    }

    gretl_matrix_free(BSB);
    gretl_matrix_free(Tmp);

    if (!err) {
	jv->Alpha = a;
    } else {
	gretl_matrix_free(a);
    }

    return err;
}

static int 
alpha_test_show_beta (JohansenInfo *jv, 
		      const gretl_matrix *M,
		      const gretl_matrix *S11a,
		      const gretl_matrix *S01a,
		      const DATAINFO *pdinfo,
		      PRN *prn)
{
    int err = 0;

    gretl_matrix_copy_values(jv->Beta, M);

    if (jv->rank) { 
	/* and if r > 1? */
	double den = jv->Beta->val[0];

	if (!floateq(den, 0.0)) {
	    gretl_matrix_divide_by_scalar(jv->Beta, den);
	}
    }

    if (!err) {
	err = alpha_compute_alpha(jv, S11a, S01a);
    }

    if (!err) {
	print_beta_alpha_Pi(jv, pdinfo, prn);
    }

    return err;
}

int vecm_alpha_test (GRETL_VAR *jvar, 
		     const gretl_restriction_set *rset,
		     const DATAINFO *pdinfo, 
		     gretlopt opt,
		     PRN *prn)
{
    const gretl_matrix *R = rset_get_R_matrix(rset);
    const gretl_matrix *S00 = jvar->jinfo->S00;
    const gretl_matrix *S01 = jvar->jinfo->S01;
    const gretl_matrix *S11 = jvar->jinfo->S11;

    gretl_matrix *ASA = NULL;
    gretl_matrix *C = NULL;
    gretl_matrix *Tmp = NULL;
    gretl_matrix *S00a = NULL;
    gretl_matrix *S11a = NULL;
    gretl_matrix *S01a = NULL;

    int rank = jvar->jinfo->rank;
    int n = jvar->neqns;
    int m = S11->rows;
    int err = 0;

    if (!simple_restriction(jvar, rset)) {
	err = E_NOTIMP;
    }

    ASA = gretl_matrix_alloc(R->rows, R->rows);
    C = gretl_matrix_alloc(n, n);
    Tmp = gretl_matrix_alloc(m, m);
    S00a = gretl_zero_matrix_new(n, n);
    S11a = gretl_zero_matrix_new(m, m);
    S01a = gretl_zero_matrix_new(n, m);

    if (ASA == NULL || C == NULL || Tmp == NULL || 
	S00a == NULL || S11a == NULL || S01a == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    gretl_matrix_qform(R, GRETL_MOD_NONE, S00,
		       ASA, GRETL_MOD_NONE);

    gretl_invert_symmetric_matrix(ASA);

    gretl_matrix_qform(R, GRETL_MOD_TRANSPOSE, ASA,
		       C, GRETL_MOD_NONE);

    /* Johansen, p. 124: compute the A_{\perp}-transformed 
       moment matrices */

    /* S_{00.A_{\perp}} */
    gretl_matrix_qform(S00, GRETL_MOD_TRANSPOSE, C, 
		       S00a, GRETL_MOD_NONE);
    gretl_matrix_subtract_reversed(S00, S00a);

    /* S_{11.A_{\perp}} */
    gretl_matrix_qform(S01, GRETL_MOD_TRANSPOSE, C, 
		       S11a, GRETL_MOD_NONE);
    gretl_matrix_subtract_reversed(S11, S11a);

    /* S_{01.A_{\perp}} */
    gretl_matrix_reuse(Tmp, n, n);
    gretl_matrix_multiply(S00, C, Tmp);
    gretl_matrix_multiply(Tmp, S01, S01a);
    gretl_matrix_subtract_reversed(S01, S01a);

#if ADEBUG
    gretl_matrix_print(S00a, "S00a");
    gretl_matrix_print(S11a, "S11a");
    gretl_matrix_print(S01a, "S01a");
#endif

    if (!err) {
	/* do the eigenvalue thingy */
	gretl_matrix *A = NULL;
	gretl_matrix *AS00 = NULL;
	gretl_matrix *M = NULL;
	gretl_matrix *evals = NULL;

	AS00 = make_AS00(R, S00a, &A, &err);

	if (!err) {
	    M = gretl_matrix_alloc(m, m);
	    if (M == NULL) {
		err = E_ALLOC;
	    }
	}

	if (!err) {
	    gretl_matrix_reuse(Tmp, m, m);
	    err = alt_get_eigenvalues(AS00, S01a, S11a,
				      M, &evals, Tmp, 
				      rank);
	}
    
	if (!err) {
	    /* is the df right, in general? */
	    johansen_LR_calc(jvar, evals, A, V_ALPHA, prn);
	} 

	if (!err && (opt & OPT_V)) {
	    alpha_test_show_beta(jvar->jinfo, M, S11a, S01a,
				 pdinfo, prn);
	}

	gretl_matrix_free(evals);
	gretl_matrix_free(M);
	gretl_matrix_free(AS00);
	gretl_matrix_free(A);
    }

 bailout:

    gretl_matrix_free(ASA);
    gretl_matrix_free(C);
    gretl_matrix_free(Tmp);
    gretl_matrix_free(S00a);
    gretl_matrix_free(S11a);
    gretl_matrix_free(S01a);
  
    return err;
}
