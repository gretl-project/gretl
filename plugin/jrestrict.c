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

/* General restrictions on beta in context of a VECM.  Based on an ox
   program by Jack, June 2007.  Rendered into C by Allin.
*/

#include "libgretl.h"
#include "random.h"
#include "johansen.h"
#include "var.h"
#include "gretl_restrict.h"
#include "jprivate.h"

#define JDEBUG 0

typedef struct Jwrap_ Jwrap;

struct Jwrap_ {
    int T;          /* length of time series used */
    int neqns;      
    int rank;
    int blen;       /* number of unrestricted coefficients in beta */
    int nC;         /* number of restriction matrices */
    int df;         /* degrees if freedom for LR test */
    double ldS00;   /* base component of log-likelihood */
    double ll;      /* log-likelihood */

    /* moment matrices and copies */
    const gretl_matrix *S00;
    const gretl_matrix *S01;
    const gretl_matrix *S11;

    gretl_matrix *S00i;
    gretl_matrix *S11m;

    /* restrictions on beta */
    gretl_matrix *H;
    gretl_matrix *s;
    gretl_matrix **h;

    /* coefficients and variances */
    gretl_matrix *beta;
    gretl_matrix *alpha;
    gretl_matrix *Omega;
    gretl_matrix *V;
    gretl_matrix *se;

    /* temp storage for beta calculation */
    gretl_matrix *phivec;
    gretl_matrix *bcol;

    /* temp storage for likelihood calculation */
    gretl_matrix *qf1;
    gretl_matrix *qf2;
};

static Jwrap *jwrap_new (int neqns, int rank, int T)
{
    Jwrap *J = malloc(sizeof *J);

    if (J == NULL) {
	return NULL;
    }

    J->T = T;
    J->neqns = neqns;
    J->rank = rank;
    J->blen = 0;
    J->nC = 0;
    J->df = 0;

    J->ll = NADBL;

    J->S00 = NULL;
    J->S01 = NULL;
    J->S11 = NULL;

    J->S00i = NULL;
    J->S11m = NULL;

    J->H = NULL;
    J->s = NULL;
    J->h = NULL;

    J->beta = NULL;
    J->alpha = NULL;
    J->Omega = NULL;
    J->V = NULL;
    J->se = NULL;

    J->phivec = NULL;
    J->bcol = NULL;

    J->qf1 = NULL;
    J->qf2 = NULL;

    return J;
}

static void jwrap_destroy (Jwrap *J)
{
    gretl_matrix_free(J->S00i);
    gretl_matrix_free(J->S11m);

    gretl_matrix_free(J->H);
    gretl_matrix_free(J->s);

    if (J->h != NULL) {
	gretl_matrix_array_free(J->h, J->nC);
    }

    gretl_matrix_free(J->beta);
    gretl_matrix_free(J->alpha);
    gretl_matrix_free(J->Omega);
    gretl_matrix_free(J->V);
    gretl_matrix_free(J->se);

    gretl_matrix_free(J->phivec);
    gretl_matrix_free(J->bcol);

    gretl_matrix_free(J->qf1);
    gretl_matrix_free(J->qf2);

    free(J);
}

static int make_S_matrices (Jwrap *J, const GRETL_VAR *jvar)
{
    gretl_matrix *Tmp = NULL;
    int err = 0;

    J->S00 = jvar->jinfo->Suu;
    J->S01 = jvar->jinfo->Suv;
    J->S11 = jvar->jinfo->Svv;

    J->S00i = gretl_matrix_copy(J->S00);
    J->S11m = gretl_matrix_copy(J->S11);

    if (J->S00i == NULL || J->S11m == NULL) {
	return E_ALLOC;
    }

    Tmp = gretl_matrix_alloc(J->S01->cols, J->S01->cols);

    if (Tmp == NULL) {
	return E_ALLOC;
    }

    J->ldS00 = gretl_matrix_log_determinant(J->S00i, &err);

    if (!err) {
	gretl_matrix_copy_values(J->S00i, J->S00);
	err = gretl_invert_symmetric_matrix(J->S00i);
    }

    if (!err) {
	err = gretl_matrix_qform(J->S01, GRETL_MOD_TRANSPOSE,
				 J->S00i, Tmp, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_matrix_subtract_from(J->S11m, Tmp);
    }

#if JDEBUG > 1
    gretl_matrix_print(J->S00, "S00");
    gretl_matrix_print(J->S01, "S01");
    gretl_matrix_print(J->S11, "S11");
    gretl_matrix_print(J->S00i, "S00i");
    gretl_matrix_print(J->S11m, "S11m");
#endif    

    gretl_matrix_free(Tmp);

    return err;
}

/* Convert the implicit constraint matrices R and q into explicit form.
   The matrices put in the locations *ph and *ps are newly allocated
   in this function.
*/

static int imp2exp (const gretl_matrix *Ri, const gretl_matrix *qi,
		    gretl_matrix **ph, gretl_matrix **ps)
{
    gretl_matrix *RRT = NULL;
    gretl_matrix *Tmp = NULL;
    int m = Ri->rows;
    int n = Ri->cols;
    int err = 0;

    int rr = gretl_matrix_rank(Ri, &err);

    *ph = NULL;
    *ps = NULL;

    if (rr<n) {
	/* 
	   standard case: \beta_i only partially constrained 
	*/

    RRT = gretl_matrix_alloc(m, m);
    Tmp = gretl_matrix_alloc(n, m);

    if (RRT == NULL || Tmp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

#if JDEBUG 
    gretl_matrix_print(Ri, "Ri");
    gretl_matrix_print(qi, "qi");
#endif    

    *ph = gretl_matrix_right_nullspace(Ri, &err);

#if JDEBUG 
    fprintf(stderr, "right_nullspace: err = %d\n", err);
    gretl_matrix_print(*ph, "*ph");
#endif

    if (!err) {
	err = gretl_matrix_multiply_mod(Ri, GRETL_MOD_NONE,
					Ri, GRETL_MOD_TRANSPOSE,
					RRT, GRETL_MOD_NONE);
    }
    
    if (!err) {
	err = gretl_invert_symmetric_matrix(RRT);
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(Ri, GRETL_MOD_TRANSPOSE,
					RRT, GRETL_MOD_NONE,
					Tmp, GRETL_MOD_NONE);
    }

    if (!err) {
	*ps = gretl_matrix_alloc(Tmp->rows, qi->cols);
	if (*ps == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = gretl_matrix_multiply(Tmp, qi, *ps);
    }
    } else {
	/* 
	   special case: \beta_i is completely specified 
	   by the restrictions
	*/

	*ps = gretl_matrix_copy(qi);
	err = gretl_LU_solve(Ri, *ps);
    }

 bailout:

    gretl_matrix_free(RRT);
    gretl_matrix_free(Tmp);

    return err;
}

/* Determine the number of rows in a given diagonal block of the "big
   R" restriction matrix: we count from the given starting row till we
   hit a row that is all zeros, up to the specified column.
*/

static int 
get_Ri_rows (const gretl_matrix *R, int nb, int i, int *start)
{
    int stop, nz = nb * (i + 1);
    int j, r = 0;

    for (i=*start; i<R->rows; i++) {
	stop = 1;
	for (j=0; j<nz; j++) {
	    if (gretl_matrix_get(R, i, j) != 0) {
		stop = 0;
		r++;
		break;
	    }
	}
	if (stop) {
	    break;
	}
    }

    *start += r;

    return r;
}

/* Here we construct both the big block-diagonal restrictions matrix,
   J->H, and an array holding the block-diagonal elements, J->h,
   which will be used later in computing the variance of beta.
*/

static int set_up_restrictions (Jwrap *J, GRETL_VAR *jvar,
				const gretl_restriction_set *rset)
{
    const gretl_matrix *R;
    const gretl_matrix *q;

    gretl_matrix *Ri;
    gretl_matrix *qi;
    gretl_matrix **ss;

    int nb = gretl_VECM_n_beta(jvar);
    int nC, hr = 0, hc = 0, sr = 0;
    int Rr = 0, Rc = 0, irmax = 0;
    int r0, ir;
    int i, err = 0;

    if (nb <= 0) {
	return E_DATA;
    }	

    R = rset_get_R_matrix(rset);
    q = rset_get_q_matrix(rset);

    nC = R->cols / nb;

#if JDEBUG
    gretl_matrix_print(R, "R, in set_up_restrictions");
    gretl_matrix_print(q, "q, in set_up_restrictions");
    fprintf(stderr, "nC = %d\n", nC);
#endif

    J->h = gretl_matrix_array_alloc(nC);
    if (J->h == NULL) {
	return E_ALLOC;
    }

    ss = gretl_matrix_array_alloc(nC);
    if (ss == NULL) {
	return E_ALLOC;
    }

    /* find max rows per restriction sub-matrix */
    r0 = 0;
    for (i=0; i<nC; i++) {
	ir = get_Ri_rows(R, nb, i, &r0);
	if (ir > irmax) {
	    irmax = ir;
	}
    }

    Ri = gretl_matrix_alloc(irmax, nb);
    qi = gretl_matrix_alloc(irmax, 1);

    J->nC = r0 = 0;

    for (i=0; i<nC && !err; i++) {
	ir = get_Ri_rows(R, nb, i, &r0);
	gretl_matrix_reuse(Ri, ir, 0);
	gretl_matrix_reuse(qi, ir, 0);
	gretl_matrix_extract_matrix(Ri, R, Rr, Rc, GRETL_MOD_NONE);
	gretl_matrix_extract_matrix(qi, q, Rr,  0, GRETL_MOD_NONE);
	err = imp2exp(Ri, qi, &J->h[i], &ss[i]);
#if JDEBUG
	gretl_matrix_print(J->h[i], "H[i], in set_up_restrictions");
	gretl_matrix_print(ss[i], "ss[i], in set_up_restrictions");
#endif
	if (!err) {
	    if (J->h[i] == NULL) {
		/* may happen if the i-th contegration vector
		   is fully restricted */
		hr += Ri->cols;
	    } else {
		hr += J->h[i]->rows;
		hc += J->h[i]->cols;
	    }
	    sr += ss[i]->rows;
	    Rr += ir;
	    Rc += nb;
	    J->nC += 1;
	}
    }

    gretl_matrix_free(Ri);
    gretl_matrix_free(qi);

#if JDEBUG
    fprintf(stderr, "H: hr = %d, hc = %d (sr = %d)\n", hr, hc, sr);
#endif

    if (!err) {
	J->H = gretl_zero_matrix_new(hr, hc);
	J->s = gretl_column_vector_alloc(sr);

	if (J->H == NULL || J->s == NULL) {
	    err = E_ALLOC;
	} else {
	    hr = hc = sr = 0;
	}
    }

    for (i=0; i<nC && !err; i++) {
	err = gretl_matrix_inscribe_matrix(J->s, ss[i], sr, 0, 
					   GRETL_MOD_NONE);
	if (!err) {
	    sr += ss[i]->rows;
	}
	if (!err && (J->h[i] != NULL)) {
	    err = gretl_matrix_inscribe_matrix(J->H, J->h[i], hr, hc, 
					       GRETL_MOD_NONE);
	    if (!err) {
		hr += J->h[i]->rows;
		hc += J->h[i]->cols;
	    }
	}
    }

    gretl_matrix_array_free(ss, nC);

    return err;
}

static int 
normalize_beta (Jwrap *J, const gretl_restriction_set *rset, 
		gretl_matrix *b)
{
    const gretl_matrix *R = rset_get_R_matrix(rset);
    const gretl_matrix *d = rset_get_q_matrix(rset);

    gretl_matrix *tmp_sq = NULL;
    gretl_matrix *tmp_b = NULL;
    gretl_matrix *tmp = NULL;
    gretl_matrix *X = NULL;

    double x;
    int i, j, ii;
    int br = b->rows;
    int bc = J->rank;
    int err = 0;

    tmp_sq = gretl_identity_matrix_new(bc);
    tmp_b = gretl_matrix_alloc(br, bc);
    if (tmp_sq == NULL || tmp_b == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    tmp = gretl_matrix_kronecker_product_new(tmp_sq, b);
    if (tmp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    X = gretl_matrix_alloc(R->rows, tmp->cols);
    if (X == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
	
    gretl_matrix_multiply(R, tmp, X);
    gretl_matrix_free(tmp);

    tmp = gretl_matrix_alloc(bc*bc, 1);
    if (tmp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    gretl_matrix_multi_ols(d, X, tmp, NULL);

    ii = 0;
    for (i=0; i<bc; i++) {
	for (j=0; j<bc; j++) {
	    x = gretl_matrix_get(tmp, ii++, 0);
	    gretl_matrix_set(tmp_sq, j, i, x);
	}
    }
	    
    gretl_matrix_copy_values(tmp_b, b);
    gretl_matrix_multiply(tmp_b, tmp_sq, b);

 bailout:

    gretl_matrix_free(tmp);
    gretl_matrix_free(tmp_b);
    gretl_matrix_free(tmp_sq);
    gretl_matrix_free(X);

    return err;
}

static int case0 (Jwrap *J)
{
    gretl_matrix *M = NULL;
    gretl_matrix *Tmp = NULL;
    gretl_matrix *evals = NULL;
    int n = J->S11->cols;
    int err = 0;

    Tmp = gretl_matrix_alloc(n, n);
    M = gretl_matrix_alloc(n, n);

    if (Tmp == NULL || M == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_matrix_qform(J->S01, GRETL_MOD_TRANSPOSE, 
				 J->S00i, Tmp, GRETL_MOD_NONE);
    }

    if (!err) {
	evals = gretl_gensymm_eigenvals(Tmp, J->S11, M, &err);
    }

    if (!err) {
	err = gretl_symmetric_eigen_sort(evals, M, J->rank);
    }

#if JDEBUG
    gretl_matrix_print(evals, "case0: evals");
    gretl_matrix_print(M, "case0: M");
    fprintf(stderr, "(err = %d)\n", err);
#endif

    if (!err) {
	J->beta = M;
    } else {
	gretl_matrix_free(M);
    }

    gretl_matrix_free(Tmp);
    gretl_matrix_free(evals);

    return err;
}

static int initval (Jwrap *J, const gretl_restriction_set *rset,
		    gretl_matrix **pb)
{
    gretl_matrix *HHi = NULL;
    gretl_matrix *vecb = NULL;
    gretl_matrix *tmp = NULL;
    int Hcols;
    int err;

    err = case0(J);
    if (err) {
	return err;
    }

    err = normalize_beta(J, rset, J->beta);
    if (err) {
	return err;
    }

    Hcols = J->H->cols;

    HHi = gretl_matrix_alloc(Hcols, Hcols);
    vecb = gretl_column_vector_alloc(J->H->rows);
    tmp = gretl_column_vector_alloc(Hcols);

    if (HHi == NULL || vecb == NULL || tmp == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(J->H, GRETL_MOD_TRANSPOSE,
					J->H, GRETL_MOD_NONE,
					HHi, GRETL_MOD_NONE);
    }
    
    if (!err) {
	err = gretl_invert_symmetric_matrix(HHi);
    }

    if (!err) {
	err = gretl_matrix_vectorize(vecb, J->beta);
    }

    if (!err) {
	err = gretl_matrix_subtract_from(vecb, J->s);
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(J->H, GRETL_MOD_TRANSPOSE,
					vecb, GRETL_MOD_NONE,
					tmp, GRETL_MOD_NONE);
    }

    if (!err) {
	gretl_matrix_reuse(vecb, tmp->rows, 1);
	err = gretl_matrix_multiply(HHi, tmp, vecb);
    }

    if (!err) {
#if JDEBUG
	fprintf(stderr, "*** initval: vecb->rows = %d\n", vecb->rows);
#endif
	J->blen = vecb->rows;
	*pb = vecb;
    } else {
	gretl_matrix_free(vecb);
    }

    gretl_matrix_free(HHi);
    gretl_matrix_free(tmp);
    
    return err;
}

static int make_omega (Jwrap *J,
		       const gretl_matrix *alpha, 
		       const gretl_matrix *beta)
{
    gretl_matrix *tmp = NULL;
    int err = 0;

    tmp = gretl_matrix_alloc(beta->cols, beta->cols);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    J->Omega = gretl_matrix_copy(J->S00);
    if (J->Omega == NULL) {
	gretl_matrix_free(tmp);
	return E_ALLOC;
    }

    err = gretl_matrix_qform(beta, GRETL_MOD_TRANSPOSE,
			     J->S11, tmp, GRETL_MOD_NONE);

    if (!err) {
	gretl_matrix_qform(alpha, GRETL_MOD_NONE,
			   tmp, J->Omega, GRETL_MOD_DECUMULATE);
    }

    if (!err) {
	gretl_matrix_divide_by_scalar(J->Omega, J->T);
    }

    gretl_matrix_free(tmp);

    return err;
}

static int make_beta_variance (Jwrap *J)
{
    const gretl_matrix *a = J->alpha;
    const gretl_matrix *b = J->beta;

    gretl_matrix *V = NULL;
    gretl_matrix *aiom = NULL;

    int r = b->cols;
    int n = b->rows;
    int npar = 0;
    int istart, jstart;
    int i, j, err = 0;

    for (i=0; i<r; i++) { 
	if (J->h[i] != NULL) {
	npar += J->h[i]->cols;
    }
    }

    V = gretl_zero_matrix_new(npar, npar);
    aiom = gretl_matrix_alloc(a->cols, a->cols);
    if (V == NULL || aiom == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = gretl_matrix_qform(a, GRETL_MOD_TRANSPOSE,
			     J->Omega, aiom, GRETL_MOD_NONE);

#if JDEBUG > 1
    gretl_matrix_print(aiom, "aiom");
#endif

    istart = 0;
    for (i=0; i<r && !err; i++) {
	if (J->h[i] != NULL) {
	gretl_matrix *HiS = gretl_matrix_alloc(J->h[i]->cols, J->S11->cols);

	err = gretl_matrix_multiply_mod(J->h[i], GRETL_MOD_TRANSPOSE,
					J->S11, GRETL_MOD_NONE,
					HiS, GRETL_MOD_NONE);

	jstart = 0;
	for (j=0; j<r && !err; j++) {
		if (J->h[j] != NULL) {
	    gretl_matrix *Vij = gretl_matrix_alloc(HiS->rows, J->h[j]->cols);
	    double rij = gretl_matrix_get(aiom, i, j);

	    gretl_matrix_multiply(HiS, J->h[j], Vij);
	    gretl_matrix_multiply_by_scalar(Vij, rij);
#if JDEBUG > 1
	    gretl_matrix_print(Vij, "Vij");
	    fprintf(stderr, "inscribing Vij at %d, %d\n", istart, jstart);
#endif
	    err = gretl_matrix_inscribe_matrix(V, Vij, istart, jstart, GRETL_MOD_NONE);

	    jstart += J->h[j]->cols;
	    gretl_matrix_free(Vij);
	}
	    }

	istart += J->h[i]->cols;
	gretl_matrix_free(HiS);
    }
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(V);
    }

    if (!err) {
	J->V = gretl_matrix_alloc(J->H->rows, J->H->rows);
	if (J->V == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gretl_matrix_qform(J->H, GRETL_MOD_NONE,
				     V, J->V, GRETL_MOD_NONE);
	}
    }

    if (!err) {
	double x;

	J->se = gretl_matrix_alloc(n, r);
	if (J->se == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<J->V->rows; i++) {
		x = gretl_matrix_get(J->V, i, i);
		J->se->val[i] = sqrt(x);
	    }
	}
    }

 bailout:

    gretl_matrix_free(aiom);
    gretl_matrix_free(V);

    return err;
}

static int make_beta (Jwrap *J, const double *phi)
{
    int i, nbeta = J->beta->rows * J->beta->cols;
    int err = 0;

    if (J->blen == 0) {
	fprintf(stderr, "*** make_beta: blen = 0\n");
	return E_DATA;
    }

    if (J->phivec == NULL) {
	/* temp storage not allocated yet */
	J->phivec = gretl_column_vector_alloc(J->blen);
	J->bcol = gretl_column_vector_alloc(nbeta);
	if (J->phivec == NULL || J->bcol == NULL) {
	    return E_ALLOC;
	}
    }

    if (!err) {
	for (i=0; i<J->blen; i++) {
	    J->phivec->val[i] = phi[i];
	}
	err = gretl_matrix_multiply(J->H, J->phivec, J->bcol);
    }

    if (!err) {
	gretl_matrix_add_to(J->bcol, J->s);
    }

    if (!err) {
	for (i=0; i<nbeta; i++) {
	    J->beta->val[i] = J->bcol->val[i];
	}
    }

#if JDEBUG
    gretl_matrix_print(J->phivec, "phi");
    gretl_matrix_print(J->beta, "beta");
#endif

    return err;
}

/* BFGS callback function */

static double Jloglik (const double *phi, void *data)
{
    Jwrap *J = (Jwrap *) data;
    double ret = J->ldS00;
    int bcols = J->beta->cols;
    int err = 0;

    if (J->qf1 == NULL) {
	/* temp storage not allocated yet */
	J->qf1 = gretl_matrix_alloc(bcols, bcols);
	J->qf2 = gretl_matrix_alloc(bcols, bcols);
	if (J->qf1 == NULL || J->qf2 == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = make_beta(J, phi);
    }

    if (!err) {
	gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE,
			   J->S11m, J->qf1, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE,
				 J->S11, J->qf2, GRETL_MOD_NONE);
    }

    if (!err) {
	ret += gretl_matrix_log_determinant(J->qf1, &err);
    }

    if (!err) {
	ret -= gretl_matrix_log_determinant(J->qf2, &err);
    }

    if (err) {
	J->ll = ret = NADBL;
    } else {
	J->ll = ret = -J->T * 0.5 * (J->neqns * (1.0 + LN_2_PI) + ret);
    }

    return ret;
}

/* See Johansen pp. 106-112 */

static void set_LR_df (Jwrap *J, GRETL_VAR *jvar)
{
    int p = gretl_matrix_rows(J->beta);
    int r = J->rank;
    int i, si;

    J->df = 0;
    for (i=0; i<J->nC; i++) {
	if (J->h[i] != NULL) {
	    si = J->h[i]->cols;
	    J->df += p - r - si;
	}
    }

    /* system was subject to a prior restriction */
    J->df -= jvar->jinfo->bdf;
}

#define VECM_WIDTH 13

static int printres (Jwrap *J, GRETL_VAR *jvar, const DATAINFO *pdinfo,
		     PRN *prn)
{
    JohansenInfo *jv = jvar->jinfo;
    const gretl_matrix *b = J->beta;
    const gretl_matrix *sd = J->se;
    char vname[32], s[16];
    int n = b->rows;
    int r = b->cols;
    int i, j;

    pprintf(prn, _("Unrestricted loglikelihood (lu) = %g\n"), jvar->ll);
    pprintf(prn, _("Restricted loglikelihood (lr) = %g\n"), J->ll);
    if (J->df > 0) {
	double x = 2.0 * (jvar->ll - J->ll);

	pprintf(prn, "2 * (lu - lr) = %g\n", x);
	pprintf(prn, _("P(Chi-Square(%d) > %g = %g\n"), J->df, x, 
		chisq_cdf_comp(x, J->df));
    }

    pputs(prn, "\n\n");
    pputs(prn, _("Restricted cointegrating vectors"));
    pprintf(prn, " (%s)", _("standard errors in parentheses"));
    pputs(prn, "\n\n");

    for (i=0; i<n; i++) {
	if (i < jv->list[0]) {
	    sprintf(vname, "%s(-1)", pdinfo->varname[jv->list[i+1]]);
	} else if (jv->code == J_REST_CONST) {
	    strcpy(vname, "const");
	} else if (jv->code == J_REST_TREND) {
	    strcpy(vname, "trend");
	}
	pprintf(prn, "%-12s", vname); /* FIXME */

	for (j=0; j<r; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(b, i, j));
	}
	pputc(prn, '\n');

	bufspace(VECM_WIDTH, prn);

	for (j=0; j<r; j++) {
	    sprintf(s, "(%#.5g)", gretl_matrix_get(sd, i, j));
	    pprintf(prn, "%12s ", s);
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');

    return 0;
}

static int simann (Jwrap *J, gretl_matrix *b)
{
    int i, SAiter = 4096;
    double f0, f1, fbest;
    double rndu;
    int jump;

    gretl_matrix *b0 = NULL;
    gretl_matrix *b1 = NULL;
    gretl_matrix *bbest = NULL;
    gretl_matrix *d = NULL;

    double Temp = 1.0;
    double radius = 1.0;
    int err = 0;

    b0 = gretl_matrix_copy(b);
    b1 = gretl_matrix_copy(b);
    bbest = gretl_matrix_copy(b);
    d = gretl_column_vector_alloc(b->rows);

    if (b0 == NULL || b1 == NULL || bbest == NULL || d == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    f0 = fbest = Jloglik(b->val, J);

    for (i=0; i<SAiter; i++) {
	gretl_matrix_random_fill(d, D_NORMAL);
	gretl_matrix_multiply_by_scalar(d, radius);
	gretl_matrix_add_to(b1, d);
	f1 = Jloglik(b1->val, J);

	if (f1 > f0) {
	    jump = 1;
	} else {
	    rndu = ((double) rand()) / RAND_MAX;
	    jump = (Temp < rndu);
	}

	if (jump) {
	    f0 = f1;
	    gretl_matrix_copy_values(b0, b1);
	    if (f0 > fbest) {
		fbest = f0;
		gretl_matrix_copy_values(bbest, b0);
		fprintf(stderr, "i:%d\tTemp = %#g, radius = %#g, fbest = %#g\n", 
			i, Temp, radius, fbest);
	    }
	} else {
	    gretl_matrix_copy_values(b1, b0);
	    f1 = f0;
	}

	Temp *= 0.999;
	radius *= 0.9999;
    }
    
    gretl_matrix_copy_values(b, bbest);

 bailout:

    gretl_matrix_free(b0);
    gretl_matrix_free(b1);
    gretl_matrix_free(bbest);
    gretl_matrix_free(d);

    return err;
}

/* public entry point */

int general_beta_analysis (GRETL_VAR *jvar, 
			   const gretl_restriction_set *rset,
			   const DATAINFO *pdinfo,
			   gretlopt opt,
			   PRN *prn)
{
    Jwrap *J = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *S01b = NULL;

    int n = jvar->neqns;
    int rank = jrank(jvar);
    int err = 0;

    J = jwrap_new(n, rank, jvar->T);
    if (J == NULL) {
	return E_ALLOC;
    }

    err = make_S_matrices(J, jvar);

    if (!err) {
	err = set_up_restrictions(J, jvar, rset);
    }

    if (!err) {
	err = initval(J, rset, &b);
#if JDEBUG
	fprintf(stderr, "after initval: err = %d\n", err);
	gretl_matrix_print(b, "b, before BFGS");
#endif
    }

#if 1
    err = simann(J, b);
#endif

    if (!err) {
	int maxit = 4000;
	double reltol = 1.0e-11;
	int fncount = 0;
	int grcount = 0;
	int nn = b->rows;

	err = LBFGS_max(b->val, nn, maxit, reltol, 
			&fncount, &grcount, Jloglik, C_LOGLIK,
			NULL, J, opt, prn);
    }

    if (!err) {
	S01b = gretl_matrix_alloc(J->S01->rows, J->rank);
	J->alpha = gretl_matrix_alloc(J->S01->rows, J->rank);
	if (S01b == NULL || J->alpha == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = gretl_matrix_multiply(J->S01, J->beta, S01b);
    }

    if (!err) {
	err = gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE,
				 J->S11, J->qf1, GRETL_MOD_NONE);
    }

    if (!err) {
	gretl_invert_symmetric_matrix(J->qf1);
    }

    if (!err) {
	err = gretl_matrix_multiply(S01b, J->qf1, J->alpha);
    }

    if (!err) {
	err = make_omega(J, J->alpha, J->beta);
    }

    if (opt & OPT_F) {
	gretl_matrix_free(jvar->S);
	jvar->S = gretl_matrix_copy(J->Omega);
    }

    if (!err) {
	gretl_invert_symmetric_matrix(J->Omega);
    }

    if (!err) {
	err = make_beta_variance(J);
    }

    if (!err) {
	set_LR_df(J, jvar);

	if (opt & OPT_F) {
	    jvar->jinfo->ll0 = jvar->ll;
	    jvar->ll = J->ll;
	    jvar->jinfo->bdf += J->df; /* ?? */

	    gretl_matrix_free(jvar->jinfo->Beta);
	    jvar->jinfo->Beta = J->beta;
	    J->beta = NULL;

	    gretl_matrix_free(jvar->jinfo->Alpha);
	    jvar->jinfo->Alpha = J->alpha;
	    J->alpha = NULL;

	    gretl_matrix_free(jvar->jinfo->Bvar);
	    jvar->jinfo->Bvar = J->V;
	    J->V = NULL;

	    gretl_matrix_free(jvar->jinfo->Bse);
	    jvar->jinfo->Bse = J->se;
	    J->se = NULL;
	} else {
	    printres(J, jvar, pdinfo, prn);
	}
    } 

    jwrap_destroy(J);

    gretl_matrix_free(S01b);
    gretl_matrix_free(b);

    return err;
}
