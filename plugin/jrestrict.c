/* General restrictions on beta in context of a VECM.
   Based on Jack's ox program, June 2007
*/

#include "libgretl.h"
#include "johansen.h"
#include "var.h"
#include "gretl_restrict.h"

#define JDEBUG 0

typedef struct Jwrap_ Jwrap;

struct Jwrap_ {
    int T;          /* length of time series used */
    int neqns;      
    int rank;
    int blen;       /* number of unrestricted coefficients in beta */
    int nC;         /* number of restriction matrices */
    double ll;      /* log-likelihood */

    /* moment matrices and copies */
    gretl_matrix *S00;
    gretl_matrix *S01;
    gretl_matrix *S11;
    gretl_matrix *S00i;
    gretl_matrix *S11c;

    /* restrictions on beta */
    gretl_matrix *H;
    gretl_matrix *s;
    gretl_matrix **Hs;

    gretl_matrix *beta;
    gretl_matrix *alpha;
    gretl_matrix *Omega;
    gretl_matrix *V;

    /* temporary storage for likelihood calculation */
    gretl_matrix *A;
    gretl_matrix *qf1;
    gretl_matrix *qf2;
};

static double Jloglik (const double *phi, void *data);

static Jwrap *jwrap_new (int neqns, int rank)
{
    Jwrap *J = malloc(sizeof *J);

    if (J == NULL) {
	return NULL;
    }

    J->T = 0;
    J->neqns = neqns;
    J->rank = rank;
    J->blen = 0;
    J->nC = 0;

    J->ll = NADBL;

    J->S00 = NULL;
    J->S01 = NULL;
    J->S11 = NULL;
    J->S00i = NULL;
    J->S11c = NULL;

    J->H = NULL;
    J->s = NULL;
    J->Hs = NULL;

    J->beta = NULL;
    J->alpha = NULL;
    J->Omega = NULL;
    J->V = NULL;

    J->A = NULL;
    J->qf1 = NULL;
    J->qf2 = NULL;

    return J;
}

static void jwrap_destroy (Jwrap *J)
{
    gretl_matrix_free(J->S00);
    gretl_matrix_free(J->S01);
    gretl_matrix_free(J->S11);
    gretl_matrix_free(J->S00i);
    gretl_matrix_free(J->S11c);

    gretl_matrix_free(J->H);
    gretl_matrix_free(J->s);

    if (J->Hs != NULL) {
	gretl_matrix_array_free(J->Hs, J->nC);
    }

    gretl_matrix_free(J->beta);
    gretl_matrix_free(J->alpha);
    gretl_matrix_free(J->Omega);
    gretl_matrix_free(J->V);

    gretl_matrix_free(J->A);
    gretl_matrix_free(J->qf1);
    gretl_matrix_free(J->qf2);

    free(J);
}

static void matrix_fill_from_array (gretl_matrix *m,
				    const double *X)
{
    int i, n = m->rows * m->cols;

    for (i=0; i<n; i++) {
	m->val[i] = X[i];
    }
}

static int make_S_matrices (Jwrap *J, const GRETL_VAR *jvar)
{
    J->S00 = gretl_matrix_copy(jvar->jinfo->Suu);
    J->S01 = gretl_matrix_copy(jvar->jinfo->Suv);
    J->S11 = gretl_matrix_copy(jvar->jinfo->Svv);

    if (J->S00 == NULL || J->S01 == NULL || J->S11 == NULL) {
	return E_ALLOC;
    }

#if JDEBUG > 1
    gretl_matrix_print(J->S00, "S00");
    gretl_matrix_print(J->S01, "S01");
    gretl_matrix_print(J->S11, "S11");
#endif    

    return 0;
}

/* Convert implicit constraint matrices R and q into explicit form.
   The matrices put in the locations **pH and **ps are newly allocated
   in this function.
*/

static int imp2exp (const gretl_matrix *Ri, const gretl_matrix *qi,
		    gretl_matrix **pH, gretl_matrix **ps)
{
    gretl_matrix *RRT = NULL;
    gretl_matrix *Tmp = NULL;
    int n = Ri->cols;
    int m = Ri->rows;
    int err = 0;

    *pH = NULL;
    *ps = NULL;

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

    *pH = gretl_matrix_right_nullspace(Ri, &err);

#if JDEBUG 
    fprintf(stderr, "right_nullspace: err = %d\n", err);
    gretl_matrix_print(*pH, "*pH");
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

 bailout:

    gretl_matrix_free(RRT);
    gretl_matrix_free(Tmp);

    return err;
}

/* Here we construct both the big block-diagonal restrictions matrix,
   J->H, and an array holding the block-diagonal elements, J->Hs,
   which will be used later in computing the variance of beta.
*/

static int set_up_restrictions (Jwrap *J, GRETL_VAR *jvar,
				const gretl_restriction_set *r)
{
    const gretl_matrix *R;
    const gretl_matrix *q;

    gretl_matrix *Ri;
    gretl_matrix *qi;
    gretl_matrix **ss;

    int nb = gretl_VECM_n_beta(jvar);
    int nC, hr = 0, hc = 0, sr = 0;
    int Rr = 0, Rc = 0;
    int i, err = 0;

    R = rset_get_R_matrix(r);
    q = rset_get_q_matrix(r);

    nC = R->cols / nb;

    J->Hs = gretl_matrix_array_alloc(nC);
    if (J->Hs == NULL) {
	return E_ALLOC;
    }

    ss = gretl_matrix_array_alloc(nC);
    if (ss == NULL) {
	return E_ALLOC;
    }

    Ri = gretl_matrix_alloc(R->rows / nb, nb);
    qi = gretl_matrix_alloc(R->rows / nb, 1);

    J->nC = 0;

    for (i=0; i<nC && !err; i++) {
	gretl_matrix_extract_matrix(Ri, R, Rr, Rc, GRETL_MOD_NONE);
	gretl_matrix_extract_matrix(qi, q, Rr,  1, GRETL_MOD_NONE);
	err = imp2exp(Ri, qi, &J->Hs[i], &ss[i]);
	if (!err) {
	    hr += J->Hs[i]->rows;
	    hc += J->Hs[i]->cols;
	    sr += ss[i]->rows;
	    Rr += R->rows / nb;
	    Rc += nC;
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
	err = gretl_matrix_inscribe_matrix(J->H, J->Hs[i], hr, hc, GRETL_MOD_NONE);
	if (!err) {
	    err = gretl_matrix_inscribe_matrix(J->s, ss[i], sr, 0, GRETL_MOD_NONE);
	}
	if (!err) {
	    hr += J->Hs[i]->rows;
	    hc += J->Hs[i]->cols;
	    sr += ss[i]->rows;
	}
    }

    gretl_matrix_array_free(ss, nC);

    return err;
}

int normalize_beta (Jwrap *J, const gretl_matrix **C, gretl_matrix *b)
{
    gretl_matrix *R = NULL;
    gretl_vector *d = NULL;
    gretl_matrix *tmp_sq = NULL;
    gretl_matrix *tmp_b = NULL;
    gretl_matrix *tmp = NULL;
    gretl_matrix *X = NULL;

    double x;
    int i, j, ii, jj, r, c;
    int nelemb, nconst = 0;
    int br = b->rows;
    int bc = J->rank;
    int err = 0;

    nelemb = br*bc;

    for (i=0; i<J->rank; i++) {
	nconst += C[i]->rows;
    }

    R = gretl_zero_matrix_new(nconst, nelemb);
    d = gretl_column_vector_alloc(nconst);

    if (R == NULL || d == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    r = 0;
    for (i=0; i<J->rank && !err; i++) {
	for (ii=0; ii<C[i]->rows; ii++) {
	    c = i * b->rows;
	    for (jj=0; jj<br; jj++) {
		x = gretl_matrix_get(C[i], ii, jj);
		gretl_matrix_set(R, r, c++, x);
	    }
	    x = gretl_matrix_get(C[i], ii, br);
	    gretl_vector_set(d, r, x);
	    r++;
	}
    }

    tmp_sq = gretl_identity_matrix_new(bc);
    tmp_b = gretl_matrix_alloc(br,bc);
    X = gretl_matrix_alloc(nconst, bc*bc);

    if (tmp_sq == NULL || tmp_b == NULL || X == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
    
    tmp = gretl_matrix_kronecker_product_new(tmp_sq, b);
    if (tmp == NULL) {
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
    for(i=0; i<bc; i++) {
	for(j=0; j<bc; j++) {
	    x = gretl_matrix_get(tmp, ii++, 0);
	    gretl_matrix_set(tmp_sq, j, i, x);
	}
    }
	    
    gretl_matrix_copy_values(tmp_b, b);
    gretl_matrix_multiply(tmp_b, tmp_sq, b);

 bailout:

    gretl_matrix_free(R);
    gretl_matrix_free(d);
    gretl_matrix_free(tmp);
    gretl_matrix_free(tmp_b);
    gretl_matrix_free(tmp_sq);
    gretl_matrix_free(X);

    return err;
}

int PCBnorm (gretl_matrix *b)
{
    gretl_matrix *c = NULL;
    gretl_matrix *tmp = NULL;
    int err = 0;

    c = gretl_matrix_alloc(b->cols, b->cols);
    tmp = gretl_matrix_copy(b);
    if (c == NULL || tmp == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_matrix_extract_matrix(c, b, 0, 0, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_invert_general_matrix(c);
    }

    if (!err) {
	err = gretl_matrix_multiply(tmp, c, b);
    }

    /* fix the upper I */
    gretl_matrix_inscribe_I(b, 0, 0, b->cols);

    gretl_matrix_free(c);
    gretl_matrix_free(tmp);

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
    gretl_matrix_print(M, "real_case0: M, again");
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

static int initval (Jwrap *J, const gretl_matrix **Cmat, gretl_matrix **pb)
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

    normalize_beta(J, Cmat, J->beta);

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
	fprintf(stderr, "*** initval: vecb->rows = %d\n", vecb->rows);
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
    int npar = 0;
    int istart, jstart;
    int i, j, err = 0;

    for (i=0; i<r; i++) { 
	npar += J->Hs[i]->cols;
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
	gretl_matrix *HiS = gretl_matrix_alloc(J->Hs[i]->cols, J->S11->cols);

	err = gretl_matrix_multiply_mod(J->Hs[i], GRETL_MOD_TRANSPOSE,
					J->S11, GRETL_MOD_NONE,
					HiS, GRETL_MOD_NONE);

	jstart = 0;
	for (j=0; j<r && !err; j++) {
	    gretl_matrix *Vij = gretl_matrix_alloc(HiS->rows, J->Hs[j]->cols);
	    double rij = gretl_matrix_get(aiom, i, j);

	    gretl_matrix_multiply(HiS, J->Hs[j], Vij);
	    gretl_matrix_multiply_by_scalar(Vij, rij);
#if JDEBUG > 1
	    gretl_matrix_print(Vij, "Vij");
	    fprintf(stderr, "inscribing Vij at %d, %d\n", istart, jstart);
#endif
	    err = gretl_matrix_inscribe_matrix(V, Vij, istart, jstart, GRETL_MOD_NONE);

	    jstart += J->Hs[j]->cols;
	    gretl_matrix_free(Vij);
	}

	istart += J->Hs[i]->cols;
	gretl_matrix_free(HiS);
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

 bailout:

    gretl_matrix_free(aiom);
    gretl_matrix_free(V);

    return err;
}

static int make_beta (Jwrap *J, const double *phi)
{
    gretl_matrix *phivec = NULL;
    gretl_matrix *tmp = NULL;
    int i, nbeta = J->beta->rows * J->beta->cols;
    int err = 0;

    if (J->blen == 0) {
	fprintf(stderr, "*** make_beta: blen = 0\n");
	return E_DATA;
    }

    phivec = gretl_column_vector_alloc(J->blen);
    tmp = gretl_column_vector_alloc(nbeta);

    if (phivec == NULL || tmp == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	matrix_fill_from_array(phivec, phi);
	err = gretl_matrix_multiply(J->H, phivec, tmp);
    }

    if (!err) {
	gretl_matrix_add_to(tmp, J->s);
    }

    if (!err) {
	for (i=0; i<nbeta; i++) {
	    J->beta->val[i] = tmp->val[i];
	}
    }

#if JDEBUG
    gretl_matrix_print(phivec, "phi");
    gretl_matrix_print(J->beta, "beta");
#endif

    gretl_matrix_free(phivec);
    gretl_matrix_free(tmp);

    return err;
}

/* BFGS callback function */

static double Jloglik (const double *phi, void *data)
{
    Jwrap *J = (Jwrap *) data;

    double ret = NADBL;
    int bcols = J->beta->cols;
    int err = 0;

    if (J->A == NULL) {
	/* temp storage not allocated yet */
	J->A = gretl_matrix_alloc(J->S01->cols, J->S01->cols);
	J->qf1 = gretl_matrix_alloc(bcols, bcols);
	J->qf2 = gretl_matrix_alloc(bcols, bcols);
	
	if (J->A == NULL || J->qf1 == NULL || J->qf2 == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = make_beta(J, phi);
    }

    if (!err) {
	gretl_matrix_copy_values(J->S11c, J->S11);
    }

    if (!err) {
	err = gretl_matrix_qform(J->S01, GRETL_MOD_TRANSPOSE,
				 J->S00i, J->A, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_matrix_subtract_from(J->S11c, J->A);
    }

    if (!err) {
	gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE,
			   J->S11c, J->qf1, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE,
				 J->S11, J->qf2, GRETL_MOD_NONE);
    }

    if (!err) {
	ret = gretl_matrix_log_determinant(J->qf1, &err);
    }

    if (!err) {
	ret -= gretl_matrix_log_determinant(J->qf2, &err);
    }

    if (err) {
	J->ll = ret = NADBL;
    } else {
	J->ll = ret = -J->T * 0.5 * ret;
	fprintf(stderr, "J->ll = %g\n", J->ll);
    }

    return ret;
}

int 
full_beta_analysis (GRETL_VAR *jvar, 
		    const gretl_restriction_set *rset,
		    PRN *prn)
{
    Jwrap *J = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *bS11b = NULL;
    gretl_matrix *S01b = NULL;

    int n = jvar->neqns;
    int rank = jrank(jvar);
    int err = 0;

    J = jwrap_new(n, rank);
    if (J == NULL) {
	return E_ALLOC;
    }

    err = make_S_matrices(J, jvar);

    if (!err) {
	err = set_up_restrictions(J, jvar, rset);
    }

#if 0 /* FIXME: not ready yet */
    if (!err) {
	err = initval(J, Cmat, &b);
#if JDEBUG
	fprintf(stderr, "after initval: err = %d\n", err);
	gretl_matrix_print(b, "b, before BFGS");
#endif
    }
#endif

    if (!err) {
	int maxit = 4000;
	double reltol = 1.0e-11;
	int fncount = 0;
	int grcount = 0;
	int nn = b->rows;

	err = BFGS_max(b->val, nn, maxit, reltol, 
		       &fncount, &grcount, Jloglik, C_LOGLIK,
		       NULL, J, (prn == NULL)? OPT_NONE : OPT_V, 
		       prn);
    }

#if JDEBUG
    if (!err) {
	printf("after BFGS: loglikelihood = %.8g\n", J->ll);
    }
#endif

    if (!err) {
	S01b = gretl_matrix_alloc(J->S01->rows, J->rank);
	bS11b = gretl_matrix_alloc(J->rank, J->rank);
	J->alpha = gretl_matrix_alloc(J->S01->rows, J->rank);
	if (S01b == NULL || bS11b == NULL || J->alpha == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = gretl_matrix_multiply(J->S01, J->beta, S01b);
    }

    if (!err) {
	err = gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE,
				 J->S11, bS11b, GRETL_MOD_NONE);
    }

    if (!err) {
	gretl_invert_symmetric_matrix(bS11b);
    }

    if (!err) {
	err = gretl_matrix_multiply(S01b, bS11b, J->alpha);
    }

    if (!err) {
	err = make_omega(J, J->alpha, J->beta);
    }

    if (!err) {
	gretl_invert_symmetric_matrix(J->Omega);
    }

    if (!err) {
	err = make_beta_variance(J);
    }

    if (!err) {
	gretl_matrix_print(J->beta, "J->beta");
	gretl_matrix_print(J->V, "J->V");
	/* FIXME attach these to jvar and nullify local ptrs */
    } 

    jwrap_destroy(J);

    gretl_matrix_free(S01b);
    gretl_matrix_free(bS11b);
    gretl_matrix_free(b);

    return err;
}
