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

#define JDEBUG 1

typedef struct Jwrap_ Jwrap;

struct Jwrap_ {
    int T;          /* length of time series used */
    int p;          /* number of equations */
    int p1;         /* number of rows in beta (>= p) */
    int rank;       /* rank, r, of VECM */
    int blen;       /* number of unrestricted coefficients in beta */
    int alen;       /* number of unrestricted coefficients in alpha */
    int noest;      /* fully constrained: no estimation required */
    int nC;         /* number of restriction matrices */
    int df;         /* degrees if freedom for LR test */
    int jr;         /* rank of Jacobian */
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

    /* homogeneous restrictions on alpha */
    gretl_matrix *G;

    /* coefficients and variances */
    gretl_matrix *beta;
    gretl_matrix *alpha;
    gretl_matrix *Omega;
    gretl_matrix *V;
    gretl_matrix *se;

    /* Jacobian */
    gretl_matrix *Jac;

    /* temp storage for beta calculation */
    gretl_matrix *phivec;

    /* alpha calculation */
    gretl_matrix *psivec;

    /* temp storage for likelihood calculation */
    gretl_matrix *qf1;
    gretl_matrix *qf2;
};

static int compute_alpha (Jwrap *J);
static int real_compute_ll (Jwrap *J);

static int make_S_matrices (Jwrap *J, const GRETL_VAR *jvar)
{
    gretl_matrix *Tmp = NULL;
    int err = 0;

    J->S00 = jvar->jinfo->S00;
    J->S01 = jvar->jinfo->S01;
    J->S11 = jvar->jinfo->S11;

    J->S00i = gretl_matrix_copy(J->S00);
    J->S11m = gretl_matrix_copy(J->S11);

    if (J->S00i == NULL || J->S11m == NULL) {
	return E_ALLOC;
    }

    J->p1 = J->S01->cols;

    Tmp = gretl_matrix_alloc(J->p1, J->p1);

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

    if (!err) {
	/* allocate beta and alpha while we're at it */
	J->beta = gretl_matrix_alloc(J->p1, J->rank);
	J->alpha = gretl_matrix_alloc(J->p, J->rank);
	if (J->beta == NULL || J->alpha == NULL) {
	    err = E_ALLOC;
	} else {
	    J->blen = J->p1 * J->rank;
	    J->alen = J->p * J->rank;
	}
    }

    if (!err) {
	/* and Omega too */
	J->Omega = gretl_matrix_alloc(J->p, J->p);
	if (J->Omega == NULL) {
	    err = E_ALLOC;
	}
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

static int J_alloc_qfs (Jwrap *J)
{
    int r = J->rank;

    J->qf1 = gretl_matrix_alloc(r, r);
    J->qf2 = gretl_matrix_alloc(r, r);

    if (J->qf1 == NULL || J->qf2 == NULL) {
	return E_ALLOC;
    }

    return 0;
}

static void jwrap_destroy (Jwrap *J)
{
    gretl_matrix_free(J->S00i);
    gretl_matrix_free(J->S11m);

    gretl_matrix_free(J->G);
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
    gretl_matrix_free(J->Jac);

    gretl_matrix_free(J->phivec);
    gretl_matrix_free(J->psivec);

    gretl_matrix_free(J->qf1);
    gretl_matrix_free(J->qf2);

    free(J);
}

static Jwrap *jwrap_new (const GRETL_VAR *jvar, int *err)
{
    Jwrap *J = malloc(sizeof *J);

    if (J == NULL) {
	return NULL;
    }

    J->T = jvar->T;
    J->p = jvar->neqns;
    J->rank = jrank(jvar);
    J->blen = 0;
    J->alen = 0;
    J->noest = 0;
    J->nC = 0;
    J->df = 0;
    J->jr = 0;

    J->ll = NADBL;

    J->S00 = NULL;
    J->S01 = NULL;
    J->S11 = NULL;

    J->S00i = NULL;
    J->S11m = NULL;

    J->H = NULL;
    J->s = NULL;
    J->h = NULL;
    J->G = NULL;

    J->beta = NULL;
    J->alpha = NULL;
    J->Omega = NULL;
    J->V = NULL;
    J->se = NULL;

    J->Jac = NULL;

    J->phivec = NULL;
    J->psivec = NULL;

    J->qf1 = NULL;
    J->qf2 = NULL;

    *err = make_S_matrices(J, jvar);

    if (!*err) {
	*err = J_alloc_qfs(J);
    }	

    if (*err) {
	jwrap_destroy(J);
	J = NULL;
    }

    return J;
}

#if SWITCHER

typedef struct switcher_ switcher;

struct switcher_ {
    gretl_matrix *K1;     /* holds kronecker product */
    gretl_matrix *K2;     /* holds kronecker product */
    gretl_matrix *I00;    /* chunk 0,0 of information matrix */
    gretl_matrix *I11;    /* chunk 1,1 of information matrix */
    gretl_matrix *TmpL;   /* used only in \Phi calculation */
    gretl_matrix *TmpR;   /* shared temp */
    gretl_matrix *TmpR1;  /* used when \alpha is restricted */
    gretl_matrix *Tmppp;  /* p x p temp matrix */
    gretl_matrix *Tmprr;  /* r x r temp */
    gretl_matrix *Tmprp;  /* r x p temp */
    gretl_matrix *Tmprp1; /* r x p1 temp */
    gretl_matrix *HK2;    /* used in \Phi calculation */
    gretl_matrix *Pi;     /* \alpha \beta' */
    gretl_matrix *lsPi;   /* vec of unrestricted \Pi */
    gretl_matrix *iOmega; /* Omega-inverse */
};

static void switcher_free (switcher *s)
{
    gretl_matrix_free(s->K1);
    gretl_matrix_free(s->K2);
    gretl_matrix_free(s->I00);
    gretl_matrix_free(s->I11);
    gretl_matrix_free(s->TmpL);
    gretl_matrix_free(s->TmpR);
    gretl_matrix_free(s->TmpR1);
    gretl_matrix_free(s->Tmppp);
    gretl_matrix_free(s->Tmprr);
    gretl_matrix_free(s->Tmprp);
    gretl_matrix_free(s->Tmprp1);
    gretl_matrix_free(s->HK2);
    gretl_matrix_free(s->Pi);
    gretl_matrix_free(s->lsPi);
    gretl_matrix_free(s->iOmega);
}

static int switcher_init (switcher *s, Jwrap *J)
{
    gretl_matrix *S11i = NULL;
    int r = J->rank;
    int err = 0;

    s->K1 = NULL; 
    s->K2 = NULL;
    s->I00 = NULL;
    s->I11 = NULL;
    s->TmpL = NULL;
    s->TmpR = NULL;
    s->TmpR1 = NULL;
    s->Tmppp = NULL;
    s->Tmprr = NULL;
    s->Tmprp = NULL;
    s->Tmprp1 = NULL;
    s->HK2 = NULL;
    s->Pi = NULL;
    s->lsPi = NULL;
    s->iOmega = NULL;

    s->K1 = gretl_matrix_alloc(J->p1 * r, J->p1 * r);
    s->K2 = gretl_matrix_alloc(J->p1 * r, J->p1 * J->p1);
    s->I00 = gretl_matrix_alloc(J->alen, J->alen);
    s->I11 = gretl_matrix_alloc(J->blen, J->blen);
    s->TmpL = gretl_matrix_alloc(J->blen, J->p * J->p1);
    s->TmpR = gretl_matrix_alloc(J->p1 * J->p1, 1);
    s->Tmppp = gretl_matrix_alloc(J->p, J->p);
    s->Tmprr = gretl_matrix_alloc(r, r);
    s->Tmprp = gretl_matrix_alloc(r, J->p);
    s->Tmprp1 = gretl_matrix_alloc(r, J->p1);
    s->HK2 = gretl_matrix_alloc(J->blen, J->p * J->p1);
    s->Pi = gretl_matrix_alloc(J->p, J->p1);
    s->iOmega = gretl_matrix_alloc(J->p, J->p);

    if (s->K1 == NULL || s->K2 == NULL || 
	s->I00 == NULL || s->I11 == NULL ||
	s->TmpL  == NULL || s->TmpR  == NULL || 
	s->Tmppp == NULL || s->Tmprr == NULL ||
	s->Tmprp == NULL || s->Tmprp1 == NULL || 
	s->HK2 == NULL || s->Pi == NULL || s->iOmega == NULL) {
	return E_ALLOC;
    }

    if (J->G != NULL) {
	s->TmpR1 = gretl_matrix_alloc(J->alen, J->p * J->p1);
	if (s->TmpR1 == NULL) {
	    return E_ALLOC;
	}
    }

    S11i = gretl_matrix_copy(J->S11);
    s->lsPi = gretl_matrix_alloc(J->p1, J->p);
    if (S11i == NULL || s->lsPi == NULL) {
	err = E_ALLOC;
    } else {
	gretl_invert_symmetric_matrix(S11i);
	gretl_matrix_multiply_mod(S11i, GRETL_MOD_NONE,
				  J->S01, GRETL_MOD_TRANSPOSE,
				  s->lsPi, GRETL_MOD_NONE);
	/* make into vec(\Pi'_{LS}) */
	gretl_matrix_reuse(s->lsPi, J->p1 * J->p, 1);
	gretl_matrix_free(S11i);
    }

    return err;
}

/* The following functions represent an attempt at implementing the
   switching algorithm as set out in Boswijk and Doornik, 2004,
   p. 455.  This does not appear to be quite right yet, but I think
   it's not very far from being right.
*/

/* 
   Update \Psi using:

   [ G'(\Omega^{-1} \otimes \beta'S_{11}\beta)G ]^{-1} 
   \times G'(\Omega^{-1} \otimes \beta'S_{11}) vec(Pi'_{LS})
*/

static int Psifun (Jwrap *J, switcher *s)
{
    int r = J->rank;
    int i, j, k = 0;
    int err = 0;

    gretl_matrix_reuse(s->K1, J->p * r, J->p * r);
    gretl_matrix_reuse(s->K2, J->p * r, J->p * J->p1);
    gretl_matrix_reuse(s->TmpR, J->alen, 1);

    /* left-hand chunk */
    err += gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE, J->S11,
			      s->Tmprr, GRETL_MOD_NONE);
    err += gretl_matrix_kronecker_product(s->iOmega, s->Tmprr, s->K1);
    if (J->G != NULL) {
	err += gretl_matrix_qform(J->G, GRETL_MOD_TRANSPOSE, s->K1,
				  s->I00, GRETL_MOD_NONE);
    } else {
	err += gretl_matrix_copy_values(s->I00, s->K1);
    }
    err += gretl_invert_symmetric_matrix(s->I00);
    if (err) {
	return err;
    }

    /* right-hand chunk */
    err += gretl_matrix_multiply_mod(J->beta, GRETL_MOD_TRANSPOSE, 
				     J->S11, GRETL_MOD_NONE,
				     s->Tmprp1, GRETL_MOD_NONE);
    err += gretl_matrix_kronecker_product(s->iOmega, s->Tmprp1, s->K2);
    if (J->G != NULL) {
	err += gretl_matrix_multiply_mod(J->G, GRETL_MOD_TRANSPOSE,
					 s->K2, GRETL_MOD_NONE,
					 s->TmpR1, GRETL_MOD_NONE);
	err += gretl_matrix_multiply(s->TmpR1, s->lsPi, s->TmpR);
    } else {
	err += gretl_matrix_multiply(s->K2, s->lsPi, s->TmpR);
    }

    /* combine */
    err += gretl_matrix_multiply(s->I00, s->TmpR, J->psivec);

    /* update alpha' from psivec */
    if (J->G != NULL) {
	gretl_matrix_reuse(s->Tmprp, J->p * r, 1);
	gretl_matrix_multiply(J->G, J->psivec, s->Tmprp);
	for (i=0; i<J->p; i++) {
	    for (j=0; j<r; j++) {
		gretl_matrix_set(J->alpha, i, j, s->Tmprp->val[k++]);
	    }
	}
	gretl_matrix_reuse(s->Tmprp, r, J->p);
    } else {
	for (i=0; i<J->p; i++) {
	    for (j=0; j<r; j++) {
		gretl_matrix_set(J->alpha, i, j, J->psivec->val[k++]);
	    }
	}
    }	

    return err;
}

/*
    Update \Phi using:

    [H'(\alpha'\Omega^{-1}\alpha \otimes S_{11})H]^{-1} \times
       H'(\alpha'\Omega^{-1} \otimes S_{11}) \times
         [vec(Pi'_{LS}) - (\alpha \otimes I_{p1})h_0]
 */

static int Phifun (Jwrap *J, switcher *s)
{
    int r = J->rank;
    int err = 0;

    gretl_matrix_reuse(s->K1, r * J->p1, r * J->p1);
    gretl_matrix_reuse(s->K2, r * J->p1, J->p * J->p1);
    gretl_matrix_reuse(s->TmpR, J->p * J->p1, 1);

    /* first big inverse */
    err += gretl_matrix_qform(J->alpha, GRETL_MOD_TRANSPOSE, s->iOmega,
			      s->Tmprr, GRETL_MOD_NONE);
    err += gretl_matrix_kronecker_product(s->Tmprr, J->S11, s->K1);
    if (J->H != NULL) {
	err += gretl_matrix_qform(J->H, GRETL_MOD_TRANSPOSE, s->K1,
				  s->I11, GRETL_MOD_NONE);
    } else {
	err += gretl_matrix_copy_values(s->I11, s->K1);
    }
    err += gretl_invert_symmetric_matrix(s->I11);
    if (err) {
	return err;
    }

    /* second chunk */
    err += gretl_matrix_multiply_mod(J->alpha, GRETL_MOD_TRANSPOSE,
				     s->iOmega, GRETL_MOD_NONE,
				     s->Tmprp, GRETL_MOD_NONE);
    err += gretl_matrix_kronecker_product(s->Tmprp, J->S11, s->K2);
    if (J->H != NULL) {
	err += gretl_matrix_multiply_mod(J->H, GRETL_MOD_TRANSPOSE,
					 s->K2, GRETL_MOD_NONE,
					 s->HK2, GRETL_MOD_NONE);
    } else {
	err += gretl_matrix_copy_values(s->HK2, s->K2);
    }

    /* combine first and second chunks */
    err += gretl_matrix_multiply(s->I11, s->HK2, s->TmpL);

    /* right-hand chunk */
    err += gretl_matrix_copy_values(s->TmpR, s->lsPi);
    if (J->s != NULL && !gretl_is_zero_matrix(J->s)) {
	gretl_matrix_reuse(s->K2, J->p * J->p1, r * J->p1);
	err += gretl_matrix_kronecker_I(J->alpha, J->p1, s->K2);
	err += gretl_matrix_multiply_mod(s->K2, GRETL_MOD_NONE,
					 J->s, GRETL_MOD_NONE,
					 s->TmpR, GRETL_MOD_DECUMULATE);
    }

    /* combine */
    err += gretl_matrix_multiply(s->TmpL, s->TmpR, J->phivec);

    /* update beta from phivec */
    if (J->H != NULL) {
	gretl_matrix_reuse(J->beta, J->p1 * r, 1);
	err += gretl_matrix_multiply(J->H, J->phivec, J->beta);
	if (!gretl_is_zero_matrix(J->s)) {
	    err += gretl_matrix_add_to(J->beta, J->s);
	}
	gretl_matrix_reuse(J->beta, J->p1, r);
    } else {
	err += gretl_matrix_copy_values_shaped(J->beta, J->phivec);
    }

    return err;
}

/* 
    Update Omega using:

    S_{00} - S_{01} \beta\alpha' - \alpha\beta' S_{10} +
       \alpha\beta' S_{11} \beta\alpha'

    then invert into iOmega.
 */

static int Omegafun (Jwrap *J, switcher *s)
{
    int err = 0;

    err += gretl_matrix_copy_values(J->Omega, J->S00);

    err += gretl_matrix_multiply_mod(J->alpha, GRETL_MOD_NONE,
				     J->beta, GRETL_MOD_TRANSPOSE,
				     s->Pi, GRETL_MOD_NONE);

    err += gretl_matrix_multiply_mod(J->S01, GRETL_MOD_NONE,
				     s->Pi, GRETL_MOD_TRANSPOSE,
				     s->Tmppp, GRETL_MOD_NONE);

    err += gretl_matrix_add_self_transpose(s->Tmppp);
    err += gretl_matrix_subtract_from(J->Omega, s->Tmppp);

    err += gretl_matrix_qform(s->Pi, GRETL_MOD_NONE, J->S11,
			      J->Omega, GRETL_MOD_CUMULATE);

    err += gretl_matrix_copy_values(s->iOmega, J->Omega);

    err += gretl_invert_symmetric_matrix(s->iOmega);

    return err;
}

static int switcher_ll (Jwrap *J, switcher *s)
{
    int err = 0;

    gretl_matrix_copy_values(s->Tmppp, J->Omega);
    J->ll = gretl_matrix_log_determinant(s->Tmppp, &err);
    if (err) {
	return err;
    }

    J->ll *= -J->T / 2.0;
    
    return 0;
}

static int switchit (Jwrap *J)
{
    switcher s;
    double lldiff, llbak = -1.0e+200;
    double tol = 4.0e-9;
    int j, jmax = 10000;
    int conv = 0;
    int err;

    err = switcher_init(&s, J);

    if (!err) {
	/* initialize alpha */
	err = compute_alpha(J);
    }

    if (!err) {
	/* initialize Omega */
	err = Omegafun(J, &s);
    }

    for (j=0; j<jmax && !err; j++) {
	fprintf(stderr, "switcher: j = %d\n", j);
	err = Psifun(J, &s);
	fprintf(stderr, " Psifun: err = %d\n", err);
	if (!err) {
	    err = Phifun(J, &s);
	    fprintf(stderr, " Phifun: err = %d\n", err);
	}
	if (!err) {
	    err = Omegafun(J, &s);
	    fprintf(stderr, " Omegafun: err = %d\n", err);
	}
	if (!err) {
	    err = switcher_ll(J, &s);
	    fprintf(stderr, " -(T/2)log|Omega| = %.8g\n", J->ll);
	}
	if (!err && j > 1) {
	    lldiff = (J->ll - llbak) / fabs(llbak);
	    fprintf(stderr, " lldiff = %.8g\n", lldiff);
	    if (lldiff < tol) {
		conv = 1;
		break;
	    }
	}
	llbak = J->ll;
    }

    if (!err && !conv) {
	err = E_NOCONV;
    }

    if (!err) {
	real_compute_ll(J);
    }

    switcher_free(&s);

    return err;
}

#endif /* SWITCHER */

/* 
   J = [(I_p \otimes \beta)G : (\alpha \otimes I_{p1})H]

   Boswijk and Doornik (2004), equation (40), page 455.
*/

static int check_jacobian (Jwrap *J)
{
    gretl_matrix *A = NULL;
    gretl_matrix *B = NULL;
    int err = 0;

    /* form both beta = H \phi + s, and alpha, for randomized 
       \theta
    */

    /* FIXME phi/psi vs J->phivec, J->psivec */

    if (J->H != NULL) {
	gretl_matrix *phi = gretl_column_vector_alloc(J->H->cols);

	gretl_matrix_random_fill(phi, D_NORMAL);
	gretl_matrix_reuse(J->beta, J->p1 * J->rank, 1);
	gretl_matrix_multiply(J->H, phi, J->beta);
	gretl_matrix_add_to(J->beta, J->s);
	gretl_matrix_reuse(J->beta, J->p1, J->rank);
	gretl_matrix_free(phi);
    } else {
	gretl_matrix_random_fill(J->beta, D_NORMAL);
    }

    if (J->G != NULL) {
	gretl_matrix *psi = gretl_column_vector_alloc(J->G->cols);

	gretl_matrix_random_fill(psi, D_NORMAL);
	gretl_matrix_reuse(J->alpha, J->p * J->rank, 1);
	gretl_matrix_multiply(J->G, psi, J->alpha);
	gretl_matrix_reuse(J->alpha, J->p, J->rank);
	gretl_matrix_free(psi);
    } else {	
	compute_alpha(J);
    }

    if (!err) {
	A = gretl_matrix_I_kronecker_new(J->p, J->beta, &err);
    }

    if (!err) {
	B = gretl_matrix_kronecker_I_new(J->alpha, J->p1, &err);
    }

    if (!err && J->G != NULL) {
	/* alpha is restricted */
	gretl_matrix *AG = gretl_matrix_multiply_new(A, J->G, &err);

	if (!err) {
	    gretl_matrix_free(A);
	    A = AG;
	}
    }

    if (!err && J->H != NULL) {
	/* beta is restricted */
	gretl_matrix *BH = gretl_matrix_multiply_new(B, J->H, &err);

	if (!err) {
	    gretl_matrix_free(B);
	    B = BH;
	}
    }

    if (!err) {
	J->Jac = gretl_matrix_col_concat(A, B, &err);
    }

    if (!err) {
	J->jr = gretl_matrix_rank(J->Jac, &err);
#if 1
	fprintf(stderr, "rank of Jacobian = %d\n", J->jr);
	gretl_matrix_print(J->Jac, "Jacobian");
#endif
    }

    gretl_matrix_free(A);
    gretl_matrix_free(B);

    return err;
}

/* Convert the implicit constraint matrices R and q into explicit form.
   The matrices put in the locations *ph and *ps are newly allocated
   in this function.
*/

static int imp2exp (gretl_matrix *Ri, const gretl_matrix *qi,
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

    if (rr == n) {
	/* special case: \beta_i is completely specified by the
	   restrictions */
	*ps = gretl_matrix_copy(qi);
	err = gretl_LU_solve(Ri, *ps);
	goto bailout;
    }	
	
    /* the standard case: \beta_i only partially constrained */

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
    fprintf(stderr, "right_nullspace: err = %d\n\n", err);
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
	*ps = gretl_matrix_multiply_new(Tmp, qi, &err);
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

/* produce the \beta matrix in the case where it is fully
   constrained (and therefore no estimation is needed)
*/

static int solve_for_beta (Jwrap *J, 
			   const gretl_matrix *R,
			   const gretl_matrix *q)
{
    gretl_matrix *b;
    int err = 0;

    b = gretl_matrix_copy(q);
    if (b == NULL) {
	return E_ALLOC;
    }

    if (!gretl_is_identity_matrix(R)) {
	gretl_matrix *Rcpy = gretl_matrix_copy(R);

	if (Rcpy == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gretl_LU_solve(Rcpy, b);
	    gretl_matrix_free(Rcpy);
	}
    }

    if (!err) {
	err = gretl_matrix_copy_values_shaped(J->beta, b);
    }

    if (!err) {
	J->noest = 1;
    }

    gretl_matrix_free(b);
    
    return err;
}

static gretl_matrix *augmented_H (const gretl_matrix *H,
				  const gretl_matrix *s,
				  int *err)
{
    gretl_matrix *H1;

    if (H != NULL) {
	H1 = gretl_matrix_col_concat(H, s, err);
    } else {
	H1 = gretl_matrix_copy(s);
    }

    return H1;
}

static void rank_check_message (int i, int j, int *k, int r, int rmin,
				PRN *prn)
{
    i++, j++; /* convert to 1-based */

    if (k == NULL) {
	pprintf(prn, "Rank of R%d * H%d = %d", i, j, r);
    } else {
	int p;

	pprintf(prn, "Rank of R%d * (H%d", i, j);
	for (p=1; p<rmin; p++) {
	    pprintf(prn, ":H%d", k[p] + 1);
	}
	pprintf(prn, ") = %d", r);
    }

    pprintf(prn, ", should be >= %d\n", rmin);
}

static int rank_check (const gretl_matrix *R, const gretl_matrix *H,
		       int i, int j, int *k, int rmin, PRN *prn)
{
    gretl_matrix *RH;
    int r, err = 0;

    RH = gretl_matrix_multiply_new(R, H, &err);

    if (!err) {
	r = gretl_matrix_rank(RH, &err);
	if (r < rmin) {
	    rank_check_message(i, j, k, r, rmin, prn);
	    err = E_NOIDENT;
	}
#if JDEBUG
	if (r >= rmin) {
	    rank_check_message(i, j, k, r, rmin, prn);
	}
#endif
	gretl_matrix_free(RH);
    } 

    return err;
}

/* Recursive routine to construct the rank tests for (R_i * H+), where
   H+ is composed of the columns of 2 or more of the per-restriction
   H_j matrices -- e.g. R_1 * (H_2~H_3).
*/

static int extra_check (int i, int j, int *k, int d, int dmax,
			gretl_matrix **R, gretl_matrix **H,
			PRN *prn)
{
    gretl_matrix *Htmp = NULL;
    int p, kd0 = k[d-1] + 1;
    int err = 0;

    for (k[d]=kd0; k[d]<=dmax && !err; k[d] += 1) {
	if (k[d] == i) continue;

	Htmp = gretl_matrix_copy(H[j]);

	if (Htmp == NULL) {
	    err = E_ALLOC;
	} else {
	    for (p=1; p<=d && !err; p++) {
#if JDEBUG > 1
		fprintf(stderr, "extra_check: i=%d, j=%d, d=%d, dmax=%d,"
			" Htmp=%p, p=%d, k[p]=%d\n", i, j, d, dmax,
			(void *) Htmp, p, k[p]);
#endif
		err = gretl_matrix_inplace_colcat(Htmp, H[k[p]], NULL);
	    }	
	    if (!err) {
		err = rank_check(R[i], Htmp, i, j, k, d + 1, prn);
	    }
	    gretl_matrix_free(Htmp);
	}

	if (!err && d < dmax - 1) {
	    err = extra_check(i, j, k, d + 1, dmax, R, H, prn);
	}
    }

    return err;
}

#if SWITCHER

/* Doornik's approach */

static int 
identification_check (Jwrap *J, gretl_matrix **R, 
		      gretl_matrix **ss, PRN *prn)
{
    int npar = J->H->cols + J->alen;
    int err;

    err = check_jacobian(J);
    if (!err && J->jr < npar) {
	fprintf(stderr, "Rank(Jacobian) = %d < %d\n", J->jr, npar);
    }

    return err;
}

#else

/* See Johansen, Journal of Econometrics, 1995 */

static int 
identification_check (Jwrap *J, gretl_matrix **R, 
		      gretl_matrix **ss, PRN *prn)
{
    gretl_matrix **Rtmp = NULL;
    gretl_matrix **Htmp = NULL;
    int i, j, *k = NULL;
    int err = 0;

    if (J->nC < 2) {
	return 0;
    }

    if (J->nC > 2) {
	/* extra indices are needed */
	k = malloc((J->nC - 1) * sizeof *k);
	if (k == NULL) {
	    return E_ALLOC;
	}
    }

    /* construct arrays of (possibly modified) R and H matrices for
       testing */
    
    Rtmp = gretl_matrix_array_alloc(J->nC);
    Htmp = gretl_matrix_array_alloc(J->nC);
    if (Rtmp == NULL || Htmp == NULL) {
	err = E_ALLOC;
    }

#if JDEBUG
    fprintf(stderr, "identification_check: Rtmp and Htmp: size %d, err = %d\n",
	    J->nC, err);
#endif

    for (i=0; i<J->nC && !err; i++) {
	if (gretl_is_zero_matrix(ss[i])) {
	    /* restriction is homogeneous */
	    if (J->h[i] == NULL) {
		err = E_DATA;
	    } else {
		Rtmp[i] = gretl_matrix_copy(R[i]);
		Htmp[i] = gretl_matrix_copy(J->h[i]);
		if (Htmp[i] == NULL || Rtmp[i] == NULL) {
		    err = E_ALLOC;
		}
	    }
	} else {
	    /* non-homogeneous: augment H with s */
	    Htmp[i] = augmented_H(J->h[i], ss[i], &err);
	    if (!err) {
		Rtmp[i] = gretl_matrix_left_nullspace(Htmp[i], 
						      GRETL_MOD_TRANSPOSE, 
						      &err);
		if (err == E_DATA) {
		    /* augmented H is of full rank: give up */
		    err = E_NOIDENT;
		    break;
		}
	    }
	}
    }	

    if (!err) {
	/* conduct the Johansen rank tests on R_i * H_j, etc. */
	for (i=0; i<J->nC && !err; i++) {
	    for (j=0; j<J->nC && !err; j++) {
		if (i == j || Rtmp[i] == NULL || Htmp[j] == NULL) {
		    continue;
		}
		err = rank_check(Rtmp[i], Htmp[j], i, j, NULL, 1, prn);
		if (!err && k != NULL) {
		    k[0] = j;
		    err = extra_check(i, j, k, 1, J->nC - 1, 
				      Rtmp, Htmp, prn);
		}
	    }
	}
    } 

    gretl_matrix_array_free(Rtmp, J->nC);
    gretl_matrix_array_free(Htmp, J->nC);
    free(k);

    return err;
}

#endif

static int count_nC (Jwrap *J, const gretl_matrix *R, int nb)
{
    int nC = 0, r0 = 0;
    int i, ir;

    for (i=0; i<J->rank; i++) {
	ir = get_Ri_rows(R, nb, i, &r0);
	if (ir > 0) {
	    nC++;
	}
    }

    return nC;
}

/* Here we construct both the big block-diagonal restrictions matrix,
   J->H, and an array holding the block-diagonal elements, J->h,
   which will be used later in computing the variance of beta.
*/

static int set_up_H_and_s (Jwrap *J, GRETL_VAR *jvar,
			   const gretl_restriction_set *rset,
			   PRN *prn)
{
    const gretl_matrix *R;
    const gretl_matrix *q;

    gretl_matrix *Ri;
    gretl_matrix *qi;
    gretl_matrix **Rarr;
    gretl_matrix **ss;

    int nC, hr = 0, hc = 0, sr = 0;
    int Rr = 0, Rc = 0, irmax = 0;
    int r0, ir;
    int i, err = 0;

    R = rset_get_R_matrix(rset);
    q = rset_get_q_matrix(rset);

#if JDEBUG
    gretl_matrix_print(R, "R, in set_up_H_and_s");
    gretl_matrix_print(q, "q, in set_up_H_and_s");
#endif

    if (R->rows == J->p1 * J->rank) {
	/* number of restrictions = total betas */
	J->nC = J->rank;
	return solve_for_beta(J, R, q);
    } else if (R->rows < J->rank * J->rank) {
	pprintf(prn, "R->rows = %d, should be >= %d\n",
		R->rows, J->rank * J->rank);
	return E_NOIDENT;
    } else if ((nC = count_nC(J, R, J->p1)) < J->rank) {
	pprintf(prn, "R blocks = %d, should be %d\n", nC, J->rank);
	return E_NOIDENT;
    }

    nC = J->rank;

    J->h = gretl_matrix_array_alloc(nC);
    if (J->h == NULL) {
	return E_ALLOC;
    }

    ss = gretl_matrix_array_alloc(nC);
    Rarr = gretl_matrix_array_alloc(nC);
    if (ss == NULL || Rarr == NULL) {
	free(ss);
	free(Rarr);
	return E_ALLOC;
    }

    /* find max rows per restriction sub-matrix */
    r0 = 0;
    for (i=0; i<nC; i++) {
	ir = get_Ri_rows(R, J->p1, i, &r0);
	if (ir > irmax) {
	    irmax = ir;
	}
    }

    Ri = gretl_matrix_alloc(irmax, J->p1);
    qi = gretl_matrix_alloc(irmax, 1);

    J->nC = r0 = 0;

    for (i=0; i<nC && !err; i++) {
	ir = get_Ri_rows(R, J->p1, i, &r0);
	gretl_matrix_reuse(Ri, ir, 0);
	gretl_matrix_reuse(qi, ir, 0);
	gretl_matrix_extract_matrix(Ri, R, Rr, Rc, GRETL_MOD_NONE);
	gretl_matrix_extract_matrix(qi, q, Rr,  0, GRETL_MOD_NONE);
	err = imp2exp(Ri, qi, &J->h[i], &ss[i]);
#if JDEBUG
	gretl_matrix_print(J->h[i], "H[i], in set_up_H_and_s");
	gretl_matrix_print(ss[i], "ss[i], in set_up_H_and_s");
#endif
	if (!err) {
	    if (J->h[i] == NULL) {
		/* may happen if the i-th cointegration vector
		   is fully restricted */
		hr += Ri->cols;
	    } else {
		Rarr[i] = gretl_matrix_copy(Ri);
		hr += J->h[i]->rows;
		hc += J->h[i]->cols;
	    }
	    sr += ss[i]->rows;
	    Rr += ir;
	    Rc += J->p1;
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
	    J->blen = hc;
	    hr = hc = sr = 0;
	}
    }
    
    for (i=0; i<nC && !err; i++) {
	err = gretl_matrix_inscribe_matrix(J->s, ss[i], sr, 0, 
					   GRETL_MOD_NONE);
	sr += ss[i]->rows;
	if (!err && J->h[i] != NULL) {
	    err = gretl_matrix_inscribe_matrix(J->H, J->h[i], hr, hc, 
					       GRETL_MOD_NONE);
	    if (!err) {
		hr += J->h[i]->rows;
		hc += J->h[i]->cols;
	    }
	}
    }

    if (!err) {
	err = identification_check(J, Rarr, ss, prn);
    }

#if JDEBUG
    if (!err) {
	gretl_matrix_print(J->H, "J->H");
	gretl_matrix_print(J->s, "J->s");
    }
#endif

    gretl_matrix_array_free(ss, nC);
    gretl_matrix_array_free(Rarr, nC);

    return err;
}

static int set_up_G (Jwrap *J, GRETL_VAR *jvar,
		     const gretl_restriction_set *rset,
		     PRN *prn)
{
    const gretl_matrix *Ra = rset_get_Ra_matrix(rset);
    const gretl_matrix *qa = rset_get_qa_matrix(rset);
    int err = 0;

#if JDEBUG
    gretl_matrix_print(Ra, "Ra, in set_up_G");
    gretl_matrix_print(qa, "qa, in set_up_G");
#endif

    if (!gretl_is_zero_matrix(qa)) {
	pprintf(prn, "alpha restriction is not homogeneous: not supported\n");
	return E_NOTIMP;
    }

    if (J->rank > 1 && Ra->cols == J->p) {
	/* common alpha restriction */
	gretl_matrix *Rtmp = gretl_matrix_I_kronecker_new(J->rank, Ra, &err);

	if (!err) {
	    J->G = gretl_matrix_right_nullspace(Rtmp, &err);
	    gretl_matrix_free(Rtmp);
	}
    } else {
	J->G = gretl_matrix_right_nullspace(Ra, &err);
    }

    if (!err) {
	J->alen = J->G->cols;
    }

    return err;
}

static int 
normalize_initial_beta (Jwrap *J, const gretl_restriction_set *rset, 
			gretl_matrix *b)
{
    const gretl_matrix *R = rset_get_R_matrix(rset);
    const gretl_matrix *d = rset_get_q_matrix(rset);

    gretl_matrix *tmp = NULL;
    gretl_matrix *tmp_sq = NULL;
    gretl_matrix *tmp_b = NULL;
    gretl_matrix *X = NULL;

    double x;
    int i, j, ii;
    int br = b->rows;
    int bc = J->rank;
    int bc2 = bc * bc;
    int err = 0;

    if (bc2 > d->rows) {
	fprintf(stderr, "*** normalize_initial_beta: df = %d\n", d->rows - bc2);
	return 0;
    }

    tmp_b = gretl_matrix_alloc(br, bc);
    tmp_sq = gretl_matrix_alloc(bc, bc);
    X = gretl_matrix_alloc(R->rows, bc2);
    if (tmp_b == NULL || tmp_sq == NULL || X == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    tmp = gretl_matrix_I_kronecker_new(bc, b, &err);
    if (err) {
	goto bailout;
    }

    gretl_matrix_multiply(R, tmp, X);
    gretl_matrix_reuse(tmp, bc2, 1);

    err = gretl_matrix_multi_ols(d, X, tmp, NULL);
    if (err) {
	fprintf(stderr, "beta initialization: gretl_matrix_multi_ols failed\n");
	err = 0;
	goto bailout;
    }

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
    gretl_matrix_free(tmp_sq);
    gretl_matrix_free(tmp_b);
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
	gretl_matrix_copy_values(J->beta, M);
    }

    gretl_matrix_free(Tmp);
    gretl_matrix_free(M);
    gretl_matrix_free(evals);

    return err;
}

static int initval (Jwrap *J, const gretl_restriction_set *rset,
		    gretl_matrix **pb)
{
    gretl_matrix *HHi = NULL;
    gretl_matrix *vecb = NULL;
    gretl_matrix *tmp = NULL;
    int err;

    err = case0(J);
    if (err) {
	return err;
    }

    if (rset_get_R_matrix(rset) == NULL) {
	return 0;
    }

    err = normalize_initial_beta(J, rset, J->beta);
    if (err) {
	return err;
    }

#if JDEBUG
    gretl_matrix_print(J->beta, "initval: 'normalized' beta");
#endif

    HHi = gretl_matrix_alloc(J->blen, J->blen);
    vecb = gretl_column_vector_alloc(J->H->rows);
    tmp = gretl_column_vector_alloc(J->blen);

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
	gretl_matrix_print(vecb, "initval: final vecb");
#endif
	*pb = vecb;
    } else {
	gretl_matrix_free(vecb);
    }

    gretl_matrix_free(HHi);
    gretl_matrix_free(tmp);
    
    return err;
}

#if !SWITCHER

static int make_omega (Jwrap *J,
		       const gretl_matrix *alpha, 
		       const gretl_matrix *beta)
{
    gretl_matrix *tmp = NULL;
    int err = 0;

    tmp = gretl_matrix_alloc(J->rank, J->rank);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_copy_values(J->Omega, J->S00);

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

#endif

static int make_zero_variance (Jwrap *J)
{
    int npar = J->p1 * J->rank;

    J->V = gretl_zero_matrix_new(npar, npar);
    J->se = gretl_zero_matrix_new(J->p1, J->rank);

    if (J->V == NULL || J->se == NULL) {
	return E_ALLOC;
    }

    return 0;
}

static int make_beta_variance (Jwrap *J)
{
    const gretl_matrix *a = J->alpha;

    gretl_matrix *V = NULL;
    gretl_matrix *aiom = NULL;

    int r = J->rank;
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
	gretl_matrix *HiS;

	if (J->h[i] == NULL) {
	    continue;
	}

	HiS = gretl_matrix_alloc(J->h[i]->cols, J->S11->cols);
	if (HiS == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gretl_matrix_multiply_mod(J->h[i], GRETL_MOD_TRANSPOSE,
					    J->S11, GRETL_MOD_NONE,
					    HiS, GRETL_MOD_NONE);
	}

	jstart = 0;
	for (j=0; j<r && !err; j++) {
	    if (J->h[j] != NULL) {
		double rij = gretl_matrix_get(aiom, i, j);
		gretl_matrix *Vij;

		Vij = gretl_matrix_multiply_new(HiS, J->h[j], &err);
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

	J->se = gretl_matrix_alloc(J->p1, r);
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
    int i, err = 0;

    gretl_matrix_reuse(J->beta, J->p1 * J->rank, 1);

    for (i=0; i<J->blen; i++) {
	J->phivec->val[i] = phi[i];
    }
    err = gretl_matrix_multiply(J->H, J->phivec, J->beta);

    if (!err) {
	gretl_matrix_add_to(J->beta, J->s);
    }

    gretl_matrix_reuse(J->beta, J->p1, J->rank);

#if JDEBUG > 1
    gretl_matrix_print(J->phivec, "phi");
    gretl_matrix_print(J->beta, "beta");
#endif

    return err;
}

static int real_compute_ll (Jwrap *J)
{
    double ll0 = J->ldS00;
    int err = 0;

    err = gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE,
			     J->S11m, J->qf1, GRETL_MOD_NONE);

    if (!err) {
	err = gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE,
				 J->S11, J->qf2, GRETL_MOD_NONE);
    }

    if (!err) {
	ll0 += gretl_matrix_log_determinant(J->qf1, &err);
    }

    if (!err) {
	ll0 -= gretl_matrix_log_determinant(J->qf2, &err);
    }

    if (err) {
	J->ll = NADBL;
    } else {
	J->ll = -J->T * 0.5 * (J->p * (1.0 + LN_2_PI) + ll0);
    }

    return err;
}

/* BFGS callback function */

static double Jloglik (const double *phi, void *data)
{
    Jwrap *J = (Jwrap *) data;
    int err;

    err = make_beta(J, phi);

    if (!err) {
	err = real_compute_ll(J);
    }

    return J->ll;
}

/* See Johansen pp. 106-112 */

static void set_LR_df (Jwrap *J, GRETL_VAR *jvar)
{
    int p = gretl_matrix_rows(J->beta);
    int r = J->rank;
    int i, si;

    J->df = 0;

    for (i=0; i<J->nC; i++) {
	si = 0;
	if (J->h != NULL && J->h[i] != NULL) {
	    si = J->h[i]->cols;
	}
	J->df += p - r - si;
    }

    /* system was subject to a prior restriction */
    J->df -= jvar->jinfo->bdf;

    if (J->df < 0) {
	fprintf(stderr, "*** warning: set_LR_df gave df = %d\n", J->df);
    }
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

static int simann (Jwrap *J, gretl_matrix *b, gretlopt opt, PRN *prn)
{
    int i, SAiter = 4096;
    double f0, f1;
    double fbest, fworst;
    double rndu;
    int jump;

    gretl_matrix *b0 = NULL;
    gretl_matrix *b1 = NULL;
    gretl_matrix *d = NULL;

    double Temp = 1.0;
    double radius = 1.0;
    int hdr = 0, err = 0;

    b0 = gretl_matrix_copy(b);
    b1 = gretl_matrix_copy(b);
    d = gretl_column_vector_alloc(b->rows);

    if (b0 == NULL || b1 == NULL || d == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    f0 = fbest = fworst = Jloglik(b->val, J);

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
		gretl_matrix_copy_values(b, b0);
		if (opt & OPT_V) {
		    if (!hdr) {
			pputs(prn, "\nSimulated annealing:\n");
			pprintf(prn, "%6s %12s %12s %12s\n",
				"iter", "temp", "radius", "fbest");
			hdr = 1;
		    }
		    pprintf(prn, "%6d %#12.6g %#12.6g %#12.6g\n", 
			    i, Temp, radius, fbest);
		}
	    } else if (f0 < fworst) {
		fworst = f0;
	    }
	} else {
	    gretl_matrix_copy_values(b1, b0);
	    f1 = f0;
	}

	Temp *= 0.999;
	radius *= 0.9999;
    }

    if (hdr) {
	pputc(prn, '\n');
    }
    
    if (fbest - fworst < 1.0e-9) {
	pprintf(prn, "Warning: likelihood seems to be flat\n");
    }

 bailout:

    gretl_matrix_free(b0);
    gretl_matrix_free(b1);
    gretl_matrix_free(d);

    return err;
}

static int compute_alpha (Jwrap *J)
{
    gretl_matrix *S01b = NULL;
    int err = 0;

    S01b = gretl_matrix_multiply_new(J->S01, J->beta, &err);

    if (!err) {
	err = gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE,
				 J->S11, J->qf1, GRETL_MOD_NONE);
    }

    if (!err) {
	gretl_invert_symmetric_matrix(J->qf1);
    }

    if (!err) {
	gretl_matrix_multiply(S01b, J->qf1, J->alpha);
    }

    gretl_matrix_free(S01b);

    return err;
}

static int allocate_phivec (Jwrap *J)
{
    if (J->H != NULL) {
	J->phivec = gretl_column_vector_alloc(J->H->cols);
    } else {
	J->phivec = gretl_column_vector_alloc(J->p1 * J->rank);
    }

    return (J->phivec == NULL)? E_ALLOC : 0;
}

static int allocate_psivec (Jwrap *J)
{
    if (J->G != NULL) {
	J->psivec = gretl_column_vector_alloc(J->G->cols);
    } else {
	J->psivec = gretl_column_vector_alloc(J->p * J->rank);
    }

    return (J->psivec == NULL)? E_ALLOC : 0;
}

/* public entry point */

int general_vecm_analysis (GRETL_VAR *jvar, 
			   const gretl_restriction_set *rset,
			   const DATAINFO *pdinfo,
			   gretlopt opt,
			   PRN *prn)
{
    Jwrap *J = NULL;
    gretl_matrix *b = NULL;
    int err = 0;

    J = jwrap_new(jvar, &err);
    if (err) {
	return err;
    }

    if (rset_VECM_bcols(rset) > 0) {
	err = set_up_H_and_s(J, jvar, rset, prn);
    }

    if (!err) {
	err = allocate_phivec(J);
    }

    if (rset_VECM_acols(rset) > 0) {
	err = set_up_G(J, jvar, rset, prn);
    }

    if (!err) {
	err = allocate_psivec(J);
    }    

    if (!err && J->noest) {
	/* nothing to be estimated */
	err = real_compute_ll(J);
	goto skipest;
    }

    if (!err) {
	err = initval(J, rset, &b);
#if JDEBUG
	fprintf(stderr, "after initval: err = %d\n", err);
	gretl_matrix_print(b, "b, before BFGS");
#endif
    }

#if SWITCHER
    if (!err) {
	err = switchit(J);
    }
#else
    if (!err) {
	err = simann(J, b, opt, prn);
    }    
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
#endif

 skipest:

#if !SWITCHER
    if (!err) {
	err = compute_alpha(J);
    }

    if (!err) {
	err = make_omega(J, J->alpha, J->beta);
    }
#endif

    if (opt & OPT_F) {
	gretl_matrix_free(jvar->S);
	jvar->S = gretl_matrix_copy(J->Omega);
    }

    /* FIXME: below here, when SWITCHER is non-zero */

    if (!err) {
	gretl_invert_symmetric_matrix(J->Omega);
    }

    if (!err) {
	if (J->noest) {
	    err = make_zero_variance(J);
	} else {
	    err = make_beta_variance(J);
	}
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
    gretl_matrix_free(b);

    return err;
}
