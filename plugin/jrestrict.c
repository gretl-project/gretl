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
#include "libset.h"
#include "jprivate.h"

#define JDEBUG 1

typedef struct gradhelper_ gradhelper;
typedef struct Jwrap_ Jwrap;

struct Jwrap_ {
    int bfgs;       /* use LBFGS */    
    int T;          /* length of time series used */
    int p;          /* number of equations */
    int p1;         /* number of rows in beta (>= p) */
    int r;          /* cointegrating rank of VECM */
    int blen;       /* number of unrestricted coefficients in beta */
    int alen;       /* number of unrestricted coefficients in alpha */
    int df;         /* degrees if freedom for LR test */
    int jr;         /* rank of Jacobian */
    double llk;     /* constant in log-likelihood */
    double ll;      /* log-likelihood */

    /* moment matrices and copies */
    const gretl_matrix *S00;
    const gretl_matrix *S01;
    const gretl_matrix *S11;

    gretl_matrix *S00i;

    /* restrictions on beta */
    gretl_matrix *H;
    gretl_matrix *h0;

    /* homogeneous restrictions on alpha */
    gretl_matrix *G;

    /* coefficients and variances */
    gretl_matrix *beta;
    gretl_matrix *alpha;
    gretl_matrix *Pi;
    gretl_matrix *lsPi;
    gretl_matrix *Omega;
    gretl_matrix *iOmega;
    gretl_matrix *Vb;
    gretl_matrix *Va;
    gretl_matrix *bse;
    gretl_matrix *ase;

    /* information matrix diagonal blocks */
    gretl_matrix *I00;
    gretl_matrix *I11;

    /* free parameter vector (BFGS or jittering only) */
    gretl_matrix *theta;

    /* temp storage for beta calculation */
    gretl_matrix *phi;

    /* alpha calculation */
    gretl_matrix *psi;

    /* temporary storage */
    gretl_matrix *qf1;
    gretl_matrix *qf2;
    gretl_matrix *Tmppp;
    gretl_matrix *Tmprp;

    /* for analytical derivatives in LBFGS */
    gradhelper *ghelper;
};

static int J_compute_alpha (Jwrap *J);
static int make_beta_se (Jwrap *J);
static int make_alpha_se (Jwrap *J);
static int phi_from_beta (Jwrap *J);
static void gradhelper_free (gradhelper *g);

static int make_S_matrices (Jwrap *J, const GRETL_VAR *jvar)
{
    gretl_matrix *Tmp = NULL;
    int err = 0;

    J->S00 = jvar->jinfo->S00;
    J->S01 = jvar->jinfo->S01;
    J->S11 = jvar->jinfo->S11;

    J->S00i = gretl_matrix_copy(J->S00);

    if (J->S00i == NULL) {
	return E_ALLOC;
    }

    J->p1 = J->S01->cols;

    Tmp = gretl_matrix_alloc(J->p1, J->p1);

    if (Tmp == NULL) {
	return E_ALLOC;
    }

    if (!err) {
	gretl_matrix_copy_values(J->S00i, J->S00);
	err = gretl_invert_symmetric_matrix(J->S00i);
    }

    if (!err) {
	err = gretl_matrix_qform(J->S01, GRETL_MOD_TRANSPOSE,
				 J->S00i, Tmp, GRETL_MOD_NONE);
    }

    if (!err) {
	/* allocate beta, alpha, etc. while we're at it */
	J->beta = gretl_matrix_alloc(J->p1, J->r);
	J->alpha = gretl_matrix_alloc(J->p, J->r);
	J->Pi = gretl_matrix_alloc(J->p, J->p1);
	J->Omega = gretl_matrix_alloc(J->p, J->p);
	if (J->beta == NULL || J->alpha == NULL || 
	    J->Pi == NULL || J->Omega == NULL) {
	    err = E_ALLOC;
	} else {
	    J->blen = J->p1 * J->r;
	    J->alen = J->p * J->r;
	}
    }

#if JDEBUG > 1
    gretl_matrix_print(J->S00, "S00");
    gretl_matrix_print(J->S01, "S01");
    gretl_matrix_print(J->S11, "S11");
    gretl_matrix_print(J->S00i, "S00i");
#endif    

    gretl_matrix_free(Tmp);

    return err;
}

static int J_alloc_aux (Jwrap *J)
{
    clear_gretl_matrix_err();

    J->qf1 = gretl_matrix_alloc(J->r, J->r);
    J->qf2 = gretl_matrix_alloc(J->r, J->r);
    J->Tmppp = gretl_matrix_alloc(J->p, J->p);
    J->Tmprp = gretl_matrix_alloc(J->r, J->p);

    return get_gretl_matrix_err();
}

static int allocate_info_blocks_etc (Jwrap *J)
{
    int err = 0;

    if (J->G != NULL) {
	J->I00 = gretl_matrix_alloc(J->alen, J->alen);
	if (J->I00 == NULL) {
	    err = E_ALLOC;
	}
    }

    if (J->blen > 0 && !err) {
	J->I11 = gretl_matrix_alloc(J->blen, J->blen);
	if (J->I11 == NULL) {
	    err = E_ALLOC;
	}	
    }

    if (J->iOmega == NULL && !err) {
	J->iOmega = gretl_matrix_alloc(J->p, J->p);
	if (J->iOmega == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static void jwrap_destroy (Jwrap *J)
{
    gretl_matrix_free(J->S00i);

    gretl_matrix_free(J->G);
    gretl_matrix_free(J->H);
    gretl_matrix_free(J->h0);

    gretl_matrix_free(J->beta);
    gretl_matrix_free(J->alpha);
    gretl_matrix_free(J->Pi);
    gretl_matrix_free(J->lsPi);
    gretl_matrix_free(J->Omega);
    gretl_matrix_free(J->iOmega);
    gretl_matrix_free(J->Vb);
    gretl_matrix_free(J->Va);
    gretl_matrix_free(J->bse);
    gretl_matrix_free(J->ase);

    gretl_matrix_free(J->I00);
    gretl_matrix_free(J->I11);

    gretl_matrix_free(J->phi);
    gretl_matrix_free(J->psi);
    gretl_matrix_free(J->theta);

    gretl_matrix_free(J->qf1);
    gretl_matrix_free(J->qf2);
    gretl_matrix_free(J->Tmppp);
    gretl_matrix_free(J->Tmprp);

    gradhelper_free(J->ghelper);

    free(J);
}

enum {
    BFGS_NUMERIC = 1,
    BFGS_ANALYTIC 
};

static int use_lbfgs (gretlopt opt)
{
    int ret = 0;

    /* could use BFGS_NUMERIC instead */

    if (opt & OPT_L) {
	ret = BFGS_ANALYTIC;
    } else if (getenv("GRETL_VECM_LBFGS") != NULL) {
	ret = BFGS_ANALYTIC;
    }

    return ret;
}

static Jwrap *jwrap_new (const GRETL_VAR *jvar, gretlopt opt, int *err)
{
    Jwrap *J = malloc(sizeof *J);

    if (J == NULL) {
	return NULL;
    }

    J->bfgs = use_lbfgs(opt);

    J->T = jvar->T;
    J->p = jvar->neqns;
    J->r = jrank(jvar);
    J->blen = 0;
    J->alen = 0;
    J->df = 0;
    J->jr = 0;

    J->llk = J->T * 0.5 * (J->p * (1.0 + LN_2_PI));
    J->ll = NADBL;

    J->S00 = NULL;
    J->S01 = NULL;
    J->S11 = NULL;

    J->S00i = NULL;

    J->H = NULL;
    J->h0 = NULL;
    J->G = NULL;

    J->beta = NULL;
    J->alpha = NULL;
    J->Pi = NULL;
    J->lsPi = NULL;
    J->Omega = NULL;
    J->iOmega = NULL;
    J->Vb = NULL;
    J->Va = NULL;
    J->bse = NULL;
    J->ase = NULL;

    J->I00 = NULL;
    J->I11 = NULL;

    J->phi = NULL;
    J->psi = NULL;
    J->theta = NULL;

    J->qf1 = NULL;
    J->qf2 = NULL;
    J->Tmppp = NULL;
    J->Tmprp = NULL;

    J->ghelper = NULL;

    *err = make_S_matrices(J, jvar);

    if (!*err) {
	*err = J_alloc_aux(J);
    }	

    if (*err) {
	jwrap_destroy(J);
	J = NULL;
    }

    return J;
}

/* vectorize src into targ */

static void 
vec_simple (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, n = src->rows * src->cols;

    for (i=0; i<n; i++) {
	targ->val[i] = src->val[i];
    }
}

/* vectorize the transpose of src into targ */

static void 
vec_transpose (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, j, k = 0;

    for (i=0; i<src->rows; i++) {
	for (j=0; j<src->cols; j++) {
	    targ->val[k++] = gretl_matrix_get(src, i, j);
	}
    }
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
	J->blen = 0;
    }

    gretl_matrix_free(b);
    
    return err;
}

static int psi_from_alpha (Jwrap *J)
{
    gretl_matrix *GG = NULL;
    gretl_matrix *GGG = NULL;
    int err = 0;

    if (J->G == NULL) {
	/* just vectorize alpha' into psi */
	vec_transpose(J->psi, J->alpha);
	return 0;
    }	

    GG = gretl_matrix_alloc(J->G->cols, J->G->cols);
    GGG = gretl_matrix_alloc(J->G->cols, J->G->rows);

    if (GG == NULL || GGG == NULL) {
	gretl_matrix_free(GG);
	gretl_matrix_free(GGG);
	return E_ALLOC;
    }

    gretl_matrix_multiply_mod(J->G, GRETL_MOD_TRANSPOSE,
			      J->G, GRETL_MOD_NONE,
			      GG, GRETL_MOD_NONE);

    err = gretl_invert_symmetric_matrix(GG);
    if (err) {
	goto bailout;
    }

    gretl_matrix_multiply_mod(GG, GRETL_MOD_NONE,
			      J->G, GRETL_MOD_TRANSPOSE,
			      GGG, GRETL_MOD_NONE);

    gretl_matrix_reuse(J->Tmprp, J->r * J->p, 1);

    /* vectorize alpha' into Tmprp */
    vec_transpose(J->Tmprp, J->alpha);

    gretl_matrix_multiply(GGG, J->Tmprp, J->psi);
    gretl_matrix_reuse(J->Tmprp, J->r, J->p);

 bailout:

    gretl_matrix_free(GG);
    gretl_matrix_free(GGG);

    return err;
}

typedef struct switcher_ switcher;

struct switcher_ {
    gretl_matrix *K1;     /* holds kronecker product */
    gretl_matrix *K2;     /* holds kronecker product */
    gretl_matrix *TmpL;   /* used only in \phi calculation */
    gretl_matrix *TmpR;   /* shared temp */
    gretl_matrix *TmpR1;  /* used when \alpha is restricted */
    gretl_matrix *Tmprp1; /* r x p1 temp */
    gretl_matrix *HK2;    /* used in \phi calculation */
};

static void switcher_free (switcher *s)
{
    gretl_matrix_free(s->K1);
    gretl_matrix_free(s->K2);
    gretl_matrix_free(s->TmpL);
    gretl_matrix_free(s->TmpR);
    gretl_matrix_free(s->TmpR1);
    gretl_matrix_free(s->Tmprp1);
    gretl_matrix_free(s->HK2);
}

static int make_lsPi (Jwrap *J)
{
    gretl_matrix *S11i;

    if (J->lsPi != NULL) {
	/* already done */
	return 0;
    }

    S11i = gretl_matrix_copy(J->S11);
    if (S11i == NULL) {
	return E_ALLOC;
    }

    J->lsPi = gretl_matrix_alloc(J->p1, J->p);
    if (J->lsPi == NULL) {
	gretl_matrix_free(S11i);
	return E_ALLOC;
    }

    gretl_invert_symmetric_matrix(S11i);
    gretl_matrix_multiply_mod(S11i, GRETL_MOD_NONE,
			      J->S01, GRETL_MOD_TRANSPOSE,
			      J->lsPi, GRETL_MOD_NONE);
    /* make into vec(\Pi'_{LS}) */
    gretl_matrix_reuse(J->lsPi, J->p1 * J->p, 1);
    gretl_matrix_free(S11i);

    return 0;
}

static int switcher_init (switcher *s, Jwrap *J)
{
    int err = 0;

    s->K1 = NULL; 
    s->K2 = NULL;
    s->TmpL = NULL;
    s->TmpR = NULL;
    s->TmpR1 = NULL;
    s->Tmprp1 = NULL;
    s->HK2 = NULL;

    clear_gretl_matrix_err();

    s->K1 = gretl_matrix_alloc(J->p1 * J->r, J->p1 * J->r);
    s->K2 = gretl_matrix_alloc(J->p1 * J->r, J->p1 * J->p1);
    s->TmpR = gretl_matrix_alloc(J->p1 * J->p1, 1);
    s->Tmprp1 = gretl_matrix_alloc(J->r, J->p1);

    if (J->blen > 0) {
	s->TmpL = gretl_matrix_alloc(J->blen, J->p * J->p1);
	s->HK2 = gretl_matrix_alloc(J->blen, J->p * J->p1);
    }

    err = get_gretl_matrix_err();
    if (err) {
	return err;
    }

    if (J->G != NULL) {
	s->TmpR1 = gretl_matrix_alloc(J->alen, J->p * J->p1);
	if (s->TmpR1 == NULL) {
	    return E_ALLOC;
	}
    }

    err = make_lsPi(J);

    if (!err) {
	err = allocate_info_blocks_etc(J);
    }

    return err;
}

/* The following functions implement the switching algorithm 
   as set out in Boswijk and Doornik, 2004, p. 455.
*/

static void alpha_from_psi (Jwrap *J)
{
    int i, j, k = 0;

    if (J->G != NULL) {
	gretl_matrix_reuse(J->Tmprp, J->p * J->r, 1);
	gretl_matrix_multiply(J->G, J->psi, J->Tmprp);
	for (i=0; i<J->p; i++) {
	    for (j=0; j<J->r; j++) {
		gretl_matrix_set(J->alpha, i, j, J->Tmprp->val[k++]);
	    }
	}
	gretl_matrix_reuse(J->Tmprp, J->r, J->p);
    } else {
	for (i=0; i<J->p; i++) {
	    for (j=0; j<J->r; j++) {
		gretl_matrix_set(J->alpha, i, j, J->psi->val[k++]);
	    }
	}
    }
}

/* 
   Update \psi using:

   [ G'(\Omega^{-1} \otimes \beta'S_{11}\beta)G ]^{-1} 
   \times G'(\Omega^{-1} \otimes \beta'S_{11}) vec(Pi'_{LS})

   or, in case \alpha is not restricted, just compute the
   ML alpha conditional on \beta in the usual Johansen
   manner.
*/

static int update_psi (Jwrap *J, switcher *s)
{
    int err = 0;

    if (J->G == NULL) {
	J_compute_alpha(J);
	return 0;
    }

    gretl_matrix_reuse(s->K1, J->p * J->r, J->p * J->r);
    gretl_matrix_reuse(s->K2, J->p * J->r, J->p * J->p1);
    gretl_matrix_reuse(s->TmpR, J->alen, 1);

    /* left-hand chunk */
    gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE, J->S11,
		       J->qf1, GRETL_MOD_NONE);
    gretl_matrix_kronecker_product(J->iOmega, J->qf1, s->K1);
    gretl_matrix_qform(J->G, GRETL_MOD_TRANSPOSE, s->K1,
		       J->I00, GRETL_MOD_NONE);

    err = gretl_invert_symmetric_matrix(J->I00);
    if (err) {
	return err;
    }

    /* right-hand chunk */
    gretl_matrix_multiply_mod(J->beta, GRETL_MOD_TRANSPOSE, 
			      J->S11, GRETL_MOD_NONE,
			      s->Tmprp1, GRETL_MOD_NONE);
    gretl_matrix_kronecker_product(J->iOmega, s->Tmprp1, s->K2);
    gretl_matrix_multiply_mod(J->G, GRETL_MOD_TRANSPOSE,
			      s->K2, GRETL_MOD_NONE,
			      s->TmpR1, GRETL_MOD_NONE);
    gretl_matrix_multiply(s->TmpR1, J->lsPi, s->TmpR);

    /* combine */
    gretl_matrix_multiply(J->I00, s->TmpR, J->psi);

    /* update alpha */
    alpha_from_psi(J);

    return err;
}

/* update \beta based on \phi vector */

static void beta_from_phi (Jwrap *J)
{
    if (J->blen == 0) {
	return;
    }

    if (J->H != NULL) {
	gretl_matrix_reuse(J->beta, J->p1 * J->r, 1);
	gretl_matrix_multiply(J->H, J->phi, J->beta);
	if (!gretl_is_zero_matrix(J->h0)) {
	    gretl_matrix_add_to(J->beta, J->h0);
	}
	gretl_matrix_reuse(J->beta, J->p1, J->r);
    } else {
	gretl_matrix_copy_values_shaped(J->beta, J->phi);
    }
}

/*
    Update \phi using:

    [H'(\alpha'\Omega^{-1}\alpha \otimes S_{11})H]^{-1} \times
       H'(\alpha'\Omega^{-1} \otimes S_{11}) \times
         [vec(Pi'_{LS}) - (\alpha \otimes I_{p1})h_0]
 */

static int update_phi (Jwrap *J, switcher *s)
{
    int err = 0;

    if (J->blen == 0) {
	return 0;
    }

    gretl_matrix_reuse(s->K1, J->r * J->p1, J->r * J->p1);
    gretl_matrix_reuse(s->K2, J->r * J->p1, J->p * J->p1);
    gretl_matrix_reuse(s->TmpR, J->p * J->p1, 1);

    /* first big inverse */
    gretl_matrix_qform(J->alpha, GRETL_MOD_TRANSPOSE, J->iOmega,
		       J->qf1, GRETL_MOD_NONE);
    gretl_matrix_kronecker_product(J->qf1, J->S11, s->K1);
    if (J->H != NULL) {
	gretl_matrix_qform(J->H, GRETL_MOD_TRANSPOSE, s->K1,
			   J->I11, GRETL_MOD_NONE);
    } else {
	gretl_matrix_copy_values(J->I11, s->K1);
    }

    err = gretl_invert_symmetric_matrix(J->I11);
    if (err) {
	return err;
    }

    /* second chunk */
    gretl_matrix_multiply_mod(J->alpha, GRETL_MOD_TRANSPOSE,
			      J->iOmega, GRETL_MOD_NONE,
			      J->Tmprp, GRETL_MOD_NONE);
    gretl_matrix_kronecker_product(J->Tmprp, J->S11, s->K2);
    if (J->H != NULL) {
	gretl_matrix_multiply_mod(J->H, GRETL_MOD_TRANSPOSE,
				  s->K2, GRETL_MOD_NONE,
				  s->HK2, GRETL_MOD_NONE);
    } else {
	gretl_matrix_copy_values(s->HK2, s->K2);
    }

    /* combine first and second chunks */
    gretl_matrix_multiply(J->I11, s->HK2, s->TmpL);

    /* right-hand chunk */
    gretl_matrix_copy_values(s->TmpR, J->lsPi);
    if (J->h0 != NULL && !gretl_is_zero_matrix(J->h0)) {
	gretl_matrix_reuse(s->K2, J->p * J->p1, J->r * J->p1);
	gretl_matrix_kronecker_I(J->alpha, J->p1, s->K2);
	gretl_matrix_multiply_mod(s->K2, GRETL_MOD_NONE,
				  J->h0, GRETL_MOD_NONE,
				  s->TmpR, GRETL_MOD_DECUMULATE);
    }

    /* combine */
    gretl_matrix_multiply(s->TmpL, s->TmpR, J->phi);

    /* update beta */
    beta_from_phi(J);

    return err;
}

enum {
    OMEGA_ONLY,
    OMEGA_PLUS
};

/* 
    Update Omega using:

    S_{00} - S_{01} \beta\alpha' - \alpha\beta' S_{10} +
       \alpha\beta' S_{11} \beta\alpha'

    then invert into iOmega.
 */

static int make_Omega (Jwrap *J, int code)
{
    int err = 0;

    gretl_matrix_copy_values(J->Omega, J->S00);

    gretl_matrix_multiply_mod(J->alpha, GRETL_MOD_NONE,
			      J->beta, GRETL_MOD_TRANSPOSE,
			      J->Pi, GRETL_MOD_NONE);

    gretl_matrix_multiply_mod(J->S01, GRETL_MOD_NONE,
			      J->Pi, GRETL_MOD_TRANSPOSE,
			      J->Tmppp, GRETL_MOD_NONE);

    gretl_matrix_add_self_transpose(J->Tmppp);
    gretl_matrix_subtract_from(J->Omega, J->Tmppp);

    gretl_matrix_qform(J->Pi, GRETL_MOD_NONE, J->S11,
		       J->Omega, GRETL_MOD_CUMULATE);

    if (code == OMEGA_PLUS) {
	gretl_matrix_copy_values(J->iOmega, J->Omega);
	err = gretl_invert_symmetric_matrix(J->iOmega);
    }

    return err;
}

/* reduced version of log-likelihood calculation for
   switching algorithm loop */

static int switcher_ll (Jwrap *J)
{
    int err = 0;

    gretl_matrix_copy_values(J->Tmppp, J->Omega);
    J->ll = gretl_matrix_log_determinant(J->Tmppp, &err);
    if (!err) {
	J->ll *= -J->T / 2.0;
    }
    
    return err;
}

static int form_I00 (Jwrap *J)
{
    gretl_matrix *K = NULL;
    int err = 0;

    J->I00 = gretl_matrix_alloc(J->alen, J->alen);
    if (J->I00 == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE, J->S11,
		       J->qf1, GRETL_MOD_NONE);

    K = gretl_matrix_kronecker_product_new(J->iOmega, J->qf1, &err);

    if (!err) {
	if (J->G != NULL) {
	    gretl_matrix_qform(J->G, GRETL_MOD_TRANSPOSE, K,
			       J->I00, GRETL_MOD_NONE);
	} else {
	    gretl_matrix_copy_values(J->I00, K);
	}
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(J->I00);
    }

    gretl_matrix_free(K);

    return err;
}

static int form_I11 (Jwrap *J)
{
    gretl_matrix *K = NULL;
    int err = 0;

    if (J->blen == 0) {
	return 0;
    }

    J->I11 = gretl_matrix_alloc(J->blen, J->blen);
    if (J->I11 == NULL) {
	return E_ALLOC;
    }    

    gretl_matrix_qform(J->alpha, GRETL_MOD_TRANSPOSE, J->iOmega,
		       J->qf1, GRETL_MOD_NONE);

    K = gretl_matrix_kronecker_product_new(J->qf1, J->S11, &err);

    if (!err) {
	if (J->H != NULL) {
	    gretl_matrix_qform(J->H, GRETL_MOD_TRANSPOSE, K,
			       J->I11, GRETL_MOD_NONE);
	} else {
	    gretl_matrix_copy_values(J->I11, K);
	}
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(J->I11);
    }

    gretl_matrix_free(K);

    return err;
}

static int add_Omega_inverse (Jwrap *J)
{
    int err = 0;

    J->iOmega = gretl_matrix_copy(J->Omega);

    if (J->iOmega == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_invert_symmetric_matrix(J->iOmega);
    }

    return err;
}

/* Use the information matrix to compute standard errors for
   beta and alpha.
*/

static int variance_from_info_matrix (Jwrap *J)
{
    int err = 0;

    if (J->jr < J->alen + J->blen) {
	/* model is not fully identified */
	return 0;
    }

    /* preliminary */

    if (J->iOmega == NULL) {
	err = add_Omega_inverse(J);
	if (err) return err;
    }

    /* variance of beta */

    if (J->blen == 0) {
	int nb = J->p1 * J->r;

	J->Vb = gretl_zero_matrix_new(nb, nb);
	J->bse = gretl_zero_matrix_new(J->p1, J->r);

	if (J->Vb == NULL || J->bse == NULL) {
	    err = E_ALLOC;
	}

	goto beta_done;
    }	

    if (J->I11 == NULL) {
	err = form_I11(J);
	if (err) return err;
    }

    gretl_matrix_divide_by_scalar(J->I11, J->T);

    if (J->H != NULL) {
	int nb = J->r * J->p1;

	J->Vb = gretl_matrix_alloc(nb, nb);
	if (J->Vb == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_matrix_qform(J->H, GRETL_MOD_NONE, J->I11,
			       J->Vb, GRETL_MOD_NONE);
	}
    } else {
	J->Vb = J->I11;
	J->I11 = NULL;
    }

    if (!err) {
	err = make_beta_se(J);
    }

 beta_done:

    if (err) {
	return err;
    }

    /* variance of alpha */

    if (J->I00 == NULL) {
	err = form_I00(J);
	if (err) return err;
    }    

    gretl_matrix_divide_by_scalar(J->I00, J->T);

    if (J->G != NULL) {
	int na = J->r * J->p;

	J->Va = gretl_matrix_alloc(na, na);
	if (J->Va == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_matrix_qform(J->G, GRETL_MOD_NONE, J->I00,
			       J->Va, GRETL_MOD_NONE);
	}
    } else {
	J->Va = J->I00;
	J->I00 = NULL;
    }

    if (!err) {
	err = make_alpha_se(J);
    }  

#if JDEBUG 
    gretl_matrix_print(J->bse, "J->bse");
    gretl_matrix_print(J->ase, "J->ase");
#endif

    return err;
}

static int user_switch_init (Jwrap *J, int *uinit)
{
    const gretl_matrix *U;
    int psilen;
    int ulen, i, k;
    int err = 0;

    U = get_init_vals();
    if (U == NULL) {
	return 0;
    }

    psilen = (J->G == NULL)? 0 : J->alen;

    ulen = gretl_vector_get_length(U);

    if (ulen == J->blen + psilen) {
	/* we're given phi and psi? */
	k = 0;
	for (i=0; i<J->blen; i++) {
	    J->phi->val[i] = U->val[k++];
	}
	if (psilen == 0) {
	    beta_from_phi(J);
	    J_compute_alpha(J);
	} else {
	    for (i=0; i<psilen; i++) {
		J->psi->val[i] = U->val[k++];
	    }
	}
    } else if (ulen == J->blen + J->alen) {
	/* we're given phi and (redundant) psi? */
	k = 0;
	for (i=0; i<J->blen; i++) {
	    J->phi->val[i] = U->val[k++];
	}
	for (i=0; i<J->alen; i++) {
	    J->psi->val[i] = U->val[k++];
	}
	alpha_from_psi(J);
    } else if (ulen == (J->p1 + J->p) * J->r) {
	/* we're given vec(beta), vec(alpha)? */
	int n = J->p1 * J->r;

	k = 0;
	for (i=0; i<n; i++) {
	    J->beta->val[i] = U->val[k++];
	}

	n = J->p * J->r;
	for (i=0; i<n; i++) {
	    J->alpha->val[i] = U->val[k++];
	}
	
	err = phi_from_beta(J);
	if (!err) {
	    err = psi_from_alpha(J);
	}
    } else {	
	gretl_errmsg_set("Invalid initialization given");
	err = E_DATA;
    } 

    *uinit = 1;

    return err;
}

/* driver function for switching algorithm */

static int switchit (Jwrap *J, PRN *prn)
{
    switcher s;
    double lldiff = NADBL;
    double llbak = -1.0e+200;
    double tol = 4.0e-11;
    int j, jmax = 50000;
    int uinit = 0;
    int conv = 0;
    int err;

    err = user_switch_init(J, &uinit);
    if (err) {
	return err;
    }

    err = switcher_init(&s, J);

    if (!err && J->H != NULL) {
	beta_from_phi(J);
    }

    if (!err && J->G != NULL) {
	alpha_from_psi(J);
    }

#if JDEBUG
    gretl_matrix_print(J->beta, "switchit: initial beta"); 
    gretl_matrix_print(J->alpha, "switchit: initial alpha");
#endif

    if (!err) {
	/* initialize Omega */
	err = make_Omega(J, OMEGA_PLUS);
    }

    for (j=0; j<jmax && !err; j++) {
	err = update_phi(J, &s);
	if (!err) {
	    err = update_psi(J, &s);
	}
	if (!err) {
	    err = make_Omega(J, OMEGA_PLUS);
	}
	if (!err) {
	    err = switcher_ll(J);
	}
	if (!err && j > 1) {
	    lldiff = (J->ll - llbak) / fabs(llbak);
	    if (lldiff < tol) {
		conv = 1;
		break;
	    }
	}
	llbak = J->ll;
    }

    pprintf(prn, "Switching algorithm: %d iterations", j);
    if (uinit) {
	pputs(prn, " (user-supplied initial values)");
    }
    pputc(prn, '\n');
    pprintf(prn, " -(T/2)log|Omega| = %.8g, lldiff = %g\n\n", J->ll, lldiff);

    if (!err && !conv) {
	err = E_NOCONV;
    }

    if (!err) {
	J->ll -= J->llk;
    }

    switcher_free(&s);

    return err;
}

/* 
   J = [(I_p \otimes \beta)G : (\alpha \otimes I_{p1})H]

   Boswijk and Doornik (2004), equation (40), page 455.
*/

static int check_jacobian (Jwrap *J)
{
    gretl_matrix *A = NULL;
    gretl_matrix *B = NULL;
    gretl_matrix *Jac = NULL;
    int err = 0;

    /* form both beta = H \phi + s, and alpha, for randomized 
       \theta
    */

    if (J->H != NULL) {
	gretl_matrix *phi = gretl_column_vector_alloc(J->H->cols);

	gretl_matrix_random_fill(phi, D_NORMAL);
	gretl_matrix_reuse(J->beta, J->p1 * J->r, 1);
	gretl_matrix_multiply(J->H, phi, J->beta);
	gretl_matrix_add_to(J->beta, J->h0);
	gretl_matrix_reuse(J->beta, J->p1, J->r);
	gretl_matrix_free(phi);
    } else {
	gretl_matrix_random_fill(J->beta, D_NORMAL);
    }

    if (J->G != NULL) {
	gretl_matrix *psi = gretl_column_vector_alloc(J->G->cols);

	gretl_matrix_random_fill(psi, D_NORMAL);
	gretl_matrix_reuse(J->alpha, J->p * J->r, 1);
	gretl_matrix_multiply(J->G, psi, J->alpha);
	gretl_matrix_reuse(J->alpha, J->p, J->r);
	gretl_matrix_free(psi);
    } else {	
	J_compute_alpha(J);
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
	Jac = gretl_matrix_col_concat(A, B, &err);
    }

    if (!err) {
	J->jr = gretl_matrix_rank(Jac, &err);
#if JDEBUG
	gretl_matrix_print(Jac, "Jacobian");
#endif
    }

    gretl_matrix_free(A);
    gretl_matrix_free(B);
    gretl_matrix_free(Jac);

    return err;
}

/* Doornik's approach to checking identification */

static int 
vecm_id_check (Jwrap *J, GRETL_VAR *jvar, PRN *prn)
{
    int npar = J->blen + J->alen;
    int err;

    err = check_jacobian(J);

    if (!err) {
	pprintf(prn, "Rank of Jacobian = %d, number of free "
		"parameters = %d\n", J->jr, npar);
	if (J->jr < npar) {
	    pputs(prn, "Model is not fully identified\n");
	    if (J->bfgs) {
		err = E_NOIDENT;
	    }
	} else {
	    pputs(prn, "Model is fully identified\n");
	}

	J->df = (J->p + J->p1 - J->r) * J->r - J->jr;
	pprintf(prn, "Based on Jacobian, df = %d\n", J->df);

	/* system was subject to a prior restriction? */
	if (jvar->jinfo->lrdf > 0) {
	    J->df -= jvar->jinfo->lrdf;
	    pprintf(prn, "Allowing for prior restriction, df = %d\n", 
		    J->df);
	}
    }

    return err;
}

/* common alpha restriction needs to be replicated */

static int G_from_expanded_R (Jwrap *J, const gretl_matrix *R)
{
    gretl_matrix *Rtmp = NULL;
    int err = 0;

    Rtmp = gretl_matrix_I_kronecker_new(J->r, R, &err);

    if (!err) {
	J->G = gretl_matrix_right_nullspace(Rtmp, &err);
	gretl_matrix_free(Rtmp);
    }

    return err;
}

/* set up restriction matrix, G, for alpha */

static int set_up_G (Jwrap *J, const gretl_restriction *rset,
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

    if (J->r > 1 && Ra->cols == J->p) {
	/* got a common alpha restriction */
	err = G_from_expanded_R(J, Ra);
    } else {
	J->G = gretl_matrix_right_nullspace(Ra, &err);
    }

    if (!err) {
	J->alen = J->G->cols;
    }

#if JDEBUG
    gretl_matrix_print(J->G, "G, in set_up_G");
#endif

    return err;
}

/* set up restriction matrices, H and h_0, for beta */

static int set_up_H_h0 (Jwrap *J, const gretl_restriction *rset)
{
    const gretl_matrix *R = rset_get_R_matrix(rset);
    const gretl_matrix *q = rset_get_q_matrix(rset);
    gretl_matrix *RRT = NULL;
    gretl_matrix *Tmp = NULL;
    int err = 0;

#if JDEBUG
    gretl_matrix_print(R, "R, in set_up_H_h0");
    gretl_matrix_print(q, "q, in set_up_H_h0");
#endif

    if (R->rows == J->p1 * J->r) {
	/* number of restrictions = total betas */
	return solve_for_beta(J, R, q);
    }

    if (J->r > 1 && R->cols == J->p1) {
	/* common beta restriction */
	gretl_matrix *Rtmp = gretl_matrix_I_kronecker_new(J->r, R, &err);

	if (!err) {
	    J->H = gretl_matrix_right_nullspace(Rtmp, &err);
	    gretl_matrix_free(Rtmp);
	}
    } else {
	J->H = gretl_matrix_right_nullspace(R, &err);
    }

    if (err) {
	return err;
    }

    J->blen = J->H->cols;

    if (q == NULL || gretl_is_zero_matrix(q)) {
	/* implicitly or explicitly zero */
	J->h0 = gretl_zero_matrix_new(R->cols, 1);
	return (J->h0 == NULL)? E_ALLOC : 0; 
    }

    /* now for h_0, for non-zero q */

    RRT = gretl_matrix_alloc(R->rows, R->rows);
    Tmp = gretl_matrix_alloc(R->cols, R->rows);
    if (RRT == NULL || Tmp == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(R, GRETL_MOD_NONE,
					R, GRETL_MOD_TRANSPOSE,
					RRT, GRETL_MOD_NONE);
    }
    
    if (!err) {
	err = gretl_invert_symmetric_matrix(RRT);
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(R, GRETL_MOD_TRANSPOSE,
					RRT, GRETL_MOD_NONE,
					Tmp, GRETL_MOD_NONE);
    }

    if (!err) {
	J->h0 = gretl_matrix_multiply_new(Tmp, q, &err);
    }

#if JDEBUG
    gretl_matrix_print(J->H, "H, in set_up_H_h0");
    gretl_matrix_print(J->h0, "h_0, in set_up_H_h0");
#endif

    gretl_matrix_free(RRT);
    gretl_matrix_free(Tmp);

    return err;
}

/* See H. Peter Boswijk, "Identifiability of Cointegrated Systems",
   http://www.ase.uva.nl/pp/bin/258fulltext.pdf
   We're assuming here that vec(\beta) = H\psi + h_0, for
   non-zero h_0.
*/

static int phi_init_nonhomog (Jwrap *J)
{
    gretl_matrix *BB = NULL;
    gretl_matrix *BP = NULL;
    gretl_matrix *IBP = NULL;
    gretl_matrix *IBPH = NULL;
    gretl_matrix *IBPh = NULL;
    gretl_matrix *E = NULL;
    int ocols = J->p1 - J->r;
    int err = 0;

    if (J->h0 == NULL || 
	gretl_is_zero_matrix(J->h0) ||
	ocols == 0 || J->blen == 0) {
	/* shouldn't be here! */
	return 0;
    }

    BB = gretl_matrix_alloc(J->p1, J->p1);
    BP = gretl_matrix_alloc(J->p1, ocols);
    IBPH = gretl_matrix_alloc(J->r * ocols, J->blen);
    IBPh = gretl_column_vector_alloc(J->r * ocols);
    
    if (BB == NULL || BP == NULL || 
	IBPH == NULL || IBPh == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    gretl_matrix_multiply_mod(J->beta, GRETL_MOD_NONE,
			      J->beta, GRETL_MOD_TRANSPOSE,
			      BB, GRETL_MOD_NONE);

    E = gretl_symmetric_matrix_eigenvals(BB, 1, &err);

    if (!err) {
	err = gretl_symmetric_eigen_sort(E, BB, 0);
    } 

    if (!err) {
	err = gretl_matrix_extract_matrix(BP, BB, 0, J->r, GRETL_MOD_NONE);
    }

    if (!err) {
	IBP = gretl_matrix_I_kronecker_new(J->r, BP, &err);
    }

    if (!err) {
	gretl_matrix_multiply_mod(IBP, GRETL_MOD_TRANSPOSE,
				  J->H, GRETL_MOD_NONE,
				  IBPH, GRETL_MOD_NONE);
	gretl_matrix_multiply_mod(IBP, GRETL_MOD_TRANSPOSE,
				  J->h0, GRETL_MOD_NONE,
				  IBPh, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_matrix_moore_penrose(IBPH);
#if JDEBUG
	fprintf(stderr, "moore-penrose on IBPH, err = %d\n", err);
	gretl_matrix_print(IBPH, "IBPH (pseudo)inverse");
#endif
    }
    
    if (!err) {
	gretl_matrix_multiply(IBPH, IBPh, J->phi);
	gretl_matrix_switch_sign(J->phi);
    }

 bailout:

    gretl_matrix_free(BB);
    gretl_matrix_free(BP);
    gretl_matrix_free(IBP);
    gretl_matrix_free(IBPH);
    gretl_matrix_free(IBPh);
    gretl_matrix_free(E);

    return err;
}

/* initialization of \phi in case the restriction on beta
   is homogeneous */

static int phi_init_homog (Jwrap *J) 
{
    gretl_matrix *HH = NULL;
    gretl_matrix *Hb = NULL;
    gretl_matrix *b = NULL;
    int nbeta = J->p1 * J->r;
    int i, err = 0;

    b = gretl_matrix_copy(J->beta);
    HH = gretl_matrix_alloc(J->blen, J->blen);
    Hb = gretl_column_vector_alloc(J->blen);

    if (b == NULL || HH == NULL || Hb == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    gretl_matrix_reuse(b, nbeta, 1);

    err = gretl_matrix_multiply_mod(J->H, GRETL_MOD_TRANSPOSE,
				    J->H, GRETL_MOD_NONE,
				    HH, GRETL_MOD_NONE);
    
    if (!err) {
	err = gretl_invert_symmetric_matrix(HH);
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(J->H, GRETL_MOD_TRANSPOSE,
					b, GRETL_MOD_NONE,
					Hb, GRETL_MOD_NONE);
    }

    if (!err) {
	gretl_matrix_reuse(b, Hb->rows, 1);
	err = gretl_matrix_multiply(HH, Hb, b);
    }

    if (!err) {
	for (i=0; i<b->rows; i++) {
	    J->phi->val[i] = b->val[i];
	}
    } 

 bailout:

    gretl_matrix_free(HH);
    gretl_matrix_free(Hb);
    gretl_matrix_free(b);

    return err;
}

static int phi_from_beta (Jwrap *J)
{
    int err = 0;

    if (J->H == NULL) {
	/* just vectorize beta into phi */
	vec_simple(J->phi, J->beta);
    } else if (gretl_is_zero_matrix(J->h0)) {
	err = phi_init_homog(J);
    } else {
	/* solve for \phi a la Boswijk */
	err = phi_init_nonhomog(J);
    }

    return err;
}

/* solution to unrestricted eigenvalue problem */

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
	err = gretl_symmetric_eigen_sort(evals, M, J->r);
    }

#if JDEBUG
    gretl_matrix_print(M, "case0: raw beta");
#endif

    if (!err) {
	gretl_matrix_copy_values(J->beta, M);
    }

    gretl_matrix_free(Tmp);
    gretl_matrix_free(M);
    gretl_matrix_free(evals);

    return err;
}

/* 
   Initialize beta, or the free parameters associated with
   beta in case beta is restricted. 
*/

static int beta_init (Jwrap *J)
{
    int err;

    err = case0(J);

    if (!err) {
	err = phi_from_beta(J);
    }

    if (!err && J->H != NULL) {
	beta_from_phi(J);
    }

#if JDEBUG
    gretl_matrix_print(J->phi, "beta_init: final phi vector");
#endif
    
    return err;
}

/* do initial computation of alpha based on beta; then, 
   if relevant, impose alpha restriction and make phi
*/

static int alpha_init (Jwrap *J)
{
    int err = 0;

    err = J_compute_alpha(J);

#if JDEBUG
    gretl_matrix_print(J->alpha, "alpha_init: result from compute_alpha");
    gretl_matrix_print(J->beta, "(using this beta)");
#endif

    if (!err) {
	err = psi_from_alpha(J);
    }

#if JDEBUG
    gretl_matrix_print(J->psi, "alpha_init: final psi");
#endif

    return err;
}

static int make_beta_se (Jwrap *J)
{
    double x;
    int i;

    J->bse = gretl_matrix_alloc(J->p1, J->r);
    if (J->bse == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<J->Vb->rows; i++) {
	x = gretl_matrix_get(J->Vb, i, i);
	J->bse->val[i] = sqrt(x);
    }

    return 0;
}

static int make_alpha_se (Jwrap *J)
{
    double x;
    int i, j, k;

    J->ase = gretl_matrix_alloc(J->p, J->r);
    if (J->ase == NULL) {
	return E_ALLOC;
    }

    k = 0;
    for (i=0; i<J->ase->rows; i++) {
	for (j=0; j<J->ase->cols; j++) {
	    x = gretl_matrix_get(J->Va, k, k);
	    gretl_matrix_set(J->ase, i, j, sqrt(x));
	    k++;
	}
    }

    return 0;
}

/* called by BFGS/simann callback (only) */

static int real_compute_ll (Jwrap *J)
{
    int err = 0;

    if (J->G == NULL) {
	err = J_compute_alpha(J);
    } else {
	alpha_from_psi(J);
    }

    if (!err) {
	err = make_Omega(J, OMEGA_ONLY);
    }

    if (!err) {
	gretl_matrix_copy_values(J->Tmppp, J->Omega);
	J->ll = gretl_matrix_log_determinant(J->Tmppp, &err);
    }

    if (!err) {
	J->ll *= -J->T * 0.5;
	J->ll -= J->llk;
    }    

    return err;
}

static void 
sync_with_theta (Jwrap *J, const double *theta)
{
    int i, k = 0;

    if (J->H != NULL) {
	for (i=0; i<J->blen; i++) {
	    J->phi->val[i] = theta[k++];
	}
    } 

    if (J->G != NULL) {
	for (i=0; i<J->alen; i++) {
	    J->psi->val[i] = theta[k++];
	}
    }
}

/* BFGS (and simann) callback function */

static double Jloglik (const double *theta, void *data)
{
    Jwrap *J = (Jwrap *) data;

    sync_with_theta(J, theta);

    if (J->blen > 0) {
	beta_from_phi(J);
    }

    real_compute_ll(J);

    return J->ll;
}

struct gradhelper_ {
    gretl_matrix *vPiT;
    gretl_matrix *dPi;
    gretl_matrix *BS;
    gretl_matrix *K;
    gretl_matrix *LK;
    gretl_matrix *psigrad;
    gretl_matrix *phigrad;
};

static void gradhelper_free (gradhelper *g)
{
    if (g == NULL) {
	return;
    }

    gretl_matrix_free(g->vPiT);
    gretl_matrix_free(g->dPi);
    gretl_matrix_free(g->BS);
    gretl_matrix_free(g->K);
    gretl_matrix_free(g->LK);
    gretl_matrix_free(g->psigrad);
    gretl_matrix_free(g->phigrad);

    free(g);
};    

static int add_gradhelper (Jwrap *J)
{
    gradhelper *g = malloc(sizeof *g);
    int LKr, err;

    if (g == NULL) {
	return E_ALLOC;
    }

    err = make_lsPi(J);
    if (err) {
	free(g);
	return err;
    }

    J->iOmega = gretl_matrix_alloc(J->p, J->p);
    if (J->iOmega == NULL) {
	free(g);
	return E_ALLOC;
    }    

    g->vPiT = NULL;
    g->dPi = NULL;
    g->K = NULL;
    g->LK = NULL;
    g->BS = NULL;
    g->psigrad = NULL;
    g->phigrad = NULL;

    if (J->blen == 0) {
	LKr = J->alen;
    } else if (J->G == NULL) {
	LKr = J->blen;
    } else if (J->H == NULL) {
	LKr = J->alen;
    } else {
	LKr = (J->blen > J->alen)? J->blen : J->alen;
    }

    clear_gretl_matrix_err();

    g->vPiT = gretl_matrix_alloc(J->p * J->p1, 1);
    g->dPi = gretl_matrix_alloc(J->p * J->p1, 1);
    g->K = gretl_matrix_alloc(J->p1 * J->r, J->p * J->p1);
    g->LK = gretl_matrix_alloc(LKr, J->p * J->p1);

    if (J->G != NULL) {
	g->BS = gretl_matrix_alloc(J->r, J->p1);
	g->psigrad = gretl_matrix_alloc(J->alen, 1);
    }

    if (J->blen > 0) {
	g->phigrad = gretl_matrix_alloc(J->blen, 1);
    }    

    err = get_gretl_matrix_err();
    if (err) {
	gradhelper_free(g);
    } else {
	J->ghelper = g;
    }

    return err;
}   

/* optional BFGS callback */ 

static int Jgradient (double *theta, double *gr, int n,
		      BFGS_CRIT_FUNC f, void *data)
{
    Jwrap *J = (Jwrap *) data;
    gradhelper *g = J->ghelper;
    int i, k, err = 0;

    sync_with_theta(J, theta);

    if (J->blen > 0) {
	beta_from_phi(J);
    }

    if (J->G == NULL) {
	err = J_compute_alpha(J);
    } else {
	alpha_from_psi(J);
    }

    if (!err) {
	err = make_Omega(J, OMEGA_PLUS);
    }

    /* form vec(\Pi_{LS}' - \Pi') */

    vec_transpose(g->vPiT, J->Pi);
    gretl_matrix_copy_values(g->dPi, J->lsPi);
    gretl_matrix_subtract_from(g->dPi, g->vPiT);

    /* form psi component of gradient, unless alpha is unrestricted */

    if (J->G != NULL) {
	gretl_matrix_multiply_mod(J->beta, GRETL_MOD_TRANSPOSE,
				  J->S11, GRETL_MOD_NONE,
				  g->BS, GRETL_MOD_NONE);
	gretl_matrix_reuse(g->K, J->p * J->r, -1);
	gretl_matrix_kronecker_product(J->iOmega, g->BS, g->K);

	gretl_matrix_reuse(g->LK, J->alen, -1);
	gretl_matrix_multiply_mod(J->G, GRETL_MOD_TRANSPOSE,
				  g->K, GRETL_MOD_NONE,
				  g->LK, GRETL_MOD_NONE);
	gretl_matrix_multiply(g->LK, g->dPi, g->psigrad);
	gretl_matrix_multiply_by_scalar(g->psigrad, J->T);
    }

    /* form phi component, unless beta is fixed */

    if (J->blen > 0) {
	gretl_matrix_multiply_mod(J->alpha, GRETL_MOD_TRANSPOSE,
				  J->iOmega, GRETL_MOD_NONE,
				  J->Tmprp, GRETL_MOD_NONE);
	gretl_matrix_reuse(g->K, J->p1 * J->r, -1);
	gretl_matrix_kronecker_product(J->Tmprp, J->S11, g->K);

	if (J->H != NULL) {
	    gretl_matrix_reuse(g->LK, J->blen, -1);
	    gretl_matrix_multiply_mod(J->H, GRETL_MOD_TRANSPOSE,
				      g->K, GRETL_MOD_NONE,
				      g->LK, GRETL_MOD_NONE);
	    gretl_matrix_multiply(g->LK, g->dPi, g->phigrad);
	} else {
	    gretl_matrix_multiply(g->K, g->dPi, g->phigrad);
	}

	gretl_matrix_multiply_by_scalar(g->phigrad, J->T);
    }

    /* transcribe results to the gr array */
    k = 0;
    for (i=0; i<J->blen; i++) {
	gr[k++] = g->phigrad->val[i];
    }
    if (J->G != NULL) {
	for (i=0; i<J->alen; i++) {
	    gr[k++] = g->psigrad->val[i];
	}
    }

    return err;
}

#define VECM_WIDTH 13

static int printres (Jwrap *J, GRETL_VAR *jvar, const DATAINFO *pdinfo,
		     PRN *prn)
{
    JohansenInfo *jv = jvar->jinfo;
    const gretl_matrix *c = J->beta;
    const gretl_matrix *sd = J->bse;
    char vname[32], s[16];
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
    pputs(prn, _("Cointegrating vectors"));
    if (sd != NULL) {
	pprintf(prn, " (%s)", _("standard errors in parentheses"));
    }
    pputs(prn, "\n\n");

    for (i=0; i<J->p1; i++) {
	if (i < jvar->ylist[0]) {
	    sprintf(vname, "%s(-1)", pdinfo->varname[jvar->ylist[i+1]]);
	} else if (jv->code == J_REST_CONST) {
	    strcpy(vname, "const");
	} else if (jv->code == J_REST_TREND) {
	    strcpy(vname, "trend");
	}
	pprintf(prn, "%-12s", vname);

	for (j=0; j<J->r; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(c, i, j));
	}
	pputc(prn, '\n');

	if (sd != NULL) {
	    bufspace(VECM_WIDTH, prn);
	    for (j=0; j<J->r; j++) {
		sprintf(s, "(%#.5g)", gretl_matrix_get(sd, i, j));
		pprintf(prn, "%12s ", s);
	    }
	    pputc(prn, '\n');
	}
    }

    c = J->alpha;
    sd = J->ase;

    pputc(prn, '\n');
    pputs(prn, _("alpha (adjustment vectors)"));
    if (sd != NULL) {
	pprintf(prn, " (%s)", _("standard errors in parentheses"));
    }
    pputs(prn, "\n\n");

    for (i=0; i<J->p; i++) {
	sprintf(vname, "%s", pdinfo->varname[jvar->ylist[i+1]]);
	pprintf(prn, "%-12s", vname);
	for (j=0; j<J->r; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(c, i, j));
	}
	pputc(prn, '\n');

	if (sd != NULL) {
	    bufspace(VECM_WIDTH, prn);
	    for (j=0; j<J->r; j++) {
		sprintf(s, "(%#.5g)", gretl_matrix_get(sd, i, j));
		pprintf(prn, "%12s ", s);
	    }
	    pputc(prn, '\n');
	}
    }

    pputc(prn, '\n');

    return 0;
}

/* simulated annealing: can help when the standard deterministic
   initialization (theta_0) leads to a local maximum
*/

/* simann method variants --

   1: Use ll-best point from annealing (including theta_0 in the
      evaluation)
   2: Use the last value from annealing, regardless.
   3: hybrid: use ll-best point in case of improvement, else
      last value.
*/

static int simann (Jwrap *J, gretlopt opt, PRN *prn)
{
    gretl_matrix *b = J->theta;
    int i, SAiter = 4096;
    double f0, f1;
    double fbest, fworst;
    double rndu;
    int jump;

    gretl_matrix *b0 = NULL;
    gretl_matrix *b1 = NULL;
    gretl_matrix *bstar = NULL;
    gretl_matrix *d = NULL;

    double Temp = 1.0;
    double radius = 1.0;
    int improved = 0;
    int method = 3; /* hybrid */
    int err = 0;

    b0 = gretl_matrix_copy(b);
    b1 = gretl_matrix_copy(b);
    bstar = gretl_matrix_copy(b);
    d = gretl_column_vector_alloc(b->rows);

    if (b0 == NULL || b1 == NULL || bstar == NULL || d == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    f0 = fbest = fworst = Jloglik(b0->val, J);

    if (opt & OPT_V) {
	pprintf(prn, "\nSimulated annealing: initial function value = %.8g\n",
		f0);
    }    

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
		gretl_matrix_copy_values(bstar, b0);
		if (opt & OPT_V) {
		    if (!improved) {
			pprintf(prn, "\n%6s %12s %12s %12s\n",
				"iter", "temp", "radius", "fbest");
		    }
		    pprintf(prn, "%6d %#12.6g %#12.6g %#12.6g\n", 
			    i, Temp, radius, fbest);
		}
		improved = 1;
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

    if (method == 1) {
	/* force to "all-time best" */
	gretl_matrix_copy_values(J->theta, bstar);
    } else if (method == 2) {
	/* take last value from annealing */
	gretl_matrix_copy_values(J->theta, b0);
    } else {
	/* hybrid method */
	if (improved) {
	    gretl_matrix_copy_values(J->theta, bstar);
	} else {
	    gretl_matrix_copy_values(J->theta, b0);
	}
    }

    sync_with_theta(J, b->val);

    if (improved) {
	if (opt & OPT_V) pputc(prn, '\n');
    } else {
	pprintf(prn, "No improvement found in %d iterations\n\n", SAiter);
    }
    
    if (fbest - fworst < 1.0e-9) {
	pprintf(prn, "Warning: likelihood seems to be flat\n");
    }

 bailout:

    gretl_matrix_free(b0);
    gretl_matrix_free(b1);
    gretl_matrix_free(bstar);
    gretl_matrix_free(d);

    return err;
}

/* solve for unrestricted alpha conditional on beta */

static int J_compute_alpha (Jwrap *J)
{
    gretl_matrix *S01b;
    int err = 0;

    S01b = gretl_matrix_reuse(J->Tmprp, J->p, J->r);
    gretl_matrix_multiply(J->S01, J->beta, S01b);

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

    gretl_matrix_reuse(J->Tmprp, J->r, J->p);

    return err;
}

static int allocate_phi (Jwrap *J)
{
    int n = (J->H != NULL)? J->H->cols : J->p1 * J->r;

    J->phi = gretl_zero_matrix_new(n, 1);

    return (J->phi == NULL)? E_ALLOC : 0;
}

static int allocate_psi (Jwrap *J)
{
    int n = (J->G != NULL)? J->G->cols : J->p * J->r;

    J->psi = gretl_zero_matrix_new(n, 1);

    return (J->psi == NULL)? E_ALLOC : 0;
}

/* Needed only when using BFGS (or when initializing using simulated
   annealing): concatenate phi and psi to form the full vector of free
   parameters, theta.  But note: if alpha is unrestricted, we'll just
   compute alpha conditional on beta, and so we don't append phi to
   the free params in theta.
*/

static int make_theta (Jwrap *J)
{
    int i, k = 0;

    if (J->H != NULL) {
	k += J->blen;
    }

    if (J->G != NULL) {
	k += J->alen;
    }

    if (k == 0) {
	/* can't happen */
	return 0;
    }

    J->theta = gretl_column_vector_alloc(k);
    if (J->theta == NULL) {
	return E_ALLOC;
    }

    k = 0;

    if (J->H != NULL) {
	for (i=0; i<J->blen; i++) {
	    J->theta->val[k++] = J->phi->val[i];
	}
    }

    if (J->G != NULL) {
	for (i=0; i<J->alen; i++) {
	    J->theta->val[k++] = J->psi->val[i];
	}
    }    

    return 0;
}

/* write info from the temporary Jwrap structure into
   the "permanent" GRETL_VAR structure */

static void transcribe_to_jvar (Jwrap *J, GRETL_VAR *jvar)
{
    jvar->jinfo->ll0 = jvar->ll;
    jvar->ll = J->ll;
    jvar->jinfo->lrdf += J->df; /* ? */

    gretl_matrix_free(jvar->S);
    jvar->S = J->Omega;
    J->Omega = NULL;

    gretl_matrix_free(jvar->jinfo->Beta);
    jvar->jinfo->Beta = J->beta;
    J->beta = NULL;

    gretl_matrix_free(jvar->jinfo->Alpha);
    jvar->jinfo->Alpha = J->alpha;
    J->alpha = NULL;

    gretl_matrix_free(jvar->jinfo->Bvar);
    jvar->jinfo->Bvar = J->Vb;
    J->Vb = NULL;

    gretl_matrix_free(jvar->jinfo->Bse);
    jvar->jinfo->Bse = J->bse;
    J->bse = NULL;

    gretl_matrix_free(jvar->jinfo->Ase);
    jvar->jinfo->Ase = J->ase;
    J->ase = NULL;    
}

static int do_bfgs (Jwrap *J, gretlopt opt, PRN *prn)
{
    int maxit = 4000;
    double reltol = 1.0e-11;
    int fncount = 0;
    int grcount = 0;
    int nn = J->theta->rows;
    int err = 0;

    if (J->bfgs == BFGS_ANALYTIC) {
	pputs(prn, "LBFGS: using analytical derivatives\n\n");
	err = add_gradhelper(J);
    } else {
	pputs(prn, "LBFGS: using numerical derivatives\n\n");
    }

    if (!err) {
	err = LBFGS_max(J->theta->val, nn, maxit, reltol, 
			&fncount, &grcount, Jloglik, C_LOGLIK,
			(J->bfgs == BFGS_ANALYTIC)? Jgradient : NULL,
			J, opt, prn);
    }

    return err;
}

/* 
   OPT_L: use LBFGS approach instead of switching algorithm.
   OPT_F: doing full estimation of restricted system, not just 
          testing the restriction.
   OPT_J: ("jitter") use simulated annealing in initialization.
*/

int general_vecm_analysis (GRETL_VAR *jvar, 
			   const gretl_restriction *rset,
			   const DATAINFO *pdinfo,
			   PRN *prn)
{
    Jwrap *J = NULL;
    gretlopt opt = gretl_restriction_get_options(rset);
    int full = (opt & OPT_F);
    int do_simann = (opt & OPT_J);
    int err = 0;

    J = jwrap_new(jvar, opt, &err);
    if (err) {
	return err;
    }

    if (rset_VECM_bcols(rset) > 0) {
	err = set_up_H_h0(J, rset);
    }

    if (!err) {
	err = allocate_phi(J);
    }

    if (!err && rset_VECM_acols(rset) > 0) {
	err = set_up_G(J, rset, prn);
    }

    if (!err) {
	err = allocate_psi(J);
    }   

    if (!err && J->blen == 0 && J->G == NULL) {
	/* beta doesn't need to be estimated */
	err = real_compute_ll(J);
	goto skipest;
    }

    if (!err) {
	err = vecm_id_check(J, jvar, prn);
    }

    if (!err) {
	err = beta_init(J);
    }

    if (!err && (!J->bfgs || J->G != NULL)) {
	err = alpha_init(J);
    }

    if (!err && (J->bfgs || do_simann)) {
	err = make_theta(J);
    }

    if (!err && do_simann) {
	err = simann(J, opt, prn);
    }

    if (!err) {
	if (J->bfgs) {
	    err = do_bfgs(J, opt, prn);
	} else {
	    err = switchit(J, prn);
	}
    }

 skipest:

    if (!err) {
	err = variance_from_info_matrix(J);
    }

    if (!err) {
	if (full) {
	    transcribe_to_jvar(J, jvar);
	} else {
	    printres(J, jvar, pdinfo, prn);
	}
    } 

    jwrap_destroy(J);

    return err;
}
