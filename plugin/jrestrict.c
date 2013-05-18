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

/* General restrictions on beta in context of a VECM.  Based on an ox
   program by Jack, June 2007.  Rendered into C by Allin.
*/

#include "libgretl.h"
#include "random.h"
#include "johansen.h"
#include "var.h"
#include "varprint.h"
#include "gretl_bfgs.h"
#include "gretl_restrict.h"
#include "libset.h"
#include "jprivate.h"

#define JDEBUG 0

typedef struct gradhelper_ gradhelper;
typedef struct Jwrap_ Jwrap;

enum {
    J_USE_LBFGS   = 1 << 0,
    J_CROSS_ALPHA = 1 << 1
};

struct Jwrap_ {
    char flags;     /* status bitflags */ 
    int T;          /* length of time series used */
    int p;          /* number of equations */
    int p1;         /* number of rows in beta (>= p) */
    int r;          /* cointegrating rank of VECM */
    int blen;       /* number of unrestricted coefficients in beta */
    int alen;       /* number of unrestricted coefficients in alpha */
    int df;         /* degrees of freedom for LR test */
    int jr;         /* rank of Jacobian */
    double llk;     /* constant in log-likelihood */
    double ll;      /* log-likelihood */

    /* moment matrices and copies */
    const gretl_matrix *R0;
    const gretl_matrix *R1;
    const gretl_matrix *S00;
    const gretl_matrix *S01;
    const gretl_matrix *S11;

    gretl_matrix *S00i;

    /* restrictions on beta */
    gretl_matrix *H;
    gretl_matrix *h;

    /* restrictions on beta, incl. normalizations */
    gretl_matrix *H0;

    /* augmented beta restriction matrices, if needed */
    gretl_matrix *bigR;
    gretl_matrix *bigq;

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

    /* for keeping track of column scaling */
    int *normrow;
    int *normcol;
    double *normval;
};

#define do_scaling(J) (J->normrow != NULL)
#define using_bfgs(J) (J->flags & J_USE_LBFGS)
#define cross_alpha(J) (J->flags & J_CROSS_ALPHA)

static int J_compute_alpha (Jwrap *J);
static int make_beta_se (Jwrap *J);
static int make_alpha_se (Jwrap *J);
static int phi_from_beta (Jwrap *J);
static void gradhelper_free (gradhelper *g);
static int real_set_up_H (Jwrap *J, const gretl_matrix *R,
			  const gretl_matrix *q);

static int make_moment_matrices (Jwrap *J, const GRETL_VAR *jvar)
{
    gretl_matrix *Tmp = NULL;
    int err = 0;

    J->R0 = jvar->jinfo->R0;
    J->R1 = jvar->jinfo->R1;

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
	/* allocate beta, alpha, etc. while we're at it:
	   but are these sizes always right? */
	J->beta = gretl_zero_matrix_new(J->p1, J->r);
	J->alpha = gretl_zero_matrix_new(J->p, J->r);
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

    if (J->H0 == J->H) {
	J->H0 = NULL;
    }

    gretl_matrix_free(J->H);
    gretl_matrix_free(J->h);
    gretl_matrix_free(J->H0);

    gretl_matrix_free(J->bigR);
    gretl_matrix_free(J->bigq);

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

    free(J->normrow);
    free(J->normcol);
    free(J->normval);

    free(J);
}

static Jwrap *jwrap_new (const GRETL_VAR *jvar, gretlopt opt, int *err)
{
    Jwrap *J = malloc(sizeof *J);

    if (J == NULL) {
	return NULL;
    }

    J->flags = 0;
    if (opt & OPT_L) {
	J->flags |= J_USE_LBFGS;
    }

    J->T = jvar->T;
    J->p = jvar->neqns;
    J->r = jrank(jvar);
    J->blen = 0;
    J->alen = 0;
    J->df = 0;
    J->jr = 0;

    J->llk = J->T * 0.5 * (J->p * (1.0 + LN_2_PI));
    J->ll = NADBL;

    J->R0 = J->R1 = NULL;
    
    J->S00 = NULL;
    J->S01 = NULL;
    J->S11 = NULL;

    J->S00i = NULL;

    J->H = NULL;
    J->h = NULL;
    J->H0 = NULL;

    J->bigR = NULL;
    J->bigq = NULL;

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

    J->normrow = NULL;
    J->normcol = NULL;
    J->normval = NULL;

    *err = make_moment_matrices(J, jvar);

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

/* unvectorize src into targ, transposing */

static void 
unvec_transpose (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, j, k = 0;

    for (i=0; i<targ->rows; i++) {
	for (j=0; j<targ->cols; j++) {
	    gretl_matrix_set(targ, i, j, src->val[k++]);
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
    gretl_matrix *b = gretl_matrix_copy(q);
    int err = 0;

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
    gretl_matrix_block *B; /* holder for the following */
    gretl_matrix *K1;      /* holds kronecker product */
    gretl_matrix *K2;      /* holds kronecker product */
    gretl_matrix *TmpL;    /* used only in \phi calculation */
    gretl_matrix *TmpR;    /* used only in \phi calculation */
    gretl_matrix *TmpR1;   /* used when \alpha is restricted */
    gretl_matrix *Tmprp1;  /* r x p1 temp */
};

static void switcher_free (switcher *s)
{
    gretl_matrix_block_destroy(s->B);

    gretl_matrix_free(s->TmpL);
    gretl_matrix_free(s->TmpR1);
}

static int make_lsPi (Jwrap *J)
{
    gretl_matrix *S11i;
    int err = 0;

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

    err = gretl_invert_symmetric_matrix(S11i);

    if (!err) {
	gretl_matrix_multiply_mod(S11i, GRETL_MOD_NONE,
				  J->S01, GRETL_MOD_TRANSPOSE,
				  J->lsPi, GRETL_MOD_NONE);
	/* make into vec(\Pi'_{LS}) */
	gretl_matrix_reuse(J->lsPi, J->p1 * J->p, 1);
    }

    gretl_matrix_free(S11i);

    return err;
}

static int switcher_init (switcher *s, Jwrap *J)
{
    int err = 0;

    s->TmpL = NULL;
    s->TmpR1 = NULL;

    s->B = gretl_matrix_block_new(&s->K1,   J->p1 * J->r, J->p1 * J->r,
				  &s->K2,   J->p1 * J->r, J->p1 * J->p1,
				  &s->TmpR, J->p * J->p1, 1,
				  &s->Tmprp1, J->r, J->p1,
				  NULL);
    if (s->B == NULL) {
	return E_ALLOC;
    }

    if (J->blen > 0) {
	s->TmpL = gretl_matrix_alloc(J->blen, J->p * J->p1);
	if (s->TmpL == NULL) {
	    return E_ALLOC;
	}
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

    /* left-hand chunk */
    gretl_matrix_qform(J->beta, GRETL_MOD_TRANSPOSE, J->S11,
		       J->qf1, GRETL_MOD_NONE);
    gretl_matrix_kronecker_product(J->iOmega, J->qf1, s->K1);
    gretl_matrix_qform(J->G, GRETL_MOD_TRANSPOSE, s->K1,
		       J->I00, GRETL_MOD_NONE);

    /* right-hand chunk */
    gretl_matrix_multiply_mod(J->beta, GRETL_MOD_TRANSPOSE, 
			      J->S11, GRETL_MOD_NONE,
			      s->Tmprp1, GRETL_MOD_NONE);
    gretl_matrix_kronecker_product(J->iOmega, s->Tmprp1, s->K2);
    gretl_matrix_multiply_mod(J->G, GRETL_MOD_TRANSPOSE,
			      s->K2, GRETL_MOD_NONE,
			      s->TmpR1, GRETL_MOD_NONE);
    gretl_matrix_multiply(s->TmpR1, J->lsPi, J->psi);

    /* combine */
    err = gretl_cholesky_decomp_solve(J->I00, J->psi);

    /* update alpha */
    if (!err) {
	alpha_from_psi(J);
    }

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
	if (!gretl_is_zero_matrix(J->h)) {
	    gretl_matrix_add_to(J->beta, J->h);
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

    /* second chunk */
    gretl_matrix_multiply_mod(J->alpha, GRETL_MOD_TRANSPOSE,
			      J->iOmega, GRETL_MOD_NONE,
			      J->Tmprp, GRETL_MOD_NONE);
    gretl_matrix_kronecker_product(J->Tmprp, J->S11, s->K2);
    if (J->H != NULL) {
	gretl_matrix_multiply_mod(J->H, GRETL_MOD_TRANSPOSE,
				  s->K2, GRETL_MOD_NONE,
				  s->TmpL, GRETL_MOD_NONE);
    } else {
	gretl_matrix_copy_values(s->TmpL, s->K2);
    }

    /* combine first and second chunks */
    err = gretl_cholesky_decomp_solve(J->I11, s->TmpL);
    if (err) {
	fprintf(stderr, "cholesky decomp failed in update_phi\n");
    }

    if (!err) {
	/* right-hand chunk */
	gretl_matrix_copy_values(s->TmpR, J->lsPi);
	if (J->h != NULL && !gretl_is_zero_matrix(J->h)) {
	    gretl_matrix_reuse(s->K2, J->p * J->p1, J->r * J->p1);
	    gretl_matrix_kronecker_I(J->alpha, J->p1, s->K2);
	    gretl_matrix_multiply_mod(s->K2, GRETL_MOD_NONE,
				      J->h, GRETL_MOD_NONE,
				      s->TmpR, GRETL_MOD_DECREMENT);
	}

	/* combine */
	gretl_matrix_multiply(s->TmpL, s->TmpR, J->phi);

	/* update beta */
	beta_from_phi(J);
    }

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

    then invert into iOmega if wanted.
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

static int form_I00_inverse (Jwrap *J)
{
    gretl_matrix *K = NULL;
    int err = 0;

    if (J->I00 == NULL) {
	J->I00 = gretl_matrix_alloc(J->alen, J->alen);
	if (J->I00 == NULL) {
	    return E_ALLOC;
	}
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

static int form_I11_inverse (Jwrap *J)
{
    gretl_matrix *K = NULL;
    int err = 0;

    if (J->blen == 0) {
	return 0;
    }

    if (J->I11 == NULL) {
	J->I11 = gretl_matrix_alloc(J->blen, J->blen);
	if (J->I11 == NULL) {
	    return E_ALLOC;
	} 
    }   

    gretl_matrix_qform(J->alpha, GRETL_MOD_TRANSPOSE, J->iOmega,
		       J->qf1, GRETL_MOD_NONE);

    K = gretl_matrix_kronecker_product_new(J->qf1, J->S11, &err);

    if (!err) {
	const gretl_matrix *H = (do_scaling(J))? J->H0 : J->H;

	if (H != NULL) {
	    gretl_matrix_qform(H, GRETL_MOD_TRANSPOSE, K,
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

static void restricted_vecm_set_df (Jwrap *J, GRETL_VAR *v)
{
    int K = 0;
    double c;

    /* lagged differences */
    K += v->order * J->p;

    /* deterministic stuff */
    K += v->jinfo->seasonals + v->ifc;
    if (jcode(v) == J_UNREST_TREND) {
	K++;
    }

    /* exogenous vars? */
    if (v->xlist != NULL) {
	K += v->xlist[0];
    }

    K *= J->p;

    /* add number of free beta and alpha terms */
    K += J->blen + J->alen;

    c = K / (double) J->p;
    v->df = v->T - floor(c);

#if JDEBUG
    fprintf(stderr, "T = %d, K = %d, c = %g, df = %d\n", v->T, K, c, v->df);
#endif
}

/* Use the information matrix to compute standard errors for
   beta and alpha.
*/

static int variance_from_info_matrix (Jwrap *J, GRETL_VAR *v)
{
    int npar = J->alen + J->blen;
    int err = 0;

    if (J->jr < npar) {
	/* model is not fully identified */
	return 0;
    }

    /* preliminary */

    restricted_vecm_set_df(J, v);

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

    if (do_scaling(J)) {
	gretl_matrix_free(J->I11);
	J->I11 = NULL;
    }

    err = form_I11_inverse(J);
    if (err) return err;

    gretl_matrix_divide_by_scalar(J->I11, v->df);

    if (J->H0 != NULL) {
	int nb = J->r * J->p1;

	J->Vb = gretl_matrix_alloc(nb, nb);
	if (J->Vb == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_matrix_qform(J->H0, GRETL_MOD_NONE, J->I11,
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

    err = form_I00_inverse(J);
    if (err) return err;

    gretl_matrix_divide_by_scalar(J->I00, v->df);

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

/* preliminary test for whether we're using manual
   initialization with the switching algorithm: if so,
   we won't try to extract column scale factors
*/

static int doing_manual_init (Jwrap *J)
{
    const gretl_matrix *U = get_init_vals();
    int ret = 0;

    if (U != NULL) {
	int psilen = (J->G == NULL)? 0 : J->alen;
	int ulen = gretl_vector_get_length(U);

	if (ulen == J->blen + psilen) {
	    ret = 1;
	} else if (ulen == J->blen + J->alen) {
	    ret = 1;
	} else if (ulen == (J->p1 + J->p) * J->r) {
	    ret = 1;
	} 
    }

    return ret;
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

#if JDEBUG > 1

static void print_switch_iter (Jwrap *J, int i)
{
    double wi, wmin = NADBL, wmax = 0;

    fprintf(stderr, "Switcher, iteration %d, ll = %.10g\n", i, J->ll);
    gretl_matrix_print(J->Omega, "J->Omega");
    gretl_matrix_print(J->iOmega, "J->iOmega");
    gretl_matrix_print(J->beta, "J->beta"); 
    gretl_matrix_print(J->alpha, "J->alpha");

    for (i=0; i<J->p; i++) {
	wi = fabs(gretl_matrix_get(J->Omega, i, i));
	if (wi > wmax) {
	    wmax = wi;
	}
	if (wi < wmin) {
	    wmin = wi;
	}
    }

    fprintf(stderr, "Largest/smallest diagonal elements of Omega = "
	    " %g / %g = %g\n", wmax, wmin, wmax / wmin);
}

#endif

/* driver function for switching algorithm */

static int switchit (Jwrap *J, PRN *prn)
{
    switcher s;
    double lldiff = NADBL;
    double llbak = -1.0e+200;
    double eps1 = 0.0001;
    double eps2 = 0.005;
    double stol = 0.01 * eps1;
    double wtol = 0.01 * eps2;
    int j, jmax = 50000;
    int wcount = 0;
    int uinit = 0;
    int conv = 0;
    int err;

    /* old: tol = 2.0e-11 */

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

#if JDEBUG > 1
	print_switch_iter(J, j);
#endif

	if (!err && j > 1) {
	    lldiff = J->ll - llbak;
	    if (lldiff < stol) {
		conv = 1;
		break;
	    } else if (lldiff < wtol) {
		if (++wcount == 4) {
		    conv = 2;
		    break;
		}
	    } else {
		wcount = 0;
	    }
	}

	llbak = J->ll;
    }

    if (err) {
	pputs(prn, _("Switching algorithm: failed"));
	pputc(prn, '\n');
    } else {
	pprintf(prn, _("Switching algorithm: %d iterations"), j);

	if (uinit) {
	    pprintf(prn, " (%s)", _("user-supplied initial values"));
	}
	pputc(prn, '\n');

	if (!na(J->ll) && !na(lldiff)) {
	    pprintf(prn, " -(T/2)log|Omega| = %.8g, lldiff = %g\n", J->ll, lldiff);
	}

	if (conv == 1) {
	    pprintf(prn, "%s\n", _("Strong convergence")); 
	} else if (conv == 2) {
	    pprintf(prn, "%s\n", _("Weak convergence")); 
	}
    }

    if (!err && !conv) {
	err = E_NOCONV;
    } else if (j > 1500) {
	pputs(prn, "*** warning: large number of iterations may "
	      "indicate a problem");
	pputc(prn, '\n');
    }

    pputc(prn, '\n');

    if (!err) {
	J->ll -= J->llk;
    }

    switcher_free(&s);

    return err;
}

/* re-scale beta and alpha based on the saved scaling
   information
*/

static int replace_col_scaling (Jwrap *J)
{
    double x, s;
    int i, j, k, ii;
    int err = 0;

    for (i=0; i<J->normcol[0]; i++) {
	j = J->normcol[i+1];
	if (j < 0) {
	    continue;
	}
	k = j / J->p1;
	s = gretl_matrix_get(J->beta, j % J->p1, k) / J->normval[i];
	for (ii=0; ii<J->p1; ii++) {
	    x = gretl_matrix_get(J->beta, ii, k);
	    if (x != 0) {
		gretl_matrix_set(J->beta, ii, k, x / s);
	    }
	}
	for (ii=0; ii<J->p; ii++) {
	    x = gretl_matrix_get(J->alpha, ii, k);
	    if (x != 0) {
		gretl_matrix_set(J->alpha, ii, k, x * s);
	    }
	}
	J->blen -= 1;
    }

    return err;
}

static int reset_null_H (Jwrap *J)
{
    fprintf(stderr, "** Doing reset_null_H\n");

    J->H0 = J->H;
    gretl_matrix_free(J->h);
    J->H = J->h = NULL;
    
    J->blen = J->r * J->p1;

    gretl_matrix_free(J->phi);
    J->phi = gretl_zero_matrix_new(J->blen, 1);
    if (J->phi == NULL) {
	return E_ALLOC;
    }

    return 0;
}

/* recalculate the H matrix after removing column
   scalings */

static int reset_H (Jwrap *J, const gretl_matrix *R,
		    const gretl_matrix *q)
{
    int err;

    J->H0 = J->H;
    gretl_matrix_free(J->h);
    J->h = NULL;

    err = real_set_up_H(J, R, q);

    if (!err) {
	gretl_matrix_free(J->phi);
	J->phi = gretl_zero_matrix_new(J->blen, 1);
	if (J->phi == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static const gretl_matrix *
J_get_R_matrix (Jwrap *J, const gretl_restriction *rset)
{
    if (J->bigR != NULL) {
	return J->bigR;
    } else {
	return rset_get_R_matrix(rset);
    }
}

static const gretl_matrix *
J_get_q_matrix (Jwrap *J, const gretl_restriction *rset)
{
    if (J->bigq != NULL) {
	return J->bigq;
    } else {
	return rset_get_q_matrix(rset);
    }
}

static int get_master_col (Jwrap *J, int pos, int *mpos)
{
    int i, k = -J->normcol[pos];
    int ci, bcol = k / J->p1;

    for (i=0; i<J->normcol[0]; i++) {
	ci = J->normcol[i+1];
	if (ci >= 0 && ci / J->p1 == bcol) {
	    *mpos = i + 1;
	    return ci;
	}
    }

    return -1;
}

/* Take a row from the original R matrix and re-write it as a
   homogenous restriction referencing the "master" value which 
   was used as a basis for scaling the relevant column of beta. 
   For example, if we found the scaling b[2,1] = 1, and then the 
   following restriction b[2,2] = -1, we remove the first 
   restriction (temporarily), and we replace the second with 
   b[1,1] + b[2,1] = 0.
*/

static int homogenize_R_line (gretl_matrix *R,
			      gretl_matrix *q,
			      Jwrap *J, int pos, 
			      int k)
{
    int j, j0, j1, mpos;
    double x;

    j0 = get_master_col(J, pos, &mpos);
    if (j0 < 0) {
	return E_DATA;
    }

    j1 = -J->normcol[pos];
    x = -J->normval[mpos-1] / J->normval[pos-1];

    for (j=0; j<R->cols; j++) {
	if (j == j0) {
	    gretl_matrix_set(R, k, j, 1);
	} else if (j == j1) {
	    gretl_matrix_set(R, k, j, x);
	} else {
	    gretl_matrix_set(R, k, j, 0);
	}
    }

    gretl_vector_set(q, k, 0);

    return 0;
}

static void norm_destroy (Jwrap *J)
{
    free(J->normrow);
    free(J->normcol);
    free(J->normval);

    J->normrow = NULL;
    J->normcol = NULL;
    J->normval = NULL;
}

static int norm_allocate (Jwrap *J, const gretl_matrix *R)
{
    int n = R->rows;
    int err = 0;

    J->normrow = gretl_list_new(n);
    if (J->normrow == NULL) err = E_ALLOC;

    if (!err) {
	J->normcol = gretl_list_new(n);
	if (J->normcol == NULL) err = E_ALLOC;
    }

    if (!err) {
	J->normval = malloc(n * sizeof *J->normval);
	if (J->normval == NULL) err = E_ALLOC;
    }

    if (err) {
	norm_destroy(J);
    } else {
	J->normrow[0] = 0;
	J->normcol[0] = 0;
    }  

    return err;
}

static void list_push (int *list, int k)
{
    list[0] += 1;
    list[list[0]] = k;
}

/* test whether or not a row of the R matrix (taking
   into account the corresponding q value) represents
   a scaling of a column of beta
*/

static int is_scaling_row (const gretl_matrix *R, 
			   const gretl_matrix *q,
			   int i, int *col,
			   double *sval)
{
    int j, c = 0;
    double rij;

    for (j=0; j<R->cols; j++) {
	rij = gretl_matrix_get(R, i, j);
	if (rij != 0 && q->val[i] != 0) {
	    c++;
	    if (c > 1) {
		break;
	    } else {
		*sval = q->val[i] / rij;
		*col = j;
	    }
	} 
    }

    return c == 1;
}

/* Check for restrictions which make scale removal impossible,
   for each column of beta.  Record a -1 in the appropriate
   place in the @sc array if scaling won't work for that column.

   We consider scale removal infeasible for a given beta column
   if (a) we find a non-homogeneous restriction involving two
   or more beta elements (e.g. b[1,1] - b[1,3] = 2), or (b) we
   find any cross-column restrictions.

   Return 1 if there are no scaleable columns on this test, else 0.
*/

static int screen_beta_cols (Jwrap *J,
			     const gretl_matrix *R,
			     const gretl_matrix *q,
			     char *sc)
{
    double x;
    int ns = 0;
    int i, j, k, m, c;

    for (i=0; i<R->rows; i++) {
        k = -1;
	c = 0;
	for (j=0; j<R->cols; j++) {
	    x = gretl_matrix_get(R, i, j);
	    if (x != 0) {
		c++;
		if (k < 0) {
		    k = j / J->p1;
		} else {
		    m = j / J->p1;
		    if (m != k) {
			/* cross-column restriction */
			sc[k] = sc[m] = -1;
		    } 
		}
	    }
	}
	if (c > 1 && q->val[i] != 0) { 
	    sc[k] = -1;
	}
    }

    for (i=0; i<J->r; i++) {
	if (sc[i] == 0) ns++;
    }

    return (ns == 0);
}

/* Check for restrictions (rows of R) which are just scalings of a
   given column of beta.  There can only be one such scaling per
   beta column.
*/

static int check_for_scaling (Jwrap *J,
			      const gretl_matrix *R,
			      const gretl_matrix *q)
{
    char *sc = calloc(J->r, 1);
    double sval = 0;
    int i, k, rcol = 0;

    if (sc == NULL) {
	return E_ALLOC;
    }

    /* first pass: screen for conditions which make
       scale-removal impossible */
    if (screen_beta_cols(J, R, q, sc)) {
	free(sc);
	return 0;
    }

    if (norm_allocate(J, R)) {
	return E_ALLOC;
    }

#if 0 /* Sven's Problem #1 */
    sc[2] = -1;
#endif

    for (i=0; i<R->rows; i++) {
	if (is_scaling_row(R, q, i, &rcol, &sval)) {
	    k = rcol / J->p1;
	    if (sc[k] >= 0) {
		/* column is scaleable */
		J->normval[J->normrow[0]] = sval;
		list_push(J->normrow, i);
		if (sc[k] == 0) {
		    /* beta column k is not yet scaled */
		    list_push(J->normcol, rcol);
		} else {
		    /* negative rcol entry indicates a "follower" */
		    list_push(J->normcol, -rcol);
		}
		sc[k] += 1;
	    }
	}
    }

#if JDEBUG
    printlist(J->normrow, "norm rows");
    printlist(J->normcol, "norm cols");
    for (i=0; i<J->normrow[0]; i++) {
	fprintf(stderr, "normval[%d] = %g\n", i, J->normval[i]);
    }
#endif

    if (J->normrow[0] == 0) {
	/* no (practicable) scalings found */
	norm_destroy(J);
    }

    free(sc);

    return 0;
}

static int just_identified (Jwrap *J)
{
    int npar = J->alen + J->blen;

    return (J->jr >= npar && J->df == 0);
}

static int 
maybe_remove_col_scaling (Jwrap *J,
			  const gretl_restriction *rset)
{
    const gretl_matrix *R = J_get_R_matrix(J, rset);
    const gretl_matrix *q = J_get_q_matrix(J, rset);
    gretl_matrix *Rr = NULL;
    gretl_matrix *qr = NULL;
    double x;
    int i, j, k, rr, pos;
    int err = 0;

#if JDEBUG
    fprintf(stderr, "maybe_remove_col_scaling...\n");
#endif

    if (R == NULL || q == NULL || gretl_is_zero_matrix(q)) {
	return 0;
    }

    if (using_bfgs(J) || cross_alpha(J) || doing_manual_init(J)) {
	return 0;
    }

#if 1 /* ?? */
    if (just_identified(J)) {
	return 0;
    }
#endif

    err = check_for_scaling(J, R, q);
    if (err || !do_scaling(J)) {
	return err;
    } 
    
    rr = R->rows;
    for (i=0; i<J->normcol[0]; i++) {
	if (J->normcol[i+1] >= 0) {
	    rr--;
	}
    }

    if (rr == 0) {
	/* all the restrictions were column scalings */
	return reset_null_H(J);
    }

    Rr = gretl_matrix_alloc(rr, R->cols);
    qr = gretl_column_vector_alloc(rr);

    if (Rr == NULL || qr == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    k = 0;
    for (i=0; i<R->rows; i++) {
	pos = in_gretl_list(J->normrow, i);
	if (pos > 0 && J->normcol[pos] < 0) {
	    /* convert to homogeneous restriction */
	    homogenize_R_line(Rr, qr, J, pos, k);
	    k++;
	} else if (pos == 0) {
	    /* pass the row through unmodified */
	    for (j=0; j<R->cols; j++) {
		x = gretl_matrix_get(R, i, j);
		gretl_matrix_set(Rr, k, j, x);
		x = gretl_vector_get(q, i);
		gretl_vector_set(qr, k, x);
	    }
	    k++;
	} 
    }

#if JDEBUG
    gretl_matrix_print(Rr, "Rr");
    gretl_matrix_print(qr, "qr");
#endif

    err = reset_H(J, Rr, qr);

 bailout:

    gretl_matrix_free(Rr);
    gretl_matrix_free(qr);

    return err;
}

/* 
   J = [(I_p \otimes \beta)G : (\alpha \otimes I_{p1})H]

   Boswijk and Doornik (2004), equation (40), page 455.
*/

static int check_jacobian (Jwrap *J)
{
    gretl_matrix *L = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *Jac = NULL;
    int err = 0;

    /* form beta and alpha, based on randomized theta */

    if (J->H != NULL) {
	gretl_matrix *phi;

	if (J->H->rows != J->p1 * J->r) {
	    fprintf(stderr, "?? matrices wrongly sized in check_jacobian\n");
	    return E_NONCONF;
	} 

	phi = gretl_column_vector_alloc(J->H->cols);
	if (phi == NULL) {
	    return E_ALLOC;
	}

	gretl_matrix_random_fill(phi, D_NORMAL);
	gretl_matrix_reuse(J->beta, J->p1 * J->r, 1);
	err = gretl_matrix_multiply(J->H, phi, J->beta);
	if (!err && J->h != NULL) {
	    err = gretl_matrix_add_to(J->beta, J->h);
	}
	gretl_matrix_reuse(J->beta, J->p1, J->r);
	gretl_matrix_free(phi);
    } else if (J->blen > 0) {
	gretl_matrix_random_fill(J->beta, D_NORMAL);
    }

    if (!err && J->G != NULL) {
	gretl_matrix *psi = gretl_column_vector_alloc(J->G->cols);
	gretl_matrix *avec = gretl_column_vector_alloc(J->p * J->r);

	if (psi == NULL || avec == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}	

	gretl_matrix_random_fill(psi, D_NORMAL);
	gretl_matrix_multiply(J->G, psi, avec);
	unvec_transpose(J->alpha, avec);
	gretl_matrix_free(psi);
	gretl_matrix_free(avec);
    } else if (!err) {	
	J_compute_alpha(J);
    }

    if (!err) {
	L = gretl_matrix_I_kronecker_new(J->p, J->beta, &err);
    }

    if (!err) {
	R = gretl_matrix_kronecker_I_new(J->alpha, J->p1, &err);
    }

    if (!err && J->G != NULL) {
	/* alpha is restricted */
	gretl_matrix *LG = gretl_matrix_multiply_new(L, J->G, &err);

	if (!err) {
	    gretl_matrix_replace(&L, LG);
	}
    }

    if (!err && J->H != NULL) {
	/* beta is restricted */
	gretl_matrix *RH = gretl_matrix_multiply_new(R, J->H, &err);

	if (!err) {
	    gretl_matrix_replace(&R, RH);
	}
    }

    if (!err) {
	Jac = gretl_matrix_col_concat(L, R, &err);
    }

    if (!err) {
	J->jr = gretl_matrix_rank(Jac, &err);
#if JDEBUG
	gretl_matrix_print(Jac, "Jacobian");
#endif
    }

 bailout:

    gretl_matrix_free(L);
    gretl_matrix_free(R);
    gretl_matrix_free(Jac);

#if JDEBUG
    fprintf(stderr, "check_jacobian: returning %d\n", err);
#endif

    return err;
}

/* Doornik's approach to checking identification */

static int vecm_id_check (Jwrap *J, GRETL_VAR *jvar, 
			  gretlopt opt, PRN *prn)
{
    int npar = J->alen + J->blen;
    int err = 0;

    err = check_jacobian(J);

#if JDEBUG
    fprintf(stderr, "vecm_id_check: check_jacobian ret'd %d\n", err);
#endif

    if (!err) {
	pprintf(prn, _("Rank of Jacobian = %d, number of free "
		"parameters = %d\n"), J->jr, npar);
	if (J->jr < npar) {
	    pputs(prn, _("Model is not fully identified\n"));
	    if (using_bfgs(J)) {
		err = E_NOIDENT;
	    }
	} else {
	    pputs(prn, _("Model is fully identified\n"));
	}

	J->df = (J->p + J->p1 - J->r) * J->r - J->jr;
	pprintf(prn, _("Based on Jacobian, df = %d\n"), J->df);

#if JDEBUG
	fprintf(stderr, "Jacobian: got df = %d\n", J->df);
#endif

	if (!(opt & OPT_B)) {
	    /* system was subject to a prior restriction? */
	    if (jvar->jinfo->lrdf > 0) {
		J->df -= jvar->jinfo->lrdf;
		pprintf(prn, _("Allowing for prior restriction, df = %d\n"), 
			J->df);
	    }
	}
    }

    return err;
}

/* common alpha restriction needs to be replicated */

static int G_from_expanded_R (Jwrap *J, const gretl_matrix *R)
{
    gretl_matrix *Rtmp = NULL;
    gretl_matrix *Rcpy = NULL;
    double x;
    int c0, c1 = 0;
    int i, j, k;
    int err = 0;

    Rtmp = gretl_matrix_I_kronecker_new(J->r, R, &err);

    if (!err) {
	Rcpy = gretl_matrix_copy(Rtmp);
	if (Rcpy == NULL) {
	    err = E_ALLOC;
	}
    }

#if JDEBUG
    gretl_matrix_print(Rtmp, "expanded R");
#endif

    for (i=0; i<J->p; i++) {
	for (j=0; j<J->r; j++) {
	    /* write col of orig R into col of remapped R */
	    c0 = i + j * J->p;
	    for (k=0; k<Rcpy->rows; k++) {
		x = gretl_matrix_get(Rcpy, k, c0);
		gretl_matrix_set(Rtmp, k, c1, x);
	    }
	    c1++;
	}
    }

#if JDEBUG
    gretl_matrix_print(Rtmp, "remapped expanded R");
#endif

    if (!err) {
	J->G = gretl_matrix_right_nullspace(Rtmp, &err);
    }

    gretl_matrix_free(Rtmp);
    gretl_matrix_free(Rcpy);

    return err;
}

static void nullspace_normalize (gretl_matrix *A)
{
    int i, j, nz, nm1, pos;
    double aij;

    for (i=0; i<A->rows; i++) {
	nz = nm1 = pos = 0;
	for (j=0; j<A->cols; j++) {
	    aij = gretl_matrix_get(A, i, j);
	    if (aij == 0.0) {
		nz++;
	    } else if (aij == -1.0) {
		nm1++;
		pos = j;
	    } else {
		break;
	    }
	}
	if (nz == A->cols - 1 && nm1 == 1) {
	    /* got a solitary -1 on row i */
	    gretl_matrix_set(A, i, pos, 1.0);
	}
    }
}

/* user-supplied R needs to be remapped for conformity with 
   vec(\alpha') = G\psi */

static int G_from_remapped_R (Jwrap *J, const gretl_matrix *R)
{
    gretl_matrix *Rtmp = NULL;
    double x;
    int c0, c1 = 0;
    int i, j, k;
    int err = 0;

    Rtmp = gretl_matrix_copy(R);
    if (Rtmp == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<J->p; i++) {
	for (j=0; j<J->r; j++) {
	    /* write col of orig R into col of remapped R */
	    c0 = i + j * J->p;
	    for (k=0; k<R->rows; k++) {
		x = gretl_matrix_get(R, k, c0);
		gretl_matrix_set(Rtmp, k, c1, x);
	    }
	    c1++;
	}
    }

#if JDEBUG
    gretl_matrix_print(Rtmp, "remapped R");
#endif

    J->G = gretl_matrix_right_nullspace(Rtmp, &err);
    gretl_matrix_free(Rtmp);

    return err;
}

/* detect the presence of a cross-equation alpha
   restriction */

static int 
cross_alpha_check (Jwrap *J, const gretl_matrix *R)
{
    double x;
    int i, j, k = -1;

    for (i=0; i<R->rows; i++) {
	for (j=0; j<R->cols; j++) {
	    x = gretl_matrix_get(R, i, j);
	    if (x != 0) {
		if (k < 0) {
		    k = j / J->p;
		} else if (j / J->p != k) {
		    fprintf(stderr, "Got a cross-alpha restriction\n");
		    J->flags |= J_CROSS_ALPHA;
		    return 1;
		}
	    }
	}
    }

    return 0;
}

/* set up restriction matrix, G, for alpha */

static int set_up_G (Jwrap *J, const gretl_restriction *rset)
{
    const gretl_matrix *Ra = rset_get_Ra_matrix(rset);
    const gretl_matrix *qa = rset_get_qa_matrix(rset);
    int err = 0;

#if JDEBUG
    gretl_matrix_print(Ra, "Ra, in set_up_G");
    gretl_matrix_print(qa, "qa, in set_up_G");
#endif

    if (qa != NULL && !gretl_is_zero_matrix(qa)) {
	gretl_errmsg_set("alpha restriction is not homogeneous: "
			 "not supported");
	return E_NOTIMP;
    }

    if (J->r > 1) {
	if (Ra->cols == J->p) {
	    /* got a common alpha restriction */
	    err = G_from_expanded_R(J, Ra);
	} else if (Ra->cols > J->p) {
	    err = G_from_remapped_R(J, Ra);
	}
    } else {
	J->G = gretl_matrix_right_nullspace(Ra, &err);
    }

    if (!err) {
	J->alen = J->G->cols;
    }

    if (!err) {
	nullspace_normalize(J->G);
    }

    if (!err) {
	cross_alpha_check(J, Ra);
    }

#if JDEBUG
    gretl_matrix_print(J->G, "G, in set_up_G");
#endif

    return err;
}

static int real_set_up_H (Jwrap *J, const gretl_matrix *R,
			  const gretl_matrix *q)
{
    gretl_matrix *RRT = NULL;
    gretl_matrix *Tmp = NULL;
    int err = 0;

#if JDEBUG
    gretl_matrix_print(R, "R, in real_set_up_H");
    gretl_matrix_print(q, "q, in real_set_up_H");
#endif

    J->H = gretl_matrix_right_nullspace(R, &err);

    if (err) {
	return err;
    }

    J->blen = J->H->cols;

    if (q == NULL || gretl_is_zero_matrix(q)) {
	/* implicitly or explicitly zero */
	J->h = gretl_zero_matrix_new(R->cols, 1);
#if JDEBUG
	gretl_matrix_print(J->H, "H, in real_set_up_H");
	gretl_matrix_print(J->h, "h, in real_set_up_H");
#endif
	return (J->h == NULL)? E_ALLOC : 0; 
    }

    /* now for h, for non-zero q */

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
	J->h = gretl_matrix_multiply_new(Tmp, q, &err);
    }

#if JDEBUG
    gretl_matrix_print(J->H, "H, in real_set_up_H");
    gretl_matrix_print(J->h, "h, in real_set_up_H");
#endif

    gretl_matrix_free(RRT);
    gretl_matrix_free(Tmp);

    return err;
}

static gretl_matrix *replicate_q (int r, const gretl_matrix *q0,
				  int *err)
{
    gretl_matrix *q;
    int k = q0->rows * r;
    int i, j;

    q = gretl_column_vector_alloc(k);

    if (q == NULL) {
	*err = E_ALLOC;
    } else {
	k = 0;
	for (j=0; j<r; j++) {
	    for (i=0; i<q0->rows; i++) {
		q->val[k++] = q0->val[i];
	    }
	}
    }

    return q;
}

/* set up restriction matrices, H and h, for beta */

static int set_up_H (Jwrap *J, const gretl_restriction *rset)
{
    const gretl_matrix *R = rset_get_R_matrix(rset);
    const gretl_matrix *q = rset_get_q_matrix(rset);
    int qzero, err = 0;

#if JDEBUG
    gretl_matrix_print(R, "R, in set_up_H");
    gretl_matrix_print(q, "q, in set_up_H");
#endif

    if (R->rows == J->p1 * J->r) {
	/* number of restrictions = total betas */
	return solve_for_beta(J, R, q);
    }

    qzero = (q == NULL || gretl_is_zero_matrix(q));

    if (J->r > 1 && R->cols == J->p1) {
	/* common beta restriction */
	J->bigR = gretl_matrix_I_kronecker_new(J->r, R, &err);
	if (!err) {
	    R = J->bigR;
	    if (!qzero) {
		J->bigq = replicate_q(J->r, q, &err);
		if (!err) {
		    q = J->bigq;
		}
	    }
	}
    } 

    if (!err) {
	err = real_set_up_H(J, R, q);
    }

    if (!err) {
	J->H0 = J->H;
    }

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

    if (J->h == NULL || 
	gretl_is_zero_matrix(J->h) ||
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

    E = gretl_symm_matrix_eigenvals_descending(BB, 1, &err);

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
				  J->h, GRETL_MOD_NONE,
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
	if (gretl_is_zero_matrix(J->phi)) {
	    fprintf(stderr, "Got a zero vector for phi\n");
	    gretl_matrix_print(IBPh, "IBPh");
	} else {
	    gretl_matrix_switch_sign(J->phi);
	}
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
    } else if (gretl_is_zero_matrix(J->h)) {
	err = phi_init_homog(J);
    } else {
	/* solve for \phi a la Boswijk */
	err = phi_init_nonhomog(J);
    }

    return err;
}

/* solve the unrestricted eigenvalue problem for \beta */

static int case0 (Jwrap *J)
{
    gretl_matrix *B = NULL;
    int err = 0;

    B = gretl_matrix_alloc(J->p1, J->p);
    if (B == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_matrix_SVD_johansen_solve(J->R0, J->R1,
					      NULL, B, NULL,
					      J->r);
    }

#if JDEBUG
    gretl_matrix_print(B, "case0: raw beta");
#endif

    if (!err) {
	gretl_matrix_copy_values(J->beta, B);
    }

    gretl_matrix_free(B);

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

/* do initial computation of alpha based on beta (unless this is
   already done); then, if relevant, impose alpha restriction and make
   psi.
*/

static int alpha_init (Jwrap *J)
{
    int err = 0;

    err = J_compute_alpha(J);
    if (err) {
	return err;
    }

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
    int i, j;

    J->bse = gretl_matrix_alloc(J->p1, J->r);
    if (J->bse == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<J->Vb->rows; i++) {
	for (j=0; j<J->Vb->cols; j++) {
	    x = gretl_matrix_get(J->Vb, i, j);
	    if (fabs(x) < 1.0e-20) {
		gretl_matrix_set(J->Vb, i, j, 0);
		x = 0;
	    }
	    if (j == i) {
		J->bse->val[i] = sqrt(x);
	    }
	}
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

static int max_beta_vname (char *tmp,
			   GRETL_VAR *v, 
			   const DATASET *dset)
{
    int r = gretl_matrix_rows(v->jinfo->Beta);
    const char *s;
    int i, ni, n = 0;

    for (i=0; i<r; i++) {    
	s = vecm_beta_varname(tmp, v, dset, i);
	ni = strlen(s);
	if (ni > n) {
	    n = ni;
	}
    }

    return n;
}

#define VECM_WIDTH 13

static int printres (Jwrap *J, GRETL_VAR *jvar, 
		     gretl_restriction *rset,
		     const DATASET *dset,
		     PRN *prn)
{
    const gretl_matrix *c = J->beta;
    const gretl_matrix *sd = J->bse;
    char s[16], vname[NAMETRUNC];
    char namefmt[8];
    int nwid, sdshow;
    int i, j;

    if (J->df > 0) {
	pprintf(prn, _("Unrestricted loglikelihood (lu) = %.8g\n"), jvar->ll);
	pprintf(prn, _("Restricted loglikelihood (lr) = %.8g\n"), J->ll);
    } else {
	char ll0[32], ll1[32];

	pprintf(prn, "%s = %.8g\n", _("loglikelihood"), J->ll);
	sprintf(ll0, "%.8g", jvar->ll);
	sprintf(ll1, "%.8g", J->ll);
	if (strcmp(ll0, ll1)) {
	    pprintf(prn, "*** warning: should equal %.8g\n", jvar->ll);
	}
    }

    if (J->df > 0) {
	double x = 2.0 * (jvar->ll - J->ll);
	double pv = chisq_cdf_comp(J->df, x);

	pprintf(prn, "2 * (lu - lr) = %g\n", x);
	pprintf(prn, "P(%s(%d) > %g) = %g\n", _("Chi-square"), J->df, x, pv);
	rset_add_results(rset, x, pv, J->ll);
    }

    sdshow = (sd != NULL && !gretl_is_zero_matrix(sd));

    pputc(prn, '\n');
    pputs(prn, _("Cointegrating vectors"));
    if (sdshow) {
	pprintf(prn, " (%s)", _("standard errors in parentheses"));
    }
    pputs(prn, "\n\n");

    nwid = max_beta_vname(vname, jvar, dset) + 1;
    sprintf(namefmt, "%%-%ds", nwid);

    for (i=0; i<J->p1; i++) {
	vecm_beta_varname(vname, jvar, dset, i);
	pprintf(prn, namefmt, vname);

	for (j=0; j<J->r; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(c, i, j));
	}
	pputc(prn, '\n');

	if (sdshow) {
	    bufspace(nwid + 1, prn);
	    for (j=0; j<J->r; j++) {
		sprintf(s, "(%#.5g)", gretl_matrix_get(sd, i, j));
		pprintf(prn, "%12s ", s);
	    }
	    pputc(prn, '\n');
	}
    }

    c = J->alpha;
    sd = J->ase;

    sdshow = (sd != NULL && !gretl_is_zero_matrix(sd));

    pputc(prn, '\n');
    pputs(prn, _("alpha (adjustment vectors)"));
    if (sdshow) {
	pprintf(prn, " (%s)", _("standard errors in parentheses"));
    }
    pputs(prn, "\n\n");

    for (i=0; i<J->p; i++) {
	maybe_trim_varname(vname, dset->varname[jvar->ylist[i+1]]);
	pprintf(prn, namefmt, vname);
	for (j=0; j<J->r; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(c, i, j));
	}
	pputc(prn, '\n');

	if (sdshow) {
	    bufspace(nwid + 1, prn);
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
	err = gretl_invert_symmetric_matrix(J->qf1);
	if (err) {
	    gretl_matrix_print(J->qf1, "J->qf1: couldn't invert");
	}
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
	fprintf(stderr, "jrestrict: make_theta: k = 0\n");
	return E_DATA;
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

static void transcribe_to_jvar (Jwrap *J, GRETL_VAR *jvar,
				gretlopt opt)
{
    if (!(opt & OPT_B)) {
	jvar->jinfo->ll0 = jvar->ll;
	jvar->ll = J->ll;
	jvar->jinfo->lrdf += J->df; /* ? */
    }

    gretl_matrix_replace(&jvar->S, J->Omega);
    J->Omega = NULL;

    gretl_matrix_replace(&jvar->jinfo->Beta, J->beta);
    J->beta = NULL;

    gretl_matrix_replace(&jvar->jinfo->Alpha, J->alpha);
    J->alpha = NULL;

    if (opt & OPT_B) {
	/* bootstrapping: we don't need the variances */
	return;
    }

    gretl_matrix_replace(&jvar->jinfo->Bvar, J->Vb);
    J->Vb = NULL;

    gretl_matrix_replace(&jvar->jinfo->Bse, J->bse);
    J->bse = NULL;

    gretl_matrix_replace(&jvar->jinfo->Ase, J->ase);
    J->ase = NULL; 
}

static int do_bfgs (Jwrap *J, gretlopt opt, PRN *prn)
{
    int err = add_gradhelper(J);

    if (!err) {
	int maxit = 4000, fncount = 0, grcount = 0;
	int nn = J->theta->rows;
	double toler = 1.0e-11;

	pputs(prn, "LBFGS: using analytical derivatives\n\n");

	err = LBFGS_max(J->theta->val, nn, maxit, toler, 
			&fncount, &grcount, Jloglik, C_LOGLIK,
			Jgradient, J, opt, prn);
    }

    return err;
}

/* 
   OPT_L: use LBFGS approach instead of switching algorithm.
   OPT_F: doing full estimation of restricted system, not just 
          testing the restriction.
   OPT_J: ("jitter") use simulated annealing in initialization.
   OPT_B: bootstrapping impulse response
*/

int general_vecm_analysis (GRETL_VAR *jvar, 
			   gretl_restriction *rset,
			   const DATASET *dset,
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
	err = set_up_H(J, rset);
    }

    if (!err) {
	err = allocate_phi(J);
    }

    if (!err && rset_VECM_acols(rset) > 0) {
	err = set_up_G(J, rset);
    }

    if (!err) {
	err = allocate_psi(J);
    }   

    if (!err) {
	err = vecm_id_check(J, jvar, opt, prn);
	if (J->df == 0) {
	    fprintf(stderr, "warning: test df = 0\n");
	    J->flags &= ~J_USE_LBFGS;
	    do_simann = 0;
	}
    }

    if (!err && J->blen == 0 && J->G == NULL) {
	/* beta doesn't need to be estimated */
	J->df = (J->p1 - J->r) * J->r;
	err = real_compute_ll(J);
	goto skipest;
    }

    if (!err && !(opt & OPT_N)) {
	err = maybe_remove_col_scaling(J, rset);
    }

    if (!err) {
	err = beta_init(J);
    }

    if (!err && (!using_bfgs(J) || J->G != NULL)) {
	err = alpha_init(J);
    }

    if (!err && (using_bfgs(J) || do_simann)) {
	err = make_theta(J);
    }

    if (!err && do_simann) {
	int n = gretl_vector_get_length(J->theta);

	err = gretl_simann(J->theta->val, n, 4096,
			   Jloglik, J, opt, prn);
	if (!err) {
	    sync_with_theta(J, J->theta->val);
	}
    }

    if (!err) {
	if (using_bfgs(J)) {
	    err = do_bfgs(J, opt, prn);
	} else {
	    err = switchit(J, prn);
	}
    }

    if (!err && do_scaling(J)) {
	err = replace_col_scaling(J);
    }     

 skipest:

    if (!err && !(opt & OPT_B)) {
	err = variance_from_info_matrix(J, jvar);
    }

    if (!err) {
	if (full) {
	    transcribe_to_jvar(J, jvar, OPT_F);
	} else if (opt & OPT_B) {
	    transcribe_to_jvar(J, jvar, OPT_B);
	} else {
	    printres(J, jvar, rset, dset, prn);
	}
    } 

    jwrap_destroy(J);

    return err;
}
