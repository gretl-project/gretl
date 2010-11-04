/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2010 Allin Cottrell and Riccardo "Jack" Lucchetti
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
#include "matrix_extra.h"
#include "gretl_bfgs.h"
#include "missing_private.h"

#define BIPDEBUG 0

typedef struct bp_container_ bp_container;

struct bp_container_ {
    const int *list;         /* incoming model specification */
    int t1, t2;              /* start and end of sample range */
    double ll0;		     /* log-likelihood for rho = 0*/
    double ll;		     /* log-likelihood */
    gretl_matrix *score;     /* score matrix */
    gretl_vector *sscore;    /* score vector (sum) */

    int ntot; 	             /* total obs */
    int k1, k2;              /* number of regressors for each eq. */
    int npar;                /* number of parameters */
    char *mask;              /* mask for missing values */
    int depvar1;	     /* location of y1 in array Z */
    int depvar2;	     /* location of y2 in array Z */
    int *X1list;	     /* regressor list for the first eq. */
    int *X2list;	     /* regressor list for the second eq. */

    int *s1;	     	     /* first dependent var */
    int *s2;		     /* second dependent var */
    gretl_matrix *reg1;	     /* first eq. regressors */
    gretl_matrix *reg2;	     /* second eq. regressors */
    gretl_vector *fitted1;   /* x_1'\beta_1 */
    gretl_vector *fitted2;   /* x_2'\beta_2 */
    gretl_vector *u1;
    gretl_vector *u2;

    gretl_vector *beta;	     /* first eq. parameters */
    gretl_vector *gama;      /* second eq. parameters */
    double arho;             /* atan(rho) */

    gretl_matrix *H11;       /* workspace for the analytical Hessian */
    gretl_matrix *H12;
    gretl_matrix *H13;
    gretl_matrix *H22;
    gretl_matrix *H23;

    gretl_matrix *vcv;	     /* Variance-covariance matrix */
};

/* ---------------------------------------------------------------------- */
/* container-bookkeeping functions                                        */
/* ---------------------------------------------------------------------- */

static void bp_container_destroy (bp_container *BPC)
{
    if (BPC == NULL) {
	return;
    }

    free(BPC->mask);

    free(BPC->X1list);
    free(BPC->X2list);

    free(BPC->s1);
    free(BPC->s2);
    gretl_matrix_free(BPC->reg1);
    gretl_matrix_free(BPC->reg2);
    gretl_vector_free(BPC->fitted1);
    gretl_vector_free(BPC->fitted2);
    gretl_vector_free(BPC->u1);
    gretl_vector_free(BPC->u2);
    gretl_matrix_free(BPC->score);
    gretl_vector_free(BPC->sscore);

    gretl_vector_free(BPC->beta);
    gretl_vector_free(BPC->gama);

    gretl_matrix_free(BPC->vcv);

    gretl_matrix_free(BPC->H11);
    gretl_matrix_free(BPC->H12);
    gretl_matrix_free(BPC->H13);
    gretl_matrix_free(BPC->H22);
    gretl_matrix_free(BPC->H23);

    free(BPC);
}

static bp_container *bp_container_new (const int *list)
{
    bp_container *BPC = malloc(sizeof *BPC);

    if (BPC == NULL) {
	return NULL;
    }

    BPC->mask = NULL;

    BPC->list = list;
    BPC->t1 = BPC->t2 = 0;
    BPC->ll = NADBL;

    BPC->X1list = NULL;
    BPC->X2list = NULL;

    BPC->reg1 = NULL;
    BPC->reg2 = NULL;
    BPC->fitted1 = NULL;
    BPC->fitted2 = NULL;
    BPC->u1 = NULL;
    BPC->u2 = NULL;
    BPC->score = NULL;
    BPC->sscore = NULL;

    BPC->beta = NULL;
    BPC->gama = NULL;
    BPC->arho = 0;

    BPC->vcv = NULL;

    return BPC;
}

/* ---------------------------------------------------------------------- */
/* Initialization-related functions                                       */
/* ---------------------------------------------------------------------- */

static int bp_cont_fill_data (bp_container *BPC, const double **Z)
{
    int i, T, err = 0;
    char *mask = BPC->mask;
    int t1 = BPC->t1;
    int t2 = BPC->t2;

    T = t2 - t1 + 1;

    int nmissing = 0;
    if (mask != NULL) {
	for (i=0; i<T; i++) {
	    if (mask[i]) nmissing++;
	}
    } 

    if (T == nmissing) {
	return E_DATA;
    }

    BPC->ntot = T - nmissing;

    int *s1;
    int *s2;
    s1 = malloc(BPC->ntot * sizeof *s1);
    s2 = malloc(BPC->ntot * sizeof *s2);

    if (s1==NULL || s2==NULL) {
	return E_ALLOC;
    }

    if (nmissing == 0) {
	for (i=t1; i<=t2; i++) {
	    s1[i] = (Z[BPC->depvar1][i] == 1.0);
	    s2[i] = (Z[BPC->depvar2][i] == 1.0);
	}
    } else {
	int k = 0;
	for (i=t1; i<=t2; i++) {
	    if (!mask[i]) {
		s1[k] = (Z[BPC->depvar1][i] == 1.0);
		s2[k] = (Z[BPC->depvar2][i] == 1.0);
		k++;
	    }
	}
    } 

    BPC->s1 = s1;
    BPC->s2 = s2;
    BPC->reg1 = gretl_matrix_data_subset_masked(BPC->X1list, Z, t1, t2, 
						mask, &err);
    BPC->reg2 = gretl_matrix_data_subset_masked(BPC->X2list, Z, t1, t2, 
						mask, &err);
    return err;
}

/*
  Records into the container a few items:

  * sample beginning and end
  * the lists
  * the number of parameters
*/

static int bp_base_info (bp_container *BPC, const double **Z, DATAINFO *pdinfo)
{
    int t1, t2, T, i, err = 0;

    t1 = BPC->t1 = pdinfo->t1;
    t2 = BPC->t2 = pdinfo->t2;

    int nelem = BPC->list[0];
    if (nelem < 3) {
	/* we need at least two dep. vars plus one regressor! */
	return E_PARSE;
    }

    BPC->depvar1 = BPC->list[1];
    BPC->depvar2 = BPC->list[2];

    int seppos = gretl_list_separator_position(BPC->list);
    if (seppos == 0) {
	/* same regressors for both equations */
	BPC->k1 = BPC->k2 = nelem - 2;
	BPC->X1list = gretl_list_new(BPC->k1);
	BPC->X2list = gretl_list_new(BPC->k2);
	for (i=1; i<=BPC->k1; i++) {
	    BPC->X1list[i] = BPC->X2list[i] = BPC->list[i+2];
	}
    } else {
	BPC->k1 = seppos - 3;
	BPC->k2 = nelem - seppos;
	BPC->X1list = gretl_list_new(BPC->k1);
	BPC->X2list = gretl_list_new(BPC->k2);
	for (i=1; i<=BPC->k1; i++) {
	    BPC->X1list[i] = BPC->list[i+2];
	}
	for (i=1; i<=BPC->k2; i++) {
	    BPC->X2list[i] = BPC->list[i+seppos];
	}
    }

    int k1, k2;
    k1 = BPC->k1;
    k2 = BPC->k2;

    BPC->npar = k1 + k2 + 1;

    return err;
}

static double match_residuals(int T, double *u1, double *u2, char *mask)
{
    int t;
    double u1t, u2t, v1, v2, cv;

    v1 = 0;
    v2 = 0;
    cv = 0;

    for (t=0; t<T; t++) {
	u1t = u1[t];
	u2t = u2[t];
	mask[t] = (na(u1t) || na(u2t));
	if (!mask[t]) {
	    v1 += u1t * u1t; 
	    v2 += u2t * u2t; 
	    cv += u1t * u2t; 
	}
    }

    double rho = cv / sqrt(v1 * v2);

    if (rho>0.99) {
	rho = 0.99;
    } else if (rho<-0.99) {
	rho = -0.99;
    }

    return rho;
}

/* 
   The following function runs the initial probit models, which 
   serves several purposes: 

   1) the parameter vectors inside the container are initialised
   2) an initial estimate of rho is computed on the basis of
      the correlation between the residuals of the single equation 
      models
   3) Several other quantities are stored in the container (ll0, 
      for one)
   4) the missing values mask is computed by ANDing the residuals;
      they should match anyway because of the preliminary OLS done 
      in run_preliminary_ols(), but hey, checking is always good.

*/

static int biprobit_first_pass(bp_container *BPC, MODEL *olsmod, 
			       double ***pZ, DATAINFO *pdinfo)
{
#if BIPDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
#endif

    MODEL probmod;
    int *list1 = NULL;
    int *list2 = NULL;
    int t1, t2, T, i, err = 0;
    int k1, k2;

    t1 = pdinfo->t1;
    t2 = pdinfo->t2;

    k1 = BPC->k1;
    k2 = BPC->k2;

    T = t2 - t1 + 1;

    /* Initial probit for the first eq. */

#if 0
    gretl_model_init(&probmod);
#endif
    set_reference_missmask_from_model(olsmod);

    list1 = gretl_list_new(k1 + 1);
    list1[1] = BPC->depvar1;
    for (i=1; i<=k1; i++) {
	list1[i+1] = BPC->X1list[i];
    }

    probmod = binary_probit(list1, pZ, pdinfo, OPT_A, NULL);
    if (probmod.errcode) {
	goto bailout;
    }
#if BIPDEBUG
    printmodel(&probmod, pdinfo, OPT_NONE, prn);
#endif

    BPC->ll0 = probmod.lnL;
    BPC->beta = gretl_vector_alloc(k1);
    for (i=0; i<k1; i++) {
	gretl_vector_set(BPC->beta, i, probmod.coeff[i]);
    }

    double *u1;
    u1 = malloc(T * sizeof *u1);
    for (i=0; i<T; i++) {
	u1[i] = probmod.uhat[i];
    }

    /* Initial probit for the second eq. */

    set_reference_missmask_from_model(&probmod);

    list2 = gretl_list_new(k2 + 1);
    list2[1] = BPC->depvar2;
    for (i=1; i<=k2; i++) {
	list2[i+1] = BPC->X2list[i];
    }

    probmod = binary_probit(list2, pZ, pdinfo, OPT_A, NULL);
    if (probmod.errcode) {
	goto bailout;
    }
#if BIPDEBUG
    printmodel(&probmod, pdinfo, OPT_NONE, prn);
#endif

    BPC->ll0 += probmod.lnL;
    BPC->gama = gretl_vector_alloc(k2);
    for (i=0; i<k2; i++) {
	gretl_vector_set(BPC->gama, i, probmod.coeff[i]);
    }

    double *u2;
    u2 = malloc(T * sizeof *u2);
    for (i=0; i<T; i++) {
	u2[i] = probmod.uhat[i];
    }

    BPC->mask = calloc(T, 1);
    double rho = match_residuals(T, u1, u2, BPC->mask);
    BPC->arho = atanh(rho);
    
 bailout:

    clear_model(&probmod);

    free(list1);
    free(list2);
    free(u1);
    free(u2);

    return err;

}

/*
  This function fills up the container: first some base info, 
  then the output from the separate univariate probit models, 
  then the rest (matrices etc)
*/

static int bp_container_fill (bp_container *BPC, MODEL *olsmod, 
			       double ***pZ, DATAINFO *pdinfo)
{
    int t1, t2, T, i, err = 0;
    const double **Z;

    Z = (const double **) *pZ;

    err = bp_base_info(BPC, Z, pdinfo);

    if (!err) {
	err = biprobit_first_pass(BPC, olsmod, pZ, pdinfo);
    }

    if (!err) {
	err = bp_cont_fill_data(BPC, Z);
    }

    if (!err) {
	BPC->fitted1 = gretl_vector_alloc(BPC->ntot);
	BPC->fitted2 = gretl_vector_alloc(BPC->ntot);

	BPC->score = gretl_matrix_alloc(BPC->ntot, BPC->npar);
	BPC->sscore = gretl_vector_alloc(BPC->npar);

	if (BPC->fitted1 == NULL || BPC->fitted1 == NULL ||
	    BPC->score == NULL || BPC->sscore == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	int k1 = BPC->k1;
	int k2 = BPC->k2;

	BPC->H11 = gretl_matrix_alloc(k1, k1);
	BPC->H12 = gretl_matrix_alloc(k1, k2);
	BPC->H13 = gretl_matrix_alloc(k1, 1);
	BPC->H22 = gretl_matrix_alloc(k2, k2);
	BPC->H23 = gretl_matrix_alloc(k2, 1);

	if (BPC->H11 == NULL || BPC->H12 == NULL ||
	    BPC->H13 == NULL || BPC->H22 == NULL || 
	    BPC->H23 == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static int biprobit_init (bp_container *BPC, double ***pZ, DATAINFO *pdinfo, 
			  PRN *prn)
{
    MODEL probmod;
    int *list1 = NULL;
    int *list2 = NULL;
    int i, err = 0;
    int k1, k2;

    k1 = BPC->k1;
    k2 = BPC->k2;

    /* Initial probit for the first eq. */

    gretl_model_init(&probmod);

    list1 = gretl_list_new(k1 + 1);
    list1[1] = BPC->depvar1;
    for (i=1; i<=k1; i++) {
	list1[i+1] = BPC->X1list[i];
    }

    probmod = binary_probit(list1, pZ, pdinfo, OPT_A, NULL);
    if (probmod.errcode) {
	goto bailout;
    }

    BPC->ll0 = probmod.lnL;
    BPC->beta = gretl_vector_alloc(k1);
    for (i=0; i<k1; i++) {
	gretl_vector_set(BPC->beta, i, probmod.coeff[i]);
    }

    clear_model(&probmod);

    /* Initial probit for the second eq. */

    gretl_model_init(&probmod);

    list2 = gretl_list_new(k2 + 1);
    list2[1] = BPC->depvar2;
    for (i=1; i<=k2; i++) {
	list2[i+1] = BPC->X2list[i];
    }

    probmod = binary_probit(list2, pZ, pdinfo, OPT_A, NULL);
    if (probmod.errcode) {
	goto bailout;
    }

    BPC->ll0 += probmod.lnL;
    BPC->gama = gretl_vector_alloc(k2);
    for (i=0; i<k2; i++) {
	gretl_vector_set(BPC->gama, i, probmod.coeff[i]);
    }

 bailout:

    clear_model(&probmod);

    free(list1);
    free(list2);

    return err;
}

int *prelim_list(const int *list)
{
    int i, j, n, m;
    int *invalid;
    n = m = list[0];

    invalid = malloc(n * sizeof *invalid);
    
    for(i=1; i<=n; i++) {
	if(list[i] == LISTSEP) {
	    invalid[i-1] = 1;
	    m--;
	} else {
	    invalid[i-1] = 0;
	    for(j=2; j<i; j++) {
		if (list[i] == list[j]) {
		    invalid[i-1] = 1;
		    m--;
		    break;
		}
	    }
	}
    }

    int *ret;
    ret = gretl_list_new(m);

    j = 1;
    for(i=1; i<=n; i++) {
	if (!invalid[i-1]) {
	    ret[j++] = list[i];
	}
    }

    free(invalid);
    return ret;
}

MODEL run_preliminary_ols (const int *list, double ***pZ, DATAINFO *pdinfo) 
{
#if BIPDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
#endif

    MODEL mod;
    int *tmplist;
    gretl_model_init(&mod);

    tmplist = prelim_list(list);

    mod = lsq(tmplist, (double **) *pZ, pdinfo, OLS, OPT_A);
#if BIPDEBUG
    printmodel(&mod, pdinfo, OPT_NONE, prn);
#endif
    set_reference_missmask_from_model(&mod);

    mod.ci = BIPROBIT;
    return mod;
}

/* ---------------------------------------------------------------------- */
/* likelihood-related functions                                           */
/* ---------------------------------------------------------------------- */

static int biprob_prelim(const double *param, bp_container *bp)
{
    int k1 = bp->k1;
    int k2 = bp->k2;
    int npar = k1 + k2 + 1;
    int err = 0;

    bp->arho = param[npar-1];
    if (fabs(bp->arho) > 16) {
	return 1;
    }

    double rho = tanh(bp->arho);
    int i;

    for (i=0; i<k1; i++) {
	gretl_vector_set(bp->beta, i, param[i]);
    }

    for (i=0; i<k2; i++) {
	gretl_vector_set(bp->gama, i, param[k1+i]);
    }

    err = gretl_matrix_multiply_mod(bp->beta, GRETL_MOD_NONE, 
				    bp->reg1, GRETL_MOD_TRANSPOSE, 
				    bp->fitted1, GRETL_MOD_NONE);

    if (!err) {
    err = gretl_matrix_multiply_mod(bp->gama, GRETL_MOD_NONE, 
				    bp->reg2, GRETL_MOD_TRANSPOSE, 
				    bp->fitted2, GRETL_MOD_NONE);

    }

    return err;

}

static double biprob_loglik (const double *param, void *ptr)
{
    bp_container *bp = (bp_container *) ptr;
    int k1 = bp->k1;
    int k2 = bp->k2;
    int npar = k1 + k2 + 1;
    int i, err;
    double ll = NADBL;

    err = biprob_prelim(param, bp);

    if (err) {
	return ll;
    }

    ll = 0.0;
    int eqt;
    double rho, llt, a, b, P;
    rho = tanh(bp->arho);
    
    for (i=0; i<bp->ntot; i++) {
	a = gretl_vector_get(bp->fitted1, i);
	b = gretl_vector_get(bp->fitted2, i);
	
	a = bp->s1[i] ? a : -a;
	b = bp->s2[i] ? b : -b;
	eqt = (bp->s1[i] == bp->s2[i]);
	
	P = bvnorm_cdf((eqt ? rho : -rho), a, b);
	llt = log(P);
	ll += llt;
    }

    bp->ll = ll;
    return ll;
}

static int biprob_score (double *theta, double *s, int npar, BFGS_CRIT_FUNC ll, 
			 void *ptr)
{
    bp_container *bp = (bp_container *) ptr;
    int k1 = bp->k1;
    int k2 = bp->k2;
    int err;

    err = biprob_prelim(theta, bp);

    if (err) {
	return 0;
    }

    double ca, sa, ssa;

    ca = cosh(bp->arho);
    sa = sinh(bp->arho);
    gretl_matrix_zero(bp->sscore);

    int i, j, eqt;
    double a, b, P, f, d1, d2, da, tmp, u_ab, u_ba;

    for (i=0; i<bp->ntot; i++) {
	a = gretl_vector_get(bp->fitted1, i);
	b = gretl_vector_get(bp->fitted2, i);
	
	a = bp->s1[i] ? a : -a;
	b = bp->s2[i] ? b : -b;
	eqt = (bp->s1[i] == bp->s2[i]);
	ssa = eqt ? sa : -sa;
	
	P = bvnorm_cdf(ssa/ca, a, b);
	
	/* score */
	
	u_ba = (ca*b - ssa*a);
	u_ab = (ca*a - ssa*b);
	f = ca/M_2PI * exp(-0.5* (a*a + u_ba*u_ba));

	d1 = exp(-0.5*a*a) * normal_cdf(u_ba) / (P * SQRT_2_PI);
	d2 = exp(-0.5*b*b) * normal_cdf(u_ab) / (P * SQRT_2_PI);
	da = f / (P*ca*ca);
	da = eqt ? da : -da;

	for(j=0; j<bp->k1; j++) {
	    if (bp->s1[i]) {
		tmp = gretl_matrix_get(bp->reg1, i, j) * d1;
	    } else {
		tmp = -gretl_matrix_get(bp->reg1, i, j) * d1;
	    }
	    
	    gretl_matrix_set(bp->score, i, j, tmp);
	    tmp += gretl_vector_get(bp->sscore, j);
	    gretl_vector_set(bp->sscore, j, tmp);
	}
	
	for(j=0; j<bp->k2; j++) {
	    if (bp->s2[i]) {
		tmp = gretl_matrix_get(bp->reg2, i, j) * d2;
	    } else {
		tmp = -gretl_matrix_get(bp->reg2, i, j) * d2;
	    }

	    gretl_matrix_set(bp->score, i, bp->k1 + j, tmp);
	    tmp += gretl_vector_get(bp->sscore, bp->k1 + j);
	    gretl_vector_set(bp->sscore, bp->k1 + j, tmp);
	}
	
	gretl_matrix_set(bp->score, i, bp->npar - 1, da);
	tmp = da + gretl_vector_get(bp->sscore, bp->npar - 1);
	gretl_vector_set(bp->sscore, bp->npar - 1, tmp);

    }

    if (s != NULL) {
	for (i=0; i<npar; i++) {
	    s[i] = gretl_vector_get(bp->sscore,i);
	}
    }

    return 0;
}

int biprobit_ahessian (double *theta, gretl_matrix *H, void *ptr)
{
    bp_container *bp = (bp_container *) ptr;
    int k1 = bp->k1;
    int k2 = bp->k2;
    int err;

    err = biprob_prelim(theta, bp);

    if (err) {
	return 0;
    }

    double ca, sa, ssa;

    ca = cosh(bp->arho);
    sa = sinh(bp->arho);
    gretl_matrix_zero(bp->sscore);

    int t, i, j, eqt;
    double a, b, P, f, d1, d2, da, tmp, u_ab, u_ba;

    double h11 = 0;
    double h12 = 0;
    double h13 = 0;
    double h22 = 0;
    double h23 = 0;
    double h33 = 0;

    gretl_matrix *H11;
    gretl_matrix *H12;
    gretl_matrix *H13;
    gretl_matrix *H22;
    gretl_matrix *H23;
 
    H11 = bp->H11; 
    H12 = bp->H12; 
    H13 = bp->H13; 
    H22 = bp->H22; 
    H23 = bp->H23; 

    gretl_matrix_zero(H11);
    gretl_matrix_zero(H12);
    gretl_matrix_zero(H13);
    gretl_matrix_zero(H22);
    gretl_matrix_zero(H23);

    err = gretl_matrix_multiply_mod(bp->score, GRETL_MOD_TRANSPOSE,
				    bp->score, GRETL_MOD_NONE,
				    H, GRETL_MOD_NONE);

    double x;
    for (t=0; t<bp->ntot; t++) {
	a = gretl_vector_get(bp->fitted1, t);
	b = gretl_vector_get(bp->fitted2, t);
	
	a = bp->s1[t] ? a : -a;
	b = bp->s2[t] ? b : -b;
	eqt = (bp->s1[t] == bp->s2[t]);
	ssa = eqt ? sa : -sa;
	
	P = bvnorm_cdf(ssa/ca, a, b);
	
	/* score */
	
	u_ba = (ca*b - ssa*a);
	u_ab = (ca*a - ssa*b);
	f = ca/M_2PI * exp(-0.5* (a*a + u_ba*u_ba));

	d1 = exp(-0.5*a*a) * normal_cdf(u_ba) / (P * SQRT_2_PI);
	d2 = exp(-0.5*b*b) * normal_cdf(u_ab) / (P * SQRT_2_PI);
	da = f / (P*ca*ca);
	da = eqt ? da : -da;


	h11 = -(a * d1 + sa*ca*da);

	tmp = f/P;
	h12 = (eqt ? tmp : -tmp);

	tmp = -da * (ca * u_ab);
	h13 = (bp->s1[t] ? tmp : -tmp);

	h22 = -(b*d2 + sa*ca*da);

	tmp = -da * (ca * u_ba);
	h23 = (bp->s2[t] ? tmp : -tmp);

	tmp = (eqt ? u_ba*u_ab : -u_ba*u_ab);
	h33 += da * (ca * tmp - sa);


	for(i=0; i<bp->k1; i++) {
	    x = gretl_matrix_get(bp->reg1, t, i);

	    for(j=i; j<bp->k1; j++) {
		tmp = gretl_matrix_get(H11, i, j);
		tmp += h11 * x * gretl_matrix_get(bp->reg1, t, j);
		gretl_matrix_set(H11, i, j, tmp);
		gretl_matrix_set(H11, j, i, tmp);
	    }

	    for(j=0; j<bp->k2; j++) {
		tmp = gretl_matrix_get(H12, i, j);
		tmp += h12 * x * gretl_matrix_get(bp->reg2, t, j);
		gretl_matrix_set(H12, i, j, tmp);
	    }

	    tmp = gretl_matrix_get(H13, i, 0);
	    tmp += h13 * x;
	    gretl_matrix_set(H13, i, 0, tmp);
	}

	for(i=0; i<bp->k2; i++) {
	    x = gretl_matrix_get(bp->reg2, t, i);

	    for(j=i; j<bp->k2; j++) {
		tmp = gretl_matrix_get(H22, i, j);
		tmp += h22 * x * gretl_matrix_get(bp->reg2, t, j);
		gretl_matrix_set(H22, i, j, tmp);
		gretl_matrix_set(H22, j, i, tmp);
	    }

	    tmp = gretl_matrix_get(H23, i, 0);
	    tmp += h23 * x;
	    gretl_matrix_set(H23, i, 0, tmp);
	}

    }
#if 0
    gretl_matrix_print(H11, "H11");
    gretl_matrix_print(H12, "H12");
    gretl_matrix_print(H13, "H13");
    gretl_matrix_print(H22, "H22");
    gretl_matrix_print(H23, "H23");
    fprintf(stderr, "h33 = %12.6f\n", h33);
    gretl_matrix_print(H, "OPG");
#endif

    for(i=0; i<bp->k1; i++) {
	for(j=i; j<bp->k1; j++) {
	    tmp = gretl_matrix_get(H, i, j) - 
		gretl_matrix_get(H11, i, j);
	    gretl_matrix_set(H, i, j, tmp);
	    gretl_matrix_set(H, j, i, tmp);
	}
	
	for(j=0; j<bp->k2; j++) {
	    tmp = gretl_matrix_get(H, i, j+k1) - 
		gretl_matrix_get(H12, i, j);
	    gretl_matrix_set(H, i, j+k1, tmp);
	    gretl_matrix_set(H, j+k1, i, tmp);
	}
	
	tmp = gretl_matrix_get(H, i, k1+k2) - 
	    gretl_matrix_get(H13, i, 0);
	gretl_matrix_set(H, i, k1+k2, tmp);
	gretl_matrix_set(H, k1+k2, i, tmp);
    }

    for(i=0; i<bp->k2; i++) {
	for(j=i; j<bp->k2; j++) {
	    tmp = gretl_matrix_get(H, i+k1, j+k1) - 
		gretl_matrix_get(H22, i, j);
	    gretl_matrix_set(H, i+k1, j+k1, tmp);
	    gretl_matrix_set(H, j+k1, i+k1, tmp);
	}
	
	tmp = gretl_matrix_get(H, i+k1, k1+k2) - 
	    gretl_matrix_get(H23, i, 0);
	gretl_matrix_set(H, i+k1, k1+k2, tmp);
	gretl_matrix_set(H, k1+k2, i+k1, tmp);
    }

    tmp = gretl_matrix_get(H, k1+k2, k1+k2) - h33;
    gretl_matrix_set(H, k1+k2, k1+k2, tmp);
    gretl_invert_symmetric_matrix(H);


    return 0;
}

#if 0
/* 
   shouldn't be needed anymore, but let's keep it around, 
   just in case
*/

int biprobit_nhessian (double *b, gretl_matrix *H, void *ptr)
{
    bp_container *bpc = (bp_container *) ptr;
    double x, eps = 1.0e-06;
    gretl_matrix *splus = NULL;
    gretl_matrix *sminus = NULL;
    double *theta;
    int n = bpc->npar;
    int m = n*(n+1)/2;
    int i, j, k, err = 0;
    
    theta  = malloc(n * sizeof *theta);
    splus  = gretl_matrix_alloc(1, n);
    sminus = gretl_matrix_alloc(1, n);
    
    if (theta == NULL || H == NULL ||
	splus == NULL || sminus == NULL) {

	err = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<n; i++) {
	theta[i] = b[i];
    }

    for (i=0; i<n; i++) {
	theta[i] += eps;

	biprob_score(theta, NULL, n, NULL, bpc);
	for (j=0; j<n; j++) {
	    x = gretl_vector_get(bpc->sscore, j);
	    gretl_vector_set(splus, j, x);
	}

	theta[i] -= 2*eps;
	biprob_score(theta, NULL, n, NULL, bpc);
	for (j=0; j<n; j++) {
	    x = gretl_vector_get(bpc->sscore, j);
	    gretl_vector_set(sminus, j, x);
	}

	theta[i] += eps;
	for (j=0; j<n; j++) {
	    x = gretl_vector_get(splus, j);
	    x -= gretl_vector_get(sminus, j);
	    gretl_matrix_set(H, i, j, -x/(2*eps));
	}
    }

    gretl_matrix_xtr_symmetric(H);
    gretl_matrix_print(H, "Negative Hessian (not inverted)");
    gretl_invert_symmetric_matrix(H);

 bailout:

    gretl_matrix_free(splus);
    gretl_matrix_free(sminus);
    free(theta);

    return err;
}

#endif

/* ---------------------------------------------------------------------- */
/* Actual likelihood maximization and post-estimation                     */
/* ---------------------------------------------------------------------- */

#define USE_NEWTON 1

static int do_maxlik(bp_container *bp, gretlopt opt, PRN *prn)
{
    int fncount, grcount;
    int np = bp->npar;
    int err = 0;

    /* prova */
    double *theta;
    theta = malloc(np * sizeof *theta);
    int i;
    int k1 = bp->k1;
    int k2 = bp->k2;
    int maxit = 1000;

    for (i=0; i<k1; i++) {
	theta[i] = gretl_vector_get(bp->beta, i);
    }
    for (i=0; i<k2; i++) {
	theta[k1 + i] = gretl_vector_get(bp->gama, i);
    }
    theta[k1 + k2] = bp->arho;

#if USE_NEWTON
    double crittol = 1.0e-06;
    double gradtol = 1.0e-03;
    err = newton_raphson_max (theta, np, maxit, crittol, gradtol,
			      &fncount, C_LOGLIK, biprob_loglik, 
			      biprob_score, biprobit_ahessian,
			      bp, opt & OPT_V, prn);
#else
    double toler = 1.0e-11;
    err = BFGS_max(theta, np, maxit, toler, &fncount, 
		   &grcount, biprob_loglik, C_LOGLIK,
		   biprob_score, bp, NULL,
		   opt & OPT_V, prn);
#endif

    free(theta);
    return err;
}

static int biprobit_vcv(bp_container *bp, gretl_matrix **V, gretlopt opt)
{
    int err;
    int do_hessian = (opt & OPT_R) || (opt & OPT_H);
    int do_opg = !(opt & OPT_H);
    gretl_matrix *Hess = NULL;
    gretl_matrix *OPG = NULL;

    if (do_hessian) {
	double *theta;
	int i, np = bp->npar;
	int k1 = bp->k1;
	int k2 = bp->k2;

	theta = malloc(np * sizeof *theta);

	for (i=0; i<k1; i++) {
	    theta[i] = gretl_vector_get(bp->beta, i);
	}

	for (i=0; i<k2; i++) {
	    theta[k1 + i] = gretl_vector_get(bp->gama, i);
	}
	theta[k1 + k2] = bp->arho;

	int nh = np * (np + 1) / 2;
	Hess   = gretl_matrix_alloc(np, np);

	err = biprobit_ahessian(theta, Hess, bp);
#if 0
	gretl_matrix_print(Hess, "iHess");
#endif
	free(theta);
    }

    if (do_opg) {
	OPG = gretl_matrix_XTX_new(bp->score);
#if 0
	gretl_matrix_print(OPG, "OPG");
#endif
    }

    if (opt & OPT_H) {
	*V = Hess;
    } else if (opt & OPT_R) {
	err = gretl_matrix_qform(Hess, GRETL_MOD_NONE, OPG, 
				 *V, GRETL_MOD_NONE);
#if 0
	gretl_matrix_print(*V, "V");
#endif
	gretl_matrix_free(OPG);
	gretl_matrix_free(Hess);
    } else {
	err = gretl_invert_symmetric_matrix(OPG);
	*V = OPG;
    }

    return err;
}

static int add_indep_LR_test (MODEL *pmod, double LR)
{
    ModelTest *test = model_test_new(GRETL_TEST_INDEP);
    int err = 0;

    if (test != NULL) {
        model_test_set_teststat(test, GRETL_STAT_LR);
        model_test_set_dfn(test, 1);
        model_test_set_value(test, LR);
        model_test_set_pvalue(test, chisq_cdf_comp(1, LR));
        maybe_add_test_to_model(pmod, test);
    } else {
        err = 1;
    }

    return err;
}

static int biprobit_fill_model(MODEL *m, bp_container *b, DATAINFO *pdinfo,
			       gretlopt opt)
{
    double *fullcoeff;
    double *se;
    double x;
    int k1 = b->k1;
    int k2 = b->k2;
    int npar = b->npar - 1;
    int i, j, nvc, err = 0;

    nvc = (npar * npar + npar) / 2;

    fullcoeff = malloc(npar * sizeof *fullcoeff);
    se = malloc(npar * sizeof *fullcoeff);
    m->vcv = malloc(nvc * sizeof *m->vcv);

    if (fullcoeff == NULL || se == NULL || m->vcv == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_model_allocate_params(m, npar);
    }

    if (err) {
	return err;
    }

    for (i=0; i<k1; i++) {
	strcpy(m->params[i], pdinfo->varname[b->X1list[i+1]]);
	fullcoeff[i] = gretl_vector_get(b->beta, i);
    }
	
    for (i=0; i<k2; i++) {
	strcpy(m->params[i+k1], pdinfo->varname[b->X2list[i+1]]);
	fullcoeff[i+k1] = gretl_vector_get(b->gama, i);
    }

    m->lnL = b->ll;
    mle_criteria(m, 0); /* FIXME? */
    m->rho = tanh(b->arho);

    free(m->list);
    m->list = gretl_list_copy(b->list);
    free(m->coeff);
    m->coeff = fullcoeff;
    m->ncoeff = npar;
    m->t1 = b->t1;
    m->t2 = b->t2;
    m->nobs = b->ntot;
    gretl_model_set_coeff_separator(m, pdinfo->varname[b->depvar2], k1);

    gretl_matrix *V;
    V = gretl_matrix_alloc(npar+1, npar+1);
    err = biprobit_vcv(b, &V, opt);

    int k = 0;

    for (i=0; i<npar; i++) {
	for (j=i; j<npar; j++) {
	    x = gretl_matrix_get(V, i, j);
	    m->vcv[k++] = x;
	    if (i==j) {
		se[i] = sqrt(x);
	    }
	}
    }

    m->sderr = se;
    gretl_matrix_free(V);

    err = add_indep_LR_test (m, 2*(b->ll - b->ll0));
    
    return err;
}

/* ---------------------------------------------------------------------- 
   The driver function for the plugin.

   The function biprobit_init estimates the two probit equations 
   separately. Then, do_maxlik performs the estimation of the full
   model.
   ---------------------------------------------------------------------- */

MODEL biprobit_estimate (const int *list, double ***pZ, DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn) 
{
    int err;
    bp_container *bpc;
    MODEL mod;

    bpc = bp_container_new(list);

    mod = run_preliminary_ols(list, pZ, pdinfo);

    if (err) {
	mod.errcode = err;
	return mod;
    } 

    err = bp_container_fill(bpc, &mod, pZ, pdinfo);

    if (err) {
	mod.errcode = err;
	return mod;
    } 

    mod.errcode = do_maxlik(bpc, opt, prn);

    if (!mod.errcode) {
	biprobit_fill_model(&mod, bpc, pdinfo, opt);
    }

    bp_container_destroy(bpc);
    
    return mod;
}
