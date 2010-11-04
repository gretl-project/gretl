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
    double ll0;		     /* log-likelihood for rho = 0 */
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

    gretl_matrix_block *B;   /* workspace for the analytical Hessian */
    gretl_matrix *H11; 
    gretl_matrix *H12;
    gretl_matrix *H13;
    gretl_matrix *H22;
    gretl_matrix *H23;

    gretl_matrix *vcv;	     /* Variance-covariance matrix */
};

/* ---------------------------------------------------------------------- */
/* container-bookkeeping functions                                        */
/* ---------------------------------------------------------------------- */

static void bp_container_destroy (bp_container *bp)
{
    if (bp == NULL) {
	return;
    }

    gretl_matrix_free(bp->score);
    gretl_vector_free(bp->sscore);
    free(bp->mask);

    free(bp->X1list);
    free(bp->X2list);

    free(bp->s1);
    free(bp->s2);

    gretl_matrix_free(bp->reg1);
    gretl_matrix_free(bp->reg2);
    gretl_vector_free(bp->fitted1);
    gretl_vector_free(bp->fitted2);
    gretl_vector_free(bp->u1);
    gretl_vector_free(bp->u2);

    gretl_vector_free(bp->beta);
    gretl_vector_free(bp->gama);
    gretl_matrix_free(bp->vcv);

    gretl_matrix_block_destroy(bp->B);

    free(bp);
}

static bp_container *bp_container_new (const int *list)
{
    bp_container *bp = malloc(sizeof *bp);

    if (bp == NULL) {
	return NULL;
    }

    bp->list = list;

    bp->t1 = bp->t2 = 0;
    bp->ll0 = bp->ll = NADBL;

    bp->score = NULL;
    bp->sscore = NULL;

    bp->mask = NULL;

    bp->X1list = NULL;
    bp->X2list = NULL;
    bp->s1 = NULL;
    bp->s2 = NULL;
    bp->reg1 = NULL;
    bp->reg2 = NULL;
    bp->fitted1 = NULL;
    bp->fitted2 = NULL;
    bp->u1 = NULL;
    bp->u2 = NULL;

    bp->beta = NULL;
    bp->gama = NULL;
    bp->B = NULL;
    bp->arho = 0;

    bp->vcv = NULL;

    return bp;
}

/* ---------------------------------------------------------------------- */
/* Initialization-related functions                                       */
/* ---------------------------------------------------------------------- */

static int bp_cont_fill_data (bp_container *bp, const double **Z)
{
    int t1 = bp->t1;
    int t2 = bp->t2;
    int nmissing = 0;
    int i, T, err = 0;

    T = t2 - t1 + 1;

    if (bp->mask != NULL) {
	for (i=0; i<T; i++) {
	    if (bp->mask[i]) nmissing++;
	}
    } 

    if (T == nmissing) {
	return E_DATA;
    }

    bp->ntot = T - nmissing;

    bp->s1 = malloc(bp->ntot * sizeof *bp->s1);
    bp->s2 = malloc(bp->ntot * sizeof *bp->s2);

    if (bp->s1 == NULL || bp->s2 == NULL) {
	return E_ALLOC;
    }

    if (nmissing == 0) {
	for (i=t1; i<=t2; i++) {
	    bp->s1[i] = (Z[bp->depvar1][i] == 1.0);
	    bp->s2[i] = (Z[bp->depvar2][i] == 1.0);
	}
    } else {
	int k = 0;

	for (i=t1; i<=t2; i++) {
	    if (!bp->mask[i]) {
		bp->s1[k] = (Z[bp->depvar1][i] == 1.0);
		bp->s2[k] = (Z[bp->depvar2][i] == 1.0);
		k++;
	    }
	}
    } 

    bp->reg1 = gretl_matrix_data_subset_masked(bp->X1list, Z, t1, t2, 
					       bp->mask, &err);
    bp->reg2 = gretl_matrix_data_subset_masked(bp->X2list, Z, t1, t2, 
					       bp->mask, &err);

    return err;
}

/*
  Records into the container a few items:

  * sample beginning and end
  * the lists
  * the number of parameters
*/

static int bp_base_info (bp_container *bp, const double **Z, DATAINFO *pdinfo)
{
    int nelem = bp->list[0];
    int i, seppos;
    int err = 0;

    bp->depvar1 = bp->list[1];
    bp->depvar2 = bp->list[2];

    seppos = gretl_list_separator_position(bp->list);

    if (seppos == 0) {
	/* same regressors for both equations */
	bp->k1 = bp->k2 = nelem - 2;
    } else {
	bp->k1 = seppos - 3;
	bp->k2 = nelem - seppos;
    }

    bp->X1list = gretl_list_new(bp->k1);
    bp->X2list = gretl_list_new(bp->k2);

    if (bp->X1list == NULL || bp->X2list == NULL) {
	err = E_ALLOC;
    } else {
	bp->npar = bp->k1 + bp->k2 + 1;
	if (seppos == 0) {
	    for (i=1; i<=bp->k1; i++) {
		bp->X1list[i] = bp->X2list[i] = bp->list[i+2];
	    }
	} else {
	    for (i=1; i<=bp->k1; i++) {
		bp->X1list[i] = bp->list[i+2];
	    }
	    for (i=1; i<=bp->k2; i++) {
		bp->X2list[i] = bp->list[i+seppos];
	    }
	}
    }

    return err;
}

static double match_residuals (int T, double *u1, double *u2, char *mask)
{
    double v1 = 0, v2 = 0, cv = 0;
    double u1t, u2t, rho;
    int t;

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

    rho = cv / sqrt(v1 * v2);

    if (rho > 0.99) {
	rho = 0.99;
    } else if (rho < -0.99) {
	rho = -0.99;
    }

    return rho;
}

static int bp_make_lists (bp_container *bp, int **plist1, int **plist2)
{
    int *list1, *list2;
    int i;

    list1 = gretl_list_new(bp->k1 + 1);
    list2 = gretl_list_new(bp->k2 + 1);

    if (list1 == NULL || list2 == NULL) {
	free(list1);
	free(list2);
	return E_ALLOC;
    }

    list1[1] = bp->depvar1;
    for (i=1; i<=bp->k1; i++) {
	list1[i+1] = bp->X1list[i];
    }

    list2[1] = bp->depvar2;
    for (i=1; i<=bp->k2; i++) {
	list2[i+1] = bp->X2list[i];
    }   

    *plist1 = list1;
    *plist2 = list2;

    return 0;
}

static int bp_add_coeff_vecs_and_mask (bp_container *bp, int T)
{
    bp->beta = gretl_vector_alloc(bp->k1);
    bp->gama = gretl_vector_alloc(bp->k2);
    bp->mask = calloc(T, 1);

    if (bp->beta == NULL || bp->gama == NULL || bp->mask == NULL) {
	return E_ALLOC;
    }

    return 0;
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
      in bp_preliminary_ols(), but hey, checking is always good.

*/

static int biprobit_first_pass (bp_container *bp, MODEL *olsmod, 
				double ***pZ, DATAINFO *pdinfo)
{
#if BIPDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
#endif
    MODEL probmod;
    double rho;
    double *u2, *u1 = NULL;
    int *list1 = NULL, *list2 = NULL;
    int T, i, err;

    err = bp_make_lists(bp, &list1, &list2);
    if (err) {
	return err;
    }

    T = bp->t2 - bp->t1 + 1;

    err = bp_add_coeff_vecs_and_mask(bp, T);
    if (err) {
	goto bailout;
    }

    u1 = malloc(2 * T * sizeof *u1);
    if (u1 == NULL) {
	err = E_ALLOC;
	goto bailout;
    } else {
	u2 = u1 + T;
    }

    /* Initial probit for the first eq. */

    set_reference_missmask_from_model(olsmod);

    probmod = binary_probit(list1, pZ, pdinfo, OPT_A, NULL);
    if (probmod.errcode) {
	goto bailout;
    }

#if BIPDEBUG
    printmodel(&probmod, pdinfo, OPT_NONE, prn);
#endif

    bp->ll0 = probmod.lnL;

    for (i=0; i<bp->k1; i++) {
	gretl_vector_set(bp->beta, i, probmod.coeff[i]);
    }

    for (i=0; i<T; i++) {
	u1[i] = probmod.uhat[probmod.t1 + i];
    }

    /* Initial probit for the second eq. */

    set_reference_missmask_from_model(&probmod);

    clear_model(&probmod);
    probmod = binary_probit(list2, pZ, pdinfo, OPT_A, NULL);
    if (probmod.errcode) {
	goto bailout;
    }

#if BIPDEBUG
    printmodel(&probmod, pdinfo, OPT_NONE, prn);
#endif

    bp->ll0 += probmod.lnL;

    for (i=0; i<bp->k2; i++) {
	gretl_vector_set(bp->gama, i, probmod.coeff[i]);
    }

    for (i=0; i<T; i++) {
	u2[i] = probmod.uhat[probmod.t1 + i];
    }

    rho = match_residuals(T, u1, u2, bp->mask);
    bp->arho = atanh(rho);
    
 bailout:

    clear_model(&probmod);

    free(list1);
    free(list2);
    free(u1);

    return err;
}

/*
  This function fills up the container: first some base info, 
  then the output from the separate univariate probit models, 
  then the rest (matrices, etc.)
*/

static int bp_container_fill (bp_container *bp, MODEL *olsmod, 
			      double ***pZ, DATAINFO *pdinfo)
{
    const double **Z = (const double **) *pZ;
    int err;

    bp->t1 = olsmod->t1;
    bp->t2 = olsmod->t2;

    err = bp_base_info(bp, Z, pdinfo);

    if (!err) {
	err = biprobit_first_pass(bp, olsmod, pZ, pdinfo);
    }

    if (!err) {
	err = bp_cont_fill_data(bp, Z);
    }

    if (!err) {
	bp->fitted1 = gretl_vector_alloc(bp->ntot);
	bp->fitted2 = gretl_vector_alloc(bp->ntot);

	bp->score = gretl_matrix_alloc(bp->ntot, bp->npar);
	bp->sscore = gretl_vector_alloc(bp->npar);

	if (bp->fitted1 == NULL || bp->fitted1 == NULL ||
	    bp->score == NULL || bp->sscore == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	bp->B = gretl_matrix_block_new(&bp->H11, bp->k1, bp->k1,
				       &bp->H12, bp->k1, bp->k2,
				       &bp->H13, bp->k1, 1,
				       &bp->H22, bp->k2, bp->k2,
				       &bp->H23, bp->k2, 1,
				       NULL);

	if (bp->B == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

/* returns a list containing all the distinct elements in the
   input @list, skipping the list separator and any duplicated
   members
*/

static int *prelim_list (const int *list)
{
    int i, j, n, m;
    int *skip;
    int *ret;

    n = m = list[0];

    skip = malloc(n * sizeof *skip);
    if (skip == NULL) {
	return NULL;
    }
    
    for (i=1; i<=n; i++) {
	if (list[i] == LISTSEP) {
	    skip[i-1] = 1;
	    m--;
	} else {
	    skip[i-1] = 0;
	    for (j=2; j<i; j++) {
		if (list[i] == list[j]) {
		    skip[i-1] = 1;
		    m--;
		    break;
		}
	    }
	}
    }

    ret = gretl_list_new(m);

    if (ret != NULL) {
	j = 1;
	for (i=1; i<=n; i++) {
	    if (!skip[i-1]) {
		ret[j++] = list[i];
	    }
	}
    }

    free(skip);

    return ret;
}

/* run an initial OLS to determine the sample that is wanted for
   the two probit models to be estimated subsequently
*/

MODEL bp_preliminary_ols (const int *list, double **Z, DATAINFO *pdinfo) 
{
#if BIPDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
#endif
    MODEL mod;
    int *tmplist;

    tmplist = prelim_list(list);

    if (tmplist == NULL) {
	gretl_model_init(&mod);
	mod.errcode = E_ALLOC;
	return mod;
    }

    mod = lsq(tmplist, Z, pdinfo, OLS, OPT_A);
    free(tmplist);

#if BIPDEBUG
    printmodel(&mod, pdinfo, OPT_NONE, prn);
#endif

    if (!mod.errcode) {
	set_reference_missmask_from_model(&mod);
	mod.ci = BIPROBIT;
    }

    return mod;
}

/* ---------------------------------------------------------------------- */
/* likelihood-related functions                                           */
/* ---------------------------------------------------------------------- */

static int biprob_prelim (const double *theta, bp_container *bp)
{
    int i, err = 0;

    bp->arho = theta[bp->npar - 1];

    if (fabs(bp->arho) > 16) {
	return 1; /* suitable error code? */
    }

    for (i=0; i<bp->k1; i++) {
	gretl_vector_set(bp->beta, i, theta[i]);
    }

    for (i=0; i<bp->k2; i++) {
	gretl_vector_set(bp->gama, i, theta[bp->k1+i]);
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

static double biprob_loglik (const double *theta, void *ptr)
{
    bp_container *bp = (bp_container *) ptr;
    double rho, llt, a, b, P;
    double ll = NADBL;
    int i, eqt, err;

    err = biprob_prelim(theta, bp);

    if (err) {
	return ll;
    }

    ll = 0.0;
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
    double ca, sa, ssa;
    double a, b, P, f, d1, d2, da, tmp, u_ab, u_ba;
    int i, j, eqt;
    int err;

    err = biprob_prelim(theta, bp);

    if (err) {
	return 0;
    }

    ca = cosh(bp->arho);
    sa = sinh(bp->arho);
    gretl_matrix_zero(bp->sscore);

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

	for (j=0; j<bp->k1; j++) {
	    if (bp->s1[i]) {
		tmp = gretl_matrix_get(bp->reg1, i, j) * d1;
	    } else {
		tmp = -gretl_matrix_get(bp->reg1, i, j) * d1;
	    }
	    
	    gretl_matrix_set(bp->score, i, j, tmp);
	    tmp += gretl_vector_get(bp->sscore, j);
	    gretl_vector_set(bp->sscore, j, tmp);
	}
	
	for (j=0; j<bp->k2; j++) {
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

/* analytical hessian */

int biprobit_ahessian (double *theta, gretl_matrix *H, void *ptr)
{
    bp_container *bp = (bp_container *) ptr;
    double a, b, P, f, d1, d2, da, tmp, u_ab, u_ba;
    double ca, sa, ssa, x;
    double h11 = 0;
    double h12 = 0;
    double h13 = 0;
    double h22 = 0;
    double h23 = 0;
    double h33 = 0;
    int k1 = bp->k1;
    int k2 = bp->k2;
    int t, i, j, eqt;
    int err;

    err = biprob_prelim(theta, bp);

    if (err) {
	return 0;
    }

    ca = cosh(bp->arho);
    sa = sinh(bp->arho);
    gretl_matrix_zero(bp->sscore);

    gretl_matrix_zero(bp->H11);
    gretl_matrix_zero(bp->H12);
    gretl_matrix_zero(bp->H13);
    gretl_matrix_zero(bp->H22);
    gretl_matrix_zero(bp->H23);

    err = gretl_matrix_multiply_mod(bp->score, GRETL_MOD_TRANSPOSE,
				    bp->score, GRETL_MOD_NONE,
				    H, GRETL_MOD_NONE);

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

	for (i=0; i<bp->k1; i++) {
	    x = gretl_matrix_get(bp->reg1, t, i);

	    for (j=i; j<bp->k1; j++) {
		tmp = gretl_matrix_get(bp->H11, i, j);
		tmp += h11 * x * gretl_matrix_get(bp->reg1, t, j);
		gretl_matrix_set(bp->H11, i, j, tmp);
		gretl_matrix_set(bp->H11, j, i, tmp);
	    }

	    for (j=0; j<bp->k2; j++) {
		tmp = gretl_matrix_get(bp->H12, i, j);
		tmp += h12 * x * gretl_matrix_get(bp->reg2, t, j);
		gretl_matrix_set(bp->H12, i, j, tmp);
	    }

	    tmp = gretl_matrix_get(bp->H13, i, 0);
	    tmp += h13 * x;
	    gretl_matrix_set(bp->H13, i, 0, tmp);
	}

	for (i=0; i<bp->k2; i++) {
	    x = gretl_matrix_get(bp->reg2, t, i);

	    for (j=i; j<bp->k2; j++) {
		tmp = gretl_matrix_get(bp->H22, i, j);
		tmp += h22 * x * gretl_matrix_get(bp->reg2, t, j);
		gretl_matrix_set(bp->H22, i, j, tmp);
		gretl_matrix_set(bp->H22, j, i, tmp);
	    }

	    tmp = gretl_matrix_get(bp->H23, i, 0);
	    tmp += h23 * x;
	    gretl_matrix_set(bp->H23, i, 0, tmp);
	}

    }

#if 0
    gretl_matrix_print(bp->H11, "H11");
    gretl_matrix_print(bp->H12, "H12");
    gretl_matrix_print(bp->H13, "H13");
    gretl_matrix_print(bp->H22, "H22");
    gretl_matrix_print(bp->H23, "H23");
    fprintf(stderr, "h33 = %12.6f\n", h33);
    gretl_matrix_print(H, "OPG");
#endif

    for (i=0; i<bp->k1; i++) {
	for (j=i; j<bp->k1; j++) {
	    tmp = gretl_matrix_get(H, i, j) - 
		gretl_matrix_get(bp->H11, i, j);
	    gretl_matrix_set(H, i, j, tmp);
	    gretl_matrix_set(H, j, i, tmp);
	}
	
	for (j=0; j<bp->k2; j++) {
	    tmp = gretl_matrix_get(H, i, j+k1) - 
		gretl_matrix_get(bp->H12, i, j);
	    gretl_matrix_set(H, i, j+k1, tmp);
	    gretl_matrix_set(H, j+k1, i, tmp);
	}
	
	tmp = gretl_matrix_get(H, i, k1+k2) - 
	    gretl_matrix_get(bp->H13, i, 0);
	gretl_matrix_set(H, i, k1+k2, tmp);
	gretl_matrix_set(H, k1+k2, i, tmp);
    }

    for (i=0; i<bp->k2; i++) {
	for (j=i; j<bp->k2; j++) {
	    tmp = gretl_matrix_get(H, i+k1, j+k1) - 
		gretl_matrix_get(bp->H22, i, j);
	    gretl_matrix_set(H, i+k1, j+k1, tmp);
	    gretl_matrix_set(H, j+k1, i+k1, tmp);
	}
	
	tmp = gretl_matrix_get(H, i+k1, k1+k2) - 
	    gretl_matrix_get(bp->H23, i, 0);
	gretl_matrix_set(H, i+k1, k1+k2, tmp);
	gretl_matrix_set(H, k1+k2, i+k1, tmp);
    }

    tmp = gretl_matrix_get(H, k1+k2, k1+k2) - h33;
    gretl_matrix_set(H, k1+k2, k1+k2, tmp);

    err = gretl_invert_symmetric_matrix(H);

    return err;
}

#if 0

/* 
   numerical hessian: shouldn't be needed anymore, but let's keep it
   around, just in case 
*/

int biprobit_nhessian (double *b, gretl_matrix *H, void *ptr)
{
    bp_container *bp = (bp_container *) ptr;
    double x, eps = 1.0e-06;
    gretl_matrix *splus = NULL;
    gretl_matrix *sminus = NULL;
    double *theta;
    int n = bp->npar;
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

	biprob_score(theta, NULL, n, NULL, bp);
	for (j=0; j<n; j++) {
	    x = gretl_vector_get(bp->sscore, j);
	    gretl_vector_set(splus, j, x);
	}

	theta[i] -= 2*eps;
	biprob_score(theta, NULL, n, NULL, bp);
	for (j=0; j<n; j++) {
	    x = gretl_vector_get(bp->sscore, j);
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

static double *make_bp_theta (bp_container *bp, int *err)
{
    double *theta = malloc(bp->npar * sizeof *theta);
    int i;

    if (theta == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<bp->k1; i++) {
	theta[i] = gretl_vector_get(bp->beta, i);
    }

    for (i=0; i<bp->k2; i++) {
	theta[bp->k1 + i] = gretl_vector_get(bp->gama, i);
    }

    theta[bp->k1 + bp->k2] = bp->arho;

    return theta;
}

/* ---------------------------------------------------------------------- */
/* Actual likelihood maximization and post-estimation                     */
/* ---------------------------------------------------------------------- */

#define USE_NEWTON 1

static int bp_do_maxlik (bp_container *bp, gretlopt opt, PRN *prn)
{
    double crittol = 1.0e-06;
    double gradtol = 1.0e-03;
    int fncount, maxit = 1000;
    double *theta;
    int err = 0;

    theta = make_bp_theta(bp, &err);
    if (err) {
	return err;
    }

#if USE_NEWTON
    err = newton_raphson_max(theta, bp->npar, maxit, crittol, gradtol,
			     &fncount, C_LOGLIK, biprob_loglik, 
			     biprob_score, biprobit_ahessian,
			     bp, opt & OPT_V, prn);
#else
    int grcount;

    crittol = 1.0e-11;
    err = BFGS_max(theta, bp->npar, maxit, crittol, &fncount, 
		   &grcount, biprob_loglik, C_LOGLIK,
		   biprob_score, bp, NULL,
		   opt & OPT_V, prn);
#endif

    free(theta);

    return err;
}

static gretl_matrix *biprobit_vcv (bp_container *bp, gretlopt opt, int *err)
{
    int do_hessian = (opt & OPT_R) || (opt & OPT_H);
    int do_opg = !(opt & OPT_H);
    gretl_matrix *Hess = NULL;
    gretl_matrix *OPG = NULL;
    gretl_matrix *V = NULL;

    if (do_hessian) {
	double *theta = make_bp_theta(bp, err);

	if (!*err) {
	    Hess = gretl_matrix_alloc(bp->npar, bp->npar);
	    if (Hess == NULL) {
		*err = E_ALLOC;
	    }
	}

	if (!*err) {
	    *err = biprobit_ahessian(theta, Hess, bp);
	}
#if 0
	gretl_matrix_print(Hess, "iHess");
#endif
	free(theta);
    }

    if (!*err && do_opg) {
	OPG = gretl_matrix_XTX_new(bp->score);
	if (OPG == NULL) {
	    *err = E_ALLOC;
	}
#if 0
	gretl_matrix_print(OPG, "OPG");
#endif
    }

    if (!*err) {
	if (opt & OPT_H) {
	    V = Hess;
	    Hess = NULL;
	} else if (opt & OPT_R) {
	    V = gretl_matrix_alloc(bp->npar, bp->npar);
	    if (V == NULL) {
		*err = E_ALLOC;
	    } else {
		*err = gretl_matrix_qform(Hess, GRETL_MOD_NONE, OPG, 
					  V, GRETL_MOD_NONE);
	    }
#if 0
	    gretl_matrix_print(V, "V");
#endif
	} else {
	    *err = gretl_invert_symmetric_matrix(OPG);
	    if (!*err) {
		V = OPG;
		OPG = NULL;
	    }
	}
    }

    gretl_matrix_free(OPG);
    gretl_matrix_free(Hess);

    return V;
}

static int add_indep_LR_test (MODEL *pmod, double LR)
{
    ModelTest *test = model_test_new(GRETL_TEST_INDEP);
    int err = 0;

    if (test == NULL) {
	err = E_ALLOC;
    } else {
        model_test_set_teststat(test, GRETL_STAT_LR);
        model_test_set_dfn(test, 1);
        model_test_set_value(test, LR);
        model_test_set_pvalue(test, chisq_cdf_comp(1, LR));
        maybe_add_test_to_model(pmod, test);
    } 

    return err;
}

static int biprobit_fill_model (MODEL *pmod, bp_container *bp, DATAINFO *pdinfo,
				gretlopt opt)
{
    gretl_matrix *V;
    double *coeff, *sderr, *vcv;
    int nvc, npar = bp->npar - 1;
    int i, j, err = 0;

    nvc = (npar * npar + npar) / 2;

    coeff = realloc(pmod->coeff, npar * sizeof *coeff);
    sderr = realloc(pmod->sderr, npar * sizeof *sderr);
    vcv = realloc(pmod->vcv, nvc * sizeof *vcv);

    if (coeff == NULL || sderr == NULL || vcv == NULL) {
	err = E_ALLOC;
    } else {
	pmod->coeff = coeff;
	pmod->sderr = sderr;
	pmod->vcv = vcv;
	pmod->ncoeff = npar;
    }

    if (!err) {
	err = gretl_model_allocate_params(pmod, npar);
    }

    if (err) {
	return err;
    }

    for (i=0; i<bp->k1; i++) {
	strcpy(pmod->params[i], pdinfo->varname[bp->X1list[i+1]]);
	pmod->coeff[i] = gretl_vector_get(bp->beta, i);
    }
	
    for (i=0; i<bp->k2; i++) {
	strcpy(pmod->params[i+bp->k1], pdinfo->varname[bp->X2list[i+1]]);
	pmod->coeff[i + bp->k1] = gretl_vector_get(bp->gama, i);
    }

    pmod->lnL = bp->ll;
    mle_criteria(pmod, 0); /* FIXME? */
    pmod->rho = tanh(bp->arho);

    free(pmod->list);
    pmod->list = gretl_list_copy(bp->list);

    pmod->t1 = bp->t1;
    pmod->t2 = bp->t2;
    pmod->nobs = bp->ntot;

    gretl_model_set_coeff_separator(pmod, pdinfo->varname[bp->depvar2], bp->k1);

    V = biprobit_vcv(bp, opt, &err);

    if (!err) {
	double x;
	int k = 0;

	for (i=0; i<npar; i++) {
	    for (j=i; j<npar; j++) {
		x = gretl_matrix_get(V, i, j);
		pmod->vcv[k++] = x;
		if (i == j) {
		    pmod->sderr[i] = sqrt(x);
		}
	    }
	}
    }

    gretl_matrix_free(V);

    err = add_indep_LR_test(pmod, 2 * (bp->ll - bp->ll0));
    
    return err;
}

/* ---------------------------------------------------------------------- 
   The driver function for the plugin.

   The function biprobit_init estimates the two probit equations 
   separately. Then, bp_do_maxlik performs the estimation of the full
   model.
   ---------------------------------------------------------------------- */

MODEL biprobit_estimate (const int *list, double ***pZ, DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn) 
{
    bp_container *bp;
    MODEL mod;
    int err = 0;

    if (list[0] < 3) {
	/* we need at least two dep. vars plus one regressor! */
	err = E_PARSE;
    } else {
	bp = bp_container_new(list);
	if (bp == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	gretl_model_init(&mod);
	mod.errcode = err;
	return mod;
    }

    mod = bp_preliminary_ols(list, *pZ, pdinfo);
    if (mod.errcode) {
	return mod;
    } 

    mod.errcode = bp_container_fill(bp, &mod, pZ, pdinfo);
    if (mod.errcode) {
	return mod;
    } 

    mod.errcode = bp_do_maxlik(bp, opt, prn);

    if (!mod.errcode) {
	biprobit_fill_model(&mod, bp, pdinfo, opt);
    }

    bp_container_destroy(bp);
    
    return mod;
}
