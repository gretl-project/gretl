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
#include "version.h"
#include "matrix_extra.h"
#include "gretl_bfgs.h"
#include "libset.h"
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

    int nobs; 	             /* number of usable observations */
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
    gretl_vector *fitted1;   /* x_1'\beta */
    gretl_vector *fitted2;   /* x_2'\gamma */

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

static int bp_cont_fill_data (bp_container *bp, DATASET *dset)
{
    int nmissing = 0;
    int i, t, T, err = 0;

    T = bp->t2 - bp->t1 + 1;

    for (i=0; i<T; i++) {
	if (bp->mask[i]) nmissing++;
    }

    if (T == nmissing) {
	return E_DATA;
    }

    bp->nobs = T - nmissing;

    bp->s1 = malloc(bp->nobs * sizeof *bp->s1);
    bp->s2 = malloc(bp->nobs * sizeof *bp->s2);

    if (bp->s1 == NULL || bp->s2 == NULL) {
	return E_ALLOC;
    }

    i = 0;
    for (t=bp->t1; t<=bp->t2; t++) {
	if (!bp->mask[t - bp->t1]) {
	    bp->s1[i] = (dset->Z[bp->depvar1][t] == 1.0);
	    bp->s2[i] = (dset->Z[bp->depvar2][t] == 1.0);
	    i++;
	}
    }	    

    bp->reg1 = gretl_matrix_data_subset_masked(bp->X1list, dset,
					       bp->t1, bp->t2, 
					       bp->mask, &err);
    bp->reg2 = gretl_matrix_data_subset_masked(bp->X2list, dset,
					       bp->t1, bp->t2, 
					       bp->mask, &err);

    return err;
}

/*
  Records into the container a few items:

  * the lists
  * the number of parameters
*/

static int bp_base_info (bp_container *bp, DATASET *dset)
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
				DATASET *dset, PRN *prn)
{
    MODEL probmod;
    double *u2, *u1 = NULL;
    int *list1 = NULL, *list2 = NULL;
    gretlopt popt = OPT_A | OPT_P;
    int T, i, err;

    err = bp_make_lists(bp, &list1, &list2);
    if (err) {
	return err;
    }

    T = bp->t2 - bp->t1 + 1;

    gretl_model_init(&probmod, NULL); /* in case of failure below */

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

    if (prn != NULL) {
	popt |= OPT_V;
	pputc(prn, '\n');
    }

    /* Initial probit for the first equation */

    set_reference_missmask_from_model(olsmod);

    probmod = binary_probit(list1, dset, popt, prn);
    if (probmod.errcode) {
	goto bailout;
    }

    if (prn != NULL) {
	probmod.aux = AUX_BIPROB;
	printmodel(&probmod, dset, OPT_NONE, prn);
    }

    bp->ll0 = probmod.lnL;

    for (i=0; i<bp->k1; i++) {
	gretl_vector_set(bp->beta, i, probmod.coeff[i]);
    }

    for (i=0; i<T; i++) {
	u1[i] = probmod.uhat[probmod.t1 + i];
    }

    /* Initial probit for the second equation */

    set_reference_missmask_from_model(&probmod);

    clear_model(&probmod);
    probmod = binary_probit(list2, dset, popt, prn);
    if (probmod.errcode) {
	goto bailout;
    }

    if (prn != NULL) {
	probmod.aux = AUX_BIPROB;
	printmodel(&probmod, dset, OPT_NONE, prn);
    }

    bp->ll0 += probmod.lnL;

    for (i=0; i<bp->k2; i++) {
	gretl_vector_set(bp->gama, i, probmod.coeff[i]);
    }

    for (i=0; i<T; i++) {
	u2[i] = probmod.uhat[probmod.t1 + i];
	bp->mask[i] = (na(u1[i]) || na(u2[i]));
    }

    bp->arho = 0;
    
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
			      DATASET *dset, PRN *prn)
{
    int err;

    /* sample beginning and end */
    bp->t1 = olsmod->t1;
    bp->t2 = olsmod->t2;

    /* lists and number of parameters per equation */
    err = bp_base_info(bp, dset);

    if (!err) {
	err = biprobit_first_pass(bp, olsmod, dset, prn);
    }

    if (!err) {
	err = bp_cont_fill_data(bp, dset);
    }

    if (!err) {
	bp->fitted1 = gretl_vector_alloc(bp->nobs);
	bp->fitted2 = gretl_vector_alloc(bp->nobs);

	bp->score = gretl_matrix_alloc(bp->nobs, bp->npar);
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
    int m = gretl_list_n_distinct_members(list);
    int *ret;

    ret = gretl_list_new(m);

    if (ret != NULL) {
	int i, j, k = 1;
	int unique;

	for (i=1; i<=list[0]; i++) {
	    if (list[i] != LISTSEP) {
		unique = 1;
		for (j=2; j<i; j++) {
		    if (list[i] == list[j]) {
			unique = 0;
			break;
		    }
		}
		if (unique) {
		    ret[k++] = list[i];
		}
	    }
	}
    }

    return ret;
}

/* run an initial OLS to determine the sample that is wanted for
   the two probit models to be estimated subsequently
*/

MODEL bp_preliminary_ols (const int *list, DATASET *dset) 
{
#if BIPDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
#endif
    MODEL mod;
    int *tmplist;

    tmplist = prelim_list(list);

    if (tmplist == NULL) {
	gretl_model_init(&mod, NULL);
	mod.errcode = E_ALLOC;
	return mod;
    }

    mod = lsq(tmplist, dset, OLS, OPT_A);
    if (gretl_model_get_data(&mod, "droplist")) {
	gretl_model_destroy_data_item(&mod, "droplist");
    }
    free(tmplist);

#if BIPDEBUG
    printmodel(&mod, dset, OPT_NONE, prn);
    gretl_print_destroy(prn);
#endif

    if (!mod.errcode) {
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

    if (fabs(bp->arho) > 18) {
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
    
    for (i=0; i<bp->nobs; i++) {
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
	return err;
    }

    ca = cosh(bp->arho);
    sa = sinh(bp->arho);
    gretl_matrix_zero(bp->sscore);

    for (i=0; i<bp->nobs; i++) {
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

	d1 = bp->s1[i] ? d1 : -d1;
	d2 = bp->s2[i] ? d2 : -d2;
	da = eqt ? da : -da;

	for (j=0; j<bp->k1; j++) {
	    tmp = gretl_matrix_get(bp->reg1, i, j) * d1;
	    gretl_matrix_set(bp->score, i, j, tmp);
	    tmp += gretl_vector_get(bp->sscore, j);
	    gretl_vector_set(bp->sscore, j, tmp);
	}
	
	for (j=0; j<bp->k2; j++) {
	    tmp = gretl_matrix_get(bp->reg2, i, j) * d2;
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

/* 
   Analytical Hessian for the biprobit model; it is assumed that 
   the score matrix has already been computed and is safely stored 
   in the container. This is important because we use a neat trick:

   \frac{\partial^2 ln f(x)}{\partial x \partial x'} =
   \frac{\partial^2 f(x)}{\partial x \partial x'} - gg'

   where g = \frac{\partial ln f(x)}{\partial x}.
   In practice, we just compute explicitly the difference between
   the OPG and the actual Hessian, and then reconstruct H from 
   there.
 */

int biprobit_hessian (double *theta, gretl_matrix *H, void *ptr)
{
    bp_container *bp = (bp_container *) ptr;
    double a, b, P, f, d1, d2, da, tmp, u_ab, u_ba;
    double ca, sa, ssa, rho, x;
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

    if (!err) {
	ca = cosh(bp->arho);
	sa = sinh(bp->arho);

	gretl_matrix_zero(bp->sscore);
	gretl_matrix_block_zero(bp->B);

	/* first, put the OPG matrix into H */

	err = gretl_matrix_multiply_mod(bp->score, GRETL_MOD_TRANSPOSE,
					bp->score, GRETL_MOD_NONE,
					H, GRETL_MOD_NONE);
    }

    if (err) {
	return err;
    }

    for (t=0; t<bp->nobs; t++) {

	/* Warning: it is important to handle properly the four
	   cases (00, 01, 10, 11) and the associated sign switches;
	   this may be very tricky. The code below should be bug-free
	   from this point of view, but you never know. Kudos to
	   Claudia Pigini for going through this with incredible 
	   perseverance.
	*/

	a = gretl_vector_get(bp->fitted1, t);
	b = gretl_vector_get(bp->fitted2, t);
	
	a = bp->s1[t] ? a : -a;
	b = bp->s2[t] ? b : -b;
	eqt = (bp->s1[t] == bp->s2[t]);
	ssa = eqt ? sa : -sa;
	rho = ssa/ca;
	P = bvnorm_cdf(rho, a, b);
	
	/* score (for atan(rho) we use the precomputed one) */
	
	u_ba = (ca*b - ssa*a);
	u_ab = (ca*a - ssa*b);

	d1 = exp(-0.5*a*a) * normal_cdf(u_ba) / (P * SQRT_2_PI);
	d2 = exp(-0.5*b*b) * normal_cdf(u_ab) / (P * SQRT_2_PI);       

	f = ca/M_2PI * exp(-0.5* (a*a + u_ba*u_ba));

	da = gretl_matrix_get(bp->score,t, k1+k2);

	/* hessian + opg */

	h11 = -(a * d1 + sa*ca*da);
	h12 = (eqt ? f/P : -f/P);
	h13 = (bp->s1[t] ? -da * (ca * u_ab) : da * (ca * u_ab));
	h22 = -(b*d2 + sa*ca*da);
	h23 = (bp->s2[t] ? -da * (ca * u_ba) : da * (ca * u_ba));

	tmp = (eqt ? u_ba*u_ab : -u_ba*u_ab);
	tmp = da * (ca * tmp - sa) /ca;
	h33 += tmp;

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
    fprintf(stderr, "ca = %g, sa = %g, h33 = %12.6f\n", ca, sa, h33);
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

    return err;
}

static gretl_matrix *biprobit_hessian_inverse (double *theta, 
					       void *ptr,
					       int *err)
{
    bp_container *bp = ptr;
    gretl_matrix *H;

    H = gretl_matrix_alloc(bp->npar, bp->npar);
    
    if (H == NULL) {
	*err = E_ALLOC;
    } else {
	*err = biprobit_hessian(theta, H, ptr);
    }

    if (!*err) {
	*err = gretl_invert_symmetric_matrix(H);
    }

    return H;
}

#if 0

/* 
   numerical hessian: shouldn't be needed, but let's keep it
   around, just in case 
*/

int biprobit_numerical_hessian (double *b, gretl_matrix *H, void *ptr)
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

static int bp_do_maxlik (bp_container *bp, gretlopt opt, PRN *prn)
{
    double crittol = 1.0e-06;
    double gradtol = 1.0e-05;
    int fncount, maxit = 1000;
    gretlopt maxopt = opt & OPT_V;
    int quiet = opt & OPT_Q;
    double *theta;
    int err = 0;

    theta = make_bp_theta(bp, &err);
    if (err) {
	return err;
    }

    if (libset_get_int(GRETL_OPTIM) == OPTIM_BFGS) {
	int grcount = 0;

	crittol = 1.0e-11;
	err = BFGS_max(theta, bp->npar, maxit, crittol, &fncount, 
		       &grcount, biprob_loglik, C_LOGLIK,
		       biprob_score, bp, NULL,
		       maxopt, quiet ? NULL : prn);
    } else {	
	err = newton_raphson_max(theta, bp->npar, maxit, crittol, gradtol,
				 &fncount, C_LOGLIK, biprob_loglik, 
				 biprob_score, biprobit_hessian,
				 bp, maxopt, quiet ? NULL : prn);
    }

    free(theta);

    return err;
}

static int biprobit_vcv (MODEL *pmod, bp_container *bp, 
			 const DATASET *dset, gretlopt opt)
{
    gretl_matrix *H = NULL;
    int err = 0;

    if (opt & OPT_G) {
	err = gretl_model_add_OPG_vcv(pmod, bp->score);
    } else {
	double *theta = make_bp_theta(bp, &err);

	if (!err) {
	    H = biprobit_hessian_inverse(theta, bp, &err);
	}
	if (!err) {
	    if (opt & OPT_R) {
		err = gretl_model_add_QML_vcv(pmod, BIPROBIT, 
					      H, bp->score,
					      dset, opt);
	    } else {
		err = gretl_model_add_hessian_vcv(pmod, H);
	    }
	}
	free(theta);
    }

    gretl_matrix_free(H);

    return err;
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

static int bp_add_hat_matrices (MODEL *pmod, bp_container *bp,
				gretlopt opt)
{
    gretl_matrix *Uh, *Yh;
    int T = bp->t2 - bp->t1 + 1;
    int Yhcols = (opt & OPT_X)? 2: 4;
    int err = 0;

    Uh = gretl_matrix_alloc(T, 2);
    Yh = gretl_matrix_alloc(T, Yhcols);

    if (Uh == NULL || Yh == NULL) {
	gretl_matrix_free(Uh);
	gretl_matrix_free(Yh);
	err = E_ALLOC;
    } else {
	double rho, ca, sa, ssa;
	double a, b, u_ab, u_ba, P, im;
	double p11, p10, p01, p00;
	int t, i = 0;

	ca = cosh(bp->arho);
	sa = sinh(bp->arho);
	rho = sa/ca;

	for (t=0; t<T; t++) {

	    if (na(pmod->uhat[t + pmod->t1])) {
		gretl_matrix_set(Uh, t, 0, M_NA);
		gretl_matrix_set(Uh, t, 1, M_NA);
		gretl_matrix_set(Yh, t, 0, M_NA);
		gretl_matrix_set(Yh, t, 1, M_NA);
		if (Yhcols == 4) {
		    gretl_matrix_set(Yh, t, 2, M_NA);
		    gretl_matrix_set(Yh, t, 3, M_NA);
		}		    
		continue;
	    }

	    a = gretl_vector_get(bp->fitted1, i);
	    b = gretl_vector_get(bp->fitted2, i);

	    if (opt & OPT_X) {
		/* Yh should contain the index function values --
		   corresponds to the flag --save-xbeta */
		gretl_matrix_set(Yh, t, 0, a);
		gretl_matrix_set(Yh, t, 1, b);
	    } 

	    p11 = bvnorm_cdf( rho,  a,  b);
	    p10 = bvnorm_cdf(-rho,  a, -b);
	    p01 = bvnorm_cdf(-rho, -a,  b);
	    p00 = bvnorm_cdf( rho, -a, -b);

	    if (!(opt & OPT_X)) {
		/* Let the columns of Yh hold the estimated probabilities
		   of the outcomes (1,1), (1,0), (0,1) and (0,0)
		   respectively.
		*/
		gretl_matrix_set(Yh, t, 0, p11);
		gretl_matrix_set(Yh, t, 1, p10);
		gretl_matrix_set(Yh, t, 2, p01);
		gretl_matrix_set(Yh, t, 3, p00);
	    }

	    /* Put the generalized residuals into uhat;
	       The generalized residuals are defined as

	       E(u_1| y_1, y_2, x_1, x_2)
	       E(u_2| y_1, y_2, x_1, x_2)

	       They should be equal to the score for the constant
	       (if any)
	    */

	    if (bp->s1[i] && bp->s2[i]) {
		P = p11;
	    } else if (bp->s1[i] && !bp->s2[i]) {
		P = p10;
	    } else if (!bp->s1[i] && bp->s2[i]) {
		P = p01;
	    } else {
		P = p00;
	    }

	    a = bp->s1[i] ? a : -a;
	    b = bp->s2[i] ? b : -b;
	    ssa = (bp->s1[i] == bp->s2[i]) ? sa : -sa;
	    u_ba = (ca*b - ssa*a);
	    u_ab = (ca*a - ssa*b);

	    im = exp(-0.5*a*a) * normal_cdf(u_ba) / (P * SQRT_2_PI);
	    gretl_matrix_set(Uh, t, 0, bp->s1[i] ? im : -im);

	    im = exp(-0.5*b*b) * normal_cdf(u_ab) / (P * SQRT_2_PI);
	    gretl_matrix_set(Uh, t, 1, bp->s2[i] ? im : -im);
	    
	    i++;
	}

	gretl_matrix_set_t1(Uh, bp->t1);
	gretl_matrix_set_t2(Uh, bp->t2);
	gretl_matrix_set_t1(Yh, bp->t1);
	gretl_matrix_set_t2(Yh, bp->t2);

	gretl_model_set_matrix_as_data(pmod, "bp_uhat", Uh);
	gretl_model_set_matrix_as_data(pmod, "bp_yhat", Yh);
    }

    return err;
}

static int biprobit_fill_model (MODEL *pmod, bp_container *bp, 
				const DATASET *dset,
				gretlopt opt)
{
    double *tmp;
    int npar = bp->npar - 1;
    int i, xi, err = 0;

    tmp = realloc(pmod->coeff, npar * sizeof *tmp);

    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	pmod->coeff = tmp;
	pmod->ncoeff = npar;
	err = gretl_model_allocate_param_names(pmod, npar);
    }

    if (err) {
	return err;
    }

    for (i=0; i<bp->k1; i++) {
	xi = bp->X1list[i+1];
	gretl_model_set_param_name(pmod, i, dset->varname[xi]);
	pmod->coeff[i] = gretl_vector_get(bp->beta, i);
    }
	
    for (i=0; i<bp->k2; i++) {
	xi = bp->X2list[i+1];
	gretl_model_set_param_name(pmod, i+bp->k1, dset->varname[xi]);
	pmod->coeff[i + bp->k1] = gretl_vector_get(bp->gama, i);
    }

    /* scrub some invalid statistics */
    pmod->rsq = pmod->adjrsq = NADBL;
    pmod->ybar = pmod->sdy = NADBL;
    pmod->fstt = pmod->sigma = NADBL;

    free(pmod->yhat);
    pmod->yhat = NULL;

    pmod->lnL = bp->ll;
    mle_criteria(pmod, 1); /* allow for estimation of rho */
    pmod->rho = tanh(bp->arho);

    free(pmod->list);
    pmod->list = gretl_list_copy(bp->list);

    pmod->t1 = bp->t1;
    pmod->t2 = bp->t2;
    pmod->nobs = bp->nobs;

    gretl_model_set_coeff_separator(pmod, dset->varname[bp->depvar2], bp->k1);

    err = biprobit_vcv(pmod, bp, dset, opt);

    if (!err) {
	err = bp_add_hat_matrices(pmod, bp, opt);
    }

    if (!err) {
	err = add_indep_LR_test(pmod, 2 * (bp->ll - bp->ll0));
    }
    
    return err;
}

/* 
   The driver function for the plugin.
*/

MODEL biprobit_estimate (const int *list, DATASET *dset, 
			 gretlopt opt, PRN *prn) 
{
    PRN *vprn = NULL;
    bp_container *bp;
    MODEL mod;
    int err = 0;

    if (list[0] < 3 || (gretl_list_has_separator(list) && list[0] < 5)) {
	/* we need at least two dep. vars plus one regressor */
	err = E_ARGS;
    } else {
	bp = bp_container_new(list);
	if (bp == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	gretl_model_init(&mod, NULL);
	mod.errcode = err;
	return mod;
    }

    mod = bp_preliminary_ols(list, dset);
    if (mod.errcode) {
	return mod;
    } 

    if (opt & OPT_C) {
	/* cluster implies robust */
	opt |= OPT_R;
    }

    if (opt & OPT_V) {
	vprn = prn;
    }

    mod.errcode = bp_container_fill(bp, &mod, dset, vprn);
    if (mod.errcode) {
	return mod;
    } 

    mod.errcode = bp_do_maxlik(bp, opt, prn);

    if (!mod.errcode) {
	biprobit_fill_model(&mod, bp, dset, opt);
    }

    bp_container_destroy(bp);
    
    return mod;
}
