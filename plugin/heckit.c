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

#include "libgretl.h"
#include "matrix_extra.h"
#include "libset.h"
#include "missing_private.h"
#include "gretl_bfgs.h"

#define HDEBUG 0

typedef struct h_container_ h_container;

struct h_container_ {
    const int *list;         /* incoming model specification */
    int t1, t2;              /* start and end of sample range */
    int kmain;		     /* no. of params in the main eq. */
    int ksel;		     /* no. of params in the selection eq. */
    double ll;		     /* log-likelihood */

    int ntot, nunc;	     /* total and uncensored obs */
    int depvar;		     /* location of y in array Z */
    int selvar;		     /* location of selection variable in array Z */
    int millsvar;            /* location of Mills ratios in array Z */
    int *Xlist;		     /* regressor list for the main eq. */
    int *Zlist;		     /* regressor list for the selection eq. */

    gretl_matrix *y;	     /* dependent var */
    gretl_matrix *reg;	     /* main eq. regressors */
    gretl_matrix *mills;     /* Mills ratios from selection eq */
    gretl_matrix *delta;     /* used in 2-step vcv calculations */
    gretl_matrix *d;	     /* selection dummy variable */
    gretl_matrix *selreg;    /* selection eq. regressors */ 
    gretl_matrix *selreg_u;  /* selection eq. regressors (subsample d==1) */
    gretl_vector *fitted;
    gretl_vector *u;
    gretl_vector *ndx;

    gretl_vector *beta;	     /* main eq. parameters */
    gretl_vector *gama;	     /* selection eq. parameters */
    double sigma;
    double rho;
    double lambda;	     /* rho*sigma by definition */

    gretl_matrix *vcv;	     /* Variance-covariance matrix */
    gretl_matrix *VProbit;   /* 1st stage probit covariance matrix */

    char *probmask;	     /* mask NAs for initial probit */
    char *fullmask;	     /* mask NAs */
    char *uncmask;	     /* mask NAs and (d==0) */
};

static void h_container_destroy (h_container *HC)
{
    if (HC == NULL) {
	return;
    }

    free(HC->Xlist);
    free(HC->Zlist);

    gretl_matrix_free(HC->y);
    gretl_matrix_free(HC->reg);
    gretl_matrix_free(HC->mills);
    gretl_matrix_free(HC->delta);
    gretl_matrix_free(HC->d);
    gretl_matrix_free(HC->selreg);
    gretl_matrix_free(HC->selreg_u);
    gretl_vector_free(HC->fitted);
    gretl_vector_free(HC->u);
    gretl_vector_free(HC->ndx);

    gretl_vector_free(HC->beta);
    gretl_vector_free(HC->gama);

    gretl_matrix_free(HC->vcv);
    gretl_matrix_free(HC->VProbit);

    free(HC->probmask);
    free(HC->fullmask);
    free(HC->uncmask);

    free(HC);
}

static h_container *h_container_new (const int *list)
{
    h_container *HC = malloc(sizeof *HC);

    if (HC == NULL) {
	return NULL;
    }

    HC->list = list;
    HC->t1 = HC->t2 = 0;
    HC->ll = NADBL;

    HC->Xlist = NULL;
    HC->Zlist = NULL;

    HC->y = NULL;
    HC->reg = NULL;
    HC->mills = NULL;
    HC->delta = NULL;
    HC->d = NULL;
    HC->selreg = NULL;
    HC->selreg_u = NULL;
    HC->fitted = NULL;
    HC->u = NULL;
    HC->ndx = NULL;

    HC->beta = NULL;
    HC->gama = NULL;
    HC->sigma = NADBL;
    HC->rho = NADBL;
    HC->lambda = NADBL;

    HC->vcv = NULL;
    HC->VProbit = NULL;

    HC->probmask = NULL;
    HC->fullmask = NULL;
    HC->uncmask = NULL;

    return HC;
}

/* The following is complicated: for the purposes of the Heckit model
   as such, we need to mask out any observations that (a) are
   "selected" and (b) have missing values for any variables in Xlist.
   With regard to the initial probit, we have to mask in addition any
   observations that have missing values for any of the variables in
   Zlist.
*/

static int 
make_heckit_NA_mask (h_container *HC, const int *reglist, const int *sellist,
		     const double **Z, const DATAINFO *pdinfo)
{
    int T = sample_size(pdinfo);
    int totmiss = 0;
    int i, t, err = 0;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	int miss = 0;

	/* screen the primary equation variables */
	if (Z[HC->selvar][t] == 1.0) {
	    for (i=1; i<=reglist[0] && !miss; i++) {
		if (na(Z[reglist[i]][t])) {
		    miss = 1;
		    totmiss++;
		}
	    }
	}

	/* screen the selection equation variables */
	for (i=1; i<=sellist[0] && !miss; i++) {
	    if (na(Z[sellist[i]][t])) {
		miss = 1;
		totmiss++;
	    }
	}
    }

    if (totmiss == T) {
	return E_MISSDATA;
    }

    if (totmiss > 0) {
	int s;

	HC->probmask = malloc(pdinfo->n + 1);
	HC->fullmask = calloc(T, 1);

	if (HC->probmask == NULL ||
	    HC->fullmask == NULL) {
	    return E_ALLOC;
	}

	memset(HC->probmask, '0', pdinfo->n);
	HC->probmask[pdinfo->n] = 0;

	for (s=0, t=pdinfo->t1; t<=pdinfo->t2; t++, s++) {
	    int miss = 0;

	    if (Z[HC->selvar][t] == 1.0) {
		for (i=1; i<=reglist[0] && !miss; i++) {
		    if (na(Z[reglist[i]][t])) {
			HC->fullmask[s] = 1;
			HC->probmask[t] = '1';
			miss = 1;
		    }
		}
	    }

	    for (i=1; i<=sellist[0] && !miss; i++) {
		if (na(Z[sellist[i]][t])) {
		    HC->fullmask[s] = 1;
		    HC->probmask[t] = '1';
		    miss = 1;
		}
	    }
	}
    }

    if (HC->probmask != NULL) {
#if HDEBUG
	fprintf(stderr, "Doing copy_to_reference_missmask\n");
#endif
	err = copy_to_reference_missmask(HC->probmask);
    }

    return err;
}

static int make_uncens_mask (h_container *HC, const double **Z, 
			     DATAINFO *pdinfo)
{
    int T = sample_size(pdinfo);
    int t, s = 0;

    HC->uncmask = calloc(T, 1);
    if (HC->uncmask == NULL) {
	return E_ALLOC;
    }

    if (HC->fullmask == NULL) {
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (Z[HC->selvar][t] == 0.0) {
		HC->uncmask[s] = 1;
	    }
	    s++;
	}
    } else {
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (Z[HC->selvar][t] == 0.0 || HC->fullmask[s]) {
		HC->uncmask[s] = 1;
	    }
	    s++;
	}
    }

    return 0;
}

#if HDEBUG > 1

static void debug_print_masks (h_container *HC,
			       const int *Xl,
			       const double **Z)
{
    double x;
    int t, i;

    for (t=HC->t1; t<=HC->t2; t++) {
	if (HC->fullmask[t] == 1) {
	    fputc('F', stderr);
	}
	if (HC->uncmask[t] == 1) {
	    fputc('U', stderr);
	}
	fputc('\t', stderr);
	x = Z[HC->selvar][t];
	if (na(x)) {
	    fprintf(stderr, "          NA");
	} else {
	    fprintf(stderr, "%12.4f", x);
	}
	for (i=1; i<=Xl[0]; i++) {
	    x = Z[Xl[i]][t];
	    if (na(x)) {
		fprintf(stderr, "          NA");
	    } else {
		fprintf(stderr, "%12.4f", x);
	    }
	}
	fputc('\n', stderr);
    }
}

#endif

/* make copies of the primary and selection lists with the first
   terms (i.e. the dependent variable and the selection variable,
   respectively) removed.
*/

static int h_make_data_lists (h_container *HC,
			      const int *reglist,
			      const int *sellist)
{
    int i;

    HC->Xlist = gretl_list_new(HC->kmain);
    HC->Zlist = gretl_list_new(HC->ksel);

    if (HC->Xlist == NULL || HC->Zlist == NULL) {
	return E_ALLOC;
    }

    for (i=1; i<=HC->Xlist[0]; i++) {
	HC->Xlist[i] = reglist[i+1];
    }

    for (i=1; i<=HC->Zlist[0]; i++) {
	HC->Zlist[i] = sellist[i+1];
    }

    return 0;
}

static int h_container_fill (h_container *HC, const int *Xl, 
			     const int *Zl, const double **Z, 
			     DATAINFO *pdinfo, MODEL *probmod, 
			     MODEL *olsmod)
{
    gretl_vector *tmp = NULL;
    double bmills, s2, mdelta;
    int tmplist[2];
    int t1, t2;
    int i, err = 0;

    /* X does NOT include the Mills ratios: hence the "-1" */
    HC->kmain = olsmod->ncoeff - 1;
    HC->ksel = probmod->ncoeff;

    if (HC->kmain < 1) {
	/* FIXME? */
	gretl_errmsg_set(_("No regression function has been specified"));
	return E_DATA;
    }

    HC->beta = gretl_column_vector_alloc(HC->kmain);
    if (HC->beta == NULL) {
	return E_ALLOC;
    }

    HC->gama = gretl_coeff_vector_from_model(probmod, NULL, &err);
    if (err) {
	gretl_matrix_free(HC->beta);
	HC->beta = NULL;
	return err;
    }

    for (i=0; i<HC->kmain; i++) {
	gretl_vector_set(HC->beta, i, olsmod->coeff[i]);
    }

    HC->lambda = olsmod->coeff[HC->kmain];

    /*
      we'll be working with two masks: fullmask just masks unusable
      observations, and may be NULL; uncmasks skips censored
      observations too.
    */
    err = make_uncens_mask(HC, Z, pdinfo);

#if HDEBUG > 1
    if (!err && HC->fullmask != NULL) {
	debug_print_masks(HC, Xl, Z);
    }
#endif

    t1 = pdinfo->t1;
    t2 = pdinfo->t2;

    tmplist[0] = 1;

    if (!err) {
	/* dependent variable vector */
	tmplist[1] = HC->depvar;
	HC->y = gretl_matrix_data_subset_masked(tmplist, Z, t1, t2, 
						HC->uncmask, &err);
    }

    if (!err) {
	HC->nunc = gretl_matrix_rows(HC->y);
	/* selection variable vector */
	tmplist[1] = HC->selvar;
	HC->d = gretl_matrix_data_subset_masked(tmplist, Z, t1, t2, 
						HC->fullmask, &err);
    }

    if (!err) {
	HC->ntot = gretl_matrix_rows(HC->d);
	/* Mills ratio vector */
	tmplist[1] = HC->millsvar;
	HC->mills = gretl_matrix_data_subset_masked(tmplist, Z, t1, t2, 
						    HC->uncmask, &err);
    }

    if (!err) {
	err = h_make_data_lists(HC, Xl, Zl);
    }

    if (err) {
	return err;
    }

    /* regressors in primary equation */
    HC->reg = gretl_matrix_data_subset_masked(HC->Xlist, Z, t1, t2, 
					      HC->uncmask, &err);

    if (!err) {
	/* regressors in selection equation */
	HC->selreg = 
	    gretl_matrix_data_subset_masked(HC->Zlist, Z, t1, t2, 
					    HC->fullmask, &err);
    }

    if (!err) {
	/* regressors in selection equation, uncensored sample */
	HC->selreg_u = 
	    gretl_matrix_data_subset_masked(HC->Zlist, Z, t1, t2, 
					    HC->uncmask, &err);
    }

    if (err) {
	return err;
    }    

    tmp = gretl_matrix_multiply_new(HC->selreg_u, HC->gama, &err);

    if (!err) {
	gretl_matrix_add_to(tmp, HC->mills);
	HC->delta = gretl_matrix_dot_op(tmp, HC->mills, '*', &err);
	gretl_matrix_free(tmp);
    }

    if (!err) {
	bmills = olsmod->coeff[olsmod->ncoeff-1];
	mdelta = gretl_vector_mean(HC->delta);
	s2 = olsmod->ess / olsmod->nobs + mdelta * bmills * bmills;
	HC->sigma = sqrt(s2);
	HC->rho = bmills / HC->sigma;
    }

    if (!err && fabs(HC->rho) > 1.0) {
	HC->rho = (HC->rho < 0)? -1 : 1;
	/* ensure consistency */
	HC->sigma = HC->lambda / HC->rho;
    }

    if (!err) {
	HC->VProbit = gretl_vcv_matrix_from_model(probmod, NULL, &err);
    }

    if (!err) {
	HC->fitted = gretl_matrix_alloc(HC->nunc, 1);
	HC->u      = gretl_matrix_alloc(HC->nunc, 1);
	HC->ndx    = gretl_matrix_alloc(HC->ntot, 1);

	if (HC->fitted == NULL || HC->u == NULL || HC->ndx == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static double h_loglik (const double *param, void *ptr)
{
    h_container *HC = (h_container *) ptr;
    double lnsig, x, ll = NADBL;
    int kmax = HC->kmain + HC->ksel;
    int i, j, err = 0;
    
    for (i=0; i<HC->kmain; i++) {
	gretl_vector_set(HC->beta, i, param[i]);
    }

    j = 0;
    for (i=HC->kmain; i<kmax; i++) {
	gretl_vector_set(HC->gama, j++, param[i]);
    }

    HC->sigma = param[kmax];
    lnsig = log(HC->sigma);

    HC->rho = param[kmax+1];

#if HDEBUG > 1
    gretl_matrix_print(HC->beta, "beta");
    gretl_matrix_print(HC->gama, "gama");
    fprintf(stderr, "sigma = %12.6f, rho = %12.6f", HC->sigma, HC->rho);
    fputc('\n', stderr);
#endif

    if (HC->sigma <= 0 || fabs(HC->rho) >= 1) {
	return NADBL;
    } 

    err = gretl_matrix_multiply(HC->reg, HC->beta, HC->fitted);
    
    if (!err) {
	gretl_matrix_copy_values(HC->u, HC->y);
	err = gretl_matrix_subtract_from(HC->u, HC->fitted);
    }

    if (!err) {
	err = gretl_matrix_divide_by_scalar(HC->u, HC->sigma);
    }

    if (!err) {
	err = gretl_matrix_multiply(HC->selreg, HC->gama, HC->ndx);
    }

    if (!err) {
	double rhofunc = 1 / sqrt(1 - HC->rho * HC->rho);
	double ll0 = 0, ll1 = 0, ll2 = 0;
	double ut, ndxt;
	int sel;

	/* i goes through all obs, while j keeps track of the uncensored
	   ones */
	j = 0;
	for (i=0; i<HC->ntot; i++) {
	    sel = (1.0 == gretl_vector_get(HC->d, i));
	    ndxt = gretl_vector_get(HC->ndx, i);
	    if (sel) {
		ut = gretl_vector_get(HC->u, j++);
		x = (ndxt + HC->rho*ut) * rhofunc;
		ll1 += log(normal_pdf(ut)) - lnsig;
		ll2 += log(normal_cdf(x));
	    } else {
		ll0 += log(normal_cdf(-ndxt));
	    }
	}
	ll = ll0 + ll1 + ll2;
#if HDEBUG
	fprintf(stderr, "ll0 = %g, ll1 = %g, ll2 = %g: ll = %g (ntot = %d)\n", 
		ll0, ll1, ll2, ll, HC->ntot);
#endif
    }

    return ll;
}

/*
  What we should do here is not entirely clear: we set yhat to the
  linear predictor, computed over uncensored AND censored
  observations.  

  On the contrary, uhat is defined as 
  y - yhat	for uncensored observations, and as
  E(u|d=0)	for censored observations
*/

static void heckit_yhat_uhat (MODEL *hm, h_container *HC, 
			      const double **Z, DATAINFO *pdinfo)
{
    double x, xb, zg, u;
    double f, F, lam = HC->lambda;
    double rho = HC->rho;
    double rhofunc = 1 / sqrt(1 - rho * rho);
    double lnsig = log(HC->sigma);
    int kb = HC->kmain;
    int kg = HC->ksel;
    int i, t, v;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	xb = 0.0;
	for (i=0; i<kb; i++) {
	    v = HC->Xlist[i+1];
	    x = Z[v][t];
	    if (na(x)) {
		xb = NADBL;
		break;
	    } else {
		xb += x * gretl_vector_get(HC->beta, i);
	    }
	}

	if (!na(xb)) {
	    hm->yhat[t] = xb;
	}

	zg = 0.0;
	for (i=0; i<kg; i++) {
	    v = HC->Zlist[i+1];
	    x = Z[v][t];
	    if (na(x)) {
		zg = NADBL;
		break;
	    } else {
		zg += x * gretl_vector_get(HC->gama, i);
	    }
	}

	if (na(zg)) {
	    hm->uhat[t] = NADBL;
	    continue;
	}

	if (Z[HC->selvar][t] == 0) {
	    /* censored */
	    f = normal_pdf(zg);
	    F = normal_cdf(-zg);
	    hm->uhat[t] = -lam * f / F;
	    hm->llt[t] = log(F);
	} else if (!na(xb) && Z[HC->selvar][t] == 1) {
	    /* uncensored */
	    hm->uhat[t] = u = Z[HC->depvar][t] - xb;
	    u /= HC->sigma;
	    x = (zg + rho*u) * rhofunc;
	    hm->llt[t] = -LN_SQRT_2_PI - 0.5*u*u - lnsig + log(normal_cdf(x));
	}
    }

#if 0
    /* R-squared based on correlation? */
    x = gretl_corr(hm->t1, hm->t2, Z[depvar], hm->yhat, NULL);
    if (!na(x)) {
	hm->rsq = x * x;
    }
#endif
}

/*
  This function works the same way for the 2-step and the ML
  estimators: all the relevant items are taken from the container HC
  anyway.
*/

static int transcribe_heckit_params (MODEL *hm, h_container *HC, DATAINFO *pdinfo)
{
    double *fullcoeff;
    int kb = HC->kmain;
    int npar = kb + HC->ksel + 1;
    int i, err = 0;

    fullcoeff = malloc(npar * sizeof *fullcoeff);
    if (fullcoeff == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_model_allocate_params(hm, npar);
    }

    if (!err) {
	for (i=0; i<kb; i++) {
	    strcpy(hm->params[i], pdinfo->varname[hm->list[i+2]]);
	    fullcoeff[i] = gretl_vector_get(HC->beta, i);
	}
	
	strcpy(hm->params[kb], "lambda");
	fullcoeff[kb] = HC->lambda;

	for (i=0; i<HC->ksel; i++) {
	    strcpy(hm->params[i+kb+1], pdinfo->varname[HC->Zlist[i+1]]);
	    fullcoeff[i+kb+1] = gretl_vector_get(HC->gama, i);
	}
    }

    hm->lnL = HC->ll;
    mle_criteria(hm, 0); /* FIXME? */
    hm->sigma = HC->sigma;
    hm->rho = HC->rho;

    if (!err) {
	free(hm->list);
	hm->list = gretl_list_copy(HC->list);
	free(hm->coeff);
	hm->coeff = fullcoeff;
	hm->ncoeff = npar;
	hm->t1 = HC->t1;
	hm->t2 = HC->t2;
	gretl_model_set_coeff_separator(hm, N_("Selection equation"), kb + 1);
	gretl_model_set_int(hm, "base-coeffs", kb);
    }
    
    return err;
}

static int transcribe_2step_vcv (MODEL *pmod, const gretl_matrix *S, 
				  const gretl_matrix *Vp)
{
    int m = S->rows;
    int n = Vp->rows;
    int nvc, tot = m + n;
    int i, j, k;
    double vij;

    if (pmod->vcv != NULL) {
	free(pmod->vcv);
    }

    if (pmod->sderr != NULL) {
	free(pmod->sderr);
    }
    
    nvc = (tot * tot + tot) / 2;

    pmod->vcv = malloc(nvc * sizeof *pmod->vcv);
    pmod->sderr = malloc(tot * sizeof *pmod->sderr);

    if (pmod->vcv == NULL || pmod->sderr == NULL) {
	return E_ALLOC;
    }

    /* build combined covariance matrix */

    k = 0;
    for (i=0; i<tot; i++) {
	for (j=i; j<tot; j++) {
	    if (i < m) {
		/* How to handle the "irrelevant" block?  In the line below
		   these elements were being set to NADBL; for the moment
		   I'm setting them to zero -- AC, 2009-04-11.
		*/
		vij = (j < m)? gretl_matrix_get(S, i, j) : 0.0;
	    } else {
		vij = gretl_matrix_get(Vp, i-m, j-m);
	    }
	    pmod->vcv[k++] = vij;
	    if (i == j) {
		pmod->sderr[i] = sqrt(vij);
	    }
	}
    }

    return 0;
}

/*
  This function computes the adjusted VCV matrix for the 2-step estimator
  and copies it into the model struct
*/

static int heckit_2step_vcv (h_container *HC, MODEL *olsmod)
{
    gretl_matrix *W = HC->selreg_u;
    gretl_matrix *w = HC->delta;
    gretl_matrix *Vp = HC->VProbit;
    gretl_matrix_block *B = NULL;
    gretl_matrix *X, *XX, *XXi;
    gretl_matrix *XXw, *XwZ, *S;
    gretl_matrix *Xw = NULL;
    int nX = HC->kmain + 1;
    int nZ = HC->ksel;
    int err = 0;

    X = gretl_matrix_col_concat(HC->reg, HC->mills, &err);

    if (!err) {
	Xw = gretl_matrix_dot_op(w, X, '*', &err);
    }

    if (!err) {
	B = gretl_matrix_block_new(&XX, nX, nX,
				   &XXi, nX, nX,
				   &XXw, nX, nX,
				   &XwZ, nX, nZ,
				   &S, nX, nX,
				   NULL);
	if (B == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	goto bailout;
    }

    gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE, 
			      X, GRETL_MOD_NONE, 
			      XX, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE, 
			      Xw, GRETL_MOD_NONE, 
			      XXw, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(Xw, GRETL_MOD_TRANSPOSE, 
			      W, GRETL_MOD_NONE, 
			      XwZ, GRETL_MOD_NONE);

    gretl_matrix_copy_values(XXi, XX);
    err = gretl_invert_symmetric_matrix(XXi); 
    if (err) {
	fprintf(stderr, "heckit_2step_vcv: error inverting X'X\n");
	goto bailout;
    }

    gretl_matrix_qform(XwZ, GRETL_MOD_NONE,
		       Vp, S, GRETL_MOD_NONE);
    gretl_matrix_subtract_from(XXw, S);
    gretl_matrix_multiply_by_scalar(XXw, HC->rho * HC->rho);
    gretl_matrix_subtract_from(XX, XXw);
    gretl_matrix_qform(XXi, GRETL_MOD_NONE,
		       XX, S, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(S, HC->sigma * HC->sigma);

    err = transcribe_2step_vcv(olsmod, S, Vp);

 bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(Xw);
    gretl_matrix_block_destroy(B);

    return err;
}

/*
  since lambda is not a parameter that is directly estimated by ML, 
  we compute the augmented vcv matrix via the delta method
*/

int add_lambda_to_ml_vcv (h_container *HC)
{
    gretl_matrix *J = NULL;
    gretl_matrix *tmp = NULL;
    int npar = HC->vcv->rows;
    int kb = HC->kmain;
    int i, err = 0;

    tmp = gretl_matrix_alloc(npar+1, npar+1);
    J = gretl_zero_matrix_new(npar+1, npar);

    if (tmp == NULL || J == NULL) {
	gretl_matrix_free(tmp);
	gretl_matrix_free(J);
	return E_ALLOC;
    }

    for (i=0; i<kb; i++) {
	gretl_matrix_set(J, i, i, 1);
    }

    gretl_matrix_set(J, kb, npar-2, HC->rho);
    gretl_matrix_set(J, kb, npar-1, HC->sigma);

    for (i=kb+1; i<=npar; i++) {
	gretl_matrix_set(J, i, i-1, 1);
    }

    gretl_matrix_qform(J, GRETL_MOD_NONE, HC->vcv, tmp, GRETL_MOD_NONE);

    gretl_matrix_free(J);
    gretl_matrix_free(HC->vcv);
    HC->vcv = tmp;
    
    return err;
}

int heckit_ml (MODEL *hm, h_container *HC, PRN *prn)
{
    int maxit, fncount, grcount;
    double hij, rho, toler;
    double *hess = NULL;
    double *theta = NULL;
    int i, j, k, np = HC->kmain + HC->ksel + 2;
    int err = 0;

    theta = malloc(np * sizeof *theta);
    if (theta == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<HC->kmain; i++) {
	theta[i] = gretl_vector_get(HC->beta, i);
    }

    j = 0;
    for (i=HC->kmain; i<np-2; i++) {
	theta[i] = gretl_vector_get(HC->gama, j++);
    }

    theta[np-2] = HC->sigma;
    rho = HC->rho;

    if (fabs(rho) > 0.99) {
	rho = (rho > 0)? 0.99 : -0.99;
    }

    theta[np-1] = rho;

    BFGS_defaults(&maxit, &toler, HECKIT);

    err = BFGS_max(theta, np, maxit, toler, 
		   &fncount, &grcount, h_loglik, C_LOGLIK,
		   NULL, HC, (prn != NULL)? OPT_V : OPT_NONE, prn);

    if (!err) {
	HC->ll = hm->lnL = h_loglik(theta, HC);
	gretl_model_set_int(hm, "fncount", fncount);	
	gretl_model_set_int(hm, "grcount", grcount);	
	HC->lambda = HC->sigma * HC->rho;
	hess = numerical_hessian(theta, np, h_loglik, HC, &err);
    }

    if (!err) {
	HC->vcv = gretl_matrix_alloc(np, np);
	if (HC->vcv == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	k = 0;
	for (i=0; i<np; i++) {
	    for (j=i; j<np; j++) {
		hij = hess[k++];
		gretl_matrix_set(HC->vcv, i, j, hij);
		if (i != j) {
		    gretl_matrix_set(HC->vcv, j, i, hij);
		}
	    }
	}
	add_lambda_to_ml_vcv(HC);

#if HDEBUG
	for (i=0; i<np; i++) {
	    hij = gretl_matrix_get(HC->vcv, i, i);
	    fprintf(stderr, "theta[%d] = %12.6f, (%12.6f)\n", i, theta[i], sqrt(hij));
	}
#endif
    }

    free(hess);
    free(theta);

    return err;
}

/*
  This function just copies the VCV matrix for the ML estimator into
  the model struct.  Note that we don't transcribe the variances for
  sigma and rho into the model.
*/

static int transcribe_ml_vcv (MODEL *pmod, h_container *HC)
{
    int nvc, npar = HC->vcv->rows - 2;
    int i, j, k = 0;
    double vij;

    if (pmod->vcv != NULL) {
	free(pmod->vcv);
    }

    if (pmod->sderr != NULL) {
	free(pmod->sderr);
    }
    
    nvc = (npar * npar + npar) / 2;

    pmod->vcv = malloc(nvc * sizeof *pmod->vcv);
    pmod->sderr = malloc(npar * sizeof *pmod->sderr);

    if (pmod->vcv == NULL || pmod->sderr == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<npar; i++) {
	for (j=i; j<npar; j++) {
	    vij = gretl_matrix_get(HC->vcv, i, j);
	    pmod->vcv[k++] = vij;
	    if (i == j) {
		pmod->sderr[i] = sqrt(vij);
	    }
	}
    }

    return 0;
}

static MODEL heckit_init (h_container *HC, double ***pZ, DATAINFO *pdinfo)
{
#if HDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
#endif
    MODEL hm;
    MODEL probmod;
    int *reglist = NULL;
    int *sellist = NULL;
    int t, err = 0;

    gretl_model_init(&hm);
    gretl_model_init(&probmod);

    err = gretl_list_split_on_separator(HC->list, &reglist, &sellist);
    if (!err && (reglist == NULL || sellist == NULL)) {
	err = E_ARGS;
    }

    if (!err) {
	HC->depvar = reglist[1];
	HC->selvar = sellist[1];
	err = make_heckit_NA_mask(HC, reglist, sellist, (const double **) *pZ, 
				  pdinfo);
    }

    if (err) {
	goto bailout;
    }

    /* run initial auxiliary probit */
    probmod = logit_probit(sellist, pZ, pdinfo, PROBIT, OPT_A, NULL);
    if (probmod.errcode) {
	hm.errcode = probmod.errcode;
	goto bailout;
    }

    HC->t1 = probmod.t1;
    HC->t2 = probmod.t2;

#if HDEBUG
    printmodel(&probmod, pdinfo, OPT_NONE, prn);
#endif

    /* add the inverse Mills ratio to the dataset */

    err = dataset_add_series(1, pZ, pdinfo);
    if (err) {
	goto bailout;
    } 

    HC->millsvar = pdinfo->v - 1;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if ((*pZ)[HC->selvar][t] == 1.0) {
	    (*pZ)[HC->millsvar][t] = probmod.uhat[t];
	} else {
	    (*pZ)[HC->millsvar][t] = NADBL;
	}
    }

    strcpy(pdinfo->varname[HC->millsvar], "lambda");
    gretl_list_append_term(&reglist, HC->millsvar);

    /* run OLS including Mills variable */
    hm = lsq(reglist, pZ, pdinfo, OLS, OPT_A);
    if (!hm.errcode) {
	hm.ci = HECKIT;
	gretl_model_set_int(&hm, "totobs", probmod.nobs); 
    }

#if HDEBUG
    printmodel(&hm, pdinfo, OPT_NONE, prn);
    gretl_print_destroy(prn);
#endif

    if (!hm.errcode) {
	/* fill the container with all the relevant data */
	err = h_container_fill(HC, reglist, sellist, (const double **) *pZ, 
			       pdinfo, &probmod, &hm);
    }

 bailout:

    clear_model(&probmod);
    free(reglist);
    free(sellist);

    if (err && !hm.errcode) {
	hm.errcode = err;
    }

    return hm;
}

/* 
   The driver function for the plugin. 

   The function heckit_init runs the 2-step estimation and fills up
   the container HC with all the various items to be used later. Then,
   if we just want a 2-step estimate we adjust the vcv matrix and
   transcribe the result into the returned model. Otherwise, we go for
   ML.
*/

MODEL heckit_estimate (const int *list, double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn) 
{
    h_container *HC = NULL;
    MODEL hm;
    int oldv = pdinfo->v;
    int err = 0;

    gretl_model_init(&hm);

    HC = h_container_new(list);
    if (HC == NULL) {
	hm.errcode = E_ALLOC;
	return hm;
    }

    hm = heckit_init(HC, pZ, pdinfo);
    if (hm.errcode) {
	h_container_destroy(HC);
	return hm;
    }

    if (opt & OPT_T) {
	/* two-step: compute appropriate correction to covariances */
	err = heckit_2step_vcv(HC, &hm);
    } else {
	/* use MLE */
	err = heckit_ml(&hm, HC, prn);
	if (!err) {
	    err = transcribe_ml_vcv(&hm, HC);
	}
    } 

    if (err) {
	hm.errcode = err;
    } else {
	err = transcribe_heckit_params(&hm, HC, pdinfo);
	heckit_yhat_uhat(&hm, HC, (const double **) *pZ, pdinfo);
    }

    if (err && hm.errcode == 0) {
	hm.errcode = err;
    }

    if (hm.errcode == 0 && (opt & OPT_T)) {
	hm.opt |= OPT_T;
    }

    if (pdinfo->v > oldv) {
	dataset_drop_last_variables(1, pZ, pdinfo);
    }

    h_container_destroy(HC);

#if 0
    int t;

    for (t=hm.t1; t<=hm.t2; t++) {
	fprintf(stderr, "%4d: u = %12.6f, llt = %12.6f, missing = %d\n", 
		t, hm.uhat[t], hm.llt[t], model_missing(&hm, t));
    }
#endif

    return hm;
}
