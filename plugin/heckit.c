/* 
 * Copyright (C) Riccardo "Jack" Lucchetti and Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "libgretl.h"
#include "matrix_extra.h"
#include "libset.h"
#include "missing_private.h"

#define HDEBUG 0
#define HTOL 1.0E-09

typedef struct h_container_ h_container;

struct h_container_ {
    int kmain;			/* no. of params in the main eq. */
    int ksel;			/* no. of params in the selection eq. */
    double ll;			/* log-likelihood */

    int ntot, nunc;		/* total & uncensored obs */
    int *Xlist;			/* regressor list for the main eq. */
    int *Zlist;			/* regressor list for the selection eq. */

    gretl_matrix *y;		/* dependent var */
    gretl_matrix *reg;		/* main eq. regressors */
    gretl_matrix *mills;	/* Mills ratios from selection eq */
    gretl_matrix *delta;	/* used in 2-step vcv calculations */
    gretl_matrix *d;		/* selection dummy variable */
    gretl_matrix *selreg;	/* selection eq. regressors */ 
    gretl_matrix *selreg_u;	/* selection eq. regressors (subsample d==1) */
    gretl_vector *fitted;
    gretl_vector *u;
    gretl_vector *ndx;

    gretl_vector *beta;		/* main eq. parameters */
    gretl_vector *gama;		/* selection eq. parameters */
    double sigma;
    double rho;
    double lambda;		/* rho*sigma by definiton */

    gretl_matrix *vcv;		/* Variance-covariance matrix */
    gretl_matrix *VProbit;	/* 1st stage probit covariance matrix */

    char *fullmask;		/* mask NAs */
    char *uncmask;		/* mask NAs and (d==0) */

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

    free(HC->fullmask);
    free(HC->uncmask);

    free(HC);
}

static h_container *h_container_new (void)
{
    h_container *HC = malloc(sizeof *HC);

    if (HC == NULL) {
	return NULL;
    }

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

    HC->fullmask = NULL;
    HC->uncmask = NULL;

    return HC;
}

static char *heckit_suppress_mask (int sel, const int *list, const double **Z,
				   const DATAINFO *pdinfo, int *err)
{
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    char *mask = NULL;
    int i, t, s = 0;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (Z[sel][t]) {
	    for (i=1; i<=list[0]; i++) {
		if (na(Z[list[i]][t])) {
		    s++;
		    break;
		}
	    }
	}
    }

    if (s == T) {
	*err = E_DATA;
	return NULL;
    } else if (s > 0) {
	mask = calloc(T, 1);
	if (mask == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	}

	s = 0;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (Z[sel][t]) {
		for (i=1; i<=list[0]; i++) {
		    if (na(Z[list[i]][t])) {
			mask[s] = 1;
			break;
		    }
		}
	    }
	    s++;
	}
    }

    return mask;
}

static int make_uncens_mask (h_container *HC, const int *Xl, 
			     const int *Zl, const double **Z, 
			     DATAINFO *pdinfo)
{
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int t, s = 0;

    HC->uncmask = calloc(T, 1);
    if (HC->uncmask == NULL) {
	return E_ALLOC;
    }

    if (HC->fullmask == NULL) {
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (Z[Zl[1]][t] == 0.0) {
		HC->uncmask[s] = 1;
	    }
	    s++;
	}
    } else {
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (Z[Zl[1]][t] == 0.0 || HC->fullmask[s]) {
		HC->uncmask[s] = 1;
	    }
	    s++;
	}
    }

    return 0;
}


static int h_container_fill (h_container *HC, const int *Xl, 
			     const int *Zl, const double **Z, 
			     DATAINFO *pdinfo, MODEL *probmod, 
			     MODEL *olsmod)
{
    gretl_vector *gama = NULL;
    gretl_vector *tmp = NULL;
    double bmills, s2, mdelta;
    int tmplist[2];
    int t1 = pdinfo->t1;
    int t2 = pdinfo->t2;
    int v, i, err = 0;

    /* it is assumed that the Mills ratios have just been added to the dataset;
       hence, they're the last variable
    */
    v = pdinfo->v - 1;

    /* X does NOT include the Mills ratios: hence the "-2" */
    HC->kmain = Xl[0] - 2;
    HC->ksel = Zl[0] - 1;

    HC->beta = gretl_vector_alloc(HC->kmain);
    HC->gama = gretl_vector_alloc(HC->ksel);

    for (i=0; i<HC->kmain; i++) {
	gretl_vector_set(HC->beta, i, olsmod->coeff[i]);
    }

    HC->lambda = olsmod->coeff[HC->kmain];

    for (i=0; i<HC->ksel; i++) {
	gretl_vector_set(HC->gama, i, probmod->coeff[i]);
    }

    /*
      we'll be working with 2 masks: fullmask just masks unusable
      observations, and may be NULL; uncmasks skips censored
      observations too.
    */
    err = make_uncens_mask(HC, Xl, Zl, Z, pdinfo);

#if HDEBUG
    double x;
    int t;

    if (HC->fullmask != NULL) {
	for (t=t1; t<=t2; t++) {
	    if (HC->fullmask[t] == '1') {
		fputc('F', stderr);
	    }
	    if (HC->uncmask[t] == '1') {
		fputs('U', stderr);
	    }
	    fputc('\t', stderr);
	    x = Z[Zl[1]][t];
	    fprintf(stderr, "%12.4f", x);
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

    tmplist[0] = 1;
    tmplist[1] = Xl[1];

    HC->y = gretl_matrix_data_subset(tmplist, Z, t1, t2, HC->uncmask);
    HC->nunc = gretl_matrix_rows(HC->y);

    tmplist[1] = Zl[1];
    HC->d = gretl_matrix_data_subset(tmplist, Z, t1, t2, HC->fullmask);
    HC->ntot = gretl_matrix_rows(HC->d);

    tmplist[1] = v;
    HC->mills = gretl_matrix_data_subset(tmplist, Z, t1, t2, HC->uncmask);

    HC->Xlist = gretl_list_new(HC->kmain);
    HC->Zlist = gretl_list_new(HC->ksel);

    if (HC->Xlist == NULL || HC->Zlist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (i=1; i<=HC->Xlist[0]; i++) {
	HC->Xlist[i] = Xl[i+1];
    }

    for (i=1; i<=HC->Zlist[0]; i++) {
	HC->Zlist[i] = Zl[i+1];
    }

    HC->reg = gretl_matrix_data_subset(HC->Xlist, Z, t1, t2, HC->uncmask);
    HC->selreg = gretl_matrix_data_subset(HC->Zlist, Z, t1, t2, HC->fullmask);
    HC->selreg_u = gretl_matrix_data_subset(HC->Zlist, Z, t1, t2, HC->uncmask);

    gama = gretl_coeff_vector_from_model(probmod, NULL);

    HC->delta = gretl_column_vector_alloc(HC->nunc);
    err = gretl_matrix_multiply(HC->selreg_u, gama, HC->delta);
    gretl_matrix_add_to(HC->delta, HC->mills);
    
    tmp = gretl_matrix_dot_op(HC->delta, HC->mills, '*', &err);
    gretl_matrix_free(HC->delta);
    HC->delta = tmp;

    bmills = olsmod->coeff[olsmod->ncoeff-1];
    mdelta  = gretl_vector_mean(HC->delta);
    s2 = olsmod->ess / olsmod->nobs + mdelta * bmills * bmills;
    HC->sigma = sqrt(s2);
    HC->rho = bmills / HC->sigma;

    HC->fitted = gretl_matrix_alloc(HC->nunc,1);
    HC->u = gretl_matrix_alloc(HC->nunc,1);
    HC->ndx = gretl_matrix_alloc(HC->ntot,1);

    HC->VProbit = gretl_vcv_matrix_from_model(probmod, NULL);

 bailout:

    gretl_vector_free(gama);

    return err;
}

static double h_loglik(const double *param, void *ptr)
{
    h_container *HC = (h_container *) ptr;
    double lnsig, isqrtrhoc, x, ll = NADBL;
    int k1 = HC->kmain;
    int k2 = HC->ksel;
    int i, j, err = 0;
    
    for (i=0; i<k1; i++) {
	gretl_vector_set(HC->beta, i, param[i]);
    }

    j = 0;
    for (i=k1; i<k1+k2; i++) {
	gretl_vector_set(HC->gama, j++, param[i]);
    }

    HC->sigma = param[k1+k2];
    lnsig = log(HC->sigma);

    HC->rho = param[k1+k2+1];

#if HDEBUG
    gretl_matrix_print(HC->beta, "beta");
    gretl_matrix_print(HC->gama, "gama");
    fprintf(stderr, "sigma = %12.6f, rho = %12.6f", HC->sigma, HC->rho);
    fputc('\n', stderr);
#endif

    if (HC->sigma <= 0 || fabs(HC->rho) >= 1) {
	return ll;
    } else {
	isqrtrhoc = 1 / sqrt(1 - HC->rho * HC->rho);
    }

    err = gretl_matrix_multiply_mod(HC->reg, GRETL_MOD_NONE,
				    HC->beta, GRETL_MOD_TRANSPOSE,
				    HC->fitted, GRETL_MOD_NONE);
    
    if(!err) {
	gretl_matrix_copy_values(HC->u, HC->y);
	err = gretl_matrix_subtract_from(HC->u, HC->fitted);
    }

    if(!err) {
	err = gretl_matrix_divide_by_scalar(HC->u, HC->sigma);
    }

    if(!err) {
	err = gretl_matrix_multiply_mod(HC->selreg, GRETL_MOD_NONE,
					HC->gama, GRETL_MOD_TRANSPOSE,
					HC->ndx, GRETL_MOD_NONE);
    }

    if (!err) {
	double ut, ndxt;
	double ll0 = 0;
	double ll1 = 0;
	double ll2 = 0;
	int sel;
	/* 
	   i goes through all obs, while j keeps track of the uncensored
	   ones
	*/
	j = 0;
	for (i=0; i<HC->ntot; i++) {
	    sel = (1.0 == gretl_vector_get(HC->d, i));
	    ndxt = gretl_vector_get(HC->ndx,i);
	    if (sel) {
		ut = gretl_vector_get(HC->u, j++);
		x = (ndxt + HC->rho*ut) * isqrtrhoc;
		ll1 += log(normal_pdf(ut)) - lnsig;
		ll2 += log(normal_cdf(x));
	    } else {
		ll0 += log(normal_cdf(-ndxt));
	    }
	}
	
#if HDEBUG
	fprintf(stderr, "ll0 = %g, ll1 = %g, ll2 = %g\n", ll0, ll1, ll2);
#endif
	ll = ll0 + ll1 + ll2;
    }

    return ll;
}

#if 0
/* remove the "lambda" effect from the fitted values */
/* TODO: this may be wrong and is obsolete anyway */

static void fix_heckit_resids (MODEL *hm, h_container *HC)
{
    int t, k = hm->ncoeff - 1;
    double x;

    for (t=hm->t1; t<=hm->t2; t++) {
	if (!na(hm->uhat[t])) {
	    x = hm->coeff[k] * lam[t];
	    hm->uhat[t] += x;
	    hm->yhat[t] -= x;
	}
    }
}
#endif

/*
  This function works the same way for the 2-step and the ML
  estimators: all the relevant items are taken from the container HC
  anyway.
*/

static int transcribe_heckit_params (MODEL *hm, h_container *HC, DATAINFO *pdinfo)
{
    double *fullcoeff;
    int ko = HC->kmain;
    int kp = HC->ksel;
    int k = ko + kp + 1;
    int i, err = 0;

    fullcoeff = malloc(k * sizeof *fullcoeff);
    if (fullcoeff == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_model_allocate_params(hm, k);
    }

    if (!err) {
	for (i=0; i<ko; i++) {
	    strcpy(hm->params[i], pdinfo->varname[hm->list[i+2]]);
	    fullcoeff[i] = gretl_vector_get(HC->beta, i);
    
	}
	
	strcpy(hm->params[ko], "lambda");
	fullcoeff[ko] = HC->lambda;

	for (i=0; i<kp; i++) {
	    strcpy(hm->params[i+ko+1], pdinfo->varname[HC->Zlist[i+1]]);
	    fullcoeff[i+ko+1] = gretl_vector_get(HC->gama, i);
	}

    }

    hm->lnL = HC->ll;
    hm->sigma = HC->sigma;
    hm->rho = HC->rho;

    if (!err) {
	free(hm->coeff);
	hm->coeff = fullcoeff;
	hm->ncoeff = k;
	gretl_model_set_coeff_separator(hm, N_("Selection equation"), ko + 1);
	gretl_model_set_int(hm, "base-coeffs", ko);
	hm->list[0] -= 1;
    }
    
    return err;
}

static int transcribe_2step_vcv (MODEL *pmod, const gretl_matrix *S, 
				  const gretl_matrix *Vp)
{
    int m = S->rows;
    int n = Vp->rows;
    int nvc, tot = m + n;
    int i, j, k = 0;
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

    for (i=0; i<tot; i++) {
	for (j=i; j<tot; j++) {
	    if (i < m) {
		vij = (j < m)? gretl_matrix_get(S, i, j) : NADBL;
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
    gretl_matrix *X = NULL;
    gretl_matrix *Xw = NULL;
    gretl_matrix *W = HC->selreg_u;
    gretl_matrix *w = HC->delta;

    gretl_matrix *XX = NULL;
    gretl_matrix *XXi = NULL;
    gretl_matrix *XXw = NULL;
    gretl_matrix *XwZ = NULL;
    gretl_matrix *S = NULL;

    gretl_matrix *Vp = HC->VProbit;

    int nX = HC->kmain + 1;
    int nZ = HC->ksel;
    int err = 0;

    double s2, sigma, rho;

    sigma = HC->sigma;
    rho = HC->rho;
    s2 = sigma * sigma;

    X = gretl_matrix_col_concat(HC->reg, HC->mills, &err);
    if (err) {
	goto bailout;
    }

    Xw = gretl_matrix_dot_op(w, X, '*', &err);
    if (err) {
	goto bailout;
    }

    XX = gretl_matrix_alloc(nX, nX);
    XXi = gretl_matrix_alloc(nX, nX);
    XXw = gretl_matrix_alloc(nX, nX);
    XwZ = gretl_matrix_alloc(nX, nZ);
    S = gretl_matrix_alloc(nX, nX);

    if (X == NULL || w == NULL || W == NULL || Xw == NULL ||
	XX == NULL || XXi == NULL || XXw == NULL ||
	XwZ == NULL || S == NULL) {
	err = E_ALLOC;
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
    gretl_matrix_multiply_by_scalar(XXw, rho * rho);
    gretl_matrix_subtract_from(XX, XXw);
    gretl_matrix_qform(XXi, GRETL_MOD_NONE,
		       XX, S, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(S, s2);

    err = transcribe_2step_vcv(olsmod, S, Vp);

 bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(Xw);
    gretl_matrix_free(XX);
    gretl_matrix_free(XXi);
    gretl_matrix_free(XXw);
    gretl_matrix_free(XwZ);
    gretl_matrix_free(S);

    return err;
}

/*
  since lambda is not a parameter that is directly estimated by ML, 
  we compute the augmented vcv matrix via the delta method
*/

int add_lambda_to_ml_vcv (h_container *HC)
{
    int i, k, npar = HC->vcv->rows;
    gretl_matrix *J = NULL;
    gretl_matrix *tmp = NULL;
    int err = 0;

    tmp = gretl_matrix_alloc(npar+1, npar+1);
    J = gretl_zero_matrix_new(npar+1, npar);

    if (tmp == NULL || J == NULL) {
	gretl_matrix_free(tmp);
	gretl_matrix_free(J);
	return E_ALLOC;
    }

    k = HC->kmain;

    for (i=0; i<k; i++) {
	gretl_matrix_set(J, i, i, 1);
    }

    gretl_matrix_set(J, k, npar-2, HC->rho);
    gretl_matrix_set(J, k, npar-1, HC->sigma);

    for (i=k+1; i<=npar; i++) {
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
    int fncount, grcount;
    double hij, rho;
    double *hess = NULL;
    double *theta = NULL;
    int i, j, k, np = HC->kmain+HC->ksel+2;
    int err = 0;

    theta = malloc(np * sizeof *theta);
    if (theta == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<HC->kmain; i++) {
	theta[i] = gretl_vector_get(HC->beta,i);
    }

    j = 0;
    for (i=HC->kmain; i<np-2; i++) {
	theta[i] = gretl_vector_get(HC->gama,j++);
    }

    theta[np-2] = HC->sigma;
    rho = HC->rho;

    if (fabs(rho) > 0.99) {
	rho = (rho > 0)? 0.99 : -0.99;
    }

    theta[np-1] = rho;

    err = BFGS_max(theta, np, 1000, HTOL, 
		   &fncount, &grcount, h_loglik, C_LOGLIK,
		   NULL, HC, (prn != NULL)? OPT_V : OPT_NONE, prn);

    if (!err) {
	HC->ll = hm->lnL = h_loglik(theta, HC);
	gretl_model_set_int(hm, "fncount", fncount);	
	gretl_model_set_int(hm, "grcount", grcount);	
	HC->lambda = HC->sigma*HC->rho;
	hess = numerical_hessian(theta, np, h_loglik, HC, &err);
    }

    if (!err) {
        HC->vcv = gretl_matrix_alloc(np,np);
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
  the model struct
*/
static int transcribe_ml_vcv (MODEL *pmod, h_container *HC)
{
    int nvc, npar;

    /*
      We don't transcribe the variances for sigma and rho into the
      model
    */
    npar = HC->vcv->rows - 2;

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

static int translate_to_reference_missmask (const char *shortmask,
					    const DATAINFO *pdinfo)
{
    char *longmask = malloc(pdinfo->n + 1);
    int t, s, err = 0;

    if (longmask == NULL) {
	return E_ALLOC;
    }

    memset(longmask, '0', pdinfo->n);
    longmask[pdinfo->n] = 0;

    s = 0;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (shortmask[s++]) {
	    longmask[t] = '1';
	}
    }

    err = copy_to_reference_missmask(longmask);
    free(longmask);

    return err;
}

static MODEL heckit_init (int *list, double ***pZ, DATAINFO *pdinfo,
			  h_container *HC) 
{

    MODEL hm;
    MODEL probmod;
    int *Xlist = NULL;
    int *Zlist = NULL;
    int v = pdinfo->v;
    int t, sel;
    int err = 0;

    gretl_model_init(&hm);
    gretl_model_init(&probmod);

    err = gretl_list_split_on_separator(list, &Xlist, &Zlist);
    if (err) {
	hm.errcode = err;
	goto bailout;
    } 

    HC->fullmask = heckit_suppress_mask(Zlist[1], Xlist, (const double **) *pZ, 
				pdinfo, &err);

    if (HC->fullmask != NULL) {
	translate_to_reference_missmask(HC->fullmask, pdinfo);
    }

    /* run initial auxiliary probit */
    probmod = logit_probit(Zlist, pZ, pdinfo, PROBIT, OPT_A, NULL);
    if (probmod.errcode) {
	free(Xlist);
	free(Zlist);
	goto bailout;
    }

    /* add the inverse Mills ratio to the dataset */

    err = dataset_add_series(1, pZ, pdinfo);
    if (err) {
	hm.errcode = err;
	clear_model(&probmod);
	goto bailout;
    } 

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	sel = ((*pZ)[Zlist[1]][t] == 1.0);
	(*pZ)[v][t] = (sel) ? probmod.uhat[t] : NADBL;
    }

    strcpy(pdinfo->varname[v], "lambda");
    gretl_list_append_term(&Xlist, v);

    /* run OLS */
    hm = lsq(Xlist, pZ, pdinfo, OLS, OPT_A);
    hm.ci = HECKIT;

    /* FIXME: this definitely doesn't belong here */
    gretl_model_set_int(&hm, "totobs", probmod.nobs);

    /* fill the container with all the relevant data */
    err = h_container_fill(HC, Xlist, Zlist, (const double **) *pZ, 
			   pdinfo, &probmod, &hm);

    clear_model(&probmod);

 bailout:

    free(Xlist);
    free(Zlist);

    return hm;
}

/* 
   The driver function for the plugin. 

   The function heckit_init runs the 2-step estimation and fills up
   the container HC with all the various items to be used later. Then,
   if we just want a 2-step estimate, then we just adjust the vcv
   matrix and transcribe the result into the returned
   model. Otherwise, we go for ML.
*/

MODEL heckit_estimate (int *list, double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn) 
{
    h_container *HC = NULL;
    MODEL hm;
    int err = 0;

    gretl_model_init(&hm);

    HC = h_container_new();
    if (HC == NULL) {
	hm.errcode = E_ALLOC;
	return hm;
    }

    hm = heckit_init(list, pZ, pdinfo, HC);
    if (hm.errcode) {
	h_container_destroy(HC);
	return hm;
    }

    if (opt & OPT_T) {
	/* two-step: compute appropriate correction to covariances */
	err = heckit_2step_vcv(HC, &hm);
    } else {
	err = heckit_ml(&hm, HC, prn);
	if(!err) {
	    err = transcribe_ml_vcv(&hm, HC);
	}
    } 

    if (err) {
	hm.errcode = err;
    } else {
	err = transcribe_heckit_params(&hm, HC, pdinfo);
    }

    if (err) {
	hm.errcode = err;
    }

    dataset_drop_last_variables(1, pZ, pdinfo);
    h_container_destroy(HC);

    return hm;
}
