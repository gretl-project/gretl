/* 
 * Copyright (C) 2004 Riccardo Lucchetti and Allin Cottrell
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

/* The algorithm here was contributed by Riccardo "Jack" Lucchetti; 
   C coding by Allin Cottrell. 
*/

#ifdef STANDALONE
# include <gretl/libgretl.h>
#else
# include "libgretl.h"
# include "internal.h"
#endif

/* #define DEBUG */

#define DEFAULT_MAX_ITER 1000

typedef struct _tobit_info tobit_info;

struct _tobit_info {
    int k;
    int n;
    double ll;
    double ll2;
    int *list;
    double *theta;
    double *delta;
    double *deltmp;
    double *ystar;
    double *P;
    double *e;
    double *f;
};

static int get_maxiter (void)
{
    char *mistr = getenv("GRETL_MAX_ITER");
    int mi = DEFAULT_MAX_ITER;

    if (mistr != NULL) {
	if (!sscanf(mistr, "%d", &mi)) {
	    mi = DEFAULT_MAX_ITER;
	}
    }

    return mi;
}

/* Below: we are buying ourselves a considerable simplification when it comes
   to the tobit_ll function.  That function needs access to the orginal y
   and X data.  But the sample used for estimation may be at an offset into
   the full dataset, and the variables chosen for the analysis may not
   be at contiguous locations in the main dataset.  So we are constructing
   a "virtual dataset" in the form of a set of const pointers into the real
   dataset.  These pointers start at the correct sample offset, and are
   contiguous, so the indexing is a lot easier.
*/

static const double **make_tobit_X (const MODEL *pmod, const double **Z)
{
    const double **X;
    int nv = pmod->list[0];
    int offset = pmod->t1;
    int v, i;

    X = malloc(nv * sizeof *X);
    if (X == NULL) return NULL;

    X[0] = Z[0] + offset;
    X[1] = Z[pmod->list[1]] + offset;

    for (i=2; i<nv; i++) {
	v = pmod->list[i + 1];
#ifdef DEBUG
	fprintf(stderr, "setting X[%d] -> Z[%d]\n", i, v);
#endif
	X[i] = Z[v] + offset;
    }

    return X;
}

/* Below: construct the regression list for the OPG regression, with
   the appropriate indices into the temporary artificial dataset.
*/

static int *make_tobit_list (int k)
{
    int *list;
    int i;

    list = malloc((k + 2) * sizeof *list);
    if (list == NULL) return NULL;

    list[0] = k + 1;
    list[1] = 0;  /* dep var is the constant */
    for (i=0; i<k; i++) {
	list[i+2] = i + 1; 
    }

#ifdef DEBUG
    printlist(list, "OPG regression list");
#endif

    return list;
}

static void tobit_free (tobit_info *tobit)
{
    free(tobit->P);
    free(tobit->e);
    free(tobit->f);
    free(tobit->ystar);
    free(tobit->theta);
    free(tobit->delta);
    free(tobit->deltmp);
    free(tobit->list);
}

/* Collect together much of the required memory allocation, under
   the aegis of the "tobit" struct.
*/

static int tobit_init (tobit_info *tobit, const MODEL *pmod)
{
    int i, k, err = 0;

    tobit->n = pmod->nobs;
    k = tobit->k = pmod->ncoeff;  

    tobit->theta = tobit->delta = tobit->deltmp = NULL;
    tobit->P = tobit->e = tobit->f = tobit->ystar = NULL;
    tobit->list = NULL;

    tobit->list = make_tobit_list(k + 1);
    tobit->theta = malloc((tobit->k + 1) * sizeof *tobit->theta);
    tobit->delta = malloc((tobit->k + 1) * sizeof *tobit->delta);
    tobit->deltmp = malloc((tobit->k + 1) * sizeof *tobit->deltmp);
    tobit->P = malloc(tobit->n * sizeof *tobit->P);
    tobit->e = malloc(tobit->n * sizeof *tobit->e);
    tobit->f = malloc(tobit->n * sizeof *tobit->f);
    tobit->ystar = malloc(tobit->n * sizeof *tobit->ystar);

    if (tobit->P == NULL || tobit->e == NULL || tobit->f == NULL ||
	tobit->ystar == NULL || tobit->theta == NULL || 
	tobit->delta == NULL || tobit->deltmp == NULL ||
	tobit->list == NULL) {
	err = 1;
	tobit_free(tobit);
    }

    /* initialize beta from the OLS coefficients */
    for (i=0; i<k; i++) {
	tobit->theta[i] = pmod->coeff[i];
#ifdef DEBUG
	fprintf(stderr, "initial OLS beta[%d] = %#.5g\n", i, pmod->coeff[i]);
#endif
    }

    /* initialize variance */
    tobit->theta[k] = 1.0;

    return err;
}

/* Compute log likelihood, and the score matrix if do_score is non-zero */

static int tobit_ll (double *theta, const double **X, double **Z, tobit_info *tobit, 
		     int do_score)
{
    const double *y = X[1];
    double siginv = theta[tobit->k]; /* inverse of variance */
    double x, llt;
    int i, t;

    if (siginv < 0.0) {
	fprintf(stderr, "tobit_ll: got a negative variance\n");
	if (!do_score) tobit->ll2 = -1.0e10;
	return 1;
    }

    /* calculate ystar, e, f, and P vectors */
    for (t=0; t<tobit->n; t++) {
	tobit->ystar[t] = theta[0];
	for (i=1; i<tobit->k; i++) {
	    tobit->ystar[t] += theta[i] * X[i+1][t]; /* coeff * orig data */
	}
	tobit->e[t] = y[t] * siginv - tobit->ystar[t];
	tobit->f[t] = siginv * normal_pdf(tobit->e[t]);
	tobit->P[t] = normal_cdf(tobit->ystar[t]);
    }

    /* compute loglikelihood for each obs, cumulate into ll */
    x = 0.0;
    for (t=0; t<tobit->n; t++) {
	if (y[t] == 0.0) {
	    llt = 1.0 - tobit->P[t];
	} else {
	    llt = tobit->f[t];
	}
	x += log(llt);
    }

    /* assign ll to the appropriate record */
    if (do_score) {
	tobit->ll = x;
    } else {
	tobit->ll2 = x;
    }

    if (do_score) {
	int i, gi, xi;
	double den, tail;

	for (t=0; t<tobit->n; t++) {
	    den = normal_pdf(tobit->ystar[t]);
	    tail = 1.0 - tobit->P[t];
	    
	    for (i=0; i<=tobit->k; i++) {

		/* set the indices into the data arrays */
		gi = i + 1;
		xi = (i == 0)? 0 : i + 1;

		if (y[t] == 0.0) {
		    /* score if y is censored */
		    if (xi <= tobit->k) {
			Z[gi][t] = -den / tail * X[xi][t];
		    } else {
			Z[gi][t] = 0.0;
		    } 
#ifdef DEBUG
		    printf("%#12.5g  ", Z[gi][t]);
#endif
		} else {
		    /* score if y is not censored */
		    if (xi <= tobit->k) {
			x = X[xi][t];
		    } else {
			x = -y[t];
		    }
		    Z[gi][t] = tobit->e[t] * x;
		    if (xi == tobit->k + 1) {
			Z[gi][t] += 1.0 / siginv;
		    }
#ifdef DEBUG
		    printf("%#12.5g  ", Z[gi][t]);
#endif
		}
	    }
#ifdef DEBUG
	    fputc('\n', stdout);
#endif
	}
    }

    return 0;
}

static void print_iter_info (int iter, double *theta, int m, double ll,
			     PRN *prn)
{
    int i;

    pprintf(prn, "\n*** iteration %d: theta and ll ***\n", iter);
    for (i=0; i<m; i++) {
	if (i && i % 5 == 0) pputc(prn, '\n');
	pprintf(prn, "%#12.5g ", theta[i]);
    }
    pprintf(prn, "\n    ll = %g\n\n", ll);
}

#ifdef STANDALONE /* the test program, not the plugin */

static void print_tobit_stats (double *beta, int k, double sigma, 
			       gretl_matrix *VCV, double ll)
{
    int i;

    for (i=0; i<k; i++) {
	printf("%#12.5g%#12.5g\n", beta[i],
	       sqrt(gretl_matrix_get(VCV, i, i)));
    }
    printf("%#12.5g%#12.5g\n", sigma, 
	   sqrt(gretl_matrix_get(VCV, k, k)));

    gretl_matrix_print(VCV, NULL, NULL);

    printf("Log-likelihood = %#g\n", ll);    
}

#else /* the plugin */

/* Transcribe the VCV matrix into packed triangular form */

static int make_vcv (MODEL *pmod, gretl_matrix *v)
{
    const int nv = pmod->ncoeff;
    const int nterms = nv * (nv + 1) / 2;
    double x;
    int i, j, k;

    pmod->vcv = malloc(nterms * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) return 1;  

    for (i=0; i<nv; i++) {
	for (j=0; j<=i; j++) {
	    k = ijton(i+1, j+1, nv);
	    x = gretl_matrix_get(v, i, j);
	    pmod->vcv[k] = x;
	}
    }

    return 0;
}

static int add_norm_test_to_model (MODEL *pmod, double chi2)
{
    pmod->tests = malloc(sizeof *pmod->tests);
    if (pmod->tests == NULL) return 1;

    strcpy(pmod->tests[0].type, N_("Test for normality of residual"));
    strcpy(pmod->tests[0].h_0, N_("error is normally distributed"));
    pmod->tests[0].param[0] = '\0';
    pmod->tests[0].teststat = GRETL_TEST_NORMAL_CHISQ;
    pmod->tests[0].value = chi2;
    pmod->tests[0].dfn = 2;
    pmod->tests[0].dfd = 0;
    pmod->tests[0].pvalue = chisq(chi2, 2);

    pmod->ntests = 1;

    return 0;
}

/* Taking the original OLS model as a basis, re-write the statistics
   to reflect the Tobit results.
*/

static int write_tobit_stats (MODEL *pmod, tobit_info *tobit, const double **X, 
			      double sigma, gretl_matrix *VCV)
{
    int i, t, cenc = 0;
    int offset = pmod->t1;
    const double *y = X[1];
    double chi2, ubar, udev, skew, kurt;

    for (i=0; i<tobit->k; i++) {
	pmod->coeff[i] = tobit->theta[i];
	pmod->sderr[i] = sqrt(gretl_matrix_get(VCV, i, i));
    }

    pmod->sigma = sigma;
    pmod->lnL = tobit->ll;

    pmod->ess = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	double yhat = pmod->coeff[0];

	for (i=1; i<tobit->k; i++) {
	    yhat += pmod->coeff[i] * X[i + 1][t - offset];
	}

	pmod->yhat[t] = yhat;
	pmod->uhat[t] = y[t - offset] - yhat;

	pmod->ess += pmod->uhat[t] * pmod->uhat[t]; /* Is this at all valid? */

	if (y[t - offset] == 0.0) cenc++;
    }

    /* run normality test on the untruncated uhat */
    moments(pmod->t1, pmod->t2, pmod->uhat, &ubar, &udev, &skew, &kurt, pmod->ncoeff);
    chi2 = doornik_chisq(skew, kurt, pmod->nobs); 
    add_norm_test_to_model(pmod, chi2);

    /* now truncate reported yhat, uhat */
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (pmod->yhat[t] < 0.0) {
	    pmod->yhat[t] = 0.0;
	    pmod->uhat[t] = y[t - offset];
	}
    }

    pmod->fstt = pmod->rsq = pmod->adjrsq = NADBL;

    make_vcv(pmod, VCV);

    pmod->ci = TOBIT;

    gretl_model_set_int(pmod, "censobs", cenc);

    return 0;
}

#endif /* STANDALONE */

/* Main Tobit iterative loop */

static int do_tobit (const double **Z, DATAINFO *pdinfo, MODEL *pmod,
		     PRN *prn)
{
    MODEL tmod;
    tobit_info tobit;

    /* temporary artificial dataset */
    double **tZ = NULL;
    DATAINFO *tinfo = NULL;

    const double **X;
    double *beta = NULL;
    double sigma;
    int i, j, k;
    int err = 0;

    /* for VCV manipulation */
    gretl_matrix *G = NULL;
    gretl_matrix *VCV = NULL;
    gretl_matrix *J = NULL; 
    gretl_matrix *tmp = NULL; 

    /* convergence-related stuff */
    int iters, itermax;
    double tol = 1.0e-09;
    double smallstep = 1.0e-06;
    double crit = 1.0e20;
    double stepsize = 0.25;

    err = tobit_init(&tobit, pmod);
    if (err) return 1;

    k = tobit.k;
    beta = tobit.theta;   /* point beta at the start of the theta array */

    /* set of pointers into original data */
    X = make_tobit_X(pmod, Z);
    if (X == NULL) {
	err = 1;
	goto bailout;
    }

    /* create temp dataset for the OPG regression: k+1 vars plus constant */    
    tinfo = create_new_dataset(&tZ, k + 2, tobit.n, 0);
    if (tinfo == NULL) {
	err = 1;
	goto bailout;
    }

    iters = 0;
    itermax = get_maxiter();

    while (crit > tol && iters++ < itermax && !err) {

	/* compute loglikelihood and score matrix */
	err = tobit_ll(tobit.theta, X, tZ, &tobit, 1); 

	/* BHHH via OPG regression (OPT_A -> "this is an auxiliary regression") */
	tmod = lsq(tobit.list, &tZ, tinfo, OLS, OPT_A, 0.0);

	for (i=0; i<=k; i++) {
#ifdef DEBUG
	    printf("OPG coeff[%d] = %g\n", i, tmod.coeff[i]);
#endif
	    tobit.delta[i] = tmod.coeff[i] * stepsize;
	    tobit.deltmp[i] = tobit.theta[i] + tobit.delta[i];
	} 
	
	clear_model(&tmod, NULL);

	/* see if we've gone up... (0 means "don't compute score") */
	tobit_ll(tobit.deltmp, X, tZ, &tobit, 0); 

	while (tobit.ll2 < tobit.ll && stepsize > smallstep && !err) { 
	    /* ...if not, halve steplength, as with ARMA models */
	    stepsize *= 0.5;
	    for (i=0; i<=k; i++) {
		tobit.delta[i] *= 0.5;
		tobit.deltmp[i] = tobit.theta[i] + tobit.delta[i];
	    }
	    tobit_ll(tobit.deltmp, X, tZ, &tobit, 0);
	}

	/* double the steplength? */
	if (stepsize < 4.0) stepsize *= 2.0;

	/* actually update parameter estimates */
	for (i=0; i<=k; i++) {
	    tobit.theta[i] += tobit.delta[i];
	}

	print_iter_info(iters, tobit.theta, k+1, tobit.ll, prn);

	if (tobit.theta[k] < 0.0) {
	    err = 1;
	    break;
	}

	crit = tobit.ll2 - tobit.ll;  
    }

    if (crit > tol || err != 0) {
	pmod->errcode = E_NOCONV;
	err = 1;
	goto bailout;
    }

    /* recover estimate of variance */
    sigma = 1.0 / tobit.theta[k]; 

    /* recover slope estimates */
    for (i=0; i<k; i++) {
	beta[i] *= sigma;
    }

    /* get estimate of variance matrix for Olsen parameters */
    G = gretl_matrix_from_2d_array((const double **) tZ + 1, tobit.n, tobit.k + 1);
    VCV = gretl_matrix_vcv(G);
    gretl_invert_symmetric_matrix(VCV);
    gretl_matrix_divide_by_scalar(VCV, tobit.n);
    gretl_matrix_free(G);

    /* Jacobian mat. for transforming VCV from Olsen to slopes + variance */
    J = gretl_matrix_alloc(tobit.k + 1, tobit.k + 1);
    gretl_matrix_zero(J);
    for (i=0; i<=tobit.k; i++) {
	for (j=0; j<=tobit.k; j++) {
	    if (i == j && i < tobit.k) {
		/* upper left diagonal component */
		gretl_matrix_set(J, i, j, sigma);
	    } else if (j == tobit.k && i < j) {
		/* right-hand column */
		gretl_matrix_set(J, i, j, -sigma * beta[i]);
	    } else if (j == tobit.k && i == j) {
		/* bottom right-hand element */
		gretl_matrix_set(J, i, j, -sigma * sigma);
	    }
	}
    }

    /* VCV matrix transformation */
    tmp = gretl_matrix_alloc(tobit.k + 1, tobit.k + 1);
    gretl_matrix_multiply(J, VCV, tmp);
    gretl_matrix_multiply_mod(tmp, GRETL_MOD_NONE,
			      J, GRETL_MOD_TRANSPOSE,
			      VCV);
    gretl_matrix_free(tmp);
    gretl_matrix_free(J);

#ifdef STANDALONE
    print_tobit_stats(beta, k, sigma, VCV, tobit.ll);
#else
    write_tobit_stats(pmod, &tobit, X, sigma, VCV);
#endif

    gretl_matrix_free(VCV);

 bailout:

    free_Z(tZ, tinfo);
    free_datainfo(tinfo);
    tobit_free(&tobit);
    free(X);

    return err;
}

#ifdef STANDALONE

static int *make_ols_list (int k)
{
    int *list;
    int i;

    list = malloc((k + 1) * sizeof *list);
    if (list == NULL) return NULL;

    list[0] = k;
    list[1] = 1;
    list[2] = 0;
    for (i=3; i<=list[0]; i++) list[i] = i - 1;

#ifdef DEBUG
    printlist(list, "input list");
#endif

    return list;
}

static double 
censored_frac (const double *y, const DATAINFO *pdinfo)
{
    int t, cenc = 0;

    for (t=0; t<pdinfo->n; t++) {
	if (y[t] == 0.0) cenc++;
    }
    return (double) cenc / pdinfo->n;
}

int main (int argc, char *argv[])
{
    MODEL model;
    PRN *prn;
    const char *datafile;
    int err;
    int *list;

    /* printing mechanism */
    prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);

    /* basic dataset apparatus */
    DATAINFO *pdinfo = NULL;
    double **Z = NULL;

    /* print libgretl version */
    lib_logo();

    if (argc < 2) {
	fputs("Please supply the name of a data file\n", stderr);
	exit(EXIT_FAILURE);
    }

    datafile = argv[1];

    err = import_csv(&Z, &pdinfo, datafile, NULL);
    if (err) exit(EXIT_FAILURE);

    pprintf(prn, "Perc. of censored observations = %.3f\n", 
	    censored_frac(Z[1], pdinfo));

    list = make_ols_list(pdinfo->v);
    if (list == NULL) exit(EXIT_FAILURE);

    /* run initial OLS */
    model = lsq(list, &Z, pdinfo, OLS, OPT_A, 0.0);
    if (model.errcode) {
	err = 1;
	goto bailout;
    }    

    /* do the actual analysis */
    err = do_tobit((const double **) Z, pdinfo, &model, prn); 

    /* clean up -- not really needed in a standalone program, but
       we want to be able to check for memory leaks.
    */

 bailout:

    clear_model(&model, NULL);
    free(list);
    free_Z(Z, pdinfo);
    free_datainfo(pdinfo);

    gretl_print_destroy(prn);

    return err;
}

#else

/* the driver function for the plugin */

MODEL tobit_estimate (int *list, double ***pZ, DATAINFO *pdinfo,
		      PRN *prn) 
{
    MODEL model;

    /* run initial OLS */
    model = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (model.errcode) {
	return model;
    }    

    /* do the actual Tobit analysis */
    do_tobit((const double **) *pZ, pdinfo, &model, prn); 

    return model;
}

#endif /* STANDALONE */
