/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
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

/* The algorithm here was contributed by Riccardo "Jack" Lucchetti. */

#ifdef STANDALONE
# include <gretl/libgretl.h>
#else
# include "libgretl.h"
# include "internal.h"
#endif

/* #define DEBUG */

typedef struct _tobit_info tobit_info;

struct _tobit_info {
    int k;
    int n;
    double ll;
    double ll2;
    double *theta;
    double *delta;
    double *deltmp;
    double *ystar;
    double *P;
    double *e;
    double *f;
};

static int tobit_init (tobit_info *tobit, const MODEL *pmod)
{
    int i, k, err = 0;

    tobit->n = pmod->nobs;
    k = tobit->k = pmod->ncoeff;    

    tobit->theta = NULL;
    tobit->P = tobit->e = tobit->f = tobit->ystar = NULL;

    tobit->theta = malloc((tobit->k + 1) * sizeof *tobit->theta);
    tobit->delta = malloc((tobit->k + 1) * sizeof *tobit->delta);
    tobit->deltmp = malloc((tobit->k + 1) * sizeof *tobit->deltmp);
    tobit->P = malloc(tobit->n * sizeof *tobit->P);
    tobit->e = malloc(tobit->n * sizeof *tobit->e);
    tobit->f = malloc(tobit->n * sizeof *tobit->f);
    tobit->ystar = malloc(tobit->n * sizeof *tobit->ystar);

    if (tobit->P == NULL || tobit->e == NULL || tobit->f == NULL ||
	tobit->ystar == NULL || tobit->theta == NULL || 
	tobit->delta == NULL || tobit->deltmp == NULL) {
	err = 1;
	free(tobit->P);
	free(tobit->e);
	free(tobit->f);
	free(tobit->ystar);
	free(tobit->theta);
	free(tobit->delta);
	free(tobit->deltmp);
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

static void tobit_free (tobit_info *tobit)
{
    free(tobit->P);
    free(tobit->e);
    free(tobit->f);
    free(tobit->ystar);
    free(tobit->theta);
    free(tobit->delta);
    free(tobit->deltmp);
}

static int tobit_ll (double *theta, const double **X, double **Z, tobit_info *tobit, 
		     int do_score)
{
    const double *y = X[1];
    double siginv = theta[tobit->k]; /* inverse of variance */
    double x, llt;
    int i, t;

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

		/* indices into the data arrays */
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

/* construct a set of pointers into the original data set */

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

    return list;
}

#ifdef STANDALONE

static void print_iter_info (double *theta, int m, double ll)
{
    int i;

    printf("\n*** theta and ll ***\n");
    for (i=0; i<m; i++) {
	printf("%#12.5g ", theta[i]);
    }
    printf("%#12.5g\n\n", ll);
}

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

#else

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

static int write_tobit_stats (MODEL *pmod, tobit_info *tobit, const double *y,
			      double sigma, gretl_matrix *VCV)
{
    int i, t;
    int offset = pmod->t1;

    for (i=0; i<tobit->k; i++) {
	pmod->coeff[i] = tobit->theta[i];
	pmod->sderr[i] = sqrt(gretl_matrix_get(VCV, i, i));
    }

    pmod->sigma = sigma;
    pmod->lnL = tobit->ll;

    pmod->ess = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->yhat[t] = tobit->ystar[t - offset];
	pmod->uhat[t] = y[t - offset] - pmod->yhat[t];
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }

    if (pmod->tss > pmod->ess) {
	pmod->fstt = pmod->dfd * (pmod->tss - pmod->ess) / 
	    (pmod->dfn * pmod->ess);
    } else {
	pmod->fstt = NADBL;
    }

    pmod->rsq = pmod->adjrsq = NADBL;

    if (pmod->tss > 0) {
	pmod->rsq = 1.0 - (pmod->ess / pmod->tss);
	if (pmod->dfd > 0) {
	    double den = pmod->tss * pmod->dfd;

	    pmod->adjrsq = 1.0 - (pmod->ess * (pmod->nobs - 1) / den);
	}
    }

    make_vcv(pmod, VCV);

    gretl_aic_etc(pmod);

    pmod->ci = TOBIT;

    return 0;
}

#endif /* STANDALONE */

static int do_tobit (const double **Z, DATAINFO *pdinfo, MODEL *pmod)
{
    MODEL tmod;
    tobit_info tobit;

    double **tZ = NULL;
    DATAINFO *tinfo = NULL;

    const double **X;
    int *tlist = NULL;
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
    double small = 1.0e-09;
    double smallstep = 1.0e-06;
    double convcrit = 1.0e20;
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

    /* create a temp dataset for the OPG regression: 
       k+1 vars plus constant */    
    tinfo = create_new_dataset(&tZ, k + 2, tobit.n, 0);
    if (tinfo == NULL) {
	err = 1;
	goto bailout;
    }

    /* OPG regression list */
    tlist = make_tobit_list(k + 1);
    if (tlist == NULL) {
	err = 1;
	goto bailout;
    }

    printlist(tlist, "OPG regression list");

    while (convcrit > small) {

	/* compute loglikelihood and score matrix */
	tobit_ll(tobit.theta, X, tZ, &tobit, 1); 

	/* BHHH via OPG regression (OPT_A -> "this is an auxiliary regression") */
	tmod = lsq(tlist, &tZ, tinfo, OLS, OPT_A, 0.0);

	for (i=0; i<=k; i++) {
#ifdef DEBUG
	    printf("OPG coeff[%d] = %g\n", i, tmod.coeff[i]);
#endif
	    tobit.delta[i] = tmod.coeff[i] * stepsize;
	    tobit.deltmp[i] = tobit.theta[i] + tobit.delta[i];
	} 
	
	clear_model(&tmod, NULL);

	/* see if we've gone up... (0 -> don't compute score) */
	tobit_ll(tobit.deltmp, X, tZ, &tobit, 0); 

	while (tobit.ll2 < tobit.ll && stepsize > smallstep) { 
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

#ifdef STANDALONE                   
	print_iter_info(tobit.theta, k+1, tobit.ll);
#endif

	convcrit = tobit.ll2 - tobit.ll;                 
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
    write_tobit_stats(pmod, &tobit, X[1], sigma, VCV);
#endif

    gretl_matrix_free(VCV);

 bailout:

    free(tlist);
    free_Z(tZ, tinfo);
    free_datainfo(tinfo);
    tobit_free(&tobit);
    free(X);

    return err;
}

#ifdef STANDALONE

static int *make_list (int k)
{
    int *list;
    int i;

    list = malloc((k + 1) * sizeof *list);
    if (list == NULL) return NULL;

    list[0] = k;
    list[1] = 1;
    list[2] = 0;
    for (i=3; i<=list[0]; i++) list[i] = i - 1;

    printlist(list, "input list");

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
    const char *datafile;
    int err;
    int *list;

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

    printf("Perc. of censored observations = %.3f\n", 
	   censored_frac(Z[1], pdinfo));

    list = make_list(pdinfo->v);
    if (list == NULL) exit(EXIT_FAILURE);

    /* run initial OLS */
    model = lsq(list, &Z, pdinfo, OLS, OPT_A, 0.0);
    if (model.errcode) {
	err = 1;
	goto bailout;
    }    

    /* do the actual analysis */
    err = do_tobit((const double **) Z, pdinfo, &model); 

    /* clean up (not really needed in a standalone program, but
       necessary to avoid memory leaks when this whole function is
       embedded in a larger program).
    */

 bailout:

    free(list);
    free_Z(Z, pdinfo);
    free_datainfo(pdinfo);

    return err;
}

#else

MODEL tobit_estimate (int *list, double ***pZ, DATAINFO *pdinfo) 
{
    MODEL model;

    /* run initial OLS */
    model = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (model.errcode) {
	return model;
    }    

    /* do the actual analysis */
    do_tobit((const double **) *pZ, pdinfo, &model); 

    return model;
}

#endif /* STANDALONE */
