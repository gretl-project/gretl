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
#endif

/* #define DEBUG */

typedef struct _tobit_info tobit_info;

struct _tobit_info {
    int k;
    int n;
    double ll;
    double ll2;
    double *theta;
    double *ystar;
    double *P;
    double *e;
    double *f;
};

static int tobit_allocate (tobit_info *tobit)
{
    int err = 0;

    tobit->theta = NULL;
    tobit->P = tobit->e = tobit->f = tobit->ystar = NULL;

    tobit->theta = malloc((tobit->k + 1) * sizeof *tobit->theta);
    tobit->P = malloc(tobit->n * sizeof *tobit->P);
    tobit->e = malloc(tobit->n * sizeof *tobit->e);
    tobit->f = malloc(tobit->n * sizeof *tobit->f);
    tobit->ystar = malloc(tobit->n * sizeof *tobit->ystar);

    if (tobit->P == NULL || tobit->e == NULL || tobit->f == NULL ||
	tobit->ystar == NULL || tobit->theta == NULL) {
	err = 1;
	free(tobit->P);
	free(tobit->e);
	free(tobit->f);
	free(tobit->ystar);
	free(tobit->theta);
    }

    return err;
}

static void tobit_free (tobit_info *tobit)
{
    free(tobit->P);
    free(tobit->e);
    free(tobit->f);
    free(tobit->ystar);
    free(tobit->theta);
}

static int tobit_ll (double *theta, double **Z, tobit_info *tobit, 
		     int do_score)
{
    double *y = Z[1];
    double siginv = theta[tobit->k]; /* inverse of variance */
    double x, llt;
    int i, t;

    /* calculate ystar, e, f, and P vectors */
    for (t=0; t<tobit->n; t++) {
	tobit->ystar[t] = theta[0];
	for (i=1; i<tobit->k; i++) {
	    tobit->ystar[t] += theta[i] * Z[i+1][t];
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

		/* indices into the Z data array */
		gi = i + tobit->k + 1;
		xi = (i == 0)? 0 : i + 1;

		if (y[t] == 0.0) {
		    /* score if y is censored */
		    if (xi <= tobit->k) {
			Z[gi][t] = -den / tail * Z[xi][t];
		    } else {
			Z[gi][t] = 0.0;
		    } 
#ifdef DEBUG
		    printf("%#12.5g  ", Z[gi][t]);
#endif
		} else {
		    /* score if y is not censored */
		    if (xi <= tobit->k) {
			x = Z[xi][t];
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

static int *make_regression_list (int ncoeff)
{
    int *list;
    int i;

    list = malloc((ncoeff + 3) * sizeof *list);
    if (list == NULL) return NULL;

    list[0] = ncoeff + 1;
    list[1] = 1; /* dep var */
    list[2] = 0; /* intercept */

    for (i=3; i<=list[0]; i++) {
	list[i] = i - 1; /* indep vars */
    }

#ifdef DEBUG
    printlist(list, "initial OLS regression list");
#endif

    return list;
}

static void print_iter_info (double *theta, int m, double ll)
{
    int i;

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

static int do_tobit (double ***pZ, DATAINFO *pdinfo, tobit_info *tobit)
{
    MODEL tmod;
    int *list = NULL;
    double *beta = NULL, *delta = NULL, *deltmp = NULL;
    double sigma;
    int i, j, k = pdinfo->v - 1; /* minus the dep var */
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

    tobit->n = pdinfo->n;
    tobit->k = k;

    err = tobit_allocate(tobit);
    if (err) return 1;

    delta = malloc((k + 1) * sizeof *delta);
    deltmp = malloc((k + 1) * sizeof *deltmp);
    if (delta == NULL || deltmp == NULL) {
	err = 1;
	goto bailout;
    }

    list = make_regression_list(k);
    if (list == NULL) {
	err = 1;
	goto bailout;
    }

    beta = tobit->theta;   /* point beta at the start of the theta array */
    tobit->theta[k] = 1.0; /* initialize variance */

    /* initialize beta as OLS */
    tmod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (tmod.errcode) {
	err = 1;
	goto bailout;
    }
    for (i=0; i<k; i++) {
	beta[i] = tmod.coeff[i];
#ifdef DEBUG
	fprintf(stderr, "initial OLS beta[%d] = %#.5g\n", i, beta[i]);
#endif
    }
    clear_model(&tmod, NULL);

    /* add k + 1 variables for the OPG regression */
    if (dataset_add_vars(k + 1, pZ, pdinfo)) {
	err = 1;
	goto bailout;
    }

    /* OPG regression list */
    list[0] += 1; /* one element longer than for initial OLS */
    list[1] = 0;  /* dep var is the constant */
    for (i=0; i<=k; i++) {
	list[i+2] = i + k + 1;
    }

#ifdef DEBUG
    printlist(list, "OPG regression list");
#endif

    while (convcrit > small) {

	/* compute loglikelihood and score matrix */
	tobit_ll(tobit->theta, *pZ, tobit, 1); 

	/* BHHH via OPG regression (OPT_A -> "this is an auxiliary regression") */
	tmod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);

	for (i=0; i<=k; i++) {
	    delta[i] = tmod.coeff[i] * stepsize;
	    deltmp[i] = tobit->theta[i] + delta[i];
	} 
	
	clear_model(&tmod, NULL);

	/* see if we've gone up... (0 -> don't compute score) */
	tobit_ll(deltmp, *pZ, tobit, 0); 

	while (tobit->ll2 < tobit->ll && stepsize > smallstep) { 
	    /* ...if not, halve steplength, as with ARMA models */
	    stepsize *= 0.5;
	    for (i=0; i<=k; i++) {
		delta[i] *= 0.5;
		deltmp[i] = tobit->theta[i] + delta[i];
	    }
	    tobit_ll(deltmp, *pZ, tobit, 0);
	}

	/* double the steplength? */
	if (stepsize < 4.0) stepsize *= 2.0;

	/* actually update parameter estimates */
	for (i=0; i<=k; i++) {
	    tobit->theta[i] += delta[i];
	}  
                   
	print_iter_info(tobit->theta, k+1, tobit->ll);

	convcrit = tobit->ll2 - tobit->ll;                 
    }

    /* recover estimate of variance */
    sigma = 1.0 / tobit->theta[k]; 

    /* recover slope estimates */
    for (i=0; i<k; i++) {
	beta[i] *= sigma;
    }

    /* get estimate of variance matrix for Olsen parameters */
    G = gretl_matrix_from_2d_array((const double **) &((*pZ)[tobit->k+1]), 
				   tobit->n, tobit->k + 1);
    VCV = gretl_matrix_vcv(G);
    gretl_invert_symmetric_matrix(VCV);
    gretl_matrix_divide_by_scalar(VCV, tobit->n);
    gretl_matrix_free(G);

    /* Jacobian mat. for transforming VCV from Olsen to slopes + variance */
    J = gretl_matrix_alloc(tobit->k + 1, tobit->k + 1);
    gretl_matrix_zero(J);
    for (i=0; i<=tobit->k; i++) {
	for (j=0; j<=tobit->k; j++) {
	    if (i == j && i < tobit->k) {
		/* upper left diagonal component */
		gretl_matrix_set(J, i, j, sigma);
	    } else if (j == tobit->k && i < j) {
		/* right-hand column */
		gretl_matrix_set(J, i, j, -sigma * beta[i]);
	    } else if (j == tobit->k && i == j) {
		/* bottom right-hand element */
		gretl_matrix_set(J, i, j, -sigma * sigma);
	    }
	}
    }

    /* VCV matrix transformation */
    tmp = gretl_matrix_alloc(tobit->k + 1, tobit->k + 1);
    gretl_matrix_multiply(J, VCV, tmp);
    gretl_matrix_multiply_mod(tmp, GRETL_MOD_NONE,
			      J, GRETL_MOD_TRANSPOSE,
			      VCV);
    gretl_matrix_free(tmp);
    gretl_matrix_free(J);

    print_tobit_stats(beta, k, sigma, VCV, tobit->ll);

    gretl_matrix_free(VCV);

 bailout:

    free(delta);
    free(deltmp);
    free(list);

    return err;
}

#ifdef STANDALONE

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
    const char *datafile;
    int err;

    /* basic dataset apparatus */
    DATAINFO *pdinfo = NULL;
    double **Z = NULL;

    tobit_info tobit;

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

    /* do the actual analysis */
    err = do_tobit(&Z, pdinfo, &tobit); 

    /* clean up (not really needed in a standalone program, but
       necessary to avoid memory leaks when this whole function is
       embedded in a larger program).
    */

    tobit_free(&tobit);
    free_Z(Z, pdinfo);
    free_datainfo(pdinfo);

    return err;
}

#else

MODEL tobit_model (int *list, double ***pZ, DATAINFO *pdinfo, 
		   PRN *prn)
{
    MODEL tmod; 
    tobit_info tobit;
    int err;

    err = do_tobit(pZ, pdinfo, &tobit); 

    tobit_free(&tobit);

    /* to be written!! */

    gretl_model_init(&tmod, NULL);

    return tmod;
}

#endif /* STANDALONE */
