#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gretl/libgretl.h>
#include <gsl/gsl_linalg.h>

static void gretl_gsl_matrix_print (gsl_matrix *X, int rows, int cols,
				    int triangle)
{
    int i, j, jmax;
    double x;
    char numstr[16];

    jmax = (triangle)? 1 : cols;

    for (i=0; i<rows; i++) {
	for (j=0; j<jmax; j++) {
	    printf("%#10.5g ", gsl_matrix_get(X, i, j));
	}
	for (j=jmax; j<cols; j++) {
	    x = gsl_matrix_get(X, i, i) * gsl_matrix_get(X, j, j);
	    x = sqrt(x);
	    x = gsl_matrix_get(X, i, j) / x;
	    sprintf(numstr,"(%.3f)", x); 
	    printf("%11s", numstr);
	}
	printf("\n");
	if (triangle && jmax < cols) jmax++;
    }
}

static void kronecker_place (gsl_matrix *X, 
			     const gsl_matrix *M,
			     int startrow, int startcol,
			     int k, double scale)
{
    int i, j;
    int row, col;
    double x;
    
    for (i=0; i<k; i++) {
	row = startrow * k + i;
	for (j=0; j<k; j++) {
	    col = startcol * k + j;
	    x = gsl_matrix_get(M, i, j);
	    gsl_matrix_set(X, row, col, x * scale);
	}
    }
}

static void make_Xi_from_Z (gsl_matrix *X, double **Z, int *list, int T)
{
    int i, t;

    for (i=2; i<=list[0]; i++) {
	for (t=0; t<T; t++) {
	    gsl_matrix_set(X, t, i-2, Z[list[i]][t]);
	}
    }
}

static int
gls_sigma_from_uhat (gsl_matrix *sigma, const gsl_matrix *e, int m, int T)
{
    int i, j, t;
    double xx;

    /* construct sigma: s_{ij} = e'_i * e_j / T  */
    for (i=0; i<m; i++) {
	for (j=0; j<m; j++) {
	    xx = 0.0;
	    for (t=0; t<T; t++) {
		xx += gsl_matrix_get(e, i, t) * gsl_matrix_get(e, j, t);
	    }
	    gsl_matrix_set (sigma, i, j, xx / T);
	}
    }

    return 0;
}

static gsl_matrix *
gls_sigma_inverse_from_uhat (const gsl_matrix *e, int m, int T)
{
    int i, j, t;
    double xx;
    gsl_matrix *sigma, *inverse;
    gsl_permutation *p;
    int sign;

    sigma = gsl_matrix_alloc (m, m);
    inverse = gsl_matrix_alloc (m, m);
    p = gsl_permutation_alloc (m);  

    /* construct sigma: s_{ij} = e'_i * e_j / T  */
    for (i=0; i<m; i++) {
	for (j=0; j<m; j++) {
	    xx = 0.0;
	    for (t=0; t<T; t++) {
		xx += gsl_matrix_get(e, i, t) * gsl_matrix_get(e, j, t);
	    }
	    gsl_matrix_set (sigma, i, j, xx / T);
	}
    }

    gsl_linalg_LU_decomp (sigma, p, &sign);
    gsl_linalg_LU_invert (sigma, p, inverse);

    gsl_matrix_free (sigma);
    gsl_permutation_free (p); 

    return inverse;
}

/* m = number of equations 
   k = number of indep vars per equation 
*/

static int ApB (const gsl_matrix * A, const gsl_matrix * B,
		gsl_matrix *C) 
{
    int err;

    err = gsl_linalg_matmult_mod (A, GSL_LINALG_MOD_TRANSPOSE,
				  B, GSL_LINALG_MOD_NONE,
				  C);
    return err;
}

static void sur_resids (MODEL *pmod, double **Z, gsl_matrix *uhat)
{
    int i, t;
    int k = pmod->ncoeff, T = pmod->nobs;
    double fit;

    for (t=0; t<T; t++) {
	fit = 0.0;
	for (i=0; i<k; i++) {
	    fit += pmod->coeff[i+1] * Z[pmod->list[i+2]][t];
	}
	pmod->yhat[t] = fit;
	pmod->uhat[t] = Z[pmod->list[1]][t] - fit;
	/* for cross-equation vcv */
	gsl_matrix_set(uhat, pmod->ID, t, pmod->uhat[t]);
    }

    pmod->ess = 0.0;
    for (t=0; t<T; t++) {
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }
    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    /* pmod->rsq = 1.0 - (pmod->ess / pmod->tss); hmm... */
}

static int calculate_coefficients (MODEL **models, double **Z,
				   gsl_matrix *X, gsl_matrix *uhat,
				   double *tmp_y, int m, int k)
{
    gsl_vector *y;
    gsl_vector *coeff;
    gsl_permutation *p; 
    gsl_matrix *vcv;
    int i, j, sign;
    int ncoeff = m * k;

    y = gsl_vector_alloc (ncoeff);
    coeff = gsl_vector_calloc (ncoeff);
    p = gsl_permutation_alloc (ncoeff);
    vcv = gsl_matrix_alloc (ncoeff, ncoeff);

    for (i=0; i<ncoeff; i++) {
	gsl_vector_set(y, i, tmp_y[i]);
    }

    gsl_linalg_LU_decomp (X, p, &sign); 
    gsl_linalg_LU_solve (X, p, y, coeff);
    gsl_linalg_LU_invert (X, p, vcv);

    for (i=0; i<m; i++) {
	for (j=0; j<k; j++) {
	    (models[i])->coeff[j+1] = gsl_vector_get(coeff, i * k + j);
	    (models[i])->sderr[j+1] = 
		sqrt(gsl_matrix_get(vcv, i * k + j, i * k + j));
	}
	sur_resids(models[i], Z, uhat);
    }

    gsl_vector_free (y);
    gsl_vector_free (coeff);
    gsl_permutation_free (p);
    gsl_matrix_free (vcv);

    return 0;
}

int sur (int m, int k, int T, int **lists, double ***pZ,
	 DATAINFO *pdinfo)
{
    int i, j, t, l;
    gsl_matrix *X, *Xi, *Xj, *M;
    gsl_matrix *uhat, *sigma;
    double *tmp_y, *y;
    int v, bigrows = m * k;
    MODEL **models;
    PRN *prn;

    prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);

    models = malloc(m * sizeof *models);
    for (i=0; i<m; i++) {
	models[i] = gretl_model_new(pdinfo);
    }

    X = gsl_matrix_alloc (bigrows, bigrows);
    Xi = gsl_matrix_alloc (T, k);
    Xj = gsl_matrix_alloc (T, k);
    M = gsl_matrix_alloc(k, k);
    uhat = gsl_matrix_alloc(m, T);

    /* first grab the OLS residuals */
    for (i=0; i<m; i++) {
	*models[i] = lsq(lists[i], pZ, pdinfo, OLS, 1, 0.0);
	if ((models[i])->errcode) {
	    fprintf(stderr, "model failed on lists[%d], code=%d\n",
		    i, (models[i])->errcode);
	    return 1;
	}
	(models[i])->ID = i;
	(models[i])->aux = AUX_SUR;
	for (t=0; t<T; t++) {
	    gsl_matrix_set(uhat, i, t, (models[i])->uhat[t]);
	}
    }

    sigma = gls_sigma_inverse_from_uhat (uhat, m, T);

    /* Xi = data matrix for equation i, specified in lists[i] */
    for (i=0; i<m; i++) {
	fprintf(stderr, "doing make_Xi_from_Z(), i=%d\n", i);
	make_Xi_from_Z(Xi, *pZ, lists[i], T);
	for (j=0; j<m; j++) { 
	    if (i != j) {
		make_Xi_from_Z(Xj, *pZ, lists[j], T);
	    }
	    ApB ((const gsl_matrix *) Xi, 
		 (i == j)? (const gsl_matrix *) Xi : 
		 (const gsl_matrix *) Xj, M);
	    kronecker_place (X, (const gsl_matrix *) M,
			     i, j, k, 
			     gsl_matrix_get(sigma, i, j)); 
	}
    }

    tmp_y = malloc((m * k) * sizeof *tmp_y);

    /* form Y column vector (m x k) */
    v = 0;
    for (i=0; i<m; i++) { /* loop over the m vertically arranged
			     blocks in the final column vector */
	double xx;

	fprintf(stderr, "working on block %d\n", i);
	make_Xi_from_Z(Xi, *pZ, lists[i], T);
	for (j=0; j<k; j++) { /* loop over the k rows within each of 
				 the m blocks */
	    fprintf(stderr, " working on row %d\n", i * k + j);
	    tmp_y[v] = 0.0;
	    for (l=0; l<m; l++) { /* loop over the m components that
				     must be added to form each element */
		fprintf(stderr, "  component %d of row %d\n", 
		       l+1, i * k + j + 1);
		fprintf(stderr, "    sigma(%d, %d) * ", i, l);
		fprintf(stderr, "X'_%d[%d] * ", i, j);
		fprintf(stderr, "y_%d\n", l);
		y = (*pZ)[lists[l][1]];
		/* multiply X'[l] into y */
		xx = 0.0;
		for (t=0; t<T; t++) {
		    xx += gsl_matrix_get(Xi, t, j) * y[t];
		}
		xx *= gsl_matrix_get(sigma, i, l);
		tmp_y[v] += xx;
	    }
	    fprintf(stderr, " finished row %d\n", i * k + j);
	    v++;
	}
	fprintf(stderr, "finished block %d\n", i);	
    }

    calculate_coefficients (models, *pZ, X, uhat, tmp_y, m, k);
    gls_sigma_from_uhat (sigma, uhat, m, T);

    for (i=0; i<m; i++) {
	printmodel(models[i], pdinfo, prn);
	free_model(models[i]);
    }

    printf("cross-equation VCV for residuals\n"
	   "(correlations above the diagonal)\n\n");
    gretl_gsl_matrix_print(sigma, m, m, 1);

    gsl_matrix_free(X);
    gsl_matrix_free(Xi);
    gsl_matrix_free(Xj);
    gsl_matrix_free(M);
    gsl_matrix_free(sigma);
    gsl_matrix_free(uhat);

    free(tmp_y);
    free(models);

    return 0;
}

static int read_grunfeld (double **Z, DATAINFO *pdinfo)
{
    FILE *fp;
    char line[128];
    int i, t, year, firm, lastfirm = 0;
    double I, F, C;

    fp = fopen("grunfeld.data", "r");
    if (fp == NULL) return 1;

    while (fgets(line, 127, fp)) {
	if (strstr(line, "Year")) continue;
	if (sscanf(line, "%d %d %lf %lf %lf", 
		   &year, &firm, &I, &F, &C) != 5) {
	    break;
	}
	if (firm != lastfirm) {
	    t = 0;
	}
	Z[1 + ((firm - 1) * 3) + 0][t] = I;
	Z[1 + ((firm - 1) * 3) + 1][t] = F;
	Z[1 + ((firm - 1) * 3) + 2][t] = C;
	lastfirm = firm;
	t++;
    }

    for (i=0; i<5; i++) {
	sprintf(pdinfo->varname[1 + 3 * i], "I_%d", i+1);
	sprintf(pdinfo->varname[2 + 3 * i], "F_%d", i+1);
	sprintf(pdinfo->varname[3 + 3 * i], "C_%d", i+1);
    }

    return 0;
}

static void print_grunfeld (double **Z) 
{
    int i, t;

    for (i=0; i<16; i++) {
	for (t=0; t<20; t++) {
	    printf("Z[%d][%d] = %7.2f ", i, t, Z[i][t]);
	}
	printf("\n");
    }
}

int main (void) 
{
    double **Z = NULL;
    int **lists = NULL;
    int m = 5;         /* number of equations */
    int k = 3;         /* indep vars per equation */
    int T = 20;        /* observations per series */
    int n = m * k + 1; /* add one for constant */
    int i, j;
    DATAINFO *pdinfo;

    pdinfo = create_new_dataset(&Z, n, T, 0);

    read_grunfeld(Z, pdinfo);
    if (0) {
	print_grunfeld(Z);
    }

    lists = malloc(m * sizeof *lists);
    for (i=0; i<m; i++) {
	lists[i] = malloc((k + 2) * sizeof **lists);
    }

    for (i=0; i<m; i++) {
	lists[i][0] = k + 1;
	for (j=1; j<=k; j++) {
	    lists[i][j] = i * 3 + j;
	}
	lists[i][k+1] = 0;
    }

    sur(m, k, T, lists, &Z, pdinfo);

    for (i=0; i<n; i++) {
	free(Z[i]);
    }
    free(Z);

    return 0;
}
