#include <stdio.h>
#include <stdlib.h>

#include <clapack.h>
#include <f2c.h>

static double *lapack_eigenvals (double **X, int k) 
{
    int n = k;
    int info, sdim;
    int one = 1;
    int i, j, p;
    int m = k * k;
    int lwork = m;

    double *a;
    double *work;
    double *wr, *wi;
    char job = 'N', sort = 'N';

    a = malloc(m * sizeof *a);
    if (a == NULL) return NULL;

    work = malloc(m * sizeof *work);
    if (work == NULL) {
	free(a);
	return NULL;
    }

    wr = malloc(n * sizeof *wr);
    if (wr == NULL) {
	free(a);
	free(work);
	return NULL;
    }

    wi = malloc(n * sizeof *wi);
    if (wi == NULL) {
	free(a);
	free(work);
	free(wr);
	return NULL;
    }

    p = 0;
    for (j=0; j<k; j++) {
	for (i=0; i<k; i++) {
	    a[p++] = X[i][j];
	}
    } 

    dgees_(&job, &sort, NULL, &n, 
	   a, &n, &sdim, wr, wi, NULL, &one, 
	   work, &lwork, NULL, &info);

    if (info != 0) {
	free(wr);
	wr = NULL;
    }

    free(a);
    free(wi);
    free(work);

    return wr;
}


