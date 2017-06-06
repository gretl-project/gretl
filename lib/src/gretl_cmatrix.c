/* complex matrices */

#include "clapack_complex.h"

gretl_matrix *gretl_zheev (gretl_array *A, gretl_array *V, int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *mr, *mi;
    integer n, info, lwork;
    double *w = NULL;
    double *rwork = NULL;
    cmplx *a = NULL;
    cmplx *work = NULL;
    char jobz = V != NULL ? 'V' : 'N';
    char uplo = 'U';
    int i, j, k;

    *err = 0;

    mr = gretl_array_get_element(A, 0, NULL, err);
    mi = gretl_array_get_element(A, 1, NULL, err);
    n = mr->rows;

    ret = gretl_matrix_alloc(n, 1);
    w = ret->val;
    work = malloc(sizeof *work);
    a = malloc(n * n * sizeof *a);

    /* write upper triangle of complex matrix */
    for (i=0; i<n; i++) {
	k = i * n;
	for (j=0; j<=i; j++) {
	    a[k].r = gretl_matrix_get(mr, i, j);
	    a[k].i = gretl_matrix_get(mi, i, j);
	    k++;
	}
    }

    lwork = -1;
    zheev_(&jobz, &uplo, &n, a, &n, w, work, &lwork, rwork, &info);

    lwork = (integer) work[0].r;
    work = realloc(work, lwork * sizeof *work);
    rwork = malloc((3 * n - 2) * sizeof *rwork);

    zheev_(&jobz, &uplo, &n, a, &n, w, work, &lwork, rwork, &info);
    if (info != 0) {
	fprintf(stderr, "zheev: info = %d\n", info);
	*err = 1;
    } else if (V != NULL) {
	mr = gretl_matrix_alloc(n, n);
	mi = gretl_matrix_alloc(n, n);
	
	/* read out the eigenvectors */
	for (i=0; i<n; i++) {
	    k = i * n;
	    for (j=0; j<n; j++) {
		gretl_matrix_set(mr, j, i, a[k].r);
		gretl_matrix_set(mi, j, i, a[k].i);
		k++;
	    }
	}

	gretl_array_set_matrix(V, 0, mr, 0);
	gretl_array_set_matrix(V, 1, mi, 0);
    }	

    free(rwork);
    free(work);
    free(a);

    return ret;
}
