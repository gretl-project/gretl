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
    if (!*err) {
	mi = gretl_array_get_element(A, 1, NULL, err);
    }

    if (!*err) {
	n = mr->rows;
	if (n == 0 || mr->cols != n || mi->rows != n || mi->cols != n) {
	    *err = E_NONCONF;
	} else {
	    ret = gretl_matrix_alloc(n, 1);
	    work = malloc(sizeof *work);
	    a = malloc(n * n * sizeof *a);
	    if (ret == NULL || work == NULL || a == NULL) {
		*err = E_ALLOC;
	    }
	}
    }

    if (*err) {
	goto bailout;
    }

    w = ret->val;

    /* write upper triangle of complex matrix into @a */
    for (i=0; i<n; i++) {
	k = i * n;
	for (j=0; j<=i; j++) {
	    a[k].r = gretl_matrix_get(mr, i, j);
	    a[k].i = gretl_matrix_get(mi, i, j);
	    k++;
	}
    }

    /* get optimal workspace size */
    lwork = -1;
    zheev_(&jobz, &uplo, &n, a, &n, w, work, &lwork, rwork, &info);

    lwork = (integer) work[0].r;
    work = realloc(work, lwork * sizeof *work);
    rwork = malloc((3 * n - 2) * sizeof *rwork);
    if (work == NULL || rwork == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    /* do the actual eigen decomposition */
    zheev_(&jobz, &uplo, &n, a, &n, w, work, &lwork, rwork, &info);
    if (info != 0) {
	fprintf(stderr, "zheev: info = %d\n", info);
	*err = E_DATA;
    } else if (V != NULL) {
	mr = gretl_matrix_alloc(n, n);
	mi = gretl_matrix_alloc(n, n);

	if (mr == NULL || mi == NULL) {
	    *err = E_ALLOC;
	    goto bailout;
	}
	
	/* read out the eigenvectors */
	for (i=0; i<n; i++) {
	    k = i * n;
	    for (j=0; j<n; j++) {
		gretl_matrix_set(mr, j, i, a[k].r);
		gretl_matrix_set(mi, j, i, a[k].i);
		k++;
	    }
	}

	/* note: array takes ownership of mr, mi */
	gretl_array_set_matrix(V, 0, mr, 0);
	gretl_array_set_matrix(V, 1, mi, 0);
    }

 bailout:

    free(rwork);
    free(work);
    free(a);

    if (*err) {
	gretl_matrix_free(ret);
	ret = NULL;
    }

    return ret;
}
