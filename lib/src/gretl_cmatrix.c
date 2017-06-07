/* complex matrices */

#include "clapack_complex.h"
#include "gretl_cmatrix.h"

static int get_two_matrices (gretl_array *A,
			     gretl_matrix **pmr,
			     gretl_matrix **pmi,
			     int square)
{
    gretl_matrix *mr = NULL, *mi = NULL;
    int err = 0;
    
    mr = gretl_array_get_element(A, 0, NULL, &err);
    if (!err) {
	mi = gretl_array_get_element(A, 1, NULL, &err);
    }

    if (!err) {
	int r = mr->rows;
	int c = mr->cols;

	if (r == 0 || c == 0 || mi->rows != r || mi->cols != c) {
	    err = E_NONCONF;
	}
	if (!err && square && r != c) {
	    err = E_NONCONF;
	}
    }

    if (!err) {
	*pmr = mr;
	*pmi = mi;
    }

    return err;
}

static int complex_mat_into_array (cmplx *a, int n,
				   gretl_array *A)
{
    gretl_matrix *mr = gretl_matrix_alloc(n, n);
    gretl_matrix *mi = gretl_matrix_alloc(n, n);
    int i, j, k;
    int err = 0;

    if (mr == NULL || mi == NULL) {
	return E_ALLOC;
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
    gretl_array_set_matrix(A, 0, mr, 0);
    gretl_array_set_matrix(A, 1, mi, 0);    

    return err;
}

/* eigenvalues and optionally eigenvectors of a Hermitian matrix */

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

    *err = get_two_matrices(A, &mr, &mi, 1);

    if (!*err) {
	n = mr->rows;
	ret = gretl_matrix_alloc(n, 1);
	work = malloc(sizeof *work);
	a = malloc(n * n * sizeof *a);
	if (ret == NULL || work == NULL || a == NULL) {
	    *err = E_ALLOC;
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
	*err = complex_mat_into_array(a, n, V);
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

/* inverse of complex matrix via LU decomposition;
   uses zgetrf(), zgetri() */

gretl_array *gretl_zgetri (gretl_array *A, int *err)
{
    gretl_array *Ainv = NULL;
    gretl_matrix *mr = NULL;
    gretl_matrix *mi = NULL;
    integer lwork;
    integer *ipiv;
    cmplx *work, *a;
    integer n, info;
    int i, j, k;

    *err = get_two_matrices(A, &mr, &mi, 1); /* square? */
    if (*err) {
	return NULL;
    }

    n = mr->rows;
    lwork = 10 * n;

    a = malloc(n * n *sizeof *a);
    ipiv = malloc(2 * n *sizeof *ipiv);
    work = malloc(lwork * sizeof *work);
    if (a == NULL || ipiv == NULL || work == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    /* write complex matrix into @a */
    for (i=0; i<n; i++) {
	k = i * n;
	for (j=0; j<n; j++) {
	    a[k].r = gretl_matrix_get(mr, i, j);
	    a[k].i = gretl_matrix_get(mi, i, j);
	    k++;
	}
    }    

    zgetrf_(&n, &n, a, &n, ipiv, &info);
    if (info != 0) {
	printf("zgetrf: info = %d\n", info);
	*err = E_DATA;
    }

    if (!*err) {
	zgetri_(&n, a, &n, ipiv, work, &lwork, &info);
	if (info != 0) {
	    printf("zgetri: info = %d\n", info);
	    *err = E_DATA;
	}
    }

    if (!*err) {
	Ainv = gretl_array_new(GRETL_TYPE_MATRICES, 2, err);
	if (!*err) {
	    *err = complex_mat_into_array(a, n, Ainv);
	}
    }

 bailout:
    
    free(work);
    free(ipiv);
    free(a);

    if (*err && Ainv != NULL) {
	gretl_array_destroy(Ainv);
	Ainv = NULL;
    }

    return Ainv;
}
