#include <stdio.h>
#include <stdlib.h>

#include <f2c.h> 
#include <cblas.h>
#include <clapack_double.h>

enum {
    GRETL_MATRIX_OK = 0,
    GRETL_MATRIX_NOMEM,
    GRETL_MATRIX_NON_CONFORM,
    GRETL_MATRIX_RANGE
} gretl_matrix_errors;

enum {
    GRETL_MMOD_NONE = 0,
    GRETL_MMOD_TRANSPOSE = 1
} gretl_matrix_mods;

typedef struct _gretl_matrix gretl_matrix;

struct _gretl_matrix {
    int rows;
    int cols;
    double *data;
};

static const char *wspace_fail = "Workspace query failed\n";

static gretl_matrix *gretl_matrix_new (int rows, int cols)
{
    gretl_matrix *m;

    m = malloc(sizeof *m);
    if (m == NULL) return m;

    m->data = malloc(rows * cols * sizeof *m->data);

    if (m->data == NULL) {
	free(m);
	return NULL;
    }

    m->rows = rows;
    m->cols = cols;

    return m;
}

static void gretl_matrix_free (gretl_matrix *m)
{
    if (m == NULL) return;

    free(m->data);
    free(m);
}

static double gretl_matrix_get (gretl_matrix *m, int i, int j)
{
    if (m == NULL || m->data == NULL) return -999.0;

    if (i >= m->rows || j >= m->cols) return -999.0;

    return m->data[i * m->cols + j];
}

static gretl_matrix *gretl_matrix_from_square_array (const double **X, int k)
{
    int i, j, p;
    gretl_matrix *m;

    m = gretl_matrix_new(k, k);
    if (m == NULL) return m;

    p = 0;
    for (j=0; j<k; j++) {
	for (i=0; i<k; i++) {
	    m->data[p++] = X[i][j];
	}
    } 

    return m;
}

static int gretl_matrix_multiply (gretl_matrix *a, gretl_matrix *b, 
				  int bflag,
				  gretl_matrix *c)
{
    int i, j, k;

    if (a->cols != b->rows) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    if (c->cols != a->cols || c->rows != b->rows) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    for (i=0; i<c->rows; i++) {
	for (j=0; j<c->cols; j++) {
	    c->data[i * c->cols + j] = 0.0;
	    for (k=0; k<c->cols; k++) {
		if (bflag == GRETL_MMOD_TRANSPOSE) {
		    c->data[i * c->cols + j] += 
			a->data[i * a->cols + k] * b->data[j * b->rows + k];
		} else {
		    c->data[i * c->cols + j] += 
			a->data[i * a->cols + k] * b->data[k * b->cols + j];
		}
	    }
	}
    }

    return GRETL_MATRIX_OK;
}

static int lapack_invert_gretl_matrix (gretl_matrix *m)
{
    integer n = m->rows;    
    integer info;
    integer lwork;
    integer *ipiv;

    double *work;

    ipiv = malloc(n * sizeof *ipiv);
    if (ipiv == NULL) {
	return 1;
    }

    work = malloc(sizeof *work);
    if (work == NULL) {
	free(ipiv);
	return 1;
    }    

    dgetrf_(&n, &n, m->data, &n, ipiv, &info);   

    if (info != 0) {
	free(ipiv);
	return info;
    }

    lwork = -1;
    dgetri_(&n, m->data, &n, ipiv, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	fputs(wspace_fail, stderr);
	free(ipiv);
	return 1;
    }    

    work = realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	free(ipiv);
	return 1;
    }  

    dgetri_(&n, m->data, &n, ipiv, work, &lwork, &info);

    return info;
}

static double *gretl_lapack_eigenvals (gretl_matrix *m) 
{
    integer n = m->rows;
    integer info, sdim;
    integer lwork;
    integer one = 1;

    double *work;
    double *wr, *wi;

    char job = 'N', sort = 'N';

    work = malloc(sizeof *work);
    if (work == NULL) {
	return NULL;
    }

    wr = malloc(n * sizeof *wr);
    if (wr == NULL) {
	free(work);
	return NULL;
    }

    wi = malloc(n * sizeof *wi);
    if (wi == NULL) {
	free(work);
	free(wr);
	return NULL;
    }

    lwork = -1; /* find optimal workspace size */
    dgees_(&job, &sort, NULL, &n, 
	   m->data, &n, &sdim, wr, wi, NULL, &one, 
	   work, &lwork, NULL, &info);

    if (info != 0 || work[0] <= 0.0) {
	fputs(wspace_fail, stderr);
	free(work);
	free(wr);
	free(wi);
	return NULL;
    }	

    lwork = (integer) work[0];

    work = realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	free(wr);
	free(wi);
	return NULL;
    }    

    dgees_(&job, &sort, NULL, &n, 
	   m->data, &n, &sdim, wr, wi, NULL, &one, 
	   work, &lwork, NULL, &info);

    if (info != 0) {
	free(wr);
	wr = NULL;
    }

    free(wi);
    free(work);

    return wr;
}

static void print_M (gretl_matrix *M, int k)
{
    int i, j;
    double x;

    for (i=0; i<k; i++) {
	for (j=0; j<k; j++) {
	   x = gretl_matrix_get(M, i, j);
	   printf("%#13.7g ", x);
	}
	putchar('\n');
    } 
    putchar('\n');
}

int johansen_eigenvals (const double **X, const double **Y, const double **Z, 
			int k, double *evals)
{
    gretl_matrix *Suu, *Svv, *Suv;
    gretl_matrix *Inv, *TmpL, *TmpR, *M;
    int err = 0;

    Suu = gretl_matrix_from_square_array(X, k);
    Svv = gretl_matrix_from_square_array(Y, k);
    Suv = gretl_matrix_from_square_array(Z, k);

    Inv = gretl_matrix_new(k, k);
    TmpL = gretl_matrix_new(k, k);
    TmpR = gretl_matrix_new(k, k);
    M = gretl_matrix_new(k, k);

    /* calculate Suu^{-1} Suv */
    lapack_invert_gretl_matrix(Suu);
    gretl_matrix_multiply(Suu, Suv, GRETL_MMOD_NONE, TmpR);

    /* calculate Svv^{-1} Suv' */
    lapack_invert_gretl_matrix(Svv);
    gretl_matrix_multiply(Svv, Suv, GRETL_MMOD_TRANSPOSE, TmpR);

    gretl_matrix_multiply(TmpL, TmpR, GRETL_MMOD_NONE, M);

    print_M(M, k);

    evals = gretl_lapack_eigenvals(M);

    /* free stuff */
    gretl_matrix_free(Svv);
    gretl_matrix_free(Suu);
    gretl_matrix_free(Suv);

    gretl_matrix_free(Inv);
    gretl_matrix_free(TmpL);
    gretl_matrix_free(TmpR);
    gretl_matrix_free(M);

    return err;
}


#if 1

int main (void)
{
    double **X;
    double *evals;
    gretl_matrix *m;
    gretl_matrix *a, *b, *c;
    int i;

    X = malloc(3 * sizeof *X);
    for (i=0; i<3; i++) {
	X[i] = malloc(3 * sizeof **X);
    }

    X[0][0] = -0.2482051;
    X[0][1] = 0.3580177;
    X[0][2] = -0.6679180;
    X[1][0] = 0.03675934;
    X[1][1] = 0.006330003;
    X[1][2] = 0.07293056;
    X[2][0] = 0.1751283;
    X[2][1] = -0.2075322;
    X[2][2] = 0.4387614;

    m = gretl_matrix_from_square_array((const double **) X, 3);

    evals = gretl_lapack_eigenvals(m);

    if (evals == NULL) {
	fprintf(stderr, "eigenvals returned NULL\n");
    } else {
        for (i=0; i<3; i++) {
            printf("lambda[%d] = %g\n", i + 1, evals[i]);
        }
    }

    gretl_matrix_free(m);

    a = gretl_matrix_new (3, 3);
    b = gretl_matrix_new (3, 3);
    c = gretl_matrix_new (3, 3);

    gretl_matrix_multiply(a, b, GRETL_MMOD_NONE, c);

    return 0;
}

#endif


