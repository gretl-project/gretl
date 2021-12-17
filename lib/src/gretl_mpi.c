/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "libgretl.h"
#include "gretl_mpi.h"
#include "gretl_typemap.h"
#include <mpi.h>

#ifdef WIN32
# include <windows.h>
#else
# include <dlfcn.h>
#endif

/* Support for MPI in libgretl. On systems other than MS Windows
   We get the MPI symbols that we need from the address space of
   the calling program to avoid introducing a hard dependency on
   libmpi; on Windows we call LoadLibrary() on msmpi.dll to the
   same effect.

   To use functions in this translation unit from elsewhere in
   libgretl one must first call gretl_MPI_init(), and then
   guard subsequent calls with "if (gretl_mpi_initialized())".
   Otherwise you'll get a segfault.

   For future reference, the MPI variants define the following
   specific symbols in mpi.h:

     Open MPI: OMPI_MAJOR_VERSION
     MPICH:    MPICH_VERSION
     MS-MPI:   MSMPI_VER
*/

#ifdef OMPI_MAJOR_VERSION
/* external constants that need to be loaded from the
   Open MPI library */
static MPI_Comm mpi_comm_world;
static MPI_Datatype mpi_double;
static MPI_Datatype mpi_int;
static MPI_Datatype mpi_unsigned;
static MPI_Datatype mpi_byte;
static MPI_Op mpi_sum;
static MPI_Op mpi_prod;
static MPI_Op mpi_max;
static MPI_Op mpi_min;
#else
/* It seems that MPICH and MS-MPI just define these symbols
   as integer values in the header */
# define mpi_comm_world MPI_COMM_WORLD
# define mpi_double     MPI_DOUBLE
# define mpi_int        MPI_INT
# define mpi_unsigned   MPI_UNSIGNED
# define mpi_byte       MPI_BYTE
# define mpi_sum        MPI_SUM
# define mpi_prod       MPI_PROD
# define mpi_max        MPI_MAX
# define mpi_min        MPI_MIN
#endif

enum {
    TAG_MATRIX_INFO = 1,
    TAG_MATRIX_VAL,
    TAG_SCALAR_VAL,
    TAG_INT_VAL,
    TAG_ARRAY_LEN,
    TAG_SERIES_LEN,
    TAG_SERIES_VAL,
    TAG_LIST_LEN,
    TAG_LIST_VAL,
    TAG_STR_LEN,
    TAG_STR_VAL,
    TAG_BMEMB_INFO,
    TAG_ARRAY_INFO,
    TAG_BUNDLE_SIZE
};

#define MI_LEN 5 /* matrix info length */

#define MPI_DEBUG 0

static void *MPIhandle;       /* handle to the MPI library */
static int gretl_MPI_err;     /* initialization error record */
static int gretl_MPI_initted; /* are we initialized or not? */

/* renamed, pointerized versions of the MPI functions we need */
static int (*mpi_comm_rank) (MPI_Comm, int *);
static int (*mpi_comm_size) (MPI_Comm, int *);
static int (*mpi_error_class) (int, int *);
static int (*mpi_error_string) (int, char *, int *);
static int (*mpi_reduce) (void *, void *, int, MPI_Datatype, MPI_Op,
			  int, MPI_Comm);
static int (*mpi_allreduce) (void *, void *, int, MPI_Datatype, MPI_Op,
			     MPI_Comm);
static int (*mpi_bcast) (void *, int, MPI_Datatype, int, MPI_Comm);
static int (*mpi_send) (void *, int, MPI_Datatype, int, int, MPI_Comm);
static int (*mpi_recv) (void *, int, MPI_Datatype, int, int, MPI_Comm,
			MPI_Status *);
static int (*mpi_barrier) (MPI_Comm);
static int (*mpi_probe) (int, int, MPI_Comm, MPI_Status *);
static double (*mpi_wtime) (void);
static int (*mpi_initialized) (int *);

static int gretl_matrix_bcast (gretl_matrix **pm, int id, int root);
static int gretl_array_bcast (gretl_array **pa, int id, int root);
static int gretl_bundle_send (gretl_bundle *b, int dest);
static gretl_bundle *gretl_bundle_receive (int source, int *err);

static void *mpi_receive_element (int source, GretlType etype,
				  int *err);

static void *mpiget (void *handle, const char *name, int *err)
{
#ifdef WIN32
    void *p = GetProcAddress(handle, name);
#else
    void *p = dlsym(handle, name);
#endif

    if (p == NULL) {
	printf("mpiget: couldn't find '%s'\n", name);
	*err += 1;
    }

#if MPI_DEBUG
    else {
	printf("mpiget: '%s' -> %p\n", name, p);
    }
#endif

    return p;
}

int gretl_MPI_init (void)
{
    int err = 0;

    if (gretl_MPI_initted) {
	return gretl_MPI_err;
    }

#if MPI_DEBUG
    printf("Loading MPI symbols...\n");
#endif

#ifdef WIN32
    MPIhandle = LoadLibrary("msmpi.dll");
#else
    MPIhandle = dlopen(NULL, RTLD_NOW);
#endif

    if (MPIhandle == NULL) {
	err = E_EXTERNAL;
	goto bailout;
    }

    mpi_comm_rank    = mpiget(MPIhandle, "MPI_Comm_rank", &err);
    mpi_comm_size    = mpiget(MPIhandle, "MPI_Comm_size", &err);
    mpi_error_class  = mpiget(MPIhandle, "MPI_Error_class", &err);
    mpi_error_string = mpiget(MPIhandle, "MPI_Error_string", &err);
    mpi_reduce       = mpiget(MPIhandle, "MPI_Reduce", &err);
    mpi_allreduce    = mpiget(MPIhandle, "MPI_Allreduce", &err);
    mpi_bcast        = mpiget(MPIhandle, "MPI_Bcast", &err);
    mpi_send         = mpiget(MPIhandle, "MPI_Send", &err);
    mpi_recv         = mpiget(MPIhandle, "MPI_Recv", &err);
    mpi_probe        = mpiget(MPIhandle, "MPI_Probe", &err);
    mpi_barrier      = mpiget(MPIhandle, "MPI_Barrier", &err);
    mpi_wtime        = mpiget(MPIhandle, "MPI_Wtime", &err);
    mpi_initialized  = mpiget(MPIhandle, "MPI_Initialized", &err);

#ifdef OMPI_MAJOR_VERSION
    if (!err) {
	mpi_comm_world = (MPI_Comm) mpiget(MPIhandle, "ompi_mpi_comm_world", &err);
	mpi_double     = (MPI_Datatype) mpiget(MPIhandle, "ompi_mpi_double", &err);
	mpi_int        = (MPI_Datatype) mpiget(MPIhandle, "ompi_mpi_int", &err);
	mpi_unsigned   = (MPI_Datatype) mpiget(MPIhandle, "ompi_mpi_unsigned", &err);
	mpi_byte       = (MPI_Datatype) mpiget(MPIhandle, "ompi_mpi_byte", &err);
	mpi_sum        = (MPI_Op) mpiget(MPIhandle, "ompi_mpi_op_sum", &err);
	mpi_prod       = (MPI_Op) mpiget(MPIhandle, "ompi_mpi_op_prod", &err);
	mpi_max        = (MPI_Op) mpiget(MPIhandle, "ompi_mpi_op_max", &err);
	mpi_min        = (MPI_Op) mpiget(MPIhandle, "ompi_mpi_op_min", &err);
    }
#endif

    if (err) {
	close_plugin(MPIhandle);
	MPIhandle = NULL;
	err = E_EXTERNAL;
    }

 bailout:

#if MPI_DEBUG
    printf("load_MPI_symbols: returning %d\n", err);
#endif

    gretl_MPI_err = err;
    gretl_MPI_initted = 1;

    return err;
}

static void gretl_mpi_error (int *err)
{
    char msg[BUFSIZ];
    int len, id = 0;

    mpi_comm_rank(mpi_comm_world, &id);
    mpi_error_string(*err, msg, &len);
    gretl_errmsg_sprintf("%3d: %s", id, msg);

    *err = E_EXTERNAL;
}

static int dim_error (int *dims, int n)
{
    int i, d0 = 0;
    int err = 0;

    /* allow for the possibility that some nodes
       carry null matrices: d0 will be the common
       non-zero dimension */

    for (i=0; i<n; i++) {
	if (dims[i] > 0) {
	    d0 = dims[i];
	    break;
	}
    }

    /* check that all non-null matrices are conformable */

    for (i=0; i<n; i++) {
	if (dims[i] != d0 && dims[i] != 0) {
	    err = E_NONCONF;
	    break;
	}
    }

    return err;
}

static void set_dim_max (int *dims, int *dtest, int n)
{
    int i;

    /* find the max row dimension of matrices with
       cols > 0, or vice versa */

    for (i=0; i<n; i++) {
	if (dtest[i] > 0) {
	    if (dims[i] > dims[n]) {
		dims[n] = dims[i];
	    }
	}
    }
}

static int matrix_dims_check (int *rows, int *cols, int n,
			      Gretl_MPI_Op op)
{
    int match_rows =
	(op == GRETL_MPI_SUM || op == GRETL_MPI_PROD ||
	 op == GRETL_MPI_HCAT);
    int match_cols =
	(op == GRETL_MPI_SUM || op == GRETL_MPI_PROD ||
	 op == GRETL_MPI_VCAT);
    int err = 0;

    if (match_rows) {
	err = dim_error(rows, n);
    }

    if (!err && match_cols) {
	err = dim_error(cols, n);
    }

    if (!err) {
	set_dim_max(rows, cols, n);
	set_dim_max(cols, rows, n);
	/* for now we'll consider a null matrix reduction
	   as an error (though _maybe_ we want to permit
	   this?)
	*/
	if (rows[n] == 0 || cols[n] == 0) {
	    err = E_DATA;
	}
    }

    return err;
}

static void fill_matrix_info (int *rc, const gretl_matrix *m)
{
    rc[0] = m->rows;
    rc[1] = m->cols;
    rc[2] = m->is_complex;
    rc[3] = gretl_matrix_get_t1(m);
    rc[4] = gretl_matrix_get_t2(m);
}

static int matrix_reduce_alloc (int *rows, int *cols,
				int n, Gretl_MPI_Op op,
				gretl_matrix **pm,
				double **px)
{
    int rtotal = 0, ctotal = 0;
    int i, err = 0;

    /* Note: in the arrays @rows and @cols the elements 0 to
       n-1 are the values for each node, and the elements n
       are the maxima.

       The tasks here are (a) to allocate a matrix of the right
       size to hold the result of the "reduce" operation and (b)
       to allocate workspace big enough to hold the values of
       the largest individual matrix.
    */

    if (op == GRETL_MPI_SUM || op == GRETL_MPI_PROD) {
	rtotal = rows[n];
	ctotal = cols[n];
    } else if (op == GRETL_MPI_HCAT) {
	/* there's a known common number of rows, but
	   we need to determine the total number of
	   columns in the result
	*/
	rtotal = rows[n];
	for (i=0; i<n; i++) {
	    if (rows[i] > 0) {
		ctotal += cols[i];
	    }
	}
    } else if (op == GRETL_MPI_VCAT) {
	/* there's a known common number of columns, but
	   we need to determine the total number of rows
	   in the result
	*/
	ctotal = cols[n];
	for (i=0; i<n; i++) {
	    if (cols[i] > 0) {
		rtotal += rows[i];
	    }
	}
    }

    if (rtotal == 0 || ctotal == 0) {
	err = E_DATA;
    } else if (op == GRETL_MPI_PROD) {
	*pm = gretl_unit_matrix_new(rtotal, ctotal);
    } else {
	*pm = gretl_zero_matrix_new(rtotal, ctotal);
    }

    if (*pm == NULL) {
	err = E_ALLOC;
    } else {
	size_t xsize = rows[n] * cols[n];

	*px = malloc(xsize * sizeof **px);
	if (*px == NULL) {
	    gretl_matrix_free(*pm);
	    *pm = NULL;
	    err = E_ALLOC;
	}
    }

    return err;
}

static int matrix_reduce_step (gretl_matrix *mtarg,
			       double * restrict src,
			       int n, Gretl_MPI_Op op,
			       int *offset)
{
    double * restrict targ = mtarg->val;
    int i;

    if (op == GRETL_MPI_SUM) {
	for (i=0; i<n; i++) {
	    targ[i] += src[i];
	}
    } else if (op == GRETL_MPI_PROD) {
	for (i=0; i<n; i++) {
	    targ[i] *= src[i];
	}
    } else if (op == GRETL_MPI_HCAT) {
	int k = *offset;

	for (i=0; i<n; i++) {
	    targ[k++] = src[i];
	}
	*offset = k;
    } else if (op == GRETL_MPI_VCAT) {
	int rmin = *offset;
	int nrows = n / mtarg->cols;
	int rmax = rmin + nrows;
	int j, k = 0;

	for (j=0; j<mtarg->cols; j++) {
	    for (i=rmin; i<rmax; i++) {
		gretl_matrix_set(mtarg, i, j, src[k++]);
	    }
	}
	*offset += nrows;
    }

    return 0;
}

static int invalid_rank_error (int r)
{
    gretl_errmsg_sprintf(_("Invalid MPI rank %d"), r);
    return E_DATA;
}

static int gretl_comm_check (int root, int *idp, int *npp)
{
    int np;

    mpi_comm_size(mpi_comm_world, &np);
    if (root < 0 || root >= np) {
	return invalid_rank_error(root);
    }

    if (npp != NULL) {
	*npp = np;
    }
    if (idp != NULL) {
	mpi_comm_rank(mpi_comm_world, idp);
    }

    return 0;
}

int gretl_matrix_mpi_reduce (gretl_matrix *sm,
			     gretl_matrix **pm,
			     Gretl_MPI_Op op,
			     int root,
			     gretlopt opt)
{
    gretl_matrix *rm = NULL;
    double *val = NULL;
    int *rows = NULL;
    int *cols = NULL;
    int rc[MI_LEN] = {0};
    int id, np;
    int err = 0;

    if (op != GRETL_MPI_SUM &&
	op != GRETL_MPI_PROD &&
	op != GRETL_MPI_HCAT &&
	op != GRETL_MPI_VCAT) {
	return E_DATA;
    }

    err = gretl_comm_check(root, &id, &np);
    if (err) {
	return err;
    }

    if (id != root) {
	/* send matrix dimensions to root */
	fill_matrix_info(rc, sm);
	err = mpi_send(rc, MI_LEN, mpi_int, root, TAG_MATRIX_INFO,
		       mpi_comm_world);
    } else {
	/* root: gather dimensions from other processes,
	   check for conformability and allocate storage
	*/
	int i;

	rows = malloc((np+1) * sizeof *rows);
	cols = malloc((np+1) * sizeof *cols);
	if (rows == NULL || cols == NULL) {
	    err = E_ALLOC;
	}

	/* initialize record of row/col maxima */
	rows[np] = 0;
	cols[np] = 0;

	for (i=0; i<np && !err; i++) {
	    if (i == root) {
		if (sm != NULL) {
		    rows[i] = sm->rows;
		    cols[i] = sm->cols;
		} else {
		    rows[i] = cols[i] = 0;
		}
	    } else {
		err = mpi_recv(rc, MI_LEN, mpi_int, i, TAG_MATRIX_INFO,
			       mpi_comm_world, MPI_STATUS_IGNORE);
		if (!err) {
		    rows[i] = rc[0];
		    cols[i] = rc[1];
		}
	    }
	}

	if (!err) {
	    err = matrix_dims_check(rows, cols, np, op);
	}

	if (!err) {
	    err = matrix_reduce_alloc(rows, cols, np, op,
				      &rm, &val);
	}
    }

    if (!err) {
	if (id != root) {
	    /* send data to root */
	    int sendsize = rc[0] * rc[1];

	    if (sendsize > 0) {
		err = mpi_send(sm->val, sendsize, mpi_double, root,
			       TAG_MATRIX_VAL, mpi_comm_world);
	    }
	} else {
	    /* root gathers and processes data */
	    int i, recvsize, offset = 0;

	    for (i=0; i<np && !err; i++) {
		recvsize = rows[i] * cols[i];
		if (recvsize == 0) {
		    continue;
		}
		if (i == root) {
		    if (sm != NULL) {
			err = matrix_reduce_step(rm, sm->val, recvsize, op,
						 &offset);
		    }
		} else {
		    err = mpi_recv(val, recvsize, mpi_double, i,
				   TAG_MATRIX_VAL,  mpi_comm_world,
				   MPI_STATUS_IGNORE);
		    if (!err) {
			err = matrix_reduce_step(rm, val, recvsize, op,
						 &offset);
		    }
		}
	    }
	}
    }

    if (id == root) {
	/* clean up */
	free(val);
	free(rows);
	free(cols);
	/* handle return value */
	if (!err) {
	    *pm = rm;
	} else {
	    gretl_matrix_free(rm);
	}
    }

    if (!err && (opt & OPT_A)) {
	err = gretl_matrix_bcast(pm, id, root);
    }

    return err;
}

int gretl_array_mpi_reduce (gretl_array *sa,
			    gretl_array **pa,
			    Gretl_MPI_Op op,
			    int root)
{
    gretl_array *a = NULL;
    GretlType atype, etype;
    int *nmvec = NULL;
    int id, np, nm;
    int i, j;
    int err = 0;

    if (op != GRETL_MPI_ACAT) {
	return E_DATA;
    }

    err = gretl_comm_check(root, &id, &np);
    if (err) {
	return err;
    }

    atype = gretl_array_get_type(sa);
    etype = gretl_array_get_content_type(sa);
    nm = gretl_array_get_length(sa);

    if (id != root) {
	/* send size of our array to root */
	err = mpi_send(&nm, 1, mpi_int, root, TAG_ARRAY_LEN,
		       mpi_comm_world);
    } else {
	/* root: gather and record array sizes from other processes;
	   allocate composite array
	*/
	int ntotal = 0;

	nmvec = malloc(np * sizeof *nmvec);
	if (nmvec == NULL) {
	    err = E_ALLOC;
	}
	for (i=0; i<np && !err; i++) {
	    if (i == root) {
		nmvec[i] = nm;
	    } else {
		err = mpi_recv(&nmvec[i], 1, mpi_int, i, TAG_ARRAY_LEN,
			       mpi_comm_world, MPI_STATUS_IGNORE);
	    }
	    ntotal += nmvec[i];
	}
	if (!err) {
	    a = gretl_array_new(atype, ntotal, &err);
	}
    }

    if (!err) {
	void *data;

	if (id != root) {
	    /* send our elements to root */
	    for (j=0; j<nm; j++) {
		data = gretl_array_get_data(sa, j);
		err = gretl_mpi_send(data, etype, root);
	    }
	} else {
	    /* root: gather elements from other processes
	       and pack into big array
	    */
	    int k = 0;

	    for (i=0; i<np && !err; i++) {
		for (j=0; j<nmvec[i] && !err; j++) {
		    if (i == root) {
			data = gretl_array_get_element(sa, j, NULL, &err);
			err = gretl_array_set_element(a, k++, data, etype, 1);
		    } else {
			data = mpi_receive_element(i, etype, &err);
			if (!err) {
			    err = gretl_array_set_data(a, k++, data);
			}
		    }
		}
	    }
	}
    }

    if (id == root) {
	free(nmvec);
	/* handle return value */
	if (!err) {
	    *pa = a;
	} else {
	    gretl_array_destroy(a);
	}
    }

    return err;
}

int gretl_scalar_mpi_reduce (double x,
			     double *xp,
			     Gretl_MPI_Op op,
			     int root,
			     gretlopt opt)
{
    MPI_Op mpi_op;
    int np, ret;

    mpi_comm_size(mpi_comm_world, &np);
    if (root < 0 || root >= np) {
	return invalid_rank_error(root);
    }

    if (op == GRETL_MPI_SUM) {
	mpi_op = mpi_sum;
    } else if (op == GRETL_MPI_PROD) {
	mpi_op = mpi_prod;
    } else if (op == GRETL_MPI_MAX) {
	mpi_op = mpi_max;
    } else if (op == GRETL_MPI_MIN) {
	mpi_op = mpi_min;
    } else {
	return E_DATA;
    }

    if (opt & OPT_A) {
	ret = mpi_allreduce(&x, xp, 1, mpi_double, mpi_op,
			    mpi_comm_world);
    } else {
	ret = mpi_reduce(&x, xp, 1, mpi_double, mpi_op,
			 root, mpi_comm_world);
    }

    return ret;
}

static void maybe_date_matrix (gretl_matrix *m, int *rc)
{
    if (rc[3] >= 0 && rc[4] >= rc[3]) {
	gretl_matrix_set_t1(m, rc[3]);
	gretl_matrix_set_t2(m, rc[4]);
    }
}

static int gretl_matrix_bcast (gretl_matrix **pm, int id, int root)
{
    gretl_matrix *m = NULL;
    int rc[MI_LEN];
    int err = 0;

    if (id == root) {
	m = *pm;
	fill_matrix_info(rc, m);
    }

    /* broadcast the matrix dimensions first */
    err = mpi_bcast(rc, MI_LEN, mpi_int, root, mpi_comm_world);

    if (!err && id != root) {
	int r = rc[0];
	int c = rc[1];
	int cmplx = rc[2];

	/* everyone but root needs to allocate space */
	if (cmplx) {
	    *pm = m = gretl_cmatrix_new(r, c);
	} else {
	    *pm = m = gretl_matrix_alloc(r, c);
	}
	if (m == NULL) {
	    return E_ALLOC;
	} else {
	    maybe_date_matrix(m, rc);
	}
    }

    if (!err) {
	/* broadcast the matrix content */
	int n = rc[0] * rc[1];

	if (rc[2]) {
	    n *= 2;
	}

	/* FIXME we can get a hang here with 100% CPU if
	   a worker bombs out on bcast(); in that case
	   it seems that root's call never returns --
	   or maybe not before some looong time-out.
	*/
	err = mpi_bcast(m->val, n, mpi_double, root,
			mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

static int gretl_series_bcast (double **px, int len,
			       int id, int root)
{
    double *x = NULL;
    int err = 0;

    if (id == root) {
	x = *px;
	if (x == NULL) {
	    return E_DATA;
	}
    } else {
	/* everyone else needs to allocate space */
	x = malloc(len * sizeof *x);
	if (x == NULL) {
	    return E_ALLOC;
	}
    }

    if (!err) {
	/* broadcast the series */
	err = mpi_bcast(x, len, mpi_double, root, mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
    } else if (id != root) {
	*px = x;
    }

    return err;
}

static int gretl_list_bcast (int **plist, int id, int root)
{
    int *list = NULL;
    int len = 0;
    int err = 0;

    if (id == root) {
	list = *plist;
	if (list == NULL) {
	    return E_DATA;
	} else {
	    len = list[0];
	}
    }

    /* broadcast the list length first */
    err = mpi_bcast(&len, 1, mpi_int, root, mpi_comm_world);

    if (!err && id != root) {
	/* everyone but root needs to allocate space */
	list = gretl_list_new(len);
	if (list == NULL) {
	    return E_ALLOC;
	}
    }

    if (!err) {
	/* broadcast the list */
	err = mpi_bcast(list, len + 1, mpi_int, root, mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
    } else if (id != root) {
	*plist = list;
    }

    return err;
}

static int gretl_string_bcast (char **pbuf, int id, int root)
{
    char *buf = NULL;
    int bytes = 0;
    int err = 0;

    if (id == root) {
	buf = *pbuf;
	if (buf == NULL) {
	    return E_DATA;
	} else {
	    bytes = strlen(buf) + 1;
	}
    }

    /* broadcast the buffer size first */
    err = mpi_bcast(&bytes, 1, mpi_int, root, mpi_comm_world);

    if (!err && id != root) {
	/* everyone but root needs to allocate space */
	buf = calloc(bytes, 1);
	if (buf == NULL) {
	    return E_ALLOC;
	}
    }

    if (!err) {
	/* broadcast the buffer */
	err = mpi_bcast(buf, bytes, mpi_byte, root, mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
    } else if (id != root) {
	*pbuf = buf;
    }

    return err;
}

static int gretl_scalar_bcast (double *px, int root)
{
    int err;

    err = mpi_bcast(px, 1, mpi_double, root, mpi_comm_world);

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

static int gretl_int_bcast (int *pi, int root)
{
    int err;

    err = mpi_bcast(pi, 1, mpi_int, root, mpi_comm_world);

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

static int gretl_unsigned_bcast (unsigned int *pu, int root)
{
    int err;

    err = mpi_bcast(pu, 1, mpi_unsigned, root, mpi_comm_world);

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

/* Compose a message indicating the type (and size, if
   applicable) of a bundle member to be passed via MPI.
*/

static int compose_msgbuf (char *buf, GretlType type,
			   const char *key, int *size,
			   void *data)
{
    int n = 0;

    if (type == GRETL_TYPE_DOUBLE) {
	n = sizeof(double);
    } else if (type == GRETL_TYPE_INT) {
	n = sizeof(int);
    } else if (type == GRETL_TYPE_UNSIGNED) {
	n = sizeof(unsigned int);
    } else if (type == GRETL_TYPE_SERIES) {
	n = *size;
    } else if (type == GRETL_TYPE_STRING) {
	char *s = data;
	n = strlen(s) + 1;
    } else if (type == GRETL_TYPE_LIST) {
	int *list = data;
	n = (list[0] + 1) * sizeof(int);
    } else if (type == GRETL_TYPE_MATRIX ||
	       type == GRETL_TYPE_BUNDLE ||
	       type == GRETL_TYPE_ARRAY) {
	n = 0;
    } else {
	return E_DATA;
    }

    sprintf(buf, "%d %s %d", type, key, n);
    *size = n;

    return 0;
}

/* Parse the information sent in @buf by root, regarding a
   bundle member.
*/

static int parse_msgbuf (const char *buf, char *key,
			 GretlType *type, int *size)
{
    int t;

    if (sscanf(buf, "%d %s %d", &t, key, size) != 3) {
	return E_DATA;
    } else {
	*type = t;
	return 0;
    }
}

static int gretl_bundle_bcast (gretl_bundle **pb,
			       int id, int root)
{
    gretl_bundle *b = NULL;
    gretl_matrix *m = NULL;
    gretl_array *a = NULL;
    gretl_bundle *bsub = NULL;
    gretl_array *keys = NULL;
    double *x = NULL;
    GretlType type;
    char key[32];
    void *data;
    char msgbuf[64];
    int msglen;
    int size;
    int i, nk;
    int err = 0;

    if (id == root) {
	b = *pb;
	nk = gretl_bundle_get_n_keys(b);
	keys = gretl_bundle_get_keys(b, &err);
	if (err) {
	    return err;
	}
    }

    /* broadcast the number of keys first */
    err = mpi_bcast(&nk, 1, mpi_int, root, mpi_comm_world);

    if (!err && id != root) {
	/* everyone but root needs to start a bundle */
	b = gretl_bundle_new();
	if (b == NULL) {
	    return E_ALLOC;
	}
    }

    msglen = sizeof msgbuf;

    for (i=0; i<nk && !err; i++) {
	/* loop across bundle keys */
	memset(msgbuf, 0, msglen);
	size = 0;
	data = NULL;
	m = NULL;
	a = NULL;
	bsub = NULL;
	x = NULL;

	if (id == root) {
	    const char *rkey = gretl_array_get_data(keys, i);

	    data = gretl_bundle_get_data(b, rkey, &type, &size, &err);
	    if (!err) {
		compose_msgbuf(msgbuf, type, rkey, &size, data);
		if (type == GRETL_TYPE_MATRIX) {
		    m = (gretl_matrix *) data;
		} else if (type == GRETL_TYPE_ARRAY) {
		    a = (gretl_array *) data;
		} else if (type == GRETL_TYPE_BUNDLE) {
		    bsub = (gretl_bundle *) data;
		} else if (type == GRETL_TYPE_SERIES) {
		    x = (double *) data;
		}
	    }
	}
	if (!err) {
	    /* broadcast info on the bundle member */
	    err = mpi_bcast(msgbuf, msglen, mpi_byte, root, mpi_comm_world);
	}
	if (!err) {
	    err = parse_msgbuf(msgbuf, key, &type, &size);
	}
	if (err) {
	    break;
	}
	if (type == GRETL_TYPE_MATRIX) {
	    err = gretl_matrix_bcast(&m, id, root);
	    if (!err && id != root) {
		err = gretl_bundle_donate_data(b, key, m, type, 0);
	    }
	} else if (type == GRETL_TYPE_ARRAY) {
	    err = gretl_array_bcast(&a, id, root);
	    if (!err && id != root) {
		err = gretl_bundle_donate_data(b, key, a, type, 0);
	    }
	} else if (type == GRETL_TYPE_BUNDLE) {
	    err = gretl_bundle_bcast(&bsub, id, root);
	    if (!err && id != root) {
		err = gretl_bundle_donate_data(b, key, bsub, type, 0);
	    }
	} else if (type == GRETL_TYPE_SERIES) {
	    err = gretl_series_bcast(&x, size, id, root);
	    if (!err && id != root) {
		err = gretl_bundle_donate_data(b, key, x, type, size);
	    }
	} else {
	    /* scalar, string or list */
	    if (id != root) {
		data = calloc(size, 1);
		if (data == NULL) {
		    err = E_ALLOC;
		}
	    }
	    err = mpi_bcast(data, size, mpi_byte, root, mpi_comm_world);
	    if (!err && id != root && data != NULL) {
		if (gretl_is_scalar_type(type)) {
		    err = gretl_bundle_set_data(b, key, data, type, 0);
		    free(data);
		} else {
		    err = gretl_bundle_donate_data(b, key, data, type, 0);
		}
	    }
	}
    }

    gretl_array_destroy(keys);

    mpi_barrier(mpi_comm_world);

    if (err) {
	gretl_mpi_error(&err);
    } else if (id != root) {
	*pb = b;
    }

    return err;
}

static int gretl_array_bcast (gretl_array **pa, int id, int root)
{
    gretl_array *a = NULL;
    GretlType type = 0;
    int nelem = 0;
    int st[2];
    int i, err = 0;

    if (id == root) {
	a = *pa;
	st[0] = nelem = gretl_array_get_length(a);
	st[1] = type = gretl_array_get_type(a);
    }

    /* broadcast the array size and type */
    err = mpi_bcast(st, 2, mpi_int, root, mpi_comm_world);

    if (!err && id != root) {
	/* everyone but root needs to allocate an array */
	nelem = st[0];
	type = st[1];
	*pa = a = gretl_array_new(type, nelem, &err);
    }

    for (i=0; i<nelem && !err; i++) {
	void *ptr, *data = NULL;

	if (id == root) {
	    data = gretl_array_get_data(a, i);
	}
	ptr = &data;
	if (type == GRETL_TYPE_MATRICES) {
	    err = gretl_matrix_bcast(ptr, id, root);
	    if (id != root) {
		data = *(gretl_matrix **) ptr;
	    }
	} else if (type == GRETL_TYPE_STRINGS) {
	    err = gretl_string_bcast(ptr, id, root);
	    if (id != root) {
		data = *(char **) ptr;
	    }
	} else if (type == GRETL_TYPE_BUNDLES) {
	    err = gretl_bundle_bcast(ptr, id, root);
	    if (id != root) {
		data = *(gretl_bundle **) ptr;
	    }
	} else if (type == GRETL_TYPE_LISTS) {
	    err = gretl_list_bcast(ptr, id, root);
	    if (id != root) {
		data = *(int **) ptr;
	    }
	} else if (type == GRETL_TYPE_ARRAYS) {
	    err = gretl_array_bcast(ptr, id, root);
	    if (id != root) {
		data = *(gretl_array **) ptr;
	    }
	} else {
	    /* ?? */
	    err = E_TYPES;
	}
	if (!err && id != root) {
	    gretl_array_set_data(a, i, data);
	}
    }

    mpi_barrier(mpi_comm_world);

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

int gretl_mpi_barrier (void)
{
    return mpi_barrier(mpi_comm_world) != MPI_SUCCESS;
}

/**
 * gretl_mpi_bcast:
 * @p: the location of the object to be broadcast.
 * @type: the type of the object to which @p points.
 * @root: the rank of the root process.
 *
 * Broadcasts the value referenced by @p to MPI_COMM_WORLD.
 * @type must be GRETL_TYPE_MATRIX, in which case
 * @p should be a (**gretl_matrix) pointer; GRETL_TYPE_BUNDLE,
 * in which case @p should be a (**gretl_bundle) pointer;
 * or GRETL_TYPE_DOUBLE, in which case @p should be a (*double)
 * pointer.
 *
 * Returns: 0 on successful completion, non-zero code otherwise.
 **/

int gretl_mpi_bcast (void *p, GretlType type, int root)
{
    int id, err;

    err = gretl_comm_check(root, &id, NULL);
    if (err) {
	return err;
    }

    if (type == GRETL_TYPE_DOUBLE) {
	return gretl_scalar_bcast((double *) p, root);
    } else if (type == GRETL_TYPE_INT) {
	return gretl_int_bcast((int *) p, root);
    } else if (type == GRETL_TYPE_UNSIGNED) {
	return gretl_unsigned_bcast((unsigned int *) p, root);
    } else if (type == GRETL_TYPE_MATRIX) {
	return gretl_matrix_bcast((gretl_matrix **) p, id, root);
    } else if (type == GRETL_TYPE_BUNDLE) {
	return gretl_bundle_bcast((gretl_bundle **) p, id, root);
    } else if (type == GRETL_TYPE_ARRAY) {
	return gretl_array_bcast((gretl_array **) p, id, root);
    } else if (type == GRETL_TYPE_STRING) {
	return gretl_string_bcast((char **) p, id, root);
    } else if (type == GRETL_TYPE_LIST) {
	return gretl_list_bcast((int **) p, id, root);
    } else {
	return E_DATA;
    }
}

/* this "send" function is public, for convenience of callers
   such as the svm plugin
*/

int gretl_matrix_mpi_send (const gretl_matrix *m, int dest)
{
    int rc[MI_LEN];
    int err;

    fill_matrix_info(rc, m);

    err = mpi_send(rc, MI_LEN, mpi_int, dest, TAG_MATRIX_INFO,
		   mpi_comm_world);

    if (!err) {
	int n = m->rows * m->cols;

	if (m->is_complex) {
	    n *= 2;
	}
	err = mpi_send(m->val, n, mpi_double, dest, TAG_MATRIX_VAL,
		       mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

static int gretl_string_send (char *s, int dest)
{
    int n = strlen(s) + 1;
    int err;

    err = mpi_send(&n, 1, mpi_int, dest, TAG_STR_LEN,
		   mpi_comm_world);

    if (!err) {
	err = mpi_send(s, n, mpi_byte, dest, TAG_STR_VAL,
		       mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

static int gretl_list_send (int *list, int dest)
{
    int n = list[0];
    int err;

    err = mpi_send(&n, 1, mpi_int, dest, TAG_LIST_LEN,
		   mpi_comm_world);

    if (!err) {
	err = mpi_send(list, n+1, mpi_int, dest, TAG_LIST_VAL,
		       mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

static int gretl_scalar_send (double *px, int dest)
{
    int err;

    err = mpi_send(px, 1, mpi_double, dest, TAG_SCALAR_VAL,
		   mpi_comm_world);

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

static int gretl_int_send (int *pi, int dest)
{
    int err;

    err = mpi_send(pi, 1, mpi_int, dest, TAG_INT_VAL,
		   mpi_comm_world);

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

static int gretl_array_send (gretl_array *a, int dest)
{
    GretlType type = 0;
    int nelem = 0;
    int st[2];
    int i, err = 0;

    st[0] = nelem = gretl_array_get_length(a);
    st[1] = type = gretl_array_get_type(a);

    /* send the array size and type */
    err = mpi_send(st, 2, mpi_int, dest, TAG_ARRAY_INFO,
		   mpi_comm_world);

    for (i=0; i<nelem && !err; i++) {
	void *data = gretl_array_get_data(a, i);

	if (type == GRETL_TYPE_MATRICES) {
	    err = gretl_matrix_mpi_send(data, dest);
	} else if (type == GRETL_TYPE_STRINGS) {
	    err = gretl_string_send(data, dest);
	} else if (type == GRETL_TYPE_BUNDLES) {
	    err = gretl_bundle_send(data, dest);
	} else if (type == GRETL_TYPE_LISTS) {
	    err = gretl_list_send(data, dest);
	} else if (type == GRETL_TYPE_ARRAYS) {
	    err = gretl_array_send(data, dest);
	}
    }

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

/**
 * gretl_mpi_send:
 * @p: pointer to the object to be sent.
 * @type: the type of the object.
 * @dest: the MPI rank of the destination.
 *
 * Sends the value referenced by @p, of gretl type @type, to the
 * MPI process with rank @dest.
 *
 * Returns: 0 on successful completion, non-zero code otherwise.
 **/

int gretl_mpi_send (void *p, GretlType type, int dest)
{
    int np;

    mpi_comm_size(mpi_comm_world, &np);
    if (dest < 0 || dest >= np) {
	return invalid_rank_error(dest);
    }

    if (type == GRETL_TYPE_DOUBLE) {
	return gretl_scalar_send((double *) p, dest);
    } else if (type == GRETL_TYPE_INT) {
	return gretl_int_send((int *) p, dest);
    } else if (type == GRETL_TYPE_MATRIX) {
	return gretl_matrix_mpi_send((gretl_matrix *) p, dest);
    } else if (type == GRETL_TYPE_BUNDLE) {
	return gretl_bundle_send((gretl_bundle *) p, dest);
    } else if (type == GRETL_TYPE_ARRAY) {
	return gretl_array_send((gretl_array *) p, dest);
    } else if (type == GRETL_TYPE_STRING) {
	return gretl_string_send((char *) p, dest);
    } else if (type == GRETL_TYPE_LIST) {
	return gretl_list_send((int *) p, dest);
    } else {
	return E_TYPES;
    }
}

gretl_matrix *gretl_matrix_mpi_receive (int source,
					int *err)
{
    gretl_matrix *m = NULL;
    int rc[MI_LEN];

    *err = mpi_recv(rc, MI_LEN, mpi_int, source, TAG_MATRIX_INFO,
		    mpi_comm_world, MPI_STATUS_IGNORE);

    if (!*err) {
	int r = rc[0];
	int c = rc[1];
	int cmplx = rc[2];
	int n = r * c;

	if (cmplx) {
	    m = gretl_cmatrix_new(r, c);
	    n *= 2;
	} else {
	    m = gretl_matrix_alloc(r, c);
	}

	if (m == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	} else {
	    *err = mpi_recv(m->val, n, mpi_double, source,
			    TAG_MATRIX_VAL, mpi_comm_world,
			    MPI_STATUS_IGNORE);
	    if (*err) {
		maybe_date_matrix(m, rc);
	    }
	}
    }

    if (*err) {
	gretl_mpi_error(err);
    }

    return m;
}

int gretl_matrix_mpi_fill (gretl_matrix **pm, int source)
{
    int rc[MI_LEN];
    int err;

    if (pm == NULL) {
	return E_DATA;
    }

    err = mpi_recv(rc, MI_LEN, mpi_int, source, TAG_MATRIX_INFO,
		   mpi_comm_world, MPI_STATUS_IGNORE);

    if (!err) {
	gretl_matrix *m = *pm;
	int r = rc[0];
	int c = rc[1];
	int cmplx = rc[2];
	int n = r * c;

	if (m == NULL) {
	    if (cmplx) {
		m = gretl_cmatrix_new(r, c);
		n *= 2;
	    } else {
		m = gretl_matrix_alloc(r, c);
	    }
	    if (m == NULL) {
		err = E_ALLOC;
	    } else {
		*pm = m;
	    }
	} else if (m->rows != r || m->cols != c ||
		   m->is_complex != cmplx) {
	    err = E_NONCONF;
	}

	if (!err) {
	    err = mpi_recv(m->val, n, mpi_double, source,
			   TAG_MATRIX_VAL, mpi_comm_world,
			   MPI_STATUS_IGNORE);
	    if (err) {
		maybe_date_matrix(m, rc);
	    }
	}
    }

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

static int *gretl_list_receive (int source, int *err)
{
    int *list = NULL;
    int n;

    *err = mpi_recv(&n, 1, mpi_int, source, TAG_LIST_LEN,
		    mpi_comm_world, MPI_STATUS_IGNORE);

    if (!*err) {
	list = gretl_list_new(n);
	if (list == NULL) {
	    *err = E_ALLOC;
	} else {
	    *err = mpi_recv(list, n+1, mpi_int, source,
			    TAG_LIST_VAL, mpi_comm_world,
			    MPI_STATUS_IGNORE);
	}
    }

    return list;
}

static char *gretl_string_receive (int source, int *err)
{
    char *s = NULL;
    int n;

    *err = mpi_recv(&n, 1, mpi_int, source, TAG_STR_LEN,
		    mpi_comm_world, MPI_STATUS_IGNORE);

    if (!*err) {
	s = calloc(n, 1);
	if (s == NULL) {
	    *err = E_ALLOC;
	} else {
	    *err = mpi_recv(s, n, mpi_byte, source,
			    TAG_STR_VAL, mpi_comm_world,
			    MPI_STATUS_IGNORE);
	}
    }

    return s;
}

static double gretl_scalar_receive (int source, int *err)
{
    double x = NADBL;

    *err = mpi_recv(&x, 1, mpi_double, source, TAG_SCALAR_VAL,
		    mpi_comm_world, MPI_STATUS_IGNORE);
    if (*err) {
	gretl_mpi_error(err);
    }

    return x;
}

static double gretl_int_receive (int source, int *err)
{
    int i = 0;

    *err = mpi_recv(&i, 1, mpi_int, source, TAG_INT_VAL,
		    mpi_comm_world, MPI_STATUS_IGNORE);
    if (*err) {
	gretl_mpi_error(err);
    }

    return i;
}

static gretl_array *gretl_array_receive (int source, int *err)
{
    gretl_array *a = NULL;
    GretlType type = 0;
    int nelem = 0;
    int st[2];
    int i;

    /* get the array size and type */
    *err = mpi_recv(st, 2, mpi_int, source, TAG_ARRAY_INFO,
		    mpi_comm_world, MPI_STATUS_IGNORE);

    if (!*err) {
	nelem = st[0];
	type = st[1];
	a = gretl_array_new(type, nelem, err);
    }

    for (i=0; i<nelem && !*err; i++) {
	void *data = NULL;

	if (type == GRETL_TYPE_MATRICES) {
	    data = gretl_matrix_mpi_receive(source, err);
	} else if (type == GRETL_TYPE_STRINGS) {
	    data = gretl_string_receive(source, err);
	} else if (type == GRETL_TYPE_BUNDLES) {
	    data = gretl_bundle_receive(source, err);
	} else if (type == GRETL_TYPE_LISTS) {
	    data = gretl_list_receive(source, err);
	} else if (type == GRETL_TYPE_ARRAYS) {
	    data = gretl_array_receive(source, err);
	}
	if (data != NULL) {
	    gretl_array_set_data(a, i, data);
	}
    }

    return a;
}

static GretlType type_from_status (MPI_Status *status)
{
    if (status->MPI_TAG == TAG_MATRIX_INFO) {
	return GRETL_TYPE_MATRIX;
    } else if (status->MPI_TAG == TAG_BUNDLE_SIZE) {
	return GRETL_TYPE_BUNDLE;
    } else if (status->MPI_TAG == TAG_ARRAY_INFO) {
	return GRETL_TYPE_ARRAY;
    } else if (status->MPI_TAG == TAG_SCALAR_VAL) {
	return GRETL_TYPE_DOUBLE;
    } else if (status->MPI_TAG == TAG_INT_VAL) {
	return GRETL_TYPE_INT;
    } else if (status->MPI_TAG == TAG_STR_LEN) {
	return GRETL_TYPE_STRING;
    } else if (status->MPI_TAG == TAG_LIST_LEN) {
	return GRETL_TYPE_LIST;
    } else {
	return GRETL_TYPE_NONE;
    }
}

void *gretl_mpi_receive (int source, GretlType *ptype, int *err)
{
    static double x;
    static int k;
    void *ret = NULL;
    MPI_Status status;
    int np;

    mpi_comm_size(mpi_comm_world, &np);
    if (source < 0 || source >= np) {
	*err = invalid_rank_error(source);
	return NULL;
    }

    /* check for the type of thing of offer from @source */
    mpi_probe(source, MPI_ANY_TAG, mpi_comm_world, &status);
    *ptype = type_from_status(&status);

    if (*ptype == GRETL_TYPE_DOUBLE) {
	x = gretl_scalar_receive(source, err);
	ret = &x;
    } else if (*ptype == GRETL_TYPE_INT) {
	k = gretl_int_receive(source, err);
	ret = &k;
    } else if (*ptype == GRETL_TYPE_MATRIX) {
	ret = gretl_matrix_mpi_receive(source, err);
    } else if (*ptype == GRETL_TYPE_BUNDLE) {
	ret = gretl_bundle_receive(source, err);
    } else if (*ptype == GRETL_TYPE_ARRAY) {
	ret = gretl_array_receive(source, err);
    } else if (*ptype == GRETL_TYPE_STRING) {
	ret = gretl_string_receive(source, err);
    } else if (*ptype == GRETL_TYPE_LIST) {
	ret = gretl_list_receive(source, err);
    } else {
	*err = E_TYPES;
    }

    return ret;
}

static void *mpi_receive_element (int source, GretlType etype,
				  int *err)
{
    void *ret = NULL;
    MPI_Status status;
    GretlType srctype;

    mpi_probe(source, MPI_ANY_TAG, mpi_comm_world, &status);
    srctype = type_from_status(&status);

    if (srctype != etype) {
	*err = E_TYPES;
    } else if (etype == GRETL_TYPE_MATRIX) {
	ret = gretl_matrix_mpi_receive(source, err);
    } else if (etype == GRETL_TYPE_BUNDLE) {
	ret = gretl_bundle_receive(source, err);
    } else if (etype == GRETL_TYPE_ARRAY) {
	ret = gretl_array_receive(source, err);
    } else if (etype == GRETL_TYPE_STRING) {
	ret = gretl_string_receive(source, err);
    } else if (etype == GRETL_TYPE_LIST) {
	ret = gretl_list_receive(source, err);
    }

    return ret;
}

/* "new"-style send/receive for bundles: we do everything
   in memory */

static int gretl_series_send (double *x, int n, int dest)
{
    int err;

    err = mpi_send(&n, 1, mpi_int, dest, TAG_SERIES_LEN,
		   mpi_comm_world);

    if (!err) {
	err = mpi_send(x, n, mpi_double, dest, TAG_SERIES_VAL,
		       mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

static double *gretl_series_receive (int source, int *err)
{
    double *x = NULL;
    int n;

    *err = mpi_recv(&n, 1, mpi_int, source, TAG_SERIES_LEN,
		    mpi_comm_world, MPI_STATUS_IGNORE);

    if (!*err) {
	x = malloc(n * sizeof *x);
	if (x == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	} else {
	    *err = mpi_recv(x, n, mpi_double, source,
			    TAG_SERIES_VAL, mpi_comm_world,
			    MPI_STATUS_IGNORE);
	}
    }

    if (*err) {
	gretl_mpi_error(err);
    }

    return x;
}

static int gretl_bundle_send (gretl_bundle *b, int dest)
{
    gretl_array *keys = NULL;
    GretlType type;
    void *data;
    char msgbuf[64];
    int msglen = sizeof msgbuf;
    int size;
    int i, nk;
    int err = 0;

    nk = gretl_bundle_get_n_keys(b);
    if (nk > 0) {
	keys = gretl_bundle_get_keys(b, &err);
	if (err) {
	    return err;
	}
    }

    /* send the number of keys first */
    err = mpi_send(&nk, 1, mpi_int, dest, TAG_BUNDLE_SIZE,
		   mpi_comm_world);

    if (nk == 0) {
	/* the bundle is empty */
	return 0;
    }

    for (i=0; i<nk && !err; i++) {
	/* loop across bundle keys */
	const char *key = gretl_array_get_data(keys, i);

	data = gretl_bundle_get_data(b, key, &type, &size, &err);
	if (!err) {
	    /* send info on the bundle member */
	    compose_msgbuf(msgbuf, type, key, &size, data);
	    err = mpi_send(msgbuf, msglen, mpi_byte, dest,
			   TAG_BMEMB_INFO, mpi_comm_world);
	}
	if (err) {
	    break;
	}
	if (type == GRETL_TYPE_MATRIX) {
	    err = gretl_matrix_mpi_send(data, dest);
	} else if (type == GRETL_TYPE_ARRAY) {
	    err = gretl_array_send(data, dest);
	} else if (type == GRETL_TYPE_BUNDLE) {
	    err = gretl_bundle_send(data, dest);
	} else if (type == GRETL_TYPE_SERIES) {
	    err = gretl_series_send(data, size, dest);
	} else if (type == GRETL_TYPE_STRING) {
	    err = gretl_string_send(data, dest);
	} else if (type == GRETL_TYPE_LIST) {
	    err = gretl_list_send(data, dest);
	} else if (type == GRETL_TYPE_DOUBLE) {
	    err = gretl_scalar_send(data, dest);
	} else if (type == GRETL_TYPE_INT) {
	    err = gretl_int_send(data, dest);
	}
    }

    gretl_array_destroy(keys);

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

gretl_bundle *gretl_bundle_receive (int source, int *err)
{
    gretl_bundle *b = NULL;
    char key[32];
    char msgbuf[64];
    int msglen = sizeof msgbuf;
    int size;
    int i, nk;

    /* get the number of keys first */
    *err = mpi_recv(&nk, 1, mpi_int, source, TAG_BUNDLE_SIZE,
		    mpi_comm_world, MPI_STATUS_IGNORE);

    if (!*err) {
	b = gretl_bundle_new();
	if (b == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	}
    }

    if (nk == 0) {
	/* the bundle is empty */
	return b;
    }

    for (i=0; i<nk && !*err; i++) {
	/* loop across bundle keys */
	GretlType type = 0;
	void *data = NULL;
	double x;

	memset(msgbuf, 0, msglen);
	size = 0;

	/* get info on bundle member @i */
	*err = mpi_recv(msgbuf, msglen, mpi_byte, source, TAG_BMEMB_INFO,
			mpi_comm_world, MPI_STATUS_IGNORE);
	if (!*err) {
	    *err = parse_msgbuf(msgbuf, key, &type, &size);
	}
	if (*err) {
	    break;
	}
	if (type == GRETL_TYPE_DOUBLE) {
	    x = gretl_scalar_receive(source, err);
	    data = &x;
	} else if (type == GRETL_TYPE_INT) {
	    x = gretl_int_receive(source, err);
	    data = &x;
	} else if (type == GRETL_TYPE_MATRIX) {
	    data = gretl_matrix_mpi_receive(source, err);
	} else if (type == GRETL_TYPE_ARRAY) {
	    data = gretl_array_receive(source, err);
	} else if (type == GRETL_TYPE_BUNDLE) {
	    data = gretl_bundle_receive(source, err);
	} else if (type == GRETL_TYPE_SERIES) {
	    data = gretl_series_receive(source, err);
	} else if (type == GRETL_TYPE_LIST) {
	    data = gretl_list_receive(source, err);
	} else if (type == GRETL_TYPE_STRING) {
	    data = gretl_string_receive(source, err);
	}
	if (!*err) {
	    *err = gretl_bundle_donate_data(b, key, data, type, size);
	}
    }

    if (*err) {
	gretl_mpi_error(err);
    }

    return b;
}

static void fill_tmp (double * restrict tmp,
		      const gretl_matrix *m,
		      int nr, int cmplx,
		      int *offset)
{
    double x;
    double complex z;
    int imin = *offset;
    int imax = imin + nr;
    int i, j, k = 0;

    for (j=0; j<m->cols; j++) {
	for (i=imin; i<imax; i++) {
	    if (cmplx) {
		z = gretl_cmatrix_get(m, i, j);
		tmp[k++] = creal(z);
		tmp[k++] = cimag(z);
	    } else {
		x = gretl_matrix_get(m, i, j);
		tmp[k++] = x;
	    }
	}
    }

    *offset += nr;
}

static int scatter_to_self (int *rc, double *val,
			    gretl_matrix **pm)
{
    int cmplx = rc[2];
    int err = 0;

    if (cmplx) {
	*pm = gretl_cmatrix_new(rc[0], rc[1]);
    } else {
	*pm = gretl_matrix_alloc(rc[0], rc[1]);
    }

    if (*pm == NULL) {
	err = E_ALLOC;
    } else {
	size_t n = rc[0] * rc[1] * sizeof *val;

	if (cmplx) {
	    n *= 2;
	}
	memcpy((*pm)->val, val, n);
    }

    return err;
}

static void matsplit_rule (int n, int np, int *k1, int *n1,
			   int *k2, int *n2)
{
    if (n <= np) {
	*k1 = n;
	*n1 = 1;
	*k2 = np - n;
	*n2 = 0;
    } else if (n % np == 0) {
	*k1 = np;
	*n1 = n / np;
	*k2 = 0;
	*n2 = 0;
    } else {
	int a = ceil(1.0e-12 + n / (double) np);
	int head = n % np;

        *k1 = head;
        *n1 = a;
        *k2 = np - head;
        *n2 = a - 1;
    }
}

int gretl_matrix_mpi_scatter (const gretl_matrix *m,
                              gretl_matrix **recvm,
                              Gretl_MPI_Op op,
                              int root)
{
    double *tmp = NULL;
    int id, np;
    int rc[MI_LEN] = {0};
    int err = 0;

    mpi_comm_rank(mpi_comm_world, &id);
    mpi_comm_size(mpi_comm_world, &np);

    if (root < 0 || root >= np) {
        return invalid_rank_error(root);
    }

    if (id == root) {
        int cmplx = m->is_complex;
        int k1, n1, k2, n2;
        int i, n;

        if (op == GRETL_MPI_VSPLIT) {
            /* scatter by rows */
            int offset = 0;

            matsplit_rule(m->rows, np, &k1, &n1, &k2, &n2);

            n = n1 * m->cols;
            rc[0] = n1;
            rc[1] = m->cols;
            if (cmplx) {
                rc[2] = 1;
                n *= 2;
            }

            /* we'll need a working buffer */
            tmp = malloc(n * sizeof *tmp);
            if (tmp == NULL) {
                err = E_ALLOC;
            }

            for (i=0; i<np; i++) {
                if (i == k1) {
                    rc[0] = n2;
                    n = n2 * m->cols;
                    if (cmplx) n *= 2;
                }
                if (rc[0] > 0) {
                    fill_tmp(tmp, m, rc[0], cmplx, &offset);
                }
                if (i == root) {
                    err = scatter_to_self(rc, tmp, recvm);
                } else {
                    err = mpi_send(rc, MI_LEN, mpi_int, i, TAG_MATRIX_INFO,
                                   mpi_comm_world);
                    err = mpi_send(tmp, n, mpi_double, i, TAG_MATRIX_VAL,
                                   mpi_comm_world);
                }
            }
        } else {
            /* scatter by columns */
            double *val = m->val;

            matsplit_rule(m->cols, np, &k1, &n1, &k2, &n2);

            n = n1 * m->rows;
            if (cmplx) {
                n *= 2;
            }
            fill_matrix_info(rc, m);
            rc[1] = n1;

            for (i=0; i<np; i++) {
                if (i == k1) {
                    rc[1] = n2;
                    n = n2 * m->rows;
                    if (cmplx) n *= 2;
                }
                if (i == root) {
                    err = scatter_to_self(rc, val, recvm);
                } else {
                    err = mpi_send(rc, MI_LEN, mpi_int, i, TAG_MATRIX_INFO,
                                   mpi_comm_world);
                    err = mpi_send(val, n, mpi_double, i, TAG_MATRIX_VAL,
                                   mpi_comm_world);
                }
                /* advance the read pointer */
                val += n;
            }
        }
    } else {
        /* non-root processes get their share-out */
        *recvm = gretl_matrix_mpi_receive(root, &err);
    }

    if (id == root) {
        free(tmp);
    }

    return err;
}

/* MPI timer */

double gretl_mpi_time (void)
{
    return mpi_wtime();
}

/* end MPI timer */

int gretl_mpi_rank (void)
{
    int id = -1;

    if (gretl_mpi_initialized()) {
	mpi_comm_rank(mpi_comm_world, &id);
    }

    return id;
}

int gretl_mpi_n_processes (void)
{
    int np = 0;

    if (gretl_mpi_initialized()) {
	mpi_comm_size(mpi_comm_world, &np);
    }

    return np;
}

int gretl_mpi_initialized (void)
{
    static int ret = -1;

    if (!gretl_MPI_initted) {
	return 0;
    }

    if (ret < 0) {
	mpi_initialized(&ret);
	ret = ret > 0;
    }

    return ret;
}
