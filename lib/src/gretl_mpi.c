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
# define mpi_sum        MPI_SUM
# define mpi_prod       MPI_PROD
# define mpi_max        MPI_MAX
# define mpi_min        MPI_MIN
#endif

enum {
    TAG_MATRIX_DIM = 1,
    TAG_MATRIX_VAL,
    TAG_SCALAR_VAL
};

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
static int (*mpi_probe) (int, int, MPI_Comm, MPI_Status *);
static double (*mpi_wtime) (void);
static int (*mpi_initialized) (int *);

static void *mpiget (void *handle, const char *name, int *err)
{
#ifdef WIN32
    void *p = GetProcAddress(handle, name);
#else
    void *p = dlsym(handle, name);
#endif
    
    if (p == NULL) {
	printf("mpi_dlget: couldn't find '%s'\n", name);
	*err += 1;
    }

#if MPI_DEBUG
    else {
	printf("mpi_dlget: '%s' -> %p\n", name, p);
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
    mpi_wtime        = mpiget(MPIhandle, "MPI_Wtime", &err);
    mpi_initialized  = mpiget(MPIhandle, "MPI_Initialized", &err);

#ifdef OMPI_MAJOR_VERSION
    if (!err) {
	mpi_comm_world = (MPI_Comm) mpiget(MPIhandle, "ompi_mpi_comm_world", &err);
	mpi_double     = (MPI_Datatype) mpiget(MPIhandle, "ompi_mpi_double", &err);
	mpi_int        = (MPI_Datatype) mpiget(MPIhandle, "ompi_mpi_int", &err);
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
    char msg1[BUFSIZ], msg2[BUFSIZ];
    int len, errclass;
    int id = 0;

    mpi_comm_rank(mpi_comm_world, &id);
    mpi_error_class(*err, &errclass);
    mpi_error_string(errclass, msg1, &len);
    mpi_error_string(*err, msg2, &len);
    gretl_errmsg_sprintf("%3d: %s %s\n", id, msg1, msg2);

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

static int matrix_reduce_step (gretl_matrix *targ, double *src,
			       int n, Gretl_MPI_Op op,
			       int *offset)
{
    int i;

    if (op == GRETL_MPI_SUM) {
	for (i=0; i<n; i++) {
	    targ->val[i] += src[i];
	}
    } else if (op == GRETL_MPI_PROD) {
	for (i=0; i<n; i++) {
	    targ->val[i] *= src[i];
	}
    } else if (op == GRETL_MPI_HCAT) {
	int k = *offset;

	for (i=0; i<n; i++) {
	    targ->val[k++] = src[i];
	}
	*offset = k;
    } else if (op == GRETL_MPI_VCAT) {
	int rmin = *offset;
	int nrows = n / targ->cols;
	int rmax = rmin + nrows;
	int j, k = 0;

	for (j=0; j<targ->cols; j++) {
	    for (i=rmin; i<rmax; i++) {
		gretl_matrix_set(targ, i, j, src[k++]);
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
    int rc[2] = {0};
    int id, np;
    int err = 0;

    if (op != GRETL_MPI_SUM &&
	op != GRETL_MPI_PROD &&
	op != GRETL_MPI_HCAT &&
	op != GRETL_MPI_VCAT) {
	return E_DATA;
    }

    mpi_comm_rank(mpi_comm_world, &id);
    mpi_comm_size(mpi_comm_world, &np);

    if (root < 0 || root >= np) {
	return invalid_rank_error(root);
    }

    if (id != root) {
	/* send matrix dimensions to root */
	rc[0] = sm->rows;
	rc[1] = sm->cols;
	err = mpi_send(rc, 2, mpi_int, root, TAG_MATRIX_DIM, 
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
		rows[i] = sm->rows;
		cols[i] = sm->cols;
	    } else {
		err = mpi_recv(rc, 2, mpi_int, i, TAG_MATRIX_DIM, 
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
		    err = matrix_reduce_step(rm, sm->val, recvsize, op,
					     &offset);
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
	err = gretl_matrix_mpi_bcast(pm, root);
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

    /* convert NA to NaN for use with MPI's built-in
       reduction functions */
    x = na(x) ? M_NA : x;

    if (opt & OPT_A) {
	ret = mpi_allreduce(&x, xp, 1, mpi_double, mpi_op, 
			    mpi_comm_world);
    } else {
	ret = mpi_reduce(&x, xp, 1, mpi_double, mpi_op, 
			 root, mpi_comm_world);
    }

    return ret;
}

int gretl_matrix_mpi_bcast (gretl_matrix **pm, int root)
{
    gretl_matrix *m = NULL;
    int rc[2];
    int id, np, err;

    mpi_comm_rank(mpi_comm_world, &id);
    mpi_comm_size(mpi_comm_world, &np);

    if (root < 0 || root >= np) {
	return invalid_rank_error(root);
    }

    if (id == root) {
	m = *pm;
	rc[0] = m->rows;
	rc[1] = m->cols;
    }

    /* broadcast the matrix dimensions first */
    err = mpi_bcast(rc, 2, mpi_int, root, mpi_comm_world);

    if (!err && id != root) {
	/* everyone but root needs to allocate space */
	*pm = m = gretl_matrix_alloc(rc[0], rc[1]);
	if (m == NULL) {
	    return E_ALLOC;
	}
    }

    if (!err) {
	/* broadcast the matrix content */
	int n = rc[0] * rc[1];
	
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

int gretl_scalar_mpi_bcast (double *px, int root)
{
    int np, err;

    mpi_comm_size(mpi_comm_world, &np);

    if (root < 0 || root >= np) {
	return invalid_rank_error(root);
    }

    err = mpi_bcast(px, 1, mpi_double, root, mpi_comm_world);

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

/**
 * gretl_mpi_bcast:
 * @p: the location of the object to be broadcast.
 * @type: the type of the object to which @p points.
 * @root: the rank of the root process.
 *
 * Broadcasts the value referenced by @p to MPI_COMM_WORLD.
 * At present @type must be GRETL_TYPE_MATRIX, in which case
 * @p should be a (**gretl_matrix) pointer, or GRETL_TYPE_DOUBLE,
 * in which case @p should be a (*double) pointer.
 *
 * Returns: 0 on successful completion, non-zero code otherwise.
 **/

int gretl_mpi_bcast (void *p, GretlType type, int root)
{
    if (type == GRETL_TYPE_DOUBLE) {
	return gretl_scalar_mpi_bcast((double *) p, root);
    } else if (type == GRETL_TYPE_MATRIX) {
	return gretl_matrix_mpi_bcast((gretl_matrix **) p, root);
    } else {
	return E_DATA;
    }
}

int gretl_matrix_mpi_send (const gretl_matrix *m, int dest)
{
    int rc[2] = {m->rows, m->cols};
    int err;

    err = mpi_send(rc, 2, mpi_int, dest, TAG_MATRIX_DIM, 
		   mpi_comm_world);

    if (!err) {
	int n = m->rows * m->cols;

	err = mpi_send(m->val, n, mpi_double, dest, TAG_MATRIX_VAL,
		       mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
    }    

    return err;
}

static int gretl_scalar_mpi_send (double *px, int dest)
{
    int err;

    err = mpi_send(px, 1, mpi_double, dest, TAG_SCALAR_VAL, 
		   mpi_comm_world);

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
 * Sends the value referenced by @p to the MPI process with
 * rank @dest. At present @type must be GRETL_TYPE_MATRIX, 
 * in which case  @p should be a (*gretl_matrix) pointer, or
 * GRETL_TYPE_DOUBLE, in which case @p should be a (*double)
 * pointer.
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
	return gretl_scalar_mpi_send((double *) p, dest);
    } else if (type == GRETL_TYPE_MATRIX) {
	return gretl_matrix_mpi_send((gretl_matrix *) p, dest);
    } else {
	return E_DATA;
    }
}

gretl_matrix *gretl_matrix_mpi_receive (int source, 
					int *err)
{
    gretl_matrix *m = NULL;
    int rc[2];

    *err = mpi_recv(rc, 2, mpi_int, source, TAG_MATRIX_DIM, 
		    mpi_comm_world, MPI_STATUS_IGNORE);

    if (!*err) {
	m = gretl_matrix_alloc(rc[0], rc[1]);
	if (m == NULL) {
	    *err = E_ALLOC;
	} else {
	    int n = rc[0] * rc[1];

	    *err = mpi_recv(m->val, n, mpi_double, source,
			    TAG_MATRIX_VAL, mpi_comm_world,
			    MPI_STATUS_IGNORE);
	    if (*err) {
		gretl_mpi_error(err);
	    }
	}
    }

    return m;
}

double gretl_scalar_mpi_receive (int source, int *err)
{
    double x = NADBL;

    *err = mpi_recv(&x, 1, mpi_double, source, TAG_SCALAR_VAL,
		    mpi_comm_world, MPI_STATUS_IGNORE);
    if (*err) {
	gretl_mpi_error(err);
    }

    return x;
}

int gretl_mpi_receive (int source, GretlType *type, 
		       gretl_matrix **pm,
		       double *px)
{
    MPI_Status status;
    int np, err = 0;

    mpi_comm_size(mpi_comm_world, &np);
    if (source < 0 || source >= np) {
	return invalid_rank_error(source);
    }

    /* check for the type of thing of offer from @source */
    mpi_probe(source, MPI_ANY_TAG, mpi_comm_world, &status);

    if (status.MPI_TAG == TAG_SCALAR_VAL) {
	*px = gretl_scalar_mpi_receive(source, &err);
	*type = GRETL_TYPE_DOUBLE;
    } else if (status.MPI_TAG == TAG_MATRIX_DIM) {
	*pm = gretl_matrix_mpi_receive(source, &err);
	*type = GRETL_TYPE_MATRIX;
    } else {
	err = E_DATA;
    }

    return err;
}

static void fill_tmp (const gretl_matrix *m, double *tmp,
		      int nr, int *offset)
{
    double x;
    int imin = *offset;
    int imax = imin + nr;
    int i, j, k = 0;

    for (j=0; j<m->cols; j++) {
	for (i=imin; i<imax; i++) {
	    x = gretl_matrix_get(m, i, j);
	    tmp[k++] = x;
	}
    }

    *offset += nr;
}

static int scatter_to_self (int *rc, double *val, 
			    gretl_matrix **pm)
{
    int err = 0;

    *pm = gretl_matrix_alloc(rc[0], rc[1]);

    if (*pm == NULL) {
	err = E_ALLOC;
    } else {
	size_t n = rc[0] * rc[1] * sizeof *val;

	memcpy((*pm)->val, val, n);
    }

    return err;
}

int gretl_matrix_mpi_scatter (const gretl_matrix *m,
			      gretl_matrix **recvm,
			      Gretl_MPI_Op op,
			      int root)
{
    double *tmp = NULL;
    int id, np;
    int rc[2];
    int err = 0;

    mpi_comm_rank(mpi_comm_world, &id);
    mpi_comm_size(mpi_comm_world, &np);

    if (root < 0 || root >= np) {
	return invalid_rank_error(root);
    }

    if (id == root) {
	int i, n;
	
	if (op == GRETL_MPI_VSPLIT) {
	    /* scatter by rows */
	    int nr = m->rows / np;
	    int rem = m->rows % np;
	    int offset = 0;

	    /* we'll need a working buffer */
	    tmp = malloc(m->cols * (nr + rem) * sizeof *tmp);
	    if (tmp == NULL) {
		err = E_ALLOC;
	    }

	    n = nr * m->cols;
	    rc[0] = nr;
	    rc[1] = m->cols;
	    
	    for (i=0; i<np; i++) {
		if (i == np - 1 && rem > 0) {
		    rc[0] += rem;
		    n += m->cols * rem;
		}
		fill_tmp(m, tmp, rc[0], &offset);
		if (i == root) {
		    err = scatter_to_self(rc, tmp, recvm);
		} else {
		    err = mpi_send(rc, 2, mpi_int, i, TAG_MATRIX_DIM, 
				   mpi_comm_world);
		    err = mpi_send(tmp, n, mpi_double, i, TAG_MATRIX_VAL, 
				   mpi_comm_world);
		}
	    }
	} else {
	    /* scatter by columns */
	    int nc = m->cols / np;
	    int rem = m->cols % np;
	    double *val = m->val;

	    n = m->rows * nc;
	    rc[0] = m->rows;
	    rc[1] = nc;

	    for (i=0; i<np; i++) {
		if (i == np - 1 && rem > 0) {
		    rc[1] += rem;
		    n += m->rows * rem;
		}
		if (i == root) {
		    err = scatter_to_self(rc, val, recvm);
		} else {
		    err = mpi_send(rc, 2, mpi_int, i, TAG_MATRIX_DIM, 
				   mpi_comm_world);
		    err = mpi_send(val, n, mpi_double, i, TAG_MATRIX_VAL, 
				   mpi_comm_world);
		}
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

static double mpi_dt0;

void gretl_mpi_stopwatch_init (void)
{
    mpi_dt0 = mpi_wtime();
}

double gretl_mpi_stopwatch (void)
{
    double dt1 = mpi_wtime();
    double x = dt1 - mpi_dt0;

    mpi_dt0 = dt1;

    return x;
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
