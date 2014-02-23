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

#ifndef WIN32
# include <dlfcn.h>
#endif

/* Support for MPI in libgretl. On systems other than MS Windows
   We get the MPI symbols that we need from the address space of
   the calling program to avoid introducing a hard dependency on 
   libmpi.

   To use functions in this translation unit from elsewhere in
   libgretl one must first call gretl_MPI_init(), and then
   guard subsequent calls with "if (gretl_mpi_initialized())".
   Otherwise you'll get a segfault or (on Windows) a fatal MPI
   error.

   For future reference, the MPI variants define the following
   specific symbols in mpi.h:

     Open MPI: OMPI_MAJOR_VERSION
     MPICH: MPICH_VERSION
     MS-MPI: MSMPI_VER
*/

#ifdef OMPI_MAJOR_VERSION
/* external constants that need to be loaded from the 
   Open MPI library */
static MPI_Comm mpi_comm_world;
static MPI_Datatype mpi_double;
static MPI_Datatype mpi_int;
static MPI_Op mpi_max;
static MPI_Op mpi_sum;
#else
/* It seems that MPICH and MS-MPI just define these symbols
   as integer values in the header */
# define mpi_comm_world MPI_COMM_WORLD
# define mpi_double     MPI_DOUBLE
# define mpi_int        MPI_INT
# define mpi_max        MPI_MAX
# define mpi_sum        MPI_SUM
#endif

#ifdef WIN32 /* msmpi.dll loaded */

#define mpi_comm_rank MPI_Comm_rank
#define mpi_error_class MPI_Error_class
#define mpi_error_string MPI_Error_string
#define mpi_reduce MPI_Reduce
#define mpi_bcast MPI_Bcast
#define mpi_send MPI_Send
#define mpi_recv MPI_Recv
#define mpi_wtime MPI_Wtime
#define mpi_initialized MPI_Initialized
	       
int gretl_MPI_init (void)
{
    int initted;

    /* this is a no-op, provided that MPI itself
       is initialized */

    mpi_initialized(&initted);

    return initted ? 0 : E_EXTERNAL;
}

#else /* !WIN32, use dlsym() */

#define MPI_DEBUG 0

static void *MPIhandle;       /* handle to the MPI library */
static int gretl_MPI_err;     /* initialization error record */
static int gretl_MPI_initted; /* are we initialized or not? */

/* renamed, pointerized versions of the MPI functions we need */
static int (*mpi_comm_rank) (MPI_Comm, int *);
static int (*mpi_error_class) (int, int *);
static int (*mpi_error_string) (int, char *, int *);
static int (*mpi_reduce) (void *, void *, int, MPI_Datatype, MPI_Op,
			  int, MPI_Comm);
static int (*mpi_bcast) (void *, int, MPI_Datatype, int, MPI_Comm);
static int (*mpi_send) (void *, int, MPI_Datatype, int, int, MPI_Comm);	       
static int (*mpi_recv) (void *, int, MPI_Datatype, int, int, MPI_Comm,
			MPI_Status *);
static double (*mpi_wtime) (void);
static int (*mpi_initialized) (int *);

static void *mpiget (void *handle, const char *name, int *err)
{
    void *p = dlsym(handle, name);
    
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

    MPIhandle = dlopen(NULL, RTLD_NOW);
    if (MPIhandle == NULL) {
	err = E_EXTERNAL;
	goto bailout;
    } 

    mpi_comm_rank    = mpiget(MPIhandle, "MPI_Comm_rank", &err);
    mpi_error_class  = mpiget(MPIhandle, "MPI_Error_class", &err);
    mpi_error_string = mpiget(MPIhandle, "MPI_Error_string", &err);
    mpi_reduce       = mpiget(MPIhandle, "MPI_Reduce", &err);
    mpi_bcast        = mpiget(MPIhandle, "MPI_Bcast", &err);
    mpi_send         = mpiget(MPIhandle, "MPI_Send", &err);
    mpi_recv         = mpiget(MPIhandle, "MPI_Recv", &err);
    mpi_wtime        = mpiget(MPIhandle, "MPI_Wtime", &err);
    mpi_initialized  = mpiget(MPIhandle, "MPI_Initialized", &err);

#ifdef OMPI_MAJOR_VERSION
    if (!err) {
	mpi_comm_world = (MPI_Comm) mpiget(MPIhandle, "ompi_mpi_comm_world", &err);
	mpi_double     = (MPI_Datatype) mpiget(MPIhandle, "ompi_mpi_double", &err);
	mpi_int        = (MPI_Datatype) mpiget(MPIhandle, "ompi_mpi_int", &err);
	mpi_max        = (MPI_Op) mpiget(MPIhandle, "ompi_mpi_op_max", &err);
	mpi_sum        = (MPI_Op) mpiget(MPIhandle, "ompi_mpi_op_sum", &err);
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

#endif /* WIN32 or not */

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

int gretl_matrix_mpi_reduce (gretl_matrix *m, Gretl_MPI_Op op,
			     double *global_x, int id)
{
    MPI_Op mpi_op;
    double local_x = 0.0;

    if (op == GRETL_MPI_SUM) {
	mpi_op = mpi_sum;
    } else if (op == GRETL_MPI_MAX) {
	mpi_op = mpi_max;
    } else {
	return E_DATA;
    }

    if (id > 0) {
	int i, n = m->rows * m->cols;

	if (mpi_op == mpi_sum) {
	    for (i=0; i<n; i++) {
		local_x += m->val[i];
	    }
	} else if (mpi_op == mpi_max) {
	    local_x = m->val[0];
	    for (i=1; i<n; i++) {
		if (m->val[i] > local_x) {
		    local_x = m->val[i];
		}
	    }
	}
	/* implement other ops here */
    }

    mpi_reduce(&local_x, global_x, 1, mpi_double, mpi_op, 
	       0, mpi_comm_world);

    return 0;
}

int gretl_matrix_mpi_bcast (gretl_matrix **pm, int id)
{
    gretl_matrix *m = NULL;
    int rc[2];
    int err;

    if (id == 0) {
	/* root is doing the broadcasting */
	m = *pm;
	rc[0] = m->rows;
	rc[1] = m->cols;
    }

    /* broadcast the matrix dimensions first */
    err = mpi_bcast(rc, 2, mpi_int, 0, mpi_comm_world);

    if (!err && id > 0) {
	/* everyone but root needs to allocate space */
	*pm = m = gretl_matrix_alloc(rc[0], rc[1]);
	if (m == NULL) {
	    return E_ALLOC;
	}
    }

    if (!err) {
	/* broadcast the matrix content */
	int n = rc[0] * rc[1];

	err = mpi_bcast(m->val, n, mpi_double, 0, mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}

int gretl_matrix_mpi_send (const gretl_matrix *m, int dest)
{
    int rc[2] = {m->rows, m->cols};
    int err;

    err = mpi_send(rc, 2, mpi_int, dest, 0, mpi_comm_world);

    if (!err) {
	int n = m->rows * m->cols;

	err = mpi_send(m->val, n, mpi_double, dest, 0, mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
    }    

    return err;
}

int gretl_matrix_mpi_send_size_known (const gretl_matrix *m, 
				      int dest)
{
    int n = m->rows * m->cols;
    int err;

    err = mpi_send(m->val, n, mpi_double, dest, 0, mpi_comm_world);

    if (err) {
	gretl_mpi_error(&err);
    }    

    return err;
}

gretl_matrix *gretl_matrix_mpi_receive (int source, 
					int *err)
{
    gretl_matrix *m = NULL;
    MPI_Status status;
    int rc[2];

    *err = mpi_recv(rc, 2, mpi_int, source, 0, mpi_comm_world, &status);

    if (!*err) {
	m = gretl_matrix_alloc(rc[0], rc[1]);
	if (m == NULL) {
	    *err = E_ALLOC;
	} else {
	    int n = rc[0] * rc[1];

	    *err = mpi_recv(m->val, n, mpi_double, source, 0, 
			    mpi_comm_world, &status);
	    if (*err) {
		gretl_mpi_error(err);
	    }
	}
    }

    return m;
}

gretl_matrix *
gretl_matrix_mpi_receive_size_known (int source, int rows, int cols,
				     int *err)
{
    gretl_matrix *m = NULL;
    MPI_Status status;

    m = gretl_matrix_alloc(rows, cols);
    if (m == NULL) {
	*err = E_ALLOC;
    } else {
	*err = mpi_recv(m->val, rows * cols, mpi_double, source, 0, 
			mpi_comm_world, &status);
	if (*err) {
	    gretl_mpi_error(err);
	}	
    }

    return m;
}

int gretl_matrix_mpi_send_cols (const gretl_matrix *m, 
				int j, int ncols)
{
    int rc[2] = {m->rows, ncols};
    int n = m->rows * ncols;
    int err;

    err = mpi_send(rc, 2, mpi_int, j, 0, mpi_comm_world);

    if (!err) {
	void *ptr = m->val + (j-1)*n;

	err = mpi_send(ptr, n, mpi_double, j, 0, mpi_comm_world);
    }

    if (err) {
	gretl_mpi_error(&err);
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
    int id;

    mpi_comm_rank(mpi_comm_world, &id);
    return id;
}

int gretl_mpi_initialized (void)
{
    static int ret = -1;

#ifndef WIN32
    if (!gretl_MPI_initted) {
	return 0;
    }
#endif

    if (ret < 0) {
	mpi_initialized(&ret);
	ret = ret > 0;
    }

    return ret;
}


