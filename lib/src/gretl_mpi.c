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

static void gretl_mpi_error (int *err)
{
    char msg1[BUFSIZ], msg2[BUFSIZ];
    int len, errclass;
    int id = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Error_class(*err, &errclass);
    MPI_Error_string(errclass, msg1, &len);
    MPI_Error_string(*err, msg2, &len);
    gretl_errmsg_sprintf("%3d: %s %s\n", id, msg1, msg2);

    *err = E_EXTERNAL;
}

int gretl_matrix_mpi_reduce (gretl_matrix *m, MPI_Op op,
			     double *global_x, int id)
{
    double local_x = 0.0;

    if (id > 0) {
	int i, n = m->rows * m->cols;

	if (op == MPI_SUM) {
	    for (i=0; i<n; i++) {
		local_x += m->val[i];
	    }
	} else if (op == MPI_MAX) {
	    local_x = m->val[0];
	    for (i=1; i<n; i++) {
		if (m->val[i] > local_x) {
		    local_x = m->val[i];
		}
	    }
	}
	/* implement other ops here */
    }

    MPI_Reduce(&local_x, global_x, 1, MPI_DOUBLE, op, 0, MPI_COMM_WORLD);

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
    err = MPI_Bcast(rc, 2, MPI_INT, 0, MPI_COMM_WORLD);

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

	err = MPI_Bcast(m->val, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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

    err = MPI_Send(rc, 2, MPI_INT, dest, 0, MPI_COMM_WORLD);

    if (!err) {
	int n = m->rows * m->cols;

	err = MPI_Send(m->val, n, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
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

    err = MPI_Send(m->val, n, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);

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

    *err = MPI_Recv(rc, 2, MPI_INT, source, 0, MPI_COMM_WORLD, &status);

    if (!*err) {
	m = gretl_matrix_alloc(rc[0], rc[1]);
	if (m == NULL) {
	    *err = E_ALLOC;
	} else {
	    int n = rc[0] * rc[1];

	    *err = MPI_Recv(m->val, n, MPI_DOUBLE, source, 0, 
			    MPI_COMM_WORLD, &status);
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
	*err = MPI_Recv(m->val, rows * cols, MPI_DOUBLE, source, 0, 
			MPI_COMM_WORLD, &status);
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

    err = MPI_Send(rc, 2, MPI_INT, j, 0, MPI_COMM_WORLD);

    if (!err) {
	void *ptr = m->val + (j-1)*n;

	err = MPI_Send(ptr, n, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
    }

    if (err) {
	gretl_mpi_error(&err);
    }

    return err;
}
