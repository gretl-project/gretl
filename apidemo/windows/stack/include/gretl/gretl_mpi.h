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

#ifndef GRETL_MPI_H
#define GRETL_MPI_H

typedef enum {
    GRETL_MPI_SUM = 1,
    GRETL_MPI_PROD,
    GRETL_MPI_MAX,
    GRETL_MPI_MIN,
    GRETL_MPI_HCAT,
    GRETL_MPI_VCAT,
    GRETL_MPI_ACAT,
    GRETL_MPI_HSPLIT,
    GRETL_MPI_VSPLIT
} Gretl_MPI_Op;

int gretl_MPI_init (void);

int gretl_mpi_initialized (void);

int gretl_mpi_rank (void);

int gretl_mpi_n_processes (void);

int gretl_mpi_barrier (void);

int gretl_mpi_bcast (void *p, GretlType type, int root);

int gretl_mpi_send (void *p, GretlType type, int dest);

int gretl_matrix_mpi_send (const gretl_matrix *m, int dest);

gretl_matrix *gretl_matrix_mpi_receive (int source, int *err);

int gretl_matrix_mpi_fill (gretl_matrix **pm, int source);

int gretl_matrix_mpi_reduce (gretl_matrix *mymat,
			     gretl_matrix **pm,
			     Gretl_MPI_Op op,
			     int root,
			     gretlopt opt);

int gretl_scalar_mpi_reduce (double x, 
			     double *xp,
			     Gretl_MPI_Op op,
			     int root,
			     gretlopt opt);

int gretl_matrix_mpi_scatter (const gretl_matrix *m,
			      gretl_matrix **recvm,
			      Gretl_MPI_Op op,
			      int root);

int gretl_array_mpi_reduce (gretl_array *sa,
			    gretl_array **pa,
			    Gretl_MPI_Op op,
			    int root);

void *gretl_mpi_receive (int source, GretlType *type,
			 int *err);

double gretl_mpi_time (void);

#endif /* GRETL_MPI_H */
