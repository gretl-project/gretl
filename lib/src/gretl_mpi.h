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
    GRETL_MPI_VCAT
} Gretl_MPI_Op;

int gretl_MPI_init (void);

int gretl_mpi_initialized (void);

int gretl_mpi_rank (void);

int gretl_mpi_reduce (void *sendp, void *recvp,
		      GretlType type, Gretl_MPI_Op op, 
		      int id);

int gretl_matrix_mpi_bcast (gretl_matrix **pm, int id);

int gretl_matrix_mpi_send (const gretl_matrix *m, int dest);

int gretl_matrix_mpi_send_size_known (const gretl_matrix *m, 
				      int dest);

gretl_matrix *gretl_matrix_mpi_receive (int source, 
					int *err);

gretl_matrix *
gretl_matrix_mpi_receive_size_known (int source, 
				     int rows, 
				     int cols,
				     int *err);

int gretl_matrix_mpi_send_cols (const gretl_matrix *m, 
				int j, int ncols);

void gretl_mpi_stopwatch_init (void);

double gretl_mpi_stopwatch (void);

#endif /* GRETL_MPI_H */
