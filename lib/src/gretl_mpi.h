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

#include <mpi.h>

int gretl_matrix_mpi_reduce (gretl_matrix *m, MPI_Op op,
			     double *global_x, int id);

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

#endif /* GRETL_MPI_H */
