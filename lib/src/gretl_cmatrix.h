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

#ifndef GRETL_CMATRIX_H
#define GRETL_CMATRIX_H

gretl_matrix *gretl_zheev (gretl_array *A, gretl_array *V, int *err);

gretl_array *gretl_zgetri (gretl_array *A, int *err);

gretl_array *gretl_zgemm (gretl_array *A, gretl_array *B, int *err);

int complex_matrix_print (gretl_array *A, PRN *prn);

gretl_array *gretl_complex_fft (gretl_array *A, int *err);

#endif /* GRETL_CMATRIX_H */
