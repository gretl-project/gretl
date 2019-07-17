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

gretl_matrix *gretl_matrix_fft (const gretl_matrix *y, int *err);

gretl_matrix *gretl_matrix_ffti (const gretl_matrix *y, int *err);

gretl_matrix *gretl_zgemm (const gretl_matrix *A,
			   GretlMatrixMod amod,
			   const gretl_matrix *B,
			   GretlMatrixMod bmod,
			   int *err);

gretl_matrix *gretl_zgetri (const gretl_matrix *A, int *err);

gretl_matrix *gretl_zheev (const gretl_matrix *A, gretl_matrix *V,
			   int *err);

gretl_matrix *gretl_zgeev (const gretl_matrix *A,
			   gretl_matrix *VL,
			   gretl_matrix *VR,
			   int *err);

gretl_matrix *gretl_complex_fft (const gretl_matrix *A, int inverse,
				 int *err);

gretl_matrix *gretl_complex_hprod (const gretl_matrix *A,
				   const gretl_matrix *B,
				   int *err);

gretl_matrix *gretl_cmatrix (const gretl_matrix *Re,
			     const gretl_matrix *Im,
			     int *err);

gretl_matrix *gretl_cxtract (const gretl_matrix *A, int im,
			    int *err);

gretl_matrix *gretl_ctran (const gretl_matrix *A, int *err);

int gretl_ctran_in_place (gretl_matrix *A);

gretl_matrix *gretl_cexp (const gretl_matrix *A, int *err);

int complex_matrix_print (gretl_matrix *A, const char *name, PRN *prn);

int complex_matrix_printf (gretl_matrix *A, const char *fmt, PRN *prn);

#endif /* GRETL_CMATRIX_H */
