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

#include <complex.h>

#define complex_scalar(m) (m->is_complex && m->rows == 1 && m->cols == 1)

gretl_matrix *gretl_matrix_fft (const gretl_matrix *y, int *err);

gretl_matrix *gretl_matrix_ffti (const gretl_matrix *y, int *err);

gretl_matrix *gretl_cmatrix_multiply (const gretl_matrix *A,
				      const gretl_matrix *B,
				      int *err);

gretl_matrix *gretl_cmatrix_AHB (const gretl_matrix *A,
				 const gretl_matrix *B,
				 int *err);

gretl_matrix *gretl_cmatrix_inverse (const gretl_matrix *A, int *err);

gretl_matrix *gretl_cmatrix_ginv (const gretl_matrix *A, int *err);

gretl_matrix *gretl_zheev (const gretl_matrix *A, int eigenvecs,
			   int *err);

gretl_matrix *gretl_zgeev (const gretl_matrix *A,
			   gretl_matrix *VL,
			   gretl_matrix *VR,
			   int *err);

int gretl_cmatrix_SVD (const gretl_matrix *x, gretl_matrix **pu,
		       gretl_vector **ps, gretl_matrix **pvt,
		       int smod);

gretl_matrix *gretl_complex_fft (const gretl_matrix *A, int inverse,
				 int *err);

gretl_matrix *gretl_cmatrix (const gretl_matrix *Re,
			     const gretl_matrix *Im,
			     double ival, int *err);

gretl_matrix *gretl_cxtract (const gretl_matrix *A, int im,
			    int *err);

gretl_matrix *gretl_ctrans (const gretl_matrix *A,
			    int conjugate, int *err);

int gretl_ctrans_in_place (gretl_matrix *A);

gretl_matrix *cmatrix_add_sub (const gretl_matrix *A,
			       const gretl_matrix *B,
			       int sgn, int *err);

int apply_cmatrix_dfunc (gretl_matrix *targ,
			 const gretl_matrix *src,
			 double (*dfunc) (double complex));

int apply_cmatrix_cfunc (gretl_matrix *targ,
			 const gretl_matrix *src,
			 double complex (*cfunc) (double complex));

int apply_cmatrix_unary_op (gretl_matrix *targ,
			    const gretl_matrix *src,
			    int op);

gretl_matrix *gretl_cmatrix_determinant (const gretl_matrix *X,
					 int log, int *err);

gretl_matrix *gretl_cmatrix_trace (const gretl_matrix *X,
				   int *err);

int gretl_cmatrix_set_diagonal (gretl_matrix *targ,
				const gretl_matrix *src,
				double x);

gretl_matrix *gretl_cmatrix_hdprod (const gretl_matrix *A,
				    const gretl_matrix *B,
				    int *err);

gretl_matrix *gretl_cmatrix_kronecker (const gretl_matrix *A,
				       const gretl_matrix *B,
				       int *err);

gretl_matrix *gretl_cmatrix_switch (const gretl_matrix *m,
				    int to_new, int *err);

gretl_matrix *gretl_cmatrix_vector_stat (const gretl_matrix *m,
					 GretlVecStat vs, int rowwise,
					 int *err);

int gretl_cmatrix_fill (gretl_matrix *m, double complex z);

gretl_matrix *scalar_to_complex (double x, int *err);

gretl_matrix *two_scalars_to_complex (double xr, double xi,
				      int *err);

gretl_matrix *complex_scalar_to_mat (double complex z, int *err);

int complex_matrix_print_range (const gretl_matrix *A,
				const char *name,
				int rmin, int rmax,
				PRN *prn);

int complex_matrix_print (const gretl_matrix *A,
			  const char *name,
			  PRN *prn);

int complex_matrix_printf (const gretl_matrix *A,
			   const char *fmt,
			   PRN *prn);

gretl_matrix *gretl_cmatrix_dot_op (const gretl_matrix *a,
				    const gretl_matrix *b,
				    int op, int *err);

gretl_matrix *cmatrix_get_element (const gretl_matrix *M,
				   int i, int *err);

#endif /* GRETL_CMATRIX_H */
