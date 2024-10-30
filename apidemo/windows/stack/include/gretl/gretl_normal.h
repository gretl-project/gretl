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

#ifndef GRETL_NORMAL_H
#define GRETL_NORMAL_H

double invmills (double x);

double bvnorm_cdf (double rho, double a, double b);

gretl_matrix *gretl_GHK (const gretl_matrix *C,
			 const gretl_matrix *A,
			 const gretl_matrix *B,
			 const gretl_matrix *U,
			 int *err);

gretl_matrix *gretl_GHK2 (const gretl_matrix *C,
			  const gretl_matrix *A,
			  const gretl_matrix *B,
			  const gretl_matrix *U,
			  gretl_matrix *dP,
			  int *err);

#endif /* GRETL_NORMAL_H */
