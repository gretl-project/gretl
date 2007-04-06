/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2007 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H
 
int bootstrap_analysis (MODEL *pmod, int p, int B, const double **Z,
			const DATAINFO *pdinfo, gretlopt opt,
			PRN *prn);

int bootstrap_test_restriction (MODEL *pmod, gretl_matrix *R, 
				gretl_matrix *q, double test, int g,
				const double **Z, const DATAINFO *pdinfo, 
				PRN *prn);

int bootstrap_ok (int ci);

int bootstrap_save_data (const char *fname);

#endif /* BOOTSTRAP_H */
