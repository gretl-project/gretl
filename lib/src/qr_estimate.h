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

#ifndef QR_ESTIMATE_H
#define QR_ESTIMATE_H

int gretl_qr_regress (MODEL *pmod, DATASET *dset, gretlopt opt);

int qr_tsls_vcv (MODEL *pmod, const DATASET *dset, gretlopt opt);

double hac_weight (int kern, int h, int i);

double qs_hac_weight (double bt, int i);

int newey_west_bandwidth (const gretl_matrix *f, int kern, int *h, double *bt);

gretl_matrix *HAC_XOX (const gretl_matrix *uhat, const gretl_matrix *X,
		       VCVInfo *vi, int *err);

#endif  /* QR_ESTIMATE_H */
