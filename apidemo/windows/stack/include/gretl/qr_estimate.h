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

int lapack_cholesky_regress (MODEL *pmod, const DATASET *dset,
			     gretlopt opt);

int qr_tsls_vcv (MODEL *pmod, DATASET *dset, gretlopt opt);

int qr_matrix_hccme (const gretl_matrix *X,
		     const gretl_matrix *h,
		     const gretl_matrix *XTXi,
		     gretl_matrix *d,
		     gretl_matrix *VCV,
		     int hc_version);

double hac_weight (int kern, int h, int i);

double qs_hac_weight (double bt, int i);

int maybe_limit_VAR_coeffs (gretl_matrix *A,
			    gretl_matrix *Y,
			    gretl_matrix *X,
			    gretl_matrix *E);

int newey_west_bandwidth (const gretl_matrix *f, 
			  const gretl_matrix *w,
			  int kern, int prewhitened,
			  int *h, double *bt);

gretl_matrix *HAC_XOX (const gretl_matrix *X,
		       const gretl_matrix *uhat,
		       VCVInfo *vi, int use_prior,
		       int *err);

gretl_matrix *newey_west_OPG (const gretl_matrix *G,
			      int *err);

gretl_matrix *long_run_covariance (const gretl_matrix *X,
				   int demean, int *err);

int set_cluster_vcv_ci (int ci);

#endif  /* QR_ESTIMATE_H */
