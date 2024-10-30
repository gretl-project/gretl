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

#ifndef KALMAN_H_
#define KALMAN_H_

typedef struct kalman_ kalman;

kalman *kalman_new_minimal (gretl_matrix *M[], int copy[],
			    int nmat, int dkvar, int *err);

void kalman_free (kalman *K);

/* the "raw" C API */

kalman *kalman_new (gretl_matrix *a, gretl_matrix *P,
		    gretl_matrix *T, gretl_matrix *BT,
		    gretl_matrix *ZT, gretl_matrix *HH,
		    gretl_matrix *GG, gretl_matrix *y,
		    gretl_matrix *x, gretl_matrix *mu,
		    gretl_matrix *V, int *err);

int kalman_forecast (kalman *K, PRN *prn);

int kalman_run (kalman *K, PRN *prn, int *errp);

int is_kalman_bundle (gretl_bundle *b);

gretl_matrix *kalman_smooth (kalman *K, gretlopt opt,
			     PRN *prn, int *err);

double kalman_get_loglik (const kalman *K);

double kalman_get_arma_variance (const kalman *K);

int kalman_set_initial_state_vector (kalman *K, const gretl_vector *a);

int kalman_set_initial_MSE_matrix (kalman *K, const gretl_matrix *P);

void kalman_set_arma_ll (kalman *K);

void kalman_attach_data (kalman *K, void *data);

void *kalman_get_data (const kalman *K);

void kalman_attach_printer (kalman *K, PRN *prn);

PRN *kalman_get_printer (const kalman *K);

/* end "raw" C API */

#ifndef __GTK_DOC_IGNORE__

int kalman_bundle_filter (gretl_bundle *b, PRN *prn, int *errp);

int kalman_bundle_smooth (gretl_bundle *b, int dist, PRN *prn);

gretl_matrix *kalman_bundle_simulate (gretl_bundle *b,
				      const gretl_matrix *U, 
				      int get_state,
				      PRN *prn, int *err);

gretl_matrix *kalman_bundle_simdata (gretl_bundle *b,
				     const gretl_matrix *U,
				     PRN *prn, int *err);

int maybe_set_kalman_element (void *kptr,
			      const char *key,
			      void *vptr,
			      GretlType vtype,
			      int copy,
			      int *err);

void *maybe_retrieve_kalman_element (void *kptr,
				     const char *key,
				     GretlType *type,
				     int *reserved,
				     int *ownit,
				     int *err);

int maybe_delete_kalman_element (void *kptr,
				 const char *key,
				 int *err);

int print_kalman_bundle_info (void *kptr, PRN *prn);

gretl_bundle *kalman_bundle_copy (const gretl_bundle *src,
				  int *err);

int kalman_serialize (void *kptr, PRN *prn);

gretl_bundle *kalman_deserialize (void *p1, void *p2,
				  int *err);

char **kalman_bundle_get_matrix_names (kalman *K, int *ns);

char **kalman_bundle_get_scalar_names (kalman *K, int *ns);

int kalman_bundle_n_members (gretl_bundle *b);

#endif /* __GTK_DOC_IGNORE__ */

#endif /* KALMAN_H_ */
