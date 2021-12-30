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

enum {
    KALMAN_USER    = 1 << 0, /* user-defined filter? */
    KALMAN_DIFFUSE = 1 << 1, /* using diffuse P_{1|0} */
    KALMAN_FORWARD = 1 << 2, /* running forward filtering pass */
    KALMAN_SMOOTH  = 1 << 3, /* preparing for smoothing pass */
    KALMAN_SIM     = 1 << 4, /* running simulation */
    KALMAN_CROSS   = 1 << 5, /* cross-correlated disturbances */
    KALMAN_CHECK   = 1 << 6, /* checking user-defined matrices */
    KALMAN_BUNDLE  = 1 << 7, /* kalman is inside a bundle */
    KALMAN_SSFSIM  = 1 << 8  /* on simulation, emulate SsfPack */
};

typedef struct kalman_ kalman;

kalman *kalman_new_minimal (gretl_matrix *M[], int copy[],
			    int nmat, int *err);

void kalman_free (kalman *K);

int kalman_forecast (kalman *K, PRN *prn);

gretl_matrix *kalman_smooth (kalman *K,
			     gretl_matrix **pP,
			     gretl_matrix **pU,
			     int *err);

#ifndef __GTK_DOC_IGNORE__

int kalman_bundle_run (gretl_bundle *b, PRN *prn, int *errp);

int kalman_bundle_smooth (gretl_bundle *b, int dist, PRN *prn);

gretl_matrix *kalman_bundle_simulate (gretl_bundle *b,
				      const gretl_matrix *U, 
				      int get_state,
				      PRN *prn, int *err);

gretl_matrix *kalman_bundle_simdata (gretl_bundle *b,
				     const gretl_matrix *U,
				     PRN *prn, int *err);

/* for interfacing with gretl bundle mechanism */

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
