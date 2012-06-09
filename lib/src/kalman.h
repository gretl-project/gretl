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
    KALMAN_ARMA_LL = 1 << 0, /* is the filter being used for ARMA estimation ? */
    KALMAN_AVG_LL  = 1 << 1, /* store total likelihood or average? */
    KALMAN_USER    = 1 << 2, /* user-defined filter? */
    KALMAN_DIFFUSE = 1 << 3, /* using diffuse P_{1|0} */
    KALMAN_FORWARD = 1 << 4, /* running forward filtering pass */
    KALMAN_SMOOTH  = 1 << 5, /* preparing for smoothing pass */
    KALMAN_SIM     = 1 << 6, /* running simulation */
    KALMAN_CROSS   = 1 << 7, /* cross-correlated disturbances */
    KALMAN_ETT     = 1 << 8, /* ARMA: producing \epsilon{t|t} estimates */
    KALMAN_CHECK   = 1 << 9  /* checking user-defined matrices */
};

typedef struct kalman_ kalman;

kalman *kalman_new (const gretl_matrix *S, const gretl_matrix *P,
		    const gretl_matrix *F, const gretl_matrix *A,
		    const gretl_matrix *H, const gretl_matrix *Q,
		    const gretl_matrix *R, const gretl_matrix *y,
		    const gretl_matrix *x, const gretl_matrix *m,
		    gretl_matrix *E, int *err);

void kalman_free (kalman *K);

int kalman_forecast (kalman *K, PRN *prn);

double kalman_get_loglik (const kalman *K);

double kalman_get_arma_variance (const kalman *K);

gretl_matrix *kalman_arma_smooth (kalman *K, int *err);

gretl_matrix *kalman_smooth (kalman *K,
			     gretl_matrix **pP,
			     gretl_matrix **pU,
			     int *err);

int kalman_set_initial_state_vector (kalman *K, const gretl_matrix *S);

int kalman_set_initial_MSE_matrix (kalman *K, const gretl_matrix *P);

void kalman_set_nonshift (kalman *K, int n);

void kalman_set_options (kalman *K, int opts);

int kalman_get_options (kalman *K);

void kalman_attach_data (kalman *K, void *data);

void *kalman_get_data (const kalman *K);

void kalman_attach_printer (kalman *K, PRN *prn);

PRN *kalman_get_printer (const kalman *K);

#ifndef __GTK_DOC_IGNORE__

int kalman_parse_line (const char *line, const DATASET *dset, 
		       gretlopt opt);

double user_kalman_get_loglik (void);

gretl_matrix *user_kalman_get_matrix (int idx, int *err);

double user_kalman_get_s2 (void);

int user_kalman_get_time_step (void);

int user_kalman_run (const char *E, const char *V, const char *S,
		     const char *P, const char *G, const DATASET *dset, 
		     PRN *prn, int *errp);

gretl_matrix *user_kalman_smooth (const char *Pname, const char *Uname,
				  int *err);

gretl_matrix *user_kalman_simulate (const gretl_matrix *V, 
				    const gretl_matrix *W,
				    const char *Sname, 
				    PRN *prn, int *err);

void kalman_cleanup (void);

int delete_kalman (PRN *prn);

#endif /* __GTK_DOC_IGNORE__ */

#endif /* KALMAN_H_ */
