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
    KALMAN_USER    = 1 << 2  /* user-defined filter? */
};

typedef struct kalman_ kalman;

kalman *kalman_new (const gretl_matrix *S, const gretl_matrix *P,
		    const gretl_matrix *F, const gretl_matrix *A,
		    const gretl_matrix *H, const gretl_matrix *Q,
		    const gretl_matrix *R, const gretl_matrix *y,
		    const gretl_matrix *x, gretl_matrix *E,
		    int *err);

void kalman_free (kalman *K);

int kalman_forecast (kalman *K);

double kalman_get_loglik (const kalman *K);

double user_kalman_get_loglik (void);

double kalman_get_arma_variance (const kalman *K);

int kalman_set_initial_state_vector (kalman *K, const gretl_matrix *S);

int kalman_set_initial_MSE_matrix (kalman *K, const gretl_matrix *P);

void kalman_set_nonshift (kalman *K, int n);

void kalman_set_options (kalman *K, int opts);

int kalman_get_options (kalman *K);

int kalman_parse_line (const char *line, gretlopt opt);

int user_kalman_run (const char *E, const char *S, const char *P,
		     const char *L, int *err);

gretl_matrix *user_kalman_smooth (const char *Pname, int *err);

void destroy_user_kalman (void);

int delete_kalman (PRN *prn);

#endif /* KALMAN_H_ */
