/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef KALMAN_H_
#define KALMAN_H_

typedef struct kalman_ kalman;

kalman *kalman_new (const gretl_matrix *S, const gretl_matrix *P,
		    const gretl_matrix *F, const gretl_matrix *A,
		    const gretl_matrix *H, const gretl_matrix *Q,
		    const gretl_matrix *R, const gretl_matrix *y,
		    int *err);

void kalman_free (kalman *K);

int kalman_forecast (kalman *K, const gretl_matrix *x, gretl_matrix *E);

double kalman_get_loglik (const kalman *K);

double kalman_get_s2 (const kalman *K);

void kalman_set_s2 (kalman *K, double s2);

int kalman_set_initial_state_vector (kalman *K, const gretl_matrix *S);

int kalman_set_initial_MSE_matrix (kalman *K, const gretl_matrix *P);

#endif /* KALMAN_H_ */


