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

#ifndef RANDOM_H
#define RANDOM_H

void gretl_rand_set_seed (unsigned int seed);

double gretl_one_snormal (void);

void gretl_rand_uniform (double *a, int t1, int t2);

int gretl_rand_uniform_minmax (double *a, int t1, int t2,
			       double min, double max);

void gretl_rand_normal (double *a, int t1, int t2);

int gretl_rand_normal_full (double *a, int t1, int t2,
			    double mean, double sd);

int gretl_rand_chisq (double *a, int t1, int t2, int v);

int gretl_rand_student (double *a, int t1, int t2, int v);

int gretl_rand_F (double *a, int t1, int t2, int v1, int v2);

int gretl_rand_binomial (double *a, int t1, int t2, int n, double p);

void gretl_rand_poisson (double *a, int t1, int t2, double *m,
			 int vec);

int gretl_rand_gamma (double *a, int t1, int t2,  
		      double shape, double scale);

unsigned int gretl_rand_int_max (unsigned int max);

unsigned int gretl_rand_int (void);

void gretl_rand_init (void);

void gretl_rand_free (void);

unsigned int gretl_rand_get_seed (void);

#endif /* RANDOM_H */

