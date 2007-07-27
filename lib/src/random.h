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

#ifndef RANDOM_H
#define RANDOM_H

void gretl_rand_set_seed (unsigned int seed);

double gretl_one_snormal (void);

void gretl_uniform_dist (double *a, int t1, int t2);

int gretl_uniform_dist_minmax (double *a, int t1, int t2,
			       double min, double max);

void gretl_normal_dist (double *a, int t1, int t2);

int gretl_normal_dist_with_params (double *a, int t1, int t2,
				   double mean, double sd);

int gretl_chisq_dist (double *a, int t1, int t2, int v);

int gretl_t_dist (double *a, int t1, int t2, int v);

int gretl_binomial_dist (double *a, int t1, int t2, int n, double p);

void gretl_poisson_dist (double *a, int t1, int t2, double *m,
			 int vec);

int gretl_gamma_dist (double *a, int t1, int t2,  
		      double shape, double scale);

unsigned int gretl_rand_int_max (unsigned int max);

unsigned int gretl_rand_int (void);

void gretl_rand_init (void);

void gretl_rand_free (void);

unsigned int get_gretl_random_seed (void);

#endif /* RANDOM_H */

