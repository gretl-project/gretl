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

void gretl_rand_init (void);

void gretl_mpi_rand_init (int n, int self, int single_rng);

void gretl_rand_set_seed (guint64 seed);

void gretl_alt_rand_set_seed (guint64 seed);

guint32 gretl_rand_int (void);

guint32 gretl_alt_rand_int (void);

guint32 gretl_rand_int_max (guint32 max);

int gretl_rand_int_minmax (int *a, int n, int min, int max);

int gretl_alt_rand_int_minmax (int *a, int n, int min, int max);

double gretl_rand_01 (void);

double gretl_one_snormal (void);

void gretl_rand_uniform (double *a, int t1, int t2);

int gretl_rand_uniform_minmax (double *a, int t1, int t2,
			       double min, double max);

int gretl_rand_uniform_int_minmax (double *a, int t1, int t2,
				   int min, int max,
				   gretlopt opt);

void gretl_rand_normal (double *a, int t1, int t2);

int gretl_rand_normal_full (double *a, int t1, int t2,
			    double mean, double sd);

int gretl_rand_chisq (double *a, int t1, int t2, int v);

int gretl_rand_student (double *a, int t1, int t2, double v);

int gretl_rand_F (double *a, int t1, int t2, int v1, int v2);

int gretl_rand_binomial (double *a, int t1, int t2, int n, double p);

int gretl_rand_poisson (double *a, int t1, int t2, const double *m,
			 int vec);

int gretl_rand_exponential (double *a, int t1, int t2, double mu);

int gretl_rand_weibull (double *a, int t1, int t2, double shape,
			double scale);

int gretl_rand_logistic (double *a, int t1, int t2,
			 double loc, double scale);

int gretl_rand_gamma (double *a, int t1, int t2,  
		      double shape, double scale);

double gretl_rand_gamma_one (double shape, double scale);

int gretl_rand_GED (double *a, int t1, int t2, double nu);

int gretl_rand_laplace (double *a, int t1, int t2,
			double mu, double b);

int gretl_rand_beta (double *x, int t1, int t2, 
		     double s1, double s2);

int gretl_rand_beta_binomial (double *x, int t1, int t2, 
			      int n, double s1, double s2);

int gretl_rand_discrete (double *x, int t1, int t2,
                         const gretl_vector *p);

gretl_matrix *gretl_rand_dirichlet (const gretl_vector *a,
				    int n, int *err);

gretl_matrix *halton_matrix (int m, int r, int offset, int *err);

gretl_matrix *inverse_wishart_matrix (const gretl_matrix *S,
				      int v, int *err);

gretl_matrix *inverse_wishart_sequence (const gretl_matrix *S,
					int v, int n, int *err);

guint64 gretl_rand_get_seed (void);

char *gretl_rand_hex_string (int len, int *err);

#endif /* RANDOM_H */

