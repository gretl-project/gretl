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

#ifndef PVALUES_H
#define PVALUES_H

typedef enum {
    D_NONE = 0,
    D_UNIFORM,
    D_UDISCRT,
    D_NORMAL,
    D_STUDENT,
    D_CHISQ,
    D_SNEDECOR,
    D_BINOMIAL,
    D_POISSON,
    D_EXPON,
    D_WEIBULL,
    D_GAMMA,
    D_GED,
    D_LAPLACE,
    D_BETA,
    D_DW,
    D_BINORM,
    D_JOHANSEN,
    D_BETABIN,
    D_NC_CHISQ,
    D_NC_F,
    D_NC_T,
    D_LOGISTIC,
    D_DIRICHLET,
    D_DISCRETE
} DistCode;

double gammafun (double x);

double lngamma (double x);

double digamma (double x);

double trigamma (double x);

double hypergeo (double a, double b, double c, double x);

double beta_cdf (double a, double b, double x);

double binomial_cdf (double p, int n, int k);

double binomial_cdf_comp (double p, int n, int k);

double binomial_pmf (double p, int n, int k);

double poisson_pmf (double lambda, int k);

double x_factorial (double x);

double log_x_factorial (double x);

double normal_pvalue_2 (double x);

double normal_pvalue_1 (double x);

double student_pvalue_2 (double df, double x);

double student_pvalue_1 (double df, double x);

double chisq_cdf (double df, double x);

double chisq_cdf_comp (double df, double x);

double nc_chisq_cdf (double df, double delta, double x);

double snedecor_cdf (double dfn, double dfd, double x);

double snedecor_cdf_comp (double dfn, double dfd, double x);

double snedecor_critval (double dfn, double dfd, double a);

double snedecor_pvalue_2 (double dfn, double dfd, double x);

double nc_snedecor_cdf (double dfn, double dfd, double delta, double x);

double normal_cdf (double x);

double normal_cdf_inverse (double x);

double normal_cdf_comp (double x);

double student_cdf (double df, double x);

double student_cdf_inverse (double df, double a);

double nc_student_cdf (double df, double delta, double x);

double nc_student_pdf (double df, double delta, double x);

double normal_pdf (double x);

double normal_critval (double a);

double student_critval (double df, double a);

double log_normal_pdf (double x);

double gamma_cdf (double s1, double s2, double x, int control);

double gamma_cdf_comp (double s1, double s2, double x, int control);

double gamma_cdf_inverse (double shape, double scale, double p);

double GED_pdf (double nu, double x);

double GED_cdf (double nu, double x);

double GED_cdf_comp (double nu, double x);

double GED_cdf_inverse (double nu, double a);

double laplace_pdf (double mu, double b, double x);

double laplace_cdf (double mu, double b, double x);

double laplace_cdf_comp (double mu, double b, double x);

double laplace_cdf_inverse (double mu, double b, double a);

double tcrit95 (int df);

double rhocrit (int n, double alpha);

double cephes_gamma (double x);

double cephes_lgamma (double x); 

double gretl_get_pvalue (int dist, const double *parm, double x);

double gretl_get_pdf (int dist, const double *parm, double x);

int gretl_fill_pdf_array (int dist, const double *parm, double *x, int n);

double gretl_get_cdf (int dist, const double *parm, double x);

double gretl_get_cdf_inverse (int dist, const double *parm, double a);

double gretl_get_critval (int dist, const double *parm, double a);

int gretl_fill_random_series (double *x, int dist, 
			      const double *parm,
			      const double *vecp1, 
			      const double *vecp2, 
			      const DATASET *dset);

gretl_matrix *gretl_get_random_matrix (int dist,
				       const double *parm,
				       const double *vecp1,
				       const double *vecp2,
				       int rows, int cols, 
				       int *err);

double gretl_get_random_scalar (int dist, const double *parm,
				int *err);

int batch_pvalue (const char *str, DATASET *dset, PRN *prn);

void print_pvalue (int dist, const double *parm, double x, double pv, PRN *prn);

void print_critval (int dist, const double *parm, double a, double c, PRN *prn);

gretl_matrix *gretl_get_DW (int n, int k, int *err);

int dist_code_from_string (const char *s);

#endif /* PVALUES_H */
