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

double gamma_function (double x);

double log_gamma_function (double x);

double binomial_cdf (double p, int n, int k);

double binomial_cdf_comp (double p, int n, int k);

double x_factorial (double x);

double log_x_factorial (double x);

double normal_pvalue_2 (double x);

double normal_pvalue_1 (double x);

double student_pvalue_2 (int df, double x);

double chisq_cdf (int df, double x);

double chisq_cdf_comp (int df, double x);

double snedecor_cdf_comp (int dfn, int dfd, double x);

double snedecor_critval (int dfn, int dfd, double a);

double normal_cdf (double x);

double normal_cdf_inverse (double x);

double normal_cdf_comp (double x);

double student_cdf_inverse (double df, double a);

double normal_pdf (double x);

double normal_critval (double a);

double log_normal_pdf (double x);

double bvnorm_cdf (double a, double b, double rho);

double gamma_cdf_comp (double s1, double s2, double x, int control);

double tcrit95 (int df);

double rhocrit95 (int n);

double cephes_gamma (double x);

double cephes_lgamma (double x); 

double gretl_get_pvalue (char st, const double *p);

double gretl_get_pdf (char st, double *p);

int gretl_fill_pdf_array (char st, double *p, double *x, int n);

double gretl_get_cdf (char st, double *p);

double gretl_get_cdf_inverse (char st, double *p);

double gretl_get_critval (char st, double *p);

double *gretl_get_random_series (char st, const double *p,
				 const DATAINFO *pdinfo,
				 int *err);

int batch_pvalue (const char *str, 
		  double ***pZ, DATAINFO *pdinfo, 
		  PRN *prn);

void print_pvalue (char st, double *p, double pv, PRN *prn);

void print_critval (char st, double *x, double c, PRN *prn);

gretl_matrix *gretl_get_DW (int n, int k, int *err);

#endif /* PVALUES_H */


