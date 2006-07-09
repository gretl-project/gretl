/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef PVALUES_H
#define PVALUES_H

double binomial_cdf (int k, int n, double p);

double binomial_pvalue (int k, int n, double p);

double x_factorial (double x);

double log_x_factorial (double x);

double normal_pvalue_2 (double x);

double normal_pvalue_1 (double x);

double t_pvalue_2 (double x, int df);

double fdist (double x, int dfn, int dfd);

double chisq (double x, int df);

double normal_cdf (double x);

double normal_pdf (double x);

double log_normal_pdf (double x);

double gamma_dist (double s1, double s2, double x, int control);

double tcrit95 (int df);

double rhocrit95 (int n);

double batch_pvalue (const char *str, 
		     const double **Z, const DATAINFO *pdinfo, 
                     PRN *prn, gretlopt opt);

void interact_pvalue (void);

double f_crit_a (double a, int df1, int df2);

int print_critical (const char *line, PRN *prn);

double genr_get_critical (const char *line, const double **Z, 
			  const DATAINFO *pdinfo);

#endif /* PVALUES_H */


