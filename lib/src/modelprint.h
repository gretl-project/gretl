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

/* modelprint.h for gretl */

#ifndef MODELPRINT_H
#define MODELPRINT_H

#define MC_NAMELEN 32

typedef struct model_coeff_ model_coeff;

struct model_coeff_ {
    double b;
    double se;
    double tval;
    double pval;
    double slope;
    double lo;
    double hi;
    int show_tval;
    int show_pval;
    int df_pval;
    int multi;
    char name[MC_NAMELEN];
};

int printmodel (MODEL *pmod, const DATASET *dset, gretlopt opt,
		PRN *prn);

const char *estimator_string (const MODEL *pmod, PRN *prn);

void arma_extra_info (const MODEL *pmod, PRN *prn);

void print_model_vcv_info (const MODEL *pmod, const DATASET *dset,
			   PRN *prn);

int ols_print_anova (const MODEL *pmod, PRN *prn);

int print_coeffs (const double *b,
		  const double *se,
		  const char **names,
		  int nc, int df, int ci,
		  PRN *prn);

int print_model_from_matrices (const gretl_matrix *cs, 
			       const gretl_matrix *adds, 
			       gretl_array *names, int df,
			       gretlopt opt, PRN *prn);

gretlopt get_printmodel_opt (const MODEL *pmod,
			     gretlopt opt);

int print_coeff_intervals (const CoeffIntervals *cf, PRN *prn);

#endif /* MODELPRINT_H */
