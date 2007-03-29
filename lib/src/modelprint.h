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

/* modelprint.h for gretl */

#ifndef MODELPRINT_H
#define MODELPRINT_H

typedef struct model_coeff_ model_coeff;

struct model_coeff_ {
    double b;
    double se;
    double tval;
    double pval;
    double slope;
    int show_pval;
    int df_pval;
    char name[32];
};

enum {
    COEFF_HEADING_VARNAME,
    COEFF_HEADING_PARAM
};

int printmodel (MODEL *pmod, const DATAINFO *pdinfo, gretlopt opt,
		PRN *prn);

const char *estimator_string (const MODEL *pmod, PRN *prn);

void print_model_vcv_info (const MODEL *pmod, PRN *prn);

int ols_print_anova (const MODEL *pmod, PRN *prn);

void print_coeff_heading (int mode, PRN *prn);

void model_coeff_init (model_coeff *mc);

void print_coeff (const model_coeff *mc, PRN *prn);

void print_arch_coeffs (const double *a, const double *se,
			int T, int order, PRN *prn, int aux);

#endif /* MODELPRINT_H */
