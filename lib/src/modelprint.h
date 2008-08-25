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

typedef struct coeff_fmt_ coeff_fmt;

struct coeff_fmt_ {
    int namelen;
    int blmax;
    int brmax;
    int slmax;
    int srmax;
};

typedef struct model_coeff_ model_coeff;

struct model_coeff_ {
    double b;
    double se;
    double tval;
    double pval;
    double slope;
    double lo;
    double hi;
    int show_pval;
    int df_pval;
    char name[32];
    coeff_fmt *fmt;
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
