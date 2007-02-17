/* gretl - The Gnu Regression, Econometrics and Time-series Library
 *
 * Copyright (C) Allin Cottrell
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

#ifndef TEXPRINT_H
#define TEXPRINT_H

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

int tex_print_equation (const MODEL *pmod, const DATAINFO *pdinfo, 
			gretlopt opt, PRN *prn);

int tex_print_model (MODEL *pmod, const DATAINFO *pdinfo, 
		     gretlopt opt, PRN *prn);

void tex_coeff_table_start (const char *col1, const char *col2,
			    int binary, PRN *prn);

void tex_coeff_table_end (PRN *prn);

void tex_print_coeff (const model_coeff *mc, PRN *prn);

void make_tex_coeff_name (const MODEL *pmod, const DATAINFO *pdinfo,
			  int i, char *name);

void tex_print_VECM_omega (GRETL_VAR *vecm, const DATAINFO *pdinfo, PRN *prn);

void tex_print_VECM_coint_eqns (GRETL_VAR *vecm, const DATAINFO *pdinfo, PRN *prn);

void tex_print_VAR_ll_stats (GRETL_VAR *var, PRN *prn);

int texprint (MODEL *pmod, const DATAINFO *pdinfo, char *texfile, 
	      gretlopt opt);

char *tex_escape (char *targ, const char *src);

void tex_dcolumn_double (double x, char *numstr);

char *tex_float_string (double x, int prec, char *targ);

void set_gretl_tex_preamble (void);

void set_tex_use_utf (int s);

void gretl_tex_preamble (PRN *prn, int ams);

void tex_print_obs_marker (int t, const DATAINFO *pdinfo, PRN *prn);

void set_tex_param_format (const char *s);

int tex_using_custom_tabular (void);

const char *tex_column_format (int i);

#endif /* TEXPRINT_H */
