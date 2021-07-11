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

#ifndef TEXPRINT_H
#define TEXPRINT_H

int tex_print_equation (const MODEL *pmod, const DATASET *dset, 
			gretlopt opt, PRN *prn);

int tex_print_model (MODEL *pmod, const DATASET *dset, 
		     gretlopt opt, PRN *prn);

int tex_coeff_table_start (const char **cols, gretlopt opt, PRN *prn);

void tex_coeff_table_end (PRN *prn);

void tex_print_coeff (const model_coeff *mc, PRN *prn);

void make_tex_coeff_name (const MODEL *pmod, const DATASET *dset,
			  int i, char *name);

void tex_print_VECM_omega (GRETL_VAR *vecm, const DATASET *dset, PRN *prn);

void tex_print_VECM_coint_eqns (GRETL_VAR *vecm, const DATASET *dset, PRN *prn);

void tex_print_VAR_ll_stats (GRETL_VAR *var, PRN *prn);

int texprint (MODEL *pmod, const DATASET *dset,
	      const char *fname, gretlopt opt);

int rtfprint (MODEL *pmod, const DATASET *dset,
	      const char *fname, gretlopt opt);

int csvprint (MODEL *pmod, const DATASET *dset,
	      const char *fname, gretlopt opt);

char *tex_escape (char *targ, const char *src);

char *tex_escape_new (const char *src);

char *tex_escape_special (char *targ, const char *src);

char *tex_rl_float (double x, char *s, int dig);

char *tex_rl_double (double x, char *s);

void tex_print_double (double x, PRN *prn);

char *tex_sprint_double (double x, char *s);

char *tex_sprint_double_digits (double x, char *s, int dig);

char *tex_modify_exponent (char *s);

void set_gretl_tex_preamble (void);

void set_tex_use_pdf (const char *prog);

int get_tex_use_pdf (void);

void gretl_tex_preamble (PRN *prn, int fmt);

void tex_print_obs_marker (int t, const DATASET *dset, PRN *prn);

int set_tex_param_format (const char *s);

int tex_using_custom_tabular (void);

const char *tex_column_format (int i);

#endif /* TEXPRINT_H */
