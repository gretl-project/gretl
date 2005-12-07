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

#ifndef GENERATE_H
#define GENERATE_H

typedef struct _GENERATOR GENERATOR;

int generate (const char *line, double ***pZ, DATAINFO *pdinfo, gretlopt opt); 

GENERATOR *
genr_compile (const char *line, double ***pZ, DATAINFO *pdinfo, gretlopt opt);

int execute_genr (GENERATOR *genr, int oldv);

void destroy_genr (GENERATOR *genr);

int genr_get_varnum (const GENERATOR *genr);

int genr_get_err (const GENERATOR *genr);

int varindex (const DATAINFO *pdinfo, const char *varname);

int genr_fit_resid (const MODEL *pmod, 
		    double ***pZ, DATAINFO *pdinfo,
		    int code, int undo);

int get_generated_value (const char *argv, double *val,
			 double ***pZ, DATAINFO *pdinfo,
			 int t);

int gretl_reserved_word (const char *str);

int genr_function_from_string (const char *s);

int get_t_from_obs_string (char *s, const double **Z, 
			   const DATAINFO *pdinfo);

#endif /* GENERATE_H */

