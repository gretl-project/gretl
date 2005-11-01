/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef GRETL_NLS_H
#define GRETL_NLS_H

typedef struct _nls_spec nls_spec;

nls_spec *nls_spec_new (int ci, const DATAINFO *pdinfo);

void nls_spec_destroy (nls_spec *spec);

int 
nls_spec_add_param_with_deriv (nls_spec *spec, 
			       const char *dstr,
			       const double **Z, 
			       const DATAINFO *pdinfo);

int 
nls_spec_set_regression_function (nls_spec *spec, 
				  const char *fnstr, 
				  const DATAINFO *pdinfo);

void nls_spec_set_t1_t2 (nls_spec *spec, int t1, int t2);

int nls_parse_line (int ci, const char *line, const double **Z,
		    const DATAINFO *pdinfo, PRN *prn);

MODEL nls (double ***pZ, DATAINFO *pdinfo, gretlopt opt, PRN *prn);

MODEL model_from_nls_spec (nls_spec *spec, double ***pZ, 
			   DATAINFO *pdinfo, gretlopt opt, 
			   PRN *prn);

double get_default_nls_toler (void);

#endif /* GRETL_NLS_H */
