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

#ifndef GRETL_NLS_H
#define GRETL_NLS_H

enum {
    C_LOGLIK,
    C_GMM,
    C_OTHER
};

typedef struct _nlspec nlspec;

nlspec *nlspec_new (int ci, const DATAINFO *pdinfo);

void nlspec_destroy (nlspec *spec);

int 
nlspec_add_param_with_deriv (nlspec *spec, 
			     const char *dstr,
			     const double **Z, 
			     const DATAINFO *pdinfo);

int nlspec_add_param_list (nlspec *spec, int np, double *vals,
			   char **names, double ***pZ,
			   DATAINFO *pdinfo);

int 
nlspec_set_regression_function (nlspec *spec, 
				const char *fnstr, 
				const DATAINFO *pdinfo);

void nlspec_set_t1_t2 (nlspec *spec, int t1, int t2);

int nl_parse_line (int ci, const char *line, const double **Z,
		   const DATAINFO *pdinfo, PRN *prn);

MODEL nl_model (double ***pZ, DATAINFO *pdinfo, gretlopt opt, PRN *prn);

MODEL model_from_nlspec (nlspec *spec, double ***pZ, 
			 DATAINFO *pdinfo, gretlopt opt, 
			 PRN *prn);

MODEL ivreg_via_gmm (const int *list, double ***pZ,
		     DATAINFO *pdinfo, gretlopt opt);

double get_default_nls_toler (void);

#endif /* GRETL_NLS_H */
