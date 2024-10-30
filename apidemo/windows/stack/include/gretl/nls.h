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
    C_SSR,
    C_OTHER
};

/**
 * nlspec:
 *
 * An opaque structure handled only via accessor functions.
 */

typedef struct nlspec_ nlspec;

nlspec *nlspec_new (int ci, const DATASET *dset);

void nlspec_destroy (nlspec *spec);

int nlspec_add_param_with_deriv (nlspec *spec, const char *s);

int nlspec_add_param_list (nlspec *spec, int np, double *vals,
			   char **names);

int aux_nlspec_add_param_list (nlspec *spec, int np, double *vals,
			       char **names);

int 
nlspec_set_regression_function (nlspec *spec, 
				const char *fnstr, 
				const DATASET *dset);

void nlspec_set_t1_t2 (nlspec *spec, int t1, int t2);

int nl_parse_line (int ci, const char *line, 
		   const DATASET *dset, PRN *prn);

int nl_set_smallstep (void);

MODEL nl_model (DATASET *dset, gretlopt opt, PRN *prn);

MODEL model_from_nlspec (nlspec *spec, DATASET *dset, 
			 gretlopt opt, PRN *prn);

MODEL GNR (int *glist, DATASET *gdset, gretlopt opt, PRN *prn);

int finalize_nls_model (MODEL *pmod, nlspec *spec,
			int perfect, int *glist);

int nls_boot_calc (const MODEL *pmod, DATASET *dset,
		   int ft1, int ft2, double *fcerr);

int nl_model_run_aux_genrs (const MODEL *pmod, 
			    DATASET *dset);

double get_default_nls_toler (void);

#endif /* GRETL_NLS_H */
