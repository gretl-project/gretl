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

/* compare.h for gretl */

#include <stdio.h>

enum aux_codes {
    NONE,
    AUX_SQ,
    AUX_LOG,
    AUX_CHOW,
    AUX_ADD,
    AUX_AR,
    AUX_WHITE,
    AUX_COINT,
    AUX_ARCH,
    AUX_ADF,
    AUX_OMIT
};

typedef struct {
    int m1;        /* ID for first model */
    int m2;        /* ID for second model */
    int ols;       /* was the first model estimated via OLS? */
    int discrete;  /* logit or probit model? */
    int dfn;
    int dfd;
    double F;
    double chisq;
    double trsq;
    int score;
} COMPARE;

/* functions follow */
 
int auxreg (int *addvars, 
	    MODEL *orig, MODEL *new, int *model_count, 
	    double **pZ, DATAINFO *pdinfo, 
	    const int aux_code, 
	    print_t *prn, GRETLTEST *test);

int handle_omit (int *omitvars, MODEL *orig, MODEL *new, 
		 int *model_count, 
		 double **pZ, DATAINFO *pdinfo, 
		 print_t *prn);

int autocorr_test (MODEL *pmod, 
		   double **pZ, DATAINFO *pdinfo, 
		   print_t *prn, GRETLTEST *test);

int chow_test (const char *line, MODEL *pmod, 
	       double **pZ, DATAINFO *pdinfo, 
	       print_t *prn, char *msg, 
	       GRETLTEST *test);

int cusum_test (MODEL *pmod, 
		double **pZ, DATAINFO *pdinfo, 
		print_t *prn, char *msg, 
		const PATHS *ppaths, 
		GRETLTEST *test);
