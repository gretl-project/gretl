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

typedef enum {
    NONE,      /* not an auxiliary regression */
    AUX_SQ,    /* aux. regression for nonlinearity (squared terms) */
    AUX_LOG,   /* aux. regression for nonlinearity (log terms) */
    AUX_CHOW,  /* aux. regression for Chow test */
    AUX_ADD,   /* aux. regression for adding variables */
    AUX_AR,    /* aux. regression for autocorrelation test */
    AUX_WHITE, /* aux. regression for heteroskedasticity (White's test) */
    AUX_COINT, /* aux. regression for cointegreation test */
    AUX_ARCH,  /* aux. regression for ARCH test */
    AUX_ADF,   /* aux. regression for augmented Dickey-Fuller test */
    AUX_OMIT   /* aux. regression for omitting variables */
} aux_codes;

typedef struct {
    int m1;        /* ID for first model */
    int m2;        /* ID for second model */
    int ols;       /* was the first model estimated via OLS? */
    int discrete;  /* logit or probit model? */
    int dfn;       /* numerator degrees of freedom */
    int dfd;       /* denominator degrees of freedom */ 
    double F;      /* F test statistic */
    double chisq;  /* Chi-square test statistic */
    double trsq;   /* T*R^2 test statistic */
    int score;     /* "cases correct" for discrete models */
} COMPARE;

/* functions follow */
 
int auxreg (LIST addvars, 
	    MODEL *orig, MODEL *new, int *model_count, 
	    double ***pZ, DATAINFO *pdinfo, 
	    int aux_code, 
	    PRN *prn, GRETLTEST *test);

int omit_test (LIST omitvars, MODEL *orig, MODEL *new, 
	       int *model_count, 
	       double ***pZ, DATAINFO *pdinfo, 
	       PRN *prn);

int autocorr_test (MODEL *pmod, int order,
		   double ***pZ, DATAINFO *pdinfo, 
		   PRN *prn, GRETLTEST *test);

int chow_test (const char *line, MODEL *pmod, 
	       double ***pZ, DATAINFO *pdinfo, 
	       PRN *prn, GRETLTEST *test);

int cusum_test (MODEL *pmod, 
		double ***pZ, DATAINFO *pdinfo, 
		PRN *prn, 
		PATHS *ppaths, 
		GRETLTEST *test);

int hausman_test (MODEL *pmod, 
		  double ***pZ, DATAINFO *pdinfo, 
		  const PATHS *ppaths, PRN *prn);

int mp_ols (const LIST list, double ***pZ, DATAINFO *pdinfo, 
	    const PATHS *ppaths, PRN *prn); 
