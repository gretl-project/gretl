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
#include "gretl_matrix.h"

enum series_save_opts {
    SAVE_LEVERAGE  = 1 << 0,
    SAVE_INFLUENCE = 1 << 1,
    SAVE_DFFITS    = 1 << 2,
    SAVE_TREND     = 1 << 3,
    SAVE_CYCLE     = 1 << 4
};

enum aux_codes {
    AUX_NONE,  /* not an auxiliary regression */
    AUX_SQ,    /* aux. regression for nonlinearity (squared terms) */
    AUX_LOG,   /* aux. regression for nonlinearity (log terms) */
    AUX_CHOW,  /* aux. regression for Chow test */
    AUX_ADD,   /* aux. regression for adding variables */
    AUX_AR,    /* aux. regression for autocorrelation test */
    AUX_SCR,   /* regression showing serial correlation-robust std errs */
    AUX_WHITE, /* aux. regression for heteroskedasticity (White's test) */
    AUX_COINT, /* aux. regression for cointegreation test */
    AUX_ARCH,  /* aux. regression for ARCH test */
    AUX_DF,    /* aux. regression for Dickey-Fuller test */
    AUX_ADF,   /* aux. regression for augmented Dickey-Fuller test */
    AUX_OMIT,  /* aux. regression for omitting variables */
    AUX_RESET, /* aux. regression for Ramsey's RESET */
    AUX_SYS,   /* single equation from multivariate system */
    AUX_AUX    /* auxiliary regression not otherwise specified */
};

typedef struct {
    int m1;        /* ID for first model */
    int m2;        /* ID for second model */
    int ci;        /* estimator code for the first */
    int dfn;       /* numerator degrees of freedom */
    int dfd;       /* denominator degrees of freedom */ 
    double F;      /* F test statistic */
    double chisq;  /* Chi-square test statistic */
    double trsq;   /* T*R^2 test statistic */
    int score;     /* "cases correct" for discrete models */
    int robust;    /* = 1 when robust vcv is in use, else 0 */
} COMPARE;

/* functions follow */
 
int auxreg (LIST addvars, 
	    MODEL *orig, MODEL *new,
	    double ***pZ, DATAINFO *pdinfo, 
	    int aux_code, GRETLTEST *test, 
	    gretlopt opt, PRN *prn);

double robust_omit_F (const int *list, MODEL *pmod);

int omit_test (LIST omitvars, MODEL *orig, MODEL *new, 
	       double ***pZ, DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn);

int reset_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		PRN *prn, GRETLTEST *test);

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
		  PRN *prn);

int vif_test (MODEL *pmod, 
	      double ***pZ, DATAINFO *pdinfo, 
	      PRN *prn);

int leverage_test (MODEL *pmod, 
		   double ***pZ, DATAINFO *pdinfo, 
		   PRN *prn, PATHS *ppaths, 
		   gretlopt oflag);

int add_leverage_values_to_dataset (double ***pZ, DATAINFO *pdinfo,
				    gretl_matrix *m, unsigned char opt);

int mp_ols (const LIST list, const char *pos,
	    double ***pZ, DATAINFO *pdinfo, 
	    PRN *prn); 

int sum_test (LIST sumvars, MODEL *pmod, 
	      double ***pZ, DATAINFO *pdinfo, 
	      PRN *prn);

