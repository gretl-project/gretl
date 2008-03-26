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

/* compare.h for gretl */

#ifndef COMPARE_H
#define COMPARE_H

#include <stdio.h>
#include "gretl_matrix.h"

typedef enum {
    SAVE_LEVERAGE  = 1 << 0,
    SAVE_INFLUENCE = 1 << 1,
    SAVE_DFFITS    = 1 << 2,
    SAVE_TREND     = 1 << 3,
    SAVE_CYCLE     = 1 << 4
} SeriesSaveCode;

typedef enum {
    AUX_NONE,  /* not an auxiliary regression */
    AUX_SQ,    /* aux. regression for nonlinearity (squared terms) */
    AUX_LOG,   /* aux. regression for nonlinearity (log terms) */
    AUX_CHOW,  /* aux. regression for Chow test */
    AUX_ADD,   /* aux. regression for adding variables */
    AUX_AR,    /* aux. regression for autocorrelation test */
    AUX_ARCH,  /* aux. regression for ARCH test */
    AUX_WHITE, /* aux. regression for heteroskedasticity (White's test) */
    AUX_COINT, /* aux. regression for cointegreation test */
    AUX_DF,    /* aux. regression for Dickey-Fuller test */
    AUX_ADF,   /* aux. regression for augmented Dickey-Fuller test */
    AUX_KPSS,  /* aux. regression for KPSS test */
    AUX_OMIT,  /* aux. regression for omitting variables */
    AUX_RESET, /* aux. regression for Ramsey's RESET */
    AUX_SYS,   /* single equation from multivariate system */
    AUX_VAR,   /* single equation from VAR system */
    AUX_VECM,  /* single equation from VECM system */
    AUX_JOHANSEN,  /* Johansen cointegration test */
    AUX_GROUPWISE, /* testing for groupwise heteroskedasticity */
    AUX_HET_1, /* aux. regression for Pesaran-Taylor HET_1 test */
    AUX_BP,    /* aux. regression for Breusch-Pagan heterosked. test */
    AUX_AUX    /* auxiliary regression not otherwise specified */
} ModelAuxCode;

/* functions follow */
 
double robust_omit_F (const int *list, MODEL *pmod);

int add_test (const int *addvars,  MODEL *orig, MODEL *new,
	      double ***pZ, DATAINFO *pdinfo, 
	      gretlopt opt, PRN *prn);

int omit_test (const int *omitvars, MODEL *orig, MODEL *new, 
	       double ***pZ, DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn);

int nonlinearity_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		       int aux_code, gretlopt opt, PRN *prn); 

int reset_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn);

int autocorr_test (MODEL *pmod, int order,
		   double ***pZ, DATAINFO *pdinfo, 
		   gretlopt opt, PRN *prn);

double ljung_box (int m, int t1, int t2, const double *y, int *err);

int chow_test (const char *line, MODEL *pmod, 
	       double ***pZ, DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn);

int cusum_test (MODEL *pmod, 
		double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn);

int panel_hausman_test (MODEL *pmod, 
			double ***pZ, DATAINFO *pdinfo, 
			gretlopt opt, PRN *prn);

int vif_test (MODEL *pmod, 
	      double ***pZ, DATAINFO *pdinfo, 
	      PRN *prn);

int leverage_test (MODEL *pmod, 
		   double ***pZ, DATAINFO *pdinfo, 
		   gretlopt opt, PRN *prn);

int add_leverage_values_to_dataset (double ***pZ, DATAINFO *pdinfo,
				    gretl_matrix *m, unsigned char flags);

int lmtest_driver (const char *param,
		   double ***pZ, DATAINFO *pdinfo, 
		   gretlopt opt, PRN *prn);

#endif /* COMPARE_H */
