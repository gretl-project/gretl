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

/**
 * ModelAuxCode:
 * @AUX_NONE: not an auxiliary regression
 * @AUX_SQ: nonlinearity test (squared terms)
 * @AUX_LOG: nonlinearity test (log terms)
 * @AUX_CHOW: Chow test
 * @AUX_ADD: LM test regression for added variables
 * @AUX_AR: autocorrelation test
 * @AUX_ARCH: ARCH test
 * @AUX_WHITE: heteroskedasticity (White's test)
 * @AUX_COINT: cointegration test
 * @AUX_DF: Dickey-Fuller test
 * @AUX_ADF: augmented Dickey-Fuller test
 * @AUX_KPSS: KPSS unit-root test
 * @AUX_OMIT: unused
 * @AUX_RESET: Ramsey's RESET
 * @AUX_SYS: single equation from multivariate system
 * @AUX_VAR: single equation from VAR system
 * @AUX_VECM: single equation from VECM system
 * @AUX_JOHANSEN: Johansen cointegration test
 * @AUX_GROUPWISE: test for groupwise heteroskedasticity
 * @AUX_HET_1: Pesaran-Taylor HET_1 test
 * @AUX_BP: Breusch-Pagan heteroskedastcity test
 * @AUX_AUX: auxiliary regression not otherwise specified
 * @AUX_COMFAC: common factor test
 * @AUX_BIPROB: biprobit initializer
 *
 * Symbolic names to keep track of auxiliary regression models,
 * which are estimated either for the purpose of carrying out
 * some sort of diagnostic test or which form part of a
 * multi-equation system.
 */

typedef enum {
    AUX_NONE,
    AUX_SQ,
    AUX_LOG,
    AUX_CHOW,
    AUX_ADD,
    AUX_AR,
    AUX_ARCH,
    AUX_WHITE,
    AUX_COINT,
    AUX_DF,
    AUX_ADF,
    AUX_KPSS,
    AUX_OMIT,
    AUX_RESET,
    AUX_SYS,
    AUX_VAR,
    AUX_VECM,
    AUX_JOHANSEN,
    AUX_GROUPWISE,
    AUX_HET_1,
    AUX_BP,
    AUX_AUX,
    AUX_COMFAC,
    AUX_BIPROB,
} ModelAuxCode;

double wald_omit_F (const int *list, MODEL *pmod);

double wald_omit_chisq (const int *list, MODEL *pmod);

int add_test (MODEL *pmod, const int *addvars,
	      DATASET *dset, gretlopt opt,
	      PRN *prn);

int add_test_full (MODEL *orig, MODEL *pmod,
		   const int *addvars, DATASET *dset,
		   gretlopt opt, PRN *prn);

int omit_test (MODEL *pmod, const int *omitvars,
	       DATASET *dset, gretlopt opt,
	       PRN *prn);

int omit_test_full (MODEL *orig, MODEL *pmod,
		    const int *omitvars, DATASET *dset,
		    gretlopt opt, PRN *prn);

int nonlinearity_test (MODEL *pmod, DATASET *dset,
		       ModelAuxCode aux, gretlopt opt, PRN *prn);

int reset_test (MODEL *pmod, DATASET *dset,
		gretlopt opt, PRN *prn);

int autocorr_test (MODEL *pmod, int order, DATASET *dset,
		   gretlopt opt, PRN *prn);

int comfac_test (MODEL *pmod, DATASET *dset,
		 gretlopt opt, PRN *prn);

double get_DW_pvalue_for_model (MODEL *pmod,
				DATASET *dset,
				int *err);

int chow_test (int splitobs, MODEL *pmod, DATASET *dset,
	       gretlopt opt, PRN *prn);

int chow_test_from_dummy (int splitvar, MODEL *pmod, DATASET *dset,
			  gretlopt opt, PRN *prn);

int QLR_test (MODEL *pmod, DATASET *dset,
	      gretlopt opt, PRN *prn);

double QLR_pval (double X2, int df, double p1, double p2);

int cusum_test (MODEL *pmod, DATASET *dset,
		gretlopt opt, PRN *prn);

int panel_specification_test (MODEL *pmod, DATASET *dset,
			      gretlopt opt, PRN *prn);

int panel_hausman_test (MODEL *pmod, DATASET *dset,
			gretlopt opt, PRN *prn);

int vif_test (MODEL *pmod, DATASET *dset,
	      gretlopt opt, PRN *prn);

int bkw_test (MODEL *pmod, DATASET *dset,
	      gretlopt opt, PRN *prn);

int leverage_test (MODEL *pmod, DATASET *dset,
		   gretlopt opt, PRN *prn);

int add_leverage_values_to_dataset (DATASET *dset, gretl_matrix *m,
				    gretlopt opt, int flags);

void print_add_omit_null (const int *list, const DATASET *dset,
			  gretlopt opt, PRN *prn);

#endif /* COMPARE_H */
