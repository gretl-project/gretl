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

#ifndef FORECAST_H
#define FORECAST_H

#ifdef  __cplusplus
extern "C" {
#endif

typedef enum {
    FC_AUTO_OK      = 1 << 0,
    FC_DYNAMIC_OK   = 1 << 1,
    FC_ADDOBS_OK    = 1 << 2,
    FC_INTEGRATE_OK = 1 << 3,
    FC_MEAN_OK      = 1 << 4
} FcastFlags;

struct FITRESID_ {
    int model_ID;   /* ID of model on which forecast is based */
    int asymp;      /* 0/1 flag for asymptotic estimator */
    int std;        /* 0/1 flag for standardized residuals */
    int model_t1;   /* start of model estimation range */
    int method;     /* one of the ForecastMethod options */
    double *actual; /* array of values of dependent variable */
    double *fitted; /* array of fitted values */
    double *resid;  /* array of residuals */
    double *sderr;  /* array of forecast standard errors (or NULL) */
    double sigma;   /* standard error of regression */
    double alpha;   /* for confidence intervals */
    int pmax;       /* if positive, suggested number of decimal places
                       for use in printing */
    int df;         /* degrees of freedom for model */
    int t0;         /* start of pre-forecast data range */
    int t1;         /* start of forecast range */
    int t2;         /* end of forecast range */
    int k;          /* number of steps ahead (method = FC_KSTEP only) */
    int nobs;       /* length of the arrays actual, fitted, resid */
    char depvar[VNAMELEN]; /* name of dependent variable */
};

void free_fit_resid (FITRESID *fr);

FITRESID *get_fit_resid (const MODEL *pmod, const DATASET *dset,
			 int *err);

FITRESID *get_forecast (MODEL *pmod, int t1, int t2, int pre_n,
			DATASET *dset, gretlopt opt, int *err);

FITRESID *get_system_forecast (void *p, int ci, int i,
			       int t1, int t2, int pre_n,
			       DATASET *dset, gretlopt opt,
			       int *err);

gretl_matrix *matrix_forecast (MODEL *pmod,
                               const gretl_matrix *X,
                               int *err);

int do_forecast (const char *str, DATASET *dset,
		 gretlopt opt, PRN *prn);

void forecast_options_for_model (MODEL *pmod, const DATASET *dset,
				 FcastFlags *flags, int *dt2max,
				 int *st2max);

gretl_matrix *get_forecast_matrix (int idx, int *err);

FITRESID *
recursive_OLS_k_step_fcast (MODEL *pmod, DATASET *dset,
			    int t1, int t2, int k, int pre_n,
			    int *err);

void fcast_get_continuous_range (const FITRESID *fr, int *pt1, int *pt2);

void forecast_matrix_cleanup (void);

#ifdef  __cplusplus
}
#endif

#endif /* FORECAST_H */
