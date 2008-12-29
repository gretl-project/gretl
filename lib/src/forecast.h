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

typedef enum {
    FC_STATIC,
    FC_DYNAMIC,
    FC_AUTO,
    FC_KSTEP
} ForecastMethod;

typedef enum {
    FC_AUTO_OK    = 1 << 0,
    FC_DYNAMIC_OK = 1 << 1,
    FC_ADDOBS_OK  = 1 << 2
} ForecastFlags;

struct FITRESID_ {
    int model_ID;
    int asymp;
    int model_t1;
    int method;
    double *actual;
    double *fitted;
    double *resid;
    double *sderr;
    double sigma;
    double tval;
    int pmax;
    int df;
    int t0, t1, t2;
    int k; 
    int nobs;
    char depvar[VNAMELEN];
};

void free_fit_resid (FITRESID *fr);

FITRESID *get_fit_resid (const MODEL *pmod, const double **Z, 
			 const DATAINFO *pdinfo, int *err);

FITRESID *get_forecast (MODEL *pmod, int t1, int t2, int pre_n,
			double ***pZ, DATAINFO *pdinfo,
			gretlopt opt, int *err);

FITRESID *get_system_forecast (void *p, int ci, int i, 
			       int t1, int t2, int pre_n,
			       const double **Z, DATAINFO *pdinfo,
			       gretlopt opt, int *err);

int do_forecast (const char *str, double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn);

void forecast_options_for_model (MODEL *pmod, const double **Z,
				 const DATAINFO *pdinfo, int *flags,
				 int *dt2max, int *st2max);

gretl_matrix *get_forecast_matrix (int idx, int *err);

FITRESID *
rolling_OLS_k_step_fcast (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
			  int t1, int t2, int k, int pre_n, int *err);

void forecast_matrix_cleanup (void);

#endif /* FORECAST_H */


