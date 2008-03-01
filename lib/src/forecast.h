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
    FC_ONESTEP
} ForecastMethod;

struct FITRESID_ {
    int model_ID;
    int model_ci;
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
    int nobs;
    int err;
    char depvar[VNAMELEN];
};

void free_fit_resid (FITRESID *fr);

FITRESID *get_fit_resid (const MODEL *pmod, const double **Z, 
			 const DATAINFO *pdinfo);

FITRESID *get_forecast (MODEL *pmod, int t0, int t1, int t2,
			double ***pZ, DATAINFO *pdinfo,
			gretlopt opt);

FITRESID *get_system_forecast (void *p, int ci, int i, 
			       int t0, int t1, int t2,
			       const double **Z, DATAINFO *pdinfo,
			       gretlopt opt);

int display_forecast (const char *str, double ***pZ, DATAINFO *pdinfo, 
		      gretlopt opt, PRN *prn);

int add_forecast (const char *str, MODEL *pmod, 
		  double ***pZ, DATAINFO *pdinfo,
		  gretlopt opt);

void forecast_options_for_model (MODEL *pmod, const double **Z,
				 const DATAINFO *pdinfo,
				 int *dyn_ok, int *add_obs_ok,
				 int *dt2max, int *st2max);

FITRESID *
rolling_OLS_one_step_fcast (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
			    int t0, int t1, int t2);

#endif /* FORECAST_H */


