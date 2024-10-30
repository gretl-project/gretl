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

#ifndef GRETL_MIDAS_H
#define GRETL_MIDAS_H

typedef enum {
    MIDAS_U,        /* "unrestricted MIDAS" */
    MIDAS_NEALMON,  /* normalized exponential Almon */
    MIDAS_BETA0,    /* normalized beta, last lag 0 */
    MIDAS_BETAN,    /* normalized beta, non-zero last lag */
    MIDAS_ALMONP,   /* plain Almon polynomial */
    MIDAS_BETA1,    /* as MIDAS_BETA0, but 1st coeff clamped */
    MIDAS_MAX       /* sentinel */
} MidasType;

int midas_days_per_period (int days_per_week, int pd);

int is_valid_midas_frequency_ratio (const DATASET *dset, int m);

int get_midas_frequency (const DATASET *dset, int m);

int midas_m_from_pd (const DATASET *dset, int freq);

const char *midas_pdstr (const DATASET *dset, int cfac);

DATASET *midas_aux_dataset (const int *list,
			    const DATASET *dset,
			    int *err);

gretl_matrix *midas_list_to_vector (const int *list,
				    const DATASET *dset,
				    int *err);

int midas_forecast_setup (const MODEL *pmod,
			  DATASET *dset,
			  ForecastMethod method,
			  char **pformula);

MODEL midas_model (const int *list, const char *param,
		   DATASET *dset, gretlopt opt,
		   PRN *prn);

char *get_midas_term_line (const MODEL *pmod, int i);

#endif /* GRETL_MIDAS_H */
