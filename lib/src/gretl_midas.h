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
    MIDAS_U,
    MIDAS_NEALMON,
    MIDAS_BETA0,
    MIDAS_BETAN,
    MIDAS_ALMONP,
    MIDAS_MAX /* sentinel */
} MidasType;

int midas_days_per_period (int days_per_week, int pd);

DATASET *midas_aux_dataset (const int *list,
			    const DATASET *dset,
			    int *err);

MODEL midas_model (const int *list, const char *param,
		   DATASET *dset, gretlopt opt,
		   PRN *prn);

#endif /* GRETL_MIDAS_H */
