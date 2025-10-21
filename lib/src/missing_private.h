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

#ifndef MISSING_PRIVATE_H
#define MISSING_PRIVATE_H

#define masked(m,t) (m != NULL && m[t] == '1')

void set_reference_missmask_from_model (const MODEL *pmod);

int copy_to_reference_missmask (const char *mask);

int apply_reference_missmask (MODEL *pmod);

int reference_missmask_present (void);

int model_adjust_sample (MODEL *pmod, const DATASET *dset,
                         int *misst);

char *model_sample_mask (const MODEL *pmod);

#endif /* MISSING_PRIVATE_H */
