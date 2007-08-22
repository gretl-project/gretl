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

#ifndef MODELSPEC_H
#define MODELSPEC_H

#undef MSPEC_DEBUG

int model_sample_problem (const MODEL *pmod, const DATAINFO *pdinfo);

int modelspec_test_check (int test_ci, gretlopt opt, int model_id, 
			  DATAINFO *pdinfo, PRN *prn);

char *modelspec_get_command_by_id (int ID);

int modelspec_save (MODEL *pmod);

void free_modelspec (void);

#endif /* MODELSPEC_H */
