/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2004 Ramu Ramanathan and Allin Cottrell
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
