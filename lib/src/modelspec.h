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

typedef struct MODELSPEC_ MODELSPEC;

int model_ci_from_modelspec (const MODELSPEC *spec, int i);

int model_sample_problem (const MODEL *pmod, const DATAINFO *pdinfo);

int modelspec_sample_problem (MODELSPEC *spec, int i, 
			      const DATAINFO *pdinfo);

int modelspec_last_index (const MODELSPEC *spec);

int modelspec_index_from_model_id (const MODELSPEC *spec, int ID);

char *modelspec_get_command_by_id (MODELSPEC *spec, int ID);

int modelspec_save (MODEL *pmod, MODELSPEC **pmspec);

void free_modelspec (MODELSPEC *spec);


#endif /* MODELSPEC_H */
