/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
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

/* subsample.h for gretl */

#ifndef SUBSAMPLE_H
#define SUBSAMPLE_H

/* functions follow */

int attach_subsample_to_model (MODEL *pmod, 
			       double ***fullZ, const DATAINFO *fulldinfo);

int model_sample_issue (const MODEL *pmod, MODELSPEC *spec,
			double **Z, const DATAINFO *pdinfo);

int set_sample_dummy (const char *line, 
		      double ***oldZ, double ***newZ,
		      DATAINFO *oldinfo, DATAINFO *newinfo,
		      const int opt);

int set_sample (const char *line, DATAINFO *pdinfo);

int restore_full_sample (double ***subZ, double ***fullZ, double ***Z,
			 DATAINFO **subinfo, DATAINFO **fullinfo,
			 DATAINFO **datainfo); 

int count_missing_values (double ***pZ, DATAINFO *pdinfo, PRN *prn);

int allocate_case_markers (char ***S, int n);

#endif /* SUBSAMPLE_H */
