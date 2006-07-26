/*
 *  Copyright (c) by Allin Cottrell and Riccardo Lucchetti
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef GENRFUNCS_H
#define GENRFUNCS_H

/* private functions, used only in generate.c */

int get_fracdiff (const double *y, double *diffvec, double d,
		  const DATAINFO *pdinfo);

int genrunit (double ***pZ, DATAINFO *pdinfo);

double genr_cov_corr (const char *s, double ***pZ, 
		      const DATAINFO *pdinfo, int fn);

double genr_vcv (const char *s, const DATAINFO *pdinfo, MODEL *pmod);

int get_observation_number (const char *s, const DATAINFO *pdinfo);

double 
get_model_data_element (MODEL *pmod, int idx, const char *s,
			const DATAINFO *pdinfo, int *err);

#endif /* GENRFUNCS_H */
