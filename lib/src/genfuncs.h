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

#ifndef GENFUNCS_H
#define GENFUNCS_H

int sort_series (const double *x, double *y, int f, 
		 const DATAINFO *pdinfo);

int gretl_sort_by (const double *x, const double *y, 
		   double *z, const DATAINFO *pdinfo);

int diff_series (const double *x, double *y, int f, 
		 const DATAINFO *pdinfo);

int orthdev_series (const double *x, double *y, const DATAINFO *pdinfo);

int cum_series (const double *x, double *y, 
		const DATAINFO *pdinfo);

int resample_series (const double *x, double *y, 
		     const DATAINFO *pdinfo);

int fracdiff_series (const double *x, double *y, double d,
		     const DATAINFO *pdinfo);

int panel_mean_series (const double *x, double *y, const DATAINFO *pdinfo);

int panel_sd_series (const double *x, double *y, const DATAINFO *pdinfo);

int hp_filter (const double *x, double *hp, const DATAINFO *pdinfo,
	       gretlopt opt);

int bkbp_filter (const double *y, double *bk, 
		 const DATAINFO *pdinfo);

int dummy (double ***pZ, DATAINFO *pdinfo, int center);

int panel_dummies (double ***pZ, DATAINFO *pdinfo, gretlopt opt);

int gen_unit (double ***pZ, DATAINFO *pdinfo);

int panel_unit_first_obs (int t, const DATAINFO *pdinfo);

int gen_time (double ***pZ, DATAINFO *pdinfo, int tm);

int gen_wkday (double ***pZ, DATAINFO *pdinfo);

int plotvar_code (const DATAINFO *pdinfo);

const double *gretl_plotx (const DATAINFO *pdinfo);

int genr_fit_resid (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		    int code, int undo);

int get_observation_number (const char *s, const DATAINFO *pdinfo);

int get_t_from_obs_string (const char *s, const double **Z, 
			   const DATAINFO *pdinfo);

#endif /* GENFUNCS_H */
