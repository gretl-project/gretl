/*
 *   Copyright (c) by Allin Cottrell
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

#ifndef GRETL_PANEL_H_
#define GRETL_PANEL_H_

int panel_diagnostics (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn);

MODEL real_panel_model (const int *list, double ***pZ, DATAINFO *pdinfo,
			gretlopt opt, PRN *prn);

MODEL panel_wls_by_unit (const int *list, double ***pZ, DATAINFO *pdinfo,
			 gretlopt opt, PRN *prn);

int panel_autocorr_test (MODEL *pmod, int order, 
			 double **Z, DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn);

int set_panel_structure_from_vars (int uv, int tv, 
				   double **Z, 
				   DATAINFO *pdinfo);

int set_panel_structure_from_line (const char *line, 
				   double **Z, 
				   DATAINFO *pdinfo);

int switch_panel_orientation (double **Z, DATAINFO *pdinfo);

int balanced_panel (const DATAINFO *pdinfo);

int *panel_list_omit (const MODEL *orig, const int *drop, int *err);

int *panel_list_add (const MODEL *orig, const int *add, int *err);

int panel_variance_info (const double *x, const DATAINFO *pdinfo,
			 double xbar, double *psw, double *psb);

int panel_obs_info (const int *list, const double **Z, const DATAINFO *pdinfo,
		    PRN *prn);

#endif /* GRETL_PANEL_H_ */
