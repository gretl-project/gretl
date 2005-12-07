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

int hp_filter (const double *x, double *hp, const DATAINFO *pdinfo);

int bkbp_filter (const double *y, double *bk, const DATAINFO *pdinfo);

int get_fracdiff (const double *y, double *diffvec, double d,
		  const DATAINFO *pdinfo);

int dummy (double ***pZ, DATAINFO *pdinfo, int center);

int panel_unit_dummies (double ***pZ, DATAINFO *pdinfo);

int panel_unit_first_obs (int t, const DATAINFO *pdinfo);

int paneldum (double ***pZ, DATAINFO *pdinfo);

int genrunit (double ***pZ, DATAINFO *pdinfo);

int genrtime (double ***pZ, DATAINFO *pdinfo, int tm);

int plotvar (double ***pZ, DATAINFO *pdinfo, const char *period);

#endif /* GENRFUNCS_H */
