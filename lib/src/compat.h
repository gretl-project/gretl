/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef COMPAT_H
#define COMPAT_H

void graphyzx (const int *list, const double *y1, const double *y2, 
	       const double *x, int n, const char *yname, 
	       const char *xname, const DATAINFO *pdinfo, 
	       gretlopt oflag, PRN *prn);

int ascii_graph (const int *list, const double **Z, const DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn);

int rhodiff (char *param, const int *list, double ***pZ, DATAINFO *pdinfo);

int gretl_multiply (char *s, int *list, char *sfx, double ***pZ,
		    DATAINFO *pdinfo);

#endif /* COMPAT_H */
