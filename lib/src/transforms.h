/*
 *  Copyright (c) 2004 by Allin Cottrell
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

#ifndef TRANSFORMS_H
#define TRANSFORMS_H

int rhodiff (char *param, const LIST list, 
	     double ***pZ, DATAINFO *pdinfo);

int diffgenr (int v, double ***pZ, DATAINFO *pdinfo, int ldiff);

int laggenr (int v, int lag, double ***pZ, DATAINFO *pdinfo);

int loggenr (int v, double ***pZ, DATAINFO *pdinfo);

int xpxgenr (int vi, int vj, double ***pZ, DATAINFO *pdinfo);

int list_diffgenr (const LIST list, double ***pZ, DATAINFO *pdinfo);

int list_ldiffgenr (const LIST list, double ***pZ, DATAINFO *pdinfo);

int list_laggenr (const LIST list, double ***pZ, DATAINFO *pdinfo);

int list_loggenr (const LIST list, double ***pZ, DATAINFO *pdinfo);

int list_xpxgenr (const LIST list, double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt);

#endif /* TRANSFORMS_H */
