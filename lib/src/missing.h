/*
 *  Copyright (c) by Allin Cottrell
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

#ifndef MISSING_H
#define MISSING_H

#include <float.h>

#define NADBL DBL_MAX
#define na(x) (x == NADBL)

#define missing_masked(m,t,t1) (m != NULL && m[t-t1] != 0)
#define model_missing(m,t)     ((m)->missmask != NULL && \
                                (m)->missmask[t - (m)->t1] != 0)
#define has_missing_obs(m)     ((m)->missmask != NULL)

int model_missval_count (const MODEL *pmod);

void set_miss (LIST list, const char *param, double **Z,
	       DATAINFO *pdinfo, PRN *prn);

#endif /* MISSING_H */
