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

#ifndef BOXPLOTS_H
#define BOXPLOTS_H

extern char boxplottmp[MAXLEN];

int boxplots (int *list, char **bools, 
	      double ***pZ, const DATAINFO *pdinfo, 
	      int notches);

int boolean_boxplots (const char *str, double ***pZ, 
		      DATAINFO *pdinfo, int notches);

int retrieve_boxplot (const char *fname);

int augment_boxplot_count (void);

void zero_boxplot_count (void);

#endif /* BOXPLOTS_H */


