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

int boxplots (int *list, char **bools, 
	      double ***pZ, const DATAINFO *pdinfo, 
	      gretlopt opt);

int boolean_boxplots (const char *str, double ***pZ, 
		      DATAINFO *pdinfo, gretlopt opt);

int retrieve_boxplot (const char *fname);

int ps_print_plots (const char *fname, int flag, gpointer data);

#ifndef G_OS_WIN32
int plot_to_xpm (const char *fname, gpointer data);
#endif

const char *get_boxplots_string (void);

const char *get_boxdump_name (void);

#endif /* BOXPLOTS_H */


