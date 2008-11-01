/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#ifndef BOXPLOTS_H
#define BOXPLOTS_H

int boxplots (int *list, double ***pZ, const DATAINFO *pdinfo, 
	      gretlopt opt);

int boolean_boxplots (const char *str, double ***pZ, 
		      DATAINFO *pdinfo, gretlopt opt);

int gnuplot_from_boxplot (const char *fname);

int boxplot_numerical_summary (const char *fname, PRN *prn);

const char *get_last_boxplots_string (void);

#endif /* BOXPLOTS_H */
