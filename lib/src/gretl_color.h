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

#ifndef GRETL_COLOR_H
#define GRETL_COLOR_H

gretl_array *colormix_array (gretlRGB c1, gretlRGB c2,
			     const double *f, int nf,
			     int do_plot, int *err);

void decompose_argb (guint32 u,
		     guint8 *a,
		     guint8 *r,
		     guint8 *g,
		     guint8 *b);

#endif /* GRETL_COLOR_H */
