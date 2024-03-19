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

#ifndef GRETL_MULTIPLOT_H_
#define GRETL_MULTIPLOT_H_

int gretl_gridplot_collecting (void);

int gretl_gridplot_start (const char *param, gretlopt opt,
			  DATASET *dset);

int gretl_gridplot_add_plot (gchar *buf);

int gretl_gridplot_finalize (gretlopt opt);

int gretl_gridplot_from_array (const char *param, gretlopt opt);

int check_gridplot_options (gretlopt opt);

void gretl_gridplot_destroy (void);

#endif /* GRETL_MULTIPLOT_H_ */
