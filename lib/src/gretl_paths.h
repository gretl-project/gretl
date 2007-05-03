/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2004 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef GRETL_PATHS_H
#define GRETL_PATHS_H

void set_string_table_written (void);

int gretl_string_table_written (void);

int gretl_path_prepend (char *file, const char *path);

FILE *gretl_fopen (const char *filename, const char *mode);

gzFile gretl_gzopen (const char *filename, const char *mode);

int gretl_is_xml_file (const char *fname);

char *addpath (char *fname, PATHS *ppaths, int script);

int getopenfile (const char *line, char *fname, PATHS *ppaths,
		 gretlopt opt);

int gretl_path_is_absolute (const char *fname);

void show_paths (const PATHS *ppaths);

int set_paths (PATHS *ppaths, gretlopt opt);

const char *gretl_lib_path (void);

const char *gretl_user_dir (void);

void set_gretl_user_dir (const char *path, PATHS *ppaths);

void gretl_maybe_switch_dir (const char *fname);

const char *gretl_gnuplot_path (void);

const char *gretl_plotfile (void);

char *set_gretl_plotfile (const char *fname);

const char *gretl_x12_arima (void);

const char *gretl_x12_arima_dir (void);

const char *gretl_png_font (void);

void set_gretl_png_font (const char *s, PATHS *ppaths);

#endif /* GRETL_PATHS_H */
