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

enum paths_status {
    GRETL_USING_GUI      = 1 << 0,
    STRING_TABLE_WRITTEN = 1 << 1
};

void set_string_table_written (PATHS *ppaths);

int gretl_string_table_written (PATHS *ppaths);

int gretl_using_gui (const PATHS *ppaths);

char *addpath (char *fname, PATHS *ppaths, int script);

int getopenfile (const char *line, char *fname, PATHS *ppaths,
		 int setpath, int script);

void show_paths (PATHS *ppaths);

int set_paths (PATHS *ppaths, int defaults, int gui);

const char *fetch_gretl_lib_path (void);

#endif /* GRETL_PATHS_H */
