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

#ifndef GRETL_PATHS_H
#define GRETL_PATHS_H

typedef enum {
    CURRENT_DIR,
    DATA_SEARCH,
    SCRIPT_SEARCH,
    FUNCS_SEARCH,
    USER_SEARCH
} SearchLocation;

void set_string_table_written (void);

int gretl_string_table_written (void);

int gretl_path_prepend (char *file, const char *path);

int gretl_normalize_path (char *path);

void set_fopen_use_utf8 (void);

FILE *gretl_fopen (const char *filename, const char *mode);

gzFile gretl_gzopen (const char *filename, const char *mode);

int gretl_mkdir (const char *path);

int gretl_deltree (const char *path);

int gretl_write_access (char *fname);

int gretl_is_xml_file (const char *fname);

char *addpath (char *fname, PATHS *ppaths, int script);

int getopenfile (const char *line, char *fname, PATHS *ppaths,
		 gretlopt opt);

int gretl_path_is_absolute (const char *fname);

int has_system_prefix (const char *fname, const PATHS *ppaths,
		       int locus);

void show_paths (const PATHS *ppaths);

int gretl_set_paths (PATHS *ppaths, gretlopt opt);

const char *gretl_lib_path (void);

const char *gretl_dot_dir (void);

const char *gretl_work_dir (void);

char *gretl_default_workdir (PATHS *paths);

int set_gretl_work_dir (const char *path, PATHS *ppaths);

const char *gretl_maybe_switch_dir (const char *fname);

const char *gretl_gnuplot_path (void);

const char *gretl_plotfile (void);

char *set_gretl_plotfile (const char *fname);

const char *gretl_tramo (void);

const char *gretl_x12_arima (void);

const char *gretl_x12_arima_dir (void);

const char *gretl_png_font (void);

void set_gretl_png_font (const char *s, PATHS *ppaths);

#ifndef WIN32

int cli_read_rc (PATHS *paths);

#endif

#endif /* GRETL_PATHS_H */
