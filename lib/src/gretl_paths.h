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

#include <sys/stat.h>
#include <dirent.h>

typedef enum {
    DATA_SEARCH,
    DB_SEARCH,
    FUNCS_SEARCH,
    SCRIPT_SEARCH,
    USER_SEARCH
} SearchType;

typedef enum {
    GRETL_HELPFILE,
    GRETL_CMD_HELPFILE,
    GRETL_CLI_HELPFILE
} HelpPaths;

typedef enum {
    PKG_ALL,
    PKG_SUBDIR,
    PKG_TOPLEV
} PkgType;

typedef struct ConfigPaths_ ConfigPaths;

/* these are all the gretl paths which are recorded in
   the gretl config file or Windows registry entries 
*/

struct ConfigPaths_ {
    char gretldir[MAXLEN];
    char workdir[MAXLEN];
#ifndef WIN32
    char gnuplot[MAXLEN];
#endif
    char x12a[MAXLEN];
    char tramo[MAXLEN];
    char rbinpath[MAXLEN];
    char rlibpath[MAXLEN];
    char oxlpath[MAXLEN];
    char octpath[MAXLEN];
    char statapath[MAXLEN];
    char pypath[MAXLEN];
    char mpiexec[MAXLEN];
    char mpi_hosts[MAXLEN];
    char dbhost[64];
    char pngfont[128];
};

void set_string_table_written (void);

int gretl_string_table_written (void);

int gretl_path_prepend (char *file, const char *path);

int gretl_normalize_path (char *path);

int slash_terminate (char *path);

void set_stdio_use_utf8 (void);

int get_stdio_use_utf8 (void);

int string_is_utf8 (const unsigned char *s);

FILE *gretl_fopen (const char *fname, const char *mode);

FILE *gretl_fopen_with_recode (const char *fname, const char *mode,
			       char **recoded_fname);

int gretl_test_fopen (const char *fname, const char *mode);

FILE *gretl_read_user_file (const char *fname);

FILE *gretl_mktemp (char *pattern, const char *mode);

int gretl_open (const char *pathname, int flags);

int gretl_rename (const char *oldpath, const char *newpath);

int gretl_remove (const char *path);

gzFile gretl_gzopen (const char *fname, const char *mode);

int gretl_stat (const char *fname, struct stat *buf);

int gretl_mkdir (const char *path);

int gretl_chdir (const char *path);

DIR *gretl_opendir (const char *name);

int gretl_deltree (const char *path);

int gretl_setenv (const char *name, const char *value);

int gretl_write_access (char *fname);

int gretl_is_xml_file (const char *fname);

int gretl_isdir (const char *path);

char *gretl_addpath (char *fname, int script);

int getopenfile (const char *line, char *fname, gretlopt opt);

int fname_has_path (const char *fname);

int has_system_prefix (const char *fname, SearchType stype);

void show_paths (void);

int gretl_set_paths (ConfigPaths *paths);

int gretl_update_paths (ConfigPaths *cpaths, gretlopt opt);

char **get_plausible_search_dirs (SearchType stype, int *n_dirs);

char *gretl_function_package_get_path (const char *name,
				       PkgType type);

void set_gretl_plugin_path (const char *path);

const char *helpfile_path (int id);

const char *gretl_home (void);

const char *gretl_lib_path (void);

const char *gretl_dotdir (void);

const char *gretl_workdir (void);

const char *gretl_default_workdir (void);

const char *maybe_get_default_workdir (void);

char *gretl_make_dotpath (const char *basename);

int set_gretl_work_dir (const char *path);

const char *gretl_maybe_switch_dir (const char *fname);

char *gretl_maybe_prepend_dir (char *fname);

const char *gretl_gnuplot_path (void);

const char *gretl_plotfile (void);

char *set_gretl_plotfile (const char *fname);

void report_plot_written (PRN *prn);

const char *gretl_binbase (void);

const char *gretl_tramo (void);

const char *gretl_tramo_dir (void);

const char *gretl_x12_arima (void);

const char *gretl_x12_arima_dir (void);

int gretl_x12_is_x13 (void);

const char *gretl_rbin_path (void);

const char *gretl_rlib_path (void);

const char *gretl_png_font (void);

const char *gretl_oxl_path (void);

const char *gretl_octave_path (void);

const char *gretl_stata_path (void);

const char *gretl_python_path (void);

const char *gretl_mpi_hosts (void);

const char *gretl_mpiexec (void);

const char *gretl_current_dir (void);

void gretl_set_current_dir (const char *s);

void set_gretl_png_font (const char *s);

void get_gretl_config_from_file (FILE *fp, ConfigPaths *cpaths,
				 char *dbproxy, int *use_proxy);

#ifdef WIN32

void win32_set_gretldir (const char *progname);

#else

void get_gretl_rc_path (char *rcfile);

int cli_read_rc (void);

#endif

#ifdef OS_OSX

const char *gretl_app_support_dir (void);

#endif

#endif /* GRETL_PATHS_H */
