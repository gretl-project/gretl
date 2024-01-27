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
    GRETL_CMDREF,
    GRETL_FUNCREF,
    GRETL_GUI_HELP,
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
#if !defined(WIN32) || !defined(PKGBUILD)
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
    char jlpath[MAXLEN];
    char lppath[MAXLEN];
    char mpiexec[MAXLEN];
    char mpi_hosts[MAXLEN];
    char pngfont[128];
    int no_dotdir;
};

int gretl_path_prepend (char *file, const char *path);

int gretl_normalize_path (char *path);

int slash_terminate (char *path);

int utf8_encoded (const char *s);

FILE *gretl_fopen (const char *fname, const char *mode);

int gretl_test_fopen (const char *fname, const char *mode);

FILE *gretl_read_user_file (const char *fname);

FILE *gretl_mktemp (char *pattern, const char *mode);

int gretl_open (const char *pathname, int flags, int mode);

int gretl_rename (const char *oldpath, const char *newpath);

int gretl_remove (const char *path);

gzFile gretl_gzopen (const char *fname, const char *mode);

int gretl_stat (const char *fname, struct stat *buf);

int gretl_file_exists (const char *fname);

int gretl_mkdir (const char *path);

int gretl_chdir (const char *path);

GDir *gretl_opendir (const char *name);

int gretl_deltree (const char *path);

int gretl_setenv (const char *name, const char *value);

int gretl_write_access (char *fname);

int gretl_is_xml_file (const char *fname);

int gretl_isdir (const char *path);

char *gretl_addpath (char *fname, int script);

int get_full_filename (const char *fname, char *fullname, 
		       gretlopt opt);

int fname_has_path (const char *fname);

int has_system_prefix (const char *fname, SearchType stype);

void show_paths (void);

int gretl_set_paths (ConfigPaths *paths);

int gretl_update_paths (ConfigPaths *cpaths, gretlopt opt);

int gretl_set_path_by_name (const char *name, const char *path);

char **get_plausible_search_dirs (SearchType stype, int *n_dirs);

char *gretl_function_package_get_path (const char *name,
				       PkgType type);

int get_package_data_path (int ci, const char *fname, char *fullname);

void set_gretl_plugin_path (const char *path);

const char *helpfile_path (int id, int cli, int en);

int using_translated_helpfile (int id);

const char *gretl_home (void);

const char *gretl_bindir (void);

const char *gretl_plugin_path (void);

const char *gretl_dotdir (void);

const char *gretl_workdir (void);

const char *maybe_get_default_workdir (void);

gchar *gretl_make_dotpath (const char *basename);

const char *gretl_maybe_switch_dir (const char *fname);

char *gretl_maybe_prepend_dir (char *fname);

const char *gretl_gnuplot_path (void);

const char *gretl_plotfile (void);

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

const char *gretl_julia_path (void);

const char *gretl_lpsolve_path (void);

const char *gretl_mpi_hosts (void);

const char *gretl_mpiexec (void);

const char *gretl_function_package_path (void);

void gretl_set_script_dir (const char *s);

void gretl_script_dirs_cleanup (void);

char *gretl_prepend_homedir (const char *fname, int *err);

void set_gretl_png_font (const char *s);

void get_gretl_config_from_file (FILE *fp, ConfigPaths *cpaths,
				 char *dbproxy, int *use_proxy,
				 int *updated, gchar **gptheme);

int gretl_path_compose (char *targ, int len,
			const char *s1,
			const char *s2);

char *gretl_build_path (char *targ,
			const gchar *first_element,
			...);

gretl_bundle *foreign_info (void);

gchar *get_download_path (const char *dlname, int *err);

#ifdef WIN32

void win32_set_gretldir (void);

#else

void get_gretl_rc_path (char *rcfile);

int cli_read_rc (void);

#endif

#ifdef OS_OSX

const char *gretl_app_support_dir (void);

#endif

#endif /* GRETL_PATHS_H */
