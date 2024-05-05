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

#ifndef GRETL_WIN32_H
#define GRETL_WIN32_H

#ifdef WIN32

#include <windows.h>

enum {
    TO_BACKSLASH,
    FROM_BACKSLASH
};

enum {
    RGUI,
    RTERM,
    REXE,
    RLIB,
    RBASE
};

void set_windebug (int s);

int read_reg_val (HKEY tree, const char *base,
		  char *keyname, char *keyval);

const char *get_gretlnet_filename (void);

int set_gretlnet_filename (const char *prog);

void win32_cli_read_rc (void);

void win_show_last_error (void);

int win32_ensure_path (void);

int win_run_sync (const char *cmdline, const char *currdir);

int win_run_async (const char *cmdline, const char *currdir);

int gretl_spawn (const char *cmdline);

int gnuplot_make_image (const char *input_fname);

int gretl_shell (const char *arg, gretlopt opt, PRN *prn);

int gretl_win32_pipe_output (const char *cmdline,
			     const char *currdir,
			     PRN *prn);

int gretl_win32_grab_stdout (const char *cmdline,
			     const char *currdir,
			     char **sout);

char *slash_convert (char *str, int which);

char *desktop_path (void);

char *appdata_path (void);

char *mydocs_path (void);

char *program_files_path (void);

char *program_files_x86_path (void);

int win32_write_access (const char *path);

int win32_remove (const char *path);

int R_home_from_registry (char *s);

int win32_R_path (char *s, int which);

int win32_check_for_program (const char *prog);

char *strptime (const char *buf, const char *format,
		struct tm *timeptr);

double win32_mktime (struct tm *tm);

double win32_fscan_nonfinite (FILE *fp, int *err);

double win32_sscan_nonfinite (const char *s, int *err);

void win32_pprint_nonfinite (PRN *prn, double x, char c);

double win32_get_time (void);

int try_for_CP_65001 (void);

int win32_get_core_count (void);

int win32_get_stack_size (void);

#endif /* WIN32 */

#endif /* GRETL_WIN32_H */
