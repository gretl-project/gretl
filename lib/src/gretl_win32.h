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

#include <sys/types.h>
#include <dirent.h>
#include <windows.h>

enum {
    TO_BACKSLASH,
    FROM_BACKSLASH
};

DIR *win32_opendir (const char *dname);

int read_reg_val (HKEY tree, const char *base,
		  char *keyname, char *keyval);

int read_reg_val_with_fallback (HKEY tree0, HKEY tree1, 
				const char *base, char *keyname, 
				char *keyval);

int write_reg_val (HKEY tree, const char *base,
		   const char *keyname, const char *keyval);

void cli_read_registry (char *callname, PATHS *ppaths);

void win_show_last_error (void);

int winfork (char *cmdline, const char *dir, int wshow,
	     DWORD flags);

int gretl_spawn (char *cmdline);

int gretl_shell (const char *arg, PRN *prn);

char *desktop_path (void);

char *appdata_path (void);

char *mydocs_path (void);

int win32_write_access (char *path);

int win32_rename (const char *oldpath, const char *newpath);

int win32_delete_dir (const char *path);

char *R_path_from_registry (void);

#endif /* WIN32 */

#endif /* GRETL_WIN32_H */
