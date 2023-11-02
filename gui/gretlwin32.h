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

#ifndef GRETLWIN32_H
#define GRETLWIN32_H

#include <sys/types.h>
#include <dirent.h>

#include "plotspec.h"

enum {
    WIN32_TO_CLIPBOARD,
    WIN32_TO_PRINTER
};

int win32_run_async (const char *prog, const char *arg);

void win32_start_R_async (void);

void gretl_win32_debug_init (int debug);

void record_win32_locale (char *s);

int win32_lang_id_from_locale (void);

void gretl_win32_init (int debug, int ignore_rc);

int prn_to_clipboard (PRN *prn, int fmt);

int buf_to_clipboard (const char *buf);

int emf_to_clipboard (char *emfname);

int browser_open (const char *url);

int win32_open_file (const char *fname);

int win32_open_pdf (const char *fname, const char *dest);

void win32_font_selector (char *fontname, int flag);

int win32_rename_dir (const char *oldname, const char *newname);

void get_default_windows_app_font (char *target);

#endif /* GRETLWIN32_H */
