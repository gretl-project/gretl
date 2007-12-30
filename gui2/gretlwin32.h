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
    TO_BACKSLASH,
    FROM_BACKSLASH
};

enum {
    WIN32_TO_CLIPBOARD,
    WIN32_TO_PRINTER
};

extern int use_wimp;

int create_child_process (char *prog);

void startR (char *Rcommand);

char *slash_convert (char *str, int which);

int unmangle (const char *dosname, char *longname);

void set_up_windows_look (void);

void menu_font_option_off (void);

void win_help (void);

void gretl_win32_init (const char *progname);

void gretl_win32_debug (void);

const char *get_network_cfg_filename (void);

int prn_to_clipboard (PRN *prn, int copycode);

int win_buf_to_clipboard (const char *buf);

int send_file (char *fullname);

void win32_process_graph (GPT_SPEC *spec, int color, int dest);

int browser_open (const char *url);

int win32_open_file (const char *fname);

void win32_raise_window (GtkWidget *w);

#endif /* GRETLWIN32_H */
