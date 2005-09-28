/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef GRETLWIN32_H
#define GRETLWIN32_H

enum {
    TO_BACKSLASH,
    FROM_BACKSLASH
};

int create_child_process (char *prog, char *env);

void startR (const char *Rcommand);

char *slash_convert (char *str, int which);

int unmangle (const char *dosname, char *longname);

void set_up_windows_look (void);

void win_help (void);

void gretl_win32_init (const char *progname);

const char *get_network_cfg_filename (void);

void win32_make_user_dirs (void);

int prn_to_clipboard (PRN *prn, int copycode);

int win_buf_to_clipboard (const char *buf);

int send_file (char *fullname);

char *desktop_path (void);

#endif /* GRETLWIN32_H */
