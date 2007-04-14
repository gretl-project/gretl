/*
 *  Copyright (c) 2004 by Allin Cottrell
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

#ifndef OPTIONS_H
#define OPTIONS_H

gretlopt get_gretl_options (char *line, int *err);

const char *print_flags (gretlopt oflags, int ci);

const char **get_opts_for_command (int ci, int *nopt);

int check_for_loop_only_options (int ci, gretlopt opt, PRN *prn);

char **get_all_option_strings (int *pn);

int incompatible_options (gretlopt opt, gretlopt test);

#endif /* OPTIONS_H */
