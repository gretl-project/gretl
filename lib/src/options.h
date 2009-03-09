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

#ifndef OPTIONS_H
#define OPTIONS_H

gretlopt get_gretl_options (char *line, int *err);

gretlopt opt_from_flag (unsigned char c);

const char *print_flags (gretlopt oflags, int ci);

const char **get_opts_for_command (int ci, int *nopt);

int check_for_loop_only_options (int ci, gretlopt opt, PRN *prn);

char **get_all_option_strings (int *pn);

gretlopt transcribe_option_flags (gretlopt *targ, gretlopt src,
				  gretlopt test);

gretlopt delete_option_flags (gretlopt *targ, gretlopt test);

int incompatible_options (gretlopt opt, gretlopt test);

int option_prereq_missing (gretlopt opt, gretlopt test,
			   gretlopt prereq);

int inapplicable_option_error (int ci, gretlopt opt);

double get_optval_double (int ci, gretlopt opt);

void set_optval_double (int ci, gretlopt opt, double x);

const char *get_optval_string (int ci, gretlopt opt);

void clear_option_params (void);

void option_flags_cleanup (void);

#endif /* OPTIONS_H */
