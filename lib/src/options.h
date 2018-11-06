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

typedef enum {
    OPT_NO_PARM = 0,
    OPT_ACCEPTS_PARM,
    OPT_NEEDS_PARM,
    OPT_AMBIGUOUS
} OptStatus;

gretlopt opt_from_flag (unsigned char c);

const char *print_flags (gretlopt oflags, int ci);

const char **get_opts_for_command (int ci, int *nopt);

int set_options_for_command (const char *cmdword,
			     const char *param,
			     gretlopt opt);

void maybe_get_stored_options (int ci, gretlopt *popt);

int check_for_loop_only_options (int ci, gretlopt opt, PRN *prn);

int cluster_option_ok (int ci);

int matrix_data_option (int ci, gretlopt opt);

char **get_all_option_strings (int *pn);

gretlopt transcribe_option_flags (gretlopt *targ, gretlopt src,
				  gretlopt test);

gretlopt delete_option_flags (gretlopt *targ, gretlopt test);

int incompatible_options (gretlopt opt, gretlopt test);

int options_incompatible_with (gretlopt opt, gretlopt base,
			       gretlopt test);

int option_prereq_missing (gretlopt opt, gretlopt test,
			   gretlopt prereq);

int inapplicable_option_error (int ci, gretlopt opt);

int push_option_param (int ci, gretlopt opt, char *val);

double get_optval_double (int ci, gretlopt opt, int *err);

int set_optval_double (int ci, gretlopt opt, double x);

const char *get_optval_string (int ci, gretlopt opt);

int set_optval_string (int ci, gretlopt opt, const char *s);

int set_optval_int (int ci, gretlopt opt, int k);

int get_optval_int (int ci, gretlopt opt, int *err);

int get_compression_option (int ci);

void destroy_option_params_at_level (int level);

void clear_stored_options_for_command (int ci);

void stored_options_cleanup (void);

void setopt_cleanup (void);

void option_printing_cleanup (void);

gretlopt valid_long_opt (int ci, const char *s, OptStatus *status);

gretlopt valid_short_opt (int ci, char c);

void debug_print_option_flags (const char *msg, gretlopt opt);

#endif /* OPTIONS_H */
