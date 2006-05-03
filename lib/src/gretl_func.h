/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2004 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef GRETL_FUNC_H
#define GRETL_FUNC_H

int n_user_functions (void);

const char *user_function_name_by_index (int i);

int user_function_index_by_name (const char *name);

int gretl_compiling_function (void);

int gretl_executing_function (void);

int gretl_start_compiling_function (const char *line, PRN *prn);

int gretl_function_append_line (const char *line);

int gretl_is_user_function (const char *line);

int gretl_is_public_user_function (const char *name);

int gretl_get_user_function (const char *line, char **fnname);

int is_user_matrix_function (const char *word);

int gretl_function_start_exec (const char *line, const char *fname,
			       double ***pZ, DATAINFO *pdinfo);

char *gretl_function_get_line (char *line, int len,
			       double ***pZ, DATAINFO **ppdinfo,
			       int *err);

int gretl_function_stack_depth (void);

void gretl_function_stop_on_error (double ***pZ, DATAINFO **ppdinfo, PRN *prn);

int gretl_function_flagged_error (const char *s, PRN *prn);

int gretl_function_set_info (int i, const char *help);

int gretl_function_get_info (int i, const char *key, char const **value);

void gretl_function_set_private (int i, int priv);

int write_selected_user_functions (const int *privlist, 
				   const int *publist, 
				   const char *author,
				   const char *version,
				   const char *date,
				   const char *descrip,
				   const char *fname);

int function_package_get_info (const char *fname,
			       int **privlist, 
			       int **publist,
			       char **author,
			       char **version,
			       char **date,
			       char **descrip);

int write_user_function_file (const char *fname);

int user_function_file_is_loaded (const char *fname);

int load_user_function_file (const char *fname);

int get_function_file_info (const char *fname, PRN *prn, char **pname);

int get_function_file_code (const char *fname, PRN *prn, char **pname);

char *get_function_file_header (const char *fname, int *err);

int user_function_help (const char *fnname, PRN *prn);

void gretl_functions_cleanup (void);

#endif /* GRETL_FUNC_H */
