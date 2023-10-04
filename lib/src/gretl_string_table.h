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

#ifndef GRETL_STRING_TABLE_H
#define GRETL_STRING_TABLE_H

gretl_string_table *gretl_string_table_new (const int *list);

int gretl_string_table_index (gretl_string_table *gst, const char *s, 
			      int col, int addcol, PRN *prn);

int gretl_string_table_finalize (gretl_string_table *gst, DATASET *dset);

int gretl_string_table_validate (gretl_string_table *gst,
				 gretlopt opt);

int gretl_string_table_save (gretl_string_table *gst, DATASET *dset);

void gretl_string_table_destroy (gretl_string_table *gst);

void gretl_string_table_add_extra (gretl_string_table *gst, PRN *prn);

int gretl_string_table_reset_column_id (gretl_string_table *gst, 
					int oldid, int newid);

series_table *gretl_string_table_detach_col (gretl_string_table *gst,
					     int col);

int in_string_table (gretl_string_table *gst, int id);

int *string_table_copy_list (gretl_string_table *gst);

int string_table_replace_list (gretl_string_table *gst,
			       int *newlist);

series_table *series_table_new (char **strs, int n_strs, int *err);

series_table *series_table_copy (series_table *st);

void series_table_destroy (series_table *st);

void series_table_free_shallow (series_table *st);

double series_table_get_value (series_table *st, const char *s);

const char *series_table_get_string (series_table *st, double val);

int series_table_add_string (series_table *st, const char *s);

int series_table_add_strings (series_table *st, const char **S,
			      int ns);

char **series_table_get_strings (series_table *st, int *n_strs);

int series_table_get_n_strings (series_table *st);

int *series_table_map (series_table *st_from, series_table *st_to);

void series_table_print (DATASET *dset, int i, PRN *prn);

void gretl_insert_builtin_string (const char *name, const char *s);

char *get_built_in_string_by_name (const char *name);

void builtin_strings_cleanup (void);

int process_string_command (const char *line, 
			    double ***pZ, DATASET *dset,
			    PRN *prn);

int substitute_named_strings (char *line, int *subst);

char *gretl_getenv (const char *key, int *defined, int *err);

char *retrieve_date_string (int t, const DATASET *dset, int *err);

gretl_array *retrieve_date_strings (const gretl_vector *v,
				    const DATASET *dset,
				    int *err);

char *retrieve_file_content (const char *fname, const char *codset,
			     int *err);

char *gretl_backtick (const char *arg, int *err);

#endif /* GRETL_STRING_TABLE_H */
