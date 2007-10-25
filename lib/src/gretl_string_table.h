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

typedef struct _gretl_string_table gretl_string_table;

gretl_string_table *gretl_string_table_new (void);

gretl_string_table *string_table_new_from_cols_list (int *list);

int gretl_string_table_index (gretl_string_table *st, const char *s, int col,
			      int addcol, PRN *prn);

int gretl_string_table_print (gretl_string_table *st, DATAINFO *pdinfo,
			      const char *fname, PRN *prn);

void gretl_string_table_destroy (gretl_string_table *st);

void gretl_insert_builtin_string (const char *name, const char *s);

int save_named_string (const char *name, const char *s, PRN *prn);

char *get_named_string (const char *name);

int process_string_command (const char *line, PRN *prn);

int substitute_named_strings (char *line);

int string_is_defined (const char *sname);

void saved_strings_cleanup (void);

int destroy_saved_strings_at_level (int d);

int is_codevar (const char *s);

int set_codevars (const char *s);

#endif /* GRETL_STRING_TABLE_H */
