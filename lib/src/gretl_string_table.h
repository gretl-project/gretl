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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
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

int process_string_command (const char *line, PRN *prn);

int substitute_named_strings (char *line);

int string_is_defined (const char *sname);

void saved_strings_cleanup (void);

#endif /* GRETL_STRING_TABLE_H */
