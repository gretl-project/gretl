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

int gretl_string_table_index (gretl_string_table *st, const char *s, int col,
			      int addcol, PRN *prn);

int gretl_string_table_print (gretl_string_table *st, PATHS *ppaths, PRN *prn);

void gretl_string_table_destroy (gretl_string_table *st);

#endif /* GRETL_STRING_TABLE_H */
