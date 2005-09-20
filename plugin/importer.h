/*
 *  Copyright (c) by Allin Cottrell 2002
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

typedef struct wbook_ wbook;
typedef struct wsheet_ wsheet;

typedef enum {
    FIRST_COL_DATE_FORMAT = 1 << 0,
    DATE_BASE_1904        = 1 << 1
} BookFlag;

#define book_numeric_dates(b) (b.flags & FIRST_COL_DATE_FORMAT)
#define book_base_1904(b) (b.flags & DATE_BASE_1904)

struct wbook_ {
    int version;
    int nsheets;
    int selected;
    int col_offset, row_offset;
    char **sheetnames;
    guint32 *byte_offsets;
    void *colspin, *rowspin;
    int *xf_list;
    BookFlag flags;
    int totmiss;
    char *missmask;
    int debug;
};

struct wsheet_ {
    int maxcol, maxrow;
    int text_cols, text_rows;
    int col_offset, row_offset;
    int ID;
    BookFlag flags;
    char *name;
    double **Z;
    char **varname;
    char **label;
};

