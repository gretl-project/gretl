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

typedef struct wbook_ wbook;
typedef struct wsheet_ wsheet;

typedef enum {
    BOOK_NUMERIC_DATES   = 1 << 0,
    BOOK_DATE_BASE_1904  = 1 << 1,
    BOOK_AUTO_VARNAMES   = 1 << 2,
    BOOK_TIME_SERIES     = 1 << 3,
    BOOK_OBS_LABELS      = 1 << 4,
    BOOK_OBS_BLANK       = 1 << 5,
    BOOK_DEBUG           = 1 << 6,
    BOOK_DATA_REVERSED   = 1 << 7,
    BOOK_TOP_LEFT_EMPTY  = 1 << 8
} BookFlag;

enum {
    COL_OFFSET,
    ROW_OFFSET
};

#define book_numeric_dates(b) ((b)->flags & BOOK_NUMERIC_DATES)
#define book_base_1904(b)     ((b)->flags & BOOK_DATE_BASE_1904)
#define book_auto_varnames(b) ((b)->flags & BOOK_AUTO_VARNAMES)
#define book_time_series(b)   ((b)->flags & BOOK_TIME_SERIES)
#define book_obs_labels(b)    ((b)->flags & BOOK_OBS_LABELS)
#define book_debugging(b)     ((b)->flags & BOOK_DEBUG)

#define book_set_numeric_dates(b) ((b)->flags |= BOOK_NUMERIC_DATES)
#define book_set_base_1904(b)     ((b)->flags |= BOOK_DATE_BASE_1904)
#define book_set_auto_varnames(b) ((b)->flags |= BOOK_AUTO_VARNAMES)
#define book_set_time_series(b)   ((b)->flags |= BOOK_TIME_SERIES)
#define book_set_obs_labels(b)    ((b)->flags |= BOOK_OBS_LABELS)
#define book_set_debug(b)         ((b)->flags |= BOOK_DEBUG)

#define book_unset_obs_labels(b)  ((b)->flags &= ~BOOK_OBS_LABELS)
#define book_unset_time_series(b) ((b)->flags &= ~BOOK_TIME_SERIES)

struct wbook_ {
    int version;
    int nsheets;
    int selected;
    int col_offset, row_offset;
    char *targname;
    char **sheetnames;
    guint32 *byte_offsets;
    void *colspin, *rowspin;
    int *xf_list;
    BookFlag flags;
    int (*get_min_offset)();
    void *data;
};

struct wsheet_ {
    int maxcol, maxrow;
    int text_cols, text_rows;
    int col_offset, row_offset;
    int colheads;
    int ID;
    BookFlag flags;
    char *name;
    double **Z;
    char **varname;
    char **label;
};

