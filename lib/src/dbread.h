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

#ifndef DBREAD_H
#define DBREAD_H

#define DB_DESCRIP_LEN 72  /* size of array to hold "# description" */

typedef enum {
    DB_MISSING_DATA = E_MAX + 1,
    DB_NO_SUCH_SERIES,
    DB_PARSE_ERROR
} db_error_codes;

typedef enum {
    COMPACT_NONE,
    COMPACT_SUM,
    COMPACT_AVG,
    COMPACT_SOP,
    COMPACT_EOP,
    COMPACT_WDAY,
    COMPACT_MAX
} CompactMethod; 

typedef float dbnumber;
typedef struct dbwrapper_ dbwrapper;
typedef struct SERIESINFO_ SERIESINFO;
typedef struct ODBC_info_ ODBC_info;

struct SERIESINFO_ {
    int t1, t2, v;
    char varname[VNAMELEN];
    char descrip[MAXLABEL];
    int nobs;
    char stobs[OBSLEN];
    char endobs[OBSLEN];
    int pd;
    int offset;
    int err;
    int undated;
};

struct dbwrapper_ {
    int nv;
    int nalloc;
    SERIESINFO *sinfo;
};

#define ODBC_MAXCOLS 4

struct ODBC_info_ {
    char *dsn;
    char *username;
    char *password;
    char *query;
    char *fmt;
    char coltypes[ODBC_MAXCOLS];
    double *x;
    char **S;
    int nrows;
    int ncols;
};

int get_native_db_data (const char *dbbase, SERIESINFO *sinfo, 
			double **Z);

int get_remote_db_data (const char *dbbase, SERIESINFO *sinfo, double **Z);

int get_pcgive_db_data (const char *dbbase, SERIESINFO *sinfo, 
			double **Z);

int get_rats_db_data (const char *fname, SERIESINFO *sinfo, double **Z);

dbwrapper *read_rats_db (FILE *fp);

dbwrapper *read_pcgive_db (FILE *fp);

dbwrapper *dbwrapper_new (int n);

void dbwrapper_destroy (dbwrapper *dw);

double *compact_db_series (const double *src, SERIESINFO *sinfo,
			   int target_pd, CompactMethod method);

double *expand_db_series (const double *src, SERIESINFO *sinfo,
			  int target_pd);

int set_db_name (const char *fname, int filetype, const PATHS *ppaths, 
		 PRN *prn);

const char *get_db_name (void);

int set_odbc_dsn (const char *line, PRN *prn);

int db_set_sample (const char *s, DATAINFO *pdinfo);

int db_get_series (char *line, double ***pZ, DATAINFO *datainfo, 
		   gretlopt opt, PRN *prn);

int db_delete_series_by_name (char *line, PRN *prn);

int db_delete_series_by_number (const int *list, const char *fname);

void get_db_padding (SERIESINFO *sinfo, DATAINFO *pdinfo, 
		     int *pad1, int *pad2);

int db_range_check (SERIESINFO *sinfo, DATAINFO *pdinfo);

int check_db_import (SERIESINFO *sinfo, DATAINFO *pdinfo);

int compact_data_set (double ***pZ, DATAINFO *pdinfo, int newpd,
		      CompactMethod default_method, int monstart,
		      int repday);

int expand_data_set (double ***pZ, DATAINFO *pdinfo, int newpd);

#endif /* DBREAD_H */
