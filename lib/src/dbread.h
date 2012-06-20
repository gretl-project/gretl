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
} DBError;

/**
 * dbnumber: 
 *
 * The type used for representing primary data in gretl databases.
 */

typedef float dbnumber;

typedef struct SERIESINFO_ {
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
} SERIESINFO;

typedef struct dbwrapper_ {
    char *fname;
    int dbtype;
    int nv;
    int nalloc;
    SERIESINFO *sinfo;
} dbwrapper;

#define ODBC_OBSCOLS 3

typedef struct ODBC_info_ {
    char *dsn;
    char *username;
    char *password;
    char *query;
    char **fmts;
    char coltypes[ODBC_OBSCOLS];
    double **X;
    char **S;
    int nrows;
    int obscols;
    int nvars;
} ODBC_info;

#if G_BYTE_ORDER == G_BIG_ENDIAN
typedef struct {
    long frac;
    short exp;
} netfloat;

float retrieve_float (netfloat nf);
#endif

int get_native_db_data (const char *dbbase, SERIESINFO *sinfo, 
			double **Z);

int get_remote_db_data (const char *dbbase, SERIESINFO *sinfo, 
			double **Z);

int get_pcgive_db_data (const char *dbbase, SERIESINFO *sinfo, 
			double **Z);

int get_rats_db_data (const char *fname, SERIESINFO *sinfo, double **Z);

dbwrapper *read_rats_db (const char *fname, FILE *fp);

dbwrapper *read_pcgive_db (const char *fname, FILE *fp);

dbwrapper *dbwrapper_new (int n, const char *fname, int dbtype);

void dbwrapper_destroy (dbwrapper *dw);

double *compact_db_series (const double *src, SERIESINFO *sinfo,
			   int target_pd, CompactMethod method);

double *expand_db_series (const double *src, SERIESINFO *sinfo,
			  int target_pd, int interpol);

int set_db_name (const char *fname, int filetype, PRN *prn);

const char *get_db_name (void);

int set_odbc_dsn (const char *line, PRN *prn);

int db_set_sample (const char *s, DATASET *dset);

int db_get_series (char *line, DATASET *datainfo, 
		   gretlopt opt, PRN *prn);

int db_delete_series_by_name (char *line, PRN *prn);

int db_delete_series_by_number (const int *list, const char *fname);

void get_db_padding (SERIESINFO *sinfo, DATASET *dset, 
		     int *pad1, int *pad2);

int db_range_check (SERIESINFO *sinfo, DATASET *dset);

int check_db_import (SERIESINFO *sinfo, DATASET *dset);

int compact_data_set (DATASET *dset, int newpd,
		      CompactMethod default_method, 
		      int monstart, int repday);

int expand_data_set (DATASET *dset, int newpd,
		     int interpol);

#endif /* DBREAD_H */
