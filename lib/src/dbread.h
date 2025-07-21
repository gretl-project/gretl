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

#include "gretl_string_table.h"

#define DB_DESCRIP_LEN 72  /* size of array to hold "# description" */
#define DB_INIT_ROWS 32

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
    char *descrip;
    int nobs;
    char stobs[OBSLEN];
    char endobs[OBSLEN];
    int pd;
    int offset;
    int err;
    int undated;
    void *data;
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
    gretl_string_table *gst;
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

void series_info_init (SERIESINFO *sinfo);

void series_info_set_description (SERIESINFO *sinfo,
				  const char *s);

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

int set_db_name (const char *fname, int filetype, PRN *prn);

const char *get_db_name (void);

int dbnomics_is_open (void);

int set_odbc_dsn (const char *line, PRN *prn);

int db_set_sample (const char *star, const char *stop, DATASET *dset);

int db_get_series (const char *line, DATASET *datainfo,
		   gretlopt opt, PRN *prn);

int db_delete_series_by_name (const char *line, PRN *prn);

int db_delete_series_by_number (const int *list, const char *fname);

int db_range_check (int db_pd,
		    const char *db_stobs,
		    const char *db_endobs,
		    const char *varname,
		    DATASET *dset);

int check_db_import_conversion (int pd, DATASET *dset);

int transcribe_db_data (DATASET *dset, int targv,
			const double *src, int pd,
			int nobs, char *stobs,
			CompactMethod cmethod);

int lib_spread_db_data (double **dbZ, SERIESINFO *sinfo,
			DATASET *dset, PRN *prn);

int lib_spread_dbnomics_data (DATASET *dset, DATASET *dbset,
			      PRN *prn);

int compact_dataset (DATASET *dset, int newpd,
		     CompactMethod default_method,
		     int wkstart, int repday);

int expand_dataset (DATASET *dset, int newpd);

int midas_days_per_period (int days_per_week, int pd);

#endif /* DBREAD_H */
