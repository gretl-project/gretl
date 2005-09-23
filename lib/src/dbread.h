/*
 *   Copyright (c) by Allin Cottrell
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

#ifndef DBREAD_H
#define DBREAD_H

#define DB_DESCRIP_LEN 72  /* size of array to hold "# description" */

typedef enum {
    DB_OK = 0,
    DB_MISSING_DATA,
    DB_NO_SUCH_SERIES,
    DB_PARSE_ERROR,
    DB_NOT_FOUND
} db_error_codes;

typedef enum {
    COMPACT_NONE,
    COMPACT_SUM,
    COMPACT_AVG,
    COMPACT_SOP,
    COMPACT_EOP,
    COMPACT_MAX
} compaction_methods; 

typedef float dbnumber;

typedef struct _db_table_row db_table_row;
typedef struct _db_table db_table;

struct _db_table_row {
    char *varname;
    char *comment;
    char *obsinfo;
};

struct _db_table {
    int nrows;
    int nalloc;
    db_table_row *rows;
};

typedef struct _SERIESINFO SERIESINFO;

struct _SERIESINFO {
    char varname[16];
    char descrip[MAXLABEL];
    int nobs;
    char stobs[OBSLEN];
    char endobs[OBSLEN];
    int pd;
    int offset;
    int err;
    int undated;
};

int get_native_db_data (const char *dbbase, SERIESINFO *sinfo, 
			double **Z);

db_table *read_rats_db (FILE *fp);

int get_rats_data_by_series_number (const char *fname, 
				    int series_number,
				    SERIESINFO *sinfo, 
				    double **Z);

double *compact_db_series (const double *src, SERIESINFO *sinfo,
			   int target_pd, int method);

double *expand_db_series (const double *src, SERIESINFO *sinfo,
			  int target_pd);

int set_db_name (const char *fname, int filetype, const PATHS *ppaths, 
		 PRN *prn);

int db_set_sample (const char *line, DATAINFO *pdinfo);

int db_get_series (const char *line, double ***pZ, DATAINFO *datainfo, 
		   PRN *prn);

void get_db_padding (SERIESINFO *sinfo, DATAINFO *pdinfo, 
		     int *pad1, int *pad2);

int check_db_import (SERIESINFO *sinfo, DATAINFO *pdinfo);

int compact_data_set (double ***pZ, DATAINFO *pdinfo, int newpd,
		      int default_method, int monstart);

int expand_data_set (double ***pZ, DATAINFO *pdinfo, int newpd);

#endif /* DBREAD_H */
