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

#ifndef DATAIO_H
#define DATAIO_H

#include <stdio.h>
#include <string.h>

typedef enum {
    GRETL_NATIVE_DATA,    /* gretl native format data file */
    GRETL_XML_DATA,       /* gretl xml format data file */
    GRETL_CSV,            /* comma-separated data file */
    GRETL_OCTAVE,         /* GNU octave ascii data file */
    GRETL_GNUMERIC,       /* gnumeric workbook data */
    GRETL_XLS,            /* MS Excel spreadsheet data */
    GRETL_ODS,            /* Open Document Spreadsheet data */
    GRETL_WF1,            /* Eviews workfile data */
    GRETL_DTA,            /* Stata .dta data */
    GRETL_SCRIPT,         /* file containing gretl commands */
    GRETL_SESSION,        /* zipped session file */
    GRETL_NATIVE_DB,      /* gretl database */
    GRETL_NATIVE_DB_WWW,  /* gretl database, accessed via internet */
    GRETL_RATS_DB,        /* RATS 4.0 database */
    GRETL_PCGIVE_DB,      /* PcGive bn7/in7 pair */
    GRETL_JMULTI,         /* JMulTi data file */
    GRETL_UNRECOGNIZED    /* none of the above */
} GretlFileType;

typedef enum {
    CLEAR_FULL,           /* fully clear the dataset */
    CLEAR_SUBSAMPLE       /* dataset is sub-sampled: clear partially */
} DataClearCode;

typedef enum {
    VARNAME_RESERVED = 1, /* vername is a gretl reserved name */
    VARNAME_FIRSTCHAR,    /* first character is not alphabetical */
    VARNAME_BADCHAR       /* illegal character in second or subsequent place */
} GretlVarnameError;

#define SPREADSHEET_IMPORT(f) (f == GRETL_XLS || \
			       f == GRETL_GNUMERIC || \
			       f == GRETL_ODS)

#define OTHER_IMPORT(f) (f == GRETL_DTA || \
                         f == GRETL_JMULTI || \
                         f == GRETL_OCTAVE || \
			 f == GRETL_WF1)

#define free_datainfo(p) do { if (p != NULL) { clear_datainfo(p, 0); free(p); } \
                            } while (0);

#define DBNA  -999.0 /* missing value code for gretl databases */

#define GRETL_SCALAR_DIGITS 12

/* functions follow */

int dateton (const char *date, const DATAINFO *pdinfo);

char *ntodate (char *datestr, int nt, const DATAINFO *pdinfo);

char *ntodate_full (char *datestr, int t, const DATAINFO *pdinfo);

int get_subperiod (int t, const DATAINFO *pdinfo, int *err);

int get_info (const char *hdrfile, PRN *prn);

int get_precision (const double *x, int n, int placemax);

double get_date_x (int pd, const char *obs);

int write_data (const char *fname, const int *list, 
		const double **Z, const DATAINFO *pdinfo, 
	        gretlopt opt, PATHS *ppaths);

int data_report (const DATAINFO *pdinfo, PATHS *ppaths, PRN *prn);

int is_gzipped (const char *fname);

int gretl_is_pkzip_file (const char *fname);

void gz_switch_ext (char *targ, char *src, char *ext);

int merge_or_replace_data (double ***pZ0, DATAINFO *pdinfo0,
			   double ***pZ1, DATAINFO **ppdinfo1,
			   gretlopt opt, PRN *prn);

int gretl_get_data (char *datfile, PATHS *ppaths,
		    double ***pZ, DATAINFO *pdinfo, 
		    gretlopt opt, PRN *prn);

int open_nulldata (double ***pZ, DATAINFO *pdinfo, 
		   int data_status, int length,
		   PRN *prn);

int import_csv (const char *fname, double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn);

int import_spreadsheet (const char *fname, int ftype,
			int *list, char *sheetname,
			double ***pZ, DATAINFO *pdinfo, 
			gretlopt opt, PRN *prn);

int import_other (const char *fname, int ftype,
		  double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn);

int add_obs_markers_from_file (DATAINFO *pdinfo, const char *fname);

GretlFileType detect_filetype (char *fname, PATHS *ppaths, PRN *prn);

gretlopt data_save_opt_from_suffix (const char *fname);

int check_varname (const char *varname);

int check_atof (const char *numstr);

int check_atoi (const char *numstr);

int transpose_data (double ***pZ, DATAINFO *pdinfo);

#endif /* DATAIO_H */
