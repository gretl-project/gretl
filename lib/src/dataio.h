/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <stdio.h>
#include <string.h>

typedef enum {
    GRETL_DATA_FLOAT = 1, /* single-precision binary data */
    GRETL_DATA_DOUBLE,    /* double-precision binary data */
    GRETL_DATA_OCTAVE,    /* data in Gnu Octave format */
    GRETL_DATA_CSV,       /* data in Comma Separated Values format */
    GRETL_DATA_R,         /* data in Gnu R format */
    GRETL_DATA_R_ALT,     /* data in alternate Gnu R format */
    GRETL_DATA_GZIPPED,   /* gzipped data */
    GRETL_DATA_TRAD       /* traditional (ESL-style) data */
} data_options;

typedef enum {
    GRETL_NATIVE_DATA,    /* gretl native format data file */
    GRETL_XML_DATA,       /* gretl xml format data file */
    GRETL_CSV_DATA,       /* comma-separated data file */
    GRETL_BOX_DATA,       /* BOX1 format data file */
    GRETL_GNUMERIC,       /* gnumeric workbook data */
    GRETL_EXCEL,          /* MS Excel spreadsheet data */
    GRETL_SCRIPT,         /* file containing gretl commands */
    GRETL_NATIVE_DB,      /* gretl database */
    GRETL_RATS_DB,        /* RATS 4.0 database */
    GRETL_UNRECOGNIZED    /* none of the above */
} gretl_filetypes;

typedef enum {
    CLEAR_FULL,
    CLEAR_SUBSAMPLE
} clear_codes;

typedef enum {
    DATA_NONE,
    DATA_CLEAR,
    DATA_APPEND
} data_open_codes;


#define free_datainfo(p) do { if (p != NULL) { clear_datainfo(p, 0); free(p); } \
                            } while (0);

/* functions follow */

void free_Z (double **Z, DATAINFO *pdinfo);

DATAINFO *datainfo_new (void);

DATAINFO *create_new_dataset (double ***pZ, /* data matrix */
			      int nvar,     /* number of variables */
			      int nobs,     /* observations per variable */
			      int markers   /* case markers or not? */
			      );

void clear_datainfo (DATAINFO *pdinfo, int code);

int start_new_Z (double ***pZ, DATAINFO *pdinfo, int resample);

int dateton (const char *date, const DATAINFO *pdinfo);

char *ntodate (char *datestr, int nt, const DATAINFO *pdinfo);

char *ntodate_full (char *datestr, int t, const DATAINFO *pdinfo);

int get_info (const char *hdrfile, PRN *prn);

int get_precision (double *x, int n, int placemax);

double get_date_x (int pd, const char *obs);

int write_data (const char *fname, const int *list, 
		double **Z, const DATAINFO *pdinfo, 
	        int opt, PATHS *ppaths);

int data_report (const DATAINFO *pdinfo, PATHS *ppaths, PRN *prn);

int is_gzipped (const char *fname);

int has_gz_suffix (const char *fname);

void gz_switch_ext (char *targ, char *src, char *ext);

int merge_data (double ***pZ, DATAINFO *pdinfo,
		double **addZ, DATAINFO *addinfo,
		PRN *prn);

int gretl_get_data (double ***pZ, DATAINFO **ppdinfo, 
		    char *datfile, PATHS *ppaths, 
		    int data_status, PRN *prn);

int open_nulldata (double ***pZ, DATAINFO *pdinfo, 
		   int data_status, int length,
		   PRN *prn);

int import_csv (double ***pZ, DATAINFO **ppdinfo, 
                const char *fname, PRN *prn);

int import_box (double ***pZ, DATAINFO **ppdinfo, 
		const char *fname, PRN *prn);

int add_case_markers (DATAINFO *pdinfo, const char *fname);

int detect_filetype (char *fname, PATHS *ppaths, PRN *prn);

int check_varname (const char *varname);

int get_xmldata (double ***pZ, DATAINFO **ppdinfo, char *fname,
		 PATHS *ppaths, int data_status, PRN *prn, int gui); 

char *get_xml_description (const char *fname);

int check_atof (const char *numstr);
