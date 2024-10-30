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
    GRETL_XML_DATA,       /* gretl XML data file (.gdt) */
    GRETL_BINARY_DATA,    /* gretl_binary data file (.gdtb) */
    GRETL_CSV,            /* comma-separated or other plain text data */
    GRETL_OCTAVE,         /* GNU octave ascii data file */
    GRETL_GNUMERIC,       /* gnumeric workbook data */
    GRETL_XLS,            /* MS Excel spreadsheet data */
    GRETL_XLSX,           /* MS Office Open XML spreadsheet data */
    GRETL_ODS,            /* Open Document Spreadsheet data */
    GRETL_WF1,            /* Eviews workfile data */
    GRETL_DTA,            /* Stata .dta data */
    GRETL_SAV,            /* SPSS .sav data */
    GRETL_SAS,            /* SAS xport data file */
    GRETL_JMULTI,         /* JMulTi data file */
    GRETL_DATA_MAX,       /* -- place marker -- */
    GRETL_SCRIPT,         /* file containing gretl commands */
    GRETL_SESSION,        /* zipped session file */
    GRETL_NATIVE_DB,      /* gretl database */
    GRETL_NATIVE_DB_WWW,  /* gretl database, accessed via internet */
    GRETL_RATS_DB,        /* RATS 4.0 database */
    GRETL_PCGIVE_DB,      /* PcGive bn7/in7 pair */
    GRETL_ODBC,           /* Open DataBase Connectivity */
    GRETL_DBNOMICS,       /* DB.NOMICS access */
    GRETL_MAP,            /* shapefile or GeoJSON */
    GRETL_UNRECOGNIZED    /* none of the above */
} GretlFileType;

typedef enum {
    CLEAR_FULL,           /* fully clear the dataset */
    CLEAR_SUBSAMPLE       /* dataset is sub-sampled: clear partially */
} DataClearCode;

#define SPREADSHEET_IMPORT(f) (f == GRETL_XLS ||	\
			       f == GRETL_XLSX ||	\
			       f == GRETL_GNUMERIC ||	\
			       f == GRETL_ODS)

#define OTHER_IMPORT(f) (f == GRETL_DTA ||	\
                         f == GRETL_SAV ||	\
			 f == GRETL_SAS ||	\
                         f == GRETL_JMULTI ||	\
                         f == GRETL_OCTAVE ||	\
			 f == GRETL_WF1 ||	\
			 f == GRETL_MAP)

#define free_datainfo(p) do { if (p != NULL) { clear_datainfo(p, 0); free(p); } \
                            } while (0);

#define DBNA  -999.0 /* missing value code for gretl databases */

#define GRETL_SCALAR_DIGITS 12

/* functions follow */

int dateton (const char *date, const DATASET *dset);

int merge_dateton (const char *date, const DATASET *dset);

char *ntolabel (char *datestr, int t, const DATASET *dset);

char *ntolabel_8601 (char *datestr, int t, const DATASET *dset);

int get_subperiod (int t, const DATASET *dset, int *err);

int get_precision (const double *x, int n, int placemax);

double get_date_x (int pd, const char *obs);

void date_maj_min (int t, const DATASET *dset, int *maj, int *min);

int write_data (const char *fname, int *list, const DATASET *dset,
		gretlopt opt, PRN *prn);

int gui_write_data (const char *fname, int *list, const DATASET *dset,
		    gretlopt opt);

int is_gzipped (const char *fname);

int gretl_is_pkzip_file (const char *fname);

gretlopt get_merge_opts (gretlopt opt);

int merge_or_replace_data (DATASET *dset0, DATASET **pdset1,
			   gretlopt opt, PRN *prn);

int gretl_seek_data (char *fname, DATASET *dset,
		     gretlopt opt, PRN *prn);

int gretl_get_data (const char *fname, DATASET *dset,
		    gretlopt opt, PRN *prn);

int open_nulldata (DATASET *dset, int data_status,
		   int length, gretlopt opt, PRN *prn);

int import_csv (const char *fname, DATASET *dset,
	        gretlopt opt, PRN *prn);

int import_spreadsheet (const char *fname, GretlFileType ftype,
			int *list, char *sheetname,
			DATASET *dset, gretlopt opt, PRN *prn);

int import_other (const char *fname, GretlFileType ftype,
		  DATASET *dset, gretlopt opt, PRN *prn);

int peek_at_csv (const char *fname, int n_lines, PRN *prn);

int gretl_read_purebin (const char *fname, DATASET *dset,
			gretlopt opt, PRN *prn);

int add_obs_markers_from_file (DATASET *dset, const char *fname);

int add_var_labels_from_file (DATASET *dset, const char *fname);

int save_var_labels_to_file (const DATASET *dset, const char *fname);

int dataset_has_var_labels (const DATASET *dset);

int read_or_write_var_labels (gretlopt opt, DATASET *dset, PRN *prn);

int read_or_write_obs_markers (gretlopt opt, DATASET *dset, PRN *prn);

int read_or_write_dset_description (gretlopt opt, DATASET *dset, PRN *prn);

GretlFileType data_file_type_from_name (const char *fname);

GretlFileType detect_filetype (char *fname, gretlopt opt);

int check_varname (const char *varname);

int check_identifier (const char *varname);

int check_atof (const char *numstr);

int check_atoi (const char *numstr);

int transpose_data (DATASET *dset);

void dataset_add_import_info (DATASET *dset, const char *fname,
			      GretlFileType type);

int analyse_daily_import (const DATASET *dset, PRN *prn);

void set_dset_matrix_target (gretl_matrix **pm);

void *get_dset_matrix_target (void);

#endif /* DATAIO_H */
