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

#define free_datainfo(p) clear_datainfo(p, 0); free(p);

/* functions follow */

DATAINFO *create_new_dataset (double **pZ,  /* data matrix */
			      int nvar,     /* number of variables */
			      int nobs,     /* observations per variable */
			      int markers   /* case markers or not? */
			      );

void clear_datainfo (DATAINFO *pdinfo, int resample);

int start_new_Z (double **pZ, DATAINFO *pdinfo, 
		 int resample);

int dateton (const char *date, const int pd, const char *startdate,
	     char *msg);

void ntodate (char *datestr, const int nt, const DATAINFO *pdinfo);

int get_info (const char *hdrfile, print_t *prn);

int write_data (const char *fname, const int *list, 
		double *Z, const DATAINFO *pdinfo, 
	        int opt);

int has_gz_suffix (const char *fname);

void gz_switch_ext (char *targ, char *src, char *ext);

int get_data (double **pZ, DATAINFO *pdinfo, 
	      PATHS *ppaths, 
	      const int data_file_open, char *msg, 
	      FILE *fp);

int open_nulldata (double **pZ, DATAINFO *pdinfo, 
		   const int data_file_open, const int length,
		   print_t *prn);

int import_csv (double **pZ, DATAINFO *pdinfo, 
                const char *fname, print_t *prn);

int add_case_markers (DATAINFO *pdinfo, const char *fname);

int import_box (double **pZ, DATAINFO *pdinfo, 
		const char *fname, print_t *prn);

int detect_filetype (char *fname, PATHS *ppaths, print_t *prn);



