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

/* subsample.h for gretl */

#ifndef SUBSAMPLE_H
#define SUBSAMPLE_H

typedef enum {
    SUBSAMPLE_NONE,
    SUBSAMPLE_DROP_MISSING,
    SUBSAMPLE_USE_DUMMY,
    SUBSAMPLE_BOOLEAN,
    SUBSAMPLE_RANDOM,
    SUBSAMPLE_UNKNOWN
} SubsampleMode;

double ***fetch_full_Z (void);

void reset_full_Z (double ***pZ);

DATAINFO *fetch_full_datainfo (void);

char *copy_subsample_mask (const char *src);

char *copy_datainfo_submask (const DATAINFO *pdinfo);

int write_datainfo_submask (const DATAINFO *pdinfo, FILE *fp);

int write_model_submask (const MODEL *pmod, FILE *fp);

int submask_cmp (const char *m1, const char *m2);

int attach_subsample_to_model (MODEL *pmod, const DATAINFO *pdinfo);

int restrict_sample (const char *line, const int *list,  
		     double ***pZ, DATAINFO **ppdinfo,
		     ExecState *state, gretlopt opt, 
		     PRN *prn);

int 
restrict_sample_from_mask (char *mask, int mode, 
			   double ***pZ, DATAINFO **ppdinfo,
			   ExecState *state);

int complex_subsampled (void);

int get_full_length_n (void);

int set_sample (const char *line, const double **Z, DATAINFO *pdinfo);

int restore_full_sample (double ***pZ, DATAINFO **ppdinfo,
			 ExecState *state); 

int count_missing_values (double ***pZ, DATAINFO *pdinfo, PRN *prn);

int add_dataset_to_model (MODEL *pmod, const DATAINFO *pdinfo);

void free_model_dataset (MODEL *pmod);

void maybe_free_full_dataset (const DATAINFO *pdinfo);

#endif /* SUBSAMPLE_H */
