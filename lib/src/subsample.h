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

#define RESAMPLED ((char *) 0xdeadbeef)

double ***fetch_full_Z (void);

void reset_full_Z (double ***pZ);

DATAINFO *fetch_full_datainfo (void);

void free_subsample_mask (char *s);

char *copy_subsample_mask (const char *src, int *err);

char *copy_datainfo_submask (const DATAINFO *pdinfo, int *err);

int write_datainfo_submask (const DATAINFO *pdinfo, FILE *fp);

int write_model_submask (const MODEL *pmod, FILE *fp);

int submask_cmp (const char *m1, const char *m2);

int attach_subsample_to_model (MODEL *pmod, const DATAINFO *pdinfo);

int restrict_sample (const char *line, const int *list,  
		     double ***pZ, DATAINFO *pdinfo,
		     ExecState *state, gretlopt opt, 
		     PRN *prn);

int 
restrict_sample_from_mask (char *mask, double ***pZ, DATAINFO *pdinfo,
			   gretlopt opt);

int complex_subsampled (void);

int get_full_length_n (void);

void set_dataset_resampled (DATAINFO *pdinfo);

int dataset_is_resampled (const DATAINFO *pdinfo);

int set_sample (const char *line, double ***pZ, DATAINFO *pdinfo);

int restore_full_sample (double ***pZ, DATAINFO *pdinfo, ExecState *state);

int backup_full_dataset (double **Z, DATAINFO *pdinfo);

int count_missing_values (const double **Z, const DATAINFO *pdinfo, 
			  gretlopt opt, PRN *prn, int *err);

int add_dataset_to_model (MODEL *pmod, const double **Z, 
			  const DATAINFO *pdinfo,
			  gretlopt opt);

void free_model_dataset (MODEL *pmod);

void maybe_free_full_dataset (const DATAINFO *pdinfo);

int model_sample_problem (MODEL *pmod, const DATAINFO *pdinfo);

#endif /* SUBSAMPLE_H */
