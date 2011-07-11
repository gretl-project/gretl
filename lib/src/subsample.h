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

DATASET *fetch_full_dataset (void);

void free_subsample_mask (char *s);

char *copy_subsample_mask (const char *src, int *err);

char *copy_datainfo_submask (const DATASET *dset, int *err);

int write_datainfo_submask (const DATASET *dset, FILE *fp);

int write_model_submask (const MODEL *pmod, FILE *fp);

int submask_cmp (const char *m1, const char *m2);

int attach_subsample_to_model (MODEL *pmod, const DATASET *dset);

int add_dataset_to_model (MODEL *pmod, const DATASET *dset,
			  gretlopt opt);

int restrict_sample (const char *line, const int *list,  
		     DATASET *dset, ExecState *state, 
		     gretlopt opt, PRN *prn);

int 
restrict_sample_from_mask (char *mask, DATASET *dset, gretlopt opt);

int complex_subsampled (void);

int get_full_length_n (void);

void set_dataset_resampled (DATASET *dset);

int dataset_is_resampled (const DATASET *dset);

int set_sample (const char *line, DATASET *dset);

int restore_full_sample (DATASET *dset, ExecState *state);

int backup_full_dataset (DATASET *dset);

int count_missing_values (const DATASET *dset, gretlopt opt, 
			  PRN *prn, int *err);

void maybe_free_full_dataset (const DATASET *dset);

int model_sample_problem (MODEL *pmod, const DATASET *dset);

void print_sample_obs (const DATASET *dset, PRN *prn);

void print_sample_status (const DATASET *dset, PRN *prn);

int data_report (const DATASET *dset, const char *fname, PRN *prn);

#endif /* SUBSAMPLE_H */
