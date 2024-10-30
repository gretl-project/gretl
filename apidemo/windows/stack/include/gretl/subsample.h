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

#define RESAMPLED ((char *) 0xdeadbeef)

DATASET *fetch_full_dataset (void);

void sync_dataset_shared_members (const DATASET *dset);

void free_subsample_mask (char *s);

char *copy_subsample_mask (const char *src, int *err);

char *copy_dataset_submask (const DATASET *dset, int *err);

int write_dataset_submask (const DATASET *dset, PRN *prn);

int write_model_submask (const MODEL *pmod, PRN *prn);

int get_dataset_submask_size (const DATASET *dset);

int get_model_submask_size (const MODEL *pmod);

int subsample_check_model (MODEL *pmod, char *mask);

int remove_model_subsample_info (MODEL *pmod);

int submask_cmp (const char *m1, const char *m2);

int attach_subsample_to_model (MODEL *pmod, const DATASET *dset);

int add_dataset_to_model (MODEL *pmod, const DATASET *dset,
			  gretlopt opt);

int restrict_sample (const char *param, const int *list,  
		     DATASET *dset, ExecState *state, 
		     gretlopt opt, PRN *prn,
		     int *n_dropped);

int restrict_sample_from_mask (char *mask, DATASET *dset, 
			       gretlopt opt);

int perma_sample (DATASET *dset, gretlopt opt, PRN *prn,
		  int *n_dropped);

int recode_strvals (DATASET *dset, gretlopt opt);

int complex_subsampled (void);

int dataset_is_subsampled (const DATASET *dset);

int get_full_length_n (void);

void set_dataset_resampled (DATASET *dset, unsigned int seed);

int dataset_is_resampled (const DATASET *dset);

int set_sample (const char *start, const char *stop,
		DATASET *dset, gretlopt opt);

int set_panel_sample (const char *start, const char *stop,
		      gretlopt opt, DATASET *dset,
		      ExecState *s, PRN *prn);

int restore_full_sample (DATASET *dset, ExecState *state);

int backup_full_dataset (DATASET *dset);

int count_missing_values (const DATASET *dset, gretlopt opt, 
			  PRN *prn, int *err);

void maybe_free_full_dataset (DATASET *dset);

int model_sample_problem (const MODEL *pmod, const DATASET *dset);

int fcast_not_feasible (const MODEL *pmod, const DATASET *dset);

int same_dataset (const MODEL *pmod, const DATASET *dset);

void print_sample_obs (const DATASET *dset, PRN *prn);

void print_sample_status (const DATASET *dset, PRN *prn);

int data_report (const DATASET *dset, const char *fname, PRN *prn);

#endif /* SUBSAMPLE_H */
