/* 
 * Copyright (C) 2004 Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
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

enum {
    PRESERVE_OPG_MODEL = 1 << 0,
    FULL_VCV_MATRIX    = 1 << 1
} bhhh_opts;

typedef struct _model_info model_info;

void model_info_free (model_info *model);

model_info *model_info_new (void);

MODEL *model_info_capture_OPG_model (model_info *model);

gretl_matrix *model_info_get_VCV (model_info *model);

double *model_info_get_theta (model_info *model);

int model_info_get_t1 (const model_info *model);

int model_info_get_t2 (const model_info *model);

int model_info_get_n (const model_info *model);

void model_info_get_pqr (const model_info *model, 
			 int *p, int *q, int *r);

double model_info_get_ll (const model_info *model);

double **model_info_get_series (const model_info *model);

void model_info_set_pqr (model_info *model, int p, int q, int r);

void model_info_set_n_series (model_info *model, int n);

void model_info_set_k (model_info *model, int k);

int model_info_get_k (model_info *model);

void model_info_set_t1_t2 (model_info *model, int t1, int t2);

void model_info_set_opts (model_info *model, unsigned char opts);

void model_info_set_tol (model_info *model, double tol);

void model_info_set_ll (model_info *model, double ll, int do_score);

void model_info_set_s2 (model_info *model, double s2);

int bhhh_max (int (*loglik) (double *, const double **, double **,
			     model_info *, int), 
	      const double **X, const double *init_coeff,
	      model_info *model, PRN *prn);
