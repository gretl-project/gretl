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

typedef enum {
    PRESERVE_OPG_MODEL = 1 << 0,  
    FULL_VCV_MATRIX    = 1 << 1
} BHHH_opts;

typedef struct _model_info model_info;

typedef int (*LL_FUNC) (double *, 
			const double **, 
			double **, 
			model_info *, 
			int);

void model_info_free (model_info *minfo);

model_info *model_info_new (int k, int t1, int t2, int bign, double tol);

MODEL *model_info_capture_OPG_model (model_info *minfo);

gretl_matrix *model_info_get_VCV (model_info *minfo);

double *model_info_get_theta (model_info *minfo);

int model_info_get_t1 (const model_info *minfo);

int model_info_get_t2 (const model_info *minfo);

int model_info_get_n (const model_info *minfo);

int model_info_get_iters (const model_info *minfo);

void *model_info_get_extra_info (model_info *minfo);

double model_info_get_ll (const model_info *minfo);

double **model_info_get_series (const model_info *minfo);

void model_info_set_extra_info (model_info *minfo, void *extra);

void model_info_set_n_series (model_info *minfo, int n);

int model_info_get_k (model_info *minfo);

void model_info_set_opts (model_info *minfo, unsigned char opts);

void model_info_set_ll (model_info *minfo, double ll, int do_score);

int bhhh_max (LL_FUNC loglik, 
	      const double **X, 
	      const double *init_coeff,
	      model_info *minfo, 
	      PRN *prn);
