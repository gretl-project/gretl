/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2004 Ramu Ramanathan and Allin Cottrell
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

#ifndef GRETL_MODEL_H
#define GRETL_MODEL_H

/**
 * free_model:
 * @p: pointer to #MODEL.
 *
 * Free allocated content of MODEL then the pointer itself.
 */

#define free_model(p) if (p != NULL) { \
                             clear_model(p); \
                             free(p); \
                          }

#define AR_MODEL(c) (c == AR || c == CORC || c == PWE || c == HILU)


MODEL *gretl_model_new (void);

void gretl_model_init (MODEL *pmod);

void gretl_model_smpl_init (MODEL *pmod, const DATAINFO *pdinfo);

void gretl_model_set_auxiliary (MODEL *pmod, int aux);

void exchange_smpl (MODEL *pmod, DATAINFO *pdinfo);

void clear_model (MODEL *pmod);

int gretl_model_set_data (MODEL *pmod, const char *key, void *ptr, size_t size);

int gretl_model_set_int (MODEL *pmod, const char *key, int val);

int gretl_model_set_double (MODEL *pmod, const char *key, double val);

void *gretl_model_get_data (const MODEL *pmod, const char *key);

int gretl_model_get_int (const MODEL *pmod, const char *key);

double gretl_model_get_double (const MODEL *pmod, const char *key);

void debug_print_model_info (const MODEL *pmod, const char *msg);

int copy_model (MODEL *targ, const MODEL *src, const DATAINFO *pdinfo);

int swap_models (MODEL **targ, MODEL **src);

int command_ok_for_model (int test_ci, int model_ci);

int get_model_count (void);

void reset_model_count (void);

int model_count_plus (void);

void model_count_minus (void);

void set_model_id (MODEL *pmod);

void gretl_test_init (GRETLTEST *test);

int add_test_to_model (MODEL *pmod, const GRETLTEST *test);

void model_list_to_string (int *list, char *buf);

int highest_numbered_var_in_model (const MODEL *pmod, 
				   const DATAINFO *pdinfo);

int mle_aic_bic (MODEL *pmod, int addk);

#endif /* GRETL_MODEL_H */
