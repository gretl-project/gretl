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

#ifndef LIBSET_H
#define LIBSET_H

typedef enum {
    VCV_UNSET,
    VCV_HESSIAN,
    VCV_IM,
    VCV_OP,
    VCV_QML,
    VCV_BW
} VcvType;

typedef enum {
    KERNEL_BARTLETT = 1,
    KERNEL_PARZEN,
    KERNEL_QS
} HACKernel;

typedef enum {
    NORM_PHILLIPS,
    NORM_DIAG,
    NORM_FIRST
} VECMnorm;

typedef int (*ITER_PRINT_FUNC) (PRN *);

int libset_init (void);
void libset_cleanup (void);
int libset_restore_state_zero (DATAINFO *pdinfo);

int push_program_state (void);
int pop_program_state (void);

void set_use_qr (int set);
int get_use_qr (void);

void set_use_lbfgs (int set);
int get_use_lbfgs (void);

void set_use_cwd (int set);
int get_use_cwd (void);

void set_shell_ok (int set);
int get_shell_ok (void);

void set_xsect_hccme (const char *s);
void set_tseries_hccme (const char *s);
void set_panel_hccme (const char *s);

void set_garch_robust_vcv (const char *s);
int get_garch_vcv_version (void);
int get_garch_robust_vcv_version (void);

int get_force_hc (void);
int get_hc_version (void);
int get_hac_lag (int T);
int get_hac_kernel (void);
void set_hac_kernel (int k);

int get_hac_prewhiten (void);
void set_hac_prewhiten (int w);

int get_panel_beck_katz (void);
void set_panel_beck_katz (int b);

double get_qs_bandwidth (void);
void set_qs_bandwidth (double w);

int data_based_hac_bandwidth (void);

int get_halt_on_error (void);

double get_hp_lambda (void);
int set_hp_lambda (double d);

int get_bkbp_k (const DATAINFO *pdinfo);
void get_bkbp_periods (const DATAINFO *pdinfo, int *l, int *u);
int set_bkbp_k (int k);
int set_bkbp_periods (int bkl, int bku);
void unset_bkbp_k (void);
void unset_bkbp_periods (void);

int gretl_get_text_pause (void);

double get_bhhh_toler (void);
int get_bhhh_maxiter (void);
int set_bhhh_toler (double tol);
int set_bhhh_maxiter (int n);

double get_bfgs_toler (void);
int get_bfgs_maxiter (void);
int set_bfgs_toler (double tol);
int set_bfgs_maxiter (int n);

const gretl_matrix *get_init_vals (void);
int n_init_vals (void);
void free_init_vals (void);

int get_VAR_horizon (void);

int get_bootstrap_replications (void);

double get_nls_toler (void);
int set_nls_toler (double tol);

void set_loop_on (void);
void set_loop_off (void);
int gretl_looping (void);

void gretl_set_batch_mode (int b);
int gretl_in_batch_mode (void);

void gretl_set_gui_mode (int g);
int gretl_in_gui_mode (void);

void set_gretl_echo (int e);
int gretl_echo_on (void);

void set_gretl_messages (int e);
int gretl_messages_on (void);

int set_long_digits (int n);
int get_long_digits (void);

int set_max_verbose (int n);
int get_max_verbose (void);

int get_vecm_norm (void);

void shelldir_init (void);
char *get_shelldir (void);

char get_csv_delim (const DATAINFO *pdinfo);

int execute_set_line (const char *line, double **Z, DATAINFO *pdinfo, 
		      PRN *prn);

void set_iter_print_func (ITER_PRINT_FUNC func);
int iter_print_callback (PRN *prn);

#endif /* LIBSET_H */

