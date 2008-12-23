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
    NORM_PHILLIPS,
    NORM_DIAG,
    NORM_FIRST,
    NORM_NONE,
    NORM_MAX
} VECMnorm;

/* guard against consequences of typos */

#define BFGS_MAXITER     "bfgs_maxiter"
#define BFGS_TOLER       "bfgs_toler"
#define BHHH_MAXITER     "bhhh_maxiter"
#define BHHH_TOLER       "bhhh_toler"
#define BKBP_K           "bkbp_k"
#define BOOTREP          "bootrep"
#define FORCE_DECP       "force_decpoint"
#define FORCE_HC         "force_hc"
#define GARCH_VCV        "garch_vcv"
#define GARCH_ROBUST_VCV "garch_robust_vcv"
#define HAC_KERNEL       "hac_kernel"
#define HAC_LAG          "hac_lag"
#define HALT_ON_ERR      "halt_on_error"
#define HC_VERSION       "hc_version"
#define HORIZON          "horizon"
#define HP_LAMBDA        "hp_lambda"
#define USE_LBFGS        "lbfgs"
#define LONGDIGITS       "longdigits"
#define LOOP_MAXITER     "loop_maxiter"
#define RQ_MAXITER       "rq_maxiter"
#define MAX_VERBOSE      "max_verbose"
#define NLS_TOLER        "nls_toler"
#define PCSE             "pcse"
#define PREWHITEN        "hac_prewhiten"
#define QS_BANDWIDTH     "qs_bandwidth"
#define SHELL_OK         "shell_ok"
#define USE_CWD          "use_cwd"
#define USE_QR           "qr"
#define USE_FCP          "fcp"
#define VECM_NORM        "vecm_norm"
#define ARMA_VCV         "arma_vcv"
#define VERBOSE_INCLUDE  "verbose_include"

typedef int (*ITER_PRINT_FUNC) (int, PRN *);

#define set_nls_toler(x) (libset_set_double(NLS_TOLER, x))

int libset_init (void);
void libset_cleanup (void);
int libset_restore_state_zero (DATAINFO *pdinfo);

int push_program_state (void);
int pop_program_state (void);

int libset_get_bool (const char *s);
void libset_set_bool (const char *s, int set);

double libset_get_double (const char *s);
int libset_set_double (const char *s, double x);

int libset_get_int (const char *s);
int libset_set_int (const char *s, int k);

/* GUI setter functions */
void set_xsect_hccme (const char *s);
void set_tseries_hccme (const char *s);
void set_panel_hccme (const char *s);
void set_garch_robust_vcv (const char *s);

int get_hac_lag (int T);
int data_based_hac_bandwidth (void);

int get_bkbp_k (const DATAINFO *pdinfo);
void get_bkbp_periods (const DATAINFO *pdinfo, int *l, int *u);
int set_bkbp_k (int k);
int set_bkbp_periods (int bkl, int bku);
void unset_bkbp_k (void);
void unset_bkbp_periods (void);

void set_mp_bits (int b);
int get_mp_bits (void);

const gretl_matrix *get_init_vals (void);
int n_init_vals (void);
void free_init_vals (void);

void set_loop_on (int quiet);
void set_loop_off (void);
int gretl_looping (void);
int gretl_looping_quietly (void);

void gretl_set_batch_mode (int b);
int gretl_in_batch_mode (void);

void gretl_set_gui_mode (int g);
int gretl_in_gui_mode (void);

void set_gretl_echo (int e);
int gretl_echo_on (void);

void set_gretl_messages (int e);
int gretl_messages_on (void);

int gretl_warnings_on (void);
int gretl_debugging_on (void);
int lists_protected (void);

void shelldir_init (const char *s);
char *get_shelldir (void);

char get_csv_delim (const DATAINFO *pdinfo);

int execute_set_line (const char *line, double **Z, DATAINFO *pdinfo, 
		      PRN *prn);

void set_iter_print_func (ITER_PRINT_FUNC func);
int iter_print_callback (int i, PRN *prn);
int iter_print_func_installed (void);

void set_workdir_callback (int (*callback)());

#endif /* LIBSET_H */

