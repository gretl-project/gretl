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

typedef enum {
    OPTIM_AUTO,
    OPTIM_BFGS,
    OPTIM_NEWTON,
    OPTIM_MAX
} OptimCode;

typedef enum {
    STEPLEN_POWER,
    STEPLEN_QUAD,
    STEPLEN_MAX
} SteplenCode;

/* guard against consequences of typos */

#define BFGS_MAXITER     "bfgs_maxiter"
#define BFGS_TOLER       "bfgs_toler"
#define BFGS_MAXGRAD     "bfgs_maxgrad"
#define BFGS_VERBSKIP    "bfgs_verbskip"
#define BFGS_RSTEP       "bfgs_richardson"
#define OPTIM_STEPLEN    "optim_steplen"
#define BHHH_MAXITER     "bhhh_maxiter"
#define BHHH_TOLER       "bhhh_toler"
#define LBFGS_MEM        "lbfgs_mem"
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
#define USE_LBFGS        "lbfgs"
#define LOOP_MAXITER     "loop_maxiter"
#define RQ_MAXITER       "rq_maxiter"
#define MAX_VERBOSE      "max_verbose"
#define NLS_TOLER        "nls_toler"
#define PCSE             "pcse"
#define PREWHITEN        "hac_prewhiten"
#define QS_BANDWIDTH     "qs_bandwidth"
#define SHELL_OK         "shell_ok"
#define USE_CWD          "use_cwd"
#define USE_SVD          "svd"
#define USE_FCP          "fcp"
#define VECM_NORM        "vecm_norm"
#define GRETL_OPTIM      "optimizer"
#define ARMA_VCV         "arma_vcv"
#define SKIP_MISSING     "skip_missing"
#define R_FUNCTIONS      "R_functions"
#define R_LIB            "R_lib"
#define NORMAL_RAND      "normal_rand"
#define DPDSTYLE         "dpdstyle"
#define GMM_MAXITER      "gmm_maxiter"
#define CSV_DIGITS       "csv_digits"
#define NADARWAT_TRIM    "nadarwat_trim"
#define CONV_HUGE        "huge"
#define USE_OPENMP       "openmp"
#define FDJAC_QUAL       "fdjac_quality"

typedef int (*ITER_PRINT_FUNC) (int, PRN *);
typedef void (*SHOW_ACTIVITY_FUNC) (void);
typedef int (*DEBUG_READLINE) (void *);
typedef int (*DEBUG_OUTPUT) (void *);
typedef int (*QUERY_STOP) (void);

#define set_nls_toler(x) (libset_set_double(NLS_TOLER, x))

int libset_init (void);
void libset_cleanup (void);
int libset_restore_state_zero (DATASET *dset);

int push_program_state (void);
int pop_program_state (void);

int libset_get_bool (const char *key);
int libset_set_bool (const char *key, int val);

double libset_get_double (const char *key);
int libset_set_double (const char *key, double val);

double libset_get_user_tolerance (const char *key);

int libset_get_int (const char *key);
int libset_set_int (const char *key, int val);

int is_libset_var (const char *s);

int libset_use_openmp (guint64 n);

int gretl_n_processors (void);

int get_omp_n_threads (void);

/* GUI setter functions */
void set_xsect_hccme (const char *s);
void set_tseries_hccme (const char *s);
void set_panel_hccme (const char *s);
void set_garch_robust_vcv (const char *s);

int get_hac_lag (int T);
int data_based_hac_bandwidth (void);

int get_bkbp_k (const DATASET *dset);
void get_bkbp_periods (const DATASET *dset, int *l, int *u);

void set_mp_bits (int b);
int get_mp_bits (void);

const gretl_matrix *get_init_vals (void);
int n_init_vals (void);
void free_init_vals (void);

const gretl_matrix *get_matrix_mask (void);
int get_matrix_mask_nobs (void);
void destroy_matrix_mask (void);

void set_loop_on (int quiet, int progressive);
void set_loop_off (void);
int gretl_looping (void);
int gretl_looping_currently (void);
int gretl_looping_quietly (void);
int gretl_looping_progressive (void);

void gretl_set_batch_mode (int b);
int gretl_in_batch_mode (void);

void gretl_set_gui_mode (void);
int gretl_in_gui_mode (void);

void set_gretl_echo (int e);
int gretl_echo_on (void);

void set_gretl_messages (int e);
int gretl_messages_on (void);

int gretl_warnings_on (void);
int gretl_debugging_on (void);

void shelldir_init (const char *s);
char *get_shelldir (void);

void set_data_export_decimal_comma (int s);
char get_data_export_decpoint (void);

void set_data_export_delimiter (char c);
char get_data_export_delimiter (void);

const char *get_csv_na_write_string (void);
int set_csv_na_write_string (const char *s);

const char *get_csv_na_read_string (void);
int set_csv_na_read_string (const char *s);

int execute_set_line (const char *line, DATASET *dset, 
		      gretlopt opt, PRN *prn);

void set_iter_print_func (ITER_PRINT_FUNC func);
int iter_print_callback (int i, PRN *prn);
int iter_print_func_installed (void);

void set_show_activity_func (SHOW_ACTIVITY_FUNC func);
void show_activity_callback (void);
int show_activity_func_installed (void);

void set_debug_read_func (DEBUG_READLINE dfunc);
DEBUG_READLINE get_debug_read_func (void);

void set_debug_output_func (DEBUG_OUTPUT dout);
DEBUG_OUTPUT get_debug_output_func (void);

void set_query_stop_func (QUERY_STOP query);
int check_for_stop (void);

void set_workdir_callback (int (*callback)());

int libset_write_script (const char *fname);
int libset_read_script (const char *fname);

#endif /* LIBSET_H */

