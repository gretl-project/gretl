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

#ifdef  __cplusplus
extern "C" {
#endif

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
    HAC_REFUSE,
    HAC_ES,
    HAC_AM
} HACopts;

typedef enum {
    ARELLANO,
    BECK_KATZ,
    DRISCOLL_KRAAY
} PanelRobust;

typedef enum {
    USE_CWD         = 1 << 0,  /* store: use current dir as default */
    ECHO_ON         = 1 << 1,  /* echoing commands or not */
    MSGS_ON         = 1 << 2,  /* emitting non-error messages or not */
    FORCE_DECPOINT  = 1 << 3,  /* override locale decimal character */
    USE_SVD         = 1 << 4,  /* SVD decomposition is matrix OLS default */
    PREWHITEN       = 1 << 5,  /* HAC pre-whitening? */
    FORCE_HC        = 1 << 6,  /* don't use HAC for time series */
    USE_LBFGS       = 1 << 7,  /* prefer LBFGS to BFGS? */
    SHELL_OK        = 1 << 8,  /* "shell" facility is approved? */
    WARNINGS        = 1 << 9,  /* print numerical warning messages */
    SKIP_MISSING    = 1 << 10, /* skip NAs when building matrix from series */
    BFGS_RSTEP      = 1 << 11, /* use Richardson in BFGS numerical gradient */
    ROBUST_Z        = 1 << 12, /* use z- not t-score with HCCM/HAC */
    MWRITE_G        = 1 << 13, /* use %g format with mwrite() */
    MPI_USE_SMT     = 1 << 14, /* MPI: use hyperthreads by default */
    STATE_FLAG_MAX  = 1 << 15, /* separator */
    /* state small int (but non-boolean) vars */
    GRETL_OPTIM,
    VECM_NORM,
    GARCH_VCV,
    GARCH_ALT_VCV,
    ARMA_VCV,
    WILDBOOT_DIST,
    FDJAC_QUAL,
    USE_QR,
    MAX_VERBOSE,
    HC_VERSION,
    PANEL_ROBUST,
    HAC_KERNEL,
    HAC_LAG,
    USER_HAC_LAG,
    LBFGS_MEM,
    QUANTILE_TYPE,
    STATE_SMALL_INT_MAX, /* separator: start state int vars */
    HORIZON,
    BOOTREP,
    LOOP_MAXITER,
    BFGS_MAXITER,
    BFGS_VERBSKIP,
    BOOT_ITERS,
    BHHH_MAXITER,
    RQ_MAXITER,
    GMM_MAXITER,
    STATE_INT_MAX, /* separator: end state int vars */
    CONV_HUGE,
    NLS_TOLER,
    BFGS_TOLER,
    BFGS_MAXGRAD,
    BHHH_TOLER,
    QS_BANDWIDTH,
    NADARWAT_TRIM,
    STATE_FLOAT_MAX, /* separator: end state floats */
    CSV_WRITE_NA,
    CSV_READ_NA,
    INITVALS,
    INITCURV,
    MATMASK,
    STATE_VARS_MAX, /* separator */
    /* non-state vars follow */
    GRETL_DEBUG,
    GRETL_ASSERT,
    DATACOLS,
    PLOT_COLLECT,
    R_FUNCTIONS,
    R_LIB,
    LOGLEVEL,
    LOGSTAMP,
    CSV_DIGITS,
    HAC_MISSVALS,
    NS_SMALL_INT_MAX, /* separator */
    GMP_BITS,
    NS_MAX, /* separator */
    BLAS_MNK_MIN,
    OMP_MNK_MIN,
    OMP_N_THREADS,
    SIMD_K_MAX,
    SIMD_MN_MIN,
    USE_DCMT,
    NS_INT_MAX, /* separator */
    SEED,
    CSV_DELIM,
    STOPWATCH,
    VERBOSE,
    SV_WORKDIR,
    SV_LOGFILE,
    GRAPH_THEME,
    DISP_DIGITS,
    TEX_PLOT_OPTS,
    SETVAR_MAX /* sentinel */
} SetKey;

typedef void (*SHOW_ACTIVITY_FUNC) (void);

#define set_nls_toler(x) (libset_set_double(NLS_TOLER, x))

int libset_init (void);
void libset_cleanup (void);

int push_program_state (void);
int pop_program_state (void);

int libset_get_bool (SetKey key);
int libset_set_bool (SetKey key, int val);

double libset_get_double (SetKey key);
int libset_set_double (SetKey key, double val);

double libset_get_user_tolerance (SetKey key);

int libset_get_int (SetKey key);
int libset_set_int (SetKey key, int val);

int libset_help_available (const char *s);

/* GUI setter functions */
void set_xsect_hccme (const char *s);
void set_tseries_hccme (const char *s);
void set_panel_hccme (const char *s);
void set_garch_alt_vcv (const char *s);

int get_hac_lag (int T);
int data_based_hac_bandwidth (void);

int get_bkbp_k (const DATASET *dset);
void get_bkbp_periods (const DATASET *dset, int *l, int *u);

/* convenience accessor */
int get_mp_bits (void);

gretl_matrix *get_initvals (void);
int n_initvals (void);

gretl_matrix *get_initcurv (void);
int n_initcurv (void);

const gretl_matrix *get_matrix_mask (void);
int get_matrix_mask_nobs (void);
void destroy_matrix_mask (void);

void set_loop_on (void);
void set_loop_off (void);

int gretl_looping (void);
int gretl_looping_currently (void);

void gretl_iteration_push (void);
void gretl_iteration_pop (void);
int gretl_iteration_depth (void);

int gretl_iterating (void);

void gretl_set_batch_mode (int b);
int gretl_in_batch_mode (void);

void gretl_set_gui_mode (void);
int gretl_in_gui_mode (void);

void gretl_set_tool_mode (void);
int gretl_in_tool_mode (void);

void set_gretl_echo (int e);
int gretl_echo_on (void);

void set_gretl_messages (int e);
int gretl_messages_on (void);

int gretl_comments_on (void);

int gretl_warnings_on (void);
int gretl_debugging_on (void);

void set_data_export_decimal_comma (int s);
char get_data_export_decpoint (void);

void set_data_export_delimiter (char c);
char get_data_export_delimiter (void);

const char *get_csv_na_write_string (void);
int set_csv_na_write_string (const char *s);

const char *get_csv_na_read_string (void);
int set_csv_na_read_string (const char *s);

const char *get_tex_plot_opts (void);

int execute_set (const char *setobj, const char *setarg,
		 DATASET *dset, gretlopt opt, PRN *prn);

void set_show_activity_func (SHOW_ACTIVITY_FUNC func);
void show_activity_callback (void);
int show_activity_func_installed (void);

void set_user_stop (int s);
int get_user_stop (void);

void set_workdir_callback (int (*callback)(char *s));

int libset_write_script (const char *fname);
int libset_read_script (const char *fname);

#ifdef  __cplusplus
}
#endif

#endif /* LIBSET_H */
