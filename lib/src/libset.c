/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* libset.c for gretl */

#include "libgretl.h"
#include "libset.h"
#include "usermat.h"

#include <unistd.h>
#include <errno.h>

#define PDEBUG 0

enum {
    AUTO_LAG_STOCK_WATSON,
    AUTO_LAG_WOOLDRIDGE,
    AUTO_LAG_NEWEYWEST
};

enum {
    STATE_USE_CWD        = 1 << 0,  /* store: use current dir as default */
    STATE_ECHO_ON        = 1 << 1,  /* echoing commands or not */
    STATE_MSGS_ON        = 1 << 2,  /* emitting non-error messages or not */
    STATE_FORCE_DECPOINT = 1 << 3   /* override locale decimal separator */
};    

/* for values that really want a non-negative integer */
#define UNSET_INT -1
#define is_unset(i) (i == UNSET_INT)

typedef struct set_vars_ set_vars;

struct robust_opts {
    int auto_lag;
    int user_lag;
    int hc_version;
    int force_hc;
    int hkern;
    int prewhite;
    int pcse;
    double qsband;
};

struct garch_opts {
    int vcv_variant;
    int robust_vcv_variant;
};

struct bkbp_opts {
    int k;
    int periods[2];
};

struct max_opts {
    double toler;
    int maxiter;
};

struct set_vars_ {
    int flags;
    int use_qr;                 /* use QR decomposition by default? */
    int use_lbfgs;              /* use LBFGS instead of plain BFGS */
    int halt_on_err;            /* halt cli program on script error? */
    unsigned int seed;          /* for PRNG */
    int shell_ok;               /* shell commands permitted? */
    double hp_lambda;           /* for Hodrick-Prescott filter */
    int horizon;                /* for VAR impulse responses */ 
    int bootrep;                /* bootstrap replications */
    double nls_toler;           /* NLS convergence criterion */
    char delim;                 /* delimiter for CSV data export */
    int longdigits;             /* digits for printing data in long form */
    int max_verbose;            /* verbose output from maximizer? */
    int vecm_norm;              /* VECM beta normalization */
    gretl_matrix *initvals;     /* for parameter initialization */
    struct robust_opts ropts;   /* robust standard error options */
    struct garch_opts gopts;    /* GARCH covariance matrix */
    struct bkbp_opts bkopts;    /* Baxter-King filter */
    struct max_opts bhhh_opts;  /* options for BHHH maximization */
    struct max_opts bfgs_opts;  /* options for BFGS maximization */
    char shelldir[MAXLEN];      /* working dir for shell commands */
};

/* global state */
set_vars *state;

static void robust_opts_init (struct robust_opts *opts)
{
    opts->auto_lag = AUTO_LAG_STOCK_WATSON;
    opts->user_lag = UNSET_INT;
    opts->hc_version = 0;
    opts->force_hc = 0; 
    opts->hkern = KERNEL_BARTLETT;
    opts->prewhite = 0;
    opts->pcse = 0;
    opts->qsband = NADBL;
}

static void robust_opts_copy (struct robust_opts *opts)
{
    opts->auto_lag = state->ropts.auto_lag;
    opts->user_lag = state->ropts.user_lag;
    opts->hc_version = state->ropts.hc_version;
    opts->force_hc = state->ropts.force_hc; 
    opts->hkern = state->ropts.hkern; 
    opts->prewhite = state->ropts.prewhite;
    opts->pcse = state->ropts.pcse;
    opts->qsband = state->ropts.qsband;
}

static const char *hac_kernel_string (void)
{
    switch (state->ropts.hkern) {
    case KERNEL_BARTLETT:
	return "bartlett";
    case KERNEL_PARZEN:
	return "parzen";
    case KERNEL_QS:
	return "quadratic-spectral";
    default:
	return "";
    }
}

static void garch_opts_init (struct garch_opts *opts)
{
    opts->vcv_variant = VCV_UNSET;
    opts->robust_vcv_variant = VCV_UNSET;
}

static void garch_opts_copy (struct garch_opts *opts)
{
    opts->vcv_variant = state->gopts.vcv_variant;
    opts->robust_vcv_variant = state->gopts.robust_vcv_variant;
}

static void bkbp_opts_init (struct bkbp_opts *opts)
{
    opts->k = UNSET_INT;
    opts->periods[0] = UNSET_INT;
    opts->periods[1] = UNSET_INT;
}

static void bkbp_opts_copy (struct bkbp_opts *opts)
{
    opts->k = state->bkopts.k;
    opts->periods[0] = state->bkopts.periods[0];
    opts->periods[1] = state->bkopts.periods[1];
}

static void max_opts_init (struct max_opts *opts)
{
    opts->toler = NADBL;
    opts->maxiter = 500;
}

static void bhhh_opts_copy (struct max_opts *opts)
{
    opts->toler = state->bhhh_opts.toler;
    opts->maxiter = state->bhhh_opts.maxiter;
}

static void bfgs_opts_copy (struct max_opts *opts)
{
    opts->toler = state->bfgs_opts.toler;
    opts->maxiter = state->bfgs_opts.maxiter;
}

static void print_initvals (const gretl_matrix *ivals, PRN *prn)
{
    if (ivals == NULL) {
	pputs(prn, " initvals = auto\n");
    } else {
	gretl_matrix_print_to_prn(ivals, " initvals =", prn);
    }
}

/* check_for_state() returns non-zero if the program options
   state is not readable */

static int check_for_state (void) 
{
    if (state == NULL) { /* was "!=" ?? */
	return libset_init();
    } else {
#if PDEBUG > 1
	fprintf(stderr, "check_for_state: state = %p\n", (void *) state);
#endif
	return 0;
    }
}

static int flag_to_bool (set_vars *sv, int flag)
{
    return (sv->flags & flag)? 1 : 0;
}

static void state_vars_copy (set_vars *sv)
{
#if PDEBUG
    fprintf(stderr, "state_vars_copy() called\n");
#endif
    sv->flags = state->flags;
    sv->use_qr = state->use_qr;
    sv->use_lbfgs = state->use_lbfgs;
    sv->halt_on_err = state->halt_on_err;
    sv->seed = state->seed;
    sv->shell_ok = state->shell_ok;
    sv->hp_lambda = state->hp_lambda;
    sv->horizon = state->horizon;
    sv->bootrep = state->bootrep;
    sv->nls_toler = state->nls_toler;
    sv->delim = state->delim; 
    sv->longdigits = state->longdigits; 
    sv->max_verbose = state->max_verbose;
    sv->vecm_norm = state->vecm_norm;
    sv->initvals = gretl_matrix_copy(state->initvals);
    strcpy(sv->shelldir, state->shelldir);

    robust_opts_copy(&sv->ropts);
    garch_opts_copy(&sv->gopts);
    bkbp_opts_copy(&sv->bkopts);
    bhhh_opts_copy(&sv->bhhh_opts);
    bfgs_opts_copy(&sv->bfgs_opts);
}

static void state_vars_init (set_vars *sv)
{
#if PDEBUG
    fprintf(stderr, "state_vars_init called\n");
#endif
    sv->flags = STATE_ECHO_ON | STATE_MSGS_ON;
    sv->use_qr = UNSET_INT; 
    sv->use_lbfgs = UNSET_INT; 
    sv->seed = 0;
    sv->halt_on_err = UNSET_INT;
    sv->shell_ok = 0;
    sv->hp_lambda = NADBL;
    sv->horizon = UNSET_INT;
    sv->bootrep = 1000;
    sv->nls_toler = NADBL;
    sv->delim = UNSET_INT;
    sv->longdigits = 10;
    sv->max_verbose = 0;
    sv->vecm_norm = NORM_PHILLIPS;
    sv->initvals = NULL;
    *sv->shelldir = '\0';

    robust_opts_init(&sv->ropts);
    garch_opts_init(&sv->gopts);
    bkbp_opts_init(&sv->bkopts);
    max_opts_init(&sv->bhhh_opts);
    max_opts_init(&sv->bfgs_opts);
}

int get_hc_version (void)
{
    if (check_for_state()) {
	return 0;
    } else {
	return state->ropts.hc_version;
    }
}

double get_hp_lambda (void)
{
    if (check_for_state()) {
	return 0.0;
    } else {
	return state->hp_lambda;
    }
}

int set_hp_lambda (double d)
{
    if (check_for_state()) {
	return E_ALLOC;
    }

    if (d > 0) {
	state->hp_lambda = d;
	return 0;
    } else {
	return 1;
    }
}

int get_bkbp_k (const DATAINFO *pdinfo)
{
    if (check_for_state()) {
	return 0;
    }

    if (is_unset(state->bkopts.k)) {
	if (pdinfo->pd == 1) {
	    return 3;
	} else if (pdinfo->pd == 4) {
	    return 12;
	} else if (pdinfo->pd == 12) {
	    return 36;
	} else {
	    return 3;
	}
    } else {
	return state->bkopts.k;
    }
}

int set_bkbp_k (int k)
{
    if (check_for_state()) {
	return E_ALLOC;
    }

    if (k > 0) {
	state->bkopts.k = k;
	return 0;
    } else {
	return 1;
    }
}

void unset_bkbp_k (void)
{
    if (check_for_state()) {
	return;
    }

    state->bkopts.k = UNSET_INT;
}

void get_bkbp_periods (const DATAINFO *pdinfo, int *l, int *u)
{
    if (check_for_state()) {
	return;
    }

    if (is_unset(state->bkopts.periods[0])) {
	*l = (pdinfo->pd == 4)? 6 :
	    (pdinfo->pd == 12)? 18 : 2;
    } else {
	*l = state->bkopts.periods[0];
    }

    if (is_unset(state->bkopts.periods[1])) {
	*u = (pdinfo->pd == 4)? 32 :
	    (pdinfo->pd == 12)? 96 : 8;
    } else {
	*u = state->bkopts.periods[1];
    }
}

int set_bkbp_periods (int l, int u)
{
    if (check_for_state()) {
	return E_ALLOC;
    }

    if (l > 0 && u > l) {
	state->bkopts.periods[0] = l;
	state->bkopts.periods[1] = u;
	return 0;
    } else {
	return 1;
    }
}

void unset_bkbp_periods (void)
{
    if (check_for_state()) {
	return;
    }

    state->bkopts.periods[0] = UNSET_INT;
    state->bkopts.periods[1] = UNSET_INT;
}

double get_bhhh_toler (void)
{
    if (check_for_state()) {
	return 1.0;
    }

    return state->bhhh_opts.toler;
}

int get_bhhh_maxiter (void)
{
    if (check_for_state()) {
	return 0;
    }

    return state->bhhh_opts.maxiter;
}

double get_bfgs_toler (void)
{
    if (check_for_state()) {
	return 1.0;
    }

    if (na(state->bfgs_opts.toler)) {
	state->bfgs_opts.toler = get_default_nls_toler();
    }

    return state->bfgs_opts.toler;
}

int get_bfgs_maxiter (void)
{
    if (check_for_state()) {
	return 0;
    }

    return state->bfgs_opts.maxiter;
}

enum {
    BHHH,
    BFGS
};

static int real_set_toler (double tol, int c)
{
    int err = 0;

    if (check_for_state()) {
	return E_ALLOC;
    }

    if (tol <= 0.0) {
	err = 1;
    } else if (c == BFGS) {
	state->bfgs_opts.toler = tol;
    } else {
	state->bhhh_opts.toler = tol;
    }

    return err;
}

int set_bhhh_toler (double tol)
{
    return real_set_toler(tol, BHHH);
}

int set_bfgs_toler (double tol)
{
    return real_set_toler(tol, BFGS);
}

static int real_set_maxiter (int n, int c)
{
    int err = 0;

    if (check_for_state()) {
	return E_ALLOC;
    }

    if (n < 1) {
	err = 1;
    } else if (c == BFGS) {
	state->bfgs_opts.maxiter = n;
    } else {
	state->bhhh_opts.maxiter = n;
    }

    return err;
}

int set_bhhh_maxiter (int n)
{
    return real_set_maxiter(n, BHHH);
}

int set_bfgs_maxiter (int n)
{
    return real_set_maxiter(n, BFGS);
}

int set_long_digits (int n)
{
    if (check_for_state()) {
	return E_ALLOC;
    }

    if (n < 1 || n > 20) {
	return 1;
    } else {
	state->longdigits = n;
	return 0;
    }
}

int set_max_verbose (int n)
{
    if (check_for_state()) {
	return E_ALLOC;
    }

    if (n < 0) {
	return 1;
    } else {
	state->max_verbose = n;
	return 0;
    }
}

int get_VAR_horizon (void)
{
    if (check_for_state()) {
	return 0;
    } else {
	return state->horizon;
    }
}

int get_bootstrap_replications (void)
{
    if (check_for_state()) {
	return 0;
    } else {
	return state->bootrep;
    }
}

double get_nls_toler (void)
{
    if (check_for_state()) {
	return 1.0;
    }

    if (na(state->nls_toler)) {
	state->nls_toler = get_default_nls_toler();
    }

    return state->nls_toler;
}

int set_nls_toler (double tol)
{
    int err = 0;

    if (check_for_state()) {
	return E_ALLOC;
    }

    if (tol <= 0.0) {
	err = 1;
    } else {
	state->nls_toler = tol;
    }

    return err;
}

static int get_or_set_force_hc (int f)
{
    if (check_for_state()) {
	return 0;
    }

    if (f >= 0) {
	state->ropts.force_hc = f;
    }

    return state->ropts.force_hc;
}

static int get_or_set_garch_vcv (int v)
{
    if (check_for_state()) {
	return 0;
    }

    if (v >= 0) {
	state->gopts.vcv_variant = v;
    }

    return state->gopts.vcv_variant;
}

static int get_or_set_garch_robust_vcv (int v)
{
    if (check_for_state()) {
	return 0;
    }

    if (v >= 0) {
	state->gopts.robust_vcv_variant = v;
    }

    return state->gopts.robust_vcv_variant;
}

struct g_vcv {
    int opt;
    char *str;
};

static struct g_vcv g_vcv_types[] = {
    { VCV_UNSET,   "unset" },
    { VCV_HESSIAN, "hessian" },
    { VCV_IM,      "im" },
    { VCV_OP,      "op" },
    { VCV_QML,     "qml" },
    { VCV_BW,      "bw" }
};

static const char *garch_vcv_string (void)
{
    int i, n = sizeof g_vcv_types / sizeof g_vcv_types[0];
    const char *ret = "";

    for (i=0; i<n; i++ ) {
	if (state->gopts.vcv_variant == g_vcv_types[i].opt) {
	    ret = g_vcv_types[i].str;
	    break;
	}
    }

    return ret;
}

static int set_garch_vcv_variant (const char *s)
{
    int i, vopt = VCV_UNSET;
    int n = sizeof g_vcv_types / sizeof g_vcv_types[0];

    for (i=0; i<n; i++ ) {
	if (!strcmp(s, g_vcv_types[i].str)) {
	    vopt = g_vcv_types[i].opt;
	    break;
	}
    }	

    get_or_set_garch_vcv(vopt);
    
    return (vopt == VCV_UNSET);
}

int get_garch_vcv_version (void)
{
    return get_or_set_garch_vcv(-1);
}

int get_garch_robust_vcv_version (void)
{
    return get_or_set_garch_robust_vcv(-1);
}

static void set_force_hc (int f)
{
    get_or_set_force_hc(f);
}

int get_force_hc (void)
{
    return get_or_set_force_hc(-1);
}

void set_gretl_echo (int e)
{
    if (check_for_state()) return;

    if (e) {
	state->flags |= STATE_ECHO_ON;
    } else {
	state->flags &= ~STATE_ECHO_ON;
    }
}

int gretl_echo_on (void)
{
    if (check_for_state()) return 1;
    return flag_to_bool(state, STATE_ECHO_ON);
}

void set_gretl_messages (int e)
{
    if (check_for_state()) return;

    if (e) {
	state->flags |= STATE_MSGS_ON;
    } else {
	state->flags &= ~STATE_MSGS_ON;
    }
}

int gretl_messages_on (void)
{
    if (check_for_state()) return 1;
    return flag_to_bool(state, STATE_MSGS_ON);
}

char get_csv_delim (const DATAINFO *pdinfo)
{
    check_for_state();
    if (state->delim > 0) {
	return state->delim;
    } else {
	return pdinfo->delim;
    }
}

int get_long_digits (void)
{
    check_for_state();
    return state->longdigits;
}

int get_max_verbose (void)
{
    check_for_state();
    return state->max_verbose;
}

int get_vecm_norm (void)
{
    check_for_state();
    return state->vecm_norm;
}

const gretl_matrix *get_init_vals (void)
{
    check_for_state();
    return state->initvals;
}

void free_init_vals (void)
{
    if (state->initvals != NULL) {
	gretl_matrix_free(state->initvals);
	state->initvals = NULL;
    }
}

int n_init_vals (void)
{
    check_for_state();
    if (state->initvals != NULL) {
	return gretl_vector_get_length(state->initvals);
    } else {
	return 0;
    }
}

char *get_shelldir (void)
{
    check_for_state();

    if (state != NULL && *state->shelldir != '\0') {
	return state->shelldir;
    } else {
	return NULL;
    }
} 

int get_hac_lag (int T)
{
    check_for_state();

    /* Variants of Newey-West */

    if (state->ropts.user_lag >= 0 && state->ropts.user_lag < T - 2) {
	/* FIXME upper limit? */
	return state->ropts.user_lag;
    } else if (state->ropts.auto_lag == AUTO_LAG_WOOLDRIDGE) {
	return 4.0 * pow(T / 100.0, 2.0 / 9.0);
    } else {
	/* Stock-Watson default */
	return 0.75 * pow(T, 1.0 / 3.0);
    }
}

/* prewhitening implies nw3, but not vice versa */

int data_based_hac_bandwidth (void)
{
    if (is_unset(state->ropts.user_lag)) {
	if (state->ropts.auto_lag == AUTO_LAG_NEWEYWEST ||
	    state->ropts.prewhite) {
	    return 1;
	}
    }

    return 0;
}

int get_hac_kernel (void)
{
    check_for_state();

    return state->ropts.hkern;
}

void set_hac_kernel (int k)
{
    check_for_state();

    if (k >= KERNEL_BARTLETT && k <= KERNEL_QS) {
	state->ropts.hkern = k;
    }
}

int get_hac_prewhiten (void)
{
    check_for_state();

    return state->ropts.prewhite;
}

void set_hac_prewhiten (int w)
{
    check_for_state();
    state->ropts.prewhite = (w != 0);
}

int get_panel_beck_katz (void)
{
    check_for_state();

    return state->ropts.pcse;
}

void set_panel_beck_katz (int b)
{
    check_for_state();
    state->ropts.pcse = (b != 0);
}

double get_qs_bandwidth (void)
{
    check_for_state();

    if (!na(state->ropts.qsband) && state->ropts.qsband > 0) {
	return state->ropts.qsband;
    } else {
	/* what's a sensible default here? */
	return 2.0;
    }
}

void set_qs_bandwidth (double w)
{
    check_for_state();
    if (w >= 0.0) {
	state->ropts.qsband = w;
    }
}

static const char *get_hac_lag_string (void)
{
    check_for_state();

    if (state->ropts.user_lag >= 0 && state->ropts.user_lag < 1000) {
	static char lagstr[6];

	sprintf(lagstr, "%d", state->ropts.user_lag);
	return lagstr;
    } else if (state->ropts.auto_lag == AUTO_LAG_STOCK_WATSON) {
	return "nw1";
    } else {
	return "nm2";
    }
}

static int 
libset_numeric_string (const char *s, int *pi, double *px, int *err)
{
    char *test;
    int ret = 1;

    if (s == NULL || *s == '\0' ||
	!strcmp(s, "inf") || !strcmp(s, "nan")) {
	return 0;
    }

    errno = 0;

    gretl_push_c_numeric_locale();

    if (px != NULL) {
	*px = strtod(s, &test);
	if (*test != '\0') {
	    ret = 0;
	} else if (errno == ERANGE) {
	    gretl_errmsg_set(strerror(errno));
	    *err = 1;
	}
    } else {
	long li = strtol(s, &test, 10);

	if (*test != '\0') {
	    ret = 0;
	} else if (errno == ERANGE) {
	    gretl_errmsg_set(strerror(errno));
	    *err = 1;
	} else {
	    *pi = (int) li;
	}
    }

    gretl_pop_c_numeric_locale();

    return ret;
}

/* all the libset ints/doubles are non-negative */

static int libset_get_scalar (const char *s, double **Z,
			      DATAINFO *pdinfo,
			      int *pi, double *px)
{
    double x = NADBL;
    int v, err = 0;

    if (libset_numeric_string(s, pi, px, &err)) {
	if (err) {
	    err = E_DATA;
	} else if (pi != NULL && *pi < 0) {
	    err = E_DATA;
	} else if (px != NULL && *px < 0.0) {
	    err = E_DATA;
	}
	return err;
    }

    v = varindex(pdinfo, s);

    if (v >= pdinfo->v) {
	sprintf(gretl_errmsg, "'%s': unrecognized name", s);
	return E_UNKVAR;
    } else if (var_is_series(pdinfo, v)) {
	return E_DATATYPE;
    } else {
	x = Z[v][0];
    }

    if (x < 0.0) {
	return E_DATA;
    }

    if (px != NULL) {
	*px = x;
    } else if (pi != NULL) {
	if (na(x) || fabs(x) > (double) INT_MAX) {
	    err = E_DATA;
	} else {
	    *pi = (int) x;
	}
    }

    return err;
}

static int parse_hc_variant (const char *s)
{
    int err = 1;

    check_for_state();

    if (!strcmp(s, "0") || !strcmp(s, "1") ||
	!strcmp(s, "2") || !strcmp(s, "3")) {
	state->ropts.hc_version = atoi(s);
	err = 0;
    } else if (!strcmp(s, "3a")) {
	state->ropts.hc_version = 4;
	err = 0;
    }

    if (err) {
	int hcv;

	if (!strcmp(s, "hc3a")) {
	    state->ropts.hc_version = 4;
	    err = 0;
	} else if (sscanf(s, "hc%d", &hcv)) {
	    if (hcv >= 0 && hcv <= 4) {
		state->ropts.hc_version = hcv;
		err = 0;
	    }
	}
    }

    return err;
}

void set_xsect_hccme (const char *s)
{
    char *scpy = gretl_strdup(s);

    if (scpy == NULL) return;

    lower(scpy);
    parse_hc_variant(scpy);
    free(scpy);
}

void set_tseries_hccme (const char *s)
{
    char *scpy = gretl_strdup(s);

    if (scpy == NULL) return;

    lower(scpy);
    if (parse_hc_variant(scpy) == 0) {
	set_force_hc(1);
    } else {
	set_force_hc(0);
    }
    free(scpy);
}

void set_panel_hccme (const char *s)
{
    check_for_state();

    if (!strcmp(s, "Arellano")) {
	state->ropts.pcse = 0;
    } else if (!strcmp(s, "PCSE")) {
	state->ropts.pcse = 1;
    }
}

void set_garch_robust_vcv (const char *s)
{
    char *scpy = gretl_strdup(s);

    if (scpy == NULL) return;

    lower(scpy);
    if (!strcmp(s, "qml")) {
	get_or_set_garch_robust_vcv(VCV_QML);
    } else if (!strcmp(s, "bw")) {
	get_or_set_garch_robust_vcv(VCV_BW);
    }
    free(scpy);
}

static int set_line_width (const char *s0, const char *s1,
			   DATAINFO *pdinfo, PRN *prn)
{
    int v, w, err = 0;

    if (!isdigit((unsigned char) *s1)) {
	return 1;
    }

    if (isdigit((unsigned char) *s0)) {
	v = atoi(s0);
    } else {
	v = varindex(pdinfo, s0);
    }

    if (v < 1 || v >= pdinfo->v) {
	return E_DATA;
    }

    w = atoi(s1);

    if (w < 0 || w > 32) {
	err = E_DATA;
    } else {
	var_set_linewidth(pdinfo, v, w);
	pprintf(prn, _("Line width for %s = %d\n"), 
		pdinfo->varname[v], w);
    }

    return err;
}

static int set_bkbp_limits (const char *s0, const char *s1,
			    double **Z, DATAINFO *pdinfo,
			    PRN *prn)
{
    int p0, p1;
    int err = 0;

    err = libset_get_scalar(s0, Z, pdinfo, &p0, NULL);
    if (!err) {
	err = libset_get_scalar(s1, Z, pdinfo, &p1, NULL);
    }

    if (err) {
	return err;
    }

    if (p1 < p0) {
	/* 2nd entry should be bigger than 1st one */
	int tmp = p1;

	p1 = p0;
	p0 = tmp;
    }

    pprintf(prn, _("Baxter-King band = %d-%d periods\n"), p0, p1);

    state->bkopts.periods[0] = p0;
    state->bkopts.periods[1] = p1;

    return 0;
}

static int set_initvals (const char *s, const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *m;
    char mname[VNAMELEN];
    int err = 0;

    /* skip past "set initvals" */
    s += 12;

    if (sscanf(s, "%15s", mname) != 1 || !strcmp(mname, "auto")) {
	gretl_matrix_free(state->initvals);
	state->initvals = NULL;
    } else {
	m = get_matrix_by_name(mname);
	if (m == NULL) {
	    pprintf(prn, _("'%s': no such matrix\n"), mname);
	    err = E_DATA;
	} else {
	    state->initvals = gretl_matrix_copy(m);
	    if (state->initvals == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    return err;
}

void shelldir_init (void)
{
    char *test = getcwd(state->shelldir, MAXLEN);

    if (test == NULL) {
	*state->shelldir = '\0';
    }
}

static int set_shelldir (const char *s)
{
    int err = 0;

    /* skip past "set shelldir" and space */
    s += 12;
    s += strspn(s, " ");

    if (*s == '\0') {
	*state->shelldir = '\0';
    } else if (*s == '"') {
	int len = haschar('"', s + 1);

	if (len <= 0 || len >= MAXLEN) {
	    err = E_PARSE;
	} else {
	    *state->shelldir = '\0';
	    strncat(state->shelldir, s + 1, len);
	}
    } else {
	*state->shelldir = '\0';
	strncat(state->shelldir, s, MAXLEN - 1);
    }

    return err;
}

static int parse_set_plotfile (const char *s)
{
    char *fname;
    int err = 0;

    while (isspace((unsigned char) *s)) {
	s++;
    }

    /* now skip two words, "set" and "plotfile" */
    s += strcspn(s, " ");
    s += strspn(s, " ");
    s += strcspn(s, " ");
    s += strspn(s, " ");

    fname = gretl_strdup(s);
    if (fname != NULL) {
	tailstrip(fname);
	set_gretl_plotfile(fname);
	free(fname);
    } else {
	err = E_ALLOC;
    }

    return err;
}

static char delim_from_arg (const char *s)
{
    char ret = 0;

    if (!strcmp(s, "comma")) {
	ret = ',';
    } else if (!strcmp(s, "space")) {
	ret = ' ';
    } else if (!strcmp(s, "tab")) {
	ret = '\t';
    } else if (!strcmp(s, "semicolon")) {
	ret = ';';
    }

    return ret;
}

static const char *arg_from_delim (char c)
{
    const char *ret = "unset";

    if (c == ',') {
	ret = "comma";
    } else if (c == ' ') {
	ret = "space";
    } else if (c == '\t') {
	ret = "tab";
    } else if (c == ';') {
	ret = "semicolon";
    }

    return ret;
}

static int display_settings (PRN *prn)
{
    double dval;
    unsigned int uval;
    int ival;

    pputs(prn, "Variables that can be set using \"set\" (do \"help set\""
	  " for details):\n");

    pprintf(prn, " echo = %d\n", flag_to_bool(state, STATE_ECHO_ON));
    pprintf(prn, " messages = %d\n", flag_to_bool(state, STATE_MSGS_ON));

    ival = get_use_qr(); /* checks env */
    pprintf(prn, " qr = %d\n", state->use_qr);

    ival = get_use_lbfgs(); /* checks env */
    pprintf(prn, " lbfgs = %d\n", state->use_lbfgs);

    pprintf(prn, " use_cwd = %d\n", flag_to_bool(state, STATE_USE_CWD));
    pprintf(prn, " force_decpoint = %d\n", flag_to_bool(state, STATE_FORCE_DECPOINT));

    uval = get_gretl_random_seed();
    pprintf(prn, " seed = %u\n", uval);

    pprintf(prn, " hac_lag = %s\n", get_hac_lag_string());
    pprintf(prn, " hc_version = %d\n", state->ropts.hc_version);
    pprintf(prn, " force_hc = %d\n", state->ropts.force_hc);
    pprintf(prn, " hac_kernel = %s\n", hac_kernel_string());
    pprintf(prn, " hac_prewhiten = %d\n", state->ropts.prewhite);
    if (na(state->ropts.qsband)) {
	pputs(prn, " qs_bandwidth: auto\n");
    } else {
	pprintf(prn, " qs_bandwidth = %g\n", state->ropts.qsband);
    }
    pprintf(prn, " pcse = %d\n", state->ropts.pcse);

    pprintf(prn, " garch_vcv = %s\n", garch_vcv_string());

    if (na(state->hp_lambda)) {
	pputs(prn, " hp_lambda: auto\n");
    } else {
	pprintf(prn, " hp_lambda = %g\n", state->hp_lambda);
    }

    if (is_unset(state->bkopts.periods[0]) ||
	is_unset(state->bkopts.periods[1])) {
	pputs(prn, " bkbp_limits: auto\n");
    } else {
	pprintf(prn, " bkbp_limits = (%d, %d)\n", state->bkopts.periods[0], 
		state->bkopts.periods[1]);
    }

    if (is_unset(state->bkopts.k)) {
	pputs(prn, " bkbp_k: auto\n");
    } else {
	pprintf(prn, " bkbp_k = %d\n", state->bkopts.k);
    }

    if (is_unset(state->horizon)) {
	pputs(prn, " horizon: auto\n");
    } else {
	pprintf(prn, " horizon = %d\n", state->horizon);
    }

    pprintf(prn, " bootrep = %d\n", state->bootrep);

    pprintf(prn, " nls_toler = %g\n", get_nls_toler());

    dval = get_bhhh_toler();
    if (na(dval)) {
	pputs(prn, " bhhh_toler = default\n");
    } else {
	pprintf(prn, " bhhh_toler = %g\n", dval);
    }
    pprintf(prn, " bhhh_maxiter = %d\n", get_bhhh_maxiter());

    dval = get_bfgs_toler();
    if (na(dval)) {
	pputs(prn, " bfgs_toler = default\n");
    } else {
	pprintf(prn, " bfgs_toler = %g\n", dval);
    }
    pprintf(prn, " bfgs_maxiter = %d\n", get_bfgs_maxiter());

    ival = get_halt_on_error(); /* checks env */
    pprintf(prn, " halt_on_error = %d\n", state->halt_on_err);

    pprintf(prn, " shell_ok = %d\n", get_shell_ok());
    pprintf(prn, " csv_delim = %s\n", arg_from_delim(state->delim));
    pprintf(prn, " longdigits = %d\n", state->longdigits);
    pprintf(prn, " max_verbose = %d\n", state->max_verbose);

    if (state->vecm_norm == NORM_DIAG) {
	pputs(prn, " vecm_norm = diag\n");
    } else {
	pputs(prn, " vecm_norm = phillips\n");
    }
    
    print_initvals(state->initvals, prn);

    if (*state->shelldir) {
	pprintf(prn, " shelldir = '%s'\n", state->shelldir);
    } else {
	pputs(prn, " shelldir = unset\n");
    }

    return 0;
}

#define boolean_on(s) (!strcmp(s, "on") || !strcmp(s, "1") || \
                       !strcmp(s, "true"))

#define boolean_off(s) (!strcmp(s, "off") || !strcmp(s, "0") || \
                        !strcmp(s, "false"))

int execute_set_line (const char *line, double **Z, DATAINFO *pdinfo, 
		      PRN *prn)
{
    char setobj[16], setarg[16], setarg2[16];
    int k, nw, err = E_PARSE;
    double x;

    check_for_state();

    *setobj = *setarg = *setarg2 = '\0';

    nw = sscanf(line, "%*s %15s %15s %15s", setobj, setarg, setarg2);

    if (nw <= 0) {
	return display_settings(prn);
    }

    /* specials which need the whole line */
    if (nw > 1) {
	if (!strcmp(setobj, "plotfile")) {
	    return parse_set_plotfile(line);
	} else if (!strcmp(setobj, "initvals")) {
	    return set_initvals(line, pdinfo, prn);
	} else if (!strcmp(setobj, "shelldir")) {
	    return set_shelldir(line);
	}
    }

    if (nw == 1) {
	if (!strcmp(setobj, "echo")) {
	    state->flags |= STATE_ECHO_ON;
	    err = 0;
	} else if (!strcmp(setobj, "stopwatch")) {
	    gretl_stopwatch();
	    err = 0;
	}
    } else if (nw == 2) {
	lower(setarg);

	/* set command echo on/off */
	if (!strcmp(setobj, "echo")) {
	    if (boolean_off(setarg)) {
		state->flags &= ~STATE_ECHO_ON;
		err = 0;
	    } else if (boolean_on(setarg)) {
		state->flags |= STATE_ECHO_ON;
		err = 0;
	    }
	} else if (!strcmp(setobj, "messages")) {
	    if (boolean_off(setarg)) {
		state->flags &= ~STATE_MSGS_ON;
		err = 0;
	    } else if (boolean_on(setarg)) {
		state->flags |= STATE_MSGS_ON;
		err = 0;
	    }
	} else if (!strcmp(setobj, "csv_delim")) {
	    char c = delim_from_arg(setarg);

	    if (c > 0) {
		state->delim = c;
		err = 0;
	    }
	} else if (!strcmp(setobj, "hac_lag")) {
	    /* set max lag for HAC estimation */
	    if (!strcmp(setarg, "nw1")) {
		state->ropts.auto_lag = AUTO_LAG_STOCK_WATSON;
		state->ropts.user_lag = UNSET_INT;
		err = 0;
	    } else if (!strcmp(setarg, "nw2")) {
		state->ropts.auto_lag = AUTO_LAG_WOOLDRIDGE;
		state->ropts.user_lag = UNSET_INT;
		err = 0;
	    } else if (!strcmp(setarg, "nw3") ||
		       !strcmp(setarg, "auto")) {
		state->ropts.auto_lag = AUTO_LAG_NEWEYWEST;
		state->ropts.user_lag = UNSET_INT;
		err = 0;
	    } else if (isdigit(*setarg)) {
		state->ropts.user_lag = atoi(setarg);
		err = 0;
	    }
	} else if (!strcmp(setobj, "hc_version")) {
	    /* set HCCM variant */
	    err = parse_hc_variant(setarg);
	} else if (!strcmp(setobj, "force_hc")) {
	    /* use HCCM, not HAC, even for time series */
	    if (boolean_on(setarg)) { 
		set_force_hc(1);
		err = 0;
	    } else if (boolean_off(setarg)) { 
		set_force_hc(0);
		err = 0;
	    }
	} else if (!strcmp(setobj, "hac_kernel")) {
	    if (!strcmp(setarg, "bartlett")) {
		state->ropts.hkern = KERNEL_BARTLETT;
		err = 0;
	    } else if (!strcmp(setarg, "parzen")) {
		state->ropts.hkern = KERNEL_PARZEN;
		err = 0;
	    } else if (!strcmp(setarg, "qs")) {
		state->ropts.hkern = KERNEL_QS;
		err = 0;
	    }
	} else if (!strcmp(setobj, "hac_prewhiten")) {
	    if (boolean_on(setarg)) { 
		set_hac_prewhiten(1);
		err = 0;
	    } else if (boolean_off(setarg)) { 
		set_hac_prewhiten(0);
		err = 0;
	    }
	} else if (!strcmp(setobj, "vecm_norm")) {
	    if (!strcmp(setarg, "phillips")) {
		state->vecm_norm = NORM_PHILLIPS;
		err = 0;
	    } else if (!strcmp(setarg, "diag")) {
		state->vecm_norm = NORM_DIAG;
		err = 0;
	    }
	} else if (!strcmp(setobj, "pcse")) {
	    if (boolean_on(setarg)) { 
		set_panel_beck_katz(1);
		err = 0;
	    } else if (boolean_off(setarg)) { 
		set_panel_beck_katz(0);
		err = 0;
	    }	    
	} else if (!strcmp(setobj, "qs_bandwidth")) {
	    err = libset_get_scalar(setarg, Z, pdinfo, NULL, &x);
	    if (!err) {
		state->ropts.qsband = x;
	    }	    
	} else if (!strcmp(setobj, "garch_vcv")) {
	    /* set GARCH VCV variant */
	    err = set_garch_vcv_variant(setarg);
	} else if (!strcmp(setobj, "qr")) {
	    /* switch QR vs Cholesky decomposition */
	    if (boolean_on(setarg)) {
		state->use_qr = 1;
		err = 0;
	    } else if (boolean_off(setarg)) {
		state->use_qr = 0;
		err = 0;
	    }
	} else if (!strcmp(setobj, "lbfgs")) {
	    /* switch LBFGS vs plain BFGS */
	    if (boolean_on(setarg)) {
		state->use_lbfgs = 1;
		err = 0;
	    } else if (boolean_off(setarg)) {
		state->use_lbfgs = 0;
		err = 0;
	    }
	} else if (!strcmp(setobj, "use_cwd")) {
	    if (boolean_on(setarg)) {
		state->flags |= STATE_USE_CWD;
		err = 0;
	    } else if (boolean_off(setarg)) {
		state->flags &= ~STATE_USE_CWD;
		err = 0;
	    }
	} else if (!strcmp(setobj, "force_decpoint")) {
	    if (boolean_on(setarg)) {
		state->flags |= STATE_FORCE_DECPOINT;
		err = 0;
	    } else if (boolean_off(setarg)) {
		state->flags &= ~STATE_FORCE_DECPOINT;
		err = 0;
	    }
	} else if (!strcmp(setobj, "halt_on_error")) {
	    if (boolean_on(setarg)) {
		state->halt_on_err = 1;
		err = 0;
	    } else if (boolean_off(setarg)) {
		state->halt_on_err = 0;
		err = 0;
	    }
	} else if (!strcmp(setobj, "shell_ok")) {
	    pprintf(prn, "You can only set this variable via the gretl GUI\n");
	} else if (!strcmp(setobj, "seed")) {
	    /* seed for PRNG */
	    err = libset_get_scalar(setarg, Z, pdinfo, &k, NULL);
	    if (!err) {
		gretl_rand_set_seed((unsigned int) k);
		pprintf(prn, 
			_("Pseudo-random number generator seeded with %d\n"), k);
		state->seed = k;
	    }
	} else if (!strcmp(setobj, "hp_lambda")) {
	    /* Hodrick-Prescott filter parameter */
	    if (!strcmp(setarg, "auto")) {
		state->hp_lambda = NADBL;
		err = 0;
	    } else {
		err = libset_get_scalar(setarg, Z, pdinfo, NULL, &x);
		if (!err) {
		    state->hp_lambda = x;
		}
	    }
	} else if (!strcmp(setobj, "bkbp_k")) {
	    /* Baxter-King approximation order */
	    err = libset_get_scalar(setarg, Z, pdinfo, &k, NULL);
	    if (!err) {
		state->bkopts.k = k;
		pprintf(prn, 
			_("Baxter-King approximation = %d\n"), state->bkopts.k);
	    }
	} else if (!strcmp(setobj, "horizon")) {
	    /* horizon for VAR impulse responses */
	    if (!strcmp(setarg, "auto")) {
		state->horizon = UNSET_INT;
		err = 0;
	    } else {
		err = libset_get_scalar(setarg, Z, pdinfo, &k, NULL);
		if (!err) {
		    state->horizon = k;
		} else {
		    state->horizon = UNSET_INT;
		}
	    }
	} else if (!strcmp(setobj, "bootrep")) {
	    err = libset_get_scalar(setarg, Z, pdinfo, &k, NULL);
	    if (!err) {
		err = (k <= 0);
	    }
	    if (!err) {
		state->bootrep = k;
	    } 
	} else if (!strcmp(setobj, "nls_toler")) {
	    err = libset_get_scalar(setarg, Z, pdinfo, NULL, &x);
	    if (!err) {
		err = set_nls_toler(x);
	    }
	} else if (!strcmp(setobj, "bhhh_toler")) {
	    /* Tolerance for BHHH (ARMA, Tobit) */
	    if (!strcmp(setarg, "default")) {
		set_bhhh_toler(NADBL);
		err = 0;
	    } else {
		err = libset_get_scalar(setarg, Z, pdinfo, NULL, &x);
		if (!err) {
		    err = set_bhhh_toler(x);
		}
	    }
	} else if (!strcmp(setobj, "bhhh_maxiter")) {
	    /* Maximum iterations for BHHH (ARMA, Tobit) */
	    if (isdigit(*setarg)) {
		err = set_bhhh_maxiter(atoi(setarg));
	    }
	} else if (!strcmp(setobj, "bfgs_toler")) {
	    if (!strcmp(setarg, "default")) {
		set_bfgs_toler(NADBL);
		err = 0;
	    } else {
		err = libset_get_scalar(setarg, Z, pdinfo, NULL, &x);
		if (!err) {
		    err = set_bfgs_toler(x);
		}
	    }
	} else if (!strcmp(setobj, "bfgs_maxiter")) {
	    err = libset_get_scalar(setarg, Z, pdinfo, &k, NULL);
	    if (!err) {
		err = set_bfgs_maxiter(k);
	    }
	} else if (!strcmp(setobj, "longdigits")) {
	    err = libset_get_scalar(setarg, Z, pdinfo, &k, NULL);
	    if (!err) {
		err = set_long_digits(k);
	    }
	} else if (!strcmp(setobj, "max_verbose")) {
	    err = libset_get_scalar(setarg, Z, pdinfo, &k, NULL);
	    if (!err) {
		err = set_max_verbose(k);
	    }
	} 	    
    } else if (nw == 3) {
	if (!strcmp(setobj, "bkbp_limits")) {
	    err = set_bkbp_limits(setarg, setarg2, Z, pdinfo, prn);
	} else if (!strcmp(setobj, "linewidth")) {
	    err = set_line_width(setarg, setarg2, pdinfo, prn);
	}
    }
		    
    return err;
}

/* use Cholesky or QR for regression? */

void set_use_qr (int set)
{
    check_for_state();

    if (state != NULL) {
	state->use_qr = set;
    }
}

int get_use_qr (void)
{
    if (check_for_state()) {
	return 0;
    }

#if PDEBUG
    fprintf(stderr, "get_use_qr: state = %p\n", (void *) state);
#endif

    if (is_unset(state->use_qr)) {
	char *s = getenv("GRETL_USE_QR");

	if (s != NULL && *s != '\0' && *s != '0') {
	    state->use_qr = 1;
	} else {
	    state->use_qr = 0;
	}
    } 

    return state->use_qr;
}

/* use limited-memory or plain BFGS? */

void set_use_lbfgs (int set)
{
    check_for_state();

    if (state != NULL) {
	state->use_lbfgs = set;
    }
}

int get_use_lbfgs (void)
{
    if (check_for_state()) {
	return 0;
    }

    if (is_unset(state->use_lbfgs)) {
	char *s = getenv("GRETL_USE_LBFGS");

	if (s != NULL && *s != '\0' && *s != '0') {
	    state->use_lbfgs = 1;
	} else {
	    state->use_lbfgs = 0;
	}
    } 

    return state->use_lbfgs;
}

void set_use_cwd (int set)
{
    check_for_state();

    if (state != NULL) {
	if (set) {
	    state->flags |= STATE_USE_CWD;
	} else {
	    state->flags &= ~STATE_USE_CWD;
	}
    }
}

int get_use_cwd (void)
{
    if (check_for_state()) {
	return 0;
    }

    return flag_to_bool(state, STATE_USE_CWD);
}

int get_halt_on_error (void)
{
    if (check_for_state()) {
	return 1;
    }

#if PDEBUG
    fprintf(stderr, "get_halt_on_error: state = %p, halt_on_err = %d\n",
	    (void *) state, state->halt_on_err);
#endif

    if (is_unset(state->halt_on_err)) {
	char *s = getenv("GRETL_KEEP_GOING");

	if (s != NULL && *s != '\0' && *s != '0') {
	    state->halt_on_err = 0;
	} else {
	    state->halt_on_err = 1;
	}
    } 

    return state->halt_on_err;
}

#ifndef WIN32
static int read_cli_shell_status (void)
{
    char shellstamp[FILENAME_MAX];
    FILE *fp;
    int ok = 0;

    sprintf(shellstamp, "%s.gretl_shell_stamp", gretl_user_dir());
    fp = fopen(shellstamp, "r");
    if (fp != NULL) {
	ok = 1;
	fclose(fp);
    }

    return ok;
}
#endif

void set_shell_ok (int set)
{
    check_for_state();
    state->shell_ok = set;
}

int get_shell_ok (void)
{
    check_for_state();

#ifndef WIN32
    if (!gretl_in_gui_mode()) {
	state->shell_ok = read_cli_shell_status();
    }
#endif

    return state->shell_ok;
}

/* Mechanism for pushing and popping program state for user-defined
   functions. push_program_state() is used when a function starts
   execution: the function gets a copy of the current program state,
   while that state is pushed onto the stack for restoration when the
   function exits.
*/

static int n_states;
static set_vars **state_stack;

int push_program_state (void)
{
    set_vars **sstack;
    set_vars *newstate;
    int ns = n_states;
    int err = 0;

#if PDEBUG
    fprintf(stderr, "push_program_state: n_states = %d\n", ns);
#endif

    newstate = malloc(sizeof *newstate);

    if (newstate == NULL) {
	err = 1;
    } else {
	sstack = realloc(state_stack, (ns + 1) * sizeof *sstack);
	if (sstack == NULL) {
	    free(newstate);
	    err = 1;
	}
    }

    if (!err) {
	if (ns == 0) {
	    /* set all defaults */
	    state_vars_init(newstate);
	} else {
	    /* copy existing state */
	    state_vars_copy(newstate);
	}
	state_stack = sstack;
	state = state_stack[ns] = newstate;
	n_states++;
    }

#if PDEBUG
    if (!err) {
	fprintf(stderr, " state is now state_stack[%d]\n", ns);
    }
#endif

    return err;
}

static void free_state (set_vars *sv)
{
    if (sv->initvals != NULL) {
	gretl_matrix_free(sv->initvals);
    }

    free(sv);
}

/* Called when a user-defined function exits: restores the program
   state that was in force when the function started executing.
*/

int pop_program_state (void)
{
    set_vars **sstack;
    int ns = n_states;
    int err = 0;

#if PDEBUG
    fprintf(stderr, "pop_program_state called: ns=%d\n", ns);
#endif

    if (ns < 2) {
	err = 1;
    } else {
	free_state(state_stack[ns - 1]);
	state_stack[ns - 1] = NULL;
	sstack = realloc(state_stack, (ns - 1) * sizeof *sstack);
	if (sstack == NULL) {
	    err = 1;
	}	
    }

    if (!err) {
	state_stack = sstack;
	state = state_stack[ns - 2];
	n_states--;
    }

#if PDEBUG
    fprintf(stderr, " state is now state_stack[%d]\n", ns - 2);
#endif

    return err;
}

/* initialization of all user-settable settings */

int libset_init (void)
{
    static int done;
    int err = 0;

#if PDEBUG
    fprintf(stderr, "libset_init called, done=%d\n", done);
#endif

    if (!done) {
	err = push_program_state();
	done = 1;
    }

    return err;
}

void libset_cleanup (void)
{
    int i;

#if PDEBUG
    fprintf(stderr, "libset_cleanup called\n");
#endif

    for (i=0; i<n_states; i++) {
	free_state(state_stack[i]);
    }

    free(state_stack);
    state_stack = NULL;
    n_states = 0;
}

/* switches for looping, batch mode, and pausing between screens of
   output: these depend on the state of the program calling libgretl,
   but are not user-settable
*/

static int gretl_text_pause;
static int loop_on;
static int batch_mode;
static int gui_mode;

int gretl_get_text_pause (void)
{
    return gretl_text_pause;
}

void set_loop_on (void)
{
    loop_on = 1;
    gretl_text_pause = 0;
}

void set_loop_off (void)
{
    loop_on = 0;
    if (!batch_mode && !gui_mode) {
	gretl_text_pause = 1;
    }
}

int gretl_looping (void)
{
    return loop_on;
}

void gretl_set_batch_mode (int b)
{
    batch_mode = b;
    if (batch_mode) {
	gretl_text_pause = 0;
    }	
}

int gretl_in_batch_mode (void)
{
    return batch_mode;
}

void gretl_set_gui_mode (int g)
{
    gui_mode = g;
}

int gretl_in_gui_mode (void)
{
    return gui_mode;
}


