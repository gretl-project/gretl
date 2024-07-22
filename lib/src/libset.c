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

/* libset.c for gretl */

#include "libgretl.h"
#include "libset.h"
#include "usermat.h"
#include "uservar.h"
#include "matrix_extra.h"
#include "gretl_func.h"
#include "gretl_string_table.h"
#include "gretl_mt.h"
#include "gretl_foreign.h"

#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

#include <stddef.h>
#include <unistd.h>
#include <errno.h>

#include <glib.h>

#define PDEBUG 0
#define SVDEBUG 0

enum {
    AUTO_LAG_STOCK_WATSON,
    AUTO_LAG_WOOLDRIDGE,
    AUTO_LAG_NEWEYWEST
};

enum {
    CAT_BEHAVE = 1,
    CAT_NUMERIC,
    CAT_RNG,
    CAT_ROBUST,
    CAT_TS,
    CAT_SPECIAL
};

typedef enum {
    SV_ALL,
    SV_INT,
    SV_DOUBLE
} SVType;

/* for values that really want a non-negative integer */
#define UNSET_INT -9
#define is_unset(i) (i == UNSET_INT)

typedef struct set_state_ set_state;

struct set_state_ {
    int flags;
    /* small integer values */
    gint8 optim;                /* code for preferred optimizer */
    gint8 vecm_norm;            /* VECM beta normalization */
    gint8 garch_vcv;            /* GARCH vcv variant */
    gint8 garch_alt_vcv;        /* GARCH vcv variant, robust estimation */
    gint8 arma_vcv;             /* ARMA vcv variant */
    gint8 wildboot_d;           /* distribution for wild bootstrap */
    gint8 fdjac_qual;           /* quality of "fdjac" function */
    gint8 use_qr;               /* off, on or pivot */
    gint8 max_verbose;          /* optimizer verbosity level */
    gint8 hc_version;           /* HCCME version */
    gint8 panel_robust;         /* panel robust vcv estimator */
    gint8 hac_kernel;           /* HAC kernel type */
    gint8 auto_hac_lag;         /* HAC automatic lag-length formula */
    gint8 user_hac_lag;         /* fixed user-set HAC lag length */
    gint8 lbfgs_mem;            /* memory span for L-BFGS-B */
    gint8 quantile_type;        /* Formula for computing quantiles */
    /* potentially larger integers */
    int horizon;                /* for VAR impulse responses */
    int bootrep;                /* bootstrap replications */
    int loop_maxiter;           /* max no. of iterations in non-for loops */
    int bfgs_maxiter;           /* max iterations, BFGS */
    int bfgs_verbskip;          /* BFGS: show one in n iterations  */
    int boot_iters;             /* max iterations, IRF bootstrap */
    int bhhh_maxiter;           /* max iterations, BHHH */
    int rq_maxiter;             /* max iterations for quantreg, simplex */
    int gmm_maxiter;            /* max iterations for iterated GMM */
    /* floating-point values */
    double conv_huge;           /* conventional value for $huge */
    double nls_toler;           /* NLS convergence criterion */
    double bfgs_toler;          /* convergence tolerance, BFGS */
    double bfgs_maxgrad;        /* max acceptable gradient norm, BFGS */
    double bhhh_toler;          /* convergence tolerance, BHHH */
    double qs_bandwidth;        /* bandwidth for QS HAC kernel */
    double nadarwat_trim;       /* multiple of h to use in nadarwat() for trimming */
    /* strings */
    char csv_write_na[8];       /* representation of NA in CSV output */
    char csv_read_na[8];        /* representation of NA in CSV input */
    /* matrices */
    gretl_matrix *initvals;     /* parameter initializer */
    gretl_matrix *initcurv;     /* initial curvature matrix for BFGS */
    gretl_matrix *matmask;      /* mask for series -> matrix conversion */
};

typedef struct global_vars_ global_vars;

static struct global_vars_ {
    gint8 gretl_debug;
    gint8 gretl_assert;
    gint8 datacols;
    gint8 plot_collect;
    gint8 R_functions;
    gint8 R_lib;
    gint8 loglevel;
    gint8 logstamp;
    gint8 csv_digits;
    gint8 hac_missvals;
    int gmp_bits;
} globals = {0, 0, 5, 0, 0, 1, 2, 0, UNSET_INT, HAC_ES, 256};

/* globals for internal use */
static int seed_is_set;
static int comments_on;
static char data_delim = ',';
static char export_decpoint = '.';

typedef struct setvar_ setvar;

struct setvar_ {
    SetKey key;       /* internal integer key */
    const char *name; /* userspace name */
    gint8 category;   /* for printing purposes */
    size_t offset;    /* byte offset into state or globals struct,
			 where applicable */
};

setvar setvars[] = {
    /* booleans (bitflags) */
    { USE_CWD,      "use_cwd",   CAT_BEHAVE },
    { ECHO_ON,      "echo",      CAT_BEHAVE },
    { MSGS_ON,      "messages",  CAT_BEHAVE },
    { FORCE_DECPOINT, "force_decpoint", CAT_BEHAVE },
    { USE_SVD,      "svd",       CAT_NUMERIC },
    { PREWHITEN,    "hac_prewhiten", CAT_ROBUST },
    { FORCE_HC,     "force_hc",      CAT_ROBUST },
    { USE_LBFGS,    "lbfgs",        CAT_NUMERIC },
    { SHELL_OK,     "shell_ok",     CAT_SPECIAL },
    { WARNINGS,     "warnings",     CAT_BEHAVE },
    { SKIP_MISSING, "skip_missing", CAT_BEHAVE },
    { BFGS_RSTEP,   "bfgs_richardson", CAT_NUMERIC },
    { ROBUST_Z,     "robust_z", CAT_ROBUST },
    { MWRITE_G,     "mwrite_g", CAT_BEHAVE },
    { MPI_USE_SMT,  "mpi_use_smt", CAT_BEHAVE },
    { STATE_FLAG_MAX, NULL },
    /* small integers */
    { GRETL_OPTIM,  "optimizer", CAT_NUMERIC, offsetof(set_state,optim) },
    { VECM_NORM,    "vecm_norm", CAT_TS,      offsetof(set_state,vecm_norm) },
    { GARCH_VCV,    "garch_vcv", CAT_ROBUST,  offsetof(set_state,garch_vcv) },
    { GARCH_ALT_VCV, "garch_alt_vcv", CAT_ROBUST, offsetof(set_state,garch_alt_vcv) },
    { ARMA_VCV,      "arma_vcv", CAT_ROBUST, offsetof(set_state,arma_vcv) },
    { WILDBOOT_DIST, "wildboot", CAT_BEHAVE, offsetof(set_state,wildboot_d) },
    { FDJAC_QUAL,    "fdjac_quality", CAT_NUMERIC, offsetof(set_state,fdjac_qual) },
    { USE_QR,        "force_qr", CAT_NUMERIC, offsetof(set_state,use_qr) },
    { MAX_VERBOSE,   "max_verbose", CAT_BEHAVE, offsetof(set_state,max_verbose) },
    { HC_VERSION,    "hc_version",  CAT_ROBUST, offsetof(set_state,hc_version) },
    { PANEL_ROBUST,  "panel_robust", CAT_ROBUST, offsetof(set_state,panel_robust) },
    { HAC_KERNEL,    "hac_kernel",  CAT_ROBUST, offsetof(set_state,hac_kernel) },
    { HAC_LAG,       "hac_lag",     CAT_ROBUST },
    { USER_HAC_LAG,  NULL },
    { LBFGS_MEM,     "lbfgs_mem",     CAT_NUMERIC, offsetof(set_state,lbfgs_mem) },
    { QUANTILE_TYPE, "quantile_type", CAT_BEHAVE, offsetof(set_state,quantile_type) },
    { STATE_SMALL_INT_MAX, NULL },
    /* larger integers */
    { HORIZON,       "horizon", CAT_TS,     offsetof(set_state,horizon) },
    { BOOTREP,       "bootrep", CAT_BEHAVE, offsetof(set_state,bootrep) },
    { LOOP_MAXITER,  "loop_maxiter",  CAT_BEHAVE,  offsetof(set_state,loop_maxiter) },
    { BFGS_MAXITER,  "bfgs_maxiter",  CAT_NUMERIC, offsetof(set_state,bfgs_maxiter) },
    { BFGS_VERBSKIP, "bfgs_verbskip", CAT_BEHAVE,  offsetof(set_state,bfgs_verbskip) },
    { BOOT_ITERS,    "boot_iters",    CAT_TS,      offsetof(set_state,boot_iters) },
    { BHHH_MAXITER,  "bhhh_maxiter",  CAT_NUMERIC, offsetof(set_state,bhhh_maxiter) },
    { RQ_MAXITER,    "rq_maxiter",    CAT_NUMERIC, offsetof(set_state,rq_maxiter) },
    { GMM_MAXITER,   "gmm_maxiter",   CAT_NUMERIC, offsetof(set_state,gmm_maxiter) },
    { STATE_INT_MAX, NULL },
    /* doubles */
    { CONV_HUGE,     "huge",         CAT_BEHAVE,  offsetof(set_state,conv_huge) },
    { NLS_TOLER,     "nls_toler",    CAT_NUMERIC, offsetof(set_state,nls_toler) },
    { BFGS_TOLER,    "bfgs_toler",   CAT_NUMERIC, offsetof(set_state,bfgs_toler) },
    { BFGS_MAXGRAD,  "bfgs_maxgrad", CAT_NUMERIC, offsetof(set_state,bfgs_maxgrad) },
    { BHHH_TOLER,    "bhhh_toler",   CAT_NUMERIC, offsetof(set_state,bhhh_toler) },
    { QS_BANDWIDTH,  "qs_bandwidth", CAT_ROBUST,  offsetof(set_state,qs_bandwidth) },
    { NADARWAT_TRIM, "nadarwat_trim", CAT_NUMERIC, offsetof(set_state,nadarwat_trim) },
    { STATE_FLOAT_MAX, NULL },
    /* strings */
    { CSV_WRITE_NA,  "csv_write_na", CAT_SPECIAL },
    { CSV_READ_NA,   "csv_read_na",  CAT_SPECIAL },
    /* matrices */
    { INITVALS,      "initvals",    CAT_NUMERIC },
    { INITCURV,      "initcurv",    CAT_NUMERIC },
    { MATMASK,       "matrix_mask", CAT_BEHAVE },
    { STATE_VARS_MAX, NULL },
    /* global ints */
    { GRETL_DEBUG,   "debug",     CAT_BEHAVE, offsetof(global_vars,gretl_debug) },
    { GRETL_ASSERT,  "assert",    CAT_BEHAVE, offsetof(global_vars,gretl_assert) },
    { DATACOLS,      "datacols",  CAT_BEHAVE, offsetof(global_vars,datacols) },
    { PLOT_COLLECT,  "plot_collection", CAT_BEHAVE, offsetof(global_vars,plot_collect) },
    { R_FUNCTIONS,   "R_functions", CAT_BEHAVE, offsetof(global_vars,R_functions) },
    { R_LIB,         "R_lib",       CAT_BEHAVE, offsetof(global_vars,R_lib) },
    { LOGLEVEL,      "loglevel",    CAT_BEHAVE, offsetof(global_vars,loglevel) },
    { LOGSTAMP,      "logstamp",    CAT_BEHAVE, offsetof(global_vars,logstamp) },
    { CSV_DIGITS,    "csv_digits",  CAT_BEHAVE, offsetof(global_vars,csv_digits) },
    { HAC_MISSVALS,  "hac_missvals", CAT_BEHAVE, offsetof(global_vars,hac_missvals) },
    { NS_SMALL_INT_MAX, NULL },
    { GMP_BITS,      "gmp_bits",    CAT_BEHAVE, offsetof(global_vars,gmp_bits) },
    { NS_MAX, NULL },
    /* delegated ints */
    { BLAS_MNK_MIN,  "blas_mnk_min", CAT_BEHAVE },
    { OMP_MNK_MIN,   "omp_mnk_min",  CAT_BEHAVE },
    { OMP_N_THREADS, "omp_num_threads", CAT_SPECIAL },
    { SIMD_K_MAX,    "simd_k_max",  CAT_BEHAVE },
    { SIMD_MN_MIN,   "simd_mn_min", CAT_BEHAVE },
    { USE_DCMT,      "use_dcmt",    CAT_RNG },
    /* specials */
    { SEED,          "seed",      CAT_RNG },
    { CSV_DELIM,     "csv_delim", CAT_SPECIAL },
    { STOPWATCH,     "stopwatch", CAT_SPECIAL },
    { VERBOSE,       "verbose",   CAT_SPECIAL },
    { SV_WORKDIR,    "workdir",   CAT_SPECIAL },
    { SV_LOGFILE,    "logfile",   CAT_SPECIAL },
    { GRAPH_THEME,   "graph_theme", CAT_SPECIAL },
    { DISP_DIGITS,   "display_digits", CAT_SPECIAL },
    { TEX_PLOT_OPTS, "tex_plot_opts", CAT_SPECIAL }
};

#define libset_boolvar(k) (k < STATE_FLAG_MAX || k==R_FUNCTIONS || \
			   k==R_LIB || k==LOGSTAMP || k==USE_DCMT)
#define libset_double(k) (k > STATE_INT_MAX && k < STATE_FLOAT_MAX)
#define libset_int(k) ((k > STATE_FLAG_MAX && k < STATE_INT_MAX) || \
		       (k > STATE_VARS_MAX && k < NS_INT_MAX))

#define libset_small_int(k) (k < STATE_SMALL_INT_MAX || \
			     (k > STATE_VARS_MAX && k < NS_SMALL_INT_MAX))

#define coded_intvar(k) (k == GARCH_VCV || \
			 k == GARCH_ALT_VCV || \
			 k == ARMA_VCV || \
			 k == HAC_LAG || \
			 k == HAC_KERNEL || \
                         k == HC_VERSION || \
			 k == PANEL_ROBUST || \
			 k == USE_QR || \
			 k == VECM_NORM || \
			 k == GRETL_OPTIM || \
			 k == MAX_VERBOSE || \
			 k == WILDBOOT_DIST || \
			 k == QUANTILE_TYPE || \
			 k == GRETL_ASSERT || \
			 k == PLOT_COLLECT || \
			 k == LOGLEVEL || \
			 k == HAC_MISSVALS)

/* the current set of state variables */
set_state *state;

static const char *hac_lag_string (void);
static int real_libset_read_script (const char *fname,
				    PRN *prn);
static int libset_get_scalar (SetKey key, const char *arg,
			      int *pi, double *px);
static int libset_int_min (SetKey key);

/* In case we have a SetKey and want to print the
   associated userspace keyword. */

static const char *setkey_get_name (SetKey key)
{
    int i;

    for (i=0; i<G_N_ELEMENTS(setvars); i++) {
	if (setvars[i].key == key) {
	    return setvars[i].name;
	}
    }
    return NULL;
}

/* Set up a hash table to map from userspace keywords
   to setvar structs.
*/

static GHashTable *libset_hash_init (void)
{
    GHashTable *ht = g_hash_table_new(g_str_hash, g_str_equal);
    int i;

    for (i=0; i<G_N_ELEMENTS(setvars); i++) {
	if (setvars[i].name != NULL) {
	    g_hash_table_insert(ht, (gpointer) setvars[i].name, &setvars[i]);
	}
    }

    return ht;
}

static GHashTable *svht;

static void libset_hash_cleanup (void)
{
    if (svht != NULL) {
        g_hash_table_destroy(svht);
        svht = NULL;
    }
}

static setvar *get_setvar_by_name (const char *name)
{
    setvar *ret = NULL;

    if (svht == NULL) {
	svht = libset_hash_init();
    }

    ret = g_hash_table_lookup(svht, name);

    if (ret == NULL) {
	/* backward compatibility */
	if (!strcmp(name, "csv_na")) {
	    ret = g_hash_table_lookup(svht, "csv_write_na");
	} else if (!strcmp(name, "mp_mnk_min")) {
	    ret = g_hash_table_lookup(svht, "omp_mnk_min");
	}
    }

    return ret;
}

static void *setvar_get_target (setvar *sv)
{
    void *p;

    if (sv->offset == 0 || sv->key > GMP_BITS) {
	/* FIXME criterion? (0 offset is legit for "debug") */
	if (sv->key != GRETL_DEBUG) {
	    return NULL;
	}
    }

    p = (sv->key < STATE_VARS_MAX)? (void *) state : (void *) &globals;
#if SVDEBUG
    fprintf(stderr, "setvar_get_target: '%s': %s=%p, offset=%lu, ret=%p\n",
	    sv->name, sv->key < STATE_VARS_MAX ? "state" : "globals",
	    (void *) p, sv->offset, p + sv->offset);
#endif
    return p + sv->offset;
}

#define INTS_OFFSET (1 + log2(STATE_FLAG_MAX))

static void *setkey_get_target (SetKey key, SVType t)
{
    int i = INTS_OFFSET + key - GRETL_OPTIM;
    setvar *sv = &setvars[i];

    if (sv->key != key) {
	fprintf(stderr, "*** internal error, looking for %s, found %s ***\n",
		setkey_get_name(key), sv->name);
	return NULL;
    } else if ((t == SV_INT && !libset_int(key)) ||
	       (t == SV_DOUBLE && !libset_double(key))) {
	fprintf(stderr, "*** type mismatch in setkey_get_target for %s ***\n",
		sv->name);
	return NULL;
    } else {
	return setvar_get_target(sv);
    }
}

/* value strings for integer-coded variables */

static const char *gvc_strs[] = {"unset", "hessian", "im", "op", "qml", "bw", NULL};
static const char *gvr_strs[] = {"qml", "bw", NULL};
static const char *avc_strs[] = {"hessian", "op", NULL};
static const char *hkn_strs[] = {"bartlett", "parzen", "qs", NULL};
static const char *hcv_strs[] = {"0", "1", "2", "3", "3a", NULL};
static const char *vnm_strs[] = {"phillips", "diag", "first", "none", NULL};
static const char *opt_strs[] = {"auto", "BFGS", "newton", NULL};
static const char *mxv_strs[] = {"off", "on", "full", NULL};
static const char *wbt_strs[] = {"rademacher", "mammen", NULL};
static const char *qnt_strs[] = {"Q6", "Q7", "Q8", NULL};
static const char *ast_strs[] = {"off", "warn", "stop", "fatal", NULL};
static const char *plc_strs[] = {"off", "auto", "on", NULL};
static const char *csv_strs[] = {"comma", "space", "tab", "semicolon", "pipe", NULL};
static const char *ahl_strs[] = {"nw1", "nw2", "nw3", NULL};
static const char *llv_strs[] = {"debug", "info", "warn", "error", "critical", NULL};
static const char *qrp_strs[] = {"off", "on", "pivot", NULL};
static const char *hmv_strs[] = {"off", "es", "am", NULL};
static const char *pnr_strs[] = {"arellano", "pcse", "scc", NULL};

struct codevar_info {
    SetKey key;
    const char **strvals;
};

/* look-up table for sets of value strings */

struct codevar_info coded[] = {
    { GARCH_VCV,     gvc_strs },
    { GARCH_ALT_VCV, gvr_strs },
    { ARMA_VCV,      avc_strs },
    { HAC_KERNEL,    hkn_strs },
    { HC_VERSION,    hcv_strs },
    { PANEL_ROBUST,  pnr_strs },
    { VECM_NORM,     vnm_strs },
    { GRETL_OPTIM,   opt_strs },
    { MAX_VERBOSE,   mxv_strs },
    { WILDBOOT_DIST, wbt_strs },
    { CSV_DELIM,     csv_strs },
    { QUANTILE_TYPE, qnt_strs },
    { GRETL_ASSERT,  ast_strs },
    { PLOT_COLLECT,  plc_strs },
    { HAC_LAG,       ahl_strs },
    { LOGLEVEL,      llv_strs },
    { USE_QR,        qrp_strs },
    { HAC_MISSVALS,  hmv_strs }
};

static const char **libset_option_strings (SetKey key)
{
    int i;

    for (i=0; i<G_N_ELEMENTS(coded); i++) {
	if (coded[i].key == key) {
	    return coded[i].strvals;
	}
    }
    return NULL;
}

static void coded_var_show_opts (SetKey key, PRN *prn)
{
    const char **S = libset_option_strings(key);

    if (S != NULL) {
	pputs(prn, "valid settings:");
	while (*S != NULL) {
	    pprintf(prn, " %s", *S);
	    S++;
	}
	pputc(prn, '\n');
    }
}

static const char *garch_alt_vcv_string (void)
{
    if (state->garch_alt_vcv == ML_QML) {
	return gvr_strs[0];
    } else if (state->garch_alt_vcv == ML_BW) {
	return gvr_strs[1];
    } else {
	return "unset";
    }
}

static const char *arma_vcv_string (void)
{
    if (state->arma_vcv == ML_HESSIAN) {
	return avc_strs[0];
    } else if (state->arma_vcv == ML_OP) {
	return avc_strs[1];
    } else {
	return "unset";
    }
}

static const char *libset_option_string (SetKey key)
{
    if (key == HAC_LAG) {
	return hac_lag_string();          /* special */
    } else if (key == GARCH_ALT_VCV) {
	return garch_alt_vcv_string();    /* special */
    } else if (key == ARMA_VCV) {
	return arma_vcv_string();         /* special */
    } else {
	const char **strs = libset_option_strings(key);
	void *valp = setkey_get_target(key, SV_INT);

	return strs[*(gint8 *) valp];
    }
}

static void print_state_matrix (SetKey key, PRN *prn, gretlopt opt)
{
    gretl_matrix *m;
    char *name;

    if (key == INITVALS) {
	name = "initvals";
	m = state->initvals;
    } else {
	name = "initcurv";
	m = state->initcurv;
    }

    if (opt & OPT_D) {
	if (m == NULL) {
	    pprintf(prn, " %s = auto\n", name);
	} else {
	    gretl_matrix_print_to_prn(m, name, prn);
	}
    }
}

/* check_for_state() returns non-zero if the program options
   state is not readable */

static int check_for_state (void)
{
    if (state == NULL) {
	return libset_init();
    } else {
#if PDEBUG > 1
	fprintf(stderr, "check_for_state: state = %p\n", (void *) state);
#endif
	return 0;
    }
}

static int flag_to_bool (set_state *sv, SetKey key)
{
    if (sv == NULL) {
	return 0;
    } else {
	return (sv->flags & key)? 1 : 0;
    }
}

static void state_vars_copy (set_state *sv)
{
#if PDEBUG
    fprintf(stderr, "state_vars_copy() called\n");
#endif
    /* copy everything */
    *sv = *state;
    /* but set matrix pointers to NULL */
    sv->initvals = NULL;
    sv->initcurv = NULL;
    sv->matmask = NULL;
}

static set_state default_state = {
    ECHO_ON | MSGS_ON | WARNINGS | SKIP_MISSING, /* .flags */
    OPTIM_AUTO,     /* .optim */
    NORM_PHILLIPS,  /* .vecm_norm */
    ML_UNSET,       /* .garch_vcv */
    ML_UNSET,       /* .garch_alt_vcv */
    ML_HESSIAN,     /* .arma_vcv */
    0,              /* .wildboot_dist */
    0,              /* .fdjac_qual */
    0,              /* .use_qr */
    0,              /* .max_verbose */
    1,              /* .hc_version */
    0,              /* .panel_robust */
    KERNEL_BARTLETT,       /* .hac_kernel */
    AUTO_LAG_STOCK_WATSON, /* .auto_hac_lag */
    UNSET_INT,             /* .user_hac_lag */
    8,             /* .lbfgs_mem */
    0,             /* .quantile_type */
    UNSET_INT,     /* .horizon */
    1000,          /* .bootrep */
    100000,        /* .loop_maxiter */
    UNSET_INT,     /* .bfgs_maxiter */
    1,       /* .bfgs_verbskip */
    1999,    /* .boot_iters */
    500,     /* .bhhh_maxiter */
    1000,    /* .rq_maxiter */
    250,     /* .gmm_maxiter */
    1.0e100, /* .conv_huge */
    NADBL,   /* .nls_toler */
    NADBL,   /* .bfgs_toler */
    5.0,     /* .bfgs_maxgrad */
    1.0e-6,  /* .bhhh_toler */
    2.0,     /* .qs_bandwidth */
    4.0,     /* .nadarwat_trim */
    "NA",      /* .csv_write_na */
    "default", /* .csv_read_na */
    NULL,  /* .initvals */
    NULL,  /* .initcurv */
    NULL   /* .matmask */
};

static void state_vars_init (set_state *sv)
{
    *sv = default_state;
}

int get_bkbp_k (const DATASET *dset)
{
    if (dset->pd == 1) {
	return 3;
    } else if (dset->pd == 4) {
	return 12;
    } else if (dset->pd == 12) {
	return 36;
    } else {
	return 3;
    }
}

void get_bkbp_periods (const DATASET *dset, int *l, int *u)
{
    *l = (dset->pd == 4)? 6 :
	(dset->pd == 12)? 18 : 2;

    *u = (dset->pd == 4)? 32 :
	(dset->pd == 12)? 96 : 8;
}

void set_gretl_echo (int e)
{
    if (check_for_state()) return;

    if (e) {
	state->flags |= ECHO_ON;
    } else {
	state->flags &= ~ECHO_ON;
    }
}

int gretl_echo_on (void)
{
    if (check_for_state()) {
	return 1;
    } else {
	return flag_to_bool(state, ECHO_ON);
    }
}

int gretl_comments_on (void)
{
    if (gretl_function_depth() > 0) {
	return 0;
    } else {
	return comments_on;
    }
}

void set_gretl_messages (int e)
{
    if (check_for_state()) return;

    if (e) {
	state->flags |= MSGS_ON;
    } else {
	state->flags &= ~MSGS_ON;
    }
}

int gretl_messages_on (void)
{
    if (check_for_state()) {
	return 1;
    } else {
	return flag_to_bool(state, MSGS_ON);
    }
}

int gretl_warnings_on (void)
{
    if (check_for_state()) return 1;
    return flag_to_bool(state, WARNINGS);
}

int gretl_debugging_on (void)
{
    return globals.gretl_debug;
}

#define DEFAULT_MP_BITS 256
#define mp_bits_ok(b) (b >= 256 && b <= 8192)

/* Called from the mp_ols plugin, also gretl_utils.c and
   the GUI model specification dialog.
*/

int get_mp_bits (void)
{
    if (globals.gmp_bits > DEFAULT_MP_BITS) {
	return globals.gmp_bits;
    } else {
	char *s = getenv("GRETL_MP_BITS");
	int b;

	if (s != NULL) {
	    b = atoi(s);
	    if (mp_bits_ok(b)) {
		return b;
	    }
	}
    }

    return DEFAULT_MP_BITS;
}

/* start accessors for libset matrices */

gretl_matrix *get_initvals (void)
{
    gretl_matrix *iv;

    /* note: we nullify initvals after first use */
    check_for_state();
    iv = state->initvals;
    state->initvals = NULL;
    return iv;
}

int n_initvals (void)
{
    check_for_state();
    if (state->initvals != NULL) {
	return gretl_vector_get_length(state->initvals);
    } else {
	return 0;
    }
}

gretl_matrix *get_initcurv (void)
{
    gretl_matrix *ic;

    /* note: like initvals, we nullify initcurv after first use */
    check_for_state();
    ic = state->initcurv;
    state->initcurv = NULL;
    return ic;
}

int n_initcurv (void)
{
    check_for_state();
    if (state->initcurv != NULL) {
	return gretl_vector_get_length(state->initcurv);
    } else {
	return 0;
    }
}

const gretl_matrix *get_matrix_mask (void)
{
    check_for_state();
    return state->matmask;
}

int get_matrix_mask_nobs (void)
{
    int n = 0;

    check_for_state();

    if (state->matmask != NULL) {
	int i;

	for (i=0; i<state->matmask->rows; i++) {
	    if (state->matmask->val[i] != 0.0) {
		n++;
	    }
	}
    }

    return n;
}

/* end accessors for libset matrices */

int get_hac_lag (int T)
{
    int h = 0;

    check_for_state();

    /* Variants of Newey-West */

    if (state->user_hac_lag >= 0 && state->user_hac_lag < T - 2) {
	/* FIXME upper limit? */
	h = state->user_hac_lag;
    } else if (state->auto_hac_lag == AUTO_LAG_WOOLDRIDGE) {
	h = (int) floor(4 * pow(T/100.0, 2.0/9));
    } else {
	/* Stock-Watson default */
	h = (int) floor(0.75 * pow(T, 1.0/3));
    }

    return h;
}

/* prewhitening implies nw3, but not vice versa */

int data_based_hac_bandwidth (void)
{
    if (is_unset(state->user_hac_lag)) {
	if (state->auto_hac_lag == AUTO_LAG_NEWEYWEST ||
	    (state->flags & PREWHITEN)) {
	    return 1;
	}
    }

    return 0;
}

static const char *hac_lag_string (void)
{
    check_for_state();

    if (state->user_hac_lag >= 0 && state->user_hac_lag < 127) {
	static char lagstr[6];

	sprintf(lagstr, "%d", state->user_hac_lag);
	return lagstr;
    } else {
	return ahl_strs[state->auto_hac_lag];
    }
}

/* set max lag for HAC estimation */

static int parse_hac_lag_variant (const char *s)
{
    int i, err = 0;

    for (i=0; ahl_strs[i] != NULL; i++) {
	if (!strcmp(s, ahl_strs[i])) {
	    state->auto_hac_lag = i;
	    state->user_hac_lag = UNSET_INT;
	    return 0;
	}
    }

    err = libset_get_scalar(HAC_LAG, s, &i, NULL);
    if (!err) {
	state->user_hac_lag = i;
    }

    return err;
}

enum {
    NUMERIC_OK,
    NUMERIC_BAD,
    NON_NUMERIC
};

/* Test @s for being a string representation of a numeric value
   in the C locale, using strtod() and/or strtol(). This can be
   used to retrieve a floating-point value (if @px is non-NULL)
   or an integer value (@pi non-NULL); exactly one of these
   pointers should be non-NULL.

   A return value of NUMERIC_OK means that @s is indeed numeric
   and the converted value is within range for a double (if @px
   is non-NULL) or a 32-bit integer (if @pi is non-NULL).

   A return of NUMERIC_BAD means that @s is numeric but out of
   range for the target type.

   A return of NON_NUMERIC means that @s is not a numeric string;
   one may then proceed to test whether it's the name of a scalar
   variable.
*/

static int
libset_numeric_test (const char *s, int *pi, double *px)
{
    int ret = NUMERIC_OK;
    char *test;

    if (!strcmp(s, "inf") || !strcmp(s, "nan")) {
	return NUMERIC_BAD;
    } else if (isalpha(*s)) {
	return NON_NUMERIC;
    }

    errno = 0;
    gretl_push_c_numeric_locale();

    if (px != NULL) {
	/* looking for a floating-point value */
	*px = strtod(s, &test);
	if (*test != '\0') {
	    ret = NON_NUMERIC;
	} else if (errno == ERANGE) {
	    gretl_errmsg_set_from_errno(s, errno);
	    ret = NUMERIC_BAD;
	}
    } else {
	/* looking for an integer value */
	long li = strtol(s, &test, 10);

	if (*test != '\0') {
	    /* try for a floating-point value that's also a valid int? */
	    char *testx;
	    double x;

	    errno = 0;
	    x = strtod(s, &testx);

	    if (*testx != '\0') {
		ret = NON_NUMERIC;
	    } else if (errno == ERANGE) {
		ret = NUMERIC_BAD;
	    } else {
		/* numeric, but does it work as an int? */
		if (x == floor(x) && fabs(x) <= INT_MAX) {
		    *pi = (int) x;
		} else {
		    ret = NUMERIC_BAD;
		}
	    }
	} else if (errno == ERANGE) {
	    gretl_errmsg_set_from_errno(s, errno);
	    ret = NUMERIC_BAD;
	} else if (labs(li) > INT_MAX) {
	    /* OK as a long but too big for 32-bit int */
	    ret = NUMERIC_BAD;
	} else {
	    *pi = (int) li;
	}
    }

    gretl_pop_c_numeric_locale();

    return ret;
}

static int negval_invalid (SetKey key)
{
    int ret = 1; /* presume invalid */

    if (key > 0) {
	if (key == BLAS_MNK_MIN || key == OMP_MNK_MIN ||
	    key == SIMD_K_MAX || key == SIMD_MN_MIN) {
	    /* these can all be set to -1 */
	    ret = 0;
	}
    }

    return ret;
}

static int libset_get_scalar (SetKey key, const char *arg,
			      int *pi, double *px)
{
    double x = NADBL;
    int nstatus, err = 0;

    if (arg == NULL || *arg == '\0') {
	return E_ARGS;
    }

    nstatus = libset_numeric_test(arg, pi, px);

    if (nstatus == NUMERIC_BAD) {
	return E_INVARG; /* handled */
    } else if (nstatus == NUMERIC_OK) {
	if (pi != NULL && negval_invalid(key) && *pi < 0) {
	    err = E_INVARG;
	} else if (px != NULL && *px < 0.0) {
	    err = E_INVARG;
	}
	return err; /* handled */
    }

    /* handle the non-numeric case */
    x = get_scalar_value_by_name(arg, &err);

    if (!err) {
	if (negval_invalid(key) && x < 0.0) {
	    err = E_INVARG;
	} else if (px != NULL) {
	    *px = x;
	} else if (pi != NULL) {
	    if (na(x) || fabs(x) > (double) INT_MAX) {
		err = E_INVARG;
	    } else {
		*pi = (int) x;
	    }
	}
    }

    return err;
}

/* Called only when setting the PRNG seed */

static int libset_get_unsigned (const char *arg, unsigned int *pu,
				int *automatic)
{
    unsigned long lu = 0;
    char *test = NULL;
    double x = NADBL;
    int err = 0;

    if (!strcmp(arg, "auto")) {
	*automatic = 1;
	*pu = time(NULL);
	return 0;
    }

    errno = 0;
    lu = strtoul(arg, &test, 10);

    if (*test == '\0' && errno == 0) {
	if (lu <= UINT_MAX) {
	    *pu = (unsigned) lu;
	    return 0;
	} else {
	    return E_DATA;
	}
    }

    x = get_scalar_value_by_name(arg, &err);
    if (err) {
	return err;
    }

    if (x < 0.0 || na(x) || x > (double) UINT_MAX) {
	err = E_DATA;
    } else {
	*pu = (unsigned) x;
    }

    return err;
}

static int n_strvals (const char **s)
{
    int n = 0;

    while (*s != NULL) {
	n++; s++;
    }
    return n;
}

static int parse_libset_int_code (SetKey key, const char *val)
{
    int i, err = E_DATA;

    if (key == HAC_LAG) {
	err = parse_hac_lag_variant(val);
    } else if (coded_intvar(key)) {
	const char **strs = libset_option_strings(key);
	void *valp = setkey_get_target(key, SV_INT);
	int ival = -1;

	for (i=0; strs[i] != NULL; i++) {
	    if (!g_ascii_strcasecmp(val, strs[i])) {
		ival = i;
		break;
	    }
	}
	if (ival >= 0) {
	    void *valp = setkey_get_target(key, SV_INT);

	    err = 0;
	    if (key == GARCH_ALT_VCV) {
		ival = (ival == 1)? ML_BW : ML_QML;
	    } else if (key == ARMA_VCV) {
		ival = (ival == 1)? ML_OP : ML_HESSIAN;
	    }
	    *(gint8 *) valp = ival;
	} else if (key == MAX_VERBOSE || key == LOGLEVEL) {
	    /* special: bare integers allowed? */
	    int n = n_strvals(strs);

	    for (i=0; i<n; i++) {
		if (val[0] == i + 48 && val[1] == '\0') {
		    *(gint8 *) valp = i;
		    err = 0;
		    break;
		}
	    }
	}
    }

#if 0
    fprintf(stderr, "parse_libset_int_code: %s, %s, err = %d\n",
	    setkey_get_name(key), val, err);
#endif

    if (err) {
	gretl_errmsg_sprintf(_("%s: invalid value '%s'"),
			     setkey_get_name(key), val);
    }

    return err;
}

/* start public functions called from gui/settings.c */

void set_xsect_hccme (const char *s)
{
    if (check_for_state()) return;

    if (!strncmp(s, "HC", 2)) {
	s += 2;
    }
    parse_libset_int_code(HC_VERSION, s);
}

void set_tseries_hccme (const char *s)
{
    if (check_for_state()) return;

    if (!strcmp(s, "HAC")) {
	libset_set_bool(FORCE_HC, 0);
    } else {
	if (!strncmp(s, "HC", 2)) {
	    s += 2;
	}
	if (parse_libset_int_code(HC_VERSION, s) == 0) {
	    /* non-HAC variant chosen */
	    libset_set_bool(FORCE_HC, 1);
	}
    }
}

void set_panel_hccme (const char *s)
{
    if (check_for_state()) return;

    parse_libset_int_code(PANEL_ROBUST, s);
}

void set_garch_alt_vcv (const char *s)
{
    if (check_for_state()) return;

    parse_libset_int_code(GARCH_ALT_VCV, s);
}

/* end public functions called from gui/settings.c */

static int set_init_matrix (SetKey key, const char *name,
			    PRN *prn)
{
    gretl_matrix **targ;

    targ = (key == INITVALS)? &state->initvals : &state->initcurv;

    gretl_matrix_free(*targ);
    *targ = NULL;

    if (strcmp(name, "auto")) {
	gretl_matrix *m = get_matrix_by_name(name);

	if (m == NULL) {
	    pprintf(prn, _("'%s': no such matrix"), name);
	    pputc(prn, '\n');
	    return E_DATA;
	}
	*targ = gretl_matrix_copy(m);
	if (*targ == NULL) {
	    return E_ALLOC;
	}
    }

    return 0;
}

static int set_matrix_mask (const char *name, DATASET *dset)
{
    gretl_matrix_free(state->matmask);
    state->matmask = NULL;

    if (strcmp(name, "null")) {
	int t, v = current_series_index(dset, name);

	if (v < 0) {
	    return E_UNKVAR;
	}
	state->matmask = gretl_column_vector_alloc(dset->n);
	if (state->matmask == NULL) {
	    return E_ALLOC;
	}
	for (t=0; t<dset->n; t++) {
	    state->matmask->val[t] = dset->Z[v][t];
	}
    }

    return 0;
}

void destroy_matrix_mask (void)
{
    check_for_state();
    gretl_matrix_free(state->matmask);
    state->matmask = NULL;
}

static int (*workdir_callback)(char *s);

void set_workdir_callback (int (*callback)(char *s))
{
    workdir_callback = callback;
}

static int set_workdir (const char *s)
{
    int err = 0;

    if (gretl_function_depth() > 0) {
	gretl_errmsg_set("set workdir: cannot be done inside a function");
	return 1;
    } else if (*s == '\0') {
	return E_DATA;
    } else {
	char workdir[MAXLEN];

	*workdir = '\0';
	strncat(workdir, s, MAXLEN - 1);
	if (!err && workdir_callback != NULL) {
	    err = (*workdir_callback)(workdir);
	} else if (!err) {
	    err = gretl_set_path_by_name("workdir", workdir);
	}
    }

    return err;
}

/* Apparatus for handling tex_plot_opts: we allow setting of this
   variable within a function, but only if it hasn't already been set
   at a "higher" level of execution. If setting within a function
   succeeds, we destroy the setting on exit so that it doesn't
   propagate "upwards".
*/

#define TPO_UNSET 1024
static char *tex_plot_opts;
static int tpo_depth = TPO_UNSET;

static int set_tex_plot_opts (const char *arg)
{
    int fd = gretl_function_depth();

    if (fd > tpo_depth) {
	gretl_errmsg_sprintf("'%s': cannot be overwritten here",
			     "tex_plot_opts");
	return E_INVARG;
    } else {
	free(tex_plot_opts);
	if (!strcmp(arg, "null") || !strcmp(arg, "none")) {
	    tex_plot_opts = NULL;
	} else {
	    tex_plot_opts = gretl_strdup(arg);
	    g_strstrip(tex_plot_opts);
	}
	tpo_depth = fd;
    }

    return 0;
}

/* if tex_plot_opts was set within a function, scrub it on exit */

static void maybe_destroy_tpo (void)
{
    if (gretl_function_depth() == tpo_depth) {
	free(tex_plot_opts);
	tex_plot_opts = NULL;
	tpo_depth = TPO_UNSET;
    }
}

/* accessor, called from graphing.c */

const char *get_tex_plot_opts (void)
{
    return tex_plot_opts;
}

/* end tex_plot_opts apparatus */

static int set_logfile (const char *s)
{
    int err = 0;

    if (gretl_function_depth() > 0) {
	gretl_errmsg_set("set logfile: cannot be done inside a function");
	return E_DATA;
    } else if (*s == '\0' || !strcmp(s, "null")) {
	gretl_insert_builtin_string("logfile", "");
    } else if (!strcmp(s, "stdout") || !strcmp(s, "stderr")) {
	gretl_insert_builtin_string("logfile", s);
    } else {
        char outname[FILENAME_MAX];

        /* switch to workdir if needed */
        strcpy(outname, s);
        gretl_maybe_prepend_dir(outname);
	err = gretl_test_fopen(outname, "w");
	if (!err) {
	    gretl_insert_builtin_string("logfile", outname);
	} else {
	    gretl_errmsg_sprintf("Couldn't write to %s", outname);
	    err = E_FOPEN;
	}
    }

    return err;
}

const char *csv_delims = ", \t;|";

static char delim_from_arg (const char *s)
{
    int i;

    for (i=0; csv_strs[i] != NULL; i++) {
	if (!strcmp(s, csv_strs[i])) {
	    return csv_delims[i];
	}
    }

    return 0;
}

static const char *arg_from_delim (char c)
{
    int i;

    for (i=0; csv_delims[i] != '\0'; i++) {
	if (c == csv_delims[i]) {
	    return csv_strs[i];
	}
    }

    return "unset";
}

static void libset_print_bool (SetKey key, const char *s,
			       PRN *prn, gretlopt opt)
{
    int v = libset_get_bool(key);

    if (gretl_function_depth() > 0 &&
	(key == ECHO_ON || key == MSGS_ON || key == VERBOSE)) {
	return;
    }

    if (s == NULL) {
	s = setkey_get_name(key);
    }

    if (opt & OPT_D) {
	pprintf(prn, " %s = %d\n", s, v);
    } else {
	pprintf(prn, "set %s %s\n", s, v? "on" : "off");
    }
}

static void libset_print_int (SetKey key, const char *s,
			      PRN *prn, gretlopt opt)
{
    if (s == NULL) {
	s = setkey_get_name(key);
    }

    if (coded_intvar(key)) {
	if (opt & OPT_D) {
	    pprintf(prn, " %s = %s\n", s, libset_option_string(key));
	} else {
	    pprintf(prn, "set %s %s\n", s, libset_option_string(key));
	}
    } else {
	int k = libset_get_int(key);

	if (opt & OPT_D) {
	    if (is_unset(k)) {
		pprintf(prn, " %s = auto\n", s);
	    } else {
		pprintf(prn, " %s = %d\n", s, k);
	    }
	} else if (!is_unset(k)) {
	    pprintf(prn, "set %s %d\n", s, k);
	}
    }
}

static void libset_print_double (SetKey key, const char *s,
				 PRN *prn, gretlopt opt)
{
    double x = libset_get_double(key);

    if (s == NULL) {
	s = setkey_get_name(key);
    }

    if (opt & OPT_D) {
	if (na(x)) {
	    pprintf(prn, " %s = auto\n", s);
	} else {
	    pprintf(prn, " %s = %.15g\n", s, x);
	}
    } else if (!na(x)) {
	pprintf(prn, "set %s %.15g\n", s, x);
    }
}

static void libset_header (char *s, PRN *prn, gretlopt opt)
{
    if (opt & OPT_D) {
	pputs(prn, "\n --- ");
	pputs(prn, _(s));
	pputs(prn, " ---\n");
    } else {
	pprintf(prn, "# %s\n", s);
    }
}

static void print_vars_for_category (int category, PRN *prn,
				     gretlopt opt)
{
    setvar *v;
    int i;

    for (i=0; i<G_N_ELEMENTS(setvars); i++) {
	if (setvars[i].category != category) {
	    continue;
	}
	v = &setvars[i];
	if (libset_boolvar(v->key)) {
	    libset_print_bool(v->key, v->name, prn, opt);
	} else if (libset_int(setvars[i].key)) {
	    libset_print_int(v->key, v->name, prn, opt);
	} else if (libset_double(v->key)) {
	    libset_print_double(v->key, v->name, prn, opt);
	}
    }
}

/* print_settings: use OPT_D for "display", otherwise
   this gives script-type output */

static int print_settings (PRN *prn, gretlopt opt)
{
    const char *workdir = gretl_workdir();

    gretl_push_c_numeric_locale();

    if (opt & OPT_D) {
	pputs(prn, _("Variables that can be set using \"set\""));
	pputs(prn, " (");
	pputs(prn, _("\"help set\" for details"));
	pputs(prn, "):\n");
    }

    libset_header(N_("Program interaction and behavior"), prn, opt);

    if (opt & OPT_D) {
	pprintf(prn, " workdir = \"%s\"\n", workdir);
    } else if (0) {
	/* do we want this? */
	pprintf(prn, "set workdir \"%s\"\n", workdir);
    }

    if (opt & OPT_D) {
	pprintf(prn, " csv_delim = %s\n", arg_from_delim(data_delim));
	pprintf(prn, " csv_write_na = %s\n", get_csv_na_write_string());
	pprintf(prn, " csv_read_na = %s\n", get_csv_na_read_string());
	pprintf(prn, " display_digits = %d\n", get_gretl_digits());
	pprintf(prn, " graph_theme = %s\n", get_plotstyle());
	pprintf(prn, " tex_plot_opts = \"%s\"\n",
		tex_plot_opts == NULL ? "null" : tex_plot_opts);
    } else {
	const char *dl = arg_from_delim(data_delim);

	if (strcmp(dl, "unset")) {
	    pprintf(prn, "set csv_delim %s\n", arg_from_delim(data_delim));
	}
	pprintf(prn, "set csv_write_na %s\n", get_csv_na_write_string());
	pprintf(prn, "set csv_read_na %s\n", get_csv_na_read_string());
	pprintf(prn, "set graph_theme %s\n", get_plotstyle());
	pprintf(prn, "set tex_plot_opts \"%s\"\n",
		tex_plot_opts == NULL ? "null" : tex_plot_opts);
    }

    print_vars_for_category(CAT_BEHAVE, prn, opt);
    if (opt & OPT_D) {
	/* display only */
	libset_print_bool(SHELL_OK, NULL, prn, opt);
    }

    libset_header(N_("Numerical methods"), prn, opt);
    print_vars_for_category(CAT_NUMERIC, prn, opt);
    if ((opt & OPT_D) && gretl_function_depth() == 0) {
	/* script version of this? */
	print_state_matrix(INITVALS, prn, opt);
	print_state_matrix(INITCURV, prn, opt);
    }

    libset_header(N_("Random number generation"), prn, opt);
    if (opt & OPT_D) {
	pprintf(prn, " seed = %u\n", gretl_rand_get_seed());
    } else {
	if (seed_is_set) {
	    pprintf(prn, "set seed %u\n", gretl_rand_get_seed());
	}
    }
    if (gretl_mpi_initialized()) {
	libset_print_bool(USE_DCMT, NULL, prn, opt);
    }

    libset_header(N_("Robust estimation"), prn, opt);
    print_vars_for_category(CAT_ROBUST, prn, opt);

    libset_header(N_("Time series"), prn, opt);
    print_vars_for_category(CAT_TS, prn, opt);

    pputc(prn, '\n');

    gretl_pop_c_numeric_locale();

    return 0;
}

static int libset_query_settings (setvar *sv, PRN *prn)
{
    int err = 0;

    if (libset_boolvar(sv->key)) {
	pprintf(prn, "%s: boolean (on/off), currently %s\n",
		sv->name, libset_get_bool(sv->key)? "on" : "off");
    } else if (coded_intvar(sv->key)) {
	pprintf(prn, "%s: code, currently \"%s\"\n", sv->name,
		libset_option_string(sv->key));
	coded_var_show_opts(sv->key, prn);
    } else if (libset_int(sv->key)) {
	int k = libset_get_int(sv->key);

	if (is_unset(k)) {
	    pprintf(prn, "%s: positive integer, currently unset\n", sv->name);
	} else if (libset_int_min(sv->key) == 0) {
	    pprintf(prn, "%s: non-negative integer, currently %d\n", sv->name, k);
	} else {
	    pprintf(prn, "%s: positive integer, currently %d\n", sv->name, k);
	}
    } else if (libset_double(sv->key)) {
	double x = libset_get_double(sv->key);

	if (na(x)) {
	    pprintf(prn, "%s: positive floating-point value, "
		    "currently automatic\n", sv->name);
	} else {
	    pprintf(prn, "%s: positive floating-point value, "
		    "currently %g\n", sv->name, x);
	}
    } else if (sv->key == INITVALS || sv->key == INITCURV ||
	       sv->key == MATMASK) {
	gretl_matrix *m =
	    (sv->key == INITVALS)? state->initvals :
	    (sv->key == INITCURV)? state->initcurv : state->matmask;

	if (m != NULL) {
	    pprintf(prn, "%s: matrix, currently\n", sv->name);
	    gretl_matrix_print_to_prn(m, NULL, prn);
	} else {
	    pprintf(prn, "%s: matrix, currently null\n", sv->name);
	}
    } else if (sv->key == SEED) {
	pprintf(prn, "%s: unsigned int, currently %u (%s)\n",
		sv->name, gretl_rand_get_seed(),
		seed_is_set ? "set by user" : "automatic");
    } else if (sv->key == CSV_DELIM) {
	pprintf(prn, "%s: named character, currently \"%s\"\n", sv->name,
		arg_from_delim(data_delim));
	coded_var_show_opts(sv->key, prn);
    } else if (sv->key == SV_WORKDIR) {
	pprintf(prn, "%s: string, currently \"%s\"\n", sv->name,
		gretl_workdir());
    } else if (sv->key == CSV_WRITE_NA) {
	pprintf(prn, "%s: string, currently \"%s\"\n", sv->name,
		state->csv_write_na);
    } else if (sv->key == CSV_READ_NA) {
	pprintf(prn, "%s: string, currently \"%s\"\n", sv->name,
		state->csv_read_na);
    } else if (sv->key == DISP_DIGITS) {
	pprintf(prn, "%s: integer, currently %d\n", sv->name,
		get_gretl_digits());
    } else if (sv->key == GRAPH_THEME) {
	pprintf(prn, "%s: keyword, currently \"%s\"\n", sv->name,
		get_plotstyle());
    } else if (sv->key == TEX_PLOT_OPTS) {
	pprintf(prn, "%s: string, currently \"%s\"\n", sv->name,
		tex_plot_opts == NULL ? "null" : tex_plot_opts);
    } else if (sv->key == VERBOSE) {
	pprintf(prn, "%s: boolean (on/off), currently %s\n", sv->name,
		(libset_get_bool(ECHO_ON) || libset_get_bool(MSGS_ON)) ?
		"on" : "off");
	err = 0;
    } else {
	err = 1;
    }

    return err;
}

/* determine whether help text is available pertaining to
   putative libset variable @s.
*/

int libset_help_available (const char *s)
{
    setvar *sv = get_setvar_by_name(s);
    int err = 1;

    if (sv != NULL) {
	err = libset_query_settings(sv, NULL);
    }

    return (err == 0);
}

/* return the enumeration value corresponding to the putative
   libset variable @s, or 0 on failure
*/

SetKey get_libset_key (const char *s)
{
    setvar *sv = get_setvar_by_name(s);

    return sv != NULL ? sv->key : 0;
}

#define default_ok(k) (k == BFGS_TOLER || k == BHHH_TOLER || k == NLS_TOLER)

#define default_str(s) (!strcmp(s, "auto") || !strcmp(s, "default"))

#define boolean_on(s) (!strcmp(s, "on") || !strcmp(s, "1") || \
                       !strcmp(s, "true") || !strcmp(s, "TRUE"))

#define boolean_off(s) (!strcmp(s, "off") || !strcmp(s, "0") || \
                        !strcmp(s, "false") || !strcmp(s, "FALSE"))

static int write_or_read_settings (gretlopt opt, PRN *prn)
{
    int err = incompatible_options(opt, (OPT_T | OPT_F));

    if (!err) {
	const char *fname = get_optval_string(SET, opt);

	if (fname == NULL) {
	    err = E_DATA;
	} else if (opt == OPT_T) {
	    err = libset_write_script(fname);
	} else {
	    err = real_libset_read_script(fname, prn);
	}
    }

    return err;
}

static int check_set_bool (SetKey key, const char *name,
			   const char *arg)
{
    if (boolean_on(arg)) {
	return libset_set_bool(key, 1);
    } else if (boolean_off(arg)) {
	return libset_set_bool(key, 0);
    } else {
	gretl_errmsg_sprintf(_("%s: invalid value '%s'"), name, arg);
	return E_PARSE;
    }
}

static int legacy_set_pcse (const char *arg)
{
    int err = 0;

    if (arg == NULL || *arg == '\0') {
	err = E_PARSE;
    } else if (boolean_on(arg)) {
	libset_set_int(PANEL_ROBUST, BECK_KATZ);
    } else if (boolean_off(arg)) {
	libset_set_int(PANEL_ROBUST, ARELLANO);
    } else {
	gretl_errmsg_sprintf(_("%s: invalid value '%s'"), "pcse", arg);
	err = E_PARSE;
    }

    return err;
}

static int set_display_digits (const char *arg)
{
    if (gretl_function_depth() > 0) {
	gretl_errmsg_sprintf("'%s': cannot be set inside a function",
			     "display_digits");
	return E_INVARG;
    } else {
	return set_gretl_digits(atoi(arg));
    }
}

static int set_verbosity (const char *arg)
{
    int err = 0;

    if (!strcmp(arg, "on")) {
	set_gretl_messages(1);
	set_gretl_echo(1);
    } else if (!strcmp(arg, "off")) {
	set_gretl_messages(0);
	set_gretl_echo(0);
	comments_on = 0;
    } else if (!strcmp(arg, "comments")) {
	set_gretl_messages(0);
	set_gretl_echo(0);
	comments_on = 1;
    } else {
	err = E_INVARG;
    }

    return err;
}

int execute_set (const char *setobj, const char *setarg,
		 DATASET *dset, gretlopt opt, PRN *prn)
{
    setvar *sv = NULL;
    int k, argc, err;
    unsigned int u;

#if 0
    fprintf(stderr, "execute_set: setobj = '%s', setarg='%s'\n",
	    setobj, setarg);
#endif

    check_for_state();

    if (opt != OPT_NONE) {
	return write_or_read_settings(opt, prn);
    }

    argc = (setobj != NULL) + (setarg != NULL);
    if (argc == 0) {
	return print_settings(prn, OPT_D);
    }

    if (!strcmp(setobj, "pcse")) {
	return legacy_set_pcse(setarg);
    } else {
	sv = get_setvar_by_name(setobj);
    }

    if (sv == NULL) {
	gretl_errmsg_sprintf(_("set: unknown variable '%s'"), setobj);
	return E_DATA;
    }

    /* set error default */
    err = E_PARSE;

    if (argc == 1) {
	if (sv->key == STOPWATCH) {
	    gretl_stopwatch();
	    return 0;
	} else {
	    return libset_query_settings(sv, prn);
	}
    } else if (argc == 2) {
	/* specials first */
	if (sv->key == CSV_WRITE_NA) {
	    return set_csv_na_write_string(setarg);
	} else if (sv->key == CSV_READ_NA) {
	    return set_csv_na_read_string(setarg);
	} else if (sv->key == INITVALS || sv->key == INITCURV) {
	    return set_init_matrix(sv->key, setarg, prn);
	} else if (sv->key == MATMASK) {
	    return set_matrix_mask(setarg, dset);
	} else if (sv->key == SV_WORKDIR) {
	    return set_workdir(setarg);
	} else if (sv->key == SV_LOGFILE) {
	    return set_logfile(setarg);
	} else if (sv->key == GRAPH_THEME) {
	    return set_plotstyle(setarg);
	} else if (sv->key == DISP_DIGITS) {
	    return set_display_digits(setarg);
	} else if (sv->key == VERBOSE) {
	    return set_verbosity(setarg);
	} else if (sv->key == TEX_PLOT_OPTS) {
	    return set_tex_plot_opts(setarg);
	} else if (sv->key == OMP_MNK_MIN) {
#if defined(_OPENMP)
	    return set_omp_mnk_min(atoi(setarg));
#else
	    pprintf(prn, "Warning: openmp not supported\n");
#endif
	}

	if (libset_boolvar(sv->key)) {
	    if (sv->key == SHELL_OK) {
		pprintf(prn, "'%s': this must be set via the gretl GUI\n", setobj);
		err = E_DATA;
	    } else {
		err = check_set_bool(sv->key, setobj, setarg);
	    }
	} else if (libset_double(sv->key)) {
	    if (default_ok(sv->key) && default_str(setarg)) {
		libset_set_double(sv->key, NADBL);
		err = 0;
	    } else {
		double x;

		err = libset_get_scalar(sv->key, setarg, NULL, &x);
		if (!err) {
		    err = libset_set_double(sv->key, x);
		}
	    }
	} else if (sv->key == CSV_DELIM) {
	    char c = delim_from_arg(setarg);

	    if (c > 0) {
		data_delim = c;
		err = 0;
	    }
	} else if (sv->key == SEED) {
	    int automatic = 0;

	    err = libset_get_unsigned(setarg, &u, &automatic);
	    if (!err) {
		gretl_rand_set_seed(u);
		if (gretl_messages_on()) {
		    pprintf(prn,
			    _("Pseudo-random number generator seeded with %u\n"), u);
		}
		seed_is_set = !automatic;
	    }
	} else if (sv->key == HORIZON) {
	    /* horizon for VAR impulse responses */
	    if (!strcmp(setarg, "auto")) {
		state->horizon = UNSET_INT;
		err = 0;
	    } else {
		err = libset_get_scalar(sv->key, setarg, &k, NULL);
		if (!err) {
		    state->horizon = k;
		} else {
		    state->horizon = UNSET_INT;
		}
	    }
	} else if (coded_intvar(sv->key)) {
	    err = parse_libset_int_code(sv->key, setarg);
	} else if (libset_int(sv->key)) {
	    err = libset_get_scalar(sv->key, setarg, &k, NULL);
	    if (!err) {
		err = libset_set_int(sv->key, k);
	    }
	} else {
	    gretl_errmsg_sprintf(_("set: unknown variable '%s'"), setobj);
	    err = E_UNKVAR;
	}
    }

    return err;
}

double libset_get_double (SetKey key)
{
    void *valp;

    if (check_for_state()) {
	return NADBL;
    }

    valp = setkey_get_target(key, SV_DOUBLE);
    if (valp != NULL) {
	double x = *(double *) valp;

	if (na(x) && (key == NLS_TOLER || key == BFGS_TOLER)) {
	    x = get_default_nls_toler();
	}
	return x;
    } else {
	fprintf(stderr, "libset_get_double: unrecognized "
		"key %d\n", key);
	return 0;
    }
}

double libset_get_user_tolerance (SetKey key)
{
    if (key >= NLS_TOLER && key <= BHHH_TOLER) {
	void *valp = setkey_get_target(key, SV_ALL);

	return *(double *) valp;
    } else {
	return NADBL;
    }
}

int libset_set_double (SetKey key, double val)
{
    void *valp;
    int err = 0;

    if (check_for_state()) {
	return 1;
    }

    /* all the libset double vals must be positive */
    if (val <= 0.0 || na(val)) {
	return E_DATA;
    }

    valp = setkey_get_target(key, SV_DOUBLE);
    if (valp != NULL) {
	*(double *) valp = val;
    } else {
	fprintf(stderr, "libset_set_double: unrecognized key %d (%s)\n",
		key, setkey_get_name(key));
	err = E_UNKVAR;
    }

    return err;
}

int libset_get_int (SetKey key)
{
    void *valp;

    if (check_for_state()) {
	return 0;
    }

    valp = setkey_get_target(key, SV_INT);

    if (valp != NULL) {
#if SVDEBUG
	fprintf(stderr, "libset_get_int: valp %p\n", valp);
#endif
	if (libset_small_int(key)) {
	    return *(gint8 *) valp;
	} else {
	    return *(int *) valp;
	}
    } else if (key == BLAS_MNK_MIN) {
	return get_blas_mnk_min();
    } else if (key == OMP_N_THREADS) {
	return get_omp_n_threads();
    } else if (key == OMP_MNK_MIN) {
	return get_omp_mnk_min();
    } else if (key == SIMD_K_MAX) {
	return get_simd_k_max();
    } else if (key == SIMD_MN_MIN) {
	return get_simd_mn_min();
    } else {
	fprintf(stderr, "libset_get_int: unrecognized "
		"key %d\n", key);
	return 0;
    }
}

struct int_limits {
    SetKey key;
    int min;
    int max;
};

static int get_int_limits (SetKey key, int *min, int *max)
{
    static struct int_limits ilims[] = {
	{ HC_VERSION, 0, 4 },
	{ PANEL_ROBUST, 0, 2 },
	{ FDJAC_QUAL, 0, 2 },
	{ USE_QR, 0, 2 },
	{ LBFGS_MEM,  3, 20 },
	{ GRETL_DEBUG, 0, 4 },
	{ DATACOLS,    1, 15 },
	{ PLOT_COLLECT, 0, 2 },
	{ CSV_DIGITS, 1, 25 },
	{ BOOT_ITERS, 499, 999999 },
	{ BFGS_VERBSKIP, 0, 1000 },
	{ BOOTREP, 1, 99999 },
	{ HORIZON, 1, 1000 },
	{ LOOP_MAXITER, 0, INT_MAX - 1 },
	{ GMP_BITS, 256, 8192 }
    };
    int i;

    for (i=0; i<G_N_ELEMENTS(ilims); i++) {
	if (ilims[i].key == key) {
	    *min = ilims[i].min;
	    *max = ilims[i].max;
	    return 1;
	    break;
	}
    }

    return 0;
}

static int libset_int_min (SetKey key)
{
    int m1, m0 = 1;

    get_int_limits(key, &m0, &m1);
    return m0;
}

/* Called from within libset.c and also from various places in
   libgretl. It's primarily designed for "real" integer variables
   (not int-coded categories), but for now we make an exception
   for HC_VERSION and PLOT_COLLECT, to support existing calls
   from lib/src/estimate.c and gui/settings.c.
*/

int libset_set_int (SetKey key, int val)
{
    int err = 0;

    if (check_for_state()) {
	return 1;
    }

    if (key == BLAS_MNK_MIN) {
	set_blas_mnk_min(val);
    } else if (key == SIMD_K_MAX) {
	set_simd_k_max(val);
    } else if (key == SIMD_MN_MIN) {
	set_simd_mn_min(val);
    } else if (key == OMP_N_THREADS) {
	err = set_omp_n_threads(val);
    } else {
	int min = 1, max = 100000;
	void *valp;

	get_int_limits(key, &min, &max);
	if (val < min || val > max) {
	    err = E_DATA;
	} else {
	    valp = setkey_get_target(key, SV_INT);
	    if (valp == NULL) {
		err = E_DATA;
	    } else if (libset_small_int(key)) {
		*(gint8 *) valp = val;
	    } else {
		*(int *) valp = val;
	    }
	}
    }

    return err;
}

static void set_flag_from_env (SetKey flag, const char *s)
{
    char *e = getenv(s);

    if (e != NULL) {
	if (*e == '\0' || *e == '0') {
	    state->flags &= ~flag;
	} else {
	    state->flags |= flag;
	}
    }
}

static void set_int_from_env (SetKey key, const char *s)
{
    char *e = getenv(s);

    if (e != NULL) {
	if (*e == '\0' || *e == '0') {
	    libset_set_int(key, 0);
	} else {
	    libset_set_int(key, atoi(e));
	}
    }
}

static void maybe_check_env (SetKey key)
{
    if (key == USE_SVD) {
	set_flag_from_env(USE_SVD, "GRETL_USE_SVD");
    } else if (key == USE_LBFGS) {
	set_flag_from_env(USE_LBFGS, "GRETL_USE_LBFGS");
    } else if (key == USE_QR) {
	set_int_from_env(USE_QR, "GRETL_USE_QR");
    }
}

int libset_get_bool (SetKey key)
{
    /* global specials */
    if (key == R_FUNCTIONS) {
	return globals.R_functions;
    } else if (key == R_LIB) {
	return globals.R_lib;
    } else if (key == USE_DCMT) {
        return gretl_rand_get_dcmt();
    } else if (key == LOGSTAMP) {
	return globals.logstamp;
    }

    if (check_for_state()) {
	return 0;
    } else {
	maybe_check_env(key);
	return flag_to_bool(state, key);
    }
}

static void libset_set_decpoint (int on)
{
#ifdef ENABLE_NLS
    if (on) {
	/* force use of the decimal dot */
	setlocale(LC_NUMERIC, "C");
    } else {
	/* revert to whatever is the local default */
	char *current = get_built_in_string_by_name("lang");

	if (current != NULL && strcmp(current, "unknown")) {
	    setlocale(LC_NUMERIC, current);
	} else {
	    setlocale(LC_NUMERIC, "");
	}
    }

    reset_local_decpoint();
#endif
}

void set_data_export_decimal_comma (int s)
{
    if (s) {
	export_decpoint = ',';
    } else {
	export_decpoint = '.';
    }
}

char get_data_export_decpoint (void)
{
    char c = export_decpoint;

    /* revert to '.' on access */
    export_decpoint = '.';
    return c;
}

void set_data_export_delimiter (char c)
{
    data_delim = c;
}

char get_data_export_delimiter (void)
{
    return data_delim;
}

static int check_R_setting (gint8 *var, SetKey key, int val)
{
    int err = 0;

#ifdef USE_RLIB
    if (key == R_FUNCTIONS && val != 0) {
	/* R_functions depends on having R_lib on, so in case
	   the latter is off we should turn it on too.
	*/
	err = check_set_R_home();
	if (err) {
	    gretl_errmsg_set(_("R_functions could not be enabled\n"
                               "Please see https://gretl.sourceforge.net/gretl_and_R.html"));
	} else {
	    libset_set_bool(R_LIB, val);
	    *var = val;
	}
    } else {
	*var = val;
    }
#else
    if (val) {
	const char *s = (key == R_FUNCTIONS)? "R_functions" : "R_lib";

	gretl_errmsg_sprintf(_("%s: not supported in this build of gretl"), s);
	err = E_EXTERNAL;
    }
#endif

    return err;
}

int libset_set_bool (SetKey key, int val)
{
    if (check_for_state()) {
	return E_ALLOC;
    }

    /* global specials */
    if (key == R_FUNCTIONS) {
	return check_R_setting(&globals.R_functions, key, val);
    } else if (key == R_LIB) {
	return check_R_setting(&globals.R_lib, key, val);
    } else if (key == USE_DCMT) {
	return gretl_rand_set_dcmt(val);
    } else if (key == LOGSTAMP) {
	globals.logstamp = val;
	return 0;
    }

    if (val) {
	state->flags |= key;
    } else {
	state->flags &= ~key;
    }

    if (key == FORCE_DECPOINT) {
	libset_set_decpoint(val);
    }

    return 0;
}

/* Mechanism for pushing and popping program state for user-defined
   functions. push_program_state() is used when a function starts
   execution: the function gets a copy of the current program state,
   while that state is pushed onto the stack for restoration when the
   function exits.
*/

#define PPDEBUG 0

static int n_states;
static GPtrArray *state_stack;
static int state_idx = -1;

#if PPDEBUG
static void print_state_stack (int pop)
{
    set_state *sv;
    int i;

    fputs(pop ? "\nafter pop:\n" : "\nafter push:\n", stderr);
    for (i=0; i<n_states; i++) {
	sv = g_ptr_array_index(state_stack, i);
	fprintf(stderr, "%d: %p", i, (void *) sv);
	fputs(sv == state ? " *\n" : "\n", stderr);
    }
}
#endif

int push_program_state (void)
{
    set_state *newstate;
    int err = 0;

#if SVDEBUG
    fprintf(stderr, "push_program_state: n_states = %d\n", n_states);
#endif

    if (n_states == 0) {
	state_stack = g_ptr_array_new();
    }

    state_idx++;

    if (state_idx < n_states) {
	newstate = g_ptr_array_index(state_stack, state_idx);
    } else {
	newstate = malloc(sizeof *newstate);
	if (newstate == NULL) {
	    err = E_ALLOC;
	} else {
	    g_ptr_array_add(state_stack, newstate);
	    n_states++;
	}
    }

    if (newstate != NULL) {
	if (n_states == 1) {
	    state_vars_init(newstate);
	} else {
	    state_vars_copy(newstate);
	}
	state = newstate;
    }

#if SVDEBUG
    fprintf(stderr, " newstate = %p, state = %p\n",
	    (void *) newstate, (void *) state);
    fprintf(stderr, " gmm_maxiter at %p, value %d\n",
	    (void *) &(state->gmm_maxiter), state->gmm_maxiter);
#endif

#if PPDEBUG
    print_state_stack(0);
#endif

    return err;
}

static void free_state (set_state *sv)
{
    if (sv != NULL) {
	gretl_matrix_free(sv->initvals);
	gretl_matrix_free(sv->initcurv);
	gretl_matrix_free(sv->matmask);
	free(sv);
    }
}

/* Called when a user-defined function exits: restores the program
   state that was in force when the function started executing.
*/

int pop_program_state (void)
{
    int err = 0;

    if (n_states < 2) {
	err = 1;
    } else {
	int fdp0, fdp = state->flags & FORCE_DECPOINT;

	/* restore prior stack level */
	state = g_ptr_array_index(state_stack, --state_idx);

	fdp0 = state->flags & FORCE_DECPOINT;
	if (fdp0 != fdp) {
	    libset_set_decpoint(fdp0);
	}
	if (tex_plot_opts != NULL) {
	    maybe_destroy_tpo();
	}
    }

#if PPDEBUG
    print_state_stack(1);
#endif

    return err;
}

/* initialization of all user-settable settings */

int libset_init (void)
{
    static int done;
    int err = 0;

    if (!done) {
	err = push_program_state();
	done = 1;
    }

    return err;
}

/* state variables for looping */
static char *looping;
static int looplen;

void libset_cleanup (void)
{
    int i;

#if PDEBUG
    fprintf(stderr, "libset_cleanup called\n");
#endif

    for (i=0; i<n_states; i++) {
	free_state(g_ptr_array_index(state_stack, i));
    }

    g_ptr_array_free(state_stack, TRUE);
    state_stack = NULL;
    n_states = 0;
    state_idx = -1;

    libset_hash_cleanup();

    free(tex_plot_opts);
    tex_plot_opts = NULL;
    tpo_depth = TPO_UNSET;

    free(looping);
    looping = NULL;
    looplen = 0;
}

/* switches for looping and batch mode: output: these depend on the
   state of the program calling libgretl, they are not user-settable
*/

#define LDEBUG 0

#if LDEBUG

static void print_looping (int on)
{
    char c;
    int i;

    fputs(on ? "loop on:  " : "loop off: ", stderr);
    for (i=0; i<looplen; i++) {
	c = looping[i] ? '1' : '0';
	if (i == gretl_function_depth()) {
	    fprintf(stderr, "[%c]", c);
	} else {
	    fputc(c, stderr);
	}
    }
    fputc('\n', stderr);
}

#endif

#define LMIN 8

void set_loop_on (void)
{
    int fd = gretl_function_depth();

    if (looping == NULL) {
	looplen = fd+1 < LMIN ? LMIN : fd+1;
	looping = calloc(1, looplen);
    } else if (looplen < fd+1) {
	int n = looplen * 2;

	n = n < fd+1 ? fd+1 : n;
	looping = realloc(looping, n);
	memset(looping + looplen, 0, n - looplen);
	looplen = n;
    }
    looping[fd] = 1;
#if LDEBUG
    print_looping(1);
#endif
}

void set_loop_off (void)
{
    if (looping != NULL) {
	looping[gretl_function_depth()] = 0;
    }
#if LDEBUG
    print_looping(0);
#endif
}

/* returns 1 if there's a loop going on anywhere in the "caller
   ancestry" of the current execution level, else 0.
*/

int gretl_looping (void)
{
    int i, fd = gretl_function_depth();
    int n = MIN(fd+1, looplen);

    for (i=0; i<n; i++) {
	if (looping[i]) {
	    return 1;
	}
    }
    return 0;
}

/* returns 1 if there's a loop going on at the current execution
   stack level, else 0.
*/

int gretl_looping_currently (void)
{
    int fd = gretl_function_depth();

    return looplen >= fd+1 && looping[fd];
}

static int iter_depth;

void gretl_iteration_push (void)
{
    iter_depth++;
}

void gretl_iteration_pop (void)
{
    if (iter_depth > 0) {
	iter_depth--;
    }
}

int gretl_iteration_depth (void)
{
    return iter_depth;
}

/* portmanteau: a "loop" is active or (internal)
   iteration is in progress (e.g. in BFGS)
*/

int gretl_iterating (void)
{
    return gretl_looping() || iter_depth > 0;
}

static int batch_mode_switch (int set, int val)
{
    static int bmode;

    if (set) {
	bmode = val;
    }

    return bmode;
}

void gretl_set_batch_mode (int b)
{
    batch_mode_switch(1, b);
}

/* Returns 1 if we're running a script, otherwise 0.
   Note: a 0 return indicates that we're in an interactive
   mode, whether GUI or CLI.
*/

int gretl_in_batch_mode (void)
{
    return batch_mode_switch(0, 0);
}

static int gui_mode;

/* set by the GUI program at start-up */

void gretl_set_gui_mode (void)
{
    gui_mode = 1;
}

/* Returns 1 if we're running the GUI program. The current
   usage may be interactive (menu-driven or typing at the
   GUI "console") or script/batch. See also
   gretl_in_batch_mode().
*/

int gretl_in_gui_mode (void)
{
    return gui_mode;
}

/* "tool_mode" is set when gretlcli is being used as
   a build tool */

static int tool_mode;

void gretl_set_tool_mode (void)
{
    set_gretl_echo(0);
    set_gretl_messages(0);
    tool_mode = 1;
}

int gretl_in_tool_mode (void)
{
    return tool_mode;
}

/* mechanism to support callback for representing ongoing
   activity in the GUI */

static SHOW_ACTIVITY_FUNC sfunc;

void set_show_activity_func (SHOW_ACTIVITY_FUNC func)
{
    sfunc = func;
}

int show_activity_func_installed (void)
{
    return sfunc != NULL;
}

void show_activity_callback (void)
{
    if (sfunc != NULL) {
	(*sfunc)();
    }
}

/* Support for the GUI "Stop" button -- could in principle be extended
   to gretlcli via some Ctrl+key combination?
*/

static int user_stop;

void set_user_stop (int s)
{
    user_stop = s;
}

int get_user_stop (void)
{
    return user_stop;
}

/* internal: set the value of a string-valued "setvar" */

static int set_string_setvar (char *targ, const char *s, int len)
{
    *targ = '\0';

    if (*s == '"') {
	const char *p = strchr(s+1, '"');

	if (p == NULL) {
	    return E_PARSE;
	} else {
	    strncat(targ, s+1, p-s-1);
	}
    } else {
	strncat(targ, s, len);
    }

    return 0;
}

/* for setting what we print for NAs on CSV output */

const char *get_csv_na_write_string (void)
{
    if (check_for_state()) {
	return "NA";
    } else {
	return state->csv_write_na;
    }
}

int set_csv_na_write_string (const char *s)
{
    if (check_for_state()) {
	return E_DATA;
    } else {
	return set_string_setvar(state->csv_write_na, s, 7);
    }
}

/* and for setting what we read as NA on CSV input */

const char *get_csv_na_read_string (void)
{
    if (check_for_state()) {
	return "default";
    } else {
	return state->csv_read_na;
    }
}

int set_csv_na_read_string (const char *s)
{
    if (check_for_state()) {
	return E_DATA;
    } else {
	return set_string_setvar(state->csv_read_na, s, 7);
    }
}

int libset_write_script (const char *fname)
{
    PRN *prn;
    int err = 0;

    /* FIXME maybe adjust path for fname? */

    prn = gretl_print_new_with_filename(fname, &err);

    if (!err) {
	print_settings(prn, OPT_NONE);
	gretl_print_destroy(prn);
    }

    return err;
}

static char *get_quoted_arg (const char *s, int *err)
{
    const char *p;
    int n = 0, matched = 0;
    char *ret = NULL;

    p = s = strchr(s, '"') + 1;

    while (*p) {
	if (*p == '"') {
	    matched = 1;
	    break;
	} else {
	    n++;
	}
	p++;
    }

    if (!matched) {
	*err = E_PARSE;
    } else {
	ret = gretl_strndup(s, n);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

/* If @prn is non-NULL we're called via the "set" command.
   Otherwise we're called by the gretl session apparatus.
*/

static int real_libset_read_script (const char *fname,
				    PRN *prn)
{
    FILE *fp;
    int err = 0;

    fp = gretl_fopen(fname, "r");

    if (fp == NULL) {
	char fullname[FILENAME_MAX];

	strcpy(fullname, fname);
	gretl_addpath(fullname, 0);
	fp = gretl_fopen(fname, "r");
	if (fp == NULL) {
	    err = E_FOPEN;
	}
    }

    if (!err) {
	char setobj[32], setarg[32], line[1024];
	int nf;

	while (fgets(line, sizeof line, fp)) {
	    if (*line == '#' || string_is_blank(line)) {
		continue;
	    }
	    tailstrip(line);
	    nf = sscanf(line, "%*s %31s %31s", setobj, setarg);
	    if (nf == 1) {
		err = execute_set(setobj, NULL, NULL, OPT_NONE, prn);
	    } else if (nf == 2) {
		if (*setarg == '"') {
		    char *q = get_quoted_arg(line, &err);

		    if (!err) {
			err = execute_set(setobj, q, NULL, OPT_NONE, prn);
			free(q);
		    }
		} else {
		    err = execute_set(setobj, setarg, NULL, OPT_NONE, prn);
		}
	    }
	    if (err && prn != NULL) {
		break;
	    }
	}

	fclose(fp);
    }

    return err;
}

int libset_read_script (const char *fname)
{
    return real_libset_read_script(fname, NULL);
}
