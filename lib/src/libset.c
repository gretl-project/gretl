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
#include "gretl_string_table.h"
#include "gretl_func.h"

#ifdef _OPENMP
# include <omp.h>
#endif

#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

#include <unistd.h>
#include <errno.h>

#include <glib.h>

#define PDEBUG 0

enum {
    AUTO_LAG_STOCK_WATSON,
    AUTO_LAG_WOOLDRIDGE,
    AUTO_LAG_NEWEYWEST
};

/* state flags */

enum {
    STATE_USE_CWD         = 1 << 0,  /* store: use current dir as default */
    STATE_ECHO_ON         = 1 << 1,  /* echoing commands or not */
    STATE_MSGS_ON         = 1 << 2,  /* emitting non-error messages or not */
    STATE_FORCE_DECPOINT  = 1 << 3,  /* override locale decimal separator */
    STATE_USE_PCSE        = 1 << 4,  /* Beck-Katz panel-corrected std errs */
    STATE_USE_SVD         = 1 << 5,  /* SVD decomposition is matrix OLS default */
    STATE_PREWHITEN       = 1 << 6,  /* HAC pre-whitening? */
    STATE_FORCE_HC        = 1 << 7,  /* don't use HAC for time series */
    STATE_HALT_ON_ERR     = 1 << 8,  /* errors fatal in batch mode */
    STATE_USE_LBFGS       = 1 << 9,  /* prefer LBFGS to BFGS? */
    STATE_SHELL_OK        = 1 << 10, /* "shell" facility is approved? */
    STATE_MAX_VERBOSE     = 1 << 11, /* verbose output from maximizer? */
    STATE_USE_FCP         = 1 << 12, /* use FCP garch code */
    STATE_WARN_ON         = 1 << 13, /* print numerical warning messages */
    STATE_SKIP_MISSING    = 1 << 14, /* skip NAs when building matrix from series */
    STATE_LOOPING         = 1 << 15, /* loop is in progress at this level */
    STATE_LOOP_QUIET      = 1 << 16, /* loop commands should be quiet */
    STATE_LOOP_PROG       = 1 << 17, /* progressive loop is in progress */
    STATE_BFGS_RSTEP      = 1 << 18, /* use Richardson method in BFGS numerical
					gradient */
    STATE_DPDSTYLE_ON     = 1 << 19, /* emulate dpd in dynamic panel data models */
    STATE_OPENMP_ON       = 1 << 20  /* using openmp */
};    

/* for values that really want a non-negative integer */
#define UNSET_INT -9
#define is_unset(i) (i == UNSET_INT)

typedef struct set_vars_ set_vars;

struct robust_opts {
    int auto_lag;
    int user_lag;
    int hc_version;
    int hkern;
    double qsband;
};

struct set_vars_ {
    int flags;
    unsigned int seed;          /* for PRNG */
    double conv_huge;           /* conventional value for $huge */
    int horizon;                /* for VAR impulse responses */ 
    int bootrep;                /* bootstrap replications */
    double nls_toler;           /* NLS convergence criterion */
    int loop_maxiter;           /* max no. of iterations in non-for loops */
    int vecm_norm;              /* VECM beta normalization */
    int optim;                  /* code for preferred optimizer */
    int bfgs_maxiter;           /* max iterations, BFGS */         
    double bfgs_toler;          /* convergence tolerance, BFGS */
    double bfgs_maxgrad;        /* max acceptable gradient norm, BFGS */
    int bfgs_verbskip;          /* BFGS: show one in n iterations  */
    int optim_steplen;          /* step length algorithm (only BFGS for now) */
    int bhhh_maxiter;           /* max iterations, BHHH */          
    double bhhh_toler;          /* convergence tolerance, BHHH */
    int lbfgs_mem;              /* memory span for L-BFGS-B */
    int garch_vcv;              /* GARCH vcv variant */
    int garch_robust_vcv;       /* GARCH vcv variant, robust estimation */
    int arma_vcv;               /* ARMA vcv variant */
    int rq_maxiter;             /* max iterations for quantreg, simplex */
    int gmm_maxiter;            /* max iterations for iterated GMM */
    gretl_matrix *initvals;     /* parameter initializer */
    gretl_matrix *matmask;      /* mask for series -> matrix conversion */
    struct robust_opts ropts;   /* robust standard error options */
    char shelldir[MAXLEN];      /* working dir for shell commands */
    char csv_write_na[8];       /* representation of NA in CSV output */
    char csv_read_na[8];        /* representation of NA in CSV input */
    double nadarwat_trim;       /* multiple of h to use in nadarwat() for trimming */
    int fdjac_qual;             /* quality of "fdjac" function */
};

#define ECHO "echo"
#define MESSAGES "messages"
#define WARNINGS "warnings"
#define GRETL_DEBUG "debug"
#define USE_DCMT "use_dcmt"

#define BLAS_MNK_MIN "blas_mnk_min"
#define SIMD_K_MAX "simd_k_max"
#define SIMD_MN_MIN "simd_mn_min"
#define MP_MNK_MIN "mp_mnk_min"
/* but now preferred */
#define OMP_MNK_MIN "omp_mnk_min"
#define OMP_N_THREADS "omp_num_threads"

#define libset_boolvar(s) (!strcmp(s, ECHO) || \
                           !strcmp(s, MESSAGES) || \
                           !strcmp(s, WARNINGS) || \
                           !strcmp(s, FORCE_DECP) || \
			   !strcmp(s, FORCE_HC) || \
			   !strcmp(s, HALT_ON_ERR) || \
                           !strcmp(s, MAX_VERBOSE) || \
			   !strcmp(s, USE_LBFGS) || \
			   !strcmp(s, PCSE) || \
			   !strcmp(s, PREWHITEN) || \
			   !strcmp(s, USE_SVD) || \
			   !strcmp(s, SHELL_OK) || \
			   !strcmp(s, USE_CWD) || \
			   !strcmp(s, USE_FCP) || \
                           !strcmp(s, SKIP_MISSING) || \
			   !strcmp(s, R_FUNCTIONS) || \
			   !strcmp(s, R_LIB) || \
			   !strcmp(s, BFGS_RSTEP) || \
			   !strcmp(s, DPDSTYLE) || \
			   !strcmp(s, USE_DCMT) || \
			   !strcmp(s, USE_OPENMP))

#define libset_double(s) (!strcmp(s, CONV_HUGE) || \
			  !strcmp(s, BFGS_TOLER) || \
			  !strcmp(s, BFGS_MAXGRAD) || \
			  !strcmp(s, BHHH_TOLER) || \
			  !strcmp(s, NLS_TOLER) || \
			  !strcmp(s, QS_BANDWIDTH) || \
			  !strcmp(s, NADARWAT_TRIM))

#define libset_int(s) (!strcmp(s, BFGS_MAXITER) || \
		       !strcmp(s, BFGS_VERBSKIP) || \
		       !strcmp(s, OPTIM_STEPLEN) || \
		       !strcmp(s, BHHH_MAXITER) || \
		       !strcmp(s, GMM_MAXITER) || \
		       !strcmp(s, LBFGS_MEM) || \
		       !strcmp(s, BOOTREP) || \
		       !strcmp(s, HAC_KERNEL) || \
                       !strcmp(s, HC_VERSION) || \
		       !strcmp(s, HORIZON) || \
		       !strcmp(s, LOOP_MAXITER) || \
                       !strcmp(s, RQ_MAXITER) || \
		       !strcmp(s, VECM_NORM) || \
		       !strcmp(s, GRETL_OPTIM) || \
		       !strcmp(s, GRETL_DEBUG) || \
		       !strcmp(s, BLAS_MNK_MIN) || \
		       !strcmp(s, OMP_MNK_MIN) || \
		       !strcmp(s, MP_MNK_MIN) || \
		       !strcmp(s, OMP_N_THREADS) || \
		       !strcmp(s, SIMD_K_MAX) || \
		       !strcmp(s, SIMD_MN_MIN) || \
		       !strcmp(s, FDJAC_QUAL))

/* global state */
set_vars *state;
static int gretl_debug;
static int user_mp_bits;
static int R_functions;
static int R_lib = 1;
static int csv_digits = UNSET_INT;
static char data_delim = ',';
static char data_export_decpoint = '.';
#ifdef OS_OSX
static int omp_mnk_min = -1;
#else
static int omp_mnk_min = 80000; /* was 65535 */
#endif
static int omp_n_threads;

static int boolvar_get_flag (const char *s);
static const char *hac_lag_string (void);
static int real_libset_read_script (const char *fname,
				    PRN *prn);
static int set_csv_digits (const char *s);

static void robust_opts_init (struct robust_opts *r)
{
    r->auto_lag = AUTO_LAG_STOCK_WATSON;
    r->user_lag = UNSET_INT;
    r->hc_version = 0;
    r->hkern = KERNEL_BARTLETT;
    r->qsband = NADBL;
}

static void robust_opts_copy (struct robust_opts *r)
{
    r->auto_lag = state->ropts.auto_lag;
    r->user_lag = state->ropts.user_lag;
    r->hc_version = state->ropts.hc_version;
    r->hkern = state->ropts.hkern; 
    r->qsband = state->ropts.qsband;
}

static const char *csv_delim_args[] = {
    "comma",
    "space",
    "tab",
    "semicolon",
    NULL
};

static const char *garch_vcv_strs[] = {
    "unset",
    "hessian",
    "im",
    "op",
    "qml",
    "bw",
    NULL
};

static const char *arma_vcv_strs[] = {
    "hessian",
    "op",
    NULL
};

static const char *hac_kernel_strs[] = {
    "bartlett", 
    "parzen", 
    "qs",
    NULL
};

static const char *hc_version_strs[] = {
    "0", "1", "2", "3", "3a", NULL
};

static const char *vecm_norm_strs[] = {
    "phillips",
    "diag",
    "first",
    "none",
    NULL
};

static const char *optim_strs[] = {
    "auto",
    "BFGS",
    "newton",
    NULL
};

static const char *steplen_strs[] = {
    "power",
    "quadratic",
    NULL
};

static const char *normal_rand_strs[] = {
    "ziggurat",
    "box-muller",
    NULL
};

static const char **libset_option_strings (const char *s)
{
    if (!strcmp(s, GARCH_VCV)) {
	return garch_vcv_strs;
    } else if (!strcmp(s, ARMA_VCV)) {
	return arma_vcv_strs;
    } else if (!strcmp(s, HAC_KERNEL)) {
	return hac_kernel_strs;
    } else if (!strcmp(s, HC_VERSION)) {
	return hc_version_strs;
    } else if (!strcmp(s, VECM_NORM)) {
	return vecm_norm_strs;
    } else if (!strcmp(s, GRETL_OPTIM)) {
	return optim_strs;
    } else if (!strcmp(s, NORMAL_RAND)) {
	return normal_rand_strs;
    } else if (!strcmp(s, "csv_delim")) {
	return csv_delim_args;
    } else if (!strcmp(s, OPTIM_STEPLEN)) {
	return steplen_strs;
    } else {
	return NULL;
    }
}

static void coded_var_show_opts (const char *s, PRN *prn)
{
    const char **S = libset_option_strings(s);

    if (S != NULL) {
	pputs(prn, "valid settings:");
	while (*S != NULL) {
	    pprintf(prn, " %s", *S);
	    S++;
	}
	pputc(prn, '\n');
    }
}

static const char *get_arma_vcv_str (int v)
{
    if (v == ML_HESSIAN) {
	return arma_vcv_strs[0];
    } else if (v == ML_OP) {
	return arma_vcv_strs[1];
    } else {
	return "unknown";
    }
}

static const char *libset_option_string (const char *s)
{
    if (!strcmp(s, HAC_LAG)) {
	return hac_lag_string(); /* special */
    } else if (!strcmp(s, GARCH_VCV)) {
	return garch_vcv_strs[state->garch_vcv];
    } else if (!strcmp(s, ARMA_VCV)) {
	return get_arma_vcv_str(state->arma_vcv);
    } else if (!strcmp(s, HAC_KERNEL)) {
	return hac_kernel_strs[state->ropts.hkern];
    } else if (!strcmp(s, HC_VERSION)) {
	return hc_version_strs[state->ropts.hc_version];
    } else if (!strcmp(s, VECM_NORM)) {
	return vecm_norm_strs[state->vecm_norm];
    } else if (!strcmp(s, GRETL_OPTIM)) {
	return optim_strs[state->optim];
    } else if (!strcmp(s, NORMAL_RAND)) {
	return normal_rand_strs[gretl_rand_get_box_muller()];
    } else if (!strcmp(s, OPTIM_STEPLEN)) {
	return steplen_strs[state->optim_steplen];
    } else {
	return "?";
    }
}

static void print_initvals (const gretl_matrix *ivals, PRN *prn,
			    gretlopt opt)
{
    if (opt & OPT_D) {
	if (ivals == NULL) {
	    pputs(prn, " initvals = auto\n");
	} else {
	    gretl_matrix_print_to_prn(ivals, " initvals =", prn);
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

static int flag_to_bool (set_vars *sv, int flag)
{
    if (!sv) {
	return 0;
    } else {
	return (sv->flags & flag)? 1 : 0;
    }
}

static void state_vars_copy (set_vars *sv)
{
#if PDEBUG
    fprintf(stderr, "state_vars_copy() called\n");
#endif
    sv->flags = state->flags;
    /* We're not (yet) looping at the current level of execution (but
       note that the STATE_LOOP_QUIET flag should be inherited).
    */
    sv->flags &= ~STATE_LOOPING;

    sv->seed = state->seed;
    sv->conv_huge = state->conv_huge;
    sv->horizon = state->horizon;
    sv->bootrep = state->bootrep;
    sv->loop_maxiter = state->loop_maxiter;
    sv->rq_maxiter = state->rq_maxiter;
    sv->gmm_maxiter = state->gmm_maxiter;
    sv->nls_toler = state->nls_toler;
    sv->vecm_norm = state->vecm_norm;
    sv->optim = state->optim;
    sv->bfgs_maxiter = state->bfgs_maxiter;
    sv->bfgs_toler = state->bfgs_toler;
    sv->bfgs_maxgrad = state->bfgs_maxgrad;
    sv->bfgs_verbskip = state->bfgs_verbskip;
    sv->optim_steplen = state->optim_steplen;
    sv->bhhh_maxiter = state->bhhh_maxiter;
    sv->bhhh_toler = state->bhhh_toler;
    sv->lbfgs_mem = state->lbfgs_mem;
    sv->garch_vcv = state->garch_vcv;
    sv->arma_vcv = state->arma_vcv;
    sv->garch_robust_vcv = state->garch_robust_vcv;
    sv->nadarwat_trim = state->nadarwat_trim;
    sv->fdjac_qual = state->fdjac_qual;

    sv->initvals = gretl_matrix_copy(state->initvals);
    sv->matmask = gretl_matrix_copy(state->matmask);
    strcpy(sv->shelldir, state->shelldir);
    strcpy(sv->csv_write_na, state->csv_write_na);
    strcpy(sv->csv_read_na, state->csv_read_na);

    robust_opts_copy(&sv->ropts);
}

/* for processors count */
#if defined(WIN32)
# include <windows.h>
#elif defined(OS_OSX)
# include <sys/param.h>
# include <sys/sysctl.h>
#endif

int gretl_n_processors (void)
{
    static int n_proc = -1;

    if (n_proc >= 1) {
	return n_proc;
    }

    n_proc = 1;

#if defined(WIN32)
    SYSTEM_INFO sysinfo;

    GetSystemInfo(&sysinfo);
    n_proc = sysinfo.dwNumberOfProcessors;
#elif defined(OS_OSX)
    int mib[2] = {CTL_HW, HW_NCPU}; /* or CTL_KERN, KERN_MAXPROC ? */
    size_t len = sizeof n_proc;

    if (sysctl(mib, 2, &n_proc, &len, NULL, 0) == -1) {
	perror("could not determine number of CPUs available");
	n_proc = 1;
    }    
#else
    n_proc = sysconf(_SC_NPROCESSORS_ONLN);
#endif

    return n_proc;
}

#define OMP_SHOW 0

int libset_use_openmp (guint64 n)
{
#if defined(_OPENMP)
    if (state == NULL || !(state->flags & STATE_OPENMP_ON)) {
	return 0;
    } else if (omp_mnk_min >= 0 && n >= (guint64) omp_mnk_min) {
# if OMP_SHOW > 1
	fprintf(stderr, "libset_use_openmp: yes\n");
# endif
	return 1;
    }
#endif

    return 0;
}

static int set_omp_n_threads (int n)
{
#if defined(_OPENMP)
    if (n < 1 || n > gretl_n_processors()) {
	gretl_errmsg_sprintf("omp_num_threads: must be >= 1 and <= %d",
			     gretl_n_processors());
	return E_DATA;
    } else {
	omp_n_threads = n;
	omp_set_num_threads(n);
	return 0;
    }
#else
    gretl_errmsg_set("OpenMP is not enabled");
    return E_DATA;
#endif
}

int get_omp_n_threads (void)
{
    return omp_n_threads;
}

#if defined(_OPENMP)

static int openmp_by_default (void)
{
    static int called = 0;
    int num_cores = gretl_n_processors();
    int ret = num_cores > 1;

    if (ret) {
	/* one can use the environment to turn this off */
	char *envstr = getenv("GRETL_USE_OPENMP");

	if (envstr != NULL && !strcmp(envstr, "0")) {
	    ret = 0;
	}
    }

    if (!called && ret) {
	char *s = getenv("OMP_NUM_THREADS");

	if (s != NULL && *s != '\0') {
	    omp_n_threads = atoi(s);
	} else {
	    omp_n_threads = num_cores;
	}
	called = 1;
    }

# if OMP_SHOW
    if (1) {
	fprintf(stderr, "number of cores detected = %d\n", num_cores);
	fprintf(stderr, "use OpenMP by default? %s\n", ret ? "yes" : "no");
	fprintf(stderr, "omp_num_threads = %d\n", omp_n_threads);
    }
# endif	

    return ret;
}

#endif /* _OPENMP defined */

static void state_vars_init (set_vars *sv)
{
#if PDEBUG
    fprintf(stderr, "state_vars_init called\n");
#endif
    sv->flags = STATE_ECHO_ON | STATE_MSGS_ON | STATE_WARN_ON | 
	STATE_HALT_ON_ERR | STATE_SKIP_MISSING;
#if defined(_OPENMP)
    if (openmp_by_default()) {
	sv->flags |= STATE_OPENMP_ON;
    }
#endif
    sv->seed = 0;
    sv->conv_huge = 1.0e100;
    sv->horizon = UNSET_INT;
    sv->bootrep = 1000;
    sv->nls_toler = NADBL;
    sv->loop_maxiter = 100000;
    sv->rq_maxiter = 1000;
    sv->gmm_maxiter = 250;
    sv->vecm_norm = NORM_PHILLIPS;
    sv->optim = OPTIM_AUTO;
    sv->initvals = NULL;
    sv->matmask = NULL;

    sv->bfgs_maxiter = UNSET_INT;
    sv->bfgs_toler = NADBL;
    sv->bfgs_maxgrad = 5.0;
    sv->bfgs_verbskip = 1;
    sv->optim_steplen = STEPLEN_POWER;
    sv->bhhh_maxiter = 500;
    sv->bhhh_toler = NADBL;
    sv->lbfgs_mem = 8;
    sv->garch_vcv = ML_UNSET;
    sv->arma_vcv = ML_HESSIAN;
    sv->garch_robust_vcv = ML_UNSET;
    sv->nadarwat_trim = 4.0;
    sv->fdjac_qual = 0;

    *sv->shelldir = '\0';

    strcpy(sv->csv_write_na, "NA");
    strcpy(sv->csv_read_na, "default");

    robust_opts_init(&sv->ropts);
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

int gretl_warnings_on (void)
{
    if (check_for_state()) return 1;
    return flag_to_bool(state, STATE_WARN_ON);
}

int gretl_debugging_on (void)
{
    return gretl_debug;
}

#define DEFAULT_MP_BITS 256
#define mp_bits_ok(b) (b >= 256 && b <= 8192)

void set_mp_bits (int b)
{
    if (mp_bits_ok(b)) {
	user_mp_bits = b;
    }
}

int get_mp_bits (void)
{
    if (user_mp_bits > DEFAULT_MP_BITS) {
	return user_mp_bits;
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
    int h = 0;

    check_for_state();

    /* Variants of Newey-West */

    if (state->ropts.user_lag >= 0 && state->ropts.user_lag < T - 2) {
	/* FIXME upper limit? */
	h = state->ropts.user_lag;
    } else if (state->ropts.auto_lag == AUTO_LAG_WOOLDRIDGE) {
	h = 4.0 * pow(T / 100.0, 2.0 / 9.0);
    } else {
	/* Stock-Watson default */
	h = 0.75 * pow(T, 1.0 / 3.0);
    }

    return h;
}

/* prewhitening implies nw3, but not vice versa */

int data_based_hac_bandwidth (void)
{
    if (is_unset(state->ropts.user_lag)) {
	if (state->ropts.auto_lag == AUTO_LAG_NEWEYWEST ||
	    (state->flags & STATE_PREWHITEN)) {
	    return 1;
	}
    }

    return 0;
}

static const char *hac_lag_string (void)
{
    check_for_state();

    if (state->ropts.user_lag >= 0 && state->ropts.user_lag < 1000) {
	static char lagstr[6];

	sprintf(lagstr, "%d", state->ropts.user_lag);
	return lagstr;
    } else if (state->ropts.auto_lag == AUTO_LAG_STOCK_WATSON) {
	return "nw1";
    } else {
	return "nw2";
    }
}

/* set max lag for HAC estimation */

static int parse_hac_lag_variant (const char *s)
{
    int err = E_DATA;

    if (!strcmp(s, "nw1")) {
	state->ropts.auto_lag = AUTO_LAG_STOCK_WATSON;
	state->ropts.user_lag = UNSET_INT;
	err = 0;
    } else if (!strcmp(s, "nw2")) {
	state->ropts.auto_lag = AUTO_LAG_WOOLDRIDGE;
	state->ropts.user_lag = UNSET_INT;
	err = 0;
    } else if (!strcmp(s, "nw3") ||
	       !strcmp(s, "auto")) {
	state->ropts.auto_lag = AUTO_LAG_NEWEYWEST;
	state->ropts.user_lag = UNSET_INT;
	err = 0;
    } else if (isdigit(*s)) {
	state->ropts.user_lag = atoi(s);
	err = 0;
    }

    return err;
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
	    gretl_errmsg_set_from_errno(s);
	    *err = 1;
	}
    } else {
	long li = strtol(s, &test, 10);

	if (*test != '\0') {
	    ret = 0;
	} else if (errno == ERANGE) {
	    gretl_errmsg_set_from_errno(s);
	    *err = 1;
	} else {
	    *pi = (int) li;
	}
    }

    gretl_pop_c_numeric_locale();

    return ret;
}

static int negval_invalid (const char *var)
{
    int ret = 1; /* presume invalid */

    if (var != NULL) {
	if (!strcmp(var, BLAS_MNK_MIN) ||
	    !strcmp(var, OMP_MNK_MIN) ||
	    !strcmp(var, MP_MNK_MIN) ||
	    !strcmp(var, SIMD_K_MAX) ||
	    !strcmp(var, SIMD_MN_MIN)) {
	    /* these can all be set to -1 */
	    ret = 0;
	}
    }

    return ret;
}

static int libset_get_scalar (const char *var, const char *arg, 
			      int *pi, double *px)
{
    double x = NADBL;
    int err = 0;

    if (libset_numeric_string(arg, pi, px, &err)) {
	if (err) {
	    err = E_DATA;
	} else if (pi != NULL && negval_invalid(var) && *pi < 0) {
	    err = E_DATA;
	} else if (px != NULL && *px < 0.0) {
	    err = E_DATA;
	}
	return err;
    }

    if (gretl_is_scalar(arg)) {
	x = gretl_scalar_get_value(arg, NULL);
    } else {
	gretl_errmsg_sprintf("'%s': not a scalar", arg);
	return E_UNKVAR;
    }

    if (negval_invalid(var) && x < 0.0) {
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
    int i;

    check_for_state();

    if (!strncmp(s, "hc", 2)) {
	s += 2;
    }

    for (i=0; hc_version_strs[i] != NULL; i++) {
	if (!strcmp(s, hc_version_strs[i])) {
	    state->ropts.hc_version = i;
	    return 0;
	}
    }

    if (!strcmp(s, "4")) {
	state->ropts.hc_version = 4;
	return 0;
    }

    return 1;
}

static int parse_libset_int_code (const char *key, 
				  const char *val)
{
    int i, err = E_DATA;

    if (!g_ascii_strcasecmp(key, HC_VERSION)) {
	err = parse_hc_variant(val);
    } else if (!g_ascii_strcasecmp(key, HAC_LAG)) {
	err = parse_hac_lag_variant(val);
    } else if (!g_ascii_strcasecmp(key, GARCH_VCV)) {
	for (i=0; i<ML_VCVMAX; i++) {
	    if (!g_ascii_strcasecmp(val, garch_vcv_strs[i])) {
		state->garch_vcv = i;
		err = 0;
		break;
	    }
	}
    } else if (!g_ascii_strcasecmp(key, ARMA_VCV)) {
	if (!g_ascii_strcasecmp(val, "op")) {
	    state->arma_vcv = ML_OP;
	    err = 0;
	} else if (!g_ascii_strcasecmp(val, "hessian")) {
	    state->arma_vcv = ML_HESSIAN;
	    err = 0;
	}
    } else if (!g_ascii_strcasecmp(key, HAC_KERNEL)) {
	for (i=0; i<KERNEL_MAX; i++) {
	    if (!g_ascii_strcasecmp(val, hac_kernel_strs[i])) {
		state->ropts.hkern = i;
		err = 0;
		break;
	    }
	}
    } else if (!g_ascii_strcasecmp(key, VECM_NORM)) {
	for (i=0; i<NORM_MAX; i++) {
	    if (!g_ascii_strcasecmp(val, vecm_norm_strs[i])) {
		state->vecm_norm = i;
		err = 0;
		break;
	    }
	}
    } else if (!g_ascii_strcasecmp(key, GRETL_OPTIM)) {
	for (i=0; i<OPTIM_MAX; i++) {
	    if (!g_ascii_strcasecmp(val, optim_strs[i])) {
		state->optim = i;
		err = 0;
		break;
	    }
	}	
    } else if (!g_ascii_strcasecmp(key, NORMAL_RAND)) {
	for (i=0; normal_rand_strs[i] != NULL; i++) {
	    if (!g_ascii_strcasecmp(val, normal_rand_strs[i])) {
		gretl_rand_set_box_muller(i);
		err = 0;
		break;
	    }
	}
    } else if (!g_ascii_strcasecmp(key, OPTIM_STEPLEN)) {
	for (i=0; i<STEPLEN_MAX; i++) {
	    if (!g_ascii_strcasecmp(val, steplen_strs[i])) {
		state->optim_steplen = i;
		err = 0;
		break;
	    }
	}	
    }

    if (err) {
	gretl_errmsg_sprintf(_("%s: invalid value '%s'"), key, val);
    }

    return err;
}	

void set_xsect_hccme (const char *s)
{
    char *scpy;

    if (check_for_state()) return;

    scpy = gretl_strdup(s);

    if (scpy != NULL) {
	gretl_lower(scpy);
	parse_hc_variant(scpy);
	free(scpy);
    }
}

void set_tseries_hccme (const char *s)
{
    char *scpy;

    if (check_for_state()) return;

    scpy = gretl_strdup(s);

    if (scpy != NULL) {
	gretl_lower(scpy);
	if (parse_hc_variant(scpy) == 0) {
	    libset_set_bool(FORCE_HC, 1);
	} else {
	    libset_set_bool(FORCE_HC, 0);
	}
	free(scpy);
    }
}

void set_panel_hccme (const char *s)
{
    if (check_for_state()) return;

    if (!strcmp(s, "Arellano")) {
	state->flags &= ~STATE_USE_PCSE;
    } else if (!strcmp(s, "PCSE")) {
	state->flags |= STATE_USE_PCSE;
    }
}

void set_garch_robust_vcv (const char *s)
{
    char *scpy;

    if (check_for_state()) return;

    scpy = gretl_strdup(s);

    if (scpy != NULL) {
	gretl_lower(scpy);
	if (!strcmp(s, "qml")) {
	    state->garch_robust_vcv = ML_QML;
	} else if (!strcmp(s, "bw")) {
	    state->garch_robust_vcv = ML_BW;
	}
	free(scpy);
    }
}

static int set_initvals (const char *mname, PRN *prn)
{
    gretl_matrix *m;
    int err = 0;

    if (!strcmp(mname, "auto")) {
	gretl_matrix_free(state->initvals);
	state->initvals = NULL;
    } else {
	m = get_matrix_by_name(mname);
	if (m == NULL) {
	    pprintf(prn, _("'%s': no such matrix"), mname);
	    pputc(prn, '\n');
	    err = E_DATA;
	} else {
	    if (state->initvals != NULL) {
		gretl_matrix_free(state->initvals);
	    }
	    state->initvals = gretl_matrix_copy(m);
	    if (state->initvals == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    return err;
}

static int set_matmask (const char *vname, const DATASET *dset,
			PRN *prn)
{
    int err = 0;

    if (!strcmp(vname, "null")) {
	gretl_matrix_free(state->matmask);
	state->matmask = NULL;
    } else {
	int t, v = current_series_index(dset, vname);

	if (v < 0) {
	    err = E_UNKVAR;
	} else {
	    if (state->matmask != NULL) {
		gretl_matrix_free(state->matmask);
	    }
	    state->matmask = gretl_column_vector_alloc(dset->n);
	    if (state->matmask == NULL) {
		err = E_ALLOC;
	    } else {
		for (t=0; t<dset->n; t++) {
		    state->matmask->val[t] = dset->Z[v][t];
		}
	    }
	}
    }

    return err;
}

void destroy_matrix_mask (void)
{
    check_for_state();
    gretl_matrix_free(state->matmask);
    state->matmask = NULL;
}

void shelldir_init (const char *s)
{
    if (s != NULL) {
	int n;

	*state->shelldir = '\0';
	strncat(state->shelldir, s, MAXLEN - 1);
	n = strlen(state->shelldir);
	if (n > 0 && (state->shelldir[n-1] == '\\' ||
		      state->shelldir[n-1] == '/')) {
	    state->shelldir[n-1] = '\0';
	}	
    } else {
	char *test = getcwd(state->shelldir, MAXLEN);

	if (test == NULL) {
	    *state->shelldir = '\0';
	} 
    }

    gretl_insert_builtin_string("shelldir", state->shelldir);
}

static int set_shelldir (const char *s)
{
    int len = 0, err = 0;

    /* skip past "set shelldir" and space */
    s += 12;
    s += strspn(s, " ");

    if (*s == '\0') {
	*state->shelldir = '\0';
	gretl_insert_builtin_string("shelldir", state->shelldir);
    } else if (*s == '"') {
	s++;
	len = gretl_charpos('"', s);
	if (len <= 0) {
	    err = E_PARSE;
	} 
    } else {
	len = strlen(s);
    }

    if (!err && len > 0) {
	char test[MAXLEN];
	char *home = NULL;
	int slen = len;

	if (*s == '~') {
	    home = getenv("HOME");
	    if (home != NULL) {
		s++;
		slen--;
		len = slen + strlen(home);
	    }
	} 

	*test = '\0';
    
	if (len >= MAXLEN) {
	    gretl_errmsg_set("shelldir: string is too long");
	    err = E_DATA;
	} else if (home != NULL) {
	    strcat(test, home);
	    strncat(test, s, slen);
	} else {
	    strncat(test, s, len);
	}

	if (!gretl_isdir(test)) {
	    gretl_errmsg_sprintf("shelldir: '%s' no such directory", test);
	    err = E_DATA;
	}

	if (!err) {
	    strcpy(state->shelldir, test);
	    gretl_insert_builtin_string("shelldir", state->shelldir);
	}
    }

    return err;
}

static int (*workdir_callback)();

void set_workdir_callback (int (*callback)())
{
    workdir_callback = callback;
}

static int set_workdir (const char *s)
{
    int err = 0;

    if (gretl_function_depth() > 0) {
	gretl_errmsg_set("set workdir: cannot be done inside a function");
	return 1;
    }

    /* skip past "set workdir" and space */
    s += 11;
    s += strspn(s, " ");

    if (*s == '\0') {
	err = E_DATA;
    } else {
	char workdir[MAXLEN];

	if (*s == '"') {
	    sscanf(s+1, "%511[^\"]", workdir);
	} else {
	    sscanf(s, "%511s", workdir);
	}
	if (workdir_callback != NULL) {
	    err = (*workdir_callback)(workdir);
	} else {
	    err = set_gretl_work_dir(workdir);
	}
    } 

    return err;
}

const char *csv_delims = ", \t;";

static char delim_from_arg (const char *s)
{
    int i;

    for (i=0; csv_delim_args[i] != NULL; i++) {
	if (!strcmp(s, csv_delim_args[i])) {
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
	    return csv_delim_args[i];
	}
    }

    return "unset";
}

static void libset_print_bool (const char *s, PRN *prn,
			       gretlopt opt)
{
    int v = libset_get_bool(s);

    if (opt & OPT_D) {
	pprintf(prn, " %s = %d\n", s, v);
    } else {
	pprintf(prn, "set %s %s\n", s, v? "on" : "off");
    }
}

#define coded_intvar(s) (!strcmp(s, GARCH_VCV) || \
			 !strcmp(s, ARMA_VCV) || \
			 !strcmp(s, HAC_LAG) || \
			 !strcmp(s, HAC_KERNEL) || \
                         !strcmp(s, HC_VERSION) || \
			 !strcmp(s, VECM_NORM) || \
			 !strcmp(s, GRETL_OPTIM) || \
			 !strcmp(s, NORMAL_RAND) || \
			 !strcmp(s, OPTIM_STEPLEN))

const char *intvar_code_string (const char *s)
{
    if (!strcmp(s, HAC_LAG)) {
	return hac_lag_string(); /* special */
    } else {
	return libset_option_string(s);
    }
}

static void libset_print_int (const char *s, PRN *prn,
			      gretlopt opt)
{
    if (coded_intvar(s)) {
	if (opt & OPT_D) {
	    pprintf(prn, " %s = %s\n", s, intvar_code_string(s));
	} else {
	    pprintf(prn, "set %s %s\n", s, intvar_code_string(s));
	}
    } else {
	int k = libset_get_int(s);

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

static void libset_print_double (const char *s, PRN *prn,
				 gretlopt opt)
{
    double x = libset_get_double(s);

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

/* print_settings: use OPT_D for "display", otherwise
   this gives script-type output */

static int print_settings (PRN *prn, gretlopt opt)
{
    if (opt & OPT_D) {
	pputs(prn, _("Variables that can be set using \"set\""));
	pputs(prn, " (");
	pputs(prn, _("\"help set\" for details"));
	pputs(prn, "):\n");
    }

    libset_header(N_("Program interaction and behavior"), prn, opt);

    if (opt & OPT_D) {
	pprintf(prn, " csv_delim = %s\n", arg_from_delim(data_delim));
	pprintf(prn, " csv_write_na = %s\n", get_csv_na_write_string());
	pprintf(prn, " csv_read_na = %s\n", get_csv_na_read_string());
    } else {
	const char *dl = arg_from_delim(data_delim);

	if (strcmp(dl, "unset")) {
	    pprintf(prn, "set csv_delim %s\n", arg_from_delim(data_delim));
	}
	pprintf(prn, "set csv_write_na %s\n", get_csv_na_write_string());
	pprintf(prn, "set csv_read_na %s\n", get_csv_na_read_string());
    }

    libset_print_int(CSV_DIGITS, prn, opt);	    
    libset_print_bool(ECHO, prn, opt);
    libset_print_bool(FORCE_DECP, prn, opt);
    libset_print_bool(HALT_ON_ERR, prn, opt);
    libset_print_int(LOOP_MAXITER, prn, opt);
    libset_print_bool(MAX_VERBOSE, prn, opt);
    libset_print_int(BFGS_VERBSKIP, prn, opt);
    libset_print_double(CONV_HUGE, prn, opt);
    libset_print_bool(MESSAGES, prn, opt);
    libset_print_bool(WARNINGS, prn, opt);
    libset_print_int(GRETL_DEBUG, prn, opt);
    libset_print_int(BLAS_MNK_MIN, prn, opt);
    libset_print_int(OMP_MNK_MIN, prn, opt);
    libset_print_int(OMP_N_THREADS, prn, opt);
    libset_print_int(SIMD_K_MAX, prn, opt);
    libset_print_int(SIMD_MN_MIN, prn, opt);

    if (opt & OPT_D) {
	libset_print_bool(SHELL_OK, prn, opt);
	if (*state->shelldir) {
	    pprintf(prn, " shelldir = '%s'\n", state->shelldir);
	} else {
	    pputs(prn, " shelldir = unset\n");
	}
    } 

    libset_print_bool(USE_CWD, prn, opt);
    libset_print_bool(SKIP_MISSING, prn, opt);

    libset_print_bool(R_LIB, prn, opt);
    libset_print_bool(R_FUNCTIONS, prn, opt);

    libset_header(N_("Numerical methods"), prn, opt);

    libset_print_int(GRETL_OPTIM, prn, opt);
    libset_print_int(BFGS_MAXITER, prn, opt);
    libset_print_double(BFGS_TOLER, prn, opt);
    libset_print_double(BFGS_MAXGRAD, prn, opt);
    libset_print_int(OPTIM_STEPLEN, prn, opt);
    libset_print_int(BHHH_MAXITER, prn, opt);
    libset_print_double(BHHH_TOLER, prn, opt);
    libset_print_int(RQ_MAXITER, prn, opt);
    libset_print_int(GMM_MAXITER, prn, opt);
    print_initvals(state->initvals, prn, opt);
    libset_print_bool(BFGS_RSTEP, prn, opt);
    libset_print_bool(USE_LBFGS, prn, opt);
    libset_print_int(LBFGS_MEM, prn, opt);
    libset_print_double(NLS_TOLER, prn, opt);
    libset_print_bool(USE_SVD, prn, opt);
    libset_print_bool(USE_FCP, prn, opt);
    libset_print_bool(DPDSTYLE, prn, opt);
    libset_print_double(NADARWAT_TRIM, prn, opt);
    libset_print_int(FDJAC_QUAL, prn, opt);

    libset_header(N_("Random number generation"), prn, opt);

    if (opt & OPT_D) {
	pprintf(prn, " seed = %u\n", gretl_rand_get_seed());
	pprintf(prn, " normal_rand = %s\n", libset_option_string(NORMAL_RAND));
    } else {
	pprintf(prn, "set seed %u\n", gretl_rand_get_seed());
	pprintf(prn, "set normal_rand %s\n", libset_option_string(NORMAL_RAND));
    }
    if (gretl_mpi_initialized()) {
	libset_print_bool(USE_DCMT, prn, opt);
    }

    libset_header(N_("Robust estimation"), prn, opt);

    libset_print_int(BOOTREP, prn, opt);
    libset_print_int(GARCH_VCV, prn, opt);
    libset_print_int(ARMA_VCV, prn, opt);
    libset_print_bool(FORCE_HC, prn, opt);
    libset_print_int(HAC_LAG, prn, opt);
    libset_print_int(HAC_KERNEL, prn, opt);
    libset_print_bool(PREWHITEN, prn, opt);
    libset_print_int(HC_VERSION, prn, opt);
    libset_print_bool(PCSE, prn, opt);
    libset_print_double(QS_BANDWIDTH, prn, opt);

    libset_header(N_("Time series"), prn, opt);

    libset_print_int(HORIZON, prn, opt);
    libset_print_int(VECM_NORM, prn, opt);

    pputc(prn, '\n');
    
    return 0;
}

static int libset_query_settings (const char *s, PRN *prn)
{
    int err = 0;

    if (libset_boolvar(s)) {
	pprintf(prn, "%s: boolean (on/off), currently %s\n", 
		s, libset_get_bool(s)? "on" : "off");
    } else if (coded_intvar(s)) {
	pprintf(prn, "%s: code, currently \"%s\"\n", s, intvar_code_string(s));
	coded_var_show_opts(s, prn);
    } else if (libset_int(s)) {
	int k = libset_get_int(s);
	
	if (is_unset(k)) {
	    pprintf(prn, "%s: positive integer, currently unset\n", s);
	} else {
	    pprintf(prn, "%s: positive integer, currently %d\n", s, k);
	}	    
    } else if (libset_double(s)) {
	double x = libset_get_double(s);

	if (na(x)) {
	    pprintf(prn, "%s: positive floating-point value, "
		    "currently automatic\n", s);
	} else {
	    pprintf(prn, "%s: positive floating-point value, "
		    "currently %g\n", s, x);
	}
    } else if (!strcmp(s, "initvals")) {
	if (state->initvals != NULL) {
	    pprintf(prn, "%s: matrix, currently\n", s);
	    gretl_matrix_print_to_prn(state->initvals, NULL, prn);
	} else {
	    pprintf(prn, "%s: matrix, currently null\n", s);
	}
    } else if (!strcmp(s, "matrix_mask")) {
	if (state->matmask != NULL) {
	    pprintf(prn, "%s: matrix, currently\n", s);
	    gretl_matrix_print_to_prn(state->matmask, NULL, prn);
	} else {
	    pprintf(prn, "%s: matrix, currently null\n", s);
	}
    } else if (!strcmp(s, "seed")) {
	pprintf(prn, "%s: unsigned int, currently %u\n",
		s, state->seed ? state->seed : gretl_rand_get_seed());
    } else if (!strcmp(s, "csv_delim")) {
	pprintf(prn, "%s: named character, currently \"%s\"\n", s,
		arg_from_delim(data_delim));
	coded_var_show_opts(s, prn);
    } else if (!strcmp(s, "shelldir")) {
	pprintf(prn, "%s: string, currently \"%s\"\n", s,
		state->shelldir);
    } else if (!strcmp(s, "workdir")) {
	pprintf(prn, "%s: string, currently \"%s\"\n", s,
		gretl_workdir());
    } else if (!strcmp(s, "csv_write_na")) {
	pprintf(prn, "%s: string, currently \"%s\"\n", s,
		state->csv_write_na);
    } else if (!strcmp(s, "csv_read_na")) {
	pprintf(prn, "%s: string, currently \"%s\"\n", s,
		state->csv_read_na);
    } else if (!strcmp(s, "stopwatch")) {
	err = 0;
    } else {
	err = 1;
    }

    return err;
}

int is_libset_var (const char *s)
{
    int err = libset_query_settings(s, NULL);

    return (err == 0);
}

#define default_ok(s) (!strcmp(s, BFGS_TOLER) || \
                       !strcmp(s, BHHH_TOLER))

#define default_str(s) (!strcmp(s, "auto") || !strcmp(s, "default"))

#define boolean_on(s) (!strcmp(s, "on") || !strcmp(s, "1") || \
                       !strcmp(s, "true"))

#define boolean_off(s) (!strcmp(s, "off") || !strcmp(s, "0") || \
                        !strcmp(s, "false"))

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

static int check_set_bool (const char *setobj, const char *setarg)
{
    if (boolean_on(setarg)) {
	return libset_set_bool(setobj, 1);
    } else if (boolean_off(setarg)) {
	return libset_set_bool(setobj, 0);
    } else {
	gretl_errmsg_sprintf(_("%s: invalid value '%s'"), setobj, setarg);
	return E_PARSE;
    }
}

int execute_set_line (const char *line, DATASET *dset, 
		      gretlopt opt, PRN *prn)
{
    char setobj[32], setarg[32], junk[8];
    int k, argc, err = E_PARSE;
    double x;

    check_for_state();

    if (opt != OPT_NONE) {
	return write_or_read_settings(opt, prn);
    }

    *setobj = *setarg = *junk = '\0';

    argc = sscanf(line, "%*s %31s %31s %7s", setobj, setarg, junk);

    if (argc <= 0) {
	return print_settings(prn, OPT_D);
    }

    if (argc > 1) {
	/* specials which need the whole line (FIXME) */
	if (!strcmp(setobj, "shelldir")) {
	    return set_shelldir(line);
	} else if (!strcmp(setobj, "workdir")) {
	    return set_workdir(line);
	}
    }

    if (argc == 3) {
	/* got some extraneous stuff */
	return err;
    }

    if (argc == 1) {
	if (!strcmp(setobj, ECHO)) {
	    state->flags |= STATE_ECHO_ON;
	    err = 0;
	} else if (!strcmp(setobj, "stopwatch")) {
	    gretl_stopwatch();
	    err = 0;
	} else {
	    return libset_query_settings(setobj, prn);
	}
    } else if (argc == 2) {
	if (!strcmp(setobj, "csv_write_na") || !strcmp(setobj, "csv_na")) {
	    return set_csv_na_write_string(setarg);
	} else if (!strcmp(setobj, "csv_read_na")) {
	    return set_csv_na_read_string(setarg);
	} else if (!strcmp(setobj, CSV_DIGITS)) {
	    return set_csv_digits(setarg);
	} else if (!strcmp(setobj, "initvals")) {
	    return set_initvals(setarg, prn);
	} else if (!strcmp(setobj, "matrix_mask")) {
	    return set_matmask(setarg, dset, prn);
	}

	if (libset_boolvar(setobj)) {
	    if (!strcmp(setobj, SHELL_OK)) {
		pprintf(prn, "You can only set this variable "
			"via the gretl GUI\n");
	    } else if (!strcmp(setobj, USE_OPENMP)) {
#if defined(_OPENMP)
		err = check_set_bool(setobj, setarg);
#else
		pprintf(prn, "Warning: openmp not supported\n");
#endif
	    } else {
		err = check_set_bool(setobj, setarg);
		if (!err && !strcmp(setobj, HALT_ON_ERR) &&
		    boolean_off(setarg)) {
		    pputs(prn, "Warning: \"set halt_on_error off\" is "
			  "deprecated and will be removed.\nPlease use "
			  "\"catch\" to trap errors instead.\n");
		}
	    }
	} else if (libset_double(setobj)) {
	    if (default_ok(setobj) && default_str(setarg)) {
		libset_set_double(setobj, NADBL);
		err = 0;
	    } else {
		err = libset_get_scalar(NULL, setarg, NULL, &x);
		if (!err) {
		    err = libset_set_double(setobj, x);
		}
	    }
	} else if (!strcmp(setobj, "csv_delim")) {
	    char c = delim_from_arg(setarg);

	    if (c > 0) {
		data_delim = c;
		err = 0;
	    }
	} else if (!strcmp(setobj, "seed")) {
	    err = libset_get_scalar(NULL, setarg, &k, NULL);
	    if (!err) {
		gretl_rand_set_seed((unsigned int) k);
		if (gretl_messages_on() && !gretl_looping_quietly()) {
		    pprintf(prn, 
			    _("Pseudo-random number generator seeded with %d\n"), k);
		}
		state->seed = k;
	    }
	} else if (!strcmp(setobj, HORIZON)) {
	    /* horizon for VAR impulse responses */
	    if (!strcmp(setarg, "auto")) {
		state->horizon = UNSET_INT;
		err = 0;
	    } else {
		err = libset_get_scalar(NULL, setarg, &k, NULL);
		if (!err) {
		    state->horizon = k;
		} else {
		    state->horizon = UNSET_INT;
		}
	    }
	} else if (coded_intvar(setobj)) {
	    err = parse_libset_int_code(setobj, setarg);
	} else if (libset_int(setobj)) {
	    err = libset_get_scalar(setobj, setarg, &k, NULL);
	    if (!err) {
		err = libset_set_int(setobj, k);
	    }
	} else {
	    gretl_errmsg_sprintf(_("set: unknown variable '%s'"), setobj);
	    err = E_UNKVAR;
	}
    }
		    
    return err;
}

double libset_get_double (const char *key)
{
    if (check_for_state()) {
	return NADBL;
    }

    if (!strcmp(key, QS_BANDWIDTH)) {
	if (!na(state->ropts.qsband) && state->ropts.qsband > 0) {
	    return state->ropts.qsband;
	} else {
	    /* what's a sensible default here? */
	    return 2.0;
	}
    } else if (!strcmp(key, NLS_TOLER)) {
	if (na(state->nls_toler)) {
	    return get_default_nls_toler();
	} else {
	    return state->nls_toler;
	}
    } else if (!strcmp(key, BHHH_TOLER)) {
	if (na(state->bhhh_toler)) {
	    return 1.0e-6;
	} else {
	    return state->bhhh_toler;
	}
    } else if (!strcmp(key, BFGS_TOLER)) {
	if (na(state->bfgs_toler)) {
	    return get_default_nls_toler();
	} else {
	    return state->bfgs_toler;
	}
    } else if (!strcmp(key, BFGS_MAXGRAD)) {
	return state->bfgs_maxgrad;
    } else if (!strcmp(key, NADARWAT_TRIM)) {
	if (na(state->nadarwat_trim)) {
	    return 4.0;
	} else {
	    return state->nadarwat_trim;
	}
    } else if (!strcmp(key, CONV_HUGE)) {
	if (na(state->conv_huge)) {
	    return 1.0e100;
	} else {
	    return state->conv_huge;
	}
    } else {
	fprintf(stderr, "libset_get_double: unrecognized "
		"variable '%s'\n", key);	
	return 0;
    }
}

double libset_get_user_tolerance (const char *key)
{
    if (!strcmp(key, NLS_TOLER)) {
	return state->nls_toler;
    } else if (!strcmp(key, BHHH_TOLER)) {
	return state->bhhh_toler;
    } else if (!strcmp(key, BFGS_TOLER)) {
	return state->bfgs_toler;
    } else if (!strcmp(key, BFGS_MAXGRAD)) {
	return state->bfgs_maxgrad;
    } else {
	return NADBL;
    }
}

int libset_set_double (const char *key, double val)
{
    int err = 0;

    if (check_for_state()) {
	return 1;
    }

    /* all the libset double vals must be positive */
    if (val <= 0.0) {
	return E_DATA;
    }

    if (!strcmp(key, QS_BANDWIDTH)) {
	state->ropts.qsband = val;
    } else if (!strcmp(key, NLS_TOLER)) {
	state->nls_toler = val;
    } else if (!strcmp(key, BHHH_TOLER)) {
	state->bhhh_toler = val;
    } else if (!strcmp(key, BFGS_TOLER)) {
	state->bfgs_toler = val;
    } else if (!strcmp(key, BFGS_MAXGRAD)) {
	state->bfgs_maxgrad = val;
    } else if (!strcmp(key, NADARWAT_TRIM)) {
	state->nadarwat_trim = val;
    } else if (!strcmp(key, CONV_HUGE)) {
	state->conv_huge = val;
    } else {
	fprintf(stderr, "libset_set_double: unrecognized "
		"variable '%s'\n", key);	
	err = E_UNKVAR;
    }

    return err;
}

int libset_get_int (const char *key)
{
    if (check_for_state()) {
	return 0;
    }

    if (!strcmp(key, BFGS_MAXITER)) {
	return state->bfgs_maxiter;
    } else if (!strcmp(key, OPTIM_STEPLEN)) {
	return state->optim_steplen;
    } else if (!strcmp(key, BHHH_MAXITER)) {
	return state->bhhh_maxiter;
    } else if (!strcmp(key, RQ_MAXITER)) {
	return state->rq_maxiter;
    } else if (!strcmp(key, GMM_MAXITER)) {
	return state->gmm_maxiter;
    } else if (!strcmp(key, LBFGS_MEM)) {
	return state->lbfgs_mem;
    } else if (!strcmp(key, BOOTREP)) {
	return state->bootrep;
    } else if (!strcmp(key, GARCH_VCV)) {
	return state->garch_vcv;
    } else if (!strcmp(key, GARCH_ROBUST_VCV)) {
	return state->garch_robust_vcv;
    } else if (!strcmp(key, ARMA_VCV)) {
	return state->arma_vcv;
    } else if (!strcmp(key, HAC_KERNEL)) {
	return state->ropts.hkern;
    } else if (!strcmp(key, HC_VERSION)) {
	return state->ropts.hc_version;
    } else if (!strcmp(key, HORIZON)) {
	return state->horizon;
    } else if (!strcmp(key, LOOP_MAXITER)) {
	return state->loop_maxiter;
    } else if (!strcmp(key, VECM_NORM)) {
	return state->vecm_norm;
    } else if (!strcmp(key, GRETL_OPTIM)) {
	return state->optim;
    } else if (!strcmp(key, GRETL_DEBUG)) {
	return gretl_debug;
    } else if (!strcmp(key, BLAS_MNK_MIN)) {
	return get_blas_mnk_min();
    } else if (!strcmp(key, OMP_MNK_MIN) || !strcmp(key, MP_MNK_MIN)) {
	return omp_mnk_min;
    } else if (!strcmp(key, OMP_N_THREADS)) {
	return omp_n_threads;
    } else if (!strcmp(key, SIMD_K_MAX)) {
	return get_simd_k_max();
    } else if (!strcmp(key, SIMD_MN_MIN)) {
	return get_simd_mn_min();
    } else if (!strcmp(key, BFGS_VERBSKIP)) {
	return state->bfgs_verbskip;
    } else if (!strcmp(key, CSV_DIGITS)) {
	return csv_digits;
    } else if (!strcmp(key, FDJAC_QUAL)) {
	return state->fdjac_qual;
    } else {
	fprintf(stderr, "libset_get_int: unrecognized "
		"variable '%s'\n", key);	
	return 0;
    }
}

static int intvar_min_max (const char *s, int *min, int *max,
			   int **var)
{
    *max = 100000;

    if (!strcmp(s, BFGS_MAXITER)) {
	*min = 0;
	*var = &state->bfgs_maxiter;
    } else if (!strcmp(s, BFGS_VERBSKIP)) {
	*min = 1;
	*var = &state->bfgs_verbskip;
    } else if (!strcmp(s, OPTIM_STEPLEN)) {
	*min = 0;
	*max = STEPLEN_MAX;
	*var = &state->optim_steplen;
    } else if (!strcmp(s, BHHH_MAXITER)) {
	*min = 1;
	*var = &state->bhhh_maxiter;
    } else if (!strcmp(s, RQ_MAXITER)) {
	*min = 1;
	*var = &state->rq_maxiter;
    } else if (!strcmp(s, GMM_MAXITER)) {
	*min = 1;
	*var = &state->gmm_maxiter;
    } else if (!strcmp(s, LBFGS_MEM)) {
	*min = 3;
	*max = 20;
	*var = &state->lbfgs_mem;
    } else if (!strcmp(s, BOOTREP)) {
	*min = 1;
	*var = &state->bootrep;
    } else if (!strcmp(s, HAC_KERNEL)) {
	*min = 0;
	*max = KERNEL_MAX;
    } else if (!strcmp(s, HC_VERSION)) {
	*min = 0;
	*max = 4 + 1;
	*var = &state->ropts.hc_version;
    } else if (!strcmp(s, HORIZON)) {
	*min = 1;
	*var = &state->horizon;
    } else if (!strcmp(s, LOOP_MAXITER)) {
	*min = 0;
	*var = &state->loop_maxiter;
    } else if (!strcmp(s, VECM_NORM)) {
	*min = 0;
	*max = NORM_MAX;
	*var = &state->vecm_norm;
    } else if (!strcmp(s, GRETL_OPTIM)) {
	*min = 0;
	*max = OPTIM_MAX;
	*var = &state->optim;
    } else if (!strcmp(s, GRETL_DEBUG)) {
	*min = 0;
	*var = &gretl_debug;
    } else if (!strcmp(s, FDJAC_QUAL)) {
	*min = 0;
	*max = 3;
	*var = &state->fdjac_qual;
    } else {
	fprintf(stderr, "libset_set_int: unrecognized "
		"variable '%s'\n", s);	
	return E_UNKVAR;
    }

    return 0;
}

int libset_set_int (const char *key, int val)
{
    if (check_for_state()) {
	return 1;
    }

    if (!strcmp(key, BLAS_MNK_MIN)) {
	set_blas_mnk_min(val);
	return 0;
    } else if (!strcmp(key, SIMD_K_MAX)) {
	set_simd_k_max(val);
	return 0;
    } else if (!strcmp(key, SIMD_MN_MIN)) {
	set_simd_mn_min(val);
	return 0;
    } else if (!strcmp(key, OMP_MNK_MIN) || !strcmp(key, MP_MNK_MIN)) {
	omp_mnk_min = val;
	return 0;
    } else if (!strcmp(key, OMP_N_THREADS)) {
	return set_omp_n_threads(val);
    } else {
	int min = 0, max = 0;
	int *ivar = NULL;
	int err = 0;

	err = intvar_min_max(key, &min, &max, &ivar);

	if (!err) {
	    if (val < min || val >= max || ivar == NULL) {
		err = E_DATA;
	    } else {
		*ivar = val;
	    }
	}

	return err;
    }
}

static int boolvar_get_flag (const char *s)
{
    if (!strcmp(s, ECHO)) {
	return STATE_ECHO_ON;
    } else if (!strcmp(s, MESSAGES)) {
	return STATE_MSGS_ON;
    } else if (!strcmp(s, WARNINGS)) {
	return STATE_WARN_ON;
    } else if (!strcmp(s, USE_SVD)) {
	return STATE_USE_SVD;
    } else if (!strcmp(s, USE_LBFGS)) {
	return STATE_USE_LBFGS;
    } else if (!strcmp(s, FORCE_DECP)) {
	return STATE_FORCE_DECPOINT;
    } else if (!strcmp(s, USE_CWD)) {
	return STATE_USE_CWD;
    } else if (!strcmp(s, USE_FCP)) {
	return STATE_USE_FCP;
    } else if (!strcmp(s, HALT_ON_ERR)) {
	return STATE_HALT_ON_ERR;
    } else if (!strcmp(s, MAX_VERBOSE)) {
	return STATE_MAX_VERBOSE;
    } else if (!strcmp(s, SHELL_OK)) {
	return STATE_SHELL_OK;
    } else if (!strcmp(s, FORCE_HC)) {
	return STATE_FORCE_HC;
    } else if (!strcmp(s, PREWHITEN)) {
	return STATE_PREWHITEN;
    } else if (!strcmp(s, PCSE)) {
	return STATE_USE_PCSE;
    } else if (!strcmp(s, SKIP_MISSING)) {
	return STATE_SKIP_MISSING;
    } else if (!strcmp(s, BFGS_RSTEP)) {
	return STATE_BFGS_RSTEP;
    } else if (!strcmp(s, DPDSTYLE)) {
	return STATE_DPDSTYLE_ON;
    } else if (!strcmp(s, USE_OPENMP)) {
	return STATE_OPENMP_ON;
    } else {
	fprintf(stderr, "libset_get_bool: unrecognized "
		"variable '%s'\n", s);	
	return 0;
    }
}

static void set_flag_from_env (int flag, const char *s, int neg)
{
    char *e = getenv(s);
    int action = 0;

    if (e != NULL) {
	if (*e != '\0' && *e != '0') {
	    action = (neg)? -1 : 1;
	} else {
	    action = (neg)? 1 : -1;
	}
    }

    if (action > 0) {
	state->flags |= flag;
    } else if (action < 0) {
	state->flags &= ~flag;
    }
}

static void maybe_check_env (const char *s)
{
    if (!strcmp(s, USE_SVD)) {
	set_flag_from_env(STATE_USE_SVD, "GRETL_USE_SVD", 0);
    } else if (!strcmp(s, USE_LBFGS)) {
	set_flag_from_env(STATE_USE_LBFGS, "GRETL_USE_LBFGS", 0);
    } else if (!strcmp(s, HALT_ON_ERR)) {
	set_flag_from_env(STATE_HALT_ON_ERR, "GRETL_KEEP_GOING", 1);
    }
}

int libset_get_bool (const char *key)
{
    int flag, ret = 0;

    /* global specials */

    if (!strcmp(key, R_FUNCTIONS)) {
	return R_functions;
    } else if (!strcmp(key, R_LIB)) {
	return R_lib;
    } else if (!strcmp(key, USE_DCMT)) {
        return gretl_rand_get_dcmt();
    }

    if (!strcmp(key, MAX_VERBOSE) && gretl_debug > 1) {
	/* strong debugging turns on max_verbose */
	return 1;
    }

    if (check_for_state()) {
	return 0;
    }

    maybe_check_env(key);

    flag = boolvar_get_flag(key);
    if (flag == 0) {
	fprintf(stderr, "libset_get_bool: unrecognized "
		"variable '%s'\n", key);
	ret = 0;
    } else {
	ret = flag_to_bool(state, flag);
    }

    return ret;
}

static void libset_set_decpoint (int on)
{
#ifdef ENABLE_NLS
    static char num_locale[32];

    if (on) {
	char *orig = setlocale(LC_NUMERIC, "");

	*num_locale = '\0';
	strncat(num_locale, orig, 31);
	setlocale(LC_NUMERIC, "C");
    } else {
	setlocale(LC_NUMERIC, num_locale);
    }

    reset_local_decpoint();
#endif
}

void set_data_export_decimal_comma (int s)
{
    if (s) {
	data_export_decpoint = ',';
    } else {
	data_export_decpoint = '.';
    }
}

char get_data_export_decpoint (void)
{
    char c = data_export_decpoint;

    /* revert to '.' on access */
    data_export_decpoint = '.';
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

static int check_R_setting (int *var, int val, const char *key)
{
    int err = 0;

#ifdef USE_RLIB
    *var = val;
#else
    if (val) {
	gretl_errmsg_sprintf("%s: not supported.", key);
	err = E_EXTERNAL;
    }
#endif

    return err;
}

int libset_set_bool (const char *key, int val)
{
    int flag, err = 0;

    if (check_for_state()) {
	return E_ALLOC;
    }

    /* global specials */

    if (!strcmp(key, R_FUNCTIONS)) {
	return check_R_setting(&R_functions, val, key);
    } else if (!strcmp(key, R_LIB)) {
	return check_R_setting(&R_lib, val, key);
    } else if (!strcmp(key, USE_DCMT)) {
	return gretl_rand_set_dcmt(val);
    }

    flag = boolvar_get_flag(key);

    if (flag == 0) {
	fprintf(stderr, "libset_set_bool: unrecognized "
		"variable '%s'\n", key);
	err = E_UNKVAR;
    } else if (val) {
	state->flags |= flag;
    } else {
	state->flags &= ~flag;
    }

    if (flag == STATE_FORCE_DECPOINT) {
	libset_set_decpoint(val);
    }

    return err;
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
	err = E_ALLOC;
    } else {
	sstack = realloc(state_stack, (ns + 1) * sizeof *sstack);
	if (sstack == NULL) {
	    free(newstate);
	    err = E_ALLOC;
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
    gretl_matrix_free(sv->initvals);
    gretl_matrix_free(sv->matmask);

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

/* switches for looping and batch mode: output: these depend on the
   state of the program calling libgretl, they are not user-settable
*/

void set_loop_on (int quiet, int progressive)
{
    state->flags |= STATE_LOOPING;
    if (quiet) {
	state->flags |= STATE_LOOP_QUIET;
    }
    if (progressive) {
	state->flags |= STATE_LOOP_PROG;
    }    
}

void set_loop_off (void)
{
    state->flags &= ~STATE_LOOPING;
    
    /* If we're not currently governed by "loop quietness" at
       caller level, turn such quietness off too 
    */
    if (state->flags & STATE_LOOP_QUIET) {
	int i = n_states - 1;

	if (i <= 0 || !(state_stack[i-1]->flags & STATE_LOOP_QUIET)) {
	    state->flags ^= STATE_LOOP_QUIET;
	}
    }
    
    /* and similarly for progressiveness */
    if (state->flags & STATE_LOOP_PROG) {
	int i = n_states - 1;

	if (i <= 0 || !(state_stack[i-1]->flags & STATE_LOOP_PROG)) {
	    state->flags ^= STATE_LOOP_PROG;
	}
    }    
}

/* returns 1 if there's a loop going on anywhere in the "caller
   ancestry" of the current execution level, else 0.
*/

int gretl_looping (void)
{
    int i, ns = n_states;

    for (i=0; i<ns; i++) {
	if (state_stack[i]->flags & STATE_LOOPING) {
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
    return (state->flags & STATE_LOOPING)? 1 : 0;
}

int gretl_looping_quietly (void)
{
    return (state->flags & STATE_LOOP_QUIET)? 1 : 0;
}

int gretl_looping_progressive (void)
{
    return (state->flags & STATE_LOOP_PROG)? 1 : 0;
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

/* mechanism to support callback for printing iteration info */

static ITER_PRINT_FUNC ifunc;

void set_iter_print_func (ITER_PRINT_FUNC func)
{
    ifunc = func;
}

int iter_print_func_installed (void)
{
    return ifunc != NULL;
}

int iter_print_callback (int i, PRN *prn)
{
    int ret = 0;

    if (ifunc != NULL) {
	ret = (*ifunc)(i, prn);
    }

    return ret;
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

/* mechanism for interactive debugging */

static DEBUG_READLINE dbg_readline;

void set_debug_read_func (DEBUG_READLINE dfunc) 
{
    dbg_readline = dfunc;
}

DEBUG_READLINE get_debug_read_func (void) 
{
    return dbg_readline;
}

static DEBUG_OUTPUT dbg_output;

void set_debug_output_func (DEBUG_OUTPUT dout) 
{
    dbg_output = dout;
}

DEBUG_OUTPUT get_debug_output_func (void) 
{
    return dbg_output;
}

/* support for GUI Stop button */

static QUERY_STOP query_stop;

void set_query_stop_func (QUERY_STOP query) 
{
    query_stop = query;
}

int check_for_stop (void)
{
    if (query_stop != NULL) {
	return (*query_stop)();
    } else {
	return 0;
    }
}

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

static int set_csv_digits (const char *s)
{
    int k = atoi(s);

    if (k > 0 && k < 26) {
	csv_digits = k;
	return 0;
    } else {
	return E_DATA;
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
	char line[1024];

	while (fgets(line, sizeof line, fp)) {
	    if (*line == '#' || string_is_blank(line)) {
		continue;
	    }
	    tailstrip(line);
	    err = execute_set_line(line, NULL, OPT_NONE, prn);
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
