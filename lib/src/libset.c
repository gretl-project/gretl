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

#ifdef _OPENMP
# include <omp.h>
#endif

#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

#ifdef WIN32
# include "gretl_win32.h"
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

enum {
    INIT_VALS,
    INIT_CURV
};

/* state flags */

enum {
    STATE_USE_CWD         = 1 << 0,  /* store: use current dir as default */
    STATE_ECHO_ON         = 1 << 1,  /* echoing commands or not */
    STATE_MSGS_ON         = 1 << 2,  /* emitting non-error messages or not */
    STATE_FORCE_DECPOINT  = 1 << 3,  /* override locale decimal separator */
    STATE_USE_PCSE        = 1 << 4,  /* Beck-Katz panel-corrected std errs */
    STATE_USE_SVD         = 1 << 5,  /* SVD decomposition is matrix OLS default */
    STATE_USE_QR          = 1 << 6,  /* QR decomp is least-squares command default */
    STATE_PREWHITEN       = 1 << 7,  /* HAC pre-whitening? */
    STATE_FORCE_HC        = 1 << 8,  /* don't use HAC for time series */
    STATE_USE_LBFGS       = 1 << 9,  /* prefer LBFGS to BFGS? */
    STATE_SHELL_OK        = 1 << 10, /* "shell" facility is approved? */
    STATE_WARN_ON         = 1 << 11, /* print numerical warning messages */
    STATE_SKIP_MISSING    = 1 << 12, /* skip NAs when building matrix from series */
    STATE_BFGS_RSTEP      = 1 << 13, /* use Richardson in BFGS numerical gradient */
    STATE_DPDSTYLE_ON     = 1 << 14, /* emulate dpd in dynamic panel data models */
    STATE_OPENMP_ON       = 1 << 15, /* using openmp */
    STATE_ROBUST_Z        = 1 << 16, /* use z- not t-score with HCCM/HAC */
    STATE_MWRITE_G        = 1 << 17, /* use %g format with mwrite() */
    STATE_ECHO_SPACE      = 1 << 18, /* preserve vertical space in output */
    STATE_MPI_SMT         = 1 << 19  /* MPI: use hyperthreads by default */
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
    int max_verbose;            /* optimizer verbosity level */
    int boot_iters;             /* max iterations, IRF bootstrap */
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
    gretl_matrix *initcurv;     /* initial curvature matrix for BFGS */
    gretl_matrix *matmask;      /* mask for series -> matrix conversion */
    struct robust_opts ropts;   /* robust standard error options */
    char csv_write_na[8];       /* representation of NA in CSV output */
    char csv_read_na[8];        /* representation of NA in CSV input */
    double nadarwat_trim;       /* multiple of h to use in nadarwat() for trimming */
    int fdjac_qual;             /* quality of "fdjac" function */
    double fdjac_eps;           /* finite increment for "fdjac" function */
    int wildboot_dist;          /* distribution for wild bootstrap */
};

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

#define libset_boolvar(s) (!strcmp(s, MESSAGES) || \
                           !strcmp(s, WARNINGS) || \
                           !strcmp(s, FORCE_DECP) || \
			   !strcmp(s, FORCE_HC) || \
			   !strcmp(s, USE_LBFGS) || \
			   !strcmp(s, PCSE) || \
			   !strcmp(s, PREWHITEN) || \
			   !strcmp(s, USE_SVD) || \
			   !strcmp(s, USE_QR) || \
			   !strcmp(s, SHELL_OK) || \
			   !strcmp(s, USE_CWD) || \
                           !strcmp(s, SKIP_MISSING) || \
			   !strcmp(s, R_FUNCTIONS) || \
			   !strcmp(s, R_LIB) || \
			   !strcmp(s, BFGS_RSTEP) || \
			   !strcmp(s, DPDSTYLE) || \
			   !strcmp(s, USE_DCMT) || \
			   !strcmp(s, ROBUST_Z) || \
			   !strcmp(s, MWRITE_G) || \
			   !strcmp(s, MPI_USE_SMT) || \
			   !strcmp(s, USE_OPENMP))

#define libset_double(s) (!strcmp(s, CONV_HUGE) || \
			  !strcmp(s, BFGS_TOLER) || \
			  !strcmp(s, BFGS_MAXGRAD) || \
			  !strcmp(s, BHHH_TOLER) || \
			  !strcmp(s, NLS_TOLER) || \
			  !strcmp(s, QS_BANDWIDTH) || \
			  !strcmp(s, NADARWAT_TRIM) || \
			  !strcmp(s, FDJAC_EPS))

#define libset_int(s) (!strcmp(s, BFGS_MAXITER) || \
		       !strcmp(s, MAX_VERBOSE) || \
		       !strcmp(s, BOOT_ITERS) || \
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
		       !strcmp(s, GRETL_ASSERT) || \
		       !strcmp(s, BLAS_MNK_MIN) || \
		       !strcmp(s, OMP_MNK_MIN) || \
		       !strcmp(s, MP_MNK_MIN) || \
		       !strcmp(s, OMP_N_THREADS) || \
		       !strcmp(s, SIMD_K_MAX) || \
		       !strcmp(s, SIMD_MN_MIN) || \
		       !strcmp(s, FDJAC_QUAL) || \
		       !strcmp(s, WILDBOOT_DIST) || \
		       !strcmp(s, QUANTILE_TYPE) || \
		       !strcmp(s, PLOT_COLLECTION))

/* global state */
set_vars *state;
static int seed_is_set;
static int gretl_debug;
static int user_mp_bits;
static int R_functions;
static int R_lib = 1;
static int csv_digits = UNSET_INT;
static int comments_on = 0;
static int gretl_assert = 0;
static int plot_collection = 0;
static int Qtype = 0;
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
static int libset_get_scalar (const char *var, const char *arg,
			      int *pi, double *px);

static void robust_opts_init (struct robust_opts *r)
{
    r->auto_lag = AUTO_LAG_STOCK_WATSON;
    r->user_lag = UNSET_INT;
    r->hc_version = 0;
    r->hkern = KERNEL_BARTLETT;
    r->qsband = NADBL;
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

static const char *maxverb_strs[] = {
    "off",
    "on",
    "full",
    NULL
};

static const char *steplen_strs[] = {
    "power",
    "quadratic",
    NULL
};

static const char *wildboot_strs[] = {
    "rademacher",
    "mammen",
    NULL
};

static const char *qtype_strs[] = {
    "Q6",
    "Q7",
    "Q8",
    NULL
};

static const char *assert_strs[] = {
    "off",
    "warn",
    "stop",
    NULL
};

static const char *plotcoll_strs[] = {
    "off",
    "auto",
    "on",
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
    } else if (!strcmp(s, MAX_VERBOSE)) {
	return maxverb_strs;
    } else if (!strcmp(s, "csv_delim")) {
	return csv_delim_args;
    } else if (!strcmp(s, OPTIM_STEPLEN)) {
	return steplen_strs;
    } else if (!strcmp(s, WILDBOOT_DIST)) {
	return wildboot_strs;
    } else if (!strcmp(s, QUANTILE_TYPE)) {
	return qtype_strs;
    } else if (!strcmp(s, GRETL_ASSERT)) {
	return assert_strs;
    } else if (!strcmp(s, PLOT_COLLECTION)) {
	return plotcoll_strs;
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
	if (!strcmp(s, "csv_delim")) {
	    pputs(prn, " or quoted punctuation character");
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
    } else if (!strcmp(s, OPTIM_STEPLEN)) {
	return steplen_strs[state->optim_steplen];
    } else if (!strcmp(s, MAX_VERBOSE)) {
	return maxverb_strs[state->max_verbose];
    } else if (!strcmp(s, WILDBOOT_DIST)) {
	return wildboot_strs[state->wildboot_dist];
    } else if (!strcmp(s, QUANTILE_TYPE)) {
	return qtype_strs[Qtype];
    } else if (!strcmp(s, GRETL_ASSERT)) {
	return assert_strs[gretl_assert];
    } else if (!strcmp(s, PLOT_COLLECTION)) {
	return plotcoll_strs[plot_collection];
    } else {
	return "?";
    }
}

static void print_initmat (const gretl_matrix *imat,
			   int type, PRN *prn,
			   gretlopt opt)
{
    char *name;

    if (type == INIT_VALS) {
	name = "initvals";
    } else if (type == INIT_CURV) {
	name = "initcurv";
    }

    if (opt & OPT_D) {
	if (imat == NULL) {
	    pprintf(prn, " %s = auto\n", name);
	} else {
	    gretl_matrix_print_to_prn(imat, name, prn);
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
    /* copy everything */
    *sv = *state;
    /* but set matrix pointers to NULL */
    sv->initvals = NULL;
    sv->initcurv = NULL;
    sv->matmask = NULL;
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
    int mib[2] = {CTL_HW, HW_NCPU};
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

int gretl_n_physical_cores (void)
{
    static int n_cores = -1;

    if (n_cores >= 1) {
	return n_cores;
    }

    /* this may well not be what we want, but it's a
       starting point and fallback
    */
    n_cores = gretl_n_processors();

#if defined(WIN32)
    int nc = win32_get_core_count();

    if (nc > 0) {
	n_cores = nc;
    }
#elif defined(OS_OSX)
    if (n_cores > 1) {
	int nc = 0;
	size_t len = sizeof nc;

	if (sysctlbyname("hw.physicalcpu", &nc, &len, NULL, 0) == -1) {
	    perror("could not determine number of physical cores available");
	} else {
	    n_cores = nc;
	}
    }
#else
    if (n_cores > 1) {
	/* check SMT status */
	FILE *fp = fopen("/sys/devices/system/cpu/smt/active", "r");
	char line[2];
	int smt = 0;

	if (fp != NULL) {
	    if (fgets(line, sizeof line, fp)) {
		smt = atoi(line);
	    }
	    fclose(fp);
	}
	if (smt) {
	    n_cores /= 2;
	}
    }
#endif

    return n_cores;
}

#define OMP_SHOW 0

int libset_use_openmp (guint64 n)
{
#if defined(_OPENMP)
    if (state == NULL || !(state->flags & STATE_OPENMP_ON) || omp_n_threads < 2) {
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
	gretl_errmsg_sprintf(_("omp_num_threads: must be >= 1 and <= %d"),
			     gretl_n_processors());
	return E_DATA;
    } else {
	omp_n_threads = n;
	omp_set_num_threads(n);
	if (blas_is_openblas()) {
	    blas_set_num_threads(n);
	}
    }
#else
    gretl_warnmsg_set(_("set_omp_n_threads: OpenMP is not enabled"));
#endif

    return 0;
}

void num_threads_init (int blas_type)
{
    int nc = gretl_n_physical_cores();

#if defined(_OPENMP)
    omp_n_threads = nc;
    omp_set_num_threads(nc);
#endif
    if (blas_type == BLAS_OPENBLAS) {
	blas_set_num_threads(nc);
    }
    if (blas_type > BLAS_NETLIB) {
	set_blas_mnk_min(90000);
    }
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

#define LOOP_MAXITER_DEFAULT 100000

static void state_vars_init (set_vars *sv)
{
#if PDEBUG
    fprintf(stderr, "state_vars_init called\n");
#endif
    sv->flags = STATE_ECHO_ON | STATE_MSGS_ON | STATE_WARN_ON |
	STATE_SKIP_MISSING;
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
    sv->loop_maxiter = LOOP_MAXITER_DEFAULT;
    sv->rq_maxiter = 1000;
    sv->gmm_maxiter = 250;
    sv->vecm_norm = NORM_PHILLIPS;
    sv->optim = OPTIM_AUTO;
    sv->max_verbose = 0;
    sv->initvals = NULL;
    sv->initcurv = NULL;
    sv->matmask = NULL;

    sv->bfgs_maxiter = UNSET_INT;
    sv->bfgs_toler = NADBL;
    sv->bfgs_maxgrad = 5.0;
    sv->bfgs_verbskip = 1;
    sv->optim_steplen = STEPLEN_POWER;
    sv->bhhh_maxiter = 500;
    sv->bhhh_toler = NADBL;
    sv->boot_iters = 1999;
    sv->lbfgs_mem = 8;
    sv->garch_vcv = ML_UNSET;
    sv->arma_vcv = ML_HESSIAN;
    sv->garch_robust_vcv = ML_UNSET;
    sv->nadarwat_trim = 4.0;
    sv->fdjac_qual = 0;
    sv->fdjac_eps = 0.0;
    sv->wildboot_dist = 0;

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

int gretl_comments_on (void)
{
    if (gretl_function_depth() > 0) {
	return 0;
    } else {
	return comments_on;
    }
}

int gretl_echo_space (void)
{
    if (check_for_state()) return 0;
    return flag_to_bool(state, STATE_ECHO_SPACE);
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
    if (check_for_state()) {
	return 1;
    } else {
	return flag_to_bool(state, STATE_MSGS_ON);
    }
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
	return state->initcurv->rows;
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
    } else {
	int k = 0;

	err = libset_get_scalar(HAC_LAG, s, &k, NULL);
	if (!err) {
	    state->ropts.user_lag = k;
	}
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
    int nstatus, err = 0;

    if (arg == NULL || *arg == '\0') {
	return E_ARGS;
    }

    nstatus = libset_numeric_test(arg, pi, px);

    if (nstatus == NUMERIC_BAD) {
	return E_INVARG; /* handled */
    } else if (nstatus == NUMERIC_OK) {
	if (pi != NULL && negval_invalid(var) && *pi < 0) {
	    err = E_INVARG;
	} else if (px != NULL && *px < 0.0) {
	    err = E_INVARG;
	}
	return err; /* handled */
    }

    /* handle the non-numeric case */
    x = get_scalar_value_by_name(arg, &err);

    if (!err) {
	if (negval_invalid(var) && x < 0.0) {
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

static int libset_get_unsigned (const char *arg, unsigned int *pu)
{
    unsigned long lu = 0;
    char *test = NULL;
    double x = NADBL;
    int err = 0;

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
    } else if (!g_ascii_strcasecmp(key, MAX_VERBOSE)) {
	for (i=0; maxverb_strs[i] != NULL; i++) {
	    if (!g_ascii_strcasecmp(val, maxverb_strs[i])) {
		state->max_verbose = i;
		err = 0;
		break;
	    }
	}
	if (err && (strcmp(val, "0") == 0 || strcmp(val, "1") == 0)) {
	    state->max_verbose = atoi(val);
	    err = 0;
	}
    } else if (!g_ascii_strcasecmp(key, WILDBOOT_DIST)) {
	for (i=0; wildboot_strs[i] != NULL; i++) {
	    if (!g_ascii_strcasecmp(val, wildboot_strs[i])) {
		state->wildboot_dist = i;
		err = 0;
		break;
	    }
	}
    } else if (!g_ascii_strcasecmp(key, QUANTILE_TYPE)) {
	for (i=0; qtype_strs[i] != NULL; i++) {
	    if (!g_ascii_strcasecmp(val, qtype_strs[i])) {
		Qtype = i;
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
    } else if (!g_ascii_strcasecmp(key, GRETL_ASSERT)) {
	for (i=0; assert_strs[i] != NULL; i++) {
	    if (!g_ascii_strcasecmp(val, assert_strs[i])) {
		gretl_assert = i;
		err = 0;
		break;
	    }
	}
    } else if (!g_ascii_strcasecmp(key, PLOT_COLLECTION)) {
	for (i=0; plotcoll_strs[i] != NULL; i++) {
	    if (!g_ascii_strcasecmp(val, plotcoll_strs[i])) {
		plot_collection = i;
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
	if (!strcmp(scpy, "qml")) {
	    state->garch_robust_vcv = ML_QML;
	} else if (!strcmp(scpy, "bw")) {
	    state->garch_robust_vcv = ML_BW;
	}
	free(scpy);
    }
}

static int set_initmat (const char *mname, int type,
			PRN *prn)
{
    gretl_matrix *m;
    int err = 0;

    if (!strcmp(mname, "auto")) {
	if (type == INIT_VALS) {
	    gretl_matrix_free(state->initvals);
	    state->initvals = NULL;
	} else if (type == INIT_CURV) {
	    gretl_matrix_free(state->initcurv);
	    state->initcurv = NULL;
	}
    } else {
	m = get_matrix_by_name(mname);
	if (m == NULL) {
	    pprintf(prn, _("'%s': no such matrix"), mname);
	    pputc(prn, '\n');
	    err = E_DATA;
	} else {
	    if (type == INIT_VALS) {
		state->initvals = gretl_matrix_copy(m);
		if (state->initvals == NULL) {
		    err = E_ALLOC;
		}
	    } else if (type == INIT_CURV) {
		state->initcurv = gretl_matrix_copy(m);
		if (state->initcurv == NULL) {
		    err = E_ALLOC;
		}
	    }
	}
    }

    return err;
}

static int set_echo_status (const char *arg)
{
    int err = 0;

    if (check_for_state()) {
	return 1;
    }

    if (!strcmp(arg, "on")) {
	state->flags |= STATE_ECHO_ON;
    } else if (!strcmp(arg, "off")) {
	state->flags &= ~STATE_ECHO_ON;
	state->flags &= ~STATE_ECHO_SPACE;
    } else if (!strcmp(arg, "space")) {
	state->flags |= STATE_ECHO_SPACE;
    } else if (!strcmp(arg, "full")) {
	state->flags |= STATE_ECHO_ON;
	state->flags |= STATE_ECHO_SPACE;
    } else {
	err = E_INVARG;
    }

    return err;
}

static const char *get_echo_status (void)
{
    if (check_for_state()) {
	return "on";
    } else if ((state->flags & STATE_ECHO_ON) &&
	       (state->flags & STATE_ECHO_SPACE)) {
	return "full";
    } else if (state->flags & STATE_ECHO_ON) {
	return "on";
    } else if (state->flags & STATE_ECHO_SPACE) {
	return "space";
    } else {
	return "off";
    }
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

const char *csv_delims = ", \t;";

static char delim_from_arg (const char *s)
{
    int i;

    for (i=0; csv_delim_args[i] != NULL; i++) {
	if (!strcmp(s, csv_delim_args[i])) {
	    return csv_delims[i];
	}
    }

    if (strlen(s) == 1 && ispunct(*s)) {
	return s[0];
    }

    return 0;
}

static const char *arg_from_delim (char c)
{
    static char d[2];
    int i;

    for (i=0; csv_delims[i] != '\0'; i++) {
	if (c == csv_delims[i]) {
	    return csv_delim_args[i];
	}
    }

    if (ispunct(c)) {
	sprintf(d, "%c", c);
	return d;
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
			 !strcmp(s, OPTIM_STEPLEN) || \
			 !strcmp(s, MAX_VERBOSE) || \
			 !strcmp(s, WILDBOOT_DIST) || \
			 !strcmp(s, QUANTILE_TYPE) || \
			 !strcmp(s, GRETL_ASSERT) || \
			 !strcmp(s, PLOT_COLLECTION))

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
	if (na(x) || (x == 0.0 && !strcmp(s, FDJAC_EPS))) {
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
    const char *workdir = gretl_workdir();

    if (opt & OPT_D) {
	pputs(prn, _("Variables that can be set using \"set\""));
	pputs(prn, " (");
	pputs(prn, _("\"help set\" for details"));
	pputs(prn, "):\n");
    }

    libset_header(N_("Program interaction and behavior"), prn, opt);

    if (opt & OPT_D) {
	pprintf(prn, " workdir = '%s'\n", workdir);
    } else if (0) {
	/* non-portable? */
	if (strchr(workdir, ' ')) {
	    pprintf(prn, "set workdir \"%s\"\n", workdir);
	} else {
	    pprintf(prn, "set workdir %s\n", workdir);
	}
    }

    if (opt & OPT_D) {
	pprintf(prn, " csv_delim = %s\n", arg_from_delim(data_delim));
	pprintf(prn, " csv_write_na = %s\n", get_csv_na_write_string());
	pprintf(prn, " csv_read_na = %s\n", get_csv_na_read_string());
	pprintf(prn, " display_digits = %d\n", get_gretl_digits());
    } else {
	const char *dl = arg_from_delim(data_delim);

	if (strcmp(dl, "unset")) {
	    pprintf(prn, "set csv_delim %s\n", arg_from_delim(data_delim));
	}
	pprintf(prn, "set csv_write_na %s\n", get_csv_na_write_string());
	pprintf(prn, "set csv_read_na %s\n", get_csv_na_read_string());
    }

    libset_print_int(CSV_DIGITS, prn, opt);

    if (opt & OPT_D) {
	pprintf(prn, " echo = %s\n", get_echo_status());
    } else {
	pprintf(prn, "set echo %s\n", get_echo_status());
    }

    libset_print_bool(FORCE_DECP, prn, opt);
    libset_print_int(LOOP_MAXITER, prn, opt);
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
	/* display only */
	libset_print_bool(SHELL_OK, prn, opt);
    }

    libset_print_bool(USE_CWD, prn, opt);
    libset_print_bool(SKIP_MISSING, prn, opt);

    libset_print_bool(R_LIB, prn, opt);
    libset_print_bool(R_FUNCTIONS, prn, opt);

    libset_header(N_("Numerical methods"), prn, opt);

    libset_print_int(GRETL_OPTIM, prn, opt);
    libset_print_int(MAX_VERBOSE, prn, opt);
    libset_print_int(BFGS_MAXITER, prn, opt);
    libset_print_double(BFGS_TOLER, prn, opt);
    libset_print_double(BFGS_MAXGRAD, prn, opt);
    libset_print_int(OPTIM_STEPLEN, prn, opt);
    libset_print_int(BHHH_MAXITER, prn, opt);
    libset_print_double(BHHH_TOLER, prn, opt);
    libset_print_int(RQ_MAXITER, prn, opt);
    libset_print_int(GMM_MAXITER, prn, opt);
    print_initmat(state->initvals, INIT_VALS, prn, opt);
    print_initmat(state->initcurv, INIT_CURV, prn, opt);
    libset_print_bool(BFGS_RSTEP, prn, opt);
    libset_print_bool(USE_LBFGS, prn, opt);
    libset_print_int(LBFGS_MEM, prn, opt);
    libset_print_double(NLS_TOLER, prn, opt);
    libset_print_bool(USE_SVD, prn, opt);
    libset_print_bool(USE_QR, prn, opt);
    libset_print_bool(DPDSTYLE, prn, opt);
    libset_print_double(NADARWAT_TRIM, prn, opt);
    libset_print_int(FDJAC_QUAL, prn, opt);
    libset_print_double(FDJAC_EPS, prn, opt);

    libset_header(N_("Random number generation"), prn, opt);

    if (opt & OPT_D) {
	pprintf(prn, " seed = %u\n", gretl_rand_get_seed());
    } else {
	if (seed_is_set) {
	    pprintf(prn, "set seed %u\n", gretl_rand_get_seed());
	}
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
    libset_print_int(BOOT_ITERS, prn, opt);

    pputc(prn, '\n');

    return 0;
}

static int libset_query_settings (const char *s, PRN *prn)
{
    int err = 0;

    if (!strcmp(s, "echo")) {
	pprintf(prn, "%s: code, currently '%s'\n", s, get_echo_status());
    } else if (libset_boolvar(s)) {
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
    } else if (!strcmp(s, "initcurv")) {
	if (state->initcurv != NULL) {
	    pprintf(prn, "%s: matrix, currently\n", s);
	    gretl_matrix_print_to_prn(state->initcurv, NULL, prn);
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
	pprintf(prn, "%s: unsigned int, currently %u (%s)\n",
		s, state->seed ? state->seed : gretl_rand_get_seed(),
		seed_is_set ? "set by user" : "automatic");
    } else if (!strcmp(s, "csv_delim")) {
	pprintf(prn, "%s: named character, currently \"%s\"\n", s,
		arg_from_delim(data_delim));
	coded_var_show_opts(s, prn);
    } else if (!strcmp(s, "workdir")) {
	pprintf(prn, "%s: string, currently \"%s\"\n", s,
		gretl_workdir());
    } else if (!strcmp(s, "csv_write_na")) {
	pprintf(prn, "%s: string, currently \"%s\"\n", s,
		state->csv_write_na);
    } else if (!strcmp(s, "csv_read_na")) {
	pprintf(prn, "%s: string, currently \"%s\"\n", s,
		state->csv_read_na);
    } else if (!strcmp(s, "display_digits")) {
	pprintf(prn, "%s: integer, currently %d\n", s,
		get_gretl_digits());
    } else if (!strcmp(s, "stopwatch")) {
	err = 0;
    } else if (!strcmp(s, "verbose")) {
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
                       !strcmp(s, BHHH_TOLER) || \
		       !strcmp(s, NLS_TOLER))

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

int execute_set (const char *setobj, const char *setarg,
		 DATASET *dset, gretlopt opt, PRN *prn)
{
    int k, argc, err;
    unsigned int u;

    check_for_state();

    if (opt != OPT_NONE) {
	return write_or_read_settings(opt, prn);
    }

    argc = (setobj != NULL) + (setarg != NULL);

    if (argc == 0) {
	return print_settings(prn, OPT_D);
    }

    /* set error default */
    err = E_PARSE;

    if (argc == 1) {
	if (!strcmp(setobj, "stopwatch")) {
	    gretl_stopwatch();
	    return 0;
	} else {
	    return libset_query_settings(setobj, prn);
	}
    } else if (argc == 2) {
	if (!strcmp(setobj, "shelldir")) {
	    pputs(prn, "'shelldir' is obsolete, please use 'workdir'\n");
	    return 0;
	} else if (!strcmp(setobj, "workdir")) {
	    return set_workdir(setarg);
	} else if (!strcmp(setobj, "csv_write_na") || !strcmp(setobj, "csv_na")) {
	    return set_csv_na_write_string(setarg);
	} else if (!strcmp(setobj, "csv_read_na")) {
	    return set_csv_na_read_string(setarg);
	} else if (!strcmp(setobj, CSV_DIGITS)) {
	    return set_csv_digits(setarg);
	} else if (!strcmp(setobj, "initvals")) {
	    return set_initmat(setarg, INIT_VALS, prn);
	} else if (!strcmp(setobj, "initcurv")) {
	    return set_initmat(setarg, INIT_CURV, prn);
	} else if (!strcmp(setobj, "matrix_mask")) {
	    return set_matmask(setarg, dset, prn);
	} else if (!strcmp(setobj, "graph_theme")) {
	    return set_plotstyle(setarg);
	} else if (!strcmp(setobj, "display_digits")) {
	    if (gretl_function_depth() > 0) {
		pprintf(prn, "'%s': cannot be set inside a function\n");
		return E_INVARG;
	    } else {
		return set_gretl_digits(atoi(setarg));
	    }
	} else if (!strcmp(setobj, "echo")) {
	    return set_echo_status(setarg);
	} else if (!strcmp(setobj, "verbose")) {
	    err = 0;
	    if (!strcmp(setarg, "on")) {
		set_gretl_messages(1);
		set_echo_status(setarg);
	    } else if (!strcmp(setarg, "off")) {
		set_gretl_messages(0);
		set_echo_status(setarg);
		comments_on = 0;
	    } else if (!strcmp(setarg, "comments")) {
		set_gretl_messages(0);
		set_echo_status("off");
		comments_on = 1;
	    } else {
		err = E_INVARG;
	    }
	    return err;
	}

	if (libset_boolvar(setobj)) {
	    if (!strcmp(setobj, SHELL_OK)) {
		pprintf(prn, "'%s': this must be set via the gretl GUI\n", setobj);
		err = E_DATA;
	    } else if (!strcmp(setobj, USE_OPENMP)) {
#if defined(_OPENMP)
		err = check_set_bool(setobj, setarg);
#else
		pprintf(prn, "Warning: openmp not supported\n");
#endif
	    } else {
		err = check_set_bool(setobj, setarg);
	    }
	} else if (libset_double(setobj)) {
	    if (default_ok(setobj) && default_str(setarg)) {
		libset_set_double(setobj, NADBL);
		err = 0;
	    } else {
		double x;

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
	    err = libset_get_unsigned(setarg, &u);
	    if (!err) {
		gretl_rand_set_seed(u);
		if (gretl_messages_on()) {
		    pprintf(prn,
			    _("Pseudo-random number generator seeded with %u\n"), u);
		}
		state->seed = u;
		seed_is_set = 1;
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
    } else if (!strcmp(key, FDJAC_EPS)) {
	if (na(state->fdjac_eps)) {
	    return 0.0;
	} else {
	    return state->fdjac_eps;
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

    /* all the libset double vals must be positive, except for
       FDJAC_EPS, where 0.0 means "auto"
    */
    if (val < 0.0) {
	return E_DATA;
    } else if (val == 0.0 && strcmp(key, FDJAC_EPS)) {
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
    } else if (!strcmp(key, FDJAC_EPS)) {
	state->fdjac_eps = val;
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
    } else if (!strcmp(key, MAX_VERBOSE)) {
	return state->max_verbose;
    } else if (!strcmp(key, BOOT_ITERS)) {
	return state->boot_iters;
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
    } else if (!strcmp(key, GRETL_ASSERT)) {
	return gretl_assert;
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
    } else if (!strcmp(key, WILDBOOT_DIST)) {
	return state->wildboot_dist;
    } else if (!strcmp(key, QUANTILE_TYPE)) {
	return Qtype;
    } else if (!strcmp(key, PLOT_COLLECTION)) {
	return plot_collection;
    } else if (!strcmp(key, "loop_maxiter_default")) {
	return LOOP_MAXITER_DEFAULT; /* for internal use */
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
    } else if (!strcmp(s, BOOT_ITERS)) {
	*max = 999999;
	*min = 499;
	*var = &state->boot_iters;
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
	*max = INT_MAX - 1;
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
	*max = 4;
	*var = &state->fdjac_qual;
    } else if (!strcmp(s, WILDBOOT_DIST)) {
	*min = 0;
	*max = 1;
	*var = &state->wildboot_dist;
    } else if (!strcmp(s, PLOT_COLLECTION)) {
	*min = 0;
	*max = 2;
	*var = &plot_collection;
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
    if (!strcmp(s, MESSAGES)) {
	return STATE_MSGS_ON;
    } else if (!strcmp(s, WARNINGS)) {
	return STATE_WARN_ON;
    } else if (!strcmp(s, USE_SVD)) {
	return STATE_USE_SVD;
    } else if (!strcmp(s, USE_QR)) {
        return STATE_USE_QR;
    } else if (!strcmp(s, USE_LBFGS)) {
	return STATE_USE_LBFGS;
    } else if (!strcmp(s, FORCE_DECP)) {
	return STATE_FORCE_DECPOINT;
    } else if (!strcmp(s, USE_CWD)) {
	return STATE_USE_CWD;
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
    } else if (!strcmp(s, ROBUST_Z)) {
	return STATE_ROBUST_Z;
    } else if (!strcmp(s, MWRITE_G)) {
	return STATE_MWRITE_G;
    } else if (!strcmp(s, MPI_USE_SMT)) {
	return STATE_MPI_SMT;
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
    } else if (!strcmp(s, USE_QR)) {
        set_flag_from_env(STATE_USE_QR, "GRETL_USE_QR", 0);
    } else if (!strcmp(s, USE_LBFGS)) {
	set_flag_from_env(STATE_USE_LBFGS, "GRETL_USE_LBFGS", 0);
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
    if (!strcmp(key, R_FUNCTIONS) && val != 0) {
	/* this depends on having R_lib on, so in
	   case it's off we should turn it on too
	*/
	libset_set_bool(R_LIB, val);
    }
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

#define PPDEBUG 0

static int n_states;
static GPtrArray *state_stack;
static int state_idx = -1;

#if PPDEBUG
static void print_state_stack (int pop)
{
    set_vars *sv;
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
    set_vars *newstate;
    int err = 0;

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

#if PPDEBUG
    print_state_stack(0);
#endif

    return err;
}

static void free_state (set_vars *sv)
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
	int fdp = state->flags & STATE_FORCE_DECPOINT;

	state_idx--;
	state = g_ptr_array_index(state_stack, state_idx);

	if (fdp && !(state->flags & STATE_FORCE_DECPOINT)) {
	    libset_set_decpoint(0);
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

/* support for GUI "Stop" button */

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
