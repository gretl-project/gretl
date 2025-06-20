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

#ifndef LIBGRETL_H
#define LIBGRETL_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef HAVE_VASPRINTF
# define _GNU_SOURCE
# define _ISOC99_SOURCE
# include <stdio.h>
# undef _GNU_SOURCE
#else
# include <stdio.h>
#endif

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>

#include <zlib.h>

#define GLIB_VERSION_MIN_REQUIRED GLIB_VERSION_2_28
#define GLIB_VERSION_MAX_ALLOWED (GLIB_VERSION_CUR_STABLE)
#include <glib.h>

#ifdef G_OS_WIN32
/* set our non-standard Windows identifier, if
   it's not already defined by the compiler
*/
# ifndef WIN32
#  define WIN32
# endif
#endif

#ifdef FULL_XML_HEADERS
# include <libxml/xmlmemory.h>
# include <libxml/parser.h>
#endif

#ifdef __APPLE__
# include <TargetConditionals.h>
#endif

#ifdef  __cplusplus
extern "C" {
#endif

#ifdef WIN32
# ifndef isnan
#  define isnan(x) ((x) != (x))
# endif
#endif

#define YMD_READ_FMT     "%d-%d-%d" 
#define YMD_WRITE_FMT    "%d-%02d-%02d"
#define YMD_WRITE_Y2_FMT "%02d-%02d-%02d"
#define YMD_WRITE_Y4_FMT "%04d-%02d-%02d"

#ifdef ENABLE_NLS
# include "libintl.h"
# include "locale.h"
# define gettext_noop(String) String
# define _(String) gettext (String)
# define N_(String) gettext_noop (String)
#else /* no NLS */
# define _(String)  ((char *) String)
# define N_(String) String
#endif /* NLS or not */

#define MAXLINE 131072 /* maximum length of command line */
#define MAXLABEL   128 /* maximum length of descriptive labels for variables */
#define MAXLEN     512 /* max length of regular "long" strings */
#define MAXDISP     32 /* max length of "display names" for variables */
#define VNAMELEN    32 /* space allocated for var names (including termination) */
#define OBSLEN      16 /* space allocated for obs strings (including termination) */

#ifndef M_PI
# define M_PI 3.1415926535897932384626432
#endif

#ifndef M_2PI
# define M_2PI 6.2831853071795864769252864
#endif

#define SQRT_2_PI     2.506628274631000502415765284811
#define LN_PI         1.144729885849400174143427351353
#define LN_2_PI       1.837877066409345483560659472811
#define LN_SQRT_2_PI  0.918938533204672741780329736406

#define PMAX_NOT_AVAILABLE 666

/* numbers smaller than the given limit will print as zero */
#define screen_zero(x)  ((fabs(x) > 1.0e-13)? x : 0.0)

typedef enum {
    GRETL_TYPE_NONE,
    GRETL_TYPE_BOOL,
    GRETL_TYPE_INT,
    GRETL_TYPE_UNSIGNED,
    GRETL_TYPE_OBS,
    GRETL_TYPE_LIST,
    GRETL_TYPE_DOUBLE,
    GRETL_TYPE_INT_ARRAY,
    GRETL_TYPE_DOUBLE_ARRAY,
    GRETL_TYPE_STRING,
    GRETL_TYPE_CMPLX_ARRAY,
    GRETL_TYPE_SERIES,
    GRETL_TYPE_MATRIX,
    GRETL_TYPE_STRUCT,
    GRETL_TYPE_SCALAR_REF,
    GRETL_TYPE_SERIES_REF,
    GRETL_TYPE_MATRIX_REF,
    GRETL_TYPE_STRING_REF,
    GRETL_TYPE_LIST_REF,
    GRETL_TYPE_USERIES,
    GRETL_TYPE_DATE,
    GRETL_TYPE_BUNDLE,
    GRETL_TYPE_BUNDLE_REF,
    GRETL_TYPE_ARRAY,
    GRETL_TYPE_ARRAY_REF,
    GRETL_TYPE_STRINGS,
    GRETL_TYPE_MATRICES,
    GRETL_TYPE_BUNDLES,
    GRETL_TYPE_LISTS,
    GRETL_TYPE_ARRAYS,
    GRETL_TYPE_STRINGS_REF,
    GRETL_TYPE_MATRICES_REF,
    GRETL_TYPE_BUNDLES_REF,
    GRETL_TYPE_LISTS_REF,
    GRETL_TYPE_ARRAYS_REF,
    GRETL_TYPE_VOID,
    GRETL_TYPE_NUMERIC,
    GRETL_TYPE_ANY
} GretlType;

#define gretl_scalar_type(t) (t == GRETL_TYPE_BOOL ||	\
                              t == GRETL_TYPE_INT ||	\
			      t == GRETL_TYPE_UNSIGNED || \
			      t == GRETL_TYPE_OBS ||	\
			      t == GRETL_TYPE_DOUBLE)

#define gretl_ref_type(t) (t == GRETL_TYPE_SCALAR_REF ||   \
	                   t == GRETL_TYPE_SERIES_REF ||   \
	                   t == GRETL_TYPE_MATRIX_REF ||   \
	                   t == GRETL_TYPE_BUNDLE_REF ||   \
			   t == GRETL_TYPE_STRING_REF ||   \
			   t == GRETL_TYPE_STRINGS_REF ||  \
			   t == GRETL_TYPE_MATRICES_REF || \
			   t == GRETL_TYPE_BUNDLES_REF || \
			   t == GRETL_TYPE_LISTS_REF || \
			   t == GRETL_TYPE_ARRAYS_REF)

#define gretl_array_type(t) (t == GRETL_TYPE_STRINGS ||  \
			     t == GRETL_TYPE_MATRICES || \
			     t == GRETL_TYPE_BUNDLES || \
			     t == GRETL_TYPE_LISTS || \
			     t == GRETL_TYPE_ARRAYS)

enum ts_codes {
    CROSS_SECTION,
    TIME_SERIES,
    STACKED_TIME_SERIES,
    STACKED_CROSS_SECTION,
    PANEL_UNKNOWN,
    PANEL_SIDE_BY_SIDE,
    SPECIAL_TIME_SERIES,
    STRUCTURE_UNKNOWN
};

enum progress_flags {
    SP_NONE,
    SP_TOTAL,
    SP_LOAD_INIT,
    SP_SAVE_INIT,
    SP_FONT_INIT,
    SP_UPDATER_INIT,
    SP_FINISH 
};

enum progress_return_flags {
    SP_RETURN_OK,
    SP_RETURN_DONE,
    SP_RETURN_CANCELED
};

enum test_stats {
    GRETL_STAT_NONE,
    GRETL_STAT_NORMAL_CHISQ,
    GRETL_STAT_LM,
    GRETL_STAT_F,
    GRETL_STAT_LMF,
    GRETL_STAT_HARVEY_COLLIER,
    GRETL_STAT_RESET,
    GRETL_STAT_LR,
    GRETL_STAT_WALD_CHISQ,
    GRETL_STAT_SUP_WALD,
    GRETL_STAT_Z,
    GRETL_STAT_STUDENT,
    GRETL_STAT_LB_CHISQ,
    GRETL_STAT_WF
};

typedef enum {
    OPT_NONE = 0,
    OPT_A = 1 <<  0,
    OPT_B = 1 <<  1,
    OPT_C = 1 <<  2,
    OPT_D = 1 <<  3,
    OPT_E = 1 <<  4,
    OPT_F = 1 <<  5,
    OPT_G = 1 <<  6,
    OPT_H = 1 <<  7,
    OPT_I = 1 <<  8,
    OPT_J = 1 <<  9,
    OPT_K = 1 << 10,
    OPT_L = 1 << 11,
    OPT_M = 1 << 12,
    OPT_N = 1 << 13,
    OPT_O = 1 << 14,
    OPT_P = 1 << 15,
    OPT_Q = 1 << 16,
    OPT_R = 1 << 17,
    OPT_S = 1 << 18,
    OPT_T = 1 << 19,
    OPT_U = 1 << 20,
    OPT_V = 1 << 21,
    OPT_W = 1 << 22,
    OPT_X = 1 << 23,
    OPT_Y = 1 << 24,
    OPT_Z = 1 << 25,
    OPT_a = 1 << 26,
    OPT_b = 1 << 27,
    OPT_i = 1 << 28,
    OPT_UNSET = 1 << 30
} gretlopt;

typedef enum {
    OP_EQ  = '=',
    OP_GT  = '>',
    OP_LT  = '<',
    OP_NEQ = 21,
    OP_GTE = 22,
    OP_LTE = 23
} GretlOp;

typedef enum {
    C_AIC,
    C_BIC,
    C_HQC,
    C_MAX
} ModelSelCriteria;

typedef enum {
    FC_STATIC,
    FC_DYNAMIC,
    FC_AUTO,
    FC_KSTEP
} ForecastMethod;

#ifndef CMPLX
typedef struct _cmplx {
    double r;
    double i;
} cmplx;
#endif

typedef struct GRETL_VAR_ GRETL_VAR;
typedef struct FITRESID_ FITRESID;

typedef struct model_data_item_ model_data_item;
typedef struct ModelTest_ ModelTest;
typedef struct equation_system_ equation_system;

typedef struct gretl_bundle_ gretl_bundle;
typedef struct gretl_array_ gretl_array;
typedef struct gretl_string_table_ gretl_string_table;

/**
 * VARINFO:
 *
 * Holds extended metadata on an individual data series.
 */

typedef struct VARINFO_ VARINFO;

/* gretl data set */
typedef struct DATASET_ { 
    int v;              /* number of data series, including const */
    int n;              /* number of observations */
    int pd;             /* periodicity or frequency of data */
    int structure;      /* time series, cross section or whatever */
    double sd0;         /* floating-point representation of stobs */
    int t1, t2;         /* start and end of current sample */
    char stobs[OBSLEN];  /* string representation of starting obs */
    char endobs[OBSLEN]; /* string representation of ending obs */
    double **Z;         /* two-dimensional data array */
    char **varname;     /* array of names of series */
    VARINFO **varinfo;  /* array of metadata per series */
    char markers;       /* NO_MARKERS (0), REGULAR MARKERS or DAILY_DATE_STRINGS */
    char modflag;       /* binary flag for dataset modified or not */
    char **S;           /* to hold observation markers */
    char *descrip;      /* to hold info on data sources, etc. */
    char *submask;      /* subsampling mask */
    char *restriction;  /* record of sub-sampling restriction */
    char *padmask;      /* record of padding to re-balance panel data */
    char *mapfile;      /* name of associated map (polygons) file, if any */
    unsigned int rseed; /* resampling seed */
    int auxiliary;      /* 0 for regular dataset, 1 for auxiliary dataset */
    int n_varinfo;      /* # of named series with metadata, if < @v above */
    char *pangrps;      /* panel-only: name of series holding group names */
    int panel_pd;       /* panel-only: panel time-series frequency */
    double panel_sd0;   /* panel-only: time-series start */
} DATASET;

typedef struct VMatrix_ {
    int ci;
    int dim;
    int t1, t2;
    int nmin, nmax;
    int ncrit; /* for assessing critical values */
    char **names;
    double *vec;
    double *xbar;
    double *ssx;
    int *list;
    int missing;
} VMatrix;

typedef struct SAMPLE_ {
    int t1;
    int t2;
    unsigned int rseed;
} SAMPLE;

typedef struct ARINFO_ {
    int *arlist;          /* list of autoregressive lags */
    double *rho;          /* array of autoreg. coeffs. */
    double *sderr;        /* and their standard errors */
} ARINFO;

/* struct to hold model results */
typedef struct MODEL_ {
    int ID;                      /* ID number for model */
    int refcount;                /* for saving/deleting */
    int ci;                      /* "command index" -- estimation method */
    gretlopt opt;                /* record of options */
    int t1, t2, nobs;            /* starting observation, ending
                                    observation, and number of obs */
    char *submask;               /* keep track of sub-sample in force
                                    when model was estimated */
    char *missmask;              /* missing observations mask */
    SAMPLE smpl;                 /* numeric start and end of current sample
                                    when model was estimated */
    int full_n;                  /* full length of dataset on estimation */
    int ncoeff, dfn, dfd;        /* number of coefficents; degrees of
                                    freedom in numerator and denominator */
    int *list;                   /* list of variables by ID number */
    int ifc;                     /* = 1 if the equation includes a constant,
				    else = 0 */
    int nwt;                     /* ID number of the weight variable (WLS) */
    int aux;                     /* code representing the sort of
				    auxiliary regression this is (or not) */
    double *coeff;               /* array of coefficient estimates */
    double *sderr;               /* array of estimated std. errors */
    double *uhat;                /* regression residuals */
    double *yhat;                /* fitted values from regression */
    double *xpx;                 /* X'X matrix, in packed form */
    double *vcv;                 /* VCV matrix for coefficient estimates */
    double ess, tss;             /* Error and Total Sums of Squares */
    double sigma;                /* Standard error of regression */
    double rsq, adjrsq;          /* Unadjusted and adjusted R^2 */     
    double fstt;                 /* overall F-statistic */
    double chisq;                /* overall chi-square statistic */
    double lnL;                  /* log-likelihood */
    double ybar, sdy;            /* mean and std. dev. of dependent var. */
    double criterion[C_MAX];     /* array of model selection statistics */
    double dw, rho;              /* Durbin-Watson stat. and estimated 1st
				    order autocorrelation coefficient */
    ARINFO *arinfo;              /* pointer to struct to hold special info for 
				    autoregressive model */ 
    int errcode;                 /* Error code in case of failure */
    char *name;                  /* for use in GUI */
    char *depvar;                /* name of dependent var in special cases */
    int nparams;                 /* number of named model parameters */
    char **params;               /* for named model parameters */
    gint64 esttime;              /* time of estimation */
    int ntests;                  /* number of attached test results */
    ModelTest *tests;            /* attached hypothesis test results */
    DATASET *dataset;            /* for handling models estimated on a
				    sub-sampled portion of the dataset */
    int n_data_items;            /* number of extra data items */
    model_data_item **data_items; /* pointer to additional data */
} MODEL;

#ifdef __ARM_ARCH_ISA_A64
# include <complex.h>
#endif

#include "gretl_commands.h"
#include "gretl_prn.h"
#include "gretl_errors.h"
#include "interact.h"
#include "dataset.h"
#include "estimate.h"
#include "genmain.h"
#include "genfuncs.h"
#include "compare.h"
#include "gretl_bundle.h"
#include "gretl_array.h"    
#include "gretl_intl.h"
#include "gretl_list.h"
#include "gretl_paths.h"
#include "gretl_utils.h"
#include "gretl_model.h"
#include "pvalues.h"
#include "dataio.h"
#include "gretl_data_io.h"
#include "strutils.h"
#include "describe.h"
#include "printout.h"
#include "printscan.h"
#include "modelprint.h"
#include "graphing.h"
#include "random.h"
#include "nonparam.h"
#include "options.h"
#include "discrete.h"
#include "adf_kpss.h"
#include "subsample.h"
#include "calendar.h"
#include "plugins.h"
#include "nls.h"
#include "missing.h"
#include "transforms.h"

#ifdef  __cplusplus
}
#endif

#endif /* LIBGRETL_H */
