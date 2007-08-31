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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>

#include <zlib.h>

#ifndef GRETLCLI
# include <libxml/xmlmemory.h>
# include <libxml/parser.h>
#endif

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef WIN32
# include "winconfig.h"
#endif

#ifdef WIN32
# ifndef isnan
#  define isnan(x) ((x) != (x))
# endif
#endif

#ifdef ENABLE_NLS
# define G_(String) gp_gettext (String)
# define I_(String) iso_gettext (String) 
# define M_(String) maybe_iso_gettext (String)
#else
# define G_(String) String
# define I_(String) String
# define M_(String) String
#endif /* ENABLE_NLS */

#ifndef __GNOME_I18N_H__
# ifdef ENABLE_NLS
#  include "libintl.h"
#  include "locale.h"
#  define gettext_noop(String) String
#  define _(String) gettext (String)
#  define N_(String) gettext_noop (String)
# else
#  define _(String) String
#  define N_(String) String
# endif /* ENABLE_NLS */
#endif /* __GNOME_I18N_H__ */

#define MAXLINE 8192  /* maximum length of command line */
#define MAXLABEL 128  /* maximum length of descriptive labels for variables */
#define MAXLEN   512  /* max length of regular "long" strings */
#define ERRLEN   256  /* max length of libgretl error messages */
#define MAXDISP   32  /* max length of "display names" for variables */
#define VNAMELEN  16  /* space allocated for var names (including termination) */
#define OBSLEN    16  /* space allocated for obs strings (including termination) */

#ifndef M_PI
# define M_PI 3.1415926535897932384626432
#endif

#ifndef M_2PI
# define M_2PI 6.2831853071795864769252864
#endif

#define LN_2_PI       1.837877066409345
#define LN_SQRT_2_PI  0.9189385332056725

#define LISTSEP            999
#define PMAX_NOT_AVAILABLE 666

/* numbers smaller than the given limit will print as zero */
#define screen_zero(x)  ((fabs(x) > 1.0e-13)? x : 0.0)

typedef enum {
    GRETL_TYPE_NONE,
    GRETL_TYPE_BOOL,
    GRETL_TYPE_INT,
    GRETL_TYPE_LIST,
    GRETL_TYPE_DOUBLE,
    GRETL_TYPE_INT_ARRAY,
    GRETL_TYPE_DOUBLE_ARRAY,
    GRETL_TYPE_STRING,
    GRETL_TYPE_CHAR_ARRAY,
    GRETL_TYPE_CMPLX_ARRAY,
    GRETL_TYPE_SERIES,
    GRETL_TYPE_MATRIX,
    GRETL_TYPE_STRUCT,
    GRETL_TYPE_SCALAR_REF,
    GRETL_TYPE_SERIES_REF,
    GRETL_TYPE_MATRIX_REF
} GretlType;

#define gretl_scalar_type(t) (t == GRETL_TYPE_BOOL || \
                              t == GRETL_TYPE_INT || \
			      t == GRETL_TYPE_DOUBLE)

#define gretl_ref_type(t) (t == GRETL_TYPE_SCALAR_REF || \
			   t == GRETL_TYPE_SERIES_REF || \
			   t == GRETL_TYPE_MATRIX_REF)

enum ts_codes {
    CROSS_SECTION,
    TIME_SERIES,
    STACKED_TIME_SERIES,
    STACKED_CROSS_SECTION,
    PANEL_UNKNOWN,
    SPECIAL_TIME_SERIES,
    STRUCTURE_UNKNOWN
};

enum auto_genr {
    GENR_RESID,
    GENR_FITTED,
    GENR_RESID2,
    GENR_H,
    GENR_AHAT
};

enum progress_flags {
    SP_NONE, 
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
    GRETL_STAT_TR2,
    GRETL_STAT_F,
    GRETL_STAT_LMF,
    GRETL_STAT_HARVEY_COLLIER,
    GRETL_STAT_RESET,
    GRETL_STAT_LR,
    GRETL_STAT_WALD_CHISQ,
    GRETL_STAT_SUP_WALD
};

enum gretl_opt_flags {
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
    OPT_Z = 1 << 24,
    OPT_UNSET = 1 << 30
};

typedef enum {
    OP_EQ  = '=',
    OP_GT  = '>',
    OP_LT  = '<',
    OP_NEQ = 21,
    OP_GTE = 22,
    OP_LTE = 23
} GretlOp;

typedef enum {
    D_NONE = 0,
    D_UNIFORM,
    D_NORMAL,
    D_GAMMA
} DistCode;

typedef enum {
    C_AIC,
    C_BIC,
    C_HQC,
    C_MAX
} ModelSelCriteria;

#ifndef CMPLX
typedef struct _cmplx cmplx;
struct _cmplx {
    double r;
    double i;
};
#endif

typedef unsigned long gretlopt;

typedef struct VARINFO_ VARINFO;
typedef struct PANINFO_ PANINFO;
typedef struct DATAINFO_ DATAINFO;
typedef struct PATHS_ PATHS;
typedef struct VMatrix_ VMatrix;
typedef struct SAMPLE_ SAMPLE;
typedef struct ARINFO_ ARINFO;
typedef struct MODEL_ MODEL;
typedef struct PRN_ PRN;
typedef struct FITRESID_ FITRESID;
typedef struct DATASET_ DATASET;
typedef struct GRETL_VAR_ GRETL_VAR;

typedef struct model_data_item_ model_data_item;
typedef struct ModelTest_ ModelTest;
typedef struct gretl_equation_system_ gretl_equation_system;

/* information on individual variable */
struct VARINFO_ {
    char label[MAXLABEL];
    char display_name[MAXDISP];
    char gp_display_name[MAXDISP];
    int flags;
    char compact_method;
    char stack_level;
    char line_width;
    char **sorted_markers;
};

/* panel indexing information */
struct PANINFO_ {
    int nunits;    /* number of cross-sectional units */
    int Tmin;      /* min. number of time-series observations per unit */
    int Tmax;      /* max. number of time-series observations per unit */
    int olen;      /* length in digits of highest time-series index value */
    int *unit;     /* index array, cross-sectional units */
    int *period;   /* index array, time periods */
    char *padmask; /* mask recording padding, when subsampled */
};

/* information on data set */
struct DATAINFO_ { 
    int v;              /* number of variables */
    int n;              /* number of observations */
    int pd;             /* periodicity or frequency of data */
    int structure;      /* time series, cross section or whatever */
    double sd0;         /* float representation of stobs */
    int t1, t2;         /* start and end of current sample */
    char stobs[OBSLEN];  /* string representation of starting obs (date) */
    char endobs[OBSLEN]; /* string representation of ending obs */
    char **varname;     /* array of names of variables */
    VARINFO **varinfo;  /* array of specific info on vars */
    PANINFO *paninfo;   /* additional info for panel data */
    char markers;       /* whether (1) or not (0) the data file has
			   observation markers */
    char delim;         /* default delimiter for "CSV" files */
    char decpoint;      /* character used to represent decimal point */
    char submode;       /* mode of sub-sampling in force, if any */
    char **S;           /* to hold observation markers */
    char *descrip;      /* to hold info on data sources etc. */
    char *submask;      /* subsampling mask */
    void *data;         /* all-purpose pointer */
};

/* wrapper for the two main elements of a gretl data set */
struct DATASET_ {
    DATAINFO *dinfo;
    double **Z;
};

struct PATHS_ {
    char currdir[MAXLEN];
    char userdir[MAXLEN];
    char gretldir[MAXLEN];
    char datadir[MAXLEN];
    char scriptdir[MAXLEN];
    char helpfile[MAXLEN];
    char cmd_helpfile[MAXLEN];
    char cli_helpfile[MAXLEN];
    char datfile[MAXLEN];
    char binbase[MAXLEN];
    char ratsbase[MAXLEN];
    char gnuplot[MAXLEN];
    char x12a[MAXLEN];
    char x12adir[MAXLEN];
    char tramo[MAXLEN];
    char tramodir[MAXLEN];
    char dbhost[32];
    char pngfont[128];
};

struct VMatrix_ {
    int ci;
    int dim;
    int t1, t2, n;
    char **names;
    double *vec;
    int *list;
    int missing;
};

struct SAMPLE_ {
    int t1;
    int t2;
};

struct ARINFO_ {
    int *arlist;          /* list of autoreg lags */
    double *rho;          /* array of autoreg. coeffs. */
    double *sderr;        /* and their standard errors */
};

/* struct to hold model results */
struct MODEL_ {
    int ID;                      /* ID number for model */
    int refcount;                /* for saving/deleting */
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
    int ci;                      /* "command index" -- depends on 
				    estimation method */
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
    double fstt;                 /* F-statistic */
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
    int ntests;                  /* number of attached test results */
    ModelTest *tests;            /* attached hypothesis test results */
    DATASET *dataset;            /* for handling models estimated on a
				    sub-sampled portion of the dataset */
    int n_data_items;            /* number of extra data items */
    model_data_item **data_items; /* pointer to additional data */
};

#define VARLABEL(p,i)        ((p->varinfo[i])->label)
#define DISPLAYNAME(p,i)     ((p->varinfo[i])->display_name)
#define COMPACT_METHOD(p,i)  ((p->varinfo[i])->compact_method)
#define STACK_LEVEL(p,i)     ((p->varinfo[i])->stack_level)
#define SORTED_MARKER(p,i,t) ((p->varinfo[i])->sorted_markers[t])

#include "gretl_commands.h"
#include "gretl_errors.h"
#include "interact.h"
#include "dataset.h"
#include "estimate.h"
#include "genmain.h"
#include "genfuncs.h"
#include "compare.h"
#include "gretl_intl.h"
#include "gretl_list.h"
#include "gretl_paths.h"
#include "gretl_utils.h"
#include "gretl_model.h"
#include "gretl_prn.h"
#include "pvalues.h"
#include "dataio.h"
#include "dbread.h"
#include "strutils.h"
#include "describe.h"
#include "printout.h"
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
#include "tsls.h"
#ifdef WIN32
# include "gretl_win32.h"
#endif

#endif /* LIBGRETL_H */
