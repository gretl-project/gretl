/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef LIBGRETL_H
#define LIBGRETL_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

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
# ifdef USE_GTK2
#  define I_(String) iso_gettext (String) 
#  define M_(String) maybe_iso_gettext (String)
# else
#  define I_(String) _(String)
#  define M_(String) _(String)
# endif /* USE_GTK2 */
#else
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

#define MAXLABEL 128  /* maximum length of descriptive labels for variables */
#define MAXLEN   512  /* max length of "long" strings */
#define ERRLEN   256  /* max length of libgretl error messages */
#define MAXDISP   20  /* max length of "display names" for variables */
#define VNAMELEN   9  /* space allocated for var names (including termination) */
#define OBSLEN    11  /* space allocated for obs strings (including termination) */

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

#define LISTSEP            999
#define PMAX_NOT_AVAILABLE 666

/* numbers smaller than the given limit will print as zero */
#define screen_zero(x)  ((fabs(x) > 1.0e-13)? x : 0.0)

typedef enum {
    GRETL_PRINT_STDOUT,
    GRETL_PRINT_STDERR,
    GRETL_PRINT_FILE,
    GRETL_PRINT_BUFFER,
    GRETL_PRINT_NULL
} prn_codes;

typedef enum {
    GRETL_PRINT_FORMAT_PLAIN,
    GRETL_PRINT_FORMAT_TEX,
    GRETL_PRINT_FORMAT_TEX_DOC,
    GRETL_PRINT_FORMAT_RTF,
    GRETL_PRINT_FORMAT_FIXED
} gretl_print_formats;

typedef enum {
    TIME_SERIES = 1,
    STACKED_TIME_SERIES,
    STACKED_CROSS_SECTION
} ts_codes;

typedef enum {
    GENR_RESID,
    GENR_FITTED,
    GENR_RESID2
} auto_genr;

typedef enum {
    SP_NONE, 
    SP_LOAD_INIT,
    SP_SAVE_INIT,
    SP_FONT_INIT,
    SP_UPDATER_INIT,
    SP_FINISH 
} progress_flags;

typedef enum {
    SP_RETURN_OK,
    SP_RETURN_DONE,
    SP_RETURN_CANCELED
} progress_return_flags;

typedef enum {
    GRETL_TEST_NORMAL_CHISQ,
    GRETL_TEST_TR2,
    GRETL_TEST_F,
    GRETL_TEST_LMF,
    GRETL_TEST_HARVEY_COLLIER,
    GRETL_TEST_RESET
} test_stats;

enum {
    OPT_A = 1 <<  0,
    OPT_B = 1 <<  1,
    OPT_C = 1 <<  2,
    OPT_D = 1 <<  3,
    OPT_F = 1 <<  4,
    OPT_I = 1 <<  5,
    OPT_L = 1 <<  6,
    OPT_M = 1 <<  7,
    OPT_N = 1 <<  8,
    OPT_O = 1 <<  9,
    OPT_Q = 1 << 10,
    OPT_R = 1 << 11,
    OPT_S = 1 << 12,
    OPT_T = 1 << 13,
    OPT_V = 1 << 14,
    OPT_W = 1 << 15,
    OPT_X = 1 << 16,
    OPT_Z = 1 << 17
} gretl_opt_flags;

typedef int *LIST; 

typedef struct _VARINFO VARINFO;
typedef struct _DATAINFO DATAINFO;
typedef struct _PATHS PATHS;
typedef struct _GRETLTEST GRETLTEST;
typedef struct _GRETLSUMMARY GRETLSUMMARY;
typedef struct _CORRMAT CORRMAT;
typedef struct _SAMPLE SAMPLE;
typedef struct _ARINFO ARINFO;
typedef struct _MODEL MODEL;
typedef struct _MODELSPEC MODELSPEC;
typedef struct _GRAPHT GRAPHT;
typedef struct _PRN PRN;
typedef struct _FITRESID FITRESID;
typedef struct _CONFINT CONFINT;
typedef struct _VCV VCV;
typedef struct _DATASET DATASET;
typedef struct _GRETL_VAR GRETL_VAR;

typedef struct _mp_results mp_results;
typedef struct _model_data_item model_data_item;


/* information on individual variable */
struct _VARINFO {
    char label[MAXLABEL];
    char display_name[MAXDISP];
    int compact_method;
};

/* information on data set */
struct _DATAINFO { 
    int v;              /* number of variables */
    int n;              /* number of observations */
    int pd;             /* periodicity or frequency of data */
    double sd0;         /* float representation of stobs */
    int t1, t2;         /* start and end of current sample */
    char stobs[OBSLEN];  /* string representation of starting obs (date) */
    char endobs[OBSLEN]; /* string representation of ending obs */
    char **varname;     /* array of names of variables */
    VARINFO **varinfo;  /* array of specific info on vars */
    char markers;       /* whether (1) or not (0) the data file has
			   observation markers */
    char delim;         /* default delimiter for "CSV" files */
    char time_series;   /* bit field to record time-series/panel nature
			   of the data set */
    char decpoint;      /* character used to represent decimal point */
    char **S;           /* to hold observation markers */
    char *descrip;      /* to hold info on data sources etc. */
    unsigned char *vector; /* hold info on vars: vector versus scalar */
    void *data;         /* all-purpose pointer */
};

/* wrapper for the two main elements of a gretl data set */
struct _DATASET {
    DATAINFO *dinfo;
    double **Z;
};

struct _PATHS {
    char currdir[MAXLEN];
    char userdir[MAXLEN];
    char gretldir[MAXLEN];
    char datadir[MAXLEN];
    char scriptdir[MAXLEN];
    char helpfile[MAXLEN];
    char cmd_helpfile[MAXLEN];
    char datfile[MAXLEN];
    char plotfile[MAXLEN];
    char binbase[MAXLEN];
    char ratsbase[MAXLEN];
    char gnuplot[MAXLEN];
#ifdef HAVE_X12A
    char x12a[MAXLEN];
    char x12adir[MAXLEN];
#endif
    char dbhost_ip[16];
    char pngfont[16];
};

struct _GRETLTEST {
    char type[72];
    char h_0[64];
    char param[VNAMELEN];
    unsigned char teststat;
    int dfn, dfd;
    double value;
    double pvalue;
};

struct _GRETLSUMMARY {
    int n;
    int *list;
    double *xskew; 
    double *xkurt;
    double *xmedian;
    double *coeff;
    double *sderr;
    double *xpx;
    double *xpy;
};

struct _CORRMAT {
    int n, t1, t2;
    int *list;
    double *xpx;
};

struct _SAMPLE {
    int t1;
    int t2;
};

struct _ARINFO {
    int *arlist;                /* list of autoreg lags */
    double *rho;                /* array of autoreg. coeffs. */
    double *sderr;              /* and their standard errors */
};

/* struct to hold model results */
struct _MODEL {
    int ID;                      /* ID number for model */
    int t1, t2, nobs;            /* starting observation, ending
                                    observation, and number of obs */
    double *subdum;              /* keep track of sub-sample in force
                                    when model was estimated */
    SAMPLE smpl;                 /* numeric start and end of current sample
                                    when model was estimated */
    int ncoeff, dfn, dfd;        /* number of coefficents; degrees of
                                    freedom in numerator and denominator */
    int *list;                   /* list of variables by ID number */
    int ifc;                     /* = 1 if the equation includes a constant,
				    else = 0 */
    int ci;                      /* "command index" -- depends on 
				    estimation method */
    int nwt;                     /* ID of the weight variable (WLS) */
    int order;                   /* lag order (e.g. for ARCH) */
    int aux;                     /* code representing the sort of
				    auxiliary regression this is (or not) */
    double *coeff;               /* array of coefficient estimates */
    double *sderr;               /* array of estimated std. errors */
    double *uhat;                /* regression residuals */
    double *yhat;                /* fitted values from regression */
    double *xpx;
    double *vcv;                 /* VCV matrix for coeff. estimates */
    double ess, tss;             /* Error and Total Sums of Squares */
    double sigma;                /* Standard error of regression */
    double ess_wt;               /* ESS using weighted data (WLS) */
    double sigma_wt;             /* same thing for std.error */
    double rsq, adjrsq;          /* Unadjusted and adjusted R^2 */     
    double fstt;                 /* F-statistic */
    double lnL;                  /* log-likelihood */
    double ybar, sdy;            /* mean and std. dev. of dependent var. */
    double criterion[8];         /* array of model selection statistics */
    double dw, rho;              /* Durbin-Watson stat. and estimated 1st
				    order autocorrelation coefficient */
    ARINFO *arinfo;              /* struct to hold special info for 
				    autoregressive model */ 
    int errcode;                 /* Error code in case of failure */
    char *name;
    int nparams;                 /* number of named model parameters */
    char **params;               /* for named model parameters */
    int ntests;
    GRETLTEST *tests;
    void *data;                  /* pointer for use in re. missing data */
    DATASET *dataset;            /* for handling models estimated on a
				    sub-sampled portion of the dataset */
    int n_data_items;            /* number of extra data items */
    model_data_item **data_items; /* pointer to additional data */
};

struct _MODELSPEC {
    char *cmd;
    double *subdum;
};

struct _GRAPHT {
    int ID;
    int sort;
    char name[32];
    char fname[MAXLEN];
}; 

struct _PRN {
    FILE *fp;
    FILE *fpaux;    
    char *buf;
    size_t bufsize;
    int format;
};

struct _mp_results {
    int ncoeff;
    int t1, t2, ifc;
    int dfn, dfd;
    int *varlist;
    char **varnames;
    double *coeff;
    double *sderr;
    double sigma;
    double ess;
    double rsq, adjrsq;
    double fstt;
};

struct _FITRESID {
    double *actual;
    double *fitted;
    double *sderr;
    double sigma;
    double tval;
    int pmax;
    int df;
    int t1, t2;
    int nobs;
    int err;
    char depvar[VNAMELEN];
};

struct _CONFINT {
    int *list;
    double *coeff;
    double *maxerr;
    int df;
    int ifc;
};

struct _VCV {
    int ci;
    int *list;
    double *vec;
};

#define VARLABEL(p,i)  ((p->varinfo[i])->label)
#define DISPLAYNAME(p,i)  ((p->varinfo[i])->display_name)
#define COMPACT_METHOD(p,i) ((p->varinfo[i])->compact_method)

#include "gretl_commands.h"
#include "gretl_errors.h"
#include "estimate.h"
#include "generate.h"
#include "gretl_paths.h"
#include "gretl_utils.h"
#include "gretl_model.h"
#include "compare.h"
#include "pvalues.h"
#include "dataio.h"
#include "dbread.h"
#include "strutils.h"
#include "describe.h"
#include "printout.h"
#include "modelprint.h"
#include "texprint.h"
#include "interact.h"
#include "graphing.h"
#include "random.h"
#include "monte_carlo.h"
#include "nonparam.h"
#include "discrete.h"
#include "subsample.h"
#include "calendar.h"
#include "system.h"
#include "plugins.h"
#include "nls.h"
#ifdef WIN32
# include "gretl_win32.h"
#endif

int parse_set_line (const char *line, int *echo_off); /* libset.c */
void set_use_qr (int set);

#endif /* LIBGRETL_H */
