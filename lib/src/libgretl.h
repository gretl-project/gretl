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
# include "../../config.h"
#endif

#ifdef OS_WIN32
# include "winconfig.h"
#endif

#ifdef OS_WIN32
# ifndef isnan
#  define isnan(x) ((x) != (x))
# endif
#endif

#ifdef ENABLE_NLS
# ifdef USE_GTK2
#  define I_(String) iso_gettext (String) 
# else
#  define I_(String) _(String)
# endif /* USE_GTK2 */
#else
# define I_(String) String
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

#define MAXLABEL 128  /* maximum length of decsriptive labels for variables */
#define MAXLEN   512  /* max length of "long" strings */
#define ERRLEN   256  /* max length of libgretl error messages */

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

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
    GRETL_PRINT_FORMAT_RTF
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
    SP_FINISH 
} progress_flags;

typedef enum {
    GRETL_TEST_NORMAL_CHISQ,
    GRETL_TEST_TR2,
    GRETL_TEST_F,
    GRETL_TEST_LMF,
    GRETL_TEST_HARVEY_COLLIER,
    GRETL_TEST_RESET
} test_stats;

typedef int *LIST;  

/* information on data set */
typedef struct { 
    int v;              /* number of variables */
    int n;              /* number of observations */
    int pd;             /* periodicity or frequency of data */
    int bin;            /* data file is binary? (1 or 0) */
    int extra;
    double sd0;         /* float representation of stobs */
    int t1, t2;         /* start and end of current sample */
    char stobs[9];      /* string representation of starting obs (date) */
    char endobs[9];     /* string representation of ending obs */
    char **varname;     /* array of names of variables */
    char **label;       /* array of descriptive labels for vars */
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
} DATAINFO;

typedef struct {
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
    char dbhost_ip[16];
} PATHS;

typedef struct {
    char type[72];
    char h_0[64];
    char param[9];
    unsigned char teststat;
    int dfn, dfd;
    double value;
    double pvalue;
} GRETLTEST;

typedef struct {
    int n;
    int *list;
    double *xskew, *xkurt, *xmedian, *coeff, *sderr, *xpx, *xpy;
} GRETLSUMMARY;

typedef struct {
    int n, t1, t2;
    int *list;
    double *xpx;
} CORRMAT;

typedef struct {
    int t1;
    int t2;
} SAMPLE;

typedef struct {
    int *arlist;                /* list of autoreg lags */
    double *rho;                /* array of autoreg. coeffs. */
    double *sderr;              /* and their standard errors */
} ARINFO;

/* struct to hold model results */
typedef struct {
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
    int wt_dummy;                /* Is the weight var a 0/1 dummy? */
    int order;                   /* lag order (e.g. for ARCH) */
    int aux;                     /* code representing the sort of
				    auxiliary regression this is (or not) */
    int ldepvar;                 /* = 1 if there's a lag of the
				    dependent variable on the RHS, else 0 */
    int correct;                 /* cases 'correct' (binary depvar) */
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
    double chisq;                /* Chi-square */
    double ybar, sdy;            /* mean and std. dev. of dependent var. */
    double criterion[8];         /* array of model selection statistics */
    double dw, rho;              /* Durbin-Watson stat. and estimated 1st
				    order autocorrelation coefficient */
    double rho_in;               /* the rho value input into the 
				    regression (e.g. CORC) */
    ARINFO *arinfo;              /* struct to hold special info for 
				    autoregressive model */ 
    double *slope;               /* for nonlinear models */
    int errcode;                 /* Error code in case of failure */
    char *name;
    int ntests;
    GRETLTEST *tests;
    void *data;                   /* pointer for use in re. missing data */
} MODEL;

typedef struct {
    char *cmd;
    double *subdum;
} MODELSPEC;

typedef struct {
    int ID;
    int sort;
    char name[24];
    char fname[MAXLEN];
} GRAPHT; 

typedef struct {
    char name[32];
    int nmodels;
    int ngraphs;
    MODEL **models;
    GRAPHT **graphs;
    char *notes;
} SESSION; 

typedef struct {
    int nmodels;
    int *model_ID;
    char **model_name;
} SESSIONBUILD;

typedef struct {
    FILE *fp;
    char *buf;
    size_t bufsize;
    int format;
} PRN;

typedef struct {
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
} mp_results;

typedef struct {
    double *actual;
    double *fitted;
    double *sderr;
    double sigma;
    double tval;
    int pmax;
    int df;
    int t1, t2;
    int nobs;
    char depvar[9];
} FITRESID;

typedef struct {
    int *list;
    double *coeff;
    double *maxerr;
    int df;
    int ifc;
} CONFINT;

typedef struct {
    int ci;
    int *list;
    double *vec;
} VCV;

#include "commands.h"
#include "errors.h"
#include "estimate.h"
#include "generate.h"
#include "gretl_utils.h"
#include "compare.h"
#include "pvalues.h"
#include "dataio.h"
#include "strutils.h"
#include "describe.h"
#include "printout.h"
#include "modelprint.h"
#include "texprint.h"
#include "var.h"
#include "interact.h"
#include "graphing.h"
#include "random.h"
#include "monte_carlo.h"
#include "nonparam.h"
#include "discrete.h"
#include "subsample.h"
#include "calendar.h"
#include "system.h"
#ifdef OS_WIN32
# include "gretl_win32.h"
#endif

#endif /* LIBGRETL_H */
