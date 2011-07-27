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

#include "libgretl.h"
#include "gretl_scalar.h"

#include <errno.h>

#define OPTDEBUG 0

#define vcv_opt_ok(c) (MODEL_COMMAND(c) || c == ADD || c == OMIT)

#define quiet_opt_ok(c) (MODEL_COMMAND(c) ||	\
			 c == ADD ||		\
			 c == ADF ||		\
			 c == ANOVA ||		\
			 c == APPEND ||		\
			 c == CHOW ||		\
			 c == COINT2 ||		\
			 c == CORRGM ||		\
			 c == CUSUM ||		\
			 c == DATA ||		\
			 c == DPANEL ||		\
			 c == ESTIMATE ||	\
			 c == FCAST ||		\
			 c == FOREIGN ||	\
                         c == FRACTINT ||       \
			 c == FREQ ||		\
			 c == KPSS ||		\
			 c == MODTEST ||	\
                         c == LEVINLIN ||       \
			 c == LOOP ||		\
                         c == MAHAL ||          \
			 c == NORMTEST ||	\
			 c == OLS ||		\
			 c == OMIT ||		\
			 c == OPEN ||		\
			 c == RESET ||		\
			 c == RESTRICT ||	\
			 c == RMPLOT ||		\
			 c == SYSTEM ||		\
			 c == VAR ||		\
			 c == VECM ||		\
			 c == XCORRGM)

struct gretl_option {
    int ci;              /* command index (gives context) */
    gretlopt o;          /* index of integer type */
    const char *longopt; /* "--"-style string representation of option */
    char parminfo;       /* 0 = option can never take a parameter,
                            1 = option may take a parameter,
                            2 = option requires a parameter 
			 */
};

struct flag_match {
    gretlopt o;
    unsigned char c;
};

/* Below: This is used as a one-way mapping from the long form
   to the index (e.g. OPT_Q), so a given index can have more than 
   one long-form counterpart, depending on context.  The last field
   flags whether the given option accepts (1), or requires (2), an
   accompanying parameter value.
*/

struct gretl_option gretl_opts[] = {
    { ADD,      OPT_B, "both", 0 },
    { ADD,      OPT_I, "silent", 0 },
    { ADD,      OPT_L, "lm", 0 },
    { ADF,      OPT_N, "nc", 0 }, 
    { ADF,      OPT_C, "c", 0 }, 
    { ADF,      OPT_D, "seasonals", 0 },
    { ADF,      OPT_R, "ctt", 0 },     
    { ADF,      OPT_T, "ct", 0 }, 
    { ADF,      OPT_V, "verbose", 0 },
    { ADF,      OPT_F, "difference", 0 },
    { ADF,      OPT_E, "test-down", 0 },
    { ADF,      OPT_G, "gls", 0 },
    { AR1,      OPT_B, "no-corc", 0 },
    { AR1,      OPT_H, "hilu", 0 },
    { AR1,      OPT_P, "pwe", 0 },
    { APPEND,   OPT_T, "time-series", 0 },
    { ARBOND,   OPT_A, "asymptotic", 0 },
    { ARBOND,   OPT_D, "time-dummies", 0 },
    { ARBOND,   OPT_H, "orthdev", 0 },
    { ARBOND,   OPT_T, "two-step", 0 },
    { ARMA,     OPT_C, "conditional", 0 },
    { ARMA,     OPT_E, "save-ehat", 0 },    
    { ARMA,     OPT_G, "opg", 0 },
    { ARMA,     OPT_H, "hessian", 0 },
    { ARMA,     OPT_L, "lbfgs", 0 },
    { ARMA,     OPT_N, "nc", 0 },    
    { ARMA,     OPT_V, "verbose", 0 },
    { ARMA,     OPT_X, "x-12-arima", 0 },
    { ARMA,     OPT_Y, "y-diff-only", 0 },
    { BIPROBIT, OPT_G, "opg", 0 },
    { BIPROBIT, OPT_R, "robust", 0 },
    { BIPROBIT, OPT_V, "verbose", 0 },
    { BIPROBIT, OPT_X, "save-xbeta", 0 },
    { BXPLOT,   OPT_O, "notches", 0 },
    { BXPLOT,   OPT_L, "outliers", 2 },
    { BXPLOT,   OPT_P, "panel", 0 },
    { BXPLOT,   OPT_U, "output", 2 },
    { BXPLOT,   OPT_Z, "factorized", 0 },
    { CHOW,     OPT_D, "dummy", 0 },
    { CLEAR,    OPT_D, "dataset", 0 },
    { CLEAR,    OPT_O, "other", 0 },
    { COINT,    OPT_E, "test-down", 0 },
    { COINT,    OPT_N, "nc", 0 },
    { COINT,    OPT_R, "ctt", 0 },     
    { COINT,    OPT_S, "skip-df", 0 },
    { COINT,    OPT_T, "ct", 0 },
    { COINT2,   OPT_A, "crt", 0 },
    { COINT2,   OPT_D, "seasonals", 0 },
    { COINT2,   OPT_N, "nc", 0 },
    { COINT2,   OPT_R, "rc", 0 },
    { COINT2,   OPT_C, "uc", 0 },
    { COINT2,   OPT_T, "ct", 0 },
    { COINT2,   OPT_S, "silent", 0 },
    { COINT2,   OPT_V, "verbose", 0 },
    { COINT2,   OPT_Y, "asy", 0 },
    { CORR,     OPT_K, "kendall", 0 },
    { CORR,     OPT_S, "spearman", 0 },
    { CORR,     OPT_U, "uniform", 0 },
    { CORR,     OPT_V, "verbose", 0 },
    { CUSUM,    OPT_R, "squares", 0 },
    { DATA,     OPT_O, "odbc", 0 },
    { DATAMOD,  OPT_P, "preserve", 0 },
    { DELEET,   OPT_D, "db", 0 },
    { DELEET,   OPT_L, "list", 0 },
    { DIFFTEST, OPT_G, "sign", 0 },
    { DIFFTEST, OPT_R, "rank-sum", 0 },
    { DIFFTEST, OPT_I, "signed-rank", 0 },
    { DIFFTEST, OPT_V, "verbose", 0 },
    { DISCRETE, OPT_R, "reverse", 0 },
    { DPANEL,   OPT_A, "asymptotic", 0 },
    { DPANEL,   OPT_D, "time-dummies", 0 },    
    { DPANEL,   OPT_L, "system", 0 },
    { DPANEL,   OPT_T, "two-step", 0 },
    { DPANEL,   OPT_V, "verbose", 0 },
    { DPANEL,   OPT_X, "dpdstyle", 0 },
    { DUMMIFY,  OPT_F, "drop-first", 0 },
    { DUMMIFY,  OPT_L, "drop-last", 0 },
    { DURATION, OPT_W, "weibull", 0 },
    { DURATION, OPT_E, "exponential", 0 },
    { DURATION, OPT_L, "loglogistic", 0 },
    { DURATION, OPT_Z, "lognormal", 0 },
    { DURATION, OPT_M, "medians", 0 },    
    { DURATION, OPT_G, "opg", 0 },
    { DURATION, OPT_R, "robust", 0 },
    { DURATION, OPT_V, "verbose", 0 },
#if 0 /* not yet */
    { EQNPRINT, OPT_F, "filename", 0 }, /* backward compatibility */
#endif
    { EQNPRINT, OPT_O, "complete", 0 },
    { EQNPRINT, OPT_T, "t-ratios", 0 },
    { TABPRINT, OPT_O, "complete", 0 },
    { TABPRINT, OPT_R, "rtf", 0 },
#if 0 /* not yet */
    { TABPRINT, OPT_F, "filename", 0 }, /* backward compatibility */
#endif
    { EQUATION, OPT_M, "multi", 0 },
    { ESTIMATE, OPT_I, "iterate", 0 },
    { ESTIMATE, OPT_M, "geomean", 0 },
    { ESTIMATE, OPT_N, "no-df-corr", 0 },
    { ESTIMATE, OPT_V, "verbose", 0 },
    { FCAST,    OPT_D, "dynamic", 0 },
    { FCAST,    OPT_M, "mean-y", 0 },
    { FCAST,    OPT_N, "no-stats", 0 },
    { FCAST,    OPT_S, "static", 0 },
    { FCAST,    OPT_R, "rolling", 0 },
    { FCAST,    OPT_O, "out-of-sample", 0 },
    { FCAST,    OPT_I, "integrate", 0 },
    { FCAST,    OPT_U, "plot", 1 },
    { FOREIGN,  OPT_D, "send-data", 0 },
    { FOREIGN,  OPT_V, "verbose", 0 },
    { FRACTINT, OPT_G, "gph", 0 },
    { FRACTINT, OPT_A, "all", 0 },
    { FREQ,     OPT_O, "gamma", 0 },
    { FREQ,     OPT_S, "silent", 0 },
    { FREQ,     OPT_Z, "normal", 0 },
    { FUNDEBUG, OPT_C, "continue", 0 },
    { FUNDEBUG, OPT_N, "next", 0 },
    { FUNDEBUG, OPT_Q, "quit", 0 },
    { GARCH,    OPT_A, "arma-init", 0 },
    { GARCH,    OPT_F, "fcp", 0 }, 
    { GARCH,    OPT_N, "nc", 0 }, 
    { GARCH,    OPT_R, "robust", 0 },
    { GARCH,    OPT_V, "verbose", 0 },
    { GARCH,    OPT_Z, "stdresid", 0 },
    { GMM,      OPT_I, "iterate", 0 },
    { GMM,      OPT_L, "lbfgs", 0 },
    { GMM,      OPT_T, "two-step", 0 },
    { GMM,      OPT_V, "verbose", 0 },
    { GNUPLOT,  OPT_D, "input", 2 },
    { GNUPLOT,  OPT_O, "with-lines", 0 },
    { GNUPLOT,  OPT_I, "inverse-fit", 0 },
    { GNUPLOT,  OPT_L, "loess-fit", 0 },
    { GNUPLOT,  OPT_Q, "quadratic-fit", 0 },
    { GNUPLOT,  OPT_N, "linear-fit", 0 },
    { GNUPLOT,  OPT_B, "cubic-fit", 0 },
    { GNUPLOT,  OPT_M, "with-impulses", 0 },
    { GNUPLOT,  OPT_S, "suppress-fitted", 0 },
    { GNUPLOT,  OPT_T, "time-series", 0 },
    { GNUPLOT,  OPT_Z, "dummy", 0 },
    { GNUPLOT,  OPT_C, "control", 0 },
    { GNUPLOT,  OPT_U, "output", 2 },
    { GNUPLOT,  OPT_Y, "single-yaxis", 0 },
    { GNUPLOT,  OPT_X, "matrix", 2 },
    { HECKIT,   OPT_M, "ml", 0 },
    { HECKIT,   OPT_R, "robust", 0 },
    { HECKIT,   OPT_T, "two-step", 0 },
    { HECKIT,   OPT_V, "verbose", 0 },
    { HELP,     OPT_F, "func", 0 },
    { INTREG,   OPT_R, "robust", 0 },
    { INTREG,   OPT_V, "verbose", 0 },
    { IVREG,    OPT_G, "gmm", 0 },
    { IVREG,    OPT_I, "iterate", 0 },
    { IVREG,    OPT_L, "liml", 0 },
    { IVREG,    OPT_R, "robust", 0 },  
    { IVREG,    OPT_S, "save", 0 },
    { IVREG,    OPT_T, "two-step", 0 },
    { IVREG,    OPT_W, "weights", 2 },
    { KALMAN,   OPT_C, "cross", 0 },
    { KALMAN,   OPT_D, "diffuse", 0 },
    { KPSS,     OPT_T, "trend", 0 },
    { KPSS,     OPT_D, "seasonals", 0 },
    { KPSS,     OPT_V, "verbose", 0 },
    { KPSS,     OPT_F, "difference", 0 },
    { LEVERAGE, OPT_S, "save", 0 },
    { LEVINLIN, OPT_N, "nc", 0 },
    { LEVINLIN, OPT_T, "ct", 0 },
    { MAKEPKG,  OPT_I, "index", 0 },
    { MAKEPKG,  OPT_T, "translations", 0 },
    { MODTEST,  OPT_A, "autocorr", 0 },
    { MODTEST,  OPT_B, "breusch-pagan", 0 },
    { MODTEST,  OPT_C, "comfac", 0 },
    { MODTEST,  OPT_H, "arch", 0 },
    { MODTEST,  OPT_L, "logs", 0 },
    { MODTEST,  OPT_N, "normality", 0 },
    { MODTEST,  OPT_S, "squares", 0 }, 
    { MODTEST,  OPT_P, "panel", 0 },
    { MODTEST,  OPT_R, "robust", 0 },
    { MODTEST,  OPT_W, "white", 0 },
    { MODTEST,  OPT_X, "white-nocross", 0 },
    { LABELS,   OPT_F, "from-file", 2 },
    { LABELS,   OPT_T, "to_file", 2 },
    { LOGIT,    OPT_M, "multinomial", 0 },
    { LOGIT,    OPT_P, "p-values", 0 },
    { LOGIT,    OPT_R, "robust", 0 },
    { LOGIT,    OPT_V, "verbose", 0 },
    { LOOP,     OPT_P, "progressive", 0 },
    { LOOP,     OPT_V, "verbose", 0 },
    { MAHAL,    OPT_S, "save", 0 },
    { MAHAL,    OPT_V, "vcv", 0 },
    { MEANTEST, OPT_O, "unequal-vars", 0 },
    { MLE,      OPT_A, "auxiliary", 0 },
    { MLE,      OPT_H, "hessian", 0 },
    { MLE,      OPT_G, "no-gradient-check", 0 },
    { MLE,      OPT_L, "lbfgs", 0 },
    { MLE,      OPT_N, "numerical", 0 },
    { MLE,      OPT_R, "robust", 0 },
    { MLE,      OPT_V, "verbose", 0 },
    { MODPRINT, OPT_A, "addstats", 2 },
    { MODPRINT, OPT_C, "csv", 0 },
    { MODPRINT, OPT_O, "complete", 0 },
    { MODPRINT, OPT_R, "rtf", 0 },
    { MODPRINT, OPT_T, "tex", 0 },
    { MPOLS,    OPT_O, "vcv", 0 },
    { MPOLS,    OPT_S, "simple-print", 0 },
    { NEGBIN,   OPT_G, "opg", 0 },
    { NEGBIN,   OPT_M, "model1", 0 },
    { NEGBIN,   OPT_R, "robust", 0 },
    { NEGBIN,   OPT_V, "verbose", 0 },
    { NLS,      OPT_N, "numerical", 0 },
    { NLS,      OPT_O, "vcv", 0 },
    { NLS,      OPT_R, "robust", 0 },
    { NLS,      OPT_V, "verbose", 0 },
    { NORMTEST, OPT_A, "all", 0 },
    { NORMTEST, OPT_D, "dhansen", 0 },
    { NORMTEST, OPT_W, "swilk", 0 },
    { NORMTEST, OPT_J, "jbera", 0 },
    { NORMTEST, OPT_L, "lillie", 0 },
    { NULLDATA, OPT_P, "preserve", 0 },
    { OLS,      OPT_F, "print-final", 0 },
    { OLS,      OPT_J, "jackknife", 0 },
    { OLS,      OPT_N, "no-df-corr", 0 },
    { OLS,      OPT_O, "vcv", 0 }, 
    { OLS,      OPT_R, "robust", 0 },
    { OLS,      OPT_Q, "quiet", 0 },
    { OLS,      OPT_S, "simple-print", 0 },
    { OLS,      OPT_V, "anova", 0 },
    { OMIT,     OPT_A, "auto", 1 },
    { OMIT,     OPT_B, "both", 0 },
    { OMIT,     OPT_X, "chi-square", 0 },
    { OMIT,     OPT_I, "silent", 0 },
    { OMIT,     OPT_W, "test-only", 0 },
    { OPEN,     OPT_B, "progress-bar", 0 },
    { OPEN,     OPT_C, "coded", 0 },
    { OPEN,     OPT_D, "drop-empty", 0 },
    { OPEN,     OPT_F, "cols", 2 },    
    { OPEN,     OPT_O, "odbc", 0 },
    { OPEN,     OPT_P, "preserve", 0 },
    { OPEN,     OPT_W, "www", 0 },
    { OUTFILE,  OPT_A, "append", 0 },
    { OUTFILE,  OPT_C, "close", 0 },
    { OUTFILE,  OPT_W, "write", 0 },
    { PANEL,    OPT_B, "between", 0 },
    { PANEL,    OPT_D, "time-dummies", 0 },
    { PANEL,    OPT_F, "fixed-effects", 0 },
    { PANEL,    OPT_H, "hausman-reg", 0 }, /* backward compatibility */
    { PANEL,    OPT_I, "iterate", 0 },
    { PANEL,    OPT_M, "matrix-diff", 0 },
    { PANEL,    OPT_O, "vcv", 0 },
    { PANEL,    OPT_P, "pooled", 0 },
    { PANEL,    OPT_R, "robust", 0 },
    { PANEL,    OPT_S, "silent", 0 },
    { PANEL,    OPT_U, "random-effects", 0 },
    { PANEL,    OPT_V, "verbose", 0 },
    { PANEL,    OPT_W, "unit-weights", 0 },
    { POISSON,  OPT_R, "robust", 0 },    
    { POISSON,  OPT_V, "verbose", 0 },
    { PCA,      OPT_C, "covariance", 0 },
    { PCA,      OPT_A, "save-all", 0 },
    { PCA,      OPT_O, "save", 0 },
    { PERGM,    OPT_O, "bartlett", 0 },
    { PERGM,    OPT_L, "log", 0 },
    { PERGM,    OPT_R, "radians", 0 },
    { PERGM,    OPT_D, "degrees", 0 },
    { PLOT,     OPT_O, "one-scale", 0 },
    { PLOT,     OPT_S, "time-series", 0 },
    { PLOT,     OPT_T, "tall", 0 },
    { PRINT,    OPT_O, "byobs", 0 },
    { PRINT,    OPT_N, "no-dates", 0 },
    { PROBIT,   OPT_P, "p-values", 0 },
    { PROBIT,   OPT_R, "robust", 0 },
    { PROBIT,   OPT_V, "verbose", 0 },
    { QQPLOT,   OPT_R, "raw", 0 },
    { QQPLOT,   OPT_Z, "z-scores", 0 },      
    { QUANTREG, OPT_I, "intervals", 1 },
    { QUANTREG, OPT_N, "no-df-corr", 0 },
    { QUANTREG, OPT_R, "robust", 0 },
    { QUIT,     OPT_X, "exit", 0 },
    { RESET,    OPT_C, "cubes-only", 0 },
    { RESET,    OPT_R, "squares-only", 0 },
    { RESTRICT, OPT_B, "bootstrap", 0 },
    { RESTRICT, OPT_F, "full", 0 },
    { RESTRICT, OPT_J, "jitter", 0 },
    { RESTRICT, OPT_V, "verbose", 0 },
    { RESTRICT, OPT_L, "lbfgs", 0 },
    { RESTRICT, OPT_N, "no-scaling", 0 },
    { RESTRICT, OPT_S, "silent", 0 },
    { RESTRICT, OPT_W, "wald", 0 },
    { RUNS,     OPT_D, "difference", 0 },
    { RUNS,     OPT_E, "equal", 0 },
    { SCATTERS, OPT_L, "with-lines", 0 },
    { SCATTERS, OPT_U, "output", 2 },
    { SET,      OPT_F, "from-file", 2 },
    { SET,      OPT_T, "to-file", 2 },
    { SETINFO,  OPT_C, "continuous", 0 },
    { SETINFO,  OPT_D, "discrete", 0 },
    { SETOBS,   OPT_C, "stacked-cross-section", 0 },
    { SETOBS,   OPT_P, "panel-vars", 0 },
    { SETOBS,   OPT_R, "restructure", 0 },
    { SETOBS,   OPT_S, "stacked-time-series", 0 },
    { SETOBS,   OPT_T, "time-series", 0 },
    { SETOBS,   OPT_X, "cross-section", 0 },
    { SETOBS,   OPT_N, "special-time-series", 0 },
    { SETOBS,   OPT_M, "labels", 2 },
    { SMPL,     OPT_B, "balanced", 0 },
    { SMPL,     OPT_C, "contiguous", 0 },
    { SMPL,     OPT_F, "full", 0 },
    { SMPL,     OPT_O, "dummy", 0 },
    { SMPL,     OPT_M, "no-missing", 0 },
    { SMPL,     OPT_N, "random", 0 },
    { SMPL,     OPT_P, "replace", 0 }, 
    { SMPL,     OPT_R, "restrict", 0 },
    { SPEARMAN, OPT_V, "verbose", 0 },
    { SQUARE,   OPT_O, "cross", 0 },
    { STORE,    OPT_C, "csv", 0 },
    { STORE,    OPT_D, "database", 0 },
    { STORE,    OPT_F, "overwrite", 0 },
    { STORE,    OPT_G, "dat", 0 },
    { STORE,    OPT_J, "jmulti", 0 },
    { STORE,    OPT_M, "gnu-octave", 0 },
    { STORE,    OPT_N, "no-header", 0 },
    { STORE,    OPT_R, "gnu-R", 0 },
    { STORE,    OPT_T, "traditional", 0 },
    { STORE,    OPT_Z, "gzipped", 0 },
    { STORE,    OPT_X, "omit-obs", 0 },
    { STORE,    OPT_E, "comment", 2 },
    { SUMMARY,  OPT_B, "by", 2 },
    { SUMMARY,  OPT_M, "matrix", 2 },
    { SUMMARY,  OPT_S, "simple", 0 },
    { SYSTEM,   OPT_I, "iterate", 0 },
    { SYSTEM,   OPT_V, "verbose", 0 },
#if 0
    { TABPRINT, OPT_F, "filename", 0 }, /* backward compatibility */
#endif
    { TOBIT,    OPT_L, "llimit", 2 },
    { TOBIT,    OPT_M, "rlimit", 2 },
    { TOBIT,    OPT_R, "robust", 0 },
    { TOBIT,    OPT_V, "verbose", 0 },
    { VAR,      OPT_D, "seasonals", 0 },
    { VAR,      OPT_F, "variance-decomp", 0 },
    { VAR,      OPT_H, "robust-hac", 0 },
    { VAR,      OPT_I, "impulse-responses", 0 },
    { VAR,      OPT_L, "lagselect", 0 },
    { VAR,      OPT_N, "nc", 0 },
    { VAR,      OPT_R, "robust", 0 }, 
    { VAR,      OPT_T, "trend", 0 }, 
    { VAR,      OPT_S, "lags", 2 },
    { VAR,      OPT_X, "silent", 0 },
    { VARLIST,  OPT_A, "accessors", 0 },
    { VARLIST,  OPT_S, "scalars", 0 },    
    { VECM,     OPT_A, "crt", 0 },
    { VECM,     OPT_D, "seasonals", 0 },
    { VECM,     OPT_F, "variance-decomp", 0 },
    { VECM,     OPT_I, "impulse-responses", 0 },
    { VECM,     OPT_N, "nc", 0 },
    { VECM,     OPT_R, "rc", 0 },
    { VECM,     OPT_T, "ct", 0 },
    { VECM,     OPT_V, "verbose", 0 },
    { VECM,     OPT_X, "silent", 0 },
    { WLS,      OPT_R, "robust", 0 },    
    { XTAB,     OPT_C, "column", 0 },
    { XTAB,     OPT_M, "matrix", 2 },
    { XTAB,     OPT_R, "row", 0 },
    { XTAB,     OPT_Z, "zeros", 0 },
    { 0,        0L,    NULL, 0 }
};

/* this function is used in compiling the gretl reference
   manual */

char **get_all_option_strings (int *pn)
{
    char **optstrs;
    int i, n = 0;

    for (i=0; gretl_opts[i].ci != 0; i++) {
	n++;
    }

    optstrs = strings_array_new(n);

    if (optstrs != NULL) {
	for (i=0; i<n; i++) {
	    optstrs[i] = gretl_strdup(gretl_opts[i].longopt);
	    if (optstrs[i] == NULL) {
		free_strings_array(optstrs, n);
		optstrs = NULL;
		break;
	    }
	}
    }

    if (optstrs != NULL) {
	strings_array_sort(&optstrs, &n, OPT_U);
	*pn = n;
    }
    
    return optstrs;
}

const char **get_opts_for_command (int ci, int *nopt)
{
    int i, j, n = 0;
    const char **ret = NULL;

    if (ci != OLS) {
	/* widely applicable options which are "attached" to OLS */
	n += vcv_opt_ok(ci);
	n += quiet_opt_ok(ci);
    }

    for (i=0; gretl_opts[i].ci != 0; i++) {
	n += (gretl_opts[i].ci == ci);
    }

    if (n == 0) {
	*nopt = 0;
	return NULL;
    }

    ret = malloc(n * sizeof *ret);
    if (ret == NULL) return NULL;

    j = 0;
    for (i=0; gretl_opts[i].ci != 0; i++) {
	if (gretl_opts[i].ci == ci) {
	    ret[j++] = gretl_opts[i].longopt;
	}
    }

    if (ci != OLS) {
	if (vcv_opt_ok(ci)) {
	    ret[j++] = "vcv";
	}
	if (quiet_opt_ok(ci)) {
	    ret[j++] = "quiet";
	}
    }

    *nopt = n;

    return ret;
}

struct flag_match flag_matches[] = {
    { OPT_A, 'a' },
    { OPT_B, 'b' },
    { OPT_C, 'c' },
    { OPT_D, 'd' },
    { OPT_E, 'e' },
    { OPT_F, 'f' },
    { OPT_G, 'g' },
    { OPT_H, 'h' },
    { OPT_I, 'i' },
    { OPT_J, 'j' },
    { OPT_K, 'k' },
    { OPT_L, 'l' },
    { OPT_M, 'm' },
    { OPT_N, 'n' },
    { OPT_O, 'o' },
    { OPT_P, 'p' },
    { OPT_Q, 'q' },
    { OPT_R, 'r' },
    { OPT_S, 's' },
    { OPT_T, 't' },
    { OPT_U, 'u' },
    { OPT_V, 'v' },
    { OPT_W, 'w' },
    { OPT_X, 'x' },
    { OPT_Y, 'y' },
    { OPT_Z, 'z' },
    { 0L,   '\0' }
};

static const char *ok_flags = "abcdefghijklmnopqrstuvwxyz";

#define isflag(c) (c && (strchr(ok_flags, c) != NULL))

/**
 * opt_from_flag:
 *
 * Returns: the gretl option value associated with a given
 * single-character flag.  Gives 0 if there is no associated
 * option.
 */

gretlopt opt_from_flag (unsigned char c)
{
    int i;

    for (i=0; flag_matches[i].c != '\0'; i++) {
	if (c == flag_matches[i].c) {
	    return flag_matches[i].o;
	}
    }

    return 0L;
}

static int opt_is_valid (gretlopt opt, int ci, char c)
{
    int i;

    if (opt == OPT_O && vcv_opt_ok(ci)) {
	return 1;
    }

    if (opt == OPT_Q && quiet_opt_ok(ci)) {
	return 1;
    }    

    for (i=0; gretl_opts[i].ci != 0; i++) {
	if (ci == gretl_opts[i].ci && opt == gretl_opts[i].o) {
	    return 1;
	}
    }

    if (c != 0) {
	gretl_errmsg_sprintf("Invalid option '-%c'", c);
    } 

    return 0;
}

gretlopt valid_short_opt (int ci, char c)
{
    gretlopt opt = 0;
    int i, ok;

    for (i=0; flag_matches[i].c != '\0'; i++) {
	if (c == flag_matches[i].c) {
	    opt = flag_matches[i].o;
	    break;
	}
    }

    if (opt) {
	ok = opt_is_valid(opt, ci, c);
	if (!ok) {
	    opt = 0;
	}
    }
	
    return opt;
}

/* See if at point @p (at which we've found '-') in string @s we
   might be at the start of an option flag: the previous character
   must be a space, and we must not be inside a quoted string.
*/

static int maybe_opt_start (char *s, char *p)
{
    int i, n = p - s;
    int quoted = 0;

    if (n > 0 && !isspace(*(p-1))) {
	return 0;
    }

    for (i=0; i<n; i++) {
	if (s[i] == '"' && (i == 0 || s[i-1] != '\\')) {
	    quoted = !quoted;
	}
    }

    return !quoted;
}

static gretlopt get_short_opts (char *line, int ci, int *err)
{
    char *s = line;
    gretlopt opt, ret = 0L;

    while ((s = strchr(s, '-')) != NULL) {
	char *p = s + 1;
	int i, n = 0;

	if (maybe_opt_start(line, s)) {
	    n = strspn(p, ok_flags);
	    if (n > 0) {
		if (isspace(p[n]) || p[n] == '\0') {
		    for (i=0; i<n; i++) {
			opt = opt_from_flag(p[i]);
			if (!opt_is_valid(opt, ci, p[i])) {
			    if (err != NULL) {
				*err = 1;
			    }
			    return 0L;
			}
			ret |= opt;
		    }
		} else {
		    n = 0;
		}
	    }
	}

	if (n > 0) {
	    gretl_delete(s, 0, n + 1);
	} else {
	    s++;
	}
    }

    return ret;
}

/* Apparatus for setting and retrieving parameters associated
   with command options, as in --opt=val.  
*/

typedef struct optparm_ optparm;

struct optparm_ {
    int ci;
    gretlopt opt;
    char *val;
};

static optparm *optparms;
static int n_parms;

/* Note: the following is called at the start of parse_command_line(),
   via gretl_cmd_clear().  Hopefully this ensures that we won't get an
   accumulation of stale data in the record of option parameters.
*/

/**
 * clear_option_params:
 *
 * Clears any ancillary parameter values currently associated 
 * with gretl command/option pairs.
 */

void clear_option_params (void)
{
    int i;

#if OPTDEBUG > 1
    fprintf(stderr, "clearing option params\n");
#endif

    for (i=0; i<n_parms; i++) {
	free(optparms[i].val);
    }
    free(optparms);
    optparms = NULL;
    n_parms = 0;
}

static optparm *matching_optparm (int ci, gretlopt opt)
{
    int i;

    for (i=0; i<n_parms; i++) {
	if (optparms[i].ci == ci && optparms[i].opt == opt) {
	    return &optparms[i];
	}
    }

    return NULL;
}

/* for a given @ci, @opt pair, determine its status with regard
   to a parameter value (not allowed, allowed, or required)
*/

static int option_parm_status (int ci, gretlopt opt)
{
    int i;

    for (i=0; gretl_opts[i].ci != 0; i++) {
	if (gretl_opts[i].ci == ci && gretl_opts[i].o == opt) {
	    return gretl_opts[i].parminfo;
	}
    }

    return 0;
}

static int real_push_option_param (int ci, gretlopt opt, char *val,
				   int checked)
{
    optparm *op;
    int n, err = 0;

    if (!checked && option_parm_status(ci, opt) == OPT_NO_PARM)  {
	return E_DATA;
    }

    n = n_parms + 1;

    op = matching_optparm(ci, opt);
    if (op != NULL) {
	/* got a match for the ci, opt pair already */
	free(op->val);
	op->val = val;
	return 0;
    }

    op = realloc(optparms, n * sizeof *op);

    if (op == NULL) {
	err = E_ALLOC;
    } else {
	optparms = op;
	op = &optparms[n-1];
	op->ci = ci;
	op->opt = opt;
	op->val = val;
	n_parms = n;
    }

#if OPTDEBUG 
    fprintf(stderr, "push_option_param: val='%s', n_parms=%d, "
	    "err = %d\n", op->val, n_parms, err);
#endif

    return err;
}

/**
 * push_option_param:
 * @ci: gretl command index.
 * @opt: gretl option flag.
 * @val: parameter value as string.
 *
 * Pushes onto an internal stack a record of the @val
 * to be associated with @opt for the current @ci.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int push_option_param (int ci, gretlopt opt, char *val)
{
    return real_push_option_param(ci, opt, val, 0);
}

/**
 * get_optval_string:
 * @ci: gretl command index.
 * @opt: gretl option value.
 *
 * Returns: the ancillary string value currently 
 * associated with option @opt for command @ci, if any,
 * otherwise %NULL.
 */

const char *get_optval_string (int ci, gretlopt opt)
{
    optparm *op = matching_optparm(ci, opt);

    return (op != NULL)? op->val : NULL;
}

/**
 * get_optval_double:
 * @ci: gretl command index.
 * @opt: gretl option value.
 *
 * Returns: the double-precision ancillary value currently 
 * associated with option @opt for command @ci, if any,
 * otherwise #NADBL.
 */

double get_optval_double (int ci, gretlopt opt)
{
    optparm *op = matching_optparm(ci, opt);
    double ret = NADBL;

    if (op != NULL && op->val != NULL) {
	if (numeric_string(op->val)) {
	    ret = dot_atof(op->val);
	} else {
	    ret = gretl_scalar_get_value(op->val);
	}
    }

    return ret;
}

/* called via GUI */

/**
 * set_optval_double:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @x: value to set.
 *
 * Sets a double-precision ancillary value to be associated
 * with option @opt for command @ci.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int set_optval_double (int ci, gretlopt opt, double x)
{
    char s[32];

    gretl_push_c_numeric_locale();
    sprintf(s, "%g", x);
    gretl_pop_c_numeric_locale();

    return real_push_option_param(ci, opt, gretl_strdup(s), 1);
}

gretlopt valid_long_opt (int ci, const char *lopt, OptStatus *status)
{
    int opt = OPT_NONE;
    int i;

    *status = 0;

    if (*lopt == '\0') {
	return 0;
    }

    if (vcv_opt_ok(ci) && !strcmp(lopt, "vcv")) {
	return OPT_O;
    }

    if (quiet_opt_ok(ci) && !strcmp(lopt, "quiet")) {
	return OPT_Q;
    }

    /* start by looking for an exact match */
    for (i=0; gretl_opts[i].o != 0; i++) {
	if (ci == gretl_opts[i].ci && 
	    !strcmp(lopt, gretl_opts[i].longopt)) {
	    opt = gretl_opts[i].o;
	    *status = gretl_opts[i].parminfo;
	    break;
	}
    }

    /* if this failed, try for an abbreviation or extension */
    if (opt == OPT_NONE) {
	int len, len1, len2;

	len1 = strlen(lopt);
	for (i=0; gretl_opts[i].o != 0; i++) {
	    len2 = strlen(gretl_opts[i].longopt);
	    len = (len2 > len1)? len1 : len2;
	    if (ci == gretl_opts[i].ci && 
		!strncmp(lopt, gretl_opts[i].longopt, len)) {
		opt = gretl_opts[i].o;
		*status = gretl_opts[i].parminfo;
		break;
	    }
	} 
    }  

    /* backward compatibility */
    if (opt == OPT_NONE && !strcmp(lopt, "wald")) {
	opt = OPT_W;
	*status = 0;
    }

    return opt;
}

#define data_open_special(s) (!strcmp(s, "sheet") || \
                              !strcmp(s, "coloffset") || \
			      !strcmp(s, "rowoffset"))

/* extract an option parameter value following '=' */

static int handle_optval (char *s, int ci, gretlopt opt, int status)
{
    char *p = s + 1; /* skip '=' */
    char *val = NULL;
    int quoted = 0;
    int n, err = 0;

    if (*p == '"') {
	/* handle a quoted value (e.g. a filename) */
	quoted = 1;
	p++;
	n = strcspn(p, "\"");
	if (n > 0 && *(p + n) == '"') {
	    val = gretl_strndup(p, n);
	} else {
	    err = E_PARSE;
	}
    } else if (*p != '\0') {
	/* plain unquoted value */
	n = strcspn(p, " =");
	if (n > 0) {
	    val = gretl_strndup(p, n);
	} else {
	    err = E_PARSE;
	}
    } else {
	err = E_PARSE;
    }

    if (val == NULL && !err) {
	/* allocation must have failed */
	err = E_ALLOC;
    } 

    if (!err) {
	if (status == OPT_NEEDS_PARM || status == OPT_ACCEPTS_PARM) {
	    err = real_push_option_param(ci, opt, val, 1);
	} else {
	    err = E_DATA;
	}
    }

    if (err) {
	free(val);
    } else {
	gretl_delete(s, 0, 1 + strlen(val) + 2 * quoted);
    }

    return err;
}

/* Crawl along @line looking for long-style option flags, possibly
   with associated parameter values; check options for validity
   in context.
*/
  
static gretlopt get_long_opts (char *line, int ci, int *err)
{
    char *s = line;
    char longopt[32];
    OptStatus status = 0;
    gretlopt match, ret = 0L;

    while ((s = strstr(s, "--")) != NULL) {
	match = 0;
	*longopt = '\0';

#if OPTDEBUG
	fprintf(stderr, "get_long_opts: s = '%s'\n", s);
#endif

	if (maybe_opt_start(line, s)) {
	    sscanf(s + 2, "%31[^ =]", longopt);
	    match = valid_long_opt(ci, longopt, &status);
	    if (match > 0) {
		/* recognized an acceptable option flag */
		ret |= match;
	    } else if ((ci == OPEN || ci == APPEND) && 
		       data_open_special(longopt)) {
		; /* no-op here: handled elsewhere (FIXME) */
	    } else {
		/* not a valid flag, or not applicable in context */
		gretl_errmsg_sprintf(_("Invalid option '--%s'"), longopt);
		fprintf(stderr, " line='%s', ci = %d\n", line, ci);
		*err = 1;
		return 0L;
	    } 
	}

	if (match > 0) {
	    gretl_delete(s, 0, 2 + strlen(longopt));
	    if (*s == '=') {
		/* got a param value: is it OK? */
		if (status == OPT_NO_PARM) {
		    *err = E_PARSE;
		} else {
		    *err = handle_optval(s, ci, match, status);
		}
	    } else if (status == OPT_NEEDS_PARM) {
		/* we need a param value but there's none */
		gretl_errmsg_sprintf(_("The option '--%s' requires a parameter"), 
				     longopt);
		*err = E_PARSE;
	    }
	    if (*err) {
		return 0L;
	    }	    
	} else {
	    s += 2;
	}
    }

    return ret;
}

/* Try to get the command word, skipping a possible assignment to a
   named object as in "<objectname> <- command", and bearing in mind
   that <objectname> may have embedded spaces (in which case it will
   be wrapped in quotes).
*/

static void get_cmdword (const char *line, char *word)
{
    if (*line == '"') {
	/* skip to closing quote */
	line++;
	line += strcspn(line, "\"");
	if (*line) line++;
	line += strspn(line, " ");
	if (!strncmp(line, "<-", 2)) {
	    line += 2;
	}
	sscanf(line, "%8s", word);
    } else if (!sscanf(line, "%*s <- %8s", word)) {
	sscanf(line, "%8s", word);
    }
}

static void tail_strip (char *s)
{
    int i, n = strlen(s);

    for (i=n-1; i>0; i--) {
	if (s[i] == ' ') {
	    s[i] = '\0';
	}
	else break;
    }
}

static int end_block_ci (const char *s)
{
    if (strlen(s) > 4) {
	char word[9];

	sscanf(s + 4, "%8s", word);
	if (!strcmp(word, "nls")) {
	    return NLS;
	} else if (!strcmp(word, "mle")) {
	    return MLE;
	} else if (!strcmp(word, "gmm")) {
	    return GMM;
	} else if (!strcmp(word, "restrict")) {
	    return RESTRICT;
	} else if (!strcmp(word, "kalman")) {
	    return KALMAN;
	} else if (!strcmp(word, "foreign")) {
	    return FOREIGN;
	} else if (!strcmp(word, "system")) {
	    return SYSTEM;
	}
    }

    return END;
}

#define ar1_alias(s) (!strcmp(s, "corc") || \
		      !strcmp(s, "hilu") || \
		      !strcmp(s, "pwe"))

#define ols_alias(s) (!strcmp(s, "hccm"))

#define smpl_alias(s) (!strcmp(s, "sample"))

#define modtest_alias(s) (!strcmp(s, "lmtest"))

/**
 * get_gretl_options:
 * @line: command line to parse.
 * @err: location for error code, which is set to 1 in case any 
 * invalid options are found, else set to 0.
 * 
 * Check for option flags in @line: if found, chop them out and set
 * the return value accordingly. 
 *
 * Returns: the options found in the line.
 */

gretlopt get_gretl_options (char *line, int *err)
{
    gretlopt oflags = 0L;
    gretlopt opt;
    char cmdword[9] = {0};
    int endblock = 0;
    int ci, myerr = 0;

    if (err != NULL) {
	*err = 0;
    }

    if (strlen(line) < 2 || *line == '#') {
	return oflags;
    }

    get_cmdword(line, cmdword);

    if (!strcmp(cmdword, "end")) {
	endblock = 1;
	ci = end_block_ci(line);
    } else if (ar1_alias(cmdword)) {
	ci = AR1;
    } else if (smpl_alias(cmdword)) {
	ci = SMPL;
    } else if (ols_alias(cmdword)) {
	ci = OLS;
    } else if (modtest_alias(cmdword)) {
	ci = MODTEST;
    } else {
	ci = gretl_command_number(cmdword);
    }

#if OPTDEBUG
    fprintf(stderr, "get gretl_options: '%s' -> ci = %d\n", cmdword, ci);
#endif

    /* some commands do not take a "flag", and "-%c" may have
       some other meaning */
    if (ci == 0 || ci == GENR || ci == PRINTF) {
	return oflags;
    } else if (!endblock && (ci == GMM || ci == MLE || ci == NLS)) {
	return oflags;
    }

    /* smpl: in some contexts options don't make sense */
    if (ci == SMPL && strchr(line, ';')) {
	return oflags;
    }

    if (ci != SET && ci != SETINFO && ci != SMPL &&
	ci != TABPRINT && ci != EQNPRINT) {
	/* try for short-form options (e.g. "-o"): but note that
	   with some commands there's the possibility of collision
	   between short options and other syntactical elements,
	   so we only recognize long-form options.
	*/
	opt = get_short_opts(line, ci, (ci == SMPL)? NULL : &myerr);
	if (!myerr && opt) {
	    oflags |= opt;
	}
    }

    /* try for long-form options (e.g. "--vcv") */
    if (!myerr) {
	opt = get_long_opts(line, ci, &myerr);
	if (!myerr && opt) {
	    oflags |= opt;
	}
    }

    /* strip trailing whitespace after processing */
    tail_strip(line);

    if (err != NULL) {
	*err = myerr;
    }

    return oflags;
}

static PRN *flagprn;

/**
 * print_flags:
 * @oflags: options.
 * @ci: command index, for context.
 * 
 * Returns: a string representation of the options in @oflags,
 * or an empty string if no options are found.  The returned
 * value should not be modified in any way.
 */

const char *print_flags (gretlopt oflags, int ci)
{
    const char *parm;
    gretlopt opt;
    int i;

    if (ci == FUNDEBUG) {
	/* options are for internal use only */
	return "";
    }

    if (flagprn == NULL) {
	int err = 0;

	flagprn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
	if (err) {
	    return "";
	}
    } else {
	gretl_print_reset_buffer(flagprn);
    }

    if (oflags == OPT_NONE || ci == QUIT || ci == GENR) {
	/* no options, or only hidden ones */
	return "";
    }

    /* special: -o (--vcv) can be used with several model
       commands */
    if ((oflags & OPT_O) && vcv_opt_ok(ci)) {
	pputs(flagprn, " --vcv");
	oflags &= ~OPT_O; /* handled */
    }

    if ((oflags & OPT_Q) && quiet_opt_ok(ci)) {
	pputs(flagprn, " --quiet");
	oflags &= ~OPT_Q; /* handled */
    }

    for (i=0; gretl_opts[i].ci != 0; i++) {
	opt = gretl_opts[i].o;
	if (ci == gretl_opts[i].ci && (oflags & opt)) {
	    pprintf(flagprn, " --%s", gretl_opts[i].longopt);
	    if (gretl_opts[i].parminfo) {
		parm = get_optval_string(ci, opt);
		if (parm != NULL && *parm != '\0') {
		    pprintf(flagprn, "=%s", parm);
		}
	    }
	}
    }

    return gretl_print_get_buffer(flagprn);
}

void option_flags_cleanup (void)
{
    gretl_print_destroy(flagprn);
    flagprn = NULL;
}

/**
 * check_for_loop_only_options:
 * @ci: gretl command index.
 * @opt: option flag to be tested.
 * @prn: gretl printing struct.
 * 
 * Returns: 1 if option @opt is applicable for command @ci only
 * in the context of a command loop (in which case a warning is
 * printed to @prn), otherwise 0.  
 */

int check_for_loop_only_options (int ci, gretlopt opt, PRN *prn)
{
    int ret = 0;

    if (ci == OLS && (opt & OPT_P)) {
	const char *flagstr = print_flags(OPT_P, OLS);

	pprintf(prn, _("Warning: option%s ignored outside of loop"), 
		flagstr);
	pputc(prn, '\n');
	ret = 1;
    }

    return ret;
}

/**
 * transcribe_options:
 * @targ: pointer to target option flags.
 * @src: source option flags.
 * @test: bitwise OR of flags that are of interest in context.
 *
 * If the intersection of the flags in @src and @test is non-
 * empty, set the corresponding flags in @targ.
 * 
 * Returns: the (possibly modified) @targ.
 */

gretlopt transcribe_option_flags (gretlopt *targ, gretlopt src,
				  gretlopt test)
{
    *targ |= (src & test);

    return *targ;
}

/**
 * delete_option_flags:
 * @targ: pointer to target option flags.
 * @test: bitwise OR of flags that to be deleted.
 *
 * If the intersection of the flags in @targ and @test is non-
 * empty, unset the corresponding flags in @targ.
 * 
 * Returns: the (possibly modified) @targ.
 */

gretlopt delete_option_flags (gretlopt *targ, gretlopt test)
{
    *targ &= ~(*targ & test);

    return *targ;
}

/**
 * incompatible_options:
 * @opt: option flags to be tested.
 * @test: bitwise OR of flags that are incompatible in context.
 * 
 * Returns: %E_BADOPT if @opt contains more than one of the flags
 * in @test, otherwise 0.
 */

int incompatible_options (gretlopt opt, gretlopt test)
{
    int optcount = 0;
    gretlopt o;

    for (o=OPT_A; o<=OPT_Z; o=o<<1) {
	if ((opt & o) && (test & o)) {
	    optcount++;
	    if (optcount > 1) {
		return E_BADOPT;
	    }
	}
    }

    return 0;
}

/**
 * option_prereq_missing:
 * @opt: option flags to be tested.
 * @test: bitwise OR of flags that have a deinite prequisite.
 * @prereq: bitwise OR of prequisite flags.
 * 
 * Returns: %E_BADOPT if @opt contains at least one element of
 * @test but no elements of @prereq, otherwise 0.
 */

int option_prereq_missing (gretlopt opt, gretlopt test,
			   gretlopt prereq)
{
    if (opt & test) {
	if (!(opt & prereq)) {
	    return E_BADOPT;
	}
    }

    return 0;
}

/**
 * inapplicable_option_error:
 * @ci: command index.
 * @opt: bad option flag.
 *
 * Flags an error: to be used when @opt is not applicable in the
 * context of command @ci, in context.
 *
 * Returns: %E_BADOPT.
 */

int inapplicable_option_error (int ci, gretlopt opt)
{
    const char *s = print_flags(opt, ci);

    gretl_errmsg_sprintf("%s: inapplicable option", s);
    return E_BADOPT;
}
