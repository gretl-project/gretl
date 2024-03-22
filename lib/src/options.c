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
#include "uservar.h"
#include "gretl_func.h"

#include <errno.h>

#define OPTDEBUG 0

/* commands for which --vcv (OPT_O) is applicable */
#define vcv_opt_ok(c) (MODEL_COMMAND(c) || c == ADD || c == OMIT)

/* commands for which --window (OPT_W) is applicable */
#define window_opt_ok(c) (MODEL_COMMAND(c) || c == VAR || c == VECM)

/* commands which support the $result accessor */
#define yields_result(c) (c == CORR || c == FREQ || c == SUMMARY)

/* commands for which --quiet (OPT_Q) is applicable */
#define quiet_opt_ok(c) (MODEL_COMMAND(c) ||    \
                         yields_result(c) ||    \
                         c == ADD ||            \
                         c == ADF ||            \
                         c == ANOVA ||          \
                         c == APPEND ||         \
			 c == BDS ||		\
                         c == BKW ||            \
                         c == COEFFSUM ||       \
                         c == CHOW ||           \
                         c == COINT2 ||         \
                         c == CORRGM ||         \
                         c == CUSUM ||          \
                         c == DATA ||           \
                         c == DIFFTEST ||       \
                         c == ESTIMATE ||       \
                         c == FCAST ||          \
                         c == FOREIGN ||        \
                         c == FRACTINT ||       \
                         c == KPSS ||           \
                         c == MAKEPKG ||        \
                         c == MODTEST ||        \
                         c == LEVERAGE ||       \
                         c == LEVINLIN ||       \
                         c == LOOP ||           \
                         c == MAHAL ||          \
                         c == NORMTEST ||       \
                         c == OLS ||            \
                         c == OMIT ||           \
                         c == OPEN ||           \
			 c == PANSPEC ||	\
                         c == PKG ||            \
                         c == QLRTEST ||        \
                         c == RENAME ||         \
                         c == RESET ||          \
                         c == RESTRICT ||       \
                         c == RMPLOT ||         \
                         c == SMPL ||           \
                         c == SYSTEM ||         \
                         c == VAR ||            \
                         c == VECM ||           \
                         c == VIF ||            \
                         c == XCORRGM ||        \
                         c == XTAB)

/* --output (OPT_U) as attached to GNUPLOT */
#define plot_output_opt_ok(c) (c == GNUPLOT ||	\
			       c == PLOT ||	\
			       c == BXPLOT ||	\
			       c == HFPLOT ||	\
			       c == PANPLOT ||	\
			       c == QQPLOT ||	\
			       c == KDPLOT ||   \
			       c == RMPLOT ||	\
			       c == SCATTERS || \
			       c == GRIDPLOT)

/* --plot (OPT_U) as attached to CORR */
#define cmd_plot_opt_ok(c) (c == CORR ||	\
			    c == CORRGM ||	\
			    c == CUSUM ||	\
			    c == FCAST ||	\
			    c == FREQ ||	\
			    c == HURST ||	\
			    c == LEVERAGE ||	\
			    c == PERGM ||	\
			    c == QLRTEST ||	\
			    c == XCORRGM)

/* --outbuf (OPT_b) as attached to GNUPLOT */
#define plot_outbuf_opt_ok(c) (plot_output_opt_ok(c) || \
			       cmd_plot_opt_ok(c))

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

/* Below: This is used as a one-way mapping from the long form to the
   index (e.g. OPT_Q), so a given index can have more than one
   long-form counterpart depending on the context.  The last field
   indicates whether the given option does not accept (0), accepts
   (1), or requires (2) an associated parameter value.
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
    { ADF,      OPT_V, "verbose", 1 },
    { ADF,      OPT_F, "difference", 0 },
    { ADF,      OPT_E, "test-down", 1 },
    { ADF,      OPT_G, "gls", 0 },
    { ADF,      OPT_U, "perron-qu", 0 },
    { AR1,      OPT_B, "no-corc", 0 },
    { AR1,      OPT_H, "hilu", 0 },
    { AR1,      OPT_P, "pwe", 0 },
    { AR1,      OPT_L, "loose", 0 },
    { APPEND,   OPT_A, "all-cols", 0 },
    { APPEND,   OPT_T, "time-series", 0 },
    { APPEND,   OPT_R, "rowoffset", 2 },
    { APPEND,   OPT_C, "coloffset", 2 },
    { APPEND,   OPT_F, "fixed-cols", 2 },
    { APPEND,   OPT_M, "rowmask", 2 },
    { APPEND,   OPT_L, "cols", 2 },
    { APPEND,   OPT_S, "sheet", 2 },
    { APPEND,   OPT_V, "verbose", 0 },
    { APPEND,   OPT_U, "update-overlap", 0 },
    { APPEND,   OPT_X, "fixed-sample", 0 },
    { APPEND,   OPT_K, "frompkg", 2 },
    { ARMA,     OPT_A, "as154", 0 },
    { ARMA,     OPT_C, "conditional", 0 },
    { ARMA,     OPT_E, "save-ehat", 0 },
    { ARMA,     OPT_G, "opg", 0 },
    { ARMA,     OPT_H, "hessian", 0 },
    { ARMA,     OPT_K, "kalman", 0 },
    { ARMA,     OPT_L, "lbfgs", 0 },
    { ARMA,     OPT_N, "nc", 0 },
    { ARMA,     OPT_V, "verbose", 0 },
    { ARMA,     OPT_X, "x-13arima", 0 },
    { ARMA,     OPT_X, "x-12-arima", 0 }, /* compatibility alias */
    { ARMA,     OPT_Y, "y-diff-only", 0 },
    { ARMA,     OPT_R, "robust", 0 },
    { ARMA,     OPT_S, "stdx", 0 },
    { ARMA,     OPT_Z, "lagselect", 0 },
    { BDS,      OPT_B, "boot", 2 },
    { BDS,      OPT_C, "corr1", 2 },
    { BDS,      OPT_S, "sdcrit", 2 },
    { BDS,      OPT_X, "matrix", 2 },
    { BIPROBIT, OPT_G, "opg", 0 },
    { BIPROBIT, OPT_R, "robust", 0 },
    { BIPROBIT, OPT_V, "verbose", 0 },
    { BIPROBIT, OPT_C, "cluster", 2 },
    { BIPROBIT, OPT_X, "save-xbeta", 0 },
    { BXPLOT,   OPT_O, "notches", 0 },
    { BXPLOT,   OPT_K, "tweaks", 2 },
    { BXPLOT,   OPT_L, "outliers", 2 },
    { BXPLOT,   OPT_P, "panel", 0 },
    { BXPLOT,   OPT_X, "matrix", 2 },
    { BXPLOT,   OPT_Z, "factorized", 0 },
    { BXPLOT,   OPT_B, "whiskerbars", 0 },
    { CHOW,     OPT_D, "dummy", 0 },
    { CHOW,     OPT_L, "limit-to", 2 },
    { CLEAR,    OPT_D, "dataset", 0 },
    { CLEAR,    OPT_F, "functions", 0 },
    { COINT,    OPT_D, "seasonals", 0 },
    { COINT,    OPT_E, "test-down", 1 },
    { COINT,    OPT_N, "nc", 0 },
    { COINT,    OPT_R, "ctt", 0 },
    { COINT,    OPT_S, "skip-df", 0 },
    { COINT,    OPT_T, "ct", 0 },
    { COINT,    OPT_V, "verbose", 0 },
    { COINT,    OPT_I, "silent", 0 },
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
    { CORR,     OPT_N, "uniform", 0 },
    { CORR,     OPT_V, "verbose", 0 },
    { CORR,     OPT_U, "plot", 2 },
    { CORR,     OPT_X, "matrix", 2 },
    { CORR,     OPT_T, "triangle", 0 },
    { CORRGM,   OPT_A, "acf-only", 0 },
    { CORRGM,   OPT_B, "bartlett", 0 },
    { CORRGM,   OPT_S, "silent", 0 },
    { CUSUM,    OPT_R, "squares", 0 },
    { DATA,     OPT_C, "compact", 2 },
    { DATA,     OPT_O, "odbc", 0 },
    { DATA,     OPT_N, "name", 2 },
    { DATA,     OPT_V, "verbose", 0 },
    { DATA,     OPT_F, "no-align", 0 },
    { DATAMOD,  OPT_P, "preserve", 0 },
    { DATAMOD,  OPT_T, "panel-time", 0 },
    { DATAMOD,  OPT_W, "weekstart", 2 },
    { DATAMOD,  OPT_R, "repday", 2 },
    { DELEET,   OPT_D, "db", 0 },
    { DELEET,   OPT_F, "force", 0 },
    { DELEET,   OPT_L, "list", 0 },
    { DELEET,   OPT_T, "type", 2 },
    { DIFFTEST, OPT_G, "sign", 0 },
    { DIFFTEST, OPT_R, "rank-sum", 0 },
    { DIFFTEST, OPT_I, "signed-rank", 0 },
    { DIFFTEST, OPT_V, "verbose", 0 },
    { DISCRETE, OPT_R, "reverse", 0 },
    { DPANEL,   OPT_A, "asymptotic", 0 },
    { DPANEL,   OPT_D, "time-dummies", 1 },
    { DPANEL,   OPT_K, "keep-extra", 0 },
    { DPANEL,   OPT_L, "system", 0 },
    { DPANEL,   OPT_T, "two-step", 0 },
    { DPANEL,   OPT_V, "verbose", 0 },
    { DPANEL,   OPT_X, "dpdstyle", 0 },
    { DPANEL,   OPT_C, "collapse", 0 },
    { DUMMIFY,  OPT_F, "drop-first", 0 },
    { DUMMIFY,  OPT_L, "drop-last", 0 },
    { DURATION, OPT_B, "weibull", 0 },
    { DURATION, OPT_E, "exponential", 0 },
    { DURATION, OPT_L, "loglogistic", 0 },
    { DURATION, OPT_Z, "lognormal", 0 },
    { DURATION, OPT_M, "medians", 0 },
    { DURATION, OPT_G, "opg", 0 },
    { DURATION, OPT_R, "robust", 0 },
    { DURATION, OPT_C, "cluster", 2 },
    { DURATION, OPT_V, "verbose", 0 },
    { EQNPRINT, OPT_O, "complete", 0 },
    { EQNPRINT, OPT_T, "t-ratios", 0 },
    { EQNPRINT, OPT_U, "output", 2 },
    { TABPRINT, OPT_O, "complete", 0 },
    { TABPRINT, OPT_C, "csv", 0 },
    { TABPRINT, OPT_R, "rtf", 0 },
    { TABPRINT, OPT_T, "format", 2 },
    { TABPRINT, OPT_U, "output", 2 },
    { EQUATION, OPT_M, "multi", 0 },
    { ESTIMATE, OPT_I, "iterate", 0 },
    { ESTIMATE, OPT_M, "geomean", 0 },
    { ESTIMATE, OPT_N, "no-df-corr", 0 },
    { ESTIMATE, OPT_U, "unrestrict-init", 0 },
    { ESTIMATE, OPT_V, "verbose", 0 },
    { ESTIMATE, OPT_W, "window", 0 },
    { FCAST,    OPT_L, "all-probs", 0 },
    { FCAST,    OPT_D, "dynamic", 0 },
    { FCAST,    OPT_M, "mean-y", 0 },
    { FCAST,    OPT_N, "no-stats", 0 },
    { FCAST,    OPT_T, "stats-only", 0 },
    { FCAST,    OPT_S, "static", 0 },
    { FCAST,    OPT_R, "recursive", 0 },
    { FCAST,    OPT_R, "rolling", 0 }, /* legacy alias */
    { FCAST,    OPT_O, "out-of-sample", 0 },
    { FCAST,    OPT_I, "integrate", 0 },
    { FOREIGN,  OPT_D, "send-data", 1 },
    { FOREIGN,  OPT_V, "verbose", 0 },
    { FOREIGN,  OPT_F, "frame", 0 },
    { FOREIGN,  OPT_N, "no-compile", 0 },
    { FOREIGN,  OPT_I, "io-funcs", 2 },
    { FRACTINT, OPT_G, "gph", 0 },
    { FRACTINT, OPT_A, "all", 0 },
    { FREQ,     OPT_Q, "quiet", 0 }, /* legacy alias */
    { FREQ,     OPT_O, "gamma", 0 },
    { FREQ,     OPT_S, "silent", 0 },
    { FREQ,     OPT_Z, "normal", 0 },
    { FREQ,     OPT_N, "nbins", 2 },
    { FREQ,     OPT_M, "min", 2 },
    { FREQ,     OPT_W, "binwidth", 2 },
    { FREQ,     OPT_X, "matrix", 2 },
    { FREQ,     OPT_K, "tweaks", 2 },
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
    { GNUPLOT,  OPT_I, "input", 2 },
    { GNUPLOT,  OPT_O, "with-lines", 1 },
    { GNUPLOT,  OPT_F, "fit", 2 },
    { GNUPLOT,  OPT_K, "tweaks", 2 },
    { GNUPLOT,  OPT_M, "with-impulses", 1 },
    { GNUPLOT,  OPT_P, "with-lp", 1 },
    { GNUPLOT,  OPT_B, "with-boxes", 1 },
    { GNUPLOT,  OPT_Q, "with-steps", 1 },
    { GNUPLOT,  OPT_T, "time-series", 0 },
    { GNUPLOT,  OPT_Z, "dummy", 0 },
    { GNUPLOT,  OPT_C, "control", 0 },
    { GNUPLOT,  OPT_U, "output", 2 },
    { GNUPLOT,  OPT_b, "outbuf", 2 },
    { GNUPLOT,  OPT_b, "buffer", 2 }, /* compatibility alias */
    { GNUPLOT,  OPT_i, "inbuf", 2 },
    { GNUPLOT,  OPT_Y, "single-yaxis", 0 },
    { GNUPLOT,  OPT_X, "matrix", 2 },
    { GNUPLOT,  OPT_N, "band", 2 },
    { GNUPLOT,  OPT_J, "band-style", 2 },
    { GNUPLOT,  OPT_a, "bands", 2 },
    { GNUPLOT,  OPT_W, "font", 2 },
    { GNUPLOT,  OPT_L, "ylogscale", 1 },
    { GNUPLOT,  OPT_D, "y2axis", 2 },
    { GRAPHPG,  OPT_M, "monochrome", 0 },
    { GRAPHPG,  OPT_O, "output", 2 },
    { GRIDPLOT, OPT_F, "fontsize", 2 },
    { GRIDPLOT, OPT_W, "width", 2 },
    { GRIDPLOT, OPT_H, "height", 2 },
    { GRIDPLOT, OPT_R, "rows", 2 },
    { GRIDPLOT, OPT_C, "cols", 2 },
    { GRIDPLOT, OPT_L, "layout", 2 },
    { GRIDPLOT, OPT_U, "output", 2 },
    { HECKIT,   OPT_M, "ml", 0 },
    { HECKIT,   OPT_G, "opg", 0 },
    { HECKIT,   OPT_R, "robust", 0 },
    { HECKIT,   OPT_C, "cluster", 2 },
    { HECKIT,   OPT_T, "two-step", 0 },
    { HECKIT,   OPT_V, "verbose", 0 },
    { HELP,     OPT_F, "func", 0 },
    { HFPLOT,   OPT_O, "with-lines", 0 },
    { HFPLOT,   OPT_T, "time-series", 0 },
    { HSK,      OPT_N, "no-squares", 0 },
    { INCLUDE,  OPT_F, "force", 0 },
    { INFO,     OPT_F, "from-file", 2},
    { INFO,     OPT_T, "to-file", 2},
    { INTREG,   OPT_G, "opg", 0 },
    { INTREG,   OPT_R, "robust", 0 },
    { INTREG,   OPT_C, "cluster", 2 },
    { INTREG,   OPT_V, "verbose", 0 },
    { IVREG,    OPT_G, "gmm", 0 },
    { IVREG,    OPT_I, "iterate", 0 },
    { IVREG,    OPT_L, "liml", 0 },
    { IVREG,    OPT_N, "no-df-corr", 0 },
    { IVREG,    OPT_R, "robust", 0 },
    { IVREG,    OPT_S, "save", 0 },
    { IVREG,    OPT_T, "two-step", 0 },
    { IVREG,    OPT_H, "weights", 2 },
    { IVREG,    OPT_X, "no-tests", 0 },
    { IVREG,    OPT_C, "cluster", 2 },
    { IVREG,    OPT_M, "matrix-diff", 0 },
    { JOIN,     OPT_I, "ikey", 2 },
    { JOIN,     OPT_O, "okey", 2 },
    { JOIN,     OPT_F, "filter", 2 },
    { JOIN,     OPT_A, "aggr", 2 },
    { JOIN,     OPT_D, "data", 2 },
    { JOIN,     OPT_T, "tconv-fmt", 2 },
    { JOIN,     OPT_K, "tkey", 2 },
    { JOIN,     OPT_X, "tconvert", 2 },
    { JOIN,     OPT_H, "no-header", 0 },
    { JOIN,     OPT_V, "verbose", 0 },
    { JOIN,     OPT_P, "pd", 2 }, /* undocumented: is it wanted? */
    { JOIN,     OPT_R, "frompkg", 2 },
    { KDPLOT,   OPT_O, "alt", 0 },
    { KDPLOT,   OPT_S, "scale", 2},
    { KPSS,     OPT_T, "trend", 0 },
    { KPSS,     OPT_D, "seasonals", 0 },
    { KPSS,     OPT_V, "verbose", 0 },
    { KPSS,     OPT_F, "difference", 0 },
    { LAGS,     OPT_L, "bylag", 0 },
    { LEVERAGE, OPT_S, "save", 0 },
    { LEVERAGE, OPT_O, "overwrite", 0 },
    { LEVINLIN, OPT_N, "nc", 0 },
    { LEVINLIN, OPT_T, "ct", 0 },
    { LEVINLIN, OPT_V, "verbose", 0 },
    { MAKEPKG,  OPT_D, "dtd", 2 },
    { MAKEPKG,  OPT_I, "index", 0 },
    { MAKEPKG,  OPT_T, "translations", 0 },
    { MARKERS,  OPT_D, "delete", 0 },
    { MARKERS,  OPT_F, "from-file", 2 },
    { MARKERS,  OPT_T, "to-file", 2 },
    { MARKERS,  OPT_A, "to-array", 2 },
    { MARKERS,  OPT_R, "from-array", 2 },
    { MARKERS,  OPT_S, "from-series", 2 },
    { MODTEST,  OPT_A, "autocorr", 0 },
    { MODTEST,  OPT_B, "breusch-pagan", 0 },
    { MODTEST,  OPT_C, "comfac", 0 },
    { MODTEST,  OPT_D, "xdepend", 0 },
    { MODTEST,  OPT_H, "arch", 0 },
    { MODTEST,  OPT_L, "logs", 0 },
    { MODTEST,  OPT_N, "normality", 0 },
    { MODTEST,  OPT_S, "squares", 0 },
    { MODTEST,  OPT_P, "panel", 0 },
    { MODTEST,  OPT_R, "robust", 0 },
    { MODTEST,  OPT_W, "white", 0 },
    { MODTEST,  OPT_X, "white-nocross", 0 },
    { MODTEST,  OPT_I, "silent", 0 },
    { MODTEST,  OPT_U, "univariate", 0 },
    { MPI,      OPT_D, "send-data", 1 },
    { MPI,      OPT_M, "send-metadata", 0 },
    { MPI,      OPT_R, "full-range", 0 },
    { MPI,      OPT_F, "send-functions", 0 },
    { MPI,      OPT_L, "local", 0 },
    { MPI,      OPT_N, "np", 2 },
    { MPI,      OPT_T, "omp-threads", 2 },
    { MPI,      OPT_Q, "quiet", 0},
    { MPI,      OPT_V, "verbose", 1 },
    { MPI,      OPT_S, "single-rng", 0},
    { LABELS,   OPT_D, "delete", 0 },
    { LABELS,   OPT_F, "from-file", 2 },
    { LABELS,   OPT_T, "to-file", 2 },
    { LABELS,   OPT_A, "from-array", 2 },
    { LABELS,   OPT_R, "to-array", 2 },
    { LAD,      OPT_N, "no-vcv", 0 },
    { LOGISTIC, OPT_M, "ymax", 2 },
    { LOGISTIC, OPT_R, "robust", 0 },
    { LOGISTIC, OPT_C, "cluster", 2 },
    { LOGISTIC, OPT_F, "fixed-effects", 0 },
    { LOGIT,    OPT_M, "multinomial", 0 },
    { LOGIT,    OPT_P, "p-values", 0 },
    { LOGIT,    OPT_R, "robust", 0 },
    { LOGIT,    OPT_C, "cluster", 2 },
    { LOGIT,    OPT_V, "verbose", 0 },
    { LOGIT,    OPT_S, "estrella", 0 },
    { LOOP,     OPT_D, "decr", 0 },
    { LOOP,     OPT_P, "progressive", 0 },
    { LOOP,     OPT_V, "verbose", 0 },
    { MAHAL,    OPT_S, "save", 0 },
    { MAHAL,    OPT_V, "vcv", 0 },
    { MEANTEST, OPT_O, "unequal-vars", 0 },
    { MIDASREG, OPT_R, "robust", 0 },
    { MIDASREG, OPT_V, "verbose", 0 },
    { MIDASREG, OPT_L, "levenberg", 0 },
    { MIDASREG, OPT_P, "print-spec", 0 },
    { MIDASREG, OPT_B, "breaktest", 0 },
    { MIDASREG, OPT_C, "clamp-beta", 0 }, /* legacy */
    { MLE,      OPT_A, "auxiliary", 0 },
    { MLE,      OPT_G, "opg", 0 },
    { MLE,      OPT_H, "hessian", 0 },
    { MLE,      OPT_S, "no-gradient-check", 0 },
    { MLE,      OPT_L, "lbfgs", 0 },
    { MLE,      OPT_N, "numerical", 0 },
    { MLE,      OPT_R, "robust", 1 },
    { MLE,      OPT_C, "cluster", 2 },
    { MLE,      OPT_V, "verbose", 0 },
    { MODPRINT, OPT_A, "addstats", 2 },
    { MODPRINT, OPT_O, "output", 2 },
    { MODPRINT, OPT_C, "complete", 0 },
    { MODELTAB, OPT_O, "output", 2 },
    { MODELTAB, OPT_C, "complete", 0 },
    { MODELTAB, OPT_B, "options", 2 },
    { MPOLS,    OPT_S, "simple-print", 0 },
    { NEGBIN,   OPT_G, "opg", 0 },
    { NEGBIN,   OPT_M, "model1", 0 },
    { NEGBIN,   OPT_R, "robust", 0 },
    { NEGBIN,   OPT_C, "cluster", 2 },
    { NEGBIN,   OPT_V, "verbose", 0 },
    { NLS,      OPT_N, "numerical", 0 },
    { NLS,      OPT_R, "robust", 0 },
    { NLS,      OPT_V, "verbose", 0 },
    { NLS,      OPT_S, "no-gradient-check", 0 },
    { NORMTEST, OPT_A, "all", 0 },
    { NORMTEST, OPT_D, "dhansen", 0 },
    { NORMTEST, OPT_W, "swilk", 0 },
    { NORMTEST, OPT_J, "jbera", 0 },
    { NORMTEST, OPT_L, "lillie", 0 },
    { NULLDATA, OPT_N, "no-index", 0 },
    { NULLDATA, OPT_P, "preserve", 0 },
    { OLS,      OPT_F, "print-final", 0 },
    { OLS,      OPT_J, "jackknife", 0 },
    { OLS,      OPT_N, "no-df-corr", 0 },
    { OLS,      OPT_O, "vcv", 0 },
    { OLS,      OPT_R, "robust", 0 },
    { OLS,      OPT_Q, "quiet", 0 }, /* note: for the sake of documentation */
    { OLS,      OPT_S, "simple-print", 0 },
    { OLS,      OPT_V, "anova", 0 },
    { OLS,      OPT_C, "cluster", 2 },
    { OLS,      OPT_W, "window", 0 },
    { OMIT,     OPT_A, "auto", 1 },
    { OMIT,     OPT_B, "both", 0 },
    { OMIT,     OPT_X, "chi-square", 0 },
    { OMIT,     OPT_I, "silent", 0 },
    { OMIT,     OPT_W, "test-only", 0 },
    { OMIT,     OPT_T, "trend", 0 },      /* omit auto-trend: VAR only */
    { OMIT,     OPT_E, "seasonals", 0 },  /* omit auto-seasonals: VAR only */
    { OPEN,     OPT_A, "all-cols", 0 },
    { OPEN,     OPT_B, "progress-bar", 0 },
    { OPEN,     OPT_D, "drop-empty", 0 },
    { OPEN,     OPT_E, "select", 2 },
    { OPEN,     OPT_F, "fixed-cols", 2 },
    { OPEN,     OPT_O, "odbc", 0 },
    { OPEN,     OPT_P, "preserve", 0 },
    { OPEN,     OPT_R, "rowoffset", 2 },
    { OPEN,     OPT_C, "coloffset", 2 },
    { OPEN,     OPT_S, "sheet", 2 },
    { OPEN,     OPT_W, "www", 0 },
    { OPEN,     OPT_L, "cols", 2 },
    { OPEN,     OPT_M, "rowmask", 2 },
    { OPEN,     OPT_V, "verbose", 0 },
    { OPEN,     OPT_K, "frompkg", 2 },
    { OPEN,     OPT_H, "no-header", 0 },
    { OPEN,     OPT_I, "ignore-quotes", 0 },
    { OPEN,     OPT_U, "bundle", 2 },
    { OUTFILE,  OPT_A, "append", 0 },
    { OUTFILE,  OPT_Q, "quiet", 0 },
    { OUTFILE,  OPT_D, "decpoint", 0 },
    { OUTFILE,  OPT_B, "buffer", 1 }, /* note: 1 is just for backward compat */
    { OUTFILE,  OPT_T, "tempfile", 2 },
    { PANEL,    OPT_B, "between", 0 },
    { PANEL,    OPT_C, "cluster", 2 },
    { PANEL,    OPT_D, "time-dummies", 1 },
    { PANEL,    OPT_E, "nerlove", 0 },
    { PANEL,    OPT_F, "fixed-effects", 0 },
    { PANEL,    OPT_I, "iterate", 0 },
    { PANEL,    OPT_M, "matrix-diff", 0 },
    { PANEL,    OPT_N, "no-df-corr", 0 },
    { PANEL,    OPT_P, "pooled", 0 },
    { PANEL,    OPT_R, "robust", 0 },
    { PANEL,    OPT_U, "random-effects", 0 },
    { PANEL,    OPT_V, "verbose", 0 },
    { PANEL,    OPT_H, "unit-weights", 0 },
    { PANEL,    OPT_X, "unbalanced", 1 },
    { PANPLOT,  OPT_M, "means", 0 },
    { PANPLOT,  OPT_V, "overlay", 0 },
    { PANPLOT,  OPT_S, "sequence", 0 },
    { PANPLOT,  OPT_D, "grid", 0 },
    { PANPLOT,  OPT_A, "stack", 0 },
    { PANPLOT,  OPT_B, "boxplots", 0 },
    { PANPLOT,  OPT_C, "boxplot", 0 },
    { PANPLOT,  OPT_Y, "single-yaxis", 0 },
    { PANSPEC,  OPT_M, "matrix-diff", 0 },
    { PANSPEC,  OPT_N, "nerlove", 0 },
    { POISSON,  OPT_R, "robust", 0 },
    { POISSON,  OPT_C, "cluster", 2 },
    { POISSON,  OPT_V, "verbose", 0 },
    { PCA,      OPT_C, "covariance", 0 },
    { PCA,      OPT_A, "save-all", 0 },
    { PCA,      OPT_O, "save", 1 },
    { PCA,      OPT_Q, "quiet", 0 },
    { PERGM,    OPT_O, "bartlett", 0 },
    { PERGM,    OPT_L, "log", 0 },
    { PERGM,    OPT_R, "radians", 0 },
    { PERGM,    OPT_D, "degrees", 0 },
    { PERGM,    OPT_S, "silent", 0 },
    { PKG,      OPT_L, "local", 0 },
    { PKG,      OPT_V, "verbose", 0 },
    { PLOT,     OPT_C, "control", 0 },
    { PLOT,     OPT_O, "with-lines", 1 },
    { PLOT,     OPT_F, "fit", 2 },
    { PLOT,     OPT_B, "with-boxes", 1 },
    { PLOT,     OPT_Q, "with-steps", 1 },
    { PLOT,     OPT_M, "with-impulses", 1 },
    { PLOT,     OPT_P, "with-lp", 1 },
    { PLOT,     OPT_T, "time-series", 0 },
    { PLOT,     OPT_Y, "single-yaxis", 0 },
    { PLOT,     OPT_Z, "dummy", 0 },
    { PLOT,     OPT_N, "band", 2 },
    { PLOT,     OPT_a, "bands", 2 },
    { PLOT,     OPT_J, "band-style", 2 },
    { PLOT,     OPT_W, "font", 2 },
    { PLOT,     OPT_L, "ylogscale", 1 },
    { PRINT,    OPT_O, "byobs", 0 },
    { PRINT,    OPT_L, "list", 0 },
    { PRINT,    OPT_D, "no-dates", 0 },
    { PRINT,    OPT_U, "numeric", 0 },
    { PRINT,    OPT_M, "midas", 0 },
    { PRINT,    OPT_C, "complex", 0 },
    { PRINT,    OPT_T, "tree", 0 },
    { PRINT,    OPT_R, "range", 2 },
    { PRINT,    OPT_X, "data-only", 0 },
    { PROBIT,   OPT_P, "p-values", 0 },
    { PROBIT,   OPT_R, "robust", 0 },
    { PROBIT,   OPT_C, "cluster", 2 },
    { PROBIT,   OPT_V, "verbose", 0 },
    { PROBIT,   OPT_E, "random-effects", 0 },
    { PROBIT,   OPT_G, "quadpoints", 2 },
    { PROBIT,   OPT_B, "bootstrap", 1 },
    { PROBIT,   OPT_S, "estrella", 0 },
    { QLRTEST,  OPT_L, "limit-to", 2 },
    { QQPLOT,   OPT_R, "raw", 0 },
    { QQPLOT,   OPT_Z, "z-scores", 0 },
    { QUANTREG, OPT_I, "intervals", 1 },
    { QUANTREG, OPT_N, "no-df-corr", 0 },
    { QUANTREG, OPT_R, "robust", 0 },
    { QUIT,     OPT_X, "exit", 0 },
    { RESET,    OPT_C, "cubes-only", 0 },
    { RESET,    OPT_R, "squares-only", 0 },
    { RESET,    OPT_I, "silent", 0 },
    { RESTRICT, OPT_B, "bootstrap", 1 },
    { RESTRICT, OPT_F, "full", 0 },
    { RESTRICT, OPT_J, "jitter", 0 },
    { RESTRICT, OPT_V, "verbose", 0 },
    { RESTRICT, OPT_L, "lbfgs", 0 },
    { RESTRICT, OPT_N, "no-scaling", 0 },
    { RESTRICT, OPT_S, "silent", 0 },
    { RESTRICT, OPT_W, "wald", 0 },
    { RMPLOT,   OPT_T, "trim", 0 },
    { RUNS,     OPT_D, "difference", 0 },
    { RUNS,     OPT_E, "equal", 0 },
    { SCATTERS, OPT_O, "with-lines", 0 },
    { SCATTERS, OPT_T, "time-series", 0 },
    { SCATTERS, OPT_U, "output", 2 },
    { SCATTERS, OPT_X, "matrix", 2 },
    { SCATTERS, OPT_K, "tweaks", 2 },
    { SET,      OPT_F, "from-file", 2 },
    { SET,      OPT_T, "to-file", 2 },
    { SETINFO,  OPT_C, "continuous", 0 },
    { SETINFO,  OPT_D, "discrete", 0 },
    { SETINFO,  OPT_I, "description", 2 },
    { SETINFO,  OPT_G, "graph-name", 2 },
    { SETINFO,  OPT_M, "midas", 0 },
    { SETINFO,  OPT_F, "coded", 0 },
    { SETINFO,  OPT_N, "numeric", 0 },
    { SETOBS,   OPT_C, "stacked-cross-section", 0 },
    { SETOBS,   OPT_P, "panel-vars", 0 },
    { SETOBS,   OPT_R, "restructure", 0 },
    { SETOBS,   OPT_S, "stacked-time-series", 0 },
    { SETOBS,   OPT_T, "time-series", 0 },
    { SETOBS,   OPT_X, "cross-section", 0 },
    { SETOBS,   OPT_N, "special-time-series", 0 },
    { SETOBS,   OPT_G, "panel-groups", 0 },
    { SETOBS,   OPT_I, "panel-time", 0 },
    { SMPL,     OPT_A, "no-all-missing", 0 },
    { SMPL,     OPT_B, "preserve-panel", 0 },
    { SMPL,     OPT_B, "balanced", 0 }, /* alias */
    { SMPL,     OPT_C, "contiguous", 0 },
    { SMPL,     OPT_D, "dates", 0 },
    { SMPL,     OPT_F, "full", 0 },
    { SMPL,     OPT_O, "dummy", 0 },
    { SMPL,     OPT_M, "no-missing", 0 },
    { SMPL,     OPT_N, "random", 0 },
    { SMPL,     OPT_P, "replace", 0 },
    { SMPL,     OPT_R, "restrict", 0 },
    { SMPL,     OPT_T, "permanent", 0 },
    { SMPL,     OPT_U, "unit", 0 },
    { SMPL,     OPT_X, "time", 0 },
    { SPEARMAN, OPT_V, "verbose", 0 },
    { SQUARE,   OPT_O, "cross", 0 },
    { STDIZE,   OPT_C, "center-only", 0 },
    { STDIZE,   OPT_N, "no-df-corr", 0 },
    { STORE,    OPT_A, "matrix", 2 },
    { STORE,    OPT_D, "database", 0 },
    { STORE,    OPT_E, "comment", 2 },
    { STORE,    OPT_F, "overwrite", 0 },
    { STORE,    OPT_G, "dat", 0 },
    { STORE,    OPT_I, "decimal-comma", 0 },
    { STORE,    OPT_J, "jmulti", 0 },
    { STORE,    OPT_L, "lcnames", 0 },
    { STORE,    OPT_M, "gnu-octave", 0 },
    { STORE,    OPT_N, "no-header", 0 },
    { STORE,    OPT_P, "preserve-strvals", 0 },
    { STORE,    OPT_R, "gnu-R", 0 },
    { STORE,    OPT_X, "omit-obs", 0 },
    { STORE,    OPT_Z, "gzipped", 1 },
    { SUMMARY,  OPT_B, "by", 2 },
    { SUMMARY,  OPT_S, "simple", 0 },
    { SUMMARY,  OPT_W, "weights", 2 },
    { SUMMARY,  OPT_X, "matrix", 2 },
    { SYSTEM,   OPT_I, "iterate", 0 },
    { SYSTEM,   OPT_V, "verbose", 0 },
    { SYSTEM,   OPT_R, "robust", 0 },
    { SYSTEM,   OPT_N, "no-df-corr", 0 },
    { TEXTPLOT, OPT_O, "one-scale", 0 },
    { TEXTPLOT, OPT_S, "time-series", 0 },
    { TEXTPLOT, OPT_T, "tall", 0 },
    { TOBIT,    OPT_L, "llimit", 2 },
    { TOBIT,    OPT_M, "rlimit", 2 },
    { TOBIT,    OPT_G, "opg", 0 },
    { TOBIT,    OPT_R, "robust", 0 },
    { TOBIT,    OPT_C, "cluster", 2 },
    { TOBIT,    OPT_V, "verbose", 0 },
    { VAR,      OPT_D, "seasonals", 0 },
    { VAR,      OPT_F, "variance-decomp", 0 },
    { VAR,      OPT_H, "robust-hac", 0 },
    { VAR,      OPT_I, "impulse-responses", 0 },
    { VAR,      OPT_L, "lagselect", 0 },
    { VAR,      OPT_N, "nc", 0 },
    { VAR,      OPT_R, "robust", 0 },
    { VAR,      OPT_T, "trend", 0 },
    { VAR,      OPT_S, "silent", 0 },
    { VAR,      OPT_M, "minlag", 2 },
    { VARLIST,  OPT_A, "accessors", 0 },
    { VARLIST,  OPT_S, "scalars", 0 },
    { VARLIST,  OPT_T, "type", 2 },
    { VARLIST,  OPT_D, "debug" },
    { VECM,     OPT_A, "crt", 0 },
    { VECM,     OPT_D, "seasonals", 0 },
    { VECM,     OPT_F, "variance-decomp", 0 },
    { VECM,     OPT_I, "impulse-responses", 1 },
    { VECM,     OPT_N, "nc", 0 },
    { VECM,     OPT_R, "rc", 0 },
    { VECM,     OPT_C, "uc", 0 },
    { VECM,     OPT_T, "ct", 0 },
    { VECM,     OPT_V, "verbose", 0 },
    { VECM,     OPT_S, "silent", 0 },
    { WLS,      OPT_R, "robust", 0 },
    { WLS,      OPT_C, "cluster", 2 },
    { WLS,      OPT_Z, "allow-zeros", 0 },
    { XTAB,     OPT_C, "column", 0 },
    { XTAB,     OPT_X, "matrix", 2 },
    { XTAB,     OPT_R, "row", 0 },
    { XTAB,     OPT_Z, "zeros", 0 },
    { XTAB,     OPT_T, "tex", 1 },
    { XTAB,     OPT_N, "no-totals", 0 },
    { XTAB,     OPT_E, "equal", 0 },
    { XTAB,     OPT_F, "no-fisher", 0 },
    { 0,        0L,    NULL, 0 }
};

static const char *get_longopt (int ci, gretlopt opt)
{
    int i, got_ci = 0;

    for (i=0; gretl_opts[i].ci; i++) {
        if (gretl_opts[i].ci == ci) {
            if (gretl_opts[i].o == opt) {
                return gretl_opts[i].longopt;
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    return "??";
}

int cluster_option_ok (int ci)
{
    int i, got_ci = 0;

    for (i=0; gretl_opts[i].ci; i++) {
        if (gretl_opts[i].ci == ci) {
            if (gretl_opts[i].o == OPT_C &&
                !strcmp(gretl_opts[i].longopt, "cluster")) {
                return 1;
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    return 0;
}

/* used in tokenize.c: detects the case where a command
   has an option of the form --matrix=foo, indicating that
   the columns of a named matrix replace the dataset
   series that the comand would otherwise require as
   arguments
*/

int matrix_data_option (int ci, gretlopt opt)
{
    if (opt & OPT_X) {
        int i, got_ci = 0;

        for (i=0; gretl_opts[i].ci; i++) {
            if (gretl_opts[i].ci == ci) {
                if (gretl_opts[i].o == OPT_X &&
                    !strcmp(gretl_opts[i].longopt, "matrix")) {
                    return 1;
                }
                got_ci = 1;
            } else if (got_ci) {
                break;
            }
        }
    }

    return 0;
}

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
                strings_array_free(optstrs, n);
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

/* used in checking command documentation */

const char **get_opts_for_command (int ci, int *nopt)
{
    int i, j, got_ci, n = 0;
    const char **ret = NULL;

    if (ci != OLS) {
        /* widely applicable options attached to OLS */
        n += vcv_opt_ok(ci);
        n += quiet_opt_ok(ci);
        n += window_opt_ok(ci);
    }

    if (ci != GNUPLOT) {
	/* common plotting options attached to GNUPLOT */
	n += plot_output_opt_ok(ci);
	n += plot_outbuf_opt_ok(ci);
    }

    if (ci != CORR) {
	/* auxiliary plotting options attached to CORR */
	n += cmd_plot_opt_ok(ci);
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

    j = got_ci = 0;
    for (i=0; gretl_opts[i].ci != 0; i++) {
        if (gretl_opts[i].ci == ci) {
            ret[j++] = gretl_opts[i].longopt;
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    if (ci != OLS) {
        if (vcv_opt_ok(ci)) {
            ret[j++] = "vcv";
        }
        if (quiet_opt_ok(ci)) {
            ret[j++] = "quiet";
        }
        if (window_opt_ok(ci)) {
            ret[j++] = "window";
        }
    } else if (ci != GNUPLOT) {
	if (plot_output_opt_ok(ci)) {
	    ret[j++] = "output";
	}
	if (plot_outbuf_opt_ok(ci)) {
	    ret[j++] = "outbuf";
	}
    } else if (ci != CORR) {
	if (cmd_plot_opt_ok(ci)) {
	    ret[j++] = "plot";
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

    return OPT_NONE;
}

static int opt_is_valid (gretlopt opt, int ci, char c)
{
    int i, got_ci = 0;

    if (opt == OPT_O && vcv_opt_ok(ci)) {
        return 1;
    } else if (opt == OPT_Q && quiet_opt_ok(ci)) {
        return 1;
    } else if (opt == OPT_W && window_opt_ok(ci)) {
        return 1;
    } else if (opt == OPT_U && plot_output_opt_ok(ci)) {
	return 1;
    } else if (opt == OPT_b && plot_outbuf_opt_ok(ci)) {
	return 1;
    } else if (opt == OPT_U && cmd_plot_opt_ok(ci)) {
	return 1;
    }

    for (i=0; gretl_opts[i].ci != 0; i++) {
        if (ci == gretl_opts[i].ci) {
            if (opt == gretl_opts[i].o) {
                return 1;
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    if (c != 0) {
        gretl_errmsg_sprintf("Invalid option '-%c'", c);
    }

    return 0;
}

enum {
    OPT_SETOPT  = 1 << 0,
    OPT_PERSIST = 1 << 1
};

/* The following apparatus is used for

   (a) setting and retrieving parameters associated with
       command options, as in --opt=val, and

   (b) storing options for a specified command via the
       "setopt" command (with or without parameters).
*/

typedef struct stored_opt_ stored_opt;

struct stored_opt_ {
    int ci;       /* index of the associated command */
    gretlopt opt; /* option flag */
    char *val;    /* option parameter value, or NULL */
    int flags;    /* may include OPT_SETOPT, OPT_PERSIST */
    int fd;       /* "function depth" at which registered */
};

static stored_opt *optinfo;
static int n_stored_opts;

static void clear_one_option (stored_opt *so)
{
    free(so->val);
    so->val = NULL;
    so->ci = 0;
    so->opt = 0;
    so->flags = 0;
}

/**
 * clear_stored_options_for_command:
 * @ci: target command index.
 *
 * Clears any (non-persistent) option information currently
 * associated with the given command, identified by its
 * index.
 */

void clear_stored_options_for_command (int ci)
{
    int i, fd = gretl_function_depth();

#if OPTDEBUG > 1
    fprintf(stderr, "clearing stored options for %s\n",
            gretl_command_word(ci));
#endif

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].fd == fd && optinfo[i].ci == ci &&
            !(optinfo[i].flags & OPT_PERSIST)) {
            clear_one_option(&optinfo[i]);
        }
    }
}

/* called on exiting a user-defined function */

void destroy_option_params_at_level (int level)
{
    int i, n = n_stored_opts;

#if OPTDEBUG
    fprintf(stderr, "destroy_option_params_at_level: %d\n", level);
#endif

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].fd == level) {
            clear_one_option(&optinfo[i]);
            n--;
        }
    }

    if (n == 0) {
        free(optinfo);
        optinfo = NULL;
    }

    n_stored_opts = n;
}

/* unconditionally clean up the entire optinfo stack */

void stored_options_cleanup (void)
{
    int i;

#if OPTDEBUG
    fprintf(stderr, "stored_options_cleanup\n");
#endif

    for (i=0; i<n_stored_opts; i++) {
        free(optinfo[i].val);
    }

    free(optinfo);
    optinfo = NULL;
    n_stored_opts = 0;
}

/* scrub just those options set via "setopt" */

void setopt_cleanup (void)
{
    int i, n = n_stored_opts;

#if OPTDEBUG
    fprintf(stderr, "setopt_cleanup\n");
#endif

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].flags & OPT_SETOPT) {
            clear_one_option(&optinfo[i]);
            n--;
        }
    }

    if (n == 0) {
        free(optinfo);
        optinfo = NULL;
    }

    n_stored_opts = n;
}

static stored_opt *matching_stored_opt (int ci, gretlopt opt)
{
    int i, fd = gretl_function_depth();

#if OPTDEBUG
    fprintf(stderr, "matching_stored_opt? ci=%d, fd=%d, opt=%d, n_stored=%d\n",
            ci, fd, opt, n_stored_opts);
#endif

    for (i=0; i<n_stored_opts; i++) {
	stored_opt *so = &optinfo[i];

        if (so->ci == ci && so->opt == opt && so->fd == fd) {
            return so;
        }
    }

    return NULL;
}

static stored_opt *empty_stored_opt_slot (void)
{
    int i;

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].ci == 0) {
            return &optinfo[i];
        }
    }

    return NULL;
}

/* for a given (@ci, @opt) pair, determine its status with regard
   to a parameter value: 0 = not allowed, 1 = optional, 2 = required
*/

static int option_parm_status (int ci, gretlopt opt)
{
    int i, got_ci = 0;

    if (opt == OPT_U) {
	if (plot_output_opt_ok(ci)) {
	    ci = GNUPLOT;
	} else if (cmd_plot_opt_ok(ci)) {
	    ci = CORR;
	}
    } else if (opt == OPT_b && plot_outbuf_opt_ok(ci)) {
        ci = GNUPLOT;
    }

    for (i=0; gretl_opts[i].ci != 0; i++) {
        if (gretl_opts[i].ci == ci) {
            if (gretl_opts[i].o == opt) {
                return gretl_opts[i].parminfo;
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    return 0;
}

/* Handle the case where we got input on the pattern

   --option-with-param=eval(foo)

   where "foo" is the name of a string variable and the
   option parameter should be set to the value of that
   variable.
*/

static int maybe_evaluate_optval (stored_opt *so)
{
    char *s = so->val;
    int n = strlen(s);
    int err = 0;

    if (!strncmp(s, "eval(", 5) && s[n-1] == ')') {
        char *sname = gretl_strndup(s + 5, n - 6);

        if (sname != NULL) {
            char *tmp = generate_string(sname, NULL, &err);

            if (tmp != NULL) {
                free(so->val);
                so->val = tmp;
            }
            free(sname);
        }
    }

    return err;
}

static int real_push_option (int ci, gretlopt opt, char *val,
                             int checked, int flags)
{
    int fd = gretl_function_depth();
    int n, err = 0;
    int val_set = 0;
    stored_opt *so;

    if (!checked && ci > 0 && option_parm_status(ci, opt) == OPT_NO_PARM)  {
        return E_DATA;
    }

#if OPTDEBUG
    fprintf(stderr, "push_option_param: ci=%d (%s), fd=%d, opt=%d,"
            " val='%s', SETOPT %d\n", ci, gretl_command_word(ci),
	    fd, opt, val, flags & OPT_SETOPT ? 1 : 0);
#endif

    so = matching_stored_opt(ci, opt);

    if (so != NULL) {
        /* got a match for the (ci, opt) pair already */
#if OPTDEBUG
        fprintf(stderr, " push_option_param: replacing\n");
#endif
        if (!(flags & OPT_SETOPT)) {
            free(so->val);
            so->val = val;
            val_set = 1;
        }
        so->flags = flags;
        goto finish;
    }

    so = empty_stored_opt_slot();

    if (so != NULL) {
        /* re-use a vacant slot */
#if OPTDEBUG
        fprintf(stderr, " push_option_param: reusing empty\n");
#endif
        so->ci = ci;
        so->opt = opt;
        so->val = val;
        val_set = 1;
        so->flags = flags;
        so->fd = fd;
        goto finish;
    }

    /* so we have to extend the array */
    n = n_stored_opts + 1;
    so = realloc(optinfo, n * sizeof *so);

    if (so == NULL) {
        err = E_ALLOC;
    } else {
#if OPTDEBUG
        fprintf(stderr, " push_option_param: appending\n");
#endif
        optinfo = so;
        so = &optinfo[n-1];
        so->ci = ci;
        so->opt = opt;
        so->val = val;
        val_set = 1;
        so->flags = flags;
        so->fd = fd;
        n_stored_opts = n;
    }

 finish:

    if (val_set && so->val != NULL) {
        maybe_evaluate_optval(so);
    }

#if OPTDEBUG
    fprintf(stderr, " push_option_param: returning err = %d\n", err);
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
 * to be associated with @opt for the current @ci. Note
 * that the command option apparatus takes ownership of
 * @val, so the value passed in should be copied if need
 * be.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int push_option_param (int ci, gretlopt opt, char *val)
{
    return real_push_option(ci, opt, val, 0, 0);
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
    stored_opt *so = matching_stored_opt(ci, opt);

    return (so != NULL)? so->val : NULL;
}

/**
 * get_optval_double:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @err: location to receive error code.
 *
 * Returns: the double-precision ancillary value currently
 * associated with option @opt for command @ci, if any,
 * otherwise #NADBL. If @opt is an active option for
 * @ci but the parameter for this option cannot be
 * interpreted as a numerical value, E_INVARG is written
 * into @err.
 */

double get_optval_double (int ci, gretlopt opt, int *err)
{
    stored_opt *so = matching_stored_opt(ci, opt);
    double ret = NADBL;

    if (so != NULL && so->val != NULL) {
        ret = gretl_double_from_string(so->val, err);
        if (err) {
            ret = generate_scalar(so->val, NULL, err);
        }
        if (*err) {
            gretl_errmsg_sprintf(_("%s: invalid option argument"), so->val);
            *err = E_INVARG;
        }
    }

    return ret;
}

/**
 * get_optval_int:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @err: location to receive error code.
 *
 * Returns: the integer ancillary value currently
 * associated with option @opt for command @ci, if any,
 * otherwise 0. A non-zero value is written to @err if
 * such a value is required for the option in question
 * but is not present.
 */

int get_optval_int (int ci, gretlopt opt, int *err)
{
    stored_opt *so = matching_stored_opt(ci, opt);
    int status = option_parm_status(ci, opt);
    int ret = 0;

    if (so != NULL && so->val != NULL) {
        ret = gretl_int_from_string(so->val, err);
        if (*err) {
            ret = generate_int(so->val, NULL, err);
        }
        if (*err) {
            gretl_errmsg_sprintf(_("%s: invalid option argument"), so->val);
            *err = E_INVARG;
        }
    } else if (status == 2 && err != NULL) {
        const char *longopt = get_longopt(ci, opt);

        gretl_errmsg_sprintf(_("The option '--%s' requires a parameter"),
                             longopt);
        *err = E_DATA;
    }

    return ret;
}

static void clear_options_for_command (int ci)
{
    int i, fd = gretl_function_depth();

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].fd == fd && optinfo[i].ci == ci) {
            clear_one_option(&optinfo[i]);
        }
    }
}

static void set_stored_options (int ci, gretlopt opt, int flags)
{
    int i, got_ci = 0;

#if OPTDEBUG
    fprintf(stderr, "setting stored options for %s (%d): ",
            gretl_command_word(ci), ci);
    debug_print_option_flags(NULL, opt);
#endif

    if (vcv_opt_ok(ci) && (opt & OPT_O)) {
        real_push_option(ci, OPT_O, NULL, 1, flags);
        opt &= ~OPT_O; /* handled */
    }

    if (quiet_opt_ok(ci) && (opt & OPT_Q)) {
        real_push_option(ci, OPT_Q, NULL, 1, flags);
        opt &= ~OPT_Q; /* handled */
    }

    if (window_opt_ok(ci) && (opt & OPT_W)) {
        real_push_option(ci, OPT_W, NULL, 1, flags);
        opt &= ~OPT_W; /* handled */
    }

    if (opt == 0) {
        return;
    }

    for (i=0; gretl_opts[i].o != 0; i++) {
	struct gretl_option *gopt = &gretl_opts[i];

        if (ci == gopt->ci) {
            if (opt & gopt->o) {
		real_push_option(ci, gopt->o, NULL, 1, flags);
            }
            got_ci = ci;
        } else if (got_ci > 0 && ci != got_ci) {
            break;
        }
    }
}

/* Apparatus for pre-selecting options for a specified command,
   using the "setopt" command.
*/

int set_options_for_command (const char *cmdword,
                             const char *param,
                             gretlopt opt)
{
    int target_ci = gretl_command_number(cmdword);
    int flags = OPT_SETOPT;
    int clear = 0;
    int err = 0;

    if (target_ci == 0 || target_ci == SETOPT) {
        gretl_errmsg_sprintf(_("field '%s' in command is invalid"),
                             cmdword);
        return E_DATA;
    }

    if (param != NULL && *param != '\0') {
        if (!strcmp(param, "persist")) {
            flags |= OPT_PERSIST;
        } else if (!strcmp(param, "clear")) {
            clear = 1;
        } else {
            gretl_errmsg_sprintf(_("field '%s' in command is invalid"),
                                 param);
            return E_DATA;
        }
    }

    if (clear) {
        clear_options_for_command(target_ci);
    } else if (opt == 0) {
        err = E_ARGS;
    } else {
        set_stored_options(target_ci, opt, flags);
    }

    return err;
}

#if OPTDEBUG

static char option_flag_char (gretlopt opt)
{
    char c = 'A';
    int i;

    for (i=OPT_A; i<=OPT_Y; i*=2) {
        if (opt == i) {
            return c;
        }
        c++;
    }

    return '?';
}

#endif

void maybe_get_stored_options (int ci, gretlopt *popt)
{
    int i, fd = gretl_function_depth();

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].fd == fd && optinfo[i].ci == ci &&
            (optinfo[i].flags & OPT_SETOPT)) {
#if OPTDEBUG
            fprintf(stderr, "ci %d: got stored OPT_%c\n",
                    ci, option_flag_char(optinfo[i].opt));
#endif
            *popt |= optinfo[i].opt;
        }
    }
}

/* When writing data as gdt this is called conditionally on OPT_Z
   being given, but when writing gdtb it is called unconditionally
   (since the format of the latter is a zipfile in all cases).
   The effects are then:

   gdt: no compression applied by default; if --gzipped given
   with no parameter, zlib level 1 used; if --gzipped given
   with parameter, param value respected within range 0 to 9.

   gdtb: zlib level 1 used by default, levels 0 to 9 can be
   specified via --gzipped + param.
*/

int get_compression_option (int ci)
{
    stored_opt *so = matching_stored_opt(ci, OPT_Z);
    int level = 0;

    if (so == NULL || so->val == NULL) {
        /* use zlib level 1 by default */
        level = 1;
    } else {
        level = atoi(so->val);
        if (level < 0) {
            level = 0;
        } else if (level > 9) {
            level = 9;
        }
    }

    return level;
}

/* below: called via GUI */

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

    return real_push_option(ci, opt, gretl_strdup(s), 1, 0);
}

/**
 * set_optval_int:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @k: value to set.
 *
 * Sets a integer ancillary value to be associated
 * with option @opt for command @ci.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int set_optval_int (int ci, gretlopt opt, int k)
{
    char *s = gretl_strdup_printf("%d", k);

    return real_push_option(ci, opt, s, 1, 0);
}

/**
 * set_optval_string:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @s: value to set.
 *
 * Sets a ancillary string value to be associated
 * with option @opt for command @ci.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int set_optval_string (int ci, gretlopt opt, const char *s)
{
    if (s == NULL) {
	gretl_errmsg_set("set_optval_string: no string was supplied");
	return E_DATA;
    } else {
	char *scpy = gretl_strdup(s);

	if (scpy == NULL) {
	    return E_ALLOC;
	} else {
	    return real_push_option(ci, opt, scpy, 1, 0);
	}
    }
}

/* valid_long_opt: this is (semi-) public because we have need of
   it in the command tokenizer
*/

gretlopt valid_long_opt (int ci, const char *s, OptStatus *status)
{
    gretlopt opt = OPT_NONE;
    int i, got_ci = 0;

    *status = 0;

    if (*s == '\0') {
        return 0;
    }

#if OPTDEBUG
    fprintf(stderr, "valid_long_opt, s = '%s'\n", s);
#endif

    /* common options without parameter */
    if (vcv_opt_ok(ci) && !strcmp(s, "vcv")) {
        return OPT_O;
    }
    if (quiet_opt_ok(ci) && !strcmp(s, "quiet")) {
        return OPT_Q;
    }
    if (window_opt_ok(ci) && !strcmp(s, "window")) {
        return OPT_W;
    }

    /* common options with parameter: switch @ci to the
       command that "owns" the option in question
    */
    if (plot_output_opt_ok(ci) && !strcmp(s, "output")) {
        ci = GNUPLOT;
    }
    if (plot_outbuf_opt_ok(ci) && (!strcmp(s, "outbuf") ||
				   !strcmp(s, "buffer"))) {
        ci = GNUPLOT;
    }
    if (cmd_plot_opt_ok(ci) && !strcmp(s, "plot")) {
        ci = CORR;
    }

    /* start by looking for an exact match */
    for (i=0; gretl_opts[i].o != 0; i++) {
        if (ci == gretl_opts[i].ci) {
            if (!strcmp(s, gretl_opts[i].longopt)) {
                opt = gretl_opts[i].o;
                *status = gretl_opts[i].parminfo;
                break;
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    /* if this failed, try for a unique abbreviation */
    if (opt == OPT_NONE) {
        int optlen, slen = strlen(s);
        int nmatch = 0;

        got_ci = 0;
        for (i=0; gretl_opts[i].o != 0; i++) {
            if (ci == gretl_opts[i].ci) {
                optlen = strlen(gretl_opts[i].longopt);
                if (optlen > slen && !strncmp(s, gretl_opts[i].longopt, slen)) {
                    opt = gretl_opts[i].o;
                    *status = gretl_opts[i].parminfo;
                    nmatch++;
                }
                got_ci = 1;
            } else if (got_ci) {
                break;
            }
        }
        if (nmatch > 1) {
            if (ci == SETOBS && !strcmp(s, "panel")) {
                /* backward compatibility: --panel-vars was there first */
                return OPT_P;
            }
            *status = OPT_AMBIGUOUS;
            return OPT_NONE;
        }
    }

    /* backward compatibility */
    if (opt == OPT_NONE && !strcmp(s, "wald")) {
        opt = OPT_W;
        *status = 0;
    }

#if OPTDEBUG
    fprintf(stderr, "valid_long_opt, returning %d (status %d)\n",
	    opt, *status);
#endif

    return opt;
}

/* see valid_long_opt() above */

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

static void print_option_param (const char *s, PRN *prn)
{
    const char *qchars = "=%, ";
    const char *p = s;
    int wrap = 0;
    int escape = 0;

    while (*p) {
        if (strspn(p, qchars)) {
            wrap = 1;
        } else if (*p == '"') {
            escape = 1;
        }
        p++;
    }

    if (wrap) {
        if (escape) {
            pputs(prn, "=\"");
            while (*s) {
                if (*s == '"') {
                    pputs(prn, "\\\"");
                } else {
                    pputc(prn, *s);
                }
                s++;
            }
            pputs(prn, "\"\n");
        } else {
            pprintf(prn, "=\"%s\"", s);
        }
    } else {
        pprintf(prn, "=%s", s);
    }
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
    int i, got_ci;

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

    if ((oflags & OPT_O) && vcv_opt_ok(ci)) {
        pputs(flagprn, " --vcv");
        oflags &= ~OPT_O; /* handled */
    }
    if ((oflags & OPT_Q) && quiet_opt_ok(ci)) {
        pputs(flagprn, " --quiet");
        oflags &= ~OPT_Q; /* handled */
    }
    if ((oflags & OPT_W) && window_opt_ok(ci)) {
        pputs(flagprn, " --window");
        oflags &= ~OPT_W; /* handled */
    }

    if (plot_outbuf_opt_ok(ci)) {
	if (oflags & OPT_b) {
	    parm = get_optval_string(ci, OPT_b);
	    pprintf(flagprn, " --outbuf=%s\n", parm);
	    oflags &= ~OPT_b; /* handled */
	} else if (oflags & OPT_U) {
	    parm = get_optval_string(ci, OPT_U);
	    if (plot_output_opt_ok(ci)) {
		pprintf(flagprn, " --output=%s\n", parm);
	    } else {
		pprintf(flagprn, " --plot=%s\n", parm);
	    }
	    oflags &= ~OPT_U; /* handled */
	}
    }

    got_ci = 0;
    for (i=0; gretl_opts[i].ci != 0; i++) {
        if (ci == gretl_opts[i].ci) {
            opt = gretl_opts[i].o;
            if (oflags & opt) {
                pprintf(flagprn, " --%s", gretl_opts[i].longopt);
                if (gretl_opts[i].parminfo) {
                    parm = get_optval_string(ci, opt);
                    if (parm != NULL && *parm != '\0') {
                        print_option_param(parm, flagprn);
                    }
                }
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    return gretl_print_get_buffer(flagprn);
}

void option_printing_cleanup (void)
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

    for (o=OPT_A; o<=OPT_i; o=o<<1) {
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
 * options_incompatible_with:
 * @opt: option flags to be tested.
 * @base: "base" option.
 * @test: bitwise OR of flags that are incompatible with @base.
 *
 * Returns: %E_BADOPT if @opt contains both @base and one or more of
 * the flags in @test, otherwise 0.
 */

int options_incompatible_with (gretlopt opt, gretlopt base,
                               gretlopt test)
{
    if (opt & base) {
        if (opt & test) {
            return E_BADOPT;
        }
    }

    return 0;
}

/**
 * option_prereq_missing:
 * @opt: option flags to be tested.
 * @test: bitwise OR of flags that have a definite prequisite.
 * @prereq: bitwise OR of prerequisite flags.
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
 * Flags an error: to be used when @opt is not applicable
 * for command @ci, in context.
 *
 * Returns: %E_BADOPT.
 */

int inapplicable_option_error (int ci, gretlopt opt)
{
    const char *s = print_flags(opt, ci);

    gretl_errmsg_sprintf(_("%s: inapplicable option"), s);
    return E_BADOPT;
}

void debug_print_option_flags (const char *msg, gretlopt opt)
{
    if (msg != NULL && *msg != '\0') {
        fprintf(stderr, "%s: ", msg);
    }

    if (opt == 0) {
        fprintf(stderr, "opt=0\n");
    } else {
        char c = 'A';
        int i, started = 0;

        fprintf(stderr, "opt=%d (", opt);

        for (i=OPT_A; i<=OPT_Y; i*=2) {
            if (opt & i) {
                if (started) {
                    fputc('|', stderr);
                }
                fprintf(stderr, "OPT_%c", c);
                started = 1;
            }
            c++;
        }

        fputs(")\n", stderr);
    }
}
