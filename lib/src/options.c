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
#include <errno.h>

/* model commands plus ADD and OMIT */
#define vcv_opt_ok(c) (c == ADD || \
                       c == AR || \
		       c == AR1 || \
                       c == ARBOND || \
                       c == ARMA || \
                       c == HECKIT || \
                       c == GARCH || \
                       c == GMM || \
                       c == HCCM || \
                       c == HSK || \
                       c == LAD || \
                       c == LOGISTIC || \
                       c == LOGIT || \
                       c == MPOLS || \
                       c == OLS || \
                       c == OMIT || \
                       c == MLE || \
                       c == NLS || \
                       c == PANEL || \
                       c == POISSON || \
                       c == PROBIT || \
		       c == QUANTREG || \
                       c == TOBIT || \
                       c == TSLS || \
                       c == WLS)

struct gretl_option {
    int ci;              /* command index (gives context) */
    gretlopt o;          /* index of integer type */
    const char *longopt; /* "--"-style string representation of option */
};

struct flag_match {
    gretlopt o;
    unsigned char c;
};

/* Below: This is used as a one-way mapping from the long form
   to the index (e.g. OPT_Q), so a given index can have more than 
   one long-form counterpart, depending on context. 
*/

struct gretl_option gretl_opts[] = {
    { ADD,      OPT_B, "both" },
    { ADD,      OPT_I, "silent" },
    { ADD,      OPT_Q, "quiet" },
    { ADD,      OPT_T, "inst" },
    { ADF,      OPT_N, "nc" }, 
    { ADF,      OPT_C, "c" }, 
    { ADF,      OPT_D, "seasonals" },
    { ADF,      OPT_R, "ctt" },     
    { ADF,      OPT_T, "ct" }, 
    { ADF,      OPT_V, "verbose" },
    { ADF,      OPT_Q, "quiet" },
    { ADF,      OPT_F, "difference" },
    { ADF,      OPT_E, "test-down" },
    { ADF,      OPT_G, "gls" },
    { AR1,      OPT_B, "no-corc"},
    { AR1,      OPT_H, "hilu"},
    { AR1,      OPT_P, "pwe"},
    { AR1,      OPT_Q, "quiet"},
    { APPEND,   OPT_T, "time-series" },
    { APPEND,   OPT_Q, "quiet" },
    { ARBOND,   OPT_A, "asymptotic" },
    { ARBOND,   OPT_D, "time-dummies" },
    { ARBOND,   OPT_H, "orthdev" },
    { ARBOND,   OPT_T, "two-step" },
    { ARMA,     OPT_C, "conditional" },
    { ARMA,     OPT_G, "opg" },
    { ARMA,     OPT_H, "hessian" },
    { ARMA,     OPT_N, "nc" },    
    { ARMA,     OPT_Q, "quiet" },
    { ARMA,     OPT_V, "verbose" },
    { ARMA,     OPT_X, "x-12-arima" },
    { BXPLOT,   OPT_O, "notches" },
    { CHOW,     OPT_Q, "quiet" },
    { COINT,    OPT_N, "nc" },
    { COINT,    OPT_R, "ctt" },     
    { COINT,    OPT_S, "skip-df" },
    { COINT,    OPT_T, "ct" },
    { COINT2,   OPT_A, "crt" },
    { COINT2,   OPT_D, "seasonals" },
    { COINT2,   OPT_N, "nc" },
    { COINT2,   OPT_Q, "quiet" },
    { COINT2,   OPT_R, "rc" },
    { COINT2,   OPT_T, "ct" },
    { COINT2,   OPT_V, "verbose" },
    { CORR,     OPT_K, "kendall" },
    { CORR,     OPT_S, "spearman" },
    { CORR,     OPT_U, "uniform" },
    { CORR,     OPT_V, "verbose" },
    { CORRGM,   OPT_Q, "quiet" },
    { CUSUM,    OPT_Q, "quiet" },
    { CUSUM,    OPT_R, "squares" },
    { DATA,     OPT_O, "odbc" },
    { DATAMOD,  OPT_P, "preserve" },
    { DELEET,   OPT_D, "db" },
    { DIFFTEST, OPT_G, "sign" },
    { DIFFTEST, OPT_R, "rank-sum" },
    { DIFFTEST, OPT_I, "signed-rank" },
    { DIFFTEST, OPT_V, "verbose" },
    { DISCRETE, OPT_R, "reverse" },
    { DUMMIFY,  OPT_F, "drop-first" },
    { DUMMIFY,  OPT_L, "drop-last" },
    { EQNPRINT, OPT_O, "complete" },
    { EQNPRINT, OPT_T, "t-ratios" },
    { TABPRINT, OPT_O, "complete" },
    { TABPRINT, OPT_R, "rtf" },
    { ESTIMATE, OPT_M, "geomean" },
    { ESTIMATE, OPT_N, "no-df-corr" },
    { ESTIMATE, OPT_Q, "quiet" },
    { ESTIMATE, OPT_T, "iterate" },
    { ESTIMATE, OPT_V, "verbose" },
    { FCAST,    OPT_D, "dynamic" },
    { FCAST,    OPT_S, "static" },
    { FCAST,    OPT_Q, "quiet" },
    { FCAST,    OPT_O, "out-of-sample" },
    { FOREIGN,  OPT_D, "send-data" },
    { FOREIGN,  OPT_Q, "quiet" },
    { FREQ,     OPT_O, "gamma" },
    { FREQ,     OPT_Q, "quiet" },
    { FREQ,     OPT_S, "silent" },
    { FREQ,     OPT_Z, "normal" },
    { GARCH,    OPT_A, "arma-init" },
    { GARCH,    OPT_F, "fcp" }, 
    { GARCH,    OPT_N, "nc" }, 
    { GARCH,    OPT_R, "robust" },
    { GARCH,    OPT_V, "verbose" },
    { GMM,      OPT_I, "iterate" },
    { GMM,      OPT_Q, "quiet" },
    { GMM,      OPT_T, "two-step" },
    { GMM,      OPT_V, "verbose" },
    { GNUPLOT,  OPT_O, "with-lines" },
    { GNUPLOT,  OPT_I, "inverse-fit" },
    { GNUPLOT,  OPT_L, "loess-fit" },
    { GNUPLOT,  OPT_Q, "quadratic-fit" },
    { GNUPLOT,  OPT_N, "linear-fit" },
    { GNUPLOT,  OPT_M, "with-impulses" },
    { GNUPLOT,  OPT_S, "suppress-fitted" },
    { GNUPLOT,  OPT_T, "time-series" },
    { GNUPLOT,  OPT_Z, "dummy" },
    { GNUPLOT,  OPT_C, "control" },
    { GRAPH,    OPT_O, "tall" },
    { HECKIT,   OPT_M, "ml" },
    { HECKIT,   OPT_T, "two-step" },
    { HECKIT,   OPT_V, "verbose" },
    { HELP,     OPT_F, "func" },
    { HSK,      OPT_Q, "quiet" },
    { KPSS,     OPT_T, "trend" },
    { KPSS,     OPT_V, "verbose" },
    { KPSS,     OPT_Q, "quiet" },
    { KPSS,     OPT_F, "difference" },
    { LEVERAGE, OPT_S, "save" },
    { LMTEST,   OPT_A, "autocorr" },
    { LMTEST,   OPT_B, "breusch-pagan" },
    { LMTEST,   OPT_H, "arch" },
    { LMTEST,   OPT_L, "logs" },
    { LMTEST,   OPT_S, "squares" }, 
    { LMTEST,   OPT_P, "panel" },
    { LMTEST,   OPT_R, "robust" },
    { LMTEST,   OPT_Q, "quiet" },
    { LMTEST,   OPT_W, "white" },
    { LMTEST,   OPT_X, "white-nocross" },
    { LOGIT,    OPT_P, "p-values" },
    { LOGIT,    OPT_Q, "quiet" },
    { LOGIT,    OPT_R, "robust" },
    { LOGIT,    OPT_V, "verbose" },
    { LOOP,     OPT_P, "progressive" },
    { LOOP,     OPT_Q, "quiet" },
    { LOOP,     OPT_V, "verbose" },
    { MAHAL,    OPT_S, "save" },
    { MAHAL,    OPT_V, "vcv" },
    { MEANTEST, OPT_O, "unequal-vars" },
    { MLE,      OPT_H, "hessian" },
    { MLE,      OPT_N, "numerical" },
    { MLE,      OPT_Q, "quiet" },
    { MLE,      OPT_R, "robust" },
    { MLE,      OPT_V, "verbose" },
    { MODPRINT, OPT_C, "csv" },
    { MODPRINT, OPT_O, "complete" },
    { MODPRINT, OPT_R, "rtf" },
    { MODPRINT, OPT_T, "tex" },
    { MPOLS,    OPT_O, "vcv" },
    { MPOLS,    OPT_Q, "quiet" },
    { MPOLS,    OPT_S, "simple-print" },
    { NLS,      OPT_N, "numerical" },
    { NLS,      OPT_O, "vcv" },
    { NLS,      OPT_Q, "quiet" },
    { NLS,      OPT_R, "robust" },
    { NLS,      OPT_V, "verbose" },
    { NORMTEST, OPT_A, "all" },
    { NORMTEST, OPT_D, "dhansen" },
    { NORMTEST, OPT_W, "swilk" },
    { NORMTEST, OPT_J, "jbera" },
    { NORMTEST, OPT_L, "lillie" },
    { NORMTEST, OPT_Q, "quiet" },
    { NULLDATA, OPT_P, "preserve" },
    { OLS,      OPT_F, "print-final" },
    { OLS,      OPT_N, "no-df-corr" },
    { OLS,      OPT_O, "vcv" }, 
    { OLS,      OPT_R, "robust" },
    { OLS,      OPT_Q, "quiet" },
    { OLS,      OPT_S, "simple-print" },
    { OLS,      OPT_V, "anova" },
    { OMIT,     OPT_A, "auto" },
    { OMIT,     OPT_B, "both" },
    { OMIT,     OPT_I, "silent" },
    { OMIT,     OPT_P, "bootstrap" },
    { OMIT,     OPT_Q, "quiet" },
    { OMIT,     OPT_T, "inst" },
    { OMIT,     OPT_W, "wald" },
    { OPEN,     OPT_C, "coded" },
    { OPEN,     OPT_D, "drop-empty" },
    { OPEN,     OPT_O, "odbc" },
    { OPEN,     OPT_P, "preserve" },
    { OPEN,     OPT_W, "www" },
    { OPEN,     OPT_Q, "quiet" },
    { OUTFILE,  OPT_A, "append" },
    { OUTFILE,  OPT_C, "close" },
    { OUTFILE,  OPT_W, "write" },
    { PANEL,    OPT_B, "between" },
    { PANEL,    OPT_D, "time-dummies" },
    { PANEL,    OPT_F, "fixed-effects" },
    { PANEL,    OPT_H, "hausman-reg" },
    { PANEL,    OPT_O, "vcv" },
    { PANEL,    OPT_P, "pooled" },
    { PANEL,    OPT_Q, "quiet" },
    { PANEL,    OPT_R, "robust" },
    { PANEL,    OPT_S, "silent" },
    { PANEL,    OPT_T, "iterate" },
    { PANEL,    OPT_U, "random-effects" },
    { PANEL,    OPT_V, "verbose" },
    { PANEL,    OPT_W, "unit-weights" },
    { POISSON,  OPT_V, "verbose" },
    { PCA,      OPT_C, "covariance" },
    { PCA,      OPT_A, "save-all" },
    { PCA,      OPT_O, "save" },
    { PERGM,    OPT_O, "bartlett" },
    { PERGM,    OPT_L, "log" },
    { PLOT,     OPT_O, "one-scale" },
    { PRINT,    OPT_O, "byobs" },
    { PRINT,    OPT_L, "long" },
    { PRINT,    OPT_N, "no-dates" },
    { PRINT,    OPT_T, "ten" },
    { PROBIT,   OPT_P, "p-values" },
    { PROBIT,   OPT_Q, "quiet" },
    { PROBIT,   OPT_R, "robust" },
    { PROBIT,   OPT_V, "verbose" },
    { QUANTREG, OPT_I, "intervals" },
    { QUANTREG, OPT_N, "no-df-corr" },
    { QUANTREG, OPT_R, "robust" },
    { QUIT,     OPT_X, "exit" },
    { RESET,    OPT_C, "cubes-only" },
    { RESET,    OPT_Q, "quiet" },
    { RESET,    OPT_R, "squares-only" },
    { RESTRICT, OPT_B, "bootstrap" },
    { RESTRICT, OPT_F, "full" },
    { RESTRICT, OPT_J, "jitter" },
    { RESTRICT, OPT_Q, "quiet" },
    { RESTRICT, OPT_V, "verbose" },
    { RESTRICT, OPT_L, "lbfgs" },
    { RESTRICT, OPT_N, "no-scaling" },
    { RUNS,     OPT_D, "difference" },
    { RUNS,     OPT_E, "equal" },
    { SCATTERS, OPT_L, "with-lines" },
    { SETINFO,  OPT_C, "continuous" },
    { SETINFO,  OPT_D, "discrete" },
    { SETOBS,   OPT_C, "stacked-cross-section" },
    { SETOBS,   OPT_P, "panel-vars" },
    { SETOBS,   OPT_R, "restructure" },
    { SETOBS,   OPT_S, "stacked-time-series" },
    { SETOBS,   OPT_T, "time-series" },
    { SETOBS,   OPT_X, "cross-section" },
    { SETOBS,   OPT_N, "special-time-series" },
    { SMPL,     OPT_F, "full" },
    { SMPL,     OPT_O, "dummy" },
    { SMPL,     OPT_M, "no-missing" },
    { SMPL,     OPT_N, "random" },
    { SMPL,     OPT_P, "replace" }, 
    { SMPL,     OPT_R, "restrict" },
    { SPEARMAN, OPT_V, "verbose" },
    { SQUARE,   OPT_O, "cross" },
    { STORE,    OPT_C, "csv" },
    { STORE,    OPT_D, "database" },
    { STORE,    OPT_F, "overwrite" },
    { STORE,    OPT_G, "dat" },
    { STORE,    OPT_J, "jmulti" },
    { STORE,    OPT_M, "gnu-octave" },
    { STORE,    OPT_R, "gnu-R" },
    { STORE,    OPT_T, "traditional" },
    { STORE,    OPT_Z, "gzipped" },
    { STORE,    OPT_X, "omit-obs" },
    { SYSTEM,   OPT_T, "iterate" },
    { SYSTEM,   OPT_V, "verbose" },
    { TOBIT,    OPT_V, "verbose" },
    { TSLS,     OPT_Q, "quiet" },
    { TSLS,     OPT_R, "robust" },  
    { TSLS,     OPT_S, "save" },
    { VAR,      OPT_D, "seasonals" },
    { VAR,      OPT_F, "variance-decomp" },
    { VAR,      OPT_I, "impulse-responses" },
    { VAR,      OPT_L, "lagselect" },
    { VAR,      OPT_N, "nc" },
    { VAR,      OPT_Q, "quiet" }, 
    { VAR,      OPT_R, "robust" }, 
    { VAR,      OPT_T, "trend" }, 
    { VECM,     OPT_A, "crt" },
    { VECM,     OPT_D, "seasonals" },
    { VECM,     OPT_F, "variance-decomp" },
    { VECM,     OPT_I, "impulse-responses" },
    { VECM,     OPT_N, "nc" },
    { VECM,     OPT_Q, "quiet" },
    { VECM,     OPT_R, "rc" },
    { VECM,     OPT_T, "ct" },
    { VECM,     OPT_V, "verbose" },
    { WLS,      OPT_R, "robust" },    
    { WLS,      OPT_Q, "quiet" },
    { XCORRGM,  OPT_Q, "quiet" },
    { XTAB,     OPT_C, "column" },
    { XTAB,     OPT_R, "row" },
    { XTAB,     OPT_Z, "zeros" },
    { 0,        0L,    NULL }
};

static int compare_strings (const void *a, const void *b)
{
    const char **sa = (const char **) a;
    const char **sb = (const char **) b;
     
    return strcmp(*sa, *sb);
}

char **get_all_option_strings (int *pn)
{
    char **optstrs;
    int i, j, m, n = 0;

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
	qsort(optstrs, n, sizeof *optstrs, compare_strings);
	m = n;
	for (i=0; i<m-1; i++) {
	    if (!strcmp(optstrs[i], optstrs[i+1])) {
		free(optstrs[i+1]);
		for (j=i+1; j<m-1; j++) {
		    optstrs[j] = optstrs[j+1];
		}
		optstrs[m-1] = NULL;
		i--;
		m--;
	    }
	}
	if (m < n) {
	    optstrs = realloc(optstrs, m * sizeof *optstrs);
	}
	*pn = m;
    }
    
    return optstrs;
}

const char **get_opts_for_command (int ci, int *nopt)
{
    int i, j, n = 0;
    const char **ret = NULL;

    if (vcv_opt_ok(ci) && ci != OLS) {
	n++; /* vcv */
    }

    for (i=0; gretl_opts[i].ci != 0; i++) {
	if (gretl_opts[i].ci == ci) n++;
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

    if (vcv_opt_ok(ci) && ci != OLS) {
	ret[j++] = "vcv";
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
    { OPT_Z, 'z' },
    { 0L,   '\0' }
};

static const char *ok_flags = "abcdefghijklmnopqrstuvwxz";

#define isflag(c) (c && (strchr(ok_flags, c) != NULL))

gretlopt opt_from_flag (unsigned char c)
{
    int i;

    for (i=0; flag_matches[i].c != '\0'; i++) {
	if (c == flag_matches[i].c) return flag_matches[i].o;
    }

    return 0L;
}

static int opt_is_valid (gretlopt opt, int ci, char c)
{
    int i;

    if (opt == OPT_O && vcv_opt_ok(ci)) {
	return 1;
    }

    for (i=0; gretl_opts[i].ci != 0; i++) {
	if (ci == gretl_opts[i].ci && opt == gretl_opts[i].o) {
	    return 1;
	}
    }

    if (c != 0) {
	sprintf(gretl_errmsg, "Invalid option '-%c'", c);
    } 

    return 0;
}

/* See if at point "p" (at which we've found '-') in string "s" we
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

static int is_long_opt (const char *lopt)
{
    int i, ret = 0;

    for (i=0; gretl_opts[i].o != 0; i++) {
	if (!strcmp(lopt, gretl_opts[i].longopt)) {
	    ret = 1;
	    break;
	}
    }

    return ret;
}

static int valid_long_opt (int ci, const char *lopt)
{
    int opt = OPT_NONE;
    int i;

    if (vcv_opt_ok(ci) && !strcmp(lopt, "vcv")) {
	return OPT_O;
    }

    /* start by looking for an exact match */
    for (i=0; gretl_opts[i].o != 0; i++) {
	if (ci == gretl_opts[i].ci && 
	    !strcmp(lopt, gretl_opts[i].longopt)) {
	    opt = gretl_opts[i].o;
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
		break;
	    }
	} 
    }   

    return opt;
}

/* The following apparatus is rather funky, and will probably have to
   be completely redone at some point.  For now, it allows passing a
   numerical value (two-sided p-value) in connection with the --auto
   option to the omit command.  We may want to make this mechanism
   available in connection with other option flags for other commands:
   in that case a redesign will be needed.  AC, 2008-01-08.
*/

static int optval_ci;
static gretlopt optval_opt;
static double optval_double = NADBL;
static char optval_str[32];

double get_optval_double (int ci, gretlopt opt)
{
    double x = NADBL;

    if (ci == optval_ci && optval_opt == opt) {
	x = optval_double;
    }

    optval_double = NADBL;

    return x;
}

void set_optval_double (int ci, gretlopt opt, double x)
{
    optval_ci = ci;
    optval_opt = opt;
    optval_double = x;
    gretl_push_c_numeric_locale();
    sprintf(optval_str, "%g", x);
    gretl_pop_c_numeric_locale();
}

static void unset_optval_double (void)
{
    optval_ci = 0;
    optval_opt = OPT_NONE;
    optval_double = NADBL;
    *optval_str = '\0';
}

static int valid_optval (int ci, gretlopt opt, const char *val)
{
    double x;

    if ((ci == OMIT && opt == OPT_A) ||
	(ci == QUANTREG && opt == OPT_I)) {
	if (numeric_string(val)) {
	    x = dot_atof(val);
	    if (x > 0.0 && x < 1.0) {
		set_optval_double(ci, opt, x);
		return 1;
	    } else {
		unset_optval_double();
	    }
	}
    }

    return 0;
}
  
static gretlopt get_long_opts (char *line, int ci, int *err)
{
    char *s = line;
    char longopt[32];
    char optval[32];
    gretlopt match, ret = 0L;

    while ((s = strstr(s, "--")) != NULL) {
	match = 0;
	*longopt = '\0';
	if (maybe_opt_start(line, s)) {
	    sscanf(s + 2, "%31[^ =]", longopt);
	    /* sscanf(s + 2, "%31s", longopt); */
	    match = valid_long_opt(ci, longopt);
	    if (match > 0) {
		/* recognized an acceptable option flag */
		ret |= match;
	    } else if (is_long_opt(longopt)) {
		/* recognized option, but not valid for the command */
		sprintf(gretl_errmsg, "Invalid option '--%s'", longopt);
		fprintf(stderr, " line='%s', ci = %d\n", line, ci);
		*err = 1;
		return 0L;
	    } 
	}

	if (match > 0) {
	    gretl_delete(s, 0, 2 + strlen(longopt));
	    if (*s == '=' && sscanf(s + 1, "%31[^ =]", optval) == 1) {
		if (valid_optval(ci, match, optval)) {
		    gretl_delete(s, 0, 1 + strlen(optval));
		}
	    }
	} else {
	    s += 2;
	}
    }

    return ret;
}

static void get_cmdword (const char *line, char *word)
{
    if (!sscanf(line, "%*s <- %8s", word)) {
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

#define ar1_alias(s) (!strcmp(s, "corc") || \
		      !strcmp(s, "hilu") || \
		      !strcmp(s, "pwe"))

#define smpl_alias(s) (!strcmp(s, "sample"))

/**
 * get_gretl_options:
 * @line: command line to parse.
 * @err: location for error code, which is set to 1 in case any 
 * invalid options are found, else set to 0.
 * 
 * Check for option flags in @line: if found, chop them out and set
 * the return value accordingly. Strip any trailing semicolon from
 * @line while we're at it.
 *
 * Returns: the options found in the line.
 */

gretlopt get_gretl_options (char *line, int *err)
{
    gretlopt oflags = 0L;
    int n = strlen(line);
    gretlopt opt;
    char cmdword[9] = {0};
    int ci, myerr = 0;

    gretl_error_clear();

    if (err != NULL) {
	*err = 0;
    }

    if (n < 2 || *line == '#') {
	return oflags;
    }

    /* to enable reading of trad. ESL input files */
    if (line[n-2] == ';' && isspace(line[n-1])) {
	line[n-2] = '\0';
    } else if (line[n-1] == ';') {
	line[n-1] = '\0';
    }

    get_cmdword(line, cmdword);

    if (strstr(line, "end nls")) {
	ci = NLS;
    } else if (strstr(line, "end mle")) {
	ci = MLE;
    } else if (strstr(line, "end gmm")) {
	ci = GMM;
    } else if (strstr(line, "end restrict")) {
	ci = RESTRICT;
    } else if (ar1_alias(cmdword)) {
	ci = AR1;
    } else if (smpl_alias(cmdword)) {
	ci = SMPL;
    } else {
	ci = gretl_command_number(cmdword);
    }

    /* some commands do not take a "flag", and "-%c" may have
       some other meaning */
    if (ci == 0 || ci == GENR || ci == PRINTF) {
	return oflags;
    }

    /* some commands are basically aliases */
    if (ci == OMITFROM) {
	ci = OMIT;
    } else if (ci == ADDTO) {
	ci = ADD;
    }

    /* smpl: in some contexts options don't make sense */
    if (ci == SMPL && strchr(line, ';')) {
	return oflags;
    }

    if (ci != SETINFO && ci != TABPRINT && ci != EQNPRINT) {
	/* try for short-form options (e.g. "-o") */
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

/**
 * print_flags:
 * @oflags: options.
 * @ci: command index, for context.
 * 
 * Constructs a string representation of the options in @oflags.
 *
 * Returns: pointer to static string (do not free!).
 */

const char *print_flags (gretlopt oflags, int ci)
{
    static char flagstr[512];
    char fbit[32];
    int i;

    flagstr[0] = '\0';

    if (oflags == OPT_NONE) {
	return flagstr;
    }

    if (ci == QUIT || ci == GENR) {
	/* any option flags are "hidden" */
	return flagstr;
    }

    if (ci == OMITFROM) {
	ci = OMIT;
    } else if (ci == ADDTO) {
	ci = ADD;
    }

    /* special: -o (--vcv) can be used with several model
       commands */
    if ((oflags & OPT_O) && vcv_opt_ok(ci)) {
	strcat(flagstr, " --vcv");
	oflags &= ~OPT_O;
    }

    for (i=0; gretl_opts[i].ci != 0; i++) {
	if (ci == gretl_opts[i].ci && (oflags & gretl_opts[i].o)) {
	    sprintf(fbit, " --%s", gretl_opts[i].longopt);
	    strcat(flagstr, fbit);
	    if (*optval_str != '\0' && gretl_opts[i].o == optval_opt 
		&& ci == optval_ci) {
		sprintf(fbit, "=%s", optval_str);
		strcat(flagstr, fbit);
		*optval_str = '\0';
	    }
	}
    }

    return flagstr;
}

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
