/*
 *  Copyright (c) 2004 by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "libgretl.h"
#include "gretl_private.h"

#define is_model_ci(c) (c == OLS || c == CORC || c == HILU || \
                        c == WLS || c == PWE || c == POOLED || c == HCCM || \
                        c == HSK || c == ADD || c == LAD || \
                        c == OMIT || c == TSLS || c == LOGIT || \
                        c == PROBIT || c == TOBIT || c == ARMA || \
                        c == AR || c == LOGISTIC || c == NLS || \
                        c == GARCH)

struct gretl_option {
    int ci;              /* command index (context) */
    gretlopt o;          /* index of integer type */
    const char *longopt; /* -- string representation of option */
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
    { ADD,      OPT_Q, "quiet" },
    { ADDTO,    OPT_Q, "quiet" },
    { ADF,      OPT_N, "nc" }, 
    { ADF,      OPT_C, "c" }, 
    { ADF,      OPT_R, "ctt" },     
    { ADF,      OPT_T, "ct" }, 
    { ADF,      OPT_V, "verbose" },
    { ARMA,     OPT_N, "native" },
    { ARMA,     OPT_V, "verbose" },
    { ARMA,     OPT_X, "x-12-arima" },
    { BXPLOT,   OPT_O, "notches" },
    { COINT,    OPT_N, "nc" },
    { COINT2,   OPT_O, "verbose" },
    { EQNPRINT, OPT_O, "complete" },
    { TABPRINT, OPT_O, "complete" },
    { ESTIMATE, OPT_M, "geomean" },
    { ESTIMATE, OPT_N, "no-df-corr" },
    { ESTIMATE, OPT_T, "iterate" },
    { FCASTERR, OPT_O, "plot" },
    { FREQ,     OPT_O, "gamma" },
    { FREQ,     OPT_Q, "quiet" },
    { GARCH,    OPT_A, "arma-init" },    
    { GARCH,    OPT_R, "robust" },
    { GARCH,    OPT_V, "verbose" },
    { GNUPLOT,  OPT_O, "with-lines" },
    { GNUPLOT,  OPT_M, "with-impulses" },
    { GNUPLOT,  OPT_S, "suppress-fitted" },
    { GNUPLOT,  OPT_Z, "dummy" },
    { GRAPH,    OPT_O, "tall" },
    { IMPORT,   OPT_O, "box1" },
    { KPSS,     OPT_T, "trend" },
    { KPSS,     OPT_V, "verbose" },
    { LEVERAGE, OPT_O, "save" },
    { LMTEST,   OPT_L, "logs" },
    { LMTEST,   OPT_M, "autocorr" },
    { LMTEST,   OPT_O, "autocorr" },
    { LMTEST,   OPT_S, "squares" }, 
    { LMTEST,   OPT_P, "panel" },
    { LMTEST,   OPT_W, "white" },
    { LOOP,     OPT_P, "progressive" },
    { LOOP,     OPT_V, "verbose" },
    { MEANTEST, OPT_O, "unequal-vars" },
    { OLS,      OPT_N, "no-df-corr" },
    { OLS,      OPT_O, "vcv" }, 
    { OLS,      OPT_P, "print-final" },
    { OLS,      OPT_R, "robust" },
    { OLS,      OPT_Q, "quiet" },
    { OMIT,     OPT_Q, "quiet" },
    { OMITFROM, OPT_Q, "quiet" },
    { OUTFILE,  OPT_A, "append" },
    { OUTFILE,  OPT_C, "close" },
    { OUTFILE,  OPT_W, "write" },
    { PANEL,    OPT_C, "cross-section" },
    { PANEL,    OPT_S, "time-series" },
    { POOLED,   OPT_T, "iterate" },
    { POOLED,   OPT_V, "verbose" },
    { POOLED,   OPT_W, "unit-weights" },
    { PCA,      OPT_A, "save-all" },
    { PCA,      OPT_O, "save" },
    { PERGM,    OPT_O, "bartlett" },
    { PLOT,     OPT_O, "one-scale" },
    { PRINT,    OPT_O, "byobs" },
    { PRINT,    OPT_T, "ten" },
    { SETOBS,   OPT_C, "stacked-cross-section" },
    { SETOBS,   OPT_S, "stacked-time-series" },
    { SETOBS,   OPT_T, "time-series" },
    { SETOBS,   OPT_X, "cross-section" },
    { SETOBS,   OPT_N, "special-time-series" },
    { SMPL,     OPT_C, "replace" },    
    { SMPL,     OPT_O, "dummy" },
    { SMPL,     OPT_M, "no-missing" },
    { SMPL,     OPT_N, "random" },
    { SMPL,     OPT_R, "restrict" },
    { SPEARMAN, OPT_O, "verbose" },
    { SQUARE,   OPT_O, "cross" },
    { STORE,    OPT_C, "csv" },
    { STORE,    OPT_D, "dat" },
    { STORE,    OPT_M, "gnu-octave" },
    { STORE,    OPT_R, "gnu-R" },
    { STORE,    OPT_T, "traditional" },
    { STORE,    OPT_Z, "gzipped" },
    { TOBIT,    OPT_V, "verbose" },
    { TSLS,     OPT_R, "robust" },  
    { TSLS,     OPT_S, "save" },
    { VAR,      OPT_V, "impulse-responses" },
    { VAR,      OPT_Q, "quiet" },  
    { VAR,      OPT_R, "robust" },  
    { 0,        0L,    NULL }
};

const char **get_opts_for_command (int ci, int *nopt)
{
    int i, j, n = 0;
    const char **ret = NULL;

    if (is_model_ci(ci) && ci != OLS && ci != LAD) {
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

    if (is_model_ci(ci) && ci != OLS && ci != LAD) {
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
    { OPT_I, 'i' },
    { OPT_L, 'l' },
    { OPT_M, 'm' },
    { OPT_N, 'n' },
    { OPT_O, 'o' },
    { OPT_P, 'p' },
    { OPT_Q, 'q' },
    { OPT_R, 'r' },
    { OPT_S, 's' },
    { OPT_T, 't' },
    { OPT_V, 'v' },
    { OPT_W, 'w' },
    { OPT_X, 'x' },
    { OPT_Z, 'z' },
    { 0L,   '\0' }
};

/* note: 'f' is not treated as an option flag for now */

#define isflag(c) (c == 'a' || c == 'b' || c == 'c' || c == 'd' || \
                   c == 'i' || c == 'l' || c == 'm' || c == 'p' || \
                   c == 'n' || c == 'o' || c == 'q' || c == 'r' || \
                   c == 's' || c == 't' || c == 'v' || c == 'w' || \
                   c == 'x' || c == 'z')

static gretlopt opt_from_flag (unsigned char c)
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

    if (opt == OPT_O && is_model_ci(ci)) {
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

static gretlopt get_short_opts (char *line, int ci, int *err)
{
    char *p = strchr(line, '-');
    gretlopt opt, ret = 0L;

    while (p != NULL) {
	unsigned char c, prev;
	int match = 0;
	size_t n = strlen(p);

	c = *(p + 1);
	prev = *(p - 1);
	
	if (isspace(prev) && isflag(c) && (n == 2 || isspace(*(p + 2)))) {
	    opt = opt_from_flag(c);
	    if (!opt_is_valid(opt, ci, c)) {
		*err = 1;
		return 0L;
	    }
	    ret |= opt;
	    gretl_delete(p, 0, 2);
	    match = 1;
	}
	if (!match) p++;
	p = strchr(p, '-');
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

    if (is_model_ci(ci) && ci != LAD && !strcmp(lopt, "vcv")) {
	/* VCV is available for all but LAD models */
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
  
static gretlopt get_long_opts (char *line, int ci, int *err)
{
    char *p = strstr(line, "--");
    gretlopt match, ret = 0L;

    while (p != NULL) {
	char longopt[32];

	sscanf(p + 2, "%31s", longopt);
	match = valid_long_opt(ci, longopt);

	if (match > 0) {
	    ret |= match;
	    gretl_delete(p, 0, 2 + strlen(longopt));
	} else if (is_long_opt(longopt)) {
	    /* recognized option, but not valid for the command */
	    sprintf(gretl_errmsg, "Invalid option '--%s'", longopt);
	    *err = 1;
	    return 0L;
	} else {
	    p += 2;
	}
	p = strstr(p, "--");
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

/**
 * get_gretl_options:
 * @line: command line to parse.
 * @err: address for error code, which is set to 1 in case any 
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

    *gretl_errmsg = '\0';

    if (err != NULL) {
	*err = 0;
    }

    if (n < 2 || *line == '#') return oflags;

    /* to enable reading of trad. esl input files */
    if (line[n-2] == ';' && isspace(line[n-1])) {
	line[n-2] = '\0';
    } else if (line[n-1] == ';') {
	line[n-1] = '\0';
    }

    /* some commands do not take a "flag", and "-%c" may have
       some other meaning */
    get_cmdword(line, cmdword);

    if (!strcmp(cmdword, "genr") || !strcmp(cmdword, "sim") ||
	!strcmp(cmdword, "label")) return oflags;

    if (strstr(line, "end nls")) {
	ci = NLS;
    } else {
	ci = gretl_command_number(cmdword);
    }

    if (ci == 0) return oflags;

    /* try for short-form options (e.g. "-o") */
    opt = get_short_opts(line, ci, &myerr);
    if (!myerr && opt) {
	oflags |= opt;
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
    static char flagstr[64];
    char fbit[20];
    int i;

    flagstr[0] = '\0';

    if (oflags == 0L) return flagstr;

    /* special: -o (--vcv) can be used with several model
       commands */
    if ((oflags & OPT_O) && is_model_ci(ci)) {
	strcat(flagstr, " --vcv");
	oflags &= ~OPT_O;
    }

    for (i=0; gretl_opts[i].ci != 0; i++) {
	if (ci == gretl_opts[i].ci && (oflags & gretl_opts[i].o)) {
	    sprintf(fbit, " --%s", gretl_opts[i].longopt);
	    strcat(flagstr, fbit);
	}
    }

    return flagstr;
}

int check_for_loop_only_options (int ci, gretlopt opt, PRN *prn)
{
    if (ci == OLS && (opt & OPT_P)) {
	const char *flagstr = print_flags(OPT_P, OLS);

	pprintf(prn, _("Warning: option%s ignored outside of loop"), 
		flagstr);
	pputc(prn, '\n');
	return 1;
    }

    return 0;
}
