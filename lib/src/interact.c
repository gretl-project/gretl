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

/* interact.c for gretl */

#include "libgretl.h"
#include "monte_carlo.h"
#include "var.h"
#include "johansen.h"
#include "gretl_func.h"
#include "compat.h"
#include "system.h"
#include "forecast.h"
#include "cmd_private.h"
#include "libset.h"
#include "usermat.h"
#include "modelspec.h"
#include "gretl_panel.h"
#include "texprint.h"
#include "gretl_xml.h"
#include "gretl_string_table.h"
#include "dbread.h"

#include <glib.h>

/* for the "shell" command */
#ifndef WIN32
# ifdef HAVE_PATHS_H
#  include <paths.h>
# endif
#endif

#define CMD_DEBUG 0
#define ARMA_DBG 0

#include "laginfo.c"

typedef struct {
    int v;
    char vname[VNAMELEN];
    int lmin;
    int lmax;
    int *laglist;
} LAGVAR;

#define cmd_set_nolist(c) (c->flags |= CMD_NOLIST)
#define cmd_unset_nolist(c) (c->flags &= ~CMD_NOLIST)

static void get_optional_filename_etc (const char *line, CMD *cmd);
static int if_eval (const char *s, double ***pZ, DATAINFO *pdinfo, int *err);
static int set_if_state (int code);
static int get_if_state (int code);

#define bare_quote(p,s)   (*p == '"' && (p-s==0 || *(p-1) != '\\'))
#define starts_comment(p) (*p == '/' && *(p+1) == '*')
#define ends_comment(p)   (*p == '*' && *(p+1) == '/')

static int strip_inline_comments (char *s)
{
    int ret = 0;

    if (*s == '#') {
	/* the entire line is a comment */
	ret = 1;
    } else if (strstr(s, "#")) {
	int quoted = 0;
	int braced = 0;
	char *p = s;

	while (*p) {
	    if (bare_quote(p, s)) {
		quoted = !quoted;
	    } else if (!quoted) {
		if (*p == '{') {
		    braced++;
		} else if (*p == '}') {
		    braced--;
		}
	    }
	    if (!quoted && !braced) {
		if (*p == '#') {
		    *p = '\0';
		    break;
		}
	    }
	    p++;
	}
    }

    return ret;
}

/* filter_comments: strip comments out of line; return non-zero if
   the whole line is a comment */

static int filter_comments (char *s, CMD *cmd)
{
    char tmp[MAXLINE];
    char *p = s;
    int quoted = 0;
    int ignore = (cmd->flags & CMD_IGNORE);
    int j = 0, filt = 0;

    if (strlen(s) >= MAXLINE) {
	cmd->err = E_TOOLONG;
	return 0;
    }

    while (*p) {
	if (!quoted && !ignore && *p == '#') {
	    break;
	}
	if (!ignore && bare_quote(p, s)) {
	    quoted = !quoted;
	}
	if (!quoted) {
	    if (starts_comment(p)) {
		ignore = 1;
		p += 2;
	    } else if (ends_comment(p)) {
		if (!ignore) {
		    cmd->err = E_PARSE;
		    return 0;
		}
		ignore = 0;
		p += 2;
		p += strspn(p, " ");
	    }
	}
	if (!ignore && *p != '\r') {
	    tmp[j++] = *p;
	}
	if (*p) {
	    p++;
	}
    }

    tmp[j] = '\0';
    strcpy(s, tmp);
    tailstrip(s);

    if (*s == '\0') { 
	filt = 1;
    } else if (!ignore) {
	/* '#' comments */
	filt = strip_inline_comments(s);
	tailstrip(s);
    }

    if (filt) {
	/* the whole line is a comment */
	cmd_set_nolist(cmd);
	cmd->ci = CMD_COMMENT;
    }

    if (ignore) {
	/* the line ends in multi-line comment mode */
	cmd->flags |= CMD_IGNORE;
    } else {
	cmd->flags &= ~CMD_IGNORE;
    }

    return filt;
}

static int get_rhodiff_or_lags_param (char *s, CMD *cmd)
{
    int k = haschar(';', s);
    int ret = 0;

    if (k > 0) {
	free(cmd->param);
	cmd->param = gretl_strndup(s, k);
	shift_string_left(s, k + 1);
	ret = 1;
    }

    return ret;
}

static void cmd_set_param_direct (CMD *cmd, const char *s)
{
    free(cmd->param);
    cmd->param = gretl_strdup(s);
    if (cmd->param == NULL) {
	cmd->err = E_ALLOC;
    }
}

/* catch aliased command words and assign ci; return
   ci if alias caught, else 0. */

static int catch_command_alias (char *line, CMD *cmd)
{
    char *s = cmd->word;

    cmd->ci = 0;

    if (!strcmp(line, "exit")) {
	strcpy(s, "quit");
	cmd->ci = QUIT;
	cmd->opt = OPT_X;
    } else if (!strcmp(s, "ls")) {
	cmd->ci = VARLIST;
    } else if (!strcmp(s, "boxplots")) { 
	cmd->ci = BXPLOT;
    } else if (!strcmp(s, "man")) {
	cmd->ci = HELP;
    } else if (!strcmp(s, "pooled")) {
	cmd->ci = OLS;
    } else if (!strcmp(s, "import")) {
	cmd->ci = OPEN;
    } else if (!strcmp(line, "smpl full")) {
	strcpy(line, "smpl");
	cmd->opt = OPT_F;
    } else if (!strcmp(s, "sample")) {
	cmd->ci = SMPL;
    } else if (!strcmp(s, "list")) {
	char lname[VNAMELEN];

	if (gretl_string_ends_with(line, "print")) {
	    if (sscanf(line, "list %15s", lname)) {
		strcpy(line, lname);
	    }
	} 
	cmd->ci = GENR;
    } else if (*s == '!' || !strcmp(s, "launch")) {
	cmd->ci = SHELL;
    } else if (!strcmp(s, "funcerr")) {
	cmd->ci = FUNCERR;
    } else if (!strcmp(line, "end if")) {
	strcpy(s, "endif");
	strcpy(line, "endif");
	cmd->ci = ENDIF;
    } else if (!strcmp(s, "elif")) {
	cmd->ci = ELSE;
	cmd->opt = OPT_I;
    } else if (!strcmp(s, "addobs")) {
	cmd->ci = DATAMOD;
	cmd_set_param_direct(cmd, "addobs");
    } else if (!strcmp(s, "transpos")) {
	cmd->ci = DATAMOD;
	cmd_set_param_direct(cmd, "transpose");
    } else if (!strcmp(s, "fit")) {
	cmd->ci = FCAST;
	cmd->opt |= OPT_F;
	strcpy(line, "fcast autofit");
    } else if (!strcmp(s, "fcasterr")) {
	cmd->ci = FCAST;
	cmd->opt |= OPT_E;
    } else if (!strcmp(s, "corc")) {
	strcpy(s, "ar1");
	cmd->ci = AR1;
    } else if (!strcmp(s, "hilu")) {
	strcpy(s, "ar1");
	cmd->ci = AR1;
	cmd->opt |= OPT_H;
    } else if (!strcmp(s, "pwe")) {
	strcpy(s, "ar1");
	cmd->ci = AR1;
	cmd->opt |= OPT_P;
    }	

    return cmd->ci;
}

#define REQUIRES_PARAM(c) (c == ADDTO || \
                           c == DATAMOD || \
                           c == FUNC || \
                           c == LOOP ||  \
                           c == NULLDATA || \
                           c == OMITFROM || \
                           c == SETMISS)

#define REQUIRES_ORDER(c) (c == ADF || \
                           c == ARCH || \
                           c == COINT || \
                           c == COINT2 || \
                           c == KPSS || \
                           c == VAR || \
                           c == VECM)

#define NO_VARLIST(c) (c == APPEND || \
                       c == BREAK || \
                       c == CHOW || \
	               c == CRITERIA || \
	               c == CUSUM || \
                       c == DATA || \
                       c == END || \
	               c == ENDLOOP || \
                       c == ESTIMATE || \
	               c == EQNPRINT || \
	               c == FCAST || \
                       c == FUNC || \
                       c == FUNCERR || \
	               c == GENR || \
                       c == GMM || \
	               c == HAUSMAN || \
                       c == HELP || \
                       c == INCLUDE || \
    	               c == INFO || \
 	               c == LABELS || \
                       c == LEVERAGE || \
                       c == LMTEST || \
                       c == LOOP || \
                       c == MLE || \
                       c == MODELTAB || \
                       c == NLS || \
                       c == NULLDATA || \
		       c == OPEN || \
                       c == OUTFILE || \
                       c == PRINTF || \
	               c == PVALUE || \
                       c == QLRTEST || \
	               c == QUIT || \
                       c == RENAME || \
                       c == RESET || \
                       c == RESTRICT || \
	               c == RUN || \
                       c == SET || \
                       c == SETINFO || \
	               c == SETOBS || \
	               c == SHELL || \
                       c == SPRINTF || \
		       c == SSCANF || \
                       c == SYSTEM || \
                       c == TABPRINT || \
                       c == TESTUHAT || \
                       c == VARLIST || \
                       c == VIF)

#define USES_LISTSEP(c) (c == AR || \
                         c == ARBOND || \
                         c == ARMA || \
                         c == COINT2 || \
                         c == EQUATION || \
                         c == HECKIT || \
                         c == GARCH || \
                         c == MPOLS || \
                         c == POISSON || \
                         c == PRINT || \
                         c == SCATTERS || \
                         c == TSLS || \
                         c == VAR || \
                         c == VECM || \
                         c == XTAB)

#define DOUBLE_SEP_OK(c) (c == ARBOND || \
                          c == ARMA || \
                          c == COINT2 || \
			  c == VECM) 

#define NEEDS_LISTSEP(c) (c == AR || \
                          c == ARBOND || \
                          c == ARMA || \
                          c == HECKIT || \
                          c == GARCH || \
                          c == TSLS)

#define DEFAULTS_TO_FULL_LIST(c) (c == CORR || \
                                  c == DIFF || \
                                  c == LDIFF || \
                                  c == LAGS || \
                                  c == LOGS || \
                                  c == PCA || \
                                  c == PRINT || \
                                  c == SDIFF || \
                                  c == SMPL || \
                                  c == SQUARE || \
                                  c == ORTHDEV || \
                                  c == STORE || \
                                  c == SUMMARY)

#define SCALARS_OK_IN_LIST(c) (c == DELEET || \
                               c == PRINT || \
                               c == STORE)

#define MODIFIES_LIST(c) (c == DIFF || \
			  c == DUMMIFY || \
			  c == LDIFF || \
			  c == SDIFF || \
			  c == LAGS || \
			  c == LOGS || \
			  c == SQUARE || \
	                  c == ORTHDEV)

enum {
    SET_FALSE,
    SET_TRUE,
    SET_ELSE,
    SET_ELIF,
    SET_ENDIF,
    IS_FALSE,
    IS_TRUE,
    UNINDENT,
    GETINDENT,
    RELAX
};

#define IFDEBUG 0

static int flow_control (const char *line, double ***pZ, 
			 DATAINFO *pdinfo, CMD *cmd)
{
    int ci = cmd->ci;
    int blocked, ok, err = 0;

    blocked = get_if_state(IS_FALSE);

    if (ci != IF && ci != ELSE && ci != ENDIF) {
	return blocked;
    }

    if (ci == IF) {
	if (blocked) {
	    err = set_if_state(SET_FALSE);
	} else {
	    ok = if_eval(line, pZ, pdinfo, &err);
	    if (!err) {
		err = set_if_state(ok? SET_TRUE : SET_FALSE);
	    }
	}
    } else if (ci == ENDIF) {
	err = set_if_state(SET_ENDIF);
    } else if (ci == ELSE && (cmd->opt & OPT_I)) {
	err = set_if_state(SET_ELIF);
	if (!err && get_if_state(IS_TRUE)) {
	    set_if_state(UNINDENT);
	    ok = if_eval(line, pZ, pdinfo, &err);
	    if (!err) {
		err = set_if_state(ok? SET_TRUE : SET_FALSE);
	    }
	}
    } else if (ci == ELSE) {
	err = set_if_state(SET_ELSE);
    }

    if (err) {
	set_if_state(RELAX);
	cmd->err = err;
    } 

    return 1;
}

static char cmd_savename[MAXSAVENAME];

static void maybe_extract_savename (char *s, CMD *cmd)
{
    *cmd->savename = 0;

    if (strncmp(s, "genr ", 5) && strstr(s, "<-")) {
	int n, len, quoted;

	quoted = (*s == '"');
	len = strcspn(s, "<");
	if (len - quoted == 0) {
	    return;
	}
	n = len - quoted;
	if (n > MAXSAVENAME - 1) {
	    n = MAXSAVENAME - 1;
	}
	strncat(cmd->savename, s + quoted, n);
	tailstrip(cmd->savename);
	n = strlen(cmd->savename);
	if (cmd->savename[n-1] == '"') {
	    cmd->savename[n-1] = 0;
	}
	strcpy(cmd_savename, cmd->savename);
	shift_string_left(s, len + 2);
	while (*s == ' ') {
	    shift_string_left(s, 1);
	}
    }
}

/* grab a filename, possibly prepending userdir */

static int filename_to_param (CMD *cmd, char *s, int *len,
			      int *quoted, int *nsp)
{
    char *fname;

    while (isspace(*s)) s++;

    if (*s == '"' || *s == '\'') {
	char *p = strchr(s + 1, *s);

	if (p == NULL) {
	    return E_PARSE;
	}
	*len = p - s - 1;
	if (*len == 0) {
	    return E_PARSE;
	}
	*quoted = 1;
    } else {
	*len = strcspn(s, " ");
    }

    free(cmd->param);
    cmd->param = NULL;

    fname = gretl_strndup(s + *quoted, *len);
    if (fname == NULL) {
	return E_ALLOC;
    }

    if (nsp != NULL) {
	/* count spaces in supplied filename */
	int i;

	*nsp = 0;
	for (i=0; i<*len; i++) {
	    if (fname[i] == ' ') *nsp += 1;
	}
    }

    if (libset_get_bool(USE_CWD) || gretl_path_is_absolute(fname)) {
	cmd->param = fname;
    } else if (cmd->ci == OUTFILE && !strcmp(fname, "null")) {
	cmd->param = fname;
    } else {
	cmd->param = gretl_strdup_printf("%s%s", gretl_work_dir(), fname);
	free(fname);
	if (cmd->param == NULL) {
	    return E_ALLOC;
	}
    }

    return 0;
}

static int 
get_maybe_quoted_filename (CMD *cmd, char *s, int *nf)
{
    int err, len = 0;
    int quoted = 0;
    int nsp = 0;

    err = filename_to_param(cmd, s, &len, &quoted, &nsp);
    if (err) {
	return err;
    }

    if (nsp) {
	*nf -= nsp;
    }

    shift_string_left(s, len + 2 * quoted);

    return 0;
} 

static void grab_gnuplot_literal_block (char *s, CMD *cmd)
{
    s = strchr(s, '{');
    if (s != NULL) {
	free(cmd->param);
	cmd->param = gretl_strdup(s);
	*s = 0;
    }
}

static int arma_special_field (const char *s)
{
    int ret = 0;

    if (*s == '{') {
	ret = 1;
    } else if (isalpha(*s) && get_matrix_by_name(s)) {
	ret = 1;
    }

    return ret;
}

/* Check that list separator and any "{xxx}" or matrix fields 
   are in appropriate positions; also see if there _are_ any
   of the latter.
 */

static int arma_basic_check (char **S, int n, int *specials)
{
    int i, spos = 0;

    for (i=0; i<n; i++) {
	if (S[i][0] == ';') {
	    spos = i;
	    break;
	}
    }

    if (spos != 2 && spos != 3) {
	return E_PARSE;
    }

    for (i=0; i<n; i++) {
	if (arma_special_field(S[i])) {
	    if (i != 0 && i != spos - 1) {
		return E_PARSE;
	    } else if (i == 0) {
		specials[0] = 1;
	    } else {
		specials[1] = i;
	    }
	}
    }

    return 0;
}

/* push onto array *pS a string of length len starting at s */

static int push_field (char ***pS, const char *s, int len, int *nf)
{
    char *chunk;
    int err = 0;

    if (len == 0) {
	return 0;
    }

    chunk = gretl_strndup(s, len);
    if (chunk == NULL) {
	err = E_ALLOC;
    } else {
	err = strings_array_add(pS, nf, chunk);
	free(chunk);
    }

    return err;
}

#define no_specials(s) (s[0] == 0 && s[1] == 0)

/* split fields in ar(i)ma command, allowing for specific
   non-seasonal AR and/or MA lags given within '{' and '}',
   or specified via a named matrix.
*/

static char **arma_split_fields (const char *s, int *nf, 
				 int *specials, int *err)
{
    const char *q, *p = s;
    char **F = NULL;
    int n;

    while (*p && !*err) {
	while (*p == ' ') p++;
	if (*p == '{') {
	    q = strchr(p, '}');
	    if (q == NULL) {
		*err = E_PARSE;
	    } else {
		n = strcspn(p, "}");
		*err = push_field(&F, p, n + 1, nf);
		p = q;
	    }
	} else if (*p == ';') {
	    *err = push_field(&F, p, 1, nf);
	} else {
	    n = strcspn(p, " {};");
	    *err = push_field(&F, p, n, nf);
	    p += n - 1;
	}
	p++;
    }

    if (!*err) {
	*err = arma_basic_check(F, *nf, specials);
    }

    if (*err || no_specials(specials)) {
	free_strings_array(F, *nf);
	F = NULL;
    }

    return F;
}

/* matrix must be vector of integers */

static int max_lag_from_matrix (char *s, int *err)
{
    const gretl_matrix *m;
    int k, kmax = 0, kmin = 99999;
    int i, n = 0;

    m = get_matrix_by_name(s);

    if (m == NULL) {
	*err = E_UNKVAR;
    } else if ((n = gretl_vector_get_length(m)) == 0) {
	*err = E_NONCONF;
    }

    for (i=0; i<n && !*err; i++) {
	k = (int) m->val[i];
	if (k > kmax) kmax = k;
	if (k < kmin) kmin = k;
    }

    if (!*err && kmin < 1) {
	*err = E_DATA;
    }

    return kmax;
}

static int max_lag_from_numerics (char *s, int *err)
{
    int k, kmax = 0, kmin = 99999;
    char *rem;

    while (*s != '\0' && *s != '}') {
	k = (int) strtol(s, &rem, 10);
	if (k > kmax) kmax = k;
	if (k < kmin) kmin = k;
	s = rem;
    }

    if (kmin < 1) {
	*err = E_DATA;
    }

    return kmax;
}

/* e.g. "matrix" or "{1 4}" */

static int max_lag_from_field (char *s, int *err)
{
    int mat = 1, kmax = 0;

    if (*s == '{') {
	mat = 0;
	s++;
    }

    if (mat) {
	kmax = max_lag_from_matrix(s, err);
    } else {
	kmax = max_lag_from_numerics(s, err);
    } 

    return kmax;
}

static int arma_maybe_rewrite (char *s, CMD *cmd)
{
    char **S = NULL;
    char *orig = s;
    char chunk[16];
    int pq[2] = {0};
    int i, k, n = 0;

#if ARMA_DBG
    fprintf(stderr, "arma_rewrite: s = '%s'\n", s);
#endif

    if (!strncmp(s, "arma ", 5)) {
	s += 5;
    } else if (!strncmp(s, "arima ", 6)) {
	s += 6;
    }

    S = arma_split_fields(s, &n, pq, &cmd->err);
    if (S == NULL) {
	return cmd->err;
    }

#if ARMA_DBG
    for (i=0; i<n; i++) {
	fprintf(stderr, "S[%d] = '%s'\n", i, S[i]);
    }
#endif

    /* preserve original command line for echo */
    free(cmd->extra);
    cmd->extra = gretl_strdup(orig);

    cmd->param = malloc(strlen(s) + 1);
    if (cmd->param == NULL) {
	cmd->err = E_ALLOC;
    } else {
	*s = *cmd->param = '\0';
    }

    for (i=0; i<n && !cmd->err; i++) {
	if ((i == 0 && pq[0]) || (i > 0 && i == pq[1])) {
	    k = max_lag_from_field(S[i], &cmd->err);
	    if (!cmd->err) {
		sprintf(chunk, "%d", k);
		strcat(s, chunk);
		if (i == 0) {
		    strcat(cmd->param, "p=");
		    strcat(cmd->param, S[i]);
		} else {
		    if (*cmd->param) {
			strcat(cmd->param, ",");
		    }
		    strcat(cmd->param, "q=");
		    strcat(cmd->param, S[i]);
		}		    
	    }
	} else {
	    strcat(s, S[i]);
	}
	if (i < n - 1) {
	    strcat(s, " ");
	}
    } 

#if ARMA_DBG
    fprintf(stderr, "new s = '%s'\n", s);
    fprintf(stderr, "cmd->param = '%s'\n", cmd->param);
#endif

    free_strings_array(S, n);

    return cmd->err;
}

/* pluck the specification for "block-diagonal" instruments out of an
   arbond command line, and put it in the command's "param" field for
   subsequent special processing in arbond.c */

static void grab_arbond_diag (char *s, CMD *cmd)
{
    char *param = NULL;
    char *s0, *p, *q;
    int k;

    s0 = s = strrchr(s, ';');

    while ((s = strstr(s, "GMM("))) {
	p = strchr(s, ')');
	if (p == NULL) {
	    cmd->err = E_PARSE;
	} else {
	    p++;
	    k = p - s;
	    q = gretl_strndup(s, k);
	    param = gretl_str_expand(&param, q, " ");
	    if (param == NULL) {
		cmd->err = E_ALLOC;
	    }
	    free(q);
	    while (*p == ' ') {
		p++; k++;
	    }
	    shift_string_left(s, k);
	}
	if (cmd->err) {
	    break;
	}
    }

    if (param != NULL) {
	free(cmd->param);
	cmd->param = param;
	tailstrip(s0);
    }
}

#define LAG_DEBUG 0

static int lag_from_lstr (const char *s, 
			  const double **Z,
			  const DATAINFO *pdinfo,
			  int *err)
{
    int lsign = 1, lag = 0;

    *err = 0;

    if (!strcmp(s, "0")) {
	/* lag zero = contemp. value */
	lsign = 1;
    } else if (isalpha(*s)) {
	lsign = -1;
    } else if (*s == '-') {
	lsign = 1;
	s++;
    } else if (*s == '+') {
	lsign = -1;
	s++;
    } else {
	*err = 1;
    }

    if (!*err) {
	if (isdigit(*s)) {
	    lag = atoi(s);
	} else {
	    int v = varindex(pdinfo, s);

	    if (v < pdinfo->v) {
		lag = Z[v][0];
	    } else {
		*err = 1;
	    }
	}
    }

    if (!*err) {
	lag = lsign * lag;
    }

    return lag;
}

static int get_contiguous_lags (LAGVAR *lv,
				const char *l1str, const char *l2str,
				const double **Z, const DATAINFO *pdinfo)
{
    int err = 0;

    lv->lmin = lag_from_lstr(l1str, Z, pdinfo, &err);

    if (!err) {
	lv->lmax = lag_from_lstr(l2str, Z, pdinfo, &err);
    }

    return err;
}

static int parse_lagvar (const char *s, LAGVAR *lv, 
			 const double **Z, const DATAINFO *pdinfo)
{
    char l1str[16], l2str[16];
    int i, err = 1;

    lv->v = 0;
    *lv->vname = 0;
    lv->lmin = 0;
    lv->lmax = 0;
    lv->laglist = NULL;

    if (sscanf(s, "%15[^(](%15s", lv->vname, l1str) != 2) {
	return err;
    }

    lv->v = varindex(pdinfo, lv->vname);
    if (lv->v == 0 || lv->v >= pdinfo->v) {
	return err;
    }

    if (sscanf(s, "%15[^(](%15s to %15[^)])", lv->vname, 
	       l1str, l2str) == 3) {
	err = get_contiguous_lags(lv, l1str, l2str, Z, pdinfo);
    } else if (strchr(l1str, ',') != NULL) {
	lv->laglist = gretl_list_from_string(strchr(s, '(') + 1);
	if (lv->laglist != NULL) {
	    for (i=1; i<=lv->laglist[0]; i++) {
		lv->laglist[i] = -lv->laglist[i];
	    }
	    err = 0;
	}
    } else {
	sscanf(s, "%15[^(](%15[^ )]", lv->vname, l1str);
	lv->lmin = lv->lmax = lag_from_lstr(l1str, Z, pdinfo, &err);
    }

#if LAG_DEBUG
    fprintf(stderr, "parse_lagvar: s = '%s'\n", s);
    fprintf(stderr, " lmin = %d, lmax = %d\n",
	    lv->lmin, lv->lmax);
    if (lv->laglist != NULL) {
	printlist(lv->laglist, "lv->laglist");
    }
#endif

    return err;
}

static int cmd_full_list (const DATAINFO *pdinfo, CMD *cmd)
{
    int nv = 0, err = 0;
    int *list;

    if (cmd->ci == PRINT && cmd->extra != NULL &&
	*cmd->extra != '\0') {
	/* no-op */
	return 0;
    }

    if (cmd->flags & CMD_NULLIST) {
	/* no-op */
	cmd->flags ^= CMD_NULLIST;
	return 0;
    }

    list = full_var_list(pdinfo, &nv);

    if (list == NULL) {
	if (nv > 0) {
	    err = E_ALLOC;
	}
    } else {
	free(cmd->list);
	cmd->list = list;
    }

    return err;
}

static int expand_command_list (CMD *cmd, int add)
{
    int i, oldn = cmd->list[0];
    int *list;

    list = realloc(cmd->list, (oldn + add) * sizeof *list);
    if (list == NULL) {
	cmd->err = E_ALLOC;
	strcpy (gretl_errmsg, 
		_("Memory allocation failed for command list"));
	return 1;
    }

    /* one of the added vars was "already assumed" */
    list[0] += (add - 1);

    /* avoid uninitialized values */
    for (i=oldn+1; i<=list[0]; i++) {
	list[i] = 0;
    }
    
    cmd->list = list;

    return 0;
}

/* Get the total number of lags and set the increment for
   generating successive lags.  Allows for mixed leads
   and lags. */

static int get_n_lags (LAGVAR *lv, int *incr)
{
    int nl = 0;

    if (lv->laglist != NULL) {
	nl = lv->laglist[0];
	*incr = 0;
    } else if (lv->lmax >= lv->lmin) {
	nl = lv->lmax - lv->lmin + 1;
	*incr = 1;
    } else {
	nl = lv->lmin - lv->lmax + 1;
	*incr = -1;
    }

    return nl;
}

/* see if we have a valid specification for automatically adding
   lags of a given variable to the command list */

int auto_lag_ok (const char *s, int *lpos,
		 double ***pZ, DATAINFO *pdinfo,
		 CMD *cmd)
{
    LAGVAR lagvar;
    int nlags, i;
    int llen = *lpos;
    int lincr = 1;
    int ok = 1;
	
    if (parse_lagvar(s, &lagvar, (const double **) *pZ, pdinfo)) {
	ok = 0;
	goto bailout;
    }

    nlags = get_n_lags(&lagvar, &lincr);

#if LAG_DEBUG
    if (lagvar.laglist != NULL) {
	fprintf(stderr, "auto lags: n=%d, incr=%d\n", nlags, lincr);
    } else {
	fprintf(stderr, "auto lags: last=%d, first=%d, n=%d, incr=%d\n",
		lagvar.lmax, lagvar.lmin, nlags, lincr);
    }
#endif

    if (nlags <= 0) {
	cmd->err = E_PARSE;
	ok = 0;
	goto bailout;
    }

    if (nlags > 1 && expand_command_list(cmd, nlags)) {
	ok = 0;
	goto bailout;
    }

    for (i=0; i<nlags && ok; i++) {
	int order, lv;

	if (lagvar.laglist != NULL) {
	    order = lagvar.laglist[i+1];
	} else {
	    order = lagvar.lmin + i * lincr;
	}

	lv = laggenr(lagvar.v, order, pZ, pdinfo);

#if LAG_DEBUG
	fprintf(stderr, "laggenr for var %d (%s), lag %d, gave lv = %d\n",
		lagvar.v, pdinfo->varname[lagvar.v], order, lv);
#endif
	if (lv < 0) {
	    cmd->err = 1;
	    sprintf(gretl_errmsg, _("generation of lag variable failed"));
	    ok = 0;
	} else {
	    /* Record info regarding the auto-generation of lags
	       so that we'll be able to echo the command properly --
	       see the echo_cmd() function.  Note: 'lagvar.v' is the
	       "source" variable; 'lv' is the ID number of the
	       generated lag.
	    */
	    cmd->list[llen++] = lv;
	    cmd->err = list_lag_info_add(lagvar.v, order, lv, llen - 1, cmd);
	    if (cmd->err) {
		ok = 0;
	    }
	}
    }

    if (ok) {
	*lpos = llen;
    }

 bailout:

    if (lagvar.laglist != NULL) {
	free(lagvar.laglist);
    }

    return ok;
} 

static void parse_laglist_spec (const char *s, int *order, char **lname,
				int *vnum, const double **Z,
				const DATAINFO *pdinfo)
{
    int len = strcspn(s, ",;");

    if (len < strlen(s)) {
	char ostr[VNAMELEN] = {0};
	char word[32] = {0};
	int v;

	sscanf(s, "%15[^ ,;]", ostr);
	if (isdigit(*ostr)) {
	    *order = atoi(ostr);
	} else {
	    v = varindex(pdinfo, ostr);
	    if (v < pdinfo->v) {
		*order = Z[v][0];
	    }
	}
	sscanf(s + len + 1, "%31[^ )]", word);
	v = varindex(pdinfo, word);
	if (v < pdinfo->v) {
	    *vnum = v;
	} else {
	    *lname = gretl_word_strdup(s + len + 1, NULL);
	}
    } else {
	*lname = gretl_word_strdup(s, NULL);
    }
}

static int auto_transform_ok (const char *s, int *lpos,
			      double ***pZ, DATAINFO *pdinfo,
			      CMD *cmd)
{
    char fword[9];
    int *genlist = NULL;
    int trans = 0;
    int order = 0;
    gretlopt opt = OPT_NONE;
    int ok = 1;

    if (sscanf(s, "%8[^(](", fword)) {
	char *param = NULL;
	int *gotlist;
	int vnum = 0;

	if (!strcmp(fword, "cross")) {
	    strcpy(fword, "square");
	    opt = OPT_O;
	} else if (!strcmp(fword, "log")) {
	    strcpy(fword, "logs");
	}

	trans = gretl_command_number(fword);
	if (!MODIFIES_LIST(trans)) {
	    trans = 0;
	}

	if (trans > 0) {
	    s = strchr(s, '(') + 1;

	    if (trans == LAGS) {
		parse_laglist_spec(s, &order, &param, &vnum,
				   (const double **) *pZ, pdinfo);
	    } else {
		param = gretl_word_strdup(s, NULL);
	    }

	    if (param != NULL) {
		/* try for a named list */
		gotlist = get_list_by_name(param);
		if (gotlist != NULL) {
		    genlist = gretl_list_copy(gotlist);
		} else {
		    vnum = varindex(pdinfo, param);
		    if (vnum == pdinfo->v) {
			vnum = 0;
		    }
		}
		free(param);
	    } 

	    if (genlist == NULL && vnum > 0) {
		/* try for a single variable */
		genlist = gretl_list_new(1);
		if (genlist != NULL) {
		    genlist[1] = vnum;
		}
	    }
	}
    }

    if (genlist == NULL) {
	cmd->err = E_PARSE;
	return 0;
    }	

    if (trans == LOGS) {
	cmd->err = list_loggenr(genlist, pZ, pdinfo);
    } else if (trans == DIFF || trans == LDIFF || trans == SDIFF) {
	cmd->err = list_diffgenr(genlist, trans, pZ, pdinfo);
    } else if (trans == ORTHDEV) {
	cmd->err = list_orthdev(genlist, pZ, pdinfo);
    } else if (trans == SQUARE) {
	cmd->err = list_xpxgenr(&genlist, pZ, pdinfo, opt);
    } else if (trans == LAGS) {
	cmd->err = list_laggenr(&genlist, order, pZ, pdinfo);
    } else if (trans == DUMMIFY) {
	cmd->err = list_dumgenr(&genlist, pZ, pdinfo, OPT_F);
    }

    if (!cmd->err) {
	cmd->list[0] -= 1;
	cmd->err = gretl_list_insert_list(&cmd->list, genlist, *lpos);
	if (!cmd->err) {
	    *lpos += genlist[0];
	}
    }

    if (cmd->err) {
	ok = 0;
    }

    free(genlist);

    return ok;
} 

static int add_time_ok (const char *s, int *lpos,
			double ***pZ, DATAINFO *pdinfo,
			CMD *cmd)
{
    int ok = 0;

    if (!strcmp(s, "time")) {
	if (cmd->ci == GNUPLOT) {
	    cmd->list[0] -= 1;
	    cmd->opt |= OPT_T;
	    ok = 1; /* handled */
	} else {
	    cmd->err = gen_time(pZ, pdinfo, 1);
	    if (!cmd->err) {
		cmd->list[*lpos] = varindex(pdinfo, "time");
		*lpos += 1;
		ok = 1; /* handled */
	    }
	}
    }

    return ok;
}

static int truncate_varname (const char *s, int *k, const DATAINFO *pdinfo,
			     CMD *cmd)
{
    int ok = 0;

    if (strlen(s) > 8) {
	char test[9];
	int v;

	*test = 0;
	strncat(test, s, 8);
	if ((v = varindex(pdinfo, test)) < pdinfo->v) {
	    cmd->list[*k] = v;
	    *k += 1;
	    ok = 1;
	} 
    } 

    return ok;
}

static int wildcard_expand (const char *s, int *lpos,
			    const DATAINFO *pdinfo, CMD *cmd)
{
    int ok = 0;

    if (strchr(s, '*') != NULL) {
	int *wildlist = varname_match_list(pdinfo, s);

	if (wildlist != NULL) {
	    int k, nw = wildlist[0];
	    int llen = *lpos;

	    if (expand_command_list(cmd, nw)) {
		return 0;
	    }
	    for (k=1; k<=nw; k++) {
		cmd->list[llen++] = wildlist[k];
	    }
	    free(wildlist);
	    *lpos = llen;
	    ok = 1;
	}
    }

    return ok;
}

static int matrix_name_ok (const char *s, CMD *cmd)
{
    int ok = 0;

    if (cmd->ci == PRINT) {
	if (get_matrix_by_name(s)) {
	    cmd->extra = gretl_str_expand(&cmd->extra, s, " ");
	    cmd->list[0] -= 1;
	    ok = 1;
	}
    }

    return ok;
}

static void parse_rename_cmd (const char *line, CMD *cmd, 
			      const DATAINFO *pdinfo)
{
    int vtest, vtarg;
    char targ[VNAMELEN];
    char newname[VNAMELEN];
    char numstr[8];

    line += strlen(cmd->word);

    if (sscanf(line, "%15s %15s", targ, newname) != 2) {
	cmd->err = E_DATA;
	sprintf(gretl_errmsg, "rename: %s", 
		_("requires a variable number and a new name"));
	return;
    }

    if (isdigit(*targ)) {
	vtarg = atoi(targ);
	if (vtarg >= pdinfo->v || vtarg < 1) {
	    cmd->err = E_DATA;
	    sprintf(gretl_errmsg, _("Variable number %d is out of bounds"), vtarg);
	    return;
	}
    } else {
	/* we're given the name of a variable? */
	vtarg = varindex(pdinfo, targ);
	if (vtarg >= pdinfo->v) {
	    cmd->err = E_UNKVAR;
	    return;
	}
    } 

    vtest = varindex(pdinfo, newname);
    if (vtest < pdinfo->v && vtest != vtarg) {
	sprintf(gretl_errmsg, _("'%s': there is already a variable "
				"of this name"), newname);
	cmd->err = E_DATA;
	return;
    }

    if (check_varname(newname)) {
	cmd->err = E_DATA;
	return;
    }

    free(cmd->param);
    cmd->param = gretl_strdup(newname);
    if (cmd->param == NULL) {
	cmd->err = E_ALLOC;
	return;
    }

    sprintf(numstr, "%d", vtarg);

    free(cmd->extra);
    cmd->extra = gretl_strdup(numstr);
}

static void parse_outfile_cmd (char *s, CMD *cmd)
{
    int len = 0, quoted = 0;

    s += strlen(cmd->word);

    while (isspace(*s)) {
	s++;
    }

    if (*s) {
	cmd->err = filename_to_param(cmd, s, &len, &quoted, NULL);
    }
}

static void parse_logistic_ymax (char *line, CMD *cmd)
{
    char *p = strstr(line, "ymax");

    if (p != NULL) {
	char *q = p + 4;
	char numstr[12];

	while (*q == ' ' || *q == '=') {
	    q++;
	}
	if (sscanf(q, "%11s", numstr)) {
	    cmd->param = realloc(cmd->param, 6 + strlen(numstr));
	    if (cmd->param == NULL) {
		cmd->err = E_ALLOC;
	    } else {
		sprintf(cmd->param, "ymax=%s", numstr);
	    }
	    *p = '\0';
	}
    }
}

static int read_dash_param (const char **ps, CMD *cmd)
{
    const char *s = *ps;
    char *test;
    int i, parm = -1;
    int ok = 0;

    if (!strncmp(s, "--sheet=", 8)) {
	s += 8;
	if (*s == '"' || *s == '\'') {
	    free(cmd->extra);
	    cmd->extra = gretl_quoted_string_strdup(s, (const char **) &test);
	    parm = 0;
	} else {
	    parm = strtol(s, &test, 10);
	}
	i = 1;
    } else if (!strncmp(s, "--coloffset=", 12)) {
	s += 12;
	parm = strtol(s, &test, 10);
	i = 2;
    } else if (!strncmp(s, "--rowoffset=", 12)) {
	s += 12;
	parm = strtol(s, &test, 10);
	i = 3;
    }

    if (parm >= 0 && (*test == '\0' || *test == ' ') && 
	!(cmd->list[0] == 3 && cmd->list[i] > 0)) { 
	ok = 1;
	*ps = test;
	if (cmd->list[0] < 3) {
	    free(cmd->list);
	    cmd->list = gretl_list_new(3);
	    if (cmd->list == NULL) {
		cmd->err = E_ALLOC;
	    }
	}
	if (!cmd->err) {
	    cmd->list[i] = parm;
	}
    }

    return ok;
}

static void parse_spreadsheet_params (const char *s, CMD *cmd)
{
    int i, ok = 1;

    /* skip command word */
    s += strcspn(s, " ");
    s += strspn(s, " ");

    if (*s == '"') {
	s = strchr(s+1, '"');
	if (s == NULL) {
	    cmd->err = E_PARSE;
	    return;
	}
	s++;
    } else {
	s += strcspn(s, " ");
    }

    s += strspn(s, " ");

    for (i=0; *s && i<3 && ok; i++) {
	ok = read_dash_param(&s, cmd);
	if (cmd->err) {
	    break;
	} else if (!ok) {
	    sprintf(gretl_errmsg, _("Parse error in '%s'"), s);
	    cmd->err = E_PARSE;
	} else {
	    s += strspn(s, " ");
	}
    }
}

#define FIELDLEN 64

static int get_field_length (const char *s)
{
    int inparen = 0;
    int len = 0;

    while (*s) {
	if (*s == '(') {
	    inparen++;
	} else if (*s == ')') {
	    inparen--;
	}
	if (!inparen && *s == ' ') {
	    break;
	}
	s++;
	len++;
    }

    if (len >= FIELDLEN) {
	fprintf(stderr, "list field in command is too long\n");
	len = -1;
    }

    return len;
}

static int get_next_field (char *field, const char *s)
{
    int len, err = 0;

    *field = '\0';

#if CMD_DEBUG
    fprintf(stderr, "get_next_field: input = '%s'\n", s);
#endif

    while (*s == ' ') s++;

    len = get_field_length(s);

    if (len >= 0) {
	strncat(field, s, len);
    } else {
	err = 1;
    }

#if CMD_DEBUG
    fprintf(stderr, "get_next_field: got '%s'\n", field);
#endif

    return err;
}

static void accommodate_obsolete_commands (char *line, CMD *cmd)
{
    if (!strcmp(cmd->word, "noecho")) {
	strcpy(cmd->word, "set");
	strcpy(line, "set echo off");
    } else if (!strcmp(cmd->word, "seed")) {
	char seedstr[16];

	strcpy(cmd->word, "set");
	if (sscanf(line, "%*s %15s", seedstr)) {
	    sprintf(line, "set seed %s", seedstr);
	} else {
	    strcpy(line, "set seed");
	}
    } else if (!strcmp(cmd->word, "list") &&
	       string_is_blank(line + 4)) {
	strcpy(cmd->word, "varlist");
	strcpy(line, "varlist");
    }
}

/* look for a line with an "implicit genr", such as
   y = 3*x, x += 10, etc. */

int plausible_genr_start (const char *s, const DATAINFO *pdinfo)
{
    int ret = 0;

    if (strchr(s, '=') || strstr(s, "++") || strstr(s, "--")) {
	const char *ok = "+-*/=[";
	char word[VNAMELEN];

	if (sscanf(s, "%15[^[ +-*/=]", word)) {
	    s += strlen(word);
	    while (*s == ' ') s++;
	    if (strspn(s, ok) && check_varname(word) == 0) {
		ret = 1;
	    }
	}
    } else if (pdinfo != NULL && varindex(pdinfo, s) < pdinfo->v) {
	ret = 1;
    } else if (get_matrix_by_name(s)) {
	ret = 1;
    } else if (get_list_by_name(s)) {
	ret = 1;
    } else if (get_string_by_name(s)) {
	ret = 1;
    } 

    return ret;
}

/* if we find a semicolon without a preceding or following space,
   insert a space so that we can count the fields in the line
   correctly */

static int fix_semicolon_separation (char *s, CMD *cmd)
{
    int len = strlen(s);
    int i, j;

    for (i=0; i<len-1; i++) {
	if ((s[i] != ' ' && s[i+1] == ';') ||
	    (s[i] == ';' && s[i+1] && s[i+1] != ' ')) {
	    if (len < MAXLINE - 1) {
		for (j=len; j>i+1; j--) {
		    s[j] = s[j-1];
		}
		s[i+1] = ' ';
		s[len + 1] = '\0';
		len++;
	    } else {
		cmd->err = E_TOOLONG;
		break;
	    }
	} 
    }

    return len;
}

static int check_datamod_command (CMD *cmd, const char *s)
{
    cmd->aux = dataset_op_from_string(cmd->param);

    if (cmd->aux == DS_NONE) {
	cmd->err = E_PARSE;
    } else {
	/* skip param word and space */
	s += strcspn(s, " ");
	s += strspn(s, " ");
	free(cmd->param);
	cmd->param = gretl_strdup(s);
	if (cmd->param == NULL) {
	    cmd->err = E_ALLOC;
	} 
    }

#if CMD_DEBUG
    fprintf(stderr, "check_datamod_command: param='%s', aux = %d\n", 
	    cmd->param, cmd->aux);
#endif

    return cmd->err;
}

/* apparatus for checking that the "end" command is valid */

#define COMMAND_CAN_END(c) (c == FUNC || \
                            c == GMM || \
                            c == MLE || \
                            c == NLS || \
			    c == RESTRICT || \
			    c == SYSTEM)

static int check_end_command (CMD *cmd)
{
    if (cmd->param != NULL && *cmd->param != 0) {
	int cmdcode = gretl_command_number(cmd->param);

	if (cmdcode == LOOP) {
	    cmd->ci = ENDLOOP;
	} else if (!COMMAND_CAN_END(cmdcode)) {
	    cmd->err = 1;
	    sprintf(gretl_errmsg, _("command 'end %s' not recognized"), 
		    cmd->param);
	}
    } else {
	cmd->err = 1;
	strcpy(gretl_errmsg, _("end: nothing to end")); 
    }

    return cmd->err;
}

static void cmd_param_grab_string (CMD *cmd, const char *s)
{
    free(cmd->param);
    cmd->param = gretl_strdup(s);
    if (cmd->param == NULL) {
	cmd->err = E_ALLOC;
    }
}

static void cmd_param_grab_word (CMD *cmd, const char *s, int *nf)
{
    int n = strcspn(s, " =\n\t");

    if (n > 0) {
	free(cmd->param);
	cmd->param = gretl_strndup(s, n);
	if (cmd->param == NULL) {
	    cmd->err = E_ALLOC;
	} 
    }
}

/* Capture the next 'word' found following the initial command word
   (or the whole remainder of the line in some cases) as the parameter
   for cmd.  Flag an error if the command requires a parameter but
   none is found.
*/

static int capture_param (const char *s, CMD *cmd,
			  int *nf, const double **Z, 
			  const DATAINFO *pdinfo)
{
    /* if param has already been written by some special
       routine, don't overwrite it */
    if (*cmd->param != '\0') {
	if (cmd->ci == DATAMOD) {
	    check_datamod_command(cmd, s);
	}
	return cmd->err;
    }

    /* skip past leading word on line */
    s += strcspn(s, " ");
    s += strspn(s, " ");

    if (string_is_blank(s)) {
	if (REQUIRES_PARAM(cmd->ci) || REQUIRES_ORDER(cmd->ci)) {
	    cmd->err = E_PARSE;
	    sprintf(gretl_errmsg, _("%s: required parameter is missing"),
		    cmd->word);
	}
    } else {
	if (cmd->ci == PRINT || cmd->ci == FUNCERR || cmd->ci == DELEET) {
	    /* grab the whole remainder of line */
	    cmd_param_grab_string(cmd, s);
	} else {
	    /* grab one 'word' */
	    cmd_param_grab_word(cmd, s, nf);
	}
#if CMD_DEBUG
	fprintf(stderr, "capture_param: s='%s', param='%s'\n",
		s, cmd->param);
#endif
	if (REQUIRES_ORDER(cmd->ci)) {
	    cmd->order = gretl_int_from_string(cmd->param,
					       Z, pdinfo, 
					       &cmd->err);
	}
    }

    if (cmd->ci == DATAMOD) {
	check_datamod_command(cmd, s);
    } else if (cmd->ci == END) {
	check_end_command(cmd);
    }

#if CMD_DEBUG
    fprintf(stderr, "capture_param: returning %d\n", cmd->err);
#endif

    return cmd->err;
}

static int gretl_cmd_clear (CMD *cmd)
{
    cmd->ci = 0;
    cmd->err = 0;
    *cmd->word = '\0';
    cmd->opt = OPT_NONE;

    cmd_unset_nolist(cmd);

    if (cmd->list == NULL || cmd->param == NULL || cmd->extra == NULL) {
	cmd->err = E_ALLOC;
    } else {
	cmd->list[0] = 0;
	*cmd->param = '\0';
	*cmd->extra = '\0';
    }

    cmd_lag_info_destroy(cmd);

    return cmd->err;
}

static int resize_command_list (CMD *cmd, int nf)
{
    int *list;
    int i;

    if (nf < 0) {
	return 0;
    }

    list = realloc(cmd->list, (1 + nf) * sizeof *cmd->list);

    if (list == NULL) {
	cmd->err = E_ALLOC;
	strcpy (gretl_errmsg, _("Memory allocation failed for command list"));
    } else {
	list[0] = nf;
	for (i=1; i<=nf; i++) {
	    list[i] = 0;
	}
	cmd->list = list;
    }

    return cmd->err;
}

/* below: count fields, considering space as the field separator but
   only in case the material is not 'glued together' with parentheses
*/

int count_free_fields (const char *s)
{
    int inparen = 0;
    int nf = 0;

#if CMD_DEBUG
    fprintf(stderr, "count_free_fields: looking at '%s'\n", s);
#endif

    if (s != NULL && *s != '\0') {
	/* step past any leading spaces */
	while (*s == ' ') {
	    s++;
	}

	if (*s) {
	    s++;
	    nf++;
	}

	while (*s) {
	    if (*s == '(') {
		inparen++;
	    } else if (*s == ')') {
		inparen--;
	    } 
	    if (!inparen && *s == ' ') {
		while (*s == ' ') s++;
		if (*s) {
		    nf++;
		} else {
		    break;
		}
	    }
	    s++;
	}
    }

#if CMD_DEBUG
    fprintf(stderr, "count_free_fields: nf = %d\n", nf);
#endif
	    
    return nf;
}

static int get_sepcount (const char *s)
{
    int c = 0;

    while (*s++) {
	if (*s == ';') c++;
    }

    return c;
}

static char *copy_remainder (const char *line, int pos)
{
    char *rem;

    if (*(line + pos) == '\0') {
	rem = gretl_strdup(line + pos);
    } else {
	rem = gretl_strdup(line + pos + 1);
    }

    return rem;
}

static int process_cmd_extra (CMD *cmd, const double **Z,
			      const DATAINFO *pdinfo)
{
    if (cmd->ci == VECM) {
	cmd->aux = gretl_int_from_string(cmd->extra, Z, pdinfo, 
					 &cmd->err);
    }

    return cmd->err;
}

static int handle_semicolon (int *k, int *ints_ok, int *poly, 
			     int *sepcount, CMD *cmd)
{
    int ok = 0;

    if (USES_LISTSEP(cmd->ci)) {
	cmd->list[*k] = LISTSEP;
	*k += 1;
	*sepcount -= 1;
	if (*ints_ok) {
	    if (*sepcount == 0 || (cmd->ci == ARBOND && *sepcount == 1)) {
		*ints_ok = 0;
	    }
	}	
	if (cmd->ci == MPOLS) { 	 
	    *poly = 1; 	 
	}
	ok = 1;
    } 

    return ok;
}

static int get_id_or_int (const char *s, int *k, int ints_ok, int poly,
			  const DATAINFO *pdinfo, CMD *cmd)
{
    int v = atoi(s);
    int ok = 0;

    if (!ints_ok && !poly && v >= pdinfo->v) {
	cmd->err = 1;
	sprintf(gretl_errmsg, _("%d is not a valid variable number"), v);
    } else {
	cmd->list[*k] = v;
	*k += 1;
	ok = 1;
    }

    return ok;
}

static int parse_alpha_list_field (const char *s, int *pk, int ints_ok,
				   double ***pZ, DATAINFO *pdinfo, 
				   CMD *cmd)
{
    int *xlist;
    int v, k = *pk;
    int ok = 0;

#if CMD_DEBUG
    fprintf(stderr, "parse_alpha_list_field: s = '%s'\n", s);
#endif

    if (ints_ok) {
	v = gretl_int_from_string(s, (const double **) *pZ,
				  pdinfo, &cmd->err);
	if (!cmd->err) {
	    cmd->list[k++] = v;
	    ok = 1;
	}
    } else if ((v = varindex(pdinfo, s)) < pdinfo->v) {
	cmd->list[k++] = v;
	ok = 1;
    } else if ((xlist = get_list_by_name(s)) != NULL) {
	if (cmd->list[0] == 1 && xlist[0] == 0) {
	    cmd->list[0] = 0;
	    cmd->flags |= CMD_NULLIST;
	    ok = 1;
	} else {
	    cmd->list[0] -= 1;
	    cmd->err = gretl_list_insert_list(&cmd->list, xlist, k);
	    if (!cmd->err) { 
		k += xlist[0];
		ok = 1;
	    }
	}
    } else if (strchr(s, '(') != NULL) {
	if (auto_lag_ok(s, &k, pZ, pdinfo, cmd)) {
	    /* lag specification, e.g. 'var(-1)' */
	    ok = 1;
	} else if (auto_transform_ok(s, &k, pZ, pdinfo, cmd)) {
	    /* automated transformations such as 'logs(list)' */
	    ok = 1;	
	}
    } else if (add_time_ok(s, &k, pZ, pdinfo, cmd)) {
	ok = 1;	
    } else if (wildcard_expand(s, &k, pdinfo, cmd)) {
	ok = 1;
    } else if (truncate_varname(s, &k, pdinfo, cmd)) {
	ok = 1;
    } else if (cmd->ci == PRINT && matrix_name_ok(s, cmd)) {
	ok = 1;
    } 

    *pk = k;

    if (!ok && cmd->err == 0) {
	sprintf(gretl_errmsg, _("'%s' is not the name of a variable"), s);
	cmd->err = E_UNKVAR;
    }

    return ok;
}

static int sepcount_error (int ci, int nsep)
{
    int err = 0;

    if (NEEDS_LISTSEP(ci) && nsep == 0) {
	err = E_ARGS;
    } else if (!USES_LISTSEP(ci) && nsep > 0) {
	err = E_PARSE;
    } else if (!DOUBLE_SEP_OK(ci) && nsep == 2) {
	err = E_PARSE;
    } else if (nsep > 2) {
	err = E_PARSE;
    }

    return err;
}

/* Get the first word out of line.  In general this should be a
   command word (starting with a alphabetical character), but there
   are a few special case: shell commands start with the "!"  escape;
   restriction specifications may start with "-" (as in "-b1 + b2 =
   0") or a numerical multiplier.
*/

static int get_command_word (const char *line, CMD *cmd)
{
    int n = gretl_varchar_spn(line);
    int ret = 0;

    if (cmd->context == RESTRICT && n == 0) {
	ret = !string_is_blank(line);
    } else if (*line == '!') {
	strcpy(cmd->word, "!");
	ret = 1;
    } else if (n > 0) {
	if (n > FN_NAMELEN - 1) {
	    n = FN_NAMELEN - 1;
	}
	*cmd->word = '\0';
	strncat(cmd->word, line, n);
	ret = 1;
    } else if (!string_is_blank(line)) {
	cmd->err = E_PARSE;
    }

    return ret;
}

/**
 * parse_command_line:
 * @line: the command line.
 * @cmd: pointer to command struct.
 * @pZ: pointer to data matrix.
 * @pdinfo: pointer to data information struct.
 *
 * Parses @line and fills out @cmd accordingly. 
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int parse_command_line (char *line, CMD *cmd, double ***pZ, DATAINFO *pdinfo) 
{
    int j, k, nf, linelen, pos;
    int poly = 0;
    int sepcount = 0;
    int ints_ok = 0;
    int subst = 0;
    char *rem = NULL;
    char s[FIELDLEN] = {0};

    if (gretl_cmd_clear(cmd)) {
	return cmd->err;
    }

    gretl_error_clear();

    cmd->err = substitute_named_strings(line, &subst);
    if (cmd->err) {
	return cmd->err;
    }

    if (subst) {
	/* record the fact that substitution has been done */
	cmd->flags |= CMD_SUBST;
    } else {
	cmd->flags &= ~CMD_SUBST;
    }

    compress_spaces(line);

#if CMD_DEBUG
    fprintf(stderr, "parsing '%s'\n", line);
#endif

    /* trap lines that are nothing but comments */
    if (filter_comments(line, cmd)) {
	return 0;
    }

    /* catch errors associated with comment syntax */
    if (cmd->err) {
	return cmd->err;
    }

    /* extract any options */
    cmd->opt = get_gretl_options(line, &cmd->err);
    if (cmd->err) {
	return cmd->err;
    }    

    /* extract "savename" for storing an object? */
    maybe_extract_savename(line, cmd);

    /* no command here? */
    if (!get_command_word(line, cmd)) {
	cmd_set_nolist(cmd);
	cmd->ci = CMD_NULL;
	return cmd->err;
    }

#if CMD_DEBUG
    fprintf(stderr, "cmd->word = '%s'\n", cmd->word);
#endif

    /* backwards compatibility */
    accommodate_obsolete_commands(line, cmd);

    /* replace simple aliases and a few specials */
    catch_command_alias(line, cmd);

    /* subsetted commands (e.g. "deriv" in relation to "nls") */
    if (!strcmp(cmd->word, "end")) {
	cmd->context = 0;
    } else if (cmd->context && strcmp(cmd->word, "equation")) {
	/* "equation" occurs in the SYSTEM context, but it is
	   a command in its own right */
	cmd->ci = cmd->context;
    }

    if (cmd->ci == 0) {
	cmd->ci = gretl_command_number(cmd->word);
	if (cmd->ci == 0) {
	    if (plausible_genr_start(line, pdinfo)) {
		cmd->ci = GENR;
	    } else if (get_user_function_by_name(cmd->word)) {
		cmd->ci = GENR;
		cmd->opt = OPT_U;
	    } else if (is_gretl_function_call(line)) {
		cmd->ci = GENR;
		cmd->opt = OPT_U;
	    } else if (get_if_state(IS_FALSE)) {
		cmd_set_nolist(cmd);
		cmd->ci = CMD_NULL;
		return 0;
	    } else {
		cmd->err = 1;
		sprintf(gretl_errmsg, _("command '%s' not recognized"), 
			cmd->word);
		goto cmd_exit;
	    }
	}
    }

#if CMD_DEBUG
    fprintf(stderr, "cmd->ci = %d\n", cmd->ci);
#endif

    /* if, else, endif controls: should this come earlier? */
    if (flow_control(line, pZ, pdinfo, cmd)) {
	cmd_set_nolist(cmd);
	cmd->ci = CMD_NULL;
	return cmd->err;
    }

    if (cmd->ci == EQNPRINT || cmd->ci == TABPRINT) {
	/* TeX printing commands can take a filename parameter, and
	   possibly a format string -- but that's all
	*/
	get_optional_filename_etc(line, cmd);
	return cmd->err;
    } else if (cmd->ci == OUTFILE) {
	/* the "outfile" command may have a filename */
	parse_outfile_cmd(line, cmd);
    } else if (cmd->ci == RENAME) {
	/* the "rename" command calls for a variable number and a
	   new name */
	parse_rename_cmd(line, cmd, pdinfo);
    } else if (cmd->ci == OPEN || cmd->ci == APPEND) {
	/* "open" and "append" may have spreadsheet parameters */
	if (!(cmd->opt & OPT_O)) {
	    parse_spreadsheet_params(line, cmd);
	    if (cmd->err) {
		return cmd->err;
	    }
	}
    } 

    /* commands that never take a list of variables */
    if (NO_VARLIST(cmd->ci) || 
	(cmd->ci == DELEET && (cmd->opt & OPT_D))) { 
	cmd_set_nolist(cmd);
	capture_param(line, cmd, NULL, (const double **) *pZ, pdinfo);
	return cmd->err;
    }

    /** now for a few commands which may or may not take a list **/

    if (cmd->ci == PRINT && strstr(line, "\"")) {
	/* no list in string literal variant */
	cmd_set_nolist(cmd);
	capture_param(line, cmd, NULL, NULL, NULL);
	return cmd->err;
    }

    /* SMPL can take a list, but only in case of OPT_M
       "--no-missing" */
    if (cmd->ci == SMPL && !(cmd->opt & OPT_M)) {
	cmd_set_nolist(cmd);
	return cmd->err;
    }	

    /* boxplots take a list, but if there are Boolean conditions
       embedded, the line has to be parsed specially */
    if (cmd->ci == BXPLOT && strchr(line, '(')) {
	cmd_set_nolist(cmd);
	return cmd->err;
    }

    /* OMIT typically takes a list, but can be given without args
       to omit the last variable */
    if (cmd->ci == OMIT && string_is_blank(line + 4)) {
	cmd_set_nolist(cmd);
	return cmd->err;
    } 

    /* dataset-modifying commands */
    if (cmd->ci == DATAMOD) {
	capture_param(line, cmd, NULL, NULL, NULL);
	cmd_set_nolist(cmd);
	return cmd->err;
    }

    /** OK, now we're definitely doing a list-oriented command,
	We begin by taking care of a few specials **/

    if (cmd->ci == GNUPLOT) {
	/* we may have a block of stuff to pass literally
	   to gnuplot */
	grab_gnuplot_literal_block(line, cmd);
    } else if (cmd->ci == LOGISTIC) {
	/* we may have a "ymax" parameter */
	parse_logistic_ymax(line, cmd);
    } else if (cmd->ci == ARMA) {
	/* allow for specific "gappy" lags */
	arma_maybe_rewrite(line, cmd);
    }

    /* fix lines that contain a semicolon stuck to another element */
    linelen = fix_semicolon_separation(line, cmd);
    if (cmd->err) {
	return cmd->err;
    }

    /* arbond special: if there's a block-diagonal instruments
       portion to the command, grab that in literal form for
       later processing */
    if (cmd->ci == ARBOND && get_sepcount(line) == 2) {
	grab_arbond_diag(line, cmd);
	if (cmd->err) {
	    return cmd->err;
	}
    }

    /* find number of space-separated fields remaining in line,
       record our reading position, and make a copy of the
       remainder of the line
    */
    nf = count_free_fields(line) - 1;
    pos = strlen(cmd->word);
    rem = copy_remainder(line, pos);
    if (rem == NULL) {
	cmd->err = E_ALLOC;
	goto cmd_exit;
    }

#if CMD_DEBUG
    fprintf(stderr, "nf=%d, remainder='%s'\n", nf, rem);
#endif

    if (cmd->ci == DELEET) {
	if (nf == 1 && (get_matrix_by_name(rem) || get_string_by_name(rem))) {
	    /* special for deleting a named matrix or string */
	    cmd_param_grab_string(cmd, rem);
	    goto cmd_exit;
	}
    }

    /* specials where there's something that goes into "param",
       before the first semicolon */
    if (cmd->ci == RHODIFF || cmd->ci == LAGS) {
	if (get_rhodiff_or_lags_param(rem, cmd)) {
	    strcpy(line, rem);
	    linelen = strlen(line);
	    nf = count_fields(line);
	    pos = 0;
	} else if (cmd->ci == RHODIFF) {
	    /* rhodiff: param field is not optional */
	    cmd->err = E_ARGS;
	    goto cmd_exit;
	} else {
	    /* lags: param is optional */
	    *rem = '\0';
	}
    }    

    /* "store" is a special case since the filename that comes
       before the list may be quoted, and have spaces in it.  Groan */
    if (cmd->ci == STORE && nf > 0) {
	cmd->err = get_maybe_quoted_filename(cmd, rem, &nf);
	if (cmd->err) {
	    goto cmd_exit;
	} else {
	    pos = 0;
	    if (--nf > 0) {
		strcpy(line, rem);	
		linelen = strlen(line);
	    } 
	}
    }

    /* 
       "omitfrom" and "addto" take the ID of a previous model;
       "setmiss" takes a value to be interpreted as missing;
       these are captured in cmd->param
    */
    if (REQUIRES_ORDER(cmd->ci) ||
	cmd->ci == ADDTO ||
	cmd->ci == OMITFROM ||
	cmd->ci == SETMISS) {
	capture_param(line, cmd, &nf, (const double **) *pZ, pdinfo);
	if (cmd->err) {
	    goto cmd_exit;
	} else {
	    strcpy(rem, line + pos + 1 + strlen(cmd->param));
	    pos = 0;
	    if (--nf > 0) {
		strcpy(line, rem);
		linelen = strlen(line);
	    }
	} 
    }

    /* quantreg requires a tau specification */
    if (cmd->ci == QUANTREG) {
	capture_param(line, cmd, &nf, NULL, NULL);
	if (cmd->err) {
	    goto cmd_exit;
	} else {
	    strcpy(rem, line + pos + 1 + strlen(cmd->param));
	    pos = 0;
	    if (--nf > 0) {
		strcpy(line, rem);
		linelen = strlen(line);
	    }
	} 
    }    

    if (cmd->ci == OMITFROM && nf == 0) {
	cmd_set_nolist(cmd);
	return cmd->err;
    }

    if (cmd->ci == VECM) { 
	free(cmd->extra);
	cmd->extra = gretl_word_strdup(line, NULL);
	shift_string_left(line, strlen(cmd->extra));
	nf--;
	pos = 0;
	linelen = strlen(line);
	if (process_cmd_extra(cmd, (const double **) *pZ, pdinfo)) {
	    goto cmd_exit;
	}
    } 

    /* get a count of ';' separators in line */
    sepcount = get_sepcount(line);
    cmd->err = sepcount_error(cmd->ci, sepcount);
    if (cmd->err) {
	goto cmd_exit;
    }

    if (cmd->ci == AR || cmd->ci == ARBOND ||
	cmd->ci == ARMA || cmd->ci == GARCH) {
	/* flag acceptance of plain ints in list */
	ints_ok = 1;
    }

    /* allocate space for the command list */
    if (resize_command_list(cmd, nf)) {
	goto cmd_exit;
    }

    /* now assemble the command list */

    for (j=1, k=1; j<=nf; j++) {
	int ok = 0;

	strcpy(rem, line + pos + 1);

	/* special: optional width for correlogram, periodogram */
	if ((cmd->ci == CORRGM || cmd->ci == PERGM) && j == 2) {
	    cmd->list[0] = 1;
	    cmd_param_grab_word(cmd, rem, NULL);
	    break;
	} else if (cmd->ci == XCORRGM && j == 3) {
	    cmd->list[0] = 2;
	    cmd_param_grab_word(cmd, rem, NULL);
	    break;
	}	    

	cmd->err = get_next_field(s, rem);
	if (cmd->err) {
	    break;
	}

	if (isalpha((unsigned char) *s)) {
	    ok = parse_alpha_list_field(s, &k, ints_ok, pZ, pdinfo, cmd);
	} else if (*s == '*') {
	    ok = wildcard_expand(s, &k, pdinfo, cmd);
	} else if (isdigit(*s)) {
	    ok = get_id_or_int(s, &k, ints_ok, poly, pdinfo, cmd);
	} else if (*s == ';') {
	    ok = handle_semicolon(&k, &ints_ok, &poly, &sepcount, cmd);
	} 

	if (!ok) {
	    if (cmd->err == 0) {
		cmd->err = 1;
	    } 
	    if (*gretl_errmsg == '\0') {
		sprintf(gretl_errmsg, _("field '%s' in command is invalid"), s);
	    }
	    break;
	}

	/* check for bad scalars */
	if (cmd->list[k-1] != LISTSEP) {
	    if (!ints_ok && !poly && !(SCALARS_OK_IN_LIST(cmd->ci))) {
		if (var_is_scalar(pdinfo, cmd->list[k-1])) {
		    cmd->err = 1;
		    sprintf(gretl_errmsg, _("variable %s is a scalar"), s);
		    break;
		}
	    }
	}

	/* advance for next read */
	pos += strlen(s) + 1;
    }

    if (cmd->err) {
	goto cmd_exit;
    }

    /* By now we're looking at a command that takes a list,
       which either has been specified already or needs to
       be filled out automatically */

    /* commands that can take a specified list, but where if the
       list is null or just ";" we want to operate on all variables
    */    
    if (DEFAULTS_TO_FULL_LIST(cmd->ci) && !gretl_looping()) {
	if (cmd->list[0] == 0) {
	    cmd_full_list(pdinfo, cmd);
	    /* suppress echo of the list -- may be too long */
	    cmd_set_nolist(cmd);
	}
    } else if (cmd->ci != SETMISS && 
	       cmd->ci != PRINT &&
	       cmd->ci != DELEET) {
	/* the command needs a list but doesn't have one */
	if (cmd->list[0] == 0) {
	    cmd->err = E_ARGS;
	}
    }

    if (NEEDS_TWO_VARS(cmd->ci) && cmd->list[0] == 1) {
	cmd->err = E_ARGS;
    }

    if (cmd->ci == GNUPLOT && !(cmd->opt & OPT_T) && cmd->list[0] < 2) {
	cmd->err = E_ARGS;
    }

    /* check list for duplicated variables? */
    if (!cmd->err && !cmd_nolist(cmd)) {
	int dupv = gretl_list_duplicates(cmd->list, cmd->ci);

	if (dupv >= 0) {
	    printlist(cmd->list, "cmd->list with duplicate(s)");
	    cmd->err = E_UNSPEC;
	    sprintf(gretl_errmsg, 
		    _("var number %d duplicated in the command list."),
		    dupv);
	} 
    }

 cmd_exit:

    /* double-check that allocation hasn't failed */
    if (cmd->err == 0 && (cmd->list == NULL || cmd->param == NULL || 
			  cmd->extra == NULL)) {
	cmd->err = E_ALLOC;
    }

#if CMD_DEBUG
    printlist(cmd->list, "cmd->list");
    fprintf(stderr, "cmd->err = %d, context=%d\n", cmd->err,
	    cmd->context);
#endif

    free(rem);

    if (cmd->err) {
	cmd->context = 0;
    }

    return cmd->err;
}

static int maybe_need_recode (void)
{
    const gchar *cset = NULL;
    int utf = g_get_charset(&cset);

    return !utf;
}

static int recode_print_line (const char *s, PRN *prn)
{
    gchar *trs;
    gsize bytes;
    GError *err = NULL;

    trs = g_locale_from_utf8(s, -1, NULL, &bytes, &err);  

    if (err != NULL) {
	pprintf(prn, "%s\n", err->message);
	g_error_free(err);
    } else {
	pputs(prn, trs);
    }

    if (trs != NULL) {
	g_free(trs);
    }

    return 0;
}

/* list the topics available in the functions help file */

static int func_help_topics (const char *helpfile, PRN *prn)
{
    char line[128], word[12];
    FILE *fp;
    int j, k;

    if ((fp = gretl_fopen(helpfile, "r")) == NULL) {
	printf(_("Unable to access the file %s.\n"), helpfile);
	return E_FOPEN;
    } 

    j = 1;
    k = 0;
    while (fgets(line, sizeof line, fp) != NULL) {
	if (!strncmp(line, "## ", 3)) {
	    /* sub-heading */
	    tailstrip(line);
	    if (k++ > 0) {
		pputc(prn, '\n');
	    }
	    pprintf(prn, "\n%s:\n", line + 3);
	    j = 1;
	} else if (*line == '#') {
	    /* actual entry */
	    sscanf(line + 2, "%10s", word);
	    pprintf(prn, "%-10s", word);
	    if (j % 7 == 0) {
		pputc(prn, '\n');
	    } else {
		pputc(prn, ' ');
	    }
	    j++;
	}
    } 

    pputs(prn, _("\n\nFor help on a specific function, type: help funname"));
    pputs(prn, _(" (e.g. help qrdecomp)\n"));

    fclose(fp);
    
    return 0;
}

/**
 * cli_help:
 * @cmdword: the command on which help is wanted.
 * @paths: pointer to gretl paths structure.
 * @opt: may include %OPT_F to give priority to functions
 * rather than commands.
 * @prn: pointer to gretl printing struct.
 *
 * Searches in the gretl helpfile for help on @cmdword and, 
 * if help is found, prints it to @prn.  If @cmdword is %NULL, 
 * lists the valid commands.
 *
 * Returns: 0 on success, 1 if the helpfile was not found or the
 * requested topic was not found.
 */

int cli_help (const char *cmdword, PATHS *paths, gretlopt opt, PRN *prn)
{
    static int recode = -1;
    char helpfile[FILENAME_MAX];
    FILE *fp;
    int noword, funhelp = (opt & OPT_F);
    char word[9];
    char line[128];
    int i, j, ok = 0;

    noword = (cmdword == NULL || *cmdword == '\0');

    if (noword && !funhelp) {
	pputs(prn, _("\nValid gretl commands are:\n"));
	j = 1;
	for (i=1; i<NC; i++) {
	    if (HIDDEN_COMMAND(i)) {
		continue;
	    }
	    pprintf(prn, "%-9s", gretl_command_word(i));
	    if (j % 8 == 0) {
		pputc(prn, '\n');
	    } else {
		pputc(prn, ' ');
	    }
	    j++;
	}

	pputs(prn, _("\n\nFor help on a specific command, type: help cmdname"));
	pputs(prn, _(" (e.g. help smpl)\n"));
	pputs(prn, _("You can also do 'help functions' for a list of functions\n"));

	return 0;
    }

    if (!strcmp(cmdword, "functions") || (noword && funhelp)) {
	sprintf(helpfile, "%s%s", paths->gretldir, _("genrcli.hlp"));
	return func_help_topics(helpfile, prn);
    }

    if (!funhelp) {
	ok = gretl_command_number(cmdword) > 0;
    } 

    if (ok) {
	strcpy(helpfile, paths->cli_helpfile);
    } else {
	if (genr_function_word(cmdword)) {
	    sprintf(helpfile, "%sgenrcli.hlp", paths->gretldir);
	} else if (gretl_is_public_user_function(cmdword)) {
	    return user_function_help(cmdword, prn);
	} else {
	    pprintf(prn, _("\"%s\" is not a gretl command.\n"), cmdword);
	    return 1;
	}
    }

    if ((fp = gretl_fopen(helpfile, "r")) == NULL) {
	printf(_("Unable to access the file %s.\n"), helpfile);
	return 1;
    } 

    if (!gretl_in_gui_mode() && recode < 0) {
	recode = maybe_need_recode();
    }

    ok = 0;
    while (fgets(line, sizeof line, fp) != NULL) {
	if (*line != '#') {
	    continue;
	}
	sscanf(line + 2, "%8s", word);
	if (!strcmp(cmdword, word)) {
	    ok = 1;
	    pprintf(prn, "\n%s\n", word);
	    while (fgets(line, sizeof line, fp)) {
		if (*line == '#') {
		    break;
		}
		if (recode > 0) {
		    recode_print_line(line, prn);
		} else {
		    pputs(prn, line);
		}
	    }
	    break;
	}
    }

    if (!ok) {
	pprintf(prn, _("%s: sorry, no help available.\n"), cmdword);
    }

    fclose(fp);

    return 0;
}

static int parse_criteria (const char *line, const double **Z,
			   const DATAINFO *pdinfo, PRN *prn)
{
    double ess;
    int i, T, k;
    char cmd[9], essstr[32], Tstr[9], kstr[9];
    
    if (sscanf(line, "%s %s %s %s", cmd, essstr, Tstr, kstr) != 4) {
	return 1;
    }

    if (isalpha((unsigned char) *essstr) && 
	(i = varindex(pdinfo, essstr)) < pdinfo->v) {
	ess = get_xvalue(i, Z, pdinfo);
    } else if (isdigit(*essstr)) {
	ess = atof(essstr);
    } else {
	return 1;
    }

    if (ess < 0) {
	pputs(prn, _("ess: negative value is out of bounds.\n"));
	return 1;
    }

    if (isalpha((unsigned char) *Tstr) &&
	(i = varindex(pdinfo, Tstr)) < pdinfo->v) { 
	T = (int) get_xvalue(i, Z, pdinfo);
    } else if (isdigit(*Tstr)) {
	T = atoi(Tstr);
    } else {
	return 1;
    }

    if (T < 0) {
	pputs(prn, _("T: negative value is out of bounds.\n"));
	return 1;
    }

    if (isalpha((unsigned char) *kstr) &&
	(i = varindex(pdinfo, kstr)) < pdinfo->v) {
	k = (int) get_xvalue(i, Z, pdinfo);
    } else if (isdigit(*kstr)) {
	k = atoi(kstr);
    } else {
	return 1;
    }

    if (k < 0) {
	pputs(prn, _("k: negative value is out of bounds.\n"));
	return 1;
    }   
 
    gretl_print_criteria(ess, T, k, prn);

    return 0;
}

/**
 * parseopt:
 * @argv: command-line argument array.
 * @argc: argument count.
 * @fname: optional filename argument.
 * @force_lang: pointer to store result of "force language" option.
 *
 * Returns: the gretl option code corresponding to the first "real"
 * option flag, or 0 if the option flag is not recognized.
 */

int parseopt (const char **argv, int argc, char *fname, int *force_lang)
{
    int opt = 0, extra = 0;

    *fname = '\0';
    *force_lang = 0;

    while (*++argv) {
	const char *s = *argv;

#ifdef ENABLE_NLS
	if (!strcmp(s, "-e") || !strncmp(s, "--english", 9)) { 
	    *force_lang = ENGLISH;
	    continue;
	} else if (!strcmp(s, "-q") || !strncmp(s, "--basque", 8)) { 
	    *force_lang = BASQUE;
	    continue;
	}
#endif

	if (!strcmp(s, "-b") || !strncmp(s, "--batch", 7)) { 
	    opt = OPT_BATCH;
	} else if (!strcmp(s, "-h") || !strcmp(s, "--help")) { 
	    opt = OPT_HELP;
	} else if (!strcmp(s, "-v") || !strcmp(s, "--version")) { 
	    opt = OPT_VERSION;
	} else if (!strcmp(s, "-r") || !strncmp(s, "--run", 5)) { 
	    opt = OPT_RUNIT;
	} else if (!strcmp(s, "-d") || !strncmp(s, "--db", 4)) { 
	    opt = OPT_DBOPEN;
	} else if (!strcmp(s, "-w") || !strncmp(s, "--webdb", 7)) { 
	    opt = OPT_WEBDB;
	} else if (!strcmp(s, "-c") || !strncmp(s, "--dump", 6)) {
	    opt = OPT_DUMP;
	} else if (!strcmp(s, "-g") || !strncmp(s, "--debug", 7)) {
	    extra = OPT_DEBUG;
	} else {
	    strncat(fname, s, MAXLEN - 1);
	}
    }

    opt |= extra;

    return opt;
}

#ifndef WIN32

static int gretl_shell_async (const char *arg, PRN *prn)
{
    GError *gerr = NULL;
    int err = 0;

    g_spawn_command_line_async(arg, &gerr);

    if (gerr != NULL) {
	pprintf(prn, "%s\n", gerr->message);
	g_error_free(gerr);
	err = 1;
    }    

    return err;
}

static int gretl_shell_sync (const char *arg, gchar **psout,
			     PRN *prn)
{
    gchar *sout = NULL;
    gchar *serr = NULL;
    GError *gerr = NULL;
    int status;
    gchar *argv[5];
    const char *theshell = getenv("SHELL");
    const char *namep;
    char shellnam[40];
    int err = 0;

    if (theshell == NULL) {
#ifdef HAVE_PATHS_H
	theshell =_PATH_BSHELL;
#else
	theshell = "/bin/sh"; 
#endif
    }

    namep = strrchr(theshell, '/');
    if (namep == NULL) {
	namep = theshell;
    }

    strcpy(shellnam, "-");
    strcat(shellnam, ++namep);
    if (strcmp(namep, "sh") != 0) {
	shellnam[0] = '+';
    }

    argv[0] = g_strdup(theshell);
    argv[1] = shellnam;
    argv[2] = g_strdup("-c");
    argv[3] = g_strdup(arg);
    argv[4] = NULL;

    g_spawn_sync(get_shelldir(), argv, NULL, 0, NULL, NULL,
		 &sout, &serr, &status, &gerr); 

    g_free(argv[0]);
    g_free(argv[2]);
    g_free(argv[3]);

    if (gerr != NULL) {
	if (prn != NULL) {
	    pprintf(prn, "%s\n", gerr->message);
	} else {
	    gretl_errmsg_set(gerr->message);
	}
	g_error_free(gerr);
	err = 1;
    }

    if (psout != NULL) {
	*psout = sout;
    } else if (sout != NULL) {
	pputs(prn, sout);
	g_free(sout);
    }

    if (serr != NULL) {
	pputs(prn, serr);
	g_free(serr);
    }

    return err;
}

/**
 * gretl_shell_grab:
 * @arg: command line to be executed.
 * @sout: location to receive output from command.
 *
 * Calls the shell to execute @arg syncronously and captures the
 * standard output, if any, in @sout.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_shell_grab (const char *arg, char **sout)
{
    return gretl_shell_sync(arg, sout, NULL);
}

static int gretl_shell (const char *arg, PRN *prn)
{
    int async = 0;
    int err = 0;
    
    if (arg == NULL || *arg == '\0') {
	return 0;
    }

    if (!libset_get_bool(SHELL_OK)) {
	strcpy(gretl_errmsg, _("The shell command is not activated."));
	return 1;
    }

    if (!strncmp(arg, "launch ", 7)) {
	async = 1;
	arg += 7;
    } else if (*arg == '!') {
	arg++;
    }

    arg += strspn(arg, " \t");

    if (async) {
	err = gretl_shell_async(arg, prn);
    } else {
	err = gretl_shell_sync(arg, NULL, prn);
    }

    return err;
}

#endif /* ! WIN32 */

#define SAFELEN 78

static void trim_to_length (char *s)
{
    int i, n = strlen(s);

    if (n < SAFELEN - 1) return;

    for (i=n-1; i>0; i--) {
	if (s[i] == ' ') {
	    s[i] = '\0';
	    break;
	}
    }
}

void safe_print_line (const char *line, int *plen, PRN *prn)
{
    char tmp[SAFELEN];
    const char *q, *p = line;
    int n, m, rem, out = 0;
    int len0 = *plen;

    rem = n = strlen(line);

    while (out < n) {
	*tmp = 0;
	q = p;
	strncat(tmp, p, SAFELEN - 1);
	len0 = 0;
	trim_to_length(tmp - len0);
	len0 = 0;
	m = strlen(tmp);
	out += m;
	rem = n - out;
	p = q + m;
	if (rem > 0) {
	    pprintf(prn, "%s \\\n ", tmp);
	    *plen = 1;
	} else {
	    pprintf(prn, "%s", tmp);
	    *plen += m;
	}
    }
}

static int
print_maybe_quoted_str (const char *s, PRN *prn)
{
    int ret = 0;

    if (strchr(s, ' ') != NULL) {
	ret += pprintf(prn, " \"%s\"", s);
    } else {
	ret += pprintf(prn, " %s", s);
    }

    return ret;
}

static int is_real_model_cmd (const CMD *cmd)
{
    if (is_model_cmd(cmd->word) &&
	strcmp(cmd->word, "omit") &&
	strcmp(cmd->word, "add")) {
	return 1;
    } else {
	return 0;
    }
}

static int 
cmd_list_print_var (const CMD *cmd, int i, const DATAINFO *pdinfo,
		    int gotsep, PRN *prn)
{
    int src, v = cmd->list[i];
    int imin = 0;
    int bytes = 0;

    if (is_real_model_cmd(cmd)) {
	imin = 1;
    }

    if (v > 0 && i > imin && is_auto_generated_lag(i, cmd->list, cmd->linfo)) {
	if (is_first_lag(i, cmd->list, gotsep, cmd->linfo, &src)) {
	    bytes += print_lags_by_varnum(src, cmd->linfo, pdinfo, 
					  gotsep, prn);
	} 
    } else {
	pputc(prn, ' ');
	bytes += 1 + pputs(prn, pdinfo->varname[v]);
    }

    return bytes;
}

static int more_coming (const CMD *cmd, int i, int gotsep)
{
    if (cmd->opt) {
	return 1;
    } else if (cmd->linfo == NULL) {
	return (i < cmd->list[0]);
    } else {
	int j;

	for (j=i+1; j<=cmd->list[0]; j++) {
	    if (!is_auto_generated_lag(i, cmd->list, cmd->linfo) ||
		is_first_lag(i, cmd->list, gotsep, cmd->linfo, NULL)) {
		return 1;
	    }
	}
    }

    return 0;
}

static int n_separators (const int *list)
{
    int i, nsep = 0;

    for (i=2; i<list[0]; i++) {
	if (list[i] == LISTSEP) {
	    nsep++;
	}
    }

    return nsep;
}

static int effective_ci (const CMD *cmd)
{
    int ci = cmd->ci;

    if (ci == END) {
	if (!strcmp(cmd->param, "nls")) {
	    ci = NLS;
	} else if (!strcmp(cmd->param, "mle")) {
	    ci = MLE;
	} else if (!strcmp(cmd->param, "gmm")) {
	    ci = GMM;
	} else if (!strcmp(cmd->param, "restrict")) {
	    ci = RESTRICT;
	}
    }

    return ci;
}

#define listsep_switch(c) (c == AR || c == MPOLS)

#define hold_param(c) (c == TSLS || c == AR || c == ARBOND || c == ARMA || \
                       c == CORRGM || c == PERGM || c == SCATTERS || c == MPOLS || \
                       c == GNUPLOT || c == LOGISTIC || c == GARCH || \
                       c == EQUATION || c == POISSON || c == XCORRGM || \
                       c == HECKIT)

#define TESTLEN 62
#define LINELEN 78

static void
cmd_print_list (const CMD *cmd, const DATAINFO *pdinfo,  
		int *plen, PRN *prn)
{
    char numstr[12];
    int use_varnames = (cmd->ci != AR && cmd->ci != DELEET);
    int nsep, gotsep, i;

    if (cmd->list == NULL || cmd->list[0] == 0) {
	return;
    }

    nsep = n_separators(cmd->list);

    if (cmd->ci == RHODIFF) {
	*plen += pprintf(prn, " %s;", cmd->param);
    } else if (cmd->ci == LAGS) {
	if (cmd->param[0] != '\0') {
	    *plen += pprintf(prn, " %s;", cmd->param);
	}
    } else if (cmd->param[0] != '\0' && !hold_param(cmd->ci)) {
	*plen += print_maybe_quoted_str(cmd->param, prn);
    }

    if (cmd->ci == VECM && cmd->extra != NULL) {
	*plen += pprintf(prn, " %s", cmd->extra);
    }

    gotsep = 0;

    for (i=1; i<=cmd->list[0]; i++) {

	if (cmd->list[i] == LISTSEP) {
	    *plen += pputs(prn, " ;");
	    gotsep++;
	    if (listsep_switch(cmd->ci) && gotsep == nsep) {
		use_varnames = !use_varnames;
	    }
	    continue;
	}

	if (use_varnames) {
	    *plen += cmd_list_print_var(cmd, i, pdinfo, gotsep, prn);
	} else {
	    sprintf(numstr, " %d", cmd->list[i]);
	    *plen += pputs(prn, numstr);
	}

	if (*plen > TESTLEN && more_coming(cmd, i, gotsep)) {
	    pputs(prn, " \\\n "); 
	    *plen = 1;
	}
    }
}

#define ECHO_DEBUG 0

static int command_is_silent (const CMD *cmd, const char *line)
{
    if (cmd->ci == FUNCERR || cmd->ci == PRINTF ||
	(cmd->ci == PRINT && strchr(line, '"'))) {
	return 1;
    }

    if (!strncmp(line, "quit", 4) && string_is_blank(line + 4)) {
	return 1;
    }

    if (cmd->ci == SET && !strcmp(cmd->param, "echo") &&
	gretl_function_depth() > 0) {
	return 1;
    }

    if (*line == '!') {
	return 1;
    }

    return 0;
}

#define rewritten_arma(c) (c->ci == ARMA && \
                           c->extra != NULL && \
			   *c->extra != '\0')

/* these commands have sub-lists that may contain either
   numerical values or the names of scalar variables:
   this can't be handled properly by the list-printing
   mechanism */

#define dont_print_list(c) ((c->flags & CMD_NOLIST) || \
                             c->ci == ARBOND || \
                             c->ci == ARMA || \
			     c->ci == GARCH || \
                             c->ci == OPEN)

#define print_param_last(c) (c == ARBOND || \
			     c == DELEET || \
			     c == LOGISTIC)

/**
 * echo_cmd:
 * @cmd: pointer to #CMD struct.
 * @pdinfo: pointer to data information struct.
 * @line: "raw" command line associated with @cmd.
 * @flags: bitwise OR of elements from #CmdEchoFlags.
 * @prn: pointer to gretl printing struct (or %NULL).
 *
 * Echoes the user command represented by @cmd and @line, to
 * %stdout (if @prn is %NULL) or @prn.  This is used for two
 * distinct purposes: to give visual feedback on the
 * command supplied, and (in some contexts) to record a
 * command that was executed interactively.
 */

void echo_cmd (const CMD *cmd, const DATAINFO *pdinfo, const char *line, 
	       unsigned char flags, PRN *prn)
{
    int batch = (flags & CMD_BATCH_MODE);
    int recording = (flags & CMD_RECORDING);
    int len, llen = 0;

    if (line == NULL || prn == NULL || cmd->ci > NC) {
	return;
    }

#if ECHO_DEBUG
    fprintf(stderr, "echo_cmd:\n line='%s'\n param='%s'\n extra='%s'\n", 
	    line, cmd->param, cmd->extra);
    fprintf(stderr, " cmd->opt=%ld, batch=%d, nolist=%d\n",
	    cmd->opt, flags & CMD_BATCH_MODE, cmd_nolist(cmd));
    fprintf(stderr, " prn=%p\n", (void *) prn);
    fprintf(stderr, " cmd->word='%s'\n", cmd->word);
    if (!cmd_nolist(cmd)) {
	printlist(cmd->list, "cmd->list");
    }
#endif

    /* certain things don't get echoed at all, if not recording */
    if (!recording && command_is_silent(cmd, line)) {
	return;
    }

    /* special case: "store" command: record as comment */
    if (recording && cmd->ci == STORE) {
	pprintf(prn, "# store '%s'", cmd->param);
	if (cmd->opt) { 
	    const char *flagstr = print_flags(cmd->opt, cmd->ci);

	    pputs(prn, flagstr);
	}
	pputc(prn, '\n');
	return;
    }

    /* print leading string before echo: none if recording */
    if (!recording) {
	if ((flags & CMD_STACKING) || gretl_compiling_function()) {
	    llen += pputs(prn, "> ");
	} else if (batch) {
	    llen += pputs(prn, "? ");
	} else if (flags & CMD_CLI) {
	    llen += pputc(prn, ' ');
	}
    }

    /* special: printing a list */
    if (cmd->ci == PRINT && !strcmp(cmd->word, "list")) {
	pprintf(prn, "list %s print\n", cmd->extra);
	return;
    }

    if (*line == '\0') {
	return;
    }

    /* command is preceded by a "savename" to which an object will
       be assigned */
    if (*cmd->savename) {
	if (haschar(' ', cmd->savename) >= 0) {
	    pprintf(prn, "\"%s\" <- ", cmd->savename);
	    llen += strlen(cmd->savename) + 6;
	} else {
	    pprintf(prn, "%s <- ", cmd->savename);
	    llen += strlen(cmd->savename) + 4;
	}
    }

    if (dont_print_list(cmd)) {
	const char *s = line;
	
	if (rewritten_arma(cmd)) {
	    s = cmd->extra;
	}
	if (strlen(s) > SAFELEN - 2) {
	    safe_print_line(s, &llen, prn);
	} else {
	    llen += pputs(prn, s);
	}
    } else {
	llen += pprintf(prn, "%s", cmd->word);
	cmd_print_list(cmd, pdinfo, &llen, prn);
    } 

    /* print parameter after list, if wanted */
    if (print_param_last(cmd->ci) && *cmd->param != '\0') {
	len = strlen(cmd->param) + 1;
	if (llen + len > LINELEN) {
	    pputs(prn, " \\\n ");
	    llen = 0;
	}	    
	pputc(prn, ' ');
	pputs(prn, cmd->param);
	llen += len;
    }

    /* add printout of any options to the command */
    if (cmd->opt) {
	const char *flagstr;

	flagstr = print_flags(cmd->opt, effective_ci(cmd));
	len = strlen(flagstr);
	if (llen + len > LINELEN) {
	    pputs(prn, " \\\n ");
	}	    
	pputs(prn, flagstr);
    }

    pputc(prn, '\n');
    gretl_print_flush_stream(prn);
}

void echo_function_call (const char *line, unsigned char flags, PRN *prn)
{
    char leadchar = (flags & CMD_STACKING)? '>' : '?';

    if (gretl_echo_on()) {
	pprintf(prn, "%c %s\n", leadchar, line);
    }
}

/* Look for a flag of the form " -x" which occurs outside of any
   quotes: if found, return a pointer to the flag.
*/

static const char *flag_present (const char *s, char f, int *quoted)
{
    int inquote = 0;

#if CMD_DEBUG
    fprintf(stderr, "flag_present: looking at '%s'\n", s);
#endif

    while (*s) {
	if (*s == '"') {
	    inquote = !inquote;
	}
	if (!inquote) {
	    if (*s == ' ' && strlen(s) >= 4 && *(s+1) == '-' && *(s+2) == f) {
		s += 3;
		while (*s) {
		    if (isspace(*s)) s++;
		    else break;
		}
		if (*s == '"' && *(s+1)) {
		    *quoted = 1;
		    return s + 1;
		}
		if (*s != '"' && *(s+1)) {
		    *quoted = 0;
		    return s;
		}
	    } 
	}
	s++;
    }

#if CMD_DEBUG
    fprintf(stderr, "flag_present: returning '%s'\n", s);
#endif

    return NULL;
}

static char *get_flag_field  (const char *s, char f)
{
    char *ret = NULL;
    int quoted = 0;

    if ((s = flag_present(s, f, &quoted)) != NULL) {
	const char *p = s;
	int len = 0;

	while (*p) {
	    if (quoted && *p == '"') break;
	    if (!quoted && isspace(*p)) break;
	    p++;
	    len++;
	}

	ret = gretl_strndup(s, len);
    }

    return ret;
}

/* grab filename for TeX output; while we're at it, if the command is
   TABPRINT, get an optional format string too
*/

static void get_optional_filename_etc (const char *line, CMD *cmd)
{
    char *p = get_flag_field(line + 8, 'f');

    if (p != NULL && *p != '\0') {
	free(cmd->param);
	if (libset_get_bool(USE_CWD) || gretl_path_is_absolute(p)) {
	    cmd->param = p;
	} else {
	    cmd->param = gretl_strdup_printf("%s%s", gretl_work_dir(), p);
	    free(p);
	}
    }

    if (cmd->ci == TABPRINT) {
	const char *q = strstr(line, "--format=");
	int len;

	if (q != NULL && !strncmp(q + 9, "default", 7)) {
	    set_tex_param_format(NULL);
	} else if (q != NULL && *(q + 9) == '"') {
	    q += 10;
	    len = haschar('"', q);
	    if (len > 0) {
		p = gretl_strndup(q, len);
		if (p != NULL) {
		    set_tex_param_format(p);
		    free(p);
		}
	    }
	}
    }
}

static int set_var_info (const char *line, gretlopt opt, 
			 DATAINFO *pdinfo, PRN *prn)
{
    char *p, vname[VNAMELEN];
    int v, setstuff = 0;

    if (pdinfo->varinfo == NULL) {
	return 1;
    }

    /* skip command word plus space */
    line += strcspn(line, " ");
    line += strspn(line, " ");

    if (sscanf(line, "%15s", vname) != 1) {
	return E_PARSE;
    }

    /* skip varname, but not following space */
    line += strcspn(line, " ");

    v = varindex(pdinfo, vname);
    if (v == pdinfo->v) {
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), vname);
	return E_UNKVAR;
    }

    if (opt & OPT_D) {
	set_var_discrete(pdinfo, v, 1);
    } else if (opt & OPT_C) {
	set_var_discrete(pdinfo, v, 0);
    }

    p = get_flag_field(line, 'd');
    if (p != NULL) {
	setstuff = 1;
	*VARLABEL(pdinfo, v) = 0;
	strncat(VARLABEL(pdinfo, v), p, MAXLABEL - 1);
	free(p);
    }

    p = get_flag_field(line, 'n');
    if (p != NULL) {
	setstuff = 1;
	var_set_display_name(pdinfo, v, p);
	free(p);
    } 

    return 0;
}

static void showlabels (const DATAINFO *pdinfo, PRN *prn)
{
    const char *label;
    int i;

    pprintf(prn, _("Listing labels for variables:\n"));

    for (i=0; i<pdinfo->v; i++) {
	label = VARLABEL(pdinfo, i);
	if (strlen(label) > 2) {
	    pprintf(prn, " %s: %s\n", pdinfo->varname[i], label);
	}
    }
}

static void do_print_string (char *str, PRN *prn)
{
    size_t len;

    if (*str == '"') str++;

    len = strlen(str);

    if (str[len-1] == '"') {
	str[len-1] = 0;
    }

    pprintf(prn, "%s\n", str);
}

static int 
do_outfile_command (gretlopt flag, const char *fname, PRN *prn)
{
    static char outname[MAXLEN];
    int diverted = 0;
    int err = 0;

    if (prn == NULL) {
	return 0;
    }

    if (flag != OPT_W && flag != OPT_A && flag != OPT_C) {
	return E_ARGS;
    }

    diverted = printing_is_redirected(prn);

    /* command to close outfile */
    if (flag == OPT_C) {
	if (!diverted) {
	    pputs(prn, _("Output is not currently diverted to file\n"));
	    return 1;
	} else {
	    print_end_redirection(prn);
	    if (gretl_messages_on()) {
		pprintf(prn, _("Closed output file '%s'\n"), outname);
	    }
	    return 0;
	}
    }

    /* command to divert output to file */
    if (diverted) {
	fprintf(stderr, _("Output is already diverted to '%s'\n"),
		outname);
	return 1;
    } else if (*fname == '\0') {
	return E_ARGS;
    } else if (!strcmp(fname, "null")) {
	if (gretl_messages_on()) {
	   pputs(prn, _("Now discarding output\n")); 
	}
	print_start_redirection(prn, NULL);
	strcpy(outname, fname);
    } else {
	FILE *fp;

	fname = gretl_maybe_switch_dir(fname);

	if (flag == OPT_W) {
	    fp = gretl_fopen(fname, "w");
	} else {
	    fp = gretl_fopen(fname, "a");
	}

	if (fp == NULL) {
	    pprintf(prn, _("Couldn't open %s for writing\n"), fname);
	    return 1;
	}

	if (gretl_messages_on()) {
	    if (flag == OPT_W) {
		pprintf(prn, _("Now writing output to '%s'\n"), fname);
	    } else {
		pprintf(prn, _("Now appending output to '%s'\n"), fname);
	    }
	    
	}

	print_start_redirection(prn, fp);
	strcpy(outname, fname);
    }

    return err;
}

int call_pca_plugin (VMatrix *cmat, double ***pZ,
		     DATAINFO *pdinfo, gretlopt opt,
		     PRN *prn)
{
    void *handle = NULL;
    int (*pca_from_cmatrix) (VMatrix *, double ***, DATAINFO *,
			     gretlopt, PRN *);
    int err = 0;

    gretl_error_clear();
    
    pca_from_cmatrix = get_plugin_function("pca_from_cmatrix", &handle);
    if (pca_from_cmatrix == NULL) {
        return 1;
    }
        
    err = (* pca_from_cmatrix) (cmat, pZ, pdinfo, opt, prn);
    close_plugin(handle);
    
    return err;
}

static int do_pca (int *list, double ***pZ, DATAINFO *pdinfo,
		   gretlopt opt, PRN *prn)
{
    int err = 0;

    if (list[0] > 0) {
	VMatrix *cmat;

	/* adding OPT_U ensures a uniform sample for the
	   correlation or covariance matrix */
	cmat = corrlist(list, (const double **) *pZ, pdinfo, 
			opt | OPT_U, &err);
	if (!err) {
	    err = call_pca_plugin(cmat, pZ, pdinfo, opt, prn);
	    if (!err && (opt & (OPT_O | OPT_A))) {
		/* results saved as series */
		maybe_list_vars(pdinfo, prn);
	    }
	    free_vmatrix(cmat);
	}
    }

    return err;
}

static void print_info (gretlopt opt, DATAINFO *pdinfo, PRN *prn)
{
    if (pdinfo->descrip != NULL) {
	pprintf(prn, "%s\n", pdinfo->descrip);
    } else {
	pputs(prn, _("No data information is available.\n"));
    }
}

static int maybe_print_model (MODEL *pmod, DATAINFO *pdinfo,
			      PRN *prn, gretlopt opt)
{
    if (pmod->errcode == 0 && !gretl_looping()) {
	printmodel(pmod, pdinfo, opt, prn);
    }

    return pmod->errcode;
}

static int model_test_check (CMD *cmd, DATAINFO *pdinfo, PRN *prn)
{
    return last_model_test_ok(cmd->ci, cmd->opt, pdinfo, prn);
}

static int get_line_continuation (char *line, FILE *fp, PRN *prn)
{
    char tmp[MAXLINE];
    int err = 0;

    if (!strncmp(line, "quit", 4)) {
	return 0;
    }

    while (top_n_tail(line, &err)) {
	if (err) {
	    break;
	}
	*tmp = '\0';
	fgets(tmp, sizeof tmp, fp);
	if (*tmp != '\0') {
	    if (strlen(line) + strlen(tmp) > MAXLINE - 1) {
		pprintf(prn, _("Maximum length of command line "
			       "(%d bytes) exceeded\n"), MAXLINE);
		err = E_TOOLONG;
		break;
	    } else {
		strcat(line, tmp);
		compress_spaces(line);
	    }
	}
    }

    return err;
}

static int run_script (const char *fname, ExecState *s, 
		       double ***pZ, DATAINFO *pdinfo,
		       PRN *prn)
{
    FILE *fp;
    int indent = get_if_state(GETINDENT);
    int err = 0;
    
    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	sprintf(gretl_errmsg, _("Couldn't open %s"), fname);
	return E_FOPEN;
    }

    strcpy(s->runfile, fname);

    if (gretl_echo_on()) {
	pprintf(prn, "run \"%s\"\n", fname);
    }

    while (fgets(s->line, MAXLINE - 1, fp) && !err) {
	err = get_line_continuation(s->line, fp, prn);
	if (!err) {
	    err = maybe_exec_line(s, pZ, pdinfo);
	}
    }

    fclose(fp);

    if (get_if_state(GETINDENT) != indent) {
	set_if_state(RELAX);
	pputc(prn, '\n');
	pprintf(prn, _("Unmatched \"%s\""), "if");
	pputc(prn, '\n');
	if (!err) {
	    err = E_PARSE;
	}
    }

    return err;
}

static int append_data (const char *line, int *list, 
			char *sheetname,
			double ***pZ, DATAINFO *pdinfo, 
			gretlopt opt, PRN *prn)
{
    char fname[MAXLEN] = {0};
    int ftype, err = 0;

    err = getopenfile(line, fname, NULL, OPT_NONE);
    if (err) {
	errmsg(err, prn);
	return err;
    }

    ftype = detect_filetype(fname, NULL, prn);

    if (ftype == GRETL_CSV) {
	err = import_csv(fname, pZ, pdinfo, opt, prn);
    } else if (SPREADSHEET_IMPORT(ftype)) {
	err = import_spreadsheet(fname, ftype, list, sheetname, 
				 pZ, pdinfo, opt, prn);
    } else if (OTHER_IMPORT(ftype)) {
	err = import_other(fname, ftype, pZ, pdinfo, opt, prn);
    } else if (ftype == GRETL_XML_DATA) {
	err = gretl_read_gdt(fname, NULL, pZ, pdinfo, opt, prn);
    } else {
	err = gretl_get_data(fname, NULL, pZ, pdinfo, opt, prn);
    }

    return err;
}

static int do_end_restrict (ExecState *s, double ***pZ, DATAINFO *pdinfo)
{
    gretlopt ropt = gretl_restriction_get_options(s->rset);
    CMD *cmd = s->cmd;
    PRN *prn = s->prn;
    int err = 0;

    if ((cmd->opt & OPT_F) || (ropt & OPT_F)) {
	/* FIXME non-vecm case */
	s->var = gretl_restricted_vecm(s->rset, (const double **) *pZ, 
				       pdinfo, cmd->opt, prn, &err);
	if (s->var != NULL) {
	    if (s->callback != NULL) {
		s->callback(s, pZ, pdinfo);
	    }
	} 
    } else {
	err = gretl_restriction_finalize(s->rset, (const double **) *pZ, 
					 pdinfo, cmd->opt, prn);
    }

    s->rset = NULL;

    return err;
}

int gretl_cmd_exec (ExecState *s, double ***pZ, DATAINFO *pdinfo)
{
    CMD *cmd = s->cmd;
    char *line = s->line;
    MODEL **models = s->models;
    PRN *prn = s->prn;
    double rho;
    char runfile[MAXLEN];
    int *listcpy = NULL;
    int k, order = 0;
    int err = 0;

    if (NEEDS_MODEL_CHECK(cmd->ci)) {
	err = model_test_check(cmd, pdinfo, prn);
	if (err) {
	    return err;
	}
    }

    if (MODIFIES_LIST(cmd->ci)) {
	if (cmd->list[0] == 0) {
	    /* no-op */
	    return 0;
	} else {
	    /* list is potentially modified -> make a copy */
	    listcpy = gretl_list_copy(cmd->list);
	    if (listcpy == NULL) {
		return E_ALLOC;
	    }
	}
    }

    if (cmd->ci == OLS && dataset_is_panel(pdinfo)) {
	/* FIXME is this too funky? */
	cmd->ci = PANEL;
	cmd->opt |= OPT_P; /* panel pooled flag */
    }

    gretl_error_clear();

    s->alt_model = 0;

    switch (cmd->ci) {

    case APPEND:
	err = append_data(line, cmd->list, cmd->extra, pZ, pdinfo, 
			  cmd->opt, prn);
	break;

    case ADF:
	err = adf_test(cmd->order, cmd->list, pZ, pdinfo, cmd->opt, prn);
	break;

    case KPSS:
	err = kpss_test(cmd->order, cmd->list, pZ, pdinfo, cmd->opt, prn);
	break;

    case COINT:
	err = coint(cmd->order, cmd->list, pZ, pdinfo, cmd->opt, prn);
	break;

    case COINT2:
	err = johansen_test_simple(cmd->order, cmd->list, (const double **) *pZ, 
				   pdinfo, cmd->opt, prn);
	break;

    case CORR:
	err = incompatible_options(cmd->opt, OPT_U | OPT_S | OPT_K);
	if (err) {
	    break;
	}
	if (cmd->opt & OPT_K) {
	    err = kendall(cmd->list, (const double **) *pZ, pdinfo,
			  cmd->opt, prn);
	} else if (cmd->opt & OPT_S) {
	    err = spearman(cmd->list, (const double **) *pZ, pdinfo,
			   cmd->opt, prn);
	} else {
	    err = gretl_corrmx(cmd->list, (const double **) *pZ, pdinfo, 
			       cmd->opt, prn);
	}
	break;

    case CORRGM:
	order = atoi(cmd->param);
	err = corrgram(cmd->list[1], order, 0, (const double **) *pZ, pdinfo, 
		       prn, cmd->opt | OPT_A);
	break;

    case XCORRGM:
	order = atoi(cmd->param);
	err = xcorrgram(cmd->list, order, (const double **) *pZ, pdinfo, 
			prn, cmd->opt | OPT_A);
	break;

    case PERGM:
	order = atoi(cmd->param);
	err = periodogram(cmd->list[1], order, (const double **) *pZ, pdinfo, 
			  cmd->opt | OPT_N, prn);
	break;

    case BREAK:
    case ENDLOOP:
	pprintf(prn, _("You can't end a loop here, "
		       "you haven't started one\n"));
	err = 1;
	break;

    case FCAST:
	err = do_forecast(line, pZ, pdinfo, cmd->opt, prn);
	break;

    case FREQ:
	err = freqdist(cmd->list[1], (const double **) *pZ, 
		       pdinfo, (s->flags == CONSOLE_EXEC),
		       cmd->opt, prn);
	if (!err && !(cmd->opt & OPT_Q) && s->callback != NULL) {
	    s->callback(s, pZ, pdinfo);
	}
	break;

    case DISCRETE:
	err = list_makediscrete(cmd->list, pdinfo, cmd->opt);
	break;

    case ESTIMATE:
	err = estimate_named_system(line, pZ, pdinfo, cmd->opt, prn);
	break;

    case FUNC:
	err = gretl_start_compiling_function(line, prn);
	break;

    case GENR:
	err = generate(line, pZ, pdinfo, cmd->opt, prn);
	break;

    case PCA:
	err = do_pca(cmd->list, pZ, pdinfo, cmd->opt, prn);
	break;

    case CRITERIA:
	err = parse_criteria(line, (const double **) *pZ, pdinfo, prn);
	break;

    case DATA:
	err = db_get_series(line, pZ, pdinfo, cmd->opt, prn);
	break;

    case DATAMOD:
	err = modify_dataset(cmd->aux, cmd->list, cmd->param, pZ, 
			     pdinfo, prn);
	if (!err) { 
	    print_smpl(pdinfo, get_full_length_n(), prn);
	    if (s->callback != NULL) {
		s->callback(s, pZ, pdinfo);
	    }
	} 
	break;

    case DIFF:
    case LDIFF:
    case SDIFF:
	err = list_diffgenr(listcpy, cmd->ci, pZ, pdinfo);
	if (!err) {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case ORTHDEV:
	err = list_orthdev(listcpy, pZ, pdinfo);
	if (!err) {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case DUMMIFY:
	err = list_dumgenr(&listcpy, pZ, pdinfo, cmd->opt);
	if (!err) {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case LAGS:
	order = atoi(cmd->param);
	err = list_laggenr(&listcpy, order, pZ, pdinfo); 
	if (!err) {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case LOGS:
	err = list_loggenr(listcpy, pZ, pdinfo);
	if (!err) {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case SQUARE:
	err = list_xpxgenr(&listcpy, pZ, pdinfo, cmd->opt);
	if (!err) {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case GRAPH:
	err = ascii_graph(cmd->list, (const double **) *pZ, pdinfo, 
			  cmd->opt, prn);
	break;

    case PLOT:
	err = ascii_graph(cmd->list, (const double **) *pZ, pdinfo, 
			  (cmd->opt | OPT_T), prn);
	break;

    case RMPLOT:
    case HURST:
	if (cmd->list[0] != 1) {
	    pputs(prn, _("This command requires one variable.\n"));
	    err = 1;
	} else {
	    if (cmd->ci == RMPLOT) {
		err = rmplot(cmd->list, (const double **) *pZ, 
			     pdinfo, prn);
	    } else {
		err = hurstplot(cmd->list, (const double **) *pZ, 
				pdinfo, prn);
	    }
	}
	break;

    case INFO:
	print_info(cmd->opt, pdinfo, prn);
	break;

    case RENAME:
	err = rename_var_by_id(cmd->extra, cmd->param, pdinfo);
	if (!err) {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case SET:
	err = execute_set_line(line, *pZ, pdinfo, prn);
	break;

    case SETINFO:
	err = set_var_info(line, cmd->opt, pdinfo, prn);
	break;

    case SETMISS:
        set_miss(cmd->list, cmd->param, *pZ, pdinfo, prn);
        break;

    case LABELS:
	showlabels(pdinfo, prn);
	break;

    case VARLIST:
	varlist(pdinfo, prn);
	break;

    case PRINT:
	if (*cmd->param != '\0') {
	    do_print_string(cmd->param, prn);
	} else {
	    printdata(cmd->list, cmd->extra, (const double **) *pZ, pdinfo, 
		      cmd->opt, prn);
	}
	break;

    case PRINTF:
    case SPRINTF:
	err = do_printf(line, pZ, pdinfo, prn);
	break;

    case SSCANF:
	err = do_sscanf(line, pZ, pdinfo, prn);
	break;

    case PVALUE:
	err = batch_pvalue(line, pZ, pdinfo, prn);
	break;

    case RHODIFF:
	err = rhodiff(cmd->param, cmd->list, pZ, pdinfo);
	if (!err) {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case SUMMARY:
	err = list_summary(cmd->list, (const double **) *pZ, pdinfo, prn);
	break; 

    case XTAB:
	err = crosstab(cmd->list, (const double **) *pZ, 
		       pdinfo, cmd->opt, prn);
	break;

    case MAHAL:
	err = mahalanobis_distance(cmd->list, pZ, pdinfo, 
				   cmd->opt, prn);
	break;

    case MEANTEST:
	err = means_test(cmd->list, (const double **) *pZ, pdinfo, 
			 cmd->opt, prn);
	break;	

    case VARTEST:
	err = vars_test(cmd->list, (const double **) *pZ, pdinfo, 
			prn);
	break;

    case RUNS:
	err = runs_test(cmd->list[1], (const double **) *pZ, pdinfo, 
			cmd->opt, prn);
	break;

    case SPEARMAN:
	err = spearman(cmd->list, (const double **) *pZ, pdinfo, 
		       cmd->opt, prn);
	break;

    case DIFFTEST:
	err = diff_test(cmd->list, (const double **) *pZ, pdinfo, 
			cmd->opt, prn);
	break;

    case OUTFILE:
	err = do_outfile_command(cmd->opt, cmd->param, prn);
	break;

    case SETOBS:
	err = set_obs(line, *pZ, pdinfo, cmd->opt);
	if (!err) {
	    if (pdinfo->n > 0) {
		print_smpl(pdinfo, 0, prn);
		if (s->callback != NULL) {
		    s->callback(s, pZ, pdinfo);
		}
	    } else {
		pprintf(prn, _("setting data frequency = %d\n"), pdinfo->pd);
	    }
	}
	break;

    case SMPL:
	if (cmd->opt == OPT_F) {
	    err = restore_full_sample(pZ, pdinfo, s);
	} else if (cmd->opt) {
	    err = restrict_sample(line, cmd->list, pZ, pdinfo, 
				  s, cmd->opt, prn);
	} else { 
	    err = set_sample(line, pZ, pdinfo);
	}
	if (!err) {
	    print_smpl(pdinfo, get_full_length_n(), prn);
	}	
	break;

    case STORE:
	if (*cmd->param == '\0') {
	    pputs(prn, _("store: no filename given\n"));
	    break;
	}
	if (gretl_messages_on()) {
	    pprintf(prn, _("store: using filename %s\n"), cmd->param);
	}
	err = write_data(cmd->param, cmd->list, (const double **) *pZ,
			 pdinfo, cmd->opt, NULL);
	if (!err && gretl_messages_on()) {
	    pprintf(prn, _("Data written OK.\n"));
	}
	break;

    case SHELL:
	err = gretl_shell(line, prn);
	break;

    case OLS:
    case WLS:
    case HCCM:
	clear_model(models[0]);
	*models[0] = lsq(cmd->list, pZ, pdinfo, cmd->ci, cmd->opt);
	err = maybe_print_model(models[0], pdinfo, prn, cmd->opt);
	break;
	
#ifdef ENABLE_GMP
    case MPOLS:
	clear_model(models[0]);
	*models[0] = mp_ols(cmd->list, (const double **) *pZ, pdinfo);
	err = maybe_print_model(models[0], pdinfo, prn, cmd->opt);
	break;
#endif

    case AR:
    case ARMA:
    case ARCH:
	clear_model(models[0]);
	if (cmd->ci == AR) {
	    *models[0] = ar_func(cmd->list, pZ, pdinfo, cmd->opt, prn);
	} else if (cmd->ci == ARMA) {
	    *models[0] = arma(cmd->list, cmd->param, (const double **) *pZ, 
			      pdinfo, cmd->opt, prn);
	} else {
	    *models[0] = arch_model(cmd->list, cmd->order, pZ, pdinfo,
				    cmd->opt, prn);
	}
	err = maybe_print_model(models[0], pdinfo, prn, cmd->opt);
	break;

    case AR1:
	rho = estimate_rho(cmd->list, pZ, pdinfo, cmd->opt, prn, &err);
	if (!err) {
	    clear_model(models[0]);
	    *models[0] = ar1_lsq(cmd->list, pZ, pdinfo, cmd->ci, cmd->opt, rho);
	    err = maybe_print_model(models[0], pdinfo, prn, cmd->opt);
	}
	break;

    case ARBOND:
    case PANEL:	
	if (!dataset_is_panel(pdinfo)) {
	    strcpy(gretl_errmsg, _("This estimator requires panel data"));
	    err = E_DATA;
	    break;
	}
    case GARCH:
    case HECKIT:
    case HSK:
    case LAD:
    case LOGISTIC:
    case LOGIT:
    case POISSON:
    case PROBIT:
    case QUANTREG:
    case TOBIT:
    case TSLS:
	clear_model(models[0]);
	if (cmd->ci == LOGIT || cmd->ci == PROBIT) {
	    *models[0] = logit_probit(cmd->list, pZ, pdinfo, cmd->ci, cmd->opt, prn);
	} else if (cmd->ci == HSK) {
	    *models[0] = hsk_func(cmd->list, pZ, pdinfo);
	} else if (cmd->ci == LOGISTIC) {
	    *models[0] = logistic_model(cmd->list, pZ, pdinfo, cmd->param);
	} else if (cmd->ci == TOBIT) {
	    *models[0] = tobit_model(cmd->list, pZ, pdinfo,
				     (cmd->opt & OPT_V)? prn : NULL);
	} else if (cmd->ci == POISSON) {
	    *models[0] = poisson_model(cmd->list, pZ, pdinfo,
				       (cmd->opt & OPT_V)? prn : NULL);
	} else if (cmd->ci == HECKIT) {
	    *models[0] = heckit_model(cmd->list, pZ, pdinfo, cmd->opt, 
				     (cmd->opt & OPT_V)? prn : NULL);
	} else if (cmd->ci == TSLS) {
	    *models[0] = tsls_func(cmd->list, TSLS, pZ, pdinfo, cmd->opt);
	} else if (cmd->ci == LAD) {
	    *models[0] = lad(cmd->list, pZ, pdinfo);
	} else if (cmd->ci == QUANTREG) {
	    *models[0] = quantreg(cmd->param, cmd->list, pZ, pdinfo,
				  cmd->opt, prn);
	} else if (cmd->ci == GARCH) {
	    *models[0] = garch(cmd->list, pZ, pdinfo, cmd->opt, prn);
	} else if (cmd->ci == PANEL) {
	    *models[0] = panel_model(cmd->list, pZ, pdinfo, cmd->opt, prn);
	} else if (cmd->ci == ARBOND) {
	    *models[0] = arbond_model(cmd->list, cmd->param, (const double **) *pZ, 
				      pdinfo, cmd->opt, prn);
	} else {
	    /* can't happen */
	    err = 1;
	    break;
	}
	err = maybe_print_model(models[0], pdinfo, prn, cmd->opt);
	break;

    case GMM:
    case MLE:
    case NLS:
	err = nls_parse_line(cmd->ci, line, (const double **) *pZ, 
			     pdinfo, prn);
	if (!err) {
	    gretl_cmd_set_context(cmd, cmd->ci);
	} 
	break;

    case ADD:
    case OMIT:
    plain_add_omit:
	clear_model(models[1]);
	if (cmd->ci == ADD || cmd->ci == ADDTO) {
	    err = add_test(cmd->list, models[0], models[1], 
			   pZ, pdinfo, cmd->opt, prn);
	} else {
	    err = omit_test(cmd->list, models[0], models[1],
			    pZ, pdinfo, cmd->opt, prn);
	}
	if (!err && !(cmd->opt & OPT_Q) && !(cmd->opt & OPT_W)) {
	    /* for command-line use, we keep a stack of 
	       two models, and recycle the places */
	    swap_models(models[0], models[1]);
	}
	if (!(cmd->opt & OPT_W)) {
	    clear_model(models[1]);
	}
	if ((cmd->opt & OPT_A) && err == E_NOOMIT) {
	    /* auto-omit was a no-op */
	    err = 0;
	}
	break;	

    case ADDTO:
    case OMITFROM:
	k = atoi(cmd->param);
	if ((err = modelspec_test_check(cmd->ci, 0, k, pdinfo, prn))) {
	    break;
	}
	if (k == (models[0])->ID) {
	    goto plain_add_omit;
	} else {
	    MODEL tmpmod;

	    err = re_estimate(modelspec_get_command_by_id(k), 
			      &tmpmod, pZ, pdinfo);
	    if (err) {
		pprintf(prn, _("Failed to reconstruct model %d\n"), k);
		break;
	    } 
	    clear_model(models[1]);
	    tmpmod.ID = k;
	    if (cmd->ci == ADDTO) {
		err = add_test(cmd->list, &tmpmod, models[1], 
			       pZ, pdinfo, cmd->opt, prn);
	    } else {
		err = omit_test(cmd->list, &tmpmod, models[1],
				pZ, pdinfo, cmd->opt, prn);
	    }
	    if (err) {
		clear_model(models[1]);
		break;
	    } else {
		if (!(cmd->opt & OPT_Q)) {
		    swap_models(models[0], models[1]);
		}
		clear_model(models[1]);
	    }
	    clear_model(&tmpmod);
	}
	break;

    case COEFFSUM:
    case CUSUM:
    case RESET:
    case CHOW:
    case QLRTEST:
    case VIF:
	if (cmd->ci == COEFFSUM) {
	    err = gretl_sum_test(cmd->list, models[0], pdinfo, prn);
	} else if (cmd->ci == CUSUM) {
	    err = cusum_test(models[0], pZ, pdinfo, cmd->opt, prn);
	} else if (cmd->ci == RESET) {
	    err = reset_test(models[0], pZ, pdinfo, OPT_NONE, prn);
	} else if (cmd->ci == CHOW || cmd->ci == QLRTEST) {
	    err = chow_test(line, models[0], pZ, pdinfo, OPT_NONE, prn);
	} else if (cmd->ci == VIF) { 
	    err = vif_test(models[0], pZ, pdinfo, prn);
	} 
	break;

    case TESTUHAT:
	err = last_model_test_uhat(pZ, pdinfo, prn);
	break;

    case HAUSMAN:
	err = panel_hausman_test(models[0], pZ, pdinfo, cmd->opt, prn);
	break;

    case LMTEST:
	err = lmtest_driver(cmd->param, pZ, pdinfo, cmd->opt, prn);
	break;

    case LEVERAGE:
	err = leverage_test(models[0], pZ, pdinfo, cmd->opt, prn);
	if (!err && (cmd->opt & OPT_S)) {
	    /* FIXME gui notification? */
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case EQNPRINT:
    case TABPRINT:
	if ((models[0])->errcode == E_NAN) {
	    pprintf(prn, _("Couldn't format model\n"));
	} else {
	    char fname[FILENAME_MAX];

	    strcpy(fname, cmd->param);

	    if (cmd->opt & OPT_R) {
		err = rtfprint(models[0], pdinfo, fname, cmd->opt);
	    } else {
		err = texprint(models[0], pdinfo, fname, 
			       (cmd->ci == EQNPRINT)? (cmd->opt | OPT_E) :
			       cmd->opt);
	    }
	    if (!err) {
		pprintf(prn, _("Model printed to %s\n"), fname);
	    }
	}
	break;

    case RESTRICT:
	/* joint hypothesis test on model or system */
	if (s->rset == NULL) {
	    if (*cmd->param == '\0') {
		/* if param is non-blank, we're restricting a named system */
		err = model_test_check(cmd, pdinfo, prn);
		if (err) break;
	    }
	    s->rset = restriction_set_start(line, cmd->opt, &err);
	    if (!err) {
		gretl_cmd_set_context(cmd, RESTRICT);
	    }
	} else {
	    err = restriction_set_parse_line(s->rset, line, pdinfo);
	    if (err) {
		s->rset = NULL;
	    }	
	}
	break;

    case SYSTEM:
	/* system of equations */
	if (s->sys == NULL) {
	    s->sys = equation_system_start(line, cmd->opt, &err);
	    if (!err) {
		gretl_cmd_set_context(cmd, SYSTEM);
	    }
	} else {
	    err = system_parse_line(s->sys, line, pZ, pdinfo);
	    if (err) {
		s->sys = NULL;
	    } 
	}
	break;

    case EQUATION:
	err = equation_system_append(s->sys, cmd->list);
	if (err) {
	    s->sys = NULL;
	}
	break;

    case END:
	if (!strcmp(cmd->param, "system")) {
	    err = equation_system_finalize(s->sys, pZ, pdinfo, cmd->opt, prn);
	    if (err || s->sys->name == NULL) {
		s->sys = NULL;
	    } else {
		system_set_save_flag(s->sys);
	    }
	} else if (!strcmp(cmd->param, "mle") || 
		   !strcmp(cmd->param, "nls") ||
		   !strcmp(cmd->param, "gmm")) {
	    clear_model(models[0]);
	    *models[0] = nls(pZ, pdinfo, cmd->opt, prn);
	    err = maybe_print_model(models[0], pdinfo, prn, cmd->opt);
	    if (!err) {
		s->alt_model = 1;
	    }
	} else if (!strcmp(cmd->param, "restrict")) {
	    err = do_end_restrict(s, pZ, pdinfo);
	} else {
	    err = 1;
	}
	break;

    case VAR:
    case VECM:
	if (cmd->ci == VAR) {
	    s->var = gretl_VAR(cmd->order, cmd->list, (const double **) *pZ, 
			       pdinfo, cmd->opt, prn, &err);
	} else {
	    s->var = gretl_VECM(cmd->order, cmd->aux, cmd->list, (const double **) *pZ,
				pdinfo, cmd->opt, prn, &err);
	}
	if (s->var != NULL) {
	    if (s->callback != NULL) {
		s->callback(s, pZ, pdinfo);
	    }
	}
	break;

    case RUN:
    case INCLUDE:
	err = getopenfile(line, runfile, NULL, OPT_S);
	if (err) { 
	    break;
	} 
	if (gretl_messages_on()) {
	    pprintf(prn, " %s\n", runfile);
	}
	if (cmd->ci == INCLUDE && gretl_is_xml_file(runfile)) {
	    err = load_user_matrix_file(runfile);
	    break;
	}
	if (!strcmp(runfile, s->runfile)) { 
	    pprintf(prn, _("Infinite loop detected in script\n"));
	    err = 1;
	    break;
	}
	err = run_script(runfile, s, pZ, pdinfo, prn);
	break;

    case FUNCERR:
	err = s->funcerr = 1;
	break;

    default:
	pprintf(prn, _("Sorry, the %s command is not yet implemented "
		       "in libgretl\n"), cmd->word);
	err = 1;
	break;
    }

    if (listcpy != NULL) {
	free(listcpy);
    }

    if (err == E_OK) {
	err = 0;
    }

    if (err) {
	gretl_cmd_destroy_context(cmd);
	if (!s->funcerr) {
	    errmsg(err, prn);
	}
    }

    return err;
}


/* called by functions, and by scripts executed from within
   functions */

int maybe_exec_line (ExecState *s, double ***pZ, DATAINFO *pdinfo)
{
    int err = 0;

    if (string_is_blank(s->line)) {
	return 0;
    }

    if (gretl_compiling_loop()) { 
	err = get_command_index(s->line, s->cmd, pdinfo);
    } else {
	err = parse_command_line(s->line, s->cmd, pZ, pdinfo);
    }

    if (err) {
        errmsg(err, s->prn);
        return 1;
    }
    
    s->in_comment = cmd_ignore(s->cmd)? 1 : 0;

    if (s->cmd->ci < 0) {
	return 0; /* nothing there, or a comment */
    }

    if (s->cmd->ci == LOOP || gretl_compiling_loop()) {  
	/* accumulating loop commands */
	if (!ok_in_loop(s->cmd->ci)) {
            pputs(s->prn, _("Sorry, this command is not available in loop mode\n"));
            return 1;
        }
	err = gretl_loop_append_line(s, pZ, pdinfo);
	if (err) {
	    errmsg(err, s->prn);
	    return 1;
	} 
	return 0;
    } 

    if (s->cmd->ci == FUNCERR) {
	s->funcerr = err = 1;
    } else {
	err = gretl_cmd_exec(s, pZ, pdinfo);
    }

    /* FIXME: problem with "set_as_last_model" below */

    if (!err && (is_model_cmd(s->cmd->word) || s->alt_model)
	&& !is_quiet_model_test(s->cmd->ci, s->cmd->opt)) {
	attach_subsample_to_model(s->models[0], pdinfo);
	set_as_last_model(s->models[0], GRETL_OBJ_EQN);
    }

    if (system_save_flag_is_set(s->sys)) {
	/* only warrants action in GUI program */
	system_unset_save_flag(s->sys);
	s->sys = NULL;
    }

    return err;
}

static int could_be_varname (const char *s)
{
    int n = gretl_varchar_spn(s);
    char word[VNAMELEN];

    if (n > 0 && n < VNAMELEN) {
	*word = '\0';
	strncat(word, s, n);
	if (check_varname(word) == 0) {
	    return 1;
	}
    }

    return 0;
}

/**
 * get_command_index:
 * @line: command line.
 * @cmd: pointer to gretl command struct.
 * @pdinfo: dataset information.
 *
 * Parse @line and assign to the %ci field of @cmd the index number of
 * the command embedded in @line.  Note: this is a "lite" version of
 * parse_command_line().  It is used when commands are being stacked
 * for execution within a loop.  Command options are not parsed out of
 * @line.
 *
 * Returns: 1 on error, otherwise 0.
 */

int get_command_index (char *line, CMD *cmd, const DATAINFO *pdinfo)
{
    static int context;
    int done = 0;

    cmd->ci = 0;
    cmd->opt = OPT_NONE;

    while (isspace(*line)) {
	line++;
    }

#if CMD_DEBUG
    fprintf(stderr, "get_command_index: line='%s'\n", line);
#endif

    if (filter_comments(line, cmd)) {
	return 0;
    }

    if (!get_command_word(line, cmd)) {
	if (*line == '$') {
	    /* most plausible possibility? */
	    strcpy(cmd->word, "genr");
	    cmd->ci = GENR;
	} else {
	    cmd_set_nolist(cmd);
	    cmd->ci = CMD_NULL;
#if CMD_DEBUG
	    fprintf(stderr, "get_command_index: got nothing, returning 0\n");
#endif
	    return 0;
	}
    }

#if CMD_DEBUG
    fprintf(stderr, " got command word = '%s'\n", cmd->word);
#endif

    /* subsetted commands (e.g. "deriv" in relation to "nls") */
    if (!strcmp(cmd->word, "end")) {
	context = 0;
	cmd->ci = END;
	done = 1;
    } else if (context && strcmp(cmd->word, "equation")) {
	/* "equation" occurs in the SYSTEM context, but it is
	   a command in its own right, so we don't set cmd->ci
	   to the context value */
	cmd->ci = context;
#if CMD_DEBUG
	fprintf(stderr, " context (static) = %d, ci = %d\n", context, cmd->ci);
#endif
	done = 1;
    } else if (catch_command_alias(line, cmd)) {
#if CMD_DEBUG
	fprintf(stderr, " caught command alias, ci = %d\n", cmd->ci);
#endif
	done = 1;
    } 

    if (!done) {
	cmd->ci = gretl_command_number(cmd->word);
#if CMD_DEBUG
	fprintf(stderr, " gretl_command_number(%s) gave %d\n", cmd->word, cmd->ci);
#endif
	if (cmd->ci == 0) {
	    if (could_be_varname(line)) {
		cmd->ci = GENR;
	    } else if (get_user_function_by_name(cmd->word)) {
		cmd->ci = GENR;
		cmd->opt = OPT_U;
	    } else {
		cmd->err = 1;
		sprintf(gretl_errmsg, _("command '%s' not recognized"), 
			cmd->word);
		return 1;
	    }
	}
    }	

    if (cmd->ci == NLS) {
	context = NLS;
    } else if (cmd->ci == MLE) {
	context = MLE;
    } else if (cmd->ci == GMM) {
	context = GMM;
    }

    if (!strcmp(line, "end loop")) {
	cmd->ci = ENDLOOP;
    }

#if CMD_DEBUG
    fprintf(stderr, " cmd->ci set to %d\n", cmd->ci);
#endif

    return 0;
}

/* which commands can we run without having first opened 
   a data file?
*/

#define startup_ci(c) (c == OPEN || \
		       c == RUN || \
		       c == INCLUDE || \
		       c == NULLDATA || \
		       c == GENR || \
		       c == PVALUE || \
		       c == PRINT || \
		       c == PRINTF || \
		       c == HELP || \
		       c == SET || \
		       c == FUNC || \
		       c == QUIT || \
		       c == LOOP || \
		       c == ENDLOOP || \
		       c == IF || \
		       c == ELIF || \
		       c == ELSE || \
		       c == ENDIF)

int ready_for_command (const char *line)
{
    const char *ok_cmds[] = {
	"import", 
	"delete",
	"eval",
	"!",
	"launch",
	"matrix",
	"scalar",
	"string",
	NULL 
    };
    int i, ok = 0;

    if (string_is_blank(line) || gretl_compiling_function()) {
	ok = 1;
    } else if (*line == '#') {
	ok = 1;
    } else if (*line == '/' && *(line+1) == '*') {
	ok = 1;
    } else {
	char word[9];
	int ci;

	sscanf(line, "%8s", word);
	ci = gretl_command_number(word);
	if (ci == 0 || startup_ci(ci)) {
	    ok = 1;
	} else {
	    for (i=0; ok_cmds[i] != NULL && !ok; i++) {
		if (!strncmp(line, ok_cmds[i], strlen(ok_cmds[i]))) {
		    ok = 1;
		}
	    }
	}
    }

    return ok;
}

int gretl_cmd_init (CMD *cmd)
{
    cmd->ci = 0;
    cmd->err = 0;
    cmd->context = 0;
    cmd->order = 0;
    cmd->aux = 0;
    cmd->flags = 0;
    *cmd->word = '\0';

    cmd->list = NULL;
    cmd->param = NULL;
    cmd->extra = NULL;
    cmd->linfo = NULL;

    /* make 'list', 'param' and 'extra' blank rather than NULL
       for safety (in case they are deferenced) */

    cmd->list = gretl_null_list();
    if (cmd->list == NULL) {
	cmd->err = E_ALLOC;
    }

    if (cmd->err == 0) {
	cmd->param = calloc(1, 1);
	if (cmd->param == NULL) {
	    cmd->err = E_ALLOC;
	}
    }

    if (cmd->err == 0) {
	cmd->extra = calloc(1, 1);
	if (cmd->extra == NULL) {
	    free(cmd->param);
	    cmd->param = NULL;
	    cmd->err = E_ALLOC;
	}
    }    

    return cmd->err;
}

void gretl_cmd_free (CMD *cmd)
{
    free(cmd->list);
    free(cmd->param);
    free(cmd->extra);

    cmd_lag_info_destroy(cmd);
}

void gretl_cmd_destroy (CMD *cmd)
{
    gretl_cmd_free(cmd);
    free(cmd);
}

CMD *gretl_cmd_new (void)
{
    CMD *cmd = malloc(sizeof *cmd);

    if (cmd != NULL) {
	gretl_cmd_init(cmd);
    }

    return cmd;
}

void gretl_cmd_set_context (CMD *cmd, int ci)
{
    cmd->context = ci;
}

void gretl_cmd_destroy_context (CMD *cmd)
{
    cmd->context = 0;
}

gretlopt gretl_cmd_get_opt (const CMD *cmd)
{
    return cmd->opt;
}

void gretl_cmd_set_opt (CMD *cmd, gretlopt opt)
{
    cmd->opt = opt;
}

char *gretl_cmd_get_savename (char *sname)
{
    strcpy(sname, cmd_savename);
    *cmd_savename = 0;

    return sname;
}

void gretl_exec_state_init (ExecState *s,
			    ExecFlags flags,
			    char *line,
			    CMD *cmd,
			    MODEL **models, 
			    PRN *prn)
{
    s->flags = flags;
    s->line = line;
    s->cmd = cmd;

    *s->runfile = '\0';

    s->models = models;
    s->prn = prn;

    s->sys = NULL;
    s->rset = NULL;
    s->var = NULL;
    s->alt_model = 0;
    s->in_comment = 0;
    s->funcerr = 0;

    s->submask = NULL;
    s->callback = NULL;
}

void gretl_exec_state_clear (ExecState *s)
{
    gretl_cmd_free(s->cmd);
    destroy_working_models(s->models, 2);
    set_as_last_model(NULL, GRETL_OBJ_NULL);
    free_subsample_mask(s->submask);
    s->funcerr = 0;
}

void gretl_exec_state_uncomment (ExecState *s)
{
    s->in_comment = 0;
    s->cmd->flags &= ~CMD_IGNORE;
}

/* if-then stuff - conditional execution */

static int if_eval (const char *s, double ***pZ, DATAINFO *pdinfo, int *err)
{
    double val = NADBL;
    int ret = -1;

#if IFDEBUG
    fprintf(stderr, "if_eval: line = '%s'\n", s);
#endif

    if (!strncmp(s, "if", 2)) {
	s += 2;
    } else if (!strncmp(s, "elif", 4)) {
	s += 4;
    }

    while (*s == ' ') s++;

    val = generate_scalar(s, pZ, pdinfo, err);

#if IFDEBUG
    if (err) {
	fprintf(stderr, "if_eval: generate returned %d\n", *err);
    }
#endif

    if (*err) {
	gretl_errmsg_set(_("error evaluating 'if'"));
    } else if (na(val)) {
	*err = 1;
	strcpy(gretl_errmsg, _("indeterminate condition for 'if'"));
    } else {
	ret = (int) val;
    }

#if IFDEBUG
    fprintf(stderr, "if_eval: returning %d\n", ret);
#endif

    return ret;
}

#if IFDEBUG
static const char *ifstr (int c)
{
    if (c == SET_FALSE) return "SET_FALSE";
    if (c == SET_TRUE)  return "SET_TRUE";
    if (c == SET_ELSE)  return "SET_ELSE";
    if (c == SET_ELIF)  return "SET_ELIF";
    if (c == SET_ENDIF) return "SET_ENDIF";
    if (c == IS_FALSE)  return "IS_FALSE";
    if (c == IS_TRUE)   return "IS_TRUE";
    if (c == UNINDENT)  return "UNINDENT";
    if (c == GETINDENT) return "GETINDENT";
    if (c == RELAX)     return "RELAX";
    return "UNKNOWN";
}
#endif

static void unmatched_message (int code)
{
    sprintf(gretl_errmsg, _("Unmatched \"%s\""),
	    (code == SET_ELSE)? "else" : 
	    (code == SET_ELIF)? "elif": "endif");
}

#define IF_DEPTH 32

static int ifstate (int code, int *err)
{
    static unsigned char T[IF_DEPTH];
    static unsigned char got_if[IF_DEPTH];
    static unsigned char got_else[IF_DEPTH];
    static unsigned char got_T[IF_DEPTH];
    static unsigned char indent;
    int i, ret = 0;

#if IFDEBUG
    fprintf(stderr, "ifstate: code = %s\n", ifstr(code));
#endif

    if (code == RELAX) {
	indent = 0;
    } else if (code == GETINDENT) {
	ret = indent;
    } else if (code == UNINDENT) {
	ret = --indent;
    } else if (code == SET_FALSE || code == SET_TRUE) {
	indent++;
	if (indent >= IF_DEPTH) {
	    sprintf(gretl_errmsg, "IF depth (%d) exceeded\n", IF_DEPTH);
	    *err = E_DATA;
	} else {
	    T[indent] = got_T[indent] = (code == SET_TRUE);
	    got_if[indent] = 1;
	    got_else[indent] = 0;
	}
    } else if (code == SET_ELSE || code == SET_ELIF) {
	if (got_else[indent] || !got_if[indent]) {
	    unmatched_message(code);
	    *err = E_PARSE;
	} else {
	    got_else[indent] = (code == SET_ELSE);
	    if (T[indent]) {
		T[indent] = 0;
	    } else if (!got_T[indent]) {
		T[indent] = 1;
	    }
	}
    } else if (code == SET_ENDIF) {
	if (!got_if[indent] || indent == 0) {
	    unmatched_message(code);
	    *err = E_PARSE;
	} else {
	    got_if[indent] = 0;
	    got_else[indent] = 0;
	    got_T[indent] = 0;
	    indent--;
	}
    } else if (code == IS_FALSE || code == IS_TRUE) {
	for (i=1; i<=indent; i++) {
	    if (T[i] == 0) {
		ret = 1;
		break;
	    }
	}
	if (code == IS_TRUE) {
	    ret = !ret;
	}
    } 

#if IFDEBUG
    fprintf(stderr, "ifstate: returning %d (indent %d, err %d)\n", 
	    ret, indent, (err == NULL)? 0 : *err);
#endif

    return ret;
}

static int set_if_state (int code)
{
    int err = 0;

    ifstate(code, &err);
    return err;
}

static int get_if_state (int code)
{
    int err = 0;

    return ifstate(code, &err);
}

void gretl_if_state_clear (void)
{
    ifstate(RELAX, NULL);
}

int gretl_if_state_finalize (void)
{
    int ret, err = 0;

    ret = ifstate(IS_TRUE, NULL);

    if (!ret) {
	ifstate(RELAX, NULL);
	err = E_PARSE;
    }

    return err;
}

int gretl_if_state_record (void)
{
    return ifstate(GETINDENT, NULL);
}

int gretl_if_state_check (int indent0)
{
    int indent, err = 0;

    indent = ifstate(GETINDENT, NULL);

    if (indent != indent0) {
	sprintf(gretl_errmsg, _("Unmatched \"%s\""), "if");
	ifstate(RELAX, NULL);
	err = E_PARSE;
    }

    return err;
}
