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
#include "uservar.h"
#include "gretl_panel.h"
#include "texprint.h"
#include "gretl_xml.h"
#include "gretl_string_table.h"
#include "dbread.h"
#include "gretl_foreign.h"
#include "boxplots.h"
#include "kalman.h"
#include "flow_control.h"
#include "libglue.h"
#include "csvdata.h"
#ifdef USE_CURL
# include "gretl_www.h"
#endif

#include <errno.h>
#include <glib.h>

/* for the "shell" command */
#ifndef WIN32
# ifdef HAVE_PATHS_H
#  include <paths.h>
# endif
#endif

#define CMD_DEBUG 0
#define LAGS_DBG 0

#include "laginfo.c"

typedef struct {
    int v;                /* ID number of series */
    int *vlist;           /* or list of series */
    char name[VNAMELEN];  /* name of series of list */
    int lmin;             /* minimum lag */
    int lmax;             /* maximum lag */
    int *laglist;         /* list of specific lags */
} LAGVAR;

#define cmd_set_nolist(c) (c->flags |= CMD_NOLIST)
#define cmd_unset_nolist(c) (c->flags &= ~CMD_NOLIST)

static void get_optional_filename (const char *s, CMD *cmd);

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

int filter_comments (char *s, CMD *cmd)
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

/* as in, e.g., "lags 4 ; <varlist>" but we allow the leading
   parameter to be given as a pre-defined scalar variable
   instead of a numeric constant
*/

static int get_lags_param (CMD *cmd, char **ps)
{
    int k = gretl_charpos(';', *ps);
    int ret = 0;

    if (k > 0) {
	char *tmp = gretl_strndup(*ps, k-1);
	int lag = positive_int_from_string(tmp);

	if (lag < 0 && gretl_is_scalar(tmp)) {
	    lag = gretl_scalar_get_value(tmp, NULL);
	    free(tmp);
	    tmp = g_strdup_printf("%d", lag);
	}

	free(cmd->param);
	cmd->param = tmp;
	*ps = *ps + k + 1;
	ret = 1;
    }

    return ret;
}

static void deprecate_alias (const char *bad, const char *good, int cmd)
{
    const char *tag = cmd ? "command" : "construction";

    gretl_warnmsg_sprintf("\"%s\": obsolete %s; please use \"%s\"",
			  bad, tag, good);
}

/* catch aliased command words and assign cmd->ci; return
   cmd->ci if alias caught, else 0. */

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
    } else if (!strcmp(s, "pooled")) {
	deprecate_alias("pooled", "ols", 1);
	cmd->ci = OLS;
    } else if (!strcmp(line, "smpl full")) {
	strcpy(line, "smpl");
	cmd->opt = OPT_F;
    } else if (!strcmp(s, "equations")) {
	cmd->ci = EQUATION;
	cmd->opt |= OPT_M;
    } else if (!strcmp(s, "graph")) { 
	deprecate_alias("graph", "textplot", 1);
	cmd->ci = PLOT; 	 
    } else if (!strcmp(s, "plot")) {
	deprecate_alias("plot", "textplot", 1);
	cmd->ci = PLOT; 	 
	cmd->opt = OPT_S;
    } else if (!strcmp(s, "list")) {
	char lname[VNAMELEN];
	char fmt[24];

	if (string_is_blank(line + 4)) {
	    cmd->ci = VARLIST;
	    strcpy(line, "varlist");
	} else if (gretl_string_ends_with(line, "delete")) {
	    sprintf(fmt, "list %%%ds delete", VNAMELEN - 1);
	    if (sscanf(line, fmt, lname)) {
		free(cmd->parm2);
		cmd->parm2 = gretl_strdup(lname);
		cmd->ci = DELEET;
	    }
	} else {
	    if (gretl_string_ends_with(line, "print")) {
		sprintf(fmt, "list %%%ds", VNAMELEN - 1);
		if (sscanf(line, fmt, lname)) {
		    strcpy(line, lname);
		}
	    } 
	    cmd->ci = GENR;
	}
    } else if (*s == '!' || !strcmp(s, "launch")) {
	cmd->ci = SHELL;
    } else if (!strcmp(s, "addobs")) { 	 
	char *tmp = gretl_strdup(line); 	 

	deprecate_alias("addobs", "dataset addobs", 1);
	strcpy(line, "dataset "); 	 
	strcat(line, tmp); 	 
	cmd->ci = DATAMOD; 	 
	free(tmp);
    } else if (!strcmp(s, "fcasterr")) {
	deprecate_alias("fcasterr", "fcast", 1);
	cmd->ci = FCAST; 	 
    } else if (!strcmp(s, "continue")) {
	cmd->ci = FUNDEBUG;
	cmd->opt |= OPT_C;
    } else if (!strcmp(s, "next")) {
	cmd->ci = FUNDEBUG;
	cmd->opt |= OPT_N;
    } else if (!strcmp(s, "undebug")) {
	cmd->ci = FUNDEBUG;
	cmd->opt |= OPT_Q;
    }

    return cmd->ci;
}

static int catch_system_alias (CMD *cmd)
{
    cmd->ci = 0;
    
    if (!strcmp(cmd->word, "equation")) {
	cmd->ci = EQUATION;
    } else if (!strcmp(cmd->word, "equations")) {
	cmd->ci = EQUATION;
	cmd->opt |= OPT_M;
    }

    return cmd->ci;
}

#define REQUIRES_PARAM(c) (c == DATAMOD ||	\
                           c == FUNC ||		\
                           c == LOOP ||		\
                           c == MAKEPKG ||	\
			   c == MODPRINT ||	\
                           c == NULLDATA ||	\
                           c == SETMISS)

#define REQUIRES_ORDER(c) (c == ADF ||		\
                           c == ARCH ||		\
                           c == COINT ||	\
                           c == COINT2 ||	\
                           c == KPSS ||		\
			   c == LEVINLIN ||	\
                           c == NULLDATA ||	\
                           c == VAR ||		\
                           c == VECM)

#define NO_VARLIST(c) (c == APPEND ||		\
                       c == BREAK ||		\
                       c == CHOW ||		\
		       c == CLEAR ||		\
	               c == CUSUM ||		\
                       c == DATA ||		\
                       c == END ||		\
	               c == ENDLOOP ||		\
                       c == ESTIMATE ||		\
	               c == EQNPRINT ||		\
	               c == FCAST ||		\
		       c == FLUSH ||		\
		       c == FOREIGN ||		\
                       c == FUNC ||		\
                       c == FUNCERR ||		\
                       c == FUNCRET ||		\
		       c == FUNDEBUG ||		\
	               c == GENR ||		\
                       c == GMM ||		\
		       c == GRAPHPG ||		\
	               c == HAUSMAN ||		\
                       c == HELP ||		\
                       c == INCLUDE ||		\
    	               c == INFO ||		\
		       c == JOIN ||		\
                       c == KALMAN ||		\
                       c == LEVERAGE ||		\
                       c == LOOP ||		\
		       c == MAKEPKG ||		\
		       c == MARKERS ||		\
                       c == MLE ||		\
                       c == MODELTAB ||		\
                       c == MODPRINT ||		\
                       c == MODTEST ||		\
		       c == MPI ||		\
                       c == NLS ||		\
                       c == NULLDATA ||		\
		       c == OPEN ||		\
                       c == OUTFILE ||		\
                       c == PRINTF ||		\
	               c == PVAL ||		\
                       c == QLRTEST ||		\
	               c == QUIT ||		\
                       c == RENAME ||		\
                       c == RESET ||		\
                       c == RESTRICT ||		\
	               c == RUN ||		\
                       c == SET ||		\
                       c == SETINFO ||		\
	               c == SETOBS ||		\
	               c == SHELL ||		\
                       c == SPRINTF ||		\
		       c == SSCANF ||		\
                       c == SYSTEM ||		\
                       c == TABPRINT ||		\
                       c == VARLIST ||		\
                       c == VIF)

#define USES_LISTSEP(c) (c == AR ||		\
                         c == ARBOND ||		\
                         c == ARMA ||		\
                         c == BIPROBIT ||	\
                         c == COINT2 ||		\
			 c == DPANEL ||		\
                         c == EQUATION ||	\
                         c == GARCH ||		\
                         c == HECKIT ||		\
                         c == IVREG ||		\
                         c == MPOLS ||		\
                         c == POISSON ||	\
			 c == NEGBIN ||		\
			 c == DURATION ||	\
                         c == PRINT ||		\
                         c == SCATTERS ||	\
                         c == VAR ||		\
                         c == VECM ||		\
                         c == XTAB)

#define DOUBLE_SEP_OK(c) (c == ARBOND ||	\
                          c == DPANEL ||	\
                          c == ARMA ||		\
                          c == COINT2 ||	\
			  c == VECM) 

#define NEEDS_LISTSEP(c) (c == AR ||		\
                          c == ARBOND ||	\
                          c == ARMA ||		\
			  c == DPANEL ||	\
                          c == GARCH ||		\
                          c == HECKIT ||	\
                          c == IVREG)

#define DEFAULTS_TO_FULL_LIST(c) (c == CORR ||		\
                                  c == DIFF ||		\
                                  c == LDIFF ||		\
                                  c == LABELS ||	\
                                  c == LAGS ||		\
                                  c == LOGS ||		\
                                  c == PCA ||		\
                                  c == PRINT ||		\
                                  c == SDIFF ||		\
                                  c == SMPL ||		\
                                  c == SQUARE ||	\
                                  c == ORTHDEV ||	\
                                  c == STORE ||		\
                                  c == SUMMARY)

#define MODIFIES_LIST(c) (c == DIFF ||		\
			  c == DUMMIFY ||	\
			  c == LDIFF ||		\
			  c == SDIFF ||		\
			  c == LAGS ||		\
			  c == LOGS ||		\
			  c == SQUARE ||	\
	                  c == ORTHDEV)

/* given an assignment such as "foo <- command", extract
   the first field and record it in the "savename"
   member of CMD.
*/

static void maybe_extract_savename (char *s, CMD *cmd)
{
    char savename[32], test[4];
    int n;

    *cmd->savename = '\0';

#if CMD_DEBUG > 1
    fprintf(stderr, "** testing for savename: '%s'\n", s);
#endif

    if (*s == '"') {
	n = sscanf(s, "\"%31[^\"]\" <- %3s", savename, test);
    } else {
	n = sscanf(s, "%31s <- %3s", savename, test);
    } 

    if (n == 2) {
	char *p = strstr(s + strlen(savename), " <- ");

	strcpy(cmd->savename, savename);
	p += 4;
	p += strspn(p, " ");
	n = p - s;
	shift_string_left(s, n);
    }

#if CMD_DEBUG > 1
    fprintf(stderr, "** after test for savename: '%s'\n", s);
#endif
}

static const char *maybe_skip_savename (const char *s)
{
    char savename[32], test[4];
    int n;

    if (*s == '"') {
	n = sscanf(s, "\"%31[^\"]\" <- %3s", savename, test);
    } else {
	n = sscanf(s, "%31s <- %3s", savename, test);
    } 

    if (n == 2) {
	s = strstr(s + strlen(savename), " <- ");
	s += 4;
    }

    return s;
}

static inline void maybe_set_catch_flag (char *s, CMD *cmd)
{
    if (strncmp(s, "catch ", 6) == 0) {
	/* the next two lines added 2012-12-27 */
	set_gretl_errno(0);
	gretl_error_clear();
	cmd->flags |= CMD_CATCH;
	shift_string_left(s, 6);
    } else if (!cmd->context) {
	cmd->flags &= ~CMD_CATCH;
    }
}

/* grab a filename, possibly prepending userdir */

static int filename_to_param (CMD *cmd, const char *s, 
			      int *len, int *quoted)
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

    if (libset_get_bool(USE_CWD) || fname_has_path(fname)) {
	cmd->param = fname;
    } else if (cmd->ci == OUTFILE && !strcmp(fname, "null")) {
	cmd->param = fname;
    } else {
	cmd->param = gretl_strdup_printf("%s%s", gretl_workdir(), fname);
	free(fname);
	if (cmd->param == NULL) {
	    return E_ALLOC;
	}
    }

    return 0;
}

static int get_maybe_quoted_filename (CMD *cmd, char **ps)
{
    int err, len = 0;
    int quoted = 0;

    err = filename_to_param(cmd, *ps, &len, &quoted);
    if (err) {
	return err;
    }

    *ps = *ps + len + 2 * quoted;

    return 0;
} 

static void grab_gnuplot_literal_block (char *s, CMD *cmd)
{
    s = strchr(s, '{');
    if (s != NULL) {
	free(cmd->param);
	cmd->param = gretl_strdup(s);
	*s = '\0';
    }
}

static int boxplot_booleans_present (const char *s)
{
    int ret = 0;

    /* Note: the '(' character might appear inside a gnuplot
       literal block, "{...}", in which case it should not be 
       taken as indicating that boolean conditions are present
       on the command line.
    */

    if (strchr(s, '(')) {
	int braced = 0;

	while (*s && !ret) {
	    if (*s == '{') {
		braced++;
	    } else if (*s == '}') {
		braced--;
	    } else if (!braced && *s == '(') {
		ret = 1;
	    }
	    s++;
	}
    }

    return ret;
}

static int is_special_lag_field (const char *s)
{
    if (*s == '{') {
	return 1;
    } else if (isalpha(*s) && gretl_is_matrix(s)) {
	return 1;
    } else {
	return 0;
    }
}

/* Check that any "{xxx}" or matrix fields are in appropriate 
   positions; also see if there _are_ any of the latter.
*/

static int gappy_lags_check (char **S, int ns, int *specials, int ci)
{
    int i, err = 0;

    if (ci == ARMA && ns != 2 && ns != 3) {
	err = E_PARSE;
    } else if (ci == DPANEL && ns != 1) {
	err = E_PARSE;
    } else if (ci == VAR && ns != 1) {
	err = E_PARSE;
    }

    for (i=0; i<ns && !err; i++) {
	if (is_special_lag_field(S[i])) {
	    if (i != 0 && i != ns - 1) {
		err = E_PARSE;
	    } else if (i == 0) {
		specials[0] = 1;
	    } else {
		specials[1] = i;
	    }
	}
    }

    return err;
}

/* push onto array *pS a string of length len starting at s */

static int push_lag_field (char ***pS, const char *s, int len, int *ns)
{
    char *chunk = gretl_strndup(s, len);
    int err = 0;

    if (chunk == NULL) {
	err = E_ALLOC;
    } else {
	err = strings_array_add(pS, ns, chunk);
	free(chunk);
    }

    return err;
}

#define no_specials(s) (s[0] == 0 && s[1] == 0)

/* For some commands -- notably ar(i)ma -- we allow for
   a lag specification taking the form of integers in
   braces or a named matrix, permitting "gappy" lags.
   Here we split such specifications into their
   components.
*/

static char **split_lag_fields (char *s, int *ns, 
				int *specials, CMD *cmd,
				char **rem)
{
    char *q, *p = s;
    char **S = NULL;
    int n;

    if (cmd->ci == VAR) {
	while (*p == ' ') p++;
    }

    while (*p && !cmd->err) {
	if (cmd->ci == VAR && *p == ' ') {
	    *rem = p;
	    break;
	}		
	while (*p == ' ') p++;
	if (*p == ';') {
	    /* reached the end of the portion of the command line
	       subject to special treatment */
	    *rem = p;
	    break;
	} else if (*p == '{') {
	    q = strchr(p, '}');
	    if (q == NULL) {
		cmd->err = E_PARSE;
	    } else {
		n = strcspn(p, "}");
		cmd->err = push_lag_field(&S, p, n + 1, ns);
		if (cmd->ci == VAR) {
		    *rem = q + 1;
		    break;
		}
		p = q;
	    }
	} else {
	    n = strcspn(p, " {};");
	    if (n == 0) {
		cmd->err = E_PARSE;
	    } else {
		cmd->err = push_lag_field(&S, p, n, ns);
		p += n - 1;
	    }
	}
	p++;
    }

    if (!cmd->err) {
	cmd->err = gappy_lags_check(S, *ns, specials, cmd->ci);
    }

    if (cmd->err || no_specials(specials)) {
	strings_array_free(S, *ns);
	S = NULL;
    }

    return S;
}

/* here we have only one field to worry about, specifying the 
   lag pattern for the dependent variable(s) */

static void handle_single_lagvec (CMD *cmd, const char *lspec,
				  char *line, const char *rem)
{
    if (*lspec == '{') {
	cmd->auxlist = gretl_list_from_string(lspec, &cmd->err);
    } else {
	gretl_matrix *m = get_matrix_by_name(lspec);

	if (m == NULL) {
	    cmd->err = E_UNKVAR;
	} else {
	    cmd->auxlist = gretl_list_from_vector(m, &cmd->err);
	}
    }

    if (!cmd->err) {
	int lmin = 0, lmax = 0;

	cmd->err = gretl_list_min_max(cmd->auxlist, &lmin, &lmax);
	if (!cmd->err && lmin < 1) {
	    cmd->err = E_DATA;
	}

	if (!cmd->err) {
	    char *tmp, numstr[16];

	    while (*rem == ' ') rem++;
	    /* we have to copy here since @rem is actually
	       a pointer into @line */
	    tmp = gretl_strdup(rem);

	    if (tmp == NULL) {
		cmd->err = E_ALLOC;
	    } else {
		sprintf(numstr, "%d ", lmax);
		*line = '\0';
		strcat(line, numstr);
		strcat(line, tmp);
		free(tmp);
	    }
	}
    }
}

/* we may have a "special" (gappy) lag spec for the AR term(s) 
   and/or the MA terms
*/

static void handle_arma_lags (CMD *cmd, char **S,
			      int ns, int *specials,
			      char *line, char *rem)
{
    int *plist = NULL, *qlist = NULL;
    const char *lspec;
    char *tmp = NULL;

    if (specials[0]) {
	lspec = S[0];
	if (*lspec == '{') {
	    plist = gretl_list_from_string(lspec, &cmd->err);
	} else {
	    gretl_matrix *m = get_matrix_by_name(lspec);

	    if (m == NULL) {
		cmd->err = E_UNKVAR;
	    } else {
		plist = gretl_list_from_vector(m, &cmd->err);
	    }
	}
    }

    if (specials[1] && !cmd->err) {
	lspec = S[ns-1];
	if (*lspec == '{') {
	    qlist = gretl_list_from_string(lspec, &cmd->err);
	} else {
	    gretl_matrix *m = get_matrix_by_name(lspec);

	    if (m == NULL) {
		cmd->err = E_UNKVAR;
	    } else {
		qlist = gretl_list_from_vector(m, &cmd->err);
	    }
	}
    }

    if (!cmd->err) {
	/* form the full list to pass to arma */
	if (plist != NULL && qlist == NULL) {
	    cmd->auxlist = plist;
	} else if (qlist != NULL) {
	    cmd->auxlist = gretl_lists_join_with_separator(plist, qlist);
	    if (cmd->auxlist == NULL) {
		cmd->err = E_ALLOC;
	    }	    
	} 
    }

    if (!cmd->err) {
	tmp = gretl_strdup(rem);
	if (tmp == NULL) {
	    cmd->err = E_ALLOC;
	} else {
	    *line = '\0';
	}
    }	

    if (!cmd->err) {
	int lmin = 0, lmax = 0;
	char numstr[16];

	if (plist != NULL) {
	    cmd->err = gretl_list_min_max(plist, &lmin, &lmax);
	    if (!cmd->err && lmin < 1) {
		cmd->err = E_DATA;
	    }
	    if (!cmd->err) {
		sprintf(numstr, "%d ", lmax);
		strcat(line, numstr);
	    }
	} else {
	    sprintf(numstr, "%s ", S[0]);
	    strcat(line, numstr);
	}

	if (ns == 3) {
	    /* ARIMA d spec */
	    sprintf(numstr, "%s ", S[1]);
	    strcat(line, numstr);
	}

	if (qlist != NULL) {
	    cmd->err = gretl_list_min_max(qlist, &lmin, &lmax);
	    if (!cmd->err && lmin < 1) {
		cmd->err = E_DATA;
	    }	    
	    if (!cmd->err) {
		sprintf(numstr, "%d ", lmax);
		strcat(line, numstr);
	    }
	} else {
	    sprintf(numstr, "%s ", S[ns-1]);
	    strcat(line, numstr);
	}

	strcat(line, tmp);
    }

    free(tmp);

    if (plist != cmd->auxlist) {
	free(plist);
    }

    if (qlist != cmd->auxlist) {
	free(qlist);
    }    
}

static int maybe_rewrite_lags (char *s, CMD *cmd)
{
    char **S = NULL;
    char *line = s;   /* save the starting point */
    char *rem = NULL; /* will point to remainder of line */
    int specials[2] = {0};
    int ns = 0;

    if (!strncmp(line, "arma ", 5)) {
	s += 5;
    } else if (!strncmp(s, "arima ", 6)) {
	s += 6;
    } else if (!strncmp(s, "dpanel ", 7)) {
	s += 7;
    } else if (!strncmp(s, "var ", 4)) {
	s += 4;
    }

#if LAGS_DBG
    fprintf(stderr, "looking at '%s'\n", s);
#endif

    S = split_lag_fields(s, &ns, specials, cmd, &rem);
    if (S == NULL) {
	return cmd->err;
    }

#if LAGS_DBG
    int i;
    for (i=0; i<ns; i++) {
	fprintf(stderr, "S[%d] = '%s'\n", i, S[i]);
    }
#endif

    /* save original command line for echo, before modifying */
    free(cmd->parm2);
    cmd->parm2 = gretl_strdup(line);

    if (cmd->ci == ARMA) {
	handle_arma_lags(cmd, S, ns, specials, s, rem);
    } else {
	handle_single_lagvec(cmd, S[0], s, rem);
    }

    strings_array_free(S, ns);

#if LAGS_DBG
    fprintf(stderr, "revised line = '%s'\n", line);
#endif

    return cmd->err;
}

static char *got_gmm_spec (char *s)
{
    /* return whichever variant is found first */
    char *p1 = strstr(s, "GMM(");
    char *p2 = strstr(s, "GMMlevel(");

    if (p1 != NULL && p2 == NULL) {
	return p1;
    } else if (p2 != NULL && p1 == NULL) {
	return p2;
    } else if (p1 != NULL && p2 != NULL) {
	return (p2 - p1 > 0)? p1 : p2;
    } else {
	return NULL;
    }
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

    while ((s = got_gmm_spec(s)) != NULL) {
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

static int lag_from_lstr (const char *s, int *err)
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
	} else if (gretl_is_scalar(s)) {
	    lag = gretl_scalar_get_value(s, NULL);
	} else {
	    *err = 1;
	}
    }

    if (!*err) {
	lag = lsign * lag;
    }

    return lag;
}

static int get_contiguous_lags (LAGVAR *lv,
				const char *l1str, 
				const char *l2str)
{
    int err = 0;

    lv->lmin = lag_from_lstr(l1str, &err);

    if (!err) {
	lv->lmax = lag_from_lstr(l2str, &err);
    }

    return err;
}

static int parse_lagvar (const char *s, LAGVAR *lv, 
			 const DATASET *dset)
{
    char l1str[16], l2str[16];
    char fmt[32];
    int i, err = 1;

    lv->v = 0;
    lv->vlist = NULL;
    *lv->name = '\0';
    lv->lmin = 0;
    lv->lmax = 0;
    lv->laglist = NULL;

    sprintf(fmt, "%%%d[^(](%%%ds", VNAMELEN - 1, 15);

    if (sscanf(s, fmt, lv->name, l1str) != 2) {
	return err;
    }

#if LAG_DEBUG
    fprintf(stderr, "parse_lagvar: name = '%s'\n", lv->name);
#endif

    lv->v = current_series_index(dset, lv->name);
    if (lv->v <= 0) {
	lv->vlist = get_list_by_name(lv->name);
	if (lv->vlist == NULL) {
	    return err;
	} else {
	    lv->v = 0;
	}
    }

    sprintf(fmt, "%%%d[^(](%%%ds to %%%d[^)])", VNAMELEN - 1, 15, 15);

    if (sscanf(s, fmt, lv->name, l1str, l2str) == 3) {
	err = get_contiguous_lags(lv, l1str, l2str);
    } else if (strchr(l1str, ',') != NULL) {
	lv->laglist = gretl_list_from_string(strchr(s, '('), &err);
	if (lv->laglist != NULL) {
	    for (i=1; i<=lv->laglist[0]; i++) {
		lv->laglist[i] = -lv->laglist[i];
	    }
	    err = 0;
	}
    } else {
	sprintf(fmt, "%%%d[^(](%%%d[^ )]", VNAMELEN - 1, 15);
	sscanf(s, fmt, lv->name, l1str);
	lv->lmin = lv->lmax = lag_from_lstr(l1str, &err);
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

static int cmd_full_list (const DATASET *dset, CMD *cmd)
{
    int nv = 0, err = 0;
    int *list;

    if (cmd->ci == PRINT && cmd->parm2 != NULL &&
	*cmd->parm2 != '\0') {
	/* no-op */
	return 0;
    }

    if (cmd->flags & CMD_NULLIST) {
	/* no-op */
	cmd->flags ^= CMD_NULLIST;
	return 0;
    }

    list = full_var_list(dset, &nv);

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
   lags of a given series (or named list of series) to the 
   command list 
*/

int auto_lag_ok (const char *s, int *lpos,
		 DATASET *dset, CMD *cmd)
{
    LAGVAR lagvar;
    int i, j, nlags;
    int *vlist = NULL;
    int unilist[2] = {1, 0};
    int llen = *lpos;
    int lincr = 1;
    int ok = 1;
	
    if (parse_lagvar(s, &lagvar, dset)) {
	ok = 0;
	goto bailout;
    }

    if (lagvar.vlist != NULL) {
	/* we got a list of series to process */
	if (lagvar.vlist[0] <= 0) {
	    cmd->err = E_DATA;
	    ok = 0;
	    goto bailout;
	} 
	vlist = lagvar.vlist;
    } else {
	/* we got a single series */
	unilist[1] = lagvar.v;
	vlist = unilist;
    }

    /* number of lags per series */
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
    } else {
	int totlags = vlist[0] * nlags;

	if (totlags > 1 && expand_command_list(cmd, totlags)) {
	    ok = 0;
	}
    }

    for (j=1; j<=vlist[0] && ok; j++) {
	int vj = vlist[j];

	for (i=0; i<nlags && ok; i++) {
	    int order, lv;

	    if (lagvar.laglist != NULL) {
		order = lagvar.laglist[i+1];
	    } else {
		order = lagvar.lmin + i * lincr;
	    }

	    lv = laggenr(vj, order, dset);

#if LAG_DEBUG
	    fprintf(stderr, "laggenr for var %d (%s), lag %d, gave lv = %d\n",
		    vj, dset->varname[vj], order, lv);
#endif
	    if (lv < 0) {
		cmd->err = 1;
		gretl_errmsg_set(_("generation of lag variable failed"));
		ok = 0;
	    } else {
		/* Record info regarding the auto-generation of lags
		   so that we'll be able to echo the command properly --
		   see the echo_cmd() function.  Note: 'lagvar.v' is the
		   "source" variable; 'lv' is the ID number of the
		   generated lag.
		*/
		cmd->list[llen++] = lv;
		cmd->err = list_lag_info_add(vj, order, lv, llen - 1, cmd);
		if (cmd->err) {
		    ok = 0;
		}
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
				int *vnum, const DATASET *dset)
{
    int len = strcspn(s, ",;");
    int err = 0;

    if (len < strlen(s)) {
	char ostr[VNAMELEN] = {0};
	char word[32] = {0};
	char fmt[12];
	int v;

	sprintf(fmt, "%%%d[^ ,;]", VNAMELEN - 1);

	sscanf(s, fmt, ostr);
	if (isdigit(*ostr)) {
	    *order = atoi(ostr);
	} else if (gretl_is_scalar(ostr)) {
	    *order = gretl_scalar_get_value(ostr, NULL);
	} else {
	    ; /* FIXME error condition */
	}
	sscanf(s + len + 1, "%31[^ )]", word);
	v = series_index(dset, word);
	if (v < dset->v) {
	    *vnum = v;
	} else {
	    *lname = gretl_word_strdup(s + len + 1, NULL, 
				       OPT_NONE, &err);
	}
    } else {
	*lname = gretl_word_strdup(s, NULL, OPT_NONE, &err);
    }
}

static char *get_transform_param (const char *s, int *err)
{
    const char *p = s;
    int inparen = 1;

    while (*p) {
	if (*p == '(') {
	    inparen++;
	} else if (*p == ')') {
	    inparen--;
	}
	if (inparen == 0) {
	    break;
	}
	p++;
    }

    if (inparen == 0 && p - s > 0) {
	return gretl_strndup(s, p - s);
    } else {
	*err = E_PARSE;
	return NULL;
    }
}

static int auto_transform_ok (const char *s, int *lpos,
			      DATASET *dset, CMD *cmd)
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
				   dset);
	    } else {
		param = get_transform_param(s, &cmd->err);
	    }

	    if (param != NULL) {
		/* try for a named list */
		gotlist = get_list_by_name(param);
		if (gotlist != NULL) {
		    genlist = gretl_list_copy(gotlist);
		} else {
		    vnum = series_index(dset, param);
		    if (vnum == dset->v) {
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
	cmd->err = list_loggenr(genlist, dset);
    } else if (trans == DIFF || trans == LDIFF || trans == SDIFF) {
	cmd->err = list_diffgenr(genlist, trans, dset);
    } else if (trans == ORTHDEV) {
	cmd->err = list_orthdev(genlist, dset);
    } else if (trans == SQUARE) {
	cmd->err = list_xpxgenr(&genlist, dset, opt);
    } else if (trans == LAGS) {
	cmd->err = list_laggenr(&genlist, order, dset, OPT_NONE);
    } else if (trans == DUMMIFY) {
	cmd->err = list_dumgenr(&genlist, dset, OPT_F);
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
			DATASET *dset, CMD *cmd)
{
    int ok = 0;

    if (!strcmp(s, "time")) {
	if (cmd->ci == GNUPLOT) {
	    cmd->list[0] -= 1;
	    cmd->opt |= OPT_T;
	    ok = 1; /* handled */
	} else {
	    cmd->err = gen_time(dset, 1);
	    if (!cmd->err) {
		cmd->list[*lpos] = series_index(dset, "time");
		*lpos += 1;
		ok = 1; /* handled */
	    }
	}
    }

    return ok;
}

static int wildcard_expand (const char *s, int *lpos,
			    const DATASET *dset, CMD *cmd)
{
    int err = 0, ok = 0;

    if (strchr(s, '*') != NULL) {
	int *wildlist = varname_match_list(dset, s, &err);

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

static int print_name_ok (const char *s, CMD *cmd)
{
    int ok = 0;

    if (cmd->ci == PRINT) {
	GretlType t = user_var_get_type_by_name(s);

	if (t == GRETL_TYPE_MATRIX ||
	    t == GRETL_TYPE_DOUBLE ||
	    t == GRETL_TYPE_STRING ||
	    t == GRETL_TYPE_BUNDLE ||
	    !strcmp(s, "scalars")) {
	    cmd->parm2 = gretl_str_expand(&cmd->parm2, s, " ");
	    cmd->list[0] -= 1;
	    ok = 1;
	}
    }

    return ok;
}

static int delete_name_ok (const char *s, CMD *cmd)
{
    char bname[VNAMELEN];
    char fmt[10];
    int ok = 0;

    sprintf(fmt, "%%%d[^[.]", VNAMELEN - 1);

    if (sscanf(s, fmt, bname) == 1 &&
	gretl_is_bundle(bname)) {
	free(cmd->param);
	cmd->param = gretl_strdup(s);
	cmd->list[0] -= 1;
	ok = 1;
    }

    return ok;
}

static void parse_rename_cmd (const char *s, CMD *cmd, 
			      const DATASET *dset)
{
    int vtest, vtarg;
    char targ[VNAMELEN];
    char newname[VNAMELEN];
    char fmt[10], numstr[8];

    sprintf(fmt, "%%%ds %%%ds", VNAMELEN-1, VNAMELEN-1);

    if (sscanf(s, fmt, targ, newname) != 2) {
	cmd->err = E_PARSE;
	return;
    }

    if (isdigit(*targ)) {
	vtarg = atoi(targ);
	if (vtarg >= dset->v || vtarg < 1) {
	    cmd->err = E_DATA;
	    gretl_errmsg_sprintf(_("Variable number %d is out of bounds"), vtarg);
	    return;
	}
    } else {
	/* we're given the name of a variable? */
	vtarg = series_index(dset, targ);
	if (vtarg >= dset->v) {
	    cmd->err = E_UNKVAR;
	    return;
	}
    } 

    vtest = series_index(dset, newname);
    if (vtest == vtarg) {
	; /* no-op */
    } else if (vtest < dset->v) {
	gretl_errmsg_sprintf(_("A series named %s already exists"), newname);
	cmd->err = E_DATA;
	return;
    }

    if (vtest != vtarg) {
	if (check_varname(newname)) {
	    cmd->err = E_DATA;
	    return;
	}
	if (gretl_type_from_name(newname, dset)) {
	    cmd->err = E_TYPES;
	    return;
	}
    }

    /* write newname into cmd->param */
    free(cmd->param);
    cmd->param = gretl_strdup(newname);

    /* write target ID into cmd->parm2 */
    sprintf(numstr, "%d", vtarg);
    free(cmd->parm2);
    cmd->parm2 = gretl_strdup(numstr);

    if (cmd->param == NULL || cmd->parm2 == NULL) {
	cmd->err = E_ALLOC;
    }
}

static void parse_outfile_cmd (const char *s, CMD *cmd)
{
    int len = 0, quoted = 0;

    while (isspace(*s)) {
	s++;
    }

    if (*s) {
	cmd->err = filename_to_param(cmd, s, &len, &quoted);
    }
}

static int small_positive_int (const char *s)
{
    if (integer_string(s)) {
	int k = atoi(s);

	if (k > 0 && k <= 10) {
	    return 1;
	}
    }

    return 0;
}

static void handle_spreadsheet_params (const char *rem, CMD *cmd)
{
    int err = 0;

    if (cmd->opt & OPT_O) {
	/* odbc: spreadsheet-specific options not acceptable */
	err = incompatible_options(cmd->opt, OPT_O | OPT_C |
				   OPT_R | OPT_S);
    } else if (cmd->opt & OPT_W) {
	/* web database: ditto */
	err = incompatible_options(cmd->opt, OPT_W | OPT_C |
				   OPT_R | OPT_S);
    }

    if (!err) {
	err = incompatible_options(cmd->opt, OPT_O | OPT_W);
    }

    if (!err && (cmd->opt & (OPT_R | OPT_C | OPT_S))) {
	/* row offset, column offset, sheet name/number */
	int r0 = 0, c0 = 0;
	const char *s = NULL;

	if (cmd->opt & OPT_R) {
	    /* --rowoffset */
	    r0 = get_optval_int(cmd->ci, OPT_R, &err);
	}

	if (!err && (cmd->opt & OPT_C)) {
	    /* --coloffset */
	    c0 = get_optval_int(cmd->ci, OPT_C, &err);
	}

	if (!err && (cmd->opt & OPT_S)) {
	    /* --sheet */
	    s = get_optval_string(cmd->ci, OPT_S);
	    if (s == NULL) {
		err = E_DATA;
	    } 
	}

	if (!err) {
	    int slist[4] = {3, 0, c0, r0};

	    free(cmd->list);
	    cmd->list = gretl_list_copy(slist);
	    if (cmd->list == NULL) {
		err = E_ALLOC;
	    } else {
		if (small_positive_int(s)) {
		    /* take the --sheet spec as giving a sheet
		       number (1-based) */
		    cmd->list[1] = atoi(s);
		} else if (s != NULL) {
		    /* take it as giving a sheet name */
		    free(cmd->parm2);
		    cmd->parm2 = gretl_strdup(s);
		    if (cmd->parm2 == NULL) {
			err = E_ALLOC;
		    }
		}
	    }
	}
    }

    cmd->err = err;
}

#define FIELDLEN 512

static int get_field_length (const char *s)
{
    const char *p = s;
    int inparen = 0;
    int len = 0;

    while (*p) {
	if (*p == '(') {
	    inparen++;
	} else if (*p == ')') {
	    inparen--;
	}
	if (!inparen && *p == ' ') {
	    break;
	}
	p++;
	len++;
    }

    if (len >= FIELDLEN) {
	fprintf(stderr, "list field in command is too long "
		"(len = %d, max = %d)\n", len, FIELDLEN);
	fprintf(stderr, "s = '%s'\n", s);
	gretl_errmsg_set("Overflow in command list field");
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
	err = E_DATA;
    }

#if CMD_DEBUG
    fprintf(stderr, "get_next_field: got '%s'\n", field);
#endif

    return err;
}

/* look for a line with an "implicit genr", such as
   y = 3*x, x += 10, etc. */

int plausible_genr_start (const char *s, const DATASET *dset)
{
    int ret = 0;

    if (strchr(s, '=') || strstr(s, "++") || strstr(s, "--")) {
	const char *ok = ".+-*/%^~|=[";
	char word[VNAMELEN] = {0};
	char fmt[20];

	sprintf(fmt, "%%%d[^[ .+*/%%^~|=-]", VNAMELEN - 1);

	if (sscanf(s, fmt, word)) {
	    s += strlen(word);
	    while (*s == ' ') s++;
	    if (strspn(s, ok) > 0 && check_varname(word) == 0) {
		ret = 1;
	    }
	}
    } else if (gretl_type_from_name(s, dset) != 0) {
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
    } else if (cmd->aux != DS_SORTBY && cmd->aux != DS_DSORTBY) {
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

#define COMMAND_CAN_END(c) (c == FOREIGN ||	\
			    c == FUNC ||	\
                            c == GMM ||		\
                            c == KALMAN ||	\
                            c == MLE ||		\
			    c == MPI ||		\
                            c == NLS ||		\
			    c == RESTRICT ||	\
			    c == SYSTEM)

static int check_end_command (CMD *cmd)
{
    if (cmd->param != NULL && *cmd->param != 0) {
	int cmdcode = gretl_command_number(cmd->param);

	if (cmdcode == LOOP) {
	    cmd->ci = ENDLOOP;
	} else if (!COMMAND_CAN_END(cmdcode)) {
	    cmd->err = 1;
	    gretl_errmsg_sprintf(_("command 'end %s' not recognized"), 
				 cmd->param);
	}
    } else {
	cmd->err = 1;
	gretl_errmsg_set(_("end: nothing to end")); 
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

static void cmd_param_grab_word (CMD *cmd, const char *s)
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

static void param_grab_braced (CMD *cmd, const char *s)
{
    if (*s == '{') {
	const char *p = strchr(s, '}');

	if (p == NULL) {
	    cmd->err = E_PARSE;
	} else {
	    int n = p - s + 1;

	    free(cmd->param);
	    cmd->param = gretl_strndup(s, n);
	    if (cmd->param == NULL) {
		cmd->err = E_ALLOC;
	    } 
	}	    
    } else {
	cmd_param_grab_word(cmd, s);
    }
}

static void param_grab_quoted (CMD *cmd, const char *s)
{
    if (*s == '"') {
	const char *p = strchr(s+1, '"');

	if (p == NULL) {
	    cmd->err = E_PARSE;
	} else {
	    int n = p - s - 1;

	    free(cmd->param);
	    cmd->param = gretl_strndup(s+1, n);
	    if (cmd->param == NULL) {
		cmd->err = E_ALLOC;
	    } 
	}	    
    } else {
	cmd_param_grab_word(cmd, s);
    }
}

/* Capture the next 'word' found following the initial command word
   (or the whole remainder of the line in some cases) as the parameter
   for @cmd.  Flag an error if the command requires a parameter but
   none is found.
*/

static int capture_param (CMD *cmd, const char *s)
{
    /* if param has already been written by some special
       routine, don't overwrite it */
    if (*cmd->param != '\0') {
	if (cmd->ci == DATAMOD) {
	    check_datamod_command(cmd, s);
	}
	return cmd->err;
    }

    s += strspn(s, " ");

    if (string_is_blank(s)) {
	if (REQUIRES_PARAM(cmd->ci) || REQUIRES_ORDER(cmd->ci)) {
	    cmd->err = E_PARSE;
	    gretl_errmsg_sprintf(_("%s: required parameter is missing"),
				 cmd->word);
	}
    } else {
	if (cmd->ci == PRINT || cmd->ci == FUNCERR || 
	    cmd->ci == DELEET || cmd->ci == HELP ||
	    cmd->ci == EQUATION) {
	    /* grab the whole remainder of line */
	    cmd_param_grab_string(cmd, s);
	} else if (cmd->ci == QUANTREG || cmd->ci == LEVINLIN) {
	    param_grab_braced(cmd, s);
	} else if (cmd->ci == OPEN || cmd->ci == APPEND ||
		   cmd->ci == JOIN) {
	    param_grab_quoted(cmd, s);
	} else {
	    /* grab one 'word' */
	    cmd_param_grab_word(cmd, s);
	}
#if CMD_DEBUG
	fprintf(stderr, "capture_param: s='%s', param='%s'\n",
		s, cmd->param);
#endif
	if (REQUIRES_ORDER(cmd->ci) && cmd->ci != LEVINLIN) {
	    cmd->order = gretl_int_from_string(cmd->param, &cmd->err);
	    if (cmd->err) {
		gretl_errmsg_sprintf(_("%s: expected an integer order"),
				     cmd->word);
		cmd->err = E_PARSE;
	    }
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

    cmd_unset_nolist(cmd);

    if (cmd->list == NULL || cmd->param == NULL || cmd->parm2 == NULL) {
	cmd->err = E_ALLOC;
    } else {
	cmd->list[0] = 0;
	*cmd->param = '\0';
	*cmd->parm2 = '\0';
    }

    free(cmd->auxlist);
    cmd->auxlist = NULL;

    cmd_lag_info_destroy(cmd);
    clear_option_params();

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
    const char *p = s;
    int inparen = 0;
    int nf = 0;

#if CMD_DEBUG
    fprintf(stderr, "count_free_fields: looking at '%s'\n", s);
#endif

    if (s != NULL) {
	while (*s) {
	    if (!inparen && *s != ' ' && (s == p || *(s-1) == ' ')) {
		/* non-space preceded by space, or at start */
		nf++;
	    }
	    if (*s == '(') {
		inparen++;
	    } else if (*s == ')') {
		inparen--;
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

#define semi_special(c) (c == ARBOND || c == DPANEL)

static int handle_semicolon (int *k, int *ints_ok, int *poly, 
			     int *sepcount, CMD *cmd)
{
    int ok = 0;

    if (USES_LISTSEP(cmd->ci)) {
	cmd->list[*k] = LISTSEP;
	*k += 1;
	*sepcount -= 1;
	if (*ints_ok) {
	    if (*sepcount == 0 || (*sepcount == 1 && semi_special(cmd->ci))) {
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
			  const DATASET *dset, CMD *cmd)
{
    char *test;
    int v, ok = 0;

    errno = 0;

    v = strtol(s, &test, 10);
    if (*test != '\0' || errno == ERANGE) {
	return 0;
    } 

    if (!ints_ok && !poly && v >= dset->v) {
	cmd->err = E_UNKVAR;
	gretl_errmsg_sprintf(_("%d is not a valid variable number"), v);
    } else {
	cmd->list[*k] = v;
	*k += 1;
	ok = 1;
    }

    return ok;
}

static int parse_alpha_list_field (const char *s, int *pk, int ints_ok,
				   DATASET *dset, CMD *cmd)
{
    int *xlist;
    int v, k = *pk;
    int ok = 0;

#if CMD_DEBUG
    fprintf(stderr, "parse_alpha_list_field: s = '%s', ci = %d (%s)\n", 
	    s, cmd->ci, cmd->word);
#endif

    if (ints_ok) {
	v = gretl_int_from_string(s, &cmd->err);
	if (!cmd->err) {
	    cmd->list[k++] = v;
	    ok = 1;
	}
    } else if ((v = series_index(dset, s)) < dset->v) {
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
	if (auto_lag_ok(s, &k, dset, cmd)) {
	    /* lag specification, e.g. 'var(-1)' */
	    ok = 1;
	} else if (auto_transform_ok(s, &k, dset, cmd)) {
	    /* automated transformations such as 'logs(list)' */
	    ok = 1;	
	}
    } else if (add_time_ok(s, &k, dset, cmd)) {
	ok = 1;	
    } else if (wildcard_expand(s, &k, dset, cmd)) {
	ok = 1;
    } else if (cmd->ci == PRINT && print_name_ok(s, cmd)) {
	ok = 1;
    } else if (cmd->ci == DELEET && delete_name_ok(s, cmd)) {
	ok = 1;
    }

    *pk = k;

    if (!ok && cmd->err == 0) {
	if (user_var_get_type_by_name(s)) {
	    gretl_errmsg_sprintf(_("'%s' is not the name of a series"), s);
	    cmd->err = E_DATATYPE;
	} else {
	    gretl_errmsg_sprintf(_("'%s' is not the name of a variable"), s);
	    cmd->err = E_UNKVAR;
	}
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

static int ends_block (const char *s, const char *blocktype)
{
    if (!strncmp(s, "end ", 4)) {
	s += 3;
	s += strspn(s, " \t");
	if (!strncmp(s, blocktype, strlen(blocktype))) {
	    return 1;
	}
    }

    return 0;
}

/* Get the first word out of line.  In general this should be a
   command word (starting with a alphabetical character), but there
   are a few special case: shell commands start with the "!" escape;
   restriction specifications may start with "-" (as in "-b1 + b2 =
   0") or a numerical multiplier. 

   In addition the beginning of the line may take the form of
   a save to a named object as in "foo <- command".
*/

static int get_first_word (const char *line, char *cnext, CMD *cmd)
{
    int n, ret = 0;

    *cmd->word = '\0';

    if (!cmd->context && gretl_function_depth() == 0) {
	if (strstr(line, " <- ") != NULL) {
	    line = maybe_skip_savename(line);
	}
    } 

    n = gretl_namechar_spn(line);

    if (cmd->context == RESTRICT && n == 0) {
	/* non-alpha may be OK */
	ret = !string_is_blank(line);
    } else if (*line == '!') {
	/* shell escape */
	strcpy(cmd->word, "!");
	ret = 1;
    } else if (n > 0) {
	/* got some alphabetical stuff */
	if (n > FN_NAMELEN - 1) {
	    n = FN_NAMELEN - 1;
	}
	strncat(cmd->word, line, n);
	*cnext = line[n];
	ret = 1;
    } else if (!string_is_blank(line)) {
	/* must be garbage? */
	cmd->err = E_PARSE;
    }

    return ret;
}

/* For commands that need a list: @line is the portion of the command
   line we're reading from and @nf is the number of fields to be
   processed.  
*/

static int process_command_list (CMD *cmd, const char *line, int nf,
				 DATASET *dset)
{
    char field[FIELDLEN] = {0};
    int poly = 0, ints_ok = 0;
    int sepcount;
    int j, k;

    /* get a count of ';' separators in line */
    sepcount = get_sepcount(line);
    cmd->err = sepcount_error(cmd->ci, sepcount);
    if (cmd->err) {
	return cmd->err;
    }

    /* allocate space for the command list */
    if (resize_command_list(cmd, nf)) {
	return cmd->err;
    }    

    if (cmd->ci == AR || cmd->ci == ARBOND || cmd->ci == DPANEL ||
	cmd->ci == ARMA || cmd->ci == GARCH) {
	/* flag acceptance of plain ints in list */
	ints_ok = 1;
    } else if (matrix_data_option(cmd->ci, cmd->opt)) {
	/* the list refers to columns of a matrix */
	ints_ok = 1;
    }

    for (j=1, k=1; j<=nf; j++) {
	int ok = 0;

	/* special: optional width for correlogram, periodogram */
	if ((cmd->ci == CORRGM || cmd->ci == PERGM ||
	     cmd->ci == FRACTINT) && j == 2) {
	    cmd->list[0] = 1;
	    cmd_param_grab_word(cmd, line);
	    break;
	} else if (cmd->ci == XCORRGM && j == 3) {
	    cmd->list[0] = 2;
	    cmd_param_grab_word(cmd, line);
	    break;
	}	    

	cmd->err = get_next_field(field, line);
	if (cmd->err) {
	    break;
	}

	if (isalpha((unsigned char) *field)) {
	    ok = parse_alpha_list_field(field, &k, ints_ok, dset, cmd);
	} else if (*field == '*') {
	    ok = wildcard_expand(field, &k, dset, cmd);
	} else if (isdigit(*field)) {
	    ok = get_id_or_int(field, &k, ints_ok, poly, dset, cmd);
	} else if (*field == ';') {
	    ok = handle_semicolon(&k, &ints_ok, &poly, &sepcount, cmd);
	} 

	if (!ok) {
	    if (cmd->err == 0) {
		cmd->err = 1;
	    } 
	    if (!gretl_errmsg_is_set()) {
		if (*field == '=' && cmd->ci != GENR) {
		    gretl_errmsg_sprintf(_("'%s' may not be used as a "
					   "variable name"), cmd->word);
		} else {
		    gretl_errmsg_sprintf(_("field '%s' in command is invalid"), 
					 field);
		}
	    }
	    break;
	}

	/* advance for next read */
	line += strlen(field) + 1;
    }

    if (cmd->err) {
	return cmd->err;
    }

    /* commands that can take a specified list, but where if the
       list is null or just ";" we want to operate on all variables
    */    
    if (DEFAULTS_TO_FULL_LIST(cmd->ci)) {
	if (cmd->list[0] == 0) {
	    if (cmd->ci == SMPL) {
		/* "smpl" accepts an empty list as "all vars", with
		   the --no-missing or --contiguous options, so
		   leave well alone */
		;
	    } else if (cmd->ci == SUMMARY && (cmd->opt & OPT_X)) {
		/* summary with --matrix option: leave alone */
		;
	    } else {
		cmd_full_list(dset, cmd);
	    }
	    /* suppress echo of the list -- may be too long */
	    cmd_set_nolist(cmd);
	}
    } else if (cmd->ci != SETMISS && 
	       cmd->ci != PRINT &&
	       cmd->ci != GNUPLOT &&
	       cmd->ci != SCATTERS &&
	       cmd->ci != BXPLOT &&
	       cmd->ci != DELEET) {
	/* the command needs a list but doesn't have one */
	if (cmd->list[0] == 0) {
	    cmd->err = E_ARGS;
	}
    }

    if (NEEDS_TWO_VARS(cmd->ci) && cmd->list[0] == 1) {
	cmd->err = E_ARGS;
    }

    if (!cmd->err && cmd->ci == GNUPLOT && cmd->list[0] < 2) {
	/* check the list for the gnuplot command */
	if (cmd->opt & (OPT_D | OPT_X)) {
	    ; /* OK: non-empty list not required */
	} else if (cmd->opt & OPT_T) {
	    /* time-series otion: only one series needed */
	    if (cmd->list[0] < 1) {
		cmd->err = E_ARGS;
	    }
	} else {
	    /* all other cases: we need at least two series */
	    cmd->err = E_ARGS;
	}
    }

    if (!cmd->err && cmd->ci == SCATTERS && cmd->list[0] < 2) {
	/* check the list for the scatters command */
	if ((cmd->opt & OPT_X) && (cmd->opt & OPT_T)) {
	    ; /* non-empty list not required */
	} else {
	    /* all other cases: we need at least two series */
	    cmd->err = E_ARGS;
	}
    }

    if (!cmd->err && cmd->ci == BXPLOT && cmd->list[0] < 1) {
	/* check the list for the boxplot command */
	if (cmd->opt & OPT_X) {
	    ; /* matrix: non-empty list not required */
	} else {
	    /* all other cases: we need at least one series */
	    cmd->err = E_ARGS;
	}
    }    

    /* check list for duplicated variables? */
    if (!cmd->err && !cmd_nolist(cmd)) {
	int dupv = gretl_list_duplicates(cmd->list, cmd->ci);

	if (dupv >= 0) {
	    printlist(cmd->list, "cmd->list with duplicate(s)");
	    cmd->err = E_DATA;
	    gretl_errmsg_sprintf(_("variable %d duplicated in the "
				   "command list."), dupv);
	} 
    }

    return cmd->err;
}

/**
 * parse_command_line:
 * @line: the command line.
 * @cmd: pointer to command struct.
 * @dset: dataset struct.
 * @ptr: pointer for use with "compilation" of
 * conditionals in loops.
 *
 * Parses @line and fills out @cmd accordingly. 
 *
 * Returns: 0 on success, non-zero code on error.
 */

int parse_command_line (char *line, CMD *cmd, DATASET *dset, void *ptr) 
{
    char *rem = NULL;
    char cnext = 0;
    int nf;

    if (gretl_cmd_clear(cmd)) {
	return cmd->err;
    }

    gretl_error_clear();

#if CMD_DEBUG
    fprintf(stderr, "parse_command_line: '%s'\n", line);
#endif

    if (!cmd_nosub(cmd)) {
	int subst = 0;

	cmd->err = substitute_named_strings(line, &subst);
	if (cmd->err) {
	    return cmd->err;
	} else if (subst) {
	    /* record the fact that substitution has been done */
	    cmd->flags |= CMD_SUBST;
	} else {
	    cmd->flags &= ~CMD_SUBST;
	}
    }

#if CMD_DEBUG
    if (cmd->flags & CMD_SUBST) {
	fprintf(stderr, "after substitution: '%s'\n", line);
    }
#endif

    if (cmd->context == FOREIGN && 
	!ends_block(line, "foreign") &&
	!ends_block(line, "mpi")) {
	cmd_set_nolist(cmd);
	cmd->opt = OPT_NONE;
	cmd->ci = FOREIGN;
	return 0;
    }

    if ((cmd->flags & CMD_SUBST) || !gretl_looping_currently()) {
	/* normalize line spaces */
	compress_spaces(line);
	
	/* trap lines that are nothing but comments */
	if (filter_comments(line, cmd)) {
	    return 0;
	}

	/* catch errors associated with comment syntax */
	if (cmd->err) {
	    return cmd->err;
	}
    }

#if 0
    test_tokenize(line, cmd, dset, ptr);
#endif

    /* check for "catch" */
    maybe_set_catch_flag(line, cmd);

    if (!cmd_noopt(cmd)) {
	/* extract any options */
	cmd->opt = get_gretl_options(line, &cmd->err);
	if (cmd->err) {
	    return cmd->err;
	}
    }

    if (!cmd->context) {
	if (gretl_function_depth() == 0 && strstr(line, " <- ") != NULL) {
	    /* extract "savename" for storing an object? */
	    maybe_extract_savename(line, cmd);
	} else {
	    *cmd->savename = '\0';
	}
    } 

    /* maybe there's no command here? */
    if (!get_first_word(line, &cnext, cmd)) {
	cmd_set_nolist(cmd);
	cmd->ci = CMD_NULL;
	return cmd->err;
    }

    if (!cmd->context) {
	/* replace simple aliases and a few specials */
	catch_command_alias(line, cmd);
    } else if (cmd->context == SYSTEM) {
	catch_system_alias(cmd);
    }

    /* subsetted commands (e.g. "deriv" in relation to "nls") */
    if (!strcmp(cmd->word, "end")) {
	cmd->context = 0;
    } else if (cmd->context && cmd->ci != EQUATION) {
	/* "equation" occurs in the SYSTEM context, but it is
	   a command in its own right, so don't overwrite 
	   cmd->ci with cmd->context
	*/
	cmd->ci = cmd->context;
    }

    if (cmd->ci == 0) {
	if (!strcmp(cmd->word, "elif") && cnext == '(') {
	    /* bodge: temporary reprieve for SVAR */
	    cmd->ci = ELIF;
	} else if (cnext != '(') {
	    /* regular command, not a function call */
	    cmd->ci = gretl_command_number(cmd->word);
	}
	if (cmd->ci == 0) {
	    if (plausible_genr_start(line, dset)) {
		cmd->ci = GENR;
	    } else if (function_lookup(cmd->word)) {
		cmd->ci = GENR;
		cmd->opt = OPT_O;
	    } else if (get_user_function_by_name(cmd->word)) {
		cmd->ci = GENR;
		cmd->opt = OPT_O;
	    } else if (gretl_if_state_false()) {
		cmd_set_nolist(cmd);
		cmd->ci = CMD_MASKED;
		return 0;
	    } else {
		cmd->err = 1;
		if (gretl_command_number(cmd->word) == 0) {
		    gretl_errmsg_sprintf(_("command '%s' not recognized"), 
					 cmd->word);
		}
		goto cmd_exit;
	    }
	}
    }

#if CMD_DEBUG
    fprintf(stderr, "cmd->ci = %d\n", cmd->ci);
#endif

    /* if, else, endif controls: should this come earlier? */
    if (flow_control(line, dset, cmd, ptr)) {
	cmd_set_nolist(cmd);
	cmd->ci = CMD_MASKED;
	return cmd->err;
    }

    /* special: list <listname> delete */
    if (cmd->ci == DELEET && *cmd->parm2 != '\0') {
	cmd_set_nolist(cmd);
	return cmd->err;
    }

    /* advance beyond the first 'word' on the line (contiguous
       non-space characters), skip a following space if there is
       one, and record our read position as 'rem'
    */
    rem = line + strcspn(line, " ");
    if (*rem != '\0') {
	rem++;
    }

    if (cmd->ci == EQNPRINT || cmd->ci == TABPRINT) {
	/* TeX printing commands can take a filename parameter,
	   but that's all
	*/
	get_optional_filename(rem, cmd);
	return cmd->err;
    } else if (cmd->ci == OUTFILE) {
	/* the "outfile" command may have a filename */
	parse_outfile_cmd(rem, cmd);
    } else if (cmd->ci == RENAME) {
	/* the "rename" command calls for a variable number and a
	   new name */
	parse_rename_cmd(rem, cmd, dset);
    } else if (cmd->ci == OPEN || cmd->ci == APPEND) {
	/* "open" and "append" may have spreadsheet parameters */
	handle_spreadsheet_params(rem, cmd);
	if (cmd->err) {
	    return cmd->err;
	}
    } 

    /* commands that never take a list of variables */
    if (NO_VARLIST(cmd->ci) || 
	(cmd->ci == DELEET && (cmd->opt & (OPT_D | OPT_T))) ||
	(cmd->ci == EQUATION && (cmd->opt & OPT_M))) { 
	cmd_set_nolist(cmd);
	if (cmd->ci != GENR) {
	    capture_param(cmd, rem);
	}
	return cmd->err;
    } 

    /* now for some commands which may or may not take a list:
       we can return early in the no-list cases 
    */

    if (cmd->ci == PRINT) {
	/* no list in string literal variant */
	if (strstr(line, "\"")) {
	    cmd_set_nolist(cmd);
	    capture_param(cmd, rem);
	    return cmd->err;
	} else if (cmd->flags & CMD_PROG) {
	    /* print in progressive loop */
	    free(cmd->parm2);
	    cmd->parm2 = gretl_strdup(rem);
	    if (cmd->parm2 == NULL) {
		cmd->err = E_ALLOC;
	    }
	    return cmd->err;
	}
    } else if (cmd->ci == SMPL) {
	/* SMPL may take a list, but only in case of OPT_M,
	   "--no-missing", or OPT_C, "--contiguous" */
	if (!(cmd->opt & (OPT_M | OPT_C))) {
	    cmd_set_nolist(cmd);
	    return cmd->err;
	}
    } else if (cmd->ci == BXPLOT) {
	/* boxplots take a list, but if there are Boolean conditions
	   embedded, the line has to be parsed specially */
	if (boxplot_booleans_present(line)) {
	    cmd_set_nolist(cmd);
	    return cmd->err;
	}
    } else if (cmd->ci == OMIT) {
	/* OMIT typically takes a list, but can be given without args
	   to omit the last variable */
	if (string_is_blank(line + 4)) {
	    cmd_set_nolist(cmd);
	    return cmd->err;
	} 
    } else if (cmd->ci == XTAB) {
	/* XTAB generally takes a list, but not with the --matrix option */
	if (cmd->opt & OPT_M) {
	    cmd_set_nolist(cmd);
	    return cmd->err;
	} 
    } else if (cmd->ci == DATAMOD) {
	/* dataset-modifying commands */
	capture_param(cmd, rem);
	if (cmd->aux != DS_SORTBY && 
	    cmd->aux != DS_DSORTBY) {
	    cmd_set_nolist(cmd);
	    return cmd->err;
	}
    } else if (cmd->ci == STORE && (cmd->flags & CMD_PROG)) {
	/* store in progressive loop */
	cmd->err = get_maybe_quoted_filename(cmd, &rem);
	if (!cmd->err) {
	    free(cmd->parm2);
	    cmd->parm2 = gretl_strdup(rem);
	    if (cmd->parm2 == NULL) {
		cmd->err = E_ALLOC;
	    }
	}
	return cmd->err;
    }	

    /* OK, now we're definitely doing a list-oriented command;
       we begin by taking care of a few specials 
    */

    if (cmd->ci == GNUPLOT || cmd->ci == BXPLOT) {
	/* we may have a block of stuff to pass literally
	   to gnuplot */
	grab_gnuplot_literal_block(rem, cmd);
    } else if (cmd->ci == ARMA || cmd->ci == DPANEL || cmd->ci == VAR) {
	/* allow for specific "gappy" lags */
	maybe_rewrite_lags(line, cmd);
    } 

    /* fix lines that contain a semicolon stuck to another element */
    fix_semicolon_separation(rem, cmd);
    if (cmd->err) {
	return cmd->err;
    }

    /* arbond special: if there's a block-diagonal instruments
       portion to the command, grab that in literal form for
       later processing. Note that this modifies @line, cutting
       out the special GMM() bits and storing them in cmd->param.
    */
    if ((cmd->ci == ARBOND || cmd->ci == DPANEL) && get_sepcount(line) == 2) {
	grab_arbond_diag(line, cmd);
	if (cmd->err) {
	    return cmd->err;
	}
    } 

    /* find the number of space-separated fields remaining
       in the command line
    */
    nf = count_free_fields(rem);

#if CMD_DEBUG
    fprintf(stderr, "nf=%d, remainder='%s'\n", nf, rem);
#endif

    if (cmd->ci == DELEET && nf == 1) {
	GretlType t = user_var_get_type_by_name(rem);

	if (t == GRETL_TYPE_DOUBLE ||
	    t == GRETL_TYPE_MATRIX ||
	    t == GRETL_TYPE_BUNDLE ||
	    t == GRETL_TYPE_STRING ||
	    !strcmp(rem, "kalman")) {
	    /* special for deleting a named matrix, string, ... */
	    cmd_param_grab_string(cmd, rem);
	    goto cmd_exit;
	}
    }

    /* specials where there's something that goes into "param",
       before the first semicolon */
    if (cmd->ci == LAGS) {
	if (get_lags_param(cmd, &rem)) {
	    nf = count_fields(rem, NULL);
	} 
    }    

    /* "store" is a special case since the filename that comes
       before the list may be quoted, and have spaces in it */
    if (cmd->ci == STORE && nf > 0) {
	cmd->err = get_maybe_quoted_filename(cmd, &rem);
	if (cmd->err) {
	    goto cmd_exit;
	} else {
	    nf = count_free_fields(rem);
	}
    }

    /* "setmiss" takes a value to be interpreted as missing;
       this are captured in cmd->param, as is the 'order' for
       a command that needs same
    */
    if (REQUIRES_ORDER(cmd->ci) || cmd->ci == SETMISS) {
	capture_param(cmd, rem);
	if (cmd->err) {
	    goto cmd_exit;
	} else {
	    rem += strlen(cmd->param) + 1;
	    nf--;
	} 
    }

    if (cmd->ci == QUANTREG) {
	/* quantreg requires a tau specification */
	capture_param(cmd, rem);
	if (cmd->err) {
	    goto cmd_exit;
	} else {
	    rem += strlen(cmd->param) + 1;
	    nf = count_free_fields(rem);
	} 
    } else if (cmd->ci == DATAMOD) {
	/* at this point, must be doing a dataset operation that
	   requires a list argument (e.g. sorting) */
	rem += strspn(rem, " ");
	rem += strcspn(rem, " ");
	nf--;
    } else if (cmd->ci == VECM) { 
	free(cmd->parm2);
	cmd->parm2 = gretl_word_strdup(rem, NULL, OPT_NONE, &cmd->err);
	rem += strlen(cmd->parm2) + 1;
	nf--;
    }

    /* By now we're looking at a command that takes a list,
       which either has been specified already or needs to
       be filled out automatically */

    cmd->err = process_command_list(cmd, rem, nf, dset);

 cmd_exit:

    /* double-check that allocation hasn't failed */
    if (cmd->err == 0 && (cmd->list == NULL || cmd->param == NULL || 
			  cmd->parm2 == NULL)) {
	cmd->err = E_ALLOC;
    }

#if CMD_DEBUG
    printlist(cmd->list, "cmd->list");
    fprintf(stderr, "cmd->err = %d, context=%d\n", cmd->err,
	    cmd->context);
#endif

    if (cmd->err) {
	cmd->context = 0;
    }

#if CMD_DEBUG
    fprintf(stderr, "parsed: '%s'\n", line);
#endif

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

static void output_help_line (const char *line, PRN *prn, int recode)
{
    if (recode > 0) {
	recode_print_line(line, prn);
    } else {
	pputs(prn, line);
    }
}

/* check in the CLI helpfile for a line of the form "  @s: ...", 
   where @s has been recognized as a libset variable */

static int got_setvar_line (const char *s, int n, char *line)
{
    if (!strncmp(line, "  ", 2) &&
	!strncmp(line + 2, s, n) &&
	line[2 + n] == ':') {
	return 1;
    } else {
	return 0;
    }
}

#define HELPLEN 128

/* special case: the user has done "help set @setvar" */	

static int do_set_help (const char *setvar, FILE *fp, 
			char *line, PRN *prn, 
			int recode)
{
    char s[9];
    int n = strlen(setvar);
    int count = 0;

    while (count < 2 && fgets(line, HELPLEN, fp) != NULL) {
	if (*line != '#') {
	    continue;
	}
	sscanf(line + 2, "%8s", s);
	if (!strcmp(s, "set")) {
	    while (count < 2 && fgets(line, HELPLEN, fp)) {
		if (got_setvar_line(setvar, n, line)) {
		    pputc(prn, '\n');
		    output_help_line(line + 2, prn, recode);
		    count++;
		} else if (count > 0) {
		    if (string_is_blank(line)) {
			/* reached the end of the item */
			pputc(prn, '\n');
			count++;
		    } else {
			output_help_line(line + 2, prn, recode);
		    }
		} else if (*line == '#') {
		    /* overshot: not found */
		    count = 999;
		}
	    }
	}
    }

    return (count > 0 && count < 999);
}

/* Is @s "set ", and if so, is @word the name of a variable that can be
   set via the 'set' command?  We try looking it up in libset.c.
   FIXME: there are some "special case" set variables that are not
   found via the function is_libset_var() at present.
*/

static int is_set_item_help (char *s, char *word, PRN *prn)
{
    if (!strncmp(s, "set ", 4)) {
	if (sscanf(s + 4, "%31s", word) == 1) {
	    if (!strcmp(word, "stopwatch")) {
		strcpy(s, "$stopwatch");
		*word = '\0';
	    } else if (is_libset_var(word)) {
		s[3] = '\0';
		return 1;
	    } else {
		pprintf(prn, "'%s' is not a settable variable\n", word);
		return -1;
	    }
	}
    }

    return 0;
}

static int is_help_alias (char *s)
{
    int ret = 0;

    if (!strcmp(s, "addobs")) {
	strcpy(s, "dataset");
	ret = 1;
    }

    return ret;
}

/**
 * cli_help:
 * @cmdword: the command on which help is wanted.
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

int cli_help (const char *cmdword, gretlopt opt, PRN *prn)
{
    static int recode = -1;
    char helpfile[FILENAME_MAX];
    FILE *fp;
    int noword, funhelp = (opt & OPT_F);
    char word[12], needle[32]; 
    char setvar[32], line[HELPLEN];
    int i, j, ok = 0;

    noword = (cmdword == NULL || *cmdword == '\0');

    *needle = *setvar = '\0';

    if (!noword) {
	strncat(needle, cmdword, 31);
    }

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

    if ((noword && funhelp) || !strcmp(needle, "functions")) {
	sprintf(helpfile, "%s%s", gretl_home(), _("genrcli.hlp"));
	return func_help_topics(helpfile, prn);
    }

    if (!funhelp) {
	ok = gretl_command_number(needle) > 0;
	if (!ok) {
	    ok = is_help_alias(needle);
	}
	if (!ok) {
	    ok = is_set_item_help(needle, setvar, prn);
	    if (ok < 0) {
		/* unrecognized "help set foo" */
		return 1;
	    } 
	}
    } 

    if (ok) {
	strcpy(helpfile, helpfile_path(GRETL_CLI_HELPFILE));
    } else if (genr_function_word(needle)) {
	sprintf(helpfile, "%sgenrcli.hlp", gretl_home());
    } else if (gretl_is_public_user_function(needle)) {
	return user_function_help(needle, OPT_NONE, prn);
    } else {
	pprintf(prn, _("\"%s\" is not a gretl command.\n"), needle);
	return 1;
    }

    if ((fp = gretl_fopen(helpfile, "r")) == NULL) {
	printf(_("Unable to access the file %s.\n"), helpfile);
	return 1;
    } 

    if (!gretl_in_gui_mode() && recode < 0) {
	recode = maybe_need_recode();
    }

    if (*setvar != '\0') {
	ok = do_set_help(setvar, fp, line, prn, recode);
	if (!ok) {
	    pprintf(prn, _("%s: sorry, no help available.\n"), cmdword);
	}
	fclose(fp);
	return 0;
    }

    ok = 0;
    while (fgets(line, sizeof line, fp) != NULL) {
	if (*line != '#') {
	    continue;
	}
	sscanf(line + 2, "%10s", word);
	if (!strcmp(needle, word)) {
	    ok = 1;
	    pprintf(prn, "\n%s\n", word);
	    while (fgets(line, sizeof line, fp)) {
		if (*line == '#') {
		    break;
		}
		output_help_line(line, prn, recode);
	    }
	    break;
	}
    }

    if (!ok) {
	pprintf(prn, _("%s: sorry, no help available.\n"), needle);
    }

    fclose(fp);

    return 0;
}

/**
 * parseopt:
 * @pargc: pointer to count of arguments.
 * @pargv: pointer to command-line argument array.
 * @popt: location to receive option(s).
 * @scriptval: location to receive numerical option value
 * to be passed to script.
 * @fname: optional filename argument.
 *
 * Parses options out of the command line into @popt and
 * fills out @fname if applicable.
 *
 * Returns: 0 on success, non-zero in case of bad options.
 */

int parseopt (int *pargc, char ***pargv, gretlopt *popt, 
	      double *scriptval, char *fname)
{
    char **argv;
    int argc, gotfile = 0;
    gretlopt opt = OPT_NONE;
    int err = 0;

    *fname = '\0';

    if (pargv == NULL) {
	return 0;
    }

    argc = *pargc;
    argv = *pargv;

    while (*++argv) {
	const char *s = *argv;

	if (!strcmp(s, "-e") || !strncmp(s, "--english", 9)) { 
	    opt |= OPT_ENGLISH;
	} else if (!strcmp(s, "-b") || !strncmp(s, "--batch", 7)) {
	    opt |= OPT_BATCH;
	} else if (!strcmp(s, "-h") || !strcmp(s, "--help")) { 
	    opt |= OPT_HELP;
	} else if (!strcmp(s, "-v") || !strcmp(s, "--version")) { 
	    opt |= OPT_VERSION;
	} else if (!strcmp(s, "-r") || !strncmp(s, "--run", 5)) { 
	    opt |= OPT_RUNIT;
	} else if (!strcmp(s, "-d") || !strncmp(s, "--db", 4)) { 
	    opt |= OPT_DBOPEN;
	} else if (!strcmp(s, "-w") || !strncmp(s, "--webdb", 7)) { 
	    opt |= OPT_WEBDB;
	} else if (!strcmp(s, "-c") || !strncmp(s, "--dump", 6)) {
	    opt |= OPT_DUMP;
	} else if (!strcmp(s, "-q") || !strcmp(s, "--quiet")) { 
	    opt |= OPT_QUIET;
	} else if (!strcmp(s, "-m") || !strcmp(s, "--makepkg")) { 
	    opt |= OPT_MAKEPKG;
	} else if (!strncmp(s, "--scriptopt=", 12)) {
	    *scriptval = atof(s + 12);
	} else if (*s == '-') {
	    /* not a valid option */
	    err = E_DATA;
	    break;
	} else if (!gotfile) {
	    strncat(fname, s, MAXLEN - 1);
	    gotfile = 1;
	}

	argc--;
    }

    if (!err) {
	err = incompatible_options(opt, OPT_BATCH | OPT_RUNIT | 
				   OPT_DBOPEN | OPT_WEBDB | OPT_MAKEPKG);
	if (!err) {
	    err = incompatible_options(opt, OPT_ENGLISH | OPT_BASQUE);
	}
    }

    *pargc = argc;
    *pargv = argv;
    *popt = opt;

    return err;
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
	gretl_errmsg_set(_("The shell command is not activated."));
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

static int print_command_param (const char *s, PRN *prn)
{
    int ret = 0;

    if (*s != '{' && strchr(s, ' ') != NULL) {
	ret += pprintf(prn, " \"%s\"", s);
    } else {
	ret += pprintf(prn, " %s", s);
    }

    return ret;
}

static int 
cmd_list_print_var (const CMD *cmd, int i, const DATASET *dset,
		    int gotsep, PRN *prn)
{
    int src, v = cmd->list[i];
    int imin = (MODEL_COMMAND(cmd->ci))? 1 : 0;
    int bytes = 0;

    if (v > 0 && i > imin && is_auto_generated_lag(i, cmd->list, cmd->linfo)) {
	if (is_first_lag(i, cmd->list, gotsep, cmd->linfo, &src)) {
	    bytes += print_lags_by_varnum(src, cmd->linfo, dset, 
					  gotsep, prn);
	} else if (cmd->ci == EQUATION && i == 1) {
	    pputc(prn, ' ');
	    bytes += 1 + pputs(prn, dset->varname[v]);
	}
    } else {
	pputc(prn, ' ');
	bytes += 1 + pputs(prn, dset->varname[v]);
    }

    return bytes;
}

static int more_coming (const CMD *cmd, int i, int gotsep)
{
    int ret = 0;

    if (cmd->opt) {
	ret = 1;
    } else if (cmd->linfo == NULL) {
	ret = (i < cmd->list[0]);
    } else {
	int j;

	for (j=i+1; j<=cmd->list[0] && !ret; j++) {
	    if (!is_auto_generated_lag(j, cmd->list, cmd->linfo) ||
		is_first_lag(j, cmd->list, gotsep, cmd->linfo, NULL)) {
		ret = 1;
	    }
	}
    }

    return ret;
}

static int n_separators (const int *list)
{
    int i, nsep = 0;

    for (i=2; i<=list[0]; i++) {
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
	} else if (!strcmp(cmd->param, "foreign")) {
	    ci = FOREIGN;
	} else if (!strcmp(cmd->param, "kalman")) {
	    ci = KALMAN;
	} else if (!strcmp(cmd->param, "mpi")) {
	    ci = MPI;
	}
    }

    return ci;
}

#define listsep_switch(c) (c == AR || c == MPOLS)

#define hold_param(c) (c == IVREG || c == AR || c == ARBOND ||		\
		       c == DPANEL || c == ARMA || c == CORRGM ||	\
		       c == PERGM || c == SCATTERS || c == MPOLS ||	\
                       c == GNUPLOT || c == GARCH || c == EQUATION ||	\
		       c == POISSON || c == XCORRGM || c == HECKIT ||	\
		       c == NEGBIN || c == DURATION || c == FRACTINT)

#define TESTLEN 62
#define LINELEN 78

static void
cmd_print_list (const CMD *cmd, const DATASET *dset,  
		int *plen, PRN *prn)
{
    char numstr[12];
    int use_varnames = (cmd->ci != AR && cmd->ci != DELEET);
    int nsep, gotsep, i;

    if (cmd->list == NULL || cmd->list[0] == 0) {
	return;
    }
    
    if (dset == NULL) {
	use_varnames = 0;
    } else if (matrix_data_option(cmd->ci, cmd->opt)) {
	/* using columns of a matrix */
	use_varnames = 0;
    }	

    nsep = n_separators(cmd->list);

    if (cmd->ci == LAGS) {
	if (cmd->param[0] != '\0') {
	    *plen += pprintf(prn, " %s;", cmd->param);
	}
    } else if (cmd->param[0] != '\0' && !hold_param(cmd->ci)) {
	*plen += print_command_param(cmd->param, prn);
    }

    if (cmd->ci == VECM && cmd->parm2 != NULL) {
	*plen += pprintf(prn, " %s", cmd->parm2);
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
	    *plen += cmd_list_print_var(cmd, i, dset, gotsep, prn);
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

    if (!strcmp(line, "set echo off") || !strcmp(line, "flush")) {
	return 1;
    }

    if (!strncmp(line, "quit", 4) && string_is_blank(line + 4)) {
	return 1;
    }

    if (cmd->ci == SET && !strcmp(cmd->param, "echo") &&
	gretl_function_depth() > 0) {
	return 1;
    }

    if (cmd->ci == OUTFILE && cmd->opt == OPT_C) {
	return 1;
    }

    if (*line == '!') {
	return 1;
    }

    return 0;
}

#define rewritten_lags(c) ((c->ci == ARMA || c->ci == DPANEL || c->ci == VAR) && \
                           c->parm2 != NULL &&				\
			   *c->parm2 != '\0')

/* these commands have sub-lists that may contain either
   numerical values or the names of scalar variables:
   this can't be handled properly by the list-printing
   mechanism */

#define dont_print_list(c) ((c->flags & CMD_NOLIST) ||	\
			    c->ci == ARBOND ||		\
			    c->ci == ARMA ||		\
			    c->ci == DPANEL ||		\
			    c->ci == GARCH ||		\
			    c->ci == OPEN)

#define print_param_last(c) (c == ARBOND ||	\
			     c == DPANEL ||	\
			     c == DELEET ||	\
	                     c == CORRGM ||	\
                             c == PERGM ||	\
	                     c == FRACTINT ||	\
                             c == XCORRGM)

/*
 * real_echo_cmd:
 * @cmd: pointer to #CMD struct.
 * @dset: pointer to dataset.
 * @line: "raw" command line associated with @cmd.
 * @recording: echo is going to command log (0/1).
 * @prn: pointer to gretl printing struct.
 *
 * Echoes the user command represented by @cmd and @line to
 * @prn.  This is used for two distinct purposes: to give 
 * visual feedback on the command supplied, and (in some 
 * contexts) to record a command that was executed interactively.
 */

static void real_echo_cmd (const CMD *cmd, const DATASET *dset, 
			   const char *line, int recording, 
			   PRN *prn)
{
    int compiling = 0;
    int skiplist = 0;
    int len, llen = 0;

    if (line == NULL || prn == NULL || cmd->ci > NC) {
	return;
    }

    if (gretl_compiling_function() || gretl_compiling_loop()) {
	compiling = 1;
    }

#if ECHO_DEBUG
    fprintf(stderr, "echo_cmd:\n*** line='%s'\n param='%s' extra='%s'\n", 
	    line, cmd->param, cmd->parm2);
    fprintf(stderr, " cmd->opt=%d, recording=%d, compiling=%d, nolist=%d\n",
	    cmd->opt, recording, compiling, cmd_nolist(cmd));
    fprintf(stderr, " cmd->word = '%s'\n", cmd->word);
    fprintf(stderr, " cmd->ci = %d, context = %d\n", cmd->ci, cmd->context);
    fprintf(stderr, " cmd->savename = '%s'\n", cmd->savename);
    if (!cmd_nolist(cmd)) {
	printlist(cmd->list, "cmd->list");
    }
#endif

    /* certain things don't get echoed at all, if not recording or
       compiling a function or loop */
    if (!recording && !compiling && command_is_silent(cmd, line)) {
	return;
    }

    /* in a progressive loop, do not apply the usual echo procedure
       for commands whose list may pertain to a temporary loop-special
       dataset */
    if (gretl_looping_progressive()) {
	if (cmd->ci == PRINT || cmd->ci == STORE) {
	    pprintf(prn, "? %s\n", line);
	    return;
	} 
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

    /* print leading string before echo? only if not recording */
    if (!recording) {
	if (compiling) {
	    llen += pputs(prn, "> ");
	} else {
	    llen += pputs(prn, "? ");
	}
    }

    /* special: printing a list */
    if (cmd->ci == PRINT && !strcmp(cmd->word, "list")) {
	pprintf(prn, "list %s print\n", cmd->parm2);
	return;
    }

    if (*line == '\0') {
	return;
    }

    /* command is preceded by a "savename" to which an object will
       be assigned */
    if (*cmd->savename && !cmd->context && cmd->ci != END) {
	if (strchr(cmd->savename, ' ') != NULL) {
	    pprintf(prn, "\"%s\" <- ", cmd->savename);
	    llen += strlen(cmd->savename) + 6;
	} else {
	    pprintf(prn, "%s <- ", cmd->savename);
	    llen += strlen(cmd->savename) + 4;
	}
    }

    skiplist = dont_print_list(cmd) || rewritten_lags(cmd);

    if (skiplist) {
	const char *s = line;
	
	if (rewritten_lags(cmd)) {
	    s = cmd->parm2;
	}
	if (strlen(s) > SAFELEN - 2) {
	    safe_print_line(s, &llen, prn);
	} else {
	    llen += pputs(prn, s);
	}
    } else {
	if (cmd->ci == EQUATION) {
	    llen += pprintf(prn, " %s", cmd->word);
	} else {
	    llen += pprintf(prn, "%s", cmd->word);
	}
	cmd_print_list(cmd, dset, &llen, prn);
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
	if (flagstr != NULL) {
	    len = strlen(flagstr);
	    if (llen + len > LINELEN) {
		pputs(prn, " \\\n ");
	    }	    
	    pputs(prn, flagstr);
	}
    }

    pputc(prn, '\n');
    gretl_print_flush_stream(prn);
}

void echo_command (const CMD *cmd, const DATASET *dset, 
		   const char *line, PRN *prn)
{
    real_echo_cmd(cmd, dset, line, 0, prn);
}

void gretl_record_command (const CMD *cmd, const DATASET *dset, 
			   const char *line, PRN *prn)
{
    real_echo_cmd(cmd, dset, line, 1, prn);
}

/* Look for a flag of the form " -x" which occurs outside of any
   quotes: if found, return a pointer to the flag.
*/

static const char *flag_present (const char *s, char f, int *quoted)
{
    const char *p = s;
    int inquote = 0;
    const char *ret = NULL;

#if CMD_DEBUG
    fprintf(stderr, "flag_present: looking at '%s'\n", s);
#endif

    while (*s) {
	if (*s == '"') {
	    inquote = !inquote;
	}
	if (!inquote) {
	    /* we're looking for, e.g., "-f", either at the
	       start of the input string or preceded by a
	       space -- and followed by something
	    */
	    if ((s == p || *(s-1) == ' ') && strlen(s) >= 3 &&
		*s == '-' && *(s+1) == f) {
		s += 2;
		s += strspn(s, " ");
		/* have we got some following text? */
		if (*s == '"' && *(s+1)) {
		    *quoted = 1;
		    ret = s + 1;
		}
		if (*s != '"' && *s != '\0') {
		    *quoted = 0;
		    ret = s;
		}
	    } 
	}
	if (ret != NULL) {
	    break;
	}
	s++;
    }

#if CMD_DEBUG
    fprintf(stderr, "flag_present: returning '%s'\n", ret);
#endif

    return ret;
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

/* grab filename for TeX output, if applicable */

static void get_optional_filename (const char *s, CMD *cmd)
{
    char *p = get_flag_field(s, 'f');

    if (p != NULL && *p != '\0') {
	free(cmd->param);
	if (libset_get_bool(USE_CWD) || fname_has_path(p)) {
	    cmd->param = p;
	} else {
	    cmd->param = gretl_strdup_printf("%s%s", gretl_workdir(), p);
	    free(p);
	}
    }
}

static int set_var_info (const char *line, gretlopt opt, 
			 DATASET *dset)
{
    char *p, vname[VNAMELEN];
    int v, err = 0;

    if (dset == NULL || dset->varinfo == NULL) {
	return E_NODATA;
    }

    /* skip command word plus space */
    line += strcspn(line, " ");
    line += strspn(line, " ");

    if (gretl_scan_varname(line, vname) != 1) {
	return E_PARSE;
    }

    /* skip varname, but not following space */
    line += strcspn(line, " ");

    if (gretl_is_scalar(vname)) {
	v = gretl_scalar_get_value(vname, NULL);
	if (v < 0 || v >= dset->v) {
	    return E_UNKVAR;
	}
    } else if (integer_string(vname)) {
	v = atoi(vname);
	if (v < 0 || v >= dset->v) {
	    return E_UNKVAR;
	}	
    } else {
	v = series_index(dset, vname);
	if (v < 0 || v >= dset->v) {
	    gretl_errmsg_sprintf(_("Unknown variable '%s'"), vname);
	    return E_UNKVAR;
	}
    }

    if (opt & OPT_D) {
	series_set_discrete(dset, v, 1);
    } else if (opt & OPT_C) {
	series_set_discrete(dset, v, 0);
    }

    if (opt & OPT_I) {
	const char *s = get_optval_string(SETINFO, OPT_I);

	if (s == NULL) {
	    err = E_ARGS;
	} else {
	    series_record_label(dset, v, s);
	}
    } else if (strstr(line, " -d ")) {
	/* backward compatibility */
	p = get_flag_field(line, 'd');
	if (p != NULL) {
	    series_record_label(dset, v, p);
	    free(p);
	}
    }

    if (opt & OPT_G) {
	const char *s = get_optval_string(SETINFO, OPT_G);

	if (s == NULL) {
	    err = E_ARGS;
	} else {
	    series_record_display_name(dset, v, s);
	}
    } else if (strstr(line, " -n ")) {
	/* backward compatibility */
	p = get_flag_field(line, 'n');
	if (p != NULL) {
	    series_record_display_name(dset, v, p);
	    free(p);
	}
    }    

    return err;
}

static void showlabels (const int *list, const DATASET *dset, PRN *prn)
{
    const char *label;
    int i, v, nl = 0;

    if (dset == NULL || dset->v == 0) {
	pprintf(prn, _("No series are defined\n"));
	return;
    }

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v >= 0 && v < dset->v) {
	    label = series_get_label(dset, v);
	    if (*label != '\0') {
		nl++;
	    }
	}
    } 

    if (nl == 0) {
	pprintf(prn, "No labels\n");
	return;
    } 

    pprintf(prn, _("Listing labels for variables:\n"));

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v >= 0 && v < dset->v) {
	    label = series_get_label(dset, v);
	    if (*label != '\0') {
		pprintf(prn, " %s: %s\n", dset->varname[v], label);
	    }
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

    pputs(prn, str);
    pputc(prn, '\n');
}

static void outfile_redirect (PRN *prn, FILE *fp, gretlopt opt,
			      int *parms)
{
    if (opt & OPT_Q) {
	parms[0] = gretl_echo_on();
	parms[1] = gretl_messages_on();
	set_gretl_echo(0);
	set_gretl_messages(0);
    } else {
	parms[0] = parms[1] = -1;
    }
    print_start_redirection(prn, fp);
}

static void maybe_restore_vparms (int *parms)
{
    if (parms[0] == 1) {
	set_gretl_echo(1);
    }
    if (parms[1] == 1) {
	set_gretl_messages(1);
    }    
    parms[0] = parms[1] = -1;
}

static int 
do_outfile_command (gretlopt opt, const char *fname, PRN *prn)
{
    static char outname[MAXLEN];
    static int vparms[2];
    int diverted = 0;
    int err = 0;

    if (prn == NULL) {
	return 0;
    }

    if (!(opt & (OPT_W | OPT_A | OPT_C))) {
	return E_ARGS;
    }

    diverted = printing_is_redirected(prn);

    /* command to close outfile */
    if (opt & OPT_C) {
	if (!diverted) {
	    pputs(prn, _("Output is not currently diverted to file\n"));
	    return 1;
	} else {
	    print_end_redirection(prn);
	    maybe_restore_vparms(vparms);
	    if (gretl_messages_on() && *outname != '\0') {
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
	outfile_redirect(prn, NULL, opt, vparms);
	*outname = '\0';
    } else if (!strcmp(fname, "stderr")) {
	outfile_redirect(prn, stderr, opt, vparms);
	*outname = '\0';
    } else if (!strcmp(fname, "stdout")) {
	outfile_redirect(prn, stdout, opt, vparms);
	*outname = '\0';
    } else {
	FILE *fp;

	fname = gretl_maybe_switch_dir(fname);

	if (opt & OPT_W) {
	    fp = gretl_fopen(fname, "w");
	} else {
	    fp = gretl_fopen(fname, "a");
	}

	if (fp == NULL) {
	    pprintf(prn, _("Couldn't open %s for writing\n"), fname);
	    return 1;
	}

	if (gretl_messages_on()) {
	    if (opt == OPT_W) {
		pprintf(prn, _("Now writing output to '%s'\n"), fname);
	    } else {
		pprintf(prn, _("Now appending output to '%s'\n"), fname);
	    }
	    
	}

	outfile_redirect(prn, fp, opt, vparms);
	strcpy(outname, fname);
    }

    return err;
}

int call_pca_plugin (VMatrix *cmat, DATASET *dset, 
		     gretlopt opt, PRN *prn)
{
    void *handle = NULL;
    int (*pca_from_cmatrix) (VMatrix *, DATASET *,
			     gretlopt, PRN *);
    int err = 0;

    gretl_error_clear();
    
    pca_from_cmatrix = get_plugin_function("pca_from_cmatrix", &handle);
    if (pca_from_cmatrix == NULL) {
        return 1;
    }
        
    err = (*pca_from_cmatrix) (cmat, dset, opt, prn);
    close_plugin(handle);
    
    return err;
}

static int do_pca (int *list, DATASET *dset,
		   gretlopt opt, PRN *prn)
{
    int err = 0;

    if (list[0] > 0) {
	VMatrix *cmat;

	/* adding OPT_U ensures a uniform sample for the
	   correlation or covariance matrix */
	cmat = corrlist(list, dset, opt | OPT_U, &err);
	if (!err) {
	    err = call_pca_plugin(cmat, dset, opt, prn);
	    if (!err && (opt & (OPT_O | OPT_A))) {
		/* results saved as series */
		maybe_list_vars(dset, prn);
	    }
	    free_vmatrix(cmat);
	}
    }

    return err;
}

static void print_info (gretlopt opt, DATASET *dset, PRN *prn)
{
    if (dset != NULL && dset->descrip != NULL) {
	pprintf(prn, "%s\n", dset->descrip);
    } else {
	pputs(prn, _("No data information is available.\n"));
    }
}

/* After estimating a model, check its errcode member to see
   if anything went wrong, and reset gretl_errno to zero.

   If we're looping (that is, if a loop is in progress at the
   current level of function execution), that's all, but if
   not then:

   (a) print the model (this requires special handling inside
   loops); 

   (b) if the user has employed the "name <- command" mechanism,
   attach the supplied name to the model; 

   (c) conditionally add the model to the stack in objstack.c,
   and if this is done, signal the fact by setting the 'pmod'
   member of @ExecState.

   (d) if we're called by the GUI program and the model has
   been assigned a name, activate the callback that adds the 
   model to the GUI session.
*/

static int print_save_model (MODEL *pmod, DATASET *dset,
			     gretlopt opt, PRN *prn, 
			     ExecState *s)
{
    int err = pmod->errcode;

    if (!err) {
	set_gretl_errno(0);
	if (!gretl_looping_currently()) {
	    int havename = *s->cmd->savename != '\0';

	    if (havename) {
		gretl_model_set_name(pmod, s->cmd->savename);
	    }
	    printmodel(pmod, dset, opt, prn);
	    attach_subsample_to_model(pmod, dset);
	    s->pmod = maybe_stack_model(pmod, s->cmd, prn, &err);
	    if (!err && s->callback != NULL && havename && 
		gretl_in_gui_mode()) {
		s->callback(s, s->pmod, GRETL_OBJ_EQN);
	    }
	}
    } 

    return err;
}

static void save_var_vecm (ExecState *s)
{
    maybe_stack_var(s->var, s->cmd);

    if (s->callback != NULL && *s->cmd->savename != '\0' &&
	gretl_in_gui_mode()) {
	s->callback(s, s->var, GRETL_OBJ_VAR);
    }    
}

static void gui_save_system (ExecState *s)
{
    /* note: with GRETL_OBJ_SYS, the business of calling
       "maybe_stack" is handled within system.c, so here
       all we have to do is invoke the GUI callback, if
       appropriate
    */
    if (gretl_in_gui_mode() && s->callback != NULL && 
	*s->cmd->savename != '\0') {
	s->callback(s, s->sys, GRETL_OBJ_SYS);
    }    
}

static int model_test_check (CMD *cmd, DATASET *dset, PRN *prn)
{
    int err = last_model_test_ok(cmd->ci, cmd->opt, dset, prn);

    if (err == E_DATA && cmd->ci == RESTRICT && *cmd->param == '\0') {
	/* try for a not-yet estimated anonymous system */
	if (get_anonymous_equation_system() != NULL) {
	    gretl_error_clear();
	    err = 0;
	}
    }

    return err;
}

static int get_line_continuation (char *line, FILE *fp, PRN *prn)
{
    char tmp[MAXLINE];
    int err = 0;

    if (!strncmp(line, "quit", 4)) {
	return 0;
    }

    while (top_n_tail(line, MAXLINE, &err)) {
	if (err) {
	    break;
	}
	*tmp = '\0';
	if (fgets(tmp, sizeof tmp, fp) && *tmp != '\0') {
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
		       DATASET *dset, gretlopt opt, 
		       PRN *prn)
{
    int indent = gretl_if_state_record();
    int echo = gretl_echo_on();
    int messages = gretl_messages_on();
    FILE *fp;
    int iferr, err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	gretl_errmsg_sprintf(_("Couldn't open %s"), fname);
	return E_FOPEN;
    }

    strcpy(s->runfile, fname);

    if (opt & OPT_Q) {
	set_gretl_echo(0);
	set_gretl_messages(0);
    }

    if (gretl_echo_on()) {
	pprintf(prn, "run \"%s\"\n", fname);
    }

    while (fgets(s->line, MAXLINE - 1, fp) && !err) {
	err = get_line_continuation(s->line, fp, prn);
	if (!err) {
	    err = maybe_exec_line(s, dset);
	}
    }

    fclose(fp);

    if (opt & OPT_Q) {
	set_gretl_echo(echo);
	set_gretl_messages(messages);
    }    

    iferr = gretl_if_state_check(indent);
    if (iferr && !err) {
	err = iferr;
    }

    return err;
}

static int lib_try_http (const char *s, char *fname, int *http)
{
    int err = 0;

    /* skip past command word */
    s += strcspn(s, " ");
    s += strspn(s, " ");

    if (strncmp(s, "http://", 7) == 0) {
#ifdef USE_CURL
	err = retrieve_public_file(s, fname);
	if (!err) {
	    *http = 1;
	}
#else
	gretl_errmsg_set(_("Internet access not supported"));
	err = E_DATA;
#endif
    }

    return err;
}

static int lib_clear_data (ExecState *s, DATASET *dset)
{
    int err = 0;

    if (dset->Z != NULL) {
	err = restore_full_sample(dset, NULL); 
	free_Z(dset); 
    }

    clear_model(s->model);
    clear_datainfo(dset, CLEAR_FULL);
    libgretl_session_cleanup(SESSION_CLEAR_DATASET);
    set_model_count(0);
    gretl_cmd_destroy_context(s->cmd);

    return err;
}

static int join_aggregation_method (const char *s, int *seqval,
				    char **auxname, int *err)
{
    int ret = -1;

    if (!strncmp(s, "seq:", 4)) {
	char *endptr;

	*seqval = (int) strtol(s + 4, &endptr, 10);
	if (*endptr == '\0' && *seqval != 0) {
	    ret = AGGR_SEQ;
	} else {
	    gretl_errmsg_sprintf(_("%s: invalid input '%s'\n"), "--seq", s + 4);
	    *err = E_DATA;
	}
    } else if (!strcmp(s, "count")) {
	ret = AGGR_COUNT;
    } else if (!strcmp(s, "avg")) {
	ret = AGGR_AVG;
    } else if (!strcmp(s, "sum")) {
	ret = AGGR_SUM;
    } else if (!strcmp(s, "min")) {
	ret = AGGR_MIN;
    } else if (!strcmp(s, "max")) {
	ret = AGGR_MAX;
    } else if (!strcmp(s, "none")) {
	ret = AGGR_NONE;
    } else if (!strncmp(s, "min(", 4) ||
	       !strncmp(s, "max(", 4)) {
	const char *p = strchr(s + 4, ')');

	if (p != NULL && strlen(p) == 1) {
	    int len = p - (s + 4);

	    if (len > 0) {
		*auxname = gretl_strndup(s + 4, len);
		if (*auxname == NULL) {
		    *err = E_ALLOC;
		} else {
		    ret = (s[1] == 'a')? AGGR_MAX : AGGR_MIN;
		}
	    }
	} else {
	    *err = E_PARSE;
	}
    } else {
	*err = E_PARSE;
    }

    return ret;
}

static int get_inner_key_id (const char *s, int n,
			     const DATASET *dset,
			     int *err)
{
    char vname[VNAMELEN];
    int id = -1;

    if (n == 0 || n >= VNAMELEN) {
	*err = E_PARSE;
    } else {
	*vname = '\0';
	strncat(vname, s, n);
	if (gretl_namechar_spn(vname) != n) {
	    gretl_errmsg_sprintf(_("field '%s' in command is invalid"), vname);
	    *err = E_PARSE;
	} else {
	    id = current_series_index(dset, vname);
	    if (id < 0) {
		*err = E_UNKVAR;
	    }
	}
    }

    return id;
}

static int *get_inner_keys (const char *s, DATASET *dset,
			    int *err)
{
    int *klist = NULL;
    int ikey1 = -1, ikey2 = -1;
    int nkeys = 0;

    if (strchr(s, ',') == NULL) {
	/* just one key, fine */
	ikey1 = current_series_index(dset, s);
	if (ikey1 < 0) {
	    *err = E_UNKVAR;
	} else {
	    nkeys = 1;
	}
    } else {
	/* we should have a double key */
	int n = strcspn(s, ",");

	ikey1 = get_inner_key_id(s, n, dset, err);

	if (!*err) {
	    s += n + 1;
	    n = strlen(s);
	    ikey2 = get_inner_key_id(s, n, dset, err);
	}

	if (!*err) {
	    nkeys = 2;
	}
    }

    if (!*err) {
	klist = gretl_list_new(nkeys);
	if (klist == NULL) {
	    *err = E_ALLOC;
	} else {
	    klist[1] = ikey1;
	    if (nkeys == 2) {
		klist[2] = ikey2;
	    }
	}
    }

    return klist;
}

static int lib_join_data (ExecState *s,
			  char *newfile,
			  DATASET *dset,
			  gretlopt opt,
			  PRN *prn)
{
    gretlopt opts[] = { 
	OPT_I, /* ikey: inner key(s) */
	OPT_O, /* okey: outer key(s) */
	OPT_F, /* filter: filter expression */
	OPT_A, /* aggr: aggregation method */
	OPT_D, /* data: "payload" spec */
	OPT_K, /* tkey: outer time-key name,format */
	OPT_X, /* tconvert: date columns for conversion */
	OPT_T, /* tconv-fmt: format for "tconvert" */ 
	0 
    };
    const char *param;
    char *p, *okey = NULL, *filter = NULL;
    char *varname = NULL, *dataname = NULL;
    char *auxname = NULL;
    char *tconvstr = NULL;
    char *tconvfmt = NULL;
    int *ikeyvars = NULL;
    int aggr = 0, seqval = 0;
    int tseries = 0;
    int i, err = 0;

    if (opt & OPT_K) {
	/* --tkey implies special handling of keys */
	if (opt & (OPT_I | OPT_O)) {
	    return E_BADOPT;
	} else if (!dataset_is_time_series(dset)) {
	    return E_PDWRONG;
	}
    }

    p = strstr(s->line, s->cmd->param);
    if (p == NULL) {
	return E_DATA;
    }

    p += strlen(s->cmd->param);
    p += strspn(p, " \"");
    if (*p == '\0') {
	return E_ARGS;
    } else if (strlen(p) != gretl_namechar_spn(p)) {
	return E_PARSE;
    }

    tseries = dataset_is_time_series(dset);

    varname = gretl_strdup(p);
    
    if (current_series_index(dset, varname) < 0) {
	err = check_varname(varname);
	if (!err && gretl_type_from_name(varname, NULL)) {
	    err = E_TYPES;
	}
    }

    for (i=0; opts[i] && !err; i++) {
	gretlopt jopt = opts[i];

	if (opt & jopt) {
	    param = get_optval_string(JOIN, jopt);
	    if (param == NULL) {
		err = E_DATA;
	    } else if (jopt == OPT_I) {		
		/* the inner key(s) string */
		ikeyvars = get_inner_keys(param, dset, &err);
	    } else if (jopt == OPT_O) {
		/* the outer key(s) string */
		okey = gretl_strdup(param);
	    } else if (jopt == OPT_F) {
		/* string specifying a row filter */
		filter = gretl_strdup(param);
	    } else if (jopt == OPT_A) {
		/* aggregation */
		aggr = join_aggregation_method(param, &seqval, 
					       &auxname, &err);
	    } else if (jopt == OPT_D) {
		/* string specifying the wanted data series */
		dataname = gretl_strdup(param);
	    } else if (jopt == OPT_K) {
		/* string specifying outer time key */
		okey = gretl_strdup(param);
	    } else if (jopt == OPT_X) {
		/* string holding list of time/date cols */
		tconvstr = gretl_strdup(param);
	    } else if (jopt == OPT_T) {
		/* format for tconvert columns */
		tconvfmt = gretl_strdup(param);
	    }
	}
    }

    if (!err && okey != NULL && ikeyvars == NULL && !(opt & OPT_K)) {
	/* We can't have an outer key but no inner one, unless
	   we're matching by the time-series structure of the
	   left-hand dataset (implied by OPT_K)
	*/
	gretl_errmsg_set(_("Inner key is missing"));
	err = E_PARSE;
    }

    if (!err && aggr != 0 && ikeyvars == NULL && !tseries) {
	/* aggregation requires ikeyvars, unless there's
	   an implicit time-series inner key
	*/
	gretl_errmsg_set(_("Inner key is missing"));
	err = E_ARGS;
    }

    if (!err) {
	PRN *vprn = (opt & OPT_V)? prn : NULL;

	err = join_from_csv(newfile, varname, dset, 
			    ikeyvars, okey, filter,
			    dataname, aggr, seqval, 
			    auxname, tconvstr,
			    tconvfmt, opt, vprn);
    }	

    free(ikeyvars);
    free(okey);
    free(filter);
    free(dataname);
    free(varname);
    free(auxname);
    free(tconvstr);
    free(tconvfmt);

    return err;
}

#define ALLOW_GUI_OPEN 1

static int lib_open_append (ExecState *s, 
			    DATASET *dset, 
			    char *newfile,
			    PRN *prn)
{
    CMD *cmd = s->cmd;
    gretlopt opt = cmd->opt;
    char *line = s->line;
    PRN *vprn = prn;
    int quiet = (opt & OPT_Q);
    int http = 0;
    int dbdata = 0;
    int odbc = 0;
    int ftype;
    int err = 0;
    
    if (cmd->ci == JOIN && (dset == NULL || dset->v == 0)) {
	return E_NODATA;
    }

#if ALLOW_GUI_OPEN
    if (cmd->ci == OPEN && gretl_function_depth() > 0) {
	gretl_errmsg_sprintf(_("The \"%s\" command cannot be used in this context"),
			     gretl_command_word(cmd->ci));
	return E_DATA;
    }
#else
    if (cmd->ci == OPEN && (gretl_in_gui_mode() || gretl_function_depth() > 0)) {
	gretl_errmsg_sprintf(_("The \"%s\" command cannot be used in this context"),
			     gretl_command_word(cmd->ci));
	return E_DATA;
    }
#endif

    if (cmd->ci != JOIN && (opt & OPT_O)) {
	odbc = 1;
    }

    err = lib_try_http(line, newfile, &http);
    if (err) {
	errmsg(err, prn);
	return err;
    }

    if (!http && !odbc) {
	/* not using http or ODBC */
	err = getopenfile(line, newfile, (opt & OPT_W)?
			  OPT_W : OPT_NONE);
	if (err) {
	    errmsg(err, prn);
	    return err;
	}
    }

    if (opt & OPT_W) {
	ftype = GRETL_NATIVE_DB_WWW;
    } else if (odbc) {
	ftype = GRETL_ODBC;
    } else {
	ftype = detect_filetype(newfile, OPT_P);
    }

    if (cmd->ci == JOIN) {
	if (ftype == GRETL_CSV) {
	    err = lib_join_data(s, newfile, dset, opt, prn);
	} else {
	    /* only CSV allowed for now */
	    err = E_NOTIMP;
	}
	if (err) {
	    errmsg(err, prn);
	}
	return err;
    }

    dbdata = (ftype == GRETL_NATIVE_DB || ftype == GRETL_NATIVE_DB_WWW ||
	      ftype == GRETL_RATS_DB || ftype == GRETL_PCGIVE_DB ||
	      ftype == GRETL_ODBC);

    if (cmd->ci == OPEN && !dbdata) {
	lib_clear_data(s, dset);
    } 

    if (quiet) {
	/* in case we hit any problems below... */
	vprn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
    } 

    if (ftype == GRETL_CSV) {
	err = import_csv(newfile, dset, opt, vprn);
    } else if (SPREADSHEET_IMPORT(ftype)) {
	err = import_spreadsheet(newfile, ftype, cmd->list, cmd->parm2,
				 dset, opt, vprn);
    } else if (OTHER_IMPORT(ftype)) {
	err = import_other(newfile, ftype, dset, opt, vprn);
    } else if (ftype == GRETL_XML_DATA) {
	err = gretl_read_gdt(newfile, dset, opt, vprn);
    } else if (ftype == GRETL_ODBC) {
	err = set_odbc_dsn(line, vprn);
    } else if (dbdata) {
	err = set_db_name(newfile, ftype, vprn);
    } else {
	err = gretl_get_data(newfile, dset, opt, vprn);
    }

    if (vprn != prn) {
	if (err) {
	    /* The user asked for --quiet operation, but something
	       went wrong so let's print any info we got on
	       vprn.
	    */
	    const char *buf = gretl_print_get_buffer(vprn);

	    if (buf != NULL && *buf != '\0') {
		pputs(prn, buf);
	    }
	} else {
	    /* print minimal success message */
	    pprintf(prn, _("Read datafile %s\n"), newfile);
	}
	gretl_print_destroy(vprn);
    }

    if (err) {
	errmsg(err, prn);
	return err;
    }

    if (dset->v > 0 && !dbdata && !quiet) {
	varlist(dset, prn);
    }

    if (http) {
	remove(newfile);
    }

    if (dbdata || http || cmd->ci == JOIN) {
	/* signal to the gretlcli callback that we didn't do
	   a regular datafile open 
	*/
	*newfile = '\0';
    }

    return err;
}

static int check_clear_data (void)
{
#if ALLOW_GUI_OPEN
    if (gretl_function_depth() > 0) {
	gretl_errmsg_sprintf(_("The \"%s\" command cannot be used in this context"),
			     gretl_command_word(CLEAR));
	return E_DATA;
    }
#else
    if (gretl_in_gui_mode() || gretl_function_depth() > 0) {
	gretl_errmsg_sprintf(_("The \"%s\" command cannot be used in this context"),
			     gretl_command_word(CLEAR));
	return E_DATA;
    }
#endif

    return 0;
}

static void schedule_callback (ExecState *s)
{
    if (s->callback != NULL) {
	s->flags |= CALLBACK_EXEC;
    } 
}

static int callback_scheduled (ExecState *s)
{
    return (s->flags & CALLBACK_EXEC) ? 1 : 0;
}

static void callback_exec (ExecState *s, char *fname, int err)
{
    if (!err && s->callback != NULL) {

	if (s->cmd->ci == OPEN) {
	    s->callback(s, fname, 0);
	} else {
	    s->callback(s, NULL, 0);
	}
    }

    s->flags &= ~CALLBACK_EXEC;
    *s->cmd->savename = '\0';
}

static int do_end_restrict (ExecState *s, DATASET *dset)
{
    GretlObjType otype = gretl_restriction_get_type(s->rset);
    gretlopt ropt = gretl_restriction_get_options(s->rset);
    gretlopt opt = s->cmd->opt | ropt;
    int err = 0;

    if (opt & OPT_F) {
	/* restrict --full */
	if (otype == GRETL_OBJ_VAR) {
	    s->var = gretl_restricted_vecm(s->rset, dset, 
					   opt, s->prn, &err);
	    if (!err && s->var != NULL) {
		save_var_vecm(s);
	    }
	} else if (otype == GRETL_OBJ_EQN) {
	    err = gretl_restriction_finalize_full(s, s->rset, dset, 
						  opt, s->prn);
	    if (!err) {
		gretlopt printopt = OPT_NONE;

		if (opt & (OPT_Q | OPT_S)) {
		    printopt = OPT_Q;
		}
		print_save_model(s->pmod, dset, printopt, s->prn, s);
	    }
	}
    } else {
	err = gretl_restriction_finalize(s->rset, dset, 
					 opt, s->prn);
    }

    s->rset = NULL;

    return err;
}

static int do_debug_command (ExecState *state, const char *param, 
			     gretlopt opt)
{
    int err = incompatible_options(opt, OPT_C | OPT_N | OPT_Q);

    if (err) {
	return err;
    }

    if (opt & (OPT_C | OPT_N)) {
	/* continue, next */
	if (!(state->flags & DEBUG_EXEC)) {
	    gretl_errmsg_set("Debugging is not in progress");
	    return E_DATA;
	} else {
	    /* handled in debug_command_loop */
	    return 0;
	}
    } else {
	/* OPT_Q quits debugging of the given function */
	return user_function_set_debug(param, !(opt & OPT_Q));
    } 
}

/* Given the name of a discrete variable, perform a command for each
   value of the discrete variable. Note that at present the only
   command supported in this way is SUMMARY.  
*/

static int do_command_by (CMD *cmd, DATASET *dset, PRN *prn)
{
    const char *byvar = get_optval_string(cmd->ci, OPT_B);
    const char **labels = NULL;
    gretl_matrix *xvals = NULL;
    const double *x;
    int i, v, nvals = 0;
    int single, err = 0;

    if (dset == NULL || byvar == NULL) {
	return E_DATA;
    }

    /* FIXME accept "unit" and "time"/"period" in place of actual
       variables for panel data? */

    v = current_series_index(dset, byvar);
    if (v < 0) {
	return E_UNKVAR;
    }

    x = (const double *) dset->Z[v];

    if (!series_is_discrete(dset, v) && !gretl_isdiscrete(dset->t1, dset->t2, x)) {
	gretl_errmsg_sprintf(_("The variable '%s' is not discrete"), byvar);
	return E_DATA;
    }

    single = cmd->list[0] == 1;

    xvals = gretl_matrix_values(x + dset->t1, dset->t2 - dset->t1 + 1, 
				OPT_S, &err);

    if (!err) {
	nvals = gretl_vector_get_length(xvals);
	if (nvals == 0) {
	    err = E_DATA;
	} else if (is_string_valued(dset, v)) {
	    int n_labels;

	    labels = series_get_string_vals(dset, v, &n_labels);
	    if (n_labels != nvals) {
		labels = NULL;
	    }
	}
    }

    if (!err && single) {
	pputc(prn, '\n');
	pprintf(prn, _("Summary statistics for %s, by value of %s"),
		dset->varname[cmd->list[1]], dset->varname[v]);
	pputc(prn, '\n');
    }

    for (i=0; i<nvals && !err; i++) {
	Summary *summ = NULL;
	char genline[64];
	double xi = gretl_vector_get(xvals, i);
	double *rv = NULL;

	gretl_push_c_numeric_locale();
	sprintf(genline, "%s == %g", byvar, xi);
	gretl_pop_c_numeric_locale();
	rv = generate_series(genline, dset, prn, &err);

	if (!err) {
	    summ = get_summary_restricted(cmd->list, dset, rv,
					  cmd->opt, prn, &err);
	}

	if (!err) {
	    if (i == 0) {
		pputc(prn, '\n');
	    }
	    if (single) {
		bufspace(2, prn);
	    }
	    if (labels != NULL) {
		pprintf(prn, "%s = %s (n = %d):\n", byvar, labels[i], summ->n);
	    } else {
		pprintf(prn, "%s = %g (n = %d):\n", byvar, xi, summ->n);
	    }
	    print_summary(summ, dset, prn);
	    free_summary(summ);
	}

	free(rv);
    }

    gretl_matrix_free(xvals);

    return err;
}

static void exec_state_prep (ExecState *s)
{
    s->flags &= ~CALLBACK_EXEC;
    s->pmod = NULL;
}

static int param_to_order (const char *s)
{
    if (s == NULL || *s == '\0') {
	/* giving an order is optional */
	return 0;
    } else if (integer_string(s)) {
	return atoi(s);
    } else if (gretl_is_scalar(s)) {
	return (int) gretl_scalar_get_value(s, NULL);
    } else {
	return -1;
    }
}

static GretlType get_type_for_deletion (const char *param, int *err)
{
    GretlType type = 0;

    if (param != NULL && *param != '\0') {
	/* not compatible with the type-deletion option */
	*err = E_DATA;
    } else {
	const char *s = get_optval_string(DELEET, OPT_T);

	if (s == NULL) {
	    *err = E_ARGS;
	} else {
	    type = gretl_type_from_string(s);
	}
    }

    return type;
}

#define can_continue(c) (c == ARMA || c == GARCH || c == GMM || \
                         c == MLE || c == NLS)

#define want_param_to_order(c) (c == CORRGM || c == XCORRGM ||	\
				c == PERGM || c == LAGS ||	\
				c == FRACTINT)

/* OMIT and ADD: if we're estimating a revised model, should
   we be saving it as the "last model", or are we just treating 
   the command as a stand-alone test? 
*/

static int add_omit_save (CMD *cmd)
{
    if (cmd->ci == ADD) {
	/* not saving if given the --lm option */
	return !(cmd->opt & OPT_L);
    } else {
	/* omit: not saving if given the --wald option */
	return !(cmd->opt & OPT_W);
    }
}

static int VAR_omit_driver (CMD *cmd, DATASET *dset, PRN *prn)
{
    GRETL_VAR *var = get_last_model(NULL);
    int err = 0;

    if (cmd->opt & OPT_W) {
	/* Wald test using VCV */
	err = gretl_VAR_wald_omit_test(var, cmd->list, dset, 
				       cmd->opt, prn);
    } else {
	/* the full deal: estimate reduced system */
	GRETL_VAR *vnew;

	vnew = gretl_VAR_omit_test(var, cmd->list, dset, cmd->opt,
				   prn, &err);
	if (!err) {
	    err = maybe_stack_var(vnew, cmd);
	}
    }

    return err;
}

static int model_print_driver (MODEL *pmod, DATASET *dset,
			       int ci, const char *param,
			       gretlopt opt, PRN *prn)
{
    int err = incompatible_options(opt, OPT_R | OPT_C);

    if (!err) {
	char fname[FILENAME_MAX];

	strcpy(fname, param);

	if (opt & OPT_R) {
	    err = rtfprint(pmod, dset, fname, opt);
	} else if (opt & OPT_C) {
	    err = csvprint(pmod, dset, fname, opt);
	} else {
	    gretlopt texopt = opt;

	    if (ci == EQNPRINT) {
		texopt |= OPT_E;
	    }		
	    err = texprint(pmod, dset, fname, texopt);
	}
	if (!err) {
	    pprintf(prn, _("Model printed to %s\n"), fname);
	}
    }

    return err;
}

static void abort_execution (ExecState *s)
{
    *s->cmd->savename = '\0';
    gretl_cmd_destroy_context(s->cmd);
    errmsg(E_STOP, s->prn);
}

static int plot_ok;

void set_plot_produced (void)
{
    plot_ok = 1;
}

static void maybe_schedule_graph_callback (ExecState *s)
{
    int gui_mode = gretl_in_gui_mode();

    if (graph_written_to_file()) {
	if (gui_mode && *s->cmd->savename != '\0') {
	    pprintf(s->prn, "Warning: ignoring \"%s <-\"\n", s->cmd->savename);
	}	
	report_plot_written(s->prn);
    } else if (gui_mode) {
	schedule_callback(s);
    }
}

static void maybe_print_error_message (CMD *cmd, int err, PRN *prn)
{
    if (gretl_function_depth() > 0) {
	; /* defer printing */
    } else if (cmd->flags & CMD_CATCH) {
	/* print only if messages on */
	if (gretl_messages_on()) {
	    errmsg(err, prn);
	}
    } else {
	/* otherwise go ahead and print */
	errmsg(err, prn);
    }
}

int gretl_cmd_exec (ExecState *s, DATASET *dset)
{
    CMD *cmd = s->cmd;
    char *line = s->line;
    MODEL *model = s->model;
    PRN *prn = s->prn;
    char readfile[MAXLEN];
    int *listcpy = NULL;
    int err = 0;

    exec_state_prep(s);
    plot_ok = 0;

    if (gretl_in_gui_mode() && check_for_stop()) {
	/* the GUI user clicked the "Stop" button */
	abort_execution(s);
	return E_STOP;
    }

    if (NEEDS_MODEL_CHECK(cmd->ci)) {
	err = model_test_check(cmd, dset, prn);
    } else if (MODIFIES_LIST(cmd->ci)) {
	if (cmd->list[0] == 0) {
	    /* no-op */
	    return 0;
	} else {
	    /* list is potentially modified -> make a copy */
	    listcpy = gretl_list_copy(cmd->list);
	    if (listcpy == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (err) {
	goto bailout;
    }

    *readfile = '\0';

    if (cmd->ci == OLS && dataset_is_panel(dset)) {
	cmd->ci = PANEL;
	cmd->opt |= OPT_P; /* panel pooled OLS flag */
    }

    if (want_param_to_order(cmd->ci)) {
	cmd->order = param_to_order(cmd->param);
    }

#if 0
    fprintf(stderr, "gretl_cmd_exec: '%s'\n", line);
#endif

    switch (cmd->ci) {

    case APPEND:
    case JOIN:
    case OPEN:
	err = lib_open_append(s, dset, readfile, prn);
	if (!err && cmd->ci == OPEN) {
	    schedule_callback(s);
	}
	break;

    case CLEAR:
	err = check_clear_data();
	if (!err) {
	    if (gretl_in_gui_mode()) {
		schedule_callback(s);
	    } else {
		lib_clear_data(s, dset);
	    }
	}
	break;

    case FLUSH:
	if (gretl_in_gui_mode()) {
	    schedule_callback(s);
	}
	break;

    case ANOVA:
	err = anova(cmd->list, dset, cmd->opt, prn);
	break;

    case ADF:
	err = adf_test(cmd->order, cmd->list, dset, cmd->opt, prn);
	break;

    case KPSS:
	err = kpss_test(cmd->order, cmd->list, dset, cmd->opt, prn);
	break;

    case LEVINLIN:
	err = llc_test_driver(cmd->param, cmd->list, dset, 
			      cmd->opt, prn);
	break;

    case COINT:
	err = engle_granger_test(cmd->order, cmd->list, dset, 
				 cmd->opt, prn);
	break;

    case COINT2:
	err = johansen_test_simple(cmd->order, cmd->list, 
				   dset, cmd->opt, prn);
	break;

    case CORR:
	err = incompatible_options(cmd->opt, OPT_U | OPT_S | OPT_K);
	if (err) {
	    break;
	}
	if (cmd->opt & OPT_K) {
	    err = kendall_tau(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->opt & OPT_S) {
	    err = spearman_rho(cmd->list, dset, cmd->opt, prn);
	} else {
	    err = gretl_corrmx(cmd->list, dset, cmd->opt, prn);
	}
	break;

    case CORRGM:
	err = corrgram(cmd->list[1], cmd->order, 0, dset, 
		       cmd->opt, prn);
	break;

    case XCORRGM:
	err = xcorrgram(cmd->list, cmd->order, dset, 
			cmd->opt, prn);
	break;

    case PERGM:
	err = periodogram(cmd->list[1], cmd->order, dset, 
			  cmd->opt, prn);
	break;

    case FRACTINT:
	err = fractint(cmd->list[1], cmd->order, dset, 
		       cmd->opt, prn);
	break;	

    case FUNDEBUG:
	err = do_debug_command(s, cmd->param, cmd->opt);
	break;

    case BREAK:
    case ENDLOOP:
	pprintf(prn, _("You can't end a loop here, "
		       "you haven't started one\n"));
	err = 1;
	break;

    case FCAST:
	err = do_forecast(line, dset, cmd->opt, prn);
	break;

    case FREQ:
	{
	    int graph = (s->flags == CONSOLE_EXEC ||
			 (cmd->opt & OPT_G));

	    if (cmd->opt & OPT_X) {
		err = matrix_freq_driver(cmd->list, &graph, 
					 cmd->opt, prn);
	    } else {
		err = freqdist(cmd->list[1], dset, &graph, 
			       cmd->opt, prn);
	    }
	    if (!err && graph) {
		maybe_schedule_graph_callback(s);
	    }
	}
	break;

    case DISCRETE:
	err = list_makediscrete(cmd->list, dset, cmd->opt);
	break;

    case ESTIMATE:
	err = estimate_named_system(line, dset, cmd->opt, prn);
	break;

    case FUNC:
	err = gretl_start_compiling_function(line, prn);
	break;

    case GENR:
	err = generate(line, dset, cmd->opt, prn);
	break;

    case PCA:
	err = do_pca(cmd->list, dset, cmd->opt, prn);
	break;

    case DATA:
	err = db_get_series(line, dset, cmd->opt, prn);
	break;

    case DATAMOD:
	err = modify_dataset(dset, cmd->aux, cmd->list, cmd->param, 
			     prn);
	if (!err) { 
	    schedule_callback(s);
	} 
	break;

    case DIFF:
    case LDIFF:
    case SDIFF:
	err = list_diffgenr(listcpy, cmd->ci, dset);
	if (!err) {
	    maybe_list_vars(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case ORTHDEV:
	err = list_orthdev(listcpy, dset);
	if (!err) {
	    maybe_list_vars(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case DUMMIFY:
	err = list_dumgenr(&listcpy, dset, cmd->opt);
	if (!err) {
	    maybe_list_vars(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case LAGS:
	err = list_laggenr(&listcpy, cmd->order, dset, OPT_NONE); 
	if (!err) {
	    maybe_list_vars(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case LOGS:
	err = list_loggenr(listcpy, dset);
	if (!err) {
	    maybe_list_vars(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case SQUARE:
	err = list_xpxgenr(&listcpy, dset, cmd->opt);
	if (!err) {
	    maybe_list_vars(dset, prn);
	    set_dataset_is_changed();
	}
	break;

    case PLOT:
	err = textplot(cmd->list, dset, cmd->opt, prn);
	break;

    case RMPLOT:
    case HURST:
	if (cmd->list[0] != 1) {
	    pputs(prn, _("This command requires one variable.\n"));
	    err = E_DATA;
	} else {
	    if (cmd->ci == RMPLOT) {
		err = rmplot(cmd->list, dset, cmd->opt, prn);
	    } else {
		err = hurstplot(cmd->list, dset, cmd->opt, prn);
	    }
	}
	break;

    case QQPLOT:
	err = qq_plot(cmd->list, dset, cmd->opt);
	break;	

    case INFO:
	print_info(cmd->opt, dset, prn);
	break;

    case RENAME:
	err = dataset_rename_series(dset, atoi(cmd->parm2), cmd->param);
	if (!err) {
	    maybe_list_vars(dset, prn);
	}
	break;

    case SET:
	err = execute_set_line(line, dset, cmd->opt, prn);
	break;

    case SETINFO:
	err = set_var_info(line, cmd->opt, dset);
	break;

    case SETMISS:
	if (dset == NULL || dset->Z == NULL) {
	    err = E_DATA;
	} else {
	    set_miss(cmd->list, cmd->param, dset, prn);
	}
        break;

    case LABELS:
	if (cmd->opt) {
	    err = read_or_write_var_labels(cmd->opt, dset, prn);
	    if (!err && (cmd->opt & (OPT_D | OPT_F))) {
		schedule_callback(s);
	    }
	} else {
	    showlabels(cmd->list, dset, prn);
	}
	break;

    case MARKERS:
	err = read_or_write_obs_markers(cmd->opt, dset, prn);
	if (!err && (cmd->opt & (OPT_D | OPT_F))) {
	    schedule_callback(s);
	}
	break;

    case VARLIST:
	if (cmd->opt & OPT_S) {
	    print_scalars(prn);
	} else if (cmd->opt & OPT_A) {
	    list_ok_dollar_vars(dset, prn);
	} else {
	    varlist(dset, prn);
	}
	break;

    case PRINT:
	if (cmd->param != NULL && *cmd->param != '\0') {
	    do_print_string(cmd->param, prn);
	} else {
	    printdata(cmd->list, cmd->parm2, dset, cmd->opt, prn);
	}
	break;

    case PRINTF:
    case SPRINTF:
    case SSCANF:
	err = do_printscan_command(line, dset, prn); 	 
	break;

    case PVAL:
	err = batch_pvalue(line, dset, prn);
	break;

    case SUMMARY:
	err = incompatible_options(cmd->opt, OPT_B | OPT_W);
	if (err) {
	    break;
	}
	if (cmd->opt & OPT_B) {
	    err = do_command_by(cmd, dset, prn);
	} else if (cmd->opt & OPT_X) {
	    err = matrix_command_driver(cmd->ci, cmd->list, cmd->param,
					dset, cmd->opt, prn);
	} else {
	    err = list_summary_driver(cmd->list, dset, cmd->opt, prn);
	}
	break; 

    case XTAB:
	if (cmd->opt & OPT_M) {
	    err = crosstab_from_matrix(cmd->opt, prn);
	} else {
	    err = crosstab(cmd->list, dset, cmd->opt, prn);
	}
	break;

    case MAHAL:
	err = mahalanobis_distance(cmd->list, dset, cmd->opt, prn);
	break;

    case MEANTEST:
	err = means_test(cmd->list, dset, cmd->opt, prn);
	break;	

    case VARTEST:
	err = vars_test(cmd->list, dset, prn);
	break;

    case RUNS:
	err = runs_test(cmd->list[1], dset, cmd->opt, prn);
	break;

    case SPEARMAN:
	err = spearman_rho(cmd->list, dset, cmd->opt, prn);
	break;

    case DIFFTEST:
	err = diff_test(cmd->list, dset, cmd->opt, prn);
	break;

    case OUTFILE:
	err = do_outfile_command(cmd->opt, cmd->param, prn);
	break;

    case SETOBS:
	err = set_obs(line, dset, cmd->opt);
	if (!err) {
	    if (dset->n > 0) {
		if (!(cmd->opt & (OPT_I | OPT_G))) {
		    print_smpl(dset, 0, prn);
		}
		schedule_callback(s);
	    } else {
		pprintf(prn, _("data frequency = %d\n"), dset->pd);
	    }
	}
	break;

    case SMPL:
	if (cmd->opt == OPT_F) {
	    err = restore_full_sample(dset, s);
	} else if (cmd->opt) {
	    err = restrict_sample(line, cmd->list, dset, 
				  s, cmd->opt, prn);
	} else { 
	    err = set_sample(line, dset);
	}
	if (!err) {
	    print_smpl(dset, get_full_length_n(), prn);
	}	
	break;

    case MAKEPKG:
	err = create_and_write_function_package(cmd->param, cmd->opt, prn);
	break;

    case STORE:
	if (dset == NULL || dset->Z == NULL) {
	    err = E_NODATA;
	} else if (*cmd->param == '\0') {
	    pputs(prn, _("store: no filename given\n"));
	    err = E_PARSE;
	}
	if (!err) {
	    err = write_data(cmd->param, cmd->list, dset, cmd->opt, prn);
	}
	break;

    case SHELL:
	err = gretl_shell(line, prn);
	break;

    case OLS:
    case WLS:
	clear_model(model);
	*model = lsq(cmd->list, dset, cmd->ci, cmd->opt);
	err = print_save_model(model, dset, cmd->opt, prn, s);
	break;
	
    case MPOLS:
	clear_model(model);
	*model = mp_ols(cmd->list, dset);
	err = print_save_model(model, dset, cmd->opt, prn, s);
	break;

    case AR:
    case AR1:
    case ARMA:
    case ARCH:
	clear_model(model);
	if (cmd->ci == AR) {
	    *model = ar_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == AR1) {
	    *model = ar1_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == ARMA) {
	    *model = arma(cmd->list, cmd->auxlist, dset, 
			  cmd->opt, prn);
	} else {
	    *model = arch_model(cmd->list, cmd->order, dset,
				cmd->opt);
	}
	err = print_save_model(model, dset, cmd->opt, prn, s);
	break;

    case ARBOND:
    case PANEL:	
    case DPANEL:
	if (!dataset_is_panel(dset)) {
	    gretl_errmsg_set(_("This estimator requires panel data"));
	    err = E_DATA;
	    break;
	}
    case GARCH:
    case HECKIT:
    case HSK:
    case INTREG:
    case IVREG:
    case LAD:
    case LOGISTIC:
    case LOGIT:
    case POISSON:
    case NEGBIN:
    case PROBIT:
    case QUANTREG:
    case TOBIT:
    case DURATION:
    case BIPROBIT:
	clear_model(model);
	if (cmd->ci == LOGIT || cmd->ci == PROBIT) {
	    *model = logit_probit(cmd->list, dset, cmd->ci, cmd->opt, prn);
	} else if (cmd->ci == HSK) {
	    *model = hsk_model(cmd->list, dset);
	} else if (cmd->ci == LOGISTIC) {
	    *model = logistic_driver(cmd->list, dset, cmd->opt);
	} else if (cmd->ci == TOBIT) {
	    *model = tobit_driver(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == POISSON || cmd->ci == NEGBIN) {
	    *model = count_model(cmd->list, cmd->ci, dset, cmd->opt, prn);
	} else if (cmd->ci == HECKIT) {
	    *model = heckit_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == IVREG) {
	    *model = ivreg(cmd->list, dset, cmd->opt);
	} else if (cmd->ci == LAD) {
	    *model = lad(cmd->list, dset);
	} else if (cmd->ci == QUANTREG) {
	    *model = quantreg_driver(cmd->param, cmd->list, dset,
				     cmd->opt, prn);
	} else if (cmd->ci == DURATION) {
	    *model = duration_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == GARCH) {
	    *model = garch(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == PANEL) {
	    *model = panel_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == ARBOND) {
	    *model = arbond_model(cmd->list, cmd->param, dset, 
				  cmd->opt, prn);
	} else if (cmd->ci == DPANEL) {
	    *model = dpd_model(cmd->list, cmd->auxlist, cmd->param, 
			       dset, cmd->opt, prn);
	} else if (cmd->ci == INTREG) {
	    *model = interval_model(cmd->list, dset, cmd->opt, prn);
	} else if (cmd->ci == BIPROBIT) {
	    *model = biprobit_model(cmd->list, dset, cmd->opt, prn);
	} else {
	    /* can't happen */
	    err = 1;
	    break;
	}
	err = print_save_model(model, dset, cmd->opt, prn, s);
	break;

    case GMM:
    case MLE:
    case NLS:
	err = nl_parse_line(cmd->ci, line, dset, prn);
	if (!err) {
	    gretl_cmd_set_context(cmd, cmd->ci);
	} 
	break;

    case FOREIGN:
    case MPI:
	err = foreign_append_line(line, cmd->opt, prn);
	if (!err) {
	    gretl_cmd_set_context(cmd, FOREIGN);
	}
	break;

    case KALMAN:
	err = kalman_parse_line(line, dset, cmd->opt);
	if (!err && (cmd->opt == OPT_NONE)) {
	    gretl_cmd_set_context(cmd, cmd->ci);
	}
	break;

    case ADD:
    case OMIT:
	if (get_last_model_type() == GRETL_OBJ_VAR) {
	    err = VAR_omit_driver(cmd, dset, prn);
	} else if (add_omit_save(cmd)) {
	    MODEL mymod;
	    
	    gretl_model_init(&mymod, dset);
	    if (cmd->ci == ADD) {
		err = add_test_full(model, &mymod, cmd->list, 
				    dset, cmd->opt, prn);
	    } else {
		err = omit_test_full(model, &mymod, cmd->list, 
				     dset, cmd->opt, prn);
	    }
	    if (!err) {
		gretlopt popt = OPT_NONE;

		if (cmd->opt & (OPT_I | OPT_Q)) {
		    popt = OPT_Q;
		} else if (cmd->opt & OPT_O) {
		    popt = OPT_O; /* --vcv printing option */
		}
		clear_model(model);
		*model = mymod;
		print_save_model(model, dset, popt, prn, s);
	    } 
	} else if (cmd->ci == ADD) {
	    err = add_test(model, cmd->list, dset, cmd->opt, prn);
	} else {
	    err = omit_test(model, cmd->list, dset, cmd->opt, prn);
	}
	if (err == E_NOOMIT) {
	    /* auto-omit was a no-op */
	    err = 0;
	}	
	break;	

    case COEFFSUM:
    case CUSUM:
    case RESET:
    case CHOW:
    case QLRTEST:
    case VIF:
	if (cmd->ci == COEFFSUM) {
	    err = gretl_sum_test(cmd->list, model, dset, prn);
	} else if (cmd->ci == CUSUM) {
	    err = cusum_test(model, dset, cmd->opt, prn);
	} else if (cmd->ci == RESET) {
	    err = reset_test(model, dset, cmd->opt, prn);
	} else if (cmd->ci == CHOW) {
	    err = chow_test_driver(line, model, dset, cmd->opt, prn);
	} else if (cmd->ci == QLRTEST) {
	    err = QLR_test(model, dset, cmd->opt, prn);
	} else if (cmd->ci == VIF) { 
	    err = vif_test(model, dset, prn);
	} 
	break;

    case NORMTEST:
	err = gretl_normality_test(cmd->list[1], dset, cmd->opt, prn);
	break;

    case HAUSMAN:
	err = panel_hausman_test(model, dset, cmd->opt, prn);
	break;

    case MODTEST:
	err = model_test_driver(cmd->param, dset, cmd->opt, prn);
	break;

    case LEVERAGE:
	err = leverage_test(model, dset, cmd->opt, prn);
	if (!err && (cmd->opt & OPT_S) && !(cmd->opt & OPT_Q)) {
	    maybe_list_vars(dset, prn);
	}
	break;

    case EQNPRINT:
    case TABPRINT:
	if (model->errcode == E_NAN) {
	    pprintf(prn, _("Couldn't format model\n"));
	} else {
	    err = model_print_driver(model, dset, cmd->ci,
				     cmd->param, cmd->opt, 
				     prn);
	}
	break;

    case RESTRICT:
	/* joint hypothesis test on model */
	if (s->rset == NULL) {
	    if (*cmd->param == '\0') {
		/* if param is non-blank, we're restricting a named system */
		err = model_test_check(cmd, dset, prn);
		if (err) break;
	    }
	    s->rset = restriction_set_start(line, cmd->opt, &err);
	    if (!err) {
		gretl_cmd_set_context(cmd, RESTRICT);
	    }
	} else {
	    err = restriction_set_parse_line(s->rset, line, dset);
	    if (err) {
		s->rset = NULL;
	    }	
	}
	break;

    case SYSTEM:
	if (s->sys == NULL) {
	    /* no equation system is defined currently */
	    s->sys = equation_system_start(line, cmd->savename, cmd->opt, &err);
	    if (!err) {
		gretl_cmd_set_context(cmd, SYSTEM);
	    }
	} else {
	    err = system_parse_line(s->sys, line, dset);
	    if (err) {
		s->sys = NULL;
	    } 
	}
	break;

    case EQUATION:
	if (cmd->opt & OPT_M) {
	    err = equation_system_append_multi(s->sys, cmd->param, dset);
	} else {
	    err = equation_system_append(s->sys, cmd->list);
	}
	if (err) {
	    s->sys = NULL;
	}
	break;

    case END:
	if (!strcmp(cmd->param, "system")) {
	    err = equation_system_finalize(s->sys, dset, cmd->opt, prn);
	    if (!err) {
		gui_save_system(s);
	    }
	    /* clear for next use */
	    s->sys = NULL;
	} else if (!strcmp(cmd->param, "mle") || 
		   !strcmp(cmd->param, "nls") ||
		   !strcmp(cmd->param, "gmm")) {
	    clear_model(model);
	    *model = nl_model(dset, cmd->opt, prn);
	    err = print_save_model(model, dset, cmd->opt, prn, s);
	} else if (!strcmp(cmd->param, "restrict")) {
	    err = do_end_restrict(s, dset);
	} else if (!strcmp(cmd->param, "foreign")) {
	    err = foreign_execute(dset, cmd->opt, prn);
	} else if (!strcmp(cmd->param, "kalman")) {
	    err = kalman_parse_line(line, dset, cmd->opt);
	} else if (!strcmp(cmd->param, "mpi")) {
	    err = foreign_execute(dset, cmd->opt, prn);
	} else {
	    err = 1;
	}
	break;

    case VAR:
    case VECM:
	if (cmd->ci == VAR) {
	    s->var = gretl_VAR(cmd->order, cmd->auxlist, cmd->list, 
			       dset, cmd->opt, prn, &err);
	} else {
	    int rank = gretl_int_from_string(cmd->parm2, &err);

	    if (!err) {
		s->var = gretl_VECM(cmd->order, rank, cmd->list, dset, 
				    cmd->opt, prn, &err);
	    }
	}
	if (!err && s->var != NULL) {
	    save_var_vecm(s);
	}
	break;

    case RUN:
    case INCLUDE:
	if (cmd->ci == RUN) {
	    err = getopenfile(line, readfile, OPT_S);
	} else {
	    err = getopenfile(line, readfile, OPT_I);
	    cmd->opt |= OPT_Q;
	}
	if (err) { 
	    break;
	} 
	if (gretl_messages_on()) {
	    pprintf(prn, " %s\n", readfile);
	}
	if (cmd->ci == INCLUDE && gretl_is_xml_file(readfile)) {
	    err = load_user_XML_file(readfile, prn);
	    break;
	}
	if (!strcmp(readfile, s->runfile)) { 
	    pprintf(prn, _("Infinite loop detected in script\n"));
	    err = 1;
	    break;
	}
	err = run_script(readfile, s, dset, cmd->opt, prn);
	break;

    case FUNCERR:
    case FUNCRET:
	if (gretl_function_depth() == 0) {
	    gretl_errmsg_sprintf("'%s': can only be used within a function",
				 gretl_command_word(cmd->ci));
	    err = 1;
	} else if (cmd->ci == FUNCERR) {
	    err = E_FUNCERR;
	} 
	break;

    case DELEET:
	if (cmd->list != NULL && cmd->list[0] > 0) {
	    if (cmd->opt & OPT_F) {
		/* got the --force option */
		err = dataset_drop_listed_variables(cmd->list, dset, 
						    NULL, prn);
	    } else {
		pputs(prn, _("You cannot delete series in this context\n"));
		err = 1;
	    }
	} else if (cmd->opt & OPT_T) {
	    GretlType type = get_type_for_deletion(cmd->param, &err);

	    if (!err) {
		err = delete_user_vars_of_type(type, prn);
	    }
	} else {
	    err = gretl_delete_var_by_name(cmd->param, prn);
	}
	break;

    case MODPRINT:
	err = do_modprint(line, cmd->opt, prn);
	break;

    case GNUPLOT:
    case BXPLOT:
    case SCATTERS:
	if (cmd->opt & OPT_X) {
	    err = matrix_command_driver(cmd->ci, cmd->list, cmd->param, 
					dset, cmd->opt, prn);
	} else if (cmd->ci == GNUPLOT) {
	    if (cmd->opt & OPT_D) {
		err = gnuplot_process_file(cmd->opt, prn);
	    } else if (cmd->opt & OPT_C) {
		err = xy_plot_with_control(cmd->list, cmd->param, 
					   dset, cmd->opt);
	    } else {
		err = gnuplot(cmd->list, cmd->param, dset, cmd->opt);
	    }
	} else if (cmd->ci == SCATTERS) {
	    err = multi_scatters(cmd->list, dset, cmd->opt);
	} else if (cmd_nolist(cmd)) {
	    err = boolean_boxplots(line, cmd->param, dset, cmd->opt);
	} else {
	    err = boxplots(cmd->list, cmd->param, dset, cmd->opt);
	}
	break;

    case MODELTAB:
    case GRAPHPG:
	if (gretl_in_gui_mode()) {
	    schedule_callback(s);
	} else {
	    pprintf(prn, _("%s: command not available\n"), cmd->word);
	}
	break;

    default:
	if (*cmd->word != '\0') {
	    pprintf(prn, _("Sorry, the %s command is not yet implemented "
			   "in libgretl\n"), cmd->word);
	} else {
	    pprintf(prn, "What?\n");
	}
	err = 1;
	break;
    }

    if (listcpy != NULL) {
	free(listcpy);
    }

    if (err == E_OK) {
	err = 0;
    }

    if (!err && plot_ok) {
	maybe_schedule_graph_callback(s);
    }

    if (callback_scheduled(s)) {
	callback_exec(s, readfile, err);
    } 

 bailout:

    if (err) {
	maybe_print_error_message(cmd, err, prn);
    }

    err = process_command_error(cmd, err);

    if (err) {
	gretl_cmd_destroy_context(cmd);
    } else {
	/* this is a no-op if there's no warning */
	warnmsg(prn);
    }

    return err;
}

/* called by functions, and by scripts executed from within
   functions */

int maybe_exec_line (ExecState *s, DATASET *dset)
{
    int err = 0;

    if (string_is_blank(s->line)) {
	return 0;
    }

    if (gretl_compiling_loop()) { 
	err = get_command_index(s->line, s->cmd);
    } else {
	/* FIXME last arg to parse_command_line() ? */
	err = parse_command_line(s->line, s->cmd, dset, NULL);
    }

    if (err) {
        errmsg(err, s->prn);
        return err;
    }

    gretl_exec_state_transcribe_flags(s, s->cmd);

    if (s->cmd->ci < 0) {
	return 0; /* nothing there, or a comment */
    }

    if (s->cmd->ci == LOOP || gretl_compiling_loop()) {  
	/* accumulating loop commands */
	err = gretl_loop_append_line(s, dset);
	if (err) {
	    errmsg(err, s->prn);
	    return err;
	} 
	return 0;
    } 

    s->pmod = NULL; /* be on the safe side */

    if (s->cmd->ci == FUNCERR) {
	err = E_FUNCERR;
    } else {
	/* note: error messages may be printed to s->prn */
	err = gretl_cmd_exec(s, dset);
    }

    return err;
}

static int could_be_varname (const char *s)
{
    int n = gretl_namechar_spn(s);
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

static int is_endloop (const char *s)
{
    char s1[4], s2[5];

    return sscanf(s, "%3s %4s", s1, s2) == 2 && !strcmp(s2, "loop");
}

/**
 * get_command_index:
 * @line: command line.
 * @cmd: pointer to gretl command struct.
 *
 * Parse @line and assign to the %ci field of @cmd the index number of
 * the command embedded in @line.  Note: this is a "lite" version of
 * parse_command_line().  It is used when commands are being stacked
 * for execution within a loop.  Command options are not parsed out of
 * @line.
 *
 * Returns: 1 on error, otherwise 0.
 */

int get_command_index (char *line, CMD *cmd)
{
    char cnext = 0;
    int done = 0;

    cmd->ci = 0;
    cmd->opt = OPT_NONE;
    *cmd->parm2 = *cmd->param = '\0';

    while (isspace(*line)) {
	line++;
    }

#if CMD_DEBUG
    fprintf(stderr, "get_command_index: line='%s'\n", line);
#endif

    if (filter_comments(line, cmd)) {
	return 0;
    }

    if (!strncmp(line, "catch ", 6)) {
	line += 6;
    }

    if (!get_first_word(line, &cnext, cmd)) {
	if (*line == '$' || *line == '@') {
	    /* most plausible possibility? */
	    strcpy(cmd->word, "genr");
	    cmd->ci = GENR;
	} else {
	    /* FIXME watch out for fallout! (2012-03-01) */
	    cmd_set_nolist(cmd);
	    cmd->ci = CMD_NULL;
#if CMD_DEBUG
	    fprintf(stderr, "get_command_word: got nothing!\n");
#endif
	    return E_PARSE; /* return 0 */
	}
    }

#if CMD_DEBUG
    fprintf(stderr, " got command word = '%s'\n", cmd->word);
#endif

    if (!strcmp(cmd->word, "end")) {
	if (is_endloop(line)) {
	    deprecate_alias("end loop", "endloop", 0);
	    cmd->ci = ENDLOOP;
	} else {
	    cmd->context = 0;
	    cmd->ci = END;
	}
	done = 1;
    } else if (cmd->context) {
	cmd->ci = cmd->context;
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
		cmd->opt = OPT_O;
	    } else {
		cmd->err = 1;
		gretl_errmsg_sprintf(_("command '%s' not recognized"), 
				     cmd->word);
		return 1;
	    }
	}
    }

    if (cmd->ci == NLS || cmd->ci == MLE ||
	cmd->ci == GMM || cmd->ci == FOREIGN ||
	cmd->ci == KALMAN) {
	cmd->context = cmd->ci;
    } else if (cmd->ci == MPI) {
	cmd->context = FOREIGN;
    }

#if CMD_DEBUG
    fprintf(stderr, " get_command_index: cmd->ci set to %d\n", cmd->ci);
#endif

    return 0;
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
    *cmd->savename = '\0';

    cmd->list = NULL;
    cmd->param = NULL;
    cmd->parm2 = NULL;
    cmd->auxlist = NULL;
    cmd->linfo = NULL;

    /* make 'list', 'param' and 'extra' blank rather than NULL
       for safety (in case they are dereferenced) */

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
	cmd->parm2 = calloc(1, 1);
	if (cmd->parm2 == NULL) {
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
    free(cmd->parm2);
    free(cmd->auxlist);

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
    *cmd->savename = '\0';
}

gretlopt gretl_cmd_get_opt (const CMD *cmd)
{
    return cmd->opt;
}

void gretl_cmd_set_opt (CMD *cmd, gretlopt opt)
{
    cmd->opt = opt;
}

const char *gretl_cmd_get_savename (CMD *cmd)
{
    return cmd->savename;
}

void gretl_exec_state_init (ExecState *s,
			    ExecFlags flags,
			    char *line,
			    CMD *cmd,
			    MODEL *model, 
			    PRN *prn)
{
    s->flags = flags;

    s->line = line;
    if (s->line != NULL) {
	*s->line = '\0';
    }    

    s->cmd = cmd;
    if (s->cmd != NULL) {
	s->cmd->ci = 0;
    }    

    *s->runfile = '\0';

    s->model = model;
    s->prn = prn;

    s->pmod = NULL;
    s->sys = NULL;
    s->rset = NULL;
    s->var = NULL;
    s->in_comment = 0;
    s->padded = 0;

    if (flags == FUNCTION_EXEC) {
	/* On entry to function execution we check if there's
	   a 'last model' in place. If so, we want to make
	   this invisible within the function, but set things
	   up so that we can restore it as last model on
	   exit from the function -- the idea being that
	   excuting a function should not change the 'last
	   model' state at caller level. To achieve this we
	   need to take out a 'private' reference to the
	   model, stored in the ExecState, and then remove
	   it from last model position for the present.
	*/
	s->prev_model = get_last_model(&s->prev_type);
	if (s->prev_model != NULL) {
	    gretl_object_ref(s->prev_model, s->prev_type);
	    set_as_last_model(NULL, GRETL_OBJ_NULL);
	}
	s->prev_model_count = get_model_count();
    } else {
	s->prev_model = NULL;
	s->prev_type = GRETL_OBJ_NULL;
	s->prev_model_count = -1;
    }

    s->submask = NULL;
    s->callback = NULL;
}

void function_state_init (CMD *cmd, ExecState *state, int *indent0)
{
    cmd->list = NULL;
    cmd->param = NULL;
    cmd->parm2 = NULL;
    cmd->linfo = NULL;

    state->cmd = NULL;
    state->model = NULL;
    state->submask = NULL;

    state->padded = 0;

    *indent0 = gretl_if_state_record();
}

static EXEC_CALLBACK gui_callback;

void gretl_exec_state_set_callback (ExecState *s, EXEC_CALLBACK callback,
				    gretlopt opt)
{
    s->callback = callback;
    s->pmod = NULL;

    if (opt & OPT_G) {
	gui_callback = callback;
    }
}

EXEC_CALLBACK get_gui_callback (void)
{
    return gui_callback;
}

void gretl_exec_state_clear (ExecState *s)
{
    gretl_cmd_free(s->cmd);

    if (s->flags & FUNCTION_EXEC) {
	/* Restore whatever was the 'last model' before 
	   function execution. Note that this includes
	   the case where there was no 'last model', in
	   which case we restore the null state. Drop
	   the extra refcount for the model we put into
	   last model position (if any), so we don't end 
	   up leaking memory.
	*/
	set_as_last_model(s->prev_model, s->prev_type);
	if (s->prev_model != NULL) {
	    gretl_object_unref(s->prev_model, s->prev_type);
	}
	/* restore the previous model count */
	if (s->prev_model_count >= 0) {
	    set_model_count(s->prev_model_count);
	}
    }

    destroy_working_model(s->model);

    s->prev_model = NULL;
    s->prev_type = GRETL_OBJ_NULL;
    s->prev_model_count = -1;

    free_subsample_mask(s->submask);
}

void gretl_exec_state_uncomment (ExecState *s)
{
    s->in_comment = 0;
    s->cmd->flags &= ~CMD_IGNORE;
}

void gretl_exec_state_transcribe_flags (ExecState *s, CMD *cmd)
{
    s->in_comment = (cmd_ignore(cmd))? 1 : 0;
}

void gretl_exec_state_set_model (ExecState *s, MODEL *pmod)
{
    s->pmod = pmod;
}

int process_command_error (CMD *cmd, int err)
{
    int ret = err;

    if (err) {
	if (gretl_compiling_function() ||
	    gretl_compiling_loop()) {
	    ; /* pass the error through */
	} else if (libset_get_bool(HALT_ON_ERR) == 0) {
	    /* global "continue on error" */
	    set_gretl_errno(err);
	    ret = 0;
	} else if (cmd->flags & CMD_CATCH) {
	    /* local "continue on error" */
	    set_gretl_errno(err);
	    cmd->flags ^= CMD_CATCH;
	    ret = 0;
	}
    }

    return ret;
}
