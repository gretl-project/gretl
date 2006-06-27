/*
 *  Copyright (c) by Allin Cottrell
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

/* interact.c for gretl */

#include "libgretl.h"
#include "var.h"
#include "gretl_func.h"
#include "loop_private.h"
#include "compat.h"
#include "system.h"
#include "forecast.h"
#include "cmd_private.h"
#include "libset.h"
#include "usermat.h"

/* equipment for the "shell" command */
#ifndef WIN32
# include <sys/wait.h>
# include <signal.h>
# include <errno.h>
# include <unistd.h>
# ifdef HAVE_PATHS_H
#  include <paths.h>
# endif
#endif

#define CMD_DEBUG 0

#include "laginfo.c"

typedef struct {
    int v;
    char vname[VNAMELEN];
    int lmin;
    int lmax;
    int *laglist;
} LAGVAR;

static void get_optional_filename (const char *line, CMD *cmd);

static int trydatafile (char *line, CMD *cmd)
{
    int i, m, n = strlen(line);
    char datfile[MAXLEN];

    *datfile = '\0';

    for (i=0; i<n; i++) {
	if ((n - i) > 4 && strncmp(line+i, "DATA", 4) == 0) {
	    sscanf(line + i, "%s", datfile);
	    m = strlen(datfile);
	    if (datfile[m-1] == ',') {
		datfile[m-1] = '\0';
	    }
	    lower(datfile);
	    i += 4;
	} else if (line[i] == '*' && line[i+1] == ')') {
	    cmd->ignore = 0;
	}
    }

    if (*datfile) {
	sprintf(line, "open %s.gdt", datfile);
	return 1;
    } 

    return 0;
}

static int filter_comments (char *line, CMD *cmd)
{
    int i, j = 0, n = strlen(line);
    char tmpstr[MAXLEN], datfile[MAXLEN];

    if (n >= MAXLEN) {
	return 0;
    }
    
    for (i=0; i<n; i++) {
	if (!strncmp(line + i, "(*", 2) || !strncmp(line + i, "/*", 2)) {
	    cmd->ignore = 1;
	    if (line[i+3] == '!') { /* special code for data file to open */
		sscanf(line + 4, "%s", datfile);
		sprintf(line, "open %s", datfile);
		cmd->ignore = 0;  /* FIXME ? */
		return 0;
	    }
	} else if (!strncmp(line + i, "*)", 2) || !strncmp(line + i, "*/", 2)) {
	    cmd->ignore = 0; 
	    i += 2;
	    while (isspace((unsigned char) line[i]) && i < n) {
		i++;
	    }
	}

	if (!cmd->ignore && line[i] != '\r') {
	    tmpstr[j] = line[i];
	    j++;
	}
    }

    tmpstr[j] = '\0';
    strcpy(line, tmpstr);

    return (*line == '\0');
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

/* catch aliased command words and assign ci; return 1
   if alias caught, else 0. */

static int catch_command_alias (char *line, CMD *cmd)
{
    char *s = cmd->word;

    cmd->ci = 0;

    if (!strcmp(line, "q")) {
	strcpy(s, "quit");
	cmd->ci = QUIT;
    } if (!strcmp(line, "x")) {
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
    } else if (!strcmp(s, "label")) {
	cmd->ci = SETINFO;
    } else if (!strcmp(line, "smpl full")) {
	strcpy(line, "smpl");
	cmd->opt = OPT_F;
    } else if (!strcmp(s, "sample")) {
	cmd->ci = SMPL;
    } else if (!strcmp(s, "list")) {
	cmd->ci = REMEMBER;
	cmd->opt = OPT_L;
    } else if (*s == '!') {
	cmd->ci = SHELL;
    } else if (!strcmp(s, "funcerr")) {
	cmd->ci = FUNCERR;
    } else if (!strcmp(line, "end if")) {
	strcpy(s, "endif");
	strcpy(line, "endif");
	cmd->ci = ENDIF;
    }

    return cmd->ci;
}

#define REQUIRES_PARAM(c) (c == ADDTO || \
                           c == FCAST || \
                           c == FUNC || \
                           c == LOOP ||  \
                           c == MULTIPLY || \
                           c == NULLDATA || \
                           c == OMITFROM || \
                           c == REMEMBER || \
                           c == SETMISS)

#define REQUIRES_ORDER(c) (c == ADF || \
                           c == ADDOBS || \
                           c == ARCH || \
                           c == COINT || \
                           c == COINT2 || \
                           c == KPSS || \
                           c == VAR || \
                           c == VECM)

#define NO_VARLIST(c) (c == ADDOBS || \
                       c == APPEND || \
                       c == BREAK || \
                       c == CHOW || \
	               c == CRITERIA || \
	               c == CRITICAL || \
	               c == CUSUM || \
                       c == DATA || \
                       c == END || \
	               c == ENDLOOP || \
                       c == ESTIMATE || \
	               c == EQNPRINT || \
	               c == FCAST || \
	               c == FCASTERR || \
	               c == FIT || \
                       c == FNCALL || \
                       c == FUNC || \
                       c == FUNCERR || \
	               c == GENR || \
	               c == HAUSMAN || \
                       c == HELP || \
	               c == IMPORT || \
                       c == INCLUDE || \
    	               c == INFO || \
 	               c == LABELS || \
                       c == LEVERAGE || \
                       c == LMTEST || \
                       c == LOOP || \
                       c == MATRIX || \
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
                       c == SYSTEM || \
                       c == TABPRINT || \
                       c == TESTUHAT || \
                       c == TRANSPOSE || \
                       c == VARLIST || \
                       c == VIF)

#define USES_LISTSEP(c) (c == AR || \
                         c == ARMA || \
                         c == EQUATION || \
                         c == GARCH || \
                         c == POISSON || \
                         c == SCATTERS || \
                         c == TSLS || \
                         c == XTAB)

#define NEEDS_LISTSEP(c) (c == AR || \
                          c == ARMA || \
                          c == GARCH || \
                          c == TSLS || \
                          c == XTAB)

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
                                  c == STORE || \
                                  c == SUMMARY)

#define SCALARS_OK_IN_LIST(c) (c == DELEET || \
                               c == PRINT || \
                               c == STORE)

#define RETURNS_LIST(c) (c == DIFF || \
                         c == DUMMIFY || \
                         c == LDIFF || \
                         c == SDIFF || \
                         c == LAGS || \
                         c == LOGS || \
                         c == SQUARE)

static int flow_control (const char *line, double ***pZ, 
			 DATAINFO *pdinfo, CMD *cmd)
{
    int ci = cmd->ci;
    int err = 0;

    /* clear to proceed? */
    if (!ifstate(IS_FALSE) && 
	ci != IF && ci != ELSE && ci != ENDIF) {
	return 0;
    }

    if (ci == IF) {
	int ok = if_eval(line, pZ, pdinfo);

	if (ok == -1) {
	    err = 1;
	} else if (ok) {
	    err = ifstate(SET_TRUE);
	} else {
	    err = ifstate(SET_FALSE);
	}
    } else if (ci == ELSE) {
	err = ifstate(SET_ELSE);
    } else if (ci == ENDIF) {
	err = ifstate(SET_ENDIF);
    }

    if (err) {
	cmd->errcode = E_SYNTAX;
    }    

    return 1;
}

static char cmd_savename[MAXSAVENAME];

static void maybe_extract_savename (char *s, CMD *cmd)
{
    *cmd->savename = 0;

    if (strncmp(s, "genr ", 5) && strstr(s, " <- ")) {
	int n, len, quote;

	quote = (*s == '"');
	len = strcspn(s, "<");
	if (len < 2) {
	    return;
	}
	n = len - 1 - quote;
	if (n > MAXSAVENAME - 1) {
	    n = MAXSAVENAME - 1;
	}
	strncat(cmd->savename, s + quote, n);
	if (cmd->savename[n-1] == '"') {
	    cmd->savename[n-1] = 0;
	}
	strcpy(cmd_savename, cmd->savename);
	shift_string_left(s, len + 3);
    }
}

static int 
get_maybe_quoted_storename (CMD *cmd, char *s, int *nf)
{
    int quoted = 0;
    int q, len;

    while (*s == ' ') s++;

    q = *s;

    if (q == '"' || q == '\'') {
	char *p = strchr(s + 1, q);

	if (p == NULL) {
	    return E_SYNTAX;
	}
	len = p - s - 1;
	if (len == 0) {
	    return E_SYNTAX;
	}
	quoted = 1;
    } else {
	len = strcspn(s, " ");
    }

    free(cmd->param);
    cmd->param = gretl_strndup(s + quoted, len);
    if (cmd->param == NULL) {
	return E_ALLOC;
    }

    if (quoted) {
	char *p = cmd->param;

	while (*p) {
	    if (*p == ' ') *nf -= 1;
	    p++;
	}
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

#define LAG_DEBUG 0

static int get_contiguous_lags (LAGVAR *lv,
				const char *l1str, const char *l2str,
				const double **Z, const DATAINFO *pdinfo)
{
    const char *p;
    int i, v, lag, lsign;
    int err = 0;

    for (i=0; i<2 && !err; i++) {
	p = (i == 0)? l1str : l2str;

	if (*p == '0') {
	    lsign = 1;
	} else if (*p == '-') {
	    lsign = 1;
	    p++;
	} else if (*p == '+') {
	    lsign = -1;
	    p++;
	} else {
	    err = 1;
	    break;
	}

	if (isdigit(*p)) {
	    lag = atoi(p);
	} else {
	    v = varindex(pdinfo, p);
	    if (v < pdinfo->v) {
		lag = Z[v][0];
	    } else {
		err = 1;
	    }
	}

	if (!err) {
	    if (i == 0) {
		lv->lmin = lsign * lag;
	    } else {
		lv->lmax = lsign * lag;
	    }
	}
    }

    return err;
}

static int parse_lagvar (const char *s, LAGVAR *lv, 
			 const double **Z, const DATAINFO *pdinfo)
{
    char l1str[16], l2str[16];
    int lag, i, err = 1;

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
    } else if (sscanf(s, "%15[^(](%d)", lv->vname, &lag) == 2) {
	lv->lmin = lv->lmax = -lag;
	err = 0;
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
    int oldn = cmd->list[0];
    int *list;

    list = realloc(cmd->list, (oldn + add) * sizeof *list);
    if (list == NULL) {
	cmd->errcode = E_ALLOC;
	strcpy (gretl_errmsg, 
		_("Memory allocation failed for command list"));
	return 1;
    }

    /* one of the vars was "already assumed" */
    list[0] += (add - 1);
    
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

int auto_lag_ok (const char *s, int *lnum,
		 double ***pZ, DATAINFO *pdinfo,
		 CMD *cmd)
{
    LAGVAR lagvar;
    int nlags, i;
    int llen = *lnum;
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
	cmd->errcode = E_PARSE;
	ok = 0;
	goto bailout;
    }

    if (nlags > 1 && expand_command_list(cmd, nlags)) {
	ok = 0;
	goto bailout;
    }

    for (i=0; i<nlags && ok; i++) {
	int laglen, vnum;

	if (lagvar.laglist != NULL) {
	    laglen = lagvar.laglist[i+1];
	} else {
	    laglen = lagvar.lmin + i * lincr;
	}

	vnum = laggenr(lagvar.v, laglen, pZ, pdinfo);

#if LAG_DEBUG
	fprintf(stderr, "laggenr for var %d (%s), lag %d, gave vnum = %d\n",
		lagvar.v, pdinfo->varname[lagvar.v], laglen, vnum);
#endif
	if (vnum < 0) {
	    cmd->errcode = 1;
	    sprintf(gretl_errmsg, _("generation of lag variable failed"));
	    ok = 0;
	} else {
	    int err;

	    /* record info regarding the auto-generation of lags,
	       so that we'll be able to echo the command properly --
	       see the echo_cmd() function. 
	    */
	    cmd->list[llen++] = vnum;
	    err = add_to_list_lag_info(lagvar.v, laglen, vnum, cmd);
	    if (err) {
		cmd->errcode = E_ALLOC;
		ok = 0;
	    }
	}
    }

    if (ok) {
	*lnum = llen;
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
	char ostr[8] = {0};
	char word[32] = {0};
	int v;

	sscanf(s, "%7[^ ,;]", ostr);
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

static int auto_transform_ok (const char *s, int *lnum,
			      double ***pZ, DATAINFO *pdinfo,
			      CMD *cmd)
{
    char fword[9];
    int *genlist = NULL;
    int trans = 0;
    int order = 0;
    gretlopt opt = OPT_NONE;
    int err = 0, ok = 1;

    if (sscanf(s, "%8[^(](", fword)) {
	char *param = NULL;
	int *gotlist;
	int vnum = 0;

	if (!strcmp(fword, "cross")) {
	    strcpy(fword, "square");
	    opt = OPT_O;
	}

	trans = gretl_command_number(fword);
	if (!RETURNS_LIST(trans)) {
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
	cmd->errcode = E_PARSE;
	return 0;
    }	

    if (trans == LOGS) {
	err = list_loggenr(genlist, pZ, pdinfo);
    } else if (trans == DIFF || trans == LDIFF || trans == SDIFF) {
	err = list_diffgenr(genlist, trans, pZ, pdinfo);
    } else if (trans == SQUARE) {
	err = list_xpxgenr(&genlist, pZ, pdinfo, opt);
    } else if (trans == LAGS) {
	err = list_laggenr(&genlist, order, pZ, pdinfo);
    } else if (trans == DUMMIFY) {
	err = list_dumgenr(&genlist, pZ, pdinfo);
    }

    if (err) {
	cmd->errcode = err;
	ok = 0;
    } else {
	cmd->list[0] -= 1;
	gretl_list_insert_list(&cmd->list, genlist, *lnum);
	*lnum = cmd->list[0];
    }

    free(genlist);

    return ok;
} 

static int plot_var_ok (const char *s, int *lnum,
			double ***pZ, DATAINFO *pdinfo,
			CMD *cmd)
{
    int pnum, ok = 1;
    
    if (strcmp(s, "qtrs") &&
	strcmp(s, "months") &&
	strcmp(s, "time")) { 
	return 0;
    }

    pnum = plotvar_from_varname(pZ, pdinfo, s);

    if (pnum < 0) {
	cmd->errcode = 1;
	ok = 0;
	sprintf(gretl_errmsg, 
		_("Failed to add plotting index variable"));
    } else {
	cmd->list[*lnum] = pnum;
	*lnum += 1;
    }

    return ok;
} 

#if defined(USE_GTK2) || defined (HAVE_FNMATCH_H)

static int wildcard_expand_ok (const char *s, int *lnum,
			       const DATAINFO *pdinfo, CMD *cmd)
{
    int ok = 0;

    if (strchr(s, '*') != NULL) {
	int *wildlist = varname_match_list(pdinfo, s);

	if (wildlist != NULL) {
	    int k, nw = wildlist[0];
	    int llen = *lnum;

	    if (expand_command_list(cmd, nw)) {
		return 0;
	    }
	    for (k=1; k<=nw; k++) {
		cmd->list[llen++] = wildlist[k];
	    }
	    free(wildlist);
	    *lnum = llen;
	    ok = 1;
	}
    }

    return ok;
}

#endif 

static void parse_rename_cmd (const char *line, CMD *cmd, 
			      const DATAINFO *pdinfo)
{
    int vtest, vtarg;
    char targ[VNAMELEN];
    char newname[VNAMELEN];
    char numstr[8];

    line += strlen(cmd->word);

    if (sscanf(line, "%15s %15s", targ, newname) != 2) {
	cmd->errcode = E_DATA;
	sprintf(gretl_errmsg, "rename: %s", 
		_("requires a variable number and a new name"));
	return;
    }

    if (isdigit(*targ)) {
	vtarg = atoi(targ);
    } else {
	/* we're given the name of a variable? */
	vtarg = varindex(pdinfo, targ);
    }

    if (vtarg >= pdinfo->v || vtarg < 1) {
	cmd->errcode = E_DATA;
	sprintf(gretl_errmsg, _("Variable number %d is out of bounds"), vtarg);
	return;
    } 

    vtest = varindex(pdinfo, newname);
    if (vtest < pdinfo->v && vtest != vtarg) {
	sprintf(gretl_errmsg, _("'%s': there is already a variable "
				"of this name"), newname);
	cmd->errcode = E_DATA;
	return;
    }

    if (check_varname(newname)) {
	cmd->errcode = E_DATA;
	return;
    }

    free(cmd->param);
    cmd->param = gretl_strdup(newname);
    if (cmd->param == NULL) {
	cmd->errcode = E_ALLOC;
	return;
    }

    sprintf(numstr, "%d", vtarg);

    free(cmd->extra);
    cmd->extra = gretl_strdup(numstr);
}

static void parse_outfile_cmd (char *s, CMD *cmd)
{
    s += strlen(cmd->word);

    while (isspace((unsigned char) *s) || *s == '"') {
	s++;
    }

    if (*s) {
	free(cmd->param);
	cmd->param = gretl_strdup(s);
	if (cmd->param == NULL) {
	    cmd->errcode = E_ALLOC;
	} else {
	    int n;

	    tailstrip(cmd->param);
	    n = strlen(cmd->param);
	    if (cmd->param[n] == '"') {
		cmd->param[n] = 0;
	    }
	}
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
		cmd->errcode = E_ALLOC;
	    } else {
		sprintf(cmd->param, "ymax=%s", numstr);
	    }
	    *p = '\0';
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

    while (*s == ' ') s++;
    len = get_field_length(s);

    if (len >= 0) {
	strncat(field, s, len);
    } else {
	err = 1;
    }

#if CMD_DEBUG
    fprintf(stderr, "field_from_line, got '%s'\n", field);
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

static int plausible_genr_start (const char *s, CMD *cmd, 
				 const DATAINFO *pdinfo)
{
    if (strchr(s, '=') || strstr(s, "++") || strstr(s, "--")) {
	const char *ok = "+-*/=[";
	char word[VNAMELEN];

	if (sscanf(s, "%15[^[ +-*/=]", word)) {
	    s += strlen(word);
	    while (*s == ' ') s++;
	    if (strspn(s, ok) && check_varname(word) == 0) {
		cmd->ci = GENR;
		if (get_matrix_by_name(word, pdinfo)) {
		    cmd->opt = OPT_M;
		}
	    }
	}
    }

    return cmd->ci;
}

/* if we find a semicolon directly after a varname, insert a space so
   that we can count the fields in the line correctly */

static int fix_semicolon_after_var (char *s)
{
    int len = strlen(s);
    int i, j;

    for (i=0; i<len-1; i++) {
	if (isalnum((unsigned char) s[i]) && s[i+1] == ';') {
	    if (len < MAXLINE - 1) {
		for (j=len; j>i+1; j--) {
		    s[j] = s[j-1];
		}
		s[i+1] = ' ';
		s[len + 1] = '\0';
		len++;
	    } else {
		break;
	    }
	}
    }

    return len;
}

/* apparatus for checking that the "end" command is valid */

#define COMMAND_CAN_END(c) (c == FUNC || \
                            c == MLE || \
                            c == NLS || \
			    c == RESTRICT || \
			    c == SYSTEM)

static int check_end_command (CMD *cmd)
{
    if (*cmd->param) {
	int cmdcode = gretl_command_number(cmd->param);

	if (cmdcode == LOOP) {
	    cmd->ci = ENDLOOP;
	} else if (!COMMAND_CAN_END(cmdcode)) {
	    cmd->errcode = 1;
	    sprintf(gretl_errmsg, _("command 'end %s' not recognized"), 
		    cmd->param);
	}
    } else {
	cmd->errcode = 1;
	strcpy(gretl_errmsg, _("end: nothing to end")); 
    }

    return cmd->errcode;
}

static int resize_cmd_param (CMD *cmd, const char *s, int inlen)
{
    char *param;
    int len;

    if (inlen > 0) {
	len = inlen;
    } else if (s != NULL) {
	if (strchr(s, ' ') == NULL) {
	    len = strlen(s);
	} else {
	    len = strcspn(s, " ");
	}
	len++;
    } else {
	return 1;
    }

    param = realloc(cmd->param, len);
    if (param == NULL) {
	return 1;
    }

    cmd->param = param;
    
    return 0;
}

/* Capture the next 'word' found following the initial command word
   (or the whole remainder of the line in some cases) as the parameter
   for cmd.  Flag an error if the command requires a parameter but
   none is found.
*/

static int capture_param (const char *s, CMD *cmd,
			  const double **Z, 
			  const DATAINFO *pdinfo)
{
    /* if param has already been written by some special
       routine, don't overwrite it */
    if (*cmd->param != '\0') {
	return cmd->errcode;
    }

    /* skip past leading word on line */
    s += strcspn(s, " ");
    s += strspn(s, " ");

    if (string_is_blank(s)) {
	if (REQUIRES_PARAM(cmd->ci) || REQUIRES_ORDER(cmd->ci)) {
	    cmd->errcode = E_PARSE;
	    sprintf(gretl_errmsg, _("%s: required parameter is missing"),
		    cmd->word);
	}
    } else {
	if (resize_cmd_param(cmd, NULL, strlen(s) + 1)) {
	    cmd->errcode = E_ALLOC;
	} else if (cmd->ci == PRINT || cmd->ci == FUNCERR) {
	    /* grab the whole remainder of line */
	    strcpy(cmd->param, s);
	} else {
	    /* grab one 'word' */
	    sscanf(s, "%s", cmd->param);
	}
#if CMD_DEBUG
	fprintf(stderr, "capture_param: s='%s', param='%s'\n",
		s, cmd->param);
#endif
	if (REQUIRES_ORDER(cmd->ci)) {
	    cmd->order = gretl_int_from_string(cmd->param,
					       Z, pdinfo, 
					       &cmd->errcode);
	}
    }

    if (cmd->ci == END) {
	/* test that param is present and valid */
	check_end_command(cmd);
    }

    return cmd->errcode;
}

static int gretl_cmd_clear (CMD *cmd)
{
    cmd->ci = 0;
    cmd->errcode = 0;
    cmd->nolist = 0;
    *cmd->word = '\0';

    if (cmd->list == NULL || cmd->param == NULL || cmd->extra == NULL) {
	cmd->errcode = E_ALLOC;
    } else {
	cmd->list[0] = 0;
	*cmd->param = '\0';
	*cmd->extra = '\0';
    }

    cmd_lag_info_destroy(cmd);

    return cmd->errcode;
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
	cmd->errcode = E_ALLOC;
	strcpy (gretl_errmsg, _("Memory allocation failed for command list"));
    } else {
	list[0] = nf;
	for (i=1; i<=nf; i++) {
	    list[i] = 0;
	}
	cmd->list = list;
    }

    return cmd->errcode;
}

static int maybe_print_object (const char *line, int nf, CMD *cmd)
{
    int ret = 0;

    if (cmd->opt & OPT_L) {
	char word[8];
	int *list;

	if (nf == 0 || (sscanf(line, "%5s", word) && !strcmp(word, "print"))) {
	    list = get_list_by_name(cmd->param);
	    if (list != NULL) {
		free(cmd->list);
		cmd->list = gretl_list_copy(list);
		if (cmd->list != NULL) {
		    cmd->opt |= OPT_P;
		    ret = 1;
		}
	    }
	}
    }

    return ret;
}

static int trap_comments (char *s, CMD *cmd)
{
    int ret = 0;

    if (*s == '#' || (*s == '/' && *(s+1) == '/')) {	
	cmd->nolist = 1;
	cmd->ci = CMD_COMMENT;
	ret = 1;
    } else if (strstr(s, " #") || strstr(s, "//")) {
	int quoted = 0;
	int braced = 0;

	while (*s) {
	    if (*s == '"') {
		quoted = !quoted;
	    } else if (!quoted) {
		if (*s == '{') {
		    braced++;
		} else if (*s == '}') {
		    braced--;
		}
	    }
	    if (!quoted && !braced) {
		if ((*s == ' ' && *(s+1) == '#') ||
		    (*s == '/' && *(s+1) == '/')) {
		    *s = '\0';
		    break;
		}
	    }
	    s++;
	}
    }

    return ret;
}

static int creating_null_list (CMD *cmd, int nf, const char *s)
{
    int ret = 0;

    if (nf == 2 && (cmd->opt & OPT_L)) {
	if (!strcmp(s, " = null")) {
	    cmd->list[0] = 0;
	    ret = 1;
	}
    }

    return ret;
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

    while (*s) {
	if (*s == ';') c++;
	s++;
    }

    return c;
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
    int j, nf, linelen, pos, v, lnum;
    int gotdata = 0;
    int sepcount = 0;
    int read_lags = 0;
    char *remainder = NULL;
    char field[FIELDLEN] = {0};

    if (gretl_cmd_clear(cmd)) {
	return cmd->errcode;
    }

    *gretl_errmsg = '\0';

    compress_spaces(line);

#if CMD_DEBUG
    fprintf(stderr, "parsing '%s'\n", line);
#endif

    /* legacy: look for Ramanathan practice files */
    if (line[0] == '(' && line[1] == '*') {
	cmd->ignore = 1;
	gotdata = trydatafile(line, cmd);
    }

    /* trap old-style comments */
    if (!gotdata) {
	if (filter_comments(line, cmd)) {
	    cmd->nolist = 1;
	    cmd->ci = CMD_COMMENT;
	    return cmd->errcode;
	}
    }

    /* also new-style comments */
    if (trap_comments(line, cmd)) {
	return cmd->errcode;
    } 

    /* extract any options */
    cmd->opt = get_gretl_options(line, &cmd->errcode);
    if (cmd->errcode) {
	return cmd->errcode;
    }    

    /* extract "savename" for storing an object? */
    maybe_extract_savename(line, cmd);

    /* no command here? */
    if (sscanf(line, "%8s", cmd->word) != 1) {
	cmd->nolist = 1;
	cmd->ci = CMD_NULL;
	return cmd->errcode;
    }

#if CMD_DEBUG
    fprintf(stderr, "cmd->word = '%s'\n", cmd->word);
#endif

    /* backwards compatibility */
    accommodate_obsolete_commands(line, cmd);

    /* replace simple aliases and a few specials */
    catch_command_alias(line, cmd);

    /* "remember": capture the line in cmd's "extra" field */
    if (cmd->ci == REMEMBER) {
	free(cmd->extra);
	cmd->extra = gretl_strdup(line);
#if CMD_DEBUG
	fprintf(stderr, "cmd->extra = '%s'\n", cmd->extra);
#endif
    }

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
	    /* call of user-defined function? */
	    if (gretl_get_user_function(line, &cmd->param)) {
		cmd->ci = FNCALL;
	    } else if (!plausible_genr_start(line, cmd, pdinfo)) {
		cmd->errcode = 1;
		sprintf(gretl_errmsg, _("command '%s' not recognized"), 
			cmd->word);
		goto bailout;
	    }
	}
    }

#if CMD_DEBUG
    fprintf(stderr, "cmd->ci = %d\n", cmd->ci);
#endif

    /* if, else, endif controls: should this come earlier? */
    if (flow_control(line, pZ, pdinfo, cmd)) {
	cmd->nolist = 1;
	cmd->ci = CMD_NULL;
	return cmd->errcode;
    }

    /* tex printing commands can take a filename parameter */
    if (cmd->ci == EQNPRINT || cmd->ci == TABPRINT) {
	get_optional_filename(line, cmd);
    } 

    /* the "outfile" command may have a filename */
    else if (cmd->ci == OUTFILE) {
	parse_outfile_cmd(line, cmd);
    }

    /* the "rename" command calls for a variable number and a
       new name */
    else if (cmd->ci == RENAME) {
	parse_rename_cmd(line, cmd, pdinfo);
    }  

    /* commands that never take a list of variables */
    if (NO_VARLIST(cmd->ci)) { 
	cmd->nolist = 1;
	capture_param(line, cmd, (const double **) *pZ, pdinfo);
	return cmd->errcode;
    }

    /** now for a few command which may or may not take a list **/

    /* PRINT can take a list, but not in its string literal variant */
    if (cmd->ci == PRINT && strstr(line, "\"")) {
	cmd->nolist = 1;
	capture_param(line, cmd, NULL, NULL);
	return cmd->errcode;
    }

    /* SMPL can take a list, but only in case of OPT_M
       "--no-missing" */
    if (cmd->ci == SMPL && !(cmd->opt & OPT_M)) {
	cmd->nolist = 1;
	return cmd->errcode;
    }	

    /* boxplots take a list, but if there are Boolean conditions
       embedded, the line has to be parsed specially */
    if (cmd->ci == BXPLOT && strchr(line, '(')) {
	cmd->nolist = 1;
	return cmd->errcode;
    }

    /* OMIT typically takes a list, but can be given without args
       to omit last var */
    if (cmd->ci == OMIT && string_is_blank(line + 4)) {
	cmd->nolist = 1;
	return cmd->errcode;
    }

    /** OK, now we're definitely doing a list-oriented command,
	We begin by taking care of a few specials **/

    /* GNUPLOT can have a block of stuff to pass literally
       to gnuplot */
    if (cmd->ci == GNUPLOT) {
	grab_gnuplot_literal_block(line, cmd);
    }

    /* "logistic" can have a "ymax" parameter */
    else if (cmd->ci == LOGISTIC) {
	parse_logistic_ymax(line, cmd);
    } 

    /* fix lines that contain a semicolon right after a var */
    linelen = fix_semicolon_after_var(line);

    /* find number of space-separated fields remaining in line,
       record our reading position, and make a copy of the
       remainder of the line
    */
    nf = count_free_fields(line) - 1;
    pos = strlen(cmd->word);
    remainder = gretl_strdup(line + pos + 1);
    if (remainder == NULL) {
	cmd->errcode = E_ALLOC;
	goto bailout;
    }

#if CMD_DEBUG
    fprintf(stderr, "nf=%d, remainder='%s'\n", nf, remainder);
#endif

    /* need to treat rhodiff specially -- put everything from
       the end of the command word to the first semicolon into
       "param", for processing later */
    if (cmd->ci == RHODIFF) { 
	if (!get_rhodiff_or_lags_param(remainder, cmd)) {
	    cmd->errcode = E_SYNTAX;
	    goto bailout;
	}
	strcpy(line, remainder);
	linelen = strlen(line);
	nf = count_fields(line);
	pos = 0;
    }

    if (cmd->ci == LAGS) {
	/* optional initial lags field */
	if (get_rhodiff_or_lags_param(remainder, cmd)) {
	    strcpy(line, remainder);
	    linelen = strlen(line);
	    nf = count_fields(line);
	    pos = 0;
	} else {
	    *remainder = '\0';
	}
    }	

    /* "store" is a special case since the filename that comes
       before the list may be quoted, and have spaces in it.  Groan */
    if (cmd->ci == STORE && nf > 0) {
	cmd->errcode = get_maybe_quoted_storename(cmd, remainder, &nf);
	if (cmd->errcode) {
	    goto bailout;
	} else {
	    pos = 0;
	    if (--nf > 0) {
		strcpy(line, remainder);	
		linelen = strlen(line);
	    }		
	}
    }

    /* 
       "multiply" takes a multiplier;
       "omitfrom" and "addto" take the ID of a previous model;
       "setmiss" takes a value to be interpreted as missing;
       "remember" takes a thing to remember;
       these are captured in cmd->param
    */
    if (REQUIRES_ORDER(cmd->ci) ||
	cmd->ci == ADDTO ||
	cmd->ci == OMITFROM ||
	cmd->ci == MULTIPLY ||
	cmd->ci == REMEMBER ||
	cmd->ci == SETMISS) {
	capture_param(line, cmd, (const double **) *pZ, pdinfo);
	if (cmd->errcode) {
	    goto bailout;
	} else {
	    strcpy(remainder, line + pos + 1 + strlen(cmd->param));
	    pos = 0;
	    if (--nf > 0) {
		strcpy(line, remainder);
		linelen = strlen(line);
	    }
	} 
    }

    if (cmd->ci == REMEMBER) {
	if (nf <= 1 && maybe_print_object(line, nf, cmd)) {
	    return cmd->errcode;
	} else if (nf == 0) {
	    /* creating empty object implicitly, OK */
	    return cmd->errcode;
	} else if (creating_null_list(cmd, nf, remainder)) {
	    /* creating empty object with explicit "null" */
	    return cmd->errcode;
	} 
	line += strspn(line, " ");
	if (*line != '=') {
	    cmd->errcode = E_PARSE;
	    return cmd->errcode;
	}
	line += strspn(line, "=");
	nf--;
	pos = 0;
	linelen = strlen(line);
    } else if (cmd->ci == MULTIPLY || cmd->ci == VECM) { 
	char extra[4];

	sscanf(line, "%3s", extra);
	free(cmd->extra);
	cmd->extra = gretl_strdup(extra);
	shift_string_left(line, strlen(extra));
	nf--;
	pos = 0;
	linelen = strlen(line);
    } 

    /* get a count of ';' separators in line */
    sepcount = get_sepcount(line);

    if (NEEDS_LISTSEP(cmd->ci) && sepcount == 0) {
	/* missing field in command */
	cmd->errcode = E_ARGS;
	free(remainder);
	return cmd->errcode;
    }

    if (cmd->ci == AR || cmd->ci == ARMA || cmd->ci == GARCH) {
	/* flag acceptance of lag orders in list */
	read_lags = 1;
    }

    /* allocate space for the command list */
    if (resize_command_list(cmd, nf)) {
	goto bailout;
    }

    /* now assemble the command list */

    for (j=1, lnum=1; j<=nf; j++) {

	strcpy(remainder, line + pos + 1);

	/* special: optional lag order for correlogram */
	if (cmd->ci == CORRGM && j == 2) {
	    cmd->list[0] = 1;
	    if (resize_cmd_param(cmd, remainder, 0)) {
		cmd->errcode = E_ALLOC;
		goto bailout;
	    }
	    sscanf(remainder, "%s", cmd->param);
	    break;
	}

	cmd->errcode = get_next_field(field, remainder);
	if (cmd->errcode) {
	    goto bailout;
	}

	if (isalpha((unsigned char) *field)) {
	    /* probably should be the name of a variable */
	    int *savedlist;

	    if (field[strlen(field) - 1] == ';') {
		/* strip any trailing semicolon */
		field[strlen(field) - 1] = '\0';
	    }

	    if (read_lags) {
		int k = gretl_int_from_string(field, (const double **) *pZ,
					      pdinfo, &cmd->errcode);
		if (!cmd->errcode) {
		    cmd->list[lnum++] = k;
		}
	    } else if ((v = varindex(pdinfo, field)) < pdinfo->v) {
		/* yes, it's an existing variable */
		cmd->list[lnum++] = v;
	    } else if ((savedlist = get_list_by_name(field)) != NULL) {
		/* or it's a pre-defined list */
		cmd->list[0] -= 1;
		cmd->errcode = gretl_list_insert_list(&cmd->list, savedlist,
						      lnum);
		lnum += savedlist[0];
	    } else { 
		/* possibly an auto-generated variable? */
		if (strchr(field, '(') != NULL) {
		    /* Case 1: automated lags: e.g. 'var(-1)' */
		    if (auto_lag_ok(field, &lnum, pZ, pdinfo, cmd)) {
			/* handled, get on with it */
			pos += strlen(field) + 1;
			continue; 
		    }
		    /* Case 2: automated transformations: e.g. 'logs(list)' */
		    else if (auto_transform_ok(field, &lnum, pZ, pdinfo, cmd)) {
			/* handled, get on with it */
			pos += strlen(field) + 1;
			continue; 
		    }
		}

		/* Case 3: special plotting variable */
		else if (!cmd->errcode && 
			 plot_var_ok(field, &lnum, pZ, pdinfo, cmd)) {
		    /* handled, get on with it */
		    pos += strlen(field) + 1;
		    continue; 
		} 

#if defined(USE_GTK2) || defined (HAVE_FNMATCH_H)
		/* wildcard expansion? */
		else if (!cmd->errcode && 
			 wildcard_expand_ok(field, &lnum, pdinfo, cmd)) {
		    /* handled, get on with it */
		    pos += strlen(field) + 1;
		    continue; 			
		}
#endif
		/* try abbreviating the varname? */
		else if (!cmd->errcode) {
		    cmd->errcode = 1; /* presume guilt at this stage */
		    if (strlen(field) > 8) {
			char test[9];

			*test = 0;
			strncat(test, field, 8);
			if ((v = varindex(pdinfo, test)) <= pdinfo->v - 1) {
			    cmd->list[lnum++] = v;
			    cmd->errcode = 0;
			} 
		    } 
		    if (cmd->errcode) {
			sprintf(gretl_errmsg, 
				_("'%s' is not the name of a variable"), field);
		    }
		}

		if (cmd->errcode) {
		    goto bailout;
		}
	    }
	} /* end if isalpha(*field) */

#if defined(USE_GTK2) || defined (HAVE_FNMATCH_H)
	else if (*field == '*') {
	    if (wildcard_expand_ok(field, &lnum, pdinfo, cmd)) {
		/* handled, get on with it */
		pos += strlen(field) + 1;
		continue; 			
	    }
	}
#endif

	else if (isdigit(*field)) {
	    /* could be the ID number of a variable */
	    v = atoi(field);
	    if (!read_lags && v > pdinfo->v - 1) {
		cmd->errcode = 1;
		sprintf(gretl_errmsg, 
                       _("%d is not a valid variable number"), v);
		goto bailout;
	    }
	    cmd->list[lnum++] = v;
	}

	else if (*field == ';') {
	    /* could be the separator between two sub-lists */
	    if (USES_LISTSEP(cmd->ci)) {
		pos += strlen(field) + 1;
		cmd->list[lnum++] = LISTSEP;
		sepcount--;
		if (read_lags && sepcount == 0) {
		    /* turn off acceptance of AR lags etc */
		    read_lags = 0;
		}
		continue;
	    } else if (cmd->ci == VAR || cmd->ci == VECM) {
		pos += strlen(field) + 1;
		cmd->list[lnum++] = LISTSEP;
		continue;
	    } else {
		cmd->list[0] -= 1;
		break;
	    }
	}

	if (!isalpha((unsigned char) *field) && !isdigit((unsigned char) *field)) {
	    cmd->errcode = 1;
	    sprintf(gretl_errmsg, _("field '%s' in command is invalid"), field);
	    goto bailout;
	}

	/* check cmd->list for scalars */
	if (!read_lags && !(SCALARS_OK_IN_LIST(cmd->ci))) {
	    if (var_is_scalar(pdinfo, cmd->list[lnum-1])) {
		cmd->errcode = 1;
		sprintf(gretl_errmsg, _("variable %s is a scalar"), field);
		goto bailout;
	    }
	}

	pos += strlen(field) + 1;
    } /* end of loop through fields in command line */

    /* By now we're looking at a command that takes a list,
       which either has been specified already or needs to
       be filled out automatically */

    /* commands that can take a specified list, but where if the
       list is null or just ";" we want to operate on all variables
    */    
    if (DEFAULTS_TO_FULL_LIST(cmd->ci)) {
	if (cmd->list[0] == 0) {
	    cmd_full_list(pdinfo, cmd);
	    /* suppress echo of the list -- may be too long */
	    cmd->nolist = 1;
	}
    } else if (cmd->ci != SETMISS && 
	       cmd->ci != DELEET &&
	       cmd->ci != REMEMBER) {
	/* the command needs a list but doesn't have one */
	if (cmd->list[0] == 0) {
	    cmd->errcode = E_ARGS;
	}
    }

    if (NEEDS_TWO_VARS(cmd->ci) && cmd->list[0] == 1) {
	cmd->errcode = E_ARGS;
    }

    /* check list for duplicated variables? */
    if (!cmd->errcode && !cmd->nolist) {
	int dupv = gretl_list_duplicates(cmd->list, cmd->ci);

	if (dupv >= 0) {
	    printlist(cmd->list, "cmd->list");
	    cmd->errcode = E_UNSPEC;
	    sprintf(gretl_errmsg, 
		    _("var number %d duplicated in the command list."),
		    dupv);
	}
    }

 bailout:

    /* double-check that allocation hasn't failed */
    if (cmd->list == NULL || cmd->param == NULL || cmd->extra == NULL) {
	cmd->errcode = E_ALLOC;
    }

#if CMD_DEBUG
    printlist(cmd->list, "cmd->list");
    fprintf(stderr, "cmd->errcode = %d, context=%d\n", cmd->errcode,
	    cmd->context);
#endif

    if (cmd->errcode) {
	cmd->context = 0;
    }

    free(remainder);

    return cmd->errcode;
}

/**
 * help:
 * @cmdword: the command on which help is wanted.
 * @helpfile: path to the gretl help file.
 * @prn: pointer to gretl printing struct.
 *
 * Searches in @helpfile for help on @cmdword and, if help is found,
 * prints it to @prn.  If @cmdword is %NULL, lists the valid
 * commands.
 *
 * Returns: 0 on success, 1 if the helpfile was not found or the
 * requested topic was not found.
 */

int help (const char *cmdword, const char *helpfile, PRN *prn)
{
    FILE *fp;
    char word[9];
    char line[128];
    int i, j, ok;

    if (cmdword == NULL || *cmdword == '\0') {
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

	return 0;
    }

    ok = gretl_command_number(cmdword) > 0;

    if (!ok) {
	if (gretl_is_public_user_function(cmdword)) {
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
		pputs(prn, line);
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
    int opt = 0;
    const char *s = argv[1];

    *fname = '\0';
    *force_lang = 0;

#ifdef ENABLE_NLS
    if (!strcmp(s, "-e") || !strncmp(s, "--english", 9)) { 
	*force_lang = ENGLISH;
    } else if (!strcmp(s, "-q") || !strncmp(s, "--basque", 8)) { 
	*force_lang = BASQUE;
    }
    if (*force_lang) {
	if (--argc < 2) {
	    return 0;
	}
	argv++;
	s = argv[1];
    }	
#endif

    if (strcmp(s, "-b") == 0 || strncmp(s, "--batch", 7) == 0) 
	opt = OPT_BATCH;
    else if (strcmp(s, "-h") == 0 || strcmp(s, "--help") == 0) 
	opt = OPT_HELP;
    else if (strcmp(s, "-p") == 0 || strcmp(s, "--pvalue") == 0) 
	opt = OPT_PVALS;
    else if (strcmp(s, "-v") == 0 || strcmp(s, "--version") == 0) 
	opt = OPT_VERSION;
    else if (strcmp(s, "-r") == 0 || strncmp(s, "--run", 5) == 0) 
	opt = OPT_RUNIT;
    else if (strcmp(s, "-d") == 0 || strncmp(s, "--db", 4) == 0) 
	opt = OPT_DBOPEN;
    else if (strcmp(s, "-w") == 0 || strncmp(s, "--webdb", 7) == 0) 
	opt = OPT_WEBDB;
    else if (strcmp(s, "-c") == 0 || strncmp(s, "--dump", 6) == 0) 
	opt = OPT_DUMP;

    if (opt != 0) {
	argv++;
	argc--;
    }

    if (argc >= 2) {
	strncat(fname, argv[1], MAXLEN - 1);
    }

    return opt;
}

#ifndef WIN32

int gretl_shell (const char *arg)
{
    int pid;
    void (*old1)(int);
    void (*old2)(int);
    char shellnam[40];
    const char *theshell, *namep; 

    if (!get_shell_ok()) {
	strcpy(gretl_errmsg, "The shell command is not activated.");
	return 1;
    }

    old1 = signal(SIGINT, SIG_IGN);
    old2 = signal(SIGQUIT, SIG_IGN);

    if ((pid = fork()) == 0) {
	for (pid = 3; pid < 20; pid++) {
	    close(pid);
	}

	signal(SIGINT, SIG_DFL);
	signal(SIGQUIT, SIG_DFL);

	theshell = getenv("SHELL");
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
	if (arg) {
	    execl(theshell, shellnam, "-c", arg, NULL);
	} else {
	    execl(theshell, shellnam, NULL);
	}
	perror(theshell);
	return 1;
    }

    if (pid > 0) {
	while (wait(NULL) != pid);
    }

    signal(SIGINT, old1);
    signal(SIGQUIT, old2);

    if (pid == -1) {
	perror(_("Try again later"));
    }

    return 0;
}

#endif /* ! WIN32 */

static void trim_to_length (char *s, int oklen)
{
    int i, n = strlen(s);

    if (n < oklen - 1) return;

    for (i=n-1; i>0; i--) {
	if (s[i] == ' ') {
	    s[i] = '\0';
	    break;
	}
    }
}

#define SAFELEN 78

static void 
real_safe_print_line (const char *line, int cli, int batch, 
		      int script, PRN *prn)
{
    char tmp[SAFELEN];
    const char *split = " \\";
    const char *leader;
    const char *leaders[] = { "? ", "> " };
    const char *p, *q;
    int n, out, rem;

    if (!cli && batch) return;

    if (gretl_compiling_loop()) {
	leader = leaders[1];
    } else {
	leader = leaders[0];
    }

    if (cli) {
	printf("%s", (batch)? leader : " ");
    } else if (script) {
	pputs(prn, leader); 
    }	

    rem = n = strlen(line);

    p = line;
    out = 0;
    while (out < n) {
	*tmp = 0;
	q = p;
	strncat(tmp, p, SAFELEN - 1);
	trim_to_length(tmp, SAFELEN);
	out += strlen(tmp);
	rem = n - out;
	p = q + strlen(tmp);
	if (cli) {
	    printf("%s%s\n", tmp, (rem > 0)? split : "");
	}
	if (!batch) {
	    pprintf(prn, "%s%s\n", tmp, (rem > 0)? split : "");
	}
    }
}

void safe_print_line (const char *line, PRN *prn)
{
    real_safe_print_line(line, 0, 0, 1, prn);
}

static int
print_maybe_quoted_str (const char *s, int cli, PRN *prn)
{
    int ret = 0;

    if (strchr(s, ' ') != NULL) {
	if (cli) {
	    ret += printf(" \"%s\"", s);
	} else {
	    pprintf(prn, " \"%s\"", s);
	    ret += strlen(s) + 3;
	}
    } else {
	if (cli) {
	    ret += printf(" %s", s); 
	} else {
	    pprintf(prn, " %s", s);
	    ret += strlen(s) + 1;
	}
    }

    return ret;
}

static int 
cmd_list_print_var (const CMD *cmd, int i, const DATAINFO *pdinfo,
		    int echo_stdout, PRN *prn)
{
    int v = cmd->list[i];
    int src, genpos;
    int bytes = 0;

    if (v > 0 &&
	(genpos = is_auto_generated_lag(v, cmd->linfo)) > 0) {
	if (is_first_lag(genpos, cmd->linfo, &src)) {
	    bytes += print_lags_by_varnum(src, cmd->linfo, echo_stdout, 
					  pdinfo, prn);
	} 
    } else {
	if (echo_stdout) {
	    bytes += printf(" %s", pdinfo->varname[v]);
	} else {
	    pputc(prn, ' ');
	    bytes += 1 + pputs(prn, pdinfo->varname[v]);
	}
    }

    return bytes;
}

static int more_coming (const CMD *cmd, int i)
{
    if (cmd->opt && cmd->ci != REMEMBER) {
	return 1;
    } else if (cmd->linfo == NULL) {
	return (i < cmd->list[0]);
    } else {
	int j, v, pos;

	for (j=i+1; j<=cmd->list[0]; j++) {
	    v = cmd->list[j];
	    pos = is_auto_generated_lag(v, cmd->linfo);
	    if (pos == 0 || is_first_lag(pos, cmd->linfo, NULL)) {
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

#define listsep_switch(c) (c == AR || c == GARCH || c == ARMA)

#define hold_param(c) (c == TSLS || c == AR || c == ARMA || \
                       c == CORRGM || c == SCATTERS || \
                       c == GNUPLOT || c == LOGISTIC || c == GARCH || \
                       c == EQUATION || c == POISSON)

#define TESTLEN 62
#define LINELEN 78

static void
print_cmd_list (const CMD *cmd, const DATAINFO *pdinfo,  
		int batch, int echo_stdout, char leadchar,
		int *stdlen, int *prnlen, PRN *prn)
{
    int use_varnames = 1;
    char first[16];
    int nsep, gotsep, i;

    nsep = n_separators(cmd->list);

    if (cmd->ci == AR || cmd->ci == GARCH || cmd->ci == ARMA) {
	use_varnames = 0;	
    }

    if (cmd->ci == REMEMBER && cmd->extra != NULL) {
	sprintf(first, "# %s", cmd->word);
    } else {
	strcpy(first, cmd->word);
    }

    if (echo_stdout) {
	if (batch) {
	    if (cmd->ci != EQUATION) {
		putchar('\n');
	    }
	    *stdlen += printf("%c %s", leadchar, first);
	} else {
	    *stdlen += printf(" %s", first);
	}
	if (cmd->ci == RHODIFF) {
	    *stdlen += printf(" %s;", cmd->param);
	} else if (cmd->ci == LAGS) {
	    if (cmd->param[0] != '\0') {
		*stdlen += printf(" %s;", cmd->param);
	    }
	} else if (cmd->param[0] != '\0' && !hold_param(cmd->ci)) {
	    *stdlen += print_maybe_quoted_str(cmd->param, 1, prn);
	}
	if (cmd->ci == REMEMBER) {
	    fputs(" =", stdout);
	    *stdlen += 2;
	    if (cmd->list[0] == 0) {
		fputs(" null", stdout);
		*stdlen += 5;
	    }
	} else if (cmd->ci == VECM) {
	    *stdlen += printf(" %s", cmd->extra);
	}
    }

    if (!batch) {
	pprintf(prn, "%s", first);
	*prnlen += strlen(first);
	if (cmd->ci == RHODIFF) {
	    pprintf(prn, " %s;", cmd->param);
	    *prnlen += 2 + strlen(cmd->param);
	} else if (cmd->ci == LAGS) {
	    if (cmd->param[0] != '\0') {
		pprintf(prn, " %s;", cmd->param);
		*prnlen += 2 + strlen(cmd->param);
	    }
	} else if (cmd->param[0] != '\0' && !hold_param(cmd->ci)) {
	    *prnlen += print_maybe_quoted_str(cmd->param, 0, prn);
	}
	if (cmd->ci == REMEMBER) {
	    *prnlen += pputs(prn, " =");
	    if (cmd->list[0] == 0) {
		pputs(prn, " null");
		*prnlen += 5;
	    }
	} else if (cmd->ci == VECM) {
	    pprintf(prn, " %s", cmd->extra);
	    *prnlen += 1 + strlen(cmd->extra);
	}
    }

#if 0
    if (cmd->ci == STORE) {
	if (echo_stdout) {
	    printf(" \\\n");
	}
	if (!batch) {
	    pputs(prn, " \\\n");
	}
    }
#endif

    gotsep = 0;

    for (i=1; i<=cmd->list[0]; i++) {

	if (cmd->list[i] == LISTSEP) {
	    if (echo_stdout) {
		*stdlen += printf(" ;");
	    }
	    if (!batch) {
		*prnlen += pputs(prn, " ;");
	    }
	    gotsep++;
	    if (listsep_switch(cmd->ci) && gotsep == nsep) {
		use_varnames = !use_varnames;
	    }
	    continue;
	}

	if (echo_stdout) {
	    if (use_varnames) {
		*stdlen += cmd_list_print_var(cmd, i, pdinfo, 1, prn);
	    } else {
		*stdlen += printf(" %d", cmd->list[i]);
	    }
	    if (*stdlen > TESTLEN && more_coming(cmd, i)) {
		printf(" \\\n "); 
		*stdlen = 1;
	    }
	}

	if (!batch) {
	    if (use_varnames) {
		*prnlen += cmd_list_print_var(cmd, i, pdinfo, 0, prn);
	    } else {
		char numstr[12];
		
		sprintf(numstr, " %d", cmd->list[i]);
		*prnlen += pputs(prn, numstr);
	    }
	    if (*prnlen > TESTLEN && more_coming(cmd, i)) {
		pputs(prn, " \\\n "); 
		*prnlen = 1;
	    }
	}
    }
}

#undef ECHO_DEBUG

static int is_silent (const CMD *cmd)
{
    if (cmd->ci == FUNCERR) {
	return 1;
    }

    if (cmd->ci == SET && !strcmp(cmd->param, "echo") &&
	gretl_executing_function()) {
	return 1;
    }

    return 0;
}

/**
 * echo_cmd:
 * @cmd: pointer to #CMD struct.
 * @pdinfo: pointer to data information struct.
 * @line: "raw" command line associated with @cmd.
 * @flags: bitwise OR of elements from #CmdEchoFlags.
 * @prn: pointer to gretl printing struct (or %NULL).
 *
 * Echoes the user command represented by @pcmd and @line, to
 * %stdout and/or @prn.
 */

void echo_cmd (const CMD *cmd, const DATAINFO *pdinfo, const char *line, 
	       unsigned char flags, PRN *prn)
{
    char leadchar = (flags & CMD_STACKING)? '>' : '?';
    int echo_stdout = (flags & CMD_ECHO_TO_STDOUT);
    int batch = (flags & CMD_BATCH_MODE);
    int len, stdlen = 0, prnlen = 0;

    if (line == NULL) {
	return;
    }

#if ECHO_DEBUG
    fprintf(stderr, "echo_cmd: line='%s', echo_stdout=%d, cmd->opt=%ld, batch=%d, "
	    "param='%s', nolist=%d\n", line, echo_stdout, cmd->opt, batch, cmd->param,
	    cmd->nolist);
    fprintf(stderr, " prn=%p\n", (void *) prn);
    fprintf(stderr, " cmd->word='%s'\n", cmd->word);
#endif

    /* don't echo certain things */
    if (is_silent(cmd)) {
	return;
    }

    /* special case: gui "store" command, which could overflow the
       line length; also I'm not sure whether we should record gui
       "store" in the command script.  As a compromise we'll record it,
       but commented out. (FIXME: loop context?)
    */
    if (!echo_stdout && !batch && cmd->ci == STORE) {  
	pprintf(prn, "# store '%s'", cmd->param);
	if (cmd->opt) { 
	    const char *flagstr = print_flags(cmd->opt, cmd->ci);

	    pputs(prn, flagstr);
	}
	pputc(prn, '\n');
	return;
    }

    /* another special, REMEMBER: we should print the literal
       command line first */
    if (cmd->ci == REMEMBER && cmd->extra != NULL) {
	pputs(prn, cmd->extra);
	pputc(prn, '\n');
    }

    if (*line == '\0' || *line == '!' || !strcmp(line, "quit")) {
	return;
    }

    /* command is preceded by a "savename" to which an object will
       be assigned */
    if (*cmd->savename && !echo_stdout && !batch) {
	pprintf(prn, "%s <- ", cmd->savename);
	prnlen += strlen(cmd->savename) + 4;
    }

    if (!cmd->nolist) { 
	/* command has a list of args to be printed */
	print_cmd_list(cmd, pdinfo, batch, echo_stdout, leadchar, 
		       &stdlen, &prnlen, prn);
    } else if ((cmd->ci == GENR || cmd->ci == SMPL) && 
	       strlen(line) > SAFELEN - 2) {
	real_safe_print_line(line, echo_stdout, batch,  
			     (flags & CMD_STACKING), prn);
    } else if (strcmp(cmd->word, "quit")) {
	if (echo_stdout) {
	    if (batch) {
		stdlen += printf("%c %s", leadchar, line);
	    } else {
		stdlen += printf(" %s", line);
	    }
	}
	if (!batch) {
	    prnlen += pputs(prn, line);
	}
    }

    /* print parameter after list, if wanted */
    if (cmd->ci == LOGISTIC) {
	if (cmd->param[0] != '\0') {
	    len = strlen(cmd->param) + 1;
	    if (echo_stdout) {
		putchar(' ');
		fputs(cmd->param, stdout);
		stdlen += len;
	    }
	    if (!batch) {
		pputc(prn, ' ');
		pputs(prn, cmd->param);
		prnlen += len;
	    }
	}
    }

    /* add printout of any options to the command (note that these
       will have been stripped from line)
    */
    if (cmd->opt) {
	const char *flagstr;
	int ci = cmd->ci;

	if (cmd->ci == END) {
	    if (!strcmp(cmd->param, "nls")) {
		ci = NLS;
	    } else if (!strcmp(cmd->param, "mle")) {
		ci = MLE;
	    }
	}
	flagstr = print_flags(cmd->opt, ci);
	len = strlen(flagstr);
	if (echo_stdout) {
	    if (stdlen + len > LINELEN) {
		fputs(" \\\n ", stdout);
	    }
	    fputs(flagstr, stdout);
	}
	if (!batch) {
	    if (prnlen + len > LINELEN) {
		pputs(prn, " \\\n ");
	    }	    
	    pputs(prn, flagstr);
	}
    }

    if (echo_stdout) {
	putchar('\n');
    }

    if (!batch) {
	pputc(prn, '\n');
	gretl_print_flush_stream(prn);
    }
}

void echo_function_call (const char *line, unsigned char flags, PRN *prn)
{
    char leadchar = (flags & CMD_STACKING)? '>' : '?';

    if (gretl_echo_on()) {
	pprintf(prn, "%c %s\n", leadchar, line);
    }
}

/* Look for a flag of the form "-x" which occurs outside of any
   quotes: if found, return a pointer to the flag.
*/

static const char *flag_present (const char *s, char f, int *quoted)
{
    int inquote = 0;
    int gotdash = 0;

    while (*s) {
	if (*s == '"') {
	    inquote = !inquote;
	}
	if (!inquote) {
	    if (*s == '-') {
		gotdash = 1;
	    } else if (gotdash && *s == f && *(s+1)) {
		s++;
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
	    } else {
		gotdash = 0;
	    }
	}
	s++;
    }

    return NULL;
}

static char *get_flag_field  (const char *s, char f)
{
    const char *p;
    char *ret = NULL;
    int quoted = 0;

    if ((p = flag_present(s, f, &quoted)) != NULL) {
	const char *q = p;
	size_t len = 0;

	while (*q) {
	    if (quoted && *q == '"') break;
	    if (!quoted && isspace(*q)) break;
	    q++;
	    len++;
	}

	ret = malloc(len + 1);

	if (ret != NULL) {
	    *ret = 0;
	    strncat(ret, p, len);
	}
    }

    return ret;
}

static void get_optional_filename (const char *line, CMD *cmd)
{
    char *p;

    p = get_flag_field(line + 8, 'f');
    if (p != NULL) {
	free(cmd->param);
	cmd->param = p;
    }    
}

static int set_var_info (const char *line, gretlopt opt, 
			 DATAINFO *pdinfo, PRN *prn)
{
    char *p;
    char vname[VNAMELEN];
    int cmdlen = 0;
    int v, setstuff = 0;

    if (pdinfo->varinfo == NULL) {
	return 1;
    }

    /* skip command word */
    cmdlen = strcspn(line, " ");
    line += cmdlen++;
    line += strspn(line, " ");

    if (sscanf(line, "%15s", vname) != 1) {
	return E_PARSE;
    }

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

    p = get_flag_field(line + cmdlen, 'd');
    if (p != NULL) {
	setstuff = 1;
	*VARLABEL(pdinfo, v) = 0;
	strncat(VARLABEL(pdinfo, v), p, MAXLABEL - 1);
	free(p);
    }

    p = get_flag_field(line + cmdlen, 'n');
    if (p != NULL) {
	setstuff = 1;
	*DISPLAYNAME(pdinfo, v) = 0;
	strncat(DISPLAYNAME(pdinfo, v), p, MAXDISP - 1);
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
do_outfile_command (gretlopt flag, char *fname, PRN *prn)
{
    static char outname[MAXLEN];
    int diverted = 0;

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
	    pprintf(prn, "Closed output file '%s'\n", outname);
	    return 0;
	}
    }

    /* command to divert output to file */
    if (diverted) {
	fprintf(stderr, _("Output is already diverted to '%s'\n"),
		outname);
	return 1;
    } else {
	if (*fname == '\0') {
	    return E_ARGS;
	} else {
	    FILE *fp;

	    if (flag == OPT_W) {
		fp = gretl_fopen(fname, "w");
	    } else {
		fp = gretl_fopen(fname, "a");
	    }

	    if (fp == NULL) {
		pprintf(prn, _("Couldn't open %s for writing\n"), fname);
		return 1;
	    } else {
		if (flag == OPT_W) {
		    pprintf(prn, _("Now writing output to '%s'\n"), fname);
		} else {
		    pprintf(prn, _("Now appending output to '%s'\n"), fname);
		}
	    }

	    print_start_redirection(prn, fp);
	    strcpy(outname, fname);
	    return 0;
	}
    }

    return 1; /* not reached */
}

int call_pca_plugin (VMatrix *corrmat, double ***pZ,
		     DATAINFO *pdinfo, gretlopt *pflag,
		     PRN *prn)
{
    void *handle = NULL;
    int (*pca_from_corrmat) (VMatrix *, double ***, DATAINFO *,
			     gretlopt *, PRN *);
    int err = 0;

    *gretl_errmsg = 0;
    
    pca_from_corrmat = get_plugin_function("pca_from_corrmat", &handle);
    if (pca_from_corrmat == NULL) {
        return 1;
    }
        
    err = (* pca_from_corrmat) (corrmat, pZ, pdinfo, pflag, prn);
    close_plugin(handle);
    
    return err;
}

static int add_obs (int n, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    int err = 0;

    if (complex_subsampled()) {
	pprintf(prn, _("The data set is currently sub-sampled.\n"));
	err = E_DATA;
    } else if (n <= 0) {
	err = E_PARSE;
    } else {
	err = dataset_add_observations(n, pZ, pdinfo, OPT_A);
	if (!err) {
	    pprintf(prn, _("Dataset extended by %d observations"), n);
	    pputc(prn, '\n');
	}
    }

    return err;
}

/* common code for command-line and GUI client programs, where the
   command doesn't require special handling on the client side 
*/

int simple_commands (CMD *cmd, const char *line, 
		     double ***pZ, DATAINFO *pdinfo,
		     PRN *prn)
{
    int err = 0, order = 0;
    VMatrix *corrmat;
    Summary *summ;
    int *listcpy = NULL;

    if (RETURNS_LIST(cmd->ci)) {
	/* list is potentially modified -> make a copy */
	listcpy = gretl_list_copy(cmd->list);
	if (listcpy == NULL) {
	    return E_ALLOC;
	}
    }

    switch (cmd->ci) {

    case ADDOBS:
	err = add_obs(cmd->order, pZ, pdinfo, prn);
	break;

    case ADF:
	err = adf_test(cmd->order, cmd->list[1], pZ, pdinfo, cmd->opt, prn);
	break;

    case KPSS:
	err = kpss_test(cmd->order, cmd->list[1], pZ, pdinfo, cmd->opt, prn);
	break;

    case COINT:
	err = coint(cmd->order, cmd->list, pZ, pdinfo, cmd->opt, prn);
	break;

    case COINT2:
	err = johansen_test_simple(cmd->order, cmd->list, pZ, pdinfo, 
				   cmd->opt, prn);
	break;

    case CORR:
	if (cmd->list[0] > 3) {
	    err = gretl_corrmx(cmd->list, (const double **) *pZ, pdinfo, 
			       prn);
	    if (err) 
		pputs(prn, _("Error in generating correlation matrix\n"));
	    break;
	}
	corrmat = corrlist(cmd->list, (const double **) *pZ, pdinfo);
	if (corrmat == NULL) {
	    pputs(prn, _("Couldn't allocate memory for correlation matrix.\n"));
	} else {
	    printcorr(corrmat, prn);
	    free_vmatrix(corrmat);
	}
	break;

    case DISCRETE:
	err = list_makediscrete(cmd->list, pdinfo, cmd->opt);
	break;

    case ESTIMATE:
	err = estimate_named_system(line, pZ, pdinfo, cmd->opt, prn);
	break;

    case FNCALL:
	err = gretl_function_start_exec(line, cmd->param, pZ, pdinfo);
	break;

    case FUNC:
	err = gretl_start_compiling_function(line, prn);
	break;

    case FUNCERR:
	err = gretl_function_flagged_error(cmd->param, prn);
	break;	

    case PCA:
	corrmat = corrlist(cmd->list, (const double **) *pZ, pdinfo);
	if (corrmat == NULL) {
	    pputs(prn, _("Couldn't allocate memory for correlation matrix.\n"));
	} else {
	    err = call_pca_plugin(corrmat, pZ, pdinfo, &cmd->opt, prn);
	    if (cmd->opt && !err) {
		maybe_list_vars(pdinfo, prn);
	    }
	    free_vmatrix(corrmat);
	}
	break;

    case CRITERIA:
	err = parse_criteria(line, (const double **) *pZ, pdinfo, prn);
	if (err) { 
	    pputs(prn, _("Error in computing model selection criteria.\n"));
	}
	break;

    case CRITICAL:
	err = print_critical(line, prn);
	break;

    case DATA:
	err = db_get_series(line, pZ, pdinfo, prn);
	break;

    case DIFF:
    case LDIFF:
    case SDIFF:
	err = list_diffgenr(listcpy, cmd->ci, pZ, pdinfo);
	if (err) {
	    if (cmd->ci == LDIFF) {
		pputs(prn, _("Error adding log differences of variables.\n"));
	    } else if (cmd->ci == DIFF) {
		pputs(prn, _("Error adding first differences of variables.\n"));
	    }
	} else {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case DUMMIFY:
	err = list_dumgenr(&listcpy, pZ, pdinfo);
	if (err) {
	    pputs(prn, _("Error adding dummy variables.\n"));
	} else {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case LAGS:
	order = atoi(cmd->param);
	err = list_laggenr(&listcpy, order, pZ, pdinfo); 
	if (err) {
	    pputs(prn, _("Error adding lags of variables.\n"));
	} else {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case LOGS:
	err = list_loggenr(listcpy, pZ, pdinfo);
	if (err) {
	    pputs(prn, _("Error adding logs of variables.\n"));
	} else {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case SQUARE:
	err = list_xpxgenr(&listcpy, pZ, pdinfo, cmd->opt);
	if (err) {
	    pputs(prn, _("Failed to generate squares\n"));
	} else {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case MULTIPLY:
	err = gretl_multiply(cmd->param, cmd->list, cmd->extra, pZ, pdinfo);
	if (!err) {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case GRAPH:
	ascii_graph(cmd->list, (const double **) *pZ, pdinfo, 
		    cmd->opt, prn);
	break;

    case PLOT:
	ascii_graph(cmd->list, (const double **) *pZ, pdinfo, 
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
	if (pdinfo->descrip != NULL) {
	    pprintf(prn, "%s\n", pdinfo->descrip);
	} else {
	    pputs(prn, _("No data information is available.\n"));
	}
	break;

    case RENAME:
	err = rename_var_by_id(cmd->extra, cmd->param, pdinfo);
	if (err) {
	    errmsg(err, prn);
	} else {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case SETINFO:
	err = set_var_info(line, cmd->opt, pdinfo, prn);
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
	    printdata(cmd->list, (const double **) *pZ, pdinfo, 
		      cmd->opt, prn);
	}
	break;

    case RHODIFF:
	err = rhodiff(cmd->param, cmd->list, pZ, pdinfo);
	if (err) {
	    errmsg(err, prn);
	} else {
	    maybe_list_vars(pdinfo, prn);
	}
	break;

    case SUMMARY:
	summ = summary(cmd->list, (const double **) *pZ, pdinfo, prn);
	if (summ == NULL) {
	    pputs(prn, _("generation of summary stats failed\n"));
	} else {
	    print_summary(summ, pdinfo, prn);
	    free_summary(summ);
	}
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
			prn);
	break;

    case SPEARMAN:
	err = spearman(cmd->list, (const double **) *pZ, pdinfo, 
		       cmd->opt, prn);
	break;

    case OUTFILE:
	err = do_outfile_command(cmd->opt, cmd->param, prn);
	break;

    case STORE:
	if (*cmd->param != '\0') {
	    if ((cmd->opt & OPT_Z) && !has_suffix(cmd->param, ".gz")) {
		pprintf(prn, _("store: using filename %s.gz\n"), cmd->param);
	    } else {
		pprintf(prn, _("store: using filename %s\n"), cmd->param);
	    }
	} else {
	    pprintf(prn, _("store: no filename given\n"));
	    break;
	}
	if (write_data(cmd->param, cmd->list, (const double **) *pZ,
		       pdinfo, cmd->opt, NULL)) {
	    pprintf(prn, _("write of data file failed\n"));
	    err = 1;
	    break;
	}
	pprintf(prn, _("Data written OK.\n"));
	if (((cmd->opt & OPT_O) || (cmd->opt & OPT_S)) && pdinfo->markers) {
	    pprintf(prn, _("Warning: case markers not saved in "
			   "binary datafile\n"));
	}
	break;

    case REMEMBER:
	if (cmd->opt & OPT_P) {
	    ; /* let echo do the work? */
	} else if (cmd->opt & OPT_L) {
	    err = remember_list(cmd->list, cmd->param, prn);
	} 
	break;

    case MATRIX:
	err = matrix_command(line, pZ, pdinfo, prn);
	break;

    case TRANSPOSE:
	err = transpose_data(pZ, pdinfo);
	break;

    case SHELL:
	err = gretl_shell(line + 1);
	break;

    default:
	break;
    }

    if (listcpy != NULL) {
	free(listcpy);
    }

    if (err == E_OK) {
	err = 0;
    }

    return err;
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
 * for execution within a loop.
 *
 * Returns: 1 on error, otherwise 0.
 */

int get_command_index (char *line, CMD *cmd, const DATAINFO *pdinfo)
{
    static int context;
    int done = 0;

    while (isspace(*line)) {
	line++;
    }

    cmd->ci = 0;

#if CMD_DEBUG
    fprintf(stderr, "get_command_index: line='%s'\n", line);
#endif

    if (*line == '#' || (*line == '(' && *(line+1) == '*')) {
	cmd->nolist = 1;
	cmd->ci = CMD_COMMENT;
	return 0;
    }

    if (sscanf(line, "%8s", cmd->word) != 1) {
	cmd->nolist = 1;
	cmd->ci = CMD_NULL;
	return 0;
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
	    if (gretl_is_user_function(line)) {
		cmd->nolist = 1;
		cmd->ci = FNCALL;
	    } else if (!plausible_genr_start(line, cmd, pdinfo)) {
		cmd->errcode = 1;
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

int ready_for_command (const char *line)
{
    const char *ok_cmds[] = {
	"open", 
	"run", 
	"include",
	"nulldata", 
	"import", 
	"pvalue",
	"print",
	"printf",
	"eval",
	"!",
	"(*", 
	"man ", 
	"help", 
	"set", 
	"critical", 
	"seed", 
	"function",
	"newfunc",
	"noecho",
	NULL 
    };
    int i, ok = 0;

    if (string_is_blank(line) || gretl_compiling_function()) {
	ok = 1;
    } else if (*line == 'q' || *line == 'x' || *line == '#') {
	ok = 1;
    } else {
	for (i=0; ok_cmds[i] != NULL && !ok; i++) {
	    if (strncmp(line, ok_cmds[i], strlen(ok_cmds[i])) == 0) {
		ok = 1;
	    }
	}
    }

    return ok;
}

int gretl_cmd_init (CMD *cmd)
{
    cmd->ci = 0;
    cmd->errcode = 0;
    cmd->context = 0;
    cmd->order = 0;
    cmd->ignore = 0;
    *cmd->word = '\0';

    cmd->list = NULL;
    cmd->param = NULL;
    cmd->extra = NULL;
    cmd->linfo = NULL;

    /* make 'list', 'param' and 'extra' blank rather than NULL
       for safety (in case they are deferenced) */

    cmd->list = gretl_null_list();
    if (cmd->list == NULL) {
	cmd->errcode = E_ALLOC;
    }

    if (cmd->errcode == 0) {
	cmd->param = calloc(1, 1);
	if (cmd->param == NULL) {
	    cmd->errcode = E_ALLOC;
	}
    }

    if (cmd->errcode == 0) {
	cmd->extra = calloc(1, 1);
	if (cmd->extra == NULL) {
	    free(cmd->param);
	    cmd->param = NULL;
	    cmd->errcode = E_ALLOC;
	}
    }    

    return cmd->errcode;
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


