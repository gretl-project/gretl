/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
#include "internal.h"

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

extern int _omitfromlist (int *list, const int *omitvars, int newlist[],
			  const DATAINFO *pdinfo, const int model_count);

static int _full_list (const DATAINFO *pdinfo, CMD *command);
static void get_optional_filename (const char *line, CMD *cmd);

typedef struct {
    int lag;
    int varnum;
    char varname[VNAMELEN];
} LAGVAR;

/* ........................................................... */

static int trydatafile (char *line, int *ignore)
{
    int i, m, n = strlen(line);
    char datfile[MAXLEN];

    *datfile = '\0';
    for (i=0; i<n; i++) {
	if ((n - i) > 4 && strncmp(line+i, "DATA", 4) == 0) {
	    sscanf(line + i, "%s", datfile);
	    m = strlen(datfile);
	    if (datfile[m-1] == ',') datfile[m-1] = '\0';
	    lower(datfile);
	    i += 4;
	}
	else if (line[i] == '*' && line[i+1] == ')') *ignore = 0;
    }
    if (*datfile) {
	sprintf(line, "open %s.gdt", datfile);
	return 1;
    } 
    return 0;
}

/* ........................................................... */

static int filter_comments (char *line, int *ignore)
{
    int i, j = 0, n = strlen(line);
    char tmpstr[MAXLEN], datfile[MAXLEN];
    
    for (i=0; i<n; i++) {
	if (line[i] == '(' && line [i+1] == '*') {
	    *ignore = 1;
	    if (line[i+3] == '!') { /* special code for data file to open */
		sscanf(line + 4, "%s", datfile);
		sprintf(line, "open %s", datfile);
		*ignore = 0;  /* FIXME ? */
		return 0;
	    }
	}
	else if (line[i] == '*' && line [i+1] == ')') {
	    *ignore = 0; i += 2;
	    while (isspace((unsigned char) line[i]) && i < n) i++;
	}
	if (!(*ignore) && line[i] != '\r') {
	    tmpstr[j] = line[i];
	    j++;
	}
    }
    tmpstr[j] = '\0';
    strcpy(line, tmpstr);

    if (*line) return 0;
    else return 1;
}

/* ........................................................... */

static int get_rhodiff_param (char *str, CMD *cmd)
{
    int k;

    if ((k = haschar(';', str)) < 0) return 1;

    cmd->param = realloc(cmd->param, k+1);
    if (cmd->param == NULL) return E_ALLOC;

    *cmd->param = 0;
    strncat(cmd->param, str, k);

    _shiftleft(str, k + 1);

    return 0;
}

/* ........................................................... */

int subsetted_command (const char *cmd)
{    
    if (strcmp(cmd, "deriv") == 0) return NLS;
    return 0;
}

/* ........................................................... */

static int aliased (char *cmd)
{
    if (!strcmp(cmd, "q")) {
	strcpy(cmd, "quit");
	return 1;
    }
    else if (!strcmp(cmd, "x")) {
	strcpy(cmd, "quit");
	return 2;
    }
    else if (!strcmp(cmd, "let")) {
	strcpy(cmd, "genr");
	return 1;
    }
    else if (!strcmp(cmd, "ls")) {
	strcpy(cmd, "varlist");
	return 1;
    }
    else if (!strcmp(cmd, "list")) {
	strcpy(cmd, "varlist");
	return 1;
    }
    else if (!strcmp(cmd, "boxplots")) { 
	strcpy(cmd, "boxplot");
	return 1;
    }
    else if (!strcmp(cmd, "man")) {
	strcpy(cmd, "help");
	return 1;
    }
    else if (!strcmp(cmd, "sample")) {
	strcpy(cmd, "smpl");
	return 1;
    }
    else if (!strcmp(cmd, "eval")) {
	strcpy(cmd, "genr");
	return 1;
    }
    else if (cmd[0] == '!') {
	strcpy(cmd, "shell");
	return 1;
    }
    return 0;
}

#define NO_VARLIST(c) (c == VARLIST || \
	               c == QUIT || \
	               c == EQNPRINT || \
	               c == TABPRINT || \
	               c == FCAST || \
	               c == FCASTERR || \
	               c == FIT || \
 	               c == LABEL || \
 	               c == LABELS || \
    	               c == INFO || \
	               c == CRITERIA || \
	               c == PVALUE || \
	               c == RUN || \
	               c == SHELL || \
	               c == SETOBS || \
	               c == CHOW || \
	               c == CUSUM || \
	               c == CRITICAL || \
	               c == HAUSMAN || \
	               c == PANEL || \
	               c == OPEN || \
	               c == IMPORT || \
	               c == ENDLOOP || \
	               c == SIM || \
                       c == RENAME || \
	               c == TESTUHAT || \
                       c == RESET || \
                       c == SYSTEM || \
                       c == LEVERAGE || \
                       c == MODELTAB || \
                       c == NLS || \
                       c == DATA || \
	               c == GENR || \
                       c == SET || \
                       c == PRINTF || \
                       c == OUTFILE)

/* ........................................................... */

static int flow_control (const char *line, double ***pZ, 
			 DATAINFO *pdinfo, CMD *cmd)
{
    /* clear to proceed */
    if (!ifstate(IS_FALSE) && 
	cmd->ci != IF && cmd->ci != ELSE && cmd->ci != ENDIF)
	return 0;

    if (cmd->ci == IF) {
	int ret = if_eval(line, pZ, pdinfo);

	if (ret == -1 || ifstate((ret)? SET_TRUE : SET_FALSE)) 
	    cmd->errcode = E_SYNTAX;
    }

    else if (cmd->ci == ELSE && ifstate(SET_ELSE)) 
	cmd->errcode = E_SYNTAX;

    else if (cmd->ci == ENDIF && ifstate(SET_ENDIF)) 
	cmd->errcode = E_SYNTAX;

    return 1;
}

static void get_savename (char *s, CMD *cmd)
{
    *cmd->savename = 0;

    if (strncmp(s, "genr ", 5) && strstr(s, " <- ")) {
	int n, len, quote;

	quote = (*s == '"');
	len = strcspn(s, "<");
	if (len < 2) return;
	n = len - 1 - quote;
	if (n > MAXSAVENAME - 1) n = MAXSAVENAME - 1;
	strncat(cmd->savename, s + quote, n);
	if (cmd->savename[n-1] == '"') cmd->savename[n-1] = 0;
	_shiftleft(s, len + 3);
    }
}

static void grab_gnuplot_literal_block (char *line, CMD *command)
{
    char *p = strchr(line, '{');
    char *bl = NULL;

    if (p != NULL) {
	bl = malloc(strlen(p) + 1);
	if (bl != NULL) {
	    strcpy(bl, p);
	    free(command->param);
	    command->param = bl;
	    *p = 0;
	}
    }
}

static int parse_lagvar (const char *s, LAGVAR *plagv, DATAINFO *pdinfo)
{
    int ret = 0;

    *plagv->varname = 0;
    plagv->lag = 0;
    plagv->varnum = 0;

    if (sscanf(s, "%8[^(](-%d)", plagv->varname, &plagv->lag) == 2) {
	plagv->varnum = varindex(pdinfo, plagv->varname);
	if (plagv->varnum < pdinfo->v) {
	    ret = 1;
	} 
    }

    return ret;
}

static void parse_rename_cmd (const char *line, CMD *cmd, 
			      const DATAINFO *pdinfo)
{
    int v, vnum;
    char vname[VNAMELEN];

    line += strlen("rename ");

    if (sscanf(line, "%d %8s", &vnum, vname) != 2) {
	cmd->errcode = E_DATA;
	sprintf(gretl_errmsg, "rename: %s", 
		_("requires a variable number and a new name"));
	return;
    }

    if (vnum >= pdinfo->v || vnum < 1) {
	cmd->errcode = E_DATA;
	sprintf(gretl_errmsg, _("Variable number %d is out of bounds"), vnum);
	return;
    } 

    v = varindex(pdinfo, vname);
    if (v < pdinfo->v && v != vnum) {
	sprintf(gretl_errmsg, _("'%s': there is already a variable "
				"of this name"), vname);
	cmd->errcode = E_DATA;
	return;
    }

    if (check_varname(vname)) {
	cmd->errcode = E_DATA;
	return;
    }

    free(cmd->param);
    cmd->param = malloc(strlen(vname) + 1);
    if (cmd->param == NULL) {
	cmd->errcode = E_ALLOC;
	return;
    }
    
    strcpy(cmd->param, vname);
    sprintf(cmd->str, "%d", vnum);
}

static void parse_outfile_cmd (char *line, CMD *cmd)
{
    /* 7 = number of chars in the command word, "outfile" */
    int n = strlen(line) - 7;

    if (n > 1) {
	free(cmd->param);
	cmd->param = malloc(n);
	if (cmd->param == NULL) {
	    cmd->errcode = E_ALLOC;
	} else {
	    char *p = line + 7 + 1;

	    while (*p && (isspace(*p) || *p == '"')) p++;
	    strcpy(cmd->param, p);
	    n = strlen(cmd->param) - 1;
	    while (n > 0) {
		if (isspace(cmd->param[n]) || cmd->param[n] == '"') {
		    cmd->param[n] = 0;
		    n--;
		} else {
		    break;
		}
	    }
	}
    }
}

static void parse_logistic_ymax (char *line, CMD *cmd)
{
    char *p;

    p = strstr(line, "ymax");
    if (p != NULL) {
	char *q = p + 4;
	char numstr[12];

	while (*q == ' ' || *q == '=') q++;
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



/**
 * getcmd:
 * @line: the command line (string).
 * @pdinfo: pointer to data information struct.
 * @command: pointer to #CMD struct.
 * @ignore: pointer to int indicating whether (1) or not (0) we're
 * in comment mode, and @line should not be parsed.
 * @pZ: pointer to data matrix.
 * @cmds: pointer to gretl printing struct.
 *
 * Parses @line and fills out @command accordingly.  In case
 * of error, @command->errcode gets a non-zero value.
 */

void getcmd (char *line, DATAINFO *pdinfo, CMD *command, 
	     int *ignore, double ***pZ, PRN *cmds)
{
    int i, j, nf, linelen, n, v, gotdata = 0, ar = 0, poly = 0;
    int spacename = 0;
    char field[10], *remainder;
    LAGVAR lagvar;

    command->ci = 0;
    command->errcode = 0;
    command->nolist = 0;
    *command->param = '\0';

    *gretl_errmsg = '\0';

    /* look for ramu practice files */
    if (line[0] == '(' && line[1] == '*') {
	*ignore = 1;
	gotdata = trydatafile(line, ignore);
    }

    /* trap comments */
    if (!gotdata) {
	if (filter_comments(line, ignore)) {
	    command->nolist = 1;
	    command->ci = -2;
	    return;
	}
    }

    /* extract "savename" for storing an object? */
    get_savename(line, command);

    linelen = strlen(line);

    if (*line == '#' || sscanf(line, "%s", command->cmd) != 1) {
	command->nolist = 1;
	command->ci = -1;
	return;
    }

    /* backwards compatibility for obsolete commands */
    if (!strcmp(command->cmd, "noecho")) {
	strcpy(command->cmd, "set");
	strcpy(line, "set echo off");
    } else if (!strcmp(command->cmd, "seed")) {
	char seedstr[16];

	strcpy(command->cmd, "set");
	if (sscanf(line, "%*s %15s", seedstr)) {
	    sprintf(line, "set seed %s", seedstr);
	} else {
	    strcpy(line, "set seed");
	}
    }

    /* command aliases */
    if (aliased(command->cmd) == 2) {
	command->param = realloc(command->param, 2);
	strcpy(command->param, "x");
    }

    /* subsetted commands (e.g. "deriv" in relation to "nls") */
    command->ci = subsetted_command(command->cmd);

    /* trap bogus commands */ 
    if (command->ci == 0) {
	command->ci = gretl_command_number(command->cmd);
	if (command->ci == 0) {
	    command->errcode = 1;
	    sprintf(gretl_errmsg, _("command '%s' not recognized"), 
		    command->cmd);
	    return;
	}
    }

    /* if, else, endif controls */
    if (flow_control(line, pZ, pdinfo, command)) {
	command->nolist = 1;
	command->ci = -1;
	return;
    }

    /* tex printing commands can take a filename parameter; so
       can gnuplot command */
    if (command->ci == EQNPRINT || command->ci == TABPRINT) {
	get_optional_filename(line, command);
    } 

    /* the "outfile" command may have a filename */
    else if (command->ci == OUTFILE) {
	parse_outfile_cmd(line, command);
    }

    /* the "rename" command calls for a variable number and a
       new name */
    else if (command->ci == RENAME) {
	parse_rename_cmd(line, command, pdinfo);
    }  

    /* "logistic" can have a special parameter */
    else if (command->ci == LOGISTIC) {
	parse_logistic_ymax(line, command);
    }      

    /* commands that never take a list of variables */
    if (NO_VARLIST(command->ci)) {
	command->nolist = 1;
	return;
    }

    /* smpl can be special: only takes a list in case of OPT_M
       "--no-missing" */
    if (command->ci == SMPL && !(command->opt & OPT_M)) {
	command->nolist = 1;
	return;
    }	

    /* boxplots can be special: boolean conditions embedded in
       the list, which have to be parsed separately */
    if (command->ci == BXPLOT && strchr(line, '(')) {
	command->nolist = 1;
	return;
    }

    /* gnuplot command can have a block of stuff to pass literally
       to gnuplot */
    if (command->ci == GNUPLOT) {
	grab_gnuplot_literal_block(line, command);
    }

    /* fix lines that contain a semicolon right after a var */
    for (i=0; i<linelen-1; i++) {
	if (isalnum((unsigned char) line[i]) && line[i+1] == ';') {
	    for (j=linelen; j>i+1; j--) line[j] = line[j-1];
	    line[i+1] = ' ';
	    line[linelen + 1] = '\0';
	    linelen = strlen(line);
	}
    }

    /* now we're probably dealing with a command that wants a list... */    
    nf = _count_fields(line) - 1;
    n = strlen(command->cmd);

    /* unless it's on a short list of specials */
    if (command->ci == HELP ||
	command->ci == LOOP ||
	command->ci == END ||
	command->ci == LMTEST ||
	command->ci == NULLDATA ||
	(command->ci == PRINT && strstr(line, "\""))) {
	command->nolist = 1;
	if (!strncmp(line, "man ", 4)) n--;
	if (nf) {
	    command->param = realloc(command->param, linelen - n + 1);
	    if (command->ci == PRINT) {
		strcpy(command->param, line + n + 1);
	    } else {
		sscanf(line + n + 1, "%s", command->param);
	    }
	}
	return;
    }

    remainder = malloc(linelen - n + 1);
    if (remainder == NULL) {
	command->errcode = E_ALLOC;
	return;
    }

    /* need to treat rhodiff specially -- put everything from
       the end of the command word to the first semicolon into
       "param", for processing later */
    if (command->ci == RHODIFF) {  /* FIXME */
	strcpy(remainder, line + n + 1);
	if (get_rhodiff_param(remainder, command)) {
	    command->errcode = E_SYNTAX;
	    free(remainder);
	    return;
	}
	strcpy(line, remainder);
	linelen = strlen(line);
	nf = _count_fields(line);
	n = 0;
    }

    /* "store" is a special case since the filename that comes
       before the list may be quoted, and have spaces in it.  Groan */
    if (command->ci == STORE && nf) {
	int q, qmatch = 0;

	if ((q = line[n+1]) == '"' || (q = line[n+1]) == '\'') {
	    spacename = 1;
	    command->param = realloc(command->param, linelen - n + 1);
	    if (command->param == NULL) {
		command->errcode = E_ALLOC;
		return;
	    }
	    strcpy(remainder, line + n + 1);
	    for (i=1; i< (int) strlen(remainder); i++) {
		if (remainder[i] == q) {
		    strncpy(command->param, remainder, i + 1);
		    command->param[i+1] = '\0';
		    qmatch = 1;
		    break;
		}
		if (remainder[i] == ' ') nf--;
	    }
	    if (!qmatch) {
		command->errcode = E_SYNTAX;
		free(remainder);
		return;
	    }
	    /* fprintf(stderr, "got filename '%s'\n", command->param); */
	    _shiftleft(remainder, strlen(command->param));
	    /* unquote the filename */
	    for (i=0; i< (int) strlen(command->param) - 2; i++) {
		command->param[i] = command->param[i+1];
	    }
	    command->param[i] = '\0';
	    nf--;
	    n = 0;
	    if (nf > 0) {
		strcpy(line, remainder);	
		linelen = strlen(line);
	    }
	}
    } /* end if STORE && nf */

    /* "store" takes a filename before the list, "var" takes a 
       lag order, "adf" takes a lag order, "arch" takes a lag 
       order, "multiply" takes a multiplier.  "omitfrom" and
       "addto" take the ID of a previous model. "setmiss" takes
       a value to be interpreted as "missing". 
    */
    if ((command->ci == STORE && !spacename) ||
	command->ci == ADF ||
	command->ci == ARCH ||
	command->ci == COINT ||
	command->ci == COINT2 ||
	command->ci == ADDTO ||
	command->ci == OMITFROM ||
	command->ci == MULTIPLY ||
	command->ci == SETMISS ||
	command->ci == VAR) {
	if (nf) {
	    command->param = realloc(command->param, linelen - n + 1);
	    if (command->param == NULL) {
		command->errcode = E_ALLOC;
		return;
	    }
	    strcpy(remainder, line + n + 1);
	    sscanf(remainder, "%s", command->param);
	    _shiftleft(remainder, strlen(command->param));
	    nf--;
	    n = 0;
	    if (nf > 0) {
		strcpy(line, remainder);
		linelen = strlen(line);
	    }
	} 
    }

    if (command->ci == MULTIPLY) {  /* suffix string */
	char suffix[4];

	sscanf(line, "%3s", suffix);
	strcpy(command->str, suffix);
	_shiftleft(line, strlen(suffix));
	nf--;
	n = 0;
	linelen = strlen(line);
    }

    if (command->ci == AR || command->ci == ARMA) ar = 1;

    command->list = realloc(command->list, (1 + nf) * sizeof *command->list);
    if (command->list == NULL) {
	command->errcode = E_ALLOC;
	strcpy (gretl_errmsg, _("Memory allocation failed for command list"));
	free(remainder);
	return;
    }
    command->list[0] = nf;

    /* now assemble the command list */
    for (j=1; j<=nf; j++) {

	strcpy(remainder, line + n + 1);

	/* special: optional lag order for correlogram */
	if (command->ci == CORRGM && j == 2) {
	    command->list[0] = 1;
	    command->param = realloc(command->param, linelen - n + 1);
	    sscanf(remainder, "%s", command->param);
	    break;
	}
	
	sscanf(remainder, "%s", field);
	/* fprintf(stderr, "remainder: %s\n", remainder); */

	if (isalpha((unsigned char) *field)) {
	    /* should be the name of a variable */
	    if (field[strlen(field) - 1] == ';')
		field[strlen(field) - 1] = '\0';
	    if ((v = varindex(pdinfo, field)) <= pdinfo->v - 1) {
		/* yes, it's an existing variable */
		command->list[j] = v;
	    } else {
		/* no: an auto-generated variable? */
		/* Case 1: automated lags:  e.g. 'var(-1)' */
		if (parse_lagvar(field, &lagvar, pdinfo)) {
		    extern int newlag; /* generate.c */
		    int lnum;

		    lnum = laggenr(lagvar.varnum, lagvar.lag, 1, pZ, pdinfo);
		    if (lnum < 0) {
			command->errcode = 1;
			sprintf(gretl_errmsg, 
				_("generation of lag variable failed"));
		    } else { 
			command->list[j] = lnum;
			if (newlag && cmds != NULL) {
			    pprintf(cmds, "genr %s\n", VARLABEL(pdinfo, lnum));
			}
			/* fully handled, get on with it */
			n += strlen(field) + 1;
			continue; 
		    }
		} 
		/* Case 2: special plotting variable */
		else if (!command->errcode && (!strcmp(field, "qtrs") || 
		    !strcmp(field, "months") || !strcmp(field, "time"))) {
		    int pnum = plotvar(pZ, pdinfo, field);

		    if (pnum < 0) {
			command->errcode = 1;
			sprintf(gretl_errmsg, 
				_("Failed to add plotting index variable"));
		    } else {
			command->list[j] = pnum;
			/* fully handled, get on with it */
			n += strlen(field) + 1;
			continue; 
		    }
		} 
		/* last chance: try abbreviating the varname? */
		else if (!command->errcode) {
		    command->errcode = 1; /* presume guilt at this stage */
		    if (strlen(field) > 8) {
			char test[VNAMELEN];

			*test = 0;
			strncat(test, field, 8);
			if ((v = varindex(pdinfo, test)) <= pdinfo->v - 1) {
			    command->list[j] = v;
			    command->errcode = 0;
			} 
		    } 
		    if (command->errcode) {
			sprintf(gretl_errmsg, 
				_("'%s' is not the name of a variable"), field);
		    }
		}

		if (command->errcode) {
		    free(remainder);
		    return;
		}
	    }
	} /* end if isalpha(*field) */

	else if (isdigit(*field)) {
	    /* could be the ID number of a variable */
	    v = atoi(field);
	    if (!ar && !poly && v > pdinfo->v - 1) {
		command->errcode = 1;
		sprintf(gretl_errmsg, 
                       _("%d is not a valid variable number"), v);
		free(remainder);
		return;
	    }	
	    command->list[j] = v;
	}

	else if (*field == ';') {
	    /* could be the separator between two sub-lists */
	    if (command->ci == TSLS || command->ci == AR ||
		command->ci == MPOLS || command->ci == SCATTERS ||
		command->ci == ARMA) {
		command->param = realloc(command->param, 4);
		sprintf(command->param, "%d", j);
		n += strlen(field) + 1;
		command->list[j] = LISTSEP;
		ar = 0; /* turn off acceptance of AR lags */
		if (command->ci == MPOLS) poly = 1;
		continue;
	    }
	    else if (command->ci == VAR) {
		n += strlen(field) + 1;
		command->list[j] = LISTSEP;
		continue;
	    }
	    else {
		command->list[0] -= 1;
		break;
	    }
	}

	if (!isalpha((unsigned char) *field) && 
	    !isdigit((unsigned char) *field) &&
	    !(spacename && (*field == '"' || *field == '\''))) { 
	    command->errcode = 1;
	    sprintf(gretl_errmsg, 
		    _("field '%s' in command is invalid"), field);
	    free(remainder);
	    return;
	}

	/* check command->list for scalars */
	if (!ar && !poly && command->ci != PRINT && command->ci != STORE) {
	    if (!pdinfo->vector[command->list[j]]) {
		command->errcode = 1;
		sprintf(gretl_errmsg, 
			_("variable %s is a scalar"), field);
		free(remainder);
		return;
	    }
	}

	n += strlen(field) + 1;
    } /* end of loop through fields in command line */

    /* commands that can take a specified list, but where if the
       list is null or just ";" we want to operate on all variables
    */    
    if (command->ci == PRINT ||
	command->ci == STORE || 
	command->ci == CORR || 
	command->ci == LOGS ||
	command->ci == SQUARE ||
	command->ci == DIFF ||
	command->ci == LDIFF ||
	command->ci == SUMMARY ||
	command->ci == SMPL ||
	command->ci == PCA ||
	command->ci == LAGS) {
	if (command->list[0] == 0) {
	    _full_list(pdinfo, command);
	    /* suppress echo of the list -- may be too long */
	    command->nolist = 1;
	}
    } else if (command->ci != SETMISS && command->ci != DELEET) {
	/* command that needs a list but doesn't have one */
	if (command->list[0] == 0) command->errcode = E_ARGS;
    }

    if (NEEDS_TWO_VARS(command->ci) && command->list[0] == 1)
	command->errcode = E_ARGS;

    if ((command->ci == AR || command->ci == TSLS || 
	 command->ci == ARMA || command->ci == SCATTERS) 
	&& strlen(command->param) == 0) {
	command->errcode = E_ARGS;
    }

    free(remainder);
}

static void nl_strip (char *line)
{
    int n = strlen(line);

    if (n && line[n-1] == '\n') line[n-1] = 0;
}

/**
 * help:
 * @cmd: the command on which help is wanted.
 * @helpfile: path to the gretl help file.
 * @prn: pointer to gretl printing struct.
 *
 * Searches in @helpfile for help on @cmd and, if help is found,
 * prints it to @prn.  If @cmd is a NULL pointer, lists the valid
 * commands.
 *
 * Returns: 0 on success, 1 if the helpfile was not found or the
 * requested topic was not found.
 */

int help (const char *cmd, const char *helpfile, PRN *prn)
{
    FILE *fq;
    char line[MAXLEN], cmdcopy[9];
    int i, ok;

    if (cmd == NULL) {
	pputs(prn, _("\nValid gretl commands are:\n"));
	for (i=1; i<NC; i++) {
	    pprintf(prn, "%-9s", gretl_command_word(i));
	    if (i%8 == 0) pputs(prn, "\n");
	    else pputs(prn, " ");
	}
	pputs(prn, _("\n\nFor help on a specific command, type: help cmdname"));
	pputs(prn, _(" (e.g. help smpl)\n"));
	return 0;
    }

    *cmdcopy = 0;
    strncat(cmdcopy, cmd, 8);

    ok = 0;

    if (gretl_command_number(cmd) > 0) ok = 1;

    if (!ok && aliased(cmdcopy)) {
	if (gretl_command_number(cmdcopy) > 0) {
	    ok = 1;
	}
    }

    if (!ok) {
	pprintf(prn, _("\"%s\" is not a gretl command.\n"), cmd);
	return 1;
    }

    if ((fq = fopen(helpfile, "r")) == NULL) {
	printf(_("Unable to access the file %s.\n"), helpfile);
	return 1;
    } 

    ok = 0;
    while (fgets(line, MAXLEN, fq) != NULL) {
	nl_strip(line);
	if (strcmp(cmdcopy, line) == 0) {
	    ok = 1;
	    pputs(prn, "\n");
	    while (fgets(line, MAXLEN, fq)) {
		if (*line == '#') break;
		nl_strip(line);
		if (*line != '@') {
		    pprintf(prn, "%s\n", line);
		}		
	    }
	    break;
	}
    }

    if (!ok) {
	pprintf(prn, _("%s: sorry, no help available.\n"), cmd);
    }

    fclose(fq);
    return 0;
}

/* ........................................................... */

static int parse_criteria (const char *line, const DATAINFO *pdinfo, 
			   double ***pZ, PRN *prn)
{
    double ess;
    int i, T, k;
    char cmd[9], essstr[32], Tstr[9], kstr[9];
    
    if (sscanf(line, "%s %s %s %s", cmd, essstr, Tstr, kstr) != 4) {
	return 1;
    }

    if (isalpha((unsigned char) *essstr) && 
	(i = varindex(pdinfo, essstr)) < pdinfo->v) 
	    ess = get_xvalue(i, *pZ, pdinfo);
    else if (isdigit(*essstr)) ess = atof(essstr);
    else return 1;
    if (ess < 0) {
	pputs(prn, _("ess: negative value is out of bounds.\n"));
	return 1;
    }

    if (isalpha((unsigned char) *Tstr) &&
	(i = varindex(pdinfo, Tstr)) < pdinfo->v) 
	    T = (int) get_xvalue(i, *pZ, pdinfo);
    else if (isdigit(*Tstr)) T = atoi(Tstr);
    else return 1;
    if (T < 0) {
	pputs(prn, _("T: negative value is out of bounds.\n"));
	return 1;
    }

    if (isalpha((unsigned char) *kstr) &&
	(i = varindex(pdinfo, kstr)) < pdinfo->v) 
	    k = (int) get_xvalue(i, *pZ, pdinfo);
    else if (isdigit(*kstr)) k = atoi(kstr);
    else return 1;
    if (k < 0) {
	pputs(prn, _("k: negative value is out of bounds.\n"));
	return 1;
    }   
 
    _criteria(ess, T, k, prn);

    return 0;
}

/**
 * fcast:
 * @line: the command line, giving a starting observation, ending
 * observation, and variable name to use for the forecast values
 * (the starting and ending observations may be omitted).
 * @pmod: pointer to gretl #MODEL.
 * @pdinfo: pointer to data information struct.
 * @pZ: pointer to data matrix.
 *
 * Creates a new variable containing predicted values for the
 * dependent variable in @pmod.
 *
 * Returns: the ID number of the newly created variable containing the
 * forecast, or a negative integer on error.
 */

int fcast (const char *line, const MODEL *pmod, DATAINFO *pdinfo, 
	   double ***pZ)
{
    int t, t1, t2, vi;
    char t1str[OBSLEN], t2str[OBSLEN], varname[VNAMELEN];

    *t1str = '\0'; *t2str = '\0';

    /* the varname should either be in the 2nd or 4th position */
    if (sscanf(line, "%*s %8s %8s %8s", t1str, t2str, varname) != 3) {
	if (sscanf(line, "%*s" "%8s", varname) != 1) return -1;
    }

    if (*t1str && *t2str) {
	t1 = dateton(t1str, pdinfo);
	t2 = dateton(t2str, pdinfo);
	if (t1 < 0 || t2 < 0 || t2 < t1) return -1;
    } else {
	t1 = pdinfo->t1;
	t2 = pdinfo->t2;
    }

    if (!isalpha((unsigned char) *varname)) {
        sprintf(gretl_errmsg, _("First char of varname ('%c') is bad\n"
               "(first must be alphabetical)"), varname[0]);
	return -1;
    }

    if (_reserved(varname)) return -1;

    vi = varindex(pdinfo, varname);

    if (vi >= pdinfo->v && dataset_add_vars(1, pZ, pdinfo)) 
	return -1 * E_ALLOC;

    strcpy(pdinfo->varname[vi], varname);
    strcpy(VARLABEL(pdinfo, vi), _("predicted values"));

    for (t=0; t<pdinfo->n; t++) (*pZ)[vi][t] = NADBL;

    _forecast(t1, t2, vi, pmod, pZ);

    return vi;
}

/* ........................................................... */

static int _full_list (const DATAINFO *pdinfo, CMD *command)
/* create a gretl "list" containing all the vars in the data set,
   except for the constant and any "hidden" variables or scalars */
{
    int i, n = 1;

    command->list = realloc(command->list, pdinfo->v * sizeof(int));
    if (command->list == NULL) return E_ALLOC;

    for (i=1; i<pdinfo->v; i++) {
	if (hidden_var(i, pdinfo)) continue;
	if (pdinfo->vector[i] == 0) continue;
	command->list[n++] = i;
    }
    command->list[0] = n - 1;
    return 0;
}

/**
 * parseopt:
 * @s: option string, as supplied on the command line.
 *
 * Returns: the gretl option code correspoding to @s, or 0 if the option
 * string is not recognized.
 */

int parseopt (const char *s)
{
    if (strcmp(s, "-b") == 0 || strncmp(s, "--batch", 7) == 0) 
	return OPT_BATCH;
    if (strcmp(s, "-h") == 0 || strcmp(s, "--help") == 0) 
	return OPT_HELP;
    if (strcmp(s, "-p") == 0 || strcmp(s, "--pvalue") == 0) 
	return OPT_PVALS;
    if (strcmp(s, "-v") == 0 || strcmp(s, "--version") == 0) 
	return OPT_VERSION;
    if (strcmp(s, "-r") == 0 || strncmp(s, "--run", 5) == 0) 
	return OPT_RUNIT;
    if (strcmp(s, "-d") == 0 || strncmp(s, "--db", 4) == 0) 
	return OPT_DBOPEN;
    if (strcmp(s, "-w") == 0 || strncmp(s, "--webdb", 7) == 0) 
	return OPT_WEBDB;
    return 0;
}

#ifndef WIN32

int shell (const char *arg)
{
    int pid;
    void (*old1)(int);
    void (*old2)(int);
    char shellnam[40];
    const char *theshell, *namep; 

    old1 = signal (SIGINT, SIG_IGN);
    old2 = signal (SIGQUIT, SIG_IGN);

    if ((pid = fork()) == 0) {
	for (pid = 3; pid < 20; pid++)
	    (void) close(pid);
	(void) signal(SIGINT, SIG_DFL);
	(void) signal(SIGQUIT, SIG_DFL);
	theshell = getenv("SHELL");
	if (theshell == NULL)
#ifdef HAVE_PATHS_H
	    theshell =_PATH_BSHELL;
#else
	    theshell = "/bin/sh"; 
#endif
	namep = strrchr(theshell, '/');
	if (namep == NULL)
	    namep = theshell;
	(void) strcpy(shellnam,"-");
	(void) strcat(shellnam, ++namep);
	if (strcmp(namep, "sh") != 0)
	    shellnam[0] = '+';
	if (arg) {
	    execl(theshell, shellnam, "-c", arg, NULL);
	}
	else {
	    execl(theshell, shellnam, NULL);
	}
	perror(theshell);
	return 1;
    }

    if (pid > 0) while (wait(NULL) != pid);

    (void) signal(SIGINT, old1);
    (void) signal(SIGQUIT, old2);

    if (pid == -1) {
	perror(_("Try again later"));
    }
    return 0;
}

#endif /* ! WIN32 */

/**
 * echo_cmd:
 * @pcmd: pointer to #CMD struct.
 * @pdinfo: pointer to data information struct.
 * @line: "raw" command line to be echoed.
 * @batch: set to 1 for batch mode, 0 for interactive.
 * @gui: 1 for the gretl GUI, 0 for command-line program.
 * @prn: pointer to gretl printing struct.
 *
 * Echoes the user command represented by @pcmd and @line.
 * 
 */

#define hold_param(c) (c == TSLS || c == AR || c == ARMA || c == CORRGM || \
                       c == MPOLS || c == SCATTERS || c == GNUPLOT || \
                       c == LOGISTIC)

void echo_cmd (CMD *cmd, const DATAINFO *pdinfo, const char *line, 
	       int batch, int gui, PRN *prn)
     /* echo a given command: depending on the circumstances, either
	to stdout or to a buffer, or both */

{
    int i, err, gotsep = 1;
    int cli = !gui;

    if (line == NULL) return;

#if 0
    fprintf(stderr, "echo_cmd: line='%s', gui=%d, cmd->opt=%ld, batch=%d, "
	    "param='%s', nolist=%d\n", line, gui, cmd->opt, batch, cmd->param,
	    cmd->nolist);
    fprintf(stderr, "cmd->cmd='%s'\n", cmd->cmd);
#endif

    /* special case: gui "store" command, which could overflow the
       "line" length; also I'm not sure whether we should record
       gui "store" in the command script; we'll record it, but
       commented out.
    */
    if (gui && !batch && cmd->ci == STORE) {  /* FIXME monte carlo loop? */
	pprintf(prn, "# store '%s'", cmd->param);
	if (cmd->opt) { 
	    const char *flagstr = print_flags(cmd->opt, cmd->ci);

	    pprintf(prn, "%s", flagstr);
	}
	pputc(prn, '\n');
	return;
    }

    if (*line == '\0' || *line == '!' || !strcmp(line, "quit"))
	return;

    if (cmd->ci == AR || cmd->ci == ARMA) gotsep = 0;

    /* command is preceded by a "savename" to which a object will
       be assigned */
    if (*cmd->savename && gui && !batch) {
	pprintf(prn, "%s <- ", cmd->savename);
    }

    if (!cmd->nolist) { 
	/* command has a list of args to be printed */
	if (cli) {
	    if (batch) {
		printf("\n? %s", cmd->cmd);
	    } else {
		printf(" %s", cmd->cmd);
	    }
	    if (cmd->ci == RHODIFF) printf(" %s;", cmd->param);
	    else if (*cmd->param && !hold_param(cmd->ci)) {
		printf(" %s", cmd->param);
	    }
	}
	if (!batch) {
	    pprintf(prn, "%s", cmd->cmd);
	    if (cmd->ci == RHODIFF) pprintf(prn, " %s;", cmd->param);
	    else if (*cmd->param && !hold_param(cmd->ci)) {
		pprintf(prn, " %s", cmd->param);
	    }
	}

	/* if list is very long, break it up over lines */
	if (cmd->ci == STORE) {
	    if (cli) printf(" \\\n");
	    if (!batch) pputs(prn, " \\\n");
	}
	for (i=1; i<=cmd->list[0]; i++) {
	    if (cmd->list[i] == LISTSEP) {
		if (cli) printf(" ;");
		if (!batch) pputs(prn, " ;");
		gotsep = (cmd->ci != MPOLS)? 1 : 0;
		continue;
	    }
	    if (cli) {
		if (gotsep) {
		    printf(" %s", pdinfo->varname[cmd->list[i]]);
		} else {
		    printf(" %d", cmd->list[i]);
		}
		if (i > 1 && i < cmd->list[0] && (i+1) % 10 == 0) {
		    printf(" \\\n"); /* line continuation */
		}
	    }
	    if (!batch) {
		if (gotsep) {
		    pprintf(prn, " %s", pdinfo->varname[cmd->list[i]]);
		} else {
		    pprintf(prn, " %d", cmd->list[i]);
		}
		if (i > 1 && i < cmd->list[0] && (i+1) % 10 == 0) {
		    pputs(prn, " \\\n"); /* line continuation */
		}
	    }
	}

	/* corrgm and gnuplot: param comes last */
	if ((cmd->ci == CORRGM || cmd->ci == GNUPLOT || cmd->ci == LOGISTIC)
	    && *cmd->param) { 
	    if (cli) printf(" %s", cmd->param);
	    if (!batch) pprintf(prn, " %s", cmd->param);
	}

	/* check for duplicated vars */
	err = _list_dups(cmd->list, cmd->ci);
	if (err) {
	    printf(_("\nvar number %d duplicated in the command list.\n"),
		   err);
	    cmd->ci = VARDUP;
	}
    } /* end if !cmd->nolist */

    else if (strcmp(cmd->cmd, "quit")) {
	if (cli) {
	    if (batch) printf("? %s", line);
	    else printf(" %s", line);
	}
	if (!batch) pputs(prn, line);
    }

    if (cmd->opt) {
	const char *flagstr;
	int ci = cmd->ci;

	if (ci == END && !strcmp(cmd->param, "nls")) {
	    ci = NLS;
	}
	flagstr = print_flags(cmd->opt, ci);
	if (cli) fputs(flagstr, stdout);
	if (!batch) pputs(prn, flagstr);
    }

    if (cli) putchar('\n');

    if (!batch) {
	pputc(prn, '\n');
	if (prn != NULL && prn->fp) {
	    fflush(prn->fp);
	}
    }
}

/* .......................................................... */

/* Look for a flag of the form "-x".  Make sure it's outside of
   any quotes.  Return pointer to flag. */

static const char *flag_present (const char *s, char f, int *quoted)
{
    int inquote = 0;
    int gotdash = 0;

    while (*s) {
	if (*s == '"') inquote = !inquote;
	if (!inquote) {
	    if (*s == '-') gotdash = 1;
	    else if (gotdash && *s == f && *(s+1)) {
#if 0
		/* blank out the flag and following in the
		   original string? */
		*s = 0;
#endif
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
	    }
	    else gotdash = 0;
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

/* .......................................................... */

static int make_var_label (const char *line, const DATAINFO *pdinfo, 
			   PRN *prn)
{
    char *p;
    char vname[VNAMELEN];
    int v;

    if (pdinfo->varinfo == NULL) return 1;

    if (sscanf(line, "label %8s", vname) != 1) return E_PARSE;

    v = varindex(pdinfo, vname);
    if (v == pdinfo->v) {
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), vname);
	return E_UNKVAR;
    }

    p = get_flag_field(line + 6, 'd');
    if (p != NULL) {
	*VARLABEL(pdinfo, v) = 0;
	strncat(VARLABEL(pdinfo, v), p, MAXLABEL - 1);
	free(p);
    }

    p = get_flag_field(line + 6, 'n');
    if (p != NULL) {
	*DISPLAYNAME(pdinfo, v) = 0;
	strncat(DISPLAYNAME(pdinfo, v), p, MAXDISP - 1);
	free(p);
    }  

    return 0;
}

/* .......................................................... */

static void get_optional_filename (const char *line, CMD *cmd)
{
    char *p;

    p = get_flag_field(line + 8, 'f');
    if (p != NULL) {
	free(cmd->param);
	cmd->param = p;
    }    
}

/* .......................................................... */

static void showlabels (const DATAINFO *pdinfo, PRN *prn)
{
    int i;

    pprintf(prn, _("Listing labels for variables:\n"));
    for (i=0; i<pdinfo->v; i++) {
	if (strlen(VARLABEL(pdinfo, i)) > 2) {
	    pprintf(prn, "%3d) %-10s %s\n", i, 
		    pdinfo->varname[i], VARLABEL(pdinfo, i));
	}
    }
}

static void do_print_string (char *str, PRN *prn)
{
    size_t len;

    if (*str == '"') str++;
    len = strlen(str);
    if (str[len-1] == '"') str[len-1] = 0;

    pprintf(prn, "%s\n", str);
}

static int do_outfile_command (unsigned long flag, char *fname,
			       PRN *prn)
{
    static char outname[MAXLEN];
    int output_diverted = 0;

    if (flag != OPT_W && flag != OPT_A && flag != OPT_C) {
	return E_ARGS;
    }

    if (prn->fpaux != NULL ||
	(prn->fp != NULL && prn->buf != NULL)) {
	output_diverted = 1;
    }

    /* command to close outfile */
    if (flag == OPT_C) {
	if (!output_diverted) {
	    pputs(prn, _("Output is not currently diverted to file\n"));
	    return 1;
	} else {
	    fclose(prn->fp);
	    prn->fp = prn->fpaux;
	    prn->fpaux = NULL;
	    pprintf(prn, "Closed output file '%s'\n", outname);
	    return 0;
	}
    }

    /* command to divert output to file */
    if (output_diverted) {
	fprintf(stderr, _("Output is already diverted to '%s'\n"),
		outname);
	return 1;
    } else {
	if (*fname == 0) {
	    return E_ARGS;
	} else {
	    FILE *fp;

	    if (flag == OPT_W) {
		fp = fopen(fname, "w");
	    } else {
		fp = fopen(fname, "a");
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

	    /* save the prn stream */
	    prn->fpaux = prn->fp;
	    /* hook output to specified file */
	    prn->fp = fp;
	    strcpy(outname, fname);
	    return 0;
	}
    }

    return 1; /* not reached */
}

/* ........................................................ */

int call_pca_plugin (CORRMAT *corrmat, double ***pZ,
		     DATAINFO *pdinfo, unsigned long *pflag,
		     PRN *prn)
{
    void *handle = NULL;
    int (*pca_from_corrmat) (CORRMAT *, double ***, DATAINFO *,
			     unsigned long *, PRN *);
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

/* ........................................................ */

int simple_commands (CMD *cmd, const char *line, 
		     double ***pZ, DATAINFO *datainfo, PATHS *paths,
		     int pause, PRN *prn)
     /* common code for command-line and GUI client programs, where
	the command doesn't require special handling on the client
	side */
{
    int err = 0, order = 0;
    CORRMAT *corrmat;
    GRETLSUMMARY *summ;

    switch (cmd->ci) {

    case ADF:
	if (!isdigit(*cmd->param)) {
	    pputs(prn, _("adf: lag order must be given first\n"));
	    break;
	}
	/* flag an error if cmd->list[0] > 1? */
	order = atoi(cmd->param);
	err = adf_test(order, cmd->list[1], pZ, datainfo, prn);
	break;

    case COINT:
	order = atoi(cmd->param);
	err = coint(order, cmd->list, pZ, datainfo, prn);
	break;

    case COINT2:
	order = atoi(cmd->param);
	err = johansen_test(order, cmd->list, pZ, datainfo, cmd->opt, prn);
	break;

    case CORR:
	if (cmd->list[0] > 3) {
	    err = esl_corrmx(cmd->list, pZ, datainfo, pause, prn);
	    if (err) 
		pputs(prn, _("Error in generating correlation matrix\n"));
	    break;
	}
	corrmat = corrlist(cmd->list, pZ, datainfo);
	if (corrmat == NULL) 
	    pputs(prn, _("Couldn't allocate memory for correlation matrix.\n"));
	else printcorr(corrmat, datainfo, prn);
	free_corrmat(corrmat);
	break;

    case PCA:
	corrmat = corrlist(cmd->list, pZ, datainfo);
	if (corrmat == NULL) {
	    pputs(prn, _("Couldn't allocate memory for correlation matrix.\n"));
	} else {
	    err = call_pca_plugin(corrmat, pZ, datainfo, &cmd->opt, prn);
	    if (cmd->opt && !err) {
		varlist(datainfo, prn);
	    }
	    free_corrmat(corrmat);
	}
	break;

    case CRITERIA:
	err = parse_criteria(line, datainfo, pZ, prn);
	if (err) 
	    pputs(prn, _("Error in computing model selection criteria.\n"));
	break;

    case CRITICAL:
	err = print_critical(line, prn);
	break;

    case DATA:
	err = db_get_series(line, pZ, datainfo, prn);
	break;

    case DIFF:
	err = list_diffgenr(cmd->list, pZ, datainfo);
	if (err) 
	    pputs(prn, _("Error adding first differences of variables.\n"));
	else varlist(datainfo, prn);
	break;

    case LDIFF:
	err = list_ldiffgenr(cmd->list, pZ, datainfo);
	if (err) 
	    pputs(prn, _("Error adding log differences of variables.\n"));
	else varlist(datainfo, prn);
	break;

    case LAGS:
	err = lags(cmd->list, pZ, datainfo); 
	if (err) 
	    pputs(prn, _("Error adding lags of variables.\n"));
	else varlist(datainfo, prn);
	break;

    case LOGS:
	err = logs(cmd->list, pZ, datainfo);
	if (err < cmd->list[0]) 
	    pputs(prn, _("Error adding logs of variables.\n"));
	if (err > 0) { 
	    varlist(datainfo, prn);
	    err = 0;
	} else
	    err = 1;
	break;

    case MULTIPLY:
	err = _multiply(cmd->param, cmd->list, cmd->str, pZ, datainfo);
	if (!err) varlist(datainfo, prn);
	break;

    case GRAPH:
	graph(cmd->list, *pZ, datainfo, cmd->opt, prn);
	break;

    case PLOT:
	plot(cmd->list, *pZ, datainfo, cmd->opt, pause, prn);
	break;

    case RMPLOT:
	if (cmd->list[0] != 1) {
	    pputs(prn, _("This command requires one variable.\n"));
	    err = 1;
	} else {
	    err = rmplot(cmd->list, *pZ, datainfo, prn, paths);
	}
	break;

    case INFO:
	if (datainfo->descrip != NULL) 
	    pprintf(prn, "%s\n", datainfo->descrip);
	else 
	    pputs(prn, _("No data information is available.\n"));
	break;

    case LABEL:
	err = make_var_label(line, datainfo, prn);
	break;

    case LABELS:
	showlabels(datainfo, prn);
	break;

    case VARLIST:
	varlist(datainfo, prn);
	break;

    case PRINT:
	if (strlen(cmd->param)) {
	    do_print_string(cmd->param, prn);
	} else {
	    printdata(cmd->list, pZ, datainfo, pause, cmd->opt, prn);
	}
	break;

    case SUMMARY:
	summ = summary(cmd->list, pZ, datainfo, prn);
	if (summ == NULL) 
	    pputs(prn, _("generation of summary stats failed\n"));
	else {
	    print_summary(summ, datainfo, pause, prn);
	    free_summary(summ);
	}
	break; 

    case MEANTEST:
	err = means_test(cmd->list, *pZ, datainfo, cmd->opt, prn);
	break;	

    case VARTEST:
	err = vars_test(cmd->list, *pZ, datainfo, prn);
	break;

    case RUNS:
	err = runs_test(cmd->list[1], *pZ, datainfo, prn);
	break;

    case SPEARMAN:
	err = spearman(cmd->list, *pZ, datainfo, cmd->opt, prn);
	break;

    case OUTFILE:
	err = do_outfile_command(cmd->opt, cmd->param, prn);
	break;

    default:
	break;
    }

    return err;
}

int ready_for_command (const char *line)
{
    const char *ok_cmds[] = {
	"open", 
	"run", 
	"nulldata", 
	"import", 
	"pvalue",
	"eval",
	"!",
	"(*", 
	"man ", 
	"help", 
	"set", 
	"critical", 
	"seed", 
	"genr",
	NULL 
    };
    const char **p = ok_cmds;

    if (string_is_blank(line)) return 1;

    if (*line == 'q' || *line == 'x' || 
	*line == '\0' || *line == '#') return 1;

    while (*p) {
	if (strncmp(line, *p, strlen(*p)) == 0)
	    return 1;
	p++;
    }

    return 0;
}




