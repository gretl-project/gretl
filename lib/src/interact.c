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
#include "internal.h"

/* equipment for the "shell" command */
#ifndef OS_WIN32
# include <sys/wait.h>
# include <signal.h>
# include <errno.h>
# include <unistd.h>
# ifdef HAVE_PATHS_H
#  include <paths.h>
# endif
#endif

#include "cmdlist.h"

extern int _omitfromlist (int *list, const int *omitvars, int newlist[],
			  const DATAINFO *pdinfo, const int model_count);
extern int _parse_lagvar (const char *varname, LAGVAR *plagv, 
			  DATAINFO *pdinfo);

static int _full_list (const DATAINFO *pdinfo, CMD *command);


/* ........................................................... */

static int trydatafile (char *line, int *ignore)
{
    int i, m, n = strlen(line);
    char datfile[MAXLEN];

    datfile[0] = '\0';
    for (i=0; i<n; i++) {
	if ((n - i) > 4 && strncmp(line+i, "DATA", 4) == 0) {
	    sscanf(line+i, "%s", datfile);
	    m = strlen(datfile);
	    if (datfile[m-1] == ',') datfile[m-1] = '\0';
	    lower(datfile);
	    i += 4;
	}
	else if (line[i] == '*' && line[i+1] == ')') *ignore = 0;
    }
    if (datfile[0] != '\0') {
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
		sscanf(line+4, "%s", datfile);
		sprintf(line, "open %s", datfile);
		*ignore = 0;  /* FIXME ? */
		return 0;
	    }
	}
	else if (line[i] == '*' && line [i+1] == ')') {
	    *ignore = 0; i += 2;
	    while (isspace((unsigned char) line[i]) && i < n) i++;
	}
	if (!(*ignore) && line[i] != 13) {
	    tmpstr[j] = line[i];
	    j++;
	}
    }
    tmpstr[j] = '\0';
    strcpy(line, tmpstr);
    if (strlen(line)) return 0;
    else return 1;
}

/* ........................................................... */

static int get_rhodiff_param (char *str, CMD *cmd)
{
    int k;

    if ((k = haschar(';', str)) < 0) return 1;
    cmd->param = realloc(cmd->param, k+1);
    if (cmd->param == NULL) return E_ALLOC;
    strncpy(cmd->param, str, k);
    cmd->param[k] = '\0';
    _shiftleft(str, k + 1);
    return 0;
}

/* ........................................................... */

int command_number (const char *cmd)
{    
    int i;

    for (i=0; i<NC; i++) 
	if (strcmp(cmd, commands[i]) == 0) 
	    return i;
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
    else if (cmd[0] == '!') {
	strcpy(cmd, "shell");
	return 1;
    }
    return 0;
}

#define NO_VARLIST(c) (c == VARLIST || \
	               c == NOECHO || \
	               c == QUIT || \
	               c == SMPL || \
	               c == EQNPRINT || \
	               c == TABPRINT || \
	               c == FCAST || \
	               c == FCASTERR || \
	               c == FIT || \
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
	               c == DELEET || \
	               c == TESTUHAT || \
	               c == GENR)

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

    command->errcode = 0;
    gretl_errmsg[0] = '\0';
    command->nolist = 0;
    command->param[0] = '\0';

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

    linelen = strlen(line);

    if (line[0] == '#' || sscanf(line, "%s", command->cmd) != 1) {
	command->nolist = 1;
	command->ci = -1;
	return;
    }

    /* command aliases */
    if (aliased(command->cmd) == 2) {
	command->param = realloc(command->param, 2);
	strcpy(command->param, "x");
    }

    /* trap bogus commands */    
    if ((command->ci = command_number(command->cmd)) == 0) {
	command->errcode = 1;
	sprintf(gretl_errmsg, _("command '%s' not recognized"), 
		command->cmd);
	return;
    }

    /* if, else, endif controls */
    if (flow_control(line, pZ, pdinfo, command)) {
	command->nolist = 1;
	command->ci = -1;
	return;
    }

    /* commands that never take a list of variables */
    if (NO_VARLIST(command->ci)) {
	command->nolist = 1;
	return;
    }

    /* boxplots can be special: boolean conditions embedded in
       the list, which have to be parsed separately */
    if (command->ci == BXPLOT && strchr(line, '(')) {
	command->nolist = 1;
	return;
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
	command->ci == SEED ||
	command->ci == LMTEST ||
	command->ci == NULLDATA ||
	(command->ci == PRINT && strstr(line, "\""))) {
	command->nolist = 1;
	if (!strncmp(line, "man ", 4)) n--;
	if (nf) {
	    command->param = realloc(command->param, linelen - n + 1);
	    if (command->ci == PRINT)
		strcpy(command->param, line + n + 1);
	    else
		sscanf(line + n + 1, "%s", command->param);
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
	    for (i=1; i<strlen(remainder); i++) {
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
	    for (i=0; i<strlen(command->param) - 2; i++)
		command->param[i] = command->param[i+1];
	    command->param[i] = '\0';
	    strcpy(line, remainder);
	    nf--;
	    n = 0;
	    linelen = strlen(line);
	}
    } /* end if STORE && nf */

    /* "store" takes a filename before the list, "var" takes a 
       lag order, "adf" takes a lag order, "arch" takes a lag 
       order, "multiply" takes a multiplier.  "omitfrom" and
       "addto" take the ID of a previous model. "setmiss" takes
       a value to be interpreted as "missing"
    */
    if ((command->ci == STORE && !spacename) ||
	command->ci == ADF ||
	command->ci == ARCH ||
	command->ci == COINT ||
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
	    strcpy(line, remainder);
	    nf--;
	    n = 0;
	    linelen = strlen(line);
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

    if (command->ci == AR) ar = 1;

    command->list = realloc(command->list, (1 + nf) * sizeof(int));
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

	if (isalpha((unsigned char) field[0])) {
	    /* fprintf(stderr, "field: %s\n", field); */
	    if (field[strlen(field)-1] == ';')
		field[strlen(field)-1] = '\0';
	    if ((v = varindex(pdinfo, field)) <= pdinfo->v - 1) 
		command->list[j] = v;
	    else {
		/* Automated lags:
		   could be the string is like "varname(-2)", meaning
		   the second lag of varname: in this case auto-
		   generate the lag variable. */
		if (_parse_lagvar(field, &lagvar, pdinfo)) {
		    if (_laggenr(lagvar.varnum, lagvar.lag, 0, pZ, pdinfo)) {
			command->errcode = 1;
			sprintf(gretl_errmsg, 
				_("generation of lag variable failed"));
			free(remainder);
			return;
		    } else { 
			command->list[j] = pdinfo->v - 1;
			if (cmds) 
			    pprintf(cmds, "genr %s\n", 
				    pdinfo->label[pdinfo->v - 1]);
			n += strlen(field) + 1;
			continue; 
		    }
		} 
		/* auto-generation of plotting variables */
		if (!strcmp(field, "qtrs") || 
		    !strcmp(field, "months") || !strcmp(field, "time")) {
		    plotvar(pZ, pdinfo, field);
		    command->list[j] = pdinfo->v - 1;
		} else {
		    command->errcode = 1;
		    sprintf(gretl_errmsg, 
			    _("'%s' is not the name of a variable"), field);
		    free(remainder);
		    return;
		}
	    }
	} /* end if isalpha(field[0]) */

	if (isdigit(field[0])) {
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

	else if (field[0] == ';') {
	    if (command->ci == TSLS || command->ci == AR ||
		command->ci == MPOLS || command->ci == SCATTERS) {
		command->param = realloc(command->param, 4);
		sprintf(command->param, "%d", j);
		n += strlen(field) + 1;
		command->list[j] = 999;
		ar = 0; /* turn off acceptance of AR lags */
		if (command->ci == MPOLS) poly = 1;
		continue;
	    }
	    else {
		command->list[0] -= 1;
		break;
	    }
	}

	if (!isalpha((unsigned char) field[0]) && 
	    !isdigit((unsigned char) field[0]) &&
	    !(spacename && (field[0] == '"' || field[0] == '\''))) { 
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
    }

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
	command->ci == LAGS) {
	if (command->list[0] == 0) {
	    _full_list(pdinfo, command);
	    /* suppress echo of the list -- may be too long */
	    command->nolist = 1;
	}
    } else if (command->ci != SETMISS) {
	/* command that needs a list but doesn't have one */
	if (command->list[0] == 0) command->errcode = E_ARGS;
    }

    if (NEEDS_TWO_VARS(command->ci) && command->list[0] == 1)
	command->errcode = E_ARGS;

    if ((command->ci == AR || command->ci == TSLS || 
	 command->ci == SCATTERS) && strlen(command->param) == 0) {
	command->errcode = E_ARGS;
    }

    free(remainder);
    return;
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
    char line[MAXLEN], tmp[MAXLEN], cmdcopy[9];
    int i, ok;

    if (cmd == NULL) {
	pprintf(prn, _("\nValid gretl commands are:\n"));
	for (i=1; i<NC; i++) {
	    pprintf(prn, "%-9s", commands[i]);
	    if (i%8 == 0) pprintf(prn, "\n");
	    else pprintf(prn, " ");
	}
	pprintf(prn, "\n");
	pprintf(prn, _("\nFor help on a specific command, type: help cmdname"));
	pprintf(prn, _(" (e.g. help smpl)\n"));
	return 0;
    }

    strncpy(cmdcopy, cmd, 8);
    cmdcopy[8] = '\0';

    ok = 0;
    for (i=1; i<NC; i++) {
	if (!strcmp(commands[i], cmd)) {
	    ok = 1;
	    break;
	}
    }
    if (!ok && aliased(cmdcopy)) {
	for (i=1; i<NC; i++) {
	    if (!strcmp(commands[i], cmdcopy)) {
		ok = 1;
		break;
	    }
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

    while (fgets(line, MAXLEN, fq) != NULL) {
	delchar('\n', line);
	ok = !strcmp(cmdcopy, line);
	if (!ok) continue;
	pprintf(prn, "\n");
	do {
	    if (fgets(tmp, MAXLEN, fq) == NULL) {
		fclose(fq);
		return 0;
	    }
	    delchar('\n', tmp);
	    i = strcmp(tmp, "#");
	    if (!i) {
		fclose(fq);
		return 0;
	    }
	    if (*tmp != '@')
		pprintf(prn, "%s\n", tmp);
	} while (i);
	if (ok) {
	    fclose(fq);
	    return 0;
	}
    }
    pprintf(prn, _("%s: sorry, no help available.\n"), cmd);
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
    if (isalpha((unsigned char) essstr[0]) && 
	(i = varindex(pdinfo, essstr)) < pdinfo->v) 
	    ess = get_xvalue(i, *pZ, pdinfo);
    else if (isdigit(essstr[0])) ess = atof(essstr);
    else return 1;
    if (ess < 0) {
	pprintf(prn, _("ess: negative value is out of bounds.\n"));
	return 1;
    }
    if (isalpha((unsigned char) Tstr[0]) &&
	(i = varindex(pdinfo, Tstr)) < pdinfo->v) 
	    T = (int) get_xvalue(i, *pZ, pdinfo);
    else if (isdigit(Tstr[0])) T = atoi(Tstr);
    else return 1;
    if (T < 0) {
	pprintf(prn, _("T: negative value is out of bounds.\n"));
	return 1;
    }
    if (isalpha((unsigned char) kstr[0]) &&
	(i = varindex(pdinfo, kstr)) < pdinfo->v) 
	    k = (int) get_xvalue(i, *pZ, pdinfo);
    else if (isdigit(kstr[0])) k = atoi(kstr);
    else return 1;
    if (k < 0) {
	pprintf(prn, _("k: negative value is out of bounds.\n"));
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
    char t1str[9], t2str[9], varname[9];

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

    if (!isalpha((unsigned char) varname[0])) return -1;
    varname[8] = 0;
    vi = varindex(pdinfo, varname);

    if (vi >= pdinfo->v && dataset_add_vars(1, pZ, pdinfo)) 
	return -1 * E_ALLOC;

    strcpy(pdinfo->varname[vi], varname);
    strcpy(pdinfo->label[vi], _("predicted values"));

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

#ifndef OS_WIN32

/* ........................................................ */

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

#endif

/**
 * echo_cmd:
 * @pcmd: pointer to #CMD struct.
 * @pdinfo: pointer to data information struct.
 * @line: "raw" command line to be echoed.
 * @batch: set to 1 for batch mode, 0 for interactive.
 * @gui: 1 for the gretl GUI, 0 for command-line program.
 * @oflag:
 * @prn: pointer to gretl printing struct.
 *
 * Echoes the use command represented by @pcmd and @line.
 * 
 */

void echo_cmd (CMD *pcmd, const DATAINFO *pdinfo, const char *line, 
	       int batch, int gui, int oflag, PRN *prn)
     /* echo a given command: depending on the circumstances, either
	to stdout or to a buffer, or both */

{
    int i, err, got999 = 1;
    char flagc;

    if (strcmp(line, "quit") == 0 || line[0] == '!' ||
	strlen(line) == 0) return;

    if (pcmd->ci == AR) got999 = 0;

    if (!pcmd->nolist) { /* print list of params to command */
	if (!gui) {
	    if (!batch) printf(" %s", pcmd->cmd);
	    else printf("\n? %s", pcmd->cmd);
	    if (pcmd->ci == RHODIFF) printf(" %s;", pcmd->param);
	    else if (strlen(pcmd->param) && pcmd->ci != TSLS 
		     && pcmd->ci != AR && pcmd->ci != CORRGM
		     && pcmd->ci != MPOLS && pcmd->ci != SCATTERS) 
		printf(" %s", pcmd->param);
	}
	if (!batch) {
	    pprintf(prn, "%s", pcmd->cmd);
	    if (pcmd->ci == RHODIFF) pprintf(prn, " %s;", pcmd->param);
	    else if (strlen(pcmd->param) && pcmd->ci != TSLS 
		     && pcmd->ci != AR && pcmd->ci != CORRGM
		     && pcmd->ci != MPOLS && pcmd->ci != SCATTERS) 
		pprintf(prn, " %s", pcmd->param);
	}
	/* if list is very long, break it up over lines */
	if (pcmd->ci == STORE) {
	    if (!gui) printf(" \\\n");
	    if (!batch) pprintf(prn, " \\\n");
	}
	for (i=1; i<=pcmd->list[0]; i++) {
	    if (pcmd->list[i] == 999) {
		if (!gui) printf(" ;");
		if (!batch) pprintf(prn, " ;");
		got999 = (pcmd->ci != MPOLS)? 1 : 0;
		continue;
	    }
	    if (!gui) {
		if (got999) 
		    printf(" %s", pdinfo->varname[pcmd->list[i]]);
		else printf(" %d", pcmd->list[i]);
		if (i > 1 && i < pcmd->list[0] && (i+1) % 10 == 0) 
		    printf(" \\\n"); /* break line */
	    }
	    if (!batch) {
		if (got999) 
		    pprintf(prn, " %s", pdinfo->varname[pcmd->list[i]]);
		else pprintf(prn, " %d", pcmd->list[i]);
		if (i > 1 && i < pcmd->list[0] && (i+1) % 10 == 0) 
		    pprintf(prn, " \\\n"); /* break line */
	    }
	}
	/* corrgm: param comes last */
	if (pcmd->ci == CORRGM && strlen(pcmd->param)) { 
	    if (!gui) printf(" %s", pcmd->param);
	    if (!batch) pprintf(prn, " %s", pcmd->param);
	}
	err = _list_dups(pcmd->list, pcmd->ci);
	if (err) {
	    printf(_("\nvar number %d duplicated in the command list.\n"),
		   err);
	    pcmd->ci = 999;
	}
    } /* end if !pcmd->nolist */
    else if (strcmp (pcmd->cmd, "quit")) {
	if (!gui) {
	    if (batch) printf("? %s", line);
	    else printf(" %s", line);
	}
	if (!batch) pprintf(prn, "%s", line);
    }
    if (oflag) { 
	flagc = getflag(oflag);
	if (!gui) printf(" -%c", flagc);
	if (!batch) pprintf(prn, " -%c", flagc);
    }
    if (!gui) putchar('\n');
    if (!batch) {
	pprintf(prn, "\n");
	if (prn != NULL && prn->fp) fflush(prn->fp);
    }
}

/* .......................................................... */

static void showlabels (const DATAINFO *pdinfo)
{
    int i;

    printf(_("Listing labels for variables:\n"));
    for (i=0; i<pdinfo->v; i++) {
	if (strlen(pdinfo->label[i]) > 2) {
	    printf("%3d) %-10s %s\n", i, 
		   pdinfo->varname[i], pdinfo->label[i]);
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

/* ........................................................ */

int simple_commands (CMD *cmd, const char *line, 
		     double ***pZ, DATAINFO *datainfo, PATHS *paths,
		     int pause, int oflag, PRN *prn)
     /* common code for command-line and GUI client programs, where
	the command doesn't require special handling on the client
	side */
{
    int err = 0, order = 0;
    CORRMAT *corrmat;
    GRETLSUMMARY *summ;

    switch (cmd->ci) {

    case ADF:
	if (!isdigit(cmd->param[0])) {
	    pprintf(prn, _("adf: lag order must be given first\n"));
	    break;
	}
	order = atoi(cmd->param);
	err = adf_test(order, cmd->list[1], pZ, datainfo, prn);
	if (err) errmsg(err, prn);
	break;

    case COINT:
	order = atoi(cmd->param);
	err = coint(order, cmd->list, pZ, datainfo, prn);
	break;

    case CORR:
	if (cmd->list[0] > 3) {
	    err = esl_corrmx(cmd->list, pZ, datainfo, pause, prn);
	    if (err) 
		pprintf(prn, _("Error in generating correlation matrix\n"));
	    break;
	}
	corrmat = corrlist(cmd->list, pZ, datainfo);
	if (corrmat == NULL) 
	    pprintf(prn, _("Couldn't allocate memory for correlation matrix.\n"));
	else printcorr(corrmat, datainfo, prn);
	free_corrmat(corrmat);
	break;

    case CRITERIA:
	err = parse_criteria(line, datainfo, pZ, prn);
	if (err) 
	    pprintf(prn, _("Error in computing model selection criteria.\n"));
	break;

    case CRITICAL:
	err = print_critical(line, prn);
	break;

    case DIFF:
	err = list_diffgenr(cmd->list, pZ, datainfo);
	if (err) 
	    pprintf(prn, _("Error adding first differences of variables.\n"));
	else varlist(datainfo, prn);
	break;

    case LDIFF:
	err = list_ldiffgenr(cmd->list, pZ, datainfo);
	if (err) 
	    pprintf(prn, _("Error adding log differences of variables.\n"));
	else varlist(datainfo, prn);
	break;

    case LAGS:
	err = lags(cmd->list, pZ, datainfo); 
	if (err) 
	    pprintf(prn, _("Error adding lags of variables.\n"));
	else varlist(datainfo, prn);
	break;

    case LOGS:
	err = logs(cmd->list, pZ, datainfo);
	if (err < cmd->list[0]) 
	    pprintf(prn, _("Error adding logs of variables.\n"));
	if (err > 0) { 
	    varlist(datainfo, prn);
	    err = 0;
	} else
	    err = 1;
	break;

    case MULTIPLY:
	err = _multiply(cmd->param, cmd->list, cmd->str, pZ, datainfo);
	if (err) errmsg(err, prn);
	else varlist(datainfo, prn);
	break;

    case GRAPH:
	graph(cmd->list, *pZ, datainfo, oflag, prn);
	break;

    case PLOT:
	plot(cmd->list, *pZ, datainfo, oflag, pause, prn);
	break;

    case RMPLOT:
	if (cmd->list[0] != 1) {
	    pprintf(prn, _("This command requires one variable.\n"));
	    err = 1;
	} else {
	    err = rmplot(cmd->list, *pZ, datainfo, prn, paths);
	    if (err) errmsg(err, prn);
	}
	break;

    case INFO:
	if (datainfo->descrip != NULL) 
	    pprintf(prn, "%s\n", datainfo->descrip);
	else 
	    pprintf(prn, _("No data information is available.\n"));
	break;

    case LABELS:
	showlabels(datainfo);
	break;

    case VARLIST:
	varlist(datainfo, prn);
	break;

    case PRINT:
	if (strlen(cmd->param)) 
	    do_print_string(cmd->param, prn);
	else 
	    printdata(cmd->list, pZ, datainfo, pause, oflag, prn);
	break;

    case SUMMARY:
	summ = summary(cmd->list, pZ, datainfo, prn);
	if (summ == NULL) 
	    pprintf(prn, _("generation of summary stats failed\n"));
	else {
	    print_summary(summ, datainfo, pause, prn);
	    free_summary(summ);
	}
	break; 

    case MEANTEST:
	err = means_test(cmd->list, *pZ, datainfo, oflag, prn);
	if (err) errmsg(err, prn);
	break;	

    case VARTEST:
	err = vars_test(cmd->list, *pZ, datainfo, prn);
	if (err) errmsg(err, prn);
	break;

    case RUNS:
	err = runs_test(cmd->list[1], *pZ, datainfo, prn);
	break;

    case SPEARMAN:
	err = spearman(cmd->list, *pZ, datainfo, oflag, prn);
	break;

    default:
	break;
    }

    return err;
}





