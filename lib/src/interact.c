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
#include <sys/wait.h>
#include <signal.h>
#include <errno.h>
#include <unistd.h>
#ifdef HAVE_PATHS_H
#include <paths.h>
#endif
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
	sprintf(line, "open %s", datfile);
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

    /*  printf("get_rhodiff_param: str = %s\n", str); */
    if ((k = haschar(';', str)) < 0) return 1;
    /*  printf("get_rhodiff_param: k = %d\n", k); */
    cmd->param = realloc(cmd->param, k+1);
    if (cmd->param == NULL) return E_ALLOC;
    strncpy(cmd->param, str, k);
    cmd->param[k] = '\0';
    /*  printf("get_rhodiff_param: param = %s\n", cmd->param); */
    shiftleft(str, k + 1);
    return 0;
}

/* ........................................................... */

int _command_number (const char *cmd)
{    
    int i;

    for (i=0; i<NC; i++) 
	if (strcmp(cmd, commands[i]) == 0) 
	    return i;
    return 0;
}

/* ........................................................... */

void getcmd (char *line, DATAINFO *pdinfo, CMD *command, 
	     int *ignore, double **pZ, print_t *cmds)
{
    int i, j, nf, linelen, n, v, gotdata = 0, ar = 0;
    char field[10], *remainder, *tmpstr;
    LAGVAR lagvar;

    command->errcode = 0;
    strcpy(command->errmsg, "");
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
    if (strcmp(command->cmd, "q") == 0) 
	strcpy(command->cmd, "quit");
    if (strcmp(command->cmd, "x") == 0) {
	strcpy(command->cmd, "quit");
	command->param = realloc(command->param, 2);
	strcpy(command->param, "x");
    }
    if (strcmp(command->cmd, "ls") == 0) 
	strcpy(command->cmd, "list");
    if (strcmp(command->cmd, "man") == 0) 
	strcpy(command->cmd, "help");
    if (strcmp(command->cmd, "sample") == 0) 
	strcpy(command->cmd, "smpl");
    if (command->cmd[0] == '!') 
	strcpy(command->cmd, "shell");

    /* trap bogus commands */    
    if ((command->ci = _command_number(command->cmd)) == 0) {
	command->errcode = 1;
	sprintf(command->errmsg, "command '%s' not recognized", 
		command->cmd);
	return;
    }

    /* commands that never take a list of variables */
    if (command->ci == LIST ||
	command->ci == QUIT || 
	command->ci == SMPL ||
	command->ci == EQNPRINT ||
	command->ci == TABPRINT ||
	command->ci == FCAST ||
	command->ci == FCASTERR ||
	command->ci == FIT ||
	command->ci == LABELS ||
	command->ci == INFO ||
	command->ci == LMTEST ||
	command->ci == CRITERIA ||
	command->ci == PVALUE ||
	command->ci == RUN ||
	command->ci == SHELL ||
	command->ci == SETOBS ||
	command->ci == CHOW ||
	command->ci == CUSUM ||
	command->ci == OPEN ||
	command->ci == IMPORT ||
	command->ci == ENDLOOP ||
	command->ci == SIM ||
	command->ci == DELEET ||
	command->ci == TESTUHAT ||
	command->ci == GENR) {
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
    nf = count_fields(line) - 1;
    n = strlen(command->cmd);
    remainder = malloc(linelen-n+1);

    /* ...unless it's "help", "loop", or "nulldata" 
       which are special */
    if (command->ci == HELP ||
	command->ci == LOOP ||
	command->ci == SEED ||
	command->ci == NULLDATA) {
	command->nolist = 1;
	if (nf) {
	    command->param = realloc(command->param, linelen - n + 1);
	    for (i=0; i<=linelen-n; i++) 
		remainder[i] = line[i+n]; 
	    remainder[linelen-n] = '\0'; 
	    sscanf(remainder, "%s", command->param);
	    free(remainder);    
	    return;
	} else command->param[0] = '\0';
	free(remainder);    
	return;
    }

    /* need to treat rhodiff specially -- put everything from
       the end of the command word to the first semicolon into
       "param", for processing later */
    if (command->ci == RHODIFF) {  /* FIXME */
	for (i=0; i<linelen-n; i++) 
	    remainder[i] = line[i+n+1];
	remainder[linelen-n] = '\0'; 
	if (get_rhodiff_param(remainder, command)) {
	    command->errcode = E_SYNTAX;
	    free(remainder);
	    return;
	}
	strcpy(line, remainder);
	linelen = strlen(line);
	nf = count_fields(line);
	n = 0;
    }

    /* "store" takes a filename before the list, "var" takes a 
       lag order, "adf" takes a lag order, "arch" takes a lag 
       order, "multiply" takes a multiplier.  "omitfrom" and
       "addto" take the ID of a previous model */
    if (command->ci == STORE || 
	command->ci == ADF ||
	command->ci == ARCH ||
	command->ci == COINT ||
	command->ci == ADDTO ||
	command->ci == OMITFROM ||
	command->ci == MULTIPLY ||
	command->ci == VAR) {
	if (nf) {
	    command->param = realloc(command->param, linelen - n + 1);
	    for (i=0; i<linelen-n; i++) 
		remainder[i] = line[i+n+1];
	    remainder[linelen-n] = '\0'; 
	    sscanf(remainder, "%s", command->param);
	    shiftleft(remainder, strlen(command->param));
	    strcpy(line, remainder);
	    nf -= 1;
	    n = 0;
	    linelen = strlen(line);
	} 
    }
    if (command->ci == MULTIPLY) {  /* suffix string */
	tmpstr = malloc(32);
	if (tmpstr == NULL) {
	    command->errcode = E_ALLOC;
	    return;
	}
	sscanf(line, "%s", tmpstr);
	strncpy(command->str, tmpstr, 3);
	shiftleft(line, strlen(tmpstr));
	nf -= 1;
	n = 0;
	linelen = strlen(line);
	free(tmpstr);
    }
    if (command->ci == AR) ar = 1;

    command->list = realloc(command->list, (1+nf) * sizeof(int));
    if (command->list == NULL) {
	command->errcode = E_ALLOC;
	strcpy (command->errmsg, "Memory allocation failed for command list");
	free(remainder);
	return;
    }
    command->list[0] = nf;

    /* now assemble the command list */
    for (j=1; j<=nf; j++) {

	for (i=0; i<linelen-n; i++) 
	    remainder[i] = line[i+n+1];
	remainder[linelen-n] = '\0'; 

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
		    if (laggenr(lagvar.varnum, lagvar.lag, 0, pZ, pdinfo)) {
			command->errcode = 1;
			sprintf(command->errmsg, 
				"generation of lag variable failed");
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
		    sprintf(command->errmsg, 
			    "'%s' is not the name of a variable", field);
		    free(remainder);
		    return;
		}
	    }
	}

	if (isdigit(field[0])) {
	    v = atoi(field);
	    if (!ar && v > pdinfo->v - 1) {
		command->errcode = 1;
		sprintf(command->errmsg, 
                       "%d is not a valid variable number", v);
		free(remainder);
		return;
	    }	
	    command->list[j] = v;
	}
	else if (field[0] == ';') {
	    if (command->ci == TSLS || command->ci == AR ||
		command->ci == SCATTERS) {
		command->param = realloc(command->param, 4);
		sprintf(command->param, "%d", j);
		n += strlen(field) + 1;
		command->list[j] = 999;
		ar = 0;
		continue;
	    }
	    else {
		command->list[0] -= 1;
		break;
	    }
	}

	if (!isalpha((unsigned char) field[0]) && 
	    !isdigit((unsigned char) field[0])) { 
	    command->errcode = 1;
	    sprintf(command->errmsg, 
		    "field '%s' in command is invalid", field);
	    free(remainder);
	    return;
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
    } else {
	/* command that needs a list but doesn't have one */
	if (command->list[0] == 0) command->errcode = E_ARGS;
    }

    if ((command->ci == AR || command->ci == TSLS || 
	 command->ci == SCATTERS) && strlen(command->param) == 0) {
	command->errcode = E_ARGS;
    }

    free(remainder);
    return;
}

/* ........................................................... */

int help (const char *cmd, const char *helpfile, print_t *prn)
{
    FILE *fq;
    char line[MAXLEN], tmp[MAXLEN];
    int ls, i, ok;

    if (cmd == NULL) {
	pprintf(prn, "\nValid gretl commands are:\n");
	for (i=1; i<NC; i++) {
	    pprintf(prn, "%-9s", commands[i]);
	    if (i%8 == 0) pprintf(prn, "\n");
	    else pprintf(prn, " ");
	}
	pprintf(prn, "\n");
	pprintf(prn, "\nFor help on a specific command, type: help cmdname");
	pprintf(prn, " (e.g. help smpl)\n");
	return 0;
    }

    ok = 0;
    for (i=1; i<NC; i++) {
	if (strcmp(commands[i], cmd) == 0) {
	    ok = 1;
	    break;
	}
    }
    if (!ok) {
	pprintf(prn, "\"%s\" is not a gretl command.\n", cmd);
	return 1;
    }

    if ((fq = fopen(helpfile, "r")) == NULL) {
	printf("Unable to access the file %s.\n", helpfile);
	return 1;
    } 

    while (fgets(line, MAXLEN, fq) != NULL) {
	delchar('\n', line);
	ls = strcmp(cmd, line);
	if (ls) continue;
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
	    pprintf(prn, "%s\n", tmp);
	} while (i);
	if (!ls) {
	    fclose(fq);
	    return 0;
	}
    }
    pprintf(prn, "%s: sorry, no help available.\n", cmd);
    fclose(fq);
    return 0;
}

/* ........................................................... */

static int parse_criteria (const char *line, const DATAINFO *pdinfo, 
			   double **pZ, print_t *prn)
{
    int i, n = pdinfo->n;
    double ess;
    int T, k;
    char cmd[9], essstr[32], Tstr[9], kstr[9];
    
    if (sscanf(line, "%s %s %s %s", cmd, essstr, Tstr, kstr) != 4) {
	return 1;
    }
    if (isalpha((unsigned char) essstr[0]) && 
	(i = varindex(pdinfo, essstr)) < pdinfo->v) 
	    ess = (*pZ)[n*i + pdinfo->t1];
    else if (isdigit(essstr[0])) ess = atof(essstr);
    else return 1;
    if (ess < 0) {
	pprintf(prn, "ess: negative value is out of bounds.\n");
	return 1;
    }
    if (isalpha((unsigned char) Tstr[0]) &&
	(i = varindex(pdinfo, Tstr)) < pdinfo->v) 
	    T = (int) (*pZ)[n*i + pdinfo->t1];
    else if (isdigit(Tstr[0])) T = atoi(Tstr);
    else return 1;
    if (T < 0) {
	pprintf(prn, "T: negative value is out of bounds.\n");
	return 1;
    }
    if (isalpha((unsigned char) kstr[0]) &&
	(i = varindex(pdinfo, kstr)) < pdinfo->v) 
	    k = (int) (*pZ)[n*i + pdinfo->t1];
    else if (isdigit(kstr[0])) k = atoi(kstr);
    else return 1;
    if (k < 0) {
	pprintf(prn, "k: negative value is out of bounds.\n");
	return 1;
    }    
    criteria(ess, T, k, prn);

    return 0;
}

/* ........................................................... */

int fcast (const char *line, const MODEL *pmod, DATAINFO *pdinfo, 
	   double **pZ, char *msg)
     /* return ID of var containing the forecast, or negative int on 
	error */
{
    int t, t1, t2, vi, n = pdinfo->n;
    char t1str[8], t2str[8], varname[9];

    *t1str = '\0'; *t2str = '\0';

    /* the varname should either be in the 2nd or 4th position */
    if (sscanf(line, "%*s %7s %7s %8s", t1str, t2str, varname) != 3) {
	if (sscanf(line, "%*s" "%8s", varname) != 1) return -1;
    }

    if (*t1str && *t2str) {
	t1 = dateton(t1str, pdinfo->pd, pdinfo->stobs, msg);
	t2 = dateton(t2str, pdinfo->pd, pdinfo->stobs, msg);
	if (t1 < 0 || t2 < 0 || t2 < t1) return -1;
    } else {
	t1 = pdinfo->t1;
	t2 = pdinfo->t2;
    }

    if (!isalpha((unsigned char) varname[0])) return -1;
    varname[8] = 0;
    vi = varindex(pdinfo, varname);

    if (vi >= pdinfo->v && grow_Z(1, pZ, pdinfo)) return -1 * E_ALLOC;

    strcpy(pdinfo->varname[vi], varname);
    strcpy(pdinfo->label[vi], "predicted values");

    for (t=0; t<pdinfo->n; t++) (*pZ)[n*vi + t] = NADBL;

    forecast(t1, t2, vi, pmod, pdinfo, pZ);

    return vi;
}

/* ........................................................... */
    
int add_new_var (DATAINFO *pdinfo, double **pZ, GENERATE *genr)
{
    int t, isconst = 1, n = pdinfo->n, v = genr->varnum;
    double xx;

    if (genr->special) return 0;
    /* is the new variable an addition to data set? */
    if (v >= pdinfo->v) {
	if (grow_Z(1, pZ, pdinfo)) return E_ALLOC;
	strcpy(pdinfo->varname[v], genr->varname);
	strcpy(pdinfo->label[v], genr->label);
    } else {
	strcpy(pdinfo->label[v], genr->label);
    }
    xx = genr->xvec[pdinfo->t1];
    for (t=pdinfo->t1+1; t<=pdinfo->t2; t++) {
	if (xx != genr->xvec[t]) {
	    isconst = 0;
	    break;
	}
    }
    if (isconst) {
	for (t=0; t<n; t++) (*pZ)[n*v + t] = xx;
    } else {
	for (t=0; t<n; t++) (*pZ)[n*v + t] = NADBL;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) 
	    (*pZ)[n*v + t] = genr->xvec[t];
    }
    if (genr->xvec != NULL) free(genr->xvec);
    return 0;
}

/* ........................................................... */

static int _full_list (const DATAINFO *pdinfo, CMD *command)
/* create a gretl "list" containing all the vars in the data set,
   except for the constant and any "hidden" variables */
{
    int i, n = 1;

    command->list = realloc(command->list, pdinfo->v * sizeof(int));
    if (command->list == NULL) return E_ALLOC;
    for (i=1; i<pdinfo->v; i++) {
	if (hidden_var(i, pdinfo)) continue;
	command->list[n++] = i;
    }
    command->list[0] = n - 1;
    return 0;
}

/* ........................................................... */

int parseopt (const char *s)
{
    if (strcmp(s, "-b") == 0 || strcmp(s, "--batch") == 0) 
	return OPT_BATCH;
    if (strcmp(s, "-h") == 0 || strcmp(s, "--help") == 0) 
	return OPT_HELP;
    if (strcmp(s, "-p") == 0 || strcmp(s, "--pvalue") == 0) 
	return OPT_PVALS;
    if (strcmp(s, "-v") == 0 || strcmp(s, "--version") == 0) 
	return OPT_VERSION;
    if (strcmp(s, "-r") == 0 || strcmp(s, "--run") == 0) 
	return OPT_RUNIT;
    if (strcmp(s, "-d") == 0 || strcmp(s, "--db") == 0) 
	return OPT_DBOPEN;
    if (strcmp(s, "-w") == 0 || strcmp(s, "--webdb") == 0) 
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
	perror("Try again later");
    }
    return 0;
}

#endif

/* ........................................................ */

void echo_cmd (CMD *pcmd, const DATAINFO *pdinfo, const char *line, 
	       const int batch, const int gui, const int oflag, 
	       print_t *prn)
     /* echo a given command: depending on the circumstances, either
	to stdout or to a buffer, or both */

{
    int i, err, got999 = 1;
    char flagc;

    if (strcmp(line, "quit") == 0 || line[0] == '!') return;

    if (pcmd->ci == AR) got999 = 0;

    if (!pcmd->nolist) { /* print list of params to command */
	if (!gui) {
	    if (!batch) printf(" %s", pcmd->cmd);
	    else printf("\n? %s", pcmd->cmd);
	    if (pcmd->ci == RHODIFF) printf(" %s;", pcmd->param);
	    else if (strlen(pcmd->param) && pcmd->ci != TSLS 
		     && pcmd->ci != AR && pcmd->ci != CORRGM
		     && pcmd->ci != SCATTERS) 
		printf(" %s", pcmd->param);
	}
	if (!batch) {
	    pprintf(prn, "%s", pcmd->cmd);
	    if (pcmd->ci == RHODIFF) pprintf(prn, " %s;", pcmd->param);
	    else if (strlen(pcmd->param) && pcmd->ci != TSLS 
		     && pcmd->ci != AR && pcmd->ci != CORRGM
		     && pcmd->ci != SCATTERS) 
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
		got999 = 1;
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
	err = list_dups(pcmd->list, pcmd->ci);
	if (err) {
	    printf("\nvar number %d duplicated in the command list.\n",
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

    printf("Listing labels for variables:\n");
    for (i=0; i<pdinfo->v; i++) {
	if (strlen(pdinfo->label[i]) > 2) {
	    printf("%3d) %-10s %s\n", i, 
		   pdinfo->varname[i], pdinfo->label[i]);
	}
    }
}

/* ........................................................ */

int simple_commands (CMD *cmd, const char *line, 
		     double **pZ, DATAINFO *datainfo, PATHS *paths,
		     const int batch, const int oflag, 
		     print_t *prn)
     /* common code for command-line and GUI client programs, where
	the command doesn't require special handling on the client
	side */
{
    int err = 0, order = 0;
    CORRMAT corrmat;

    switch (cmd->ci) {

    case ADF:
	if (!isdigit(cmd->param[0])) {
	    pprintf(prn, "adf: lag order must be given first\n");
	    break;
	}
	order = atoi(cmd->param);
	err = adf_test(order, cmd->list, pZ, datainfo, prn);
	if (err) errmsg(err, NULL, prn);
	break;

    case COINT:
	order = atoi(cmd->param);
	err = coint(order, cmd->list, pZ, datainfo, prn);
	break;

    case CORR:
	if (cmd->list[0] > 3) {
	    err = esl_corrmx(cmd->list, pZ, datainfo, 
			     batch, prn);
	    if (err) 
		pprintf(prn, "Error in generating correlation matrix\n");
	    break;
	}
	corrmat = corrlist(cmd->list, *pZ, datainfo);
	if ((err = corrmat.errcode)) 
	    pprintf(prn, "Couldn't allocate memory for correlation matrix.\n");
	else printcorr(cmd->list, corrmat, datainfo, prn);
	if (corrmat.r != NULL) free(corrmat.r);
	break;

    case CRITERIA:
	err = parse_criteria(line, datainfo, pZ, prn);
	if (err) 
	    pprintf(prn, "Error in computing model selection criteria.\n");
	break;

    case DIFF:
	err = list_diffgenr(cmd->list, pZ, datainfo);
	if (err) 
	    pprintf(prn, "Error adding first differences of variables.\n");
	else varlist(datainfo, prn);
	break;

    case LDIFF:
	err = list_ldiffgenr(cmd->list, pZ, datainfo);
	if (err) 
	    pprintf(prn, "Error adding log differences of variables.\n");
	else varlist(datainfo, prn);
	break;

    case LAGS:
	err = lags(cmd->list, pZ, datainfo); 
	if (err) 
	    pprintf(prn, "Error adding lags of variables.\n");
	else varlist(datainfo, prn);
	break;

    case LOGS:
	err = logs(cmd->list, pZ, datainfo, NULL);
	if (err < cmd->list[0]) 
	    pprintf(prn, "Error adding logs of variables.\n");
	if (err > 0) { 
	    varlist(datainfo, prn);
	    err = 0;
	} else
	    err = 1;
	break;

    case MULTIPLY:
	err = multiply(cmd->param, cmd->list, cmd->str, pZ, datainfo);
	if (err) errmsg(err, NULL, prn);
	else varlist(datainfo, prn);
	break;

    case GRAPH:
	graph(cmd->list, *pZ, datainfo, oflag, prn);
	break;

    case PLOT:
	plot(cmd->list, *pZ, datainfo, oflag, batch, prn);
	break;

    case INFO:
	err = get_info(paths->hdrfile, prn);
	if (err) 
	    pprintf(prn, "Error reading data header file.\n");
	break;

    case LABELS:
	showlabels(datainfo);
	break;

    case LIST:
	varlist(datainfo, prn);
	break;

    case PRINT:
	printdata(cmd->list, pZ, datainfo, batch, oflag, prn);
	break;

    case SUMMARY:
	err = summary(cmd->list, pZ, datainfo, batch, prn);
	if (err) 
	    pprintf(prn, "generation of summary stats failed\n");
	break; 

    case MEANTEST:
	err = means_test(cmd->list, *pZ, datainfo, oflag, prn);
	if (err) errmsg(err, NULL, prn);
	break;	

    case VARTEST:
	err = vars_test(cmd->list, *pZ, datainfo, prn);
	if (err) errmsg(err, NULL, prn);
	break;

    case RUNS:
	err = runs_test(cmd->list, *pZ, datainfo, prn);
	break;

    case SPEARMAN:
	err = spearman(cmd->list, *pZ, datainfo, prn, oflag);
	break;

    default:
	break;
    }

    return err;
}





