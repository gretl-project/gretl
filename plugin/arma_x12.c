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

#include "libgretl.h"
#include <gtk/gtk.h>

#ifdef WIN32
# include <windows.h>
#else
# if GLIB_CHECK_VERSION(2,0,0)
#  define GLIB2
#  include <signal.h>
# endif /* GLIB_CHECK_VERSION */
#endif

#ifdef GLIB2

static int tramo_x12a_spawn (const char *workdir, const char *fmt, ...)
{
    va_list ap;
    int i, nargs;
    int ok;
    int status = 0, ret = 0;
    GError *error = NULL;
    gchar **argv = NULL;
    gchar *sout = NULL, *serr = NULL;
    char *s;

    argv = malloc(2 * sizeof *argv);
    if (argv == NULL) return 1;
    argv[0] = g_strdup(fmt);
    argv[1] = NULL;
    i = nargs = 1;

    va_start(ap, fmt);
    while ((s = va_arg(ap, char *))) {
	i++;
	argv = realloc(argv, (i+1) * sizeof *argv);
	if (argv == NULL) {
	    status = 1;
	    break;
	}
	argv[i-1] = g_strdup(s);
	argv[i] = NULL;
    }
    va_end(ap);

    if (status == 1) return 1;

    nargs = i;

    signal(SIGCHLD, SIG_DFL);

    ok = g_spawn_sync (workdir,
		       argv,
		       NULL,
		       0,
		       NULL,
		       NULL,
		       &sout,
		       &serr,
		       &status,
		       &error);

    if (!ok) {
	fprintf(stderr, "spawn: '%s'\n", error->message);
	g_error_free(error);
	ret = 1;
    } else if (serr && *serr) {
	fprintf(stderr, "stderr: '%s'\n", serr);
	ret = 1;
    } else if (status != 0) {
	fprintf(stderr, "status=%d: stdout: '%s'\n", status, sout);
	ret = 1;
    }

    if (serr != NULL) g_free(serr);
    if (sout != NULL) g_free(sout);

    if (ret != 0) fputc(' ', stderr);
    for (i=0; i<nargs; i++) {
	if (ret != 0) fprintf(stderr, "%s ", argv[i]);
	free(argv[i]);
    }
    free(argv);
    if (ret != 0) fputc('\n', stderr);
    
    return ret;
}

#endif

static int write_spc_file (const char *fname, 
			   const double **Z, const DATAINFO *pdinfo, 
			   int depvar, int p, int q) 
{
    double x;
    FILE *fp;
    int i, t;
    int startyr, startper;
    char *s, tmp[8];

    fp = fopen(fname, "w");
    if (fp == NULL) return 1;    

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif 

    x = date(pdinfo->t1, pdinfo->pd, pdinfo->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    s = strchr(tmp, '.');
    if (s != NULL) startper = atoi(s + 1);
    else startper = 1;

    fprintf(fp, "series{\n period=%d\n title=\"%s\"\n", pdinfo->pd, 
	    pdinfo->varname[depvar]);
    fprintf(fp, " start=%d.%d\n", startyr, startper);
    fputs(" data=(\n", fp);

    i = 0;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (na(Z[depvar][t])) {
	    fputs("-99999 ", fp); /* FIXME? */
	} else {
	    fprintf(fp, "%g ", Z[depvar][t]);
	}
	if ((i + 1) % 7 == 0) fputc('\n', fp);
	i++;
    }
    fputs(" )\n}\n", fp);

    fprintf(fp, "arima{\nmodel = (%d 0 %d)(0 0 0)\n}\n", p, q);
    fputs("estimate{\nprint = (acm itr lkf lks mdl est rts rcm)\n}\n", fp);
    /* adding 'rsd' gives print of residuals */

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    return 0;
}

static int check_arma_list (const int *list)
{
    int err = 0;

    if (list[0] != 4) err = 1;

    /* for now we'll accept ARMA (4,4) at max */
    else if (list[1] < 0 || list[1] > 4) err = 1;
    else if (list[2] < 0 || list[2] > 4) err = 1;
    else if (list[1] + list[2] == 0) err = 1;

    if (err) {
	gretl_errmsg_set(_("Syntax error in arma command"));
    }
    
    return err;
}

int write_x12_arma (int *list, const double **Z, 
		    const DATAINFO *pdinfo, PATHS *paths, 
		    const char *prog, const char *workdir,
		    char *fname)
{
    int err = 0;
    char varname[VNAMELEN];
#ifndef GLIB2
    char cmd[MAXLEN];
#endif
    int depvar, p, q;

    err = check_arma_list(list);
    if (err) return 1;

    p = list[1];
    q = list[2];
    depvar = list[4];

    /* sanity check */
    if (!pdinfo->vector[depvar]) {
	char msg[48];

	sprintf(msg, "%s %s", pdinfo->varname[depvar], 
		_("is a scalar"));
	gretl_errmsg_set(msg);
	return 1;
    }

    sprintf(varname, pdinfo->varname[depvar]);

    /* write out an .spc file */
    sprintf(fname, "%s%c%s.spc", workdir, SLASH, varname);
    write_spc_file(fname, Z, pdinfo, depvar, p, q);

    /* run the program */
#if defined(WIN32)
    sprintf(cmd, "\"%s\" %s -r -p -q", prog, varname);
    err = winfork(cmd, workdir, SW_SHOWMINIMIZED, 
		  CREATE_NEW_CONSOLE | HIGH_PRIORITY_CLASS);
#elif defined(GLIB2)
    err = tramo_x12a_spawn(workdir, prog, varname, "-r", "-p", "-q", "-n", NULL);
#else
    sprintf(cmd, "cd \"%s\" && \"%s\" %s -r -p -q -n >/dev/null", 
	    workdir, prog, varname);
    err = gretl_spawn(cmd);
#endif

    if (!err) {
	sprintf(fname, "%s%c%s.out", workdir, SLASH, varname); 
    }

    return err;
}

