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
#include "bhhh_max.h"
#include "libset.h"

#include "../cephes/mconf.h"

#include <glib.h>

#ifdef WIN32
# include <windows.h>
# include <io.h>
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
		       G_SPAWN_SEARCH_PATH,
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

#include "arma_common.c"

static int add_unique_output_file (MODEL *pmod, const char *path)
{
    char fname[FILENAME_MAX];
    char unique[FILENAME_MAX];
    int err;

    sprintf(fname, "%s.out", path);
    sprintf(unique, "%s.XXXXXX", fname);
    if (mktemp(unique) == NULL) return 1;

    err = rename(fname, unique);
    if (!err) {
	gretl_model_set_data(pmod, "x12a_output", g_strdup(unique),
			     strlen(fname) + 1);
    } 

    return err;
}

static int print_iterations (const char *path, PRN *prn)
{
    char fname[MAXLEN];
    FILE *fp;
    char line[129];
    int print = 0;
    int err = 0;

    sprintf(fname, "%s.out", path);
    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	err = 1;
    } else {
	while (fgets(line, sizeof line, fp)) {
	    if (!strncmp(line, " MODEL EST", 10)) print = 1;
	    if (print) pputs(prn, line);
	    if (!strncmp(line, " Estimatio", 10)) break;
	}
	fclose(fp);
    }
    
    return err;
}

static int x12_date_to_n (const char *s, const DATAINFO *pdinfo)
{
    char date[12];

    *date = 0;
    strncat(date, s, 4);
    if (pdinfo->pd > 1) {
	strcat(date, ":");
	strncat(date, s + 4, 4);
    }

    return dateton(date, pdinfo);
}

/* Parse the statistics from the X12ARIMA output file foo.lks */

static int get_ll_stats (const char *fname, MODEL *pmod)
{
    FILE *fp;
    char line[80], statname[12];
    double x;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	return 1;
    }

    gretl_push_c_numeric_locale();

    while (fgets(line, sizeof line, fp)) {
	if (sscanf(line, "%11s %lf", statname, &x) == 2) {
	    if (!strcmp(statname, "nobs")) pmod->nobs = (int) x;
	    else if (!strcmp(statname, "lnlkhd")) pmod->lnL = x;
	    else if (!strcmp(statname, "aic")) pmod->criterion[C_AIC] = x;
	    else if (!strcmp(statname, "bic")) pmod->criterion[C_BIC] = x;
	    else if (!strcmp(statname, "hnquin")) pmod->criterion[C_HQC] = x;
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

/* Parse the roots information from the X12ARIMA output file foo.rts */

static int get_roots (const char *fname, MODEL *pmod, 
		      struct arma_info *ainfo)
{
    FILE *fp;
    char line[132];
    int i, nr, err = 0;
    cmplx *roots;

    nr = ainfo->p + ainfo->q + ainfo->P + ainfo->Q;

    roots = malloc(nr * sizeof *roots);
    if (roots == NULL) {
	return E_ALLOC;
    }

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	free(roots);
	return 1;
    }

    gretl_push_c_numeric_locale();

    i = 0;
    while (fgets(line, sizeof line, fp) && i < nr) {
	double re, im;

	if (!strncmp(line, "AR", 2) || !strncmp(line, "MA", 2)) {
	    if (sscanf(line, "%*s %*s %*s %lf %lf", &re, &im) == 2) {
		roots[i].r = re;
		roots[i].i = im;
		i++;
	    }
	}
    }
    
    gretl_pop_c_numeric_locale();

    fclose(fp);

    if (i != nr) {
	free(roots);
	roots = NULL;
	err = 1;
    }

    if (roots != NULL) {
	gretl_model_set_data(pmod, "roots", roots, nr * sizeof *roots);
    }

    return err;
}

/* Note: X12ARIMA does not give the full covariance matrix: it gives
   it only for the ARMA terms, and not for the constant.  Also the
   signs of off-diagonal elements are hard to disentangle.
*/

#if 0
static int get_x12a_vcv (const char *fname, MODEL *pmod, int nc)
{
    FILE *fp;
    char line[1024], valstr[24];
    double x;
    int i, j, k, nt = (nc * nc + nc) / 2;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) return 1;

    pmod->vcv = malloc(nt * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) {
	fclose(fp);
	return 1;
    }

    for (i=0; i<nt; i++) {
	pmod->vcv[i] = NADBL;
    }

    gretl_push_c_numeric_locale();

    j = 1;
    while (fgets(line, sizeof line, fp)) {
	if (!strncmp(line, "Nonseas", 7)) {
	    char *p = line + strcspn(line, "+-");

	    for (i=1; i<nc; i++) {
		sscanf(p, "%22s", valstr);
		p += 22;
		if (i >= j) {
		    x = atof(valstr);
		    k = ijton(i, j, nc);
		    pmod->vcv[k] = x;
		}
	    }
	    j++;
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return err;
}
#endif

/* Below: parse the coefficient estimates and standard errors from
   the X12ARIMA output file foo.est
*/

static int 
get_estimates (const char *fname, MODEL *pmod, struct arma_info *ainfo)
{
    double *ar_coeff = pmod->coeff + ainfo->ifc;
    double *ma_coeff = ar_coeff + ainfo->p + ainfo->P;
    double *x_coeff = ma_coeff + ainfo->q + ainfo->Q;

    double *ar_sderr = pmod->sderr + ainfo->ifc;
    double *ma_sderr = ar_sderr + ainfo->p + ainfo->P;
    double *x_sderr = ma_sderr + ainfo->q + ainfo->Q;

    FILE *fp;
    char line[128], word[16];
    double b, se;
    int i, j, k;

    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	return 1;
    }

    for (i=0; i<ainfo->nc; i++) {
	pmod->coeff[i] = pmod->sderr[i] = NADBL;
    }

    gretl_push_c_numeric_locale();

    i = j = k = 0;

    while (fgets(line, sizeof line, fp) && i < ainfo->nc) {
	if (sscanf(line, "%15s", word) == 1) {
	    if (!strcmp(word, "Constant")) {
		if (sscanf(line, "%*s %*s %lf %lf", &b, &se) == 2) {
		    pmod->coeff[0] = b;
		    pmod->sderr[0] = se;
		}
	    } else if (!strcmp(word, "User-defined")) {
		if (sscanf(line, "%*s %*s %lf %lf", &b, &se) == 2) {
		    x_coeff[i] = b;
		    x_sderr[i] = se;
		    i++;
		}		
	    } else if (!strcmp(word, "AR")) {
		if (sscanf(line, "%*s %*s %*s %*s %lf %lf", &b, &se) == 2) {
		    ar_coeff[j] = b;
		    ar_sderr[j] = se;
		    j++;
		}
	    } else if (!strcmp(word, "MA")) {
		if (sscanf(line, "%*s %*s %*s %*s %lf %lf", &b, &se) == 2) {
		    ma_coeff[k] = -b;  /* MA sign conventions */
		    ma_sderr[k] = se;
		    k++;
		}
	    }
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    for (i=0; i<ainfo->nc; i++) {
	if (na(pmod->coeff[i]) || na(pmod->sderr[i])) {
	    err = 1;
	    break;
	}
    }

    /* FIXME seasonal case!! */

    if (!err) {
	double narfac = 1.0;
	double sarfac = 1.0;
	double arfac;

	for (i=1; i<=ainfo->p; i++) {
	    narfac -= pmod->coeff[i];
	}

	for (i=1; i<=ainfo->P; i++) {
	    sarfac -= pmod->coeff[i+ainfo->p];
	}

	arfac = narfac * sarfac;

	pmod->coeff[0] *= arfac;
	pmod->sderr[0] *= arfac;

	for (i=0; i<ainfo->r; i++) {
	    x_coeff[i] *= arfac;
	    x_sderr[i] *= arfac;
	}
    }

    return err;
}

/* Parse the residuals from the X12ARIMA output file foo.rsd */

static int 
get_uhat (const char *fname, MODEL *pmod, const DATAINFO *pdinfo)
{
    FILE *fp;
    char line[64], date[9];
    double x;
    int t, start = 0, nobs = 0;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	return E_FOPEN;
    }

    gretl_push_c_numeric_locale();

    while (fgets(line, sizeof line, fp)) {
	if (*line == '-') {
	    start = 1;
	    continue;
	}
	if (start && sscanf(line, "%s %lf", date, &x) == 2) {
	    t = x12_date_to_n(date, pdinfo);
	    if (t >= 0 && t < pdinfo->n) {
		pmod->uhat[t] = x;
		nobs++;
	    }
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    if (nobs == 0) {
	err = E_DATA;
    }

    return err;
}

static void 
populate_arma_model (MODEL *pmod, const int *list, const char *path, 
		     const double **Z, const DATAINFO *pdinfo, 
		     struct arma_info *ainfo)
{
    char fname[MAXLEN];
    int t, err = 0;

    pmod->uhat = malloc(pdinfo->n * sizeof *pmod->uhat);
    pmod->yhat = malloc(pdinfo->n * sizeof *pmod->yhat);
    pmod->coeff = malloc(ainfo->nc * sizeof *pmod->coeff);
    pmod->sderr = malloc(ainfo->nc * sizeof *pmod->sderr);

    if (pmod->uhat == NULL || pmod->yhat == NULL || 
	pmod->coeff == NULL || pmod->sderr == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    pmod->full_n = pdinfo->n;

    for (t=0; t<pdinfo->n; t++) {
	pmod->uhat[t] = pmod->yhat[t] = NADBL;
    }

    sprintf(fname, "%s.est", path);
    err = get_estimates(fname, pmod, ainfo);

    if (!err) {
	sprintf(fname, "%s.rsd", path);
	err = get_uhat(fname, pmod, pdinfo);
    }

    if (!err) {
	sprintf(fname, "%s.lks", path);
	err = get_ll_stats(fname, pmod);
    }

    if (!err) {
	sprintf(fname, "%s.rts", path);
	err = get_roots(fname, pmod, ainfo);
    }

#if 0
    if (!err) {
	sprintf(fname, "%s.acm", path);
	err = get_x12a_vcv(fname, pmod, nc);
	/* also .rcm */
    }
#endif

    if (err) {
	fprintf(stderr, "problem getting model info\n");
	pmod->errcode = E_FOPEN;
    } else {
	write_arma_model_stats(pmod, NULL, list, Z, NULL, ainfo);
	gretl_model_add_arma_varnames(pmod, pdinfo, ainfo->yno,
				      ainfo->p, ainfo->q, 
				      ainfo->P, ainfo->Q,
				      ainfo->r);
    }
}

static void 
output_series_to_spc (const int *list, const double **Z, 
		      int t1, int t2, FILE *fp)
{
    int i, t;

    fputs(" data = (\n", fp);

    for (t=t1; t<=t2; t++) {
	for (i=1; i<=list[0]; i++) {
	    fprintf(fp, "%g ", Z[list[i]][t]);
	}
	fputc('\n', fp);
    }

    fputs(" )\n", fp);
}

static int *
arma_info_get_x_list (struct arma_info *ainfo, const int *alist)
{
    int *xlist = NULL;
    int start = (ainfo->seasonal)? 7 : 4;
    int i;

    xlist = gretl_list_new(ainfo->r);

    if (xlist != NULL) {
	for (i=1; i<=xlist[0]; i++) {
	    xlist[i] = alist[i + start];
	}
    }

    return xlist;
}

static int write_spc_file (const char *fname, 
			   const double **Z, const DATAINFO *pdinfo,
			   struct arma_info *ainfo, const int *alist,
			   int verbose)
{
    int ylist[2];
    int *xlist = NULL;
    double dx;
    FILE *fp;
    int startyr, startper;
    char *s, tmp[8];
    int i;

    if (ainfo->r > 0) {
	xlist = arma_info_get_x_list(ainfo, alist);
	if (xlist == NULL) {
	    return E_ALLOC;
	}
    }

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't write to '%s'\n", fname);
	return 1;  
    } 

    gretl_push_c_numeric_locale();

    dx = date(ainfo->t1, pdinfo->pd, pdinfo->sd0);
    startyr = (int) dx;
    sprintf(tmp, "%g", dx);
    s = strchr(tmp, '.');
    if (s != NULL) {
	startper = atoi(s + 1);
    } else if (pdinfo->pd > 1) {
	startper = 1;
    } else {
	startper = 0;
    }

    fprintf(fp, "series {\n period = %d\n title = \"%s\"\n", pdinfo->pd, 
	    pdinfo->varname[ainfo->yno]);
    if (startper > 0) {
	fprintf(fp, " start = %d.%d\n", startyr, startper);
    } else {
	fprintf(fp, " start = %d\n", startyr);
    }

    ylist[0] = 1;
    ylist[1] = ainfo->yno;
    output_series_to_spc(ylist, Z, ainfo->t1, ainfo->t2, fp);
    fputs("}\n", fp);

    /* regression specification */
    fputs("Regression {\n", fp);
    if (ainfo->ifc) {
	fputs(" variables = (const)\n", fp);
    }
    if (ainfo->r > 0) {
	fputs(" user = ( ", fp);
	for (i=1; i<=xlist[0]; i++) {
	    fprintf(fp, "%s ", pdinfo->varname[xlist[i]]);
	}
	fputs(")\n", fp);
	output_series_to_spc(xlist, Z, ainfo->t1, ainfo->t2, fp);
	free(xlist);
    }
    fputs("\n}\n", fp);

    /* arima specification */
    if (ainfo->P > 0 || ainfo->Q > 0) {
	fprintf(fp, "arima {\n model = (%d 0 %d)(%d 0 %d)\n}\n", 
		ainfo->p, ainfo->q, ainfo->P, ainfo->Q);
    } else {
	fprintf(fp, "arima {\n model = (%d 0 %d)\n}\n", ainfo->p, ainfo->q); 
    }

    if (verbose) {
	fputs("estimate {\n print = (acm itr lkf lks mdl est rts rcm)\n", fp);
    } else {
	fputs("estimate {\n print = (acm lkf lks mdl est rts rcm)\n", fp);
    }

    fputs(" save = (rsd est lks acm rts rcm)\n}\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

MODEL arma_x12_model (const int *list, const double **Z, const DATAINFO *pdinfo,
		      gretlopt opt, PRN *prn)
{
    int verbose = (opt & OPT_V);
    const char *prog = gretl_x12_arima();
    const char *workdir = gretl_x12_arima_dir();
    char yname[VNAMELEN], path[MAXLEN];
    int *alist = NULL;
    PRN *aprn = NULL;
#ifndef GLIB2
    char cmd[MAXLEN];
#endif
    MODEL armod;
    struct arma_info ainfo;
    int err = 0;

    if (verbose) {
	aprn = prn;
    } 

    ainfo.atype = ARMA_X12A;
    ainfo.dx = NULL;
    ainfo.T = pdinfo->n;
    ainfo.pd = pdinfo->pd;

    gretl_model_init(&armod); 
    gretl_model_smpl_init(&armod, pdinfo);

    alist = gretl_list_copy(list);
    if (alist == NULL) {
	armod.errcode = E_ALLOC;
	goto bailout;
    }

    if (opt & OPT_I) {
	err = check_arima_list(alist, opt, Z, pdinfo, &ainfo);
    } else {
	err = check_arma_list(alist, opt, Z, pdinfo, &ainfo);
    }

    if (err) {
	armod.errcode = err;
	goto bailout;
    }

    /* dependent variable */
    if (opt & OPT_I) {
	ainfo.yno = (ainfo.seasonal)? alist[9] : alist[5];
    } else {
	ainfo.yno = (ainfo.seasonal)? alist[7] : alist[4];
    }

    /* calculate maximum lag */
    calc_max_lag(&ainfo);

    /* adjust sample range if need be */
    if (arma_adjust_sample(pdinfo, Z, alist, &ainfo)) {
        armod.errcode = E_DATA;
	goto bailout;
    }	

    sprintf(yname, pdinfo->varname[ainfo.yno]);

    /* write out an .spc file */
    sprintf(path, "%s%c%s.spc", workdir, SLASH, yname);
    write_spc_file(path, Z, pdinfo, &ainfo, alist, verbose);

    /* run the program */
#if defined(WIN32)
    sprintf(cmd, "\"%s\" %s -r -p -q", prog, yname);
    err = winfork(cmd, workdir, SW_SHOWMINIMIZED, 
		  CREATE_NEW_CONSOLE | HIGH_PRIORITY_CLASS);
#elif defined(GLIB2)
    err = tramo_x12a_spawn(workdir, prog, yname, "-r", "-p", "-q", "-n", NULL);
#else
    sprintf(cmd, "cd \"%s\" && \"%s\" %s -r -p -q -n >/dev/null", 
	    workdir, prog, yname);
    err = gretl_spawn(cmd);
#endif

    if (!err) {
	sprintf(path, "%s%c%s", workdir, SLASH, yname); 
	armod.t1 = ainfo.t1;
	armod.t2 = ainfo.t2;
	armod.nobs = armod.t2 - armod.t1 + 1;
	populate_arma_model(&armod, alist, path, Z, pdinfo, &ainfo);
	if (verbose && !armod.errcode) {
	    print_iterations(path, aprn);
	}
	if (!armod.errcode && gretl_in_gui_mode()) {
	    add_unique_output_file(&armod, path);
	    gretl_model_set_int(&armod, "arma_by_x12a", 1);
	}	
    } else {
	armod.errcode = E_UNSPEC;
	gretl_errmsg_set(_("Failed to execute x12arima"));
    }

 bailout:

    free(alist);

    return armod;
}

