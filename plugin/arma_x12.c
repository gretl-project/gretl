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

struct x12_arma_info {
    int yno; /* dependent variable */
    int p;   /* non-seasonal AR order */
    int q;   /* non-seasonal MA order */
    int P;   /* seasonal AR order */
    int Q;   /* seasonal MA order */
    int t1;  /* starting obs */
    int t2;  /* ending obs */
};

static void write_arma_model_stats (MODEL *pmod, const int *list, 
				    const double *y, 
				    const DATAINFO *pdinfo,
				    struct x12_arma_info *aset)
{
    double mean_error;
    int t;

    pmod->ci = ARMA;
    pmod->ifc = 1;
    pmod->nobs = pmod->t2 - pmod->t1 + 1; 
    pmod->dfn = aset->p + aset->q + aset->P + aset->Q;
    pmod->dfd = pmod->nobs - pmod->dfn;
    pmod->ncoeff = pmod->dfn + 1;

    pmod->list = gretl_list_copy(list);

    pmod->ybar = gretl_mean(pmod->t1, pmod->t2, y);
    pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, y);

    mean_error = pmod->ess = 0.0;
    for (t=0; t<pdinfo->n; t++) {
	if (!na(pmod->uhat[t])) {
	    pmod->yhat[t] = y[t] - pmod->uhat[t];
	    pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	    mean_error += pmod->uhat[t];
	} else {
	    pmod->yhat[t] = NADBL;
	}
    }

    mean_error /= pmod->nobs;
    gretl_model_set_double(pmod, "mean_error", mean_error);

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    pmod->tss = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->tss += (y[t] - pmod->ybar) * (y[t] - pmod->ybar);
    }

    pmod->fstt = pmod->dfd * (pmod->tss - pmod->ess) / (pmod->dfn * pmod->ess);

    pmod->rsq = pmod->adjrsq = NADBL;

    if (pmod->tss > 0) {
	pmod->rsq = 1.0 - (pmod->ess / pmod->tss);
	if (pmod->dfd > 0) {
	    double den = pmod->tss * pmod->dfd;

	    pmod->adjrsq = 1.0 - (pmod->ess * (pmod->nobs - 1) / den);
	}
    }
}

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

    sprintf(fname, "%s.out", path);
    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	return 1;
    }

    while (fgets(line, sizeof line, fp)) {
	if (!strncmp(line, " MODEL EST", 10)) print = 1;
	if (print) pputs(prn, line);
	if (!strncmp(line, " Estimatio", 10)) break;
    }

    fclose(fp);
    
    return 0;
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

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif 

    while (fgets(line, sizeof line, fp)) {
	if (sscanf(line, "%11s %lf", statname, &x) == 2) {
	    if (!strcmp(statname, "nobs")) pmod->nobs = (int) x;
	    else if (!strcmp(statname, "lnlkhd")) pmod->lnL = x;
	    else if (!strcmp(statname, "aic")) pmod->criterion[C_AIC] = x;
	    else if (!strcmp(statname, "bic")) pmod->criterion[C_BIC] = x;
	}
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    return 0;
}

/* Parse the roots information from the X12ARIMA output file foo.rts */

static int get_roots (const char *fname, MODEL *pmod, 
		      struct x12_arma_info *aset)
{
    FILE *fp;
    char line[132];
    int i, nr, err = 0;
    cmplx *roots;

    nr = aset->p + aset->q + aset->P + aset->Q;

    roots = malloc(nr * sizeof *roots);
    if (roots == NULL) return E_ALLOC;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	free(roots);
	return 1;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif 

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

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

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

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif 

    for (i=0; i<nt; i++) {
	pmod->vcv[i] = NADBL;
    }

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

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    return err;
}
#endif

/* Below: parse the coefficient estimates and standard errors from
   the X12ARIMA output file foo.est
*/

static int get_estimates (const char *fname, double *coeff, double *sderr,
			  struct x12_arma_info *aset)
{
    FILE *fp;
    char line[128], word[16];
    double b, se;
    int ptot = aset->p + aset->P;
    int qtot = aset->q + aset->Q;
    int nc = ptot + qtot + 1;
    int i, j;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	return 1;
    }

    for (i=0; i<nc; i++) {
	coeff[i] = sderr[i] = NADBL;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif 

    i = 1;
    j = ptot + 1;

    while (fgets(line, sizeof line, fp) && i < nc) {
	if (sscanf(line, "%15s", word) == 1) {
	    if (!strcmp(word, "Constant")) {
		if (sscanf(line, "%*s %*s %lf %lf", &b, &se) == 2) {
		    coeff[0] = b;
		    sderr[0] = se;
		}
	    } else if (!strcmp(word, "AR")) {
		if (sscanf(line, "%*s %*s %*s %*s %lf %lf", &b, &se) == 2) {
		    coeff[i] = b;
		    sderr[i] = se;
		    i++;
		}
	    } else if (!strcmp(word, "MA")) {
		if (sscanf(line, "%*s %*s %*s %*s %lf %lf", &b, &se) == 2) {
		    coeff[j] = -b;  /* MA sign conventions */
		    sderr[j] = se;
		    j++;
		}
	    }
	}
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    for (i=0; i<nc; i++) {
	if (na(coeff[i]) || na(sderr[i])) {
	    err = 1;
	    break;
	}
    }

    if (!err) {
	double narfac = 1.0;
	double sarfac = 1.0;
	double arfac;

	for (i=1; i<=aset->p; i++) {
	    narfac -= coeff[i];
	}

	for (i=1; i<=aset->P; i++) {
	    sarfac -= coeff[i+aset->p];
	}

	arfac = narfac * sarfac;

	coeff[0] *= arfac;
	sderr[0] *= arfac;
    }

    return err;
}

/* Parse the residuals from the X12ARIMA output file foo.rsd */

static double *get_uhat (const char *fname, const DATAINFO *pdinfo)
{
    FILE *fp;
    char line[64], date[9];
    double x, *uhat;
    int t, start = 0, nobs = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	return NULL;
    }

    uhat = malloc(pdinfo->n * sizeof *uhat);
    if (uhat == NULL) return NULL;

    for (t=0; t<pdinfo->n; t++) uhat[t] = NADBL;

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif 

    while (fgets(line, sizeof line, fp)) {
	if (*line == '-') {
	    start = 1;
	    continue;
	}
	if (start && sscanf(line, "%s %lf", date, &x) == 2) {
	    t = x12_date_to_n(date, pdinfo);
	    if (t >= 0 && t < pdinfo->n) {
		uhat[t] = x;
		nobs++;
	    }
	}
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    if (nobs == 0) {
	free(uhat);
	uhat = NULL;
    }

    return uhat;
}

static void 
populate_arma_model (MODEL *pmod, const int *list, const char *path, 
		     const double *y, const DATAINFO *pdinfo, 
		     struct x12_arma_info *aset)
{
    double *uhat = NULL, *yhat = NULL;
    double *coeff = NULL, *sderr = NULL;
    char fname[MAXLEN];
    int nc, err = 0;

    sprintf(fname, "%s.rsd", path);
    uhat = get_uhat(fname, pdinfo);
    if (uhat == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    nc = aset->p + aset->q + aset->P + aset->Q + 1;

    yhat = malloc(pdinfo->n * sizeof *yhat);
    coeff = malloc(nc * sizeof *coeff);
    sderr = malloc(nc * sizeof *sderr);
    if (yhat == NULL || coeff == NULL || sderr == NULL) {
	free(yhat);
	free(coeff);
	free(uhat);
	pmod->errcode = E_ALLOC;
	return;
    }

    coeff[0] = sderr[0] = 0.0;
    sprintf(fname, "%s.est", path);
    err = get_estimates(fname, coeff, sderr, aset);

    if (!err) {
	sprintf(fname, "%s.lks", path);
	err = get_ll_stats(fname, pmod);
    }

    if (!err) {
	sprintf(fname, "%s.rts", path);
	err = get_roots(fname, pmod, aset);
    }

#if 0
    if (!err) {
	sprintf(fname, "%s.acm", path);
	err = get_x12a_vcv(fname, pmod, nc);
    }
#endif

    if (err) {
	fprintf(stderr, "problem getting model info\n");
	pmod->errcode = E_FOPEN;
    } else {
	pmod->uhat = uhat;
	pmod->yhat = yhat;
	pmod->coeff = coeff;
	pmod->sderr = sderr;
	write_arma_model_stats(pmod, list, y, pdinfo, aset);
	gretl_model_add_arma_varnames(pmod, pdinfo, aset->yno,
				      aset->p, aset->q, 
				      aset->P, aset->Q,
				      0);
	if (aset->P || aset->Q) {
	    gretl_model_set_int(pmod, "arma_P", aset->P);
	    gretl_model_set_int(pmod, "arma_Q", aset->Q);
	}
    }
}

static void output_series_to_spc (const double *x, int t1, int t2, 
				  FILE *fp)
{
    int i, t;

    fputs(" data = (\n", fp);

    i = 0;
    for (t=t1; t<=t2; t++) {
	fprintf(fp, "%g ", x[t]);
	if ((i + 1) % 7 == 0) fputc('\n', fp);
	i++;
    }
    fputs(" )\n", fp);
}

static int 
arma_missobs_check (const double **Z, const DATAINFO *pdinfo,
		    struct x12_arma_info *aset)
{
    int misst = 0;
    int list[2];

    list[0] = 1;
    list[1] = aset->yno;

    aset->t1 = pdinfo->t1;
    aset->t2 = pdinfo->t2;

    if (check_for_missing_obs(list, &aset->t1, &aset->t2, Z, &misst)) {
	gchar *msg;

	msg = g_strdup_printf(_("Missing value encountered for "
				"variable %d, obs %d"), aset->yno, misst);
	gretl_errmsg_set(msg);
	g_free(msg);
	return 1;
    }       

    return 0;
}

static int write_spc_file (const char *fname, 
			   const double **Z, const DATAINFO *pdinfo,
			   struct x12_arma_info *aset,
			   int verbose)
{
    double x;
    FILE *fp;
    int startyr, startper;
    char *s, tmp[8];

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't write to '%s'\n", fname);
	return 1;  
    }  

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif 

    x = date(aset->t1, pdinfo->pd, pdinfo->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    s = strchr(tmp, '.');
    if (s != NULL) {
	startper = atoi(s + 1);
    } else if (pdinfo->pd > 1) {
	startper = 1;
    } else {
	startper = 0;
    }

    fprintf(fp, "series {\n period = %d\n title = \"%s\"\n", pdinfo->pd, 
	    pdinfo->varname[aset->yno]);
    if (startper > 0) {
	fprintf(fp, " start = %d.%d\n", startyr, startper);
    } else {
	fprintf(fp, " start = %d\n", startyr);
    }

    output_series_to_spc(Z[aset->yno], aset->t1, aset->t2, fp);
    fputs("}\n", fp);

    fputs("Regression {\n Variables = (const)\n}\n", fp);

    if (aset->P > 0 || aset->Q > 0) {
	fprintf(fp, "arima {\n model = (%d 0 %d)(%d 0 %d)\n}\n", 
		aset->p, aset->q, aset->P, aset->Q);
    } else {
	fprintf(fp, "arima {\n model = (%d 0 %d)\n}\n", aset->p, aset->q); 
    }

    if (verbose) {
	fputs("estimate {\n print = (acm itr lkf lks mdl est rts rcm)\n", fp);
    } else {
	fputs("estimate {\n print = (acm lkf lks mdl est rts rcm)\n", fp);
    }

    fputs(" save = (rsd est lks acm rts rcm)\n}\n", fp);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    return 0;
}

#define has_seasonals(l) (l[0] > 5 && l[3] == LISTSEP && l[6] == LISTSEP)

static int check_arma_list (const int *list)
{
    int err = 0;

    if (list[0] < 4) err = 1;

    /* for now we'll accept ARMA (4,4) at max */
    else if (list[1] < 0 || list[1] > 4) err = 1;
    else if (list[2] < 0 || list[2] > 4) err = 1;
    else if (list[1] + list[2] == 0) err = 1;

    if (!err && has_seasonals(list)) {
	if (list[4] < 0 || list[4] > 4) err = 1;
	else if (list[5] < 0 || list[5] > 4) err = 1;
	else if (list[4] + list[5] == 0) err = 1;
    }

    if (err) {
	gretl_errmsg_set(_("Syntax error in arma command"));
    }
    
    return err;
}

MODEL arma_x12_model (const int *list, const double **Z, const DATAINFO *pdinfo, 
		      const char *prog, const char *workdir,
		      gretlopt opt, int gui, PRN *prn)
{
    int err = 0;
    int verbose = (prn != NULL);
    char varname[VNAMELEN], path[MAXLEN];
    PRN *aprn = NULL;
#ifndef GLIB2
    char cmd[MAXLEN];
#endif
    struct x12_arma_info aset;
    MODEL armod;

    if (opt & OPT_V) {
	aprn = prn;
    }

    gretl_model_init(&armod);  
    gretl_model_smpl_init(&armod, pdinfo);

    if (check_arma_list(list)) {
	armod.errcode = E_UNSPEC;
	return armod;
    }

    aset.p = list[1];
    aset.q = list[2];

    if (has_seasonals(list)) {
	aset.P = list[4];
	aset.Q = list[5];
	aset.yno = list[7];
    } else {
	aset.P = aset.Q = 0;
	aset.yno = list[4];
    }

    /* sanity check */
    if (!pdinfo->vector[aset.yno]) {
	char msg[48];

	sprintf(msg, "%s %s", pdinfo->varname[aset.yno], 
		_("is a scalar"));
	gretl_errmsg_set(msg);
	armod.errcode = E_DATA;
	return armod;
    }

    if (arma_missobs_check(Z, pdinfo, &aset)) {
	armod.errcode = E_MISSDATA;
	return armod;
    }	

    sprintf(varname, pdinfo->varname[aset.yno]);

    /* write out an .spc file */
    sprintf(path, "%s%c%s.spc", workdir, SLASH, varname);
    write_spc_file(path, Z, pdinfo, &aset, verbose);

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
	const double *y = Z[aset.yno];

	sprintf(path, "%s%c%s", workdir, SLASH, varname); 
	armod.t1 = aset.t1;
	armod.t2 = aset.t2;
	populate_arma_model(&armod, list, path, y, pdinfo, &aset);
	if (verbose && !armod.errcode) {
	    print_iterations(path, aprn);
	}
	if (!armod.errcode && gui) {
	    add_unique_output_file(&armod, path);
	    gretl_model_set_int(&armod, "arma_by_x12a", 1);
	}	
    } else {
	armod.errcode = E_UNSPEC;
	gretl_errmsg_set(_("Failed to execute x12arima"));
    }

    return armod;
}

