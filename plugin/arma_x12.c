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
 */

#include "../cephes/mconf.h"

#include "libgretl.h"
#include "version.h"
#include "bhhh_max.h"
#include "libset.h"
#include "arma_priv.h"

#include <glib.h>

#ifdef WIN32
# include <windows.h>
# include <io.h>
#endif

#ifndef WIN32

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

    ok = g_spawn_sync(workdir,
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

#endif /* !WIN32 */

#define X12A_CODE
#include "arma_common.c"
#undef X12A_CODE

static int add_unique_output_file (MODEL *pmod, const char *path)
{
    char outname[FILENAME_MAX];
    char unique[FILENAME_MAX];
    char line[256];
    FILE *f0, *f1;

    sprintf(outname, "%s.out", path);
    f0 = gretl_fopen(outname, "r");
    if (f0 == NULL) {
	return E_FOPEN;
    }

    sprintf(unique, "%s.XXXXXX", outname);
    f1 = gretl_mktemp(unique, "w");
    if (f1 == NULL) {
	fclose(f0);
	return E_FOPEN;
    }

    while (fgets(line, sizeof line, f0)) {
	fputs(line, f1);
    }

    fclose(f0);
    fclose(f1);
    gretl_remove(outname);

    gretl_model_set_data(pmod, "x12a_output", g_strdup(unique),
			 GRETL_TYPE_STRING, strlen(unique) + 1);

    return 0;
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

#define non_yearly_frequency(f) (f != 1 && f != 4 && f != 12)

static int x12_date_to_n (const char *s, const DATASET *dset)
{
    char date[12] = {0};

    if (non_yearly_frequency(dset->pd)) {
	int t, maj = 0, min = 0, n = strlen(s);
	char fmt[16];

	sprintf(fmt, "%%%dd%%2d", n - 2);
	if (sscanf(s, fmt, &maj, &min) == 2) {
	    t = (maj - 1) * dset->pd + min - 1;
	} else {
	    t = -1;
	}

	return t;
    }

    if (dset->pd > 1) {
	int len = strlen(s);

	if (len <= 4) {
	    strncat(date, s, len - 2);
	    strcat(date, ":");
	    strcat(date, s + len - 2);
	} else {
	    strncat(date, s, 4);
	    strcat(date, ":");
	    strncat(date, s + 4, 4);
	}
    } else {
	strncat(date, s, 4);
    }

    return dateton(date, dset);
}

/* Parse the statistics from the X12ARIMA output file foo.lks */

static int get_ll_stats (const char *fname, MODEL *pmod)
{
    FILE *fp;
    char line[80], statname[12];
    int nefobs = 0;
    double x;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	return E_FOPEN;
    }

    pmod->sigma = NADBL;

    gretl_push_c_numeric_locale();

    while (fgets(line, sizeof line, fp)) {
	if (sscanf(line, "%11s %lf", statname, &x) == 2) {
	    if (!strcmp(statname, "nefobs")) nefobs = (int) x;
	    else if (!strcmp(statname, "var")) pmod->sigma = sqrt(x);
	    else if (!strcmp(statname, "lnlkhd")) pmod->lnL = x;
	    else if (!strcmp(statname, "aic")) pmod->criterion[C_AIC] = x;
	    else if (!strcmp(statname, "bic")) pmod->criterion[C_BIC] = x;
	    else if (!strcmp(statname, "hnquin")) pmod->criterion[C_HQC] = x;
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    if (nefobs > 0) {
	pmod->nobs = nefobs;
	pmod->t1 = pmod->t2 - nefobs + 1;
    }

    return 0;
}

/* Parse the roots information from the X12ARIMA output file foo.rts */

static int get_roots (const char *fname, MODEL *pmod, 
		      arma_info *ainfo)
{
    FILE *fp;
    char line[132];
    int i, nr, err = 0;
    cmplx *roots;

    nr = ainfo->p + ainfo->q + ainfo->P + ainfo->Q;
    if (nr == 0) {
	return 0;
    }

    roots = malloc(nr * sizeof *roots);
    if (roots == NULL) {
	return E_ALLOC;
    }

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	free(roots);
	return E_FOPEN;
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
	fprintf(stderr, "Error reading '%s'\n", fname);
	free(roots);
	roots = NULL;
	err = E_DATA;
    }

    if (roots != NULL) {
	gretl_model_set_data(pmod, "roots", roots, GRETL_TYPE_CMPLX_ARRAY,
			     nr * sizeof *roots);
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
   the X-12-ARIMA output file foo.est
*/

static int 
get_estimates (const char *fname, MODEL *pmod, arma_info *ainfo)
{
    double *ar_coeff = pmod->coeff + ainfo->ifc;
    double *ma_coeff = ar_coeff + ainfo->np + ainfo->P;
    double *x_coeff = ma_coeff + ainfo->nq + ainfo->Q;
    double *ar_sderr = pmod->sderr + ainfo->ifc;
    double *ma_sderr = ar_sderr + ainfo->np + ainfo->P;
    double *x_sderr = ma_sderr + ainfo->nq + ainfo->Q;
    FILE *fp;
    char line[128], word[16];
    double b, se;
    int i, j, k;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read from '%s'\n", fname);
	return E_FOPEN;
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
	    fprintf(stderr, "Error reading '%s'\n", fname);
	    err = E_DATA;
	    break;
	}
    }

    return err;
}

/* Parse the residuals from the x12a output file foo.rsd */

static int 
get_uhat (const char *fname, MODEL *pmod, const DATASET *dset)
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
	    t = x12_date_to_n(date, dset);
	    if (t >= 0 && t < dset->n) {
		pmod->uhat[t] = x;
		nobs++;
	    } 
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    if (nobs == 0) {
	fprintf(stderr, "Error reading '%s'\n", fname);
	err = E_DATA;
    }

    return err;
}

static void 
populate_x12a_arma_model (MODEL *pmod, const char *path, 
			  const DATASET *dset, 
			  arma_info *ainfo)
{
    char fname[MAXLEN];
    int err;

    pmod->t1 = ainfo->t1;
    pmod->t2 = ainfo->t2;
    pmod->ncoeff = ainfo->nc;
    pmod->full_n = dset->n;

    err = gretl_model_allocate_storage(pmod);
    if (err) {
	pmod->errcode = err;
	return;
    }

    sprintf(fname, "%s.est", path);
    err = get_estimates(fname, pmod, ainfo);

    if (!err) {
	sprintf(fname, "%s.rsd", path);
	err = get_uhat(fname, pmod, dset);
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
	fprintf(stderr, "problem reading X-12-ARIMA model info\n");
	pmod->errcode = err;
    } else {
	write_arma_model_stats(pmod, ainfo, dset);
    }
}

static void 
output_series_to_spc (const int *list, const DATASET *dset,
		      int t1, int t2, FILE *fp)
{
    int i, t, done = 0;

    for (t=t1; t<=t2 && !done; t++) {
	for (i=1; i<=list[0] && !done; i++) {
	    if (na(dset->Z[list[i]][t])) {
		fputs(" missingcode=-99999\n", fp);
		done = 1;
	    }
	}
    }

    fputs(" data=(\n", fp);

    for (t=t1; t<=t2; t++) {
	for (i=1; i<=list[0]; i++) {
	    if (na(dset->Z[list[i]][t])) {
		fputs("-99999 ", fp);
	    } else {
		fprintf(fp, "%.15g ", dset->Z[list[i]][t]);
	    }
	}
	fputc('\n', fp);
    }

    fputs(" )\n", fp);
}

static int *arma_info_get_x_list (arma_info *ainfo)
{
    int *xlist = NULL;
    int start = arma_list_y_position(ainfo);
    int i;

    xlist = gretl_list_new(ainfo->nexo);

    if (xlist != NULL) {
	for (i=1; i<=xlist[0]; i++) {
	    xlist[i] = ainfo->alist[i + start];
	}
    }

    return xlist;
}

static void 
make_x12a_date_string (int t, const DATASET *dset, char *str)
{
    double dx;
    int yr, subper = 0;
    char *s;

    if (non_yearly_frequency(dset->pd)) {
	int maj = t / dset->pd + 1;
	int min = t % dset->pd + 1;

	sprintf(str, "%d.%d", maj, min);
	return;
    } 

    dx = date_as_double(t, dset->pd, dset->sd0);
    yr = (int) dx;

    sprintf(str, "%g", dx);
    s = strchr(str, '.');

    if (s != NULL) {
	subper = atoi(s + 1);
    } else if (dset->pd > 1) {
	subper = 1;
    } 

    if (subper > 0) {
	sprintf(str, "%d.%d", yr, subper);
    } else {
	sprintf(str, "%d", yr);
    }    
}

static void x12_pdq_string (arma_info *ainfo, FILE *fp)
{
    fputc('(', fp);
    
    if (ainfo->pmask == NULL) {
	fprintf(fp, "%d", ainfo->p);
    } else {
	int i;

	fputc('[', fp);
	for (i=0; i<ainfo->p; i++) {
	    if (AR_included(ainfo, i)) {
		fprintf(fp, "%d", i+1);
		if (i < ainfo->p - 1) {
		    fputc(' ', fp);
		}
	    }
	}
	fputc(']', fp);
    }

    fprintf(fp, " %d ", ainfo->d);

    if (ainfo->qmask == NULL) {
	fprintf(fp, "%d", ainfo->q);
    } else {
	int i;

	fputc('[', fp);
	for (i=0; i<ainfo->q; i++) {
	    if (MA_included(ainfo, i)) {
		fprintf(fp, "%d", i+1);
		if (i < ainfo->q - 1) {
		    fputc(' ', fp);
		}
	    }
	}
	fputc(']', fp);
    }

    fputc(')', fp);
}

static int write_arma_spc_file (const char *fname, 
				const DATASET *dset,
				arma_info *ainfo, int pdmax, 
				gretlopt opt)
{
    int maxobs = pdmax * 50;
    int maxfcast = pdmax * 5;
    int ylist[2];
    int *xlist = NULL;
    FILE *fp;
    char datestr[12];
    int nfcast = 0;
    int t1 = ainfo->t1;
    int tmax;
    int i;

    if (ainfo->nexo > 0) {
	xlist = arma_info_get_x_list(ainfo);
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

    fprintf(fp, "series{\n title=\"%s\"\n period=%d\n", 
	    dset->varname[ainfo->yno], dset->pd);

    t1 -= ainfo->maxlag;
    
    make_x12a_date_string(t1, dset, datestr);
    fprintf(fp, " start=%s\n", datestr);

    ylist[0] = 1;
    ylist[1] = ainfo->yno;

    tmax = ainfo->t2;

    if ((opt & OPT_F) && ainfo->t2 < dset->n - 1) {
	int nobs;

	tmax = dset->n - 1;
	nobs = tmax - ainfo->t1 + 1; /* t1? */
	if (nobs > maxobs) {
	    tmax -= nobs - maxobs;
	}
	nfcast = tmax - ainfo->t2;
	if (nfcast > maxfcast) {
	    tmax -= nfcast - maxfcast;
	    nfcast -= nfcast - maxfcast;
	}
#if 0
	fprintf(stderr, "x12a: doing forecast: nfcast = %d\n", nfcast);
#endif
    } 

    output_series_to_spc(ylist, dset, t1, tmax, fp);

    if (tmax > ainfo->t2) {
	make_x12a_date_string(ainfo->t2, dset, datestr);
	fprintf(fp, " span = ( , %s)\n", datestr);
    }

    fputs("}\n", fp);

    /* regression specification */
    fputs("Regression {\n", fp);
    if (ainfo->ifc) {
	fputs(" variables = (const)\n", fp);
    }
    if (ainfo->nexo > 0) {
	fputs(" user = ( ", fp);
	for (i=1; i<=xlist[0]; i++) {
	    fprintf(fp, "%s ", dset->varname[xlist[i]]);
	}
	fputs(")\n", fp);
	output_series_to_spc(xlist, dset, t1, tmax, fp);
	free(xlist);
    }
    fputs("}\n", fp);

    /* arima specification */
    fputs("arima {\n model = ", fp);
    x12_pdq_string(ainfo, fp);
    if (ainfo->P > 0 || ainfo->Q > 0) {
	fprintf(fp, "(%d %d %d)\n}\n", ainfo->P, ainfo->D, ainfo->Q);
    } else {
	fputs("\n}\n", fp);
    }

    fputs("estimate {\n", fp);

    /* could enable here: "tol = XX, maxiter = NN"  
       the default is tol = 1.0e-5, maxiter = 200 
    */

    if (opt & OPT_V) {
	fputs(" print = (acm itr lkf lks mdl est rts rcm)\n", fp);
    } else {
	fputs(" print = (acm lkf lks mdl est rts rcm)\n", fp);
    }

    fputs(" save = (rsd est lks acm rts rcm)\n", fp);

    if (opt & OPT_C) {
	fputs(" exact = none\n", fp);
    }

    fputs("}\n", fp);

    if (nfcast > 0) {
	fputs("forecast {\n save = (ftr)\n", fp);
	fprintf(fp, " maxlead = %d\n}\n", nfcast);
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static void print_x12a_error_file (const char *fname, PRN *prn)
{
    FILE *fp = gretl_fopen(fname, "r");

    if (fp != NULL) {
	char line[512];

	while (fgets(line, sizeof line, fp)) {
	    pputs(prn, line);
	}
	fclose(fp);
    }
}

static void delete_old_files (const char *path)
{
    const char *exts[] = {
	"acm", "itr", "lkf", "lks", "mdl", "est", "rts",
	"rcm", "rsd", "ftr", "err", "log", "out", NULL
    };
    char old[MAXLEN];
    int i, n = strlen(path);

    for (i=0; exts[i] != NULL; i++) {
	*old = '\0';
	strncat(old, path, n - 3);
	strcat(old, exts[i]);
	gretl_remove(old);
    }
}

static void x12a_maybe_allow_missvals (arma_info *ainfo)
{
    if (arma_exact_ml(ainfo)) {
	ainfo->pflags |= ARMA_NAOK;
    }
}

MODEL arma_x12_model (const int *list, const int *pqspec,
		      const DATASET *dset, int pdmax, 
		      gretlopt opt, PRN *prn)
{
    int verbose = (opt & OPT_V);
    const char *prog = gretl_x12_arima();
    const char *workdir = gretl_x12_arima_dir();
    char yname[VNAMELEN], path[MAXLEN];
    MODEL armod;
    arma_info ainfo_s, *ainfo;
    int missv = 0, misst = 0;
#ifdef WIN32
    char *cmd;
#endif
    int err = 0;

    if (dset->t2 < dset->n - 1) {
	/* FIXME this is temporary (OPT_F -> generate forecast) */
	opt |= OPT_F;
    }

    ainfo = &ainfo_s;
    arma_info_init(ainfo, opt | OPT_X, pqspec, dset);
    ainfo->prn = set_up_verbose_printer(opt, prn);
    gretl_model_init(&armod, dset); 

    ainfo->alist = gretl_list_copy(list);
    if (ainfo->alist == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = arma_check_list(ainfo, dset, opt);
    }

    if (err) {
	armod.errcode = err;
	goto bailout;
    }     

    /* calculate maximum lag */
    calc_max_lag(ainfo);

    x12a_maybe_allow_missvals(ainfo);

    /* adjust sample range if need be */
    armod.errcode = arma_adjust_sample(ainfo, dset, &missv, &misst);
    if (armod.errcode) {
	goto bailout;
    } else if (missv > 0) {
	set_arma_missvals(ainfo);
    }

    ainfo->y = (double *) dset->Z[ainfo->yno]; /* it's really const */
    strcpy(yname, dset->varname[ainfo->yno]);

    /* write out an .spc file */
    sprintf(path, "%s%c%s.spc", workdir, SLASH, yname);
    write_arma_spc_file(path, dset, ainfo, pdmax, opt);

    /* remove any files from an old run, in case of error */
    delete_old_files(path);

    /* run the program */
#ifdef WIN32
    cmd = g_strdup_printf("\"%s\" %s -r -p -q", prog, yname);
    err = win_run_sync(cmd, workdir);
    g_free(cmd);
#else
    err = tramo_x12a_spawn(workdir, prog, yname, "-r", "-p", "-q", "-n", NULL);
#endif

    if (!err) {
	sprintf(path, "%s%c%s", workdir, SLASH, yname); 
	populate_x12a_arma_model(&armod, path, dset, ainfo);
	if (verbose && !armod.errcode) {
	    print_iterations(path, ainfo->prn);
	}
	if (!armod.errcode) {
	    if (gretl_in_gui_mode()) {
		add_unique_output_file(&armod, path);
	    }
	    gretl_model_set_int(&armod, "arma_flags", (int) ainfo->flags);
	}	
    } else {
	armod.errcode = E_UNSPEC;
	gretl_errmsg_set(_("Failed to execute x12arima"));
    }

    if (armod.errcode && ainfo->prn != NULL) {
	sprintf(path, "%s%c%s.err", workdir, SLASH, yname);
	print_x12a_error_file(path, ainfo->prn);
    }

    if (ainfo->prn != NULL && ainfo->prn != prn) {
	iter_print_callback(0, ainfo->prn);
	close_down_verbose_printer(ainfo->prn);
    }

 bailout:

    arma_info_cleanup(ainfo);

    return armod;
}
