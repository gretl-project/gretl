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

/* describe.c - gretl descriptive statistics */

#include "libgretl.h"
#include "internal.h"
#include <unistd.h>

#ifdef OS_WIN32
# include <windows.h>
#endif

extern void _mxout (const double *rr, const int *list, int ci,
		    const DATAINFO *pdinfo, int pause, PRN *prn);

/* ........................................................... */

static int missvals (double *x, int n)
{
    int t;
    
    for (t=0; t<n; t++)
	if (na(x[t])) return 1;
	    
    return 0;
}

/**
 * esl_median:
 * @zx: data series (which should be pre-sorted).
 * @n: length of the series.
 *
 * Returns: the median value of the given series.
 *
 */

double esl_median (const double *zx, int n)
{
    double xx;
    int n2, n2p;

    n2p = (n2 = n/2) + 1;
    xx = (n % 2)? zx[n2p - 1] : 0.5 * (zx[n2 - 1] + zx[n2p - 1]);
    return xx;
}

/* ........................................................... */

static void moments (int t1, int t2, const double *zx, 
		     double *xbar, double *std, 
		     double *skew, double *kurt, int k)
     /* k is the "degrees of freedom loss": it will generally be one,
	other than when dealing with a regression residual */
{
    int t, n = t2 - t1 + 1;
    double sum = 0.0, p, var;

    if (_isconst(t1, t2, zx)) {
	*xbar = zx[t1];
	*std = 0.0;
	*skew = *kurt = NADBL;
	return;
    }
    for (t=t1; t<=t2; t++) sum += zx[t];
    *xbar = sum/n;
    var = *skew = *kurt = 0.0;
    for (t=t1; t<=t2; t++) {
	sum = zx[t] - *xbar;
	p = sum * sum;
	var += p;
	p *= sum;
	*skew = *skew + p;
	p *= sum;
	*kurt = *kurt + p;
    }
    var /= (n - k);
    if (var < 0.0) {
	*std = *skew = *kurt = NADBL;
	return;
    }
    *std = sqrt(var);
    if (var > 0.0) {
	*skew /= (n * var * (*std));
	*kurt = (*kurt)/(n * var * var) - 3.0;
    } else *skew = *kurt = NADBL;
}

/**
 * free_freq:
 * @freq: gretl frequency distribution struct
 *
 * Frees all malloced elements of the struct.
 *
 */

void free_freq (FREQDIST *freq)
{
    free(freq->midpt);
    free(freq->endpt);
    free(freq->f);
    free(freq);
}

/**
 * freqdist:
 * @pZ: pointer to data matrix
 * @pdinfo: information on the data set.
 * @varno: ID number of variable to process.
 * @params: degrees of freedom loss (generally = 1 unless we're dealing
 * with the residual from a regression)
 *
 * Calculates the frequency distribution for the specified variable.
 *
 * Returns: struct containing the distribution.
 *
 */

FREQDIST *freqdist (double ***pZ, const DATAINFO *pdinfo, 
		    int varno, int params)
{
    FREQDIST *freq;
    double *x = NULL;
    double xx, xmin, xmax, xbar, sdx;
    double skew, kurt;
    int i, k, n, t;
    int maxend = 16;

    freq = malloc(sizeof *freq);
    if (freq == NULL) return NULL;

    gretl_errno = 0;
    gretl_errmsg[0] = '\0';
    freq->midpt = NULL;
    freq->endpt = NULL;
    freq->f = NULL;

    x = malloc((pdinfo->t2 - pdinfo->t1 + 1) * sizeof *x);
    if (x == NULL) {
	sprintf(gretl_errmsg, _("Out of memory in frequency distribution"));
	free(freq);
	return NULL;
    }
    n = ztox(varno, x, *pZ, pdinfo);
    if (n < 3) {
	gretl_errno = E_DATA;
	sprintf(gretl_errmsg, _("Insufficient data to build frequency "
		"distribution for variable %s"), pdinfo->varname[varno]);
	free(x);
	return freq;
    }
    freq->t1 = pdinfo->t1; 
    freq->t2 = pdinfo->t2;

    strcpy(freq->varname, pdinfo->varname[varno]);

    if (_isconst(0, n-1, x)) {
	gretl_errno = 1;
	sprintf(gretl_errmsg, _("%s is a constant"), freq->varname);
	return freq;
    }    
    
    moments(0, n-1, x, &freq->xbar, &freq->sdx, &skew, &kurt, params);
    xbar = freq->xbar;
    sdx = freq->sdx;

    freq->endpt = malloc((maxend + 1) * sizeof(double));
    freq->midpt = malloc(maxend * sizeof(double));
    freq->f = malloc(maxend * sizeof(int));
    if (freq->endpt == NULL || freq->midpt == NULL ||
	freq->f == NULL) {
	gretl_errno = E_ALLOC;
	strcpy(gretl_errmsg, _("Out of memory for frequency distribution"));
	free(x);
	return freq;
    }
    
    _minmax(0, n-1, x, &xmin, &xmax);
    freq->n = n;
    freq->endpt[0] = xmin;

    xx = xbar - 3.75 * sdx;
    sdx /= 2.0; 
    while (xx < xmin) xx += sdx;
    freq->endpt[1] = xx;
    freq->endpt[maxend] = xmax;
    for (t=2; t<maxend; t++) {
	xx += sdx;
	if (xx > xmax) {
	    freq->endpt[t] = xmax;
	    break;
	}
	freq->endpt[t] = xx;
    }
    freq->numbins = t;

    for (k=0; k<freq->numbins; k++) {
	freq->f[k] = 0;
	freq->midpt[k] = (freq->endpt[k] + freq->endpt[k+1])/2;
    }
    for (i=0; i<n; i++) {
	xx = x[i];
	if (xx < freq->endpt[1]) {
	    freq->f[0] += 1;
	    continue;
	}
	if (xx >= freq->endpt[freq->numbins]) {
	    freq->f[freq->numbins-1] += 1;
	    continue;
	}
	for (k=1; k<freq->numbins; k++) 
	    if (freq->endpt[k] <= xx && xx <freq->endpt[k+1])
		freq->f[k] += 1;
    }

    freq->chisqu = freq->n * (skew * skew/6.0 + kurt * kurt/24.0); 

    free(x);

    return freq;
}

/* ...................................................... */

static int get_pacf (double *pacf, int *maxlag, int varnum, 
		     double ***pZ, DATAINFO *pdinfo)
{
    int i, j, err = 0, *laglist, *list;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int nobs = t2 - t1 + 1;
    int v = pdinfo->v;
    MODEL tmp;

    pdinfo->t1 = 0;

    *maxlag = 14;
    if (*maxlag > nobs / 2 - 1) *maxlag = nobs / 2 - 1;

    list = malloc((*maxlag + 3) * sizeof *list);
    laglist = malloc(*maxlag * sizeof *laglist);
    if (list == NULL || laglist == NULL) {
	pdinfo->t1 = t1;
	return 1;
    }

    /* add appropriate number of lags to data set */
    for (i=1; i<=*maxlag; i++) {
	_laggenr(varnum, i, 0, pZ, pdinfo);
	/* the lagvar may already exist */
	laglist[i-1] = _lagvarnum(varnum, i, pdinfo); 
    }

    _init_model(&tmp, pdinfo);
    pdinfo->t1 = t1;

    list[1] = varnum;
    for (i=2; i<=*maxlag; i++) {
	list[0] = i + 2;
	list[i+2] = 0;
	for (j=2; j<i+2; j++) list[j] = laglist[j-2];
	tmp = lsq(list, pZ, pdinfo, OLS, 0, 0);
	if ((err = tmp.errcode)) break;
	pacf[i-1] = tmp.coeff[i];
	if (i < *maxlag) clear_model(&tmp, pdinfo);
    }

    clear_model(&tmp, pdinfo);
    dataset_drop_vars(pdinfo->v - v, pZ, pdinfo);
    free(laglist);
    free(list);

    return err;
}

/**
 * corrgram:
 * @varno: ID number of variable to process.
 * @order: integer order for autocorrelation function.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @ppaths: struct containing path information.
 * @batch: if = 1, use ASCII graphic rather than gnuplot graph.
 * @prn: gretl printing struct.
 *
 * Computes autocorrelation function and plots the correlogram for
 * the variable specified in @list.
 *
 * Returns: 0 on successful completion, error code on error.
 *
 */

int corrgram (int varno, int order, double ***pZ, 
	      DATAINFO *pdinfo, PATHS *ppaths, 
	      int batch, PRN *prn)
{
    double *x, *y, *acf, *xl, box;
    double *pacf = NULL;
    int err = 0, k, l, m, nobs, n = pdinfo->n; 
    int maxlag = 0, t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int list[2];
    FILE *fq = NULL;

    list[0] = 1;
    list[1] = varno;
    _adjust_t1t2(NULL, list, &t1, &t2, *pZ, NULL);
    nobs = t2 - t1 + 1;

    if (missvals(&(*pZ)[varno][t1], nobs)) {
	pprintf(prn, "\n%s",
		_("Missing values within sample -- can't do correlogram"));
	return 1;
    }

    if (nobs < 4) {
	pprintf(prn, _("\nInsufficient observations for correlogram"));
	return 1;
    }
    if (_isconst(t1, t2, &(*pZ)[varno][0])) {
	sprintf(gretl_tmp_str, _("%s is a constant"), 
		pdinfo->varname[varno]);
	pprintf(prn, "\n%s\n", gretl_tmp_str);
	return 1;
    }

    switch (pdinfo->pd) {
    case 4: 
	m = (nobs <= 20)? nobs - 5 : 16;
	break;
    case 12: 
    case 52: 
	m = (nobs <= 40)? nobs - 13 : 36;
	break;
    case 24: 
	m = (nobs <= 100)? nobs - 25 : 96;
	break;
    default:  
	m = (nobs <= 18)? nobs - 5 : 14;
	break;
    }

    if (order && m > order) m = order;

    x = malloc(n * sizeof *x);
    y = malloc(n * sizeof *y);
    acf = malloc((m + 1) * sizeof *acf);
    if (x == NULL || y == NULL || acf == NULL)
	return E_ALLOC;    

    for (l=1; l<=m; l++) {
	for (t=t1+l; t<=t2; t++) {
	    k = t - (t1+l);
	    x[k] = (*pZ)[varno][t];
	    y[k] = (*pZ)[varno][t-l];
	}
	acf[l] = _corr(nobs-l, x, y);
    }

    sprintf(gretl_tmp_str, _("Autocorrelation function for %s"), 
	    pdinfo->varname[varno]);
    pprintf(prn, "\n%s\n\n", gretl_tmp_str);

    /* add Box-Pierce statistic */
    box = 0;
    for (t=1; t<=m; t++) 
	box += acf[t] * acf[t];
    box *= m;
    pprintf(prn, _("Box-Pierce Q statistic = %.4f\n"), box);
    pprintf(prn, _("Degrees of freedom = %d, significance level = %.4f\n\n"),
	    m, chisq(box, m));

    for (t=1; t<=m; t++) {
	pprintf(prn, "%5d)%7.3f", t, acf[t]);
	if (t%5 == 0) pprintf(prn, "\n");
    }
    pprintf(prn, "\n");

    if (batch) { /* use ASCII graphics, not gnuplot */
	xl = malloc(m * sizeof *xl);
	if (xl == NULL) return E_ALLOC;
	for (l=0; l<m; l++) xl[l] = l + 1.0;
        pprintf(prn, "\n\n%s\n\n", _("Correlogram"));
	_graphyzx(NULL, acf + 1, NULL, xl, m, pdinfo->varname[varno], 
		  _("lag"), NULL, 0, prn);
	free(x);
	free(xl);
	free(y);
	free(acf);
	return 0;	
    }

    /* generate partial autocorrelation function */
    pacf = malloc(14 * sizeof *pacf);
    if (pacf == NULL) {
	err = E_ALLOC;
	goto getout;
    }
    err = get_pacf(pacf, &maxlag, varno, pZ, pdinfo);
    pacf[0] = acf[1];
    if (!err) {
	pprintf(prn, "\n%s", _("Partial autocorrelations"));
	if (maxlag < m) 
	    pprintf(prn, " (%s %d):\n\n", _("to lag"), maxlag);
	else
	    pprintf(prn, ":\n\n");
	for (l=1; l<=maxlag; l++) {
	    pprintf(prn, "%5d)%7.3f", l, pacf[l-1]);
	    if (l%5 == 0) pprintf(prn, "\n");
	}
    }
    pprintf(prn, "\n");
    if (maxlag % 5 > 0) pprintf(prn, "\n");

    if (gnuplot_init(ppaths, &fq)) return E_FOPEN;

    fprintf(fq, "# correlogram\n");
    fprintf(fq, "set xlabel \"%s\"\n", _("lag"));
    fprintf(fq, "set xzeroaxis\n");
    fprintf(fq, "set title \"%s %s\"\n", I_("Correlogram for"), 
	    pdinfo->varname[varno]);
    if (maxlag) {
	fprintf(fq, "plot '-' using 1:2 title '%s' "
		"w impulses, \\\n'-' using 1:2 title '%s",
		_("autocorrelations"), _("partial autocorrelations"));
	if (maxlag < m)
	    fprintf(fq, "(%s %d)' w impulses\n", I_("to lag"), maxlag);
	else
	    fprintf(fq, "' w impulses\n");
    } else {
	fprintf(fq, "set nokey\n");
	fprintf(fq, "plot '-' using 1:2 w impulses\n");
    }

    /* send data inline */
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif
    for (l=1; l<=m; l++) 
	fprintf(fq, "%d %f\n", l, acf[l]);
    fprintf(fq, "e\n");
    for (l=1; l<=maxlag; l++) 
	fprintf(fq, "%f %f\n", l + .1, pacf[l-1]);
    fprintf(fq, "e\n");
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

#if defined(OS_WIN32) && !defined(GNUPLOT_PNG)
    fprintf(fq, "pause -1\n");
#endif
    fclose(fq);
    err = gnuplot_display(ppaths);

 getout:
    free(x);
    free(y);
    free(acf);
    if (pacf) free(pacf);

    return err;
}

/* ...................................................... */

static int roundup_half (int i)
{
    return (int) ceil((double) i / 2.0);
}

#ifdef notdef
/* Not sure what I was doing here 8-/ */
static int old_roundup_half (int i)
{
    int j;

    j = 10 * (int) ((i / 2) / 10.0 + .5);
    if (j < i / 2) j += 5;
    else if (j - i / 2 > 5) j-= 5;
    return j;
}
#endif

/* ...................................................... */

static int fract_int (int n, double *hhat, double *omega, PRN *prn)
{
    double xx, tstat, **tmpZ = NULL;
    DATAINFO tmpdinfo;
    MODEL tmp;
    int t, err = 0, list[4];

    tmpdinfo.vector = NULL;
    tmpdinfo.n = n;
    tmpdinfo.v = 3;
    tmpdinfo.pd = 1;
    tmpdinfo.extra = 1;
    if (start_new_Z(&tmpZ, &tmpdinfo, 1))
	return 1;
    
    for (t=0; t<n; t++) {
	tmpZ[0][t] = 1.0;
	tmpZ[1][t] = log(hhat[t]);
	xx = sin(omega[t] / 2);
	tmpZ[2][t] = log(4 * xx * xx);
    }

    list[0] = 3;
    list[1] = 1;
    list[2] = 2;
    list[3] = 0;

    _init_model(&tmp, &tmpdinfo);
    tmp = lsq(list, &tmpZ, &tmpdinfo, OLS, 0, 0);

    if (!tmp.errcode) {
	tstat = tmp.coeff[1] / tmp.sderr[1];
	pprintf(prn, "\n%s\n"
		"  %s = %f\n"
		"  %s: t(%d) = %f, %s %.4f\n",
		_("Test for fractional integration"),
		_("Estimated degree of integration"), tmp.coeff[1], 
		_("test statistic"), tmp.dfd, tstat, 
		_("with p-value"), tprob(tstat, tmp.dfd));
    } else err = tmp.errcode;

    clear_model(&tmp, &tmpdinfo);
    free_Z(tmpZ, &tmpdinfo);
    clear_datainfo(&tmpdinfo, CLEAR_FULL);

    return err;
}

/**
 * periodogram:
 * @varno: ID number of variable to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @ppaths: struct containing path information.
 * @batch: if non-zero, don't show gnuplot graph.
 * @opt: if non-zero, use Bartlett lag window for periodogram.
 * @prn: gretl printing struct.
 *
 * Computes and displays the periodogram for the variable specified in @list.
 *
 * Returns: 0 on successful completion, error code on error.
 *
 */

int periodogram (int varno, double ***pZ, const DATAINFO *pdinfo, 
		 PATHS *ppaths, int batch, 
		 int opt, PRN *prn)
{
    double *autocov, *omega, *hhat, *savexx = NULL;
    double xx, yy, varx, w;
    int err = 0, k, xmax, L, nT; 
    int nobs, t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int list[2];
    FILE *fq = NULL;

    list[0] = 1;
    list[1] = varno;
    _adjust_t1t2(NULL, list, &t1, &t2, *pZ, NULL);
    nobs = t2 - t1 + 1;

    if (missvals(&(*pZ)[varno][t1], nobs)) {
	pprintf(prn, "\n%s",
		_("Missing values within sample -- can't do periodogram"));
	return 1;
    }    

    if (nobs < 9) {
	pprintf(prn, "\n%s",
		_("Insufficient observations for periodogram"));
	return 1;
    }
    if (_isconst(t1, t2, &(*pZ)[varno][0])) {
	sprintf(gretl_tmp_str, _("'%s' is a constant"), pdinfo->varname[varno]);
	pprintf(prn, "\n%s\n", gretl_tmp_str);
	return 1;
    }

    /* Chatfield (1996); Greene 4ed, p. 772 */
    if (opt) L = (int) 2.0 * sqrt((double) nobs);
    else L = nobs - 1; 

    /* prepare for fractional integration test */
    nT = (int) sqrt((double) nobs);
    xx = sqrt((double) nobs);
    if ((double) nT < xx) nT += 1;
    
    autocov = malloc((L + 1) * sizeof *autocov);
    omega = malloc(nT * sizeof *omega);
    hhat = malloc(nT * sizeof *hhat);
    if (autocov == NULL || omega == NULL || hhat == NULL) 
	return E_ALLOC;

    xx = _esl_mean(t1, t2, (*pZ)[varno]);
    /* find autocovariances */
    for (k=1; k<=L; k++) {
	autocov[k] = 0;
	for (t=t1+k; t<=t2; t++) {
	    autocov[k] += 
		((*pZ)[varno][t] - xx) * ((*pZ)[varno][t-k] - xx);
	}
	autocov[k] /= nobs;
    }

    xmax = roundup_half(nobs);

    if (!batch && gnuplot_init(ppaths, &fq) == 0) {
	char titlestr[80];

	fprintf(fq, "# periodogram\n");
	fprintf(fq, "set xtics nomirror\n"); 
	if (pdinfo->pd == 4)
	    fprintf(fq, "set x2label '%s'\n", I_("quarters"));
	else if (pdinfo->pd == 12)
	    fprintf(fq, "set x2label '%s'\n", I_("months"));
	else if (pdinfo->pd == 1 && pdinfo->time_series == TIME_SERIES)
	    fprintf(fq, "set x2label '%s'\n", I_("years"));
	else
	    fprintf(fq, "set x2label '%s'\n", I_("periods"));
	fprintf(fq, "set x2range [0:%d]\n", xmax);
	fprintf(fq, "set x2tics(");
	k = (nobs / 2) / 6;
	for (t=1; t<=nobs/2; t += k) {
	    fprintf(fq, "\"%.1f\" %d, ", 
		    (double) (nobs / 2) / (2 * t), t);
	}
	fprintf(fq, "\"\" %d)\n", nobs);
	fprintf(fq, "set xlabel '%s'\n", I_("scaled frequency"));
	fprintf(fq, "set xzeroaxis\n");
	fprintf(fq, "set nokey\n");
	sprintf(titlestr, I_("Spectrum of %s"), pdinfo->varname[varno]);
	fprintf(fq, "set title '%s", titlestr);
	if (opt) {
	    sprintf(titlestr, I_("Bartlett window, length %d"), L);
	    fprintf(fq, " (%s)'\n", titlestr);
	}
	else
	    fprintf(fq, "'\n");
	fprintf(fq, "set xrange [0:%d]\n", xmax);
	fprintf(fq, "plot '-' using 1:2 w lines\n");
    }

    pprintf(prn, _("\nPeriodogram for %s\n"), pdinfo->varname[varno]);
    pprintf(prn, _("Number of observations = %d\n"), nobs);
    if (opt) 
	pprintf(prn, _("Using Bartlett lag window, length %d\n\n"), L);
    pprintf(prn, _(" omega  scaled frequency  periods  spectral density\n\n"));

    if (!batch && fq) savexx = malloc((1 + nobs/2) * sizeof *savexx);

    varx = _esl_variance(t1, t2, &(*pZ)[varno][0]);
    varx *= (double) (nobs - 1) / nobs;
    for (t=1; t<=nobs/2; t++) {
	yy = 2 * M_PI * t / (double) nobs;
	xx = varx; 
	for (k=1; k<=L; k++) {
	    if (opt) w = 1 - (double) k/(L + 1);
	    else w = 1.0;
	    xx += 2.0 * w * autocov[k] * cos(yy * k);
	}
	xx /= 2 * M_PI;
	pprintf(prn, " %.4f%9d%16.2f%14.4f\n", yy, t, 
		(double) (nobs / 2) / (2 * t), xx);
	if (!batch && fq && savexx) savexx[t] = xx;
	if (t <= nT) {
	    omega[t-1] = yy;
	    hhat[t-1] = xx;
	}
    }
    pprintf(prn, "\n");

    if (!batch && fq) {
	if (savexx == NULL) {
	    fclose(fq);
	} else {
#ifdef ENABLE_NLS
	    setlocale(LC_NUMERIC, "C");
#endif
	    for (t=1; t<=nobs/2; t++) fprintf(fq, "%d %f\n", t, savexx[t]);
#ifdef ENABLE_NLS
	    setlocale(LC_NUMERIC, "");
#endif
	    fprintf(fq, "e\n");
#ifdef OS_WIN32
	    fprintf(fq, "pause -1\n");
#endif
	    fclose(fq);
	    free(savexx);
	    err = gnuplot_display(ppaths);
	}
    }

    if (opt == 0 && fract_int(nT, hhat, omega, prn)) {
	pprintf(prn, "\n%s\n",
		_("Fractional integration test failed"));
	err = 1;
    }

    free(autocov);
    free(omega);
    free(hhat);

    return err;
}

/* ............................................................. */

static void printf17 (double zz, PRN *prn)
{
    if (na(zz)) pprintf(prn, "%17s", _("undefined"));
    else {
	pprintf(prn, " ");
	gretl_print_value(zz, prn);
    }
}

/* ............................................................. */

#define LINEWID 78

static void center_line (char *str, PRN *prn, int dblspc)
{
    size_t len = strlen(str);

    if (LINEWID > len) {
	size_t i, pad = (LINEWID - len) / 2;
	char cstr[84];

	for (i=0; i<pad; i++) cstr[i] = ' ';
	strcpy(cstr + i, str);
	if (dblspc) strcat(cstr, "\n");
	pprintf(prn, "%s\n", cstr);
    } else {
	if (dblspc) strcat(str, "\n");
	pprintf(prn, "%s\n", str);
    }
}

/* ............................................................... */

static void prhdr (const char *str, const DATAINFO *pdinfo, 
		   int ci, PRN *prn)
{
    char date1[9], date2[9], tmp[96];

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    pprintf(prn, "\n");

    sprintf(tmp, _("%s, using the observations %s - %s"), str, date1, date2);
    center_line(tmp, prn, 0);

    if (ci == CORR) {
	strcpy(tmp, _("(missing values denoted by -999 will be skipped)"));
	center_line(tmp, prn, 1);
    }
}

/**
 * print_summary:
 * @summ: gretl summary statistics struct.
 * @pdinfo: information on the data set.
 * @pause: if non-zero, pause after showing each screen of info.
 * @prn: gretl printing struct.
 *
 * Print the summary statistics for a given variable.
 *
 */

void print_summary (GRETLSUMMARY *summ,
		    const DATAINFO *pdinfo,
		    int pause, PRN *prn)
{
    double xbar, std, xcv;
    int lo = summ->list[0], v, lv, lineno = 4;
    char tmp[96];

    prhdr(_("Summary Statistics"), pdinfo, SUMMARY, prn);
    if (lo == 1) {
	sprintf(tmp, _("for the variable '%s' (%d valid observations)"), 
		pdinfo->varname[summ->list[1]], summ->n);
	center_line(tmp, prn, 1);
    } else {
	strcpy(tmp, _("(missing values denoted by -999 will be skipped)"));
	center_line(tmp, prn, 1);
	pprintf(prn, "\n%s  ", _("Variable"));
    }

    pprintf(prn, _("             MEAN         MEDIAN            MIN"
            "            MAX\n"));


    for (v=1; v<=lo; v++) {
	if (pause) page_break(1, &lineno, 0);
	lineno++;
	lv = summ->list[v];
	xbar = summ->coeff[v];
	if (lo > 1)
	    pprintf(prn, "%-10s", pdinfo->varname[lv]);
	else _bufspace(2, prn);
	printf17(xbar, prn);
	printf17(summ->xmedian[v], prn);
	printf17(summ->xpx[v], prn);
	printf17(summ->xpy[v], prn);
	pprintf(prn, "\n");
    }

    if (pause) page_break(lo + 2, &lineno, 0);
    lineno += 2;
    pprintf(prn, "\n");

    if (lo > 1) pprintf(prn, "\n%s  ", _("Variable"));
    pprintf(prn, _("             S.D.           C.V.           "
	 "SKEW       EXCSKURT\n"));

    for (v=1; v<=lo; v++) {
	if (pause) page_break(1, &lineno, 0);
	lineno++;
	lv = summ->list[v];
	if (lo > 1)
	    pprintf(prn, "%-10s", pdinfo->varname[lv]);
	else _bufspace(2, prn);
	xbar = summ->coeff[v];
	std = summ->sderr[v];
	if (xbar != 0.0) xcv = (xbar > 0)? std/xbar: (-1) * std/xbar;
	else xcv = -999;
	printf17(std, prn);
	printf17(xcv, prn);
	printf17(summ->xskew[v], prn);
	printf17(summ->xkurt[v], prn);
	pprintf(prn, "\n");
    }
    pprintf(prn, "\n");
}

/**
 * free_summary:
 * @summ: gretl summary statistics struct
 *
 * Frees all malloced elements of the struct.
 *
 */

void free_summary (GRETLSUMMARY *summ)
{
    free(summ->xskew);
    free(summ->xkurt);
    free(summ->xmedian);
    free(summ->coeff);
    free(summ->sderr);
    free(summ->xpx);
    free(summ->xpy); 
    free(summ->list);
    free(summ);
}

/**
 * summary:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 *
 * Calculates descriptive summary statistics for the specified variables.
 *
 * Returns: struct containing the summary statistics.
 *
 */

GRETLSUMMARY *summary (LIST list, 
		       double ***pZ, const DATAINFO *pdinfo,
		       PRN *prn) 
{
    int mm, lo;
    int v, *tmp = NULL;
    GRETLSUMMARY *summ;
    double xbar, std, low, high, skew, kurt, *x = NULL;

    summ = malloc(sizeof *summ);
    if (summ == NULL) return NULL;
    summ->list = NULL;

    lo = list[0];
    mm = lo + 1;

    if ((summ->xskew = malloc(mm * sizeof(double))) == NULL) return NULL;
    if ((summ->xkurt = malloc(mm * sizeof(double))) == NULL) return NULL;
    if ((summ->xmedian = malloc(mm * sizeof(double))) == NULL) return NULL;
    if ((summ->coeff = malloc(mm * sizeof(double))) == NULL) return NULL;
    if ((summ->sderr = malloc(mm * sizeof(double))) == NULL) return NULL;
    if ((summ->xpx = malloc(mm * sizeof(double))) == NULL) return NULL;
    if ((summ->xpy = malloc(mm * sizeof(double))) == NULL) return NULL;
    if ((x = malloc((pdinfo->t2 - pdinfo->t1 + 1) * sizeof *x)) == NULL) 
	return NULL;

    for (v=1; v<=lo; v++)  {
	summ->n = ztox(list[v], x, *pZ, pdinfo);
	if (summ->n < 2) { /* zero or one observations */
	    if (summ->n == 0)
		pprintf(prn, _("Dropping %s: sample range contains no valid "
			"observations\n"), pdinfo->varname[list[v]]);
	    else
		pprintf(prn, _("Dropping %s: sample range has only one "
			"obs, namely %g\n"), pdinfo->varname[list[v]], x[0]);
	    list_exclude(v, list);
	    if (list[0] == 0) {
		free_summary(summ);
		free(x);
		return NULL;
	    } else {
		lo--;
		v--;
		continue;
	    }
	}
	_minmax(0, summ->n-1, x, &low, &high);	
	moments(0, summ->n-1, x, &xbar, &std, &skew, &kurt, 1);
	summ->xpx[v] = low;
	summ->xpy[v] = high;
	summ->coeff[v] = xbar;
	summ->sderr[v] = std;
	summ->xskew[v] = skew;
	summ->xkurt[v] = kurt;
	qsort(x, summ->n, sizeof *x, _compare_doubles); 
	if (summ->n > 1) summ->xmedian[v] = esl_median(x, summ->n);
	else summ->xmedian[v] = x[1];
    }

    copylist(&tmp, list);
    summ->list = tmp;
    free(x);

    return summ;
}

/**
 * free_corrmat:
 * @corrmat: gretl correlation matrix struct
 *
 * Frees all malloced elements of the struct.
 *
 */

void free_corrmat (CORRMAT *corrmat)
{
    if (corrmat != NULL) {
	free(corrmat->list);
	free(corrmat->xpx);
	free(corrmat);
    }
}

/**
 * corrlist:
 * @list: list of variables to process, by ID number.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Computes pairwise correlation coefficients for the variables
 * specified in @list, skipping any constants.
 *
 * Returns: gretl correlation matrix struct.
 * 
 */

CORRMAT *corrlist (LIST list, double ***pZ, const DATAINFO *pdinfo)
{
    CORRMAT *corrmat;
    int *p = NULL;
    int i, j, lo, nij, mm;
    int t1 = pdinfo->t1, t2 = pdinfo->t2; 

    corrmat = malloc(sizeof *corrmat);
    if (corrmat == NULL) return NULL;

    copylist(&p, list);
    if (p == NULL) {
	free(corrmat);
	return NULL;
    }	

    /* drop any constants from list */
    for (i=1; i<=p[0]; i++) {
	if (_isconst(t1, t2, (*pZ)[p[i]])) {
	    list_exclude(i, p);
	    i--;
	}
    }
    corrmat->list = p;

    lo = corrmat->list[0];  
    corrmat->n = t2 - t1 + 1;
    mm = (lo * (lo + 1))/2;
    corrmat->xpx = malloc(mm * sizeof(double));
    if (corrmat->xpx == NULL) {
	free_corrmat(corrmat);
	return NULL;
    }
    for (i=1; i<=lo; i++) {   
	for (j=i; j<=lo; j++)  {
	    nij = ijton(i, j, lo);
	    if (i == j) {
		corrmat->xpx[nij] = 1.0;
		continue;
	    }
	    corrmat->xpx[nij] = _corr(corrmat->n, 
				      &(*pZ)[corrmat->list[i]][t1],
				      &(*pZ)[corrmat->list[j]][t1]);
	}
    }

    corrmat->t1 = t1;
    corrmat->t2 = t2;

    return corrmat;
}

/**
 * matrix_print_corr:
 * @corr: gretl correlation matrix
 * @pdinfo: data information struct.
 * @pause: = 1 to pause after showing each screen of info.
 * @prn: gretl printing struct.
 *
 * Prints a gretl correlation matrix.
 *
 */

void matrix_print_corr (CORRMAT *corr, const DATAINFO *pdinfo,
			int pause, PRN *prn)
{
    char tmp[96];

    prhdr(_("Correlation Coefficients"), pdinfo, CORR, prn);
    sprintf(tmp, _("5%% critical value (two-tailed) = "
	    "%.4f for n = %d"), rhocrit95(corr->n), corr->n);
    center_line(tmp, prn, 1);
    _mxout(corr->xpx, corr->list, CORR, pdinfo, pause, prn);
}

/**
 * esl_corrmx:
 * @list: gives the ID numbers of the variables to process.
 * @pZ: pointer to the data matrix.
 * @pdinfo: data information struct.
 * @pause: if non-zero, pause after showing each screen of info.
 * @prn: gretl printing struct.
 *
 * Computes and prints the correlation matrix for the specified list
 * of variables.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int esl_corrmx (LIST list, double ***pZ, const DATAINFO *pdinfo, 
		int pause, PRN *prn)
{
    CORRMAT *corr;

    corr = corrlist(list, pZ, pdinfo);
    if (corr == NULL) return 1;
    matrix_print_corr(corr, pdinfo, pause, prn);
    free_corrmat(corr);
    return 0;
}

/**
 * means_test:
 * @list: gives the ID numbers of the variables to compare.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @vareq: assume population variances are equal (1) or not (0).
 * @prn: gretl printing struct.
 *
 * Carries out test of the null hypothesis that the means of two
 * variables are equal.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int means_test (LIST list, double **Z, const DATAINFO *pdinfo, 
		int vareq, PRN *prn)
{
    double m1, m2, s1, s2, skew, kurt, se, mdiff, t, pval;
    double *x = NULL, *y = NULL;
    int df, n1, n2, n = pdinfo->n;

    if (list[0] < 2) return E_ARGS;

    if ((x = malloc(n * sizeof *x)) == NULL) return E_ALLOC;
    if ((y = malloc(n * sizeof *y)) == NULL) return E_ALLOC;

    n1 = ztox(list[1], x, Z, pdinfo);
    n2 = ztox(list[2], y, Z, pdinfo);
    if (n1 == 0 || n2 == 0) {
	pprintf(prn, _("Sample range has no valid observations."));
	free(x); free(y);
	return 1;
    }
    if (n1 == 1 || n2 == 1) {
	pprintf(prn, _("Sample range has only one observation."));
	free(x); free(y);
	return 1;
    }
    df = n1 + n2 - 2;

    moments(0, n1-1, x, &m1, &s1, &skew, &kurt, 1);
    moments(0, n2-1, y, &m2, &s2, &skew, &kurt, 1);
    mdiff = m1 - m2;

    if (vareq) {
	double sp2;

	sp2 = ((n1-1)*s1*s1 + (n2-1)*s2*s2) / df;
	se = sqrt(sp2/n1 + sp2/n2);
    } else 
	se = sqrt((s1*s1/n1) + (s2*s2/n2));

    t = mdiff / se;
    pval = tprob(t, df);

    pprintf(prn, _("\nEquality of means test "
	    "(assuming %s variances)\n\n"), (vareq)? _("equal") : _("unequal"));
    pprintf(prn, _("   Difference between sample means = %g - %g = %g\n"), 
	    m1, m2, mdiff);
    pprintf(prn, _("   Null hypothesis: The two population means are the same.\n"));
    pprintf(prn, _("   Estimated standard error = %g\n"), se);
    pprintf(prn, _("   Test statistic: t(%d) = %g\n"), df, t);
    pprintf(prn, _("   p-value (two-tailed) = %g\n\n"), pval);
    if (pval > .10)
	pprintf(prn, _("   The difference is not statistically significant.\n\n"));

    free(x);
    free(y);

    return 0;
}

/**
 * vars_test:
 * @list: gives the ID numbers of the variables to compare.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Carries out test of the null hypothesis that the variances of two
 * variables are equal.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int vars_test (LIST list, double **Z, const DATAINFO *pdinfo, 
	       PRN *prn)
{
    double m, skew, kurt, s1, s2, var1, var2, F;
    double *x = NULL, *y = NULL;
    int dfn, dfd, n1, n2, n = pdinfo->n;

    if (list[0] < 2) return E_ARGS;

    if ((x = malloc(n * sizeof *x)) == NULL) return E_ALLOC;
    if ((y = malloc(n * sizeof *y)) == NULL) return E_ALLOC;

    n1 = ztox(list[1], x, Z, pdinfo);
    n2 = ztox(list[2], y, Z, pdinfo);
    if (n1 == 0 || n2 == 0) {
	pprintf(prn, _("Sample range has no valid observations."));
	free(x); free(y);
	return 1;
    }
    if (n1 == 1 || n2 == 1) {
	pprintf(prn, _("Sample range has only one observation."));
	free(x); free(y);
	return 1;
    }
    
    moments(0, n1-1, x, &m, &s1, &skew, &kurt, 1);
    moments(0, n2-1, y, &m, &s2, &skew, &kurt, 1);

    var1 = s1*s1;
    var2 = s2*s2;
    if (var1 > var2) { 
	F = var1/var2;
	dfn = n1 - 1;
	dfd = n2 - 1;
    } else {
	F = var2/var1;
	dfn = n2 - 1;
	dfd = n1 - 1;
    }

    pprintf(prn, _("\nEquality of variances test\n\n"));
    pprintf(prn, _("   Ratio of sample variances = %g\n"), F);
    pprintf(prn, _("   Null hypothesis: The two population variances are equal.\n"));
    pprintf(prn, _("   Test statistic: F(%d,%d) = %g\n"), dfn, dfd, F);
    pprintf(prn, _("   p-value (two-tailed) = %g\n\n"), fdist(F, dfn, dfd));
    if (fdist(F, dfn, dfd) > .10)
	pprintf(prn, _("   The difference is not statistically significant.\n\n"));

    free(x);
    free(y);

    return 0;

}
