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

#ifndef PI
#define PI 3.14159265358979323846
#endif

extern void _mxout (const double *rr, const int *list, const int ci,
		    const DATAINFO *pdinfo, const int batch, print_t *prn);

/* ........................................................... */

static int missvals (double *x, int n)
{
    int t;
    
    for (t=0; t<n; t++)
	if (na(x[t])) return 1;
	    
    return 0;
}

/* ........................................................... */

double esl_median (const double *zx, const int n)
/* if zx is in increasing order, finds median for observations 0
   to n-1 */
{
    double xx;
    int n2, n2p;

    n2p = (n2 = n/2) + 1;
    xx = (n % 2)? zx[n2p - 1] : 0.5 * (zx[n2 - 1] + zx[n2p - 1]);
    return xx;
}

/* ........................................................... */

static void moments (const int t1, const int t2, const double *zx, 
		     double *xbar, double *std, 
		     double *skew, double *kurt, const int k)
     /* k is the "degrees of freedom loss": it will generally be one,
	other than when dealing with a regression residual */
{
    int t, n = t2 - t1 + 1;
    double sum = 0.0, p, var;

    if (isconst(t1, t2, zx)) {
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

/* ........................................................ */

void free_freq (FREQDIST *freq)
{
    free(freq->midpt);
    free(freq->endpt);
    free(freq->f);
    free(freq);
}

/* ........................................................ */

FREQDIST *freq_func (double **pZ, const DATAINFO *pdinfo, double *zz,
		     const int nzz, const char *varname, const int params)
     /* generates frequency distribution:
	params is the "degrees of freedom loss", which will
	genrally be 1 unless dealing with a regression residual. 
        if pZ is NULL, use the passed-in data zz, otherwise look
        up the variable by name and find its values in pZ.
     */
{
    FREQDIST *freq;
    double *x = NULL;
    double xx, xmin, xmax, xbar, sdx;
    double skew, kurt;
    int i, k, n, t;
    int maxend = 16;

    freq = malloc(sizeof *freq);
    if (freq == NULL) return NULL;

    x = malloc((pdinfo->t2 - pdinfo->t1 + 1) * sizeof *x);
    if (x == NULL) {
	free(freq);
	return NULL;
    }

    freq->errcode = 0;
    freq->errmsg[0] = '\0';
    freq->midpt = NULL;
    freq->endpt = NULL;
    freq->f = NULL;

    if (pZ != NULL) {
	i = varindex(pdinfo, varname);
	if (i > pdinfo->v - 1) {
	    freq->errcode = E_DATA;
	    sprintf(freq->errmsg, "'%s' is not in the data set", varname);
	    free(x);
	    return freq;
	}	
	n = ztox(i, x, pdinfo, *pZ);
	if (n < 3) {
	    freq->errcode = E_DATA;
	    sprintf(freq->errmsg, "Insufficient data to build frequency "
		    "distribution for variable %s", varname);
	    free(x);
	    return freq;
	}
    } else {
	n = nzz;
	x = zz;
    }

    strcpy(freq->varname, varname);
    freq->t1 = pdinfo->t1; 
    freq->t2 = pdinfo->t2;

    if (isconst(0, n-1, x)) {
	freq->errcode = 1;
	sprintf(freq->errmsg, "%s is a constant", freq->varname);
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
	freq->errcode = E_ALLOC;
	strcpy(freq->errmsg, "Out of memory for frequency distribution");
	free(x);
	return freq;
    }
    
    minmax(0, n-1, x, &xmin, &xmax);
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

static int get_pacf (double *pacf, int *maxlag, const int varnum, 
		     double **pZ, DATAINFO *pdinfo)
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
	laggenr(varnum, i, 0, pZ, pdinfo);
	/* the lagvar may already exist */
	laglist[i-1] = lagvarnum(varnum, i, pdinfo); 
    }

    init_model(&tmp);
    pdinfo->t1 = t1;

    list[1] = varnum;
    for (i=2; i<=*maxlag; i++) {
	list[0] = i + 2;
	list[i+2] = 0;
	for (j=2; j<i+2; j++) list[j] = laglist[j-2];
	tmp = lsq(list, *pZ, pdinfo, OLS, 0, 0);
	if ((err = tmp.errcode)) break;
	pacf[i-1] = tmp.coeff[i];
	if (i < *maxlag) clear_model(&tmp, NULL, NULL);
    }

    clear_model(&tmp, NULL, NULL);
    shrink_Z(pdinfo->v - v, pZ, pdinfo);
    free(laglist);
    free(list);

    return err;
}

/* ...................................................... */

int corrgram (const int *list, const int order, double **pZ, 
	      DATAINFO *pdinfo, const PATHS *ppaths, 
	      const int batch, print_t *prn)
/* computes values of autocorrelation function and plots correlogram */
{
    double *x, *y, *acf, *xl, box;
    double *pacf = NULL;
    int err = 0, k, l, m, v = list[1], nobs, n = pdinfo->n; 
    int maxlag = 0, t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    FILE *fq;

    adjust_t1t2(NULL, list, &t1, &t2, *pZ, pdinfo->n, NULL);
    nobs = t2 - t1 + 1;

    if (missvals(&(*pZ)[v*n + t1], nobs)) {
	pprintf(prn, "\nMissing values within sample -- can't do correlogram");
	return 1;
    }

    if (nobs < 4) {
	pprintf(prn, "\nInsufficient observations for correlogram");
	return 1;
    }
    if (isconst(t1, t2, &(*pZ)[v*n])) {
	pprintf(prn, "\n'%s' is a constant\n", pdinfo->varname[v]);
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
	    x[k] = (*pZ)[v*n + t];
	    y[k] = (*pZ)[v*n + t-l];
	}
	acf[l] = corr(nobs-l, x, y);
    }
    pprintf(prn, "\nAutocorrelation function for %s\n\n", pdinfo->varname[v]);

    /* add Box-Pierce statistic */
    box = 0;
    for (t=1; t<=m; t++) 
	box += acf[t] * acf[t];
    box *= m;
    pprintf(prn, "Box-Pierce Q statistic = %.4f\n", box);
    pprintf(prn, "Degrees of freedom = %d, significance level = %.4f\n\n",
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
        pprintf(prn, "\n\nCorrelogram\n\n");
	graphyzx(NULL, acf + 1, NULL, xl, m, pdinfo->varname[v], 
		 "lag", NULL, 0, prn);
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
    err = get_pacf(pacf, &maxlag, v, pZ, pdinfo);
    pacf[0] = acf[1];
    if (!err) {
	pprintf(prn, "\nPartial autocorrelations");
	if (maxlag < m) 
	    pprintf(prn, " (to lag %d):\n\n", maxlag);
	else
	    pprintf(prn, ":\n\n");
	for (l=1; l<=maxlag; l++) {
	    pprintf(prn, "%5d)%7.3f", l, pacf[l-1]);
	    if (l%5 == 0) pprintf(prn, "\n");
	}
    }
    pprintf(prn, "\n");
    if (maxlag%5 > 0) pprintf(prn, "\n");

    fq = fopen(ppaths->plotfile, "w");
    /*  fq = popen("gnuplot -persist", "w"); */
    if (fq == NULL) return E_FOPEN;
    fprintf(fq, "# correlogram\n");
    fprintf(fq, "set xlabel \"lag\"\n");
    fprintf(fq, "set xzeroaxis\n");
    fprintf(fq, "set title \"Correlogram for %s\"\n", pdinfo->varname[v]);
    if (maxlag) {
	fprintf(fq, "plot '-' using 1:2 title 'autocorrelations' "
		"w impulses, \\\n"
		"'-' using 1:2 title 'partial autocorrelations");
	if (maxlag < m)
	    fprintf(fq, "(to lag %d)' w impulses\n", maxlag);
	else
	    fprintf(fq, "' w impulses\n");
    } else {
	fprintf(fq, "set nokey\n");
	fprintf(fq, "plot '-' using 1:2 w impulses\n");
    }
    /* send data inline */
    for (l=1; l<=m; l++) 
	fprintf(fq, "%d %f\n", l, acf[l]);
    fprintf(fq, "e\n");
    for (l=1; l<=maxlag; l++) 
	fprintf(fq, "%f %f\n", l + .1, pacf[l-1]);
    fprintf(fq, "e\n");

#ifdef OS_WIN32
    fprintf(fq, "pause -1\n");
#endif
    fclose(fq);
    err = gnuplot_display(ppaths->gnuplot, ppaths->plotfile);
    /*  pclose(fq); */

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
    int j;

    j = 10 * (int) ((i / 2) / 10.0 + .5);
    if (j < i / 2) j += 5;
    else if (j - i / 2 > 5) j-= 5;
    return j;
}

/* ...................................................... */

static int fract_int (int n, double *hhat, double *omega, print_t *prn)
{
    double xx, tstat, *tmpZ = NULL;
    DATAINFO tmpdinfo;
    MODEL tmp;
    int t, err = 0, list[4];

    tmpdinfo.S = NULL;
    tmpdinfo.n = n;
    tmpdinfo.v = 3;
    tmpdinfo.pd = 1;
    tmpdinfo.extra = 1;
    if (start_new_Z(&tmpZ, &tmpdinfo, 1))
	return 1;
    
    for (t=0; t<n; t++) {
	tmpZ[t] = 1.0;
	tmpZ[n + t] = log(hhat[t]);
	xx = sin(omega[t] / 2);
	tmpZ[2*n + t] = log(4 * xx * xx);
    }

    list[0] = 3;
    list[1] = 1;
    list[2] = 2;
    list[3] = 0;

    init_model(&tmp);
    tmp = lsq(list, tmpZ, &tmpdinfo, OLS, 0, 0);

    if (!tmp.errcode) {
	tstat = tmp.coeff[1] / tmp.sderr[1];
	pprintf(prn, "\nTest for fractional integration\n"
		"  Estimated degree of integration = %f\n"
		"  test statistic: t(%d) = %f, with p-value %.4f\n",
		tmp.coeff[1], tmp.dfd, tstat, tprob(tstat, tmp.dfd));
    } else err = tmp.errcode;

    clear_model(&tmp, NULL, NULL);
    free(tmpZ);
    clear_datainfo(&tmpdinfo, 0);

    return err;
}

/* ...................................................... */

int periodogram (const int *list, double **pZ, const DATAINFO *pdinfo, 
		 const PATHS *ppaths, const int batch, 
		 const int opt, print_t *prn)
/* compute and display periodogram */
{
    double *autocov, *omega, *hhat;
    double xx, yy, varx, w;
    int err = 0, k, xmax, L, v = list[1], n = pdinfo->n, nT; 
    int nobs, t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    FILE *fq = NULL;
    
    adjust_t1t2(NULL, list, &t1, &t2, *pZ, pdinfo->n, NULL);
    nobs = t2 - t1 + 1;

    if (missvals(&(*pZ)[v*n + t1], nobs)) {
	pprintf(prn, "\nMissing values within sample -- can't do periodogram");
	return 1;
    }    

    if (nobs < 9) {
	pprintf(prn, "\nInsufficient observations for periodogram");
	return 1;
    }
    if (isconst(t1, t2, &(*pZ)[v*n])) {
	pprintf(prn, "\n'%s' is a constant\n", pdinfo->varname[v]);
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

    xx = esl_mean(t1, t2, &(*pZ)[v*n]);
    /* find autocovariances */
    for (k=1; k<=L; k++) {
	autocov[k] = 0;
	for (t=t1+k; t<=t2; t++) {
	    autocov[k] += 
		((*pZ)[v*n + t] - xx) * ((*pZ)[v*n + t-k] - xx);
	}
	autocov[k] /= nobs;
    }

    xmax = roundup_half(nobs);

    if (!batch) {
	fq = fopen(ppaths->plotfile, "w");
	fprintf(fq, "# periodogram\n");
	fprintf(fq, "set xtics nomirror\n"); 
	if (pdinfo->pd == 4)
	    fprintf(fq, "set x2label 'quarters'\n");
	else if (pdinfo->pd == 12)
	    fprintf(fq, "set x2label 'months'\n");
	else if (pdinfo->pd == 1 && pdinfo->time_series == 1)
	    fprintf(fq, "set x2label 'years'\n");
	else
	    fprintf(fq, "set x2label 'periods'\n");
	fprintf(fq, "set x2range [0:%d]\n", xmax);
	fprintf(fq, "set x2tics(");
	k = (nobs / 2) / 6;
	for (t=1; t<=nobs/2; t += k) {
	    fprintf(fq, "\"%.1f\" %d, ", 
		    (double) (nobs / 2) / (2 * t), t);
	}
	fprintf(fq, "\"\" %d)\n", nobs);
	fprintf(fq, "set xlabel 'scaled frequency'\n");
	fprintf(fq, "set xzeroaxis\n");
	fprintf(fq, "set nokey\n");
	fprintf(fq, "set title 'Spectrum of %s", pdinfo->varname[v]);
	if (opt) 
	    fprintf(fq, " (Bartlett window, length %d)'\n", L);
	else
	    fprintf(fq, "'\n");
	fprintf(fq, "set xrange [0:%d]\n", xmax);
	fprintf(fq, "plot '-' using 1:2 w lines\n");
    }

    pprintf(prn, "\nPeriodogram for %s\n", pdinfo->varname[v]);
    pprintf(prn, "Number of observations = %d\n", nobs);
    if (opt) 
	pprintf(prn, "Using Bartlett lag window, length %d\n\n", L);
    pprintf(prn, " omega  scaled frequency  periods  spectral density\n\n");

    varx = esl_variance(t1, t2, &(*pZ)[v*n]);
    varx *= (double) (nobs - 1) / nobs;
    for (t=1; t<=nobs/2; t++) {
	yy = 2 * PI * t / (double) nobs;
	xx = varx; 
	for (k=1; k<=L; k++) {
	    if (opt) w = 1 - (double) k/(L + 1);
	    else w = 1.0;
	    xx += 2.0 * w * autocov[k] * cos(yy * k);
	}
	xx /= 2 * PI;
	pprintf(prn, " %.4f%9d%16.2f%14.4f\n", yy, t, 
		(double) (nobs / 2) / (2 * t), xx);
	if (!batch) fprintf(fq, "%d %f\n", t, xx);
	if (t <= nT) {
	    omega[t-1] = yy;
	    hhat[t-1] = xx;
	}
    }
    pprintf(prn, "\n");

    if (!batch) {
#ifdef OS_WIN32
	fprintf(fq, "pause -1\n");
#endif
	fclose(fq);
	err = gnuplot_display(ppaths->gnuplot, ppaths->plotfile);
    }

    if (opt == 0 && fract_int(nT, hhat, omega, prn)) {
	pprintf(prn, "\nFractional integration test failed\n");
	err = 1;
    }

    free(autocov);
    free(omega);
    free(hhat);

    return err;
}

/* ............................................................. */

static void printf17 (const double zz, print_t *prn)
{
    if (na(zz)) pprintf(prn, "      undefined");
    else printxs(zz, 17, SUMMARY, prn);
}

/* ............................................................... */

static void prhdr (const char *str, const DATAINFO *pdinfo, 
		   const int ci, print_t *prn)
{
    char date1[9], date2[9];

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    pprintf(prn, "\n");
    if (pdinfo->pd != 1) space((ci == SUMMARY)? 10 : 7, prn);
    else {
	if (pdinfo->sd0 > 1900) space((ci == SUMMARY)? 12 : 9, prn);
	else space((ci == SUMMARY)? 15 : 12, prn);
    } 
    
    pprintf(prn, "%s, using the observations %s - %s\n", str, date1, date2);
    if (ci == CORR) {
	pprintf(prn, "               "
		"(missing values denoted by -999 will be skipped)\n\n");
    }
}

/* ............................................................. */

void print_summary (GRETLSUMMARY *summ,
		    const DATAINFO *pdinfo,
		    print_t *prn, int batch)
{
    double xbar, std, xcv;
    int lo = summ->list[0], v, lv, lineno = 4;

    prhdr("Summary Statistics", pdinfo, SUMMARY, prn);
    if (lo == 1) {
	space(16, prn);
	pprintf(prn, "for the variable '%s' (%d valid observations)\n\n", 
		pdinfo->varname[summ->list[1]], summ->n);
    } else {
	pprintf(prn, "               "
		"(missing values denoted by -999 will be skipped)\n\n");
	pprintf(prn, "\nVariable    ");
    }

    pprintf(prn, "             MEAN         MEDIAN            MIN"
            "            MAX\n");
    for (v=1; v<=lo; v++) {
	_pgbreak(1, &lineno, batch);
	lineno++;
	lv = summ->list[v];
	xbar = summ->coeff[v];
	if (lo > 1)
	    pprintf(prn, "%-14s", pdinfo->varname[lv]);
	else space(2, prn);
	printf17(xbar, prn);
	printf17(summ->xmedian[v], prn);
	printf17(summ->xpx[v], prn);
	printf17(summ->xpy[v], prn);
	pprintf(prn, "\n");
    }
    _pgbreak(lo + 2, &lineno, batch);
    lineno += 2;
    pprintf(prn, "\n");
    if (lo > 1) pprintf(prn, "\nVariable    ");
    pprintf(prn, "             S.D.           C.V.           "
	 "SKEW       EXCSKURT\n");
    for (v=1; v<=lo; v++) {
	_pgbreak(1, &lineno, batch);
	lineno++;
	lv = summ->list[v];
	if (lo > 1)
	    pprintf(prn, "%-14s", pdinfo->varname[lv]);
	else space(2, prn);
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

/* ............................................................. */

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

/* ............................................................. */

GRETLSUMMARY *summary (int *list, 
		       double **pZ, const DATAINFO *pdinfo,
		       print_t *prn) 
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
	summ->n = ztox(list[v], x, pdinfo, *pZ);
	if (summ->n < 2) { /* zero or one observations */
	    if (summ->n == 0)
		pprintf(prn, "Dropping %s: sample range contains no valid "
			"observations\n", pdinfo->varname[list[v]]);
	    else
		pprintf(prn, "Dropping %s: sample range has only one "
			"obs, namely %g\n", pdinfo->varname[list[v]], x[0]);
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
	minmax(0, summ->n-1, x, &low, &high);	
	moments(0, summ->n-1, x, &xbar, &std, &skew, &kurt, 1);
	summ->xpx[v] = low;
	summ->xpy[v] = high;
	summ->coeff[v] = xbar;
	summ->sderr[v] = std;
	summ->xskew[v] = skew;
	summ->xkurt[v] = kurt;
	qsort(x, summ->n, sizeof *x, compare_doubles); 
	if (summ->n > 1) summ->xmedian[v] = esl_median(x, summ->n);
	else summ->xmedian[v] = x[1];
    }

    copylist(&tmp, list);
    summ->list = tmp;
    free(x);

    return summ;
}

/* ....................................................... */

void free_corrmat (CORRMAT *corrmat)
{
    if (corrmat != NULL) {
	free(corrmat->list);
	free(corrmat->xpx);
	free(corrmat);
    }
}

/* ....................................................... */

CORRMAT *corrlist (int *list, double **pZ, const DATAINFO *pdinfo)
/* computes pairwise correlation coefficients for 
   variables in list, skipping any constants, from
   observation pdinfo->t1 to pdinfo->t2.  
*/
{
    CORRMAT *corrmat;
    int *p = NULL;
    int i, j, lo, nij, mm;
    int t1 = pdinfo->t1, t2 = pdinfo->t2, n = pdinfo->n; 

    corrmat = malloc(sizeof *corrmat);
    if (corrmat == NULL) return NULL;

    copylist(&p, list);
    if (p == NULL) {
	free(corrmat);
	return NULL;
    }	

    /* drop any constants from list */
    for (i=1; i<=p[0]; i++) {
	if (isconst(t1, t2, &(*pZ)[n * p[i]])) {
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
	    corrmat->xpx[nij] = corr(corrmat->n, 
				     &(*pZ)[n * corrmat->list[i] + t1],
				     &(*pZ)[n * corrmat->list[j] + t1]);
	}
    }

    corrmat->t1 = t1;
    corrmat->t2 = t2;

    return corrmat;
}

/* ............................................................ */

void matrix_print_corr (CORRMAT *corr, const DATAINFO *pdinfo,
			const int batch, print_t *prn)
{
    prhdr("Correlation Coefficients", pdinfo, CORR, prn);
    pprintf(prn, "              5%% critical value (two-tailed) = "
	    "%.3f for n = %d\n\n", rhocrit95(corr->n), corr->n);

    _mxout(corr->xpx, corr->list, CORR, pdinfo, batch, prn);
}

/* ............................................................ */

int esl_corrmx (int *list, double **pZ, const DATAINFO *pdinfo, 
		const int batch, print_t *prn)
{
    CORRMAT *corr;

    corr = corrlist(list, pZ, pdinfo);
    if (corr == NULL) return 1;
    matrix_print_corr(corr, pdinfo, batch, prn);
    free_corrmat(corr);
    return 0;
}

/* ............................................................ */

int means_test (int *list, double *Z, const DATAINFO *pdinfo, 
		const int vareq, print_t *prn)
{
    double m1, m2, s1, s2, skew, kurt, se, mdiff, t, pval;
    double *x = NULL, *y = NULL;
    int df, n1, n2, n = pdinfo->n;

    if (list[0] < 2) return E_ARGS;

    if ((x = malloc(n * sizeof *x)) == NULL) return E_ALLOC;
    if ((y = malloc(n * sizeof *y)) == NULL) return E_ALLOC;

    n1 = ztox(list[1], x, pdinfo, Z);
    n2 = ztox(list[2], y, pdinfo, Z);
    if (n1 == 0 || n2 == 0) {
	pprintf(prn, "Sample range has no valid observations.");
	free(x); free(y);
	return 1;
    }
    if (n1 == 1 || n2 == 1) {
	pprintf(prn, "Sample range has only one observation.");
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

    pprintf(prn, "\nEquality of means test "
	    "(assuming %s variances)\n\n", (vareq)? "equal" : "unequal");
    pprintf(prn, "   Difference between sample means = %g - %g = %g\n", 
	    m1, m2, mdiff);
    pprintf(prn, "   Null hypothesis: The two population means are the same.\n");
    pprintf(prn, "   Estimated standard error = %g\n", se);
    pprintf(prn, "   Test statistic: t(%d) = %g\n", df, t);
    pprintf(prn, "   p-value (two-tailed) = %g\n\n", pval);
    if (pval > .10)
	pprintf(prn, "   The difference is not statistically significant.\n\n");

    free(x);
    free(y);

    return 0;
}

/* ............................................................ */

int vars_test (int *list, double *Z, const DATAINFO *pdinfo, 
	       print_t *prn)
{
    double m, skew, kurt, s1, s2, var1, var2, F;
    double *x = NULL, *y = NULL;
    int dfn, dfd, n1, n2, n = pdinfo->n;

    if (list[0] < 2) return E_ARGS;

    if ((x = malloc(n * sizeof *x)) == NULL) return E_ALLOC;
    if ((y = malloc(n * sizeof *y)) == NULL) return E_ALLOC;

    n1 = ztox(list[1], x, pdinfo, Z);
    n2 = ztox(list[2], y, pdinfo, Z);
    if (n1 == 0 || n2 == 0) {
	pprintf(prn, "Sample range has no valid observations.");
	free(x); free(y);
	return 1;
    }
    if (n1 == 1 || n2 == 1) {
	pprintf(prn, "Sample range has only one observation.");
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

    pprintf(prn, "\nEquality of variances test\n\n");
    pprintf(prn, "   Ratio of sample variances = %g\n", F);
    pprintf(prn, "   Null hypothesis: The two population variances are equal.\n");
    pprintf(prn, "   Test statistic: F(%d,%d) = %g\n", dfn, dfd, F);
    pprintf(prn, "   p-value (two-tailed) = %g\n\n", fdist(F, dfn, dfd));
    if (fdist(F, dfn, dfd) > .10)
	pprintf(prn, "   The difference is not statistically significant.\n\n");

    free(x);
    free(y);

    return 0;

}
