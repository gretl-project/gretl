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
#include "gretl_private.h"
#include "gretl_matrix.h"

#include <unistd.h>

#ifdef WIN32
# include <windows.h>
#endif

#if defined(ENABLE_NLS) && defined(USE_GTK2)
#include <glib.h>
#endif

/* ........................................................... */

static int missvals (double *x, int n)
{
    int t;
    
    for (t=0; t<n; t++)
	if (na(x[t])) return 1;
	    
    return 0;
}

/**
 * gretl_median:
 * @x: data series.
 * @n: length of the series.
 *
 * Returns: the median value of the given series.
 *
 */

double gretl_median (const double *x, int n)
{
    double *sx, med;
    int t, n2, n2p;

    sx = malloc(n * sizeof *sx);
    if (sx == NULL) return NADBL;

    for (t=0; t<n; t++) {
	sx[t] = x[t];
    }

    qsort(sx, n, sizeof *sx, gretl_compare_doubles); 

    n2p = (n2 = n/2) + 1;
    med = (n % 2)? sx[n2p - 1] : 0.5 * (sx[n2 - 1] + sx[n2p - 1]);

    free(sx);

    return med;
}

/* ........................................................... */

int moments (int t1, int t2, const double *zx, 
	     double *xbar, double *std, 
	     double *skew, double *kurt, int k)
     /* k is the "degrees of freedom loss": it will generally be one,
	other than when dealing with a regression residual */
{
    int t, n = t2 - t1 + 1;
    double dev, var;
    double s, s2, s3, s4;
    int allstats = 1;

    if (skew == NULL && kurt == NULL) allstats = 0;

    if (gretl_isconst(t1, t2, zx)) {
	*xbar = zx[t1];
	*std = 0.0;
	if (allstats) {
	    *skew = *kurt = NADBL;
	}
	return 1;
    }

    s = 0.0;
    for (t=t1; t<=t2; t++) s += zx[t];

    *xbar = s / n;
    var = 0.0;
    if (allstats) *skew = *kurt = 0.0;

    s2 = s3 = s4 = 0.0;
    for (t=t1; t<=t2; t++) {
	dev = zx[t] - *xbar;
	s2 += dev * dev;
	if (allstats) {
	    s3 += pow(dev, 3);
	    s4 += pow(dev, 4);
	}
    }

    var = s2 / (n-k);

    if (var < 0.0) {
	*std = NADBL;
	if (allstats) *skew = *kurt = NADBL;
	return 1;
    }

    *std = sqrt(var);

    if (allstats) {
	if (var > 0.0) {
	    *skew = (s3 / n) / pow(s2 / n, 1.5);
	    *kurt = ((s4 / n) / pow(s2 / n, 2)) - 3.0; /* excess kurtosis */
	} else {
	    *skew = *kurt = NADBL;
	}
    }

    return 0;
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

static double rb1_to_z1 (double rb1, double n)
{
    double b, w2, d, y, z1;

    b = 3.0 * (n*n + 27*n - 70) * (n+1) * (n+3) /
	((n-2) * (n+5) * (n+7) * (n+9));

    w2 = -1.0 + sqrt(2 * (b-1));
    
    d = 1.0 / sqrt(log(sqrt(w2)));

    y = rb1 * sqrt(((w2-1.0)/2.0) * ((n+1.0)*(n+3.0))/(6.0*(n-2)));

    z1 = d * log(y + sqrt(y*y + 1));

    return z1;
}

static double b2_to_z2 (double b1, double b2, double n)
{
    double d, a, c, k, alpha, chi, z2;
    double n2 = n * n;

    d = (n-3) * (n+1) * (n2 + 15*n - 4.0);

    a = ((n-2) * (n+5) * (n+7) * (n2 + 27*n - 70.0)) / (6.0 * d);

    c = ((n-7) * (n+5) * (n+7) * (n2 + 2*n - 5.0)) / (6.0 * d);

    k = ((n+5) * (n+7) * (n2*n + 37*n2 + 11*n - 313.0)) / (12.0 * d);

    alpha = a + b1 * c;

    chi = (b2 - 1.0 - b1) * 2.0 * k;

    z2 = (pow(chi/(2*alpha), 1.0/3.0) - 1.0 + (1.0 / (9.0*alpha))) *
	sqrt(9.0*alpha);

    return z2;
}

double doornik_chisq (double skew, double kurt, int n)
     /* Bowman-Shenton as modified by Doornik & Hansen */
	
{
    double rb1, b1, b2, z1, z2;

    rb1 = skew;
    b1 = skew * skew;
    b2 = kurt + 3.0; /* convert from "excess" to regular */

    z1 = rb1_to_z1 (rb1, (double) n);
    z2 = b2_to_z2 (b1, b2, (double) n);

    return z1*z1 + z2*z2;
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

#if 0

FREQDIST *old_freqdist (double ***pZ, const DATAINFO *pdinfo, 
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
    if (n < 8) {
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

    freq->endpt = malloc((maxend + 1) * sizeof *freq->endpt);
    freq->midpt = malloc(maxend * sizeof *freq->midpt);
    freq->f = malloc(maxend * sizeof *freq->f);
    if (freq->endpt == NULL || freq->midpt == NULL ||
	freq->f == NULL) {
	gretl_errno = E_ALLOC;
	strcpy(gretl_errmsg, _("Out of memory for frequency distribution"));
	free(x);
	return freq;
    }
    
    gretl_minmax(0, n-1, x, &xmin, &xmax);
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
	freq->midpt[k] = (freq->endpt[k] + freq->endpt[k+1]) / 2.0;
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
	for (k=1; k<freq->numbins; k++) {
	    if (freq->endpt[k] <= xx && xx <freq->endpt[k+1]) {
		freq->f[k] += 1;
	    }
	}
    }

    if (freq->n > 7) {
	freq->chisqu = doornik_chisq(skew, kurt, freq->n); 
    } else {
	freq->chisqu = NADBL;
    }

    free(x);

    return freq;
}

#endif

FREQDIST *freqdist (double ***pZ, const DATAINFO *pdinfo, 
		    int varno, int params)
{
    FREQDIST *freq;
    double *x = NULL;
    double xx, xmin, xmax, xbar, sdx;
    double skew, kurt;
    double range, binwidth;
    int i, k, n;
    int nbins;

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
    if (n < 8) {
	gretl_errno = E_DATA;
	sprintf(gretl_errmsg, _("Insufficient data to build frequency "
		"distribution for variable %s"), pdinfo->varname[varno]);
	free(x);
	return freq;
    }

    freq->t1 = pdinfo->t1; 
    freq->t2 = pdinfo->t2;
    freq->n = n;

    strcpy(freq->varname, pdinfo->varname[varno]);

    if (gretl_isconst(0, n-1, x)) {
	gretl_errno = 1;
	sprintf(gretl_errmsg, _("%s is a constant"), freq->varname);
	return freq;
    }    
    
    moments(0, n-1, x, &freq->xbar, &freq->sdx, &skew, &kurt, params);
    xbar = freq->xbar;
    sdx = freq->sdx;

    gretl_minmax(0, n-1, x, &xmin, &xmax);
    range = xmax - xmin;
    
    if (n < 16) {
	nbins = 5; 
    } else if (n < 50) {
	nbins = 7;
    } else if (n > 850) {
	nbins = 29;
    } else {
	nbins = (int) sqrt((double) n);
	if (nbins % 2 == 0) nbins++;
    }

    freq->numbins = nbins;
    binwidth = range / (nbins - 1);

    freq->endpt = malloc((nbins + 1) * sizeof *freq->endpt);
    freq->midpt = malloc(nbins * sizeof *freq->midpt);
    freq->f = malloc(nbins * sizeof *freq->f);
    if (freq->endpt == NULL || freq->midpt == NULL ||
	freq->f == NULL) {
	gretl_errno = E_ALLOC;
	strcpy(gretl_errmsg, _("Out of memory for frequency distribution"));
	free(x);
	return freq;
    }
    
    freq->endpt[0] = xmin - .5 * binwidth;
    if (xmin > 0.0 && freq->endpt[0] < 0.0) {
	double rshift;

	freq->endpt[0] = 0.0;
	rshift = 1.0 - xmin / binwidth;
	freq->endpt[freq->numbins] = xmax + rshift * binwidth;
    } else {
	freq->endpt[freq->numbins] = xmax + .5 * binwidth;
    }
    
    for (k=0; k<freq->numbins; k++) {
	freq->f[k] = 0;
	if (k > 0) {
	    freq->endpt[k] = freq->endpt[k-1] + binwidth;
	}
	freq->midpt[k] = freq->endpt[k] + .5 * binwidth;
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
	for (k=1; k<freq->numbins; k++) {
	    if (freq->endpt[k] <= xx && xx < freq->endpt[k+1]) {
		freq->f[k] += 1;
		break;
	    }
	}
    }

    if (freq->n > 7) {
	freq->chisqu = doornik_chisq(skew, kurt, freq->n); 
    } else {
	freq->chisqu = NADBL;
    }

    free(x);

    return freq;
}

/* ...................................................... */

#ifdef PACF_BY_OLS

static int get_pacf (double *pacf, int m, int varnum, 
		     double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    int i, j, err = 0;
    int *laglist, *list;
    int t1 = pdinfo->t1;
    int v = pdinfo->v;
    MODEL tmp;

    pdinfo->t1 = 0; 

    list = malloc((m + 3) * sizeof *list);
    laglist = malloc(m * sizeof *laglist);
    if (list == NULL || laglist == NULL) {
	pdinfo->t1 = t1;
	return 1;
    }

    /* add appropriate number of lags to data set */
    for (i=1; i<=m; i++) {
	int lnum = laggenr(varnum, i, 0, pZ, pdinfo);

	if (lnum < 0) {
	    free(list);
	    free(laglist);
	    return 1;
	}
	laglist[i-1] = lnum; 
    }

    gretl_model_init(&tmp, pdinfo);

    pdinfo->t1 = t1;

    list[1] = varnum;
    for (i=2; i<=m; i++) {
	list[0] = i + 2;
	list[2] = 0;
	for (j=0; j<i; j++) {
	    list[j+3] = laglist[j];
	}
	tmp = lsq(list, pZ, pdinfo, OLS, OPT_A, 0);
	if ((err = tmp.errcode)) {
	    fprintf(stderr, "error estimating model for pacf\n");
	    break;
	}
	pacf[i-1] = tmp.coeff[i];
	clear_model(&tmp, pdinfo);
    }

    dataset_drop_vars(pdinfo->v - v, pZ, pdinfo);
    free(laglist);
    free(list);

    pdinfo->t1 = t1;

    return err;
}

#else

/* Durbin-Levinson algorithm */

static int get_pacf (double *pacf, const double *acf, int m)
{
    int i, j;
    gretl_matrix *phi;
    double x, num, den;

    phi = gretl_matrix_alloc(m, m);
    if (phi == NULL) return 1;

    pacf[0] = acf[0];
    gretl_matrix_set(phi, 0, 0, acf[0]);

    for (i=1; i<m; i++) {
	num = acf[i];
	for (j=0; j<i; j++) {
	    num -= gretl_matrix_get(phi, i-1, j) * acf[i-j-1];
	}
	den = 1.0;
	for (j=0; j<i; j++) {
	    den -= gretl_matrix_get(phi, i-1, j) * acf[j];
	}
	pacf[i] = num / den;
	gretl_matrix_set(phi, i, i, pacf[i]);
	for (j=0; j<i; j++) {
	    x = gretl_matrix_get(phi, i-1, j);
	    x -= pacf[i] * gretl_matrix_get(phi, i-1, i-j-1);
	    gretl_matrix_set(phi, i, j, x);
	}
    }
	
    gretl_matrix_free(phi);

    return 0;
}

#endif

static int auto_acf_order (int pd, int nobs)
{
    int acf_m;

    switch (pd) {
    case 4: 
	acf_m = (nobs <= 20)? nobs - 5 : 14; 
	break;
    case 12: 
    case 52: 
	acf_m = (nobs <= 40)? nobs - 13 : 28;
	break;
    case 24: 
	acf_m = (nobs <= 100)? nobs - 25 : 96;
	break;
    default:  
	acf_m = (nobs <= 18)? nobs - 5 : 14;
	break;
    }

    return acf_m;
}

static double gretl_acf (int n, int k, const double *y)
{
    int t;
    double z, ybar, num = 0.0, den = 0.0;

    if (n == 0 || gretl_isconst(0, n-1, y)) 
	return NADBL;

    ybar = gretl_mean(0, n-1, y);

    for (t=0; t<n; t++) {
	z = y[t] - ybar;
	den += z * z;
	if (t >= k) {
	    num += z * (y[t-k] - ybar);
	}
    }

    return num / den;
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
    double *acf, box, pm;
    double *pacf = NULL;
    int k, l, acf_m, pacf_m, nobs; 
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int list[2];
    FILE *fq = NULL;
    int err = 0, pacf_err = 0;

    list[0] = 1;
    list[1] = varno;
    adjust_t1t2(NULL, list, &t1, &t2, (const double **) *pZ, NULL);
    nobs = t2 - t1 + 1;

    if (missvals(&(*pZ)[varno][t1], nobs)) {
	pprintf(prn, "\n%s",
		_("Missing values within sample -- can't do correlogram"));
	return 1;
    }

    if (nobs < 4) {
	pputs(prn, _("\nInsufficient observations for correlogram"));
	return 1;
    }

    if (gretl_isconst(t1, t2, &(*pZ)[varno][0])) {
	sprintf(gretl_tmp_str, _("%s is a constant"), 
		pdinfo->varname[varno]);
	pprintf(prn, "\n%s\n", gretl_tmp_str);
	return 1;
    }

    acf_m = order;
    if (acf_m == 0) {
	acf_m = auto_acf_order(pdinfo->pd, nobs);
    } else if (acf_m > nobs - pdinfo->pd) {
	acf_m = nobs - 1;
    }

    acf = malloc(acf_m * sizeof *acf);
    if (acf == NULL) {
	return E_ALLOC;    
    }

    /* calculate acf, with lag order m */
    for (l=1; l<=acf_m; l++) {
	acf[l-1] = gretl_acf(nobs, l, &(*pZ)[varno][t1]);
    }

    sprintf(gretl_tmp_str, _("Autocorrelation function for %s"), 
	    pdinfo->varname[varno]);
    pprintf(prn, "\n%s\n\n", gretl_tmp_str);

    /* add Ljung-Box statistic */
    box = 0;
    for (t=0; t<acf_m; t++) { 
	box += acf[t] * acf[t] / (nobs - t + 1);
    }
    box *= nobs * (nobs + 2.0);

    pprintf(prn, "Ljung-Box Q' = %.4f\n", box);
    pprintf(prn, _("Degrees of freedom = %d, p-value = %.4f\n\n"),
	    acf_m, chisq(box, acf_m));

    /* print acf */
    for (t=0; t<acf_m; t++) {
	pprintf(prn, "%5d)%8.4f", t + 1, acf[t]);
	if ((t + 1) % 5 == 0) pputc(prn, '\n');
    }
    pputc(prn, '\n');

    if (batch) { 
	/* batch mode: use ASCII graphics, not gnuplot */
	double *xl = malloc(acf_m * sizeof *xl);

	if (xl == NULL) {
	    err = E_ALLOC;
	    goto acf_getout;
	}
	for (l=0; l<acf_m; l++) {
	    xl[l] = l + 1.0;
	}
        pprintf(prn, "\n\n%s\n\n", _("Correlogram"));
	graphyzx(NULL, acf, NULL, xl, acf_m, pdinfo->varname[varno], 
		 _("lag"), NULL, 0, prn);
	free(xl);
    } 

    /* determine lag order for pacf (may have to be shorter than acf_m) */
    if (acf_m > nobs / 2 - 1) pacf_m = nobs / 2 - 1;
    else pacf_m = acf_m;

    /* generate (and if not in batch mode) plot partial 
       autocorrelation function */

    pacf = malloc(pacf_m * sizeof *pacf); 
    if (pacf == NULL) {
	err = E_ALLOC;
	goto acf_getout;
    }

#ifdef PACF_BY_OLS
    err = pacf_err = get_pacf(pacf, pacf_m, varno, pZ, pdinfo, prn);
#else
    err = pacf_err = get_pacf(pacf, acf, pacf_m);
#endif

    if (!err) {
	pacf[0] = acf[0];
	pprintf(prn, "\n%s", _("Partial autocorrelations"));
	if (pacf_m < acf_m) {
	    pprintf(prn, " (%s %d):\n\n", _("to lag"), pacf_m);
	} else {
	    pputs(prn, ":\n\n");
	}
	for (k=0; k<pacf_m; k++) {
	    pprintf(prn, "%5d)%8.4f", k+1, pacf[k]);
	    if ((k + 1) % 5 == 0) pputc(prn, '\n');
	}
    }
    pputc(prn, '\n');
    if (pacf_m % 5 > 0) pputc(prn, '\n');

    if (batch) {
	goto acf_getout;
    } else if (gnuplot_init(ppaths, PLOT_CORRELOGRAM, &fq)) {
	err = E_FOPEN;
	goto acf_getout;
    }

    /* for confidence bands */
    pm = 1.0 / sqrt((double) nobs);

    fprintf(fq, "# correlogram\n");

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    /* create two separate plots, if both are OK */
    if (!pacf_err) {
	fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.48\n", fq);
    }
    fputs("set xzeroaxis\n", fq);
    fputs("set key top right\n", fq); 
    fprintf(fq, "set xlabel \"%s\"\n", _("lag"));
    fputs("set yrange [-1.1:1.1]\n", fq);

    /* upper plot: Autocorrelation Function or ACF */
    if (!pacf_err) {
	fputs("set origin 0.0,0.50\n", fq);
    }
    fprintf(fq, "set title \"%s %s\"\n", I_("ACF for"), 
	    pdinfo->varname[varno]);
    fprintf(fq, "set xrange [0:%d]\n", acf_m + 1);
    fprintf(fq, "plot \\\n"
	    "'-' using 1:2 notitle w impulses, \\\n"
	    "%g title '%s' lt 2, \\\n"
	    "%g notitle lt 2\n", pm, 
	    "+- 1.96/T^0.5",
	    -pm);
    for (k=0; k<acf_m; k++) {
	fprintf(fq, "%d %g\n", k + 1, acf[k]);
    }
    fputs("e\n", fq);

    if (!pacf_err) {
	/* lower plot: Partial Autocorrelation Function or PACF */
	fputs("set origin 0.0,0.0\n", fq);
	fprintf(fq, "set title \"%s %s\"\n", I_("PACF for"), 
		pdinfo->varname[varno]);
	fprintf(fq, "set xrange [0:%d]\n", pacf_m + 1);
	fprintf(fq, "plot \\\n"
		"'-' using 1:2 notitle w impulses, \\\n"
		"%g title '%s' lt 2, \\\n"
		"%g notitle lt 2\n", pm,
		"+- 1.96/T^0.5",
		-pm);
	for (k=0; k<pacf_m; k++) {
	    fprintf(fq, "%d %g\n", k + 1, pacf[k]);
	}
	fputs("e\n", fq);
    }

    if (!pacf_err) fputs("set nomultiplot\n", fq);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fq);

    err = gnuplot_display(ppaths);

 acf_getout:
    free(acf);
    free(pacf);

    return err;
}

/* ...................................................... */

static int roundup_half (int i)
{
    return (int) ceil((double) i / 2.0);
}

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
    if (start_new_Z(&tmpZ, &tmpdinfo, 1))
	return 1;

    /* Test from Geweke and Porter-Hudak, as set out in
       Greene, Econometric Analysis 4e, p. 787 */
    
    for (t=0; t<n; t++) {
	tmpZ[0][t] = 1.0;
	tmpZ[1][t] = log(hhat[t]);
	xx = sin(omega[t] / 2);
	tmpZ[2][t] = log(4 * xx * xx);
    }

    list[0] = 3;
    list[1] = 1;
    list[2] = 0;
    list[3] = 2;

    gretl_model_init(&tmp, &tmpdinfo);
    tmp = lsq(list, &tmpZ, &tmpdinfo, OLS, OPT_A, 0);

    if (!tmp.errcode) {
	tstat = -tmp.coeff[1] / tmp.sderr[1];
	pprintf(prn, "\n%s\n"
		"  %s = %g\n"
		"  %s: t(%d) = %g, %s %.4f\n",
		_("Test for fractional integration"),
		_("Estimated degree of integration"), -tmp.coeff[1], 
		_("test statistic"), tmp.dfd, tstat, 
		_("with p-value"), tprob(tstat, tmp.dfd));
    } else {
	err = tmp.errcode;
    }

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
    int do_graph = !batch;
    FILE *fq = NULL;

    *gretl_errmsg = 0;

    list[0] = 1;
    list[1] = varno;
    adjust_t1t2(NULL, list, &t1, &t2, (const double **) *pZ, NULL);
    nobs = t2 - t1 + 1;

    if (missvals(&(*pZ)[varno][t1], nobs)) {
	strcpy(gretl_errmsg, 
	       _("Missing values within sample -- can't do periodogram"));
	return 1;
    }    

    if (nobs < 12) {
	strcpy(gretl_errmsg,
	       _("Insufficient observations for periodogram"));
	return 1;
    }

    if (gretl_isconst(t1, t2, &(*pZ)[varno][0])) {
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

    xx = gretl_mean(t1, t2, (*pZ)[varno]);

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

    if (do_graph && gnuplot_init(ppaths, PLOT_PERIODOGRAM, &fq) == 0) {
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
	else {
	    fprintf(fq, "'\n");
	}
	fprintf(fq, "set xrange [0:%d]\n", xmax);
	fprintf(fq, "plot '-' using 1:2 w lines\n");
    }

    if (do_graph && fq == NULL) {
	do_graph = 0;
	err = 1;
    }

    pprintf(prn, _("\nPeriodogram for %s\n"), pdinfo->varname[varno]);
    pprintf(prn, _("Number of observations = %d\n"), nobs);
    if (opt) 
	pprintf(prn, _("Using Bartlett lag window, length %d\n\n"), L);
    pputs(prn, _(" omega  scaled frequency  periods  spectral density\n\n"));

    if (do_graph) { 
	savexx = malloc((1 + nobs/2) * sizeof *savexx);
	if (savexx == NULL) {
	    err = 1;
	    fclose(fq);
	    do_graph = 0;
	}
    }

    varx = gretl_variance(t1, t2, &(*pZ)[varno][0]);
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
	pprintf(prn, " %.4f%9d%16.2f%16.5f\n", yy, t, 
		(double) (nobs / 2) / (2 * t), xx);
	if (savexx != NULL) savexx[t] = xx;
	if (t <= nT) {
	    omega[t-1] = yy;
	    hhat[t-1] = xx;
	}
    }
    pputc(prn, '\n');

    if (do_graph) {
#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "C");
#endif
	for (t=1; t<=nobs/2; t++) fprintf(fq, "%d %f\n", t, savexx[t]);
#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "");
#endif
	fprintf(fq, "e\n");

	fclose(fq);
	free(savexx);
	err = gnuplot_display(ppaths);
    }

    if (opt == 0 && fract_int(nT, hhat, omega, prn)) {
	pprintf(prn, "\n%s\n",
		_("Fractional integration test failed"));
    }

    free(autocov);
    free(omega);
    free(hhat);

    return err;
}

/* ............................................................. */

static void printf15 (double zz, PRN *prn)
{
    if (na(zz)) pprintf(prn, "%*s", UTF_WIDTH(_("undefined"), 15), 
			_("undefined"));
    else {
	pputc(prn, ' ');
	gretl_print_fullwidth_double(zz, 5, prn);	
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
    char date1[OBSLEN], date2[OBSLEN], tmp[96];

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    pputc(prn, '\n');

    sprintf(tmp, _("%s, using the observations %s - %s"), str, date1, date2);
    center_line(tmp, prn, 0);

    if (ci == CORR) {
	strcpy(tmp, _("(missing values denoted by -999 will be skipped)"));
	center_line(tmp, prn, 1);
    }
}

static void print_summary_single (GRETLSUMMARY *summ,
				  const DATAINFO *pdinfo,
				  PRN *prn)
{
    char obs1[OBSLEN], obs2[OBSLEN], tmp[128];
    double vals[8];
    double xbar, sd, cv = NADBL;
    const char *labels[] = {
	N_("Mean"),
	N_("Median"),
	N_("Minimum"),
	N_("Maximum"),
	N_("Standard deviation"),
	N_("C.V."),
	N_("Skewness"),
	N_("Ex. kurtosis")
    };
    int slen = 0, i = 0;

    ntodate(obs1, pdinfo->t1, pdinfo);
    ntodate(obs2, pdinfo->t2, pdinfo);

    prhdr(_("Summary Statistics"), pdinfo, SUMMARY, prn);
    sprintf(tmp, _("for the variable '%s' (%d valid observations)"), 
	    pdinfo->varname[summ->list[1]], summ->n);
    center_line(tmp, prn, 1);

    xbar = summ->coeff[0];
    sd = summ->sderr[0];
    if (xbar != 0.0 && !na(sd)) cv = fabs(sd/xbar);

    vals[0] = xbar;
    vals[1] = summ->xmedian[0];
    vals[2] = summ->xpx[0];
    vals[3] = summ->xpy[0];
    vals[4] = sd;
    vals[5] = cv;
    vals[6] = summ->xskew[0];
    vals[7] = summ->xkurt[0];

    for (i=0; i<8; i++) {
	if (strlen(_(labels[i])) > slen) {
#if defined(ENABLE_NLS) && defined(USE_GTK2)
	    slen = g_utf8_strlen(_(labels[i]), -1);	    
#else
	    slen = strlen(_(labels[i]));
#endif
	}
    }
    slen++;

    for (i=0; i<8; i++) {
	pprintf(prn, "  %-*s", UTF_WIDTH(_(labels[i]), slen), _(labels[i]));
	printf15(vals[i], prn);
	pputc(prn, '\n');
    }

    pputs(prn, "\n\n");    
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
    char tmp[128];

    if (lo == 1) {
	print_summary_single(summ, pdinfo, prn);
	return;
    }

    prhdr(_("Summary Statistics"), pdinfo, SUMMARY, prn);
    strcpy(tmp, _("(missing values denoted by -999 will be skipped)"));
    center_line(tmp, prn, 1);
    pprintf(prn, "\n%s  ", _("Variable"));

    pputs(prn, _("      MEAN           MEDIAN           MIN"
            "             MAX\n\n"));

    for (v=0; v<lo; v++) {
	if (pause) page_break(1, &lineno, 0);
	lineno++;
	lv = summ->list[v+1];
	pprintf(prn, "%-10s", pdinfo->varname[lv]);
	xbar = summ->coeff[v];
	printf15(xbar, prn);
	printf15(summ->xmedian[v], prn);
	printf15(summ->xpx[v], prn);
	printf15(summ->xpy[v], prn);
	pputc(prn, '\n');
    }

    if (pause) page_break(lo + 2, &lineno, 0);
    lineno += 2;
    pputc(prn, '\n');

    pprintf(prn, "\n%s  ", _("Variable"));
    pputs(prn, _("      S.D.            C.V.           "
	 " SKEW          EXCSKURT\n\n"));

    for (v=0; v<lo; v++) {
	if (pause) page_break(1, &lineno, 0);
	lineno++;
	lv = summ->list[v+1];
	pprintf(prn, "%-10s", pdinfo->varname[lv]);

	xbar = summ->coeff[v];
	std = summ->sderr[v];
	if (xbar != 0.0) xcv = (xbar > 0)? std/xbar: (-1) * std/xbar;
	else xcv = NADBL;

	printf15(std, prn);
	printf15(xcv, prn);
	printf15(summ->xskew[v], prn);
	printf15(summ->xkurt[v], prn);
	pputc(prn, '\n');
    }
    pputc(prn, '\n');
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
    int lo = list[0];
    int v, *tmp = NULL;
    GRETLSUMMARY *summ;
    double xbar, std, low, high, skew, kurt, *x = NULL;

    summ = malloc(sizeof *summ);
    if (summ == NULL) return NULL;
    summ->list = NULL;

    if ((summ->xskew = malloc(lo * sizeof(double))) == NULL) return NULL;
    if ((summ->xkurt = malloc(lo * sizeof(double))) == NULL) return NULL;
    if ((summ->xmedian = malloc(lo * sizeof(double))) == NULL) return NULL;
    if ((summ->coeff = malloc(lo * sizeof(double))) == NULL) return NULL;
    if ((summ->sderr = malloc(lo * sizeof(double))) == NULL) return NULL;
    if ((summ->xpx = malloc(lo * sizeof(double))) == NULL) return NULL;
    if ((summ->xpy = malloc(lo * sizeof(double))) == NULL) return NULL;
    if ((x = malloc((pdinfo->t2 - pdinfo->t1 + 1) * sizeof *x)) == NULL) 
	return NULL;

    for (v=0; v<lo; v++)  {
	summ->n = ztox(list[v+1], x, *pZ, pdinfo);
	if (summ->n < 2) { /* zero or one observations */
	    if (summ->n == 0)
		pprintf(prn, _("Dropping %s: sample range contains no valid "
			"observations\n"), pdinfo->varname[list[v+1]]);
	    else
		pprintf(prn, _("Dropping %s: sample range has only one "
			"obs, namely %g\n"), pdinfo->varname[list[v+1]], x[0]);
	    list_exclude(v+1, list);
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

	gretl_minmax(0, summ->n-1, x, &low, &high);	
	moments(0, summ->n-1, x, &xbar, &std, &skew, &kurt, 1);

	summ->xpx[v] = low;
	summ->xpy[v] = high;
	summ->coeff[v] = xbar;
	summ->sderr[v] = std;
	summ->xskew[v] = skew;
	summ->xkurt[v] = kurt;

	if (summ->n > 1) {
	    summ->xmedian[v] = gretl_median(x, summ->n);
	} else {
	    summ->xmedian[v] = x[1];
	}
    } /* end loop over variables in list */

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
	if (gretl_isconst(t1, t2, (*pZ)[p[i]])) {
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
	    corrmat->xpx[nij] = 
		gretl_corr(corrmat->n, 
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
    text_print_matrix(corr->xpx, corr->list, NULL, pdinfo, pause, prn);
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
 * @vardiff: if non-zero, assume population variances are different.
 * @prn: gretl printing struct.
 *
 * Carries out test of the null hypothesis that the means of two
 * variables are equal.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int means_test (LIST list, double **Z, const DATAINFO *pdinfo, 
		unsigned long vardiff, PRN *prn)
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
	pputs(prn, _("Sample range has no valid observations."));
	free(x); free(y);
	return 1;
    }
    if (n1 == 1 || n2 == 1) {
	pputs(prn, _("Sample range has only one observation."));
	free(x); free(y);
	return 1;
    }
    df = n1 + n2 - 2;

    moments(0, n1-1, x, &m1, &s1, &skew, &kurt, 1);
    moments(0, n2-1, y, &m2, &s2, &skew, &kurt, 1);
    mdiff = m1 - m2;

    if (vardiff) {
	se = sqrt((s1 * s1 / n1) + (s2 * s2 / n2));
    } else {
	double sp2;

	sp2 = ((n1-1) * s1 * s1 + (n2-1) * s2 * s2) / df;
	se = sqrt(sp2 / n1 + sp2 / n2);
    }

    t = mdiff / se;
    pval = tprob(t, df);

    pprintf(prn, _("\nEquality of means test "
	    "(assuming %s variances)\n\n"), (vardiff)? _("unequal") : _("equal"));
    pprintf(prn, "   %s: ", pdinfo->varname[list[1]]);
    pprintf(prn, _("Number of observations = %d\n"), n1);
    pprintf(prn, "   %s: ", pdinfo->varname[list[2]]);
    pprintf(prn, _("Number of observations = %d\n"), n2);
    pprintf(prn, _("   Difference between sample means = %g - %g = %g\n"), 
	    m1, m2, mdiff);
    pputs(prn, _("   Null hypothesis: The two population means are the same.\n"));
    pprintf(prn, _("   Estimated standard error = %g\n"), se);
    pprintf(prn, _("   Test statistic: t(%d) = %g\n"), df, t);
    pprintf(prn, _("   p-value (two-tailed) = %g\n\n"), pval);
    if (pval > .10)
	pputs(prn, _("   The difference is not statistically significant.\n\n"));

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
	pputs(prn, _("Sample range has no valid observations."));
	free(x); free(y);
	return 1;
    }
    if (n1 == 1 || n2 == 1) {
	pputs(prn, _("Sample range has only one observation."));
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

    pputs(prn, _("\nEquality of variances test\n\n"));
    pprintf(prn, _("   Ratio of sample variances = %g\n"), F);
    pprintf(prn, "   %s: %s\n", _("Null hypothesis"), 
	    _("The two population variances are equal"));
    pprintf(prn, "   %s: F(%d,%d) = %g\n", _("Test statistic"), dfn, dfd, F);
    pprintf(prn, _("   p-value (two-tailed) = %g\n\n"), fdist(F, dfn, dfd));
    if (fdist(F, dfn, dfd) > .10)
	pputs(prn, _("   The difference is not statistically significant.\n\n"));

    free(x);
    free(y);

    return 0;
}
