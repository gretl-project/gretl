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
 *
 */

#include "libgretl.h"
#include "gretl_matrix.h"
#include "gretl_cmatrix.h"
#include "matrix_extra.h"
#include "gretl_panel.h"
#include "gretl_string_table.h"
#include "libset.h"
#include "plotspec.h"
#include "usermat.h"

#include <unistd.h>

#ifdef WIN32
# include <windows.h>
#endif

/**
 * SECTION:describe
 * @short_description: descriptive statistics plus some tests
 * @title: Descriptive statistics
 * @include: libgretl.h
 *
 * Computation and printing of numerous descriptive statistics
 * along with some hypothesis tests, for example regarding the
 * normality of a data series.
 */

/* return the number of "good" (non-missing) observations in series x;
   also (optionally) write to x0 the first non-missing value */

static int good_obs (const double *x, int n, double *x0)
{
    int t, good = n;

    if (x0 != NULL) {
	*x0 = NADBL;
    }

    for (t=0; t<n; t++) {
	if (na(x[t])) {
	    good--;
	} else if (x0 != NULL && na(*x0)) {
	    *x0 = x[t];
	}
    }

    return good;
}

/**
 * gretl_minmax:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @min: location to receive minimum value.
 * @max: location to receive maximum value.
 *
 * Puts the minimum and maximum values of the series @x,
 * from obs @t1 to obs @t2, into the variables @min and @max.
 *
 * Returns: the number of valid observations in the given
 * data range.
 */

int gretl_minmax (int t1, int t2, const double *x,
		  double *min, double *max)
{
    int t, n = 0;

    *max = *min = NADBL;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    if (n == 0) {
		*max = *min = x[t];
	    } else {
		if (x[t] > *max) *max = x[t];
		if (x[t] < *min) *min = x[t];
	    }
	    n++;
	}
    }

    return n;
}

static int summary_minmax (int t1, int t2, const double *x,
			   Summary *summ, int i)
{
    int n = gretl_minmax(t1, t2, x, &summ->low[i], &summ->high[i]);

    if (summ->high[i] > 1.0e5 && summ->high[i] < 1.0e10) {
	/* 2024-05-29: add check for big integers */
	summ->ival[i] = gretl_isint(t1, t2, x);
    }

    return n;
}

/**
 * gretl_min:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the minimum value of @x over the given range,
 * or #NADBL if no valid vaues are found.
 */

double gretl_min (int t1, int t2, const double *x)
{
    double min, max;

    gretl_minmax(t1, t2, x, &min, &max);
    return min;
}

/**
 * gretl_max:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the maximum value of @x over the given range,
 * or #NADBL if no valid vaues are found.
 */

double gretl_max (int t1, int t2, const double *x)
{
    double min, max;

    gretl_minmax(t1, t2, x, &min, &max);
    return max;
}

/**
 * gretl_sum:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the sum of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * in case there are no valid observations.
 */

double gretl_sum (int t1, int t2, const double *x)
{
    double sum = 0.0;
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
	if (!(na(x[t]))) {
	    sum += x[t];
	    n++;
	}
    }

    return (n == 0)? NADBL : sum;
}

/**
 * gretl_mean:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the arithmetic mean of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * in case there are no valid observations.
 */

double gretl_mean (int t1, int t2, const double *x)
{
    double xbar, sum = 0.0;
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
	if (!(na(x[t]))) {
	    sum += x[t];
	    n++;
	}
    }

    if (n == 0) {
	return NADBL;
    }

    xbar = sum / n;
    sum = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!(na(x[t]))) {
	    sum += (x[t] - xbar);
	}
    }

    return xbar + sum / n;
}

/**
 * eval_ytest:
 * @y: reference numerical value.
 * @op: operator.
 * @test: numerical test value.
 *
 * Returns: 1 if the expression @y @yop @test (for example
 * "y = 2" or "y <= 45") evaluates as true, else 0.
 */

int eval_ytest (double y, GretlOp op, double test)
{
    int ret = 0;

    switch (op) {
    case OP_EQ:
	if (na(test)) {
	    ret = na(y);
	} else {
	    ret = (y == test);
	}
	break;
    case OP_GT:
	ret = (y > test);
	break;
    case OP_LT:
	ret = (y < test);
	break;
    case OP_NEQ:
	if (na(test)) {
	    ret = !na(y);
	} else {
	    ret = (y != test);
	}
	break;
    case OP_GTE:
	ret = (y >= test);
	break;
    case OP_LTE:
	ret = (y <= test);
	break;
    default:
	break;
    }

    return ret;
}

/**
 * gretl_restricted_mean:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: criterion series.
 * @yop: criterion operator.
 * @yval: criterion value.
 *
 * Returns: the arithmetic mean of the series @x in the
 * range @t1 to @t2 (inclusive), but including only
 * observations where the criterion variable @y bears the
 * relationship @yop to the value @yval -- or #NADBL in case
 * there are no observations that satisfy the restriction.
 */

double gretl_restricted_mean (int t1, int t2, const double *x,
			      const double *y, GretlOp yop,
			      double yval)
{
    double xbar, sum = 0.0;
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t]) && eval_ytest(y[t], yop, yval)) {
	    sum += x[t];
	    n++;
	}
    }

    if (n == 0) {
	return NADBL;
    }

    xbar = sum / n;
    sum = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t]) && eval_ytest(y[t], yop, yval)) {
	    sum += (x[t] - xbar);
	}
    }

    return xbar + sum / n;
}

/*
   For an array a with n elements, find the element which would be a[h]
   if the array were sorted from smallest to largest (without the need
   to do a full sort).  The algorithm is due to C. A. R. Hoare; the C
   implementation is shamelessly copied, with minor adaptations, from
   http://www.astro.northwestern.edu/~ato/um/programming/sorting.html
   (Atakan G\"urkan: link now dead).
*/

static double find_hoare (double *a, int n, int h)
{
    double w, x;
    int l = 0;
    int r = n - 1;
    int i, j;

    while (l < r) {
	x = a[h];
	i = l;
	j = r;
	while (j >= i) {
	    while (a[i] < x) i++;
	    while (x < a[j]) j--;
	    if (i <= j) {
		w = a[i];
		a[i++] = a[j];
		a[j--] = w;
	    }
	}
	if (j < h) l = i;
	if (h < i) r = j;
    }

    return a[h];
}

static double find_hoare_inexact (double *a, double p,
				  double xmin, double xmax,
				  double frac, int n,
				  int hf, int hc)
{
    double tmp, high, low;
    int i;

    if (p < 0.5) {
	high = find_hoare(a, n, hc);
	tmp = xmin;
	for (i=0; i<hc; i++) {
	    if (a[i] > tmp) {
		tmp = a[i];
	    }
	}
	low = tmp;
    } else {
	low = find_hoare(a, n, hf);
	tmp = xmax;
	for (i=hc; i<n; i++) {
	    if (a[i] < tmp) {
		tmp = a[i];
	    }
	}
	high = tmp;
    }

    return low + frac * (high - low);
}

/* See https://en.wikipedia.org/wiki/Quantile, also
   Hyndman and Fan, "Sample Quantiles in Statistical
   Packages" (The American Statistician, 1996).
   Return the "index" of the @p-quantile for a sample
   of size @n (which may not be an integer).
*/

static double quantile_index (int n, double p)
{
    int Qtype = libset_get_int(QUANTILE_TYPE);
    double h;

    if (Qtype == 0) {
	/* \hat{Q}_6 */
	h = (n + 1) * p;
    } else if (Qtype == 1) {
	/* \hat{Q}_7 */
	h = (n - 1) * p + 1;
    } else {
	/* \hat{Q}_8 */
	h = (n + 1.0/3) * p + 1.0/3;
    }

    return h - 1; /* 0-based */
}

/**
 * gretl_quantile:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @p: probability.
 * @opt: may include OPT_Q to hush warning when sample is
 * too small.
 * @err: location to receive error code.
 *
 * Returns: the @p quantile of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_quantile (int t1, int t2, const double *x, double p,
		       gretlopt opt, int *err)
{
    double *a = NULL;
    double xmin, xmax;
    double h, ret;
    int hf, hc;
    int t, n = 0;

    if (*err) {
	/* don't compound a prior error */
	return NADBL;
    }

    if (p <= 0 || p >= 1) {
	/* sanity check */
	*err = E_DATA;
    } else {
	n = gretl_minmax(t1, t2, x, &xmin, &xmax);
	if (n == 0) {
	    *err = E_DATA;
	}
    }

    if (!*err) {
	h = quantile_index(n, p);
	hf = floor(h);
	hc = ceil(h);

	if (hc == 0 || hc == n) {
	    /* too few usable observations for the specified
	       quantile; don't treat as fatal error */
	    return NADBL;
	} else {
	    a = malloc(n * sizeof *a);
	    if (a == NULL) {
		*err = E_ALLOC;
	    }
	}
    }

    if (*err) {
	return NADBL;
    }

    n = 0;
    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    a[n++] = x[t];
	}
    }

    if (hf == hc) {
	/* "exact" */
	ret = find_hoare(a, n, hf);
    } else {
	ret = find_hoare_inexact(a, p, xmin, xmax, h - hf,
				 n, hf, hc);
    }

    free(a);

    return ret;
}

/**
 * gretl_array_quantiles:
 * @a: data array (this gets re-ordered).
 * @n: length of array.
 * @p: array of probabilities (over-written by quantiles).
 * @k: number of probabilities.
 *
 * Computes @k quantiles (given by the elements of @p) for the
 * first n elements of the array @a, which is re-ordered in
 * the process. On successful exit, @p contains the quantiles.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_array_quantiles (double *a, int n, double *p, int k)
{
    double h, xmin = 0, xmax = NADBL;
    int hf, hc, i;
    int err = 0;

    if (n <= 0 || k <= 0) {
	return E_DATA;
    }

    for (i=0; i<k; i++) {
	if (p[i] <= 0.0 || p[i] >= 1.0) {
	    p[i] = NADBL;
	    err = E_INVARG;
	    break;
	}

	h = quantile_index(n, p[i]);
	hf = floor(h);
	hc = ceil(h);

	if (hc == 0 || hc == n) {
	    p[i] = NADBL;
	} else if (hf == hc) {
	    p[i] = find_hoare(a, n, hf);
	} else {
	    if (na(xmax)) {
		gretl_minmax(0, n-1, a, &xmin, &xmax);
	    }
	    p[i] = find_hoare_inexact(a, p[i], xmin, xmax, h - hf,
				      n, hf, hc);
	}
    }

    return err;
}

/**
 * gretl_array_quantile:
 * @a: array on which to operate.
 * @n: number of elements in @a.
 * @p: probability.
 *
 * Returns: the @p quantile of the first @n elements in @a,
 * which is re-ordered in the process, or #NADBL on
 * failure.
 */

double gretl_array_quantile (double *a, int n, double p)
{
    gretl_array_quantiles(a, n, &p, 1);
    return p;
}

/**
 * gretl_median:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the median value of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_median (int t1, int t2, const double *x)
{
    int err = 0;

    return gretl_quantile(t1, t2, x, 0.5, OPT_NONE, &err);
}

/**
 * gretl_sst:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the sum of squared deviations from the mean for
 * the series @x from obs @t1 to obs @t2, skipping any missing
 * values, or #NADBL on failure.
 */

double gretl_sst (int t1, int t2, const double *x)
{
    int t, n = t2 - t1 + 1;
    double sumsq, xx, xbar;

    if (n == 0) {
	return NADBL;
    }

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
	return NADBL;
    }

    sumsq = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    xx = x[t] - xbar;
	    sumsq += xx * xx;
	}
    }

    return sumsq;
}

/* allows for computing MLE rather than using df correction,
   via the @asy argument */

double real_gretl_variance (int t1, int t2, const double *x,
			    int asy)
{
    int t, n = t2 - t1 + 1;
    double v, dx, xbar;

    if (n == 0) {
	/* null sample */
	return NADBL;
    }

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
	return NADBL;
    }

    v = 0.0;
    n = 0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    dx = x[t] - xbar;
	    v += dx * dx;
	    n++;
	}
    }

    if (asy && n > 0) {
	v /= n;
    } else if (n > 1) {
	v /= (n - 1);
    } else {
	v = 0.0;
    }

    return (v >= 0)? v : NADBL;
}

/**
 * gretl_variance:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the variance of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_variance (int t1, int t2, const double *x)
{
    return real_gretl_variance(t1, t2, x, 0);
}

/**
 * gretl_restricted_variance:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: criterion series.
 * @yop: criterion operator.
 * @yval: criterion value.
 *
 * Returns: the variance of the series @x from obs
 * @t1 to obs @t2, skipping any missing values and
 * observations where the series @y does not bear the
 * relationship @yop to the value @yval, or #NADBL on
 * failure.
 */

double gretl_restricted_variance (int t1, int t2, const double *x,
				  const double *y, GretlOp yop,
				  double yval)
{
    double sumsq, xx, xbar;
    int t, n = 0;

    xbar = gretl_restricted_mean(t1, t2, x, y, yop, yval);
    if (na(xbar)) {
	return NADBL;
    }

    sumsq = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t]) && eval_ytest(y[t], yop, yval)) {
	    xx = x[t] - xbar;
	    sumsq += xx * xx;
	    n++;
	}
    }

    sumsq = (n > 1)? sumsq / (n - 1) : 0.0;

    return (sumsq >= 0)? sumsq : NADBL;
}

/**
 * gretl_stddev:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the standard deviation of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_stddev (int t1, int t2, const double *x)
{
    double v = gretl_variance(t1, t2, x);

    return (na(v))? v : sqrt(v);
}

/**
 * gretl_restricted_stddev:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: criterion series.
 * @yop: criterion operator.
 * @yval: criterion value.
 *
 * Returns: the standard deviation of the series @x from obs
 * @t1 to obs @t2, skipping any missing values and observations
 * where the series @y does not bear the relationship @yop to
 * the value @yval, or #NADBL on failure.
 */

double gretl_restricted_stddev (int t1, int t2, const double *x,
				const double *y, GretlOp yop,
				double yval)
{
    double xx = gretl_restricted_variance(t1, t2, x, y, yop, yval);

    return (na(xx))? xx : sqrt(xx);
}

/**
* gretl_long_run_variance:
* @t1: starting observation.
* @t2: ending observation.
* @x: data series.
* @m: bandwidth (<= 0 for automatic).
* @mu: mean (or NADBL to use sample mean).
*
* Returns: the long-run variance of the series @x from obs
* @t1 to obs @t2, using Bartlett kernel weights, or #NADBL
* on failure (which includes encountering missing values).
*/

double gretl_long_run_variance (int t1, int t2, const double *x,
				int m, double mu)
{
    double zt, wi, xbar, s2 = 0.0;
    int use_mu = !na(mu);
    int i, t, n, order;

    if (series_adjust_sample(x, &t1, &t2) != 0) {
	return NADBL;
    }

    n = t2 - t1 + 1;

    if (n < 2) {
	return NADBL;
    }

    if (m < 0) {
	order = (int) exp(log(n) / 3.0);
    } else {
	order = m;
    }

    if (use_mu) {
	xbar = mu;
    }

    if (use_mu && mu == 0.0) {
	for (t=t1; t<=t2; t++) {
	    s2 += x[t] * x[t];
	}
    } else {
	if (!use_mu) {
	    xbar = 0.0;
	    for (t=t1; t<=t2; t++) {
		xbar += x[t];
	    }
	    xbar /= n;
	}
	for (t=t1; t<=t2; t++) {
	    zt = x[t] - xbar;
	    s2 += zt * zt;
	}
    }

    for (i=1; i<=order; i++) {
	wi = 1.0 - i /((double) order + 1);
	for (t=t1+i; t<=t2; t++) {
	    zt = x[t] - xbar;
	    s2 += 2 * wi * zt * (x[t-i] - xbar);
	}
    }

    return s2 / n;
}

/**
 * gretl_covar:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: data series.
 * @missing: location to receive information on the number
 * of missing observations that were skipped, or %NULL.
 *
 * Returns: the covariance of the series @x and @y from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_covar (int t1, int t2, const double *x, const double *y,
		    int *missing)
{
    double sx, sy, sxy, xbar, ybar;
    int t, nn, n = t2 - t1 + 1;

    if (n == 0) {
	/* void sample */
	return NADBL;
    }

    nn = 0;
    sx = sy = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && !na(y[t])) {
	    sx += x[t];
	    sy += y[t];
	    nn++;
	}
    }

    if (nn < 2) {
	return NADBL;
    }

    xbar = sx / nn;
    ybar = sy / nn;
    sxy = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && !na(y[t])) {
	    sx = x[t] - xbar;
	    sy = y[t] - ybar;
	    sxy = sxy + (sx * sy);
	}
    }

    if (missing != NULL) {
	/* record number of missing obs */
	*missing = n - nn;
    }

    return sxy / (nn - 1);
}

/**
 * gretl_corr:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: data series.
 * @missing: location to receive information on the number
 * of missing observations that were skipped, or %NULL.
 *
 * Returns: the correlation coefficient for the series @x and @y
 * from obs @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_corr (int t1, int t2, const double *x, const double *y,
		   int *missing)
{
    int t, nn, n = t2 - t1 + 1;
    double sx, sy, sxx, syy, sxy, den, xbar, ybar;
    double cval = 0.0;

    if (n == 0) {
	/* void sample */
	return NADBL;
    }

    if (gretl_isconst(t1, t2, x) || gretl_isconst(t1, t2, y)) {
	/* correlation between any x and a constant is
	   undefined */
	return NADBL;
    }

    nn = 0;
    sx = sy = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && !na(y[t])) {
  	    sx += x[t];
	    sy += y[t];
	    nn++;
	}
    }

    if (nn < 2) {
	/* zero or 1 observations: no go! */
	return NADBL;
    }

    xbar = sx / nn;
    ybar = sy / nn;
    sxx = syy = sxy = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && !na(y[t])) {
	    sx = x[t] - xbar;
	    sy = y[t] - ybar;
	    sxx += sx * sx;
	    syy += sy * sy;
	    sxy += sx * sy;
	}
    }

    if (sxy != 0.0) {
        den = sxx * syy;
        if (den > 0.0) {
	    cval = sxy / sqrt(den);
        } else {
	    cval = NADBL;
	}
    }

    if (missing != NULL) {
	/* record number of missing observations */
	*missing = n - nn;
    }

    return cval;
}

/**
 * gretl_corr_rsq:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: data series.
 *
 * Returns: the square of the correlation coefficient for the series
 * @x and @y from obs @t1 to obs @t2, skipping any missing values,
 * or #NADBL on failure.  Used as alternative value for R^2 in a
 * regression without an intercept.
 */

double gretl_corr_rsq (int t1, int t2, const double *x, const double *y)
{
    double r = gretl_corr(t1, t2, x, y, NULL);

    return (na(r))? r : r * r;
}

/* we're supposing a variance smaller than this is just noise */
#define TINYVAR 1.0e-36

double gretl_skewness (int t1, int t2, const double *x)
{
    double s, sd, xx, xbar, ret;
    int t, n = 0;

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
	return NADBL;
    }

    s = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    xx = x[t] - xbar;
	    s += xx * xx;
	    n++;
	}
    }

    s /= n;

    if (s <= TINYVAR) {
	ret = NADBL;
    } else {
	sd = sqrt(s);
	s = 0.0;

	for (t=t1; t<=t2; t++) {
	    if (!na(x[t])) {
		xx = (x[t] - xbar) / sd;
		s += xx * xx * xx;
	    }
	}

	ret = s/n;
    }

    return ret;
}

/* note: the following computes EXCESS kurtosis */

double gretl_kurtosis (int t1, int t2, const double *x)
{
    double s, sd, xx, xbar, ret;
    int t, n = 0;

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
	return NADBL;
    }

    s = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    xx = x[t] - xbar;
	    s += xx * xx;
	    n++;
	}
    }

    s /= n;

    if (s <= TINYVAR) {
	ret = NADBL;
    } else {
	sd = sqrt(s);
	s = 0.0;

	for (t=t1; t<=t2; t++) {
	    if (!na(x[t])) {
		xx = (x[t] - xbar) / sd;
		s += xx * xx * xx * xx;
	    }
	}

	ret = s/n - 3.0;
    }

    return ret;
}

/**
 * gretl_moments:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @wts: weights (may be NULL).
 * @xbar: pointer to receive mean.
 * @sd: pointer to receive standard deviation.
 * @skew: pointer to receive skewness (may be NULL).
 * @kurt: pointer to receive excess kurtosis (may be NULL).
 * @k: degrees of freedom loss (generally 1).
 *
 * Calculates sample moments for series @x from obs @t1 to obs
 * @t2.
 *
 * Returns: 0 on success, 1 on error.
 */

int gretl_moments (int t1, int t2, const double *x,
		   const double *wts,
		   double *xbar, double *std,
		   double *skew, double *kurt, int k)
{
    int t, n;
    double dev, var;
    double s, s2, s3, s4;
    int allstats = 1;
    int weighted = (wts != NULL);
    double wt, wn = 0;

    if (skew == NULL && kurt == NULL) {
	allstats = 0;
    }

    while (na(x[t1]) && t1 <= t2) {
	/* skip missing values at start of sample */
	t1++;
    }

    if (gretl_isconst(t1, t2, x)) {
	/* no variation in x */
	*xbar = x[t1];
	*std = 0.0;
	if (allstats) {
	    *skew = *kurt = NADBL;
	}
	return 1;
    }

    /* get number and sum of valid observations */
    s = 0.0;
    n = 0;
    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    if (weighted) {
		wt = wts[t];
		if (!na(wt) && wt != 0.0) {
		    s += wt * x[t];
		    wn += wt;
		    n++;
		}
	    } else {
		s += x[t];
		n++;
	    }
	}
    }

    if (n == 0) {
	/* no valid data */
	*xbar = *std = NADBL;
	if (allstats) {
	    *skew = *kurt = 0.0;
	}
	return 1;
    }

    if (!weighted) {
	wn = n;
    }

    /* calculate mean and initialize other stats */
    *xbar = s / wn;
    var = 0.0;
    if (allstats) {
	*skew = *kurt = 0.0;
    }

    /* compute quantities needed for higher moments */
    s2 = s3 = s4 = 0.0;
    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    dev = x[t] - *xbar;
	    if (weighted) {
		wt = wts[t];
		if (!na(wt) && wt != 0.0) {
		    s2 += wt * dev * dev;
		    if (allstats) {
			s3 += wt * pow(dev, 3);
			s4 += wt * pow(dev, 4);
		    }
		}
	    } else {
		s2 += dev * dev;
		if (allstats) {
		    s3 += pow(dev, 3);
		    s4 += pow(dev, 4);
		}
	    }
	}
    }


    /* variance */
    if (weighted) {
	/* as per Stata's -summary- with aweights */
	var = n * s2 / (wn * (n-1));
    } else {
	var = s2 / (wn - k);
    }

    if (var < 0.0) {
	/* in principle impossible, but this is a computer */
	*std = NADBL;
	if (allstats) {
	    *skew = *kurt = NADBL;
	}
	return 1;
    }

    if (var > TINYVAR) {
	*std = sqrt(var);
    } else {
	/* screen out minuscule variance */
	*std = 0.0;
    }

    if (allstats) {
	if (var > TINYVAR) {
	    *skew = (s3 / wn) / pow(s2 / wn, 1.5);
	    *kurt = ((s4 / wn) / pow(s2 / wn, 2)) - 3.0; /* excess kurtosis */
	} else {
	    /* if variance is effectively zero, these should be undef'd */
	    *skew = *kurt = NADBL;
	}
    }

    return 0;
}

/**
 * gretl_sorted_series:
 * @v: ID number of input series.
 * @dset: dataset struct.
 * @opt: may include OPT_M to flag an error in case
 * missing values are found.
 * @n: on input, the minimum acceptable number of
 * non-missing observations; on output, the number
 * of valid observations.
 * @err: location to receive error code.
 *
 * Returns: an array containing the valid values of the
 * input series over the sample range given in @dset,
 * sorted from smallest to largest, or NULL on error.
 * An error is flagged if the number of valid observations
 * is less than that given in @n on input, or if OPT_M
 * is given and the input contains missing values.
 */

double *gretl_sorted_series (int v, const DATASET *dset,
			     gretlopt opt, int *n,
			     int *err)
{
    double *y = NULL;
    int t, k = 0;

    if (*n < 1) {
	*n = 1;
    }

    for (t=dset->t1; t<=dset->t2; t++) {
	if (!na(dset->Z[v][t])) {
	    k++;
	} else if (opt & OPT_M) {
	    *err = E_MISSDATA;
	    return NULL;
	}
    }

    if (k < *n) {
	gretl_errmsg_set(_("Insufficient data"));
	*err = E_DATA;
	return NULL;
    }

    y = malloc(k * sizeof *y);
    if (y == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    *n = k;

    k = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
	if (!na(dset->Z[v][t])) {
	    y[k++] = dset->Z[v][t];
	}
    }

    qsort(y, k, sizeof *y, gretl_compare_doubles);

    return y;
}

/**
 * free_freq:
 * @freq: pointer to gretl frequency distribution struct
 *
 * Frees all resources associated with @freq, and the
 * pointer itself.
 */

void free_freq (FreqDist *freq)
{
    if (freq == NULL) {
	return;
    }

    free(freq->midpt);
    free(freq->endpt);
    free(freq->f);

    if (freq->S != NULL) {
	strings_array_free(freq->S, freq->numbins);
    }

    free(freq);
}

static FreqDist *freq_new (const DATASET *dset, int v)
{
    FreqDist *freq;

    freq = malloc(sizeof *freq);
    if (freq == NULL) return NULL;

    strcpy(freq->varname, dset->varname[v]);
    strcpy(freq->gname, plotname(dset, v, 1));

    freq->midpt = NULL;
    freq->endpt = NULL;
    freq->S = NULL;
    freq->f = NULL;

    freq->dist = 0;
    freq->discrete = 0;
    freq->strvals = 0;

    freq->xbar = NADBL;
    freq->sdx = NADBL;
    freq->test = NADBL;

    return freq;
}

/*
 * dh_root_b1_to_z1:
 * @rb1: square root b1, skewness.
 * @n: number of observations.
 *
 * Performs the transformation from skewness, root b1, to the
 * normal score, z1, as set out in Doornik and Hansen, "An Omnibus
 * Test for Normality", 1994.  The transformation is originally
 * due to D'Agostino (Biometrika, 1970).
 *
 * Returns: the z1 value.
 */

static double dh_root_b1_to_z1 (double rb1, double n)
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

/*
 * dh_b2_to_z2:
 * @b1: skewness.
 * @b2: kurtosis.
 * @n: number of observations.
 *
 * Performs the transformation from kurtosis, b2, to the
 * normal score, z2, as set out in Doornik and Hansen, "An Omnibus
 * Test for Normality", 1994.
 *
 * Returns: the z2 value, or #NADBL on failure.
 */

static double dh_b2_to_z2 (double b1, double b2, double n)
{
    double d, a, c, k, alpha, chi, z2;
    double n2 = n * n;

    d = (n-3) * (n+1) * (n2 + 15*n - 4.0);
    a = ((n-2) * (n+5) * (n+7) * (n2 + 27*n - 70.0)) / (6.0 * d);
    c = ((n-7) * (n+5) * (n+7) * (n2 + 2*n - 5.0)) / (6.0 * d);
    k = ((n+5) * (n+7) * (n2*n + 37*n2 + 11*n - 313.0)) / (12.0 * d);

    alpha = a + b1 * c;
    z2 =  (1.0 / (9.0 * alpha)) - 1.0;
    chi = (b2 - 1.0 - b1) * 2.0 * k;

    if (chi > 0.0) {
       z2 += pow(chi/(2*alpha), 1.0 / 3.0);
    }

    z2 *= sqrt(9.0*alpha);

    return z2;
}

/**
 * doornik_chisq:
 * @skew: skewness.
 * @xkurt: excess kurtosis.
 * @n: number of observations.
 *
 * Calculates the Chi-square test for normality as set out by
 * Doornik and Hansen, "An Omnibus Test for Normality", 1994.
 * This is a modified version of the test proposed by Bowman and
 * Shenton (Biometrika, 1975).
 *
 * Returns: the Chi-square value, which has 2 degrees of freedom.
 */

double doornik_chisq (double skew, double xkurt, int n)
{
    double z1, z2, x2 = NADBL;

    z1 = dh_root_b1_to_z1(skew, (double) n);
    z2 = dh_b2_to_z2(skew * skew, xkurt + 3.0, (double) n);

    if (!na(z2)) {
	x2 = z1*z1 + z2*z2;
    }

    return x2;
}

static int
series_get_moments (int t1, int t2, const double *x,
		    double *skew, double *xkurt,
		    int *pn)
{
    double dev, s[4] = {0.0};
    int t, n = 0;
    int err = 0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    s[0] += x[t];
	    n++;
	}
    }

    *pn = n;
    if (n < 3) {
	return E_DATA;
    }

    s[0] /= n;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    dev = x[t] - s[0];
	    s[1] += dev * dev;
	    s[2] += pow(dev, 3);
	    s[3] += pow(dev, 4);
	}
    }

    s[1] /= n;

    if (s[1] > 0.0) {
	*skew = (s[2] / n) / pow(s[1], 1.5);
	*xkurt = ((s[3] / n) / pow(s[1], 2)) - 3.0;
    } else {
	/* if variance is effectively zero, these should be undef'd */
	*skew = *xkurt = NADBL;
	err = 1;
    }

    return err;
}

static int
get_moments (const gretl_matrix *M, int row, double *skew, double *kurt)
{
    int j, n = gretl_matrix_cols(M);
    double dev, s[4] = {0.0};
    int err = 0;

    for (j=0; j<n; j++) {
	s[0] += gretl_matrix_get(M, row, j);
    }

    s[0] /= n;

    for (j=0; j<n; j++) {
	dev = gretl_matrix_get(M, row, j) - s[0];
	s[1] += dev * dev;
	s[2] += pow(dev, 3);
	s[3] += pow(dev, 4);
    }

    s[1] /= n;

    if (s[1] > 0.0) {
	*skew = (s[2] / n) / pow(s[1], 1.5);
	*kurt = ((s[3] / n) / pow(s[1], 2));
    } else {
	/* if variance is effectively zero, these should be undef'd */
	*skew = *kurt = NADBL;
	err = 1;
    }

    return err;
}

int
multivariate_normality_test (const gretl_matrix *E,
			     const gretl_matrix *Sigma,
			     gretlopt opt, PRN *prn)
{
    gretl_matrix *S = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *C = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *evals = NULL;
    gretl_matrix *tmp = NULL;
    /* convenience pointers: do not free! */
    gretl_matrix *H;
    gretl_vector *Z1;
    gretl_vector *Z2;
    double x, skew, kurt;
    double X2 = NADBL;
    int n, p;
    int i, j;
    int err = 0;

    if (E == NULL || Sigma == NULL) {
	err = E_DATA;
	goto bailout;
    }

    p = gretl_matrix_cols(E);
    n = gretl_matrix_rows(E);

    clear_gretl_matrix_err();

    S = gretl_matrix_copy(Sigma);
    V = gretl_vector_alloc(p);
    C = gretl_matrix_alloc(p, p);
    X = gretl_matrix_copy_transpose(E);
    R = gretl_matrix_alloc(p, n);
    tmp = gretl_matrix_alloc(p, p);

    err = get_gretl_matrix_err();
    if (err) {
	goto bailout;
    }

    for (i=0; i<p; i++) {
	x = gretl_matrix_get(S, i, i);
	gretl_vector_set(V, i, 1.0 / sqrt(x));
    }

    err = gretl_matrix_diagonal_sandwich(V, S, C);
    if (err) {
	goto bailout;
    }

    if (!(opt & OPT_Q)) {
	pputc(prn, '\n');
	gretl_matrix_print_to_prn(C, _("Residual correlation matrix, C"), prn);
    }

    evals = gretl_symmetric_matrix_eigenvals(C, 1, &err);
    if (err) {
	goto bailout;
    }

    if (!(opt & OPT_Q)) {
	pprintf(prn, "%s\n\n", _("Eigenvalues of C"));
	for (i=0; i<p; i++) {
	    pprintf(prn, " %10g\n", evals->val[i]);
	}
	pputc(prn, '\n');
    }

    /* C should now contain eigenvectors of the original C:
       relabel as 'H' for perspicuity */
    H = C;
#if 0
    gretl_matrix_print_to_prn(H, "Eigenvectors, H", prn);
#endif
    gretl_matrix_copy_values(tmp, H);

    /* make "tmp" into $H \Lambda^{-1/2}$ */
    for (i=0; i<p; i++) {
	for (j=0; j<p; j++) {
	    x = gretl_matrix_get(tmp, i, j);
	    x *= 1.0 / sqrt(evals->val[j]);
	    gretl_matrix_set(tmp, i, j, x);
	}
    }

    /* make S into $H \Lambda^{-1/2} H'$ */
    gretl_matrix_multiply_mod(tmp, GRETL_MOD_NONE,
			      H, GRETL_MOD_TRANSPOSE,
			      S, GRETL_MOD_NONE);

#if 1
    gretl_matrix_demean_by_row(X);
#endif

    /* compute VX', in X (don't need to subtract means, because these
       are OLS residuals)
    */
    for (i=0; i<p; i++) {
	for (j=0; j<n; j++) {
	    x = gretl_matrix_get(X, i, j);
	    x *= gretl_vector_get(V, i);
	    gretl_matrix_set(X, i, j, x);
	}
    }

    /* finally, compute $R' = H  \Lambda^{-1/2} H' V X'$
       Doornik and Hansen, 1994, section 3 */
    gretl_matrix_multiply(S, X, R);

    /* Z_1 and Z_2 are p-vectors: use existing storage */
    Z1 = V;
    Z2 = gretl_matrix_reuse(tmp, p, 1);

    for (i=0; i<p && !err; i++) {
	get_moments(R, i, &skew, &kurt);
	if (na(skew) || na(kurt)) {
	    err = 1;
	} else {
	    double z1i = dh_root_b1_to_z1(skew, n);
	    double z2i = dh_b2_to_z2(skew * skew, kurt, n);

	    if (na(z2i)) {
		X2 = NADBL;
		err = E_NAN;
	    } else {
		gretl_vector_set(Z1, i, z1i);
		gretl_vector_set(Z2, i, z2i);
	    }
	}
    }

    if (!err) {
	X2 = gretl_vector_dot_product(Z1, Z1, &err);
	X2 += gretl_vector_dot_product(Z2, Z2, &err);
    }

    if (na(X2)) {
	pputs(prn, "Calculation of test statistic failed\n");
    } else {
	/* print and record result */
	double pv = chisq_cdf_comp(2 * p, X2);

	if (!(opt & OPT_I)) {
	    pputs(prn, _("Doornik-Hansen test"));
	    pprintf(prn, "\n %s(%d) = %g [%.4f]\n\n", _("Chi-square"), 2 * p,
		    X2, pv);
	}
	record_test_result(X2, pv);
    }

 bailout:

    gretl_matrix_free(S);
    gretl_matrix_free(V);
    gretl_matrix_free(C);
    gretl_matrix_free(X);
    gretl_matrix_free(R);
    gretl_matrix_free(evals);
    gretl_matrix_free(tmp);

    return err;
}

static int freq_add_arrays (FreqDist *freq, int n)
{
    int err = 0;

    if (!freq->discrete) {
	freq->endpt = malloc((n + 1) * sizeof *freq->endpt);
	if (freq->endpt == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	if (freq->strvals) {
	    freq->S = strings_array_new(n);
	    if (freq->S == NULL) {
		err = E_ALLOC;
	    }
	} else {
	    freq->midpt = malloc(n * sizeof *freq->midpt);
	    if (freq->midpt == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (!err) {
	freq->f = malloc(n * sizeof *freq->f);
	if (freq->f == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	freq->numbins = n;
    }

    return err;
}

static int freq_n_bins (int k0, int n, double s, double r)
{
    int k;

#if 0 /* maybe enable this? */
    /* Scott's bin width rule */
    double h = 3.49 * s * pow((double) n, -1.0/3.0);
    fprintf(stderr, "s = %g, n = %d, Scott's h = %g\n", s, n, h);
    fprintf(stderr, " => nbins ~= %d\n", (int) (r/h));
#endif

    if (k0 > 0) {
	/* got a pre-specified value */
	k = k0;
    } else if (n < 16) {
	k = 5;
    } else if (n < 50) {
	k = 7;
    } else if (n > 850) {
	k = 29;
    } else {
	k = (int) sqrt((double) n);
	if (k > 99) {
	    k = 99;
	}
    }
    if (k % 2 == 0) {
	k++;
    }

    return k;
}

int freq_setup (int v, const DATASET *dset, int *pn,
		double *pxmax, double *pxmin, int *nbins,
		double *binwidth)
{
    const double *x = dset->Z[v];
    gretl_matrix *xvals;
    int all_n = sample_size(dset);
    int nv, missvals = 0;
    int n = 0;
    int err = 0;

    xvals = gretl_matrix_values_full(x + dset->t1, all_n, OPT_S,
				     &missvals, &err);
    if (!err) {
	nv = xvals->rows;
	n = all_n - missvals;
	if (nv < 2) {
	    gretl_errmsg_sprintf(_("%s is a constant"), dset->varname[v]);
	    err = E_DATA;
	} else if (n < 8) {
	    gretl_errmsg_sprintf(_("Insufficient data to build frequency "
				   "distribution for variable %s"),
				 dset->varname[v]);
	    err = E_TOOFEW;
	}
    }

    if (!err) {
	double xmin = xvals->val[0];
	double xmax = xvals->val[nv-1];
	double xrange = xmax - xmin;
	double s = gretl_stddev(dset->t1, dset->t2, x);
	int k = freq_n_bins(*nbins, n, s, xrange);

	if (nv < k) {
	    /* too few distinct values: don't do binning */
	    *nbins = 0;
	} else {
	    *nbins = k;
	    *binwidth = xrange / (k - 1);
	}
	*pn = n;
	*pxmax = xmax;
	*pxmin = xmin;
    }

    gretl_matrix_free(xvals);

    return err;
}

/* calculate test stat for distribution, if the sample
   is big enough */

static void
freq_dist_stat (FreqDist *freq, const double *x, gretlopt opt, int k)
{
    double skew, kurt;

    freq->test = NADBL;
    freq->dist = 0;

    gretl_moments(freq->t1, freq->t2, x, NULL,
		  &freq->xbar, &freq->sdx,
		  &skew, &kurt, k);

    if (freq->n > 7) {
	if (opt & OPT_O) {
            int err = 0;

	    if (freq->n > 500) {
		freq->test = vge_gamma_test(x, freq->t1, freq->t2, &err);
	    } else {
		freq->test = lockes_test(x, freq->t1, freq->t2, &err);
	    }
	    freq->dist = D_GAMMA;
	} else if (opt & OPT_Z) {
	    freq->test = doornik_chisq(skew, kurt, freq->n);
	    freq->dist = D_NORMAL;
	}
    }
}

struct strval_sorter {
    char *s;
    int n;
};

int compare_strvals (const void *a, const void *b)
{
    const struct strval_sorter *sa = a;
    const struct strval_sorter *sb = b;

    return strcmp(sa->s, sb->s);
}

FreqDist *get_string_freq (int v, const DATASET *dset,
			   int *err)
{
    FreqDist *freq;
    struct strval_sorter *ss;
    const double *x = dset->Z[v];
    series_table *st;
    char **S;
    int i, t, n, ns;

    freq = freq_new(dset, v);
    if (freq == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    st = series_get_string_table(dset, v);
    S = series_table_get_strings(st, &ns);

    ss = malloc(ns * sizeof *ss);
    if (ss == NULL) {
	*err = E_ALLOC;
	free(freq);
	return NULL;
    }

    for (i=0; i<ns; i++) {
	ss[i].s = S[i];
	ss[i].n = 0;
    }

    n = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
	if (!na(x[t])) {
	    i = x[t] - 1;
	    ss[i].n += 1;
	    n++;
	}
    }

    qsort(ss, ns, sizeof *ss, compare_strvals);

    freq->t1 = dset->t1;
    freq->t2 = dset->t2;
    freq->n = n;
    freq->discrete = 1;
    freq->strvals = 1;

    if (freq_add_arrays(freq, ns)) {
	*err = E_ALLOC;
    } else {
	for (i=0; i<ns; i++) {
	    freq->S[i] = gretl_strdup(ss[i].s);
	    freq->f[i] = ss[i].n;
	}
    }

    free(ss);

    if (*err && freq != NULL) {
	free_freq(freq);
	freq = NULL;
    }

    return freq;
}

FreqDist *get_discrete_freq (int v, const DATASET *dset,
			     gretlopt opt, int *err)
{
    FreqDist *freq;
    const double *x = dset->Z[v];
    int *ifreq = NULL;
    double *ivals = NULL;
    double *sorted = NULL;
    double last;
    int i, t, nv;

    freq = freq_new(dset, v);
    if (freq == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    freq->t1 = dset->t1;
    freq->t2 = dset->t2;

    freq->n = 0;
    for (t=freq->t1; t<=freq->t2; t++) {
	if (!na(x[t])) {
	    freq->n += 1;
	}
    }

    if (freq->n < 3) {
	gretl_errmsg_sprintf(_("Insufficient data to build frequency "
			       "distribution for variable %s"),
			     dset->varname[v]);
	*err = E_TOOFEW;
	goto bailout;
    }

    freq->discrete = 1;
    freq->test = NADBL;
    freq->dist = 0;

    sorted = malloc(freq->n * sizeof *sorted);
    if (sorted == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    i = 0;
    for (t=freq->t1; t<=freq->t2; t++) {
	if (!na(x[t])) {
	    sorted[i++] = x[t];
	}
    }

    qsort(sorted, freq->n, sizeof *sorted, gretl_compare_doubles);
    nv = count_distinct_values(sorted, freq->n);

    if (nv >= 10 && !(opt & OPT_X)) {
	freq_dist_stat(freq, x, opt, 1);
    } else if (opt & (OPT_Z | OPT_O)) {
	freq->xbar = gretl_mean(freq->t1, freq->t2, x);
	freq->sdx = gretl_stddev(freq->t1, freq->t2, x);
    }

    ifreq = malloc(nv * sizeof *ifreq);
    ivals = malloc(nv * sizeof *ivals);
    if (ifreq == NULL || ivals == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    ivals[0] = last = sorted[0];
    ifreq[0] = i = 1;

    for (t=1; t<freq->n; t++) {
	if (sorted[t] != last) {
	    last = sorted[t];
	    ifreq[i] = 1;
	    ivals[i++] = last;
	} else {
	    ifreq[i-1] += 1;
	}
    }

    if (freq_add_arrays(freq, nv)) {
	*err = E_ALLOC;
    } else {
	int allints = 1;

	for (i=0; i<nv; i++) {
	    if (allints && ivals[i] != floor(ivals[i])) {
		allints = 0;
	    }
	    freq->midpt[i] = ivals[i];
	    freq->f[i] = ifreq[i];
	}

	if (allints) {
	    freq->discrete = 2;
	}
    }

 bailout:

    free(sorted);
    free(ivals);
    free(ifreq);

    if (*err && freq != NULL) {
	free_freq(freq);
	freq = NULL;
    }

    return freq;
}

/**
 * get_freq:
 * @varno: ID number of variable to process.
 * @dset: dataset struct.
 * @fmin: lower limit of left-most bin (or #NADBL for automatic).
 * @fwid: bin width (or #NADBL for automatic).
 * @nbins: number of bins to use (or 0 for automatic).
 * @params: degrees of freedom loss (generally = 1 unless we're dealing
 * with the residual from a regression).
 * @opt: if includes %OPT_Z, set up for comparison with normal dist;
 * if includes %OPT_O, compare with gamma distribution;
 * if includes %OPT_S we're not printing results; if includes %OPT_D,
 * treat the variable as discrete; %OPT_X indicates that this function
 * is called as part of a cross-tabulation; if includes %OPT_G,
 * add the arrays that are needed for a plot even when we're not
 * printing (%OPT_S).
 * @err: location to receive error code.
 *
 * Calculates the frequency distribution for the specified variable.
 *
 * Returns: pointer to struct containing the distribution.
 */

FreqDist *get_freq (int varno, const DATASET *dset,
		    double fmin, double fwid,
		    int nbins, int params,
		    gretlopt opt, int *err)
{
    FreqDist *freq;
    const double *x;
    double xx, xmin, xmax;
    double binwidth = fwid;
    int t, k, n;

    if (is_string_valued(dset, varno)) {
	return get_string_freq(varno, dset, err);
    } else if (series_is_discrete(dset, varno) || (opt & OPT_D)) {
	return get_discrete_freq(varno, dset, opt, err);
    } else if (gretl_isdiscrete(dset->t1, dset->t2, dset->Z[varno]) > 1) {
	return get_discrete_freq(varno, dset, opt, err);
    }

    freq = freq_new(dset, varno);
    if (freq == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    *err = freq_setup(varno, dset, &n, &xmax, &xmin, &nbins, &binwidth);

    if (!*err && nbins == 0) {
	/* no binning wanted: switch to discrete variant */
	free_freq(freq);
	return get_discrete_freq(varno, dset, opt, err);
    }

    if (*err) {
	goto bailout;
    }

    if (!na(fmin) && !na(fwid)) {
	/* endogenous implied number of bins */
	nbins = (int) ceil((xmax - fmin) / fwid);
	if (nbins <= 0 || nbins > 5000) {
	    *err = E_INVARG;
	    goto bailout;
	} else {
	    binwidth = fwid;
	}
    }

    freq->t1 = dset->t1;
    freq->t2 = dset->t2;
    freq->n = n;

    x = dset->Z[varno];
    freq_dist_stat(freq, x, opt, params);

    if (freq_add_arrays(freq, nbins)) {
	*err = E_ALLOC;
	goto bailout;
    }

    if (na(fmin) || na(fwid)) {
	freq->endpt[0] = xmin - .5 * binwidth;
	if (xmin > 0.0 && freq->endpt[0] < 0.0) {
	    double rshift;

	    freq->endpt[0] = 0.0;
	    rshift = 1.0 - xmin / binwidth;
	    freq->endpt[freq->numbins] = xmax + rshift * binwidth;
	} else {
	    freq->endpt[freq->numbins] = xmax + .5 * binwidth;
	}
    } else {
	freq->endpt[0] = fmin;
	freq->endpt[freq->numbins] = fmin + nbins * fwid;
    }

    for (k=0; k<freq->numbins; k++) {
	freq->f[k] = 0;
	if (k > 0) {
	    freq->endpt[k] = freq->endpt[k-1] + binwidth;
	}
	freq->midpt[k] = freq->endpt[k] + .5 * binwidth;
    }

    for (t=dset->t1; t<=dset->t2; t++) {
	xx = x[t];
	if (na(xx)) {
	    continue;
	}
	if (xx < freq->endpt[1]) {
	    freq->f[0] += 1;
	} else if (xx >= freq->endpt[freq->numbins]) {
	    freq->f[freq->numbins - 1] += 1;
	    continue;
	} else {
	    for (k=1; k<freq->numbins; k++) {
		if (freq->endpt[k] <= xx && xx < freq->endpt[k+1]) {
		    freq->f[k] += 1;
		    break;
		}
	    }
	}
    }

 bailout:

    if (*err && freq != NULL) {
	free_freq(freq);
	freq = NULL;
    }

    return freq;
}

static void record_freq_test (const FreqDist *freq)
{
    double pval = NADBL;

    if (freq->dist == D_NORMAL) {
	pval = chisq_cdf_comp(2, freq->test);
    } else if (freq->dist == D_GAMMA) {
	pval = normal_pvalue_2(freq->test);
    }

    if (!na(pval)) {
	record_test_result(freq->test, pval);
    }
}

static int check_freq_opts (gretlopt opt, int *n_bins,
			    double *fmin, double *fwid)
{
    double x;
    int err = 0;

    if (opt & OPT_N) {
	/* number of bins given */
	if (opt & (OPT_M | OPT_W)) {
	    /* can't give min, width spec */
	    return E_BADOPT;
	} else {
	    int n = get_optval_int(FREQ, OPT_N, &err);

	    if (!err) {
		if (n < 2 || n > 10000) {
		    err = E_INVARG;
		} else {
		    *n_bins = n;
		}
	    }
	    if (err) {
		return err;
	    }
	}
    }

    if (opt & OPT_M) {
	/* minimum specified */
	if (!(opt & OPT_W)) {
	    /* but no width given */
	    return E_ARGS;
	} else {
	    x = get_optval_double(FREQ, OPT_M, &err);
	    if (err) {
		return err;
	    } else if (na(x)) {
		return E_ARGS;
	    } else {
		*fmin = x;
	    }
	}
    }

    if (opt & OPT_W) {
	/* width specified */
	if (!(opt & OPT_M)) {
	    /* but no min given */
	    return E_ARGS;
	} else {
	    x = get_optval_double(FREQ, OPT_W, &err);
	    if (err) {
		return err;
	    } else if (na(x)) {
		return E_ARGS;
	    } else if (x <= 0) {
		return E_INVARG;
	    } else {
		*fwid = x;
	    }
	}
    }

    return 0;
}

static void record_freq_matrix (FreqDist *fd)
{
    gretl_matrix *m = NULL;
    char **Sr = NULL;
    int i, n = fd->numbins;

    if (fd->S != NULL) {
	m = gretl_matrix_alloc(n, 1);
	Sr = strings_array_dup(fd->S, n);
    } else {
	m = gretl_matrix_alloc(n, 2);
    }

    if (m != NULL) {
	char **Sc = NULL;

	if (fd->S == NULL) {
	    Sc = strings_array_new(2);
	    if (Sc != NULL) {
		Sc[0] = gretl_strdup("midpoint");
		Sc[1] = gretl_strdup("count");
	    }
	}
	for (i=0; i<n; i++) {
	    if (fd->S != NULL) {
		gretl_matrix_set(m, i, 0, fd->f[i]);
	    } else {
		gretl_matrix_set(m, i, 0, fd->midpt[i]);
		gretl_matrix_set(m, i, 1, fd->f[i]);
	    }
	}
	if (Sc != NULL) {
	    gretl_matrix_set_colnames(m, Sc);
	}
	if (Sr != NULL) {
	    gretl_matrix_set_rownames(m, Sr);
	}
	set_last_result_data(m, GRETL_TYPE_MATRIX);
    }
}

static int freq_plot_wanted (PlotType ptype, gretlopt opt,
			     int *err)
{
    if (opt & OPT_Q) {
	/* handle legacy option: --quiet = no plot */
	return 0;
    } else {
	return gnuplot_graph_wanted(ptype, opt, err);
    }
}

/* Wrapper function: get the distribution, print it if
   wanted, graph it if wanted, then free stuff.
*/

int freqdist (int varno, const DATASET *dset,
	      gretlopt opt, PRN *prn)
{
    FreqDist *freq = NULL;
    DistCode dist = D_NONE;
    PlotType ptype = 0;
    double fmin = NADBL;
    double fwid = NADBL;
    int n_bins = 0;
    int do_graph = 0;
    int err = 0;

    if (opt & OPT_O) {
	dist = D_GAMMA;
	ptype = PLOT_FREQ_GAMMA;
    } else if (opt & OPT_Z) {
	dist = D_NORMAL;
	ptype = PLOT_FREQ_NORMAL;
    } else {
	ptype = PLOT_FREQ_SIMPLE;
    }

    do_graph = freq_plot_wanted(ptype, opt, &err);

    if (!err) {
	err = check_freq_opts(opt, &n_bins, &fmin, &fwid);
    }

    if (!err) {
	gretlopt fopt = opt;

	if (do_graph) {
	    fopt |= OPT_G;
	}
	freq = get_freq(varno, dset, fmin, fwid, n_bins, 1, fopt, &err);
    }

    if (!err) {
	if (!(opt & OPT_S)) {
	    print_freq(freq, varno, dset, prn);
	} else if (dist) {
	    record_freq_test(freq);
	}
	record_freq_matrix(freq);

	if (do_graph && freq->numbins < 2) {
	    do_graph = 0;
	}
	if (do_graph) {
	    int gerr = plot_freq(freq, dist, opt);

	    if (gerr) {
		pputs(prn, _("gnuplot command failed\n"));
		do_graph = 0;
	    }
	}
	free_freq(freq);
    }

    return err;
}

/* Wrapper function: get the frequency distribution and write it into
   a gretl_matrix: the first column holds the mid-points of the bins
   and the second holds the frequencies.
*/

gretl_matrix *freqdist_matrix (const double *x, int t1, int t2, int *err)
{
    DATASET *dset = NULL;
    FreqDist *freq = NULL;
    gretl_matrix *m = NULL;
    int i, t, T = t2 - t1 + 1;

    dset = create_auxiliary_dataset(1, T, 0);
    if (dset == NULL) {
	*err = E_ALLOC;
    }

    if (!*err) {
	i = 0;
	for (t=t1; t<=t2; t++) {
	    dset->Z[0][i++] = x[t];
	}
	freq = get_freq(0, dset, NADBL, NADBL, 0, 1, OPT_NONE, err);
    }

    if (!*err) {
	m = gretl_matrix_alloc(freq->numbins, 2);
	if (m == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	for (i=0; i<freq->numbins; i++) {
	    gretl_matrix_set(m, i, 0, freq->midpt[i]);
	    gretl_matrix_set(m, i, 1, freq->f[i]);
	}
    }

    destroy_dataset(dset);
    free_freq(freq);

    return m;
}

static inline double alt_mget (const gretl_matrix *m,
			       int i, int j, int tr)
{
    return tr ? m->val[i*m->rows+j] : m->val[j*m->rows+i];
}

static int corresp_status (const gretl_matrix *H, int tr)
{
    int dim1 = tr ? H->rows : H->cols;
    int dim2 = tr ? H->cols : H->rows;
    int sum1, max1 = 0;
    int sum2, max2 = 0;
    int i, j;
    int ret = 0;

    /* find the max sum in dimension 1 */
    for (j=0; j<dim1 && max1 < 2; j++) {
	sum1 = 0;
	for (i=0; i<dim2 && sum1 < 2; i++) {
	    sum1 += alt_mget(H, i, j, tr) > 0;
	}
	if (sum1 > max1) {
	    max1 = sum1;
	}
    }
    if (max1 == 1) {
	ret = tr ? -1 : 1;
	/* find the max sum in dimension 2 */
	for (i=0; i<dim2 && max2 < 2; i++) {
	    sum2 = 0;
	    for (j=0; j<dim1 && sum2 < 2; j++) {
		sum2 += alt_mget(H, i, j, tr) > 0;
	    }
	    if (sum2 > max2) {
		max2 = sum2;
	    }
	}
	if (max2 == 1) {
	    ret = 2;
	}
    }

    return ret;
}

int correspondence (const double *x, const double *y,
		    int n, int *err)
{
    gretl_matrix *H = NULL;
    int status = 0;

    /* status: 2 means 1-to-1 x:y relationship
               1 means 1-to-n x:y relationship
	      -1 means n-to-1 x:y relationship
	       0 means no relationship
    */

    H = gretl_matrix_xtab(x, y, n, err);

    if (!*err) {
	status = corresp_status(H, 0);
	if (status == 0) {
	    status = corresp_status(H, 1);
	}
    }

    gretl_matrix_free(H);

    return status;
}

gretl_matrix *xtab_to_matrix (const Xtab *tab)
{
    gretl_matrix *m;
    double x;
    int i, j;

    if (tab == NULL) {
	return NULL;
    }

    m = gretl_matrix_alloc(tab->rows, tab->cols);
    if (m == NULL) {
	return NULL;
    }

    for (j=0; j<tab->cols; j++) {
	for (i=0; i<tab->rows; i++) {
	    x = (double) tab->f[i][j];
	    gretl_matrix_set(m, i, j, x);
	}
    }

    return m;
}

static void *last_result;
static GretlType last_result_type;

void *get_last_result_data (GretlType *type, int *err)
{
    void *ret = NULL;

    if (last_result == NULL) {
	*type = GRETL_TYPE_NONE;
	*err = E_BADSTAT;
    } else {
	*type = last_result_type;
	if (*type == GRETL_TYPE_MATRIX) {
	    ret = gretl_matrix_copy(last_result);
	} else {
	    ret = gretl_bundle_copy(last_result, err);
	}
    }

    return ret;
}

void set_last_result_data (void *data, GretlType type)
{
    if (last_result != NULL) {
	if (last_result_type == GRETL_TYPE_MATRIX) {
	    gretl_matrix_free(last_result);
	} else if (last_result_type == GRETL_TYPE_BUNDLE) {
	    gretl_bundle_destroy(last_result);
	}
    }

    last_result = data;
    last_result_type = type;
}

void last_result_cleanup (void)
{
    if (last_result != NULL) {
	if (last_result_type == GRETL_TYPE_MATRIX) {
	    gretl_matrix_free(last_result);
	} else if (last_result_type == GRETL_TYPE_BUNDLE) {
	    gretl_bundle_destroy(last_result);
	}
    }

    last_result = NULL;
    last_result_type = 0;
}

/**
 * free_xtab:
 * @tab: pointer to gretl crosstab struct.
 *
 * Frees all resources associated with @tab, and the
 * pointer itself.
 */

void free_xtab (Xtab *tab)
{
    int i;

    if (tab == NULL) {
	return;
    }

    free(tab->rtotal);
    free(tab->ctotal);
    free(tab->rval);
    free(tab->cval);

    if (tab->f != NULL) {
	for (i=0; i<tab->rows; i++) {
	    free(tab->f[i]);
	}
	free(tab->f);
    }

    if (tab->Sr != NULL) {
	strings_array_free(tab->Sr, tab->rows);
    }
    if (tab->Sc != NULL) {
	strings_array_free(tab->Sc, tab->cols);
    }

    free(tab);
}

static Xtab *xtab_new (int n, int t1, int t2)
{
    Xtab *tab = malloc(sizeof *tab);

    if (tab == NULL) return NULL;

    tab->rtotal = NULL;
    tab->ctotal = NULL;
    tab->rval = NULL;
    tab->cval = NULL;
    tab->f = NULL;

    tab->n = n;
    tab->t1 = t1;
    tab->t2 = t2;
    tab->missing = 0;

    *tab->rvarname = '\0';
    *tab->cvarname = '\0';
    tab->Sr = NULL;
    tab->Sc = NULL;
    tab->rstrs = 0;
    tab->cstrs = 0;

    return tab;
}

/* allocate the required arrays (apart from S, which must be
   handled separately) and initialize counters to zero
*/

static int xtab_allocate_arrays (Xtab *tab)
{
    int i, j, err = 0;

    if (!tab->rstrs && tab->rval == NULL) {
	tab->rval = malloc(tab->rows * sizeof *tab->rval);
	if (tab->rval == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && !tab->cstrs && tab->cval == NULL) {
	tab->cval = malloc(tab->cols * sizeof *tab->cval);
	if (tab->cval == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	tab->rtotal = malloc(tab->rows * sizeof *tab->rtotal);
	tab->ctotal = malloc(tab->cols * sizeof *tab->ctotal);
	if (tab->rtotal == NULL || tab->ctotal == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<tab->rows; i++) {
		tab->rtotal[i] = 0;
	    }
	    for (j=0; j<tab->cols; j++) {
		tab->ctotal[j] = 0;
	    }
	}
    }

    if (!err) {
	tab->f = malloc(tab->rows * sizeof *tab->f);
	if (tab->f == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<tab->rows; i++) {
		tab->f[i] = NULL;
	    }
	}
    }

    for (i=0; i<tab->rows && !err; i++) {
	tab->f[i] = malloc(tab->cols * sizeof *tab->f[i]);
	if (tab->f[i] == NULL) {
	    err = E_ALLOC;
	} else {
	    for (j=0; j<tab->cols; j++) {
		tab->f[i][j] = 0;
	    }
	}
    }

    return err;
}

/* also used in gretl_matrix.c */

int compare_xtab_rows (const void *a, const void *b)
{
    const double **da = (const double **) a;
    const double **db = (const double **) b;
    double ret = da[0][0] - db[0][0];

    if (ret == 0) {
	ret = da[0][1] - db[0][1];
    }

    return ret < 0 ? -1 : ret > 0 ? 1: 0;
}

static int xtab_get_data (Xtab *tab, int v, int j,
			  const DATASET *dset,
			  series_table **pst)
{
    double **xtarg = (j == 1)? &tab->cval : &tab->rval;
    char ***Starg = (j == 1)? &tab->Sc : &tab->Sr;
    int *itarg = (j == 1)? &tab->cols : &tab->rows;
    int *ttarg = (j == 1)? &tab->cstrs : &tab->rstrs;
    double *x = dset->Z[v] + dset->t1;
    int n = sample_size(dset);
    gretl_matrix *u;
    int err = 0;

    u = gretl_matrix_values(x, n, OPT_S, &err);

    if (!err && is_string_valued(dset, v)) {
	series_table *st;
	int ns, nv = u->rows;
	char **S0, **S = NULL;

	*pst = st = series_get_string_table(dset, v);
	S0 = series_table_get_strings(st, &ns);

	if (ns > nv) {
	    int i, k;

	    S = strings_array_new(nv);
	    if (S != NULL) {
		for (i=0; i<nv; i++) {
		    k = u->val[i] - 1;
		    S[i] = gretl_strdup(S0[k]);
		}
	    }
	} else {
	    S = strings_array_dup(S0, ns);
	}
	if (S == NULL) {
	    err = E_ALLOC;
	} else {
	    *ttarg = 1;
	    *itarg = nv;
	    *Starg = S;
	    strings_array_sort(Starg, &nv, OPT_NONE);
	}
    } else if (!err) {
	*itarg = u->rows;
	*xtarg = u->val;
	u->val = NULL;
    }

    gretl_matrix_free(u);

    return err;
}

static int xtab_row_match (Xtab *tab, int i, const char *s, double x)
{
    if (tab->rstrs) {
	return strcmp(s, tab->Sr[i]) == 0;
    } else {
	return x == tab->rval[i];
    }
}

static int xtab_col_match (Xtab *tab, int j, const char *s, double x)
{
    if (tab->cstrs) {
	return strcmp(s, tab->Sc[j]) == 0;
    } else {
	return x == tab->cval[j];
    }
}

#define complete_obs(x,y,t) (!na(x[t]) && !na(y[t]))

/* crosstab struct creation functions */

static Xtab *get_new_xtab (int rv, int cv, const DATASET *dset,
			   int *err)
{
    series_table *sti = NULL;
    series_table *stj = NULL;
    const char *s1 = NULL;
    const char *s2 = NULL;
    double x1 = 0, x2 = 0;
    Xtab *tab = NULL;
    int imatch, jmatch;
    int i, j, t, n = 0;

    for (t=dset->t1; t<=dset->t2; t++) {
	if (complete_obs(dset->Z[rv], dset->Z[cv], t)) {
	    n++;
	}
    }

    if (n == 0) {
	*err = E_MISSDATA;
    } else {
	tab = xtab_new(n, dset->t1, dset->t2);
	if (tab == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	/* assemble row data */
	*err = xtab_get_data(tab, rv, 0, dset, &sti);
    }
    if (!*err) {
	/* assemble column data */
	*err = xtab_get_data(tab, cv, 1, dset, &stj);
    }
    if (!*err) {
	*err = xtab_allocate_arrays(tab);
    }

    if (*err) goto bailout;

    tab->missing = (dset->t2 - dset->t1 + 1) - n;
    strcpy(tab->rvarname, dset->varname[rv]);
    strcpy(tab->cvarname, dset->varname[cv]);

    /* The following could be made more efficient by substituting
       sorted arrays for dset->Z[rv] and dset->Z[cv] but I'm not
       sure if the fixed cost would be recouped.
    */

    for (t=dset->t1; t<=dset->t2; t++) {
	if (!complete_obs(dset->Z[rv], dset->Z[cv], t)) {
	    continue;
	}
	x1 = dset->Z[rv][t];
	if (tab->rstrs) {
	    s1 = series_table_get_string(sti, x1);
	}
	x2 = dset->Z[cv][t];
	if (tab->cstrs) {
	    s2 = series_table_get_string(stj, x2);
	}
	for (i=0; i<tab->rows; i++) {
	    imatch = xtab_row_match(tab, i, s1, x1);
	    if (imatch) {
		jmatch = 0;
		for (j=0; j<tab->cols && !jmatch; j++) {
		    jmatch = xtab_col_match(tab, j, s2, x2);
		    if (jmatch) {
			tab->f[i][j] += 1;
			tab->rtotal[i] += 1;
			tab->ctotal[j] += 1;
		    }
		}
		break;
	    }
	}
    }

 bailout:

    if (*err) {
	free_xtab(tab);
	tab = NULL;
    }

    return tab;
}

static void record_xtab (const Xtab *tab, const DATASET *dset,
			 gretlopt opt)
{
    gretl_matrix *X = NULL;
    char **Sc = NULL, **Sr = NULL;
    double xij, cj, ri;
    int totals, rows, cols;
    int i, j;

    totals = (opt & OPT_N)? 0 : 1;
    rows = tab->rows + totals;
    cols = tab->cols + totals;

    X = gretl_zero_matrix_new(rows, cols);
    if (X == NULL) {
	return;
    }

    /* column labels */
    Sc = strings_array_new(cols);
    if (Sc != NULL) {
	for (j=0; j<tab->cols; j++) {
	    if (tab->cstrs) {
		Sc[j] = gretl_strdup(tab->Sc[j]);
	    } else {
		cj = tab->cval[j];
		Sc[j] = gretl_strdup_printf("%4g", cj);
	    }
	}
	if (totals) {
	    Sc[cols-1] = gretl_strdup("TOTAL");
	}
    }

    /* row labels */
    Sr = strings_array_new(rows);
    if (Sr != NULL) {
	for (i=0; i<tab->rows; i++) {
	    if (tab->rstrs) {
		Sr[i] = gretl_strdup(tab->Sr[i]);
	    } else {
		ri = tab->rval[i];
		Sr[i] = gretl_strdup_printf("%4g", ri);
	    }
	}
	if (totals) {
	    Sr[rows-1] = gretl_strdup("TOTAL");
	}
    }

    /* body of table */

    for (i=0; i<tab->rows; i++) {
	if (tab->rtotal[i] > 0) {
	    /* row counts */
	    for (j=0; j<tab->cols; j++) {
		if (tab->ctotal[j] > 0) {
		    if (opt & (OPT_C | OPT_R)) {
			if (opt & OPT_C) {
			    xij = 100.0 * tab->f[i][j] / tab->ctotal[j];
			} else {
			    xij = 100.0 * tab->f[i][j] / tab->rtotal[i];
			}
		    } else {
			xij = tab->f[i][j];
		    }
		    gretl_matrix_set(X, i, j, xij);
		}
	    }
	    if (totals) {
		/* row totals */
		if (opt & OPT_C) {
		    xij = 100.0 * tab->rtotal[i] / tab->n;
		} else {
		    xij = tab->rtotal[i];
		}
		gretl_matrix_set(X, i, cols-1, xij);
	    }
	}
    }

    if (totals) {
	/* column totals */
	for (j=0; j<tab->cols; j++) {
	    if (opt & OPT_R) {
		xij = 100.0 * tab->ctotal[j] / tab->n;
	    } else {
		xij = tab->ctotal[j];
	    }
	    gretl_matrix_set(X, rows-1, j, xij);
	}
	gretl_matrix_set(X, rows-1, cols-1, tab->n);
    }

    gretl_matrix_set_colnames(X, Sc);
    gretl_matrix_set_rownames(X, Sr);
    set_last_result_data(X, GRETL_TYPE_MATRIX);
}

/* For use in the context of "xtab" with --quiet option:
   compute and record the Pearson chi-square value and its
   p-value.
*/

static int xtab_do_pearson (const Xtab *tab)
{
    double x, y, ymin = 1.0e-7;
    double pearson = 0.0;
    int i, j, err = 0;

    for (i=0; i<tab->rows && !err; i++) {
	if (tab->rtotal[i] > 0) {
	    for (j=0; j<tab->cols && !err; j++) {
		y = ((double) tab->rtotal[i] * tab->ctotal[j]) / tab->n;
		if (y < ymin) {
		    err = E_DATA;
		} else {
		    x = (double) tab->f[i][j] - y;
		    pearson += x * x / y;
		}
	    }
	}
    }

    if (!err) {
	int df = (tab->rows - 1) * (tab->cols - 1);
	double pval = chisq_cdf_comp(df, pearson);

	if (!na(pval)) {
	    record_test_result(pearson, pval);
	} else {
	    err = E_DATA;
	}
    }

    if (err) {
	record_test_result(NADBL, NADBL);
    }

    return err;
}

int crosstab_from_matrix (gretlopt opt, PRN *prn)
{
    const char *mname;
    gretl_matrix *m;
    Xtab *tab = NULL;
    int i, j, nvals, n = 0;
    int free_m = 0;
    double x;
    int err = 0;

    mname = get_optval_string(XTAB, OPT_X);
    if (mname == NULL) {
	return E_DATA;
    }

    m = get_matrix_by_name(mname);
    if (m == NULL) {
        m = generate_matrix(mname, NULL, &err);
        if (m == NULL) {
            return E_UNKVAR;
        } else {
            free_m = 1;
        }
    }

    if (m->rows < 2 || m->cols < 2) {
	err = E_DATA;
    }

    nvals = m->rows * m->cols;

    for (i=0; i<nvals && !err; i++) {
	x = m->val[i];
	if (x < 0 || x != floor(x) || x > INT_MAX) {
	    err = E_DATA;
	}
	n += x;
    }

    if (err) {
	gretl_errmsg_sprintf(_("Matrix %s does not represent a "
			       "contingency table"), mname);
    } else {
        tab = xtab_new(n, 0, 0);
        if (tab == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        goto bailout;
    }

    tab->rows = m->rows;
    tab->cols = m->cols;
    tab->Sc = (char **) gretl_matrix_get_colnames(m);
    tab->Sr = (char **) gretl_matrix_get_rownames(m);

    if (xtab_allocate_arrays(tab)) {
        err = E_ALLOC;
        goto bailout;
    }

    for (i=0; i<m->rows; i++) {
	tab->rval[i] = i + 1;
	tab->rtotal[i] = 0.0;
	for (j=0; j<m->cols; j++) {
	    tab->f[i][j] = gretl_matrix_get(m, i, j);
	    tab->rtotal[i] += tab->f[i][j];
	}
    }

    for (j=0; j<m->cols; j++) {
	tab->cval[j] = j + 1;
	tab->ctotal[j] = 0.0;
	for (i=0; i<m->rows; i++) {
	    tab->ctotal[j] += tab->f[i][j];
	}
    }

    if (opt & OPT_Q) {
	xtab_do_pearson(tab);
    } else {
	print_xtab(tab, NULL, opt | OPT_S, prn);
    }

 bailout:

    if (tab != NULL) {
        tab->Sc = NULL;
        tab->Sr = NULL;
        free_xtab(tab);
    }
    if (free_m) {
        gretl_matrix_free(m);
    }

    return err;
}

int crosstab (const int *list, const DATASET *dset,
	      gretlopt opt, PRN *prn)
{
    Xtab *tab;
    int *rowvar = NULL;
    int *colvar = NULL;
    int i, j, vi, vj, k;
    int pos, onelist;
    int nrv, ncv;
    int err = 0;

    pos = gretl_list_separator_position(list);
    onelist = (pos == 0);

    if (pos == 0) {
	/* single list case */
	nrv = list[0];
	ncv = nrv - 1;
    } else {
	/* double list case */
	nrv = pos - 1;
	ncv = list[0] - pos;
    }

    if (nrv == 0 || ncv == 0) {
	return E_PARSE;
    }

    rowvar = gretl_list_new(nrv);
    if (rowvar == NULL) {
	return E_ALLOC;
    }

    j = 1;
    for (i=1; i<=nrv; i++) {
	k = list[i];
	if (accept_as_discrete(dset, k, 0)) {
	    rowvar[j++] = k;
	} else {
	    pprintf(prn, _("dropping %s: not a discrete variable\n"),
		    dset->varname[k]);
	    rowvar[0] -= 1;
	}
    }

    if (rowvar[0] == 0 || (onelist && rowvar[0] == 1)) {
	gretl_errmsg_set("xtab: variables must be discrete");
	free(rowvar);
	return E_TYPES;
    }

    if (onelist && rowvar[0] == 2) {
	/* the bivariate case */
	tab = get_new_xtab(rowvar[1], rowvar[2], dset, &err);
	if (!err) {
	    /* make $result matrix available */
	    record_xtab(tab, dset, opt);
	    if (opt & OPT_Q) {
		/* quiet: run and record the Pearson test */
		xtab_do_pearson(tab);
	    } else {
		/* print, and record Pearson test */
		print_xtab(tab, dset, opt | OPT_S, prn);
	    }
	    free_xtab(tab);
	}
	goto finish;
    }

    if (!onelist) {
	/* construct the second list */
	colvar = gretl_list_new(ncv);
	if (colvar == NULL) {
	    err = E_ALLOC;
	} else {
	    j = 1;
	    for (i=1; i<=ncv; i++) {
		k = list[pos+i];
		if (accept_as_discrete(dset, k, 0)) {
		    colvar[j++] = k;
		} else {
		    colvar[0] -= 1;
		}
	    }
	    if (colvar[0] == 0) {
		err = E_TYPES;
	    }
	}
    }

    for (i=1; i<=rowvar[0] && !err; i++) {
	vi = rowvar[i];
	if (onelist) {
	    /* single list case */
	    for (j=1; j<i && !err; j++) {
		vj = rowvar[j];
		tab = get_new_xtab(vj, vi, dset, &err);
		if (!err) {
		    print_xtab(tab, dset, opt, prn);
		    free_xtab(tab);
		}
	    }
	} else {
	    /* double list case */
	    for (j=1; j<=colvar[0] && !err; j++) {
		vj = colvar[j];
		tab = get_new_xtab(vi, vj, dset, &err);
		if (!err) {
		    print_xtab(tab, dset, opt, prn);
		    free_xtab(tab);
		}
	    }
	}
    }

 finish:

    free(rowvar);
    free(colvar);

    return err;
}

Xtab *single_crosstab (const int *list, const DATASET *dset,
		       gretlopt opt, PRN *prn, int *err)
{
    Xtab *tab = NULL;
    int rv, cv;

    if (list[0] != 2) {
	*err = E_DATA;
	return NULL;
    }

    rv = list[1];
    cv = list[2];

    if (accept_as_discrete(dset, rv, 0) &&
	accept_as_discrete(dset, cv, 0)) {
	tab = get_new_xtab(rv, cv, dset, err);
    } else {
	*err = E_TYPES;
    }

    if (!*err) {
	print_xtab(tab, dset, opt, prn);
    }

    return tab;
}

static int just_record_freq_test (const FreqDist *freq)
{
    double pval = NADBL;

    if (freq->dist == D_NORMAL) {
	pval = chisq_cdf_comp(2, freq->test);
    } else if (freq->dist == D_GAMMA) {
	pval = normal_pvalue_2(freq->test);
    }

    if (na(pval)) {
	return E_NAN;
    } else {
	record_test_result(freq->test, pval);
	return 0;
    }
}

int model_error_dist (const MODEL *pmod, DATASET *dset,
		      gretlopt opt, PRN *prn)
{
    FreqDist *freq = NULL;
    gretlopt fopt = OPT_Z; /* show normality test */
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int err = 0;

    if (pmod == NULL || pmod->uhat == NULL) {
	return E_DATA;
    }

    err = gretl_model_get_normality_test(pmod, prn);

    if (!err) {
	return 0;
    } else if (LIMDEP(pmod->ci)) {
	return err;
    } else {
	err = 0;
    }

    if (exact_fit_check(pmod, prn)) {
	return 0;
    }

    if (genr_fit_resid(pmod, dset, M_UHAT)) {
	return E_ALLOC;
    }

    if (!err) {
	dset->t1 = pmod->t1;
	dset->t2 = pmod->t2;
	freq = get_freq(dset->v - 1, dset, NADBL, NADBL, 0,
			pmod->ncoeff, fopt, &err);
    }

    if (!err) {
	if (opt & OPT_I) {
	    err = just_record_freq_test(freq);
	} else if (opt & OPT_Q) {
	    print_freq_test(freq, prn);
	} else {
	    print_freq(freq, 0, NULL, prn);
	}
	free_freq(freq);
    }

    dataset_drop_last_variables(dset, 1);

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return err;
}

/* Compute PACF via Durbin-Levinson algorithm and write
   the values into @A following the ACF values already
   present.
*/

static int get_pacf (gretl_matrix *A)
{
    int m = gretl_matrix_rows(A);
    gretl_matrix *phi;
    double *acf = A->val;
    double *pacf = acf + m;
    double x, num, den;
    int i, j;

    phi = gretl_matrix_alloc(m, m);
    if (phi == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_set(A, 0, 1, acf[0]);
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

/* also used in GUI, library.c */

int auto_acf_order (int T)
{
    int p = 10 * log10(T);

    if (p > T / 5) {
	/* restrict to 20 percent of data (Tadeusz) */
	p = T / 5;
    }

    return p;
}

/**
 * gretl_acf:
 * @k: lag order.
 * @t1: starting observation.
 * @t2: ending observation.
 * @y: data series.
 * @ybar: mean of @y over range @t1 to @t2.
 *
 * Returns: the autocorrelation at lag @k for the series @y over
 * the range @t1 to @t2, or #NADBL on failure.
 */

static double gretl_acf (int k, int t1, int t2, const double *y,
			 double ybar)
{
    double z, num = 0, den = 0;
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
	if (na(y[t])) {
	    return NADBL;
	}
	z = y[t] - ybar;
	den += z * z;
	if (t - k >= t1) {
	    if (na(y[t-k])) {
		return NADBL;
	    }
	    num += z * (y[t-k] - ybar);
	    n++;
	}
    }

    if (n == 0) {
	return NADBL;
    }

    return num / den;
}

/**
 * ljung_box:
 * @m: maximum lag.
 * @t1: starting observation.
 * @t2: ending observation.
 * @y: data series.
 * @err: location to receive error code.
 *
 * Returns: the Ljung-Box statistic for lag order @m for
 * the series @y over the sample @t1 to @t2, or #NADBL
 * on failure.
 */

double ljung_box (int m, int t1, int t2, const double *y, int *err)
{
    double acf, ybar = 0.0, LB = 0.0;
    int k, n = t2 - t1 + 1;

    *err = 0;

    if (n == 0 || gretl_isconst(t1, t2, y)) {
	*err = E_DATA;
    } else if (m <= 0) {
	gretl_errmsg_sprintf(_("Invalid lag order %d"), m);
	*err = E_DATA;
    } else {
	ybar = gretl_mean(t1, t2, y);
	if (na(ybar)) {
	    *err = E_DATA;
	}
    }

    if (*err) {
	return NADBL;
    }

    /* calculate acf up to lag m, cumulating LB */
    for (k=1; k<=m; k++) {
	acf = gretl_acf(k, t1, t2, y, ybar);
	if (na(acf)) {
	    *err = E_MISSDATA;
	    break;
	}
	LB += acf * acf / (n - k);
    }

    if (*err) {
	LB = NADBL;
    } else {
	LB *= n * (n + 2.0);
    }

    return LB;
}

/**
 * gretl_xcf:
 * @k: lag order (or lead order if < 0).
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: first data series.
 * @y: second data series.
 * @xbar: mean of first series.
 * @ybar: mean of second series.
 *
 * Returns: the cross-correlation at lag (or lead) @k for the
 * series @x and @y over the range @t1 to @t2, or #NADBL on failure.
 */

static double
gretl_xcf (int k, int t1, int t2, const double *x, const double *y,
	   double xbar, double ybar)
{
    double num = 0, den1 = 0, den2 = 0;
    double zx, zy;
    int t;

    for (t=t1; t<=t2; t++) {
	if (na(x[t]) || na(y[t])) {
	    return NADBL;
	}
	zx = x[t] - xbar;
	zy = y[t] - ybar;
	den1 += zx * zx;
	den2 += zy * zy;
	if ((k >= 0 && t - k >= t1) || (k < 0 && t - k <= t2)) {
	    num += zx * (y[t-k] - ybar);
	}
    }

    return num / sqrt(den1 * den2);
}

static int corrgm_plot_wanted (PlotType ptype, gretlopt opt, int *err)
{
    if (opt & OPT_Q) {
	/* handle legacy option: --quiet = no plot */
	return 0;
    } else if (opt & OPT_U) {
	/* ignore legacy option: --plot=ascii */
	int ci = ptype == PLOT_CORRELOGRAM ? CORRGM : XCORRGM;
	const char *s = get_optval_string(ci, OPT_U);

	if (s != NULL && !strcmp(s, "ascii")) {
	    return 0;
	}
    }

    return gnuplot_graph_wanted(ptype, opt, err);
}

/* Full Bartlett calculation for three z-values, called from
   corrgm_plot() if we're printing ACF output. BPM will have three
   columns.
*/

static void do_acf_bartlett (gretl_matrix *BPM, int i,
			     int T, double ssr,
			     const double *z)
{
    double c = sqrt((1.0/T) * (1 + 2*ssr));
    int j;

    for (j=0; j<BPM->cols; j++) {
	gretl_matrix_set(BPM, i, j, z[j] * c);
    }
}

/* Basic Bartlett if we're not printing ACF output. BPM will have a
   single column.
*/

static void corrgm_basic_bartlett (const double *acf,
				   gretl_matrix *BPM,
				   int T)
{
    double Tm1 = 1.0 / T;
    double z = 1.96;
    double ssr = 0;
    int k;

    for (k=0; k<BPM->rows; k++) {
	if (na(acf[k])) {
	    continue;
	}
	if (k > 0 && !na(acf[k-1])) {
	    ssr += acf[k-1] * acf[k-1];
	}
	BPM->val[k] = z * sqrt(Tm1 * (1 + 2*ssr));
    }
}

static void corrgm_do_ljung_box (const double *acf,
				 int m, int T,
				 int nparam)
{
    double lbox = 0.0;
    int k, dfQ = 0;

    for (k=0; k<m; k++) {
	if (!na(acf[k])) {
	    lbox += (T * (T + 2.0)) * acf[k] * acf[k] / (T - (k + 1));
	    dfQ++;
	}
    }

    dfQ -= nparam;
    if (lbox > 0 && dfQ > 0) {
	double pval = chisq_cdf_comp(dfQ, lbox);

	if (!na(pval)) {
	    record_test_result(lbox, pval);
	}
    }
}

/* The following is called if the --silent option is not given to
   "corrgm". Note that it embeds Bartlett calculations for three
   confidence levels. We make do with a simpler Bartlett variant if
   we're not printing the ACF table.
*/

static void corrgm_print (const char *vname,
			  const double *acf,
			  const double *pacf,
			  gretl_matrix *BPM,
			  int m, int T,
			  int nparam,
			  gretlopt opt,
			  PRN *prn)
{
    const double z[] = {1.65, 1.96, 2.58};
    double pm[3] = {0};
    double pval = NADBL;
    double lbox = 0.0;
    double ssr = 0.0;
    int i, k, dfQ = 1;

    if (opt & OPT_R) {
	pprintf(prn, "\n%s\n", _("Residual autocorrelation function"));
    } else {
	pputc(prn, '\n');
	pprintf(prn, _("Autocorrelation function for %s"), vname);
	pputc(prn, '\n');
    }

    pputs(prn, _("***, **, * indicate significance at the 1%, 5%, 10% levels\n"));

    if (opt & OPT_B) {
	pputs(prn, _("using Bartlett standard errors for ACF"));
    } else {
	pprintf(prn, _("using standard error 1/T^%.1f"), 0.5);
    }

    pputs(prn, "\n\n");
    if (pacf != NULL) {
	pputs(prn, _("  LAG      ACF          PACF         Q-stat. [p-value]"));
    } else {
	pputs(prn, _("  LAG      ACF          Q-stat. [p-value]"));
    }
    pputs(prn, "\n\n");

    if (BPM == NULL || pacf != NULL) {
	/* if we're not doing Bartlett confidence bands, or we're
	   showing PACF as well as ACF
	*/
	for (i=0; i<3; i++) {
	    pm[i] = z[i] / sqrt((double) T);
	    if (pm[i] > 0.5) {
		pm[i] = 0.5;
	    }
	}
    }

    for (k=0; k<m; k++) {
	double pm0, pm1, pm2;

	if (na(acf[k])) {
	    pprintf(prn, "%5d\n", k + 1);
	    continue;
	}

	/* ACF */
	pprintf(prn, "%5d%9.4f ", k + 1, acf[k]);
	if (BPM != NULL) {
	    if (k > 0 && !na(acf[k-1])) {
		ssr += acf[k-1] * acf[k-1];
	    }
	    do_acf_bartlett(BPM, k, T, ssr, z);
	    pm0 = gretl_matrix_get(BPM, k, 0);
	    pm1 = gretl_matrix_get(BPM, k, 1);
	    pm2 = gretl_matrix_get(BPM, k, 2);
	} else {
	    pm0 = pm[0];
	    pm1 = pm[1];
	    pm2 = pm[2];
	}
	if (fabs(acf[k]) > pm2) {
	    pputs(prn, " ***");
	} else if (fabs(acf[k]) > pm1) {
	    pputs(prn, " ** ");
	} else if (fabs(acf[k]) > pm0) {
	    pputs(prn, " *  ");
	} else {
	    pputs(prn, "    ");
	}

	/* PACF, if present */
	if (pacf != NULL) {
	    if (na(pacf[k])) {
		bufspace(13, prn);
	    } else {
		/* PACF */
		pprintf(prn, "%9.4f", pacf[k]);
		if (fabs(pacf[k]) > pm[2]) {
		    pputs(prn, " ***");
		} else if (fabs(pacf[k]) > pm[1]) {
		    pputs(prn, " ** ");
		} else if (fabs(pacf[k]) > pm[0]) {
		    pputs(prn, " *  ");
		} else {
		    pputs(prn, "    ");
		}
	    }
	}

	/* Ljung-Box Q */
	lbox += (T * (T + 2.0)) * acf[k] * acf[k] / (T - (k + 1));

	if (k >= nparam) {
	    /* i.e., if the real df is > 0 */
	    pprintf(prn, "%12.4f", lbox);
	    pval = chisq_cdf_comp(dfQ++, lbox);
	    pprintf(prn, "  [%5.3f]", pval);
	}

	pputc(prn, '\n');
    }
    pputc(prn, '\n');

    if (lbox > 0 && !na(pval)) {
	record_test_result(lbox, pval);
    }
}

/**
 * corrgram:
 * @varno: ID number of variable to process.
 * @order: integer order for autocorrelation function.
 * @nparam: number of estimated parameters (e.g. for the
 * case of ARMA), used to correct the degrees of freedom
 * for Q test.
 * @dset: dataset struct.
 * @opt: if includes OPT_R, variable in question is a model
 * residual generated "on the fly"; OPT_U can be used to
 * specify a plot option; OPT_S says not to print anything;
 * OPT_A means ACF only (no PACF); OPT_B means Bartlett.
 * @prn: gretl printing struct.
 *
 * Computes the autocorrelation function and plots the correlogram for
 * the variable specified by @varno.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int corrgram (int varno, int order, int nparam, DATASET *dset,
	      gretlopt opt, PRN *prn)
{
    gretl_matrix *BPM = NULL;
    gretl_matrix *A = NULL;
    double ybar;
    double *acf = NULL;
    double *pacf = NULL;
    const char *vname;
    int BPMcols = 3;
    int k, m, T;
    int do_plot = 0;
    int do_print = 1;
    int do_pacf = 1;
    int t1 = dset->t1;
    int t2 = dset->t2;
    int err = 0;

    gretl_error_clear();

    if (order < 0) {
	gretl_errmsg_sprintf(_("Invalid lag order %d"), order);
	return E_DATA;
    }

    err = series_adjust_sample(dset->Z[varno], &t1, &t2);
    if (err) {
	return err;
    }

    if ((T = t2 - t1 + 1) < 4) {
	return E_TOOFEW;
    }

    if (gretl_isconst(t1, t2, dset->Z[varno])) {
	gretl_errmsg_sprintf(_("%s is a constant"), dset->varname[varno]);
	return E_DATA;
    }

    ybar = gretl_mean(t1, t2, dset->Z[varno]);
    if (na(ybar)) {
	return E_DATA;
    }

    if (opt & OPT_S) {
	do_print = 0;
	BPMcols = 1;
    }
    if (opt & OPT_A) {
	do_pacf = 0;
    }

    /* lag order for acf */
    m = order;
    if (m == 0) {
	m = auto_acf_order(T);
    } else if (m > T - dset->pd) {
	int mmax = T - 1;

	if (m > mmax) {
	    m = mmax;
	}
    }

    /* allocate space for ACF, and perhaps PACF */
    A = gretl_matrix_alloc(m, do_pacf ? 2 : 1);
    if (A == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (opt & OPT_B) {
	/* Bartlett */
	BPM = gretl_matrix_alloc(m, BPMcols);
	if (BPM == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    /* graphing? */
    do_plot = corrgm_plot_wanted(PLOT_CORRELOGRAM, opt, &err);
    if (err) {
	goto bailout;
    }

    /* convenience pointer into @A */
    acf = A->val;

    /* calculate acf up to order @m */
    for (k=0; k<m; k++) {
	acf[k] = gretl_acf(k+1, t1, t2, dset->Z[varno], ybar);
    }

    /* try adding pacf into @A if wanted */
    if (do_pacf) {
	int pacf_err = get_pacf(A);

	pacf = pacf_err ? NULL : acf + m;
    }

    vname = plotname(dset, varno, 1);

    if (do_print) {
	corrgm_print(vname, acf, pacf, BPM, m, T, nparam, opt, prn);
    } else {
	if (opt & OPT_B) {
	    corrgm_basic_bartlett(acf, BPM, T);
	}
	corrgm_do_ljung_box(acf, m, T, nparam);
    }

    if (do_plot) {
	double pm = 1.96 / sqrt((double) T);

	if (pm > 0.5) pm = 0.5;
	err = correlogram_plot(vname, acf, pacf, BPM, m, pm, opt);
    }

 bailout:

    gretl_matrix_free(A);
    gretl_matrix_free(BPM);

    return err;
}

/**
 * acf_matrix:
 * @x: series to analyse.
 * @order: maximum lag for autocorrelation function.
 * @dset: information on the data set, or %NULL.
 * @n: length of series (required if @dset is %NULL).
 * @err: location to receive error code.
 *
 * Computes the autocorrelation function for series @x with
 * maximum lag @order.
 *
 * Returns: two-column matrix containing the values of the
 * ACF and PACF at the successive lags, or %NULL on error.
 */

gretl_matrix *acf_matrix (const double *x, int order,
			  const DATASET *dset, int n,
			  int *err)
{
    gretl_matrix *A = NULL;
    double xbar;
    int m, k, t, T;
    int t1, t2;

    if (dset != NULL) {
	t1 = dset->t1;
	t2 = dset->t2;

	while (na(x[t1])) t1++;
	while (na(x[t2])) t2--;

	T = t2 - t1 + 1;
    } else {
	t1 = 0;
	t2 = n - 1;
	T = n;
    }

    if (T < 4) {
	*err = E_TOOFEW;
	return NULL;
    }

    for (t=t1; t<=t2; t++) {
	if (na(x[t])) {
	    *err = E_MISSDATA;
	    return NULL;
	}
    }

    if (gretl_isconst(t1, t2, x)) {
	gretl_errmsg_set(_("Argument is a constant"));
	*err = E_DATA;
	return NULL;
    }

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
	*err = E_DATA;
	return NULL;
    }

    m = order;

    if (dset == NULL) {
	if (m < 1 || m > T) {
	    *err = E_DATA;
	    return NULL;
	}
    } else {
	if (m == 0) {
	    m = auto_acf_order(T);
	} else if (m > T - dset->pd) {
	    int mmax = T - 1;

	    if (m > mmax) {
		m = mmax;
	    }
	}
    }

    A = gretl_matrix_alloc(m, 2);
    if (A == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* calculate ACF up to order m */
    for (k=0; k<m && !*err; k++) {
	A->val[k] = gretl_acf(k+1, t1, t2, x, xbar);
	if (na(A->val[k])) {
	    *err = E_DATA;
	}
    }

    /* add PACF */
    if (!*err) {
	*err = get_pacf(A);
    }

    if (*err) {
	gretl_matrix_free(A);
	A = NULL;
    }

    return A;
}

static int xcorrgm_graph (const char *xname, const char *yname,
			  double *xcf, int m, double *pm,
			  int allpos)
{
    char crit_string[16];
    gchar *title;
    FILE *fp;
    int k, err = 0;

    fp = open_plot_input_file(PLOT_XCORRELOGRAM, 0, &err);
    if (err) {
	return err;
    }

    sprintf(crit_string, "%.2f/T^%.1f", 1.96, 0.5);

    gretl_push_c_numeric_locale();

    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);
    print_keypos_string(GP_KEY_RIGHT_TOP, fp);
    fprintf(fp, "set xlabel '%s'\n", _("lag"));
    if (allpos) {
	fputs("set yrange [-0.1:1.1]\n", fp);
    } else {
	fputs("set yrange [-1.1:1.1]\n", fp);
    }
    title = g_strdup_printf(_("Correlations of %s and lagged %s"),
			    xname, yname);
    fprintf(fp, "set title '%s'\n", title);
    g_free(title);
    fprintf(fp, "set xrange [%d:%d]\n", -(m + 1), m + 1);
    if (allpos) {
	fprintf(fp, "plot \\\n"
		"'-' using 1:2 notitle w impulses lw 5, \\\n"
		"%g title '%s' lt 2\n", pm[1], crit_string);
    } else {
	fprintf(fp, "plot \\\n"
		"'-' using 1:2 notitle w impulses lw 5, \\\n"
		"%g title '+- %s' lt 2, \\\n"
		"%g notitle lt 2\n", pm[1], crit_string, -pm[1]);
    }

    for (k=-m; k<=m; k++) {
	fprintf(fp, "%d %g\n", k, xcf[k+m]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

/* We assume here that all data issues have already been
   assessed (lag length, missing values etc.) and we just
   get on with the job.
*/

static gretl_matrix *real_xcf_vec (const double *x, const double *y,
				   int p, int T, int *err)
{
    gretl_matrix *xcf;
    double xbar, ybar;
    int i;

    xbar = gretl_mean(0, T-1, x);
    if (na(xbar)) {
	*err = E_DATA;
	return NULL;
    }

    ybar = gretl_mean(0, T-1, y);
    if (na(ybar)) {
	*err = E_DATA;
	return NULL;
    }

    xcf = gretl_column_vector_alloc(p * 2 + 1);
    if (xcf == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=-p; i<=p; i++) {
	xcf->val[i+p] = gretl_xcf(i, 0, T - 1, x, y, xbar, ybar);
    }

    return xcf;
}

/* for arrays @x and @y, check that there are no missing values
   and that neither series has a constant value
*/

static int xcf_data_check (const double *x, const double *y, int T,
			   int *badvar)
{
    int xconst = 1, yconst = 1;
    int t;

    if (T < 5) {
	return E_TOOFEW;
    }

    for (t=0; t<T; t++) {
	if (na(x[t]) || na(y[t])) {
	    return E_MISSDATA;
	}
	if (t > 0 && x[t] != x[0]) {
	    xconst = 0;
	}
	if (t > 0 && y[t] != y[0]) {
	    yconst = 0;
	}
    }

    if (xconst) {
	*badvar = 1;
	return E_DATA;
    } else if (yconst) {
	*badvar = 2;
	return E_DATA;
    }

    return 0;
}

/**
 * xcf_vec:
 * @x: first series.
 * @y: second series.
 * @p: maximum lag for cross-correlation function.
 * @dset: information on the data set, or NULL.
 * @n: length of series (required only if @dset is NULL).
 * @err: location to receive error code.
 *
 * Computes the cross-correlation function for series @x with
 * series @y up to maximum lag @order.
 *
 * Returns: column vector containing the values of the
 * cross-correlation function, or NULL on error.
 */

gretl_matrix *xcf_vec (const double *x, const double *y,
		       int p, const DATASET *dset,
		       int n, int *err)
{
    gretl_matrix *xcf = NULL;
    int t1, t2;
    int T, badvar = 0;

    if (p <= 0) {
	*err = E_DATA;
	return NULL;
    }

    if (dset != NULL) {
	int yt1, yt2;

	t1 = yt1 = dset->t1;
	t2 = yt2 = dset->t2;

	while (na(x[t1])) t1++;
	while (na(y[yt1])) yt1++;

	while (na(x[t2])) t2--;
	while (na(y[yt2])) yt2--;

	t1 = (yt1 > t1)? yt1 : t1;
	t2 = (yt2 < t2)? yt2 : t2;

	T = t2 - t1 + 1;
    } else {
	t1 = 0;
	t2 = n - 1;
	T = n;
    }

#if 0
    fprintf(stderr, "t1=%d, t2=%d, T=%d\n", t1, t2, T);
    fprintf(stderr, "x[t1]=%g, y[t1]=%g\n\n", x[t1], y[t1]);
#endif

    if (dset != NULL) {
	if (2 * p > T - dset->pd) { /* ?? */
	    *err = E_DATA;
	}
    } else if (2 * p > T) {
	*err = E_DATA;
    }

    if (!*err) {
	*err = xcf_data_check(x + t1, y + t1, T, &badvar);
	if (badvar) {
	    gretl_errmsg_sprintf(_("Argument %d is a constant"),
				 badvar);
	}
    }

    if (!*err) {
	xcf = real_xcf_vec(x + t1, y + t1, p, T, err);
    }

    return xcf;
}

/**
 * xcorrgram:
 * @list: should contain ID numbers of two variables.
 * @order: integer order for autocorrelation function.
 * @dset: dataset struct.
 * @opt: may include OPT_U for plot options.
 * @prn: gretl printing struct.
 *
 * Computes the cross-correlation function and plots the
 * cross-correlogram for the specified variables.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int xcorrgram (const int *list, int order, DATASET *dset,
	       gretlopt opt, PRN *prn)
{
    gretl_matrix *xcf = NULL;
    double pm[3];
    const char *xname, *yname;
    const double *x, *y;
    int t1 = dset->t1, t2 = dset->t2;
    int do_plot = 0;
    int k, p, badvar = 0;
    int T, err = 0;

    gretl_error_clear();

    if (order < 0) {
	gretl_errmsg_sprintf(_("Invalid lag order %d"), order);
	return E_DATA;
    }

    if (list[0] != 2) {
	return E_DATA;
    }

    x = dset->Z[list[1]];
    y = dset->Z[list[2]];

    err = list_adjust_sample(list, &t1, &t2, dset, NULL);

    if (!err) {
	T = t2 - t1 + 1;
	err = xcf_data_check(x + t1, y + t1, T, &badvar);
    }

    if (err) {
	if (badvar) {
	    gretl_errmsg_sprintf(_("%s is a constant"),
				 dset->varname[list[badvar]]);
	}
	return err;
    }

    /* graphing? */
    do_plot = corrgm_plot_wanted(PLOT_XCORRELOGRAM, opt, &err);
    if (err) {
	return err;
    }

    xname = dset->varname[list[1]];
    yname = dset->varname[list[2]];

    p = order;
    if (p == 0) {
	p = auto_acf_order(T) / 2;
    } else if (2 * p > T - dset->pd) {
	p = (T - 1) / 2; /* ?? */
    }

    xcf = real_xcf_vec(x + t1, y + t1, p, T, &err);
    if (err) {
	return err;
    }

    /* for confidence bands */
    pm[0] = 1.65 / sqrt((double) T);
    pm[1] = 1.96 / sqrt((double) T);
    pm[2] = 2.58 / sqrt((double) T);

    pputc(prn, '\n');
    pprintf(prn, _("Cross-correlation function for %s and %s"),
	    xname, yname);
    pputs(prn, "\n\n");
    pputs(prn, _("  LAG      XCF"));
    pputs(prn, "\n\n");

    for (k=-p; k<=p; k++) {
	double x = xcf->val[k + p];

	pprintf(prn, "%5d%9.4f", k, x);
	if (fabs(x) > pm[2]) {
	    pputs(prn, " ***");
	} else if (fabs(x) > pm[1]) {
	    pputs(prn, " **");
	} else if (fabs(x) > pm[0]) {
	    pputs(prn, " *");
	}
	pputc(prn, '\n');
    }
    pputc(prn, '\n');

    if (do_plot) {
	int allpos = 1;

	for (k=-p; k<=p; k++) {
	    if (xcf->val[k+p] < 0) {
		allpos = 0;
		break;
	    }
	}
	err = xcorrgm_graph(xname, yname, xcf->val, p, pm, allpos);
    }

    gretl_matrix_free(xcf);

    return err;
}

struct fractint_test {
    double d;    /* estimated degree of integration */
    double se;   /* standard error of the above */
};

static int fract_int_GPH (int T, int m, const double *dens,
			  struct fractint_test *ft)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *V = NULL;
    double x, w;
    int t, err = 0;

    ft->d = ft->se = NADBL;

    y = gretl_column_vector_alloc(m);
    X = gretl_unit_matrix_new(m, 2);
    b = gretl_column_vector_alloc(2);
    V = gretl_matrix_alloc(2, 2);

    if (y == NULL || X == NULL || b == NULL || V == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* Test from Geweke and Porter-Hudak, as set out in
       Greene, Econometric Analysis 4e, p. 787 */

    for (t=0; t<m; t++) {
	y->val[t] = log(dens[t]);
	w = M_2PI * (t + 1) / (double) T;
	x = sin(w / 2);
	gretl_matrix_set(X, t, 1, log(4 * x * x));
    }

    err = gretl_matrix_ols(y, X, b, V, NULL, &x);

    if (!err) {
	ft->d = -b->val[1];
	ft->se = sqrt(gretl_matrix_get(V, 1, 1));
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(V);

    return err;
}

static gretl_matrix *gretl_matrix_pergm (const gretl_matrix *x, int m,
					 int *err)
{
    gretl_matrix *p = NULL;
    gretl_matrix *f = NULL;

    f = gretl_matrix_fft(x, err);
    if (*err) {
	return NULL;
    }

    p = gretl_column_vector_alloc(m);

    if (p == NULL) {
	*err = E_ALLOC;
    } else {
	int T = gretl_vector_get_length(x);
	double re, im, scale = M_2PI * T;
	double complex z;
	int i;

	for (i=0; i<m; i++) {
	    z = gretl_cmatrix_get(f, i+1, 0);
	    re = creal(z);
	    im = cimag(z);
	    p->val[i] = (re*re + im*im) / scale;
	}
    }

    gretl_matrix_free(f);

    return p;
}

struct LWE_helper {
    gretl_matrix *lambda;
    gretl_matrix *lpow;
    gretl_matrix *I1;
    gretl_matrix *I2;
    double lcm;
};

static void LWE_free (struct LWE_helper *L)
{
    gretl_matrix_free(L->lambda);
    gretl_matrix_free(L->lpow);
    gretl_matrix_free(L->I1);
    gretl_matrix_free(L->I2);
}

static double LWE_obj_func (struct LWE_helper *L, double d)
{
    double dd = 2.0 * d;
    int i;

    gretl_matrix_copy_values(L->lpow, L->lambda);
    gretl_matrix_raise(L->lpow, dd);

    for (i=0; i<L->I1->rows; i++) {
	L->I2->val[i] = L->I1->val[i] * L->lpow->val[i];
    }

    return -(log(gretl_vector_mean(L->I2)) - dd * L->lcm);
}

static gretl_matrix *LWE_lambda (const gretl_matrix *I1, int n)
{
    gretl_matrix *lambda;
    int i, m = gretl_vector_get_length(I1);

    lambda = gretl_column_vector_alloc(m);

    if (lambda != NULL) {
	for (i=0; i<m; i++) {
	    gretl_vector_set(lambda, i, (M_2PI / n) * (i + 1));
	}
    }

    return lambda;
}

static int
LWE_init (struct LWE_helper *L, const gretl_matrix *X, int m)
{
    int i, err = 0;

    L->I2 = L->lpow = NULL;

    L->I1 = gretl_matrix_pergm(X, m, &err);
    if (err) {
	return err;
    }

    L->lambda = LWE_lambda(L->I1, X->rows);
    if (L->lambda == NULL) {
	gretl_matrix_free(L->I1);
	return E_ALLOC;
    }

    L->lpow = gretl_matrix_copy(L->lambda);
    L->I2 = gretl_matrix_copy(L->I1);

    if (L->lpow == NULL || L->I2 == NULL) {
	err = E_ALLOC;
	LWE_free(L);
    } else {
	L->lcm = 0.0;
	for (i=0; i<m; i++) {
	    L->lcm += log(L->lambda->val[i]);
	}
	L->lcm /= m;
    }

    return err;
}

#define LWE_MAXITER 100

static double LWE_calc (const gretl_matrix *X, int m, int *err)
{
    struct LWE_helper L = {0};
    double d = 0, dd = 1.0;
    double eps = 1.0e-05;
    double f, incr, incl, deriv, h;
    int iter = 0;

    *err = LWE_init(&L, X, m);
    if (*err) {
	return NADBL;
    }

    while (fabs(dd) > 1.0e-06 && iter++ < LWE_MAXITER) {
	f = LWE_obj_func(&L, d);
	incr = LWE_obj_func(&L, d + eps) / eps;
	incl = LWE_obj_func(&L, d - eps) / eps;

	deriv = (incr - incl) / 2.0;
	h = (0.5 * (incr + incl) - f / eps) / eps;
	dd = (h < 0)? (-deriv / h) : deriv;

	if (fabs(dd) > 1) {
	    dd = (dd > 0)? 1 : -1;
	}

	d += 0.5 * dd;
    }

    if (*err) {
	d = NADBL;
    } else if (iter == LWE_MAXITER) {
	fprintf(stderr, "LWE: max iterations reached\n");
	d = NADBL;
    }

    LWE_free(&L);

    return d;
}

int auto_spectrum_order (int T, gretlopt opt)
{
    int m;

    if (opt & OPT_O) {
	/* Bartlett */
	m = (int) 2.0 * sqrt((double) T);
    } else {
	/* fractional integration test */
	double m1 = floor((double) T / 2.0);
	double m2 = floor(pow((double) T, 0.6));

	m = (m1 < m2)? m1 : m2;
	m--;
    }

    return m;
}

static int fract_int_LWE (const double *x, int m, int t1, int t2,
			  struct fractint_test *ft)
{
    gretl_matrix *X;
    int err = 0;

    X = gretl_vector_from_series(x, t1, t2);

    if (X == NULL) {
	err = E_ALLOC;
    } else {
	ft->d = LWE_calc(X, m, &err);
	if (!err) {
	    ft->se = 1.0 / (2 * sqrt((double) m));
	}
	gretl_matrix_free(X);
    }

    return err;
}

/* called by pergm_or_fractint(), when responding to the
   "pergm" command only, without its --silent option
*/

static void pergm_print (const char *vname, const double *d,
			 int T, int L, gretlopt opt, PRN *prn)
{
    char xstr[32];
    double dt, yt;
    int t;

    if (vname == NULL) {
	pprintf(prn, "\n%s\n", _("Residual periodogram"));
    } else {
	pprintf(prn, _("\nPeriodogram for %s\n"), vname);
    }

    pprintf(prn, _("Number of observations = %d\n"), T);
    if (opt & OPT_O) {
	pprintf(prn, _("Using Bartlett lag window, length %d\n\n"), L);
    } else {
	pputc(prn, '\n');
    }

    if (opt & OPT_L) {
	pputs(prn, _(" omega  scaled frequency  periods  log spectral density\n\n"));
    } else {
	pputs(prn, _(" omega  scaled frequency  periods  spectral density\n\n"));
    }

    for (t=1; t<=T/2; t++) {
	dt = (opt & OPT_L)? log(d[t-1]) : d[t-1];
	if (opt & OPT_R) {
	    yt = 2 * M_PI * (double) t / T;
	    pprintf(prn, " %.5f%8d%16.2f", yt, t, (double) T / t);
	} else if (opt & OPT_D) {
	    yt = 360 * t / (double) T;
	    pprintf(prn, " %7.3f%8d%16.2f", yt, t, (double) T / t);
	} else {
	    yt = M_2PI * t / (double) T;
	    pprintf(prn, " %.5f%8d%16.2f", yt, t, (double) T / t);
	}
	dt = (opt & OPT_L)? log(d[t-1]) : d[t-1];
	sprintf(xstr, "%#.5g", dt);
	gretl_fix_exponent(xstr);
	pprintf(prn, "%16s\n", xstr);
    }

    pputc(prn, '\n');
}

static int finalize_fractint (const double *x,
			      const double *dens,
			      int t1, int t2, int width,
			      const char *vname,
			      gretlopt opt,
			      PRN *prn)
{
    struct fractint_test ft;
    gretl_matrix *result = NULL;
    int do_GPH = opt & (OPT_G | OPT_A);
    int do_LWE = !(opt & OPT_G);
    int do_print = !(opt & OPT_Q);
    int T = t2 - t1 + 1;
    int m, err = 0;

    /* order for test */
    if (width <= 0) {
	m = auto_spectrum_order(T, OPT_NONE);
    } else {
	m = (width > T / 2)? T / 2 : width;
    }

    if (do_print && vname != NULL) {
	pprintf(prn, "\n%s, T = %d\n\n", vname, T);
    }

    if (do_GPH && do_LWE) {
	result = gretl_matrix_alloc(2, 2);
    } else {
	result = gretl_matrix_alloc(1, 2);
    }

    if (do_LWE) {
	/* do Local Whittle if wanted */
	err = fract_int_LWE(x, m, t1, t2, &ft);
	if (!err) {
	    double z = ft.d / ft.se;
	    double pv = normal_pvalue_2(z);

	    record_test_result(z, pv);
	    gretl_matrix_set(result, 0, 0, ft.d);
	    gretl_matrix_set(result, 0, 1, ft.se);

	    if (do_print) {
		pprintf(prn, "%s (m = %d)\n"
			"  %s = %g (%g)\n"
			"  %s: z = %g, %s %.4f\n\n",
			_("Local Whittle Estimator"), m,
			_("Estimated degree of integration"), ft.d, ft.se,
			_("test statistic"), z,
			_("with p-value"), pv);
	    }
	}
    }

    if (!err && do_GPH) {
	/* do GPH if wanted (options --all or --gph) */
	int row = do_LWE ? 1 : 0;

	err = fract_int_GPH(T, m, dens, &ft);
	gretl_matrix_set(result, row, 0, ft.d);
	gretl_matrix_set(result, row, 1, ft.se);

	if (!err && (!do_LWE || do_print)) {
	    double tval = ft.d / ft.se;
	    int df = m - 2;
	    double pv = student_pvalue_2(df, tval);

	    if (!do_LWE) {
		record_test_result(tval, pv);
	    }
	    if (do_print) {
		pprintf(prn, "%s (m = %d)\n"
			"  %s = %g (%g)\n"
			"  %s: t(%d) = %g, %s %.4f\n\n",
			_("GPH test"), m,
			_("Estimated degree of integration"), ft.d, ft.se,
			_("test statistic"), df, tval,
			_("with p-value"), pv);
	    }
	}
    }

    if (err) {
	gretl_matrix_free(result);
    } else {
	char **colnames = strings_array_new(2);

	colnames[0] = gretl_strdup("d");
	colnames[1] = gretl_strdup("se");
	gretl_matrix_set_colnames(result, colnames);
	set_last_result_data(result, GRETL_TYPE_MATRIX);
    }

    return err;
}

static double *pergm_bartlett_density (const double *x, int t1, int t2,
				       int L, int *err)
{
    double *sdy, *acov, *dens;
    double xx, yy, vx, sx, w;
    int k, t, T = t2 - t1 + 1;

    sdy = malloc(T * sizeof *sdy);
    acov = malloc((L + 1) * sizeof *acov);
    dens = malloc(T/2 * sizeof *dens);

    if (sdy == NULL || acov == NULL || dens == NULL) {
	*err = E_ALLOC;
	free(sdy);
	free(acov);
	free(dens);
	return NULL;
    }

    xx = gretl_mean(t1, t2, x);
    vx = real_gretl_variance(t1, t2, x, 1);
    sx = sqrt(vx);

    for (t=t1; t<=t2; t++) {
	sdy[t-t1] = (x[t] - xx) / sx;
    }

    /* autocovariances */
    for (k=1; k<=L; k++) {
	acov[k] = 0.0;
	for (t=k; t<T; t++) {
	    acov[k] += sdy[t] * sdy[t-k];
	}
	acov[k] /= T;
    }

    vx /= M_2PI;

    for (t=0; t<T/2; t++) {
	yy = M_2PI * (t+1) / (double) T;
	xx = 1.0;
	for (k=1; k<=L; k++) {
	    w = 1.0 - (double) k/(L + 1);
	    xx += 2.0 * w * acov[k] * cos(yy * k);
	}
	dens[t] = xx * vx;
    }

    free(sdy);
    free(acov);

    return dens;
}

enum {
    PERGM_CMD,
    PERGM_FUNC,
    FRACTINT_CMD
};

/* The following is a little complicated because we're supporting
   three use cases:

   1) Implementing the "pergm" command, with or without the
   option of using a Bartlett window (OPT_O); in this case the
   final matrix-location argument will be NULL.

   2) Implementing the user-space pergm function: in this case
   the final matrix-location argument will be non-NULL, and the
   @vname argument will be NULL.

   3) Implementing the fractint command, computing the Local
   Whittle Estimator and/or the GPH test -- this is flagged by
   opt & OPT_F. If we're doing Whittle only, we don't need to compute
   the spectral density via the autocovariances approach; that is
   handled via FFT in the LWE functions above.
*/

static int
pergm_or_fractint (int usage, const double *x, int t1, int t2,
		   int width, const char *vname, gretlopt opt,
		   PRN *prn, gretl_matrix **pmat)
{
    double *dens = NULL;
    int bartlett = (opt & OPT_O);
    int whittle_only = 0;
    int t, T, L = 0;
    int do_plot = 0;
    int err = 0;

    gretl_error_clear();

    /* common to all uses: check for data problems */

    err = series_adjust_sample(x, &t1, &t2);
    if (!err && (T = t2 - t1 + 1) < 12) {
	err = E_TOOFEW;
    }
    if (!err && gretl_isconst(t1, t2, x)) {
	if (vname != NULL) {
	    gretl_errmsg_sprintf(_("%s is a constant"), vname);
	}
	err = E_DATA;
    }
    if (err) {
	return err;
    }

    if (usage == PERGM_CMD) {
	do_plot = gnuplot_graph_wanted(PLOT_PERIODOGRAM, opt, &err);
	if (err) {
	    return err;
	}
    }

    if (usage == FRACTINT_CMD) {
	/* doing fractint */
	if (!(opt & (OPT_A | OPT_G))) {
	    /* not --all, not --gph */
	    whittle_only = 1;
	}
    }

    if (!whittle_only) {
	/* Chatfield (1996); William Greene, 4e, p. 772 */
	if (bartlett) {
	    if (width <= 0) {
		L = auto_spectrum_order(T, opt);
	    } else {
		L = (width > T / 2)? T / 2 : width;
	    }
	} else {
	    /* use full sample */
	    L = T - 1;
	}

	if (bartlett) {
	    dens = pergm_bartlett_density(x, t1, t2, L, &err);
	} else {
	    gretl_matrix *pg;
	    gretl_vector vt;

	    gretl_matrix_init(&vt);
	    vt.rows = t2 - t1 + 1;
	    vt.cols = 1;
	    vt.val = (double *) x + t1;

	    pg = gretl_matrix_pergm(&vt, T/2, &err);
	    if (!err) {
		dens = gretl_matrix_steal_data(pg);
	    }
	    gretl_matrix_free(pg);
	}
	if (err) {
	    return err;
	}
    }

    if (usage == PERGM_FUNC) {
	/* make matrix for pergm() function */
	int T2 = T / 2;
	gretl_matrix *pm;

	*pmat = pm = gretl_matrix_alloc(T2, 2);

	if (pm == NULL) {
	    err = E_ALLOC;
	} else {
	    for (t=0; t<T2; t++) {
		gretl_matrix_set(pm, t, 0, M_2PI * (t+1) / (double) T);
		gretl_matrix_set(pm, t, 1, dens[t]);
	    }
	}
    } else if (usage == FRACTINT_CMD) {
	/* supporting the "fractint" command */
	err = finalize_fractint(x, dens, t1, t2, width,
				vname, opt, prn);
    } else {
	/* supporting the "pergm" command */
	if (!(opt & OPT_S)) {
	    /* not --silent */
	    pergm_print(vname, dens, T, L, opt, prn);
	}
	if (do_plot) {
	    err = periodogram_plot(vname, T, L, dens, opt);
	}
    }

    free(dens);

    return err;
}

/**
 * periodogram:
 * @varno: ID number of variable to process.
 * @width: width of window.
 * @dset: dataset struct.
 * @opt: if includes OPT_O, use Bartlett lag window for periodogram;
 * OPT_L, use log scale; OPT_S, don't print anything (just do plot).
 * @prn: gretl printing struct.
 *
 * Computes and displays the periodogram for the series specified
 * by @varno.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int periodogram (int varno, int width, const DATASET *dset,
		 gretlopt opt, PRN *prn)
{
    return pergm_or_fractint(PERGM_CMD, dset->Z[varno],
			     dset->t1, dset->t2,
			     width, dset->varname[varno],
			     opt, prn, NULL);
}

/**
 * residual_periodogram:
 * @x: series to process.
 * @width: width of window.
 * @dset: dataset struct.
 * @opt: if includes OPT_O, use Bartlett lag window for periodogram;
 * OPT_L, use log scale.
 * @prn: gretl printing struct.
 *
 * Computes and displays the periodogram for @x, which is presumed to
 * be a model residual series.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int residual_periodogram (const double *x, int width, const DATASET *dset,
			  gretlopt opt, PRN *prn)
{
    return pergm_or_fractint(PERGM_CMD, x,
			     dset->t1, dset->t2,
			     width, NULL,
			     opt, prn, NULL);
}

/**
 * periodogram_matrix:
 * @x: the series to process.
 * @t1: starting observation in @x.
 * @t2: ending observation in @x.
 * @width: width of Bartlett window, or -1 for plain sample
 * periodogram.
 * @err: location to receive error code.
 *
 * Implements the userspace gretl pergm function, which can
 * be used on either a series from the dataset or a gretl
 * vector.
 *
 * Returns: allocated matrix on success, NULL on failure.
 */

gretl_matrix *periodogram_matrix (const double *x, int t1, int t2,
				  int width, int *err)
{
    gretlopt opt = (width < 0)? OPT_NONE : OPT_O;
    gretl_matrix *m = NULL;

    *err = pergm_or_fractint(PERGM_FUNC, x, t1, t2, width,
			     NULL, opt, NULL, &m);

    return m;
}

/**
 * fractint:
 * @varno: ID number of variable to process.
 * @order: lag order / window size.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Computes and prints a test for fractional integration of the
 * series specified by @varno. By default the test uses the
 * Local Whittle Estimator but if @opt includes OPT_G then
 * the Geweke and Porter-Hudak test is done instead, or if
 * OPT_A then both tests are shown. If OPT_Q is given the
 * test results are not printed, just recorded (with
 * preference given to the LWE in case of OPT_A).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int fractint (int varno, int order, const DATASET *dset,
	      gretlopt opt, PRN *prn)
{
    int err = incompatible_options(opt, (OPT_G | OPT_A));

    if (!err) {
	err = pergm_or_fractint(FRACTINT_CMD, dset->Z[varno],
				dset->t1, dset->t2,
				order,
				dset->varname[varno],
				opt, prn, NULL);
    }

    return err;
}

static void printfw (int w, double x, int d, int ival, PRN *prn)
{
    w = w - 1; /* allow for leading space */
    pputc(prn, ' ');
    if (na(x)) {
	pprintf(prn, "%*s", UTF_WIDTH(_("NA"), w), _("NA"));
    } else if (ival) {
	pprintf(prn, "%*.0f", w, x);
    } else if (x > 999 && x < 100000) {
	int p = 1 + floor(log10(x));

	d -= p;
	if (d < 0) d = 0;
	pprintf(prn, "%*.*f", w, d, x);
    } else {
	pprintf(prn, "%#*.*g", w, d, x);
    }
}

static void output_line (char *str, PRN *prn, int dblspc)
{
    pputs(prn, str);
    if (dblspc) {
	pputs(prn, "\n\n");
    } else {
	pputc(prn, '\n');
    }
}

static void prhdr (const char *str, const DATASET *dset,
		   int missing, PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN];
    gchar *tmp;

    ntolabel(date1, dset->t1, dset);
    ntolabel(date2, dset->t2, dset);

    pputc(prn, '\n');

    tmp = g_strdup_printf(_("%s, using the observations %s - %s"),
			  str, date1, date2);
    output_line(tmp, prn, 0);
    g_free(tmp);

    if (missing) {
	output_line(_("(missing values were skipped)"), prn, 1);
    }
}

static void summary_print_val (int w, double x, int digits,
			       int places, PRN *prn)
{
    pputc(prn, ' ');
    if (na(x)) {
	pprintf(prn, "%*s", UTF_WIDTH(_("NA"), w), _("NA"));
    } else if (digits < 0) {
	pprintf(prn, "%*.0f", w, x);
    } else if (digits > 0) {
	pprintf(prn, "%#*.*g", w, digits, x);
    } else if (places > 0) {
	pprintf(prn, "%*.*f", w, places, x);
    } else {
	/* the default */
	pprintf(prn, "%#*.5g", w, x);
    }
}

#define NSUMM 12

void print_summary_single (const Summary *s,
			   int digits, int places,
			   const DATASET *dset,
			   PRN *prn)
{
    const char *labels[NSUMM] = {
	N_("Mean"),
	N_("Median"),
	N_("Minimum"),
	N_("Maximum"),
	N_("Standard deviation"),
	N_("C.V."),
	N_("Skewness"),
	N_("Ex. kurtosis"),
	/* xgettext:no-c-format */
	N_("5% percentile"),
	/* xgettext:no-c-format */
	N_("95% percentile"),
	N_("Interquartile range"),
	N_("Missing obs.")
    };
    const char *wstr = N_("Within s.d.");
    const char *bstr = N_("Between s.d.");
    double vals[NSUMM];
    int simple_skip[NSUMM] = {0,1,0,0,0,1,1,1,1,1,1,0};
    int skip0595 = 0;
    int offset = 2;
    int slen = 0, i = 0;

    if (s->opt & OPT_B) {
	offset = 4;
    } else {
	const char *vname = dset->varname[s->list[1]];
	char obs1[OBSLEN], obs2[OBSLEN];
	gchar *tmp = NULL;

	ntolabel(obs1, dset->t1, dset);
	ntolabel(obs2, dset->t2, dset);

	prhdr(_("Summary statistics"), dset, 0, prn);

	if (isdigit(*vname)) {
	    const char *mname = dataset_get_matrix_name(dset);

	    if (mname != NULL) {
		tmp = g_strdup_printf(_("for column %d of %s (%d valid observations)"),
				      atoi(vname), mname, s->n);
	    } else {
		tmp = g_strdup_printf(_("for column %d (%d valid observations)"),
				      atoi(vname), s->n);
	    }
	} else {
	    tmp = g_strdup_printf(_("for the variable '%s' (%d valid observations)"),
				  dset->varname[s->list[1]], s->n);
	}
	output_line(tmp, prn, 1);
	g_free(tmp);
    }

    vals[0]  = s->mean[0];
    vals[1]  = s->median[0];
    vals[2]  = s->low[0];
    vals[3]  = s->high[0];
    vals[4]  = s->sd[0];
    vals[5]  = s->cv[0];
    vals[6]  = s->skew[0];
    vals[7]  = s->xkurt[0];
    vals[8]  = s->perc05[0];
    vals[9]  = s->perc95[0];
    vals[10] = s->iqr[0];
    vals[11] = s->misscount[0];

    if (na(vals[8]) && na(vals[9])) {
	skip0595 = 1;
    }

    for (i=0; i<NSUMM; i++) {
	if ((s->opt & OPT_S) && simple_skip[i]) {
	    continue;
	} else if (skip0595 && (i == 8 || i == 9)) {
	    continue;
	}
	if (strlen(_(labels[i])) > slen) {
	    slen = g_utf8_strlen(_(labels[i]), -1);
	}
    }
    slen++;

    for (i=0; i<NSUMM; i++) {
	if ((s->opt & OPT_S) && simple_skip[i]) {
	    continue;
	} else if (skip0595 && (i == 8 || i == 9)) {
	    continue;
	}
	bufspace(offset, prn);
	if (i == 8 || i == 9) {
	    /* the strings contain '%' */
	    gchar *pcstr = g_strdup(_(labels[i]));
	    int n = slen - g_utf8_strlen(pcstr, -1);

	    pputs(prn, pcstr);
	    if (n > 0) {
		bufspace(n, prn);
	    }
	    g_free(pcstr);
	} else {
	    pprintf(prn, "%-*s", UTF_WIDTH(_(labels[i]), slen), _(labels[i]));
	}
	if (i == NSUMM - 1 || s->ival[0]) {
	    summary_print_val(14, vals[i], -1, places, prn);
	} else {
	    summary_print_val(14, vals[i], digits, places, prn);
	}
	pputc(prn, '\n');
    }

    if (!na(s->sw) && !na(s->sb)) {
	pputc(prn, '\n');
	bufspace(offset, prn);
	pprintf(prn, "%-*s", UTF_WIDTH(_(wstr), slen), _(wstr));
	summary_print_val(14, s->sw, digits, places, prn);
	pputc(prn, '\n');
	bufspace(offset, prn);
	pprintf(prn, "%-*s", UTF_WIDTH(_(bstr), slen), _(bstr));
	summary_print_val(14, s->sb, digits, places, prn);
    }

    pputc(prn, '\n');
}

static void summary_print_varname (const char *src,
				   int len, PRN *prn)
{
    char vname[NAMETRUNC];

    maybe_trim_varname(vname, src);
    pprintf(prn, "%-*s", len, vname);
}

/**
 * print_summary:
 * @summ: pointer to gretl summary statistics struct.
 * @dset: information on the data set.
 * @prn: gretl printing struct.
 *
 * Prints the summary statistics for a given variable.
 */

void print_summary (const Summary *summ,
		    const DATASET *dset,
		    PRN *prn)
{
    int dmax, d = get_gretl_digits();
    int nv = summ->list[0];
    int len = 0;
    int i, vi;

    if (summ->list == NULL || summ->list[0] == 0) {
	return;
    }

    if (summ->weight_var > 0) {
	pputc(prn, '\n');
	pprintf(prn, _("Weighting variable: %s\n"),
		       dset->varname[summ->weight_var]);
    }

    if (nv == 1 && !(summ->opt & OPT_B)) {
	print_summary_single(summ, 6, 0, dset, prn);
	return;
    }

    /* number of significant figures to use */
    dmax = 5; // was: (summ->opt & OPT_S)? 4 : 5;
    d = d > dmax ? dmax : d;

    if (nv > 1) {
	int maxlen = max_namelen_in_list(summ->list, dset);

	len = maxlen + 1;
    }

    pputc(prn, '\n');

    if (summ->opt & OPT_S) {
	/* the "simple" option */
	const char *h[] = {
	    N_("Mean"),
	    N_("Median"),
	    /* TRANSLATORS: 'S.D.' means Standard Deviation */
	    N_("S.D."),
	    N_("Min"),
	    N_("Max"),
	};
	int w = 12;
	int ival;

	if (len > 0) {
	    pprintf(prn, "%*s", len, " ");
	}
	pprintf(prn, "%*s%*s%*s%*s%*s\n",
		UTF_WIDTH(_(h[0]), w), _(h[0]),
		UTF_WIDTH(_(h[1]), w), _(h[1]),
		UTF_WIDTH(_(h[2]), w), _(h[2]),
		UTF_WIDTH(_(h[3]), w), _(h[3]),
		UTF_WIDTH(_(h[4]), w), _(h[4]));

	for (i=0; i<summ->list[0]; i++) {
	    vi = summ->list[i+1];
	    ival = summ->ival[i];
	    if (nv > 1) {
		summary_print_varname(dset->varname[vi], len, prn);
	    }
	    printfw(w, summ->mean[i], d, ival, prn);
	    printfw(w, summ->median[i], d, ival, prn);
	    printfw(w, summ->sd[i], d, ival, prn);
	    printfw(w, summ->low[i], d, ival, prn);
	    printfw(w, summ->high[i], d, ival, prn);
	    pputc(prn, '\n');
	}
	pputc(prn, '\n');
    } else {
	/* print all available stats */
	const char *ha[] = {
	    N_("Mean"),
	    N_("Median"),
	    N_("Minimum"),
	    N_("Maximum"),
	};
	const char *hb[] = {
	    N_("Std. Dev."),
	    N_("C.V."),
	    N_("Skewness"),
	    N_("Ex. kurtosis")
	};
	const char *hc[] = {
	    /* xgettext:no-c-format */
	    N_("5% perc."),
	    /* xgettext:no-c-format */
	    N_("95% perc."),
	    N_("IQ range"),
	    N_("Missing obs.")
	};
	/* cases where 0.05 and 0.95 quantiles are OK */
	int npct = 0;
	int w = 15;
	int ival;

	if (len > 0) {
	    pprintf(prn, "%*s", len, " ");
	}
	pprintf(prn, "%*s%*s%*s%*s\n",
		UTF_WIDTH(_(ha[0]), w), _(ha[0]),
		UTF_WIDTH(_(ha[1]), w), _(ha[1]),
		UTF_WIDTH(_(ha[2]), w), _(ha[2]),
		UTF_WIDTH(_(ha[3]), w), _(ha[3]));

	for (i=0; i<summ->list[0]; i++) {
	    vi = summ->list[i+1];
	    ival = summ->ival[i];
	    if (nv > 1) {
		summary_print_varname(dset->varname[vi], len, prn);
	    }
	    printfw(w, summ->mean[i], d, ival, prn);
	    printfw(w, summ->median[i], d, ival, prn);
	    printfw(w, summ->low[i], d, ival, prn);
	    printfw(w, summ->high[i], d, ival, prn);
	    pputc(prn, '\n');
	    /* while we're at it, register cases where we can
	       show the 0.05 and 0.95 quantiles
	    */
	    if (!na(summ->perc05[i]) && !na(summ->perc95[i])) {
		npct++;
	    }
	}
	pputc(prn, '\n');

	if (len > 0) {
	    pprintf(prn, "%*s", len, " ");
	}
	pprintf(prn, "%*s%*s%*s%*s\n",
		UTF_WIDTH(_(hb[0]), w), _(hb[0]),
		UTF_WIDTH(_(hb[1]), w), _(hb[1]),
		UTF_WIDTH(_(hb[2]), w), _(hb[2]),
		UTF_WIDTH(_(hb[3]), w), _(hb[3]));

	for (i=0; i<summ->list[0]; i++) {
	    double cv;

	    vi = summ->list[i+1];
	    ival = summ->ival[i];
	    if (nv > 1) {
		summary_print_varname(dset->varname[vi], len, prn);
	    }
	    if (floateq(summ->mean[i], 0.0)) {
		cv = NADBL;
	    } else if (floateq(summ->sd[i], 0.0)) {
		cv = 0.0;
	    } else {
		cv = fabs(summ->sd[i] / summ->mean[i]);
	    }
	    printfw(w, summ->sd[i], d, ival, prn);
	    printfw(w, cv, d, 0, prn);
	    printfw(w, summ->skew[i], d, 0, prn);
	    printfw(w, summ->xkurt[i], d, 0, prn);
	    pputc(prn, '\n');
	}
	pputc(prn, '\n');

	if (npct > 0) {
	    /* note: use pputs for strings containing literal '%' */
	    gchar *hc0 = g_strdup(_(hc[0]));
	    gchar *hc1 = g_strdup(_(hc[1]));
	    int n;

	    if (len > 0) {
		pprintf(prn, "%*s", len, " ");
	    }
	    n = w - g_utf8_strlen(hc0, -1);
	    if (n > 0) bufspace(n, prn);
	    pputs(prn, hc0);
	    n = w - g_utf8_strlen(hc1, -1);
	    if (n > 0) bufspace(n, prn);
	    pputs(prn, hc1);
	    g_free(hc0); g_free(hc1);

	    pprintf(prn, "%*s%*s\n",
		    UTF_WIDTH(_(ha[2]), w), _(hc[2]),
		    UTF_WIDTH(_(ha[3]), w), _(hc[3]));
	} else {
	    /* not showing any 0.05, 0.95 quantiles */
	    if (len > 0) {
		pprintf(prn, "%*s", len, " ");
	    }
	    pprintf(prn, "%*s%*s\n",
		    UTF_WIDTH(_(ha[2]), w), _(hc[2]),
		    UTF_WIDTH(_(ha[3]), w), _(hc[3]));
	}

	for (i=0; i<summ->list[0]; i++) {
	    vi = summ->list[i+1];
	    ival = summ->ival[i];
	    if (nv > 1) {
		summary_print_varname(dset->varname[vi], len, prn);
	    }
	    if (!na(summ->perc05[i]) && !na(summ->perc95[i])) {
		printfw(w, summ->perc05[i], d, ival, prn);
		printfw(w, summ->perc95[i], d, ival, prn);
	    } else if (npct > 0) {
		pprintf(prn, "%*s", w, "NA");
		pprintf(prn, "%*s", w, "NA");
	    }
	    printfw(w, summ->iqr[i], d, ival, prn);
	    pprintf(prn, "%*d", w, (int) summ->misscount[i]);
	    pputc(prn, '\n');
	}
	pputc(prn, '\n');
    }
}

/**
 * free_summary:
 * @summ: pointer to gretl summary statistics struct
 *
 * Frees all resources associated with @summ, and
 * the pointer itself.
 */

void free_summary (Summary *summ)
{
    free(summ->list);
    free(summ->misscount);
    free(summ->stats);
    free(summ->ival);

    free(summ);
}

/***
 * summary_new:
 * @list: list of variables we want summary statistics for
 * @wv: weighting variable (0=no weights)
 * @opt: options
 *
 * Allocates a new "Summary" struct and initializes a few
 * things inside it
 */

static Summary *summary_new (const int *list, int wv,
			     gretlopt opt, int *err)
{
    Summary *s;
    int nv;

    if (list == NULL) {
	fprintf(stderr, "summary_new: list is NULL\n");
	*err = E_DATA;
	return NULL;
    }

    nv = list[0];

    s = malloc(sizeof *s);
    if (s == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    s->list = gretl_list_copy(list);
    if (s->list == NULL) {
	*err = E_ALLOC;
	free(s);
	return NULL;
    }

    s->weight_var = wv;
    s->opt = opt;
    s->n = 0;
    s->misscount = malloc(nv * sizeof *s->misscount);
    s->ival = calloc(nv, sizeof *s->ival);

    s->stats = malloc(11 * nv * sizeof *s->stats);
    if (s->stats == NULL) {
	*err = E_ALLOC;
	free_summary(s);
	return NULL;
    }

    s->mean   = s->stats;
    s->median = s->mean + nv;
    s->sd     = s->median + nv;
    s->skew   = s->sd + nv;
    s->xkurt  = s->skew + nv;
    s->low    = s->xkurt + nv;
    s->high   = s->low + nv;
    s->cv     = s->high + nv;
    s->perc05 = s->cv + nv;
    s->perc95 = s->perc05 + nv;
    s->iqr    = s->perc95 + nv;

    s->sb = s->sw = NADBL;

    return s;
}

int summary_has_missing_values (const Summary *summ)
{
    if (summ->misscount != NULL) {
	int i, nv = summ->list[0];

	for (i=0; i<nv; i++) {
	    if (summ->misscount[i] > 0) {
		return 1;
	    }
	}
    }

    return 0;
}

static int compare_wgt_ord_rows (const void *a, const void *b)
{
    const double **da = (const double **) a;
    const double **db = (const double **) b;

    return (da[0][0] > db[0][0]) - (da[0][0] < db[0][0]);
}

static int weighted_order_stats (const double *y, const double *w,
				 double *ostats, int n, int k)
{
    /* on input "ostats" contains the quantiles to compute;
       on output it holds the results
    */
    double p, q, wsum = 0.0;
    double **X;
    int t, i, err = 0;

    X = doubles_array_new(n, 2);
    if (X == NULL) {
	return E_ALLOC;
    }

    i = 0;
    for (t=0; t<n; t++) {
	if (na(y[t]) || na(w[t]) || w[t] == 0.0) {
	    continue;
	} else {
	    X[i][0] = y[t];
	    X[i][1] = w[t];
	    wsum += w[t];
	    i++;
	}
    }

    qsort(X, i, sizeof *X, compare_wgt_ord_rows);

    for (i=0; i<k; i++) {
	p = ostats[i] * wsum;
	q = X[0][1];
	t = 0;
	while (q <= p) {
	    q += X[t][1];
	    if (q < p) t++;
	}

	if (t == 0 || t >= n - 1) {
	    ostats[i] = NADBL;
	    continue;
	}

	if (X[t-1][0] == X[t][0]) {
	    ostats[i] = X[t][0];
	} else {
	    double a = (q - p) / X[t][1];

	    ostats[i] = a * X[t+1][0] + (1-a) * X[t][0];
	}
    }

    doubles_array_free(X, n);

    return err;
}

/**
 * get_summary_weighted:
 * @list: list of variables to process.
 * @dset: dataset struct.
 * @wgt: series to use as weights.
 * @opt: may include OPT_S for "simple" version.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Calculates descriptive summary statistics for the specified
 * variables, weighting the observations by @rv. The series @rv must
 * be of full length (dset->n).
 *
 * Returns: #Summary object containing the summary statistics, or NULL
 * on failure.
 */

Summary *get_summary_weighted (const int *list, const DATASET *dset,
			       int wtvar, gretlopt opt,
			       PRN *prn, int *err)
{
    const double *wts;
    double *x, ostats[5];
    int t1 = dset->t1;
    int t2 = dset->t2;
    Summary *s;
    int i, t;

    s = summary_new(list, wtvar, opt, err);
    if (s == NULL) {
	return NULL;
    }

    x = malloc(dset->n * sizeof *x);
    if (x == NULL) {
	*err = E_ALLOC;
	free_summary(s);
	return NULL;
    }

    wts = dset->Z[wtvar];

    for (i=0; i<list[0]; i++)  {
	double *pskew = NULL, *pkurt = NULL;
	int vi = s->list[i+1];
	int ni = 0;
        int ntot = 0;
        double wt;

	/* prepare the series for weighting: substitute NAs for values
	   at which the weights are invalid or zero
	*/
	for (t=t1; t<=t2; t++) {
	    wt = wts[t];
	    if (!na(wt) && wt != 0.0) {
                ntot++;
		x[t] = dset->Z[vi][t];
		if (!na(x[t])) {
		    ni++;
		}
	    } else {
		x[t] = NADBL;
	    }
	}

	s->misscount[i] = ntot - ni;

	if (ni > s->n) {
	    s->n = ni;
	}

	if (ni == 0) {
	    pprintf(prn, _("Dropping %s: sample range contains no valid "
			   "observations\n"), dset->varname[vi]);
	    gretl_list_delete_at_pos(s->list, i + 1);
	    if (s->list[0] == 0) {
		return s;
	    } else {
		i--;
		continue;
	    }
	}

	if (opt & OPT_S) {
	    s->skew[i] = NADBL;
	    s->xkurt[i] = NADBL;
	    s->cv[i] = NADBL;
	    s->median[i] = NADBL;
	} else {
	    pskew = &s->skew[i];
	    pkurt = &s->xkurt[i];
	}

	summary_minmax(t1, t2, x, s, i);
	gretl_moments(t1, t2, x, wts, &s->mean[i], &s->sd[i], pskew, pkurt, 1);

	if (!(opt & OPT_S)) {
	    int err;

	    if (floateq(s->mean[i], 0.0)) {
		s->cv[i] = NADBL;
	    } else if (floateq(s->sd[i], 0.0)) {
		s->cv[i] = 0.0;
	    } else {
		s->cv[i] = fabs(s->sd[i] / s->mean[i]);
	    }

	    ostats[0] = 0.05;
	    ostats[1] = 0.25;
	    ostats[2] = 0.5;
	    ostats[3] = 0.75;
	    ostats[4] = 0.95;

	    err = weighted_order_stats(x, wts, ostats, ni, 5);

	    if (!err) {
		s->median[i] = ostats[2];
		s->perc05[i] = ostats[0];
		s->perc95[i] = ostats[4];
		if (!na(ostats[1]) && !na(ostats[3])) {
		    s->iqr[i] = ostats[3] - ostats[1];
		}
	    }
	}
#if 0
	if (dataset_is_panel(dset) && list[0] == 1) {
	    panel_variance_info(x, dset, s->mean[0], &s->sw, &s->sb);
	}
#endif
    }

    free(x);

    return s;
}

/* Get the additional statistics wanted if the --simple
   flag was not given to the summary command.
*/

static int get_extra_stats (Summary *s, int i,
			    int t1, int t2,
			    const double *x)
{
    int err = 0;

    if (floateq(s->mean[i], 0.0)) {
	s->cv[i] = NADBL;
    } else if (floateq(s->sd[i], 0.0)) {
	s->cv[i] = 0.0;
    } else {
	s->cv[i] = fabs(s->sd[i] / s->mean[i]);
    }

    s->perc05[i] = gretl_quantile(t1, t2, x, 0.05, OPT_Q, &err);
    s->perc95[i] = gretl_quantile(t1, t2, x, 0.95, OPT_Q, &err);
    s->iqr[i] = gretl_quantile(t1, t2, x, 0.75, OPT_NONE, &err);
    if (!na(s->iqr[i])) {
	s->iqr[i] -= gretl_quantile(t1, t2, x, 0.25, OPT_NONE, &err);
    }

    return err;
}

/**
 * get_summary_restricted:
 * @list: list of variables to process.
 * @dset: dataset struct.
 * @rv: series to use as restriction dummy.
 * @opt: may include OPT_S for "simple" version.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Calculates descriptive summary statistics for the specified variables,
 * with the observations restricted to those for which @rv has a non-zero
 * (and non-missing) value. The series @rv must be of full length (dset->n).
 *
 * Returns: #Summary object containing the summary statistics, or NULL
 * on failure.
 */

Summary *get_summary_restricted (const int *list, const DATASET *dset,
				 const double *rv, gretlopt opt,
				 PRN *prn, int *err)
{
    int t1 = dset->t1;
    int t2 = dset->t2;
    Summary *s;
    double *x;
    int i, t;

    s = summary_new(list, 0, opt, err);
    if (s == NULL) {
	return NULL;
    }

    x = malloc(dset->n * sizeof *x);
    if (x == NULL) {
	*err = E_ALLOC;
	free_summary(s);
	return NULL;
    }

    for (i=0; i<s->list[0]; i++)  {
	double *pskew = NULL, *pkurt = NULL;
	int vi = s->list[i+1];
	int strvals;
	int ni = 0;
        int ntot = 0;

	strvals = is_string_valued(dset, vi);

	if (!strvals) {
	    /* create the restricted series: substitute NAs
	       for values at which the restriction dummy is
	       invalid or zero
	    */
	    for (t=t1; t<=t2; t++) {
		if (!na(rv[t]) && rv[t] != 0.0) {
		    ntot++;
		    x[t] = dset->Z[vi][t];
		    if (!na(x[t])) {
			ni++;
		    }
		} else {
		    x[t] = NADBL;
		}
	    }
	    s->misscount[i] = ntot - ni;
	    if (ni > s->n) {
		s->n = ni;
	    }
	}

	if (ni == 0) {
	    if (strvals) {
		pprintf(prn, _("Dropping %s: string-valued series\n"),
			dset->varname[vi]);
	    } else {
		pprintf(prn, _("Dropping %s: sample range contains no valid "
			       "observations\n"), dset->varname[vi]);
	    }
	    gretl_list_delete_at_pos(s->list, i + 1);
	    if (s->list[0] == 0) {
		return s;
	    } else {
		i--;
		continue;
	    }
	}

	if (opt & OPT_S) {
	    s->skew[i] = NADBL;
	    s->xkurt[i] = NADBL;
	    s->cv[i] = NADBL;
	} else {
	    pskew = &s->skew[i];
	    pkurt = &s->xkurt[i];
	}

	summary_minmax(t1, t2, x, s, i);
	gretl_moments(t1, t2, x, NULL, &s->mean[i], &s->sd[i], pskew, pkurt, 1);
	s->median[i] = gretl_median(t1, t2, x);

	if (!(opt & OPT_S)) {
	    *err = get_extra_stats(s, i, t1, t2, x);
	}

	if (dataset_is_panel(dset) && list[0] == 1) {
	    panel_variance_info(x, dset, s->mean[0], &s->sw, &s->sb);
	}
    }

    free(x);

    return s;
}

/**
 * get_summary:
 * @list: list of variables to process.
 * @dset: dataset struct.
 * @opt: may include OPT_S for "simple" version.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Calculates descriptive summary statistics for the specified variables.
 *
 * Returns: #Summary object containing the summary statistics, or NULL
 * on failure.
 */

Summary *get_summary (const int *list, const DATASET *dset,
		      gretlopt opt, PRN *prn, int *err)
{
    int t1 = dset->t1;
    int t2 = dset->t2;
    Summary *s;
    int i, nmax;

    s = summary_new(list, 0, opt, err);
    if (s == NULL) {
	return NULL;
    }

    nmax = sample_size(dset);

    for (i=0; i<s->list[0]; i++)  {
	double *pskew = NULL, *pkurt = NULL;
	const double *x;
	double x0;
	int vi = s->list[i+1];
	int strvals;
	int ni = 0;

	strvals = is_string_valued(dset, vi);
	if (!strvals) {
	    x = dset->Z[vi];
	    ni = good_obs(x + t1, nmax, &x0);
	    s->misscount[i] = nmax - ni;
	    if (ni > s->n) {
		s->n = ni;
	    }
	}
	if (ni == 0) {
	    if (strvals) {
		pprintf(prn, _("Dropping %s: string-valued series\n"),
			dset->varname[vi]);
	    } else {
		pprintf(prn, _("Dropping %s: sample range contains no valid "
			       "observations\n"), dset->varname[vi]);
	    }
	    gretl_list_delete_at_pos(s->list, i + 1);
	    if (s->list[0] == 0) {
		return s;
	    } else {
		i--;
		continue;
	    }
	}

	if (opt & OPT_S) {
	    s->skew[i] = NADBL;
	    s->xkurt[i] = NADBL;
	    s->cv[i] = NADBL;
	    s->median[i] = NADBL;
	} else {
	    pskew = &s->skew[i];
	    pkurt = &s->xkurt[i];
	}

	summary_minmax(t1, t2, x, s, i);

	if (s->weight_var == 0) {
	    gretl_moments(t1, t2, x, NULL, &s->mean[i], &s->sd[i],
			  pskew, pkurt, 1);
	} else {
	    gretl_moments(t1, t2, x, dset->Z[s->weight_var], &s->mean[i], &s->sd[i],
			  pskew, pkurt, 0);
	}

	/* included in both simple and full variants */
	s->median[i] = gretl_median(t1, t2, x);

	if (!(opt & OPT_S)) {
	    *err = get_extra_stats(s, i, t1, t2, x);
	}

	if (dataset_is_panel(dset) && list[0] == 1) {
	    panel_variance_info(x, dset, s->mean[0], &s->sw, &s->sb);
	}
    }

    return s;
}

/* basic summary stats column heads */
static const char *hb[] = {
    N_("Byval"), N_("Mean"), N_("Median"), N_("S.D."), N_("Min"), N_("Max")
};

/* full summary stats column heads */
static const char *hf[] = {
    N_("Byval"), N_("Mean"), N_("Median"), N_("Minimum"), N_("Maximum"),
    N_("Std. Dev."), N_("C.V."), N_("Skewness"), N_("Ex. kurtosis"),
    /* xgettext:no-c-format */
    N_("5% perc."),
    /* xgettext:no-c-format */
    N_("95% perc."),
    N_("IQ range"),
    N_("Missing obs.")
};

static int add_summary_colnames (gretl_matrix *m,
				 gretlopt opt)
{
    char **S = strings_array_new(m->cols);

    if (S == NULL) {
	return E_ALLOC;
    } else {
	const char **h = (opt & OPT_S)? hb : hf;
	int k = (opt & OPT_B) ? 0 : 1;
	int i;

	for (i=0; i<m->cols; i++) {
	    S[i] = gretl_strdup(_(h[k++]));
	}
	gretl_matrix_set_colnames(m, S);
    }

    return 0;
}

/* the factorized statistics case */

static int add_summary_by_rownames (gretl_matrix *m,
				    const gretl_matrix *xvals,
				    const DATASET *dset,
				    series_table *st)
{
    char **S = strings_array_new(m->rows);
    gchar *tmp;
    double xvi;
    int i;

    if (S == NULL) {
	return E_ALLOC;
    } else {
	for (i=0; i<m->rows; i++) {
	    xvi = xvals->val[i];
	    if (st != NULL) {
		S[i] = gretl_strdup(series_table_get_string(st, xvi));
	    } else {
		tmp = gretl_strdup_printf("%g", xvi);
		S[i] = gretl_strdup(tmp);
		g_free(tmp);
	    }
	}
	gretl_matrix_set_rownames(m, S);
    }

    return 0;
}

/* the (relatively) simple case */

static int add_summary_rownames (gretl_matrix *m,
				 const int *list,
				 const DATASET *dset)
{
    char **S = strings_array_new(m->rows);
    int i, vi;

    if (S == NULL) {
	return E_ALLOC;
    } else {
	for (i=0; i<m->rows; i++) {
	    vi = list[i+1];
	    S[i] = gretl_strdup(dset->varname[vi]);
	}
	gretl_matrix_set_rownames(m, S);
    }

    return 0;
}

static gretl_matrix *make_summary_recorder (const int *list,
					    gretlopt opt,
					    const DATASET *dset,
					    const gretl_matrix *xvals,
					    series_table *st,
					    int *err)
{
    gretl_matrix *m = NULL;
    int rows, cols;

    cols = (opt & OPT_S)? 5 : 12;

    if (opt & OPT_B) {
	/* add an extra column for factor value */
	cols++;
	/* number of distinct factor values */
	rows = gretl_vector_get_length(xvals);
    } else {
	rows = list[0];
    }

    m = gretl_zero_matrix_new(rows, cols);

    if (m == NULL) {
	*err = E_ALLOC;
    } else {
	*err = add_summary_colnames(m, opt);
	if (!*err) {
	    if (xvals != NULL) {
		*err = add_summary_by_rownames(m, xvals, dset, st);
	    } else {
		*err = add_summary_rownames(m, list, dset);
	    }
	}
    }

    if (*err && m != NULL) {
	gretl_matrix_free(m);
	m = NULL;
    }

    return m;
}

/* If @ival is NULL, write the content of @summ into matrix @m,
   otherwise
*/

static void transcribe_summary (Summary *summ, gretl_matrix *m,
				int *ival, double byval)
{
    int factorized = ival != NULL;
    int nv = summ->list[0];
    int offset = 0;
    int i, r;

    if (factorized) {
	/* writing to a single row */
	gretl_matrix_set(m, *ival, 0, byval);
	offset = 1;
    }

    if (summ->opt & OPT_S) {
	/* the "simple" case */
	for (i=0; i<nv; i++) {
	    r = factorized ? *ival : i;
	    gretl_matrix_set(m, r, 0 + offset, summ->mean[i]);
	    gretl_matrix_set(m, r, 1 + offset, summ->median[i]);
	    gretl_matrix_set(m, r, 2 + offset, summ->sd[i]);
	    gretl_matrix_set(m, r, 3 + offset, summ->low[i]);
	    gretl_matrix_set(m, r, 4 + offset, summ->high[i]);
	}
    } else {
	/* record all available stats */
	double cv;

	for (i=0; i<nv; i++) {
	    r = factorized ? *ival : i;
	    gretl_matrix_set(m, r, 0 + offset, summ->mean[i]);
	    gretl_matrix_set(m, r, 1 + offset, summ->median[i]);
	    gretl_matrix_set(m, r, 2 + offset, summ->low[i]);
	    gretl_matrix_set(m, r, 3 + offset, summ->high[i]);
	    gretl_matrix_set(m, r, 4 + offset, summ->sd[i]);
	    if (floateq(summ->mean[i], 0.0)) {
		cv = NADBL;
	    } else if (floateq(summ->sd[i], 0.0)) {
		cv = 0.0;
	    } else {
		cv = fabs(summ->sd[i] / summ->mean[i]);
	    }
	    gretl_matrix_set(m, r,  5 + offset, cv);
	    gretl_matrix_set(m, r,  6 + offset, summ->skew[i]);
	    gretl_matrix_set(m, r,  7 + offset, summ->xkurt[i]);
	    gretl_matrix_set(m, r,  8 + offset, summ->perc05[i]);
	    gretl_matrix_set(m, r,  9 + offset, summ->perc95[i]);
	    gretl_matrix_set(m, r, 10 + offset, summ->iqr[i]);
	    gretl_matrix_set(m, r, 11 + offset, summ->misscount[i]);
	}
    }
}

static void record_summary (Summary *summ, const DATASET *dset)
{
    gretl_matrix *m;
    int err = 0;

    m = make_summary_recorder(summ->list, summ->opt, dset,
			      NULL, NULL, &err);
    if (!err) {
	transcribe_summary(summ, m, NULL, 0);
	set_last_result_data(m, GRETL_TYPE_MATRIX);
    }
}

/**
 * list_summary:
 * @list: list of series to process.
 * @dset: dataset struct.
 * @opt: may include %OPT_S for "simple" version.
 * @prn: gretl printing struct.
 *
 * Prints descriptive statistics for the listed variables.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int list_summary (const int *list, int wgtvar,
		  const DATASET *dset,
		  gretlopt opt, PRN *prn)
{
    Summary *summ = NULL;
    int err = 0;

    if (list != NULL) {
	if (list[0] == 0) {
	    return 0;
	}
	if (wgtvar == 0) {
	    /* no weights */
	    summ = get_summary(list, dset, opt, prn, &err);
	} else {
	    summ = get_summary_weighted(list, dset, wgtvar,
					opt, prn, &err);
	}
    } else {
	int *tmplist = full_var_list(dset, NULL);

	if (tmplist != NULL) {
	    if (wgtvar == 0) {
		/* no weights */
		summ = get_summary(tmplist, dset, opt, prn, &err);
	    } else {
		summ = get_summary_weighted(tmplist, dset, wgtvar,
					    opt, prn, &err);
	    }
	    free(tmplist);
	}
    }

    if (summ != NULL) {
	if (!(opt & OPT_Q)) {
	    print_summary(summ, dset, prn);
	}
	record_summary(summ, dset);
	free_summary(summ);
    }

    return err;
}

int summary_statistics_by (const int *list, DATASET *dset,
			   gretlopt opt, PRN *prn)
{
    const char *byname = get_optval_string(SUMMARY, OPT_B);
    int quiet = (opt & OPT_Q);
    int byvar;
    series_table *st = NULL;
    gretl_matrix *xvals = NULL;
    gretl_matrix *m = NULL;
    int *mylist = NULL;
    const double *x;
    int i, nvals = 0;
    int single, err = 0;

    if (dset == NULL || byname == NULL) {
        return E_DATA;
    }

    /* FIXME accept "unit" and "time"/"period" in place of actual
       variables for panel data? */

    byvar = current_series_index(dset, byname);
    if (byvar < 0) {
        return E_UNKVAR;
    } else if (!accept_as_discrete(dset, byvar, 0)) {
        gretl_errmsg_sprintf(_("The variable '%s' is not discrete"), byname);
        return E_DATA;
    }

    x = (const double *) dset->Z[byvar];

    if (list == NULL) {
        /* compose full series list, minus the "by" variable */
        int pos;

        mylist = full_var_list(dset, NULL);
        if (mylist == NULL) {
            return E_ALLOC;
        }
        pos = in_gretl_list(mylist, byvar);
        if (pos > 0) {
            gretl_list_delete_at_pos(mylist, pos);
        }
        if (mylist[0] == 0) {
            free(mylist);
            return E_DATA;
        }
	list = mylist;
    }

    single = (list[0] == 1);

    xvals = gretl_matrix_values(x + dset->t1, dset->t2 - dset->t1 + 1,
                                OPT_S, &err);

    if (!err) {
        nvals = gretl_vector_get_length(xvals);
        if (nvals == 0) {
            err = E_DATA;
        } else {
            st = series_get_string_table(dset, byvar);
        }
    }

    if (!err && single) {
	m = make_summary_recorder(list, opt, dset, xvals, st, &err);
	if (!quiet) {
	    pputc(prn, '\n');
	    pprintf(prn, _("Summary statistics for %s, by value of %s"),
		    dset->varname[list[1]], byname);
	    pputc(prn, '\n');
	}
    }

    for (i=0; i<nvals && !err; i++) {
        Summary *summ = NULL;
        char genline[64];
        double xi = gretl_vector_get(xvals, i);
        double *rv = NULL;

        gretl_push_c_numeric_locale();
        sprintf(genline, "%s == %g", byname, xi);
        gretl_pop_c_numeric_locale();
        rv = generate_series(genline, dset, prn, &err);

        if (!err) {
            summ = get_summary_restricted(list, dset, rv,
                                          opt, prn, &err);
        }
        if (!err && !quiet) {
            if (i == 0) {
                pputc(prn, '\n');
            }
#if 0 /* 2024-05-29 */
            if (single) {
                bufspace(2, prn);
            }
#endif
            if (st != NULL) {
                const char *si = series_table_get_string(st, xi);

                pprintf(prn, "%s (n = %d):", si, summ->n);
            } else {
                pprintf(prn, "%s = %g (n = %d):\n", byname, xi, summ->n);
            }
            print_summary(summ, dset, prn);
        }
	if (!err && m != NULL) {
	    transcribe_summary(summ, m, &i, xi);
	}
	if (summ != NULL) {
	    free_summary(summ);
	}
        free(rv);
    }

    if (!err && m != NULL) {
	set_last_result_data(m, GRETL_TYPE_MATRIX);
    }

    gretl_matrix_free(xvals);
    if (mylist != NULL) {
        free(mylist);
    }

    return err;
}

/**
 * vmatrix_new:
 *
 * Returns: an allocated and initialized #VMatrix, or
 * %NULL on failure.
 */

VMatrix *vmatrix_new (void)
{
    VMatrix *v = malloc(sizeof *v);

    if (v != NULL) {
	v->vec = NULL;
	v->xbar = NULL;
	v->ssx = NULL;
	v->list = NULL;
	v->names = NULL;
	v->ci = 0;
	v->dim = 0;
	v->t1 = 0;
	v->t2 = 0;
	v->nmin = 0;
	v->nmax = 0;
	v->ncrit = 0;
	v->missing = 0;
    }

    return v;
}

/**
 * free_vmatrix:
 * @vmat: pointer to gretl correlation matrix struct
 *
 * Frees all resources associated with @vmat, and
 * the pointer itself.
 */

void free_vmatrix (VMatrix *vmat)
{
    if (vmat != NULL) {
	if (vmat->names != NULL) {
	    strings_array_free(vmat->names, vmat->dim);
	}
	if (vmat->vec != NULL) {
	    free(vmat->vec);
	}
	if (vmat->xbar != NULL) {
	    free(vmat->xbar);
	}
	if (vmat->ssx != NULL) {
	    free(vmat->ssx);
	}
	if (vmat->list != NULL) {
	    free(vmat->list);
	}
	free(vmat);
    }
}

enum {
    CORRMAT = 0,
    COVMAT
};

/* Compute correlation or covariance matrix, using the maximum
   available sample for each coefficient.
*/

static int max_corrcov_matrix (VMatrix *v, const DATASET *dset,
			       int flag)
{
    double **Z = dset->Z;
    int i, j, vi, vj, idx;
    int vn = v->t2 - v->t1 + 1;
    int m = v->dim;
    int nij, nmiss;

    v->nmin = vn;
    v->nmax = 0;

    for (i=0; i<m; i++) {
	vi = v->list[i+1];
	for (j=i; j<m; j++)  {
	    vj = v->list[j+1];
	    idx = ijton(i, j, m);
	    nij = vn;
	    if (i == j && flag == CORRMAT) {
		v->vec[idx] = 1.0;
	    } else {
		nmiss = 0;
		if (flag == COVMAT) {
		    v->vec[idx] = gretl_covar(v->t1, v->t2, Z[vi], Z[vj],
					      &nmiss);
		} else {
		    v->vec[idx] = gretl_corr(v->t1, v->t2, Z[vi], Z[vj],
					     &nmiss);
		}
		if (nmiss > 0) {
		    v->missing += 1;
		    nij -= nmiss;
		}
		if (nij > v->nmax) {
		    v->nmax = nij;
		}
		if (nij < v->nmin) {
		    v->nmin = nij;
		}
	    }
	}
    }

    /* We'll record an "ncrit" value if there's something resembling
       a common number of observations across the coefficients.
       Specifically, we require that the difference between the max
       and min number of observations is less than 10 percent of the
       maximum; in that case we set "ncrit" to the minimum.
    */
    if (v->missing > 0) {
	v->ncrit = 0;
	if (v->nmax > 0) {
	    double d = (v->nmax - v->nmin) / (double) v->nmax;

	    if (d < 0.10) {
		v->ncrit = v->nmin;
	    }
	}
    } else {
	v->ncrit = v->nmax;
    }

    return 0;
}

/* Compute correlation or covariance matrix, ensuring we use the same
   sample for all coefficients. We may be doing this in the context of
   the "corr" command or in the context of "pca". In the latter case
   we want to save the "xbar" and "ssx" values for later use when
   calculating the component series.
*/

static int uniform_corrcov_matrix (VMatrix *v, const DATASET *dset,
				   int ci, int flag)
{
    double **Z = dset->Z;
    double *xbar = NULL, *ssx = NULL;
    double *x;
    int m = v->dim;
    int mm = (m * (m + 1)) / 2;
    double d1, d2;
    int i, j, jmin, nij, t;
    int miss, n = 0;

    xbar = malloc(m * sizeof *xbar);
    if (xbar == NULL) {
	return E_ALLOC;
    }

    if (ci == PCA || flag == CORRMAT) {
	ssx = malloc(m * sizeof *ssx);
	if (ssx == NULL) {
	    free(xbar);
	    return E_ALLOC;
	}
    }

    for (i=0; i<m; i++) {
	xbar[i] = 0.0;
    }

    if (ssx != NULL) {
	for (i=0; i<m; i++) {
	    ssx[i] = 0.0;
	}
    }

    v->missing = 0;

    /* first pass: get sample size and sums */

    for (t=v->t1; t<=v->t2; t++) {
	miss = 0;
	for (i=0; i<m; i++) {
	    if (na(Z[v->list[i+1]][t])) {
		miss = 1;
		v->missing += 1;
		break;
	    }
	}
	if (!miss) {
	    n++;
	    for (i=0; i<m; i++) {
		xbar[i] += Z[v->list[i+1]][t];
	    }
	}
    }

    if (n < 2) {
	free(xbar);
	free(ssx);
	return E_MISSDATA;
    }

    for (i=0; i<m; i++) {
	xbar[i] /= n;
    }

    for (i=0; i<mm; i++) {
	v->vec[i] = 0.0;
    }

    /* second pass: get deviations from means and cumulate */

    for (t=v->t1; t<=v->t2; t++) {
	miss = 0;
	for (i=0; i<m; i++) {
	    x = Z[v->list[i+1]];
	    if (na(x[t])) {
		miss = 1;
		break;
	    }
	}
	if (!miss) {
	    for (i=0; i<m; i++) {
		x = Z[v->list[i+1]];
		d1 = x[t] - xbar[i];
		if (ssx != NULL) {
		    ssx[i] += d1 * d1;
		}
		jmin = (flag == COVMAT)? i : i + 1;
		for (j=jmin; j<m; j++) {
		    x = Z[v->list[j+1]];
		    nij = ijton(i, j, m);
		    d2 = x[t] - xbar[j];
		    v->vec[nij] += d1 * d2;
		}
	    }
	}
    }

    /* finalize: compute correlations or covariances */

    if (flag == CORRMAT) {
	/* correlations */
	for (i=0; i<m; i++) {
	    for (j=i; j<m; j++)  {
		nij = ijton(i, j, m);
		if (i == j) {
		    v->vec[nij] = 1.0;
		} else if (ssx[i] == 0.0 || ssx[j] == 0.0) {
		    v->vec[nij] = NADBL;
		} else {
		    v->vec[nij] /= sqrt(ssx[i] * ssx[j]);
		}
	    }
	}
    } else {
	/* covariances */
	for (i=0; i<mm; i++) {
	    v->vec[i] /= (n - 1);
	}
    }

    v->nmax = v->nmin = v->ncrit = n;

    if (ci == PCA) {
	v->xbar = xbar;
	v->ssx = ssx;
    } else {
	free(xbar);
	free(ssx);
    }

    return 0;
}

static int corr_skip_obs (int k, int n_ok, int uniform)
{
    if (uniform) {
	return n_ok < k;
    } else {
	return n_ok < 2;
    }
}

static void corrlist_adjust_sample (const int *list,
				    int *t1, int *t2,
				    gretlopt opt,
				    const DATASET *dset)
{
    int uniform = (opt & OPT_N)? 1 : 0;
    int t1min = *t1;
    int t2max = *t2;
    int k = list[0];
    int n_ok;
    int i, t;

    /* Note: In this context what counts as a "usable" observation
       depends on whether the --uniform option (OPT_N) is given. If
       so, an observation is usable only if it's "complete" (no NAs);
       otherwise it's potentially usable so long as there at least two
       non-missing values among the members of @list.
    */

    /* advance start of sample range to skip any unusable obs */
    for (t=t1min; t<t2max; t++) {
	n_ok = k;
	for (i=1; i<=k; i++) {
	    if (na(dset->Z[list[i]][t])) {
		n_ok--;
	    }
	}
	if (corr_skip_obs(k, n_ok, uniform)) {
	    t1min++;
	} else {
	    break;
	}
    }

    /* retard end of sample range to skip any unusable obs */
    for (t=t2max; t>t1min; t--) {
	n_ok = k;
	for (i=1; i<=k; i++) {
	    if (na(dset->Z[list[i]][t])) {
		n_ok--;
	    }
	}
	if (corr_skip_obs(k, n_ok, uniform)) {
	    t2max--;
	} else {
	    break;
	}
    }

    *t1 = t1min;
    *t2 = t2max;
}

/**
 * corrlist:
 * @ci: command index.
 * @list: list of variables to process, by ID number.
 * @dset: dataset struct.
 * @opt: option flags.
 * @err: location to receive error code.
 *
 * Computes pairwise correlation coefficients for the variables
 * specified in @list, skipping any constants.  If the option
 * flags contain OPT_N, a uniform sample is ensured: only those
 * observations for which all the listed variables have valid
 * values are used.  If OPT_C is included, we actually calculate
 * covariances rather than correlations.
 *
 * Returns: gretl correlation matrix struct, or NULL on failure.
 */

VMatrix *corrlist (int ci, int *list, const DATASET *dset,
		   gretlopt opt, int *err)
{
    int flag = (opt & OPT_C)? COVMAT : CORRMAT;
    VMatrix *v = NULL;
    int i, m, mm;

    if (list == NULL) {
	fprintf(stderr, "corrlist: list is NULL\n");
	*err = E_DATA;
	return NULL;
    }

    v = vmatrix_new();
    if (v == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* drop any constants from list */
    for (i=list[0]; i>0; i--) {
	if (gretl_isconst(dset->t1, dset->t2, dset->Z[list[i]])) {
	    gretl_list_delete_at_pos(list, i);
	}
    }

    if (list[0] < 2) {
	gretl_errmsg_set(_("corr: needs at least two non-constant arguments"));
	*err = E_DATA;
	goto bailout;
    }

    v->t1 = dset->t1;
    v->t2 = dset->t2;

    /* adjust sample range if need be */
    corrlist_adjust_sample(list, &v->t1, &v->t2, opt, dset);

    if (v->t2 - v->t1 + 1 < 3) {
	*err = E_TOOFEW;
	goto bailout;
    }

    v->dim = m = list[0];
    mm = (m * (m + 1)) / 2;

    v->names = strings_array_new(m);
    v->vec = malloc(mm * sizeof *v->vec);
    v->list = gretl_list_copy(list);

    if (v->names == NULL || v->vec == NULL || v->list == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    if (ci == PCA) {
	opt |= OPT_N;
    }

    if (opt & OPT_N) {
	/* impose uniform sample size */
	*err = uniform_corrcov_matrix(v, dset, ci, flag);
    } else {
	/* sample sizes may differ */
	*err = max_corrcov_matrix(v, dset, flag);
    }

    if (!*err) {
	for (i=0; i<m; i++) {
	    v->names[i] = gretl_strdup(dset->varname[list[i+1]]);
	    if (v->names[i] == NULL) {
		*err = E_ALLOC;
		goto bailout;
	    }
	}
    }

    v->ci = (flag == CORRMAT)? CORR : 0; /* FIXME? */

 bailout:

    if (*err) {
	free_vmatrix(v);
	v = NULL;
    }

    return v;
}

static void printcorr (const VMatrix *v, PRN *prn)
{
    double r = v->vec[1];

    pprintf(prn, "\ncorr(%s, %s)", v->names[0], v->names[1]);

    if (na(r)) {
	pprintf(prn, ": %s\n\n", _("undefined"));
    } else {
	pprintf(prn, " = %.8f\n", r);
	if (fabs(r) < 1.0) {
	    int nm2 = v->ncrit - 2;
	    double tval = r * sqrt(nm2 / (1 - r*r));

	    pputs(prn, _("Under the null hypothesis of no correlation:\n "));
	    pprintf(prn, _("t(%d) = %g, with two-tailed p-value %.4f\n"),
		    nm2, tval, student_pvalue_2(nm2, tval));
	    pputc(prn, '\n');
	} else {
	    pprintf(prn, _("Two-tailed critical values for n = %d"), v->ncrit);
	    pprintf(prn, ": 5%% %.4f, 1%% %.4f\n", rhocrit(v->ncrit, 0.05),
		    rhocrit(v->ncrit, 0.01));
	    pputc(prn, '\n');
	}
    }
}

/**
 * print_corrmat:
 * @corr: gretl correlation matrix.
 * @dset: dataset information.
 * @prn: gretl printing struct.
 *
 * Prints a gretl correlation matrix to @prn.
 */

void print_corrmat (VMatrix *corr, const DATASET *dset, PRN *prn)
{
    if (corr->dim == 2) {
	printcorr(corr, prn);
    } else {
	char date1[OBSLEN], date2[OBSLEN];
	gchar *tmp = NULL;

	ntolabel(date1, corr->t1, dset);
	ntolabel(date2, corr->t2, dset);

	pputc(prn, '\n');

	if (corr->nmin == corr->t2 - corr->t1 + 1) {
	    tmp = g_strdup_printf(_("%s, using the %d observations %s - %s"),
				  _("Correlation Coefficients"), corr->ncrit,
				  date1, date2);
	} else if (corr->nmin == corr->nmax) {
	    tmp = g_strdup_printf(_("%s, using %d observations from %s - %s"),
				  _("Correlation Coefficients"), corr->ncrit,
				  date1, date2);
	} else {
	    tmp = g_strdup_printf(_("%s, using samples of size %d to %d"),
				  _("Correlation Coefficients"),
				  corr->nmin, corr->nmax);
	}
	output_line(tmp, prn, 0);
	g_free(tmp);

	if (corr->missing > 0) {
	    output_line(_("(missing values were skipped)"), prn, 1);
	}

	if (corr->ncrit > 0) {
	    pprintf(prn, _("Two-tailed critical values for n = %d"), corr->ncrit);
	    pprintf(prn, ": 5%% %.4f, 1%% %.4f\n", rhocrit(corr->ncrit, 0.05),
		    rhocrit(corr->ncrit, 0.01));
	    pputc(prn, '\n');
	}

	text_print_vmatrix(corr, prn);
    }
}

static void record_corr_matrix (VMatrix *c)
{
    gretl_matrix *m = NULL;
    int i, j, k, n = c->dim;
    double cij;

    m = gretl_matrix_alloc(n, n);

    if (m != NULL) {
	k = 0;
	for (i=0; i<n; i++) {
	    for (j=i; j<n; j++) {
		cij = c->vec[k];
		gretl_matrix_set(m, i, j, cij);
		gretl_matrix_set(m, j, i, cij);
		k++;
	    }
	}

	if (c->names) {
	    char **rlab = strings_array_new(n);
	    char **clab = strings_array_new(n);

	    if (rlab != NULL && clab != NULL) {
		for (i=0; i<n; i++) {
		    rlab[i] = gretl_strdup(c->names[i]);
		    clab[i] = gretl_strdup(c->names[i]);
		}
		gretl_matrix_set_rownames(m, rlab);
		gretl_matrix_set_colnames(m, clab);
	    }
	}

	set_last_result_data(m, GRETL_TYPE_MATRIX);
    }
}

/**
 * gretl_corrmx:
 * @list: gives the ID numbers of the variables to process.
 * @dset: dataset struct.
 * @opt: option flags: OPT_N means use uniform sample size,
 * OPT_U controls plotting options.
 * @prn: gretl printing struct.
 *
 * Computes and prints the correlation matrix for the specified list
 * of variables.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int gretl_corrmx (int *list, const DATASET *dset,
		  gretlopt opt, PRN *prn)
{
    VMatrix *corr = NULL;
    int err = 0;

    if (list != NULL) {
	if (list[0] == 0) {
	    return 0;
	} else {
	    corr = corrlist(CORR, list, dset, opt, &err);
	}
    } else {
	int *tmplist = full_var_list(dset, NULL);

	if (tmplist != NULL) {
	    corr = corrlist(CORR, tmplist, dset, opt, &err);
	    free(tmplist);
	}
    }

    if (corr != NULL) {
	if (!(opt & OPT_Q)) {
	    print_corrmat(corr, dset, prn);
	}
	if (corr->dim > 2 && gnuplot_graph_wanted(PLOT_HEATMAP, opt, NULL)) {
	    err = plot_corrmat(corr, opt);
	}

	record_corr_matrix(corr);
	free_vmatrix(corr);
    }

    return err;
}

/*
  Satterthwaite, F.E. (1946). "An Approximate Distribution of Estimates
  of Variance Components". Biometrics Bulletin, 2, 6, pp. 110–114.

  v = (s^2_1/n_1 + s^2_2/n_2)^2 /
        [(s^2_1/n1)^2 / (n_1 - 1) + (s^2_2/n_2)^2 / (n_2 - 1)]
*/

int satterthwaite_df (double v1, int n1,
		      double v2, int n2)
{
    double rtnum = v1 / n1 + v2 / n2;
    double d1 = (v1/n1) * (v1/n1) / (n1 - 1);
    double d2 = (v2/n2) * (v2/n2) / (n2 - 1);
    double v = (rtnum * rtnum) / (d1 + d2);

    return (int) floor(v);
}

static int test_data_from_dummy (const DATASET *dset,
                                 int v1, int v2,
                                 double **px,
                                 double **py,
                                 int *pn1, int *pn2)
{
    double *x = NULL;
    double *y = NULL;
    int n1 = 0;
    int n2 = 0;
    int t, ix, iy;

    if (!gretl_isdummy(dset->t1, dset->t2, dset->Z[v2])) {
        return E_INVARG;
    }

    for (t=dset->t1; t<=dset->t2; t++) {
        if (na(dset->Z[v1][t])) {
            continue;
        } else if (dset->Z[v2][t] == 0.0) {
            n1++;
        } else if (dset->Z[v2][t] == 1.0) {
            n2++;
        }
    }
    if (n1 < 2 || n2 < 2) {
        return E_TOOFEW;
    }

    x = malloc(n1 * sizeof *x);
    y = malloc(n2 * sizeof *y);
    if (x == NULL || y == NULL) {
        free(x); free(y);
        return E_ALLOC;
    }

    ix = iy = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
        if (na(dset->Z[v1][t])) {
            continue;
        } else if (dset->Z[v2][t] == 0.0) {
            x[ix++] = dset->Z[v1][t];
        } else if (dset->Z[v2][t] == 1.0) {
            y[iy++] = dset->Z[v1][t];
        }
    }

    *px = x;
    *py = y;
    *pn1 = n1;
    *pn2 = n2;

    return 0;
}

static double *transcribe_ok_obs (const DATASET *dset, int v,
                                  int *pn, int *err)
{
    const double *src = dset->Z[v];
    double *ret = NULL;
    int t, n = 0;

    for (t=dset->t1; t<=dset->t2; t++) {
        if (!na(src[t])) {
            n++;
        }
    }

    if (n < 2) {
        *err = E_TOOFEW;
    } else {
        ret = malloc(n * sizeof *ret);
        if (ret == NULL) {
            *err = E_ALLOC;
        }
    }
    if (!*err) {
        *pn = n;
        n = 0;
        for (t=dset->t1; t<=dset->t2; t++) {
            if (!na(src[t])) {
                ret[n++] = src[t];
            }
        }
    }

    return ret;
}

static double *transcribe_ok_pairs (const DATASET *dset,
                                    int v1, int v2,
                                    int *pn, int *err)
{
    const double *x = dset->Z[v1];
    const double *y = dset->Z[v2];
    double *ret = NULL;
    int t, n = 0;

    for (t=dset->t1; t<=dset->t2; t++) {
        if (!na(x[t]) && !na(y[t])) {
            n++;
        }
    }

    if (n < 2) {
        *err = E_TOOFEW;
    } else {
        ret = malloc(n * sizeof *ret);
        if (ret == NULL) {
            *err = E_ALLOC;
        }
    }
    if (!*err) {
        *pn = n;
        n = 0;
        for (t=dset->t1; t<=dset->t2; t++) {
            if (!na(x[t]) && !na(y[t])) {
                ret[n++] = x[t] - y[t];
            }
        }
    }

    return ret;
}

/* Note: @y may be NULL; in that case @x holds paired differences, and
   @n2 should be 0.
*/

static int ols_means_test (const double *x, int n1,
                           const double *y, int n2,
                           double *m1, double *m2,
                           double *se, int *df,
                           PRN *prn)
{
    MODEL mod;
    DATASET *aset;
    int list[4] = {3, 1, 0, 2};
    int nobs = n1 + n2;
    int i, err = 0;

    if (y == NULL) {
        /* using paired data */
        list[0] = 2;
    }

    aset = create_auxiliary_dataset(list[0], nobs, OPT_NONE);
    if (aset == NULL) {
        return E_ALLOC;
    }

    /* pairs, or stacked dependent variable, @x | @y */
    memcpy(aset->Z[1], x, n1 * sizeof *x);
    if (y != NULL) {
        memcpy(aset->Z[1] + n1, y, n2 * sizeof *y);
    }
    strcpy(aset->varname[1], "x_y");

    if (y != NULL) {
        /* dummy to indicate @y observations */
        for (i=0; i<nobs; i++) {
            aset->Z[2][i] = i >= n1;
        }
        strcpy(aset->varname[2], "dum_y");
    }

    gretl_model_init(&mod, aset);
    mod = lsq(list, aset, OLS, OPT_A | OPT_R);
    err = mod.errcode;

    if (!err) {
        if (prn != NULL) {
            printmodel(&mod, aset, OPT_S, prn);
        }
        *m1 = mod.coeff[0];
        if (y != NULL) {
            *m2 = mod.coeff[0] + mod.coeff[1];
            *se = mod.sderr[1];
        } else {
            *se = mod.sderr[0];
        }
        *df = mod.dfd;
    }

    clear_model(&mod);
    destroy_dataset(aset);

    return err;
}

static int test_data_from_list (const DATASET *dset,
                                int v1, int v2,
                                double **px,
                                double **py,
                                int *pn1, int *pn2)
{
    double *x, *y;
    int n1, n2;
    int err = 0;

    x = transcribe_ok_obs(dset, v1, &n1, &err);
    if (!err) {
        y = transcribe_ok_obs(dset, v2, &n2, &err);
    }
    if (!err) {
        *px = x;
        *py = y;
        *pn1 = n1;
        *pn2 = n2;
    }

    return err;
}

static int paired_data_from_list (const DATASET *dset,
                                  int v1, int v2,
                                  double **px, int *pn)
{
    double *x;
    int n, err = 0;

    x = transcribe_ok_pairs(dset, v1, v2, &n, &err);
    if (!err) {
        *px = x;
        *pn = n;
    }

    return err;
}

static void paired_diff_print (int n,
                               double mdiff,
                               double se,
                               int df,
                               double t,
                               double pval,
                               PRN *prn)
{
    pputc(prn, '\n');
    pputs(prn, _("Paired difference test"));
    pputs(prn, "\n\n");
    pputs(prn, "   ");
    pprintf(prn, _("Number of observations = %d\n"), n);
    pputs(prn, "   ");
    pprintf(prn, _("Sample mean paired difference = %g"), mdiff);
    pputc(prn, '\n');
    pputs(prn, _("   Null hypothesis: The population mean difference is zero."));
    pputc(prn, '\n');
    pprintf(prn, _("   Estimated standard error = %g\n"), se);
    pprintf(prn, _("   Test statistic: t(%d) = %g\n"), df, t);
    pprintf(prn, _("   p-value (two-tailed) = %g\n\n"), pval);
}

/**
 * means_test:
 * @list: gives the ID numbers of the variables to compare.
 * @dset: dataset struct.
 * @opt: if OPT_O, assume that the population variances differ;
 * if OPT_D interpret the second series as a factor (dummy);
 * if OPT_R, carry out a robust, regression-based test; if OPT_P,
 * carry out a paired-difference test; if OPT_Q, don't print
 * the resuults.
 * @prn: gretl printing struct.
 *
 * Carries out test of the null hypothesis that the means of two
 * variables are equal.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int means_test (const int *list, const DATASET *dset,
		gretlopt opt, PRN *prn)
{
    double m1, m2, s1, s2, se, mdiff, t, pval;
    double *x = NULL, *y = NULL;
    int vardiff = (opt & OPT_O)? 1 : 0;
    int usedum = (opt & OPT_D)? 1 : 0;
    int robust = (opt & OPT_R)? 1 : 0;
    int paired = (opt & OPT_P)? 1 : 0;
    int v1, v2;
    int n1 = 0;
    int n2 = 0;
    int df, err = 0;

    if (list[0] < 2) {
	err = E_ARGS;
    } else {
        err = incompatible_options(opt, OPT_O | OPT_P);
    }
    if (err) {
        return err;
    }

    v1 = list[1];
    v2 = list[2];

    if (paired) {
        /* get @x and @y as portions of dset->Z[v1] */
        err = paired_data_from_list(dset, v1, v2, &x, &n1);
    } else if (usedum) {
        /* get @x and @y as portions of dset->Z[v1] */
        err = test_data_from_dummy(dset, v1, v2, &x, &y, &n1, &n2);
    } else {
        /* get @x and @y as dset->Z[v1], dset->Z[v2] */
        err = test_data_from_list(dset, v1, v2, &x, &y, &n1, &n2);
    }

    if (err) {
        free(x); free(y);
        return err;
    }

    if (robust || paired) {
        /* use regression method */
        err = ols_means_test(x, n1, y, n2, &m1, &m2, &se, &df, prn);
        if (err) {
            goto bailout;
        }
        mdiff = paired ? m1 : m1 - m2;
    } else {
        double var1, var2;

        df = n1 + n2 - 2;
        gretl_moments(0, n1-1, x, NULL, &m1, &s1, NULL, NULL, 1);
        gretl_moments(0, n2-1, y, NULL, &m2, &s2, NULL, NULL, 1);
        mdiff = m1 - m2;
        var1 = s1 * s1;
        var2 = s2 * s2;

        if (vardiff) {
            /* assuming unequal variances */
            se = sqrt((var1 / n1) + (var2 / n2));
            df = satterthwaite_df(var1, n1, var2, n2);
        } else {
            /* form pooled estimate of variance */
            double s2p;

            s2p = ((n1-1) * var1 + (n2-1) * var2) / df;
            se = sqrt(s2p / n1 + s2p / n2);
            df = n1 + n2 - 2;
        }
    }

    t = mdiff / se;
    pval = student_pvalue_2(df, t);
    record_test_result(t, pval);

    if (opt & OPT_Q) {
        /* quiet */
        goto bailout;
    }

    if (paired) {
        paired_diff_print(n1, mdiff, se, df, t, pval, prn);
    } else {
        if (robust) {
            pputs(prn, _("Robust equality of means test"));
            pputs(prn, "\n\n");
        } else {
            pprintf(prn, _("\nEquality of means test "
                           "(assuming %s variances)\n\n"),
                    vardiff ? _("unequal") : _("equal"));
        }
        if (usedum) {
            pprintf(prn, "   %s (%s==0): ", dset->varname[v1], dset->varname[v2]);
            pprintf(prn, _("Number of observations = %d\n"), n1);
            pprintf(prn, "   %s (%s==1): ", dset->varname[v1], dset->varname[v2]);
            pprintf(prn, _("Number of observations = %d\n"), n2);
        } else {
            pprintf(prn, "   %s: ", dset->varname[v1]);
            pprintf(prn, _("Number of observations = %d\n"), n1);
            pprintf(prn, "   %s: ", dset->varname[v2]);
            pprintf(prn, _("Number of observations = %d\n"), n2);
        }
        pprintf(prn, _("   Difference between sample means = %g - %g = %g\n"),
                m1, m2, mdiff);
        pputs(prn, _("   Null hypothesis: The two population means are the same.\n"));
        pprintf(prn, _("   Estimated standard error = %g\n"), se);
        pprintf(prn, _("   Test statistic: t(%d) = %g\n"), df, t);
        pprintf(prn, _("   p-value (two-tailed) = %g\n\n"), pval);
    }

 bailout:

    free(x);
    free(y);

    return 0;
}

static double trimmed_mean (const double *x, double p,
                            int n, int *err)
{
    double ret;
    double *sx;
    int skip;

    sx = malloc(n * sizeof *sx);
    memcpy(sx, x, n * sizeof *sx);
    qsort(sx, n, sizeof *sx, gretl_compare_doubles);
    skip = floor(p * n);
    ret = gretl_mean(skip, n-skip-1, sx);
    free(sx);

    return ret;
}

static double levene_center (const double *x, int n,
                             int f, double p,
                             int *err)
{
    if (f == 1) {
        return gretl_mean(0, n-1, x);
    } else if (f == 2) {
        return gretl_median(0, n-1, x);
    } else {
        return trimmed_mean(x, p, n, err);
    }
}

/* Levene's test for difference of variance, with the options added by
   Brown and Forsythe. Refs:

   Levene, H. (1960) Robust Tests for Equality of Variances. In: Olkin,
   I. (ed), Contributions to Probability and Statistics, Stanford
   University Press. pp. 278-292.

   Brown, Morton B and Forsythe, Alan B. (1974). "Robust tests for the
   equality of variances". Journal of the American Statistical
   Association 69 (346), pp. 364–367.
*/

static double levene_test (const double *x, int n1,
                           const double *y, int n2,
                           int f, double p,
                           int *err)
{
    double zxbar, zybar, zbar, lc;
    double num, den, W;
    double d1, d2;
    double *zx, *zy;
    int i, dfd;

    zx = malloc(n1 * sizeof *zx);
    zy = malloc(n2 * sizeof *zy);
    if (zx == NULL || zy == NULL) {
        *err = E_ALLOC;
        return NADBL;
    }

    /* absolute deviations, x */
    lc = levene_center(x, n1, f, p, err);
    zxbar = 0;
    for (i=0; i<n1; i++) {
        zx[i] = fabs(x[i] - lc);
        zxbar += zx[i];
    }
    zxbar /= n1;

    /* absolute deviations, y */
    lc = levene_center(y, n2, f, p, err);
    zybar = 0;
    for (i=0; i<n2; i++) {
        zy[i] = fabs(y[i] - lc);
        zybar += zy[i];
    }
    zybar /= n2;

    /* grand mean */
    zbar = (n1 * zxbar + n2 * zybar) / (n1+n2);
    dfd = n1 + n2 - 2;

    /* numerator of W */
    d1 = zxbar - zbar;
    d2 = zybar - zbar;
    num = n1 * d1 * d1 + n2 * d2 * d2;

    /* denominator of W */
    den = 0;
    for (i=0; i<n1; i++) {
        d1 = zx[i] - zxbar;
        den += d1 * d1;
    }
    for (i=0; i<n2; i++) {
        d2 = zy[i] - zybar;
        den += d2 * d2;
    }

    W = dfd * num / den;

    free(zx);
    free(zy);

    return W;
}

static int handle_levene_options (double *p, int *err)
{
    const char *s = get_optval_string(VARTEST, OPT_R);
    double default_trim = 0.10;
    int f = 1;

    if (s != NULL) {
        if (!strcmp(s, "mean")) {
            ; /* OK */
        } else if (!strcmp(s, "median")) {
            f = 2;
        } else if (!strcmp(s, "trimmed")) {
            f = 3;
            *p = default_trim;
        } else if (!strncmp(s, "trimmed,", 8)) {
            int k = positive_int_from_string(s + 8);

            if (k > 0 && k < 50) {
                *p = k / 100.0;
                f = 3;
            } else {
                *err = E_INVARG;
            }
        }
    }

    return f;
}

static void levene_print (double W, double pval,
                          int dfn, int dfd, int f,
                          double p, PRN *prn)
{
    pputs(prn, _("\nRobust equality of variances test\n"));
    pputs(prn, _("Centering: "));
    if (f == 1) {
        pputs(prn, _("mean"));
    } else if (f == 2) {
        pputs(prn, _("median"));
    } else {
        pprintf(prn, _("%d percent trimmed mean"), (int) (p*100));
    }
    pputc(prn, '\n');
    pprintf(prn, "F(%d,%d) = %g ", dfn, dfd, W);
    pprintf(prn, "%s %g\n", _("with p-value"), pval);
}

static void vartest_print (double F, double pval,
                           int dfn, int dfd,
                           const DATASET *dset,
                           int v1, int v2,
                           int n1, int n2,
                           int usedum, PRN *prn)
{
    pputs(prn, _("\nEquality of variances test\n\n"));
    if (usedum) {
        pprintf(prn, "   %s (%s==0): ", dset->varname[v1], dset->varname[v2]);
        pprintf(prn, _("Number of observations = %d\n"), n1);
        pprintf(prn, "   %s (%s==1): ", dset->varname[v1], dset->varname[v2]);
        pprintf(prn, _("Number of observations = %d\n"), n2);
    } else {
        pprintf(prn, "   %s: ", dset->varname[v1]);
        pprintf(prn, _("Number of observations = %d\n"), n1);
        pprintf(prn, "   %s: ", dset->varname[v2]);
        pprintf(prn, _("Number of observations = %d\n"), n2);
    }
    pprintf(prn, _("   Ratio of sample variances = %g\n"), F);
    pprintf(prn, "   %s: %s\n", _("Null hypothesis"),
            _("The two population variances are equal"));
    pprintf(prn, "   %s: F(%d,%d) = %g\n", _("Test statistic"), dfn, dfd, F);
    pprintf(prn, _("   p-value (two-tailed) = %g\n\n"), pval);
}

/**
 * vars_test:
 * @list: gives the ID numbers of the variables to compare.
 * @dset: pointer to dataset.
 * @opt: if includes OPT_R, carry out a robust (Levene-type) test;
 * if OPT_D, interpret the second series as a factor (dummy); if
 * OPT_Q (quiet), don't print the results.
 * @prn: gretl printing struct.
 *
 * Carries out test of the null hypothesis that the variances of two
 * variables are equal.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int vars_test (const int *list, const DATASET *dset,
               gretlopt opt, PRN *prn)
{
    double F, pval;
    double p = NADBL;
    double *x = NULL;
    double *y = NULL;
    int robust = (opt & OPT_R)? 1 : 0;
    int usedum = (opt & OPT_D)? 1 : 0;
    int dfn, dfd;
    int v1, v2;
    int n1 = 0;
    int n2 = 0;
    int f = 0;
    int err = 0;

    if (list[0] < 2) {
        return E_ARGS;
    }

    if (robust) {
        f = handle_levene_options(&p, &err);
        if (err) {
            return err;
        }
    }

    v1 = list[1];
    v2 = list[2];

    if (usedum) {
        /* get @x and @y as portions of dset->Z[v1] */
        err = test_data_from_dummy(dset, v1, v2, &x, &y, &n1, &n2);
    } else {
        /* get @x and @y as dset->Z[v1], dset->Z[v2] */
        err = test_data_from_list(dset, v1, v2, &x, &y, &n1, &n2);
    }

    if (err) {
        free(x); free(y);
        return err;
    }

    if (robust) {
        F = levene_test(x, n1, y, n2, f, p, &err);
        dfn = 1; /* = 2 - 1 */
        dfd = n1 + n2 - 2;
    } else {
        double mx, my, sx, sy;

        gretl_moments(0, n1-1, x, NULL, &mx, &sx, NULL, NULL, 1);
        gretl_moments(0, n2-1, y, NULL, &my, &sy, NULL, NULL, 1);
        F = (sx * sx) / (sy * sy);
        dfn = n1 - 1;
        dfd = n2 - 1;
    }

    pval = snedecor_cdf_comp(dfn, dfd, F);
    record_test_result(F, pval);

    if (!(opt & OPT_Q)) {
        /* not quiet */
        if (robust) {
            levene_print(F, pval, dfn, dfd, f, p, prn);
        } else {
            vartest_print(F, pval, dfn, dfd, dset, v1, v2,
                          n1, n2, usedum, prn);
        }
    }

    free(x);
    free(y);

    return err;
}

struct MahalDist_ {
    int *list;
    int n;
    gretl_matrix *d;
};

const double *mahal_dist_get_distances (const MahalDist *md)
{
    return md->d->val;
}

int mahal_dist_get_n (const MahalDist *md)
{
    return md->n;
}

const int *mahal_dist_get_varlist(const MahalDist *md)
{
    return md->list;
}

void free_mahal_dist (MahalDist *md)
{
    if (md != NULL) {
	free(md->list);
	gretl_matrix_free(md->d);
	free(md);
    }
}

static MahalDist *mahal_dist_new (const int *list, int n)
{
    MahalDist *md = malloc(sizeof *md);

    if (md != NULL) {
	md->d = gretl_column_vector_alloc(n);
	if (md->d == NULL) {
	    free(md);
	    md = NULL;
	} else {
	    md->list = gretl_list_copy(list);
	    if (md->list == NULL) {
		gretl_matrix_free(md->d);
		free(md);
		md = NULL;
	    } else {
		md->n = n;
	    }
	}
    }

    if (md != NULL) {
	int t;

	for (t=0; t<n; t++) {
	    md->d->val[t] = NADBL;
	}
    }

    return md;
}

static int mdist_saver (DATASET *dset)
{
    int t, v, err;

    err = dataset_add_series(dset, 1);

    if (err) {
	v = 0;
    } else {
	v = dset->v - 1;

	for (t=0; t<dset->n; t++) {
	    dset->Z[v][t] = NADBL;
	}

	strcpy(dset->varname[v], "mdist");
	make_varname_unique(dset->varname[v], v, dset);
	series_set_label(dset, v, _("Mahalanobis distances"));
    }

    return v;
}

static int
real_mahalanobis_distance (const int *list, DATASET *dset,
			   gretlopt opt, MahalDist *md,
			   PRN *prn)
{
    gretl_matrix *S = NULL;
    gretl_vector *means = NULL;
    gretl_vector *xdiff;
    int orig_t1 = dset->t1;
    int orig_t2 = dset->t2;
    int n, err = 0;

    list_adjust_sample(list, &dset->t1, &dset->t2, dset, NULL);

    n = sample_size(dset);
    if (n < 2) {
	err = E_DATA;
	goto bailout;
    }

    xdiff = gretl_column_vector_alloc(list[0]);
    if (xdiff == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    S = gretl_covariance_matrix_from_varlist(list,
					     dset,
					     &means,
					     &err);

    if (!err) {
	if (opt & OPT_V) {
	    gretl_matrix_print_to_prn(S, _("Covariance matrix"), prn);
	}
	err = gretl_invert_symmetric_matrix(S);
	if (err) {
	    fprintf(stderr, "error inverting covariance matrix\n");
	} else if (opt & OPT_V) {
	    gretl_matrix_print_to_prn(S, _("Inverse of covariance matrix"), prn);
	}
    }

    if (!err) {
	int k = gretl_vector_get_length(means);
	int miss, obslen = 0, savevar = 0;
	double m, x, xbar;
	int i, t, vi;

	if (opt & OPT_S) {
	    /* save the results to a data series */
	    savevar = mdist_saver(dset);
	}

	if (prn != NULL) {
	    pprintf(prn, "%s\n", _("Mahalanobis distances from the centroid"));
	    pprintf(prn, "%s\n", _("using the variables:"));
	    for (i=1; i<=list[0]; i++) {
		pprintf(prn, " %s\n", dset->varname[list[i]]);
	    }
	    pputc(prn, '\n');

	    obslen = max_obs_marker_length(dset);
	    if (obslen < 8) {
		obslen = 8;
	    }
	}

	for (t=dset->t1; t<=dset->t2; t++) {
	    miss = 0;

	    /* write vector of deviations from centroid for
	       observation t */
	    for (i=0; i<k; i++) {
		vi = list[i+1];
		xbar = gretl_vector_get(means, i);
		x = dset->Z[vi][t];
		miss += na(x);
		if (!miss) {
		    gretl_vector_set(xdiff, i,  x - xbar);
		}
	    }

	    if (miss) {
		m = NADBL;
	    } else {
		m = gretl_scalar_qform(xdiff, S, &err);
		if (!err) {
		    m = sqrt(m);
		}
	    }

	    if (prn != NULL) {
		print_obs_marker(t, dset, obslen, prn);
		if (err || miss) {
		    pprintf(prn, "NA\n");
		} else {
		    pprintf(prn, "%9.6f\n", m);
		}
	    }

	    if (savevar > 0) {
		dset->Z[savevar][t] = m;
	    } else if (md != NULL) {
		md->d->val[t] = m;
	    }
	}

	if (savevar > 0 && prn != NULL) {
	    pputc(prn, '\n');
	    pprintf(prn, _("Distances saved as '%s'"),
		    dset->varname[savevar]);
	    pputc(prn, '\n');
	}

	if (prn != NULL) {
	    pputc(prn, '\n');
	}
    }

    gretl_matrix_free(xdiff);
    gretl_matrix_free(means);
    gretl_matrix_free(S);

 bailout:

    dset->t1 = orig_t1;
    dset->t2 = orig_t2;

    return err;
}

int mahalanobis_distance (const int *list, DATASET *dset,
			  gretlopt opt, PRN *prn)
{
    return real_mahalanobis_distance(list, dset, opt, NULL,
				     (opt & OPT_Q) ? NULL : prn);
}

MahalDist *get_mahal_distances (const int *list, DATASET *dset,
				gretlopt opt, PRN *prn, int *err)
{
    MahalDist *md = mahal_dist_new(list, dset->n);

    if (md == NULL) {
	*err = E_ALLOC;
    } else {
	*err = real_mahalanobis_distance(list, dset, opt,
					 md, prn);
	if (*err) {
	    free_mahal_dist(md);
	    md = NULL;
	}
    }

    return md;
}

enum {
    EUCLIDEAN,
    MANHATTAN,
    HAMMING,
    CHEBYSHEV,
    COSINE,
    MAHALANOBIS
};

static int distance_type (const char *s)
{
    int n = strlen(s);

    if (n == 0 || (n == 1 && (*s == 'c' || *s == 'm'))) {
	return -1;
    }

    if (!strncmp(s, "euclidean", n)) {
	return EUCLIDEAN;
    } else if (!strncmp(s, "manhattan", n)) {
	return MANHATTAN;
    } else if (!strncmp(s, "hamming", n)) {
	return HAMMING;
    } else if (!strncmp(s, "chebyshev", n) ||
	       !strncmp(s, "chebychev", n)) {
	return CHEBYSHEV;
    } else if (!strncmp(s, "cosine", n)) {
	return COSINE;
    } else if (!strncmp(s, "mahalanobis", n)) {
	return MAHALANOBIS;
    } else {
	return -1;
    }
}

/* Convert from matrix @X supplied to distance() to
   dataset, to use the Mahalanobis distance code
   above.
*/

static gretl_matrix *mahal_bridge (const gretl_matrix *X,
				   const gretl_matrix *Y,
				   int *err)
{
    gretl_matrix *d = NULL;
    MahalDist *md = NULL;
    DATASET *dset = NULL;

    if (Y != NULL) {
	/* we're not going to handle a second matrix */
	*err = E_INVARG;
	return NULL;
    }

    dset = gretl_dataset_from_matrix(X, NULL, OPT_B, err);
    if (!*err) {
	int *list = gretl_consecutive_list_new(1, X->cols);

	md = get_mahal_distances(list, dset, OPT_NONE, NULL, err);
	free(list);
    }
    if (!*err) {
	d = md->d;
	md->d = NULL;
    }

    destroy_dataset(dset);
    free_mahal_dist(md);

    return d;
}

/* distance(): produces a set of pairwise distances, either between
   the rows of @X, or if @Y is non-NULL, between the rows of @X and
   the rows of @Y.

   If @X is m x n and @Y is NULL, the result is a vector &v of length
   m * (m - 1) / 2, comprising the non-redundant pairwise distances
   with the zeros on the diagonal suppressed. This can be turned into
   a square matrix via unvech(v, 0).

   If @X is m x n and @Y is p x n, the result is a matrix with m rows
   and p columns.

   If @X and @Y have different column dimensions an error is flagged.
*/

gretl_matrix *distance (const gretl_matrix *X,
			const gretl_matrix *Y,
			const char *type, int *err)
{
    gretl_matrix *ret;
    double d, dij, x, y;
    double den1, den2;
    int dtype, vlen, pos;
    int r1, r2, c, jmin;
    int i, j, k;
    int nothirdarg = (Y == NULL);

    if (gretl_is_null_matrix(X) || gretl_is_complex(X)) {
	*err = E_INVARG;
	return NULL;
    }

    if (Y != NULL && (gretl_is_complex(Y) || Y->cols != X->cols)) {
	*err = E_INVARG;
	return NULL;
    }

    dtype = distance_type(type);
    if (dtype < 0) {
	*err = E_INVARG;
	return NULL;
    }

    if (dtype == MAHALANOBIS) {
	return mahal_bridge(X, Y, err);
    }

    r1 = X->rows;
    c  = X->cols;

    if (nothirdarg) {
	r2 = r1;
	vlen = r1 * (r1 - 1) / 2;
	/* column vector for results */
	ret = gretl_matrix_alloc(vlen, 1);
    } else {
	r2 = Y->rows;
	/* matrix for results */
	ret = gretl_matrix_alloc(r1, r2);
    }

    if (ret == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    pos = 0;
    for (i=0; i<r1; i++) {
	jmin = nothirdarg ? i + 1 : 0;
	for (j=jmin; j<r2; j++) {
	    den1 = den2 = dij = 0;
	    for (k=0; k<c; k++) {
		x = gretl_matrix_get(X, i, k);
		if (Y == NULL) {
		    y = gretl_matrix_get(X, j, k);
		} else {
		    y = gretl_matrix_get(Y, j, k);
		}
		if (dtype == MANHATTAN) {
		    dij += fabs(x - y);
		} else if (dtype == HAMMING) {
		    dij += (x - y != 0);
		} else if (dtype == CHEBYSHEV) {
		    d = fabs(x - y);
		    if (d > dij) {
			dij = d;
		    }
		} else if (dtype == COSINE) {
		    dij += x * y;
		    den1 += x * x;
		    den2 += y * y;
		} else {
		    /* euclidean */
		    d = x - y;
		    dij += d * d;
		}
	    }
	    if (dtype == EUCLIDEAN) {
		dij = sqrt(dij);
	    } else if (dtype == HAMMING) {
		dij /= c;
	    } else if (dtype == COSINE) {
		dij = 1.0 - dij / sqrt(den1 * den2);
	    }
	    if (nothirdarg) {
		ret->val[pos++] = dij;
	    } else {
		gretl_matrix_set(ret, i, j, dij);
	    }
	}
    }

    return ret;
}

/*
   G = [(2 * \sum_i^n (i * x_i)) / (n * \sum_i^n x_i )] - [(n + 1)/n]

   where x is sorted: x_i <= x_{i+1}
*/

static double gini_coeff (const double *x, int t1, int t2, double **plz,
			  int *pn, int *err)
{
    int m = t2 - t1 + 1;
    double *sx = NULL;
    double csx = 0.0, sumx = 0.0, sisx = 0.0;
    double G;
    int t, n = 0;

    gretl_error_clear();

    sx = malloc(m * sizeof *sx);

    if (sx == NULL) {
	*err = E_ALLOC;
	return NADBL;
    }

    n = 0;
    for (t=t1; t<=t2; t++) {
	if (na(x[t])) {
	    continue;
	} else if (x[t] < 0.0) {
	    gretl_errmsg_set(_("illegal negative value"));
	    *err = E_DATA;
	    break;
	} else {
	    sx[n++] = x[t];
	    sumx += x[t];
	}
    }

    if (!*err && (n == 0 || sumx == 0.0)) {
	*err = E_DATA;
    }

    if (*err) {
	free(sx);
	return NADBL;
    }

    qsort(sx, n, sizeof *sx, gretl_compare_doubles);

#if 0
    if (1) {
	/* just testing alternative calculation for equivalence */
	double num, S[2];
	int s, fyi;

	num = S[0] = S[1] = 0.0;

	for (t=0; t<n; t++) {
	    fyi = 0;
	    s = t;
	    while (s < n && sx[s] == sx[t]) {
		fyi++;
		s++;
	    }
	    S[1] += (double) fyi/n * sx[t];
	    num += (double) fyi/n * (S[0] + S[1]);
	    t = s - 1;
	    S[0] = S[1];
	}

	G = 1.0 - num / S[1];
	fprintf(stderr, "alt G = %g\n", G);
    }
#endif

    for (t=0; t<n; t++) {
	csx += sx[t];
	sisx += (t + 1) * sx[t];
	sx[t] = csx / sumx; /* sx now = Lorenz curve */
    }

    G = 2.0 * sisx / (n * sumx) - ((double) n + 1) / n;

    if (plz != NULL) {
	*plz = sx;
	*pn = n;
    } else {
	free(sx);
    }

    return G;
}

/**
 * gretl_gini:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the Gini coefficient for the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_gini (int t1, int t2, const double *x)
{
    double g;
    int err = 0;

    g = gini_coeff(x, t1, t2, NULL, NULL, &err);

    if (err) {
	g = NADBL;
    }

    return g;
}

static int lorenz_graph (const char *vname, double *lz, int n)
{
    FILE *fp;
    double idx;
    int downsample = 0;
    int t, err = 0;

    fp = open_plot_input_file(PLOT_REGULAR, 0, &err);
    if (err) {
	return err;
    }

    print_keypos_string(GP_KEY_LEFT_TOP, fp);

    fprintf(fp, "set title '%s'\n", vname);
    fprintf(fp, "plot \\\n"
	    "'-' using 1:2 title '%s' w lines, \\\n"
	    "'-' using 1:2 notitle w lines\n",
	    _("Lorenz curve"));

    gretl_push_c_numeric_locale();

    if (n > 4000) {
	downsample = (int) (n / 1000.0);
    }

    for (t=0; t<n; t++) {
	if (downsample && t % downsample) {
	    continue;
	}
	idx = t + 1;
	fprintf(fp, "%g %g\n", idx / n, lz[t]);
    }
    fputs("e\n", fp);

    for (t=0; t<n; t++) {
	if (downsample && t % downsample) {
	    continue;
	}
	idx = ((double) t + 1) / n;
	fprintf(fp, "%g %g\n", idx, idx);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

/**
 * gini:
 * @varno: ID number of variable to examine.
 * @dset: dataset struct.
 * @opt: unused at present.
 * @prn: gretl printing struct.
 *
 * Graphs the Lorenz curve for variable @vnum and prints the
 * Gini coefficient.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gini (int varno, DATASET *dset, gretlopt opt, PRN *prn)
{
    double *lz = NULL;
    double gini;
    int n, fulln;
    int err = 0;

    gini = gini_coeff(dset->Z[varno], dset->t1, dset->t2,
		      &lz, &n, &err);

    if (err) {
	return err;
    }

    fulln = dset->t2 - dset->t1 - 1;
    pprintf(prn, "%s\n", dset->varname[varno], n);
    pprintf(prn, _("Number of observations = %d\n"), n);

    if (n < fulln) {
	pputs(prn, _("Warning: there were missing values\n"));
    }

    pputc(prn, '\n');

    pprintf(prn, "%s = %g\n", _("Sample Gini coefficient"), gini);
    pprintf(prn, "%s = %g\n", _("Estimate of population value"),
	    gini * (double) n / (n - 1));

    err = lorenz_graph(dset->varname[varno], lz, n);

    free(lz);

    return err;
}

/* Following: apparatus for Shapiro-Wilk normality test.
   Thanks to Marcin Blazejowski for the contribution.
*/

#ifndef min
# define min(a, b) ((a) > (b) ? (b) : (a))
#endif

static int sign (int a)
{
    return (a == 0)? 0 : ((a < 0)? -1 : 1);
}

static int compare_floats (const void *a, const void *b)
{
    const float *fa = (const float *) a;
    const float *fb = (const float *) b;

    return (*fa > *fb) - (*fa < *fb);
}

/* Algorithm AS 181.2 Appl. Statist. (1982) Vol. 31, No. 2
   Calculates the algebraic polynomial of order n-1 with
   array of coefficients cc.  Zero-order coefficient is
   cc(1) = cc[0]
*/

static double poly (const float *cc, int n, float x) {
    double p, ret = cc[0]; /* preserve precision! */
    int j;

    if (n > 1) {
	p = x * cc[n-1];
	for (j=n-2; j>0; j--) {
	    p = (p + cc[j]) * x;
	}
	ret += p;
    }

    return ret;
}

/* Calculate coefficients for the Shapiro-Wilk test */

static void sw_coeff (int n, int n2, float *a)
{
    const float c1[6] = { 0.f,.221157f,-.147981f,-2.07119f, 4.434685f, -2.706056f };
    const float c2[6] = { 0.f,.042981f,-.293762f,-1.752461f,5.682633f, -3.582633f };
    float summ2, ssumm2;
    float a0, a1, an = (float) n;
    float fac, an25, rsn;
    int i, i1, nn2 = n / 2;

    if (n == 3) {
	a[0] = sqrt(0.5);
    } else {
	an25 = an + .25;
	summ2 = 0.0;
	for (i=0; i<n2; i++) {
	    a[i] = (float) normal_cdf_inverse((i + 1 - .375f) / an25);
	    summ2 += a[i] * a[i];
	}
	summ2 *= 2.0;
	ssumm2 = sqrt(summ2);
	rsn = 1.0 / sqrt(an);
	a0 = poly(c1, 6, rsn) - a[0] / ssumm2;

	/* Normalize array a */
	if (n > 5) {
	    i1 = 2;
	    a1 = -a[1] / ssumm2 + poly(c2, 6, rsn);
	    fac = sqrt((summ2 - 2.0 * (a[0] * a[0]) - 2.0 * (a[1] * a[1])) /
		       (1.0 - 2.0 * (a0 * a0) - 2.0 * (a1 * a1)));
	    a[1] = a1;
	} else {
	    i1 = 1;
	    fac = sqrt((summ2 - 2.0 * (a[0] * a[0])) / (1.0 - 2.0 * (a0 * a0)));
	}
	a[0] = a0;
	for (i=i1; i<nn2; i++) {
	    a[i] /= -fac;
	}
    }
}

/* ALGORITHM AS R94 APPL. STATIST. (1995) vol. 44, no. 4, 547-551.
   Calculates the Shapiro-Wilk W test and its significance level,
   Rendered into C by Marcin Blazejowski, May 2008.
*/

static int sw_w (float *x, int n, int n1, float *a,
		 double *W, double *pval)
{
    const float z90 = 1.2816f;
    const float z95 = 1.6449f;
    const float z99 = 2.3263f;
    const float zm = 1.7509f;
    const float zss = .56268f;
    const float bf1 = .8378f;
    const double xx90 = .556;
    const double xx95 = .622;
    const float little = 1e-19f; /* "small" is reserved on win32 */
    const float pi6 = 1.909859f;
    const float stqr = 1.047198f;

    /* polynomial coefficients */
    const float  g[2] = { -2.273f, .459f };
    const float c3[4] = { .544f, -.39978f, .025054f, -6.714e-4f };
    const float c4[4] = { 1.3822f, -.77857f, .062767f, -.0020322f };
    const float c5[4] = { -1.5861f, -.31082f, -.083751f, .0038915f };
    const float c6[3] = { -.4803f, -.082676f, .0030302f };
    const float c7[2] = { .164f, .533f };
    const float c8[2] = { .1736f, .315f };
    const float c9[2] = { .256f, -.00635f };

    float r1, zbar, ssassx, gamma, range;
    float bf, ld, m, s, sa, xi, sx, xx, y, w1;
    float asa, ssa, z90f, sax, zfm, z95f, zsd, z99f, ssx, xsx;
    float an = (float) n;
    int i, j, ncens = n - n1;

    /* Check for zero range */
    range = x[n1 - 1] - x[0];
    if (range < little) {
	fprintf(stderr, "sw_w: range is too small\n");
	return 1;
    }

    /* Check for correct sort order on range-scaled X */
    /* In fact we don't need do this, since we use qsort() from libc
       and we can suppose, that qsort() can sort... */
    xx = x[0] / range;
    sx = xx;
    sa = -a[0];
    j = n - 1;
    for (i = 1; i < n1; --j) {
	xi = x[i] / range;
	sx += xi;
	++i;
	if (i != j) {
	    sa += sign(i - j) * a[min(i,j) - 1];
	}
	xx = xi;
    }

    /* Calculate W statistic as squared correlation between data and
       coefficients */
    sa /= n1;
    sx /= n1;
    ssa = ssx = sax = 0.0;
    j = n - 1;
    for (i=0; i<n1; i++, j--) {
	if (i != j) {
	    asa = sign(i - j) * a[min(i,j)] - sa;
	} else {
	    asa = -sa;
	}
	xsx = x[i] / range - sx;
	ssa += asa * asa;
	ssx += xsx * xsx;
	sax += asa * xsx;
    }

    /* W1 equals (1-W) calculated to avoid excessive rounding error
       for W very near 1 (a potential problem in very large samples)
    */
    ssassx = sqrt(ssa * ssx);
    w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
    *W = 1 - w1;

    /* Calculate significance level for W */
    if (n == 3) {
	/* exact P-value */
	*pval = pi6 * (asin(sqrt(*W)) - stqr);
	return 0;
    }

    y = log(w1);
    xx = log(an);

    if (n <= 11) {
	gamma = poly(g, 2, an);
	if (y >= gamma) {
	    /* FIXME: rather use an even smaller value, or NA? */
	    *pval = little;
	    return 0;
	}
	y = -log(gamma - y);
	m = poly(c3, 4, an);
	s = exp(poly(c4, 4, an));
    } else {
	/* n >= 12 */
	m = poly(c5, 4, xx);
	s = exp(poly(c6, 3, xx));
    }

    if (ncens > 0) {
	/* Censoring by proportion ncens/n.
	   Calculate mean and sd of normal equivalent deviate of W
	*/
	float delta = (float) ncens / an;

	ld = -log(delta);
	bf = 1.0 + xx * bf1;
	r1 = pow(xx90, (double) xx);
	z90f = z90 + bf * pow(poly(c7, 2, r1), (double) ld);
	r1 = pow(xx95, (double) xx);
	z95f = z95 + bf * pow(poly(c8, 2, r1), (double) ld);
	z99f = z99 + bf * pow(poly(c9, 2, xx), (double) ld);

	/* Regress Z90F,...,Z99F on normal deviates Z90,...,Z99 to get
	   pseudo-mean and pseudo-sd of z as the slope and intercept
	*/
	zfm = (z90f + z95f + z99f) / 3.0;
	zsd = (z90 * (z90f - zfm) + z95 * (z95f - zfm) + z99 * (z99f - zfm)) / zss;
	zbar = zfm - zsd * zm;
	m += zbar * s;
	s *= zsd;
    }

    *pval = normal_cdf_comp(((double) y - (double) m) / (double) s);

    return 0;
}

static int sw_sample_check (int n, int n1)
{
    int ncens = n - n1;
    float delta, an = n;

    if (n < 3 || n1 < 3) {
	fprintf(stderr, "There is no way to run SW test for less then 3 obs\n");
	return E_DATA;
    }

    if (ncens < 0 || (ncens > 0 && n < 20)) {
	fprintf(stderr, "sw_w: not enough uncensored obserations\n");
	return E_DATA;
    }

    delta = (float) ncens / an;
    if (delta > .8f) {
	fprintf(stderr, "sw_w: too many censored obserations\n");
	return E_DATA;
    }

    return 0;
}

/**
 * shapiro_wilk:
 * @x: data array.
 * @t1: starting observation.
 * @t2: ending observation.
 * @W: location to receive test statistic.
 * @pval: location to receive p-value.
 *
 * Computes the Shapiro-Wilk W statistic as a test for
 * normality of the data @x, and also the p-value for
 * the test. These are written into the pointer
 * arguments @W and @pval.
 *
 * Returns: 0 on success, non-zero on failure.
*/

int shapiro_wilk (const double *x, int t1, int t2, double *W, double *pval)
{
    int n1 = 0;         /* number of uncensored obs */
    float *a = NULL;	/* array of coefficients */
    float *xf = NULL;   /* copy of x in float format */
    int i, t, n2, n = 0;
    int err = 0;

    *W = *pval = NADBL;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) n++;
    }

    /* for now we assume all obs are uncensored */
    n1 = n;

    err = sw_sample_check(n, n1);
    if (err) {
	return err;
    }

    /* How many coeffs should be computed? */
    n2 = fmod(n, 2.0);
    n2 = (n2 == 0)?  n / 2 : (n - 1) / 2;

    xf = malloc(n * sizeof *xf);
    a = malloc(n2 * sizeof *a);

    if (xf == NULL || a == NULL) {
	err = E_ALLOC;
    } else {
	i = 0;
	for (t=t1; t<=t2; t++) {
	    if (!na(x[t])) {
		xf[i++] = x[t];
	    }
	}
	qsort(xf, n, sizeof *xf, compare_floats);
	/* Main job: compute W stat and its p-value */
	sw_coeff(n, n2, a);
	err = sw_w(xf, n, n1, a, W, pval);
    }

    free(a);
    free(xf);

    return err;
}

static double round2 (double x)
{
    x *= 100.0;
    x = (x - floor(x) < 0.5)? floor(x) : ceil(x);
    return x / 100;
}

/* Polynomial approximation to the p-value for the Lilliefors test,
   said to be good to 2 decimal places.  See Abdi and Molin,
   "Lilliefors/Van Soest's test of normality", in Salkind (ed.)
   Encyclopedia of Measurement and Statistics, Sage, 2007; also
   http://www.utd.edu/~herve/Abdi-Lillie2007-pretty.pdf
 */

static double lilliefors_pval (double L, int N)
{
    double b0 = 0.37872256037043;
    double b1 = 1.30748185078790;
    double b2 = 0.08861783849346;
    double b1pN = b1 + N;
    double Lm2 = 1.0 / (L * L);
    double A, pv;

    A = (-b1pN + sqrt(b1pN * b1pN - 4 * b2 * (b0 - Lm2))) / (2 * b2);

    pv = -0.37782822932809 + 1.67819837908004*A
	- 3.02959249450445*A*A + 2.80015798142101*pow(A, 3.)
	- 1.39874347510845*pow(A, 4.) + 0.40466213484419*pow(A, 5.)
	- 0.06353440854207*pow(A, 6.) + 0.00287462087623*pow(A, 7.)
	+ 0.00069650013110*pow(A, 8.) - 0.00011872227037*pow(A, 9.)
	+ 0.00000575586834*pow(A, 10.);

    if (pv < 0) {
	pv = 0;
    } else if (pv > 1) {
	pv = 1;
    } else {
	/* round to 2 digits */
	pv = round2(pv);
    }

    return pv;
}

static int lilliefors_test (const double *x, int t1, int t2,
			    double *L, double *pval)
{
    double *zx;   /* copy of x for sorting, z-scoring */
    int i, t, n = 0;
    int err = 0;

    *L = *pval = NADBL;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) n++;
    }

    if (n < 5) {
	/* we need more than 4 data points */
	return E_DATA;
    }

    zx = malloc(n * sizeof *zx);

    if (zx == NULL) {
	err = E_ALLOC;
    } else {
	double dx, mx = 0.0, sx = 0.0;
	double Phi, Dp = 0.0, Dm = 0.0;
	double Dmax = 0.0;

	/* compute sample mean and standard deviation plus
	   sorted z-scores */

	i = 0;
	for (t=t1; t<=t2; t++) {
	    if (!na(x[t])) {
		zx[i++] = x[t];
		mx += x[t];
	    }
	}
	mx /= n;

	for (t=t1; t<=t2; t++) {
	    if (!na(x[t])) {
		dx = x[t] - mx;
		sx += dx * dx;
	    }
	}
	sx = sqrt(sx / (n - 1));

	qsort(zx, n, sizeof *zx, gretl_compare_doubles);

	for (i=0; i<n; i++) {
	    zx[i] = (zx[i] - mx) / sx;
	}

	/* The following is based on the formula given in the
	   "nortest" package for GNU R, function lillie.test (written
	   by Juergen Gross).  It differs from the algorithm given by
	   Abdi and Molin (whose p-value approximation is used above).
	   The A & M account of the statistic itself seems to be
	   wrong: it consistently over-rejects (i.e. around 8 percent
	   of the time for a nominal 5 percent significance level).
	*/

	for (i=0; i<n; i++) {
	    Phi = normal_cdf(zx[i]);
	    Dp = (double) (i+1) / n - Phi;
	    Dm = Phi - (double) i / n;
	    if (Dp > Dmax) Dmax = Dp;
	    if (Dm > Dmax) Dmax = Dm;
	}

	*L = Dmax;
	*pval = lilliefors_pval(Dmax, n);
    }

    free(zx);

    return err;
}

static int skew_kurt_test (const double *x, int t1, int t2,
			   double *test, double *pval,
			   gretlopt opt)
{
    double skew, xkurt;
    int n, err = 0;

    *test = *pval = NADBL;

    err = series_get_moments(t1, t2, x, &skew, &xkurt, &n);

    if (!err) {
	if (opt & OPT_J) {
	    /* Jarque-Bera */
	    *test = (n / 6.0) * (skew * skew + xkurt * xkurt/4.0);
	} else {
	    /* Doornik-Hansen */
	    *test = doornik_chisq(skew, xkurt, n);
	}
    }

    if (!err && na(*test)) {
	err = E_NAN;
    }

    if (!na(*test)) {
	*pval = chisq_cdf_comp(2, *test);
    }

    return err;
}

static void print_normality_stat (double test, double pval,
				  gretlopt opt, PRN *prn)
{
    const char *tstrs[] = {
	N_("Shapiro-Wilk W"),
	N_("Jarque-Bera test"),
	N_("Lilliefors test"),
	N_("Doornik-Hansen test")
    };
    int i = (opt & OPT_W)? 0 : (opt & OPT_J)? 1 :
	(opt & OPT_L)? 2 : 3;

    if (na(test) || na(pval)) {
	return;
    }

    if (opt & OPT_L) {
	pprintf(prn, " %s = %g, %s ~= %g\n\n",
		_(tstrs[i]), test, _("with p-value"), pval);
    } else {
	pprintf(prn, " %s = %g, %s %g\n\n",
		_(tstrs[i]), test, _("with p-value"), pval);
    }
}

/**
 * gretl_normality_test:
 * @varno: ID number of the variable to process.
 * @dset: dataset struct.
 * @opt: OPT_A: all tests; OPT_D: Doornik-Hansen; OPT_W: Shapiro-Wilk;
 * OPT_J: Jarque-Bera; OPT_L: Lilliefors; default is Doornik-Hansen.
 * Also use OPT_Q for quiet.
 * @prn: gretl printing struct.
 *
 * Performs, and prints the results of, the specified test(s) randomness
 * for the variable specified by @v.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_normality_test (int varno, const DATASET *dset,
			  gretlopt opt, PRN *prn)
{
    gretlopt alltests = OPT_D | OPT_W | OPT_J | OPT_L;
    double test = NADBL;
    double pval = NADBL;
    double trec = NADBL;
    double pvrec = NADBL;
    int err;

    if (varno < 0 || varno >= dset->v) {
	return E_DATA;
    }

    err = incompatible_options(opt, OPT_A | alltests);
    if (err) {
	return err;
    }

    if (opt & OPT_A) {
	/* show all tests */
	opt |= alltests;
    }

    if (!(opt & alltests)) {
	/* no method selected: use Doornik-Hansen */
	opt |= OPT_D;
    }

    if (!(opt & OPT_Q)) {
	pprintf(prn, _("Test for normality of %s:"), dset->varname[varno]);
	if (opt & OPT_A) {
	    pputs(prn, "\n\n");
	} else {
	    pputc(prn, '\n');
	}
    }

    if (opt & OPT_D) {
	err = skew_kurt_test(dset->Z[varno], dset->t1, dset->t2,
			     &test, &pval, OPT_D);
	if (!err && !(opt & OPT_Q)) {
	    print_normality_stat(test, pval, OPT_D, prn);
	}
	if (!err) {
	    /* record Doornik-Hansen result by default */
	    trec = test;
	    pvrec = pval;
	}
    }

    if (opt & OPT_W) {
	err = shapiro_wilk(dset->Z[varno], dset->t1, dset->t2,
			   &test, &pval);
	if (!err && !(opt & OPT_Q)) {
	    print_normality_stat(test, pval, OPT_W, prn);
	}
    }

    if (opt & OPT_L) {
	err = lilliefors_test(dset->Z[varno], dset->t1, dset->t2,
			      &test, &pval);
	if (!err && !(opt & OPT_Q)) {
	    print_normality_stat(test, pval, OPT_L, prn);
	}
    }

    if (opt & OPT_J) {
	err = skew_kurt_test(dset->Z[varno], dset->t1, dset->t2,
			     &test, &pval, OPT_J);
	if (!err && !(opt & OPT_Q)) {
	    print_normality_stat(test, pval, OPT_J, prn);
	}
    }

    if (na(trec) && !na(test)) {
	trec = test;
    }

    if (na(pvrec) && !na(pval)) {
	pvrec = pval;
    }

    if (!na(trec) && !na(pvrec)) {
	record_test_result(trec, pvrec);
    }

    return err;
}

gretl_matrix *gretl_normtest_matrix (const double *y,
				     int t1, int t2,
				     gretlopt opt,
				     int *err)
{
    gretl_matrix *ret = NULL;
    double test = NADBL;
    double pval = NADBL;
    int do_all = 0;

    if (opt & OPT_J) {
	/* Jarque-Bera */
	*err = skew_kurt_test(y, t1, t2, &test, &pval, opt);
    } else if (opt & OPT_W) {
	/* Shapiro-Wilk */
	*err = shapiro_wilk(y, t1, t2, &test, &pval);
    } else if (opt & OPT_L) {
	/* Lilliefors */
	*err = lilliefors_test(y, t1, t2, &test, &pval);
    } else if (opt & OPT_A) {
	/* all tests */
	do_all = 1;
    } else {
	/* default: Doornik-Hansen */
	*err = skew_kurt_test(y, t1, t2, &test, &pval, opt);
    }

    if (do_all) {
	ret = gretl_matrix_alloc(4, 2);
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else {
	    char **Sc, **Sr;

	    Sc = strings_array_new(2);
	    Sr = strings_array_new(4);
	    if (Sc != NULL) {
		Sc[0] = gretl_strdup("test");
		Sc[1] = gretl_strdup("p-value");
		gretl_matrix_set_colnames(ret, Sc);
	    }
	    if (Sr != NULL) {
		Sr[0] = gretl_strdup("Doornik-Hansen");
		Sr[1] = gretl_strdup("Shapiro-Wilk");
		Sr[2] = gretl_strdup("Jarque-Bera");
		Sr[3] = gretl_strdup("Lilliefors");
		gretl_matrix_set_rownames(ret, Sr);
	    }

	    /* Doornik-Hansen */
	    skew_kurt_test(y, t1, t2, &test, &pval, OPT_NONE);
	    gretl_matrix_set(ret, 0, 0, test);
	    gretl_matrix_set(ret, 0, 1, pval);
	    /* Shapiro-Wilk */
	    shapiro_wilk(y, t1, t2, &test, &pval);
	    gretl_matrix_set(ret, 1, 0, test);
	    gretl_matrix_set(ret, 1, 1, pval);
	    /* Jarque-Bera */
	    skew_kurt_test(y, t1, t2, &test, &pval, OPT_J);
	    gretl_matrix_set(ret, 2, 0, test);
	    gretl_matrix_set(ret, 2, 1, pval);
	    /* Lilliefors */
	    lilliefors_test(y, t1, t2, &test, &pval);
	    gretl_matrix_set(ret, 3, 0, test);
	    gretl_matrix_set(ret, 3, 1, pval);
	}
	return ret;
    }

    if (!*err) {
	ret = gretl_matrix_alloc(1, 2);
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else {
	    ret->val[0] = test;
	    ret->val[1] = pval;
	}
    }

    return ret;
}
