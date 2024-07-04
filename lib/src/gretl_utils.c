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
#include "version.h"
#include "monte_carlo.h"
#include "gretl_func.h"
#include "objstack.h"
#include "cmd_private.h"
#include "libset.h"
#include "gretl_mt.h"
#include "uservar.h"
#include "gretl_panel.h"
#include "gretl_string_table.h"
#include "gretl_xml.h"
#include "forecast.h"
#include "gretl_typemap.h"

#ifdef USE_CURL
# include "gretl_www.h"
#endif

#ifdef WIN32
# include "gretl_win32.h"
#endif

#ifdef _OPENMP
# include <omp.h>
#endif

#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

#if defined(HAVE_MPI) || defined(USE_RLIB)
# include "gretl_foreign.h"
#endif

#if HAVE_GMP
# include <gmp.h>
#endif

#include <errno.h>

#if defined(_OPENMP) || defined(HAVE_MPI) || defined(WIN32)
# define WANT_XTIMER
#endif
#if !defined(_OPENMP) && !defined(WIN32)
# define NEED_ITIMER
#endif

static void gretl_tests_cleanup (void);

/**
 * date_as_double:
 * @t: observation number (zero-based).
 * @pd: data periodicity or frequency.
 * @sd0: floating point representation of starting date.
 *
 * Returns: the date corresponding to @t, as a double-precision number.
 */

double date_as_double (int t, int pd, double sd0)
{
    int ysd, yy, pp, yp;
    int p10 = 10;

    ysd = (int) sd0;

    if (pd == 1) {
        return (double) (ysd + t);
    }

    pp = pd;
    while ((pp = pp / 10)) {
        p10 *= 10;
    }

    pp = t % pd + p10 * (sd0 - ysd) + .5;
    if (pp != pd)  {
        yy = ysd + t/pd  + pp/pd + .5;
        yp = pp % pd;
    }  else {
        yy = ysd + t/pd + .5;
        yp = pp;
    }

    return yy + (double) yp / p10;
}

/**
 * ijton:
 * @i: row number (0-based)
 * @j: column number (0-based)
 * @nrows: number of rows (and columns) in symmetric matrix.
 *
 * Given a (row, column) reference into a symmetric matrix A, finds
 * the index into a 1-dimensional array x composed of the
 * non-redundant (lower) elements of A.
 *
 * E.g. for the 3 x 3 case with 6 non-redundant elements, 0 to 5,
 *
 *    A(0,0) = x[0]  A(0,1) = x[1]  A(0,2) = x[2]
 *    A(1,0) = x[1]  A(1,1) = x[3]  A(1,2) = x[4]
 *    A(2,0) = x[2]  A(2,1) = x[4]  A(2,1) = x[5]
 *
 * Returns: 0-based index into flat array.
 */

int ijton (int i, int j, int nrows)
{
    if (i > j) {
        int tmp = i;

        i = j;
        j = tmp;
    }

    return nrows * i + j - i - ((i - 1) * i / 2);
}

/**
 * transcribe_array:
 * @targ: arrat to which to write.
 * @src: array from which to read.
 * @dset: data information struct.
 *
 * Copy from @src to @targ, skipping any missing values,
 * over the sample range defined in @dset.
 *
 * Returns: the number of valid observations put into @targ.
 */

int transcribe_array (double *targ, const double *src,
                      const DATASET *dset)
{
    int t, n = 0;

    for (t=dset->t1; t<=dset->t2; t++) {
        if (!na(src[t])) {
            targ[n++] = src[t];
        }
    }

    return n;
}

/**
 * gretl_isdummy:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series to examine.
 *
 * Check whether series @x has only 0 or 1 values over the
 * given sample range (or possibly missing values).
 *
 * Returns: 0 if the variable is not a 0/1 dummy, otherwise the
 * number of 1s in the series.
 */

int gretl_isdummy (int t1, int t2, const double *x)
{
    int t, m = 0, goodobs = 0;

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        }
        if (x[t] != 0.0 && x[t] != 1.0) {
            return 0;
        }
        if (x[t] == 1.0) {
            m++;
        }
        goodobs++;
    }

    if (m < goodobs) {
        return m;
    }

    return 0;
}

/**
 * gretl_iszero:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x has only zero (or missing) values
 * over the given sample range.
 *
 * Returns: 0 if the series contains any valid non-zero values,
 * otherwise 1.
 */

int gretl_iszero (int t1, int t2, const double *x)
{
    int t;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && x[t] != 0) {
            return 0;
        }
    }

    return 1;
}

/**
 * gretl_isconst:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x is constant over the
 * given sample range (aside from any missing values).
 *
 * Returns: 1 if the variable is constant, otherwise 0.
 */

int gretl_isconst (int t1, int t2, const double *x)
{
    int t, ret = 1;

    while (na(x[t1]) && t1 <= t2) {
        t1++;
    }

    if (t1 >= t2) {
        return 0;
    }

    for (t=t1+1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        }
        if (floatneq(x[t], x[t1])) {
            ret = 0;
            break;
        }
    }

    return ret;
}

/**
 * gretl_isstoch:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x is stochastic, and contains a
 * contiguous set of valid values within the given sample span.
 * The simple and fallible heuristic for a stochastic series is
 * that the differences between successive values are all the
 * same.
 *
 * Returns: 1 if the variable appears to be stochastic on
 * the specified criterion, otherwise 0.
 */

int gretl_isstoch (int t1, int t2, const double *x)
{
    double dx0;
    int na_count = 0;
    int multidiff = 0;
    int s1 = -1, s2 = -1;
    int t;

    if (t1 >= t2) {
        return 0;
    }

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            s1 = t;
            break;
        }
    }
    for (t=t2; t>=s1; t--) {
        if (!na(x[t])) {
            s2 = t;
            break;
        }
    }
    if (s1 < 0 || s2 < 0 || s1 >= s2) {
        return 0;
    }

    dx0 = x[s1+1] - x[s1];

    for (t=s1; t<=s2; t++) {
        if (na(x[t])) {
            na_count++;
            break;
        }
        if (t > s1 && !multidiff) {
            if (x[t] - x[t-1] != dx0) {
                multidiff = 1;
            }
        }
    }

    return (na_count == 0 && multidiff);
}

/**
 * gretl_isunits:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x equals 1 over the
 * given sample range (aside from any missing values).
 *
 * Returns: 1 if so, otherwise 0.
 */

int gretl_isunits (int t1, int t2, const double *x)
{
    int t, ret = 1;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && x[t] != 1.0) {
            ret = 0;
            break;
        }
    }

    return ret;
}

/**
 * gretl_iscount:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x contains nothing but non-negative
 * integer values (some of which are > 1) over the
 * given sample range.
 *
 * Returns: 1 if so, otherwise 0.
 */

int gretl_iscount (int t1, int t2, const double *x)
{
    int t, xi;
    int g1 = 0;
    int ret = 1;

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        }
        if (x[t] < 0.0) {
            ret = 0;
            break;
        }
        xi = x[t];
        if (x[t] != (double) xi) {
            ret = 0;
            break;
        }
        if (x[t] > 1.0) {
            g1 = 1;
        }
    }

    if (g1 == 0) {
        ret = 0;
    }

    return ret;
}

/**
 * gretl_isint:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x contains nothing but integer
 * values (disregarding NAs) over the given sample range.
 *
 * Returns: 1 if so, otherwise 0.
 */

int gretl_isint (int t1, int t2, const double *x)
{
    int t;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t]) && x[t] != floor(x[t])) {
	    return 0;
	}
    }

    return 1;
}

#define FEW_VALS 32
#define FEWER_VALS 8

static int few_vals (int t1, int t2, const double *x, double *ratio)
{
    double test[FEW_VALS];
    int match;
    int i, t, n = 0, nv = 0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            match = 0;
            for (i=0; i<nv; i++) {
                if (x[t] == test[i]) {
                    /* we've already seen this value */
                    match = 1;
                    break;
                }
            }
            if (!match) {
                /* x[t] is a "new" value */
                if (nv == FEW_VALS) {
                    /* hit the limit */
                    nv++;
                    break;
                }
                test[nv++] = x[t];
            }
            n++;
        }
    }

    /* ratio of distinct values to observations */
    *ratio = nv / (double) n;

    return nv;
}

/**
 * gretl_isdiscrete:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Checks the variable @x over the range @t1 to @t2 for discreteness.
 * This is a heuristic whose components are (a) whether the values
 * are "fairly round" (multiples of 0.25) or not, and, if test (a) is
 * passed, (b) whether the variable takes on only "few" distinct
 * values.
 *
 * Returns: 0 if test (a) is not passed or the number of distinct values
 * is > 32; else 1 if the number of distinct values is <= 32; else 2 if
 * the number of distinct values is <= 8.  A return of 1 is supposed
 * to indicate that it's "reasonable" to treat @x as discrete, while
 * a return of 2 indicates that it's probably unreasonable _not_ to
 * treat @x as discrete for the purpose of drawing up a frequency
 * distribution.
 */

int gretl_isdiscrete (int t1, int t2, const double *x)
{
    int t, n = 0, disc = 1;
    int allints = 1;
    double r = 0;

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        }
        n++;
        if (!ok_int(x[t])) {
            allints = disc = 0;
            break;
        }
        r = x[t] - floor(x[t]);
        if (allints && r != 0) {
            allints = 0;
        }
        if (r != 0.0 && r != 0.25 && r != 0.5 && r != 0.75) {
            disc = 0;
            break;
        }
    }

    if (n == 0) {
        disc = 0;
    }

    if (disc) {
        n = few_vals(t1, t2, x, &r);
        if (allints) {
            if (n <= FEW_VALS && r < 0.2) {
                disc = 2;
            } else {
                disc = (n <= FEWER_VALS)? 2 : 1;
            }
        } else if (n > FEW_VALS) {
            disc = 0;
        } else if (r > 0.9 && n > 30) {
            /* somewhat arbitrary: but if r (= ratio of distinct
               values to number of cases) is "too high", and the
               number of observations is not tiny, perhaps we should
               not automatically take the var as discrete
            */
            disc = 0;
        } else if (n <= FEWER_VALS) {
            disc = 2;
        }
    }

    return disc;
}

int accept_as_discrete (const DATASET *dset, int v, int strict)
{
    if (series_is_discrete(dset, v)) {
        /* the series has been explicitly marked as discrete */
        return 1;
    } else {
        /* check for plausibility of discreteness */
        int d = gretl_isdiscrete(dset->t1, dset->t2, dset->Z[v]);

        return strict ? d > 1 : d > 0;
    }
}

/**
 * gretl_ispositive:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 * @strict: boolean, strict inequality.
 *
 * Returns: 1 if all non-missing values in @x, over the range @t1 to
 * @t2, are positive (if strict==1) or non-negative (if strict==0),
 * otherwise 0.
 */

int gretl_ispositive (int t1, int t2, const double *x, int strict)
{
    int t;

    if (strict) {
        for (t=t1; t<=t2; t++) {
            if (x[t] <= 0) {
                return 0;
            }
        }
    } else {
        for (t=t1; t<=t2; t++) {
            if (x[t] < 0) {
                return 0;
            }
        }
    }

    return 1;
}

/**
 * gretl_is_oprobit_ok:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Checks the variable @x over the range @t1 to @t2 for its
 * suitability as the dependent variable in an ordered probit
 * analysis.  The criterion used is that the variable has
 * only non-negative integer values.
 *
 * Returns: 1 if the test succeeds, otherwise 0.
 */

int gretl_is_oprobit_ok (int t1, int t2, const double *x)
{
    int t, n = 0, d = 1;

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        }
        n++;
        if (x[t] != floor(x[t]) || x[t] < 0) {
            d = 0;
            break;
        }
    }

    return (d > 0 && n > 0);
}

/**
 * true_const:
 * @v: index number of variable to test.
 * @dset: dataset struct.
 *
 * Check whether variable Z[v] equals 1 over the sample
 * range given in @dset (aside from any missing values).
 *
 * Returns: 1 if so, otherwise 0.
 */

int true_const (int v, const DATASET *dset)
{
    if (v < 0 || v >= dset->v) {
        return 0;
    }

    return gretl_isunits(dset->t1, dset->t2, dset->Z[v]);
}

/**
 * gretl_compare_doubles:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with qsort.  Sorts doubles in
 * ascending order.
 *
 * Returns: appropriate value for qsort.
 */

int gretl_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    if (isnan(*da) || isnan(*db)) {
        if (!isnan(*da)) {
            return -1;
        } else if (!isnan(*db)) {
            return 1;
        } else {
            return 0;
        }
    } else {
        return (*da > *db) - (*da < *db);
    }
}

/**
 * gretl_inverse_compare_doubles:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with qsort.  Sorts doubles in
 * descending order.
 *
 * Returns: appropriate value for qsort.
 */

int gretl_inverse_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    if (isnan(*da) || isnan(*db)) {
        if (!isnan(*da)) {
            return 1;
        } else if (!isnan(*db)) {
            return -1;
        } else {
            return 0;
        }
    } else {
        return (*da < *db) - (*da > *db);
    }
}

/**
 * gretl_compare_ints:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with qsort.  Sorts integers in
 * ascending order.
 *
 * Returns: appropriate value for qsort.
 */

int gretl_compare_ints (const void *a, const void *b)
{
    const int *ia = (const int *) a;
    const int *ib = (const int *) b;

    return *ia - *ib;
}

/**
 * gretl_compare_strings:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with qsort.  Sorts strings in
 * alphabetical order.
 *
 * Returns: appropriate value for qsort.
 */

int gretl_compare_strings (const void *a, const void *b)
{
    const char **sa = (const char **) a;
    const char **sb = (const char **) b;

    return g_strcmp0(*sa, *sb);
}

/**
 * count_distinct_values:
 * @x: sorted array of doubles.
 * @n: number of elements in array.
 *
 * Returns: the number of distinct values in array @x,
 * provided that @x is already sorted.
 */

int count_distinct_values (const double *x, int n)
{
    int i, c = 1;

    for (i=1; i<n; i++) {
        if (x[i] != x[i-1]) {
            c++;
        }
    }

    return c;
}

/**
 * count_distinct_int_values:
 * @x: sorted array of ints.
 * @n: number of elements in array.
 *
 * Returns: the number of distinct values in array @x,
 * provided that @x is already sorted.
 */

int count_distinct_int_values (const int *x, int n)
{
    int i, c = 1;

    for (i=1; i<n; i++) {
        if (x[i] != x[i-1]) c++;
    }

    return c;
}

/**
 * rearrange_id_array:
 * @x: sorted array of doubles.
 * @m: number of distinct values in array.
 * @n: number of elements in array.
 *
 * Rearranges the sorted array @x such that the first @m
 * elements contain the @m distinct values in sorted order.
 *
 * Returns: 0 on success, 1 on error (in case @m is greater
 * than @n).
 */

int rearrange_id_array (double *x, int m, int n)
{
    int i, k = 1;

    if (m >= n || m == 1) {
        return 1;
    }

    for (i=1; i<n && k<m; i++) {
        if (x[i] != x[i-1]) {
            x[k++] = x[i];
        }
    }

    return 0;
}

/**
 * printlist:
 * @list: array of integers.
 * @msg: message to print along with @list (or NULL).
 *
 * Prints to stderr the given @list of integers along with a message.
 */

void printlist (const int *list, const char *msg)
{
    int i;

    if (msg) {
        fprintf(stderr, "%s:\n", msg);
    } else {
        fprintf(stderr, "list: ");
    }

    if (list == NULL) {
        fputs("list is NULL", stderr);
    } else {
        fprintf(stderr, "%d : ", list[0]);
        for (i=1; i<=list[0]; i++) {
            if (list[i] == LISTSEP) {
                fputs("; ", stderr);
            } else {
                fprintf(stderr, "%d ", list[i]);
            }
        }
    }

    fputc('\n', stderr);
}

/**
 * gretl_calculate_criteria:
 * @ess: error sum of squares.
 * @n: number of observations.
 * @k: number of parameters estimated.
 * @ll: pointer to recieve loglikelihood.
 * @aic: pointer to recieve Akaike criterion.
 * @bic: pointer to recieve Schwarz Bayesian criterion.
 * @hqc: pointer to recieve Hannan-Quinn criterion.
 *
 * Calculates model selection criteria based on @ess, @nobs and
 * @k, for a model estimated via least squares.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_calculate_criteria (double ess, int n, int k,
                              double *ll, double *aic, double *bic,
                              double *hqc)
{
    double lnl, c[3];
    int err = 0;

    if (na(ess) || ess <= 0.0 || n <= k) {
        err = 1;
    } else {
        errno = 0;
        lnl = -.5 * n * log(ess);
        if (errno == EDOM || errno == ERANGE) {
            err = 1;
        } else {
            lnl += -.5 * n * (1 + LN_2_PI - log(n));
            c[0] = -2.0 * lnl + 2 * k;
            c[1] = -2.0 * lnl + k * log(n);
            c[2] = -2.0 * lnl + 2 * k * log(log(n));
        }
    }

    if (err) {
        *ll = *aic = *bic = *hqc = NADBL;
    } else {
        *ll = lnl;
        *aic = c[0];
        *bic = c[1];
        *hqc = c[2];
    }

    return err;
}

/**
 * ls_criteria:
 * @pmod: pointer to gretl model structure.
 *
 * Fills out the model selection criteria members of @pmod, using
 * gretl_calculate_criteria().
 *
 * Returns: 0 on success, non-zero on error.
 */

int ls_criteria (MODEL *pmod)
{
    double ll, aic, bic, hqc;
    int err;

    err = gretl_calculate_criteria(pmod->ess, pmod->nobs, pmod->ncoeff,
                                   &ll, &aic, &bic, &hqc);

    pmod->lnL = ll;
    pmod->criterion[C_AIC] = aic;
    pmod->criterion[C_BIC] = bic;
    pmod->criterion[C_HQC] = hqc;

    return err;
}

static char *
real_format_obs (char *obs, int maj, int min, int pd, char sep)
{
    if (pd >= 10) {
        int pdp = pd / 10, minlen = 2;
        char fmt[18];

        while ((pdp = pdp / 10)) minlen++;
        sprintf(fmt, "%%d%c%%0%dd", sep, minlen);
        sprintf(obs, fmt, maj, min);
    } else {
        sprintf(obs, "%d%c%d", maj, sep, min);
    }

    return obs;
}

/**
 * format_obs:
 * @obs: target string (should be of length #OBSLEN).
 * @maj: major period (e.g. year).
 * @min: minor period (e.g. quarter, month).
 * @pd: data frequency.
 *
 * Prints to @obs the gretl-type date string representing
 * the observation given by @maj, @min and @pd.
 *
 * Returns: @obs.
 */

char *format_obs (char *obs, int maj, int min, int pd)
{
    return real_format_obs(obs, maj, min, pd, ':');
}

static int get_stobs_maj_min (char *stobs, int structure,
                              int *maj, int *min)
{
    int dotc = 0;
    char *p = stobs;
    int err = 0;

    while (*p) {
        if (*p == ':') {
            *p = '.';
            dotc++;
        } else if (*p == '.') {
            dotc++;
        } else if (!isdigit((unsigned char) *p)) {
            err = 1;
            break;
        }
        p++;
    }

    if (!err) {
        if (dotc > 1 || *stobs == '.' ||
            stobs[strlen(stobs) - 1] == '.') {
            err = 1;
        }
    }

    if (!err) {
        if (dotc > 0) {
            sscanf(stobs, "%d.%d", maj, min);
            if (*maj <= 0 || *min <= 0) {
                err = 1;
            }
        } else {
            sscanf(stobs, "%d", maj);
            if (*maj <= 0 && structure != SPECIAL_TIME_SERIES) {
                err = 1;
            }
        }
    }

    return err;
}

static int
catch_setobs_errors (const char *stobs, int pd, int min, int structure)
{
    int panel = structure == STACKED_TIME_SERIES ||
        structure == STACKED_CROSS_SECTION;
    int err = 0;

    if (pd == 1) {
        if (min > 0) {
            gretl_errmsg_set(_("no ':' allowed in starting obs with "
                               "frequency 1"));
            err = 1;
        } else if (panel) {
            gretl_errmsg_set(_("panel data must have frequency > 1"));
            err = 1;
        }
    } else {
        if (min == 0) {
            gretl_errmsg_set(_("starting obs must contain a ':' with "
                               "frequency > 1"));
            err = 1;
        } else if (min > pd) {
            gretl_errmsg_sprintf(_("starting obs '%s' is incompatible with frequency"),
                                 stobs);
            err = 1;
        } else if (structure == CROSS_SECTION) {
            gretl_errmsg_set(_("cross-sectional data: frequency must be 1"));
            err = 1;
        }
    }

    return err;
}

static int invalid_stobs (const char *s)
{
    gretl_errmsg_sprintf(_("starting obs '%s' is invalid"), s);
    return E_DATA;
}

static void maybe_fix_daily_start (guint32 *ed, int pd)
{
    int wday = weekday_from_epoch_day(*ed);
    int fix = 0;

    if (wday == G_DATE_SUNDAY) {
        /* 5- or 6-day data: sunday not valid */
        fix = 1;
    } else if (pd == 5 && wday == G_DATE_SATURDAY) {
        /* 5-day data: saturday not valid */
        fix = 2;
    }

    if (fix) {
        char *fixed, *msg;

        *ed += fix;
        fixed = ymd_extended_from_epoch_day(*ed, 0, NULL);
        msg = gretl_strdup_printf(_("the starting date was corrected to Monday %s"),
                                  fixed);
        gretl_warnmsg_set(msg);
        free(msg);
        free(fixed);
    }
}

#define likely_calendar_obs_string(s) (strchr(s, '-') || strchr(s, '/'))

#define recognized_ts_frequency(f) (f == 4 || f == 12 || f == 24)

static int process_starting_obs (const char *stobs_in, int pd,
                                 int *pstructure, double *psd0,
                                 guint32 *ped0, gretlopt opt)
{
    char stobs[OBSLEN];
    int structure = *pstructure;
    double sd0 = 0.0;
    int maybe_tseries = 1;
    int dated = 0;
    int err = 0;

    if (structure == CROSS_SECTION ||
        structure == STACKED_TIME_SERIES ||
        structure == STACKED_CROSS_SECTION) {
        maybe_tseries = 0;
    }

    *stobs = '\0';
    strncat(stobs, stobs_in, OBSLEN - 1);

    /* check for possible calendar date */
    if (likely_calendar_obs_string(stobs)) {
        if (maybe_tseries) {
            dated = 1;
        } else {
            return invalid_stobs(stobs);
        }
    } else {
        ; /* stobs[8] = '\0'; */
    }

    /* 2023-01-04: if the likely_calendar_obs_string() test failed, we
       were truncating @stobs to 8 bytes. Now I don't understand
       why. That truncation prevents interpretation of, e.g., "setobs
       24 726468:01" as calling for hourly data starting in epoch day
       726468, since it knocks off the trailing '1'.
    */

    if (dated) {
        if (pd == 5 || pd == 6 || pd == 7 || pd == 52) {
            /* calendar-dated data, daily or weekly */
            guint32 ed0 = get_epoch_day(stobs);

            if (ed0 <= 0) {
                return invalid_stobs(stobs);
            } else {
                if (pd < 7) {
                    maybe_fix_daily_start(&ed0, pd);
                }
                sd0 = ed0;
                *ped0 = ed0;
                structure = TIME_SERIES;
            }
        } else {
            return invalid_stobs(stobs);
        }
    } else if (structure == TIME_SERIES && pd == 10) {
        /* decennial data */
        sd0 = (double) atoi(stobs);
    } else {
        int maj = 0, min = 0;

        if (get_stobs_maj_min(stobs, structure, &maj, &min)) {
            return invalid_stobs(stobs);
        }

        if ((pd == 5 || pd == 6 || pd == 7 || pd == 52) &&
            min == 0 && maybe_tseries) {
            /* catch undated daily or weekly data */
            structure = TIME_SERIES;
        } else {
            if (catch_setobs_errors(stobs, pd, min, structure)) {
                return E_DATA;
            } else if (pd == 1) {
                sprintf(stobs, "%d", maj);
                if (structure == STRUCTURE_UNKNOWN) {
                    if (maj > 1) {
                        structure = TIME_SERIES; /* annual? */
                    } else {
                        structure = CROSS_SECTION;
                    }
                }
            } else {
                if (structure == TIME_SERIES && min > 0 &&
                    !recognized_ts_frequency(pd)) {
                    if (opt & OPT_I) {
                        return invalid_stobs(stobs);
                    } else {
                        structure = SPECIAL_TIME_SERIES;
                    }
                }
                real_format_obs(stobs, maj, min, pd, '.');
                if (structure == STRUCTURE_UNKNOWN &&
                    recognized_ts_frequency(pd)) {
                    structure = TIME_SERIES;
                }
            }
        }

        /* for non-calendar data */
        sd0 = dot_atof(stobs);
    }

    if (!err) {
        *pstructure = structure;
        *psd0 = sd0;
    }

    return err;
}

/**
 * set_obs:
 * @parm1: first parameter.
 * @parm2: second parameter.
 * @dset: dataset struct.
 * @opt: %OPT_S for stacked time-series, %OPT_C for stacked cross-section,
 * %OPT_T for time series, %OPT_X for cross section, %OPT_P to set
 * panel structure via two variables representing unit and period
 * respectively. For data already set as panel, %OPT_G to set panel
 * group names or %OPT_I to set panel time-dimension information.
 *
 * Set the frequency and initial observation for a dataset, or
 * in the case of a panel dataset, extra group or time information.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int set_obs (const char *parm1, const char *parm2,
             DATASET *dset, gretlopt opt)
{
    const char *stobs = NULL;
    int structure = STRUCTURE_UNKNOWN;
    double sd0 = dset->sd0;
    guint32 ed0 = 0;
    int pd, panel = 0;
    int err = 0;

    if (dset == NULL) {
        return E_NODATA;
    }

    if ((opt & (OPT_R | OPT_P)) && dset->Z == NULL) {
        return E_NODATA;
    }

    if ((opt & (OPT_G | OPT_I)) && !dataset_is_panel(dset)) {
        return E_DATA;
    }

    /* prevent substantive reorganization of the dataset
       within a function */
    if ((opt & (OPT_P | OPT_C)) && gretl_function_depth() > 0) {
        gretl_errmsg_set(_("You cannot do this within a function"));
        return E_DATA;
    }

    gretl_error_clear();

    if (opt & OPT_R) {
        /* restructure panel: "hidden" option */
        return switch_panel_orientation(dset);
    }

    if (opt == OPT_NONE && parm1 == NULL && parm2 == NULL) {
        /* we'll just print current obs info */
        return 0;
    }

    if (parm1 == NULL || (!(opt & OPT_G) && parm2 == NULL)) {
        /* all cases need at least one param, most need two */
        return E_ARGS;
    }

    if (opt & OPT_P) {
        return set_panel_structure_from_varnames(parm1, parm2, dset);
    } else if (opt & OPT_G) {
        /* --panel-groups */
        return set_panel_group_strings(parm1, parm2, dset);
    }

    /* now we get down to business */

    pd = gretl_int_from_string(parm1, &err);
    if (!err && pd < 1) {
        gretl_errmsg_sprintf(_("frequency (%d) does not seem to make sense"), pd);
        return E_DATA;
    }

    stobs = parm2;

    /* if an explicit structure option was passed in, respect it */
    if (opt == OPT_X) {
        structure = CROSS_SECTION;
    } else if (opt == OPT_T) {
        structure = TIME_SERIES;
    } else if (opt == OPT_S) {
        structure = STACKED_TIME_SERIES;
        panel = 1;
    } else if (opt == OPT_C) {
        structure = STACKED_CROSS_SECTION;
        panel = 1;
    } else if (opt == OPT_N) {
        structure = SPECIAL_TIME_SERIES;
    } else if (opt == OPT_I) {
        /* --panel-time */
        structure = TIME_SERIES;
    }

    if (panel && dset->n > 0 && pd > dset->n) {
        gretl_errmsg_sprintf(_("frequency (%d) does not seem to make sense"), pd);
        return 1;
    }

    err = process_starting_obs(stobs, pd, &structure, &sd0, &ed0, opt);

    if (err) {
        return err;
    }

    if (opt == OPT_I) {
        /* --panel-time */
        dset->panel_pd = pd;
        dset->panel_sd0 = sd0;
        return 0;
    }

    if (panel && dset->n % pd != 0) {
        int sts = structure == STACKED_TIME_SERIES;

        gretl_errmsg_sprintf(_("Panel datasets must be balanced.\n"
                               "The number of observations (%d) is not a multiple\n"
                               "of the number of %s (%d)."),
                             dset->n, sts ? _("periods") : _("units"), pd);
        return E_DATA;
    }

    if (ed0 > 0) {
        /* replace any existing markers with date strings */
        dataset_destroy_obs_markers(dset);
    } else if (structure == TIME_SERIES && (pd == 1 || pd == 4 || pd == 12)) {
        /* force use of regular time-series obs labels */
        dataset_destroy_obs_markers(dset);
    }

    dset->pd = pd;
    dset->structure = structure;
    dset->sd0 = sd0;

    if (ed0 > 0) {
        calendar_date_string(dset->stobs, 0, dset);
        calendar_date_string(dset->endobs, dset->n - 1, dset);
    } else {
        ntolabel(dset->stobs, 0, dset);
        ntolabel(dset->endobs, dset->n - 1, dset);
    }

    /* pre-process stacked cross-sectional panels: put into canonical
       stacked time series form
    */
    if (dset->structure == STACKED_CROSS_SECTION) {
        err = switch_panel_orientation(dset);
    }

#if 0
    fprintf(stderr, "setobs: pd=%d, stobs=%s, sd0=%g, markers=%d, S=%p\n",
            dset->pd, dset->stobs, dset->sd0, dset->markers, (void *) dset->S);
#endif

    return err;
}

/**
 * simple_set_obs:
 * @dset: pointer to data information struct.
 * @pd: data frequency.
 * @stobs: string representation of starting observation.
 * @opt: options flags.
 *
 * See the "setobs" command.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int simple_set_obs (DATASET *dset, int pd, const char *stobs,
                    gretlopt opt)
{
    int structure = STRUCTURE_UNKNOWN;
    double sd0 = dset->sd0;
    guint32 ed0 = 0;
    int panel = 0;
    int err = 0;

    if (opt == OPT_X) {
        structure = CROSS_SECTION;
    } else if (opt == OPT_T) {
        structure = TIME_SERIES;
    } else if (opt == OPT_S) {
        structure = STACKED_TIME_SERIES;
        panel = 1;
    } else if (opt == OPT_C) {
        structure = STACKED_CROSS_SECTION;
        panel = 1;
    } else if (opt == OPT_N) {
        structure = SPECIAL_TIME_SERIES;
    } else if (opt == OPT_I) {
        /* --panel-time */
        structure = TIME_SERIES;
    }

    err = process_starting_obs(stobs, pd, &structure, &sd0, &ed0, OPT_NONE);

    if (err) {
        return err;
    }

    if (panel && dset->n % pd != 0) {
        int sts = structure == STACKED_TIME_SERIES;

        gretl_errmsg_sprintf(_("Panel datasets must be balanced.\n"
                               "The number of observations (%d) is not a multiple\n"
                               "of the number of %s (%d)."),
                             dset->n, sts ? _("periods") : _("units"), pd);
        return E_DATA;
    }

    if (ed0 > 0) {
        /* replace any existing markers with date strings */
        dataset_destroy_obs_markers(dset);
    } else if (structure == TIME_SERIES && (pd == 1 || pd == 4 || pd == 12)) {
        /* force use of regular time-series obs labels */
        dataset_destroy_obs_markers(dset);
    }

    dset->pd = pd;
    dset->structure = structure;
    dset->sd0 = sd0;

    if (ed0 > 0) {
        calendar_date_string(dset->stobs, 0, dset);
        calendar_date_string(dset->endobs, dset->n - 1, dset);
    } else {
        ntolabel(dset->stobs, 0, dset);
        ntolabel(dset->endobs, dset->n - 1, dset);
    }

    if (dset->structure == STACKED_CROSS_SECTION) {
        err = switch_panel_orientation(dset);
    }

    return err;
}

/**
 * gretl_double_from_string:
 * @s: string to examine.
 * @err: location to receive error code.
 *
 * If @s is a valid string representation of a double,
 * return its value, otherwise if @s is the name of a
 * scalar variable, return the value of that variable,
 * otherwise set the content of @err to a non-zero value.
 *
 * Returns: double value.
 */

double gretl_double_from_string (const char *s, int *err)
{
    char *test;
    double x;

    if (s == NULL || *s == '\0') {
        *err = E_DATA;
        return NADBL;
    }

    if (isalpha(*s)) {
        return get_scalar_value_by_name(s, err);
    }

    gretl_push_c_numeric_locale();
    errno = 0;
    x = strtod(s, &test);
    gretl_pop_c_numeric_locale();

    if (*test != '\0' || errno == ERANGE) {
        *err = E_DATA;
        errno = 0;
        return NADBL;
    }

    return x;
}

/**
 * gretl_int_from_string:
 * @s: string to examine.
 * @err: location to receive error code.
 *
 * If @s is a valid string representation of an integer,
 * return that integer, otherwise if @s is the name of a
 * scalar variable, return the value of that variable,
 * provided it can be converted to an integer, otherwise
 * set the content of @err to a non-zero value.
 *
 * Returns: integer value.
 */

int gretl_int_from_string (const char *s, int *err)
{
    char *test;
    int n = 0;

    if (s == NULL || *s == '\0') {
        *err = E_DATA;
        return 0;
    }

    if (isalpha(*s)) {
        double x = get_scalar_value_by_name(s, err);

        if (!*err) {
            n = gretl_int_from_double(x, err);
        }
        return n;
    }

    errno = 0;
    n = strtol(s, &test, 10);

    if (*test != '\0' || errno == ERANGE) {
        *err = E_DATA;
        errno = 0;
        return 0;
    }

    return n;
}

/**
 * positive_int_from_string:
 * @s: string to examine.
 *
 * If @s is a valid string representation of a positive integer
 * that can be represented as a 32-bit signed value return that
 * integer, otherwise return -1.
 *
 * Returns: integer value.
 */

int positive_int_from_string (const char *s)
{
    int ret = -1;

    if (s != NULL && *s != '\0') {
        char *test;

        errno = 0;

        ret = strtol(s, &test, 10);
        if (*test != '\0' || !strcmp(s, test) ||
	    errno == ERANGE || ret <= 0) {
            ret = -1;
        }
    }

    return ret;
}

#if 0 /* not yet */

/**
 * natural_number_from_string:
 * @s: string to examine.
 *
 * If @s is a valid string representation of a natural number that
 * can be represented as a 32-bit signed integer return that number,
 * otherwise return -1.
 *
 * Returns: integer value.
 */

int natural_number_from_string (const char *s)
{
    int ret = -1;

    if (s != NULL && *s != '\0') {
	long lval;
        char *test;

        errno = 0;

        lval = strtol(s, &test, 10);
        if (*test != '\0' || !strcmp(s, test) || errno == ERANGE ||
	    lval < 0 || lval > INT_MAX) {
            ret = -1;
        } else {
	    ret = (int) lval;
    }

    return ret;
}

#endif

static int letter_to_int (char c)
{
    const char *s = "abcdefghij";
    int i = 0;

    while (*s) {
        if (c == *s) {
            return i;
        }
        s++;
        i++;
    }

    return 0;
}

/**
 * gretl_version_number:
 * @version: gretl program version in string form.
 *
 * Returns: the integer gretl version number.
 */

int gretl_version_number (const char *version)
{
    int vnum = 0;

    if (!strcmp(version, "@VERSION@")) {
	version = GRETL_VERSION;
    }

    if (atoi(version) >= 2015) {
        /* as in "2015d" and subsequent releases */
        int Y;
        char c;

        sscanf(version, "%d%c", &Y, &c);
        vnum = 10 * Y + letter_to_int(c);
    } else {
        /* old style: "1.9.14" or whatever */
        int x, y, z;

        sscanf(version, "%d.%d.%d", &x, &y, &z);
        vnum = 10000 * x + 100 * y + z;
    }

    return vnum;
}

/* Mapping from new-style version numbers, based on
   year, to old style based on maj.min.pl. This
   covers gretl 1.9.4 to 1.10.2.
*/

static int old_style_gretl_version_number (int v)
{
    const int vtrans[18][2] = {
        {10904, 20110},
        {10905, 20111},
        {10906, 20112},
        {10907, 20113},
        {10908, 20120},
        {10909, 20121},
        {10910, 20122},
        {10911, 20123},
        {10912, 20130},
        {10913, 20131},
        {10914, 20132},
        {10990, 20140},
        {10991, 20141},
        {10992, 20142},
        {11000, 20150},
        {11001, 20151},
        {11002, 20152},
        {11003, 20153}
    };
    int i;

    for (i=0; i<17; i++) {
        if (v == vtrans[i][1] || v < vtrans[i+1][1]) {
            return vtrans[i][0];
        }
    }

    /* default to 1.9.4 */

    return vtrans[0][0];
}

/**
 * gretl_version_string:
 * @targ: string into which to write (9 bytes minimum).
 * @vnum: integer program version number.
 *
 * Returns: the string representation of @vnum, in @targ.
 */

char *gretl_version_string (char *targ, int vnum)
{
    if (vnum >= 20153) {
        /* "2105d" or higher */
        const char *s = "abcdefghij";
        int y, idx;
        char c;

        y = vnum / 10;
        idx = vnum - 10 * y;

        if (idx >= 0 && idx < 10) {
            c = s[idx];
        } else {
            c = 'a';
        }

        sprintf(targ, "%d%c", y, c);
    } else {
        int x, y, z;

        if (vnum > 20000) {
            /* translate back to old-style */
            vnum = old_style_gretl_version_number(vnum);
        }

        x = vnum / 10000;
        y = (vnum - x * 10000) / 100;
        z = vnum % 100;

        sprintf(targ, "%d.%d.%d", x, y, z);
    }

    return targ;
}

/**
 * varnum_from_string:
 * @str: string representation of an integer ID number.
 * @dset: dataset information.
 *
 * Returns: integer ID number, or -1 on failure.
 */

int varnum_from_string (const char *str, DATASET *dset)
{
    int v = positive_int_from_string(str);

    if (v <= 0 || v >= dset->v) {
        v = -1;
    }

    return v;
}

/**
 * gretl_int_from_double:
 * @x: double-precision floating point value
 * @err: location to receive error code.
 *
 * Returns: the value of @x converted to an integer, if
 * possible. Otherwise returns -1 with @err set to a
 * non-zero value. Note that it is considered an
 * error if @x is "too far" from the nearest integer;
 * it must be "almost integral", with tolerance 0.001.
 */

int gretl_int_from_double (double x, int *err)
{
    int k = -1;

    if (na(x) || fabs(x) > INT_MAX || fabs(x - nearbyint(x)) > 0.001) {
        *err = E_INVARG;
    } else {
        k = (int) lrint(x);
    }

    return k;
}

/**
 * gretl_unsigned_from_double:
 * @x: double-precision floating point value
 * @err: location to receive error code.
 *
 * Returns: the value of @x converted to an unsigned
 * 32-bit integer, if possible. Otherwise returns -1
 * with @err set to a non-zero value. Note that it is
 * considered an error if @x is "too far" from the
 * nearest integer; it must be "almost integral", with
 * tolerance 1.0e-6.
 */

guint32 gretl_unsigned_from_double (double x, int *err)
{
    guint32 u = 0;

    if (na(x) || x < 0 || fabs(x) > G_MAXUINT32) {
        *err = E_INVARG;
    } else {
        double f = floor(x);
        double c = ceil(x);

        if (x - f < 1e-6) {
            u = f;
        } else if (c - x < 1e-6) {
            u = c;
        } else {
            *err = E_INVARG;
        }
    }

    return u;
}

/**
 * gretl_uint53_from_double:
 * @x: double-precision floating point value
 * @err: location to receive error code.
 *
 * Returns: the value of @x converted to an unsigned 64-bit
 * integer, if possible. Otherwise returns -1 with @err set
 * to a non-zero value. It is considered an error if @x is
 * not "almost integral", with tolerance 1.0e-6.
 */

guint64 gretl_uint53_from_double (double x, int *err)
{
    guint64 k = 0;

    if (na(x) || x > 9007199254740992 /* 2^53 */) {
        *err = E_INVARG;
    } else {
        double f = floor(x);
        double c = ceil(x);

        if (x - f < 1e-6) {
            k = f;
        } else if (c - x < 1e-6) {
            k = c;
        } else {
            *err = E_INVARG;
        }
    }

    return k;
}

/**
 * gretl_type_from_name:
 * @s: the name to check.
 * @dset: dataset information.
 *
 * Returns: non-zero if @s is the name of an existing
 * series, matrix, scalar, list, string or bundle,
 * otherwise 0.
 */

GretlType gretl_type_from_name (const char *s, const DATASET *dset)
{
    if (dset != NULL && gretl_is_series(s, dset)) {
        return GRETL_TYPE_SERIES;
    } else {
        return user_var_get_type_by_name(s);
    }
}

/**
 * copyvec:
 * @src: array of doubles.
 * @n: number of elements to copy.
 *
 * Returns: an allocated copy of the first @n elements of
 * array @src, or %NULL on failure.
 */

double *copyvec (const double *src, int n)
{
    double *targ = NULL;
    size_t sz = n * sizeof *targ;

    if (n > 0 && src != NULL) {
        targ = malloc(sz);
    }

    if (targ != NULL) {
        memcpy(targ, src, sz);
    }

    return targ;
}

/**
 * doubles_array_free:
 * @X: 2-dimensional array of doubles.
 * @m: number of sub-arrays.
 *
 * Frees a 2-dimensional array of doubles, first freeing
 * each sub-array.
 */

void doubles_array_free (double **X, int m)
{
    if (X != NULL) {
        int i;

        for (i=0; i<m; i++) {
            free(X[i]);
        }
        free(X);
    }
}

/**
 * doubles_array_new:
 * @m: number of sub-arrays.
 * @n: length of each sub-array.
 *
 * Allocates a 2-dimensional array of doubles, that is,
 * @m arrays each containing @n elements.  If @n is
 * zero the sub-arrays are just set to %NULL.
 *
 * Returns: the allocated array, or %NULL on failure.
 */

double **doubles_array_new (int m, int n)
{
    double **X;
    int i;

    if (m == 0) {
        return NULL;
    }

    X = malloc(m * sizeof *X);
    if (X == NULL) {
        return X;
    }

    for (i=0; i<m; i++) {
        X[i] = NULL;
    }

    if (n > 0) {
        for (i=0; i<m; i++) {
            X[i] = malloc(n * sizeof **X);
            if (X[i] == NULL) {
                doubles_array_free(X, m);
                X = NULL;
                break;
            }
        }
    }

    return X;
}

/**
 * doubles_array_new0:
 * @m: number of sub-arrays.
 * @n: length of each sub-array.
 *
 * Works just as doubles_array_new(), except that on
 * successful allocation all values in the arrays are
 * set to zero.
 *
 * Returns: the allocated array, or %NULL on failure.
 */

double **doubles_array_new0 (int m, int n)
{
    double **X = doubles_array_new(m, n);

    if (X != NULL && n > 0) {
        int i, j;

        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                X[i][j] = 0.0;
            }
        }
    }

    return X;
}

/**
 * doubles_array_adjust_length:
 * @X: original two-dimensional array.
 * @m: number of sub-arrays in @X.
 * @new_n: new length of each sub-array.
 *
 * For each of the @m sub-arrays in @X, reallocate to
 * a length of @new_n.
 *
 * Returns: 0 on success, non-zero on error.
 */

int doubles_array_adjust_length (double **X, int m, int new_n)
{
    int err = 0;

    if (X == NULL || m == 0) {
        ; /* no-op */
    } else {
        double *xi;
        int i;

        for (i=0; i<m && !err; i++) {
            if (new_n == 0) {
                free(X[i]);
                X[i] = NULL;
            } else {
                xi = realloc(X[i], new_n * sizeof *xi);
                if (xi == NULL) {
                    err = E_ALLOC;
                } else {
                    X[i] = xi;
                }
            }
        }
    }

    return err;
}

/**
 * data_array_from_model:
 * @pmod: reference model.
 * @Z: main data array.
 * @missv: should equal 1 if there are missing values to be
 * skipped, else 0.
 *
 * Constructs a dataset containing all the variables referenced in
 * @pmod.  The arrays start at the correct sample offset for @pmod,
 * and are contiguous.  If @missvals equals 0, this is done by creating
 * a set of pointers into the main dataset, but if there are missing
 * values to be handled, the sub-arrays are newly allocated and purged
 * of NAs.
 *
 * Returns: two-dimensional array, or %NULL on failure.
 */

double **data_array_from_model (const MODEL *pmod, double **Z, int missv)
{
    double **X;
    int nv = pmod->list[0];
    int offset = pmod->t1;
    int v, i;

    if (missv) {
        X = doubles_array_new(nv, pmod->nobs);
    } else {
        X = malloc(nv * sizeof *X);
    }

    if (X == NULL) return NULL;

    if (missv) {
        int t, s;

        for (t=0; t<pmod->nobs; t++) {
            X[0][t] = 1.0;
        }

        for (i=1; i<nv; i++) {
            v = (i == 1)? pmod->list[1] : pmod->list[i + 1];
            s = 0;
            for (t=pmod->t1; t<=pmod->t2; t++) {
                if (!na(pmod->uhat[t])) {
                    X[i][s++] = Z[v][t];
                }
            }
        }
    } else {
        /* constant in slot 0 */
        X[0] = Z[0] + offset;

        /* dependent var in slot 1 */
        X[1] = Z[pmod->list[1]] + offset;

        /* independent vars in slots 2, 3, ... */
        for (i=2; i<nv; i++) {
            v = pmod->list[i + 1];
            X[i] = Z[v] + offset;
        }
    }

    return X;
}

#ifndef WIN32

int gretl_spawn (char *cmdline)
{
    GError *error = NULL;
    gchar *errout = NULL;
    gchar *sout = NULL;
    int ok, status;
    int ret = 0;

    gretl_error_clear();

    ok = g_spawn_command_line_sync(cmdline,
                                   &sout,   /* standard output */
                                   &errout, /* standard error */
                                   &status, /* exit status */
                                   &error);

    if (!ok) {
        gretl_errmsg_set(error->message);
        fprintf(stderr, "gretl_spawn: '%s'\n", error->message);
        g_error_free(error);
        if (errout != NULL) {
            fprintf(stderr, " stderr = '%s'\n", errout);
        }
        ret = 1;
    } else if (status != 0) {
        if (errout != NULL && *errout) {
            gretl_errmsg_set(errout);
            fprintf(stderr, "gretl_spawn: status = %d: '%s'\n", status, errout);
        } else if (sout != NULL && *sout) {
            gretl_errmsg_set(sout);
            fprintf(stderr, "gretl_spawn: status = %d: '%s'\n", status, sout);
        } else {
            gretl_errmsg_set(_("Command failed"));
            fprintf(stderr, "gretl_spawn: status = %d\n", status);
        }
        ret = 1;
    }

    if (errout != NULL) g_free(errout);
    if (sout != NULL) g_free(sout);

    if (ret) {
        fprintf(stderr, "Failed command: '%s'\n", cmdline);
    }

    return ret;
}

int gnuplot_make_image (const char *input_fname)
{
    gchar *cmdline;
    int err;

    cmdline = g_strdup_printf("\"%s\" \"%s\"",
                              gretl_gnuplot_path(),
                              input_fname);
    err = gretl_spawn(cmdline);
    g_free(cmdline);

    return err;
}

int gretl_pipe_output (gchar **argv, gchar **envp,
                       const char *currdir, PRN *prn,
                       gchar **errp)
{
    GError *error = NULL;
    gchar *errout = NULL;
    gchar *sout = NULL;
    int ok, status;
    int err = 0;

    gretl_error_clear();

    ok = g_spawn_sync(currdir,
                      argv,
                      envp,    /* may be NULL to inherit env */
                      G_SPAWN_SEARCH_PATH, /* ? */
                      NULL,    /* child_setup */
                      NULL,    /* user_data */
                      &sout,   /* standard output */
                      &errout, /* standard error */
                      &status, /* exit status */
                      &error);

    if (!ok) {
        gretl_errmsg_set(error->message);
        fprintf(stderr, "gretl_pipe_output: '%s'\n", error->message);
        g_error_free(error);
        if (errout != NULL) {
            fprintf(stderr, " stderr = '%s'\n", errout);
        }
        err = 1;
    } else if (status != 0) {
        if (errout != NULL && *errout) {
            gretl_errmsg_set(errout);
            fprintf(stderr, "gretl_pipe_output: status = %d: '%s'\n", status, errout);
        } else {
            gretl_errmsg_set(_("Command failed"));
            fprintf(stderr, "gretl_pipe_output: status = %d\n", status);
        }
        err = 1;
    }

    if (sout != NULL) {
        if (*sout != '\0') {
            pputs(prn, sout);
        }
        g_free(sout);
    }
    if (errout != NULL) {
        if (*errout != '\0') {
            if (errp != NULL) {
                *errp = errout;
                errout = NULL;
            } else {
                fputs(errout, stderr);
            }
        }
        g_free(errout);
    }

    return err;
}

#endif /* !WIN32 */

/* file copying */

int gretl_copy_file (const char *src, const char *dest)
{
    FILE *srcfd, *destfd;
    char buf[8192];
    size_t n;

    if (!strcmp(src, dest)) {
        /* @src and @dest are the same: no-op */
        return 0;
    }

    if ((srcfd = gretl_fopen(src, "rb")) == NULL) {
        gretl_errmsg_sprintf(_("Couldn't open %s"), src);
        return E_FOPEN;
    }

    if ((destfd = gretl_fopen(dest, "wb")) == NULL) {
        gretl_errmsg_sprintf(_("Couldn't write to %s"), dest);
        fclose(srcfd);
        return E_FOPEN;
    }

    while ((n = fread(buf, 1, sizeof buf, srcfd)) > 0) {
        fwrite(buf, 1, n, destfd);
    }

    fclose(srcfd);
    fclose(destfd);

    return 0;
}

static int maybe_unload_function_package (const char *s,
                                          PRN *prn)
{
    char *p, pkgname[FN_NAMELEN+4];
    const char *path;
    fnpkg *pkg;
    int done = 0;

    *pkgname = '\0';
    strncat(pkgname, s, FN_NAMELEN+3);
    p = strrchr(pkgname, '.');
    *p = '\0';

    pkg = get_function_package_by_name(pkgname);

    if (pkg != NULL) {
        path = get_function_package_path_by_name(pkgname);
        if (path != NULL) {
            function_package_unload_full_by_filename(path);
            done = 1;
        }
    }

    if (done) {
        pprintf(prn, _("Unloaded package %s\n"), pkgname);
    } else {
        pprintf(prn, _("Package %s was not loaded\n"), pkgname);
    }

    return 0;
}

/* note: this is not used for series */

int gretl_delete_var_by_name (const char *s, PRN *prn)
{
    int err = 0;

    if (s == NULL || *s == '\0') {
        return E_PARSE;
    }

    if (object_is_function_arg(s)) {
        gretl_errmsg_sprintf(_("The variable %s is read-only"), s);
        return E_DATA;
    }

    if (has_suffix(s, ".gfn")) {
        err = maybe_unload_function_package(s, prn);
    } else if (gretl_is_user_var(s)) {
        err = user_var_delete_by_name(s, prn);
    } else {
        /* try for a bundle member? */
        gchar *genstr = g_strdup_printf("%s=null", s);

        err = generate(genstr, NULL, GRETL_TYPE_ANY, OPT_P, prn);
        g_free(genstr);
    }

    return err;
}

/* apparatus to support gretl's stopwatch */

#ifdef WANT_XTIMER

/* Both OpenMP and MPI offer timers that return the
   time as a double (seconds and fractions thereof).
   If we may be using either of these we'll define
   the "xtimers".
*/

struct xtimer {
    int level;
    double t0;
};

static GPtrArray *xtimers;

static double get_xtime (void)
{
#if defined(HAVE_MPI)
    if (gretl_mpi_initialized()) {
        return gretl_mpi_time();
    }
#endif
#if defined(WIN32)
    return win32_get_time();
#elif defined(_OPENMP)
    return omp_get_wtime();
#endif
    return 0; /* not reached */
}

static struct xtimer *new_xtimer (int level)
{
    struct xtimer *xt = g_malloc(sizeof *xt);

    xt->level = level;
    xt->t0 = get_xtime();
    g_ptr_array_insert(xtimers, xtimers->len, xt);
    return xt;
}

static void xtimer_init (void)
{
    if (xtimers == NULL) {
        xtimers = g_ptr_array_sized_new(1);
        new_xtimer(gretl_function_depth());
    }
}

static struct xtimer *get_xtimer (void)
{
    struct xtimer *xt;
    int i, lev = gretl_function_depth();

    if (xtimers == NULL) {
        xtimers = g_ptr_array_sized_new(1);
        return new_xtimer(lev);
    }

    for (i=0; i<xtimers->len; i++) {
        xt = g_ptr_array_index(xtimers, i);
        if (xt->level == lev) {
            return xt;
        }
    }

    return new_xtimer(lev);
}

static double xtimer_stopwatch (void)
{
    struct xtimer *xt = get_xtimer();
    double t1 = get_xtime();
    double ret = t1 - xt->t0;

    xt->t0 = t1;

    return ret;
}

#endif /* WANT_XTIMER */

#ifdef NEED_ITIMER

/* Alternative timer, giving time in microseconds as a 64-bit int.
   If we don't have OpenMP we'll need this when not under MPI.
*/

struct itimer {
    int level;
    gint64 t0;
};

static GPtrArray *itimers;

static struct itimer *new_itimer (int level)
{
    struct itimer *it = g_malloc(sizeof *it);

    it->level = level;
    it->t0 = g_get_monotonic_time();
    g_ptr_array_insert(itimers, itimers->len, it);
    return it;
}

static void itimer_init (void)
{
    if (itimers == NULL) {
        itimers = g_ptr_array_sized_new(1);
        new_itimer(gretl_function_depth());
    }
}

static struct itimer *get_itimer (void)
{
    struct itimer *it;
    int i, lev = gretl_function_depth();

    if (itimers == NULL) {
        itimers = g_ptr_array_sized_new(1);
        return new_itimer(lev);
    }

    for (i=0; i<itimers->len; i++) {
        it = g_ptr_array_index(itimers, i);
        if (it->level == lev) {
            return it;
        }
    }

    return new_itimer(lev);
}

static double itimer_stopwatch (void)
{
    struct itimer *it = get_itimer();
    gint64 t1 = g_get_monotonic_time();
    double ret = (t1 - it->t0) * 1.0e-6;

    it->t0 = t1;

    return ret;
}

#endif /* NEED_ITIMER */

static void gretl_stopwatch_cleanup (void)
{
    int i;

#ifdef WANT_XTIMER
    if (xtimers != NULL) {
        for (i=0; i<xtimers->len; i++) {
            g_free(xtimers->pdata[i]);
        }
        g_ptr_array_free(xtimers, TRUE);
        xtimers = NULL;
    }
#endif

#ifdef NEED_ITIMER
    if (itimers != NULL) {
        for (i=0; i<itimers->len; i++) {
            g_free(itimers->pdata[i]);
        }
        g_ptr_array_free(itimers, TRUE);
        itimers = NULL;
    }
#endif
}

double gretl_stopwatch (void)
{
#if defined(HAVE_MPI)
    if (gretl_mpi_initialized()) {
        return xtimer_stopwatch();
    }
#endif
#if defined(_OPENMP)
    return xtimer_stopwatch();
#else
    return itimer_stopwatch();
#endif
}

static void gretl_stopwatch_init (void)
{
#if defined(HAVE_MPI)
    if (gretl_mpi_initialized()) {
        xtimer_init();
        return;
    }
#endif
#if defined(_OPENMP)
    xtimer_init();
#else
    itimer_init();
#endif
}

/* BLAS detection and analysis */

static char blas_core[64];
static char blas_parallel[32];
static char blas_version[32];

static int blas_variant;

static int parse_ldd_output (const char *s)
{
    char found[6] = {0};
    char line[512];
    int i = 0;
    int ret = BLAS_UNKNOWN;

    *line = '\0';

    while (*s) {
        if (*s == '\n') {
            /* got to the end of a line */
            line[i] = '\0';
	    if (strstr(line, "libopenblas")) {
		found[5] = 5;
	    } else if (strstr(line, "libblis")) {
		found[4] = 1;
	    } else if (strstr(line, "Accelerate.frame")) {
		found[3] = 1;
	    } else if (strstr(line, "libmkl")) {
		found[2] = 1;
	    } else if (strstr(line, "atlas")) {
		found[1] = 1;
	    } else if (strstr(line, "libblas")) {
		found[0] = 1;
	    }
	    *line = '\0';
	    i = 0;
	} else {
	    line[i++] = *s;
	}
	s++;
    }

    if (found[5]) {
	ret = BLAS_OPENBLAS;
    } else if (found[4]) {
        ret = BLAS_BLIS;
    } else if (found[3]) {
        ret = BLAS_VECLIB;
    } else if (found[2]) {
        ret = BLAS_MKL;
    } else if (found[1]) {
        ret = BLAS_ATLAS;
    } else if (found[0]) {
        ret = BLAS_NETLIB;
    } else {
        fputs("detect blas: found no relevant libs!\n", stderr);
    }

    return ret;
}

#if !defined(WIN32)

static int detect_blas_via_ldd (void)
{
    gchar *targ;
    gchar *sout = NULL;
    gchar *errout = NULL;
    gint status = 0;
    GError *gerr = NULL;
    int variant = 0;

#ifdef OS_OSX
    gchar *argv[4];
    targ = g_strdup_printf("%sgretlcli", gretl_bindir());
    argv[0] = "otool";
    argv[1] = "-L";
    argv[2] = targ;
    argv[3] = NULL;
#else
    gchar *argv[3];
    targ = g_strdup(GRETL_PREFIX "/lib/libgretl-1.0.so");
    argv[0] = "ldd";
    argv[1] = targ;
    argv[2] = NULL;
#endif

    g_spawn_sync(NULL, argv, NULL, G_SPAWN_SEARCH_PATH,
                 NULL, NULL, &sout, &errout,
                 &status, &gerr);

    if (gerr != NULL) {
        fprintf(stderr, "%s\n", gerr->message);
        g_error_free(gerr);
    } else if (status != 0) {
        fprintf(stderr, "%s exited with status %d\n", argv[0], status);
    } else if (sout != NULL) {
        variant = parse_ldd_output(sout);
    } else {
        fprintf(stderr, "%s: %s\n", argv[0], "Got no output");
    }

    g_free(sout);
    g_free(errout);
    g_free(targ);

    return variant;
}

#endif /* neither Windows nor Mac */

#include <dlfcn.h>

static void (*OB_set_num_threads) (int);
static int (*OB_get_num_threads) (void);
static void (*BLIS_set_num_threads) (int);
static int (*BLIS_get_num_threads) (void);
static void (*BLIS_init) (void);
static void (*BLIS_finalize) (void);
static void (*MKL_finalize) (void);
static void (*MKL_domain_set_num_threads) (int, int);
static int (*MKL_domain_get_max_threads) (int);
static void (*FLAME_init) (void);
static int (*FLAME_initialized) (void);
static void (*FLAME_finalize) (void);

static void register_openblas_details (void *handle)
{
    char *(*OB_get_corename) (void);
    int (*OB_get_parallel) (void);
    char *(*OB_get_config) (void);

    OB_get_corename = dlsym(handle, "openblas_get_corename");
    OB_get_parallel = dlsym(handle, "openblas_get_parallel");
    OB_get_config = dlsym(handle, "openblas_get_config");

    if (OB_get_corename != NULL) {
        char *s = OB_get_corename();

        if (s != NULL) {
            *blas_core = '\0';
            strncat(blas_core, s, 31);
        }
    } else {
        fprintf(stderr, "Couldn't find openblas_get_corename()\n");
    }

    if (OB_get_parallel != NULL) {
        int p = OB_get_parallel();

        if (p == 0) {
            strcpy(blas_parallel, "none");
        } else if (p == 1) {
            strcpy(blas_parallel, "pthreads");
        } else if (p == 2) {
            strcpy(blas_parallel, "OpenMP");
        }
    } else {
        fprintf(stderr, "Couldn't find openblas_get_parallel()\n");
    }

    if (OB_get_config != NULL) {
        char *ver = OB_get_config();

        if (ver != NULL) {
            *blas_version = '\0';
            sscanf(ver, "OpenBLAS %16[^ ]", blas_version);
        }
    } else {
        fprintf(stderr, "Couldn't find openblas_get_config()\n");
    }
}

static void register_blis_details (void *handle)
{
    /* The last element on arch_t enum in libblis =>
       => must be updated whenever new architecture/cpu model appears*/
    /* Shouldn't this come from a header? (Allin) */
    /* -> Well, this enum is defined in blis.h which we do not include.
       And we can't retrive it via dlopen(), can we? (Marcin) */
    const int BLIS_NUM_ARCHS = 26;
    char *buf = NULL;
    int id;

    /* Functions from libblis we need. */
    char *(*BLIS_get_version_str) (void);
    int (*BLIS_get_enable_threading) (void);
    int (*BLIS_get_enable_openmp) (void);
    int (*BLIS_get_enable_pthreads) (void);
    char *(*BLIS_arch_string) (int);
    int (*BLIS_arch_query_id) (void);

    BLIS_get_version_str = dlsym(handle, "bli_info_get_version_str");
    BLIS_get_enable_threading = dlsym(handle, "bli_info_get_enable_threading");
    BLIS_get_enable_openmp = dlsym(handle, "bli_info_get_enable_openmp");
    BLIS_get_enable_pthreads = dlsym(handle, "bli_info_get_enable_pthreads");
    BLIS_arch_string = dlsym(handle, "bli_arch_string");
    BLIS_arch_query_id = dlsym(handle, "bli_arch_query_id");

    if (BLIS_get_version_str != NULL) {
        buf = BLIS_get_version_str();
        if (buf != NULL) {
            *blas_version = '\0';
            strncat(blas_version, buf, 31);
            buf = NULL;
        }
    } else {
        fprintf(stderr, "Couldn't find bli_info_get_version_str()\n");
    }

    /* Model we have: threaded or sequential */
    if (!BLIS_get_enable_threading()) {
        buf = "sequential (non-threading)";
    } else {
        if (BLIS_get_enable_openmp()) {
            buf = "OpenMP";
        } else if (BLIS_get_enable_pthreads()) {
            buf = "pthreads";
        } else {
            buf = "unrecognized";
        }
    }
    if (buf != NULL) {
        *blas_parallel = '\0';
        strncat(blas_parallel, buf, 31);
        buf = NULL;
    }

    /* BLIS core in use */
    if (BLIS_arch_query_id != NULL) {
        id = BLIS_arch_query_id();
    } else {
        fprintf(stderr, "Couldn't find bli_arch_query_id()\n");
        id = BLIS_NUM_ARCHS-1;
    }
    if (BLIS_arch_string != NULL) {
        buf = BLIS_arch_string(id);
    } else {
        fprintf(stderr, "Couldn't find bli_arch_string()\n");
        buf = "unrecognized";
    }
    if (buf != NULL) {
        *blas_core = '\0';
        strncat(blas_core, buf, 31);
    }
}

static void register_mkl_details (void *handle)
{
    typedef struct {
        int MajorVersion;
        int MinorVersion;
        int UpdateVersion;
        char *ProductStatus;
        char *Build;
        char *Processor;
        char *Platform;
    } MKL_version;

    MKL_version *version;
    char *buf = NULL;
    int nt = 0, blas_nt = 0, id = 0;

    /* Functions from libmkl_intel_lp64.so we need. */
    void (*MKL_get_version) (MKL_version *);
    int (*MKL_get_max_threads) (void);
    int (*MKL_cbwr_get_auto_branch) (void);

    MKL_get_version = dlsym(handle, "MKL_Get_Version");
    MKL_get_max_threads = dlsym(handle, "MKL_Get_Max_Threads");
    MKL_cbwr_get_auto_branch = dlsym(handle, "MKL_CBWR_Get_Auto_Branch");

    /* MKL version */
    if (MKL_get_version != NULL) {
        version = malloc(sizeof *version);
        if (version == NULL) {
            return;
        } else {
            MKL_get_version(version);
            if (version != NULL) {
                *blas_version = '\0';
                snprintf(blas_version, 31, "%d.%d build %s",
                         version->MajorVersion,
                         version->UpdateVersion,
                         version->Build);
            }
        }
        free(version);
    } else {
        fprintf(stderr, "Couldn't find MKL_Get_Version()\n");
    }

    /* Model we have: threaded or sequential */
    if (MKL_get_max_threads != NULL && MKL_domain_get_max_threads != NULL) {
        nt = MKL_get_max_threads();
        blas_nt = MKL_domain_get_max_threads(1); /* MKL_DOMAIN_BLAS */
        if (nt == 1)  {
            buf = "sequential (non-threading)";
        } else {
            if (blas_nt == 1) {
                buf = "TBB";
            } else if (blas_nt > 1) {
                buf = "OpenMP";
            } else {
                buf = "unrecognized";
            }
        }
        strncat(blas_parallel, buf, 31);
        buf = NULL;
    } else {
        fprintf(stderr, "Couldn't find MKL_get_max_threads() and/or MKL_domain_get_max_threads()\n");
    }

    /* MKL core in use */
    if (MKL_cbwr_get_auto_branch != NULL) {
        id = MKL_cbwr_get_auto_branch();
        switch (id) {
        case 3: /* MKL_CBWR_COMPATIBLE */
            buf = "SSE2 without rcpps/rsqrtps instructions";
            break;
        case 4: /* MKL_CBWR_SSE2 */
            buf = "SSE2";
            break;
        case 5: /* MKL_CBWR_SSE3 */
            buf = "SSE2 (deprecated SSE3)";
            break;
        case 6: /* MKL_CBWR_SSSE3 */
            buf = "SSSE3";
            break;
        case 7: /* MKL_CBWR_SSE4_1 */
            buf = "SSE4-1";
            break;
        case 8: /* MKL_CBWR_SSE4_2 */
            buf = "SSE4-2";
            break;
        case 9: /* MKL_CBWR_AVX */
            buf = "AVX";
            break;
        case 10: /* MKL_CBWR_AVX2 */
            buf = "AVX2";
            break;
        case 11: /* MKL_CBWR_AVX512_MIC */
            buf = "AVX2 (deprecated AVX-512MIC)";
            break;
        case 12: /* MKL_CBWR_AVX512 */
            buf = "AVX-512";
            break;
        case 13: /* MKL_CBWR_AVX512_MIC_E1 */
            buf = "AVX2 (deprecated AVX-512MICE1)";
            break;
        case 14: /* MKL_CBWR_AVX512_E1 */
            buf = "AVX-512 with support of Vector Neural Network Instructions";
            break;
        default:
            buf = "unknown";
        }
        strncat(blas_core, buf, 63);
    }
}

static void libflame_utils (void *handle)
{
    /* Functions form libflame we need */
    FLAME_initialized = dlsym(handle, "FLA_Initialized");
    FLAME_finalize = dlsym(handle, "FLA_Finalize");
}

/* below: called in creating $sysinfo bundle */

int get_blas_details (char **s1, char **s2, char **s3)
{
    if (*blas_core == '\0' || *blas_parallel == '\0') {
        return 0;
    } else {
        *s1 = blas_core;
        *s2 = blas_parallel;
        if (s3 != NULL && *blas_version != '\0') {
            *s3 = blas_version;
        }
        return 1;
    }
}

/* for $sysinfo bundle */

#ifndef WIN32

#include <sys/resource.h>

int get_stack_size (void)
{
    struct rlimit rl;

    if (getrlimit(RLIMIT_STACK, &rl) == 0) {
        return (int) rl.rlim_cur;
    }

    return 0;
}

#endif

const char *blas_variant_string (void)
{
    if (blas_variant == BLAS_NETLIB) {
        return "netlib";
    } else if (blas_variant == BLAS_ATLAS) {
        return "atlas";
    } else if (blas_variant == BLAS_OPENBLAS) {
        return "openblas";
    } else if (blas_variant == BLAS_MKL) {
        return "mkl";
    } else if (blas_variant == BLAS_VECLIB) {
        return "veclib";
    } else if (blas_variant == BLAS_BLIS) {
        return "blis";
    } else {
        return "unknown";
    }
}

/* for gretl_matrix.c, gretl_mt.c */

int blas_is_threaded (void)
{
    return blas_variant == BLAS_OPENBLAS ||
        blas_variant == BLAS_BLIS ||
        blas_variant == BLAS_MKL;
}

void blas_set_num_threads (int nt)
{
    if (OB_set_num_threads != NULL) {
        OB_set_num_threads(nt);
    } else if (BLIS_set_num_threads != NULL) {
        BLIS_set_num_threads(nt);
    } else if (MKL_domain_set_num_threads != NULL) {
        MKL_domain_set_num_threads(nt, 1); /* MKL_DOMAIN_BLAS */
    }
}

int blas_get_num_threads (void)
{
    if (OB_get_num_threads != NULL) {
        return OB_get_num_threads();
    } else if (BLIS_get_num_threads != NULL) {
        return BLIS_get_num_threads();
    } else if (MKL_domain_get_max_threads != NULL) {
        return MKL_domain_get_max_threads(1); /* MKL_DOMAIN_BLAS */
    } else {
        return 0;
    }
}

static void blas_init (void)
{
    void *ptr = NULL;

#if defined(OS_OSX) && defined(PKGBUILD)
    blas_variant = BLAS_VECLIB; /* the default */
    return;
#else
    blas_variant = BLAS_UNKNOWN;
#endif

    ptr = dlopen(NULL, RTLD_NOW);

    if (ptr != NULL) {
        OB_set_num_threads = dlsym(ptr, "openblas_set_num_threads");
        OB_get_num_threads = dlsym(ptr, "openblas_get_num_threads");
        if (OB_set_num_threads != NULL) {
            blas_variant = BLAS_OPENBLAS;
            register_openblas_details(ptr);
        }
    }

    if (ptr != NULL && blas_variant == BLAS_UNKNOWN) {
        BLIS_init = dlsym(ptr, "bli_init");
        BLIS_finalize = dlsym(ptr, "bli_finalize");
        BLIS_set_num_threads = dlsym(ptr, "bli_thread_set_num_threads");
        BLIS_get_num_threads = dlsym(ptr, "bli_thread_get_num_threads");
        if (BLIS_init != NULL) {
            blas_variant = BLAS_BLIS;
            BLIS_init(); /* This is only to be sure that BLIS was initialized */
            register_blis_details(ptr);
        }
    }

    if (ptr != NULL && blas_variant == BLAS_UNKNOWN) {
        MKL_finalize = dlsym(ptr, "mkl_finalize");
        MKL_domain_get_max_threads = dlsym(ptr, "MKL_Domain_Get_Max_Threads");
        MKL_domain_set_num_threads = dlsym(ptr, "MKL_Domain_Set_Num_Threads");
        if (MKL_domain_set_num_threads != NULL) {
            blas_variant = BLAS_MKL;
            register_mkl_details(ptr);
        }
    }

    if (ptr != NULL && FLAME_init == NULL) {
        FLAME_init = dlsym(ptr, "FLA_Init");
        if (FLAME_init != NULL) {
            FLAME_init();
            libflame_utils(ptr);
        }
    }

    if (blas_variant != BLAS_VECLIB &&
        blas_variant != BLAS_OPENBLAS &&
        blas_variant != BLAS_BLIS &&
        blas_variant != BLAS_MKL) {
#ifdef WIN32
        blas_variant = BLAS_NETLIB; /* ?? */
#else
        blas_variant = detect_blas_via_ldd();
#endif
    }
}

void blas_cleanup (void)
{
    if (blas_variant == BLAS_BLIS && BLIS_finalize != NULL) {
        BLIS_finalize();
    } else if (blas_variant == BLAS_MKL && MKL_finalize != NULL) {
        MKL_finalize();
    }

    if (FLAME_initialized != NULL && FLAME_finalize != NULL) {
        if (FLAME_initialized()) {
            FLAME_finalize();
        }
    }
}

/* CPUID detection and analysis (for $sysinfo bundle) */

#if defined(OS_OSX)
# include <sys/types.h>
# include <sys/sysctl.h>
# define CPU_IDENT
#elif (defined(__x86_64__) || defined(__i386__))
# include <cpuid.h>
# define CPU_IDENT
#endif

#ifdef CPU_IDENT

char *get_cpu_details (void)
{
    char *ret = NULL;
    char vendor_buf[16] = {0};
    char brand_buf[64] = {0};

#if defined(OS_OSX)
    size_t bsz = sizeof brand_buf;

    sysctlbyname("machdep.cpu.brand_string", &brand_buf, &bsz, NULL, 0);
#else
    guint32 i, j, data[4];
    int n_bytes = 4;

    __cpuid(0, data[0], data[1], data[2], data[3]);
    sprintf(vendor_buf, "%.*s%.*s%.*s",
            n_bytes, (char *) &data[1],
            n_bytes, (char *) &data[3],
            n_bytes, (char *) &data[2]);

    __cpuid(0x80000000, data[0], data[1], data[2], data[3]);
    if (data[0] & 0x80000000) {
        if (data[0] >= 0x80000004) {
            for (i=0x80000002; i<=0x80000004; i++) {
                __cpuid(i, data[0], data[1], data[2], data[3]);
                for (j=0; j<4; j++) {
                    strncat(brand_buf, (char *) &data[j], n_bytes);
                }
            }
        }
    }
#endif

    if (*brand_buf) {
        ret = gretl_strdup(g_strstrip(brand_buf));
    } else if (*vendor_buf) {
        ret = gretl_strdup(g_strstrip(vendor_buf));
    } else {
        ret = gretl_strdup("unknown");
    }

    return ret;
}

#else

char *get_cpu_details (void)
{
    return gretl_strdup("unknown");
}

#endif /* CPU_IDENT defined or not */

/* AVX support detection */
#if defined(__x86_64__)
#ifndef __CPUID_H
# include <cpuid.h>
#endif
#define CPU_AVX_DETECT
#endif

int avx_support (void)
{
#if defined(CPU_AVX_DETECT)
    guint32 eax, ebx, ecx, edx;
    int avx = 0, avx2 = 0, avx512 = 0;

    __cpuid(0x80000000, eax, ebx, ecx, edx);
    __cpuid_count(0x00000001, 0, eax, ebx, ecx, edx);
    avx = (ecx & bit_AVX) != 0;

    __cpuid_count(0x00000007, 0, eax, ebx, ecx, edx);
    avx2 = (ebx & bit_AVX2) != 0;
    avx512 = (ebx & bit_AVX512F) != 0;

    if (avx512) {
        return 512;
    } else if (avx2) {
        return 2;
    } else if (avx) {
        return 1;
    } else {
        return 0;
    }
#else
    return 0;
#endif
}

/* library init and cleanup functions */

/**
 * libgretl_init:
 *
 * In a program that uses libgretl, this function should be
 * called once, before any other libgretl functions are
 * used. See also libgretl_cleanup(), and libgretl_mpi_init().
 **/

void libgretl_init (void)
{
    libset_init();
    gretl_rand_init();
    gretl_xml_init();
    gretl_stopwatch_init();
    if (!gretl_in_tool_mode()) {
        blas_init();
        num_threads_init(blas_variant);
    }
#if HAVE_GMP
    mpf_set_default_prec(get_mp_bits());
#endif
#ifdef WIN32
    win32_ensure_path();
#endif
}

#ifdef HAVE_MPI

/**
 * libgretl_mpi_init:
 * @self: the MPI rank of the calling process.
 * @np: the number of MPI processes.
 * @dcmt: if non-zero, set up per-process RNG using DCMT.
 *
 * This function provides an alternative to libgretl_init()
 * which should be used when a libgretl program is to be run in
 * MPI mode.
 **/

int libgretl_mpi_init (int self, int np, int dcmt)
{
    int err;

    libset_init();

    /* let geneval know */
    set_mpi_rank_and_size(self, np);

    /* load MPI library functions */
    err = gretl_MPI_init();
    if (err) {
        return err;
    }

    if (dcmt && np > 1) {
        /* use DCMT for multiple RNGs */
        gretl_dcmt_init(np, self, 0);
    } else {
        /* just use one RNG */
        gretl_rand_init();
    }

    gretl_xml_init();
    gretl_stopwatch_init();
#if HAVE_GMP
    mpf_set_default_prec(get_mp_bits());
#endif

    /* be relatively quiet by default */
    set_gretl_echo(0);
    set_gretl_messages(0);

    return 0;
}

#else

int gretl_mpi_initialized (void)
{
    /* stub for non-MPI build */
    return 0;
}

#endif

static int other_gretl_running (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "rb");
    int ret = 0;

    if (fp != NULL) {
        /* Each entry in gretl.pid will contain two
           integers: PID and sequence number. So if
           we can read more than two values that must
           mean there's more than one gretl process
           running.
        */
        long lv[4];
        int np;

        np = fscanf(fp, "%ld %ld %ld %ld", &lv[0],
                    &lv[1], &lv[2], &lv[3]);
        if (np > 2) {
            ret = 1;
        }
        fclose(fp);
    }

    return ret;
}

static void dotdir_cleanup (void)
{
    const char *dotdir = gretl_dotdir();
    const char *workdir = gretl_workdir();
    int err;

    if (!strcmp(dotdir, workdir)) {
        return;
    }

    err = gretl_chdir(dotdir);

    if (!err) {
        GDir *dir = gretl_opendir(".");
        const gchar *fname;
        int skipit = 0;

        if (dir != NULL) {
            while ((fname = g_dir_read_name(dir)) != NULL) {
                if (!strcmp(fname, "gretl.pid")) {
                    skipit = other_gretl_running(fname);
                    break;
                }
            }
            if (!skipit) {
                g_dir_rewind(dir);
                while ((fname = g_dir_read_name(dir)) != NULL) {
                    if (gretl_isdir(fname)) {
                        if (*fname == '.' && strlen(fname) > 3) {
                            /* failed auto-dot directory? */
                            gretl_deltree(fname);
                        }
                    } else if (!strncmp(fname, "prntmp.", 7) &&
                               !gretl_in_gui_mode()) {
                        ; /* don't let gretlcli delete a GUI file */
                    } else if (strcmp(fname, "..") &&
                               strcmp(fname, ".") &&
                               strcmp(fname, ".gretl2rc") &&
                               strcmp(fname, "gretl.pid") &&
                               strcmp(fname, "addons.idx") &&
                               strcmp(fname, "mail.dat")) {
                        gretl_remove(fname);
                    }
                }
            }
            g_dir_close(dir);
        }
    }
}

void libgretl_session_cleanup (int mode)
{
    /* "last model" for accessors */
    gretl_saved_objects_cleanup();

    /* trash dataset-related items */
    gretl_transforms_cleanup();
    gretl_lists_cleanup();
    gretl_tests_cleanup();
    gretl_plotx(NULL, OPT_NONE);

    /* scrub options set via "setopt" */
    setopt_cleanup();

    if (mode != SESSION_CLEAR_DATASET) {
        destroy_user_vars();
    }
}

/**
 * libgretl_cleanup:
 *
 * In a program that uses libgretl, this function may be
 * called to free various chunks of memory after the program
 * is finished with libgretl. Do not attempt to call any
 * other libgretl functions after invoking this cleanup.
 *
 * See also libgretl_init().
 **/

void libgretl_cleanup (void)
{
    libgretl_session_cleanup(SESSION_CLEAR_ALL);

    gretl_rand_free();
    gretl_functions_cleanup();
    libset_cleanup();
    gretl_command_hash_cleanup();
    gretl_function_hash_cleanup();
    lapack_mem_free();
    forecast_matrix_cleanup();
    stored_options_cleanup();
    option_printing_cleanup();
    gnuplot_cleanup();
    bufgets_cleanup();
    plugins_cleanup();
    gretl_bundle_cleanup();
    gretl_typemap_cleanup();
    gretl_stopwatch_cleanup();
#ifdef USE_CURL
    gretl_www_cleanup();
#endif
    builtin_strings_cleanup();
    last_result_cleanup();

#ifdef HAVE_MPI
    if (!gretl_mpi_initialized()) {
        dotdir_cleanup();
    }
#else
    dotdir_cleanup();
#endif

#ifdef USE_RLIB
    gretl_R_cleanup();
#endif

    gretl_script_dirs_cleanup();
    gretl_xml_cleanup();
    blas_cleanup();
}

/* record and retrieve hypothesis test results */

enum {
    SET_TEST_STAT,
    GET_TEST_STAT,
    GET_TEST_PVAL,
    GET_TEST_LNL,
    GET_TEST_BRK,
    TESTS_CLEANUP
};

#define getcode(c) (c == GET_TEST_STAT || c == GET_TEST_PVAL || \
                    c == GET_TEST_LNL || c == GET_TEST_BRK)

static GretlType last_test_type;

static double record_or_get_test_result (double teststat,
                                         double pval,
                                         double lnl,
                                         double inbrk,
                                         int code)
{
    static double val = NADBL;
    static double pv = NADBL;
    static double ll = NADBL;
    static double brk = NADBL;
    double ret = NADBL;

    if (code == SET_TEST_STAT) {
        last_test_type = GRETL_TYPE_DOUBLE;
        val = teststat;
        pv = pval;
        ll = lnl;
        brk = inbrk;
    } else if (getcode(code) && last_test_type == GRETL_TYPE_DOUBLE) {
        if (code == GET_TEST_STAT) {
            ret = val;
        } else if (code == GET_TEST_PVAL) {
            ret = pv;
        } else if (code == GET_TEST_LNL) {
            ret = ll;
        } else if (code == GET_TEST_BRK) {
            ret = brk;
        }
    }

    return ret;
}

static gretl_matrix *
record_or_get_test_matrix (gretl_matrix *tests,
                           gretl_matrix *pvals,
                           int code, int *err)
{
    static gretl_matrix *vals = NULL;
    static gretl_matrix *pvs = NULL;
    gretl_matrix *ret = NULL;

    if (code == TESTS_CLEANUP) {
        gretl_matrix_free(vals);
        gretl_matrix_free(pvs);
        vals = pvs = NULL;
        last_test_type = GRETL_TYPE_NONE;
        return NULL;
    }

    if (code == SET_TEST_STAT) {
        last_test_type = GRETL_TYPE_MATRIX;
        gretl_matrix_free(vals);
        vals = tests;
        gretl_matrix_free(pvs);
        pvs = pvals;
    } else if (getcode(code)) {
        gretl_matrix *src = (code == GET_TEST_STAT)? vals : pvs;

        if (src != NULL) {
            ret = gretl_matrix_copy(src);
            if (ret == NULL) {
                *err = E_ALLOC;
            }
        } else {
            *err = E_BADSTAT;
        }
    }

    return ret;
}

static void gretl_tests_cleanup (void)
{
    record_or_get_test_matrix(NULL, NULL, TESTS_CLEANUP, NULL);
}

/* Returns the type (scalar, matrix or none) of the last recorded
   test statistic/p-value pair.
*/

int get_last_test_type (void)
{
    return last_test_type;
}

void record_test_result (double teststat, double pval)
{
    record_or_get_test_result(teststat, pval, NADBL, NADBL, SET_TEST_STAT);
}

void record_LR_test_result (double teststat, double pval, double lnl)
{
    record_or_get_test_result(teststat, pval, lnl, NADBL, SET_TEST_STAT);
}

void record_QLR_test_result (double teststat, double pval, double brk)
{
    record_or_get_test_result(teststat, pval, NADBL, brk, SET_TEST_STAT);
}

/* Note: the hypothesis-test recorder "takes ownership" of the
   two input matrices, @tests and @pvals, which should therefore
   NOT be freed by the caller.
*/

void record_matrix_test_result (gretl_matrix *tests,
                                gretl_matrix *pvals)
{
    record_or_get_test_matrix(tests, pvals, SET_TEST_STAT, NULL);
}

double get_last_test_statistic (void)
{
    return record_or_get_test_result(0, 0, 0, 0, GET_TEST_STAT);
}

double get_last_pvalue (void)
{
    return record_or_get_test_result(0, 0, 0, 0, GET_TEST_PVAL);
}

double get_last_lnl (void)
{
    return record_or_get_test_result(0, 0, 0, 0, GET_TEST_LNL);
}

double get_last_break (void)
{
    return record_or_get_test_result(0, 0, 0, 0, GET_TEST_BRK);
}

gretl_matrix *get_last_test_matrix (int *err)
{
    return record_or_get_test_matrix(NULL, NULL, GET_TEST_STAT, err);
}

gretl_matrix *get_last_pvals_matrix (int *err)
{
    return record_or_get_test_matrix(NULL, NULL, GET_TEST_PVAL, err);
}

#define NEED_ALIGNED_MALLOC 0

#if NEED_ALIGNED_MALLOC

/* We needed this -- in the absence of posix_memalign -- when
   experimenting with the array variant of SFMT (random.c). Since that
   did not seem to be faster than the "one at a time" variant when we
   last experimented, we junked the array-using code, which makes the
   following functionality redundant. However, we may want to try
   using SFMT arrays again. In that case, see git commit
   65f3395bbdb01450496fb20fcf0f8f15092f9eeb for the last revision of
   random.c that included the array code.
*/

/*
  malloc and free for alignments greater than that guaranteed by the C
  library, based on Steven G. Johnson's public domand code at
  http://ab-initio.mit.edu/~stevenj/align.c
*/

#ifdef HAVE_STDINT_H
# include <stdint.h> /* for uintptr_t */
#else
# define uintptr_t size_t
#endif

#define NOT_POWER_OF_TWO(n) (((n) & ((n) - 1)))
#define UI(p) ((uintptr_t) (p))

#define PTR_ALIGN(p0, align)                            \
    ((void *) (((UI(p0) + (align + sizeof(void*)))      \
                & (~UI(align - 1)))))

/* pointer must sometimes be aligned; assume sizeof(void*) is a power of two */
#define ORIG_PTR(p) (*(((void **) (UI(p) & (~UI(sizeof(void*) - 1)))) - 1))

static void *real_aligned_malloc (size_t size, size_t align)
{
    void *p0, *p;

    if (NOT_POWER_OF_TWO(align)) {
        errno = EINVAL;
        return NULL;
    }

    if (align < sizeof(void *)) {
        align = sizeof(void *);
    }

    /* including the extra sizeof(void*) is overkill on a 32-bit
       machine, since malloc is already 8-byte aligned, as long
       as we enforce alignment >= 8 ...but oh well */

    p0 = malloc(size + align + sizeof(void *));
    if (!p0) {
        return NULL;
    }

    p = PTR_ALIGN(p0, align);
    ORIG_PTR(p) = p0;

    return p;
}

void *gretl_aligned_malloc (size_t size, size_t alignment)
{
    if (size < 1) {
        return NULL;
    } else {
        return real_aligned_malloc(size, alignment);
    }
}

void gretl_aligned_free (void *mem)
{
    if (mem != NULL) {
        free(ORIG_PTR(mem));
    }
}

#endif /* NEED_ALIGNED_MALLOC */

#ifdef WIN32

int check_for_program (const char *prog)
{
    int ret = 0;

    if (prog == NULL || *prog == '\0') {
        return 0;
    }

    if (has_suffix(prog, ".exe")) {
        ret = win32_check_for_program(prog);
    } else {
        gchar *test = g_strdup_printf("%s.exe", prog);

        ret = win32_check_for_program(test);
        g_free(test);
    }

    return ret;
}

#else /* !WIN32 */

#include <sys/stat.h>
#include <unistd.h>

static int is_executable (const char *s, uid_t myid, gid_t mygrp)
{
    struct stat buf;
    int ok = 0;

    if (gretl_stat(s, &buf) != 0) {
        return 0;
    }

    if (buf.st_mode & (S_IFREG|S_IFLNK)) {
        if (buf.st_uid == myid && (buf.st_mode & S_IXUSR)) {
            ok = 1;
        } else if (buf.st_gid == mygrp && (buf.st_mode & S_IXGRP)) {
            ok = 1;
        } else if (buf.st_uid != myid && buf.st_gid != mygrp &&
                   (buf.st_mode & S_IXOTH)) {
            ok = 1;
        }
    }

    return ok;
}

/**
 * check_for_program:
 * @prog: name of program.
 *
 * Returns: 1 if @prog is found in the PATH, otherwise 0.
 */

int check_for_program (const char *prog)
{
    uid_t myid = getuid();
    gid_t mygrp = getgid();
    char *path;
    char *pathcpy;
    char **dirs;
    char *fullpath;
    char *p;
    int max_dlen = 0;
    int found = 0;
    int i, ndirs;

    if (prog == NULL || *prog == '\0') {
        return 0;
    }

    if (*prog == '/') {
        return is_executable(prog, myid, mygrp);
    }

    path = getenv("PATH");
    if (path == NULL || *path == '\0') {
        return 0;
    }

    pathcpy = gretl_strdup(path);
    if (pathcpy == NULL) {
        return 0;
    }

    ndirs = 1;
    p = pathcpy;
    while (*p) {
        if (*p == ':') ndirs++;
        p++;
    }

    dirs = malloc(ndirs * sizeof *dirs);
    if (dirs == NULL) {
        free(pathcpy);
        return 0;
    }

    if (ndirs == 1) {
        dirs[0] = pathcpy;
        max_dlen = strlen(pathcpy);
    } else {
        for (i=0; i<ndirs; i++) {
            int dlen;

            dirs[i] = strtok((i == 0)? pathcpy : NULL, ":");
            if (dirs[i] == NULL) {
                ndirs = i;
                break;
            }
            dlen = strlen(dirs[i]);
            if (dlen > max_dlen) {
                max_dlen = dlen;
            }
        }
    }

    if (ndirs == 0 ||
        (fullpath = malloc(max_dlen + strlen(prog) + 2)) == NULL) {
        free(dirs);
        free(pathcpy);
        return 0;
    }

    for (i=0; i<ndirs && !found; i++) {
        sprintf(fullpath, "%s/%s", dirs[i], prog);
        found = is_executable(fullpath, myid, mygrp);
    }

    free(dirs);
    free(pathcpy);
    free(fullpath);

    return found;
}

#endif /* WIN32 or not */

gint64 gretl_monotonic_time (void)
{
    return g_get_monotonic_time();
}

/* implementations of gzip and gunzip for unitary files */

#define GRETL_ZBUFSIZ 131072 /* 128k bytes */

int gretl_gzip (char *fname, char *zname)
{
    char buf[GRETL_ZBUFSIZ];
    size_t len;
    FILE *fp;
    gzFile fz;

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
        return E_FOPEN;
    }

    fz = gretl_gzopen(zname, "wb");
    if (fz == Z_NULL) {
        fclose(fp);
        return E_FOPEN;
    }

    while ((len = fread(buf, 1, GRETL_ZBUFSIZ, fp)) > 0) {
        gzwrite(fz, buf, len);
    }

    fclose(fp);
    gzclose(fz);

    return 0;
}

int gretl_gunzip (char *zname, char *fname)
{
    char buf[GRETL_ZBUFSIZ];
    size_t len;
    FILE *fp;
    gzFile fz;

    fz = gretl_gzopen(zname, "rb");
    if (fz == Z_NULL) {
        return E_FOPEN;
    }

    fp = gretl_fopen(fname, "wb");
    if (fp == NULL) {
        gzclose(fz);
        return E_FOPEN;
    }

    while ((len = gzread(fz, buf, GRETL_ZBUFSIZ)) > 0) {
        fwrite(buf, 1, len, fp);
    }

    gzclose(fz);
    fclose(fp);

    return 0;
}

/* Check if we're able to launch an automatic MPI routine, on the
   local machine, from within some time-consuming function.  This
   requires that MPI is available but not already running, and
   the local machine has at least 2 processors.
*/

int auto_mpi_ok (void)
{
    int ret = 0;

#ifdef HAVE_MPI
    if (gretl_mpi_initialized()) {
        ; /* No: can't run MPI under MPI */
    } else if (gretl_n_processors() < 2) {
        ; /* No: can't do local MPI */
    } else {
        /* Yes, if mpiexec is installed */
        ret = check_for_mpiexec();
    }
#endif

    return ret;
}

gretl_matrix *dec2bin (double x, const gretl_matrix *v, int *err)
{
    gretl_matrix *ret = NULL;
    double *val;
    guint32 ui;
    int n = 1;
    int i, j;

    if (v != NULL) {
        n = gretl_vector_get_length(v);
        if (n == 0) {
            *err = E_INVARG;
            return NULL;
        }
        val = v->val;
    } else {
        val = &x;
    }

    ret = gretl_zero_matrix_new(n, 32);
    if (ret == NULL) {
        *err = E_ALLOC;
    }

    for (i=0; i<n; i++) {
        ui = gretl_unsigned_from_double(val[i], err);
        if (*err) {
            break;
        }
        j = 0;
        while (ui > 0) {
            gretl_matrix_set(ret, i, j++, ui % 2);
            ui = ui >> 1;
        }
    }

    if (*err) {
        gretl_matrix_free(ret);
        ret = NULL;
    }

    return ret;
}

gretl_matrix *bin2dec (const gretl_matrix *m, int *err)
{
    int r = gretl_matrix_rows(m); /* must be 1 or more */
    int c = gretl_matrix_cols(m); /* must be >= 1 && <= 32 */
    gretl_matrix *ret = NULL;

    if (r == 0 || c == 0 || c > 32) {
        *err = E_INVARG;
    } else {
        /* allocate column vector for return */
        ret = gretl_zero_matrix_new(r, 1);
        if (ret == NULL) {
            *err = E_ALLOC;
        }
    }

    if (!*err) {
        guint32 ui, k;
        double mij;
        int i, j;

        for (i=0; i<r && !*err; i++) {
            ui = 0;
            k = 0x01;
            for (j=0; j<c && !*err; j++) {
                mij = gretl_matrix_get(m, i, j);
                if (mij == 1) {
                    ui += k;
                } else if (mij != 0) {
                    *err = E_INVARG;
                    break;
                }
                k = k << 1;
            }
            ret->val[i] = (double) ui;
        }
    }

    return ret;
}

/* apparatus for combining categorical indices */

typedef struct {
    int a;
    int b;
    int pos;
} tri_sorter;

static int compare_ts (const void *a, const void *b)
{
    const tri_sorter *ts1 = a;
    const tri_sorter *ts2 = b;
    int ret = ts1->a - ts2->a;

    if (ret == 0) {
	ret = ts1->b - ts2->b;
    }

    return ret;
}

/* Given discrete series @v0 and @v1, fill series @v2 with
   indices of the combination of values of the first two
   series.
*/

int combine_categories (DATASET *dset, int v0, int v1, int v2)
{
    tri_sorter *tsr;
    int i, k;

    tsr = malloc(dset->n * sizeof *tsr);
    if (tsr == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<dset->n; i++) {
	tsr[i].a = dset->Z[v0][i];
	tsr[i].b = dset->Z[v1][i];
	tsr[i].pos = i;
    }

    qsort(tsr, dset->n, sizeof *tsr, compare_ts);

    dset->Z[v2][tsr[0].pos] = k = 1;
    for (i=1; i<dset->n; i++) {
	if (tsr[i].b != tsr[i-1].b || tsr[i].a != tsr[i-1].a) {
	    k++;
	}
	dset->Z[v2][tsr[i].pos] = k;
    }

#if 0 /* debugging */
    int list[4] = {3, v0, v1, v2};
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
    printdata(list, NULL, dset, OPT_O, prn);
    gretl_print_destroy(prn);
#endif

    free(tsr);

    return 0;
}
