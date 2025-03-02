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

/* various functions called by 'genr' */

#include "libgretl.h"
#include "genparse.h"
#include "gen_public.h"
#include "libset.h"
#include "monte_carlo.h"
#include "gretl_panel.h"
#include "gretl_cmatrix.h"
#include "gretl_typemap.h"
#include "gretl_midas.h"
#include "estim_private.h"
#include "matrix_extra.h"
#include "gretl_btree.h"
#include "gretl_func.h"
#include "kalman.h"
#include "../../cephes/cephes.h"

#include <errno.h>

int sort_series (const double *x, double *y, int f,
                 const DATASET *dset)
{
    double *z = NULL;
    int n = sample_size(dset);
    int i, t;

    z = malloc(n * sizeof *z);
    if (z == NULL) {
        return E_ALLOC;
    }

    i = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
        if (!na(x[t])) {
            z[i++] = x[t];
        }
    }

    if (f == F_DSORT) {
        qsort(z, i, sizeof *z, gretl_inverse_compare_doubles);
    } else {
        qsort(z, i, sizeof *z, gretl_compare_doubles);
    }

    i = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
        if (na(x[t])) {
            y[t] = NADBL;
        } else {
            y[t] = z[i++];
        }
    }

    free(z);

    return 0;
}

struct pair_sorter {
    double x;
    double y;
};

/* sort the series y by the values of x, putting the result
   into z */

int gretl_sort_by (const double *x, const double *y,
                   double *z, const DATASET *dset)
{
    struct pair_sorter *xy;
    int n = sample_size(dset);
    int i, t;

    for (t=dset->t1; t<=dset->t2; t++) {
        if (na(x[t])) {
            return E_MISSDATA;
        }
    }

    xy = malloc(n * sizeof *xy);
    if (xy == NULL) {
        return E_ALLOC;
    }

    i = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
        xy[i].x = x[t];
        xy[i].y = y[t];
        i++;
    }

    qsort(xy, n, sizeof *xy, gretl_compare_doubles);

    i = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
        z[t] = xy[i++].y;
    }

    free(xy);

    return 0;
}

static void genrank (const double *sz, int m,
                     const double *z, int n,
                     double *rz)
{
    int cases, k, i, j;
    double avg, r = 1;

    for (i=0; i<m; i++) {
        /* scan sorted z */
        cases = k = 0;

        if (i > 0 && sz[i] == sz[i-1]) {
            continue;
        }

        for (j=0; j<n; j++) {
            /* scan raw z for matches */
            if (!na(z[j])) {
                if (z[j] == sz[i]) {
                    rz[k] = r;
                    cases++;
                }
                k++;
            }
        }

        if (cases > 1) {
            avg = (r + r + cases - 1.0) / 2.0;
            for (j=0; j<m; j++) {
                if (rz[j] == r) {
                    rz[j] = avg;
                }
            }
        }

        r += cases;
    }
}

/* implements both rank_series() and rank_vector() */

static int rank_array (const double *x, double *y, int f, int n)
{
    double *sx = NULL;
    double *rx = NULL;
    int m = n;
    int i, t;

    for (t=0; t<n; t++) {
        if (na(x[t])) m--;
    }

    sx = malloc(m * sizeof *sx);
    rx = malloc(m * sizeof *rx);

    if (sx == NULL || rx == NULL) {
        free(sx);
        free(rx);
        return E_ALLOC;
    }

    i = 0;
    for (t=0; t<n; t++) {
        if (!na(x[t])) {
            sx[i] = x[t];
            rx[i] = 0.0;
            i++;
        }
    }

    if (f == F_DSORT) {
        qsort(sx, m, sizeof *sx, gretl_inverse_compare_doubles);
    } else {
        qsort(sx, m, sizeof *sx, gretl_compare_doubles);
    }

    genrank(sx, m, x, n, rx);

    i = 0;
    for (t=0; t<n; t++) {
        if (na(x[t])) {
            y[t] = NADBL;
        } else {
            y[t] = rx[i++];
        }
    }

    free(sx);
    free(rx);

    return 0;
}

int rank_series (const double *x, double *y, int f,
                 const DATASET *dset)
{
    return rank_array(x + dset->t1, y + dset->t1,
                      f, sample_size(dset));
}

gretl_matrix *rank_vector (const gretl_matrix *x, int f, int *err)
{
    gretl_matrix *y = NULL;
    int n = gretl_vector_get_length(x);

    if (n == 0) {
        *err = E_DATA;
    } else {
        y = gretl_matrix_alloc(x->rows, x->cols);
        if (y == NULL) {
            *err = E_ALLOC;
        }
    }

    if (!*err) {
        *err = rank_array(x->val, y->val, f, n);
        if (*err) {
            gretl_matrix_free(y);
            y = NULL;
        }
    }

    return y;
}

/**
 * diff_series:
 * @x: array of original data.
 * @y: array into which to write the result.
 * @f: function, F_DIFF, F_SDIFF or F_LDIFF.
 * @dset: data set information.
 *
 * Calculates the differenced counterpart to the input
 * series @x.  If @f = F_SDIFF, the seasonal difference is
 * computed; if @f = F_LDIFF, the log difference, and if
 * @f = F_DIFF, the ordinary first difference.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int diff_series (const double *x, double *y, int f,
                 const DATASET *dset)
{
    int k = (f == F_SDIFF)? dset->pd : 1;
    int t, t1 = dset->t1;

    if (t1 < k) {
        t1 = k;
    }

    for (t=t1; t<=dset->t2; t++) {
        if (dataset_is_panel(dset) && t % dset->pd == 0) {
            /* skip first observation in panel unit */
            continue;
        }
        if (na(x[t]) || na(x[t-k])) {
            continue;
        }
        if (f == F_DIFF || f == F_SDIFF) {
            y[t] = x[t] - x[t-k];
        } else if (f == F_LDIFF) {
            if (x[t] > 0.0 && x[t-k] > 0.0) {
                y[t] = log(x[t]) - log(x[t-k]);
            }
        }
    }

    return 0;
}

int two_series_minmax (const double *x, const double *y,
                       double *z, int f, const DATASET *dset)
{
    int t;

    for (t=dset->t1; t<=dset->t2; t++) {
        if (na(x[t]) || na(y[t])) {
            continue;
        }
        if (f == F_MIN) {
            z[t] = MIN(x[t], y[t]);
        } else {
            z[t] = MAX(x[t], y[t]);
        }
    }

    return 0;
}

gretl_matrix *two_matrices_minmax (const gretl_matrix *a,
                                   const gretl_matrix *b,
                                   int f, int *err)
{
    gretl_matrix *ret = NULL;

    if (a->rows != b->rows || a->cols != b->cols) {
        *err = E_NONCONF;
        return NULL;
    }

    ret = gretl_matrix_alloc(a->rows, a->cols);

    if (ret == NULL) {
        *err = E_ALLOC;
    } else {
        int i, n = a->rows * a->cols;
        double ai, bi;

        for (i=0; i<n; i++) {
            ai = a->val[i];
            bi = b->val[i];
            if (na(ai) || na(bi)) {
                ret->val[i] = NADBL;
            } else if (f == F_MIN) {
                ret->val[i] = MIN(ai, bi);
            } else {
                ret->val[i] = MAX(ai, bi);
            }
        }
    }

    return ret;
}

/* Special case of interpolate_series() for panel data. */

static int panel_interpolate (const double *x, double *y,
                              const DATASET *dset)
{
    int i, t, s, dt;
    int i1 = dset->t1 / dset->pd;
    int i2 = dset->t2 / dset->pd;
    int t1, t2, tt, tt0, tt2, tnext;
    double b;

    for (i=i1; i<=i2; i++) {
        /* process the data for unit @i */
        tt0 = i * dset->pd;
        t1 = 0;
        t2 = dset->pd - 1;
        for (t=0; t<=t2; t++) {
            s = tt0 + t;
            y[s] = x[s];
            if (na(x[s])) {
                t1++;
            } else {
                break;
            }
        }
        for (t=t2; t>=t1; t--) {
            s = tt0 + t;
            y[s] = x[s];
            if (na(x[s])) {
                t2--;
            } else {
                break;
            }
        }
        /* having found the first and last non-missing
           observations, perform interpolation
        */
        tt2 = tt0 + t2;
        for (t=t1; t<=t2; t++) {
            tt = tt0 + t;
            dt = 0;
            for (s=tt+1; s<=tt2; s++) {
                if (na(x[s])) {
                    dt++;
                } else {
                    break;
                }
            }
            y[tt] = x[tt];
            if (dt > 0) {
                tnext = s - 1 - tt0;
                b = (x[s] - x[tt]) / (dt+1);
                for (s=1; s<=dt; s++) {
                    y[tt+s] = x[tt] + b*s;
                }
                t = tnext;
            }
        }
    }

    return 0;
}

/**
 * interpolate_series:
 * @x: array of original data.
 * @y: array into which to write the result.
 * @dset: data set information.
 *
 * Performs linear interpolation of missing values in @x,
 * writing the result into @y. Panel data are handled:
 * interpolation is strictly in the time-series dimension.
 * No extrapolation is performed.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int interpolate_series (const double *x, double *y,
                        const DATASET *dset)
{
    int t, s, dt;
    int t1 = dset->t1;
    int t2 = dset->t2;
    int tnext;
    double b;

    if (dataset_is_panel(dset)) {
        return panel_interpolate(x, y, dset);
    } else if (!dataset_is_time_series(dset)) {
        gretl_errmsg_set(_("This function requires time-series or panel data"));
        return E_PDWRONG;
    }

    /* first determine the first and last non-missing values */
    for (t=dset->t1; t<=dset->t2; t++) {
        y[t] = x[t];
        if (na(x[t])) {
            t1++;
        } else {
            break;
        }
    }
    for (t=dset->t2; t>=t1; t--) {
        y[t] = x[t];
        if (na(x[t])) {
            t2--;
        } else {
            break;
        }
    }

    /* perform linear interpolation */
    for (t=t1; t<=t2; t++) {
        dt = 0;
        for (s=t+1; s<=t2; s++) {
            if (na(x[s])) {
                dt++;
            } else {
                break;
            }
        }
        y[t] = x[t];
        if (dt > 0) {
            tnext = s-1;
            b = (x[s] - x[t]) / (dt+1);
            for (s=1; s<=dt; s++) {
                y[t+s] = x[t] + b*s;
            }
            t = tnext;
        }
    }

    return 0;
}

/**
 * standardize_series:
 * @x: array of original data.
 * @y: array into which to write the result.
 * @dfc: degrees of freedom correction.
 * @dset: data set information.
 *
 * By default calculates the standardized counterpart to the input
 * series @x, but if @dfc < 0 the result is just centered.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int standardize_series (const double *x, double *y, int dfc,
                        const DATASET *dset)
{
    double d, xbar = 0;
    int t, n = 0;

    for (t=dset->t1; t<=dset->t2; t++) {
        if (!na(x[t])) {
            xbar += x[t];
            n++;
        }
    }

    if (dfc >= 0 && n < dfc + 1) {
        return E_TOOFEW;
    }

    xbar /= n;

    if (dfc < 0) {
        /* just centering */
        for (t=dset->t1; t<=dset->t2; t++) {
            if (na(x[t])) {
                y[t] = NADBL;
            } else {
                y[t] = x[t] - xbar;
            }
        }
    } else {
        /* dividing by s.d. */
        double sd, TSS = 0;

        for (t=dset->t1; t<=dset->t2; t++) {
            if (!na(x[t])) {
                d = x[t] - xbar;
                TSS += d * d;
            }
        }

        sd = sqrt(TSS / (n - dfc));

        for (t=dset->t1; t<=dset->t2; t++) {
            if (na(x[t])) {
                y[t] = NADBL;
            } else {
                y[t] = (x[t] - xbar) / sd;
            }
        }
    }

    return 0;
}

/**
 * orthdev_series:
 * @x: array of original data.
 * @y: array into which to write the result.
 * @dset: data set information.
 *
 * Calculates in @y the forward orthogonal deviations of the input
 * series @x.  That is, y[t] is the scaled difference between x[t]
 * and the mean of the subsequent observations on x.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int orthdev_series (const double *x, double *y, const DATASET *dset)
{
    double xbar;
    int n, s, t, Tt;

    for (t=dset->t1; t<dset->t2; t++) {

        if (na(x[t])) {
            continue;
        }

        if (dataset_is_panel(dset)) {
            Tt = dset->pd - (t % dset->pd) - 1;
        } else {
            Tt = dset->t2 - t;
        }

        xbar = 0.0;
        n = 0;
        for (s=1; s<=Tt; s++) {
            if (!na(x[t+s])) {
                xbar += x[t+s];
                n++;
            }
        }

        if (n > 0) {
            xbar /= n;
            /* Lead one period, for compatibility with first diffs.
               I.e. we lose the first observation rather than the
               last.  This is for arbond.  Cf. Doornik, Bond and
               Arellano, DPD documentation.
            */
            y[t+1] = sqrt(n / (n + 1.0)) * (x[t] - xbar);
        }
    }

    return 0;
}

/**
 * fracdiff_series:
 * @x: array of original data.
 * @y: array into which to write the result.
 * @d: fraction by which to difference.
 * @diff: boolean variable 1 for fracdiff, 0 for fraclag
 * @obs: used for autoreg calculation, -1 if whole series
 *      should be calculated otherwise just the observation for
 *      obs is calculated
 * @dset: data set information.
 *
 * Calculates the fractionally differenced or lagged
 * counterpart to the input series @x. The fractional
 * difference operator is defined as (1-L)^d, while the
 * fractional lag operator 1-(1-L)^d.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int fracdiff_series (const double *x, double *y, double d,
                     int diff, int obs, const DATASET *dset)
{
    int dd, t, T;
    const double TOL = 1.0E-12;
    int t1 = dset->t1;
    int t2 = (obs >= 0)? obs : dset->t2;
    int tmiss = 0;
    double phi = (diff)? -d : d;

    if (na(d)) {
        return E_NAN;
    }

#if 0
    fprintf(stderr, "Doing fracdiff_series, with d = %g\n", d);
#endif

    if (series_adjust_sample(x, &t1, &t2) != 0) {
        tmiss = first_missing_index(x, t1, t2);
        if (tmiss > 0 && tmiss < t2) {
            t2 = tmiss;
        }
    }

    if (obs >= 0) {
        /* doing a specific observation */
        T = obs - t1 + 1;
        for (t=0; t<dset->n; t++) {
            y[t] = NADBL;
        }
        if (obs != t1) {
            y[obs] = (diff)? x[obs] : 0;
            for (dd=1; dd<T && fabs(phi)>TOL; dd++) {
                y[obs] += phi * x[obs - dd];
                phi *= (dd - d)/(dd + 1);
            }
        } else if (diff) {
            y[obs] = x[obs];
        }
    } else {
        /* doing the whole series */
        T = t2 - t1 + 1;
        for (t=0; t<=t2; t++) {
            if (t >= t1 && t <= t2) {
                y[t] = (diff)? x[t] : 0;
            } else {
                y[t] = NADBL;
            }
        }
        for (dd=1; dd<=T && fabs(phi)>TOL; dd++) {
            for (t=t1+dd; t<=t2; t++) {
                y[t] += phi * x[t - dd];
            }
            phi *= (dd - d)/(dd + 1);
        }
    }

    return 0;
}

static int boxcox_vector (const double *x, double *y, double d,
                          int t1, int t2)
{
    int t, special_case = 999;
    double tol = 1.0e-12;

    /* 999 for general case; here we want to avoid
       the pow() function whenever possible.
    */

    if (fabs(d) < tol) {
        special_case = 0;
    } else if (fabs(1.0 - d) < tol) {
        special_case = 1;
    } else if (fabs(-1.0 - d) < tol) {
        special_case = -1;
    } else if (fabs(0.5 - d) < tol) {
        special_case = 2;
    } else if (fabs(-0.5 - d) < tol) {
        special_case = -2;
    }

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            y[t] = NADBL;
        } else if (special_case == -2) {
            y[t] = (x[t] > 0)? 2.0 * (1.0 - 1.0/sqrt(x[t])) : NADBL;
        } else if (special_case == -1) {
            y[t] = 1.0 - 1.0/x[t];
        } else if (special_case == 0) {
            y[t] = (x[t] > 0)? log(x[t]) : NADBL;
        } else if (special_case == 1) {
            y[t] = x[t] - 1.0;
        } else if (special_case == 2) {
            y[t] = (x[t] > 0)? 2.0 * (sqrt(x[t]) - 1.0) : NADBL;
        } else {
            y[t] = (pow(x[t], d) - 1) / d;
        }
    }

    return 0;
}

/**
 * boxcox_series:
 * @x: array of original data.
 * @y: array into which to write the result.
 * @d: lambda parameter.
 * @dset: data set information.
 *
 * Calculates in @y the Box-Cox transformation for the
 * input series @x.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int boxcox_series (const double *x, double *y, double d,
                   const DATASET *dset)
{
    return boxcox_vector(x, y, d, dset->t1, dset->t2);
}

/**
 * boxcox_matrix:
 * @m: matrix of original data.
 * @d: lambda parameter.
 *
 * Returns: a matrix whose columns are Box-Cox transformations
 * of the columns of @m.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

gretl_matrix *boxcox_matrix (const gretl_matrix *m, double d,
                             int *err)
{
    gretl_matrix *bc;

    bc = gretl_matrix_alloc(m->rows, m->cols);

    if (bc == NULL) {
        *err = E_ALLOC;
    } else {
        const double *x = m->val;
        double *y = bc->val;
        int j, n = m->rows;

        for (j=0; j<m->cols; j++) {
            boxcox_vector(x, y, d, 0, n-1);
            x += n;
            y += n;
        }
    }

    return bc;
}

int cum_series (const double *x, double *y, const DATASET *dset)
{
    int t, s = dset->t1;

    /* corner case: sample size of 1! */
    if (dset->t1 == dset->t2) {
        y[dset->t1] = x[dset->t1];
        return 0;
    }

    for (t=dset->t1; t<=dset->t2 && na(x[t]); t++) {
        s++;
    }

    if (s == dset->t2) {
        /* no usable data */
        return 0;
    }

    y[s] = x[s];

    if (dataset_is_panel(dset)) {
        int k;

        for (t=s+1; t<=dset->t2; t++) {
            if (t % dset->pd == 0) {
                /* first obs for panel unit */
                for (k=0; k<dset->pd; k++) {
                    if (!na(x[t+k])) {
                        t = t + k;
                        y[t] = x[t];
                        break;
                    }
                }
            } else if (!na(y[t-1]) && !na(x[t])) {
                y[t] = y[t-1] + x[t];
            }
        }
    } else {
        for (t=s+1; t<=dset->t2 && !na(x[t]); t++) {
            y[t] = y[t-1] + x[t];
        }
    }

    return 0;
}

static int panel_resample_series (const double *x, double *y,
                                  const DATASET *dset)
{
    int n = panel_sample_size(dset);
    int T = dset->pd;
    int i, j, t, s;
    int u1, u2;
    int *z;

    if (n <= 1) {
        return E_DATA;
    }

    z = malloc(n * sizeof *z);
    if (z == NULL) {
        return E_ALLOC;
    }

    u1 = dset->t1 / T;
    u2 = dset->t2 / T;

    /* select n units from [u1 .. u2] */
    gretl_rand_int_minmax(z, n, u1, u2);

    s = dset->t1;
    for (i=0; i<n; i++) {
        j = z[i] * T;
        for (t=0; t<T; t++) {
            y[s++] = x[j + t];
        }
    }

    free(z);

    return 0;
}

int resample_series (const double *x, double *y, const DATASET *dset)
{
    int *z = NULL;
    int t1, t2;
    int i, t, n;

    if (dataset_is_panel(dset)) {
        return panel_resample_series(x, y, dset);
    }

    t1 = dset->t1;
    t2 = dset->t2;
    series_adjust_sample(x, &t1, &t2);

    n = t2 - t1 + 1;
    if (n <= 1) {
        return E_DATA;
    }

    z = malloc(n * sizeof *z);
    if (z == NULL) {
        return E_ALLOC;
    }

    /* generate n uniform drawings from [t1 .. t2] */
    gretl_rand_int_minmax(z, n, t1, t2);

    /* sample from source series x based on indices z */
    for (t=t1, i=0; t<=t2; t++, i++) {
        y[t] = x[z[i]];
    }

    free(z);

    return 0;
}

int block_resample_series (const double *x, double *y, int blocklen,
                           const DATASET *dset)
{
    int t1 = dset->t1;
    int t2 = dset->t2;
    int *z = NULL;
    int m, rem, bt2, x0;
    int i, s, t, n;

    if (dataset_is_panel(dset)) {
        return E_PDWRONG;
    }

    if (blocklen <= 0) {
        return E_DATA;
    }

    if (blocklen == 1) {
        return resample_series(x, y, dset);
    }

    series_adjust_sample(x, &t1, &t2);

    n = t2 - t1 + 1;

    m = n / blocklen;
    rem = n % blocklen;

    /* Let n now represent the number of blocks of @blocklen
       contiguous observations which we need to select; the
       last of these may not be fully used.
    */
    n = m + (rem > 0);

    /* the last selectable starting point for a block */
    bt2 = t2 - blocklen + 1;

    if (bt2 < t1) {
        return E_DATA;
    }

    z = malloc(n * sizeof *z);
    if (z == NULL) {
        return E_ALLOC;
    }

    /* Generate uniform random series: we want n drawings from the
       range [t1 .. t2 - blocklen + 1], each of which will be
       interpreted as the starting point of a block to be used.
    */
    gretl_rand_int_minmax(z, n, t1, bt2);

    /* Sample from the source series using blocks given by the random
       indices: note that the last block will be incomplete if rem > 0
    */
    t = t1;
    for (i=0; i<n; i++) {
        x0 = z[i];
        for (s=0; s<blocklen; s++) {
            if (t <= t2) {
                y[t++] = x[x0+s];
            } else {
                break;
            }
        }
    }

    free(z);

    return 0;
}

/* retrieve a pre-sample @x value if it's available via the n-vector
   @x0, an optional argument to the filter() function.
*/

static double get_xlag (int lag, int t1, gretl_vector *x0, int n)
{
    int p = t1 - lag;
    int i = n - p;

    if (i >= 0 && i < n) {
        return x0->val[i];
    } else {
        return 0;
    }
}

/* implements filter_series() and filter_matrix() */

static int filter_vector (const double *x, double *y, int t1, int t2,
                          gretl_vector *A, gretl_vector *C, double y0,
                          gretl_vector *x0)
{
    int t, s, i, n;
    int amax, cmax;
    double xlag, ylag;
    double coef, *e;
    int x0len = 0;
    int err = 0;

    if (gretl_is_null_matrix(C)) {
        cmax = 0;
    } else {
        cmax = gretl_vector_get_length(C);
        if (cmax == 0) {
            /* if present, C must be a vector */
            return E_NONCONF;
        }
    }

    if (gretl_is_null_matrix(A)) {
        amax = 0;
    } else {
        amax = gretl_vector_get_length(A);
        if (amax == 0) {
            /* if present, A must be a vector */
            return E_NONCONF;
        }
    }

    n = t2 - t1 + 1;
    e = malloc(n * sizeof *e);

    if (e == NULL) {
        return E_ALLOC;
    }

    if (x0 != NULL) {
        x0len = gretl_vector_get_length(x0);
    }

    s = 0;
    if (cmax) {
        for (t=t1; t<=t2; t++) {
            e[s] = 0;
            for (i=0; i<cmax; i++) {
                if (t-i >= t1) {
                    xlag = x[t-i];
                } else if (x0 == NULL) {
                    xlag = 0;
                } else {
                    xlag = get_xlag(t-i, t1, x0, x0len);
                }
                if (na(xlag)) {
                    e[s] = NADBL;
                    break;
                } else {
                    coef = gretl_vector_get(C, i);
                    e[s] += xlag * coef;
                }
            }
            s++;
        }
    } else {
        /* implicitly MA(0) = 1 */
        for (t=t1; t<=t2; t++) {
            e[s++] = x[t];
        }
    }

    s = 0;
    if (amax) {
        for (t=t1; t<=t2; t++) {
            if (na(e[s])) {
                y[t] = NADBL;
            } else {
                y[t] = e[s];
                for (i=0; i<amax; i++) {
                    ylag = (t-i > t1)? y[t-i-1] : y0;
                    if (na(ylag)) {
                        y[t] = NADBL;
                        break;
                    } else {
                        coef = gretl_vector_get(A, i);
                        y[t] += coef * ylag;
                    }
                }
            }
            s++;
        }
    } else {
        for (t=t1; t<=t2; t++) {
            y[t] = e[s++];
        }
    }

    free(e);

    return err;
}

/**
 * filter_series:
 * @x: array of original data.
 * @y: array into which to write the result.
 * @dset: data set information.
 * @A: vector for autoregressive polynomial.
 * @C: vector for moving average polynomial.
 * @y0: initial value of output series.
 * @x0: prior values of x.
 *
 * Filters x according to y_t = C(L)/A(L) x_t.  If the intended
 * AR order is p, @A should be a vector of length p.  If the
 * intended MA order is q, @C should be vector of length (q+1),
 * the first entry giving the coefficient at lag 0.  However, if
 * @C is NULL this is taken to mean that the lag-0 MA coefficient
 * is unity (and all others are zero).
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int filter_series (const double *x, double *y, const DATASET *dset,
                   gretl_vector *A, gretl_vector *C,
                   double y0, gretl_vector *x0)
{
    int t1 = dset->t1;
    int t2 = dset->t2;
    int err;

    err = series_adjust_sample(x, &t1, &t2);

    if (!err) {
        err = filter_vector(x, y, t1, t2, A, C, y0, x0);
    }

    return err;
}

/**
 * filter_matrix:
 * @X: matrix of original data, r x c.
 * @A: vector for autoregressive polynomial.
 * @C: vector for moving average polynomial.
 * @y0: initial value of output series.
 * @x0: prior values of x.
 * @err: location to receive error code.
 *
 * Filters the columns of x according to y_t = C(L)/A(L) x_t.  If the
 * intended AR order is p, @A should be a vector of length p.  If the
 * intended MA order is q, @C should be vector of length (q+1), the
 * first entry giving the coefficient at lag 0.  However, if @C is
 * NULL this is taken to mean that the lag-0 MA coefficient is unity
 * (and all others are zero).
 *
 * Returns: r x c matrix of filtered values, or NULL on failure.
 */

gretl_matrix *filter_matrix (gretl_matrix *X, gretl_vector *A,
                             gretl_vector *C, double y0,
                             gretl_matrix *x0, int *err)
{
    int r = X->rows;
    int c = X->cols;
    gretl_matrix *Y = NULL;
    double *a = NULL, *b = NULL;
    int i, j;

    if (gretl_is_complex(X) ||
        gretl_is_complex(A) ||
        gretl_is_complex(C)) {
        *err = E_CMPLX;
        return NULL;
    }

    Y = gretl_matrix_alloc(r, c);
    a = calloc(r, sizeof *a);
    b = calloc(r, sizeof *b);

    if (Y == NULL || a == NULL || b == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    for (j=0; j<c; j++) {
        for (i=0; i<r; i++) {
            a[i] = gretl_matrix_get(X, i, j);
        }
        *err = filter_vector(a, b, 0, r-1, A, C, y0, x0);
        if (*err) {
            break;
        } else {
            for (i=0; i<r; i++) {
                gretl_matrix_set(Y, i, j, b[i]);
            }
        }
    }

    free(a);
    free(b);

    return Y;
}

static int series_goodobs (const double *x, int *t1, int *t2)
{
    int t, t1min = *t1, t2max = *t2;
    int T = 0;

    for (t=t1min; t<t2max; t++) {
        if (na(x[t])) t1min++;
        else break;
    }

    for (t=t2max; t>t1min; t--) {
        if (na(x[t])) t2max--;
        else break;
    }

    for (t=t1min; t<=t2max; t++) {
        if (!na(x[t])) {
            T++;
        }
    }

    *t1 = t1min;
    *t2 = t2max;

    return T;
}

/**
 * exponential_movavg_series:
 * @x: array of original data.
 * @y: array into which to write the result.
 * @dset: dataset information.
 * @d: coefficient on lagged @x.
 * @n: number of @x observations to average to give the
 * initial @y value.
 * @y0: optional initial value for EMA (or #NADBL).
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int exponential_movavg_series (const double *x, double *y,
                               const DATASET *dset,
                               double d, int n,
                               double y0)
{
    int t1 = dset->t1;
    int t2 = dset->t2;
    int t, T, err = 0;

    if (dataset_is_panel(dset)) {
        return E_PDWRONG;
    }

    if (na(y0) && n < 0) {
        return E_INVARG;
    }

    T = series_goodobs(x, &t1, &t2);

    if (na(y0) && T < n) {
        return E_TOOFEW;
    }

    if (na(y0) && n == 0) {
        /* signal to use full sample mean */
        n = T;
    }

    if (!na(y0)) {
        ; /* initialize on supplied value */
    } else if (n == 1) {
        /* initialize on first observation */
        y0 = x[t1];
    } else {
        /* initialize on mean of first n obs */
        y0 = 0.0;
        for (t=t1; t<t1+n; t++) {
            if (na(x[t])) {
                err = E_MISSDATA;
                break;
            }
            y0 += x[t];
        }
        if (!err) {
            y0 /= n;
        }
    }

    if (!err && na(y0)) {
        err = E_MISSDATA;
    }

    if (!err) {
        y[t1] = d * x[t1] + (1-d) * y0;
        for (t=t1+1; t<=t2; t++) {
            if (na(x[t]) || na(y[t-1])) {
                y[t] = NADBL;
            } else {
                y[t] = d * x[t] + (1-d) * y[t-1];
            }
        }
    }

    return err;
}

/**
 * movavg_series:
 * @x: array of original data.
 * @y: array into which to write the result.
 * @dset: data set information.
 * @k: number of terms in MA.
 * @center: if non-zero, produce centered MA.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int movavg_series (const double *x, double *y, const DATASET *dset,
                   int k, int center)
{
    int t1 = dset->t1;
    int t2 = dset->t2;
    int k1 = k-1, k2 = 0;
    int i, s, t, T;

    T = series_goodobs(x, &t1, &t2);

    T = t2 - t1 + 1;
    if (T < k) {
        return E_TOOFEW;
    }

    if (center) {
        k1 = k / 2;
        k2 = (k % 2 == 0)? (k1 - 1) : k1;
    }

    t1 += k1;
    t2 -= k2;

    for (t=t1; t<=t2; t++) {
        double xs, msum = 0.0;

        for (i=-k1; i<=k2; i++) {
            s = t + i;
            if (dset->structure == STACKED_TIME_SERIES) {
                if (s / dset->pd != t / dset->pd) {
                    s = -1;
                }
            }

            if (s >= 0) {
                xs = x[s];
            } else {
                xs = NADBL;
            }

            if (na(xs)) {
                msum = NADBL;
                break;
            } else {
                msum += x[s];
            }
        }

        if (!na(msum)) {
            y[t] = (k > 0)? (msum / k) : msum;
        }
    }

    if (center && k % 2 == 0) {
        /* centered, with even number of terms */
        for (t=t1; t<t2; t++) {
            if (na(y[t]) || na(y[t+1])) {
                y[t] = NADBL;
            } else {
                y[t] = (y[t] + y[t+1]) / 2.0;
            }
        }
        y[t2] = NADBL;
    }

    return 0;
}

int seasonally_adjust_series (const double *x, double *y,
                              const char *vname, DATASET *dset,
                              int tramo, gretl_bundle *b,
                              PRN *prn)
{
    int (*adjust_series) (const double *, double *,
                          const char *, const DATASET *,
                          int, gretl_bundle *, PRN *);
    int t1 = dset->t1;
    int t2 = dset->t2;
    int T, err = 0;

    if (!quarterly_or_monthly(dset)) {
        gretl_errmsg_set(_("Input must be a monthly or quarterly time series"));
        return E_PDWRONG;
    }

    series_adjust_sample(x, &t1, &t2);
    T = t2 - t1 + 1;

    if (T < dset->pd * 3) {
        return E_TOOFEW;
    } else if (tramo && T > 600) {
        gretl_errmsg_set(_("TRAMO can't handle more than 600 observations.\n"
                           "Please select a smaller sample."));
        return E_EXTERNAL;
    } else if (!tramo) {
        int pdmax = get_x13as_maxpd();

        if (T > 50 * pdmax) {
            gretl_errmsg_sprintf(_("X-13ARIMA can't handle more than %d observations.\n"
                                   "Please select a smaller sample."), 50 * pdmax);
            return E_EXTERNAL;
        }
    }

    gretl_error_clear();

    adjust_series = get_plugin_function("adjust_series");

    if (adjust_series == NULL) {
        err = E_FOPEN;
    } else {
        int save_t1 = dset->t1;
        int save_t2 = dset->t2;

        dset->t1 = t1;
        dset->t2 = t2;

        err = (*adjust_series) (x, y, vname, dset, tramo, b, prn);

        dset->t1 = save_t1;
        dset->t2 = save_t2;
    }

    return err;
}

int tramo_linearize_series (const double *x, double *y,
                            DATASET *dset)
{
    int (*linearize_series) (const double *, double *,
                             const DATASET *);
    int t1 = dset->t1;
    int t2 = dset->t2;
    int T, err = 0;

    series_adjust_sample(x, &t1, &t2);
    T = t2 - t1 + 1;

    if (T < 8) {
        return E_TOOFEW;
    } else if (T > 600) {
        gretl_errmsg_set(_("TRAMO can't handle more than 600 observations.\n"
                           "Please select a smaller sample."));
        return E_EXTERNAL;
    }

    gretl_error_clear();

    linearize_series = get_plugin_function("linearize_series");

    if (linearize_series == NULL) {
        err = E_FOPEN;
    } else {
        int save_t1 = dset->t1;
        int save_t2 = dset->t2;

        dset->t1 = t1;
        dset->t2 = t2;

        err = (*linearize_series) (x, y, dset);

        dset->t1 = save_t1;
        dset->t2 = save_t2;
    }

    return err;
}

#define panel_obs_ok(x,t,m) ((m == NULL || m[t] != 0) && !na(x[t]))
#define panel_obs_masked(m,t) (m != NULL && m[t] == 0)

#define PXSUM_SKIP_NA 1

/* Run a check for time-invariance of a panel series. If it isn't
   time-invariant that should become apparent quickly, but if it is
   we get a substantial shortcut for standard deviation.
*/

static int time_invariant (const double *x,
                           int u1, int u2, int T,
                           const double *mask)
{
    double xmin, xmax;
    int i, t, s;
    int ret = 1;

    for (i=u1; i<=u2 && ret; i++) {
        xmin = xmax = NADBL;
        for (t=0; t<T && ret; t++) {
            s = i*T + t;
            if (panel_obs_ok(x, s, mask)) {
                if (na(xmax)) {
                    xmin = xmax = x[s];
                } else if (x[s] < xmin) {
                    xmin = x[s];
                } else if (x[s] > xmax) {
                    xmax = x[s];
                }
            }
            if (xmax > xmin) {
                ret = 0;
            }
        }
        if (ret && na(xmax)) {
            ret = 0;
        }
    }

    return ret;
}

/**
 * panel_statistic:
 * @x: source data.
 * @y: target into which to write.
 * @dset: data set information.
 * @k: code representing the desired statistic: F_PNOBS,
 * F_PMIN, F_PMAX, F_PSUM, F_PMEAN, F_PXSUM, F_PSD, F_PXNOBS
 * or F_PXSUM.
 * @mask: either NULL or a series with 0s for observations
 * to be excluded from the calculations, non-zero values
 * at other observations.
 *
 * Given the data in @x, constructs in @y a series containing
 * a panel-data statistic.
 *
 * Returns: 0 on success, non-zero on error.
 */

int panel_statistic (const double *x, double *y, const DATASET *dset,
                     int k, const double *mask)
{
    int T, Ti;
    int i, s, t, u1, u2;
    int err = 0;

    if (!dataset_is_panel(dset)) {
        return E_DATA;
    }

    T = dset->pd;
    u1 = dset->t1 / T;
    u2 = dset->t2 / T;

    if (k == F_PSD && time_invariant(x, u1, u2, T, mask)) {
        for (i=u1; i<=u2; i++) {
            for (t=0; t<T; t++) {
                y[i*T + t] = 0.0;
            }
        }
        return 0;
    }

    if (k == F_PNOBS) {
	/* the number of valid observations */
        for (i=u1; i<=u2; i++) {
            Ti = 0;
            for (t=0; t<T; t++) {
                if (panel_obs_ok(x, i*T + t, mask)) {
                    Ti++;
                }
            }
            for (t=0; t<T; t++) {
                y[i*T + t] = Ti;
            }
        }
    } else if (k == F_PMIN || k == F_PMAX || k == F_PSUM) {
        double yi;

        for (i=u1; i<=u2; i++) {
            yi = NADBL;
            for (t=0; t<T; t++) {
                s = i*T + t;
                if (panel_obs_ok(x, s, mask)) {
                    if (na(yi)) {
			/* the first valid observation */
                        yi = x[s];
                    } else if (k == F_PMIN && x[s] < yi) {
                        yi = x[s];
                    } else if (k == F_PMAX && x[s] > yi) {
                        yi = x[s];
                    } else if (k == F_PSUM) {
                        yi += x[s];
                    }
                }
            }
            for (t=0; t<T; t++) {
                y[i*T + t] = yi;
            }
        }
    } else if (k == F_PMEAN || k == F_PSD) {
	/* in both of these cases we need the mean */
	double xbar;

	for (i=u1; i<=u2; i++) {
	    xbar = NADBL;
	    Ti = 0;
	    for (t=0; t<T; t++) {
		s = i*T + t;
		if (panel_obs_ok(x, s, mask)) {
		    if (na(xbar)) {
			/* the first valid observation */
			xbar = x[s];
		    } else {
			/* cumulate */
			xbar += x[s];
		    }
		    Ti++;
		}
	    }
	    if (!na(xbar)) {
		xbar /= Ti;
	    }
	    for (t=0; t<T; t++) {
		y[i*T + t] = xbar;
	    }
	}
        if (k == F_PSD) {
	    /* now we need the deviations from the mean */
	    double xbar, dev, ssx, sd;

            for (i=u1; i<=u2; i++) {
                xbar = y[i*T]; /* the mean calculated above */
		if (na(xbar)) {
		    /* no calculation needed */
		    for (t=1; t<T; t++) {
			y[i*T + t] = NADBL;
		    }
		    continue;
		}
		ssx = 0;
                Ti = 0;
		for (t=0; t<T; t++) {
		    s = i*T + t;
		    if (panel_obs_ok(x, s, mask)) {
			dev = x[s] - xbar;
			ssx += dev * dev;
			Ti++;
		    }
		}
                if (Ti == 1) {
                    sd = 0.0;
                } else {
                    sd = sqrt(ssx / (Ti-1));
                }
                for (t=0; t<T; t++) {
                    y[i*T + t] = sd;
                }
            }
        }
    } else if (k == F_PXSUM) {
        /* the sum of cross-sectional values for each period */
        double yt;
        int nt;

        for (t=0; t<T; t++) {
            yt = 0.0;
            nt = 0;
            for (i=u1; i<=u2; i++) {
                s = i*T + t;
                if (panel_obs_masked(mask, s)) {
                    continue;
                } else if (na(x[s])) {
#if PXSUM_SKIP_NA
                    continue;
#else
                    yt = NADBL;
                    break;
#endif
                } else {
                    yt += x[s];
                    nt++;
                }
            }
            if (nt == 0) {
                yt = NADBL;
            }
            for (i=u1; i<=u2; i++) {
                y[i*T + t] = yt;
            }
        }
    } else if (k == F_PXNOBS) {
        /* number of valid x-sectional obs in each period */
        int nt;

        for (t=0; t<T; t++) {
            nt = 0;
            for (i=u1; i<=u2; i++) {
                s = i*T + t;
                if (panel_obs_masked(mask, s)) {
                    continue;
                } else if (!na(x[s])) {
                    nt++;
                }
            }
            for (i=u1; i<=u2; i++) {
                y[i*T + t] = nt;
            }
        }
    } else {
        /* unsupported option */
        err = E_DATA;
    }

    return err;
}

/**
 * panel_shrink:
 * @x: panel-data source series.
 * @noskip: keep NAs in output.
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * By default, constructs a column vector holding the first
 * non-missing observation of @x for each panel unit within
 * the current sample range (hence skipping any units that
 * have no valid observations). However, if @noskip is non-zero,
 * the vector contains an NA for units with all missing
 * values.
 *
 * Returns: a new column vector, or NULL on error.
 */

gretl_matrix *panel_shrink (const double *x, int noskip,
                            const DATASET *dset, int *err)
{
    gretl_matrix *m = NULL;
    int n, T;

    if (!dataset_is_panel(dset)) {
        *err = E_DATA;
        return NULL;
    }

    n = panel_sample_size(dset);
    T = dset->pd;
    m = gretl_column_vector_alloc(n);

    if (m == NULL) {
        *err = E_ALLOC;
    } else {
        int i1 = dset->t1 / T;
        int i2 = dset->t2 / T;
        int s = dset->t1;
        int i, t, k = 0;
        int gotit;

        for (i=i1; i<=i2; i++) {
            /* loop across units */
            gotit = 0;
            for (t=0; t<T; t++) {
                /* loop inside units */
                if (!na(x[s+t])) {
                    gotit = 1;
                    m->val[k++] = x[s+t];
                    break;
                }
            }
            if (!gotit && noskip) {
                m->val[k++] = NADBL;
            }
            /* skip to the next unit */
            s += T;
        }

        if (k < n) {
            m->rows = k;
        }
    }

    return m;
}

/**
 * panel_expand:
 * @x: source vector.
 * @y target array.
 * @dset: pointer to panel dataset.
 * @opt: OPT_X to repeat across individuals instead
 * of across time.
 *
 * Constructs a panel series by repeating the values in @x.
 * By default expansion is across time, and the source vector
 * must have as many elements as there are individuals in the
 * current sample range. But if OPT_X is given, expansion is
 * across individuals and the length of @x must equal the
 * time-series length of the dataset.
 *
 * Returns: zero on success, non-zero code on error.
 */

int panel_expand (const gretl_matrix *x, double *y,
                  gretlopt opt, const DATASET *dset)
{
    int xlen, N, T;
    int xsect = 0;
    int i, t, s;

    if (x == NULL || !dataset_is_panel(dset)) {
        return E_DATA;
    }

    T = dset->pd;
    N = sample_size(dset) / T;
    xlen = gretl_vector_get_length(x);

    if (opt & OPT_X) {
        if (xlen != T) {
            return E_INVARG;
        } else {
            xsect = 1;
        }
    } else if (xlen != N) {
        return E_INVARG;
    }

    s = dset->t1;
    for (i=0; i<N; i++) {
        for (t=0; t<T; t++) {
            y[s++] = xsect? x->val[t] : x->val[i];
        }
    }

    return 0;
}

static double default_hp_lambda (const DATASET *dset)
{
    return 100 * dset->pd * dset->pd;
}

/**
 * hp_filter:
 * @x: array of original data.
 * @hp: array in which filtered series is computed.
 * @dset: pointer to dataset.
 * @lambda: smoothing parameter (or #NADBL to use the default
 * value).
 * @opt: if %OPT_T, return the trend rather than the cycle.
 *
 * Calculates the "cycle" component of the time series in
 * array @x, using the Hodrick-Prescott filter.  Adapted from the
 * original FORTRAN code by E. Prescott.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int hp_filter (const double *x, double *hp, const DATASET *dset,
               double lambda, gretlopt opt)
{
    int t1 = dset->t1, t2 = dset->t2;
    double v00 = 1.0, v11 = 1.0, v01 = 0.0;
    double det, tmp0, tmp1;
    double e0, e1, b00, b01, b11;
    double **V = NULL;
    double m[2], tmp[2];
    int i, t, T, tb;
    int err = 0;

    for (t=t1; t<=t2; t++) {
        hp[t] = NADBL;
    }

    err = series_adjust_sample(x, &t1, &t2);
    if (err) {
        goto bailout;
    }

    T = t2 - t1 + 1;
    if (T < 4) {
        err = E_TOOFEW;
        goto bailout;
    }

    if (na(lambda)) {
        lambda = default_hp_lambda(dset);
    }

    V = doubles_array_new(4, T);
    if (V == NULL) {
        return E_ALLOC;
    }

    /* adjust starting points */
    x += t1;
    hp += t1;

    /* covariance matrices for each obs */

    for (t=2; t<T; t++) {
        tmp0 = v00;
        tmp1 = v01;
        v00 = 1.0 / lambda + 4.0 * (tmp0 - tmp1) + v11;
        v01 = 2.0 * tmp0 - tmp1;
        v11 = tmp0;

        det = v00 * v11 - v01 * v01;
        V[0][t] =  v11 / det;
        V[1][t] = -v01 / det;
        V[2][t] =  v00 / det;

        tmp0 = v00 + 1.0;
        tmp1 = v00;
        v00 -= v00 * v00 / tmp0;
        v11 -= v01 * v01 / tmp0;
        v01 -= (tmp1 / tmp0) * v01;
    }

    m[0] = x[0];
    m[1] = x[1];

    /* forward pass */
    for (t=2; t<T; t++) {
        tmp[0] = m[1];
        m[1] = 2.0 * m[1] - m[0];
        m[0] = tmp[0];

        V[3][t-1] = V[0][t] * m[1] + V[1][t] * m[0];
        hp[t-1]   = V[1][t] * m[1] + V[2][t] * m[0];

        det = V[0][t] * V[2][t] - V[1][t] * V[1][t];
        v00 =  V[2][t] / det;
        v01 = -V[1][t] / det;

        tmp[1] = (x[t] - m[1]) / (v00 + 1.0);
        m[1] += v00 * tmp[1];
        m[0] += v01 * tmp[1];
    }

    V[3][T-2] = m[0];
    V[3][T-1] = m[1];
    m[0] = x[T-2];
    m[1] = x[T-1];

    /* backward pass */
    for (t=T-3; t>=0; t--) {
        t1 = t+1;
        tb = T - t - 1;

        tmp[0] = m[0];
        m[0] = 2.0 * m[0] - m[1];
        m[1] = tmp[0];

        if (t > 1) {
            /* combine info for y < i with info for y > i */
            e0 = V[2][tb] * m[1] + V[1][tb] * m[0] + V[3][t];
            e1 = V[1][tb] * m[1] + V[0][tb] * m[0] + hp[t];
            b00 = V[2][tb] + V[0][t1];
            b01 = V[1][tb] + V[1][t1];
            b11 = V[0][tb] + V[2][t1];

            det = b00 * b11 - b01 * b01;
            V[3][t] = (b00 * e1 - b01 * e0) / det;
        }

        det = V[0][tb] * V[2][tb] - V[1][tb] * V[1][tb];
        v00 =  V[2][tb] / det;
        v01 = -V[1][tb] / det;

        tmp[1] = (x[t] - m[0]) / (v00 + 1.0);
        m[1] += v01 * tmp[1];
        m[0] += v00 * tmp[1];
    }

    V[3][0] = m[0];
    V[3][1] = m[1];

    if (opt & OPT_T) {
        for (t=0; t<T; t++) {
            hp[t] = V[3][t];
        }
    } else {
        for (t=0; t<T; t++) {
            hp[t] = x[t] - V[3][t];
        }
    }

 bailout:

    if (V != NULL) {
        for (i=0; i<4; i++) {
            free(V[i]);
        }
        free(V);
    }

    return err;
}

/**
 * oshp_filter:
 * @x: array of original data.
 * @hp: array in which filtered series is computed.
 * @dset: pointer to dataset.
 * @lambda: smoothing parameter (or #NADBL to use the default
 * value).
 * @opt: if %OPT_T, return the trend rather than the cycle.
 *
 * Calculates the "cycle" component of the time series in
 * array @x, using the a one-sided Hodrick-Prescott filter.
 * The implementation uses the Kalman filter.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int oshp_filter (const double *x, double *hp, const DATASET *dset,
                 double lambda, gretlopt opt)
{
    int t1 = dset->t1, t2 = dset->t2;
    gretl_matrix *M[4] = {NULL};
    gretl_matrix *Lambda = NULL;
    gretl_matrix *a0 = NULL;
    gretl_matrix *mu = NULL;
    gretl_bundle *b = NULL;
    int copy[4] = {0};
    int T, t, err, kerr;
    double mt, sqrt_lam;

    for (t=t1; t<=t2; t++) {
        hp[t] = NADBL;
    }

    err = series_adjust_sample(x, &t1, &t2);
    if (err) {
        return err;
    }

    T = t2 - t1 + 1;
    if (T < 4) {
        return E_TOOFEW;
    }

    if (na(lambda)) {
        lambda = default_hp_lambda(dset);
    }
    sqrt_lam = sqrt(lambda);

    /* adjust starting points */
    x += t1;
    hp += t1;

    /* begin Kalman filter setup */

    /* obsy */
    M[0] = gretl_matrix_alloc(T+1, 1);
    for (t=0; t<T; t++) {
        gretl_matrix_set(M[0], t, 0, x[t]);
    }
    /* add a dummy trailing observation */
    gretl_matrix_set(M[0], T, 0, 0.0);

    /* obsymat */
    M[1] = gretl_zero_matrix_new(2, 1);
    gretl_matrix_set(M[1], 0, 0, 1);

    /* statemat */
    M[2] = gretl_zero_matrix_new(2, 2);
    gretl_matrix_set(M[2], 0, 0, 2);
    gretl_matrix_set(M[2], 0, 1, -1);
    gretl_matrix_set(M[2], 1, 0, 1);

    /* statevar */
    M[3] = gretl_zero_matrix_new(2, 2);
    gretl_matrix_set(M[3], 0, 0, 1/sqrt_lam);

    b = kalman_bundle_new(M, copy, 4, 0, &err);
    if (err) {
        goto bailout;
    }

    /* obsvar */
    Lambda = gretl_matrix_from_scalar(sqrt_lam);
    err = gretl_bundle_donate_data(b, "obsvar", Lambda,
                                   GRETL_TYPE_MATRIX, 0);

    /* diffuse prior */
    err = gretl_bundle_set_scalar(b, "diffuse", 1);

    if (err) {
        goto bailout;
    }

    /* inistate */
    a0 = gretl_matrix_alloc(2, 1);
    gretl_matrix_set(a0, 0, 0, 2*x[0] - x[1]);
    gretl_matrix_set(a0, 1, 0, 3*x[0] - 2*x[1]);
    err = gretl_bundle_donate_data(b, "inistate", a0,
                                   GRETL_TYPE_MATRIX, 0);
    if (err) {
        goto bailout;
    }

#if DEBUG
    gretl_matrix_print(M[0], "obsy");
    gretl_matrix_print(M[1], "obsymat");
    gretl_matrix_print(M[2], "statemat");
    gretl_matrix_print(M[3], "statevar");
    gretl_matrix_print(a0, "inistate");
#endif

    kerr = kalman_bundle_filter(b, NULL, &err);
    if (err || kerr) {
        goto bailout;
    }

    mu = gretl_bundle_get_matrix(b, "state", &err);
    if (err) {
        goto bailout;
    }

#if DEBUG
    printf("kerr = %d, err = %d\n", kerr, err);
    gretl_matrix_print(mu, "state");
#endif

    /* take the "lead" of the second element of the
       state estimate */
    if (opt & OPT_T) {
        for (t=0; t<T; t++) {
            mt = gretl_matrix_get(mu, t+1, 1);
            hp[t] = mt;
        }
    } else {
        for (t=0; t<T; t++) {
            mt = gretl_matrix_get(mu, t+1, 1);
            hp[t] = x[t] - mt;
        }
    }

 bailout:

    /* since all the matrices we've allocated above belong
       to the Kalman bundle, the following will suffice to
       free everything on exit
    */
    gretl_bundle_destroy(b);

    return err;
}

/**
 * bkbp_filter:
 * @x: array of original data.
 * @bk: array into which to write the filtered series.
 * @dset: data set information.
 * @bkl: lower frequency bound (or 0 for automatic).
 * @bku: upper frequency bound (or 0 for automatic).
 * @k: approximation order (or 0 for automatic).
 *
 * Calculates the Baxter and King bandpass filter. If bku exceeds the
 * number of available observations, then the "low-pass" version will
 * be computed (weights sum to 1). Otherwise, the ordinary bandpass
 * filter will be applied (weights sum to 0).
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int bkbp_filter (const double *x, double *bk, const DATASET *dset,
                 int bkl, int bku, int k)
{
    int t1 = dset->t1, t2 = dset->t2;
    double omubar, omlbar;
    double avg_a;
    double *a;
    int i, t, err = 0;
    int lowpass = 0;

    if (bkl <= 0 || bku <= 0) {
        /* get user settings if available (or the defaults) */
        get_bkbp_periods(dset, &bkl, &bku);
    }

    if (k <= 0) {
        k = get_bkbp_k(dset);
    }

#if BK_DEBUG
    fprintf(stderr, "lower limit = %d, upper limit = %d, \n",
            bkl, bku);
#endif

    if (bkl >= bku) {
        gretl_errmsg_set(_("Error in Baxter-King frequencies"));
        return 1;
    }

    err = series_adjust_sample(x, &t1, &t2);
    if (err) {
        return err;
    }

    if (2 * k >= t2 - t1 + 1) {
        gretl_errmsg_set(_("Insufficient observations"));
        return E_DATA;
    }

    if (bku >= t2 - t1 + 1) {
        lowpass = 1;
    }

    a = malloc((k + 1) * sizeof *a);
    if (a == NULL) {
        return E_ALLOC;
    }

    omubar = M_2PI / bkl;
    omlbar = lowpass? 0 : M_2PI / bku;

    /* first we compute the coefficients */

    avg_a = a[0] = (omubar - omlbar) / M_PI;

    if (lowpass) {
        avg_a -= 1.0;
    }

    for (i=1; i<=k; i++) {
        a[i] = (sin(i * omubar) - sin(i * omlbar)) / (i * M_PI);
        avg_a += 2 * a[i];
    }

    avg_a /= (2 * k + 1);

    for (i=0; i<=k; i++) {
        a[i] -= avg_a;
#if BK_DEBUG
        fprintf(stderr, "a[%d] = %#9.6g\n", i, a[i]);
#endif
    }

    /* now we filter the series, skipping the first
       and last k observations */

    for (t=0; t<dset->n; t++) {
        if (t < t1 + k || t > t2 - k) {
            bk[t] = NADBL;
        } else {
            bk[t] = a[0] * x[t];
            for (i=1; i<=k; i++) {
                bk[t] += a[i] * (x[t-i] + x[t+i]);
            }
        }
    }

    free(a);

    return err;
}

/* following: material relating to the Butterworth filter */

#define Min(x,y) (((x) > (y))? (y) : (x))
#define NEAR_ZERO 1e-35

static double safe_pow (double x, int n)
{
    double lp;

    x = fabs(x);

    if (x < NEAR_ZERO) {
        return 0.0;
    } else {
        lp = n * log(x);
        if (lp < -80) {
            return 0.0;
        } else if (lp > 80) {
            return exp(80.0);
        } else {
            return exp(lp);
        }
    }
}

static double cotan (double theta)
{
    double s = sin(theta);

    if (fabs(s) < NEAR_ZERO) {
        s = copysign(NEAR_ZERO, s);
    }

    return cos(theta) / s;
}

/**
 * hp_gain:
 * @lambda: H-P parameter.
 * @hipass: 1 for high-pass filter, 0 for low-pass.
 *
 * Returns: a matrix holding omega values from 0 to \pi in
 * column 0, and the corresponding filter gain in column 1.
 */

gretl_matrix *hp_gain (double lambda, int hipass)
{
    gretl_matrix *G;
    int i, width = 180;
    double inc = M_PI / width;
    double x, gain, omega = 0.0;

    G = gretl_matrix_alloc(width + 1, 2);
    if (G == NULL) {
        return NULL;
    }

    for (i=0; i<=width; i++) {
        x = 2 * sin(omega / 2);
        gain = 1 / (1 + lambda * pow(x, 4));
        if (hipass) {
            gain = 1.0 - gain;
        }
        gretl_matrix_set(G, i, 0, omega);
        gretl_matrix_set(G, i, 1, gain);
        omega += inc;
    }

    return G;
}

/**
 * butterworth_gain:
 * @n: order of the filter.
 * @cutoff: angular cutoff in radians.
 * @hipass: 1 for high-pass filter, 0 for low-pass.
 *
 * Returns: a matrix holding omega values from 0 to \pi in
 * column 0, and the corresponding filter gain in column 1.
 */

gretl_matrix *butterworth_gain (int n, double cutoff, int hipass)
{
    gretl_matrix *G;
    int i, width = 180;
    double inc = M_PI / width;
    double x, gain, omega = 0.0;

    G = gretl_matrix_alloc(width + 1, 2);
    if (G == NULL) {
        return NULL;
    }

    for (i=0; i<=width; i++) {
        if (hipass) {
            x = cotan(omega / 2) * cotan((M_PI - cutoff) / 2);
        } else {
            x = tan(omega / 2) * cotan(cutoff / 2);
        }
        gain = 1 / (1 + pow(x, 2 * n));
        gretl_matrix_set(G, i, 0, omega);
        gretl_matrix_set(G, i, 1, gain);
        omega += inc;
    }

    return G;
}

/* Toeplitz methods:

   1 = Stephen Pollock's symmetric Toeplitz solver
   2 = gretl's netlib-based general Toeplitz solver
   3 = build full Toeplitz matrix and apply SVD

*/

#define TOEPLITZ_METHOD 1

#if (TOEPLITZ_METHOD == 3)

/* form the complete Toeplitz matrix represented in compressed
   form by the coefficients in g
*/

static gretl_matrix *toeplize (const double *g, int T, int k)
{
    gretl_matrix *X;
    double x;
    int i, j;

    X = gretl_zero_matrix_new(T, T);
    if (X == NULL) {
        return NULL;
    }

    for (i=0; i<T; i++) {
        gretl_matrix_set(X, i, i, g[0]);
    }

    for (i=1; i<k; i++) {
        x = g[i];
        for (j=i; j<T; j++) {
            gretl_matrix_set(X, j, j-i, x);
            gretl_matrix_set(X, j-i, j, x);
        }
    }

    return X;
}

#include "matrix_extra.h"

static int toeplitz_solve (const double *g, double *dy, int T, int q)
{
    gretl_matrix *y, *X;
    int err = 0;

    fprintf(stderr, "svd_toeplitz_solve...\n");

    X = toeplize(g, T, q + 1);
    if (X == NULL) {
        return E_ALLOC;
    }

    err = gretl_SVD_invert_matrix(X);

    if (!err) {
        y = gretl_vector_from_array(dy, T, GRETL_MOD_NONE);
        if (y == NULL) {
            err = E_ALLOC;
        } else {
            gretl_matrix s;

            s.t1 = s.t2 = 0;
            s.rows = T;
            s.cols = 1;
            s.val = dy;

            err = gretl_matrix_multiply(X, y, &s);
            gretl_matrix_free(y);
        }
    }

    gretl_matrix_free(X);

    return err;
}

#elif (TOEPLITZ_METHOD == 2)

static int toeplitz_solve (double *g, double *y, int T, int q)
{
    gretl_vector *mg = NULL;
    gretl_vector *my = NULL;
    gretl_vector *mx = NULL;
    int i, err = 0;

    fprintf(stderr, "netlib toeplitz solve...\n");

    mg = gretl_vector_alloc(T);
    my = gretl_vector_alloc(T);

    if (mg == NULL || my == NULL) {
        return E_ALLOC;
    }

    for (i=0; i<T; i++) {
        mg->val[i] = (i <= q) ? g[i] : 0;
        my->val[i] = y[i];
    }

    mx = gretl_toeplitz_solve(mg, mg, my, &err);

    if (err) {
        fprintf(stderr, "symm_toeplitz: err = %d\n", err);
        for (i=0; i<T; i++) {
            y[i] = NADBL;
        }
    } else {
        for (i=0; i<T; i++) {
            y[i] = mx->val[i];
        }
    }

    gretl_vector_free(mg);
    gretl_vector_free(my);
    gretl_vector_free(mx);

    return err;
}

#else /* TOEPLITZ_METHOD == 1 */

/* This function uses a Cholesky decomposition to find the solution of
   the equation Gx = y, where G is a symmetric Toeplitz matrix of order
   T with q supra-diagonal and q sub-diagonal bands. The coefficients of
   G are contained in g. The RHS vector y contains the elements of
   the solution vector x on completion.
*/

static int toeplitz_solve (double *g, double *y, int T, int q)
{
    double **mu;
    int t, j, k, jmax;

    mu = doubles_array_new(q+1, T);
    if (mu == NULL) {
        return E_ALLOC;
    }

    /* initialize */
    for (t=0; t<q; t++) {
        for (j=t+1; j<=q; j++) {
            mu[j][t] = 0.0;
        }
    }

    /* factorize */
    for (t=0; t<T; t++) {
        for (k=Min(q, t); k>=0; k--) {
            mu[k][t] = g[k];
            jmax = q - k;
            for (j=1; j<=jmax && t-k-j >= 0; j++) {
                mu[k][t] -= mu[j][t-k] * mu[j+k][t] * mu[0][t-k-j];
            }
            if (k > 0) {
                mu[k][t] /= mu[0][t-k];
            }
        }
    }

    /* forward solve */
    for (t=0; t<T; t++) {
        jmax = Min(t, q);
        for (j=1; j<=jmax; j++) {
            y[t] -= mu[j][t] * y[t-j];
        }
    }

    /* divide by the diagonal */
    for (t=0; t<T; t++) {
        y[t] /= mu[0][t];
    }

    /* backsolve */
    for (t=T-1; t>=0; t--) {
        jmax = Min(q, T - 1 - t);
        for (j=1; j<=jmax; j++) {
            y[t] -= mu[j][t+j] * y[t+j];
        }
    }

    doubles_array_free(mu, q+1);

    return 0;
}

#endif /* alternate TOEPLITZ_METHODs */


/* Premultiply vector y by a symmetric banded Toeplitz matrix
   Gamma with n nonzero sub-diagonal bands and n nonzero
   supra-diagonal bands. The elements of Gamma are contained
   in g; tmp is used as workspace.
*/

static int GammaY (double *g, double *y, double *tmp,
                   int T, int n)
{
    double lx, rx;
    int t, j;

    for (j=0; j<=n; j++) {
        tmp[j] = 0.0;
    }

    for (t=0; t<T; t++) {
        for (j=n; j>0; j--) {
            tmp[j] = tmp[j-1];
        }
        tmp[0] = g[0] * y[t];
        for (j=1; j<=n; j++) {
            lx = (t - j < 0) ? 0 : y[t-j];
            rx = (t + j >= T) ? 0 : y[t+j];
            tmp[0] += g[j] * (lx + rx);
        }
        if (t >= n) {
            y[t-n] = tmp[n];
        }
    }

    for (j=0; j<n; j++) {
        y[T-j-1] = tmp[j];
    }

    return 0;
}

#undef Min
#undef NEAR_ZERO

/* Form the autocovariances of an MA process of order q. */

static void form_gamma (double *g, double *mu, int q)
{
    int j, k;

    for (j=0; j<=q; j++) {
        g[j] = 0.0;
        for (k=0; k<=q-j; k++) {
            g[j] += mu[k] * mu[k+j];
        }
    }
}

/* Find the coefficients of the nth power of the summation
   operator (if sign = 1) or of the difference operator (if
   sign = -1).
*/

static void sum_or_diff (double *theta, int n, int sign)
{
    int j, q;

    theta[0] = 1.0;

    for (q=1; q<=n; q++) {
        theta[q] = 0.0;
        for (j=q; j>0; j--) {
            theta[j] += sign * theta[j-1];
        }
    }
}

/* g is the target, mu is used as workspace for the
   summation coefficients */

static void form_mvec (double *g, double *mu, int n)
{
    sum_or_diff(mu, n, 1);
    form_gamma(g, mu, n);
}

/* g is the target, mu is used as workspace for the
   differencing coefficients
   svec = Q'SQ where Q is the 2nd-diff operator
*/

static void form_svec (double *g, double *mu, int n)
{
    sum_or_diff(mu, n, -1);
    form_gamma(g, mu, n);
}

/* g is the target, mu and tmp are used as workspace */

static void form_wvec (double *g, double *mu, double *tmp,
                       int n, double lam1, double lam2)
{
    int i;

    form_svec(tmp, mu, n);
    for (i=0; i<=n; i++) {
        g[i] = lam1 * tmp[i];
    }

    form_mvec(tmp, mu, n);
    for (i=0; i<=n; i++) {
        g[i] += lam2 * tmp[i];
    }
}

/* Find the second differences of the elements of a vector.
   Notice that the first two elements of the vector are lost
   in the process. However, we index the leading element of
   the differenced vector by t = 0.
*/

static void QprimeY (double *y, int T)
{
    int t;

    for (t=0; t<T-2; t++) {
        y[t] += y[t+2] - 2 * y[t+1];
    }
}

/* Multiply the vector y by matrix Q of order T x T-2,
   where Q' is the matrix which finds the second
   differences of a vector.
*/

static void form_Qy (double *y, int T)
{
    double tmp, lag1 = 0.0, lag2 = 0.0;
    int t;

    for (t=0; t<T-2; t++) {
        tmp = y[t];
        y[t] += lag2 - 2 * lag1;
        lag2 = lag1;
        lag1 = tmp;
    }

    y[T-2] = lag2 - 2 * lag1;
    y[T-1] = lag1;
}

#if HAVE_GMP

static int mp_butterworth (const double *x, double *bw, int T,
                           int order, double cutoff)
{
    int (*mpfun) (const double *, double *, int, int, double);

    mpfun = get_plugin_function("mp_bw_filter");

    if (mpfun == NULL) {
        return E_FOPEN;
    }

    return (*mpfun) (x, bw, T, order, cutoff);
}

#endif

#if 0

/* note: cut should be in degrees */

static int bw_numeric_check (double cut, int n)
{
    const double b0 = -8.56644923648263;
    const double b1 = -0.66992007228155;
    const double b2 =  2.82851981714539;
    double x;

    /* let's be conservative here */
    x = b0 + b1 * (cut - 1.0) + b2 * (n + 1);
    x = 1/(1+exp(-x));

    return (x > 0.479717591568)
}

#endif

/* maximum modulus among the poles of the filter:
   too close to 1.0 is a problem
*/

static double bw_max_mod (double cut, int n)
{
    double tc = tan(cut / 2);
    double theta = M_PI / (2 * n);
    double delta = 1 + tc * (2 * sin(theta) + tc);
    double x = (1 - tc * tc) / delta;
    double y = 2 * tc * cos(theta) / delta;

    return sqrt(x * x + y * y);
}

#define MAX_LAM 10000

/* Calculate the Butterworth lambda based on cutoff and
   order. Return non-zero if it seems that lambda is too
   extreme. FIXME: conditioning on lambda alone is not an
   adequate approach.
*/

static int
set_bw_lambda (double cutoff, int n, double *lam1, double *lam2)
{
    int ret = 0;

    *lam1 = 1 / tan(cutoff / 2);
    *lam1 = safe_pow(*lam1, n * 2);

    fprintf(stderr, "for cutoff %g, order %d: lambda=%g, maxmod=%g\n",
            cutoff, n, *lam1, bw_max_mod(cutoff, n));

    if (*lam1 > 1.0e18) {
        /* can't cope, even with multiple precision? */
        ret = 2;
    } else if (*lam1 > 1.0e6) {
        /* may be OK with multiple precision? */
        ret = 1;
    } else {
        /* OK with regular double-precision? */
        if (*lam1 > MAX_LAM) {
            /* note: AC thinks the following doesn't help */
            *lam1 = sqrt(*lam1);
            *lam2 = 1 / *lam1;
        }
    }

    return ret;
}

/**
 * butterworth_filter:
 * @x: array of original data.
 * @bw: array into which to write the filtered series.
 * @dset: data set information.
 * @n: desired lag order.
 * @cutoff: desired angular cutoff in degrees (0, 180).
 *
 * Calculates the Butterworth filter. The code that does this
 * is based on D.S.G. Pollock's IDEOLOG -- see
 * http://www.le.ac.uk/users/dsgp1/
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int butterworth_filter (const double *x, double *bw, const DATASET *dset,
                        int n, double cutoff)
{
    double *g, *ds, *tmp, *y;
    double lam1, lam2 = 1.0;
    int t1 = dset->t1, t2 = dset->t2;
    int T, t, m;
    int bad_lambda = 0;
    int err = 0;

    if (!dataset_is_time_series(dset)) {
        gretl_errmsg_set(_("This function requires time-series data"));
        return E_DATA;
    }

    err = series_adjust_sample(x, &t1, &t2);
    if (err) {
        return err;
    }

    if (2 * n >= t2 - t1) {
        gretl_errmsg_set(_("Insufficient observations"));
        return E_DATA;
    }

    if (cutoff <= 0.0 || cutoff >= 180.0) {
        return E_INVARG;
    }

    T = t2 - t1 + 1;
    m = 3 * (n+1);

    /* sample offset */
    y = bw + t1;
    x = x + t1;

    /* the cutoff is expressed in radians internally */
    cutoff *= M_PI / 180.0;

    bad_lambda = set_bw_lambda(cutoff, n, &lam1, &lam2);

#if HAVE_GMP
    if (bad_lambda > 1) {
        gretl_errmsg_set(_("Butterworth: infeasible lambda value"));
        return E_DATA;
    } else if (bad_lambda) {
        return mp_butterworth(x, y, T, n, cutoff);
    }
#else
    if (bad_lambda) {
        gretl_errmsg_set(_("Butterworth: infeasible lambda value"));
        return E_DATA;
    }
#endif

    /* the workspace we need for everything except the
       Toeplitz solver */
    g = malloc(m * sizeof *g);
    if (g == NULL) {
        return E_ALLOC;
    }

    ds = g + n + 1;
    tmp = ds + n + 1;

    /* place a copy of the data in y */
    memcpy(y, x, T * sizeof *y);

    form_wvec(g, ds, tmp, n, lam1, lam2); /* W = M + lambda * Q'SQ */

    QprimeY(y, T);

    /* solve (M + lambda*Q'SQ)x = d for x */
    err = toeplitz_solve(g, y, T-2, n);

    if (!err) {
        form_Qy(y, T);
        form_svec(g, ds, n-2);
        GammaY(g, y, tmp, T, n-2);   /* Form SQy in y */
        /* write the trend into y (low-pass) */
        for (t=0; t<T; t++) {
            y[t] = x[t] - lam1 * y[t];
        }
    }

    free(g);

    return err;
}

/* common code for both plain and weighted polynomial trend
   estimation */

static int real_poly_trend (const double *x, double *fx, double *w,
                            int T, int order)
{
    double tmp, denom, lagdenom = 1;
    double alpha, gamma, delta;
    double *phi, *philag;
    int i, t;

    phi = malloc(2 * T * sizeof *phi);
    if (phi == NULL) {
        return E_ALLOC;
    }

    philag = phi + T;

    for (t=0; t<T; t++) {
        phi[t] = 1;
        philag[t] = 0;
        fx[t] = 0;
    }

    for (i=0; i<=order; i++) {
        double xadd;

        alpha = gamma = denom = 0.0;

        for (t=0; t<T; t++) {
            xadd = x[t] * phi[t];
            if (w != NULL) xadd *= w[t];
            alpha += xadd;
            xadd = phi[t] * phi[t] * t;
            if (w != NULL) xadd *= w[t];
            gamma += xadd;
            xadd = phi[t] * phi[t];
            if (w != NULL) xadd *= w[t];
            denom += xadd;
        }

        alpha /= denom;
        gamma /= denom;
        delta = denom / lagdenom;
        lagdenom = denom;

        for (t=0; t<T; t++) {
            fx[t] += alpha * phi[t];
            tmp = phi[t];
            phi[t] = (t - gamma) * phi[t] - delta * philag[t];
            philag[t] = tmp;
        }
    }

    free(phi);

    return 0;
}

/* Special case of poly_trend() for panel data. */

static int panel_poly_trend (const double *x, double *fx,
                             const DATASET *dset,
                             int order)
{
    int i1 = dset->t1 / dset->pd;
    int i2 = dset->t2 / dset->pd;
    int t0, t1, t2, T;
    int i, t, err = 0;

    for (i=i1; i<=i2 && !err; i++) {
        /* process the data for unit @i */
        t0 = i * dset->pd;
        t1 = 0;
        t2 = dset->pd - 1;
        for (t=0; t<=t2; t++) {
            if (na(x[t0 + t])) {
                t1++;
            } else {
                break;
            }
        }
        for (t=t2; t>=t1; t--) {
            if (na(x[t0 + t])) {
                t2--;
            } else {
                break;
            }
        }
        for (t=t1; t<=t2; t++) {
            if (na(x[t0 + t])) {
                err = E_MISSDATA;
                break;
            }
        }
        /* having found the first and last non-missing
           observations, perform the fitting if possible
        */
        T = t2 - t1 + 1;
        if (order < T) {
            t1 = t0 + t1;
            err = real_poly_trend(x + t1, fx + t1, NULL, T, order);
        } else {
            err = E_DF;
        }
    }

    return err;
}

/**
 * poly_trend:
 * @x: array of original data.
 * @fx: array into which to write the fitted series.
 * @dset: data set information.
 * @order: desired polynomial order.
 *
 * Calculates a trend via the method of orthogonal polynomials.
 * Based on C code for D.S.G. Pollock's DETREND program.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int poly_trend (const double *x, double *fx, const DATASET *dset, int order)
{
    int t1 = dset->t1, t2 = dset->t2;
    int T, err;

    if (dataset_is_panel(dset)) {
        return panel_poly_trend(x, fx, dset, order);
    }

    err = series_adjust_sample(x, &t1, &t2);
    if (err) {
        return err;
    }

    T = t2 - t1 + 1;

    if (order > T) {
        return E_DF;
    }

    return real_poly_trend(x + t1, fx + t1, NULL, T, order);
}

/**
 * poly_weights:
 * @w: array into which the weights will be written.
 * @T: the length of @w.
 * @wmax: the ratio of maximum to minimum weight.
 * @midfrac: the size of the central section that should be given
 * the minimum weight.
 * @opt: weighting scheme option (OPT_Q = quadratic, OPT_B = cosine bell,
 * OPT_C = crenelated).
 *
 * Calculates a set of weights; intended for use with polynomial
 * trend fitting.
 */

void poly_weights (double *w, int T, double wmax,
                   double midfrac, gretlopt opt)
{
    double wt, a = 0, b = 0;
    int t, cut, T2 = T / 2;

    if (midfrac > 0) {
        cut = round(T * (1.0 - midfrac) / 2.0);
    } else {
        cut = T2;
    }

    if (opt == OPT_Q) {
        /* quadratic */
        double z1, z2, det;

        if (midfrac > 0) {
            z1 = cut;
            z2 = 2 * cut;
        } else {
            z1 = (T-1) / 2.0;
            z2 = T-1;
        }

        det = 1.0 / (z1 * (z1 * z2 - z2 * z2));
        a = z2 * (1 - wmax) * det;
        b = -z2 * a;
    } else if (opt == OPT_B) {
        /* cosine-bell */
        b = (wmax - 1) / 2.0;
    }

    for (t=0; t<=T2; t++) {
        if (t <= cut) {
            if (opt == OPT_Q) {
                wt = a * t + b;
                wt = wt * t + wmax;
            } else if (opt == OPT_B) {
                a = (t * M_PI) / cut;
                wt = 1 + b * (1 + cos(a));
            } else {
                /* "crenelated" */
                wt = wmax;
            }
        } else {
            wt = 1;
        }
        w[t] = w[T-1-t] = wt;
    }
}

/**
 * weighted_poly_trend:
 * @x: array of original data.
 * @fx: array into which to write the fitted series.
 * @dset: data set information.
 * @order: desired polynomial order.
 * @opt: weighting option (OPT_Q = quadratic, OPT_B = cosine bell,
 * OPT_C = crenelated).
 * @wratio: ratio of maximum to minimum weight.
 * @midfrac: proportion of the data to be treated specially, in the
 * middle.
 *
 * Calculates a trend via the method of orthogonal polynomials, using
 * the specified weighting scheme.
 * Based on C code for D.S.G. Pollock's DETREND program.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int weighted_poly_trend (const double *x, double *fx, const DATASET *dset,
                         int order, gretlopt opt, double wratio,
                         double midfrac)
{
    double *w = NULL;
    int t1 = dset->t1, t2 = dset->t2;
    int T, err;

    err = series_adjust_sample(x, &t1, &t2);
    if (err) {
        return err;
    }

    T = t2 - t1 + 1;

    if (order > T) {
        return E_DF;
    }

    w = malloc(T * sizeof *w);
    if (w == NULL) {
        err = E_ALLOC;
    } else {
        poly_weights(w, T, wratio, midfrac, opt);
        err = real_poly_trend(x + t1, fx + t1, w, T, order);
        free(w);
    }

    return err;
}

static int n_new_dummies (const DATASET *dset,
                          int nunits, int nperiods)
{
    char dname[VNAMELEN];
    int i, nnew = nunits + nperiods;

    for (i=0; i<nunits; i++) {
        sprintf(dname, "du_%d", i + 1);
        if (gretl_is_series(dname, dset)) {
            nnew--;
        }
    }

    for (i=0; i<nperiods; i++) {
        sprintf(dname, "dt_%d", i + 1);
        if (gretl_is_series(dname, dset)) {
            nnew--;
        }
    }

    return nnew;
}

static const char *dayname (int i)
{
    const char *days[] = {
        N_("Monday"), N_("Tuesday"), N_("Wednesday"),
        N_("Thursday"), N_("Friday"), N_("Saturday"),
        N_("Sunday")
    };

    if (i >= 1 && i <= 7) {
        return days[i-1];
    } else {
        return "??";
    }
}

static gchar *seas_name_and_label (int i, const DATASET *dset,
                                   gretlopt opt, char *vname)
{
    int pd = dataset_is_panel(dset) ? dset->panel_pd : dset->pd;
    int ts = dset->structure == TIME_SERIES || dset->panel_pd > 1;
    gchar *ret = NULL;

    if (opt & OPT_C) {
        sprintf(vname, "S%d", i);
        ret = g_strdup(_("centered periodic dummy"));
    } else if (opt & OPT_S) {
        sprintf(vname, "S%d", i);
        ret = g_strdup(_("uncentered periodic dummy"));
    } else if (pd == 4 && ts) {
        sprintf(vname, "dq%d", i);
        ret = g_strdup_printf(_("= 1 if quarter = %d, 0 otherwise"), i);
    } else if (pd == 12 && ts) {
        sprintf(vname, "dm%d", i);
        ret = g_strdup_printf(_("= 1 if month = %d, 0 otherwise"), i);
    } else if (dated_daily_data(dset)) {
#if 0 /* wait: backward-incompatible */
        sprintf(vname, "dd%d", i);
#else
        sprintf(vname, "dummy_%d", i);
#endif
        ret = g_strdup_printf(_("= 1 if day = %s, 0 otherwise"), _(dayname(i)));
    } else {
        char dumstr[8] = "dummy_";
        char numstr[12];
        int len;

        sprintf(numstr, "%d", i);
        len = strlen(numstr);
        dumstr[8 - len] = '\0';
        sprintf(vname, "%s%d", dumstr, i);
        ret = g_strdup_printf(_("%s = 1 if period is %d, 0 otherwise"), vname, i);
    }

    return ret;
}

static int get_first_panel_period (DATASET *dset)
{
    double x = date_as_double(0, dset->panel_pd, dset->panel_sd0);
    int i, d = ceil(log10(dset->panel_pd));
    int ret = 0;

    x -= floor(x);
    for (i=0; i<d; i++) {
        x *= 10;
    }

    x = (x-floor(x) > 0.5)? ceil(x) : floor(x);
    ret = x - 1;

    return ret;
}

static int real_seasonals (DATASET *dset, int ref, int center,
                           int snames, int **plist)
{
    char vname[VNAMELEN];
    gchar *vlabel = NULL;
    gretlopt opt = 0;
    int *list = NULL;
    int i, vi, k, t, pp;
    int pd, ndums, nnew = 0;
    int pantime = 0;
    int nfix = 0;
    int err = 0;

    if (dset == NULL || dset->n == 0) {
        return E_NODATA;
    }

    if (dataset_is_seasonal_panel(dset)) {
        pd = dset->panel_pd;
        if (pd != 4 && pd != 12 && pd != 24) {
            /* not handled yet */
            return E_PDWRONG;
        } else {
            pantime = 1;
        }
    }

    if (!pantime) {
        pd = dset->pd;
        if (pd < 2 || pd > 99999) {
            return E_PDWRONG;
        }
    }

    if (dated_weekly_data(dset)) {
        /* allow for up to 53 ISO 8601 weeks */
        pd = 53;
    }

    if (ref < 0 || ref > pd) {
        return E_INVARG;
    }

    ndums = pd - (ref > 0);
    list = gretl_list_new(ndums);
    if (list == NULL) {
        return E_ALLOC;
    }

    if (center) opt |= OPT_C;
    if (snames) opt |= OPT_S;

    /* check for appropriate series already in place */

    for (k=1, i=1; i<=pd; i++) {
        if (i == ref) {
            /* dummy not wanted */
            continue;
        }
        vlabel = seas_name_and_label(i, dset, opt, vname);
        vi = series_index(dset, vname);
        if (vi < dset->v) {
            const char *orig = series_get_label(dset, vi);

            list[k] = vi;
            if (orig == NULL || strcmp(vlabel, orig)) {
                series_set_label(dset, vi, vlabel);
                nfix++;
            }
        } else {
            nnew++;
        }
        g_free(vlabel);
        k++;
    }

    if (nnew == 0 && nfix == 0) {
        goto finish;
    }

    /* starting point for new dummy IDs */
    vi = dset->v;

    if (nnew > 0 && dataset_add_series(dset, nnew)) {
        return E_ALLOC;
    }

    for (k=1, i=1; i<=pd; i++) {
        if (i == ref) {
            continue;
        }
        if (list[k] == 0) {
            /* no pre-existing series */
            vlabel = seas_name_and_label(i, dset, opt, vname);
            strcpy(dset->varname[vi], vname);
            series_set_label(dset, vi, vlabel);
            g_free(vlabel);
            list[k] = vi++;
        }
        k++;
    }

    if (pantime) {
        int p0 = get_first_panel_period(dset);
        int T = dset->pd;

        pp = 0;

        for (t=0; t<dset->n; t++) {
            if (t % T == 0) {
                pp = p0;
            } else {
                pp = (t + p0) % pd;
            }
            for (k=0, i=1; i<=list[0]; i++) {
                vi = list[i];
                if (k+1 == ref) k++;
                dset->Z[vi][t] = (pp == k)? 1 : 0;
                k++;
            }
        }
    } else if (dated_daily_data(dset)) {
        char datestr[OBSLEN];
        int wkday;

        for (t=0; t<dset->n; t++) {
            ntolabel(datestr, t, dset);
            wkday = weekday_from_date(datestr);
            for (k=1, i=1; i<=list[0]; i++) {
                vi = list[i];
                if (k == ref) k++;
                dset->Z[vi][t] = (wkday == k)? 1 : 0;
                k++;
            }
        }
    } else if (dataset_is_daily(dset)) {
        int yy, mm = 10;
        double xx;

        pp = pd;
        while ((pp = pp / 10)) {
            mm *= 10;
        }

        for (k=1, i=1; i<=list[0]; i++) {
            vi = list[i];
            if (k == ref) k++;
            for (t=0; t<dset->n; t++) {
                xx = date_as_double(t, pd, dset->sd0) + .1;
                yy = (int) xx;
                pp = (int) (mm * (xx - yy) + 0.5);
                dset->Z[vi][t] = (pp == k)? 1 : 0;
            }
            k++;
        }
    } else if (dated_weekly_data(dset)) {
        char datestr[OBSLEN];
        int wknum;

        for (t=0; t<dset->n; t++) {
            ntolabel(datestr, t, dset);
            wknum = iso_week_from_date(datestr);
            for (k=1, i=1; i<=list[0]; i++) {
                vi = list[i];
                if (k+1 == ref) k++;
                dset->Z[vi][t] = (wknum == k)? 1 : 0;
                k++;
            }
        }
    } else {
        int p0 = get_subperiod(0, dset, NULL);

        for (t=0; t<dset->n; t++) {
            pp = (t + p0) % pd;
            for (k=0, i=1; i<=list[0]; i++) {
                vi = list[i];
                if (k+1 == ref) k++;
                dset->Z[vi][t] = (pp == k)? 1 : 0;
                k++;
            }
        }
    }

    if (center) {
        double cx = 1.0 / pd;

        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            for (t=0; t<dset->n; t++) {
                dset->Z[vi][t] -= cx;
            }
        }
    }

 finish:

    if (plist != NULL) {
        *plist = list;
    } else {
        free(list);
    }

    return err;
}

/**
 * gen_seasonal_dummies:
 * @dset: dataset struct.
 * @ref: 1-based reference period, or 0 for none.
 * @center: if non-zero subtract the population mean from
 * each of the generated dummies.
 *
 * Adds to the data set (if these variables are not already
 * present) a set of seasonal dummy variables.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gen_seasonal_dummies (DATASET *dset, int ref, int center)
{
    return real_seasonals(dset, ref, center, 0, NULL);
}

/**
 * seasonals_list:
 * @dset: dataset struct.
 * @ref: 1-based reference period, or 0 for none.
 * @center: if non-zero, subtract the population mean from
 * each of the generated dummies.
 * @err: location to receive error code.
 *
 * Adds to the data set (if these variables are not already
 * present) a set of seasonal dummy variables.
 *
 * Returns: list holding the ID numbers of the seasonal
 * dummies, or NULL on error.
 */

int *seasonals_list (DATASET *dset, int ref, int center, int *err)
{
    int *list = NULL;

    *err = real_seasonals(dset, ref, center, 1, &list);

    return list;
}

/**
 * gen_panel_dummies:
 * @dset: dataset struct.
 * @opt: %OPT_T for time dummies, otherwise unit dummies.
 * @prn: printer for warning, or NULL.
 *
 * Adds to the data set a set of dummy variables corresponding
 * to either the cross-sectional units in a panel, or the
 * time periods.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gen_panel_dummies (DATASET *dset, gretlopt opt, PRN *prn)
{
    char vname[16], label[MAXLABEL];
    int vi, t, yy, pp, mm;
    int orig_v = dset->v;
    int ndum, nnew;
    int n_unitdum = 0;
    int n_timedum = 0;
    int allnew = 0;
    int newvnum;
    double xx;

    if (opt & OPT_T) {
        ndum = n_timedum = dset->pd;
    } else {
        n_unitdum = dset->n / dset->pd;
        if (dset->n % dset->pd) {
            n_unitdum++;
        }
        ndum = n_unitdum;
    }

    if (ndum == 1) {
        return E_PDWRONG;
    }

    if (complex_subsampled()) {
        nnew = ndum;
        allnew = 1;
    } else {
        nnew = n_new_dummies(dset, n_unitdum, n_timedum);
    }

    if (nnew > 0 && prn != NULL) {
        double sz = nnew * dset->n * 8 / (1024.0 * 1024.0);

        if (sz > 1024) {
            pprintf(prn, _("warning: requested %gMb of storage\n"), sz);
        }
    }

    if (nnew > 0 && dataset_add_series(dset, nnew)) {
        return E_ALLOC;
    }

    pp = dset->pd;
    mm = 10;
    while ((pp = pp / 10)) {
        mm *= 10;
    }

    newvnum = orig_v;

    /* generate time-based dummies, if wanted */

    for (vi=1; vi<=n_timedum; vi++) {
        int dnum;

        sprintf(vname, "dt_%d", vi);

        if (allnew) {
            dnum = newvnum++;
        } else {
            dnum = series_index(dset, vname);
            if (dnum >= orig_v) {
                dnum = newvnum++;
            }
        }

        strcpy(dset->varname[dnum], vname);
        sprintf(label, _("%s = 1 if %s is %d, 0 otherwise"), vname,
                _("period"), vi);
        series_set_label(dset, dnum, label);

        for (t=0; t<dset->n; t++) {
            xx = date_as_double(t, dset->pd, dset->sd0);
            yy = (int) xx;
            pp = (int) (mm * (xx - yy) + 0.5);
            dset->Z[dnum][t] = (pp == vi)? 1.0 : 0.0;
        }
    }

    /* generate unit-based dummies, if wanted */

    for (vi=1; vi<=n_unitdum; vi++) {
        int dmin = (vi - 1) * dset->pd;
        int dmax = vi * dset->pd;
        int dnum;

        sprintf(vname, "du_%d", vi);

        if (allnew) {
            dnum = newvnum++;
        } else {
            dnum = series_index(dset, vname);
            if (dnum >= orig_v) {
                dnum = newvnum++;
            }
        }

        strcpy(dset->varname[dnum], vname);
        sprintf(label, _("%s = 1 if %s is %d, 0 otherwise"), vname,
                _("unit"), vi);
        series_set_label(dset, dnum, label);

        for (t=0; t<dset->n; t++) {
            dset->Z[dnum][t] = (t >= dmin && t < dmax)? 1 : 0;
        }
    }

    return 0;
}

/**
 * gen_unit:
 * @dset: dataset struct.
 * @vnum: location to receive ID number of series, or NULL.
 *
 * (For panel data only) adds to the data set an index variable
 * that uniquely identifies the cross-sectional units.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gen_unit (DATASET *dset, int *vnum)
{
    int xt = 0;
    int i, t;

    if (dset->structure != STACKED_TIME_SERIES) {
        gretl_errmsg_set(_("'genr unit' can be used only with "
                         "panel data"));
        return E_DATA;
    }

    i = series_index(dset, "unit");

    if (i == dset->v && dataset_add_series(dset, 1)) {
        return E_ALLOC;
    }

    strcpy(dset->varname[i], "unit");
    series_set_label(dset, i, _("cross-sectional unit index"));

    for (t=0; t<dset->n; t++) {
        if (t % dset->pd == 0) {
            xt++;
        }
        dset->Z[i][t] = (double) xt;
    }

    if (vnum != NULL) {
        *vnum = i;
    }

    return 0;
}

/**
 * panel_unit_first_obs:
 * @t: zero-based observation number.
 * @dset: data information struct.
 *
 * Returns: 1 if observation @t is the first time-series
 * observation on a given cross-sectional unit in a
 * panel dataset, else 0.
 */

int panel_unit_first_obs (int t, const DATASET *dset)
{
    char *p, obs[OBSLEN];
    int ret = 0;

    ntolabel(obs, t, dset);
    p = strchr(obs, ':');
    if (p != NULL && atoi(p + 1) == 1) {
        ret = 1;
    }

    return ret;
}

/* make special time variable for panel data */

void make_panel_time_var (double *x, const DATASET *dset)
{
    int t, xt = 0;

    for (t=0; t<dset->n; t++) {
        if (t % dset->pd == 0) {
            xt = 1;
        }
        x[t] = (double) xt++;
    }
}

/**
 * gen_time:
 * @dset: dataset struct.
 * @tm: if non-zero, an actual time trend is wanted,
 * otherwise just an index of observations.
 * @vnum: location to receive ID number of series,
 * or NULL.
 *
 * Generates (and adds to the dataset, if it's not already
 * present) a time-trend or index variable.  This function
 * is panel-data aware: if the dataset is a panel and
 * @tm is non-zero, the trend will not simply run
 * consecutively over the entire range of the data, but
 * will correctly represent the location in time of
 * each observation.  The index is 1-based.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gen_time (DATASET *dset, int tm, int *vnum)
{
    int v, t;

    v = series_index(dset, (tm)? "time" : "index");

    if (v == dset->v && dataset_add_series(dset, 1)) {
        return E_ALLOC;
    }

    if (tm) {
        strcpy(dset->varname[v], "time");
        series_set_label(dset, v, _("time trend variable"));
    } else {
        strcpy(dset->varname[v], "index");
        series_set_label(dset, v, _("data index variable"));
    }

    if (tm && dset->structure == STACKED_TIME_SERIES) {
        make_panel_time_var(dset->Z[v], dset);
    } else {
        for (t=0; t<dset->n; t++) {
            dset->Z[v][t] = (double) (t + 1);
        }
    }

    if (vnum != NULL) {
        *vnum = v;
    }

    return 0;
}

/**
 * genr_wkday:
 * @dset: dataset struct.
 * @vnum: location to receive ID number of series, or NULL.
 *
 * Generates (and adds to the dataset, if it's not already
 * present) an index representing the day of the week for
 * each observation (for dated daily data only).
 * The index has value 0 for Sunday, 1 for Monday, and
 * so on.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gen_wkday (DATASET *dset, int *vnum)
{
    char datestr[OBSLEN];
    int i, t, wd;

    if (!dated_daily_data(dset)) {
        return E_PDWRONG;
    }

    i = series_index(dset, "weekday");
    if (i == dset->v && dataset_add_series(dset, 1)) {
        return E_ALLOC;
    }

    strcpy(dset->varname[i], "weekday");
    series_set_label(dset, i, _("day of week (1 = Monday)"));

    for (t=0; t<dset->n; t++) {
        ntolabel(datestr, t, dset);
        wd = weekday_from_date(datestr);
        /* traditional userspace weekdays */
        wd = (wd == G_DATE_SUNDAY) ? 0 : wd;
        dset->Z[i][t] = (double) wd;
    }

    if (vnum != NULL) {
        *vnum = i;
    }

    return 0;
}

typedef enum {
    PLOTVAR_INDEX = 1,
    PLOTVAR_TIME,
    PLOTVAR_ANNUAL,
    PLOTVAR_QUARTERS,
    PLOTVAR_MONTHS,
    PLOTVAR_CALENDAR,
    PLOTVAR_GPTIME,
    PLOTVAR_DECADES,
    PLOTVAR_HOURLY,
    PLOTVAR_PANEL,
    PLOTVAR_MAX
} plotvar_type;

static int plotvar_code (const DATASET *dset, gretlopt opt)
{
    if ((opt & OPT_S) && calendar_data(dset)) {
        return PLOTVAR_TIME;
    } else if (!dataset_is_time_series(dset)) {
        return PLOTVAR_INDEX;
    } else if (dset->pd == 1) {
        return PLOTVAR_ANNUAL;
    } else if (dset->pd == 4) {
        return PLOTVAR_QUARTERS;
    } else if (dset->pd == 12) {
        return PLOTVAR_MONTHS;
    } else if (dset->pd == 24) {
        return PLOTVAR_HOURLY;
    } else if (calendar_data(dset)) {
        return (opt & OPT_T)? PLOTVAR_GPTIME : PLOTVAR_CALENDAR;
    } else if (dataset_is_decennial(dset)) {
        return PLOTVAR_DECADES;
    } else {
        return PLOTVAR_TIME;
    }
}

static int panel_plotvar_code (const DATASET *dset)
{
    int pd = dset->panel_pd;

    if (pd == 0) {
        return 0;
    } else if (pd == 1) {
        return PLOTVAR_ANNUAL;
    } else if (pd == 4) {
        return PLOTVAR_QUARTERS;
    } else if (pd == 12) {
        return PLOTVAR_MONTHS;
    } else if (pd == 24) {
        return PLOTVAR_HOURLY;
    } else if ((pd == 5 || pd == 6 || pd == 7 || pd == 52) &&
               dset->panel_sd0 > 10000.0) {
        return PLOTVAR_CALENDAR;
    } else if (pd == 10) {
        return PLOTVAR_DECADES;
    } else {
        return 0;
    }
}

static double correct_to_int (double x)
{
    if (x - floor(x) < 1.0e-6) {
        return floor(x);
    } else if (ceil(x) - x < 1.0e-6) {
        return ceil(x);
    }
    return x;
}

static double unix_time_from_iso_date (const char *s)
{
    double x, ret = NADBL;
    char *test;

    test = gretl_strptime(s, "%Y-%m-%d", &x);
    if (test != NULL && *test == '\0') {
        ret = x;
    }

    return ret;
}

/**
 * gretl_plotx:
 * @dset: data information struct.
 * @opt: can include OPT_P for panel time-series plot;
 * OPT_T to use gnuplot time (seconds of unix epoch);
 * OPT_S to indicate that the context is the "scatters"
 * command.
 *
 * Finds or creates a special dummy variable for use on the
 * x-axis in plotting; this will have the full length of the
 * data series as given in @dset, and will be appropriately
 * configured for the data frequency.  Do not try to free this
 * variable.
 *
 * Returns: pointer to plot x-variable, or NULL on failure.
 */

const double *gretl_plotx (const DATASET *dset, gretlopt opt)
{
    static double *x;
    static int ptype;
    static int Tbak;
    static double sd0bak;
    char datestr[OBSLEN];
    char *tstr = NULL;
    double mul, frac;
    int t, y1, T = 0;
    int maj, min;
    int new_ptype = 0;
    int panvar = 0;
    int failed = 0;
    double sd0 = 0;
    float rm;

    if (dset == NULL) {
        /* cleanup signal */
        free(x);
        x = NULL;
        ptype = 0;
        T = 0;
        sd0 = 0;
        return NULL;
    }

    if (dataset_is_panel(dset) && ((opt & OPT_P) ||
                                   sample_size(dset) == dset->pd)) {
        /* we're plotting a single time-series from a panel */
        new_ptype = panel_plotvar_code(dset);
        if (new_ptype != 0) {
            sd0 = dset->panel_sd0;
        } else {
            panvar = plausible_panel_time_var(dset);
            if (panvar > 0) {
                new_ptype = PLOTVAR_PANEL;
            } else {
                new_ptype = PLOTVAR_INDEX;
            }
        }
        T = dset->pd;
    }

    if (new_ptype == 0) {
        /* not already determined */
        new_ptype = plotvar_code(dset, opt);
    }
    if (T == 0) {
        /* not already determined */
        T = dset->n;
    }
    if (sd0 == 0.0) {
        /* not already determined */
        sd0 = dset->sd0;
    }

    if (x != NULL && new_ptype == ptype && Tbak == T && sd0 == sd0bak) {
        /* a suitable array is already at hand */
        return x;
    }

    if (x != NULL) {
        free(x);
    }

    if (new_ptype == PLOTVAR_PANEL) {
        x = copyvec(dset->Z[panvar], dset->n);
    } else {
        x = malloc(T * sizeof *x);
    }

    if (x == NULL || new_ptype == PLOTVAR_PANEL) {
        return x;
    }

 try_again:

    Tbak = T;
    ptype = new_ptype;
    sd0bak = sd0;

    y1 = (int) sd0;
    rm = sd0 - y1;

    switch (ptype) {
    case PLOTVAR_ANNUAL:
        for (t=0; t<T; t++) {
            x[t] = sd0 + t;
        }
        break;
    case PLOTVAR_QUARTERS:
    case PLOTVAR_MONTHS:
    case PLOTVAR_HOURLY:
        frac = 1.0 / dset->pd;
        if (sscanf(dset->stobs, "%d:%d", &maj, &min) == 2) {
            for (t=0; t<T; t++) {
                x[t] = maj + frac * (min - 1);
                if (min < dset->pd) {
                    min++;
                } else {
                    min = 1;
                    maj++;
                }
            }
        } else {
            mul = ptype == PLOTVAR_QUARTERS ? 10 : 100;
            x[0] = correct_to_int(y1 + (mul * rm - 1.0) / dset->pd);
            for (t=1; t<T; t++) {
                x[t] = correct_to_int(x[t-1] + frac);
            }
        }
        break;
    case PLOTVAR_CALENDAR:
    case PLOTVAR_GPTIME:
        for (t=0; t<T; t++) {
            if (dset->S != NULL) {
                tstr = dset->S[t];
            } else {
                calendar_date_string(datestr, t, dset);
                tstr = datestr;
            }
            if (ptype == PLOTVAR_GPTIME) {
                x[t] = unix_time_from_iso_date(tstr);
            } else {
                x[t] = get_dec_date(tstr);
            }
            if (na(x[t])) {
                failed++;
                break;
            }
        }
        break;
    case PLOTVAR_DECADES:
        for (t=0; t<T; t++) {
            x[t] = dset->sd0 + 10 * t;
        }
        break;
    case PLOTVAR_INDEX:
    case PLOTVAR_TIME:
        for (t=0; t<T; t++) {
            x[t] = (double) (t + 1);
        }
        break;
    default:
        break;
    }

    if (failed == 1) {
        /* calendar dating failed */
        failed++; /* ensure we don't loop */
        new_ptype = PLOTVAR_TIME;
        goto try_again;
    }

    return x;
}

/**
 * get_fit_or_resid:
 * @pmod: pointer to source model.
 * @dset: information on the data set.
 * @idx: %M_UHAT, %M_UHAT2, %M_YHAT, %M_AHAT or %M_H.
 * @vname: location to write series name (length %VNAMELEN)
 * @pdesc: location to receive copy of series description.
 * @err: location to receive error code.
 *
 * Creates a full-length array holding the specified model
 * data, and writes name and description into the @vname and
 * @vlabel.
 *
 * Returns: allocated array on success or NULL on failure.
 */

double *get_fit_or_resid (const MODEL *pmod, DATASET *dset,
                          ModelDataIndex idx, char *vname,
                          gchar **pdesc, int *err)
{
    const double *src = NULL;
    double *ret = NULL;
    int t;

    /* Note: this is called only from the GUI, and in that
       context @vname will be shown as the default name
       for the saved series, though the user can change
       it if she chooses.
    */

    if (idx == M_H) {
        src = gretl_model_get_data(pmod, "garch_h");
    } else if (idx == M_AHAT) {
        src = gretl_model_get_data(pmod, "ahat");
    } else if (idx == M_UHAT || idx == M_UHAT2) {
        src = pmod->uhat;
    } else if (idx == M_YHAT) {
        src = pmod->yhat;
    }

    if (src == NULL) {
        if (*err == 0) {
            *err = E_BADSTAT;
        }
        return NULL;
    }

    ret = malloc(dset->n * sizeof *ret);
    if (ret == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    for (t=0; t<dset->n; t++) {
        if (t >= pmod->t1 && t <= pmod->t2) {
            if (idx == M_UHAT2) {
                ret[t] = na(src[t]) ? NADBL : (src[t] * src[t]);
            } else {
                ret[t] = src[t];
            }
        } else {
            ret[t] = NADBL;
        }
    }

    if (idx == M_UHAT) {
        sprintf(vname, "uhat%d", pmod->ID);
        if (pmod->ci == GARCH && (pmod->opt & OPT_Z)) {
            *pdesc = g_strdup_printf(_("standardized residual from model %d"), pmod->ID);
        } else {
            *pdesc = g_strdup_printf(_("residual from model %d"), pmod->ID);
        }
    } else if (idx == M_YHAT) {
        sprintf(vname, "yhat%d", pmod->ID);
        *pdesc = g_strdup_printf(_("fitted value from model %d"), pmod->ID);
    } else if (idx == M_UHAT2) {
        /* squared residuals */
        sprintf(vname, "usq%d", pmod->ID);
        if (pmod->ci == GARCH && (pmod->opt & OPT_Z)) {
            *pdesc = g_strdup_printf(_("squared standardized residual from model %d"), pmod->ID);
        } else {
            *pdesc = g_strdup_printf(_("squared residual from model %d"), pmod->ID);
        }
    } else if (idx == M_H) {
        /* garch variance */
        sprintf(vname, "h%d", pmod->ID);
        *pdesc = g_strdup_printf(_("fitted variance from model %d"), pmod->ID);
    } else if (idx == M_AHAT) {
        sprintf(vname, "ahat%d", pmod->ID);
        if (pmod->opt & OPT_U) {
            *pdesc = g_strdup_printf(_("individual effects from model %d"), pmod->ID);
        } else {
            *pdesc = g_strdup_printf(_("per-unit constants from model %d"), pmod->ID);
        }
    }

    return ret;
}

/**
 * genr_fit_resid:
 * @pmod: pointer to source model.
 * @dset: dataset struct.
 * @idx: %M_UHAT, %M_UHAT2, %M_YHAT, %M_AHAT or %M_H.
 *
 * Adds residuals or fitted values or squared residuals from a
 * given model to the data set.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int genr_fit_resid (const MODEL *pmod, DATASET *dset,
                    ModelDataIndex idx)
{
    char vname[VNAMELEN];
    gchar *vlabel = NULL;
    double *x;
    int err = 0;

    x = get_fit_or_resid(pmod, dset, idx, vname, &vlabel, &err);

    if (!err) {
        err = dataset_add_allocated_series(dset, x);
    }

    if (err) {
        free(x);
    } else {
        int v = dset->v - 1;

        strcpy(dset->varname[v], vname);
        series_set_label(dset, v, vlabel);
    }

    g_free(vlabel);

    return err;
}

int get_observation_number (const char *s, const DATASET *dset)
{
    char test[OBSLEN];
    size_t n;
    int t;

    *test = 0;
    strncat(test, (*s == '"')? s + 1 : s, OBSLEN - 1);

    n = strlen(test);
    if (test[n-1] == '"') {
        test[n-1] = '\0';
    }

    if (dataset_has_markers(dset)) {
        for (t=0; t<dset->n; t++) {
            if (!strcmp(test, dset->S[t])) {
                return t + 1;
            }
        }
        if (calendar_data(dset)) {
            for (t=0; t<dset->n; t++) {
                if (!strcmp(test, dset->S[t]) ||
                    !strcmp(test, dset->S[t] + 2)) {
                    return t + 1;
                }
            }
        }
    }

    if (dset->structure == TIME_SERIES) {
        t = dateton(test, dset);
        if (t >= 0) {
            return t + 1;
        }
    }

    if (calendar_data(dset)) {
        /* never reached? */
        char datestr[OBSLEN];

        for (t=0; t<dset->n; t++) {
            calendar_date_string(datestr, t, dset);
            if (!strcmp(test, datestr) ||
                !strcmp(test, datestr + 2)) {
                return t + 1;
            }
        }
        return 0;
    }

    return 0;
}

#define OBS_DEBUG 0

static int plain_obs_number (const char *obs, const DATASET *dset)
{
    char *test;
    int t = -1;

    errno = 0;

    strtol(obs, &test, 10);

    if (errno == 0 && *test == '\0') {
        t = atoi(obs) - 1; /* convert from 1-based to 0-based */
        if (t >= dset->n) {
            t = -1;
        }
    }

    return t;
}

/* Given what looks like an observation number or date within "[" and
   "]", try to determine the observation number.  This is quite tricky
   since we try to handle both dates and plain observation numbers
   (and in addition, variables representing the latter); and we may
   have to deal with the translation from 1-based indexing in user
   space to 0-based indexing for internal purposes.
*/

int get_t_from_obs_string (const char *s, const DATASET *dset)
{
    int t;

    if (*s == '"') {
        char obs[OBSLEN+2];
        int err = 0;

        *obs = '\0';
        strncat(obs, s, OBSLEN+1);
        gretl_unquote(obs, &err);
        t = dateton(obs, dset);
    } else {
        t = dateton(s, dset);
    }

#if OBS_DEBUG
    fprintf(stderr, "\nget_t_from_obs_string: s ='%s', dateton gives t = %d\n",
            s, t);
#endif

    if (t < 0) {
        if (isdigit((unsigned char) *s)) {
            t = plain_obs_number(s, dset);
#if OBS_DEBUG
            fprintf(stderr, " plain_obs_number gives t = %d\n", t);
#endif
        } else {
            if (gretl_is_scalar(s)) {
                t = gretl_scalar_get_value(s, NULL);
            }
            if (t > dset->n) {
                /* e.g. annual dates */
                char try[16];

                sprintf(try, "%d", t);
                t = dateton(try, dset);
#if OBS_DEBUG
                fprintf(stderr, " revised via dateton: t = %d\n", t);
#endif
            } else {
                /* convert to 0-based */
                t--;
            }
        }
    }

    if (t < 0) {
        gretl_errmsg_set(_("Observation number out of bounds"));
    }

#if OBS_DEBUG
    fprintf(stderr, " return value: t = %d\n", t);
#endif

    return t;
}

int check_declarations (char ***pS, parser *p)
{
    char **S;
    const char *s;
    int exists = 0;
    int badname = 0;
    int i, n = 1;

    gretl_error_clear();

    if (p->lh.expr == NULL) {
        p->err = E_ALLOC;
        return 0;
    }

    s = p->lh.expr;
    s += strspn(s, " ");

    while (*s) {
        if (*s == ',' || *s == ' ') {
            n++;
            s++;
            s += strspn(s, " ");
        } else {
            s++;
        }
    }

    S = strings_array_new(n);
    if (S == NULL) {
        p->err = E_ALLOC;
        return 0;
    }

    s = p->lh.expr;
    for (i=0; i<n && !p->err; i++) {
        S[i] = gretl_word_strdup(s, &s, OPT_S | OPT_U, &p->err);
    }

    if (!p->err && *s != '\0') {
        p->err = E_DATA;
    }

    for (i=0; i<n && !p->err; i++) {
        if (gretl_type_from_name(S[i], p->dset)) {
            /* variable already exists */
            exists = 1;
            p->err = E_DATA;
        } else if (check_identifier(S[i])) {
            /* invalid name */
            badname = 1;
            p->err = E_DATA;
        }
    }

    if (p->err) {
        if (exists) {
            gretl_errmsg_set(_("Invalid declaration: maybe you need "
                               "the \"clear\" command?"));
        } else if (!badname) {
            gretl_errmsg_set(_("Invalid declaration"));
        }
        strings_array_free(S, n);
    } else {
        *pS = S;
    }

    return n;
}

/* cross-sectional (list-based) statistics */

static double xsect_minmax (double *x, int n, int f)
{
    double ret = x[0];
    int i;

    for (i=1; i<n; i++) {
        if (f == F_MIN) {
            if (x[i] < ret) {
                ret = x[i];
            }
        } else if (x[i] > ret) {
            ret = x[i];
        }
    }

    return ret;
}

/* Fill @targ with data (from @xlist) and associated
   weights (from @wlist) for the cross-section at obs @t.
   If @partial_ok is 0 (false) we return E_MISSDATA if
   any values or weights are missing at @t, otherwise
   we skip missing terms and return the number of valid
   terms in @pn.
*/

static int data_and_weight_at_obs (double **targ, int *pn,
                                   const int *xlist,
                                   const int *wlist,
                                   const DATASET *dset,
                                   int t, int partial_ok)
{
    int i, j = 0;
    double xi, wi;

    for (i=1; i<=xlist[0]; i++) {
        xi = dset->Z[xlist[i]][t];
        wi = dset->Z[wlist[i]][t];
        if (wi < 0) {
            return E_INVARG;
        } else if (wi == 0) {
            continue;
        } else if (wi > 0 && na(xi)) {
            if (!partial_ok) {
                return E_MISSDATA;
            }
        } else if (!na(xi) && !na(wi)) {
            targ[0][j] = xi;
            targ[1][j] = wi;
            j++;
        } else if (!partial_ok) {
            return E_MISSDATA;
        }
    }

    *pn = j;
    return (j == 0)? E_MISSDATA : 0;
}

/* data_at_obs: returns a count of the non-missing values of
   the members of @list at observation @t. If @targ is non-NULL,
   enter valid values into this array.
*/

static int data_at_obs (const int *list, int t,
                        const DATASET *dset,
			double *targ)
{
    int i, j = 0;
    double xit;

    for (i=1; i<=list[0]; i++) {
        xit = dset->Z[list[i]][t];
	if (!na(xit)) {
	    if (targ != NULL) {
		targ[j] = xit;
	    }
	    j++;
	}
    }

    return j;
}

int cross_sectional_stat (double *y, const int *list,
                          const DATASET *dset, int f,
                          int partial_ok)
{
    double *x = NULL;
    int nv = list[0];
    int i, t, n;

    if (f != F_NOBS) {
	/* storage for the cross section, if needed */
	x = malloc(nv * sizeof *x);
	if (x == NULL) {
	    return E_ALLOC;
	}
    }

    for (t=dset->t1; t<=dset->t2; t++) {
        n = data_at_obs(list, t, dset, x);
	if (f == F_NOBS) {
	    y[t] = (double) n;
        } else if (n == 0 || (n < nv && !partial_ok)) {
            y[t] = NADBL;
        } else if ((f == F_VCE || f == F_SD) && n == 1) {
            y[t] = NADBL;
        } else if (f == F_MIN || f == F_MAX) {
            y[t] = xsect_minmax(x, n, f);
        } else if (f == F_MEDIAN) {
            y[t] = gretl_median(0, n-1, x);
        } else {
            double xsum = 0;

            for (i=0; i<n; i++) {
                xsum += x[i];
            }
            if (f == F_SUM) {
                y[t] = xsum;
            } else if (f == F_MEAN) {
                y[t] = xsum / n;
            } else {
                /* F_VCE or F_SD */
                double d, s2 = 0, xbar = xsum / n;

                for (i=0; i<n; i++) {
                    d = x[i] - xbar;
                    s2 += d * d;
                }
                if (f == F_VCE) {
                    y[t] = s2 / (n-1);
                } else {
                    y[t] = sqrt(s2 / (n-1));
                }
            }
        }
    }

    free(x);

    return 0;
}

int x_sectional_weighted_stat (double *y, const int *xlist,
                               const int *wlist,
                               const DATASET *dset,
                               int f, int partial_ok)
{
    double wxbar, wsum, d, ws2, V2, adj;
    double *x, *w, **X;
    int i, t, n, missing;
    int err = 0;

    if (wlist[0] != xlist[0]) {
        gretl_errmsg_sprintf(_("Weighted stats: data list has %d members but weight "
                             "list has %d"), xlist[0], wlist[0]);
        return E_DATA;
    }

    /* storage for values and weights, per obs */
    X = doubles_array_new(2, xlist[0]);
    if (X == NULL) {
        return E_ALLOC;
    }

    x = X[0];
    w = X[1];

    for (t=dset->t1; t<=dset->t2; t++) {
        missing = data_and_weight_at_obs(X, &n, xlist, wlist,
                                         dset, t, partial_ok);
        if (missing == E_INVARG) {
            err = E_INVARG;
            break;
        } else if (missing) {
            y[t] = NADBL;
        } else if (f != F_WMEAN && n == 1) {
            y[t] = NADBL;
        } else {
            wxbar = wsum = 0;
            for (i=0; i<n; i++) {
                wxbar += x[i] * w[i];
                wsum += w[i];
            }
            wxbar /= wsum;
            if (f == F_WMEAN) {
                y[t] = wxbar;
            } else {
                ws2 = V2 = 0;
                for (i=0; i<n; i++) {
                    d = (x[i] - wxbar);
                    ws2 += w[i] * d * d;
                    V2 += w[i] * w[i];
                }
                /* "frequency" weights */
                adj = (n-1) * wsum / n;
                /* for "reliability" weights? */
                // adj = wsum - V2 / wsum;
                if (f == F_WVAR) {
                    y[t] = ws2 / adj;
                } else {
                    y[t] = sqrt(ws2 / adj);
                }
            }
        }
    }

    doubles_array_free(X, 2);

    return err;
}

/* writes to the series @y a linear combination of the variables given
   in @list, using the coefficients given in the vector @b.
*/

int list_linear_combo (double *y, const int *list,
                       const gretl_vector *b,
                       const DATASET *dset)
{
    int nb = gretl_vector_get_length(b);
    int nl = list[0];
    int err = 0;

    if (nb != nl) {
        gretl_errmsg_sprintf(_("List has %d members, but length "
                               "of vector b is %d"), nl, nb);
        err = E_DATA;
    } else {
        int i, t;
        double xit, yt;

        for (t=dset->t1; t<=dset->t2; t++) {
            yt = 0;
            for (i=0; i<nl; i++) {
                xit = dset->Z[list[i+1]][t];
                if (na(xit)) {
                    yt = NADBL;
                    break;
                } else {
                    yt += xit * gretl_vector_get(b, i);
                }
            }
            y[t] = yt;
        }
    }

    return err;
}

#define MDEBUG 0

#if HAVE_GMP

static int try_mp_midas_weights (const double *theta, int k,
                                 gretl_matrix *w, int method)
{
    int (*mpfun) (const double *, int, gretl_matrix *, int);
    int err = 0;

    mpfun = get_plugin_function("mp_midas_weights");

    if (mpfun == NULL) {
        return E_FOPEN;
    }

    err = (*mpfun) (theta, k, w, method);

    if (!err) {
        int i, n = gretl_vector_get_length(w);

        for (i=0; i<n; i++) {
            if (!isfinite(w->val[i])) {
                err = E_NAN;
                break;
            }
        }
    }

    return err;
}

#endif /* HAVE_GMP */

#define INF_CHECK_VERBOSE 0

static int inf_check (double val, const gretl_matrix *theta,
                      const char *context, int *err)
{
    if (isinf(val)) {
#if INF_CHECK_VERBOSE
        fprintf(stderr, "range error in %s\n", context);
        gretl_matrix_print(theta, "theta");
#endif
        *err = E_NAN;
    } else {
        /* tolerate underflow */
        errno = 0;
    }

    return *err;
}

#define BETA_METHOD(m) (m==MIDAS_BETA0 || m==MIDAS_BETA1 || m==MIDAS_BETAN)

static int check_beta_params (int method, double *theta,
                              int k, double eps)
{
    int err = 0;

    if ((method == MIDAS_BETA0 || method == MIDAS_BETA1) && k != 2) {
        gretl_errmsg_set(_("theta must be a 2-vector"));
        err = E_INVARG;
    } else if (method == MIDAS_BETAN && k != 3) {
        gretl_errmsg_set(_("theta must be a 3-vector"));
        err = E_INVARG;
    } else if (theta[0] < eps || theta[1] < eps) {
        if (theta[0] < 0.0 || theta[1] < 0.0) {
            gretl_errmsg_set(_("beta: parameters must be positive"));
            fprintf(stderr, "beta: theta1=%g, theta2=%g\n", theta[0], theta[1]);
            err = E_INVARG;
        } else {
            /* bodge! */
            if (theta[0] < eps) theta[0] = eps;
            if (theta[1] < eps) theta[1] = eps;
        }
    }

    return err;
}

/* Computes a column m-vector holding weights for use with MIDAS */

gretl_matrix *midas_weights (int p, const gretl_matrix *m,
                             int method, int *err)
{
    double *theta;
    double eps = pow(2.0, -52);
    double wsum = 0.0;
    gretl_matrix *w = NULL;
    int using_mp = 0;
    int k, i, j;

#if MDEBUG
    gretl_matrix_print(m, "m, in midas_weights");
#endif

    /* @p = lag order
       @k = number of hyper-parameters
    */

    if (method < MIDAS_NEALMON || method >= MIDAS_MAX || p < 1) {
        *err = E_INVARG;
        return NULL;
    }

    k = gretl_vector_get_length(m);
    if (k == 0) {
        *err = E_INVARG;
        return NULL;
    }

    theta = m->val;

    if (BETA_METHOD(method)) {
        *err = check_beta_params(method, theta, k, eps);
        if (*err) {
            return NULL;
        }
    }

    w = gretl_zero_matrix_new(p, 1);
    if (w == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    errno = 0;

    if (method == MIDAS_NEALMON) {
        for (i=0; i<p; i++) {
            w->val[i] = (i+1) * theta[0];
            for (j=1; j<k; j++) {
                w->val[i] += pow(i+1, j+1) * theta[j];
            }
            w->val[i] = exp(w->val[i]);
            wsum += w->val[i];
            if (errno && inf_check(wsum, m, "nealmon weights", err)) {
                break;
            }
        }
    } else if (BETA_METHOD(method)) {
        double si, ai, bi;

        for (i=0; i<p; i++) {
            si = i / (p - 1.0);
            if (i == 0) {
                si += eps;
            } else if (i == p-1) {
                si -= eps;
            }
            ai = pow(si, theta[0] - 1.0);
            bi = pow(1.0 - si, theta[1] - 1.0);
            w->val[i] = ai * bi;
            wsum += w->val[i];
            if (errno && inf_check(wsum, m, "beta weights", err)) {
                break;
            }
        }
    } else if (method == MIDAS_ALMONP) {
        /* straight Almon ploynomial */
        for (i=0; i<p; i++) {
            w->val[i] = theta[0];
            for (j=1; j<k; j++) {
                w->val[i] += pow(i+1, j) * theta[j];
            }
        }
    }

#if HAVE_GMP
    if (*err == E_NAN || (method != MIDAS_ALMONP && wsum < eps)) {
        /* attempt a fix-up using multiple precision */
        int save_errno = errno;

        errno = 0;
        *err = try_mp_midas_weights(theta, k, w, method);
        if (*err) {
            errno = save_errno;
        } else {
            using_mp = 1;
        }
    }
#endif

    if (errno) {
        gretl_errmsg_sprintf(_("Failed to calculate MIDAS weights: %s"),
                             gretl_strerror(errno));
        if (*err == 0) {
            *err = E_INVARG;
        }
        gretl_matrix_free(w);
        errno = 0;
        return NULL;
    }

    if (!using_mp && method != MIDAS_ALMONP) {
        /* normalize the weights */
        for (i=0; i<p; i++) {
            w->val[i] /= wsum;
        }
    }

    if (!using_mp && method == MIDAS_BETAN) {
        /* beta with third param, not zero-terminated:
           add theta[2] and renormalize
        */
        wsum = 1 + p * theta[2];
        for (i=0; i<p; i++) {
            w->val[i] = (w->val[i] + theta[2]) / wsum;
        }
    }

#if MDEBUG
    gretl_matrix_print(w, "w: midas weights");
#endif

    return w;
}

#if HAVE_GMP

static int try_mp_midas_grad (const double *theta,
                              gretl_matrix *G,
                              int method)
{
    int (*mpfun) (const double *, gretl_matrix *, int);
    int err = 0;

    mpfun = get_plugin_function("mp_midas_gradient");

    if (mpfun == NULL) {
        return E_FOPEN;
    }

    err = (*mpfun) (theta, G, method);

    if (!err) {
        int i, n = G->rows * G->cols;

        for (i=0; i<n; i++) {
            if (!isfinite(G->val[i])) {
                err = E_NAN;
                break;
            }
        }
    }

    return err;
}

static int mgrad_zero (const gretl_matrix *G)
{
    int i, j, colzero;

    for (j=0; j<G->cols; j++) {
        colzero = 1;
        for (i=0; i<G->rows; i++) {
            if (gretl_matrix_get(G, i, j) != 0.0) {
                colzero = 0;
                break;
            }
        }
        if (colzero) {
            return 1;
        }
    }

    return 0;
}

#endif /* HAVE_GMP */

gretl_matrix *midas_gradient (int p, const gretl_matrix *m,
                              int method, int *err)
{
    double *theta;
    double eps = pow(2.0, -52);
    double ws2, wsum = 0.0;
    gretl_matrix *w = NULL;
    gretl_matrix *G = NULL;
    int k, i, j;

#if MDEBUG
    gretl_matrix_print(m, "m, in midas_gradient");
#endif

    if (method < MIDAS_NEALMON || method >= MIDAS_MAX || p < 1) {
        *err = E_INVARG;
        return NULL;
    }

    /* @p = lag order
       @k = number of hyper-parameters
    */

    k = gretl_vector_get_length(m);
    if (k == 0) {
        *err = E_INVARG;
        return NULL;
    }

    theta = m->val;
    errno = 0;

    if (BETA_METHOD(method)) {
        *err = check_beta_params(method, theta, k, eps);
        if (*err) {
            return NULL;
        }
    }

    if (method != MIDAS_ALMONP) {
        w = gretl_zero_matrix_new(p, 1);
        if (w == NULL) {
            *err = E_ALLOC;
            return NULL;
        }
    }

    G = gretl_matrix_alloc(p, k);
    if (G == NULL) {
        gretl_matrix_free(w);
        *err = E_ALLOC;
        return NULL;
    }

    if (method == MIDAS_NEALMON) {
        double *dsum = malloc(k * sizeof *dsum);
        double gij;

        for (j=0; j<k; j++) {
            dsum[j] = 0.0;
        }
        for (i=0; i<p; i++) {
            w->val[i] = (i+1) * theta[0];
            for (j=1; j<k; j++) {
                w->val[i] += pow(i+1, j+1) * theta[j];
            }
            w->val[i] = exp(w->val[i]);
            wsum += w->val[i];
            if (errno && inf_check(wsum, m, "nealmon gradient", err)) {
                free(dsum);
                goto range_error;
            }
        }
        for (i=0; i<p; i++) {
            for (j=0; j<k; j++) {
                dsum[j] += w->val[i] * pow(i+1, j+1);
            }
        }
        for (j=0; j<k; j++) {
            dsum[j] /= wsum;
        }
        for (i=0; i<p; i++) {
            w->val[i] /= wsum;
            for (j=0; j<k; j++) {
                gij = w->val[i] * (pow(i+1, j+1) - dsum[j]);
                gretl_matrix_set(G, i, j, gij);
            }
        }
        free(dsum);
    } else if (BETA_METHOD(method)) {
        double si, ai, bi;
        double g1sum = 0;
        double g2sum = 0;

        /* loop 1: form raw weights and their sum */
        for (i=0; i<p; i++) {
            si = i / (p - 1.0);
            if (i == 0) {
                si += eps;
            } else if (i == p - 1) {
                si -= eps;
            }
            ai = pow(si, theta[0] - 1.0);
            bi = pow(1.0 - si, theta[1] - 1.0);
            w->val[i] = ai * bi;
            wsum += w->val[i];
            if (errno && inf_check(wsum, m, "beta gradient", err)) {
                goto range_error;
            }
        }
        if (wsum <= eps) {
            /* should we just set G to zero in this case? */
#if 0
            fprintf(stderr, "sum of weights = %g\n", wsum);
#endif
            *err = E_NAN;
            goto range_error;
        }
        ws2 = wsum * wsum;
        if (errno && inf_check(ws2, m, "beta gradient", err)) {
            goto range_error;
        }
        /* loop 2: form first component of derivative
           and second cumulant */
        for (i=0; i<p; i++) {
            si = i / (p - 1.0);
            if (i == 0) {
                si += eps;
            } else if (i == p - 1) {
                si -= eps;
            }
            ai = w->val[i] * log(si);
            g1sum += ai;
            gretl_matrix_set(G, i, 0, ai/wsum);
            bi = w->val[i] * log(1-si);
            g2sum += bi;
            gretl_matrix_set(G, i, 1, bi/wsum);
        }
        /* loop 3: form second component and subtract */
        for (i=0; i<p; i++) {
            ai = gretl_matrix_get(G, i, 0);
            ai -= w->val[i] * g1sum/ws2;
            gretl_matrix_set(G, i, 0, ai);
            bi = gretl_matrix_get(G, i, 1);
            bi -= w->val[i] * g2sum/ws2;
            gretl_matrix_set(G, i, 1, bi);
        }
        if (k == 3) {
            /* not zero-terminated */
            double c3 = theta[2];
            double m3 = 1 / (1 + p * c3);
            double g;

            for (i=0; i<2*p; i++) {
                /* scale the first two columns */
                G->val[i] *= m3;
            }
            for (i=0; i<p; i++) {
                /* compute the third-col derivative */
                g = 1 - p * w->val[i]/wsum;
                gretl_matrix_set(G, i, 2, m3 * m3 * g);
            }
        }
    } else if (method == MIDAS_ALMONP) {
        /* straight Almon polynomial */
        for (i=0; i<p; i++) {
            gretl_matrix_set(G, i, 0, 1.0);
            for (j=1; j<k; j++) {
                gretl_matrix_set(G, i, j, pow(i + 1.0, j));
            }
        }
    }

 range_error:

#if HAVE_GMP
    if (*err == E_NAN || mgrad_zero(G)) {
        /* attempt a fix-up using multiple precision */
        int save_errno = errno;

        errno = 0;
        *err = try_mp_midas_grad(theta, G, method);
        if (*err) {
            errno = save_errno;
        }
    }
#endif

    gretl_matrix_free(w);

    if (errno) {
        gretl_errmsg_sprintf(_("Failed to calculate MIDAS gradient: %s"),
                             gretl_strerror(errno));
        if (*err == 0) {
            *err = E_INVARG;
        }
        gretl_matrix_free(G);
        G = NULL;
        errno = 0;
    }

#if MDEBUG
    gretl_matrix_print(G, "G: midas gradient");
#endif

    return G;
}

/* Retrieve from a midasreg $model bundle all the info
   needed to construct the array of per-lag multipliers
   and their standard errors.
*/

static int process_midas_bundle (gretl_bundle *mb, int idx,
                                 gretl_matrix **ptheta,
                                 gretl_matrix **pV,
                                 int *pmtype, int *ph,
                                 int *pminlag,
                                 const char **plname)
{
    gretl_array *mt;
    gretl_matrix *b, *vcv;
    gretl_matrix *theta = NULL;
    gretl_matrix *V = NULL;
    int *xlist = NULL;
    int mtype, np;
    int i0, i1;
    int i, j, k, nm;
    int err = 0;

    mt = gretl_bundle_get_array(mb, "midas_info", &err);
    if (!err) {
        xlist = gretl_bundle_get_list(mb, "xlist", &err);
    }
    if (!err) {
        b = gretl_bundle_get_matrix(mb, "coeff", &err);
    }
    if (!err) {
        vcv = gretl_bundle_get_matrix(mb, "vcv", &err);
    }
    if (err) {
        return err;
    }

    i0 = xlist[0];
    nm = gretl_array_get_length(mt);
    if (idx < 1 || idx > nm) {
        gretl_errmsg_set(_("Invalid MIDAS term index"));
        err = E_DATA;
    } else {
        idx--; /* convert index to 0-based */
    }

    for (i=0; i<=idx && !err; i++) {
        gretl_bundle *mti;
        int minlag, maxlag;
        int nl, leader = 1;

        mti = gretl_array_get_data(mt, i);
        minlag = gretl_bundle_get_int(mti, "minlag", &err);
        maxlag = gretl_bundle_get_int(mti, "maxlag", &err);
        mtype = gretl_bundle_get_int(mti, "type", &err);
        np = gretl_bundle_get_int(mti, "nparm", &err);
        if (err) {
            break;
        }

        if (mtype == MIDAS_U || mtype == MIDAS_ALMONP) {
            leader = 0;
        }
        nl = maxlag - minlag + 1;
        i1 = i0 + np - 1 + leader;
#if 1
        fprintf(stderr, "midas term %d, type %d, np %d, nlags %d\n",
                i, mtype, np, nl);
#endif
        if (i == idx) {
            double vjk;
            int r = i1 - i0 + 1;
            int jj;

            theta = gretl_matrix_alloc(r, 1);
            V = gretl_matrix_alloc(r, r);
            for (j=0; j<r; j++) {
                jj = j + i0;
                theta->val[j] = b->val[jj];
                for (k=0; k<r; k++) {
                    vjk = gretl_matrix_get(vcv, jj, k+i0);
                    gretl_matrix_set(V, j, k, vjk);
                }
            }
            *pmtype = mtype;
            *ph = nl;
            *pminlag = minlag;
            *plname = gretl_bundle_get_string(mti, "lname", NULL);
        }
        i0 = i1 + 1;
    }

    if (err) {
        gretl_matrix_free(theta);
        gretl_matrix_free(V);
    } else {
        *ptheta = theta;
        *pV = V;
    }

    return err;
}

gretl_matrix *midas_multipliers (gretl_bundle *mb, int cumulate,
                                 int idx, int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *theta = NULL;
    gretl_matrix *mult = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *J = NULL;
    gretl_matrix *vcv = NULL;
    const char *lname = NULL;
    int minlag = 0;
    int mtype = 0, h = 0;
    int i, k;

    *err = process_midas_bundle(mb, idx, &theta, &V, &mtype,
                                &h, &minlag, &lname);
    if (*err) {
        gretl_errmsg_set(_("Not a valid midasreg bundle"));
        return NULL;
    }

    k = theta->rows;

    if (mtype == MIDAS_U) {
        mult = gretl_matrix_copy(theta);
        J = gretl_zero_matrix_new(h, k);
        gretl_matrix_inscribe_I(J, 0, 0, MIN(k,h));
    } else if (mtype == MIDAS_ALMONP) {
        mult = midas_weights(h, theta, mtype, err);
        J = midas_gradient(h, theta, mtype, err);
    } else {
        gretl_matrix *tmp1, *tmp2;
        double mg, th0 = theta->val[0];
        int j;

        tmp1 = gretl_matrix_alloc(k-1, 1);
        for (i=1; i<k; i++) {
            tmp1->val[i-1] = theta->val[i];
        }
        mult = midas_weights(h, tmp1, mtype, err);
        tmp2 = midas_gradient(h, tmp1, mtype, err);
        J = gretl_matrix_alloc(h, k);
        for (i=0; i<h; i++) {
            gretl_matrix_set(J, i, 0, mult->val[i]);
            for (j=1; j<k; j++) {
                mg = gretl_matrix_get(tmp2, i, j-1);
                gretl_matrix_set(J, i, j, th0*mg);
            }
        }
        gretl_matrix_multiply_by_scalar(mult, th0);
        gretl_matrix_free(tmp1);
        gretl_matrix_free(tmp2);
    }

    if (*err || mult == NULL || J == NULL) {
        if (!*err) {
            *err = E_ALLOC;
        }
        return NULL;
    }

    vcv = gretl_matrix_alloc(J->rows, J->rows);
    if (vcv == NULL) {
        *err = E_ALLOC;
    } else {
        *err = gretl_matrix_qform(J, GRETL_MOD_NONE,
                                  V, vcv, GRETL_MOD_NONE);
    }

    if (!*err && cumulate) {
        gretl_matrix *tmp;

        for (i=1; i<h; i++) {
            mult->val[i] += mult->val[i-1];
        }
        gretl_matrix_free(J);
        J = gretl_unit_matrix_new(h, h);
        gretl_matrix_set_triangle(J, NULL, 0, 1);
        tmp = gretl_matrix_alloc(h, h);
        *err = gretl_matrix_qform(J, GRETL_MOD_NONE,
                                  vcv, tmp, GRETL_MOD_NONE);
        gretl_matrix_free(vcv);
        vcv = tmp;
    }

    if (!*err) {
        double vi;

        ret = gretl_matrix_alloc(h, 2);
        for (i=0; i<h; i++) {
            gretl_matrix_set(ret, i, 0, mult->val[i]);
            vi = gretl_matrix_get(vcv, i, i);
            gretl_matrix_set(ret, i, 1, sqrt(vi));
        }
    }

    if (ret != NULL) {
        char **S = strings_array_new(h);
        char rname[32];

        if (S != NULL) {
            for (i=0; i<h; i++) {
                if (lname != NULL) {
                    if (cumulate) {
                        sprintf(rname, "c%s_%d", lname, i + minlag);
                    } else {
                        sprintf(rname, "%s_%d", lname, i + minlag);
                    }
                } else {
                    sprintf(rname, "%s_%d", cumulate ? "cmult" : "mult",
                            i + minlag);
                }
                S[i] = gretl_strdup(rname);
            }
            gretl_matrix_set_rownames(ret, S);
        }
    }

    gretl_matrix_free(theta);
    gretl_matrix_free(V);
    gretl_matrix_free(mult);
    gretl_matrix_free(J);
    gretl_matrix_free(vcv);

    return ret;
}

int midas_linear_combo (double *y, const int *list,
                        const gretl_matrix *theta,
                        int method, const DATASET *dset)
{
    gretl_matrix *w = NULL;
    int m = list[0];
    int err = 0;

    w = midas_weights(m, theta, method, &err);

    if (BETA_METHOD(method) && !err && w == NULL) {
        int t;

        for (t=dset->t1; t<=dset->t2; t++) {
            y[t] = NADBL;
        }
        return 0;
    }

    if (!err) {
        err = list_linear_combo(y, list, w, dset);
    }

    gretl_matrix_free(w);

    return err;
}

int *vector_to_midas_list (const gretl_matrix *v,
                           int f_ratio,
                           const char *prefix,
                           DATASET *dset,
                           int *err)
{
    char vname[VNAMELEN];
    int *list = NULL;
    int origv = dset->v;
    int i;

    /* double-check to avoid crashing below */
    if (gretl_vector_get_length(v) != sample_size(dset) * f_ratio) {
        *err = E_DATA;
        return NULL;
    }

    /* check names for collisions first */
    for (i=0; i<f_ratio && !*err; i++) {
        sprintf(vname, "%s%d", prefix, i+1);
        if (current_series_index(dset, vname) >= 1 ||
            get_user_var_by_name(vname) != NULL) {
            gretl_errmsg_set(_("The constructed series names would "
                             "collide with those of existing objects"));
            *err = E_INVARG;
        }
    }

    if (!*err) {
        /* try adding the required number of series */
        *err = dataset_add_series(dset, f_ratio);
        if (!*err) {
            list = gretl_list_new(f_ratio);
            if (list == NULL) {
                *err = E_ALLOC;
            }
        }
    }

    if (!*err) {
        /* actually construct the series */
        char label[MAXLABEL];
        int pos, pos0 = f_ratio - 1;
        int t, k = origv;

        for (i=0; i<f_ratio; i++) {
            sprintf(dset->varname[k], "%s%d", prefix, f_ratio - i);
            sprintf(label, "%s in sub-period %d", prefix, f_ratio - i);
            series_record_label(dset, k, label);
            list[i+1] = k;
            k++;
        }

        for (t=dset->t1; t<=dset->t2; t++) {
            pos = pos0;
            for (k=origv; k<dset->v; k++) {
                dset->Z[k][t] = v->val[pos--];
            }
            pos0 += f_ratio;
        }

        gretl_list_set_midas(list, dset);
    }

    return list;
}

/* Imhof: draws on the RATS code in IMHOF.SRC from Estima, 2004.

   Imhof Procedure for computing P(u'Au < x) for a quadratic form in
   Normal(0,1) variables. This can be used for ratios of quadratic
   forms as well, since P((u'Au/u'Bu) < x) = P(u'(A-xB)u < 0).

   In the Durbin-Watson context, 'B' is the identity matrix and
   'x' is the Durbin-Watson statistic.

   References:

   Imhof, J.P (1961), Computing the Distribution of Quadratic Forms of
   Normal Variables, Biometrika 48, 419-426.

   Abrahamse, A.P.J and Koerts, J. (1969), On the Theory and
   Application of the General Linear Model, Rotterdam University
   Press.
*/

static double imhof_bound (const gretl_matrix *lambda, int *err)
{
    double e1 = 0.0001; /* Max truncation error due to finite
                           upper bound on domain */
    double e2 = 0.0001; /* Cutoff for deciding whether an
                           eigenvalue is effectively zero */
    double absl, bound;
    double nl = 0.0, sum = 0.0;
    int i;

    for (i=0; i<lambda->rows; i++) {
        absl = fabs(lambda->val[i]);
        if (absl > e2) {
            nl += 1.0;
            sum += log(absl);
        }
    }

    if (nl == 0.0) {
        fprintf(stderr, "imhof_bound: got no non-zero eigenvalues\n");
        *err = E_DATA;
        return NADBL;
    }

    /* The key factor in the integrand is the product of
       (1+(lambda(i)*x)^2)^(1/4) across i. Since, for those
       factors with small |lambda(i)|, this won't go to zero very
       quickly, we count only the terms on the bigger eigenvalues.
    */

    nl *= 0.5;
    sum = 0.5 * sum + log(M_PI * nl);
    bound = exp(-(sum + log(e1)) / nl);
    bound += 5.0 / nl;

    if (bound < 0) {
        fprintf(stderr, "imhof_bound: got negative result\n");
        *err = E_DATA;
        bound = NADBL;
    }

    return bound;
}

static double imhof_f (double u, const gretl_matrix *lambda, double arg)
{
    double ul, rho = 0.0;
    double theta = -u * arg;
    int i, k = lambda->rows;

    /* The value at zero isn't directly computable as
       it produces 0/0. The limit is computed below.
    */
    if (u == 0.0) {
        double lsum = 0;

        for (i=0; i<k; i++) {
            lsum += lambda->val[i];
        }
        return 0.5 * (-arg + lsum);
    }

    for (i=0; i<k; i++) {
        ul = u * lambda->val[i];
        theta += atan(ul);
        rho += log(1.0 + ul * ul);
    }

    return sin(0.5 * theta) / (u * exp(0.25 * rho));
}

#define gridlimit 2048

/*
  Adaptation of Abrahamse and Koert's Pascal code.  Evaluates the
  integral by Simpson's rule with grid size halving each time until
  the change from moving to a tighter grid is negligible. By halving
  the grid, only the odd terms need to be computed, as the old ones
  are already in the sum. Points entering get x 4 weights, which then
  get reduced to x 2 on the next iteration.
*/

static double imhof_integral (double arg, const gretl_matrix *lambda,
                              double bound, int *err)
{
    double e3 = 0.0001;
    double base, step, sum1;
    double int0 = 0.0, int1 = 0.0;
    double eps4 = 3.0 * M_PI * e3;
    double sum4 = 0.0;
    double ret = NADBL;
    int j, n = 2;

    base = imhof_f(0, lambda, arg);
    base += imhof_f(bound, lambda, arg);

    while (n < gridlimit) {
        step = bound / n;
        sum1 = base + sum4 * 2.0;
        base = sum1;
        sum4 = 0.0;
        for (j=1; j<=n; j+=2) {
            sum4 += imhof_f(j * step, lambda, arg);
        }
        int1 = (sum1 + 4 * sum4) * step;
        if (n > 8 && fabs(int1 - int0) < eps4) {
            break;
        }
        int0 = int1;
        n *= 2;
    }

    if (n > gridlimit) {
        fprintf(stderr, "n = %d, Imhof integral failed to converge\n", n);
        *err = E_NOCONV;
    } else {
        ret = 0.5 - int1 / (3.0 * M_PI);
        if (ret < 0) {
            fprintf(stderr, "n = %d, Imhof integral gave negative value %g\n", n, ret);
        }
    }

    return ret;
}

/* Implements the "imhof" function in genr: computes the probability
   P(u'Au < arg) for a quadratic form in Normal(0,1) variables.  The
   argument @m may be either the square matrix A or a vector
   containing the precomputed eigenvalues of A.
*/

double imhof (const gretl_matrix *m, double arg, int *err)
{
    gretl_matrix *lambda = NULL;
    double bound, ret = NADBL;
    int free_lambda = 0;

    errno = 0;

    if (m->cols == 1 || m->rows == 1) {
        /* we'll assume m is a vector of eigenvalues */
        lambda = (gretl_matrix *) m;
    } else if (m->rows == m->cols) {
        /* we'll assume m is the 'A' matrix */
        lambda = gretl_general_matrix_eigenvals(m, err);
        free_lambda = 1;
    } else {
        /* huh? */
        *err = E_INVARG;
    }

    if (!*err) {
        bound = imhof_bound(lambda, err);
    }

    if (!*err) {
        ret = imhof_integral(arg, lambda, bound, err);
    }

    if (errno != 0) {
        fprintf(stderr, "imhof: %s\n", gretl_strerror(errno));
        if (!*err) {
            *err = E_NOCONV;
        }
        ret = NADBL;
        errno = 0;
    } else if (!*err && ret < 0 && ret > -1.0e-14) {
        /* just approximation error? */
        ret = 0;
    }

    if (free_lambda) {
        gretl_matrix_free(lambda);
    }

    return ret;
}

/* Implements the $dwpval accessor: given the residual vector
   @u and the matrix of regressors, @X, calculates the Durbin-Watson
   statistic then finds its p-value via the Imhof/Koerts/Abrahamse
   procedure. If @pDW is non-NULL, the Durbin-Watson statstic is
   written to that location.
*/

double dw_pval (const gretl_matrix *u, const gretl_matrix *X,
                double *pDW, int *perr)
{
    gretl_matrix *M = NULL;
    gretl_matrix *A = NULL;
    gretl_matrix *MA = NULL;
    gretl_matrix *XX = NULL;
    gretl_matrix *E = NULL;
    double uu, DW = 0;
    double pv = NADBL;
    int k = X->cols;
    int n = X->rows;
    int i, err = 0;

    M = gretl_identity_matrix_new(n);
    A = gretl_DW_matrix_new(n);
    MA = gretl_matrix_alloc(n, n);
    XX = gretl_matrix_alloc(k, k);

    if (M == NULL || A == NULL || MA == NULL || XX == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
                              X, GRETL_MOD_NONE,
                              XX, GRETL_MOD_NONE);
    err = gretl_invert_symmetric_matrix(XX);

    if (!err) {
        /* M = I - X(X'X)^{-1}X' */
        err = gretl_matrix_qform(X, GRETL_MOD_NONE,
                                 XX, M, GRETL_MOD_DECREMENT);
    }
    if (!err) {
        err = gretl_matrix_multiply(M, A, MA);
    }
    if (!err) {
        uu = gretl_matrix_dot_product(u, GRETL_MOD_TRANSPOSE,
                                      u, GRETL_MOD_NONE,
                                      &err);
    }
    if (!err) {
        DW = gretl_scalar_qform(u, A, &err);
    }
    if (!err) {
        DW /= uu;
        E = gretl_general_matrix_eigenvals(MA, &err);
#if 0
        fprintf(stderr, "DW = %g\n", DW);
        gretl_matrix_print(E, "eigenvals of M*A");
#endif
    }

    if (!err) {
        k = n - k;
        for (i=0; i<k; i++) {
            E->val[i] -= DW;
        }
        gretl_matrix_reuse(E, k, 1);
        pv = imhof(E, 0.0, &err);
        if (!err && pDW != NULL) {
            *pDW = DW;
        }
    }

 bailout:

    gretl_matrix_free(M);
    gretl_matrix_free(A);
    gretl_matrix_free(MA);
    gretl_matrix_free(XX);
    gretl_matrix_free(E);

    *perr = err;

    return pv;
}

/* create a matrix containing ACF and PACF values for each
   column of the input matrix, @m, with lag order @p.
*/

gretl_matrix *multi_acf (const gretl_matrix *m,
                         const int *list,
                         const DATASET *dset,
                         int p, int *err)
{
    gretl_matrix *a, *A = NULL;
    const double *x;
    double xa;
    int nv, T, acol, pcol;
    int i, j;

    if (list == NULL && gretl_is_null_matrix(m)) {
        *err = E_DATA;
        return NULL;
    }

    if (m != NULL) {
        nv = m->cols;
    } else {
        nv = list[0];
    }

    A = gretl_matrix_alloc(p, 2 * nv);
    if (A == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    if (m != NULL) {
        x = m->val;
        T = m->rows;
    } else {
        x = dset->Z[list[1]] + dset->t1;
        T = sample_size(dset);
    }

    acol = 0;
    pcol = nv;

    for (j=0; j<nv; j++) {
        /* get ACF/PACF for column/series */
        a = acf_matrix(x, p, NULL, T, err);
        if (*err) {
            gretl_matrix_free(a);
            gretl_matrix_free(A);
            return NULL;
        }

        /* transcribe into A matrix and free */
        for (i=0; i<p; i++) {
            xa = gretl_matrix_get(a, i, 0);
            gretl_matrix_set(A, i, acol, xa);
            xa = gretl_matrix_get(a, i, 1);
            gretl_matrix_set(A, i, pcol, xa);
        }
        acol++;
        pcol++;

        gretl_matrix_free(a);

        /* move to next data read position */
        if (j < nv - 1) {
            if (m != NULL) {
                x += m->rows;
            } else {
                x = dset->Z[list[j+2]] + dset->t1;
            }
        }
    }

    return A;
}

gretl_matrix *multi_xcf (const void *px, int xtype,
                         const void *py, int ytype,
                         const DATASET *dset,
                         int p, int *err)
{
    const int *xlist = NULL;
    const gretl_matrix *Xmat = NULL;
    const double *xvec = NULL;
    const double *yvec = NULL;
    gretl_matrix *xj, *XCF = NULL;
    int T = sample_size(dset);
    int np = 2 * p + 1;
    int Ty, nx = 1;
    int i, j;

    if (xtype == LIST) {
        xlist = px;
        nx = xlist[0];
        if (nx < 1) {
            *err = E_DATA;
            return NULL;
        }
        xvec = dset->Z[xlist[1]] + dset->t1;
    } else if (xtype == MAT) {
        Xmat = px;
        if (gretl_is_null_matrix(Xmat)) {
            *err = E_DATA;
            return NULL;
        }
        nx = Xmat->cols;
        T = Xmat->rows;
        xvec = Xmat->val;
    } else {
        /* VEC: note: px is of type *void */
        xvec = px;
        xvec += dset->t1;
    }

    if (ytype == MAT) {
        const gretl_matrix *ymat = py;

        if (ymat->cols != 1) {
            *err = E_NONCONF;
            return NULL;
        }
        yvec = ymat->val;
        Ty = ymat->rows;
    } else {
        yvec = (const double *) py + dset->t1;
        Ty = sample_size(dset);
    }

    if (Ty != T) {
        *err = E_NONCONF;
        return NULL;
    }

    if (nx > 1) {
        XCF = gretl_matrix_alloc(np, nx);
        if (XCF == NULL) {
            *err = E_ALLOC;
            return NULL;
        }
    }

    for (j=0; j<nx; j++) {
        /* get XCF for left-hand column/series and y */
        xj = xcf_vec(xvec, yvec, p, NULL, T, err);
        if (*err) {
            gretl_matrix_free(XCF);
            return NULL;
        }

        if (nx == 1) {
            XCF = xj;
            break;
        }

        /* transcribe into big XCF matrix and free */
        for (i=0; i<np; i++) {
            gretl_matrix_set(XCF, i, j, xj->val[i]);
        }
        gretl_matrix_free(xj);

        /* move to next data read position */
        if (j < nx - 1) {
            if (Xmat != NULL) {
                xvec += Xmat->rows;
            } else {
                xvec = dset->Z[xlist[j+2]] + dset->t1;
            }
        }
    }

    return XCF;
}

static int theil_decomp (double *m, double MSE,
                         const double *y,
                         const double *f,
                         int T)
{
    double da, dp;
    double Abar, Pbar;
    double sa, sp, r;
    int t, err = 0;

    if (MSE <= 0.0) {
        m[0] = m[1] = m[2] = NADBL;
        return E_DATA;
    }

    Abar = Pbar = 0.0;

    for (t=0; t<T; t++) {
        Abar += y[t];
        Pbar += f[t];
    }

    Abar /= T;
    Pbar /= T;

    sa = sp = r = 0.0;

    for (t=0; t<T; t++) {
        da = y[t] - Abar;
        sa += da * da;
        dp = f[t] - Pbar;
        sp += dp * dp;
        r += da * dp;
    }

    sa = sqrt(sa / T);
    sp = sqrt(sp / T);
    r /= T * sa * sp;

    if (sa == 0.0 || sp == 0.0) {
        err = E_DATA;
        m[0] = m[1] = m[2] = NADBL;
    } else {
        m[0] = (Abar - Pbar) * (Abar - Pbar) / MSE; /* U^M */
        m[1] = (sp - r * sa) * (sp - r * sa) / MSE; /* U^R */
        m[2] = (1.0 - r * r) * sa * sa / MSE;       /* U^D */

        if (m[2] > 0.99999999999999) {
            /* U^M and U^R are just machine noise? */
            m[2] = 1.0;
            m[0] = m[1] = 0.0;
        }
    }

    return err;
}

/* OPT_O allows embedded missing values (which will be skipped);
   OPT_D requests Theil decomposition; OPT_T means time series.
*/

static int fill_fcstats_column (gretl_matrix *m,
                                const double *y,
                                const double *f,
                                int T,
                                gretlopt opt,
                                int col)
{
    double ME, MSE, MAE, MPE, MAPE, U;
    double Unum, Uden1, Uden2;
    double fe, u[2];
    int do_theil = 0;
    int do_U2 = (opt & OPT_T);
    int ok_T = T;
    int t, err = 0;

    ME = MSE = MAE = MPE = MAPE = U = 0.0;
    Unum = Uden1 = Uden2 = 0.0;
    u[0] = u[1] = 0.0;

    for (t=0; t<T; t++) {
        if (na(y[t]) || na(f[t])) {
            if (opt & OPT_O) {
                ok_T--;
                continue;
            } else {
                err = E_MISSDATA;
                break;
            }
        }
        fe = y[t] - f[t];
        ME += fe;
        MSE += fe * fe;
        MAE += fabs(fe);
        if (floateq(fe, 0)) {
            ; /* OK, zero contribution to MPE, MAPE */
        } else if (y[t] == 0.0) {
            /* can't calculate percentage */
            MPE = MAPE = NADBL;
        } else {
            MPE += 100 * fe / y[t];
            MAPE += 100 * fabs(fe / y[t]);
        }
        if (do_U2) {
            /* let U = U2 */
            if (t < T-1 && !na(U)) {
                if (na(f[t+1]) || na(y[t+1])) {
                    U = NADBL;
                } else {
                    fe = f[t+1] - y[t+1];
                    if (floatneq(fe, 0)) {
                        if (y[t] == 0.0) {
                            U = NADBL;
                        } else {
                            fe /= y[t];
                            u[0] += fe * fe;
                        }
                    }
                    fe = y[t+1] - y[t];
                    if (floatneq(fe, 0)) {
                        if (y[t] == 0.0) {
                            U = NADBL;
                        } else {
                            fe /= y[t];
                            u[1] += fe * fe;
                        }
                    }
                }
            }
        } else {
            /* let U = U1 */
            Unum += fe * fe;
            Uden1 += y[t] * y[t];
            Uden2 += f[t] * f[t];
        }
    }

    if (!err) {
        if (ok_T == 0) {
            err = E_MISSDATA;
        } else if (ok_T < T) {
            T = ok_T;
        } else if (opt & OPT_D) {
            do_theil = 1;
        }
    }

    if (!err) {
        fnpkg *pkg = get_active_function_package(0);
        const char *pkgname = NULL;
        int show_MSE = 0;

        if (pkg != NULL) {
            pkgname = function_package_get_name(pkg);
        }
        if (pkgname != NULL && !strcmp(pkgname, "fcModels") &&
            function_package_get_version(pkg) <= 1.1) {
            show_MSE = 1;
        }

        ME /= T;
        MSE /= T;
        MAE /= T;
        if (!isnan(MPE)) {
            MPE /= T;
        }
        if (!isnan(MAPE)) {
            MAPE /= T;
        }
        if (do_U2) {
            if (!isnan(U) && u[1] > 0.0) {
                U = sqrt(u[0] / T) / sqrt(u[1] / T);
            }
        } else {
            Unum = sqrt(Unum/T);
            Uden1 = sqrt(Uden1/T);
            Uden2 = sqrt(Uden2/T);
            U = Unum / (Uden1+ Uden2);
        }
        gretl_matrix_set(m, 0, col, ME);
        if (show_MSE) {
            gretl_matrix_set(m, 1, col, MSE);
        } else {
            gretl_matrix_set(m, 1, col, sqrt(MSE));
        }
        gretl_matrix_set(m, 2, col, MAE);
        gretl_matrix_set(m, 3, col, MPE);
        gretl_matrix_set(m, 4, col, MAPE);
        gretl_matrix_set(m, 5, col, U);

        if (do_theil) {
            double *targ = m->val + col * m->rows + 6;

            theil_decomp(targ, MSE, y, f, T);
        }
    }

    return err;
}

static void add_fcstats_rownames (gretl_matrix *m,
                                  gretlopt opt)
{
    const char *S[] = {
        "ME", "RMSE", "MAE", "MPE", "MAPE",
        "U1", "UM", "UR", "UD"
    };
    int i, ns = m->rows;
    char **rownames = strings_array_new(ns);

    if (rownames != NULL) {
        for (i=0; i<ns; i++) {
            if (i == 5 && (opt & OPT_T)) {
                rownames[5] = gretl_strdup("U2");
            } else {
                rownames[i] = gretl_strdup(S[i]);
            }
        }
        gretl_matrix_set_rownames(m, rownames);
    }
}

static int fcstats_sample_check (const double *y,
                                 const double *f,
                                 int *pt1,
                                 int *pt2,
                                 int *nmiss)
{
    int t1 = *pt1;
    int t2 = *pt2;
    int t, err = 0;

    for (t=t1; t<=t2; t++) {
        if (na(y[t]) || na(f[t])) {
            t1++;
        } else {
            break;
        }
    }

    for (t=t2; t>=t1; t--) {
        if (na(y[t]) || na(f[t])) {
            t2--;
        } else {
            break;
        }
    }

    if (t2 - t1 + 1 < 1) {
        err = E_MISSDATA;
    } else if (nmiss != NULL) {
        /* allow internal NAs */
        for (t=t1; t<=t2; t++) {
            if (na(y[t]) || na(f[t])) {
                *nmiss += 1;
            }
        }
    } else {
        /* flag error on internal NAs */
        for (t=t1; t<=t2; t++) {
            if (na(y[t]) || na(f[t])) {
                err = E_MISSDATA;
                break;
            }
        }
    }

    *pt1 = t1;
    *pt2 = t2;

    return err;
}

/*
   Forecast evaluation statistics: @y is the data series, @f the
   forecast.

   cf. http://www.economicsnetwork.ac.uk/showcase/cook_forecast
   by Steven Cook of Swansea University <s.cook@Swansea.ac.uk>

   OPT_D indicates that we should include the Theil decomposition.
   OPT_O allows missing values (which will be skipped).
   OPT_T indicates that the data are time series.
*/

gretl_matrix *forecast_stats (const double *y, const double *f,
                              int t1, int t2, int *n_used,
                              gretlopt opt, int *err)
{
    gretl_matrix *m = NULL;
    int rows = 6;
    int nmiss = 0;

    if (opt & OPT_O) {
        *err = fcstats_sample_check(y, f, &t1, &t2, &nmiss);
    } else {
        *err = fcstats_sample_check(y, f, &t1, &t2, NULL);
    }
    if (*err) {
        return NULL;
    }

    if (opt & OPT_D) {
        if (nmiss == 0) {
            /* extra rows for Theil decomp */
            rows = 9;
        } else {
            /* scrub the Theil decomp option */
            opt &= ~OPT_D;
        }
    }

    m = gretl_column_vector_alloc(rows);

    if (m == NULL) {
        *err = E_ALLOC;
    } else {
        int T = t2 - t1 + 1;

        *err = fill_fcstats_column(m, y + t1, f + t1, T, opt, 0);
    }

    if (*err) {
        gretl_matrix_free(m);
        m = NULL;
    } else {
        if (n_used != NULL) {
            *n_used = t2 - t1 + 1 - nmiss;
        }
        add_fcstats_rownames(m, opt);
    }

    return m;
}

gretl_matrix *matrix_fc_stats (const double *y,
                               const gretl_matrix *F,
                               gretlopt opt,
                               int *err)
{
    int ns = (opt & OPT_D)? 9 : 6;
    gretl_matrix *m = NULL;
    const double *f;
    int j, T = F->rows;

    m = gretl_matrix_alloc(ns, F->cols);
    if (m == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    for (j=0; j<F->cols; j++) {
        f = F->val + T * j;
        *err = fill_fcstats_column(m, y, f, T, opt, j);
        if (*err) {
            break;
        }
    }

    if (*err) {
        gretl_matrix_free(m);
        m = NULL;
    } else {
        add_fcstats_rownames(m, opt);
    }

    return m;
}

/* Remove "repeats" from Nelson-Aalen or Kaplan-Meier output matrix:
   these will be identified by all-zero rows, since we defer entering
   the calculated values until we reach the end of a sequence of
   repeated duration values.
*/

static int durations_squeeze (gretl_matrix **pA, gretl_matrix **pTmp,
                              int n, int repeats)
{
    gretl_matrix *T, *A = *pA;
    gretl_matrix *Tmp = *pTmp;
    double aij;
    int m = n - repeats;
    int i, j, err;

    err = gretl_matrix_realloc(Tmp, m, 3);
    if (err) {
        return err;
    }

    j = 0;
    for (i=0; i<n; i++) {
        aij = gretl_matrix_get(A, i, 0);
        if (aij > 0) {
            /* got some actual values: transcribe */
            gretl_matrix_set(Tmp, j, 0, aij);
            aij = gretl_matrix_get(A, i, 1);
            gretl_matrix_set(Tmp, j, 1, aij);
            aij = gretl_matrix_get(A, i, 2);
            gretl_matrix_set(Tmp, j, 2, aij);
            j++;
        }
    }

    /* swap matrix pointers */
    T = *pA;
    *pA = Tmp;
    *pTmp = T;

    return err;
}

#define cmissing(c,t) (c != NULL && na(c[t]))

/*
  Given a series of duration values @y and (possibly) a series of
  censoring values @cens (zero for uncensored, non-zero for censored),
  computes either the Nelson-Aalen estimator of the cumulative hazard
  function OR (given OPT_K) the Kaplan-Meier estimate of the survival
  function.

  Th return value is a matrix with the (unique) duration value in the
  first column, the estimator in the second and its standard error
  in the third. Any repeated duration values in the input data are
  squeezed out.

  On Nelson-Aalen, the piece by O. Borgan at
  http://www.med.mcgill.ca/epidemiology/hanley/c609/material/NelsonAalenEstimator.pdf
  is quite helpful, and on Kaplan-Meier see
  www.public.iastate.edu/~kkoehler/stat565/kaplanmeier.4page.pdf
*/

gretl_matrix *duration_func (const double *y, const double *cens,
                             int t1, int t2, gretlopt opt,
                             int *err)
{
    gretl_matrix *Tmp = NULL;
    gretl_matrix *M = NULL;
    gretl_matrix *A = NULL;
    double yi, ai, aibak;
    double vaibak, vai = 0;
    int kmeier = (opt & OPT_K);
    int nmax = t2 - t1 + 1;
    int i, t, n = nmax;
    int repeats = 0;
    int d, r, ni, ibak;

    for (t=t1; t<=t2; t++) {
        if (na(y[t]) || cmissing(cens, t)) {
            n--;
        }
    }

    if (n < 1) {
        *err = E_TOOFEW;
        return NULL;
    }

    Tmp = gretl_zero_matrix_new(n, 2);
    A   = gretl_zero_matrix_new(n, 3);
    if (Tmp == NULL || A == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    i = 0;
    for (t=t1; t<=t2; t++) {
        if (n < nmax && (na(y[t]) || cmissing(cens, t))) {
            continue;
        }
        gretl_matrix_set(Tmp, i, 0, y[t]);
        if (cens != NULL) {
            gretl_matrix_set(Tmp, i, 1, cens[t]);
        }
        i++;
    }

    M = gretl_matrix_sort_by_column(Tmp, 0, err);
    if (*err) {
        goto bailout;
    }

    r = n;
    ibak = -1;

    for (i=0; i<n; i++) {
        /* if not censored, increment d count */
        d = gretl_matrix_get(M, i, 1) == 0;
        /* start count of "drop-outs" */
        ni = 1;
        /* record current y value */
        yi = gretl_matrix_get(M, i, 0);
        while (i < n-1) {
            if (gretl_matrix_get(M, i+1, 0) != yi) {
                break;
            }
            i++;
            d += gretl_matrix_get(M, i, 1) == 0;
            ni++;
            repeats++;
        }
        if (kmeier) {
            /* survival */
            ai = (r - d) / (double) r;
            vai = d / (r * ((double) r - d));
        } else {
            /* hazard */
            ai = d / (double) r;
            vai = ((r - d) * d) / ((r - 1) * (double) r*r);
        }
#if 0
        fprintf(stderr, "i=%d, y=%g, d=%d, r=%d\n", i, yi, d, r);
#endif
        gretl_matrix_set(A, i, 0, yi);
        if (ibak >= 0) {
            aibak = gretl_matrix_get(A, ibak, 1);
            vaibak = gretl_matrix_get(A, ibak, 2);
            if (kmeier) {
                gretl_matrix_set(A, i, 1, aibak * ai);
                gretl_matrix_set(A, i, 2, vaibak + vai);
            } else {
                gretl_matrix_set(A, i, 1, aibak + ai);
                gretl_matrix_set(A, i, 2, vaibak + vai);
            }
        } else {
            gretl_matrix_set(A, i, 1, ai);
            gretl_matrix_set(A, i, 2, vai);
        }
        /* decrement number "at risk" */
        r -= ni;
        ibak = i;
    }

    if (kmeier) {
        /* compute variance using @vai cumulator */
        for (i=0; i<n; i++) {
            ai  = gretl_matrix_get(A, i, 1);
            vai = gretl_matrix_get(A, i, 2);
            gretl_matrix_set(A, i, 2, sqrt(ai * ai * vai));
        }
    } else {
        /* Nelson-Aalen: convert to std error */
        for (i=0; i<n; i++) {
            vai = gretl_matrix_get(A, i, 2);
            gretl_matrix_set(A, i, 2, sqrt(vai));
        }
    }

#if 0
    gretl_matrix_print(A, "A");
#endif

    if (repeats > 0) {
        *err = durations_squeeze(&A, &Tmp, n, repeats);
    }

 bailout:

    if (*err) {
        gretl_matrix_free(A);
        A = NULL;
    }
    gretl_matrix_free(Tmp);
    gretl_matrix_free(M);

    return A;
}

double gretl_round (double x)
{
    double fx = floor(x);

    if (x < 0) {
        return (x - fx <= 0.5)? fx : ceil(x);
    } else {
        return (x - fx < 0.5)? fx : ceil(x);
    }
}

double gretl_sgn (double x)
{
    if (x == 0) {
        return 0;
    } else if (x > 0) {
        return 1;
    } else if (x < 0) {
        return -1;
    } else {
        return NADBL;
    }
}

#define neg_x_real_v_err(t) (t == 'J' || t == 'I')

/* evaluates bessel function for scalar nu and x */

double gretl_bessel (char type, double v, double x, int *err)
{
    if (na(x) || na(v)) {
        return NADBL;
    }

    if (x < 0) {
        /* catch invalid cases for x < 0 */
        if (type == 'K') {
            *err = E_INVARG;
            return NADBL;
        } else if (v != floor(v) && (type == 'J' || type == 'I')) {
            *err = E_INVARG;
            return NADBL;
        }
    }

    switch (type) {
    case 'J':
        return cephes_bessel_Jv(v, x);
    case 'Y':
        return cephes_bessel_Yv(v, x);
    case 'I':
        if (v == 0) {
            return cephes_bessel_I0(x);
        } else if (v == 1) {
            return cephes_bessel_I1(x);
        } else if (v > 0) {
            return cephes_bessel_Iv(v, x);
        } else {
            /* cephes_bessel_Iv is not right for v < 0 */
            double b1 = netlib_bessel_K(-v, x, 1);
            double b2 = cephes_bessel_Iv(-v, x);

            return (2*b1*sin(-v*M_PI)) / M_PI + b2;
        }
        break;
    case 'K':
        /* bessel K is symmetric around v = 0 */
        v = fabs(v);
        if (v == 0) {
            return cephes_bessel_K0(x);
        } else if (v == 1) {
            return cephes_bessel_K1(x);
        } else if (v == floor(v) && v <= 30.0) {
            /* cephes doesn't do non-integer v, and also loses
               accuracy beyond |v| = 30
            */
            return cephes_bessel_Kn(v, x);
        } else {
            /* accurate but expensive */
            return netlib_bessel_K(v, x, 1);
        }
        break;
    default:
        /* unknown function type */
        return NADBL;
    }
}

/* Net Present Value: note that the first value is taken as
   occurring "now" and is not discounted */

double gretl_npv (int t1, int t2, const double *x, double r,
                  int pd, int *err)
{
    double d, PV = 0.0;
    int i, n = 0;

    if (pd != 1 && pd != 4 && pd != 12) {
        *err = E_PDWRONG;
        return NADBL;
    }

    if (pd == 1) {
        d = 1 + r;
    } else if (r < -1.0) {
        *err = E_NAN;
        return 0.0 / 0.0;
    } else {
        d = pow(1 + r, 1.0 / pd);
    }

    for (i=t1; i<=t2; i++) {
        if (!na(x[i])) {
            PV += x[i] / (pow(d, i-t1));
            n++;
        }
    }

    if (n == 0) {
        PV = NADBL;
    }

    return PV;
}

/* Internal Rate of Return */

double gretl_irr (const double *x, int n, int pd, int *err)
{
    double PV, r, r1 = 0.02, r0 = -0.02;
    int gotplus = 0, gotminus = 0;
    int i, m = n;

    for (i=0; i<n; i++) {
        if (na(x[i])) {
            m--;
        } else if (x[i] > 0) {
            gotplus = 1;
        } else if (x[i] < 0) {
            gotminus = 1;
        }
    }

    if (!gotplus && !gotminus) {
        /* null payment stream */
        return (m > 0)? 0 : NADBL;
    }

    if (gotplus && !gotminus) {
        /* payments all positive */
        return (x[0] > 0)? 0.0/0.0 : 1.0/0.0;
    } else if (gotminus && !gotplus) {
        /* payments all negative */
        return (x[0] < 0)? 0.0/0.0 : -1.0/0.0;
    }

    /* find (r0, r1) bracket for solution, if possible */

    while ((PV = gretl_npv(0, n-1, x, r0, pd, err)) < 0 && !*err) {
        if (r0 < -DBL_MAX / 2.0) {
            return -1.0/0.0;
        }
        r1 = r0;
        r0 *= 2.0;
    }

    while ((PV = gretl_npv(0, n-1, x, r1, pd, err)) > 0 && !*err) {
        if (r1 > DBL_MAX / 2.0) {
            return 1.0/0.0;
        }
        r0 = r1;
        r1 *= 2.0;
    }

#if 0
    fprintf(stderr, "initial bracket for r: %g to %g\n", r0, r1);
#endif

    r = r1;

    /* now do binary search */

    for (i=0; i<32 && !*err; i++) {
        if (floateq(PV, 0.0)) {
            break;
        }
        if (PV < 0) {
            /* r is too high */
            if (r < r1) {
                r1 = r;
            }
            r = (r + r0) / 2.0;
        } else {
            /* r too low */
            if (r > r0) {
                r0 = r;
            }
            r = (r + r1) / 2.0;
        }
        PV = gretl_npv(0, n-1, x, r, pd, err);
#if 0
        fprintf(stderr, "binary search: r = %.9g, PV = %g\n", r, PV);
#endif
    }

    if (*err) {
        r = NADBL;
    }

    return r;
}

double logistic_cdf (double x)
{
    double emx, ret;

    errno = 0;

    emx = exp(-x);
    if (errno == ERANGE) {
        ret = (x > 0)? 1.0 : 0.0;
        errno = 0;
    } else {
        ret = 1.0 / (1.0 + emx);
    }

    return ret;
}

int list_ok_dollar_vars (DATASET *dset, PRN *prn)
{
    int nm = 0;
    int i;

    pprintf(prn, "\n%s\n", _("model-related"));

    for (i=R_MAX+1; i<M_MAX; i++) {
        GretlType type = GRETL_TYPE_NONE;
        double x = NADBL;
        gretl_matrix *m = NULL;
        int *list = NULL;
        int err = 0;

        if (i < M_SCALAR_MAX) {
            x = saved_object_get_scalar(NULL, i, dset, &err);
            if (!na(x)) {
                type = GRETL_TYPE_DOUBLE;
            }
        } else if (i > M_SCALAR_MAX && i < M_SERIES_MAX) {
            type = GRETL_TYPE_SERIES;
            err = saved_object_get_series(NULL, NULL, i, dset);
            if (err) {
                if (i == M_UHAT || i == M_YHAT || i == M_SIGMA) {
                    /* maybe the result is a matrix? */
                    type = saved_object_get_data_type(NULL, i);
                    if (type == GRETL_TYPE_MATRIX) {
                        m = saved_object_get_matrix(NULL, i, &err);
                    }
                }
            }
        } else if (i > M_SERIES_MAX && i < M_MATRIX_MAX) {
            type = GRETL_TYPE_MATRIX;
            m = saved_object_get_matrix(NULL, i, &err);
        } else if (i > M_MATRIX_MAX && i < M_MBUILD_MAX) {
            type = GRETL_TYPE_MATRIX;
            m = saved_object_build_matrix(NULL, i,
                                          dset, &err);
        } else {
            type = GRETL_TYPE_LIST;
            list = saved_object_get_list(NULL, i, &err);
        }

        if (!err && type != GRETL_TYPE_NONE) {
            const char *typestr = gretl_type_get_name(type);

            if (!na(x)) {
                pprintf(prn, " %s (%s: %g)\n", mvarname(i), typestr, x);
            } else {
                pprintf(prn, " %s (%s)\n", mvarname(i), typestr);
            }
            gretl_matrix_free(m);
            free(list);
            nm++;
        }
    }

    if (nm == 0) {
        pprintf(prn, " %s\n", _("none"));
    }

    pprintf(prn, "\n%s\n", _("other"));

    for (i=1; i<R_SCALAR_MAX; i++) {
        double x = dvar_get_scalar(i, dset);

        if (!na(x)) {
            pprintf(prn, " %s (scalar: %g)\n", dvarname(i), x);
        }
    }

    pputc(prn, '\n');

    return 0;
}

static double nw_kernel (double x)
{
    /*
       eventually, a libset variable will be used to choose among
       various kernels; for now, the Normal density will have to do.
    */
    double ret, x2 = x*x;

#if 0
    ret = exp(-0.5*x2)/SQRT_2_PI;
#else
    ret = exp(-0.5*x2);
#endif

    return ret;
}

/**
 * nadaraya_watson:
 * @y: array with "dependent variable"
 * @x: array with "explanatory variable"
 * @h: double, bandwidth (may be negative; see below)
 * @dset: data set information.
 * @LOO: Boolean flag (see below)
 * @trim: trim parameter (see below)
 * @m: array to hold results
 *
 * Implements the Nadaraya-Watson nonparametric estimator for the
 * conditional mean of @y given @x via the formula
 *
 * \widehat{m}_h(x)=\frac{\sum_{i=1}^n K_h(x-X_i)
 * Y_i}{\sum_{i=1}^nK_h(x-X_i)}
 *
 * and computes it for all elements of @x. Note that, in principle,
 * the kernel K(X_i-X_j) must be computed for every combination of i
 * and j, but since the function K() is assumed to be symmetric, we
 * compute it once to save time.
 *
 * The scalar @h holds the kernel bandwidth; if @LOO is non-zero the
 * "leave-one-out" variant of the estimator (essentially a jackknife
 * estimator; see Pagan and Ullah, Nonparametric Econometrics, p. 119)
 * is computed.
 *
 * A rudimentary form of trimming is implemented: the kernel function
 * is set to 0 when the product of the @trim parameter times the
 * bandwidth @h exceeds |X_i - X_j|, so as to speed computation and
 * enhance numerical stability.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int nadaraya_watson (const double *y, const double *x, double h,
                     DATASET *dset, int LOO, double trim,
                     double *m)
{
    int t, s, err = 0;
    int t1 = dset->t1, t2 = dset->t2;
    double xt, xs, ys, yt, k;
    double *num, *den;
    int n = t2 + 1; /* really? */

    if (h < 0.0) {
        return E_DATA;
    } else if (h == 0.0) {
        /* automatic data-based bandwidth */
        const double *sampx = x + dset->t1;
        int n = sample_size(dset);

        h = kernel_bandwidth(sampx, n);
    }

    num = malloc(2 * n * sizeof *num);
    if (num == NULL) {
        return E_ALLOC;
    }

    den = num + n;
    trim *= h;

    /* here we initialize numerator and denominator; we use the
       "diagonal" in the standard case and 0 in the leave-one-out
       case
    */

    if (LOO) {
        for (t=t1; t<=t2; t++) {
            num[t] = den[t] = 0.0;
        }
    } else {
        k = nw_kernel(0);
        for (t=t1; t<=t2; t++) {
            if (!na(y[t])) {
                num[t] = k * y[t];
                den[t] = k;
            } else {
                num[t] = 0;
                den[t] = 0;
            }
        }
    }

    for (t=t1; t<=t2; t++) {
        xt = x[t];
        if (!na(xt)) {
            yt = y[t];
            for (s=t+1; s<=t2; s++) {
                xs = x[s];
                if (!na(xs) && fabs(xs-xt) < trim) {
                    k = nw_kernel((xt - xs)/h);
                    if (!na(yt)) {
                        num[s] += k * yt;
                        den[s] += k;
                    }
                    ys = y[s];
                    if (!na(ys)) {
                        num[t] += k * ys;
                        den[t] += k;
                    }
                }
            }
            m[t] = num[t] / den[t];
        } else {
            m[t] = NADBL;
        }
    }

    free(num);

    return err;
}

static int xy_get_sample (const double *y, const double *x,
                          int *t1, int *t2, int *n)
{
    int t, nxy, err = 0;

    for (t=*t1; t<=*t2; t++) {
        if (na(x[t])) {
            *t1 += 1;
        } else {
            break;
        }
    }

    for (t=*t2; t>=*t1; t--) {
        if (na(x[t])) {
            *t2 -= 1;
        } else {
            break;
        }
    }

    /* nxy = the number of points where we have valid values
       for both x and y; n = the number of valid x-values
    */

    nxy = *n = 0;

    for (t=*t1; t<=*t2; t++) {
        if (!na(x[t])) {
            *n += 1;
            if (!na(y[t])) {
                nxy++;
            }
        }
    }

    if (nxy < 16) {
        err = E_TOOFEW;
    }

    return err;
}

static int get_dataset_t (const double *x, int pos, int t1)
{
    int t, k = 0;

    for (t=t1; k<=pos; t++) {
        if (!na(x[t])) {
            if (k == pos) {
                return t;
            }
            k++;
        }
    }

    return 0;
}

/* note: the following requires that @y and @x are truly dataset
   series; vectors cannot be accepted given the dependence on
   dset->t1 and dset->t2 for the range of data to be used
*/

int gretl_loess (const double *y, const double *x, int poly_order,
                 double bandwidth, gretlopt opt, DATASET *dset,
                 double *m)
{
    gretl_matrix *my, *mx;
    gretl_matrix *yh = NULL;
    int *s_order = NULL;
    int t1 = dset->t1;
    int t2 = dset->t2;
    int s, t, n;
    int err = 0;

    if (poly_order < 0 || poly_order > 2 ||
        bandwidth <= 0 || bandwidth >= 1) {
        return E_DATA;
    }

    err = xy_get_sample(y, x, &t1, &t2, &n);
    if (err) {
        return err;
    }

    /* note: n holds the number of non-missing observations
       on x; the associated y-value may or may not be
       missing
    */

    my = gretl_column_vector_alloc(n);
    mx = gretl_column_vector_alloc(n);

    if (my == NULL || mx == NULL) {
        err = E_ALLOC;
    } else {
        s = 0;
        for (t=t1; t<=t2; t++) {
            if (!na(x[t])) {
                my->val[s] = y[t];
                mx->val[s] = x[t];
                s++;
            }
        }
        /* sort the points by the value of x */
        err = sort_pairs_by_x(mx, my, &s_order, NULL);
    }

    if (!err) {
        yh = loess_fit(mx, my, poly_order, bandwidth, opt, &err);
    }

    if (!err) {
        /* put the yh values into dataset order */
        for (s=0; s<n; s++) {
            t = get_dataset_t(x, s_order[s], t1);
            m[t] = yh->val[s];
        }
    }

    gretl_matrix_free(my);
    gretl_matrix_free(mx);
    gretl_matrix_free(yh);
    free(s_order);

    return err;
}

double series_get_nobs (int t1, int t2, const double *x)
{
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) n++;
    }

    return n;
}

double series_sum_all (int t1, int t2, const double *x)
{
    double xsum = 0.0;
    int t;

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            xsum = NADBL;
            break;
        } else {
            xsum += x[t];
        }
    }

    return xsum;
}

static gretl_matrix *delete_null_cases (gretl_matrix *m,
                                        int keeprows,
                                        int *err)
{
    gretl_matrix *ret;

    ret = gretl_matrix_alloc(keeprows, m->cols);
    if (ret == NULL) {
        *err = E_ALLOC;
    } else {
        double x;
        int i, j;

        for (i=0; i<keeprows; i++) {
            for (j=0; j<m->cols; j++) {
                x = gretl_matrix_get(m, i, j);
                gretl_matrix_set(ret, i, j, x);
            }
        }
    }

    gretl_matrix_free(m);

    return ret;
}

static gretl_matrix *real_aggregate_by (const double *x,
                                        const double *y,
                                        const int *xlist,
                                        const int *ylist,
                                        DATASET *dset,
                                        double *tmp,
                                        double (*dbuiltin)
                                        (int, int, const double *),
                                        int (*ibuiltin)
                                        (int, int, const double *),
                                        fncall *fc,
                                        int just_count,
                                        int *err)
{
    gretl_matrix *m = NULL;
    gretl_matrix **listvals;
    gretl_matrix *yvals;
    int n = sample_size(dset);
    int match;
    int skipnull = 0;
    int countcol = 1;
    int maxcases;
    int ny, nx, mcols;
    int *idx;
    double *valvec;
    double fx;
    int i, j, k, t, ni, ii;

    /* note:
       - to skip null cases in the output matrix, set skipnull = 1
       - to eliminate the case-count column, set countcol = 0
    */

    ny = ylist == NULL ? 1 : ylist[0];

    valvec = malloc(ny * sizeof *valvec);
    idx = malloc(ny * sizeof *idx);
    listvals = malloc(ny * sizeof *listvals);

    if (valvec == NULL || idx == NULL || listvals == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    for (j=0; j<ny; j++) {
        listvals[j] = NULL;
    }

    /* For @y (or each member of @ylist), create a vector holding
       its distinct values. As we go, count the combinations of
       y-values and initialize the @valvec and @idx arrays.
    */

    maxcases = 1;
    for (j=0; j<ny && !*err; j++) {
        if (ylist != NULL) {
            y = dset->Z[ylist[j+1]] + dset->t1;
        }
        listvals[j] = gretl_matrix_values(y, n, OPT_S, err);
        if (!*err) {
            maxcases *= listvals[j]->rows;
            valvec[j] = listvals[j]->val[0];
            idx[j] = 0;
        }
    }

    if (just_count) {
        x = NULL;
        xlist = NULL;
        skipnull = 0;
        countcol = 0;
        nx = 1;
    } else {
        nx = xlist == NULL ? 1 : xlist[0];
    }

    mcols = ny + nx + countcol;

    /* Allocate a matrix with enough rows to hold all the y-value
       combinations (maxcases) and enough columns to hold a record of
       the @y values, a count of matching cases, and the value(s) of
       f(x).
    */

    m = gretl_zero_matrix_new(maxcases, mcols);
    if (m == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    ii = 0;
    for (i=0; i<maxcases; i++) {
        ni = 0;
        for (k=0; k<nx && !*err; k++) {
            /* fill tmp with x for {y} == {yi} */
            if (xlist != NULL) {
                x = dset->Z[xlist[k+1]] + dset->t1;
            }
            ni = 0;
            for (t=0; t<n; t++) {
                match = 1;
                for (j=0; j<ny; j++) {
                    if (ylist != NULL) {
                        y = dset->Z[ylist[j+1]] + dset->t1;
                    }
                    if (y[t] != valvec[j]) {
                        match = 0;
                        break;
                    }
                }
                if (match) {
                    if (x != NULL) {
                        tmp[ni] = x[t];
                    }
                    ni++;
                }
            }
            if (just_count) {
                gretl_matrix_set(m, ii, ny, ni);
                break;
            } else if (ni == 0 && skipnull) {
                /* exclude cases where the obs count is 0 */
                break;
            } else {
                /* aggregate x at current y values */
                if (dbuiltin != NULL) {
                    fx = (*dbuiltin)(0, ni-1, tmp);
		} else if (ibuiltin != NULL) {
		    fx = (double) (*ibuiltin)(0, ni-1, tmp);
                } else {
                    double *pfx = &fx;

                    dset->t2 = ni-1;
                    *err = gretl_function_exec(fc, GRETL_TYPE_DOUBLE,
                                               dset, &pfx, NULL);
                    fx = *pfx;
                }
                gretl_matrix_set(m, ii, ny+k+countcol, fx);
            }
        }

        if (ni > 0 || !skipnull) {
            for (j=0; j<ny; j++) {
                gretl_matrix_set(m, ii, j, valvec[j]);
            }
            if (countcol) {
                gretl_matrix_set(m, ii, ny, ni);
            }
            ii++;
        }

        if (*err) {
            break;
        }

        /* set up the next y-values array */
        for (j=ny-1; j>=0; j--) {
            yvals = listvals[j];
            if (idx[j] == yvals->rows - 1) {
                /* restart the index at this level and pass
                   the buck */
                idx[j] = 0;
                valvec[j] = yvals->val[0];
            } else {
                /* increment the index at this level and
                   we're done */
                idx[j] += 1;
                valvec[j] = yvals->val[idx[j]];
                break;
            }
        }
    }

    if (skipnull && ii < maxcases) {
        /* the matrix contains some null cases: shrink it */
        m = delete_null_cases(m, ii, err);
    }

 bailout:

    free(valvec);
    free(idx);

    if (listvals != NULL) {
        for (j=0; j<ny; j++) {
            gretl_matrix_free(listvals[j]);
        }
        free(listvals);
    }

    return m;
}

/* Add suitable column names to the matrix to be returned
   by aggregate()
*/

static void aggr_add_colnames (gretl_matrix *m,
                               const int *ylist,
                               const int *xlist,
                               const DATASET *dset,
                               int just_count)
{
    char **S = NULL;
    int i, j, n = m->cols;
    int err = 0;

    S = strings_array_new(n);
    if (S == NULL) {
        return;
    }

    j = 0;

    if (ylist != NULL && ylist[0] > 0) {
        char **Sy = gretl_list_get_names_array(ylist, dset, &err);

        if (!err) {
            for (i=0; i<ylist[0]; i++) {
                S[j++] = Sy[i];
                Sy[i] = NULL;
            }
            free(Sy);
        }
    } else {
        S[j++] = gretl_strdup("byvar");
    }

    if (!err) {
        S[j++] = gretl_strdup("count");
    }

    if (!err && j < n) {
        if (xlist != NULL && xlist[0] > 0) {
            char **Sx = gretl_list_get_names_array(xlist, dset, &err);

            if (!err) {
                for (i=0; i<xlist[0]; i++) {
                    S[j++] = Sx[i];
                    Sx[i] = NULL;
                }
                free(Sx);
            }
        } else {
            S[j] = gretl_strdup("f(x)");
        }
    }

    if (!err) {
        gretl_matrix_set_colnames(m, S);
    } else {
        strings_array_free(S, n);
    }
}

/* Set up a call structure so we can call a user-defined
   function in the context of aggregate() directly, without
   having to go via "genr". Note that in this context the
   function is likely to be called repeatedly.
*/

static fncall *get_user_aggrby_call (const char *s,
                                     int nparam,
                                     double *tmp,
                                     int *err)
{
    ufunc *uf = get_user_function_by_name(s);
    fncall *fc = NULL;

    if (uf == NULL) {
	gretl_errmsg_sprintf("aggregate: '%s' is not a supported function", s);
        *err = E_INVARG;
    } else {
        GretlType rt = user_func_get_return_type(uf);

        if (rt != GRETL_TYPE_DOUBLE || fn_n_params(uf) != nparam) {
            *err = E_INVARG;
        }
    }

    if (!*err) {
        fc = user_func_get_fncall(uf);
        if (fc == NULL) {
            *err = E_ALLOC;
        } else if (tmp != NULL) {
            *err = push_function_arg(fc, NULL, NULL,
                                     GRETL_TYPE_SERIES,
                                     tmp);
        } else {
            ; /* not ready yet! */
        }
    }

    return fc;
}

static int get_aggregator (const char *fname,
                           double (**dbuiltin) (int, int, const double *),
                           int (**ibuiltin) (int, int, const double *))
{
    int f = function_lookup(fname);

    switch (f) {
    case F_SUM:
        *dbuiltin = gretl_sum;
        break;
    case F_SUMALL:
        *dbuiltin = series_sum_all;
        break;
    case F_MEAN:
        *dbuiltin = gretl_mean;
        break;
    case F_SD:
        *dbuiltin = gretl_stddev;
        break;
    case F_VCE:
        *dbuiltin = gretl_variance;
        break;
    case F_SST:
        *dbuiltin = gretl_sst;
        break;
    case F_SKEWNESS:
        *dbuiltin = gretl_skewness;
        break;
    case F_KURTOSIS:
        *dbuiltin = gretl_kurtosis;
        break;
    case F_MIN:
        *dbuiltin = gretl_min;
        break;
    case F_MAX:
        *dbuiltin = gretl_max;
        break;
    case F_MEDIAN:
        *dbuiltin = gretl_median;
        break;
    case F_GINI:
        *dbuiltin = gretl_gini;
        break;
    case F_NOBS:
        *dbuiltin = series_get_nobs;
        break;
    case F_ISCONST:
        *ibuiltin = gretl_isconst;
        break;
    case F_ISDUMMY:
        *ibuiltin = gretl_isdummy;
        break;
    default:
        break;
    }

    return *dbuiltin != NULL || *ibuiltin != NULL;
}

/**
 * aggregate_by:
 * @x: data array.
 * @y: discrete variable.
 * @xlist: list of x series or NULL.
 * @ylist: list of y series or NULL.
 * @fncall: the name of the aggregation function.
 * @dset: data set information.
 * @err: location to receive error code.
 *
 * Aggregates one or more data series (x) "by" the values of
 * one or more discrete series (y). In general either @x or
 * @xlist should be non-NULL, and one of @y or @ylist should
 * be non-NULL. (If @xlist is non-NULL then @x is ignored,
 * and similarly for @ylist and @y). For an account of the
 * matrix that is returned, see the help for gretl's
 * "aggregate" command.
 *
 * Returns: allocated matrix, or NULL on failure.
 */

gretl_matrix *aggregate_by (const double *x,
                            const double *y,
                            const int *xlist,
                            const int *ylist,
                            const char *fname,
                            DATASET *dset,
                            int *err)
{
    gretl_matrix *m = NULL;
    double *tmp = NULL;
    double (*dbuiltin) (int, int, const double *) = NULL;
    int (*ibuiltin) (int, int, const double *) = NULL;
    fncall *fc = NULL;
    int just_count = 0;
    int n;

    if (fname == NULL || !strcmp(fname, "null")) {
        just_count = 1;
    }

    if ((y == NULL && ylist == NULL) ||
        (!just_count && x == NULL && xlist == NULL)) {
        *err = E_DATA;
        return NULL;
    }

    if (!just_count) {
        get_aggregator(fname, &dbuiltin, &ibuiltin);
    }

    n = sample_size(dset);

    if (just_count) {
        ; /* nothing to do here */
    } else {
        tmp = malloc(n * sizeof *tmp);
        if (tmp == NULL) {
            *err = E_ALLOC;
        } else if (dbuiltin == NULL && ibuiltin == NULL) {
            fc = get_user_aggrby_call(fname, 1, tmp, err);
        }
    }

    if (!*err) {
        int save_t2 = dset->t2;

        x = (x == NULL)? NULL : x + dset->t1;
        y = (y == NULL)? NULL : y + dset->t1;
        m = real_aggregate_by(x, y, xlist, ylist, dset, tmp,
                              dbuiltin, ibuiltin, fc, just_count,
			      err);
        dset->t2 = save_t2;
    }

    if (m != NULL && *err) {
        gretl_matrix_free(m);
        m = NULL;
    }

    if (m != NULL) {
        aggr_add_colnames(m, ylist, xlist, dset, just_count);
    }

    free(tmp);

    return m;
}

/* Remedial function for the case where we try the matrix variant
   of aggregate but the aggregator function turns out not to be
   a built-in: we fall back to the "convert to series" case.
*/

static gretl_matrix *series_aggregate (const gretl_matrix *X,
                                       const gretl_matrix *y,
                                       const char *func,
                                       int *err)
{
    gretl_matrix *A = NULL;
    DATASET *dset = NULL;
    int *xlist = NULL;
    int *ylist = NULL;

    dset = matrix_dset_plus_lists(X, y, &xlist, &ylist, err);

    if (!*err) {
        A = aggregate_by(NULL, NULL, xlist, ylist, func, dset, err);
    }

    destroy_dataset(dset);
    free(xlist);
    free(ylist);

    return A;
}

#define AGGR_SORT_ALL 1

#if AGGR_SORT_ALL

/* Design: start by forming (y ~ X), where y is a vector, and sorting
   this matrix by its first column. This is an expensive operation,
   but once it's done the rest of the calculation is very
   straightforward.
*/

gretl_matrix *matrix_aggregate (const gretl_matrix *X,
                                const gretl_matrix *y,
                                const char *func,
                                int *err)
{
    double (*dbuiltin) (int, int, const double *) = NULL;
    int (*ibuiltin) (int, int, const double *) = NULL;
    gretl_matrix *Tmp;
    gretl_matrix *M;
    gretl_matrix *ret;
    double *z, *zsave;
    double rkj;
    int *counts = NULL;
    int nx = X->cols;
    int nby = 0;
    int i, j, k, nk;

    if (X == NULL || y == NULL || X->rows != y->rows) {
        *err = E_NONCONF;
    } else if (!get_aggregator(func, &dbuiltin, &ibuiltin)) {
        /* can't handle a user-defined function here */
        return series_aggregate(X, y, func, err);
    }

    if (*err) {
        return NULL;
    }

    /* allocate space for y ~ X */
    Tmp = gretl_matrix_alloc(X->rows, 1 + X->cols);
    if (Tmp == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    /* make M = (y ~ X), sorted by y */
    memcpy(Tmp->val, y->val, X->rows * sizeof(double));
    memcpy(Tmp->val + X->rows, X->val, X->rows * nx * sizeof(double));
    M = gretl_matrix_sort_by_column(Tmp, 0, err);
    gretl_matrix_free(Tmp);
    if (*err) {
        return NULL;
    }

    /* determine the number of 'by' values */
    nby = 0;
    for (i=0; i<M->rows; i++) {
        if (i == 0 || M->val[i] > M->val[i-1]) {
            nby++;
        }
    }

    /* allocate the return matrix */
    ret = gretl_zero_matrix_new(nby, 2 + nx);
    if (ret == NULL) {
        *err = E_ALLOC;
        gretl_matrix_free(M);
        return NULL;
    }

    /* get the per-value counts, filling the first two columns
       of @ret as we go
    */
    counts = calloc(nby, sizeof *counts);
    counts[0] = 1;
    for (i=1, k=0; i<M->rows; i++) {
        if (M->val[i] > M->val[i-1]) {
            gretl_matrix_set(ret, k, 0, M->val[i-1]);
            gretl_matrix_set(ret, k, 1, (double) counts[k]);
            counts[++k] = 1;
            if (i == M->rows - 1) {
                gretl_matrix_set(ret, k, 0, M->val[i]);
                gretl_matrix_set(ret, k, 1, 1.0);
            }
        } else {
            counts[k] += 1;
            if (i == M->rows - 1) {
                gretl_matrix_set(ret, k, 0, M->val[i]);
                gretl_matrix_set(ret, k, 1, (double) counts[k]);
            }
        }
    }

    z = M->val + M->rows; /* the @X portion of @M */

    /* do the actual aggregation */
    for (k=0; k<nby; k++) {
        zsave = z;
        nk = counts[k] - 1;
        for (j=0; j<nx; j++) {
             if (dbuiltin) {
                rkj = (*dbuiltin)(0, nk, z);
            } else {
                rkj = (double) (*ibuiltin)(0, nk, z);
            }
            gretl_matrix_set(ret, k, j+2, rkj);
            /* skip to the next column */
            z += M->rows;
        }
        /* move to the next block of rows */
        z = zsave + counts[k];
    }

    gretl_matrix_free(M);
    free(counts);

    return ret;
}

#else /* not AGGR_SORT_ALL */

gretl_matrix *matrix_aggregate (const gretl_matrix *X,
                                const gretl_matrix *y,
                                const char *func,
                                int *err)
{
    double (*dbuiltin) (int, int, const double *) = NULL;
    int (*ibuiltin) (int, int, const double *) = NULL;
    gretl_matrix *yvals = NULL;
    gretl_matrix *sel = NULL;
    gretl_matrix *xi = NULL;
    gretl_matrix *ret = NULL;
    double by, rkj;
    double *z;
    int n, ni, nx, nby;
    int i, j, k;

    if (X == NULL || y == NULL || X->rows != y->rows) {
        *err = E_NONCONF;
    } else if (!get_aggregator(func, &dbuiltin, &ibuiltin)) {
        /* can't handle a user-defined function here */
        return series_aggregate(X, y, func, err);
    }
    if (!*err) {
        yvals = gretl_matrix_values(y->val, y->rows, OPT_S, err);
    }

    if (*err) {
        return NULL;
    }

    n = y->rows;
    nby = yvals->rows;
    nx = X->cols;

    ret = gretl_matrix_alloc(nby, 2 + nx);
    sel = gretl_matrix_alloc(n, 1);

    if (ret == NULL || sel == NULL) {
        *err = E_ALLOC;
        gretl_matrix_free(yvals);
        gretl_matrix_free(ret);
        gretl_matrix_free(sel);
        return NULL;
    }

    for (i=0; i<nby && !*err; i++) {
        by = yvals->val[i];
        ni = 0;
        for (k=0; k<n; k++) {
            if (y->val[k] == by) {
                sel->val[k] = 1.0;
                ni++;
            } else {
                sel->val[k] = 0.0;
            }
        }
        xi = gretl_matrix_bool_sel(X, sel, 1, err);
        if (*err) {
            break;
        }
        gretl_matrix_set(ret, i, 0, by);
        gretl_matrix_set(ret, i, 1, (double) ni);
        z = xi->val;
        for (j=0; j<nx; j++) {
            if (dbuiltin) {
                rkj = (*dbuiltin)(0, ni-1, z);
            } else {
                rkj = (double) (*ibuiltin)(0, ni-1, z);
            }
            gretl_matrix_set(ret, i, j+2, rkj);
            z += ni;
        }
        gretl_matrix_free(xi);
    }

    gretl_matrix_free(yvals);
    gretl_matrix_free(sel);

    return ret;
}

#endif /* AGGR_SORT_ALL or not */

static int monthly_or_quarterly_dates (const DATASET *dset,
                                       int pd, double sd0, int T,
                                       double *x)
{
    char obstr[8];
    int mo, yr = floor(sd0);
    int t, mmax, dm;

    if (pd == 4) {
        sprintf(obstr, "%.1f", sd0 - yr);
        mo = 1 + (atoi(obstr + 2) - 1) * 3;
        mmax = 10;
        dm = 3;
    } else {
        sprintf(obstr, "%.2f", sd0 - yr);
        mo = atoi(obstr + 2);
        mmax = 12;
        dm = 1;
    }

    for (t=0; t<T; t++) {
        x[t] = 10000*yr + 100*mo + 1;
        if (mo == mmax) {
            yr++;
            mo = 1;
        } else {
            mo += dm;
        }
    }

    return 0;
}

static int annual_or_decennial_dates (const DATASET *dset,
                                      int pd, double sd0, int T,
                                      double *x)
{
    int t, yr = (int) sd0;

    for (t=0; t<T; t++) {
        x[t] = 10000*yr + 101;
        yr += pd;
    }

    return 0;
}

static int panel_daily_or_weekly (const DATASET *dset, double *x)
{
    DATASET tsset = {0};
    char datestr[16];
    int t, y, m, d, n;
    int err = 0;

    tsset.structure = TIME_SERIES;
    tsset.pd = dset->panel_pd;
    tsset.sd0 = dset->panel_sd0;
    tsset.n = dset->pd;
    calendar_date_string(tsset.stobs, 0, &tsset);

    for (t=0; t<dset->pd && !err; t++) {
        ntolabel(datestr, t, &tsset);
        n = sscanf(datestr, "%d-%d-%d", &y, &m, &d);
        if (n != 3) {
            err = E_DATA;
        } else {
            x[t] = 10000*y + 100*m + d;
        }
    }

    return err;
}

static int regular_daily_or_weekly (const DATASET *dset, double *x)
{
    char datestr[16];
    int t, y, m, d, n;
    int err = 0;

    for (t=0; t<dset->n && !err; t++) {
        ntolabel(datestr, t, dset);
#if 0
        fprintf(stderr, "regular: datestr = '%s'\n", datestr);
#endif
        n = sscanf(datestr, "%d-%d-%d", &y, &m, &d);
        if (n != 3) {
            err = E_DATA;
        } else {
            x[t] = 10000*y + 100*m + d;
        }
    }

    return err;
}

int fill_dataset_dates_series (const DATASET *dset, double *x)
{
    double sd0;
    int pd, T;
    int err = 0;

    if (dataset_is_panel(dset)) {
        /* looking at time dimension of panel */
        pd = dset->panel_pd;
        sd0 = dset->panel_sd0;
        T = dset->pd;
    } else {
        /* regular time series */
        pd = dset->pd;
        sd0 = dset->sd0;
        T = dset->n;
    }

    if (dataset_has_panel_time(dset)) {
        if (pd == 4 || pd == 12) {
            err = monthly_or_quarterly_dates(dset, pd, sd0, T, x);
        } else if (pd == 1 || pd == 10) {
            err = annual_or_decennial_dates(dset, pd, sd0, T, x);
        } else if (pd == 5 || pd == 6 || pd == 7 || pd == 52) {
            err = panel_daily_or_weekly(dset, x);
        }
    } else if (calendar_data(dset)) {
        err = regular_daily_or_weekly(dset, x);
    } else if (quarterly_or_monthly(dset)) {
        err = monthly_or_quarterly_dates(dset, pd, sd0, T, x);
    } else if (annual_data(dset) || decennial_data(dset)) {
        err = annual_or_decennial_dates(dset, pd, sd0, T, x);
    } else {
        err = E_PDWRONG;
    }

    if (!err && dataset_is_panel(dset)) {
        /* expand the date series for all groups */
        int i, N = dset->n / dset->pd;
        size_t bytes = dset->pd * sizeof *x;
        double *dest = x + dset->pd;

        for (i=1; i<N; i++) {
            memcpy(dest, x, bytes);
            dest += dset->pd;
        }
    }

    return err;
}

int fill_day_of_week_array (double *dow,
                            const double *y,
                            const double *m,
                            const double *d,
                            const DATASET *dset)
{
    char tstr[16];
    int n, yt, mt, dt;
    int julian;
    int t, err = 0;

    for (t=dset->t1; t<=dset->t2 && !err; t++) {
        julian = 0;
        if (m == NULL && d == NULL) {
            /* ISO 8601 basic input */
            sprintf(tstr, "%d", (int) y[t]);
            n = sscanf(tstr, "%4d%2d%2d", &yt, &mt, &dt);
            if (n != 3) {
                err = E_INVARG;
                break;
            }
        } else {
            /* full broken-down input */
            yt = (int) y[t];
            mt = (int) m[t];
            dt = (int) d[t];
        }
        if (yt < 0) {
            yt = -yt;
            julian = 1;
        }
        dow[t] = legacy_day_of_week(yt, mt, dt, julian, &err);
    }

    return err;
}

int fill_isoweek_array (double *wknum,
                        const double *y,
                        const double *m,
                        const double *d,
                        const DATASET *dset)
{
    char tstr[16];
    int n, yt, mt, dt;
    int t, err = 0;

    for (t=dset->t1; t<=dset->t2 && !err; t++) {
        if (m == NULL && d == NULL) {
            sprintf(tstr, "%d", (int) y[t]);
            n = sscanf(tstr, "%4d%2d%2d", &yt, &mt, &dt);
            if (n != 3) {
                err = E_INVARG;
                break;
            }
        } else {
            yt = (int) y[t];
            mt = (int) m[t];
            dt = (int) d[t];
        }
        wknum[t] = iso_week_number(yt, mt, dt, &err);
    }

    return err;
}

/**
 * empirical_cdf:
 * @y: array to process
 * @n: length of @y.
 * @err: location to receive error code on failure.
 *
 * Calculates the empirical CDF of @y, skipping any missing values.
 *
 * Returns: on successful completion, a matrix with row
 * dimension equal to the number of unique values in @y, and two
 * columns, the first containing the unique values of @y and the
 * second the cumulative relative frequency. On failure, returns
 * NULL.
 */

gretl_matrix *empirical_cdf (const double *y, int n, int *err)
{
    gretl_matrix *m = NULL;
    double *sy = NULL;
    double *crf = NULL;
    int n_le, n_u, n_ok = 0;
    int percentages = 0;
    int i, j, k, kbak;

    /* count non-missing values */
    for (i=0; i<n; i++) {
        if (!na(y[i])) {
            n_ok += 1;
        }
    }

    if (n_ok == 0) {
        *err = E_MISSDATA;
        return NULL;
    }

    /* allocate full array for sorting */
    sy = malloc(n_ok * sizeof *sy);
    if (sy == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    j = 0;
    for (i=0; i<n; i++) {
        if (!na(y[i])) {
            sy[j++] = y[i];
        }
    }

    qsort(sy, n_ok, sizeof *sy, gretl_compare_doubles);

    /* count unique values */
    n_u = 1;
    for (i=1; i<n_ok; i++) {
        if (sy[i] != sy[i-1]) {
            n_u++;
        }
    }

    m = gretl_matrix_alloc(n_u, 2);
    if (m == NULL) {
        *err = E_ALLOC;
        free(sy);
        return NULL;
    }

    crf = m->val + n_u;

    /* The first column of @m holds unique values,
       the second cumulative relative frequency
    */
    j = kbak = n_le = 0;
    for (i=0; i<n_ok; i++) {
        if (i == 0 || sy[i] != sy[i-1]) {
            gretl_matrix_set(m, j, 0, sy[i]);
            for (k=kbak; k<n_ok; k++) {
                if (sy[k] <= sy[i]) {
                    n_le++;
                } else {
                    break;
                }
            }
            if (percentages) {
                crf[j++] = 100 * n_le / (double) n_ok;
            } else {
                crf[j++] = n_le / (double) n_ok;
            }
            kbak = k;
        }
    }

    free(sy);

    return m;
}

/* Convert from ISO 8601 (daily) dates to year, quarter or month
   (discarding excess precision), if possible.
*/

static int handle_excess_precision (int *ymd1, int *ymd2, int pd,
                                    char *obs1, char *obs2)
{
    int err = 0;

    if (ymd1[2] != 1 || ymd2[2] != 1) {
        /* day-of-month must be 1 in all cases */
        return E_INVARG;
    }

    *obs1 = *obs2 = '\0';

    if (pd == 1) {
        /* annual: the month must be 1 */
        if (ymd1[1] != 1 || ymd1[2] != 1) {
            err = E_INVARG;
        } else {
            sprintf(obs1, "%d", ymd1[0]);
            sprintf(obs2, "%d", ymd2[0]);
        }
    } else if (pd == 4) {
        /* quarterly: the month must start a quarter */
        int m1 = ymd1[1], m2 = ymd2[1];

        if ((m1 != 1 && m1 != 4 && m1 != 7 && m1 != 10) ||
            (m2 != 1 && m2 != 4 && m2 != 7 && m2 != 10)) {
            err = E_INVARG;
        } else {
            sprintf(obs1, "%d:%d", ymd1[0], 1 + (m1-1)/3);
            sprintf(obs2, "%d:%d", ymd2[0], 1 + (m2-1)/3);
        }
    } else {
        /* pd == 12: monthly */
        int m1 = ymd1[1], m2 = ymd2[1];

        if (m1 < 1 || m1 > 12 || m2 < 1 || m2 > 12) {
            err = E_INVARG;
        } else {
            sprintf(obs1, "%d:%02d", ymd1[0], m1);
            sprintf(obs2, "%d:%02d", ymd2[0], m2);
        }
    }

    return err;
}

int sample_span (const char *stobs, const char *endobs,
                 int pd, int *err)
{
    char obs1[16], obs2[16];
    DATASET dset = {0};
    int ymd1[3], ymd2[3];
    int n, nf, t2, span = -1;

    if (pd == 1 || pd == 12 || (pd >= 4 && pd <= 7) || pd == 52) {
        ; /* OK, supported frequency */
    } else {
        *err = E_INVARG;
        return -1;
    }

    n = strlen(stobs);

    if (n > 10 || strlen(endobs) != n) {
        /* invalid and/or inconsistent observation strings */
        *err = E_INVARG;
        return -1;
    }

    strcpy(obs1, stobs);
    strcpy(obs2, endobs);

    if (n == 10) {
        /* should be ISO 8601 dates */
        nf = sscanf(obs1, "%d-%d-%d", &ymd1[0], &ymd1[1], &ymd1[2]);
        if (nf != 3) {
            *err = E_INVARG;
        } else {
            nf = sscanf(obs2, "%d-%d-%d", &ymd2[0], &ymd2[1], &ymd2[2]);
            if (nf != 3) {
                *err = E_INVARG;
            }
        }
        if (!*err && (pd == 1 || pd == 4 || pd == 12)) {
            *err = handle_excess_precision(ymd1, ymd2, pd, obs1, obs2);
        }
    }

    if (*err) {
        return -1;
    }

    strcpy(dset.stobs, obs1);
    dset.pd = pd;
    dset.structure = TIME_SERIES;

    if (n == 10 && ((pd >= 5 && pd <= 7) || pd == 52)) {
        /* validate daily data */
        int ed1, ed2;

        ed1 = epoch_day_from_ymd(ymd1[0], ymd1[1], ymd1[2]);
        ed2 = epoch_day_from_ymd(ymd2[0], ymd2[1], ymd2[2]);

        if (ed1 > 0 && ed2 > 0 && (pd == 5 || pd == 6)) {
            /* validate days-of-week */
            int wd1 = weekday_from_epoch_day(ed1);
            int wd2 = weekday_from_epoch_day(ed2);

            if (wd1 > pd) {
                ed1 = 0;
            } else if (wd2 > pd) {
                ed2 = 0;
            }
        }
        if (ed1 == 0) {
            gretl_errmsg_sprintf(_("Invalid observation %s"), stobs);
            *err = E_INVARG;
        } else if (ed2 == 0) {
            gretl_errmsg_sprintf(_("Invalid observation %s"), endobs);
            *err = E_INVARG;
        } else {
            dset.sd0 = (double) ed1;
        }
    }

    if (*err) {
        span = -1;
    } else {
        t2 = dateton(obs2, &dset);
        if (t2 >= 0) {
            span = t2 + 1; /* t2 is 0-based */
        } else {
            *err = E_INVARG;
            span = -1;
        }
    }

    return span;
}

/* simple driver function */

gretl_matrix *gretl_matrix_vector_stat (const gretl_matrix *m,
                                        GretlVecStat vs, int rowwise,
                                        int skip_na, int *err)
{
    if (m->is_complex) {
        return gretl_cmatrix_vector_stat(m, vs, rowwise, skip_na, err);
    } else {
        return gretl_rmatrix_vector_stat(m, vs, rowwise, skip_na, err);
    }
}

/* Fill @v with unique random draws from {1,2,...,@n} */

int fill_permutation_vector (gretl_vector *v, int n)
{
    int *pool = NULL;
    int n_tail;
    unsigned u;
    int i, k;

    if (n <= 0) {
        return E_INVARG;
    }

    k = gretl_vector_get_length(v);
    if (k <= 0 || k > n) {
        gretl_errmsg_sprintf(_("Invalid number of draws %d"), k);
        return E_INVARG;
    }

    pool = malloc(n * sizeof *pool);
    if (pool == NULL) {
        return E_ALLOC;
    }

    /* initialize selection pool */
    for (i=0; i<n; i++) {
        pool[i] = i + 1;
    }

    for (i=0; i<k; i++) {
        u = gretl_rand_int_max(n);
        v->val[i] = pool[u];
        /* remove pool[u] from the selectable set */
        if (u < n - 1) {
            n_tail = n - u - 1;
            memmove(pool + u, pool + u + 1, n_tail * sizeof *pool);
        }
        n--;
    }

    free(pool);

    return 0;
}

enum {
    SUBST_SIMPLE,
    SUBST_BTREE
};

static double subst_val_via_tree (double x, const double *x0, int n0,
                                  const double *x1, int n1)
{
    static BTree *tree;

    if (x0 == NULL) {
        /* cleanup */
        gretl_btree_destroy(tree);
        tree = NULL;
        return 0;
    }

    if (tree == NULL) {
        /* allocate and populate tree */
        double x1val;
        int i;

        tree = gretl_btree_new();
        for (i=0; i<n0; i++) {
            x1val = n1 == 1 ? *x1 : x1[i];
            gretl_btree_insert(tree, x0[i], x1val);
        }
    }

    /* do the actual lookup */
    x = gretl_btree_lookup(tree, x);

    return x;
}

static int select_subst_method (int n, int nfind)
{
#if 1
    /* for testing */
    if (getenv("REPLACE_USE_BTREE") != NULL) {
        return SUBST_BTREE;
    } else if (getenv("REPLACE_NAIVE") != NULL) {
        return SUBST_SIMPLE;
    }
#endif
    /* The idea here is that it's worth using a binary
       tree mapping "find" to "replace" values only if
       the problem is big enough. Testing seems to
       indicate that the following condition for "big
       enough" may be roughly appropriate.
    */
    if (nfind > 11 && n >= 80) {
        return SUBST_BTREE;
    } else {
        return SUBST_SIMPLE;
    }
}

/* Given an original value @x, see if it matches any of the @n0 values
   in @x0.  If so, return the substitute value from @x1, otherwise
   return the original.
*/

static double subst_val (double x, const double *x0, int n0,
                         const double *x1, int n1,
                         int method)
{
    if (method & SUBST_BTREE) {
        x = subst_val_via_tree(x, x0, n0, x1, n1);
    } else {
        int i;

        for (i=0; i<n0; i++) {
            if (x == x0[i]) {
                return (n1 == 1)? *x1 : x1[i];
            } else if (isnan(x) && isnan(x0[i])) {
                /* we'll count this as a match */
                return (n1 == 1)? *x1 : x1[i];
            }
        }
    }

    return x;
}

/**
 * substitute_values:
 * @dest: target array.
 * @src: source array.
 * @n: length of @src and @dest.
 * @v0: array of "find" values.
 * @n0: length of @v0.
 * @v1: array of "replace" values.
 * @n1: length of @v1.
 *
 * For each of the @n elements in @src, determine if it appears
 * in @v0: if so, write to @targ the corresponding element of
 * @v1, otherwise write the @src element to @targ. The method
 * employed is either "naive" lookup, or if the problem is of
 * sufficient size a binary tree.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int substitute_values (double *dest, const double *src, int n,
                       const double *v0, int n0,
                       const double *v1, int n1)
{
    int i, method = select_subst_method(n, n1);

    for (i=0; i<n; i++) {
        dest[i] = subst_val(src[i], v0, n0, v1, n1, method);
    }
    if (method & SUBST_BTREE) {
        /* cleanup call */
        subst_val(0, NULL, 0, NULL, 0, method);
    }

    return 0;
}

/* Write into @vname the putative name for a series to be made from
   column @i of a matrix. Check for validity of the name and return
   non-zero in case of error.

   The naming procedure is as follows.

   * If @prefix is non-NULL, name series @i as <prefix><i>.
   * Otherwise if @S (matrix column names) is non-NULL, name
     the series using S[i].
   * Otherwise name the series as column<i>.
*/

static int get_column_name (char *vname, int i,
                            const char *prefix,
                            int numlen,
                            const char **S)
{
    int done = 0;
    int err = 0;

    if (prefix != NULL) {
        if (numlen + strlen(prefix) >= VNAMELEN) {
            err = E_INVARG;
        } else {
            err = check_varname(prefix);
        }
        if (!err) {
            sprintf(vname, "%s%d", prefix, i+1);
            done = 1;
        }
    } else if (S != NULL && check_varname(S[i]) == 0) {
        strcpy(vname, S[i]);
        done = 1;
    }

    if (!err && !done) {
        sprintf(vname, "column%d", i+1);
    }

    return err;
}

static int name_collision (const char *vname)
{
    GretlType t = user_var_get_type_by_name(vname);
    gchar *msg = name_conflict_message(vname, t);

    gretl_errmsg_set(msg);
    g_free(msg);

    return E_TYPES;
}

/**
 * list_from_matrix:
 * @m: source matrix.
 * @prefix: prefix for series names or NULL.
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * Constructs a series from each column of @m and returns
 * a list of the series.
 *
 * Returns: list on successful completion, NULL otherwise.
 */

int *list_from_matrix (const gretl_matrix *m,
                       const char *prefix,
                       DATASET *dset, int *err)
{
    int *list = NULL;
    int n, k;

    if (m != NULL && m->is_complex) {
        *err = E_CMPLX;
        return NULL;
    } else if (gretl_is_null_matrix(m)) {
        return gretl_null_list();
    }

    n = m->rows;
    k = m->cols;

    if (n != dset->n && n != sample_size(dset)) {
        *err = E_NONCONF;
    } else {
        const char **S = gretl_matrix_get_colnames(m);
        char vname[VNAMELEN];
        int numlen = ceil(log10(k));
        int orig_v = dset->v;
        int i, vi, nadd = 0;

        list = gretl_list_new(k);

        /* check the names to be assigned to series */
        for (i=0; i<k && !*err; i++) {
            *err = get_column_name(vname, i, prefix, numlen, S);
            vi = current_series_index(dset, vname);
            if (vi > 0) {
                list[i+1] = vi;
            } else if (gretl_is_user_var(vname)) {
                *err = name_collision(vname);
            } else {
                nadd++;
            }
        }

        if (!*err && nadd > 0) {
            *err = dataset_add_NA_series(dset, nadd);
        }

        if (!*err) {
            double *dest, *src = m->val;
            size_t csize = n * sizeof *src;
            int vnew = orig_v;

            if (n == dset->n) {
                src += dset->t1;
            }

            for (i=0; i<k; i++) {
                vi = list[i+1];
                if (vi == 0) {
                    /* a new series */
                    get_column_name(vname, i, prefix, numlen, S);
                    vi = vnew++;
                    list[i+1] = vi;
                    strcpy(dset->varname[vi], vname);
                }
                dest = dset->Z[vi] + dset->t1;
                memcpy(dest, src, csize);
                src += n;
            }
        }

        if (*err && list != NULL) {
            free(list);
            list = NULL;
        }
    }

    return list;
}

#define SPHCORR_DBG 0

static int sphc_unitvec (double *f,
                         const double *omega, int h,
                         gretl_matrix *J,
                         double *work)
{
    double *s, *k, *c;
    double x, ki = 1.0;
    int i, j;

    s = work;
    k = s + h;
    c = k + h + 1;

    for (i=0; i<h; i++) {
        s[i] = sin(omega[i]);
        k[i] = ki;
        ki *= s[i];
        c[i] = cos(omega[i]);
    }

    k[h] = ki;
    c[h] = 1.0;

    for (i=0; i<=h; i++) {
        f[i] = c[i] * k[i];
#if SPHCORR_DBG
        fprintf(stderr, "f[%d] = %g\n", i, f[i]);
#endif
    }

    if (J != NULL) {
        gretl_matrix_zero(J);
        /* we reuse @s for the tangents */
        for (i=0; i<h; i++) {
            s[i] = tan(omega[i]);
        }
        for (i=1; i<=h; i++) {
            for (j=0; j<i; j++) {
                gretl_matrix_set(J, i, j, k[i] / s[j]);
            }
        }
        for (i=0; i<h; i++) {
            for (j=0; j<i; j++) {
                x = gretl_matrix_get(J, i, j);
                gretl_matrix_set(J, i, j, c[i] * x);
#if 0
                fprintf(stderr, "J[%d,%d] = c[%d]*x\n", i, j, i);
#endif
            }
            gretl_matrix_set(J, i, i, -k[i+1]);
        }
    }

    return 0;
}

/**
 * omega_from_R:
 * @R: correlation matrix.
 * @err: location to receive error code.
 *
 * Expresses the (n x n) correlation matrix R in spherical coordinates,
 * via m angles, where m = n*(n-1)/2.
 *
 * Returns: an m-element vector if successful, NULL on failure.
 */

gretl_matrix *omega_from_R (const gretl_matrix *R, int *err)
{
    gretl_matrix *K, *omega = NULL;
    int n = R->rows;
    int m = n * (n-1) / 2;

    if (R->cols != n) {
        *err = E_NONCONF;
        return omega;
    }

    K = gretl_matrix_copy(R);
    if (K == NULL) {
        *err = E_ALLOC;
        return omega;
    }

    *err = gretl_matrix_psd_root(K, 0);
    if (!*err) {
        *err = gretl_matrix_transpose_in_place(K);
    }

    if (!*err) {
        omega = gretl_matrix_alloc(m, 1);
        if (omega == NULL) {
            *err = E_ALLOC;
        }
    }

    if (!*err) {
        int i, j, k, l = 0;
        double theta, x, y, prod;
        int rank_def;

        for (i=1; i<n; i++) {
            k = i * n;
            prod = 1.0;
            rank_def = 0;
            for (j=0; j<i; j++) {
                if (prod < 1.0e-20) {
                    rank_def = 1;
                } else {
                    x = K->val[k++] / prod;
                    if (1.0 - fabs(x) < 1.0e-20) {
                        rank_def = x > 0 ? 1 : 2;
                    }
                }
                if (rank_def) {
                    omega->val[l++] = rank_def == 2 ? M_PI : 0;
                } else {
                    theta = acos(x);
                    omega->val[l++] = theta;
                    y = sin(theta);
                    prod *= y;
                }
            }
        }
    }

    gretl_matrix_free(K);

    return omega;
}

/* build a vector holding the indices for the vech Jacobian */

static int *funky_index (int n, int m)
{
    int i, j, l = 0;
    int sum;
    int *index;

    index = malloc((n + m) * sizeof *index);

    for (i=0; i<n; i++){
        index[l] = i;
        l++;
        sum = 0;
        for (j=1; j<=i; j++) {
            sum += n-j;
            index[l] = i + sum;
#if SPHCORR_DBG
            printf("funky_index: %d\n", index[l]);
#endif
            l++;
        }
    }

    return index;
}

static gretl_matrix *JR_from_Jk (gretl_matrix *Jk,
                                 gretl_matrix *K,
                                 int *err)
{
    /*
      This function computes the Jacobian of vech(R) wrt the spherical
      coordinates omega. Surprisingly, omega does not feature as an
      argument; the reason is that we use K instead (the Cholesky
      factor of R).

      Since we have the Jacobian of vech(K) wrt omega Jk, we use the
      relationship R = KK' to define the full Jacobian as JJ * Jk,
      where JJ is an appropriate subset of the rows and columns of

      (I + K_nn) (I \otimes K)
    */
    gretl_matrix *IkronK = NULL, *tmp = NULL;
    gretl_matrix *JJ = NULL, *ret = NULL;
    int i, j, n = K->rows;
    int h = n * (n+1)/2, m = n * (n-1)/2, k, l;
    int *ndxc = NULL, *ndxr = NULL;
    int ndxi;
    double x;

    /* allocate stuff */
    ndxc = malloc(2 * h * sizeof *ndxc);
    if (ndxc == NULL) {
        *err = E_ALLOC;
        return NULL;
    }
    ndxr = ndxc + h;

    l = 0;
    k = 0;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            if (i <= j) {
                ndxr[l] = k;
                ndxc[l++] = i + j*n;
            }
            k++;
        }
    }

    *err = 0;
    IkronK = gretl_matrix_I_kronecker_new(n, K, err);
    if (*err) {
        goto bailout;
    }

    JJ = gretl_matrix_commute(IkronK, n, n, 1, 1, err);
    if (*err) {
        goto bailout;
    }

    tmp = gretl_matrix_alloc(h, h);

    for (i=0; i<h; i++) {
        ndxi = ndxr[i];
        for (j=0; j<h; j++) {
            x = gretl_matrix_get(JJ, ndxi, ndxc[j]);
            gretl_matrix_set(tmp, i, j, x);
        }
    }

    ret = gretl_matrix_alloc(h, m);
    gretl_matrix_multiply(tmp, Jk, ret);

 bailout:

    free(ndxc);
    gretl_matrix_free(IkronK);
    gretl_matrix_free(JJ);
    gretl_matrix_free(tmp);

    return ret;
}

/**
 * R_from_omega:
 * @omega: angle vector.
 * @cholesky: boolean, compute cholesky factor instead of correlation matrix.
 * @J: location to receive Jacobian matrix if wanted, or NULL.
 * @err: location to receive error code.
 *
 * Returns an (n x n) correlation matrix R given its
 * spherical-coordinates representation omega, which should contain
 * m numbers between 0 and M_PI, where m = n*(n-1)/2.
 *
 * Returns: an n x n symmetric matrix if successful, NULL on failure.
 */

gretl_matrix *R_from_omega (const gretl_matrix *omega, int cholesky,
                            gretl_matrix *J, int *err)
{
    gretl_matrix *K, *R = NULL;
    int m = omega->rows;
    double der, tmp = 0.5 * (1.0 + sqrt(1 + 8*m));
    int n = nearbyint(tmp);
    int *Jindex = NULL;
    gretl_matrix *localJ = NULL;
    gretl_matrix *tmpJac = NULL;
    double *om, *f = NULL;
    double *work = NULL;
    int ii, jj, r_i, c_i;
    int i, j, h;

    if (fabs(tmp - n) > 1.0e-3) {
        *err = E_NONCONF;
        return NULL;
    }

    K = gretl_zero_matrix_new(n, n);
    gretl_matrix_set(K, 0, 0, 1.0);

    if (J != NULL) {
        /* pre-allocate Jacobian matrix */
        localJ = gretl_zero_matrix_new(m + n, m);
        Jindex = funky_index(n, m);
        tmpJac = gretl_matrix_alloc(n+1, n);
    }

#if SPHCORR_DBG
    if (J != NULL) {
        for (i=0; i<m+n; i++) {
            fprintf(stderr, "Jindex[%d] = %d\n", i, Jindex[i]);
        }
    }
#endif

    om = omega->val;
    work = malloc((n + 2*(n+1)) * sizeof *work);
    f = malloc(m * sizeof *f);

    r_i = 1;
    c_i = 0;

    for (i=0; i<n-1; i++) {
        h = i + 1;
        *err = sphc_unitvec(f, om, h, tmpJac, work);
        om += h;
        for (j=0; j<=h; j++) {
            gretl_matrix_set(K, h, j, f[j]);
        }
        if (J != NULL) {
            for (ii=0; ii<=h; ii++) {
                for (jj=0; jj<=i; jj++) {
                    der = gretl_matrix_get(tmpJac, ii, jj);
#if SPHCORR_DBG
                    printf("d[%d,%d] = %f\n", ii, jj, der);
                    printf("The place is [%d, %d]\n", Jindex[r_i], c_i + jj);
#endif
                    gretl_matrix_set(localJ, Jindex[r_i], c_i + jj, der);
                }
                r_i++;
            }
            c_i += h;
        }
    }

#if SPHCORR_DBG
    gretl_matrix_print(K, "K via f");
#endif

    free(work);
    free(f);

    if (J != NULL) {
        gretl_matrix_replace_content(J, localJ);
        gretl_matrix_free(localJ);
        gretl_matrix_free(tmpJac);
        free(Jindex);
    }

#if SPHCORR_DBG
    gretl_matrix_print(K, "reconstructed cholesky");
    if (do_jacobian) {
        gretl_matrix_print(*J, "Jacobian");
    }
#endif

    if (!cholesky) {
        gretl_matrix *foo;
        int err2;

        R = gretl_matrix_alloc(n, n);
        if (R == NULL) {
            *err = E_ALLOC;
        } else {
            gretl_matrix_multiply_mod(K, GRETL_MOD_NONE,
                                      K, GRETL_MOD_TRANSPOSE,
                                      R, GRETL_MOD_NONE);
        }
        if (J != NULL) {
            foo = JR_from_Jk(J, K, &err2);
            gretl_matrix_replace_content(J, foo);
            gretl_matrix_free(foo);
        }
        gretl_matrix_free(K);
    }

    return cholesky ? K : R;
}

/**
 * felogit_rec_loglik:
 *
 * @t: int, number of ones
 * @T: int, number of observations
 * @U: matrix with the index function
 * @X: matrix with the covariates
 *
 * Computes recursively the denominator of the fixed-effect logit
 * likelihood for one unit with @T observations and @t ones, as well
 * as its first derivative wrt the parameter vector beta. The matrix
 * @U is a vector containing the values of the index function for each
 * time period and @X is the corresponding matrix of covariates, so
 * that U[i] = exp(X[i,]*beta).
 *
 * The first element of the output vector is the likelihood
 * denominator (B in Gail et al's notation) and the subsequent k
 * elements contain the derivative of B wrt beta.
 *
 * The algorithm is described in Gail et al. (1981), "Likelihood
 * Calculations for Matched Case-Control Studies and Survival Studies
 * with Tied Death Times", Biometrika, Vol. 68 (3), pp. 703-707
 *
 * TODO: the function could be enhanced by also including the
 * recursive formula for the Hessian (unnumbered equation at the end
 * of page 706 in Gail et al).
 *
 * Returns: a (k+1)-vector if successful, NULL on failure.
 *
 */

gretl_matrix *felogit_rec_loglik (int t, int T,
                                  const gretl_matrix *U,
                                  const gretl_matrix *xi)
{
    gretl_matrix *ret;
    int k = xi->cols;

    ret = gretl_zero_matrix_new(1, k+1);

    if (t > T) {
        ; /* nothing to be done */
    } else if (t == 0) {
        ret->val[0] = 1;
    } else if (t == T) {
        double B, rj;
        int i, j;

        /* scalar B = prodc(U[1:T]) */
        B = 1.0;
        for (i=0; i<T; i++) {
            B *= U->val[i];
        }
        /* ret = (1 ~ sumc(xi[1:T,])) .* B */
        ret->val[0] = B;
        for (j=0; j<k; j++) {
            rj = 0;
            for (i=0; i<T; i++) {
                rj += gretl_matrix_get(xi, i, j);
            }
            ret->val[j+1] = B * rj;
        }
    } else {
        gretl_matrix *tmp1 = felogit_rec_loglik(t, T-1, U, xi);
        gretl_matrix *tmp2 = felogit_rec_loglik(t-1, T-1, U, xi);
        double xij, UT = U->val[T-1];
        double tmp20 = tmp2->val[0];
        int j;

        /* B = tmp[1,1] + UT * tmp[2,1] */
        ret->val[0] = tmp1->val[0] + UT * tmp20;

        for (j=1; j<=k; j++) {
            xij = gretl_matrix_get(xi, T-1, j-1);
            ret->val[j] = tmp1->val[j] + UT * (tmp2->val[j] + xij * tmp20);
        }

        gretl_matrix_free(tmp1);
        gretl_matrix_free(tmp2);
    }

    return ret;
}

static int val_map_search (double needle, const double *haystack,
                           int n, int offset)
{
    int m = n/2;

    if (needle == haystack[m]) {
        return m + offset;
    } else if (needle < haystack[0] || needle > haystack[n-1]) {
        return -1;
    } else if (needle < haystack[m]) {
        return val_map_search(needle, haystack, m, offset);
    } else {
        return val_map_search(needle, haystack + m, n - m, offset + m);
    }
}

/* Helper function for strftime() as applied to a series argument.  If
   the numeric input series @x contains any duplicated values we need
   a mapping from position in @x to position in an array of unique
   values for the purpose of stringifying the output series.
*/

int *maybe_get_values_map (const double *x, int n, int *pnv, int *err)
{
    gretl_vector *v = NULL;
    int *ret = NULL;
    int mono = 0;
    int n_ok = n;
    int i, j = -1;

    *pnv = 0;

    /* Is @x increasing or decreasing monotonically?
       This is a fairly cheap check for a sufficient
       condition for no repetition. We count missing
       values and check for bogus ones as we go.
    */
    for (i=0; i<n; i++) {
        if (na(x[i])) {
            n_ok--;
        } else if (x[i] < -6.2e10 || x[i] > 2.5e11) {
            /* outside of plausible time_t range */
            *err = E_INVARG;
            break;
        } else if (i > 0 && !na(x[i-1])) {
            if (j < 0) {
                mono = (x[i] > x[i-1])? 1 : (x[i] < x[i-1])? 2 : 0;
                j = i;
            } else if ((mono == 1 && x[i] <= x[i-1]) ||
                       (mono == 2 && x[i] >= x[i-1])) {
                mono = 0;
            }
        }
    }

    if (*err) {
        return NULL;
    } else if (n_ok == 0) {
        *err = E_MISSDATA;
        return NULL;
    } else if (mono) {
        /* no repetition: no map needed */
        *pnv = n_ok;
        return NULL;
    }

    v = gretl_matrix_values(x, n, OPT_S, err);

    if (!*err) {
        /* allocate the integer map */
        ret = malloc(n_ok * sizeof *ret);
        if (ret == NULL) {
            *err = E_ALLOC;
        }
    }

    if (!*err) {
        int nv = gretl_vector_get_length(v);

        for (i=0, j=0; i<n; i++) {
            if (na(x[i])) {
                continue;
            }
            /* find the position of x[i] in @v */
            ret[j++] = val_map_search(x[i], v->val, nv, 0);
        }
        *pnv = nv;
    }

    gretl_matrix_free(v);

    return ret;
}

/* Here we allow that @m1 or @m2 can be NULL. Optionally,
   the arguments @plist1 and/or @plist2 can be filled out
   with lists corresponding to the positions of the
   associated matrix columns in the returned dataset.
*/

DATASET *matrix_dset_plus_lists (const gretl_matrix *m1,
                                 const gretl_matrix *m2,
                                 int **plist1, int **plist2,
                                 int *err)
{
    DATASET *mdset = NULL;

    if (m1 != NULL && m2 != NULL) {
        mdset = gretl_dataset_from_matrices(m1, m2, err);
    } else if (m1 == NULL && m2 != NULL) {
        mdset = gretl_dataset_from_matrix(m2, NULL, OPT_B, err);
    } else if (m1 != NULL && m2 == NULL) {
        mdset = gretl_dataset_from_matrix(m1, NULL, OPT_B, err);
    } else {
        *err = E_DATA;
        return NULL;
    }

    if (mdset != NULL) {
        int i, j = 1;

        if (m1 != NULL && plist1 != NULL) {
            int *l1 = gretl_list_new(m1->cols);

            for (i=0; i<m1->cols; i++) {
                l1[i+1] = j++;
            }
            *plist1 = l1;
        }
        if (m2 != NULL && plist2 != NULL) {
            int *l2 = gretl_list_new(m2->cols);

            for (i=0; i<m2->cols; i++) {
                l2[i+1] = j++;
            }
            *plist2 = l2;
        }
    }

    return mdset;
}
