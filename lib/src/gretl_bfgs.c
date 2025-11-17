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

/* The BFGS and Newton optimizers plus the fdjac Jacobian function */

#include "libgretl.h"
#include "gretl_bfgs.h"
#include "libset.h"
#include "matrix_extra.h"
#include "usermat.h"
#include "uservar.h"
#include "gretl_func.h"

#include "../../minpack/minpack.h"
#include <float.h>
#include <errno.h>

#define BFGS_DEBUG 0

#define BFGS_MAXITER_DEFAULT 600

void BFGS_defaults (int *maxit, double *tol, int ci)
{
    *maxit = libset_get_int(BFGS_MAXITER);
    *tol = libset_get_user_tolerance(BFGS_TOLER);

    if (ci != MLE && ci != GMM && *maxit <= 0) {
        *maxit = 1000;
    }

    if (ci == PROBIT || ci == INTREG || ci == ARMA ||
        ci == NEGBIN || ci == DURATION) {
        if (na(*tol)) {
            *tol = 1.0e-12;
        }
    } else if (ci == TOBIT) {
        if (na(*tol)) {
            *tol = 1.0e-10; /* calibrated against Wm Greene */
        }
    } else if (ci == HECKIT) {
        if (na(*tol)) {
            *tol = 1.0e-09;
        }
    } else if (ci == GARCH) {
        if (na(*tol)) {
            *tol = 1.0e-13;
        }
    } else if (ci == MLE || ci == GMM) {
        if (*maxit <= 0) {
            *maxit = BFGS_MAXITER_DEFAULT;
        }
        if (na(*tol)) {
            *tol = libset_get_double(BFGS_TOLER);
        }
    }
}

static void free_triangular_array (double **m, int n)
{
    if (m != NULL) {
        int i;

        for (i=0; i<n; i++) {
            free(m[i]);
        }
        free(m);
    }
}

static double **triangular_array_new (int n)
{
    double **m = malloc(n * sizeof *m);
    int i;

    if (m != NULL) {
        for (i=0; i<n; i++) {
            m[i] = NULL;
        }
        for (i=0; i<n; i++) {
            m[i] = malloc((i + 1) * sizeof **m);
            if (m[i] == NULL) {
                free_triangular_array(m, n);
                return NULL;
            }
        }
    }

    return m;
}

/**
 * hessian_from_score:
 * @b: array of k parameter estimates.
 * @H: k x k matrix to receive the (negative) Hessian.
 * @gradfunc: function to compute gradient.
 * @cfunc: function to compute criterion (or NULL, see below).
 * @data: data to be passed to the @gradfunc callback.
 *
 * Uses the score function (@gradfunc) is to construct a
 * numerical approximation to the Hessian. This is primarily
 * intended for building a covariance matrix at convergence;
 * note that it may not work well at an arbitrary point in
 * the parameter space.
 *
 * Note that the only use of @cfunc within this function is
 * as an argument to be passed to @gradfunc. It is therefore
 * OK to pass NULL for @cfunc provided that @gradfunc does not
 * use its 4th argument, which corresponds to the BFGS_CRIT_FUNC
 * parameter.
 *
 * Returns: 0 on successful completion, non-zero error code
 * on error.
 */

int hessian_from_score (double *b, gretl_matrix *H,
                        BFGS_GRAD_FUNC gradfunc,
                        BFGS_CRIT_FUNC cfunc,
                        void *data)
{
    double *g, *splus, *sminus, *splus2, *sminus2;
    double x, den, b0, eps = 1.0e-05;
    int n = gretl_matrix_rows(H);
    int extra_precision = 0;
    int i, j, err = 0;

    char *s = getenv("H_EXTRA");

    if (s != NULL && *s != '\0') {
        extra_precision = 1;
        fprintf(stderr, "hessian_from_score: using extra precision\n");
    }

    if (extra_precision) {
        splus = malloc(5 * n * sizeof *splus);
        splus2 = splus + n;
        sminus = splus2 + n;
        sminus2 = sminus + n;
        g = sminus2 + n;
        den = 12*eps;
    } else {
        splus = malloc(3 * n * sizeof *splus);
        sminus = splus + n;
        g = sminus + n;
        splus2 = sminus2 = NULL;
        den = 2*eps;
    }

    if (splus == NULL) {
        return E_ALLOC;
    }

    for (i=0; i<n; i++) {
        b0 = b[i];
        b[i] = b0 + eps;
        err = gradfunc(b, g, n, cfunc, data);
        if (err) goto restore;
        for (j=0; j<n; j++) {
            splus[j] = g[j];
        }
        b[i] = b0 - eps;
        err = gradfunc(b, g, n, cfunc, data);
        if (err) goto restore;
        for (j=0; j<n; j++) {
            sminus[j] = g[j];
        }
        if (extra_precision) {
            b[i] = b0 - 2*eps;
            err = gradfunc(b, g, n, cfunc, data);
            if (err) goto restore;
            for (j=0; j<n; j++) {
                sminus2[j] = g[j];
            }
            b[i] = b0 + 2*eps;
            err = gradfunc(b, g, n, cfunc, data);
            if (err) goto restore;
            for (j=0; j<n; j++) {
                splus2[j] = g[j];
            }
        }
    restore:
        b[i] = b0;
        if (err) {
            break;
        }
        for (j=0; j<n; j++) {
            if (extra_precision) {
                x = -(splus2[j] - sminus2[j]) + 8*(splus[j] - sminus[j]);
            } else {
                x = splus[j] - sminus[j];
            }
            gretl_matrix_set(H, i, j, -x / den);
        }
    }

    if (!err) {
        gretl_matrix_xtr_symmetric(H);
    }

    free(splus);

    return err;
}

/**
 * hessian_inverse_from_score:
 * @b: array of parameter estimates.
 * @n: the number of elements in @b.
 * @gradfunc: function to compute gradient.
 * @cfunc: function to compute criterion.
 * @data: data to be passed to the @gradfunc callback.
 * @err: location to receive error code.
 *
 * A wrapper for hessian_from_score() which takes care of
 * (a) allocation of the Hessian and (b) inversion.
 *
 * Returns: the inverse of the (negative) Hessian on successful
 * completion, NULL on error.
 */

gretl_matrix *hessian_inverse_from_score (double *b, int n,
                                          BFGS_GRAD_FUNC gradfunc,
                                          BFGS_CRIT_FUNC cfunc,
                                          void *data, int *err)
{
    gretl_matrix *H = gretl_zero_matrix_new(n, n);

    if (H == NULL) {
        *err = E_ALLOC;
    } else {
        *err = hessian_from_score(b, H, gradfunc, cfunc, data);
    }

    if (!*err) {
        *err = gretl_invert_symmetric_matrix(H);
        if (*err) {
            fprintf(stderr, "hessian_inverse_from_score: failed\n");
            gretl_matrix_free(H);
            H = NULL;
        }
    }

    return H;
}

struct uhess_data {
    GENERATOR *genr;
    DATASET *dset;
};

/* Callback from numerical_hessian() for use in user_hess().
   The first argument is redundant here, since the user
   matrix referenced in the function-call in @genr will
   already contain the values in @b. But we need to
   respect the typedef for BFGS_CRIT_FUNC.
*/

static double uhess_callback (const double *b, void *data)
{
    struct uhess_data *uh = data;
    double ret = NADBL;
    int err;

    err = execute_genr(uh->genr, uh->dset, NULL);
    if (!err) {
        ret = get_scalar_value_by_name("$umax", &err);
    }

    return ret;
}

/* apparatus for constructing numerical approximation to
   the Hessian */

static void hess_h_init (double *h, double *h0, int n)
{
    memcpy(h, h0, n * sizeof *h);
}

static void hess_h_reduce (double *h, double v, int n)
{
    int i;

    for (i=0; i<n; i++) {
        h[i] /= v;
    }
}

static int numgrad_status;

int numgrad_in_progress (void)
{
    return numgrad_status;
}

static void set_numgrad_status (int s)
{
    numgrad_status = s;
}

/* number of Richardson steps */
#define RSTEPS 4

/* Default @d for numerical_hessian (2017-10-03: was 0.0001).
   Note 2018-10-05: this "new" value seems to be much too big
   in some cases.
*/
#define numhess_d 0.01

/* The algorithm below implements the method of Richardson
   Extrapolation.  It is derived from code in the gnu R package
   "numDeriv" by Paul Gilbert, which was in turn derived from code
   by Xinqiao Liu.  Turned into C and modified for gretl by
   Allin Cottrell, June 2006.  On successful completion, writes
   the Hessian (or the negative Hessian, if @neg is non-zero)
   into @H, which must be correctly sized to receive the result.
*/

int numerical_hessian (double *b, gretl_matrix *H,
                       BFGS_CRIT_FUNC func, void *data,
                       int neg, double d)
{
    double Dx[RSTEPS];
    double Hx[RSTEPS];
    double *wspace;
    double *h0, *h, *Hd, *D;
    int r = RSTEPS;
    double dsmall = 0.0001;
    double ztol, eps = 1e-4;
    double v = 2.0;    /* reduction factor for h */
    double f0, f1, f2;
    double p4m, hij;
    double bi0, bj0;
    int n = gretl_matrix_rows(H);
    int vn = (n * (n + 1)) / 2;
    int dn = vn + n;
    int i, j, k, m, u;
    int err = 0;

    if (d == 0.0) {
        d = numhess_d;
    }

    wspace = malloc((3 * n + dn) * sizeof *wspace);
    if (wspace == NULL) {
        return E_ALLOC;
    }

    h0 = wspace;
    h = h0 + n;
    Hd = h + n;
    D = Hd + n; /* D is of length dn */

#if 0
    ztol = sqrt(DBL_EPSILON / 7e-7); /* as per R */
#else
    ztol = 0.01;
#endif

    set_numgrad_status(1);

 try_again:

    /* note: numDeriv has

       h0 <- abs(d*x) + eps * (abs(x) < zero.tol)

       where the defaults are eps = 1e-4, d = 0.1,
       and zero.tol = sqrt(double.eps/7e-7)

       C translation:
       double ztol = sqrt(DBL_EPSILON / 7e-7);
       h0[i] = fabs(d*b[i]) + eps * (fabs(b[i]) < ztol);

       The above @ztol is about 1.78e-05. Below, we are
       currently using a bigger value, 0.01.
    */
    for (i=0; i<n; i++) {
        h0[i] = fabs(d*b[i]) + eps * (fabs(b[i]) < ztol);
    }

    f0 = func(b, data);

    /* first derivatives and Hessian diagonal */

    for (i=0; i<n; i++) {
        bi0 = b[i];
        hess_h_init(h, h0, n);
        for (k=0; k<r; k++) {
            b[i] = bi0 + h[i];
            f1 = func(b, data);
            if (na(f1)) {
                if (d <= dsmall) {
                    fprintf(stderr, "numerical_hessian: 1st derivative: "
                            "criterion=NA for theta[%d] = %g (d=%g)\n", i, b[i], d);
                }
                b[i] = bi0;
                err = E_NAN;
                goto end_first_try;
            }
            b[i] = bi0 - h[i];
            f2 = func(b, data);
            if (na(f2)) {
                if (d <= dsmall) {
                    fprintf(stderr, "numerical_hessian: 1st derivative: "
                            "criterion=NA for theta[%d] = %g (d=%g)\n", i, b[i], d);
                }
                b[i] = bi0;
                err = E_NAN;
                goto end_first_try;
            }
            /* F'(i) */
            Dx[k] = (f1 - f2) / (2 * h[i]);
            /* F''(i) */
            Hx[k] = (f1 - 2*f0 + f2) / (h[i] * h[i]);
            hess_h_reduce(h, v, n);
        }
        b[i] = bi0;
        p4m = 4.0;
        for (m=0; m<r-1; m++) {
            for (k=0; k<r-m-1; k++) {
                Dx[k] = (Dx[k+1] * p4m - Dx[k]) / (p4m - 1);
                Hx[k] = (Hx[k+1] * p4m - Hx[k]) / (p4m - 1);
            }
            p4m *= 4.0;
        }
        D[i] = Dx[0];
        Hd[i] = Hx[0];
    }

    /* second derivatives: lower half of Hessian only */

    u = n;
    for (i=0; i<n; i++) {
        bi0 = b[i];
        for (j=0; j<=i; j++) {
            if (i == j) {
                D[u] = Hd[i];
            } else {
                hess_h_init(h, h0, n);
                bj0 = b[j];
                for (k=0; k<r; k++) {
                    b[i] = bi0 + h[i];
                    b[j] = bj0 + h[j];
                    f1 = func(b, data);
                    if (na(f1)) {
                        if (d <= dsmall) {
                            fprintf(stderr, "numerical_hessian: 2nd derivatives (%d,%d): "
                                    "objective function gave NA\n", i, j);
                        }
                        b[i] = bi0;
                        b[j] = bj0;
                        err = E_NAN;
                        goto end_first_try;
                    }
                    b[i] = bi0 - h[i];
                    b[j] = bj0 - h[j];
                    f2 = func(b, data);
                    if (na(f2)) {
                        if (d <= dsmall) {
                            fprintf(stderr, "numerical_hessian: 2nd derivatives (%d,%d): "
                                    "objective function gave NA\n", i, j);
                        }
                        b[i] = bi0;
                        b[j] = bj0;
                        err = E_NAN;
                        goto end_first_try;
                    }
                    /* cross-partial */
                    Dx[k] = (f1 - 2*f0 + f2 - Hd[i]*h[i]*h[i]
                             - Hd[j]*h[j]*h[j]) / (2*h[i]*h[j]);
                    hess_h_reduce(h, v, n);
                }
                p4m = 4.0;
                for (m=0; m<r-1; m++) {
                    for (k=0; k<r-m-1; k++) {
                        Dx[k] = (Dx[k+1] * p4m - Dx[k]) / (p4m - 1);
                    }
                    p4m *= 4.0;
                }
                D[u] = Dx[0];
                b[j] = bj0;
            }
            u++;
        }
        b[i] = bi0;
    }

 end_first_try:
    if (err == E_NAN && d > dsmall) {
        err = 0;
        gretl_error_clear();
        d /= 10;
        goto try_again;
    }

    if (!err) {
        /* transcribe the (negative of?) the Hessian */
        u = n;
        for (i=0; i<n; i++) {
            for (j=0; j<=i; j++) {
                hij = neg ? -D[u] : D[u];
                gretl_matrix_set(H, i, j, hij);
                gretl_matrix_set(H, j, i, hij);
                u++;
            }
        }
    }

    if (neg) {
        /* internal use, not user_numhess(): we should ensure
           that the original criterion value is restored after
           calculating with a perturbed version of @b
        */
        func(b, data);
    }

    if (err && err != E_ALLOC) {
        gretl_errmsg_set(_("Failed to compute numerical Hessian"));
    }

    set_numgrad_status(0);

    free(wspace);

    return err;
}

/**
 * numerical_hessian_inverse:
 * @b: array of parameter estimates.
 * @n: the number of elements in @b.
 * @func: function to compute criterion.
 * @data: data to be passed to the @gradfunc callback.
 * @d: step size (give 0.0 for automatic).
 * @err: location to receive error code.
 *
 * A wrapper for numerical_hessian() which takes care of
 * (a) allocation of the Hessian and (b) inversion.
 *
 * Returns: the inverse of the (negative) Hessian on successful
 * completion, NULL on error.
 */

gretl_matrix *numerical_hessian_inverse (const double *b, int n,
                                         BFGS_CRIT_FUNC func,
                                         void *data, double d,
                                         int *err)
{
    gretl_matrix *H = gretl_zero_matrix_new(n, n);

    if (H == NULL) {
        *err = E_ALLOC;
    } else {
        *err = numerical_hessian((double *) b, H, func, data, 1, d);
    }

    if (!*err) {
        *err = gretl_invert_symmetric_matrix(H);
        if (*err) {
            fprintf(stderr, "numerical_hessian_inverse: failed\n");
            gretl_errmsg_set(_("Failed to compute numerical Hessian"));
            gretl_matrix_free(H);
            H = NULL;
        }
    }

    return H;
}

static int NR_fallback_hessian (double *b, gretl_matrix *H,
                                BFGS_GRAD_FUNC gradfunc,
                                BFGS_CRIT_FUNC cfunc,
                                void *data)
{
    if (gradfunc != NULL) {
        return hessian_from_score(b, H, gradfunc, cfunc, data);
    } else {
        return numerical_hessian(b, H, cfunc, data, 1, 0.0);
    }
}

/* In the case of ARMA with missing values we can end up with
   a certain number of trailing NAs in the score calculation.
   This function checks for that condition and if necessary
   trims off the rows in question.
*/

static gretl_matrix *maybe_trim_score (gretl_matrix *G)
{
    gretl_matrix *ret = NULL;
    int T = G->rows;
    int okrows = T;
    int i, t, nancount;

    for (t=T-1; t>0; t--) {
        nancount = 0;
        for (i=0; i<G->cols; i++) {
            if (isnan(gretl_matrix_get(G, t, i))) {
                nancount++;
            }
        }
        if (nancount == G->cols) {
            okrows--;
        } else {
            break;
        }
    }

    if (okrows < T) {
        ret = gretl_matrix_alloc(okrows, G->cols);

        if (ret != NULL) {
            size_t csize = okrows * sizeof(double);
            double *dest = ret->val;
            double *src = G->val;

            for (i=0; i<G->cols; i++) {
                memcpy(dest, src, csize);
                dest += okrows;
                src += G->rows;
            }
            gretl_matrix_free(G);
        }
    } else {
        ret = G;
    }

    return ret;
}

#define ALT_OPG 0

/* build the T x k matrix G, given a set of coefficient estimates,
   @b, and a function for calculating the per-observation contributions
   to the loglikelihood, @lltfun
*/

gretl_matrix *numerical_score_matrix (double *b, int T, int k,
                                      BFGS_LLT_FUNC lltfun,
                                      void *data, int *err)
{
    double h = 1e-8;
#if ALT_OPG
    double d = 1.0e-4;
#endif
    gretl_matrix *G;
    const double *x;
    double bi0, x0;
    int i, t;

    G = gretl_zero_matrix_new(T, k);
    if (G == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    set_numgrad_status(1);

    for (i=0; i<k; i++) {
        bi0 = b[i];
#if ALT_OPG
        h = d * bi0 + d * (floateq(b[i], 0.0));
#endif
        b[i] = bi0 - h;
        x = lltfun(b, i, data);
        if (x == NULL) {
            *err = E_NAN;
            goto bailout;
        }
        for (t=0; t<T; t++) {
            gretl_matrix_set(G, t, i, x[t]);
        }
        b[i] = bi0 + h;
        x = lltfun(b, i, data);
        if (x == NULL) {
            *err = E_NAN;
            goto bailout;
        }
        for (t=0; t<T; t++) {
            x0 = gretl_matrix_get(G, t, i);
            gretl_matrix_set(G, t, i, (x[t] - x0) / (2.0 * h));
        }
        b[i] = bi0;
#if NLS_DEBUG
        fprintf(stderr, "b[%d]: using %#.12g and %#.12g\n", i, bi0 - h, bi0 + h);
#endif
    }

    set_numgrad_status(0);

    /* trim missing values? */
    G = maybe_trim_score(G);

#if NLS_DEBUG
    gretl_matrix_print(G, "Numerically estimated score");
#endif

 bailout:

    if (*err) {
        gretl_matrix_free(G);
        G = NULL;
    }

    return G;
}

static int richardson_gradient (double *b, double *g, int n,
                                BFGS_CRIT_FUNC func, void *data)
{
    double df[RSTEPS];
    double eps = 1.0e-4;
    double d = 0.0001;
    double h, p4m;
    double bi0, f1, f2;
    int r = RSTEPS;
    int i, k, m;

    for (i=0; i<n; i++) {
        bi0 = b[i];
        h = fabs(d * b[i]) + eps * (floateq(b[i], 0.0));
        for (k=0; k<r; k++) {
            b[i] = bi0 - h;
            f1 = func(b, data);
            b[i] = bi0 + h;
            f2 = func(b, data);
            if (na(f1) || na(f2)) {
                b[i] = bi0;
                return 1;
            }
            df[k] = (f2 - f1) / (2 * h);
            h /= 2.0;
        }
        b[i] = bi0;
        p4m = 4.0;
        for (m=0; m<r-1; m++) {
            for (k=0; k<r-m-1; k++) {
                df[k] = (df[k+1] * p4m - df[k]) / (p4m - 1.0);
            }
            p4m *= 4.0;
        }
        g[i] = df[0];
    }

    return 0;
}

/* trigger for switch to Richardson gradient */
#define B_RELMIN 1.0e-14

static int simple_gradient (double *b, double *g, int n,
                            BFGS_CRIT_FUNC func, void *data,
                            int *redo)
{
    const double h = 1.0e-8;
    double bi0, f1, f2;
    int i;

    for (i=0; i<n; i++) {
        bi0 = b[i];
        b[i] = bi0 - h;
        if (bi0 != 0.0 && fabs((bi0 - b[i]) / bi0) < B_RELMIN) {
            fprintf(stderr, "numerical gradient: switching to Richardson\n");
            *redo = 1;
            return 0;
        }
        f1 = func(b, data);
        b[i] = bi0 + h;
        f2 = func(b, data);
        b[i] = bi0;
        if (na(f1) || na(f2)) {
            return 1;
        }
        g[i] = (f2 - f1) / (2.0 * h);
#if BFGS_DEBUG > 1
        fprintf(stderr, "g[%d] = (%.16g - %.16g) / (2.0 * %g) = %g\n",
                i, f2, f1, h, g[i]);
#endif
    }

    return 0;
}

/* default numerical calculation of gradient in context of BFGS */

int numeric_gradient (double *b, double *g, int n,
                      BFGS_CRIT_FUNC func, void *data)
{
    int err = 0;

    if (libset_get_bool(BFGS_RSTEP)) {
        err = richardson_gradient(b, g, n, func, data);
    } else {
        int redo = 0;

        err = simple_gradient(b, g, n, func, data, &redo);
        if (redo) {
            err = richardson_gradient(b, g, n, func, data);
        }
    }

#if BFGS_DEBUG
    fprintf(stderr, "numeric_gradient returning, err = %d\n", err);
#endif

    return err;
}

#define STEPFRAC        0.2
#define acctol          1.0e-7 /* alt: 0.0001 or 1.0e-7 (?) */
#define reltest         10.0

#define coeff_unchanged(a,b) (reltest + a == reltest + b)

static int broken_gradient (double *g, int n)
{
    int i;

    for (i=0; i<n; i++) {
        if (isnan(g[i])) {
            return 1;
        }
    }

    return 0;
}

/* If "set initvals" has been used, replace whatever initial values
   might have been in place with those given by the user (the customer
   is always right).  In addition, respect user settings for the maximum
   number of iterations, the convergence tolerance and so on.
*/

static void optim_get_user_values (double *b, int n, int *maxit,
                                   double *reltol, double *gradmax,
                                   gretlopt opt, PRN *prn)
{
    int umaxit;
    double utol;

    if (opt & OPT_U) {
        /* We first check to see if we've been a usable initialization
           for the parameter estimates.
        */
        gretl_matrix *uinit;
        int i, uilen;

        uinit = get_initvals();
        uilen = gretl_vector_get_length(uinit);

        if (uilen > 0) {
            /* the user has given something */
            if (uilen < n) {
                fprintf(stderr, "Only %d initial values given, but %d "
                        "are necessary\n", uilen, n);
            } else {
                for (i=0; i<n; i++) {
                    b[i] = uinit->val[i];
                }
                if ((opt & OPT_V) && !(opt & OPT_A)) {
                    /* OPT_A: arma: this is handled elsewhere */
                    pputs(prn, _("\n\n*** User-specified starting values:\n"));
                    for (i=0; i<n; i++) {
                        pprintf(prn, " %12.6f", b[i]);
                        if (i % 6 == 5) {
                            pputc(prn, '\n');
                        }
                    }
                    pputs(prn, "\n\n");
                }
            }
        }
        gretl_matrix_free(uinit);
    }

    if (reltol == NULL || gradmax == NULL) {
        /* Newton */
        return;
    }

    /* check for a setting of the maximum number of iterations */
    umaxit = libset_get_int(BFGS_MAXITER);
    if (umaxit >= 0) {
        *maxit = umaxit;
    } else if (*maxit < 0) {
        *maxit = BFGS_MAXITER_DEFAULT;
    }

    /* convergence tolerance */
    utol = libset_get_user_tolerance(BFGS_TOLER);
    if (!na(utol)) {
        /* the user has actually set a value */
        *reltol = utol;
        if (!(opt & OPT_Q)) {
            fprintf(stderr, "user-specified BFGS tolerance = %g\n", utol);
        }
    } else if (*reltol == 0) {
        /* use the generic BFGS default */
        *reltol = libset_get_double(BFGS_TOLER);
    }

    /* maximum acceptable gradient norm */
    *gradmax = libset_get_double(BFGS_MAXGRAD);
}

#define bfgs_print_iter(v,s,i) (v && (s == 1 || i % s == 0))

#define GRAD_TOLER 1.0

static int copy_initial_hessian (double **H,
                                 const gretl_matrix *A,
                                 int n)
{
    int i, j, vlen = gretl_vector_get_length(A);
    int err = 0;

#if BFGS_DEBUG > 1
    gretl_matrix_print(A, "BFGS: initial Hessian inverse");
#endif

    if (gretl_is_null_matrix(A)) {
        /* set identity matrix */
        for (i=0; i<n; i++) {
            for (j=0; j<i; j++) {
                H[i][j] = 0.0;
            }
            H[i][i] = 1.0;
        }
    } else if (vlen == n) {
        /* set the diagonal */
        for (i=0; i<n; i++) {
            for (j=0; j<i; j++) {
                H[i][j] = 0.0;
            }
            H[i][i] = A->val[i];
        }
    } else  if (A->rows == n && A->cols == n) {
        /* set the whole matrix */
        for (i=0; i<n; i++) {
            for (j=0; j<=i; j++) {
                H[i][j] = gretl_matrix_get(A, i, j);
            }
        }
    } else {
        err = E_NONCONF;
    }

    return err;
}

static double optim_fncall (BFGS_CRIT_FUNC cfunc,
                            double *b, void *data,
                            int minimize)
{
    double ret = cfunc(b, data);

    return na(ret) ? ret : minimize ? -ret : ret;
}

static int optim_gradcall (BFGS_GRAD_FUNC gradfunc,
                           double *b, double *g, int n,
                           BFGS_CRIT_FUNC cfunc,
                           void *data,
                           int minimize)
{
    int ret = gradfunc(b, g, n, cfunc, data);

    if (minimize) {
        int i;

        for (i=0; i<n; i++) {
            if (!na(g[i])) {
                g[i] = -g[i];
            }
        }
    }

    return ret;
}

static double simple_slen (int n, int *pndelta, double *b, double *X, double *t,
                           double *pf, BFGS_CRIT_FUNC cfunc, void *data,
                           double g0, double f0, int *pfcount, int minimize)
{
    double d, steplen = 1.0;
    double f1 = *pf;
    int i, crit_ok = 0, fcount = 0;
    int ndelta = *pndelta;

    /* Below: iterate so long as (a) we haven't achieved an acceptable
       value of the criterion and (b) there is still some prospect
       of doing so.
    */

    do {
        ndelta = n;
        crit_ok = 0;
        for (i=0; i<n; i++) {
            b[i] = X[i] + steplen * t[i];
            if (coeff_unchanged(b[i], X[i])) {
                ndelta--;
            }
        }
        if (ndelta > 0) {
            f1 = optim_fncall(cfunc, b, data, minimize);
            fcount++;
            d = g0 * steplen * acctol;
            crit_ok = !na(f1) && (f1 >= f0 + d);
            if (!crit_ok) {
                /* calculated criterion no good: try smaller step */
                steplen *= STEPFRAC;
            }
        }
    } while (ndelta != 0 && !crit_ok);

    *pndelta = ndelta;
    *pfcount += fcount;
    *pf = f1;

    return steplen;
}

#define GRADCHECK 0 /* 2023-10-30: for testing */

static int BFGS_orig (double *b, int n, int maxit, double reltol,
                      int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc,
                      int crittype, BFGS_GRAD_FUNC gradfunc, void *data,
                      const gretl_matrix *A0, gretlopt opt, PRN *prn)
{
    int verbskip, verbose = (opt & OPT_V);
    int minimize = (opt & OPT_I);
    double *wspace = NULL;
    double **H = NULL;
    double *g, *t, *X, *c;
    int fcount, gcount, ndelta = 0;
    int show_activity = 0;
    double sumgrad, gradmax, gradnorm = 0.0;
#if GRADCHECK
    double L2;
#endif
    double fmax, f, f0, s, steplen = 0.0;
    double fdiff, D1, D2;
    int i, j, ilast, iter, done = 0;
    int err = 0;

    optim_get_user_values(b, n, &maxit, &reltol, &gradmax, opt, prn);

    wspace = malloc(4 * n * sizeof *wspace);
    H = triangular_array_new(n);
    if (wspace == NULL || H == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    if (gradfunc == NULL) {
        gradfunc = numeric_gradient;
    }

    /* initialize curvature matrix */
    if (A0 != NULL) {
        err = copy_initial_hessian(H, A0, n);
    } else {
        gretl_matrix *A1 = get_initcurv();

        err = copy_initial_hessian(H, A1, n);
        gretl_matrix_free(A1);
    }
    if (err) {
        goto bailout;
    }

    g = wspace;
    t = g + n;
    X = t + n;
    c = X + n;

    f = optim_fncall(cfunc, b, data, minimize);

    if (na(f)) {
        gretl_errmsg_set(_("BFGS: initial value of objective function is not finite"));
        err = E_NAN;
        goto bailout;
    }

#if BFGS_DEBUG
    fprintf(stderr, "*** BFGS: first evaluation of f = %g\n", f);
#endif

    f0 = fmax = f;
    iter = ilast = fcount = gcount = 1;
    optim_gradcall(gradfunc, b, g, n, cfunc, data, minimize);

#if BFGS_DEBUG > 1
    fprintf(stderr, "initial gradient:\n");
    for (i=0; i<n; i++) {
        fprintf(stderr, " g[%d] = %g\n", i, g[i]);
    }
#endif

    if (maxit == 0) {
        goto skipcalc;
    }

    verbskip = libset_get_int(BFGS_VERBSKIP);
    show_activity = show_activity_func_installed();

    do {
        if (get_user_stop()) {
            err = E_STOP;
            break;
        }
        if (bfgs_print_iter(verbose, verbskip, iter)) {
            print_iter_info(iter, f, crittype, n, b, g, steplen, prn);
        }
        if (show_activity && (iter % 10 == 0)) {
            show_activity_callback();
        }

        if (iter > 1 && ilast == gcount) {
            /* restart: set curvature matrix to I */
            for (i=0; i<n; i++) {
                for (j=0; j<i; j++) {
                    H[i][j] = 0.0;
                }
                H[i][i] = 1.0;
            }
        }

        for (i=0; i<n; i++) {
            /* copy coefficients to X, gradient to c */
            X[i] = b[i];
            c[i] = g[i];
        }

        gradnorm = sumgrad = 0.0;
#if GRADCHECK
	L2 = 0.0;
#endif

        for (i=0; i<n; i++) {
            s = 0.0;
            for (j=0; j<=i; j++) {
                s += H[i][j] * g[j];
            }
            for (j=i+1; j<n; j++) {
                s += H[j][i] * g[j];
            }
            t[i] = s;
            sumgrad += s * g[i];
            gradnorm += fabs(b[i] * g[i]);
#if GRADCHECK
	    L2 += g[i] * g[i];
#endif
        }

        gradnorm = sqrt(gradnorm / n);

#if BFGS_DEBUG
        fprintf(stderr, "\niter %d: sumgrad=%g, gradnorm=%g\n",
                iter, sumgrad, gradnorm);
# if BFGS_DEBUG > 1
        fprintf(stderr, "H = \n");
        for (i=0; i<n; i++) {
            for (j=0; j<=i; j++) {
                fprintf(stderr, "%15.6f", H[i][j]);
            }
            fputc('\n', stderr);
        }
# endif
#endif
        if (sumgrad > 0.0) {
            /* heading in the right direction (uphill) */
            steplen = simple_slen(n, &ndelta, b, X, t, &f, cfunc, data,
                                  sumgrad, fmax, &fcount, minimize);
            fdiff = fabs(fmax - f);
            if (iter > 1 || fdiff > 0) {
                done = fdiff <= reltol * (fabs(fmax) + reltol);
#if BFGS_DEBUG
                fprintf(stderr, "convergence test: LHS=%g, RHS=%g; done=%d\n",
                        fdiff, reltol * (fabs(fmax) + reltol), done);
#endif
            }

            /* prepare to stop if relative change is small enough */
            if (done) {
                ndelta = 0;
                fmax = f;
            }

            if (ndelta > 0) {
                /* making progress */
#if BFGS_DEBUG
                fprintf(stderr, "making progress, ndelta = %d\n", ndelta);
#endif
                fmax = f;
                optim_gradcall(gradfunc, b, g, n, cfunc, data, minimize);
#if BFGS_DEBUG > 1
                fprintf(stderr, "new gradient:\n");
                for (i=0; i<n; i++) {
                    fprintf(stderr, "%15.6f", g[i]);
                }
                fputc('\n', stderr);
#endif
                gcount++;
                iter++;
                D1 = 0.0;
                for (i=0; i<n; i++) {
                    t[i] *= steplen;
                    c[i] -= g[i];
                    D1 += t[i] * c[i];
                }
                if (D1 > 0.0) {
                    D2 = 0.0;
                    for (i=0; i<n; i++) {
                        s = 0.0;
                        for (j=0; j<=i; j++) {
                            s += H[i][j] * c[j];
                        }
                        for (j=i+1; j<n; j++) {
                            s += H[j][i] * c[j];
                        }
                        X[i] = s;
                        D2 += s * c[i];
                    }
                    D2 = 1.0 + D2 / D1;
                    for (i=0; i<n; i++) {
                        for (j=0; j<=i; j++) {
                            H[i][j] += (D2 * t[i]*t[j] - X[i]*t[j] - t[i]*X[j]) / D1;
                        }
                    }
                } else {
                    /* D1 <= 0.0 */
                    ilast = gcount;
                }
            } else if (ilast < gcount) {
                ndelta = n;
                ilast = gcount;
            }
        } else if (sumgrad == 0.0) {
            break;
        } else {
            /* heading in the wrong direction */
            if (ilast == gcount) {
                /* we just did a reset, so don't reset again; instead set
                   ndelta = 0 so that we exit the main loop
                */
                ndelta = 0;
                if (gcount == 1) {
                    err = (broken_gradient(g, n))? E_NAN : E_NOCONV;
                }
            } else {
                /* reset for another attempt */
                ilast = gcount;
                ndelta = n;
            }
        }

        if (iter >= maxit) {
            break;
        }

        if (gcount - ilast > 2 * n) {
            /* periodic restart of curvature computation */
            ilast = gcount;
        }

    } while (ndelta > 0 || ilast < gcount); /* end of outer do-loop */

    if (err == E_STOP) {
        goto bailout;
    }

#if BFGS_DEBUG
    fprintf(stderr, "terminated: fmax=%g, ndelta=%d, ilast=%d, gcount=%d\n",
            fmax, ndelta, ilast, gcount);
    fprintf(stderr, "gradnorm = %g, vs gradmax = %g\n", gradnorm, gradmax);
#endif
#if GRADCHECK
    fprintf(stderr, "gradnorm = %g, L2 = %g\n", gradnorm, sqrt(L2));
#endif

    if (iter >= maxit) {
        gretl_errmsg_sprintf(_("Reached the maximum iterations (%d)"), maxit);
        err = E_NOCONV;
    } else if (gradnorm > gradmax) {
        gretl_errmsg_sprintf(_("Norm of gradient %g exceeds maximum of %g"),
                             gradnorm, gradmax);
        err = E_NOCONV;
    } else if (fmax < f0) {
        /* allow a small slop factor here? */
        double rdiff = fabs(f0 - fmax);

        if (fabs(f0) > 1.0e-16) {
            rdiff /= fabs(f0);
        }
        if (rdiff > 1.0e-12) {
            fprintf(stderr, "failed to match initial value of objective function:\n"
                    " f0=%.16g, fmax=%.16g, rdiff=%g\n", f0, fmax, rdiff);
            err = E_NOCONV;
        }
    }

    if (!err && gradnorm > GRAD_TOLER && gretl_warnings_on()) {
        gretl_warnmsg_sprintf(_("norm of gradient = %g"), gradnorm);
        set_gretl_warning(W_GRADIENT);
    }

 skipcalc:

    /* particularly relevant for iterated GMM: increment,
       don't just set, these two counts
    */
    *fncount += fcount;
    *grcount += gcount;

    if (verbose) {
        print_iter_info(-1, f, crittype, n, b, g, steplen, prn);
        /* pputc(prn, '\n'); */
    }

 bailout:

    free(wspace);
    free_triangular_array(H, n);

#if BFGS_DEBUG
    fprintf(stderr, "BFGS_max: returning %d\n", err);
#endif

    return err;
}

/* Note: we need this because the original L-BFGS-B code is
   set up as a minimizer.  We could get rid of it if anyone
   has the strength to go into lbfgsb.c (ex Fortran) and make
   the necessary adjustments.
*/

static void reverse_gradient (double *g, int n)
{
    int i;

    for (i=0; i<n; i++) {
        if (!na(g[i])) {
            g[i] = -g[i];
        }
    }
}

static int transcribe_lbfgs_bounds (const gretl_matrix *m,
                                    int nparm, int *nbd,
                                    double *l, double *u)
{
    double h = libset_get_double(CONV_HUGE);
    int i, j, c = gretl_matrix_cols(m);
    int err = 0;

    if (c != 3) {
        fprintf(stderr, "lbfgs_bounds: matrix should have 3 cols\n");
        return E_INVARG;
    }

    for (i=0; i<nparm; i++) {
        /* mark as unbounded */
        nbd[i] = 0;
    }

    for (i=0; i<m->rows && !err; i++) {
        j = (int) gretl_matrix_get(m, i, 0);
        if (j < 1 || j > nparm) {
            fprintf(stderr, "lbfgs_bounds: out-of-bounds index %d\n", j);
            err = E_INVARG;
        } else {
            j--; /* convert to zero-based */
            l[j] = gretl_matrix_get(m, i, 1);
            u[j] = gretl_matrix_get(m, i, 2);
            if (l[j] > u[j]) {
                err = E_INVARG;
            } else if (l[j] != -h && u[j] != h) {
                /* both lower and upper bounds */
                nbd[j] = 2;
            } else if (l[j] != -h) {
                /* lower bound only */
                nbd[j] = 1;
            } else if (u[j] != h) {
                /* upper bound only */
                nbd[j] = 3;
            }
        }
    }

    return err;
}

int LBFGS_max (double *b, int n,
               int maxit, double reltol,
               int *fncount, int *grcount,
               BFGS_CRIT_FUNC cfunc, int crittype,
               BFGS_GRAD_FUNC gradfunc,
               BFGS_COMBO_FUNC combfunc,
               void *data,
               const gretl_matrix *bounds,
               gretlopt opt,
               PRN *prn)
{
    double *wspace = NULL;
    int *ispace = NULL;
    double *g, *l, *u, *wa;
    int *iwa, *nbd;
    int i, m, dim;
    char task[60];
    char csave[60];
    double f, pgtol;
    double factr;
    double gradmax;
    double dsave[29];
    int isave[44];
    int lsave[4];
    int iter, ibak = 0;
    int show_activity = 0;
    int maximize;
    int verbskip, verbose = (opt & OPT_V);
    int err = 0;

    maximize = (crittype != C_SSR) && !(opt & OPT_I);

    *fncount = *grcount = 0;

    optim_get_user_values(b, n, &maxit, &reltol, &gradmax, opt, prn);

    /*
      m: the number of corrections used in the limited memory matrix.
      It is not altered by the routine.  Values of m < 3 are not
      recommended, and large values of m can result in excessive
      computing time. The range 3 <= m <= 20 is recommended.

      Was initially set to 5 (then 10, then 8; and 8 is the default).
    */
    m = libset_get_int(LBFGS_MEM);

    dim = (2*m+5)*n + 11*m*m + 8*m; /* for wa */
    dim += 3*n;                     /* for g, l and u */

    wspace = malloc(dim * sizeof *wspace);
    ispace = malloc(4*n * sizeof *ispace);

    if (wspace == NULL || ispace == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    g = wspace;
    l = g + n;
    u = l + n;
    wa = u + n;

    nbd = ispace;
    iwa = nbd + n;

    verbskip = libset_get_int(BFGS_VERBSKIP);
    show_activity = show_activity_func_installed();

    if (gradfunc == NULL && combfunc == NULL) {
        gradfunc = numeric_gradient;
    }

    /* Gradient convergence criterion (currently unused --
       we use reltol instead) */
    pgtol = 0.0;

    /* tol = (factr * macheps) => factr = tol/macheps */
    factr = reltol / pow(2.0, -52);

    if (!gretl_is_null_matrix(bounds)) {
        /* Handle specified bounds on the parameters */
        err = transcribe_lbfgs_bounds(bounds, n, nbd, l, u);
        if (err) {
            goto bailout;
        }
    } else {
        /* By default we just set all parameters to be
           less than some ridiculously large number */
        for (i=0; i<n; i++) {
            nbd[i] = 3; /* case 3: upper bound only */
            u[i] = DBL_MAX / 100;
        }
    }

    /* Start the iteration by initializing @task */
    strcpy(task, "START");

    while (1) {
        if (get_user_stop()) {
            err = E_STOP;
            break;
        }

        /* Call the L-BFGS-B code */
        setulb_(&n, &m, b, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa,
                task, csave, lsave, isave, dsave);

        iter = isave[29] + 1;

        if (!strncmp(task, "FG", 2)) {
            /* Compute function value, f */
            if (combfunc != NULL) {
                f = combfunc(b, g, n, data);
            } else {
                f = cfunc(b, data);
            }
            if (!na(f)) {
                if (maximize) f = -f;
            } else if (*fncount == 0) {
                fprintf(stderr, "initial value of f is not finite\n");
                err = E_DATA;
                break;
            }
            *fncount += 1;
            if (combfunc == NULL) {
                /* Compute gradient, g */
                gradfunc(b, g, n, cfunc, data);
            }
            if (maximize) {
                reverse_gradient(g, n);
            }
            *grcount += 1;
        } else if (!strncmp(task, "NEW_X", 5)) {
            /* The optimizer has produced a new set of parameter values */
            if (isave[33] >= maxit) {
                strcpy(task, "STOP: TOTAL NO. of f AND g "
                       "EVALUATIONS EXCEEDS LIMIT");
                err = E_NOCONV;
                break;
            }
        } else {
            if (strncmp(task, "CONVER", 6)) {
                fprintf(stderr, "%s\n", task);
            }
            break;
        }

        if (bfgs_print_iter(verbose, verbskip, iter)) {
            if (iter != ibak) {
                double steplen = (iter == 1)? NADBL : dsave[13];

                if (maximize) reverse_gradient(g, n);
                print_iter_info(iter, -f, crittype, n, b, g, steplen, prn);
                if (maximize) reverse_gradient(g, n);
            }
            ibak = iter;
        }

        if (show_activity && (iter % 10 == 0)) {
            show_activity_callback();
        }
    }

    if (err == E_STOP) {
        goto bailout;
    }

    if (!err && crittype == C_GMM) {
        /* finalize GMM computations */
        f = cfunc(b, data);
    }

    if (opt & OPT_V) {
        if (maximize) reverse_gradient(g, n);
        print_iter_info(-1, -f, crittype, n, b, g, dsave[13], prn);
        pputc(prn, '\n');
    }

 bailout:

    free(wspace);
    free(ispace);

    return err;
}

/**
 * BFGS_max:
 * @b: array of adjustable coefficients.
 * @n: number elements in array @b.
 * @maxit: the maximum number of iterations to allow.
 * @reltol: relative tolerance for terminating iteration.
 * @fncount: location to receive count of function evaluations.
 * @grcount: location to receive count of gradient evaluations.
 * @cfunc: pointer to function used to calculate maximand.
 * @crittype: code for type of the maximand/minimand: should
 * be %C_LOGLIK, %C_GMM or %C_OTHER.  Used only in printing
 * iteration info.
 * @gradfunc: pointer to function used to calculate the
 * gradient, or %NULL for default numerical calculation.
 * @data: pointer that will be passed as the last
 * parameter to the callback functions @cfunc and @gradfunc.
 * @A0: initial approximation to the inverse of the Hessian
 * (or %NULL to use identity matrix)
 * @opt: may contain %OPT_V for verbose operation, %OPT_L to
 * force use of L-BFGS-B.
 * @prn: printing struct (or %NULL).  Only used if @opt
 * includes %OPT_V.
 *
 * Obtains the set of values for @b which jointly maximize the
 * criterion value as calculated by @cfunc. By default uses the BFGS
 * variable-metric method (based on Pascal code in J. C. Nash,
 * "Compact Numerical Methods for Computers," 2nd edition, converted
 * by p2c then re-crafted by B. D. Ripley for gnu R; revised for
 * gretl by Allin Cottrell and Jack Lucchetti). Alternatively,
 * if OPT_L is given, uses the L-BFGS-B method (limited memory
 * BFGS), based on Lbfgsb.3.0 by Ciyou Zhu, Richard Byrd, Jorge
 * Nocedal and Jose Luis Morales.
 *
 * Returns: 0 on successful completion, non-zero error code
 * on error.
 */

int BFGS_max (double *b, int n, int maxit, double reltol,
              int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc,
              int crittype, BFGS_GRAD_FUNC gradfunc, void *data,
              const gretl_matrix *A0, gretlopt opt, PRN *prn)
{
    int ret, wnum;

    gretl_iteration_push();

    if ((opt & OPT_L) || libset_get_bool(USE_LBFGS)) {
        ret = LBFGS_max(b, n, maxit, reltol,
                        fncount, grcount, cfunc,
                        crittype, gradfunc, NULL, data,
                        NULL, opt, prn);
    } else {
        ret = BFGS_orig(b, n, maxit, reltol,
                        fncount, grcount, cfunc,
                        crittype, gradfunc, data,
                        A0, opt, prn);
    }

    gretl_iteration_pop();

    wnum = check_gretl_warning();
    if (wnum != W_GRADIENT) {
        /* suppress expected numerical warnings */
        set_gretl_warning(0);
    }

    return ret;
}

/* constrained L_BFGS-B */

static int BFGS_cmax (double *b, int n,
                      int maxit, double reltol,
                      int *fncount, int *grcount,
                      BFGS_CRIT_FUNC cfunc, int crittype,
                      BFGS_GRAD_FUNC gradfunc,
                      void *data,
                      const gretl_matrix *bounds,
                      gretlopt opt, PRN *prn)
{
    int ret, wnum;

    gretl_iteration_push();

    ret = LBFGS_max(b, n, maxit, reltol,
                    fncount, grcount, cfunc,
                    crittype, gradfunc, NULL, data,
                    bounds, opt, prn);

    gretl_iteration_pop();

    wnum = check_gretl_warning();
    if (wnum != W_GRADIENT) {
        /* suppress expected numerical warnings */
        set_gretl_warning(0);
    }

    return ret;
}

/* user-level access to BFGS, NR and fdjac */

typedef struct umax_ umax;

struct umax_ {
    GretlType gentype;    /* GRETL_TYPE_DOUBLE or GRETL_TYPE_MATRIX */
    char *pxname;         /* name of scalar parameter value */
    gretl_matrix *b;      /* parameter vector */
    gretl_matrix *g;      /* gradient vector */
    gretl_matrix *h;      /* hessian matrix */
    int ncoeff;           /* number of coefficients */
    GENERATOR *gf;        /* for generating scalar or matrix result */
    GENERATOR *gg;        /* for generating gradient */
    GENERATOR *gh;        /* for generating Hessian */
    double fx_out;        /* function double value */
    gretl_matrix *fm_out; /* function matrix value */
    char gmname[VNAMELEN]; /* name of user-defined gradient vector */
    char hmname[VNAMELEN]; /* name of user-defined Hessian matrix */
    DATASET *dset;        /* dataset */
    PRN *prn;             /* optional printing struct */
};

static umax *umax_new (GretlType t)
{
    umax *u = malloc(sizeof *u);

    if (u != NULL) {
        u->gentype = t;
        u->pxname = NULL;
        u->b = NULL;
        u->g = NULL;
        u->h = NULL;
        u->ncoeff = 0;
        u->gf = NULL;
        u->gg = NULL;
        u->gh = NULL;
        u->fx_out = NADBL;
        u->fm_out = NULL;
        u->gmname[0] = '\0';
        u->hmname[0] = '\0';
        u->dset = NULL;
        u->prn = NULL;
    }

    return u;
}

static void umax_unset_no_replace_flag (umax *u)
{
    user_var *uv = get_user_var_by_data(u->b);

    if (uv != NULL) {
        user_var_unset_flag(uv, UV_NOREPL);
    }
}

static void umax_destroy (umax *u)
{
    if (u->dset != NULL) {
        /* drop any private "$" series created */
        dataset_drop_listed_variables(NULL, u->dset, NULL, NULL);
    }

    umax_unset_no_replace_flag(u);
    user_var_delete_by_name("$umax", NULL);
    free(u->pxname);

    destroy_genr(u->gf);
    destroy_genr(u->gg);
    destroy_genr(u->gh);

    free(u);
}

/* user-defined optimizer: get the criterion value */

static double user_get_criterion (const double *b, void *p)
{
    umax *u = (umax *) p;
    double x = NADBL;
    int i, t, err;

    for (i=0; i<u->ncoeff; i++) {
        u->b->val[i] = b[i];
    }

    err = execute_genr(u->gf, u->dset, u->prn);

    if (err) {
        return NADBL;
    }

    t = genr_get_output_type(u->gf);

    if (t == GRETL_TYPE_DOUBLE) {
        x = genr_get_output_scalar(u->gf);
    } else if (t == GRETL_TYPE_MATRIX) {
        gretl_matrix *m = genr_get_output_matrix(u->gf);

        if (gretl_matrix_is_scalar(m)) {
            x = m->val[0];
        } else {
            err = E_TYPES;
        }
    } else {
        err = E_TYPES;
    }

    u->fx_out = x;

    return x;
}

static double user_get_criterion2 (double b, void *p)
{
    umax *u = (umax *) p;
    double x = NADBL;
    int t, err;

    gretl_scalar_set_value(u->pxname, b);
    err = execute_genr(u->gf, u->dset, u->prn);

    if (err) {
        return NADBL;
    }

    t = genr_get_output_type(u->gf);

    if (t == GRETL_TYPE_DOUBLE) {
        x = genr_get_output_scalar(u->gf);
    } else if (t == GRETL_TYPE_MATRIX) {
        gretl_matrix *m = genr_get_output_matrix(u->gf);

        if (gretl_matrix_is_scalar(m)) {
            x = m->val[0];
        } else {
            err = E_TYPES;
        }
    } else {
        err = E_TYPES;
    }

    u->fx_out = x;

    return x;
}

/* user-defined optimizer: get the gradient, if specified */

static int user_get_gradient (double *b, double *g, int k,
                              BFGS_CRIT_FUNC func, void *p)
{
    umax *u = (umax *) p;
    gretl_matrix *ug;
    int i, err;

    for (i=0; i<k; i++) {
        u->b->val[i] = b[i];
    }

    err = execute_genr(u->gg, u->dset, u->prn);

    if (err) {
        return err;
    }

    ug = get_matrix_by_name(u->gmname);

    if (ug == NULL) {
        err = E_UNKVAR;
    } else if (gretl_vector_get_length(ug) != k) {
        err = E_NONCONF;
    } else {
        for (i=0; i<k; i++) {
            g[i] = ug->val[i];
        }
    }

    return err;
}

/* user-defined optimizer: get the hessian, if specified */

static int user_get_hessian (double *b, gretl_matrix *H,
                             void *p)
{
    umax *u = (umax *) p;
    gretl_matrix *uH;
    int k = H->rows;
    int i, err;

    for (i=0; i<k; i++) {
        u->b->val[i] = b[i];
    }

    err = execute_genr(u->gh, u->dset, u->prn);

    if (err) {
        return err;
    }

    uH = get_matrix_by_name(u->hmname);

    if (uH == NULL) {
        err = E_UNKVAR;
    } else if (uH->rows != k || uH->cols != k) {
        err = E_NONCONF;
    } else {
        gretl_matrix_copy_values(H, uH);
    }

    return err;
}

/* parse the name of the user gradient matrix (vector) or
   Hessian out of the associated function call, where it
   must be the first argument, given in pointer form
*/

int optimizer_get_matrix_name (const char *fncall, char *name)
{
    const char *s = strchr(fncall, '(');
    int n, err = 0;

    if (s == NULL) {
        err = E_DATA;
    } else {
        s++;
        s += strspn(s, " ");
        if (*s != '&') {
            err = E_TYPES;
        } else {
            s++;
            n = gretl_namechar_spn(s);
            if (n >= VNAMELEN) {
                err = E_DATA;
            } else {
                strncat(name, s, n);
            }
        }
    }

    return err;
}

static int check_optimizer_scalar_parm (umax *u, const char *s)
{
    int n = gretl_namechar_spn(s);
    int err = 0;

    if (n > 0 && n < VNAMELEN) {
        char vname[VNAMELEN];

        *vname = '\0';
        strncat(vname, s, n);
        if (!gretl_is_scalar(vname)) {
            err = gretl_scalar_add(vname, NADBL);
        }
        if (!err) {
            u->pxname = gretl_strdup(vname);
        }
    } else if (n >= VNAMELEN) {
        err = E_INVARG;
    }

    return err;
}

/* Ensure that we can find the specified callback function,
   and that it cannot replace/move its param-vector argument,
   given by u->b.
*/

static int check_optimizer_callback (umax *u, const char *fncall)
{
    int n = strcspn(fncall, "(");
    int err = 0;

    if (n > 0 && n < FN_NAMELEN) {
        char fname[FN_NAMELEN];
        user_var *uvar = NULL;
        ufunc *ufun;

        *fname = '\0';
        strncat(fname, fncall, n);
        ufun = get_user_function_by_name(fname);
        if (u->b == NULL) {
            check_optimizer_scalar_parm(u, fncall + n + 1);
        } else {
            uvar = get_user_var_by_data(u->b);
        }
        if (ufun != NULL && uvar != NULL) {
            user_var_set_flag(uvar, UV_NOREPL);
        }
    } else if (n >= FN_NAMELEN) {
        err = E_INVARG;
    }

    return err;
}

static int user_gen_setup (umax *u,
                           const char *fncall,
                           const char *gradcall,
                           const char *hesscall,
                           DATASET *dset)
{
    gchar *formula;
    int err;

    err = check_optimizer_callback(u, fncall);
    if (!err && gradcall != NULL) {
        err = check_optimizer_callback(u, gradcall);
    }
    if (!err && hesscall != NULL) {
        err = check_optimizer_callback(u, hesscall);
    }
    if (err) {
        return err;
    }

    formula = g_strdup_printf("$umax=%s", fncall);
    u->gf = genr_compile(formula, dset, u->gentype, OPT_P,
                         u->prn, &err);

    if (!err && gradcall != NULL) {
        /* process gradient formula */
        err = optimizer_get_matrix_name(gradcall, u->gmname);
        if (!err) {
            u->gg = genr_compile(gradcall, dset, GRETL_TYPE_ANY,
                                 OPT_P | OPT_O, u->prn, &err);
        }
    }

    if (!err && hesscall != NULL) {
        /* process Hessian formula */
        err = optimizer_get_matrix_name(hesscall, u->hmname);
        if (!err) {
            u->gh = genr_compile(hesscall, dset, GRETL_TYPE_ANY,
                                 OPT_P | OPT_O, u->prn, &err);
        }
    }

    if (!err) {
        u->dset = dset;
        u->fm_out = genr_get_output_matrix(u->gf);
    } else {
        umax_unset_no_replace_flag(u);
        destroy_genr(u->gf);
        destroy_genr(u->gg);
        u->gf = u->gg = NULL;
    }

    g_free(formula);

    return err;
}

double user_BFGS (gretl_matrix *b,
                  const char *fncall,
                  const char *gradcall,
                  DATASET *dset,
                  const gretl_matrix *bounds,
                  int minimize, PRN *prn,
                  int *err)
{
    umax *u;
    gretlopt opt = OPT_NONE;
    int maxit = BFGS_MAXITER_DEFAULT;
    int verbose, fcount = 0, gcount = 0;
    double tol;
    double ret = NADBL;

    u = umax_new(GRETL_TYPE_DOUBLE);
    if (u == NULL) {
        *err = E_ALLOC;
        return ret;
    }

    u->ncoeff = gretl_vector_get_length(b);
    if (u->ncoeff == 0) {
        *err = E_DATA;
        goto bailout;
    }

    u->b = b;
    u->prn = prn; /* placement of this? */

    *err = user_gen_setup(u, fncall, gradcall, NULL, dset);
    if (*err) {
        return NADBL;
    }

    tol = libset_get_double(BFGS_TOLER);
    verbose = libset_get_int(MAX_VERBOSE);

    if (verbose) {
        opt |= OPT_V;
    }

    if (minimize) {
        opt |= OPT_I;
    }

    if (bounds != NULL) {
        *err = BFGS_cmax(b->val, u->ncoeff,
                         maxit, tol, &fcount, &gcount,
                         user_get_criterion, C_OTHER,
                         (u->gg == NULL)? NULL : user_get_gradient,
                         u, bounds, opt, prn);
    } else {
        *err = BFGS_max(b->val, u->ncoeff,
                        maxit, tol, &fcount, &gcount,
                        user_get_criterion, C_OTHER,
                        (u->gg == NULL)? NULL : user_get_gradient,
                        u, NULL, opt, prn);
    }

    if (fcount > 0 && (verbose || !gretl_looping())) {
        pprintf(prn, _("Function evaluations: %d\n"), fcount);
        pprintf(prn, _("Evaluations of gradient: %d\n"), gcount);
    }

    if (!*err) {
        ret = u->fx_out;
    }

 bailout:

    umax_destroy(u);

    return ret;
}

double user_NR (gretl_matrix *b,
                const char *fncall,
                const char *gradcall,
                const char *hesscall,
                DATASET *dset,
                int minimize,
                PRN *prn, int *err)
{
    umax *u;
    double crittol = 1.0e-7;
    double gradtol = 1.0e-7;
    gretlopt opt = OPT_NONE;
    int maxit = 100;
    int iters = 0;
    double ret = NADBL;

    u = umax_new(GRETL_TYPE_DOUBLE);
    if (u == NULL) {
        *err = E_ALLOC;
        return ret;
    }

    u->ncoeff = gretl_vector_get_length(b);
    if (u->ncoeff == 0) {
        *err = E_DATA;
        goto bailout;
    }

    u->b = b;

    *err = user_gen_setup(u, fncall, gradcall, hesscall, dset);
    if (*err) {
        fprintf(stderr, "user_NR: failed on user_gen_setup\n");
        return NADBL;
    }

    if (libset_get_int(MAX_VERBOSE)) {
        opt |= OPT_V;
    }

    if (minimize) {
        opt |= OPT_I;
    }

    u->prn = prn; /* 2015-03-10: this was conditional on OPT_V */

    *err = newton_raphson_max(b->val, u->ncoeff, maxit,
                              crittol, gradtol,
                              &iters, C_OTHER,
                              user_get_criterion,
                              (u->gg == NULL)? NULL : user_get_gradient,
                              (u->gh == NULL)? NULL : user_get_hessian,
                              u, opt, prn);

    if (!*err) {
        ret = u->fx_out;
    }

 bailout:

    umax_destroy(u);

    return ret;
}

double deriv_free_optimize (MaxMethod method,
                            gretl_matrix *b,
                            const char *fncall,
                            int maxit,
                            double tol,
                            int minimize,
                            DATASET *dset,
                            PRN *prn, int *err)
{
    umax *u;
    gretlopt opt = OPT_NONE;
    double ret = NADBL;

    if (b->is_complex) {
        *err = E_CMPLX;
        return ret;
    }

    u = umax_new(GRETL_TYPE_DOUBLE);
    if (u == NULL) {
        *err = E_ALLOC;
        return ret;
    }

    u->ncoeff = gretl_vector_get_length(b);
    if (u->ncoeff == 0 || (method == ROOT_FIND && u->ncoeff != 2)) {
        *err = E_INVARG;
        goto bailout;
    }

    if (method != ROOT_FIND) {
        u->b = b;
    }

    *err = user_gen_setup(u, fncall, NULL, NULL, dset);
    if (*err) {
        return NADBL;
    }

    if (libset_get_int(MAX_VERBOSE)) {
        opt = OPT_V;
        u->prn = prn;
    }

    if (minimize) {
        opt |= OPT_I;
    }

    if (method == SIMANN_MAX) {
        *err = gretl_simann(b->val, u->ncoeff, maxit,
                            user_get_criterion, u,
                            opt, prn);
    } else if (method == NM_MAX) {
        *err = gretl_amoeba(b->val, u->ncoeff, maxit,
                            user_get_criterion, u,
                            opt, prn);
    } else if (method == GSS_MAX) {
        int iters = 0;

        *err = gretl_gss(b->val, tol, &iters,
                         user_get_criterion, u,
                         opt, prn);
    } else if (method == ROOT_FIND) {
        *err = gretl_fzero(b->val, tol,
                           user_get_criterion2, u,
                           &ret, opt, prn);
    }

    if (!*err && method != ROOT_FIND) {
        ret = user_get_criterion(b->val, u);
    }

 bailout:

    umax_destroy(u);

    return ret;
}

#define JAC_DEBUG 0

static int user_calc_fvec (int m, int n, double *x, double *fvec,
                           int *iflag, void *p)
{
    umax *u = (umax *) p;
    gretl_matrix *v;
    int i, err;

    for (i=0; i<n; i++) {
        u->b->val[i] = x[i];
    }

#if JAC_DEBUG
    gretl_matrix_print(u->b, "user_calc_fvec: u->b");
#endif

    err = execute_genr(u->gf, u->dset, u->prn);
    if (err) {
        fprintf(stderr, "execute_genr: err = %d\n", err);
    }

    if (err) {
        *iflag = -1;
        return 0;
    }

    v = genr_get_output_matrix(u->gf);

#if JAC_DEBUG
    gretl_matrix_print(v, "matrix from u->f");
#endif

    if (v == NULL || gretl_vector_get_length(v) != m) {
        fprintf(stderr, "user_calc_fvec: got bad matrix\n");
        *iflag = -1;
    } else {
        for (i=0; i<m; i++) {
            fvec[i] = v->val[i];
        }
    }

    return 0;
}

static int fdjac_allocate (int m, int n,
                           gretl_matrix **J,
                           double **w, double **f)
{
    *J = gretl_matrix_alloc(m, n);
    if (*J == NULL) {
        return E_ALLOC;
    }

    *w = malloc(m * sizeof **w);
    *f = malloc(m * sizeof **f);

    if (*w == NULL || *f == NULL) {
        return E_ALLOC;
    }

    return 0;
}

gretl_matrix *user_fdjac (gretl_matrix *theta, const char *fncall,
                          double eps, DATASET *dset, int *err)
{
    umax *u;
    gretl_matrix *J = NULL;
    int m, n, nnf = 0;
    int iflag = 0;
    double *wa = NULL;
    double *fvec = NULL;
    int i;

    *err = 0;

    u = umax_new(GRETL_TYPE_MATRIX);
    if (u == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    n = gretl_vector_get_length(theta);
    if (n == 0) {
        fprintf(stderr, "fdjac: gretl_vector_get_length gave %d for theta\n",
                n);
        *err = E_DATA;
        return NULL;
    }

    u->b = theta;
    u->ncoeff = n;

    *err = user_gen_setup(u, fncall, NULL, NULL, dset);
    if (*err) {
        fprintf(stderr, "fdjac: error %d from user_gen_setup\n", *err);
        goto bailout;
    }

    if (u->fm_out == NULL) {
        fprintf(stderr, "fdjac: u.fm_out is NULL\n");
        *err = E_DATA; /* FIXME */
        goto bailout;
    }

    m = gretl_vector_get_length(u->fm_out);
    if (m == 0) {
        *err = E_DATA;
        goto bailout;
    }

    *err = fdjac_allocate(m, n, &J, &wa, &fvec);
    if (*err) {
        goto bailout;
    }

    for (i=0; i<m; i++) {
        fvec[i] = u->fm_out->val[i];
        if (na(fvec[i])) {
            nnf++;
        }
    }

    if (nnf > 0) {
        gretl_errmsg_sprintf(_("fdjac: got %d non-finite value(s) on input"),
                             nnf);
        *err = E_DATA;
    } else {
        int quality = libset_get_int(FDJAC_QUAL);

        fdjac2_(user_calc_fvec, m, n, quality, theta->val, fvec, J->val,
                m, &iflag, eps, wa, u);
    }

 bailout:

    free(wa);
    free(fvec);

    if (*err) {
        gretl_matrix_free(J);
        J = NULL;
    }

    umax_destroy(u);

    return J;
}

gretl_matrix *user_numhess (gretl_matrix *b, const char *fncall,
                            double d, DATASET *dset, int *err)
{
    gretl_matrix *H = NULL;
    struct uhess_data uh = {0};
    gchar *formula = NULL;
    int n;

    if (get_user_var_by_data(b) == NULL) {
        fprintf(stderr, "numhess: b must be a named matrix\n");
        *err = E_DATA;
        return NULL;
    }

    n = gretl_vector_get_length(b);
    if (n == 0) {
        fprintf(stderr, "numhess: gretl_vector_get_length gave %d for b\n", n);
        *err = E_DATA;
        return NULL;
    }

    if (!*err) {
        formula = g_strdup_printf("$umax=%s", fncall);
        uh.genr = genr_compile(formula, dset, GRETL_TYPE_DOUBLE, OPT_P,
                               NULL, err);
    }

    if (*err) {
        fprintf(stderr, "numhess: error %d from genr_compile\n", *err);
    } else {
        H = gretl_zero_matrix_new(n, n);
        if (H == NULL) {
            *err = E_ALLOC;
        }
    }

    if (!*err) {
        uh.dset = dset;
        *err = numerical_hessian(b->val, H, uhess_callback,
                                 &uh, 0, d);
    }

    g_free(formula);
    user_var_delete_by_name("$umax", NULL);
    destroy_genr(uh.genr);

    if (*err && H != NULL) {
        gretl_matrix_free(H);
        H = NULL;
    }

    return H;
}

/* Below: Newton-Raphson code, starting with a few
   auxiliary functions */

static int broken_matrix (const gretl_matrix *m)
{
    int i, n = m->rows * m->cols;

    for (i=0; i<n; i++) {
        if (na(m->val[i])) {
            return 1;
        }
    }

    return 0;
}

static void copy_to (double *targ, const double *src, int n)
{
    int i;

    for (i=0; i<n; i++) {
        targ[i] = src[i];
    }
}

static void copy_plus (double *targ, const double *src,
                       double step, const double *a, int n)
{
    int i;

    for (i=0; i<n; i++) {
        targ[i] = src[i] + step * a[i];
    }
}

static double scalar_xpx (const gretl_matrix *m)
{
    double xpx = 0.0;
    int i;

    for (i=0; i<m->rows; i++) {
        xpx += m->val[i] * m->val[i];
    }

    return xpx;
}

enum {
    GRADTOL_MET = 1,
    CRITTOL_MET,
    STEPMIN_MET
};

static void print_NR_status (int status, double crittol, double gradtol,
                             double sumgrad, PRN *prn)
{
    int msgs = gretl_messages_on();

    if (status == GRADTOL_MET) {
        if (msgs) {
            pprintf(prn, _("Gradient within tolerance (%g)\n"), gradtol);
        }
    } else if (status == CRITTOL_MET) {
        if (msgs) {
            pprintf(prn, _("Successive criterion values within tolerance (%g)\n"),
                    crittol);
        }
    } else if (status == STEPMIN_MET && gretl_warnings_on()) {
        if (sumgrad > 0) {
            pprintf(prn, _("Warning: couldn't improve criterion (gradient = %g)\n"),
                    sumgrad);
        } else {
            pprintf(prn, _("Warning: couldn't improve criterion\n"));
        }
    }
}

/*
   The strategy here may be simplistic, but appears to be effective in
   quite a few cases: If the (negative) Hessian is not pd, we try to
   nudge it towards positive definiteness without losing too much
   information on curvature by downsizing the off-diagonal elements of
   H. In practice, we use

   C = \lambda H + (1-lambda) diag(H)

   where lambda is reasonably close to 1. Note that, if lambda was 0,
   a NR iteration would just amount to an iteration of the simple
   gradient method (on a scaled version of the coefficients).
   Hopefully, a few of those should be able to "tow us away" from the
   non-pd region. In desperate cases, we use the a diagonal matrix
   with the absolute values of H_{i,i} plus one.
*/

static int NR_invert_hessian (gretl_matrix *H, const gretl_matrix *Hcpy)
{
    double hii, hmin = 1.0e-28; /* was 1.0e-20 */
    int i, j, n = H->rows;
    int restore = 0;
    double x;
    int err = 0;

    /* first, check if all the elements along the diagonal are
       numerically positive
    */
    for (i=0; i<n && !err; i++) {
        hii = gretl_matrix_get(H, i, i);
        if (hii < hmin) {
            err = 1;
            break;
        }
    }

    if (err) {
        fprintf(stderr, "NR_invert_hessian: non-positive "
                "diagonal (H(%d,%d) = %g)\n", i+1, i+1, hii);
    } else {
        err = gretl_invert_symmetric_matrix(H);

        if (err == E_NOTPD) {
            double lambda = 1.0;
            int s;

            for (s=0; lambda > 0.1 && err; s++) {
                lambda *= 0.8;
                fprintf(stderr, "newton hessian fixup: round %d, lambda=%g\n",
                        s, lambda);
                /* restore the original H */
                gretl_matrix_copy_values(H, Hcpy);
                for (i=0; i<n; i++) {
                    for (j=0; j<i; j++) {
                        x = lambda * gretl_matrix_get(H, i, j);
                        gretl_matrix_set(H, i, j, x);
                        gretl_matrix_set(H, j, i, x);
                    }
                }
                err = gretl_invert_symmetric_matrix(H);
            }
            restore = (err != 0);
        }
    }

    if (err) {
        fprintf(stderr, "NR_invert_hessian: major surgery needed\n");
        if (restore) {
            gretl_matrix_copy_values(H, Hcpy);
        }
        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                x = (i==j) ? 1/(1 + fabs(gretl_matrix_get(H, i, j))): 0;
                gretl_matrix_set(H, i, j, x);
            }
        }
    }

    return 0;
}

/**
 * newton_raphson_max:
 * @b: array of adjustable coefficients.
 * @n: number elements in array @b.
 * @maxit: the maximum number of iterations to allow.
 * @crittol: tolerance for terminating iteration, in terms of
 * the change in the criterion.
 * @gradtol: tolerance for terminating iteration, in terms of
 * the gradient.
 * @itercount: location to receive count of iterations.
 * @crittype: code for type of the maximand/minimand: should
 * be %C_LOGLIK, %C_GMM or %C_OTHER.  Used only in printing
 * iteration info.
 * @cfunc: pointer to function used to calculate maximand.
 * @gradfunc: pointer to function used to calculate the
 * gradient, or %NULL for default numerical calculation.
 * @hessfunc: pointer to function used to calculate the
 * Hessian.
 * @data: pointer that will be passed as the last
 * parameter to the callback functions @cfunc, @gradfunc
 * and @hessfunc.
 * @opt: may contain %OPT_V for verbose operation.
 * @prn: printing struct (or %NULL).  Only used if @opt
 * includes %OPT_V.
 *
 * The functions @cfunc (computes the criterion, usually a
 * loglikelihood) and @gradfunc (score, provides an updated
 * estimate of the gradient in its second argument) are just as in
 * BFGS_max above.
 *
 * The @hessfunc callback should compute the negative Hessian,
 * _not_ inverted; newton_raphson_max takes care of the inversion,
 * with a routine for fixing up the matrix if it's not positive
 * definite. If @hessfunc is NULL we fall back on a numerical
 * approximation to the Hessian.
 *
 * Returns: 0 on successful completion, non-zero error code
 * on error.
 */

int newton_raphson_max (double *b, int n, int maxit,
                        double crittol, double gradtol,
                        int *itercount, int crittype,
                        BFGS_CRIT_FUNC cfunc,
                        BFGS_GRAD_FUNC gradfunc,
                        HESS_FUNC hessfunc,
                        void *data, gretlopt opt,
                        PRN *prn)
{
    BFGS_GRAD_FUNC realgrad = NULL;
    int verbose = (opt & OPT_V);
    int minimize = (opt & OPT_I);
    double stepmin = 1.0e-6;
    gretl_matrix_block *B;
    gretl_matrix *H0, *H1;
    gretl_matrix *g, *a;
    double *b0, *b1;
    double f0, f1, sumgrad = 0;
    double steplen = 1.0;
    int status = 0;
    int iter = 0;
    int err = 0;

    b0 = malloc(2 * n * sizeof *b0);
    if (b0 == NULL) {
        return E_ALLOC;
    }

    B = gretl_matrix_block_new(&H0, n, n,
                               &H1, n, n,
                               &g, n, 1,
                               &a, n, 1,
                               NULL);
    if (B == NULL) {
        free(b0);
        return E_ALLOC;
    }

    if (gradfunc == NULL) {
        gradfunc = numeric_gradient;
    } else {
        realgrad = gradfunc;
    }

    /* needs some work */
    optim_get_user_values(b, n, NULL, NULL, NULL, opt, prn);

    b1 = b0 + n;
    copy_to(b1, b, n);

    f1 = optim_fncall(cfunc, b1, data, minimize);
    if (na(f1)) {
        gretl_errmsg_set(_("Initial value of objective function is not finite"));
        err = E_NAN;
    }

    if (!err) {
        err = optim_gradcall(gradfunc, b, g->val, n, cfunc, data, minimize);
        if (err) {
            fprintf(stderr, "newton_raphson_max: err = %d from %s\n", err,
                    realgrad ? "supplied gradfunc" : "numeric_gradient");
        }
    }

    if (!err) {
        if (hessfunc != NULL) {
            err = hessfunc(b, H1, data);
        } else {
            err = NR_fallback_hessian(b, H1, realgrad, cfunc, data);
        }
        if (!err) {
            if (minimize) {
                gretl_matrix_multiply_by_scalar(H1, -1.0);
            }
            gretl_matrix_copy_values(H0, H1);
            err = NR_invert_hessian(H1, H0);
        }
    }

    gretl_iteration_push();

    while (status == 0 && !err) {
        iter++;
        steplen = 1.0;
        f0 = f1;

        copy_to(b0, b1, n);

        if (broken_matrix(g)) {
            fprintf(stderr, "NA in gradient\n");
            err = E_NAN;
            break;
        }

        gretl_matrix_copy_values(H0, H1);
        if (broken_matrix(H0)) {
            fprintf(stderr, "NA in Hessian\n");
            err = E_NAN;
            break;
        }

        /* apply quadratic approximation */
        gretl_matrix_multiply(H0, g, a);
        copy_plus(b1, b0, steplen, a->val, n);
        f1 = optim_fncall(cfunc, b1, data, minimize);

        while (na(f1) || ((f1 < f0) && (steplen > stepmin))) {
            /* try smaller step */
            steplen /= 2.0;
            copy_plus(b1, b0, steplen, a->val, n);
            f1 = optim_fncall(cfunc, b1, data, minimize);
        }

        if (verbose) {
            print_iter_info(iter, f1, crittype, n, b1, g->val,
                            steplen, prn);
        }

        err = optim_gradcall(gradfunc, b1, g->val, n, cfunc, data,
                             minimize);

        if (err || broken_matrix(g)) {
            err = (err == 0)? E_NAN : err;
            break;
        }

        if (hessfunc != NULL) {
            err = hessfunc(b1, H1, data);
        } else {
            err = NR_fallback_hessian(b1, H1, realgrad, cfunc, data);
        }

        if (err || broken_matrix(H1)) {
            err = (err == 0)? E_NAN : err;
            break;
        }

        if (!err) {
            if (minimize) {
                gretl_matrix_multiply_by_scalar(H1, -1.0);
            }
            gretl_matrix_copy_values(H0, H1);
            err = NR_invert_hessian(H1, H0);
            if (err) {
                break;
            }
        }

        sumgrad = sqrt(scalar_xpx(g));

        if (steplen < stepmin) {
            status = STEPMIN_MET;
        } else if (iter > maxit) {
            err = E_NOCONV;
        } else if (sumgrad < gradtol) {
            status = GRADTOL_MET;
        } else if (f1 - f0 < crittol) {
            status = CRITTOL_MET;
        }
    }

    gretl_iteration_pop();

    if (verbose) {
        print_iter_info(-1, f1, crittype, n, b1, g->val,
                        steplen, prn);
        pputc(prn, '\n');
    }

    *itercount = iter;

    if (!err) {
        copy_to(b, b1, n);
        if (prn != NULL) {
            print_NR_status(status, crittol, gradtol, sumgrad, prn);
        }
    }

    free(b0);
    gretl_matrix_block_destroy(B);

    return err;
}

static double simann_call (BFGS_CRIT_FUNC cfunc,
                           double *b, void *data,
                           int minimize)
{
    double ret = cfunc(b, data);

    return na(ret) ? NADBL : minimize ? -ret : ret;
}

/**
 * gretl_simann:
 * @theta: parameter array.
 * @n: length of @theta.
 * @maxit: the maximum number of iterations to perform.
 * @cfunc: the function to be maximized.
 * @data: pointer to be passed to the @cfunc callback.
 * @opt: may include %OPT_V for verbose operation.
 * @prn: printing struct, or NULL.
 *
 * Simulated annealing: can help to improve the initialization
 * of @theta for numerical optimization. On exit the value of
 * @theta is set to the func-best point in case of improvement,
 * otherwise to the last point visited.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_simann (double *theta, int n, int maxit,
                  BFGS_CRIT_FUNC cfunc, void *data,
                  gretlopt opt, PRN *prn)
{
    gretl_matrix b;
    gretl_matrix *b0 = NULL;
    gretl_matrix *b1 = NULL;
    gretl_matrix *bstar = NULL;
    gretl_matrix *d = NULL;
    double f0, f1;
    double fbest, fworst;
    double Temp = 1.0;
    double radius = 1.0;
    int improved = 0;
    int minimize;
    int i, err = 0;

    if (maxit <= 0) {
        maxit = 1024;
    }

    /* maximize by default, but minimize if OPT_I is given */
    minimize = (opt & OPT_I)? 1 : 0;

    gretl_matrix_init_full(&b, n, 1, theta);

    b0 = gretl_matrix_copy(&b);
    b1 = gretl_matrix_copy(&b);
    bstar = gretl_matrix_copy(&b);
    d = gretl_column_vector_alloc(n);

    if (b0 == NULL || b1 == NULL || bstar == NULL || d == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    f0 = fbest = fworst = simann_call(cfunc, b.val, data, minimize);

    if (opt & OPT_V) {
        pprintf(prn, _("\nSimulated annealing: initial function value = %.8g\n"),
                f0);
    }

    gretl_iteration_push();

    /* Question: should the initial radius be a function of the scale
       of the initial parameter vector?
    */

    for (i=0; i<maxit; i++) {
        gretl_matrix_random_fill(d, D_NORMAL);
        gretl_matrix_multiply_by_scalar(d, radius);
        gretl_matrix_add_to(b1, d);
        f1 = simann_call(cfunc, b1->val, data, minimize);

        if (!na(f1) && (f1 > f0 || gretl_rand_01() < Temp)) {
            /* jump to the new point */
            f0 = f1;
            gretl_matrix_copy_values(b0, b1);
            if (f0 > fbest) {
                fbest = f0;
                gretl_matrix_copy_values(bstar, b0);
                if (opt & OPT_V) {
                    if (!improved) {
                        pprintf(prn, "\n%6s %12s %12s %12s\n",
                                "iter", "temp", "radius", "fbest");
                    }
                    pprintf(prn, "%6d %#12.6g %#12.6g %#12.6g\n",
                            i, Temp, radius, fbest);
                }
                improved = 1;
            } else if (f0 < fworst) {
                fworst = f0;
            }
        } else {
            /* revert to where we were */
            gretl_matrix_copy_values(b1, b0);
            f1 = f0;
        }

        Temp *= 0.999;
        radius *= 0.9999;
    }

    gretl_iteration_pop();

    if (improved) {
        /* use the best point */
        gretl_matrix_copy_values(&b, bstar);
    } else {
        /* use the last point */
        gretl_matrix_copy_values(&b, b0);
    }

    if (improved) {
        if (opt & OPT_V) pputc(prn, '\n');
    } else {
        pprintf(prn, _("No improvement found in %d iterations\n\n"), maxit);
    }

    if (fbest - fworst < 1.0e-9) {
        pprintf(prn, _("*** warning: surface seems to be flat\n"));
    }

 bailout:

    gretl_matrix_free(b0);
    gretl_matrix_free(b1);
    gretl_matrix_free(bstar);
    gretl_matrix_free(d);

    return err;
}

/*
  nelder_mead: this is based closely on the nelmin function as written
  in C by John Burkardt; see
  http://people.sc.fsu.edu/~jburkardt/c_src/asa047
  It is converted to a maximizer, and modified for use with the
  gretl library.

  Burkardt's code is distributed under the GNU LGPL license, and was
  last modified on 28 October 2010. It draws on original Fortran code
  by R. O'Neill, for which see "Algorithm AS 47: Function Minimization
  Using a Simplex Procedure", Journal of the Royal Statistical
  Society, Series C (Applied Statistics), Vol. 20, No. 3, 1971.
*/

static double nm_call (BFGS_CRIT_FUNC cfunc,
                       double *b, void *data,
                       int *ncalls, int minimize)
{
    double ret = cfunc(b, data);

    *ncalls += 1;

    return minimize ? ret : na(ret) ? ret : -ret;
}

static int
nelder_mead (BFGS_CRIT_FUNC cfunc, int n, double start[],
             double xmin[], double *ynewlo, double reqmin,
             double step[], int maxcalls, int *ncalls, int *nresets,
             void *data, gretlopt opt, PRN *prn)
{
    gretl_matrix *pmat = NULL;
    double ccoeff = 0.5;
    double ecoeff = 2.0;
    double rcoeff = 1.0;
    double eps = 0.001;
    double del = 1.0;
    double *wspace;
    double *p;
    double *p2star;
    double *pbar;
    double *pstar;
    double *y;
    double rq, x, z;
    double ylo, ystar, y2star;
    /* frequency of convergence check */
    int konvge = 10;
    int i, ihi, ilo;
    int l, j, jcount;
    int outer, inner;
    int getmin;
    int err = 0;

    /* use a gretl_matrix in case we want to print it */
    pmat = gretl_matrix_alloc(n, n + 1);
    wspace = malloc((4*n + 1) * sizeof *wspace);

    if (pmat == NULL || wspace == NULL) {
        gretl_matrix_free(pmat);
        free(wspace);
        return E_ALLOC;
    }

    p = pmat->val;
    pstar = wspace;
    p2star = pstar + n;
    pbar = p2star + n;
    y = pbar + n;

    *ncalls = *nresets = 0;
    jcount = konvge;
    rq = reqmin * n;

    /* maximize by default, but minimize if OPT_I is given */
    getmin = (opt & OPT_I)? 1 : 0;

    gretl_iteration_push();

    for (outer=1; ; outer++) {
        for (i=0; i<n; i++) {
            p[i+n*n] = start[i];
        }
        y[n] = nm_call(cfunc, start, data, ncalls, getmin);

        if (opt & OPT_V) {
            if (outer == 1) {
                pprintf(prn, _("\nOuter iteration %d: function value = %#g\n"),
                        outer, y[n]);
            } else {
                pprintf(prn, _("Outer iteration %d (reset)\n"), outer);
            }
        }

        /* construct the simplex */
        for (j=0; j<n; j++) {
            x = start[j];
            start[j] += step[j] * del;
            for (i=0; i<n; i++) {
                p[i+j*n] = start[i];
            }
            y[j] = nm_call(cfunc, start, data, ncalls, getmin);
            start[j] = x;
        }

        /* find the lowest y value */
        ylo = y[0];
        ilo = 0;
        for (i=1; i<=n; i++) {
            if (y[i] < ylo) {
                ylo = y[i];
                ilo = i;
            }
        }

        for (inner=1; *ncalls < maxcalls; inner++) {
            *ynewlo = y[0];
            ihi = 0;

            for (i=1; i<=n; i++) {
                if (*ynewlo < y[i]) {
                    *ynewlo = y[i];
                    ihi = i;
                }
            }
            /* calculate pbar, the centroid of the simplex */
            for (i=0; i<n; i++) {
                z = 0.0;
                for (j=0; j<=n; j++) {
                    z += p[i+j*n];
                }
                z -= p[i+ihi*n];
                pbar[i] = z / n;
            }

            /* reflection through the centroid */
            for (i=0; i<n; i++) {
                pstar[i] = pbar[i] + rcoeff * (pbar[i] - p[i+ihi*n]);
            }
            ystar = nm_call(cfunc, pstar, data, ncalls, getmin);

            if ((opt & OPT_V) && (inner == 1 || inner % 10 == 0)) {
                pprintf(prn, _(" inner iter %3d: function value %#g\n"),
                        inner, ystar);
            }

            if (ystar < ylo) {
                /* successful reflection, so extension */
                for (i=0; i<n; i++) {
                    p2star[i] = pbar[i] + ecoeff * (pstar[i] - pbar[i]);
                }
                y2star = nm_call(cfunc, p2star, data, ncalls, getmin);
                /* check extension */
                if (ystar < y2star) {
                    for (i=0; i<n; i++) {
                        p[i+ihi*n] = pstar[i];
                    }
                    y[ihi] = ystar;
                } else {
                    for (i=0; i<n; i++) {
                        p[i+ihi*n] = p2star[i];
                    }
                    y[ihi] = y2star;
                }
            } else {
                l = 0;
                for (i=0; i<=n; i++) {
                    if (ystar < y[i]) {
                        l++;
                    }
                }
                if (l > 1) {
                    for (i=0; i<n; i++) {
                        p[i+ihi*n] = pstar[i];
                    }
                    y[ihi] = ystar;
                } else if (l == 0) {
                    for (i=0; i<n; i++) {
                        p2star[i] = pbar[i] + ccoeff * (p[i+ihi*n] - pbar[i]);
                    }
                    y2star = nm_call(cfunc, p2star, data, ncalls, getmin);
                    /* contract the whole simplex */
                    if (y[ihi] < y2star) {
                        for (j=0; j<=n; j++) {
                            for (i=0; i<n; i++) {
                                p[i+j*n] = (p[i+j*n] + p[i+ilo*n]) * 0.5;
                                xmin[i] = p[i+j*n];
                            }
                            y[j] = nm_call(cfunc, xmin, data, ncalls, getmin);
                        }
                        ylo = y[0];
                        ilo = 0;
                        for (i=1; i<=n; i++) {
                            if (y[i] < ylo) {
                                ylo = y[i];
                                ilo = i;
                            }
                        }
                        continue;
                    } else {
                        for (i=0; i<n; i++) {
                            p[i+ihi*n] = p2star[i];
                        }
                        y[ihi] = y2star;
                    }
                } else if (l == 1) {
                    for (i=0; i<n; i++) {
                        p2star[i] = pbar[i] + ccoeff * (pstar[i] - pbar[i]);
                    }
                    y2star = nm_call(cfunc, p2star, data, ncalls, getmin);
                    /* retain reflection? */
                    if (y2star <= ystar) {
                        for (i=0; i<n; i++) {
                            p[i+ihi*n] = p2star[i];
                        }
                        y[ihi] = y2star;
                    } else {
                        for (i=0; i<n; i++) {
                            p[i+ihi*n] = pstar[i];
                        }
                        y[ihi] = ystar;
                    }
                }
            }

            /* check if ylo improved */
            if (y[ihi] < ylo) {
                ylo = y[ihi];
                ilo = ihi;
            }
            jcount--;

            if (jcount > 0) {
                continue;
            }

            /* check to see if optimum reached */
            if (*ncalls <= maxcalls) {
                jcount = konvge;
                z = 0.0;
                for (i=0; i<=n; i++) {
                    z += y[i];
                }
                x = z / (n+1);
                z = 0.0;
                for (i=0; i<=n; i++) {
                    z += (y[i] - x) * (y[i] - x);
                }
                if (z <= rq) {
                    break;
                }
            }
        }

        /* check that ynewlo is a local minimum */
        for (i=0; i<n; i++) {
            xmin[i] = p[i+ilo*n];
        }
        *ynewlo = y[ilo];

        if (*ncalls > maxcalls) {
            err = E_NOCONV;
            break;
        }

        err = 0;

        for (i=0; i<n; i++) {
            double xsave = xmin[i];
            double dx = step[i] * eps;

            xmin[i] = xsave + dx;
            z = nm_call(cfunc, xmin, data, ncalls, getmin);
            if (z < *ynewlo) {
                err = E_NOCONV;
                break;
            }
            xmin[i] = xsave - dx;
            z = nm_call(cfunc, xmin, data, ncalls, getmin);
            if (z < *ynewlo) {
                err = E_NOCONV;
                break;
            }
            xmin[i] = xsave;
        }

        if (err == 0) {
            if (opt & OPT_V) {
                double crit = getmin ? *ynewlo : -(*ynewlo);

                pprintf(prn, _("Found optimum %#g after %d function calls, "
                        "%d resets\n\n"), crit, *ncalls, *nresets);
            }
            break;
        }

        /* prepare to restart */
        for (i=0; i<n; i++) {
            start[i] = xmin[i];
        }
        del = eps;
        *nresets += 1;
    }

    gretl_iteration_pop();

    gretl_matrix_free(pmat);
    free(wspace);

    return err;
}

int gretl_amoeba (double *theta, int n, int maxit,
                  BFGS_CRIT_FUNC cfunc, void *data,
                  gretlopt opt, PRN *prn)
{
    int ncalls = 0;
    int nresets = 0;
    int maxcalls;
    double reqmin = 1.0e-12;
    double fval = 0.0;
    double *step;
    double *xmin;
    int i, err = 0;

    step = malloc(n * sizeof *step);
    xmin = malloc(n * sizeof *xmin);

    if (step == NULL || xmin == NULL) {
        free(step); free(xmin);
        return E_ALLOC;
    }

    for (i=0; i<n; i++) {
        step[i] = 1.0;
        xmin[i] = 0.0;
    }

    if (maxit < 0) {
        maxcalls = -maxit;
        opt |= OPT_F; /* fail on non-convergence */
    } else {
        maxcalls = maxit == 0 ? 2000 : maxit;
    }

    err = nelder_mead(cfunc, n, theta, xmin, &fval, reqmin,
                      step, maxcalls, &ncalls, &nresets,
                      data, opt, prn);

    fprintf(stderr, "asa047: fncalls=%d, err=%d\n", ncalls, err);

    if (err == E_NOCONV && !(opt & OPT_F)) {
        /* tolerate non-convergence */
        err = 0;
    }

    if (!err) {
        for (i=0; i<n; i++) {
            theta[i] = xmin[i];
        }
    }

    free(step);
    free(xmin);

    return err;
}

int gretl_gss (double *theta, double tol, int *ic,
               BFGS_CRIT_FUNC cfunc, void *data,
               gretlopt opt, PRN *prn)
{
    double gr = (sqrt(5.0) + 1) / 2.0;
    double a = theta[1];
    double b = theta[2];
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;
    double fc, fd;
    int minimize = (opt & OPT_I);
    int cond, iter = 1;
    int err = 0;

    if (na(tol)) {
        tol = 1.0e-4; /* ?? */
    }

    while (fabs(c - d) > tol) {
        theta[0] = c;
        fc = cfunc(theta, data);
        theta[0] = d;
        fd = cfunc(theta, data);
        if (opt & OPT_V) {
            pprintf(prn, "%d: bracket={%g,%g}, values={%g,%g}\n",
                    iter, c, d, fc, fd);
        }
        if (na(fc) || na(fd)) {
            err = E_NAN;
            break;
        }
        cond = minimize ? fc < fd : fc > fd;
        if (cond) {
            b = d;
        } else {
            a = c;
        }
        /* recompute c and d to preserve precision */
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
        iter++;
    }

    if (ic != NULL) {
        *ic = iter;
    }

    if (!err) {
        /* record the optimum and the bracket */
        theta[0] = (b + a) / 2;
        theta[1] = a;
        theta[2] = b;
    }

    return err;
}

static int sgn (double x)
{
    return x == 0 ? 0 : x < 0 ? -1 : 1;
}

/* Try to find a value for x1 which brackets a zero-crossing.
   Based on Octave's fzero.m but with a little extra stuff.
*/

static double find_x1 (double x0, double y0, double *py1,
                       ZFUNC zfunc, void *data)
{
    double srch[] = {
        -0.01, 0.025, -0.05, 0.10, -0.25, 0.50,
        -1.01, 2.5, -5, 10, -50, 100, -500, 1000
    };
    double x1, y1, a = x0;
    int warnsave, negbad = 0;
    int i;

    warnsave = libset_get_bool(WARNINGS);
    libset_set_bool(WARNINGS, 0);

    if (fabs(x0) < 0.001) {
        a = (x0 == 0)? 0.1 : sgn(x0) * 0.1;
    }

    errno = 0;
    *py1 = NADBL;

    for (i=0; i<G_N_ELEMENTS(srch); i++) {
        x1 = a + a * srch[i];
        if (negbad && x1 < 0) {
            x1 = -x1;
        }
        y1 = zfunc(x1, data);
        if (na(y1)) {
            if (errno == EDOM && x1 < 0) {
                negbad = 1;
            }
            errno = 0;
        } else if (sgn(y1) * sgn(y0) <= 0) {
            /* found a candidate */
            *py1 = y1;
            break;
        }
    }

    libset_set_bool(WARNINGS, warnsave);

    return na(*py1) ? NADBL : x1;
}

int gretl_fzero (double *bracket, double tol,
                 ZFUNC zfunc, void *data, double *px,
                 gretlopt opt, PRN *prn)
{
    int MAXITER = 80;
    double ytol = 1.0e-14;
    double xtol = 1.0e-15;
    double y, y0, y1, y2;
    double x, x0, x1, x2;
    double dx, d0, d1;
    int i, err = 0;

    x0 = bracket[0];
    x1 = bracket[1];

    if (!na(tol)) {
        ytol = tol;
    }
    if (na(x0)) {
        /* exact zero is too risky as default */
        x0 = 0.107;
    }

    y0 = zfunc(x0, data);
    if (na(y0)) {
        return E_NAN;
    }
    if (!na(x1)) {
        y1 = zfunc(x1, data);
        if (na(y1)) {
            return E_NAN;
        }
    } else {
        x1 = find_x1(x0, y0, &y1, zfunc, data);
        if (na(x1)) {
            return E_NOCONV;
        }
    }

    if (y0 == 0) {
        *px = x0;
        return 0;
    } else if (y1 == 0) {
        *px = x1;
        return 0;
    } else if (y0 * y1 > 0) {
        return E_NOCONV;
    }

    if (opt & OPT_V) {
        pprintf(prn, "Initial bracket {%#.6g, %#.6g}\n", x0, x1);
    }

    for (i=0; i<MAXITER && !err; i++) {
        x2 = (x0 + x1) / 2;
        y2 = zfunc(x2, data);
        /* exponential interpolation of x */
        x = x2 + (x2-x0) * sgn(y0-y1)*y2 / sqrt(y2*y2 - y0*y1);
        d0 = fabs(x-x0);
        d1 = fabs(x-x1);
        dx = d0 < d1 ? d0 : d1;
        y = zfunc(x, data);
        if (opt & OPT_V) {
            pprintf(prn, "Iter %2d: f(%#.9g) = %g\n", i+1, x, y);
        }
        if (dx < xtol || fabs(y) < ytol) {
            break;
        }
        if (sgn(y2) != sgn(y)) {
            x0 = x2;
            y0 = y2;
            x1 = x;
            y1 = y;
        } else if (sgn(y1) != sgn(y)) {
            x0 = x;
            y0 = y;
        } else {
            x1 = x;
            y1 = y;
        }
    }

    if (!err && fabs(y) > ytol) {
        x = NADBL;
        err = E_NOCONV;
    }

    if (!err) {
        /* record the root */
        *px = x;
    }

    return err;
}
