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

#include "../../minpack/minpack.h"
#include <float.h> 

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
    double *g, *splus, *sminus;
    double x, eps = 1.0e-05;
    int n = gretl_matrix_rows(H);
    int i, j, err = 0;
    
    splus  = malloc(n * sizeof *splus);
    sminus = malloc(n * sizeof *sminus);
    g      = malloc(n * sizeof *g);

    if (splus == NULL || sminus == NULL || g == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<n; i++) {
	double b0 = b[i];

	b[i] = b0 + eps;
	err = gradfunc(b, g, n, cfunc, data);
	if (err) break;
	for (j=0; j<n; j++) {
	    splus[j] = g[j];
	}

	b[i] = b0 - eps;
	err = gradfunc(b, g, n, cfunc, data);
	if (err) break;
	for (j=0; j<n; j++) {
	    sminus[j] = g[j];
	}

	b[i] = b0;
	for (j=0; j<n; j++) {
	    x = -(splus[j] - sminus[j]) / (2*eps);
	    gretl_matrix_set(H, i, j, x);
	}
    }

    if (!err) {
	gretl_matrix_xtr_symmetric(H);
    }

 bailout:

    free(splus);
    free(sminus);
    free(g);

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

static void hess_b_adjust_i (double *c, const double *b, double *h, int n, 
			     int i, double sgn)
{
    memcpy(c, b, n * sizeof *b);
    c[i] += sgn * h[i];
}

static void hess_b_adjust_ij (double *c, const double *b, double *h, int n, 
			      int i, int j, double sgn)
{
    memcpy(c, b, n * sizeof *b);
    c[i] += sgn * h[i];
    c[j] += sgn * h[j];
}

#define RSTEPS 4

/* The algorithm below implements the method of Richardson
   Extrapolation.  It is derived from code in the gnu R package
   "numDeriv" by Paul Gilbert, which was in turn derived from code
   by Xinqiao Liu.  Turned into C and modified for gretl by
   Allin Cottrell, June 2006.  On successful completion, writes
   the negative inverse of the Hessian into @H.
*/

static int numerical_hessian (const double *b, gretl_matrix *H,
			      BFGS_CRIT_FUNC func, void *data)
{
    double Dx[RSTEPS];
    double Hx[RSTEPS];
    double *wspace;
    double *c, *h0, *h, *Hd, *D;
    /* numerical parameters */
    int r = RSTEPS;      /* number of Richardson steps */
    double eps = 1.0e-4;
    double d = 0.0001;
    double v = 2.0;      /* reduction factor for h */
    double f0, f1, f2;
    double p4m, hij;
    int n = gretl_matrix_rows(H);
    int vn = (n * (n + 1)) / 2;
    int dn = vn + n;
    int i, j, k, m, u;
    int err = 0;

    wspace = malloc((4 * n + dn) * sizeof *wspace);
    if (wspace == NULL) {
	return E_ALLOC;
    }

    c = wspace;
    h0 = c + n;
    h = h0 + n;
    Hd = h + n;
    D = Hd + n; /* D is of length dn */

    for (i=0; i<n; i++) {
	h0[i] = (fabs(b[i]) < 0.01)? eps : d * b[i];
    }

    f0 = func(b, data);

    /* first derivatives and Hessian diagonal */

    for (i=0; i<n; i++) {
	hess_h_init(h, h0, n);
	for (k=0; k<r; k++) {
	    hess_b_adjust_i(c, b, h, n, i, 1);
	    f1 = func(c, data);
	    if (na(f1)) {
		fprintf(stderr, "numerical_hessian: 1st derivative: "
			"criterion = NA for theta[%d] = %g\n", i, c[i]);
		err = E_NAN;
		goto bailout;
	    }
	    hess_b_adjust_i(c, b, h, n, i, -1);
	    f2 = func(c, data);
	    if (na(f2)) {
		fprintf(stderr, "numerical_hessian: 1st derivative: "
			"criterion = NA for theta[%d] = %g\n", i, c[i]);
		err = E_NAN;
		goto bailout;
	    }
	    /* F'(i) */
	    Dx[k] = (f1 - f2) / (2.0 * h[i]); 
	    /* F''(i) */
	    Hx[k] = (f1 - 2.0*f0 + f2) / (h[i] * h[i]);
	    hess_h_reduce(h, v, n);
	}
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
	for (j=0; j<=i; j++) {
	    if (i == j) {
		D[u] = Hd[i];
	    } else {
		hess_h_init(h, h0, n);
		for (k=0; k<r; k++) {
		    hess_b_adjust_ij(c, b, h, n, i, j, 1);
		    f1 = func(c, data);
		    if (na(f1)) {
			fprintf(stderr, "numerical_hessian: 2nd derivatives (%d,%d): "
				"objective function gave NA\n", i, j);
			err = E_NAN;
			goto bailout;
		    }
		    hess_b_adjust_ij(c, b, h, n, i, j, -1);
		    f2 = func(c, data);
		    if (na(f2)) {
			fprintf(stderr, "numerical_hessian: 2nd derivatives (%d,%d): "
				"objective function gave NA\n", i, j);
			err = E_NAN;
			goto bailout;
		    }
		    /* cross-partial */
		    Dx[k] = (f1 - 2.0*f0 + f2 - Hd[i]*h[i]*h[i]
			     - Hd[j]*h[j]*h[j]) / (2.0*h[i]*h[j]);
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
	    }
	    u++;
	}
    }

    /* transcribe the negative of the Hessian */
    u = n;
    for (i=0; i<n; i++) {
	for (j=0; j<=i; j++) {
	    hij = -D[u++];
	    gretl_matrix_set(H, i, j, hij);
	    gretl_matrix_set(H, j, i, hij);
	}
    }

 bailout:

    if (err && err != E_ALLOC) {
	gretl_errmsg_set(_("Failed to compute numerical Hessian"));
    }

    free(wspace);

    return err;
}

/**
 * numerical_hessian_inverse:
 * @b: array of parameter estimates.
 * @n: the number of elements in @b.
 * @func: function to compute criterion.
 * @data: data to be passed to the @gradfunc callback.
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
					 void *data, int *err)
{
    gretl_matrix *H = gretl_zero_matrix_new(n, n);

    if (H == NULL) {
	*err = E_ALLOC;
    } else {
	*err = numerical_hessian(b, H, func, data);
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
	return numerical_hessian(b, H, cfunc, data);
    }
}

#define ALT_OPG 0

/* build the T x k G matrix, given a set of coefficient estimates, 
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

    for (i=0; i<k; i++) {
	bi0 = b[i];
#if ALT_OPG
	h = d * bi0 + d * (b[i] == 0.0);
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
    double v = 2.0;
    double h, p4m;
    double bi0, f1, f2;
    int r = RSTEPS;
    int i, k, m;

    for (i=0; i<n; i++) {
	bi0 = b[i];
	h = d * b[i] + eps * (b[i] == 0.0);
	for (k=0; k<r; k++) {
	    b[i] = bi0 - h;
	    f1 = func(b, data);
	    b[i] = bi0 + h;
	    f2 = func(b, data);
	    if (na(f1) || na(f2)) {
		b[i] = bi0;
		return 1;
	    }		    
	    df[k] = (f2 - f1) / (2.0 * h); 
	    h /= v;
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

#define STEPFRAC	0.2
#define acctol		1.0e-7 /* alt: 0.0001 or 1.0e-7 (?) */
#define reltest		10.0

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

/* 
   If "set initvals" has been used, replace whatever initial values
   might have been in place with those given by the user (the customer
   is always right).  In addition, respect user settings for the
   maximum number of iterations, the convergence tolerance and
   so on.
*/

static void optim_get_user_values (double *b, int n, int *maxit,
				   double *reltol, double *gradmax,
				   int *quad, gretlopt opt, PRN *prn)
{
    const gretl_matrix *uinit;
    int uilen, umaxit;
    double utol;
    int i;

    /* we first check to see if we've been a usable initialization
       for the parameter estimates */

    uinit = get_init_vals();
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
	    if (opt & OPT_V) {
		pputs(prn, _("\n\n*** User-specified starting values:\n"));
		for (i=0; i<n; i++) {
		    pprintf(prn, " %12.6f", b[i]);
		    if (i % 6 == 5) {
			pputc(prn, '\n');
		    }
		}
		pputs(prn, "\n\n");
	    }
	    free_init_vals();
	}
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

    /* step-length algorithm (BFGS only at present) */
    if (quad != NULL) {
	if (libset_get_int(OPTIM_STEPLEN) == STEPLEN_QUAD) {
	    *quad = 1;
	}
    }
}

#define bfgs_print_iter(v,s,i) (v && (s == 1 || i % s == 0))

#define GRAD_TOLER 1.0

static int copy_initial_hessian (gretl_matrix *A, double **H, 
				 int n)
{
    int i, j;

#if BFGS_DEBUG
    gretl_matrix_print(A, "BFGS: initial Hessian inverse");
#endif

    if (gretl_is_null_matrix(A)) {
	for (i=0; i<n; i++) {
	    for (j=0; j<i; j++) {
		H[i][j] = 0.0;
	    }
	    H[i][i] = 1.0;
	}
    } else if (A->rows != n || A->cols != n) {
	return E_NONCONF;
    } else {
	for (i=0; i<n; i++) {
	    for (j=0; j<=i; j++) {
		H[i][j] = gretl_matrix_get(A, i, j);
	    }
	}	
    }

    return 0;
}

/* returns number of coefficient that have actually changed */

static int coeff_at_end (double *b, const double *X, const double *t, 
			 int n, double length)
{
    int i, ndelta = n;

    for (i=0; i<n; i++) {
	b[i] = X[i] + length * t[i];
	if (coeff_unchanged(b[i], X[i])) {
	    ndelta--;
	}
    }

    return ndelta;
}

static double quad_slen (int n, int *pndelta, double *b, 
			 const double *X, const double *t, 
			 double *pf, BFGS_CRIT_FUNC cfunc, void *data, 
			 double g0, double f0, int *pfcount)
{
    double d, f1 = *pf;
    double steplen = 1.0, endpoint = 1.0;
    int ndelta, crit_ok = 0, fcount = 0, f1_done = 0;
    double safelen = 1.0e-12;
    double incredible = -1.0e12;

    /* Below: iterate so long as (a) we haven't achieved an acceptable
       value of the criterion and (b) there is still some prospect
       of doing so.
    */    

    do {
	crit_ok = 0;
	ndelta = coeff_at_end(b, X, t, n, endpoint);

	if (ndelta > 0) {
	    if (!f1_done) {
		f1 = cfunc(b, data);
		fcount++;
	    }
	    d = -g0 * endpoint * acctol;

	    /* find the optimal steplength by quadratic interpolation; 
	       inspired by Kelley (1999), "Iterative Methods for Optimization", 
	       especially section 3.2.1. 
	    */
	    if (xna(f1)) {
		/* function goes into NA zone, presumably outside the 
		   admissible parameter space; hence, try a much smaller 
		   step. FIXME execution can come back here indefinitely.
		*/
		endpoint *= STEPFRAC;
#if BFGS_DEBUG
		fprintf(stderr, "quad_slen: f1 is NA; trimming\n");
#endif
	    } else if ((f1 - f0) < incredible) {
		/* Same as above, with the exception that the objective
		   function technically computes, but returns a fishy value.
		*/
		endpoint *= STEPFRAC;
#if BFGS_DEBUG
		fprintf(stderr, "opt_slen: %g is incredible; trimming\n", 
			f1 - f0);
#endif
	    } else if (f1 < f0 + d) {
		/* function computes, but goes down: try quadratic approx */
		steplen = 0.5 * endpoint * g0 / (f0 - f1 + g0);
#if BFGS_DEBUG
		fprintf(stderr, "quad_slen, interpolate: f0 = %g, f1 = %g, g0 = %g\n", 
			f0, f1, g0);
		fprintf(stderr, "quad_slen, interpolate: endpoint = %g, "
			"steplen = %g\n", endpoint, steplen);
#endif
		    
		if (steplen < safelen) {
		    /* We have a ludicrously small steplength here,
		       most likely because the endpoint is too far out.
		       Let's trim it down and retry. 
		    */
		    endpoint *= STEPFRAC;
		} else {
		    ndelta = coeff_at_end(b, X, t, n, steplen);
		    f1 = cfunc(b, data);
		    fcount++;		    
#if BFGS_DEBUG
		    fprintf(stderr, "quad_slen, interpolate: %g is safe\n", steplen); 
#endif
		    crit_ok = !na(f1) && (f1 >= f0 + d);
		    /* if the function still goes down (or berserk), let's 
		       trim the endpoint one more time and retry */
#if BFGS_DEBUG
		    fprintf(stderr, "quad_slen, interpolate: crit_ok = %d"
			    " (f1 = %g)\n", crit_ok, f1);
#endif
		    if (!crit_ok) {
			endpoint = steplen;
			f1_done = 1;
		    } else {
			f1_done = 0;
		    }
		}
	    } else {
		crit_ok = 1;
		steplen = endpoint;
	    }
	}
    } while (ndelta > 0 && !crit_ok);

#if BFGS_DEBUG
    fprintf(stderr, "opt_slen: steplen = %g\n", steplen);
#endif
    *pndelta = ndelta;
    *pfcount += fcount;
    *pf = f1;

    return steplen;
}

static double simple_slen (int n, int *pndelta, double *b, double *X, double *t, 
			   double *pf, BFGS_CRIT_FUNC cfunc, void *data, 
			   double g0, double f0, int *pfcount)
{
    double d, f1 = *pf, steplen = 1.0;
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
	    f1 = cfunc(b, data);
	    d = -g0 * steplen * acctol;
	    fcount++;
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

static int BFGS_orig (double *b, int n, int maxit, double reltol,
		      int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
		      int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
		      gretl_matrix *A0, gretlopt opt, PRN *prn)
{
    int verbskip, verbose = (opt & OPT_V);
    double *wspace = NULL;
    double **H = NULL;
    double *g, *t, *X, *c;
    int fcount, gcount, ndelta = 0;
    int quad = 0, show_activity = 0;
    double sumgrad, gradmax, gradnorm = 0.0;
    double fmax, f, f0, s, steplen = 0.0;
    double D1, D2;
    int i, j, ilast, iter, done;
    int err = 0;

    optim_get_user_values(b, n, &maxit, &reltol, &gradmax, &quad, opt, prn);

    if (gradfunc == NULL) {
	gradfunc = numeric_gradient;
    }

    wspace = malloc(4 * n * sizeof *wspace);
    H = triangular_array_new(n);

    if (wspace == NULL || H == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* initialize curvature matrix */
    err = copy_initial_hessian(A0, H, n);
    if (err) {
	goto bailout;
    }

    g = wspace;
    t = g + n;
    X = t + n;
    c = X + n;

    f = cfunc(b, data);

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
    gradfunc(b, g, n, cfunc, data);

#if BFGS_DEBUG
    fprintf(stderr, "initial gradient:\n");
    for (i=0; i<n; i++) {
	fprintf(stderr, " g[%d] = %g\n", i, g[i]);
    }
#endif

    if (maxit == 0) {
	goto skipcalc;
    }

    verbskip = libset_get_int("bfgs_verbskip");
    show_activity = show_activity_func_installed();

    do {
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
	}

	gradnorm = sqrt(gradnorm / n);

#if BFGS_DEBUG
	fprintf(stderr, "\niter %d: sumgrad = %g\n", iter, sumgrad);
#endif

#if BFGS_DEBUG
	fprintf(stderr, "H = \n");
	for (i=0; i<n; i++) {
	    for (j=0; j<=i; j++) {
		fprintf(stderr, "%15.6f", H[i][j]);
	    }
	    fputc('\n', stderr);
	}
#endif
	if (sumgrad > 0.0) { 
	    /* heading in the right direction */
	    if (quad) {
		steplen = quad_slen(n, &ndelta, b, X, t, &f, cfunc, data, 
				    sumgrad, fmax, &fcount);
	    } else {
		steplen = simple_slen(n, &ndelta, b, X, t, &f, cfunc, data, 
				      sumgrad, fmax, &fcount);
	    }		
	    done = fabs(fmax - f) <= reltol * (fabs(fmax) + reltol);

#if BFGS_DEBUG
	    fprintf(stderr, "convergence test: LHS=%g, RHS=%g; done = %d\n",
		    fabs(fmax - f), reltol * (fabs(fmax) + reltol),
		    done);
#endif

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
		gradfunc(b, g, n, cfunc, data);
#if BFGS_DEBUG
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
#if BFGS_DEBUG
		fprintf(stderr, "D1 = %g\n", D1);
#endif
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
#if BFGS_DEBUG
		    fprintf(stderr, "D2 = %g\n", D2);
#endif
		} else {
		    /* D1 <= 0.0 */
		    ilast = gcount;
		}
	    } else if (ilast < gcount) {
		ndelta = n;
		ilast = gcount;
	    }
	} else if (sumgrad == 0.0) {
	    fprintf(stderr, "gradient is exactly zero!\n");
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

    } while (ndelta > 0 || ilast < gcount);

#if BFGS_DEBUG
    fprintf(stderr, "terminated: fmax=%g, ndelta=%d, ilast=%d, gcount=%d\n",
	    fmax, ndelta, ilast, gcount);
#endif

    if (iter >= maxit) {
	fprintf(stderr, _("stopped after %d iterations\n"), iter);
	err = E_NOCONV;
    } else if (gradnorm > gradmax) {
	err = E_NOCONV;
    } else if (fmax < f0) {
	/* allow a small sloppiness factor here? */
	double rdiff;

	rdiff = (f0 == 0.0)? -fmax : fabs((f0 - fmax) / f0);
	if (rdiff > 1.0e-12) {
	    fprintf(stderr, "failed to match initial value of objective function:\n"
		    " f0=%.18g, fmax=%.18g\n", f0, fmax);
	    err = E_NOCONV;
	}
    } 

    if (!err && gradnorm > GRAD_TOLER) {
	gretl_warnmsg_sprintf(_("norm of gradient = %g"), gradnorm);
	set_gretl_warning(W_GRADIENT);
    }

 skipcalc:

    *fncount = fcount;
    *grcount = gcount;

    if (verbose) {
	print_iter_info(-1, f, crittype, n, b, g, steplen, prn);
	pputc(prn, '\n');
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
   has the strength to go into lbfgsb.c and make the
   necessary adjustments.
*/

static void reverse_gradient (double *g, int n)
{
    int i;

    for (i=0; i<n; i++) {
	g[i] = -g[i];
    }
}

int LBFGS_max (double *b, int n, int maxit, double reltol,
	       int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	       int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	       gretlopt opt, PRN *prn)
{
    double *wspace = NULL;
    double *g, *l, *u, *wa;
    int *iwa, *nbd = NULL;
    int i, m, dim;
    char task[60];
    char csave[60];
    double f, pgtol;
    double gradmax;
    double dsave[29];
    int isave[44];
    int lsave[4];
    int iter, ibak = 0;
    int show_activity = 0;
    int verbskip, verbose = (opt & OPT_V);
    int err = 0;

    *fncount = *grcount = 0;    

    optim_get_user_values(b, n, &maxit, &reltol, &gradmax, NULL, opt, prn);

    /*
      m: the number of corrections used in the limited memory matrix.
      It is not altered by the routine.  Values of m < 3 are not
      recommended, and large values of m can result in excessive
      computing time. The range 3 <= m <= 20 is recommended.

      Was initially set to 5 (then 10, then 8; and 8 is the default).
    */
    m = libset_get_int(LBFGS_MEM); 

    dim = (2*m+4)*n + 12*m*m + 12*m; /* for wa */
    dim += 3*n; /* for g, l and u */

    wspace = malloc(dim * sizeof *wspace);
    nbd = malloc(4*n * sizeof *nbd);

    if (wspace == NULL || nbd == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    g = wspace;
    l = g + n;
    u = l + n;
    wa = u + n;
    iwa = nbd + n;

    verbskip = libset_get_int("bfgs_verbskip");
    show_activity = show_activity_func_installed();

    if (gradfunc == NULL) {
	gradfunc = numeric_gradient;
    }

    /* Gradient convergence criterion (not used -- we use reltol instead) */
    pgtol = 0;

    /* Bounds on the parameters: for now we just set them all to be
       less than some ridiculously large number */
    for (i=0; i<n; i++) {
	nbd[i] = 3; /* case 3: upper bound only */
	u[i] = NADBL / 100;
    }	

    /* Start the iteration by initializing 'task' */
    strcpy(task, "START");

    while (1) {
	/* Call the L-BFGS-B code */
	setulb_(n, m, b, l, u, nbd, &f, g, &reltol, &pgtol, wa, iwa, 
		task, csave, lsave, isave, dsave);

	iter = isave[29] + 1;

	if (!strncmp(task, "FG", 2)) {
	    /* Compute function value, f */
	    f = cfunc(b, data);
	    if (!na(f)) {
		f = -f; /* maximize, don't minimize */
	    } else if (*fncount == 0) {
		fprintf(stderr, "initial value of f is not finite\n");
		err = E_DATA;
		break;
	    }
	    *fncount += 1;
	    /* Compute gradient, g */
	    gradfunc(b, g, n, cfunc, data);
	    reverse_gradient(g, n);
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
	    break;
	}

	if (bfgs_print_iter(verbose, verbskip, iter)) {
	    if (iter != ibak) {
		double steplen = (iter == 1)? NADBL : dsave[13];

		reverse_gradient(g, n);
		print_iter_info(iter, -f, crittype, n, b, g, steplen, prn);
		reverse_gradient(g, n);
	    }
	    ibak = iter;
	}

	if (show_activity && (iter % 10 == 0)) {
	    show_activity_callback();
	}
    }

    if (!err && crittype == C_GMM) {
	/* finalize GMM computations */
	f = cfunc(b, data);
    }

    if (opt & OPT_V) {
	reverse_gradient(g, n);
	print_iter_info(-1, -f, crittype, n, b, g, dsave[13], prn);
	pputc(prn, '\n');
    }

 bailout:

    free(wspace);
    free(nbd);

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
 * criterion value as calculated by @cfunc.  Uses the BFGS
 * variable-metric method.  Based on Pascal code in J. C. Nash,
 * "Compact Numerical Methods for Computers," 2nd edition, converted
 * by p2c then re-crafted by B. D. Ripley for gnu R.  Revised for 
 * gretl by Allin Cottrell.
 * 
 * Returns: 0 on successful completion, non-zero error code
 * on error.
 */

int BFGS_max (double *b, int n, int maxit, double reltol,
	      int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	      int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	      gretl_matrix *A0, gretlopt opt, PRN *prn)
{
    int ret, wnum;

    if ((opt & OPT_L) || libset_get_bool(USE_LBFGS)) {
	ret = LBFGS_max(b, n, maxit, reltol,
			fncount, grcount, cfunc, 
			crittype, gradfunc, data, 
			opt, prn);
    } else {
	ret = BFGS_orig(b, n, maxit, reltol,
			fncount, grcount, cfunc, 
			crittype, gradfunc, data, 
			A0, opt, prn);
    }

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

static void umax_destroy (umax *u)
{
    if (u->dset != NULL) {
	/* drop any private "$" series created */
	dataset_drop_listed_variables(NULL, u->dset, NULL, NULL);
    }

    user_var_delete_by_name("$umax", NULL);

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

static int user_gen_setup (umax *u,
			   const char *fncall,
			   const char *gradcall,
			   const char *hesscall,
			   DATASET *dset)
{
    char formula[MAXLINE];
    int err = 0;

    if (u->gentype == GRETL_TYPE_MATRIX) {
	sprintf(formula, "matrix $umax=%s", fncall);
    } else {
	sprintf(formula, "$umax=%s", fncall);
    }

    u->gf = genr_compile(formula, dset, OPT_P, &err);

    if (!err) {
	/* see if the formula actually works */
	err = execute_genr(u->gf, dset, u->prn);
    }

    if (!err && gradcall != NULL) {
	/* process gradient formula */
	err = optimizer_get_matrix_name(gradcall, u->gmname);
	if (!err) {
	    u->gg = genr_compile(gradcall, dset, OPT_P | OPT_O, &err);
	    if (!err) {
		err = execute_genr(u->gg, dset, u->prn);
	    } 
	}
    }

    if (!err && hesscall != NULL) {
	/* process Hessian formula */
	err = optimizer_get_matrix_name(hesscall, u->hmname);
	if (!err) {
	    u->gh = genr_compile(hesscall, dset, OPT_P | OPT_O, &err);
	    if (!err) {
		err = execute_genr(u->gg, dset, u->prn);
	    } 
	}
    }

    if (!err) {
	u->dset = dset;
	u->fm_out = genr_get_output_matrix(u->gf);
    } else {
	destroy_genr(u->gf);
	destroy_genr(u->gg);
	u->gf = u->gg = NULL;
    }

    return err;
}

double user_BFGS (gretl_matrix *b, 
		  const char *fncall,
		  const char *gradcall, 
		  DATASET *dset,
		  PRN *prn, int *err)
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

    *err = user_gen_setup(u, fncall, gradcall, NULL, dset);
    if (*err) {
	return NADBL;
    }

    tol = libset_get_double(BFGS_TOLER);
    verbose = libset_get_bool(MAX_VERBOSE);

    if (verbose) {
	opt = OPT_V;
	u->prn = prn;
    }

    *err = BFGS_max(b->val, u->ncoeff, 
		    maxit, tol, &fcount, &gcount,
		    user_get_criterion, C_OTHER, 
		    (u->gg == NULL)? NULL : user_get_gradient, 
		    u, NULL, opt, prn);

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
	return NADBL;
    }

    if (libset_get_bool(MAX_VERBOSE)) {
	opt = OPT_V;
	u->prn = prn;
    }

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

double user_simann (gretl_matrix *b, 
		    const char *fncall,
		    int maxit,
		    DATASET *dset,
		    PRN *prn, int *err)
{
    umax *u;
    gretlopt opt = OPT_NONE;
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

    *err = user_gen_setup(u, fncall, NULL, NULL, dset);
    if (*err) {
	return NADBL;
    }

    if (libset_get_bool(MAX_VERBOSE)) {
	opt = OPT_V;
	u->prn = prn;
    }

    *err = gretl_simann(b->val, u->ncoeff, maxit, 
			user_get_criterion, u, 
			opt, prn);

    if (!*err) {
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

gretl_matrix *fdjac (gretl_matrix *theta, const char *fncall,
		     DATASET *dset, int *err)
{
    umax *u;
    gretl_matrix *J = NULL;
    int m, n;
    int iflag = 0;
    int quality = 0;
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
    }

    quality = libset_get_int(FDJAC_QUAL);

    fdjac2_(user_calc_fvec, m, n, quality, theta->val, fvec, J->val, 
	    m, &iflag, 0.0, wa, u);

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

/* Below: Newton-Raphson code, starting with a few
   auxiliary functions */

static int broken_matrix (const gretl_matrix *m)
{
    int i, n = m->rows * m->cols;

    for (i=0; i<n; i++) {
	if (xna(m->val[i])) {
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
    } else if (status == STEPMIN_MET) {
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

#define SPECTRAL 0

static int NR_invert_hessian (gretl_matrix *H, const gretl_matrix *Hcpy)
{
    double hii, hmin = 1.0e-28; /* was 1.0e-20 */
    int i, j, err = 0;
    int n = H->rows;
    int restore = 0;
    double x;

    /* first, check if all the elements along the diagonal are 
       numerically positive
    */

    for (i=0; i<n && !err; i++) {
	hii = gretl_matrix_get(H, i, i);
	err = hii < hmin;
    }

    if (err) {
	fprintf(stderr, "NR_invert_hessian: non-positive "
		"diagonal (hii=%g)\n", hii);
#if 1
	gretl_matrix_print(H, "H");
#endif
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
#if 0
		gretl_matrix_print(H, "after shrinkage");
#endif
		err = gretl_invert_symmetric_matrix(H);
	    }
	    restore = (err != 0);
	}
    }


#if SPECTRAL
    if (err) {
	gretl_matrix *evecs;
	gretl_matrix *evals;
	gretl_matrix *tmp;
	double y;

	fprintf(stderr, "newton hessian fixup: spectral method\n");

	evecs = gretl_matrix_copy(Hcpy);
	tmp = gretl_matrix_alloc(n, n);
	evals = gretl_symmetric_matrix_eigenvals(evecs, 1, &err);

	if (!err) {
	    for (i=0; i<n; i++) {
		x = 1.0 / (1.0 + fabs(gretl_vector_get(evals, i)));
		for (j=0; j<n; j++) {
		    y = x * gretl_matrix_get(evecs, j, i);
		    gretl_matrix_set(tmp, i, j, y);
		}
	    }

	    gretl_matrix_multiply(evecs, tmp, H);
	    gretl_matrix_print(H, "after spectral fixup");
	}

	gretl_matrix_free(evals);
	gretl_matrix_free(evecs);
	gretl_matrix_free(tmp);
    }
#endif

    if (err) {
	fprintf(stderr, "newton hessian fixup: err = %d -> desperation!\n",
		err);
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
    int verbose = (opt & OPT_V);
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

    /* needs some work */
    optim_get_user_values(b, n, NULL, NULL, NULL, NULL, opt, prn);

    b1 = b0 + n;
    copy_to(b1, b, n);

    f1 = cfunc(b1, data);
    if (na(f1)) {
	gretl_errmsg_set(_("Initial value of objective function is not finite"));
	err = E_NAN;
    }

    if (!err) {
	if (gradfunc != NULL) {
	    err = gradfunc(b, g->val, n, cfunc, data);
	} else {
	    err = numeric_gradient(b, g->val, n, cfunc, data);
	}
    }

    if (!err) {
	if (hessfunc != NULL) {
	    err = hessfunc(b, H1, data);
	} else {
	    err = NR_fallback_hessian(b, H1, gradfunc, cfunc, data);
	}
	if (!err) {
	    gretl_matrix_copy_values(H0, H1);
	    err = NR_invert_hessian(H1, H0);
	}
    }

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
	f1 = cfunc(b1, data);

	while (na(f1) || ((f1 < f0) && (steplen > stepmin))) {
	    /* try smaller step */
	    steplen /= 2.0;
	    copy_plus(b1, b0, steplen, a->val, n);
	    f1 = cfunc(b1, data);
	}

	if (verbose) {
	    print_iter_info(iter, f1, crittype, n, b1, g->val, 
			    steplen, prn);
	}

	if (gradfunc != NULL) {
	    err = gradfunc(b1, g->val, n, cfunc, data);
	} else {
	    err = numeric_gradient(b1, g->val, n, cfunc, data);
	}

	if (err || broken_matrix(g)) {
	    err = (err == 0)? E_NAN : err;
	    break;
	}	

	if (hessfunc != NULL) {
	    err = hessfunc(b1, H1, data);
	} else {
	    err = NR_fallback_hessian(b1, H1, gradfunc, cfunc, data);
	}

	if (err || broken_matrix(H1)) {
	    err = (err == 0)? E_NAN : err;
	    break;
	}

	if (!err) {
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

static void set_up_matrix (gretl_matrix *m, double *val, 
			   int rows, int cols)
{
    m->val = val;
    m->rows = rows;
    m->cols = cols;
    m->info = NULL;
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
    int i, err = 0;

    if (maxit <= 0) {
	maxit = 1024;
    }

    set_up_matrix(&b, theta, n, 1);

    b0 = gretl_matrix_copy(&b);
    b1 = gretl_matrix_copy(&b);
    bstar = gretl_matrix_copy(&b);
    d = gretl_column_vector_alloc(n);

    if (b0 == NULL || b1 == NULL || bstar == NULL || d == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    f0 = fbest = fworst = cfunc(b.val, data);

    if (opt & OPT_V) {
	pprintf(prn, "\nSimulated annealing: initial function value = %.8g\n",
		f0);
    }

    /* Question: should the initial radius be a function of the scale
       of the initial parameter vector?
    */

    for (i=0; i<maxit; i++) {
	gretl_matrix_random_fill(d, D_NORMAL);
	gretl_matrix_multiply_by_scalar(d, radius);
	gretl_matrix_add_to(b1, d);
	f1 = cfunc(b1->val, data);

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
	pprintf(prn, "No improvement found in %d iterations\n\n", maxit);
    }
    
    if (fbest - fworst < 1.0e-9) {
	pprintf(prn, "*** warning: surface seems to be flat\n");
    }

 bailout:

    gretl_matrix_free(b0);
    gretl_matrix_free(b1);
    gretl_matrix_free(bstar);
    gretl_matrix_free(d);

    return err;
}
