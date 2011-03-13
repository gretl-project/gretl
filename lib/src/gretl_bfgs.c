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

/* The BFGS optimizer and the fdjac Jacobian function */

#include "libgretl.h"
#include "gretl_bfgs.h"
#include "libset.h"
#include "matrix_extra.h"
#include "usermat.h"
#include "gretl_scalar.h"

#include "gretl_f2c.h"
#include "../../minpack/minpack.h"  

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

/* In case an analytical score function (@gradfun) is available,
   construct a numerical approximation to the Hessian using that
   function.  This is intended for building a covariance matrix at
   convergence; note that it will be unlikely to work at an arbitrary
   point in the parameter space.  

   Also note that the NULL pointer passed to gradfun in the 4th 
   argument corresponds to the BFGS_CRIT_FUNC parameter, so this 
   will fail horribly unless @gradfun calculates the score independently 
   of the criterion function (and so does not try to access the 4th 
   argument).
*/

gretl_matrix *hessian_from_score (double *b, int n, 
				  BFGS_GRAD_FUNC gradfun, 
				  void *data, int *err)
{
    gretl_matrix *H = NULL;
    double *g, *splus, *sminus;
    double x, eps = 1.0e-05;
    int i, j;
    
    splus  = malloc(n * sizeof *splus);
    sminus = malloc(n * sizeof *sminus);
    g      = malloc(n * sizeof *g);

    if (splus == NULL || sminus == NULL || g == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    H = gretl_matrix_alloc(n, n);
    if (H == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<n; i++) {
	double b0 = b[i];

	b[i] = b0 + eps;
	*err = gradfun(b, g, n, NULL, data);
	if (*err) break;
	for (j=0; j<n; j++) {
	    splus[j] = g[j];
	}

	b[i] = b0 - eps;
	*err = gradfun(b, g, n, NULL, data);
	if (*err) break;
	for (j=0; j<n; j++) {
	    sminus[j] = g[j];
	}

	b[i] = b0;
	for (j=0; j<n; j++) {
	    x = -(splus[j] - sminus[j]) / (2*eps);
	    gretl_matrix_set(H, i, j, x);
	}
    }

    if (!*err) {
	gretl_matrix_xtr_symmetric(H);
	*err = gretl_invert_symmetric_matrix(H);
    }

 bailout:

    if (*err) {
	gretl_matrix_free(H);
	H = NULL;
    }

    free(splus);
    free(sminus);
    free(g);

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
   "numDeriv" by Paul Gilbert, which was in turn derived from C code
   by Xinqiao Liu.  Turned back into C and modified for gretl by
   Allin Cottrell, June 2006.  On successful completion, returns
   the negative inverse of the Hessian.
*/

gretl_matrix *numerical_hessian (const double *b, int n, 
				 BFGS_CRIT_FUNC func, 
				 void *data, int *err)
{
    gretl_matrix *H = NULL;
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
    int vn = (n * (n + 1)) / 2;
    int dn = vn + n;
    int i, j, k, m, u;

    wspace = malloc((4 * n + dn) * sizeof *wspace);
    if (wspace == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    c = wspace;
    h0 = c + n;
    h = h0 + n;
    Hd = h + n;
    D = Hd + n; /* D is of length dn */

    H = gretl_matrix_alloc(n, n);
    if (H == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }	

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
		*err = E_NAN;
		goto bailout;
	    }
	    hess_b_adjust_i(c, b, h, n, i, -1);
	    f2 = func(c, data);
	    if (na(f2)) {
		fprintf(stderr, "numerical_hessian: 1st derivative: "
			"criterion = NA for theta[%d] = %g\n", i, c[i]);
		*err = E_NAN;
		goto bailout;
	    }
	    /* F'(i) */
	    Dx[k] = (f1 - f2) / (2.0 * h[i]); 
	    /* F''(i) */
	    Hx[k] = (f1 - 2.0*f0 + f2) / (h[i] * h[i]);
	    hess_h_reduce(h, v, n);
	}
	p4m = 4;
	for (m=0; m<r-1; m++) {
	    for (k=0; k<r-m; k++) {
		Dx[k] = (Dx[k+1] * p4m - Dx[k]) / (p4m - 1);
		Hx[k] = (Hx[k+1] * p4m - Hx[k]) / (p4m - 1);
	    }
	    p4m *= 4;
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
			*err = E_NAN;
			goto bailout;
		    }
		    hess_b_adjust_ij(c, b, h, n, i, j, -1);
		    f2 = func(c, data);
		    if (na(f2)) {
			fprintf(stderr, "numerical_hessian: 2nd derivatives (%d,%d): "
				"objective function gave NA\n", i, j);
			*err = E_NAN;
			goto bailout;
		    }
		    /* cross-partial */
		    Dx[k] = (f1 - 2.0*f0 + f2 - Hd[i]*h[i]*h[i]
			     - Hd[j]*h[j]*h[j]) / (2.0*h[i]*h[j]);
		    hess_h_reduce(h, v, n);
		}
		p4m = 4.0;
		for (m=0; m<r-1; m++) {
		    for (k=0; k<r-m; k++) {
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

    *err = gretl_invert_symmetric_matrix(H);

 bailout:

    if (*err != 0 && *err != E_ALLOC) {
	gretl_errmsg_set(_("Failed to compute numerical Hessian"));
    }

    free(wspace);

    return H;
}

#define ALT_OPG 0

/* build the G matrix, given a set of coefficient estimates, b, and a
   function for calculating the per-observation contributions
   to the loglikelihood, lltfun 
*/

gretl_matrix *build_score_matrix (double *b, int k, int T,
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

    G = gretl_zero_matrix_new(k, T);
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
	    gretl_matrix_set(G, i, t, x[t]);
	}
	b[i] = bi0 + h;
	x = lltfun(b, i, data);
	if (x == NULL) {
	    *err = E_NAN;
	    goto bailout;
	}
	for (t=0; t<T; t++) {
	    x0 = gretl_matrix_get(G, i, t);
	    gretl_matrix_set(G, i, t, (x[t] - x0) / (2.0 * h));
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
    int err = 0;

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
	    for (k=0; k<r-m; k++) {
		df[k] = (df[k+1] * p4m - df[k]) / (p4m - 1.0);
	    }
	    p4m *= 4.0;
	}
	g[i] = df[0];
    }

    return err;
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

int BFGS_numeric_gradient (double *b, double *g, int n,
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
    fprintf(stderr, "BFGS_numeric_gradient returning, err = %d\n", err);
#endif

    return err;
}

#define stepfrac	0.2
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
   maximum number of iterations and the convergence tolerance.
*/

static void BFGS_get_user_values (double *b, int n, int *maxit,
				  double *reltol, double *gradmax,
				  gretlopt opt, PRN *prn)
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

    /* then check for a setting of the maximum number
       of iterations */

    umaxit = libset_get_int(BFGS_MAXITER);
    if (umaxit >= 0) {
	*maxit = umaxit;
    } else if (*maxit < 0) {
	*maxit = BFGS_MAXITER_DEFAULT;
    }

    /* and then the convergence tolerance */

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

    /* and the maximum acceptable gradient norm */

    *gradmax = libset_get_double(BFGS_MAXGRAD);
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

static int BFGS_orig (double *b, int n, int maxit, double reltol,
		      int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
		      int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
		      gretl_matrix *A0, gretlopt opt, PRN *prn)
{
    int crit_ok, done;
    double *wspace = NULL;
    double **H = NULL;
    double *g, *t, *X, *c;
    int verbskip, verbose = (opt & OPT_V);
    int fcount, gcount, ndelta = 0;
    int show_activity = 0;
    double sumgrad, gradmax, gradnorm = 0.0;
    double fmax, f, f0, d = 0.0;
    double s, steplen = 0.0;
    double D1, D2;
    int i, j, ilast, iter;
    int err = 0;

    BFGS_get_user_values(b, n, &maxit, &reltol, &gradmax, opt, prn);

    if (gradfunc == NULL) {
	gradfunc = BFGS_numeric_gradient;
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

	if (sumgrad > 0.0) { 
	    /* heading in the right direction */
	    steplen = 1.0;
	    crit_ok = 0;
	    do {
		/* loop so long as (a) we haven't achieved an
		   acceptable value of the criterion and (b) there is
		   still some prospect of doing so */
		ndelta = n;
		for (i=0; i<n; i++) {
		    b[i] = X[i] + steplen * t[i];
		    if (coeff_unchanged(b[i], X[i])) {
			ndelta--;
		    }
		}
		if (ndelta > 0) {
		    f = cfunc(b, data);
		    d = -sumgrad * steplen * acctol;
		    fcount++;
		    crit_ok = !na(f) && (f >= fmax + d);
#if BFGS_DEBUG
		    fprintf(stderr, "crit_ok: f=%.10g, fmax=%.10g, d=%g, steplen=%g, ok=%d\n",
			    f, fmax, d, steplen, crit_ok);
#endif
		    if (!crit_ok) {
			/* calculated criterion no good: try smaller step */
			steplen *= stepfrac;
		    }
		}
	    } while (ndelta != 0 && !crit_ok);

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
		/* we just reset: don't reset again; set ndelta = 0 so
		   that we exit the main loop
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
    } else if (gradnorm > gradmax || gradnorm > GRAD_TOLER) {
	if (gradnorm > gradmax) {
	    err = E_NOCONV;
	} else {
	    gretl_warnmsg_sprintf(_("norm of gradient = %g"), gradnorm);
	    set_gretl_warning(W_GRADIENT);
	}
    } else if (fmax < f0) {
	/* FIXME this should never happen */
	fprintf(stderr, "failed to match initial value of objective function:\n"
		" f0=%.15g, fmax=%.15g\n", f0, fmax);
	err = E_NOCONV;
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
    double *g = NULL;
    double *l = NULL;
    double *u = NULL;
    double *wa = NULL;
    int *nbd = NULL;
    int *iwa = NULL;
    int i, m, wadim;
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

    BFGS_get_user_values(b, n, &maxit, &reltol, &gradmax, opt, prn);

    /*
      m: the number of corrections used in the limited memory matrix.
      It is not altered by the routine.  Values of m < 3 are not
      recommended, and large values of m can result in excessive
      computing time. The range 3 <= m <= 20 is recommended.

      Was initially set to 5 (then 10, then 8; and 8 is the default).
    */
    m = libset_get_int(LBFGS_MEM); 

    wadim = (2*m+4)*n + 12*m*m + 12*m;

    g = malloc(n * sizeof *g);
    l = malloc(n * sizeof *l);
    u = malloc(n * sizeof *u);
    wa = malloc(wadim * sizeof *wa);
    nbd = malloc(n * sizeof *nbd);
    iwa = malloc(3*n * sizeof *iwa);

    if (g == NULL || l == NULL || u == NULL ||
	wa == NULL || nbd == NULL || iwa == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    verbskip = libset_get_int("bfgs_verbskip");
    show_activity = show_activity_func_installed();

    if (gradfunc == NULL) {
	gradfunc = BFGS_numeric_gradient;
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
	setulb_(&n, &m, b, l, u, nbd, &f, g, &reltol, &pgtol, wa, iwa, 
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

    free(g);
    free(l);
    free(u);
    free(wa);
    free(nbd);
    free(iwa);

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

/* user-level access to BFGS and fdjac */

typedef struct umax_ umax;

struct umax_ {
    GretlType gentype;    /* GRETL_TYPE_DOUBLE or GRETL_TYPE_MATRIX */
    gretl_matrix *b;      /* parameter vector */
    gretl_matrix *g;      /* gradient vector */
    int ncoeff;           /* number of coefficients */
    GENERATOR *gf;        /* for generating scalar or matrix result */
    GENERATOR *gg;        /* for generating gradient */
    double fx_out;        /* function double value */
    gretl_matrix *fm_out; /* function matrix value */
    char gmname[VNAMELEN]; /* name of user-defined gradient vector */
    double ***Z;          /* pointer to data array */
    DATAINFO *dinfo;      /* dataset info */
    PRN *prn;             /* optional printing struct */
};

static umax *umax_new (GretlType t)
{
    umax *u = malloc(sizeof *u);

    if (u != NULL) {
	u->gentype = t;
	u->b = NULL;
	u->g = NULL;
	u->ncoeff = 0;
	u->gf = NULL;
	u->gg = NULL;
	u->fx_out = NADBL;
	u->fm_out = NULL;
	u->gmname[0] = '\0';
	u->Z = NULL;
	u->dinfo = NULL;
	u->prn = NULL;
    }

    return u;
}

static void umax_destroy (umax *u)
{
    if (u->dinfo != NULL) {
	/* drop any private "$" series created */
	dataset_drop_listed_variables(NULL, u->Z, u->dinfo, NULL, NULL);
    }

    gretl_scalar_delete("$umax", NULL);
    user_matrix_destroy_by_name("$umax", NULL);

    destroy_genr(u->gf);
    destroy_genr(u->gg);

    free(u);
}

/* user-defined BFGS: get the criterion value */

static double user_get_criterion (const double *b, void *p)
{
    umax *u = (umax *) p;
    double x = NADBL;
    int i, t, err;

    for (i=0; i<u->ncoeff; i++) {
	u->b->val[i] = b[i];
    }

    err = execute_genr(u->gf, u->Z, u->dinfo, u->prn); 

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

/* user-defined BFGS: get the gradient, if specified */

static int user_get_gradient (double *b, double *g, int k,
			      BFGS_CRIT_FUNC func, void *p)
{
    umax *u = (umax *) p;
    gretl_matrix *ug;
    int i, err;

    for (i=0; i<k; i++) {
	u->b->val[i] = b[i];
    }

    err = execute_genr(u->gg, u->Z, u->dinfo, u->prn); 

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

/* parse the name of the user gradient matrix (vector) out of
   the gradient function call, where it must be the first
   argument, given in pointer form
*/

static int get_grad_vector_name (umax *u, const char *gradcall)
{
    const char *s = strchr(gradcall, '(');
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
		strncat(u->gmname, s, n);
	    }
	}
    }

    return err;
}

static int user_gen_setup (umax *u,
			   const char *fncall,
			   const char *gradcall,
			   double ***pZ, 
			   DATAINFO *pdinfo)
{
    char formula[MAXLINE];
    int err = 0;

    if (u->gentype == GRETL_TYPE_MATRIX) {
	sprintf(formula, "matrix $umax=%s", fncall);
    } else {
	sprintf(formula, "$umax=%s", fncall);
    }

    u->gf = genr_compile(formula, pZ, pdinfo, OPT_P, &err);

    if (!err) {
	/* see if the formula actually works */
	err = execute_genr(u->gf, pZ, pdinfo, u->prn);
    }

    if (!err && gradcall != NULL) {
	/* process gradient formula */
	err = get_grad_vector_name(u, gradcall);
	if (!err) {
	    u->gg = genr_compile(gradcall, pZ, pdinfo, OPT_P | OPT_U, &err);
	    if (!err) {
		err = execute_genr(u->gg, pZ, pdinfo, u->prn);
	    } 
	}
    }

    if (!err) {
	u->Z = pZ;
	u->dinfo = pdinfo;
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
		  double ***pZ, DATAINFO *pdinfo,
		  PRN *prn, int *err)
{
    umax *u;
    double ret = NADBL;
    gretlopt opt = OPT_NONE;
    int verbose = libset_get_bool(MAX_VERBOSE);
    int maxit = BFGS_MAXITER_DEFAULT;
    int fcount = 0, gcount = 0;
    double tol;

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

    *err = user_gen_setup(u, fncall, gradcall, pZ, pdinfo);
    if (*err) {
	return NADBL;
    }

    tol = libset_get_double(BFGS_TOLER);

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

#define JAC_DEBUG 0

static int user_calc_fvec (integer *m, integer *n, double *x, double *fvec,
			   integer *iflag, void *p)
{
    umax *u = (umax *) p;
    gretl_matrix *v;
    int i, err;

    for (i=0; i<*n; i++) {
	u->b->val[i] = x[i];
    }

#if JAC_DEBUG
    gretl_matrix_print(u->b, "user_calc_fvec: u->b");
#endif

    err = execute_genr(u->gf, u->Z, u->dinfo, u->prn); 
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

    if (v == NULL || gretl_vector_get_length(v) != *m) {
	fprintf(stderr, "user_calc_fvec: got bad matrix\n"); 
	*iflag = -1;
    } else {
	for (i=0; i<*m; i++) {
	    fvec[i] = v->val[i];
	}
    }
    
    return 0;
}

static int fdjac_allocate (integer m, integer n,
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
		     double ***pZ, DATAINFO *pdinfo,
		     int *err)
{
    umax *u;
    gretl_matrix *J = NULL;
    integer m, n;
    integer iflag = 0;
    double *wa = NULL;
    double *fvec = NULL;
    double epsfcn = 0.0;
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

    *err = user_gen_setup(u, fncall, NULL, pZ, pdinfo);
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

    fdjac2_(user_calc_fvec, &m, &n, theta->val, fvec, J->val, 
	    &m, &iflag, &epsfcn, wa, u);

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

/* Below: experimental Newton-Raphson code, starting with a few
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
    GRADTOL_MET = 1 << 0,
    CRITTOL_MET = 1 << 1,
    STEPMIN_MET = 1 << 2
};

static void print_NR_status (int status, double crittol, double gradtol, 
			     PRN *prn)
{
    if (status == GRADTOL_MET) {
	pprintf(prn, "Gradient < %g; may be a solution\n", gradtol);
    } else if (status == CRITTOL_MET) {
	pprintf(prn, "Successive criterion values within tolerance (%g); "
		"may be a solution\n", crittol);
    } else if (status == STEPMIN_MET) {
	pprintf(prn, "Couldn't increase criterion; "
		"may be near a solution?\n");
    }
}

/* Newton-Raphson maximizer, loosely based on R's maxNR(). 

   The functions cfunc (computes the criterion, usually a
   loglikelihood) and gradfunc (score, provides an updated
   estimate of the gradient in its second argument) are just as in
   BFGS_max above.

   The precise working of the function hessfunc is debatable.  As
   things stand below, this function should write into its second
   argument the negative inverse of the Hessian (or some variant
   thereof).  This pushes onto hessfunc the problem of what to do if
   the observed Hessian is not negative definite.  

   The alternative would be to get hessfunc to provide the Hessian
   itself, negative definite or not, and for newton_raphson_max to
   take on the task of fixing it up if need be (which is what R's
   maxNR does). But it seemed to me that the best fix might depend on
   the particular application.  For example if the full gradient
   matrix, G, is available one might substitute (GG')^{-1}, hence
   falling back on BHHH temporarily.

   Note that unlike maxNR, there is no provision here for automatic
   substitution of numerical routines for gradfunc or hessfunc if
   these functions are not provided. My guess is that that without
   analytical functions Newton-Raphson is bound to perform worse than
   BFGS.
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
    gretl_matrix *g0, *g1;
    gretl_matrix *a;
    double *b0, *b1;
    double f0, f1;
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
			       &g0, n, 1,
			       &g1, n, 1,
			       &a, n, 1,
			       NULL);
    if (B == NULL) {
	free(b0);
	return E_ALLOC;
    }

    b1 = b0 + n;
    copy_to(b1, b, n);

    f1 = cfunc(b1, data);
    if (na(f1)) {
	gretl_errmsg_set(_("Initial value of objective function is not finite"));
	err = E_NAN;
    }

    if (!err) {
	err = gradfunc(b, g1->val, n, cfunc, data);
    }

    if (!err) {
	err = hessfunc(b, H1, data);
    }

    while (status == 0 && !err) {
	iter++;
	steplen = 1.0;
	f0 = f1;

	copy_to(b0, b1, n);

	gretl_matrix_copy_values(g0, g1);
	if (broken_matrix(g0)) {
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
	gretl_matrix_multiply(H0, g0, a);
	copy_plus(b1, b0, steplen, a->val, n);
	f1 = cfunc(b1, data);

	while (na(f1) || ((f1 < f0) && (steplen > stepmin))) {
	    /* try smaller step */
	    steplen /= 2.0;
	    copy_plus(b1, b0, steplen, a->val, n);
	    f1 = cfunc(b1, data);
	}

	if (verbose) {
	    print_iter_info(iter, f1, crittype, n, b1, g1->val, 
			    steplen, prn);
	}

	err = gradfunc(b1, g1->val, n, cfunc, data);
	if (err || broken_matrix(g1)) {
	    err = (err == 0)? E_NAN : err;
	    break;
	}	

	err = hessfunc(b1, H1, data);
	if (err || broken_matrix(H1)) {
	    err = (err == 0)? E_NAN : err;
	    break;
	}	

	if (steplen < stepmin) {
	    status = STEPMIN_MET;
	} else if (iter > maxit) {
	    err = E_NOCONV;
	} else if (sqrt(scalar_xpx(g1)) < gradtol) {
	    status = GRADTOL_MET;
	} else if (f1 - f0 < crittol) {
	    status = CRITTOL_MET;
	}
    }

    if (verbose) {
	print_iter_info(-1, f1, crittype, n, b1, g1->val, 
			steplen, prn);
	pputc(prn, '\n');
    }

    *itercount = iter;

    if (!err) {
	copy_to(b, b1, n);
	if (prn != NULL) {
	    print_NR_status(status, crittol, gradtol, prn);
	}
    }

    free(b0);
    gretl_matrix_block_destroy(B);

    return err;
}

/* Start by using BFGS with a sloppy tolerance, then switch to
   Newton-Raphson. An experimental thing.
*/

int hybrid_max (double *b, int n, int maxit, 
		double crittol, double gradtol, 
		int *itercount, int crittype, 
		BFGS_CRIT_FUNC cfunc,
		BFGS_GRAD_FUNC gradfunc, 
		HESS_FUNC hessfunc,
		void *data, gretlopt opt, 
		PRN *prn)
{
    double pre_tol = 1.0e-3;
    double pre_maxit = 100;
    int fncount = 0;
    int grcount = 0;
    int err;

    err = BFGS_max(b, n, pre_maxit, pre_tol, 
		   &fncount, &grcount, cfunc, crittype,
		   gradfunc, data, NULL, 
		   opt & OPT_V, prn);

    fprintf(stderr, "BFGS phase: fncount=%d, grcount=%d, err=%d\n",
	    fncount, grcount, err);

    if (!err) {
	err = newton_raphson_max(b, n, maxit, crittol, gradtol,
				 itercount, crittype, cfunc, 
				 gradfunc, hessfunc, data, 
				 opt & OPT_V, prn);
	fprintf(stderr, "Newton phase: itercount=%d, err=%d\n",
		*itercount, err);
    }

    *itercount += grcount;

    return err;
}
