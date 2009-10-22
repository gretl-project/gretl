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

#include "gretl_f2c.h"
#include "../../minpack/minpack.h"  

#define BFGS_DEBUG 0

void BFGS_defaults (int *maxit, double *tol, int ci)
{
    *maxit = libset_get_int(BFGS_MAXITER);
    *tol = libset_get_user_tolerance(BFGS_TOLER);

    if (ci != MLE && ci != GMM && *maxit <= 0) {
	*maxit = 1000;
    }

    if (ci == PROBIT || ci == INTREG || ci == ARMA) {
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
	    *maxit = 500;
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

static void hess_b_adjust_i (double *c, double *b, double *h, int n, 
			     int i, double sgn)
{
    int k;

    for (k=0; k<n; k++) {
	c[k] = b[k] + (k == i) * sgn * h[i];
    }
}

static void hess_b_adjust_ij (double *c, double *b, double *h, int n, 
			      int i, int j, double sgn)
{
    int k;

    for (k=0; k<n; k++) {
	c[k] = b[k] + (k == i) * sgn * h[i] +
	    (k == j) * sgn * h[j];
    }
}

#define RSTEPS 4

/* The algorithm below implements the method of Richardson
   Extrapolation.  It is derived from code in the gnu R package
   "numDeriv" by Paul Gilbert, which was in turn derived from C code
   by Xinqiao Liu.  Turned back into C and modified for gretl by
   Allin Cottrell, June 2006.
*/

double *numerical_hessian (const double *b, int n, BFGS_CRIT_FUNC func, 
			   void *data, int *err)
{
    double Dx[RSTEPS];
    double Hx[RSTEPS];
    double *bcpy, *wspace;
    double *c, *h0, *h, *Hd, *D;
    gretl_matrix *V = NULL;
    double *vcv = NULL;
    /* numerical parameters */
    int r = RSTEPS;      /* number of Richardson steps */
    double eps = 1.0e-4;
    double d = 0.0001;
    double v = 2.0;      /* reduction factor for h */
    double f0, f1, f2;
    double p4m;
    int vn = (n * (n + 1)) / 2;
    int dn = vn + n;
    int i, j, k, m, u;

    bcpy = copyvec(b, n);
    if (bcpy == NULL) {
	*err = E_ALLOC;
	return NULL;
    }	

    wspace = malloc((4 * n + dn) * sizeof *wspace);
    if (wspace == NULL) {
	free(bcpy);
	*err = E_ALLOC;
	return NULL;
    }

    c = wspace;
    h0 = c + n;
    h = h0 + n;
    Hd = h + n;
    D = Hd + n; /* D is of length dn */

    /* vech form of variance matrix */
    V = gretl_column_vector_alloc(vn);
    if (V == NULL) {
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
	    hess_b_adjust_i(c, bcpy, h, n, i, 1);
	    f1 = func(c, data);
	    if (na(f1)) {
		fprintf(stderr, "numerical_hessian: 1st derivative: "
			"criterion = NA for theta[%d] = %g\n", i, c[i]);
		*err = E_NAN;
		goto bailout;
	    }
	    hess_b_adjust_i(c, bcpy, h, n, i, -1);
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
		    hess_b_adjust_ij(c, bcpy, h, n, i, j, 1);
		    f1 = func(c, data);
		    if (na(f1)) {
			fprintf(stderr, "numerical_hessian: 2nd derivatives (%d,%d): "
				"objective function gave NA\n", i, j);
			*err = E_NAN;
			goto bailout;
		    }
		    hess_b_adjust_ij(c, bcpy, h, n, i, j, -1);
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
	    k = ijton(i, j, n);
	    V->val[k] = -D[u++];
	}
    }

    *err = gretl_invert_packed_symmetric_matrix(V);
    if (!*err) {
	vcv = gretl_matrix_steal_data(V);
    } else {
	fprintf(stderr, "numerical hessian: failed to invert V\n");
	gretl_packed_matrix_print(V, "V");
    }

    gretl_matrix_free(V);

 bailout:

    if (*err != 0 && *err != E_ALLOC) {
	gretl_errmsg_set(_("Failed to compute numerical Hessian"));
    }

    free(wspace);
    free(bcpy);

    return vcv;
}

/* default numerical calculation of gradient in context of BFGS */

int BFGS_numeric_gradient (double *b, double *g, int n,
			   BFGS_CRIT_FUNC func, void *data)
{
    double bi0, f1, f2;
    gretlopt opt = OPT_NONE;
    int i, err = 0;

    if (opt == OPT_R) {
	/* Richardson */
	double df[RSTEPS];
	double eps = 1.0e-4;
	double d = 0.0001;
	double v = 2.0;
	double h, p4m;
	int r = RSTEPS;
	int k, m;

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
		    err = 1;
		    goto bailout;
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
    } else {
	/* simple gradient calculation */
	const double h = 1.0e-8;

	for (i=0; i<n; i++) {
	    bi0 = b[i];
	    b[i] = bi0 - h;
	    f1 = func(b, data);
	    b[i] = bi0 + h;
	    f2 = func(b, data);
	    b[i] = bi0;
	    if (na(f1) || na(f2)) {
		err = 1;
		goto bailout;
	    }
	    g[i] = (f2 - f1) / (2.0 * h);
#if BFGS_DEBUG > 1
	    fprintf(stderr, "g[%d] = (%.16g - %.16g) / (2.0 * %g) = %g\n",
		    i, f2, f1, h, g[i]);
#endif
	}
    }

 bailout:

#if BFGS_DEBUG
    fprintf(stderr, "BFGS_numeric_gradient: returning %d\n", err);
#endif

    return err;
}

#define stepfrac	0.2
#define acctol		1.0e-7 /* alt: 0.0001 or 1.0e-7 (?) */
#define reltest		10.0

#define coeff_unchanged(a,b) (reltest + a == reltest + b)

static void reverse_gradient (double *g, int n)
{
    int i;

    for (i=0; i<n; i++) {
	g[i] = -g[i];
    }
}

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
				  double *reltol, gretlopt opt,
				  PRN *prn)
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
	*maxit = 500;
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
}

int BFGS_orig (double *b, int n, int maxit, double reltol,
	       int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	       int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	       gretlopt opt, PRN *prn)
{
    int crit_ok, done;
    double *wspace = NULL;
    double **H = NULL;
    double *g, *t, *X, *c;
    int verbose = (opt & OPT_V);
    int fcount, gcount, ndelta = 0;
    double fmax, f, f0, sumgrad, d = 0.0;
    int i, j, ilast, iter;
    double s, steplen = 0.0;
    double D1, D2;
    int err = 0;

    BFGS_get_user_values(b, n, &maxit, &reltol, opt, prn);

    if (gradfunc == NULL) {
	gradfunc = BFGS_numeric_gradient;
    }

    wspace = malloc(4 * n * sizeof *wspace);
    H = triangular_array_new(n);

    if (wspace == NULL || H == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    g = wspace;
    t = g + n;
    X = t + n;
    c = X + n;

    f = cfunc(b, data);

    if (na(f)) {
	gretl_errmsg_set(_("BFGS: initial value of objective function is not finite"));
	err = E_DATA;
	goto bailout;
    }

#if BFGS_DEBUG
    fprintf(stderr, "*** BFGS: first evaluation of f = %g\n", f);
#endif

    f0 = fmax = f;
    iter = ilast = fcount = gcount = 1;
    gradfunc(b, g, n, cfunc, data);
    reverse_gradient(g, n);

#if BFGS_DEBUG
    fprintf(stderr, "initial gradient:\n");
    for (i=0; i<n; i++) {
	fprintf(stderr, " g[%d] = %g\n", i, g[i]);
    }
#endif

    if (maxit == 0) {
	goto skipcalc;
    }

    do {
	if (verbose) {
	    reverse_gradient(g, n);
	    print_iter_info(iter, f, crittype, n, b, g, steplen, prn);
	    reverse_gradient(g, n);
	}

	if (ilast == gcount) {
	    /* (re-)start: initialize curvature matrix */
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

	sumgrad = 0.0;
	for (i=0; i<n; i++) {
	    s = 0.0;
	    for (j=0; j<=i; j++) {
		s -= H[i][j] * g[j];
	    }
	    for (j=i+1; j<n; j++) {
		s -= H[j][i] * g[j];
	    }
	    t[i] = s;
	    sumgrad += s * g[i];
	}

#if BFGS_DEBUG
	fprintf(stderr, "\niter %d: sumgrad = %g\n", iter, sumgrad);
#endif

	if (sumgrad < 0.0) {	
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
		fmax = f;
		gradfunc(b, g, n, cfunc, data);
		reverse_gradient(g, n);
		gcount++;
		iter++;
		D1 = 0.0;
		for (i=0; i<n; i++) {
		    t[i] = steplen * t[i];
		    c[i] = g[i] - c[i];
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
	reverse_gradient(g, n);
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
    double dsave[29];
    int isave[44];
    int lsave[4];
    int iter, ibak = 0;
    int err = 0;

    *fncount = *grcount = 0;    

    BFGS_get_user_values(b, n, &maxit, &reltol, opt, prn);

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
		f = -f;
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
	    fprintf(stderr, "%s\n", task);
	    break;
	}

	if (opt & OPT_V) {
	    if (iter != ibak) {
		double steplen = (iter == 1)? NADBL : dsave[13];

		reverse_gradient(g, n);
		print_iter_info(iter, -f, crittype, n, b, g, steplen, prn);
		reverse_gradient(g, n);
	    }
	    ibak = iter;
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
	      gretlopt opt, PRN *prn)
{
    if ((opt & OPT_L) || libset_get_bool(USE_LBFGS)) {
	return LBFGS_max(b, n, maxit, reltol,
			 fncount, grcount, cfunc, 
			 crittype, gradfunc, data, 
			 opt, prn);
    } else {
	return BFGS_orig(b, n, maxit, reltol,
			 fncount, grcount, cfunc, 
			 crittype, gradfunc, data, 
			 opt, prn);
    }
}

/* user-level access to BFGS and fdjac */

typedef struct umax_ umax;

struct umax_ {
    GretlType gentype;   /* GRETL_TYPE_DOUBLE or GRETL_TYPE_MATRIX */
    gretl_matrix *b;     /* parameter vector */
    int ncoeff;          /* number of coefficients */
    GENERATOR *g;        /* for generating scalar or matrix result */
    double x_out;        /* generated double value */
    gretl_matrix *m_out; /* generated matrix value */
    double ***Z;         /* pointer to data array */
    DATAINFO *dinfo;     /* dataset info */
    PRN *prn;            /* optional printing struct */
};

static umax *umax_new (GretlType t)
{
    umax *u = malloc(sizeof *u);

    if (u != NULL) {
	u->gentype = t;
	u->b = NULL;
	u->ncoeff = 0;
	u->g = NULL;
	u->x_out = NADBL;
	u->m_out = NULL;
	u->Z = NULL;
	u->dinfo = NULL;
	u->prn = NULL;
    }

    return u;
}

static void umax_destroy (umax *u)
{
    if (u->dinfo != NULL) {
	/* drop any private "$" variables created */
	dataset_drop_listed_variables(NULL, u->Z, u->dinfo, NULL, NULL);
    }

    destroy_genr(u->g);

    free(u);
}

static double user_get_criterion (const double *b, void *p)
{
    umax *u = (umax *) p;
    double x = NADBL;
    int i, t, err;

    for (i=0; i<u->ncoeff; i++) {
	u->b->val[i] = b[i];
    }

    err = execute_genr(u->g, u->Z, u->dinfo, OPT_S, u->prn); 

    if (err) {
	return NADBL;
    }

    t = genr_get_output_type(u->g);

    if (t == GRETL_TYPE_DOUBLE) {
	x = genr_get_output_scalar(u->g);
    } else if (t == GRETL_TYPE_MATRIX) {
	gretl_matrix *m = genr_get_output_matrix(u->g);

	if (gretl_matrix_is_scalar(m)) {
	    x = m->val[0];
	}
    } else {
	err = E_TYPES;
    }

    u->x_out = x;
    
    return x;
}

static int user_gen_setup (umax *u,
			   const char *fncall,
			   double ***pZ, 
			   DATAINFO *pdinfo)
{
    char formula[MAXLINE];
    GENERATOR *g;
    int err = 0;

    if (u->gentype == GRETL_TYPE_MATRIX) {
	sprintf(formula, "matrix $umax=%s", fncall);
    } else {
	sprintf(formula, "$umax=%s", fncall);
    }

    g = genr_compile(formula, pZ, pdinfo, OPT_P, &err);

    if (!err) {
	/* see if the formula actually works */
	err = execute_genr(g, pZ, pdinfo, OPT_S, u->prn);
    }

    if (!err) {
	u->g = g;
	u->Z = pZ;
	u->dinfo = pdinfo;
	u->m_out = genr_get_output_matrix(g);
    } else {
	destroy_genr(g);
    }

    return err;
}

double user_BFGS (gretl_matrix *b, const char *fncall,
		  double ***pZ, DATAINFO *pdinfo,
		  PRN *prn, int *err)
{
    umax *u;
    double ret = NADBL;
    gretlopt opt = OPT_NONE;
    int maxit = 500;
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

    *err = user_gen_setup(u, fncall, pZ, pdinfo);
    if (*err) {
	return NADBL;
    }

    tol = libset_get_double(BFGS_TOLER);

    if (libset_get_bool(MAX_VERBOSE)) {
	opt = OPT_V;
	u->prn = prn;
    }

    *err = BFGS_max(b->val, u->ncoeff, maxit, tol, &fcount, &gcount,
		    user_get_criterion, C_OTHER, NULL, u, opt, prn);

    if (fcount > 0) {
	pprintf(prn, _("Function evaluations: %d\n"), fcount);
	pprintf(prn, _("Evaluations of gradient: %d\n"), gcount);
    }

    if (!*err) {
	ret = u->x_out;
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

    err = execute_genr(u->g, u->Z, u->dinfo, OPT_S, u->prn); 
    if (err) {
	fprintf(stderr, "execute_genr: err = %d\n", err); 
    }

    if (err) {
	*iflag = -1;
	return 0;
    }

    v = genr_get_output_matrix(u->g);

#if JAC_DEBUG
    gretl_matrix_print(v, "matrix from u->g");
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

    *err = user_gen_setup(u, fncall, pZ, pdinfo);
    if (*err) {
	fprintf(stderr, "fdjac: error %d from user_gen_setup\n", *err);
	goto bailout;
    }

    if (u->m_out == NULL) {
	fprintf(stderr, "fdjac: u.m_out is NULL\n");
	*err = E_DATA; /* FIXME */
	goto bailout;
    }

    m = gretl_vector_get_length(u->m_out);
    if (m == 0) {
	*err = E_DATA;
	goto bailout;
    }
    
    *err = fdjac_allocate(m, n, &J, &wa, &fvec);
    if (*err) {
	goto bailout;
    }

    for (i=0; i<m; i++) {
	fvec[i] = u->m_out->val[i];
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
