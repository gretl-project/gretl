/* 
    The functions here are derived from the Fortran in VSGARCMX.FOR,
    from the Fiorentini, Calzolari and Panattoni GARCH implementation
    (Journal of Applied Econometrics, 1996), the original notice for
    which is reproduced infra.  The relevant code was extracted from
    its original Monte Carlo context, translated to C using f2c, and
    then edited extensively to turn it into more idiomatic C, and
    to convert from fixed-size arrays to dynamic memory allocation.

    I have also modified the comments at certain places in the code,
    where reference is made to the equations in the FCP paper.  The
    comments in the Fortran apparently pertained to a draft of the
    paper; I have updated them relative to the paper as published in
    JAE, 1996, pp. 399-417.

    Allin Cottrell, Wake Forest University, March 2004.
*/

/* Gabriele FIORENTINI, Giorgio CALZOLARI, Lorenzo PANATTONI
   Journal of APPLIED ECONOMETRICS, 1996 

   mixed gradient algorithm

   garch(p,q) estimates of a linear equation 

   SEE BOLLERSLEV, JOE 31(1986), 307-327. 
*/

#include "libgretl.h"
#include "libset.h"
#include "fcp.h"

#define FDEBUG 0

#define ABNUM   3   /* max number of GARCH lags */

#define SMALL_HT  1.0e-7 
#define S1MIN     1.0e-10
#define SUMGRMAX  1.0e-4

typedef struct fcpinfo_ fcpinfo;

struct fcpinfo_ {
    int nc;
    int t1, t2;
    int T;
    int nx;
    int p, q;
    int npar;
    double scale;

    const double *y;
    const double **X;
    double *e;
    double *e2;
    double *h;

    double *theta;
    double *grad;
    double *parpre;
    double **dhdp;
    double **g;

    gretl_matrix *vch;
};

static int fcp_allocate (fcpinfo *f)
{
    f->theta = malloc(f->npar * sizeof *f->theta);
    f->grad = malloc(f->npar * sizeof *f->grad);
    f->parpre = malloc(f->npar * sizeof *f->parpre);

    if (f->theta == NULL || f->grad == NULL || f->parpre) {
	return E_ALLOC;
    }

    f->dhdp = doubles_array_new(f->npar, f->T);
    if (f->dhdp == NULL) {
	return E_ALLOC;
    }

    f->g = doubles_array_new(f->nc, f->T);
    if (f->g == NULL) {
	return E_ALLOC;
    }

    f->vch = gretl_zero_matrix_new(f->npar, f->npar);
    if (f->vch == NULL) {
	return E_ALLOC;
    }  

    return 0;
}

static void fcpinfo_destroy (fcpinfo *f)
{
    free(f->theta);
    free(f->grad);
    free(f->parpre);
    doubles_array_free(f->dhdp, f->npar);
    doubles_array_free(f->g, f->nc); 
    gretl_matrix_free(f->vch);

    free(f);
}

static fcpinfo *fcpinfo_new (int nc, int q, int p, int t1, int t2, int T,
			     const double *y, const double **X, int nx,
			     double *e, double *e2, double *h,
			     double scale)
{
    fcpinfo *f = malloc(sizeof *f);

    if (f == NULL) {
	return NULL;
    }

    f->theta = NULL;
    f->grad = NULL;
    f->parpre = NULL;
    f->dhdp = NULL;
    f->g = NULL;
    f->vch = NULL;

    f->nc = nc;
    f->t1 = t1;
    f->t2 = t2;
    f->T = T,
    f->nx = nx;
    f->p = p;
    f->q = q;
    f->y = y;
    f->X = X;
    f->e = e;
    f->e2 = e2;
    f->h = h;   
    f->scale = scale;

    f->npar = nc + 1 + q + p;

    if (fcp_allocate(f)) {
	fcpinfo_destroy(f);
	f = NULL;
    }

    return f;
}

static void garch_iter_info (fcpinfo *f, int iter, double ll,
			     int hess, PRN *prn)
{
    double x;
    int i;

    pprintf(prn, "\n*** %s %d%s\n", ("iteration"), iter + 1,
	    (hess)? _(" (using Hessian)") : " (using Information Matrix)");


    pputs(prn, "Parameters:\n");

    for (i=0; i<f->npar; i++) {
	if (i && i % 5 == 0) {
	    pputc(prn, '\n');
	}
	x = f->theta[i];
	if (i < f->nc) {
	    x *= f->scale;
	} else if (i == f->nc) {
	    x *= f->scale * f->scale;
	}
	pprintf(prn, "%#12.5g ", x);
    }

    pprintf(prn, "\nll = %f\n", ll);
}

static double 
get_yhat (const double **X, int n, int t, const double *a)
{
    int i;
    double yhat = a[0];

    for (i=0; i<n; i++) {
	yhat += a[i+1] * X[i][t];
    }

    return yhat;
} 

static void copy_coeff (double *targ, const double *src, int n)
{
    int i;

    for (i=0; i<n; i++) {
	targ[i] = src[i];
    }
} 

static int g_ols (fcpinfo *f, double *c, double *amax, double *b)
{
    gretl_matrix *vc = NULL;
    double *aux = NULL;
    int i, j, k, t, err;
    double x, deltc, relinc = 0.5;
    double oldc, yh;

    vc = gretl_zero_matrix_new(f->nc, f->nc);
    if (vc == NULL) {
	return E_ALLOC;
    }

    aux = malloc(f->nc * sizeof *aux);
    if (aux == NULL) {
	gretl_matrix_free(vc);
	return E_ALLOC;
    }

    copy_coeff(b, c, f->nc);

    for (t=f->t1; t<=f->t2; t++) {
	amax[t] = get_yhat(f->X, f->nx, t, b);
    }

    for (i=0; i<f->nc; i++) {
	aux[i] = 0.0;
    }

    for (t=f->t1; t<=f->t2; t++) {
	for (k=0; k<f->nc; k++) {
	    oldc = c[k];
	    deltc = relinc;
	    if (oldc != 0.0) {
		deltc = oldc * relinc;
	    }
	    c[k] = oldc + deltc;
	    copy_coeff(b, c, f->nc);
	    yh = get_yhat(f->X, f->nx, t, b);
	    deltc = c[k] - oldc;
	    c[k] = oldc;
	    f->g[k][t] = (yh - amax[t]) / deltc;
	}
	copy_coeff(b, c, f->nc);

	/* cumulates all the w'z into diagonal blocks of vc, 
	   and w'y into elements of aux */
	for (i=0; i<f->nc; i++) {
	    aux[i] += f->g[i][t] * f->y[t];
	    for (j=0; j<f->nc; j++) {
		x = gretl_matrix_get(vc, i, j);
		x += f->g[i][t] * f->g[j][t];
		gretl_matrix_set(vc, i, j, x);
	    }
	}
    }

    err = gretl_invert_general_matrix(vc);

    if (err) {
	fputs("OLS: matrix is singular, initial coefficients are unchanged\n",
	      stderr);
    } else {
	/* compute coefficients */
	for (i=0; i<f->nc; i++) {
	    c[i] = 0.0;
	    for (j=0; j<f->nc; j++) {
		c[i] += gretl_matrix_get(vc, i, j) * aux[j];
	    }
	}
	copy_coeff(b, c, f->nc);
    } 

    gretl_matrix_free(vc);
    free(aux);

    return 0;
} 

/* Compute the GARCH log-likelihood.  Parameters are passed in f->theta;
   e, e2 and ht are be computed here (e2 holds squared residuals).
*/

static double garch_ll (fcpinfo *f, double *c)
{
    int t1 = f->t1;
    int t2 = f->t2;
    int p = f->p;
    int q = f->q;
    int nc = f->nc;
    int i, t, lag;
    int n = t2 - t1 + 1;
    double uncvar, ll;

    const double *alpha = f->theta + nc + 1;
    const double *beta = f->theta + nc + q + 1;

    for (i=0; i<nc; i++) {
	c[i] = f->theta[i];
    }

#if FDEBUG
    fprintf(stderr, "garch_ll: a0 = %.9g\n", f->theta[nc]);
    for (i=0; i<q; i++) {
	fprintf(stderr, " alpha[%d] = %.9g\n", i, alpha[i]);
    }
    for (i=0; i<p; i++) {
	fprintf(stderr, " beta[%d] = %.9g\n", i, beta[i]);
    }
#endif

    /* Compute residuals, squared residuals, and unconditional
       variance over the real estimation period 
    */
    uncvar = 0.0;
    for (t = t1; t <= t2; t++) {
	f->e[t] = f->y[t] - get_yhat(f->X, f->nx, t, c);
	f->e2[t] = f->e[t] * f->e[t];
	uncvar += f->e2[t];
    }
    uncvar /= n;

#if FDEBUG
    fprintf(stderr, "uncvar = %.9g (T=%d)\n", uncvar, n);
#endif

    /* 
       We use sample unconditional variance as the starting value 
       (at time 0, -1, -2, etc.) for the squared residuals and ht;
       we use 0 as the starting value for residuals.
    */

    lag = (p > q)? p : q;

    for (t = t1-lag; t < t1; ++t) { 
	f->e[t] = 0.0;
	f->e2[t] = f->h[t] = uncvar;
    }

    for (t=t1; t<=t2; t++) {
	f->h[t] = f->theta[nc];
	for (i=1; i<=q; i++) {
	    f->h[t] += f->e2[t-i] * alpha[i-1];
	}
	for (i=1; i<=p; i++) {
	    f->h[t] += f->h[t-i] * beta[i-1];
	}
	/* arbitrary */
	if (f->h[t] <= 0.0) {
	    f->h[t] = SMALL_HT;
	}
    }

    ll = 0.0;

    for (t=t1; t<=t2; t++) {
	double hts = f->h[t] * f->scale * f->scale;

	ll -= 0.5 * log(hts) + 0.5 * f->e2[t] / f->h[t] + LN_SQRT_2_PI;
    }

    return ll;
} 

/*
     Check that the values of the parameters of the conditional
     variance, ht, are in the set of the admissible values.  If a0 is
     less or equal than zero it is set to SMALL_HT.  If alpha and beta
     are less than zero they are set to zero; if the sum of alpha
     and beta is > 1.0, then alpha and beta are normalized (divided by
     sum).
 */

static void check_ht (double *b, int n)
{
    double sum = 0.0;
    int i;

    if (b[0] <= 0.0) {
	b[0] = SMALL_HT;
    }

    for (i=1; i<n; i++) {
	if (b[i] < 0.0) {
	    b[i] = 0.0;
	}
	sum += b[i];
    }

    if (sum > 1.0) {
	for (i=1; i<n; i++) {
	    b[i] /= sum;
	}
    }
} 

static void update_theta (fcpinfo *f, const double *gg, 
			  const double *step, double d)
{
    int i;

    for (i=0; i<f->npar; i++) {
	f->theta[i] = gg[i] + step[i] * d;
    }

    check_ht(f->theta + f->nc, f->p + f->q + 1);
}

/* combined setup for OP matrix, Information Matrix and Hessian */

static int vcv_setup (fcpinfo *f, double *c, int *count, 
		      double *zt, double ***H, gretl_matrix *V,
		      int code)
{
    int i, j, k, t, n, lag, nvpar;
    double x, *asum2;

    int t1 = f->t1;
    int t2 = f->t1;
    int p = f->p;
    int q = f->q;
    int nc = f->nc;

    double **dhdp = f->dhdp;
    double **g = f->dhdp;
    double *e = f->e;
    double *e2 = f->e2;
    double *h = f->h;

    const double *alpha = f->theta + nc + 1;
    const double *beta = f->theta + nc + q + 1;

#if FDEBUG
    fprintf(stderr, "vcv_setup: a0 = %.9g\n", a0);
    for (i=0; i<q; i++) {
	fprintf(stderr, " alpha[%d] = %.9g\n", i, alpha[i]);
    }
    for (i=0; i<p; i++) {
	fprintf(stderr, " beta[%d] = %.9g\n", i, beta[i]);
    }
#endif

    /* some useful abbreviations */
    nvpar = 1 + q + p;
    lag = (p > q)? p : q;
    n = t2 - t1 + 1;

#if FDEBUG
    fprintf(stderr, "make vcv: lag=%d, nc=%d, npar=%d\n",
	    lag, nc, npar);
#endif

    asum2 = malloc(nc * sizeof *asum2);
    if (asum2 == NULL) {
	return 1;
    }

    if (count != NULL) {
	++(*count);
    }

    for (i=0; i<nc; i++) {
	c[i] = f->theta[i];
    }

    /* Begin computation of dhtdp wrt the variance parameters; we
       start computing derivatives of starting values; for ht starting
       values are obtained from the unconditional variance of the
       residuals.
     */

    for (k=1; k<=p; k++) {
	for (i=0; i<nvpar; i++) {
	    dhdp[nc+i][t1-k] = 0.0;
	    if (H != NULL) { /* hessian only */
		for (j=0; j<nvpar; j++) {
		    H[nc+i][nc+j][k] = 0.0;
		}
	    }
	}
    }

    for (t=t1; t<=t2; t++) {
	/* fill in zt at time t (see p. 401) */
	zt[0] = 1.0;
	for (i=1; i<=q; i++) {
	    zt[i] = e2[t-i];
	}
	for (i=1; i<=p; i++) {
	    zt[q+i] = h[t-i];
	}

	/* Fill in dhtdp at time t, part relative to variance parameters
	   (eq. 7, p. 402) 
	*/
	for (i=0; i<nvpar; i++) {
	    dhdp[nc+i][t] = zt[i];
	    for (j=1; j<=p; j++) {
		dhdp[nc+i][t] += dhdp[nc+i][t-j] * beta[j-1];
	    }
	}
    }

    /* Build matrix dhtdp, block for regression coefficients
       (eq. 13). We use 0 as starting value (time 0, -1, etc.) for the
       derivatives of ht wrt the coefficients; we also use 0 as
       starting values for the residuals.
    */

    /* building blocks for pre-sample terms */
    for (i=0; i<nc; i++) {
	asum2[i] = 0.0;
	for (t=t1; t<=t2; t++) {
	    asum2[i] -= e[t] * 2.0 * g[i][t];
	}
	asum2[i] /= n; 
    }   

    /* pre-sample range */
    for (t=t1-lag; t<t1; t++) {
	for (i=0; i<nc; i++) {
	    dhdp[i][t] = asum2[i];
	}
    }

    /* actual sample range */
    for (t=t1; t<=t2; t++) {
	for (i=0; i<nc; i++) {
	    dhdp[i][t] = 0.0;
	    for (j=1; j<=q; j++) {
		if (t - q < t1) {
		    dhdp[i][t] += alpha[j-1] * asum2[i];
		} else {
		    dhdp[i][t] -= alpha[j-1] * 2.0 * g[i][t-j] * e[t-j];
		}
	    }
	    for (j=1; j<=p; j++) {
		dhdp[i][t] += dhdp[i][t-j] * beta[j-1];
	    }
	}
    }

    /* Initialize gradient and vcv */
    for (i=0; i<f->npar; i++) {
	f->grad[i] = 0.0;
    }
    gretl_matrix_zero(V);

    for (t=t1; t<=t2; t++) {
	double r_h = e[t] / h[t];
	double r2_h = e[t] * r_h;
	double aa, bb;

	/* 
	   First part, relative to regression coefficients (eq. 10, p. 402) 
 	*/
	for (i=0; i<nc; i++) {
	    aa = r_h * g[i][t] + .5 / h[t] * dhdp[i][t] * (r2_h - 1.0);
	    f->grad[i] += aa;
	    if (code == VCV_OP) {
		for (j=0; j<=i; j++) {
		    bb = r_h * g[j][t] + .5 / h[t] * dhdp[j][t] * (r2_h - 1.0);
		    x = gretl_matrix_get(V, i, j);
		    x += aa * bb;
		    gretl_matrix_set(V, i, j, x);
		    gretl_matrix_set(V, j, i, x);
		}
		for (j=0; j<nvpar; j++) {
		    int ncj = nc + j;

		    x = gretl_matrix_get(V, i, ncj);
		    x += aa * 0.5 / h[t] * dhdp[ncj][t] * (r2_h - 1.0);
		    gretl_matrix_set(V, i, ncj, x);
		    gretl_matrix_set(V, ncj, i, x);
		}
	    }
	}

	/* 
	   Second part, relative to variance parameters (eq. 6, p. 401) 
	*/
	for (i=0; i<nvpar; i++) {
	    int nci = nc + i;

	    aa = .5 / h[t] * dhdp[nci][t] * (r2_h - 1.0);
	    f->grad[nci] += aa;
	    if (code == VCV_OP) {
		for (j=0; j<=i; j++) {
		    int ncj = nc + j;

		    x = gretl_matrix_get(V, nci, ncj);
		    x += aa * 0.5 / h[t] * dhdp[ncj][t] * (r2_h - 1.0);
		    gretl_matrix_set(V, nci, ncj, x);
		    gretl_matrix_set(V, ncj, nci, x);
		}
	    }
	}
    }

    if (code == VCV_IM) {
	for (t=t1; t<=t2; t++) {
	    double ht2 = h[t] * h[t];

	    /* Part relative to the coefficients (eq. 30, p. 406).
	       Since we take the expected value, only the first two terms
	       remain.
	    */
	    for (i=0; i<nc; i++) {
		for (j=0; j<nc; j++) {
		    x = gretl_matrix_get(V, i, j);
		    x -= g[i][t] * g[j][t] / h[t] + .5 * dhdp[i][t] * dhdp[j][t] / ht2;
		    gretl_matrix_set(V, i, j, x);
		}
	    }

	    /* Part relative to the variance parameters (eq. 29, p. 406). 
	       Since we take the expected value, only the second term
	       remains.
	    */
	    for (i=nc; i<f->npar; i++) {
		for (j=nc; j<f->npar; j++) {
		    x = gretl_matrix_get(V, i, j);
		    x -= .5 * dhdp[i][t] * dhdp[j][t] / ht2;
		    gretl_matrix_set(V, i, j, x);
		}
	    }
	}
    }

    if (code != VCV_HESSIAN) {
	free(asum2);
	return 0;
    }

    /* Initial values in dhdpdp (here, "H") are 2/t x'x, and zero for
       off-diagonal (mixed) blocks. */

    for (k=0; k<lag; k++) {
	for (i=0; i<nc; i++) {
	    for (j=0; j<nc; j++) {
		H[i][j][k+1] = 0.; 
	    }
	}
	for (t=t1; t<=t2; t++) {
	    for (i=0; i<nc; i++) {
		for (j=0; j<nc; j++) {
		    H[i][j][k+1] += 2.0 * g[i][t] * g[j][t] / n;
		}
	    }
	}
	for (i=0; i<nc; i++) {
	    for (j=0; j<nvpar; j++) {
		/* mod. by AC: zero _all_ mixed entries */
		H[i][nc+j][k+1] = H[nc+j][i][k+1] = 0.0; 
	    }
	}
    }

    /* Now we fill out the full Hessian */

    for (t=t1; t<=t2; ++t) {
	double r_h = e[t] / h[t];
	double r2_h = r_h * e[t];
	double r2_h3 = r2_h / (h[t] * h[t]);
	double u_h2 = 1.0 / (h[t] * h[t]);

	for (i=0; i<f->npar; i++) {
	    for (j=0; j<f->npar; j++) {
		H[i][j][0] = 0.0; 
	    }
	}

	if (lag <= 0) {
	    goto L90;
	}

	for (k=1; k<=q; k++) {
	    for (i=0; i<nc; i++) {
		for (j=0; j<nc; j++) {
		    if (t - q < t1) {
			H[i][j][0] += H[i][j][q] * alpha[k-1];
		    } else {
			H[i][j][0] += 2.0 *
			    g[i][t-k] * g[j][t-k] * alpha[k-1];
		    }
		}
	    }
	}

	for (k=1; k<=p; k++) {
	    for (i=0; i<nc; i++) {
		for (j=0; j<nc; j++) {
		    H[i][j][0] += H[i][j][k] * beta[k-1];
		}
	    }
	}

	for (i=0; i<nc; i++) {
	    for (k=1; k<=q; k++) {
		if (t - q < t1) {
		    H[i][nc+k][0] += asum2[i];
		} else {
		    H[i][nc+k][0] -= 2.0 * g[i][t-k] * e[t-k];
		}
	    }
	    for (k=1; k<=p; k++) {
		H[i][nc+q+k][0] += dhdp[i][t-k];
	    }
	}

	for (k=1; k<=p; k++) { 
	    for (i=0; i<nc; i++) {
		for (j=0; j<nvpar; j++) {
		    H[i][nc+j][0] += H[i][nc+j][k] * beta[k-1];
		}
	    }
	}
    L90:

	/* Part relative to the coefficients (eq. 15, p. 403). 
	   Since we take the expected value, only the first two terms
	   remain.
 	*/

	for (i=0; i<nc; i++) {
	    for (j=0; j<nc; j++) {
		x = gretl_matrix_get(V, i, j);
		x = x - g[i][t] * g[j][t] / h[t] 
		    - .5 * r2_h3 * dhdp[i][t] * dhdp[j][t] 
		    - (r_h * g[j][t] * dhdp[i][t]) / h[t] 
		    - (r_h * g[i][t] * dhdp[j][t]) / h[t] 
		    + 0.5 * (r2_h - 1.0) * 
		    (H[i][j][0] / h[t] - dhdp[i][t] 
		     * dhdp[j][t] / (h[t] * h[t]));
		gretl_matrix_set(V, i, j, x);
	    }
	}

	/* Part relative to the variance parameters (eq. 14, p. 403). 
	   Since we take the expected value, only the second term
	   remains.
 	*/

	if (p > 0) {
	    for (i=0; i<nvpar; i++) {
		for (j=1; j<=p; j++) {
		    H[nc+i][nc+q+j][0] += dhdp[nc+i][t-j];
		}
	    }
	    for (i=1; i<=p; i++) {
		for (j=0; j<nvpar; j++) {
		    H[nc+q+i][nc+j][0] += dhdp[nc+j][t-i];
		}
	    }
	    for (k=1; k<=p; k++) {
		for (i=0; i<nvpar; i++) {
		    for (j=0; j<nvpar; j++) { 
			H[nc+i][nc+j][0] += H[nc+i][nc+j][k] * beta[k-1];
#if FDEBUG
			if ((t==t1 || t==t2)) {
			    fprintf(stderr, "Set H[%d][%d][0] = %.9g\n",
				    nc+i, nc+j, H[nc+i][nc+j][0]);
			    fprintf(stderr, "incr = %.9g * %.9g\n",
				    H[nc+i][nc+j][k], beta[k-1]);
			}
#endif
		    }
		}
	    }
	}

	for (i=nc; i<f->npar; i++) {
	    for (j=nc; j<f->npar; j++) {
		x = gretl_matrix_get(V, i, j);
		x = x + .5 * u_h2 * dhdp[i][t] * dhdp[j][t] 
		    - r2_h3 * dhdp[i][t] * dhdp[j][t] 
		    + .5 * (r2_h - 1.0) / h[t] * H[i][j][0];
		gretl_matrix_set(V, i, j, x);
#if FDEBUG
		if ((t==t1 || t==t2)) {
		    fprintf(stderr, "Set V(%d,%d) = %.9g (using Hval %.9g)\n",
			    i, j, x, H[i][j][0]);
		}
#endif
	    }
	}

	/* top-right mixed part (eq. 17, p. 403) */
	for (i=0; i<nc; i++) {
	    for (j=0; j<nvpar; j++) {
		x = gretl_matrix_get(V, i, nc+j);
		x = x - g[i][t] * r_h * dhdp[nc+j][t] / h[t] 
		    - .5 * (r2_h - 1.0) * dhdp[nc+j][t] * dhdp[i][t] / (h[t] * h[t]) 
		    + .5 * (r2_h - 1.0) * H[i][nc+j][0] / h[t] 
		    - .5 * r2_h * u_h2 * dhdp[i][t] * dhdp[nc+j][t];
		gretl_matrix_set(V, i, nc+j, x);
		/* and bottom left too */
		gretl_matrix_set(V, nc+j, i, x);
	    }
	}

	/* before quitting time t, tidy up dhdpdp */
	for (k=0; k<lag; k++) { 
	    for (i=0; i<f->npar; i++) {
		for (j=0; j<f->npar; j++) {
#if FDEBUG
		    if (t < 5) {
			fprintf(stderr, "t=%d: setting H(%d,%d,%d) to H(%d,%d,%d)=%g\n",
				t,i,j,lag-k,i,j,lag-k-1, H[i][j][lag-k-1]);
		    }
#endif
		    H[i][j][lag-k] = H[i][j][lag-k-1];
		}
	    }
	}
    }

    free(asum2);

    return 0;
} /* vcv_setup */


/* calculate the step for the new coefficients */

static double step_calc (fcpinfo *f,
			 double *gg, double *step,
			 const gretl_matrix *V, const double *grad,
			 double *ds, double *ps2)
{
    double s1 = 0.0, s2 = 0.0;
    double vij, stre;
    int i, j;

    for (i=0; i<f->npar; i++) {
	s1 += f->theta[i] * f->theta[i];
	gg[i] = f->theta[i];
	step[i] = 0.0;
	for (j=0; j<f->npar; j++) {
	    vij = gretl_matrix_get(V, i, j);
	    step[i] -= vij * grad[j];
	}
	s2 += step[i] * step[i];
    }

    if (s1 == 0.0) {
	s1 = S1MIN;
    }

    stre = sqrt(s2 / s1);
    s2 = sqrt(s2);

    for (i=0; i<f->npar; i++) {
	step[i] /= s2;
    }

    *ds = *ps2 = s2;

    return stre;
}

/* Block-diagonal information matrix.  Parameters are passed in the
   array "c"
*/

static int 
garch_info_matrix (fcpinfo *f, double *c, double toler, int *count, 
		   gretl_matrix *V, double *b, double *zt)
{
    static double ll1 = 0.0, fs = 0.0;

    int p = f->p;
    int q = f->q;
    int npar = f->npar;

    double d, d0, d1, d2, ll2, d3;
    double ll3, d12, d31, d23, dd;
    double di, ff, dm, ds;
    int nexp, ncall, nvpar;
    double a1s, a2s, a3s;
    double d12s, dmin, dmax, d23s, d31s;
    double s2, stre, bigd;
    double *gg = NULL, *step = NULL;
    int err;

    nvpar = q + p + 1;

    /* calculate information matrix */
    vcv_setup(f, c, count, zt, NULL, V, VCV_IM);

    err = gretl_invert_general_matrix(V);
    if (err) {
	fprintf(stderr, "garch_info_matrix: matrix inversion failed\n");
	return 1;
    }

    if (count == NULL) {
	/* just calculating vcv at convergence */
	goto vcv_exit;
    }

    gg = malloc(npar * sizeof *gg);
    if (gg == NULL) {
	return 1;
    }

    step = malloc(npar * sizeof *step);
    if (step == NULL) {
	free(gg);
	return 1;
    } 

    /* here we start iteration */

    /* calculate the step for the new coefficients */
    stre = step_calc(f, gg, step, V, f->grad, &ds, &s2);

#if FDEBUG
    fprintf(stderr, "s2 = %.9g\n", s2);
#endif

    if (stre <= toler) {
	goto L496;
    }

    ncall = 0;
    nexp = *count / 5;
    if (nexp > 5) {
	nexp = 5;
    }

    d0 = s2 / pow(2.0, nexp);
    dmin = d0 * .001;
    dmax = d0 * 4.0;

    if (*count == 1) {
	ll1 = -garch_ll(f, c);
#if FDEBUG > 1
	fprintf(stderr, "count=1, ll1=%.9g\n", ll1);
#endif
    }

    update_theta(f, gg, step, d0);

#if FDEBUG
    for (i=0; i<f->npar; i++) {
	fprintf(stderr, "theta[%d] in matinf(1) = %.9g\n", i, f->theta[i]);
    }
#endif 

    ll2 = -garch_ll(f, c);
#if FDEBUG
    fprintf(stderr, "ll2=%.9g, ll1=%.9g\n", ll2, ll1);
#endif    

    if (ll2 > ll1) {
	d1 = -d0;
	d2 = 0.0;
	d3 = d0;
	ll3 = ll2;
	ll2 = ll1;

	update_theta(f, gg, step, d1);

	ll1 = -garch_ll(f, c);
    } else {
	d1 = 0.0;
	d2 = d0;
	d3 = d0 + d0;

	update_theta(f, gg, step, d3);

	ll3 = -garch_ll(f, c);
    }

 L325:
    d23 = d2 - d3;
    d31 = d3 - d1;
    d12 = d1 - d2;
    di = d23 * ll1 + d31 * ll2 + d12 * ll3;
    bigd = di * -2.0 / (d23 * d31 * d12);

    if (bigd > 0.0) {
	goto L400;
    }

    if (ll3 <= ll1) {
	goto L341;
    }

 L329:
    d3 = d2;
    ll3 = ll2;
    d2 = d1;
    ll2 = ll1;
    d1 -= dmax;

    update_theta(f, gg, step, d1);

#if FDEBUG
    for (i=0; i<f->npar; i++) {
	fprintf(stderr, "theta[%d] in matinf(4) = %.9g\n", i, f->theta[i]);
    }
#endif 

    ll1 = -garch_ll(f, c);

    if (++ncall > 100) {
	goto endloop;
    }
    goto L325;

 L341:
    d1 = d2;
    ll1 = ll2;
    d2 = d3;
    ll2 = ll3;
    d3 += dmax;

    update_theta(f, gg, step, d3);

#if FDEBUG
    for (i=0; i<npar; i++) {
	fprintf(stderr, "theta[%d] in matinf(5) = %.9g\n", i, theta[i]);
    }
#endif

    ll3 = -garch_ll(f, c);

    if (++ncall > 100) {
	goto endloop;
    }
    goto L325;

 L400:
    d23s = d23 * (d2 + d3);
    d31s = d31 * (d3 + d1);
    d12s = d12 * (d1 + d2);
    ds = (d23s * ll1 + d31s * ll2 + d12s * ll3) * .5 / di;

    update_theta(f, gg, step, ds);

#if FDEBUG
    for (i=0; i<npar; i++) {
	fprintf(stderr, "theta[%d] in matinf(6) = %.9g\n", i, theta[i]);
    }
#endif

    fs = -garch_ll(f, c);

    if (++ncall > 100) {
	goto endloop;
    }

    a1s = (d = d1 - ds, fabs(d));
    a2s = (d = d2 - ds, fabs(d));
    a3s = (d = d3 - ds, fabs(d));

    dm = a1s;
    if (a3s < dm) {
	dm = a3s;
    }

    if (dmax < dm) {
	if (ds < d1 - dmax) {
	    goto L329;
	}
	if (ds > d3 + dmax) {
	    goto L341;
	}
    }

    if (a1s < dmin || a2s < dmin || a3s < dmin) {
	goto endloop;
    }

    if (ll1 < ll2 || ll1 < ll3) {
	if (ll2 < ll3 || ll2 < ll1) {
	    d3 = ds;
	    ll3 = fs;
	} else {
	    d2 = ds;
	    ll2 = fs;
	}
    } else {
	d1 = ds;
	ll1 = fs;
    }

 L459:
    if (d2 > d3) {
	dd = d2;
	ff = ll2;
	d2 = d3;
	ll2 = ll3;
	d3 = dd;
	ll3 = ff;
    }

    if (d1 <= d2) {
	goto L325;
    }

    dd = d1;
    ff = ll1;
    d1 = d2;
    ll1 = ll2;
    d2 = dd;
    ll2 = ff;
    goto L459;

 endloop:

    if (fs > ll1) {
	fs = ll1;
	ds = d1;
    }
    if (fs > ll2) {
	fs = ll2;
	ds = d2;
    }
    if (fs > ll3) {
	fs = ll3;
	ds = d3;
    }

 L496:
    update_theta(f, gg, step, ds);

#if FDEBUG
    for (i=0; i<npar; i++) {
	fprintf(stderr, "theta[%d] in matinf(7) = %.9g\n", i, theta[i]);
    }
#endif

    ll1 = fs;
    fs = -fs;

 vcv_exit:

    gretl_matrix_switch_sign(V);

    free(gg);
    free(step);

    return 0;
} /* garch_info_matrix */

static void free_dhdpdp (double ***H, int np)
{
    int i, j;

    for (i=0; i<np; i++) {
	for (j=0; j<np; j++) {
	    free(H[i][j]);
	}
	free(H[i]);
    }
    
    free(H);
}

static double ***allocate_dhdpdp (int np, int p, int q)
{
    double ***H;
    int i, j, lag = (p > q)? p : q;

    H = malloc(np * sizeof *H);
    if (H == NULL) {
	return NULL;
    }
    
    for (i=0; i<np; i++) {
	H[i] = malloc(np * sizeof **H);
	if (H[i] == NULL) {
	    return NULL;
	}
	for (j=0; j<np; j++) {
	    H[i][j] = malloc((lag + 1) * sizeof ***H);
	    if (H[i][j] == NULL) {
		return NULL;
	    }
	}
    }

    return H;
}

static int 
garch_full_hessian (fcpinfo *f, double *c, double toler, int *count, 
		    double *b, double *zt)
{
    static double ll1 = 0.0, fs = 0.0;

    int p = f->p;
    int q = f->q;
    int npar = f->npar;

    double d, d0, d1, d2, ll2, d3; 
    double ll3, d12, d31, d23, dd, di, ff; 
    int ncall, nexp, nvpar, lag;
    double a1s, a2s, a3s;
    double dm, ds, dmin, d12s, dmax, d23s, d31s;
    double s2, stre, bigd;
    double *gg = NULL, *step = NULL;
    double ***H = NULL;
    int err;

    gg = malloc(npar * sizeof *gg);
    if (gg == NULL) {
	return E_ALLOC;
    }

    step = malloc(npar * sizeof *step);
    if (step == NULL) {
	free(gg);
	return E_ALLOC;
    } 

    /* 3rd dimension of dhdpdp is max(p,q) + 1 */
    H = allocate_dhdpdp(npar, p, q); 
    if (H == NULL) {
	free(gg);
	free(step);
	return E_ALLOC;
    }

    nvpar = q + p + 1;
    lag = (p > q)? p : q;

    /* calculate the full Hessian */
    vcv_setup(f, c, count, zt, H, f->vch, VCV_HESSIAN);

    /* invert the Hessian */
    err = gretl_invert_symmetric_indef_matrix(f->vch);
    if (err) {
	fprintf(stderr, "garch_full_hessian: matrix inversion failed\n");
    }

    /* Start iteration here */

    /* calculate the step for the new coefficients */
    stre = step_calc(f, gg, step, f->vch, f->grad, &ds, &s2);

    if (stre <= toler) {
	goto L496;
    }

    ncall = 0;
    nexp = *count / 5;
    if (nexp > 5) {
	nexp = 5;
    }

    d0 = s2 / pow(2.0, nexp);
    dmin = d0 * .001;
    dmax = d0 * 4.0;

    if (*count == 1) {
	ll1 = -garch_ll(f, c);
#if FDEBUG
	fprintf(stderr, "hess: ll1 = %.9g\n", ll1);
#endif  
    }

    update_theta(f, gg, step, d0);

    ll2 = -garch_ll(f, c);
#if FDEBUG
    fprintf(stderr, "hess: ll2 = %.9g\n", ll2);
#endif   
 
    if (ll2 > ll1) {
	d1 = -d0;
	d2 = 0.0;
	d3 = d0;
	ll3 = ll2;
	ll2 = ll1;

	update_theta(f, gg, step, d1);

	ll1 = -garch_ll(f, c);
#if FDEBUG
	fprintf(stderr, "hess: ll1(2) = %.9g\n", ll1);
#endif    
    } else {
	d1 = 0.0;
	d2 = d0;
	d3 = d0 + d0;

	update_theta(f, gg, step, d3);

	ll3 = -garch_ll(f, c);
#if FDEBUG
	fprintf(stderr, "hess: ll3 = %.9g\n", ll3);
#endif    
    }

 L325:
    d23 = d2 - d3;
    d31 = d3 - d1;
    d12 = d1 - d2;
    di = d23 * ll1 + d31 * ll2 + d12 * ll3;
    bigd = di * -2.0 / (d23 * d31 * d12);
#if FDEBUG
    fprintf(stderr, "hess: bigd = %.9g\n", bigd);
#endif    

    if (bigd > 0.0) {
	goto L400;
    }

    if (ll3 <= ll1) {
	goto L341;
    }

 L329:
    d3 = d2;
    ll3 = ll2;
    d2 = d1;
    ll2 = ll1;
    d1 -= dmax;

    update_theta(f, gg, step, d1);

    ll1 = -garch_ll(f, c);
#if FDEBUG
    fprintf(stderr, "hess: ll1(3) = %.9g\n", ll1);
#endif    

    if (++ncall > 100) {
	goto endloop;
    }
    goto L325;

 L341:
    d1 = d2;
    ll1 = ll2;
    d2 = d3;
    ll2 = ll3;
    d3 += dmax;

    update_theta(f, gg, step, d3);

    ll3 = -garch_ll(f, c);
#if FDEBUG
    fprintf(stderr, "hess: ll3(2) = %.9g\n", ll3);
#endif    

    if (++ncall > 100) {
	goto endloop;
    }
    goto L325;

 L400:
    d23s = d23 * (d2 + d3);
    d31s = d31 * (d3 + d1);
    d12s = d12 * (d1 + d2);
    ds = (d23s * ll1 + d31s * ll2 + d12s * ll3) * .5 / di;

    update_theta(f, gg, step, ds);

    fs = -garch_ll(f, c);
#if FDEBUG
    fprintf(stderr, "hess: fs = %.9g, ds = %.9g\n", fs, ds);
#endif    

    if (++ncall > 100) {
	goto endloop;
    }

    a1s = (d = d1 - ds, fabs(d));
    a2s = (d = d2 - ds, fabs(d));
    a3s = (d = d3 - ds, fabs(d));

    dm = a1s;
    if (a3s < dm) {
	dm = a3s;
    }

    if (dmax < dm) {
	if (ds < d1 - dmax) {
	    goto L329;
	}
	if (ds > d3 + dmax) {
	    goto L341;
	}
    }

    if (a1s < dmin || a2s < dmin || a3s < dmin) {
	goto endloop;
    }

    if (ll1 < ll2 || ll1 < ll3) {
	if (ll2 < ll3 || ll2 < ll1) {
	    d3 = ds;
	    ll3 = fs;
	} else {
	    d2 = ds;
	    ll2 = fs;
	}
    } else {
	d1 = ds;
	ll1 = fs;
    }

 L459:
    if (d2 > d3) {
	dd = d2;
	ff = ll2;
	d2 = d3;
	ll2 = ll3;
	d3 = dd;
	ll3 = ff;
    }

    if (d1 <= d2) {
	goto L325;
    }
    dd = d1;
    ff = ll1;
    d1 = d2;
    ll1 = ll2;
    d2 = dd;
    ll2 = ff;
    goto L459;

 endloop:

    if (fs > ll1) {
	fs = ll1;
	ds = d1;
    }
    if (fs > ll2) {
	fs = ll2;
	ds = d2;
    }
    if (fs > ll3) {
	fs = ll3;
	ds = d3;
    }

 L496:
    update_theta(f, gg, step, ds);

    ll1 = fs;
    fs = -fs;

    gretl_matrix_switch_sign(f->vch);

    free(gg);
    free(step);
    free_dhdpdp(H, npar);

    return 0;
} /* garch_full_hessian */

static int 
make_garch_vcv (fcpinfo *f, double *c, double *b, 
		double *zt, const gretl_matrix *ihess,
		gretl_matrix *V, int vopt)
{
    gretl_matrix *OP = NULL;
    gretl_matrix *iinfo = NULL;
    int k = f->npar;
    int err = 0;

    /* OP and robust variants need OP matrix */
    if (vopt == VCV_OP || vopt == VCV_QML || vopt == VCV_BW) {
	OP = gretl_matrix_alloc(k, k);
	if (OP == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}

	vcv_setup(f, c, NULL, zt, NULL, OP, VCV_OP);

	if (vopt == VCV_OP) {
	    gretl_matrix_copy_values(V, OP);
	    err = gretl_invert_symmetric_matrix(V);
	}
    }

    /* IM and BW variants need the info matrix */
    if (vopt == VCV_IM || vopt == VCV_BW) {
	iinfo = gretl_matrix_alloc(k, k);
	if (iinfo == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}

	garch_info_matrix(f, c, 0.0, NULL, iinfo, b, zt);

	if (vopt == VCV_IM) {
	    gretl_matrix_copy_values(V, iinfo);
	} else {
	    /* Bollerslev-Wooldridge */
	    gretl_matrix_qform(iinfo, GRETL_MOD_NONE,
			       OP, V, GRETL_MOD_NONE);
	}
    } else if (vopt == VCV_QML) {
	gretl_matrix_qform(ihess, GRETL_MOD_NONE,
			   OP, V, GRETL_MOD_NONE);
    } else if (vopt == VCV_HESSIAN) {
	gretl_matrix_copy_values(V, ihess);
    }	

 bailout:

    gretl_matrix_free(OP);
    gretl_matrix_free(iinfo);

    return err;
}

/*
   Parameters to garch_estimate()

   t1:    beginning of sample in auxiliary database
   t2:    end of sample in auxiliary database
   nobs:  total number of observations in auxiliary database
   X:     data matrix for auxiliary database (regressors, not needed on
          output)
   nx:    number of columns of X
   coeff: vector of coefficient for the conditional mean, normally
          initialised by OLS on input (not needed on output)
   nc:    number of elements in coeff
   V:     covariance matrix of coeffs (all 0 on input)
   e:     vector of 0's on input, resids on output
   e2:    storage for squared residuals
   res:   vector of 0's on input, resids on output
   h:     null pointer on input, conditional variances on output
   y:     on input, vector with dep. var.; not needed on output
   amax:  vector; element 0 holds the garch intercept; 1 and 2 the
          arch & garch orders; from 3 onwards, the arch & garch 
          parameters
   b:     0 on input, holds vector of coefficient for the conditional
          mean on output
   scale: double used to scale dependent variable
   iters: int, 0 on input, holds number of iterations on output
   prn:   print handle for info on iterations and other diagnostic output

*/

int garch_estimate (int t1, int t2, int nobs, 
		    const double **X, int nx, double *coeff, int nc, 
		    gretl_matrix *V, double *e, double *e2, double *h,
		    const double *y, double *amax, double *b, 
		    double scale, int *iters, PRN *prn, int vopt)
{
    fcpinfo *f;
    int it1, it2, ittot;
    int count, conv = 0;
    int i, q, p, npar;

    double pdiff, tol1, tol2;
    double sumgra; 
    double ll, s1, s2;

    double *c = NULL;

    double zt[6];   /* 6 == max value of (1 + q + p) */

    int err = 0;

    q = (int) amax[1];
    p = (int) amax[2];
    npar = nc + 1 + q + p;

    c = malloc(nc * sizeof *c);
    if (c == NULL) {
	return E_ALLOC;
    }

    f = fcpinfo_new(nc, q, p, t1, t2, nobs, y, X, nx,
		    e, e2, h, scale);
    if (f == NULL) {
	free(c);
	return E_ALLOC;
    }

    for (i=0; i<q; i++) {
	f->theta[nc + i + 1] = amax[3 + i];
    }

    for (i=0; i<p; i++) {
	f->theta[nc + q + i + 1] = amax[3 + q + i];
    }

    f->theta[nc] = amax[0];

    tol1 = .05; /* tolerance before switching to Hessian */
    tol2 = 1e-8;

    for (i=0; i<nc; i++) {
	c[i] = f->theta[i] = coeff[i];
    }

#if FDEBUG
    fprintf(stderr, "before g_ols:\n");
    for (i=0; i<nc; i++) {
	fprintf(stderr, "c[%d]=%g, b[%d]=%g\n", i, c[i], i, b[i]);
    }
#endif

    /* this is only to calculate matrix of regressors (g) */
    g_ols(f, c, amax, b);

#if FDEBUG
    fprintf(stderr, "after g_ols:\n");
    for (i=0; i<nc; i++) {
	fprintf(stderr, "c[%d]=%g, b[%d]=%g\n", i, c[i], i, b[i]);
    }
#endif    

    /* iterative estimation */

    count = 0;

    for (it1=0; it1<100; it1++) {

	ll = garch_ll(f, c);

	garch_iter_info(f, it1, ll, 0, prn);

	/* store previous coefficients */
	for (i=0; i<npar; i++) {
	    f->parpre[i] = f->theta[i];
	}

#if FDEBUG
	fprintf(stderr, "*** Calling garch_info_matrix, round %d\n", it1);	    
#endif

	err = garch_info_matrix(f, c, tol1, &count, f->vch, b, zt);

	if (err) {
	    return E_NOCONV;
	}

	s1 = s2 = 0.0;
	for (i=0; i<npar; i++) {
	    s1 += f->parpre[i] * f->parpre[i];
	    pdiff = f->theta[i] - f->parpre[i];
	    s2 += pdiff * pdiff;
	}

	if (s1 == 0.0) {
	    s1 = S1MIN;
	}

	if (s2 / s1 <= tol1 * tol1) {
	    break;
	}
    }

#if FDEBUG
    fprintf(stderr, "\n\n*** Now going to Hessian ***\n\n");
#endif

    ittot = it1 + 1;

    for (it2=0; it2<100 && !conv; it2++) {

	/* compute residuals for covariance matrix */
	ll = garch_ll(f, c);

	/* store previous coefficients */
	for (i=0; i<npar; i++) {
	    f->parpre[i] = f->theta[i];
	}

#if FDEBUG
	fprintf(stderr, "*** Calling garch_full_hessian, round %d\n", it2);	    
#endif

	garch_full_hessian(f, c, tol2, &it2, b, zt);

	garch_iter_info(f, ittot++, ll, 1, prn);

	s1 = s2 = 0.0;

	for (i=0; i<npar; i++) {
	    s1 += f->parpre[i] * f->parpre[i];
	    pdiff = f->theta[i] - f->parpre[i];
	    s2 += pdiff * pdiff;
	}

	if (s1 == 0.0) {
	    s1 = S1MIN;
	}

	if (s2 / s1 > tol2 * tol2) {
	    /* not converged yet */
	    continue;
	}

	sumgra = 0.0;
	for (i=0; i<npar; i++) {
	    sumgra += f->grad[i] * f->grad[i];
	}

        if (sumgra >= SUMGRMAX) {
	    pprintf(prn, "\nParameters and gradients at iteration %d:\n\n", 
		    ittot);

	    for (i=0; i<f->npar; i++) {
		pprintf(prn, "%12.6f (%9.6f)\n", f->theta[i], f->grad[i]);
	    }

            pprintf(prn, "\nSum of squared gradients = %.9g (should be less " 
		    "than %g)\n", sumgra, SUMGRMAX);
            err = E_NOCONV;
            goto garch_exit;
	} else {
	    pprintf(prn, "\nFull Hessian convergence at iteration %d, "
		"tol = %.9g\n\n", ittot, tol2);
	    conv = 1;
	    *iters = ittot;
	    amax[0] = ll;
	}
    }

    if (!conv && !err) {
	err = E_NOCONV;
    }

    if (!err) {
	double sderr;

	/* build the desired VCV variant */
	err = make_garch_vcv(f, c, b, zt, f->vch, V, vopt); /* FIXME vch? */

	if (!err) {
	    /* transcribe coefficients and standard errors */
	    for (i=0; i<npar; i++) {
		double vii = gretl_matrix_get(V, i, i);

		sderr = (vii > 0.0)? sqrt(vii) : 0.0;
		amax[i+1] = f->theta[i];
		amax[i+1+npar] = sderr;
	    }
	}
    }

 garch_exit:

    fcpinfo_destroy(f);
    free(c);

    return err;
}

