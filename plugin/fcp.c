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

#include "f2c.h"
#include "clapack_double.h"

#undef FDEBUG

#define NLL    50   /* number of iterative results to store */
#define ABNUM   3   /* max number of GARCH lags */

#define SMALL_HT      1.0e-7 

int global_np;
double gscale;
int gncoeff;

#define vix(i,j) ((i) + global_np * (j))

/* private functions */

static int invert (double *g, int dim);

static  
double garch_ll (double *c, int nc, double *res2, 
		 double *res, double *yhat, 
		 const double *y, const double **X, 
		 int nexo, int t1, int t2, double *param, 
		 double *b, double *a0, 
		 double *alfa, double *beta, int q, 
		 int p, double *h);

static int vcv_setup (int t1, int t2, double *c, int nc, 
		      double *res2, double *res, int *count, 
		      const double **g, double *aux3, 
		      double *param, int nparam, 
		      double *a0, double *alfa, double *beta, 
		      int q, int p, double *h, 
		      double **dhdp, double *zt, 
		      double ***H, double *vco, int code);

static int 
garch_info_matrix (int t1, int t2, 
		   const double **X, int nexo, 
		   double *yhat, double *c, int nc, double *res2,
		   double *res, const double *y, double toler, 
		   int *count, double *vc5, const double **g,
		   double *aux3, double *param, int nparam,
		   double *b, double *a0, double *alfa,
		   double *beta, int q, int p, double *h,
		   double **dhdp, double *zt);

static int 
garch_full_hessian (int t1, int t2, 
		    const double **X, int nexo, 
		    double *yhat, double *c, int nc, double *res2,
		    double *res, const double *y, double toler, 
		    int *count, double *vc5, const double **g,
		    double *aux3, double *param, int nparam,
		    double *b, double *a0, double *alfa,
		    double *beta, int q, int p, double *h,
		    double **dhdp, double *zt);

static void free_2d_array (double **A, int k)
{
    int i;

    if (A == NULL) return;

    for (i=0; i<k; i++) free(A[i]);
    free(A);
}

static double **allocate_2d_array (int k, int T)
{
    double **A;
    int i, j;

    A = malloc(k * sizeof *A);
    if (A == NULL) return NULL;

    for (i=0; i<k; i++) {
	A[i] = malloc(T * sizeof **A);
	if (A[i] == NULL) {
	    for (j=0; j<i; j++) free(A[j]);
	    free(A);
	    return NULL;
	}
    }

    return A;
}

static int vs_allocate (double ***pdhdp, double ***pg, 
			double **pparam, double **paux3, double **psvc5,
			double **pc, double **paux, double **pvc5,
			double **pparpre, double ***ppartrc,
			int np, int nrc, int T, int nll)
{
    double *param = NULL, *aux3 = NULL, *svc5 = NULL;
    double *c = NULL, *aux = NULL, *vc5 = NULL, *parpre = NULL;
    double **D = NULL, **G = NULL, **P = NULL;

    param = malloc(np * sizeof *param);
    aux3 = malloc(np * sizeof *aux3);
    svc5 = malloc(np * sizeof *svc5);
    parpre = malloc(np * sizeof *parpre);
    if (param == NULL || aux3 == NULL || 
	svc5 == NULL || parpre == NULL) {
	goto bailout;
    }

    c = malloc(nrc * sizeof *c);
    aux = malloc(nrc * sizeof *aux);
    if (c == NULL || aux == NULL) {
	goto bailout;
    }

    vc5 = malloc(np * np * sizeof *vc5);
    if (vc5 == NULL) {
	goto bailout;
    }

    D = allocate_2d_array(np, T);
    if (D == NULL) {
	goto bailout;
    }

    G = allocate_2d_array(nrc, T);
    if (G == NULL) {
	goto bailout;
    }

    P = allocate_2d_array(np, nll);
    if (P == NULL) {
	goto bailout;
    }

    *pdhdp = D;
    *pg = G;
    *pparam = param;
    *paux3 = aux3;
    *psvc5 = svc5;
    *pc = c;
    *paux = aux;
    *pvc5 = vc5;
    *pparpre = parpre;
    *ppartrc = P;

    return 0;

 bailout:

    free(param);
    free(aux3);
    free(svc5);
    free(c);
    free(aux);
    free(vc5);
    free(parpre);
    free_2d_array(D, np);
    free_2d_array(G, nrc);
    free_2d_array(P, np);

    return 1;
}

static void vs_free (double **dhdp, int np, double **g, int nrc, 
		     double *param, double *aux3, double *svc5,
		     double *c, double *aux, double *vcv,
		     double *parpre, double **partrc)
{
    free_2d_array(dhdp, np);
    free_2d_array(g, nrc);
    free_2d_array(partrc, np);
    free(param);
    free(aux3);
    free(svc5);
    free(c);
    free(aux);
    free(vcv);
    free(parpre);
}

static void print_iter_info (int iter, double *theta, int m, double ll,
			     int hess, PRN *prn)
{
    int i;
    double x;

    pprintf(prn, "\n*** %s %d%s: theta, ll ***\n", ("iteration"), iter,
	    (hess)? _(" (using Hessian)") : "");
    for (i=0; i<m; i++) {
	if (i && i % 5 == 0) {
	    pputc(prn, '\n');
	}
	x = theta[i];
	if (i < gncoeff - 1) {
	    x *= gscale;
	} else if (i == gncoeff - 1) {
	    x *= gscale * gscale;
	}
	pprintf(prn, "%#12.5g ", x);
    }
    pprintf(prn, "\n    ll = %f\n", ll);
}

static double 
get_yhat (const double **X, int nx, int t, const double *a)
{
    int j;
    double yhat = a[0];

    for (j = 0; j < nx; j++) {
	yhat += a[j+1] * X[j][t];
    }

    return yhat;
} 

static void copy_coeff (const double *c, int nc, double *b)
{
    int i;

    for (i = 0; i < nc; ++i) {
	b[i] = c[i];
    }
} 

int ols_ (int t1, int t2, const double **X, int nx, 
	  double *c, int nc, const double *y, 
	  double *amax, double *aux, double *b, double **g)
{
    int i, j, k, t, err;
    double deltc, deriv, relinc = 0.5;
    double oldc, yy;
    double *vc;

    vc = malloc(nc * nc * sizeof *vc);
    if (vc == NULL) return 1;

    copy_coeff(c, nc, b);

    for (t = t1; t <= t2; ++t) {
	amax[t] = get_yhat(X, nx, t, b);
    }

    for (i = 0; i < nc; ++i) {
	aux[i] = 0.0;
	for (j = 0; j < nc; ++j) {
	    vc[i + j * nc] = 0.0;
	}
    }

    for (t = t1; t <= t2; ++t) {
	for (k = 0; k < nc; ++k) {
	    oldc = c[k];
	    deltc = relinc;
	    if (oldc != 0.0) {
		deltc = oldc * relinc;
	    }
	    c[k] = oldc + deltc;
	    copy_coeff(c, nc, b);
	    yy = get_yhat(X, nx, t, b);
	    deltc = c[k] - oldc;
	    deriv = (yy - amax[t]) / deltc;
	    c[k] = oldc;
	    g[k][t] = deriv;
	}
	copy_coeff(c, nc, b);

	/* cumulates all the w'z into diagonal blocks of vc 
	   and w'y into elements of aux */
	for (i = 0; i < nc; ++i) {
	    aux[i] += g[i][t] * y[t];
	    for (j = 0; j < nc; ++j) {
		vc[i + j * nc] += g[i][t] * g[j][t];
	    }
	}
    }

    err = invert(vc, nc);

    if (err == 0) {
	/* compute coefficients */
	for (i = 0; i < nc; ++i) {
	    c[i] = 0.0;
	}
	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nc; ++j) {
		c[i] += vc[i + j * nc] * aux[j];
	    }
	}
	copy_coeff(c, nc, b);
    } else {
	fputs("OLS: matrix is singular, initial coefficients are unchanged\n",
	      stderr);

	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nc; ++j) {
		vc[i + j * nc] = 0.0;
	    }
	}
    }

    free(vc);

    return 0;
} 

/* Calculate robust VCV.  If the first input (vci) is the inverse of
   the Hessian, the output is the QML (White) estimator; if vci is the
   inverse of the information matrix, the output is the
   Bollerslev-Wooldridge estimator. The second input matrix (vco) is
   the (uninverted) outer product (OP) matrix.
*/

static double *robust_vcv (const double *vci, const double *vco,
			   int nparam) 
{
    int i, j, k;
    int np2 = nparam * nparam;
    double *vv = NULL, *vcr = NULL;

    vv = malloc(np2 * sizeof *vv); /* temporary workspace */
    vcr = malloc(np2 * sizeof *vcr);

    if (vv == NULL || vcr == NULL) {
	free(vv);
	free(vcr);
	return NULL;     
    }

    /* multiply H^{-1} (or I^{-1}) into OP */
    for (i=0; i<nparam; i++) { 
	for (j=0; j<nparam; j++) {
	    vv[vix(i,j)] = 0.0;
	    for (k=0; k<nparam; k++) {
		vv[vix(i,j)] += vci[vix(i,k)] * vco[vix(k,j)];
	    }
	}
    }

    /* post-multiply by H^{-1} (or I^{-1}) */
    for (i=0; i<nparam; i++) { 
	for (j=0; j<nparam; j++) {
	    vcr[vix(i,j)] = 0.0;
	    for (k=0; k<nparam; k++) {
		vcr[vix(i,j)] += vv[vix(i,k)] * vci[vix(k,j)];
	    }
	}
    }

    free(vv);

    return vcr;
}

static int 
make_garch_vcv (int t1, int t2, 
		const double **X, int nx,
		double *yhat, const double *y,
		double *c, int nc, 
		double *res2, double *res, 
		const double **g, double *aux3,
		double *param, int nparam,
		double *b, double *a0, 
		double *alfa, double *beta, int q, int p, 
		double *h, double **dhdp, double *zt,
		const double *vch, double *vcv,
		int vopt)
{
    double *vco = NULL, *vcr = NULL, *vci = NULL;
    int np2 = nparam * nparam;
    int k, err = 0;

    /* OP and robust variants need OP matrix */
    if (vopt == VCV_OP || vopt == VCV_QML || vopt == VCV_BW) {
	vco = malloc(nparam * nparam * sizeof *vco);
	if (vco == NULL) {
	    err = 1;
	    goto bailout;
	}

	vcv_setup(t1, t2, c, nc, 
		  res2, res, 
		  NULL, g, aux3, 
		  param, nparam, 
		  a0, alfa, beta, q, p, h,
		  dhdp, zt, NULL, vco, 
		  VCV_OP);

	if (vopt == VCV_OP) {
	    err = invert(vco, nparam);
	    for (k=0; k<np2; k++) {
		vcv[k] = vco[k];
	    }
	}
    }

    /* IM and BW variants need the info matrix */
    if (vopt == VCV_IM || vopt == VCV_BW) {
	vci = malloc(nparam * nparam * sizeof *vci);
	if (vci == NULL) {
	    err = 1;
	    goto bailout;
	}

	garch_info_matrix(t1, t2, X, nx, yhat, c, 
			  nc, res2, res, y,
			  0.0, NULL, vci, g, 
			  aux3, param, nparam, b, 
			  a0, alfa, beta, q,
			  p, h, dhdp, zt);

	if (vopt == VCV_IM) {
	    for (k=0; k<np2; k++) {
		vcv[k] = vci[k];
	    }
	} else {
	    /* Bollerslev-Wooldridge */
            vcr = robust_vcv(vci, vco, nparam);
	    if (vcr == NULL) {
		err = 1;
		goto bailout;
	    }
	    for (k=0; k<np2; k++) {
		vcv[k] = vcr[k];
	    }
	}
    }

    else if (vopt == VCV_QML) {
	vcr = robust_vcv(vch, vco, nparam);
	if (vcr == NULL) {
	    err = 1;
	    goto bailout;
	}
	for (k=0; k<np2; k++) {
	    vcv[k] = vcr[k];
	}	
    }	
	
    else if (vopt == VCV_HESSIAN) {
	for (k=0; k<np2; k++) {
	    vcv[k] = vch[k];
	}
    }	

 bailout:

    free(vco);
    free(vci);
    free(vcr);

    return err;
}

int garch_estimate (int t1, int t2, int nobs, 
		    const double **X, int nx, double *yhat, 
		    double *coeff, int nc, double *vcv, 
		    double *res2, double *res, double *h,
		    const double *y, double *amax, double *b, 
		    double scale, int *iters, PRN *prn, int vopt)
{
    int i, j;

    int izo, nzo, nzo1;
    int count, count2;
    int q, p, nparam;

    double alfa[ABNUM], beta[ABNUM];
    double zt[6];   /* max value of (1 + q + p) */

    double pappo, toler1, toler2;
    double reldis, rellog, tollog, sumgra, totdis; 
    double ll, s_2, a0, s_1;

    double *param = NULL, *aux3 = NULL, *svc5 = NULL;
    double **dhdp = NULL, **g = NULL;
    double *c = NULL, *aux = NULL;
    double *parpre, **partrc;
    double *vc5 = NULL;

    int err = 0;

    q = (int) amax[1];
    p = (int) amax[2];

    /* "export" scale as file-scope global */
    gscale = scale;

    /* number of parameters of unconcentrated likelihood */
    nparam = nc + 1 + q + p;
    global_np = nparam;
    gncoeff = nc + 1;

    if (vs_allocate(&dhdp, &g, &param, &aux3, &svc5, 
		    &c, &aux, &vc5, &parpre, &partrc,
		    nparam, nc, nobs, NLL)) {
	pprintf(prn, "Out of memory\n");
	return E_ALLOC;
    }

    for (i = 0; i < q; ++i) {
	param[nc + i + 1] = amax[3 + i];
    }

    for (i = 0; i < p; ++i) {
	param[nc + q + i + 1] = amax[3 + q + i];
    }

    param[nc] = amax[0];

    toler1 = .05;
    toler2 = 1e-8;
    tollog = log10(toler2);

    for (i = 0; i < nc; ++i) {
	c[i] = coeff[i];
    }

    for (i = 0; i < nparam; ++i) {
	svc5[i] = 0.0;
    }

    for (i = 0; i < nc; ++i) {
	param[i] = coeff[i];
    }

    /* this is only to calculate matrix of regressors (g) */
    ols_(t1, t2, X, nx, c, nc, y, amax, aux, b, g);

#ifdef FFDEBUG
    for (i = 0; i < nc; i++) {
	fprintf(stderr, "after ols g[%d] = %.9g\n", i, g[i]);
    } 
#endif    

    /* iterative estimation */

    count = 0;
    nzo = 0;

    for (izo = 1; izo <= 100; ++izo) {

	ll = garch_ll(c, nc, res2, res, yhat, y, 
		      X, nx, t1, t2, param, b, &a0, alfa, 
		      beta, q, p, h);

	print_iter_info(izo, param, nparam, ll, 0, prn);

	if (++nzo > NLL) {
	    nzo = NLL;
	}
	
	/* store previous coefficients */
	for (i = 0; i < nparam; ++i) {
	    parpre[i] = param[i];
	    partrc[i][nzo-1] = param[i];
	}

#ifdef FDEBUG
	fprintf(stderr, "*** Calling garch_info_matrix, count=%d\n", count);	    
#endif

	err = garch_info_matrix(t1, t2, X, nx, yhat, c, 
				nc, res2, res, y,
				toler1, &count, vc5, (const double **) g, 
				aux3, param, nparam, b, 
				&a0, alfa, beta, q,
				p, h, dhdp, zt);

	if (err) {
	    return E_NOCONV;
	}

	/* if relative euclidean distance is used as converg. */
	s_1 = s_2 = 0.0;
	for (i = 0; i < nparam; ++i) {
	    s_1 += parpre[i] * parpre[i];
	    pappo = param[i] - parpre[i];
	    s_2 += pappo * pappo;
	}

	if (s_1 == 0.0) {
	    s_1 = 1e-10;
	}

	if (s_2 / s_1 <= toler1 * toler1) {
	    break;
	}
    }

#ifdef FDEBUG
    fprintf(stderr, "\n\n*** Now going to Hessian ***\n\n");
#endif

    count2 = 0;
    nzo1 = 0;

    for (izo = 1; izo <= 100; ++izo) {

	/* compute residuals for covariance matrix */
	ll = garch_ll(c, nc, res2, res, yhat, y, 
		      X, nx, t1, t2, param, b, &a0, alfa, 
		      beta, q, p, h);

	print_iter_info(nzo + 1, param, nparam, ll, 1, prn);

	if (++nzo > NLL) {
	    nzo = NLL;
	}

	/* store previous coefficients */
	for (i = 0; i < nparam; ++i) {
	    parpre[i] = param[i];
	    partrc[i][nzo-1] = param[i]; 
	}

#ifdef FDEBUG
	fprintf(stderr, "*** Calling garch_full_hessian, count2=%d\n", count2);	    
#endif

	garch_full_hessian(t1, t2, X, nx, yhat, c, nc, 
			   res2, res, y, toler2, &count2, vc5, 
			   (const double **) g, 
			   aux3, param, nparam, b, &a0, alfa, beta, 
			   q, p, h, dhdp, zt);

	/* if relative Euclidean distance is used as converg. */
	s_1 = 0.0;
	s_2 = 0.0;

	for (i = 0; i < nparam; ++i) {
	    s_1 += parpre[i] * parpre[i];
	    pappo = param[i] - parpre[i];
	    s_2 += pappo * pappo;
	}

	if (s_1 == 0.0) s_1 = 1e-10;

	if (s_2 / s_1 > toler2 * toler2) {
	    goto L6765;
	}

	sumgra = 0.0;
	for (i = 0; i < nparam; ++i) {
	    sumgra += aux3[i] * aux3[i];
	}
	if (sumgra >= 1.0e-4) {
	    fprintf(stderr, "Sum of gradients = %.9g\n", (double) sumgra);
	    err = E_NOCONV;
	    goto L999;
	}

	pprintf(prn, "\nFull Hessian convergence at iteration %d, "
		"tol = %.9g\n\n", nzo, toler2);

	*iters = nzo;
	amax[0] = ll;

	tollog = 0.0;
	totdis = 0.0;
	for (j = 0; j < nparam; ++j) {
	    /* Compute 2nd power */
	    double d = param[j] - partrc[j][0];

	    totdis += d * d;
	}
	totdis = sqrt(totdis);

	for (i = 0; i < nzo; ++i) {
	    /* if euclidean distance or distance in log-likel. is used */
	    s_2 = 0.0;
	    for (j = 0; j < nparam; ++j) {
		/* Compute 2nd power */
		double d = param[j] - partrc[j][i];

		s_2 += d * d;
	    }
	    s_2 = sqrt(s_2);

	    reldis = (totdis != 0.0)? s_2 / totdis : 0.0;
	    rellog = (reldis != 0.0)? reldis : tollog;
	}

	if (++nzo1 > NLL) {
	    goto L8767;
	}
    L8767:
	goto L6766;
    L6765:
	;
    }

    err = E_NOCONV;
    goto L999;

 L6766:
    if (!err) {
	double sderr;

	/* build the desired VCV variant */
	err = make_garch_vcv(t1, t2, X, nx,
			     yhat, y, c, nc, 
			     res2, res, 
			     (const double **) g, aux3,
			     param, nparam,
			     b, &a0, alfa, beta, q, p, 
			     h, dhdp, zt,
			     (const double *) vc5, vcv,
			     vopt);

	if (!err) {
	    /* transcribe coefficients and standard errors */
	    for (i = 0; i < nparam; ++i) {
		if (vcv[vix(i,i)] > 0.0) {
		    sderr = sqrt(vcv[vix(i,i)]);
		} else {
		    sderr = 0.0;
		}
		amax[i+1] = param[i];
		amax[i+1+nparam] = sderr;
	    }
	}
    }

 L999:
    vs_free(dhdp, nparam, g, nc, param, aux3, svc5, c, aux, vc5,
	    parpre, partrc);

    return err;
}

/* Compute the log-likelihood.  Parameters are passed in the
   param(nparam) vector; a0, alfa and beta get extracted from param in
   garch_ll; res, res2 and ht must be computed inside garch_ll (res2
   holds squared residuals).
*/

static double 
garch_ll (double *c, int nc, double *res2, 
	  double *res, double *yhat, 
	  const double *y, const double **X, int nx,
	  int t1, int t2, double *param, double *b,  
	  double *a0, double *alfa, double *beta, int q,
	  int p, double *h)
{
    int i, t, lag;
    int n = t2 - t1 + 1;
    double uncvar, ll;

    for (i = 0; i < nc; ++i) {
	c[i] = param[i];
    }

    *a0 = param[nc];
#ifdef FDEBUG
    fprintf(stderr, "garch_ll: a0 = %.9g\n", *a0);
#endif

    for (i = 0; i < q; ++i) {
	alfa[i] = param[nc + 1 + i];
#ifdef FDEBUG
	fprintf(stderr, " alfa[%d] = %.9g\n", i, alfa[i]);
#endif	
    }

    for (i = 0; i < p; ++i) {
	beta[i] = param[nc + q + 1 + i];
#ifdef FDEBUG
	fprintf(stderr, " beta[%d] = %.9g\n", i, beta[i]);
#endif	
    }

    /* Compute residuals, squared residuals, and unconditional
       variance over the real estimation period */
    copy_coeff(c, nc, b);
    uncvar = 0.0;
    for (t = t1; t <= t2; ++t) {
	yhat[t] = get_yhat(X, nx, t, b);
	res[t] = y[t] - yhat[t];
	res2[t] = res[t] * res[t];
	uncvar += res2[t];
    }
    uncvar /= n;

#ifdef FDEBUG
    fprintf(stderr, "uncvar = %.9g (T=%d)\n", uncvar, n);
#endif

    /* 
       We use sample unconditional variance as the starting value 
       (at time 0, -1, -2, etc.) for the squared residuals and ht;
       we use 0 as the starting value for residuals.
    */

    lag = (p > q)? p : q;
    for (t = t1-lag; t < t1; ++t) { 
	res[t] = 0.0;
	res2[t] = h[t] = uncvar;
    }

    for (t = t1; t <= t2; ++t) {
	h[t] = *a0;
	for (i = 1; i <= q; ++i) {
	    h[t] += res2[t-i] * alfa[i-1];
	}
	for (i = 1; i <= p; ++i) {
	    h[t] += h[t-i] * beta[i-1];
	}
	/* arbitrary */
	if (h[t] <= 0.0) h[t] = SMALL_HT;
    }

    ll = 0.0;
    for (t = t1; t <= t2; ++t) {
	double hts = h[t] * gscale * gscale;

	ll -= 0.5 * log(hts) + 0.5 * res2[t] / h[t] + LN_SQRT_2_PI;
    }

    return ll;
} 

static void check_ht (double *param, int np)
{
    /*
     Check that the values of the parameters of the conditional
     variance, ht, are in the set of the admissible values.  If a0 is
     less or equal than zero it is set to SMALL_HT.  If alfa and beta
     are less than zero they are set to zero; also if the sum of alfa
     and beta is > 1.0, then alfa and beta are normalized (divided by
     sum).
    */

    int i;
    double sum = 0.;

    if (param[0] <= 0.) {
	param[0] = SMALL_HT;
    }

    for (i = 1; i < np; i++) {
	if (param[i] < 0.) {
	    param[i] = 0.0;
	}
	sum += param[i];
    }

    if (sum > 1.0) {
	for (i = 1; i < np; ++i) {
	    param[i] /= sum;
	}
    }
} 

/* combined setup for OP matrix, Information Matrix and Hessian */

static int vcv_setup (int t1, int t2, double *c, int nc, 
		      double *res2, double *res, int *count, 
		      const double **g, double *aux3, 
		      double *param, int nparam, 
		      double *a0, double *alfa, double *beta, 
		      int q, int p, double *h, 
		      double **dhdp, double *zt, 
		      double ***H, double *vcv, int code)
{
    int i, j, k, t, n, lag, nvparm;
    double *asum2;

    /* some useful abbreviations */
    nvparm = 1 + q + p;
    lag = (p > q)? p : q;
    n = t2 - t1 + 1;

#ifdef FDEBUG
    fprintf(stderr, "make vcv: lag=%d, nc=%d, nparam=%d\n",
	    lag, nc, nparam);
#endif

    asum2 = malloc(nc * sizeof *asum2);
    if (asum2 == NULL) {
	return 1;
    }

    if (count != NULL) {
	++(*count);
    }

    for (i = 0; i < nc; ++i) {
	c[i] = param[i];
    }

    *a0 = param[nc];
#ifdef FDEBUG
    fprintf(stderr, "a0=%.9g\n", *a0);
#endif

    for (i = 0; i < q; ++i) {
	alfa[i] = param[nc + 1 + i];
#ifdef FDEBUG
	fprintf(stderr, "alfa[%d]=%.9g\n", i, alfa[i]);
#endif
    }

    for (i = 0; i < p; ++i) {
	beta[i] = param[nc + 1 + q + i];
#ifdef FDEBUG
	fprintf(stderr, "beta[%d]=%.9g\n", i, beta[i]);
#endif
    }

    /* Begin computation of dhtdp wrt the variance parameters; we
       start computing derivatives of starting values; for ht starting
       values are obtained from the unconditional variance of the
       residuals.
     */

    for (k = 1; k <= p; k++) {
	for (i = 0; i < nvparm; ++i) {
	    dhdp[nc+i][t1-k] = 0.0;
	    if (H != NULL) { /* hessian only */
		for (j = 0; j < nvparm; ++j) {
		    H[nc+i][nc+j][k] = 0.0;
		}
	    }
	}
    }

    for (t = t1; t <= t2; ++t) {
	/* fill in zt at time t (see p. 401) */
	zt[0] = 1.0;
	for (i = 1; i <= q; ++i) {
	    zt[i] = res2[t-i];
	}
	for (i = 1; i <= p; ++i) {
	    zt[q+i] = h[t-i];
	}

	/* Fill in dhtdp at time t, part relative to variance parameters
	   (eq. 7, p. 402) 
	*/
	for (i = 0; i < nvparm; ++i) {
	    dhdp[nc+i][t] = zt[i];
	    for (j = 1; j <= p; ++j) {
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
    for (i = 0; i < nc; ++i) {
	asum2[i] = 0.0;
	for (t = t1; t <= t2; t++) {
	    asum2[i] -= res[t] * 2.0 * g[i][t];
	}
	asum2[i] /= n; 
    }   

    /* pre-sample range */
    for (t = t1-lag; t < t1; ++t) {
	for (i = 0; i < nc; ++i) {
	    dhdp[i][t] = asum2[i];
	}
    }

    /* actual sample range */
    for (t = t1; t <= t2; ++t) {
	for (i = 0; i < nc; i++) {
	    dhdp[i][t] = 0.0;
	    for (j = 1; j <= q; ++j) {
		if (t - q < t1) {
		    dhdp[i][t] += alfa[j-1] * asum2[i];
		} else {
		    dhdp[i][t] -= alfa[j-1] * 2.0 * g[i][t-j] * res[t-j];
		}
	    }
	    for (j = 1; j <= p; ++j) {
		dhdp[i][t] += dhdp[i][t-j] * beta[j-1];
	    }
	}
    }

    /* Initialize gradient, aux3, and vcv */
    for (i = 0; i < nparam; i++) {
	aux3[i] = 0.0;
	for (j = 0; j < nparam; j++) {
	    vcv[vix(i,j)] = 0.0;
	}
    }

    for (t = t1; t <= t2; ++t) {
	double r_h = res[t] / h[t];
	double r2_h = res[t] * r_h;
	double aa, bb;

	/* 
	   First part, relative to regression coefficients (eq. 10, p. 402) 
 	*/
	for (i = 0; i < nc; ++i) {
	    aa = r_h * g[i][t] + .5 / h[t] * dhdp[i][t] * (r2_h - 1.0);
	    aux3[i] += aa;
	    if (code == VCV_OP) {
		for (j=0; j<=i; j++) {
		    bb = r_h * g[j][t] + .5 / h[t] * dhdp[j][t] * (r2_h - 1.0);
		    vcv[vix(i,j)] += aa * bb;
		    vcv[vix(j,i)] = vcv[vix(i,j)];
		}
		for (j=0; j<nvparm; j++) {
		    int ncj = nc + j;

		    vcv[vix(i,ncj)] += aa * 0.5 / h[t] * dhdp[ncj][t] * (r2_h - 1.0);
		    vcv[vix(ncj,i)] = vcv[vix(i,ncj)];
		}
	    }
	}

	/* 
	   Second part, relative to variance parameters (eq. 6, p. 401) 
	*/
	for (i = 0; i < nvparm; ++i) {
	    int nci = nc + i;

	    aa = .5 / h[t] * dhdp[nci][t] * (r2_h - 1.0);
	    aux3[nci] += aa;
	    if (code == VCV_OP) {
		for (j=0; j<=i; j++) {
		    int ncj = nc + j;

		    vcv[vix(nci,ncj)] += aa * 0.5 / h[t] * dhdp[ncj][t] * (r2_h - 1.0);
		    vcv[vix(ncj,nci)] = vcv[vix(nci,ncj)];
		}
	    }
	}
    }

    if (code == VCV_IM) {
	for (t = t1; t <= t2; ++t) {
	    double ht2 = h[t] * h[t];

	    /* Part relative to the coefficients (eq. 30, p. 406).
	       Since we take the expected value, only the first two terms
	       remain.
	    */
	    for (i = 0; i < nc; ++i) {
		for (j = 0; j < nc; ++j) {
		    vcv[vix(i,j)] = vcv[vix(i,j)] 
			- g[i][t] * g[j][t] / h[t] 
			- .5 * dhdp[i][t] * dhdp[j][t] / ht2;
		}
	    }

	    /* Part relative to the variance parameters (eq. 29, p. 406). 
	       Since we take the expected value, only the second term
	       remains.
	    */
	    for (i = nc; i < nparam; ++i) {
		for (j = nc; j < nparam; ++j) {
		    vcv[vix(i,j)] -= .5 * dhdp[i][t] * dhdp[j][t] / ht2;
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
    for (k = 0; k < lag; ++k) {
	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nc; ++j) {
		H[i][j][k+1] = 0.; 
	    }
	}
	for (t = t1; t <= t2; ++t) {
	    for (i = 0; i < nc; ++i) {
		for (j = 0; j < nc; ++j) {
		    H[i][j][k+1] += 2.0 * g[i][t] * g[j][t] / n;
		}
	    }
	}
	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nvparm; ++j) {
		/* mod. by AC: zero _all_ mixed entries */
		H[i][nc+j][k+1] = H[nc+j][i][k+1] = 0.0; 
	    }
	}
    }

    /* Now we fill out the full Hessian */

    for (t = t1; t <= t2; ++t) {
	double r_h = res[t] / h[t];
	double r2_h = r_h * res[t];
	double r2_h3 = r2_h / (h[t] * h[t]);
	double u_h2 = 1.0 / (h[t] * h[t]);

	for (i = 0; i < nparam; ++i) {
	    for (j = 0; j < nparam; ++j) {
		H[i][j][0] = 0.0; 
	    }
	}

	if (lag <= 0) {
	    goto L90;
	}

	for (k = 1; k <= q; ++k) {
	    for (i = 0; i < nc; ++i) {
		for (j = 0; j < nc; ++j) {
		    if (t - q < t1) {
			H[i][j][0] += H[i][j][q] * alfa[k-1];
		    } else {
			H[i][j][0] += 2.0 *
			    g[i][t-k] * g[j][t-k] * alfa[k-1];
		    }
		}
	    }
	}

	for (k = 1; k <= p; ++k) {
	    for (i = 0; i < nc; ++i) {
		for (j = 0; j < nc; ++j) {
		    H[i][j][0] += H[i][j][k] * beta[k-1];
		}
	    }
	}

	for (i = 0; i < nc; ++i) {
	    for (k = 1; k <= q; ++k) {
		if (t - q < t1) {
		    H[i][nc+k][0] += asum2[i];
		} else {
		    H[i][nc+k][0] -= 2.0 * g[i][t-k] * res[t-k];
		}
	    }
	    for (k = 1; k <= p; ++k) {
		H[i][nc+q+k][0] += dhdp[i][t-k];
	    }
	}

	for (k = 1; k <= p; ++k) { 
	    for (i = 0; i < nc; ++i) {
		for (j = 0; j < nvparm; ++j) {
		    H[i][nc+j][0] += H[i][nc+j][k] * beta[k-1];
		}
	    }
	}
    L90:

	/* Part relative to the coefficients (eq. 15, p. 403). 
	   Since we take the expected value, only the first two terms
	   remain.
 	*/
	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nc; ++j) {
		vcv[vix(i,j)] = vcv[vix(i,j)] 
		    - g[i][t] * g[j][t] / h[t] 
		    - .5 * r2_h3 * dhdp[i][t] * dhdp[j][t] 
		    - (r_h * g[j][t] * dhdp[i][t]) / h[t] 
		    - (r_h * g[i][t] * dhdp[j][t]) / h[t] 
		    + 0.5 * (r2_h - 1.0) * 
		    (H[i][j][0] / h[t] - dhdp[i][t] 
		     * dhdp[j][t] / (h[t] * h[t]));
	    }
	}

	/* Part relative to the variance parameters (eq. 14, p. 403). 
	   Since we take the expected value, only the second term
	   remains.
 	*/
	if (p > 0) {
	    for (i = 0; i < nvparm; ++i) {
		for (j = 1; j <= p; ++j) {
		    H[nc+i][nc+q+j][0] += dhdp[nc+i][t-j];
		}
	    }
	    for (i = 1; i <= p; ++i) {
		for (j = 0; j < nvparm; ++j) {
		    H[nc+q+i][nc+j][0] += dhdp[nc+j][t-i];
		}
	    }
	    for (k = 1; k <= p; ++k) {
		for (i = 0; i < nvparm; ++i) {
		    for (j = 0; j < nvparm; ++j) { 
			H[nc+i][nc+j][0] += H[nc+i][nc+j][k] * beta[k-1];
#ifdef FDEBUG
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

	for (i = nc; i < nparam; ++i) {
	    for (j = nc; j < nparam; ++j) {
		vcv[vix(i,j)] = vcv[vix(i,j)]
		    + .5 * u_h2 * dhdp[i][t] * dhdp[j][t] 
		    - r2_h3 * dhdp[i][t] * dhdp[j][t] 
		    + .5 * (r2_h - 1.0) / h[t] * H[i][j][0];
#ifdef FDEBUG
		    if ((t==t1 || t==t2)) {
			fprintf(stderr, "Set vcv[%d] = %.9g (using Hval %.9g)\n",
				vix(i,j), vcv[vix(i,j)], H[i][j][0]);
		    }
#endif
	    }
	}

	/* top-right mixed part (eq. 17, p. 403) */
	for (i = 0; i < nc; ++i) {
	    for (j = 0; j < nvparm; ++j) {
		vcv[vix(i,nc+j)] = vcv[vix(i,nc+j)]
		    - g[i][t] * r_h * dhdp[nc+j][t] / h[t] 
		    - .5 * (r2_h - 1.0) * dhdp[nc+j][t] * dhdp[i][t] / (h[t] * h[t]) 
		    + .5 * (r2_h - 1.0) * H[i][nc+j][0] / h[t] 
		    - .5 * r2_h * u_h2 * dhdp[i][t] * dhdp[nc+j][t];
		/* and bottom left too */
		vcv[vix(nc+j, i)] = vcv[vix(i, nc+j)];
	    }
	}

	/* before quitting time t, tidy up dhdpdp */
	for (k = 0; k < lag; ++k) { 
	    for (i = 0; i < nparam; ++i) {
		for (j = 0; j < nparam; ++j) {
#ifdef FDEBUG
		    if (t<5) fprintf(stderr, "t=%d: setting H(%d,%d,%d) to H(%d,%d,%d)=%g\n",
				     t,i,j,lag-k,i,j,lag-k-1, H[i][j][lag-k-1]);
#endif
		    H[i][j][lag-k] = H[i][j][lag-k-1];
		}
	    }
	}
    }

    free(asum2);

    return 0;
} /* vcv_setup */


/* block-diagonal information matrix */

/* parameters are passed within the vector param
   c, alfa and beta are fetched from param at the start
*/

static int 
garch_info_matrix (int t1, int t2, 
		   const double **X, int nx, double *yhat, double *c, 
		   int nc, double *res2, double *res, const double *y,
		   double toler, int *count, double *vcv, 
		   const double **g, double *aux3, 
		   double *param, int nparam, double *b,  
		   double *a0, double *alfa, double *beta, int q,
		   int p, double *h, double **dhdp, double *zt)
{
    int i, j;
    double d, d0, d1, d2, ll2, d3;
    double ll3, d12, d31, d23, dd;
    double di, ff, dm, ds;
    int iv, nexp;
    double a1s, a2s, a3s;
    int it1, it2, it3, it4, it5;
    double d12s, dac, dub, d23s, d31s;
    int err;
    double s_1, s_2, stre, bigd;
    double cappa;
    int ncall;
    int nvparm;
    double oldstp = 9.0e+39;
    static double ll1 = 0.0, fs = 0.0;

    double *gg, *step;

    gg = malloc(nparam * sizeof *gg);
    if (gg == NULL) {
	return 1;
    }

    step = malloc(nparam * sizeof *step);
    if (step == NULL) {
	free(gg);
	return 1;
    } 

    iv = 0;
    it1 = it2 = it3 = it4 = it5 = 0;
    nvparm = q + p + 1;

    /* calculate information matrix */
    vcv_setup(t1, t2, c, nc, res2, res, 
	      count, g, aux3, param, nparam, 
	      a0, alfa, beta, q, p, h,
	      dhdp, zt, NULL, vcv, 
	      VCV_IM);

    /* invert the information matrix */
    err = invert(vcv, nparam);
    if (err) {
	fprintf(stderr, "garch_info_matrix: matrix inversion failed\n");
	return 1;
    }

    if (count == NULL) {
	/* just calculating vcv at convergence */
	goto vcv_exit;
    }

    /* here we start iteration */

    /* calculate the step for the new coefficients */
    s_2 = 0.0;
    for (i = 0; i < nparam; ++i) {
	gg[i] = param[i];
	step[i] = 0.0;
	for (j = 0; j < nparam; ++j) {
	    step[i] -= vcv[vix(i,j)] * aux3[j];
	}
	s_2 += step[i] * step[i];
    }

    /* if rel. Euclidean distance is used as converg. */
    s_1 = 0.0;
    for (i = 0; i < nparam; ++i) {
	s_1 += param[i] * param[i];
    }
    if (s_1 == 0.0) s_1 = 1e-10;
    stre = s_2 / s_1;
    s_2 = sqrt(s_2);
    stre = sqrt(stre);

#ifdef FDEBUG
    fprintf(stderr, "s_2 = %.9g\n", s_2);
#endif

    oldstp = s_2;
    for (i = 0; i < nparam; ++i) {
	step[i] /= s_2;
    }

    it4 += iv;
    ds = s_2;
    if (stre <= toler) {
	goto L496;
    }

    ncall = 0;
    nexp = *count / 5;
    if (nexp > 5) {
	nexp = 5;
    }

    cappa = pow(2.0, nexp);
    d0 = s_2;
    d0 /= cappa;
    dac = d0 * .001;
    dub = d0 * 4.0;

    if (*count == 1) {
	ll1 = -garch_ll(c, nc, res2, res, yhat, y, 
			X, nx, t1, t2, param, b, a0, alfa, 
			beta, q, p, h);
#ifdef FFDEBUG
	fprintf(stderr, "count=1, ll1=%.9g\n", ll1);
#endif
    }
    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d0;
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(1) = %.9g\n", i, param[i]);
    }
#endif 

    ll2 = -garch_ll(c, nc, res2, res, yhat, y, 
		    X, nx, t1, t2, param, b, a0, alfa, 
		    beta, q, p, h);
#ifdef FDEBUG
    fprintf(stderr, "ll2=%.9g, ll1=%.9g\n", ll2, ll1);
#endif    
    if (ll2 > ll1) {
	goto L307;
    }

    d1 = 0.0;
    d2 = d0;
    d3 = d0 + d0;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d3;
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(2) = %.9g\n", i, param[i]);
    }
#endif    

    ll3 = -garch_ll(c, nc, res2, res, yhat, y, 
		    X, nx, t1, t2, param, b, a0, alfa, 
		    beta, q, p, h);
    goto L325;

 L307:
    d1 = -d0;
    d2 = 0.0;
    d3 = d0;
    ll3 = ll2;
    ll2 = ll1;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d1;
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(3) = %.9g\n", i, param[i]);
    }
#endif 

    ll1 = -garch_ll(c, nc, res2, res, yhat, y, 
		    X, nx, t1, t2, param, b, a0, alfa, 
		    beta, q, p, h);

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
    d1 -= dub;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + d1 * step[i];
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(4) = %.9g\n", i, param[i]);
    }
#endif 

    ll1 = -garch_ll(c, nc, res2, res, yhat, y, 
		    X, nx, t1, t2, param, b, a0, alfa, 
		    beta, q, p, h);

    if (++ncall > 100) {
	goto L490;
    }
    goto L325;

 L341:
    d1 = d2;
    ll1 = ll2;
    d2 = d3;
    ll2 = ll3;
    d3 += dub;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + d3 * step[i];
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(5) = %.9g\n", i, param[i]);
    }
#endif

    ll3 = -garch_ll(c, nc, res2, res, yhat, y, 
		    X, nx, t1, t2, param, b, a0, alfa, 
		    beta, q, p, h);

    if (++ncall > 100) {
	goto L490;
    }
    goto L325;

 L400:
    d23s = d23 * (d2 + d3);
    d31s = d31 * (d3 + d1);
    d12s = d12 * (d1 + d2);
    ds = (d23s * ll1 + d31s * ll2 + d12s * ll3) * .5 / di;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * ds;
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(6) = %.9g\n", i, param[i]);
    }
#endif

    fs = -garch_ll(c, nc, res2, res, yhat, y, 
		   X, nx, t1, t2, param, b, a0, alfa, 
		   beta, q, p, h);

    if (++ncall > 100) {
	goto L490;
    }

    a1s = (d = d1 - ds, fabs(d));
    a2s = (d = d2 - ds, fabs(d));
    a3s = (d = d3 - ds, fabs(d));

    dm = a1s;
    if (a3s < dm) {
	dm = a3s;
    }
    if (dub >= dm) {
	goto L422;
    }
    if (ds < d1 - dub) {
	goto L329;
    }
    if (ds > d3 + dub) {
	goto L341;
    }

 L422:
    if (a1s < dac || a2s < dac || a3s < dac) {
	goto L490;
    }
    if (ll1 < ll2 || ll1 < ll3) {
	goto L434;
    }
    d1 = ds;
    ll1 = fs;
    goto L459;

 L434:
    if (ll2 < ll3 || ll2 < ll1) {
	goto L447;
    }
    d2 = ds;
    ll2 = fs;
    goto L459;

 L447:
    d3 = ds;
    ll3 = fs;

 L459:
    if (d2 <= d3) {
	goto L463;
    }
    dd = d2;
    ff = ll2;
    d2 = d3;
    ll2 = ll3;
    d3 = dd;
    ll3 = ff;

 L463:
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

 L490:
    if (fs <= ll1) {
	goto L491;
    }
    fs = ll1;
    ds = d1;

 L491:
    if (fs <= ll2) {
	goto L492;
    }
    fs = ll2;
    ds = d2;

 L492:
    if (fs <= ll3) {
	goto L496;
    }
    fs = ll3;
    ds = d3;

 L496:
    for (i = 0; i <nparam; ++i) {
	param[i] = gg[i] + ds * step[i];
    }
    check_ht(param + nc, nvparm);

#ifdef FDEBUG
    for (i=0; i<nparam; i++) {
	fprintf(stderr, "param[%d] in matinf(7) = %.9g\n", i, param[i]);
    }
#endif

    ll1 = fs;
    fs = -fs;
    it5 += iv;

 vcv_exit:

    /* change the sign of the matrix */
    for (i = 0; i < nparam; ++i) {
	for (j = 0; j < nparam; ++j) {
	    vcv[vix(i,j)] *= -1.0;
	}
    }

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
    
    for (i=0; i<np; i++) {
	H[i] = malloc(np * sizeof **H);
	for (j=0; j<np; j++) {
	    H[i][j] = malloc((lag + 1) * sizeof ***H);
	}
    }

    return H;
}

static int 
garch_full_hessian (int t1, int t2, 
		    const double **X, int nx, double *yhat, double *c, 
		    int nc, double *res2, double *res, const double *y,
		    double toler, int *count, double *vcv, 
		    const double **g, double *aux3, 
		    double *param, int nparam, double *b, 
		    double *a0, double *alfa, double *beta, int q,
		    int p, double *h, double **dhdp, double *zt)
{
    int i, j;
    double d;
    double d0, d1, d2, ll2, d3; 
    double ll3, d12, d31, d23, dd, di, ff; 
    int iv, ncall, nexp, nvparm;
    double dm, ds; 
    double a1s, a2s, a3s;
    int it1, it2, it3, it4, it5;
    double dac, d12s, dub, d23s, d31s;
    int err;
    double s_1, s_2, stre, bigd;
    double cappa;
    int lag;
    static double ll1 = 0.0, fs = 0.0;
    double oldstp = 9.0e+39;

    double *gg, *step;
    double ***H;

    gg = malloc(nparam * sizeof *gg);
    if (gg == NULL) {
	return 1;
    }

    step = malloc(nparam * sizeof *step);
    if (step == NULL) {
	free(gg);
	return 1;
    } 

    /* 3rd dimension of dhdpdp is max(p,q) + 1 */
    H = allocate_dhdpdp(nparam, p, q); 
    if (H == NULL) {
	free(gg);
	free(step);
	return 1;
    }

    iv = 0;
    it1 = it2 = it3 = it4 = it5 = 0;
    nvparm = q + p + 1;
    lag = (p > q)? p : q;

    /* calculate the full Hessian */
    vcv_setup(t1, t2, c, nc, res2, res, 
	      count, g, aux3, param, nparam, 
	      a0, alfa, beta, q, p, h,
	      dhdp, zt, H, vcv, 
	      VCV_HESSIAN);

    /* invert the Hessian */
    err = invert(vcv, nparam);
    if (err) {
	fprintf(stderr, "garch_full_hessian: matrix inversion failed\n");
    }

    /* Start iteration here */

    /* calculate the step for the new coefficients */
    s_2 = 0.0;
    for (i = 0; i < nparam; ++i) {
	gg[i] = param[i];
	step[i] = 0.0;
	for (j = 0; j < nparam; ++j) {
	    step[i] -= vcv[vix(i,j)] * aux3[j];
	}
	s_2 += step[i] * step[i];
    }

    /* if rel. Euclidean distance is used as converg. */
    s_1 = 0.0;
    for (i = 0; i < nparam; ++i) {
	s_1 += param[i] * param[i];
    }
    if (s_1 == 0.0) s_1 = 1e-10;
    stre = s_2 / s_1;
    s_2 = sqrt(s_2);
    stre = sqrt(stre);

    oldstp = s_2;
    for (i = 0; i < nparam; ++i) {
	step[i] /= s_2;
    }

    it4 += iv;
    ds = s_2;
    if (stre <= toler) {
	goto L496;
    }

    ncall = 0;
    nexp = *count / 5;
    if (nexp > 5) {
	nexp = 5;
    }

    cappa = pow(2.0, nexp);
    d0 = s_2;
    d0 /= cappa;
    dac = d0 * .001;
    dub = d0 * 4.0;

    if (*count == 1) {
	ll1 = -garch_ll(c, nc, res2, res, yhat, y, 
			X, nx, t1, t2, param, b, a0, alfa, 
			beta, q, p, h);
#ifdef FDEBUG
	fprintf(stderr, "hess: ll1 = %.9g\n", ll1);
#endif  
    }

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d0;
    }
    check_ht(param + nc, nvparm);

    ll2 = -garch_ll(c, nc, res2, res, yhat, y, 
		    X, nx, t1, t2, param, b, a0, alfa, 
		    beta, q, p, h);
#ifdef FDEBUG
    fprintf(stderr, "hess: ll2 = %.9g\n", ll2);
#endif    
    if (ll2 > ll1) {
	goto L307;
    }

    d1 = 0.0;
    d2 = d0;
    d3 = d0 + d0;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d3;
    }
    check_ht(param + nc, nvparm);

    ll3 = -garch_ll(c, nc, res2, res, yhat, y, 
		    X, nx, t1, t2, param, b, a0, alfa, 
		    beta, q, p, h);
#ifdef FDEBUG
    fprintf(stderr, "hess: ll3 = %.9g\n", ll3);
#endif    
    goto L325;

 L307:
    d1 = -d0;
    d2 = 0.0;
    d3 = d0;
    ll3 = ll2;
    ll2 = ll1;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * d1;
    }
    check_ht(param + nc, nvparm);

    ll1 = -garch_ll(c, nc, res2, res, yhat, y, 
		    X, nx, t1, t2, param, b, a0, alfa, 
		    beta, q, p, h);
#ifdef FDEBUG
    fprintf(stderr, "hess: ll1(2) = %.9g\n", ll1);
#endif    


 L325:
    d23 = d2 - d3;
    d31 = d3 - d1;
    d12 = d1 - d2;
    di = d23 * ll1 + d31 * ll2 + d12 * ll3;
    bigd = di * -2.0 / (d23 * d31 * d12);
#ifdef FDEBUG
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
    d1 -= dub;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + d1 * step[i];
    }
    check_ht(param + nc, nvparm);

    ll1 = -garch_ll(c, nc, res2, res, yhat, y, 
		    X, nx, t1, t2, param, b, a0, alfa, 
		    beta, q, p, h);
#ifdef FDEBUG
    fprintf(stderr, "hess: ll1(3) = %.9g\n", ll1);
#endif    

    if (++ncall > 100) {
	goto L490;
    }
    goto L325;

 L341:
    d1 = d2;
    ll1 = ll2;
    d2 = d3;
    ll2 = ll3;
    d3 += dub;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + d3 * step[i];
    }
    check_ht(param + nc, nvparm);

    ll3 = -garch_ll(c, nc, res2, res, yhat, y, 
		    X, nx, t1, t2, param, b, a0, alfa, 
		    beta, q, p, h);
#ifdef FDEBUG
    fprintf(stderr, "hess: ll3(2) = %.9g\n", ll3);
#endif    

    if (++ncall > 100) {
	goto L490;
    }
    goto L325;

 L400:
    d23s = d23 * (d2 + d3);
    d31s = d31 * (d3 + d1);
    d12s = d12 * (d1 + d2);
    ds = (d23s * ll1 + d31s * ll2 + d12s * ll3) * .5 / di;

    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + step[i] * ds;
    }
    check_ht(param + nc, nvparm);

    fs = -garch_ll(c, nc, res2, res, yhat, y, 
		   X, nx, t1, t2, param, b, a0, alfa, 
		   beta, q, p, h);
#ifdef FDEBUG
    fprintf(stderr, "hess: fs = %.9g, ds = %.9g\n", fs, ds);
#endif    

    if (++ncall > 100) {
	goto L490;
    }

    a1s = (d = d1 - ds, fabs(d));
    a2s = (d = d2 - ds, fabs(d));
    a3s = (d = d3 - ds, fabs(d));

    dm = a1s;
    if (a3s < dm) {
	dm = a3s;
    }
    if (dub >= dm) {
	goto L422;
    }
    if (ds < d1 - dub) {
	goto L329;
    }
    if (ds > d3 + dub) {
	goto L341;
    }

 L422:
    if (a1s < dac || a2s < dac || a3s < dac) {
	goto L490;
    }
    if (ll1 < ll2 || ll1 < ll3) {
	goto L434;
    }
    d1 = ds;
    ll1 = fs;
    goto L459;

 L434:
    if (ll2 < ll3 || ll2 < ll1) {
	goto L447;
    }
    d2 = ds;
    ll2 = fs;
    goto L459;

 L447:
    d3 = ds;
    ll3 = fs;

 L459:
    if (d2 <= d3) {
	goto L463;
    }
    dd = d2;
    ff = ll2;
    d2 = d3;
    ll2 = ll3;
    d3 = dd;
    ll3 = ff;

 L463:
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
 
L490:
    if (fs <= ll1) {
	goto L491;
    }
    fs = ll1;
    ds = d1;

 L491:
    if (fs <= ll2) {
	goto L492;
    }
    fs = ll2;
    ds = d2;

 L492:
    if (fs <= ll3) {
	goto L496;
    }
    fs = ll3;
    ds = d3;

 L496:
    for (i = 0; i < nparam; ++i) {
	param[i] = gg[i] + ds * step[i];
    }

    check_ht(param + nc, nvparm);

    ll1 = fs;
    fs = -fs;
    it5 += iv;

    /* change the sign of the inverse */
    for (i = 0; i < nparam; ++i) {
	for (j = 0; j < nparam; ++j) {
	    vcv[vix(i,j)] *= -1.0;
	}
    }

    free(gg);
    free(step);
    free_dhdpdp(H, nparam);

    return 0;
} /* garch_full_hessian */

static int invert (double *g, int dim)
{
    integer m = dim;
    integer info, lwork;
    integer *ipiv;
    double *work;

#ifdef FDEBUG
    int i, j;
    static int k;

    fprintf(stderr, "Matrix inversion %d, dim=%d\n", k, dim);

    for (i=0; i<dim; i++) {
	for (j=0; j<dim; j++) {
	    fprintf(stderr, "%14.9g ", g[i + dim * j]);
	}
	fputc('\n', stderr);
    }
    k++;
#endif

    ipiv = malloc(m * sizeof *ipiv);
    if (ipiv == NULL) {
	return 1;
    }

    work = malloc(sizeof *work);
    if (work == NULL) {
	free(ipiv);
	return 1;
    }  

    dgetrf_(&m, &m, g, &m, ipiv, &info);   

    if (info != 0) {
	free(ipiv);
	free(work);
	return 1;
    }

    lwork = -1;
    dgetri_(&m, g, &m, ipiv, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	free(ipiv);
	free(work);
	return 1;
    }

    lwork = (integer) work[0];

    work = realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	free(ipiv);
	return 1;
    }  

    dgetri_(&m, g, &m, ipiv, work, &lwork, &info);

    free(work);
    free(ipiv);

    return (int) info;
}

