static int tw_acf (const double *phi, int p,
                   const double *theta, int q,
		   double *acf, int ma,
		   double *cvli, int mxpqp1,
		   double *alpha, int mxpq);

#define max0(i,j) ((i)>(j) ? (i) : (j))
#define min0(i,j) ((i)<(j) ? (i) : (j))

#define missing(x) (na(x) || isnan(x))

#define DSHRINK 0.0625
#define DEXPAND 16.0
#define DCDELTA 4.0

/*
  Algorithm AS 197 Applied Statistics (1984) Vol. 33, No. 1

  Computes a quantity monotonically related to the likelihood
  function of an autoregressive-moving average process,
  expressed as fact * sumsq.

  Sign of the MA coefficients switched to agree with the
  convention followed by gretl.
*/

int flikam (const double *phi, int p,
	    const double *theta, int q,
	    const double *W, double *E, int n,
	    double *sumsq, double *fact,
	    double *vw, double *vl, int rp1,
	    double *vk, int r, double toler)
{
    double eps1 = 1.0e-10;
    double detman = 1.0;
    double detcar = 0.0;
    double R, A, aor, wi;
    double vw0, vl0, alf, flj;
    int mxpq = max0(p, q);
    int mnpq = min0(p, q);
    int mxpqp1 = mxpq + 1;
    int do_quick = 0;
    int last, loop, jfrom, nexti;
    int i, j, k, ret = 0;

    *sumsq = *fact = 0;

    /* calculation of the autocovariance function of a process with
       unit innovation variance (vw) and the covariance between the
       variable and the lagged innovations (vl).
    */
    ret = tw_acf(phi, p, theta, q, vw, mxpqp1, vl, mxpqp1, vk, mxpq);
    if (ret > 0) {
	return ret;
    }

    /* computation of vk */
    vk[0] = vw[0];
    for (k=1; k<r; k++) {
	vk[k] = 0;
	for (j=k; j<p; j++) {
	    vk[k] += phi[j] * vw[j+1-k];
	}
	for (j=k; j<=q; j++) {
	    vk[k] += theta[j-1] * vl[j-k];
	}
    }

    /* computation of the initial vectors (vl, vk) */
    R = vk[0];
    if (R <= eps1) {
	return 8; /* fault code */
    }
    vl[r-1] = 0;
    for (j=0; j<r; j++) {
        vw[j] = 0;
	if (j < r-1) {
	    vl[j] = vk[j+1];
	}
	if (j < p) {
	    vl[j] += phi[j] * R;
	}
	vk[j] = vl[j];
    }

    /* exit if no observations */
    if (n <= 0) {
	return 9; /* fault code */
    }

    /* initialization for loop */
    last = p - q;
    loop = p;
    jfrom = p;
    vw[p] = 0;
    vl[mxpq] = 0;

    /* loop on time */
    for (i=0; i<n; i++) {
	/* test for skipping updating */
        if (i == last) {
	    loop = mnpq;
	    jfrom = loop;
	    /* test for switching */
	    if (q == 0) {
		/* pure AR */
		do_quick = 1;
		break;
	    }
	}
        if (i > mxpq && fabs(R - 1.0) < toler) {
	    do_quick = 1;
	    break;
	}
	/* update scalars */
        detman *= R;
	while (fabs(detman) >= 1.0) {
	    detman *= DSHRINK;
	    detcar += DCDELTA;
	}
	while (fabs(detman) < DSHRINK) {
	    detman *= DEXPAND;
	    detcar -= DCDELTA;
	}
	vw0 = vw[0];
        if (missing(W[i])) {
	    /* FIXME?! */
	    wi = vw0;
	    aor = A = E[i] = 0.0;
	} else {
	    wi = W[i];
	    A = W[i] - vw0;
	    E[i] = A / sqrt(R);
	    aor = A / R;
	    *sumsq += A * aor;
	}
        vl0 = vl[0];
        alf = vl0 / R;
        R -= alf * vl0;
	if (R <= eps1) {
	    return 8; /* fault code */
	}
 	/* update vectors */
        for (j=0; j<loop; j++) {
	    flj = vl[j+1] + phi[j] * vl0;
	    vw[j] = vw[j+1] + phi[j] * vw0 + aor * vk[j];
	    vl[j] = flj - alf * vk[j];
	    vk[j] -= alf * flj;
	}
	for (j=jfrom; j<q; j++) {
	    vw[j] = vw[j+1] + aor * vk[j];
	    vl[j] = vl[j+1] - alf * vk[j];
	    vk[j] -= alf * vl[j+1];
	}
	for (j=jfrom; j<p; j++) {
	    vw[j] = vw[j+1] + phi[j] * wi;
	}
    }

    if (do_quick) {
	/* quick recursions */
	nexti = i;
	ret = -nexti;
	for (i=nexti; i<n; i++) {
	    if (missing(W[i])) {
		continue;
	    }
	    E[i] = W[i];
	    for (j=0; j<p; j++) {
		E[i] -= phi[j] * W[i-j-1];
	    }
	    for (j=0; j<q; j++) {
		E[i] -= theta[j] * E[i-j-1];
	    }
	    *sumsq += E[i] * E[i];
	}
    }

    *fact = pow(detman, 1.0 / n) * pow(2.0, detcar / n);

    return ret;
}

/*
  Algorithm AS 197.1 Applied Statistics (1984) VOL. 33, NO. 1

  Implementation of the algorithm of G. Tunnicliffe-Wilson
  (J. Statist. Comput. Simul. 8, 1979, 301-309) for computation of the
  autocovariance function of an ARMA process of order (p, q) and unit
  innovation variance. The autoregressive and moving average
  coefficients are stored in vectors phi and theta, using Box and
  Jenkins notation. On output, vector cvli contains the covariances
  between the variable and the (K-1)-lagged innovation for K = 1, ...,
  q + 1.

  Dimensions: phi(p), theta(q), acf(ma), cvli(mxpqp1), alpha(mxpq)
*/

static int tw_acf (const double *phi, int p,
                   const double *theta, int q,
		   double *acf, int ma,
		   double *cvli, int mxpqp1,
		   double *alpha, int mxpq)
{
    double div, eps2 = 1.0e-10;
    int i, j, k, kc, j1;
    int kp1, mikp, miim1p;

    if (ma < mxpqp1) {
	return 4;
    }

    /* Initialization, and return if p = q = 0 */
    acf[0] = cvli[0] = 1.0;
    if (ma == 1) {
	return 0;
    }
    for (i=1; i<ma; i++) {
	acf[i] = 0;
    }
    if (mxpqp1 == 1) {
	return 0;
    }
    for (i=1; i<mxpqp1; i++) {
	cvli[i] = 0;
    }
    for (k=0; k<mxpq; k++) {
	alpha[k] = 0;
    }

    /* computation of the ACF of the moving average part,
       stored in @acf
    */
    for (k=0; k<q; k++) {
        cvli[k+1] = acf[k+1] = theta[k];
        kc = q - k - 1;
        for (j=0; j<kc; j++) {
	    acf[k+1] += theta[j] * theta[j+k+1];
	}
	acf[0] += theta[k] * theta[k];
    }

    if (p == 0) {
	return 0;
    }

    /* initialization of cvli */
    for (k=0; k<p; k++) {
        alpha[k] = cvli[k] = phi[k];
    }

    /* computation of T-W alpha and delta (delta stored in ACF
       which is gradually overwritten)
    */
    for (k=0; k<mxpq; k++) {
        kc = mxpq - k - 1;
        if (kc < p) {
	    div = 1.0 - alpha[kc] * alpha[kc];
	    if (div < eps2) {
		return 5;
	    }
	    if (kc == 0) {
		break;
	    }
	    for (j=0; j<kc; j++) {
		alpha[j] = (cvli[j] + alpha[kc] * cvli[kc-j-1]) / div;
	    }
	}
	if (kc < q) {
	    j1 = max0(kc - p, 1);
	    for (j=j1; j<=kc; j++) {
		acf[j] += acf[kc+1] * alpha[kc-j];
	    }
	}
	if (kc < p) {
	    for (j=0; j<kc; j++) {
		cvli[j] = alpha[j];
	    }
	}
    }

    /* computation of T-W Nu (Nu is stored in cvli,
       copied into acf)
    */
    acf[0] *= 0.5;
    for (k=0; k<mxpq; k++) {
        if (k < p) {
	    kp1 = k + 1;
	    div = 1.0 - alpha[k] * alpha[k];
	    for (j=0; j<=kp1; j++) {
		cvli[j] = (acf[j] + alpha[k] * acf[k+1-j]) / div;
	    }
	    for (j=0; j<=kp1; j++) {
		acf[j] = cvli[j];
	    }
	}
    }

    /* Computation of ACF (acf is gradually overwritten) */
    for (i=0; i<ma; i++) {
        miim1p = min0(i, p);
        for (j=0; j<miim1p; j++) {
	    acf[i] += phi[j] * acf[i-j-1];
	}
    }
    acf[0] *= 2.0;

    /* computation of cvli */
    cvli[0] = 1.0;
    for (k=0; k<q; k++) {
	cvli[k+1] = theta[k];
	if (p > 0) {
	    mikp = min0(k+1, p);
	    for (j=0; j<mikp; j++) {
		cvli[k+1] += phi[j] * cvli[k-j];
	    }
	}
    }

    return 0;
}
