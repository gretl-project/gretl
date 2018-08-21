static int inclu2_(int, int, const double *,
		   double *, double, double *,
		   double *, double *, double *,
		   double *, int *);

static int regres_(int, int, const double *, const double *,
		   double *);

int starma (int ip, int iq, int ir, int np, double *phi,
	    double *theta, double *a, double *p0, double *v,
	    double *thetab, double *xnext, double *xrow,
	    double *rbar, int nrbar)
{
    int ind, npr, ind1, ind2, npr1, indi, indj, indn;
    double vj, phii, phij;
    int i, j, irank;
    double ynext, recres, ssqerr;

    /* set A(0), v and phi */

    for (i=1; i<ir; i++) {
	a[i] = 0;
	if (i >= ip) {
	    phi[i] = 0;
	}
	v[i] = 0;
	if (i <= iq) {
	    v[i] = theta[i-1];
	}
    }
    a[0] = 0;
    if (ip == 0) {
	phi[0] = 0;
    }
    v[0] = 1;
    ind = ir-1;
    for (j=1; j<ir; j++) {
	vj = v[j];
	for (i=j; i<ir; i++) {
	    v[++ind] = v[i] * vj;
	}
    }

    if (ip == 0) {
	/* pure MA */
	goto backsub;
    }

    irank = 0;
    ssqerr = 0;
    for (i=0; i<nrbar; i++) {
	rbar[i] = 0;
    }
    for (i=0; i<np; i++) {
	p0[i] = 0;
	thetab[i] = 0;
	xnext[i] = 0;
    }

    ind = -1;
    ind1 = -1;
    npr = np - ir;
    npr1 = npr + 1;
    indj = npr1 - 1;
    ind2 = npr - 1;

    for (j=0; j<ir; j++) {
	phij = phi[j];
	xnext[indj++] = 0;
	indi = npr1 + j;
	for (i=j; i<ir; i++) {
	    ynext = v[++ind];
	    phii = phi[i];
	    if (j < ir-1) {
		xnext[indj] = -phii;
		if (i < ir-1) {
		    xnext[indi] -= phij;
		    xnext[++ind1] = -1;
		}
	    }
	    xnext[npr] = -phii * phij;
	    if (++ind2 >= np) {
		ind2 = 0;
	    }
	    xnext[ind2] += 1;
	    inclu2_(np, nrbar, xnext, xrow, ynext, p0, rbar,
		    thetab, &ssqerr, &recres, &irank);
	    xnext[ind2] = 0;
	    if (i < ir-1) {
		xnext[indi++] = 0;
		xnext[ind1] = 0;
	    }
	}
    }

    regres_(np, nrbar, rbar, thetab, p0);

    /* reorder p0 */
    ind = npr-1;
    for (i=0; i<ir; i++) {
	xnext[i] = p0[++ind];
    }
    ind = np-1;
    ind1 = npr-1;
    for (i=0; i<npr; i++) {
	p0[ind--] = p0[ind1--];
    }
    for (i=0; i<ir; i++) {
	p0[i] = xnext[i];
    }
    return 0;

 backsub:

    indn = np;
    ind = np;
    for (i=0; i<ir; i++) {
	for (j=0; j<=i; j++) {
	    ind--;
	    p0[ind] = v[ind];
	    if (j > 0) {
		indn--;
		p0[ind] += p0[indn];
	    }
	}
    }

    return 0;
}

/* Handle the required transformation of the P matrix
   in case the dependent variable is missing at the
   prior observation. In principle this could be done
   without constructing the full T and P matrices --
   exploiting the symmetry of P and the special
   structure of T -- but it would be quite fiddly!
   See the Gardner et al AS 154 paper for details.
*/

static int handle_missing_obs (int ip, int ir, int np,
			       const double *phi,
			       const double *v,
			       double *p0)
{
    gretl_matrix_block *B;
    gretl_matrix *T, *P, *TPT;
    gretl_matrix m;
    int i, j, k, err = 0;

    B = gretl_matrix_block_new(&T, ir, ir,
			       &P, ir, ir,
			       &TPT, ir, ir,
			       NULL);
    if (B == NULL) {
	return E_ALLOC;
    }

    /* construct the full T matrix using @phi */
    gretl_matrix_zero(T);
    for (i=0; i<ip; i++) {
	T->val[i] = phi[i];
    }
    gretl_matrix_inscribe_I(T, 0, 1, ir-1);

    /* construct full P matrix using @p0 */
    gretl_matrix_init(&m);
    m.rows = np; m.cols = 1; m.val = p0;
    err = gretl_matrix_unvectorize_h(P, &m);

    if (!err) {
	/* form T*P*T' */
	err = gretl_matrix_qform(T, GRETL_MOD_NONE, P,
				 TPT, GRETL_MOD_NONE);
    }

    if (!err) {
	/* p0 <- TPT' + V (= RR'), in vech mode */
	k = 0;
	for (j=0; j<ir; j++) {
	    for (i=j; i<ir; i++) {
		p0[k] = gretl_matrix_get(TPT, i, j) + v[k];
		k++;
	    }
	}
    }

    gretl_matrix_block_destroy(B);

    return err;
}

int karma (int ip, int iq, int ir, int np, double *phi,
	   double *theta, double *a, double *p0,
	   double *v, int n, double *w, double *resid,
	   double *sumlog, double *ssq, int iupd,
	   double delta, double *e, int *nit)
{
    double g, a0;
    int i, j, l, ii;
    double dt, et, ft, ut, wnext;
    int ind, indn, indw, inde = 0;
    int reset_i = 1;
    int ir1 = ir - 1;
    int skip = 0;

    for (i=0; i<ir; i++) {
	e[i] = 0;
    }

    /* for non-zero values of @nit, perform quick recursions */
    if (*nit != 0) {
	goto quick_recurse;
    }

    for (i=0; i<n; i++) {
	wnext = w[i];
	/* prediction */
	if (iupd == 0 || i > 0) {
	    dt = 0;
	    if (ir != 1) {
		dt = p0[ir];
	    }
	    if (dt < delta) {
		reset_i = 0;
		goto quick_recurse;
	    }
	    a0 = a[0];
	    if (ir != 1) {
		for (j=0; j<ir1; j++) {
		    a[j] = a[j+1];
		}
	    }
	    a[ir-1] = 0;
	    for (j=0; j<ip; j++) {
		a[j] += phi[j] * a0;
	    }
	    if (skip) {
		/* gretl addition to algorithm */
		handle_missing_obs(ip, ir, np, phi, v, p0);
	    } else {
		ind = -1;
		indn = ir-1;
		for (l=0; l<ir; l++) {
		    for (j=l; j<ir; j++) {
			ind++;
			p0[ind] = v[ind];
			if (j != ir-1) {
			    indn++;
			    p0[ind] += p0[indn];
			}
		    }
		}
	    }
	}
	if (isnan(wnext)) {
	    resid[i] = NADBL;
	    skip = 1;
	    continue;
	}
	/* updating */
	ft = p0[0];        /* MSE */
	ut = wnext - a[0]; /* unscaled forecast error */
	if (ir != 1) {
	    ind = ir-1;
	    for (j=1; j<ir; j++) {
		g = p0[j] / ft;
		a[j] += g * ut;
		for (l=j; l<ir; l++) {
		    p0[++ind] -= g * p0[l];
		}
	    }
	}
	a[0] = wnext;
	for (l=0; l<ir; l++) {
	    p0[l] = 0;
	}
#if 0	
	resid[i] = ut;
	e[inde++] = ut / sqrt(ft);
#else
	resid[i] = ut / sqrt(ft);
	e[inde++] = resid[i];
#endif
	if (inde >= iq) {
	    inde = 0;
	}
	*ssq += ut * ut / ft;
	*sumlog += log(ft);
	skip = 0;
    }

    *nit = n;
    return 0;

 quick_recurse:

    if (reset_i) {
	*nit = i = 0;
    } else {
	*nit = i - 1; /* ? */
    }

    for (ii=i; ii<n; ii++) {
	if (isnan(w[ii])) {
	    resid[ii] = 0;
	    continue;
	}
	et = w[ii];
	indw = ii;
	for (j=0; j<ip; j++) {
	    if (--indw >= 0) {
		et -= phi[j] * w[indw];
	    }
	}
	for (j=0; j<iq; ++j) {
	    if (--inde == -1) {
		inde = iq-1;
	    }
	    et -= theta[j] * e[inde];
	}
	e[inde] = resid[ii] = et;
	*ssq += et * et;
	if (++inde >= iq) {
	    inde = 0;
	}
    }

    return 0;
}

static int inclu2_ (int np, int nrbar, const double *xnext,
		    double *xrow, double y, double *p0, double *rbar,
		    double *thetab, double *ssqerr, double *recres,
		    int *irank)
{
    double di, xi, xk, dpi, cbar, sbar, rbthis;
    int i, k, i1, thisr;
    double wt = 1.0;

    for (i=0; i<np; i++) {
	xrow[i] = xnext[i];
    }

    *recres = 0.0;
    thisr = -1;

    for (i=0; i<np; i++) {
	if (xrow[i] == 0) {
	    thisr = thisr + np - i - 1;
	    continue;
	}
	xi = xrow[i];
	di = p0[i];
	dpi = di + wt * xi * xi;
	p0[i] = dpi;
	cbar = di / dpi;
	sbar = wt * xi / dpi;
	wt = cbar * wt;
	if (i < np-1) {
	    i1 = i + 1;
	    for (k=i1; k<np; ++k) {
		xk = xrow[k];
		thisr++;
		rbthis = rbar[thisr];
		xrow[k] = xk - xi * rbthis;
		rbar[thisr] = cbar * rbthis + sbar * xk;
	    }
	}
	xk = y;
	y = xk - xi * thetab[i];
	thetab[i] = cbar * thetab[i] + sbar * xk;
	if (di == 0) {
	    *irank += 1;
	    return 0;
	}
    }

    *ssqerr += wt * y * y;
    *recres = y * sqrt(wt);

    return 0;
}

static int regres_(int np, int nrbar, const double *rbar,
		   const double *thetab, double *beta)
{
    int i, j, i1, jm;
    int thisr = nrbar-1;
    int im = np-1;
    double bi;

    for (i=0; i<np; i++) {
	bi = thetab[im];
	if (im < np-1) {
	    i1 = i-1;
	    jm = np-1;
	    for (j=0; j<=i1; j++) {
		bi -= rbar[thisr--] * beta[jm--];
	    }
	}
	beta[im--] = bi;
    }

    return 0;
}
