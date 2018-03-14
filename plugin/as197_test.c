struct as197_info {
    int p;
    int P;
    int q;
    int Q;
    int pd;
    int plen;
    int qlen;
    int r;
    int rp1;
    int ifault;
    int n;
    int ifc;
    double *phi, *theta;
    double *w, *w0, *e;
    double *vw, *vl, *vk;
    double sumsq, fact;
    double toler;
    double mu;
    double loglik;
    int verbose;
    int ncalls, nbad;
    int use_loglik;
    arma_info *ai;
};

static void as197_info_init (struct as197_info *as,
			     arma_info *ai,
			     int verbose,
			     double delta,
			     int use_loglik)
{
    as->p = ai->p;
    as->P = ai->P;
    as->q = ai->q;
    as->Q = ai->Q;
    as->pd = ai->pd;
    as->n = ai->T;
    as->ifc = ai->ifc;

    as->plen = as->p + as->pd * as->P;
    as->qlen = as->q + as->pd * as->Q;
    as->r = (as->plen > as->qlen + 1)? as->plen : as->qlen + 1;
    as->rp1 = as->r + 1;

    as->phi = as->theta = NULL;
    if (as->plen > 0) {
	as->phi = malloc(as->plen * sizeof *as->phi);
    }
    if (as->qlen > 0) {
	as->theta = malloc(as->qlen * sizeof *as->theta);
    }

    as->e =  malloc(as->n * sizeof *as->e);
    as->vw = malloc(as->rp1 * sizeof *as->vw);
    as->vl = malloc(as->rp1 * sizeof *as->vl);
    as->vk = malloc(as->r * sizeof *as->vk);

    as->w = as->w0 = NULL; /* later! */

    as->toler = delta;
    as->mu = 0;
    as->loglik = NADBL;
    as->ifault = 0;

    as->verbose = verbose > 1;
    as->ncalls = as->nbad = 0;
    as->use_loglik = use_loglik;
}

static void as197_info_free (struct as197_info *as)
{
    free(as->phi);
    free(as->theta);
    free(as->e);
    free(as->vw);
    free(as->vl);
    free(as->vk);
    free(as->w0);
}

static void write_big_phi_197 (const double *b,
			       struct as197_info *as)
{
    const double *bs = b + as->p;
    double x, y;
    int i, j, k, ii;

    for (i=0; i<as->plen; i++) {
	as->phi[i] = 0.0;
    }

    for (j=-1; j<as->P; j++) {
	x = (j < 0)? -1 : bs[j];
	k = 0;
        for (i=-1; i<as->p; i++) {
	    if (i < 0) {
		y = -1;
	    } else {
		y = b[k++];
	    }
            ii = (j+1) * as->pd + (i+1);
	    if (ii > 0) {
		as->phi[ii-1] -= x * y;
	    }
        }
    }
}

static void write_big_theta_197 (const double *b,
				 struct as197_info *as)
{
    const double *bs = b + as->q;
    double x, y;
    int i, j, k, ii;

    for (i=0; i<as->qlen; i++) {
	as->theta[i] = 0.0;
    }

    for (j=-1; j<as->Q; j++) {
	x = (j < 0)? 1 : bs[j];
	k = 0;
        for (i=-1; i<as->q; i++) {
	    if (i < 0) {
		y = 1;
	    } else {
		y = b[k++];
	    }
            ii = (j+1) * as->pd + (i+1);
	    if (ii > 0) {
		as->theta[ii-1] = x * y;
	    }
        }
    }
}

static double as197_iteration (const double *b, void *data)
{
    struct as197_info *as = data;
    double crit = NADBL;
    int np = as->p + as->P;
    int i;

    as->ncalls += 1;

    if (as->q > 0 || as->Q > 0) {
	/* check that MA term(s) are within bounds */
	const double *theta = b + as->ifc + np;
	const double *Theta = theta + as->q;

	if (ma_out_of_bounds(as->ai, theta, Theta)) {
	    as->nbad += 1;
	    return crit;
	}
    }

    if (as->ifc) {
	/* subtract the constant */
	as->mu = b[0];
	for (i=0; i<as->n; i++) {
	    as->w[i] = as->w0[i] - as->mu;
	}
	b++;
    }
    if (as->P > 0) {
	write_big_phi_197(b, as);
    } else if (as->p > 0) {
	for (i=0; i<as->p; i++) {
	    as->phi[i] = b[i];
	}
    }
    b += np;
    if (as->Q > 0) {
	write_big_theta_197(b, as);
    } else if (as->q > 0) {
	for (i=0; i<as->q; i++) {
	    as->theta[i] = b[i];
	}
    }

    as->ifault = flikam(as->phi, as->plen, as->theta, as->qlen,
			as->w, as->e, as->n, &as->sumsq, &as->fact,
			as->vw, as->vl, as->rp1, as->vk, as->r,
			as->toler);

    if (as->ifault > 0) {
	if (as->ifault == 5) {
	    ; // fputs("flikam: (near) non-stationarity\n", stderr);
	} else {
	    fprintf(stderr, "flikam: ifault = %d\n", as->ifault);
	}
	as->nbad += 1;
	return NADBL;
    }

    if (isnan(as->sumsq) || isnan(as->fact)) {
	as->nbad += 1; /* leave crit as NA */
    } else {
	/* The criterion used by Melard may work better than
	   the full loglikelihood in the context of his
	   algorithm?
	*/
	if (as->use_loglik) {
	    /* full loglikelihood */
	    double ll1 = 1.0 + LN_2_PI + log(as->sumsq / as->n);
	    double sumldet = as->n * log(as->fact);

	    crit = -0.5 * (as->n * ll1 + sumldet);
	} else {
	    /* Melard's criterion */
	    crit = -as->fact * as->sumsq;
	}
    }

    if (as->verbose) {
	printf("flikam: ssq=%#.12g, fact=%#.12g, crit=%#.12g\n",
	       as->sumsq, as->fact, crit);
    }

    return crit;
}

/* calculate the full loglikelihood on completion */

static void as197_full_loglik (struct as197_info *as)
{
    double ll1 = 1.0 + LN_2_PI + log(as->sumsq / as->n);
    double sumldet = as->n * log(as->fact);

    as->loglik = -0.5 * (as->n * ll1 + sumldet);
}

static int as197_arma_finish (MODEL *pmod,
			      arma_info *ainfo,
			      const DATASET *dset,
			      struct as197_info *as,
			      double *b, gretlopt opt,
			      PRN *prn)
{
    int i, t, k = ainfo->nc;
    double s2;
    int err;

    pmod->t1 = ainfo->t1;
    pmod->t2 = ainfo->t2;
    pmod->nobs = ainfo->T;
    pmod->ncoeff = ainfo->nc;
    pmod->full_n = dset->n;

    err = gretl_model_allocate_storage(pmod);
    if (err) {
	return err;
    }

    for (i=0; i<k; i++) {
	pmod->coeff[i] = b[i];
    }

    s2 = 0.0;
    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	s2 += as->e[i] * as->e[i]; /* ?? */
	pmod->uhat[t] = as->e[i++]; /* gretl_vector_get(kh->E, i++); */
    }

    pmod->sigma = sqrt(s2 / ainfo->T);
    pmod->lnL = as->loglik;

    if (1) { /* !do_opg */
	/* base covariance matrix on Hessian (FIXME perhaps QML) */
	gretl_matrix *Hinv;

	as->use_loglik = 1;
	Hinv = numerical_hessian_inverse(b, ainfo->nc, as197_iteration,
					 as, &err);
	if (!err) {
	    err = gretl_model_write_vcv(pmod, Hinv);
	    if (!err) {
		gretl_model_set_vcv_info(pmod, VCV_ML, ML_HESSIAN);
	    }
	}
	gretl_matrix_free(Hinv);
    }

    if (!err) {
	write_arma_model_stats(pmod, ainfo, dset);
	arma_model_add_roots(pmod, ainfo, b);
	gretl_model_set_int(pmod, "arma_flags", ARMA_EXACT);
	gretl_model_set_int(pmod, "AS197", 1);
	if (arima_ydiff_only(ainfo)) {
	    pmod->opt |= OPT_Y;
	}
    }

    return err;
}

static int as197_undo_y_scaling (arma_info *ainfo,
				 gretl_matrix *y,
				 double *b,
				 struct as197_info *as)
{
    double *beta = b + 1 + ainfo->np + ainfo->P +
	ainfo->nq + ainfo->Q;
    int i, t, T = ainfo->t2 - ainfo->t1 + 1;
    int err = 0;

    /* adjust the constant */
    b[0] /= ainfo->yscale;

    for (i=0; i<ainfo->nexo; i++) {
	beta[i] /= ainfo->yscale;
    }

    i = ainfo->t1;
    for (t=0; t<T; t++) {
	y->val[t] /= ainfo->yscale;
    }

    if (na(as197_iteration(b, as))) {
	err = 1;
    }

    return err;
}

static int as197_arma (const double *coeff,
		       const DATASET *dset,
		       arma_info *ainfo,
		       MODEL *pmod,
		       gretlopt opt)
{
    struct as197_info as;
    gretl_matrix *y = NULL;
    double *b = NULL;
    double delta = -1.0;
    int use_loglik = 0;
    int verbose = 1;
    int err = 0;

    as197_info_init(&as, ainfo, verbose, delta, use_loglik);
    as.ai = ainfo; /* link */

    b = copyvec(coeff, ainfo->nc);
    if (b == NULL) {
	return E_ALLOC;
    }

    /* FIXME scaling of y */
    y = form_arma_y_vector(ainfo, &err);

    if (!err) {
	as.w = y->val;
	if (as.ifc) {
	    as.w0 = copyvec(as.w, as.n);
	    if (as.w0 == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (!err) {
	/* maximize loglikelihood via BFGS */
	int maxit;
	double toler;
	int fncount = 0;
	int grcount = 0;
	int nparam = as.p + as.q + as.P + as.Q + as.ifc;

	if (as.n > 2000) {
	    /* try to avoid slowdown on big samples */
	    as.toler = 0.0001;
	    as.use_loglik = 1; /* ? */
	}

	BFGS_defaults(&maxit, &toler, ARMA);

	err = BFGS_max(b, nparam, maxit, toler,
		       &fncount, &grcount, as197_iteration, C_LOGLIK,
		       NULL, &as, NULL, opt, ainfo->prn);
	if (!err) {
	    if (ainfo->yscale != 1.0) {
		/* broken? */
		as.use_loglik = 1;
		as197_undo_y_scaling(ainfo, y, b, &as);
	    } else if (!as.use_loglik) {
		as197_full_loglik(&as);
	    }
	    gretl_model_set_int(pmod, "fncount", fncount);
	    gretl_model_set_int(pmod, "grcount", grcount);
	    err = as197_arma_finish(pmod, ainfo, dset, &as, b,
				    opt, ainfo->prn);
	}
    }

    if (err && !pmod->errcode) {
	pmod->errcode = err;
    }

    as197_info_free(&as);
    gretl_matrix_free(y);
    free(b);

    return err;
}

/* As of 2018-03-13, the AS197 implementation for gretl
   can't handle missing values, "gappy" non-seasonal
   AR or MA specifications, or exogenous variables
   (ARMAX). So we need to screen out these conditions
   before saying OK to using the testing code.
   Added: y-scaling doesn't work yet either?
*/

static int as197_ok (arma_info *ainfo)
{
    if (ainfo->nexo > 0) {
	/* exogenous vars included */
	return 0;
    } else if (arma_missvals(ainfo)) {
	/* NAs in sample range */
	return 0;
    } else if (ainfo->pqspec != NULL && *ainfo->pqspec != '\0') {
	/* marker for "gappy" case */
	return 0;
    } else if (ainfo->yscale != 1.0) {
	/* not working yet! */
	return 0;
    } else {
	return 1;
    }
}
