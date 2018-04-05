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
    int ok_n;
    int ifc;
    double *phi, *theta;
    double *y, *y0, *e;   /* dependent var and forecast errors */
    double *vw, *vl, *vk; /* workspace */
    double sumsq, fact;   /* components of likelihood */
    double toler;         /* tolerance for switching to fast iterations */
    double loglik;
    int ma_check;
    int ncalls;
    int use_loglik;
    arma_info *ai;
    gretl_matrix *X;
    int free_X;
};

static int as197_info_init (struct as197_info *as,
			    arma_info *ai,
			    double toler,
			    int use_loglik)
{
    int err = 0;

    as->ai = ai; /* create accessor */

    /* convenience copies of @ai integer values */
    as->p = ai->p;
    as->P = ai->P;
    as->q = ai->q;
    as->Q = ai->Q;
    as->pd = ai->pd;
    as->n = ai->fullT;
    as->ok_n = ai->T;
    as->ifc = ai->ifc;

    as->plen = as->p + as->pd * as->P;
    as->qlen = as->q + as->pd * as->Q;
    as->r = (as->plen > as->qlen + 1)? as->plen : as->qlen + 1;
    as->rp1 = as->r + 1;

    as->y = as->y0 = NULL; /* later! */
    as->X = NULL; /* later too */
    as->free_X = 0;

    as->phi = as->theta = NULL;
    if (as->plen > 0) {
	as->phi = malloc(as->plen * sizeof *as->phi);
	if (as->phi == NULL) {
	    err = E_ALLOC;
	}
    }
    if (!err && as->qlen > 0) {
	as->theta = malloc(as->qlen * sizeof *as->theta);
	if (as->theta == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	int worklen = as->n + 2*as->rp1 + as->r;

	as->e =  malloc(worklen * sizeof *as->e);
	if (as->e == NULL) {
	    err = E_ALLOC;
	} else {
	    as->vw = as->e + as->n;
	    as->vl = as->vw + as->rp1;
	    as->vk = as->vl + as->rp1;
	}
    }

    if (!err) {
	as->toler = toler;
	as->loglik = NADBL;
	as->ifault = 0;
	as->ma_check = 0;
	as->ncalls = 0;
	as->use_loglik = use_loglik;
    }

    return err;
}

static void as197_info_free (struct as197_info *as)
{
    free(as->phi);
    free(as->theta);
    free(as->e);
    free(as->y0);

    if (as->free_X) {
	gretl_matrix_free(as->X);
    }
}

static void write_big_phi_197 (const double *b,
			       struct as197_info *as)
{
    const double *bs = b + as->ai->np;
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
	    } else if (AR_included(as->ai, i)) {
		y = b[k++];
	    } else {
		y = 0.0;
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
    const double *bs = b + as->ai->nq;
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
	    } else if (MA_included(as->ai, i)) {
		y = b[k++];
	    } else {
		y = 0.0;
	    }
            ii = (j+1) * as->pd + (i+1);
	    if (ii > 0) {
		as->theta[ii-1] += x * y;
	    }
        }
    }
}

static void as197_fill_arrays (struct as197_info *as,
			       const double *b)
{
    int np = as->ai->np + as->P;
    int nq = as->ai->nq + as->Q;
    double mu = 0.0;
    int i, j;

    if (as->ifc) {
	mu = b[0];
	if (as->ai->nexo == 0) {
	    /* just subtract the constant */
	    for (i=0; i<as->n; i++) {
		as->y[i] = as->y0[i];
		if (!isnan(as->y0[i])) {
		    as->y[i] -= mu;
		}
	    }
	}
	b++;
    }

    if (as->P > 0) {
	write_big_phi_197(b, as);
    } else if (as->p > 0) {
	j = 0;
	for (i=0; i<as->p; i++) {
	    if (AR_included(as->ai, i)) {
		as->phi[i] = b[j++];
	    } else {
		as->phi[i] = 0.0;
	    }
	}
    }
    b += np;

    if (as->Q > 0) {
	write_big_theta_197(b, as);
    } else if (as->q > 0) {
	j = 0;
	for (i=0; i<as->q; i++) {
	    if (MA_included(as->ai, i)) {
		as->theta[i] = b[j++];
	    } else {
		as->theta[i] = 0.0;
	    }
	}
    }
    b += nq;

    if (as->ai->nexo > 0) {
	/* subtract the regression effect */
	double xij;

	for (i=0; i<as->n; i++) {
	    as->y[i] = as->y0[i];
	    if (!isnan(as->y[i])) {
		if (as->ifc) {
		    as->y[i] -= mu;
		}
		for (j=0; j<as->ai->nexo; j++) {
		    xij = gretl_matrix_get(as->X, i, j);
		    as->y[i] -= xij * b[j];
		}
	    }
	}
    }
}

static double as197_loglikelihood (const struct as197_info *as)
{
    /* full ARMA loglikelihood */
    double ll1 = 1.0 + LN_2_PI + log(as->sumsq / as->ok_n);

    return -0.5 * as->ok_n * (ll1 + log(as->fact));
}

static double as197_iteration (const double *b, void *data)
{
    struct as197_info *as = data;
    double crit = NADBL;
    /* number of actually included AR terms */
    int np = as->ai->np + as->P;

    as->ncalls += 1;

    if (as->ma_check) {
	/* check that MA term(s) are within bounds */
	const double *theta = b + as->ifc + np;
	const double *Theta = theta + as->ai->nq;

	if (ma_out_of_bounds(as->ai, theta, Theta)) {
	    return crit;
	}
    }

    as197_fill_arrays(as, b);

    as->ifault = flikam(as->phi, as->plen, as->theta, as->qlen,
			as->y, as->e, as->n, &as->sumsq, &as->fact,
			as->vw, as->vl, as->rp1, as->vk, as->r,
			as->toler);

    if (as->ifault > 0) {
	if (as->ifault == 5) {
	    ; // fputs("flikam: (near) non-stationarity\n", stderr);
	} else {
	    fprintf(stderr, "flikam: ifault = %d\n", as->ifault);
	}
	return NADBL;
    }

    if (isnan(as->sumsq) || isnan(as->fact)) {
	; /* leave crit as NA */
    } else {
	/* The criterion used by Melard may work better than
	   the full loglikelihood in the context of his
	   algorithm? But if we're on the first iteration
	   and the sum of squares is too massive, switch to
	   the loglikelihood.
	*/
	if (!as->use_loglik) {
	    /* Melard's criterion */
	    crit = -as->fact * as->sumsq;
	    if (as->ncalls == 1 && crit < -5000) {
		as->use_loglik = 1;
	    }
	}
	if (as->use_loglik) {
	    as->loglik = crit = as197_loglikelihood(as);
	}
    }

    return crit;
}

static const double *as197_llt_callback (const double *b, int i,
					 void *data)
{
    struct as197_info *as = data;
    int err;

    as197_fill_arrays(as, b);
    err = flikam(as->phi, as->plen, as->theta, as->qlen,
		 as->y, as->e, as->n, &as->sumsq, &as->fact,
		 as->vw, as->vl, as->rp1, as->vk, as->r,
		 as->toler);

    return (err)? NULL : as->e;
}

/* FIXME consolidate this with arma_OPG_vcv? */

static int as197_OPG_vcv (MODEL *pmod,
			  struct as197_info *as,
			  double *b, double s2,
			  int k, int T,
			  PRN *prn)
{
    gretl_matrix *G = NULL;
    gretl_matrix *V = NULL;
    int err = 0;

    G = numerical_score_matrix(b, T, k, as197_llt_callback,
			       as, &err);

    if (!err) {
	V = gretl_matrix_XTX_new(G);
	if (V == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	double rcond = gretl_symmetric_matrix_rcond(V, &err);

	if (!err && rcond < 1.0E-10) {
	    pprintf(prn, "OPG: rcond = %g; will try Hessian\n", rcond);
	    err = 1;
	}
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(V);
    }

    if (!err) {
	gretl_matrix_multiply_by_scalar(V, s2);
	err = gretl_model_write_vcv(pmod, V);
    }

    gretl_matrix_free(G);
    gretl_matrix_free(V);

    return err;
}

static int as197_undo_y_scaling (arma_info *ainfo,
				 gretl_matrix *y,
				 double *b,
				 struct as197_info *as)
{
    double *beta = b + ainfo->ifc + ainfo->np + ainfo->P +
	ainfo->nq + ainfo->Q;
    int i, t, err = 0;

    if (ainfo->ifc) {
	b[0] /= ainfo->yscale;
    }

    for (i=0; i<ainfo->nexo; i++) {
	beta[i] /= ainfo->yscale;
    }

    i = ainfo->t1;
    for (t=0; t<ainfo->fullT; t++) {
	if (!isnan(as->y[t])) {
	    as->y[t] /= ainfo->yscale;
	    if (as->y0 != NULL) {
		as->y0[t] /= ainfo->yscale;
	    }
	}
    }

    as->use_loglik = 1;

    if (na(as197_iteration(b, as))) {
	err = 1;
    }

    return err;
}

static int as197_arma_finish (MODEL *pmod,
			      arma_info *ainfo,
			      const DATASET *dset,
			      struct as197_info *as,
			      double *b, gretlopt opt,
			      PRN *prn)
{
    int i, t, k = ainfo->nc;
    int vcv_err = 0;
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
	s2 += as->e[i] * as->e[i];
	pmod->uhat[t] = as->e[i++];
    }

    s2 /= ainfo->T;
    pmod->sigma = sqrt(s2);
    pmod->lnL = as->loglik;

    /* configure for computing variance matrix */
    as->use_loglik = 1;
    as->ma_check = 0;

    if (arma_use_opg(opt)) {
	vcv_err = as197_OPG_vcv(pmod, as, b, s2, k, ainfo->T, prn);
	if (!vcv_err) {
	    gretl_model_set_vcv_info(pmod, VCV_ML, ML_OP);
	    pmod->opt |= OPT_G;
	}
    } else {
	/* base covariance matrix on Hessian (FIXME perhaps QML);
	   for now we'll not fail entirely if we can't come up
	   with a Hessian-based covariance matrix
	*/
	gretl_matrix *Hinv;

	Hinv = numerical_hessian_inverse(b, ainfo->nc, as197_iteration,
					 as, &vcv_err);
	if (!vcv_err) {
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
	gretl_model_set_int(pmod, "as_algo", 197);
	if (arima_ydiff_only(ainfo)) {
	    pmod->opt |= OPT_Y;
	}
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
    double toler = -1.0;
    int use_loglik = 0;
    int err;

    err = as197_info_init(&as, ainfo, toler, use_loglik);
    if (err) {
	return err;
    }

    b = copyvec(coeff, ainfo->nc);
    if (b == NULL) {
	return E_ALLOC;
    }

    y = form_arma_y_vector(ainfo, &err);
    if (!err) {
	as.y = y->val;
	if (as.ifc || ainfo->nexo > 0) {
	    as.y0 = copyvec(as.y, as.n);
	    if (as.y0 == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (!err && ainfo->nexo > 0) {
	if (ainfo->dX != NULL) {
	    as.X = ainfo->dX;
	} else {
	    as.X = form_arma_X_matrix(ainfo, dset, &err);
	    as.free_X = 1;
	}
    }

    if (!err) {
	/* maximize loglikelihood via BFGS */
	int maxit;
	double toler;

	if (as.n > 2000) {
	    /* try to avoid slowdown on big samples */
	    as.toler = 0.0001;
	    as.use_loglik = 1; /* ? */
	} else if (!as.ifc) {
	    as.use_loglik = 1;
	}

	if (as.q > 0 || as.Q > 0) {
	    as.ma_check = 1; /* ? */
	}

	BFGS_defaults(&maxit, &toler, ARMA);

	err = BFGS_max(b, ainfo->nc, maxit, toler,
		       &ainfo->fncount, &ainfo->grcount,
		       as197_iteration, C_LOGLIK,
		       NULL, &as, NULL, opt, ainfo->prn);
	if (!err) {
	    if (ainfo->yscale != 1.0) {
		/* note: this implies recalculation of loglik */
		as197_undo_y_scaling(ainfo, y, b, &as);
	    } else if (!as.use_loglik) {
		/* we haven't already computed this */
		as.loglik = as197_loglikelihood(&as);
	    }
	    gretl_model_set_int(pmod, "fncount", ainfo->fncount);
	    gretl_model_set_int(pmod, "grcount", ainfo->grcount);
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

/* As of 2018-03-17, the AS197 implementation for gretl
   can't handle missing values within the sample range;
   all other "special cases" should be OK.
*/

static int as197_ok (arma_info *ainfo)
{
    if (arma_missvals(ainfo)) {
	return 0;
    } else {
	return 1;
    }
}
