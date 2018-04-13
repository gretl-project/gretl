/* shared between AS 154 and AS 197 */

struct as_info {
    int algo;
    int p;
    int P;
    int q;
    int Q;
    int pd;
    int plen;
    int qlen;
    int r;
    int rp1;       /* AS 197 */
    int np, nrbar; /* AS 154 */
    int ifault;
    int n;
    int ok_n;
    int ifc;
    /* AR and MA coeffs */
    double *phi, *theta;
     /* dependent var and forecast errors */
    double *y, *y0, *e;
    /* AS 197 workspace */
    double *vw, *vl, *vk;
    /* AS 154 workspace */
    double *A, *P0, *V;
    double *thetab;
    double *xnext, *xrow, *rbar;
    double *evec;
    /* components of likelihood */
    double sumsq, fact, sumlog;
    double toler;  /* tolerance for switching to fast iterations */
    double loglik;
    int ma_check;
    int iupd; /* specific to AS 154 */
    int ncalls;
    int use_loglik;
    arma_info *ai;
    gretl_matrix *X;
    int free_X;
};

static int as_154_alloc (struct as_info *as)
{
    int err = 0;

    /* unused pointers specific to AS 197 */
    as->vw = as->vl = as->vk = NULL;
    
    as->phi =   malloc(as->r * sizeof *as->phi);
    as->theta = malloc(as->r * sizeof *as->theta);
    as->A =     malloc(as->r * sizeof *as->A);
    as->P0 =    malloc(as->np * sizeof *as->P0);
    as->V =     malloc(as->np * sizeof *as->V);
    as->e =     malloc(as->n * sizeof *as->e);
    as->evec =  malloc(as->r * sizeof *as->evec);

    if (as->phi == NULL || as->theta == NULL || as->A == NULL ||
	as->P0 == NULL || as->V == NULL || as->e == NULL ||
	as->evec == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	int worklen = 3 * as->np + as->nrbar;

	as->thetab = malloc(worklen * sizeof *as->thetab);
	if (as->thetab == NULL) {
	    err = E_ALLOC;
	} else {
	    as->xnext = as->thetab + as->np;
	    as->xrow = as->xnext + as->np;
	    as->rbar = as->xrow + as->np;
	}
    }

    return err;
}

static int as_197_alloc (struct as_info *as)
{
    int err = 0;

    /* unused pointers specific to AS 154 */
    as->A = as->P0 = as->V = as->evec = NULL;
    as->thetab = as->xnext = as->xrow = as->rbar = NULL;
		
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

    return err;
}

static int as_info_init (struct as_info *as,
			 int algo,
			 arma_info *ai,
			 double toler,
			 int use_loglik)
{
    int err = 0;

    as->algo = algo;
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

    if (algo == 154) {
	as->np = as->r * (as->r + 1)/2;
	as->nrbar = as->np * (as->np - 1)/2;
	as->rp1 = 0;
    } else {
	as->rp1 = as->r + 1;
	as->np = as->nrbar = 0;
    }

    as->y = as->y0 = NULL; /* later! */
    as->X = NULL; /* later too */
    as->free_X = 0;

    if (algo == 154) {
	err = as_154_alloc(as);
    } else {
	err = as_197_alloc(as);
    }

    if (!err) {
	as->toler = toler;
	as->loglik = NADBL;
	as->ifault = 0;
	as->ma_check = 0;
	as->iupd = 1; /* AS 154: FIXME AR(1) */
	as->ncalls = 0;
	as->use_loglik = use_loglik;
    }

    return err;
}

static void as_info_free (struct as_info *as)
{
    free(as->phi);
    free(as->theta);
    free(as->e);
    free(as->y0);

    if (as->algo == 154) {
	free(as->A);
	free(as->P0);
	free(as->V);
	free(as->evec);
	free(as->thetab);
    }

    if (as->free_X) {
	gretl_matrix_free(as->X);
    }
}

static void as_write_big_phi (const double *b,
			      struct as_info *as)
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

static void as_write_big_theta (const double *b,
				struct as_info *as)
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

static void as_fill_arrays (struct as_info *as,
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
	as_write_big_phi(b, as);
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
	as_write_big_theta(b, as);
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

/* full ARMA loglikelihood */

static double as_loglikelihood (const struct as_info *as)
{
    double ll1 = 1.0 + LN_2_PI + log(as->sumsq / as->ok_n);

    if (as->algo == 154) {
	return -0.5 * (as->ok_n * ll1 + as->sumlog);
    } else {
	return -0.5 * as->ok_n * (ll1 + log(as->fact));
    }
}

static double as197_iteration (const double *b, void *data)
{
    struct as_info *as = data;
    double crit = NADBL;
    /* number of actually included AR terms */
    int np = as->ai->np + as->P;

    as->ncalls += 1;

    if (as->ma_check) {
	/* check that MA term(s) are within bounds */
	double *theta = (double *) b + as->ifc + np;
	double *Theta = theta + as->ai->nq;

	if (maybe_correct_MA(as->ai, theta, Theta)) {
	    return NADBL;
	}
    }

    as_fill_arrays(as, b);

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
	    as->loglik = crit = as_loglikelihood(as);
	}
    }

    return crit;
}

static double as154_iteration (const double *b, void *data)
{
    struct as_info *as = data;
    double crit = NADBL;
    /* number of actually included AR terms */
    int np = as->ai->np + as->P;
    int nit = 0;

    if (as->ma_check) {
	/* check that MA term(s) are within bounds */
	double *theta = (double *) b + as->ifc + np;
	double *Theta = theta + as->ai->nq;

	if (maybe_correct_MA(as->ai, theta, Theta)) {
	    return NADBL;
	}
    }

    as_fill_arrays(as, b);

    as->ifault = starma(as->plen, as->qlen, as->r, as->np,
			as->phi, as->theta, as->A, as->P0, as->V,
			as->thetab, as->xnext, as->xrow,
			as->rbar, as->nrbar);

    if (as->ifault) {
	fprintf(stderr, "starma: ifault = %d\n", as->ifault);
	return NADBL;
    }

    /* initialization required */
    as->sumlog = as->sumsq = 0;

    karma(as->plen, as->qlen, as->r, as->np,
	  as->phi, as->theta, as->A, as->P0, as->V,
	  as->n, as->y, as->e,
	  &as->sumlog, &as->sumsq, as->iupd,
	  as->toler, as->evec, &nit);

    if (isnan(as->sumlog) || isnan(as->sumsq) || as->sumsq <= 0) {
	; // fprintf(stderr, "karma: got NaNs, nit = %d\n", nit);
    } else {
	/* The criterion used by Gardner at al may work
	   better than the full loglikelihood in the
	   context of their algorithm?
	*/
	if (as->use_loglik) {
	    as->loglik = crit = as_loglikelihood(as);
	} else {
	    /* Gardner et al criterion */
	    crit = -(as->ok_n * log(as->sumsq) + as->sumlog);
	}
    }

    return crit;
}

static const double *as197_llt_callback (const double *b, int i,
					 void *data)
{
    struct as_info *as = data;
    int err;

    as_fill_arrays(as, b);
    err = flikam(as->phi, as->plen, as->theta, as->qlen,
		 as->y, as->e, as->n, &as->sumsq, &as->fact,
		 as->vw, as->vl, as->rp1, as->vk, as->r,
		 as->toler);

    return (err)? NULL : as->e;
}

static const double *as154_llt_callback (const double *b, int i,
					 void *data)
{
    struct as_info *as = data;
    int err = 0, nit = 0;

    as_fill_arrays(as, b);
    as->ifault = starma(as->plen, as->qlen, as->r, as->np,
			as->phi, as->theta, as->A, as->P0, as->V,
			as->thetab, as->xnext, as->xrow,
			as->rbar, as->nrbar);

    as->sumlog = as->sumsq = 0;
    karma(as->plen, as->qlen, as->r, as->np,
	  as->phi, as->theta, as->A, as->P0, as->V,
	  as->n, as->y, as->e,
	  &as->sumlog, &as->sumsq, as->iupd,
	  as->toler, as->evec, &nit);

    if (isnan(as->sumlog) || isnan(as->sumsq) || as->sumsq <= 0) {
	fprintf(stderr, "as154_llt_callback: failed\n");
	err = E_NAN;
    }

    return (err)? NULL : as->e;
}

static int as_undo_y_scaling (arma_info *ainfo,
			      gretl_matrix *y,
			      double *b,
			      struct as_info *as)
{
    double *beta = b + ainfo->ifc + ainfo->np + ainfo->P +
	ainfo->nq + ainfo->Q;
    double lnl;
    int i, t, err = 0;

    if (ainfo->ifc) {
	b[0] /= ainfo->yscale;
    }

    for (i=0; i<ainfo->nexo; i++) {
	beta[i] /= ainfo->yscale;
    }

    for (t=0; t<ainfo->fullT; t++) {
	if (!isnan(as->y[t])) {
	    as->y[t] /= ainfo->yscale;
	    if (as->y0 != NULL) {
		as->y0[t] /= ainfo->yscale;
	    }
	}
    }

    as->use_loglik = 1;

    if (as->algo == 154) {
	lnl = as154_iteration(b, as);
    } else {
	lnl = as197_iteration(b, as);
    }

    if (na(lnl)) {
	err = 1;
    }

    return err;
}

static int as_arma_finish (MODEL *pmod,
			   arma_info *ainfo,
			   const DATASET *dset,
			   struct as_info *as,
			   double *b, gretlopt opt,
			   PRN *prn)
{
    int i, t, k = ainfo->nc;
    int do_opg = arma_use_opg(opt);
    int QML = (opt & OPT_R);
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
	if (isnan(as->e[i])) {
	    pmod->uhat[t] = NADBL;
	} else {
	    s2 += as->e[i] * as->e[i];
	    pmod->uhat[t] = as->e[i];
	}
	i++;
    }

    s2 /= ainfo->T;
    pmod->sigma = sqrt(s2);
    pmod->lnL = as->loglik;

    /* configure for computing variance matrix */
    as->use_loglik = 1;
    as->ma_check = 0;

    if (!do_opg) {
	/* base covariance matrix on Hessian (perhaps QML) */
	gretl_matrix *Hinv;

	if (as->algo == 154) {
	    Hinv = numerical_hessian_inverse(b, ainfo->nc, as154_iteration,
					     as, &vcv_err);
	} else {
	    Hinv = numerical_hessian_inverse(b, ainfo->nc, as197_iteration,
					     as, &vcv_err);
	}
	if (!vcv_err) {
	    if (QML) {
		vcv_err = arma_QML_vcv(pmod, Hinv, as, as->algo, b, s2,
				       k, ainfo->T, prn);
	    } else {
		err = gretl_model_write_vcv(pmod, Hinv);
		if (!err) {
		    gretl_model_set_vcv_info(pmod, VCV_ML, ML_HESSIAN);
		}
	    }
	} else if (!(opt & OPT_H)) {
	    /* fallback when Hessian not explicitly requested */
	    do_opg = 1;
	    gretl_model_set_int(pmod, "hess-error", 1);
	}
	gretl_matrix_free(Hinv);
    }

    if (do_opg) {
	vcv_err = arma_OPG_vcv(pmod, as, as->algo, b, s2, k, ainfo->T, prn);
	if (!vcv_err) {
	    gretl_model_set_vcv_info(pmod, VCV_ML, ML_OP);
	    pmod->opt |= OPT_G;
	}
    }

    if (!err) {
	write_arma_model_stats(pmod, ainfo, dset);
	arma_model_add_roots(pmod, ainfo, b);
	gretl_model_set_int(pmod, "arma_flags", ARMA_EXACT);
	gretl_model_set_int(pmod, "as_algo", as->algo);
	if (arima_ydiff_only(ainfo)) {
	    pmod->opt |= OPT_Y;
	}
    }

    return err;
}

/* As of 2018-04, the AS 197 implementation for gretl
   can't properly handle missing values within the sample
   range; all other "special cases" should be OK.
*/

static int as197_ok (arma_info *ainfo)
{
    return arma_missvals(ainfo) ? 0 : 1;
}

static int as_arma (const double *coeff,
		    const DATASET *dset,
		    arma_info *ainfo,
		    MODEL *pmod,
		    gretlopt opt)
{
    struct as_info as = {0};
    gretl_matrix *y = NULL;
    double *b = NULL;
    double toler = -1.0;
    int use_loglik = 0;
    int algo, err = 0;

    algo = get_optval_int(ARMA, OPT_A, &err);
    if (algo == 0) {
	/* user didn't specify: prefer 197 unless there are
	   missing values to be handled */
	algo = as197_ok(ainfo) ? 197 : 154;
    }
    if (err || (algo != 154 && algo != 197)) {
	return E_BADOPT;
    }

    err = as_info_init(&as, algo, ainfo, toler, use_loglik);
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

	if (as.algo == 197) {
	    if (as.n > 2000) {
		/* try to avoid slowdown on big samples */
		as.toler = 0.0001;
		as.use_loglik = 1; /* ? */
	    } else if (!as.ifc) {
		as.use_loglik = 1;
	    }
	} else {
	    /* AS 154 */
	    as.use_loglik = 1; /* generally better? */
	}

	if (as.q > 0 || as.Q > 0) {
	    as.ma_check = 1;
	}

	BFGS_defaults(&maxit, &toler, ARMA);

	if (as.algo == 154) {
	    err = BFGS_max(b, ainfo->nc, maxit, toler,
			   &ainfo->fncount, &ainfo->grcount,
			   as154_iteration, C_LOGLIK,
			   NULL, &as, NULL, opt, ainfo->prn);
	} else {
	    err = BFGS_max(b, ainfo->nc, maxit, toler,
			   &ainfo->fncount, &ainfo->grcount,
			   as197_iteration, C_LOGLIK,
			   NULL, &as, NULL, opt, ainfo->prn);
	}
	if (!err) {
	    if (ainfo->yscale != 1.0) {
		/* note: this implies recalculation of loglik */
		as_undo_y_scaling(ainfo, y, b, &as);
	    } else if (!as.use_loglik) {
		/* we haven't already computed this */
		as.loglik = as_loglikelihood(&as);
	    }
	    gretl_model_set_int(pmod, "fncount", ainfo->fncount);
	    gretl_model_set_int(pmod, "grcount", ainfo->grcount);
	    err = as_arma_finish(pmod, ainfo, dset, &as, b,
				 opt, ainfo->prn);
	}
    }

    if (err && !pmod->errcode) {
	pmod->errcode = err;
    }

    as_info_free(&as);
    gretl_matrix_free(y);
    free(b);

    return err;
}

