#define MAX_ARMA_ORDER 6

enum {
    ARMA_NATIVE,
    ARMA_X12A
};

struct arma_info {
    int atype;    /* code for native or x12arima */
    int yno;      /* ID of dependent variable */
    int p;        /* non-seasonal AR order */
    int q;        /* non-seasonal MA order */
    int P;        /* seasonal AR order */
    int Q;        /* seasonal MA order */
    int maxlag;   /* longest lag in model */
    int r;        /* number of other regressors (ARMAX) */
    int ifc;      /* 1 for intercept included, otherwise 0 */
    int nc;       /* total number of coefficients */
    int t1;       /* starting observation */
    int t2;       /* ending observation */
    int seasonal; /* 1 if any seasonal terms, otherwise 0 */
    int pd;       /* periodicity of data */
};

/* write the various statistics from ARMA estimation into
   a gretl MODEL struct */

static void write_arma_model_stats (MODEL *pmod, model_info *arma,
				    const int *list, const double *y, 
				    const double *theta, 
				    struct arma_info *ainfo)
{
    double **series = NULL;
    const double *e = NULL;
    double mean_error;
    int i, t;

    if (arma != NULL) {
	series = model_info_get_series(arma);
	e = series[0];
	pmod->lnL = model_info_get_ll(arma);
    }

    pmod->ci = ARMA;
    pmod->ifc = ainfo->ifc;

    pmod->dfn = ainfo->nc - ainfo->ifc;
    pmod->dfd = pmod->nobs - pmod->dfn;
    pmod->ncoeff = ainfo->nc;

    if (theta != NULL) {
	for (i=0; i<pmod->ncoeff; i++) {
	    pmod->coeff[i] = theta[i];
	}
    }

    free(pmod->list);
    pmod->list = gretl_list_copy(list);

    pmod->ybar = gretl_mean(pmod->t1, pmod->t2, y);
    pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, y);

    mean_error = pmod->ess = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (e != NULL) {
	    pmod->uhat[t] = e[t];
	}
	pmod->yhat[t] = y[t] - pmod->uhat[t];
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	mean_error += pmod->uhat[t];
    }

    mean_error /= pmod->nobs;
    gretl_model_set_double(pmod, "mean_error", mean_error);

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    pmod->tss = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->tss += (y[t] - pmod->ybar) * (y[t] - pmod->ybar);
    }

    if (!pmod->ifc) {
	double syh2 = 0.0;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    syh2 += pmod->yhat[t] * pmod->yhat[t];
	}
	pmod->fstt = pmod->dfd * syh2 / (pmod->dfn * pmod->ess);
    } else if (pmod->tss > pmod->ess) {
	pmod->fstt = pmod->dfd * (pmod->tss - pmod->ess) / (pmod->dfn * pmod->ess);
    } else {
	pmod->fstt = NADBL;
    }

    pmod->rsq = pmod->adjrsq = NADBL;

    if (!pmod->ifc) {
	double r2 = gretl_corr_rsq(pmod->t1, pmod->t2, y, pmod->yhat);

	pmod->rsq = r2;
	pmod->adjrsq = 1.0 - ((1.0 - r2) * (pmod->nobs - 1.0) / pmod->dfd);
    } else if (pmod->tss > 0) {
	pmod->rsq = 1.0 - (pmod->ess / pmod->tss);
	if (pmod->dfd > 0) {
	    double den = pmod->tss * pmod->dfd;

	    pmod->adjrsq = 1.0 - (pmod->ess * (pmod->nobs - 1) / den);
	}
    }

    if (arma != NULL) {
	mle_aic_bic(pmod, 1);
    }

    if (ainfo->seasonal) {
	gretl_model_set_int(pmod, "arma_P", ainfo->P);
	gretl_model_set_int(pmod, "arma_Q", ainfo->Q);
	gretl_model_set_int(pmod, "arma_pd", ainfo->pd);	
    }

    if (ainfo->r > 0) {
	gretl_model_set_int(pmod, "armax", 1);
    }
}

static void calc_max_lag (struct arma_info *ainfo)
{
    int pmax = ainfo->p;
    int qmax = ainfo->q;

    if (ainfo->seasonal) {
	pmax += ainfo->P * ainfo->pd;
	qmax += ainfo->Q * ainfo->pd;
    }

    ainfo->maxlag = (pmax > qmax)? pmax : qmax;
}

static int 
arma_adjust_sample (const DATAINFO *pdinfo, const double **Z, const int *list,
		    struct arma_info *ainfo)
{
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int an, i, v, t, t1min = 0;
    int vstart, pmax, anymiss;

    vstart = (ainfo->seasonal)? 7 : 4;

    pmax = ainfo->p;
    if (ainfo->P > 0) {
	pmax += ainfo->P * ainfo->pd;
    }   

    for (t=0; t<=pdinfo->t2; t++) {
	anymiss = 0;
	for (i=vstart; i<=list[0]; i++) {
	    v = list[i];
	    if (na(Z[v][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (anymiss) {
	    t1min++;
        } else {
	    break;
	}
    }

    if (ainfo->atype == ARMA_NATIVE) {
	/* is this right for X-12-ARIMA? */
	t1min += ainfo->maxlag;
    }

    if (t1 < t1min) {
	t1 = t1min;
    }

    for (t=pdinfo->t2; t>=t1; t--) {
	anymiss = 0;
	for (i=vstart; i<=list[0]; i++) {
	    v = list[i];
	    if (na(Z[v][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (anymiss) {
	    t2--;
        } else {
	    break;
	}
    }

    for (t=t1-pmax; t<t2; t++) {
	for (i=vstart; i<=list[0]; i++) {
	    if (t < t1 && i > vstart) {
		continue;
	    }
	    v = list[i];
	    if (na(Z[v][t])) {
		char msg[64];

		sprintf(msg, _("Missing value encountered for "
			       "variable %d, obs %d"), v, t + 1);
		gretl_errmsg_set(msg);
		return 1;
	    }
	}
    }

    an = t2 - t1 + 1;
    if (an <= ainfo->nc) {
	return 1; 
    }

    ainfo->t1 = t1;
    ainfo->t2 = t2;

    return 0;
}

/* remove the intercept from list of regressors */

static int remove_const (int *list, int seasonal)
{
    int xstart, ret = 0;
    int i, j;

    xstart = (seasonal)? 8 : 5;

    for (i=xstart; i<=list[0]; i++) {
	if (list[i] == 0) {
	    for (j=i; j<list[0]; j++) {
		list[j] = list[j+1];
	    }
	    list[0] -= 1;
	    ret = 1;
	    break;
	}
    }

    return ret;
}

#define has_seasonals(l) (l[0] > 5 && l[3] == LISTSEP && l[6] == LISTSEP)

static int 
check_arma_list (int *list, gretlopt opt, struct arma_info *ainfo)
{
    int armax = 0;
    int hadconst = 0;
    int err = 0;

    /* FIXME need more checks for seasonal case, e.g. non-
       seasonal order must be less than pd */

    ainfo->seasonal = has_seasonals(list);
    ainfo->p = ainfo->q = 0;
    ainfo->P = ainfo->Q = 0;
    ainfo->r = ainfo->nc = 0;

    if (ainfo->seasonal) {
	armax = (list[0] > 7);
    } else {
	armax = (list[0] > 4);
    }

    if (list[0] < 4) {
	err = 1;
    } else if (list[1] < 0 || list[1] > MAX_ARMA_ORDER) {
	err = 1;
    } else if (list[2] < 0 || list[2] > MAX_ARMA_ORDER) {
	err = 1;
    } else if (list[1] + list[2] == 0 && !ainfo->seasonal) {
	err = 1;
    }

    if (!err) {
	ainfo->p = list[1];
	ainfo->q = list[2];
    }

    if (!err && ainfo->seasonal) {
	if (list[0] < 7) {
	    err = 1;
	} else if (list[4] < 0 || list[4] > MAX_ARMA_ORDER) {
	    err = 1;
	} else if (list[5] < 0 || list[5] > MAX_ARMA_ORDER) {
	    err = 1;
	} else if (list[4] + list[5] == 0) {
	    err = 1;
	}
    }

    if (!err && ainfo->seasonal) {
	ainfo->P = list[4];
	ainfo->Q = list[5];
    }

    /* If there's an explicit constant in the list here, we'll remove
       it, since it is added implicitly later.  But if we're supplied
       with OPT_S (meaning: suppress the intercept) we'll flag this by
       setting ifc = 0.  Also, if the user gave an armax list
       (specifying regressors) we'll respect the absence of a constant
       from that list by setting ifc = 0.
    */

    if (!err) {
	if (armax) {
	    hadconst = remove_const(list, ainfo->seasonal);
	}
	if ((opt & OPT_S) || (armax && !hadconst)) {
	    ainfo->ifc = 0;
	} else {
	    ainfo->ifc = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("Error in arma command"));
    } else {
	ainfo->r = list[0] - ((ainfo->seasonal)? 7 : 4);
	ainfo->nc = ainfo->p + ainfo->q + ainfo->P + ainfo->Q
	    + ainfo->r + ainfo->ifc;
    }

    return err;
}

