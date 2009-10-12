#include "usermat.h"

#define MAX_ARMA_ORDER 128
#define MAX_ARIMA_DIFF 2

struct arma_info {
    int yno;         /* ID of dependent variable */
    ArmaFlags flags; /* specification flags */
    int ifc;         /* specification includes a constant? */
    int p;           /* max non-seasonal AR order */
    int d;           /* non-seasonal difference */
    int q;           /* max non-seasonal MA order */
    int P;           /* seasonal AR order */
    int D;           /* seasonal difference */
    int Q;           /* seasonal MA order */
    char *pmask;     /* specific AR lags included */
    char *qmask;     /* specific MA lags included */
    int np;          /* total non-seasonal AR lags */
    int nq;          /* total non-seasonal MA lags */
    int maxlag;      /* longest lag in model */
    int nexo;        /* number of other regressors (ARMAX) */
    int nc;          /* total number of coefficients */
    int t1;          /* starting observation */
    int t2;          /* ending observation */
    int pd;          /* periodicity of data */
    int T;           /* full length of data series */
    int *xlist;      /* list of armax regressors */
    double *dy;      /* differenced dependent variable */
    double dyscale;  /* scale factor for dy */
    gretl_matrix *dX;   /* differenced regressors (if applicable) */
    const char *pqspec; /* input string with specific AR, MA lags */
};

#define arma_has_seasonal(a)   ((a)->flags & ARMA_SEAS)
#define arma_is_arima(a)       ((a)->flags & ARMA_DSPEC)
#define arma_by_x12a(a)        ((a)->flags & ARMA_X12A)
#define arma_exact_ml(a)       ((a)->flags & ARMA_EXACT)
#define arma_using_vech(a)     ((a)->flags & ARMA_VECH)
#define arma_least_squares(a)  ((a)->flags & ARMA_LS)

#define set_arma_has_seasonal(a)  ((a)->flags |= ARMA_SEAS)
#define set_arma_is_arima(a)      ((a)->flags |= ARMA_DSPEC)
#define unset_arma_is_arima(a)    ((a)->flags &= ~ARMA_DSPEC)
#define set_arma_use_vech(a)      ((a)->flags |= ARMA_VECH)
#define set_arma_least_squares(a) ((a)->flags |= ARMA_LS)

#define AR_included(a,i) (a->pmask == NULL || a->pmask[i] == '1')
#define MA_included(a,i) (a->qmask == NULL || a->qmask[i] == '1')

static void 
arma_info_init (struct arma_info *ainfo, char flags, 
		const char *pqspec, const DATAINFO *pdinfo)
{
    ainfo->yno = 0;
    ainfo->flags = flags;

    ainfo->p = 0;
    ainfo->d = 0;
    ainfo->q = 0;
    ainfo->P = 0;
    ainfo->D = 0;
    ainfo->Q = 0; 

    ainfo->pmask = NULL;
    ainfo->qmask = NULL;
    
    ainfo->np = 0;
    ainfo->nq = 0;

    ainfo->maxlag = 0;
    ainfo->ifc = 0;
    ainfo->nexo = 0;
    ainfo->nc = 0;

    ainfo->t1 = pdinfo->t1;
    ainfo->t2 = pdinfo->t2;
    ainfo->pd = pdinfo->pd;
    ainfo->T = pdinfo->n;

    ainfo->xlist = NULL;
    ainfo->dy = NULL;
    ainfo->dyscale = 1.0;
    ainfo->dX = NULL;
    ainfo->pqspec = pqspec;
}

static void arma_info_cleanup (struct arma_info *ainfo)
{
    free(ainfo->pmask);
    free(ainfo->qmask);
    free(ainfo->xlist);
    free(ainfo->dy);
    gretl_matrix_free(ainfo->dX);
}

enum {
    AR_MASK,
    MA_MASK
};

/* Create a mask for skipping certain intermediate lags, 
   AR or MA.  This function also sets ainfo->np and ainfo->nq,
   which record the actual number of non-seasonal AR and MA
   lags used.
*/

static char *mask_from_vec (const gretl_vector *v, 
			    struct arma_info *ainfo,
			    int m, int *err)
{
    int vlen = gretl_vector_get_length(v);
    int mlen = (m == AR_MASK)? ainfo->p : ainfo->q;
    int nv = 0, nmax = 0;
    char *mask;
    int i, k;

    mask = malloc(mlen + 1);
    if (mask == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<mlen; i++) {
	mask[i] = '0';
    }
    mask[mlen] = '\0';

    for (i=0; i<vlen; i++) {
	k = v->val[i] - 1; /* convert to zero-based */
	if (k >= 0 && k < mlen) {
	    mask[k] = '1'; 
	    nv++;
	    if (k + 1 > nmax) {
		nmax = k + 1;
	    }
	}
    }

    if (m == AR_MASK) {
	ainfo->p = nmax;
	ainfo->np = nv;
    } else {
	ainfo->q = nmax;
	ainfo->nq = nv;
    }

    if (nv == 0) {
	free(mask);
	mask = NULL;
    }

    return mask;
}

/* construct vector of specific lags from string spec */

static gretl_matrix *
matrix_from_spec (const char *s, int *tmp, int *err)
{
    gretl_matrix *m = NULL;
    const char *p = s;
    char *rem;
    int i, k, n = 0;

    while (*s != '\0' && *s != '}' && n < 20) {
	strtol(s, &rem, 10);
	n++;
	s = rem;
    }

    m = gretl_vector_alloc(n);
    if (m == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    s = p;
    i = 0;
    while (*s != '\0' && *s != '}' && i < 20) {
	k = (int) strtol(s, &rem, 10);
	m->val[i++] = k;
	s = rem;
    }

    *tmp = 1;

    return m;
}

static gretl_matrix *get_arma_pq_vec (struct arma_info *ainfo,
				      int t, int *tmp, int *err)
{
    gretl_matrix *m = NULL;
    const char *test = (t == AR_MASK)? "p=" : "q=";
    const char *s = strstr(ainfo->pqspec, test);

    *tmp = 0;

    if (s != NULL) {
	s += 2;
	if (*s == '{') {
	    m = matrix_from_spec(s + 1, tmp, err);
	} else {
	    char *p, mname[VNAMELEN];

	    *mname = '\0';
	    strncat(mname, s, VNAMELEN - 1);
	    p = strchr(mname, ',');
	    if (p != NULL) {
		*p = '\0';
	    }
	    m = get_matrix_by_name(mname);
	    if (m == NULL) {
		*err = E_UNKVAR;
	    }
	} 
    }

    return m;
}

static int 
arma_make_masks (struct arma_info *ainfo, int *list)
{
    gretl_matrix *m;
    int tmp, err = 0;

    if (ainfo->p > 0) {
	ainfo->np = ainfo->p;
	if (ainfo->pqspec != NULL && *ainfo->pqspec != '\0') {
	    m = get_arma_pq_vec(ainfo, AR_MASK, &tmp, &err);
	    if (m != NULL) {
		ainfo->pmask = mask_from_vec(m, ainfo, AR_MASK, &err);
		if (tmp) {
		    gretl_matrix_free(m);
		}
	    }
	}
    }

    if (ainfo->q > 0 && !err) {
	ainfo->nq = ainfo->q;
	if (ainfo->pqspec != NULL && *ainfo->pqspec != '\0') {
	    m = get_arma_pq_vec(ainfo, MA_MASK, &tmp, &err);
	    if (m != NULL) {
		ainfo->qmask = mask_from_vec(m, ainfo, MA_MASK, &err);
		if (tmp) {
		    gretl_matrix_free(m);
		}
	    }
	}
    }

    return err;
}

static int arma_list_y_position (struct arma_info *ainfo)
{
    int ypos;

    if (arma_is_arima(ainfo)) {
	ypos = (arma_has_seasonal(ainfo))? 9 : 5;
    } else {
	ypos = (arma_has_seasonal(ainfo))? 7 : 4;
    }

    return ypos;
}

#define INT_DEBUG 0

static int arima_integrate (double *dx, const double *x,
			    int t1, int t2, int d, int D, int s)
{
    double *ix;
    int t;

#if INT_DEBUG
    fprintf(stderr, "arima_integrate: t1=%d, t2=%d, d=%d, D=%d, s=%d\n",
	    t1, t2, d, D, s);
#endif

    ix = malloc((t2 + 1) * sizeof *ix);
    if (ix == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<t1; t++) {
	ix[t] = 0.0;
    }

    for (t=t1; t<=t2; t++) {
	ix[t] = dx[t];
	if (d > 0) {
	    ix[t] += x[t-1];
	} 
	if (d == 2) {
	    ix[t] += x[t-1] - x[t-2];
	}
	if (D > 0) {
	    ix[t] += x[t-s];
	    if (d > 0) {
		ix[t] -= x[t-s-1];
	    }
	    if (d == 2) {
		ix[t] -= x[t-s-1] - x[t-s-2];
	    }	    
	} 
	if (D == 2) {
	    ix[t] += x[t-s] - x[t-2*s];
	    if (d > 0) {
		ix[t] -= x[t-s-1] - x[t-2*s-1];
	    }
	    if (d == 2) {
		ix[t] -= x[t-s-1] - x[t-2*s-1];
		ix[t] += x[t-s-2] - x[t-2*s-2];
	    }
	}
    }

#if INT_DEBUG
    for (t=0; t<=t2; t++) {
	fprintf(stderr, "%2d: %12.5g %12.5g %12.5g\n",
		t, x[t], dx[t], ix[t]);
    }
#endif

    /* transcribe integrated result back into "dx" */
    for (t=0; t<=t2; t++) {
	if (t < t1) {
	    dx[t] = NADBL;
	} else {
	    dx[t] = ix[t];
	}
    }

    free(ix);

    return 0;
}

static void ainfo_data_to_model (struct arma_info *ainfo, MODEL *pmod)
{
    pmod->ifc = ainfo->ifc;
    pmod->dfn = ainfo->nc - pmod->ifc;
    pmod->dfd = pmod->nobs - pmod->dfn;
    pmod->ncoeff = ainfo->nc;

    if (arma_has_seasonal(ainfo)) {
	gretl_model_set_int(pmod, "arma_P", ainfo->P);
	gretl_model_set_int(pmod, "arma_Q", ainfo->Q);
	gretl_model_set_int(pmod, "arma_pd", ainfo->pd);	
    }

    if (ainfo->d > 0 || ainfo->D > 0) {
	gretl_model_set_int(pmod, "arima_d", ainfo->d);
	gretl_model_set_int(pmod, "arima_D", ainfo->D);
    }

    if (ainfo->nexo > 0) {
	gretl_model_set_int(pmod, "armax", 1);
    }

    if (ainfo->pmask != NULL) {
	gretl_model_set_string_as_data(pmod, "pmask", 
				       gretl_strdup(ainfo->pmask));
    }

    if (ainfo->qmask != NULL) {
	gretl_model_set_string_as_data(pmod, "qmask", 
				       gretl_strdup(ainfo->qmask));
    }
}

/* write the various statistics from ARMA estimation into
   a gretl MODEL struct */

static void write_arma_model_stats (MODEL *pmod, const int *list, 
				    struct arma_info *ainfo,
				    const double **Z, 
				    const DATAINFO *pdinfo)
{
    const double *y = NULL;
    double mean_error;
    int do_criteria = 1;
    int t;

    pmod->ci = ARMA;

    ainfo_data_to_model(ainfo, pmod);

    free(pmod->list);
    pmod->list = gretl_list_copy(list);

    if (arma_is_arima(ainfo)) {
	y = ainfo->dy;
    } else {
	y = Z[ainfo->yno];
    }

    if (!arma_least_squares(ainfo)) {
	pmod->ybar = gretl_mean(pmod->t1, pmod->t2, y);
	pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, y);
    }

    mean_error = pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(y[t]) && !na(pmod->uhat[t])) {
	    pmod->yhat[t] = y[t] - pmod->uhat[t];
	    pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	    mean_error += pmod->uhat[t];
	}
    }

    if (arma_is_arima(ainfo)) {
	arima_integrate(pmod->yhat, Z[ainfo->yno], pmod->t1, pmod->t2, 
			ainfo->d, ainfo->D, ainfo->pd);
    }

    mean_error /= pmod->nobs;
    gretl_model_set_double(pmod, "mean_error", mean_error);

    if (na(pmod->sigma)) {
	/* in X12A or native exact cases this is already done */
	pmod->sigma = sqrt(pmod->ess / pmod->nobs);
    } 

    pmod->rsq = pmod->adjrsq = pmod->fstt = pmod->chisq = NADBL;
    pmod->tss = NADBL;

    if (arma_least_squares(ainfo)) {
	/* not applicable */
	do_criteria = 0;
    } else if (arma_by_x12a(ainfo) && !na(pmod->criterion[C_AIC])) {
	/* already given by x12a */
	do_criteria = 0;
    }

    if (do_criteria) {
	mle_criteria(pmod, 1);
    }

    gretl_model_add_arma_varnames(pmod, pdinfo, ainfo->yno,
				  ainfo->p, ainfo->q, 
				  ainfo->pmask, ainfo->qmask,
				  ainfo->P, ainfo->Q,
				  ainfo->nexo);
}

static void calc_max_lag (struct arma_info *ainfo)
{
    int pmax = ainfo->p;
    int dmax = ainfo->d;

    if (arma_has_seasonal(ainfo)) {
	pmax += ainfo->P * ainfo->pd;
	dmax += ainfo->D * ainfo->pd;
    }

    ainfo->maxlag = pmax + dmax;

#if ARMA_DEBUG
    fprintf(stderr, "calc_max_lag: ainfo->maxlag = %d\n", ainfo->maxlag);
#endif
}

#define SAMPLE_DEBUG 1

static int 
arma_adjust_sample (const DATAINFO *pdinfo, const double **Z, const int *list,
		    struct arma_info *ainfo)
{
    int vstart = arma_list_y_position(ainfo);
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int an, i, v, t, t1min;
    int anymiss;

#if SAMPLE_DEBUG
    fprintf(stderr, "arma_adjust_sample: at start, t1=%d, t2=%d, maxlag = %d\n",
	    t1, t2, ainfo->maxlag);
#endif

    /* determine the starting point for valid data, t1min */
    t1min = 0;
    for (t=0; t<=pdinfo->t2; t++) {
	anymiss = 0;
	for (i=vstart; i<=list[0]; i++) {
	    if (na(Z[list[i]][t])) {
		anymiss = 1;
		t1min++;
		break;
	    }
	}
	if (!anymiss) {
	    break;
	}
    }

#if SAMPLE_DEBUG
    fprintf(stderr, " phase 1: t1min = %d\n", t1min);
#endif

    /* if the notional starting point, t1, is before the start of
       valid data, t1min, advance t1 */
    if (t1 < t1min) {
	t1 = t1min;
    }
    
    if (!arma_exact_ml(ainfo)) {
	/* conditional ML: ensure that the sample start allows for
	   the required lags of y */
	int t0;

	if (t1 < ainfo->maxlag) {
	    t1 = ainfo->maxlag;
	}
	t0 = t1 - ainfo->maxlag;
	t1min = t1;
	v = list[vstart];
	for (t=t0; t<t1min; t++) {
	    if (na(Z[v][t])) {
		t1 = t + ainfo->maxlag + 1;
	    }
	}
    }

    /* trim any missing obs from the end of the specified sample
       range 
    */
    for (t=pdinfo->t2; t>=t1; t--) {
	anymiss = 0;
	for (i=vstart; i<=list[0]; i++) {
	    if (na(Z[list[i]][t])) {
		anymiss = 1;
		t2--;
		break;
	    }
	}
	if (!anymiss) {
	    break;
	}
    }

    /* check for missing obs within the sample range */
    for (t=t1; t<t2; t++) {
	for (i=vstart; i<=list[0]; i++) {
	    v = list[i];
	    if (na(Z[v][t])) {
		gretl_errmsg_sprintf(_("Missing value encountered for "
				       "variable %d, obs %d"), v, t + 1);
		return 1;
	    }
	}
    }

    an = t2 - t1 + 1;
    if (an <= ainfo->nc) {
	/* insufficient observations */
	return E_DF; 
    }

#if SAMPLE_DEBUG
    fprintf(stderr, "arma_adjust_sample: at end, t1=%d, t2=%d\n",
	    t1, t2);
#endif

    ainfo->t1 = t1;
    ainfo->t2 = t2;

    return 0;
}

#define ARIMA_DEBUG 0

/* remove the intercept from list of regressors */

static int arma_remove_const (int *list, int seasonal, int diffs,
			      const double **Z, const DATAINFO *pdinfo)
{
    int xstart, ret = 0;
    int i, j;

    if (diffs) {
	xstart = (seasonal)? 10 : 6;
    } else {
	xstart = (seasonal)? 8 : 5;
    }

    for (i=xstart; i<=list[0]; i++) {
	if (list[i] == 0 || true_const(list[i], Z, pdinfo)) {
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

static int check_arma_sep (int *list, int sep1, struct arma_info *ainfo)
{
    int sep2 = (sep1 == 3)? 6 : 8;
    int i, err = 0;

    for (i=sep1+1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    if (i == sep2) {
		/* there's a second list separator in the right place:
		   we've got a seasonal specification */
		set_arma_has_seasonal(ainfo);
	    } else {
		err = 1;
	    }
	}
    }

    if (!err && sep1 == 4) {
	/* check for apparent but not "real" arima spec */
	if (arma_has_seasonal(ainfo)) {
	    if (list[2] == 0 && list[6] == 0) {
		gretl_list_delete_at_pos(list, 2);
		gretl_list_delete_at_pos(list, 5);
		unset_arma_is_arima(ainfo);
	    }
	} else {
	    if (list[2] == 0) {
		gretl_list_delete_at_pos(list, 2);
		unset_arma_is_arima(ainfo);
	    }
	}
    }

#if ARIMA_DEBUG
    fprintf(stderr, "check_arma_sep: returning %d\n", err);
#endif

    return err;
}

/* add list of regressors in the ar(i)max case */

static int ainfo_add_xlist (struct arma_info *ainfo, const int *list,
			    int ypos)
{
    ainfo->xlist = gretl_list_new(ainfo->nexo);

    if (ainfo->xlist == NULL) {
	return E_ALLOC;
    } else {
	int i;

	for (i=1; i<=ainfo->xlist[0]; i++) {
	    ainfo->xlist[i] = list[ypos + i];
	}
	
	return 0;
    }
}

#define count_arma_coeffs(a) (a->ifc + a->np + a->nq + a->P + a->Q + a->nexo)

static int check_arma_list (int *list, gretlopt opt, 
			    const double **Z, const DATAINFO *pdinfo,
			    struct arma_info *ainfo)
{
    int armax = 0;
    int hadconst = 0;
    int err = 0;

    if (arma_has_seasonal(ainfo)) {
	armax = (list[0] > 7);
    } else {
	armax = (list[0] > 4);
    }

    if (list[1] < 0 || list[1] > MAX_ARMA_ORDER) {
	err = 1;
    } else if (list[2] < 0 || list[2] > MAX_ARMA_ORDER) {
	err = 1;
    } 

    if (!err) {
	ainfo->p = list[1];
	ainfo->q = list[2];
    }

    if (!err && arma_has_seasonal(ainfo)) {
	if (list[0] < 7) {
	    err = 1;
	} else if (list[4] < 0 || list[4] > MAX_ARMA_ORDER) {
	    err = 1;
	} else if (list[5] < 0 || list[5] > MAX_ARMA_ORDER) {
	    err = 1;
	} 
    }

    if (!err && arma_has_seasonal(ainfo)) {
	ainfo->P = list[4];
	ainfo->Q = list[5];
    }

    /* now that we have p and q we can check for masked lags */

    if (!err) {
	err = arma_make_masks(ainfo, list);
    }

    /* If there's an explicit constant in the list here, we'll remove
       it, since it is added implicitly later.  But if we're supplied
       with OPT_N (meaning: no intercept) we'll flag this by
       setting ifc = 0.  Also, if the user gave an armax list
       (specifying regressors) we'll respect the absence of a constant
       from that list by setting ifc = 0.
    */

    if (!err) {
	if (armax) {
	    hadconst = arma_remove_const(list, arma_has_seasonal(ainfo),
					 0, Z, pdinfo);
	}
	if ((opt & OPT_N) || (armax && !hadconst)) {
	    ;
	} else {
	    ainfo->ifc = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("Error in arma command"));
    } else {
	int ypos = arma_has_seasonal(ainfo) ? 7 : 4;

	ainfo->nexo = list[0] - ypos;
	ainfo->nc = count_arma_coeffs(ainfo);
	ainfo->yno = list[ypos];
	if (ainfo->nexo > 0) {
	    err = ainfo_add_xlist(ainfo, list, ypos);
	}
    }

    return err;
}

static int check_arima_list (int *list, gretlopt opt, 
			     const double **Z, const DATAINFO *pdinfo,
			     struct arma_info *ainfo)
{
    int ypos, armax = 0;
    int hadconst = 0;
    int err = 0;

#if ARIMA_DEBUG
    printlist(list, "check_arima_list");
#endif

    ypos = arma_has_seasonal(ainfo)? 9 : 5;
    armax = list[0] > ypos;

    if (list[1] < 0 || list[1] > MAX_ARMA_ORDER) {
	err = 1;
    } else if (list[2] < 0 || list[2] > MAX_ARIMA_DIFF) {
	err = 1;
    } else if (list[3] < 0 || list[3] > MAX_ARMA_ORDER) {
	err = 1;
    } 

    if (!err) {
	ainfo->p = list[1];
	ainfo->d = list[2];
	ainfo->q = list[3];
    }

    if (!err && arma_has_seasonal(ainfo)) {
	if (list[0] < 9) {
	    err = 1;
	} else if (list[5] < 0 || list[5] > MAX_ARMA_ORDER) {
	    err = 1;
	} else if (list[6] < 0 || list[6] > MAX_ARIMA_DIFF) {
	    err = 1;
	} else if (list[7] < 0 || list[7] > MAX_ARMA_ORDER) {
	    err = 1;
	} 
    }

    if (!err && arma_has_seasonal(ainfo)) {
	ainfo->P = list[5];
	ainfo->D = list[6];
	ainfo->Q = list[7];
    }

    /* now that we have p and q we can check for masked lags */

    if (!err) {
	err = arma_make_masks(ainfo, list);
    }

    /* If there's an explicit constant in the list here, we'll remove
       it, since it is added implicitly later.  But if we're supplied
       with OPT_N (meaning: no intercept) we'll flag this by
       setting ifc = 0.  Also, if the user gave an armax list
       (specifying regressors) we'll respect the absence of a constant
       from that list by setting ifc = 0.
    */

    if (!err) {
	if (armax) {
	    hadconst = arma_remove_const(list, arma_has_seasonal(ainfo),
					 1, Z, pdinfo);
	}
	if ((opt & OPT_N) || (armax && !hadconst)) {
	    ;
	} else {
	    ainfo->ifc = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("Error in arma command"));
    } else {
	ainfo->nexo = list[0] - ypos;
	ainfo->nc = count_arma_coeffs(ainfo);
	ainfo->yno = list[ypos];
	if (ainfo->nexo > 0) {
	    err = ainfo_add_xlist(ainfo, list, ypos);
	}
    }

    return err;
}

static int arma_check_list (int *list, gretlopt opt,
			    const double **Z, const DATAINFO *pdinfo,
			    struct arma_info *ainfo)
{
    int sep1 = gretl_list_separator_position(list);
    int err = 0;

#if ARIMA_DEBUG
    fprintf(stderr, "arma_check_list: sep1 = %d\n", sep1);
    printlist(list, "incoming list");
#endif

    if (sep1 == 3) {
	if (list[0] < 4) {
	    err = E_PARSE;
	}
    } else if (sep1 == 4) {
	if (list[0] < 5) {
	    err = E_PARSE;
	} else {
	    set_arma_is_arima(ainfo);
	}
    } else {
	err = E_PARSE;
    }

    if (!err) {
	err = check_arma_sep(list, sep1, ainfo);
    }

    if (!err) {
	if (arma_is_arima(ainfo)) {
	    /* check for arima spec */
	    err = check_arima_list(list, opt, Z, pdinfo, ainfo);
	} else {	    
	    /* check for simple arma spec */
	    err = check_arma_list(list, opt, Z, pdinfo, ainfo);
	} 
    }

    /* catch null model */
    if (ainfo->nc == 0) {
	err = E_ARGS;
    }

#if ARIMA_DEBUG
    printlist(list, "ar(i)ma list after checking");
    fprintf(stderr, "err = %d\n", err);
#endif

    return err;
}

static void real_difference_series (double *dx, const double *x,
				    struct arma_info *ainfo,
				    int t1)
{
    int t, s = ainfo->pd;
    
    for (t=t1; t<ainfo->T; t++) {
	dx[t] = x[t];
	if (ainfo->d > 0) {
	    dx[t] -= x[t-1];
	} 
	if (ainfo->d == 2) {
	    dx[t] -= x[t-1] - x[t-2];
	}
	if (ainfo->D > 0) {
	    dx[t] -= x[t-s];
	    if (ainfo->d > 0) {
		dx[t] += x[t-s-1];
	    }
	    if (ainfo->d == 2) {
		dx[t] += x[t-s-1] - x[t-s-2];
	    }	    
	} 
	if (ainfo->D == 2) {
	    dx[t] -= x[t-s] - x[t-2*s];
	    if (ainfo->d > 0) {
		dx[t] += x[t-s-1] - x[t-2*s-1];
	    }
	    if (ainfo->d == 2) {
		dx[t] += x[t-s-1] - x[t-2*s-1];
		dx[t] -= x[t-s-2] - x[t-2*s-2];
	    }
	}
    }
}

static int
arima_difference (struct arma_info *ainfo, const double **Z)
{
    const double *y = Z[ainfo->yno];
    double *dy;
    int t, t1 = 0;

#if ARMA_DEBUG
    fprintf(stderr, "doing arima_difference: d = %d, D = %d\n",
	    ainfo->d, ainfo->D);
#endif

    dy = malloc(ainfo->T * sizeof *dy);
    if (dy == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<ainfo->T; t++) {
	if (na(y[t])) {
	    t1++;
	} else {
	    break;
	}
    }

    t1 += ainfo->d + ainfo->D * ainfo->pd;

    for (t=0; t<t1; t++) {
	dy[t] = NADBL;
    }

    real_difference_series(dy, y, ainfo, t1);

#if ARMA_DEBUG > 1
    for (t=0; t<ainfo->T; t++) {
	fprintf(stderr, "dy[%d] = % 12.5f\n", t, dy[t]);
    }
#endif    

    ainfo->dy = dy;

#if 0
    if (ainfo->xlist != NULL) {
	int i, vi, T = ainfo->T - t1;
	double *val;

	ainfo->dX = gretl_matrix_alloc(T, ainfo->nexo);
	if (ainfo->dX != NULL) {
	    val = ainfo->dX;
	    for (i=0; i<ainfo->nexo; i++) {
		vi = ainfo->xlist[i+1];
		real_difference_series(val, Z[vi], ainfo, t1);
		val += T;
	    }
	}
    }
#endif

    return 0;
}

