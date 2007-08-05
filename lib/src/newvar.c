#define VDEBUG 0

typedef struct VARspec_ VARspec;

struct VARspec_ {
    int n;       /* number of equations */
    int g;       /* total coefficients per equation */
    int t1;      /* starting obs */
    int t2;      /* ending obs */
    int T;       /* total observations */
    int p;       /* lag order */
    int robust;  /* robust std errors? */
    int qr;      /* use QR decomp for OLS? */
    int *ylist;
    int *xlist;
    int detflags;
    gretl_matrix *Y;
    gretl_matrix *X;
    gretl_matrix *B;
    gretl_matrix *E;
    gretl_matrix *XTX;
};

static void fill_VAR_X (VARspec *v, int p, const double **Z, 
			const DATAINFO *pdinfo);

/* apparatus for selecting the optimal lag length for a VAR */

static int 
alt_VAR_do_lagsel (GRETL_VAR *var, VARspec *vspec, 
		   const double **Z, const DATAINFO *pdinfo, 
		   PRN *prn)
{
    gretl_matrix *crittab = NULL;
    gretl_matrix *lltab = NULL;

    int p = vspec->p;
    int r = p - 1;
    int T = vspec->T;
    int n = vspec->n;

    /* initialize the "best" at the longest lag */
    double best[N_IVALS] = { var->AIC, var->BIC, var->HQC };
    int best_row[N_IVALS] = { r, r, r };
    double crit[N_IVALS];
    double LRtest;
    double ldet = NADBL;
    int cols0;
    int j, m = 0;
    int err = 0;

    if (p < 2) {
	return 0;
    }

    if (var->F != NULL) {
	gretl_matrix_free(var->F);
    }

    var->F = gretl_matrix_alloc(T, n);
    if (var->F == NULL) {
	return E_ALLOC;
    }

    crittab = gretl_matrix_alloc(p, N_IVALS);
    lltab = gretl_matrix_alloc(p, 2);
    if (crittab == NULL || lltab == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* # of cols in X that are not Y lags */
    cols0 = vspec->g - p * n; 

    for (j=1; j<p && !err; j++) {
	int jxcols = cols0 + j * n;

	fill_VAR_X(vspec, j, Z, pdinfo);

	gretl_matrix_reuse(vspec->X, T, jxcols);
	gretl_matrix_reuse(vspec->B, jxcols, n);

	err = gretl_matrix_multi_ols(vspec->Y, vspec->X, vspec->B, 
				     var->F, NULL);

	if (!err) {
	    ldet = gretl_VAR_ldet(var, &err);
	}

	if (!err) {
	    double ll;
	    int q = vspec->g - (n * (p - j));
	    int c, k = n * q;

	    ll = -(n * T / 2.0) * (LN_2_PI + 1) - (T / 2.0) * ldet;
	    crit[0] = (-2.0 * ll + 2.0 * k) / T;               /* AIC */
	    crit[1] = (-2.0 * ll + k * log(T)) / T;            /* BIC */
	    crit[2] = (-2.0 * ll + 2.0 * k * log(log(T))) / T; /* HQC */

	    gretl_matrix_set(lltab, m, 0, ll);
	    if (j == 1) {
		gretl_matrix_set(lltab, m, 1, 0);
	    } else {
		LRtest = 2.0 * (ll - gretl_matrix_get(lltab, m-1, 0));
		gretl_matrix_set(lltab, m, 1, chisq_cdf_comp(LRtest, n * n));
	    }	
	    
	    for (c=0; c<N_IVALS; c++) {
		gretl_matrix_set(crittab, m, c, crit[c]);
		if (crit[c] < best[c]) {
		    best[c] = crit[c];
		    best_row[c] = m;
		}
	    }
	
	    m++;
	}
    }

    if (!err) {
	gretl_matrix_set(lltab, m, 0, var->ll);
	LRtest = 2.0 * (var->ll - gretl_matrix_get(lltab, m - 1, 0));
	gretl_matrix_set(lltab, m, 1, chisq_cdf_comp(LRtest, n * n));
	gretl_matrix_set(crittab, m, 0, var->AIC);
	gretl_matrix_set(crittab, m, 1, var->BIC);
	gretl_matrix_set(crittab, m, 2, var->HQC);
	gretl_VAR_print_lagsel(lltab, crittab, best_row, prn);
    }

    bailout:

    gretl_matrix_free(crittab);
    gretl_matrix_free(lltab);

    gretl_matrix_free(var->F);
    var->F = NULL;

    return err;
}

static int 
gretl_matrix_delete_columns (gretl_matrix *X, int *list)
{
    size_t csz = X->rows * sizeof *X->val;
    void *dest, *src;
    int i, j, n, col;

    for (i=1; i<=list[0]; i++) {
	col = list[i];
	if (col < 0 || col >= X->cols) {
	    return E_NONCONF;
	}
    }

    for (i=1; i<=list[0]; i++) {
	col = list[i];
	dest = X->val + col * X->rows;
	src = X->val + (col + 1) * X->rows;
	n = X->cols - col - 1;
	if (n > 0) {
	    memmove(dest, src, n * csz);
	}
	for (j=i+1; j<=list[0]; j++) {
	    list[j] -= 1;
	}
    }

    X->cols -= list[0];

    return 0;
}

static void varspec_free (VARspec *v)
{
    if (v == NULL) {
	return;
    }

    free(v->ylist);
    free(v->xlist);

    gretl_matrix_free(v->Y);
    gretl_matrix_free(v->X);
    gretl_matrix_free(v->B);
    gretl_matrix_free(v->E);
    gretl_matrix_free(v->XTX);

    free(v);
}

static VARspec *varspec_new (const int *list, int order, 
			     const double **Z, DATAINFO *pdinfo,
			     gretlopt opt, int *err)
{
    VARspec *v = malloc(sizeof *v);
    int i, t, vi;

    if (v == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    v->p = order;

    v->ylist = v->xlist = NULL;
    v->Y = v->X = NULL;
    v->B = v->E = NULL;
    v->XTX = NULL;
    v->g = 0;
    v->detflags = 0;

    v->robust = (opt & OPT_R)? 1 : 0;
    v->qr = get_use_qr();

    if (gretl_list_has_separator(list)) {
	*err = gretl_list_split_on_separator(list, &v->ylist, &v->xlist);
    } else {
	v->ylist = gretl_list_copy(list);
	if (v->ylist == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	gretl_list_purge_const(v->ylist, Z, pdinfo);
	if (v->xlist != NULL) {
	    gretl_list_purge_const(v->xlist, Z, pdinfo);
	}
    }

#if VDEBUG
    printlist(v->ylist, "v->ylist");
    printlist(v->xlist, "v->xlist");
#endif

    if (!*err) {
	v->n = v->ylist[0];
	v->t1 = pdinfo->t1;
	v->t2 = pdinfo->t2;
    }

    /* advance t1 if needed */

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	int miss = 0, s = t - v->p;

	for (i=1; i<=v->ylist[0] && !miss; i++) {
	    vi = v->ylist[i];
	    if (na(Z[vi][t]) || s < 0 || na(Z[vi][s])) {
		v->t1 += 1;
		miss = 1;
	    }
	}
	if (v->xlist != NULL && !miss) {
	    for (i=1; i<=v->xlist[0] && !miss; i++) {
		if (na(Z[v->xlist[i]][t])) {
		    v->t1 += 1;
		    miss = 1;
		}
	    }
	}
	if (!miss) {
	    break;
	}
    }

    /* retard t2 if needed */

    for (t=pdinfo->t2; t>=v->t1; t--) {
	int miss = 0;

	for (i=1; i<=v->ylist[0] && !miss; i++) {
	    if (na(Z[v->ylist[i]][t])) {
		v->t2 -= 1;
		miss = 1;
	    }
	}
	if (v->xlist != NULL && !miss) {
	    for (i=1; i<=v->xlist[0] && !miss; i++) {
		if (na(Z[v->xlist[i]][t])) {
		    v->t2 -= 1;
		    miss = 1;
		}
	    }
	}
	if (!miss) {
	    break;
	}
    }

    /* reject sample in case of internal missing values */

    for (t=v->t1; t<=v->t2 && !*err; t++) {
	for (i=1; i<=v->ylist[0] && !err; i++) {
	    if (na(Z[v->ylist[i]][t])) {
		*err = E_MISSDATA;
	    }
	}
	if (v->xlist != NULL && !*err) {
	    for (i=1; i<=v->xlist[0] && !err; i++) {
		if (na(Z[v->xlist[i]][t])) {
		    *err = E_MISSDATA;
		}
	    }
	}
    }

    /* account for deterministic terms and check for
       non-negative degrees of freedom */

    if (!*err) {
	v->T = v->t2 - v->t1 + 1;
	v->g = v->p * v->n;

	if (v->xlist != NULL) {
	    v->g += v->xlist[0];
	}

	if (!(opt & OPT_N)) {
	    v->detflags |= DET_CONST;
	    v->g += 1;
	}
	if ((opt & OPT_D) && pdinfo->pd != 1) {
	    v->detflags |= DET_SEAS;
	    v->g += pdinfo->pd - 1;
	}
	if (opt & OPT_T) {
	    v->detflags |= DET_TREND;
	    v->g += 1;
	}
	if (v->T < v->g) {
	    *err = E_DF;
	}
    }

#if VDEBUG
    fprintf(stderr, "varspec: n=%d, p=%d, g=%d, t1=%d, t2=%d\n",
	    v->n, v->p, v->g, v->t1, v->t2);
#endif

    return v;
}

/* get the starting sub-period for a given t, so as to construct
   correctly aligned seasonal dummies
*/

static int startp (int t, const DATAINFO *pdinfo)
{
    int yy, pp = pdinfo->pd;
    int mm = 10;
    double xx;

    while ((pp = pp / 10)) {
	mm *= 10;
    }

    xx = date(t, pdinfo->pd, pdinfo->sd0);
    if (dataset_is_daily(pdinfo)) {
	xx += .1;
    }

    yy = (int) xx;
    pp = (int) (mm * (xx - yy) + 0.5);

    return pp;
}

static void fill_VAR_X (VARspec *v, int p, const double **Z, 
			const DATAINFO *pdinfo)
{
    int i, j, s, t, vi;
    int k = 0; /* X column */

    /* construct the X matrix: const first */

    if (v->detflags & DET_CONST) {
	s = 0;
	for (t=v->t1; t<=v->t2; t++) {
	    gretl_matrix_set(v->X, s++, k, 1.0);
	}
	k++;
    }    

    /* add lagged Ys */

    for (i=0; i<v->n; i++) {
	vi = v->ylist[i+1];
	for (j=1; j<=p; j++) {
	    s = 0;
	    for (t=v->t1; t<=v->t2; t++) {
		gretl_matrix_set(v->X, s++, k, Z[vi][t-j]);
	    }
	    k++;
	}
    }

    /* add any exogenous vars */

    if (v->xlist != NULL) {
	for (i=1; i<=v->xlist[0]; i++) {
	    vi = v->xlist[i];
	    s = 0;
	    for (t=v->t1; t<=v->t2; t++) {
		gretl_matrix_set(v->X, s++, k, Z[vi][t]);
	    }
	    k++;
	}
    }

    /* add other deterministics */

    if (v->detflags & DET_SEAS) {
	int per, per0 = startp(v->t1, pdinfo);

	for (i=0; i<pdinfo->pd-1; i++) {
	    per = per0;
	    s = 0;
	    for (t=v->t1; t<=v->t2; t++) {
		gretl_matrix_set(v->X, s++, k, (per == i+1));
		if (per < pdinfo->pd) {
		    per++;
		} else {
		    per = 1;
		}
	    }
	    k++;
	}	
    }

    if (v->detflags & DET_TREND) {
	s = 0;
	for (t=v->t1; t<=v->t2; t++) {
	    gretl_matrix_set(v->X, s++, k, (double) (t + 1));
	}
	k++;
    }

#if VDEBUG
    gretl_matrix_print(v->X, "X");
#endif
}

static int 
make_VAR_matrices (VARspec *v, const double **Z, const DATAINFO *pdinfo)
{
    int i, s, t, vi;

    v->Y = gretl_matrix_alloc(v->T, v->n);
    v->X = gretl_matrix_alloc(v->T, v->g);
    v->B = gretl_matrix_alloc(v->g, v->n);
    v->E = gretl_matrix_alloc(v->T, v->n);

    if (v->Y == NULL || v->X == NULL || v->B == NULL || 
	v->E == NULL) {
	return E_ALLOC;
    }    

    /* construct the Y matrix */

    for (i=0; i<v->n; i++) {
	vi = v->ylist[i+1];
	s = 0;
	for (t=v->t1; t<=v->t2; t++) {
	    gretl_matrix_set(v->Y, s++, i, Z[vi][t]);
	}
    }

#if VDEBUG
    gretl_matrix_print(v->Y, "Y");
#endif

    fill_VAR_X(v, v->p, Z, pdinfo);

    return 0;
}

static void 
set_VAR_param_names (VARspec *v, char **params, const DATAINFO *pdinfo)
{
    char lagstr[8];
    int i, j, n, k = 0;

    if (v->detflags & DET_CONST) {
	strcpy(params[k++], pdinfo->varname[0]);
    }     

    for (i=1; i<=v->ylist[0]; i++) {
	for (j=1; j<=v->p; j++) {
	    sprintf(lagstr, "_%d", j);
	    n = strlen(lagstr);
	    strncat(params[k], pdinfo->varname[v->ylist[i]],
		    VNAMELEN - n - 1);
	    strncat(params[k], lagstr, n);
	    k++;
	}
    }

    if (v->xlist != NULL) {
	for (i=1; i<=v->xlist[0]; i++) {
	    strcpy(params[k++], pdinfo->varname[v->xlist[i]]);
	}
    }

    if (v->detflags & DET_SEAS) {
	for (i=1; i<pdinfo->pd; i++) {
	    sprintf(params[k++], "S%d", i);
	}	
    }

    if (v->detflags & DET_TREND) {
	strcpy(params[k++], "time");
    }
}

/* (X'X)^{-1} * X'\Omega X * (X'X)^{-1} : right now
   we're supporting only HC0 and HC1 */

static int VAR_robust_vcv (gretl_matrix *V, VARspec *vspec, 
			   MODEL *pmod, int hcv, int k)
{
    gretl_matrix *XOX = NULL;
    double xij, xti, xtj, utk;
    int T = vspec->T;
    int g = vspec->g;
    int i, j, t;

    XOX = gretl_matrix_alloc(g, g);
    if (XOX == NULL) {
	return E_ALLOC;
    }

    /* form X' \Omega X */
    for (i=0; i<g; i++) {
	for (j=i; j<g; j++) {
	    xij = 0.0;
	    for (t=0; t<vspec->T; t++) {
		xti = gretl_matrix_get(vspec->X, t, i);
		xtj = gretl_matrix_get(vspec->X, t, j);
		utk = gretl_matrix_get(vspec->E, t, k);
		xij += utk * utk * xti * xtj;
	    }
	    if (hcv > 0) {
		/* cheating here, for now */
		xij *= (double) T / (T - g);
	    }
	    gretl_matrix_set(XOX, i, j, xij);
	    if (i != j) {
		gretl_matrix_set(XOX, j, i, xij);
	    }
	}
    }

    gretl_matrix_qform(vspec->XTX, GRETL_MOD_TRANSPOSE, XOX,
		       V, GRETL_MOD_NONE);

    gretl_model_set_int(pmod, "hc", 1);
    if (hcv > 0) {
	gretl_model_set_int(pmod, "hc_version", 1);
    }

    gretl_matrix_free(XOX);

    return 0;
}

/* Run the various per-equation omit tests (all lags of each var in
   turn, last lag of all vars) using the Wald method.  We also
   add the standard errors to the models here, since we have the
   covariance matrix to hand.
*/

static int VAR_wald_omit_tests (GRETL_VAR *var, VARspec *vspec, int ifc)
{
    gretl_matrix *V = NULL;
    gretl_matrix *C = NULL;
    gretl_vector *b = NULL;
    int hcv = get_hc_version();
    int p = vspec->p;
    int n = vspec->n;
    int g = vspec->g;
    int dim = (p > n)? p : n;
    int i, j, k, m = 0;
    int err = 0;

    if (ifc && vspec->robust && g - 1 > dim) {
	/* need bigger arrays for robust overall F-test */
	dim = g - 1;
    }

    V = gretl_matrix_alloc(g, g);
    C = gretl_matrix_alloc(dim, dim);
    b = gretl_column_vector_alloc(dim);

    if (V == NULL || C == NULL || b == NULL) {
	return E_ALLOC;
    }     

    for (i=0; i<n && !err; i++) {
	MODEL *pmod = var->models[i];
	int ii, jj, jpos, ipos = ifc;
	double w, vij;

	gretl_matrix_reuse(V, g, g);

	if (vspec->robust) {
	    err = VAR_robust_vcv(V, vspec, pmod, hcv, i);
	} else {
	    gretl_matrix_copy_values(V, vspec->XTX);
	    gretl_matrix_multiply_by_scalar(V, pmod->sigma * pmod->sigma);
	}
	
	if (!err) {
	    /* set (possibly robust) standard errors */
	    for (j=0; j<g; j++) {
		vij = gretl_matrix_get(V, j, j);
		pmod->sderr[j] = sqrt(vij);
	    }
	}

	/* exclusion of each var, all lags */

	gretl_matrix_reuse(C, p, p);
	gretl_matrix_reuse(b, p, 1);

	for (j=0; j<n && !err; j++) {
	    double w = NADBL;

	    gretl_matrix_extract_matrix(C, V, ipos, ipos, GRETL_MOD_NONE);
	    for (k=0; k<p; k++) {
		b->val[k] = pmod->coeff[k + ipos];
	    }
	    err = gretl_invert_symmetric_matrix(C);
	    if (!err) {
		w = gretl_scalar_qform(b, C, &err);
	    }
	    if (!err) {
		var->Fvals[m++] = w / p;
	    }

	    ipos += p;
	}

	/* exclusion of last lag, all vars? */

	if (p > 1) {
	    gretl_matrix_reuse(C, n, n);
	    gretl_matrix_reuse(b, n, 1);

	    ipos = ifc + p - 1;
	    for (ii=0; ii<n; ii++) {
		jpos = ifc + p - 1;
		for (jj=0; jj<n; jj++) {
		    vij = gretl_matrix_get(V, ipos, jpos);
		    gretl_matrix_set(C, ii, jj, vij);
		    jpos += p;
		}
		b->val[ii] = pmod->coeff[ipos];
		ipos += p;
	    }

	    err = gretl_invert_symmetric_matrix(C);
	    if (!err) {
		w = gretl_scalar_qform(b, C, &err);
	    }
	    if (!err) {
		var->Fvals[m++] = w / n;
	    }
	}

	/* exclusion of all but const? */

	if (ifc && vspec->robust) {
	    gretl_matrix_reuse(C, g-1, g-1);
	    gretl_matrix_reuse(b, g-1, 1);

	    gretl_matrix_extract_matrix(C, V, 1, 1, GRETL_MOD_NONE);
	    for (k=0; k<g-1; k++) {
		b->val[k] = pmod->coeff[k+1];
	    }
	    err = gretl_invert_symmetric_matrix(C);
	    if (!err) {
		w = gretl_scalar_qform(b, C, &err);
	    }
	    if (!err) {
		pmod->fstt = w / (g-1);
	    }
	}
    }

    gretl_matrix_free(V);
    gretl_matrix_free(C);
    gretl_matrix_free(b);

    return err;
}

/* make and record residuals for LR test, etc. */

static int 
last_lag_LR_prep (GRETL_VAR *var, VARspec *vspec, int ifc)
{
    int *collist = NULL;
    int g = vspec->g - vspec->n;
    int i, err = 0;

    if (var->F == NULL) {
	var->F = gretl_matrix_alloc(var->T, vspec->n);
	if (var->F == NULL) {
	    return E_ALLOC;
	}
    }   

    collist = gretl_list_new(vspec->n);
    if (collist == NULL) {
	return E_ALLOC;
    }

    collist[1] = ifc + vspec->p - 1;
    for (i=2; i<=collist[0]; i++) {
	collist[i] = collist[i-1] + vspec->p;
    }

    gretl_matrix_delete_columns(vspec->X, collist);
    gretl_matrix_reuse(vspec->B, g, vspec->n);
    err = gretl_matrix_multi_ols(vspec->Y, vspec->X, 
				 vspec->B, var->F,
				 NULL);

    free(collist);

    return err;
}

static int make_A_matrix (GRETL_VAR *var, VARspec *vspec, int ifc)
{
    int i, j, v, lag;
    int dim = vspec->n * vspec->p;
    double bij;

    for (j=0; j<vspec->n; j++) {
	v = lag = 0;
	for (i=0; i<dim; i++) {
	    bij = gretl_matrix_get(vspec->B, i+ifc, j);
	    gretl_matrix_set(var->A, j, vspec->n * lag + v, bij);
	    if (lag < vspec->p - 1) {
		lag++;
	    } else {
		lag = 0;
		v++;
	    }
	}
    }

    return 0;
}

static int transcribe_var_from_vspec (GRETL_VAR *var, VARspec *vspec,
				      int lagsel)
{
    int err = 0;

    var->ncoeff = vspec->g;
    var->t1 = vspec->t1;
    var->t2 = vspec->t2;
    var->T = var->t2 - var->t1 + 1;
    var->detflags = vspec->detflags;
    var->ylist = vspec->ylist;
    var->xlist = vspec->xlist;
    var->E = vspec->E;

    if (vspec->detflags & DET_CONST) {
	var->ifc = 1;
    }

    if (!lagsel) {
	/* not needed if we're just doing lag selection */

	err = make_A_matrix(var, vspec, var->ifc);

	if (!err) {
	    err = VAR_wald_omit_tests(var, vspec, var->ifc);
	}

	if (!err && vspec->p > 1) {
	    err = last_lag_LR_prep(var, vspec, var->ifc);
	}
    }

    if (!err) {
	err = VAR_add_stats(var);
    }

    if (!err && !lagsel) {
	err = gretl_VAR_do_error_decomp(var->S, var->C);
    }

    /* relinquish some pointers to var struct */
    vspec->E = NULL;
    vspec->ylist = NULL;
    vspec->xlist = NULL;

    return err;
}

static int set_up_VAR_models (GRETL_VAR *var, VARspec *vspec,
			      const double **Z,
			      const DATAINFO *pdinfo)
{
    MODEL *pmod;
    char **params = NULL;
    int yno, N = pdinfo->n;
    const double *y;
    int i, j;
    int err = 0;

    params = strings_array_new_with_length(vspec->g, VNAMELEN);
    if (params == NULL) {
	return E_ALLOC;
    }

    set_VAR_param_names(vspec, params, pdinfo);

    for (i=0; i<vspec->n && !err; i++) {
	yno = vspec->ylist[i+1];
	y = Z[yno];

	pmod = var->models[i];
	pmod->ID = i + 1;
	pmod->ci = VAR;
	pmod->aux = AUX_VAR;

	pmod->full_n = N;
	pmod->nobs = vspec->T;
	pmod->t1 = vspec->t1;
	pmod->t2 = vspec->t2;
	pmod->ncoeff = vspec->g;
	pmod->dfd = pmod->nobs - pmod->ncoeff;
	pmod->ifc = (vspec->detflags & DET_CONST)? 1 : 0;
	pmod->dfn = vspec->g - pmod->ifc;

	err = gretl_model_allocate_storage(pmod);
	pmod->depvar = gretl_strdup(pdinfo->varname[yno]);

	if (i == 0) {
	    pmod->params = params;
	} else {
	    pmod->params = strings_array_dup(params, vspec->g);
	}
	pmod->nparams = vspec->g;

	pmod->list = gretl_list_new(1);
	pmod->list[1] = yno;

	set_VAR_model_stats(pmod, vspec->E, y, i);

	for (j=0; j<vspec->g; j++) {
	    pmod->coeff[j] = gretl_matrix_get(vspec->B, j, i);
	}
    }

    return err;
}

/* Alternative mechanism for creating a GRETL_VAR using matrix OLS
   with multiple LHSs for a given RHS matrix.  As of now this almost
   duplicates the functionality of the old, clunky VAR apparatus.  The
   one thing not handled yet is production of robust standard errors.
*/

GRETL_VAR *alt_VAR (int order, int *list, double ***pZ, DATAINFO *pdinfo,
		    gretlopt opt, PRN *prn, int *errp)
{
    const double **Z = (const double **) *pZ;
    GRETL_VAR *var = NULL;
    VARspec *vspec = NULL;
    int lagsel = (opt & OPT_L)? 1 : 0;
    int err = 0;

    vspec = varspec_new(list, order, Z, pdinfo, opt, &err);

    if (!err) {
	err = make_VAR_matrices(vspec, Z, pdinfo);
    }

    /* run the regressions: use QR or Cholesky */
    if (vspec->qr) {
	err = gretl_matrix_QR_ols(vspec->Y, vspec->X, 
				  vspec->B, vspec->E,
				  &vspec->XTX, NULL);
    } else {
	err = gretl_matrix_multi_ols(vspec->Y, vspec->X, 
				     vspec->B, vspec->E,
				     &vspec->XTX);
    }

    if (!err) {
	var = gretl_VAR_new(VAR, vspec->n, vspec->p);
	if (var == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	if (lagsel) {
	    err = transcribe_var_from_vspec(var, vspec, 1);
	    if (!err) {
		err = alt_VAR_do_lagsel(var, vspec, Z, pdinfo, prn);
	    }
	} else {
	    err = set_up_VAR_models(var, vspec, Z, pdinfo);
	    if (!err) {
		err = transcribe_var_from_vspec(var, vspec, 0);
	    }
	    if (!err) {
		gretl_VAR_print(var, pdinfo, opt, prn);
	    }
	}
    }

    if (lagsel || (err && var != NULL)) {
	gretl_VAR_free(var);
	var = NULL;
    }

    varspec_free(vspec);

    *errp = err;

    return var;
}
