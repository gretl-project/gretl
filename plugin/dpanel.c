/* 
   Holder (possibly temporary) for the code to support the new
   "dpanel" command, which is intended to generalize the old "arbond"
   to handle system GMM.

   Added by Allin, 2010-07-09; contents are mostly from Jack's 
   prototype code, newmask.c.
*/

#define DPDSTYLE 0

#if 0 /* FIXME this needs to be hooked up or somehow
	 replaced */

static int dpd_add_residual_vector (dpdinfo *dpd, const double *y,
				    const double **Z, const DATAINFO *pdinfo)
{
    const double *b = dpd->beta->val;
    const double *x;
    double ut, xd;
    int i, j, s, t;

    if (dpd->uhat == NULL) {
	dpd->uhat = gretl_column_vector_alloc(dpd->nobs);
	if (dpd->uhat == NULL) {
	    return E_ALLOC;
	}
    }

    dpd->SSR = 0.0;
    s = 0;

    for (t=0; t<pdinfo->n; t++) {
	if (dpd->used[t]) {
	    ut = y[t] - y[t-1];
	    j = 0;
	    for (i=0; i<dpd->p; i++) {
		xd = y[t-i-1] - y[t-i-2];
		ut -= b[j++] * xd;
	    }
	    /* FIXME automatic time dummies */
	    if (dpd->xlist != NULL) {
		for (i=1; i<=dpd->xlist[0]; i++) {
		    x = Z[dpd->xlist[i]];
		    xd = x[t] - x[t-1];
		    ut -= b[j++] * xd;
		}
	    }
	    gretl_vector_set(dpd->uhat, s++, ut);
	    dpd->SSR += ut * ut;
	}
    }

    dpd->s2 = dpd->SSR / (dpd->nobs - dpd->k);

    return 0;
}

#endif /* 0 : FIXME */

/* 
   This is the first thing we do for each panel unit: construct a list
   holding the usable observations. A usable observation is defined
   as one for which we have all the data for the equation _in levels_.

   The list of good observations is a standard gretl list: the first
   element gives the number of good observations and the following
   elements give their positions, i.e. the 0-based places relative to
   the start of the time-series for the unit.

   Parameter @t0 gives the offset in the full dataset of the start of 
   the unit's data.

   We return the number of "good observations" minus one, to allow
   for differencing; the return value (if positive) will be the
   max number of observations that are actually usable.
*/

static int check_unit_obs (dpdinfo *dpd, int *goodobs, const double **Z, 
			   int t0)
{
    int i, s, t, ok;

    goodobs[0] = 0;

    for (t=0; t<dpd->T; t++) {
	/* do we have the dependent variable? */
	ok = !na(Z[dpd->yno][t+t0]);

	/* lags of dependent variable? */
	for (i=1; i<=dpd->laglist[0] && ok; i++) {
	    s = t - dpd->laglist[i];
	    if (s < 0) {
		ok = 0;
	    } else {
		ok &= !na(Z[dpd->yno][s+t0]);
	    }
	}

	if (ok && dpd->xlist != NULL) {
	    /* regressors */
	    for (i=1; i<=dpd->xlist[0] && ok; i++) {
		ok &= !na(Z[dpd->xlist[i]][t+t0]);
	    }
	}

	if (ok) {
	    goodobs[0] += 1;
	    goodobs[goodobs[0]] = t;
	}
    }
    
    /* allow for differencing */
    return goodobs[0] - 1;
}

/* Based on the accounting of good observations for a unit, recorded
   in the @goodobs list, fill matrix D, which will be used to
   construct H.
*/

static void build_unit_D_matrix (dpdinfo *dpd, int *goodobs, gretl_matrix *D)
{
    int usable = goodobs[0] - 1;
    int maxlag = dpd->p;
    int i, j, i0, i1;    

    /* zero all elements */
    gretl_matrix_zero(D);

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	j = i1 - 1 - maxlag;
	gretl_matrix_set(D, i0, j, -1);
	gretl_matrix_set(D, i1, j,  1);

    }

    /* levels */
    if (dpd->flags & DPD_SYSTEM) {
	for (i=1; i<=goodobs[0]; i++) {
	    i1 = goodobs[i];
	    j = i1 - 1 - maxlag;
	    gretl_matrix_set(D, i1, j + dpd->T - maxlag, 1);
	}
    }

#if ADEBUG > 2
    gretl_matrix_print_to_prn(D, "D", prn);
#endif
}

/* Variant computations of H: for "dpdstyle" we emulate what
   the Ox DPD package does */

static void compute_H (const gretl_matrix *D, gretl_matrix *H, 
		       int ni_d, int dpdstyle)
{
    if (!dpdstyle) {
	gretl_matrix_multiply_mod(D, GRETL_MOD_TRANSPOSE, 
				  D, GRETL_MOD_NONE, 
				  H, GRETL_MOD_NONE);
    } else {
	int i;

	gretl_matrix_set(H, 0, 0, 2);
	for (i=1; i<ni_d; i++) {
	    gretl_matrix_set(H, i, i, 2);
	    gretl_matrix_set(H, i-1, i, -1);
	    gretl_matrix_set(H, i, i-1, -1);
	}
	gretl_matrix_multiply_by_scalar(H, 1); /* hello? */
    }
}

/* Build row vector of dependent variable values in @dep using
   differences, followed by levels if wanted.
*/

static void build_Y (dpdinfo *dpd, int *goodobs, const double **Z, 
		     int t, gretl_matrix *dep)
{
    int i, usable = goodobs[0] - 1;
    int maxlag = dpd->p;
    int T = dpd->T;
    int t0, t1, i0, i1;
    double y0, y1;

    /* zero all elements */
    gretl_matrix_zero(dep);

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	t0 = t + i0;
	t1 = t + i1;
	y0 = Z[dpd->yno][t0];
	y1 = Z[dpd->yno][t1];
	gretl_vector_set(dep, i1-1-maxlag, y1 - y0);
    }
    
    if (dpd->flags & DPD_SYSTEM) {
	/* levels */
	for (i=0; i<=usable; i++) {
	    i1 = goodobs[i+1];
	    t1 = t + i1;
	    y1 = Z[dpd->yno][t1];
	    gretl_vector_set(dep, (T-1-maxlag) + (i1-maxlag), y1);
	}
    }
}

/* Build matrix of right-hand side variable values in @ind using
   differences, followed by levels if wanted.
*/

static void build_X (dpdinfo *dpd, int *goodobs, const double **Z, 
		     int t, gretl_matrix *ind)
{
    int i, j, usable = goodobs[0] - 1;
    int nlags = dpd->laglist[0];
    int maxlag = dpd->p;
    int t0, t1, i0, i1;
    double x0, x1;

    /* zero all elements */
    gretl_matrix_zero(ind);

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	t0 = t+i0;
	t1 = t+i1;

	for (j=1; j<=nlags; j++) {
	    x0 = Z[dpd->yno][t0 - dpd->laglist[j]];
	    x1 = Z[dpd->yno][t1 - dpd->laglist[j]];
	    gretl_matrix_set(ind, j-1, (i1-1-maxlag), x1 - x0);
	}

	for (j=1; j<=dpd->nx; j++) {
	    x0 = Z[dpd->xlist[j]][t0];
	    x1 = Z[dpd->xlist[j]][t1];
	    gretl_matrix_set(ind, j+nlags-1, (i1-1-maxlag), x1 - x0);
	}
    }
    
    if (dpd->flags & DPD_SYSTEM) {
	/* levels */
	int col;

	for (i=0; i<=usable; i++) {
	    i1 = goodobs[i+1];
	    t1 = t+i1;
	    col = (dpd->T-1-maxlag) + (i1-maxlag);

	    for (j=1; j<=nlags; j++) {
		x1 = Z[dpd->yno][t1 - dpd->laglist[j]];
		gretl_matrix_set(ind, j-1, col, x1);
	    }

	    for (j=1; j<=dpd->nx; j++) {
		x1 = Z[dpd->xlist[j]][t1];
		gretl_matrix_set(ind, j+nlags-1, col, x1);
	    }
	}
    }
}

/* Build matrix of instrument values in @ind: instruments
   for the equations in differences come first, followed by
   instruments for equations in levels if wanted.
*/

static void build_Z (dpdinfo *dpd, int *goodobs, const double **Z, 
		     int t, int ni_d, int ni_l, 
		     gretl_matrix *inst)
{
    int usable = goodobs[0] - 1;
    int maxlag = dpd->p;
    int T = dpd->T;
    double x0, x1;
    int t0, t1, i0, i1;
    int k, k2 = ni_d + ni_l;
    int i, j;

    /* zero all elements */
    gretl_matrix_zero(inst);

    /* equations in differences -- lagged levels of y */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	k = i0 * (i0-1) / 2;
	for (j=0; j<i0; j++) {
	    x0 = Z[dpd->yno][t+j];
	    if (!na(x0)) {
		gretl_matrix_set(inst, k, i1-1-maxlag, x0);
	    }
	    k++;
	}
    }

    /* equations in differences -- differenced exog vars */
    if (dpd->nx > 0) {
	for (i=0; i<usable; i++) {
	    i0 = goodobs[i+1];
	    i1 = goodobs[i+2];
	    t0 = t + i0;
	    t1 = t + i1;
	    for (j=0; j<dpd->nx; j++) {
		x0 = Z[dpd->xlist[j+1]][t0];
		x1 = Z[dpd->xlist[j+1]][t1];
		gretl_matrix_set(inst, k2 + j, (i1-1-maxlag), x1 - x0);
	    }
	}
    }

    if (dpd->flags & DPD_SYSTEM) {
	int offset, col, row;
	int lastdiff = 0;

	/* equations in levels -- lagged differences of y */

	row = (T-1)*(T-2)/2;
	offset = T-maxlag-1;

	for (i=0; i<=usable; i++) {
	    i0 = goodobs[i+1];
	    col = offset + i0 - maxlag;
	    for (j=lastdiff; j<i0; j++) {
		x0 = (j < 1) ? NADBL : Z[dpd->yno][t+j-1];
		x1 = Z[dpd->yno][t+j];
		if (!na(x1) && !na(x0)) {
		    gretl_matrix_set(inst, row, col, x1 - x0);
		}
		row++;
	    }
	    lastdiff = j;
	}

	/* equations in levels -- exog vars */
	if (dpd->nx > 0) {
	    for (i=0; i<=usable; i++) {
		i1 = goodobs[i+1];
		t1 = t+i1;
		col = (T-1-maxlag) + (i1-maxlag);
		for (j=0; j<dpd->nx; j++) {
		    x1 = Z[dpd->xlist[j+1]][t1];
		    gretl_matrix_set(inst, k2 + j, col, x1);
		}
	    }
	}
    }
}

/*
  here we construct \hat{\beta} from the moment matrices;
  no attempt is made to compute the covariance matrix
*/

static int do_estimator (dpdinfo *dpd)
{
    gretl_matrix *iZZ, *M, *M1, *M2;
    int k, m, err = 0;
    
    trim_zero_inst(dpd->XZ, dpd->ZZ, dpd->ZY);

#if ADEBUG > 0
    gretl_matrix_print(dpd->XZ, "XZ (trimmed)");
    gretl_matrix_print(dpd->ZZ, "ZZ (trimmed)");
    gretl_matrix_print(dpd->ZY, "ZY (trimmed)");
#endif

    k = dpd->XZ->rows;
    m = dpd->XZ->cols;

    iZZ = gretl_matrix_copy(dpd->ZZ);
    M = gretl_matrix_alloc(k, m);
    M1 = gretl_matrix_alloc(k, k);
    M2 = gretl_matrix_alloc(k, 1);

    if (iZZ == NULL || M == NULL || M1 == NULL || M2 == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	/* symmetrize iZZ; you never know */
	gretl_matrix_xtr_symmetric(iZZ);
	err = gretl_invert_symmetric_matrix(iZZ);
    }

    if (!err) {
	gretl_matrix_multiply(dpd->XZ, iZZ, M);
	gretl_matrix_multiply_mod(M, GRETL_MOD_NONE, 
				  dpd->XZ, GRETL_MOD_TRANSPOSE,
				  M1, GRETL_MOD_NONE);
	gretl_matrix_xtr_symmetric(M1);
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(M1);
    }

    if (!err) {
	gretl_matrix_multiply(M, dpd->ZY, M2);
	gretl_matrix_multiply(M1, M2, dpd->beta);
    }

#if ADEBUG > 0
    gretl_matrix_print(M1, "M1");
    gretl_matrix_print(M2, "M2");
#endif

    gretl_matrix_free(iZZ);
    gretl_matrix_free(M);
    gretl_matrix_free(M1);
    gretl_matrix_free(M2);

    return err;
}

/* allocate storage needed in do_units() */

static int dpd_allocate_new (dpdinfo *dpd, int nz, int north,
			     gretl_matrix **D, gretl_matrix **H,
			     gretl_matrix **dep, gretl_matrix **ind,
			     gretl_matrix **inst)
{
    dpd->XZ = gretl_zero_matrix_new(dpd->k, nz);
    dpd->ZZ = gretl_zero_matrix_new(nz, nz);
    dpd->ZY = gretl_zero_matrix_new(nz, 1);

    if (dpd->XZ == NULL || dpd->ZZ == NULL || dpd->ZY == NULL) {
	return E_ALLOC;
    }
    
    *D = gretl_matrix_alloc(dpd->T, north);
    *H = gretl_identity_matrix_new(north);
    *dep = gretl_matrix_alloc(1, north);
    *ind = gretl_matrix_alloc(dpd->k, north);
    *inst = gretl_matrix_alloc(nz, north);

    if (*D == NULL || *H == NULL || *dep == NULL ||
	*ind == NULL || *inst == NULL) {
	return E_ALLOC;
    }

    return 0;
}

/* Main driver for system GMM: the core is a loop across
   the panel units to build the moment matrices.
*/

static int do_units (dpdinfo *dpd, const double **Z, 
		     const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *D = NULL;
    gretl_matrix *dep = NULL;
    gretl_matrix *ind = NULL;
    gretl_matrix *inst = NULL;
    gretl_matrix *H = NULL;
    int *goodobs;
    int T = dpd->T;
    int t, unit = 0;
    int maxlag = dpd->p;
    int north, nz, nz_d, nz_l;
    int err = 0;

    /* do H as per Ox/dpd? */
    int dpdstyle = DPDSTYLE;

    if (dpd->flags & DPD_SYSTEM) {
	north = (T-maxlag-1) * 2 + 1;
    } else {
	north = T-maxlag-1;
    }

    /* number of instruments */
    nz = nz_d = (T-1)*(T-2)/2;
    if (dpd->flags & DPD_SYSTEM) {
	nz_l = T-1;
	nz += nz_l;
    } else {
	nz_l = 0;
    }
    nz += dpd->nx;

#if ADEBUG > 0
    pprintf(prn, "nz = %4d:\n", nz);
#endif

    goodobs = gretl_list_new(T);
    if (goodobs == NULL) {
	return E_ALLOC;
    }

    err = dpd_allocate_new(dpd, nz, north, &D, &H,
			   &dep, &ind, &inst);
    if (err) {
	goto bailout;
    }

    /* initialize observation counts */
    dpd->effN = dpd->nobs = 0;

    for (t=pdinfo->t1; t<pdinfo->t2; t+=T) {
	int uT = check_unit_obs(dpd, goodobs, Z, t);

#if ADEBUG > 0
	pprintf(prn, "\n\nUnit %4d:", unit);
	pprintf(prn, "usable obs = %4d:\n", uT);
#endif

	if (uT > 0) {
	    /* this unit is usable */
	    dpd->effN += 1;
	    dpd->nobs += uT;
	    build_unit_D_matrix(dpd, goodobs, D);
	    compute_H(D, H, T-maxlag-1, dpdstyle);
#if ADEBUG > 2
	    if (t == pdinfo->t1) {
		gretl_matrix_print_to_prn(H, "H", prn);
	    }
#endif
	    build_Y(dpd, goodobs, Z, t, dep);
	    build_X(dpd, goodobs, Z, t, ind);
	    build_Z(dpd, goodobs, Z, t, nz_d, nz_l, inst);
#if ADEBUG > 1
	    gretl_matrix_print_to_prn(dep, "Y", prn);
	    gretl_matrix_print_to_prn(ind, "X", prn);
	    gretl_matrix_print_to_prn(inst, "Z", prn);
#endif
	    gretl_matrix_multiply_mod(ind, GRETL_MOD_NONE,
				      inst, GRETL_MOD_TRANSPOSE,
				      dpd->XZ, GRETL_MOD_CUMULATE);
	    gretl_matrix_qform(inst, GRETL_MOD_NONE,
			       H, dpd->ZZ, GRETL_MOD_CUMULATE);
	    gretl_matrix_multiply_mod(inst, GRETL_MOD_NONE,
				      dep, GRETL_MOD_TRANSPOSE,
				      dpd->ZY, GRETL_MOD_CUMULATE);
	}
	unit++;
    }

    pprintf(prn, "Total obs = %d\n", dpd->nobs);

 bailout:

    gretl_matrix_free(D);
    gretl_matrix_free(H);
    gretl_matrix_free(dep);
    gretl_matrix_free(ind);
    gretl_matrix_free(inst);

    free(goodobs);

    return err;
}
