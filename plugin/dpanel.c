/* 
   Holder (possibly temporary) for the code to support the new
   "dpanel" command, which is intended to generalize the old "arbond"
   to handle system GMM.

   Added by Allin, 2010-07-09; initial content was from Jack's 
   prototype code, newmask.c.
*/

#define DPDEBUG 0

#define use_levels(d) (d->flags & DPD_SYSTEM)

static int dpanel_step1_variance (dpdinfo *dpd)
{
    gretl_matrix_block *B;
    gretl_matrix *kk, *kz, *V;
    gretl_matrix *ui, *uZ;
    int i, s, t;
    int Zcol = 0;
    int err = 0;

    /* Unlike the other matrices used here, we want to preserve
       V for use in the second step (if applicable), so we'll
       not make V part of matrix block B. 
    */
    V = gretl_zero_matrix_new(dpd->nz, dpd->nz);
    if (V == NULL) {
	return E_ALLOC;
    }

    /* max number of stacked observations per unit */
    dpd->max_ni = dpd->maxTi;
    if (dpd->flags & DPD_SYSTEM) {
	dpd->max_ni += dpd->maxTi + 1;
    }

#if DPDEBUG
    fprintf(stderr, "nz = %d, max_ni = %d\n", nz, dpd->max_ni);
#endif

    B = gretl_matrix_block_new(&kk, dpd->k, dpd->k,
			       &kz, dpd->k, dpd->nz,
			       &ui, dpd->max_ni, 1,
			       &uZ, 1, dpd->nz,
			       NULL);
    if (B == NULL) {
	gretl_matrix_free(V);
	return E_ALLOC;
    }

    s = 0;

    for (i=0; i<dpd->N; i++) {
	int ni = dpd->ui[i].nobs;

	if (ni == 0) {
	    continue;
	}

	gretl_matrix_reuse(dpd->Zi, ni, dpd->nz);
	gretl_matrix_reuse(ui, ni, -1);

	/* extract ui */
	for (t=0; t<ni; t++) {
	    ui->val[t] = dpd->uhat->val[s++];
	}	

	/* extract Zi */
	gretl_matrix_extract_matrix(dpd->Zi, dpd->ZT, 0, Zcol,
				    GRETL_MOD_TRANSPOSE);
	Zcol += dpd->Zi->rows;

#if DPDEBUG
	gretl_matrix_print(ui, "ui (vcalc)");
	gretl_matrix_print(Zi, "Zi (vcalc)");
#endif

	err = gretl_matrix_multiply_mod(ui, GRETL_MOD_TRANSPOSE,
					dpd->Zi, GRETL_MOD_NONE,
					uZ, GRETL_MOD_NONE);
	if (!err) {
	    err = gretl_matrix_multiply_mod(uZ, GRETL_MOD_TRANSPOSE,
					    uZ, GRETL_MOD_NONE,
					    V, GRETL_MOD_CUMULATE);
	}

	if (err) {
	    fprintf(stderr, "dpanel_step1_variance: error at unit %d (ni = %d)\n", 
		    i, ni);
	    break;
	}
    }

    if (!err) {
	gretl_matrix_divide_by_scalar(V, dpd->effN);
	err = gretl_matrix_multiply(dpd->XZ, dpd->A, kz);
	if (!err) {
	    err = gretl_matrix_qform(kz, GRETL_MOD_NONE, V,
				     kk, GRETL_MOD_NONE);
	}
    }

#if DPDEBUG > 1
    gretl_matrix_print(V, "V");
    gretl_matrix_print(dpd->XZ, "XZ");
    gretl_matrix_print(dpd->A, "A");
    gretl_matrix_print(dpd->den, "den");
#endif

    if (!err) {
	/* pre- and post-multiply by den */
	err = gretl_matrix_qform(dpd->den, GRETL_MOD_NONE, kk, 
				 dpd->vbeta, GRETL_MOD_NONE);
    }

    if (!err) {
	gretl_matrix_multiply_by_scalar(dpd->vbeta, dpd->effN);
    }

    gretl_matrix_block_destroy(B);

    if (!err && (dpd->flags & DPD_TWOSTEP)) {
	/* preserve V for second stage */
	dpd->V = V;
    } else {
	gretl_matrix_free(V);
    }

    if (err) {
	fprintf(stderr, "dpanel_step1_variance: err = %d (nz = %d)\n", 
		err, dpd->nz);
    }

    return err;
}

static int dpanel_step2_variance (dpdinfo *dpd, 
				  gretl_matrix *u1,
				  gretl_matrix *V1)
{
    int err;

    err = gretl_matrix_qform(dpd->XZ, GRETL_MOD_NONE, dpd->V,
			     dpd->vbeta, GRETL_MOD_NONE);

    if (!err) {
	err = gretl_invert_symmetric_matrix(dpd->vbeta);
    }

    if (!err) {
	gretl_matrix_multiply_by_scalar(dpd->vbeta, dpd->effN);
    }

    if (!err && u1 != NULL && V1 != NULL) {
	dpd->kmtmp = gretl_matrix_alloc(dpd->k, dpd->nz);
	dpd->R1 = gretl_matrix_alloc(dpd->nz, 1);

	if (dpd->kmtmp == NULL || dpd->R1 == NULL) {
	    err = E_ALLOC;
	} else {
	    err = windmeijer_correct(dpd, u1, V1);
	}
    }    

    return err;
}

/* Populate the residual vector, dpd->uhat. In the system case
   we stack the residuals in levels under the residuals in
   differences, per unit. Calculate SSR and \sigma^2 while
   we're at it.
*/

static int dpanel_add_residuals (dpdinfo *dpd)
{
    const double *b = dpd->beta->val;
    double SSRd = 0.0, SSRl = 0.0;
    double ut;
    int i, j, s, t;

    /* allocate if we haven't already done so */

    if (dpd->uhat == NULL) {
	int totobs = dpd->ndiff + dpd->nlev;

	dpd->uhat = gretl_column_vector_alloc(totobs);
	if (dpd->uhat == NULL) {
	    return E_ALLOC;
	}
    }

    s = 0;

    for (i=0; i<dpd->N; i++) {
	for (t=0; t<dpd->ui[i].t1; t++) {
	    /* differences */
	    ut = dpd->Y->val[s];
	    for (j=0; j<dpd->k; j++) {
		ut -= b[j] * gretl_matrix_get(dpd->X, s, j);
	    }
	    gretl_vector_set(dpd->uhat, s, ut);
	    SSRd += ut * ut;
	    s++;
	}
	for (t=0; t<dpd->ui[i].t2; t++) {
	    /* levels, it applicable */
	    ut = dpd->Y->val[s];
	    for (j=0; j<dpd->k; j++) {
		ut -= b[j] * gretl_matrix_get(dpd->X, s, j);
	    }
	    gretl_vector_set(dpd->uhat, s, ut);
	    SSRl += ut * ut;
	    s++;
	}
    }

    if (use_levels(dpd)) {
	/* we could attach these to the final model */
	fprintf(stderr, "nobs (diffs) = %d\n", dpd->ndiff);
	fprintf(stderr, "SSR (diffs) = %g\n", SSRd);
	dpd->nobs = dpd->nlev;
	dpd->SSR = SSRl;
    } else {
	dpd->nobs = dpd->ndiff;
	dpd->SSR = SSRd;
    }

    dpd->s2 = dpd->SSR / (dpd->nobs - dpd->k);

    return 0;
}

/* 
   This is the first thing we do for each panel unit: construct a list
   holding the usable observations. A usable observation is defined
   as one for which we have all the data for the equation _in levels_.

   The list of good observations takes this form: the first element
   gives the number of good observations and the following elements
   give their positions, i.e. the 0-based places relative to the start
   of the time-series for the unit. While we're at it we record the
   positions of the good observations relative to the full dataset,
   in the big array dpd->used, so that we can place the residuals
   correctly later on.

   (Note that @t0 gives the offset in the full dataset of the start of 
   the unit's data.)

   We return the number of "good observations" minus one, to allow
   for differencing; the return value (if positive) will be the
   max number of observations that are actually usable.
*/

static int check_unit_obs (dpdinfo *dpd, int *goodobs, const double **Z, 
			   int t0)
{
    const double *y = Z[dpd->yno];
    int i, s, t, ok;

    goodobs[0] = 0;

    for (t=0; t<dpd->T; t++) {
	int big_t = t + t0;

	/* do we have the dependent variable? */
	ok = !na(y[big_t]);

	/* lags of dependent variable? */
	for (i=1; i<=dpd->laglist[0] && ok; i++) {
	    s = t - dpd->laglist[i];
	    if (s < 0) {
		ok = 0;
	    } else {
		ok &= !na(y[s+t0]);
	    }
	}

	if (ok && dpd->xlist != NULL) {
	    /* regressors */
	    for (i=1; i<=dpd->xlist[0] && ok; i++) {
		ok &= !na(Z[dpd->xlist[i]][big_t]);
	    }
	}

	if (ok) {
	    goodobs[0] += 1;
	    goodobs[goodobs[0]] = t;
	    if (goodobs[0] > 1 || use_levels(dpd)) {
		dpd->used[big_t] = 1;
	    } 
	}
    }
    
    /* allow for differencing */
    return goodobs[0] - 1;
}

/* Based on the accounting of good observations for a unit recorded
   in the @goodobs list, fill matrix D (which will be used to
   construct H unless we're doing things "dpdstyle").
*/

static void build_unit_D_matrix (dpdinfo *dpd, int *goodobs, gretl_matrix *D)
{
    int usable = goodobs[0] - 1;
    int maxlag = dpd->p;
    int i, j, i0, i1;    

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
    if (use_levels(dpd)) {
	for (i=1; i<=goodobs[0]; i++) {
	    i1 = goodobs[i];
	    j = i1 - 1 - maxlag;
	    gretl_matrix_set(D, i1, j + dpd->T - maxlag, 1);
	}
    }

#if DPDEBUG > 2
    gretl_matrix_print(D, "D");
#endif
}

/* Variant computations of H: for "dpdstyle" we emulate what
   the Ox DPD package does, otherwise we use the pre-formed
   D matrix.
*/

static void compute_H (gretl_matrix *H, int nh,
		       const gretl_matrix *D)
{
    if (D != NULL) {
	gretl_matrix_multiply_mod(D, GRETL_MOD_TRANSPOSE, 
				  D, GRETL_MOD_NONE, 
				  H, GRETL_MOD_NONE);
    } else {	
	int i;

	gretl_matrix_set(H, 0, 0, 2);
	for (i=1; i<nh; i++) {
	    gretl_matrix_set(H, i, i, 2);
	    gretl_matrix_set(H, i-1, i, -1);
	    gretl_matrix_set(H, i, i-1, -1);
	}
#if 0
	gretl_matrix_multiply_by_scalar(H, scale);
#endif
    } 
}

/* Build row vector of dependent variable values in @Yi using
   differences, followed by levels if wanted. 
*/

static void build_Y (dpdinfo *dpd, int *goodobs, const double **Z, 
		     int t, gretl_matrix *Yi)
{
    const double *y = Z[dpd->yno];
    int i, usable = goodobs[0] - 1;
    int maxlag = dpd->p;
    int T = dpd->T;
    int Tshort = T - maxlag - 1;
    int t0, t1, i0, i1;
    double dy;

    gretl_matrix_zero(Yi);

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	t0 = t + i0;
	t1 = t + i1;
	dy = y[t1] - y[t0];
	gretl_vector_set(Yi, i1-1-maxlag, dy);
    }
    
    if (use_levels(dpd)) {
	/* levels */
	for (i=0; i<=usable; i++) {
	    i1 = goodobs[i+1];
	    t1 = t + i1;
	    gretl_vector_set(Yi, Tshort + (i1-maxlag), y[t1]);
	}
    }
}

/* Build matrix of right-hand side variable values in @Xi using
   differences, followed by levels if wanted.
*/

static void build_X (dpdinfo *dpd, int *goodobs, const double **Z, 
		     int t, gretl_matrix *Xi)
{
    const double *y = Z[dpd->yno];
    int usable = goodobs[0] - 1;
    int nlags = dpd->laglist[0];
    int maxlag = dpd->p;
    int Tshort = dpd->T - maxlag - 1;
    const double *xj;
    int t0, t1, i0, i1;
    int i, j, lj;
    double dx;

    gretl_matrix_zero(Xi);

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	t0 = t + i0;
	t1 = t + i1;

	for (j=1; j<=nlags; j++) {
	    lj = dpd->laglist[j];
	    dx = y[t1-lj] - y[t0-lj];
	    gretl_matrix_set(Xi, j-1, i1-1-maxlag, dx);
	}

	for (j=1; j<=dpd->nx; j++) {
	    /* Note: for now not differencing out the constant
	       Allin 2010-07-22 */
	    if (!use_levels(dpd) && dpd->xlist[j] == 0) {
		dx = 1.0;
	    } else {
		xj = Z[dpd->xlist[j]];
		dx = xj[t1] - xj[t0];
	    }
	    gretl_matrix_set(Xi, j+nlags-1, i1-1-maxlag, dx);
	}
    }
    
    if (use_levels(dpd)) {
	/* levels */
	int col;

	for (i=0; i<=usable; i++) {
	    i1 = goodobs[i+1];
	    t1 = t+i1;
	    col = Tshort + (i1-maxlag);

	    for (j=1; j<=nlags; j++) {
		lj = dpd->laglist[j];
		gretl_matrix_set(Xi, j-1, col, y[t1-lj]);
	    }

	    for (j=1; j<=dpd->nx; j++) {
		xj = Z[dpd->xlist[j]];
		gretl_matrix_set(Xi, j+nlags-1, col, xj[t1]);
	    }
	}
    }
}

/* Build matrix of instrument values in @Zi: instruments
   for the equations in differences come first, followed by
   instruments for equations in levels if wanted.
*/

static void build_Z (dpdinfo *dpd, int *goodobs, const double **Z, 
		     int t, int nz_diff, int nz_lev, 
		     gretl_matrix *Zi)
{
    const double *y = Z[dpd->yno];
    int usable = goodobs[0] - 1;
    int maxlag = dpd->p;
    int T = dpd->T;
    const double *xj;
    double dx, y0;
    int t0, t1, i0, i1;
    int k, k2 = nz_diff + nz_lev;
    int i, j, col;
    int qmax = dpd->qmax;

    gretl_matrix_zero(Zi);

    /* equations in differences: lagged levels of y */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	col = i1 - 1 - maxlag;
	k = i0 * (i0 - 1) / 2;
#if DPDEBUG
	fprintf(stderr, "Zi: i0=%d, i1=%d, base row=%d\n", i0, i1, k);
#endif
	for (j=0; j<i0; j++) {
	    if (qmax == 0 || j > i0 - qmax) {
		y0 = y[t+j];
		if (!na(y0)) {
		    gretl_matrix_set(Zi, k, col, y0);
		} 
	    } 
	    k++;
	}
    }

    /* equations in differences: differenced exog vars */
    if (dpd->nx > 0) {
	for (i=0; i<usable; i++) {
	    i0 = goodobs[i+1];
	    i1 = goodobs[i+2];
	    col = i1 - 1 - maxlag;
	    t0 = t + i0;
	    t1 = t + i1;
	    for (j=0; j<dpd->nx; j++) {
		/* Allin 2010-07-22: not differencing out the constant */
		if (!use_levels(dpd) && dpd->xlist[j+1] == 0) {
		    dx = 1.0;
		} else {
		    xj = Z[dpd->xlist[j+1]];
		    dx = xj[t1] - xj[t0];
		}
		gretl_matrix_set(Zi, k2 + j, col, dx);
	    }
	}
    }

    if (use_levels(dpd)) {
	int offset = T - maxlag - 1;
	int row = (T-1)*(T-2)/2;
	int lastdiff = 0;
	double y1;

	/* equations in levels: lagged differences of y */
	for (i=0; i<=usable; i++) {
	    i0 = goodobs[i+1];
	    col = offset + i0 - maxlag;
	    for (j=lastdiff; j<i0; j++) {
		if (qmax == 0 || j > i0 - qmax) {
		    y0 = (j < 1)? NADBL : y[t+j-1];
		    y1 = y[t+j];
		    if (!na(y1) && !na(y0)) {
			gretl_matrix_set(Zi, row, col, y1 - y0);
		    }
		}
		row++;
	    }
	    lastdiff = j;
	}

	/* equations in levels: levels of exog vars */
	if (dpd->nx > 0) {
	    for (i=0; i<=usable; i++) {
		i1 = goodobs[i+1];
		t1 = t + i1;
		col = (T-1-maxlag) + (i1-maxlag);
		for (j=0; j<dpd->nx; j++) {
		    xj = Z[dpd->xlist[j+1]];
		    gretl_matrix_set(Zi, k2 + j, col, xj[t1]);
		}
	    }
	}
    }
}

static int trim_zero_inst (dpdinfo *dpd)
{
    int i, n = dpd->ZZ->rows;
    int trim = 0;

    for (i=0; i<n; i++) {
	if (gretl_matrix_get(dpd->ZZ, i, i) == 0.0) {
	    trim = 1;
	    break;
	}
    }

    if (trim) {
	char *mask = calloc(n, 1);

	for (i=0; i<n; i++) {
	    if (gretl_matrix_get(dpd->ZZ, i, i) == 0.0) {
		mask[i] = 1;
	    }
	}

	gretl_matrix_cut_cols(dpd->XZ, mask);
	gretl_matrix_cut_rows_cols(dpd->ZZ, mask);
	gretl_matrix_cut_rows(dpd->ZY, mask);
	gretl_matrix_cut_rows(dpd->ZT, mask);

	free(mask);
    }

    return dpd->ZZ->rows;
}

/* Here we compute \hat{\beta} from the moment matrices. */

static int do_estimator (dpdinfo *dpd)
{
    gretl_matrix *M1 = NULL, *M2 = NULL;
    int k, m, err = 0;

    if (dpd->step == 1) {
	dpd->nz = trim_zero_inst(dpd);
#if DPDEBUG > 1
	gretl_matrix_print(dpd->XZ, "XZ (trimmed)");
	gretl_matrix_print(dpd->ZZ, "ZZ (trimmed)");
	gretl_matrix_print(dpd->ZY, "ZY (trimmed)");
#endif
    }

    if (dpd->A == NULL) {
	/* first step: A will be set to ZZ^{-1} */
	dpd->A = gretl_matrix_copy(dpd->ZZ);
	if (dpd->A == NULL) {
	    return E_ALLOC;
	}
    }

    k = dpd->XZ->rows;
    m = dpd->XZ->cols;

    if (dpd->XZA == NULL) {
	dpd->XZA = gretl_matrix_alloc(k, m);
	if (dpd->XZA == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	M1 = gretl_matrix_alloc(k, k);
	M2 = gretl_matrix_alloc(k, 1);
	if (M1 == NULL || M2 == NULL) {
	    err = E_ALLOC;
	}
    }

    if (dpd->step == 1 && !err) {
	/* symmetrize A; you never know */
	gretl_matrix_xtr_symmetric(dpd->A);
	err = gretl_invert_symmetric_matrix(dpd->A);
    }

    if (!err) {
	/* M1 <- XZ * A * XZ' */
	gretl_matrix_qform(dpd->XZ, GRETL_MOD_NONE,
			   dpd->A, M1, GRETL_MOD_NONE);
    }

    if (!err) {
	/* M1 <- (XZ * A * XZ')^{1} */
	err = gretl_invert_symmetric_matrix(M1);
    }

    if (!err) {
	/* dpd->XZA ("M") <- XZ * A */
	gretl_matrix_multiply(dpd->XZ, dpd->A, dpd->XZA);
	/* M2 <- XZ * A * ZY */
	gretl_matrix_multiply(dpd->XZA, dpd->ZY, M2);
	/* beta <- (XZ * A * XZ')^{1} * (XZ * A * ZY) */
	gretl_matrix_multiply(M1, M2, dpd->beta);
    }

#if DPDEBUG > 1
    gretl_matrix_print(M1, "M1");
    gretl_matrix_print(M2, "M2");
#endif

    gretl_matrix_free(M2);

    if (!err) {
	/* wanted for computing the variance of \hat{\beta} */
	dpd->den = M1;
    } else {
	gretl_matrix_free(M1);
    }

    return err;
}

static int do_estimator_2 (dpdinfo *dpd)
{
    gretl_matrix *Vcpy;
    int err = 0;

#if DPDEBUG
    gretl_matrix_print(dpd->V, "V, in do_estimator_2");
#endif

    if (gretl_matrix_rows(dpd->V) > dpd->effN) {
	return E_DF; /* this case won't work? */
    }

    Vcpy = gretl_matrix_copy(dpd->V);
    if (Vcpy == NULL) {
	return E_ALLOC;
    }

    err = gretl_invert_symmetric_matrix(dpd->V);
    if (err) {
	/* revert the data in dpd->V */
	gretl_matrix_copy_values(dpd->V, Vcpy);
    }

    if (err) {
	err = gretl_SVD_invert_matrix(dpd->V);
	if (!err) {
	    gretl_matrix_xtr_symmetric(dpd->V);
	}
    }

    if (!err) {
	/* A <- V^{-1} */
	gretl_matrix_copy_values(dpd->A, dpd->V);
	dpd->step = 2;
	err = do_estimator(dpd);
    }

    gretl_matrix_free(Vcpy);

    if (err) {
	fprintf(stderr, "step 2: dpd_calculate returning %d\n", err);
    }

    return err;
}

/* allocate storage needed by do_units() and friends */

static int dpd_allocate_new (dpdinfo *dpd, int north,
			     gretl_matrix **D, gretl_matrix **H,
			     gretl_matrix **Yi, gretl_matrix **Xi,
			     gretl_matrix **Zi)
{
    int err = 0;

    dpd->XZ = gretl_zero_matrix_new(dpd->k, dpd->nz);
    dpd->ZZ = gretl_zero_matrix_new(dpd->nz, dpd->nz);
    dpd->ZY = gretl_zero_matrix_new(dpd->nz, 1);

    if (dpd->XZ == NULL || dpd->ZZ == NULL || dpd->ZY == NULL) {
	return E_ALLOC;
    }

    if (dpd->flags & DPD_DPDSTYLE) {
	/* Ox/DPD-style H matrix: D is not needed */
	*D = NULL;
    } else {
	*D = gretl_matrix_alloc(dpd->T, north);
	if (*D == NULL) {
	    return E_ALLOC;
	}
    }

    *H = gretl_identity_matrix_new(north);
    *Yi = gretl_matrix_alloc(1, north);
    *Xi = gretl_matrix_alloc(dpd->k, north);
    *Zi = gretl_matrix_alloc(dpd->nz, north);

    if (*H == NULL || *Yi == NULL || *Xi == NULL || *Zi == NULL) {
	err = E_ALLOC;
    }
    
    if (!err) {
	int bigN = dpd->N * north;

	/* holders for stacking the per-unit data matrices */
	dpd->Y = gretl_zero_matrix_new(bigN, 1);
	dpd->X = gretl_zero_matrix_new(bigN, dpd->k);
	dpd->ZT = gretl_zero_matrix_new(dpd->nz, bigN);

	if (dpd->Y == NULL || dpd->X == NULL || dpd->ZT == NULL) {
	    err = E_ALLOC;
	}
    } 

    return err;
}

/* Stack the per-unit data matrices from unit @unum for future use, 
   skipping unused observations and recording the numbers of 
   observations in differences and in levels.
*/

static void stack_Y_and_X (dpdinfo *dpd,
			   const gretl_matrix *Yi, 
			   const gretl_matrix *Xi,
			   int *goodobs, int unum,
			   int *row)
{
    double xij;
    int i, j, k, t = *row;

    for (i=2; i<=goodobs[0]; i++) {
	k = goodobs[i] - 2;
	gretl_vector_set(dpd->Y, t, Yi->val[k]);
	for (j=0; j<Xi->rows; j++) {
	    xij = gretl_matrix_get(Xi, j, k);
	    gretl_matrix_set(dpd->X, t, j, xij);
	}
	t++;
    }

    /* record the number of differenced obs in t1 */
    dpd->ui[unum].t1 = goodobs[0] - 1;

    if (use_levels(dpd)) {
	int k0 = dpd->T - dpd->p - 1;

	for (i=1; i<=goodobs[0]; i++) {
	    k = k0 + goodobs[i] - 1;
	    gretl_vector_set(dpd->Y, t, Yi->val[k]);
	    for (j=0; j<Xi->rows; j++) {
		xij = gretl_matrix_get(Xi, j, k);
		gretl_matrix_set(dpd->X, t, j, xij);
	    }
	    t++;
	}
	/* record the number of levels obs in t2 */
	dpd->ui[unum].t2 = goodobs[0];
    }

    dpd->ui[unum].nobs = dpd->ui[unum].t1 + dpd->ui[unum].t2;

    *row = t;
}

/* Main driver for system GMM: the core is a loop across
   the panel units to build the moment matrices.
*/

static int do_units (dpdinfo *dpd, const double **Z, 
		     const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *D = NULL;
    gretl_matrix *H = NULL;
    gretl_matrix *Yi = NULL;
    gretl_matrix *Xi = NULL;
    gretl_matrix *Zi = NULL;
    int *goodobs, Tshort;
    int Yrow, Zcol;
    int nz_diff, nz_lev = 0;
    int north, i, t;
    int err = 0;

    /* full T - maxlag - 1 */
    Tshort = dpd->T - dpd->p - 1;

    if (use_levels(dpd)) {
	north = Tshort * 2 + 1;
    } else {
	north = Tshort;
    }

    /* instruments specific to equations in differences */
    nz_diff = (dpd->T - 1) * (dpd->T - 2) / 2;

    /* and to equations in levels */
    if (use_levels(dpd)) {
	nz_lev = dpd->T - 1;
    } 

    /* total instruments */
    dpd->nz = nz_diff + nz_lev + dpd->nx;

    goodobs = gretl_list_new(dpd->T);
    if (goodobs == NULL) {
	return E_ALLOC;
    }

    err = dpd_allocate_new(dpd, north, &D, &H, &Yi, &Xi, &Zi);
    if (err) {
	goto bailout;
    }

    /* initialize observation counts */
    dpd->effN = dpd->ndiff = dpd->nlev = 0;
    dpd->minTi = dpd->T;

    /* initialize data stackers */
    Yrow = Zcol = 0;

    for (i=0, t=pdinfo->t1; i<dpd->N; i++, t+=dpd->T) {
	int Ti = check_unit_obs(dpd, goodobs, Z, t);

#if DPDEBUG > 0
	fprintf(stderr, "\n\nUnit %4d:", i + 1);
	fprintf(stderr, " usable obs = %4d:\n", Ti);
#endif

	if (Ti > 0) {
	    /* this unit is usable */
	    dpd->effN += 1;
	    dpd->ndiff += Ti;
	    if (use_levels(dpd)) {
		dpd->nlev += Ti + 1;
	    }
	    if (Ti > dpd->maxTi) {
		dpd->maxTi = Ti;
	    }
	    if (Ti < dpd->minTi) {
		dpd->minTi = Ti;
	    }
	    if (D != NULL) {
		build_unit_D_matrix(dpd, goodobs, D);
	    }
	    compute_H(H, Tshort, D);
#if DPDEBUG > 2
	    if (t == pdinfo->t1) {
		gretl_matrix_print(H, "H");
	    }
#endif
	    build_Y(dpd, goodobs, Z, t, Yi);
	    build_X(dpd, goodobs, Z, t, Xi);
	    build_Z(dpd, goodobs, Z, t, nz_diff, nz_lev, Zi);
#if DPDEBUG
	    gretl_matrix_print(Yi, "do_units: Yi");
	    gretl_matrix_print(Xi, "do_units: Xi");
	    gretl_matrix_print(Zi, "do_units: Zi");
#endif
	    gretl_matrix_multiply_mod(Xi, GRETL_MOD_NONE,
				      Zi, GRETL_MOD_TRANSPOSE,
				      dpd->XZ, GRETL_MOD_CUMULATE);
	    gretl_matrix_qform(Zi, GRETL_MOD_NONE,
			       H, dpd->ZZ, GRETL_MOD_CUMULATE);
	    gretl_matrix_multiply_mod(Zi, GRETL_MOD_NONE,
				      Yi, GRETL_MOD_TRANSPOSE,
				      dpd->ZY, GRETL_MOD_CUMULATE);

	    /* stack the individual data matrices for future use */
	    stack_Y_and_X(dpd, Yi, Xi, goodobs, i, &Yrow);
	    gretl_matrix_inscribe_matrix(dpd->ZT, Zi, 0, Zcol,
					 GRETL_MOD_NONE);
	    Zcol += Zi->cols;
	}
    }

    dpd->Y->rows = dpd->X->rows = Yrow;

#if DPDEBUG
    gretl_matrix_print(dpd->Y, "dpd->Y");
    gretl_matrix_print(dpd->X, "dpd->X");
#endif

 bailout:

    free(goodobs);

    gretl_matrix_free(D);
    gretl_matrix_free(H);
    gretl_matrix_free(Yi);
    gretl_matrix_free(Xi);

    dpd->Zi = Zi; /* useful later */

    return err;
}

/* recompute \hat{\beta} and its variance matrix */

static int dpanel_step_2 (dpdinfo *dpd)
{
    gretl_matrix *u1 = NULL;
    gretl_matrix *V1 = NULL;
    int err;

    err = do_estimator_2(dpd);

    if (!err && (dpd->flags & DPD_WINCORR)) {
	/* we'll need a copy of the step-1 residuals, and also
	   the step-1 var(\hat{\beta})
	*/
	u1 = gretl_matrix_copy(dpd->uhat);
	V1 = gretl_matrix_copy(dpd->vbeta);
	if (u1 == NULL || V1 == NULL) {
	    err = E_ALLOC;
	}
    }		
	
    if (!err) {
	err = dpanel_add_residuals(dpd);
    }

    if (!err) {
	err = dpanel_step2_variance(dpd, u1, V1);
    }

    gretl_matrix_free(u1);
    gretl_matrix_free(V1);

    return err;
}

/* "Secret" public interface for new approach, including system GMM:
   in a script use --system (OPT_L, think "levels") to get
   Blundell-Bond. To build the H matrix as per Ox/DPD use the
   option --dpdstyle (OPT_X).

   The dpdinfo struct that is used here is a hybrid. It contains some
   useful common information, filled out by dpdinfo_new; it also has a
   fair number of members that are used by arbond but not (currently)
   by dpanel.  And it contains a few members that are specific to
   dpanel. This will have to be tidied up eventually.

   When dpdinfo_new() and dpd_allocate_matrices() have returned
   successfully, we have the following information available in the
   dpd struct:

   dpd->T:       total observations per unit (= pdinfo->pd)
   dpd->N:       number of units in sample range
   dpd->yno:     ID number of dependent variable
   dpd->xlist:   list of exogenous vars
   dpd->nx:      number of exogenous variables (= dpd->xlist[0])
   dpd->laglist: list of lags of y (faked up at present)
   dpd->p:       (max) lag order for y
   dpd->k:       total number of parameters

   and the following pre-allocated matrices:

   dpd->beta:    dpd->k x 1
   dpd->vbeta:   dpd->k x dpd->k (identity matrix, so far)

   We have the following things still to be set:

   dpd->ndiff: number of observations in differences
   dpd->nlev:  number of observations in levels
   dpd->nobs:  total number of observations actually used
   dpd->effN:  number of units actually used
   dpd->maxTi: the max number of equations in differences
               for any unit
   dpd->minTi: the min number of equations in differences
               for any unit actually used
   dpd->nz:    total number of instruments
   dpd->XZ:    cross-moment matrix (not allocated yet)
   dpd->ZZ:    ditto
   dpd->ZY:    ditto
   dpd->ZT:    all Zs, stacked, not allocated yet

*/

MODEL dpd_estimate (const int *list, const char *istr, 
		    const double **Z, const DATAINFO *pdinfo, 
		    gretlopt opt, PRN *prn)
{
    struct diag_info *d = NULL;
    dpdinfo *dpd = NULL;
    int nzb = 0;
    MODEL mod;
    int err = 0;

    gretl_model_init(&mod);
    gretl_model_smpl_init(&mod, pdinfo);

    dpd = dpdinfo_new(list, Z, pdinfo, opt | OPT_B, d, nzb, &mod.errcode);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in dpd_init\n", mod.errcode);
	return mod;
    }

    dpd_allocate_matrices(dpd);

    /* additional info that should be folded into dpdinfo
       at some point */

    dpd->used = calloc(pdinfo->n, 1);
    if (dpd->used == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	/* build the moment matrices */
	err = do_units(dpd, Z, pdinfo, prn);
    }

#if DPDEBUG > 0
    fprintf(stderr, "dpd_estimate, after do_units:\n");
    gretl_matrix_print(dpd->XZ, "XZ");
    gretl_matrix_print(dpd->ZZ, "ZZ");
    gretl_matrix_print(dpd->ZY, "ZY");
#endif

    if (!err) {
	/* calculate \hat{\beta} */
	do_estimator(dpd);
    }

    if (!err) {
	/* add residuals in differences */
	err = dpanel_add_residuals(dpd);
    }

    if (!err) {
	/* calculate variance of \hat{\beta} */
	err = dpanel_step1_variance(dpd);
    }

    if (!err && (opt & OPT_T)) {
	/* second step, if wanted */
	err = dpanel_step_2(dpd);
    }	

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	/* write estimation info into model struct */
	mod.errcode = dpd_finalize_model(&mod, dpd, list, istr, 
					 Z[dpd->yno], pdinfo, opt);
    }

    /* FIXME these deallocations should be handled by dpdinfo_free(),
       eventually, but right now that's awkward because of the way
       gretl_matrix_block is used in dpdinfo.
    */
    gretl_matrix_free(dpd->uhat);
    gretl_matrix_free(dpd->XZ);
    gretl_matrix_free(dpd->ZY);
    gretl_matrix_free(dpd->ZZ);
    gretl_matrix_free(dpd->Y);
    gretl_matrix_free(dpd->X);
    gretl_matrix_free(dpd->ZT);
    gretl_matrix_free(dpd->den);
    gretl_matrix_free(dpd->A);
    gretl_matrix_free(dpd->XZA);
    gretl_matrix_free(dpd->kmtmp);
    gretl_matrix_free(dpd->R1);
    gretl_matrix_free(dpd->Zi);

    dpdinfo_free(dpd);

    return mod;
}
