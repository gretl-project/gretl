/* 
   Holder (possibly temporary) for the code to support the new
   "dpanel" command, which is intended to generalize the old "arbond"
   to handle system GMM.

   Added by Allin, 2010-07-09; contents are mostly from Jack's 
   prototype code, newmask.c.
*/

#define DPDEBUG 0

#define LEV_ONLY (-1) /* flag for an obs that's good only for levels */

#define use_levels(d) (d->flags & DPD_SYSTEM)

/* For a given unit, whose observations start at @t0: count the
   available residuals and while we're at it transcribe those
   residuals into the vector @ui. We place residuals from the
   differences equations first, followed by those from the levels
   equations in the "system" case.

   (If @ui is NULL we just count and don't transcribe.)
*/

static int n_available_residuals (dpdinfo *dpd, gretl_matrix *ui,
				  int t0, int *sd, int *sl)
{
    int t, tmax = t0 + dpd->T;
    int s, n = 0;

    /* count and transcribe the residuals from the differences
       equations */

    s = *sd;

    for (t=t0; t<tmax; t++) {
	if (dpd->used[t] > 0) {
	    if (ui != NULL) {
		gretl_vector_set(ui, n++, dpd->uhat->val[s++]);
	    } else {
		n++;
		s++;
	    }
	}
    }

    *sd = s;

    /* count and transcribe the residuals from the levels
       equations, if applicable */

    if (use_levels(dpd)) {
	s = *sl;
	for (t=t0; t<tmax; t++) {
	    if (dpd->used[t]) {
		if (ui != NULL) {
		    gretl_vector_set(ui, n++, dpd->uhl->val[s++]);
		} else {
		    n++;
		    s++;
		}
	    }
	}
	*sl = s;
    }	

    return n;
}

static int dpanel_step1_variance (dpdinfo *dpd, const DATAINFO *pdinfo,
				  const char *ZZmask)
{
    gretl_matrix_block *B;
    gretl_matrix *kk, *kz, *V;
    gretl_matrix *ui, *Zi, *uZ;
    int sd, sl, t, nz, max_ni;
    int Zcol = 0;
    int err = 0;

    nz = dpd->A->rows; /* note: may be trimmed relative to dpd->nz */

    sd = sl = max_ni = 0;

    for (t=pdinfo->t1; t<pdinfo->t2; t+=dpd->T) {
	int ni = n_available_residuals(dpd, NULL, t, &sd, &sl);

	if (ni > max_ni) {
	    max_ni = ni;
	}
    }

#if DPDEBUG
    fprintf(stderr, "nz = %d, max_ni = %d\n", nz, max_ni);
#endif

    B = gretl_matrix_block_new(&kk, dpd->k, dpd->k,
			       &kz, dpd->k, nz,
			       &V,  nz, nz,
			       &ui, max_ni, 1,
			       &Zi, max_ni, dpd->nz,
			       &uZ, 1, nz,
			       NULL);
    if (B == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_zero(V);

    sd = sl = 0;

    for (t=pdinfo->t1; t<pdinfo->t2; t+=dpd->T) {
	int ni = n_available_residuals(dpd, ui, t, &sd, &sl);

	if (ni == 0) {
	    continue;
	}

	gretl_matrix_reuse(Zi, ni, dpd->nz);
	gretl_matrix_reuse(ui, ni, -1);

	/* Extract the saved Zi from the stack in dpd->ZT:
	   it will need trimming if the cross-moment matrices
	   have been trimmed, which is recorded in ZZmask.
	*/

	gretl_matrix_extract_matrix(Zi, dpd->ZT, 0, Zcol,
				    GRETL_MOD_TRANSPOSE);
	Zcol += Zi->rows;

	if (ZZmask != NULL) {
	    gretl_matrix_cut_cols(Zi, ZZmask);
	}

#if DPDEBUG
	gretl_matrix_print(ui, "ui (vcalc)");
	gretl_matrix_print(Zi, "Zi (vcalc)");
#endif

	err = gretl_matrix_multiply_mod(ui, GRETL_MOD_TRANSPOSE,
					Zi, GRETL_MOD_NONE,
					uZ, GRETL_MOD_NONE);
	if (!err) {
	    err = gretl_matrix_multiply_mod(uZ, GRETL_MOD_TRANSPOSE,
					    uZ, GRETL_MOD_NONE,
					    V, GRETL_MOD_CUMULATE);
	}

	if (err) {
	    fprintf(stderr, "dpanel_step1_variance: error at t=%d (ni = %d)\n", 
		    t, ni);
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
	err = gretl_matrix_qform(dpd->den, GRETL_MOD_NONE, kk, dpd->vbeta,
				 GRETL_MOD_NONE);
    }

    if (!err) {
	gretl_matrix_multiply_by_scalar(dpd->vbeta, dpd->effN);
    }

    gretl_matrix_block_destroy(B);

    if (err) {
	fprintf(stderr, "dpanel_step1_variance: err = %d (nz = %d)\n", err, nz);
    }

    return err;
}

static int 
dpanel_add_residuals (dpdinfo *dpd, const char *d_order,
		      const double **Z, const DATAINFO *pdinfo)
{
    const double *b = dpd->beta->val;
    const double *x;
    double ut, dx;
    int i, j, k, li, s, t;

    if (dpd->uhat == NULL) {
	dpd->uhat = gretl_column_vector_alloc(dpd->nobs);
	if (dpd->uhat == NULL) {
	    return E_ALLOC;
	}
    }

    dpd->SSR = 0.0;
    s = 0;

    /* calculate and store the residuals in differences */

    for (t=0; t<pdinfo->n; t++) {
	if (dpd->used[t] > 0) {
	    /* look up the lag needed to create a valid difference */
	    k = d_order[t];
	    j = 0;
	    ut = dpd->y[t] - dpd->y[t-k];
	    for (i=1; i<=dpd->laglist[0]; i++) {
		li = dpd->laglist[i];
		dx = dpd->y[t-li] - dpd->y[t-li-k];
		ut -= b[j++] * dx;
	    }
	    /* FIXME automatic time dummies */
	    if (dpd->nx > 0) {
		for (i=1; i<=dpd->xlist[0]; i++) {
		    /* preserve the constant */
		    if (!use_levels(dpd) && dpd->xlist[i] == 0) {
			ut -= b[j++];
		    } else {
			x = Z[dpd->xlist[i]];
			dx = x[t] - x[t-k];
			ut -= b[j++] * dx;
		    }
		    
		}
	    }
	    gretl_vector_set(dpd->uhat, s++, ut);
	    dpd->SSR += ut * ut;
	}
    }

    dpd->s2 = dpd->SSR / (dpd->nobs - dpd->k);

    /* residuals in levels: note that the values will agree with
       Ox/DPD only if the --dpdstyle option is given.
    */

    if (use_levels(dpd)) {
	double SSR_lev = 0.0;
	int nobs_lev = 0;

	dpd->uhl = gretl_column_vector_alloc(dpd->nobs + dpd->effN);
	if (dpd->uhl == NULL) {
	    return E_ALLOC;
	}

	/* These will be over-written, print them here for
	   reference */
	fprintf(stderr, "nobs (diffs) = %d\n", dpd->nobs);
	fprintf(stderr, "SSR (diffs) = %g\n", dpd->SSR);
	fprintf(stderr, "s^2 (diffs) = %g\n", dpd->s2);

	s = 0;

 	for (t=0; t<pdinfo->n; t++) {
	    if (dpd->used[t]) {
		j = 0;
		ut = dpd->y[t];
		for (i=1; i<=dpd->laglist[0]; i++) {
		    li = dpd->laglist[i];
		    ut -= b[j++] * dpd->y[t-li];
		}
		if (dpd->nx > 0) {
		    for (i=1; i<=dpd->xlist[0]; i++) {
			x = Z[dpd->xlist[i]];
			ut -= b[j++] * x[t];
		    }
		}
		gretl_vector_set(dpd->uhl, s++, ut);
		SSR_lev += ut * ut;
		nobs_lev++;
	    }
	}

	dpd->nobs = nobs_lev;
	dpd->SSR = SSR_lev;
	dpd->s2 = dpd->SSR / (dpd->nobs - dpd->k);
    }    

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
    int i, s, t, ok;

    goodobs[0] = 0;

    for (t=0; t<dpd->T; t++) {
	int big_t = t + t0;

	/* do we have the dependent variable? */
	ok = !na(dpd->y[big_t]);

	/* lags of dependent variable? */
	for (i=1; i<=dpd->laglist[0] && ok; i++) {
	    s = t - dpd->laglist[i];
	    if (s < 0) {
		ok = 0;
	    } else {
		ok &= !na(dpd->y[s+t0]);
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
	    if (goodobs[0] > 1) {
		dpd->used[big_t] = 1;
	    } else {
		dpd->used[big_t] = LEV_ONLY;
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
   differences, followed by levels if wanted. As we go, we write into
   @d_order the order of differencing at each observation (allowing
   for the possibility that we need a difference of order greater than
   1 due to missing values). This information will be needed when
   calculating the residuals from the equations in differences.
*/

static void build_Y (dpdinfo *dpd, int *goodobs, const double **Z, 
		     int t, gretl_matrix *Yi, char *d_order)
{
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
	dy = dpd->y[t1] - dpd->y[t0];
	d_order[t1] = t1 - t0;
	gretl_vector_set(Yi, i1-1-maxlag, dy);
    }
    
    if (use_levels(dpd)) {
	/* levels */
	for (i=0; i<=usable; i++) {
	    i1 = goodobs[i+1];
	    t1 = t + i1;
	    gretl_vector_set(Yi, Tshort + (i1-maxlag), dpd->y[t1]);
	}
    }
}

/* Build matrix of right-hand side variable values in @Xi using
   differences, followed by levels if wanted.
*/

static void build_X (dpdinfo *dpd, int *goodobs, const double **Z, 
		     int t, gretl_matrix *Xi)
{
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
	    dx = dpd->y[t1-lj] - dpd->y[t0-lj];
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
		gretl_matrix_set(Xi, j-1, col, dpd->y[t1-lj]);
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
		y0 = dpd->y[t+j];
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
		    y0 = (j < 1)? NADBL : dpd->y[t+j-1];
		    y1 = dpd->y[t+j];
		    if (!na(y1) && !na(y0)) {
			gretl_matrix_set(Zi, row, col, y1 - y0);
		    }
		}
		row++;
	    }
	    lastdiff = j;
	}

	/* equations in levels: exog vars */
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

/* Here we compute \hat{\beta} from the moment matrices. */

static int do_estimator (dpdinfo *dpd, char **ZZmask)
{
    gretl_matrix *iZZ, *M, *M1, *M2;
    int k, m, err = 0;
    
    trim_zero_inst(dpd->XZ, dpd->ZZ, dpd->ZY, ZZmask);

#if DPDEBUG > 0
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
	/* symmetrize ZZ copy; you never know */
	gretl_matrix_xtr_symmetric(iZZ);
	err = gretl_invert_symmetric_matrix(iZZ);
    }

    if (!err) {
	/* M1 <- XZ * ZZ^{-1} * XZ' */
	gretl_matrix_qform(dpd->XZ, GRETL_MOD_NONE,
			   iZZ, M1, GRETL_MOD_NONE);
    }

    if (!err) {
	/* M1 <- (XZ * ZZ^{-1} * XZ')^{1} */
	err = gretl_invert_symmetric_matrix(M1);
    }

    if (!err) {
	/* M <- XZ * ZZ^{-1} */
	gretl_matrix_multiply(dpd->XZ, iZZ, M);
	/* M2 <- XZ * ZZ^{-1} * ZY */
	gretl_matrix_multiply(M, dpd->ZY, M2);
	/* beta <- (XZ * ZZ^{-1} * XZ')^{1} * (XZ * ZZ^{-1} * ZY) */
	gretl_matrix_multiply(M1, M2, dpd->beta);
    }

#if DPDEBUG > 0
    gretl_matrix_print(M1, "M1");
    gretl_matrix_print(M2, "M2");
#endif

    gretl_matrix_free(M);
    gretl_matrix_free(M2);

    if (!err) {
	/* save the matrices needed for computing the
	   variance of \hat{\beta} */
	dpd->den = M1;
	dpd->A = iZZ;
	gretl_matrix_divide_by_scalar(dpd->den, dpd->N);
	gretl_matrix_multiply_by_scalar(dpd->A, dpd->N);
    } else {
	gretl_matrix_free(M1);
	gretl_matrix_free(iZZ);
    }

    return err;
}

/* allocate storage needed by do_units() and friends */

static int dpd_allocate_new (dpdinfo *dpd, int north,
			     gretl_matrix **D, gretl_matrix **H,
			     gretl_matrix **Yi, gretl_matrix **Xi,
			     gretl_matrix **Zi)
{
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
	return E_ALLOC;
    }

    return 0;
}

/* Main driver for system GMM: the core is a loop across
   the panel units to build the moment matrices.
*/

static int do_units (dpdinfo *dpd, char *d_order,
		     const double **Z, const DATAINFO *pdinfo, 
		     PRN *prn)
{
    gretl_matrix *D = NULL;
    gretl_matrix *H = NULL;
    gretl_matrix *Yi = NULL;
    gretl_matrix *Xi = NULL;
    gretl_matrix *Zi = NULL;
    int *goodobs;
    int Tshort, Zcol = 0;
    int nz_diff, nz_lev = 0;
    int north, t, unit = 0;
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

#if DPDEBUG > 0
    pprintf(prn, "nz = %4d:\n", dpd->nz);
#endif

    goodobs = gretl_list_new(dpd->T);
    if (goodobs == NULL) {
	return E_ALLOC;
    }

    err = dpd_allocate_new(dpd, north, &D, &H, &Yi, &Xi, &Zi);
    if (err) {
	goto bailout;
    }

    /* holder for stacking the per-unit Zi matrices */
    dpd->ZT = gretl_zero_matrix_new(Zi->rows, dpd->N * north);
    if (dpd->ZT == NULL) {
	err = E_ALLOC;
	goto bailout;
    }    

    /* initialize observation counts */
    dpd->effN = dpd->nobs = 0;
    dpd->minTi = dpd->T;

    for (t=pdinfo->t1; t<pdinfo->t2; t+=dpd->T) {
	int Ti = check_unit_obs(dpd, goodobs, Z, t);

#if DPDEBUG > 0
	fprintf(stderr, "\n\nUnit %4d:", unit);
	fprintf(stderr, " usable obs = %4d:\n", Ti);
#endif

	if (Ti > 0) {
	    /* this unit is usable */
	    dpd->effN += 1;
	    dpd->nobs += Ti;
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
	    build_Y(dpd, goodobs, Z, t, Yi, d_order);
	    build_X(dpd, goodobs, Z, t, Xi);
	    build_Z(dpd, goodobs, Z, t, nz_diff, nz_lev, Zi);
#if DPDEBUG > 1
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
	    /* stack the individual Zi's for future use */
#if DPDEBUG == 1
	    gretl_matrix_print(Zi, "Zi, in do units");
#endif
	    gretl_matrix_inscribe_matrix(dpd->ZT, Zi, 0, Zcol,
					 GRETL_MOD_NONE);
	    Zcol += Zi->cols;
	}

	unit++;
    }

 bailout:

    gretl_matrix_free(D);
    gretl_matrix_free(H);
    gretl_matrix_free(Yi);
    gretl_matrix_free(Xi);
    gretl_matrix_free(Zi);

    free(goodobs);

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
   dpd->y:       convenience pointer to dep. var. in dataset
   dpd->xlist:   list of exogenous vars
   dpd->nx:      number of exogenous variables (= dpd->xlist[0])
   dpd->laglist: list of lags of y (faked up at present)
   dpd->p:       (max) lag order for y
   dpd->k:       total number of parameters

   and the following pre-allocated matrices:

   dpd->beta:    dpd->k x 1
   dpd->vbeta:   dpd->k x dpd->k (identity matrix, so far)

   We have the following things still to be set:

   dpd->nobs: total number of observations actually used
   dpd->effN: number of units actually used
   dpd->maxTi: the max number of equations in differences
               for any unit
   dpd->minTi: the min number of equations in differences
               for any unit actually used
   dpd->nz:   total number of instruments
   dpd->XZ:   cross-moment matrix (not allocated yet)
   dpd->ZZ:   ditto
   dpd->ZY:   ditto
   dpd->ZT:   all Zs, stacked, not allocated yet

*/

MODEL dpd_estimate (const int *list, const char *istr, 
		    const double **Z, const DATAINFO *pdinfo, 
		    gretlopt opt, PRN *prn)
{
    struct diag_info *d = NULL;
    dpdinfo *dpd = NULL;
    char *ZZmask = NULL;
    char *d_order = NULL;
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

    d_order = calloc(pdinfo->n, 1);
    if (d_order == NULL) {
	err = E_ALLOC;
    }     

    if (!err) {
	/* build the moment matrices */
	err = do_units(dpd, d_order, Z, pdinfo, prn);
    }

#if DPDEBUG > 0
    fprintf(stderr, "dpd_estimate, after do_units:\n");
    gretl_matrix_print(dpd->XZ, "XZ");
    gretl_matrix_print(dpd->ZZ, "ZZ");
    gretl_matrix_print(dpd->ZY, "ZY");
#endif

    if (!err) {
	/* calculate \hat{\beta} */
	do_estimator(dpd, &ZZmask);
    }

    if (!err) {
	/* add residuals in differences */
	err = dpanel_add_residuals(dpd, d_order, Z, pdinfo);
    }

    if (!err) {
	/* calculate variance of \hat{\beta} */
	err = dpanel_step1_variance(dpd, pdinfo, ZZmask);
    }

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	/* write estimation info into model struct */
	mod.errcode = dpd_finalize_model(&mod, dpd, list, istr, 
					 pdinfo, opt);
    }

    /* FIXME these deallocations should be handled by dpdinfo_free(),
       eventually, but right now that's awkward because of the way
       gretl_matrix_block is used in dpdinfo.
    */
    gretl_matrix_free(dpd->uhat);
    gretl_matrix_free(dpd->XZ);
    gretl_matrix_free(dpd->ZY);
    gretl_matrix_free(dpd->ZZ);
    gretl_matrix_free(dpd->ZT);
    gretl_matrix_free(dpd->den);
    gretl_matrix_free(dpd->A);
    gretl_matrix_free(dpd->uhl);

    dpdinfo_free(dpd);
    free(ZZmask);
    free(d_order);

    return mod;
}
