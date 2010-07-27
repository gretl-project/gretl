/* 
   Holder (possibly temporary) for the code to support the new
   "dpanel" command, which is intended to generalize the old "arbond"
   to handle system GMM.

   Added by Allin, 2010-07-09; initial content was from Jack's 
   prototype code, newmask.c.
*/

#define DPDEBUG 0

#define use_levels(d) (d->flags & DPD_SYSTEM)

/* Populate the residual vector, dpd->uhat. In the system case
   we stack the residuals in levels under the residuals in
   differences, per unit. Calculate SSR and \sigma^2 while
   we're at it.
*/

static int dpanel_residuals (dpdinfo *dpd)
{
    const double *b = dpd->beta->val;
    double SSRd = 0.0, SSRl = 0.0;
    double x, ut;
    int i, j, k, kd, t;

    if (use_levels(dpd) && dpd->uhatd == NULL) {
	/* resids in differences will be needed for AR tests */
	dpd->uhatd = gretl_column_vector_alloc(dpd->ndiff);
	if (dpd->uhatd == NULL) {
	    return E_ALLOC;
	}
    }

    k = kd = 0;

    for (i=0; i<dpd->N; i++) {
	int ndiff = dpd->ui[i].nobs - dpd->ui[i].nlev;

	for (t=0; t<ndiff; t++) {
	    /* differences */
	    ut = dpd->Y->val[k];
	    for (j=0; j<dpd->k; j++) {
		x = gretl_matrix_get(dpd->X, k, j);
		ut -= b[j] * x;
	    }
	    SSRd += ut * ut;
	    dpd->uhat->val[k++] = ut;
	    if (dpd->uhatd != NULL) {
		dpd->uhatd->val[kd++] = ut;
	    }
	}
	for (t=0; t<dpd->ui[i].nlev; t++) {
	    /* levels, it applicable */
	    ut = dpd->Y->val[k];
	    for (j=0; j<dpd->k; j++) {
		x = gretl_matrix_get(dpd->X, k, j);
		ut -= b[j] * x;
	    }
	    SSRl += ut * ut;
	    dpd->uhat->val[k++] = ut;
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
	    if (goodobs[0] > 1) {
		dpd->used[big_t] = 1;
	    } else if (use_levels(dpd)) {
		dpd->used[big_t] = 2;
	    }
	}
    }
    
    /* allow for differencing */
    return goodobs[0] - 1;
}

/* We do this accounting first so that we know the sizes of the
   various (and numerous) matrices that we'll need for the whole
   analysis at the outset; we can then allocate memory en bloc.
*/

static void do_unit_accounting (dpdinfo *dpd, const double **Z,
				int **Goodobs)
{
    int i, t;

    /* instruments specific to equations in differences */
    dpd->nzdiff = (dpd->T - 1) * (dpd->T - 2) / 2;

    /* and to equations in levels */
    if (use_levels(dpd)) {
	dpd->nzlev = dpd->T - 1;
    } 

    /* total instruments */
    dpd->nz = dpd->nzdiff + dpd->nzlev + dpd->nx;

    /* initialize observation counts */
    dpd->effN = dpd->ndiff = dpd->nlev = 0;
    dpd->minTi = dpd->T;

    for (i=0, t=dpd->t1; i<dpd->N; i++, t+=dpd->T) {
	int *goodobs = Goodobs[i];
	int Ti = check_unit_obs(dpd, goodobs, Z, t);

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
	}
    }

    /* figure the max number of stacked observations per unit */
    dpd->max_ni = dpd->maxTi;
    if (dpd->flags & DPD_SYSTEM) {
	dpd->max_ni += dpd->maxTi + 1;
    }  

    /* and the total observations overall */
    dpd->totobs = dpd->ndiff + dpd->nlev;
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

static void build_unit_H_matrix (dpdinfo *dpd, int *goodobs, 
				 gretl_matrix *D)
{
    build_unit_D_matrix(dpd, goodobs, D);
    gretl_matrix_multiply_mod(D, GRETL_MOD_TRANSPOSE, 
			      D, GRETL_MOD_NONE, 
			      dpd->H, GRETL_MOD_NONE);
}

static void make_dpdstyle_H (gretl_matrix *H, int nd)
{
    int i;

    gretl_matrix_zero(H);
    gretl_matrix_set(H, 0, 0, 2);

    for (i=1; i<H->rows; i++) {
	if (i < nd) {
	    gretl_matrix_set(H, i, i, 2);
	    gretl_matrix_set(H, i-1, i, -1);
	    gretl_matrix_set(H, i, i-1, -1);
	} else {
	    gretl_matrix_set(H, i, i, 1);
	}
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
		     int t, gretl_matrix *Zi)
{
    const double *y = Z[dpd->yno];
    int usable = goodobs[0] - 1;
    int maxlag = dpd->p;
    int T = dpd->T;
    const double *xj;
    double dx, y0;
    int t0, t1, i0, i1;
    int k, k2 = dpd->nzdiff + dpd->nzlev;
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

static void trim_zero_inst (dpdinfo *dpd)
{
    int i, n = dpd->A->rows;
    int trim = 0;

    for (i=0; i<n; i++) {
	if (gretl_matrix_get(dpd->A, i, i) == 0.0) {
	    trim = 1;
	    break;
	}
    }

    if (trim) {
	char *mask = calloc(n, 1);

	for (i=0; i<n; i++) {
	    if (gretl_matrix_get(dpd->A, i, i) == 0.0) {
		mask[i] = 1;
	    }
	}

	gretl_matrix_cut_rows_cols(dpd->A, mask);
	gretl_matrix_cut_cols(dpd->XZ, mask);
	gretl_matrix_cut_rows(dpd->ZY, mask);
	gretl_matrix_cut_rows(dpd->ZT, mask);

	dpd->nz = dpd->A->rows;

	/* adjust other matrices whose dimensions involve nz */
	gretl_matrix_reuse(dpd->Acpy,  dpd->nz, dpd->nz);
	gretl_matrix_reuse(dpd->kmtmp, -1, dpd->nz);
	gretl_matrix_reuse(dpd->L1,    -1, dpd->nz);
	gretl_matrix_reuse(dpd->XZA,   -1, dpd->nz);

	free(mask);
    }
}

/* Here we compute \hat{\beta} from the moment matrices. */

static int do_estimator (dpdinfo *dpd)
{
    gretl_matrix *M1 = NULL, *M2 = NULL;
    int k, err = 0;

    if (dpd->step == 1) {
	trim_zero_inst(dpd);
#if DPDEBUG > 1
	gretl_matrix_print(dpd->A, "A (trimmed)");
	gretl_matrix_print(dpd->XZ, "XZ (trimmed)");
	gretl_matrix_print(dpd->ZY, "ZY (trimmed)");
#endif
    }

    k = dpd->XZ->rows;

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
	/* dpd->XZA <- XZ * A */
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
	dpd->M = M1;
    } else {
	gretl_matrix_free(M1);
    }

    return err;
}

static int do_estimator_2 (dpdinfo *dpd)
{
    int err = 0;

#if DPDEBUG
    gretl_matrix_print(dpd->V, "V, in do_estimator_2");
#endif

    if (dpd->V->rows > dpd->effN) {
	/* we know this case requires special treatment */
	err = gretl_SVD_invert_matrix(dpd->V);
	if (!err) {
	    gretl_matrix_xtr_symmetric(dpd->V);
	}
    } else {
	gretl_matrix_copy_values(dpd->Acpy, dpd->V);
	err = gretl_invert_symmetric_matrix(dpd->V);
	if (err) {
	    /* revert the data in dpd->V */
	    gretl_matrix_copy_values(dpd->V, dpd->Acpy);
	    err = gretl_SVD_invert_matrix(dpd->V);
	    if (!err) {
		gretl_matrix_xtr_symmetric(dpd->V);
	    }
	}
    }

    if (!err) {
	/* A <- V^{-1} */
	gretl_matrix_copy_values(dpd->A, dpd->V);
	dpd->step = 2;
	err = do_estimator(dpd);
    }

    if (err) {
	fprintf(stderr, "step 2: do_estimator_2 returning %d\n", err);
    }

    return err;
}

/* allocate temporary storage needed by do_units() */

static int make_units_workspace (dpdinfo *dpd, gretl_matrix **D, 
				 gretl_matrix **Yi, gretl_matrix **Xi)
{
    int err = 0;

    if (dpd->flags & DPD_DPDSTYLE) {
	/* Ox/DPD-style H matrix: D is not needed */
	*D = NULL;
    } else {
	*D = gretl_matrix_alloc(dpd->T, dpd->max_ni);
	if (*D == NULL) {
	    return E_ALLOC;
	}
    }

    *Yi = gretl_matrix_alloc(1, dpd->max_ni);
    *Xi = gretl_matrix_alloc(dpd->k, dpd->max_ni);

    if (*Yi == NULL || *Xi == NULL) {
	err = E_ALLOC;
    }
    
    return err;
}

/* Stack the per-unit data matrices from unit @unum for future use, 
   skipping unused observations and recording the numbers of 
   observations in differences and in levels.
*/

static void stack_unit_data (dpdinfo *dpd,
			     const gretl_matrix *Yi, 
			     const gretl_matrix *Xi,
			     const gretl_matrix *Zi,
			     int *goodobs, int unum,
			     int *row)
{
    double x;
    int i, j, k, t = *row;

    for (i=2; i<=goodobs[0]; i++) {
	k = goodobs[i] - 2;
	gretl_vector_set(dpd->Y, t, Yi->val[k]);
	for (j=0; j<Xi->rows; j++) {
	    x = gretl_matrix_get(Xi, j, k);
	    gretl_matrix_set(dpd->X, t, j, x);
	}
	for (j=0; j<dpd->nz; j++) {
	    x = gretl_matrix_get(Zi, j, k);
	    gretl_matrix_set(dpd->ZT, j, t, x);
	}
	t++;
    }

    /* record the indices of the first and last
       differenced observations */
    dpd->ui[unum].t1 = goodobs[2];
    dpd->ui[unum].t2 = goodobs[goodobs[0]];

    /* record the number of differenced obs */
    dpd->ui[unum].nobs = goodobs[0] - 1;

    if (use_levels(dpd)) {
	int k0 = dpd->T - dpd->p - 1;

	for (i=1; i<=goodobs[0]; i++) {
	    k = k0 + goodobs[i] - 1;
	    gretl_vector_set(dpd->Y, t, Yi->val[k]);
	    for (j=0; j<Xi->rows; j++) {
		x = gretl_matrix_get(Xi, j, k);
		gretl_matrix_set(dpd->X, t, j, x);
	    }
	    for (j=0; j<dpd->nz; j++) {
		x = gretl_matrix_get(Zi, j, k);
		gretl_matrix_set(dpd->ZT, j, t, x);
	    }
	    t++;
	}
	/* record the number of levels obs and augment total */
	dpd->ui[unum].nlev = goodobs[0];
	dpd->ui[unum].nobs += dpd->ui[unum].nlev;
    }

    *row = t;
}

/* Trims the last @ntrim rows from matrix @m, in place; note 
   that you cannot do this just by adusting m->rows if @m has 
   more than one column.
*/

void matrix_trim_rows (gretl_matrix *m, int ntrim)
{
    if (m->cols > 1) {
	double *dest, *src;
	int k = m->rows - ntrim;
	int n = k * sizeof *src;
	int j;

	dest = m->val + k;
	src = m->val + m->rows;

	for (j=1; j<m->cols; j++) {
	    memmove(dest, src, n);
	    dest += k;
	    src += m->rows;
	}
    }

    m->rows -= ntrim;
}

/* Main driver for system GMM: the core is a loop across
   the panel units to build the moment matrices.
*/

static int do_units (dpdinfo *dpd, const double **Z, 
		     int **Goodobs)
{
    gretl_matrix *D = NULL;
    gretl_matrix *Yi = NULL;
    gretl_matrix *Xi = NULL;
    gretl_matrix *Zi = NULL;
    int i, t, Yrow;
    int err = 0;

    err = make_units_workspace(dpd, &D, &Yi, &Xi);
    if (err) {
	goto bailout;
    }

    Zi = dpd->Zi;
    gretl_matrix_reuse(Zi, dpd->nz, dpd->max_ni);
    
    if (D == NULL) {
	/* H will not vary by unit */
	make_dpdstyle_H(dpd->H, dpd->maxTi);
    }

    /* initialize matrices to cumulate */
    gretl_matrix_zero(dpd->XZ);
    gretl_matrix_zero(dpd->A);
    gretl_matrix_zero(dpd->ZY);

    /* initialize data stacker */
    Yrow = 0;

    for (i=0, t=dpd->t1; i<dpd->N; i++, t+=dpd->T) {
	int *goodobs = Goodobs[i];
	int Ti = goodobs[0] - 1;

#if DPDEBUG
	fprintf(stderr, "\n\nUnit %d, usable obs %d:", i + 1, Ti);
#endif

	if (Ti > 0) {
	    /* this unit is usable */
	    if (D != NULL) {
		build_unit_H_matrix(dpd, goodobs, D);
	    }
	    build_Y(dpd, goodobs, Z, t, Yi);
	    build_X(dpd, goodobs, Z, t, Xi);
	    build_Z(dpd, goodobs, Z, t, Zi);
#if DPDEBUG
	    gretl_matrix_print(Yi, "do_units: Yi");
	    gretl_matrix_print(Xi, "do_units: Xi");
	    gretl_matrix_print(Zi, "do_units: Zi");
#endif
	    gretl_matrix_multiply_mod(Xi, GRETL_MOD_NONE,
				      Zi, GRETL_MOD_TRANSPOSE,
				      dpd->XZ, GRETL_MOD_CUMULATE);
	    gretl_matrix_qform(Zi, GRETL_MOD_NONE,
			       dpd->H, dpd->A, GRETL_MOD_CUMULATE);
	    gretl_matrix_multiply_mod(Zi, GRETL_MOD_NONE,
				      Yi, GRETL_MOD_TRANSPOSE,
				      dpd->ZY, GRETL_MOD_CUMULATE);
	    /* stack the individual data matrices for future use */
	    stack_unit_data(dpd, Yi, Xi, Zi, goodobs, i, &Yrow);
	}
    }

    if (Yrow < dpd->Y->rows) {
	fprintf(stderr, "*** obs TRIMMING is needed *** \n");
	dpd->Y->rows = dpd->ZT->cols = Yrow;
	matrix_trim_rows(dpd->X, dpd->X->rows - Yrow);
    }

#if DPDEBUG
    gretl_matrix_print(dpd->Y, "dpd->Y");
    gretl_matrix_print(dpd->X, "dpd->X");
#endif

 bailout:

    gretl_matrix_free(D);
    gretl_matrix_free(Yi);
    gretl_matrix_free(Xi);

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
	err = dpanel_residuals(dpd);
    }

    if (!err) {
	err = dpd_variance_2(dpd, u1, V1);
    }

    gretl_matrix_free(u1);
    gretl_matrix_free(V1);

    if (!err) {
#if 0 /* not just yet */
	ar_test(dpd, dpd->M);
#endif
	sargan_test(dpd);
	dpd_wald_test(dpd);
    }

    return err;
}

static int dpanel_step_1 (dpdinfo *dpd)
{
    int err;

    /* calculate \hat{\beta} */
    err = do_estimator(dpd);

    if (!err) {
	/* add residuals */
	err = dpanel_residuals(dpd);
    }

    if (!err) {
	/* calculate variance of \hat{\beta} */
	err = dpd_variance_1(dpd);
    }

    if (!err && !(dpd->flags & DPD_TWOSTEP)) {
#if 0 /* not just yet */
	ar_test(dpd, dpd->M);
#endif
	sargan_test(dpd);
	dpd_wald_test(dpd);
    }

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
   dpd->A:     ditto
   dpd->ZY:    ditto
   dpd->ZT:    all Zs, stacked, not allocated yet

*/

MODEL dpd_estimate (const int *list, const char *istr, 
		    const double **Z, const DATAINFO *pdinfo, 
		    gretlopt opt, PRN *prn)
{
    struct diag_info *d = NULL;
    dpdinfo *dpd = NULL;
    int **Goodobs = NULL;
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

    Goodobs = gretl_list_array_new(dpd->N, dpd->T);
    if (Goodobs == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	do_unit_accounting(dpd, Z, Goodobs);
	err = dpd_allocate_matrices(dpd);
    }

    if (!err) {
	/* build the moment matrices */
	err = do_units(dpd, Z, Goodobs);
    }

    gretl_list_array_free(Goodobs, dpd->N);

    if (!err) {
	err = dpanel_step_1(dpd);
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

    dpdinfo_free(dpd);

    return mod;
}
