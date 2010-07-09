/* 
   Holder (possibly temporary) for the code to support the new
   "dpanel" command, which is intended to generalize the old "arbond"
   to handle system GMM.

   Added by Allin, 2010-07-09; contents are mostly from Jack's 
   prototype code, newmask.c.
*/

#define DPDSTYLE 0
#define DPDEBUG 0

#define LEV_ONLY (-1) /* flag for an obs that's good only for levels */

static int dpd_add_residual_vector (dpdinfo *dpd, const double **Z, 
				    const DATAINFO *pdinfo)
{
    const double *b = dpd->beta->val;
    const double *x;
    double ut, dx;
    int i, j, li, s, t;

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
	if (dpd->used[t] && dpd->used[t] != LEV_ONLY) {
	    j = 0;
	    ut = dpd->y[t] - dpd->y[t-1];
	    for (i=1; i<=dpd->laglist[0]; i++) {
		li = dpd->laglist[i];
		dx = dpd->y[t-li] - dpd->y[t-li-1];
		ut -= b[j++] * dx;
		
	    }
	    /* FIXME automatic time dummies */
	    if (dpd->nx > 0) {
		for (i=1; i<=dpd->xlist[0]; i++) {
		    x = Z[dpd->xlist[i]];
		    dx = x[t] - x[t-1];
		    ut -= b[j++] * dx;
		    
		}
	    }
	    gretl_vector_set(dpd->uhat, s++, ut);
	    dpd->SSR += ut * ut;
	}
    }

    dpd->s2 = dpd->SSR / (dpd->nobs - dpd->k);

    /* residuals in levels: for the present these are not used
       subsequently, we just print their summary stats for
       comparison with Ox/DPD. But note that the values will
       agree only if DPDSTYLE is turned on.
    */

    if (dpd->flags & DPD_SYSTEM) {
	double SSR_lev = 0.0;
	int nobs_lev = 0;

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
		SSR_lev += ut * ut;
		nobs_lev++;
	    }

	    /* for now, clear the "used" flag to avoid crashing 
	       when saving the residuals to the model */
	    if (dpd->used[t] == LEV_ONLY) {
		dpd->used[t] = 0;
	    }
	}

	fprintf(stderr, "SSR (levels) = %g\n", SSR_lev);
	fprintf(stderr, "s^2 (levels) = %g\n", SSR_lev / (nobs_lev - dpd->k));
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
   in the @goodobs list, fill matrix D (which may be used to
   construct H).
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
    if (dpd->flags & DPD_SYSTEM) {
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

static void compute_H (const gretl_matrix *D, gretl_matrix *H, 
		       int nh, int dpdstyle)
{
    if (dpdstyle) {
	int i;

	gretl_matrix_set(H, 0, 0, 2);
	for (i=1; i<nh; i++) {
	    gretl_matrix_set(H, i, i, 2);
	    gretl_matrix_set(H, i-1, i, -1);
	    gretl_matrix_set(H, i, i-1, -1);
	}
	gretl_matrix_multiply_by_scalar(H, 1); /* or not? */
    } else {
	gretl_matrix_multiply_mod(D, GRETL_MOD_TRANSPOSE, 
				  D, GRETL_MOD_NONE, 
				  H, GRETL_MOD_NONE);
    }
}

/* Build row vector of dependent variable values in @Yi using
   differences, followed by levels if wanted.
*/

static void build_Y (dpdinfo *dpd, int *goodobs, const double **Z, 
		     int t, gretl_matrix *Yi)
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
	gretl_vector_set(Yi, i1-1-maxlag, dy);
    }
    
    if (dpd->flags & DPD_SYSTEM) {
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
	    xj = Z[dpd->xlist[j]];
	    dx = xj[t1] - xj[t0];
	    gretl_matrix_set(Xi, j+nlags-1, i1-1-maxlag, dx);
	}
    }
    
    if (dpd->flags & DPD_SYSTEM) {
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
    int i, j;

    gretl_matrix_zero(Zi);

    /* equations in differences: lagged levels of y */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	k = i0 * (i0-1) / 2;
	for (j=0; j<i0; j++) {
	    y0 = dpd->y[t+j];
	    if (!na(y0)) {
		gretl_matrix_set(Zi, k, i1-1-maxlag, y0);
	    }
	    k++;
	}
    }

    /* equations in differences: differenced exog vars */
    if (dpd->nx > 0) {
	for (i=0; i<usable; i++) {
	    i0 = goodobs[i+1];
	    i1 = goodobs[i+2];
	    t0 = t + i0;
	    t1 = t + i1;
	    for (j=0; j<dpd->nx; j++) {
		xj = Z[dpd->xlist[j+1]];
		dx = xj[t1] - xj[t0];
		gretl_matrix_set(Zi, k2 + j, (i1-1-maxlag), dx);
	    }
	}
    }

    if (dpd->flags & DPD_SYSTEM) {
	int offset = T - maxlag - 1;
	int col, row = (T-1)*(T-2)/2;
	int lastdiff = 0;
	double y1;

	/* equations in levels: lagged differences of y */
	for (i=0; i<=usable; i++) {
	    i0 = goodobs[i+1];
	    col = offset + i0 - maxlag;
	    for (j=lastdiff; j<i0; j++) {
		y0 = (j < 1)? NADBL : dpd->y[t+j-1];
		y1 = dpd->y[t+j];
		if (!na(y1) && !na(y0)) {
		    gretl_matrix_set(Zi, row, col, y1 - y0);
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

/*
  here we construct \hat{\beta} from the moment matrices;
  no attempt is made to compute the covariance matrix
*/

static int do_estimator (dpdinfo *dpd)
{
    gretl_matrix *iZZ, *M, *M1, *M2;
    int k, m, err = 0;
    
    trim_zero_inst(dpd->XZ, dpd->ZZ, dpd->ZY);

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

#if DPDEBUG > 0
    gretl_matrix_print(M1, "M1");
    gretl_matrix_print(M2, "M2");
#endif

    gretl_matrix_free(iZZ);
    gretl_matrix_free(M);
    gretl_matrix_free(M1);
    gretl_matrix_free(M2);

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
    
    *D = gretl_matrix_alloc(dpd->T, north);
    *H = gretl_identity_matrix_new(north);
    *Yi = gretl_matrix_alloc(1, north);
    *Xi = gretl_matrix_alloc(dpd->k, north);
    *Zi = gretl_matrix_alloc(dpd->nz, north);

    if (*D == NULL || *H == NULL || *Yi == NULL ||
	*Xi == NULL || *Zi == NULL) {
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
    gretl_matrix *H = NULL;
    gretl_matrix *Yi = NULL;
    gretl_matrix *Xi = NULL;
    gretl_matrix *Zi = NULL;
    int *goodobs;
    int Tshort;
    int nz_diff, nz_lev = 0;
    int north, t, unit = 0;
    int err = 0;

    /* do H as per Ox/dpd? */
    int dpdstyle = DPDSTYLE; /* defined at top for now */

    /* full T - maxlag - 1 */
    Tshort = dpd->T - dpd->p - 1;

    if (dpd->flags & DPD_SYSTEM) {
	north = Tshort * 2 + 1;
    } else {
	north = Tshort;
    }

    /* instruments specific to equations in differences */
    nz_diff = (dpd->T-1)*(dpd->T-2)/2;

    /* and to equations in levels */
    if (dpd->flags & DPD_SYSTEM) {
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

    /* initialize observation counts */
    dpd->effN = dpd->nobs = 0;

    for (t=pdinfo->t1; t<pdinfo->t2; t+=dpd->T) {
	int uT = check_unit_obs(dpd, goodobs, Z, t);

#if DPDEBUG > 0
	pprintf(prn, "\n\nUnit %4d:", unit);
	pprintf(prn, "usable obs = %4d:\n", uT);
#endif

	if (uT > 0) {
	    /* this unit is usable */
	    dpd->effN += 1;
	    dpd->nobs += uT;
	    build_unit_D_matrix(dpd, goodobs, D);
	    compute_H(D, H, Tshort, dpdstyle);
#if DPDEBUG > 2
	    if (t == pdinfo->t1) {
		gretl_matrix_print_to_prn(H, "H", prn);
	    }
#endif
	    build_Y(dpd, goodobs, Z, t, Yi);
	    build_X(dpd, goodobs, Z, t, Xi);
	    build_Z(dpd, goodobs, Z, t, nz_diff, nz_lev, Zi);
#if DPDEBUG > 1
	    gretl_matrix_print_to_prn(Yi, "Yi", prn);
	    gretl_matrix_print_to_prn(Xi, "Xi", prn);
	    gretl_matrix_print_to_prn(Zi, "Zi", prn);
#endif
	    gretl_matrix_multiply_mod(Xi, GRETL_MOD_NONE,
				      Zi, GRETL_MOD_TRANSPOSE,
				      dpd->XZ, GRETL_MOD_CUMULATE);
	    gretl_matrix_qform(Zi, GRETL_MOD_NONE,
			       H, dpd->ZZ, GRETL_MOD_CUMULATE);
	    gretl_matrix_multiply_mod(Zi, GRETL_MOD_NONE,
				      Yi, GRETL_MOD_TRANSPOSE,
				      dpd->ZY, GRETL_MOD_CUMULATE);
	}
	unit++;
    }

    pprintf(prn, "Total obs = %d\n", dpd->nobs);

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
   Blundell-Bond.

   The dpdinfo struct that is used here is a hybrid. It contains some
   useful common information, filled out by dpdinfo_new; it also has a
   fair number of members that are used by arbond but not (currently)
   by dpanel.  And it contains a few members that are specific to
   dpanel. This will have to be tidied up eventually.

   When dpdinfo_new() and dpd_allocate_matrices() have returned
   successfully, we have the following information available in the
   dpd struct:

   dpd->T:       total observations per unit (= pdinfo->pd)
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
   dpd->nz:   total number of instruments
   dpd->XZ:   cross-moment matrix (NULL, so far)
   dpd->ZZ:   ditto
   dpd->ZY:   ditto

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
	err = dpd_add_residual_vector(dpd, Z, pdinfo);
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

    dpdinfo_free(dpd);

    return mod;
}
