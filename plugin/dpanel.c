/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#define DPDEBUG 0

/* Populate the residual vector, dpd->uhat. In the system case
   we stack the residuals in levels under the residuals in
   differences, per unit. Calculate SSR and \sigma^2 while
   we're at it.
*/

static void dpanel_residuals (dpdinfo *dpd)
{
    const double *b = dpd->beta->val;
    double SSRd = 0.0, SSRl = 0.0;
    double x, ut;
    int i, j, k, t;

    k = 0;

    for (i=0; i<dpd->N; i++) {
	unit_info *unit = &dpd->ui[i];
	int ndiff = unit->nobs - unit->nlev;

	for (t=0; t<ndiff; t++) {
	    /* differences */
	    ut = dpd->Y->val[k];
	    for (j=0; j<dpd->k; j++) {
		x = gretl_matrix_get(dpd->X, k, j);
		ut -= b[j] * x;
	    }
	    SSRd += ut * ut;
	    dpd->uhat->val[k++] = ut;
	}
	for (t=0; t<unit->nlev; t++) {
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
#if 0
	fprintf(stderr, "nobs (diffs) = %d\n", dpd->ndiff);
	fprintf(stderr, "SSR (diffs) = %g\n", SSRd);
#endif
	dpd->nobs = dpd->nlev;
	dpd->SSR = SSRl;
    } else {
	dpd->nobs = dpd->ndiff;
	dpd->SSR = SSRd;
    }

    dpd->s2 = dpd->SSR / (dpd->nobs - dpd->k);
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
		dpd->used[big_t] = LEVEL_ONLY;
	    }
	}
    }

    ok = goodobs[0];

    /* allow for differencing */
    if (ok > 0) ok--;

    return ok;
}

static int block_instrument_count (dpdinfo *dpd, int **Goodobs)
{
    /* the greatest lag we can actually support? */
    int maxlag = dpd->t2max - dpd->t1min;
    int nrows = 0;
    int i, j, t;

     for (i=0; i<dpd->nzb; i++) {
	if (dpd->d[i].minlag > maxlag) {
	    /* this spec is altogether otiose */
	    dpd->nzb -= 1;
	    for (j=i; j<dpd->nzb; j++) {
		dpd->d[j].v = dpd->d[j+1].v;
		dpd->d[j].v = dpd->d[j+1].minlag;
		dpd->d[j].v = dpd->d[j+1].maxlag;
		dpd->d[j].v = dpd->d[j+1].level;
	    }	    
	    i--;
	    continue;
	} 
	if (dpd->d[i].maxlag > maxlag) {
	    /* trim maxlag to what's feasible */
	    dpd->d[i].maxlag = maxlag;
	}
	for (t=dpd->t1min; t<=dpd->t2max; t++) {
	    /* FIXME levels equations */
	    for (j=dpd->d[i].minlag; j<=dpd->d[i].maxlag; j++) {
		if (t - j >= 0) {
		    nrows++;
		}
	    }
	}
    }

#if 1
    fprintf(stderr, "nzb = %d, extra Z rows = %d (total %d)\n", dpd->nzb,
	    nrows, dpd->nz + nrows);
#endif

    return nrows;
}

/* We do this accounting first so that we know the sizes of the
   various (and numerous) matrices that we'll need for the whole
   analysis at the outset; we can then allocate memory en bloc.
*/

static void do_unit_accounting (dpdinfo *dpd, const double **Z,
				int **Goodobs)
{
    int gmin, i, t;

    /* instruments specific to equations in differences */
    dpd->nzdiff = (dpd->T - 1) * (dpd->T - 2) / 2;

    /* and to equations in levels */
    if (use_levels(dpd)) {
	dpd->nzlev = dpd->T - 1;
    } 

    /* total instruments, so far */
    dpd->nz = dpd->nzdiff + dpd->nzlev + dpd->nzr;

    /* initialize observation counts */
    dpd->effN = dpd->ndiff = dpd->nlev = 0;
    dpd->minTi = dpd->T;

    /* initialize other accounts */
    gmin = (use_levels(dpd))? 1 : 2;
    dpd->t1min = dpd->T;
    dpd->t2max = 0;

    for (i=0, t=dpd->t1; i<dpd->N; i++, t+=dpd->T) {
	int *goodobs = Goodobs[i];
	int Ti = check_unit_obs(dpd, goodobs, Z, t);
	int gmax = goodobs[0];

#if DPDEBUG
	fprintf(stderr, "unit %d: Ti = %d\n", i+1, Ti);
#endif
	if (Ti > 0) {
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
	    if (goodobs[gmin] < dpd->t1min) {
		dpd->t1min = goodobs[gmin];
	    }
	    if (goodobs[gmax] > dpd->t2max) {
		dpd->t2max = goodobs[gmax];
	    }
	}
    }

    /* figure number of time dummies, if wanted */
    if (dpd->flags & DPD_TIMEDUM) {
	dpd->ndum = dpd->t2max - dpd->t1min;
	dpd->k += dpd->ndum;
	dpd->nz += dpd->ndum;
    }

    if (dpd->nzb > 0) {
	/* figure number of extra block-diagonal instruments */
	block_instrument_count(dpd, Goodobs);
    }

    /* figure the required number of columns for the Yi and Xi data
       matrices: this must be great enough to span the data range 
       for all units taken together
    */
    dpd->max_ni = dpd->t2max - dpd->t1min + 1;
    if (dpd->flags & DPD_SYSTEM) {
	dpd->max_ni += dpd->max_ni - 1;
    }      

    /* sum the total observations overall */
    dpd->totobs = dpd->ndiff + dpd->nlev;

#if DPDEBUG
    fprintf(stderr, "*** after dpanel accounting:\n"
	    " effN=%d, max_ni=%d, k=%d, ndum=%d, nz=%d\n", 
	    dpd->effN, dpd->max_ni, dpd->k, dpd->ndum, dpd->nz);
    fprintf(stderr, " maxTi=%d, minTi=%d\n", dpd->maxTi, dpd->minTi);
#endif    
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

#if DPDEBUG
    gretl_matrix_print(H, "dpdstyle H");
#endif
}

static void build_unit_H_matrix (dpdinfo *dpd, int *goodobs, 
				 gretl_matrix *D)
{
    if (D != NULL) {
	build_unit_D_matrix(dpd, goodobs, D);
	gretl_matrix_multiply_mod(D, GRETL_MOD_TRANSPOSE, 
				  D, GRETL_MOD_NONE, 
				  dpd->H, GRETL_MOD_NONE);
    } 
}

#define timedum_level(d,j,t) ((t == j + 1 + d->t1min)? 1 : 0)

static double timedum_diff (dpdinfo *dpd, int j, int t)
{
    int d0 = timedum_level(dpd, j, t);
    int d1 = timedum_level(dpd, j, t-1);

    return d0 - d1;
}

/* Build row vector of dependent variable values in @Yi using
   differences, followed by levels if wanted. 
*/

static int build_Y (dpdinfo *dpd, int *goodobs, const double **Z, 
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
	if (i1-1-maxlag >= Yi->cols) {
	    fprintf(stderr, "Bzzt! scribbling off the end of Yi\n"
		    " Yi->cols = %d; i1-1-maxlag = %d-1-%d = %d\n", 
		    Yi->cols, i1, maxlag, i1-1-maxlag);
	    return E_DATA;
	} else {
	    gretl_vector_set(Yi, i1-1-maxlag, dy);
	}
    }
    
    if (use_levels(dpd)) {
	/* levels */
	for (i=0; i<=usable; i++) {
	    i1 = goodobs[i+1];
	    t1 = t + i1;
	    gretl_vector_set(Yi, Tshort + (i1-maxlag), y[t1]);
	}
    }

    return 0;
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
    int row, col;
    double dx;

    gretl_matrix_zero(Xi);

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	t0 = t + i0;
	t1 = t + i1;

	row = 0;
	col = i1 - 1 - maxlag;

	for (j=1; j<=nlags; j++) {
	    lj = dpd->laglist[j];
	    dx = y[t1-lj] - y[t0-lj];
	    gretl_matrix_set(Xi, row++, col, dx);
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
	    gretl_matrix_set(Xi, row++, col, dx);
	}

	if (dpd->ndum > 0) {
	    for (j=0; j<dpd->ndum; j++) {
		if (use_levels(dpd)) {
		    dx = timedum_diff(dpd, j, i1);
		} else {
		    dx = timedum_level(dpd, j, i1);
		}
		gretl_matrix_set(Xi, row++, col, dx);
	    }
	}
    }
    
    if (use_levels(dpd)) {
	for (i=0; i<=usable; i++) {
	    i1 = goodobs[i+1];
	    t1 = t + i1;
	    row = 0;
	    col = Tshort + (i1-maxlag);

	    for (j=1; j<=nlags; j++) {
		lj = dpd->laglist[j];
		gretl_matrix_set(Xi, row++, col, y[t1-lj]);
	    }

	    for (j=1; j<=dpd->nx; j++) {
		xj = Z[dpd->xlist[j]];
		gretl_matrix_set(Xi, row++, col, xj[t1]);
	    }

	    for (j=0; j<dpd->ndum; j++) {
		dx = timedum_level(dpd, j, i1);
		gretl_matrix_set(Xi, row++, col, dx);
	    }
	}
    }
}

/* level instruments for the equations in differences */
 
static void gmm_inst_diff (const double *x, int t, int *goodobs, 
			   int maxlag, int qmax, int roffset,
			   gretl_matrix *Zi)
{
    int i, j, n = goodobs[0] - 1;
    int i0, i1, col, row;
    double x0;

    for (i=0; i<n; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	col = i1 - 1 - maxlag;
	row = roffset + i0 * (i0 - 1) / 2;
	for (j=0; j<i0; j++) {
	    if (qmax == 0 || j > i0 - qmax) {
		x0 = x[t+j];
		if (!na(x0)) {
		    gretl_matrix_set(Zi, row, col, x0);
		} 
	    } 
	    row++;
	}
    }
}

/* diff instruments for the equations in levels */
 
static void gmm_inst_lev (const double *x, int t, int *goodobs, 
			  int maxlag, int qmax, int roffset,
			  int coffset, gretl_matrix *Zi)
{
    int i, j, n = goodobs[0] - 1;
    int lastdiff = 0;
    int i0, col, row = roffset;
    double x0, x1;

    for (i=0; i<=n; i++) {
	i0 = goodobs[i+1];
	col = coffset + i0 - maxlag;
	for (j=lastdiff; j<i0; j++) {
	    if (qmax == 0 || j > i0 - qmax) {
		x0 = (j < 1)? NADBL : x[t+j-1];
		x1 = x[t+j];
		if (!na(x1) && !na(x0)) {
		    gretl_matrix_set(Zi, row, col, x1 - x0);
		}
	    }
	    row++;
	}
	lastdiff = j;
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
    double dx;
    int t0, t1, i0, i1;
    int k2 = dpd->nzdiff + dpd->nzlev;
    int k3 = k2 + dpd->nx;
    int i, j, col, row = 0;
    int qmax = dpd->qmax;

    gretl_matrix_zero(Zi);

    /* equations in differences: lagged levels of y */
    gmm_inst_diff(y, t, goodobs, maxlag, qmax, row, Zi);

    /* insert here: GMM-style instruments blocks for any
       exogenous vars selected for such treatment */
#if 0
    for (i=0; i<dpd->nzb; i++) {
	const double *x = Z[dpd->d[i].v];

	row = gmm_inst_diff(x, t, goodobs, maxlag, dpd->d[i].maxlag, row, Zi);
    }
#endif

    /* equations in differences: differenced exog vars 
       ("regular" instruments) 
    */
    if (dpd->nzr > 0) {
	for (i=0; i<usable; i++) {
	    i0 = goodobs[i+1];
	    i1 = goodobs[i+2];
	    col = i1 - 1 - maxlag;
	    t0 = t + i0;
	    t1 = t + i1;
	    for (j=0; j<dpd->nzr; j++) {
		/* Allin 2010-07-22: not differencing out the constant */
		if (!use_levels(dpd) && dpd->ilist[j+1] == 0) {
		    dx = 1.0;
		} else {
		    xj = Z[dpd->ilist[j+1]];
		    dx = xj[t1] - xj[t0];
		}
		gretl_matrix_set(Zi, k2 + j, col, dx);
	    }
	}
    }

    /* equations in differences: time dummies */
    if (dpd->ndum > 0 && !use_levels(dpd)) {
	for (i=0; i<usable; i++) {
	    i1 = goodobs[i+2];
	    col = i1 - 1 - maxlag;
	    for (j=0; j<dpd->ndum; j++) {
		dx = timedum_level(dpd, j, i1);
		gretl_matrix_set(Zi, k3 + j, col, dx);
	    }
	}	
    }

    if (use_levels(dpd)) {
	/* equations in levels: lagged differences of y */
	int roffset = (T-1) * (T-2) / 2;
	int coffset = T - maxlag - 1;

	gmm_inst_lev(y, t, goodobs, maxlag, qmax, roffset, coffset, Zi);

	/* insert here: GMMlevel-style blocks, if any */

	/* equations in levels: levels of exog vars:
	   FIXME we may need two "ilists", one for 
	   diffs and one for levels, and two "nzr"s?
	*/
	if (dpd->nzr > 0) {
	    for (i=0; i<=usable; i++) {
		i1 = goodobs[i+1];
		t1 = t + i1;
		col = (T-1-maxlag) + (i1-maxlag);
		for (j=0; j<dpd->nzr; j++) {
		    xj = Z[dpd->ilist[j+1]];
		    gretl_matrix_set(Zi, k2 + j, col, xj[t1]);
		}
	    }
	}

	/* equations in levels: time dummies */
	if (dpd->ndum > 0) {
	    for (i=0; i<=usable; i++) {
		i1 = goodobs[i+1];
		col = (T-1-maxlag) + (i1-maxlag);
		for (j=0; j<dpd->ndum; j++) {
		    dx = timedum_level(dpd, j, i1);
		    gretl_matrix_set(Zi, k3 + j, col, dx);
		}
	    }
	}
    }
}

static int trim_zero_inst (dpdinfo *dpd)
{
    char *mask;
    int err = 0;

#if WRITE_MATRICES
    gretl_matrix_write_as_text(dpd->A, "dpd-bigA.mat");
#endif

    mask = gretl_matrix_zero_diag_mask(dpd->A, &err);

    if (mask != NULL) {
	err = gretl_matrix_cut_rows_cols(dpd->A, mask);
	if (!err) {
	    dpd_shrink_matrices(dpd, mask);
	}
	free(mask);
    }

    if (!err) {
	gretl_matrix_divide_by_scalar(dpd->A, dpd->effN);
    }

#if DPDEBUG
    gretl_matrix_print(dpd->A, "dpd->A, after trim_zero_inst");
#endif

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
	gretl_matrix_free(*D);
	gretl_matrix_free(*Yi);
	gretl_matrix_free(*Xi);
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
    unit_info *unit = &dpd->ui[unum];
    double x;
    int i, j, k, s = *row;

    for (i=2; i<=goodobs[0]; i++) {
	k = goodobs[i] - 1 - dpd->p;
#if DPDEBUG
	if (i==2) {
	    fprintf(stderr, "Stacking Y for unit %d: starting at "
		    "k = %d - %d = %d (s=%d)\n", unum, goodobs[i], 
		    1 + dpd->p, k, s);
	}
#endif
	gretl_vector_set(dpd->Y, s, Yi->val[k]);
	for (j=0; j<Xi->rows; j++) {
	    x = gretl_matrix_get(Xi, j, k);
	    gretl_matrix_set(dpd->X, s, j, x);
	}
	for (j=0; j<dpd->nz; j++) {
	    x = gretl_matrix_get(Zi, j, k);
	    gretl_matrix_set(dpd->ZT, j, s, x);
	}
	s++;
    }

    /* record the indices of the first and last
       differenced observations */
    unit->t1 = goodobs[2];
    unit->t2 = goodobs[goodobs[0]];

    /* record the number of differenced obs */
    unit->nobs = (goodobs[0] > 0)? (goodobs[0] - 1) : 0;

    if (use_levels(dpd)) {
	int k0 = dpd->T - dpd->p - 1;

	for (i=1; i<=goodobs[0]; i++) {
	    k = k0 + goodobs[i] - 1;
	    gretl_vector_set(dpd->Y, s, Yi->val[k]);
	    for (j=0; j<Xi->rows; j++) {
		x = gretl_matrix_get(Xi, j, k);
		gretl_matrix_set(dpd->X, s, j, x);
	    }
	    for (j=0; j<dpd->nz; j++) {
		x = gretl_matrix_get(Zi, j, k);
		gretl_matrix_set(dpd->ZT, j, s, x);
	    }
	    s++;
	}

	/* record the number of levels obs and augment total */
	unit->nlev = goodobs[0];
	unit->nobs += unit->nlev;
    }

#if WRITE_MATRICES
    gretl_matrix_write_as_text(dpd->ZT, "dpdZT.mat");
#endif

    *row = s;
}

/* Main driver for system GMM: the core is a loop across
   the panel units to build the moment matrices. At this
   point we have already done the observations
   accounts, which are recorded in the Goodobs lists.
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
	return err;
    }

    Zi = dpd->Zi;
    gretl_matrix_reuse(Zi, dpd->nz, dpd->max_ni);
    
    if (D == NULL) {
	/* the H matrix will not vary by unit */
	int tau = dpd->t2max - dpd->t1min + 1;

	if (use_levels(dpd)) {
	    /* t1min is actually "levels-only" */
	    tau--;
	}
	make_dpdstyle_H(dpd->H, tau);
    }

    /* initialize cumulators */
    gretl_matrix_zero(dpd->XZ);
    gretl_matrix_zero(dpd->A);
    gretl_matrix_zero(dpd->ZY);

    /* initialize data stacker */
    Yrow = 0;

#if DPDEBUG
    gretl_matrix_zero(dpd->Y);
    gretl_matrix_zero(dpd->X);
#endif

    for (i=0; i<dpd->N; i++) {
	int *goodobs = Goodobs[i];
	int Ti = goodobs[0] - 1;

	if (Ti == 0) {
	    continue;
	}

	t = data_index(dpd, i);
	err = build_Y(dpd, goodobs, Z, t, Yi);
	if (err) {
	    break;
	}
	build_X(dpd, goodobs, Z, t, Xi);
	build_Z(dpd, goodobs, Z, t, Zi);
#if DPDEBUG
	gretl_matrix_print(Yi, "do_units: Yi");
	gretl_matrix_print(Xi, "do_units: Xi");
	gretl_matrix_print(Zi, "do_units: Zi");
#endif
	if (D != NULL) {
	    build_unit_H_matrix(dpd, goodobs, D);
	}
	gretl_matrix_qform(Zi, GRETL_MOD_NONE,
			   dpd->H, dpd->A, GRETL_MOD_CUMULATE);
	/* stack the individual data matrices for future use */
	stack_unit_data(dpd, Yi, Xi, Zi, goodobs, i, &Yrow);
    }

#if DPDEBUG
    gretl_matrix_print(dpd->Y, "dpd->Y");
    gretl_matrix_print(dpd->X, "dpd->X");
#endif

#if WRITE_MATRICES
    gretl_matrix_write_as_text(dpd->Y, "dpdY.mat");
    gretl_matrix_write_as_text(dpd->X, "dpdX.mat");
#endif    

    gretl_matrix_free(D);
    gretl_matrix_free(Yi);
    gretl_matrix_free(Xi);

    return err;
}

/* "Secret" public interface for new approach, including system GMM:
   in a script use --system (OPT_L, think "levels") to get
   Blundell-Bond. To build the H matrix as per Ox/DPD use the
   option --dpdstyle (OPT_X).
*/

MODEL dpd_estimate (const int *list, const char *ispec, 
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

    /* parse GMM instrument info, if present */
    if (ispec != NULL && *ispec != '\0') {
	mod.errcode = parse_GMM_instrument_spec(ispec, pdinfo, &d, &nzb);
	if (mod.errcode) {
	    return mod;
	}
    }

    dpd = dpdinfo_new(DPANEL, list, Z, pdinfo, opt, d, nzb, &mod.errcode);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in dpd_init\n", mod.errcode);
	return mod;
    }

#if 1
    if (ispec != NULL && *ispec != '\0') {
	int i;

	fprintf(stderr, "ispec = %s\n", ispec); 
	fprintf(stderr, "nzb = %d\n", dpd->nzb); 
	fprintf(stderr, "nzr = %d\n", dpd->nzr); 
	for (i=0; i<dpd->nzb; i++) {
	    fprintf(stderr, "var %d (%s): lags %d to %d (%s)\n", 
		    dpd->d[i].v, pdinfo->varname[dpd->d[i].v], 
		    dpd->d[i].minlag, dpd->d[i].maxlag, 
		    (dpd->d[i].level)? "Gmmlev" : "Gmm");
	}
	printlist(dpd->ilist, "ilist (regular Z)");
    }
#endif

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
	err = trim_zero_inst(dpd);
    }

    if (!err) {
	err = dpd_step_1(dpd);
    }

    if (!err && (opt & OPT_T)) {
	/* second step, if wanted */
	err = dpd_step_2(dpd);
    }	

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	/* write estimation info into model struct */
	mod.errcode = dpd_finalize_model(&mod, dpd, list, ispec, 
					 Z[dpd->yno], pdinfo, opt);
    }

    dpdinfo_free(dpd);

    return mod;
}
