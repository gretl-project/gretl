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

#include "libset.h"

#define DPDEBUG 0
#define IVDEBUG 0

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

    if (dpd_style(dpd)) {
	dpd->s2 = dpd->SSR / (dpd->nobs - dpd->k);
    } else {
	/* xtabond2 always uses differences for this? */
	dpd->s2 = SSRd / (dpd->ndiff - dpd->k);
    }
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

	/* and lags of dependent variable? */
	for (i=1; i<=dpd->laglist[0] && ok; i++) {
	    s = t - dpd->laglist[i];
	    if (s < 0) {
		ok = 0;
	    } else {
		ok &= !na(y[s+t0]);
	    }
	}

	if (ok && dpd->xlist != NULL) {
	    /* and regressors? */
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

    /* allow for differencing (but don't set ok < 0) */
    if (ok > 0) ok--;

    return ok;
}

static void copy_diag_info (diag_info *targ, diag_info *src)
{
    targ->v = src->v;
    targ->minlag = src->minlag;
    targ->maxlag = src->maxlag;
    targ->level = src->level;
    targ->rows = src->rows;
}

/* diff_iv_accounts:

   On input tmin should be the first available obs in levels;
   tmax should be the second-last available obs in levels
   (in each case, for any unit). These indices are based at 0
   for the first period in the unit's data, and they 
   represent the subtractive terms in the first and last
   feasible observations in differences, respectively.

   minlag and maxlag represent the minimum and maximum
   lags that have been specified for the given instrument.
   These lags are relative to the "base" of the differenced
   observation, that is, the x_t from which x_{t-k} is
   subtracted to foem a difference. With non-gappy data
   k = 1 (the differences are x_t - x_{t-1}) but with
   gappy data we may have k > 1 for some observations;
   nonethless, we "impute" a difference-base index of
   (the subtractive term's index + 1). This ensures
   that we don't use level instruments that are entangled 
   in the difference they are supposed to be instrumenting.
*/

int diff_iv_accounts (dpdinfo *dpd, int tmin, int tmax)
{
    int t, tbot, ttop;
    int k, i, nrows = 0;

    tbot = tmin + 1;
    ttop = tmax + 1;

#if IVDEBUG
    fprintf(stderr, "*** diff_iv_accounts: tbot = %d, ttop = %d\n", tbot, ttop);
#endif

    for (i=0; i<dpd->nzb; i++) {
	int minlag = dpd->d[i].minlag;
	int maxlag = dpd->d[i].maxlag;
	int usable_maxlag = 0;
	int tbase = tmax + 2;
	int ii, itot = 0;

	dpd->d[i].rows = 0;

#if IVDEBUG	
	fprintf(stderr, "GMM spec %d, incoming: minlag = %d, maxlag = %d\n", 
		i, minlag, maxlag);
#endif

	/* find tbase = the 'base' of the first differenced observation 
	   for which there can be any usable instruments */

	for (t=tbot; t<=ttop; t++) {
	    if (t - minlag >= 0) {
		tbase = t;
		break;
	    }
	}

	if (tbase > ttop) {
	    fprintf(stderr, " no usable instruments for this spec\n");
	    dpd->nzb -= 1;
	    for (k=i; k<dpd->nzb; k++) {
		copy_diag_info(&dpd->d[k], &dpd->d[k+1]);
	    }	    
	    i--;
	    continue;
	}

#if IVDEBUG
	fprintf(stderr, " tbase = %d\n", tbase);
#endif

	/* step forward, cumulating in-principle usable instruments */
	for (t=tbase; t<=ttop; t++) {
	    ii = 0;
	    for (k=minlag; k<=maxlag && t-k >= 0; k++) {
		ii++;
		if (k > usable_maxlag) {
		    usable_maxlag = k;
		}
	    }
#if IVDEBUG
	    fprintf(stderr, "  ii = max insts at t=%d = %d\n", t, ii);
#endif
	    itot += ii;
	}

#if IVDEBUG    
	fprintf(stderr, " total insts = %d\n", itot);
	fprintf(stderr, " usable maxlag = %d\n", usable_maxlag);
#endif

	dpd->d[i].tbase = tbase;
	dpd->d[i].rows = itot;
	dpd->d[i].maxlag = usable_maxlag;
	nrows += itot;
    }

    return nrows;
}

/* lev_iv_accounts:

   On input tbot should be the first available obs in levels
   and ttop should be the last available obs in levels
   (in each case, for any unit). These indices are based at 0
   for the first period in the unit's data.

   minlag and maxlag represent the minimum and maximum
   lags that have been specified for the given instrument.
*/

int lev_iv_accounts (dpdinfo *dpd, int tbot, int ttop)
{
    int i, t, k, nrows = 0;

#if IVDEBUG    
    fprintf(stderr, "*** lev_iv_accounts: tbot = %d, ttop = %d\n", tbot, ttop);
#endif

    for (i=0; i<dpd->nzb2; i++) {
	int minlag = dpd->d2[i].minlag;
	int maxlag = dpd->d2[i].maxlag;
	int usable_maxlag = 0;
	int tbase = ttop + 1;
	int ii, itot = 0;

	dpd->d2[i].rows = 0;

#if IVDEBUG  	
	fprintf(stderr, "spec %d: minlag = %d, maxlag = %d\n", 
		i, minlag, maxlag);
#endif

	/* find tbase = the first obs in levels for which there can 
	   be any usable instruments; since these instruments are
	   differences we need to go back one step beyond minlag
	*/

	for (t=tbot; t<=ttop; t++) {
	    if (t - minlag - 1 >= 0) {
		tbase = t;
		break;
	    }
	}

	if (tbase > ttop) {
	    fprintf(stderr, " no usable instruments for this spec\n");
	    dpd->nzb2 -= 1;
	    for (k=i; k<dpd->nzb2; k++) {
		copy_diag_info(&dpd->d2[k], &dpd->d2[k+1]);
	    }	    
	    i--;
	    continue;
	}

#if IVDEBUG  
	fprintf(stderr, " tbase = %d\n", tbase);
#endif

	/* step forward, cumulating in-principle usable instruments */
	for (t=tbase; t<=ttop; t++) {
	    ii = 0;
	    for (k=minlag; k<=maxlag && t-k-1 >= 0; k++) {
		ii++;
		if (k > usable_maxlag) {
		    usable_maxlag = k;
		}
	    }
#if IVDEBUG
	    fprintf(stderr, "  ii = max insts at t=%d = %d\n", t, ii);
#endif
	    itot += ii;
	}

#if IVDEBUG    
	fprintf(stderr, " total insts = %d\n", itot);
	fprintf(stderr, " usable maxlag = %d\n", usable_maxlag);
#endif

	dpd->d2[i].tbase = tbase;
	dpd->d2[i].rows = itot;
	dpd->d2[i].maxlag = usable_maxlag;
	nrows += itot;
    }

    return nrows;
}

/* Work through the array of block-diagonal instrument
   specifications. Discard any useless ones, trim the
   maxlag values to what is supported on the data,
   and for each spec, count and record the implied number 
   of instrument rows that will appear in the Z matrix.
*/

static int block_instrument_count (dpdinfo *dpd, int t1lev, int t2pen)
{
    int nrows;

    nrows = diff_iv_accounts(dpd, t1lev, t2pen);
    dpd->nzdiff = nrows;
    dpd->nz += dpd->nzdiff;

#if IVDEBUG
    fprintf(stderr, "block_instrument_count, diffs: got %d rows (total = %d)\n", 
	    dpd->nzdiff, dpd->nz);
#endif

    if (use_levels(dpd)) {
	nrows = lev_iv_accounts(dpd, dpd->t1min, dpd->t2max);
	dpd->nzlev = nrows;
	dpd->nz += dpd->nzlev;
    } else {
	nrows = 0;
    }

#if IVDEBUG
    fprintf(stderr, "block_instrument_count, levels: got %d rows (total = %d)\n", 
	    dpd->nzlev, dpd->nz);
#endif  

    return 0;
}

/* We do this accounting first so that we know the sizes of the
   various (and numerous) matrices that we'll need for the whole
   analysis at the outset; we can then allocate memory en bloc.
*/

static void do_unit_accounting (dpdinfo *dpd, const double **Z,
				int **Goodobs)
{
    /* t1lev = index of first good obs in levels, and
       t2pen = index of penultimate good obs in levels 
    */
    int t1lev = dpd->T, t2pen = 0;
    int gmin, i, t;

    /* just make sure these are zeroed */
    dpd->nzdiff = dpd->nzlev = 0;

    /* total instruments, so far */
    dpd->nz = dpd->nzr;

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

#if DPDEBUG > 1
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
	    if (goodobs[1] < t1lev) {
		t1lev = goodobs[1];
	    }
	    if (goodobs[gmax] > dpd->t2max) {
		dpd->t2max = goodobs[gmax];
	    }
	    if (goodobs[gmax-1] > t2pen) {
		t2pen = goodobs[gmax-1];
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
	block_instrument_count(dpd, t1lev, t2pen);
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
    fprintf(stderr, " t1min=%d, t2max=%d\n", dpd->t1min, dpd->t2max);
#endif    
}

/* Based on the accounting of good observations for a unit recorded
   in the @goodobs list, fill matrix D (which is used to
   construct H unless we're doing things "dpdstyle").
*/

static void build_unit_D_matrix (dpdinfo *dpd, int *goodobs, gretl_matrix *D)
{
    int usable = goodobs[0] - 1;
    int colminus, maxlag = dpd->p;
    int i, j, i0, i1;    

    gretl_matrix_zero(D);

#if 1 /* maybe this is right in general? */
    colminus = dpd->t1min;
#else /* what we had before */
    colminus = 1 + maxlag;
#endif

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	j = i1 - colminus;
	gretl_matrix_set(D, i0, j, -1);
	gretl_matrix_set(D, i1, j,  1);
    }

    /* levels */
    if (use_levels(dpd)) {
	for (i=1; i<=goodobs[0]; i++) {
	    i1 = goodobs[i];
	    j = i1 - colminus;
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

#if DPDEBUG > 1
    gretl_matrix_print(H, "dpdstyle H");
#endif
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
    int colminus, maxlag = dpd->p;
    int T = dpd->T;
    int Tshort = T - maxlag - 1;
    int t0, t1, i0, i1;
    double dy;

    gretl_matrix_zero(Yi);

#if 1 /* maybe this is right in general? */
    colminus = dpd->t1min;
#else /* what we had before */
    colminus = 1 + maxlag;
#endif

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	t0 = t + i0;
	t1 = t + i1;
	dy = y[t1] - y[t0];
	if (i1-colminus >= Yi->cols) {
	    fprintf(stderr, "Bzzt! scribbling off the end of Yi\n"
		    " Yi->cols = %d; i1-colminus = %d-%d = %d\n", 
		    Yi->cols, i1, colminus, i1 - colminus);
	    return E_DATA;
	} else {
	    gretl_vector_set(Yi, i1-colminus, dy);
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
    int colminus, maxlag = dpd->p;
    int Tshort = dpd->T - maxlag - 1;
    const double *xj;
    int t0, t1, i0, i1;
    int i, j, lj;
    int row, col;
    double dx;

    gretl_matrix_zero(Xi);

#if 1 /* maybe this is right in general? */
    colminus = dpd->t1min;
#else /* what we had before */
    colminus = 1 + maxlag;
#endif

    /* differences */
    for (i=0; i<usable; i++) {
	i0 = goodobs[i+1];
	i1 = goodobs[i+2];
	t0 = t + i0;
	t1 = t + i1;

	row = 0;
	col = i1 - colminus;

	for (j=1; j<=nlags; j++) {
	    lj = dpd->laglist[j];
	    dx = y[t1-lj] - y[t0-lj];
	    gretl_matrix_set(Xi, row++, col, dx);
	}

	for (j=1; j<=dpd->nx; j++) {
	    /* Note: we don't difference away the constant
	       here, but if dpdstyle is not in force the
	       constant will have been removed already.
	    */
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
		} else if (dpd_style(dpd)) {
		    /* as per DPD: leave dummies in levels */
		    dx = timedum_level(dpd, j, i1);
		} else {
		    /* as per xtabond2 */
		    dx = timedum_diff(dpd, j, i1);
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

/* note: this amounts to t1*(t1-1)/2 in the
   straightforward case */

static int row_increment (diag_info *d, int t1)
{
    int k1 = d->level ? 1 : 0;
    int t, k, r = 0;

    for (t=d->tbase; t<t1; t++) {
	for (k=d->minlag; k<=d->maxlag && t-k-k1 >= 0; k++) {
	    r++;
	}
    }

    return r;
}

#if IVDEBUG

/* verify that we're not writing instrument values to rows of
   the Z matrix that are incompatible with what was figured
   out by block_instrument_count (see above).
*/

static int bad_write_check (dpdinfo *dpd, int row, int lev)
{
    if (!lev && row >= dpd->nzdiff) {
	fprintf(stderr, "*** ERROR in gmm_inst_diff: writing to "
		"bad row %d (max is %d)\n", row,
		dpd->nzdiff - 1);
	return 1;
    } else if (lev && (row < dpd->nzdiff || row >= dpd->nzdiff + dpd->nzlev)) {
	fprintf(stderr, "*** ERROR in gmm_inst_lev: writing to "
		"bad row %d (min is %d, max is %d)\n", row,
		dpd->nzdiff, dpd->nzdiff + dpd->nzlev - 1);
	return 1;
    }

    return 0;
}

#endif

/* GMM-style instruments in levels for the eqns in differences */

static int gmm_inst_diff (dpdinfo *dpd, int bnum, const double *x, 
			  int s, int *goodobs, int row0, int col0, 
			  gretl_matrix *Zi)
{
    int maxlag = dpd->d[bnum].maxlag;
    int minlag = dpd->d[bnum].minlag;
    int tmax = goodobs[goodobs[0]];
    int i, t, t1, t2;
    int colminus, col, row;
    double xt;

#if 1 /* maybe this is right in general? */
    colminus = dpd->t1min;
#else /* what we had before */
    colminus = 1 + dpd->p;
#endif

    for (i=1; i<goodobs[0]; i++) {
	t1 = goodobs[i];
	t2 = goodobs[i+1];
	col = col0 + t2 - colminus;
	row = row0 + row_increment(&dpd->d[bnum], t1+1);
	for (t=0; t<tmax; t++) {
	    if (t2 - t >= minlag && t1 - t < maxlag) {
		/* the criterion here needs some thought? */
		xt = x[s+t];
		if (!na(xt)) {
#if IVDEBUG
		    bad_write_check(dpd, row, 0);
#endif
		    gretl_matrix_set(Zi, row, col, xt);
		} 
		row++;
	    } 
	}

    }

    return row0 + dpd->d[bnum].rows;
}

/* GMM-style instruments in differences for the eqns in levels */

static int gmm_inst_lev (dpdinfo *dpd, int bnum, const double *x, 
			 int s, int *goodobs, int row0, int col0, 
			 gretl_matrix *Zi)
{
    int maxlag = dpd->d2[bnum].maxlag;
    int minlag = dpd->d2[bnum].minlag;
    int tmax = goodobs[goodobs[0]];
    int i, k, t, t1;
    int col, row;
    double x0, x1;

    for (i=1; i<=goodobs[0]; i++) {
	t1 = goodobs[i];
	col = col0 + t1 - dpd->p;
	row = row0 + row_increment(&dpd->d2[bnum], t1);
	for (t=1; t<=tmax; t++) {
	    k = t1 - t;
	    if (k <= maxlag && k >= minlag) {
		x0 = x[s+t-1];
		x1 = x[s+t];
		if (!na(x1) && !na(x0)) {
#if IVDEBUG
		    bad_write_check(dpd, row, 1);
#endif
		    gretl_matrix_set(Zi, row, col, x1 - x0);
		} 
		row++;
	    } 
	}
    }

    return row0 + dpd->d2[bnum].rows;
}

#if 0

/* GMM-style instruments in differences for the eqns in levels.
   In think that handling of row numbers is not right in the 
   general case where p may be > 1 and instrument lags limited.

   Keeping this older version around for comparison, for the
   moment. AC 2010-08-09
*/
 
static int gmm_inst_lev (dpdinfo *dpd, int bnum, const double *x,  
			 int s, int *goodobs, int row0, int col0, 
			 gretl_matrix *Zi)
{
    int maxlag = dpd->d2[bnum].maxlag;
    int minlag = dpd->d2[bnum].minlag;
    int i, k, t, t1;
    /* note: lastdiff seems to have a dependency on p,
       but exactly how that works is unclear */
    int lastdiff = 1; /* was 0 */
    int col, row = row0;
    double x0, x1;

    for (i=1; i<=goodobs[0]; i++) {
	t1 = goodobs[i];
	col = col0 + t1 - dpd->p;
	for (t=lastdiff; t<t1; t++) {
	    k = t1 - t;
	    if (k <= maxlag && k >= minlag) {
		x0 = x[s+t-1];
		x1 = x[s+t];
		if (!na(x1) && !na(x0)) {
		    gretl_matrix_set(Zi, row, col, x1 - x0);
		}
	    }
	    row++; /* not always right? */
	}
	lastdiff = t;
    }

    return row0 + dpd->d2[bnum].rows;
}

#endif /* 0 */

/* Build the matrix of per-unit instrument values in @Zi, which
   has the instruments in rows and the observations in columns.

   Note that each unit's Zi is the same size, padded with zero columns
   for missing observations as needed. The number of columns in Zi
   equals the maximal span of the data for all units taken together,
   counting both observations in differences and observations in
   levels, if applicable.

   We pack the instruments in the following order: 

   1) G1: GMM-style instruments in levels for equations in
      differences

   3) G2: GMM-style instruments in differences for equations in
      levels, if present

   5) I1: "Regular" instruments, differenced exog vars for eqns 
      in differences

   6) I2: "Regular" instruments, levels of exog vars for eqns 
      in levels, if any

   7) D1: Time dummies for eqns in differences, if specified and if
      "system" estimation is not being done

   8) D2: Time dummies for eqns in levels, if specified

   The pattern for the non-system case is

        Z' = | G1 : I1 : D1 |

   and for the full system case it is

        Z' = | G1 :  0 : I1 :  0 |
             |  0 : G2 : I2 : D2 |

*/

static void build_Z (dpdinfo *dpd, int *goodobs, const double **Z, 
		     int t, gretl_matrix *Zi)
{
    const int usable = goodobs[0] - 1;
    const int maxlag = dpd->p;
    const int T = dpd->T;
    const double *x;
    double dx;
    /* k2 is the starting row for "regular" instruments */
    int k2 = dpd->nzdiff + dpd->nzlev;
    /* k3 marks the starting row for time dummies */
    int k3 = k2 + dpd->nzr;
    int t0, t1, i0, i1;
    int i, j, col, row = 0;

    gretl_matrix_zero(Zi);

    /* GMM-style instruments in levels for diffs equations */
    for (i=0; i<dpd->nzb; i++) {
	x = Z[dpd->d[i].v];
	row = gmm_inst_diff(dpd, i, x, t, goodobs, row, 0, Zi);	
    }

    col = dpd->t2max - dpd->t1min;

    /* GMM-style instruments in diffs for levels equations */
    for (i=0; i<dpd->nzb2; i++) {
	x = Z[dpd->d2[i].v];
	row = gmm_inst_lev(dpd, i, x, t, goodobs, row, col, Zi);	
    }

    /* equations in differences: differenced exog vars */
    if (dpd->nzr > 0) {
	for (i=0; i<usable; i++) {
	    i0 = goodobs[i+1];
	    i1 = goodobs[i+2];
	    col = i1 - 1 - maxlag;
	    t0 = t + i0;
	    t1 = t + i1;
	    for (j=0; j<dpd->nzr; j++) {
		/* we don't difference the constant, but unless
		   dpdstyle is in force it will have been
		   dropped by this point
		*/
		if (!use_levels(dpd) && dpd->ilist[j+1] == 0) {
		    dx = 1.0;
		} else {
		    x = Z[dpd->ilist[j+1]];
		    dx = x[t1] - x[t0];
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
	/* equations in levels: levels of exog vars */
	if (dpd->nzr > 0) {
	    for (i=0; i<=usable; i++) {
		i1 = goodobs[i+1];
		t1 = t + i1;
		col = (T-1-maxlag) + (i1-maxlag);
		for (j=0; j<dpd->nzr; j++) {
		    x = Z[dpd->ilist[j+1]];
		    gretl_matrix_set(Zi, k2 + j, col, x[t1]);
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

#if DPDEBUG
    fprintf(stderr, "before trimming, order of A = %d\n", dpd->A->rows);
#endif

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

#if DPDEBUG
    gretl_matrix_print(dpd->A, "dpd->A, after trim_zero_inst");
#endif

    if (!err) {
	gretl_matrix_divide_by_scalar(dpd->A, dpd->effN);
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
    int colminus;
    int i, j, k, s = *row;

#if 1
    colminus = dpd->t1min;
#else
    colminus = 1 + dpd->p;
#endif

    for (i=2; i<=goodobs[0]; i++) {
	k = goodobs[i] - colminus;
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
	/* the starting position for levels in the data arrays */
	int k0 = dpd->t2max - dpd->t1min;

	for (i=1; i<=goodobs[0]; i++) {
	    k = k0 + goodobs[i] - dpd->p; /* note: was - 1 */
	    if (k >= Yi->cols) {
		fprintf(stderr, "*** stack_unit_data: reading off "
			"end of Yi (k=%d, Yi->cols=%d)\n", k, Yi->cols);
		fprintf(stderr, " at goodobs[%d] = %d\n", i, goodobs[i]);
		continue;
	    }
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
   the panel units to build the data and instrument matrices
   and cumulate A = \sum_i Z_i H_i Z_i'.

   At this point we have already done the observations
   accounts, which are recorded in the Goodobs lists.
*/

static int do_units (dpdinfo *dpd, const double **Z, 
		     int **Goodobs)
{
#if DPDEBUG
    char ystr[16];
#endif
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
    /* this should not be necessary if stack_unit_data() is
       working correctly */
    gretl_matrix_zero(dpd->Y);
    gretl_matrix_zero(dpd->X);
    gretl_matrix_zero(dpd->ZT);
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
	sprintf(ystr, "do_units: Y_%d", i);
	gretl_matrix_print(Yi, ystr);
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

/* If we're not doing DPD-style, the constant, if present, 
   gets differenced away, both as a regressor and as a regular 
   instrument. This is for the case where we are not including
   any equations in levels.
*/

static void maybe_prune_const (dpdinfo *dpd)
{
    int i;

    if (dpd->xlist != NULL) {
	for (i=dpd->xlist[0]; i>0; i--) {
	    if (dpd->xlist[i] == 0) {
		if (dpd->xlist[0] == 1) {
		    free(dpd->xlist);
		    dpd->xlist = NULL;
		} else {
		    gretl_list_delete_at_pos(dpd->xlist, i);
		}
		dpd->nx -= 1;
		dpd->k -= 1;
		break;
	    }
	}
    }

    if (dpd->ilist != NULL) {
	for (i=dpd->ilist[0]; i>0; i--) {
	    if (dpd->ilist[i] == 0) {
		if (dpd->ilist[0] == 1) {
		    free(dpd->ilist);
		    dpd->ilist = NULL;
		} else {
		    gretl_list_delete_at_pos(dpd->ilist, i);
		}
		dpd->nzr -= 1;
		break;
	    }
	}
    }
}

/* the user hasn't supplied a block-diagonal spec for
   y in the differences equations: here we set up the
   default version, with unlimited lags 
*/

static int insert_default_ydiff_spec (dpdinfo *dpd)
{
    diag_info *d;

    d = realloc(dpd->d, (dpd->nzb + 1) * sizeof *d);

    if (d == NULL) {
	return E_ALLOC;
    } else {
	/* insert the y spec in first place, moving
	   any other specs up */
	int i;

	dpd->d = d;

	for (i=dpd->nzb; i>0; i--) {
	    copy_diag_info(&dpd->d[i], &dpd->d[i-1]);
	}

	d = &dpd->d[0];
	d->v = dpd->yno;
	d->minlag = 2;
	d->maxlag = 99;
	d->level = 0;
	d->rows = 0;

	dpd->nzb += 1;
    }

    return 0;
}

/* the user has specified "system" but hasn't supplied a 
   block-diagonal spec for y in the levels equations: here 
   we set up the default version, with 1 lag 
*/

static int insert_default_ylev_spec (dpdinfo *dpd)
{
    diag_info *d;

    d = realloc(dpd->d, (dpd->nzb + 1) * sizeof *d);

    if (d == NULL) {
	return E_ALLOC;
    } else {
	/* append the y level spec to the diff specs */
	dpd->d = d;
	d = &dpd->d[dpd->nzb];
	d->v = dpd->yno;
	d->minlag = 1;
	d->maxlag = 1;
	d->level = 1;
	d->rows = 0;

	dpd->nzb += 1;
	dpd->nzb2 = 1;
    }

    return 0;
}

static int compare_gmm_specs (const void *a, const void *b)
{
    const diag_info *da = a;
    const diag_info *db = b;

    return (da->level - db->level);
}

/* Given the info on instrument specification returned by the
   parser that's in common between arbond and dpanel, make
   any adjustments that may be needed in the system case.
*/

static int dpanel_adjust_GMM_spec (dpdinfo *dpd)
{
    int have_yspec = 0;
    int i, err = 0;

    /* has the user given a GMM-style spec for the 
       dependent variable? */

    for (i=0; i<dpd->nzb; i++) {
	if (dpd->d[i].v == dpd->yno && dpd->d[i].level == 0) {
	    have_yspec = 1;
	    break;
	}
    }

    if (!have_yspec) {
	/* hmm, should we be doing this? */
	err = insert_default_ydiff_spec(dpd);
	if (err) {
	    return err;
	}
    }

    /* do we have any block-diagonal specs for levels eqns? */
    for (i=0; i<dpd->nzb; i++) {
	if (dpd->d[i].level) {
	    dpd->nzb2 += 1;
	}
    }

    if (dpd->nzb2 > 0 && dpd->nzb2 < dpd->nzb) {
	/* ensure the levels-equations specs come last */
	qsort(dpd->d, dpd->nzb, sizeof *dpd->d, compare_gmm_specs);
    }

    if (dpd->nzb2 == 0 && use_levels(dpd)) {
	/* add default GMM-level spec for y */
	err = insert_default_ylev_spec(dpd);
	if (err) {
	    return err;
	}	
    }    

    if (dpd->nzb2 > 0) {
	/* henceforth dpd->nzb refers to the number of specs in
	   differences, not the total number */
	dpd->nzb -= dpd->nzb2;
	dpd->d2 = dpd->d + dpd->nzb;
	dpd->flags |= DPD_SYSTEM; /* in case it's not present */
    }

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
    diag_info *d = NULL;
    dpdinfo *dpd = NULL;
    int **Goodobs = NULL;
    int nzb = 0;
    MODEL mod;
    int err = 0;

    gretl_model_init(&mod);
    gretl_model_smpl_init(&mod, pdinfo);

    if (libset_get_bool(DPDSTYLE)) {
	opt |= OPT_X;
    }

    /* parse GMM instrument info, if present */
    if (ispec != NULL && *ispec != '\0') {
	mod.errcode = parse_GMM_instrument_spec(DPANEL, ispec, pdinfo, &d, &nzb);
	if (mod.errcode) {
	    return mod;
	}
    }

    dpd = dpdinfo_new(DPANEL, list, Z, pdinfo, opt, d, nzb, &mod.errcode);
    if (mod.errcode) {
	fprintf(stderr, "Error %d in dpd_init\n", mod.errcode);
	return mod;
    }

    dpanel_adjust_GMM_spec(dpd);

    if (!dpd_style(dpd) && !use_levels(dpd)) { 
	maybe_prune_const(dpd);
    }

#if DPDEBUG
    if (dpd->nzb > 0 || dpd->nzb2 > 0) {
	int i;

	if (ispec != NULL) {
	    fprintf(stderr, "user's ispec = '%s'\n", ispec);
	} 
	fprintf(stderr, "nzb = %d, nzb2 = %d\n", dpd->nzb, dpd->nzb2); 
	fprintf(stderr, "nzr = %d\n", dpd->nzr); 
	for (i=0; i<dpd->nzb; i++) {
	    fprintf(stderr, "var %d (%s): lags %d to %d (GMM)\n", 
		    dpd->d[i].v, pdinfo->varname[dpd->d[i].v], 
		    dpd->d[i].minlag, dpd->d[i].maxlag); 
	}
	for (i=0; i<dpd->nzb2; i++) {
	    fprintf(stderr, "var %d (%s): lags %d to %d (GMMlevel)\n", 
		    dpd->d2[i].v, pdinfo->varname[dpd->d[i].v], 
		    dpd->d2[i].minlag, dpd->d2[i].maxlag); 
	}
	printlist(dpd->ilist, "ilist (regular instruments)");
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
