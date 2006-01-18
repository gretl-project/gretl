/*
 *  Copyright (c) by Allin Cottrell and Riccardo Lucchetti
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* various functions for generating variables, mostly 
   used in generate.c */

#include "libgretl.h"
#include "libset.h"
#include "genstack.h"

#include <errno.h>

static double hp_lambda (const DATAINFO *pdinfo)
{
    double la = get_hp_lambda();

    if (na(la)) {
	la = 100 * pdinfo->pd * pdinfo->pd;
    }

    return la;
}

/**
 * hp_filter:
 * @x: array of original data.
 * @hp: array in which Hodrick-Prescott "cycle" is computed.
 * @pdinfo: data set information.
 *
 * Calculates the "cycle" component of the time series in
 * array @x, using the Hodrick-Prescott filter.  Adapted from the 
 * original FORTRAN code by E. Prescott. Very few changes.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int hp_filter (const double *x, double *hp, const DATAINFO *pdinfo)
{
    int i, t, T, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int err = 0;
    double v00 = 1.0, v11 = 1.0, v01 = 0.0;
    double det, tmp0, tmp1;
    double lambda;

    double **V = NULL;
    double m[2], tmp[2];

    int tb;
    double e0, e1, b00, b01, b11;

    for (t=t1; t<=t2; t++) {
	hp[t] = NADBL;
    }

    err = array_adjust_t1t2(x, &t1, &t2);
    if (err) {
	err = E_DATA;
	goto bailout;
    }

    T = t2 - t1 + 1;
    if (T < 4) {
	err = E_DATA;
	goto bailout;
    }

    lambda = hp_lambda(pdinfo);

    V = malloc(4 * sizeof *V);
    if (V == NULL) return E_ALLOC;

    for (i=0; i<4; i++) {
	V[i] = malloc(T * sizeof **V);
	if (V[i] == NULL) {
	    int j;
	    
	    for (j=0; j<i; j++) {
		free(V[j]);
	    }
	    free(V);
	    return E_ALLOC;
	}
    }

    /* adjust starting points */
    x += t1;
    hp += t1;

    /* covariance matrices for each obs */

    for (t=2; t<T; t++) {
	tmp0 = v00;
	tmp1 = v01;
	v00 = 1.0 / lambda + 4.0 * (tmp0 - tmp1) + v11;
	v01 = 2.0 * tmp0 - tmp1;
	v11 = tmp0;

	det = v00 * v11 - v01 * v01;

	V[0][t] =  v11 / det;
	V[1][t] = -v01 / det;
	V[2][t] =  v00 / det;

	tmp0 = v00 + 1.0;
	tmp1 = v00;
      
	v00 -= v00 * v00 / tmp0;
	v11 -= v01 * v01 / tmp0;
	v01 -= (tmp1 / tmp0) * v01;

    }

    m[0] = x[0];
    m[1] = x[1];

    /* forward pass */
    for (t=2; t<T; t++) {
	tmp[0] = m[1];
	m[1] = 2.0 * m[1] - m[0];
	m[0] = tmp[0];

	V[3][t-1] = V[0][t] * m[1] + V[1][t] * m[0];
	hp[t-1]   = V[1][t] * m[1] + V[2][t] * m[0];
	  
	det = V[0][t] * V[2][t] - V[1][t] * V[1][t];
	  
	v00 =  V[2][t] / det;
	v01 = -V[1][t] / det;
	  
	tmp[1] = (x[t] - m[1]) / (v00 + 1.0);
	m[1] += v00 * tmp[1];
	m[0] += v01 * tmp[1];
    }

    V[3][T-2] = m[0];
    V[3][T-1] = m[1];
    m[0] = x[T-2];
    m[1] = x[T-1];

    /* backward pass */
    for (t=T-3; t>=0; t--) {
	t1 = t+1;
	tb = T - t - 1;
      
	tmp[0] = m[0];
	m[0] = 2.0 * m[0] - m[1];
	m[1] = tmp[0];

	if (t > 1) {
	    /* combine info for y < i with info for y > i */
	    e0 = V[2][tb] * m[1] + V[1][tb] * m[0] + V[3][t];
	    e1 = V[1][tb] * m[1] + V[0][tb] * m[0] + hp[t];
	    b00 = V[2][tb] + V[0][t1];
	    b01 = V[1][tb] + V[1][t1];
	    b11 = V[0][tb] + V[2][t1];
	      
	    det = b00 * b11 - b01 * b01;
	      
	    V[3][t] = (b00 * e1 - b01 * e0) / det;
	}
	  
	det = V[0][tb] * V[2][tb] - V[1][tb] * V[1][tb];
	v00 =  V[2][tb] / det;
	v01 = -V[1][tb] / det;

	tmp[1] = (x[t] - m[0]) / (v00 + 1.0);
	m[1] += v01 * tmp[1];
	m[0] += v00 * tmp[1];
    }

    V[3][0] = m[0];
    V[3][1] = m[1];

    for (t=0; t<T; t++) {
	hp[t] = x[t] - V[3][t];
    }

 bailout:

    for (i=0; i<4; i++) {
	free(V[i]);
    }
    free(V);

    return err;
}

/**
 * bkbp_filter:
 * @y: array of original data.
 * @bk: array into which to write the filtered series.
 * @pdinfo: data set information.
 *
 * Calculates the Baxter & King bandpass filter.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int bkbp_filter (const double *y, double *bk, const DATAINFO *pdinfo)
{
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int periods[2];

    double omubar, omlbar;
    double avg_a;
    double *a;

    int i, k, t;
    int err = 0;

    /*
      periods[0] and periods[1]: threshold periodicities for business cycle
      k: order of the approximation
    */

    /* get user settings if available (or the defaults) */
    get_bkbp_periods(periods);
    k = get_bkbp_k();

#if BK_DEBUG
    fprintf(stderr, "lower limit = %d, upper limit = %d, \n", 
	    periods[0], periods[1]);
#endif

    if (periods[0] >= periods[1]) {
	strcpy(gretl_errmsg, "Error in Baxter-King frequencies");
	return 1;
    }

    err = array_adjust_t1t2(y, &t1, &t2);
    if (err) {
	return err;
    } 

    if (2 * k >= t2 - t1 + 1) {
	strcpy(gretl_errmsg, "Insufficient observations");
	return E_DATA;
    }

    a = malloc((k + 1) * sizeof *a);
    if (a == NULL) {
	return E_ALLOC;
    }
    
    omubar = 2.0 * M_PI / periods[0];
    omlbar = 2.0 * M_PI / periods[1];
    
    /* first we compute the coefficients */

    avg_a = a[0] = (omubar - omlbar) / M_PI;

    for (i=1; i<=k; i++) {
	a[i] = (sin(i * omubar) - sin(i * omlbar)) / (i * M_PI);
	avg_a += 2 * a[i];
    }

    avg_a /= (2 * k + 1);

    for (i=0; i<=k; i++) {
	a[i] -= avg_a;
#if BK_DEBUG
	fprintf(stderr, "a[%d] = %#9.6g\n", i, a[i]);
#endif
    }

    /* now we filter the series, skipping the first
       and last k observations */

    for (t=0; t<pdinfo->n; t++) {
	if (t < t1 + k || t >= t2 - k) {
	    bk[t] = NADBL;
	} else {
	    bk[t] = a[0] * y[t];
	    for (i=1; i<=k; i++) {
		bk[t] += a[i] * (y[t-i] + y[t+i]);
	    }
	}
    }

    free(a);

    return err;
}

/**
 * get_fracdiff:
 * @y: array of original data.
 * @diffvec: array into which to write the result.
 * @d: fraction by which to difference.
 * @pdinfo: data set information.
 *
 * Calculates the fractionally differenced counterpart
 * to the input series @y.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int get_fracdiff (const double *y, double *diffvec, double d,
		  const DATAINFO *pdinfo)
{
    int dd, t, T;
    const double TOL = 1.0E-07;
    int t1 = pdinfo->t1;
    int t2 = pdinfo->t2;
    double phi = -d;
    int err;

#if 0
    fprintf(stderr, "Doing get_fracdiff, with d = %g\n", d);
#endif

    err = array_adjust_t1t2(y, &t1, &t2);
    if (err) {
	return E_DATA;
    } 

    T = t2 - t1 + 1;

    for (t=0; t<pdinfo->n; t++) {
	if (t >= t1 && t <= t2) {
	    diffvec[t] = y[t];
	} else {
	    diffvec[t] = NADBL;
	}
    }   

    for (dd=1; dd<=T && fabs(phi)>TOL; dd++) {
	for (t=t1+dd; t<=t2; t++) {
	    diffvec[t] += phi * y[t - dd];
	}
	phi *= (dd - d)/(dd + 1);
    }

    return 0;
}

static int panel_x_offset (const DATAINFO *pdinfo, int *bad)
{
    char *p = strchr(pdinfo->stobs, ':');
    int offset = 0;

    if (p == NULL) {
	p = strchr(pdinfo->stobs, '.');
    }

    if (p == NULL) {
	*bad = 1;
    } else {
	offset = atoi(p + 1) - 1;
    }

    return offset;
}

static void 
make_x_panel_dummy (double *x, const DATAINFO *pdinfo, int i)
{
    int t, offset, bad = 0;
    int dmin, dmax;

    offset = panel_x_offset(pdinfo, &bad);

    dmin = (i - 1) * pdinfo->pd;
    dmax = i * pdinfo->pd - offset;

    if (i > 1) dmin -= offset;

    for (t=0; t<pdinfo->n; t++) {
	if (bad) {
	    x[t] = NADBL;
	} else if (t >= dmin && t < dmax) {
	    x[t] = 1.0;
	} else {
	    x[t] = 0.0;
	}
    }
}

static int n_new_dummies (const DATAINFO *pdinfo,
			  int nunits, int nperiods)
{
    char dname[VNAMELEN];
    int i, nnew = nunits + nperiods;

    for (i=0; i<nunits; i++) {
	sprintf(dname, "du_%d", i + 1);
	if (varindex(pdinfo, dname) < pdinfo->v) {
	    nnew--;
	}
    }

    for (i=0; i<nperiods; i++) {
	sprintf(dname, "dt_%d", i + 1);
	if (varindex(pdinfo, dname) < pdinfo->v) {
	    nnew--;
	}
    }

    return nnew;
}

static void 
make_dummy_name_and_label (int vi, const DATAINFO *pdinfo, int center,
			   char *vname, char *vlabel)
{
    if (center > 0) {
	sprintf(vname, "S%d", vi);
	strcpy(vlabel, "centered periodic dummy");
    } else if (center < 0) {
	sprintf(vname, "S%d", vi);
	strcpy(vlabel, "uncentered periodic dummy");
    } else if (pdinfo->pd == 4 && pdinfo->structure == TIME_SERIES) {
	sprintf(vname, "dq%d", vi);
	sprintf(vlabel, 
		_("= 1 if quarter = %d, 0 otherwise"), vi);
    } else if (pdinfo->pd == 12 && pdinfo->structure == TIME_SERIES) {
	char mname[8];

	get_month_name(mname, vi);
	sprintf(vname, "d%s", mname);
	sprintf(vlabel, _("= 1 if month is %s, 0 otherwise"), mname);
    } else {
	char dumstr[8] = "dummy_";
	char numstr[8];
	int len;

	sprintf(numstr, "%d", vi);
	len = strlen(numstr);
	dumstr[8 - len] = '\0';
	sprintf(vname, "%s%d", dumstr, vi);
	sprintf(vlabel, _("%s = 1 if period is %d, 0 otherwise"), vname, vi);
    }
}

/**
 * dummy:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @center: if greater than zero subtract the population mean from
 * each of the generated dummies; if less than zero, do not
 * subtract the mean but generate dummies with labels on the
 * same pattern as centered dummies (for internal use in VECMs).
 * Usually this argument is set to zero.
 *
 * Adds to the data set (if these variables are not already 
 * present) a set of periodic (usually seasonal) dummy variables.
 *
 * Returns: the ID number of the first dummy variable on success,
 * or 0 on error.
 */

int dummy (double ***pZ, DATAINFO *pdinfo, int center)
{
    char vname[USER_VLEN];
    char vlabel[MAXLABEL];
    int vi, t, yy, pp, mm;
    int ndums, nnew = 0;
    int di, di0 = pdinfo->v;
    double xx, dx;

    if (pdinfo->structure == STACKED_CROSS_SECTION) {
	ndums = pdinfo->n / pdinfo->pd;
	if (pdinfo->n % pdinfo->pd) {
	    ndums++;
	}
    } else {
	ndums = pdinfo->pd;
    }

    if (ndums == 1 || ndums > 99999) {
	strcpy(gretl_errmsg, _("This command won't work with the current periodicity"));
	return 0;
    }

    for (vi=0; vi<ndums; vi++) {
	make_dummy_name_and_label(vi + 1, pdinfo, center, vname, vlabel);
	di = varindex(pdinfo, vname);
	if (di >= pdinfo->v || strcmp(vlabel, VARLABEL(pdinfo, di))) {
	    nnew++;
	} else if (vi == 0) {
	    di0 = di;
	} else if (di != di0 + vi) {
	    /* dummies not consecutive: problem */
	    di0 = pdinfo->v;
	    nnew = ndums;
	    break;
	}
    }

    if (nnew == 0) {
	/* all dummies already present */
	return di0;
    } else if (pZ == NULL) {
	return -1;
    }

    if (dataset_add_series(ndums, pZ, pdinfo)) {
	strcpy(gretl_errmsg, _("Out of memory error"));
	return 0;
    }

    pp = pdinfo->pd;
    mm = 10;
    while ((pp = pp / 10)) {
	mm *= 10;
    }

    for (vi=1, di = di0; vi<=ndums; vi++, di++) {
	make_dummy_name_and_label(vi, pdinfo, center, vname, vlabel);
	strcpy(pdinfo->varname[di], vname);
	strcpy(VARLABEL(pdinfo, di), vlabel);

	if (pdinfo->structure == STACKED_CROSS_SECTION) {
	    make_x_panel_dummy((*pZ)[di], pdinfo, vi);
	} else {
	    for (t=0; t<pdinfo->n; t++) {
		xx = date(t, pdinfo->pd, pdinfo->sd0);
		if (dataset_is_daily(pdinfo)) { /* FIXME weekly? */
		    xx += .1;
		}
		yy = (int) xx;
		pp = (int) (mm * (xx - yy) + 0.5);
		dx = (pp == vi)? 1.0 : 0.0;
		(*pZ)[di][t] = dx;
	    }
	}
    }

    if (center > 0) {
	double cx = 1.0 / pdinfo->pd;
	int vimax = di0 + pdinfo->pd - 1;

	for (vi=di0; vi<=vimax; vi++) {
	    for (t=0; t<pdinfo->n; t++) {
		(*pZ)[vi][t] -= cx;
	    }
	}	
    }

    return di0;
}

/* if both == 1, generate both unit and period dummies,
   otherwise just generate unit dummies
*/

static int real_paneldum (double ***pZ, DATAINFO *pdinfo,
			  int both)
{
    char vname[16];
    int vi, t, yy, pp, mm;
    int xsect, orig_v = pdinfo->v;
    int ndum, nnew, n_blockdum = 0, n_freqdum = 0;
    int newvnum, offset, bad = 0;
    double xx;

    xsect = (pdinfo->structure == STACKED_CROSS_SECTION);

    /* in case xsect, block dummies are per time-period,
       frequency dummies are for units;
       in case of stacked time series, block dummies are
       per cross-sectional unit, frequency ones are for periods
    */

    if (both || xsect) {
	n_freqdum = pdinfo->pd;
	if (n_freqdum == 1) {
	    return E_PDWRONG;
	}
    }

    if (both || !xsect) {
	n_blockdum = pdinfo->n / pdinfo->pd;
	if (pdinfo->n % pdinfo->pd) {
	    n_blockdum++;
	}
	if (n_blockdum == 1) {
	    return E_PDWRONG;
	}
    }

    ndum = n_freqdum + n_blockdum;

    nnew = n_new_dummies(pdinfo, 
			 (xsect)? n_freqdum : n_blockdum,
			 (xsect)? n_blockdum : n_freqdum);

    if (dataset_add_series(nnew, pZ, pdinfo)) {
	return E_ALLOC;
    }

    pp = pdinfo->pd;
    mm = 10;
    while ((pp = pp / 10)) {
	mm *= 10;
    }

    newvnum = orig_v;

    /* first generate the frequency-based dummies */
    for (vi=1; vi<=n_freqdum; vi++) {
	int dnum;

	if (xsect) {
	    sprintf(vname, "du_%d", vi);
	} else {
	    sprintf(vname, "dt_%d", vi);
	}

	dnum = varindex(pdinfo, vname);
	if (dnum >= orig_v) {
	    dnum = newvnum++;
	}

	strcpy(pdinfo->varname[dnum], vname);
	sprintf(VARLABEL(pdinfo, dnum), 
		_("%s = 1 if %s is %d, 0 otherwise"), vname, 
		(xsect)? _("unit"): _("period"), vi);

	for (t=0; t<pdinfo->n; t++) {
	    xx = date(t, pdinfo->pd, pdinfo->sd0);
	    yy = (int) xx;
	    pp = (int) (mm * (xx - yy) + 0.5);
	    (*pZ)[dnum][t] = (pp == vi)? 1.0 : 0.0;
	}
    }

    offset = panel_x_offset(pdinfo, &bad);

    /* and then the block-based ones */
    for (vi=1; vi<=n_blockdum; vi++) {
	int dmin = (vi-1) * pdinfo->pd;
	int dmax = vi * pdinfo->pd - offset;
	int dnum;

	if (vi > 1) dmin -= offset;

	if (xsect) {
	    sprintf(vname, "dt_%d", vi);
	} else {
	    sprintf(vname, "du_%d", vi);
	}

	dnum = varindex(pdinfo, vname);
	if (dnum >= orig_v) {
	    dnum = newvnum++;
	}	

	strcpy(pdinfo->varname[dnum], vname);
	sprintf(VARLABEL(pdinfo, dnum), 
		_("%s = 1 if %s is %d, 0 otherwise"), vname, 
		(xsect)? _("period"): _("unit"), vi);

	for (t=0; t<pdinfo->n; t++) {
	    if (bad) {
		(*pZ)[dnum][t] = NADBL;
	    } else if (t >= dmin && t < dmax) {
		(*pZ)[dnum][t] = 1.0;
	    } else {
		(*pZ)[dnum][t] = 0.0;
	    }
	}
    }

    return 0;
}

/**
 * panel_unit_dummies:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Adds to the data set a set of dummy variables corresponding
 * to the cross-sectional units in a panel.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int panel_unit_dummies (double ***pZ, DATAINFO *pdinfo)
{
    return real_paneldum(pZ, pdinfo, 0);
}

/**
 * paneldum:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Adds to the data set a set of panel data dummy variables (for
 * both unit and period).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int paneldum (double ***pZ, DATAINFO *pdinfo)
{
    return real_paneldum(pZ, pdinfo, 1);
}

/**
 * genrunit:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * (For panel data only) adds to the data set an index variable 
 * that uniquely identifies the cross-sectional units.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int genrunit (double ***pZ, DATAINFO *pdinfo)
{
    int xt = 0;
    int i, t;

    if (pdinfo->structure != STACKED_TIME_SERIES &&
	pdinfo->structure != STACKED_CROSS_SECTION) {
	strcpy(gretl_errmsg, "'genr unit' can be used only with "
	       "panel data");
	return 1;
    }

    i = varindex(pdinfo, "unit");

    if (i == pdinfo->v) {
	if (dataset_add_series(1, pZ, pdinfo)) return E_ALLOC;
    }

    strcpy(pdinfo->varname[i], "unit");
    strcpy(VARLABEL(pdinfo, i), _("cross-sectional unit index"));

    if (pdinfo->structure == STACKED_CROSS_SECTION) {
	for (t=0; t<pdinfo->n; t++) {
	    if (t % pdinfo->pd == 0) {
		xt = 1;
	    }
	    (*pZ)[i][t] = (double) xt++;
	}
    } else {
	/* stacked time series */
	for (t=0; t<pdinfo->n; t++) {
	    if (t % pdinfo->pd == 0) {
		xt++;
	    }
	    (*pZ)[i][t] = (double) xt;
	}
    }    

    return 0;
}

/**
 * panel_unit_first_obs:
 * @t: zero-based observation number.
 * @pdinfo: data information struct.
 *
 * Returns: 1 if observation @t is the first time-series
 * observation on a given cross-sectional unit in a
 * panel dataset, else 0.
 */

int panel_unit_first_obs (int t, const DATAINFO *pdinfo)
{
    char *p, obs[OBSLEN];
    int ret = 0;

    ntodate(obs, t, pdinfo);
    p = strchr(obs, ':');
    if (p != NULL && atoi(p + 1) == 1) {
	ret = 1;
    }

    return ret;
}

/* make special time variable for panel data */

static void 
make_panel_time_var (double *x, const DATAINFO *pdinfo)
{
    int t, xt = 0;

    if (pdinfo->structure == STACKED_TIME_SERIES) {
	for (t=0; t<pdinfo->n; t++) {
	    if (t % pdinfo->pd == 0) {
		xt = 1;
	    }
	    x[t] = (double) xt++;
	}
    } else {
	/* stacked cross-sections */
	for (t=0; t<pdinfo->n; t++) {
	    if (t % pdinfo->pd == 0) {
		xt++;
	    }
	    x[t] = (double) xt;
	}
    }
}

/**
 * genrtime:
 * @pZ: pointer to data array.
 * @pdinfo: data information struct.
 * @tm: if non-zero, an actual time trend is wanted,
 * otherwise just an index of observations.
 *
 * Generates (and adds to the dataset, if it's not already
 * present) a time-trend or index variable.  The function
 * is panel-data aware: if the dataset is a panel and
 * @tm is non-zero, the trend will not simply run
 * consecutively over the entire range of the data, but
 * will correctly represent the location in time of
 * each observation.  The index is 1-based.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int genrtime (double ***pZ, DATAINFO *pdinfo, int tm)
{
    int i, t;

    i = varindex(pdinfo, (tm)? "time" : "index");

    if (i == pdinfo->v && dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    if (tm) {
	strcpy(pdinfo->varname[i], "time");
	strcpy(VARLABEL(pdinfo, i), _("time trend variable"));
    } else {
	strcpy(pdinfo->varname[i], "index");
	strcpy(VARLABEL(pdinfo, i), _("data index variable"));
    }
    
    if (tm && 
	(pdinfo->structure == STACKED_TIME_SERIES ||
	 pdinfo->structure == STACKED_CROSS_SECTION)) {
	make_panel_time_var((*pZ)[i], pdinfo);
    } else {
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[i][t] = (double) (t + 1);
	}
    }

    return 0;
}

static int plotvar_is_full_size (int v, int n, const double *x)
{
    int t, ret = 1;

    for (t=0; t<n; t++) {
	if (na(x[t])) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

/**
 * plotvar:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @period: string to identify periodicity: "annual", "qtrs",
 * "months", "decdate" (calendar data), "time" or "index".
 *
 * Adds to the data set a special dummy variable for use in plotting.
 *
 * Returns: the ID number of the variable (> 0) or -1 on failure
 */

int plotvar (double ***pZ, DATAINFO *pdinfo, const char *period)
{
    int t, vi, y1, n = pdinfo->n;
    float rm;

    vi = varindex(pdinfo, period);

    if (vi < pdinfo->v) {
	if (plotvar_is_full_size(vi, pdinfo->n, (*pZ)[vi])) {
	    return vi;
	} 
    } else if (dataset_add_series(1, pZ, pdinfo)) {
	return -1;
    }

    strcpy(pdinfo->varname[vi], period);

    y1 = (int) pdinfo->sd0;
    rm = pdinfo->sd0 - y1;

    switch (period[0]) {
    case 'a':
	strcpy(VARLABEL(pdinfo, vi), _("annual plotting variable")); 
	for (t=0; t<n; t++) {
	    (*pZ)[vi][t] = (double) (t + atoi(pdinfo->stobs));
	}
	break;
    case 'q':
	strcpy(VARLABEL(pdinfo, vi), _("quarterly plotting variable"));
	(*pZ)[vi][0] = y1 + (10.0 * rm - 1.0) / 4.0;
	for (t=1; t<n; t++) {
	    (*pZ)[vi][t] = (*pZ)[vi][t-1] + .25;
	}
	break;
    case 'm':
	strcpy(VARLABEL(pdinfo, vi), _("monthly plotting variable"));
	(*pZ)[vi][0] = y1 + (100.0 * rm - 1.0) / 12.0;
	for (t=1; t<n; t++) {
	    (*pZ)[vi][t] = (*pZ)[vi][t-1] + (1.0 / 12.0);
	}
	break;
    case 'h':
	strcpy(VARLABEL(pdinfo, vi), _("hourly plotting variable"));
	(*pZ)[vi][0] = y1 + (100.0 * rm - 1.0) / 24.0;
	for (t=1; t<n; t++) {
	    (*pZ)[vi][t] = (*pZ)[vi][t-1] + (1.0 / 24.0);
	}
	break;
    case 'd':
	if ((dated_daily_data(pdinfo) && pdinfo->n > 365) ||
	    (dated_weekly_data(pdinfo) && pdinfo->n > 52)) {
	    strcpy(VARLABEL(pdinfo, vi), _("daily plotting variable"));
	    for (t=0; t<n; t++) {
		if (pdinfo->S != NULL) {
		    (*pZ)[vi][t] = get_dec_date(pdinfo->S[t]);
		} else {
		    char datestr[OBSLEN];
		    
		    calendar_date_string(datestr, t, pdinfo);
		    (*pZ)[vi][t] = get_dec_date(datestr);
		}
	    } 
	} else if (dataset_is_decennial(pdinfo)) {
	    strcpy(pdinfo->varname[vi], "time");
	    strcpy(VARLABEL(pdinfo, vi), _("time trend variable"));
	    for (t=0; t<n; t++) {
		(*pZ)[vi][t] = pdinfo->sd0 + 10 * t;
	    }	    
	} else {
	    strcpy(pdinfo->varname[vi], "time");
	    strcpy(VARLABEL(pdinfo, vi), _("time trend variable"));
	    for (t=0; t<n; t++) {
		(*pZ)[vi][t] = (double) (t + 1);
	    }
	}
	break; 
    case 'i':
	strcpy(VARLABEL(pdinfo, vi), _("index variable"));
	for (t=0; t<n; t++) {
	    (*pZ)[vi][t] = (double) (t + 1);
	}
	break;
    case 't':
	strcpy(VARLABEL(pdinfo, vi), _("time trend variable"));
	for (t=0; t<n; t++) {
	    (*pZ)[vi][t] = (double) (t + 1);
	}
	break;
    default:
	break;
    }

    return vi;
}

static int varnames_from_arg (const char *s, char *v1str, char *v2str)
{
    int i, p, n = strlen(s);

    if (n > 17) {
	return 1;
    }

    p = haschar(',', s);
    if (p < 0 || p > 8) {
	return 1;
    }

    /* get first var name */
    for (i=0; i<p; i++) {
	v1str[i] = s[i];
    }
    v1str[p] = '\0';

    /* get second var name */
    n = n - p - 1;
    for (i=0; i<n; i++) {
	v2str[i] = s[p+1+i];
    }
    v2str[i] = '\0';

    return 0;
}

double genr_cov_corr (const char *s, double ***pZ, 
		      const DATAINFO *pdinfo, int fn)
{
    char v1str[USER_VLEN], v2str[USER_VLEN];
    int v1, v2;
    double ret = NADBL;

    if (varnames_from_arg(s, v1str, v2str)) {
	return NADBL;
    }

    v1 = varindex(pdinfo, v1str);
    v2 = varindex(pdinfo, v2str);

    if (v1 >= pdinfo->v || v2 >= pdinfo->v) {
	return NADBL;
    }

    if (fn == T_COV) {
	ret = gretl_covar(pdinfo->t1, pdinfo->t2, (*pZ)[v1], (*pZ)[v2]);
    } else if (fn == T_CORR) {
	ret = gretl_corr(pdinfo->t1, pdinfo->t2, (*pZ)[v1], (*pZ)[v2], NULL);
    }

    return ret;
}

static int 
get_model_param_number (const MODEL *pmod, const char *vname)
{
    int i, ret = 0;

    if (pmod->params == NULL) {
	return 0;
    }

    for (i=0; i<=pmod->ncoeff; i++) {
	if (!strcmp(vname, pmod->params[i])) {
	    ret = i + 1;
	    break;
	}
    }

    return ret;
}

double genr_vcv (const char *s, const DATAINFO *pdinfo, MODEL *pmod)
{
    int v1 = 0, v2 = 0;
    int i, j, k, v1l, v2l;
    char v1str[USER_VLEN], v2str[USER_VLEN];
    int gotit;
    double ret = NADBL;

    if (pmod == NULL || pmod->list == NULL) {
	return NADBL;
    }

    if (varnames_from_arg(s, v1str, v2str)) {
	return NADBL;
    }

    /* are the varnames valid? */
    if (pmod->ci != NLS && pmod->ci != ARMA) {
	v1 = varindex(pdinfo, v1str);
	v2 = varindex(pdinfo, v2str);
	if (v1 >= pdinfo->v || v2 >= pdinfo->v) {
	    return NADBL;
	}
    }

    /* check model list */
    if (pmod->ci == NLS || pmod->ci == ARMA) {
	v1l = get_model_param_number(pmod, v1str);
	v2l = get_model_param_number(pmod, v2str);
    } else {
	v1l = gretl_list_position(v1, pmod->list);
	v2l = gretl_list_position(v2, pmod->list);
    }

    if (v1l == 0 || v2l == 0) {
	return NADBL;
    }

    v1l -= 2;
    v2l -= 2;

    /* make model vcv matrix if need be */
    if (pmod->vcv == NULL && makevcv(pmod)) {
	return NADBL;
    }

    /* now find the right entry */
    if (v1l > v2l) {
	k = v1l;
	v1l = v2l;
	v2l = k;
    }

    gotit = 0;
    k = 0;
    for (i=0; i<pmod->ncoeff && !gotit; i++) {
	for (j=0; j<pmod->ncoeff; j++) {
	    if (j < i) {
		continue;
	    }
	    if (i == v1l && j == v2l) {
		ret = pmod->vcv[k];
		gotit = 1;
		break;
	    }
	    k++;
	}
    }

    return ret;
}

/**
 * genr_fit_resid:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @code: GENR_RESID or GENR_FITTED or GENR_RESID2.
 * @undo: if non-zero, don't bother labeling the variables
 *
 * Adds residuals or fitted values or squared residuals from a
 * given model to the data set.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int genr_fit_resid (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		    int code, int undo)
{
    char vname[USER_VLEN], vlabel[MAXLABEL];
    int i, t;
    double *h = NULL;

    if (code == GENR_H) {
	h = gretl_model_get_data(pmod, "garch_h");
	if (h == NULL) return E_DATA;
    }

    if (dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    i = pdinfo->v - 1;

    for (t=0; t<pdinfo->n; t++) {
	(*pZ)[i][t] = NADBL;
    }

    if (code == GENR_RESID) {
	sprintf(vname, "uhat%d", pmod->ID);
	sprintf(vlabel, _("residual from model %d"), pmod->ID);
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    (*pZ)[i][t] = pmod->uhat[t];
	}
    } else if (code == GENR_FITTED) {
	sprintf(vname, "yhat%d", pmod->ID);
	sprintf(vlabel, _("fitted value from model %d"), pmod->ID);
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    (*pZ)[i][t] = pmod->yhat[t];
	}
    } else if (code == GENR_RESID2) { 
	/* squared residuals */
	sprintf(vname, "usq%d", pmod->ID);
	sprintf(vlabel, _("squared residual from model %d"), pmod->ID);
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (na(pmod->uhat[t])) {
		(*pZ)[i][t] = NADBL;
	    } else {
		(*pZ)[i][t] = pmod->uhat[t] * pmod->uhat[t];
	    }
	}
    } else if (code == GENR_H) { 
	/* garch variance */
	sprintf(vname, "h%d", pmod->ID);
	sprintf(vlabel, _("fitted variance from model %d"), pmod->ID);
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    (*pZ)[i][t] = h[t];
	}
    }

    strcpy(pdinfo->varname[i], vname);

    if (!undo) {
	strcpy(VARLABEL(pdinfo, i), vlabel);
    }

    return 0;
}

int get_observation_number (const char *s, const DATAINFO *pdinfo)
{
    char test[OBSLEN];
    size_t n;
    int t;

    *test = 0;
    strncat(test, (*s == '"')? s + 1 : s, OBSLEN - 1);

    n = strlen(test);
    if (test[n-1] == '"') {
	test[n-1] = '\0';
    }

    if (pdinfo->markers && pdinfo->S != NULL) {
	for (t=0; t<pdinfo->n; t++) {
	    if (!strcmp(test, pdinfo->S[t])) return t + 1;
	}
	if (calendar_data(pdinfo)) {
	    charsub(test, ':', '/');
	    for (t=0; t<pdinfo->n; t++) {
		if (!strcmp(test, pdinfo->S[t]) ||
		    !strcmp(test, pdinfo->S[t] + 2)) {
		    return t + 1;
		}
	    }
	}
    }

    if (pdinfo->structure == TIME_SERIES) {
	t = dateton(test, pdinfo);
	if (t >= 0) return t + 1;
    }

    if (calendar_data(pdinfo)) {
	char datestr[OBSLEN];

	charsub(test, ':', '/');
	for (t=0; t<pdinfo->n; t++) {
	    calendar_date_string(datestr, t, pdinfo);
	    if (!strcmp(test, datestr) ||
		!strcmp(test, datestr + 2)) {
		return t + 1;
	    }
	}
    }

    return 0;
}

#define OBS_DEBUG 0

static int plain_obs_number (const char *obs, const DATAINFO *pdinfo)
{
    char *test;
    int t = -1;

#if OBS_DEBUG
    fprintf(stderr, "plain_obs_number: looking at '%s'\n", obs);
#endif

    errno = 0;

    strtol(obs, &test, 10);

    if (*test != '\0' || !strcmp(obs, test) || errno == ERANGE) {
#if OBS_DEBUG
	fprintf(stderr, "plain_obs_number: failed on '%s'\n", obs);
#endif
    } else {
	t = atoi(obs) - 1; /* convert from 1-based to 0-based */ 
    } 

    if (t >= pdinfo->n) {
	t = -1;
    }

    return t;
}

static void fix_calendar_date (char *s)
{
    while (*s) {
	if (*s == ':') *s = '/';
	s++;
    }
}

/* Given what looks like an observation number or date within "[" and
   "]", try to determine the observation number.  This is quite tricky
   since we try to handle both dates and plain observation numbers
   (and in addition, variables representing the latter); and we may
   have to deal with the translation from 1-based indexing in user
   space to 0-based indexing for internal purposes.
*/

int get_t_from_obs_string (char *s, const double **Z, 
			   const DATAINFO *pdinfo)
{
    int t;

#if OBS_DEBUG
    fprintf(stderr, "\nget_t_from_obs_string: s ='%s'\n", s);
#endif

    if (calendar_data(pdinfo)) {
	fix_calendar_date(s);
    } 

    t = dateton(s, pdinfo);

#if OBS_DEBUG
    fprintf(stderr, " dateton gives t = %d\n", t);
#endif

    if (t < 0) {
	if (isdigit((unsigned char) *s)) {
	    t = plain_obs_number(s, pdinfo);
#if OBS_DEBUG
	    fprintf(stderr, " plain_obs_number gives t = %d\n", t);
#endif
	} else {
	    int v = varindex(pdinfo, s);

	    if (v < pdinfo->v) {
		t = (int) Z[v][0];
#if OBS_DEBUG
		fprintf(stderr, " based on var %d: t = %d\n", v, t);
#endif
		if (t > pdinfo->n) {
		    char try[16];

		    sprintf(try, "%d", t);
		    t = dateton(try, pdinfo);
#if OBS_DEBUG
		    fprintf(stderr, " revised via dateton: t = %d\n", t);
#endif
		} else {
		    /* convert to 0-based */
		    t--;
		}
	    }
	}
    }

    if (t < 0) {
	sprintf(gretl_errmsg, _("Observation number out of bounds"));
    }

#if OBS_DEBUG
    fprintf(stderr, " return value: t = %d\n", t);
#endif

    return t;
}

static int arma_model_stat_pos (const char *s, const MODEL *pmod)
{
    int p = -1;

    if (numeric_string(s)) {
	p = atoi(s) - 1;
	if (p >= pmod->ncoeff) p = -1;
	return p;
    } else if (pmod->params != NULL) {
	int i;

	for (i=1; i<=pmod->ncoeff; i++) {
	    if (!strcmp(s, pmod->params[i])) {
		p = i - 1;
		break;
	    }
	}
    }

    return p;
}

/* retrieve a specific element from one of the arrays of data
   on a model */

double 
get_model_data_element (MODEL *pmod, int idx, const char *key,
			const DATAINFO *pdinfo, int *err)
{
    char s[32] = {0};
    const char *p;
    int lv, vi = 0;
    double x = NADBL;

    /* extract arg string in parentheses */
    p = strchr(key, '(');
    if (p != NULL) {
	sscanf(p + 1, "%31[^) ]", s);
    } else {
	strncat(s, key, 31);
    }

    DPRINTF(("get_model_data_element: looking at key = '%s'\n", s));

    if (idx == M_RHO_S) {
	if (!(numeric_string(s))) {
	    *err = E_INVARG;
	} else if (dot_atof(s) == 1 && AR1_MODEL(pmod->ci)) {
	    x = gretl_model_get_double(pmod, "rho_in");
	} else if (pmod->ci != AR && dot_atof(s) == 1) {
	    x = pmod->rho;
	} else if (pmod->arinfo == NULL || 
		   pmod->arinfo->arlist == NULL || 
		   pmod->arinfo->rho == NULL) {
	    *err = E_INVARG;
	} else if (!(vi = gretl_list_position(atoi(s), pmod->arinfo->arlist))) {
	    *err = E_INVARG;
	} else {
	    x = pmod->arinfo->rho[vi-1];
	}
    } else if (idx == M_VCV_S) {
	x = genr_vcv(s, pdinfo, pmod);
	if (na(x)) {
	    *err = E_INVARG;
	}
    } else if (idx == M_COEFF_S || idx == M_SE_S) {
	if (pmod == NULL || pmod->list == NULL) {
	    *err = E_INVARG;
	} else if (pmod->ci == ARMA) {
	    vi = arma_model_stat_pos(s, pmod);
	    if (vi < 0) {
		*err = E_INVARG;
	    }
	} else {
	    lv = numeric_string(s)? atoi(s) : varindex(pdinfo, s);
	    vi = gretl_list_position(lv, pmod->list);

	    if (vi < 2) {
		*err = E_INVARG;
	    } else {
		vi -= 2;
	    }
	}

	if (!*err) {
	    if (idx == M_COEFF_S && pmod->coeff != NULL) { 
		x = pmod->coeff[vi];
	    } else if (pmod->sderr != NULL) {
		x = pmod->sderr[vi];
	    } else {
		*err = E_INVARG;
	    }
	}
    } 

    if (*err) {
	gretl_errno = *err;
    }

    DPRINTF(("get_model_data_element: err = %d\n", *err));

    return x;
}
