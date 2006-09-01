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

struct model_grab {
    int idx;
    const char *mword;
};

struct model_grab grabs[] = {
    { M_ESS,        "$ess" },
    { M_T,          "$t" },
    { M_RSQ,        "$rsq" },
    { M_SIGMA,      "$sigma" },
    { M_DF,         "$df" },
    { M_NCOEFF,     "$ncoeff" },
    { M_LNL,        "$lnl" },
    { M_AIC,        "$aic" },
    { M_BIC,        "$bic" },
    { M_HQC,        "$hqc" },
    { M_TRSQ,       "$trsq" },
    { M_SCALAR_MAX,  NULL },
    { M_COEFF_S,    "$coeff" },
    { M_SE_S,       "$stderr" },
    { M_VCV_S,      "$vcv" },
    { M_RHO_S,      "$rho" },
    { M_ELEM_MAX,    NULL },
    { M_UHAT,       "$uhat" },
    { M_YHAT,       "$yhat" },
    { M_AHAT,       "$ahat" },
    { M_H,          "$h" },
    { M_SERIES_MAX,  NULL },
    { M_COEFF,      "$coeff" },
    { M_SE,         "$ess" },
    { M_VCV,        "$vcv" },
    { M_RHO,        "$rho" },
    { M_JALPHA,     "$jalpha" }, 
    { M_JBETA,      "$jbeta" },
    { M_JVBETA,     "$jvbeta" },
    { M_MAX,         NULL }
};

const char *get_model_stat_word (int idx)
{
    int i;

    for (i=0; grabs[i].idx != M_MAX; i++) {
	if (idx == grabs[i].idx) {
	    return grabs[i].mword;
	}
    }

    return NULL;
}

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
 * @hp: array in which filtered series is computed.
 * @pdinfo: data set information.
 * @opt: if %OPT_T, return the trend rather than the cycle.
 *
 * Calculates the "cycle" component of the time series in
 * array @x, using the Hodrick-Prescott filter.  Adapted from the 
 * original FORTRAN code by E. Prescott. Very few changes.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int hp_filter (const double *x, double *hp, const DATAINFO *pdinfo,
	       gretlopt opt)
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

    V = doubles_array_new(4, T);
    if (V == NULL) {
	return E_ALLOC;
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

    if (opt & OPT_T) {
	for (t=0; t<T; t++) {
	    hp[t] = V[3][t];
	}
    } else {
	for (t=0; t<T; t++) {
	    hp[t] = x[t] - V[3][t];
	}
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
    int bkl, bku;

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
    get_bkbp_periods(pdinfo, &bkl, &bku);
    k = get_bkbp_k(pdinfo);

#if BK_DEBUG
    fprintf(stderr, "lower limit = %d, upper limit = %d, \n", 
	    bkl, bku);
#endif

    if (bkl >= bku) {
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
    
    omubar = 2.0 * M_PI / bkl;
    omlbar = 2.0 * M_PI / bku;
    
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
    char vname[VNAMELEN];
    char vlabel[MAXLABEL];
    int vi, t, yy, pp, mm;
    int ndums = pdinfo->pd, nnew = 0;
    int di, di0 = pdinfo->v;
    double xx, dx;

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

/**
 * panel_dummies:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: %OPT_T for time dummies, otherwise unit dummies.
 *
 * Adds to the data set a set of dummy variables corresponding
 * to either the cross-sectional units in a panel, or the
 * time periods.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int panel_dummies (double ***pZ, DATAINFO *pdinfo, gretlopt opt)
{
    char vname[16];
    int vi, t, yy, pp, mm;
    int orig_v = pdinfo->v;
    int ndum, nnew;
    int n_unitdum = 0;
    int n_timedum = 0;
    int newvnum, offset, bad = 0;
    double xx;

    if (opt & OPT_T) {
	ndum = n_timedum = pdinfo->pd;
    } else {	
	n_unitdum = pdinfo->n / pdinfo->pd;
	if (pdinfo->n % pdinfo->pd) {
	    n_unitdum++;
	}
	ndum = n_unitdum;
    }
    
    if (ndum == 1) {
	return E_PDWRONG;
    }

    nnew = n_new_dummies(pdinfo, n_unitdum, n_timedum);

    if (dataset_add_series(nnew, pZ, pdinfo)) {
	return E_ALLOC;
    }

    pp = pdinfo->pd;
    mm = 10;
    while ((pp = pp / 10)) {
	mm *= 10;
    }

    newvnum = orig_v;

    /* generate time-based dummies, if wanted */

    for (vi=1; vi<=n_timedum; vi++) {
	int dnum;

	sprintf(vname, "dt_%d", vi);

	dnum = varindex(pdinfo, vname);
	if (dnum >= orig_v) {
	    dnum = newvnum++;
	}

	strcpy(pdinfo->varname[dnum], vname);
	sprintf(VARLABEL(pdinfo, dnum), 
		_("%s = 1 if %s is %d, 0 otherwise"), vname, 
		_("period"), vi);

	for (t=0; t<pdinfo->n; t++) {
	    xx = date(t, pdinfo->pd, pdinfo->sd0);
	    yy = (int) xx;
	    pp = (int) (mm * (xx - yy) + 0.5);
	    (*pZ)[dnum][t] = (pp == vi)? 1.0 : 0.0;
	}
    }

    offset = panel_x_offset(pdinfo, &bad);

    /* generate unit-based dummies, if wanted */

    for (vi=1; vi<=n_unitdum; vi++) {
	int dmin = (vi-1) * pdinfo->pd;
	int dmax = vi * pdinfo->pd - offset;
	int dnum;

	if (vi > 1) dmin -= offset;

	sprintf(vname, "du_%d", vi);

	dnum = varindex(pdinfo, vname);
	if (dnum >= orig_v) {
	    dnum = newvnum++;
	}	

	strcpy(pdinfo->varname[dnum], vname);
	sprintf(VARLABEL(pdinfo, dnum), 
		_("%s = 1 if %s is %d, 0 otherwise"), vname, 
		_("unit"), vi);

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

    if (pdinfo->structure != STACKED_TIME_SERIES) {
	strcpy(gretl_errmsg, "'genr unit' can be used only with "
	       "panel data");
	return 1;
    }

    i = varindex(pdinfo, "unit");

    if (i == pdinfo->v && dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    strcpy(pdinfo->varname[i], "unit");
    strcpy(VARLABEL(pdinfo, i), _("cross-sectional unit index"));

    for (t=0; t<pdinfo->n; t++) {
	if (t % pdinfo->pd == 0) {
	    xt++;
	}
	(*pZ)[i][t] = (double) xt;
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

    for (t=0; t<pdinfo->n; t++) {
	if (t % pdinfo->pd == 0) {
	    xt = 1;
	}
	x[t] = (double) xt++;
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
 * present) a time-trend or index variable.  This function
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
    
    if (tm && pdinfo->structure == STACKED_TIME_SERIES) {
	make_panel_time_var((*pZ)[i], pdinfo);
    } else {
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[i][t] = (double) (t + 1);
	}
    }

    return 0;
}

/**
 * genrwkday:
 * @pZ: pointer to data array.
 * @pdinfo: data information struct.
 *
 * Generates (and adds to the dataset, if it's not already
 * present) an index representing the day of the week for
 * each observation (for dated daily data only).
 * The index has value 0 for Sunday, 1 for Monday, and
 * so on.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int genrwkday (double ***pZ, DATAINFO *pdinfo)
{
    char datestr[OBSLEN];
    int i, t;

    if (!dated_daily_data(pdinfo)) {
	return E_PDWRONG;
    }

    i = varindex(pdinfo, "weekday");

    if (i == pdinfo->v && dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    strcpy(pdinfo->varname[i], "weekday");
    strcpy(VARLABEL(pdinfo, i), _("day of week (1 = Monday)"));
    
    for (t=0; t<pdinfo->n; t++) {
	ntodate_full(datestr, t, pdinfo);
	(*pZ)[i][t] = get_day_of_week(datestr);
    }

    return 0;
}

typedef enum {
    PLOTVAR_INDEX,
    PLOTVAR_TIME,
    PLOTVAR_ANNUAL,
    PLOTVAR_QUARTERS,
    PLOTVAR_MONTHS,
    PLOTVAR_CALENDAR,
    PLOTVAR_DECADES,
    PLOTVAR_HOURLY,
    PLOTVAR_MAX
} plotvar_type; 

int plotvar_code (const DATAINFO *pdinfo)
{
    if (!dataset_is_time_series(pdinfo)) {
	return PLOTVAR_INDEX;
    } else if (pdinfo->pd == 1) {
	return PLOTVAR_ANNUAL;
    } else if (pdinfo->pd == 4) {
	return PLOTVAR_QUARTERS;
    } else if (pdinfo->pd == 12) {
	return PLOTVAR_MONTHS;
    } else if (pdinfo->pd == 24) {
	return PLOTVAR_HOURLY;
    } else if ((dated_daily_data(pdinfo) && pdinfo->n > 365) ||
	    (dated_weekly_data(pdinfo) && pdinfo->n > 52)) {
	return PLOTVAR_CALENDAR;
    } else if (dataset_is_decennial(pdinfo)) {
	return PLOTVAR_DECADES;
    } else {
	return PLOTVAR_TIME;
    }
}

/**
 * gretl_plotx:
 * @pdinfo: data information struct.
 *
 * Finds or creates a special dummy variable for use on the
 * x-axis in plotting; this will have the full length of the
 * data series as given in @pdinfo, and will be appropriately
 * configured for the data frequency.  Do not try to free this
 * variable.
 *
 * Returns: pointer to plot x-variable, or %NULL on failure.
 */

const double *gretl_plotx (const DATAINFO *pdinfo)
{
    static double *x;
    static int ptype;
    static int T;

    int t, y1;
    int new_ptype;
    float rm;

    if (pdinfo == NULL) {
	/* cleanup signal */
	free(x);
	ptype = 0;
	T = 0;
	return NULL;
    }

    new_ptype = plotvar_code(pdinfo);

    if (x != NULL && new_ptype == ptype && T == pdinfo->n) {
	/* a suitable array is already at hand */
	return x;
    }

    if (x != NULL) {
	free(x);
    }

    x = malloc(pdinfo->n * sizeof *x);
    if (x == NULL) {
	return NULL;
    }

    T = pdinfo->n;
    ptype = new_ptype;

    y1 = (int) pdinfo->sd0;
    rm = pdinfo->sd0 - y1;

    switch (ptype) {
    case PLOTVAR_ANNUAL: 
	for (t=0; t<T; t++) {
	    x[t] = (double) (t + atoi(pdinfo->stobs));
	}
	break;
    case PLOTVAR_QUARTERS:
	x[0] = y1 + (10.0 * rm - 1.0) / 4.0;
	for (t=1; t<T; t++) {
	    x[t] = x[t-1] + .25;
	}
	break;
    case PLOTVAR_MONTHS:
	x[0] = y1 + (100.0 * rm - 1.0) / 12.0;
	for (t=1; t<T; t++) {
	    x[t] = x[t-1] + (1.0 / 12.0);
	}
	break;
    case PLOTVAR_HOURLY:
	x[0] = y1 + (100.0 * rm - 1.0) / 24.0;
	for (t=1; t<T; t++) {
	    x[t] = x[t-1] + (1.0 / 24.0);
	}
	break;
    case PLOTVAR_CALENDAR:
	for (t=0; t<T; t++) {
	    if (pdinfo->S != NULL) {
		x[t] = get_dec_date(pdinfo->S[t]);
	    } else {
		char datestr[OBSLEN];
		    
		calendar_date_string(datestr, t, pdinfo);
		x[t] = get_dec_date(datestr);
	    }
	}
	break;
    case PLOTVAR_DECADES:
	for (t=0; t<T; t++) {
	    x[t] = pdinfo->sd0 + 10 * t;
	}
	break;
    case PLOTVAR_INDEX:
	for (t=0; t<T; t++) {
	    x[t] = (double) (t + 1);
	}
	break;
    case PLOTVAR_TIME:
	for (t=0; t<T; t++) {
	    x[t] = (double) (t + 1);
	}
	break;
    default:
	break;
    }

    return x;
}

static int varnames_from_arg (const char *s, char *v1str, char *v2str)
{
    int i, p, n = strlen(s);

    if (n > 2 * VNAMELEN + 1) {
	return 1;
    }

    p = haschar(',', s);
    if (p < 0 || p > VNAMELEN - 1) {
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
    char v1str[VNAMELEN], v2str[VNAMELEN];
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

double genr_vcv (const char *s, const DATAINFO *pdinfo, MODEL *pmod)
{
    int v1 = 0, v2 = 0;
    int i, j, k;
    char v1str[VNAMELEN], v2str[VNAMELEN];
    int gotit;
    double ret = NADBL;

    if (pmod == NULL || pmod->list == NULL) {
	return NADBL;
    }

    if (varnames_from_arg(s, v1str, v2str)) {
	return NADBL;
    }

    v1 = gretl_model_get_param_number(pmod, pdinfo, v1str);
    if (v1 < 0) {
	return NADBL;
    }

    v2 = gretl_model_get_param_number(pmod, pdinfo, v2str);
    if (v2 < 0) {
	return NADBL;
    }

    /* make model vcv matrix if need be */
    if (pmod->vcv == NULL && makevcv(pmod, pmod->sigma)) {
	return NADBL;
    }

    /* now find the right entry */
    if (v1 > v2) {
	k = v1;
	v1 = v2;
	v2 = k;
    }

    gotit = 0;
    k = 0;
    for (i=0; i<pmod->ncoeff && !gotit; i++) {
	for (j=0; j<pmod->ncoeff; j++) {
	    if (j < i) {
		continue;
	    }
	    if (i == v1 && j == v2) {
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
    char vname[VNAMELEN], vlabel[MAXLABEL];
    const double *x = NULL;
    int i, t;

    if (code == GENR_H) {
	x = gretl_model_get_data(pmod, "garch_h");
	if (x == NULL) return E_DATA;
    } else if (code == GENR_AHAT) {
	x = gretl_model_get_data(pmod, "ahat");
	if (x == NULL) return E_DATA;
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
	    (*pZ)[i][t] = x[t];
	}
    } else if (code == GENR_AHAT) { 
	/* fixed-effects constants */
	sprintf(vname, "ahat%d", pmod->ID);
	sprintf(vlabel, _("per-unit constants from model %d"), pmod->ID);
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    (*pZ)[i][t] = x[t];
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
	    if (!strcmp(test, pdinfo->S[t])) {
		return t + 1;
	    }
	}
	if (calendar_data(pdinfo)) {
	    charsub(test, ':', '/'); /* FIXME */
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
	if (t >= 0) {
	    return t + 1;
	}
    }

    if (calendar_data(pdinfo)) {
	char datestr[OBSLEN];

	charsub(test, ':', '/'); /* FIXME */
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

static void invalid_stat_error (int idx, int *err)
{
    const char *mword = get_model_stat_word(idx);

    sprintf(gretl_errmsg, _("Invalid argument for %s"), 
	    (mword != NULL)? mword : _("function"));
    *err = E_INVARG;    
}

/* retrieve a specific element from one of the arrays of data
   on a model */

double 
get_model_data_element (MODEL *pmod, int idx, const char *key,
			const DATAINFO *pdinfo, int *err)
{
    char s[32] = {0};
    const char *p;
    int vi = 0;
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
	    invalid_stat_error(idx, err);
	} else if (dot_atof(s) == 1 && AR1_MODEL(pmod->ci)) {
	    x = gretl_model_get_double(pmod, "rho_in");
	} else if (pmod->ci != AR && dot_atof(s) == 1) {
	    x = pmod->rho;
	} else if (pmod->arinfo == NULL || 
		   pmod->arinfo->arlist == NULL || 
		   pmod->arinfo->rho == NULL) {
	    invalid_stat_error(idx, err);
	} else if (!(vi = gretl_list_position(atoi(s), pmod->arinfo->arlist))) {
	    invalid_stat_error(idx, err);
	} else {
	    x = pmod->arinfo->rho[vi-1];
	}
    } else if (idx == M_VCV_S) {
	x = genr_vcv(s, pdinfo, pmod);
	if (na(x)) {
	    invalid_stat_error(idx, err);
	}
    } else if (idx == M_COEFF_S || idx == M_SE_S) {
	if (pmod == NULL || pmod->list == NULL) {
	    invalid_stat_error(idx, err);
	} else {
	    vi = gretl_model_get_param_number(pmod, pdinfo, s);
	    if (vi < 0) {
		invalid_stat_error(idx, err);
	    }
	} 

	if (!*err) {
	    if (idx == M_COEFF_S && pmod->coeff != NULL) { 
		x = pmod->coeff[vi];
	    } else if (pmod->sderr != NULL) {
		x = pmod->sderr[vi];
	    } else {
		invalid_stat_error(idx, err);
	    }
	}
    } 

    if (*err) {
	gretl_errno = *err;
    }

    DPRINTF(("get_model_data_element: err = %d\n", *err));

    return x;
}
