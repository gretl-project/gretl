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

#include "libgretl.h"
#include "usermat.h"
#include "uservar.h"
#include "matrix_extra.h"
#include "gretl_midas.h"

#define MIDAS_DEBUG 0

/* days per month or quarter: maybe make this user-
   configurable? */

int midas_days_per_period (int days_per_week, int pd)
{
    int ret;
    
    if (days_per_week == 5) {
	ret = 22;
    } else if (days_per_week == 6) {
	ret = 26;
    } else {
	ret = 30;
    }

    return (pd == 12)? ret : 3 * ret;
}

/* Infer the current month from the current @qtr along
   with the number of days per period, @ndays, and the
   current index within the month-days array, @day.
*/

static int quarter_to_month (int qtr, int ndays, int day)
{
    return qtr * 3 - 2 + (day - 1) / (ndays/3);
}

/* Construct an auxiliary dataset in which the data from
   a MIDAS list are represented at their "native"
   frequency. We use this for pretty-printing a high
   frequency series.
*/

DATASET *midas_aux_dataset (const int *list,
			    const DATASET *dset,
			    int *err)
{
    DATASET *mset = NULL;
    gretlopt opt = 0;
    int mpd, pd = dset->pd;
    int T, m = list[0];
    int yr, mon;
    int daily = 0;

    if (m < 3 || gretl_list_has_separator(list)) {
	*err = E_INVARG;
    } else if (!dataset_is_time_series(dset)) {
	*err = E_INVARG;
    } else if (pd != 1 && pd != 4 && pd != 12) {
	/* host dataset should be annual, quarterly or monthly */
	*err = E_PDWRONG;
    }

    if (*err) {
	return NULL;
    }

    if (pd == 1) {
	/* annual: midas series should be quarterly or monthly */
	if (m != 4 && m != 12) {
	    *err = E_INVARG;
	} else {
	    mpd = m;
	}
    } else if (pd == 4) {
	/* quarterly: midas series should be monthly or daily */
	if (m == 3) {
	    mpd = 12;
	} else if (m == midas_days_per_period(5, 4)) {
	    mpd = 5;
	} else if (m == midas_days_per_period(6, 4)) {
	    mpd = 6;
	} else if (m == midas_days_per_period(7, 4)) {
	    mpd = 7;
	} else {
	    *err = E_INVARG;
	}
    } else {
	/* monthly: midas series should be daily */
	if (m == midas_days_per_period(5, 12)) {
	    mpd = 5;
	} else if (m == midas_days_per_period(6, 12)) {
	    mpd = 6;
	} else if (m == midas_days_per_period(7, 12)) {
	    mpd = 7;
	} else {
	    *err = E_INVARG;
	}
    }

    if (*err) {
	return NULL;
    }    

    if (!gretl_is_midas_list(list, dset)) {
	gretl_warnmsg_set("The argument does not seem to be a MIDAS list");
    }

    T = sample_size(dset) * m;

    if (mpd >= 5 && mpd <= 7) {
	/* we'll add markers for daily dates */
	daily = 1;
	opt = OPT_M;
    }

    mset = create_auxiliary_dataset(1, T, opt);
    if (mset == NULL) {
	*err = E_ALLOC;
    }

    if (!*err) {
	char *p, obs[OBSLEN];
	int nonex, qtr = 0;
	int i, t, s, m3 = 0;

	mset->pd = mpd;
	mset->structure = TIME_SERIES;
	strcpy(mset->varname[0], dset->varname[list[1]]);
	p = strrchr(mset->varname[0], '_');
	if (p != NULL) *p = '\0';

	ntodate(obs, dset->t1, dset);

	if (mpd == 4) {
	    sprintf(mset->stobs, "%d:1", atoi(obs));
	} else if (mpd == 12) {
	    sprintf(mset->stobs, "%d:01", atoi(obs));
	}

	if (daily && pd == 4) {
	    m3 = m / 3;
	}

	/* loop across observations in low-frequency dataset */

	s = 0;
	for (t=dset->t1; t<=dset->t2; t++) {
	    if (daily) {
		ntodate(obs, t, dset);
		sscanf(obs, "%d:%d", &yr, &mon);
		if (pd == 4) {
		    qtr = mon;
		}
	    }
	    /* read data right-to-left! */
	    for (i=m; i>0; i--) {
		int vi = list[i];
		
		if (daily) {
		    if (pd == 4) {
			mon = quarter_to_month(qtr, m, m-i+1);
			nonex = daily_index_to_date(mset->S[s], yr, mon,
						    (m-i) % m3, mpd);
		    } else {
			nonex = daily_index_to_date(mset->S[s], yr, mon,
						    m-i, mpd);
		    }
		    if (nonex) {
			/* skip any non-existent daily dates */
			mset->t2 -= 1;
		    } else {
			mset->Z[0][s++] = dset->Z[vi][t];
		    }
		} else {
		    mset->Z[0][s++] = dset->Z[vi][t];
		}
	    }
	}

	if (daily) {
	    strcpy(mset->stobs, mset->S[0]);
	    strcpy(mset->endobs, mset->S[mset->t2]);
	    mset->markers = DAILY_DATE_STRINGS;
	}

	mset->sd0 = get_date_x(mset->pd, mset->stobs);
	if (!daily) {
	    ntodate(mset->endobs, mset->t2, mset);
	}
    }

    return mset;
}

struct midas_info_ {
    char lname[VNAMELEN];  /* name of MIDAS list */
    char mname[VNAMELEN];  /* name of initial theta vector */
    int minlag;            /* minimum lag */
    int maxlag;            /* maximum lag */
    int type;              /* type of parameterization */
    int k;                 /* number of parameters */
    int nterms;            /* number of lag terms */
    int *laglist;          /* list of lag series */
};

typedef struct midas_info_ midas_info;

/* Parse a particular entry in the incoming array of
   "mds(foo,m1,m2,p,theta)" specifications.  We check that
   foo is an existent MIDAS list, and that m1, m2, p and
   theta all have admissible values.
*/

static int parse_midas_info (const char *s, midas_info *m,
			     const DATASET *dset)
{
    char lname[VNAMELEN];
    char mname[VNAMELEN];
    char fmt[48];
    int n, m1, m2, p;
    int umidas = 0;
    int err = 0;

    m->lname[0] = '\0';
    m->mname[0] = '\0';
    m->minlag = 0;
    m->maxlag = 0;
    m->type = 0;
    m->k = 0;
    m->nterms = 0;

    if (!strncmp(s, "mds(", 4)) {
	s += 4;
    }

    sprintf(fmt, "%%%d[^, ] , %%d , %%d , %%d, %%%d[^) ])",
	    VNAMELEN-1, VNAMELEN-1);

    n = sscanf(s, fmt, lname, &m1, &m2, &p, mname);

    if (n == 4 && p == 0) {
	umidas = 1;
    } else if (n != 5) {
	err = E_PARSE;
    }

    if (!err) {
	gretl_matrix *theta = NULL;
	int *list = get_list_by_name(lname);
	int k = 0;

	if (!umidas) {
	    theta = get_matrix_by_name(mname);
	}

	if (!gretl_is_midas_list(list, dset)) {
	    gretl_errmsg_set("mds(): the first term must be a MIDAS list");
	    err = E_INVARG;
	} else if (m1 > m2) {
	    err = E_INVARG;
	} else if (p < 0 || p >= MIDAS_MAX) {
	    err = E_INVARG;
	} else if (umidas) {
	    k = m2 - m1 + 1;
	} else {
	    k = gretl_vector_get_length(theta);
	    if (k < 1 || (p == MIDAS_BETA0 && k != 2) ||
		(p == MIDAS_BETAN && k != 3)) {
		err = E_INVARG;
	    }
	}

	if (!err) {
	    strcpy(m->lname, lname);
	    if (!umidas) {
		strcpy(m->mname, mname);
	    }
	    m->minlag = m1;
	    m->maxlag = m2;
	    m->type = p;
	    m->k = k;
	}
    }

    return err;
}

/* parse one or more strings of the form 

     mds(list, minlag, maxlag, type, theta) 
*/

static int 
parse_midas_specs (const char *spec, const DATASET *dset,
		   midas_info **pm, int *pnspec,
		   int *pumidas)
{
    midas_info *m = NULL;
    const char *s;
    int nspec = 0;
    int umidas = 0;
    int err = 0;

    /* first check: count closing parentheses */
    s = spec;
    while (*s) {
	if (*s == ')') {
	    nspec++;
	}
	s++;
    }

    if (nspec == 0) {
	/* spec is junk */
	err = E_PARSE;
    } else {
	/* allocate info structs */
	m = malloc(nspec * sizeof *m);
	if (m == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* parse and record individual MIDAS specs */
	char test[128];
	const char *p;
	int len, i = 0;

	s = spec;
	for (i=0; i<nspec && !err; i++) {
	    while (*s == ' ') s++;
	    p = s;
	    while (*p && *p != ')') p++;
	    len = p - s + 1;
	    if (len > 127) {
		err = E_PARSE;
	    } else {
		*test = '\0';
		strncat(test, s, len);
		err = parse_midas_info(test, &m[i], dset);
		if (!err && m[i].type == 0) {
		    umidas++;
		}
		s = p + 1;
	    }		
	}
    }

    if (!err && umidas > 0 && umidas < nspec) {
	/* FIXME? */
	gretl_errmsg_set("U-MIDAS specifications cannot be mixed with others");
	err = E_INVARG;
    }

    if (err) {
	free(m);
	*pnspec = 0;
    } else {
	*pm = m;
	*pnspec = nspec;
	*pumidas = umidas;
    }

    return err;
}

/* extract the list of regular (low-frequency) regressors */

static int *make_midas_xlist (const int *list, int *err)
{
    int i, nx = list[0] - 1;
    int *xlist = NULL;

    if (nx > 0) {
	xlist = gretl_list_new(nx);
	if (xlist == NULL) {
	    *err = E_ALLOC;
	} else {
	    for (i=1; i<=nx; i++) {
		xlist[i] = list[i+1];
	    }
	}
    }

    return xlist;
}

/* Given the name of an incoming MIDAS list plus minlag and
   maxlag values (at this point stored in the midas_info
   structure) build the list of required lags of the
   MIDAS series.
*/

static int *make_midas_laglist (midas_info *m,
				DATASET *dset,
				int *err)
{
    const int *list = get_list_by_name(m->lname);
    int *lcpy = gretl_list_copy(list);

    if (lcpy == NULL) {
	*err = E_ALLOC;
    } else {
	gretl_matrix *lv;

	lv = gretl_matrix_seq(m->minlag, m->maxlag, 1, err);
	if (!*err) {
	    *err = list_laggenr(&lcpy, 0, lv, dset, lcpy[0], OPT_L);
	    gretl_matrix_free(lv);
	}
    }

    return lcpy;
}

/* Build a full list of all series involved: dependent
   variable, regular regressors, and all lags of MIDAS
   terms. We want this either for setting the usable
   sample range -- or in the case of U-MIDAS, as the
   list to pass to OLS.
*/

static int *make_midas_biglist (const int *list,
				midas_info *m,
				int nmidas)
{
    int i, j, nt = list[0];
    int *biglist = NULL;

    for (j=0; j<nmidas; j++) {
	nt += m[j].nterms;
    }

    biglist = gretl_list_new(nt);

    if (biglist != NULL) {
	int k = 1;
	
	for (i=1; i<=list[0]; i++) {
	    biglist[k++] = list[i];
	}
	for (j=0; j<nmidas; j++) {
	    for (i=1; i<=m[j].nterms; i++) {
		biglist[k++] = m[j].laglist[i];
	    }
	}
    }

    return biglist;
}

/* If we're doing MIDAS via nls, set the usable
   sample range first. That way we'll know how
   big certain matrices ought to be.
*/

static int midas_set_sample (const int *list,
			     DATASET *dset,
			     midas_info *m,
			     int nmidas)
{
    int *biglist;
    int j, err = 0;

    biglist = make_midas_biglist(list, m, nmidas);
    if (biglist == NULL) {
	err = E_ALLOC;
    }
				 
    if (!err) {
	int t1 = dset->t1, t2 = dset->t2;

	err = list_adjust_sample(biglist, &t1, &t2, dset, NULL);
	if (!err) {
	    dset->t1 = t1;
	    dset->t2 = t2;
	}
    }

    for (j=0; j<nmidas; j++) {
	/* we're done with these lists now */
	free(m[j].laglist);
	m[j].laglist = NULL;
    }

    free(biglist);

    return err;
}

/* estimate U-MIDAS via OLS */

static int umidas_ols (MODEL *pmod,
		       const int *list,
		       DATASET *dset,
		       midas_info *m,
		       int nmidas,
		       gretlopt opt)
{
    int *biglist;
    int j;

    biglist = make_midas_biglist(list, m, nmidas);
    
    if (biglist == NULL) {
	return E_ALLOC;
    } else {
	*pmod = lsq(biglist, dset, OLS, opt | OPT_Z);
	free(biglist);
    }

    for (j=0; j<nmidas; j++) {
	/* we're done with these lists now */
	free(m[j].laglist);
	m[j].laglist = NULL;
    }

    return pmod->errcode;
}

int midas_model_calculate_depvar (MODEL *pmod,
				  DATASET *dset)
{
    const char *param;
    midas_info *minfo = NULL;
    int *xlist = NULL;
    int nmidas = 0;
    int umidas = 0;
    int err = 0;

    /* FIXME this is supposed to help with implementing
       forecasting, but it's not figured out yet
    */

    param = gretl_model_get_data(pmod, "midas_spec");
    err = parse_midas_specs(param, dset, &minfo,
			    &nmidas, &umidas);

    if (!err) {
	xlist = make_midas_xlist(pmod->list, &err);
	printlist(xlist, "midas xlist");
    }

    if (!err) {
	/* re-build MIDAS lag-lists */
	int i, *mlist;
	
	for (i=0; i<nmidas && !err; i++) {
	    mlist = make_midas_laglist(&minfo[i], dset, &err);
	    if (!err) {
		minfo[i].laglist = mlist;
	    }
	}	
    }

    free(xlist);
    
    return err;
}

/* Identify parameterizations which take an additional
   leading coefficient: all but U-MIDAS and the plain
   Almon polynomial specification.
*/

#define takes_coeff(t) (t != 0 && t != MIDAS_ALMON)

static int finalize_midas_model (MODEL *pmod,
				 const int *list,
				 const char *param,
				 const DATASET *dset,
				 midas_info *minfo,
				 int nmidas,
				 int hfslopes,
				 int umidas)
{
    gretl_matrix *m = NULL;
    gretl_matrix *theta = NULL;
    int i, rows = 0;
    int err = 0;

    gretl_model_set_string_as_data(pmod, "midas_spec",
				   gretl_strdup(param));

    if (umidas) {
	/* @pmod is the result of OLS estimation */
	int vi;
	
	gretl_model_allocate_param_names(pmod, pmod->ncoeff);
	for (i=0; i<pmod->ncoeff; i++) {
	    vi = pmod->list[i+2];
	    gretl_model_set_param_name(pmod, i, dset->varname[vi]);
	}
	gretl_model_set_int(pmod, "umidas", 1);
    } else {
	/* @pmod is the result of NLS estimation */
	pmod->ci = MIDASREG;
	free(pmod->depvar);
	pmod->depvar = gretl_strdup(dset->varname[list[1]]);
	free(pmod->list);
	pmod->list = gretl_list_copy(list);
    }

    for (i=0; i<nmidas; i++) {
	if (minfo[i].nterms > rows) {
	    rows = minfo[i].nterms;
	}
    }

    m = gretl_matrix_alloc(rows, nmidas);
    theta = gretl_matrix_alloc(rows, 1);

    if (m == NULL || theta == NULL) {
	err = E_ALLOC;
    } else {
	gretl_matrix *w = NULL;
	double *b = pmod->coeff;
	double wij, hfb = 0;
	int k = list[0] - 1;
	int pos = k + hfslopes;
	int ptype, j;

	for (i=0; i<nmidas && !err; i++) {
	    ptype = minfo[i].type;
	    if (ptype == 0) {
		/* U-MIDAS */
		for (j=0; j<minfo[i].k; j++) {
		    gretl_matrix_set(m, j, i, b[pos++]);
		}
		for (j=minfo[i].k; j<rows; j++) {
		    gretl_matrix_set(m, j, i, M_NA);
		}
	    } else {
		hfb = takes_coeff(ptype) ? 1.0 : b[k++];
		theta->rows = minfo[i].k;
		for (j=0; j<minfo[i].k; j++) {
		    theta->val[j] = b[pos++];
		}		
		w = midas_weights(minfo[i].nterms, theta,
				  ptype, &err);
		if (!err) {
		    for (j=0; j<minfo[i].nterms; j++) {
			wij = hfb * w->val[j];
			gretl_matrix_set(m, j, i, wij);
		    }
		    for (j=minfo[i].nterms; j<rows; j++) {
			gretl_matrix_set(m, j, i, M_NA);
		    }
		}
		gretl_matrix_free(w);
	    }
	}

	if (!err) {
	    /* save "gross" coefficients onto the model */
	    err = gretl_model_set_matrix_as_data(pmod, "midas_coeffs", m);
	} else {
	    gretl_matrix_free(m);
	}
    }

    gretl_matrix_free(theta);

    return err;
}

/* Define "private" matrices to hold the regular X data
   (MX___), the vector of coefficients on these data
   (bx___), and the vector of coefficients on the
   MIDAS list terms (high-frequency slopes).
*/

static int add_midas_matrices (const int *xlist,
			       const DATASET *dset,
			       midas_info *minfo,
			       int nmidas,
			       int *pslopes)
{
    gretl_matrix *m;
    int hfslopes = 0;
    int i, err = 0;

    if (xlist != NULL) {
	m = gretl_matrix_data_subset(xlist, dset,
				     dset->t1, dset->t2,
				     M_MISSING_ERROR,
				     &err);
	if (!err) {
	    err = private_matrix_add(m, "MX___");
	}
	if (!err) {
	    m = gretl_zero_matrix_new(xlist[0], 1);
	    if (m != NULL) {
		err = private_matrix_add(m, "bx___");
	    } else {
		err = E_ALLOC;
	    }
	}
    }
    
    if (!err) {
	for (i=0; i<nmidas; i++) {
	    if (takes_coeff(minfo[i].type)) {
		hfslopes++;
	    }
	}
	if (hfslopes > 0) {
	    m = gretl_zero_matrix_new(hfslopes, 1);
	    if (m != NULL) {
		err = private_matrix_add(m, "bmlc___");
	    } else {
		err = E_ALLOC;
	    }
	}
    }

    *pslopes = hfslopes;

    return err;
}

static int put_midas_nls_line (char *line,
			       DATASET *dset,
			       PRN *prn)
{
    int err;
    
#if MIDAS_DEBUG
    pputs(prn, line);
    pputc(prn, '\n');
#endif
    
    err = nl_parse_line(NLS, line, dset, prn);
    *line = '\0';

    return err;
}

static void append_pname (char *s, const char *pname)
{
    char c = s[strlen(s) - 1];
    
    if (c != ' ' && c != '"') {
	strncat(s, " ", 1);
    }
    strcat(s, pname);
}

MODEL midas_model (const int *list,
		   const char *param,
		   DATASET *dset,
		   gretlopt opt,
		   PRN *prn)
{
    midas_info *minfo = NULL;
    char tmp[64];
    int *xlist = NULL;
    int i, nmidas = 0;
    int hfslopes = 0;
    int origv = dset->v;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int umidas = 0;
    MODEL mod;
    int err = 0;

    gretl_model_init(&mod, dset);

    if (param == NULL || *param == '\0') {
	err = E_DATA;
    } else {
	err = parse_midas_specs(param, dset, &minfo,
				&nmidas, &umidas);
    }

    if (!err) {
	/* build list of regular regressors */
	xlist = make_midas_xlist(list, &err);
	if (xlist != NULL) {
	    err = remember_list(xlist, "XL___", NULL);
	    user_var_privatize_by_name("XL___");
	}
    }

    if (!err) {
	/* build MIDAS lag-lists */
	int *mlist;
	
	for (i=0; i<nmidas && !err; i++) {
	    mlist = make_midas_laglist(&minfo[i], dset, &err);
	    if (!err) {
		sprintf(minfo[i].lname, "ML___%d", i+1);
		err = remember_list(mlist, minfo[i].lname, NULL);
		minfo[i].nterms = mlist[0];
		minfo[i].laglist = mlist;
	    }
	}
    }

    if (!err && umidas) {
	err = umidas_ols(&mod, list, dset, minfo, nmidas, opt);
#if MIDAS_DEBUG	
	pputs(prn, "*** U-MIDAS via OLS ***\n");
#endif	
	goto umidas_finish;
    }

    if (!err) {
	/* determine usable sample range */
	err = midas_set_sample(list, dset, minfo, nmidas);
    }

    if (!err) {
	/* add the required matrices */
	err = add_midas_matrices(xlist, dset, minfo, nmidas,
				 &hfslopes);
    }

#if MIDAS_DEBUG	
    pputs(prn, "*** MIDAS via NLS ***\n\n");
#endif    

    if (!err) {
	char line[MAXLEN];
	int j = 0;

	/* initial "nls" line */
	sprintf(line, "nls %s = ", dset->varname[list[1]]);
	if (xlist != NULL) {
	    strcat(line, "lincomb(XL___, bx___)");
	}
	for (i=0; i<nmidas; i++) {
	    if (takes_coeff(minfo[i].type)) {
		sprintf(tmp, " + bmlc___[%d]*mlc___%d", ++j, i+1);
	    } else {
		sprintf(tmp, " + mlc___%d", i+1);
	    }
	    strcat(line, tmp);
	}
	err = put_midas_nls_line(line, dset, prn);

	/* MIDAS series and gradient matrices */
	for (i=0; i<nmidas && !err; i++) {
	    if (minfo[i].type == 0) {
		sprintf(line, "series mlc___%d = lincomb(%s, %s)",
			i+1, minfo[i].lname, minfo[i].mname);
	    } else {
		sprintf(line, "series mlc___%d = mlincomb(%s, %s, %d)",
			i+1, minfo[i].lname, minfo[i].mname,
			minfo[i].type);
	    }
	    err = put_midas_nls_line(line, dset, prn);
	    if (!err) {
		sprintf(line, "matrix mgr___%d = mgradient(%d, %s, %d)",
			i+1, minfo[i].nterms, minfo[i].mname, minfo[i].type);
		err = put_midas_nls_line(line, dset, prn);
	    }
	}

	/* derivatives */
	if (!err && xlist != NULL) {
	    strcpy(line, "deriv bx___ = MX___");
	    err = put_midas_nls_line(line, dset, prn);
	}
	if (!err && hfslopes > 0) {
	    strcpy(line, "deriv bmlc___ = {");
	    for (i=0; i<nmidas; i++) {
		if (takes_coeff(minfo[i].type)) {
		    sprintf(tmp, "mlc___%d", i+1);
		    strcat(line, tmp);
		}
		for (j=i+1; j<nmidas; j++) {
		    if (takes_coeff(minfo[j].type)) {
			strcat(line, ", ");
			break;
		    }
		}
	    }
	    strcat(line, "}");
	    err = put_midas_nls_line(line, dset, prn);
	}
	j = 0;
	for (i=0; i<nmidas && !err; i++) {
	    if (takes_coeff(minfo[i].type)) {
		sprintf(line, "deriv %s = bmlc___[%d] * {%s} * mgr___%d",
			minfo[i].mname, ++j, minfo[i].lname, i+1);
	    } else {
		sprintf(line, "deriv %s = {%s} * mgr___%d",
			minfo[i].mname, minfo[i].lname, i+1);
	    }		
	    err = put_midas_nls_line(line, dset, prn);
	}

	/* parameter names */
	if (!err) {
	    strcpy(line, "param_names \"");
	    if (xlist != NULL) {
		for (i=1; i<=xlist[0]; i++) {
		    strcpy(tmp, dset->varname[xlist[i]]);
		    append_pname(line, tmp);
		}
	    }
	    for (i=0; i<nmidas && !err; i++) {
		if (takes_coeff(minfo[i].type)) {
		    if (hfslopes > 1) {
			sprintf(tmp, "HF_slope%d", i+1);
		    } else {
			strcpy(tmp, "HF_slope");
		    }
		    append_pname(line, tmp);
		}
	    }
	    for (i=0; i<nmidas && !err; i++) {
		for (j=0; j<minfo[i].k; j++) {
		    sprintf(tmp, "%s[%d]", minfo[i].mname, j+1);
		    append_pname(line, tmp);
		}
	    }
	    strcat(line, "\"");
	    err = put_midas_nls_line(line, dset, prn);
	}
    }

#if MIDAS_DEBUG
    pputc(prn, '\n');
#endif

    if (!err) {
	mod = nl_model(dset, opt & OPT_G, prn);
    }

 umidas_finish:

    for (i=0; i<nmidas; i++) {
	user_var_delete_by_name(minfo[i].lname, NULL);
	if (!umidas) {
	    sprintf(tmp, "mgr___%d", i+1);
	    user_var_delete_by_name(tmp, NULL);
	}
    }

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	finalize_midas_model(&mod, list, param, dset,
			     minfo, nmidas, hfslopes,
			     umidas);
    }    

    destroy_private_uvars();
    free(minfo);
    free(xlist);

    if (dset->v > origv) {
	dataset_drop_last_variables(dset, dset->v - origv);
    }

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return mod;
}
