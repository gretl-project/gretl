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
	    /* read data right-to-left */
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
   theta have admissible values.
*/

static int parse_midas_info (const char *s, midas_info *m,
			     const DATASET *dset)
{
    char lname[VNAMELEN];
    char mname[VNAMELEN];
    char fmt[48];
    int m1, m2, p;
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

    if (sscanf(s, fmt, lname, &m1, &m2, &p, mname) != 5) {
	err = E_PARSE;
    } else {
	gretl_matrix *theta = get_matrix_by_name(mname);
	int *list = get_list_by_name(lname);
	int k = 0;

	if (!gretl_is_midas_list(list, dset)) {
	    gretl_errmsg_set("mds(): the first term must be a MIDAS list");
	    err = E_INVARG;
	} else if (m1 > m2) {
	    err = E_INVARG;
	} else if (p < 1 || p >= MIDAS_MAX) {
	    err = E_INVARG;
	} else {
	    k = gretl_vector_get_length(theta);
	    if (k < 1 || (p == MIDAS_BETA0 && k != 2) ||
		(p == MIDAS_BETAN && k != 3)) {
		err = E_INVARG;
	    }
	}

	if (!err) {
	    strcpy(m->lname, lname);
	    strcpy(m->mname, mname);
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
		   midas_info **pm, int *pnspec)
{
    midas_info *m = NULL;
    const char *s;
    int nspec = 0;
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
		s = p + 1;
	    }		
	}
    }

    if (err) {
	free(m);
	*pnspec = 0;
    } else {
	*pm = m;
	*pnspec = nspec;
    }

    return err;
}

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

static int midas_set_sample (const int *list,
			     DATASET *dset,
			     midas_info *m,
			     int nmidas)
{
    int i, j, nt = list[0];
    int *biglist = NULL;
    int err = 0;

    for (j=0; j<nmidas; j++) {
	nt += m[j].nterms;
    }

    biglist = gretl_list_new(nt);

    if (biglist == NULL) {
	err = E_ALLOC;
    } else {
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

MODEL midas_model (const int *list, const char *param,
		   DATASET *dset, gretlopt opt,
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
    MODEL mod;
    int err = 0;

    gretl_model_init(&mod, dset);

    if (param == NULL || *param == '\0') {
	err = E_DATA;
    } else {
	err = parse_midas_specs(param, dset, &minfo, &nmidas);
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

    if (!err) {
	/* determine usable sample range */
	err = midas_set_sample(list, dset, minfo, nmidas);
    }

    if (!err) {
	/* add the required matrices */
	gretl_matrix *m;

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
		if (minfo[i].type != MIDAS_ALMON) {
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
    }

    if (!err) {
	char line[MAXLEN];
	int j = 0;

	/* initial "nls" line */
	sprintf(line, "nls %s = ", dset->varname[list[1]]);
	if (xlist != NULL) {
	    strcat(line, "lincomb(XL___, bx___)");
	}
	for (i=0; i<nmidas; i++) {
	    if (minfo[i].type == MIDAS_ALMON) {
		sprintf(tmp, " + mlc___%d", i+1);
	    } else {
		sprintf(tmp, " + bmlc___[%d]*mlc___%d", ++j, i+1);
	    }
	    strcat(line, tmp);
	}
	err = nl_parse_line(NLS, line, dset, prn);

	/* MIDAS series and gradient matrices */
	for (i=0; i<nmidas && !err; i++) {
	    sprintf(line, "series mlc___%d = mlincomb(%s, %s, %d)",
		    i+1, minfo[i].lname, minfo[i].mname, minfo[i].type);
	    err = nl_parse_line(NLS, line, dset, prn);
	    if (!err) {
		sprintf(line, "matrix mgr___%d = mgradient(%d, %s, %d)",
			i+1, minfo[i].nterms, minfo[i].mname, minfo[i].type);
		err = nl_parse_line(NLS, line, dset, prn);
	    }
	}

	/* derivatives */
	if (!err && xlist != NULL) {
	    strcpy(line, "deriv bx___ = MX___");
	    err = nl_parse_line(NLS, line, dset, prn);
	}
	if (!err && hfslopes > 0) {
	    strcpy(line, "deriv bmlc___ = {");
	    for (i=0; i<nmidas; i++) {
		if (minfo[i].type != MIDAS_ALMON) {
		    sprintf(tmp, "mlc___%d", i+1);
		    strcat(line, tmp);
		}
		for (j=i; j<nmidas; j++) {
		    if (minfo[j].type != MIDAS_ALMON) {
			strcat(line, ", ");
			break;
		    }
		}
	    }
	    strcat(line, "}");
	    err = nl_parse_line(NLS, line, dset, prn);
	}
	j = 0;
	for (i=0; i<nmidas && !err; i++) {
	    if (minfo[i].type != MIDAS_ALMON) {
		sprintf(line, "deriv %s = bmlc___[%d] * {%s} * mgr___%d",
			minfo[i].mname, ++j, minfo[i].lname, i+1);
	    } else {
		sprintf(line, "deriv %s = {%s} * mgr___%d",
			minfo[i].mname, minfo[i].lname, i+1);
	    }		
	    err = nl_parse_line(NLS, line, dset, prn);
	}

	/* parameter names */
	if (!err) {
	    strcpy(line, "param_names \"");
	    if (xlist != NULL) {
		for (i=1; i<=xlist[0]; i++) {
		    strcat(line, dset->varname[xlist[i]]);
		    strcat(line, " ");
		}
	    }
	    for (i=0; i<nmidas && !err; i++) {
		if (minfo[i].type != MIDAS_ALMON) {
		    if (hfslopes > 1) {
			sprintf(tmp, " HF_slope%d", i+1);
		    } else {
			strcpy(tmp, " HF_slope");
		    }
		    strcat(line, tmp);
		}
	    }
	    for (i=0; i<nmidas && !err; i++) {
		for (j=0; j<minfo[i].k; j++) {
		    sprintf(tmp, " %s[%d]", minfo[i].mname, j+1);
		    strcat(line, tmp);
		}
	    }
	    strcat(line, "\"");
	    fprintf(stderr, "%s\n", line);
	    err = nl_parse_line(NLS, line, dset, prn);
	}
    }

    if (!err) {
	mod = nl_model(dset, opt & OPT_G, prn);
    }

    for (i=0; i<nmidas; i++) {
	user_var_delete_by_name(minfo[i].lname, NULL);
	sprintf(tmp, "mgr___%d", i+1);
	user_var_delete_by_name(tmp, NULL);
    }

    destroy_private_uvars();
    free(minfo);
    free(xlist);

    if (dset->v > origv) {
	dataset_drop_last_variables(dset, dset->v - origv);
    }

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	mod.ci = MIDASREG;
	free(mod.depvar);
	mod.depvar = gretl_strdup(dset->varname[list[1]]);
    }

    return mod;
}
