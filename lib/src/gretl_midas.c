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

/* Note, 2016-09-26: to avoid clutter and unclarity, this version of
   gretl_midas.c defaults to "conditional OLS" for estimation of MIDAS
   models that comprise one or more beta or normalized exponential
   Almon (nealmon) terms. By this we mean a combination of L-BFGS-B
   (with constraints on the beta and/or nealmon hyper-parameters) and
   OLS: all coefficients other than the hyper-parameters are estimated
   via OLS conditional on the "theta" vector of hyper-parameters as
   optimized by L-BFGS-B.

   In case this proves to be a bad idea, we could go back to commit
   8aa2c9, the step before we started working towards conditional OLS.
*/

#include "libgretl.h"
#include "usermat.h"
#include "uservar.h"
#include "matrix_extra.h"
#include "libset.h"
#include "gretl_bfgs.h"
#include "nlspec.h"
#include "gretl_midas.h"

#define MIDAS_DEBUG 0
#define FC_DEBUG 0

enum mflags {
    M_PRELAG  = 1 << 0,
    M_AUTO    = 1 << 1,
    M_TMPLIST = 1 << 2
};

enum midas_methods {
    MDS_NLS,
    MDS_OLS,
    MDS_BFGS,
    MDS_GSS
};

/* information on an individual MIDAS term */

struct midas_term_ {
    char lnam0[VNAMELEN];  /* name of MIDAS list on input */
    char lname[VNAMELEN];  /* name of MIDAS list */
    char mname[VNAMELEN];  /* name of initial theta vector */
    gretl_matrix *theta;   /* value of initial theta vector */
    char flags;            /* values from @mflags */
    int minlag;            /* minimum lag */
    int maxlag;            /* maximum lag */
    int type;              /* type of parameterization */
    int nparm;             /* number of parameters */
    int nlags;             /* number of lag terms */
    int *laglist;          /* list of lag series */
};

typedef struct midas_term_ midas_term;

/* information on a MIDAS model as a whole */

struct midas_info_ {
    const int *list;      /* incoming regression list */
    int yno;              /* ID number of dependent variable */
    int nx;               /* number of low-frequency regressors */
    midas_term *mterms;   /* array of MIDAS terms */
    int nmidas;           /* number of elements in @mterms */
    int *xlist;           /* list of low-frequency regressors */
    int *seplist;         /* list of coeff separator positions */
    gchar *pnames;        /* string holding parameter names */
    int hfslopes;         /* # of high-frequency slope coeffs */
    int nalmonp;          /* number of straight Almon poly terms */
    int nbeta1;           /* number of clamped beta terms */
    int method;           /* estimation method */
    int ldepvar;          /* position of lagged dependent variable? */
    int nobs;             /* number of observations used */
    size_t colsize;       /* matrix column size in bytes */
    DATASET *dset;        /* dataset */
    gretl_matrix *b;      /* all, or OLS, coefficients */
    gretl_matrix *g;      /* SSR gradients */
    gretl_matrix *theta;  /* MIDAS hyper-params, combined */
    gretl_matrix *bounds; /* bounds on hyper-params */
    gretl_matrix *y;      /* dependent vector */
    gretl_matrix *u;      /* residual vector */
    gretl_matrix *xi;     /* single column of X */
    gretl_matrix *X;      /* X(\beta) */
    double SSR;           /* sum of squared residuals */
    double orig_toler;    /* prior NLS_TOLER value */
};

typedef struct midas_info_ midas_info;

#define prelag(m)     (m->flags & M_PRELAG)
#define auto_theta(m) (m->flags & M_AUTO)

static midas_term *mterms_from_array (gretl_array *A,
				      int *nmidas,
				      int *err);

static double cond_ols_callback (double *theta, double *g,
				 int n, void *ptr);

/* Identify parameterizations which take an additional
   leading coefficient: all but U-MIDAS and the plain
   Almon polynomial specification.
*/
#define takes_coeff(t) (t != MIDAS_U && t != MIDAS_ALMONP)

/* convenience macro for picking out beta specifications */
#define beta_type(t) (t == MIDAS_BETA0 || t == MIDAS_BETA1 || t == MIDAS_BETAN)

/* days per month or quarter: maybe make this user-
   configurable? */

int midas_days_per_period (int days_per_week, int pd)
{
    int ret = 0;

    if (days_per_week == 5) {
	ret = 22;
    } else if (days_per_week == 6) {
	ret = 26;
    } else if (days_per_week == 7) {
	ret = 30;
    }

    return (pd == 12)? ret : 3 * ret;
}

static int midas_days_to_freq (int pd, int m)
{
    int days[] = {22, 26, 30};
    int freq[] = {5, 6, 7};
    int i;

    for (i=0; i<3; i++) {
	if ((pd == 12 && m == days[i]) ||
	    (pd == 4 && m == 3 * days[i])) {
	    return freq[i];
	}
    }

    return 0;
}

/* Could @m be a valid frequency ratio (= number of members
   of a valid "MIDAS list"), given the periodicity of @dset?
   If so, return 1; if not, return 0.
*/

int is_valid_midas_frequency_ratio (const DATASET *dset, int m)
{
    if (dset->pd == 1) {
	/* lf = annual, hf = quarterly or monthly */
	return (m == 4 || m == 12);
    } else if (dset->pd == 4 && m == 3) {
	/* lf = quarterly, hf = monthly */
	return 1;
    } else if (dset->pd == 4 || dset->pd == 12) {
	/* lf = quarterly or monthly, hf = daily */
	if (m == midas_days_per_period(5, dset->pd)) {
	    return 1;
	} else if (m == midas_days_per_period(6, dset->pd)) {
	    return 1;
	} else if (m == midas_days_per_period(7, dset->pd)) {
	    return 1;
	}
    }

    return 0;
}

int get_midas_frequency (const DATASET *dset, int m)
{
    if (dset->pd == 1 && (m == 4 || m == 12)) {
	/* lf = annual -> hf quarterly or monthly */
	return m;
    } else if (dset->pd == 4 && m == 3) {
	/* lf = quarterly -> hf monthly */
	return 12;
    } else if (dset->pd == 4 || dset->pd == 12) {
	/* hf daily 5, 6 or 7? */
	return midas_days_to_freq(dset->pd, m);
    }

    return 0;
}

const char *midas_pdstr (const DATASET *dset, int cfac)
{
    const char *pdstrs[] = {
	N_("quarterly"),
	N_("monthly"),
	N_("daily")
    };
    int f = get_midas_frequency(dset, cfac);

    if (f == 4) {
	return pdstrs[0];
    } else if (f == 12) {
	return pdstrs[1];
    } else if (f == 5 || f == 6 || f == 7) {
	return pdstrs[2];
    } else {
	return "unknown frequency";
    }
}

int midas_m_from_pd (const DATASET *dset, int pd)
{
    if (dset->pd == 1 && (pd == 4 || pd == 12)) {
	return pd;
    } else if (dset->pd == 4 && pd == 12) {
	return 3;
    } else if (dset->pd == 4 || dset->pd == 12) {
	return midas_days_per_period(pd, dset->pd);
    }

    return 0;
}

/* Infer the current month from the current @qtr along
   with the number of days per period, @ndays, and the
   current index within the month-days array, @day.
*/

static int quarter_to_month (int qtr, int ndays, int day)
{
    return qtr * 3 - 2 + (day - 1) / (ndays/3);
}

static int get_midas_pd (int pd, int m, int *err)
{
    int mpd = 0;

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

    return mpd;
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

    if (!*err) {
	mpd = get_midas_pd(pd, m, err);
    }

    if (*err) {
	return NULL;
    }

    if (!gretl_is_midas_list(list, dset)) {
	gretl_warnmsg_set(_("The argument does not seem to be a MIDAS list"));
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

	ntolabel(obs, dset->t1, dset);

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
		ntolabel(obs, t, dset);
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
	    ntolabel(mset->endobs, mset->t2, mset);
	}
    }

    return mset;
}

static int midas_list_check (const int *list,
			     const DATASET *dset)
{
    int mpd, m = list[0];
    int err = 0;

    if (m < 3 || gretl_list_has_separator(list)) {
	return E_INVARG;
    } else if (!dataset_is_time_series(dset)) {
	return E_PDWRONG;
    } else if (dset->pd != 1 && dset->pd != 4) {
	/* host dataset should be annual or quarterly */
	return E_PDWRONG;
    }

    mpd = get_midas_pd(dset->pd, m, &err);
    if (!err && mpd != 3 && mpd != 4 && mpd != 12) {
	err = E_INVARG;
    }

    return err;
}

gretl_matrix *midas_list_to_vector (const int *list,
				    const DATASET *dset,
				    int *err)
{
    gretl_matrix *mvec;
    int t, i, vi, n_tail = 0;
    int T, m = list[0];

    *err = midas_list_check(list, dset);
    if (*err) {
	return NULL;
    }

    if (dset->t2 < dset->n - 1) {
	t = dset->t2 + 1;
	for (i=m; i>0; i--) {
	    vi = list[i];
	    if (!na(dset->Z[vi][t])) {
		n_tail++;
	    } else {
		break;
	    }
	}
    }

    T = sample_size(dset) * m + n_tail;
    mvec = gretl_matrix_alloc(T, 1);

    if (mvec == NULL) {
	*err = E_ALLOC;
    } else {
	int s = 0;

	for (t=dset->t1; t<=dset->t2; t++) {
	    /* read data right-to-left! */
	    for (i=m; i>0; i--) {
		vi = list[i];
		mvec->val[s++] = dset->Z[vi][t];
	    }
	}
	if (n_tail > 0) {
	    t = dset->t2 + 1;
	    for (i=0; i<n_tail; i++) {
		vi = list[m-i];
		mvec->val[s++] = dset->Z[vi][t];
	    }
	}
    }

    return mvec;
}

static gretl_matrix *make_auto_theta (char *mstr, int i,
				      int ptype, int k,
				      int m1, int m2,
				      int *err)
{
    gretl_matrix *theta = NULL;

    if (k == 0) {
	/* we have to infer k */
	if (*mstr == '\0' || !strcmp(mstr, "null")) {
	    /* OK if we know how many parameters are needed? */
	    if (ptype == MIDAS_BETA0 || ptype == MIDAS_BETA1) {
		k = 2;
	    } else if (ptype == MIDAS_BETAN) {
		k = 3;
	    } else if (ptype == MIDAS_NEALMON) {
		/* we'll default to 2 */
		k = 2;
	    } else if (ptype == MIDAS_U) {
		k = m2 - m1 + 1;
	    } else {
		*err = E_INVARG;
	    }
	} else if (integer_string(mstr)) {
	    /* does @mstr give the number of params? */
	    int chk = atoi(mstr);

	    if (chk >= 1 && chk < 100) {
		k = chk;
	    }
	} else {
	    /* try treating @mstr as matrix definition? */
	    theta = generate_matrix(mstr, NULL, err);
	}
    }

    if (k > 0) {
	theta = gretl_zero_matrix_new(1, k);
	if (theta != NULL) {
	    if (beta_type(ptype)) {
		theta->val[0] = 1;
		theta->val[1] = 10; /* optimize somehow? */
	    }
	} else {
	    *err = E_ALLOC;
	}
    }

    if (theta != NULL) {
	sprintf(mstr, "theta___%d", i+1);
	private_matrix_add(theta, mstr);
    } else if (!*err) {
	*err = E_PARSE;
    }

#if MIDAS_DEBUG
    gretl_matrix_print(theta, "auto-generated theta");
#endif

    return theta;
}

static int lag_info_from_prelag_list (midas_term *mt,
				      const int *list,
				      const DATASET *dset)
{
    int mf = series_get_midas_freq(dset, list[1]);
    int m1 = series_get_lag(dset, list[1]);
    int p1 = series_get_midas_period(dset, list[1]);
    int i, p, maxp = 0;

    if (mf > 0) {
	maxp = midas_m_from_pd(dset, mf);
    }

    if (maxp == 0 && p1 > 0) {
	for (i=1; i<=list[0]; i++) {
	    p = series_get_midas_period(dset, list[i]);
	    if (p > maxp) {
		maxp = p;
	    }
	}
    }

    if (is_valid_midas_frequency_ratio(dset, maxp)) {
	int hfl = maxp * m1 - (p1 - 1);

	mt->minlag = hfl;
	mt->maxlag = hfl + list[0] - 1;
    } else {
	/* oof! just report 1,... */
	mt->minlag = 1;
	mt->maxlag = list[0];
    }

    return 0;
}

static void midas_term_init (midas_term *mt)
{
    mt->lnam0[0] = '\0';
    mt->lname[0] = '\0';
    mt->mname[0] = '\0';
    mt->theta = NULL;
    mt->flags = 0;
    mt->minlag = 0;
    mt->maxlag = 0;
    mt->type = 0;
    mt->nparm = 0;
    mt->nlags = 0;
}

static int get_p1_p2 (const char *s1, const char *s2,
		      int *p1, int *p2)
{
    int err = 0;

    *p1 = gretl_int_from_string(s1, &err);

    if (!err) {
	*p2 = gretl_int_from_string(s2, &err);
    }

    return err;
}

static int midas_type_from_string (const char *s)
{
    if (!strcmp(s, "\"umidas\"") || !strcmp(s, "umidas")) {
	return MIDAS_U;
    } else if (!strcmp(s, "\"nealmon\"") || !strcmp(s, "nealmon")) {
	return MIDAS_NEALMON;
    } else if (!strcmp(s, "\"beta0\"") || !strcmp(s, "beta0")) {
	return MIDAS_BETA0;
    } else if (!strcmp(s, "\"betan\"") || !strcmp(s, "betan")) {
	return MIDAS_BETAN;
    } else if (!strcmp(s, "\"almonp\"") || !strcmp(s, "almonp")) {
	return MIDAS_ALMONP;
    } else if (!strcmp(s, "\"beta1\"") || !strcmp(s, "beta1")) {
	return MIDAS_BETA1;
    } else {
	return -1;
    }
}

static int midas_type_from_arg (const char *s, int *umidas, int *err)
{
    int ret = -1;

    if (integer_string(s)) {
	ret = atoi(s);
    } else {
	ret = midas_type_from_string(s);
	if (ret < 0) {
	    char *sval = get_string_by_name(s);

	    if (sval != NULL) {
		ret = midas_type_from_string(sval);
	    }
	}
	if (ret < 0 && gretl_is_scalar(s)) {
	    ret = gretl_scalar_get_value(s, err);
	}
    }

    if (ret < MIDAS_U || ret >= MIDAS_MAX) {
	*err = E_INVARG;
    } else if (ret == 0) {
	*umidas = 1;
    }

    return ret;
}

/* Parse a particular entry in the incoming array of MIDAS
   specifications. Each entry should look like one of
   the following:

   (1) mds(list, minlag, maxlag, type, theta)
   (2) mds(list, minlag, maxlag, 0)

   (3) mdsl(list, type, theta)
   (4) mdsl(list, 0)

   Forms (1) and (2) imply that lags of the MIDAS terms
   in @list should be auto-generated, while (3) and (4)
   imply that the incoming @list already holds whatever
   lags are wanted.

   Variants (2) and (4), with @type set to 0, are for
   U-MIDAS terms: in this case @theta (a vector to
   initialize the coefficients) is not required.
*/

static int parse_midas_term (const char *s,
			     midas_term *mt,
			     int i,
			     const DATASET *dset)
{
    char lname[VNAMELEN];
    char mname[VNAMELEN];
    char p1str[VNAMELEN];
    char p2str[VNAMELEN];
    char typestr[10];
    char fmt[64];
    int ns, p1 = 0, p2 = 0;
    int type, umidas = 0;
    int err = 0;

    midas_term_init(mt);
    *typestr = *lname = *mname = '\0';

    if (!strncmp(s, "mds(", 4)) {
	/* calling for auto-generated lags */
	s += 4;
	sprintf(fmt, "%%%d[^, ] , %%%d[^, ] , %%%d[^, ] , %%%d[^,) ] , %%%d[^) ])",
		VNAMELEN-1, VNAMELEN-1, VNAMELEN-1, 9, VNAMELEN-1);
	ns = sscanf(s, fmt, lname, p1str, p2str, typestr, mname);
	err = get_p1_p2(p1str, p2str, &p1, &p2);
	if (!err) {
	    type = midas_type_from_arg(typestr, &umidas, &err);
	}
	if (!err && ns < 4) {
	    err = E_PARSE;
	}
    } else if (!strncmp(s, "mdsl(", 5)) {
	/* list already holds lags */
	mt->flags |= M_PRELAG;
	s += 5;
	sprintf(fmt, "%%%d[^, ] , %%%d[^,) ] , %%%d[^) ])",
		VNAMELEN-1, 9, VNAMELEN-1);
	ns = sscanf(s, fmt, lname, typestr, mname);
	type = midas_type_from_arg(typestr, &umidas, &err);
	if (!err && ns < 2) {
	    err = E_PARSE;
	}
    } else {
	err = E_INVARG;
    }

    if (!err) {
	int *list = get_list_by_name(lname);
	int k = 0;

	if (!umidas) {
	    gretl_matrix *theta = NULL;
	    char c = mname[0];

	    if (c != '\0' && !isdigit(c) && c != '{' && strcmp(mname, "null")) {
		theta = get_matrix_by_name(mname);
	    }

	    if (theta != NULL) {
		/* copy, don't clobber, the user's initializer */
		sprintf(mname, "theta___%d", i+1);
		mt->theta = gretl_matrix_copy(theta);
		if (mt->theta != NULL) {
		    private_matrix_add(mt->theta, mname);
		}
	    } else {
		/* create automatic initializer */
		mt->flags |= M_AUTO;
		mt->theta = make_auto_theta(mname, i, type, 0,
					    p1, p2, &err);
	    }
	}

	if (err) {
	    ; /* leave it alone */
	} else if (!umidas && mt->theta == NULL) {
	    err = E_ALLOC;
	} else if (prelag(mt) && list == NULL) {
	    err = E_INVARG;
	} else if (!prelag(mt) && !gretl_is_midas_list(list, dset)) {
	    gretl_errmsg_set("mds(): the first term must be a MIDAS list");
	    err = E_INVARG;
	} else if (p1 > p2) {
	    err = E_INVARG;
	} else if (type < 0 || type >= MIDAS_MAX) {
	    err = E_INVARG;
	} else if (umidas) {
	    if (prelag(mt)) {
		k = list[0];
	    } else {
		k = p2 - p1 + 1;
	    }
	} else {
	    k = gretl_vector_get_length(mt->theta);
	    if (k < 1 || (type == MIDAS_BETA0 && k != 2) ||
		(type == MIDAS_BETA1 && k != 2) ||
		(type == MIDAS_BETAN && k != 3)) {
		err = E_INVARG;
	    }
	}

	if (!err) {
	    strcpy(mt->lnam0, lname);
	    strcpy(mt->lname, lname);
	    if (!umidas) {
		strcpy(mt->mname, mname);
	    }
	    if (prelag(mt)) {
		/* scrounge lag info from incoming list */
		lag_info_from_prelag_list(mt, list, dset);
	    } else {
		mt->minlag = p1;
		mt->maxlag = p2;
	    }
	    mt->type = type;
	    mt->nparm = k;
	}
    }

    return err;
}

/* In case we got any U-MIDAS terms in the specification,
   check to see if we need to add any initializers for
   the associated coefficients.
*/

static int umidas_check (midas_info *mi, int n_umidas)
{
    int i, err = 0;

    if (n_umidas == mi->nmidas) {
	/* all U-MIDAS: so use OLS */
	mi->method = MDS_OLS;
	return 0;
    }

    /* mix of U-MIDAS and "unboxed" terms: we'll be using
       NLS and we need initializers */

    for (i=0; i<mi->nmidas && !err; i++) {
	midas_term *mt = &mi->mterms[i];

	if (mt->type == MIDAS_U && mt->theta == NULL) {
	    mt->flags |= M_AUTO;
	    mt->theta = make_auto_theta(mt->mname, i, MIDAS_U,
					mt->nparm, 0, 0, &err);
	}
    }

    return err;
}

/* Parse the @spec string, which should contain one or more
   MIDAS specifications. For details on what exactly we're
   looking for, see the comment on parse_midas_term() above.
   This function also potentially revises the estimation
   method: to OLS if we have nothing but UMIDAS terms; to
   MDS_BFGS if we "boxed" terms and the --levenberg option
   has not been specified.
*/

static int
parse_midas_specs (midas_info *mi, const char *spec,
		   const DATASET *dset, gretlopt opt)
{
    const char *s;
    int n_spec = 0;
    int n_umidas = 0;
    int n_boxed = 0;
    int n_beta1 = 0;
    int n_almonp = 0;
    int err = 0;

    /* first check: count closing parentheses */
    s = spec;
    while (*s) {
	if (*s == ')') {
	    n_spec++;
	}
	s++;
    }

    if (n_spec == 0) {
	/* spec is junk! */
	err = E_PARSE;
    } else {
	/* allocate info structs */
	mi->mterms = malloc(n_spec * sizeof *mi->mterms);
	if (mi->mterms == NULL) {
	    err = E_ALLOC;
	}
    }

    /* number of slope coeffs on high-frequency terms */
    mi->hfslopes = 0;

    if (!err) {
	/* parse and record individual MIDAS specs */
	char test[128];
	const char *p;
	int len, i = 0;

	s = spec;
	for (i=0; i<n_spec && !err; i++) {
	    midas_term *mt = &mi->mterms[i];

	    while (*s == ' ') s++;
	    p = s;
	    while (*p && *p != ')') p++;
	    len = p - s + 1;
	    if (len > 127) {
		err = E_PARSE;
	    } else {
		*test = '\0';
		strncat(test, s, len);
		err = parse_midas_term(test, mt, i, dset);
		if (!err) {
		    if (mt->type == MIDAS_BETA0 && (opt & OPT_C)) {
			/* support legacy option */
			mt->type = MIDAS_BETA1;
		    }
		    if (mt->type == MIDAS_U) {
			n_umidas++;
		    } else if (beta_type(mt->type) ||
			       mt->type == MIDAS_NEALMON) {
			n_boxed++;
		    } else if (mt->type == MIDAS_ALMONP) {
			n_almonp++;
		    }
		    if (takes_coeff(mt->type)) {
			mi->hfslopes += 1;
		    }
		    if (mt->type == MIDAS_BETA1) {
			n_beta1++;
		    }
		}
		s = p + 1;
	    }
	}
    }

    if (!err && n_beta1 > 0) {
	/* check validity of beta clamping */
	if (n_almonp > 0) {
	    gretl_errmsg_set("Can't combine beta1 with almonp terms");
	    err = E_INVARG;
	} else if (opt & OPT_L) {
	    /* can't combine beta1 and --levenberg */
	    err = E_BADOPT;
	}
    }

    if (!err) {
	mi->nmidas = n_spec;
	mi->nalmonp = n_almonp;
	mi->nbeta1 = n_beta1;
	if (n_boxed > 0 && n_almonp == 0) {
	    /* prefer LBFGS, but respect OPT_L if it's given */
	    if (!(opt & OPT_L)) {
		/* ! --levenberg */
		mi->method = MDS_BFGS;
	    }
	} else if (n_umidas > 0) {
	    /* note: switch to MDS_OLS if appropriate */
	    err = umidas_check(mi, n_umidas);
	}
    }

    if (err) {
	/* scrub all this stuff to avoid a crash on cleanup */
	free(mi->mterms);
	mi->mterms = NULL;
	mi->nmidas = 0;
	mi->nalmonp = 0;
	mi->nbeta1 = 0;
    }

    return err;
}

/* Extract the list of regular (low-frequency) regressors.
   While we're at it, see if any regressors are lags of the
   dependent variable and if so record this fact in the
   ldepvar member of @mi.
*/

static int make_midas_xlist (midas_info *mi)
{
    int ldvpos = 0, ldv_count = 0;
    int i, xi, err = 0;

    mi->nx = mi->list[0] - 1;

    if (mi->nx > 0) {
	mi->xlist = gretl_list_new(mi->nx);
	if (mi->xlist == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=2; i<=mi->list[0]; i++) {
		xi = mi->list[i];
		mi->xlist[i-1] = xi;
		if (standard_lag_of(xi, mi->yno, mi->dset) > 0) {
		    ldvpos = i;
		    ldv_count++;
		}
	    }
	}
    }

    if (ldv_count == 1) {
	mi->ldepvar = ldvpos;
    }

    return err;
}

/* Given the name of an incoming MIDAS list plus minlag and
   maxlag values (at this point stored in the midas_term
   structure, @mt) build the list of required lags of the
   MIDAS series. Or in case the incoming list already
   includes the required lags, just take a pointer to it.
*/

static int *make_midas_laglist (midas_term *mt,
				DATASET *dset,
				int *err)
{
    int *list = get_list_by_name(mt->lname);

    if (list == NULL) {
	fprintf(stderr, "make_midas_laglist, '%s': no list!\n",
		mt->lname);
	*err = E_DATA;
	return NULL;
    }

    if (!prelag(mt) && mt->minlag == 1 && mt->maxlag == list[0]) {
	/* the incoming list was not flagged as "mdsl", and so
	   must be a MIDAS list, but it nonetheless contains all
	   the lags we need
	*/
	mt->flags |= M_PRELAG;
    }

    if (prelag(mt)) {
	/* don't copy the list (and don't free it either!) */
	return list;
    } else {
	/* copy, because we're going to modify the list */
	int *lcpy = gretl_list_copy(list);

	if (lcpy == NULL) {
	    *err = E_ALLOC;
	} else {
	    *err = list_laggenr(&lcpy, mt->minlag, mt->maxlag,
				NULL, dset, lcpy[0], OPT_L);
	}

	return lcpy;
    }
}

static int fcast_make_midas_laglists (midas_term *mterms,
				      int nmidas,
				      DATASET *dset)
{
    midas_term *mt;
    int *mlist;
    int i, err = 0;

    for (i=0; i<nmidas && !err; i++) {
	mt = &mterms[i];
	mlist = make_midas_laglist(mt, dset, &err);
	if (!err) {
	    if (!prelag(mt)) {
		sprintf(mt->lname, "ML___%d", i+1);
		err = remember_list(mlist, mt->lname, NULL);
		mt->flags |= M_TMPLIST;
	    }
	    mt->nlags = mlist[0];
	    mt->laglist = mlist;
	}
    }

    return err;
}

static int make_midas_laglists (midas_info *mi,
				DATASET *dset)
{
    midas_term *mt;
    int *mlist;
    int i, err = 0;

    for (i=0; i<mi->nmidas && !err; i++) {
	mt = &mi->mterms[i];
	mlist = make_midas_laglist(mt, dset, &err);
	if (!err) {
	    /* In the "prelag" case the laglist is already
	       in userspace. If method is BFGS the list
	       doesn't have to be pushed into userspace,
	       but we'll record it nonetheless for model
	       printout purposes.
	    */
	    if (!prelag(mt)) {
		sprintf(mt->lname, "ML___%d", i+1);
		/* note: remember_list copies its first arg */
		err = remember_list(mlist, mt->lname, NULL);
		mt->flags |= M_TMPLIST;
	    }
	    mt->nlags = mlist[0];
	    mt->laglist = mlist;
	}
    }

    return err;
}

/* If @list is non-NULL, build a full list of all series
   involved: dependent variable, regular regressors, and
   all lags of MIDAS terms. We want this either for setting
   the usable sample range, or (in the case of U-MIDAS via
   OLS), as the list to pass to lsq().

   If @list is NULL, however, the returned list will just
   contain all the MIDAS lag terms: we use this variant
   when setting up for forecasting.
*/

static int *make_midas_biglist (const int *list,
				midas_term *mterms,
				int nmidas)
{
    int i, j, nt = 0;
    int *biglist = NULL;

    if (list != NULL) {
	nt = list[0];
    }

    for (j=0; j<nmidas; j++) {
	nt += mterms[j].nlags;
    }

    biglist = gretl_list_new(nt);

    if (biglist != NULL) {
	int k = 1;

	if (list != NULL) {
	    for (i=1; i<=list[0]; i++) {
		biglist[k++] = list[i];
	    }
	}
	for (j=0; j<nmidas; j++) {
	    for (i=1; i<=mterms[j].nlags; i++) {
		biglist[k++] = mterms[j].laglist[i];
	    }
	}
    }

    return biglist;
}

/* If we're doing MIDAS via nls, set the usable
   sample range first. That way we'll know how
   big the X data matrix ought to be.
*/

static int midas_set_sample (midas_info *mi,
			     DATASET *dset)
{
    int *biglist;
    int err = 0;

    biglist = make_midas_biglist(mi->list, mi->mterms, mi->nmidas);
    if (biglist == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	int t1 = dset->t1, t2 = dset->t2;

	err = list_adjust_sample(biglist, &t1, &t2, dset, NULL);
	if (!err) {
	    dset->t1 = t1;
	    dset->t2 = t2;
	    mi->nobs = t2 - t1 + 1;
	    mi->colsize = mi->nobs * sizeof(double);
	}
    }

    free(biglist);

#if MIDAS_DEBUG
    fprintf(stderr, "midas_set_sample: returning %d\n", err);
#endif

    return err;
}

/* estimate unrestricted MIDAS via OLS */

static int umidas_ols (MODEL *pmod,
		       midas_info *mi,
		       DATASET *dset,
		       gretlopt opt)
{
    int *biglist;

    biglist = make_midas_biglist(mi->list, mi->mterms, mi->nmidas);

    if (biglist == NULL) {
	return E_ALLOC;
    } else {
	*pmod = lsq(biglist, dset, OLS, opt | OPT_Z);
	free(biglist);
    }

    return pmod->errcode;
}

/* For forecasting, put into uservar space a "private"
   matrix containing all the coefficients on MIDAS
   lag terms.
*/

static int push_midas_coeff_array (const MODEL *pmod)
{
    gretl_matrix *mhfb, *hfb;
    int err = 0;

    mhfb = gretl_model_get_data(pmod, "hfb");
    if (mhfb == NULL) {
	fprintf(stderr, "Couldn't retrieve gross midas coeffs\n");
	err = E_DATA;
    } else {
	hfb = gretl_matrix_copy(mhfb);
	if (hfb == NULL) {
	    err = E_ALLOC;
	}
    }

#if FC_DEBUG
    gretl_matrix_print(hfb, "hfb in fcast");
#endif

    if (!err) {
	err = private_matrix_add(hfb, "hfb___");
    }

    return err;
}

/* Prepare for generating a MIDAS forecast: this
   is called from midas_fcast() in forecast.c.
*/

int midas_forecast_setup (const MODEL *pmod,
			  DATASET *dset,
			  ForecastMethod method,
			  char **pformula)
{
    gretl_array *mA;
    midas_term *mterms = NULL;
    int *xlist = NULL;
    int *hflist = NULL;
    int nmidas = 0;
    int nx = 0;
    int err = 0;

    mA = gretl_model_get_data(pmod, "midas_info");
    if (mA == NULL) {
	err = E_DATA;
    }

    if (!err && gretl_model_get_int(pmod, "no_lfx") == 0) {
	xlist = gretl_model_get_list(pmod, "lfxlist");
	if (xlist == NULL) {
	    err = E_DATA;
	} else {
	    nx = xlist[0];
	}
    }

    if (!err) {
	mterms = mterms_from_array(mA, &nmidas, &err);
    }

    if (!err) {
	/* reconstitute MIDAS lag-lists */
	err = fcast_make_midas_laglists(mterms, nmidas, dset);
    }

    if (!err) {
	/* build and push list of all MIDAS terms */
	hflist = make_midas_biglist(NULL, mterms, nmidas);
	if (hflist == NULL) {
	    err = E_ALLOC;
	} else {
	    /* note: remember_list() copies its argument */
#if FC_DEBUG
	    int j;
	    printlist(hflist, "hflist");
	    for (j=1; j<=hflist[0]; j++) {
		fprintf(stderr, "%d %s\n", hflist[j], dset->varname[hflist[j]]);
	    }
#endif
	    err = remember_list(hflist, "HFL___", NULL);
	    user_var_privatize_by_name("HFL___");
	    free(hflist);
	}
    }

    if (!err) {
	/* retrieve and push vector of all MIDAS coeffs */
	err = push_midas_coeff_array(pmod);
    }

    if (!err) {
	/* build a string that can be passed to "genr" to
	   calculate fitted values */
	char tmp[64], line[MAXLEN];
	double *b = pmod->coeff;
	int p, yno = pmod->list[1];
	int i, xi, j = 0;

	*line = '\0';
	gretl_push_c_numeric_locale();

	for (i=0; i<nx; i++) {
	    xi = xlist[i+1];
	    if (xi == 0) {
		sprintf(tmp, "%+.15g", b[j++]);
	    } else {
		/* allow for dynamic formula? */
		p = (method == FC_STATIC)? 0 :
		    standard_lag_of(xi, yno, dset);
		if (p > 0) {
		    sprintf(tmp, "%+.15g*%s(-%d)", b[j++],
			    dset->varname[yno], p);
		} else {
		    sprintf(tmp, "%+.15g*%s", b[j++],
			    dset->varname[xi]);
		}
	    }
	    strcat(line, tmp);
	}
	strcat(line, "+lincomb(HFL___,hfb___)");

	gretl_pop_c_numeric_locale();

#if FC_DEBUG
	fprintf(stderr, "formula='%s'\n", line);
#endif
	*pformula = gretl_strdup(line);
    }

    free(mterms);

    return err;
}

static int midas_beta_init (midas_info *mi)
{
    midas_term *mt = &mi->mterms[0];
    int i;

    if (!auto_theta(mt)) {
	/* the user provided a theta value, so use it */
	for (i=0; i<mt->nparm; i++) {
	    mi->theta->val[i] = mt->theta->val[i];
	}
    } else {
	/* Trial some alternative initializations for the
	   beta hyper-parameters (theta) and select the one
	   that produces the smallest SSR on running OLS.
	*/
	double theta[3][2] = {
	    {1, 1}, {1, 10}, {2, 10}
	};
	double SSR, SSRmin = 1e200;
	int imax = mt->type == MIDAS_BETA1 ? 2 : 3;
	int best_idx = -1;

	for (i=0; i<imax; i++) {
	    SSR = cond_ols_callback(theta[i], NULL, 2, mi);
	    if (SSR < SSRmin) {
		SSRmin = SSR;
		best_idx = i;
	    }
	}

	if (best_idx >= 0) {
	    mi->theta->val[0] = theta[best_idx][0];
	    mi->theta->val[1] = theta[best_idx][1];
#if MIDAS_DEBUG
	    fprintf(stderr, "midas_beta_init: best_idx = %d, SSRmin=%g\n",
		    best_idx, SSRmin);
	    gretl_matrix_print(mi->theta, "best theta");
#endif
	}
    }

    return 0;
}

/* L-BFGS-B apparatus */

static int midas_bfgs_setup (midas_info *mi, DATASET *dset,
			     gretlopt opt)
{
    double eps = pow(2.0, -52);
    double *src, *targ = NULL;
    int i, j, k, ii, vi;
    int nb, bound_rows = 0;
    int type0 = 0;
    midas_term *mt;

    if (mi->nbeta1 == 1 && mi->nmidas == 1) {
	mi->method = MDS_GSS;
    }

    mi->u = gretl_column_vector_alloc(mi->nobs);
    mi->y = gretl_column_vector_alloc(mi->nobs);
    mi->xi = gretl_column_vector_alloc(mi->nobs);

    if (mi->u == NULL || mi->y == NULL || mi->xi == NULL) {
	return E_ALLOC;
    }

    if (mi->nalmonp == 0) {
	memcpy(mi->y->val, dset->Z[mi->yno] + dset->t1, mi->colsize);
    }

    /* number of low-frequency regressors */
    nb = mi->nx;

    /* How many boxed coeffs do we have? And how many
       coeffs in conditional OLS (@nb)
    */
    for (i=0; i<mi->nmidas; i++) {
	mt = &mi->mterms[i];
	if (beta_type(mt->type)) {
	    bound_rows += 2;
	} else if (mt->type == MIDAS_NEALMON) {
	    bound_rows += mt->nparm;
	}
	if (mt->type == MIDAS_U) {
	    nb += mt->nlags;
	} else if (takes_coeff(mt->type)) {
	    nb += 1;
	}
    }

    mi->bounds = gretl_matrix_alloc(bound_rows, 3);
    if (mi->bounds == NULL) {
	return E_ALLOC;
    }

    /* build X matrix for cond. OLS */
    mi->X = gretl_matrix_alloc(mi->nobs, nb);
    if (mi->X == NULL) {
	return E_ALLOC;
    }

    targ = mi->X->val;
    for (i=0; i<mi->nx; i++) {
	vi = mi->list[i+2];
	src = dset->Z[vi] + dset->t1;
	memcpy(targ, src, mi->colsize);
	targ += mi->nobs;
    }

    k = j = 0;

    for (i=0; i<mi->nmidas; i++) {
	mt = &mi->mterms[i];
	if (mi->bounds != NULL) {
	    if (beta_type(mt->type)) {
		/* set up beta minima and maxima */
		for (ii=0; ii<2; ii++) {
		    /* columns: index, minimum, maximum */
		    gretl_matrix_set(mi->bounds, j, 0, k + ii + 1);
		    if (mt->type == MIDAS_BETA1) {
			if (ii == 0) {
			    gretl_matrix_set(mi->bounds, j, 1, mt->theta->val[0]);
			    gretl_matrix_set(mi->bounds, j, 2, mt->theta->val[0]);
			} else {
			    gretl_matrix_set(mi->bounds, j, 1, 1.0 + eps);
			    gretl_matrix_set(mi->bounds, j, 2, 30);
			}
		    } else {
			gretl_matrix_set(mi->bounds, j, 1, eps);
			gretl_matrix_set(mi->bounds, j, 2, 500);
		    }
		    j++;
		}
	    } else if (mt->type == MIDAS_NEALMON) {
		/* set up nealmon minima and maxima */
		for (ii=0; ii<mt->nparm; ii++) {
		    gretl_matrix_set(mi->bounds, j, 0, k + ii + 1);
		    gretl_matrix_set(mi->bounds, j, 1, -2.0);
		    gretl_matrix_set(mi->bounds, j, 2, +2.0);
		    j++;
		}
	    }
	}
	if (mt->type == MIDAS_U) {
	    /* transcribe HF lags data */
	    for (ii=0; ii<mt->nlags; ii++) {
		vi = mt->laglist[ii+1];
		src = dset->Z[vi] + dset->t1;
		memcpy(targ, src, mi->colsize);
		targ += mi->nobs;
	    }
	} else {
	    k += mt->nparm;
	    if (takes_coeff(mt->type)) {
		/* leave space for weighted combo */
		targ += mi->nobs;
	    }
	}
    }

    mi->g = gretl_column_vector_alloc(k);
    mi->b = gretl_column_vector_alloc(nb);
    mi->theta = gretl_column_vector_alloc(k);

    if (mi->g == NULL || mi->b == NULL || mi->theta == NULL) {
	return E_ALLOC;
    }

    type0 = mi->mterms[0].type;

    if (mi->nmidas == 1 && (type0 == MIDAS_BETA0 || type0 == MIDAS_BETA1)) {
	/* FIXME multiple beta terms? */
	midas_beta_init(mi);
    } else {
	/* initialize overall mi->theta by composition
	   of one or more individual theta vectors
	*/
	ii = 0;
	for (i=0; i<mi->nmidas; i++) {
	    mt = &mi->mterms[i];
	    if (mt->type != MIDAS_U) {
		for (j=0; j<mt->nparm; j++) {
		    mi->theta->val[ii++] = mt->theta->val[j];
		}
	    }
	}
    }

    return 0;
}

static void midas_info_destroy (midas_info *mi)
{
    gretl_matrix_free(mi->b);
    gretl_matrix_free(mi->g);
    gretl_matrix_free(mi->theta);
    gretl_matrix_free(mi->bounds);
    gretl_matrix_free(mi->u);
    gretl_matrix_free(mi->X);
    gretl_matrix_free(mi->y);
    gretl_matrix_free(mi->xi);

    free(mi->mterms);
    free(mi->xlist);
    free(mi->seplist);
    g_free(mi->pnames);

    free(mi);
}

/* combined callback that returns SSR given @theta, and also
   calculates the gradient in @g
*/

static double cond_ols_callback (double *theta, double *g,
				 int n, void *ptr)
{
    midas_info *mi = ptr;
    DATASET *dset = mi->dset;
    gretl_matrix *w;
    gretl_matrix *mg;
    int i, j, k, t, s;
    int vj, xcol;
    int err = 0;

#if MIDAS_DEBUG > 1
    fprintf(stderr, "cond_ols_callback: starting\n");
    for (i=0; i<n; i++) {
	fprintf(stderr, " theta[%d] = %g\n", i, theta[i]);
    }
#endif

    if (mi->nalmonp > 0) {
	/* reset mi->y, it may have been altered */
	memcpy(mi->y->val, dset->Z[mi->yno] + dset->t1, mi->colsize);
    }

    /* Fill relevant columns of mi->X with weighted
       linear combinations of HF data
    */

    /* initialize index of column to fill */
    xcol = mi->nx;

    /* initialize index for reading from @theta */
    k = 0;

    for (i=0; i<mi->nmidas && !err; i++) {
	midas_term *mt = &mi->mterms[i];

	if (mt->type == MIDAS_U) {
	    /* copy X data: job done already */
	    xcol += mt->nlags;
	    continue;
	}
	for (j=0; j<mt->nparm; j++) {
	    mt->theta->val[j] = theta[k++];
	}
	w = midas_weights(mt->nlags, mt->theta, mt->type, &err);
	if (!err && mt->type == MIDAS_ALMONP) {
	    /* net the estimated effect out of @y */
	    for (j=0; j<mt->nlags; j++) {
		vj = mt->laglist[j+1];
		s = 0;
		for (t=dset->t1; t<=dset->t2; t++) {
		    mi->y->val[s++] -= dset->Z[vj][t] * w->val[j];
		}
	    }
	} else if (!err) {
	    /* fill the next X column and advance */
	    gretl_matrix_zero(mi->xi);
	    for (j=0; j<mt->nlags; j++) {
		vj = mt->laglist[j+1];
		s = 0;
		for (t=dset->t1; t<=dset->t2; t++) {
		    mi->xi->val[s++] += dset->Z[vj][t] * w->val[j];
		}
	    }
	    for (t=0; t<mi->nobs; t++) {
		gretl_matrix_set(mi->X, t, xcol, mi->xi->val[t]);
	    }
	    xcol++;
	}
	gretl_matrix_free(w);
    }

    if (!err) {
	/* run OLS conditional on theta */
	err = gretl_matrix_ols(mi->y, mi->X, mi->b, NULL, mi->u, NULL);
    }

    if (g == NULL) {
	goto skip_gradient;
    }

    /* now for the gradient with respect to the hyper-parameters
       in @theta
    */

    k = 0;
    xcol = mi->nx;

    for (i=0; i<mi->nmidas && !err; i++) {
	midas_term *mt = &mi->mterms[i];
	double xit;
	int p, vp;

	if (mt->type == MIDAS_U) {
	    xcol += mt->nlags;
	    continue;
	}
	mg = midas_gradient(mt->nlags, mt->theta, mt->type, &err);
	if (err) {
	    break;
	}
	for (j=0; j<mg->cols; j++) {
	    /* loop across hyper-parameters for this term */
	    gretl_matrix_zero(mi->xi);
	    for (t=0; t<mi->nobs; t++) {
		/* loop across observations */
		s = t + dset->t1;
		xit = 0.0;
		for (p=0; p<mt->nlags; p++) {
		    /* loop across HF lags at this observation */
		    vp = mt->laglist[p+1];
		    xit += dset->Z[vp][s] * gretl_matrix_get(mg, p, j);
		}
		mi->xi->val[t] = xit;
	    }
	    if (takes_coeff(mt->type)) {
		gretl_matrix_multiply_by_scalar(mi->xi, mi->b->val[xcol]);
	    }
	    /* compute gradient of SSR wrt hyper-parameter @k */
	    g[k] = 0.0;
	    for (t=0; t<mi->nobs; t++) {
		g[k] -= 2 * mi->xi->val[t] * mi->u->val[t];
	    }
	    /* advance to next hyper-parameter */
	    k++;
	}
	if (takes_coeff(mt->type)) {
	    /* advance read position for OLS coeff */
	    xcol++;
	}
	gretl_matrix_free(mg);
    }

 skip_gradient:

    if (err) {
	const char *s = errmsg_get_with_default(err);

	fprintf(stderr, "cond_ols_callback: err=%d (%s)\n", err, s);
	return NADBL;
    }

    mi->SSR = 0.0;
    for (t=0; t<mi->nobs; t++) {
	mi->SSR += mi->u->val[t] * mi->u->val[t];
    }

    return mi->SSR;
}

/* Given mi->b (OLS coefficients) and mi->theta (hyper-
   parameters optimized by L-BFGS-B), combine them into
   the full coefficient array.
*/

static int make_full_coeff_vector (midas_info *mi, int n)
{
    gretl_matrix *c;
    midas_term *mt;
    int i, j, k, l, m;

    c = gretl_column_vector_alloc(n);

    if (c == NULL) {
	return E_ALLOC;
    }

    k = l = m = 0;

    for (i=0; i<mi->nx; i++) {
	c->val[k++] = mi->b->val[l++];
    }
    for (i=0; i<mi->nmidas; i++) {
	mt = &mi->mterms[i];
	if (mt->type == MIDAS_U) {
	    for (j=0; j<mt->nlags; j++) {
		c->val[k++] = mi->b->val[l++];
	    }
	} else {
	    if (takes_coeff(mt->type)) {
		c->val[k++] = mi->b->val[l++];
	    }
	    for (j=0; j<mt->nparm; j++) {
		c->val[k++] = mi->theta->val[m++];
	    }
	}
    }

    /* swap out the original mi->b */
    gretl_matrix_free(mi->b);
    mi->b = c;

    return 0;
}

static int transcribe_to_nlspec (nlspec *s,
				 midas_info *mi,
				 gretlopt opt)
{
    int err = 0;

    s->ncoeff = mi->b->rows + mi->theta->rows;
    err = make_full_coeff_vector(mi, s->ncoeff);
    if (err) {
	return err;
    }

    s->ci = MIDASREG;
    s->opt = opt;
    s->dv = mi->list[1];
    s->coeff = mi->b->val;
    s->t1 = mi->dset->t1;
    s->t2 = mi->dset->t2;
    s->crit = mi->SSR;
    s->parnames = mi->pnames;
    s->fvec = mi->u->val;
    s->dset = mi->dset;

    return err;
}

/* The Gauss-Newton Regression after running BFGS plus
   conditional OLS */

static int cond_ols_GNR (MODEL *pmod,
			 midas_info *mi,
			 gretlopt opt,
			 PRN *prn)
{
    DATASET *dset = mi->dset;
    DATASET *gdset = NULL;
    int *glist = NULL;
    gretl_matrix *w, *G;
    int nc, zcol, vi;
    int i, j, k, kk, t, s, v;
    int cpi, *cpos = NULL;
    int err = 0;

#if MIDAS_DEBUG
    fprintf(stderr, "cond_ols_GNR, starting\n");
#endif

    /* total number of coefficients in the GNR */
    nc = mi->b->rows + mi->theta->rows;

    v = nc + 2;
    for (i=2; i<=mi->list[0]; i++) {
	if (mi->list[i] == 0) {
	    v--;
	    break;
	}
    }

    glist = gretl_list_new(nc + 1);

    gdset = create_auxiliary_dataset(v, mi->nobs, 0);
    if (gdset == NULL) {
	return E_ALLOC;
    }

    if (dataset_is_time_series(dset)) {
	gdset->structure = dset->structure;
	gdset->pd = dset->pd;
	ntolabel(gdset->stobs, dset->t1, dset);
	gdset->sd0 = get_date_x(gdset->pd, gdset->stobs);
    }

    /* dependent var: nls residual */
    glist[0] = glist[1] = 1;
    memcpy(gdset->Z[1], mi->u->val, mi->colsize);
    strcpy(gdset->varname[1], "gy");

    zcol = 2;

    /* low-frequency regressors */
    for (i=2; i<=mi->list[0]; i++) {
	int gvi;

	vi = mi->list[i];
	if (vi == 0) {
	    gvi = 0;
	} else {
	    gvi = zcol++;
	    memcpy(gdset->Z[gvi], dset->Z[vi] + dset->t1, mi->colsize);
	    sprintf(gdset->varname[gvi], "gx%d", i-1);
	}
	glist[i] = gvi;
	glist[0] += 1;
    }

    if (mi->nbeta1 > 0) {
	cpos = gretl_list_new(mi->nbeta1);
	if (cpos == NULL) {
	    err = E_ALLOC;
	}
    }

    /* the read position for OLS coefficients */
    k = mi->nx;

    /* the write position and coeff position for beta1 terms */
    cpi = 1;
    kk = mi->nx;

    /* MIDAS terms */
    for (i=0; i<mi->nmidas && !err; i++) {
	midas_term *mt = &mi->mterms[i];
	double zt, gij, hfb = 1.0;
	int ii, pos;

	if (mt->type == MIDAS_U) {
	    /* gradient wrt U-MIDAS coeffs */
	    for (j=0; j<mt->nlags; j++) {
		vi = mt->laglist[j+1];
		memcpy(gdset->Z[zcol], dset->Z[vi] + dset->t1,
		       mi->colsize);
		glist[0] += 1;
		glist[glist[0]] = zcol;
		sprintf(gdset->varname[zcol], "mdx%d", i+1);
		zcol++;
	    }
	    k += mt->nlags;
	    kk += mt->nlags;
	    continue;
	}

	/* gradient wrt weighted linear combination */
	w = midas_weights(mt->nlags, mt->theta, mt->type, &err);
	if (!err) {
	    s = 0;
	    for (t=dset->t1; t<=dset->t2; t++) {
		zt = 0.0;
		for (j=0; j<mt->nlags; j++) {
		    vi = mt->laglist[j+1];
		    zt += dset->Z[vi][t] * w->val[j];
		}
		gdset->Z[zcol][s] = zt;
		s++;
	    }
	    gretl_matrix_free(w);
	}
	glist[0] += 1;
	glist[glist[0]] = zcol;
	sprintf(gdset->varname[zcol], "mdx%d", i+1);
	zcol++;
	kk++;

	/* gradient wrt hyper-parameters */
	pos = zcol;
	if (takes_coeff(mt->type)) {
	    hfb = mi->b->val[k++];
	}
	G = midas_gradient(mt->nlags, mt->theta, mt->type, &err);
	if (!err) {
	    s = 0;
	    for (t=dset->t1; t<=dset->t2; t++) {
		for (ii=0; ii<mt->nparm; ii++) {
		    zt = 0.0;
		    for (j=0; j<mt->nlags; j++) {
			vi = mt->laglist[j+1];
			gij = gretl_matrix_get(G, j, ii);
			zt += dset->Z[vi][t] * gij;
		    }
		    gdset->Z[pos+ii][s] = zt * hfb;
		}
		s++;
	    }
	    gretl_matrix_free(G);
	}

	if (!err) {
	    for (ii=0; ii<mt->nparm; ii++) {
		glist[0] += 1;
		glist[glist[0]] = pos + ii;
		if (ii == 0 && mt->type == MIDAS_BETA1) {
		    cpos[cpi++] = kk;
		    sprintf(gdset->varname[pos+ii], "grad%dc", ii+1);
		} else {
		    sprintf(gdset->varname[pos+ii], "grad%d", ii+1);
		}
	    }
	    zcol += mt->nparm;
	    kk += mt->nparm;
	}
    }

#if MIDAS_DEBUG
    printlist(glist, "glist, mi");
    printlist(cpos, "cpos");
# if MIDAS_DEBUG > 1
    printdata(glist, NULL, gdset, OPT_O, prn);
# endif
#endif

    *pmod = GNR(glist, gdset, opt, prn);

#if MIDAS_DEBUG
    printmodel(pmod, gdset, opt, prn);
#endif

    if (pmod->errcode == 0 || pmod->errcode == E_JACOBIAN) {
	nlspec spec = {0};

	if (pmod->errcode == 0 && (opt & OPT_B)) {
	    /* test for structural break */
	    QLR_test(pmod, gdset, OPT_Q | OPT_M, NULL);
	}
	err = transcribe_to_nlspec(&spec, mi, opt);
	if (!err) {
	    err = finalize_nls_model(pmod, &spec, 0, glist);
	    if (!err && cpos != NULL) {
		/* invalidate std errors of clamped coeffs */
		for (cpi=1; cpi<=cpos[0]; cpi++) {
		    i = cpos[cpi];
		    if (i < pmod->ncoeff) {
			/* this conditionality should NOT be needed! */
			pmod->sderr[i] = NADBL;
		    }
		}
	    }
	}
	if (err && !pmod->errcode) {
	    pmod->errcode = err;
	}
	if (!pmod->errcode) {
	    set_model_id(pmod, opt);
	}
    }

    destroy_dataset(gdset);
    free(glist);
    free(cpos);

    return err;
}

/* Golden Section search for one-dimensional optimization:
   faster (though somewhat less robust?) than using L_BFGS-B.
*/

static int midas_gss (double *theta, int n, midas_info *mi,
		      int *fncount)
{
    double gr = (sqrt(5.0) + 1) / 2.0;
    double a = gretl_matrix_get(mi->bounds, 1, 1);
    double b = gretl_matrix_get(mi->bounds, 1, 2);
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;
    double tol = 1.0e-4;
    double fc, fd;
    int err = 0;

#if MIDAS_DEBUG
    fprintf(stderr, "gss: bounds = [%g,%g], incoming value %g\n",
	    a, b, theta[1]);
#endif

    /* theta[0] will be the clamped term, and
       theta[1] the one we're optimizing
    */

    while (fabs(c - d) > tol) {
	theta[1] = c;
	fc = cond_ols_callback(theta, NULL, n, mi);
	theta[1] = d;
	fd = cond_ols_callback(theta, NULL, n, mi);
	*fncount += 2;
	if (na(fc) || na(fd)) {
	    err = E_NAN;
	    break;
	}
        if (fc < fd) {
            b = d;
        } else {
            a = c;
	}
	/* recompute c and d to preserve precision */
	c = b - (b - a) / gr;
	d = a + (b - a) / gr;
    }

    if (!err) {
	/* record the optimum */
	theta[1] = (b + a) / 2;
    }

    return err;
}

/* driver function for estimating MIDAS model using
   L-BFGS-B or Golden Section search */

static int midas_bfgs_run (MODEL *pmod, midas_info *mi,
			   gretlopt opt, PRN *prn)
{
    int n = gretl_vector_get_length(mi->theta);
    double *theta = mi->theta->val;
    int fncount = 0, grcount = 0;
    int err;

    if (theta[0] == 1.0 && mi->nmidas == 1 && mi->nbeta1 == 1) {
#if MIDAS_DEBUG
	fprintf(stderr, "midas_bfgs_run: calling midas_gss\n");
#endif
	err = midas_gss(theta, n, mi, &fncount);
    } else {
	double reltol = libset_get_double(BFGS_TOLER);

#if MIDAS_DEBUG
	fprintf(stderr, "midas_bfgs_run: calling LBFGS_max\n");
#endif
	err = LBFGS_max(theta, n, 1000, reltol,
			&fncount, &grcount, NULL,
			C_SSR, NULL, cond_ols_callback,
			mi, mi->bounds, opt, prn);
    }

    if (!err) {
	err = cond_ols_GNR(pmod, mi, opt, prn);
    }

    if (!err) {
	gretl_model_set_int(pmod, "iters", fncount);
    }

    return err;
}

/* Below: we add a column vector containing the "gross"
   MIDAS coefficients (i.e. high-frequency slope, where
   applicable, times the weights implied by the hyper-
   parameter values). This is required for forecasting.
   In addition, if there's just a single MIDAS term, we
   add to the model a two-column matrix holding the
   gross coefficients and the associated lags -- this is
   to support a plot of the coefficients in the GUI.
*/

static int add_gross_midas_coeffs (MODEL *pmod,
				   midas_info *mi)
{
    const double *b = pmod->coeff + mi->nx;
    gretl_matrix *hfb, *C = NULL;
    int j, k = 0, nt = 0;
    int row = 0;
    int err = 0;

    for (j=0; j<mi->nmidas; j++) {
	nt += mi->mterms[j].nlags;
    }

    /* vector to hold all gross midas coeffs */
    hfb = gretl_column_vector_alloc(nt);
    if (hfb == NULL) {
	err = E_ALLOC;
    }

    if (!err && mi->nmidas == 1) {
	/* matrix for plotting purposes */
	C = gretl_matrix_alloc(mi->mterms[0].nlags, 2);
    }

    for (j=0; j<mi->nmidas && !err; j++) {
	midas_term *mt = &mi->mterms[j];
	gretl_matrix *theta = NULL;

	if (mt->type != MIDAS_U) {
	    theta = gretl_matrix_alloc(mt->nparm, 1);
	    if (theta == NULL) {
		err = E_ALLOC;
	    }
	}

	if (!err) {
	    gretl_matrix *w = NULL;
	    double ci, hfslope = 0;
	    int p = mt->minlag;
	    int i;

	    if (mt->type == MIDAS_U) {
		for (i=0; i<mt->nparm; i++) {
		    gretl_vector_set(hfb, row++, b[k]);
		    if (C != NULL) {
			gretl_matrix_set(C, i, 0, b[k]);
			gretl_matrix_set(C, i, 1, p++);
		    }
		    k++;
		}
	    } else {
		/* compute slope times weights */
		hfslope = takes_coeff(mt->type) ? b[k++] : 1.0;
		for (i=0; i<mt->nparm; i++) {
		    theta->val[i] = b[k++];
		}
		w = midas_weights(mt->nlags, theta, mt->type, &err);
		if (!err) {
		    for (i=0; i<mt->nlags; i++) {
			ci = hfslope * w->val[i];
			gretl_vector_set(hfb, row++, ci);
			if (C != NULL) {
			    gretl_matrix_set(C, i, 0, ci);
			    gretl_matrix_set(C, i, 1, p++);
			}
		    }
		}
		gretl_matrix_free(w);
	    }
	}

	if (theta != NULL) {
	    gretl_matrix_free(theta);
	    theta = NULL;
	}

	if (!err && C != NULL) {
	    /* finalize plotting matrix and attach to model */
	    char *cnames[2] = {"coeff", "lag"};
	    char **S = strings_array_dup(cnames, 2);

	    if (S != NULL) {
		gretl_matrix_set_colnames(C, S);
	    }
	    gretl_model_set_matrix_as_data(pmod, "midas_coeffs", C);
	} else {
	    gretl_matrix_free(C);
	}
    }

    if (!err) {
	gretl_model_set_matrix_as_data(pmod, "hfb", hfb);
    }

    return err;
}

/* Add to @pmod an array of bundles holding information
   on the individual MIDAS terms (although in most cases
   there will be just one such bundle).
*/

static int model_add_mterms_array (MODEL *pmod,
				   midas_term *mterms,
				   int nmidas)
{
    gretl_array *A;
    int err = 0;

    A = gretl_array_new(GRETL_TYPE_BUNDLES, nmidas, &err);

    if (A != NULL) {
	midas_term *mt;
	gretl_bundle *b;
	int i;

	for (i=0; i<nmidas && !err; i++) {
	    b = gretl_bundle_new();
	    if (b == NULL) {
		err = E_ALLOC;
	    } else {
		mt = &mterms[i];
		err += gretl_bundle_set_string(b, "lname",  mt->lnam0);
		err += gretl_bundle_set_int(b, "prelag", prelag(mt) ? 1 : 0);
		err += gretl_bundle_set_int(b, "minlag", mt->minlag);
		err += gretl_bundle_set_int(b, "maxlag", mt->maxlag);
		err += gretl_bundle_set_int(b, "type",  mt->type);
		err += gretl_bundle_set_int(b, "nparm", mt->nparm);
		err += gretl_array_set_bundle(A, i, b, 0);
		if (err) {
		    err = E_DATA;
		}
	    }
	}

	if (err) {
	    gretl_array_destroy(A);
	} else {
	    gretl_model_set_array_as_data(pmod, "midas_info", A);
	}
    }

    return err;
}

/* Retrieve from array @A, which was earlier attached to
   a model, information on the individual MIDAS terms
   for the purpose of forecasting.
*/

static midas_term *mterms_from_array (gretl_array *A,
				      int *nmidas,
				      int *err)
{
    int n = gretl_array_get_length(A);
    midas_term *mt, *mterms = NULL;

    if (n == 0) {
	*err = E_DATA;
    } else {
	mterms = malloc(n * sizeof *mterms);
	if (mterms == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (mterms != NULL) {
	gretl_bundle *b;
	int i;

	for (i=0; i<n && !*err; i++) {
	    mt = &mterms[i];
	    midas_term_init(mt);
	    b = gretl_array_get_bundle(A, i);
	    if (b == NULL) {
		*err = E_DATA;
	    } else {
		strcpy(mt->lname, gretl_bundle_get_string(b, "lname", err));
		if (gretl_bundle_get_int(b, "prelag", err)) {
		    mt->flags |= M_PRELAG;
		}
		mt->minlag = gretl_bundle_get_int(b, "minlag", err);
		mt->maxlag = gretl_bundle_get_int(b, "maxlag", err);
		mt->type   = gretl_bundle_get_int(b, "type", err);
		mt->nparm  = gretl_bundle_get_int(b, "nparm", err);
		sprintf(mt->mname, "theta___%d", i+1);
	    }
	}

	if (*err) {
	    free(mterms);
	    mterms = NULL;
	}
    }

    if (!*err) {
	*nmidas = n;
    }

    return mterms;
}

/* For use when printing a midasreg model: returns a line
   to be printed before the coefficients associated with a
   MIDAS term, or NULL on failure.
*/

const char *get_midas_term_line (const MODEL *pmod, int i)
{
    static char targ[MAXLEN];
    const char *ret = NULL;
    gretl_array *A = gretl_model_get_data(pmod, "midas_info");
    gretl_bundle *b = NULL;

    *targ = '\0';

    if (A != NULL) {
	int n = gretl_array_get_length(A);

	if (i >= 0 && i < n) {
	    b = gretl_array_get_bundle(A, i);
	}
    }

    if (b != NULL) {
	const char *lname;
	int l1, l2, err = 0;

	lname = gretl_bundle_get_string(b, "lname", &err);
	l1 = gretl_bundle_get_int(b, "minlag", &err);
	l2 = gretl_bundle_get_int(b, "maxlag", &err);

	if (!err) {
	    sprintf(targ, _("MIDAS list %s, high-frequency lags %d to %d"),
		    lname, l1, l2);
	    ret = targ;
	}
    }

    return ret;
}

/* Get the MIDAS model ready for shipping out. What
   exactly we do here depends in part on whether
   estimation was done by NLS or OLS.
*/

static int finalize_midas_model (MODEL *pmod,
				 midas_info *mi,
				 const char *param,
				 gretlopt opt,
				 const DATASET *dset)
{
    int type0, mixed = 0;
    int i, err = 0;

    /* ensure that the incoming sample range is recorded */
    gretl_model_smpl_init(pmod, dset);

    gretl_model_set_string_as_data(pmod, "midas_spec",
				   gretl_strdup(param));

    if (pmod->ci == OLS) {
	/* @pmod is the result of OLS estimation */
	int vi;

	gretl_model_allocate_param_names(pmod, pmod->ncoeff);
	for (i=0; i<pmod->ncoeff; i++) {
	    vi = pmod->list[i+2];
	    gretl_model_set_param_name(pmod, i, dset->varname[vi]);
	}
	gretl_model_set_int(pmod, "umidas", 1);
    } else {
	/* @pmod is the result of some sort of NLS estimation */
	free(pmod->depvar);
	pmod->depvar = gretl_strdup(dset->varname[mi->list[1]]);
	free(pmod->list);
	pmod->list = gretl_list_copy(mi->list);
	if (mi->method == MDS_BFGS) {
	    gretl_model_set_int(pmod, "BFGS", 1);
	} else if (mi->method == MDS_GSS) {
	    gretl_model_set_int(pmod, "GSS", 1);
	}
    }

    pmod->ci = MIDASREG;

    /* record list of low-frequency regressors */
    if (mi->xlist == NULL) {
	gretl_model_set_int(pmod, "no_lfx", 1);
    } else {
	gretl_model_set_list_as_data(pmod, "lfxlist", mi->xlist);
	mi->xlist = NULL; /* donated to pmod */
    }

    /* and positions where MIDAS terms start */
    if (mi->seplist != NULL) {
	gretl_model_set_list_as_data(pmod, "seplist", mi->seplist);
	mi->seplist = NULL; /* donated to pmod */
    }

    if (mi->ldepvar) {
	gretl_model_set_int(pmod, "ldepvar", mi->ldepvar);
	gretl_model_set_int(pmod, "dynamic", 1);
    }

    /* record the (common) type of MIDAS term? */
    type0 = mi->mterms[0].type;
    if (mi->nmidas > 1) {
	for (i=1; i<mi->nmidas; i++) {
	    if (mi->mterms[i].type != type0) {
		mixed = 1;
		break;
	    }
	}
    }
    if (!mixed && type0 > 0) {
	gretl_model_set_int(pmod, "midas_type", type0);
    }

    if (!err) {
	/* Compute and record the "gross" MIDAS coefficients:
	   we may need these for forecasting, and possibly to
	   enable a plot in the GUI.
	*/
	err = add_gross_midas_coeffs(pmod, mi);
    }

    if (!err) {
	err = model_add_mterms_array(pmod, mi->mterms, mi->nmidas);
    }

    return err;
}

/* simple initialization of MIDAS slope terms for use with
   Levenberg-Marquardt */

static int midas_means_init (midas_info *mi,
			     const gretl_matrix *y,
			     const gretl_matrix *X,
			     gretl_matrix *c,
			     const DATASET *dset)
{
    gretl_matrix *XZ;
    int err = 0;

    XZ = gretl_matrix_alloc(mi->nobs, mi->nx + mi->hfslopes);

    if (XZ == NULL) {
	err = E_ALLOC;
    } else {
	midas_term *mt;
	double zti;
	int i, j, k, s, t;

	if (X != NULL) {
	    /* transcribe low-frequency regressors */
	    memcpy(XZ->val, X->val, mi->nobs * mi->nx * sizeof(double));
	}

	/* transcribe time-mean of MIDAS terms */
	k = mi->nx * mi->nobs;
	for (i=0; i<mi->nmidas; i++) {
	    mt = &mi->mterms[i];
	    if (takes_coeff(mt->type)) {
		for (t=0; t<mi->nobs; t++) {
		    s = t + dset->t1;
		    zti = 0.0;
		    for (j=0; j<mt->nlags; j++) {
			zti += dset->Z[mt->laglist[j+1]][s];
		    }
		    XZ->val[k++] = zti / mt->nlags;
		}
	    }
	}

	err = gretl_matrix_ols(y, XZ, c, NULL,
			       NULL, NULL);
	gretl_matrix_free(XZ);
    }

    return err;
}

/* Define "private" matrices to hold the regular X data
   (MX___) and the vector of coefficients on these data
   (bx___). If we have an "hfslope" terms, then also try
   some initialization.

   Note: we come here only if we're using native NLS
   (not U-MIDAS OLS, or BFGS + conditional OLS).
*/

static int add_midas_matrices (midas_info *mi,
			       const DATASET *dset)
{
    gretl_matrix *X = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *c = NULL;
    int init_err = 0;
    int err = 0;

    if (mi->nx > 0) {
	X = gretl_matrix_data_subset(mi->xlist, dset,
				     dset->t1, dset->t2,
				     M_MISSING_ERROR,
				     &err);
	if (!err) {
	    err = private_matrix_add(X, "MX___");
	}

	if (!err) {
	    b = gretl_zero_matrix_new(mi->nx, 1);
	    if (b != NULL) {
		err = private_matrix_add(b, "bx___");
	    } else {
		err = E_ALLOC;
	    }
	}
    }

    if (mi->nx == 0 && mi->hfslopes == 0) {
	return err;
    }

    if (!err) {
	/* for initialization only */
	y = gretl_column_vector_alloc(mi->nobs);
	if (y != NULL) {
	    memcpy(y->val, dset->Z[mi->yno] + dset->t1,
		   mi->colsize);
	} else {
	    init_err = 1;
	}
    }

    if (!err && !init_err) {
	/* "full-length" coeff vector */
	c = gretl_zero_matrix_new(mi->nx + mi->hfslopes, 1);
	if (c == NULL) {
	    init_err = 1;
	}
    }

    if (!err && !init_err) {
	if (mi->hfslopes > 0) {
	    init_err = midas_means_init(mi, y, X, c, dset);
	}
	if (init_err || (mi->hfslopes == 0 && mi->nx > 0)) {
	    /* simpler initialization, ignoring HF terms */
	    c->rows = mi->nx;
	    init_err = gretl_matrix_ols(y, X, c, NULL,
					NULL, NULL);
	}
    }

    if (!err) {
	int i;

	if (!init_err) {
	    /* initialize X coeffs from OLS */
	    for (i=0; i<mi->nx; i++) {
		b->val[i] = c->val[i];
	    }
	}
	if (mi->hfslopes > 0) {
	    /* initialize HF slopes, with fallback to zero */
	    int use_c = !init_err && c->rows > mi->nx;
	    char tmp[24];
	    double bzi;

	    for (i=0; i<mi->nmidas && !err; i++) {
		if (takes_coeff(mi->mterms[i].type)) {
		    sprintf(tmp, "bmlc___%d", i + 1);
		    bzi = use_c ? c->val[mi->nx+i] : 0.0;
		    err = private_scalar_add(bzi, tmp);
		}
	    }
	}
    }

    gretl_matrix_free(y);
    gretl_matrix_free(c);

    return err;
}

static int put_midas_nls_line (char *line,
			       DATASET *dset,
			       gretlopt opt,
			       PRN *prn)
{
    int err;

    if (opt & OPT_P) {
	/* display what we're passing to nls */
	pputs(prn, line);
	pputc(prn, '\n');
    }

    err = nl_parse_line(NLS, line, dset, prn);
    *line = '\0';

    return err;
}

/* Append @pname to @s, ensuring that it's preceded by
   a single space if it's not preceded by a double quote.
*/

static void append_pname (char *s, const char *pname)
{
    if (*s != '\0') {
	char c = s[strlen(s) - 1];

	if (c != ' ' && c != '"') {
	    strcat(s, " ");
	}
    }

    strcat(s, pname);
}

static void make_pname (char *targ, midas_term *mt, int i,
			const DATASET *dset)
{
    if (mt->type == MIDAS_NEALMON) {
	sprintf(targ, "Almon%d", i+1);
    } else if (beta_type(mt->type)) {
	sprintf(targ, "Beta%d", i+1);
    } else if (mt->type == MIDAS_ALMONP) {
	sprintf(targ, "Almon%d", i);
    } else {
	/* U-MIDAS */
	int *list = get_list_by_name(mt->lname);

	if (list != NULL) {
	    strcpy(targ, dset->varname[list[i+1]]);
	} else {
	    sprintf(targ, "U-MIDAS%d", i+1);
	}
    }
}

static int add_param_names (midas_info *mi,
			    const DATASET *dset)
{
    char tmp[64], str[2*MAXLEN];
    int i, j, k = 0;
    int err = 0;

    mi->seplist = gretl_list_new(mi->nmidas);
    *tmp = *str = '\0';

    if (mi->xlist != NULL) {
	for (i=1; i<=mi->xlist[0]; i++) {
	    strcpy(tmp, dset->varname[mi->xlist[i]]);
	    append_pname(str, tmp);
	}
	k += mi->xlist[0];
    }

    for (i=0; i<mi->nmidas; i++) {
	if (mi->seplist != NULL) {
	    mi->seplist[i+1] = k;
	}
	if (takes_coeff(mi->mterms[i].type)) {
	    if (mi->hfslopes > 1) {
		sprintf(tmp, "HF_slope%d", i+1);
	    } else {
		strcpy(tmp, "HF_slope");
	    }
	    append_pname(str, tmp);
	    k++;
	}
	for (j=0; j<mi->mterms[i].nparm; j++) {
	    make_pname(tmp, &mi->mterms[i], j, dset);
	    append_pname(str, tmp);
	    k++;
	}
    }

    mi->pnames = g_strdup(str);

    return err;
}

/* Allocate and initialize midas_info struct, and set
   MDS_NLS (Levenberg-Marquardt) as the default estimation
   method -- but this will be subject to revision when we
   digest the user's specification.
*/

static midas_info *midas_info_new (const int *list,
				   DATASET *dset)
{
    midas_info *mi = malloc(sizeof *mi);

    if (mi != NULL) {
	/* zero everything first */
	memset(mi, 0, sizeof *mi);
	mi->yno = list[1];
	mi->list = list;
	mi->dset = dset;
	mi->method = MDS_NLS;
	mi->orig_toler = 0;
    }

    return mi;
}

/* When using Levenberg-Marquardt NLS, check to see if
   we have any terms that call for a "small" step in
   the optimizer (to try to head off numerical
   problems). This applies to the "highly nonlinear"
   parameterizations.
*/

static int any_smallstep_terms (midas_info *mi)
{
    midas_term *mt;
    int i;

    for (i=0; i<mi->nmidas; i++) {
	mt = &mi->mterms[i];
	if (mt->type == MIDAS_NEALMON ||
	    mt->type == MIDAS_BETA0 ||
	    mt->type == MIDAS_BETAN) {
	    return 1;
	}
    }

    return 0;
}

/* General setup for NLS via Levenverg-Marquardt */

static int midas_nls_setup (midas_info *mi, DATASET *dset,
			    gretlopt opt, PRN *prn)
{
    char tmp[64], line[MAXLEN];
    midas_term *mt;
    const int *list = mi->list;
    int nmidas = mi->nmidas;
    int i, j = 0;
    int err = 0;

    if (opt & OPT_P) {
	pputs(prn, "\n=== auto-generated nls specification ===\n");
    }

    /* initial "nls" line */
    sprintf(line, "nls %s = ", dset->varname[list[1]]);
    if (mi->xlist != NULL) {
	strcat(line, "lincomb(XL___, bx___)");
    }
    for (i=0; i<nmidas; i++) {
	mt = &mi->mterms[i];
	if (takes_coeff(mt->type)) {
	    sprintf(tmp, " + bmlc___%d*mlc___%d", ++j, i+1);
	} else {
	    sprintf(tmp, " + mlc___%d", i+1);
	}
	strcat(line, tmp);
    }
    err = put_midas_nls_line(line, dset, opt, prn);

    /* MIDAS series and gradient matrices */
    for (i=0; i<nmidas && !err; i++) {
	mt = &mi->mterms[i];
	if (mt->type == MIDAS_U) {
	    sprintf(line, "series mlc___%d = lincomb(%s, %s)",
		    i+1, mt->lname, mt->mname);
	} else {
	    sprintf(line, "series mlc___%d = mlincomb(%s, %s, %d)",
		    i+1, mt->lname, mt->mname, mt->type);
	}
	err = put_midas_nls_line(line, dset, opt, prn);
	if (!err && mt->type != MIDAS_U) {
	    sprintf(line, "matrix mgr___%d = mgradient(%d, %s, %d)",
		    i+1, mt->nlags, mt->mname, mt->type);
	    err = put_midas_nls_line(line, dset, opt, prn);
	}
    }

    /* derivatives */
    if (!err && mi->xlist != NULL) {
	strcpy(line, "deriv bx___ = MX___");
	err = put_midas_nls_line(line, dset, opt, prn);
    }
    for (i=0; i<nmidas && !err; i++) {
	mt = &mi->mterms[i];
	if (takes_coeff(mt->type)) {
	    sprintf(line, "deriv bmlc___%d = {mlc___%d}", i+1, i+1);
	    err = put_midas_nls_line(line, dset, opt, prn);
	}
	if (!err) {
	    if (takes_coeff(mt->type)) {
		sprintf(line, "deriv %s = bmlc___%d * {%s} * mgr___%d",
			mt->mname, i+1, mt->lname, i+1);
	    } else if (mt->type == MIDAS_ALMONP) {
		sprintf(line, "deriv %s = {%s} * mgr___%d",
			mt->mname, mt->lname, i+1);
	    } else {
		sprintf(line, "deriv %s = {%s}", mt->mname, mt->lname);
	    }
	    err = put_midas_nls_line(line, dset, opt, prn);
	}
    }

    if (!err) {
	/* add parameter names */
	sprintf(line, "param_names \"%s\"", mi->pnames);
	err = put_midas_nls_line(line, dset, opt, prn);
    }

    if (opt & OPT_P) {
	pputs(prn, "=== end nls specification ===\n\n");
    }

    if (!err) {
	if (any_smallstep_terms(mi)) {
	    nl_set_smallstep();
	}
	if (mi->nmidas > 1) {
	    /* allow for a sloppier than usual fit? */
	    mi->orig_toler = libset_get_double(NLS_TOLER);
	    libset_set_double(NLS_TOLER, 1.0e-7);
	}
    }

    return err;
}

/* Main driver function for built-in MIDAS estimation.
   The actual engine used is NLS, LBFGS or OLS.
*/

MODEL midas_model (const int *list,
		   const char *param,
		   DATASET *dset,
		   gretlopt opt,
		   PRN *prn)
{
    midas_info *mi = NULL;
    int origv = dset->v;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    MODEL mod;
    int i, err = 0;

    gretl_model_init(&mod, dset);

    mi = midas_info_new(list, dset);
    if (mi == NULL) {
	mod.errcode = E_ALLOC;
	return mod;
    }

#if MIDAS_DEBUG
    fprintf(stderr, "\nmidas_model, starting...\n");
#endif

    if (param == NULL || *param == '\0') {
	err = E_DATA;
    } else {
	err = parse_midas_specs(mi, param, dset, opt);
    }

    if (!err) {
	/* build list of regular regressors */
	err = make_midas_xlist(mi);
	if (mi->xlist != NULL &&
	    (mi->method == MDS_NLS || mi->method == MDS_OLS)) {
	    err = remember_list(mi->xlist, "XL___", NULL);
	    user_var_privatize_by_name("XL___");
	}
    }

    if (!err) {
	/* build (or just read) MIDAS lag-lists */
	err = make_midas_laglists(mi, dset);
    }

    if (!err && mi->method == MDS_OLS) {
	err = umidas_ols(&mod, mi, dset, opt);
	goto midas_finish;
    }

    if (!err) {
	/* determine usable sample range */
	err = midas_set_sample(mi, dset);
    }

    if (!err && mi->method == MDS_NLS && mi->nx > 0) {
	/* add some required matrices */
	err = add_midas_matrices(mi, dset);
    }

    if (!err) {
	/* construct string with names of parameters */
	err = add_param_names(mi, dset);
    }

    if (!err && mi->method == MDS_NLS) {
	/* estimation using native NLS */
	err = midas_nls_setup(mi, dset, opt, prn);
	if (!err) {
	    /* OPT_S: suppress gradient check */
	    mod = nl_model(dset, (opt | OPT_S | OPT_M), prn);
	}
    } else if (!err) {
	/* estimation using L-BFGS-B or GSS */
	err = midas_bfgs_setup(mi, dset, opt);
	if (!err) {
	    err = midas_bfgs_run(&mod, mi, opt, prn);
	}
    }

 midas_finish:

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    for (i=0; i<mi->nmidas; i++) {
	midas_term *mt = &mi->mterms[i];

	if (!prelag(mt)) {
	    free(mt->laglist);
	    if (mt->flags & M_TMPLIST) {
		user_var_delete_by_name(mt->lname, NULL);
	    }
	}
	if (mi->method == MDS_NLS && mt->type != MIDAS_U) {
	    char tmp[24];

	    sprintf(tmp, "mgr___%d", i+1);
	    user_var_delete_by_name(tmp, NULL);
	}
    }

    if (err && !mod.errcode) {
	mod.errcode = err;
    }

    if (!mod.errcode) {
	finalize_midas_model(&mod, mi, param, opt, dset);
    }

    destroy_private_uvars();

    if (mi != NULL) {
	if (mi->orig_toler != 0) {
	    /* put back what we found on input */
	    libset_set_double(NLS_TOLER, mi->orig_toler);
	}
	midas_info_destroy(mi);
    }

    if (dset->v > origv) {
	/* or maybe not? */
	dataset_drop_last_variables(dset, dset->v - origv);
    }

    return mod;
}
