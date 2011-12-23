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
#include "version.h"

#include <glib.h>

/* Levin, Lin and Chu (Journal of Econometrics, 2002), Table 2:
   correction factors for mean (mu) and standard deviation (s)
   in context of panel unit-root statistic.

   T = sample size

   k = 1: no constant
   k = 2: constant included
   k = 3: constant plus trend

*/

static int get_LLC_corrections (int T, int k, double *mu, double *sigma)
{   
    const double LLCfac[] = {
      /*  T    mu1     s1     mu2     s2     mu3     s3 */
	 25, 0.004, 1.049, -0.554, 0.919, -0.703, 1.003,
	 30, 0.003, 1.035, -0.546, 0.889, -0.674, 0.949,
	 35, 0.002, 1.027, -0.541, 0.867, -0.653, 0.906,
	 40, 0.002, 1.021, -0.537, 0.850, -0.637, 0.871,
	 45, 0.001, 1.017, -0.533, 0.837, -0.624, 0.842,
	 50, 0.001, 1.014, -0.531, 0.826, -0.614, 0.818,
	 60, 0.001, 1.011, -0.527, 0.810, -0.598, 0.780,
	 70, 0.000, 1.008, -0.524, 0.798, -0.587, 0.751,
	 80, 0.000, 1.007, -0.521, 0.789, -0.578, 0.728,
	 90, 0.000, 1.006, -0.520, 0.782, -0.571, 0.710,
	100, 0.000, 1.005, -0.518, 0.776, -0.566, 0.695,
	250, 0.000, 1.001, -0.509, 0.742, -0.533, 0.603,
	  0, 0.000, 1.000, -0.500, 0.707, -0.500, 0.500 /* \infty */
    };
    int c = 0, err = 0;

    if (k > 0 && k < 4) {
	c = 2 * k - 1;
    } else {
	err = E_DATA;
    }

    if (!err) {
	int i, r = 12;

	for (i=0; i<12; i++) {
	    if (T <= LLCfac[7*i]) {
		r = i;
		break;
	    }
	}
	*mu = LLCfac[7*r+c];
	*sigma = LLCfac[7*r+c+1];
    }

    return err;
}

/* detrend \delta y for Levin-Lin-Chu case 3 */

static int LLC_detrend (gretl_matrix *dy)
{
    gretl_matrix *X, *b;
    int t, T = dy->rows;
    int err;

    X = gretl_matrix_alloc(T, 2);
    b = gretl_matrix_alloc(2, 1);

    if (X == NULL || b == NULL) {
	err = E_ALLOC;
    } else {
	for (t=0; t<T; t++) {
	    gretl_matrix_set(X, t, 0, 1.0);
	    gretl_matrix_set(X, t, 1, t+1);
	}
	err = gretl_matrix_ols(dy, X, b, NULL, NULL, NULL);
    }

    if (!err) {
	for (t=0; t<T; t++) {
	    /* replace with detrended values */
	    dy->val[t] -= (b->val[0] + b->val[1] * (t+1));
	}
    }

    gretl_matrix_free(X);
    gretl_matrix_free(b);

    return err;
}

/* We could use gretl_long_run_variance() here, except that it
   would have to be generalized to cover the cases (m == 1),
   where we're _not_ subtracting the mean (on the maintained
   hypothesis that the mean = 0), and (m == 3), where we have to
   subtract a linear trend before computing the variance.
*/

static double LLC_lrvar (gretl_matrix *vdy, int K, int m, int *err)
{
    double w, s21 = 0, s22 = 0;
    double *dy = vdy->val;
    int T = vdy->rows;
    int t, j;

    if (m == 3) {
	/* subtract linear trend */
	*err = LLC_detrend(vdy);
	if (*err) {
	    return NADBL;
	}
    } else if (m == 2) {
	/* subtract the mean */
	double dybar = 0;

	for (t=0; t<T; t++) {
	    dybar += dy[t];
	}
	dybar /= T;
	for (t=0; t<T; t++) {
	    dy[t] -= dybar;
	}
    }	

    for (t=0; t<T; t++) {
	s21 += dy[t] * dy[t];
    }

    for (j=1; j<=K; j++) {
	w = 1.0 - j /((double) K + 1);
	for (t=j; t<T; t++) {
	    s22 += w * dy[t] * dy[t-j];
	}
    }

    return (s21 + 2 * s22) / T;
}

/* In case we got a list of individual-specific ADF order terms,
   check that it makes sense and do some basic accounting.
*/

static int LLC_check_plist (const int *list, int N, int *pmax, int *pmin,
			    double *pbar)
{
    int err = 0;

    if (list == NULL || list[0] == 0) {
	err = E_DATA;
    } else if (list[0] > 1 && list[0] != N) {
	err = E_DATA;
    } else {
	int i;

	*pmax = *pmin = *pbar = 0;

	for (i=1; i<=list[0]; i++) {
	    if (list[i] < 0) {
		err = E_DATA;
		break;
	    } 
	    if (list[i] > *pmax) {
		*pmax = list[i];
	    }
	    if (i == 1 || list[i] < *pmin) {
		*pmin = list[i];
	    }
	    *pbar += list[i];
	}

	if (list[0] > 1) {
	    *pbar /= N;
	}
    }

    return err;
}

static int LLC_sample_check (int N, int t1, int t2, int m,
			     const int *plist, int *NT)
{
    int i, p, minT, T;
    int err = 0;

    *NT = 0;

    for (i=1; i<=plist[0] && !err; i++) {
	p = plist[i];
	minT = m + p + 1; /* ensure df > 0 */

	if (minT < 4) {
	    minT = 4;
	}

	/* T_i denotes the regression-usable series length, after
	   accounting for required lags */
	T = t2 - t1 + 1 - (1 + p);
	if (T < minT) {
	    err = E_DATA;
	} else if (plist[0] == 1) {
	    *NT = N * T;
	} else {
	    *NT += T;
	}
    }

    return err;
}

static const char *DF_test_spec (int m)
{
    const char *tests[] = {
	N_("test without constant"),
	N_("test with constant"),
	N_("with constant and trend"),
    };

    if (m > 0 && m < 4) {
	return tests[m-1];
    } else {
	return "";
    }
}

#define LLC_DEBUG 0

/* Levin-Lin-Chu panel unit-root test */

int real_levin_lin (int vnum, const int *plist, DATASET *dset, 
		    gretlopt opt, PRN *prn)
{
    int u0 = dset->t1 / dset->pd;
    int uN = dset->t2 / dset->pd;
    int N = uN - u0 + 1; /* units in sample range */
    gretl_matrix_block *B;
    gretl_matrix *y, *yavg, *b;
    gretl_matrix *dy, *X, *ui;
    gretl_matrix *e, *ei, *v, *vi;
    gretl_matrix *eps;
    double pbar, SN = 0;
    int t, t1, t2, T, NT;
    int s, pt1, pt2, dyT;
    int i, j, k, K, m;
    int p, pmax, pmin;
    int bigrow, p_varies = 0;
    int err;
    
    err = LLC_check_plist(plist, N, &pmax, &pmin, &pbar);

    if (err) {
	return err;
    }

    /* the 'case' (1 = no const, 2 = const, 3 = const + trend */
    m = 2; /* the default */
    if (opt & OPT_N) {
	/* --nc */
	m = 1;
    } else if (opt & OPT_T) {
	/* --ct */
	m = 3;
    }

    /* does p vary by individual? */
    if (pmax > pmin) {
	p_varies = 1;
    }
    p = pmax;

    /* the max number of regressors */
    k = m + pmax;

    t1 = t2 = 0;
    
    /* check that we have a useable common sample */
    
    for (i=0; i<N && !err; i++) {
	int pt1 = (i + u0) * dset->pd;
	int t1i, t2i;

	dset->t1 = pt1;
	dset->t2 = dset->t1 + dset->pd - 1;
	err = series_adjust_sample(dset->Z[vnum], &dset->t1, &dset->t2);
	t1i = dset->t1 - pt1;
	t2i = dset->t2 - pt1;
	if (i == 0) {
	    t1 = t1i;
	    t2 = t2i;
	} else if (t1i != t1 || t2i != t2) {
	    err = E_MISSDATA;
	    break;
	}
    }

    if (!err) {
	err = LLC_sample_check(N, t1, t2, m, plist, &NT);
    } 

    if (!err) {
	int Tbar = NT / N;

	/* the biggest T we'll need for regressions */
	T = t2 - t1 + 1 - (1 + pmin);

	/* Bartlett lag truncation (Andrews, 1991) */
	K = (int) floor(3.21 * pow(Tbar, 1.0/3));
	if (K > Tbar - 3) {
	    K = Tbar - 3;
	}	

	/* full length of dy vector */
	dyT = t2 - t1;

	B = gretl_matrix_block_new(&y, T, 1,
				   &yavg, T+1+p, 1,
				   &dy, dyT, 1,
				   &X, T, k,
				   &b, k, 1,
				   &ui, T, 1,
				   &ei, T, 1,
				   &vi, T, 1,
				   &e, NT, 1,
				   &v, NT, 1,
				   &eps, NT, 1,
				   NULL);
	if (B == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	return err;
    }

    if (m > 1) {
	/* constant in first column, if wanted */
	for (t=0; t<T; t++) {
	    gretl_matrix_set(X, t, 0, 1.0);
	}
    }

    if (m == 3) {
	/* trend in second column, if wanted */
	for (t=0; t<T; t++) {
	    gretl_matrix_set(X, t, 1, t+1);
	}
    }    

    gretl_matrix_zero(yavg);

    /* compute period sums of y for time-demeaning */

    for (i=0; i<N; i++) {
	pt1 = t1 + (i + u0) * dset->pd;
	pt2 = t2 + (i + u0) * dset->pd;
	s = 0;
	for (t=pt1; t<=pt2; t++) {
	    yavg->val[s++] += dset->Z[vnum][t];
	}
    }

    gretl_matrix_divide_by_scalar(yavg, N);
    bigrow = 0;

    for (i=0; i<N && !err; i++) {
	double yti, yti_1;
	int p_i, T_i, k_i;
	int pt0, ss;

	if (p_varies) {
	    p_i = plist[i+1];
	    T_i = t2 - t1 + 1 - (1 + p_i);
	    k_i = m + p_i;
	    gretl_matrix_reuse(y, T_i, 1);
	    gretl_matrix_reuse(X, T_i, k_i);
	    gretl_matrix_reuse(b, k_i, 1);
	    gretl_matrix_reuse(ei, T_i, 1);
	    gretl_matrix_reuse(vi, T_i, 1);
	} else {
	    p_i = p;
	    T_i = T;
	    k_i = k;
	}

	/* indices into Z array */
	pt1 = t1 + (i + u0) * dset->pd;
	pt2 = t2 + (i + u0) * dset->pd;
	pt0 = pt1 + 1 + p_i;

	/* build (full length) \delta y_t in dy */
	s = 0;
	for (t=pt1+1; t<=pt2; t++) {
	    ss = t - pt1;
	    yti = dset->Z[vnum][t] - gretl_vector_get(yavg, ss);
	    yti_1 = dset->Z[vnum][t-1] - gretl_vector_get(yavg, ss-1);
	    gretl_vector_set(dy, s++, yti - yti_1);
	}

	/* build y_{t-1} in y */
	s = 0;
	for (t=pt0; t<=pt2; t++) {
	    yti_1 = dset->Z[vnum][t-1] - gretl_vector_get(yavg, t - pt1 - 1);
	    gretl_vector_set(y, s++, yti_1);
	}	

	/* augmented case: write lags of dy into X */
	for (j=1; j<=p_i; j++) {
	    int col = m + j - 2;
	    double dylag;

	    s = 0;
	    for (t=pt0; t<=pt2; t++) {
		dylag = gretl_vector_get(dy, t - pt1 - 1 - j);
		gretl_matrix_set(X, s++, col, dylag);
	    }
	}

	/* set lagged y as last regressor */
	for (t=0; t<T_i; t++) {
	    gretl_matrix_set(X, t, k_i - 1, y->val[t]);
	}

#if LLC_DEBUG > 1
	gretl_matrix_print(dy, "dy");
	gretl_matrix_print(y, "y1");
	gretl_matrix_print(X, "X");
#endif

	if (p_i > 0) {
	    /* "virtual trimming" of dy for regressions */
	    dy->val += p_i;
	    dy->rows -= p_i;
	}

	/* run (A)DF regression */
	err = gretl_matrix_ols(dy, X, b, NULL, ui, NULL);
	if (err) {
	    break;
	}

	if (k_i > 1) {
	    /* reduced regressor matrix for auxiliary regressions:
	       omit the last column containing the lagged level of y
	    */
	    gretl_matrix_reuse(X, T_i, k_i - 1);
	    gretl_matrix_reuse(b, k_i - 1, 1);

	    err = gretl_matrix_ols(dy, X, b, NULL, ei, NULL);
	    if (!err) {
		err = gretl_matrix_ols(y, X, b, NULL, vi, NULL);
	    }

	    gretl_matrix_reuse(X, T, k);
	    gretl_matrix_reuse(b, k, 1);
	} else {
	    /* no auxiliary regressions required */
	    gretl_matrix_copy_values(ei, dy);
	    gretl_matrix_copy_values(vi, y);
	}

	if (p_i > 0) {
	    /* restore dy to full length */
	    dy->val -= p_i;
	    dy->rows += p_i;
	}

	if (!err) {
	    double sui, s2yi, s2ui = 0.0;

	    for (t=0; t<T_i; t++) {
		s2ui += ui->val[t] * ui->val[t];
	    }

	    s2ui /= (T_i - 1);
	    sui = sqrt(s2ui);

	    /* write normalized per-unit ei and vi into big matrices */
	    gretl_matrix_divide_by_scalar(ei, sui);
	    gretl_matrix_divide_by_scalar(vi, sui);
	    gretl_matrix_inscribe_matrix(e, ei, bigrow, 0, GRETL_MOD_NONE);
	    gretl_matrix_inscribe_matrix(v, vi, bigrow, 0, GRETL_MOD_NONE);
	    bigrow += T_i;

	    s2yi = LLC_lrvar(dy, K, m, &err);
	    if (!err) {
		/* cumulate ratio of LR std dev to innovation std dev */
		SN += sqrt(s2yi) / sui;
	    }

#if LLC_DEBUG
	    pprintf(prn, "s2ui = %.8f, s2yi = %.8f\n", s2ui, s2yi);
#endif
	}

	if (p_varies) {
	    gretl_matrix_reuse(y, T, 1);
	    gretl_matrix_reuse(X, T, k);
	    gretl_matrix_reuse(b, k, 1);
	    gretl_matrix_reuse(ei, T, 1);
	    gretl_matrix_reuse(vi, T, 1);
	}	    
    }

    if (!err) {
	/* the final step: full-length regression of e on v */
	double ee = 0, vv = 0;
	double delta, s2e, STD, td;
	double mstar, sstar;

	gretl_matrix_reuse(b, 1, 1);
	err = gretl_matrix_ols(e, v, b, NULL, eps, NULL);

	if (!err) {
	    for (t=0; t<NT; t++) {
		ee += eps->val[t] * eps->val[t];
		vv += v->val[t] * v->val[t];
	    }

	    SN /= N;
	    delta = b->val[0];
	    s2e = ee / NT;
	    STD = sqrt(s2e) / sqrt(vv);
	    td = delta / STD;

	    /* fetch the Levin-Lin-Chu corrections factors */
	    err = get_LLC_corrections(T, m, &mstar, &sstar);
	}

	if (!err) {
	    double z = (td - NT * (SN / s2e) * STD * mstar) / sstar;
	    double pval = normal_cdf(z);

#if LLC_DEBUG
	    pprintf(prn, "mustar = %g, sigstar = %g\n", mstar, sstar);
	    pprintf(prn, "SN = %g, se = %g, STD = %g\n", SN, sqrt(s2e), STD);
#endif

	    if (!(opt & OPT_Q)) {
		const char *heads[] = {
		    N_("coefficient"),
		    N_("t-ratio"),
		    N_("z-score")
		};
		const char *s = dset->varname[vnum];
		char NTstr[32];
		int sp[3] = {0, 3, 5};
		int w[3] = {4, 6, 0};
 
		pputc(prn, '\n');
		pprintf(prn, _("Levin-Lin-Chu pooled ADF test for %s\n"), s);
		pprintf(prn, "%s ", _(DF_test_spec(m)));

		if (p_varies) {
		    pprintf(prn, _("including %.2f lags of (1-L)%s (average)"), pbar, s);
		} else if (p == 1) {
		    pprintf(prn, _("including one lag of (1-L)%s"), s);
		} else {
		    pprintf(prn, _("including %d lags of (1-L)%s"), p, s);
		}
		pputc(prn, '\n');

		pprintf(prn, _("Bartlett truncation at %d lags\n"), K);
		sprintf(NTstr, "N,T = (%d,%d)", N, dyT + 1);
		pprintf(prn, _("%s, using %d observations"), NTstr, NT);

		pputs(prn, "\n\n");
		for (i=0; i<3; i++) {
		    pputs(prn, _(heads[i]));
		    bufspace(w[i], prn);
		    w[i] = sp[i] + g_utf8_strlen(_(heads[i]), -1);
		}
		pputc(prn, '\n');

		pprintf(prn, "%*.5g %*.3f %*.6g [%.4f]\n\n", 
			w[0], delta, w[1], td, w[2], z, pval);
	    }

	    record_test_result(z, pval, "Levin-Lin-Chu");
	}
    }

    gretl_matrix_block_destroy(B);

    return err;
}
