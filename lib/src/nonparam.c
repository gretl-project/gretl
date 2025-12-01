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

/**
 * SECTION:nonparam
 * @short_description: selected nonparametric routines
 * @title: Nonparametric
 * @include: libgretl.h
 *
 * Includes various nonparametric tests: rank correlation,
 * runs (randomness), and differences between series.  Also
 * "loess" regression a la William Cleveland.
 */

/* Spearman, one-tailed alpha = .005, .025, .05 */

static double srhocrit[18][3] = {
    { 0.8929, 0.7450, 0.6786 }, /*  n = 7 */
    { 0.8571, 0.7143, 0.6190 }, /*  8 */
    { 0.8167, 0.6833, 0.5833 }, /*  9 */
    { 0.7818, 0.6364, 0.5515 }, /* 10 */
    { 0.7545, 0.6091, 0.5273 }, /* 11 */
    { 0.7273, 0.5804, 0.4965 }, /* 12 */
    { 0.6978, 0.5549, 0.4780 }, /* 13 */
    { 0.6747, 0.5341, 0.4593 }, /* 14 */
    { 0.6536, 0.5179, 0.4429 }, /* 15 */
    { 0.6324, 0.5000, 0.4265 }, /* 16 */
    { 0.6152, 0.4853, 0.4118 }, /* 17 */
    { 0.5975, 0.4716, 0.3994 }, /* 18 */
    { 0.5825, 0.4579, 0.3895 }, /* 19 */
    { 0.5684, 0.4451, 0.3789 }, /* 20 */
    { 0.5545, 0.4351, 0.3688 }, /* 21 */
    { 0.5426, 0.4241, 0.3597 }, /* 22 */
    { 0.5306, 0.4150, 0.3518 }, /* 23 */
    { 0.5200, 0.4061, 0.3435 }  /* 24 */
};

/* The values above are from Daniel and Terrell,
   Business Statistics, 4e (Houghton Mifflin, 1986)
*/

static double spearman_signif (double rho, int n)
{
    const double *x;

    if (n < 7) {
	return 1.0;
    }

    x = srhocrit[n - 7];

    if (rho > x[0]) return .01;
    if (rho > x[1]) return .05;
    if (rho > x[2]) return .10;

    return 1.0;
}

enum {
    RANK_X,
    RANK_Y
};

/* m = number of sorted values (sz); n = number of raw values (either
   x or y); in case of missing values in x and/or y we may have m < n.
*/

static void make_ranking (const double *sz, int m,
			  const double *x, const double *y,
			  int n, double *rz, int *ties,
			  int which)
{
    const double *z;
    int cases, k, i, j;
    double avg, r = 1;

    z = (which == RANK_X)? x : y;

    for (i=0; i<m; i++) {
	/* scan sorted z */
	cases = k = 0;

	if (i > 0 && sz[i] == sz[i-1]) {
	    continue;
	}

	for (j=0; j<n; j++) {
	    /* scan raw z for matches */
	    if (!na(x[j]) && !na(y[j])) {
		if (z[j] == sz[i]) {
		    rz[k] = r;
		    cases++;
		}
		k++;
	    }
	}

	if (cases > 1) {
	    avg = (r + r + cases - 1.0) / 2.0;
	    for (j=0; j<m; j++) {
		if (rz[j] == r) {
		    rz[j] = avg;
		}
	    }
	    if (ties != NULL) {
		*ties = 1;
	    }
	}

	r += cases;
    }
}

static int rankcorr_get_rankings (const double *x, const double *y, int n,
				  double **rxout, double **ryout,
				  int *pm, int *ties)
{
    double *sx = NULL, *sy = NULL;
    double *rx = NULL, *ry = NULL;
    int i, m = 0;

    /* count non-missing pairs */
    for (i=0; i<n; i++) {
	if (!na(x[i]) && !na(y[i])) {
	    m++;
	}
    }

    if (m < 2) {
	return E_DATA;
    }

    sx = malloc(m * sizeof *sx);
    sy = malloc(m * sizeof *sy);
    rx = malloc(m * sizeof *rx);
    ry = malloc(m * sizeof *ry);

    if (sx == NULL || sy == NULL ||
	rx == NULL || ry == NULL) {
	free(sx);
	free(sy);
	free(rx);
	free(ry);
	return E_ALLOC;
    }

    /* copy non-missing x and y into sx, sy */
    m = 0;
    for (i=0; i<n; i++) {
	if (!na(x[i]) && !na(y[i])) {
	    sx[m] = x[i];
	    sy[m] = y[i];
	    m++;
	}
    }

    /* get sorted series */
    qsort(sx, m, sizeof *sx, gretl_inverse_compare_doubles);
    qsort(sy, m, sizeof *sy, gretl_inverse_compare_doubles);

    for (i=0; i<m; i++) {
	rx[i] = ry[i] = 0.0;
    }

    /* make rankings by comparing "raw" x, y with sorted */
    make_ranking(sx, m, x, y, n, rx, ties, RANK_X);
    make_ranking(sy, m, x, y, n, ry, ties, RANK_Y);

    /* save the ranks */
    *rxout = rx;
    *ryout = ry;

    if (pm != NULL) {
	*pm = m;
    }

    free(sx);
    free(sy);

    return 0;
}

static int real_spearman_rho (const double *x, const double *y, int n,
			      double *rho, double *zval,
			      double **rxout, double **ryout,
			      int *pm)
{
    double *rx = NULL, *ry = NULL;
    int m, ties = 0;
    int err = 0;

    *rho = *zval = NADBL;

    if (n < 2) {
	return E_TOOFEW;
    }

    err = rankcorr_get_rankings(x, y, n, &rx, &ry, &m, &ties);
    if (err) {
	return err;
    }

    /* Pearson correlation in ranks */
    *rho = gretl_corr(0, m - 1, rx, ry, NULL);

    /* save the ranks, if wanted */
    if (rxout != NULL && ry != NULL) {
	*rxout = rx;
	*ryout = ry;
    } else {
	free(rx);
	free(ry);
    }

    *pm = m;

    return err;
}

gretl_matrix *spearman_rho_func (const double *x,
				 const double *y,
				 int n, int *err)
{
    gretl_matrix *ret = NULL;
    double z, rho = NADBL;
    int m = 0;

    *err = real_spearman_rho(x, y, n, &rho, &z,
			     NULL, NULL, &m);

    if (!*err) {
	ret = gretl_matrix_alloc(1, 3);
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else {
	    ret->val[0] = rho;
	    if (!na(z)) {
		ret->val[1] = z;
		ret->val[2] = normal_pvalue_2(z);
	    } else if (m > 24) {
		ret->val[1] = rho * sqrt((m - 2) / (1 - rho*rho));
		ret->val[2] = student_pvalue_2(m - 2, ret->val[1]);
	    } else {
		ret->val[1] = ret->val[2] = NADBL;
	    }
	}
    }

    return ret;
}

static void print_raw_and_ranked (int vx, int vy,
				  const double *x, const double *y,
				  const double *rx, const double *ry,
				  const DATASET *dset,
				  PRN *prn)
{
    int T = dset->t2 - dset->t1 + 1;
    int obslen = max_obs_marker_length(dset);
    int n1 = strlen(dset->varname[vx]);
    int n2 = strlen(dset->varname[vy]);
    int t, i = 0;

    n1 = (n2 > n1)? n2 : n1;
    n1 = (n1 < 15)? 15 : n1;
    pputc(prn, '\n');
    bufspace(obslen + 1, prn);
    pprintf(prn, "%*s%8s%*s%8s\n\n", n1, dset->varname[vx], _("rank"),
	    n1, dset->varname[vy], _("rank"));

    for (t=0; t<T; t++) {
	print_obs_marker(t + dset->t1, dset, obslen, prn);
	if (!(na(x[t])) && !(na(y[t]))) {
	    pprintf(prn, "%#*g", n1, x[t]);
	    pprintf(prn, "%8g", rx[i]);
	    pprintf(prn, "%#*g", n1, y[t]);
	    pprintf(prn, "%8g", ry[i]);
	    i++;
	}
	pputc(prn, '\n');
    }
}

/**
 * spearman_rho:
 * @list: list of (two) variables to process.
 * @dset: dataset struct.
 * @opt: if includes %OPT_V, print both the "raw" and the ranked data.
 * @prn: gretl printing struct.
 *
 * Calculates and prints Spearman's rank correlation coefficient for the two
 * variables specified in the @list.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int spearman_rho (const int *list, const DATASET *dset,
		  gretlopt opt, PRN *prn)
{
    int T = dset->t2 - dset->t1 + 1;
    double *rx = NULL, *ry = NULL;
    const double *x, *y;
    double zval, rho = 0;
    int vx, vy, m;
    int err;

    if (list == NULL || list[0] != 2) {
	pputs(prn, _("This command requires two variables\n"));
	return 1;
    }

    vx = list[1];
    vy = list[2];

    x = dset->Z[vx] + dset->t1;
    y = dset->Z[vy] + dset->t1;

    err = real_spearman_rho(x, y, T, &rho, &zval,
			    (opt & OPT_V)? &rx : NULL,
			    (opt & OPT_V)? &ry : NULL,
			    &m);
    if (err) {
	return err;
    }

    pprintf(prn, _("\nFor the variables '%s' and '%s'\n"), dset->varname[vx],
	    dset->varname[vy]);

    if (na(rho)) {
	pputs(prn, _("Spearman's rank correlation is undefined\n"));
	return err;
    }

    pprintf(prn, _("Spearman's rank correlation coefficient (rho) = %.8f\n"), rho);

    if (rho == 0.0) {
	goto skipit;
    }

    if (!na(zval)) {
	pputs(prn, _("Under the null hypothesis of no correlation:\n "));
	pprintf(prn, _("z-score = %g, with two-tailed p-value %.4f\n"), zval,
		normal_pvalue_2(zval));
    } else if (m > 24) {
	double tval = rho * sqrt((m - 2) / (1 - rho*rho));

	pputs(prn, _("Under the null hypothesis of no correlation:\n "));
	pprintf(prn, _("t(%d) = %g, with two-tailed p-value %.4f\n"), m - 2,
		tval, student_pvalue_2(m - 2, tval));
    } else if (m >= 7) {
	double pval = spearman_signif(fabs(rho), m);

	if (pval < 1.0) {
	    pprintf(prn, _("significant at the %g%% level (two-tailed)\n"),
		    100.0 * pval);
	} else {
	    /* xgettext:no-c-format */
	    pputs(prn, _("not significant at the 10% level\n"));
	}
    } else {
	pputs(prn, _("Sample is too small to calculate a p-value based on "
		     "the normal distribution\n"));
    }

 skipit:

    if (rx != NULL && ry != NULL) {
	print_raw_and_ranked(vx, vy, x, y, rx, ry, dset, prn);
	free(rx);
	free(ry);
    }

    pputc(prn, '\n');

    return 0;
}

struct xy_pair {
    double x;
    double y;
};

static int compare_pairs_x (const void *a, const void *b)
{
    const struct xy_pair *pa = a;
    const struct xy_pair *pb = b;
    int ret;

    ret = (pa->x > pb->x) - (pa->x < pb->x);
    if (ret == 0) {
	ret = (pa->y > pb->y) - (pa->y < pb->y);
    }

    return ret;
}

static int compare_pairs_y (const void *a, const void *b)
{
    const struct xy_pair *pa = a;
    const struct xy_pair *pb = b;

    return (pa->y > pb->y) - (pa->y < pb->y);
}

static int real_kendall_tau (const double *x, const double *y,
			     int n, struct xy_pair *xy, int nn,
			     double *ptau, double *pz)
{
    double tau, nn1, s2, z;
    int tt1, tx = 0, ty = 0;
    int Tx = 0, Ty = 0;
    int Tx2 = 0, Ty2 = 0;
    int Tx25 = 0, Ty25 = 0;
    int tie, N0, N1, S;
    int i, j;

    /* populate sorter */
    j = 0;
    for (i=0; i<n; i++) {
	if (!na(x[i]) && !na(y[i])) {
	    xy[j].x = x[i];
	    xy[j].y = y[i];
	    j++;
	}
    }

    /* sort pairs by x */
    qsort(xy, nn, sizeof *xy, compare_pairs_x);

    /* make order counts */
    N0 = N1 = 0;
    for (i=0; i<nn; i++) {
	for (j=i+1; j<nn; j++) {
	    if (xy[j].y == xy[i].y) {
		Ty++;
	    } else if (xy[j].x != xy[i].x) {
		if (xy[j].y > xy[i].y) {
		    N0++;
		} else if (xy[j].y < xy[i].y) {
		    N1++;
		}
	    }
	}
	if (i > 0) {
	    /* account for ties in x */
	    tie = (xy[i].x == xy[i-1].x);
	    if (tie) {
		tx++;
	    }
	    if (tx > 0 && (!tie || i == nn - 1)) {
		tx++;
		tt1 = tx * (tx - 1);
		Tx += tt1;
		Tx2 += tt1 * (tx - 2);
		Tx25 += tt1 * (2 * tx + 5);
		tx = 0;
	    }
	}
    }

    if (Ty > 0) {
	/* account for ties in y */
	Ty = 0;
	qsort(xy, nn, sizeof *xy, compare_pairs_y);
	for (i=1; i<nn; i++) {
	    tie = (xy[i].y == xy[i-1].y);
	    if (tie) {
		ty++;
	    }
	    if (ty > 0 && (!tie || i == nn - 1)) {
		ty++;
		tt1 = ty * (ty - 1);
		Ty += tt1;
		Ty2 += tt1 * (ty - 2);
		Ty25 += tt1 * (2 * ty + 5);
		ty = 0;
	    }
	}
    }

    S = N0 - N1;

#if 0
    fprintf(stderr, "N0 = %d, N1 = %d, S = %d\n", N0, N1, S);
    fprintf(stderr, "Tx = %d, Ty = %d\n", Tx, Ty);
#endif

    nn1 = nn * (nn - 1.0);

    /* normal approximation as in Shapiro and Chen,
       Journal of Quality Technology, 2001 */

    if (Tx == 0 && Ty == 0) {
	tau = 2 * S / nn1;
	s2 = (1.0/18) * nn1 * (2 * nn + 5);
    } else {
	double den = (nn1 - Tx) * (nn1 - Ty);

	tau = 2 * S / sqrt(den);
	s2 = (1.0/18) * (nn1 * (2 * nn + 5) - Tx25 - Ty25);
	if (Tx2 != 0 && Ty2 != 0) {
	    s2 += (1.0/(9*nn1*(nn-2))) * Tx2 * Ty2;
	}
	if (Tx != 0 && Ty != 0) {
	    s2 += (1.0/(2*nn1)) * Tx * Ty;
	}
    }

    z = (S - 1) / sqrt(s2);

    if (ptau != NULL) {
	*ptau = tau;
    }

    if (pz != NULL) {
	*pz = z;
    }

    return 0;
}

gretl_matrix *kendall_tau_func (const double *x,
				const double *y,
				int n, int *err)
{
    gretl_matrix *ret = NULL;
    struct xy_pair *xy;
    double z, tau = NADBL;
    int i, nn = 0;

    /* count valid pairs */
    for (i=0; i<n; i++) {
	if (!na(x[i]) && !na(y[i])) {
	    nn++;
	}
    }

    if (nn < 2) {
	*err = E_TOOFEW;
	return NULL;
    }

    xy = malloc(nn * sizeof *xy);

    if (xy == NULL) {
	*err = E_ALLOC;
    } else {
	*err = real_kendall_tau(x, y, n, xy, nn, &tau, &z);
    }

    free(xy);

    if (!*err) {
	ret = gretl_matrix_alloc(1, 3);
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else {
	    ret->val[0] = tau;
	    ret->val[1] = z;
	    ret->val[2] = normal_pvalue_2(z);
	}
    }

    return ret;
}

/**
 * kendall_tau:
 * @list: list of (two) variables to process.
 * @dset: dataset struct.
 * @opt: if includes %OPT_V, print both the "raw" and the ranked data.
 * @prn: gretl printing struct.
 *
 * Calculates and prints Kendall's rank correlation tau statistic for
 * the two variables specified in @list.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int kendall_tau (const int *list, const DATASET *dset,
		 gretlopt opt, PRN *prn)
{
    struct xy_pair *xy = NULL;
    int T = dset->t2 - dset->t1 + 1;
    const double *x, *y;
    double tau, z;
    int vx, vy;
    int t, nn = 0;
    int err = 0;

    if (list == NULL || list[0] != 2) {
	pputs(prn, _("This command requires two variables\n"));
	return 1;
    }

    vx = list[1];
    vy = list[2];

    x = dset->Z[vx] + dset->t1;
    y = dset->Z[vy] + dset->t1;

    /* count valid pairs */
    for (t=0; t<T; t++) {
	if (!na(x[t]) && !na(y[t])) {
	    nn++;
	}
    }

    if (nn < 2) {
	return E_TOOFEW;
    }

    /* allocate */
    xy = malloc(nn * sizeof *xy);
    if (xy == NULL) {
	return E_ALLOC;
    }

    /* calculate */
    err = real_kendall_tau(x, y, T, xy, nn, &tau, &z);

    if (!err) {
	pprintf(prn, _("\nFor the variables '%s' and '%s'\n"), dset->varname[vx],
		dset->varname[vy]);
	pprintf(prn, "%s = %.8f\n", _("Kendall's tau"), tau);
	pputs(prn, _("Under the null hypothesis of no correlation:\n "));
	pprintf(prn, _("z-score = %g, with two-tailed p-value %.4f\n"), z,
		normal_pvalue_2(z));
    }

    if (opt & OPT_V) {
	double *rx = NULL, *ry = NULL;

	rankcorr_get_rankings(x, y, T, &rx, &ry, NULL, NULL);

	if (rx != NULL && ry != NULL) {
	    print_raw_and_ranked(vx, vy, x, y, rx, ry, dset, prn);
	    free(rx);
	    free(ry);
	}
    }

    pputc(prn, '\n');

    free(xy);

    return err;
}

#define LOCKE_DEBUG 0

static int randomize_doubles (const void *a, const void *b)
{
    return gretl_rand_int_max(8096) - 4097;
}

/* check for negative values, screen out missing values,
   and load x into an array for randomization */

static int locke_shuffle_init (const double *x,
			       int *n, double **psx)
{
    double *sx = NULL;
    int i, m = 0;

    for (i=0; i<*n; i++) {
	if (x[i] < 0.0) {
	    return E_DATA;
	}
	if (!na(x[i])) {
	    m++;
	}
    }

    if (m < 4) {
	return E_DATA;
    }

    sx = malloc(m * sizeof *sx);
    if (sx == NULL) {
	return E_ALLOC;
    }

    *psx = sx;

    m = 0;
    for (i=0; i<*n; i++) {
	if (!na(x[i])) {
	    sx[m++] = x[i];
	}
    }

    m = 2 * m / 2;
    *n = m;

    return 0;
}

#define NREPEAT 100

/**
 * lockes_test:
 * @x: data series.
 * @t1: start of sample range.
 * @t2: end of sample range.
 * @err: location to receive error code.
 *
 * Performs Charles Locke's nonparametric test for whether an
 * empirical distribution (namely, that of @x over the range
 * @t1 to @t2) is gamma.  See C. Locke, "A Test for the Composite
 * Hypothesis that a Population has a Gamma Distribution,"
 * Commun. Statis.-Theor. Meth. A5(4), 351-364 (1976).  Also
 * see Shapiro and Chen, Journal of Quality Technology 33(1),
 * Jan 2001.
 *
 * Returns: the z value for test, or #NADBL on error.
 */

double lockes_test (const double *x, int t1, int t2, int *err)
{
    struct xy_pair *uv = NULL;
    double *sx = NULL, *u = NULL, *v = NULL;
    double zj, z;
    int m = t2 - t1 + 1;
    int i, j, t;

    *err = locke_shuffle_init(x + t1, &m, &sx);
    if (*err) {
	return NADBL;
    }

    m /= 2;
    u = malloc(m * sizeof *u);
    v = malloc(m * sizeof *v);
    uv = malloc(m * sizeof *uv);

    if (u == NULL || v == NULL || uv == NULL) {
        *err = E_ALLOC;
	z = NADBL;
	goto bailout;
    }

    z = 0.0;

    /* repeat the shuffling of the series NREPEAT times, since the
       test statistic is sensitive to the ordering under the null
    */

    for (j=0; j<NREPEAT; j++) {
#if LOCKE_DEBUG
	fprintf(stderr, "locke's test: j = %d\n", j);
#endif
	qsort(sx, 2 * m, sizeof *sx, randomize_doubles);
	t = 0;
	for (i=0; i<m; i++) {
	    u[i] = sx[t] + sx[t+1];
	    v[i] = sx[t] / sx[t+1];
	    if (sx[t+1] / sx[t] > v[i]) {
		v[i] = sx[t+1] / sx[t];
	    }
	    t += 2;
	}
	real_kendall_tau(u, v, m, uv, m, NULL, &zj);
	z += zj;
#if LOCKE_DEBUG
	fprintf(stderr, "z[%d] = %g\n", j, zj);
#endif
    }

    z /= (double) NREPEAT;

#if LOCKE_DEBUG
    fprintf(stderr, "Kendall's tau: average z = %g\n", z);
#endif

 bailout:

    free(u);
    free(v);
    free(uv);
    free(sx);

    return z;
}

/**
 * vge_gamma_test:
 * @x: data series.
 * @t1: start of sample range.
 * @t2: end of sample range.
 * @err: location to receive error code.
 *
 * Performs the test described by J. A. Villaseñor and E. González-Estrada
 * (Statistics and Probability Letters, 96 (2015) pp. 281–286) for the null
 * hypothesis that @x is gamma-distributed over the range  @t1 to @t2.
 *
 * For the sake of compatibility with the gamma_test() function in the R
 * package named "goft" we divide by n-1 in computing the covariance
 * term @sxz, although this is not recommended by the authors of the above-
 * noted publication.
 *
 * Returns: the z value for the test, or #NADBL on error.
 */

double vge_gamma_test (const double *x, int t1, int t2, int *err)
{
    int t, n = t2 - t1 + 1;
    double xbar = 0, zbar = 0;
    double xc, zc, sxz, s2;
    double a, V, Vstar;

    for (t=t1; t<=t2; t++) {
        if (x[t] <= 0) {
            gretl_errmsg_set(_("Non-positive values encountered"));
            *err = E_DATA;
            return NADBL;
        } else if (na(x[t])) {
	    n--;
	} else {
            xbar += x[t];
            zbar += log(x[t]);
        }
    }

    if (n < 30) {
        /* minimum obs? */
        *err = E_TOOFEW;
        return NADBL;
    }

    xbar /= n;
    zbar /= n;

    sxz = s2 = 0;
    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            xc = x[t] - xbar;
            zc = log(x[t]) - zbar;
            sxz += xc * zc;
            s2 += xc * xc;
        }
    }

    sxz /= (n-1);
    s2 /= (n-1);
    V = s2 / (xbar * sxz);
    a = xbar / sxz;
    Vstar = sqrt(n*a) * (V-1);

    return fabs(Vstar) / sqrt(2.0);
}

/**
 * runs_test:
 * @v: ID number of the variable to process.
 * @dset: dataset struct.
 * @opt: %OPT_D to use first difference of variable, %OPT_E
 * to assume positive and negative are equiprobable.
 * @prn: gretl printing struct.
 *
 * Performs, and prints the results of, the runs test for randomness
 * for the variable specified by @v.  The normal approximation
 * is that given in Gary Smith, Statistical Reasoning, 2e, p. 674.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int runs_test (int v, const DATASET *dset,
	       gretlopt opt, PRN *prn)
{
    double xt, *x, mu, s2, sigma;
    double N2, z, pval;
    int Np, Nm;
    int t, n, runs = 1;

    n = dset->t2 - dset->t1 + 1;

    x = malloc(n * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    n = 0;

    if (opt & OPT_D) {
	double xt1;

	for (t=dset->t1 + 1; t<=dset->t2; t++) {
	    xt = dset->Z[v][t];
	    xt1 = dset->Z[v][t-1];
	    if (na(xt) || na(xt1)) {
		continue;
	    }
	    x[n++] = xt - xt1;
	}
    } else {
	for (t=dset->t1; t<=dset->t2; t++) {
	    xt = dset->Z[v][t];
	    if (na(xt)) {
		continue;
	    }
	    x[n++] = xt;
	}
    }

    if (n <= 1) {
	free(x);
	return E_TOOFEW;
    }

    Np = (x[0] > 0);
    Nm = 1 - Np;

    for (t=1; t<n; t++) {
	if (x[t] > 0) {
	    Np++;
	} else {
	    Nm++;
	}
	if ((x[t] > 0 && x[t-1] <= 0) || (x[t] <= 0 && x[t-1] > 0)) {
	    runs++;
	}
    }

    if (opt & OPT_E) {
	mu = 1.0 + n / 2.0;
	s2 = ((double) n - 1) / 4.0;
    } else {
	/* don't assume that + and - are equiprobable */
	N2 = 2.0 * Np * Nm;
	mu = 1.0 + N2 / n;
	s2 = (N2 * (N2 - n)) / (n*n * (n-1));
    }

    if (s2 <= 0) {
	sigma = 0;
	z = pval = NADBL;
    } else {
	sigma = sqrt(s2);
	z = (runs - mu) / sigma;
	pval = normal_pvalue_2(z);
    }

    if (opt & OPT_D) {
	pprintf(prn, "\n%s\n", _("Runs test (first difference)"));
    } else {
	pprintf(prn, "\n%s\n", _("Runs test (level)"));
    }

    pprintf(prn, _("\nNumber of runs (R) in the variable '%s' = %d\n"),
	    dset->varname[v], runs);

    if (na(z)) {
	pprintf(prn, _("Test statistic cannot be computed: try "
		       "the deviation from the median?\n"));
    } else {
	if (opt & OPT_E) {
	    pprintf(prn, _("Under the null hypothesis of independence and equal "
			   "probability of positive\nand negative values, R "
			   "follows N(%g, %g)\n"), mu, sigma);
	} else {
	    pprintf(prn, _("Under the null hypothesis of independence, R "
			   "follows N(%g, %g)\n"), mu, sigma);
	}
	pprintf(prn, _("z-score = %g, with two-tailed p-value %g\n"), z, pval);
    }

    pputc(prn, '\n');

    record_test_result(z, pval);

    free(x);

    return 0;
}

static double print_z_prob (double z, PRN *prn)
{
    double pv = NADBL;

    if (z > 0) {
	pv = normal_pvalue_1(z);
	if (!na(pv)) {
	    pprintf(prn, "  P(Z > %g) = %g\n", z, pv);
	}
    } else if (z < 0) {
	pv = normal_cdf(z);
	if (!na(pv)) {
	    pprintf(prn, "  P(Z < %g) = %g\n", z, pv);
	}
    }

    if (!na(pv) && pv > 0) {
	pprintf(prn, "  %s = %g\n", _("Two-tailed p-value"), 2 * pv);
    }

    return pv;
}

static void diff_test_header (int v1, int v2,
                              const DATASET *dset,
                              gretlopt opt,
			      PRN *prn)
{
    pputc(prn, '\n');
    if (opt & OPT_D) {
        char term1[72];
        char term2[72];

        sprintf(term1, "%s(%s==0)", dset->varname[v1], dset->varname[v2]);
        sprintf(term2, "%s(%s==1)", dset->varname[v1], dset->varname[v2]);
        pprintf(prn, _("Test for difference between %s and %s"),
                term1, term2);
    } else {
        pprintf(prn, _("Test for difference between %s and %s"),
                dset->varname[v1], dset->varname[v2]);
    }
}

static const int rank5[3][2] = {
    { 0, 2 }, /* n = 6 */
    { 2, 3 }, /* n = 7 */
    { 3, 5 }  /* n = 8 */
};

struct ranker {
    double val;
    double rank;
    char c;
};

/* Wilcoxon signed-rank test, with handling of zero-differences and
   non-zero ties, plus continuity correction, as in E. Cureton, "The
   Normal Approximation to the Signed-Rank Sampling Distribution when
   Zero Differences are Present", JASA(62), 1967, 1068-1069.
*/

int signed_rank_test (const double *x, const double *y,
                      int v1, int v2, const DATASET *dset,
                      double *result, gretlopt opt, PRN *prn)
{
    struct ranker *r;
    double d, T, wp, wm;
    double z = NADBL, pval = NADBL;
    int quiet = (opt & OPT_Q);
    int Z = 0, N = 0;
    int i, k, t, n = 0;

    for (t=dset->t1; t<=dset->t2; t++) {
	if (!na(x[t]) && !na(y[t])) {
	    if (x[t] != y[t]) {
		n++;
	    }
	    N++;
	}
    }

    if (N == 0) {
	return E_MISSDATA;
    }

    Z = N - n; /* number of zero-differences */

    r = malloc(n * sizeof *r);
    if (r == NULL) {
	return E_ALLOC;
    }

    i = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
	if (!na(x[t]) && !na(y[t]) && x[t] != y[t]) {
	    d = x[t] - y[t];
	    r[i].val = fabs(d);
	    r[i].c = (d > 0)? '+' : '-';
	    i++;
	}
    }

    qsort(r, n, sizeof *r, gretl_compare_doubles);

    T = 0.0; /* non-zero ties correction */
    k = 0;   /* number of non-zero ties */

    for (i=0; i<n; i++) {
	int m = 0;

	for (t=i+1; t<n && r[t].val == r[i].val; t++) {
	    m++;
	}

	if (m == 0) {
	    r[i].rank = Z + i + 1;
	} else {
	    double avg = (Z + i + 1 + Z + i + m + 1) / 2.0;

	    for (t=0; t<=m; t++) {
		r[i+t].rank = avg;
	    }
	    i += m;
	    T += pow((double) m, 3) - m;
	    k++;
	}
    }

    if (!quiet) {
	diff_test_header(v1, v2, dset, opt, prn);
	pprintf(prn, "\n%s\n", _("Wilcoxon Signed-Rank Test"));
	pprintf(prn, "%s: %s\n\n", _("Null hypothesis"),
		_("the median difference is zero"));
	if (opt & OPT_V) {
	    pprintf(prn, "%16s %8s %16s\n\n", _("difference"), _("rank"),
		    _("signed rank"));
	}
    }

    wp = wm = 0.0;

    for (i=0; i<n; i++) {
	d = r[i].rank;
	if (r[i].c == '-') {
	    wm += d;
	    d = -d;
	    r[i].val = -r[i].val;
	} else {
	    wp += d;
	}
	if (opt & OPT_V) {
	    pprintf(prn, "%16g %8g %16g\n", r[i].val, r[i].rank, d);
	}
    }

    if (opt & OPT_V) {
	pputc(prn, '\n');
    }

    if (!quiet) {
	pprintf(prn, "  n = %d\n", n);
	pprintf(prn, "  W+ = %g, W- = %g\n", wp, wm);
	pprintf(prn, "  (%s: %d, %s: %d)\n", _("zero differences"),
		Z, _("non-zero ties"), k);
	if (N > 0 && n == 0) {
	    pprintf(prn, "  %s\n", _("Fail to reject H0"));
	}
    }

    if (n > 8) {
	double s2, x, num;

	x = (N*(N+1) - Z*(Z+1)) / 4.0;
	s2 = (N*(N+1)*(2*N+1) - Z*(Z+1)*(2*Z+1) - T/2) / 24.0;
	if (!quiet) {
	    pprintf(prn, "  %s = %g\n", _("Expected value"), x);
	    pprintf(prn, "  %s = %g\n", _("Variance"), s2);
	}
	num = wp - x;
	if (num > 0.25) {
	    num -= .5;
	} else {
	    num += .5;
	}
	z = num / sqrt(s2);
	if (quiet) {
	    pval = print_z_prob(z, NULL);
	} else {
	    pprintf(prn, "  z = %g\n", z);
	    pval = print_z_prob(z, prn);
	}
    } else if (n > 5) {
	pprintf(prn, _("  5%% %s: %d (%s), %d (%s)\n"),
		_("critical values"),
                rank5[n-6][0], _("two-tailed"),
                rank5[n-6][1], _("one-tailed"));
    } else {
	pprintf(prn, "  %s\n",
		_("Sample too small for statistical significance"));
    }

    if (!quiet) {
	pputc(prn, '\n');
    }

    result[0] = z;
    result[1] = pval;

    free(r);

    return 0;
}

int rank_sum_test (const double *x, const double *y,
                   int v1, int v2, const DATASET *dset,
                   double *result, gretlopt opt, PRN *prn)
{
    struct ranker *r;
    double wa, z = NADBL, pval = NADBL;
    char xc = 'a', yc = 'b';
    int na = 0, nb = 0;
    int quiet = (opt & OPT_Q);
    int i, t, n = 0;

    for (t=dset->t1; t<=dset->t2; t++) {
	if (!na(x[t])) na++;
	if (!na(y[t])) nb++;
    }

    if (na == 0 || nb == 0) {
	return E_MISSDATA;
    }

    n = na + nb;
    r = malloc(n * sizeof *r);
    if (r == NULL) {
	return E_ALLOC;
    }

    if (na > nb) {
	/* Make 'a' the smaller group */
	int tmp = na;

	na = nb;
	nb = tmp;
	xc = 'b';
	yc = 'a';
    }

    i = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
	if (!na(x[t])) {
	    r[i].val = x[t];
	    r[i].c = xc;
	    i++;
	}
	if (!na(y[t])) {
	    r[i].val = y[t];
	    r[i].c = yc;
	    i++;
	}
    }

    qsort(r, n, sizeof *r, gretl_compare_doubles);

    for (i=0; i<n; i++) {
	int m = 1;

	for (t=i+1; t<n && r[t].val == r[i].val; t++) {
	    m++;
	}

	if (m == 1) {
	    r[i].rank = i + 1;
	} else {
	    for (t=0; t<m; t++) {
		r[i+t].rank = (i + 1 + i + m) / 2.0;
	    }
	    i += m - 1;
	}
    }

    if (!quiet) {
	diff_test_header(v1, v2, dset, opt, prn);
	pprintf(prn, "\n%s\n", _("Wilcoxon Rank-Sum Test"));
	pprintf(prn, "%s: %s\n\n", _("Null hypothesis"),
		_("the two medians are equal"));
    }

    wa = 0.0;
    for (i=0; i<n; i++) {
	if (r[i].c == 'a') {
	    wa += r[i].rank;
	}
    }

    if (opt & OPT_V) {
	pprintf(prn, "%10s %7s %8s\n\n", "value", "rank", "group");
	for (i=0; i<n; i++) {
	    pprintf(prn, "%10g %7g %8c\n", r[i].val, r[i].rank, r[i].c);
	}
 	pputc(prn, '\n');
    }

    if (!quiet) {
	pprintf(prn, "  n1 = %d, n2 = %d\n", na, nb);
	pprintf(prn, "  w (%s) = %g\n", _("sum of ranks, sample 1"), wa);
    }

    if (na >= 10 && nb >= 10) {
	double m, s;

	m = na * (na + nb + 1) / 2.0;
	s = sqrt(na * nb * (na + nb + 1) / 12.0);
	z = (wa - m) / s;
	if (quiet) {
	    pval = print_z_prob(z, NULL);
	} else {
	    pprintf(prn, "  z = (%g - %g) / %g = %g\n", wa, m, s, z);
	    pval = print_z_prob(z, prn);
	}
    } else if (na >= 4 && nb >= 4 && nb <= 12) {
	void (*cv) (int, int, PRN *);

	cv = get_plugin_function("rank_sum_lookup");
	if (cv != NULL) {
	    (*cv)(na, nb, prn);
	}
    } else {
	pprintf(prn, "  %s\n",
		_("Sample too small for statistical significance"));
    }

    if (!quiet) {
	pputc(prn, '\n');
    }

    result[0] = z;
    result[1] = pval;

    free(r);

    return 0;
}

int sign_test (const double *x, const double *y,
	      int v1, int v2, const DATASET *dset,
	      double *result, gretlopt opt, PRN *prn)
{
    double pv;
    int n, w, t;

    n = w = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
	if (!na(x[t]) && !na(y[t]) && x[t] != y[t]) {
	    w += (x[t] > y[t]);
	    n++;
	}
    }

    if (n == 0) {
	return E_MISSDATA;
    }

    if (w == 0) {
	pv = 1.0;
    } else {
	pv = binomial_cdf_comp(0.5, n, w - 1);
    }

    if (!(opt & OPT_Q)) {
	diff_test_header(v1, v2, dset, opt, prn);
	pprintf(prn, "\n%s\n\n", _("Sign Test"));
	pprintf(prn, _("Number of differences: n = %d\n"), n);
	pputs(prn, "  ");
        pprintf(prn, _("Number of cases with %s > %s: w = %d (%.2f%%)\n"),
                dset->varname[v1], dset->varname[v2],
                w, 100.0 * w / n);
	pputs(prn, "  ");
        pprintf(prn, _("Under the null hypothesis of no difference, W "
                       "follows B(%d, %.1f)\n"), n, 0.5);
	pprintf(prn, "  %s(W <= %d) = %g\n", _("Prob"), w,
		binomial_cdf(0.5, n, w));
	pprintf(prn, "  %s(W >= %d) = %g\n\n", _("Prob"), w, pv);
    }

    result[0] = w;
    result[1] = pv;

    return 0;
}

struct pair_sorter {
    double xi;
    double yi;
    int i;
    char *s;
};

/**
 * sort_pairs_by_x:
 * @x: data vector by which to sort.
 * @y: data vector.
 * @order: location to receive sort order, or %NULL.
 * @labels: array of strings to be sorted along with
 * the data, or %NULL.
 *
 * Orders the elements of @x and @y by increasing value
 * of @x.  Optionally, returns in @order an array of
 * integers representing the order in which the original
 * observations appear in the sorted vectors.  Also
 * optionally sorts an accomanying array of observation
 * labels.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int sort_pairs_by_x (gretl_matrix *x, gretl_matrix *y, int **order,
		     char **labels)
{
    struct pair_sorter *s = NULL;
    int t, T = gretl_vector_get_length(x);
    int err = 0;

    if (T == 0 || gretl_vector_get_length(y) != T) {
	return E_NONCONF;
    }

    s = malloc(T * sizeof *s);
    if (s == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<T; t++) {
	s[t].xi = x->val[t];
	s[t].yi = y->val[t];
	s[t].i = t;
	if (labels != NULL) {
	    s[t].s = labels[t];
	} else {
	    s[t].s = NULL;
	}
    }

    qsort(s, T, sizeof *s, gretl_compare_doubles);

    for (t=0; t<T; t++) {
	x->val[t] = s[t].xi;
	y->val[t] = s[t].yi;
	if (labels != NULL) {
	    labels[t] = s[t].s;
	}
    }

    if (order != NULL) {
	int *idx = malloc(T * sizeof *idx);

	if (idx == NULL) {
	    err = E_ALLOC;
	} else {
	    for (t=0; t<T; t++) {
		idx[t] = s[t].i;
	    }
	    *order = idx;
	}
    }

    free(s);

    return err;
}

static int data_pre_sorted (const gretl_matrix *x)
{
    int t, T = gretl_vector_get_length(x);

    for (t=1; t<T; t++) {
	if (x->val[t] < x->val[t-1]) {
	    return 0;
	}
    }

    return 1;
}

#define LDEBUG 0

/* convenience struct for passing loess data */

struct loess_info {
    const gretl_matrix *y;
    const gretl_matrix *x;
    gretl_matrix *Xi;
    gretl_matrix *yi;
    gretl_matrix *wt;
    int d;
    int n;
    int N;
    int n_ok;
    int loo;
};

/* Compute robustness weights for loess, if wanted: on input @rw
   contains the N residuals; on return it contains the robustness
   weights as a function of the residuals. For observations
   where y is missing, the weight is NADBL.
*/

static int make_robustness_weights (gretl_matrix *rw, int N)
{
    double *tmp = malloc(N * sizeof *tmp);
    double s, eis;
    int i, k, err = 0;

    if (tmp == NULL) {
	return E_ALLOC;
    }

    /* pack the absolute values of the non-missing residuals
       into tmp */
    k = 0;
    for (i=0; i<N; i++) {
	if (!na(rw->val[i])) {
	    tmp[k++] = fabs(rw->val[i]);
	}
    }

    /* find the median absolute residual */
    s = gretl_median(0, k-1, tmp);
    free(tmp);

    if (na(s)) {
	err = E_DATA;
    } else {
	s *= 6;
	for (i=0; i<N && !err; i++) {
	    if (!na(rw->val[i])) {
		eis = rw->val[i] / s;
		if (fabs(eis) < 1.0) {
		    /* apply bisquare */
		    rw->val[i] = (1 - eis * eis) * (1 - eis * eis);
		} else {
		    rw->val[i] = 0.0;
		}
	    }
	}
    }

    return err;
}

/* Multiply the robustness weights, @rw, into the regular
   weights, @wt, taking care to register the two vectors
   (rw is full-length while wt is local), and to skip any
   missing values in rw.
*/

static void adjust_weights (const gretl_matrix *rw,
			    gretl_matrix *wt, int a)
{
    int k, n = gretl_vector_get_length(wt);
    int t = a;

    for (k=0; k<n; k++) {
	wt->val[k] *= rw->val[t];
	if (k < n - 1) {
	    t++;
	    while (na(rw->val[t])) t++;
	}
    }
}

/* Apply the weights in lo->wt to the local data matrices
   lo->Xi and lo->yi. This is straightforward since all
   the matrices have lo->n rows.
*/

static void weight_local_data (struct loess_info *lo)
{
    double xkj, wk;
    int k, j;

    for (k=0; k<lo->n; k++) {
	wk = sqrt(lo->wt->val[k]);
	for (j=0; j<lo->Xi->cols; j++) {
	    xkj = gretl_matrix_get(lo->Xi, k, j);
	    gretl_matrix_set(lo->Xi, k, j, xkj * wk);
	}
	lo->yi->val[k] *= wk;
    }

#if LDEBUG > 1
    gretl_matrix_print(lo->Xi, "Xi, weighted");
    gretl_matrix_print(lo->yi, "yi, weighted");
#endif
}

static int next_ok_obs (const double *y, int i)
{
    i++;
    while (na(y[i])) i++;
    return i;
}

static int loess_get_local_data (int i, int *pa,
				 struct loess_info *lo,
				 int *xconst)
{
    const double *x = lo->x->val;
    const double *y = lo->y->val;
    double xk, xk1, xds, h = 0;
    int n = lo->n, n_ok = lo->n_ok;
    int k_skip = -1;
    int a, b, k, m, t;
    int err = 0;

    gretl_matrix_zero(lo->Xi);
    gretl_matrix_zero(lo->yi);
    gretl_matrix_zero(lo->wt);

    a = *pa;

#if LDEBUG
    fprintf(stderr, "\ni=%d, a=%d, n_ok=%d\n", i, a, n_ok);
#endif

    /* First determine where we should start reading the
       neighbors of xi: search rightward from a, so far as
       this is feasible.
    */

    while (n_ok > n) {
	/* find b, the next point that would be included if
	   we move the neighbor set one place to the right:
	   start at a+1 and proceed until we have n points
	   with valid y-values
	*/
	for (m=0, b=a+1; ; b++) {
	    if (!na(y[b]) && ++m == n) {
		break;
	    }
	}
	if (fabs(x[i] - x[a]) < fabs(x[i] - x[b])) {
	    /* the max distance increased: get out */
	    break;
	}
	/* shift one (valid) place to the right */
	a = next_ok_obs(y, a);
	n_ok--;
    }

    /* record status for next round */
    *pa = a;
    lo->n_ok = n_ok;

    /* Having found the starting index for the n nearest neighbors
       of xi, transcribe the relevant data into Xi and yi. As we
       go, check whether x is constant in this sub-sample.
    */

    *xconst = 1;
    xk1 = 0;
    t = a;

    for (k=0; k<n; k++) {
	if (lo->loo && t == i) {
	    /* leave-one-out: mark this observation */
	    k_skip = k;
	}
	lo->yi->val[k] = y[t];
	xk = x[t];
	gretl_matrix_set(lo->Xi, k, 0, 1.0);
	gretl_matrix_set(lo->Xi, k, 1, xk);
	if (lo->Xi->cols > 2) {
	    gretl_matrix_set(lo->Xi, k, 2, xk * xk);
	}
	if (*xconst && k > 0 && xk > xk1) {
	    *xconst = 0;
	}
	xk1 = xk;
	if (k < n - 1) {
	    t = next_ok_obs(y, t);
	}
    }

#if LDEBUG
    fprintf(stderr, " -> a=%d, n_ok=%d (xconst = %d)\n",
	    a, n_ok, *xconst);
#endif

    if (!err) {
	/* find the max(abs) distance from xi */
	double h0 = fabs(x[i] - gretl_matrix_get(lo->Xi, 0, 1));
	double hn = fabs(x[i] - gretl_matrix_get(lo->Xi, n-1, 1));

	h = (h0 > hn)? h0 : hn;
    }

    /* compute scaled distances and tricube weights */
    for (k=0; k<n; k++) {
	if (k == k_skip) {
	    /* exclude this obs via a zero weight? */
	    lo->wt->val[k] = 0.0;
	    continue;
	}
	xk = gretl_matrix_get(lo->Xi, k, 1);
	if (h == 0.0) {
	    lo->wt->val[k] = 1.0;
	} else {
	    xds = fabs(x[i] - xk) / h;
	    if (xds < 1.0) {
		lo->wt->val[k] = pow(1.0 - pow(xds, 3.0), 3.0);
	    }
	}
#if LDEBUG > 1
	fprintf(stderr, "y=%10g, x=%10g, dist=%10g\n",
		lo->yi->val[k], xk, fabs(x[i] - xk));
#endif
    }

#if LDEBUG > 1
    gretl_matrix_print(lo->Xi, "Xi");
    gretl_matrix_print(lo->yi, "yi");
    gretl_matrix_print(lo->wt, "wt");
#endif

    return err;
}

static int loess_count_usable_obs (const gretl_matrix *y,
				   int N, int *amin)
{
    int i, n_ok = 0;

    for (i=0; i<N; i++) {
	if (!na(y->val[i])) {
	    if (n_ok == 0) {
		/* record the first usable obs index */
		*amin = i;
	    }
	    n_ok++;
	}
    }

    return n_ok;
}

/**
 * loess_fit:
 * @x: x-axis variable (must be pre-sorted).
 * @y: response variable.
 * @d: order for polynomial fit (0 <= d <= 2).
 * @q: bandwidth (0 < q <= 1).
 * @opt: give %OPT_R for robust variant (with re-weighting based on
 * the first-stage residuals).
 * @err: location to receive error code.
 *
 * Computes loess estimates based on William Cleveland, "Robust Locally
 * Weighted Regression and Smoothing Scatterplots", Journal of the
 * American Statistical Association, Vol. 74 (1979), pp. 829-836.
 * Typically one expects that @d = 1 and @q is in the neighborhood
 * of 0.5.
 *
 * The x,y pairs must be pre-sorted by increasing value of @x; an
 * error is flagged if this is not the case.  See also
 * sort_pairs_by_x().
 *
 * Returns: allocated vector containing the loess fitted values, or
 * %NULL on failure.
 */

gretl_matrix *loess_fit (const gretl_matrix *x, const gretl_matrix *y,
			 int d, double q, gretlopt opt, int *err)
{
    struct loess_info lo;
    gretl_matrix_block *B;
    gretl_matrix *Xi, *yi, *wt, *b;
    gretl_matrix *yh = NULL;
    gretl_matrix *rw = NULL;
    int N = gretl_vector_get_length(y);
    int k, iters, Xic, amin = 0;
    int n_ok, robust = 0, loo = 0;
    int i, n;

    if (d < 0 || d > 2 || q > 1.0) {
	*err = E_DATA;
	return NULL;
    }

    if (!data_pre_sorted(x)) {
	gretl_errmsg_set(_("loess: the data must be sorted by x"));
	*err = E_DATA;
	return NULL;
    }

    /* check for usable data points */
    n_ok = loess_count_usable_obs(y, N, &amin);
    if (n_ok < 4) {
	*err = E_TOOFEW;
	return NULL;
    }

    if (opt & OPT_O) {
	/* leave one out: experimental */
	loo = 1;
    }

    /* check for q too small */
    if (q < (d + 1.0) / n_ok) {
	q = (d + 1.0) / n_ok;
    }

    /* set the local sub-sample size */
    n = (int) ceil(q * n_ok);

    /* we need a minimum of two columns in Xi, in order
       to compute the x-distance based weights */
    Xic = (d == 0)? 2 : d + 1;

    B = gretl_matrix_block_new(&Xi, n, Xic,
			       &yi, n, 1,
			       &wt, n, 1,
			       &b, d+1, 1,
			       NULL);
    if (B == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* vector to hold the fitted values */
    yh = gretl_column_vector_alloc(N);
    if (yh == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    if (opt & OPT_R) {
	/* extra storage for residuals/robustness weights */
	rw = gretl_column_vector_alloc(N);
	if (rw == NULL) {
	    *err = E_ALLOC;
	    goto bailout;
	}
	iters = 2; /* maybe 3? */
	robust = 1;
    } else {
	iters = 1;
    }

    /* fill out the convenience struct */
    lo.y = y;
    lo.x = x;
    lo.Xi = Xi;
    lo.yi = yi;
    lo.wt = wt;
    lo.d = d;
    lo.n = n;
    lo.N = N;
    lo.loo = loo;

    for (k=0; k<iters && !*err; k++) {
	/* iterations for robustness, if wanted */
	int a = amin;

	lo.n_ok = n_ok;

	for (i=0; i<N && !*err; i++) {
	    /* iterate across points in full sample */
	    double xi = x->val[i];
	    int xconst = 0;

	    /* a is the leftmost possible index for starting to
	       read the n nearest neighbors of xi.
	    */
	    *err = loess_get_local_data(i, &a, &lo, &xconst);
	    if (*err) {
		break;
	    }

	    if (k > 0) {
		/* We have robustness weights, rw, based on the residuals
		   from the last round, which should be used to adjust
		   the wt as computed in loess_get_local_data().
		*/
		adjust_weights(rw, wt, a);
	    }

	    /* apply weights to the local data */
	    weight_local_data(&lo);

	    if (d == 0 || xconst) {
		/* not using x in the regressions: mask the x column(s) */
		gretl_matrix_reuse(Xi, -1, 1);
		if (b->rows > 1) {
		    gretl_matrix_reuse(b, 1, 1);
		}
	    }

	    /* run local WLS */
	    *err = gretl_matrix_SVD_ols(yi, Xi, b, NULL, NULL, NULL);

	    if (!*err) {
		/* yh: evaluate the polynomial at xi */
		yh->val[i] = b->val[0];
		if (b->rows > 1) {
		    yh->val[i] += b->val[1] * xi;
		    if (b->rows == 3) {
			yh->val[i] += b->val[2] * xi * xi;
		    }
		}
		if (robust && k < iters - 1) {
		    /* save residual for robustness weights */
		    if (na(y->val[i])) {
			rw->val[i] = NADBL;
		    } else {
			rw->val[i] = y->val[i] - yh->val[i];
		    }
		}
	    }

	    /* ensure matrices are at full size */
	    gretl_matrix_reuse(Xi, -1, Xic);
	    gretl_matrix_reuse(b, d+1, -1);

	} /* end loop over sample */

	if (!*err && robust && k < iters - 1) {
	    *err = make_robustness_weights(rw, N);
	}
    } /* end robustness iterations */

 bailout:

    gretl_matrix_block_destroy(B);
    gretl_matrix_free(rw);

    if (*err) {
	gretl_matrix_free(yh);
	yh = NULL;
    }

    return yh;
}

/* Silverman's rule-of-thumb bandwidth formula for kernel
   density estimation. See his Density Estimation for
   Statistics and Data Analysis (Chapman and Hall, 1986),
   pp. 45-48, and also Davidson and MacKinnon, Economic
   Theory and Methods (OUP, 2004), pp. 680-681.
*/

double kernel_bandwidth (const double *x, int n)
{
    double n5 = pow((double) n, -0.20);
    double s, A, q1, q3, r;
    int err = 0;

    s = gretl_stddev(0, n-1, x);
    q1 = gretl_quantile(0, n-1, x, 0.25, OPT_Q, &err);
    q3 = gretl_quantile(0, n-1, x, 0.75, OPT_Q, &err);
    r = (q3 - q1) / 1.349;
    A = (r > 0 && r < s)? r : s;

#if 0
    fprintf(stderr, "kernel_bandwidth: s=%g, q1=%g, q3=%g, IQR=%g, w=%g\n",
	    s, q1, q3, q3 - q1, w);
#endif

    return 0.9 * A * n5;
}
