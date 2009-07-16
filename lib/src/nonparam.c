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

/* non-parametric stats for gretl */

#include "libgretl.h"

/* Spearman, one-tailed alpha = .005, .025, .05 */

static double rhocrit[18][3] = { 
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
    
    x = rhocrit[n - 7];

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
    double sd, sr;
    int i, m, ties = 0;
    int err = 0;

    *rho = *zval = NADBL;

    err = rankcorr_get_rankings(x, y, n, &rx, &ry, &m, &ties);
    if (err) {
	return err;
    }

    if (ties == 0) {
	/* calculate rho and z-score, no ties */
	sd = 0.0;
	for (i=0; i<m; i++) { 
	    sd += (rx[i] - ry[i]) * (rx[i] - ry[i]);
	}
	sr = 1.0 - 6.0 * sd / (m * (m * m - 1));
	sd = sqrt(1.0 / (m - 1.0));

	*rho = sr;
	*zval = sr / sd;
    } else {
	/* use Pearson in case of ties */
	sr = gretl_corr(0, m - 1, rx, ry, NULL);
	*rho = sr;
    }

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

static void print_raw_and_ranked (int vx, int vy,
				  const double *x, const double *y, 
				  const double *rx, const double *ry,
				  const DATAINFO *pdinfo,
				  PRN *prn)
{
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int t, i = 0;

    obs_marker_init(pdinfo);
    pprintf(prn, "\n     %s ", _("Obs"));
    pprintf(prn, "%13s%13s%13s%13s\n\n", pdinfo->varname[vx], _("rank"),
	    pdinfo->varname[vy], _("rank"));

    for (t=0; t<T; t++) {
	print_obs_marker(t + pdinfo->t1, pdinfo, prn);
	if (!(na(x[t])) && !(na(y[t]))) {
	    gretl_printxn(x[t], 15, prn);
	    pprintf(prn, "%15g", rx[i]);
	    gretl_printxn(y[t], 15, prn);
	    pprintf(prn, "%15g", ry[i]);
	    i++;
	}
	pputc(prn, '\n');
    }
}

/**
 * spearman:
 * @list: list of (two) variables to process.
 * @Z: data array.
 * @pdinfo: information on the data set.
 * @opt: if includes %OPT_V, print both the "raw" and the ranked data.
 * @prn: gretl printing struct.
 *
 * Calculates and prints Spearman's rank correlation coefficient for the two
 * variables specified in the @list.
 * 
 * Returns: 0 on successful completion, 1 on error.
 */

int spearman (const int *list, const double **Z, const DATAINFO *pdinfo, 
	      gretlopt opt, PRN *prn)
{
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    double *rx = NULL, *ry = NULL;
    const double *x, *y;
    double rho, zval;
    int vx, vy, m;
    int err;

    if (list[0] != 2) {
	pputs(prn, _("This command requires two variables\n"));
	return 1;
    }

    vx = list[1];
    vy = list[2];

    x = Z[vx] + pdinfo->t1;
    y = Z[vy] + pdinfo->t1;

    err = real_spearman_rho(x, y, T, &rho, &zval,
			    (opt & OPT_V)? &rx : NULL,
			    (opt & OPT_V)? &ry : NULL,
			    &m);
    if (err) {
	return err;
    }

    pprintf(prn, _("\nFor the variables '%s' and '%s'\n"), pdinfo->varname[vx],
	    pdinfo->varname[vy]);

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
	print_raw_and_ranked(vx, vy, x, y, rx, ry, pdinfo, prn);
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

int kendall (const int *list, const double **Z, const DATAINFO *pdinfo, 
	     gretlopt opt, PRN *prn)
{
    struct xy_pair *xy = NULL;
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    const double *x, *y;
    double tau, z;
    int vx, vy;
    int t, nn = 0;
    int err = 0;

    if (list[0] != 2) {
	pputs(prn, _("This command requires two variables\n"));
	return 1;
    }

    vx = list[1];
    vy = list[2];

    x = Z[vx] + pdinfo->t1;
    y = Z[vy] + pdinfo->t1;

    /* count valid pairs */
    for (t=0; t<T; t++) {
	if (!na(x[t]) && !na(y[t])) {
	    nn++;
	}
    } 

    if (nn < 2) {
	return E_MISSDATA;
    }

    /* allocate */
    xy = malloc(nn * sizeof *xy);
    if (xy == NULL) {
	return E_ALLOC;
    }

    /* calculate */
    err = real_kendall_tau(x, y, T, xy, nn, &tau, &z);

    if (!err) {
	pprintf(prn, _("\nFor the variables '%s' and '%s'\n"), pdinfo->varname[vx],
		pdinfo->varname[vy]);
	pprintf(prn, "Kendall's tau = %.8f\n", tau);
	pputs(prn, _("Under the null hypothesis of no correlation:\n "));
	pprintf(prn, _("z-score = %g, with two-tailed p-value %.4f\n"), z,
		normal_pvalue_2(z));
    }

    if (opt & OPT_V) {
	double *rx = NULL, *ry = NULL;

	rankcorr_get_rankings(x, y, T, &rx, &ry, NULL, NULL);

	if (rx != NULL && ry != NULL) { 
	    print_raw_and_ranked(vx, vy, x, y, rx, ry, pdinfo, prn);
	    free(rx);
	    free(ry);
	}
    }

    pputc(prn, '\n');

    /* clean up */
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

double lockes_test (const double *x, int t1, int t2)
{
    struct xy_pair *uv = NULL;
    double *sx = NULL, *u = NULL, *v = NULL;
    double zj, z;
    int m = t2 - t1 + 1;
    int i, j, t;
    int err;

    err = locke_shuffle_init(x + t1, &m, &sx);
    if (err) {
	return NADBL;
    }

    m /= 2;
	
    u = malloc(m * sizeof *u);
    v = malloc(m * sizeof *v);
    uv = malloc(m * sizeof *uv);

    if (u == NULL || v == NULL || uv == NULL) {
	free(u);
	free(v);
	free(uv);
	free(sx);
	return NADBL;
    }

    z = 0.0;

    /* repeat the shuffling of the series NREPEAT times, since the
       test statistic is sensitive to the ordering under the null
    */

    for (j=0; j<NREPEAT; j++) {
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
	printf("z[%d] = %g\n", j, zj);
#endif
    }   

    z /= (double) NREPEAT;

#if LOCKE_DEBUG
    fprintf(stderr, "Kendall's tau: average z = %g\n", z);
#endif 

    free(u);
    free(v);
    free(uv);
    free(sx);

    return z;
}

/**
 * runs_test:
 * @v: ID number of the variable to process.
 * @Z: data matrix.
 * @pdinfo: information on the data set.
 * @opt: %OPT_D to use first difference of variable, %OPT_E
 * to assume positive and negative are equiprobable.
 * @prn: gretl printing struct.
 *
 * Performs, and prints the results of, the runs test for randomness
 * for the variable specified by @varno.  The normal approximation
 * is that given in Gary Smith, Statistical Reasoning, 2e, p. 674.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 */

int runs_test (int v, const double **Z, const DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn)
{
    double xt, *x, mu, s2, sigma;
    double N2, z, pval;
    int Np, Nm;
    int t, n, runs = 1;

    n = pdinfo->t2 - pdinfo->t1 + 1;

    x = malloc(n * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    n = 0;

    if (opt & OPT_D) {
	double xt1;

	for (t=pdinfo->t1 + 1; t<=pdinfo->t2; t++) {
	    xt = Z[v][t];
	    xt1 = Z[v][t-1];
	    if (na(xt) || na(xt1)) {
		continue;
	    }
	    x[n++] = xt - xt1;
	} 	
    } else {
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    xt = Z[v][t];
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
	    pdinfo->varname[v], runs);

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

    record_test_result(z, pval, "runs");
  
    free(x);

    return 0;
}

static void print_z_prob (double z, PRN *prn)
{
    double p = NADBL;

    if (z > 0) {
	p = normal_pvalue_1(z);
	if (!na(p)) {
	    pprintf(prn, "  Prob(Z > %g) = %g\n", z, p);
	}
    } else if (z < 0) {
	p = normal_cdf(z);
	if (!na(p)) {
	    pprintf(prn, "  Prob(Z < %g) = %g\n", z, p);
	}
    }

    if (!na(p) && p > 0) {
	pprintf(prn, "  Two-tailed p-value = %g\n", 2 * p);
    }
}

static void diff_test_header (int v1, int v2, const DATAINFO *pdinfo,
			      PRN *prn)
{
    pputc(prn, '\n');
    pprintf(prn, _("Test for difference between %s and %s"),
	    pdinfo->varname[v1], pdinfo->varname[v2]);
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

static int 
signed_rank_test (const double *x, const double *y, 
		  int v1, int v2, const DATAINFO *pdinfo,
		  gretlopt opt, PRN *prn)
{
    struct ranker *r;
    double d, T, wp, wm;
    int Z = 0, N = 0;
    int i, k, t, n = 0;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (!na(x[t]) && !na(y[t])) {
	    if (x[t] != y[t]) {
		n++;
	    }
	    N++;
	}
    }

    if (n == 0) {
	return E_MISSDATA;
    }

    Z = N - n; /* number of zero-differences */

    r = malloc(n * sizeof *r);
    if (r == NULL) {
	return E_ALLOC;
    }

    i = 0;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
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

    diff_test_header(v1, v2, pdinfo, prn);
    pprintf(prn, "\n%s\n", _("Wilcoxon Signed-Rank Test"));
    pprintf(prn, "%s: %s\n\n", _("Null hypothesis"),
	    _("the median difference is zero"));

    if (opt & OPT_V) {
	pprintf(prn, "%16s %8s %16s\n\n", "difference", "rank", 
		"signed rank");
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

    pprintf(prn, "  n = %d\n", n);
    pprintf(prn, "  W+ = %g, W- = %g\n", wp, wm);
    pprintf(prn, "  (zero differences: %d, non-zero ties: %d)\n", Z, k);

    if (n > 8) {
	double s2, x, num, z;

	x = (N*(N+1) - Z*(Z+1)) / 4.0;
	pprintf(prn, "  %s = %g\n", _("Expected value"), x);
	s2 = (N*(N+1)*(2*N+1) - Z*(Z+1)*(2*Z+1) - T/2) / 24.0;
	pprintf(prn, "  %s = %g\n", _("Variance"), s2);
	num = wp - x;
	if (num > 0.25) {
	    num -= .5;
	} else {
	    num += .5;
	}
	z = num / sqrt(s2);
	pprintf(prn, "  z = %g\n", z);
	print_z_prob(z, prn);
    } else if (n > 5) {
	pprintf(prn, "  5%% critical values: %d (two-tailed), %d (one-tailed)\n",
		rank5[n-6][0], rank5[n-6][1]);
    } else {
	pprintf(prn, "  Sample too small for statistical significance\n");
    }

    pputc(prn, '\n');

    free(r);

    return 0;
}

static int rank_sum_test (const double *x, const double *y, 
			  int v1, int v2, const DATAINFO *pdinfo,
			  gretlopt opt, PRN *prn)
{
    struct ranker *r;
    double wa;
    char xc = 'a', yc = 'b';
    int na = 0, nb = 0;
    int i, t, n = 0;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
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
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
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

    diff_test_header(v1, v2, pdinfo, prn);
    pprintf(prn, "\n%s\n", _("Wilcoxon Rank-Sum Test"));
    pprintf(prn, "%s: %s\n\n", _("Null hypothesis"),
	    _("the two medians are equal"));

    wa = 0.0;

    if (opt & OPT_V) {
	pprintf(prn, "%10s %7s %8s\n\n", "value", "rank", "group");
    }

    for (i=0; i<n; i++) {
	if (opt & OPT_V) {
	    pprintf(prn, "%10g %7g %8c\n", r[i].val, r[i].rank, r[i].c);
	}
	if (r[i].c == 'a') {
	    wa += r[i].rank;
	} 
    }

    if (opt & OPT_V) {
	pputc(prn, '\n');
    }

    pprintf(prn, "  n1 = %d, n2 = %d\n", na, nb);
    pprintf(prn, "  w (%s) = %g\n", _("sum of ranks, sample 1"), wa);

    if (na >= 10 && nb >= 10) {
	double m, s, z;

	m = na * (na + nb + 1) / 2.0;
	s = sqrt(na * nb * (na + nb + 1) / 12.0);
	z = (wa - m) / s;
	pprintf(prn, "  z = (%g - %g) / %g = %g\n", wa, m, s, z);
	print_z_prob(z, prn);
    } else if (na >= 4 && nb >= 4 && nb <= 12) {
	void (*cv) (int, int, PRN *);
	void *handle;

	cv = get_plugin_function("rank_sum_lookup", &handle);
	if (cv != NULL) {
	    (*cv)(na, nb, prn);
	    close_plugin(handle);
	}
    } else {
	pprintf(prn, "  Sample too small for statistical significance\n");
    }

    pputc(prn, '\n');

    free(r);

    return 0;
}

static int sign_test (const double *x, const double *y, 
		      int v1, int v2, const DATAINFO *pdinfo,
		      gretlopt opt, PRN *prn)
{
    double pv;
    int n, w, t;

    n = w = 0;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (!na(x[t]) && !na(y[t]) && x[t] != y[t]) {
	    w += (x[t] > y[t]);
	    n++;
	}
    }

    if (n == 0) {
	return E_MISSDATA;
    }

    diff_test_header(v1, v2, pdinfo, prn);
    pprintf(prn, "\n%s\n\n", _("Sign Test"));
    pprintf(prn, _("Number of differences: n = %d\n"), n);
    pputs(prn, "  ");
    pprintf(prn, _("Number of cases with %s > %s: w = %d (%.2f%%)\n"), 
	    pdinfo->varname[v1], pdinfo->varname[v2],
	    w, 100.0 * w / n);

    pputs(prn, "  ");
    pprintf(prn, _("Under the null hypothesis of no difference, W "
		   "follows B(%d, %.1f)\n"), n, 0.5);
    pprintf(prn, "  %s(W <= %d) = %g\n", _("Prob"), w, 
	    binomial_cdf(0.5, n, w));
    if (w == 0) {
	pv = 1.0;
    } else {
	pv = binomial_cdf_comp(0.5, n, w - 1);
    }
    pprintf(prn, "  %s(W >= %d) = %g\n\n", _("Prob"), w, pv);

    return 0;
}

/**
 * diff_test:
 * @list: list containing 2 variables.
 * @Z: data array.
 * @pdinfo: information on the data set.
 * @opt: %OPT_G, sign test; %OPT_R, rank sum; %OPT_I,
 * signed rank; %OPT_V, verbose (for rank tests).
 * @prn: gretl printing struct.
 *
 * Performs, and prints the results of, a non-parametric
 * test for a difference between two variables or groups.
 * The specific test performed depends on @opt.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 */

int diff_test (const int *list, const double **Z, const DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn)
{
    const double *x, *y;
    int v1, v2;

    if (list[0] != 2) {
	pputs(prn, _("This command requires two variables\n"));
	return 1;
    }

    v1 = list[1];
    v2 = list[2];

    x = Z[v1];
    y = Z[v2];

    if (opt == OPT_NONE || opt == OPT_V) {
	opt = OPT_G;
    }

    if (opt & OPT_G) {
	return sign_test(x, y, v1, v2, pdinfo, opt, prn);
    } else if (opt & OPT_R) {
	return rank_sum_test(x, y, v1, v2, pdinfo, opt, prn);
    } else if (opt & OPT_I) {
	return signed_rank_test(x, y, v1, v2, pdinfo, opt, prn);
    }

    return 1;
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
 * @order: location to recieve sort order, or %NULL.
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

/* revised weights based on bisquare function (Cleveland) */

static void robust_weights (gretl_matrix *ei, const gretl_matrix *X, 
			    const gretl_matrix *y, const gretl_matrix *bj, 
			    gretl_matrix *wt, const gretl_matrix *wi,
			    int d)
{
    int i, n = gretl_vector_get_length(wt);
    double s, es, xi, di;

    for (i=0; i<n; i++) {
	ei->val[i] = y->val[i] - bj->val[0];
	if (d > 0) {
	    xi = gretl_matrix_get(X, i, 1);
	    ei->val[i] -= bj->val[1] * xi;
	    if (d == 2) {
		ei->val[i] -= bj->val[2] * xi * xi;
	    }
	}
	ei->val[i] = fabs(ei->val[i]);
    }

    s = gretl_median(0, n - 1, ei->val);

    for (i=0; i<n; i++) {
	es = ei->val[i] / (6.0 * s);
	if (es < 1.0) {
	    di = (1.0 - es * es) * (1.0 - es * es);
	} else {
	    di = 0.0;
	}
	wi->val[i] = di * wt->val[i];
    }
}

static void weight_x_y (const gretl_matrix *x, const gretl_matrix *y,
			gretl_matrix *Xr, gretl_matrix *yr,
			const gretl_matrix *w, int j, int d)
{
    int t, n = gretl_vector_get_length(w);
    double xrt, wt;

    for (t=0; t<n; t++) {
	wt = sqrt(w->val[t]);
	gretl_matrix_set(Xr, t, 0, wt);
	if (d > 0) {
	    xrt = x->val[t+j];
	    gretl_matrix_set(Xr, t, 1, xrt * wt);
	    if (d == 2) {
		gretl_matrix_set(Xr, t, 1, xrt * xrt * wt);
	    }
	}		
	yr->val[t] = y->val[t+j] * wt;
    }
}

/**
 * loess_fit:
 * @x: x-axis variable (must be pre-sorted).
 * @y: response variable.
 * @d: order for polynomial fit (0 <= d <= 2).
 * @q: bandwidth (0 < q < 1).
 * @opt: give %OPT_R for robust variant (with re-weighting based on
 * the first-stage residuals).
 * @err: location to recieve error code.
 *
 * Computes loess estimates based on William Cleveland, "Robust Locally 
 * Weighted Regression and Smoothing Scatterplots", Journal of the 
 * American Statistical Association, Vol. 74 (1979), pp. 829-836.
 * Generally one expects that @d = 1 and @q is in the neighborhood
 * of 0.5.
 *
 * The x,y pairs must be pre-sorted by increasing value of x; an
 * error is flagged if this is not the case.  See also
 * sort_pairs_by_x().
 * 
 * Returns: allocated vector containing the loess fitted values, or
 * %NULL on failure (in which case @err will contain a non-zero
 * error code).
 */

gretl_matrix *loess_fit (const gretl_matrix *x, const gretl_matrix *y,
			 int d, double q, gretlopt opt, int *err)
{
    gretl_matrix *Xr = NULL;
    gretl_matrix *yr = NULL;
    gretl_matrix *wt = NULL;
    gretl_matrix *bj = NULL;
    gretl_matrix *yh = NULL;
    gretl_matrix *ei = NULL;
    gretl_matrix *wi = NULL;

    int T = gretl_vector_get_length(y);
    double xi, d1, d2, dmax;
    int jmin, jmax;
    int k, kmax;
    int i, j, t, n;

    if (d < 0 || d > 2) {
	*err = E_DATA;
	return NULL;
    }

    if (!data_pre_sorted(x)) {
	strcpy(gretl_errmsg, "loess: the data must be sorted by x");
	*err = E_DATA;
	return NULL;
    }

    /* check the supplied q value */
    if (q > 1.0) {
	q = 1;
    } else if (q < (d + 1.0) / T) {
	q = (d + 1.0) / T;
    }

    n = (int) ceil(q * T);

    Xr = gretl_matrix_alloc(n, d + 1);
    yr = gretl_column_vector_alloc(n);
    wt = gretl_column_vector_alloc(n);
    bj = gretl_column_vector_alloc(d + 1);
    yh = gretl_column_vector_alloc(T);

    if (Xr == NULL || yr == NULL || wt == NULL || 
	bj == NULL || yh == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    if (opt & OPT_R) {
	ei = gretl_column_vector_alloc(n);
	wi = gretl_column_vector_alloc(n);
	if (ei == NULL || wi == NULL) {
	    *err = E_ALLOC;
	    goto bailout;
	}
	kmax = 3;
    } else {
	kmax = 1;
    }

    for (i=0; i<T && !*err; i++) {
	double xdt, xds;

	xi = x->val[i];

	/* get the search bounds */
	jmin = i - n + 1;
	if (jmin < 0) {
	    jmin = 0;
	}
	jmax = i;
	if (jmax + n > T) {
	    jmax = T - n;
	}

	/* find the n x-values closest to xi */
	for (j=jmin; j<=jmax; j++) {
	    if (j + n >= T) {
		break;
	    }
	    d1 = fabs(xi - x->val[j]);
	    d2 = fabs(xi - x->val[j + n]);
	    if (d1 <= d2) {
		break;
	    }
	}

	dmax = 0.0;
	for (t=0; t<n; t++) {
	    xdt = fabs(xi - x->val[t+j]);
	    if (xdt > dmax) {
		dmax = xdt;
	    }
	}

	for (t=0; t<n; t++) {
	    /* compute scaled distances */
	    xdt = fabs(xi - x->val[t+j]);
	    xds = xdt / dmax;
	    /* form tricube weights */
	    if (xds >= 1.0) {
		wt->val[t] = 0.0;
	    } else {
		wt->val[t] = pow(1.0 - pow(xds, 3.0), 3.0);
	    }
	}

	/* form weighted vars */
	weight_x_y(x, y, Xr, yr, wt, j, d);

	for (k=0; k<kmax && !*err; k++) {
	    /* run local WLS */
	    *err = gretl_matrix_ols(yr, Xr, bj, NULL, NULL, NULL);

	    if (!*err && k < kmax - 1) {
		/* re-weight based on the previous residuals */
		for (t=0; t<n; t++) {
		    gretl_matrix_set(Xr, t, 1, x->val[t+j]);
		    gretl_vector_set(yr, t, y->val[t+j]);
		}
		robust_weights(ei, Xr, yr, bj, wt, wi, d);
		weight_x_y(x, y, Xr, yr, wi, j, d);
	    }
	}	

	if (!*err) {
	    yh->val[i] = bj->val[0];
	    if (d > 0) {
		yh->val[i] += bj->val[1] * x->val[i];
		if (d == 2) {
		    yh->val[i] += bj->val[2] * x->val[i] * x->val[i];
		}
	    }
	}
    }

 bailout:
	
    gretl_matrix_free(Xr);
    gretl_matrix_free(yr);
    gretl_matrix_free(wt);
    gretl_matrix_free(bj);
    gretl_matrix_free(ei);
    gretl_matrix_free(wi);

    if (*err) {
	gretl_matrix_free(yh);
	yh = NULL;
    } 

    return yh;
}

