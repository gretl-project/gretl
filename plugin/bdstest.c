/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2021 Allin Cottrell and Riccardo "Jack" Lucchetti
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

/*  Apparatus for BDS nonlinearity test: implementation in C of the
    algorithm advocated by Ludwig Kanzler in his 1999 paper,
    "Very fast and correctly sized estimation of the BDS statistic",
    DOI:10.2139/ssrn.151669.
*/

#include "libgretl.h"
#include "version.h"

typedef guint32 bword;
typedef struct kinfo_ kinfo;

struct kinfo_ {
    int n;           /* number of observations */
    int maxdim;      /* max correlation dimension */
    double eps;      /* "closeness" criterion */
    double kfac;     /* Dechert's 'k' */
    double c0;       /* first correlation */
    int bits;        /* number of bits required */
    int nwords;      /* number of words required */
    bword **wmat;    /* words matrix */
    guint32 *colsum; /* column sums */
    guint32 *rowsum; /* row sums */
    guint8 *bitinfo; /* counts of bits set */
    double *c1;      /* initial statistics */
    double *dspace;  /* workspace, doubles */
    double *wsave;   /* bootstrap test stats */
    int *bcount;     /* bootstrap counts */
    double *rx;      /* resampled x for bootstrap */
    int *z;          /* resampling indices */
};

static void resample_array (kinfo *ki, const double *x)
{
    int i;

    /* generate n drawings from [0 .. n-1] */
    gretl_rand_int_minmax(ki->z, ki->n, 0, ki->n - 1);

    /* fill target with selected rows of @x */
    for (i=0; i<ki->n; i++) {
	ki->rx[i] = x[ki->z[i]];
    }
}

/* return the sum of the sequence @imin:@imax */

static guint32 sumseq (int l, int r)
{
    guint32 ret = 0;
    int i, imin, imax;

    imin = l < r ? l : r;
    imax = r > l ? r : l;

    for (i=imin; i<=imax; i++) {
	ret += i;
    }

    return ret;
}

/* write to @targ the powers of @base specified in @seq */

static void ivp (double *targ, double base, int *seq, int len)
{
    int i;

    for (i=0; i<len; i++) {
	targ[i] = pow(base, seq[i]);
    }
}

static int make_bitinfo (kinfo *ki)
{
    int pp, len = pow(2.0, ki->bits);
    double *pb;
    int i, j;

    pb = malloc(ki->bits * sizeof *pb);
    ki->bitinfo = malloc(len * sizeof *ki->bitinfo);

    if (pb == NULL || ki->bitinfo == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<ki->bits; i++) {
	pb[i] = pow(2.0, -(ki->bits-i-1));
    }
    for (i=0; i<len; i++) {
	ki->bitinfo[i] = 0;
	for (j=0; j<ki->bits; j++) {
	    pp = floor(i * pb[j]);
	    ki->bitinfo[i] += pp % 2;
	}
    }
    free(pb);

    return 0;
}

static void free_words_mat (bword **wordmat, int n)
{
    if (wordmat != NULL) {
	int i;

	for (i=0; i<n-1; i++) {
	    free(wordmat[i]);
	}
	free(wordmat);
    }
}

static int allocate_words_matrix (kinfo *ki, int n)
{
    int i, err = 0;

    ki->wmat = calloc(n-1, sizeof *ki->wmat);
    if (ki->wmat == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<n-1 && !err; i++) {
	    ki->wmat[i] = calloc(ki->nwords, sizeof **ki->wmat);
	    if (ki->wmat[i] == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    return err;
}

/* only required for bootstrap */

static void reset_words_matrix (kinfo *ki, int n)
{
    int i, j;

    for (i=0; i<n-1; i++) {
	for (j=0; j<ki->nwords; j++) {
	    ki->wmat[i][j] = 0;
	}
    }
}

static int make_words_matrix (const double *x, int n, kinfo *ki)
{
    bword **wmat = NULL;
    bword *wrow = NULL;
    guint8 *bitvec;
    int *p2b;
    int nw, bvlen;
    int i, j, k, jj;
    int err;

    if (ki->wmat == NULL) {
	/* allocate not already done */
	err = allocate_words_matrix(ki, n);
	if (err) {
	    return err;
	}
    } else {
	/* bootstrapping */
	reset_words_matrix(ki, n);
	x = ki->rx;
    }

    wmat = ki->wmat;

    bitvec = calloc(n-1, sizeof *bitvec);
    p2b = malloc(ki->bits * sizeof *p2b);
    wrow = calloc(ki->nwords, sizeof *wrow);

    if (bitvec == NULL || p2b == NULL || wrow == NULL) {
	free(bitvec);
	free(p2b);
	free(wrow);
	return E_ALLOC;
    }

    p2b[0] = 1;
    for (i=1; i<ki->bits; i++) {
	p2b[i] = 2 * p2b[i-1];
    }

    for (i=0; i<n-1; i++) {
	bvlen = n-1-i;
	k = 0;
	for (j=i+1; j<n; j++) {
	    bitvec[k] = fabs(x[j] - x[i]) <= ki->eps;
	    ki->rowsum[i] += bitvec[k];
	    ki->colsum[j] += bitvec[k];
	    k++;
	}
	nw = ceil((n-i) / (double) ki->bits);
	k = 0;
	for (j=0; j<nw; j++) {
	    wrow[j] = 0;
	    for (jj=0; jj<ki->bits && k<bvlen; jj++) {
		wrow[j] += bitvec[k++] * p2b[jj];
	    }
	    wmat[i][j] = wrow[j];
	}
    }

    free(p2b);
    free(wrow);
    free(bitvec);

    return err;
}

static int get_num_bits (int n, double *mem)
{
    const double a[] = {
	0.000016,
	0.000005,
	0.000032
    };
    double mi, minv = 1e6;
    double n2 = n * (double) n;
    int i, bits = 1;

    for (i=1; i<=52; i++) {
	mi = (a[0] * i + a[1]) * pow(2.0, i);
	mi += a[2] * n2 / i;
	if (i == 1) {
	    minv = mi;
	} else if (mi < minv) {
	    minv = mi;
	    bits = i;
	}
    }

    if (mem != NULL) {
	*mem = minv;
    }

    return bits;
}

static int compute_k_c1 (int n, int maxdim, kinfo *ki)
{
    guint32 *md = malloc(maxdim * sizeof *md);
    guint32 *den = malloc(maxdim * sizeof *den);
    double bs0, tmp, fs2, xn = n;
    guint32 rcs;
    int i, k, err = 0;

    if (md == NULL || den == NULL) {
	free(md); free(den);
	return E_ALLOC;
    }

    /* materials for c1 */
    md[0] = 0;
    for (i=maxdim-1; i<n-1; i++) {
	md[0] += ki->rowsum[i];
    }
    k = maxdim - 1;
    for (i=0; i<maxdim-1; i++) {
	md[k--] = ki->rowsum[i];
    }

    /* c1, first round */
    ki->c1[0] = md[0];
    for (i=1; i<maxdim; i++) {
	ki->c1[i] = ki->c1[i-1] + md[i];
    }
    bs0 = ki->c1[maxdim-1];

    /* the divisor vector */
    den[0] = 0;
    for (i=1; i<=n-maxdim; i++) {
	den[0] += i;
    }

    /* divide c1 by cumulated den */
    ki->c1[0] /= den[0];
    k = n - maxdim + 1;
    for (i=1; i<maxdim; i++) {
	den[i] = den[i-1] + k++;
	ki->c1[i] /= den[i];
    }

    /* reverse c1 */
    for (i=0; i<maxdim/2; i++) {
	tmp = ki->c1[i];
	ki->c1[i] = ki->c1[maxdim-i-1];
	ki->c1[maxdim-i-1] = tmp;
    }

    /* sum of squares */
    fs2 = 0;
    for (i=0; i<n; i++) {
	rcs = ki->rowsum[i] + ki->colsum[i];
	fs2 += rcs * rcs;
    }

    ki->c0 = ki->c1[0];
    ki->kfac = (fs2 + 2.0*n - 3*(2*bs0+n)) / (xn * (xn-1) * (xn-2));

    if (isnan(ki->c0) || isnan(ki->kfac)) {
	fprintf(stderr, "bdstest: got NaN for c0 and/or k\n");
	err = E_NAN;
    }

    free(md);
    free(den);

    return err;
}

static double kanzler_eps (const double *x, int n, double e)
{
    guint32 s1 = 0, s2 = 0, i1 = 0, i2 = 0;
    guint32 smax = sumseq(1, n-1);
    double *dist;
    double eps;
    int i, j, k, n1 = n-1;

    dist = calloc(smax, sizeof *dist);
    if (dist == NULL) {
	return NADBL;
    }

    for (i=1; i<n; i++) {
	s1 = sumseq(0, i-2);
	s2 = sumseq(1, i-1);
	i1 = (i-1) * n1 - s1;
	i2 = i * n1 - s2 - 1;
	k = i;
	for (j=i1; j<=i2; j++) {
	    dist[j] = fabs(x[k++] - x[i-1]);
	}
    }

    qsort(dist, smax, sizeof *dist, gretl_compare_doubles);
    i = (int) round(e * smax);
    eps = dist[i-1];
    free(dist);

    return eps;
}

static void destroy_kinfo (kinfo *ki)
{
    free(ki->colsum);
    free(ki->rowsum);
    free(ki->c1);
    free(ki->bitinfo);
    free(ki->dspace);
    free(ki->wsave);
    free(ki->bcount);
    free(ki->rx);
    free(ki->z);

    if (ki->wmat != NULL) {
	free_words_mat(ki->wmat, ki->n);
    }
}

static int kanzler_init (kinfo *ki, int n, int maxdim)
{
    int i, err = 0;

    ki->n = n;
    ki->maxdim = maxdim;
    ki->bits = get_num_bits(n, NULL);
    ki->nwords = ceil((n-1) / (double) ki->bits);
    ki->wmat = NULL;
    ki->colsum = NULL;
    ki->rowsum = NULL;
    ki->c1 = NULL;
    ki->bitinfo = NULL;
    ki->dspace = NULL;
    ki->wsave = NULL;
    ki->bcount = NULL;
    ki->rx = NULL;
    ki->z = NULL;

    if (ki->bits > 32) {
	gretl_errmsg_sprintf("bdstest: input is too big (requires %d bits)\n", ki->bits);
	err = E_DATA;
    } else {
	ki->colsum = malloc(n * sizeof *ki->colsum);
	ki->rowsum = malloc(n * sizeof *ki->rowsum);
	ki->c1 = malloc(maxdim * sizeof *ki->c1);
	if (ki->colsum == NULL ||
	    ki->rowsum == NULL ||
	    ki->c1 == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<n; i++) {
		ki->colsum[i] = 1;
		ki->rowsum[i] = 0;
	    }
	}
    }

    return err;
}

/* allocate workspace for use in the main function */

static int make_dspace (kinfo *ki, int md1,
			double **iv1, double **iv2,
			double **c, double **sigma,
			double **w)
{
    ki->dspace = malloc(5 * md1 * sizeof *ki->dspace);

    if (ki->dspace == NULL) {
	return E_ALLOC;
    } else {
	*iv1 = ki->dspace;
	*iv2 = *iv1 + md1;
	*c =   *iv2 + md1;
	*sigma = *c + md1;
	*w = *sigma + md1;
    }

    return 0;
}

/* apparatus for bootstrap */

static int bds_boot_round (kinfo *ki, const double *x,
			   double *w, int md1)
{
    int i;

    if (ki->wsave == NULL) {
	/* initialize */
	ki->rx = malloc(ki->n * sizeof *ki->rx);
	ki->wsave = malloc(md1 * sizeof *ki->wsave);
	ki->bcount = malloc(md1 * sizeof *ki->bcount);
	ki->z = malloc(ki->n * sizeof *ki->z);

	if (ki->rx == NULL || ki->wsave == NULL ||
	    ki->bcount == NULL || ki->z == NULL) {
	    return E_ALLOC;
	}
	for (i=0; i<md1; i++) {
	    ki->wsave[i] = w[i];
	    ki->bcount[i] = 0;
	}
    } else {
	/* record */
	for (i=0; i<md1; i++) {
	    if (fabs(w[i]) > fabs(ki->wsave[i])) {
		ki->bcount[i] += 1;
	    }
	}
    }

    /* re-initialize counters */
    for (i=0; i<ki->n; i++) {
	ki->colsum[i] = 1;
	ki->rowsum[i] = 0;
    }

    /* resample */
    resample_array(ki, x);

    return 0;
}

/* Matrix with normalized statistics on first row,
   p-values on second row; the p-values are either
   asymptotic or bootstrapped.
*/

gretl_matrix *make_return_matrix (kinfo *ki, int NB,
				  double *w, int mdi)
{
    gretl_matrix *m = gretl_matrix_alloc(2, mdi);

    if (m != NULL) {
	double pv;
	int i;

	for (i=0; i<mdi; i++) {
	    if (ki->wsave != NULL) {
		gretl_matrix_set(m, 0, i, ki->wsave[i]);
		pv = ki->bcount[i] / (double) NB;
	    } else {
		gretl_matrix_set(m, 0, i, w[i]);
		pv = normal_pvalue_2(w[i]);
	    }
	    gretl_matrix_set(m, 1, i, pv);
	}
    }

    return m;
}

gretl_matrix *bdstest (const double *x, int n, int maxdim,
		       double eps, int *err)
{
    kinfo ki = {0};
    gretl_matrix *ret = NULL;
    bword *wcol = NULL;
    double *iv1 = NULL, *iv2 = NULL;
    double *c = NULL, *w = NULL;
    double *sigma = NULL;
    int *seq = NULL;
    int NB = 0, iter = 0;
    int i, j, m, md1;

    *err = kanzler_init(&ki, n, maxdim);
    if (*err) {
	return NULL;
    }

    if (eps < 0) {
	ki.eps = kanzler_eps(x, n, -eps);
    } else {
	ki.eps = eps * gretl_stddev(0, n-1, x);
    }
    if (na(ki.eps) || ki.eps <= 0) {
	*err = E_INVARG;
    }

    md1 = maxdim - 1;

    /* construct bitinfo vector */
    *err = make_bitinfo(&ki);

 run_again:

    /* construct words matrix */
    if (!*err) {
	*err = make_words_matrix(x, n, &ki);
    }

    /* allocate double-precision workspace */
    if (!*err && ki.dspace == NULL) {
	*err = make_dspace(&ki, md1, &iv1, &iv2, &c, &sigma, &w);
    }

    if (!*err) {
	*err = compute_k_c1(n, maxdim, &ki);
    }

    if (!*err && iter == 0) {
	seq = malloc(md1 * sizeof *seq);
	wcol = malloc((n-1) * sizeof *wcol);
	if (seq == NULL || wcol == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err) {
	goto bailout;
    }

    for (m=1; m<maxdim; m++) {
	double a1, a2a3, a4, a5;
	int bitcount = 0;

	for (j=0; j<ki.nwords; j++) {
	    for (i=0; i<n-1; i++) {
		wcol[i] = ki.wmat[i][j];
	    }
	    for (i=m; i<n-1; i++) {
		ki.wmat[i][j] = wcol[i] & wcol[i-1];
		bitcount += ki.bitinfo[ki.wmat[i][j]];
	    }
	}

	c[m-1] = bitcount / (double) sumseq(1, n-m-1);

	/* the a's hold the terms required to compute sigma */
	a1 = pow(ki.kfac, m+1);

	/* compute a2 * a3 */
	for (i=0; i<m; i++) {
	    seq[i] = m - i;
	}
	ivp(iv1, ki.kfac, seq, m);
	for (i=0; i<m; i++) {
	    seq[i] = 2*(i+1);
	}
	ivp(iv2, ki.c0, seq, m);
	a2a3 = 0;
	for (i=0; i<m; i++) {
	    a2a3 += 2 * iv1[i] * iv2[i];
	}

	a4 = pow(ki.c0, 2*(m+1));
	a5 = (m+1) * (m+1) * ki.kfac * pow(ki.c0, 2*m);

	sigma[m-1] = 2 * sqrt(a1 + a2a3 + m*m*a4 - a5);
    }

    /* compute normalized statistics */
    for (i=0; i<md1 && !*err; i++) {
	w[i] = sqrt(n-(i+2)+1) * (c[i] - pow(ki.c1[i+1], i+2)) / sigma[i];
	if (isnan(w[i])) {
	    fprintf(stderr, "bdstest: NaN computing w[%d]\n", i);
	    *err = E_NAN;
	}
    }

    if (!*err && iter < NB) {
	/* bootstrapping */
	*err = bds_boot_round(&ki, x, w, md1);
	if (!*err) {
	    iter++;
	    goto run_again;
	}
    }

    if (!*err) {
	ret = make_return_matrix(&ki, NB, w, md1);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

 bailout:

    destroy_kinfo(&ki);
    free(seq);
    free(wcol);

    return ret;
}
