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

/*  Apparatus for BDS nonlinearity test based on code from Blake
    LeBaron.

    Original notice:

    Blake LeBaron
    Department of Economics
    University Of Wisconsin-Madison
    July 1988
    March 1990

    This software is distributed with the understanding that it is free.
    Users may distribute it to anyone as long as they don't charge anything.
    Also, the author gives this out without any support or responsibility
    for errors of any kind. I hope that the distribution of this software
    will further enhance are understanding of these new measures of
    dependence.
*/

#define STANDALONE 0

#if STANDALONE 
# include <gretl/libgretl.h>
#else
# include "libgretl.h"
# include "version.h"
#endif

/* NBITS: the number of usable bits per word entry. Only 15 bits
   are used to keep the lookup table reasonably small.
*/
#define NBITS 15
#define ALLBITS 0xffff
#define TABLEN 32767

#define DEBUG 0

struct position {
    double value;
    int pos;
};

typedef struct bds_info_ bds_info;

struct bds_info_ {
    int bits[NBITS];
    int *mask;
    gint16 *grid;
    gint16 **start;
    int *lookup;
    struct position *postab;
    struct position *poslast;
};

static void bds_info_destroy (bds_info *bi)
{
    free(bi->mask);
    free(bi->lookup);
    free(bi->postab);
    free(bi->start);
    free(bi->grid);
    free(bi);
}

static bds_info *bds_info_new (int n)
{
    bds_info *bi = calloc(1, sizeof *bi);
    int err = 0;

    if (bi == NULL) {
	err = E_ALLOC;
    } else {
	int memsize = 0;
	int i, j;

	bi->mask = calloc(2*n, sizeof *bi->mask);
	bi->lookup = calloc(TABLEN+1, sizeof *bi->lookup);
	bi->postab = calloc(n, sizeof *bi->postab);
	bi->start = calloc(n+1, sizeof *bi->start);

	if (bi->mask == NULL || bi->lookup == NULL ||
	    bi->postab == NULL || bi->start == NULL) {
	    err = 1;
	    goto bi_finish;
	}

	/* determine the required grid size */
	for (i=0; i<=n; i++) {
	    memsize += (n-i)/NBITS + 1;
	}

	bi->grid = calloc(memsize, sizeof *bi->grid);
	if (bi->grid == NULL) {
	    err = E_ALLOC;
	    goto bi_finish;
	}

	bi->start[0] = bi->grid;
	for (i=1; i<=n; i++) {
	    bi->start[i] = bi->start[i-1] + (n-i)/NBITS + 1;
	}

	/* bit vector */
	bi->bits[0] = 1;
	for (i=1; i<15; i++) {
	    bi->bits[i] = (bi->bits[i-1] << 1);
	}

	/* table for bit counting */
	for (i=0; i<=TABLEN; i++){
	    bi->lookup[i] = 0;
	    for (j=0; j<NBITS; j++) {
		if ((i & bi->bits[j]) != 0) {
		    bi->lookup[i] += 1;
		}
	    }
	}
    }

 bi_finish:

    if (err && bi != NULL) {
	bds_info_destroy(bi);
	bi = NULL;
    }

    return bi;
}

/* qsort callback */

static int comp (const void *a, const void *b)
{
    const struct position *pa = a;
    const struct position *pb = b;

    return (pa->value > pb->value) - (pa->value < pb->value);
}

/* embed bitmap grid to next higher dimension

   g(i,j) = g(i,j) & g(i+1,j+1)
*/

static void embed (int n, int dim, bds_info *bi)
{
    gint16 *i, *i2;
    int j;

    for (j=0; j<n-dim; j++) {
	i = *(bi->start+j);
	for (i2 = *(bi->start+j+1); i2 <= *(bi->start+j+2)-1; i2++) {
	    *i = (*i) & (*i2);
	    i++;
	}
	if (i != *(bi->start+j+1)) {
	    *i = 0;
	}
    }
}

/* mask pattern for row l, nbits: number of bits used
                           omit: number of bits omitted
                           mask: mask[0],mask[1] two word mask
*/

static void genmask (int l, int n, int nbits, int omit,
		     int *mask, int *bits)
{
    int i, k, j, last, itrue;

    mask[0] = mask[1] = ALLBITS;
    last = (n-l-1)/nbits;
    for (i=n-omit; i<n; i++) {
	itrue = i - l -1;
	j = itrue/nbits;
	k = nbits-1-(itrue % nbits);
	j = last-j;
	mask[j] = mask[j] ^ bits[k];
    }
}

/* count @c stats for grid: zero out parts not counted using mask;
   returns @c
*/

static double evalc (int n, bds_info *bi)
{
    long int count = 0;
    double nd = n;
    gint16 *i;
    int j;

    for (j=0; j<n; j++) {
	if ((*(bi->start+j+1) - *(bi->start+j)) > 2) {
	    for (i = *(bi->start+j); i< *(bi->start+j+1)-2; i++) {
		count += bi->lookup[*i];
	    }
	    for (i = *(bi->start+j+1)-2; i< *(bi->start+j+1); i++) {
		count += bi->lookup[(*i) & bi->mask[j*2 + *(bi->start+j+1)-i-1]];
	    }
	} else {
	    for (i = *(bi->start+j); i<*(bi->start+j+1); i++) {
		count += bi->lookup[(*i) & bi->mask[j*2 + *(bi->start+j+1)-i-1]];
	    }
	}
    }

    return 2.0 * count / (nd*(nd-1));
}

static void gridon (int ix, int iy, bds_info *bi)
{
    int tmp, ipos, ibit;

    if (ix == iy) {
	return;
    }
    if (ix > iy) {
	tmp = ix;
	ix = iy;
	iy = tmp;
    }
    iy = iy - ix - 1;
    ipos = iy / NBITS;
    ibit = NBITS - 1 - (iy % NBITS);
    *(*(bi->start+ix)+ipos) |= bi->bits[ibit];
}

static double fkc (double *x, int n, double *c,
		   int m, double eps, bds_info *bi)
{
    gint16 *ip;
    int i, nobs;
    struct position *pt;
    struct position *p;
    int remove = m - 1;
    long count, tcount;
    double dlen, k;
    long phi;

    nobs = n - remove;
    dlen = nobs;

    /* clear grid */
    for (ip=bi->grid; ip<=bi->start[n]; ip++) {
	*ip = 0;
    }

    /* initialize table from data */
    for (i=0; i<n; i++){
	bi->postab[i].value = x[i];
	bi->postab[i].pos = i;
    }
    /* and perform Thieler sort */
    qsort(bi->postab, n, sizeof *bi->postab, comp);
    bi->poslast = bi->postab + n - 1;

    /* start row by row construction */
    count = 0;
    phi = 0;
    for (p=bi->postab; p<=bi->poslast; p++) {
	tcount = 0;
	pt = p;
	/* count to right */
	while (pt->value - p->value <= eps) {
	    gridon(p->pos, pt->pos, bi);
	    if (p->pos < nobs && pt->pos < nobs) {
		tcount++;
	    }
	    if (pt == bi->poslast) {
		break;
	    } else {
		pt++;
	    }
	}
	if (p != bi->postab){
	    /* count to left: this is not necessary for
	       building the grid, but is needed to get
	       the Dechert k
	    */
	    pt = p-1;
	    while (p->value - pt->value <= eps) {
		if (p->pos < nobs && pt->pos < nobs) {
		    tcount++;
		}
		if (pt == bi->postab) {
		    break;
		} else {
		    pt--;
		}
	    }
	}
	count += tcount;
	/* Dechert speed up k */
	phi += tcount * tcount;
    }

    /* adjust @k and @c to U statistic */
    count = count - nobs;
    phi = phi - nobs - 3*count;
#if DEBUG
    printf("%ld %ld\n", count, phi);
#endif
    k = ((double) phi) / (dlen * (dlen-1) * (dlen-2));
    c[1] = ((double) count) / (dlen * (dlen-1));

    /* build mask */
    for (i=0; i<nobs; i++) {
	genmask(i, n, NBITS, remove, bi->mask+2*i, bi->bits);
    }

    for (i=2; i<=m; i++) {
	embed(n, i, bi);
	c[i] = evalc(nobs, bi);
    }

    return k;
}

static double ipow (double x, int m)
{
    double y = 1;
    int j;

    for (j=0; j<m; j++) {
	y *= x;
    }

    return y;
}

/* This function calculates the asymptotic standard error from @c
   and @k. It then returns the test statistic which is asymptotically
   N(0,1). These formulas can be found in Brock, Hsieh, LeBaron,
   page 43.
*/

static double cstat (double c, double cm, double k, int m, int n)
{
    double std, sigma = 0;
    int j;

    for (j=1; j<=m-1; j++) {
	sigma += 2 * ipow(k, m-j) * ipow(c, 2*j);
    }
    sigma += ipow(k, m) + (m-1)*(m-1) * ipow(c, 2*m)
	-m * m * k * ipow(c, 2*m-2);
    sigma *= 4.0;

    std = sqrt(sigma/n);

    return (cm - ipow(c, m)) / std;
}

#if STANDALONE

/* stand-alone program for Monte Carlo usage */

int main (int argc, char **argv)
{
    gretl_matrix *a = NULL;
    bds_info *bi = NULL;
    int *rej = NULL;
    int i, n, m, r;
    double *x;
    int iters = 0;
    FILE *ifile;
    double k, *c = NULL;
    double eps, z, pv;

    libgretl_init();

    if (argc < 2) {
	fprintf(stderr, "give a number of iterations\n");
	exit(EXIT_FAILURE);
    }

    iters = atoi(argv[1]);
    printf("monte carlo: iters = %d\n", iters);

    n = 1000;
    m = 5;
    eps = 0.50;
    a = gretl_matrix_alloc(n, 1);
    rej = malloc((m+1) * sizeof *rej);
    c = malloc((m+1) * sizeof *c);
    x = a->val;
    for (i=0; i<=m; i++) {
	rej[i] = 0;
    }

    bi = bds_info_new(n);
    if (bi == NULL) {
	fputs("Out of memory\n", stderr);
	exit(EXIT_FAILURE);
    }

    for (r=0; r<iters; r++) {
	if (a != NULL) {
	    gretl_matrix_random_fill(a, D_NORMAL);
	}
	/* calculate raw c and k statistics: this is the hard part */
	k = fkc(x, n, c, m, eps, bi);
#if DEBUG
	printf("k = %lf\n",k);
	for (i=1; i<=m; i++) {
	    printf("c(%d) %lf\n", i, c[i]);
	}
#endif
	/* calculate normalized stats: this is the easy part */
	for (i=2; i<=m; i++) {
	    z = cstat(c[1], c[i], k, i, n-m+1);
	    pv = normal_pvalue_2(z);
	    if (iters <= 10) {
		printf("dim %d: %f [%.4f]\n", i, z, pv);
	    }
	    if (pv < 0.05) {
		rej[i] += 1;
	    }
	}
	if (iters <= 10) {
	    putchar('\n');
	}
    }

    if (rej != NULL) {
	printf("Rejection rates at nominal 5 percent level:\n");
	for (i=2; i<=m; i++) {
	    printf("dim %d: %.3f\n", i, rej[i] / (double) iters);
	}
    }

    bds_info_destroy(bi);
    gretl_matrix_free(a);
    free(c);
    free(rej);

    libgretl_cleanup();

    return 0;
}

#else /* not STANDALONE */

gretl_matrix *bdstest (double *x, int n, int m, double eps, int *err)
{
    gretl_matrix *ret = NULL;
    bds_info *bi = NULL;
    double *c;
    double k;
    int i;

    c = malloc((m+1) * sizeof *c);
    ret = gretl_matrix_alloc(m-1, 1);
    bi = bds_info_new(n);
    
    if (c == NULL || ret == NULL || bi == NULL) {
	*err = E_ALLOC;
	gretl_matrix_free(ret);
	ret = NULL;
    } else {
	/* calculate raw c and k statistics */
	k = fkc(x, n, c, m, eps, bi);

	/* calculate normalized statistics */
	for (i=2; i<=m; i++) {
	    ret->val[i-2] = cstat(c[1], c[i], k, i, n-m+1);
	}
    }

    bds_info_destroy(bi);
    free(c);

    return ret;
}

#endif /* STANDALONE or not */
