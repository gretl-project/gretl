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
		if (i & bi->bits[j]) {
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
    gint16 *ip, *ip2;
    int j;

    for (j=0; j<n-dim; j++) {
	ip = bi->start[j];
	for (ip2 = bi->start[j+1]; ip2 <= bi->start[j+2]-1; ip2++) {
	    *ip = (*ip) & (*ip2);
	    ip++;
	}
	if (ip != bi->start[j+1]) {
	    *ip = 0;
	}
    }
}

/* mask pattern for row l, nbits: number of bits used
                           omit: number of bits omitted
                           mask: mask[0], mask[1] two-word mask
*/

static void genmask (int l, int n, int nbits, int omit,
		     int *mask, int *bits)
{
    int i, k, j, last, itrue;

    mask[0] = mask[1] = ALLBITS;
    last = (n-l-1)/nbits;
    for (i=n-omit; i<n; i++) {
	itrue = i - l - 1;
	j = itrue/nbits;
	k = nbits - 1 -(itrue % nbits);
	j = last - j;
	mask[j] = mask[j] ^ bits[k];
    }
}

/* count @c stats for grid: zero out parts not counted using mask;
   returns @c
*/

static double evalc (int n, bds_info *bi)
{
    gint32 nc = 0;
    double nd = n;
    gint16 *ip;
    int j;

    for (j=0; j<n; j++) {
	if (bi->start[j+1] - bi->start[j] > 2) {
	    for (ip = bi->start[j]; ip < bi->start[j+1]-2; ip++) {
		nc += bi->lookup[*ip];
	    }
	    for (ip = bi->start[j+1]-2; ip < bi->start[j+1]; ip++) {
		nc += bi->lookup[(*ip) & bi->mask[j*2 + bi->start[j+1]-ip-1]];
	    }
	} else {
	    for (ip = bi->start[j]; ip < bi->start[j+1]; ip++) {
		nc += bi->lookup[(*ip) & bi->mask[j*2 + bi->start[j+1]-ip-1]];
	    }
	}
    }

    return 2.0 * nc / (nd*(nd-1));
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
    bi->start[ix][ipos] |= bi->bits[ibit];
}

static double fkc (const double *x, int n, double *c,
		   int m, double eps, bds_info *bi)
{
    gint16 *ip;
    int i, nobs;
    struct position *pt;
    struct position *p;
    int remove = m - 1;
    gint32 count, tcount, phi;
    double dlen, k;

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
    count = phi = 0;
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
    count -= nobs;
    phi -= nobs + 3*count;
#if DEBUG
    printf("%d %d\n", count, phi);
#endif
    k = phi / (dlen * (dlen-1) * (dlen-2));
    c[1] = count / (dlen * (dlen-1));

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
    gretl_matrix *A = NULL;
    gretl_matrix *Ar = NULL;
    bds_info *bi = NULL;
    int *rej = NULL;
    int *bc = NULL;
    double *x = NULL;
    double *c = NULL;
    double *z = NULL;
    int i, j, n, m, r;
    int iters = 5000;
    int NB = 1999;
    int Ha = 0;
    double k, eps, pv;
    int err = 0;

    libgretl_init();

    if (argc == 1) {
	; /* OK, testing under the null */
    } else if (!strcmp(argv[1], "Ha")) {
	Ha = 1; /* testing under an alternative */
    } else {
	fprintf(stderr, "please give \"Ha\" or nothing\n");
	exit(EXIT_FAILURE);
    }

    n = 200;
    m = 5;
    eps = 0.50;
    A = gretl_matrix_alloc(n, 1);
    c = malloc((m+1) * sizeof *c);
    rej = malloc((m-1) * sizeof *rej);
    z = malloc((m-1) * sizeof *z);

    if (Ha) {
	x = malloc(n * sizeof *x);
    }

    /* initialize rejection counts */
    for (i=0; i<m-1; i++) {
	rej[i] = 0;
    }

    printf("monte carlo: iters=%d, n=%d, m=%d, eps=%g, NB=%d, Ha=%d\n",
	   iters, n, m, eps, NB, Ha);

    if (NB > 0) {
	bc = malloc((m-1) * sizeof *bc);
	Ar = gretl_matrix_alloc(n, 1);
    }

    bi = bds_info_new(n);
    if (bi == NULL) {
	fputs("Out of memory\n", stderr);
	exit(EXIT_FAILURE);
    }

    for (r=0; r<iters && !err; r++) {
	/* Monte Carlo iterations */
	if (NB > 0 && r > 0 && r % 100 == 0) {
	    fprintf(stderr, "mc iter %d\n", r);
	}
	gretl_matrix_random_fill(A, D_NORMAL);
	if (Ha) {
	    /* nonlinear MA alternative */
	    for (i=0; i<n; i++) {
		x[i] = A->val[i];
	    }
	    for (i=2; i<n; i++) {
		A->val[i] = x[i] + 0.8*x[i-1]*x[i-2];
	    }
	}

	/* calculate raw c and k statistics */
	k = fkc(A->val, n, c, m, eps, bi);
#if DEBUG
	printf("k = %lf\n",k);
	for (i=1; i<=m; i++) {
	    printf("c(%d) %lf\n", i, c[i]);
	}
#endif
	/* calculate normalized stats */
	for (i=2; i<=m; i++) {
	    z[i-2] = cstat(c[1], c[i], k, i, n-m+1);
	    if (NB == 0) {
		/* not doing bootstrap */
		pv = normal_pvalue_2(z[i-2]);
		if (iters <= 10) {
		    printf("dim %d: %f [%.4f]\n", i, z[i-2], pv);
		}
		if (pv < 0.05) {
		    rej[i-2] += 1;
		}
	    }
	}

	if (NB > 0) {
	    /* doing bootstrap */
	    double zji;

	    for (i=0; i<m-1; i++) {
		bc[i] = 0;
	    }
	    for (j=0; j<NB; j++) {
		err = gretl_matrix_resample2(Ar, A);
		k = fkc(Ar->val, n, c, m, eps, bi);
		for (i=2; i<=m; i++) {
		    zji = cstat(c[1], c[i], k, i, n-m+1);
		    if (fabs(zji) > fabs(z[i-2])) {
			bc[i-2] += 1;
		    }
		}
	    }
	    for (i=0; i<m-1 && !err; i++) {
		pv = bc[i] / (double) (NB + 1);
		if (pv < 0.05) {
		    rej[i] += 1;
		}
	    }
	} else if (iters <= 10) {
	    putchar('\n');
	}
    }

    if (err) {
	printf("Hit an error!\n");
    } else {
	printf("Rejection rates at nominal 5 percent level:\n");
	for (i=0; i<m-1; i++) {
	    printf("dim %d: %.3f\n", i+2, rej[i] / (double) iters);
	}
    }

    bds_info_destroy(bi);
    gretl_matrix_free(A);
    gretl_matrix_free(Ar);
    free(x);
    free(c);
    free(z);
    free(rej);
    free(bc);

    libgretl_cleanup();

    return 0;
}

#else /* not STANDALONE */

gretl_matrix *bdstest (const double *x, int n, int m, double eps, int *err)
{
    gretl_matrix *ret;
    gretl_matrix A = {0};
    gretl_matrix *Ar = NULL;
    bds_info *bi;
    double *c, *z;
    double k, pv;
    int *bc = NULL;
    int NB = 0;
    int i, j;

    if (n < 1500) {
	NB = 1999;
    }

    if (NB > 0) {
	/* apparatus for bootstrap */
	bc = malloc((m-1) * sizeof *bc);
	Ar = gretl_matrix_alloc(n, 1);
	if (bc == NULL || Ar == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	}
	gretl_matrix_init(&A);
	A.rows = n;
	A.cols = 1;
	A.val = (double *) x;
    }

    /* common allocations */
    c = malloc((m+1) * sizeof *c);
    z = malloc((m-1) * sizeof *z);
    ret = gretl_matrix_alloc(m-1, 2);
    bi = bds_info_new(n);

    if (c == NULL || z == NULL || ret == NULL || bi == NULL) {
	*err = E_ALLOC;
	gretl_matrix_free(ret);
	ret = NULL;
    } else {
	/* calculate raw c and k statistics */
	k = fkc(x, n, c, m, eps, bi);
	/* calculate normalized statistics */
	for (i=2; i<=m; i++) {
	    z[i-2] = cstat(c[1], c[i], k, i, n-m+1);
	    gretl_matrix_set(ret, i-2, 0, z[i-2]);
	    if (NB == 0) {
		/* record asymptotic p-value */
		pv = normal_pvalue_2(z[i-2]);
		gretl_matrix_set(ret, i-2, 1, pv);
	    }
	}
	if (NB > 0) {
	    double zji;

	    for (i=0; i<m-1; i++) {
		bc[i] = 0;
	    }
	    for (j=0; j<NB; j++) {
		gretl_matrix_resample2(Ar, &A);
		k = fkc(Ar->val, n, c, m, eps, bi);
		for (i=2; i<=m; i++) {
		    zji = cstat(c[1], c[i], k, i, n-m+1);
		    if (fabs(zji) >= fabs(z[i-2])) {
			bc[i-2] += 1;
		    }
		}
	    }
	    for (i=0; i<m-1; i++) {
		pv = bc[i] / (double) (NB + 1);
		gretl_matrix_set(ret, i, 1, pv);
	    }
	    gretl_matrix_free(Ar);
	    free(bc);
	}
    }

    bds_info_destroy(bi);
    free(c);
    free(z);

    return ret;
}

#endif /* STANDALONE or not */
