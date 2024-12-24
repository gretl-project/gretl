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

/*  gretl_normal.c - "advanced" routines relating to the normal
    distribution
*/  

#include "libgretl.h"
#include "libset.h"
#include "../../cephes/libprob.h"

#if defined(_OPENMP) && !defined(__APPLE__)
# include <omp.h>
#endif

/**
 * invmills:
 * @x: double-precision value.
 *
 * Adapted by putting together code from gsl and TDA (Univ. Bochum). 
 * The latter is, in turn, based on 
 * A. V. Swan, The Reciprocal of Mills's Ratio, Algorithm AS 17,
 * Journal of the Royal Statistical Society. Series C (Applied Statistics), 
 * Vol. 18, No. 1 (1969), 115-116.
 *
 * Returns: the inverse Mills ratio, that is the ratio between the
 * normal density function and the complement of the distribution 
 * function, both evaluated at @x.
 */

#define SQRT_HALF_PI 1.2533141373155002512078826424
#define MILLS_BOTTOM -22.9
#define MILLS_TOP 25
#define MILLS_EPS 1.0e-09

double invmills (double x)
{
    double a, a0, a1, a2;
    double b, b0, b1, b2;
    double r, s, t, d;
    double ret;

    if (x == 0.0) {
        return 1.0 / SQRT_HALF_PI;
    }

    if (x < MILLS_BOTTOM) {
	return 0;
    }

    if (x > MILLS_TOP) {
	a0 = 1.0/(x * x);
	a1 = 1.0 - 9.0 * a0 * (1.0 - 11.0 * a0);
	a2 = 1.0 - 5.0 * a0 * (1.0 - 7.0 * a0 * a1);
	d = 1.0 - a0 * (1.0 - 3.0 * a0 * a2);
	return x / d;
    }

    d = (x < 0.0)? -1.0 : 1.0;
    x = fabs(x);

    if (x <= 2.0) {
	s = 0.0;
	a = 1.0;
	r = t = x;
	b = x * x;
	while (fabs(s-t) > MILLS_EPS) {
	    a += 2.0;
	    s = t;
	    r *= b / a;
	    t += r;
	}
	ret = 1.0 / (SQRT_HALF_PI * exp(0.5 * b) - d * t);
    } else {
	a = 2.0;
	r = s = b1 = x;
	a1 = x * x + 1.0;
	a2 = x * (a1 + 2.0);
	b2 = a1 + 1.0;
	t  = a2 / b2;
	while (fabs(r-t) > MILLS_EPS && fabs(s-t) > MILLS_EPS) {
	    a += 1.0;
	    a0 = a1;
	    a1 = a2;
	    a2 = x * a1 + a * a0;
	    b0 = b1;
	    b1 = b2;
	    b2 = x * b1 + a * b0;
	    r  = s;
	    s  = t;
	    t  = a2 / b2;
	}
	ret = t;
	if (d < 0.0) {
	    ret /= (2.0 * SQRT_HALF_PI * exp(0.5 * x * x) * t - 1.0);
	}
    }

    return ret;
}

#define GENZ_BVN 1

#if GENZ_BVN

/**
 * genz04:
 * @rho: correlation coefficient.
 * @limx: abscissa value, first Gaussian r.v.
 * @limy: abscissa value, second Gaussian r.v.
 *
 * Based on FORTRAN code by Alan Genz, with minor adaptations.
 * Original source at 
 * http://www.math.wsu.edu/faculty/genz/software/fort77/tvpack.f
 * No apparent license.
 *
 * The algorithm is from Drezner and Wesolowsky (1989), 'On the
 * Computation of the Bivariate Normal Integral', Journal of
 * Statistical Computation and Simulation, 35 pp. 101-107, with major
 * modifications for double precision, and for |rho| close to 1.
 *
 * Returns: for (x, y) a bivariate standard Normal rv with correlation
 * coefficient @rho, the joint probability that (x < @limx) and (y < @limy),
 * or #NADBL on failure.
 */

static double genz04 (double rho, double limx, double limy)
{
    double w[10], x[10];
    double absrho = fabs(rho);
    double h, k, hk, bvn, hs, asr;
    double a, b, as, d1, bs, c, d, tmp;
    double sn, xs, rs;
    int i, lg, j;

    if (absrho < 0.3) {
	w[0] = .1713244923791705;
	w[1] = .3607615730481384;
	w[2] = .4679139345726904;

	x[0] = -.9324695142031522;
	x[1] = -.6612093864662647;
	x[2] = -.238619186083197;

	lg = 3;
    } else if (absrho < 0.75) {
	w[0] = .04717533638651177;
	w[1] = .1069393259953183;
	w[2] = .1600783285433464;
	w[3] = .2031674267230659;
	w[4] = .2334925365383547;
	w[5] = .2491470458134029;

	x[0] = -.9815606342467191;
	x[1] = -.904117256370475;
	x[2] = -.769902674194305;
	x[3] = -.5873179542866171;
	x[4] = -.3678314989981802;
	x[5] = -.1252334085114692;

	lg = 6;
    } else {
	w[0] = .01761400713915212;
	w[1] = .04060142980038694;
	w[2] = .06267204833410906;
	w[3] = .08327674157670475;
	w[4] = .1019301198172404;
	w[5] = .1181945319615184;
	w[6] = .1316886384491766;
	w[7] = .1420961093183821;
	w[8] = .1491729864726037;
	w[9] = .1527533871307259;

	x[0] = -.9931285991850949;
	x[1] = -.9639719272779138;
	x[2] = -.9122344282513259;
	x[3] = -.8391169718222188;
	x[4] = -.7463319064601508;
	x[5] = -.636053680726515;
	x[6] = -.5108670019508271;
	x[7] = -.3737060887154196;
	x[8] = -.2277858511416451;
	x[9] = -.07652652113349733;

	lg = 10;
    }

    h = -limx;
    k = -limy;
    hk = h * k;
    bvn = 0.0;

    if (absrho < 0.925) {
	hs = (h * h + k * k) / 2;
	asr = asin(rho);
	for (i=0; i<lg; i++) {
	    for (j=0; j<=1; j++) {
		sn = sin(asr * (1 + (2*j-1)*x[i]) / 2);
		bvn += w[i] * exp((sn * hk - hs) / (1 - sn * sn));
	    }
	}
	bvn = bvn * asr / (2 * M_2PI); 

	d1 = -h;
	bvn += normal_cdf(d1) * normal_cdf(-k);
    } else {
	if (rho < 0.0) {
	    k = -k;
	    hk = -hk;
	}

	as = (1 - rho) * (1 + rho);
	a = sqrt(as);
	bs = (h - k) * (h - k);
	c = (4 - hk) / 8;
	d = (12 - hk) / 16;
	asr = -(bs / as + hk) / 2;
	if (asr > -100.0) {
	    bvn = a * exp(asr) * (1 - c * (bs - as) * (1 - d * bs / 5)
				  / 3 + c * d * as * as / 5);
	}

	if (-hk < 100.0) {
	    b = sqrt(bs);
	    d1 = -b / a;
	    /*
	      Note: the condition below was not in the original
	      FORTRAN code this was ripped off from. Without it, there
	      are a few problems for rho very near -1; the version
	      below seems to work ok. Could it be a problem with
	      normal_cdf in the left tail?
	    */
	    if (d1 > -12.0) {
		bvn -= exp(-hk / 2) * SQRT_2_PI * normal_cdf(d1) * b * 
		    (1 - c * bs * (1 - d * bs / 5) / 3);
	    }
	}

	a /= 2;

	for (i=0; i<lg; i++) {
	    for (j=0; j<=1; j++) {
		d1 = a * (1 + (2*j-1)*x[i]);
		xs = d1 * d1;
		rs = sqrt(1 - xs);
		asr = -(bs / xs + hk) / 2;
		if (asr > -100.0) {
		    tmp = exp(-hk * (1 - rs) / ((rs + 1) * 2)) / 
			rs - (c * xs * (d * xs + 1) + 1);
		    bvn += a * w[i] * exp(asr) * tmp;
		}
	    }
	}
	
	bvn = -bvn / M_2PI;

	if (rho > 0.0) {
	    bvn += normal_cdf((h > k) ? -h : -k);
	} else {
	    bvn = -bvn;
	    if (k > h) {
		bvn = bvn + normal_cdf(k) - normal_cdf(h);
	    }
	}
    }

    /* sanity check */
    return (bvn < 0) ? 0 : bvn;
}

#else

/**
 * drezner78:
 * @rho: correlation coefficient.
 * @a: abscissa value, first Gaussian r.v.
 * @b: abscissa value, second Gaussian r.v.
 *
 * Ripped and adapted from Gnumeric, with a bug corrected for the case
 * (a * b < 0) && (rho < 0).
 *
 * The algorithm is from Drezner (1978), 'Computation of the Bivariate
 * Normal Integral', Mathematics of Computation, volume 32, number 141.
 *
 * Returns: for (x, y) a bivariate standard Normal rv with correlation
 * coefficient @rho, the joint probability that (x < @a) and (y < @b), or
 * #NADBL on failure.
 */

static double drezner78 (double rho, double a, double b)
{
    static const double x[] = {0.24840615, 0.39233107, 0.21141819, 
			       0.03324666, 0.00082485334};
    static const double y[] = {0.10024215, 0.48281397, 1.0609498, 
			       1.7797294, 2.6697604};
    double ret = NADBL;
    double a1, b1, den;
    int i, j;

    den = sqrt(2.0 * (1 - rho * rho));

    a1 = a / den;
    b1 = b / den;

    if (a <= 0 && b <= 0 && rho < 0) {
	/* standard case */
	double sum = 0.0;

	for (i=0; i<5; i++) {
	    for (j=0; j<5; j++) {
		sum += x[i] * x[j] * 
		    exp (a1 * (2 * y[i] - a1) + 
			 b1 * (2 * y[j] - b1) + 
			 2 * rho * (y[i] - a1) * (y[j] - b1));
	    }
	}
	ret = (sqrt(1 - (rho * rho)) / M_PI * sum);
    } else if (a <= 0 && b >= 0 && rho > 0) {
	ret = normal_cdf(a) - bvnorm_cdf(-rho, a, -b);
    } else if (a >= 0 && b <= 0 && rho > 0) {
	ret = normal_cdf(b) - bvnorm_cdf(-rho, -a, b);
    } else if (a >= 0 && b >= 0 && rho < 0) {
	ret = normal_cdf(a) + normal_cdf(b) - 1 + bvnorm_cdf(rho, -a, -b);
    } else if ((a * b * rho) > 0) {
	int sgna = (a < 0)? -1 : 1;
	int sgnb = (b < 0)? -1 : 1;
	double rho1, rho2, tmp, delta;

	tmp = sqrt((a * a) - 2 * rho * a * b + (b * b));
	rho1 = (rho * a - b) * sgna / tmp;
	rho2 = (rho * b - a) * sgnb / tmp;
	delta = (sgna * sgnb && (rho > 0))? 0 : 0.5;

	ret = (bvnorm_cdf(rho1, a, 0) + bvnorm_cdf(rho2, b, 0) - delta);
    }    

    return ret;
}

#endif /* bvnorm variants */

/**
 * bvnorm_cdf:
 * @rho: correlation coefficient.
 * @a: abscissa value, first Gaussian r.v.
 * @b: abscissa value, second Gaussian r.v.
 *
 * Returns: for (x, y) a bivariate standard Normal rv with correlation
 * coefficient @rho, the joint probability that (x < @a) and (y < @b), or
 * #NADBL on failure.
 */

double bvnorm_cdf (double rho, double a, double b)
{
    if (fabs(rho) > 1) {
	return NADBL;
    }	

    if (rho == 0.0) {
	/* joint prob is just the product of the marginals */
	return normal_cdf(a) * normal_cdf(b);
    }

    if (rho == 1.0) {
	/* the two variables are in fact the same */
	return normal_cdf(a < b ? a : b);
    }
    
    if (rho == -1.0) {
	/* the two variables are perfectly negatively correlated: 
	   P(x<a, y<b) = P((x<a) && (x>b)) = P(x \in (b,a))
	*/
	return (a <= b) ? 0 : normal_cdf(a) - normal_cdf(b);
    }

#if GENZ_BVN
    return genz04(rho, a, b);
#else 
    return drezner78(rho, a, b);
#endif
}

/* next: GHK apparatus with various helper functions */

#define GHK_DEBUG 0

static int ghk_input_check (const gretl_matrix *C,
			    const gretl_matrix *A,
			    const gretl_matrix *B,
			    const gretl_matrix *U,
			    const gretl_matrix *dP)
{
    if (gretl_is_null_matrix(C) ||
	gretl_is_null_matrix(A) ||
	gretl_is_null_matrix(B) ||
	gretl_is_null_matrix(U)) {
	return E_DATA;
    }

    if (A->rows != B->rows ||
	A->cols != B->cols ||
	C->rows != A->cols ||
	C->cols != A->cols ||
	U->rows != A->cols) {
	return E_NONCONF;
    }

    if (dP != NULL) {
	int m = C->rows;
	int np = m + m + m*(m+1)/2;

	if (dP->rows != A->rows || dP->cols != np) {
	    return E_NONCONF;
	}
    }

    return 0;
}

static void vector_diff (gretl_matrix *targ,
			 const gretl_matrix *m1,
			 const gretl_matrix *m2)
{
    int i, n = gretl_vector_get_length(targ);

    for (i=0; i<n; i++) {
	targ->val[i] = m1->val[i] - m2->val[i];
    }
}

static void vector_subtract (gretl_matrix *targ,
			     const gretl_matrix *src)
{
    int i, n = gretl_vector_get_length(targ);

    for (i=0; i<n; i++) {
	targ->val[i] -= src->val[i];
    }
}

/*
  C  Lower triangular Cholesky factor of \Sigma, m x m
  A  Lower bound of rectangle, m x 1
  B  Upper bound of rectangle, m x 1
  U  Random variates, m x r
*/

static double GHK_1 (const gretl_matrix *C, 
		     const gretl_matrix *A, 
		     const gretl_matrix *B, 
		     const gretl_matrix *U,
		     gretl_matrix *TA,
		     gretl_matrix *TB,
		     gretl_matrix *WGT,
		     gretl_matrix *TT,
		     double huge)
{
    int m = C->rows; /* Dimension of the multivariate normal */
    int r = U->cols; /* Number of repetitions */
    double P, den = gretl_matrix_get(C, 0, 0);
    double ui, x, z, cjk, tki;
    int i, j, k;

    z = A->val[0];
    TA->val[0] = (z == -huge) ? 0 : ndtr(z / den);
    z = B->val[0];
    TB->val[0] = (z == huge) ? 1 : ndtr(z / den);
    
    for (i=1; i<r; i++) {
	TA->val[i] = TA->val[0];
	TB->val[i] = TB->val[0];
    }

    /* form WGT = TB - TA */
    vector_diff(WGT, TB, TA);
    gretl_matrix_zero(TT);

    for (i=0; i<r; i++) {
	ui = gretl_matrix_get(U, 0, i);
	x = TB->val[i] - ui * (TB->val[i] - TA->val[i]);
	gretl_matrix_set(TT, 0, i, ndtri(x));
    }

    for (j=1; j<m; j++) {
	den = gretl_matrix_get(C, j, j);

	for (i=0; i<r; i++) {
	    if (WGT->val[i] == 0) {
		/* If WGT[i] ever comes to be zero, it cannot in
		   principle be modified by the code below; in fact,
		   however, running through the computations
		   regardless may produce a NaN (since 0 * NaN = NaN).
		   This becomes more likely for large dimension @m.
		*/
		continue;
	    }
	    x = 0.0;
	    for (k=0; k<j; k++) {
		cjk = gretl_matrix_get(C, j, k);
		tki = gretl_matrix_get(TT, k, i);
		x += cjk * tki;
	    }

	    if (A->val[j] == -huge) {
		TA->val[i] = 0.0;
	    } else {
		TA->val[i] = ndtr((A->val[j] - x) / den);
	    }

	    if (B->val[j] == huge) {
		TB->val[i] = 1.0;
	    } else {
		TB->val[i] = ndtr((B->val[j] - x) / den);
	    }

	    /* component j draw */
	    ui = gretl_matrix_get(U, j, i);
	    x = TB->val[i] - ui * (TB->val[i] - TA->val[i]);
	    gretl_matrix_set(TT, j, i, ndtri(x));
	}

	/* accumulate weight */
	vector_subtract(TB, TA);
	for (i=0; i<r; i++) {
	    WGT->val[i] *= TB->val[i];
	}
    }

    P = 0.0;
    for (i=0; i<r; i++) {
	P += WGT->val[i];
    }
    P /= r;

    if (P < 0.0 || P > 1.0) {
	fprintf(stderr, "*** ghk error: P = %g\n", P);
	P = 0.0/0.0; /* force a NaN */
    }	

    return P;
}

/* Should we enable OMP for GHK calculations? And if so,
   what's the threshold problem size that makes use of
   OMP worthwhile?
*/
#if defined(_OPENMP) && !defined(__APPLE__)
# define GHK_OMP 1
# define OMP_GHK_MIN 59
#endif

/**
 * gretl_GHK:
 * @C: Cholesky decomposition of covariance matrix, lower triangular,
 * m x m.
 * @A: Lower bounds, n x m; in case a lower bound is minus infinity
 * this should be represented as -1.0E10.
 * @B: Upper bounds, n x m; in case an upper bound is plus infinity
 * this should be represented as +1.0E10.
 * @U: Uniform random matrix, m x r.
 *
 * Computes the GHK (Geweke, Hajivassiliou, Keane) approximation to 
 * the multivariate normal distribution function for n observations 
 * on m variates, using r draws. 
 *
 * Returns: an n x 1 vector of probabilities.
 */

gretl_matrix *gretl_GHK (const gretl_matrix *C,
			 const gretl_matrix *A,
			 const gretl_matrix *B,
			 const gretl_matrix *U,
			 int *err)
{
#ifdef GHK_OMP
    unsigned sz = 0;
#endif
    gretl_matrix_block *Bk = NULL;
    gretl_matrix *P = NULL;
    gretl_matrix *Ai, *Bi;
    gretl_matrix *TA, *TB;
    gretl_matrix *WT, *TT;
    double huge;
    int m, n, r;
    int ierr, ghk_err = 0;
    int ABok, pzero;
    int i, j;

    *err = ghk_input_check(C, A, B, U, NULL);
    if (*err) {
	return NULL;
    }

    huge = libset_get_double(CONV_HUGE);

    m = C->rows;
    n = A->rows;
    r = U->cols;

    P = gretl_matrix_alloc(n, 1);
    if (P == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

#ifdef GHK_OMP
    if (n >= 2) {
	sz = n * m * r;
    }
#endif

    set_cephes_hush(1);

#ifdef GHK_OMP
#pragma omp parallel if (sz>OMP_GHK_MIN) private(i,j,Bk,Ai,Bi,TA,TB,WT,TT,ABok,pzero,ierr)
#endif
    {
	Bk = gretl_matrix_block_new(&Ai, m, 1,
				    &Bi, m, 1,
				    &TA, 1, r,
				    &TB, 1, r,
				    &WT, 1, r,
				    &TT, m, r,
				    NULL);
	if (Bk == NULL) {
	    ierr = E_ALLOC;
	    goto calc_end;
	} else {
	    ierr = 0;
	}

#ifdef GHK_OMP
#pragma omp for
#endif
	for (i=0; i<n; i++) {
	    ABok = 1; pzero = 0;

	    for (j=0; j<m && !ierr; j++) {
		Ai->val[j] = gretl_matrix_get(A, i, j);
		Bi->val[j] = gretl_matrix_get(B, i, j);
		ABok = !(isnan(Ai->val[j]) || isnan(Bi->val[j]));

		if (!ABok) {
		    /* If there are any NaNs in A or B, there's no
		       point in continuing
		    */
		    P->val[i] = 0.0/0.0; /* NaN */
		    break;
		} else if (Bi->val[j] < Ai->val[j]) {
		    gretl_errmsg_sprintf("ghk: inconsistent bounds: B[%d,%d] < A[%d,%d]",
					 i+1, j+1, i+1, j+1);
		    ierr = E_DATA;
		} else if (Bi->val[j] == Ai->val[j]) {
		    P->val[i] = 0.0;
		    pzero = 1;
		    break;
		}
	    }
	    if (!ierr && !pzero && ABok) {
		P->val[i] = GHK_1(C, Ai, Bi, U, TA, TB, WT, TT, huge);
	    }
	}

    calc_end:
	if (ierr) {
	    ghk_err = ierr;
	}

	gretl_matrix_block_destroy(Bk);
    } /* end (possibly) parallel section */

    set_cephes_hush(0);

    if (ghk_err) {
	*err = ghk_err;
	gretl_matrix_free(P);
	P = NULL;
    }

    return P;
}

/* below: revised version of GHK (plus score) */

static void scaled_convex_combo (gretl_matrix *targ, double w,
				 const gretl_matrix *m1,
				 const gretl_matrix *m2,
				 double scale)
{
    int i, n = gretl_vector_get_length(targ);
    double x1, x2;

    for (i=0; i<n; i++) {
	x1 = m1->val[i] * scale;
	x2 = m2->val[i] * scale;
	targ->val[i] = x2 - w * (x2 - x1);
    }
}

static void combo_2 (gretl_matrix *targ,
		     double w1, const gretl_matrix *m1,
		     double w2, const gretl_matrix *m2)
{
    int i, n = gretl_vector_get_length(targ);

    for (i=0; i<n; i++) {
	targ->val[i] = w1 * m1->val[i] + w2 * m2->val[i];
    }
}

static void vector_copy_mul (gretl_matrix *targ,
			     const gretl_matrix *src,
			     double w)
{
    int i, n = gretl_vector_get_length(targ);

    for (i=0; i<n; i++) {
	targ->val[i] = w * src->val[i];
    }
}

/* New-style GHK computation for observation t, iteration j,
   including the derivatives.
*/

static double ghk_tj (const gretl_matrix *C,
		      const gretl_matrix *a,
		      const gretl_matrix *b,
		      const double *u,
		      gretl_matrix *dWT,
		      gretl_matrix_block *Bk,
		      double huge,
		      int *err)
{
    double phi_min = 1.0e-300;
    gretl_matrix *dTA = NULL;
    gretl_matrix *dTB = NULL;
    gretl_matrix *dm = NULL;
    gretl_matrix *dx = NULL;
    gretl_matrix *tmp = NULL;
    gretl_matrix *cj = NULL;
    gretl_matrix *TT = NULL;
    gretl_matrix *dTT = NULL;
    double TA, TB, Tdiff, WT;
    double z, x, fx, den;
    int m = C->rows;
    int npar = m + m + m*(m+1)/2;
    int inicol, j, i;

    dTA = gretl_matrix_block_get_matrix(Bk, 0);
    dTB = gretl_matrix_block_get_matrix(Bk, 1);
    dm  = gretl_matrix_block_get_matrix(Bk, 2);
    dx  = gretl_matrix_block_get_matrix(Bk, 3);
    tmp = gretl_matrix_block_get_matrix(Bk, 4);
    cj  = gretl_matrix_block_get_matrix(Bk, 5);
    TT  = gretl_matrix_block_get_matrix(Bk, 6);
    dTT = gretl_matrix_block_get_matrix(Bk, 7);

    gretl_matrix_block_zero(Bk);
    den = C->val[0];
    gretl_matrix_reuse(dTT, npar, 1);

    if (a->val[0] == -huge) {
	TA = 0.0;
    } else {
	z = a->val[0] / den;
	TA = normal_cdf(z);
	x = normal_pdf(z) / den;
	dTA->val[0] = x;
	dTA->val[2*m] = -x/den;
    }

    if (b->val[0] == huge) {
	TB = 1.0;
    } else {
	z = b->val[0] / den;
	TB = normal_cdf(z);
	x = normal_pdf(z) / den;
	dTB->val[m] = x;
	dTB->val[2*m] = -x/den;
    }

    WT = TB - TA;
    x = TB - u[0] * WT;
    TT->val[0] = normal_cdf_inverse(x);

    fx = normal_pdf(TT->val[0]);
    scaled_convex_combo(dTT, u[0], dTA, dTB, 1/fx);
    vector_diff(dWT, dTB, dTA);

    /* first column of the gradient which refers to C */
    inicol = 2*m + 1;

    for (j=1; j<m; j++) {
	double mj = 0.0;
	int k = 0, flip = 0;

	/* the "flip" switch implements a numerical trick
	   that's needed to achieve acceptable precision when a[i]
	   is large: in that case, we flip the signs of a and b
	   so to exploit the greater accuracy of ndtr in the left-hand
	   tail than in the right-hand one.
	*/

	gretl_matrix_reuse(TT, j, 1);
	gretl_matrix_reuse(cj, 1, j);

	for (i=0; i<j; i++) {
	    cj->val[i] = gretl_matrix_get(C, j, i);
	    mj += cj->val[i] * TT->val[i];
	}

	gretl_matrix_zero(dx);
	gretl_matrix_zero(dm);

	for (i=inicol; i<=inicol+j-1; i++) {
	    dm->val[i] = TT->val[k++];
	}
	/* don't do threaded multiplication under an OMP thread */
	gretl_matrix_multiply_mod_single(cj, GRETL_MOD_NONE,
					 dTT, GRETL_MOD_TRANSPOSE,
					 tmp, GRETL_MOD_NONE);
	gretl_matrix_add_to(dm, tmp);

        den = gretl_matrix_get(C, j, j);

	x = (a->val[j] - mj) / den;
	if (x <= -huge) {
            TA = 0.0;
	    gretl_matrix_zero(dTA);
	} else {	    
	    if (x > 8.0) {
#if GHK_DEBUG > 1
		fprintf(stderr, "x=%.3f, flipping!\n", x);
#endif
		flip = 1;
		TA = normal_cdf(-x);
	    } else {
		TA = normal_cdf(x);
	    }
	    fx = normal_pdf(x);
	    dx->val[j] = 1;
	    vector_subtract(dx, dm);
	    dx->val[inicol+j] -= x;
	    vector_copy_mul(dTA, dx, fx/den);
	}
	
	gretl_matrix_zero(dx);

	x = (b->val[j] - mj) / den;
 	if (x >= huge) {
            TB = flip ? 0.0 : 1.0;
	    gretl_matrix_zero(dTB);
	} else {	    
	    TB = normal_cdf(flip ? -x : x);
	    fx = normal_pdf(x);
	    dx->val[m+j] = 1;
	    vector_subtract(dx, dm);
	    dx->val[inicol+j] -= x;
	    vector_copy_mul(dTB, dx, fx/den);
	}

	if (flip) {
	    Tdiff = TA - TB;
	    x = TA - u[j] * Tdiff;
	    TT->val[j] = -normal_cdf_inverse(x);
	} else {
	    Tdiff = TB - TA;
	    x = TB - u[j] * Tdiff;
	    TT->val[j] = normal_cdf_inverse(x);
	}

	if (na(TT->val[j])) {
#if GHK_DEBUG
	    fprintf(stderr, "TT is NA at j=%d (x=%g)\n", j, x);
	    fprintf(stderr, " (TA=%.16g, TB=%.16g, u[j]=%g)\n", TA, TB, u[j]);
#endif
	    fx = 0.0;
	} else {
	    fx = normal_pdf(TT->val[j]);
	}

	if (fx < phi_min) {
#if GHK_DEBUG
	    fprintf(stderr, "uh-oh, phi=%g\n", fx);
#endif
	    gretl_matrix_zero(tmp);	    
	} else {
	    scaled_convex_combo(tmp, u[j], dTA, dTB, 1/fx);
	}
	gretl_matrix_reuse(dTT, npar, j+1);
	gretl_matrix_inscribe_matrix(dTT, tmp, 0, j, GRETL_MOD_TRANSPOSE);
	vector_diff(tmp, dTB, dTA);
	combo_2(dWT, WT, tmp, Tdiff, dWT);

	if (WT > 0) {
	    WT *= Tdiff; /* accumulate weight */
	}

        inicol += j+1;
    }

    return WT;
}

/* This function translates column positions
   for the derivatives of elements of C from the
   "horizontal vech" into the "proper vech"
*/

static int *column_indices (int m)
{
    int i, j, k = 0;
    int *ndx;

    ndx = malloc((m * (m+1))/2 * sizeof *ndx);
    if (ndx == NULL) {
	return ndx;
    }

    for (i=0; i<m; i++) {
	for (j=0; j<=i; j++) {
	    ndx[k++] = j*(j+1)/2 + j*(m-j-1) + i;
	}
    }

    return ndx;
}

static int reorder_dP (gretl_matrix *dP, int m)
{
    int nc = dP->cols - 2*m;
    gretl_matrix *tmp;
    int *ndx;
    double xij;
    int i, j, k;

    ndx = column_indices(m);
    if (ndx == NULL) {
	return E_ALLOC;
    }

    tmp = gretl_matrix_alloc(dP->rows, nc);
    if (tmp == NULL) {
	free(ndx);
	return E_ALLOC;
    }

    gretl_matrix_extract_matrix(tmp, dP, 0, 2*m, GRETL_MOD_NONE);

    /* the first and last two elements of the
       vech never have to be moved */

    for (j=2; j<nc-2; j++) {
	if (ndx[j] != j) {
	    k = ndx[j] + 2*m;
	    for (i=0; i<dP->rows; i++) {
		xij = gretl_matrix_get(tmp, i, j);
		gretl_matrix_set(dP, i, k, xij);
	    }
	}
    }

    gretl_matrix_free(tmp);
    free(ndx);

    return 0;
}

/* workspace for ghk_tj */

static gretl_matrix_block *ghk_block_alloc (int m, int npar)
{
    gretl_matrix_block *B = NULL;
    gretl_matrix *M[8] = {NULL};

    B = gretl_matrix_block_new(&M[0], 1, npar, /* dTA */
			       &M[1], 1, npar, /* dTB */
			       &M[2], 1, npar, /* dm */
			       &M[3], 1, npar, /* dx */
			       &M[4], 1, npar, /* tmp */
			       &M[5], 1, m,    /* cj */
			       &M[6], m, 1,    /* TT */
			       &M[7], npar, m, /* dTT */
			       NULL);
    return B;
}

/* GHK including calculation of derivative */

gretl_matrix *gretl_GHK2 (const gretl_matrix *C,
			  const gretl_matrix *A,
			  const gretl_matrix *B,
			  const gretl_matrix *U,
			  gretl_matrix *dP,
			  int *err)
{
#ifdef GHK_OMP
    unsigned sz = 0;
#endif
    gretl_matrix_block *Bk;
    gretl_matrix_block *Bk2;
    gretl_matrix *a, *b;
    gretl_matrix *P = NULL;
    gretl_matrix *dpj = NULL;
    const double *uj;
    int r, n, m, npar;
    double huge;
    int t, i, j;

    if (gretl_is_null_matrix(dP)) {
	*err = E_DATA;
	return NULL;
    }

    *err = ghk_input_check(C, A, B, U, dP);
    if (*err) {
	return NULL;
    }

    r = U->cols;
    n = A->rows;
    m = C->rows;
    npar = m + m + m*(m+1)/2;

    P = gretl_zero_matrix_new(n, 1);
    if (P == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    gretl_matrix_zero(dP);

#ifdef GHK_OMP
    if (n >= 2) {
	sz = n * m * r;
    }
#endif

    huge = libset_get_double(CONV_HUGE);
    set_cephes_hush(1);

#ifdef GHK_OMP
#pragma omp parallel if (sz>OMP_GHK_MIN) private(i,j,t,a,b,uj,Bk,Bk2,dpj)
#endif
    {
	Bk = gretl_matrix_block_new(&a, 1, m,
				    &b, 1, m,
				    NULL);
	if (Bk == NULL) {
	    *err = E_ALLOC;
	}
	Bk2 = ghk_block_alloc(m, npar);
	if (Bk2 == NULL) {
	    *err = E_ALLOC;
	}
	dpj = gretl_matrix_alloc(1, npar);
	if (dpj == NULL) {
	    *err = E_ALLOC;
	}
	if (*err) {
	    goto calc_end;
	}

#ifdef GHK_OMP
#pragma omp for
#endif
	for (t=0; t<n; t++) {
	    /* loop across observations */
	    int err_t = 0;

	    for (i=0; i<m; i++) {
		/* transcribe and check bounds at current obs */
		a->val[i] = gretl_matrix_get(A, t, i);
		b->val[i] = gretl_matrix_get(B, t, i);
		if (isnan(a->val[i]) || isnan(b->val[i])) {
		    err_t = E_MISSDATA;
		    break;
		} else if (b->val[i] < a->val[i]) {
		    *err = err_t = E_DATA;
		    break;
		}
	    }

	    if (err_t == E_DATA) {
		gretl_errmsg_sprintf("ghk: inconsistent bounds: B[%d,%d] < A[%d,%d]",
				     t+1, i+1, t+1, i+1);
	    } else if (err_t == E_MISSDATA) {
		P->val[t] = 0.0/0.0; /* NaN */
		for (i=0; i<npar; i++) {
		    gretl_matrix_set(dP, t, i, 0.0/0.0);
		}
	    }

	    if (!err_t) {
		for (j=0; j<r && !*err; j++) {
		    /* Monte Carlo iterations, using successive columns of U */
		    uj = U->val + j * m;
		    P->val[t] += ghk_tj(C, a, b, uj, dpj, Bk2, huge, err);
		    gretl_matrix_inscribe_matrix(dP, dpj, t, 0, GRETL_MOD_CUMULATE);
		}
	    }
	}

    calc_end:
	gretl_matrix_block_destroy(Bk);
	gretl_matrix_block_destroy(Bk2);
	gretl_matrix_free(dpj);
    } /* end (possibly) parallel section */

    set_cephes_hush(0);

    if (*err) {
	gretl_matrix_free(P);
	P = NULL;
    } else {
	double x;

	for (t=0; t<n; t++) {
	    P->val[t] /= r;
	    for (i=0; i<npar; i++) {
		x = gretl_matrix_get(dP, t, i);
		gretl_matrix_set(dP, t, i, x / r);
	    }
	}
	if (m > 2) {
	    reorder_dP(dP, m);
	}
    }

    return P;
}
