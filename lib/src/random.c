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

/* random.c for gretl */

#include "libgretl.h"

#include <time.h>
#include <glib.h>

static GRand *gretl_rand;
static unsigned int useed;

#undef OLD_NORMAL

unsigned int gretl_rand_get_seed (void)
{
    return useed;
}

/**
 * gretl_rand_init:
 *
 * Initialize gretl's PRNG, using the system time as seed.
 */

void gretl_rand_init (void)
{
    useed = time(NULL);
    gretl_rand = g_rand_new();
    gretl_rand_set_seed((guint32) useed);
}

/**
 * gretl_rand_set_seed:
 * @seed: the chosen seed value.
 *
 * Set a specific (and hence reproducible) seed for gretl's PRNG.
 * But if the value 0 is given for @seed, set the seed using
 * the system time (which is the default when libgretl is
 * initialized).
 */

void gretl_rand_set_seed (unsigned int seed)
{
    if (seed == 0) {
	useed = time(NULL);
    } else {
	useed = seed;
    }
    g_rand_set_seed(gretl_rand, useed);
}

#define gretl_one_uniform() (g_rand_double_range(gretl_rand, 0, 1))

/**
 * gretl_one_snormal:
 *
 */

double gretl_one_snormal (void) 
{
    double x, y, z;

 tryagain:
    x = gretl_one_uniform();
    y = gretl_one_uniform();
    z = sqrt(-2. * log(x));
    if (isnan(z) || isinf(z)) {
	goto tryagain;
    }

    return (z * cos(M_2PI * y));
}

#define MOD_POLAR 1

static void gretl_two_snormals (double *z1, double *z2) 
{
    double x, y, z;

#if MOD_POLAR
 tryagain:
    x = 2 * gretl_one_uniform() - 1;
    y = 2 * gretl_one_uniform() - 1;
    z = x*x + y*y;
    if (z >= 1) {
	goto tryagain;
    }

    z = sqrt(-2 * log(z)/z);
    *z1 = z * x;
    *z2 = z * y;
#else
 tryagain:
    x = gretl_one_uniform();
    y = gretl_one_uniform();
    z = sqrt(-2 * log(x));
    if (isnan(z) || isinf(z)) {
	goto tryagain;
    }

    y *= M_2PI;

    *z1 = z * cos(y);
    *z2 = z * sin(y);
#endif
}

/**
 * gretl_rand_uniform_minmax:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @min: lower bound of range.
 * @max: upper bound of range.
 *
 * Fill the selected subset of array @a with pseudo-random drawings
 * from the uniform distribution on @min to @max, using the Mersenne
 * Twister.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int gretl_rand_uniform_minmax (double *a, int t1, int t2,
			       double min, double max) 
{
    int t;

    if (na(min) && na(max)) {
	min = 0.0;
	max = 1.0;
    } else if (na(min) || na(max) || max <= min) {
	return E_INVARG;
    }

    for (t=t1; t<=t2; t++) {
	a[t] = g_rand_double_range(gretl_rand, min, max);
    }

    return 0;
}

/**
 * gretl_rand_uniform:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the uniform distribution on 0-1, using the Mersenne
 * Twister.
 */

void gretl_rand_uniform (double *a, int t1, int t2) 
{
    int t;

    for (t=t1; t<=t2; t++) {
	a[t] = gretl_one_uniform();
    }
}

/**
 * gretl_rand_normal:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the standard normal distribution, using the Mersenne Twister
 * for uniform input and the Box-Muller method for converting to the
 * normal distribution.
 */

void gretl_rand_normal (double *a, int t1, int t2) 
{
    int t;

#ifdef OLD_NORMAL
    for (t=t1; t<=t2; t++) {
	a[t] = gretl_one_snormal();
    }
#else
    double z1, z2;

    for (t=t1; t<=t2; t++) {
	gretl_two_snormals(&z1, &z2);
	a[t] = z1;
	if (t < t2) {
	    a[++t] = z2;
	}
    }
#endif
}

/**
 * gretl_rand_normal_full:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 * @mean: mean of the distribution
 * @sd: standard deviation
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the normal distribution with the given mean and standard
 * deviation, using the Mersenne Twister for uniform input and the 
 * Box-Muller method for converting to the normal distribution.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int
gretl_rand_normal_full (double *a, int t1, int t2,
			double mean, double sd) 
{
    int t;

    if (na(mean) && na(sd)) {
	mean = 0.0;
	sd = 1.0;
    } else if (na(mean) || na(sd) || sd <= 0.0) {
	return E_INVARG;
    }

    gretl_rand_normal(a, t1, t2);

    if (mean != 0.0 || sd != 1.0) {
	for (t=t1; t<=t2; t++) {
	    a[t] = mean + a[t] * sd;
	}
    }
    
    return 0;
}

/**
 * gretl_rand_gamma:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @shape: shape parameter.
 * @scale: scale parameter.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the specified gamma distribution.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_gamma (double *a, int t1, int t2, 
		      double shape, double scale) 
{
    double *U = NULL;
    double e = 2.718281828459045235;
    double delta, dinv = 0, d1 = 0;
    double u, v, x, y, u0 = 0;
    int k, i, t;

    if (shape <= 0 || scale <= 0) {
	return E_DATA;
    }

    k = shape;
    if (k > 0) {
	U = malloc(k * sizeof *U);
	if (U == NULL) {
	    return E_ALLOC;
	}
    }

    delta = shape - k;

    if (delta > 0) {
	d1 = delta - 1;
	dinv = 1 / delta;
	u0 = e / (e + delta);
    }

    /* 
       Case of shape < 1 from Kundu and Gupta, "A convenient way of
       generating gamma random variables using generalized exponential
       distribution", Computational Statistics and Data Analysis, 51
       (2007).  Case of shape >= 1 taken from the Wikipedia entry on
       the gamma distribution.
    */

    for (t=t1; t<=t2; t++) {
	a[t] = 0.0;
	if (shape < 1) {
	    double ex2;

	    while (1) {
		u = gretl_one_uniform();
		v = gretl_one_uniform();
		x = -2 * log(1 - pow(u, dinv));
		ex2 = exp(-x/2);
		u0 = pow(x, d1) * ex2;
		u0 /= pow(2.0, d1) * pow((1-ex2), d1);
		if (v <= u0) {
		    a[t] = x;
		    break;
		}
	    }
	} else {
	    for (i=0; i<k; i++) {
		U[i] = gretl_one_uniform();
		while (U[i] == 0.0) {
		    U[i] = gretl_one_uniform();
		}
	    }
	    if (delta > 0) {
		while (1) {
		    u = gretl_one_uniform();
		    v = gretl_one_uniform();
		    if (u <= u0) {
			x = pow(u, dinv);
			y = v * pow(x, d1);
		    } else {
			x = 1 - log(u);
			y = v * exp(-x);
		    }
		    if (y <= pow(x, d1) * exp(-x)) {
			a[t] = x;
			break;
		    }
		}
	    } 
	    for (i=0; i<k; i++) {
		a[t] -= log(U[i]);
	    }
	}
	a[t] *= scale;	
    }

    free(U);
	
    return 0;
}

/**
 * gretl_rand_chisq:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @v: degrees of freedom.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Chi-Squared distribution with @v degrees of freedom, 
 * using the gamma r.v. generator.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_chisq (double *a, int t1, int t2, int v) 
{
    return gretl_rand_gamma(a, t1, t2, 0.5*v, 2);
}

/**
 * gretl_rand_student:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @v: degrees of freedom.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Student t distribution with @v degrees of freedom, 
 * using the Mersenne Twister for uniform input and the Box-Muller 
 * method for converting to the normal distribution.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_student (double *a, int t1, int t2, int v) 
{
    double *X2 = NULL;
    int T = t2 - t1 + 1;
    int t;

    if (v < 1) {
	return E_INVARG;
    }

    X2 = malloc(T * sizeof *X2);
    if (X2 == NULL) {
	return E_ALLOC;
    }

    gretl_rand_normal(a, t1, t2);
    gretl_rand_chisq(X2, 0, T-1, v);

    for (t=0; t<T; t++) {
	a[t + t1] /= sqrt(X2[t] / v);
    }

    free(X2);

    return 0;
}

/**
 * gretl_rand_F:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @v1: numerator degrees of freedom.
 * @v2: denominator degrees of freedom.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the F distribution with @v1 and @v2 degrees of freedom, 
 * using the Mersenne Twister for uniform input and the Box-Muller 
 * method for converting to the normal distribution.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_F (double *a, int t1, int t2, int v1, int v2) 
{
    double *b = NULL;
    int T = t2 - t1 + 1;
    int s, t;

    if (v1 < 1 || v2 < 1) {
	return E_INVARG;
    }

    b = malloc(T * sizeof *b);
    if (b == NULL) {
	return E_ALLOC;
    }

    gretl_rand_chisq(a, t1, t2, v1);
    gretl_rand_chisq(b, 0, T-1, v2);

    for (t=0; t<T; t++) {
	s = t + t1;
	a[s] = (a[s]/v1) / (b[t]/v2);
    }

    free(b);

    return 0;
}

/**
 * gretl_rand_binomial:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @n: number of trials.
 * @p: success probability per trial.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the binomial distribution with parameters @n and @p.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_binomial (double *a, int t1, int t2, int n, double p) 
{
    double *b;
    int i, t;

    if (n < 1 || p <= 0 || p >= 1) {
	return E_INVARG;
    }

    b = malloc(n * sizeof *b);
    if (b == NULL) {
	return E_ALLOC;
    }

    for (t=t1; t<=t2; t++) {
	a[t] = 0.0;
	gretl_rand_uniform(b, 0, n - 1);
	for (i=0; i<n; i++) {
	    if (b[i] <= p) {
		a[t] += 1;
	    }
	}
    }

    free(b);

    return 0;
}

/* Poisson rv with mean m */

static double genpois (const double m)
{
    double x;

    if (m > 200) {
	x = (m + 0.5) + sqrt(m) * gretl_one_snormal();
	x = floor(x);
    } else {
	int y = 0;

	x = exp(m) * gretl_one_uniform();
	while (x > 1) {
	    y++;
	    x *= gretl_one_uniform();
	}
	x = (double) y;
    }

    return x;
}

/**
 * gretl_rand_poisson:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @m: mean (see below).
 * @vec: should be 1 if @m is an array, else 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Poisson distribution with a mean determined by 
 * @m, which can either be a pointer to a scalar, or an array
 * of length greater than or equal to @t2 + 1.  
 */

void gretl_rand_poisson (double *a, int t1, int t2, const double *m,
			 int vec) 
{
    int t;

    for (t=t1; t<=t2; t++) {
	a[t] = (vec)? genpois(m[t]) : genpois(*m);
    }
}

/* f(x; k, \lambda) = (k/\lambda) (x/\lambda)^{k-1} e^{-(x/\lambda)^k} */

/**
 * gretl_rand_weibull:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @shape: shape parameter > 0.
 * @scale: scale parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Weibull distribution with shape @k and scale @lambda.
 *
 * Returns: 0 on success, non-zero if either parameter is out of
 * bounds.
 */

int gretl_rand_weibull (double *a, int t1, int t2, double shape,
			double scale) 
{
    int err = 0;

    if (shape <= 0 || scale <= 0) {
	err = E_DATA;
    } else {
	double u, kinv = 1.0 / shape;
	int t;

	for (t=t1; t<=t2; t++) {
	    u = gretl_one_uniform();
	    while (u == 0.0) {
		u = gretl_one_uniform();
	    }
	    a[t] = scale * pow(-log(u), kinv);
	}
    }

    return err;
}

/**
 * gretl_rand_int_max:
 * @max: the maximum value (open)
 *
 * Returns: a pseudo-random unsigned int in the interval [0, @max)
 * using the Mersenne Twister.
 */

unsigned int gretl_rand_int_max (unsigned int max)
{
    return g_rand_int_range(gretl_rand, 0, max);
}

/**
 * gretl_rand_int:
 *
 * Returns: a pseudo-random unsigned int on the interval
 * [0,0xffffffff] using the Mersenne Twister.
 */

unsigned int gretl_rand_int (void)
{
    return g_rand_int(gretl_rand);
}

/**
 * gretl_rand_free:
 *
 * Free the gretl_rand structure (may be called at program exit).
 */

void gretl_rand_free (void)
{
    g_rand_free(gretl_rand);
}


