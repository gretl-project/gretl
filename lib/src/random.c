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

/* random.c for gretl: RNGs */

#include "libgretl.h"
#include <time.h>
#include <fcntl.h>

#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

/* For optimizing the Ziggurat */
#if defined(i386) || defined (__i386__)
# define HAVE_X86_32 1
# define UMASK 0x80000000UL /* most significant w-r bits */
#else
# define HAVE_X86_32 0
#endif

#include "../../rng/splitmix64.c"
#include "../../rng/xoshiro256plus.c"

/**
 * SECTION:random
 * @short_description: generate pseudo-random values
 * @title: PRNG
 * @include: libgretl.h
 *
 * Libgretl uses xoshiro256+ as its underlying engine for
 * uniform random values, but offers added value in the
 * form of generators for several distributions commonly
 * used in econometrics.
 *
 * Note that before using the libgretl PRNG you must call
 * either libgretl_init() or the specific initialization
 * function shown below, gretl_rand_init(). And once you're
 * finished with it, you may call gretl_rand_free() or the
 * global libgretl function libgretl_cleanup().
 */

static uint64_t xor_seed;

/* alternate RNG */
static uint64_t alt_xor_seed;

static inline double double_from_uint64 (uint64_t u);

static inline double xor_01 (void)
{
    return double_from_uint64(xor_i64());
}

static inline double double_from_uint64 (uint64_t u)
{
    /* Set the exponent to 0x3FF (for 1.0) and the sign bit to 0;
       the remaining 52 bits are from the random uint64_t.
    */
    const uint64_t MASK = 0x3FF0000000000000ULL;
    uint64_t bits = MASK | (u >> 12);
    double ret;

    /* create a double in the range [1.0, 2.0) */
    memcpy(&ret, &bits, sizeof(double));

    return ret - 1.0;
}

static uint32_t xor_i32 (void)
{
    return (uint32_t) (xor_i64() >> 32);
}

/**
 * gretl_rand_init:
 *
 * Initialize gretl's PRNG, using the system time as seed.
 */

void gretl_rand_init (void)
{
    char *fseed = getenv("GRETL_FORCE_SEED");

    if (fseed != NULL) {
        xor_seed = (uint64_t) atoi(fseed);
    } else {
        xor_seed = (uint64_t) time(NULL);
    }

    set_xor_state(xor_seed);
}

void gretl_rand_free (void)
{
    /* FIXME */
    return;
}

#ifdef HAVE_MPI

/**
 * gretl_multi_rng_init:
 *
 * Initialize RNG per MPI process, if needed.
 */

void gretl_multi_rng_init (int n, int self, guint64 seed)
{
    int i, err;

    if (self == 0 && seed > 0) {
        set_xor_state(seed);
    } else if (self == 0) {
        /* automatic seeding */
#ifdef G_OS_WIN32
        /* consider using BCryptGenRandom() ? */
        uint64_t u = time(NULL);
#else
        int fd = open("/dev/urandom", O_RDONLY);
        uint64_t u;

        if (fd == -1) {
            u = time(NULL);
        } else {
            ssize_t sz = read(fd, &u, sizeof(uint64_t));
            if (sz == -1) {
                u = time(NULL);
            }
        }
#endif
        set_xor_state(u);
    }

    err = gretl_mpi_bcast_rng(xor_state, 0);
    if (!err) {
        for (i=0; i<self; i++) {
            xor_jump();
        }
    }
#if 0
    printf("rank %d: jumped or_state[1] = %" G_GUINT64_FORMAT "\n",
           self, xor_state[1]);
#endif
}

#endif /* HAVE_MPI */

/**
 * gretl_rand_get_seed:
 *
 * Returns: the value of the seed for gretl's PRNG.
 */

guint64 gretl_rand_get_seed (void)
{
    return xor_seed;
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

void gretl_rand_set_seed (guint64 seed)
{
    seed = (seed == 0)? time(NULL) : seed;
    set_xor_state((uint64_t) seed);
}

void gretl_alt_rand_set_seed (guint64 seed)
{
    seed = (seed == 0)? time(NULL) : seed;
    alt_xor_seed = seed;
    // FIXME
}

/**
 * gretl_rand_01:
 *
 * Returns: the next random double, equally distributed over
 * the range [0, 1).
 */

double gretl_rand_01 (void)
{
    return double_from_uint64(xor_i64());
}

/* Select which 32 bit generator to use for Ziggurat */

static inline uint32_t randi32 (void)
{
    return xor_i32();
}

#if !(HAVE_X86_32)

/* 53 bits for mantissa + 1 bit sign */

static uint64_t randi54 (void)
{
    const uint64_t u = xor_i64();
    const uint64_t mask = (1ULL << 54) - 1;

    return u & mask;
}

#endif

/* generates a uniform random double on (0,1) with 53-bit resolution */

static double randu53 (void)
{
    const uint64_t u = xor_i64() >> 11;

    return (double) u * 0x1.0p-53;
}

/* Ziggurat normal generator: this Ziggurat code here is shamelessly
   filched from GNU Octave (with minor modifications), since it turned
   out to be faster than what we had (until 2016-10-08), with no loss
   in quality of the random normal variates (in fact, perhaps a gain).

   See randmtzig.c in the Octave code-base. It appears that David
   Bateman was the main author; anyway, thanks to whoever wrote the
   implementation!
*/

#define ZIGGURAT_TABLE_SIZE 256
#define ZIGGURAT_NOR_R 3.6541528853610088
#define ZIGGURAT_NOR_INV_R 0.27366123732975828
#define NOR_SECTION_AREA 0.00492867323399

#define ZIGINT uint64_t
#define NMANTISSA 9007199254740992.0  /* 53 bit mantissa */

static ZIGINT ki[ZIGGURAT_TABLE_SIZE];
static double wi[ZIGGURAT_TABLE_SIZE];
static double fi[ZIGGURAT_TABLE_SIZE];

static int initt = 1;

static void create_ziggurat_tables (void)
{
    int i;
    double x, x1;

    x1 = ZIGGURAT_NOR_R;
    wi[255] = x1 / NMANTISSA;
    fi[255] = exp(-0.5 * x1 * x1);

    /* Index zero is special for tail strip, where Marsaglia and Tsang
       defines this as:
       k_0 = 2^31 * r * f(r) / v, w_0 = 0.5^31 * v / f(r), f_0 = 1,
       where v is the area of each strip of the ziggurat.
    */
    ki[0] = (ZIGINT) (x1 * fi[255] / NOR_SECTION_AREA * NMANTISSA);
    wi[0] = NOR_SECTION_AREA / fi[255] / NMANTISSA;
    fi[0] = 1.;

    for (i = 254; i > 0; i--) {
	/* New x is given by x = f^{-1}(v/x_{i+1} + f(x_{i+1})), thus
	   need inverse operator of y = exp(-0.5*x*x) -> x = sqrt(-2*ln(y))
	*/
	x = sqrt (-2. * log(NOR_SECTION_AREA / x1 + fi[i+1]));
	ki[i+1] = (ZIGINT) (x / x1 * NMANTISSA);
	wi[i] = x / NMANTISSA;
	fi[i] = exp(-0.5 * x * x);
	x1 = x;
    }

    ki[1] = 0;

    initt = 0;
}

/**
 * gretl_one_snormal:
 *
 * Returns: a single drawing from the standard normal distribution.
 */

double gretl_one_snormal (void)
{
    if (initt) {
	create_ziggurat_tables();
    }

    while (1) {
#if HAVE_X86_32
	/* Specialized for x86 32-bit architecture: 53-bit mantissa,
	   1-bit sign
	*/
	double x;
	int si, idx;
	uint32_t lo, hi;
	int64_t rabs;
	uint32_t *p = (uint32_t *) &rabs;

	lo = randi32();
	idx = lo & 0xFF;
	hi = randi32();
	si = hi & UMASK;
	p[0] = lo;
	p[1] = hi & 0x1FFFFF;
	x = (si ? -rabs : rabs) * wi[idx];
#else
	const uint64_t r = randi54();
	const int64_t rabs = r >> 1;
	const int idx = (int) (rabs & 0xFF);
	const double x = ((r & 1) ? -rabs : rabs) * wi[idx];
#endif

	if (rabs < (int64_t) (ki[idx])) {
	    return x; /* 99.3% of the time we return here 1st try */
	} else if (idx == 0) {
	    /* As stated in Marsaglia and Tsang:
	       For the normal tail, the method of Marsaglia provides:
	       generate x = -ln(U_1)/r, y = -ln(U_2), until y+y > x*x,
	       then return r+x. Except that r+x is always in the positive
	       tail. Anything random might be used to determine the
	       sign, but as we already have r we might as well use it.
	    */
	    double xx, yy;

	    do {
		xx = - ZIGGURAT_NOR_INV_R * log(randu53());
		yy = - log(randu53());
            } while (yy+yy <= xx*xx);
	    return (rabs & 0x100) ? -ZIGGURAT_NOR_R-xx : ZIGGURAT_NOR_R+xx;
        } else if ((fi[idx-1] - fi[idx]) * randu53() + fi[idx] < exp(-0.5*x*x)) {
	    return x;
	}
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
 * for uniform input and the Ziggurat method for converting to the
 * normal distribution.
 */

void gretl_rand_normal (double *a, int t1, int t2)
{
    int t;

    for (t=t1; t<=t2; t++) {
	a[t] = gretl_one_snormal();
    }
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
 * deviation, using the Mersenne Twister for uniform input and
 * the Ziggurat method for converting to the normal distribution.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int gretl_rand_normal_full (double *a, int t1, int t2,
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

static guint32 rand_int_range (guint32 begin,
                               guint32 end,
                               int alt)
{
    guint32 dist = end - begin;
    guint32 rval = 0;

    if (dist > 0) {
	/* maxval is set to the predecessor of the greatest
	   multiple of dist less than or equal to 2^32
	*/
	guint32 maxval;

	if (dist <= 0x80000000u) { /* 2^31 */
	    /* maxval = 2^32 - 1 - (2^32 % dist) */
	    guint32 rem = (0x80000000u % dist) * 2;

	    if (rem >= dist) rem -= dist;
	    maxval = 0xffffffffu - rem;
	} else {
	    maxval = dist - 1;
	}

	if (alt) {
	    do {
                /* FIXME */
                //rval = alt_xor_i32();
		rval = xor_i32();
	    } while (rval > maxval);
	} else {
	    do {
		rval = xor_i32();
	    } while (rval > maxval);
	}

	rval %= dist;
    }

    return begin + rval;
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
        a[t] = double_from_uint64(xor_i64()) * (max - min) + min;
    }

    return 0;
}

static int real_gretl_rand_int_minmax (int *a, int n,
				       int min, int max,
				       int alt)
{
    int i, err = 0;

    if (max < min) {
	err = E_INVARG;
    } else if (min == max) {
	for (i=0; i<n; i++) {
	    a[i] = min;
	}
    } else {
	int offset = 0;

	if (min < 0) {
	    offset = -min;
	    max += offset;
	    min += offset;
	}

	for (i=0; i<n; i++) {
	    a[i] = rand_int_range(min, max + 1, alt) - offset;
	}
    }

    return err;
}

/**
 * gretl_rand_int_minmax:
 * @a: target array.
 * @n: length of array.
 * @min: lower closed bound of range.
 * @max: upper closed bound of range.
 *
 * Fill array @a of length @n with pseudo-random drawings
 * from the uniform distribution on [@min, @max], using the
 * Mersenne Twister.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int gretl_rand_int_minmax (int *a, int n, int min, int max)
{
    return real_gretl_rand_int_minmax(a, n, min, max, 0);
}

/**
 * gretl_alt_rand_int_minmax:
 * @a: target array.
 * @n: length of array.
 * @min: lower closed bound of range.
 * @max: upper closed bound of range.
 *
 * Fill array @a of length @n with pseudo-random drawings
 * from the uniform distribution on [@min, @max], using a
 * Mersenne Twister which is independent of the main one
 * employed by libgretl.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int gretl_alt_rand_int_minmax (int *a, int n, int min, int max)
{
    return real_gretl_rand_int_minmax(a, n, min, max, 1);
}

static int already_selected (double *a, int n, double val,
			     int offset)
{
    int i;

    for (i=0; i<n; i++) {
	if (offset > 0) {
	    if (a[i] == val - offset) {
		return 1;
	    }
	} else if (a[i] == val) {
	    return 1;
	}
    }

    return 0;
}

/**
 * gretl_rand_uniform_int_minmax:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @min: lower closed bound of range.
 * @max: upper closed bound of range.
 * @opt: include OPT_O for sampling without replacement.
 *
 * Fill the selected subset of array @a with pseudo-random drawings
 * from the uniform distribution on [@min, @max], using the
 * Mersenne Twister.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int gretl_rand_uniform_int_minmax (double *a, int t1, int t2,
				   int min, int max,
				   gretlopt opt)
{
    int t, err = 0;

    if (max < min) {
	err = E_INVARG;
    } else if (max == min) {
	for (t=t1; t<=t2; t++) {
	    a[t] = min;
	}
    } else {
	int i = 0, offset = 0;
	double x;

	if (min < 0) {
	    offset = -min;
	    max += offset;
	    min += offset;
	}

	for (t=t1; t<=t2; t++) {
	    x = rand_int_range(min, max + 1, 0);
	    if (opt & OPT_O) {
		while (already_selected(a, i, x, offset)) {
		    x = rand_int_range(min, max + 1, 0);
		}
	    }
	    a[t] = x - offset;
	    i++;
	}
    }

    return err;
}

/**
 * gretl_rand_uniform:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the uniform distribution on [0-1), using the Mersenne
 * Twister.
 */

void gretl_rand_uniform (double *a, int t1, int t2)
{
    int t;

    for (t=t1; t<=t2; t++) {
        a[t] = double_from_uint64(xor_i64());
    }
}

static double gretl_rand_uniform_one (void)
{
    return double_from_uint64(xor_i64());
}

double gretl_rand_gamma_one (double shape, double scale)
{
    double k = shape;
    double d, c, x, v, u, dv;

    if (shape <= 0 || scale <= 0) {
	return NADBL;
    }

    if (shape < 1) {
	k = shape + 1.0;
    }

    d = k - 1.0/3;
    c = 1.0 / sqrt(9*d);

    while (1) {
	x = gretl_one_snormal();
	v = pow(1 + c*x, 3);
	if (v > 0.0) {
	    dv = d * v;
	    u = gretl_rand_01();
	    /* apply squeeze */
	    if (u < 1 - 0.0331 * pow(x, 4) ||
		log(u) < 0.5*x*x + d*(1-v+log(v))) {
		break;
	    }
	}
    }
    if (shape < 1) {
	u = gretl_rand_01();
	dv *= pow(u, 1/shape);
    }

    return dv * scale;
}

/* Marsaglia-Tsang, "A Simple Method for Generating Gamma Variables",
   ACM Transactions on Mathematical Software, Vol. 26, No. 3,
   September 2000, Pages 363­372.

   (1) Setup: d=a-1/3, c=1/sqrt(9*d).
   (2) Generate v=(1+c*x)^3 with x normal; repeat if v <= 0.
   (3) Generate uniform U.
   (4) If U < 1-0.0331*x^4 return d*v.
   (5) If log(U) < 0.5*x^2+d*(1-v+log(v)) return d*v.
   (6) Go to step 2.
*/

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
    double k = shape;
    double d, c, x, v, u, dv;
    int t;

    if (shape <= 0 || scale <= 0) {
	return E_DATA;
    }

    if (shape < 1) {
	k = shape + 1.0;
    }

    d = k - 1.0/3;
    c = 1.0 / sqrt(9*d);

    for (t=t1; t<=t2; t++) {
	while (1) {
	    x = gretl_one_snormal();
	    v = pow(1 + c*x, 3);
	    if (v > 0.0) {
		dv = d * v;
		u = gretl_rand_01();
		/* apply squeeze */
		if (u < 1 - 0.0331 * pow(x, 4) ||
		    log(u) < 0.5*x*x + d*(1-v+log(v))) {
		    break;
		}
	    }
	}
	if (shape < 1) {
	    u = gretl_rand_01();
	    dv *= pow(u, 1/shape);
	}
	a[t] = dv * scale;
    }

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
 * from the Chi-squared distribution with @v degrees of freedom,
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
 * using the Mersenne Twister for uniform input and the ziggurat
 * method for converting to the normal distribution.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_student (double *a, int t1, int t2, double v)
{
    double *X2 = NULL;
    int T = t2 - t1 + 1;
    int t;

    if (v <= 0) {
	return E_INVARG;
    }

    X2 = malloc(T * sizeof *X2);
    if (X2 == NULL) {
	return E_ALLOC;
    }

    gretl_rand_normal(a, t1, t2);
    gretl_rand_gamma(X2, 0, T-1, 0.5*v, 2);

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
 * using the Mersenne Twister for uniform input and the Ziggurat
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
    int t;

    if (n < 0 || p < 0 || p > 1) {
	return E_INVARG;
    }

    if (n == 0 || p == 0.0) {
	for (t=t1; t<=t2; t++) {
	    a[t] = 0.0;
	}
    } else if (p == 1.0) {
	for (t=t1; t<=t2; t++) {
	    a[t] = n;
	}
    } else {
	double *b = malloc(n * sizeof *b);
	int i;

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
    }

    return 0;
}

static double gretl_rand_binomial_one (int n, double p,
				       double *b)
{
    double ret;

    if (n < 0 || p < 0 || p > 1) {
	return NADBL;
    }

    if (n == 0 || p == 0.0) {
	ret = 0.0;
    } else if (p == 1.0) {
	ret = n;
    } else {
	int i;

	ret = 0.0;
	gretl_rand_uniform(b, 0, n - 1);
	for (i=0; i<n; i++) {
	    if (b[i] <= p) {
		ret += 1;
	    }
	}
    }

    return ret;
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

	x = exp(m) * gretl_rand_01();
	while (x > 1) {
	    y++;
	    x *= gretl_rand_01();
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
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_poisson (double *a, int t1, int t2, const double *m,
			int vec)
{
    double mt;
    int t;

    for (t=t1; t<=t2; t++) {
	mt = (vec)? m[t] : *m;
	a[t] = (mt <= 0)? NADBL : genpois(mt);
    }

    return 0;
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
	    u = gretl_rand_01();
	    while (u == 0.0) {
		u = gretl_rand_01();
	    }
	    a[t] = scale * pow(-log(u), kinv);
	}
    }

    return err;
}

/**
 * gretl_rand_exponential:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @mu: scale parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the exponential distribution with scale @mu.
 *
 * Returns: 0 on success, non-zero if @mu is out of
 * bounds.
 */

int gretl_rand_exponential (double *a, int t1, int t2, double mu)
{
    int err = 0;

    if (mu <= 0) {
	err = E_DATA;
    } else {
	double u;
	int t;

	for (t=t1; t<=t2; t++) {
	    u = gretl_rand_01();
	    while (u == 0.0) {
		u = gretl_rand_01();
	    }
	    a[t] = -mu * log(u);
	}
    }

    return err;
}

/**
 * gretl_rand_logistic:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @loc: location parameter.
 * @scale: scale parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the logistic distribution with location @loc and scale @scale.
 *
 * Returns: 0 on success, non-zero if @scale is out of
 * bounds.
 */

int gretl_rand_logistic (double *a, int t1, int t2,
			 double loc, double scale)
{
    int err = 0;

    if (scale <= 0) {
	err = E_DATA;
    } else {
	double u;
	int t;

	for (t=t1; t<=t2; t++) {
	    u = gretl_rand_01();
	    while (u == 0.0) {
		u = gretl_rand_01();
	    }
	    a[t] = loc + scale * log(u / (1 - u));
	}
    }

    return err;
}

/**
 * gretl_rand_GED:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @nu: shape parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the GED distribution with shape @nu. We exploit the fact that
 * if x ~ GED(n), then |x/k|^n is a Gamma rv.
 *
 * Returns: 0 on success, non-zero if @nu is out of bounds.
 */

int gretl_rand_GED (double *a, int t1, int t2, double nu)
{
    int err, t;
    double p, scale;

    if (nu < 0) {
	return E_INVARG;
    }

    p = 1.0/nu;
    scale = pow(0.5, p) * sqrt(gammafun(p) / gammafun(3.0*p));
    err = gretl_rand_gamma(a, t1, t2, p, 2);

    if (!err) {
	for (t=t1; t<=t2; t++) {
	    a[t] = scale * pow(a[t], p);
	    if (gretl_rand_01() < 0.5) {
		a[t] = -a[t];
	    }
	}
    }

    return err;
}

/**
 * gretl_rand_laplace:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @mu: mean.
 * @b: shape parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Laplace distribution with mean @mu and scale @b.
 *
 * Returns: 0 on success, non-zero if @b is out of bounds.
 */

int gretl_rand_laplace (double *a, int t1, int t2,
			double mu, double b)
{
    int t, sgn;
    double U;

    if (b < 0) {
	return E_INVARG;
    }

    /* uniform on [0,1) */
    gretl_rand_uniform(a, t1, t2);

    for (t=t1; t<=t2; t++) {
	/* convert to (-1/2,1/2] */
	U = 0.5 - a[t];
	sgn = U < 0 ? -1 : 1;
	a[t] = mu - b*sgn * log(1 - 2*fabs(U));
    }

    return 0;
}

/**
 * gretl_rand_beta:
 * @x: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @s1: shape parameter > 0.
 * @s2: shape parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the beta distribution with shape parameters @s1 and @s2.
 * The code here is adapted from http://www.netlib.org/random/random.f90
 * which implements the method of R.C.H. Cheng, "Generating beta
 * variates with nonintegral shape parameters", Communications of the
 * ACM, 21(4), April 1978.
 *
 * Returns: 0 on success, non-zero if @s1 or @s2 are out of bounds.
 */

int gretl_rand_beta (double *x, int t1, int t2,
		     double s1, double s2)
{
    double aln4 = 1.3862944;
    double a, b, s, u, v, y, z;
    double d, f, h, u0, c;
    int t, swap;
    double val;

    if (s1 <= 0 || s2 <= 0) {
	return E_DATA;
    }

    /* initialization */
    a = s1;
    b = s2;
    swap = b > a;
    if (swap) {
	f = b;
	b = a;
	a = f;
    }
    d = a/b;
    f = a+b;
    if (b > 1) {
	h = sqrt((2*a*b - f)/(f - 2.0));
	u0 = 1.0;
    } else {
	h = b;
	u0 = 1.0/(1.0 + pow(a/(DBL_MAX*b), b));
    }
    c = a+h;

    /* generation */
    for (t=t1; t<=t2; t++) {
	while (1) {
	    u = gretl_rand_uniform_one();
	    v = gretl_rand_uniform_one();
	    s = u * u * v;
	    if (u < DBL_MIN || s <= 0) continue;
	    if (u < u0) {
		v = log(u/(1.0 - u))/h;
		y = d*exp(v);
		z = c*v + f*log((1.0 + d)/(1.0 + y)) - aln4;
		if (s - 1.0 > z) {
		    if (s - s*z > 1.0) continue;
		    if (log(s) > z) continue;
		}
		val = y/(1.0 + y);
	    } else {
		if (4.0*s > pow(1.0 + 1.0/d, f)) continue;
		val = 1.0;
	    }
	    break;
	}
	x[t] = swap ? (1.0 - val) : val;
    }

    return 0;
}

/**
 * gretl_rand_beta_binomial:
 * @x: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @n: number of trials.
 * @s1: beta shape parameter > 0.
 * @s2: beta shape parameter > 0.
 *
 * Fill the selected range of array @x with pseudo-random drawings
 * from the binomial distribution with @n trials and success
 * probability distributed according to the beta distribution with
 * shape parameters @s1 and @s2.
 *
 * Returns: 0 on success, non-zero if @n, @s1 or @s2 are out of bounds.
 */

int gretl_rand_beta_binomial (double *x, int t1, int t2,
			      int n, double s1, double s2)
{
    int t, err;

    err = gretl_rand_beta(x, t1, t2, s1, s2);

    if (!err) {
	double *b = malloc(n * sizeof *b);

	if (b == NULL) {
	    return E_ALLOC;
	} else {
	    for (t=t1; t<=t2; t++) {
		x[t] = gretl_rand_binomial_one(n, x[t], b);
	    }
	    free(b);
	}
    }

    return err;
}

#define DEBUG 0

/**
 * gretl_rand_discrete:
 * @x: target vector.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @p: probabilities vector.
 *
 * Fill the selected range of array @x with pseudo-random drawings
 * from the discrete distribution defined by @p, with support from 1
 * to @n.
 *
 * Returns: 0 on success, non-zero if @p is not a proper probability vector.
 */

int gretl_rand_discrete (double *x, int t1, int t2,
                         const gretl_vector *p)
{
    gretl_matrix *S = NULL;
    gretl_matrix *sS = NULL;
    double *cp;
    double u;
    int n, nr = t2 - t1 + 1;
    int i, j, t, pos;
    int err = 0;

    n = gretl_vector_get_length(p);

    /* sanity check on @p */
    for (i=0; i<n; i++) {
        if (p->val[i] < 0.0 || p->val[i] > 1.0) {
            fprintf(stderr, "gretl_rand_discrete: p is outside[0,1]\n");
            return E_INVARG;
        }
    }

    cp = malloc(n * sizeof *cp);

    if (cp == NULL) {
        err = E_ALLOC;
    } else {
        cp[0] = p->val[0];
        for (i=1; i<n; i++) {
            cp[i] = cp[i-1] + p->val[i];
        }
        if (fabs(1.0 - cp[n-1]) > 1.0e-7) {
            free(cp);
            err = E_INVARG;
        }
    }

    if (err) {
        goto bailout;
    }

    /* force cp[n-1] to 1 to avoid numerical issues */
    cp[n-1] = 1.0;
#if DEBUG
    for (i=0; i<n; i++) {
        printf("cp[%d] = %11.7f\n", i, cp[i]);
    }
#endif

    /* build a matrix with @nr uniform rvs sorted ascendingly
       in column 2 and the corresponding row index in column 1
    */
    S = gretl_random_matrix_new(nr, 2, D_UNIFORM);
    if (S == NULL) {
        err = E_ALLOC;
        goto bailout;
    } else {
        for (i=0; i<nr; i++) {
            S->val[i] = i;
        }
#if DEBUG
        gretl_matrix_print(S, "S");
#endif
    }

    sS = gretl_matrix_sort_by_column(S, 1, &err);
    gretl_matrix_free(S);
    if (err) {
        goto bailout;
    }

#if DEBUG
    gretl_matrix_print(sS, "sS");
#endif

    /* go through the rows of @sS and assign the proper value
       to the elements of @x
    */
    i = pos = 0;
    for (t=t1; t<=t2; t++) {
        j = (int)(sS->val[i]);
        u = gretl_matrix_get(sS, i++, 1);
        while (u > cp[pos]) {
            pos++;
        }
#if DEBUG
        fprintf(stderr, "u = %8.4f, j = %d, pos = %d\n", u, j, pos);
#endif
        x[t1+j] = pos+1;
    }

 bailout:

    free(cp);
    gretl_matrix_free(sS);

    return err;
}

/**
 * gretl_rand_dirichlet:
 * @a: parameter vector, length k.
 * @n: number of rows (replications) in return matrix.
 * @err: location to receive error code.
 *
 * Returns: on success, an @n x @k matrix matrix containing
 * @n drawings from the Dirichlet distribution of order k.
 */

gretl_matrix *gretl_rand_dirichlet (const gretl_vector *a,
				    int n, int *err)
{
    int k = gretl_vector_get_length(a);
    gretl_matrix *D = NULL;

    if (k < 2) {
	*err = E_NONCONF;
	return NULL;
    }

    D = gretl_matrix_alloc(n, k);

    if (D == NULL) {
	*err = E_ALLOC;
    } else {
	double dij, rsum;
	int t1 = 0, t2 = n-1;
	int i, j;

	for (j=0; j<k && !*err; j++) {
	    *err = gretl_rand_gamma(D->val, t1, t2, a->val[j], 1.0);
	    t1 += n;
	    t2 += n;
	}
	if (*err) {
	    gretl_matrix_free(D);
	    D = NULL;
	} else {
	    for (i=0; i<n; i++) {
		rsum = 0.0;
		for (j=0; j<k; j++) {
		    rsum += gretl_matrix_get(D, i, j);
		}
		for (j=0; j<k; j++) {
		    dij = gretl_matrix_get(D, i, j);
		    gretl_matrix_set(D, i, j, dij / rsum);
		}
	    }
	}
    }

    return D;
}

/**
 * gretl_rand_int_max:
 * @max: the maximum value (open)
 *
 * Returns: a pseudo-random unsigned int in the interval
 * [0, max-1].
 */

guint32 gretl_rand_int_max (unsigned int max)
{
    return rand_int_range(0, max, 0);
}

/**
 * gretl_rand_int:
 *
 * Returns: a pseudo-random unsigned int on the interval
 * [0, 2^32-1].
 */

guint32 gretl_rand_int (void)
{
    return xor_i32();
}

guint32 gretl_alt_rand_int (void)
{
    // FIXME
    return 0;
    // return sfmt_alt_rand32();
}

static double halton (int i, int base)
{
    double f = 1.0 / base;
    double h = 0.0;

    while (i > 0) {
	h += f * (i % base);
	i = floor(i / base);
	f /= base;
    }

    return h;
}

/* This is simple, and no doubt could be more efficient, but it's fast
   up to the number of primes that might reasonably be expected in the
   context of Halton sequences. The incoming @p is prime, and
   therefore odd. Starting at the next odd number q and skipping
   evens, we check for odd factors (since odd numbers cannot have even
   factors) from 3 to a maximum of q/3. If none are found then q is
   prime and we return it. There's a guard against running out of
   positive 32-bit integers; in that case we return 0.
*/

static int next_prime (int p)
{
    int f, fmax, pr, q;

    for (q=p+2; q<INT_MAX-1; q+=2) {
	fmax = q/3;
	pr = 1;
	for (f=3; f<=fmax; f+=2) {
	    if (q % f == 0) {
		pr = 0;
		break;
	    }
	}
	if (pr) {
	    return q;
	}
    }

    return 0;
}

/**
 * halton_matrix:
 * @m: number of rows (sequences).
 * @r: number of columns (elements in each sequence).
 * @offset: the number of initial values to discard.
 * @err: location to receive error code.
 *
 * Returns: an @m x @r matrix containing @m Halton
 * sequences. The sequences are contructed using the first
 * @m primes, and the first @offset elements of each sequence
 * are discarded.
 */

gretl_matrix *halton_matrix (int m, int r, int offset, int *err)
{
    const int bases[50] = {
	/* the first 50 primes */
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
	53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
	109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
	173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229
    };
    gretl_matrix *H;
    double hij;
    int i, j, k, n, p = 0;

    if (offset < 0 || m <= 0 || r <= 0) {
	*err = E_INVARG;
	return NULL;
    }

    H = gretl_matrix_alloc(m, r);
    if (H == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* we'll discard the first @offset elements */
    n = r + offset;

    for (i=0; i<m; i++) {
	if (i < 50) {
	    p = bases[i];
	} else {
	    p = next_prime(p);
	    if (p == 0) {
		*err = E_INVARG;
		break;
	    }
	}
	j = 0;
	for (k=1; k<n; k++) {
	    hij = halton(k, p);
	    if (k >= offset) {
		gretl_matrix_set(H, i, j++, hij);
	    }
	}
    }

    if (*err) {
	gretl_matrix_free(H);
	H = NULL;
    }

    return H;
}

static int wishart_workspace (gretl_matrix **pW,
			      gretl_matrix **pB,
			      double **pZ,
			      int p)
{
    int err = 0;

    *pW = gretl_matrix_alloc(p, p);
    *pB = gretl_matrix_alloc(p, p);

    if (*pW == NULL || *pB == NULL) {
        err = E_ALLOC;
    } else {
	int n = p * (p + 1) / 2;

	*pZ = malloc(n * sizeof **pZ);
	if (*pZ == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_rand_normal(*pZ, 0, n - 1);
	}
    }

    if (err) {
	gretl_matrix_free(*pW);
	gretl_matrix_free(*pB);
    }

    return err;
}

static int wishart_inverse_finalize (const gretl_matrix *C,
				     gretl_matrix *B,
				     gretl_matrix *W)
{
    gretl_matrix_qform(C, GRETL_MOD_NONE,
		       B, W, GRETL_MOD_NONE);
    return gretl_invpd(W);
}

static void odell_feiveson_compute (gretl_matrix *B,
				    double *Z,
				    int v)
{
    double Xi, Zri, Zrj;
    double wii, wij;
    int p = B->rows;
    int i, j, r;

    for (i=0; i<p; i++) {
	gretl_rand_chisq(&Xi, 0, 0, v - i);
	wii = Xi;
	for (r=0; r<i; r++) {
	    Zri = Z[ijton(r, i, p)];
	    wii += Zri * Zri;
	}
	gretl_matrix_set(B, i, i, wii);
	for (j=i+1; j<p; j++) {
	    wij = Z[ijton(i, j, p)] * sqrt(Xi);
	    for (r=0; r<i; r++) {
		Zri = Z[ijton(r, i, p)];
		Zrj = Z[ijton(r, j, p)];
		wij += Zri * Zrj;
	    }
	    gretl_matrix_set(B, i, j, wij);
	    gretl_matrix_set(B, j, i, wij);
	}
    }
}

/**
 * inverse_wishart_matrix:
 * @S: p x p positive definite scale matrix.
 * @v: degrees of freedom: v >= p.
 * @err: location to receive error code.
 *
 * Computes a draw from the Inverse Wishart distribution
 * with @v degrees of freedom, using the method of Odell and
 * Feiveson, "A numerical procedure to generate a sample
 * covariance matrix", JASA 61, pp. 199­203, 1966.
 *
 * Returns: a p x p Inverse Wishart matrix, or NULL on error.
 */

gretl_matrix *inverse_wishart_matrix (const gretl_matrix *S,
				      int v, int *err)
{
    gretl_matrix *C, *B = NULL, *W = NULL;
    double *Z = NULL;

    if (S == NULL || S->cols != S->rows || v < S->rows) {
	*err = E_INVARG;
	return NULL;
    }

    C = gretl_matrix_copy(S);
    if (C == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    /* invert and decompose C */
    *err = cholesky_factor_of_inverse(C);

    if (!*err) {
	*err = wishart_workspace(&W, &B, &Z, S->rows);
    }

    if (*err) {
	gretl_matrix_free(C);
	return NULL;
    }

    odell_feiveson_compute(B, Z, v);
    *err = wishart_inverse_finalize(C, B, W);

    if (*err) {
	gretl_matrix_free(W);
	W = NULL;
    }

    gretl_matrix_free(B);
    gretl_matrix_free(C);
    free(Z);

    return W;
}

static void vech_into_row (gretl_matrix *targ, int row,
			   const gretl_matrix *m)
{
    double mij;
    int n = m->rows;
    int i, j, jj = 0;

    for (i=0; i<n; i++) {
	for (j=i; j<n; j++) {
	    mij = gretl_matrix_get(m, i, j);
	    gretl_matrix_set(targ, row, jj++, mij);
	}
    }
}

/**
 * inverse_wishart_sequence:
 * @S: p x p positive definite scale matrix.
 * @v: degrees of freedom: v >= p.
 * @n: number of replications.
 * @err: location to receive error code.
 *
 * Computes @n draws from the Inverse Wishart distribution
 * with @v degrees of freedom, using the method of Odell and
 * Feiveson, "A numerical procedure to generate a sample
 * covariance matrix", JASA 61, pp. 199­203, 1966.
 *
 * Returns: an n x k matrix, where k = p*(p+1)/2, each row
 * being the vech of an Inverse Wishart matrix, or NULL on error.
 */

gretl_matrix *inverse_wishart_sequence (const gretl_matrix *S,
					int v, int n, int *err)
{
    gretl_matrix *C, *B = NULL, *W = NULL;
    gretl_matrix *Seq = NULL;
    double *Z = NULL;
    int i, k = 0;

    if (S == NULL || S->cols != S->rows || v < S->rows) {
	*err = E_INVARG;
	return NULL;
    }

    if (n < 1) {
	*err = E_INVARG;
	return NULL;
    }

    C = gretl_matrix_copy(S);
    if (C == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    /* invert and decompose C */
    *err = cholesky_factor_of_inverse(C);

    if (!*err) {
	*err = wishart_workspace(&W, &B, &Z, S->rows);
    }

    if (!*err) {
	k = C->rows * (C->rows + 1) / 2;
	Seq = gretl_matrix_alloc(n, k);
	if (Seq == NULL) {
	    *err = E_ALLOC;
	}
    }

    for (i=0; i<n && !*err; i++) {
	odell_feiveson_compute(B, Z, v);
	*err = wishart_inverse_finalize(C, B, W);
        if (!*err) {
	    /* write vech of W into row i of Seq */
	    vech_into_row(Seq, i, W);
	    if (i < n - 1) {
		/* refresh Normal draws */
		gretl_rand_normal(Z, 0, k - 1);
	    }
	}
    }

    gretl_matrix_free(C);
    gretl_matrix_free(W);
    gretl_matrix_free(B);
    free(Z);

    if (*err && Seq != NULL) {
	gretl_matrix_free(Seq);
	Seq = NULL;
    }

    return Seq;
}

char *gretl_rand_hex_string (int len, int *err)
{
    const char *src = "0123456789abcdef";
    int *ivals = NULL;
    char *ret = NULL;
    int i;

    if (len < 0) {
	*err = E_INVARG;
	return NULL;
    } else if (len == 0) {
	return gretl_strdup("");
    }

    ivals = malloc(len * sizeof *ivals);
    ret = malloc(len + 1);

    if (ivals == NULL || ret == NULL) {
	*err = E_ALLOC;
	return ret;
    }

    real_gretl_rand_int_minmax(ivals, len, 0, 15, 0);

    for (i=0; i<len; i++) {
	ret[i] = src[ivals[i]];
    }

    ret[i] = '\0';
    free(ivals);

    return ret;
}
