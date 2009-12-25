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

/**
 * gretl_rand_get_seed:
 *
 * Returns: the value of the seed for gretl's PRNG.
 */

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

static double gretl_rand_01 (void)
{
    return g_rand_double(gretl_rand);
}

/**
 * gretl_one_snormal:
 *
 * Returns: a single drawing from the standard normal distribution.
 */

double gretl_one_snormal (void) 
{
    double x, y, z;

 tryagain:
    x = gretl_rand_01();
    y = gretl_rand_01();
    z = sqrt(-2. * log(x));
    if (isnan(z) || isinf(z)) {
	goto tryagain;
    }

    return (z * cos(M_2PI * y));
}

#define USE_ZIGGURAT 1

#if USE_ZIGGURAT

/* the following based on gauss.c - gaussian random numbers, using the Ziggurat method
 * Copyright (C) 2005 Jochen Voss (GPL'd).
 * 
 * See http://seehuhn.de/pages/ziggurat, and also
 *
 *     George Marsaglia, Wai Wan Tsang
 *     The Ziggurat Method for Generating Random Variables
 *     Journal of Statistical Software, vol. 5 (2000), no. 8
 *     http://www.jstatsoft.org/v05/i08/
 */

/* position of right-most step */
#define PARAM_R 3.44428647676

/* tabulated values for the height of the Ziggurat levels */
static const double ytab[128] = {
    1, 0.963598623011, 0.936280813353, 0.913041104253,
    0.892278506696, 0.873239356919, 0.855496407634, 0.838778928349,
    0.822902083699, 0.807732738234, 0.793171045519, 0.779139726505,
    0.765577436082, 0.752434456248, 0.739669787677, 0.727249120285,
    0.715143377413, 0.703327646455, 0.691780377035, 0.68048276891,
    0.669418297233, 0.65857233912, 0.647931876189, 0.637485254896,
    0.62722199145, 0.617132611532, 0.607208517467, 0.597441877296,
    0.587825531465, 0.578352913803, 0.569017984198, 0.559815170911,
    0.550739320877, 0.541785656682, 0.532949739145, 0.524227434628,
    0.515614886373, 0.507108489253, 0.498704867478, 0.490400854812,
    0.482193476986, 0.47407993601, 0.466057596125, 0.458123971214,
    0.450276713467, 0.442513603171, 0.434832539473, 0.427231532022,
    0.419708693379, 0.41226223212, 0.404890446548, 0.397591718955,
    0.390364510382, 0.383207355816, 0.376118859788, 0.369097692334,
    0.362142585282, 0.355252328834, 0.348425768415, 0.341661801776,
    0.334959376311, 0.328317486588, 0.321735172063, 0.31521151497,
    0.308745638367, 0.302336704338, 0.29598391232, 0.289686497571,
    0.283443729739, 0.27725491156, 0.271119377649, 0.265036493387,
    0.259005653912, 0.253026283183, 0.247097833139, 0.241219782932,
    0.235391638239, 0.229612930649, 0.223883217122, 0.218202079518,
    0.212569124201, 0.206983981709, 0.201446306496, 0.195955776745,
    0.190512094256, 0.185114984406, 0.179764196185, 0.174459502324,
    0.169200699492, 0.1639876086, 0.158820075195, 0.153697969964,
    0.148621189348, 0.143589656295, 0.138603321143, 0.133662162669,
    0.128766189309, 0.123915440582, 0.119109988745, 0.114349940703,
    0.10963544023, 0.104966670533, 0.100343857232, 0.0957672718266,
    0.0912372357329, 0.0867541250127, 0.082318375932, 0.0779304915295,
    0.0735910494266, 0.0693007111742, 0.065060233529, 0.0608704821745,
    0.056732448584, 0.05264727098, 0.0486162607163, 0.0446409359769,
    0.0407230655415, 0.0368647267386, 0.0330683839378, 0.0293369977411,
    0.0256741818288, 0.0220844372634, 0.0185735200577, 0.0151490552854,
    0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};

/* tabulated values for 2^24 times x[i]/x[i+1],
 * used to accept for U*x[i+1]<=x[i] without any floating point operations */
static const guint32 ktab[128] = {
    0, 12590644, 14272653, 14988939,
    15384584, 15635009, 15807561, 15933577,
    16029594, 16105155, 16166147, 16216399,
    16258508, 16294295, 16325078, 16351831,
    16375291, 16396026, 16414479, 16431002,
    16445880, 16459343, 16471578, 16482744,
    16492970, 16502368, 16511031, 16519039,
    16526459, 16533352, 16539769, 16545755,
    16551348, 16556584, 16561493, 16566101,
    16570433, 16574511, 16578353, 16581977,
    16585398, 16588629, 16591685, 16594575,
    16597311, 16599901, 16602354, 16604679,
    16606881, 16608968, 16610945, 16612818,
    16614592, 16616272, 16617861, 16619363,
    16620782, 16622121, 16623383, 16624570,
    16625685, 16626730, 16627708, 16628619,
    16629465, 16630248, 16630969, 16631628,
    16632228, 16632768, 16633248, 16633671,
    16634034, 16634340, 16634586, 16634774,
    16634903, 16634972, 16634980, 16634926,
    16634810, 16634628, 16634381, 16634066,
    16633680, 16633222, 16632688, 16632075,
    16631380, 16630598, 16629726, 16628757,
    16627686, 16626507, 16625212, 16623794,
    16622243, 16620548, 16618698, 16616679,
    16614476, 16612071, 16609444, 16606571,
    16603425, 16599973, 16596178, 16591995,
    16587369, 16582237, 16576520, 16570120,
    16562917, 16554758, 16545450, 16534739,
    16522287, 16507638, 16490152, 16468907,
    16442518, 16408804, 16364095, 16301683,
    16207738, 16047994, 15704248, 15472926
};

/* tabulated values of 2^{-24}*x[i] */
static const double wtab[128] = {
    1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
    3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
    3.8950989572e-08, 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
    4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
    5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
    5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
    5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
    6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
    6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08,
    6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
    7.2824062723e-08, 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
    7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
    7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
    8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
    8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08,
    8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
    9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
    9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
    9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07,
    1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
    1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
    1.09477300508e-07, 1.1042504257e-07, 1.11384564771e-07, 1.12356564007e-07,
    1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
    1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
    1.21808209468e-07, 1.2295639141e-07, 1.24129212952e-07, 1.25328445797e-07,
    1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
    1.31797105598e-07, 1.3320433736e-07, 1.34657379914e-07, 1.36160594606e-07,
    1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
    1.44635331499e-07, 1.4657889173e-07, 1.48632138436e-07, 1.50811780719e-07,
    1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
    1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
    1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};

static double ran_normal_ziggurat (void)
{
    guint32 U, sign, i, j;
    double x, y;

    while (1) {
	U = g_rand_int(gretl_rand);
	i = U & 0x0000007F;	/* 7 bits to choose the step */
	sign = U & 0x00000080;	/* 1 bit for the sign */
	j = U >> 8;		/* 24 bits for the x-value */

	x = j * wtab[i];

	if (j < ktab[i]) {
	    break;
	}

	if (i < 127) {
	    double y0 = ytab[i], y1 = ytab[i+1];

	    y = y1 + (y0 - y1) * gretl_rand_01();
	} else {
	    x = PARAM_R - log(1.0 - gretl_rand_01()) / PARAM_R;
	    y = exp(-PARAM_R * (x - 0.5 * PARAM_R)) * gretl_rand_01();
	}

	if (y < exp(-0.5 * x * x)) {
	    break;
	}
    }

    return sign ? x : -x;
}

/* end of Ziggurat code derived from Jochen Voss */

#else /* Box-Muller, as of old */

static void gretl_two_snormals (double *z1, double *z2) 
{
    double x, y, z;

 tryagain:
    x = 2 * gretl_rand_01() - 1;
    y = 2 * gretl_rand_01() - 1;
    z = x*x + y*y;
    if (z >= 1) {
	goto tryagain;
    }

    z = sqrt(-2 * log(z)/z);
    *z1 = z * x;
    *z2 = z * y;
}

#endif

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

#if USE_ZIGGURAT
    for (t=t1; t<=t2; t++) {
	a[t] = ran_normal_ziggurat();
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
 * gretl_rand_int_minmax:
 * @a: target array.
 * @n: length of array.
 * @min: lower closed bound of range.
 * @max: upper closed bound of range.
 *
 * Fill array @a of length @n with pseudo-random drawings
 * from the uniform distribution on @min to @max, using the
 * Mersenne Twister.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int gretl_rand_int_minmax (int *a, int n, int min, int max) 
{
    int i;

    if (max < min) {
	return E_INVARG;
    }

    /* note: in g_rand_int_range the upper bound is open */

    for (i=0; i<n; i++) {
	a[i] = g_rand_int_range(gretl_rand, min, max + 1);
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
	a[t] = gretl_rand_01();
    }
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
		u = gretl_rand_01();
		v = gretl_rand_01();
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
		U[i] = gretl_rand_01();
		while (U[i] == 0.0) {
		    U[i] = gretl_rand_01();
		}
	    }
	    if (delta > 0) {
		while (1) {
		    u = gretl_rand_01();
		    v = gretl_rand_01();
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


