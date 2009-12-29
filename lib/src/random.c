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
    useed = (seed == 0)? time(NULL) : seed;
    g_rand_set_seed(gretl_rand, useed);
}

static double gretl_rand_01 (void)
{
    return g_rand_double(gretl_rand);
}

/* Below: an implementation of the Marsaglia/Tsang Ziggurat method for
   generating random normal samples (Journal of Statistical Software,
   vol. 5, no. 8, 2000).

   The code is based on Jochen Voss's gauss.c, written for use with
   the Gnu Scientific Library -- see http://seehuhn.de/pages/ziggurat.
   It was modified for gretl in two main ways: (a) we use the Mersenne
   Twister uniform RNG from GLib, and (b) we use a 30-bit unsigned
   random integer for conversion via the Ziggurat where Voss uses
   a 24-bit value.  This sacrifices a little speed in exchange for
   better coverage of the real line: we get over 10^9 distinct 
   normal values as opposed to around 30 million values from the
   24-bit input.

   See also Jurgen Doornik's discussion of Ziggurat:
   http://www.doornik.com/research/ziggurat.pdf
*/

#define ZLEVELS 128

/* position of right-most step */
#define ZIG_R 3.44428647676

/* tabulated values for the height of the Ziggurat levels */

static const double z_ytab[ZLEVELS] = {
               1,  0.963598623011,  0.936280813353,  0.913041104253,
  0.892278506696,  0.873239356919,  0.855496407634,  0.838778928349,
  0.822902083699,  0.807732738234,  0.793171045519,  0.779139726505,
  0.765577436082,  0.752434456248,  0.739669787677,  0.727249120285,
  0.715143377413,  0.703327646455,  0.691780377035,  0.680482768910,
  0.669418297233,  0.658572339120,  0.647931876189,  0.637485254896,
  0.627221991450,  0.617132611532,  0.607208517467,  0.597441877296,
  0.587825531465,  0.578352913803,  0.569017984198,  0.559815170911,
  0.550739320877,  0.541785656682,  0.532949739145,  0.524227434628,
  0.515614886373,  0.507108489253,  0.498704867478,  0.490400854812,
  0.482193476986,  0.474079936010,  0.466057596125,  0.458123971214,
  0.450276713467,  0.442513603171,  0.434832539473,  0.427231532022,
  0.419708693379,  0.412262232120,  0.404890446548,  0.397591718955,
  0.390364510382,  0.383207355816,  0.376118859788,  0.369097692334,
  0.362142585282,  0.355252328834,  0.348425768415,  0.341661801776,
  0.334959376311,  0.328317486588,  0.321735172063,  0.315211514970,
  0.308745638367,  0.302336704338,  0.295983912320,  0.289686497571,
  0.283443729739,  0.277254911560,  0.271119377649,  0.265036493387,
  0.259005653912,  0.253026283183,  0.247097833139,  0.241219782932,
  0.235391638239,  0.229612930649,  0.223883217122,  0.218202079518,
  0.212569124201,  0.206983981709,  0.201446306496,  0.195955776745,
  0.190512094256,  0.185114984406,  0.179764196185,  0.174459502324,
  0.169200699492,  0.163987608600,  0.158820075195,  0.153697969964,
  0.148621189348,  0.143589656295,  0.138603321143,  0.133662162669,
  0.128766189309,  0.123915440582,  0.119109988745,  0.114349940703,
  0.109635440230,  0.104966670533,  0.100343857232,  0.0957672718266,
  0.0912372357329,  0.0867541250127,  0.082318375932,  0.0779304915295,
  0.0735910494266,  0.0693007111742,  0.065060233529,  0.0608704821745,
  0.0567324485840,  0.0526472709800,  0.0486162607163,  0.0446409359769,
  0.0407230655415,  0.0368647267386,  0.0330683839378,  0.0293369977411,
  0.0256741818288,  0.0220844372634,  0.0185735200577,  0.0151490552854,
  0.0118216532614,  0.00860719483079,  0.00553245272614,  0.00265435214565
};

/* quick acceptance check: tabulated values for 2^30 times x[i] / x[i+1], 
   used to accept for U * x[i+1] <= x[i] without any floating point 
   operations
*/

static const guint32 z_ktab[ZLEVELS] = {
             0,   805801241,   913449795,   959292101,
     984613382,  1000640624,  1011683932,  1019748961,
    1025894063,  1030729937,  1034633429,  1037849573,
    1040544531,  1042834934,  1044805040,  1046517229,
    1048018667,  1049345670,  1050526658,  1051584182,
    1052536369,  1053397953,  1054181032,  1054895635,
    1055550142,  1056151611,  1056706028,  1057218504,
    1057693425,  1058134583,  1058545264,  1058928337,
    1059286313,  1059621402,  1059935555,  1060230499,
    1060507771,  1060768742,  1061014636,  1061246555,
    1061465485,  1061672318,  1061867859,  1062052837,
    1062227911,  1062393683,  1062550699,  1062699456,
    1062840408,  1062973969,  1063100520,  1063220406,
    1063333944,  1063441425,  1063543115,  1063639256,
    1063730072,  1063815765,  1063896520,  1063972507,
    1064043878,  1064110772,  1064173314,  1064231616,
    1064285778,  1064335887,  1064382020,  1064424241,
    1064462606,  1064497159,  1064527934,  1064554954,
    1064578233,  1064597773,  1064613567,  1064625597,
    1064633833,  1064638235,  1064638750,  1064635313,
    1064627846,  1064616255,  1064600432,  1064580255,
    1064555581,  1064526248,  1064492076,  1064452858,
    1064408364,  1064358335,  1064302479,  1064240469,
    1064171940,  1064096477,  1064013619,  1063922841,
    1063823555,  1063715090,  1063596689,  1063467488,
    1063326500,  1063172593,  1063004460,  1062820593,
    1062619233,  1062398322,  1062155436,  1061887700,
    1061591677,  1061263222,  1060897294,  1060487689,
    1060026690,  1059504563,  1058908848,  1058223320,
    1057426429,  1056488880,  1055369768,  1054010073,
    1052321211,  1050163479,  1047302125,  1043307773,
    1037295293,  1027071672,  1005071925,   990267308
};

/* quick value conversion: tabulated values of 2^{-30} * x[i] */

static const double z_wtab[ZLEVELS] = {
  2.53622366902e-10,  3.37955476897e-10,  3.97259851698e-10,  4.44655509278e-10,
  4.84906285128e-10,  5.20330822255e-10,  5.52248531789e-10,  5.81488551029e-10,
  6.08609212063e-10,  6.34006194944e-10,  6.57971170180e-10,  6.80725976412e-10,
  7.02443700525e-10,  7.23262287545e-10,  7.43293665211e-10,  7.62630058214e-10,
  7.81348477163e-10,  7.99513985044e-10,  8.17182122838e-10,  8.34400743163e-10,
  8.51211418548e-10,  8.67650538408e-10,  8.83750174467e-10,  8.99538771398e-10,
  9.15041703768e-10,  9.30281729472e-10,  9.45279362150e-10,  9.60053179533e-10,
  9.74620080666e-10,  9.88995501967e-10,  1.00319359990e-09,  1.01722740637e-09,
  1.03110896164e-09,  1.04484942866e-09,  1.05845919195e-09,  1.07194794345e-09,
  1.08532475751e-09,  1.09859815657e-09,  1.11177616911e-09,  1.12486638081e-09,
  1.13787598005e-09,  1.15081179843e-09,  1.16368034712e-09,  1.17648784954e-09,
  1.18924027079e-09,  1.20194334448e-09,  1.21460259701e-09,  1.22722336991e-09,
  1.23981084030e-09,  1.25237003980e-09,  1.26490587207e-09,  1.27742312922e-09,
  1.28992650709e-09,  1.30242061978e-09,  1.31491001339e-09,  1.32739917906e-09,
  1.33989256563e-09,  1.35239459175e-09,  1.36490965774e-09,  1.37744215717e-09,
  1.38999648831e-09,  1.40257706548e-09,  1.41518833039e-09,  1.42783476359e-09,
  1.44052089605e-09,  1.45325132095e-09,  1.46603070584e-09,  1.47886380515e-09,
  1.49175547320e-09,  1.50471067776e-09,  1.51773451438e-09,  1.53083222136e-09,
  1.54400919578e-09,  1.55727101046e-09,  1.57062343216e-09,  1.58407244109e-09,
  1.59762425197e-09,  1.61128533679e-09,  1.62506244951e-09,  1.63896265305e-09,
  1.65299334865e-09,  1.66716230821e-09,  1.68147770973e-09,  1.69594817650e-09,
  1.71058282044e-09,  1.72539129015e-09,  1.74038382454e-09,  1.75557131261e-09,
  1.77096536049e-09,  1.78657836679e-09,  1.80242360761e-09,  1.81851533267e-09,
  1.83486887464e-09,  1.85150077365e-09,  1.86842891988e-09,  1.88567271732e-09,
  1.90325327294e-09,  1.92119361579e-09,  1.93951895237e-09,  1.95825696558e-09,
  1.97743816654e-09,  1.99709631119e-09,  2.01726889649e-09,  2.03799775532e-09,
  2.05932977497e-09,  2.08131777125e-09,  2.10402156115e-09,  2.12750929072e-09,
  2.15185909535e-09,  2.17716119811e-09,  2.20352059330e-09,  2.23106052399e-09,
  2.25992705467e-09,  2.29029518328e-09,  2.32237716305e-09,  2.35643407373e-09,
  2.39279230315e-09,  2.43186768824e-09,  2.47420205352e-09,  2.52052071731e-09,
  2.57182738837e-09,  2.62957023755e-09,  2.69595513584e-09,  2.77459811713e-09,
  2.87208672620e-09,  3.00259438883e-09,  3.20774174925e-09,  3.47813812333e-09
};

union wraprand {
    guint32 u;
    gchar c[4];
};

/* Split a 32-bit random value from GLib into 4 octets: each octet
   provides a 7-bit value for indexing into the Ziggurat plus a
   sign bit.  Load a new guint32 when the material is exhausted.
*/

static guint32 gretl_rand_octet (guint32 *sign)
{
    static union wraprand wr;
    static int i;

    if (i == 0) {
	wr.u = g_rand_int(gretl_rand);
	i = 4;
    }

    *sign = wr.c[--i] & 0x80;

    return wr.c[i] & 0x07F;
}

static double ran_normal_ziggurat (void)
{
    guint32 sign, i, j;
    double x, y;

    while (1) {
	j = g_rand_int(gretl_rand);
	i = gretl_rand_octet(&sign);
	j = j >> 2;

	x = j * z_wtab[i];

	if (j < z_ktab[i]) {
	    break;
	}

	if (i < 127) {
	    double y0 = z_ytab[i], y1 = z_ytab[i+1];

	    y = y1 + (y0 - y1) * gretl_rand_01();
	} else {
	    x = ZIG_R - log(1.0 - gretl_rand_01()) / ZIG_R;
	    y = exp(-ZIG_R * (x - 0.5 * ZIG_R)) * gretl_rand_01();
	}

	if (y < exp(-0.5 * x * x)) {
	    break;
	}
    }

    return sign ? x : -x;
}

/* Box-Muller polar method */

static void gretl_two_snormals (double *z1, double *z2) 
{
    double x, y, z;

 tryagain:
    x = 2 * gretl_rand_01() - 1;
    y = 2 * gretl_rand_01() - 1;
    z = x * x + y * y;
    if (z >= 1) {
	goto tryagain;
    }

    z = sqrt(-2 * log(z)/z);
    *z1 = z * x;
    *z2 = z * y;
}

static double ran_normal_box_muller (void)
{
    double x, y, z;

 tryagain:
    x = gretl_rand_01();
    y = gretl_rand_01();
    z = sqrt(-2. * log(x));
    if (isnan(z) || isinf(z)) {
	goto tryagain;
    }

    return z * cos(M_2PI * y);
}

static int use_box_muller;

void gretl_rand_set_box_muller (int s)
{
    use_box_muller = s;
}

int gretl_rand_get_box_muller (void)
{
    return use_box_muller;
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
    
    if (use_box_muller) {
	double z1, z2;

	for (t=t1; t<=t2; t++) {
	    gretl_two_snormals(&z1, &z2);
	    a[t] = z1;
	    if (t < t2) {
		a[++t] = z2;
	    }
	}
    } else {	
	for (t=t1; t<=t2; t++) {
	    a[t] = ran_normal_ziggurat();
	}
    }
}

/**
 * gretl_one_snormal:
 *
 * Returns: a single drawing from the standard normal distribution.
 */

double gretl_one_snormal (void) 
{
    if (use_box_muller) {
	return ran_normal_box_muller();
    } else {
	return ran_normal_ziggurat();
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


