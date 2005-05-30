/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* random.c for gretl */

#include "libgretl.h"

#include <time.h>

#ifdef USE_GTK2
# define HAVE_G_RAND
#endif

#ifdef HAVE_G_RAND
# include <glib.h>
#else
# include "mt19937ar.c"
#endif

#ifdef HAVE_G_RAND
static GRand *gretl_rand;
#endif

static unsigned int useed;

unsigned int get_gretl_random_seed (void)
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
#ifdef HAVE_G_RAND
    gretl_rand = g_rand_new();
    gretl_rand_set_seed((guint32) useed);
#else
    init_genrand(useed);
#endif
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
#ifdef HAVE_G_RAND
    g_rand_set_seed(gretl_rand, useed);
#else
    init_genrand(useed);
#endif
}

/**
 * gretl_uniform_dist:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the uniform distribution on 0-1, using the Mersenne
 * Twister.
 */

void gretl_uniform_dist (double *a, int t1, int t2) 
{
    int t;

    for (t=t1; t<=t2; t++) {
#ifdef HAVE_G_RAND
	a[t] = g_rand_double_range(gretl_rand, 0.0, 1.0);
#else
	a[t] = genrand_int32() * (1.0 / 4294967296.0);
#endif 
    }
}

/**
 * gretl_normal_dist:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the standard distribution, using the Mersenne Twister for
 * uniform input and the Box-Muller method for converting to the
 * normal distribution.
 */

void gretl_normal_dist (double *a, int t1, int t2) 
{
    int t;
    double xx, yy, zz;

    for (t=t1; t<=t2; t++) {
    tryagain:
#ifdef HAVE_G_RAND
	xx = g_rand_double(gretl_rand);
	yy = g_rand_double(gretl_rand);
#else
	xx = genrand_real2();
	yy = genrand_real2();
#endif
	zz = sqrt(-2. * log(xx));
	if (isnan(zz) || isinf(zz)) goto tryagain;
	a[t] = zz * cos(2. * M_PI * yy);
    }
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
#ifdef HAVE_G_RAND
    return g_rand_int_range(gretl_rand, 0, max);
#else
    return genrand_int32() * (max / 4294967296.0);
#endif
}

/**
 * gretl_rand_free:
 *
 * Free the gretl_rand structure (may be called at program exit).
 */

void gretl_rand_free (void)
{
#ifdef HAVE_G_RAND
    g_rand_free(gretl_rand);
#else
    return;
#endif
}


