/* mt19937.c */

/* A C-program for MT19937: Integer version (1999/10/28)          */
/*  genrand() generates one pseudorandom unsigned integer (32bit) */
/* which is uniformly distributed among 0 to 2^32-1  for each     */
/* call. sgenrand(seed) sets initial values to the working area   */
/* of 624 words. Before genrand(), sgenrand(seed) must be         */
/* called once. (seed is any 32-bit integer.)                     */
/*   Coded by Takuji Nishimura, considering the suggestions by    */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.              */

/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */
/* 02111-1307  USA                                                 */

/* Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura. */
/* Copyright (C) 2009 Makoto Matsumoto, Takuji Nishimura and       */
/* Mutsuo Saito.                                                   */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */

#include<stdio.h>
#include<inttypes.h>
#include"mt19937.h"

/* Period parameters */
/* #define N 624 */
#define M 397
#define MATRIX_A UINT32_C(0x9908b0df)   /* constant vector a */
#define UPPER_MASK UINT32_C(0x80000000) /* most significant w-r bits */
#define LOWER_MASK UINT32_C(0x7fffffff) /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B UINT32_C(0x9d2c5680)
#define TEMPERING_MASK_C UINT32_C(0xefc60000)
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

/* Initializing the array with a seed */
void _sgenrand_dc(_org_state *st, uint32_t seed)
{
    int i;

    for (i=0;i<N;i++) {
	st->mt[i] = seed;
        seed = (UINT32_C(1812433253) * (seed  ^ (seed >> 30))) + i + 1;
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
    }
    st->mti = N;
}


uint32_t _genrand_dc(_org_state *st)
{
    uint32_t y;
    static const uint32_t mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (st->mti >= N) { /* generate N words at one time */
        int kk;

        for (kk=0;kk<N-M;kk++) {
            y = (st->mt[kk]&UPPER_MASK)|(st->mt[kk+1]&LOWER_MASK);
            st->mt[kk] = st->mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (st->mt[kk]&UPPER_MASK)|(st->mt[kk+1]&LOWER_MASK);
            st->mt[kk] = st->mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (st->mt[N-1]&UPPER_MASK)|(st->mt[0]&LOWER_MASK);
        st->mt[N-1] = st->mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        st->mti = 0;
    }

    y = st->mt[st->mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return y;
}
