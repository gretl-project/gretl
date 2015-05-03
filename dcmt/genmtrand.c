/* genmtrand.c */

/*
  Copyright (C) 2001-2009 Makoto Matsumoto and Takuji Nishimura.
  Copyright (C) 2009 Mutsuo Saito
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include "dci.h"

#define SHIFT1 18

void sgenrand_mt(uint32_t seed, mt_struct *mts)
{
    int i;

    for (i=0; i<mts->nn; i++) {
	mts->state[i] = seed;
        seed = (UINT32_C(1812433253) * (seed  ^ (seed >> 30))) + i + 1;
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
    }
    mts->i = mts->nn;

    for (i=0; i<mts->nn; i++)
	mts->state[i] &= mts->wmask;
}

uint32_t genrand_mt(mt_struct *mts)
{
    uint32_t *st, uuu, lll, aa, x;
    int k,n,m,lim;

    if ( mts->i >= mts->nn ) {
	n = mts->nn; m = mts->mm;
	aa = mts->aaa;
	st = mts->state;
	uuu = mts->umask; lll = mts->lmask;

	lim = n - m;
	for (k=0; k<lim; k++) {
	    x = (st[k]&uuu)|(st[k+1]&lll);
	    st[k] = st[k+m] ^ (x>>1) ^ (x&1U ? aa : 0U);
	}
	lim = n - 1;
	for (; k<lim; k++) {
	    x = (st[k]&uuu)|(st[k+1]&lll);
	    st[k] = st[k+m-n] ^ (x>>1) ^ (x&1U ? aa : 0U);
	}
	x = (st[n-1]&uuu)|(st[0]&lll);
	st[n-1] = st[m-1] ^ (x>>1) ^ (x&1U ? aa : 0U);
	mts->i=0;
    }

    x = mts->state[mts->i];
    mts->i += 1;
    x ^= x >> mts->shift0;
    x ^= (x << mts->shiftB) & mts->maskB;
    x ^= (x << mts->shiftC) & mts->maskC;
    x ^= x >> mts->shift1;

    return x;
}

