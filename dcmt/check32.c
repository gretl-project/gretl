/* check32.c */

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
#include "mt19937.h"

#define LSB 0x1
#define WORDLEN 32

#if defined(DEBUG)
/* period parameters */
#define NNN 17
#define MMM 9
#define WWW 32
#define RRR 23
#define NREPEAT 1000
/********************/
static uint32_t next_A(check32_t *ck, _org_state *st, int w);

int main(int argc, char **argv)
{
    int i;
    uint32_t a;
    _org_state st;
    check32_t ck;

    _sgenrand_dc(&st, 3249);

    for(i=0; i<NREPEAT; i++) {
	a = next_A(&ck, &st, WWW);
	if (IRRED == _CheckPeriod_dc(&ck, &st, a,MMM,NNN,RRR,WWW)) {
	    printf ("%x\n",a);
	}
    }

    return 0;
}
#endif

void _InitCheck32_dc(check32_t *ck, int r, int w)
{
    int i;

    /* word_mask (least significant w bits) */
    ck->word_mask = 0xFFFFFFFF;
    ck->word_mask <<= WORDLEN - w;
    ck->word_mask >>= WORDLEN - w;
    /* lower_mask (least significant r bits) */
    for (ck->lower_mask=0,i=0; i<r; ++i) {
	ck->lower_mask <<= 1;
	ck->lower_mask |= LSB;
    }
    /* upper_mask (most significant (w-r) bits */
    ck->upper_mask = (~ck->lower_mask) & ck->word_mask;
}

#if defined(DEBUG)
static uint32_t next_A(check32_t *ck, _org_state *st, int w)
{
    uint32_t x;

    x = _genrand_dc(st);

    x &= ck->word_mask;
    x |= (LSB << (w-1));

/*  printf ("AAA %8x\n", x);
  getchar();
*/
    return x;
}
#endif

int _CheckPeriod_dc(check32_t *ck, _org_state *st,
		    uint32_t a, int m, int n, int r, int w)
{
    int i, j, p, pp;
    uint32_t y, *x, *init, mat[2];


    p = n*w-r;
    x = (uint32_t*) malloc (2*p*sizeof(uint32_t));
    if (NULL==x) {
	printf("malloc error in \"_CheckPeriod_dc()\"\n");
	exit(1);
    }

    init = (uint32_t*) malloc (n*sizeof(uint32_t));
    if (NULL==init) {
	printf("malloc error \"_CheckPeriod_dc()\"\n");
	free(x);
	exit(1);
    }

    /* set initial values */
    for (i=0; i<n; ++i)
	x[i] = init[i] = (ck->word_mask & _genrand_dc(st));
    /* it is better that LSBs of x[2] and x[3] are different */
    if ( (x[2]&LSB) == (x[3]&LSB) ) {
	x[3] ^= 1;
	init[3] ^= 1;
    }

    pp = 2*p-n;
    mat[0] = 0; mat[1] = a;
    for (j=0; j<p; ++j) {

	/* generate */
	for (i=0; i<pp; ++i){
	    y = (x[i]&ck->upper_mask) | (x[i+1]&ck->lower_mask);
	    x[i+n] = x[i+m] ^ ( (y>>1) ^ mat[y&LSB] );
	}

	/* pick up odd subscritpt elements */
	for (i=2; i<=p; ++i)
	    x[i] = x[(i<<1)-1];

	/* reverse generate */
	for (i=p-n; i>=0; --i) {
	    y = x[i+n] ^ x[i+m] ^ mat[ x[i+1]&LSB ];
	    y <<=1; y |= x[i+1]&LSB;

	    x[i+1] = (x[i+1]&ck->upper_mask) | (y&ck->lower_mask);
	    x[i] = (y&ck->upper_mask) | (x[i]&ck->lower_mask);
	}

    }

    if ((x[0]&ck->upper_mask)==(init[0]&ck->upper_mask)) {
	for (i=1; i<n; ++i) {
	    if (x[i] != init[i])
		break;
	}
	if (i==n) {
	    free(x); free(init);
	    return IRRED;
	}
    }


    free(x); free(init);
    return REDU;
}
