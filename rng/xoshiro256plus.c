/*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide.

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */

#include <stdint.h>

/* This is xoshiro256+ 1.0, our best and fastest generator for floating-point
   numbers. We suggest to use its upper bits for floating-point
   generation, as it is slightly faster than xoshiro256++/xoshiro256**. It
   passes all tests we are aware of except for the lowest three bits,
   which might fail linearity tests (and just those), so if low linear
   complexity is not considered an issue (as it is usually the case) it
   can be used to generate 64-bit outputs, too.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s.
*/

static inline uint64_t xor_rotl (const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t xor_state[4];

static void set_xor_state (uint64_t u)
{
    xor_state[0] = smx = u;
    xor_state[1] = splitmix64_next();
    xor_state[2] = splitmix64_next();
    xor_state[3] = splitmix64_next();
}

/* get the next pseudo-random uint64_t */

static uint64_t xor_i64 (void) {
    const uint64_t ret = xor_state[0] + xor_state[3];
    const uint64_t t = xor_state[1] << 17;

    xor_state[2] ^= xor_state[0];
    xor_state[3] ^= xor_state[1];
    xor_state[1] ^= xor_state[2];
    xor_state[0] ^= xor_state[3];

    xor_state[2] ^= t;
    xor_state[3] = xor_rotl(xor_state[3], 45);

    return ret;
}

/* This is the jump function for the generator. It is equivalent to
   2^128 calls to xor_i64(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations.
*/

static void xor_jump (void) {
    static const uint64_t JUMP[] = {
        0x180ec6d33cfd0aba, 0xd5a61266f0c9392c,
        0xa9582618e03fc9aa, 0x39abdc4529b1661c
    };
    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    int i, b;

    for (i = 0; i < sizeof JUMP / sizeof *JUMP; i++) {
        for (b = 0; b < 64; b++) {
            if (JUMP[i] & UINT64_C(1) << b) {
                s0 ^= xor_state[0];
                s1 ^= xor_state[1];
                s2 ^= xor_state[2];
                s3 ^= xor_state[3];
            }
            xor_i64();
        }
    }

    xor_state[0] = s0;
    xor_state[1] = s1;
    xor_state[2] = s2;
    xor_state[3] = s3;
}

#if 0 /* we might want this later? */

/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to xor_next(); it can be used to generate 2^64 starting
   points, from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations.
*/

static void xor_long_jump (void) {
    static const uint64_t LONG_JUMP[] = {
        0x76e15d3efefdcbbf, 0xc5004e441c522fb3,
        0x77710069854ee241, 0x39109bb02acbe635
    };
    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    int i, b;

    for (i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++) {
        for (b = 0; b < 64; b++) {
            if (LONG_JUMP[i] & UINT64_C(1) << b) {
                s0 ^= xor_state[0];
                s1 ^= xor_state[1];
                s2 ^= xor_state[2];
                s3 ^= xor_state[3];
            }
            xor_i64();
        }
    }

    xor_state[0] = s0;
    xor_state[1] = s1;
    xor_state[2] = s2;
    xor_state[3] = s3;
}

#endif
