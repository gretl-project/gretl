#ifndef SSE2NEON_H
#define SSE2NEON_H

/* A small subset of sse2neon.h by John W. Ratcliff et al,
   released under the MIT license, reproduced below. Should
   be sufficient for useage of SSE2 in SFMT.
*/

/*
 * sse2neon is freely redistributable under the MIT License.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <arm_neon.h>

typedef int64x2_t __m128i; /* 128-bit vector containing integers */

#ifndef likely
# define likely(x) __builtin_expect(!!(x), 1)
#endif
# ifndef unlikely
#define unlikely(x) __builtin_expect(!!(x), 0)
#endif

#define vreinterpretq_m128_u32(x) vreinterpretq_f32_u32(x)
#define vreinterpretq_u32_m128i(x) vreinterpretq_u32_s64(x)
#define vreinterpretq_m128_s8(x) vreinterpretq_f32_s8(x)
#define vreinterpretq_s8_m128i(x) vreinterpretq_s8_s64(x)
#define vreinterpretq_m128i_s32(x) vreinterpretq_s64_s32(x)
#define vreinterpretq_m128i_u32(x) vreinterpretq_s64_u32(x)
#define vreinterpretq_s32_m128i(x) vreinterpretq_s32_s64(x)
#define vreinterpretq_m128i_s8(x) vreinterpretq_s64_s8(x)

#define _mm_srli_epi32(a, imm)                                             \
    __extension__({                                                        \
        __m128i ret;                                                       \
        if (unlikely((imm) == 0)) {                                        \
            ret = a;                                                       \
        }                                                                  \
        if (likely(0 < (imm) && (imm) < 32)) {                             \
            ret = vreinterpretq_m128i_u32(                                 \
                vshlq_u32(vreinterpretq_u32_m128i(a), vdupq_n_s32(-imm))); \
        } else {                                                           \
            ret = _mm_setzero_si128();                                     \
        }                                                                  \
        ret;                                                               \
    })

#define _mm_srli_si128(a, imm)                                              \
    __extension__({                                                         \
        __m128i ret;                                                        \
        if (unlikely((imm) <= 0)) {                                         \
            ret = a;                                                        \
        }                                                                   \
        if (unlikely((imm) > 15)) {                                         \
            ret = _mm_setzero_si128();                                      \
        } else {                                                            \
            ret = vreinterpretq_m128i_s8(                                   \
                vextq_s8(vreinterpretq_s8_m128i(a), vdupq_n_s8(0), (imm))); \
        }                                                                   \
        ret;                                                                \
    })

static inline __m128i _mm_setzero_si128(void)
{
    return vreinterpretq_m128i_s32(vdupq_n_s32(0));
}

static inline __m128i _mm_slli_epi32(__m128i a, int imm)
{
    if (unlikely(imm <= 0)) /* TODO: add constant range macro: [0, 255] */
        return a;
    if (unlikely(imm > 31))
        return _mm_setzero_si128();
    return vreinterpretq_m128i_s32(
        vshlq_s32(vreinterpretq_s32_m128i(a), vdupq_n_s32(imm)));
}

static inline __m128i _mm_xor_si128(__m128i a, __m128i b)
{
    return vreinterpretq_m128i_s32(
        veorq_s32(vreinterpretq_s32_m128i(a), vreinterpretq_s32_m128i(b)));
}

#define _mm_slli_si128(a, imm)                                          \
    __extension__({                                                     \
        __m128i ret;                                                    \
        if (unlikely((imm) <= 0)) {                                     \
            ret = a;                                                    \
        }                                                               \
        if (unlikely((imm) > 15)) {                                     \
            ret = _mm_setzero_si128();                                  \
        } else {                                                        \
            ret = vreinterpretq_m128i_s8(vextq_s8(                      \
                vdupq_n_s8(0), vreinterpretq_s8_m128i(a), 16 - (imm))); \
        }                                                               \
        ret;                                                            \
    })

static inline __m128i _mm_and_si128(__m128i a, __m128i b)
{
    return vreinterpretq_m128i_s32(
        vandq_s32(vreinterpretq_s32_m128i(a), vreinterpretq_s32_m128i(b)));
}

#endif /* SSE2NEON_H */
