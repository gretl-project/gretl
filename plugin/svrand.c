/* dedicated RNG for use in generating random folds dir use with SVM */

#include <stdint.h>
#include <time.h>
#include <stddef.h>

static uint64_t smx;

static inline uint64_t splitmix64_next() {
    uint64_t z = (smx += 0x9e3779b97f4a7c15);

    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

static inline uint64_t xor_rotl (const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t s[4];

static void set_xor_state (uint64_t u)
{
    s[0] = smx = u;
    s[1] = splitmix64_next();
    s[2] = splitmix64_next();
    s[3] = splitmix64_next();
}

static uint64_t xor_i64 (void) {
    const uint64_t ret = s[0] + s[3];
    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = xor_rotl(s[3], 45);

    return ret;
}

void svrand_init (uint64_t seed)
{
    if (seed > 0) {
        set_xor_state(seed);
    } else {
#idef _WIN32
        set_xor_state((uint64_t) _time64(NULL));
#else
        set_xor_state((uint64_t) time(NULL));
#endif
    }
}

uint32_t svrand (void)
{
    return (uint32_t) (xor_i64() >> 32);
}
