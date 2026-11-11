/* Drawing the hash parameters x and lam for the KP-7 hashes (kp7.c, kp7e.c, kp7m31.c).
 *
 * The collision bound (Schwartz-Zippel over (x, lam)) assumes x and lam are INDEPENDENT and uniform. Each init function
 * therefore takes two separate sources, one for x and one for lam (round-1 adjudication L1 F1: before 2026-11-11 both
 * came from one 64-bit splitmix64 state, so lam was a function of x).
 *
 *   kp7_rng_os()        every word from the operating system (getrandom(2), falling back to /dev/urandom). This is
 *                       the default: an init called with a NULL seed uses it. Field elements are then drawn from fresh
 *                       entropy, so an F_{q^2} element gets its full ~122 bits.
 *   kp7_rng_seeded(s)   deterministic, for tests and benchmarks: a splitmix64 stream started at s ^ tag. The tag
 *                       separates the x stream from the lam stream, so equal seeds still give unrelated values. A
 *                       64-bit seed carries at most 64 bits of entropy, which is less than an F_{q^2} element holds:
 *                       the deterministic mode reproduces runs, it does not meet the bound's premise for F_{q^2}.
 *
 * kp7_uniform(r, m, lo) is uniform on [lo, m) for a Mersenne modulus m = 2^b - 1: it masks a word to b bits, which is
 * uniform on [0, 2^b), and rejects m itself and anything below lo. No modular reduction, so no bias.
 * Header-only (static inline) so that every library and tool shares one definition. */
#ifndef KP7SEED_H
#define KP7SEED_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/random.h>

#define KP7_TAG_X 0x78A5C3D2E1F00B17ULL   /* domain tags: the x stream and the lam stream never coincide */
#define KP7_TAG_LAM 0x1A3B5C7D9E0F2468ULL

typedef struct {
  int os;     /* 1: words from the OS; 0: splitmix64 stream */
  uint64_t s; /* splitmix64 state */
} kp7_rng;

static inline uint64_t kp7_splitmix64(uint64_t* s)
{
  uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

/* 8 bytes of OS entropy; aborts rather than ever return a predictable value */
static inline uint64_t kp7_os_u64(void)
{
  uint64_t v = 0;
  unsigned char* b = (unsigned char*)&v;
  size_t got = 0;
  while (got < sizeof v) {
    ssize_t n = getrandom(b + got, sizeof v - got, 0);
    if (n > 0) { got += (size_t)n; continue; }
    if (n < 0 && errno == EINTR) continue;
    break; /* ENOSYS or similar: try the device */
  }
  if (got < sizeof v) {
    FILE* f = fopen("/dev/urandom", "rb");
    if (!f || fread(b + got, 1, sizeof v - got, f) != sizeof v - got) {
      fprintf(stderr, "kp7seed: no OS randomness (getrandom and /dev/urandom both failed)\n");
      abort();
    }
    fclose(f);
  }
  return v;
}

static inline kp7_rng kp7_rng_os(void) { kp7_rng r = { 1, 0 }; return r; }
static inline kp7_rng kp7_rng_seeded(uint64_t seed, uint64_t tag) { kp7_rng r = { 0, seed ^ tag }; return r; }
/* NULL -> OS randomness, else the deterministic stream for *seed */
static inline kp7_rng kp7_rng_from(const uint64_t* seed, uint64_t tag)
{
  return seed ? kp7_rng_seeded(*seed, tag) : kp7_rng_os();
}

static inline uint64_t kp7_rng_next(kp7_rng* r) { return r->os ? kp7_os_u64() : kp7_splitmix64(&r->s); }

/* uniform on [lo, m), m = 2^b - 1 (all ones), by masking and rejection */
static inline uint64_t kp7_uniform(kp7_rng* r, uint64_t m, uint64_t lo)
{
  for (;;) {
    uint64_t v = kp7_rng_next(r) & m;
    if (v >= lo && v < m) return v;
  }
}
#endif
