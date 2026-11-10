/* Batched reimplementation of ntHash's canonical k-mer hash (ntHash2 code, BirolLab/ntHash at c26bd45), for a fair
 * comparison with kp7's batch APIs: the same rolling recurrence in a loop with the state in locals, emitting blocks.
 * It must be bit-identical to nthash::NtHash::roll()'s hashes()[0] (bench_pack checks the XOR, count and positions).
 *
 * srol/sror and the seeds are copied from ntHash src/internal.hpp, MIT licence:
 *   Copyright (c) 2018 Hamid Mohamadi. Permission is hereby granted, free of charge, to any person obtaining a copy of
 *   this software ... THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND. (full text: ntHash/LICENSE)
 * The initial hashes use srol's linearity (it permutes bits, so srol(a ^ b) = srol(a) ^ srol(b)) instead of ntHash's
 * tetramer tables:  fwd = XOR_i srol^{k-1-i}(seed[s_i]),  rev = XOR_i srol^{i}(seed[s_i & 7])  (the complement seed).
 */
#ifndef NTBATCH_H
#define NTBATCH_H
#include <stddef.h>
#include <stdint.h>
#include <initializer_list>

namespace ntbatch {

static const uint64_t SEED_A = 0x3c8bfbb395c60474ULL, SEED_C = 0x3193c18562a02b4cULL, SEED_G = 0x20323ed082572324ULL,
                      SEED_T = 0x295549f54be24456ULL;

inline uint64_t srol(uint64_t x)
{
  uint64_t m = ((x & 0x8000000000000000ULL) >> 30) | ((x & 0x100000000ULL) >> 32);
  return ((x << 1) & 0xFFFFFFFDFFFFFFFFULL) | m;
}
inline uint64_t srol(uint64_t x, unsigned d)
{
  if (d == 0) return x;
  uint64_t v = (x << d) | (x >> (64 - d));
  uint64_t y = (v ^ (v >> 33)) & (~0ULL >> (64 - d));
  return v ^ (y | (y << 33));
}
inline uint64_t sror(uint64_t x)
{
  uint64_t m = ((x & 0x200000000ULL) << 30) | ((x & 1ULL) << 32);
  return ((x >> 1) & 0xFFFFFFFEFFFFFFFFULL) | m;
}

struct Tables {
  uint64_t seed[256];  /* ntHash's SEED_TAB: ACGTU/acgtu, complements at [c & 7], 0 = invalid */
  uint64_t rolk[256];  /* srol(seed[c], k) for the outgoing forward base */
  uint64_t rolkc[256]; /* srol(seed[c & 7], k) for the incoming reverse base */
  unsigned k;
  explicit Tables(unsigned kk) : k(kk)
  {
    for (int i = 0; i < 256; i++) seed[i] = 0;
    seed[1] = SEED_T; seed[3] = SEED_G; seed[4] = SEED_A; seed[5] = SEED_A; seed[7] = SEED_C;
    for (int c : { 'A', 'a' }) seed[c] = SEED_A;
    for (int c : { 'C', 'c' }) seed[c] = SEED_C;
    for (int c : { 'G', 'g' }) seed[c] = SEED_G;
    for (int c : { 'T', 't', 'U', 'u' }) seed[c] = SEED_T;
    for (int i = 0; i < 256; i++) { rolk[i] = srol(seed[i], k); rolkc[i] = srol(seed[i & 7], k); }
  }
};

typedef void (*emit_fn)(void* ctx, size_t start, const uint64_t* h, size_t n);

/* canonical hash (fwd + rev) of every valid k-mer, in blocks of at most cap consecutive k-mers */
inline size_t hash_seq(const Tables& T, const char* seq, size_t len, uint64_t* buf, size_t cap, emit_fn emit, void* ctx)
{
  const unsigned k = T.k;
  size_t total = 0, i = 0;
  while (i < len) {
    while (i < len && !T.seed[(unsigned char)seq[i]]) i++;
    size_t a = i;
    while (i < len && T.seed[(unsigned char)seq[i]]) i++;
    if (i - a < k) continue;
    const size_t nk = i - a - k + 1;
    const unsigned char* s = (const unsigned char*)seq + a;
    uint64_t f = 0, r = 0;
    for (unsigned t = 0; t < k; t++) f = srol(f) ^ T.seed[s[t]];
    for (int t = (int)k - 1; t >= 0; t--) r = srol(r) ^ T.seed[s[t] & 7];
    size_t n = 0, start = a;
    buf[n++] = f + r;
    for (size_t q = 0; q + 1 < nk; q++) {
      f = srol(f) ^ T.seed[s[q + k]] ^ T.rolk[s[q]];
      r = sror(r ^ T.rolkc[s[q + k]] ^ T.seed[s[q] & 7]);
      if (n == cap) { emit(ctx, start, buf, n); start += n; n = 0; }
      buf[n++] = f + r;
    }
    emit(ctx, start, buf, n);
    total += nk;
  }
  return total;
}

} // namespace ntbatch
#endif
