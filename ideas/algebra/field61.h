/* Mersenne-61 and Gaussian F_{q^2} arithmetic, q = 2^61 - 1, as small inline helpers for ideas/algebra/.
 * These are copies of the static inline helpers in hash/kp7.c (red, addq, subq, mulq, fadd, fmul, fold), kept here so
 * hash/ stays untouched; the tests check bit-identity against kp7_roll64 wherever the two overlap. */
#ifndef IDEAS_FIELD61_H
#define IDEAS_FIELD61_H
#include <stdint.h>
#include "../../hash/kp7.h"

#define F61_Q ((uint64_t)0x1FFFFFFFFFFFFFFFULL)

static inline uint64_t f61_red(unsigned __int128 p) /* p < 2^124 */
{
  uint64_t r = ((uint64_t)p & F61_Q) + (uint64_t)(p >> 61);
  r = (r & F61_Q) + (r >> 61);
  return r >= F61_Q ? r - F61_Q : r;
}
static inline uint64_t f61_add(uint64_t a, uint64_t b) { uint64_t r = a + b; return r >= F61_Q ? r - F61_Q : r; }
static inline uint64_t f61_sub(uint64_t a, uint64_t b) { return a >= b ? a - b : a + F61_Q - b; }
static inline uint64_t f61_mul(uint64_t a, uint64_t b) { return f61_red((unsigned __int128)a * b); }
static inline uint64_t f61_neg(uint64_t a) { return a ? F61_Q - a : 0; }
static inline uint64_t f61_pow(uint64_t a, uint64_t e)
{
  uint64_t r = 1;
  while (e) { if (e & 1) r = f61_mul(r, a); a = f61_mul(a, a); e >>= 1; }
  return r;
}

static inline kp7_fe g_add(kp7_fe a, kp7_fe b) { kp7_fe r = { f61_add(a.re, b.re), f61_add(a.im, b.im) }; return r; }
static inline kp7_fe g_sub(kp7_fe a, kp7_fe b) { kp7_fe r = { f61_sub(a.re, b.re), f61_sub(a.im, b.im) }; return r; }
static inline kp7_fe g_neg(kp7_fe a) { kp7_fe r = { f61_neg(a.re), f61_neg(a.im) }; return r; }
static inline kp7_fe g_mul(kp7_fe a, kp7_fe b)
{
  const unsigned __int128 Q2 = (unsigned __int128)F61_Q * F61_Q;
  unsigned __int128 ac = (unsigned __int128)a.re * b.re, bd = (unsigned __int128)a.im * b.im;
  unsigned __int128 ad = (unsigned __int128)a.re * b.im, bc = (unsigned __int128)a.im * b.re;
  kp7_fe r = { f61_red(ac + (Q2 - bd)), f61_red(ad + bc) };
  return r;
}
static inline kp7_fe g_mul_real(kp7_fe a, uint64_t x) { kp7_fe r = { f61_mul(a.re, x), f61_mul(a.im, x) }; return r; }
static inline uint64_t g_fold(kp7_fe z) { return z.re ^ ((z.im << 32) | (z.im >> 32)); }

/* A = 1, C = i, G = -i, T = -1; code 0..3 = A C G T, as kp7_base */
static inline kp7_fe g_base(int c)
{
  static const kp7_fe B[4] = { { 1, 0 }, { 0, 1 }, { 0, F61_Q - 1 }, { F61_Q - 1, 0 } };
  return B[c];
}
static inline int base_code(char ch)
{
  switch (ch) { case 'A': case 'a': return 0; case 'C': case 'c': return 1; case 'G': case 'g': return 2;
                case 'T': case 't': return 3; default: return -1; }
}
#endif
