/* Cheaper evaluation of kept k-mers; see kp7kept.h. */
#include "kp7kept.h"
#include <stdlib.h>
#include <string.h>
#ifdef __BMI2__
#include <immintrin.h>
#endif

/* ---- codes and gate: identical to kp7pack.c ---- */
static int8_t PCODE[256];
static int pcode_ready = 0;
static void pcode_init(void)
{
  for (int i = 0; i < 256; i++) PCODE[i] = -1;
  PCODE['A'] = PCODE['a'] = 0; PCODE['C'] = PCODE['c'] = 1;
  PCODE['G'] = PCODE['g'] = 2; PCODE['T'] = PCODE['t'] = 3;
  pcode_ready = 1;
}
static inline uint64_t mix(uint64_t x)
{
  x ^= x >> 32;
  x *= 0xd6e8feb86659fd93ULL;
  return x ^ (x >> 32);
}
static inline uint64_t gate(uint64_t f, uint64_t r) { return mix(f) + mix(r); }
uint64_t kp7kept_gate(uint64_t fwd, uint64_t rc) { return gate(fwd, rc); }

static inline uint64_t pext64(uint64_t v, uint64_t m)
{
#ifdef __BMI2__
  return _pext_u64(v, m);
#else
  uint64_t r = 0;
  for (uint64_t bb = 1; m; bb <<= 1, m &= m - 1)
    if (v & m & -m) r |= bb;
  return r;
#endif
}
#define EVEN 0x5555555555555555ULL
#define ODD 0xAAAAAAAAAAAAAAAAULL

/* ---- F_q, q = 2^61 - 1 ---- */
#define Q KP7E_Q
static inline uint64_t red(unsigned __int128 p) /* p < 2^124 -> [0, q) */
{
  uint64_t r = ((uint64_t)p & Q) + (uint64_t)(p >> 61);
  r = (r & Q) + (r >> 61);
  return r >= Q ? r - Q : r;
}
static inline uint64_t part(uint64_t s) { return (s & Q) + (s >> 61); } /* s < 2^64 -> < q + 8 */
static inline uint64_t mulq(uint64_t a, uint64_t b) { return red((unsigned __int128)a * b); } /* a, b < 2^62 */
static inline uint64_t addq(uint64_t a, uint64_t b) { uint64_t r = a + b; return r >= Q ? r - Q : r; }
static inline uint64_t subq(uint64_t a, uint64_t b) { return a >= b ? a - b : a + Q - b; }

/* Gaussian F_{q^2} product with parts < 4q, folded as kp7_fold: identical to kp7pack.c's fq_sym_fold */
static inline uint64_t fq_red_wide(unsigned __int128 v)
{
  unsigned __int128 t = (v & Q) + (v >> 61);
  uint64_t r = (uint64_t)(t & Q) + (uint64_t)(t >> 61);
  r = (r & Q) + (r >> 61);
  return r >= Q ? r - Q : r;
}
static inline uint64_t fq_sym_fold(uint64_t a, uint64_t b, uint64_t c, uint64_t d)
{
  const unsigned __int128 Q16 = (unsigned __int128)16 * Q * Q;
  unsigned __int128 ac = (unsigned __int128)a * c, bd = (unsigned __int128)b * d;
  unsigned __int128 ad = (unsigned __int128)a * d, bc = (unsigned __int128)b * c;
  uint64_t re = fq_red_wide(ac + (Q16 - bd)), im = fq_red_wide(ad + bc);
  return re ^ ((im << 32) | (im >> 32));
}

/* ---- CM31, p = 2^31 - 1: identical to kp7pack.c's m31_red64 / m31_sym_pack ---- */
#define MP ((uint64_t)KP7M31_P)
static inline uint32_t m31_red64(uint64_t v) /* v < 2^63 -> [0, p) */
{
  uint64_t r = (v & MP) + (v >> 31);
  r = (r & MP) + (r >> 31);
  return (uint32_t)(r >= MP ? r - MP : r);
}
static inline uint64_t m31_sym_pack(uint32_t a, uint32_t b, uint32_t c, uint32_t d)
{
  const uint64_t PP = MP * MP;
  uint64_t re = m31_red64((uint64_t)a * c + (PP - (uint64_t)b * d));
  uint64_t im = m31_red64((uint64_t)a * d + (uint64_t)b * c);
  return re | (im << 31);
}
static inline uint64_t mpk(m31c z) { return (uint64_t)z.re | ((uint64_t)z.im << 32); }
static inline m31c munpk(uint64_t v) { m31c z = { (uint32_t)(v & 0xFFFFFFFFu), (uint32_t)(v >> 32) }; return z; }
static inline m31c msub(m31c a, m31c b)
{
  m31c r = { a.re >= b.re ? a.re - b.re : a.re + KP7M31_P - b.re, a.im >= b.im ? a.im - b.im : a.im + KP7M31_P - b.im };
  return r;
}

/* ---- tables ---- */
static int chunk_ok(unsigned B) { return B == 8 || B == 10 || B == 11 || B == 12 || B == 13 || B == 16; }
static int planes_ok(unsigned B) { return B == 8 || B == 11 || B == 16; }

static int alloc_tab(kp7kept* t, int kind, unsigned k, unsigned B, unsigned nbits)
{
  memset(t, 0, sizeof *t);
  if (k < 1 || k > 32) return 0;
  t->kind = kind;
  t->k = k;
  t->B = B;
  t->nch = (nbits + B - 1) / B;
  t->batch = 512;
  size_t n = (size_t)t->nch << B;
  t->T = (uint64_t*)aligned_alloc(64, n * sizeof(uint64_t));
  if (!t->T) return 0;
  memset(t->T, 0, n * sizeof(uint64_t));
  return 1;
}

/* T[j][v] = c0 (j == 0) + sum of wt[bit] over the set bits of v at global bit positions B*j + u, in F_q */
static void build_q(kp7kept* t, const uint64_t* wt, unsigned nbits, uint64_t c0)
{
  const unsigned B = t->B;
  for (unsigned j = 0; j < t->nch; j++) {
    uint64_t* T = t->T + ((size_t)j << B);
    T[0] = j == 0 ? c0 : 0;
    for (uint64_t v = 1; v < (1ULL << B); v++) {
      unsigned u = (unsigned)__builtin_ctzll(v), g = B * j + u;
      T[v] = addq(T[v & (v - 1)], g < nbits ? wt[g] : 0);
    }
  }
}
/* the same in CM31, entries packed re | im << 32 */
static void build_m(kp7kept* t, const m31c* wt, unsigned nbits, m31c c0)
{
  const unsigned B = t->B;
  const m31c zero = { 0, 0 };
  for (unsigned j = 0; j < t->nch; j++) {
    uint64_t* T = t->T + ((size_t)j << B);
    T[0] = mpk(j == 0 ? c0 : zero);
    for (uint64_t v = 1; v < (1ULL << B); v++) {
      unsigned u = (unsigned)__builtin_ctzll(v), g = B * j + u;
      T[v] = mpk(kp7m31_add(munpk(T[v & (v - 1)]), g < nbits ? wt[g] : zero));
    }
  }
}

int kp7kept_init_e_chunk(kp7kept* t, const kp7e_params* p, unsigned B)
{
  if (!chunk_ok(B) || !alloc_tab(t, KP7KEPT_E_CHUNK, p->k, B, 64)) return 0;
  /* per code bit: bit 2t has weight (val0 - val1) x^t, bit 2t+1 (val0 - val2) x^t; constant lam - val0 X */
  const uint64_t al = subq(p->val[0], p->val[1]), be = subq(p->val[0], p->val[2]);
  uint64_t wt[64] = { 0 }, xt = 1, X = 0;
  for (unsigned s = 0; s < p->k; s++) {
    wt[2 * s] = mulq(al, xt);
    wt[2 * s + 1] = mulq(be, xt);
    X = addq(X, xt);
    xt = mulq(xt, p->x);
  }
#ifdef KP7KEPT_MUTATE_CONST /* test-only: the constant left out */
  build_q(t, wt, 2 * p->k, 0);
#else
  build_q(t, wt, 2 * p->k, subq(p->lam, mulq(p->val[0], X)));
#endif
  return 1;
}

int kp7kept_init_e_planes(kp7kept* t, const kp7e_params* p, unsigned B)
{
  /* alpha P0 + beta P1 = P0 + P1 + w (P1 - P0) needs val = 1, w, -w, -1 */
  if (p->val[0] != 1 || p->val[1] != p->w || p->val[2] != subq(0, p->w)) return 0;
  if (!planes_ok(B) || !alloc_tab(t, KP7KEPT_E_PLANES, p->k, B, 32)) return 0;
  uint64_t wt[32] = { 0 }, xt = 1, X = 0;
  for (unsigned s = 0; s < p->k; s++) { wt[s] = xt; X = addq(X, xt); xt = mulq(xt, p->x); }
  build_q(t, wt, p->k, 0);
  t->c = subq(p->lam, X);
  t->w = p->w;
  return 1;
}

int kp7kept_init_g_planes(kp7kept* t, const kp7_params* p, unsigned B)
{
  if (p->variant != KP7_XREAL || p->x.im != 0) return 0;
  /* Gaussian val = 1, i, -i, -1 (kp7_base): alpha = 1 - i, beta = 1 + i */
  kp7_fe v1 = kp7_base(1), v2 = kp7_base(2), v0 = kp7_base(0);
  if (v0.re != 1 || v0.im != 0 || v1.re != 0 || v1.im != 1 || v2.re != 0 || v2.im != Q - 1) return 0;
  if (!planes_ok(B) || !alloc_tab(t, KP7KEPT_G_PLANES, p->k, B, 32)) return 0;
  uint64_t wt[32] = { 0 }, xt = 1, X = 0;
  for (unsigned s = 0; s < p->k; s++) { wt[s] = xt; X = addq(X, xt); xt = mulq(xt, p->x.re); }
  build_q(t, wt, p->k, 0);
  t->cre = subq(p->lam.re, X);
  t->cim = p->lam.im;
  return 1;
}

static int m_alpha_beta(const kp7m31_params* p, m31c* al, m31c* be)
{
  /* kp7m31's val[] is indexed A C T G; packed codes are A0 C1 G2 T3 */
  const m31c v0 = p->val[0], v1 = p->val[1], v2 = p->val[3];
  *al = msub(v0, v1);
  *be = msub(v0, v2);
  return v0.re == 1 && v0.im == 0;
}

int kp7kept_init_m_planes(kp7kept* t, const kp7m31_params* p, unsigned B)
{
  m31c al, be;
  if (!m_alpha_beta(p, &al, &be) || al.re != 1 || al.im != KP7M31_P - 1 || be.re != 1 || be.im != 1) return 0;
  if (!planes_ok(B) || !alloc_tab(t, KP7KEPT_M_PLANES, p->k, B, 32)) return 0;
  m31c wt[32], xt = { 1, 0 }, X = { 0, 0 };
  memset(wt, 0, sizeof wt);
  for (unsigned s = 0; s < p->k; s++) { wt[s] = xt; X = kp7m31_add(X, xt); xt = kp7m31_mul(xt, p->x); }
  m31c zero = { 0, 0 };
  build_m(t, wt, p->k, zero);
  t->c = mpk(msub(p->lam, X));
  return 1;
}

int kp7kept_init_m_chunk(kp7kept* t, const kp7m31_params* p, unsigned B)
{
  m31c al, be;
  if (!m_alpha_beta(p, &al, &be) || !chunk_ok(B) || !alloc_tab(t, KP7KEPT_M_CHUNK, p->k, B, 64)) return 0;
  m31c wt[64], xt = { 1, 0 }, X = { 0, 0 };
  memset(wt, 0, sizeof wt);
  for (unsigned s = 0; s < p->k; s++) {
    wt[2 * s] = kp7m31_mul(al, xt);
    wt[2 * s + 1] = kp7m31_mul(be, xt);
    X = kp7m31_add(X, xt);
    xt = kp7m31_mul(xt, p->x);
  }
  build_m(t, wt, 2 * p->k, msub(p->lam, X)); /* val0 = 1 */
  return 1;
}

/* ntHash's split rotation (ntHash src/internal.hpp, MIT; as ntbatch.h) */
static inline uint64_t srol1(uint64_t x)
{
  uint64_t m = ((x & 0x8000000000000000ULL) >> 30) | ((x & 0x100000000ULL) >> 32);
  return ((x << 1) & 0xFFFFFFFDFFFFFFFFULL) | m;
}
int kp7kept_init_nt(kp7kept* t, unsigned k, unsigned B)
{
  if ((B != 8 && B != 16) || !alloc_tab(t, KP7KEPT_NT_CHUNK, k, B, 64)) return 0;
  /* seeds by packed code A0 C1 G2 T3 (ntHash's SEED_A, _C, _G, _T); digit p contributes srol^p(seed[digit]) */
  const uint64_t seed[4] = { 0x3c8bfbb395c60474ULL, 0x3193c18562a02b4cULL, 0x20323ed082572324ULL,
                             0x295549f54be24456ULL };
  uint64_t rot[32][4];
  for (int d = 0; d < 4; d++) {
    uint64_t s = seed[d];
    for (unsigned q = 0; q < 32; q++) { rot[q][d] = s; s = srol1(s); }
  }
  const unsigned bpc = B / 2;
  for (unsigned j = 0; j < t->nch; j++)
    for (uint64_t v = 0; v < (1ULL << B); v++) {
      uint64_t h = 0;
      for (unsigned u = 0; u < bpc; u++) {
        unsigned q = bpc * j + u;
        if (q < k) h ^= rot[q][(v >> (2 * u)) & 3];
      }
      t->T[((size_t)j << B) | v] = h;
    }
  return 1;
}

int kp7kept_init_uh(kp7kept* t, unsigned k, uint64_t seed)
{
  memset(t, 0, sizeof *t);
  if (k < 1 || k > 32) return 0;
  t->kind = KP7KEPT_UH;
  t->k = k;
  t->batch = 512;
  uint64_t s = seed;
  for (int i = 0; i < 2; i++) { /* splitmix64 */
    uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    if (i == 0) t->a_lo = z | 1; else t->a_hi = z;
  }
  return 1;
}

void kp7kept_free(kp7kept* t)
{
  free(t->T);
  t->T = NULL;
}

/* ---- evaluators, B a compile-time constant after inlining ---- */
#define NCH(nbits, B) (((nbits) + (B) - 1) / (B))
#define LK(T, B, j, v) (T)[((size_t)(j) << (B)) | (((v) >> ((B) * (j))) & ((1ULL << (B)) - 1))]

static inline __attribute__((always_inline)) uint64_t ev_e_chunk(const kp7kept* t, const unsigned B, uint64_t fwd,
                                                                 uint64_t rc)
{
  const uint64_t* T = t->T;
  uint64_t g = 0, h = 0; /* <= 8 canonical terms: < 8q < 2^64 */
#pragma GCC unroll 8
  for (unsigned j = 0; j < NCH(64, B); j++) { g += LK(T, B, j, rc); h += LK(T, B, j, fwd); }
  return mulq(part(g), part(h));
}

static inline __attribute__((always_inline)) uint64_t ev_e_planes(const kp7kept* t, const unsigned B, uint64_t fwd,
                                                                  uint64_t rc)
{
  const uint64_t* T = t->T;
  const uint64_t r0 = pext64(rc, EVEN), r1 = pext64(rc, ODD), f0 = pext64(fwd, EVEN), f1 = pext64(fwd, ODD);
  uint64_t a0 = 0, a1 = 0, b0 = 0, b1 = 0; /* <= 4 canonical terms each: < 4q */
#pragma GCC unroll 4
  for (unsigned j = 0; j < NCH(32, B); j++) {
    a0 += LK(T, B, j, r0); a1 += LK(T, B, j, r1); b0 += LK(T, B, j, f0); b1 += LK(T, B, j, f1);
  }
  /* G = c + P0 + P1 + w (P1 - P0): every intermediate < 2^64 and every multiply input < 2^62 */
#ifdef KP7KEPT_MUTATE_PLANE /* test-only: the w term's sign flipped */
  const uint64_t g = t->c + part(a0 + a1) + mulq(t->w, part(a0 + 4 * Q - a1));
#else
  const uint64_t g = t->c + part(a0 + a1) + mulq(t->w, part(a1 + 4 * Q - a0));
#endif
  const uint64_t h = t->c + part(b0 + b1) + mulq(t->w, part(b1 + 4 * Q - b0));
  return mulq(part(g), part(h));
}

static inline __attribute__((always_inline)) uint64_t ev_g_planes(const kp7kept* t, const unsigned B, uint64_t fwd,
                                                                  uint64_t rc)
{
  const uint64_t* T = t->T;
  const uint64_t r0 = pext64(rc, EVEN), r1 = pext64(rc, ODD), f0 = pext64(fwd, EVEN), f1 = pext64(fwd, ODD);
  uint64_t a0 = 0, a1 = 0, b0 = 0, b1 = 0;
#pragma GCC unroll 4
  for (unsigned j = 0; j < NCH(32, B); j++) {
    a0 += LK(T, B, j, r0); a1 += LK(T, B, j, r1); b0 += LK(T, B, j, f0); b1 += LK(T, B, j, f1);
  }
  /* (1 - i) P0 + (1 + i) P1 with P real: re P0 + P1, im P1 - P0; parts < 3q for fq_sym_fold (< 4q) */
  return fq_sym_fold(t->cre + part(a0 + a1), t->cim + part(a1 + 4 * Q - a0), t->cre + part(b0 + b1),
                     t->cim + part(b1 + 4 * Q - b0));
}

static inline __attribute__((always_inline)) uint64_t ev_m_planes(const kp7kept* t, const unsigned B, uint64_t fwd,
                                                                  uint64_t rc)
{
  const uint64_t* T = t->T;
  const uint64_t r0 = pext64(rc, EVEN), r1 = pext64(rc, ODD), f0 = pext64(fwd, EVEN), f1 = pext64(fwd, ODD);
  uint64_t ar0 = 0, ai0 = 0, ar1 = 0, ai1 = 0, br0 = 0, bi0 = 0, br1 = 0, bi1 = 0; /* <= 4 terms < p each */
#pragma GCC unroll 4
  for (unsigned j = 0; j < NCH(32, B); j++) {
    const uint64_t e0 = LK(T, B, j, r0), e1 = LK(T, B, j, r1), d0 = LK(T, B, j, f0), d1 = LK(T, B, j, f1);
    ar0 += e0 & 0xFFFFFFFFu; ai0 += e0 >> 32; ar1 += e1 & 0xFFFFFFFFu; ai1 += e1 >> 32;
    br0 += d0 & 0xFFFFFFFFu; bi0 += d0 >> 32; br1 += d1 & 0xFFFFFFFFu; bi1 += d1 >> 32;
  }
  /* (1 - i)(a + bi) = (a + b) + (b - a) i,  (1 + i)(a + bi) = (a - b) + (a + b) i; sums < 17p */
  const uint64_t cr = t->c & 0xFFFFFFFFu, ci = t->c >> 32, P4 = 4 * MP;
  const uint32_t gr = m31_red64(cr + ar0 + ai0 + ar1 + (P4 - ai1)), gi = m31_red64(ci + ai0 + (P4 - ar0) + ar1 + ai1);
  const uint32_t hr = m31_red64(cr + br0 + bi0 + br1 + (P4 - bi1)), hi = m31_red64(ci + bi0 + (P4 - br0) + br1 + bi1);
  return m31_sym_pack(gr, gi, hr, hi);
}

static inline __attribute__((always_inline)) uint64_t ev_m_chunk(const kp7kept* t, const unsigned B, uint64_t fwd,
                                                                 uint64_t rc)
{
  const uint64_t* T = t->T;
  uint64_t gr = 0, gi = 0, hr = 0, hi = 0; /* <= 8 terms < p each */
#pragma GCC unroll 8
  for (unsigned j = 0; j < NCH(64, B); j++) {
    const uint64_t a = LK(T, B, j, rc), b = LK(T, B, j, fwd);
    gr += a & 0xFFFFFFFFu; gi += a >> 32; hr += b & 0xFFFFFFFFu; hi += b >> 32;
  }
  return m31_sym_pack(m31_red64(gr), m31_red64(gi), m31_red64(hr), m31_red64(hi));
}

static inline __attribute__((always_inline)) uint64_t ev_nt(const kp7kept* t, const unsigned B, uint64_t fwd,
                                                            uint64_t rc)
{
  const uint64_t* T = t->T;
  uint64_t f = 0, r = 0;
#pragma GCC unroll 8
  for (unsigned j = 0; j < NCH(64, B); j++) { f ^= LK(T, B, j, fwd); r ^= LK(T, B, j, rc); }
  return f + r;
}

static inline __attribute__((always_inline)) uint64_t ev_uh(const kp7kept* t, const unsigned B, uint64_t fwd,
                                                            uint64_t rc)
{
  (void)B;
  const uint64_t c = fwd < rc ? fwd : rc;
  return (uint64_t)(((unsigned __int128)t->a_lo * c) >> 64) + t->a_hi * c;
}

/* ---- samplers ---- */
#define RC_SHIFT(k) (2 * ((k) - 1))
#define FOR_KMERS(...)                                                                                                \
  const uint64_t mask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);                                                   \
  const unsigned rsh = RC_SHIFT(k);                                                                                  \
  size_t i = 0;                                                                                                      \
  while (i < len) {                                                                                                  \
    uint64_t fwd = 0, rc = 0;                                                                                        \
    unsigned have = 0;                                                                                               \
    for (; i < len; i++) {                                                                                           \
      int c = PCODE[(unsigned char)seq[i]];                                                                          \
      if (c < 0) { have = 0; fwd = rc = 0; continue; }                                                               \
      fwd = ((fwd << 2) | (uint64_t)c) & mask;                                                                       \
      rc = (rc >> 2) | ((uint64_t)(3 - c) << rsh);                                                                   \
      if (have < k) have++;                                                                                          \
      if (have == k) { size_t pos = i + 1 - k; __VA_ARGS__ }                                                         \
    }                                                                                                                \
  }

#define DEF_SAMPLERS(NAME, EV, B)                                                                                     \
  static size_t NAME##_bat(const kp7kept* t, const char* seq, size_t len, uint64_t thr, size_t* pos_buf,             \
                           uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx)                         \
  {                                                                                                                  \
    const unsigned k = t->k;                                                                                         \
    size_t bmax = t->batch < cap ? t->batch : cap;                                                                   \
    if (bmax > KP7KEPT_BATCH_MAX) bmax = KP7KEPT_BATCH_MAX;                                                          \
    if (bmax < 1) bmax = 1;                                                                                          \
    uint64_t bf[KP7KEPT_BATCH_MAX], br[KP7KEPT_BATCH_MAX];                                                           \
    size_t bp[KP7KEPT_BATCH_MAX], nb = 0, n = 0, kept = 0;                                                           \
    FOR_KMERS({                                                                                                      \
      bf[nb] = fwd; br[nb] = rc; bp[nb] = pos;                                                                       \
      nb += gate(fwd, rc) < thr;                                                                                     \
      if (nb == bmax) {                                                                                              \
        for (size_t q = 0; q < nb; q++) {                                                                            \
          if (n == cap) { emit(ctx, pos_buf, hash_buf, n); n = 0; }                                                  \
          pos_buf[n] = bp[q]; hash_buf[n] = EV(t, B, bf[q], br[q]); n++;                                             \
        }                                                                                                            \
        kept += nb; nb = 0;                                                                                          \
      }                                                                                                              \
    })                                                                                                               \
    for (size_t q = 0; q < nb; q++) {                                                                                \
      if (n == cap) { emit(ctx, pos_buf, hash_buf, n); n = 0; }                                                      \
      pos_buf[n] = bp[q]; hash_buf[n] = EV(t, B, bf[q], br[q]); n++;                                                 \
    }                                                                                                                \
    kept += nb;                                                                                                      \
    if (n) emit(ctx, pos_buf, hash_buf, n);                                                                          \
    return kept;                                                                                                     \
  }                                                                                                                  \
  static size_t NAME##_br(const kp7kept* t, const char* seq, size_t len, uint64_t thr, size_t* pos_buf,              \
                          uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx)                          \
  {                                                                                                                  \
    const unsigned k = t->k;                                                                                         \
    size_t n = 0, kept = 0;                                                                                          \
    FOR_KMERS({                                                                                                      \
      if (gate(fwd, rc) < thr) {                                                                                     \
        if (n == cap) { emit(ctx, pos_buf, hash_buf, n); n = 0; }                                                    \
        pos_buf[n] = pos; hash_buf[n] = EV(t, B, fwd, rc); n++; kept++;                                              \
      }                                                                                                              \
    })                                                                                                               \
    if (n) emit(ctx, pos_buf, hash_buf, n);                                                                          \
    return kept;                                                                                                     \
  }

DEF_SAMPLERS(ec8, ev_e_chunk, 8)
DEF_SAMPLERS(ec10, ev_e_chunk, 10)
DEF_SAMPLERS(ec11, ev_e_chunk, 11)
DEF_SAMPLERS(ec12, ev_e_chunk, 12)
DEF_SAMPLERS(ec13, ev_e_chunk, 13)
DEF_SAMPLERS(ec16, ev_e_chunk, 16)
DEF_SAMPLERS(ep8, ev_e_planes, 8)
DEF_SAMPLERS(ep11, ev_e_planes, 11)
DEF_SAMPLERS(ep16, ev_e_planes, 16)
DEF_SAMPLERS(gp8, ev_g_planes, 8)
DEF_SAMPLERS(gp11, ev_g_planes, 11)
DEF_SAMPLERS(gp16, ev_g_planes, 16)
DEF_SAMPLERS(mp8, ev_m_planes, 8)
DEF_SAMPLERS(mp11, ev_m_planes, 11)
DEF_SAMPLERS(mp16, ev_m_planes, 16)
DEF_SAMPLERS(mc8, ev_m_chunk, 8)
DEF_SAMPLERS(mc10, ev_m_chunk, 10)
DEF_SAMPLERS(mc11, ev_m_chunk, 11)
DEF_SAMPLERS(mc12, ev_m_chunk, 12)
DEF_SAMPLERS(mc13, ev_m_chunk, 13)
DEF_SAMPLERS(mc16, ev_m_chunk, 16)
DEF_SAMPLERS(nt8, ev_nt, 8)
DEF_SAMPLERS(nt16, ev_nt, 16)
DEF_SAMPLERS(uh, ev_uh, 0)

typedef size_t (*sampler_fn)(const kp7kept*, const char*, size_t, uint64_t, size_t*, uint64_t*, size_t,
                             kp7_emit_pairs_fn, void*);
#define PICK(NAME) return branchy ? NAME##_br : NAME##_bat
#define BSW(PFX, ...)                                                                                                 \
  switch (t->B) { __VA_ARGS__ default: return NULL; }
static sampler_fn pick(const kp7kept* t, int branchy)
{
  switch (t->kind) {
  case KP7KEPT_E_CHUNK:
    BSW(ec, case 8: PICK(ec8); case 10: PICK(ec10); case 11: PICK(ec11); case 12: PICK(ec12); case 13: PICK(ec13);
        case 16: PICK(ec16);)
  case KP7KEPT_E_PLANES: BSW(ep, case 8: PICK(ep8); case 11: PICK(ep11); case 16: PICK(ep16);)
  case KP7KEPT_G_PLANES: BSW(gp, case 8: PICK(gp8); case 11: PICK(gp11); case 16: PICK(gp16);)
  case KP7KEPT_M_PLANES: BSW(mp, case 8: PICK(mp8); case 11: PICK(mp11); case 16: PICK(mp16);)
  case KP7KEPT_M_CHUNK:
    BSW(mc, case 8: PICK(mc8); case 10: PICK(mc10); case 11: PICK(mc11); case 12: PICK(mc12); case 13: PICK(mc13);
        case 16: PICK(mc16);)
  case KP7KEPT_NT_CHUNK: BSW(nt, case 8: PICK(nt8); case 16: PICK(nt16);)
  case KP7KEPT_UH: PICK(uh);
  default: return NULL;
  }
}

size_t kp7kept_sample(const kp7kept* t, int branchy, const char* seq, size_t len, uint64_t threshold,
                      size_t* pos_buf, uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx)
{
  if (!pcode_ready) pcode_init();
  sampler_fn f = pick(t, branchy);
  return f ? f(t, seq, len, threshold, pos_buf, hash_buf, cap, emit, ctx) : 0;
}

uint64_t kp7kept_eval(const kp7kept* t, uint64_t fwd, uint64_t rc)
{
#define EVB(EV) switch (t->B) { case 8: return EV(t, 8, fwd, rc); case 10: return EV(t, 10, fwd, rc); \
    case 11: return EV(t, 11, fwd, rc); case 12: return EV(t, 12, fwd, rc); case 13: return EV(t, 13, fwd, rc); \
    case 16: return EV(t, 16, fwd, rc); default: return 0; }
  switch (t->kind) {
  case KP7KEPT_E_CHUNK: EVB(ev_e_chunk)
  case KP7KEPT_E_PLANES: EVB(ev_e_planes)
  case KP7KEPT_G_PLANES: EVB(ev_g_planes)
  case KP7KEPT_M_PLANES: EVB(ev_m_planes)
  case KP7KEPT_M_CHUNK: EVB(ev_m_chunk)
  case KP7KEPT_NT_CHUNK: EVB(ev_nt)
  case KP7KEPT_UH: return ev_uh(t, 0, fwd, rc);
  default: return 0;
  }
#undef EVB
}

size_t kp7kept_uh_all(const kp7kept* t, const char* seq, size_t len, uint64_t* buf, size_t cap, kp7_emit_fn emit,
                      void* ctx)
{
  if (!pcode_ready) pcode_init();
  const unsigned k = t->k;
  size_t total = 0, n = 0, start = 0;
  FOR_KMERS({
    if (n && (n == cap || pos != start + n)) { emit(ctx, start, buf, n); n = 0; }
    if (n == 0) start = pos;
    buf[n++] = ev_uh(t, 0, fwd, rc);
    total++;
  })
  if (n) emit(ctx, start, buf, n);
  return total;
}
