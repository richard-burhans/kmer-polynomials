/* KP-7 over CM31 (p = 2^31 - 1): scalar reference, AVX2 8-lane and AVX-512 16-lane blocks. See kp7m31.h.
 *
 * SIMD arithmetic follows the standard Mersenne-31 vector techniques (as in Plonky3's mersenne-31 AVX2/AVX-512
 * packing and the Circle-STARK CM31 extension arithmetic): 32-bit lanes, 32x32->64 products by mul_epu32 on the
 * even lanes plus a 32-bit-shifted copy for the odd lanes, Mersenne folding (v & p) + (v >> 31), and canonical
 * correction by min_epu32(t, t - p). No code was copied; the techniques are public folklore for M31.
 */
#include "kp7m31.h"
#include "kp7seed.h"
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#define P KP7M31_P

/* ---- scalar field ---- */

static inline uint32_t madd(uint32_t a, uint32_t b) { uint32_t t = a + b; return t >= P ? t - P : t; } /* a, b < p */
static inline uint32_t msub(uint32_t a, uint32_t b) { return a >= b ? a - b : a + P - b; }
static inline uint32_t mred(uint64_t v) /* v < 2^63 -> [0, p) */
{
  uint64_t r = (v & P) + (v >> 31); /* < 2^31 + 2^32 */
  r = (r & P) + (r >> 31);          /* < 2^31 + 2 */
  return (uint32_t)(r >= P ? r - P : r);
}
static inline uint32_t mmul(uint32_t a, uint32_t b) { return mred((uint64_t)a * b); }
static inline m31c cadd(m31c a, m31c b) { m31c r = { madd(a.re, b.re), madd(a.im, b.im) }; return r; }
static inline m31c cneg(m31c a) { m31c r = { msub(0, a.re), msub(0, a.im) }; return r; }
static inline m31c cmul(m31c a, m31c b)
{
  m31c r = { msub(mmul(a.re, b.re), mmul(a.im, b.im)), madd(mmul(a.re, b.im), mmul(a.im, b.re)) };
  return r;
}
static inline uint64_t pack(m31c z) { return (uint64_t)z.re | ((uint64_t)z.im << 31); }

m31c kp7m31_mul(m31c a, m31c b) { return cmul(a, b); }
m31c kp7m31_add(m31c a, m31c b) { return cadd(a, b); }
uint64_t kp7m31_pack(m31c z) { return pack(z); }

static uint32_t mpow(uint32_t a, uint32_t e)
{
  uint32_t r = 1;
  while (e) { if (e & 1) r = mmul(r, a); a = mmul(a, a); e >>= 1; }
  return r;
}
static m31c cinv(m31c a) /* p = 3 mod 4: the norm re^2 + im^2 is nonzero for a != 0 */
{
  uint32_t n = madd(mmul(a.re, a.re), mmul(a.im, a.im));
  uint32_t ni = mpow(n, P - 2);
  m31c r = { mmul(a.re, ni), mmul(msub(0, a.im), ni) };
  return r;
}

/* ---- codes: (c >> 1) & 3 for ACGT/acgt; -1 otherwise ---- */

static int8_t CODE[256];
static int code_ready = 0;
static void code_init(void)
{
  for (int i = 0; i < 256; i++) CODE[i] = -1;
  const char* b = "ACGTacgt";
  for (int i = 0; i < 8; i++) CODE[(unsigned char)b[i]] = (int8_t)(((unsigned char)b[i] >> 1) & 3);
  code_ready = 1;
}
int kp7m31_code(char c) { if (!code_ready) code_init(); return CODE[(unsigned char)c]; }

/* ---- setup ---- */

void kp7m31_init(kp7m31_params* p, unsigned k, uint64_t seed_x, uint64_t seed_lam)
{
  kp7m31_init_seeds(p, k, &seed_x, &seed_lam);
}

void kp7m31_init_seeds(kp7m31_params* p, unsigned k, const uint64_t* seed_x, const uint64_t* seed_lam)
{
  if (!code_ready) code_init();
  p->k = k;
  /* x uniform on F_{p^2} \ {0}, lam uniform on F_{p^2}, from two independent sources (kp7seed.h); each part uniform
   * on [0, p) by masking to 31 bits and rejecting p */
  kp7_rng rx = kp7_rng_from(seed_x, KP7_TAG_X);
#ifdef KP7M31_MUTATE_SEED /* test-only: lam from x's stream, the pre-2026-11-11 derivation; the seed tests must fail */
  kp7_rng* rl = &rx;
  (void)seed_lam;
#else
  kp7_rng rlam = kp7_rng_from(seed_lam, KP7_TAG_LAM), *rl = &rlam;
#endif
  do {
    p->x.re = (uint32_t)kp7_uniform(&rx, P, 0);
    p->x.im = (uint32_t)kp7_uniform(&rx, P, 0);
  } while (p->x.re == 0 && p->x.im == 0);
  p->lam.re = (uint32_t)kp7_uniform(rl, P, 0);
  p->lam.im = (uint32_t)kp7_uniform(rl, P, 0);
  p->xinv = cinv(p->x);
  m31c one = { 1, 0 };
  p->xk1 = one;
  for (unsigned t = 0; t + 1 < k; t++) p->xk1 = cmul(p->xk1, p->x);
  p->xk = cmul(p->xk1, p->x);
  /* by code: A=0, C=1, T=2, G=3 */
  m31c A = { 1, 0 }, C = { 0, 1 }, T = { P - 1, 0 }, G = { 0, P - 1 };
  p->val[0] = A; p->val[1] = C; p->val[2] = T; p->val[3] = G;
  m31c dl_f = cmul(p->lam, cadd(one, cneg(p->xinv))), dl_r = cmul(p->lam, cadd(one, cneg(p->x)));
  for (int co = 0; co < 4; co++)
    for (int ci = 0; ci < 4; ci++) {
      p->tf[co * 4 + ci] = cadd(cadd(cneg(cmul(p->val[co], p->xinv)), cmul(p->val[ci], p->xk1)), dl_f);
      p->tr[co * 4 + ci] = cadd(cadd(cmul(p->val[co], p->xk), cneg(p->val[ci])), dl_r);
    }
}

int kp7m31_direct(const kp7m31_params* p, const char* kmer, m31c* F, m31c* R)
{
  if (!code_ready) code_init();
  m31c f = { 0, 0 }, r = { 0, 0 };
  for (unsigned t = 0; t < p->k; t++) { /* F: Horner from the end */
    int c = CODE[(unsigned char)kmer[p->k - 1 - t]];
    if (c < 0) return 0;
    f = cadd(cmul(f, p->x), p->val[c]);
  }
  for (unsigned t = 0; t < p->k; t++) r = cadd(cmul(r, p->x), p->val[CODE[(unsigned char)kmer[t]]]);
  *F = f;
  *R = cneg(r);
  return 1;
}

uint64_t kp7m31_sym(const kp7m31_params* p, m31c F, m31c R) { return pack(cmul(cadd(F, p->lam), cadd(R, p->lam))); }

/* ---- reference roller (canonical, no lam state) ---- */

void kp7m31_start(kp7m31_roller* r, const kp7m31_params* p, const char* seq, size_t len)
{
  if (!code_ready) code_init();
  r->p = p; r->seq = seq; r->len = len; r->pos = 0; r->valid = 0;
}
static int m_restart(kp7m31_roller* r, size_t from)
{
  size_t run = 0, i = from;
  while (run < r->p->k) {
    if (i >= r->len) { r->valid = 0; r->pos = r->len; return 0; }
    run = CODE[(unsigned char)r->seq[i]] < 0 ? 0 : run + 1;
    i++;
  }
  r->pos = i - r->p->k;
  kp7m31_direct(r->p, r->seq + r->pos, &r->F, &r->R);
  r->valid = 1;
  return 1;
}
int kp7m31_roll(kp7m31_roller* r)
{
  if (!r->valid) return m_restart(r, r->pos);
  const kp7m31_params* p = r->p;
  size_t in = r->pos + p->k;
  if (in >= r->len) { r->valid = 0; r->pos = r->len; return 0; }
  int ci = CODE[(unsigned char)r->seq[in]];
  if (ci < 0) return m_restart(r, in + 1);
  int co = CODE[(unsigned char)r->seq[r->pos]];
  /* F' = (F - v(co)) x^{-1} + v(ci) x^{k-1};  R' = R x + v(co) x^k - v(ci) */
  r->F = cadd(cmul(cadd(r->F, cneg(p->val[co])), p->xinv), cmul(p->val[ci], p->xk1));
  r->R = cadd(cadd(cmul(r->R, p->x), cmul(p->val[co], p->xk)), cneg(p->val[ci]));
  r->pos++;
  return 1;
}

/* ---- scalar block (lam state, merged tables) ---- */

static inline void lane_start(const kp7m31_params* p, const char* s, m31c* G, m31c* H)
{
  m31c F, R;
  kp7m31_direct(p, s, &F, &R);
  *G = cadd(F, p->lam);
  *H = cadd(R, p->lam);
}
/* continue one lane from state (G, H) at offset q0 (already emitted) for `steps` more k-mers */
static inline void scalar_tail(const kp7m31_params* p, const char* s, size_t q0, size_t steps, m31c G, m31c H,
                               uint64_t* out)
{
  const size_t k = p->k;
  for (size_t j = 0, q = q0; j < steps; j++, q++) {
    int idx = CODE[(unsigned char)s[q]] * 4 + CODE[(unsigned char)s[q + k]];
    G = cadd(cmul(G, p->xinv), p->tf[idx]);
    H = cadd(cmul(H, p->x), p->tr[idx]);
    out[q + 1] = pack(cmul(G, H));
  }
}
static void block_scalar(const kp7m31_params* p, const char* s, size_t nk, uint64_t* out)
{
  m31c G, H;
  lane_start(p, s, &G, &H);
  out[0] = pack(cmul(G, H));
  scalar_tail(p, s, 0, nk - 1, G, H, out);
}

/* lane-interleaved index pre-pass: t[i*L + l] = code(s[l*seg + i])*4 + code(s[l*seg + i + k]), i < seg - 1 */
static _Thread_local uint8_t TIDX[1 << 16];

/* ---- AVX2, 8 lanes ---- */
#if defined(__x86_64__)
#define AVX2 __attribute__((target("avx2")))

/* 64-bit-lane value v < 2^63 -> canonical residue in the low dword, zero high dword */
static inline AVX2 __m256i red64_2(__m256i v, __m256i p64, __m256i p32)
{
  __m256i r = _mm256_add_epi64(_mm256_and_si256(v, p64), _mm256_srli_epi64(v, 31));
  r = _mm256_add_epi64(_mm256_and_si256(r, p64), _mm256_srli_epi64(r, 31)); /* < 2^31 + 2, high dword 0 */
  return _mm256_min_epu32(r, _mm256_sub_epi32(r, p32));                    /* high dword: min(0, -p) = 0 */
}
/* (ar + i ai)(br + i bi) on 8 lanes of canonical 31-bit values; canonical result */
static inline AVX2 void cmul8(__m256i ar, __m256i ai, __m256i br, __m256i bi, __m256i* rr, __m256i* ri)
{
  const __m256i p64 = _mm256_set1_epi64x(P), p32 = _mm256_set1_epi32((int)P);
  const __m256i K = _mm256_set1_epi64x((long long)((uint64_t)P << 31)); /* p*2^31 > every product (< 2^62) */
  /* even lanes: mul_epu32 reads the low dword of each 64-bit lane */
  __m256i re_e = red64_2(_mm256_sub_epi64(_mm256_add_epi64(_mm256_mul_epu32(ar, br), K), _mm256_mul_epu32(ai, bi)), p64, p32);
  __m256i im_e = red64_2(_mm256_add_epi64(_mm256_mul_epu32(ar, bi), _mm256_mul_epu32(ai, br)), p64, p32);
  /* odd lanes: shift the high dwords down */
#ifdef KP7M31_MUTATE_ODD /* test-only: odd lanes reuse ar's EVEN dwords; `make mutation` expects FAIL */
  __m256i aro = ar;
#else
  __m256i aro = _mm256_srli_epi64(ar, 32);
#endif
  __m256i aio = _mm256_srli_epi64(ai, 32), bro = _mm256_srli_epi64(br, 32), bio = _mm256_srli_epi64(bi, 32);
  __m256i re_o = red64_2(_mm256_sub_epi64(_mm256_add_epi64(_mm256_mul_epu32(aro, bro), K), _mm256_mul_epu32(aio, bio)), p64, p32);
  __m256i im_o = red64_2(_mm256_add_epi64(_mm256_mul_epu32(aro, bio), _mm256_mul_epu32(aio, bro)), p64, p32);
  *rr = _mm256_or_si256(re_e, _mm256_slli_epi64(re_o, 32));
  *ri = _mm256_or_si256(im_e, _mm256_slli_epi64(im_o, 32));
}
static inline AVX2 __m256i addm8(__m256i a, __m256i b) /* canonical + canonical, < 2^32 */
{
  __m256i t = _mm256_add_epi32(a, b);
  return _mm256_min_epu32(t, _mm256_sub_epi32(t, _mm256_set1_epi32((int)P)));
}
/* 16-entry lookup: two registers, vpermd each, blend on idx > 7 */
static inline AVX2 __m256i look16(__m256i lo, __m256i hi, __m256i idx, __m256i sel)
{
  return _mm256_blendv_epi8(_mm256_permutevar8x32_epi32(lo, idx), _mm256_permutevar8x32_epi32(hi, idx), sel);
}

static AVX2 void block_avx2(const kp7m31_params* p, const char* s, size_t nk, uint64_t* out)
{
  enum { L = 8 };
  const size_t k = p->k, seg = nk / L, rem = nk - seg * L;
  uint32_t gr[L], gi[L], hr[L], hi[L];
  for (int l = 0; l < L; l++) {
    m31c G, H;
    lane_start(p, s + l * seg, &G, &H);
    gr[l] = G.re; gi[l] = G.im; hr[l] = H.re; hi[l] = H.im;
  }
  for (int l = 0; l < L; l++)
    for (size_t i = 0; i + 1 < seg; i++) {
      size_t q = l * seg + i;
      TIDX[i * L + l] = (uint8_t)(CODE[(unsigned char)s[q]] * 4 + CODE[(unsigned char)s[q + k]]);
    }
  uint32_t t[4][16];
  for (int e = 0; e < 16; e++) { t[0][e] = p->tf[e].re; t[1][e] = p->tf[e].im; t[2][e] = p->tr[e].re; t[3][e] = p->tr[e].im; }
  __m256i tlo[4], thi[4];
  for (int c = 0; c < 4; c++) { tlo[c] = _mm256_loadu_si256((const __m256i*)t[c]); thi[c] = _mm256_loadu_si256((const __m256i*)(t[c] + 8)); }
  const __m256i xr = _mm256_set1_epi32((int)p->x.re), xim = _mm256_set1_epi32((int)p->x.im);
  const __m256i yr = _mm256_set1_epi32((int)p->xinv.re), yim = _mm256_set1_epi32((int)p->xinv.im);
#ifdef KP7M31_MUTATE_BLEND /* test-only: wrong blend threshold, entry 8 comes from the low table; expects FAIL */
  const __m256i seven = _mm256_set1_epi32(8);
#else
  const __m256i seven = _mm256_set1_epi32(7);
#endif
  __m256i Gr = _mm256_loadu_si256((const __m256i*)gr), Gi = _mm256_loadu_si256((const __m256i*)gi);
  __m256i Hr = _mm256_loadu_si256((const __m256i*)hr), Hi = _mm256_loadu_si256((const __m256i*)hi);
  uint32_t ore[L], oim[L];
  for (size_t i = 0; i < seg; i++) {
    __m256i Sr, Si;
    cmul8(Gr, Gi, Hr, Hi, &Sr, &Si);
    _mm256_storeu_si256((__m256i*)ore, Sr);
    _mm256_storeu_si256((__m256i*)oim, Si);
    for (int l = 0; l < L; l++) out[l * seg + i] = (uint64_t)ore[l] | ((uint64_t)oim[l] << 31);
    if (i + 1 < seg) {
      __m256i idx = _mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i*)(TIDX + i * L)));
      __m256i sel = _mm256_cmpgt_epi32(idx, seven);
      __m256i Ar, Ai, Br, Bi;
      cmul8(Gr, Gi, yr, yim, &Ar, &Ai);
      cmul8(Hr, Hi, xr, xim, &Br, &Bi);
      Gr = addm8(Ar, look16(tlo[0], thi[0], idx, sel));
      Gi = addm8(Ai, look16(tlo[1], thi[1], idx, sel));
      Hr = addm8(Br, look16(tlo[2], thi[2], idx, sel));
      Hi = addm8(Bi, look16(tlo[3], thi[3], idx, sel));
    }
  }
  if (rem) {
    _mm256_storeu_si256((__m256i*)gr, Gr); _mm256_storeu_si256((__m256i*)gi, Gi);
    _mm256_storeu_si256((__m256i*)hr, Hr); _mm256_storeu_si256((__m256i*)hi, Hi);
    m31c G = { gr[L - 1], gi[L - 1] }, H = { hr[L - 1], hi[L - 1] };
    scalar_tail(p, s, L * seg - 1, rem, G, H, out);
  }
}
#endif /* __x86_64__ */

/* ---- AVX-512, 16 lanes (built only when the compiler targets AVX-512, e.g. -march=native on laila) ---- */
#if defined(__AVX512F__)
static inline __m512i red64_5(__m512i v, __m512i p64, __m512i p32)
{
  __m512i r = _mm512_add_epi64(_mm512_and_si512(v, p64), _mm512_srli_epi64(v, 31));
  r = _mm512_add_epi64(_mm512_and_si512(r, p64), _mm512_srli_epi64(r, 31));
  return _mm512_min_epu32(r, _mm512_sub_epi32(r, p32));
}
static inline void cmul16(__m512i ar, __m512i ai, __m512i br, __m512i bi, __m512i* rr, __m512i* ri)
{
  const __m512i p64 = _mm512_set1_epi64(P), p32 = _mm512_set1_epi32((int)P);
  const __m512i K = _mm512_set1_epi64((long long)((uint64_t)P << 31));
  __m512i re_e = red64_5(_mm512_sub_epi64(_mm512_add_epi64(_mm512_mul_epu32(ar, br), K), _mm512_mul_epu32(ai, bi)), p64, p32);
  __m512i im_e = red64_5(_mm512_add_epi64(_mm512_mul_epu32(ar, bi), _mm512_mul_epu32(ai, br)), p64, p32);
#ifdef KP7M31_MUTATE_ODD
  __m512i aro = ar;
#else
  __m512i aro = _mm512_srli_epi64(ar, 32);
#endif
  __m512i aio = _mm512_srli_epi64(ai, 32), bro = _mm512_srli_epi64(br, 32), bio = _mm512_srli_epi64(bi, 32);
  __m512i re_o = red64_5(_mm512_sub_epi64(_mm512_add_epi64(_mm512_mul_epu32(aro, bro), K), _mm512_mul_epu32(aio, bio)), p64, p32);
  __m512i im_o = red64_5(_mm512_add_epi64(_mm512_mul_epu32(aro, bio), _mm512_mul_epu32(aio, bro)), p64, p32);
  *rr = _mm512_or_si512(re_e, _mm512_slli_epi64(re_o, 32));
  *ri = _mm512_or_si512(im_e, _mm512_slli_epi64(im_o, 32));
}
static inline __m512i addm16(__m512i a, __m512i b)
{
  __m512i t = _mm512_add_epi32(a, b);
  return _mm512_min_epu32(t, _mm512_sub_epi32(t, _mm512_set1_epi32((int)P)));
}
static void block_avx512(const kp7m31_params* p, const char* s, size_t nk, uint64_t* out)
{
  enum { L = 16 };
  const size_t k = p->k, seg = nk / L, rem = nk - seg * L;
  uint32_t gr[L], gi[L], hr[L], hi[L];
  for (int l = 0; l < L; l++) {
    m31c G, H;
    lane_start(p, s + l * seg, &G, &H);
    gr[l] = G.re; gi[l] = G.im; hr[l] = H.re; hi[l] = H.im;
  }
  for (int l = 0; l < L; l++)
    for (size_t i = 0; i + 1 < seg; i++) {
      size_t q = l * seg + i;
      TIDX[i * L + l] = (uint8_t)(CODE[(unsigned char)s[q]] * 4 + CODE[(unsigned char)s[q + k]]);
    }
  uint32_t t[4][16];
  for (int e = 0; e < 16; e++) { t[0][e] = p->tf[e].re; t[1][e] = p->tf[e].im; t[2][e] = p->tr[e].re; t[3][e] = p->tr[e].im; }
  __m512i tab[4];
  for (int c = 0; c < 4; c++) tab[c] = _mm512_loadu_si512(t[c]); /* 16 entries: one register, no blend */
  const __m512i xr = _mm512_set1_epi32((int)p->x.re), xim = _mm512_set1_epi32((int)p->x.im);
  const __m512i yr = _mm512_set1_epi32((int)p->xinv.re), yim = _mm512_set1_epi32((int)p->xinv.im);
  __m512i Gr = _mm512_loadu_si512(gr), Gi = _mm512_loadu_si512(gi), Hr = _mm512_loadu_si512(hr), Hi = _mm512_loadu_si512(hi);
  uint32_t ore[L], oim[L];
  for (size_t i = 0; i < seg; i++) {
    __m512i Sr, Si;
    cmul16(Gr, Gi, Hr, Hi, &Sr, &Si);
    _mm512_storeu_si512(ore, Sr);
    _mm512_storeu_si512(oim, Si);
    for (int l = 0; l < L; l++) out[l * seg + i] = (uint64_t)ore[l] | ((uint64_t)oim[l] << 31);
    if (i + 1 < seg) {
      __m512i idx = _mm512_cvtepu8_epi32(_mm_loadu_si128((const __m128i*)(TIDX + i * L)));
      __m512i Ar, Ai, Br, Bi;
      cmul16(Gr, Gi, yr, yim, &Ar, &Ai);
      cmul16(Hr, Hi, xr, xim, &Br, &Bi);
      Gr = addm16(Ar, _mm512_permutexvar_epi32(idx, tab[0]));
      Gi = addm16(Ai, _mm512_permutexvar_epi32(idx, tab[1]));
      Hr = addm16(Br, _mm512_permutexvar_epi32(idx, tab[2]));
      Hi = addm16(Bi, _mm512_permutexvar_epi32(idx, tab[3]));
    }
  }
  if (rem) {
    _mm512_storeu_si512(gr, Gr); _mm512_storeu_si512(gi, Gi); _mm512_storeu_si512(hr, Hr); _mm512_storeu_si512(hi, Hi);
    m31c G = { gr[L - 1], gi[L - 1] }, H = { hr[L - 1], hi[L - 1] };
    scalar_tail(p, s, L * seg - 1, rem, G, H, out);
  }
}
#endif

int kp7m31_have(int impl)
{
  if (impl == KP7M31_SCALAR) return 1;
#if defined(__x86_64__)
  if (impl == KP7M31_AVX2) return __builtin_cpu_supports("avx2");
#endif
#if defined(__AVX512F__)
  if (impl == KP7M31_AVX512) return __builtin_cpu_supports("avx512f");
#endif
  return 0;
}

size_t kp7m31_sym_seq(const kp7m31_params* p, const char* seq, size_t len, int impl, uint64_t* buf, size_t cap,
                      kp7m31_emit_fn emit, void* ctx)
{
  if (!code_ready) code_init();
  if (!kp7m31_have(impl)) return 0;
  if (cap > (1u << 16)) cap = 1u << 16; /* TIDX holds one block */
  const size_t k = p->k;
  const size_t lanes = impl == KP7M31_AVX512 ? 16 : impl == KP7M31_AVX2 ? 8 : 1;
  size_t total = 0, i = 0;
  while (i < len) {
    while (i < len && CODE[(unsigned char)seq[i]] < 0) i++;
    size_t a = i;
    while (i < len && CODE[(unsigned char)seq[i]] >= 0) i++;
    size_t run = i - a;
    if (run < k) continue;
    size_t nk = run - k + 1;
    for (size_t off = 0; off < nk; off += cap) {
      size_t n = nk - off < cap ? nk - off : cap;
      const char* s = seq + a + off;
      if (n < lanes || lanes == 1) block_scalar(p, s, n, buf);
#if defined(__x86_64__)
      else if (impl == KP7M31_AVX2) block_avx2(p, s, n, buf);
#endif
#if defined(__AVX512F__)
      else if (impl == KP7M31_AVX512) block_avx512(p, s, n, buf);
#endif
      else block_scalar(p, s, n, buf);
      emit(ctx, a + off, buf, n);
      total += n;
    }
  }
  return total;
}
