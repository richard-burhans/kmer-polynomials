/* Hybrid v2: packed-code gate, provable KP-7 hash for kept k-mers; see kp7pack.h. */
#include "kp7pack.h"

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
  /* one multiply, as kp7's mix64; the packed codes are dense small integers, so both xorshifts matter */
  x ^= x >> 32;
  x *= 0xd6e8feb86659fd93ULL;
  return x ^ (x >> 32);
}
#ifdef KP7PACK_MUTATE_GATE /* test-only: asymmetric gate; `make mutation_pack` expects FAIL */
static inline uint64_t gate(uint64_t f, uint64_t r) { return mix(f) + r; }
#else
static inline uint64_t gate(uint64_t f, uint64_t r) { return mix(f) + mix(r); }
#endif
#ifdef KP7PACK_MUTATE_RC /* test-only: complemented base lands one base too low */
#define RC_SHIFT(k) (2 * ((k) - 2))
#else
#define RC_SHIFT(k) (2 * ((k) - 1))
#endif

uint64_t kp7pack_gate(uint64_t fwd, uint64_t rc) { return gate(fwd, rc); }
uint64_t kp7pack_gate_seeded(uint64_t fwd, uint64_t rc, uint64_t seed) { return gate(fwd ^ seed, rc ^ seed); }

int kp7pack_codes(const char* kmer, unsigned k, uint64_t* fwd, uint64_t* rc)
{
  if (!pcode_ready) pcode_init();
  uint64_t f = 0, r = 0;
  for (unsigned t = 0; t < k; t++) {
    int c = PCODE[(unsigned char)kmer[t]];
    if (c < 0) return 0;
    f = (f << 2) | (uint64_t)c;                       /* base t at bits 2(k-1-t) */
    r |= (uint64_t)(3 - c) << (2 * t);                /* rc: complement of base t at bits 2t */
  }
  *fwd = f;
  *rc = r;
  return 1;
}

/* walk maximal N-free runs; for each run roll both strands and call BODY(pos, fwd, rc) for every k-mer */
#define PACK_FOR_KMERS(...)                                                                                           \
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
      if (have == k) { size_t pos = i + 1 - k; __VA_ARGS__ }                                                              \
    }                                                                                                                \
  }

size_t kp7pack_gate_seq(unsigned k, const char* seq, size_t len, uint64_t* buf, size_t cap, kp7_emit_fn emit,
                        void* ctx)
{
  if (!pcode_ready) pcode_init();
  size_t total = 0, n = 0, start = 0;
  PACK_FOR_KMERS({
    if (n && (n == cap || pos != start + n)) { emit(ctx, start, buf, n); n = 0; }
    if (n == 0) start = pos;
    buf[n++] = gate(fwd, rc);
    total++;
  })
  if (n) emit(ctx, start, buf, n);
  return total;
}

size_t kp7pack_sample(const kp7_params* p, const char* seq, size_t len, uint64_t threshold, size_t* pos_buf,
                      uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx)
{
  if (!pcode_ready) pcode_init();
  const unsigned k = p->k;
  size_t kept = 0, n = 0;
  PACK_FOR_KMERS({
    if (gate(fwd, rc) < threshold) {
      kp7_fe F, R;
      kp7_direct(p, seq + pos, &F, &R); /* the window is N-free */
      if (n == cap) { emit(ctx, pos_buf, hash_buf, n); n = 0; }
      pos_buf[n] = pos;
      hash_buf[n] = kp7_fold(kp7_sym(p, F, R));
      n++;
      kept++;
    }
  })
  if (n) emit(ctx, pos_buf, hash_buf, n);
  return kept;
}

/* ---- multiply-free direct evaluation (see kp7pack.h) ---- */

#define PQ KP7_Q
static inline uint64_t pq_red64(uint64_t s) /* s < 2^64 (sum of <= 8 canonical values) -> [0, q) */
{
  uint64_t r = (s & PQ) + (s >> 61);
  return r >= PQ ? r - PQ : r;
}
static inline uint64_t pq_red128(unsigned __int128 s) /* any sum of <= 64 canonical values -> [0, q) */
{
  unsigned __int128 t = (s & PQ) + (s >> 61);
  uint64_t r = (uint64_t)(t & PQ) + (uint64_t)(t >> 61);
  return r >= PQ ? r - PQ : r;
}

void kp7pack_tab_init(kp7pack_tab* tab, const kp7_params* p)
{
  const unsigned k = p->k;
  tab->k = k;
  tab->nchunk = (k + 3) / 4;
  kp7_fe xt = { 1, 0 };
  for (unsigned t = 0; t < 32; t++) {
    for (int d = 0; d < 4; d++) {
      kp7_fe v = kp7_mul(kp7_base(d), xt); /* val(d) x^t */
      kp7_fe nv = { v.re ? PQ - v.re : 0, v.im ? PQ - v.im : 0 };
      tab->N[t][d] = nv;                   /* -val(d) x^t */
    }
    xt = kp7_mul(xt, p->x);
  }
  for (unsigned j = 0; j < 8; j++)
    for (unsigned b = 0; b < 256; b++) {
      kp7_fe sum = { 0, 0 };
#ifdef KP7PACK_MUTATE_CHUNK /* test-only: the last chunk also sums positions >= k; `make mutation_pack` expects FAIL */
      for (unsigned u = 0; u < 4 && 4 * j + u < 32; u++)
#else
      for (unsigned u = 0; u < 4 && 4 * j + u < k; u++)
#endif
        sum = kp7_add(sum, tab->N[4 * j + u][(b >> (2 * u)) & 3]);
      tab->C[j][b] = sum;
    }
}

static inline void eval_chunk(const kp7pack_tab* tab, uint64_t fwd, uint64_t rc, kp7_fe* F, kp7_fe* R)
{
  uint64_t fr = 0, fi = 0, rr = 0, ri = 0;
  for (unsigned j = 0; j < tab->nchunk; j++) { /* nchunk <= 8: the sums stay below 8q < 2^64 */
    const kp7_fe a = tab->C[j][(rc >> (8 * j)) & 0xFF], b = tab->C[j][(fwd >> (8 * j)) & 0xFF];
    fr += a.re; fi += a.im; rr += b.re; ri += b.im;
  }
  F->re = pq_red64(fr); F->im = pq_red64(fi);
  R->re = pq_red64(rr); R->im = pq_red64(ri);
}
static inline void eval_pos(const kp7pack_tab* tab, uint64_t fwd, uint64_t rc, kp7_fe* F, kp7_fe* R)
{
  unsigned __int128 fr = 0, fi = 0, rr = 0, ri = 0;
  for (unsigned t = 0; t < tab->k; t++) {
    const kp7_fe a = tab->N[t][(rc >> (2 * t)) & 3], b = tab->N[t][(fwd >> (2 * t)) & 3];
    fr += a.re; fi += a.im; rr += b.re; ri += b.im;
  }
  F->re = pq_red128(fr); F->im = pq_red128(fi);
  R->re = pq_red128(rr); R->im = pq_red128(ri);
}

void kp7pack_eval(const kp7pack_tab* tab, int mode, uint64_t fwd, uint64_t rc, kp7_fe* F, kp7_fe* R)
{
  if (mode == KP7PACK_EVAL_POS) eval_pos(tab, fwd, rc, F, R);
  else eval_chunk(tab, fwd, rc, F, R);
}

size_t kp7pack_sample_eval(const kp7_params* p, const kp7pack_tab* tab, int mode, const char* seq, size_t len,
                           uint64_t threshold, size_t* pos_buf, uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit,
                           void* ctx)
{
  if (mode == KP7PACK_EVAL_DIRECT || !tab) return kp7pack_sample(p, seq, len, threshold, pos_buf, hash_buf, cap, emit, ctx);
  if (!pcode_ready) pcode_init();
  const unsigned k = p->k;
  size_t kept = 0, n = 0;
  PACK_FOR_KMERS({
    if (gate(fwd, rc) < threshold) {
      kp7_fe F, R;
      if (mode == KP7PACK_EVAL_POS) eval_pos(tab, fwd, rc, &F, &R); else eval_chunk(tab, fwd, rc, &F, &R);
      if (n == cap) { emit(ctx, pos_buf, hash_buf, n); n = 0; }
      pos_buf[n] = pos;
      hash_buf[n] = kp7_fold(kp7_sym(p, F, R));
      n++;
      kept++;
    }
  })
  if (n) emit(ctx, pos_buf, hash_buf, n);
  return kept;
}

size_t kp7pack_sym_all(const kp7_params* p, const kp7pack_tab* tab, const char* seq, size_t len, uint64_t* buf,
                       size_t cap, kp7_emit_fn emit, void* ctx)
{
  if (!pcode_ready) pcode_init();
  const unsigned k = p->k;
  size_t total = 0, n = 0, start = 0;
  PACK_FOR_KMERS({
    if (n && (n == cap || pos != start + n)) { emit(ctx, start, buf, n); n = 0; }
    if (n == 0) start = pos;
    kp7_fe F, R;
    eval_chunk(tab, fwd, rc, &F, &R);
    buf[n++] = kp7_fold(kp7_sym(p, F, R));
    total++;
  })
  if (n) emit(ctx, start, buf, n);
  return total;
}

/* ---- cheaper evaluation of kept k-mers (see kp7pack.h) ---- */

static inline uint64_t fq_red_wide(unsigned __int128 v) /* any 128-bit v -> [0, q) */
{
  unsigned __int128 t = (v & PQ) + (v >> 61);
  uint64_t r = (uint64_t)(t & PQ) + (uint64_t)(t >> 61);
  r = (r & PQ) + (r >> 61);
  return r >= PQ ? r - PQ : r;
}
/* (a + bi)(c + di) with parts < 4q, canonical result, folded as kp7_fold */
static inline uint64_t fq_sym_fold(uint64_t a, uint64_t b, uint64_t c, uint64_t d)
{
  const unsigned __int128 Q16 = (unsigned __int128)16 * PQ * PQ;
  unsigned __int128 ac = (unsigned __int128)a * c, bd = (unsigned __int128)b * d;
  unsigned __int128 ad = (unsigned __int128)a * d, bc = (unsigned __int128)b * c;
  uint64_t re = fq_red_wide(ac + (Q16 - bd)), im = fq_red_wide(ad + bc);
  return re ^ ((im << 32) | (im >> 32));
}
static inline uint64_t fq_part(uint64_t s) { return (s & PQ) + (s >> 61); } /* s < 16q -> < 2q */

void kp7pack_ftab_init(kp7pack_ftab* t, const kp7_params* p)
{
  kp7pack_tab base;
  kp7pack_tab_init(&base, p);
  const unsigned k = p->k;
  t->k = k;
  for (unsigned j = 0; j < 8; j++)
    for (unsigned b = 0; b < 256; b++) t->T8[j][b] = base.C[j][b];
  for (unsigned j = 0; j < 16; j++)
    for (unsigned b = 0; b < 16; b++) {
      kp7_fe sum = { 0, 0 };
      for (unsigned u = 0; u < 2 && 2 * j + u < k; u++) sum = kp7_add(sum, base.N[2 * j + u][(b >> (2 * u)) & 3]);
      t->T4[j][b] = sum;
    }
#ifndef KP7PACK_MUTATE_LAM /* test-only: forget lam; `make mutation_pack` expects FAIL */
  for (unsigned b = 0; b < 256; b++) t->T8[0][b] = kp7_add(t->T8[0][b], p->lam);
  for (unsigned b = 0; b < 16; b++) t->T4[0][b] = kp7_add(t->T4[0][b], p->lam);
#endif
}

#define M31P KP7M31_P
static inline uint32_t m31_red64(uint64_t v) /* v < 2^63 -> [0, p) */
{
  uint64_t r = (v & M31P) + (v >> 31);
  r = (r & M31P) + (r >> 31);
  return (uint32_t)(r >= M31P ? r - M31P : r);
}
static inline uint64_t m31_sym_pack(uint32_t a, uint32_t b, uint32_t c, uint32_t d) /* canonical inputs */
{
  const uint64_t PP = (uint64_t)M31P * M31P;
  uint64_t re = m31_red64((uint64_t)a * c + (PP - (uint64_t)b * d));
  uint64_t im = m31_red64((uint64_t)a * d + (uint64_t)b * c);
  return re | (im << 31); /* kp7m31_pack */
}
void kp7pack_mtab_init(kp7pack_mtab* t, const kp7m31_params* p)
{
  /* kp7pack codes A0 C1 G2 T3; kp7m31's val[] is indexed A C T G */
  const m31c v[4] = { p->val[0], p->val[1], p->val[3], p->val[2] };
  m31c N[32][4];
  m31c xt = { 1, 0 };
  for (unsigned s = 0; s < 32; s++) {
    for (int d = 0; d < 4; d++) {
      m31c m = kp7m31_mul(v[d], xt);
      m31c nm = { m.re ? M31P - m.re : 0, m.im ? M31P - m.im : 0 };
      N[s][d] = nm; /* -val(d) x^s */
    }
    xt = kp7m31_mul(xt, p->x);
  }
  const unsigned k = p->k;
  t->k = k;
  for (unsigned j = 0; j < 8; j++)
    for (unsigned b = 0; b < 256; b++) {
      m31c sum = { 0, 0 };
      for (unsigned u = 0; u < 4 && 4 * j + u < k; u++) sum = kp7m31_add(sum, N[4 * j + u][(b >> (2 * u)) & 3]);
      if (j == 0) sum = kp7m31_add(sum, p->lam);
      t->T8[j][b] = sum.re | ((uint64_t)sum.im << 32);
    }
  for (unsigned j = 0; j < 16; j++)
    for (unsigned b = 0; b < 16; b++) {
      m31c sum = { 0, 0 };
      for (unsigned u = 0; u < 2 && 2 * j + u < k; u++) sum = kp7m31_add(sum, N[2 * j + u][(b >> (2 * u)) & 3]);
      if (j == 0) sum = kp7m31_add(sum, p->lam);
      t->T4[j][b] = sum.re | ((uint64_t)sum.im << 32);
    }
}

/* one k-mer's hash from its packed codes */
static inline uint64_t eval_f8(const kp7pack_ftab* t, uint64_t fwd, uint64_t rc)
{
  uint64_t gr = 0, gi = 0, hr = 0, hi = 0;
#pragma GCC unroll 8
  for (unsigned j = 0; j < 8; j++) {
    const kp7_fe a = t->T8[j][(rc >> (8 * j)) & 0xFF], b = t->T8[j][(fwd >> (8 * j)) & 0xFF];
    gr += a.re; gi += a.im; hr += b.re; hi += b.im;
  }
  return fq_sym_fold(fq_part(gr), fq_part(gi), fq_part(hr), fq_part(hi));
}
static inline uint64_t eval_f4(const kp7pack_ftab* t, uint64_t fwd, uint64_t rc)
{
  /* 16 canonical terms would reach 16q > 2^64, so sum two halves of 8 terms (each < 8q) separately */
  uint64_t gr = 0, gi = 0, hr = 0, hi = 0;
  unsigned __int128 Gr = 0, Gi = 0, Hr = 0, Hi = 0;
#pragma GCC unroll 16
  for (unsigned j = 0; j < 16; j++) {
    const kp7_fe a = t->T4[j][(rc >> (4 * j)) & 0xF], b = t->T4[j][(fwd >> (4 * j)) & 0xF];
    if (j == 8) { Gr = gr; Gi = gi; Hr = hr; Hi = hi; gr = gi = hr = hi = 0; }
    gr += a.re; gi += a.im; hr += b.re; hi += b.im;
  }
  /* two halves of 8 terms each: each half < 8q; their partial reductions sum to < 4q */
  return fq_sym_fold(fq_part((uint64_t)Gr) + fq_part(gr), fq_part((uint64_t)Gi) + fq_part(gi),
                     fq_part((uint64_t)Hr) + fq_part(hr), fq_part((uint64_t)Hi) + fq_part(hi));
}
static inline uint64_t eval_m8(const kp7pack_mtab* t, uint64_t fwd, uint64_t rc)
{
  uint64_t gr = 0, gi = 0, hr = 0, hi = 0; /* 8 terms < 2^31 each */
#pragma GCC unroll 8
  for (unsigned j = 0; j < 8; j++) {
    const uint64_t a = t->T8[j][(rc >> (8 * j)) & 0xFF], b = t->T8[j][(fwd >> (8 * j)) & 0xFF];
    gr += a & 0xFFFFFFFFu; gi += a >> 32; hr += b & 0xFFFFFFFFu; hi += b >> 32;
  }
  return m31_sym_pack(m31_red64(gr), m31_red64(gi), m31_red64(hr), m31_red64(hi));
}
static inline uint64_t eval_m4(const kp7pack_mtab* t, uint64_t fwd, uint64_t rc)
{
  uint64_t gr = 0, gi = 0, hr = 0, hi = 0; /* 16 terms < 2^31 each */
#pragma GCC unroll 16
  for (unsigned j = 0; j < 16; j++) {
    const uint64_t a = t->T4[j][(rc >> (4 * j)) & 0xF], b = t->T4[j][(fwd >> (4 * j)) & 0xF];
    gr += a & 0xFFFFFFFFu; gi += a >> 32; hr += b & 0xFFFFFFFFu; hi += b >> 32;
  }
  return m31_sym_pack(m31_red64(gr), m31_red64(gi), m31_red64(hr), m31_red64(hi));
}

#define KP7PACK_DEF_ALL(NAME, TAB, EVAL)                                                                              \
  static size_t NAME(const TAB* t, const char* seq, size_t len, uint64_t* buf, size_t cap, kp7_emit_fn emit,          \
                     void* ctx)                                                                                       \
  {                                                                                                                  \
    if (!pcode_ready) pcode_init();                                                                                  \
    const unsigned k = t->k;                                                                                         \
    size_t total = 0, n = 0, start = 0;                                                                              \
    PACK_FOR_KMERS({                                                                                                 \
      if (n && (n == cap || pos != start + n)) { emit(ctx, start, buf, n); n = 0; }                                  \
      if (n == 0) start = pos;                                                                                       \
      buf[n++] = EVAL(t, fwd, rc);                                                                                   \
      total++;                                                                                                       \
    })                                                                                                               \
    if (n) emit(ctx, start, buf, n);                                                                                 \
    return total;                                                                                                    \
  }
KP7PACK_DEF_ALL(all_f8, kp7pack_ftab, eval_f8)
KP7PACK_DEF_ALL(all_f4, kp7pack_ftab, eval_f4)
KP7PACK_DEF_ALL(all_m8, kp7pack_mtab, eval_m8)
KP7PACK_DEF_ALL(all_m4, kp7pack_mtab, eval_m4)

size_t kp7pack_sym_all_f(const kp7pack_ftab* t, int W, const char* seq, size_t len, uint64_t* buf, size_t cap,
                         kp7_emit_fn emit, void* ctx)
{ return W == 4 ? all_f4(t, seq, len, buf, cap, emit, ctx) : all_f8(t, seq, len, buf, cap, emit, ctx); }
size_t kp7pack_sym_all_m(const kp7pack_mtab* t, int W, const char* seq, size_t len, uint64_t* buf, size_t cap,
                         kp7_emit_fn emit, void* ctx)
{ return W == 4 ? all_m4(t, seq, len, buf, cap, emit, ctx) : all_m8(t, seq, len, buf, cap, emit, ctx); }

/* Sampling: kept k-mers are collected WITHOUT a branch (the slot is always written, the count advances by the
 * comparison), then evaluated in a tight loop over the batch; output order is position order, as kp7pack_sample. */
#define KP7PACK_BATCH 512
#define KP7PACK_DEF_SAMPLE(NAME, TAB, EVAL)                                                                           \
  static size_t NAME(const TAB* t, const char* seq, size_t len, uint64_t thr, size_t* pos_buf, uint64_t* hash_buf,   \
                     size_t cap, kp7_emit_pairs_fn emit, void* ctx)                                                   \
  {                                                                                                                  \
    if (!pcode_ready) pcode_init();                                                                                  \
    const unsigned k = t->k;                                                                                         \
    const size_t bmax = cap < KP7PACK_BATCH ? cap : KP7PACK_BATCH;                                                   \
    uint64_t bf[KP7PACK_BATCH], br[KP7PACK_BATCH];                                                                   \
    size_t bp[KP7PACK_BATCH], nb = 0, n = 0, kept = 0;                                                               \
    PACK_FOR_KMERS({                                                                                                 \
      bf[nb] = fwd; br[nb] = rc; bp[nb] = pos;                                                                       \
      nb += gate(fwd, rc) < thr;                                                                                     \
      if (nb == bmax) {                                                                                              \
        for (size_t q = 0; q < nb; q++) {                                                                            \
          if (n == cap) { emit(ctx, pos_buf, hash_buf, n); n = 0; }                                                  \
          pos_buf[n] = bp[q]; hash_buf[n] = EVAL(t, bf[q], br[q]); n++;                                              \
        }                                                                                                            \
        kept += nb; nb = 0;                                                                                          \
      }                                                                                                              \
    })                                                                                                               \
    for (size_t q = 0; q < nb; q++) {                                                                                \
      if (n == cap) { emit(ctx, pos_buf, hash_buf, n); n = 0; }                                                      \
      pos_buf[n] = bp[q]; hash_buf[n] = EVAL(t, bf[q], br[q]); n++;                                                  \
    }                                                                                                                \
    kept += nb;                                                                                                      \
    if (n) emit(ctx, pos_buf, hash_buf, n);                                                                          \
    return kept;                                                                                                     \
  }
KP7PACK_DEF_SAMPLE(samp_f8, kp7pack_ftab, eval_f8)
KP7PACK_DEF_SAMPLE(samp_f4, kp7pack_ftab, eval_f4)
KP7PACK_DEF_SAMPLE(samp_m8, kp7pack_mtab, eval_m8)
KP7PACK_DEF_SAMPLE(samp_m4, kp7pack_mtab, eval_m4)

size_t kp7pack_sample_f(const kp7pack_ftab* t, int W, const char* seq, size_t len, uint64_t threshold, size_t* pos_buf,
                        uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx)
{
  return W == 4 ? samp_f4(t, seq, len, threshold, pos_buf, hash_buf, cap, emit, ctx)
                : samp_f8(t, seq, len, threshold, pos_buf, hash_buf, cap, emit, ctx);
}
size_t kp7pack_sample_m(const kp7pack_mtab* t, int W, const char* seq, size_t len, uint64_t threshold, size_t* pos_buf,
                        uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx)
{
  return W == 4 ? samp_m4(t, seq, len, threshold, pos_buf, hash_buf, cap, emit, ctx)
                : samp_m8(t, seq, len, threshold, pos_buf, hash_buf, cap, emit, ctx);
}
