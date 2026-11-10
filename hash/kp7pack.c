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
