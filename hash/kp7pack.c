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
