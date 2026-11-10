/* KP-7 over the Eisenstein encoding in F_q, q = 2^61 - 1; see kp7e.h. */
#include "kp7e.h"

#define Q KP7E_Q

/* ---- field ---- */
static inline uint64_t red(unsigned __int128 p) /* p < 2^124 */
{
  uint64_t r = ((uint64_t)p & Q) + (uint64_t)(p >> 61);
  r = (r & Q) + (r >> 61);
  return r >= Q ? r - Q : r;
}
static inline uint64_t addq(uint64_t a, uint64_t b) { uint64_t r = a + b; return r >= Q ? r - Q : r; }
static inline uint64_t subq(uint64_t a, uint64_t b) { return a >= b ? a - b : a + Q - b; }
static inline uint64_t mulq(uint64_t a, uint64_t b) { return red((unsigned __int128)a * b); } /* a, b < 2^63 */
static inline uint64_t part(uint64_t s) { return (s & Q) + (s >> 61); } /* s < 2^64 -> < q + 8 (canonical after canon) */

static uint64_t powq(uint64_t a, uint64_t e)
{
  uint64_t r = 1;
  while (e) { if (e & 1) r = mulq(r, a); a = mulq(a, a); e >>= 1; }
  return r;
}
uint64_t kp7e_mul(uint64_t a, uint64_t b) { return mulq(a, b); }
uint64_t kp7e_add(uint64_t a, uint64_t b) { return addq(a, b); }
uint64_t kp7e_inv(uint64_t a) { return powq(a, Q - 2); }

static int8_t CODE[256];
static int code_ready = 0;
static void code_init(void)
{
  for (int i = 0; i < 256; i++) CODE[i] = -1;
  CODE['A'] = CODE['a'] = 0; CODE['C'] = CODE['c'] = 1;
  CODE['G'] = CODE['g'] = 2; CODE['T'] = CODE['t'] = 3;
  code_ready = 1;
}
int kp7e_code(char c) { if (!code_ready) code_init(); return CODE[(unsigned char)c]; }

static uint64_t splitmix64(uint64_t* s)
{
  uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static uint64_t rand_q(uint64_t* s)
{
  for (;;) { uint64_t v = splitmix64(s) & Q; if (v > 1 && v != Q) return v; }
}

void kp7e_init(kp7e_params* p, unsigned k, uint64_t seed)
{
  if (!code_ready) code_init();
  /* w: the first g = 2, 3, ... with g^((q-1)/3) != 1 */
  uint64_t w = 1;
  for (uint64_t g = 2; w == 1; g++) w = powq(g, (Q - 1) / 3);
  p->k = k;
  p->w = w;
  p->val[0] = 1; p->val[1] = w; p->val[2] = subq(0, w); p->val[3] = Q - 1;
  uint64_t s = seed;
  p->x = rand_q(&s);
  p->lam = rand_q(&s);
  p->xinv = kp7e_inv(p->x);
  p->xk1 = powq(p->x, k - 1);
  p->xk = mulq(p->xk1, p->x);
  const uint64_t dlf = mulq(p->lam, subq(1, p->xinv)), dlr = mulq(p->lam, subq(1, p->x));
  for (int co = 0; co < 4; co++)
    for (int ci = 0; ci < 4; ci++) {
      /* TF = -s_out x^{-1} + s_in x^{k-1} + lam(1 - x^{-1});  TR = s_out x^k - s_in + lam(1 - x) */
      uint64_t tf = addq(subq(0, mulq(p->val[co], p->xinv)), mulq(p->val[ci], p->xk1));
#ifdef KP7E_MUTATE_ROLL /* test-only: the incoming base enters the reverse strand with the wrong sign */
      uint64_t tr = addq(mulq(p->val[co], p->xk), p->val[ci]);
#else
      uint64_t tr = subq(mulq(p->val[co], p->xk), p->val[ci]);
#endif
      p->tf[co * 4 + ci] = addq(tf, dlf);
      p->tr[co * 4 + ci] = addq(tr, dlr);
    }
}

int kp7e_direct(const kp7e_params* p, const char* kmer, uint64_t* F, uint64_t* R)
{
  /* F = sum s_t x^t (Horner from the end); R = -sum s_t x^{k-1-t} (Horner from the start) */
  uint64_t f = 0, r = 0;
  for (unsigned t = 0; t < p->k; t++) {
    int c = kp7e_code(kmer[p->k - 1 - t]);
    if (c < 0) return 0;
    f = addq(mulq(f, p->x), p->val[c]);
  }
  for (unsigned t = 0; t < p->k; t++) r = addq(mulq(r, p->x), p->val[kp7e_code(kmer[t])]);
  *F = f;
  *R = subq(0, r);
  return 1;
}

uint64_t kp7e_sym(const kp7e_params* p, uint64_t F, uint64_t R) { return mulq(addq(F, p->lam), addq(R, p->lam)); }

/* ---- scalar reference roller ---- */
void kp7e_start(kp7e_roller* r, const kp7e_params* p, const char* seq, size_t len)
{
  r->p = p; r->seq = seq; r->len = len; r->pos = 0; r->valid = 0;
}
static int restart(kp7e_roller* r, size_t from)
{
  const unsigned k = r->p->k;
  size_t run = 0, i = from;
  while (run < k) {
    if (i >= r->len) { r->valid = 0; return 0; }
    run = CODE[(unsigned char)r->seq[i]] < 0 ? 0 : run + 1;
    i++;
  }
  r->pos = i - k;
  kp7e_direct(r->p, r->seq + r->pos, &r->F, &r->R);
  r->valid = 1;
  return 1;
}
int kp7e_roll(kp7e_roller* r)
{
  if (!code_ready) code_init();
  if (!r->valid) return restart(r, r->pos);
  const kp7e_params* p = r->p;
  size_t in = r->pos + p->k;
  if (in >= r->len) { r->valid = 0; r->pos = r->len; return 0; }
  int ci = CODE[(unsigned char)r->seq[in]];
  if (ci < 0) return restart(r, in + 1);
  int co = CODE[(unsigned char)r->seq[r->pos]];
  /* F' = (F - s_out) x^{-1} + s_in x^{k-1};   R' = x (R + s_out x^{k-1}) - s_in */
  r->F = addq(mulq(subq(r->F, p->val[co]), p->xinv), mulq(p->val[ci], p->xk1));
  r->R = subq(mulq(addq(r->R, mulq(p->val[co], p->xk1)), p->x), p->val[ci]);
  r->pos++;
  return 1;
}

/* ---- optimized: lam state, merged tables, lazy state in [0, 2q), L lanes per block ---- */
#ifdef KP7E_MUTATE_LANE /* test-only: the last lane's remainder lands one slot early */
#define KP7E_REM_OFF -1
#else
#define KP7E_REM_OFF
#endif
static inline void lane_init(const kp7e_params* p, const char* s, uint64_t* G, uint64_t* H)
{
  uint64_t F, R;
  kp7e_direct(p, s, &F, &R);
  *G = addq(F, p->lam);
  *H = addq(R, p->lam);
}
/* G, H < 2q: mulq takes (2q) q < 2^123, the symmetric product (2q)^2 < 2^124; red returns canonical */
#define KP7E_DEF_BLOCK(L)                                                                                             \
  static void block_##L(const kp7e_params* p, const char* s, size_t nk, uint64_t* out)                               \
  {                                                                                                                  \
    const size_t k = p->k, seg = nk / L, rem = nk - seg * L;                                                         \
    const uint64_t xi = p->xinv, x = p->x;                                                                           \
    const uint64_t *tf = p->tf, *tr = p->tr;                                                                         \
    uint64_t g[L], h[L];                                                                                             \
    for (int l = 0; l < L; l++) lane_init(p, s + l * seg, &g[l], &h[l]);                                             \
    for (size_t i = 0; i < seg; i++) {                                                                               \
      for (int l = 0; l < L; l++) {                                                                                  \
        size_t q = l * seg + i;                                                                                      \
        out[q] = mulq(g[l], h[l]);                                                                                   \
        if (i + 1 < seg) {                                                                                           \
          int idx = CODE[(unsigned char)s[q]] * 4 + CODE[(unsigned char)s[q + k]];                                   \
          g[l] = mulq(g[l], xi) + tf[idx];                                                                           \
          h[l] = mulq(h[l], x) + tr[idx];                                                                            \
        }                                                                                                            \
      }                                                                                                              \
    }                                                                                                                \
    for (size_t j = 0, q = L * seg - 1; j < rem; j++, q++) {                                                         \
      int idx = CODE[(unsigned char)s[q]] * 4 + CODE[(unsigned char)s[q + k]];                                       \
      g[L - 1] = mulq(g[L - 1], xi) + tf[idx];                                                                       \
      h[L - 1] = mulq(h[L - 1], x) + tr[idx];                                                                        \
      out[q + 1 KP7E_REM_OFF] = mulq(g[L - 1], h[L - 1]);                                                            \
    }                                                                                                                \
  }
KP7E_DEF_BLOCK(1)
KP7E_DEF_BLOCK(2)

size_t kp7e_sym_seq(const kp7e_params* p, const char* seq, size_t len, int lanes, uint64_t* buf, size_t cap,
                    kp7e_emit_fn emit, void* ctx)
{
  if (!code_ready) code_init();
  const size_t k = p->k;
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
      if (lanes >= 2 && n >= 2) block_2(p, s, n, buf); else block_1(p, s, n, buf);
      emit(ctx, a + off, buf, n);
      total += n;
    }
  }
  return total;
}

/* ---- hybrid v2: packed-code gate and table evaluation ---- */
static inline uint64_t mix(uint64_t x) /* identical to kp7pack.c's mix */
{
  x ^= x >> 32;
  x *= 0xd6e8feb86659fd93ULL;
  return x ^ (x >> 32);
}
static inline uint64_t gate(uint64_t f, uint64_t r) { return mix(f) + mix(r); }
uint64_t kp7e_gate(uint64_t fwd, uint64_t rc) { return gate(fwd, rc); }

int kp7e_codes(const char* kmer, unsigned k, uint64_t* fwd, uint64_t* rc)
{
  if (!code_ready) code_init();
  uint64_t f = 0, r = 0;
  for (unsigned t = 0; t < k; t++) {
    int c = CODE[(unsigned char)kmer[t]];
    if (c < 0) return 0;
    f = (f << 2) | (uint64_t)c;          /* base t at bits 2(k-1-t) */
    r |= (uint64_t)(3 - c) << (2 * t);   /* rc: complement of base t at bits 2t */
  }
  *fwd = f;
  *rc = r;
  return 1;
}

void kp7e_tab_init(kp7e_tab* t, const kp7e_params* p)
{
  t->k = p->k;
  uint64_t N[32][4];
  uint64_t xt = 1;
  for (unsigned u = 0; u < 32; u++) {
    for (int d = 0; d < 4; d++) N[u][d] = u < p->k ? subq(0, mulq(p->val[d], xt)) : 0;
    xt = mulq(xt, p->x);
  }
  for (unsigned j = 0; j < 8; j++)
    for (unsigned b = 0; b < 256; b++) {
      uint64_t v = 0;
      for (unsigned u = 0; u < 4; u++) v = addq(v, N[4 * j + u][(b >> (2 * u)) & 3]);
#ifndef KP7E_MUTATE_TAB /* test-only: lam left out of the table */
      if (j == 0) v = addq(v, p->lam);
#endif
      t->T[j][b] = v;
    }
}

/* 8 canonical entries sum to < 8q < 2^64; part() leaves < q + 8 < 2q; the product of two such is < 2^124 */
static inline uint64_t tab_sym(const kp7e_tab* t, uint64_t fwd, uint64_t rc)
{
  uint64_t g = 0, h = 0;
  for (unsigned j = 0; j < 8; j++) {
    g += t->T[j][(rc >> (8 * j)) & 0xFF];
    h += t->T[j][(fwd >> (8 * j)) & 0xFF];
  }
  return mulq(part(g), part(h));
}
uint64_t kp7e_tab_sym(const kp7e_tab* t, uint64_t fwd, uint64_t rc) { return tab_sym(t, fwd, rc); }

#define KP7E_FOR_KMERS(...)                                                                                           \
  const uint64_t mask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);                                                   \
  const unsigned rsh = 2 * (k - 1);                                                                                  \
  size_t i = 0;                                                                                                      \
  while (i < len) {                                                                                                  \
    uint64_t fwd = 0, rc = 0;                                                                                        \
    unsigned have = 0;                                                                                               \
    for (; i < len; i++) {                                                                                           \
      int c = CODE[(unsigned char)seq[i]];                                                                           \
      if (c < 0) { have = 0; fwd = rc = 0; continue; }                                                               \
      fwd = ((fwd << 2) | (uint64_t)c) & mask;                                                                       \
      rc = (rc >> 2) | ((uint64_t)(3 - c) << rsh);                                                                   \
      if (have < k) have++;                                                                                          \
      if (have == k) { size_t pos = i + 1 - k; __VA_ARGS__ }                                                         \
    }                                                                                                                \
  }

size_t kp7e_sample_b(const kp7e_tab* t, const char* seq, size_t len, uint64_t thr, size_t* pos_buf, uint64_t* hash_buf,
                     size_t cap, kp7e_emit_pairs_fn emit, void* ctx)
{
  if (!code_ready) code_init();
  const unsigned k = t->k;
  size_t n = 0, kept = 0;
  KP7E_FOR_KMERS({
    if (gate(fwd, rc) < thr) {
      if (n == cap) { emit(ctx, pos_buf, hash_buf, n); n = 0; }
      pos_buf[n] = pos; hash_buf[n] = tab_sym(t, fwd, rc); n++; kept++;
    }
  })
  if (n) emit(ctx, pos_buf, hash_buf, n);
  return kept;
}

#define KP7E_BATCH 512
size_t kp7e_sample_batch(const kp7e_tab* t, const char* seq, size_t len, uint64_t thr, size_t* pos_buf,
                         uint64_t* hash_buf, size_t cap, kp7e_emit_pairs_fn emit, void* ctx)
{
  if (!code_ready) code_init();
  const unsigned k = t->k;
  const size_t bmax = cap < KP7E_BATCH ? cap : KP7E_BATCH;
  uint64_t bf[KP7E_BATCH], br[KP7E_BATCH];
  size_t bp[KP7E_BATCH], nb = 0, n = 0, kept = 0;
  KP7E_FOR_KMERS({
    bf[nb] = fwd; br[nb] = rc; bp[nb] = pos;
    nb += gate(fwd, rc) < thr;
    if (nb == bmax) {
      for (size_t q = 0; q < nb; q++) {
        if (n == cap) { emit(ctx, pos_buf, hash_buf, n); n = 0; }
        pos_buf[n] = bp[q]; hash_buf[n] = tab_sym(t, bf[q], br[q]); n++;
      }
      kept += nb; nb = 0;
    }
  })
  for (size_t q = 0; q < nb; q++) {
    if (n == cap) { emit(ctx, pos_buf, hash_buf, n); n = 0; }
    pos_buf[n] = bp[q]; hash_buf[n] = tab_sym(t, bf[q], br[q]); n++;
  }
  kept += nb;
  if (n) emit(ctx, pos_buf, hash_buf, n);
  return kept;
}
