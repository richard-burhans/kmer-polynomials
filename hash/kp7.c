/* KP-6 / KP-7 rolling hashes; see kp7.h. */
#include "kp7.h"
#include <string.h>

#define Q KP7_Q
static const unsigned __int128 Q2 = (unsigned __int128)KP7_Q * KP7_Q;

/* ---- field arithmetic (inline in this translation unit, exported wrappers below) ---- */

static inline uint64_t red(unsigned __int128 p) /* p < 2^124: p >> 61 < 2^63, so the first sum fits 64 bits */
{
  uint64_t r = ((uint64_t)p & Q) + (uint64_t)(p >> 61);
  r = (r & Q) + (r >> 61);
  return r >= Q ? r - Q : r;
}
/* any 128-bit p: the first fold stays in 128 bits ((p & Q) + (p >> 61) < 2^68), the second fits 64 bits */
static inline uint64_t red2(unsigned __int128 p)
{
  unsigned __int128 t = (p & Q) + (p >> 61); /* < 2^61 + 2^67 */
  uint64_t r = (uint64_t)(t & Q) + (uint64_t)(t >> 61);
  r = (r & Q) + (r >> 61);
  return r >= Q ? r - Q : r;
}
static inline uint64_t addq(uint64_t a, uint64_t b) { uint64_t r = a + b; return r >= Q ? r - Q : r; }
static inline uint64_t subq(uint64_t a, uint64_t b) { return a >= b ? a - b : a + Q - b; }
static inline uint64_t mulq(uint64_t a, uint64_t b) { return red((unsigned __int128)a * b); }

static inline kp7_fe fadd(kp7_fe a, kp7_fe b) { kp7_fe r = { addq(a.re, b.re), addq(a.im, b.im) }; return r; }
static inline kp7_fe fmul(kp7_fe a, kp7_fe b)
{
  unsigned __int128 ac = (unsigned __int128)a.re * b.re, bd = (unsigned __int128)a.im * b.im;
  unsigned __int128 ad = (unsigned __int128)a.re * b.im, bc = (unsigned __int128)a.im * b.re;
  kp7_fe r = { red(ac + (Q2 - bd)), red(ad + bc) };
  return r;
}
/* Karatsuba: 3 widening multiplies; see kp7.h */
static inline kp7_fe fmul_k(kp7_fe a, kp7_fe b)
{
  unsigned __int128 ac = (unsigned __int128)a.re * b.re, bd = (unsigned __int128)a.im * b.im;
  unsigned __int128 s = (unsigned __int128)(a.re + a.im) * (b.re + b.im); /* sums < 2^62, product < 2^124 */
#ifdef KP7_MUTATE_K /* test-only: drop one Karatsuba term; `make mutation` expects the tests to FAIL */
  kp7_fe r = { red(ac + (Q2 - bd)), red(s - ac) };
#else
  kp7_fe r = { red(ac + (Q2 - bd)), red(s - ac - bd) };                   /* s >= ac + bd: ad + bc >= 0 */
#endif
  return r;
}
static inline kp7_fe fmul_p(const kp7_params* p, kp7_fe a, kp7_fe b) { return p->karatsuba ? fmul_k(a, b) : fmul(a, b); }
static inline kp7_fe fmul_real(kp7_fe a, uint64_t x) { kp7_fe r = { mulq(a.re, x), mulq(a.im, x) }; return r; }
static inline uint64_t fold(kp7_fe z) { return z.re ^ ((z.im << 32) | (z.im >> 32)); }

static uint64_t powq(uint64_t a, uint64_t e)
{
  uint64_t r = 1;
  while (e) { if (e & 1) r = mulq(r, a); a = mulq(a, a); e >>= 1; }
  return r;
}

uint64_t kp7_addq(uint64_t a, uint64_t b) { return addq(a, b); }
uint64_t kp7_subq(uint64_t a, uint64_t b) { return subq(a, b); }
uint64_t kp7_mulq(uint64_t a, uint64_t b) { return mulq(a, b); }
kp7_fe kp7_add(kp7_fe a, kp7_fe b) { return fadd(a, b); }
kp7_fe kp7_mul(kp7_fe a, kp7_fe b) { return fmul(a, b); }
kp7_fe kp7_mul_karatsuba(kp7_fe a, kp7_fe b) { return fmul_k(a, b); }
static inline kp7_fe fmul_wide(kp7_fe a, kp7_fe b);
kp7_fe kp7_mul_wide(kp7_fe a, kp7_fe b) { return fmul_wide(a, b); } /* test hook: parts may be < 4q */
uint64_t kp7_fold(kp7_fe z) { return fold(z); }

kp7_fe kp7_inv(kp7_fe a)
{
  /* q = 3 mod 4, so -1 is a non-residue and the norm a.re^2 + a.im^2 is nonzero for a != 0 */
  uint64_t n = addq(mulq(a.re, a.re), mulq(a.im, a.im));
  uint64_t ni = powq(n, Q - 2);
  kp7_fe r = { mulq(a.re, ni), mulq(subq(0, a.im), ni) };
  return r;
}

kp7_fe kp7_base(int code)
{
  static const kp7_fe B[4] = { { 1, 0 }, { 0, 1 }, { 0, KP7_Q - 1 }, { KP7_Q - 1, 0 } }; /* A C G T */
  return B[code];
}

static int8_t CODE[256];
static int code_ready = 0;
static void code_init(void)
{
  for (int i = 0; i < 256; i++) CODE[i] = -1;
  CODE['A'] = CODE['a'] = 0; CODE['C'] = CODE['c'] = 1;
  CODE['G'] = CODE['g'] = 2; CODE['T'] = CODE['t'] = 3;
  code_ready = 1;
}
int kp7_code(char c) { if (!code_ready) code_init(); return CODE[(unsigned char)c]; }

/* ---- setup ---- */

static uint64_t splitmix64(uint64_t* s)
{
  uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static uint64_t rand_q(uint64_t* s) /* uniform in [1, q) */
{
  for (;;) { uint64_t v = splitmix64(s) & Q; if (v != 0 && v != Q) return v; }
}

static kp7_fe fneg(kp7_fe a) { kp7_fe r = { subq(0, a.re), subq(0, a.im) }; return r; }

void kp7_init(kp7_params* p, unsigned k, int variant, uint64_t seed)
{
  if (!code_ready) code_init();
  uint64_t s = seed;
  p->k = k;
  p->variant = variant;
  p->karatsuba = 0;
  p->opt = 0;
  p->x.re = rand_q(&s);
  p->x.im = variant == KP7_XREAL ? 0 : rand_q(&s);
  p->lam.re = rand_q(&s);
  p->lam.im = rand_q(&s);
  p->xinv = kp7_inv(p->x);
  kp7_fe xk1 = { 1, 0 };
  for (unsigned t = 0; t + 1 < k; t++) xk1 = fmul(xk1, p->x); /* x^{k-1} */
  kp7_fe xk = fmul(xk1, p->x);                                /* x^k */
  for (int c = 0; c < 4; c++) {
    kp7_fe b = kp7_base(c);
    p->out_f[c] = fneg(fmul(b, p->xinv));
    p->in_f[c] = fmul(b, xk1);
    p->out_r[c] = fmul(b, xk);
#ifdef KP7_MUTATE /* test-only: break one rolling term; `make mutation` expects the tests to FAIL */
    p->in_r[c] = b;
#else
    p->in_r[c] = fneg(b);
#endif
  }
  kp7_fe one = { 1, 0 };
#ifdef KP7_MUTATE_OPT /* test-only: wrong lam-state constant; `make mutation` expects the tests to FAIL */
  p->dl_f = fmul(p->lam, fadd(one, p->xinv));
#else
  p->dl_f = fmul(p->lam, fadd(one, fneg(p->xinv))); /* lam(1 - x^{-1}) */
#endif
  p->dl_r = fmul(p->lam, fadd(one, fneg(p->x)));    /* lam(1 - x) */
  for (int co = 0; co < 4; co++)
    for (int ci = 0; ci < 4; ci++) {
      p->tf[co * 4 + ci] = fadd(p->out_f[co], p->in_f[ci]);
      p->tr[co * 4 + ci] = fadd(p->out_r[co], p->in_r[ci]);
      p->tfl[co * 4 + ci] = fadd(p->tf[co * 4 + ci], p->dl_f);
      p->trl[co * 4 + ci] = fadd(p->tr[co * 4 + ci], p->dl_r);
    }
}

int kp7_direct(const kp7_params* p, const char* kmer, kp7_fe* F, kp7_fe* R)
{
  /* F = sum s_t x^t (Horner from the end); R = -sum s_t x^{k-1-t} (Horner from the start) */
  kp7_fe f = { 0, 0 }, r = { 0, 0 };
  for (unsigned t = 0; t < p->k; t++) {
    int c = kp7_code(kmer[p->k - 1 - t]);
    if (c < 0) return 0;
    f = fadd(fmul(f, p->x), kp7_base(c));
  }
  for (unsigned t = 0; t < p->k; t++) {
    int c = kp7_code(kmer[t]);
    r = fadd(fmul(r, p->x), kp7_base(c));
  }
  *F = f;
  *R = fneg(r);
  return 1;
}

kp7_fe kp7_sym(const kp7_params* p, kp7_fe F, kp7_fe R) { return fmul_p(p, fadd(F, p->lam), fadd(R, p->lam)); }

/* ---- rolling ---- */

void kp7_start(kp7_roller* r, const kp7_params* p, const char* seq, size_t len)
{
  r->p = p; r->seq = seq; r->len = len; r->pos = 0; r->valid = 0; r->resume = 0;
}

static int restart(kp7_roller* r, size_t from)
{
  const unsigned k = r->p->k;
  size_t run = 0, i = from;
  while (run < k) {
    if (i >= r->len) { r->valid = 0; return 0; }
    run = CODE[(unsigned char)r->seq[i]] < 0 ? 0 : run + 1;
    i++;
  }
  r->pos = i - k;
  kp7_direct(r->p, r->seq + r->pos, &r->F, &r->R);
  r->valid = 1;
  return 1;
}

int kp7_roll(kp7_roller* r)
{
  if (!r->valid) return restart(r, r->pos);
  const kp7_params* p = r->p;
  size_t in = r->pos + p->k;
  if (in >= r->len) { r->valid = 0; r->pos = r->len; return 0; }
  int ci = CODE[(unsigned char)r->seq[in]];
  if (ci < 0) return restart(r, in + 1);
  int co = CODE[(unsigned char)r->seq[r->pos]];
  if (p->variant == KP7_XREAL) {
    r->F = fadd(fadd(fmul_real(r->F, p->xinv.re), p->out_f[co]), p->in_f[ci]);
    r->R = fadd(fadd(fmul_real(r->R, p->x.re), p->out_r[co]), p->in_r[ci]);
  } else {
    r->F = fadd(fadd(fmul_p(p, r->F, p->xinv), p->out_f[co]), p->in_f[ci]);
    r->R = fadd(fadd(fmul_p(p, r->R, p->x), p->out_r[co]), p->in_r[ci]);
  }
  r->pos++;
  return 1;
}

int kp7_roll64(kp7_roller* r, int mode, uint64_t* out)
{
  if (!kp7_roll(r)) return 0;
  switch (mode) {
  case KP7_OUT_FWD: *out = fold(r->F); break;
  case KP7_OUT_CANON_MIN: { uint64_t a = fold(r->F), b = fold(r->R); *out = a < b ? a : b; } break;
  default: *out = fold(fmul_p(r->p, fadd(r->F, r->p->lam), fadd(r->R, r->p->lam))); break;
  }
  return 1;
}

/* ---- symmetric-only roller with addition-side options (see kp7.h) ---- */

static inline uint64_t addq_lazy(uint64_t a, uint64_t b) { return a + b; } /* a, b < q: result < 2q */
static inline kp7_fe fadd_lazy(kp7_fe a, kp7_fe b) { kp7_fe r = { addq_lazy(a.re, b.re), addq_lazy(a.im, b.im) }; return r; }
/* product of two elements with parts < 4q: ac, bd, ad, bc < 16Q^2, so ac + (16Q^2 - bd) and ad + bc are < 32Q^2 < 2^127 */
static inline kp7_fe fmul_wide(kp7_fe a, kp7_fe b)
{
  const unsigned __int128 Q2x16 = (unsigned __int128)16 * Q * Q;
  unsigned __int128 ac = (unsigned __int128)a.re * b.re, bd = (unsigned __int128)a.im * b.im;
  unsigned __int128 ad = (unsigned __int128)a.re * b.im, bc = (unsigned __int128)a.im * b.re;
  kp7_fe r = { red2(ac + (Q2x16 - bd)), red2(ad + bc) };
  return r;
}

int kp7_roll_sym(kp7_roller* r, uint64_t* out)
{
  const kp7_params* p = r->p;
  const int opt = p->opt;
  const int lam = opt & KP7_OPT_LAMSTATE, merged = opt & KP7_OPT_MERGED;
  const int lazy = (opt & KP7_OPT_LAZY) && p->variant == KP7_XREAL;
  if (!r->valid) {
    if (!restart(r, r->pos)) return 0;
    r->G = lam ? fadd(r->F, p->lam) : r->F;
    r->H = lam ? fadd(r->R, p->lam) : r->R;
  } else {
    size_t in = r->pos + p->k;
    if (in >= r->len) { r->valid = 0; r->pos = r->len; return 0; }
    int ci = CODE[(unsigned char)r->seq[in]];
    if (ci < 0) {
      if (!restart(r, in + 1)) return 0;
      r->G = lam ? fadd(r->F, p->lam) : r->F;
      r->H = lam ? fadd(r->R, p->lam) : r->R;
    } else {
      int co = CODE[(unsigned char)r->seq[r->pos]];
      kp7_fe g, h;
      if (p->variant == KP7_XREAL) { g = fmul_real(r->G, p->xinv.re); h = fmul_real(r->H, p->x.re); } /* G, H < 2q ok */
      else { g = fmul_p(p, r->G, p->xinv); h = fmul_p(p, r->H, p->x); }
      if (merged) {
        const kp7_fe* tf = lam ? p->tfl : p->tf;
        const kp7_fe* tr = lam ? p->trl : p->tr;
        if (lazy) { g = fadd_lazy(g, tf[co * 4 + ci]); h = fadd_lazy(h, tr[co * 4 + ci]); }
        else { g = fadd(g, tf[co * 4 + ci]); h = fadd(h, tr[co * 4 + ci]); }
      } else if (lazy) {
        /* two lazy additions would leave [0, 3q); reduce once in the middle to stay below 2q */
        g = fadd_lazy(fadd(g, p->out_f[co]), p->in_f[ci]);
        h = fadd_lazy(fadd(h, p->out_r[co]), p->in_r[ci]);
        if (lam) { g.re = g.re >= Q ? g.re - Q : g.re; g.im = g.im >= Q ? g.im - Q : g.im;
                   h.re = h.re >= Q ? h.re - Q : h.re; h.im = h.im >= Q ? h.im - Q : h.im;
                   g = fadd_lazy(g, p->dl_f); h = fadd_lazy(h, p->dl_r); }
      } else {
        g = fadd(fadd(g, p->out_f[co]), p->in_f[ci]);
        h = fadd(fadd(h, p->out_r[co]), p->in_r[ci]);
        if (lam) { g = fadd(g, p->dl_f); h = fadd(h, p->dl_r); }
      }
      r->G = g; r->H = h;
      r->pos++;
    }
  }
  kp7_fe a = lam ? r->G : (lazy ? fadd_lazy(r->G, p->lam) : fadd(r->G, p->lam));
  kp7_fe b = lam ? r->H : (lazy ? fadd_lazy(r->H, p->lam) : fadd(r->H, p->lam));
  /* lazy: a, b < 2q with the lam state, < 3q when lam is added lazily to a < 2q state; fmul_wide takes < 4q */
  *out = fold(lazy ? fmul_wide(a, b) : fmul_p(p, a, b));
  return 1;
}

/* ---- batch symmetric hashing (see kp7.h) ---- */

size_t kp7_sym_batch(kp7_roller* r, uint64_t* out, size_t max)
{
  const kp7_params* p = r->p;
  const char* seq = r->seq;
  const size_t len = r->len, k = p->k;
  size_t n = 0;
  if (max == 0) return 0;
  for (;;) {
    if (!r->valid) {
      /* only reached with n == 0, so the outputs of one call stay consecutive */
      if (!restart(r, r->resume)) return 0; /* true end of the sequence */
      r->G = fadd(r->F, p->lam);
      r->H = fadd(r->R, p->lam);
      out[n++] = fold(fmul_p(p, r->G, r->H));
    }
    size_t pos = r->pos;
    kp7_fe G = r->G, H = r->H;
    int stop_n = 0; /* stopped at a window with N */
    if (p->variant == KP7_XREAL) {
      const uint64_t xi = p->xinv.re, x = p->x.re;
      uint64_t gr = G.re, gi = G.im, hr = H.re, hi = H.im; /* each < 2q (lazy) */
      while (n < max) {
        size_t in = pos + k;
        if (in >= len) { r->valid = 0; r->resume = len; break; }
        int ci = CODE[(unsigned char)seq[in]];
        if (ci < 0) { r->valid = 0; r->resume = in + 1; stop_n = 1; break; }
        int idx = CODE[(unsigned char)seq[pos]] * 4 + ci;
        gr = mulq(gr, xi) + p->tfl[idx].re; gi = mulq(gi, xi) + p->tfl[idx].im; /* < q + q */
        hr = mulq(hr, x) + p->trl[idx].re;  hi = mulq(hi, x) + p->trl[idx].im;
        kp7_fe a = { gr, gi }, b = { hr, hi };
        out[n++] = fold(fmul_wide(a, b));
        pos++;
      }
      G.re = gr; G.im = gi; H.re = hr; H.im = hi;
    } else {
      while (n < max) {
        size_t in = pos + k;
        if (in >= len) { r->valid = 0; r->resume = len; break; }
        int ci = CODE[(unsigned char)seq[in]];
        if (ci < 0) { r->valid = 0; r->resume = in + 1; stop_n = 1; break; }
        int idx = CODE[(unsigned char)seq[pos]] * 4 + ci;
        G = fadd(fmul_p(p, G, p->xinv), p->tfl[idx]);
        H = fadd(fmul_p(p, H, p->x), p->trl[idx]);
        out[n++] = fold(fmul_p(p, G, H));
        pos++;
      }
    }
    /* pos is always the last k-mer written; after a stop, r->resume is where the next restart searches */
    r->pos = pos;
    r->G = G; r->H = H;
    if (stop_n && n == 0) continue; /* the first incoming base was N: restart now rather than return 0 */
    return n;
  }
}

/* ---- interleaved symmetric hashing (step 2; see kp7.h) ---- */

#ifdef KP7_MUTATE_SEQ /* test-only: the last lane's remainder lands one slot early; `make mutation` expects FAIL */
#define KP7_MUTATE_REM_OFF -1
#else
#define KP7_MUTATE_REM_OFF
#endif

/* lane start: canonical G = F + lam, H = R + lam of the k-mer at s (no N, checked by the caller) */
static inline void lane_init(const kp7_params* p, const char* s, kp7_fe* G, kp7_fe* H)
{
  kp7_fe F, R;
  kp7_direct(p, s, &F, &R);
  *G = fadd(F, p->lam);
  *H = fadd(R, p->lam);
}

/* x real, L lanes. s: first k-mer of an N-free stretch holding nk >= L k-mers. Lanes cover [l*seg, (l+1)*seg); the
 * last lane also takes the remainder. A lane advances only while its next k-mer is still in its stretch, so it
 * never reads past s[nk + k - 2]. */
#define KP7_DEF_BLOCK_REAL(L)                                                                                         \
  static void block_real_##L(const kp7_params* p, const char* s, size_t nk, uint64_t* out)                           \
  {                                                                                                                  \
    const size_t k = p->k, seg = nk / L, rem = nk - seg * L;                                                         \
    const uint64_t xi = p->xinv.re, x = p->x.re;                                                                     \
    const kp7_fe *tf = p->tfl, *tr = p->trl;                                                                         \
    uint64_t gr[L], gi[L], hr[L], hi[L];                                                                             \
    for (int l = 0; l < L; l++) {                                                                                    \
      kp7_fe G, H;                                                                                                   \
      lane_init(p, s + l * seg, &G, &H);                                                                             \
      gr[l] = G.re; gi[l] = G.im; hr[l] = H.re; hi[l] = H.im;                                                        \
    }                                                                                                                \
    for (size_t i = 0; i < seg; i++) {                                                                               \
      for (int l = 0; l < L; l++) {                                                                                  \
        size_t q = l * seg + i;                                                                                      \
        kp7_fe a = { gr[l], gi[l] }, b = { hr[l], hi[l] };                                                           \
        out[q] = fold(fmul_wide(a, b));                                                                              \
        if (i + 1 < seg) {                                                                                           \
          int idx = CODE[(unsigned char)s[q]] * 4 + CODE[(unsigned char)s[q + k]];                                   \
          gr[l] = mulq(gr[l], xi) + tf[idx].re; gi[l] = mulq(gi[l], xi) + tf[idx].im;                                \
          hr[l] = mulq(hr[l], x) + tr[idx].re;  hi[l] = mulq(hi[l], x) + tr[idx].im;                                 \
        }                                                                                                            \
      }                                                                                                              \
    }                                                                                                                \
    /* the last lane continues through the remainder */                                                              \
    for (size_t j = 0, q = L * seg - 1; j < rem; j++, q++) {                                                         \
      int idx = CODE[(unsigned char)s[q]] * 4 + CODE[(unsigned char)s[q + k]];                                       \
      gr[L - 1] = mulq(gr[L - 1], xi) + tf[idx].re; gi[L - 1] = mulq(gi[L - 1], xi) + tf[idx].im;                    \
      hr[L - 1] = mulq(hr[L - 1], x) + tr[idx].re;  hi[L - 1] = mulq(hi[L - 1], x) + tr[idx].im;                     \
      kp7_fe a = { gr[L - 1], gi[L - 1] }, b = { hr[L - 1], hi[L - 1] };                                             \
      out[q + 1 KP7_MUTATE_REM_OFF] = fold(fmul_wide(a, b));                                                          \
    }                                                                                                                \
  }
KP7_DEF_BLOCK_REAL(1)
KP7_DEF_BLOCK_REAL(2)
KP7_DEF_BLOCK_REAL(3)
KP7_DEF_BLOCK_REAL(4)

/* x complex, L lanes, canonical arithmetic (lazy reduction does not apply) */
#define KP7_DEF_BLOCK_CPLX(L)                                                                                         \
  static void block_cplx_##L(const kp7_params* p, const char* s, size_t nk, uint64_t* out)                           \
  {                                                                                                                  \
    const size_t k = p->k, seg = nk / L, rem = nk - seg * L;                                                         \
    kp7_fe G[L], H[L];                                                                                               \
    for (int l = 0; l < L; l++) lane_init(p, s + l * seg, &G[l], &H[l]);                                             \
    for (size_t i = 0; i < seg; i++) {                                                                               \
      for (int l = 0; l < L; l++) {                                                                                  \
        size_t q = l * seg + i;                                                                                      \
        out[q] = fold(fmul_p(p, G[l], H[l]));                                                                        \
        if (i + 1 < seg) {                                                                                           \
          int idx = CODE[(unsigned char)s[q]] * 4 + CODE[(unsigned char)s[q + k]];                                   \
          G[l] = fadd(fmul_p(p, G[l], p->xinv), p->tfl[idx]);                                                        \
          H[l] = fadd(fmul_p(p, H[l], p->x), p->trl[idx]);                                                           \
        }                                                                                                            \
      }                                                                                                              \
    }                                                                                                                \
    for (size_t j = 0, q = L * seg - 1; j < rem; j++, q++) {                                                         \
      int idx = CODE[(unsigned char)s[q]] * 4 + CODE[(unsigned char)s[q + k]];                                       \
      G[L - 1] = fadd(fmul_p(p, G[L - 1], p->xinv), p->tfl[idx]);                                                    \
      H[L - 1] = fadd(fmul_p(p, H[L - 1], p->x), p->trl[idx]);                                                       \
      out[q + 1] = fold(fmul_p(p, G[L - 1], H[L - 1]));                                                              \
    }                                                                                                                \
  }
KP7_DEF_BLOCK_CPLX(1)
KP7_DEF_BLOCK_CPLX(2)

size_t kp7_sym_seq(const kp7_params* p, const char* seq, size_t len, int lanes, uint64_t* buf, size_t cap,
                   kp7_emit_fn emit, void* ctx)
{
  if (!code_ready) code_init();
  const size_t k = p->k;
  const int cplx = p->variant != KP7_XREAL;
  if (lanes < 1) lanes = 1;
  if (lanes > (cplx ? 2 : 4)) lanes = cplx ? 2 : 4;
  size_t total = 0, i = 0;
  while (i < len) {
    while (i < len && CODE[(unsigned char)seq[i]] < 0) i++; /* skip N */
    size_t a = i;
    while (i < len && CODE[(unsigned char)seq[i]] >= 0) i++;
    size_t run = i - a;
    if (run < k) continue;
    size_t nk = run - k + 1;
    for (size_t off = 0; off < nk; off += cap) {
      size_t n = nk - off < cap ? nk - off : cap;
      const char* s = seq + a + off;
      int L = (size_t)lanes <= n ? lanes : 1; /* every lane needs at least one k-mer */
      if (cplx) { if (L == 2) block_cplx_2(p, s, n, buf); else block_cplx_1(p, s, n, buf); }
      else switch (L) {
        case 4: block_real_4(p, s, n, buf); break;
        case 3: block_real_3(p, s, n, buf); break;
        case 2: block_real_2(p, s, n, buf); break;
        default: block_real_1(p, s, n, buf); break;
      }
      emit(ctx, a + off, buf, n);
      total += n;
    }
  }
  return total;
}

/* ---- hybrid: cheap symmetric gate, provable hash for kept k-mers (see kp7.h) ---- */

static inline uint64_t canon_q(uint64_t v) { return v >= Q ? v - Q : v; } /* [0, 2q) -> [0, q) */
static inline uint64_t mix64(uint64_t x)
{
  x ^= x >> 32;
  x *= 0xd6e8feb86659fd93ULL; /* odd multiplier (a known good 64-bit mixing constant) */
  return x ^ (x >> 32);
}
#ifdef KP7_MUTATE_GATE /* test-only: break the symmetry of the gate; `make mutation` expects FAIL */
static inline uint64_t gate_of(uint64_t g, uint64_t h) { return mix64(g) + h; }
#else
static inline uint64_t gate_of(uint64_t g, uint64_t h) { return mix64(g) + mix64(h); }
#endif

uint64_t kp7_gate(const kp7_params* p, kp7_fe F, kp7_fe R)
{
  return gate_of(fold(fadd(F, p->lam)), fold(fadd(R, p->lam)));
}

static inline uint64_t rot32(uint64_t v) { return (v << 32) | (v >> 32); } /* fold's rotation */

/* MODE 0: out[q] = gate for every k-mer. MODE 1: keep gate < thr, with the provable hash. x real, lazy state. */
#define KP7_DEF_HYBRID(NAME, L, MODE)                                                                                 \
  static size_t NAME(const kp7_params* p, const char* s, size_t nk, size_t base, uint64_t thr, uint64_t* out,         \
                     size_t* opos)                                                                                    \
  {                                                                                                                  \
    const size_t k = p->k, seg = nk / L, rem = nk - seg * L;                                                         \
    const uint64_t xi = p->xinv.re, x = p->x.re;                                                                     \
    const kp7_fe *tf = p->tfl, *tr = p->trl;                                                                         \
    uint64_t gr[L], gi[L], hr[L], hi[L];                                                                             \
    /* MODE 1: lane l writes its kept k-mers into its own region [l*seg, ...) so each block's output stays in */     \
    /* position order; the regions are joined at the end. The last lane's region also holds the remainder. */        \
    size_t kl[L] = { 0 };                                                                                            \
    (void)thr; (void)opos; (void)base;                                                                               \
    for (int l = 0; l < L; l++) {                                                                                    \
      kp7_fe G, H;                                                                                                   \
      lane_init(p, s + l * seg, &G, &H);                                                                             \
      gr[l] = G.re; gi[l] = G.im; hr[l] = H.re; hi[l] = H.im;                                                        \
    }                                                                                                                \
    for (size_t i = 0; i < seg; i++) {                                                                               \
      for (int l = 0; l < L; l++) {                                                                                  \
        size_t q = l * seg + i;                                                                                      \
        uint64_t g = canon_q(gr[l]) ^ rot32(canon_q(gi[l])), h = canon_q(hr[l]) ^ rot32(canon_q(hi[l]));             \
        uint64_t gt = gate_of(g, h);                                                                                 \
        if (MODE == 0) out[q] = gt;                                                                                  \
        else if (gt < thr) {                                                                                         \
          kp7_fe a = { gr[l], gi[l] }, b = { hr[l], hi[l] };                                                         \
          size_t w = l * seg + kl[l]++;                                                                              \
          out[w] = fold(fmul_wide(a, b));                                                                            \
          opos[w] = base + q;                                                                                        \
        }                                                                                                            \
        if (i + 1 < seg) {                                                                                           \
          int idx = CODE[(unsigned char)s[q]] * 4 + CODE[(unsigned char)s[q + k]];                                   \
          gr[l] = mulq(gr[l], xi) + tf[idx].re; gi[l] = mulq(gi[l], xi) + tf[idx].im;                                \
          hr[l] = mulq(hr[l], x) + tr[idx].re;  hi[l] = mulq(hi[l], x) + tr[idx].im;                                 \
        }                                                                                                            \
      }                                                                                                              \
    }                                                                                                                \
    for (size_t j = 0, q = L * seg - 1; j < rem; j++, q++) {                                                         \
      int idx = CODE[(unsigned char)s[q]] * 4 + CODE[(unsigned char)s[q + k]];                                       \
      gr[L - 1] = mulq(gr[L - 1], xi) + tf[idx].re; gi[L - 1] = mulq(gi[L - 1], xi) + tf[idx].im;                    \
      hr[L - 1] = mulq(hr[L - 1], x) + tr[idx].re;  hi[L - 1] = mulq(hi[L - 1], x) + tr[idx].im;                     \
      uint64_t g = canon_q(gr[L - 1]) ^ rot32(canon_q(gi[L - 1])), h = canon_q(hr[L - 1]) ^ rot32(canon_q(hi[L - 1]));\
      uint64_t gt = gate_of(g, h);                                                                                   \
      if (MODE == 0) out[q + 1] = gt;                                                                                \
      else if (gt < thr) {                                                                                           \
        kp7_fe a = { gr[L - 1], gi[L - 1] }, b = { hr[L - 1], hi[L - 1] };                                           \
        size_t w = (L - 1) * seg + kl[L - 1]++;                                                                      \
        out[w] = fold(fmul_wide(a, b));                                                                              \
        opos[w] = base + q + 1;                                                                                      \
      }                                                                                                              \
    }                                                                                                                \
    if (MODE == 0) return nk;                                                                                        \
    size_t kept = kl[0];                                                                                             \
    for (int l = 1; l < L; l++) {                                                                                    \
      memmove(out + kept, out + l * seg, kl[l] * sizeof *out);                                                       \
      memmove(opos + kept, opos + l * seg, kl[l] * sizeof *opos);                                                    \
      kept += kl[l];                                                                                                 \
    }                                                                                                                \
    return kept;                                                                                                     \
  }
KP7_DEF_HYBRID(gate_1, 1, 0)
KP7_DEF_HYBRID(gate_2, 2, 0)
KP7_DEF_HYBRID(sample_1, 1, 1)
KP7_DEF_HYBRID(sample_2, 2, 1)

/* walk N-free runs in blocks, as kp7_sym_seq does */
#define KP7_FOR_BLOCKS(BODY)                                                                                          \
  size_t i = 0;                                                                                                      \
  while (i < len) {                                                                                                  \
    while (i < len && CODE[(unsigned char)seq[i]] < 0) i++;                                                          \
    size_t a = i;                                                                                                    \
    while (i < len && CODE[(unsigned char)seq[i]] >= 0) i++;                                                         \
    size_t run = i - a;                                                                                              \
    if (run < k) continue;                                                                                           \
    size_t nk = run - k + 1;                                                                                         \
    for (size_t off = 0; off < nk; off += cap) {                                                                     \
      size_t n = nk - off < cap ? nk - off : cap;                                                                    \
      const char* s = seq + a + off;                                                                                 \
      int L = lanes >= 2 && n >= 2 ? 2 : 1;                                                                          \
      BODY                                                                                                           \
    }                                                                                                                \
  }

size_t kp7_gate_seq(const kp7_params* p, const char* seq, size_t len, int lanes, uint64_t* buf, size_t cap,
                    kp7_emit_fn emit, void* ctx)
{
  if (!code_ready) code_init();
  if (p->variant != KP7_XREAL) return 0;
  const size_t k = p->k;
  size_t total = 0;
  KP7_FOR_BLOCKS({
    if (L == 2) gate_2(p, s, n, 0, 0, buf, NULL); else gate_1(p, s, n, 0, 0, buf, NULL);
    emit(ctx, a + off, buf, n);
    total += n;
  })
  return total;
}

size_t kp7_sym_sample(const kp7_params* p, const char* seq, size_t len, int lanes, uint64_t threshold, size_t* pos_buf,
                      uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx)
{
  if (!code_ready) code_init();
  if (p->variant != KP7_XREAL) return 0;
  const size_t k = p->k;
  size_t kept_total = 0;
  KP7_FOR_BLOCKS({
    size_t kept = L == 2 ? sample_2(p, s, n, a + off, threshold, hash_buf, pos_buf)
                         : sample_1(p, s, n, a + off, threshold, hash_buf, pos_buf);
    if (kept) emit(ctx, pos_buf, hash_buf, kept);
    kept_total += kept;
  })
  return kept_total;
}
