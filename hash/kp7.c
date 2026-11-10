/* KP-6 / KP-7 rolling hashes; see kp7.h. */
#include "kp7.h"

#define Q KP7_Q
static const unsigned __int128 Q2 = (unsigned __int128)KP7_Q * KP7_Q;

/* ---- field arithmetic (inline in this translation unit, exported wrappers below) ---- */

static inline uint64_t red(unsigned __int128 p) /* p < 2^123 */
{
  uint64_t r = ((uint64_t)p & Q) + (uint64_t)(p >> 61);
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

kp7_fe kp7_sym(const kp7_params* p, kp7_fe F, kp7_fe R) { return fmul(fadd(F, p->lam), fadd(R, p->lam)); }

/* ---- rolling ---- */

void kp7_start(kp7_roller* r, const kp7_params* p, const char* seq, size_t len)
{
  r->p = p; r->seq = seq; r->len = len; r->pos = 0; r->valid = 0;
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
    r->F = fadd(fadd(fmul(r->F, p->xinv), p->out_f[co]), p->in_f[ci]);
    r->R = fadd(fadd(fmul(r->R, p->x), p->out_r[co]), p->in_r[ci]);
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
  default: *out = fold(fmul(fadd(r->F, r->p->lam), fadd(r->R, r->p->lam))); break;
  }
  return 1;
}
