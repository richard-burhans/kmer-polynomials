/* KP-5: reverse-complement palindromes and inverted repeats as a polynomial identity.
 *
 * Encoding as hash/kp7e.h (Eisenstein, one element of F_q, q = 2^61 - 1): A = 1, C = w, G = -w, T = -1, complement is
 * negation. For an a-mer u, F(u) = sum_t u_t x^t and R(u) = -x^{a-1} p_u(1/x) = p_{rc(u)}(x), the forward polynomial of
 * rc(u). Hence:
 *   u is an RC palindrome          <=>  p_u = p_{rc u}        =>  F(u) == R(u)
 *   right arm v is the RC of left u <=> p_u = p_{rc v}        =>  F(u) == R(v)
 * and a false positive needs a nonzero polynomial of degree <= a-1 to vanish at the random x: probability <= (a-1)/q
 * per compared pair (Schwartz-Zippel). Odd-length RC palindromes cannot exist (the middle base would equal its own
 * complement), so `palin` is only meaningful for even k.
 *
 * Methods (identical loops; only the per-position key differs):
 *   hash  F and R of every a-mer, rolled with two 61-bit multiplies per position (merged 16-entry tables)
 *   pack  2-bit packed fwd and rc codes (unsigned __int128, so a <= 64), rolled with shifts
 * Modes:
 *   ir FASTA RECORD palin K METHOD          count windows of length K that equal their reverse complement
 *   ir FASTA RECORD inv A GMAX METHOD       perfect inverted repeats: arm length A, spacer 0..GMAX (j - i - A = g),
 *                                           reported as (i, g); prints the count and an order-independent checksum
 * RECORD may be "-" for the first record. METHOD is hash or pack. `pack` is the exact reference (no false positives,
 * no misses), so equal counts and checksums from `hash` mean it found exactly the same repeats.
 */
#include "../../hash/kp7e.h"
#include "fasta.h"
#include <stdint.h>
#include <time.h>

#define Q KP7E_Q
static inline uint64_t red(unsigned __int128 v)
{
  uint64_t r = ((uint64_t)v & Q) + (uint64_t)(v >> 61);
  r = (r & Q) + (r >> 61);
  return r >= Q ? r - Q : r;
}
static inline uint64_t mulq(uint64_t a, uint64_t b) { return red((unsigned __int128)a * b); }
static inline uint64_t addq(uint64_t a, uint64_t b) { uint64_t r = a + b; return r >= Q ? r - Q : r; }
static inline uint64_t subq(uint64_t a, uint64_t b) { return a >= b ? a - b : a + Q - b; }

static int8_t CODE[256];
static void code_init(void)
{
  memset(CODE, -1, sizeof CODE);
  CODE['A'] = 0; CODE['C'] = 1; CODE['G'] = 2; CODE['T'] = 3; /* complement 3 - c, as kp7e */
}

/* key arrays: kf[i], kr[i] for every a-mer start i; valid[i] = 0 if the window holds a non-ACGT base */
typedef struct { size_t n; uint8_t* valid; uint64_t* hf; uint64_t* hr; unsigned __int128* pf; unsigned __int128* pr; } keys;

static void keys_hash(const char* s, size_t len, unsigned a, keys* K)
{
  kp7e_params p;
  kp7e_init(&p, a, 0x6b70355ULL); /* seed fixed so runs are reproducible */
  uint64_t xa = mulq(p.xk1, p.x); /* x^a */
  uint64_t tf[16], tr[16];
  for (int co = 0; co < 4; co++)
    for (int ci = 0; ci < 4; ci++) {
      /* F' = F x^{-1} - s_out x^{-1} + s_in x^{a-1};  R' = R x + s_out x^a - s_in */
      tf[co * 4 + ci] = addq(subq(0, mulq(p.val[co], p.xinv)), mulq(p.val[ci], p.xk1));
#ifdef KP5_MUTATE /* test-only: wrong sign on the incoming base in the RC roll */
      tr[co * 4 + ci] = addq(mulq(p.val[co], xa), p.val[ci]);
#else
      tr[co * 4 + ci] = subq(mulq(p.val[co], xa), p.val[ci]);
#endif
    }
  size_t n = len >= a ? len - a + 1 : 0;
  K->n = n;
  size_t run = 0;
  uint64_t F = 0, R = 0;
  for (size_t j = 0; j < len; j++) { /* window [j-a+1, j] ends at j */
    int c = CODE[(unsigned char)s[j]];
    if (c < 0) { run = 0; if (j + 1 >= a) K->valid[j + 1 - a] = 0; continue; }
    run++;
    if (run < a) { if (j + 1 >= a) K->valid[j + 1 - a] = 0; continue; }
    size_t i = j + 1 - a;
    if (run == a) { /* direct start */
      kp7e_direct(&p, s + i, &F, &R);
    } else {
      int co = CODE[(unsigned char)s[i - 1]];
      F = addq(mulq(F, p.xinv), tf[co * 4 + c]);
      R = addq(mulq(R, p.x), tr[co * 4 + c]);
    }
    K->valid[i] = 1; K->hf[i] = F; K->hr[i] = R;
  }
}

static void keys_pack(const char* s, size_t len, unsigned a, keys* K)
{
  size_t n = len >= a ? len - a + 1 : 0;
  K->n = n;
  unsigned __int128 mask = a == 64 ? ~(unsigned __int128)0 : (((unsigned __int128)1 << (2 * a)) - 1);
  unsigned __int128 fw = 0, rc = 0;
  size_t run = 0;
  for (size_t j = 0; j < len; j++) {
    int c = CODE[(unsigned char)s[j]];
    if (c < 0) { run = 0; fw = rc = 0; if (j + 1 >= a) K->valid[j + 1 - a] = 0; continue; }
    run++;
    /* fw: first base most significant; rc: complement of the newest base becomes the most significant digit */
    fw = ((fw << 2) | (unsigned)c) & mask;
    rc = (rc >> 2) | ((unsigned __int128)(3 - c) << (2 * (a - 1)));
    if (run < a) { if (j + 1 >= a) K->valid[j + 1 - a] = 0; continue; }
    size_t i = j + 1 - a;
    K->valid[i] = 1; K->pf[i] = fw; K->pr[i] = rc;
  }
}

static inline uint64_t mix64(uint64_t x) { x ^= x >> 31; x *= 0x7fb5d329728ea185ULL; x ^= x >> 27; x *= 0x81dadef4bc2dd44dULL; return x ^ (x >> 33); }

int main(int argc, char** argv)
{
  if (argc < 5) { fprintf(stderr, "usage: ir FASTA RECORD palin K hash|pack | ir FASTA RECORD inv A GMAX hash|pack\n"); return 2; }
  code_init();
  size_t len;
  char* s = load_record(argv[1], strcmp(argv[2], "-") ? argv[2] : NULL, &len);
  int palin = !strcmp(argv[3], "palin");
  unsigned a = (unsigned)atoi(argv[4]);
  unsigned gmax = palin ? 0 : (unsigned)atoi(argv[5]);
  const char* method = argv[palin ? 5 : 6];
  int hash = !strcmp(method, "hash");
  if (!hash && strcmp(method, "pack")) { fprintf(stderr, "method must be hash or pack\n"); return 2; }
  if (a < 1 || a > 64) { fprintf(stderr, "arm/k must be 1..64\n"); return 2; }
  keys K = { 0 };
  K.valid = calloc(len + 1, 1);
  if (hash) { K.hf = malloc((len + 1) * 8); K.hr = malloc((len + 1) * 8); }
  else { K.pf = malloc((len + 1) * 16); K.pr = malloc((len + 1) * 16); }
  double t0 = now_s();
  if (hash) keys_hash(s, len, a, &K); else keys_pack(s, len, a, &K);
  double t1 = now_s();
  uint64_t count = 0, sum = 0, xr = 0;
  if (palin) {
    for (size_t i = 0; i < K.n; i++) {
      if (!K.valid[i]) continue;
      int eq = hash ? K.hf[i] == K.hr[i] : K.pf[i] == K.pr[i];
      if (eq) { count++; uint64_t h = mix64(i); sum += h; xr ^= h; }
    }
  } else {
    /* right arm starts at j = i + a + g, g in [0, gmax]; compare F(left) with R(right) */
    for (size_t i = 0; i < K.n; i++) {
      if (!K.valid[i]) continue;
      for (unsigned g = 0; g <= gmax; g++) {
        size_t j = i + a + g;
        if (j >= K.n) break;
        if (!K.valid[j]) continue;
        int eq = hash ? K.hf[i] == K.hr[j] : K.pf[i] == K.pr[j];
        if (eq) { count++; uint64_t h = mix64(i * 1315423911ULL + g); sum += h; xr ^= h; }
      }
    }
  }
  double t2 = now_s();
  printf("mode=%s a=%u gmax=%u method=%s len=%zu windows=%zu count=%llu checksum=%016llx%016llx keys_s=%.3f scan_s=%.3f\n",
         palin ? "palin" : "inv", a, gmax, method, len, K.n, (unsigned long long)count, (unsigned long long)sum,
         (unsigned long long)xr, t1 - t0, t2 - t1);
  return 0;
}
