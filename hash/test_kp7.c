/* Tests for kp7: rolling == direct, strand symmetry, exhaustive collision counts. Exit 0 iff all pass. */
#include "kp7.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fails = 0;
#define CHECK(cond, ...) do { if (!(cond)) { fails++; if (fails < 10) { printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } } while (0)

static int fe_eq(kp7_fe a, kp7_fe b) { return a.re == b.re && a.im == b.im; }
static void rc(const char* s, char* out, unsigned k)
{
  for (unsigned t = 0; t < k; t++) {
    char c = s[k - 1 - t];
    out[t] = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N';
  }
  out[k] = 0;
}

/* 1 + 2: rolling equals direct at every valid window, skips exactly the N windows, and sym(s) == sym(rc(s)) */
static void test_rolling(int variant, unsigned k, size_t len, uint64_t seed)
{
  kp7_params p;
  kp7_init(&p, k, variant, seed);
  char* s = malloc(len + 1);
  srand((unsigned)seed);
  for (size_t i = 0; i < len; i++) s[i] = "ACGTACGTACGTACGTacgtN"[rand() % 21]; /* ~5% N, some lowercase */
  s[len] = 0;
  size_t expect = 0;
  for (size_t i = 0; i + k <= len; i++) {
    int ok = 1;
    for (unsigned t = 0; t < k; t++) if (kp7_code(s[i + t]) < 0) { ok = 0; break; }
    expect += ok;
  }
  kp7_roller r;
  kp7_start(&r, &p, s, len);
  size_t n = 0, last = (size_t)-1;
  char buf[256], rcb[256];
  while (kp7_roll(&r)) {
    kp7_fe F, R;
    CHECK(last == (size_t)-1 || r.pos > last, "positions not increasing");
    last = r.pos;
    CHECK(kp7_direct(&p, s + r.pos, &F, &R), "direct says invalid at %zu", r.pos);
    CHECK(fe_eq(F, r.F) && fe_eq(R, r.R), "variant %d k %u: rolling != direct at %zu", variant, k, r.pos);
    for (unsigned t = 0; t < k; t++) buf[t] = s[r.pos + t] & ~0x20; /* upper-case */
    buf[k] = 0;
    rc(buf, rcb, k);
    kp7_fe F2, R2;
    kp7_direct(&p, rcb, &F2, &R2);
    CHECK(fe_eq(F2, R) && fe_eq(R2, F), "F(rc s) != R(s) at %zu", r.pos);
    CHECK(fe_eq(kp7_sym(&p, F, R), kp7_sym(&p, F2, R2)), "sym(s) != sym(rc s) at %zu", r.pos);
    n++;
  }
  CHECK(n == expect, "variant %d k %u: rolled %zu windows, expected %zu", variant, k, n, expect);
  /* kp7_roll64 agrees with kp7_roll + fold */
  kp7_roller a, b;
  kp7_start(&a, &p, s, len);
  kp7_start(&b, &p, s, len);
  uint64_t o;
  while (kp7_roll64(&a, KP7_OUT_SYM, &o)) {
    kp7_roll(&b);
    CHECK(o == kp7_fold(kp7_sym(&p, b.F, b.R)), "roll64 sym mismatch");
  }
  free(s);
}

static int cmp_u64(const void* a, const void* b)
{
  uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
  return x < y ? -1 : x > y;
}
static int cmp_fe(const void* a, const void* b)
{
  const kp7_fe *x = a, *y = b;
  if (x->re != y->re) return x->re < y->re ? -1 : 1;
  return x->im < y->im ? -1 : x->im > y->im;
}
static size_t distinct_u64(uint64_t* v, size_t n)
{
  qsort(v, n, sizeof *v, cmp_u64);
  size_t d = n ? 1 : 0;
  for (size_t i = 1; i < n; i++) d += v[i] != v[i - 1];
  return d;
}
static size_t distinct_fe(kp7_fe* v, size_t n)
{
  qsort(v, n, sizeof *v, cmp_fe);
  size_t d = n ? 1 : 0;
  for (size_t i = 1; i < n; i++) d += !fe_eq(v[i], v[i - 1]);
  return d;
}

/* 3: exhaustive at small k: F injective; sym has exactly one value per double-stranded class (also after fold) */
static void test_exhaustive(int variant, unsigned k, size_t classes)
{
  kp7_params p;
  kp7_init(&p, k, variant, 12345);
  size_t n = (size_t)1 << (2 * k);
  kp7_fe* fs = malloc(n * sizeof *fs);
  kp7_fe* ss = malloc(n * sizeof *ss);
  uint64_t* folded = malloc(n * sizeof *folded);
  char kmer[32];
  for (size_t v = 0; v < n; v++) {
    for (unsigned t = 0; t < k; t++) kmer[t] = "ACGT"[(v >> (2 * t)) & 3];
    kp7_fe F, R;
    kp7_direct(&p, kmer, &F, &R);
    fs[v] = F;
    ss[v] = kp7_sym(&p, F, R);
    folded[v] = kp7_fold(ss[v]);
  }
  size_t df = distinct_fe(fs, n), ds = distinct_fe(ss, n), dfold = distinct_u64(folded, n);
  printf("  variant %d k=%u: %zu k-mers, F distinct %zu, sym distinct %zu, folded sym distinct %zu (classes %zu)\n",
         variant, k, n, df, ds, dfold, classes);
  CHECK(df == n, "F not injective at k=%u", k);
  CHECK(ds == classes, "sym classes %zu != %zu at k=%u", ds, classes, k);
  CHECK(dfold == classes, "folded sym classes %zu != %zu at k=%u", dfold, classes, k);
  free(fs); free(ss); free(folded);
}

/* field sanity */
static void test_field(void)
{
  kp7_fe a = { 123456789012345ULL, 987654321098765ULL };
  kp7_fe ai = kp7_inv(a);
  kp7_fe one = kp7_mul(a, ai);
  CHECK(one.re == 1 && one.im == 0, "inverse");
  kp7_fe i = { 0, 1 };
  kp7_fe m1 = kp7_mul(i, i);
  CHECK(m1.re == KP7_Q - 1 && m1.im == 0, "i^2 != -1");
  CHECK(kp7_mulq(KP7_Q - 1, KP7_Q - 1) == 1, "(-1)^2 != 1");
}

int main(void)
{
  test_field();
  for (int variant = 0; variant < 2; variant++) {
    test_rolling(variant, 31, 200000, 7 + variant);
    test_rolling(variant, 8, 50000, 9 + variant);
    test_rolling(variant, 1, 1000, 11 + variant);
    test_exhaustive(variant, 7, 8192);   /* 4^7 / 2: odd k has no RC palindromes */
    test_exhaustive(variant, 8, 32896);  /* (4^8 + 4^4) / 2 */
  }
  printf(fails ? "FAILED: %d check(s)\n" : "all tests passed\n", fails);
  return fails != 0;
}
