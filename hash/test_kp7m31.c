/* Tests for kp7m31 (CM31, p = 2^31 - 1): rolling == direct, strand symmetry, exhaustive class counts, and
 * scalar / AVX2 / AVX-512 sym_seq bit-identical to the reference roller (hashes AND positions). Exit 0 iff all pass.
 * Also: x and lam independent of each other (separate seeds), OS mode, ranges.
 */
#include "kp7m31.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fails = 0;
#define CHECK(cond, ...) do { if (!(cond)) { fails++; if (fails < 10) { printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } } while (0)

static uint64_t sm(uint64_t* s)
{
  uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static void rc(const char* s, char* out, unsigned k)
{
  for (unsigned t = 0; t < k; t++) {
    char c = s[k - 1 - t];
    out[t] = c == 'A' || c == 'a' ? 'T' : c == 'C' || c == 'c' ? 'G' : c == 'G' || c == 'g' ? 'C' : c == 'T' || c == 't' ? 'A' : 'N';
  }
  out[k] = 0;
}

/* 1 + 2: the reference roller equals direct evaluation at every window; sym(s) == sym(rc s) */
static void test_rolling(unsigned k, const char* seq, size_t len, const char* label)
{
  kp7m31_params p;
  kp7m31_init(&p, k, 0xC0FFEE + k, (0xC0FFEE + k) + 1);
  kp7m31_roller r;
  kp7m31_start(&r, &p, seq, len);
  size_t n = 0, valid = 0;
  for (size_t i = 0; i + k <= len; i++) {
    int ok = 1;
    for (unsigned t = 0; t < k; t++) ok &= kp7m31_code(seq[i + t]) >= 0;
    valid += ok;
  }
  char* w = malloc(k + 1); char* wr = malloc(k + 1);
  while (kp7m31_roll(&r)) {
    m31c F, R, F2, R2;
    memcpy(w, seq + r.pos, k); w[k] = 0;
    CHECK(kp7m31_direct(&p, w, &F, &R), "%s: rolled window %zu has N", label, r.pos);
    CHECK(F.re == r.F.re && F.im == r.F.im && R.re == r.R.re && R.im == r.R.im, "%s: rolling != direct at %zu", label, r.pos);
    rc(w, wr, k);
    kp7m31_direct(&p, wr, &F2, &R2);
    CHECK(F2.re == R.re && F2.im == R.im, "%s: F(rc s) != R(s) at %zu", label, r.pos);
    CHECK(kp7m31_sym(&p, F, R) == kp7m31_sym(&p, F2, R2), "%s: sym not strand-symmetric at %zu", label, r.pos);
    n++;
  }
  CHECK(n == valid, "%s: rolled %zu windows, %zu valid", label, n, valid);
  free(w); free(wr);
}

/* 3: one sym value per double-stranded class */
static int cmpu64(const void* a, const void* b) { uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b; return x < y ? -1 : x > y; }
static void test_exhaustive(unsigned k, size_t classes)
{
  kp7m31_params p;
  kp7m31_init(&p, k, 0xBEEF, 0xBEEF + 1);
  size_t nk = (size_t)1 << (2 * k), distinct = 0;
  uint64_t* v = malloc(nk * sizeof *v);
  char w[16];
  for (size_t code = 0; code < nk; code++) {
    for (unsigned t = 0; t < k; t++) w[t] = "ACGT"[(code >> (2 * t)) & 3];
    w[k] = 0;
    m31c F, R;
    kp7m31_direct(&p, w, &F, &R);
    v[code] = kp7m31_sym(&p, F, R);
  }
  qsort(v, nk, sizeof *v, cmpu64);
  for (size_t i = 0; i < nk; i++) distinct += i == 0 || v[i] != v[i - 1];
  printf("  m31 k=%u: %zu sym values for %zu classes\n", k, distinct, classes);
  CHECK(distinct == classes, "m31 k=%u: %zu distinct for %zu classes", k, distinct, classes);
  free(v);
}

/* 4 + 5: sym_seq (every impl) == reference roller, hashes and positions */
typedef struct { kp7m31_roller r0; const kp7m31_params* p; uint64_t bad, n; } ctx_t;
static void emit(void* vc, size_t start, const uint64_t* h, size_t n)
{
  ctx_t* c = vc;
  for (size_t j = 0; j < n; j++) {
    if (!kp7m31_roll(&c->r0)) { c->bad++; return; }
    if (c->r0.pos != start + j || h[j] != kp7m31_sym(c->p, c->r0.F, c->r0.R)) c->bad++;
    c->n++;
  }
}
static void seq_case(const char* seq, size_t len, int impl, size_t cap, const char* label)
{
  kp7m31_params p;
  kp7m31_init(&p, 31, 4242, 4242 + 1);
  ctx_t c = { { 0 }, &p, 0, 0 };
  kp7m31_start(&c.r0, &p, seq, len);
  uint64_t* buf = malloc(cap * sizeof *buf);
  size_t total = kp7m31_sym_seq(&p, seq, len, impl, buf, cap, emit, &c);
  if (kp7m31_roll(&c.r0)) c.bad++;
  CHECK(!c.bad && total == c.n, "sym_seq impl %d %s cap %zu: %llu bad, %zu vs %llu", impl, label, cap,
        (unsigned long long)c.bad, total, (unsigned long long)c.n);
  free(buf);
}

/* x and lam are independent (round-1 L1 F1): x depends only on seed_x, lam only on seed_lam; equal seeds give
 * unrelated x and lam; OS mode (NULL seeds) gives fresh values; every draw is in range */
static void test_seeds(void)
{
  kp7m31_params a, b, c, d;
  kp7m31_init(&a, 31, 11, 22);
  kp7m31_init(&b, 31, 11, 23);
  kp7m31_init(&c, 31, 12, 22);
  CHECK(a.x.re == b.x.re && a.x.im == b.x.im && (a.lam.re != b.lam.re || a.lam.im != b.lam.im),
        "seeds: changing seed_lam must change lam only");
  CHECK((a.x.re != c.x.re || a.x.im != c.x.im) && a.lam.re == c.lam.re && a.lam.im == c.lam.im,
        "seeds: changing seed_x must change x only (lam derived from x?)");
  kp7m31_init(&d, 31, 5, 5);
  CHECK(d.x.re != d.lam.re || d.x.im != d.lam.im, "seeds: equal seeds gave x == lam");
  kp7m31_init_seeds(&a, 31, NULL, NULL);
  kp7m31_init_seeds(&b, 31, NULL, NULL);
  CHECK(a.x.re != b.x.re && a.lam.re != b.lam.re, "seeds: two OS-random inits agree");
  int bad = 0;
  for (uint64_t s = 0; s < 2000; s++) {
    kp7m31_init(&a, 3, s, s * 3 + 1);
    bad += !(a.x.re < KP7M31_P && a.x.im < KP7M31_P && (a.x.re | a.x.im) && a.lam.re < KP7M31_P && a.lam.im < KP7M31_P);
  }
  CHECK(!bad, "seeds: %d draws out of range", bad);
}

int main(void)
{
  test_seeds();
  uint64_t s = 31337;
  const size_t len = 30000;
  char* seqs[4] = { malloc(len), malloc(len), malloc(len), malloc(len) };
  for (size_t i = 0; i < len; i++) {
    seqs[0][i] = "ACGTN"[sm(&s) % 1000 == 0 ? 4 : sm(&s) % 4];
    seqs[1][i] = "ACGTN"[sm(&s) % 25 == 0 ? 4 : sm(&s) % 4];
    seqs[2][i] = (i / 997) % 3 == 1 && i % 997 < 40 ? 'N' : "ACGT"[sm(&s) % 4];
    seqs[3][i] = "acgtACGTn"[sm(&s) % 9]; /* lowercase and n */
  }
  seqs[1][0] = 'N';
  seqs[0][31] = 'N';
  const char* labels[4] = { "rare N", "frequent N", "N runs", "mixed case" };
  for (int q = 0; q < 4; q++) {
    test_rolling(31, seqs[q], len, labels[q]);
    test_rolling(8, seqs[q], 5000, labels[q]);
    test_rolling(1, seqs[q], 500, labels[q]);
  }
  test_exhaustive(7, 8192);
  test_exhaustive(8, 32896);
  const char* exact = "ACGTACGTACGTACGTACGTACGTACGTACG";
  const size_t caps[11] = { 1, 3, 5, 7, 8, 9, 15, 16, 17, 64, 4096 };
  for (int impl = 0; impl <= 2; impl++) {
    if (!kp7m31_have(impl)) { printf("  impl %d not available on this CPU/build: skipped\n", impl); continue; }
    for (int c = 0; c < 11; c++) {
      for (int q = 0; q < 4; q++) seq_case(seqs[q], len, impl, caps[c], labels[q]);
      seq_case(exact, 31, impl, caps[c], "length k");
      seq_case(exact, 30, impl, caps[c], "shorter than k");
      /* long single run: remainders at every block size */
      seq_case(seqs[2], 997 * 1 + 30, impl, caps[c], "one run");
    }
  }
  for (int q = 0; q < 4; q++) free(seqs[q]);
  printf(fails ? "FAILED: %d check(s)\n" : "all tests passed\n", fails);
  return fails != 0;
}
