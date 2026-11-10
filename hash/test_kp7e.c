/* Tests for kp7e (Eisenstein KP-7 in F_q, q = 2^61 - 1). Exit 0 iff all pass.
 *  1. field: w^3 = 1, w != 1, w^2 + w + 1 = 0, the four base values distinct, inverses
 *  2. the complex-embedding metric: squared distances 1 (transitions), 3 (A-C, G-T), 4 (complements)
 *  3. rolling == direct at every window, for N layouts and k = 31, 8, 1
 *  4. sym(s) == sym(rc s) at every window
 *  5. exhaustive k = 7, 8: F distinct for all 4^k k-mers, S distinct per double-stranded class (8192, 32896)
 *  6. kp7e_sym_seq (1 and 2 lanes) == the reference, hashes AND positions, block sizes 1/3/5/64/4096
 *  7. table evaluation == direct for k = 1..32; gate == kp7pack_gate; samplers == {gate < thr} with reference S */
#include "kp7e.h"
#include "kp7pack.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fails = 0;
#define CHECK(cond, ...) do { if (!(cond)) { fails++; if (fails < 12) { printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } } while (0)

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
    out[t] = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N';
  }
  out[k] = 0;
}

static void test_field(void)
{
  kp7e_params p;
  kp7e_init(&p, 31, 1);
  uint64_t w = p.w, w2 = kp7e_mul(w, w), w3 = kp7e_mul(w2, w);
  CHECK(w != 1 && w3 == 1, "w is not a primitive cube root of unity");
  CHECK(kp7e_add(kp7e_add(w2, w), 1) == 0, "w^2 + w + 1 != 0");
  for (int a = 0; a < 4; a++) for (int b = a + 1; b < 4; b++) CHECK(p.val[a] != p.val[b], "base values %d, %d equal", a, b);
  for (int c = 0; c < 4; c++) CHECK(kp7e_add(p.val[c], p.val[3 - c]) == 0, "complement of %d is not its negation", c);
  CHECK(kp7e_mul(p.x, p.xinv) == 1, "x * x^-1 != 1");
}

static void test_metric(void)
{
  const double complex w = cexp(2 * M_PI * I / 3);
  const double complex v[4] = { 1, w, -w, -1 }; /* A C G T */
  double d2[4][4];
  for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) d2[a][b] = cabs(v[a] - v[b]) * cabs(v[a] - v[b]);
  /* transitions A-G, C-T; transversions A-C, G-T; complements A-T, C-G */
  CHECK(fabs(d2[0][2] - 1) < 1e-12 && fabs(d2[1][3] - 1) < 1e-12, "transitions not at squared distance 1");
  CHECK(fabs(d2[0][1] - 3) < 1e-12 && fabs(d2[2][3] - 3) < 1e-12, "A-C/G-T not at squared distance 3");
  CHECK(fabs(d2[0][3] - 4) < 1e-12 && fabs(d2[1][2] - 4) < 1e-12, "complements not at squared distance 4");
  printf("  metric (complex embedding): transitions %.0f, A-C/G-T %.0f, complements %.0f\n", d2[0][2], d2[0][1], d2[0][3]);
}

static char* mkseq(size_t len, unsigned nrate, uint64_t* s)
{
  char* q = malloc(len);
  for (size_t i = 0; i < len; i++) q[i] = "ACGTN"[sm(s) % nrate == 0 ? 4 : sm(s) % 4];
  return q;
}

static void test_rolling(unsigned k, size_t len, uint64_t seed)
{
  kp7e_params p;
  kp7e_init(&p, k, seed);
  uint64_t s = seed * 7 + 1;
  char* seq = mkseq(len, 50, &s);
  kp7e_roller r;
  kp7e_start(&r, &p, seq, len);
  size_t n = 0, expect = 0;
  for (size_t i = 0; i + k <= len; i++) { int ok = 1; for (unsigned t = 0; t < k; t++) ok &= seq[i + t] != 'N'; expect += ok; }
  char w[40], wr[40];
  while (kp7e_roll(&r)) {
    uint64_t F, R, F2, R2;
    kp7e_direct(&p, seq + r.pos, &F, &R);
    CHECK(F == r.F && R == r.R, "k=%u rolling != direct at %zu", k, r.pos);
    memcpy(w, seq + r.pos, k); w[k] = 0; rc(w, wr, k);
    kp7e_direct(&p, wr, &F2, &R2);
    CHECK(F2 == R && R2 == F, "k=%u F(rc s) != R(s) at %zu", k, r.pos);
    CHECK(kp7e_sym(&p, F, R) == kp7e_sym(&p, F2, R2), "k=%u sym not strand-symmetric at %zu", k, r.pos);
    n++;
  }
  CHECK(n == expect && n > 0, "k=%u rolled %zu windows, expected %zu", k, n, expect);
  free(seq);
}

static int cmp64(const void* a, const void* b) { uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b; return x < y ? -1 : x > y; }
static void test_exhaustive(unsigned k, size_t classes)
{
  kp7e_params p;
  kp7e_init(&p, k, 99 + k);
  size_t nk = (size_t)1 << (2 * k), dF = 0, dS = 0;
  uint64_t* F = malloc(nk * 8);
  uint64_t* S = malloc(nk * 8);
  char w[16];
  for (size_t code = 0; code < nk; code++) {
    for (unsigned t = 0; t < k; t++) w[t] = "ACGT"[(code >> (2 * t)) & 3];
    w[k] = 0;
    uint64_t f, r;
    kp7e_direct(&p, w, &f, &r);
    F[code] = f; S[code] = kp7e_sym(&p, f, r);
  }
  qsort(F, nk, 8, cmp64); qsort(S, nk, 8, cmp64);
  for (size_t i = 0; i < nk; i++) { dF += i == 0 || F[i] != F[i - 1]; dS += i == 0 || S[i] != S[i - 1]; }
  printf("  exhaustive k=%u: F distinct %zu of %zu, sym distinct %zu (classes %zu)\n", k, dF, nk, dS, classes);
  CHECK(dF == nk, "k=%u: F collides", k);
  CHECK(dS == classes, "k=%u: %zu sym values for %zu classes", k, dS, classes);
  free(F); free(S);
}

/* sym_seq against the reference roller: hashes and positions */
typedef struct { kp7e_roller r0; const kp7e_params* p; uint64_t bad, n; } seq_ctx;
static void seq_emit(void* vc, size_t start, const uint64_t* h, size_t n)
{
  seq_ctx* c = vc;
  for (size_t j = 0; j < n; j++) {
    if (!kp7e_roll(&c->r0)) { c->bad++; return; }
    if (c->r0.pos != start + j || h[j] != kp7e_sym(c->p, c->r0.F, c->r0.R)) c->bad++;
    c->n++;
  }
}
static void seq_case(const char* seq, size_t len, int lanes, size_t cap, const char* label)
{
  kp7e_params p;
  kp7e_init(&p, 31, 777);
  seq_ctx c = { { 0 }, &p, 0, 0 };
  kp7e_start(&c.r0, &p, seq, len);
  uint64_t* buf = malloc(cap * 8);
  size_t total = kp7e_sym_seq(&p, seq, len, lanes, buf, cap, seq_emit, &c);
  if (kp7e_roll(&c.r0)) c.bad++;
  CHECK(!c.bad && total == c.n, "sym_seq differs: %s lanes %d cap %zu (%llu bad)", label, lanes, cap, (unsigned long long)c.bad);
  free(buf);
}
static void test_sym_seq(void)
{
  uint64_t s = 4242;
  const size_t len = 30000;
  char* rare = mkseq(len, 1000, &s);
  char* freq = mkseq(len, 25, &s);
  char* runs = malloc(len);
  for (size_t i = 0; i < len; i++) runs[i] = (i / 997) % 3 == 1 && i % 997 < 40 ? 'N' : "ACGT"[sm(&s) % 4];
  freq[0] = 'N'; rare[31] = 'N';
  const char* exact = "ACGTACGTACGTACGTACGTACGTACGTACG";
  const size_t caps[5] = { 1, 3, 5, 64, 4096 };
  for (int lanes = 1; lanes <= 2; lanes++)
    for (int c = 0; c < 5; c++) {
      seq_case(rare, len, lanes, caps[c], "rare N");
      seq_case(freq, len, lanes, caps[c], "frequent N");
      seq_case(runs, len, lanes, caps[c], "N runs");
      seq_case(exact, 31, lanes, caps[c], "length k");
      seq_case(exact, 30, lanes, caps[c], "shorter than k");
    }
  free(rare); free(freq); free(runs);
}

/* table evaluation and samplers */
typedef struct { size_t* pos; uint64_t* h; size_t n; } samp_ctx;
static void samp_emit(void* vc, const size_t* pos, const uint64_t* h, size_t n)
{
  samp_ctx* c = vc;
  for (size_t j = 0; j < n; j++) { c->pos[c->n] = pos[j]; c->h[c->n] = h[j]; c->n++; }
}
static void test_tables(void)
{
  uint64_t s = 31337;
  /* table == direct for every k, random k-mers including all-A and all-T */
  for (unsigned k = 1; k <= 32; k++) {
    kp7e_params p;
    kp7e_init(&p, k, 500 + k);
    static kp7e_tab t;
    kp7e_tab_init(&t, &p);
    char w[40];
    for (int n = 0; n < 2000; n++) {
      for (unsigned i = 0; i < k; i++) w[i] = n == 0 ? 'A' : n == 1 ? 'T' : "ACGT"[sm(&s) % 4];
      uint64_t f, r, F, R;
      kp7e_codes(w, k, &f, &r);
      kp7e_direct(&p, w, &F, &R);
      CHECK(kp7e_tab_sym(&t, f, r) == kp7e_sym(&p, F, R), "table != direct at k=%u", k);
      CHECK(kp7e_gate(f, r) == kp7pack_gate(f, r), "gate differs from kp7pack_gate");
    }
  }
  /* samplers: kept set == {gate < thr}, hashes == reference S, positions in order */
  const size_t len = 30000;
  char* seqs[3] = { mkseq(len, 1000, &s), mkseq(len, 25, &s), malloc(len) };
  for (size_t i = 0; i < len; i++) seqs[2][i] = (i / 997) % 3 == 1 && i % 997 < 40 ? 'N' : "ACGT"[sm(&s) % 4];
  seqs[1][0] = 'N';
  kp7e_params p;
  kp7e_init(&p, 31, 2468);
  static kp7e_tab t;
  kp7e_tab_init(&t, &p);
  const uint64_t thrs[4] = { 0, 1ULL << 60, 1ULL << 62, ~0ULL };
  const size_t caps[4] = { 1, 3, 64, 4096 };
  for (int q = 0; q < 3; q++)
    for (int th = 0; th < 4; th++)
      for (int c = 0; c < 4; c++)
        for (int which = 0; which < 2; which++) {
          samp_ctx sc = { malloc(len * sizeof(size_t)), malloc(len * 8), 0 };
          size_t* pb = malloc(caps[c] * sizeof(size_t));
          uint64_t* hb = malloc(caps[c] * 8);
          size_t kept = which ? kp7e_sample_batch(&t, seqs[q], len, thrs[th], pb, hb, caps[c], samp_emit, &sc)
                              : kp7e_sample_b(&t, seqs[q], len, thrs[th], pb, hb, caps[c], samp_emit, &sc);
          kp7e_roller r0;
          kp7e_start(&r0, &p, seqs[q], len);
          size_t j = 0, expect = 0;
          uint64_t bad = 0;
          while (kp7e_roll(&r0)) {
            uint64_t f, r;
            kp7e_codes(seqs[q] + r0.pos, 31, &f, &r);
            if (kp7e_gate(f, r) >= thrs[th]) continue;
            expect++;
            if (j >= sc.n || sc.pos[j] != r0.pos || sc.h[j] != kp7e_sym(&p, r0.F, r0.R)) bad++;
            j++;
          }
          CHECK(!bad && kept == expect && sc.n == expect, "%s seq %d thr %d cap %zu: %llu bad, kept %zu want %zu",
                which ? "sample_batch" : "sample_b", q, th, caps[c], (unsigned long long)bad, kept, expect);
          free(sc.pos); free(sc.h); free(pb); free(hb);
        }
  free(seqs[0]); free(seqs[1]); free(seqs[2]);
}

int main(void)
{
  test_field();
  test_metric();
  test_rolling(31, 200000, 7);
  test_rolling(8, 50000, 9);
  test_rolling(1, 2000, 11);
  test_exhaustive(7, 8192);
  test_exhaustive(8, 32896);
  test_sym_seq();
  test_tables();
  printf(fails ? "FAILED: %d check(s)\n" : "all tests passed\n", fails);
  return fails != 0;
}
