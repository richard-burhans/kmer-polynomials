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

/* 4: Karatsuba gives bit-identical products and hashes */
static uint64_t sm(uint64_t* s)
{
  uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static void test_karatsuba(void)
{
  const uint64_t edge[] = { 0, 1, 2, KP7_Q - 2, KP7_Q - 1, (KP7_Q - 1) / 2, 1ULL << 60 };
  const int ne = sizeof edge / sizeof edge[0];
  for (int a = 0; a < ne; a++) for (int b = 0; b < ne; b++) for (int c = 0; c < ne; c++) for (int d = 0; d < ne; d++) {
    kp7_fe x = { edge[a], edge[b] }, y = { edge[c], edge[d] };
    CHECK(fe_eq(kp7_mul(x, y), kp7_mul_karatsuba(x, y)), "karatsuba edge (%d,%d)x(%d,%d)", a, b, c, d);
  }
  uint64_t s = 42;
  for (int n = 0; n < 1000000; n++) {
    kp7_fe x = { sm(&s) % KP7_Q, sm(&s) % KP7_Q }, y = { sm(&s) % KP7_Q, sm(&s) % KP7_Q };
    CHECK(fe_eq(kp7_mul(x, y), kp7_mul_karatsuba(x, y)), "karatsuba random %d", n);
  }
  const size_t len = 100000;
  char* seq = malloc(len);
  for (size_t i = 0; i < len; i++) seq[i] = "ACGTN"[sm(&s) % 1000 == 0 ? 4 : sm(&s) % 4]; /* ~97% of 31-mer windows N-free */
  for (int variant = 0; variant < 2; variant++) {
    kp7_params p0, p1;
    kp7_init(&p0, 31, variant, 99);
    kp7_init(&p1, 31, variant, 99);
    p1.karatsuba = 1;
    kp7_roller r0, r1;
    kp7_start(&r0, &p0, seq, len);
    kp7_start(&r1, &p1, seq, len);
    uint64_t h0, h1, n = 0;
    for (;;) {
      int a = kp7_roll64(&r0, KP7_OUT_SYM, &h0), b = kp7_roll64(&r1, KP7_OUT_SYM, &h1);
      CHECK(a == b, "karatsuba roll lengths differ");
      if (!a || !b) break;
      CHECK(h0 == h1 && r0.pos == r1.pos, "karatsuba sym differs at %zu (variant %d)", r0.pos, variant);
      n++;
    }
    CHECK(n > 90000, "karatsuba roll compared only %llu windows", (unsigned long long)n);
  }
  free(seq);
}

/* 5: kp7_roll_sym gives the same hashes as kp7_roll64(SYM) for every option combination */
static void test_roll_sym(void)
{
  uint64_t s = 7;
  /* the wide product at its bound: parts up to 4q - 1 must match the product of the reduced parts */
  const uint64_t big[] = { 0, KP7_Q - 1, KP7_Q, 2 * KP7_Q - 1, 3 * KP7_Q + 5, 4 * KP7_Q - 1 };
  for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) for (int c = 0; c < 6; c++) for (int d = 0; d < 6; d++) {
    kp7_fe x = { big[a], big[b] }, y = { big[c], big[d] };
    kp7_fe xr = { big[a] % KP7_Q, big[b] % KP7_Q }, yr = { big[c] % KP7_Q, big[d] % KP7_Q };
    CHECK(fe_eq(kp7_mul_wide(x, y), kp7_mul(xr, yr)), "wide product at bound (%d,%d)x(%d,%d)", a, b, c, d);
  }
  for (int n = 0; n < 200000; n++) {
    kp7_fe x = { sm(&s) % (4 * KP7_Q), sm(&s) % (4 * KP7_Q) }, y = { sm(&s) % (4 * KP7_Q), sm(&s) % (4 * KP7_Q) };
    kp7_fe xr = { x.re % KP7_Q, x.im % KP7_Q }, yr = { y.re % KP7_Q, y.im % KP7_Q };
    CHECK(fe_eq(kp7_mul_wide(x, y), kp7_mul(xr, yr)), "wide product random %d", n);
  }
  /* two sequences: rare N (long rolls) and frequent N (many restarts) */
  const size_t len = 60000;
  char* seqs[2] = { malloc(len), malloc(len) };
  const unsigned nrate[2] = { 1000, 25 };
  for (int q = 0; q < 2; q++)
    for (size_t i = 0; i < len; i++) seqs[q][i] = "ACGTN"[sm(&s) % nrate[q] == 0 ? 4 : sm(&s) % 4];
  for (int variant = 0; variant < 2; variant++)
    for (int kar = 0; kar < 2; kar++)
      for (int opt = 0; opt < 8; opt++)
        for (int q = 0; q < 2; q++) {
          kp7_params p0, p1;
          kp7_init(&p0, 31, variant, 1234);
          kp7_init(&p1, 31, variant, 1234);
          p1.karatsuba = kar;
          p1.opt = opt;
          kp7_roller r0, r1;
          kp7_start(&r0, &p0, seqs[q], len);
          kp7_start(&r1, &p1, seqs[q], len);
          uint64_t h0, h1, n = 0, bad = 0;
          for (;;) {
            int a = kp7_roll64(&r0, KP7_OUT_SYM, &h0), b = kp7_roll_sym(&r1, &h1);
            if (a != b) { bad++; break; }
            if (!a) break;
            if (h0 != h1 || r0.pos != r1.pos) bad++;
            n++;
          }
          CHECK(!bad, "roll_sym differs: variant %d karatsuba %d opt %d seq %d (%llu bad)", variant, kar, opt, q,
                (unsigned long long)bad);
          CHECK(n > (q ? 10000u : 55000u), "roll_sym compared only %llu windows (seq %d)", (unsigned long long)n, q);
        }
  free(seqs[0]); free(seqs[1]);
}

/* 6: kp7_sym_batch: same hashes AND positions as kp7_roll64(SYM), for any buffer size and N layout */
static void batch_case(const char* seq, size_t len, int variant, int kar, size_t bufsz, const char* label)
{
  kp7_params p0, p1;
  kp7_init(&p0, 31, variant, 4321);
  kp7_init(&p1, 31, variant, 4321);
  p1.karatsuba = kar;
  kp7_roller r0, r1;
  kp7_start(&r0, &p0, seq, len);
  kp7_start(&r1, &p1, seq, len);
  uint64_t* buf = malloc(bufsz * sizeof *buf);
  uint64_t h0, ref = 0, got = 0, bad = 0;
  size_t i = 0, n;
  while ((n = kp7_sym_batch(&r1, buf, bufsz)) > 0) {
    size_t first = r1.pos - (n - 1);
    for (size_t j = 0; j < n; j++, i++) {
      if (!kp7_roll64(&r0, KP7_OUT_SYM, &h0)) { bad++; break; }
      ref++; got++;
      if (h0 != buf[j] || r0.pos != first + j) bad++;
    }
  }
  if (kp7_roll64(&r0, KP7_OUT_SYM, &h0)) bad++; /* the batch stopped before the reference did */
  CHECK(!bad, "sym_batch differs: %s variant %d karatsuba %d buf %zu (%llu bad of %llu)", label, variant, kar, bufsz,
        (unsigned long long)bad, (unsigned long long)got);
  free(buf);
}
static void test_sym_batch(void)
{
  uint64_t s = 99;
  const size_t len = 40000;
  char* rare = malloc(len); char* freq = malloc(len); char* runs = malloc(len);
  for (size_t i = 0; i < len; i++) {
    rare[i] = "ACGTN"[sm(&s) % 1000 == 0 ? 4 : sm(&s) % 4];
    freq[i] = "ACGTN"[sm(&s) % 25 == 0 ? 4 : sm(&s) % 4];
    runs[i] = (i / 997) % 3 == 1 && i % 997 < 40 ? 'N' : "ACGT"[sm(&s) % 4]; /* N runs longer than k */
  }
  freq[0] = 'N';  /* N at the very start */
  rare[31] = 'N'; /* N right after the first k-mer: the second call's first incoming base */
  const char* exact = "ACGTACGTACGTACGTACGTACGTACGTACG"; /* exactly k = 31 */
  const size_t bufs[4] = { 1, 2, 7, 4096 };
  for (int variant = 0; variant < 2; variant++)
    for (int kar = 0; kar < 2; kar++)
      for (int b = 0; b < 4; b++) {
        batch_case(rare, len, variant, kar, bufs[b], "rare N");
        batch_case(freq, len, variant, kar, bufs[b], "frequent N");
        batch_case(runs, len, variant, kar, bufs[b], "N runs");
        batch_case(exact, 31, variant, kar, bufs[b], "length k");
        batch_case(exact, 30, variant, kar, bufs[b], "shorter than k");
      }
  free(rare); free(freq); free(runs);
}

/* 7: kp7_sym_seq (interleaved lanes): same hashes and positions as kp7_roll64(SYM) */
typedef struct { const char* seq; size_t len; kp7_roller r0; uint64_t bad, n; } seq_ctx;
static void seq_emit(void* vctx, size_t start, const uint64_t* h, size_t n)
{
  seq_ctx* c = vctx;
  uint64_t h0;
  for (size_t j = 0; j < n; j++) {
    if (!kp7_roll64(&c->r0, KP7_OUT_SYM, &h0)) { c->bad++; return; }
    if (h0 != h[j] || c->r0.pos != start + j) c->bad++;
    c->n++;
  }
}
static void seq_case(const char* seq, size_t len, int variant, int lanes, size_t cap, const char* label)
{
  kp7_params p;
  kp7_init(&p, 31, variant, 777);
  seq_ctx c = { seq, len, { 0 }, 0, 0 };
  kp7_start(&c.r0, &p, seq, len);
  uint64_t* buf = malloc(cap * sizeof *buf);
  size_t total = kp7_sym_seq(&p, seq, len, lanes, buf, cap, seq_emit, &c);
  uint64_t h0;
  if (kp7_roll64(&c.r0, KP7_OUT_SYM, &h0)) c.bad++; /* the reference had more k-mers */
  CHECK(!c.bad && total == c.n, "sym_seq differs: %s variant %d lanes %d cap %zu (%llu bad, %zu vs %llu)", label, variant,
        lanes, cap, (unsigned long long)c.bad, total, (unsigned long long)c.n);
  free(buf);
}
static void test_sym_seq(void)
{
  uint64_t s = 1234;
  const size_t len = 30000;
  char* rare = malloc(len); char* freq = malloc(len); char* runs = malloc(len);
  for (size_t i = 0; i < len; i++) {
    rare[i] = "ACGTN"[sm(&s) % 1000 == 0 ? 4 : sm(&s) % 4];
    freq[i] = "ACGTN"[sm(&s) % 25 == 0 ? 4 : sm(&s) % 4];
    runs[i] = (i / 997) % 3 == 1 && i % 997 < 40 ? 'N' : "ACGT"[sm(&s) % 4];
  }
  freq[0] = 'N';
  rare[31] = 'N';
  const char* exact = "ACGTACGTACGTACGTACGTACGTACGTACG";
  const size_t caps[5] = { 1, 3, 5, 64, 4096 };
  for (int variant = 0; variant < 2; variant++)
    for (int lanes = 1; lanes <= (variant ? 2 : 4); lanes++)
      for (int c = 0; c < 5; c++) {
        seq_case(rare, len, variant, lanes, caps[c], "rare N");
        seq_case(freq, len, variant, lanes, caps[c], "frequent N");
        seq_case(runs, len, variant, lanes, caps[c], "N runs");
        seq_case(exact, 31, variant, lanes, caps[c], "length k");
        seq_case(exact, 30, variant, lanes, caps[c], "shorter than k");
      }
  free(rare); free(freq); free(runs);
}

/* 8: the hybrid. gate == reference gate everywhere; gate symmetric; exhaustive distinctness; sampling exact */
typedef struct { kp7_roller r0; const kp7_params* p; uint64_t bad, n; } gate_ctx;
static void gate_emit(void* vctx, size_t start, const uint64_t* h, size_t n)
{
  gate_ctx* c = vctx;
  uint64_t h0;
  for (size_t j = 0; j < n; j++) {
    if (!kp7_roll64(&c->r0, KP7_OUT_FWD, &h0)) { c->bad++; return; }
    if (c->r0.pos != start + j || h[j] != kp7_gate(c->p, c->r0.F, c->r0.R)) c->bad++;
    c->n++;
  }
}
typedef struct { const kp7_params* p; const char* seq; size_t len; uint64_t thr; size_t* pos; uint64_t* h; size_t n, cap_; } samp_ctx;
static void samp_emit(void* vctx, const size_t* pos, const uint64_t* h, size_t n)
{
  samp_ctx* c = vctx;
  for (size_t j = 0; j < n; j++) { c->pos[c->n] = pos[j]; c->h[c->n] = h[j]; c->n++; }
}
static void test_hybrid(void)
{
  uint64_t s = 555;
  const size_t len = 30000;
  char* seqs[3] = { malloc(len), malloc(len), malloc(len) };
  for (size_t i = 0; i < len; i++) {
    seqs[0][i] = "ACGTN"[sm(&s) % 1000 == 0 ? 4 : sm(&s) % 4];
    seqs[1][i] = "ACGTN"[sm(&s) % 25 == 0 ? 4 : sm(&s) % 4];
    seqs[2][i] = (i / 997) % 3 == 1 && i % 997 < 40 ? 'N' : "ACGT"[sm(&s) % 4];
  }
  seqs[1][0] = 'N';
  kp7_params p;
  kp7_init(&p, 31, KP7_XREAL, 2468);
  const size_t caps[4] = { 1, 3, 64, 4096 };
  /* (a) gate_seq equals the reference gate, for both lane counts and several block sizes */
  for (int q = 0; q < 3; q++)
    for (int lanes = 1; lanes <= 2; lanes++)
      for (int c = 0; c < 4; c++) {
        gate_ctx g = { { 0 }, &p, 0, 0 };
        kp7_start(&g.r0, &p, seqs[q], len);
        uint64_t* buf = malloc(caps[c] * sizeof *buf);
        size_t total = kp7_gate_seq(&p, seqs[q], len, lanes, buf, caps[c], gate_emit, &g);
        uint64_t h0;
        if (kp7_roll64(&g.r0, KP7_OUT_FWD, &h0)) g.bad++;
        CHECK(!g.bad && total == g.n && total > 0, "gate_seq differs: seq %d lanes %d cap %zu (%llu bad)", q, lanes, caps[c],
              (unsigned long long)g.bad);
        free(buf);
      }
  /* (b) symmetry: gate(s) == gate(rc s) on every window of seq 0 */
  {
    char w[32], wr[32];
    uint64_t bad = 0;
    for (size_t i = 0; i + 31 <= len; i++) {
      memcpy(w, seqs[0] + i, 31); w[31] = 0;
      kp7_fe F, R, F2, R2;
      if (!kp7_direct(&p, w, &F, &R)) continue;
      rc(w, wr, 31);
      kp7_direct(&p, wr, &F2, &R2);
      if (kp7_gate(&p, F, R) != kp7_gate(&p, F2, R2)) bad++;
    }
    CHECK(!bad, "gate not strand-symmetric on %llu windows", (unsigned long long)bad);
  }
  /* (c) exhaustive: one gate value per double-stranded class at k = 7 and 8 (heuristic quality, 64-bit gate) */
  for (unsigned k = 7; k <= 8; k++) {
    kp7_params pk;
    kp7_init(&pk, k, KP7_XREAL, 13579);
    size_t nk = (size_t)1 << (2 * k), classes = k == 7 ? 8192 : 32896, distinct = 0;
    uint64_t* v = malloc(nk * sizeof *v);
    char w[16];
    for (size_t code = 0; code < nk; code++) {
      for (unsigned t = 0; t < k; t++) w[t] = "ACGT"[(code >> (2 * t)) & 3];
      w[k] = 0;
      kp7_fe F, R;
      kp7_direct(&pk, w, &F, &R);
      v[code] = kp7_gate(&pk, F, R);
    }
    /* sort and count distinct */
    for (size_t gap = nk / 2; gap > 0; gap /= 2) /* shell sort, small n */
      for (size_t i = gap; i < nk; i++) for (size_t j = i; j >= gap && v[j - gap] > v[j]; j -= gap) { uint64_t t = v[j]; v[j] = v[j - gap]; v[j - gap] = t; }
    for (size_t i = 0; i < nk; i++) if (i == 0 || v[i] != v[i - 1]) distinct++;
    printf("  gate k=%u: distinct %zu (classes %zu)\n", k, distinct, classes);
    CHECK(distinct == classes, "gate k=%u: %zu distinct values for %zu classes", k, distinct, classes);
    free(v);
  }
  /* (d) sampling: kept set == {gate < thr} with provable hashes equal to kp7_roll64(SYM) */
  const uint64_t thrs[3] = { 0, 1ULL << 60, ~0ULL };
  for (int q = 0; q < 3; q++)
    for (int lanes = 1; lanes <= 2; lanes++)
      for (int t = 0; t < 3; t++)
        for (int c = 0; c < 4; c++) {
          samp_ctx sc = { &p, seqs[q], len, thrs[t], malloc(len * sizeof(size_t)), malloc(len * sizeof(uint64_t)), 0, 0 };
          size_t* pb = malloc(caps[c] * sizeof *pb);
          uint64_t* hb = malloc(caps[c] * sizeof *hb);
          size_t kept = kp7_sym_sample(&p, seqs[q], len, lanes, thrs[t], pb, hb, caps[c], samp_emit, &sc);
          /* reference */
          kp7_roller r0;
          kp7_start(&r0, &p, seqs[q], len);
          uint64_t hs, bad = 0;
          size_t j = 0, expect = 0;
          while (kp7_roll64(&r0, KP7_OUT_SYM, &hs)) {
            if (kp7_gate(&p, r0.F, r0.R) >= thrs[t]) continue;
            expect++;
            if (j >= sc.n || sc.pos[j] != r0.pos || sc.h[j] != hs) bad++;
            j++;
          }
          CHECK(!bad && kept == expect && sc.n == expect, "sym_sample differs: seq %d lanes %d thr %d cap %zu (%llu bad, kept %zu want %zu)",
                q, lanes, t, caps[c], (unsigned long long)bad, kept, expect);
          if (t == 1 && c == 3 && lanes == 1)
            printf("  sample seq %d at 1/16: kept %zu of %llu (%.4f)\n", q, kept, (unsigned long long)(r0.pos + 1), (double)kept / (r0.pos + 1));
          free(sc.pos); free(sc.h); free(pb); free(hb);
        }
  free(seqs[0]); free(seqs[1]); free(seqs[2]);
}

int main(void)
{
  test_field();
  test_karatsuba();
  test_roll_sym();
  test_sym_batch();
  test_sym_seq();
  test_hybrid();
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
