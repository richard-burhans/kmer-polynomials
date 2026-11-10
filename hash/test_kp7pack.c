/* Tests for kp7pack (hybrid v2). Exit 0 iff all pass. */
#include "kp7pack.h"
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
    out[t] = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N';
  }
  out[k] = 0;
}

/* sequences: rare N, frequent N (leading N), N runs longer than k, length k, shorter than k */
#define NSEQ 5
static char* SEQ[NSEQ];
static size_t LEN[NSEQ];
static const char* LABEL[NSEQ] = { "rare N", "frequent N", "N runs", "length k", "shorter than k" };
static void make_seqs(void)
{
  uint64_t s = 31337;
  const size_t len = 30000;
  for (int q = 0; q < 3; q++) { SEQ[q] = malloc(len); LEN[q] = len; }
  for (size_t i = 0; i < len; i++) {
    SEQ[0][i] = "ACGTN"[sm(&s) % 1000 == 0 ? 4 : sm(&s) % 4];
    SEQ[1][i] = "ACGTN"[sm(&s) % 25 == 0 ? 4 : sm(&s) % 4];
    SEQ[2][i] = (i / 997) % 3 == 1 && i % 997 < 40 ? 'N' : "acgt"[sm(&s) % 4]; /* lowercase too */
  }
  SEQ[1][0] = 'N';
  SEQ[0][31] = 'N';
  SEQ[3] = strdup("ACGTACGTACGTACGTACGTACGTACGTACG"); LEN[3] = 31;
  SEQ[4] = strdup("ACGTACGTACGTACGTACGTACGTACGTAC"); LEN[4] = 30;
}

/* 1: gate_seq equals the reference gate from directly built codes, at every window */
typedef struct { const char* seq; unsigned k; size_t next; uint64_t bad, n; } gctx;
static void gemit(void* v, size_t start, const uint64_t* h, size_t n)
{
  gctx* c = v;
  for (size_t j = 0; j < n; j++) {
    size_t pos = start + j;
    uint64_t f, r;
    /* every window between the previous emitted position and this one must contain N */
    for (size_t q = c->next; q < pos; q++) if (kp7pack_codes(c->seq + q, c->k, &f, &r)) c->bad++;
    if (!kp7pack_codes(c->seq + pos, c->k, &f, &r) || kp7pack_gate(f, r) != h[j]) c->bad++;
    c->next = pos + 1;
    c->n++;
  }
}
static void test_gate_seq(void)
{
  const size_t caps[4] = { 1, 3, 64, 4096 };
  for (int q = 0; q < NSEQ; q++)
    for (int c = 0; c < 4; c++) {
      gctx g = { SEQ[q], 31, 0, 0, 0 };
      uint64_t* buf = malloc(caps[c] * sizeof *buf);
      size_t total = kp7pack_gate_seq(31, SEQ[q], LEN[q], buf, caps[c], gemit, &g);
      uint64_t f, r;
      for (size_t p = g.next; p + 31 <= LEN[q]; p++) if (kp7pack_codes(SEQ[q] + p, 31, &f, &r)) g.bad++; /* missed tail */
      CHECK(!g.bad && total == g.n, "gate_seq differs: %s cap %zu (%llu bad)", LABEL[q], caps[c], (unsigned long long)g.bad);
      free(buf);
    }
}

/* 2: gate(s) == gate(rc s), and codes swap under rc */
static void test_symmetry(void)
{
  uint64_t bad = 0, n = 0;
  char w[33], wr[33];
  for (unsigned k = 1; k <= 32; k += (k < 8 ? 1 : 7)) {
    for (size_t i = 0; i + k <= LEN[0]; i += 7) {
      memcpy(w, SEQ[0] + i, k); w[k] = 0;
      uint64_t f, r, f2, r2;
      if (!kp7pack_codes(w, k, &f, &r)) continue;
      rc(w, wr, k);
      kp7pack_codes(wr, k, &f2, &r2);
      if (f2 != r || r2 != f || kp7pack_gate(f, r) != kp7pack_gate(f2, r2)) bad++;
      n++;
    }
  }
  CHECK(!bad && n > 10000, "gate symmetry: %llu bad of %llu", (unsigned long long)bad, (unsigned long long)n);
  /* the seeded gate: symmetric for any seed, and seed 0 is the plain gate */
  uint64_t sbad = 0, ss = 77;
  for (int rep = 0; rep < 20000; rep++) {
    for (unsigned t = 0; t < 31; t++) w[t] = "ACGT"[sm(&ss) % 4];
    w[31] = 0;
    uint64_t f, r, f2, r2, seed = sm(&ss);
    kp7pack_codes(w, 31, &f, &r);
    rc(w, wr, 31);
    kp7pack_codes(wr, 31, &f2, &r2);
    if (kp7pack_gate_seeded(f, r, seed) != kp7pack_gate_seeded(f2, r2, seed) || kp7pack_gate_seeded(f, r, 0) != kp7pack_gate(f, r)) sbad++;
  }
  CHECK(!sbad, "seeded gate: %llu bad", (unsigned long long)sbad);
}

/* 3: exhaustive, one gate value per double-stranded class at k = 7, 8 */
static int cmpu64(const void* a, const void* b) { uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b; return x < y ? -1 : x > y; }
static void test_exhaustive(void)
{
  for (unsigned k = 7; k <= 8; k++) {
    size_t nk = (size_t)1 << (2 * k), classes = k == 7 ? 8192 : 32896, distinct = 0;
    uint64_t* v = malloc(nk * sizeof *v);
    char w[16];
    for (size_t code = 0; code < nk; code++) {
      for (unsigned t = 0; t < k; t++) w[t] = "ACGT"[(code >> (2 * t)) & 3];
      w[k] = 0;
      uint64_t f, r;
      kp7pack_codes(w, k, &f, &r);
      v[code] = kp7pack_gate(f, r);
    }
    qsort(v, nk, sizeof *v, cmpu64);
    for (size_t i = 0; i < nk; i++) if (i == 0 || v[i] != v[i - 1]) distinct++;
    printf("  pack gate k=%u: distinct %zu (classes %zu)\n", k, distinct, classes);
    CHECK(distinct == classes, "pack gate k=%u: %zu distinct for %zu classes", k, distinct, classes);
    free(v);
  }
}

/* 4: sampling: kept == {gate < thr}, provable hashes and positions equal kp7_roll64(SYM) */
typedef struct { size_t* pos; uint64_t* h; size_t n; } sctx;
static void semit(void* v, const size_t* pos, const uint64_t* h, size_t n)
{
  sctx* c = v;
  for (size_t j = 0; j < n; j++) { c->pos[c->n] = pos[j]; c->h[c->n] = h[j]; c->n++; }
}
static void test_sample(void)
{
  const uint64_t thrs[3] = { 0, 1ULL << 60, ~0ULL };
  const size_t caps[4] = { 1, 3, 64, 4096 };
  static kp7pack_tab tab;
  for (int mode = 0; mode < 3; mode++)
  for (int variant = 0; variant < 2; variant++)
    for (int kar = 0; kar < 2; kar++) {
      kp7_params p;
      kp7_init(&p, 31, variant, 8642);
      p.karatsuba = kar;
      kp7pack_tab_init(&tab, &p);
      for (int q = 0; q < NSEQ; q++)
        for (int t = 0; t < 3; t++)
          for (int c = 0; c < 4; c++) {
            sctx sc = { malloc((LEN[q] + 1) * sizeof(size_t)), malloc((LEN[q] + 1) * sizeof(uint64_t)), 0 };
            size_t* pb = malloc(caps[c] * sizeof *pb);
            uint64_t* hb = malloc(caps[c] * sizeof *hb);
            size_t kept = kp7pack_sample_eval(&p, &tab, mode, SEQ[q], LEN[q], thrs[t], pb, hb, caps[c], semit, &sc);
            kp7_roller r0;
            kp7_start(&r0, &p, SEQ[q], LEN[q]);
            uint64_t hs, bad = 0;
            size_t j = 0, expect = 0;
            while (kp7_roll64(&r0, KP7_OUT_SYM, &hs)) {
              uint64_t f, r;
              if (!kp7pack_codes(SEQ[q] + r0.pos, 31, &f, &r)) { bad++; continue; }
              if (kp7pack_gate(f, r) >= thrs[t]) continue;
              expect++;
              if (j >= sc.n || sc.pos[j] != r0.pos || sc.h[j] != hs) bad++;
              j++;
            }
            CHECK(!bad && kept == expect && sc.n == expect,
                  "pack sample differs: mode %d variant %d kar %d %s thr %d cap %zu (%llu bad, kept %zu want %zu)", mode, variant, kar,
                  LABEL[q], t, caps[c], (unsigned long long)bad, kept, expect);
            if (mode == 0 && variant == 0 && kar == 0 && t == 1 && c == 3 && q == 0)
              printf("  pack sample at 1/16 (%s): kept %zu (%.4f of k-mers)\n", LABEL[q], kept, (double)kept / (r0.pos + 1));
            free(sc.pos); free(sc.h); free(pb); free(hb);
          }
    }
}

/* 5: table evaluation equals kp7_direct for every k (1..32), both variants */
static void test_eval(void)
{
  static kp7pack_tab tab;
  uint64_t s = 4242, bad = 0, n = 0;
  char w[33];
  for (int variant = 0; variant < 2; variant++)
    for (unsigned k = 1; k <= 32; k++) {
      kp7_params p;
      kp7_init(&p, k, variant, 1000 + k);
      kp7pack_tab_init(&tab, &p);
      for (int rep = 0; rep < 300; rep++) {
        for (unsigned t = 0; t < k; t++) w[t] = "ACGT"[sm(&s) % 4];
        w[k] = 0;
        uint64_t f, r;
        kp7pack_codes(w, k, &f, &r);
        kp7_fe F0, R0, F1, R1, F2, R2;
        kp7_direct(&p, w, &F0, &R0);
        kp7pack_eval(&tab, KP7PACK_EVAL_POS, f, r, &F1, &R1);
        kp7pack_eval(&tab, KP7PACK_EVAL_CHUNK, f, r, &F2, &R2);
        if (F0.re != F1.re || F0.im != F1.im || R0.re != R1.re || R0.im != R1.im) bad++;
        if (F0.re != F2.re || F0.im != F2.im || R0.re != R2.re || R0.im != R2.im) bad++;
        n++;
      }
    }
  CHECK(!bad && n == 2 * 32 * 300, "table evaluation differs from kp7_direct: %llu bad of %llu", (unsigned long long)bad,
        (unsigned long long)n);
}

/* 6: sym_all (provable hash for every k-mer from packed codes) equals kp7_roll64(SYM), hashes and positions */
typedef struct { kp7_roller r0; uint64_t bad, n; } actx;
static void aemit(void* v, size_t start, const uint64_t* h, size_t n)
{
  actx* c = v;
  uint64_t h0;
  for (size_t j = 0; j < n; j++) {
    if (!kp7_roll64(&c->r0, KP7_OUT_SYM, &h0)) { c->bad++; return; }
    if (c->r0.pos != start + j || h0 != h[j]) c->bad++;
    c->n++;
  }
}
static void test_sym_all(void)
{
  static kp7pack_tab tab;
  const size_t caps[4] = { 1, 3, 64, 4096 };
  for (int variant = 0; variant < 2; variant++) {
    kp7_params p;
    kp7_init(&p, 31, variant, 97531);
    kp7pack_tab_init(&tab, &p);
    for (int q = 0; q < NSEQ; q++)
      for (int c = 0; c < 4; c++) {
        actx a;
        a.bad = a.n = 0;
        kp7_start(&a.r0, &p, SEQ[q], LEN[q]);
        uint64_t* buf = malloc(caps[c] * sizeof *buf);
        size_t total = kp7pack_sym_all(&p, &tab, SEQ[q], LEN[q], buf, caps[c], aemit, &a);
        uint64_t h0;
        if (kp7_roll64(&a.r0, KP7_OUT_SYM, &h0)) a.bad++;
        CHECK(!a.bad && total == a.n, "sym_all differs: variant %d %s cap %zu (%llu bad)", variant, LABEL[q], caps[c],
              (unsigned long long)a.bad);
        free(buf);
      }
  }
}

int main(void)
{
  make_seqs();
  test_eval();
  test_sym_all();
  test_gate_seq();
  test_symmetry();
  test_exhaustive();
  test_sample();
  printf(fails ? "FAILED: %d check(s)\n" : "all tests passed\n", fails);
  return fails != 0;
}
