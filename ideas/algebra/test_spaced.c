/* Tests for KP-14 spaced seeds (spaced.h). Exit 0 iff all pass. */
#include "spaced.h"
#include <stdio.h>
#include <stdlib.h>

static int fails = 0;
#define CHECK(c, ...) do { if (!(c)) { fails++; if (fails < 12) { printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } } while (0)
static uint64_t sm(uint64_t* s)
{
  uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static void rc(const char* s, char* o, unsigned k)
{
  for (unsigned t = 0; t < k; t++) {
    char c = s[k - 1 - t];
    o[t] = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N';
  }
  o[k] = 0;
}

static const char* MASKS31[] = {
  "1111111111111111111111111111111", /* all ones: must equal kp7 */
  "1110111011101110111011101110111", /* palindromic, 8 blocks */
  "1101100111010111011100101110101", /* not palindromic */
  "0111111111111111111111111111110", /* care set starts and ends inside */
  "1000000000000001000000000000001", /* sparse, single-base blocks */
  "1010101010101010101010101010101", /* alternating, 16 blocks */
};

int main(void)
{
  uint64_t s = 2024;
  kp7_params kp;
  kp7_init(&kp, 31, KP7_XREAL, 0x5eed);
  const size_t len = 60000;
  char* seq[2] = { malloc(len), malloc(len) };
  for (size_t i = 0; i < len; i++) {
    seq[0][i] = "ACGTN"[sm(&s) % 1000 == 0 ? 4 : sm(&s) % 4];
    seq[1][i] = "ACGTN"[sm(&s) % 30 == 0 ? 4 : sm(&s) % 4];
  }
  seq[1][0] = 'N';
  const int nm = sizeof MASKS31 / sizeof *MASKS31;
  for (int m = 0; m < nm; m++) {
    sp_params p;
    sp_init(&p, &kp, MASKS31[m]);
    for (int q = 0; q < 2; q++) {
      sp_roller r;
      sp_start(&r, &p, seq[q], len);
      uint64_t h, bad = 0, n = 0, badsym = 0;
      while (sp_roll(&r, &h)) {
        kp7_fe F, R;
        if (!sp_direct(&p, seq[q] + r.pos, &F, &R)) { bad++; continue; }
        if (F.re != r.F.re || F.im != r.F.im || R.re != r.R.re || R.im != r.R.im || h != sp_sym(&p, F, R)) bad++;
        if (n % 97 == 0) { /* strand symmetry on a sample of windows */
          char w[32], wr[32];
          memcpy(w, seq[q] + r.pos, 31); w[31] = 0;
          rc(w, wr, 31);
          kp7_fe F2, R2;
          sp_direct(&p, wr, &F2, &R2);
          if (sp_sym(&p, F2, R2) != h) badsym++;
        }
        n++;
      }
      CHECK(!bad && n > 10000, "rolling != direct: mask %d seq %d (%llu bad of %llu)", m, q, (unsigned long long)bad, (unsigned long long)n);
      CHECK(!badsym, "not strand-symmetric: mask %d seq %d (%llu)", m, q, (unsigned long long)badsym);
    }
  }
  /* all-ones mask == kp7_roll64 symmetric output, every window */
  {
    sp_params p;
    sp_init(&p, &kp, MASKS31[0]);
    for (int q = 0; q < 2; q++) {
      sp_roller r; kp7_roller r0;
      sp_start(&r, &p, seq[q], len);
      kp7_start(&r0, &kp, seq[q], len);
      uint64_t h, h0, bad = 0;
      for (;;) {
        int a = sp_roll(&r, &h), b = kp7_roll64(&r0, KP7_OUT_SYM, &h0);
        if (a != b) { bad++; break; }
        if (!a) break;
        if (h != h0 || r.pos != r0.pos) bad++;
      }
      CHECK(!bad, "all-ones mask differs from kp7_roll64 SYM on seq %d (%llu)", q, (unsigned long long)bad);
    }
  }
  /* exhaustive k = 7: distinct symmetric values == exact classes of {masked(s), masked(rc s)} */
  {
    const char* masks7[] = { "1101011", "1111111", "1011101", "1100111" };
    for (int m = 0; m < 4; m++) {
      kp7_params kp7k;
      kp7_init(&kp7k, 7, KP7_XREAL, 77);
      sp_params p;
      sp_init(&p, &kp7k, masks7[m]);
      size_t n = 1u << 14;
      uint64_t* hv = malloc(n * sizeof *hv);
      uint64_t* cls = malloc(n * sizeof *cls);
      for (size_t code = 0; code < n; code++) {
        char w[8], wr[8];
        for (int t = 0; t < 7; t++) w[t] = "ACGT"[(code >> (2 * t)) & 3];
        w[7] = 0;
        rc(w, wr, 7);
        kp7_fe F, R;
        sp_direct(&p, w, &F, &R);
        hv[code] = sp_sym(&p, F, R);
        uint64_t a = 0, b = 0; /* masked care bases of s and of rc(s), as 2-bit strings */
        for (int t = 0; t < 7; t++) if (masks7[m][t] == '1') { a = a * 4 + base_code(w[t]); b = b * 4 + base_code(wr[t]); }
        cls[code] = a < b ? (a << 32) | b : (b << 32) | a;
      }
      int cmp(const void* x, const void* y);
      size_t dh = 0, dc = 0;
      qsort(hv, n, sizeof *hv, cmp); qsort(cls, n, sizeof *cls, cmp);
      for (size_t i = 0; i < n; i++) { if (!i || hv[i] != hv[i - 1]) dh++; if (!i || cls[i] != cls[i - 1]) dc++; }
      printf("  k=7 mask %s: %zu distinct hashes, %zu exact classes\n", masks7[m], dh, dc);
      CHECK(dh == dc, "mask %s: %zu hashes vs %zu classes", masks7[m], dh, dc);
      free(hv); free(cls);
    }
  }
  free(seq[0]); free(seq[1]);
  printf(fails ? "FAILED: %d check(s)\n" : "all tests passed\n", fails);
  return fails != 0;
}
int cmp(const void* x, const void* y) { uint64_t a = *(const uint64_t*)x, b = *(const uint64_t*)y; return a < b ? -1 : a > b; }
