/* KP-3: de Bruijn adjacency as a polynomial identity, tested on hashes.
 *
 *   u -> v (v = u shifted, c appended)  iff  x v(x) - u(x) = c x^k - u0
 * so, from the hash F(u) alone, the successors' hashes are among the 16 values
 *   cand(g, c) = (F(u) - b(g)) x^-1 + b(c) x^{k-1},   g = guessed first base of u, c = appended base.
 * The right guess g = u0 gives exactly F(v). A wrong guess gives the evaluation of u(x) - b(g) + ... whose constant
 * term u0 - b(g) is nonzero, so it equals F(w) for a stored w only if a nonzero polynomial of degree <= k vanishes at
 * x: probability <= k/q per (w, guess). So hash-only navigation over a set S has at most 12|S|k/q false successors
 * per node (union bound), and no false negatives. Predecessors are symmetric:
 *   cand_pred(l, p0) = x F(u) - b(l) x^k + b(p0),  l = guessed last base of u, p0 = the new first base.
 *
 * This program builds the k-mer set of a FASTA slice (k = 31, forward strand, a directed de Bruijn graph), then
 * compares hash-only navigation (fold of F looked up in a hash set) with exact navigation on 2-bit packed codes:
 * edges found, false positives and negatives, out- and in-degree agreement, non-branching node counts, and time.
 *   dbg <fasta> [max_bases=20000000]
 */
#include "field61.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define K 31
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
static inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; return x ^ (x >> 33); }

/* open-addressing set of nonzero 64-bit keys */
typedef struct { uint64_t* t; size_t mask; } hset;
static void hs_init(hset* h, size_t n) { size_t c = 1; while (c < 2 * n) c <<= 1; h->t = calloc(c, 8); h->mask = c - 1; }
static int hs_has(const hset* h, uint64_t k) { for (size_t i = mix(k) & h->mask;; i = (i + 1) & h->mask) { if (!h->t[i]) return 0; if (h->t[i] == k) return 1; } }
static int hs_add(hset* h, uint64_t k) { for (size_t i = mix(k) & h->mask;; i = (i + 1) & h->mask) { if (!h->t[i]) { h->t[i] = k; return 1; } if (h->t[i] == k) return 0; } }

int main(int argc, char** argv)
{
  if (argc < 2) { fprintf(stderr, "usage: dbg <fasta> [max_bases]\n"); return 2; }
  size_t maxb = argc > 2 ? strtoull(argv[2], 0, 10) : 20000000;
  FILE* f = fopen(argv[1], "r");
  char* seq = malloc(maxb + 1);
  size_t n = 0;
  char line[1 << 16];
  while (n < maxb && fgets(line, sizeof line, f)) {
    if (line[0] == '>') continue;
    for (char* p = line; *p && *p != '\n' && n < maxb; p++) seq[n++] = *p;
  }
  fclose(f);
  kp7_params kp;
  kp7_init(&kp, K, KP7_XREAL, 0x5eed);
  const uint64_t xinv = kp.xinv.re, x = kp.x.re, xk1 = f61_pow(x, K - 1), xk = f61_pow(x, K);
  const uint64_t cmask = (K == 32) ? ~0ULL : ((1ULL << (2 * K)) - 1);
  /* distinct k-mers: packed code (first base high) and F */
  size_t cap = n, m = 0;
  uint64_t* code = malloc(cap * 8);
  kp7_fe* F = malloc(cap * sizeof *F);
  hset cs; hs_init(&cs, n);
  kp7_roller r; kp7_start(&r, &kp, seq, n);
  uint64_t h, cur = 0; size_t curpos = (size_t)-2;
  double t0 = now();
  while (kp7_roll64(&r, KP7_OUT_FWD, &h)) {
    if (r.pos == curpos + 1) cur = ((cur << 2) | (uint64_t)base_code(seq[r.pos + K - 1])) & cmask;
    else { cur = 0; for (int t = 0; t < K; t++) cur = (cur << 2) | (uint64_t)base_code(seq[r.pos + t]); }
    curpos = r.pos;
    if (hs_add(&cs, cur + 1)) { code[m] = cur; F[m] = r.F; m++; } /* +1: keys must be nonzero */
  }
  hset fs; hs_init(&fs, m);
  uint64_t fold0 = 0;
  for (size_t i = 0; i < m; i++) { uint64_t fk = g_fold(F[i]); if (!fk) fold0++; hs_add(&fs, fk ? fk : 1); }
  printf("bases %zu, distinct %d-mers %zu (built in %.1f s); folds equal to 0: %llu\n", n, K, m, now() - t0, (unsigned long long)fold0);

  /* exact navigation on packed codes */
  uint64_t e_exact = 0, nb_exact = 0;
  double ta = now();
  unsigned char* outd = calloc(m, 1);
  unsigned char* ind = calloc(m, 1);
  for (size_t i = 0; i < m; i++) {
    for (uint64_t c = 0; c < 4; c++) if (hs_has(&cs, (((code[i] << 2) | c) & cmask) + 1)) { outd[i]++; e_exact++; }
    for (uint64_t c = 0; c < 4; c++) if (hs_has(&cs, ((code[i] >> 2) | (c << (2 * (K - 1)))) + 1)) ind[i]++;
    if (outd[i] == 1 && ind[i] == 1) nb_exact++;
  }
  double t_exact = now() - ta;

  /* hash-only navigation: 16 successor and 16 predecessor candidates from F alone */
  kp7_fe B[4], Bx[4], Bxk1[4], Bxk[4];
  for (int c = 0; c < 4; c++) { B[c] = g_base(c); Bx[c] = g_mul_real(B[c], xinv); Bxk1[c] = g_mul_real(B[c], xk1); Bxk[c] = g_mul_real(B[c], xk); }
  uint64_t tp = 0, fp_wrong = 0, fp_right = 0, fn = 0, deg_mismatch = 0, nb_hash = 0, pred_mismatch = 0;
  double tb = now();
  for (size_t i = 0; i < m; i++) {
    kp7_fe Fx = g_mul_real(F[i], xinv), Fxu = g_mul_real(F[i], x);
    int u0 = (int)(code[i] >> (2 * (K - 1))), ul = (int)(code[i] & 3);
    int od = 0, idg = 0;
    for (int g = 0; g < 4; g++) {
      kp7_fe base = g_sub(Fx, Bx[g]);
      for (int c = 0; c < 4; c++) {
        uint64_t fk = g_fold(g_add(base, Bxk1[c]));
        if (!hs_has(&fs, fk ? fk : 1)) continue;
        od++;
        if (g != u0) fp_wrong++;
        else if (hs_has(&cs, (((code[i] << 2) | (uint64_t)c) & cmask) + 1)) tp++;
        else fp_right++;
      }
      kp7_fe pb = g_sub(Fxu, Bxk[g]); /* g = guessed LAST base of u */
      for (int p0 = 0; p0 < 4; p0++) {
        uint64_t fk = g_fold(g_add(pb, B[p0]));
        if (!hs_has(&fs, fk ? fk : 1)) continue;
        idg++;
        if (g != ul) pred_mismatch++;
      }
    }
    if (od != outd[i] || idg != ind[i]) deg_mismatch++;
    if (od == 1 && idg == 1) nb_hash++;
  }
  double t_hash = now() - tb;
  fn = e_exact - tp;
  /* the provable bound: per node, 12 wrong-guess successor candidates, each hitting a stored k-mer w.p. <= |S| k / q */
  double bound = 12.0 * (double)m * (double)m * K / (double)F61_Q;
  printf("exact (packed codes): %llu edges, %llu non-branching nodes, %.2f s (%.1f ns/node)\n",
         (unsigned long long)e_exact, (unsigned long long)nb_exact, t_exact, t_exact / m * 1e9);
  printf("hash-only (16+16 candidates from F): true %llu, false via wrong guess %llu, false with right guess %llu, missed %llu,\n"
         "  predecessor hits via wrong guess %llu, nodes whose degrees differ %llu, non-branching %llu, %.2f s (%.1f ns/node)\n",
         (unsigned long long)tp, (unsigned long long)fp_wrong, (unsigned long long)fp_right, (unsigned long long)fn,
         (unsigned long long)pred_mismatch, (unsigned long long)deg_mismatch, (unsigned long long)nb_hash, t_hash, t_hash / m * 1e9);
  printf("provable bound on expected false successors over the whole set: <= %.3g; storage per node: packed %d bytes, hash-only %zu bytes (F)\n",
         bound, 8, sizeof(kp7_fe));
  return (fn || fp_wrong || fp_right || deg_mismatch || nb_hash != nb_exact) ? 1 : 0;
}
