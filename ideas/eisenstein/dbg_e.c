/* KP-3 Eisenstein revisit: de Bruijn navigation from 8-byte Eisenstein hashes (hash/kp7e.h), against packed codes.
 *
 * Same identity as ideas/algebra/dbg.c, now in F_q (q = 2^61 - 1) with A = 1, C = w, G = -w, T = -1:
 *   successor candidates   (F(u) - val(g)) x^-1 + val(c) x^{k-1}     g = guessed first base, c = appended base
 *   predecessor candidates  x F(u) - val(l) x^k + val(p0)            l = guessed last base, p0 = new first base
 * The right guess gives exactly F of the neighbour. A wrong guess leaves a nonzero constant term, so a false hit on a
 * stored k-mer needs a nonzero polynomial of degree <= k to vanish at the random x: <= k/q per (node, guess, w), and
 * no false negatives. The hash is ONE field element, 8 bytes (the Gaussian version in dbg.c needed 16).
 * Baseline: exact navigation on 2-bit packed codes, 64-bit for k <= 32 and 128-bit for k <= 64 (16 bytes per node).
 *   dbg_e <fasta> [max_bases=20000000]     (k is the compile-time DBG_K: make dbg_e31 dbg_e47 dbg_e63)
 */
#include "../../hash/kp7e.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef DBG_K
#define DBG_K 31
#endif
#define K DBG_K
#if K <= 32
typedef uint64_t pk_t;
#else
typedef unsigned __int128 pk_t;
#endif

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
static inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; return x ^ (x >> 33); }
static inline uint64_t pk_mix(pk_t v) { return mix((uint64_t)v ^ mix((uint64_t)(v >> 32 >> 32))); }

/* open addressing; key 0 means empty, so callers store key + 1 */
typedef struct { uint64_t* t; size_t mask; } hset;
static void hs_init(hset* h, size_t n) { size_t c = 1; while (c < 2 * n) c <<= 1; h->t = calloc(c, 8); h->mask = c - 1; }
static int hs_has(const hset* h, uint64_t k) { for (size_t i = mix(k) & h->mask;; i = (i + 1) & h->mask) { if (!h->t[i]) return 0; if (h->t[i] == k) return 1; } }
static int hs_add(hset* h, uint64_t k) { for (size_t i = mix(k) & h->mask;; i = (i + 1) & h->mask) { if (!h->t[i]) { h->t[i] = k; return 1; } if (h->t[i] == k) return 0; } }
typedef struct { pk_t* t; size_t mask; } pset;
static void ps_init(pset* h, size_t n) { size_t c = 1; while (c < 2 * n) c <<= 1; h->t = calloc(c, sizeof(pk_t)); h->mask = c - 1; }
static int ps_has(const pset* h, pk_t k) { for (size_t i = pk_mix(k) & h->mask;; i = (i + 1) & h->mask) { if (!h->t[i]) return 0; if (h->t[i] == k) return 1; } }
static int ps_add(pset* h, pk_t k) { for (size_t i = pk_mix(k) & h->mask;; i = (i + 1) & h->mask) { if (!h->t[i]) { h->t[i] = k; return 1; } if (h->t[i] == k) return 0; } }

#define Q KP7E_Q
static inline uint64_t sub_q(uint64_t a, uint64_t b) { return a >= b ? a - b : a + Q - b; }

int main(int argc, char** argv)
{
  if (argc < 2) { fprintf(stderr, "usage: dbg_e <fasta> [max_bases]\n"); return 2; }
  size_t maxb = argc > 2 ? strtoull(argv[2], 0, 10) : 20000000;
  FILE* f = fopen(argv[1], "r");
  if (!f) { perror(argv[1]); return 2; }
  char* seq = malloc(maxb + 1);
  size_t n = 0;
  char line[1 << 16];
  while (n < maxb && fgets(line, sizeof line, f)) {
    if (line[0] == '>') continue;
    for (char* p = line; *p && *p != '\n' && n < maxb; p++) seq[n++] = *p;
  }
  fclose(f);
  kp7e_params kp;
  kp7e_init(&kp, K, 0x5eed);
  const pk_t one = 1;
  const pk_t cmask = (K == 64) ? ~(pk_t)0 : ((one << (2 * K)) - 1);
  size_t m = 0;
  pk_t* code = malloc(n * sizeof(pk_t));
  uint64_t* F = malloc(n * 8);
  pset cs; ps_init(&cs, n);
  kp7e_roller r; kp7e_start(&r, &kp, seq, n);
  pk_t cur = 0; size_t curpos = (size_t)-2;
  double t0 = now();
  while (kp7e_roll(&r)) {
    if (r.pos == curpos + 1) cur = ((cur << 2) | (pk_t)kp7e_code(seq[r.pos + K - 1])) & cmask;
    else { cur = 0; for (int t = 0; t < K; t++) cur = (cur << 2) | (pk_t)kp7e_code(seq[r.pos + t]); }
    curpos = r.pos;
    if (ps_add(&cs, cur + 1)) { code[m] = cur; F[m] = r.F; m++; }
  }
  hset fs; hs_init(&fs, m);
  uint64_t f0 = 0;
  for (size_t i = 0; i < m; i++) { if (!F[i]) f0++; hs_add(&fs, F[i] + 1); } /* F < q < 2^61, so F + 1 never wraps */
  printf("k %d, bases %zu, distinct k-mers %zu (built in %.1f s); hashes equal to 0: %llu\n", K, n, m, now() - t0,
         (unsigned long long)f0);

  /* exact navigation on packed codes */
  uint64_t e_exact = 0, nb_exact = 0;
  unsigned char* outd = calloc(m, 1);
  unsigned char* ind = calloc(m, 1);
  double ta = now();
  for (size_t i = 0; i < m; i++) {
    for (pk_t c = 0; c < 4; c++) if (ps_has(&cs, (((code[i] << 2) | c) & cmask) + 1)) { outd[i]++; e_exact++; }
    for (pk_t c = 0; c < 4; c++) if (ps_has(&cs, ((code[i] >> 2) | (c << (2 * (K - 1)))) + 1)) ind[i]++;
    if (outd[i] == 1 && ind[i] == 1) nb_exact++;
  }
  double t_exact = now() - ta;

  /* hash-only navigation */
  uint64_t Bx[4], Bxk1[4], Bxk[4];
  for (int c = 0; c < 4; c++) { Bx[c] = kp7e_mul(kp.val[c], kp.xinv); Bxk1[c] = kp7e_mul(kp.val[c], kp.xk1); Bxk[c] = kp7e_mul(kp.val[c], kp.xk); }
  uint64_t tp = 0, fp_wrong = 0, fp_right = 0, deg_mismatch = 0, nb_hash = 0, pred_mismatch = 0;
  double tb = now();
  for (size_t i = 0; i < m; i++) {
    uint64_t Fx = kp7e_mul(F[i], kp.xinv), Fxu = kp7e_mul(F[i], kp.x);
    int u0 = (int)(code[i] >> (2 * (K - 1))), ul = (int)(code[i] & 3);
    int od = 0, idg = 0;
    for (int g = 0; g < 4; g++) {
      uint64_t base = sub_q(Fx, Bx[g]);
      for (int c = 0; c < 4; c++) {
        uint64_t fk = kp7e_add(base, Bxk1[c]);
        if (!hs_has(&fs, fk + 1)) continue;
        od++;
        if (g != u0) fp_wrong++;
        else if (ps_has(&cs, (((code[i] << 2) | (pk_t)c) & cmask) + 1)) tp++;
        else fp_right++;
      }
      uint64_t pb = sub_q(Fxu, Bxk[g]);
      for (int p0 = 0; p0 < 4; p0++) {
        uint64_t fk = kp7e_add(pb, kp.val[p0]);
        if (!hs_has(&fs, fk + 1)) continue;
        idg++;
        if (g != ul) pred_mismatch++;
      }
    }
    if (od != outd[i] || idg != ind[i]) deg_mismatch++;
    if (od == 1 && idg == 1) nb_hash++;
  }
  double t_hash = now() - tb;
  uint64_t fn = e_exact - tp;
  double bound = 12.0 * (double)m * (double)m * K / (double)Q;
  printf("exact (packed codes, %zu B/node): %llu edges, %llu non-branching, %.2f s (%.1f ns/node)\n", sizeof(pk_t),
         (unsigned long long)e_exact, (unsigned long long)nb_exact, t_exact, t_exact / m * 1e9);
  printf("hash-only (Eisenstein F, 8 B/node): true %llu, false via wrong guess %llu, false with right guess %llu, missed %llu,\n"
         "  predecessor hits via wrong guess %llu, nodes whose degrees differ %llu, non-branching %llu, %.2f s (%.1f ns/node)\n",
         (unsigned long long)tp, (unsigned long long)fp_wrong, (unsigned long long)fp_right, (unsigned long long)fn,
         (unsigned long long)pred_mismatch, (unsigned long long)deg_mismatch, (unsigned long long)nb_hash, t_hash, t_hash / m * 1e9);
  printf("provable bound on expected false successors: <= %.3g; key-table bytes per node: packed %zu, hash %zu\n",
         bound, sizeof(pk_t), sizeof(uint64_t));
  return (fn || fp_wrong || fp_right || deg_mismatch || nb_hash != nb_exact) ? 1 : 0;
}
