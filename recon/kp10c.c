/* KP-10 in C, end to end on two genomes, k-mers AS WRITTEN (no canonicalization, no 2-bit packing):
 *
 *   kp10c --a A.fa --b B.fa --k 31 --m 64 --c 2 --dest D1[,D2,...] [--seed-x X|os] [--seed-lam L|os]
 *         [--max-rounds 16] [--threads T] [--truth-amb amb.u64 --truth-bma bma.u64] [--label s]
 *         [--rule pow2|exact] [--load 0.8] [--buckets N] [--oneway 0|1] [--truth-kmc-amb amb.txt --truth-kmc-bma bma.txt]
 *
 * One JSON line per --dest value on stdout; progress on stderr.
 *
 * 1. Element sets. Each genome is hashed by the Eisenstein KP-7 hash S = (F + lam)(R + lam) over F_q, q = 2^61 - 1
 *    (hash/kp7e.c, kp7e_sym_seq: the C path the tests check bit for bit against the reference roller). S(s) = S(rc s),
 *    so the set of distinct values is the set of double-stranded k-mer classes, as long as no two classes collide
 *    (bound 2(k-1)/q per pair). The distinct values are sorted (LSD radix) and deduplicated.
 * 2. Sketch. A bucket at resolution `bits` is the range of values sharing their top `bits` bits (v >> (61 - bits)),
 *    as recon/aswritten.py keys them: nothing else is hashed. Per bucket: the size and power sums S_1..S_{m+c}.
 *    Each side computes its sketch from its own values only. Bytes: (m + c + 1) entries x 8 B per bucket per round, ONE
 *    side's message (as recon/recon61.py counts, and as one km-peeler IBLT file is one side's message).
 * 3. Decode, per bucket, ported from recon/recon61.py:decode (Newton's identities, then a Pade step of degrees
 *    dA = (mbar + delta)/2, dB = (mbar - delta)/2, free variables 0). The polynomials P (roots A\B) and Q (roots B\A)
 *    come from the two sketches alone. Each side then flags ITS OWN elements that are roots of its polynomial and
 *    checks itself with what it holds: side A requires #flagged = deg P_rev and, for j = m+1..m+c,
 *    S_j = psum_j(flagged A) - p_j(Q), where p_j(Q) is the j-th power sum of Q's roots computed from Q's
 *    coefficients by Newton's identities (no factoring); side B symmetrically. A bucket is accepted when both sides
 *    pass and deg P_rev - deg Q_rev = delta. (recon61.py checked psum(flagged A) - psum(flagged B) in one process;
 *    when both counts pass the two checks coincide, but each half here uses one side's data only, so the decision
 *    needs no oracle. Each party still recovers only ITS OWN differing k-mers; the other side's appear as roots of
 *    an unfactored polynomial.) A bucket with equal sizes and equal sketches is taken to have no difference.
 * 4. Rounds. A failed bucket is split on the next value bit and both halves are sketched and sent again.
 * 5. Back to k-mers. Each side rescans its genome and emits the k-mers whose S it flagged, then canonicalizes them
 *    (k <= 32) ONLY to compare with an independent truth (KMC 3 / recon/kmtools.c canonical set difference).
 *
 * Initial bucket count from the difference estimate --dest: the smallest bits with dest / 2^bits <= 0.8 m (this rule
 * reproduces the 2^16 and 2^8 buckets of the 2026-11-10 runs at the true D; it is fixed before any run here).
 * Timing: wall and per-process CPU (getrusage) per phase. With --threads 1 everything is one thread.
 *
 * Round-3 additions (2026-11-12, results/round3_owed_2026-11-12.md); every default reproduces the round-2 behaviour.
 *  - Buckets are value RANGES [lo, hi) of [0, 2^61). Round 1 has N0 ranges, lo_i = floor(i 2^61 / N0); a failed range
 *    is split at its midpoint. For N0 = 2^b these are exactly the top-b-bit buckets above (same slices, same bytes).
 *    --rule pow2 (default): N0 = 2^b as above. --rule exact: N0 = ceil(dest / (load m)), no power of two.
 *    --load sets the 0.8. --buckets N forces N0 = N for every dest; a dest written "D@N" forces N0 = N for that one.
 *  - --oneway 1: ONE sketch direction. The party with the larger set (sizes exchanged first, 2 x 4 B) decodes; the
 *    other sends its sketch each round. Acceptance uses the decoder's own check only (its flagged roots of P and
 *    p_j(Q) from Q's coefficients). After each round the decoder sends the failed bucket indices and, for each
 *    accepted bucket with deg Q > 0, its index and Q's eB coefficients; the sender flags its own roots of Q (a bucket
 *    where the count disagrees is reported failed and split in the next round) and returns its differing k-mers.
 *    So the decoder ends with BOTH sides' differences as k-mers, as one IBLT listing gives them.
 *  - Byte accounting, both protocols: "stored" = 8-byte entries, 4-byte indices, 8-byte (k <= 32) or 16-byte k-mers;
 *    "packed" = 61-bit field elements, a 32-bit bucket count, ceil(log2 N)-bit indices (or an N-bit bitmap if
 *    smaller), 2k-bit k-mers. Two-way = both sketches; one-way = every message in either direction.
 *  - k <= 64: recovered k-mers are canonical 128-bit codes; the truth can be KMC dump text (--truth-kmc-*).
 */
#define _GNU_SOURCE
#include "kp7e.h"
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <time.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;
#define Q KP7E_Q
#define MAXM 256

/* ---------- F_q, q = 2^61 - 1 (canonical in [0, q)) ---------- */
static inline u64 red(u128 x)
{
  u64 r = ((u64)x & Q) + (u64)(x >> 61);
  r = (r & Q) + (r >> 61);
  return r == Q ? 0 : r;
}
static inline u64 mulq(u64 a, u64 b) { return red((u128)a * b); }
static inline u64 addq(u64 a, u64 b) { u64 s = a + b; return s >= Q ? s - Q : s; }
static inline u64 subq(u64 a, u64 b) { return a >= b ? a - b : a + Q - b; }
static u64 powq(u64 a, u64 e) { u64 r = 1; while (e) { if (e & 1) r = mulq(r, a); a = mulq(a, a); e >>= 1; } return r; }
static inline u64 invq(u64 a) { return powq(a, Q - 2); }
static u64 INV[MAXM + 2]; /* INV[j] = 1/j */

/* ---------- timing ---------- */
static double wall(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
static double cpu(void)
{
  struct rusage u;
  getrusage(RUSAGE_SELF, &u);
  return u.ru_utime.tv_sec + u.ru_utime.tv_usec * 1e-6 + u.ru_stime.tv_sec + u.ru_stime.tv_usec * 1e-6;
}
typedef struct { double w, c; } tick;
static tick now(void) { tick t = {wall(), cpu()}; return t; }
static tick since(tick a) { tick b = now(); tick d = {b.w - a.w, b.c - a.c}; return d; }

/* ---------- genome ---------- */
typedef struct { char* s; size_t n; } genome;

/* whole FASTA in one buffer, records joined by one 'N' so no k-mer spans two records */
static genome read_fasta(const char* path)
{
  FILE* f = fopen(path, "r");
  if (!f) { perror(path); exit(1); }
  struct stat st;
  fstat(fileno(f), &st);
  genome g = {malloc((size_t)st.st_size + 16), 0};
  char* line = NULL;
  size_t lcap = 0;
  ssize_t l;
  while ((l = getline(&line, &lcap, f)) >= 0) {
    if (l && line[0] == '>') { g.s[g.n++] = 'N'; continue; }
    while (l && (line[l - 1] == '\n' || line[l - 1] == '\r')) l--;
    memcpy(g.s + g.n, line, (size_t)l);
    g.n += (size_t)l;
  }
  free(line);
  fclose(f);
  return g;
}

/* ---------- hashing a genome into a value array ---------- */
typedef struct { u64* v; size_t n, cap; } vec;
static void vpush_n(vec* a, const u64* x, size_t n)
{
  if (a->n + n > a->cap) {
    while (a->n + n > a->cap) a->cap = a->cap ? a->cap * 2 : 1 << 24;
    a->v = realloc(a->v, a->cap * 8);
    if (!a->v) { perror("realloc"); exit(1); }
  }
  memcpy(a->v + a->n, x, n * 8);
  a->n += n;
}
static void emit_all(void* ctx, size_t start, const u64* h, size_t n) { (void)start; vpush_n((vec*)ctx, h, n); }

static size_t sort_unique(u64* a, size_t n)
{
  u64* b = malloc(n * 8 + 8);
  size_t* cnt = malloc(65536 * sizeof(size_t));
  if (!b || !cnt) { perror("malloc"); exit(1); }
  for (int pass = 0; pass < 4; pass++) {
    int sh = 16 * pass;
    memset(cnt, 0, 65536 * sizeof(size_t));
    for (size_t i = 0; i < n; i++) cnt[(a[i] >> sh) & 0xFFFF]++;
    size_t s = 0;
    for (int d = 0; d < 65536; d++) { size_t c = cnt[d]; cnt[d] = s; s += c; }
    for (size_t i = 0; i < n; i++) b[cnt[(a[i] >> sh) & 0xFFFF]++] = a[i];
    u64* t = a; a = b; b = t;
  }
  free(b);
  free(cnt);
  size_t m = 0;
  for (size_t i = 0; i < n; i++)
    if (m == 0 || a[i] != a[m - 1]) a[m++] = a[i];
  return m;
}

/* ---------- one side ---------- */
typedef struct {
  genome g;
  u64* v;    /* distinct S values, sorted */
  size_t n;
  u64 kmers; /* all k-mer windows hashed */
  tick t_hash, t_dedup;
} side_t;

static void build_side(side_t* s, const kp7e_params* p, const char* path)
{
  s->g = read_fasta(path);
  tick t0 = now();
  vec all = {0};
  u64* hb = malloc(4096 * 8);
  s->kmers = kp7e_sym_seq(p, s->g.s, s->g.n, 2, hb, 4096, emit_all, &all);
  free(hb);
  s->t_hash = since(t0);
  t0 = now();
  s->n = sort_unique(all.v, all.n);
  s->v = realloc(all.v, s->n * 8 + 8);
  s->t_dedup = since(t0);
}

static inline size_t lower(const u64* v, size_t n, u64 key)
{
  size_t lo = 0, hi = n;
  while (lo < hi) { size_t mid = (lo + hi) / 2; if (v[mid] < key) lo = mid + 1; else hi = mid; }
  return lo;
}
/* a bucket is the value range [lo, hi) of [0, 2^61) */
typedef struct { u64 lo, hi; } range;
#define TOP ((u64)1 << 61)
static void slice(const side_t* s, range r, size_t* lo, size_t* hi)
{
  *lo = lower(s->v, s->n, r.lo);
  *hi = r.hi >= TOP ? s->n : lower(s->v, s->n, r.hi);
}
static unsigned ceil_log2(u64 n) { unsigned b = 0; while (b < 64 && ((u64)1 << b) < n) b++; return b; }

/* power sums S_1..S_M of e[0..n) into out[0..M) (out[j-1] = S_j); 4 independent chains for ILP */
static void psums(const u64* e, size_t n, unsigned M, u64* out)
{
  u64 acc[4][MAXM + 2];
  memset(acc, 0, sizeof acc);
  size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    u64 h0 = e[i], h1 = e[i + 1], h2 = e[i + 2], h3 = e[i + 3];
    u64 p0 = h0, p1 = h1, p2 = h2, p3 = h3;
    for (unsigned j = 0; j < M; j++) {
      acc[0][j] = addq(acc[0][j], p0); p0 = mulq(p0, h0);
      acc[1][j] = addq(acc[1][j], p1); p1 = mulq(p1, h1);
      acc[2][j] = addq(acc[2][j], p2); p2 = mulq(p2, h2);
      acc[3][j] = addq(acc[3][j], p3); p3 = mulq(p3, h3);
    }
  }
  for (; i < n; i++) {
    u64 h = e[i], p = h;
    for (unsigned j = 0; j < M; j++) { acc[0][j] = addq(acc[0][j], p); p = mulq(p, h); }
  }
  for (unsigned j = 0; j < M; j++) out[j] = addq(addq(acc[0][j], acc[1][j]), addq(acc[2][j], acc[3][j]));
}

/* ---------- decode ---------- */
/* Gauss-Jordan mod q, free variables 0; returns 0 if inconsistent. A is n x (cols + 1), row-major. */
static int solve(u64* A, int n, int cols, u64* x)
{
  int r = 0, piv[MAXM];
  for (int c = 0; c < cols && r < n; c++) {
    int p = -1;
    for (int i = r; i < n; i++) if (A[i * (cols + 1) + c]) { p = i; break; }
    if (p < 0) continue;
    if (p != r) for (int j = 0; j <= cols; j++) { u64 t = A[r * (cols + 1) + j]; A[r * (cols + 1) + j] = A[p * (cols + 1) + j]; A[p * (cols + 1) + j] = t; }
    u64 inv = invq(A[r * (cols + 1) + c]);
    for (int j = 0; j <= cols; j++) A[r * (cols + 1) + j] = mulq(A[r * (cols + 1) + j], inv);
    for (int i = 0; i < n; i++) {
      u64 f = A[i * (cols + 1) + c];
      if (i == r || !f) continue;
      for (int j = 0; j <= cols; j++) A[i * (cols + 1) + j] = subq(A[i * (cols + 1) + j], mulq(f, A[r * (cols + 1) + j]));
    }
    piv[r++] = c;
  }
  for (int i = r; i < n; i++) if (A[i * (cols + 1) + cols]) return 0;
  for (int c = 0; c < cols; c++) x[c] = 0;
  for (int i = 0; i < r; i++) x[piv[i]] = A[i * (cols + 1) + cols];
  return 1;
}

/* flag own elements e[0..n) that are zeros of P(x) = sum_{i=0..d} rev[i] x^(d-i) (rev = P_rev, low first);
 * writes them to out (at most cap); returns the count, or cap + 1 if more */
static size_t flag_roots(const u64* rev, int d, const u64* e, size_t n, u64* out, size_t cap)
{
  size_t f = 0;
  for (size_t i = 0; i < n; i++) {
    u64 acc = 0, x = e[i];
    for (int j = 0; j <= d; j++) acc = addq(mulq(acc, x), rev[j]);
    if (!acc) { if (f < cap) out[f] = x; f++; if (f > cap) return f; }
  }
  return f;
}

/* power sums p_1..p_M of the roots of f with f_rev = 1 + c1 t + ... + cd t^d = prod (1 - a t):
 * p_j = -j c_j - sum_{i=1}^{j-1} c_i p_{j-i} */
static void root_psums(const u64* rev, int d, unsigned M, u64* p)
{
  for (unsigned j = 1; j <= M; j++) {
    u64 s = (int)j <= d ? mulq(j % Q, rev[j]) : 0;
    for (unsigned i = 1; i < j && (int)i <= d; i++) s = addq(s, mulq(rev[i], p[j - i - 1]));
    p[j - 1] = subq(0, s);
  }
}
static u64 psum_j(const u64* e, size_t n, unsigned j)
{
  u64 s = 0;
  for (size_t i = 0; i < n; i++) s = addq(s, powq(e[i], j));
  return s;
}

/* returns 1 and fills fa/fb (own flagged elements) on success; 0 on overflow / failed check.
 * one_sided: accept on A's check alone (A decodes); B's flagging of Q's roots is then the sender's later check, and
 * *late = 1 when it disagrees (the bucket must be split). *eBout = deg Q (the number of B-side differences). */
static int decode(const u64* SA, const u64* SB, long nA, long nB, int m, int c, const u64* ownA, size_t lenA,
                  const u64* ownB, size_t lenB, u64* fa, size_t* nfa, u64* fb, size_t* nfb, int one_sided, int* late,
                  int* eBout)
{
  *late = 0;
  *eBout = 0;
  u64 S[MAXM + 2];
  for (int j = 0; j < m + c; j++) S[j] = subq(SA[j], SB[j]);
  long delta = nA - nB;
  int mbar = ((m - delta) % 2 == 0) ? m : m - 1;
  if (labs(delta) > mbar) return 0;
  int dA = (int)((mbar + delta) / 2), dB = (int)((mbar - delta) / 2);
  u64 r[MAXM + 2];
  r[0] = 1;
  for (int k = 1; k <= mbar; k++) {
    u64 acc = 0;
    for (int j = 1; j <= k; j++) acc = addq(acc, mulq(S[j - 1], r[k - j]));
    r[k] = mulq(subq(0, acc), INV[k]);
  }
  u64 qr[MAXM + 2], pr[MAXM + 2];
  qr[0] = 1;
  if (dB) {
    static __thread u64 A[MAXM * (MAXM + 1)];
    for (int row = 0; row < dB; row++) {
      int n = dA + 1 + row;
      for (int j = 1; j <= dB; j++) A[row * (dB + 1) + (j - 1)] = n - j >= 0 ? r[n - j] : 0;
      A[row * (dB + 1) + dB] = subq(0, r[n]);
    }
    if (!solve(A, dB, dB, qr + 1)) return 0;
  }
  for (int n = 0; n <= dA; n++) {
    u64 s = 0;
    for (int j = 0; j <= (n < dB ? n : dB); j++) s = addq(s, mulq(qr[j], r[n - j]));
    pr[n] = s;
  }
  int eA = dA, eB = dB; /* true degrees of P_rev and Q_rev = numbers of nonzero roots */
  while (eA > 0 && !pr[eA]) eA--;
  while (eB > 0 && !qr[eB]) eB--;
  if (eA - eB != delta) return 0;
  *eBout = eB;
  /* side A: its own roots of P (as P(x) = x^eA P_rev(1/x), descending = pr[0..eA]) */
  *nfa = flag_roots(pr, eA, ownA, lenA, fa, (size_t)m);
  if (*nfa != (size_t)eA) return 0;
  u64 pP[MAXM + 2], pQ[MAXM + 2];
  root_psums(pr, eA, (unsigned)(m + c), pP);
  root_psums(qr, eB, (unsigned)(m + c), pQ);
  /* side A's check (own flagged + Q's coefficients) */
  for (int j = m + 1; j <= m + c; j++)
    if (subq(psum_j(fa, *nfa, (unsigned)j), pQ[j - 1]) != S[j - 1]) return 0;
  *nfb = flag_roots(qr, eB, ownB, lenB, fb, (size_t)m);
  if (one_sided) { /* A has accepted; B's flagging is the sender's check after Q's coefficients arrive */
    *late = *nfb != (size_t)eB;
    return 1;
  }
  if (*nfb != (size_t)eB) return 0;
  /* side B's check (P's coefficients + own flagged) */
  for (int j = m + 1; j <= m + c; j++)
    if (subq(pP[j - 1], psum_j(fb, *nfb, (unsigned)j)) != S[j - 1]) return 0;
  return 1;
}

/* ---------- the protocol ---------- */
#define MAXR 64
typedef struct {
  u64 *amb, *bma; /* recovered values */
  size_t namb, nbma;
  int rounds, unresolved, bits0; /* bits0 = log2 N0 when N0 is a power of two, else -1 */
  u64 buckets0, splits, entries, buckets_with_difference, late_failures;
  u64 sent[MAXR], failed[MAXR];
  /* one-way replies (decoder -> sender: failed indices, Q coefficients; sender -> decoder: its k-mers) */
  u64 fail_idx_bytes, fail_idx_bits, q_buckets, q_coeffs, q_bytes, q_bits;
  tick t_sketchA, t_sketchB, t_decode;
} result;

static int cmp_u64(const void* a, const void* b) { u64 x = *(const u64*)a, y = *(const u64*)b; return x < y ? -1 : x > y; }

/* A decodes. Two-sided (oneway = 0): both sides' checks, as in round 2. One-sided: A's check, then B's flagging. */
static result reconcile(const side_t* A, const side_t* B, u64 n0, int m, int c, int max_rounds, int oneway)
{
  result R = {0};
  R.buckets0 = n0;
  R.bits0 = (n0 & (n0 - 1)) ? -1 : (int)ceil_log2(n0);
  const unsigned M = (unsigned)(m + c);
  size_t np = (size_t)n0;
  range* pending = malloc(np * sizeof(range));
  for (size_t i = 0; i < np; i++) {
    pending[i].lo = (u64)(((u128)i << 61) / n0);
    pending[i].hi = (u64)(((u128)(i + 1) << 61) / n0);
  }
  vec oa = {0}, ob = {0};
  while (np) {
    R.rounds++;
    if (R.rounds <= MAXR) R.sent[R.rounds - 1] = np;
    R.entries += (u64)np * (M + 1);
    u64* skA = malloc(np * M * 8);
    u64* skB = malloc(np * M * 8);
    long *nA = malloc(np * sizeof(long)), *nB = malloc(np * sizeof(long));
    size_t *loA = malloc(np * 8), *hiA = malloc(np * 8), *loB = malloc(np * 8), *hiB = malloc(np * 8);
    /* side A's sketch, from A's values only */
    tick t0 = now();
#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < np; i++) {
      slice(A, pending[i], &loA[i], &hiA[i]);
      nA[i] = (long)(hiA[i] - loA[i]);
      psums(A->v + loA[i], hiA[i] - loA[i], M, skA + i * M);
    }
    tick d = since(t0);
    R.t_sketchA.w += d.w; R.t_sketchA.c += d.c;
    t0 = now();
#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < np; i++) {
      slice(B, pending[i], &loB[i], &hiB[i]);
      nB[i] = (long)(hiB[i] - loB[i]);
      psums(B->v + loB[i], hiB[i] - loB[i], M, skB + i * M);
    }
    d = since(t0);
    R.t_sketchB.w += d.w; R.t_sketchB.c += d.c;
    /* decode */
    t0 = now();
    unsigned char* ok = calloc(np, 1);
    int* eB = calloc(np, sizeof(int));
    u64* fa = malloc(np * (size_t)m * 8);
    u64* fb = malloc(np * (size_t)m * 8);
    size_t *nfa = calloc(np, sizeof(size_t)), *nfb = calloc(np, sizeof(size_t));
    u64 withdiff = 0, late = 0, qb = 0, qc = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : withdiff, late, qb, qc)
    for (size_t i = 0; i < np; i++) {
      if (nA[i] == nB[i] && !memcmp(skA + i * M, skB + i * M, M * 8)) { ok[i] = 1; continue; }
      withdiff++;
      int lt = 0;
      ok[i] = (unsigned char)decode(skA + i * M, skB + i * M, nA[i], nB[i], m, c, A->v + loA[i], hiA[i] - loA[i],
                                    B->v + loB[i], hiB[i] - loB[i], fa + i * (size_t)m, &nfa[i], fb + i * (size_t)m,
                                    &nfb[i], oneway, &lt, &eB[i]);
      if (ok[i] && eB[i] > 0) { qb++; qc += (u64)eB[i]; } /* one-way: Q's coefficients go to the sender */
      if (ok[i] && lt) { late++; ok[i] = 0; }            /* the sender's flagging disagreed: split it */
    }
    if (R.rounds == 1) R.buckets_with_difference = withdiff;
    size_t nover = 0;
    for (size_t i = 0; i < np; i++) {
      if (ok[i]) {
        vpush_n(&oa, fa + i * (size_t)m, nfa[i]);
        vpush_n(&ob, fb + i * (size_t)m, nfb[i]);
      } else {
        pending[nover++] = pending[i];
      }
    }
    if (oneway) {
      unsigned ib = ceil_log2(np);
      R.late_failures += late;
      R.q_buckets += qb;
      R.q_coeffs += qc;
      R.q_bytes += qb * 4 + qc * 8;
      R.q_bits += qb * ib + qc * 61;
      /* failed indices for the sender: 4 B each stored; packed, the smaller of a list and an np-bit bitmap. A late
       * failure is reported by the sender itself, with its k-mer reply, so it is not in the decoder's list. */
      u64 nf = nover - late;
      R.fail_idx_bytes += nf * 4;
      R.fail_idx_bits += nf * ib < (u64)np ? nf * ib : (u64)np;
    }
    if (R.rounds <= MAXR) R.failed[R.rounds - 1] = nover;
    d = since(t0);
    R.t_decode.w += d.w; R.t_decode.c += d.c;
    fprintf(stderr, "  round %d: %zu buckets sent, %llu with a difference, %zu failed (%llu late)\n", R.rounds, np,
            (unsigned long long)withdiff, nover, (unsigned long long)late);
    free(skA); free(skB); free(nA); free(nB); free(loA); free(hiA); free(loB); free(hiB);
    free(ok); free(eB); free(fa); free(fb); free(nfa); free(nfb);
    if (!nover) { np = 0; break; }
    int narrow = 0;
    for (size_t i = 0; i < nover; i++) if (pending[i].hi - pending[i].lo < 2) narrow = 1;
    if (R.rounds >= max_rounds || narrow) { R.unresolved = (int)nover; break; }
    R.splits += nover;
    range* next = malloc(2 * nover * sizeof(range));
    for (size_t i = 0; i < nover; i++) {
      u64 mid = pending[i].lo + (pending[i].hi - pending[i].lo) / 2;
      next[2 * i].lo = pending[i].lo; next[2 * i].hi = mid;
      next[2 * i + 1].lo = mid; next[2 * i + 1].hi = pending[i].hi;
    }
    free(pending);
    pending = next;
    np = 2 * nover;
  }
  free(pending);
  qsort(oa.v, oa.n, 8, cmp_u64);
  qsort(ob.v, ob.n, 8, cmp_u64);
  R.amb = oa.v; R.namb = oa.n;
  R.bma = ob.v; R.nbma = ob.n;
  return R;
}

/* ---------- back to k-mers ---------- */
typedef struct { u64* t; size_t mask; } hset;
static inline u64 mixh(u64 x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; return x; }
static hset hs_make(const u64* v, size_t n)
{
  size_t cap = 16;
  while (cap < 2 * n + 16) cap <<= 1;
  hset h = {malloc(cap * 8), cap - 1};
  for (size_t i = 0; i < cap; i++) h.t[i] = ~0ULL;
  for (size_t i = 0; i < n; i++) {
    size_t j = mixh(v[i]) & h.mask;
    while (h.t[j] != ~0ULL && h.t[j] != v[i]) j = (j + 1) & h.mask;
    h.t[j] = v[i];
  }
  return h;
}
static inline int hs_has(const hset* h, u64 x)
{
  size_t j = mixh(x) & h->mask;
  while (h->t[j] != ~0ULL) { if (h->t[j] == x) return 1; j = (j + 1) & h->mask; }
  return 0;
}
/* canonical k-mer codes, k <= 64: base t at bits 2(k-1-t) of a 128-bit word (A0 C1 G2 T3), min of both strands */
typedef struct { u128* v; size_t n, cap; } vec2;
static void push2(vec2* a, u128 x)
{
  if (a->n == a->cap) {
    a->cap = a->cap ? a->cap * 2 : 1 << 16;
    a->v = realloc(a->v, a->cap * sizeof(u128));
    if (!a->v) { perror("realloc"); exit(1); }
  }
  a->v[a->n++] = x;
}
static int cmp_u128(const void* a, const void* b) { u128 x = *(const u128*)a, y = *(const u128*)b; return x < y ? -1 : x > y; }
static size_t uniq2(u128* v, size_t n)
{
  if (n) qsort(v, n, sizeof(u128), cmp_u128);
  size_t m = 0;
  for (size_t i = 0; i < n; i++) if (!m || v[i] != v[m - 1]) v[m++] = v[i];
  return m;
}
typedef struct { hset h; const char* seq; unsigned k; vec2 out; u64 hits; } rctx;
static int CODE[256];
static void emit_find(void* ctx, size_t start, const u64* hv, size_t n)
{
  rctx* r = ctx;
  for (size_t i = 0; i < n; i++) {
    if (!hs_has(&r->h, hv[i])) continue;
    r->hits++;
    const char* s = r->seq + start + i;
    u128 f = 0, rc = 0;
    for (unsigned t = 0; t < r->k; t++) {
      u128 b = (u128)CODE[(unsigned char)s[t]];
      f = (f << 2) | b;
      rc |= (u128)(3 - b) << (2 * t);
    }
    push2(&r->out, f < rc ? f : rc);
  }
}
/* canonical codes of the k-mers of g whose S is in vals (sorted unique); *hits = k-mer windows matched */
static vec2 recover(const kp7e_params* p, const genome* g, const u64* vals, size_t n, u64* hits)
{
  rctx r = {hs_make(vals, n), g->s, p->k, {0}, 0};
  u64* hb = malloc(4096 * 8);
  kp7e_sym_seq(p, g->s, g->n, 2, hb, 4096, emit_find, &r);
  free(hb);
  free(r.h.t);
  r.out.n = uniq2(r.out.v, r.out.n);
  *hits = r.hits;
  return r.out;
}

/* truths: native uint64 codes (recon/kmtools.c, k <= 32) or KMC `transform dump` text (KMER<tab>count, any k <= 64) */
static u128* read_u64_as128(const char* path, size_t* n)
{
  FILE* f = fopen(path, "rb");
  if (!f) { perror(path); exit(1); }
  struct stat st;
  fstat(fileno(f), &st);
  *n = (size_t)st.st_size / 8;
  u128* v = malloc(*n * sizeof(u128) + 16);
  u64 x;
  for (size_t i = 0; i < *n; i++) {
    if (fread(&x, 8, 1, f) != 1) { perror(path); exit(1); }
    v[i] = x;
  }
  fclose(f);
  *n = uniq2(v, *n);
  return v;
}
static u128* read_kmc_dump(const char* path, unsigned k, size_t* n)
{
  FILE* f = fopen(path, "r");
  if (!f) { perror(path); exit(1); }
  vec2 a = {0};
  char* line = NULL;
  size_t cap = 0;
  while (getline(&line, &cap, f) >= 0) {
    u128 x = 0;
    unsigned t = 0;
    for (; t < k && line[t] && line[t] != '\t'; t++) x = (x << 2) | (u128)CODE[(unsigned char)line[t]];
    if (t != k) { fprintf(stderr, "%s: bad line %s", path, line); exit(1); }
    push2(&a, x);
  }
  free(line);
  fclose(f);
  *n = uniq2(a.v, a.n);
  return a.v;
}
static size_t minus(const u128* a, size_t na, const u128* b, size_t nb)
{
  size_t i = 0, j = 0, m = 0;
  for (; i < na; i++) {
    while (j < nb && b[j] < a[i]) j++;
    if (!(j < nb && b[j] == a[i])) m++;
  }
  return m;
}

static int parse_seed(const char* s, u64* out)
{
  if (!strcmp(s, "os")) return 0;
  errno = 0;
  char* end;
  *out = strtoull(s, &end, 0);
  if (errno || *end) { fprintf(stderr, "bad seed %s\n", s); exit(2); }
  return 1;
}

int main(int argc, char** argv)
{
  const char *pa = NULL, *pb = NULL, *tamb = NULL, *tbma = NULL, *kamb = NULL, *kbma = NULL, *label = "", *dests = NULL;
  const char *sx = "os", *sl = "os", *rule = "pow2";
  unsigned k = 31;
  int m = 64, c = 2, max_rounds = 16, threads = 1, oneway = 0;
  double load = 0.8;
  u64 force_n0 = 0;
  for (int i = 1; i + 1 < argc; i += 2) {
    const char *o = argv[i], *v = argv[i + 1];
    if (!strcmp(o, "--a")) pa = v;
    else if (!strcmp(o, "--b")) pb = v;
    else if (!strcmp(o, "--k")) k = (unsigned)atoi(v);
    else if (!strcmp(o, "--m")) m = atoi(v);
    else if (!strcmp(o, "--c")) c = atoi(v);
    else if (!strcmp(o, "--dest")) dests = v;
    else if (!strcmp(o, "--seed-x")) sx = v;
    else if (!strcmp(o, "--seed-lam")) sl = v;
    else if (!strcmp(o, "--max-rounds")) max_rounds = atoi(v);
    else if (!strcmp(o, "--threads")) threads = atoi(v);
    else if (!strcmp(o, "--truth-amb")) tamb = v;
    else if (!strcmp(o, "--truth-bma")) tbma = v;
    else if (!strcmp(o, "--truth-kmc-amb")) kamb = v;
    else if (!strcmp(o, "--truth-kmc-bma")) kbma = v;
    else if (!strcmp(o, "--label")) label = v;
    else if (!strcmp(o, "--rule")) rule = v;
    else if (!strcmp(o, "--load")) load = atof(v);
    else if (!strcmp(o, "--buckets")) force_n0 = strtoull(v, NULL, 0);
    else if (!strcmp(o, "--oneway")) oneway = atoi(v);
    else { fprintf(stderr, "unknown option %s\n", o); return 2; }
  }
  int exact_rule = !strcmp(rule, "exact");
  if (!pa || !pb || !dests || k < 2 || k > 64 || m < 2 || m + c > MAXM || c < 1 || load <= 0 ||
      (!exact_rule && strcmp(rule, "pow2"))) {
    fprintf(stderr, "usage: kp10c --a A.fa --b B.fa --dest D[,D...] [--k 31 --m 64 --c 2 --seed-x X|os --seed-lam L|os "
                    "--max-rounds R --threads T --truth-amb F --truth-bma F | --truth-kmc-amb F --truth-kmc-bma F "
                    "--rule pow2|exact --load 0.8 --buckets N --oneway 0|1 --label S]\n");
    return 2;
  }
  omp_set_num_threads(threads);
  for (int j = 1; j <= MAXM + 1; j++) INV[j] = invq((u64)j);
  for (int i = 0; i < 256; i++) CODE[i] = 0;
  CODE['C'] = CODE['c'] = 1; CODE['G'] = CODE['g'] = 2; CODE['T'] = CODE['t'] = 3;
  u64 vx, vl;
  int hx = parse_seed(sx, &vx), hl = parse_seed(sl, &vl);
  kp7e_params p;
  kp7e_init_seeds(&p, k, hx ? &vx : NULL, hl ? &vl : NULL);

  side_t A = {0}, B = {0};
  build_side(&A, &p, pa);
  fprintf(stderr, "A: %" PRIu64 " k-mers, %zu distinct values; hash %.1f s, dedup %.1f s\n", A.kmers, A.n, A.t_hash.w, A.t_dedup.w);
  build_side(&B, &p, pb);
  fprintf(stderr, "B: %" PRIu64 " k-mers, %zu distinct values; hash %.1f s, dedup %.1f s\n", B.kmers, B.n, B.t_hash.w, B.t_dedup.w);
  size_t namb = 0, nbma = 0;
  u128 *Tamb = NULL, *Tbma = NULL;
  if (tamb && tbma) { Tamb = read_u64_as128(tamb, &namb); Tbma = read_u64_as128(tbma, &nbma); }
  else if (kamb && kbma) { Tamb = read_kmc_dump(kamb, k, &namb); Tbma = read_kmc_dump(kbma, k, &nbma); }
  /* one-way: the party with the larger set decodes (B decodes when |B| > |A|) */
  int swap = oneway && B.n > A.n;
  const side_t *D0 = swap ? &B : &A, *D1 = swap ? &A : &B;
  const unsigned M = (unsigned)(m + c), bucket_bits = M * 61 + 32, kmer_bytes = k <= 32 ? 8 : 16;

  char* dl = strdup(dests);
  for (char* tok = strtok(dl, ","); tok; tok = strtok(NULL, ",")) {
    double dest = atof(tok);
    u64 n0, tok_n0 = 0;
    char* at = strchr(tok, '@');
    if (at) tok_n0 = strtoull(at + 1, NULL, 0); /* "D@N": this dest with N first-round buckets */
    if (tok_n0) n0 = tok_n0;
    else if (force_n0) n0 = force_n0;
    else if (exact_rule) { n0 = (u64)ceil(dest / (load * m)); if (n0 < 1) n0 = 1; }
    else { unsigned bits = 0; while (bits < 40 && dest / (double)(1ULL << bits) > load * m) bits++; n0 = 1ULL << bits; }
    fprintf(stderr, "dest %.0f -> %" PRIu64 " buckets (%s)\n", dest, n0, tok_n0 || force_n0 ? "forced" : rule);
    tick t0 = now();
    result R = reconcile(D0, D1, n0, m, c, max_rounds, oneway);
    tick tr = since(t0);
    if (swap) { u64* t = R.amb; R.amb = R.bma; R.bma = t; size_t tn = R.namb; R.namb = R.nbma; R.nbma = tn; }
    t0 = now();
    u64 hitsA = 0, hitsB = 0;
    vec2 ka = recover(&p, &A.g, R.amb, R.namb, &hitsA);
    vec2 kb = recover(&p, &B.g, R.bma, R.nbma, &hitsB);
    tick tk = since(t0);
    long miss = -1, spur = -1;
    if (Tamb && Tbma)
      miss = (long)(minus(Tamb, namb, ka.v, ka.n) + minus(Tbma, nbma, kb.v, kb.n)),
      spur = (long)(minus(ka.v, ka.n, Tamb, namb) + minus(kb.v, kb.n, Tbma, nbma));
    u64 bytes = R.entries * 8, nbuckets_sent = R.entries / (M + 1);
    double packed_sketch = (double)nbuckets_sent * bucket_bits / 8.0;
    /* the sender's own differences, returned to the decoder as k-mers in the one-way protocol */
    size_t nsend = swap ? ka.n : kb.n;
    double ow_stored = (double)bytes + (double)R.fail_idx_bytes + (double)R.q_bytes + (double)nsend * kmer_bytes + 8;
    double ow_packed = packed_sketch + R.fail_idx_bits / 8.0 + R.q_bits / 8.0 + (double)nsend * 2 * k / 8.0 + 8;
    size_t drec = R.namb + R.nbma;
    double den = (namb + nbma) ? (double)(namb + nbma) : (double)(drec ? drec : 1);
    char rs[2048] = "", rf[2048] = "";
    for (int r = 0; r < R.rounds && r < MAXR; r++) {
      snprintf(rs + strlen(rs), sizeof rs - strlen(rs), "%s%" PRIu64, r ? ", " : "", R.sent[r]);
      snprintf(rf + strlen(rf), sizeof rf - strlen(rf), "%s%" PRIu64, r ? ", " : "", R.failed[r]);
    }
    printf("{\"tool\": \"kp10c\", \"label\": \"%s\", \"A\": \"%s\", \"B\": \"%s\", \"k\": %u, \"m\": %d, \"c\": %d, "
           "\"x\": %" PRIu64 ", \"lam\": %" PRIu64 ", \"seed_x\": \"%s\", \"seed_lam\": \"%s\", \"threads\": %d, "
           "\"kmers_A\": %" PRIu64 ", \"kmers_B\": %" PRIu64 ", \"values_A\": %zu, \"values_B\": %zu, "
           "\"dest\": %.0f, \"rule\": \"%s\", \"load\": %.4f, \"buckets0\": %" PRIu64 ", \"forced_buckets\": %s, "
           "\"oneway\": %d, \"decoder\": \"%s\", "
           "\"bits0\": %d, \"rounds\": %d, \"splits\": %" PRIu64 ", \"unresolved_buckets\": %d, "
           "\"buckets_sent_per_round\": [%s], \"buckets_failed_per_round\": [%s], \"late_failures\": %" PRIu64 ", "
           "\"buckets_with_difference_round1\": %" PRIu64 ", \"sketch_entries\": %" PRIu64 ", \"sketch_bytes\": %" PRIu64 ", "
           "\"bytes_round1\": %" PRIu64 ", \"sketch_bytes_packed\": %.0f, "
           "\"twoway_bytes\": %" PRIu64 ", \"twoway_bytes_packed\": %.0f, "
           "\"oneway_fail_index_bytes\": %" PRIu64 ", \"oneway_fail_index_bits\": %" PRIu64 ", \"oneway_q_buckets\": %" PRIu64
           ", \"oneway_q_coeffs\": %" PRIu64 ", \"oneway_reply_kmers\": %zu, \"oneway_bytes\": %.0f, \"oneway_bytes_packed\": %.0f, "
           "\"values_recovered_A_minus_B\": %zu, \"values_recovered_B_minus_A\": %zu, "
           "\"kmers_recovered_A_minus_B\": %zu, \"kmers_recovered_B_minus_A\": %zu, \"windows_matched_A\": %" PRIu64 ", "
           "\"windows_matched_B\": %" PRIu64 ", \"truth_A_minus_B\": %zu, \"truth_B_minus_A\": %zu, \"missing\": %ld, "
           "\"spurious\": %ld, \"exact\": %s, \"bytes_per_true_difference\": %.3f, "
           "\"packed_bytes_per_true_difference\": %.3f, \"twoway_bytes_per_true_difference\": %.3f, "
           "\"twoway_packed_bytes_per_true_difference\": %.3f, \"oneway_bytes_per_true_difference\": %.3f, "
           "\"oneway_packed_bytes_per_true_difference\": %.3f, "
           "\"wall_hash_A\": %.2f, \"cpu_hash_A\": %.2f, \"wall_dedup_A\": %.2f, \"cpu_dedup_A\": %.2f, "
           "\"wall_hash_B\": %.2f, \"cpu_hash_B\": %.2f, \"wall_dedup_B\": %.2f, \"cpu_dedup_B\": %.2f, "
           "\"wall_sketch_A\": %.2f, \"cpu_sketch_A\": %.2f, \"wall_sketch_B\": %.2f, \"cpu_sketch_B\": %.2f, "
           "\"wall_decode\": %.2f, \"cpu_decode\": %.2f, \"wall_reconcile\": %.2f, \"cpu_reconcile\": %.2f, "
           "\"wall_recover_kmers\": %.2f, \"cpu_recover_kmers\": %.2f}\n",
           label, pa, pb, k, m, c, p.x, p.lam, sx, sl, threads, A.kmers, B.kmers, A.n, B.n, dest, rule, load, n0,
           tok_n0 || force_n0 ? "true" : "false", oneway, oneway ? (swap ? "B" : "A") : "both", R.bits0, R.rounds,
           R.splits, R.unresolved, rs, rf, R.late_failures, R.buckets_with_difference, R.entries, bytes,
           n0 * (u64)(M + 1) * 8, packed_sketch, 2 * bytes, 2 * packed_sketch, R.fail_idx_bytes, R.fail_idx_bits,
           R.q_buckets, R.q_coeffs, oneway ? nsend : (size_t)0, oneway ? ow_stored : 0.0, oneway ? ow_packed : 0.0,
           R.namb, R.nbma, ka.n, kb.n, hitsA, hitsB, namb, nbma, miss, spur,
           (Tamb && miss == 0 && spur == 0 && !R.unresolved) ? "true" : "false",
           (double)bytes / den, packed_sketch / den, 2.0 * bytes / den, 2.0 * packed_sketch / den,
           oneway ? ow_stored / den : 0.0, oneway ? ow_packed / den : 0.0,
           A.t_hash.w, A.t_hash.c, A.t_dedup.w, A.t_dedup.c, B.t_hash.w, B.t_hash.c, B.t_dedup.w, B.t_dedup.c,
           R.t_sketchA.w, R.t_sketchA.c, R.t_sketchB.w, R.t_sketchB.c, R.t_decode.w, R.t_decode.c, tr.w, tr.c, tk.w, tk.c);
    fflush(stdout);
    free(R.amb); free(R.bma); free(ka.v); free(kb.v);
  }
  return 0;
}
