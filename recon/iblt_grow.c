/* km-peeler-compatible IBLT protocols on a known k-mer difference (round-3 owed runs 2, 3 and 5;
 * results/round3_owed_2026-11-12.md). k <= 32, canonical 2-bit codes as recon/kmtools.c writes them.
 *
 *   iblt_grow proto  --amb amb.u64 --bma bma.u64 --proto fixed|resend|grow|subtract --f F[,F...] --reps R
 *                    [--threads T] [--seed S|os] [--label L] [--k 31] [--max-rounds 8]
 *   iblt_grow write  --in set.u64 --n N --seed S --out t.iblt [--k 31]   one km-peeler table file (validation)
 *   iblt_grow stat   --in t.iblt                                         max |counter|, nonzero cells of a .iblt
 *   iblt_grow peel   --in d.iblt                                         this file's peeler on a km-peeler table
 *
 * THE TABLE IS km-peeler's (github.com/yhhshb/km-peeler b8ec39b, ibflib.c), reimplemented and checked byte for byte
 * against `ibltseq build` (the `write` subcommand; recon/round3_owed_job.sh, phase iblt_validate):
 *   r = 3 subtables of chunk = ceil((1.222f + 0.0f) * n / r + 1) cells (float arithmetic, as ibf_sketch_init);
 *   key = the k-mer 2-bit packed, base t in byte t/4 at bits 2(3 - t%4) (pack2bit), 8 bytes for k = 31;
 *   subtable j (seed root + j) holds the key at MurmurHash3_x64_128(key, 8, root + j).lo % chunk;
 *   a cell is an int64 counter and the XOR of its keys (16 bytes; file = 14-byte header, the last byte the key
 *   length k, + cells + 4 B per seed);
 *   a cell is pure when its counter is +-1 and the key it holds hashes back to it; peeling follows ibf_list_seq
 *   step for step, and listing succeeds when every counter is 0 (find_peelable_bucket's test).
 *
 * LINEARITY. T(A) - T(B) = T(A\B) - T(B\A), cell for cell (counters add, keys XOR), so the difference table is built
 * from the difference sets alone. That is exact, not an approximation, and recon/round3_owed_job.sh checks it on
 * whole genomes (km-peeler's own d.iblt from the full sets vs from the difference sets, cmp). The truth therefore
 * BUILDS the tables here; no protocol DECISION reads it: each one acts on the listing's own exit and residual.
 *
 * PROTOCOLS. n0 = round(1.05 f D), f in --f; a round's message is one party's table (or its new cells).
 *   fixed     one table at n = round(f D) (here f is the ratio n/D); success rate, listed fraction, max |counter|.
 *   resend    round 2's protocol: on a failed listing resend a fresh table (new seeds) at double n, all cells.
 *   grow      Mitzenmacher-Pagh style: subtable sizes chunk0 * 2^t with the SAME seeds, so the table of round t
 *             folds onto round t-1's (cell i of the half = cells i and i + half of the whole). Round t sends only
 *             the new upper half of each subtable; the receiver rebuilds the lower half by subtraction.
 *   subtract  on a failed listing the receiver keeps what it listed (it reinserts listed A\B and removes listed
 *             B\A from its own copy), and replies with an estimate of what is left, R = ceil(sum |counter| / 3)
 *             over the residual (exact when the residual is one-sided; 8 bytes). The sender sends a fresh table
 *             (new seeds) at n = round(1.05 R), or at least twice the last n when that table listed under half its
 *             n (it was overloaded, and sum |counter| undercounts a two-sided residual); repeated until a listing
 *             succeeds.
 * Bytes: stored = km-peeler's file bytes per table message (14 + 16 cells + 12; a grow increment = 16 B per new
 * cell), plus 8 B per subtract reply. Packed = (2k + b) bits per cell, where b is the counter width needed to carry
 * the DIFFERENCE counters modulo 2^b (b = 1 + ceil(log2(max |counter| + 1)), max over the run's tables); headers
 * and seeds are omitted from packed. Each rep draws its root seeds from a splitmix64 stream seeded by --seed (os =
 * /dev/urandom) and writes them; distinct reps never share a hash function (checked).
 */
#define _GNU_SOURCE
#include <inttypes.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "murmur3.h" /* km-peeler's copy, compiled from the pinned checkout */

typedef uint64_t u64;
#define R3 3
#define MAXROUNDS 16

static unsigned K = 31;

/* ---------- km-peeler key and hash ---------- */
static inline u64 packkey(u64 code)
{
  uint8_t b[8] = {0};
  for (unsigned t = 0; t < K; t++) {
    unsigned c = (unsigned)(code >> (2 * (K - 1 - t))) & 3;
    b[t / 4] |= (uint8_t)(c << (2 * (3 - (t % 4))));
  }
  u64 x;
  memcpy(&x, b, 8);
  return x;
}
static inline u64 unpackkey(u64 key)
{
  uint8_t b[8];
  memcpy(b, &key, 8);
  u64 code = 0;
  for (unsigned t = 0; t < K; t++) code = (code << 2) | ((b[t / 4] >> (2 * (3 - (t % 4)))) & 3);
  return code;
}
static inline u64 hpos(u64 key, uint32_t seed)
{
  u64 out[2];
  MurmurHash3_x64_128(&key, 8, seed, out);
  return out[0];
}
static u64 chunk_for(unsigned n)
{
  float ck = 1.222f, eps = 0.0f;
  unsigned char r = R3;
  return (unsigned long)ceil((ck + eps) * n / r + 1);
}

typedef struct { int64_t c; u64 k; } cell;
typedef struct { cell* t; u64 chunk; uint32_t seed; } table; /* 3 * chunk cells, subtable j seeded seed + j */

static void tinit(table* T, u64 chunk, uint32_t seed)
{
  T->chunk = chunk;
  T->seed = seed;
  T->t = calloc(R3 * chunk, sizeof(cell));
  if (!T->t) { perror("calloc"); exit(1); }
}
static inline void tadd(table* T, u64 key, int64_t s)
{
  for (unsigned j = 0; j < R3; j++) {
    u64 p = hpos(key, T->seed + j) % T->chunk + j * T->chunk;
    T->t[p].c += s;
    T->t[p].k ^= key;
  }
}
/* the half-size table: cell i of subtable j = cells i and i + chunk/2 */
static table fold(const table* T)
{
  table H;
  tinit(&H, T->chunk / 2, T->seed);
  for (unsigned j = 0; j < R3; j++)
    for (u64 i = 0; i < T->chunk; i++) {
      cell* d = &H.t[j * H.chunk + i % H.chunk];
      d->c += T->t[j * T->chunk + i].c;
      d->k ^= T->t[j * T->chunk + i].k;
    }
  return H;
}

/* ---------- peeling, km-peeler's purity test ---------- */
typedef struct { u64* v; size_t n, cap; } vec;
static void push(vec* a, u64 x)
{
  if (a->n == a->cap) {
    a->cap = a->cap ? 2 * a->cap : 1024;
    a->v = realloc(a->v, a->cap * 8);
    if (!a->v) { perror("realloc"); exit(1); }
  }
  a->v[a->n++] = x;
}
typedef struct {
  int success;           /* every counter 0 */
  u64 residual_cells;    /* nonzero cells left */
  u64 residual_abs;      /* sum |counter| left */
  int64_t max_abs;       /* max |counter| BEFORE peeling */
  int capped;            /* the peel-operation cap was hit */
} peelres;
/* km-peeler's find_peelable_bucket: the next cell with counter +-1 at or after *last, circularly; returns 1 if one was
 * found (*last = it) or if every counter is 0 (*last = blen), 0 if nonzero counters remain but none is +-1 */
static int find_peelable(const cell* t, u64 blen, u64* last)
{
  u64 start = *last;
  int empty = 1;
  for (; *last < blen; ++*last) {
    if (t[*last].c == 1 || t[*last].c == -1) return 1;
    if (t[*last].c) empty = 0;
  }
  for (*last = 0; *last < start; ++*last) {
    if (t[*last].c == 1 || t[*last].c == -1) return 1;
    if (t[*last].c) empty = 0;
  }
  *last = blen;
  return empty;
}
/* destructive; ibf_list_seq step for step (scan order, the skip of an impure +-1 cell, the peeled cell's counter
 * zeroed and its key left, MAXPASSES = 10). Listed keys with counter +1 (first table = A side) to la, -1 to lb. */
static peelres peel(table* T, vec* la, vec* lb)
{
  peelres P = {0};
  u64 blen = R3 * T->chunk, idx = 0, seen = 0;
  for (u64 i = 0; i < blen; i++) {
    int64_t c = T->t[i].c < 0 ? -T->t[i].c : T->t[i].c;
    if (c > P.max_abs) P.max_abs = c;
  }
  int peelable = find_peelable(T->t, blen, &idx);
  while (idx != blen && seen < 10 * blen) {
    u64 key = T->t[idx].k, pos[R3];
    int pure = 0;
    for (unsigned j = 0; j < R3; j++) {
      pos[j] = hpos(key, T->seed + j) % T->chunk + j * T->chunk;
      if (pos[j] == idx) pure = 1;
    }
    if (pure) {
      for (unsigned j = 0; j < R3; j++)
        if (pos[j] != idx) { T->t[pos[j]].c -= T->t[idx].c; T->t[pos[j]].k ^= key; }
      push(T->t[idx].c == 1 ? la : lb, key);
      T->t[idx].c = 0;
    } else {
      idx = (idx + 1) % blen;
    }
    seen++;
    peelable = find_peelable(T->t, blen, &idx);
  }
  P.capped = seen >= 10 * blen;
  P.success = peelable && !P.capped;
  for (u64 i = 0; i < blen; i++)
    if (T->t[i].c) {
      P.residual_cells++;
      P.residual_abs += (u64)(T->t[i].c < 0 ? -T->t[i].c : T->t[i].c);
    }
  return P;
}

/* ---------- sets ---------- */
static u64* read_u64(const char* path, size_t* n)
{
  FILE* f = fopen(path, "rb");
  if (!f) { perror(path); exit(1); }
  struct stat st;
  fstat(fileno(f), &st);
  *n = (size_t)st.st_size / 8;
  u64* v = malloc(*n * 8 + 8);
  if (*n && fread(v, 8, *n, f) != *n) { perror(path); exit(1); }
  fclose(f);
  return v;
}
static int has(const u64* v, size_t n, u64 x)
{
  size_t lo = 0, hi = n;
  while (lo < hi) { size_t m = (lo + hi) / 2; if (v[m] < x) lo = m + 1; else hi = m; }
  return lo < n && v[lo] == x;
}

static u64 splitmix(u64* s)
{
  u64 z = (*s += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static unsigned counter_bits(int64_t max_abs)
{
  unsigned b = 0;
  while (((int64_t)1 << b) < max_abs + 1) b++;
  return b + 1; /* sign */
}

/* ---------- one rep of one protocol ---------- */
typedef struct {
  const u64 *amb, *bma; /* truth codes, sorted */
  size_t na, nb;
  u64 *kamb, *kbma;     /* packed keys */
} diffset;

/* the net listing: every listed (key, sign), +1 = A side, summed per key; a key with net +1 is claimed in A\B, -1 in
 * B\A, 0 cancels (a key listed on both sides, or listed and later re-listed with the other sign by subtract), and
 * |net| > 1 is wrong. Compared element by element with the truth. */
typedef struct { u64 wrong, missing, claimed_a, claimed_b, cancelled; } score_t;
typedef struct { u64 code; int s; } signed_key;
static int cmp_sk(const void* a, const void* b)
{
  const signed_key *x = a, *y = b;
  return x->code < y->code ? -1 : x->code > y->code;
}
static score_t score(const diffset* S, const vec* la, const vec* lb)
{
  score_t r = {0};
  size_t n = la->n + lb->n, m = 0;
  signed_key* v = malloc((n + 1) * sizeof *v);
  for (size_t i = 0; i < la->n; i++) v[m++] = (signed_key){unpackkey(la->v[i]), +1};
  for (size_t i = 0; i < lb->n; i++) v[m++] = (signed_key){unpackkey(lb->v[i]), -1};
  qsort(v, m, sizeof *v, cmp_sk);
  for (size_t i = 0; i < m;) {
    size_t j = i;
    int net = 0;
    while (j < m && v[j].code == v[i].code) net += v[j++].s;
    if (net == 1) { r.claimed_a++; if (!has(S->amb, S->na, v[i].code)) r.wrong++; }
    else if (net == -1) { r.claimed_b++; if (!has(S->bma, S->nb, v[i].code)) r.wrong++; }
    else if (net == 0) r.cancelled++;
    else r.wrong++;
    i = j;
  }
  free(v);
  /* claimed sets are subsets of the truth when wrong == 0 */
  r.missing = (S->na + S->nb) - (r.claimed_a + r.claimed_b - r.wrong);
  return r;
}

/* the difference table T(A\B) - T(B\A); subtract also removes what the receiver has listed so far (reinserting
 * listed A-side keys into its own copy and deleting listed B-side keys), garbage keys included */
static void build_diff(table* T, const diffset* S, const vec* ca, const vec* cb)
{
  for (size_t i = 0; i < S->na; i++) tadd(T, S->kamb[i], +1);
  for (size_t i = 0; i < S->nb; i++) tadd(T, S->kbma[i], -1);
  if (ca) for (size_t i = 0; i < ca->n; i++) tadd(T, ca->v[i], -1);
  if (cb) for (size_t i = 0; i < cb->n; i++) tadd(T, cb->v[i], +1);
}

static double file_bytes(u64 chunk) { return 14.0 + 16.0 * R3 * chunk + 4.0 * R3; }

static void run_rep(const diffset* S, const char* proto, double f, int rep, u64 rootseed, int max_rounds,
                    const char* label, FILE* out)
{
  u64 D = S->na + S->nb, st = rootseed;
  int fixed = !strcmp(proto, "fixed"), grow = !strcmp(proto, "grow"), sub = !strcmp(proto, "subtract");
  unsigned n = (unsigned)floor((fixed ? f : 1.05 * f) * (double)D + 0.5);
  if (n < 1) n = 1;
  uint32_t seed = (uint32_t)splitmix(&st);
  int rounds = 0, success = 0;
  double bytes = 0, cells_sent = 0;
  int64_t maxabs = 0;
  char log[8192] = "";
  size_t L = 0;
  u64 allwrong = 0, allmissing = 0;
  vec ca = {0}, cb = {0}; /* subtract: everything the receiver has listed so far, by sign */
  table big = {0};
  u64 chunk0 = chunk_for(n);
  int G = 0;
  if (grow) { /* the largest table the run may need, built once; earlier rounds are its folds */
    G = max_rounds - 1;
    while (G > 0 && (chunk0 << G) > ((u64)1 << 34)) G--;
    tinit(&big, chunk0 << G, seed);
    build_diff(&big, S, NULL, NULL);
  }
  while (rounds < max_rounds) {
    rounds++;
    table T;
    double msg;
    if (grow) {
      int down = G - (rounds - 1);
      if (down < 0) break;
      T = big;
      table cur = {0};
      int owned = 0;
      for (int d = 0; d < down; d++) {
        table h = fold(owned ? &cur : &T);
        if (owned) free(cur.t);
        cur = h;
        owned = 1;
      }
      if (!owned) { tinit(&cur, big.chunk, seed); memcpy(cur.t, big.t, R3 * big.chunk * sizeof(cell)); }
      T = cur;
      msg = rounds == 1 ? file_bytes(T.chunk) : 16.0 * R3 * (T.chunk / 2);
      cells_sent += rounds == 1 ? (double)R3 * T.chunk : (double)R3 * (T.chunk / 2);
    } else {
      tinit(&T, chunk_for(n), seed);
      build_diff(&T, S, sub ? &ca : NULL, sub ? &cb : NULL);
      msg = file_bytes(T.chunk);
      cells_sent += (double)R3 * T.chunk;
    }
    bytes += msg;
    vec la = {0}, lb = {0};
    peelres P = peel(&T, &la, &lb);
    if (P.max_abs > maxabs) maxabs = P.max_abs;
    u64 listed = la.n + lb.n;
    if (sub) { for (size_t i = 0; i < la.n; i++) push(&ca, la.v[i]); for (size_t i = 0; i < lb.n; i++) push(&cb, lb.v[i]); }
    int exact = 0;
    if (P.success) {
      score_t sc = score(S, sub ? &ca : &la, sub ? &cb : &lb);
      exact = sc.wrong == 0 && sc.missing == 0;
      allwrong = sc.wrong;
      allmissing = sc.missing;
    }
    L += (size_t)snprintf(log + L, sizeof log - L,
                          "%s{\"n\": %u, \"chunk\": %" PRIu64 ", \"seed\": %u, \"msg_bytes\": %.0f, \"listed\": %" PRIu64
                          ", \"ok\": %d, \"residual_cells\": %" PRIu64 ", \"residual_abs\": %" PRIu64 ", \"max_abs_counter\": %"
                          PRId64 "}", rounds > 1 ? ", " : "", grow ? 0u : n, T.chunk, seed, msg, listed, P.success,
                          P.residual_cells, P.residual_abs, P.max_abs);
    free(T.t);
    free(la.v);
    free(lb.v);
    if (P.success) { success = exact ? 1 : -1; break; }
    if (fixed) break;
    if (grow) continue;
    if (sub) {
      u64 R = (P.residual_abs + 2) / 3;
      bytes += 8; /* the receiver's reply: its estimate of what is left */
      unsigned nr = (unsigned)floor(1.05 * (double)R + 0.5);
      /* sum |counter| / 3 undercounts when A- and B-side keys share cells; a table that listed under half its n was
       * overloaded, so its successor is also at least twice its size */
      if (listed < n / 2 && nr < 2 * n) nr = 2 * n;
      n = nr < 1 ? 1 : nr;
    } else {
      n *= 2;
    }
    seed = (uint32_t)splitmix(&st);
  }
  if (grow) free(big.t);
  free(ca.v);
  free(cb.v);
  unsigned b = counter_bits(maxabs);
  double packed = cells_sent * (2.0 * K + b) / 8.0 + (sub ? 8.0 * (rounds - 1) : 0);
  fprintf(out,
          "{\"label\": \"%s\", \"proto\": \"%s\", \"f\": %g, \"rep\": %d, \"D\": %" PRIu64 ", \"A_minus_B\": %zu, "
          "\"B_minus_A\": %zu, \"rootseed\": %" PRIu64 ", \"rounds\": %d, \"success\": %s, \"exact\": %s, "
          "\"undetected_wrong\": %s, \"wrong_keys\": %" PRIu64 ", \"missing_keys\": %" PRIu64 ", \"bytes\": %.0f, \"bytes_per_diff\": %.4f, \"cells_sent\": %.0f, \"max_abs_counter\": %" PRId64
          ", \"counter_bits\": %u, \"packed_bytes\": %.0f, \"packed_bytes_per_diff\": %.4f, \"per_round\": [%s]}\n",
          label, proto, f, rep, D, S->na, S->nb, rootseed, rounds, success ? "true" : "false",
          success == 1 ? "true" : "false", success == -1 ? "true" : "false", allwrong, allmissing, bytes,
          bytes / (double)D, cells_sent, maxabs, b, packed, packed / (double)D, log);
}

/* ---------- km-peeler file I/O (validation) ---------- */
static void put_be64(FILE* f, u64 x) { for (int i = 7; i >= 0; i--) fputc((int)(x >> (8 * i)) & 0xff, f); }
static void put_be32(FILE* f, uint32_t x) { for (int i = 3; i >= 0; i--) fputc((int)(x >> (8 * i)) & 0xff, f); }
static int cmd_write(const char* in, unsigned n, uint32_t seed, const char* outp)
{
  size_t ns;
  u64* v = read_u64(in, &ns);
  table T;
  tinit(&T, chunk_for(n), seed);
  for (size_t i = 0; i < ns; i++) tadd(&T, packkey(v[i]), +1);
  FILE* f = fopen(outp, "wb");
  if (!f) { perror(outp); return 1; }
  fputc(R3, f);
  float eps = 0.0f;
  uint32_t e;
  memcpy(&e, &eps, 4);
  put_be32(f, e);
  put_be64(f, T.chunk);
  fputc((int)K, f); /* GLEN: the key length lives in the header */
  for (u64 i = 0; i < R3 * T.chunk; i++) {
    put_be64(f, (u64)T.t[i].c);
    fwrite(&T.t[i].k, 8, 1, f);
  }
  for (unsigned j = 0; j < R3; j++) put_be32(f, seed + j);
  if (fclose(f)) { perror(outp); return 1; }
  fprintf(stderr, "{\"cmd\": \"write\", \"elements\": %zu, \"n\": %u, \"chunk\": %" PRIu64 ", \"seed\": %u}\n", ns, n,
          T.chunk, seed);
  return 0;
}
static int cmd_stat(const char* in)
{
  FILE* f = fopen(in, "rb");
  if (!f) { perror(in); return 1; }
  unsigned char hdr[14];
  if (fread(hdr, 1, 14, f) != 14) { perror(in); return 1; }
  u64 chunk = 0;
  for (int i = 0; i < 8; i++) chunk = (chunk << 8) | hdr[5 + i];
  int64_t mx = 0;
  u64 nz = 0, cells = hdr[0] * chunk;
  unsigned char c[16];
  for (u64 i = 0; i < cells; i++) {
    if (fread(c, 1, 16, f) != 16) { perror(in); return 1; }
    u64 x = 0;
    for (int b = 0; b < 8; b++) x = (x << 8) | c[b];
    int64_t v = (int64_t)x, a = v < 0 ? -v : v;
    if (a > mx) mx = a;
    if (v) nz++;
  }
  fclose(f);
  printf("{\"file\": \"%s\", \"r\": %u, \"chunk\": %" PRIu64 ", \"max_abs_counter\": %" PRId64 ", \"counter_bits\": %u, "
         "\"nonzero_cells\": %" PRIu64 "}\n", in, hdr[0], chunk, mx, counter_bits(mx), nz);
  return 0;
}

/* load a km-peeler table file (r = 3, 8-byte keys) and peel it with peel(); prints success and listed counts */
static int cmd_peel(const char* in)
{
  FILE* f = fopen(in, "rb");
  if (!f) { perror(in); return 1; }
  unsigned char hdr[14], c[16], sb[4];
  if (fread(hdr, 1, 14, f) != 14 || hdr[0] != R3) { fprintf(stderr, "%s: not an r = 3 table\n", in); return 1; }
  u64 chunk = 0;
  for (int i = 0; i < 8; i++) chunk = (chunk << 8) | hdr[5 + i];
  table T;
  tinit(&T, chunk, 0);
  for (u64 i = 0; i < R3 * chunk; i++) {
    if (fread(c, 1, 16, f) != 16) { perror(in); return 1; }
    u64 x = 0;
    for (int b = 0; b < 8; b++) x = (x << 8) | c[b];
    T.t[i].c = (int64_t)x;
    memcpy(&T.t[i].k, c + 8, 8);
  }
  uint32_t seeds[R3];
  for (unsigned j = 0; j < R3; j++) {
    if (fread(sb, 1, 4, f) != 4) { perror(in); return 1; }
    seeds[j] = (uint32_t)sb[0] << 24 | (uint32_t)sb[1] << 16 | (uint32_t)sb[2] << 8 | sb[3];
  }
  fclose(f);
  if (seeds[1] != seeds[0] + 1 || seeds[2] != seeds[0] + 2) { fprintf(stderr, "%s: seeds not consecutive\n", in); return 1; }
  T.seed = seeds[0];
  vec la = {0}, lb = {0};
  peelres P = peel(&T, &la, &lb);
  printf("{\"file\": \"%s\", \"seed\": %u, \"chunk\": %" PRIu64 ", \"success\": %d, \"listed_i\": %zu, \"listed_j\": %zu, "
         "\"residual_cells\": %" PRIu64 ", \"max_abs_counter\": %" PRId64 ", \"capped\": %d}\n", in, T.seed, chunk, P.success, la.n,
         lb.n, P.residual_cells, P.max_abs, P.capped);
  return 0;
}

int main(int argc, char** argv)
{
  if (argc < 2) goto usage;
  const char *amb = NULL, *bma = NULL, *proto = "fixed", *fs = "1.05", *label = "", *seeds = "os", *in = NULL,
             *outp = NULL;
  int reps = 1, threads = 1, max_rounds = 8;
  if (argc > 4096) goto usage;
  unsigned n = 0;
  for (int i = 2; i + 1 < argc; i += 2) {
    const char *o = argv[i], *v = argv[i + 1];
    if (!strcmp(o, "--amb")) amb = v;
    else if (!strcmp(o, "--bma")) bma = v;
    else if (!strcmp(o, "--proto")) proto = v;
    else if (!strcmp(o, "--f")) fs = v;
    else if (!strcmp(o, "--reps")) reps = atoi(v);
    else if (!strcmp(o, "--threads")) threads = atoi(v);
    else if (!strcmp(o, "--seed")) seeds = v;
    else if (!strcmp(o, "--label")) label = v;
    else if (!strcmp(o, "--k")) K = (unsigned)atoi(v);
    else if (!strcmp(o, "--max-rounds")) max_rounds = atoi(v);
    else if (!strcmp(o, "--in")) in = v;
    else if (!strcmp(o, "--n")) n = (unsigned)strtoul(v, NULL, 0);
    else if (!strcmp(o, "--out")) outp = v;
    else { fprintf(stderr, "unknown option %s\n", o); return 2; }
  }
  if (K < 1 || K > 32 || max_rounds < 1 || max_rounds > MAXROUNDS) goto usage;
  if (!strcmp(argv[1], "stat") && in) return cmd_stat(in);
  if (!strcmp(argv[1], "peel") && in) return cmd_peel(in);
  if (!strcmp(argv[1], "write") && in && outp && n) return cmd_write(in, n, (uint32_t)strtoul(seeds, NULL, 0), outp);
  if (strcmp(argv[1], "proto") || !amb || !bma) goto usage;
  if (strcmp(proto, "fixed") && strcmp(proto, "resend") && strcmp(proto, "grow") && strcmp(proto, "subtract")) goto usage;
  diffset S;
  S.amb = read_u64(amb, &S.na);
  S.bma = read_u64(bma, &S.nb);
  S.kamb = malloc(S.na * 8 + 8);
  S.kbma = malloc(S.nb * 8 + 8);
  for (size_t i = 0; i < S.na; i++) S.kamb[i] = packkey(S.amb[i]);
  for (size_t i = 0; i < S.nb; i++) S.kbma[i] = packkey(S.bma[i]);
  u64 st;
  if (!strcmp(seeds, "os")) {
    FILE* r = fopen("/dev/urandom", "rb");
    if (!r || fread(&st, 8, 1, r) != 1) { perror("/dev/urandom"); return 1; }
    fclose(r);
  } else {
    st = strtoull(seeds, NULL, 0);
  }
  double fv[64];
  int nf = 0;
  char* dup = strdup(fs);
  for (char* t = strtok(dup, ","); t && nf < 64; t = strtok(NULL, ",")) fv[nf++] = atof(t);
  /* one root seed stream per (f, rep); the rep's own stream then gives each round's km-peeler root seed */
  if (reps < 1 || reps > 100000) goto usage;
  int tot = nf * reps;
  u64* roots = malloc((size_t)tot * 8);
  for (int i = 0; i < tot; i++) roots[i] = splitmix(&st);
  omp_set_num_threads(threads);
  char** outs = calloc((size_t)tot, sizeof(char*));
#pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < tot; i++) {
    size_t sz = 0;
    FILE* mf = open_memstream(&outs[i], &sz);
    run_rep(&S, proto, fv[i / reps], i % reps, roots[i], max_rounds, label, mf);
    fclose(mf);
  }
  for (int i = 0; i < tot; i++) { fputs(outs[i], stdout); free(outs[i]); }
  return 0;
usage:
  fprintf(stderr, "usage: iblt_grow proto --amb A --bma B --proto fixed|resend|grow|subtract --f F[,F] --reps R "
                  "[--threads T --seed S|os --label L --k 31 --max-rounds 8]\n"
                  "       iblt_grow write --in set.u64 --n N --seed S --out t.iblt [--k 31]\n"
                  "       iblt_grow stat --in t.iblt\n       iblt_grow peel --in d.iblt\n");
  return 2;
}
