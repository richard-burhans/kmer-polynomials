/* k-mer set tools for the KP-10 vs Shibuya comparison (recon/kp10_vs_iblt_job.sh). k <= 32, 2-bit codes A0 C1 G2 T3,
 * base t of the k-mer at bits 2(k-1-t), so numeric order of codes = lexicographic order of the ACGT strings.
 *
 *   kmtools kset K < in.fa > set.u64             distinct CANONICAL k-mers, sorted (min of forward and reverse-complement
 *                                                 code; case-blind; windows with a non-ACGT byte skipped; records do
 *                                                 not join). The same definition as recon/kmerset.py and KMC 3.
 *   kmtools text K < set.u64 > kmers.txt         one ACGT line per code (km-peeler's `ibltseq build` input)
 *   kmtools diff A.u64 B.u64 amb.u64 bma.u64     exact A\B and B\A of two sorted sets
 *   kmtools kmcdump K < dump.txt > set.u64       KMC `transform dump` lines (KMER<tab>count) -> sorted codes
 *   kmtools same X.u64 Y.u64                     are two sorted sets identical?
 *   kmtools listscore K amb.u64 bma.u64 < list   score km-peeler `ibltseq list` output ("i,KMER" = in the first IBLT
 *                                                 only = A\B, "j,KMER" = B\A) against the truth, raw and with k-mers
 *                                                 listed on both sides cancelled
 *
 * Every subcommand prints one JSON line on stderr (listscore and same: stdout). Inputs are native uint64. */
#define _GNU_SOURCE
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

typedef uint64_t u64;

static signed char CODE[256];
static void code_init(void)
{
  memset(CODE, -1, sizeof CODE);
  CODE['A'] = CODE['a'] = 0;
  CODE['C'] = CODE['c'] = 1;
  CODE['G'] = CODE['g'] = 2;
  CODE['T'] = CODE['t'] = 3;
}

typedef struct {
  u64* v;
  size_t n, cap;
} vec;
static void push(vec* a, u64 x)
{
  if (a->n == a->cap) {
    a->cap = a->cap ? a->cap * 2 : 1 << 20;
    a->v = realloc(a->v, a->cap * 8);
    if (!a->v) { perror("realloc"); exit(1); }
  }
  a->v[a->n++] = x;
}

/* LSD radix sort, 16-bit digits; then unique. Returns the new length. */
static size_t sort_unique(u64* a, size_t n)
{
  u64* b = malloc(n * 8 + 8);
  if (!b) { perror("malloc"); exit(1); }
  size_t* cnt = malloc(65536 * sizeof(size_t));
  for (int pass = 0; pass < 4; pass++) {
    int sh = 16 * pass;
    memset(cnt, 0, 65536 * sizeof(size_t));
    for (size_t i = 0; i < n; i++) cnt[(a[i] >> sh) & 0xFFFF]++;
    size_t s = 0;
    for (int d = 0; d < 65536; d++) { size_t c = cnt[d]; cnt[d] = s; s += c; }
    for (size_t i = 0; i < n; i++) b[cnt[(a[i] >> sh) & 0xFFFF]++] = a[i];
    u64* t = a; a = b; b = t;
  }
  /* 4 passes: the sorted data is back in the caller's array */
  free(b);
  free(cnt);
  size_t m = 0;
  for (size_t i = 0; i < n; i++)
    if (m == 0 || a[i] != a[m - 1]) a[m++] = a[i];
  return m;
}

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

static void write_u64(FILE* f, const u64* v, size_t n)
{
  if (n && fwrite(v, 8, n, f) != n) { perror("fwrite"); exit(1); }
  if (fflush(f)) { perror("fflush"); exit(1); }
}

static int cmd_kset(unsigned k)
{
  code_init();
  const u64 mask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);
  const int top = 2 * (int)(k - 1);
  vec all = {0};
  char* line = NULL;
  size_t lcap = 0;
  ssize_t l;
  u64 f = 0, r = 0, run = 0, records = 0, bases = 0, kmers = 0;
  while ((l = getline(&line, &lcap, stdin)) >= 0) {
    if (l && line[0] == '>') { records++; run = 0; continue; }
    for (ssize_t i = 0; i < l; i++) {
      int c = CODE[(unsigned char)line[i]];
      if (line[i] == '\n' || line[i] == '\r') continue;
      bases++;
      if (c < 0) { run = 0; continue; }
      f = ((f << 2) | (u64)c) & mask;
      r = (r >> 2) | ((u64)(3 - c) << top);
      if (++run >= k) { push(&all, f < r ? f : r); kmers++; }
    }
  }
  size_t n = sort_unique(all.v, all.n);
  write_u64(stdout, all.v, n);
  fprintf(stderr, "{\"cmd\": \"kset\", \"k\": %u, \"records\": %" PRIu64 ", \"bases\": %" PRIu64 ", \"kmers\": %" PRIu64
          ", \"distinct_canonical\": %zu}\n", k, records, bases, kmers, n);
  return 0;
}

static int cmd_text(unsigned k)
{
  static const char ACGT[4] = {'A', 'C', 'G', 'T'};
  u64 buf[8192];
  char out[8192 * 33];
  size_t got;
  u64 lines = 0;
  while ((got = fread(buf, 8, 8192, stdin)) > 0) {
    char* p = out;
    for (size_t i = 0; i < got; i++) {
      for (unsigned t = 0; t < k; t++) p[t] = ACGT[(buf[i] >> (2 * (k - 1 - t))) & 3];
      p[k] = '\n';
      p += k + 1;
    }
    if (fwrite(out, 1, (size_t)(p - out), stdout) != (size_t)(p - out)) { perror("fwrite"); return 1; }
    lines += got;
  }
  if (fflush(stdout)) { perror("fflush"); return 1; }
  fprintf(stderr, "{\"cmd\": \"text\", \"k\": %u, \"lines\": %" PRIu64 "}\n", k, lines);
  return 0;
}

/* a \ b of sorted unique arrays into out (may be NULL: count only) */
static size_t minus(const u64* a, size_t na, const u64* b, size_t nb, u64* out)
{
  size_t i = 0, j = 0, m = 0;
  while (i < na) {
    while (j < nb && b[j] < a[i]) j++;
    if (j < nb && b[j] == a[i]) { i++; continue; }
    if (out) out[m] = a[i];
    m++;
    i++;
  }
  return m;
}

static int cmd_diff(const char* pa, const char* pb, const char* pamb, const char* pbma)
{
  size_t na, nb;
  u64 *a = read_u64(pa, &na), *b = read_u64(pb, &nb);
  u64* o = malloc((na > nb ? na : nb) * 8 + 8);
  size_t amb = minus(a, na, b, nb, o);
  FILE* f = fopen(pamb, "wb");
  write_u64(f, o, amb);
  fclose(f);
  size_t bma = minus(b, nb, a, na, o);
  f = fopen(pbma, "wb");
  write_u64(f, o, bma);
  fclose(f);
  fprintf(stderr, "{\"cmd\": \"diff\", \"A\": \"%s\", \"B\": \"%s\", \"distinct_A\": %zu, \"distinct_B\": %zu, "
          "\"A_minus_B\": %zu, \"B_minus_A\": %zu, \"D\": %zu}\n", pa, pb, na, nb, amb, bma, amb + bma);
  return 0;
}

static int cmd_kmcdump(unsigned k)
{
  code_init();
  vec all = {0};
  char* line = NULL;
  size_t lcap = 0;
  ssize_t l;
  u64 bad = 0;
  while ((l = getline(&line, &lcap, stdin)) >= 0) {
    if (l < (ssize_t)k) { bad++; continue; }
    u64 c = 0;
    int ok = 1;
    for (unsigned t = 0; t < k; t++) {
      int v = CODE[(unsigned char)line[t]];
      if (v < 0) { ok = 0; break; }
      c = (c << 2) | (u64)v;
    }
    if (!ok || (l > (ssize_t)k && line[k] != '\t' && line[k] != ' ' && line[k] != '\n')) { bad++; continue; }
    push(&all, c);
  }
  size_t n = sort_unique(all.v, all.n);
  write_u64(stdout, all.v, n);
  fprintf(stderr, "{\"cmd\": \"kmcdump\", \"k\": %u, \"lines\": %zu, \"distinct\": %zu, \"unparsed\": %" PRIu64 "}\n",
          k, all.n, n, bad);
  return 0;
}

static int cmd_same(const char* px, const char* py)
{
  size_t nx, ny;
  u64 *x = read_u64(px, &nx), *y = read_u64(py, &ny);
  int same = nx == ny && (nx == 0 || !memcmp(x, y, nx * 8));
  printf("{\"cmd\": \"same\", \"X\": \"%s\", \"Y\": \"%s\", \"n_X\": %zu, \"n_Y\": %zu, \"X_minus_Y\": %zu, "
         "\"Y_minus_X\": %zu, \"identical\": %s}\n", px, py, nx, ny, minus(x, nx, y, ny, NULL), minus(y, ny, x, nx, NULL),
         same ? "true" : "false");
  return same ? 0 : 1;
}

static int cmd_listscore(unsigned k, const char* pamb, const char* pbma)
{
  code_init();
  size_t namb, nbma;
  u64 *amb = read_u64(pamb, &namb), *bma = read_u64(pbma, &nbma);
  vec gi = {0}, gj = {0};
  char* line = NULL;
  size_t lcap = 0;
  ssize_t l;
  u64 bad = 0;
  while ((l = getline(&line, &lcap, stdin)) >= 0) {
    if (l < (ssize_t)k + 2 || line[1] != ',' || (line[0] != 'i' && line[0] != 'j')) { bad++; continue; }
    u64 c = 0;
    int ok = 1;
    for (unsigned t = 0; t < k; t++) {
      int v = CODE[(unsigned char)line[2 + t]];
      if (v < 0) { ok = 0; break; }
      c = (c << 2) | (u64)v;
    }
    if (!ok) { bad++; continue; }
    push(line[0] == 'i' ? &gi : &gj, c);
  }
  size_t ni = sort_unique(gi.v, gi.n), nj = sort_unique(gj.v, gj.n);
  /* raw: km-peeler's listing as printed */
  size_t miss_a = minus(amb, namb, gi.v, ni, NULL), spur_a = minus(gi.v, ni, amb, namb, NULL);
  size_t miss_b = minus(bma, nbma, gj.v, nj, NULL), spur_b = minus(gj.v, nj, bma, nbma, NULL);
  int exact = !miss_a && !spur_a && !miss_b && !spur_b && !bad;
  /* cancelled: a k-mer listed on BOTH sides is a net-zero pair (a false "pure" bucket peeled, then undone: km-peeler
   * stores no check hash, only re-hashes a candidate key to its bucket) and cannot be in both A\B and B\A, so both copies are dropped */
  u64* ci = malloc(ni * 8 + 8);
  u64* cj = malloc(nj * 8 + 8);
  size_t nci = minus(gi.v, ni, gj.v, nj, ci), ncj = minus(gj.v, nj, gi.v, ni, cj), both = ni - nci;
  size_t cmiss = minus(amb, namb, ci, nci, NULL) + minus(bma, nbma, cj, ncj, NULL);
  size_t cspur = minus(ci, nci, amb, namb, NULL) + minus(cj, ncj, bma, nbma, NULL);
  int cexact = !cmiss && !cspur && !bad;
  printf("{\"cmd\": \"listscore\", \"listed_A_minus_B\": %zu, \"listed_B_minus_A\": %zu, \"truth_A_minus_B\": %zu, "
         "\"truth_B_minus_A\": %zu, \"missing\": %zu, \"spurious\": %zu, \"unparsed_lines\": %" PRIu64 ", \"exact\": %s, "
         "\"listed_both_sides\": %zu, \"cancelled_missing\": %zu, \"cancelled_spurious\": %zu, \"cancelled_exact\": %s}\n",
         ni, nj, namb, nbma, miss_a + miss_b, spur_a + spur_b, bad, exact ? "true" : "false", both, cmiss, cspur,
         cexact ? "true" : "false");
  return 0;
}

int main(int argc, char** argv)
{
  if (argc >= 3 && !strcmp(argv[1], "kset")) return cmd_kset((unsigned)atoi(argv[2]));
  if (argc >= 3 && !strcmp(argv[1], "text")) return cmd_text((unsigned)atoi(argv[2]));
  if (argc >= 6 && !strcmp(argv[1], "diff")) return cmd_diff(argv[2], argv[3], argv[4], argv[5]);
  if (argc >= 3 && !strcmp(argv[1], "kmcdump")) return cmd_kmcdump((unsigned)atoi(argv[2]));
  if (argc >= 4 && !strcmp(argv[1], "same")) return cmd_same(argv[2], argv[3]);
  if (argc >= 5 && !strcmp(argv[1], "listscore")) return cmd_listscore((unsigned)atoi(argv[2]), argv[3], argv[4]);
  fprintf(stderr, "usage: kmtools kset K | text K | diff A B AMB BMA | kmcdump K | same X Y | listscore K AMB BMA\n");
  return 2;
}
