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
 *   kmtools tally K diff.u64 < genome.fa         per record (round-3 run 8): name, length, valid k-mer windows, windows
 *                                                 whose canonical k-mer is in diff, distinct diff k-mers in the record,
 *                                                 and those first seen in it (records in file order); TSV on stdout
 *   kmtools extscore K amb.u64 bma.u64 miss.txt < list   km-peeler's extended-syncmer recovery (their Section 3.5,
 *                                                 as recon/kmpeeler_run.py ext): canonical k-mers of the "i," and "j,"
 *                                                 sequences, those in both dropped, scored against amb U bma; the
 *                                                 missed k-mers are written to miss.txt (one per line)
 *   kmtools locate K kmers.txt < genome.fa       every occurrence of the given k-mers (either strand): k-mer, record,
 *                                                 0-based position, strand, and the 20 bases on either side (lowercase
 *                                                 and non-ACGT as written); TSV on stdout
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

/* ---------- round 3 (2026-11-12) ---------- */
static long find(const u64* v, size_t n, u64 x)
{
  size_t lo = 0, hi = n;
  while (lo < hi) { size_t m = (lo + hi) / 2; if (v[m] < x) lo = m + 1; else hi = m; }
  return lo < n && v[lo] == x ? (long)lo : -1;
}

typedef struct { char* name; char* s; size_t n, cap; } record;
/* next FASTA record from stdin into r (sequence as written); 0 at EOF */
static int next_record(record* r, char** line, size_t* lcap, int* have_header)
{
  ssize_t l;
  r->n = 0;
  if (!*have_header) {
    while ((l = getline(line, lcap, stdin)) >= 0 && (*line)[0] != '>') {}
    if (l < 0) return 0;
  }
  free(r->name);
  r->name = strdup(*line + 1);
  r->name[strcspn(r->name, " \t\r\n")] = 0;
  *have_header = 0;
  while ((l = getline(line, lcap, stdin)) >= 0) {
    if ((*line)[0] == '>') { *have_header = 1; break; }
    while (l && ((*line)[l - 1] == '\n' || (*line)[l - 1] == '\r')) l--;
    if (r->n + (size_t)l + 1 > r->cap) {
      while (r->n + (size_t)l + 1 > r->cap) r->cap = r->cap ? 2 * r->cap : 1 << 20;
      r->s = realloc(r->s, r->cap);
      if (!r->s) { perror("realloc"); exit(1); }
    }
    memcpy(r->s + r->n, *line, (size_t)l);
    r->n += (size_t)l;
  }
  return 1;
}

static int cmd_tally(unsigned k, const char* pdiff)
{
  code_init();
  size_t nd;
  u64* d = read_u64(pdiff, &nd);
  unsigned* seen_in = calloc(nd + 1, sizeof(unsigned)); /* 1 + last record index that held it; 0 = never */
  const u64 mask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);
  const int top = 2 * (int)(k - 1);
  record r = {0};
  char* line = NULL;
  size_t lcap = 0;
  int hh = 0;
  unsigned rec = 0;
  u64 tot_first = 0;
  printf("record\tindex\tlength\twindows\twindows_in_diff\tdistinct_in_diff\tdistinct_first_seen\n");
  while (next_record(&r, &line, &lcap, &hh)) {
    rec++;
    u64 f = 0, rv = 0, run = 0, win = 0, hit = 0, dist = 0, first = 0;
    for (size_t i = 0; i < r.n; i++) {
      int c = CODE[(unsigned char)r.s[i]];
      if (c < 0) { run = 0; continue; }
      f = ((f << 2) | (u64)c) & mask;
      rv = (rv >> 2) | ((u64)(3 - c) << top);
      if (++run < k) continue;
      win++;
      long j = find(d, nd, f < rv ? f : rv);
      if (j < 0) continue;
      hit++;
      if (seen_in[j] != rec) { if (!seen_in[j]) first++; dist++; seen_in[j] = rec; }
    }
    tot_first += first;
    printf("%s\t%u\t%zu\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n", r.name, rec, r.n, win, hit, dist, first);
  }
  fprintf(stderr, "{\"cmd\": \"tally\", \"k\": %u, \"diff\": \"%s\", \"diff_size\": %zu, \"records\": %u, "
          "\"diff_kmers_found\": %" PRIu64 "}\n", k, pdiff, nd, rec, tot_first);
  return 0;
}

static int cmd_extscore(unsigned k, const char* pamb, const char* pbma, const char* pmiss)
{
  code_init();
  size_t namb, nbma;
  u64 *amb = read_u64(pamb, &namb), *bma = read_u64(pbma, &nbma);
  const u64 mask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);
  const int top = 2 * (int)(k - 1);
  vec ka = {0}, kb = {0};
  char* line = NULL;
  size_t lcap = 0;
  ssize_t l;
  u64 li = 0, lj = 0, bad = 0;
  while ((l = getline(&line, &lcap, stdin)) >= 0) {
    if (l < 3 || line[1] != ',' || (line[0] != 'i' && line[0] != 'j')) { bad++; continue; }
    vec* o = line[0] == 'i' ? &ka : &kb;
    if (line[0] == 'i') li++; else lj++;
    u64 f = 0, rv = 0, run = 0;
    for (ssize_t i = 2; i < l; i++) {
      int c = CODE[(unsigned char)line[i]];
      if (line[i] == '\n' || line[i] == '\r') break;
      if (c < 0) { run = 0; continue; }
      f = ((f << 2) | (u64)c) & mask;
      rv = (rv >> 2) | ((u64)(3 - c) << top);
      if (++run >= k) push(o, f < rv ? f : rv);
    }
  }
  size_t na = sort_unique(ka.v, ka.n), nb = sort_unique(kb.v, kb.n);
  u64* ga = malloc(na * 8 + 8);
  u64* gb = malloc(nb * 8 + 8);
  size_t nga = minus(ka.v, na, kb.v, nb, ga), ngb = minus(kb.v, nb, ka.v, na, gb);
  /* got = ga U gb (disjoint), truth = amb U bma (disjoint) */
  vec got = {0}, want = {0};
  for (size_t i = 0; i < nga; i++) push(&got, ga[i]);
  for (size_t i = 0; i < ngb; i++) push(&got, gb[i]);
  for (size_t i = 0; i < namb; i++) push(&want, amb[i]);
  for (size_t i = 0; i < nbma; i++) push(&want, bma[i]);
  size_t ng = sort_unique(got.v, got.n), nw = sort_unique(want.v, want.n);
  u64* miss = malloc(nw * 8 + 8);
  size_t nmiss = minus(want.v, nw, got.v, ng, miss), nspur = minus(got.v, ng, want.v, nw, NULL);
  FILE* mf = fopen(pmiss, "w");
  if (!mf) { perror(pmiss); return 1; }
  static const char ACGT[4] = {'A', 'C', 'G', 'T'};
  size_t miss_a = 0;
  for (size_t i = 0; i < nmiss; i++) {
    char s[33];
    for (unsigned t = 0; t < k; t++) s[t] = ACGT[(miss[i] >> (2 * (k - 1 - t))) & 3];
    s[k] = 0;
    int in_a = find(amb, namb, miss[i]) >= 0;
    miss_a += (size_t)in_a;
    fprintf(mf, "%s\t%s\n", s, in_a ? "A_minus_B" : "B_minus_A");
  }
  fclose(mf);
  printf("{\"cmd\": \"extscore\", \"k\": %u, \"listed_i\": %" PRIu64 ", \"listed_j\": %" PRIu64 ", \"unparsed\": %" PRIu64
         ", \"kmers_i\": %zu, \"kmers_j\": %zu, \"kmers_in_both\": %zu, \"recovered\": %zu, \"truth\": %zu, "
         "\"missing\": %zu, \"missing_A_minus_B\": %zu, \"spurious\": %zu}\n", k, li, lj, bad, na, nb, na - nga, ng, nw,
         nmiss, miss_a, nspur);
  return 0;
}

static int cmd_locate(unsigned k, const char* pk)
{
  code_init();
  FILE* f = fopen(pk, "r");
  if (!f) { perror(pk); return 1; }
  vec q = {0};
  char* line = NULL;
  size_t lcap = 0;
  ssize_t l;
  while ((l = getline(&line, &lcap, f)) >= 0) {
    if (l < (ssize_t)k) continue;
    u64 c = 0;
    int ok = 1;
    for (unsigned t = 0; t < k; t++) {
      int v = CODE[(unsigned char)line[t]];
      if (v < 0) { ok = 0; break; }
      c = (c << 2) | (u64)v;
    }
    if (ok) push(&q, c);
  }
  fclose(f);
  size_t nq = sort_unique(q.v, q.n);
  const u64 mask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);
  const int top = 2 * (int)(k - 1);
  static const char ACGT[4] = {'A', 'C', 'G', 'T'};
  record r = {0};
  int hh = 0;
  u64 hits = 0;
  printf("kmer\trecord\tpos\tstrand\tleft20\twindow\tright20\n");
  while (next_record(&r, &line, &lcap, &hh)) {
    u64 fw = 0, rv = 0, run = 0;
    for (size_t i = 0; i < r.n; i++) {
      int c = CODE[(unsigned char)r.s[i]];
      if (c < 0) { run = 0; continue; }
      fw = ((fw << 2) | (u64)c) & mask;
      rv = (rv >> 2) | ((u64)(3 - c) << top);
      if (++run < k) continue;
      u64 can = fw < rv ? fw : rv;
      if (find(q.v, nq, can) < 0) continue;
      hits++;
      size_t pos = i + 1 - k, lo = pos >= 20 ? pos - 20 : 0, hi = i + 1 + 20 < r.n ? i + 1 + 20 : r.n;
      char ks[33];
      for (unsigned t = 0; t < k; t++) ks[t] = ACGT[(can >> (2 * (k - 1 - t))) & 3];
      ks[k] = 0;
      printf("%s\t%s\t%zu\t%c\t%.*s\t%.*s\t%.*s\n", ks, r.name, pos, can == fw ? '+' : '-', (int)(pos - lo), r.s + lo,
             (int)k, r.s + pos, (int)(hi - i - 1), r.s + i + 1);
    }
  }
  fprintf(stderr, "{\"cmd\": \"locate\", \"k\": %u, \"queries\": %zu, \"occurrences\": %" PRIu64 "}\n", k, nq, hits);
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
  if (argc >= 4 && !strcmp(argv[1], "tally")) return cmd_tally((unsigned)atoi(argv[2]), argv[3]);
  if (argc >= 6 && !strcmp(argv[1], "extscore")) return cmd_extscore((unsigned)atoi(argv[2]), argv[3], argv[4], argv[5]);
  if (argc >= 4 && !strcmp(argv[1], "locate")) return cmd_locate((unsigned)atoi(argv[2]), argv[3]);
  fprintf(stderr, "usage: kmtools kset K | text K | diff A B AMB BMA | kmcdump K | same X Y | listscore K AMB BMA | "
                  "tally K DIFF | extscore K AMB BMA MISS | locate K KMERS\n");
  return 2;
}
