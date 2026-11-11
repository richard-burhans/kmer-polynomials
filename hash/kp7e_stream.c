/* Stream every k-mer of a FASTA, AS WRITTEN, as Eisenstein KP-7 values. No 2-bit packing and no canonicalization:
 * each value is computed from the sequence bytes by the rolling hash, so any k works (k > 32 too).
 *
 *   kp7e_stream <k> <sym|fwd> <seed_x|os> <seed_lam|os> <part_bits> <lo> <hi>  < in.fa  > values.u64
 *
 *   sym  S = (F + lam)(R + lam), the strand-symmetric hash (kp7e_sym_seq, bit-identical to kp7e_sym over kp7e_roll)
 *   fwd  F = sum_t s_t x^t, the forward polynomial alone: strand-SPECIFIC, the control (kp7e_roll's F)
 *
 * Seeds are parsed with strtoull base 0; "os" draws that parameter from OS randomness (kp7e_init_seeds). Windows with
 * a non-ACGT byte are skipped; case is ignored. A value v (canonical, < q = 2^61 - 1) belongs to partition
 * v >> (61 - part_bits); values of partitions lo <= P < hi are written to stdout as native uint64, duplicates kept, so
 * a caller can hold a slice of the genome's values in memory (recon/aswritten.py deduplicates). Writing all
 * partitions to files instead ran at ~4 MB/s on this host's virtiofs. One line of JSON on stderr: k, mode, x, lam,
 * records, bases, kmers (all k-mers), written (those in [lo, hi)). */
#include "kp7e.h"
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFN 8192

typedef struct {
  unsigned bits, lo, hi;
  uint64_t buf[BUFN];
  size_t n;
  uint64_t kmers, written;
} parts_t;

static void put(parts_t* P, uint64_t v)
{
  unsigned b = (unsigned)(v >> (61 - P->bits));
  P->kmers++;
  if (b < P->lo || b >= P->hi) return;
  P->buf[P->n++] = v;
  P->written++;
  if (P->n == BUFN) {
    if (fwrite(P->buf, 8, BUFN, stdout) != BUFN) { perror("fwrite"); exit(1); }
    P->n = 0;
  }
}

static void emit(void* ctx, size_t start, const uint64_t* h, size_t n)
{
  (void)start;
  for (size_t i = 0; i < n; i++) put((parts_t*)ctx, h[i]);
}

static void hash_record(const kp7e_params* p, int sym, const char* seq, size_t len, parts_t* P, uint64_t* hb)
{
  if (sym) {
    kp7e_sym_seq(p, seq, len, 2, hb, 4096, emit, P);
  } else {
    kp7e_roller r;
    kp7e_start(&r, p, seq, len);
    while (kp7e_roll(&r)) put(P, r.F);
  }
}

static int parse_seed(const char* s, uint64_t* out)
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
  if (argc != 8) {
    fprintf(stderr, "usage: kp7e_stream <k> <sym|fwd> <seed_x|os> <seed_lam|os> <part_bits> <lo> <hi> < in.fa\n");
    return 2;
  }
  unsigned k = (unsigned)atoi(argv[1]);
  int sym = !strcmp(argv[2], "sym");
  if (!sym && strcmp(argv[2], "fwd")) { fprintf(stderr, "mode must be sym or fwd\n"); return 2; }
  uint64_t sx, sl;
  int hx = parse_seed(argv[3], &sx), hl = parse_seed(argv[4], &sl);
  unsigned bits = (unsigned)atoi(argv[5]), lo = (unsigned)atoi(argv[6]), hi = (unsigned)atoi(argv[7]);
  if (k < 1 || bits < 1 || bits > 16 || lo >= hi || hi > (1u << bits)) { fprintf(stderr, "bad k, part_bits or range\n"); return 2; }
  kp7e_params p;
  kp7e_init_seeds(&p, k, hx ? &sx : NULL, hl ? &sl : NULL);
  static parts_t P;
  P.bits = bits; P.lo = lo; P.hi = hi;
  uint64_t* hb = malloc(4096 * 8);

  size_t cap = 1 << 20, len = 0;
  char* seq = malloc(cap);
  char* line = NULL;
  size_t lcap = 0;
  ssize_t l;
  uint64_t records = 0, bases = 0;
  while ((l = getline(&line, &lcap, stdin)) >= 0) {
    if (l && line[0] == '>') {
      if (records) hash_record(&p, sym, seq, len, &P, hb);
      records++;
      bases += len;
      len = 0;
      continue;
    }
    while (l && (line[l - 1] == '\n' || line[l - 1] == '\r')) l--;
    if (len + (size_t)l > cap) { while (len + (size_t)l > cap) cap *= 2; seq = realloc(seq, cap); }
    memcpy(seq + len, line, (size_t)l);
    len += (size_t)l;
  }
  if (records) { hash_record(&p, sym, seq, len, &P, hb); bases += len; }
  if (P.n && fwrite(P.buf, 8, P.n, stdout) != P.n) { perror("fwrite"); return 1; }
  if (fflush(stdout)) { perror("fflush"); return 1; }
  fprintf(stderr, "{\"k\": %u, \"mode\": \"%s\", \"x\": %" PRIu64 ", \"lam\": %" PRIu64 ", \"records\": %" PRIu64
          ", \"bases\": %" PRIu64 ", \"kmers\": %" PRIu64 ", \"written\": %" PRIu64 "}\n",
          k, sym ? "sym" : "fwd", p.x, p.lam, records, bases, P.kmers, P.written);
  return 0;
}
