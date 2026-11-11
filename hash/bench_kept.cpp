// Kept-k-mer evaluation (kp7kept) against batched ntHash: per INPUT k-mer, for laila_perf counter runs.
//
//   bench_kept <fasta> <k> only=<method>     time one method once (or only=__load_only__: load the FASTA, nothing else)
//   bench_kept <fasta> <k> check             identity of every evaluator against its reference sampler, whole input
//
// Methods:
//   "ntHash batched"                         ntbatch.h, bit-identical to NtHash::roll() (hash/nteq)
//   "uh all"                                 min(fwd, rc) through multiply-shift, EVERY k-mer (round-1 critic gap 2)
//   "<ev> <B> <bat|br> <rate>[ b<N>]"        the hybrid-v2 sampler: packed-code gate, then <ev> on the kept k-mers
//       ev:   lib-e (kp7e_sample_batch / kp7e_sample_b, B = 8: the current Eisenstein sampler)
//             lib-m (kp7pack_sample_m / _mb, B = 8: the CM31 16 KB tables)
//             ec, ep, gp, mp, mc, nt, uh (kp7kept: E_CHUNK, E_PLANES, G_PLANES, M_PLANES, M_CHUNK, NT_CHUNK, UH)
//       rate: 0 (nothing kept: the base cost), 1/1000, 1/100, 1/16, 1/4 (threshold 2^64 * rate, as bench_e)
//       bN:   batch size for bat (default 512)
// The FASTA is loaded once and upper-cased (as bench_e). The sink XORs every hash and counts them; it is the same
// for every method. Prints "RESULT <method>\tkept=<n>\txor=<hex>\tsecs=<s>".
#include "kp7e.h"
#include "kp7kept.h"
#include "kp7pack.h"
#include "ntbatch.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

struct Rec { std::string name, seq; };
static std::vector<Rec> load_fasta(const char* path)
{
  std::ifstream in(path);
  std::vector<Rec> recs;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '>') { recs.push_back({ line.substr(1, line.find_first_of(" \t") - 1), "" }); continue; }
    for (char& c : line) c = (char)toupper((unsigned char)c);
    recs.back().seq += line;
  }
  return recs;
}

struct Sink { uint64_t x = 0, n = 0; };
// The sink works in locals: s->x and s->n may alias h (both uint64_t), so updating them through the pointer made the
// compiler load and store both per hash, about 10 instructions per k-mer that bench_e's lambda sink does not pay (the
// first sweep, job 3dddf48cdf81c0ca, measured batched ntHash at 57.1 instructions/k-mer here against 47.0 there).
static void pairs_sink(void* ctx, const size_t*, const uint64_t* h, size_t n)
{
  Sink* s = static_cast<Sink*>(ctx);
  uint64_t x = s->x;
  for (size_t j = 0; j < n; j++) x ^= h[j];
  s->x = x;
  s->n += n;
}
static void seq_sink(void* ctx, size_t, const uint64_t* h, size_t n) { pairs_sink(ctx, nullptr, h, n); }

struct Spec {
  std::string ev;
  unsigned B = 8;
  bool branchy = false;
  uint64_t thr = 0;
  size_t batch = 512;
};
static uint64_t parse_rate(const std::string& r)
{
  if (r == "0") return 0;
  if (r == "1") return ~0ULL;
  if (r.rfind("1/", 0) == 0) {
    double d = atof(r.c_str() + 2);
    if (d == 16) return 1ULL << 60;
    if (d == 4) return 1ULL << 62;
    return (uint64_t)(18446744073709551615.0 / d);
  }
  fprintf(stderr, "bad rate %s\n", r.c_str());
  exit(2);
}
static bool parse_spec(const std::string& m, Spec* s)
{
  char ev[16], col[8], rate[16], extra[16] = "";
  unsigned B;
  int n = sscanf(m.c_str(), "%15s %u %7s %15s %15s", ev, &B, col, rate, extra);
  if (n < 4) return false;
  s->ev = ev;
  s->B = B;
  s->branchy = !strcmp(col, "br");
  if (!s->branchy && strcmp(col, "bat")) return false;
  s->thr = parse_rate(rate);
  if (n == 5) { if (extra[0] != 'b') return false; s->batch = (size_t)atol(extra + 1); }
  return true;
}

static kp7e_params g_pe;
static kp7e_tab g_te;
static kp7_params g_pg;
static kp7m31_params g_pm;
static kp7pack_mtab g_mtab;
static kp7pack_ftab g_ftab;

static bool make_tab(const Spec& s, unsigned k, kp7kept* t)
{
  int ok = 0;
  if (s.ev == "ec") ok = kp7kept_init_e_chunk(t, &g_pe, s.B);
  else if (s.ev == "ep") ok = kp7kept_init_e_planes(t, &g_pe, s.B);
  else if (s.ev == "gp") ok = kp7kept_init_g_planes(t, &g_pg, s.B);
  else if (s.ev == "mp") ok = kp7kept_init_m_planes(t, &g_pm, s.B);
  else if (s.ev == "mc") ok = kp7kept_init_m_chunk(t, &g_pm, s.B);
  else if (s.ev == "nt") ok = kp7kept_init_nt(t, k, s.B);
  else if (s.ev == "uh") ok = kp7kept_init_uh(t, k, 0xC0FFEE);
  if (ok) t->batch = s.batch;
  return ok;
}

// run one method over every record; returns false for an unknown method
static bool run(const std::string& m, const std::vector<Rec>& recs, unsigned k, Sink* sink, double* secs)
{
  static uint64_t hbuf[4096];
  static size_t pbuf[4096];
  static ntbatch::Tables nt(k);
  if (m == "__load_only__") { *secs = 0; return true; }
  auto a = std::chrono::steady_clock::now();
  if (m == "ntHash batched") {
    for (const Rec& r : recs) ntbatch::hash_seq(nt, r.seq.data(), r.seq.size(), hbuf, 4096, seq_sink, sink);
  } else if (m == "uh all") {
    kp7kept t;
    kp7kept_init_uh(&t, k, 0xC0FFEE);
    a = std::chrono::steady_clock::now();
    for (const Rec& r : recs) kp7kept_uh_all(&t, r.seq.data(), r.seq.size(), hbuf, 4096, seq_sink, sink);
  } else {
    Spec s;
    if (!parse_spec(m, &s)) return false;
    if (s.ev == "lib-e") {
      if (s.B != 8) return false;
      for (const Rec& r : recs)
        (s.branchy ? kp7e_sample_b : kp7e_sample_batch)(&g_te, r.seq.data(), r.seq.size(), s.thr, pbuf, hbuf, 4096,
                                                         pairs_sink, sink);
    } else if (s.ev == "lib-m") {
      if (s.B != 8) return false;
      for (const Rec& r : recs) {
        if (s.branchy) kp7pack_sample_mb(&g_mtab, r.seq.data(), r.seq.size(), s.thr, pbuf, hbuf, 4096, pairs_sink, sink);
        else kp7pack_sample_m(&g_mtab, 8, r.seq.data(), r.seq.size(), s.thr, pbuf, hbuf, 4096, pairs_sink, sink);
      }
    } else if (s.ev == "lib-g") {
      if (s.B != 8 || s.branchy) return false;
      for (const Rec& r : recs)
        kp7pack_sample_f(&g_ftab, 8, r.seq.data(), r.seq.size(), s.thr, pbuf, hbuf, 4096, pairs_sink, sink);
    } else {
      kp7kept t;
      if (!make_tab(s, k, &t)) return false;
      a = std::chrono::steady_clock::now(); // table construction is not timed by the clock (it is by the counters)
      for (const Rec& r : recs)
        kp7kept_sample(&t, s.branchy, r.seq.data(), r.seq.size(), s.thr, pbuf, hbuf, 4096, pairs_sink, sink);
      kp7kept_free(&t);
    }
  }
  *secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
  return true;
}

int main(int argc, char** argv)
{
  if (argc < 4) { fprintf(stderr, "usage: bench_kept <fasta> <k> only=<method> | check\n"); return 2; }
  const unsigned k = (unsigned)atoi(argv[2]);
  if (k < 1 || k > 32) { fprintf(stderr, "k must be 1..32\n"); return 2; }
  kp7e_init(&g_pe, k, 0x5eedE, 0x5eedE + 1);
  kp7e_tab_init(&g_te, &g_pe);
  kp7_init(&g_pg, k, KP7_XREAL, 0x5eed, 0x5eed + 1);
  kp7pack_ftab_init(&g_ftab, &g_pg);
  kp7m31_init(&g_pm, k, 0x5eed31, 0x5eed31 + 1);
  kp7pack_mtab_init(&g_mtab, &g_pm);
  std::vector<Rec> recs = load_fasta(argv[1]);
  size_t bases = 0;
  for (auto& r : recs) bases += r.seq.size();
  printf("loaded %zu records, %zu bases; k=%u\n", recs.size(), bases, k);

  if (!strncmp(argv[3], "only=", 5)) {
    const std::string m = argv[3] + 5;
    Sink s;
    double secs;
    if (!run(m, recs, k, &s, &secs)) { fprintf(stderr, "unknown method '%s'\n", m.c_str()); return 2; }
    printf("RESULT %s\tkept=%llu\txor=%016llx\tsecs=%.3f\n", m.c_str(), (unsigned long long)s.n,
           (unsigned long long)s.x, secs);
    return 0;
  }
  if (strcmp(argv[3], "check")) { fprintf(stderr, "third argument: only=<method> or check\n"); return 2; }

  // identity: each new evaluator against its family's existing sampler, batched and branchy, at 1/16 and 1/1000
  bool all = true;
  auto one = [&](const std::string& m) { Sink s; double t; if (!run(m, recs, k, &s, &t)) { printf("  unknown %s\n", m.c_str()); all = false; } return s; };
  for (const char* rate : { "1/16", "1/1000" }) {
    const std::string R = rate;
    const Sink e = one("lib-e 8 bat " + R), m = one("lib-m 8 bat " + R), g = one("lib-g 8 bat " + R);
    // ntHash and UH references at the gate positions
    Sink n, u;
    {
      const uint64_t thr = parse_rate(R);
      static uint64_t hb[4096];
      kp7kept ut;
      kp7kept_init_uh(&ut, k, 0xC0FFEE);
      for (const Rec& r : recs) {
        struct C { const std::string* s; unsigned k; uint64_t thr; Sink* n; } c{ &r.seq, k, thr, &n };
        ntbatch::Tables nt(k);
        ntbatch::hash_seq(nt, r.seq.data(), r.seq.size(), hb, 4096,
                          [](void* ctx, size_t start, const uint64_t* h, size_t cnt) {
                            C* c = static_cast<C*>(ctx);
                            for (size_t j = 0; j < cnt; j++) {
                              uint64_t f, rc;
                              if (kp7pack_codes(c->s->data() + start + j, c->k, &f, &rc) && kp7pack_gate(f, rc) < c->thr) {
                                c->n->x ^= h[j]; c->n->n++;
                              }
                            }
                          }, &c);
        for (size_t pos = 0; pos + k <= r.seq.size(); pos++) {
          uint64_t f, rc;
          if (!kp7pack_codes(r.seq.data() + pos, k, &f, &rc) || kp7pack_gate(f, rc) >= thr) continue;
          u.x ^= kp7kept_eval(&ut, f, rc); // the formula itself is checked against __int128 in test_kp7kept
          u.n++;
        }
      }
    }
    printf("  %s: reference kept %llu (Eisenstein xor %016llx)\n", rate, (unsigned long long)e.n, (unsigned long long)e.x);
    struct { const char* ev; const Sink* want; std::vector<unsigned> Bs; } fams[] = {
      { "lib-e", &e, { 8 } }, { "ec", &e, { 8, 10, 11, 12, 13, 16 } }, { "ep", &e, { 8, 11, 16 } },
      { "gp", &g, { 8, 11, 16 } }, { "mp", &m, { 8, 11, 16 } }, { "mc", &m, { 8, 11, 16 } }, { "lib-m", &m, { 8 } },
      { "nt", &n, { 8, 16 } }, { "uh", &u, { 0 } } };
    for (auto& f : fams)
      for (unsigned B : f.Bs)
        for (const char* col : { "bat", "br" }) {
          if (!strcmp(f.ev, "lib-e") && !strcmp(col, "bat")) continue; // that is the reference itself
          const std::string name = std::string(f.ev) + " " + std::to_string(B) + " " + col + " " + R;
          const Sink got = one(name);
          const bool same = got.x == f.want->x && got.n == f.want->n;
          printf("  %-24s kept %llu xor %016llx: %s\n", name.c_str(), (unsigned long long)got.n,
                 (unsigned long long)got.x, same ? "identical" : "DIFFERENT");
          all = all && same;
        }
  }
  printf(all ? "check: every evaluator identical to its reference\n" : "check: FAILED\n");
  return all ? 0 : 1;
}
