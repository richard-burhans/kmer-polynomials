// KP-6/KP-7 against ntHash: throughput, minimizer density and gaps, context-free sampling gaps, bit uniformity.
//
//   bench <fasta> [k=31] [rounds=5] [throughput] [only=<name substring>]
//     `throughput` skips the stats; `only=` times just the matching method(s), for profiling (perf record).
//   bench --nthash-classes                exhaustive k=8 class count for ntHash's canonical (fwd + rev)
//
// Every method sees identical bytes: the FASTA is loaded once, upper-cased (case-blind), split per record.
// Throughput: `rounds` interleaved passes over all methods, timing only the hashing loop; median reported.
// Stats (untimed): per record, runs of consecutive k-mer positions; random minimizers (leftmost minimum of w
// consecutive k-mers) at w = 10, 20, 50; FracMinHash-style threshold sampling at rate 2/(w+1) for w = 20 to
// look for "deserts" (long gaps), compared with the geometric expectation; per-bit frequency of the outputs.
#include "kp7.h"
#include <nthash/nthash.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <deque>
#include <fstream>
#include <map>
#include <set>
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

enum Method { NTHASH = 0, KP7_FWD_REAL, KP7_MIN_REAL, KP7_SYM_REAL, KP7_SYM_COMPLEX, KP7_SYM_REAL_K, KP7_SYM_COMPLEX_K,
              SYMR_OPT0, SYMR_LAM, SYMR_MERGED, SYMR_LAZY, SYMR_ALL, SYMC_OPT0, SYMC_LAM_MERGED, NMETHODS };
static const char* NAMES[NMETHODS] = { "ntHash canonical (fwd+rev)", "kp7 forward, x real", "kp7 canonical min, x real",
                                       "kp7 symmetric, x real", "kp7 symmetric, x complex",
                                       "kp7 sym, x real, Karatsuba", "kp7 sym, x cplx, Karatsuba",
                                       "roll_sym x real opt0", "roll_sym x real lamstate", "roll_sym x real merged",
                                       "roll_sym x real lazy", "roll_sym x real all three", "roll_sym x cplx opt0",
                                       "roll_sym x cplx lam+merged" };
// roll_sym methods (kp7_roll_sym with params.opt) must equal the plain symmetric hashes; main checks the sinks.
static int sym_opt(int m) { switch (m) { case SYMR_OPT0: case SYMC_OPT0: return 0; case SYMR_LAM: return KP7_OPT_LAMSTATE;
  case SYMR_MERGED: return KP7_OPT_MERGED; case SYMR_LAZY: return KP7_OPT_LAZY; case SYMR_ALL: return 7;
  case SYMC_LAM_MERGED: return KP7_OPT_LAMSTATE | KP7_OPT_MERGED; default: return -1; } }
static bool is_cplx_opt(int m) { return m == SYMC_OPT0 || m == SYMC_LAM_MERGED; }
// Karatsuba methods give bit-identical hashes (test_kp7), so their STATS equal the plain methods'; main checks the sinks.
static bool is_k(Method m) { return m == KP7_SYM_REAL_K || m == KP7_SYM_COMPLEX_K || sym_opt(m) >= 0; }

// Calls fn(pos, hash) for every valid k-mer of one record.
static const kp7_params *g_real_k, *g_cplx_k; // Karatsuba copies of real/cplx, set in main
static kp7_params g_opt[NMETHODS];             // per roll_sym method: real or cplx with that method's opt bits
template <class Fn> static void hash_record(Method m, const std::string& s, unsigned k, const kp7_params* real,
                                            const kp7_params* cplx, Fn&& fn)
{
  if (s.size() < k) return;
  if (m == NTHASH) {
    nthash::NtHash nth(s.data(), s.size(), 1, k);
    while (nth.roll()) fn(nth.get_pos(), nth.hashes()[0]);
    return;
  }
  if (sym_opt(m) >= 0) {
    kp7_roller r;
    kp7_start(&r, &g_opt[m], s.data(), s.size());
    uint64_t h;
    while (kp7_roll_sym(&r, &h)) fn(r.pos, h);
    return;
  }
  const kp7_params* p = m == KP7_SYM_COMPLEX ? cplx : m == KP7_SYM_REAL_K ? g_real_k : m == KP7_SYM_COMPLEX_K ? g_cplx_k : real;
  int mode = m == KP7_FWD_REAL ? KP7_OUT_FWD : m == KP7_MIN_REAL ? KP7_OUT_CANON_MIN : KP7_OUT_SYM;
  kp7_roller r;
  kp7_start(&r, p, s.data(), s.size());
  uint64_t h;
  while (kp7_roll64(&r, mode, &h)) fn(r.pos, h);
}

struct Stats {
  uint64_t kmers = 0;
  uint64_t bits[64] = { 0 };
  // minimizers
  static constexpr int NW = 3;
  int ws[NW] = { 10, 20, 50 };
  uint64_t min_kmers[NW] = { 0 }, min_sel[NW] = { 0 }, min_maxgap[NW] = { 0 };
  std::map<uint64_t, uint64_t> min_gaps20;
  // context-free sampling at rate 2/(w+1), w = 20
  uint64_t cf_sel = 0, cf_gaps = 0, cf_maxgap = 0, cf_gaps_over[4] = { 0 };
};
static const uint64_t CF_OVER[4] = { 50, 100, 200, 400 };

static void run_stats(Stats& st, const std::vector<uint64_t>& h)
{
  // h is one run of consecutive positions
  size_t L = h.size();
  for (int wi = 0; wi < Stats::NW; wi++) {
    size_t w = st.ws[wi];
    if (L < w) continue;
    st.min_kmers[wi] += L;
    std::deque<size_t> dq;
    size_t last = (size_t)-1;
    for (size_t i = 0; i < L; i++) {
      while (!dq.empty() && h[dq.back()] > h[i]) dq.pop_back(); // keep earlier equal values: leftmost minimum
      dq.push_back(i);
      if (dq.front() + w <= i) dq.pop_front();
      if (i + 1 >= w) {
        size_t sel = dq.front();
        if (sel != last) {
          if (last != (size_t)-1) {
            uint64_t g = sel - last;
            st.min_maxgap[wi] = std::max(st.min_maxgap[wi], g);
            if (w == 20) st.min_gaps20[g]++;
          }
          st.min_sel[wi]++;
          last = sel;
        }
      }
    }
  }
  const double rate = 2.0 / 21.0;
  const uint64_t thr = (uint64_t)(rate * 18446744073709551616.0);
  size_t last = (size_t)-1;
  for (size_t i = 0; i < L; i++) {
    if (h[i] < thr) {
      if (last != (size_t)-1) {
        uint64_t g = i - last;
        st.cf_gaps++;
        st.cf_maxgap = std::max(st.cf_maxgap, g);
        for (int j = 0; j < 4; j++) st.cf_gaps_over[j] += g > CF_OVER[j];
      }
      st.cf_sel++;
      last = i;
    }
  }
}

static Stats stats_for(Method m, const std::vector<Rec>& recs, unsigned k, const kp7_params* real,
                       const kp7_params* cplx)
{
  Stats st;
  std::vector<uint64_t> run;
  for (const Rec& r : recs) {
    size_t prev = (size_t)-2;
    run.clear();
    hash_record(m, r.seq, k, real, cplx, [&](size_t pos, uint64_t h) {
      if (pos != prev + 1 && !run.empty()) { run_stats(st, run); run.clear(); }
      run.push_back(h);
      prev = pos;
      st.kmers++;
      for (int b = 0; b < 64; b++) st.bits[b] += (h >> b) & 1;
    });
    if (!run.empty()) run_stats(st, run);
  }
  return st;
}

static int nthash_classes()
{
  // ntHash canonical = fwd + rev (mod 2^64). Count distinct values over all 4^8 8-mers; 32,896 classes expected.
  const unsigned k = 8;
  std::set<uint64_t> vals;
  char s[9];
  for (unsigned v = 0; v < (1u << 16); v++) {
    for (unsigned t = 0; t < k; t++) s[t] = "ACGT"[(v >> (2 * t)) & 3];
    s[k] = 0;
    nthash::NtHash nth(s, k, 1, k);
    nth.roll();
    vals.insert(nth.hashes()[0]);
  }
  printf("ntHash canonical, k=8: %zu distinct values over 65536 k-mers (32896 double-stranded classes)\n", vals.size());
  return 0;
}

int main(int argc, char** argv)
{
  if (argc > 1 && !strcmp(argv[1], "--nthash-classes")) return nthash_classes();
  if (argc < 2) { fprintf(stderr, "usage: bench <fasta> [k] [rounds]\n"); return 2; }
  unsigned k = argc > 2 ? (unsigned)atoi(argv[2]) : 31;
  int rounds = argc > 3 ? atoi(argv[3]) : 5;
  bool throughput_only = argc > 4 && !strcmp(argv[4], "throughput");
  const char* only = argc > 5 && !strncmp(argv[5], "only=", 5) ? argv[5] + 5 : nullptr;
  auto t0 = std::chrono::steady_clock::now();
  std::vector<Rec> recs = load_fasta(argv[1]);
  size_t bases = 0;
  for (auto& r : recs) bases += r.seq.size();
  printf("loaded %zu records, %zu bases in %.1f s; k=%u, rounds=%d\n", recs.size(), bases,
         std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), k, rounds);

  kp7_params real, cplx;
  kp7_init(&real, k, KP7_XREAL, 0x5eed);
  kp7_init(&cplx, k, KP7_XCOMPLEX, 0x5eed);
  kp7_params real_k = real, cplx_k = cplx;
  real_k.karatsuba = cplx_k.karatsuba = 1;
  g_real_k = &real_k; g_cplx_k = &cplx_k;
  for (int m = 0; m < NMETHODS; m++)
    if (sym_opt(m) >= 0) { g_opt[m] = is_cplx_opt(m) ? cplx : real; g_opt[m].opt = sym_opt(m); }

  // throughput: interleaved rounds
  std::vector<std::vector<double>> secs(NMETHODS);
  uint64_t counts[NMETHODS] = { 0 }, sink = 0, msink[NMETHODS] = { 0 };
  for (int round = 0; round < rounds; round++) {
    for (int m = 0; m < NMETHODS; m++) {
      if (only && !strstr(NAMES[m], only)) continue;
      uint64_t n = 0, x = 0;
      auto a = std::chrono::steady_clock::now();
      for (const Rec& r : recs) hash_record((Method)m, r.seq, k, &real, &cplx, [&](size_t, uint64_t h) { x ^= h; n++; });
      double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
      secs[m].push_back(s);
      counts[m] = n;
      sink ^= x;
      msink[m] = x;
      printf("  round %d %-30s %.3f s  %.1f M k-mers/s\n", round, NAMES[m], s, n / s / 1e6);
      fflush(stdout);
    }
  }
  printf("\nTHROUGHPUT (median of %d interleaved rounds; sink %016llx)\n", rounds, (unsigned long long)sink);
  for (int m = 0; m < NMETHODS; m++) {
    if (secs[m].empty()) continue;
    std::vector<double> v = secs[m];
    std::sort(v.begin(), v.end());
    double med = v[v.size() / 2];
    printf("  %-30s k-mers %llu  median %.3f s  %.1f M k-mers/s\n", NAMES[m], (unsigned long long)counts[m], med,
           counts[m] / med / 1e6);
  }

  // the Karatsuba methods must reproduce the plain ones exactly (xor of every hash)
  if (only) return 0; // the identity check needs every method
  bool same = msink[KP7_SYM_REAL_K] == msink[KP7_SYM_REAL] && msink[KP7_SYM_COMPLEX_K] == msink[KP7_SYM_COMPLEX];
  for (int m = 0; m < NMETHODS; m++)
    if (sym_opt(m) >= 0) same = same && msink[m] == msink[is_cplx_opt(m) ? KP7_SYM_COMPLEX : KP7_SYM_REAL];
  printf("  Karatsuba and roll_sym hashes identical to plain: %s\n", same ? "yes" : "NO");
  if (!same) return 1;
  if (throughput_only) return 0;
  printf("\nSTATS\n");
  for (int m = 0; m < NMETHODS; m++) {
    if (is_k((Method)m)) continue;
    Stats st = stats_for((Method)m, recs, k, &real, &cplx);
    double maxdev = 0;
    for (int b = 0; b < 64; b++) maxdev = std::max(maxdev, std::fabs((double)st.bits[b] / st.kmers - 0.5));
    printf("%s\n  k-mers %llu, max per-bit deviation from 0.5: %.5f\n", NAMES[m], (unsigned long long)st.kmers, maxdev);
    for (int wi = 0; wi < Stats::NW; wi++) {
      double d = (double)st.min_sel[wi] / st.min_kmers[wi], e = 2.0 / (st.ws[wi] + 1);
      printf("  minimizer w=%-2d density %.5f (expected %.5f, ratio %.4f), max gap %llu\n", st.ws[wi], d, e, d / e,
             (unsigned long long)st.min_maxgap[wi]);
    }
    uint64_t g20 = 0, s20 = 0;
    for (auto& [g, c] : st.min_gaps20) { g20 += c; s20 += g * c; }
    printf("  w=20 gap mean %.3f; gap histogram:", (double)s20 / g20);
    for (auto& [g, c] : st.min_gaps20) printf(" %llu:%.4f", (unsigned long long)g, (double)c / g20);
    printf("\n");
    const double rate = 2.0 / 21.0;
    printf("  context-free sampling at rate %.5f: selected %.5f of k-mers, gaps %llu, max gap %llu\n", rate,
           (double)st.cf_sel / st.kmers, (unsigned long long)st.cf_gaps, (unsigned long long)st.cf_maxgap);
    for (int j = 0; j < 4; j++)
      printf("    gaps > %3llu: %llu (geometric expectation %.1f)\n", (unsigned long long)CF_OVER[j],
             (unsigned long long)st.cf_gaps_over[j], st.cf_gaps * std::pow(1 - rate, (double)CF_OVER[j]));
    fflush(stdout);
  }
  return 0;
}
