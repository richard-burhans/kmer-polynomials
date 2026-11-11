// Multi-seed sampling statistics for strand-symmetric k-mer hashes on a genome (KP-7 next steps, item 3).
//
//   stats_seeds <fasta> [k=31] [seeds=5] [only=<method>]
//
// Methods, each with `seeds` different hash functions:
//   ntHash   canonical (fwd + rev) from the library; seed i > 0 uses ntHash's own multi-hash extension
//            h_i = h * (i ^ k * MULTISEED); h_i ^= h_i >> 27 (seed 0 = the canonical hash itself)
//   v2gate   kp7pack_gate_seeded(fwd, rc, seed)            (packed codes, mix(fwd ^ s) + mix(rc ^ s))
//   v1gate   kp7_gate_seq with kp7_init(seed)              (gate on the rolling field values)
//   provable kp7_sym_seq with kp7_init(seed)               (the provable (F + lam)(R + lam), fourth-root over F_{q^2})
//   eprovable kp7e_sym_seq with kp7e_init(seed)            (the Eisenstein provable S over F_q, q = 2^61 - 1). S is a
//            61-bit field element, so it is fed to the statistics as (S << 3) | (S >> 58): every bit then carries a
//            uniform field bit, and the order (minimizers, thresholds) is that of S up to its low 3 bits.
// For each: max per-bit deviation from 1/2; random-minimizer density at w = 10, 20, 50 (expected 2/(w+1)); and
// context-free (FracMinHash-style) sampling at rate 2/21: number of gaps, gaps > 200 k-mers (geometric expectation
// n (1 - rate)^200), and how many of those long gaps lie OUTSIDE repeats. A gap counts as inside a repeat when at least
// half of the bases it spans (from the selected k-mer before it to the end of the one after) are soft-masked (lowercase)
// in the input. The FASTA's case is kept for that; hashing is case-blind (every method treats a/A alike).
#include "kp7.h"
#include "kp7pack.h"
#include "kp7e.h"
#include <nthash/nthash.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <deque>
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
    recs.back().seq += line; // case KEPT: lowercase = soft-masked repeat
  }
  return recs;
}

static const double RATE = 2.0 / 21.0;
static const uint64_t THR = (uint64_t)(RATE * 18446744073709551615.0);
static const uint64_t LONG = 200;

struct Stats {
  uint64_t kmers = 0, bits[64] = { 0 };
  int ws[3] = { 10, 20, 50 };
  uint64_t windows[3] = { 0 }, sel[3] = { 0 };
  uint64_t gaps = 0, long_gaps = 0, long_outside = 0, maxgap = 0;
};

// h: hashes of one run of consecutive k-mer positions starting at `start` in record `seq`
static void run_stats(Stats& st, const std::vector<uint64_t>& h, size_t start, const std::string& seq, unsigned k)
{
  for (uint64_t v : h) { st.kmers++; for (int b = 0; b < 64; b++) st.bits[b] += (v >> b) & 1; }
  for (int wi = 0; wi < 3; wi++) {
    size_t w = st.ws[wi];
    if (h.size() < w) continue;
    std::deque<size_t> dq;
    long last = -1;
    for (size_t i = 0; i < h.size(); i++) {
      while (!dq.empty() && h[dq.back()] > h[i]) dq.pop_back();
      dq.push_back(i);
      if (dq.front() + w <= i) dq.pop_front();
      if (i + 1 >= w) { st.windows[wi]++; if ((long)dq.front() != last) { st.sel[wi]++; last = (long)dq.front(); } }
    }
  }
  long prev = -1;
  for (size_t i = 0; i < h.size(); i++) {
    if (h[i] >= THR) continue;
    if (prev >= 0) {
      uint64_t g = i - prev;
      st.gaps++;
      st.maxgap = std::max(st.maxgap, g);
      if (g > LONG) {
        st.long_gaps++;
        size_t a = start + prev, b = start + i + k; // bases spanned
        size_t lower = 0;
        for (size_t q = a; q < b; q++) lower += seq[q] >= 'a';
        if (2 * lower < b - a) st.long_outside++;
      }
    }
    prev = (long)i;
  }
}

enum { NTHASH = 0, V2GATE, V1GATE, PROVABLE, EPROVABLE, NM };
static const char* NAMES[NM] = { "ntHash", "v2gate", "v1gate", "provable", "eprovable" };
static kp7e_params g_ep; // set per seed in main for EPROVABLE

// calls fn(pos, hash) for every valid k-mer of one record
template <class Fn> static void hash_record(int m, int seed, const std::string& s, unsigned k, const kp7_params* p, Fn&& fn)
{
  if (s.size() < k) return;
  static thread_local uint64_t buf[4096];
  auto tramp = [](void* ctx, size_t start, const uint64_t* h, size_t n) {
    Fn& f = *static_cast<Fn*>(ctx);
    for (size_t j = 0; j < n; j++) f(start + j, h[j]);
  };
  if (m == NTHASH) {
    const uint64_t MULTISEED = 0x90b45d39fb6da1faULL; // ntHash src/internal.hpp (MIT)
    nthash::NtHash nth(s.data(), s.size(), 1, k);
    while (nth.roll()) {
      uint64_t h = nth.hashes()[0];
      if (seed > 0) { h = h * ((uint64_t)seed ^ (uint64_t)k * MULTISEED); h ^= h >> 27; }
      fn(nth.get_pos(), h);
    }
  } else if (m == V2GATE) {
    uint64_t sd = seed == 0 ? 0 : 0x9E3779B97F4A7C15ULL * (uint64_t)seed;
    // packed codes rolled here (case-blind, A0 C1 G2 T3), gate via the library's seeded gate
    const uint64_t mask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);
    uint64_t f = 0, r = 0;
    unsigned have = 0;
    for (size_t i = 0; i < s.size(); i++) {
      int c;
      switch (s[i]) { case 'A': case 'a': c = 0; break; case 'C': case 'c': c = 1; break;
                      case 'G': case 'g': c = 2; break; case 'T': case 't': c = 3; break; default: c = -1; }
      if (c < 0) { have = 0; f = r = 0; continue; }
      f = ((f << 2) | (uint64_t)c) & mask;
      r = (r >> 2) | ((uint64_t)(3 - c) << (2 * (k - 1)));
      if (have < k) have++;
      if (have == k) fn(i + 1 - k, kp7pack_gate_seeded(f, r, sd));
    }
  } else if (m == EPROVABLE) {
    auto etramp = [](void* ctx, size_t start, const uint64_t* h, size_t n) {
      Fn& f = *static_cast<Fn*>(ctx);
      for (size_t j = 0; j < n; j++) f(start + j, (h[j] << 3) | (h[j] >> 58));
    };
    kp7e_sym_seq(&g_ep, s.data(), s.size(), 2, buf, 4096, etramp, &fn);
  } else if (m == V1GATE) {
    kp7_gate_seq(p, s.data(), s.size(), 2, buf, 4096, tramp, &fn);
  } else {
    kp7_sym_seq(p, s.data(), s.size(), 2, buf, 4096, tramp, &fn);
  }
}

int main(int argc, char** argv)
{
  if (argc < 2) { fprintf(stderr, "usage: stats_seeds <fasta> [k] [seeds] [only=<method>]\n"); return 2; }
  unsigned k = argc > 2 ? (unsigned)atoi(argv[2]) : 31;
  int seeds = argc > 3 ? atoi(argv[3]) : 5;
  const char* only = argc > 4 && !strncmp(argv[4], "only=", 5) ? argv[4] + 5 : nullptr;
  std::vector<Rec> recs = load_fasta(argv[1]);
  size_t bases = 0, lower = 0;
  for (auto& r : recs) { bases += r.seq.size(); for (char c : r.seq) lower += c >= 'a' && c <= 'z'; }
  printf("loaded %zu records, %zu bases, %.1f%% soft-masked; k=%u, seeds=%d, rate %.5f, long gap > %llu\n", recs.size(),
         bases, 100.0 * lower / bases, k, seeds, RATE, (unsigned long long)LONG);
  printf("method\tseed\tkmers\tmax_bit_dev\tdens10\tdens20\tdens50\tgaps\tlong_gaps\tlong_expected\tlong_outside_repeats\tmax_gap\n");
  for (int m = 0; m < NM; m++) {
    if (only && strcmp(only, NAMES[m])) continue;
    for (int seed = 0; seed < seeds; seed++) {
      kp7_params p;
      kp7_init(&p, k, KP7_XREAL, 0x5eed + 1000 * (uint64_t)seed);
      kp7e_init(&g_ep, k, 0x5eed + 1000 * (uint64_t)seed);
      Stats st;
      for (const Rec& r : recs) {
        std::vector<uint64_t> run;
        size_t prev = SIZE_MAX, start = 0;
        hash_record(m, seed, r.seq, k, &p, [&](size_t pos, uint64_t h) {
          if (prev != SIZE_MAX && pos != prev + 1) { run_stats(st, run, start, r.seq, k); run.clear(); }
          if (run.empty()) start = pos;
          run.push_back(h);
          prev = pos;
        });
        if (!run.empty()) run_stats(st, run, start, r.seq, k);
      }
      double maxdev = 0;
      for (int b = 0; b < 64; b++) maxdev = std::max(maxdev, std::fabs((double)st.bits[b] / st.kmers - 0.5));
      double d[3];
      for (int wi = 0; wi < 3; wi++) d[wi] = (double)st.sel[wi] / st.windows[wi] / (2.0 / (st.ws[wi] + 1));
      printf("%s\t%d\t%llu\t%.5f\t%.4f\t%.4f\t%.4f\t%llu\t%llu\t%.2f\t%llu\t%llu\n", NAMES[m], seed,
             (unsigned long long)st.kmers, maxdev, d[0], d[1], d[2], (unsigned long long)st.gaps,
             (unsigned long long)st.long_gaps, st.gaps * std::pow(1 - RATE, (double)LONG),
             (unsigned long long)st.long_outside, (unsigned long long)st.maxgap);
      fflush(stdout);
    }
  }
  return 0;
}
