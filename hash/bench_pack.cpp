// Hybrid v2 (kp7pack) against ntHash: throughput per INPUT k-mer, plus gate stats.
//
//   bench_pack <fasta> [k=31] [rounds=5] [throughput] [only=<name substring>]
//     `throughput` skips the stats; `only=` times just the matching method(s), for laila_perf counter runs.
//
// The FASTA is loaded once and upper-cased; every method sees identical bytes. Throughput: `rounds` interleaved passes,
// timing only the hashing loop, median reported. Sample methods emit only kept k-mers, so their rate is per input k-mer.
// Identity: "pack sample 1/16" emits provable KP-7 hashes, and their XOR must equal the XOR of the kp7_sym_seq (L2)
// hashes at the same positions, recomputed here with the gate as the filter.
// Stats (untimed, gate and ntHash only): random-minimizer density at w = 10, 20, 50 (expected 2/(w+1)) and per-bit
// frequency.
#include "kp7pack.h"
#include "ntbatch.h"
#include <nthash/nthash.hpp>

#include <algorithm>
#include <chrono>
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
    for (char& c : line) c = (char)toupper((unsigned char)c);
    recs.back().seq += line;
  }
  return recs;
}

enum Method { NTHASH = 0, PACK_GATE, PACK_S16, PACK_S1000, KP7_SEQ_L2, NT_BATCH, PACK_S100, CHUNK_S1000, CHUNK_S100,
              CHUNK_S16, POS_S16, CHUNK_ALL,
              // kept-k-mer evaluation (results/kp7_kept_eval_2026-11-10.md): f = 61-bit, m = CM31; 8/4 = chunk bits
              F8_ALL, F4_ALL, M8_ALL, M4_ALL, M31_REF, F8_S1000, F8_S100, F8_S16, F8_S4, F4_S16, M8_S1000, M8_S100, M8_S16,
              M8_S4, M4_S16, CHUNK_S4, NMETHODS };
static const char* NAMES[NMETHODS] = { "ntHash canonical (fwd+rev)", "pack gate (cheap sym)", "pack sample 1/16",
                                       "pack sample 1/1000", "kp7 sym_seq x real L2 (provable, all)",
                                       "ntHash batched (ntbatch)", "pack sample 1/100",
                                       "pack chunk-eval sample 1/1000", "pack chunk-eval sample 1/100",
                                       "pack chunk-eval sample 1/16", "pack pos-eval sample 1/16",
                                       "pack chunk-eval all (provable, all)",
                                       "fast f8 all", "fast f4 all", "fast m8 all (CM31)", "fast m4 all (CM31)",
                                       "m31 sym_seq scalar (CM31 reference)", "fast f8 sample 1/1000",
                                       "fast f8 sample 1/100", "fast f8 sample 1/16", "fast f8 sample 1/4",
                                       "fast f4 sample 1/16", "fast m8 sample 1/1000 (CM31)", "fast m8 sample 1/100 (CM31)",
                                       "fast m8 sample 1/16 (CM31)", "fast m8 sample 1/4 (CM31)", "fast m4 sample 1/16 (CM31)",
                                       "pack chunk-eval sample 1/4" };
static const uint64_t THR100 = (uint64_t)(18446744073709551615.0 / 100.0);
static const kp7pack_tab* g_tab;           // chunk/position tables for the bench params, set in main
static const kp7pack_ftab* g_ftab;         // fast 61-bit tables (same params), set in main
static const kp7pack_mtab* g_mtab;         // fast CM31 tables, set in main
static const kp7m31_params* g_p31;         // CM31 params, set in main
static const uint64_t THR4 = 1ULL << 62;
static const ntbatch::Tables* g_nt;        // batched ntHash tables, set in main
static const uint64_t THR16 = 1ULL << 60;
static const uint64_t THR1000 = (uint64_t)(18446744073709551615.0 / 1000.0);
static bool is_sample(int m)
{
  return m == PACK_S16 || m == PACK_S1000 || m == PACK_S100 || m == CHUNK_S1000 || m == CHUNK_S100 || m == CHUNK_S16 ||
         m == POS_S16 || (m >= F8_S1000 && m <= CHUNK_S4);
}

template <class Fn> static void hash_record(Method m, const std::string& s, const kp7_params* p, Fn&& fn)
{
  const unsigned k = p->k;
  if (s.size() < k) return;
  static thread_local uint64_t hbuf[4096];
  static thread_local size_t pbuf[4096];
  switch (m) {
  case NTHASH: {
    nthash::NtHash nth(s.data(), s.size(), 1, k);
    while (nth.roll()) fn(nth.get_pos(), nth.hashes()[0]);
    return;
  }
  case PACK_GATE: {
    auto tramp = [](void* ctx, size_t start, const uint64_t* h, size_t n) {
      Fn& f = *static_cast<Fn*>(ctx);
      for (size_t j = 0; j < n; j++) f(start + j, h[j]);
    };
    kp7pack_gate_seq(k, s.data(), s.size(), hbuf, 4096, tramp, &fn);
    return;
  }
  case PACK_S16:
  case PACK_S1000: {
    auto tramp = [](void* ctx, const size_t* pos, const uint64_t* h, size_t n) {
      Fn& f = *static_cast<Fn*>(ctx);
      for (size_t j = 0; j < n; j++) f(pos[j], h[j]);
    };
    kp7pack_sample(p, s.data(), s.size(), m == PACK_S16 ? THR16 : THR1000, pbuf, hbuf, 4096, tramp, &fn);
    return;
  }
  case NT_BATCH: {
    auto tramp = [](void* ctx, size_t start, const uint64_t* h, size_t n) {
      Fn& f = *static_cast<Fn*>(ctx);
      for (size_t j = 0; j < n; j++) f(start + j, h[j]);
    };
    ntbatch::hash_seq(*g_nt, s.data(), s.size(), hbuf, 4096, tramp, &fn);
    return;
  }
  case PACK_S100:
  case CHUNK_S1000:
  case CHUNK_S100:
  case CHUNK_S16:
  case POS_S16: {
    auto tramp = [](void* ctx, const size_t* pos, const uint64_t* h, size_t n) {
      Fn& f = *static_cast<Fn*>(ctx);
      for (size_t j = 0; j < n; j++) f(pos[j], h[j]);
    };
    uint64_t thr = m == PACK_S100 || m == CHUNK_S100 ? THR100 : m == CHUNK_S1000 ? THR1000 : THR16;
    int mode = m == PACK_S100 ? KP7PACK_EVAL_DIRECT : m == POS_S16 ? KP7PACK_EVAL_POS : KP7PACK_EVAL_CHUNK;
    kp7pack_sample_eval(p, g_tab, mode, s.data(), s.size(), thr, pbuf, hbuf, 4096, tramp, &fn);
    return;
  }
  case CHUNK_ALL: {
    auto tramp = [](void* ctx, size_t start, const uint64_t* h, size_t n) {
      Fn& f = *static_cast<Fn*>(ctx);
      for (size_t j = 0; j < n; j++) f(start + j, h[j]);
    };
    kp7pack_sym_all(p, g_tab, s.data(), s.size(), hbuf, 4096, tramp, &fn);
    return;
  }
  case F8_ALL: case F4_ALL: case M8_ALL: case M4_ALL: case M31_REF: {
    auto tramp = [](void* ctx, size_t start, const uint64_t* h, size_t n) {
      Fn& f = *static_cast<Fn*>(ctx);
      for (size_t j = 0; j < n; j++) f(start + j, h[j]);
    };
    if (m == M31_REF) kp7m31_sym_seq(g_p31, s.data(), s.size(), KP7M31_SCALAR, hbuf, 4096, tramp, &fn);
    else if (m == F8_ALL || m == F4_ALL) kp7pack_sym_all_f(g_ftab, m == F8_ALL ? 8 : 4, s.data(), s.size(), hbuf, 4096, tramp, &fn);
    else kp7pack_sym_all_m(g_mtab, m == M8_ALL ? 8 : 4, s.data(), s.size(), hbuf, 4096, tramp, &fn);
    return;
  }
  case F8_S1000: case F8_S100: case F8_S16: case F8_S4: case F4_S16:
  case M8_S1000: case M8_S100: case M8_S16: case M8_S4: case M4_S16: case CHUNK_S4: {
    auto tramp = [](void* ctx, const size_t* pos, const uint64_t* h, size_t n) {
      Fn& f = *static_cast<Fn*>(ctx);
      for (size_t j = 0; j < n; j++) f(pos[j], h[j]);
    };
    uint64_t thr = (m == F8_S1000 || m == M8_S1000) ? THR1000 : (m == F8_S100 || m == M8_S100) ? THR100
                 : (m == F8_S4 || m == M8_S4 || m == CHUNK_S4) ? THR4 : THR16;
    if (m == CHUNK_S4) kp7pack_sample_eval(p, g_tab, KP7PACK_EVAL_CHUNK, s.data(), s.size(), thr, pbuf, hbuf, 4096, tramp, &fn);
    else if (m >= M8_S1000) kp7pack_sample_m(g_mtab, m == M4_S16 ? 4 : 8, s.data(), s.size(), thr, pbuf, hbuf, 4096, tramp, &fn);
    else kp7pack_sample_f(g_ftab, m == F4_S16 ? 4 : 8, s.data(), s.size(), thr, pbuf, hbuf, 4096, tramp, &fn);
    return;
  }
  case KP7_SEQ_L2: {
    auto tramp = [](void* ctx, size_t start, const uint64_t* h, size_t n) {
      Fn& f = *static_cast<Fn*>(ctx);
      for (size_t j = 0; j < n; j++) f(start + j, h[j]);
    };
    kp7_sym_seq(p, s.data(), s.size(), 2, hbuf, 4096, tramp, &fn);
    return;
  }
  default: return;
  }
}

// minimizer density and per-bit frequency over runs of consecutive positions
struct Stats {
  uint64_t kmers = 0, bits[64] = { 0 };
  int ws[3] = { 10, 20, 50 };
  uint64_t windows[3] = { 0 }, sel[3] = { 0 };
};
static void run_stats(Stats& st, const std::vector<uint64_t>& h)
{
  for (uint64_t v : h) { st.kmers++; for (int b = 0; b < 64; b++) st.bits[b] += (v >> b) & 1; }
  for (int wi = 0; wi < 3; wi++) {
    size_t w = st.ws[wi];
    if (h.size() < w) continue;
    std::deque<size_t> dq; // leftmost minimum
    long last = -1;
    for (size_t i = 0; i < h.size(); i++) {
      while (!dq.empty() && h[dq.back()] > h[i]) dq.pop_back();
      dq.push_back(i);
      if (dq.front() + w <= i) dq.pop_front();
      if (i + 1 >= w) {
        st.windows[wi]++;
        if ((long)dq.front() != last) { st.sel[wi]++; last = (long)dq.front(); }
      }
    }
  }
}

int main(int argc, char** argv)
{
  if (argc < 2) { fprintf(stderr, "usage: bench_pack <fasta> [k] [rounds] [throughput] [only=<name>]\n"); return 2; }
  unsigned k = argc > 2 ? (unsigned)atoi(argv[2]) : 31;
  int rounds = argc > 3 ? atoi(argv[3]) : 5;
  bool throughput_only = argc > 4 && !strcmp(argv[4], "throughput");
  const char* only = argc > 5 && !strncmp(argv[5], "only=", 5) ? argv[5] + 5 : nullptr;
  // only= picks the method whose name equals the argument when one does, else every name containing it. A plain
  // substring once matched "sample 1/100" inside "sample 1/1000" and timed both in one counted run (2026-11-10).
  bool only_exact = false;
  if (only)
    for (int m = 0; m < NMETHODS; m++) only_exact = only_exact || !strcmp(NAMES[m], only);
  if (k > 32) { fprintf(stderr, "kp7pack needs k <= 32\n"); return 2; }
  auto t0 = std::chrono::steady_clock::now();
  std::vector<Rec> recs = load_fasta(argv[1]);
  size_t bases = 0;
  for (auto& r : recs) bases += r.seq.size();
  printf("loaded %zu records, %zu bases in %.1f s; k=%u, rounds=%d\n", recs.size(), bases,
         std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), k, rounds);

  kp7_params real;
  kp7_init(&real, k, KP7_XREAL, 0x5eed);
  static kp7pack_tab tab;
  kp7pack_tab_init(&tab, &real);
  g_tab = &tab;
  static kp7pack_ftab ftab;
  kp7pack_ftab_init(&ftab, &real);
  g_ftab = &ftab;
  static kp7m31_params p31;
  kp7m31_init(&p31, k, 0x5eed31);
  g_p31 = &p31;
  static kp7pack_mtab mtab;
  kp7pack_mtab_init(&mtab, &p31);
  g_mtab = &mtab;
  static ntbatch::Tables nt(k);
  g_nt = &nt;

  std::vector<std::vector<double>> secs(NMETHODS);
  uint64_t counts[NMETHODS] = { 0 }, msink[NMETHODS] = { 0 }, input_kmers = 0;
  for (int round = 0; round < rounds; round++) {
    for (int m = 0; m < NMETHODS; m++) {
      if (only && !(only_exact ? !strcmp(NAMES[m], only) : strstr(NAMES[m], only) != nullptr)) continue;
      uint64_t n = 0, x = 0;
      auto a = std::chrono::steady_clock::now();
      for (const Rec& r : recs) hash_record((Method)m, r.seq, &real, [&](size_t, uint64_t h) { x ^= h; n++; });
      double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
      secs[m].push_back(s);
      counts[m] = n;
      msink[m] = x;
      if (!is_sample(m) && !input_kmers) input_kmers = n;
      uint64_t per = is_sample(m) && input_kmers ? input_kmers : n;
      printf("  round %d %-38s %.3f s  %.1f M k-mers/s\n", round, NAMES[m], s, per / s / 1e6);
      fflush(stdout);
    }
  }
  printf("\nTHROUGHPUT (median of %d interleaved rounds)\n", rounds);
  for (int m = 0; m < NMETHODS; m++) {
    if (secs[m].empty()) continue;
    std::vector<double> v = secs[m];
    std::sort(v.begin(), v.end());
    double med = v[v.size() / 2];
    uint64_t per = is_sample(m) && input_kmers ? input_kmers : counts[m];
    printf("  %-38s %s %llu  median %.3f s  %.1f M k-mers/s\n", NAMES[m], is_sample(m) ? "kept" : "k-mers",
           (unsigned long long)counts[m], med, per / med / 1e6);
  }
  if (only) return 0;

  // identity: every alternative must reproduce its reference exactly (XOR of all hashes and the count)
  {
    const int pairs[][2] = { { NT_BATCH, NTHASH }, { CHUNK_S1000, PACK_S1000 }, { CHUNK_S100, PACK_S100 },
                             { CHUNK_S16, PACK_S16 }, { POS_S16, PACK_S16 }, { CHUNK_ALL, KP7_SEQ_L2 },
                             { F8_ALL, KP7_SEQ_L2 }, { F4_ALL, KP7_SEQ_L2 }, { M8_ALL, M31_REF }, { M4_ALL, M31_REF },
                             { F8_S1000, PACK_S1000 }, { F8_S100, PACK_S100 }, { F8_S16, PACK_S16 }, { F4_S16, PACK_S16 },
                             { F8_S4, CHUNK_S4 }, { M4_S16, M8_S16 } };
    bool all = true;
    for (auto& pr : pairs) {
      bool same = msink[pr[0]] == msink[pr[1]] && counts[pr[0]] == counts[pr[1]];
      printf("  %-38s identical to %s: %s\n", NAMES[pr[0]], NAMES[pr[1]], same ? "yes" : "NO");
      all = all && same;
    }
    // positions too: the first 2,000,000 (pos, hash) pairs of batched ntHash against the library, on the first record
    std::vector<std::pair<size_t, uint64_t>> a, b;
    hash_record(NTHASH, recs[0].seq, &real, [&](size_t pos, uint64_t h) { if (a.size() < 2000000) a.push_back({ pos, h }); });
    hash_record(NT_BATCH, recs[0].seq, &real, [&](size_t pos, uint64_t h) { if (b.size() < 2000000) b.push_back({ pos, h }); });
    bool pos_same = a == b;
    printf("  ntHash batched positions and hashes identical on the first %zu k-mers of %s: %s\n", a.size(),
           recs[0].name.c_str(), pos_same ? "yes" : "NO");
    if (!all || !pos_same) return 1;
  }
  // identity: the 1/16 sample's provable hashes equal kp7_sym_seq's at the gate-selected positions
  {
    uint64_t want = 0, nwant = 0;
    for (const Rec& r : recs) {
      const std::string& s = r.seq;
      hash_record(KP7_SEQ_L2, s, &real, [&](size_t pos, uint64_t h) {
        uint64_t f, rc;
        if (kp7pack_codes(s.data() + pos, k, &f, &rc) && kp7pack_gate(f, rc) < THR16) { want ^= h; nwant++; }
      });
    }
    bool same = want == msink[PACK_S16] && nwant == counts[PACK_S16];
    printf("  pack sample 1/16 provable hashes identical to kp7_sym_seq at the same positions: %s (%llu)\n",
           same ? "yes" : "NO", (unsigned long long)nwant);
    if (!same) return 1;
  }
  if (throughput_only) return 0;

  printf("\nSTATS (gate and ntHash)\n");
  for (int m : { (int)NTHASH, (int)PACK_GATE }) {
    Stats st;
    for (const Rec& r : recs) {
      std::vector<uint64_t> run;
      size_t prev = SIZE_MAX;
      hash_record((Method)m, r.seq, &real, [&](size_t pos, uint64_t h) {
        if (prev != SIZE_MAX && pos != prev + 1) { run_stats(st, run); run.clear(); }
        run.push_back(h);
        prev = pos;
      });
      run_stats(st, run);
    }
    double maxdev = 0;
    for (int b = 0; b < 64; b++) maxdev = std::max(maxdev, std::fabs((double)st.bits[b] / st.kmers - 0.5));
    printf("%s\n  k-mers %llu, max per-bit deviation from 0.5: %.5f\n", NAMES[m], (unsigned long long)st.kmers, maxdev);
    for (int wi = 0; wi < 3; wi++) {
      double d = (double)st.sel[wi] / st.windows[wi], e = 2.0 / (st.ws[wi] + 1);
      printf("  minimizer w=%-2d density %.5f (expected %.5f, ratio %.4f)\n", st.ws[wi], d, e, d / e);
    }
    fflush(stdout);
  }
  return 0;
}
