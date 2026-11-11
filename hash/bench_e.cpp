// Eisenstein KP-7 (kp7e) against batched ntHash, the 61-bit Gaussian kp7 and CM31: throughput per INPUT k-mer.
//
//   bench_e <fasta> [k=31] [rounds=5] [only=<name>]
//     only= times just the method whose name equals the argument (else every name containing it), for laila_perf
//     counter runs; only=__load_only__ times nothing (the load baseline).
//
// The FASTA is loaded once and upper-cased; every method sees identical bytes. `rounds` interleaved passes, timing
// only the hashing loop, median reported. Sample methods emit only kept k-mers, so their rate is per input k-mer.
// Identity (full runs only): kp7e sym_seq L1/L2 must equal the kp7e scalar roller (XOR and count); the branchy and
// batched kp7e samplers must agree with each other and with the roller's S at the gate-selected positions.
#include "kp7e.h"
#include "kp7pack.h"
#include "ntbatch.h"
#include <nthash/nthash.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
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

enum Method { NT_BATCH = 0, NT_LIB, E_REF, E_L1, E_L2, G_L1, G_L2, M31_SIMD,
              EB_S1000, EB_S100, EB_S16, EB_S4, EX_S1000, EX_S100, EX_S16, EX_S4,
              MB_S1000, MB_S100, MB_S16, MB_S4, M8_S16, M8_S4, NMETHODS };
static const char* NAMES[NMETHODS] = {
  "ntHash batched (ntbatch)", "ntHash library roll", "e scalar roller (reference)", "e sym_seq L1", "e sym_seq L2",
  "gauss61 sym_seq L1", "gauss61 sym_seq L2", "m31 simd (avx512 else avx2)",
  "e branchy 1/1000", "e branchy 1/100", "e branchy 1/16", "e branchy 1/4",
  "e batched 1/1000", "e batched 1/100", "e batched 1/16", "e batched 1/4",
  "m31 branchy 1/1000", "m31 branchy 1/100", "m31 branchy 1/16", "m31 branchy 1/4",
  "m31 batched 1/16", "m31 batched 1/4" };
static const uint64_t THR1000 = (uint64_t)(18446744073709551615.0 / 1000.0);
static const uint64_t THR100 = (uint64_t)(18446744073709551615.0 / 100.0);
static const uint64_t THR16 = 1ULL << 60, THR4 = 1ULL << 62;
static bool is_sample(int m) { return m >= EB_S1000; }
static uint64_t thr_of(int m)
{
  switch (m) {
  case EB_S1000: case EX_S1000: case MB_S1000: return THR1000;
  case EB_S100: case EX_S100: case MB_S100: return THR100;
  case EB_S16: case EX_S16: case MB_S16: case M8_S16: return THR16;
  default: return THR4;
  }
}

static const kp7e_params* g_pe;
static const kp7e_tab* g_te;
static const kp7_params* g_pg;
static const kp7m31_params* g_p31;
static const kp7pack_mtab* g_mtab;
static const ntbatch::Tables* g_nt;
static int g_m31_impl;

template <class Fn> static void hash_record(Method m, const std::string& s, unsigned k, Fn&& fn)
{
  if (s.size() < k) return;
  static thread_local uint64_t hbuf[4096];
  static thread_local size_t pbuf[4096];
  auto seq_tramp = [](void* ctx, size_t start, const uint64_t* h, size_t n) {
    Fn& f = *static_cast<Fn*>(ctx);
    for (size_t j = 0; j < n; j++) f(start + j, h[j]);
  };
  auto pair_tramp = [](void* ctx, const size_t* pos, const uint64_t* h, size_t n) {
    Fn& f = *static_cast<Fn*>(ctx);
    for (size_t j = 0; j < n; j++) f(pos[j], h[j]);
  };
  switch (m) {
  case NT_BATCH: ntbatch::hash_seq(*g_nt, s.data(), s.size(), hbuf, 4096, seq_tramp, &fn); return;
  case NT_LIB: {
    nthash::NtHash nth(s.data(), s.size(), 1, k);
    while (nth.roll()) fn(nth.get_pos(), nth.hashes()[0]);
    return;
  }
  case E_REF: {
    kp7e_roller r;
    kp7e_start(&r, g_pe, s.data(), s.size());
    while (kp7e_roll(&r)) fn(r.pos, kp7e_sym(g_pe, r.F, r.R));
    return;
  }
  case E_L1: case E_L2: kp7e_sym_seq(g_pe, s.data(), s.size(), m == E_L2 ? 2 : 1, hbuf, 4096, seq_tramp, &fn); return;
  case G_L1: case G_L2: kp7_sym_seq(g_pg, s.data(), s.size(), m == G_L2 ? 2 : 1, hbuf, 4096, seq_tramp, &fn); return;
  case M31_SIMD: kp7m31_sym_seq(g_p31, s.data(), s.size(), g_m31_impl, hbuf, 4096, seq_tramp, &fn); return;
  case EB_S1000: case EB_S100: case EB_S16: case EB_S4:
    kp7e_sample_b(g_te, s.data(), s.size(), thr_of(m), pbuf, hbuf, 4096, pair_tramp, &fn); return;
  case EX_S1000: case EX_S100: case EX_S16: case EX_S4:
    kp7e_sample_batch(g_te, s.data(), s.size(), thr_of(m), pbuf, hbuf, 4096, pair_tramp, &fn); return;
  case MB_S1000: case MB_S100: case MB_S16: case MB_S4:
    kp7pack_sample_mb(g_mtab, s.data(), s.size(), thr_of(m), pbuf, hbuf, 4096, pair_tramp, &fn); return;
  case M8_S16: case M8_S4:
    kp7pack_sample_m(g_mtab, 8, s.data(), s.size(), thr_of(m), pbuf, hbuf, 4096, pair_tramp, &fn); return;
  default: return;
  }
}

int main(int argc, char** argv)
{
  if (argc < 2) { fprintf(stderr, "usage: bench_e <fasta> [k] [rounds] [only=<name>]\n"); return 2; }
  unsigned k = argc > 2 ? (unsigned)atoi(argv[2]) : 31;
  int rounds = argc > 3 ? atoi(argv[3]) : 5;
  const char* only = argc > 4 && !strncmp(argv[4], "only=", 5) ? argv[4] + 5 : nullptr;
  bool only_exact = false;
  if (only) for (int m = 0; m < NMETHODS; m++) only_exact = only_exact || !strcmp(NAMES[m], only);
  if (k > 32) { fprintf(stderr, "the samplers need k <= 32\n"); return 2; }
  auto t0 = std::chrono::steady_clock::now();
  std::vector<Rec> recs = load_fasta(argv[1]);
  size_t bases = 0;
  for (auto& r : recs) bases += r.seq.size();
  printf("loaded %zu records, %zu bases in %.1f s; k=%u, rounds=%d\n", recs.size(), bases,
         std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), k, rounds);

  static kp7e_params pe; kp7e_init(&pe, k, 0x5eedE, 0x5eedE + 1); g_pe = &pe;
  static kp7e_tab te; kp7e_tab_init(&te, &pe); g_te = &te;
  static kp7_params pg; kp7_init(&pg, k, KP7_XREAL, 0x5eed, 0x5eed + 1); g_pg = &pg;
  static kp7m31_params p31; kp7m31_init(&p31, k, 0x5eed31, 0x5eed31 + 1); g_p31 = &p31;
  static kp7pack_mtab mtab; kp7pack_mtab_init(&mtab, &p31); g_mtab = &mtab;
  static ntbatch::Tables nt(k); g_nt = &nt;
  g_m31_impl = kp7m31_have(KP7M31_AVX512) ? KP7M31_AVX512 : kp7m31_have(KP7M31_AVX2) ? KP7M31_AVX2 : KP7M31_SCALAR;
  printf("m31 simd impl: %s\n", g_m31_impl == KP7M31_AVX512 ? "AVX-512" : g_m31_impl == KP7M31_AVX2 ? "AVX2" : "scalar");

  std::vector<std::vector<double>> secs(NMETHODS);
  uint64_t counts[NMETHODS] = { 0 }, msink[NMETHODS] = { 0 }, input_kmers = 0;
  for (int round = 0; round < rounds; round++)
    for (int m = 0; m < NMETHODS; m++) {
      if (only && !(only_exact ? !strcmp(NAMES[m], only) : strstr(NAMES[m], only) != nullptr)) continue;
      uint64_t n = 0, x = 0;
      auto a = std::chrono::steady_clock::now();
      for (const Rec& r : recs) hash_record((Method)m, r.seq, k, [&](size_t, uint64_t h) { x ^= h; n++; });
      double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
      secs[m].push_back(s); counts[m] = n; msink[m] = x;
      if (!is_sample(m) && !input_kmers) input_kmers = n;
      uint64_t per = is_sample(m) && input_kmers ? input_kmers : n;
      printf("  round %d %-30s %.3f s  %.1f M k-mers/s\n", round, NAMES[m], s, per / s / 1e6);
      fflush(stdout);
    }
  printf("\nTHROUGHPUT (median of %d interleaved rounds)\n", rounds);
  for (int m = 0; m < NMETHODS; m++) {
    if (secs[m].empty()) continue;
    std::vector<double> v = secs[m];
    std::sort(v.begin(), v.end());
    double med = v[v.size() / 2];
    uint64_t per = is_sample(m) && input_kmers ? input_kmers : counts[m];
    printf("  %-30s %s %llu  median %.3f s  %.1f M k-mers/s\n", NAMES[m], is_sample(m) ? "kept" : "k-mers",
           (unsigned long long)counts[m], med, per / med / 1e6);
  }
  if (only) return 0;

  bool all = true;
  const int pairs[][2] = { { E_L1, E_REF }, { E_L2, E_REF }, { G_L2, G_L1 }, { EX_S1000, EB_S1000 }, { EX_S100, EB_S100 },
                           { EX_S16, EB_S16 }, { EX_S4, EB_S4 }, { M8_S16, MB_S16 }, { M8_S4, MB_S4 } };
  for (auto& pr : pairs) {
    bool same = msink[pr[0]] == msink[pr[1]] && counts[pr[0]] == counts[pr[1]];
    printf("  %-30s identical to %s: %s\n", NAMES[pr[0]], NAMES[pr[1]], same ? "yes" : "NO");
    all = all && same;
  }
  // the e samplers' hashes equal the reference roller's S at the gate-selected positions (1/16)
  uint64_t want = 0, nwant = 0;
  for (const Rec& r : recs) {
    const std::string& s = r.seq;
    hash_record(E_REF, s, k, [&](size_t pos, uint64_t h) {
      uint64_t f, rc;
      if (kp7e_codes(s.data() + pos, k, &f, &rc) && kp7e_gate(f, rc) < THR16) { want ^= h; nwant++; }
    });
  }
  bool same = want == msink[EB_S16] && nwant == counts[EB_S16];
  printf("  e branchy 1/16 hashes identical to the reference roller at the gate positions: %s (%llu)\n",
         same ? "yes" : "NO", (unsigned long long)nwant);
  return all && same ? 0 : 1;
}
