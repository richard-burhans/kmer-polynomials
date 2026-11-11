// Direct reverse-complement consistency of spaced-seed hashes: is h(window) == h(rc window)?
//
//   seedrc [n=200000] [seed=7]
//
// For each mask, hash every window of a random sequence s (n bases) and of rc(s); window p of s is the reverse
// complement of window n-k-p of rc(s). Counts the windows where the two hashes agree, for
//   ntHash2's SeedNtHash (pinned ntHash c26bd45, hashes()[0], 1 hash per seed), and
//   our Eisenstein S_m = (F_m + lam)(R_m + lam) (spaced_e.h), which is strand-symmetric for any mask by construction.
// Masks: the paper's non-palindromic 1101100111010111011100101110101 (the bench_spaced_e mask 1), two of its
// non-palindromic relatives, and palindromic masks (bench_spaced_e mask 0 and all-ones). Round-1 adjudications L1 F5
// and L7 F11 (C14) asked for this test; the reproduction attacker reported 0/199,970 for the non-palindromic mask.
// Adapted from bench/data/kmer_poly_work/round1_L7/seedrc.cpp. Output: one TSV row per mask and hash.
#include <nthash/nthash.hpp>
#include "spaced_e.h"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

static bool palindromic(const std::string& m) { return std::string(m.rbegin(), m.rend()) == m; }

int main(int argc, char** argv)
{
  const size_t n = argc > 1 ? (size_t)atoll(argv[1]) : 200000;
  const unsigned long long rs = argc > 2 ? strtoull(argv[2], nullptr, 0) : 7;
  const unsigned k = 31;
  std::mt19937_64 g(rs);
  std::string s(n, 'A');
  for (char& c : s) c = "ACGT"[g() & 3];
  std::string r(s.rbegin(), s.rend());
  for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
  const char* masks[] = { "1101100111010111011100101110101", "1101100111010111011100101110111",
                          "1111111111111110000000000000001", "1110111011101110111011101110111",
                          "1111111111111111111111111111111" };
  kp7e_params kpe;
  kp7e_init(&kpe, k, 0x5eed, 0x1a3d);
  printf("# seedrc n=%zu random seed=%llu k=%u; window p of s vs window n-k-p of rc(s)\n", n, rs, k);
  printf("mask\tpalindromic\tweight\thash\twindows\th(w)==h(rc w)\n");
  for (const char* mk : masks) {
    const size_t nw = n - k + 1;
    int weight = 0;
    for (const char* c = mk; *c; c++) weight += *c == '1';
    std::vector<uint64_t> a(n, 0), b(n, 0);
    std::vector<std::string> sd(1, mk);
    { nthash::SeedNtHash h(s.data(), s.size(), sd, 1, k); while (h.roll()) a[h.get_pos()] = h.hashes()[0]; }
    { nthash::SeedNtHash h(r.data(), r.size(), sd, 1, k); while (h.roll()) b[h.get_pos()] = h.hashes()[0]; }
    size_t eq = 0;
    for (size_t p = 0; p < nw; p++) eq += a[p] == b[nw - 1 - p];
    printf("%s\t%s\t%d\tSeedNtHash\t%zu\t%zu\n", mk, palindromic(mk) ? "yes" : "no", weight, nw, eq);

    spe_params pe;
    spe_init(&pe, &kpe, mk);
    std::fill(a.begin(), a.end(), 0);
    std::fill(b.begin(), b.end(), 0);
    spe_roller ro;
    uint64_t h;
    spe_start(&ro, &pe, s.data(), s.size());
    while (spe_roll(&ro, &h)) a[ro.pos] = h;
    spe_start(&ro, &pe, r.data(), r.size());
    while (spe_roll(&ro, &h)) b[ro.pos] = h;
    eq = 0;
    for (size_t p = 0; p < nw; p++) eq += a[p] == b[nw - 1 - p];
    printf("%s\t%s\t%d\tEisenstein S_m\t%zu\t%zu\n", mk, palindromic(mk) ? "yes" : "no", weight, nw, eq);
  }
  return 0;
}
