// KP-14 Eisenstein revisit: spaced seeds over the Eisenstein hash (spaced_e.h) against ntHash2's SeedNtHash and the
// Gaussian spaced roller (ideas/algebra/spaced.h): throughput, interleaved rounds, and sensitivity to substitutions.
// Adapted from ideas/algebra/bench_spaced.cpp.
//   bench_spaced_e <fasta> [max_bases=60000000] [rounds=3]
// Same masks for both, one hash per seed, single thread, same -O3 -march=native build (Makefile).
// Sensitivity: a copy of the sequence with substitutions at rate r (no indels, so windows stay aligned); for each
// method the fraction of windows whose hash is unchanged, against the exact fraction whose masked double-stranded
// content is unchanged ({care bases of s, care bases of rc(s)} as a set).
#include "../algebra/spaced.h"
#include "spaced_e.h"
#include <nthash/nthash.hpp>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

static std::string load(const char* path, size_t maxb)
{
  std::ifstream in(path);
  std::string line, s;
  while (std::getline(in, line) && s.size() < maxb) {
    if (line.empty() || line[0] == '>') continue;
    for (char& c : line) c = (char)toupper((unsigned char)c);
    s += line;
  }
  if (s.size() > maxb) s.resize(maxb);
  return s;
}
static uint64_t sm(uint64_t* s)
{
  uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static const char* MASKS[] = { "1110111011101110111011101110111", "1101100111010111011100101110101" };

// ours: all masks for every window; fn(pos, mask_index, hash)
template <class Fn> static void ours(const std::string& s, const std::vector<sp_params>& ps, Fn&& fn)
{
  std::vector<sp_roller> rs(ps.size());
  for (size_t i = 0; i < ps.size(); i++) sp_start(&rs[i], &ps[i], s.data(), s.size());
  uint64_t h;
  for (;;) {
    bool any = false;
    for (size_t i = 0; i < ps.size(); i++) {
      if (!sp_roll(&rs[i], &h)) continue;
      any = true;
      fn(rs[i].pos, i, h);
    }
    if (!any) break;
  }
}
template <class Fn> static void ours_e(const std::string& s, const std::vector<spe_params>& ps, Fn&& fn)
{
  std::vector<spe_roller> rs(ps.size());
  for (size_t i = 0; i < ps.size(); i++) spe_start(&rs[i], &ps[i], s.data(), s.size());
  uint64_t h;
  for (;;) {
    bool any = false;
    for (size_t i = 0; i < ps.size(); i++) {
      if (!spe_roll(&rs[i], &h)) continue;
      any = true;
      fn(rs[i].pos, i, h);
    }
    if (!any) break;
  }
}
template <class Fn> static void nt(const std::string& s, const std::vector<std::string>& seeds, Fn&& fn)
{
  nthash::SeedNtHash h(s.data(), s.size(), seeds, 1, 31);
  while (h.roll())
    for (size_t i = 0; i < seeds.size(); i++) fn(h.get_pos(), i, h.hashes()[i]);
}

int main(int argc, char** argv)
{
  if (argc < 2) { fprintf(stderr, "usage: bench_spaced <fasta> [max_bases] [rounds]\n"); return 2; }
  size_t maxb = argc > 2 ? strtoull(argv[2], 0, 10) : 60000000;
  int rounds = argc > 3 ? atoi(argv[3]) : 3;
  std::string s = load(argv[1], maxb);
  kp7_params kp;
  kp7_init(&kp, 31, KP7_XREAL, 0x5eed);
  kp7e_params kpe;
  kp7e_init(&kpe, 31, 0x5eed);
  printf("bases %zu, rounds %d\n", s.size(), rounds);
  for (int nmask = 1; nmask <= 2; nmask++) {
    std::vector<sp_params> ps(nmask);
    std::vector<spe_params> pe(nmask);
    std::vector<std::string> seeds;
    for (int i = 0; i < nmask; i++) { sp_init(&ps[i], &kp, MASKS[i]); spe_init(&pe[i], &kpe, MASKS[i]); seeds.push_back(MASKS[i]); }
    std::vector<double> to, te, tn;
    uint64_t no = 0, ne = 0, nn = 0, xo = 0, xe = 0, xn = 0;
    for (int r = 0; r < rounds; r++) {
      auto a = std::chrono::steady_clock::now();
      no = 0; xo = 0;
      ours(s, ps, [&](size_t, size_t, uint64_t h) { xo ^= h; no++; });
      auto b = std::chrono::steady_clock::now();
      ne = 0; xe = 0;
      ours_e(s, pe, [&](size_t, size_t, uint64_t h) { xe ^= h; ne++; });
      auto c = std::chrono::steady_clock::now();
      nn = 0; xn = 0;
      nt(s, seeds, [&](size_t, size_t, uint64_t h) { xn ^= h; nn++; });
      auto d = std::chrono::steady_clock::now();
      to.push_back(std::chrono::duration<double>(b - a).count());
      te.push_back(std::chrono::duration<double>(c - b).count());
      tn.push_back(std::chrono::duration<double>(d - c).count());
    }
    std::sort(to.begin(), to.end()); std::sort(te.begin(), te.end()); std::sort(tn.begin(), tn.end());
    double mo = to[to.size() / 2], me = te[te.size() / 2], mn = tn[tn.size() / 2];
    printf("masks %d: Gaussian %llu hashes %.3f s (%.1f M/s)  Eisenstein %llu %.3f s (%.1f M/s)  SeedNtHash %llu %.3f s (%.1f M/s)"
           "  Gaussian/nt %.2fx  Eisenstein/nt %.2fx  (sinks %016llx %016llx %016llx)\n",
           nmask, (unsigned long long)no, mo, no / mo / 1e6, (unsigned long long)ne, me, ne / me / 1e6,
           (unsigned long long)nn, mn, nn / mn / 1e6, mn / mo, mn / me,
           (unsigned long long)xo, (unsigned long long)xe, (unsigned long long)xn);
  }
  // sensitivity
  std::vector<sp_params> ps(2);
  std::vector<std::string> seeds;
  for (int i = 0; i < 2; i++) { sp_init(&ps[i], &kp, MASKS[i]); seeds.push_back(MASKS[i]); }
  const double rates[] = { 0.001, 0.01, 0.05 };
  for (double rate : rates) {
    std::string m = s;
    uint64_t st = 99;
    for (char& c : m) if (c != 'N' && (sm(&st) % 1000000) < rate * 1e6) { const char* alt = "ACGT"; char nc; do nc = alt[sm(&st) % 4]; while (nc == c); c = nc; }
    for (int mi = 0; mi < 2; mi++) {
      std::vector<sp_params> one(1, ps[mi]);
      std::vector<std::string> sd(1, seeds[mi]);
      std::vector<uint64_t> ho(s.size(), 0), hm(s.size(), 0), no(s.size(), 0), nmu(s.size(), 0);
      std::vector<uint64_t> he(s.size(), 0), hem(s.size(), 0);
      std::vector<spe_params> onee(1);
      spe_init(&onee[0], &kpe, MASKS[mi]);
      std::vector<char> v(s.size(), 0);
      ours(s, one, [&](size_t p, size_t, uint64_t h) { ho[p] = h; v[p] |= 1; });
      ours(m, one, [&](size_t p, size_t, uint64_t h) { hm[p] = h; v[p] |= 2; });
      nt(s, sd, [&](size_t p, size_t, uint64_t h) { no[p] = h; v[p] |= 4; });
      nt(m, sd, [&](size_t p, size_t, uint64_t h) { nmu[p] = h; v[p] |= 8; });
      ours_e(s, onee, [&](size_t p, size_t, uint64_t h) { he[p] = h; v[p] |= 16; });
      ours_e(m, onee, [&](size_t p, size_t, uint64_t h) { hem[p] = h; v[p] |= 32; });
      uint64_t tot = 0, so = 0, sn = 0, se = 0, see = 0;
      const char* mk = MASKS[mi];
      for (size_t p = 0; p + 31 <= s.size(); p++) {
        if (v[p] != 63) continue;
        tot++;
        if (ho[p] == hm[p]) so++;
        if (no[p] == nmu[p]) sn++;
        if (he[p] == hem[p]) see++;
        // exact: care bases of s and of rc(s) unchanged as a set
        bool fw = true, rv = true, sw1 = true, sw2 = true;
        for (int t = 0; t < 31; t++) if (mk[t] == '1') {
          if (s[p + t] != m[p + t]) fw = false;
          if (s[p + 30 - t] != m[p + 30 - t]) rv = false; // rc(s) care base t is comp(s[p+30-t])
        }
        (void)sw1; (void)sw2;
        if (fw && rv) se++;
      }
      printf("SNP rate %.3f mask %d: windows %llu, unchanged: Gaussian %.5f  Eisenstein %.5f  SeedNtHash %.5f  exact(care positions of both strands) %.5f\n",
             rate, mi, (unsigned long long)tot, (double)so / tot, (double)see / tot, (double)sn / tot, (double)se / tot);
    }
  }
  return 0;
}
