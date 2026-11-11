// Full-genome, element-by-element identity check: the batched ntHash (hash/ntbatch.h) against the ntHash library's
// canonical hash (nthash::NtHash(seq, len, 1, k).roll(), hashes()[0]), position AND value at every k-mer.
//
//   nteq <fasta> [k=31]
//
// Adapted from the round-1 reproduction attacker's nteq.cpp (bench/data/kmer_poly_work/round1_L7). Differences: the
// sequence is NOT uppercased (the benchmarks hash the soft-masked FASTA as written, so case handling is part of what
// is checked), k is a parameter, and the output lists per-record counts. bench_pack.cpp's older check compared an XOR
// checksum over all k-mers and positions only on the first 2,000,000 (round-1 adjudications L1 F11, L4 F4, L7 F3).
// Exit 0 iff every record has the same k-mer count, positions and hashes in both.
#include <nthash/nthash.hpp>
#include "ntbatch.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

int main(int argc, char** argv)
{
  if (argc < 2) { fprintf(stderr, "usage: nteq <fasta> [k=31]\n"); return 2; }
  unsigned k = argc > 2 ? (unsigned)atoi(argv[2]) : 31;
  std::ifstream in(argv[1]);
  if (!in) { perror(argv[1]); return 2; }
  std::string line, seq, name;
  std::vector<std::pair<std::string, std::string>> recs;
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] == '>') {
      if (!name.empty()) recs.push_back({ name, seq });
      name = line.substr(1, line.find_first_of(" \t") - 1);
      seq.clear();
    } else {
      seq += line;
    }
  }
  if (!name.empty()) recs.push_back({ name, seq });
  ntbatch::Tables T(k);
  static uint64_t buf[4096];
  unsigned long long n = 0, nlib = 0, bad = 0, bases = 0;
  printf("record\tbases\tbatched_kmers\tlibrary_kmers\tmismatches\n");
  for (auto& r : recs) {
    const std::string& s = r.second;
    bases += s.size();
    unsigned long long rb = 0, rl = 0, rbad = 0;
    if (s.size() >= k) {
      std::vector<std::pair<size_t, uint64_t>> b;
      ntbatch::hash_seq(T, s.data(), s.size(), buf, 4096,
                        [](void* ctx, size_t start, const uint64_t* h, size_t m) {
                          auto v = (std::vector<std::pair<size_t, uint64_t>>*)ctx;
                          for (size_t j = 0; j < m; j++) v->push_back({ start + j, h[j] });
                        },
                        &b);
      nthash::NtHash nth(s.data(), s.size(), 1, k);
      size_t i = 0;
      while (nth.roll()) {
        rl++;
        if (i >= b.size() || b[i].first != nth.get_pos() || b[i].second != nth.hashes()[0]) {
          if (bad + rbad < 5) printf("# mismatch in %s at %zu\n", r.first.c_str(), (size_t)nth.get_pos());
          rbad++;
        }
        i++;
      }
      if (i != b.size()) { printf("# count differs in %s: library %zu batched %zu\n", r.first.c_str(), i, b.size()); rbad++; }
      rb = b.size();
    }
    printf("%s\t%zu\t%llu\t%llu\t%llu\n", r.first.c_str(), s.size(), rb, rl, rbad);
    n += rb; nlib += rl; bad += rbad;
  }
  printf("TOTAL records %zu  bases %llu  k %u  batched k-mers %llu  library k-mers %llu  mismatches %llu\n", recs.size(),
         bases, k, n, nlib, bad);
  return bad != 0;
}
