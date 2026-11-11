/* KP-11, the form that does work: k-mer mismatch-neighbourhood counts by convolution over the group (Z_4)^k.
 *
 * Why not the ordinary generating function: with F_S(y) = sum y^code(s), F_U(y^{4^|v|}) F_V(y) = sum y^code(uv), and
 * code(uv) = code(u) 4^|v| + code(v) is INJECTIVE in (u, v). No two terms ever share an exponent, so the product is
 * just the Cartesian product U x V written as a polynomial: it aggregates nothing, and an FFT cannot beat listing the
 * pairs (results/ideas_kp5_kp11_*.md).
 *
 * The group-algebra version does aggregate. Treat a k-mer as a vector in (Z_4)^k (2-bit digits). Hamming distance is
 * the weight of the difference vector, since t_j - s_j = 0 mod 4 iff t_j = s_j. So for a k-mer count table C,
 *   N_d(t) = #{ occurrences s : H(s, t) <= d } = (C * B_d)(t),   B_d = indicator of weight <= d,
 * a convolution over (Z_4)^k. Its characters are chi_m(a) = i^{a.m}: the fourth roots of unity of the original KP
 * encoding. B_d's transform depends only on wt(m) and is a sum of Krawtchouk polynomials (q = 4):
 *   hat B_d(m) = sum_{j<=d} K_j(wt m),   K_j(w) = sum_i (-1)^i 3^{j-i} C(w, i) C(k - w, j - i)
 * so N_d = IDFT( DFT(C) . K ) costs O(k 4^k) for EVERY target and ANY d. Direct ball enumeration costs
 * (#distinct k-mers) x |ball(k, d)|.
 *
 * Exactness: the transform runs mod P = 998244353 (P = 1 mod 4, so i exists: I4 = 3^((P-1)/4)). Every true N_d(t) is
 * at most the number of k-mer occurrences, which must be < P (checked), so the residue IS the integer.
 *
 *   mismatch FASTA RECORD K D fft|direct|both     RECORD "-" = first record; K <= 13; D <= 3
 *   mismatch FASTA RECORD K D check < CODES        fft result checked against "code expected" lines on stdin
 * prints a checksum over all 4^K targets and timings; `both` also checks the two arrays are identical.
 */
#include "fasta.h"
#include <stdint.h>
#include <time.h>

#define P 998244353ULL
static uint64_t pw(uint64_t b, uint64_t e) { uint64_t r = 1; b %= P; while (e) { if (e & 1) r = r * b % P; b = b * b % P; e >>= 1; } return r; }

static int8_t CODE[256];

static uint64_t ball_size(unsigned k, unsigned d)
{
  uint64_t s = 0, c = 1, p3 = 1;
  for (unsigned j = 0; j <= d; j++) { s += c * p3; c = c * (k - j) / (j + 1); p3 *= 3; }
  return s;
}

/* count every valid k-mer occurrence (forward strand) */
static uint32_t* count_kmers(const char* s, size_t len, unsigned k, uint64_t* occ)
{
  size_t n = (size_t)1 << (2 * k);
  uint32_t* C = calloc(n, sizeof *C);
  uint64_t mask = n - 1, code = 0;
  size_t run = 0;
  *occ = 0;
  for (size_t j = 0; j < len; j++) {
    int c = CODE[(unsigned char)s[j]];
    if (c < 0) { run = 0; code = 0; continue; }
    code = ((code << 2) | (unsigned)c) & mask;
    if (++run >= k) { C[code]++; (*occ)++; }
  }
  return C;
}

static void direct(const uint32_t* C, unsigned k, unsigned d, uint64_t* N)
{
  size_t n = (size_t)1 << (2 * k);
  memset(N, 0, n * sizeof *N);
  for (size_t s = 0; s < n; s++) {
    uint32_t c = C[s];
    if (!c) continue;
    N[s] += c;
    if (d < 1) continue;
    for (unsigned p1 = 0; p1 < k; p1++)
      for (unsigned a1 = 1; a1 < 4; a1++) {
        size_t t1 = s ^ ((size_t)a1 << (2 * p1)); /* XOR with 1..3 visits the three other bases at p1 */
        N[t1] += c;
        if (d < 2) continue;
        for (unsigned p2 = p1 + 1; p2 < k; p2++)
          for (unsigned a2 = 1; a2 < 4; a2++) {
            size_t t2 = t1 ^ ((size_t)a2 << (2 * p2));
            N[t2] += c;
            if (d < 3) continue;
            for (unsigned p3 = p2 + 1; p3 < k; p3++)
              for (unsigned a3 = 1; a3 < 4; a3++) N[t2 ^ ((size_t)a3 << (2 * p3))] += c;
          }
      }
  }
}

/* in-place DFT over (Z_4)^k mod P with characters i^{a.m} (inv = 0) or i^{-a.m} (inv = 1) */
static void dft4k(uint64_t* x, unsigned k, int inv, uint64_t I4)
{
  size_t n = (size_t)1 << (2 * k);
  uint64_t w = inv ? (P - I4) % P : I4; /* i or i^{-1} = -i */
  for (unsigned j = 0; j < k; j++) {
    size_t st = (size_t)1 << (2 * j);
    for (size_t base = 0; base < n; base += 4 * st)
      for (size_t r = 0; r < st; r++) {
        uint64_t* a = x + base + r;
        uint64_t x0 = a[0], x1 = a[st], x2 = a[2 * st], x3 = a[3 * st];
        uint64_t wx1 = w * x1 % P, wx3 = w * x3 % P;
        a[0] = (x0 + x1 + x2 + x3) % P;
        a[st] = (x0 + wx1 + (P - x2) + (P - wx3)) % P;      /* X1 = x0 + w x1 - x2 - w x3 */
        a[2 * st] = (x0 + (P - x1) + x2 + (P - x3)) % P;   /* X2 = x0 - x1 + x2 - x3 */
        a[3 * st] = (x0 + (P - wx1) + (P - x2) + wx3) % P; /* X3 = x0 - w x1 - x2 + w x3 */
      }
  }
}

static int64_t binom(int n, int r) { if (r < 0 || r > n) return 0; int64_t c = 1; for (int i = 0; i < r; i++) c = c * (n - i) / (i + 1); return c; }

static void via_fft(const uint32_t* C, unsigned k, unsigned d, uint64_t* N)
{
  size_t n = (size_t)1 << (2 * k);
  uint64_t I4 = pw(3, (P - 1) / 4); /* 3 is a primitive root mod P */
  for (size_t s = 0; s < n; s++) N[s] = C[s];
  dft4k(N, k, 0, I4);
  /* hat B_d by weight w = number of nonzero digits of m */
  uint64_t Kw[64];
  for (unsigned w = 0; w <= k; w++) {
    int64_t tot = 0, p3 = 1;
    for (unsigned j = 0; j <= d; j++) {
      int64_t kj = 0;
      for (unsigned i = 0; i <= j; i++) {
        int64_t term = binom((int)w, (int)i) * binom((int)(k - w), (int)(j - i));
        int64_t p = 1; for (unsigned e = 0; e < j - i; e++) p *= 3;
#ifdef KP11_MUTATE /* test-only: drop the Krawtchouk sign */
        kj += term * p;
#else
        kj += (i & 1 ? -term : term) * p;
#endif
      }
      tot += kj; (void)p3;
    }
    Kw[w] = (uint64_t)((tot % (int64_t)P + (int64_t)P) % (int64_t)P);
  }
  for (size_t m = 0; m < n; m++) {
    uint64_t nz = (m | (m >> 1)) & 0x5555555555555555ULL;
    N[m] = N[m] * Kw[__builtin_popcountll(nz)] % P;
  }
  dft4k(N, k, 1, I4);
  uint64_t invn = pw(pw(4, k), P - 2);
  for (size_t m = 0; m < n; m++) N[m] = N[m] * invn % P;
}

static inline uint64_t mix64(uint64_t x) { x ^= x >> 31; x *= 0x7fb5d329728ea185ULL; x ^= x >> 27; x *= 0x81dadef4bc2dd44dULL; return x ^ (x >> 33); }
static uint64_t checksum(const uint64_t* N, size_t n) { uint64_t h = 0; for (size_t t = 0; t < n; t++) h += mix64(t) * (N[t] + 1); return h; }

int main(int argc, char** argv)
{
  if (argc < 6) { fprintf(stderr, "usage: mismatch FASTA RECORD K D fft|direct|both\n"); return 2; }
  memset(CODE, -1, sizeof CODE);
  CODE['A'] = 0; CODE['C'] = 1; CODE['G'] = 2; CODE['T'] = 3;
  size_t len;
  char* s = load_record(argv[1], strcmp(argv[2], "-") ? argv[2] : NULL, &len);
  unsigned k = (unsigned)atoi(argv[3]), d = (unsigned)atoi(argv[4]);
  const char* m = argv[5];
  if (k < 1 || k > 13 || d > 3 || d > k) { fprintf(stderr, "need 1 <= K <= 13, D <= 3\n"); return 2; }
  size_t n = (size_t)1 << (2 * k);
  uint64_t occ;
  double t0 = now_s();
  uint32_t* C = count_kmers(s, len, k, &occ);
  double t1 = now_s();
  if (occ >= P) { fprintf(stderr, "occurrences %llu >= P: residues would not be exact\n", (unsigned long long)occ); return 3; }
  size_t distinct = 0; for (size_t t = 0; t < n; t++) distinct += C[t] != 0;
  uint64_t* A = NULL; uint64_t* B = NULL;
  double tf = 0, td = 0;
  if (!strcmp(m, "check")) {
    A = malloc(n * 8);
    via_fft(C, k, d, A);
    unsigned long long code, want;
    int bad = 0, seen = 0;
    while (scanf("%llu %llu", &code, &want) == 2) { seen++; if (code >= n || A[code] != want) { printf("BAD %llu got %llu want %llu\n", code, code < n ? (unsigned long long)A[code] : 0ULL, want); bad = 1; } }
    printf("checked %d targets%s\n", seen, bad ? " with mismatches" : "");
    return bad || !seen;
  }
  if (!strcmp(m, "fft") || !strcmp(m, "both")) { A = malloc(n * 8); double a = now_s(); via_fft(C, k, d, A); tf = now_s() - a; }
  if (!strcmp(m, "direct") || !strcmp(m, "both")) { B = malloc(n * 8); double a = now_s(); direct(C, k, d, B); td = now_s() - a; }
  printf("k=%u d=%u len=%zu occurrences=%llu distinct=%zu targets=%zu ball=%llu count_s=%.3f", k, d, len,
         (unsigned long long)occ, distinct, n, (unsigned long long)ball_size(k, d), t1 - t0);
  if (A) printf(" fft_s=%.3f fft_sum=%016llx", tf, (unsigned long long)checksum(A, n));
  if (B) printf(" direct_s=%.3f direct_sum=%016llx", td, (unsigned long long)checksum(B, n));
  if (A && B) {
    size_t bad = 0; for (size_t t = 0; t < n; t++) bad += A[t] != B[t];
    printf(" identical=%s", bad ? "NO" : "yes");
    if (bad) { printf("\n"); return 1; }
  }
  printf("\n");
  return 0;
}
