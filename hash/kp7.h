/* KP-6 / KP-7: rolling polynomial k-mer hashes over F_{q^2} = Z[i]/(q), q = 2^61 - 1.
 *
 * Encoding A = 1, C = i, G = -i, T = -1, so complement is negation and
 *   F(s)  = sum_t s_t x^t                      (forward)
 *   R(s)  = rc(p_s)(x) = -x^{k-1} p_s(1/x)     (reverse complement)
 *   S(s)  = (F + lam)(R + lam)                 (strand-symmetric; notes/kp7_strand_symmetric_hash.md)
 *
 * Rolling, both strands O(1) per base:
 *   F' = F * x^{-1} - s_out * x^{-1} + s_in * x^{k-1}
 *   R' = R * x      + s_out * x^{k}  - s_in
 *
 * Two variants of the evaluation point:
 *   KP7_XREAL     x in F_q (a real field element). Multiplying an F_{q^2} element by x costs 2 modular
 *                 multiplications instead of 4. Distinct k-mers collide with probability <= (k-1)/q
 *                 per pair (a nonzero polynomial with F_{q^2} coefficients has <= k-1 roots in F_q).
 *   KP7_XCOMPLEX  x in F_{q^2}. Collision probability <= (k-1)/q^2. Twice the multiplications.
 * lam is always a random element of F_{q^2}.
 *
 * 64-bit output: fold(z) = re ^ rotl64(im, 32). re and im are each in [0, q), i.e. 61 bits. The rotation
 * puts im's three always-zero top bits (61..63) at positions 29..31, where re's bits are live, and puts im's
 * live bits 29..31 at positions 61..63, where re's bits are always zero. So every output bit is the XOR of at
 * least one live bit, and the fold of a uniform (re, im) is uniform up to the 2^-61 bias of q < 2^61.
 */
#ifndef KP7_H
#define KP7_H
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define KP7_Q ((uint64_t)0x1FFFFFFFFFFFFFFFULL) /* 2^61 - 1 */

typedef struct { uint64_t re, im; } kp7_fe; /* element of F_{q^2}, both parts in [0, q) */

enum { KP7_XREAL = 0, KP7_XCOMPLEX = 1 };

typedef struct {
  unsigned k;
  int variant;
  kp7_fe x, xinv, lam;
  kp7_fe out_f[4], in_f[4]; /* -s*x^{-1},  s*x^{k-1} */
  kp7_fe out_r[4], in_r[4]; /*  s*x^{k},   -s        */
} kp7_params;

/* rolling state over one sequence */
typedef struct {
  const kp7_params* p;
  const char* seq;
  size_t len;
  size_t pos;      /* start of the current k-mer */
  int valid;       /* 1 when F, R hold the k-mer at pos */
  kp7_fe F, R;
} kp7_roller;

/* field */
uint64_t kp7_addq(uint64_t a, uint64_t b);
uint64_t kp7_subq(uint64_t a, uint64_t b);
uint64_t kp7_mulq(uint64_t a, uint64_t b);
kp7_fe kp7_add(kp7_fe a, kp7_fe b);
kp7_fe kp7_mul(kp7_fe a, kp7_fe b);
kp7_fe kp7_inv(kp7_fe a);
kp7_fe kp7_base(int code); /* 0..3 = A C G T */
int kp7_code(char c);      /* A C G T (either case) -> 0..3, else -1 */
uint64_t kp7_fold(kp7_fe z);

/* setup: seed chooses x and lam deterministically (splitmix64) */
void kp7_init(kp7_params* p, unsigned k, int variant, uint64_t seed);

/* direct (non-rolling) evaluation of one k-mer; returns 0 if it contains a non-ACGT base */
int kp7_direct(const kp7_params* p, const char* kmer, kp7_fe* F, kp7_fe* R);
kp7_fe kp7_sym(const kp7_params* p, kp7_fe F, kp7_fe R);

/* rolling: kp7_roll returns 1 and fills F/R/pos for each valid k-mer, skipping windows with non-ACGT */
void kp7_start(kp7_roller* r, const kp7_params* p, const char* seq, size_t len);
int kp7_roll(kp7_roller* r);

/* roll and produce the 64-bit output for one mode in the same call (what a caller would actually use) */
enum { KP7_OUT_FWD = 0, KP7_OUT_CANON_MIN = 1, KP7_OUT_SYM = 2 };
int kp7_roll64(kp7_roller* r, int mode, uint64_t* out);

#ifdef __cplusplus
}
#endif
#endif
