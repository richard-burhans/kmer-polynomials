/* KP-7 over CM31 = F_{p^2} = Z[i]/(p), p = 2^31 - 1 (p = 3 mod 4, so i is native): the strand-symmetric
 * polynomial hash with 32-bit lanes, so it vectorizes (AVX2: 8 lanes, AVX-512: 16 lanes).
 *
 * A DIFFERENT hash from kp7.h's (q = 2^61 - 1); same construction:
 *   encoding A = 1, C = i, G = -i, T = -1 (complement = negation)
 *   F = sum_t s_t x^t,  R = rc(p)(x) = -x^{k-1} p(1/x),  sym = (F + lam)(R + lam)
 * with x and lam random in F_{p^2} (x COMPLEX). Two distinct double-stranded k-mers collide with probability
 * <= about 2k / p^2 (Schwartz-Zippel, notes/kp7_strand_symmetric_hash.md); for k = 31, ~1.4e-17 per pair.
 *
 * Output: pack(z) = re | (im << 31), a 62-bit value (injective on F_{p^2}; the top two bits are zero).
 *
 * Base codes: (c >> 1) & 3 for A C T G in either case: A -> 0, C -> 1, T -> 2, G -> 3. That is vectorizable;
 * tables are built for this mapping. N (any other byte) is found by a scalar pass that splits N-free runs.
 *
 * Rolling with the lam state G = F + lam, H = R + lam and merged tables (as kp7.c):
 *   G' = G * x^{-1} + tf[co*4 + ci],  tf = -v(co) x^{-1} + v(ci) x^{k-1} + lam (1 - x^{-1})
 *   H' = H * x      + tr[co*4 + ci],  tr =  v(co) x^{k}  - v(ci)         + lam (1 - x)
 *   sym = G * H
 */
#ifndef KP7M31_H
#define KP7M31_H
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define KP7M31_P 0x7FFFFFFFu

typedef struct { uint32_t re, im; } m31c; /* element of F_{p^2}, both parts in [0, p) */

typedef struct {
  unsigned k;
  m31c x, xinv, lam, xk1, xk;
  m31c val[4];         /* by base code: A, C, T, G */
  m31c tf[16], tr[16]; /* merged lam-state tables, index code(out)*4 + code(in) */
} kp7m31_params;

void kp7m31_init(kp7m31_params* p, unsigned k, uint64_t seed);
int kp7m31_code(char c); /* A C T G either case -> 0..3, else -1 */

/* field (exported for tests) */
m31c kp7m31_mul(m31c a, m31c b);
m31c kp7m31_add(m31c a, m31c b);
uint64_t kp7m31_pack(m31c z);

/* reference: direct evaluation of one k-mer; returns 0 if it contains a non-ACGT byte */
int kp7m31_direct(const kp7m31_params* p, const char* kmer, m31c* F, m31c* R);
uint64_t kp7m31_sym(const kp7m31_params* p, m31c F, m31c R); /* pack((F + lam)(R + lam)) */

/* reference roller (canonical F, R; no lam state): kp7m31_roll returns 1 and sets pos/F/R for each valid k-mer */
typedef struct { const kp7m31_params* p; const char* seq; size_t len, pos; int valid; m31c F, R; } kp7m31_roller;
void kp7m31_start(kp7m31_roller* r, const kp7m31_params* p, const char* seq, size_t len);
int kp7m31_roll(kp7m31_roller* r);

/* Whole-sequence symmetric hashing, blocks of at most `cap` k-mers per N-free run (cap is clamped to 65536):
 * emit(ctx, start, buf, n) gets the hashes of n consecutive k-mers starting at sequence position `start`.
 *   impl 0: scalar (lam state)    impl 1: AVX2, 8 lanes    impl 2: AVX-512, 16 lanes
 * All three are bit-identical (test_kp7m31). Returns the number of k-mers hashed, or 0 if `impl` is not
 * available on this CPU / build (see kp7m31_have). Blocks shorter than the lane count run scalar. */
enum { KP7M31_SCALAR = 0, KP7M31_AVX2 = 1, KP7M31_AVX512 = 2 };
int kp7m31_have(int impl);
typedef void (*kp7m31_emit_fn)(void* ctx, size_t start, const uint64_t* hashes, size_t n);
size_t kp7m31_sym_seq(const kp7m31_params* p, const char* seq, size_t len, int impl, uint64_t* buf, size_t cap,
                      kp7m31_emit_fn emit, void* ctx);

#ifdef __cplusplus
}
#endif
#endif
