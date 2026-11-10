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
 * Karatsuba (params.karatsuba = 1): (a+bi)(c+di) has imaginary part (a+b)(c+d) - ac - bd, so 3 widening
 * multiplies instead of 4. Exact in 128 bits: a+b, c+d < 2^62, the product is < 2^124 (red() accepts it), and
 * (a+b)(c+d) >= ac + bd, so the subtraction never wraps. It gives bit-identical values (test_kp7 checks this).
 * It applies to the symmetric product and, for KP7_XCOMPLEX, to the rolling multiplies.
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
  int karatsuba;   /* 1: F_{q^2} products use 3 widening multiplies (Karatsuba) instead of 4; same values */
  int opt;         /* KP7_OPT_* bits for kp7_roll_sym; every combination gives the same hashes */
  kp7_fe x, xinv, lam;
  kp7_fe out_f[4], in_f[4]; /* -s*x^{-1},  s*x^{k-1} */
  kp7_fe out_r[4], in_r[4]; /*  s*x^{k},   -s        */
  kp7_fe tf[16], tr[16];     /* merged [co*4+ci]: out_f+in_f, out_r+in_r                    (KP7_OPT_MERGED) */
  kp7_fe tfl[16], trl[16];   /* the same plus lam(1-x^{-1}) and lam(1-x), for the G/H state (both opts) */
  kp7_fe dl_f, dl_r;         /* lam(1-x^{-1}), lam(1-x) alone                     (KP7_OPT_LAMSTATE only) */
} kp7_params;

/* rolling state over one sequence */
typedef struct {
  const kp7_params* p;
  const char* seq;
  size_t len;
  size_t pos;      /* start of the current k-mer */
  int valid;       /* 1 when F, R hold the k-mer at pos */
  kp7_fe F, R;
  kp7_fe G, H;     /* kp7_roll_sym's state: F, R, or F+lam, R+lam with KP7_OPT_LAMSTATE */
  size_t resume;   /* kp7_sym_batch: where the next restart searches (pos keeps the last k-mer written) */
} kp7_roller;

/* field */
uint64_t kp7_addq(uint64_t a, uint64_t b);
uint64_t kp7_subq(uint64_t a, uint64_t b);
uint64_t kp7_mulq(uint64_t a, uint64_t b);
kp7_fe kp7_add(kp7_fe a, kp7_fe b);
kp7_fe kp7_mul(kp7_fe a, kp7_fe b);
kp7_fe kp7_mul_karatsuba(kp7_fe a, kp7_fe b);
kp7_fe kp7_mul_wide(kp7_fe a, kp7_fe b); /* parts may be in [0, 4q); result canonical */
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

/* Symmetric hash only, with optional addition-side optimizations (params.opt bits). Bit-identical to
 * kp7_roll64(.., KP7_OUT_SYM, ..) for every combination (test_kp7):
 *   KP7_OPT_LAMSTATE  roll G = F + lam and H = R + lam directly. The update constant absorbs lam(1 - x^{-1}) and
 *                     lam(1 - x), so the two "+ lam" additions per k-mer disappear.
 *   KP7_OPT_MERGED    one 16-entry table per strand [s_out*4 + s_in] instead of two lookups and two additions.
 *   KP7_OPT_LAZY      (x real only; ignored for x complex) the rolling additions skip their conditional
 *                     subtraction, so the state lives in [0, 2q). The x-real multiply takes < 2q inputs
 *                     ((2q)*q < 2^123), and the final product takes parts < 4q with a 16q^2 offset and a reduction
 *                     (red2) valid for any 128-bit value, so the output is still canonical. */
enum { KP7_OPT_LAMSTATE = 1, KP7_OPT_MERGED = 2, KP7_OPT_LAZY = 4 };
int kp7_roll_sym(kp7_roller* r, uint64_t* out);

/* Batch symmetric hashing: writes up to `max` hashes of CONSECUTIVE k-mers to out[] and returns how many. It stops
 * early at a window containing N (the next call restarts after it) or at the end of the sequence (then returns 0).
 * The last k-mer written starts at r->pos, so the first starts at r->pos - (n - 1). Bit-identical to
 * kp7_roll64(.., KP7_OUT_SYM, ..) (test_kp7). It always uses the lam state and merged tables, plus lazy reduction for
 * KP7_XREAL, and keeps the rolling state in locals: no per-k-mer call, no per-k-mer state store. It ignores
 * params.opt; params.karatsuba still applies to KP7_XCOMPLEX. Do not mix it with kp7_roll/kp7_roll_sym on one roller. */
size_t kp7_sym_batch(kp7_roller* r, uint64_t* out, size_t max);

/* Interleaved symmetric hashing of a whole sequence (step 2): each maximal N-free run is cut into blocks of at most
 * `cap` k-mers, and each block into `lanes` contiguous stretches rolled in the SAME loop iteration, so their serial
 * multiply-reduce-add chains overlap. A lane starts from a direct evaluation (k multiplies per lane per block).
 * For every block it calls emit(ctx, start, buf, n) with the hashes of the n consecutive k-mers starting at
 * sequence position `start`. buf must hold `cap` hashes. lanes: 1..4 for KP7_XREAL (lazy reduction), 1..2 for
 * KP7_XCOMPLEX; other values are clamped. Bit-identical to kp7_roll64(.., KP7_OUT_SYM, ..) (test_kp7). Returns the
 * number of k-mers hashed. */
typedef void (*kp7_emit_fn)(void* ctx, size_t start, const uint64_t* hashes, size_t n);
size_t kp7_sym_seq(const kp7_params* p, const char* seq, size_t len, int lanes, uint64_t* buf, size_t cap,
                   kp7_emit_fn emit, void* ctx);

#ifdef __cplusplus
}
#endif
#endif
