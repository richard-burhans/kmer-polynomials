/* KP-7 over the Eisenstein encoding: a strand-symmetric k-mer hash that is ONE element of F_q, q = 2^61 - 1.
 *
 * Encoding (results/ideas_algebra_2026-11-10.md, ideas/algebra/roots.py): A = 1, C = w, G = -w, T = -1, where w is a
 * primitive cube root of unity. q = 2^61 - 1 is 1 mod 3, so w lies in F_q itself (the Gaussian encoding needs i, and
 * q is 3 mod 4, which is why kp7.c works in F_{q^2}). Complement is negation, as before: A <-> T and C <-> G.
 * In the complex embedding (w = e^{2 pi i/3}) the squared distances are 1 for transitions (A-G, C-T), 3 for the A-C and
 * G-T transversions and 4 for complements. That metric lives in C only; F_q has no conjugation or absolute value, so
 * the hash does not "know" the metric. The field version keeps only the algebra (complement = negation, w^2+w+1 = 0).
 *
 * With x and lam independent and uniform in F_q (x != 0; see kp7e_init_seeds):
 *   F(s) = sum_t s_t x^t                 R(s) = rc(p_s)(x) = -x^{k-1} p_s(1/x)
 *   S(s) = (F + lam)(R + lam)            one canonical 61-bit value
 * Collision bound. The four base values are distinct in F_q, so two different k-mers give a nonzero difference
 * polynomial of degree <= k-1 in x: Pr[F(s) = F(u)] <= (k-1)/q. S, as a polynomial in lam over F_q(x), factors into
 * (lam + F)(lam + R), so it determines the unordered pair {p_s, rc(p_s)}, i.e. the double-stranded class; by
 * Schwartz-Zippel two different classes collide with probability <= deg/q over random (x, lam). The difference of
 * two S polynomials has total degree <= 2k - 2 (the lam^2 terms cancel), so the bound is 2(k-1)/q, about 2.6e-17 at
 * k = 31 (the looser 2k/q is 2.7e-17). That matches kp7's KP7_XREAL bound and is weaker than CM31's 2(k-1)/p^2
 * (1.3e-17) only by p^2/q ~ 2^1: it is NOT the F_{q^2} x-complex bound 2(k-1)/q^2.
 *
 * Rolling, both strands, with the lam state G = F + lam, H = R + lam and merged 16-entry tables [s_out*4 + s_in]:
 *   G' = G * x^{-1} + TF[.],   TF = -s_out x^{-1} + s_in x^{k-1} + lam (1 - x^{-1})
 *   H' = H * x      + TR[.],   TR =  s_out x^k   - s_in          + lam (1 - x)
 *   S  = G * H
 * three 61-bit modular multiplies per k-mer (kp7.c with x real needs 4 plus a full F_{q^2} product).
 */
#ifndef KP7E_H
#define KP7E_H
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define KP7E_Q ((uint64_t)0x1FFFFFFFFFFFFFFFULL)

typedef struct {
  unsigned k;
  uint64_t w;            /* primitive cube root of unity in F_q */
  uint64_t val[4];       /* A C G T = 1, w, -w, -1 (codes A0 C1 G2 T3: complement 3 - c) */
  uint64_t x, xinv, lam;
  uint64_t xk1, xk;      /* x^{k-1}, x^k */
  uint64_t tf[16], tr[16]; /* merged, lam-state tables [s_out*4 + s_in] */
} kp7e_params;

/* x and lam are drawn independently (hash/kp7seed.h): x uniform on F_q \ {0}, lam uniform on F_q.
 *   kp7e_init_seeds: seed_x / seed_lam NULL -> that parameter from OS randomness (getrandom), the default for use;
 *                    non-NULL -> a deterministic stream for it (tests, benchmarks).
 *   kp7e_init:       both deterministic. */
void kp7e_init_seeds(kp7e_params* p, unsigned k, const uint64_t* seed_x, const uint64_t* seed_lam);
void kp7e_init(kp7e_params* p, unsigned k, uint64_t seed_x, uint64_t seed_lam);
/* field helpers (canonical inputs and outputs) */
uint64_t kp7e_mul(uint64_t a, uint64_t b);
uint64_t kp7e_add(uint64_t a, uint64_t b);
uint64_t kp7e_inv(uint64_t a);
int kp7e_code(char c); /* A C G T either case -> 0..3, else -1 */

/* reference: F, R of one k-mer by Horner (canonical); 0 if it holds a non-ACGT base */
int kp7e_direct(const kp7e_params* p, const char* kmer, uint64_t* F, uint64_t* R);
/* S = (F + lam)(R + lam), canonical */
uint64_t kp7e_sym(const kp7e_params* p, uint64_t F, uint64_t R);

/* scalar reference roller: canonical F, R updated by the plain rolling identities (no tables, no laziness) */
typedef struct {
  const kp7e_params* p;
  const char* seq;
  size_t len, pos;
  int valid;
  uint64_t F, R;
} kp7e_roller;
void kp7e_start(kp7e_roller* r, const kp7e_params* p, const char* seq, size_t len);
int kp7e_roll(kp7e_roller* r); /* 1 and F/R/pos for each valid k-mer, skipping windows with N */

/* optimized whole-sequence symmetric hashing: lam state, merged tables, lazy reduction, `lanes` (1 or 2) interleaved
 * stretches per block. emit(ctx, start, hashes, n) per block of n consecutive k-mers starting at `start`; buf holds
 * cap hashes. Bit-identical to kp7e_sym over kp7e_roll (test_kp7e). Returns the number of k-mers. */
typedef void (*kp7e_emit_fn)(void* ctx, size_t start, const uint64_t* hashes, size_t n);
size_t kp7e_sym_seq(const kp7e_params* p, const char* seq, size_t len, int lanes, uint64_t* buf, size_t cap,
                    kp7e_emit_fn emit, void* ctx);

/* ---- hybrid v2 over the Eisenstein hash: packed-code gate, table evaluation of kept k-mers ----
 * The gate is kp7pack's (mix(fwd) + mix(rc), same constant), so a sampler keeps the same k-mers as kp7pack's samplers
 * (test_kp7e checks the gate against kp7pack_gate). With codes A0 C1 G2 T3 the rc code holds the complement of base t
 * at digit t and the fwd code holds base k-1-p at digit p, so one shared table N[t][d] = -val(d) x^t gives
 *   F = sum_t N[t][digit t of rc],   R = sum_p N[p][digit p of fwd]
 * (the same identity as kp7pack.h). Byte tables T[j][b] = sum_u N[4j+u][(b >> 2u) & 3] over positions < k, with lam
 * folded into T[0]: 8 x 256 x 8 B = 16 KB. G and H are 8 lookups and additions each, one partial reduction, then one
 * product. k <= 32. */
typedef struct {
  unsigned k;
  uint64_t T[8][256];
} kp7e_tab;
void kp7e_tab_init(kp7e_tab* t, const kp7e_params* p);
int kp7e_codes(const char* kmer, unsigned k, uint64_t* fwd, uint64_t* rc); /* packed codes, as kp7pack_codes */
uint64_t kp7e_tab_sym(const kp7e_tab* t, uint64_t fwd, uint64_t rc);      /* S from packed codes */
uint64_t kp7e_gate(uint64_t fwd, uint64_t rc);                            /* = kp7pack_gate */
typedef void (*kp7e_emit_pairs_fn)(void* ctx, const size_t* pos, const uint64_t* hashes, size_t n);
/* keep gate < threshold; emit (position, S) for kept k-mers in position order, at most cap per emit.
 *   _b: a branch on the gate, evaluation inline (best when almost nothing is kept)
 *   _batch: kept k-mers collected without a branch, then evaluated in a tight loop */
size_t kp7e_sample_b(const kp7e_tab* t, const char* seq, size_t len, uint64_t threshold, size_t* pos_buf,
                     uint64_t* hash_buf, size_t cap, kp7e_emit_pairs_fn emit, void* ctx);
size_t kp7e_sample_batch(const kp7e_tab* t, const char* seq, size_t len, uint64_t threshold, size_t* pos_buf,
                         uint64_t* hash_buf, size_t cap, kp7e_emit_pairs_fn emit, void* ctx);

#ifdef __cplusplus
}
#endif
#endif
