/* Cheaper evaluation of KEPT k-mers in the hybrid-v2 sampler (results/kp7_kept_eval_2026-11-11.md).
 *
 * The sampler is unchanged: packed 2-bit codes (A0 C1 G2 T3, complement 3 - c) rolled by shifts, the strand-symmetric
 * gate mix(fwd) + mix(rc) (kp7pack's, so every sampler here keeps exactly the k-mers kp7pack and kp7e keep), and
 * the provable hash evaluated only for k-mers with gate < threshold. What changes is how a kept k-mer is evaluated.
 *
 * The identity everything rests on. With complement = 3 - c and val(3 - c) = -val(c) (true for every encoding in this
 * repository: Gaussian 1, i, -i, -1 and Eisenstein 1, w, -w, -1), the per-digit value is AFFINE in the digit's two
 * bits b0 (low) and b1 (high):
 *     -val(d) = -val(0) + (val(0) - val(1)) b0 + (val(0) - val(2)) b1,
 * because the b0*b1 coefficient is val(3) - val(2) - val(1) + val(0) = 0. So with F = sum_t N[t][digit t of rc] and
 * R = sum_p N[p][digit p of fwd] (kp7pack.h), N[t][d] = -val(d) x^t:
 *     F = -val(0) X + alpha P0(rc) + beta P1(rc),    X = sum_{t<k} x^t,  alpha = 1 - val(1),  beta = 1 + val(1)
 * where P0 (P1) is the subset sum of x^t over the digits whose low (high) bit is set. Two consequences:
 *   1. Chunk tables need not align to bases: any B-bit slice of the 2k-bit code has a table (CHUNK, any B).
 *   2. One table of subset sums of x^t serves both bit planes and both strands (PLANES): 4 plane sums per kept k-mer.
 *      alpha, beta = 1 -+ i are free in Z[i] ((1 - i)(a + bi) = (a + b) + (b - a)i), so Gaussian planes need only
 *      additions; for Eisenstein alpha P0 + beta P1 = P0 + P1 + w (P1 - P0), one multiply by w per strand.
 *
 * Evaluators (all bit-identical to their family's reference; test_kp7kept checks every one):
 *   E_CHUNK   Eisenstein (kp7e), B-bit chunks of the code, ceil(64/B) lookups per strand, lam folded into chunk 0;
 *             B = 8 is kp7e's table. Output S = (F + lam)(R + lam) in F_q, as kp7e_tab_sym.
 *   E_PLANES  Eisenstein, one subset-sum table, 4 * ceil(32/B) lookups, + 2 multiplies by w.
 *   G_PLANES  Gaussian 61-bit, x REAL (KP7_XREAL): real 8-byte subset sums; output fold(kp7_sym(F, R)), as
 *             kp7pack_sample_f. (x real is what makes the table real; with x complex the entries would be 16 bytes.)
 *   M_PLANES  CM31 (x complex): complex subset sums packed re | im << 32; output kp7m31_sym, as kp7pack_sample_m.
 *   M_CHUNK   CM31 with B-bit chunks of the code (B = 8 is kp7pack_mtab's 16 KB table).
 *   NT_CHUNK  BASELINE, not KP-7: ntHash's canonical hash (fwd + rev, BirolLab/ntHash c26bd45) of the kept k-mer, from
 *             tables of rotated seeds (srol is a bit permutation, so XOR-linear). B = 8 or 16 (base-aligned only: the
 *             four seeds are not affine in the code bits). Bit-identical to ntbatch.h / NtHash::roll().
 *   UH        BASELINE, not KP-7: canonical code min(fwd, rc) through multiply-shift with a random odd 128-bit a,
 *             h = (a * c mod 2^128) >> 64. Universal: Pr[collision] <= 2^-63 per pair of distinct canonical codes
 *             (Dietzfelbinger et al. 1997), k <= 32 only, no algebra.
 *
 * Collection: BATCHED (no branch on the gate; fwd, rc, pos stored every k-mer, the count advanced by the comparison,
 * then the batch evaluated in a tight loop; batch size 1..KP7KEPT_BATCH_MAX) or BRANCHY (evaluated inline when the
 * gate passes). Output order and the emit contract are kp7pack_sample's. k <= 32.
 */
#ifndef KP7KEPT_H
#define KP7KEPT_H
#include <stddef.h>
#include <stdint.h>
#include "kp7.h"
#include "kp7e.h"
#include "kp7m31.h"

#ifdef __cplusplus
extern "C" {
#endif

enum { KP7KEPT_E_CHUNK = 0, KP7KEPT_E_PLANES, KP7KEPT_G_PLANES, KP7KEPT_M_PLANES, KP7KEPT_M_CHUNK, KP7KEPT_NT_CHUNK,
       KP7KEPT_UH, KP7KEPT_NKIND };
#define KP7KEPT_BATCH_MAX 4096

typedef struct {
  int kind;
  unsigned k, B, nch; /* chunk width in bits and chunks per lookup set */
  uint64_t* T;        /* nch tables of 2^B entries, T[(j << B) | v]; NULL for UH */
  uint64_t c;         /* E_PLANES: lam - X;   M_PLANES: (lam - X) as re | im << 32 */
  uint64_t cre, cim;  /* G_PLANES: lam - X */
  uint64_t w;         /* E_PLANES: the cube root of unity */
  uint64_t a_lo, a_hi; /* UH multiplier (odd) */
  size_t batch;       /* BATCHED collection size, default 512 */
} kp7kept;

/* constructors; B in {8, 10, 11, 12, 13, 16} for CHUNK and {8, 11, 16} for PLANES; return 0 on a bad argument */
int kp7kept_init_e_chunk(kp7kept* t, const kp7e_params* p, unsigned B);
int kp7kept_init_e_planes(kp7kept* t, const kp7e_params* p, unsigned B);
int kp7kept_init_g_planes(kp7kept* t, const kp7_params* p, unsigned B); /* p->variant must be KP7_XREAL */
int kp7kept_init_m_planes(kp7kept* t, const kp7m31_params* p, unsigned B);
int kp7kept_init_m_chunk(kp7kept* t, const kp7m31_params* p, unsigned B);
int kp7kept_init_nt(kp7kept* t, unsigned k, unsigned B);
int kp7kept_init_uh(kp7kept* t, unsigned k, uint64_t seed);
void kp7kept_free(kp7kept* t);

/* the hash of one k-mer from its packed codes (as kp7pack_codes / kp7e_codes build them) */
uint64_t kp7kept_eval(const kp7kept* t, uint64_t fwd, uint64_t rc);
/* the gate (= kp7pack_gate) */
uint64_t kp7kept_gate(uint64_t fwd, uint64_t rc);

/* keep gate < threshold; emit (position, hash) in position order, at most cap per emit. Returns the number kept. */
size_t kp7kept_sample(const kp7kept* t, int branchy, const char* seq, size_t len, uint64_t threshold,
                      size_t* pos_buf, uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx);
/* UH of EVERY k-mer (no gate): the round-1 "canonical code + universal hash" row. Blocks of consecutive k-mers as
 * kp7pack_gate_seq. Returns the number of k-mers. */
size_t kp7kept_uh_all(const kp7kept* t, const char* seq, size_t len, uint64_t* buf, size_t cap, kp7_emit_fn emit,
                      void* ctx);

#ifdef __cplusplus
}
#endif
#endif
