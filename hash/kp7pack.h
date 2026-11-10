/* Hybrid v2 (the KMC idea: 2-bit packed k-mers rolled by shifts; ideas only, no KMC code).
 *
 * For k <= 32 a k-mer packs exactly into 64 bits, 2 bits per base, with codes A0 C1 G2 T3, so complement is 3 - c.
 * Both strands roll with shifts and masks, no multiplies and no table of seeds:
 *   fwd' = ((fwd << 2) | c_in) & mask                        (first base in the high bits)
 *   rc'  = (rc >> 2) | ((3 - c_in) << 2(k-1))
 * The packed codes are an exact, injective encoding of each strand.
 *
 * Gate: gate(s) = mix(fwd) + mix(rc) mod 2^64. It is strand-symmetric because reverse complementing swaps fwd and rc.
 * It is NOT polynomial and its collisions are heuristic, as for ntHash's fwd + rev; test_kp7pack checks it
 * exhaustively at k = 7 and 8.
 *
 * Sampling: a k-mer with gate < threshold (FracMinHash rate threshold / 2^64) gets the PROVABLE KP-7 symmetric hash
 * fold((F + lam)(R + lam)) by DIRECT evaluation (kp7_direct + kp7_sym, about 2k field multiplies). That is
 * bit-identical to kp7_roll64(.., KP7_OUT_SYM, ..) with the same params. At rate 1/1000 the multiplies average about
 * 0.06 per k-mer, so the rolling field multiplies of hybrid v1 disappear.
 */
#ifndef KP7PACK_H
#define KP7PACK_H
#include "kp7.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the gate of one k-mer from its packed codes (reference) */
uint64_t kp7pack_gate(uint64_t fwd, uint64_t rc);
/* a seeded gate, mix(fwd ^ seed) + mix(rc ^ seed): still strand-symmetric; seed 0 is kp7pack_gate. For multi-seed
 * sampling statistics (stats_seeds.cpp). */
uint64_t kp7pack_gate_seeded(uint64_t fwd, uint64_t rc, uint64_t seed);
/* packed codes of one k-mer by direct construction; returns 0 if it has a non-ACGT base. k <= 32. */
int kp7pack_codes(const char* kmer, unsigned k, uint64_t* fwd, uint64_t* rc);

/* gate for every k-mer, in blocks of at most cap: emit(ctx, start, buf, n) per block of consecutive k-mers. k <= 32.
 * Returns the number of k-mers. */
size_t kp7pack_gate_seq(unsigned k, const char* seq, size_t len, uint64_t* buf, size_t cap, kp7_emit_fn emit,
                        void* ctx);

/* keep gate < threshold; emit (position, provable KP-7 symmetric hash) for kept k-mers in position order, at most cap
 * per call of emit. p->k must equal k (<= 32); either variant. Returns the number kept. */
size_t kp7pack_sample(const kp7_params* p, const char* seq, size_t len, uint64_t threshold, size_t* pos_buf,
                      uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx);

/* ---- multiply-free direct evaluation of the provable hash (next steps, item 1) ----
 * Codes A0 C1 G2 T3 make complement 3 - c, and val(3 - c) = -val(c) (A=1, C=i, G=-i, T=-1). The rc code holds the
 * complement of base t at digit t, and the fwd code holds base k-1-p at digit p, so with ONE shared table
 *   N[t][d] = -val(d) * x^t                                       (canonical F_{q^2} elements)
 *   F = sum_t N[t][digit t of rc],    R = -x^{k-1} p(1/x) = sum_p N[p][digit p of fwd].
 * Multiplying by +-1 and +-i is a sign or a re/im swap, so N needs only the powers x^t; F and R are additions.
 *   KP7PACK_EVAL_POS    k lookups and additions per strand (4-entry table per position)
 *   KP7PACK_EVAL_CHUNK  ceil(k/4) lookups per strand from byte tables C[j][b] = sum_u N[4j+u][(b >> 2u) & 3]
 *                       (8 x 256 x 16 B = 32 KB for k = 31, shared by both strands; the last chunk sums only
 *                       positions < k, and code bits above 2k are zero)
 *   KP7PACK_EVAL_DIRECT kp7_direct (Horner, about 2k field multiplies), the reference
 * Each sum is accumulated without reduction and reduced once, so F and R are canonical, and the result
 * fold(kp7_sym(F, R)) is bit-identical to kp7_roll64(.., KP7_OUT_SYM, ..). Either variant; k <= 32. */
enum { KP7PACK_EVAL_DIRECT = 0, KP7PACK_EVAL_POS = 1, KP7PACK_EVAL_CHUNK = 2 };
typedef struct {
  unsigned k, nchunk;
  kp7_fe N[32][4];
  kp7_fe C[8][256];
} kp7pack_tab;
void kp7pack_tab_init(kp7pack_tab* tab, const kp7_params* p);
/* F and R of one k-mer from its packed codes */
void kp7pack_eval(const kp7pack_tab* tab, int mode, uint64_t fwd, uint64_t rc, kp7_fe* F, kp7_fe* R);
/* kp7pack_sample with a chosen evaluation (tab may be NULL for KP7PACK_EVAL_DIRECT) */
size_t kp7pack_sample_eval(const kp7_params* p, const kp7pack_tab* tab, int mode, const char* seq, size_t len,
                           uint64_t threshold, size_t* pos_buf, uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit,
                           void* ctx);
/* the provable symmetric hash of EVERY k-mer from packed codes and chunk tables (no gate), in blocks of consecutive
 * k-mers as kp7pack_gate_seq; bit-identical to kp7_roll64(.., KP7_OUT_SYM, ..) */
size_t kp7pack_sym_all(const kp7_params* p, const kp7pack_tab* tab, const char* seq, size_t len, uint64_t* buf,
                       size_t cap, kp7_emit_fn emit, void* ctx);

/* ---- cheaper evaluation of kept k-mers (results/kp7_kept_eval_2026-11-10.md) ----
 * Same identity as above (F from rc digits, R from fwd digits, one shared table), with:
 *   - lam folded into chunk 0's table, so the sums ARE G = F + lam and H = R + lam
 *   - all chunks unrolled (tables for chunks >= k are all zero, so they add nothing)
 *   - one partial reduction of the sums (61-bit: < 2q; the wide product takes parts < 4q)
 *   - the product and fold inline, and kept k-mers collected without a branch, then evaluated in a batch
 * W is the chunk width in bits: 8 = 4-base chunks (8 x 256 entries), 4 = 2-base chunks (16 x 16 entries).
 * 61-bit tables (kp7pack_ftab) give fold(kp7_sym(F, R)), bit-identical to kp7_roll64(.., KP7_OUT_SYM, ..).
 * CM31 tables (kp7pack_mtab, 8-byte entries re | im << 32) give kp7m31_sym(F, R), the CM31 hash: bit-identical to
 * kp7m31_direct + kp7m31_sym (a DIFFERENT hash from the 61-bit one; a sketch must use one of them throughout). */
#include "kp7m31.h"
typedef struct {
  unsigned k;
  kp7_fe T8[8][256];  /* W = 8, lam in T8[0] */
  kp7_fe T4[16][16];  /* W = 4, lam in T4[0] */
} kp7pack_ftab;
typedef struct {
  unsigned k;
  uint64_t T8[8][256]; /* re | im << 32, lam in T8[0] */
  uint64_t T4[16][16];
} kp7pack_mtab;
void kp7pack_ftab_init(kp7pack_ftab* t, const kp7_params* p);
void kp7pack_mtab_init(kp7pack_mtab* t, const kp7m31_params* p);
/* every k-mer (no gate), blocks of consecutive k-mers as kp7pack_sym_all; W = 8 or 4 */
size_t kp7pack_sym_all_f(const kp7pack_ftab* t, int W, const char* seq, size_t len, uint64_t* buf, size_t cap,
                         kp7_emit_fn emit, void* ctx);
size_t kp7pack_sym_all_m(const kp7pack_mtab* t, int W, const char* seq, size_t len, uint64_t* buf, size_t cap,
                         kp7_emit_fn emit, void* ctx);
/* gate < threshold, batch-evaluated; same output order and emit contract as kp7pack_sample */
size_t kp7pack_sample_f(const kp7pack_ftab* t, int W, const char* seq, size_t len, uint64_t threshold, size_t* pos_buf,
                        uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx);
size_t kp7pack_sample_m(const kp7pack_mtab* t, int W, const char* seq, size_t len, uint64_t threshold, size_t* pos_buf,
                        uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx);
/* the same with a BRANCH on the gate and the fast evaluation inline (no per-k-mer stores): cheaper when few k-mers
 * are kept, since the branch is almost never taken; W = 8 only */
size_t kp7pack_sample_fb(const kp7pack_ftab* t, const char* seq, size_t len, uint64_t threshold, size_t* pos_buf,
                         uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx);
size_t kp7pack_sample_mb(const kp7pack_mtab* t, const char* seq, size_t len, uint64_t threshold, size_t* pos_buf,
                         uint64_t* hash_buf, size_t cap, kp7_emit_pairs_fn emit, void* ctx);

#ifdef __cplusplus
}
#endif
#endif
