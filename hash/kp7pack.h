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

#ifdef __cplusplus
}
#endif
#endif
