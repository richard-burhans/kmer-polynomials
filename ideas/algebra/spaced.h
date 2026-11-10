/* KP-14: spaced seeds as coefficient masks, rolled in O(1 multiply + table adds) per k-mer.
 *
 * Mask m (length k, '1' = care) with care set M. For the k-mer starting at j,
 *   F_m(j) = sum_{t in M} s_{j+t} x^t                         (forward, masked)
 *   R_m(j) = F_m(rc(window j)) = -sum_{u in Mrev} s_{j+u} x^{k-1-u},  Mrev = {k-1-t : t in M}
 *   S_m(j) = fold((F_m + lam)(R_m + lam))                      strand-symmetric for ANY mask:
 *            rc swaps F_m and R_m, because R_m is defined as F_m of the reverse complement.
 * With the all-ones mask these are exactly kp7's F, R and symmetric hash (x real).
 *
 * ROLLING. Split M into maximal blocks [a,b). Because F_m is linear in the bases,
 *   F_m(j+1) = x^-1 F_m(j) - sum_blocks s_{j+a} x^{a-1} + sum_blocks s_{j+b} x^{b-1}
 *   R_m(j+1) = x R_m(j)    + sum_rblocks s_{j+a'} x^{k-a'} - sum_rblocks s_{j+b'} x^{k-b'}
 * (blocks [a',b') of Mrev). So every block shares ONE multiply per strand; each block adds two table entries,
 * and multiplying a base value by a precomputed power is a table lookup (no multiply). ntHash2's SeedNtHash rolls
 * blocks too (with rotations); the polynomial form amortises the multiply across all blocks of a mask.
 *
 * A window containing N in ANY of its k positions (care or not) is skipped, as SeedNtHash does. */
#ifndef IDEAS_SPACED_H
#define IDEAS_SPACED_H
#include "field61.h"
#include <stddef.h>
#include <string.h>

#define SP_MAXB 32
typedef struct {
  unsigned k;
  uint64_t x, xinv;      /* x real (kp7 KP7_XREAL params) */
  kp7_fe lam;
  char mask[64];
  int nb, nbr;           /* blocks of M and of Mrev */
  unsigned a[SP_MAXB], b[SP_MAXB], ar[SP_MAXB], br[SP_MAXB];
  kp7_fe fo[SP_MAXB][4], fi[SP_MAXB][4], ro[SP_MAXB][4], ri[SP_MAXB][4];
} sp_params;

typedef struct {
  const sp_params* p;
  const char* seq;
  size_t len, pos;
  int valid;
  kp7_fe F, R;
} sp_roller;

#ifdef __cplusplus
extern "C" {
#endif
void sp_init(sp_params* p, const kp7_params* kp, const char* mask); /* kp must be KP7_XREAL */
/* direct evaluation of window j; returns 0 if it holds an N */
int sp_direct(const sp_params* p, const char* w, kp7_fe* F, kp7_fe* R);
static inline uint64_t sp_sym(const sp_params* p, kp7_fe F, kp7_fe R) { return g_fold(g_mul(g_add(F, p->lam), g_add(R, p->lam))); }

static inline void sp_start(sp_roller* r, const sp_params* p, const char* seq, size_t len)
{
  r->p = p; r->seq = seq; r->len = len; r->pos = 0; r->valid = 0;
}

int sp_restart(sp_roller* r, size_t from); /* in spaced.c */
#ifdef __cplusplus
}
#endif

/* next window: returns 1 with r->pos and *out (the symmetric hash) set, 0 at the end */
static inline int sp_roll(sp_roller* r, uint64_t* out)
{
  const sp_params* p = r->p;
  if (!r->valid) {
    if (!sp_restart(r, r->pos)) return 0;
  } else {
    size_t in = r->pos + p->k;
    if (in >= r->len) { r->valid = 0; r->pos = r->len; return 0; }
    if (base_code(r->seq[in]) < 0) { if (!sp_restart(r, in + 1)) return 0; }
    else {
      const char* s = r->seq + r->pos;
      kp7_fe F = g_mul_real(r->F, p->xinv), R = g_mul_real(r->R, p->x);
      for (int i = 0; i < p->nb; i++) {
        F = g_add(F, p->fo[i][base_code(s[p->a[i]])]);
        F = g_add(F, p->fi[i][base_code(s[p->b[i]])]);
      }
      for (int i = 0; i < p->nbr; i++) {
        R = g_add(R, p->ro[i][base_code(s[p->ar[i]])]);
        R = g_add(R, p->ri[i][base_code(s[p->br[i]])]);
      }
      r->F = F; r->R = R;
      r->pos++;
    }
  }
  *out = sp_sym(p, r->F, r->R);
  return 1;
}
#endif
