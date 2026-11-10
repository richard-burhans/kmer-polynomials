/* KP-14 Eisenstein revisit: spaced seeds over the Eisenstein KP-7 hash (hash/kp7e.h), ONE F_q element per strand.
 *
 * Same construction as ideas/algebra/spaced.h, with A = 1, C = w, G = -w, T = -1 in F_q (q = 2^61 - 1):
 *   F_m(j) = sum_{t in M} s_{j+t} x^t,   R_m(j) = F_m(rc(window j)),   S_m(j) = (F_m + lam)(R_m + lam)
 * strand-symmetric for ANY mask (rc swaps F_m and R_m). With the all-ones mask these equal kp7e's F, R and
 * kp7e_sym. Rolling: per mask ONE 61-bit multiply per strand (the Gaussian version multiplied a pair by a real x:
 * two multiplies per strand) plus two table additions per block, then one multiply for S.
 * A window with N anywhere in its k positions is skipped (as ideas/algebra/spaced.h does). */
#ifndef IDEAS_SPACED_E_H
#define IDEAS_SPACED_E_H
#include "../../hash/kp7e.h"
#include "../algebra/field61.h" /* f61_* helpers (same field) and base_code (A0 C1 G2 T3, as kp7e) */
#include <stddef.h>
#include <string.h>

#define SPE_MAXB 32
typedef struct {
  unsigned k;
  uint64_t x, xinv, lam, val[4];
  char mask[64];
  int nb, nbr;
  unsigned a[SPE_MAXB], b[SPE_MAXB], ar[SPE_MAXB], br[SPE_MAXB];
  uint64_t fo[SPE_MAXB][4], fi[SPE_MAXB][4], ro[SPE_MAXB][4], ri[SPE_MAXB][4];
} spe_params;

typedef struct {
  const spe_params* p;
  const char* seq;
  size_t len, pos;
  int valid;
  uint64_t F, R;
} spe_roller;

#ifdef __cplusplus
extern "C" {
#endif
void spe_init(spe_params* p, const kp7e_params* kp, const char* mask);
int spe_direct(const spe_params* p, const char* w, uint64_t* F, uint64_t* R); /* 0 if the window holds N */
int spe_restart(spe_roller* r, size_t from);
#ifdef __cplusplus
}
#endif
static inline uint64_t spe_sym(const spe_params* p, uint64_t F, uint64_t R) { return f61_mul(f61_add(F, p->lam), f61_add(R, p->lam)); }
static inline void spe_start(spe_roller* r, const spe_params* p, const char* seq, size_t len)
{
  r->p = p; r->seq = seq; r->len = len; r->pos = 0; r->valid = 0;
}
/* next window: 1 with r->pos and *out (the symmetric hash) set, 0 at the end */
static inline int spe_roll(spe_roller* r, uint64_t* out)
{
  const spe_params* p = r->p;
  if (!r->valid) {
    if (!spe_restart(r, r->pos)) return 0;
  } else {
    size_t in = r->pos + p->k;
    if (in >= r->len) { r->valid = 0; r->pos = r->len; return 0; }
    if (base_code(r->seq[in]) < 0) { if (!spe_restart(r, in + 1)) return 0; }
    else {
      const char* s = r->seq + r->pos;
      uint64_t F = f61_mul(r->F, p->xinv), R = f61_mul(r->R, p->x);
      for (int i = 0; i < p->nb; i++) {
        F = f61_add(F, p->fo[i][base_code(s[p->a[i]])]);
        F = f61_add(F, p->fi[i][base_code(s[p->b[i]])]);
      }
      for (int i = 0; i < p->nbr; i++) {
        R = f61_add(R, p->ro[i][base_code(s[p->ar[i]])]);
        R = f61_add(R, p->ri[i][base_code(s[p->br[i]])]);
      }
      r->F = F; r->R = R;
      r->pos++;
    }
  }
  *out = spe_sym(p, r->F, r->R);
  return 1;
}
#endif
