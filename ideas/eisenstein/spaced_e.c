/* KP-14 Eisenstein spaced-seed rolling; see spaced_e.h. */
#include "spaced_e.h"

static uint64_t xpow(const spe_params* p, int e) { return e < 0 ? p->xinv : f61_pow(p->x, (uint64_t)e); }
static int blocks(const char* mask, unsigned k, unsigned* a, unsigned* b)
{
  int n = 0;
  for (unsigned t = 0; t < k;) {
    if (mask[t] != '1') { t++; continue; }
    unsigned s = t;
    while (t < k && mask[t] == '1') t++;
    a[n] = s; b[n] = t; n++;
  }
  return n;
}

void spe_init(spe_params* p, const kp7e_params* kp, const char* mask)
{
  memset(p, 0, sizeof *p);
  p->k = kp->k; p->x = kp->x; p->xinv = kp->xinv; p->lam = kp->lam;
  memcpy(p->val, kp->val, sizeof p->val);
  memcpy(p->mask, mask, p->k);
  char rev[64];
  for (unsigned t = 0; t < p->k; t++) rev[t] = mask[p->k - 1 - t];
  p->nb = blocks(mask, p->k, p->a, p->b);
  p->nbr = blocks(rev, p->k, p->ar, p->br);
  const unsigned k = p->k;
  for (int i = 0; i < p->nb; i++)
    for (int c = 0; c < 4; c++) {
#ifdef SPE_MUTATE /* test-only: wrong outgoing power; `make mutation` expects FAIL */
      p->fo[i][c] = f61_neg(f61_mul(p->val[c], xpow(p, (int)p->a[i])));
#else
      p->fo[i][c] = f61_neg(f61_mul(p->val[c], xpow(p, (int)p->a[i] - 1)));
#endif
      p->fi[i][c] = f61_mul(p->val[c], xpow(p, (int)p->b[i] - 1));
    }
  for (int i = 0; i < p->nbr; i++)
    for (int c = 0; c < 4; c++) {
      p->ro[i][c] = f61_mul(p->val[c], xpow(p, (int)(k - p->ar[i])));
      p->ri[i][c] = f61_neg(f61_mul(p->val[c], xpow(p, (int)(k - p->br[i]))));
    }
}

int spe_direct(const spe_params* p, const char* w, uint64_t* F, uint64_t* R)
{
  uint64_t f = 0, r = 0;
  for (unsigned t = 0; t < p->k; t++) if (base_code(w[t]) < 0) return 0;
  for (unsigned t = 0; t < p->k; t++) {
    if (p->mask[t] == '1') f = f61_add(f, f61_mul(p->val[base_code(w[t])], xpow(p, (int)t)));
    if (p->mask[p->k - 1 - t] == '1') r = f61_add(r, f61_mul(p->val[base_code(w[t])], xpow(p, (int)(p->k - 1 - t))));
  }
  *F = f;
  *R = f61_neg(r);
  return 1;
}

int spe_restart(spe_roller* r, size_t from)
{
  const unsigned k = r->p->k;
  size_t run = 0, i = from;
  while (run < k) {
    if (i >= r->len) { r->valid = 0; r->pos = r->len; return 0; }
    run = base_code(r->seq[i]) < 0 ? 0 : run + 1;
    i++;
  }
  r->pos = i - k;
  spe_direct(r->p, r->seq + r->pos, &r->F, &r->R);
  r->valid = 1;
  return 1;
}
