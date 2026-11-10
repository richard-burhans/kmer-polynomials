/* KP-14 spaced-seed rolling; see spaced.h. */
#include "spaced.h"

static uint64_t xpow(const sp_params* p, int e) /* x^e for e >= -1 */
{
  return e < 0 ? p->xinv : f61_pow(p->x, (uint64_t)e);
}

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

void sp_init(sp_params* p, const kp7_params* kp, const char* mask)
{
  memset(p, 0, sizeof *p);
  p->k = kp->k;
  p->x = kp->x.re;
  p->xinv = kp->xinv.re;
  p->lam = kp->lam;
  memcpy(p->mask, mask, p->k);
  char rev[64];
  for (unsigned t = 0; t < p->k; t++) rev[t] = mask[p->k - 1 - t];
  p->nb = blocks(mask, p->k, p->a, p->b);
  p->nbr = blocks(rev, p->k, p->ar, p->br);
  const unsigned k = p->k;
  for (int i = 0; i < p->nb; i++)
    for (int c = 0; c < 4; c++) {
#ifdef SP_MUTATE /* test-only: wrong outgoing power; `make mutation` expects FAIL */
      p->fo[i][c] = g_neg(g_mul_real(g_base(c), xpow(p, (int)p->a[i])));
#else
      p->fo[i][c] = g_neg(g_mul_real(g_base(c), xpow(p, (int)p->a[i] - 1)));
#endif
      p->fi[i][c] = g_mul_real(g_base(c), xpow(p, (int)p->b[i] - 1));
    }
  for (int i = 0; i < p->nbr; i++)
    for (int c = 0; c < 4; c++) {
      p->ro[i][c] = g_mul_real(g_base(c), xpow(p, (int)(k - p->ar[i])));
      p->ri[i][c] = g_neg(g_mul_real(g_base(c), xpow(p, (int)(k - p->br[i]))));
    }
}

int sp_direct(const sp_params* p, const char* w, kp7_fe* F, kp7_fe* R)
{
  kp7_fe f = { 0, 0 }, r = { 0, 0 };
  for (unsigned t = 0; t < p->k; t++) if (base_code(w[t]) < 0) return 0;
  for (unsigned t = 0; t < p->k; t++) {
    if (p->mask[t] == '1') f = g_add(f, g_mul_real(g_base(base_code(w[t])), xpow(p, (int)t)));
    unsigned u = t; /* Mrev membership: u in Mrev iff mask[k-1-u] == '1' */
    if (p->mask[p->k - 1 - u] == '1') r = g_add(r, g_mul_real(g_base(base_code(w[u])), xpow(p, (int)(p->k - 1 - u))));
  }
  *F = f;
  *R = g_neg(r);
  return 1;
}

int sp_restart(sp_roller* r, size_t from)
{
  const unsigned k = r->p->k;
  size_t run = 0, i = from;
  while (run < k) {
    if (i >= r->len) { r->valid = 0; r->pos = r->len; return 0; }
    run = base_code(r->seq[i]) < 0 ? 0 : run + 1;
    i++;
  }
  r->pos = i - k;
  sp_direct(r->p, r->seq + r->pos, &r->F, &r->R);
  r->valid = 1;
  return 1;
}
