"""KP-16 check: evaluating at {x, 1/x} carries the same information as KP-7's (F(x), R(x)).
Gaussian encoding over F_{q^2}, q = 2^61-1 (i^2 = -1). Exhaustive at small k; runs in seconds.
Usage: python3 literature/kp16_moebius_check.py  (literature/kp15_kp18_prior_art_2026-11-12.md)"""
import itertools, random
q = (1 << 61) - 1
def mul(a, b): return ((a[0]*b[0]-a[1]*b[1]) % q, (a[0]*b[1]+a[1]*b[0]) % q)
def add(a, b): return ((a[0]+b[0]) % q, (a[1]+b[1]) % q)
def neg(a): return ((-a[0]) % q, (-a[1]) % q)
def inv(a):
    n = (a[0]*a[0]+a[1]*a[1]) % q; ni = pow(n, q-2, q); return (a[0]*ni % q, (-a[1])*ni % q)
def pw(a, e):
    r = (1, 0)
    for _ in range(e): r = mul(r, a)
    return r
V = {'A': (1, 0), 'C': (0, 1), 'G': (0, q-1), 'T': (q-1, 0)}
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def ev(s, x):
    r = (0, 0)
    for c in reversed(s): r = add(mul(r, x), V[c])
    return r
def rc(s): return ''.join(COMP[c] for c in reversed(s))
random.seed(1)
y = (random.randrange(q), random.randrange(q)); x = mul(y, y); xi = inv(x)
for k in (4, 5, 6):
    bad = 0
    for t in itertools.product('ACGT', repeat=k):
        s = ''.join(t)
        F, R, Fi = ev(s, x), ev(rc(s), x), ev(s, xi)
        # identity: R(x) = -x^(k-1) F(1/x)
        if R != neg(mul(pw(x, k-1), Fi)): bad += 1
    # u - v and u*v vs F+R and F*R, u = y^-(k-1) F(x), v = y^(k-1) F(1/x)
    yk = pw(y, k-1); yki = inv(yk)
    cls_sum, cls_prod, cls_exact, classes = set(), set(), set(), set()
    lam = (random.randrange(q), random.randrange(q))
    for t in itertools.product('ACGT', repeat=k):
        s = ''.join(t); classes.add(min(s, rc(s)))
        u = mul(yki, ev(s, x)); v = mul(yk, ev(s, xi))
        cls_sum.add(add(u, neg(v))); cls_prod.add(mul(u, v))
        cls_exact.add(mul(add(u, lam), add(neg(v), lam)))
    print(f"k={k}: identity failures {bad}; classes {len(classes)}; distinct u-v {len(cls_sum)}; u*v {len(cls_prod)}; (u+l)(-v+l) {len(cls_exact)}")
