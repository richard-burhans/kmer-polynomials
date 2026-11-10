# KP-7: exhaustive check of strand-symmetric hash candidates over F_{q^2} = Z[i]/(q), q = 2^31-1 (q = 3 mod 4).
import itertools, random
q = 2**31 - 1
def mul(a, b): return ((a[0]*b[0] - a[1]*b[1]) % q, (a[0]*b[1] + a[1]*b[0]) % q)
def add(a, b): return ((a[0]+b[0]) % q, (a[1]+b[1]) % q)
def neg(a): return ((-a[0]) % q, (-a[1]) % q)
V = {'A': (1, 0), 'C': (0, 1), 'G': (0, q-1), 'T': (q-1, 0)}
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
random.seed(1)
x = (random.randrange(q), random.randrange(q)); lam = (random.randrange(q), random.randrange(q))
def F(s):  # p_s(x) = sum s_t x^t
    acc, p = (0, 0), (1, 0)
    for c in s: acc = add(acc, mul(V[c], p)); p = mul(p, x)
    return acc
def rc(s): return ''.join(COMP[c] for c in reversed(s))
for k in (7, 8):
    kmers = [''.join(t) for t in itertools.product('ACGT', repeat=k)]
    classes = {frozenset((s, rc(s))) for s in kmers}
    vals = {}
    for s in kmers:
        f, r = F(s), F(rc(s))
        vals[s] = (f, r)
    def distinct(h): return len({h(*vals[s]) for s in kmers})
    print(f"k={k}: {len(kmers)} k-mers, {len(classes)} double-stranded classes")
    print("  e1 = F+R            distinct values:", distinct(lambda f, r: add(f, r)))
    print("  e2 = F*R            distinct values:", distinct(lambda f, r: mul(f, r)))
    print("  (F+lam)(R+lam)      distinct values:", distinct(lambda f, r: mul(add(f, lam), add(r, lam))))
    # does e2 merge s with plain reverse(s)?
    s = 'AACGTTG' if k == 7 else 'AACGTTGC'  # not palindromes, so the reverse is a different k-mer
    print("  e2(s)==e2(reverse(s)):", mul(*vals[s]) == mul(*vals[s[::-1]]), " sym(s)==sym(reverse(s)):",
          mul(add(vals[s][0], lam), add(vals[s][1], lam)) == mul(add(vals[s[::-1]][0], lam), add(vals[s[::-1]][1], lam)))
