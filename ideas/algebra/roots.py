"""Root assignment (IDEAS.md section 7): can a base map with complement = negation separate transitions from
transversions?  Run: python3 roots.py   (exits nonzero if any check fails)

Write A, C in C (or a field), with T = -A and G = -C. Then
  transitions   A<->G : |A - G| = |A + C|,    C<->T : |C - T| = |C + A|        (always equal)
  transversions A<->C, G<->T : |A - C|;   A<->T : 2|A|;   C<->G : 2|C|
|A+C|^2 - |A-C|^2 = 4 Re(A conj C), so transitions differ from the A-C/G-T transversions iff Re(A conj C) != 0,
and are SHORTER iff the angle between A and C exceeds 90 degrees. Fourth roots (A=1, C=i) are exactly the
orthogonal case, so they cannot separate them. Sixth roots (A=1, C=w=e^{2 pi i/3}, G=-w, T=-1) give squared
distances 1 (transitions), 3 (A-C, G-T), 4 (complements): transitions are the cheapest substitution.

In finite fields: w (a primitive cube root of 1) lies in F_q itself when q = 1 (mod 3). Both Mersenne primes used
here qualify (2^61-1 and 2^31-1 are 1 mod 3), while i needs q = 1 (mod 4) and both are 3 (mod 4), which is why the
Gaussian encoding needs F_{q^2}. So the Eisenstein encoding puts every hash value in F_q: one field element, not two.
The Euclidean distances above live in the complex embedding only; F_q has no conjugation.
"""
import cmath, itertools, random, sys

fails = []
def check(cond, msg):
    if not cond: fails.append(msg)

def dists(A, C):
    T, G = -A, -C
    return {"A-G": abs(A - G), "C-T": abs(C - T), "A-C": abs(A - C), "G-T": abs(G - T), "A-T": abs(A - T), "C-G": abs(C - G)}

print("distances (squared) under complement = negation")
for name, A, C in (("4th roots (A=1, C=i)", 1, 1j), ("6th roots (A=1, C=w)", 1, cmath.exp(2j * cmath.pi / 3))):
    d = {k: round(v * v, 6) for k, v in dists(A, C).items()}
    print(f"  {name}: transitions A-G {d['A-G']}, C-T {d['C-T']}; transversions A-C {d['A-C']}, G-T {d['G-T']}, A-T {d['A-T']}, C-G {d['C-G']}")
d4 = dists(1, 1j); d6 = dists(1, cmath.exp(2j * cmath.pi / 3))
check(abs(d4["A-G"] - d4["A-C"]) < 1e-12, "4th roots should tie transitions with A-C transversions")
check(d6["A-G"] < min(d6["A-C"], d6["A-T"], d6["C-G"]) - 1e-9, "6th roots should make transitions the shortest")
# sweep the angle: separation iff Re(A conj C) != 0
for deg in range(1, 180):
    C = cmath.exp(1j * cmath.pi * deg / 180)
    d = dists(1, C)
    re = (1 * C.conjugate()).real
    sep = abs(d["A-G"] - d["A-C"]) > 1e-9
    check(sep == (abs(re) > 1e-12), f"angle {deg}: separation should match Re(A conj C) != 0")

# finite fields
def has_cube_root(q):
    return (q - 1) % 3 == 0
def cube_root_of_unity(q):
    random.seed(1)
    while True:
        g = random.randrange(2, q - 1)
        w = pow(g, (q - 1) // 3, q)
        if w != 1:
            return w
for q in (2**61 - 1, 2**31 - 1):
    print(f"q = 2^{q.bit_length()}-1: q mod 3 = {q % 3}, q mod 4 = {q % 4}")
    check(has_cube_root(q), f"q={q} should have a cube root of unity in F_q")
    w = cube_root_of_unity(q)
    check((w * w + w + 1) % q == 0, "w^2 + w + 1 = 0")
    vals = {1, w, (-w) % q, q - 1}
    check(len(vals) == 4, "the four base values must be distinct in F_q")
    print(f"  w = {w}; w^2+w+1 = 0 mod q: {(w*w+w+1) % q == 0}; base values distinct: {len(vals) == 4}")

# exhaustive: the strand-symmetric product (F+lam)(R+lam) with the Eisenstein encoding in F_q, k = 7 and 8
q = 2**61 - 1
w = cube_root_of_unity(q)
V = {"A": 1, "C": w, "G": (-w) % q, "T": q - 1}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
# x and lam from two independent, fixed streams (round-1 adjudication L1 F1: they were one stream); randrange is
# uniform by rejection. x != 0; the exhaustive counts do not depend on the draw unless a chance collision occurs.
x = random.Random(7).randrange(1, q); lam = random.Random(8).randrange(0, q)
for k in (7, 8):
    xs = [pow(x, t, q) for t in range(k)]
    seen = set(); classes = set()
    for tup in itertools.product("ACGT", repeat=k):
        s = "".join(tup)
        F = sum(V[c] * xs[t] for t, c in enumerate(s)) % q
        R = (-sum(V[c] * xs[k - 1 - t] for t, c in enumerate(s))) % q
        rcs = "".join(COMP[c] for c in reversed(s))
        Frc = sum(V[c] * xs[t] for t, c in enumerate(rcs)) % q
        check(Frc == R, "R must equal F of the reverse complement (complement = negation)")
        seen.add(((F + lam) * (R + lam)) % q)
        classes.add(min(s, rcs))
        if not (len(fails) < 5): break
    print(f"  Eisenstein in F_q, k={k}: distinct symmetric values {len(seen)}, double-stranded classes {len(classes)}")
    check(len(seen) == len(classes), f"k={k}: symmetric hash should be one value per class")

if fails:
    print("FAILED:", *fails[:10], sep="\n  ")
    sys.exit(1)
print("all checks passed")
