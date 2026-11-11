"""Synthetic data and brute-force references for ideas/kp5_kp11 (used by test.sh).

  gen_test.py fasta OUT SEED            random sequence with planted inverted repeats (arms 12/20/40/64) and N runs
  gen_test.py inv FASTA A GMAX          brute-force count of (i, g): s[i:i+A] == rc(s[i+A+g : i+2A+g]), no N
  gen_test.py palin FASTA K             brute-force count of windows equal to their reverse complement
  gen_test.py mismatch FASTA K D T      brute-force N_D(t) for T random targets, printed as "code count"
"""
import random
import sys

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s):
    return "".join(COMP[c] for c in reversed(s))


def read(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def ok(w):
    return all(c in COMP for c in w)


cmd = sys.argv[1]
if cmd == "fasta":
    rnd = random.Random(int(sys.argv[3]))
    s = [rnd.choice("ACGT") for _ in range(30000)]
    for a in (12, 20, 40, 64):
        for _ in range(8):
            i, g = rnd.randrange(0, 29000), rnd.randrange(0, 41)
            arm = "".join(rnd.choice("ACGT") for _ in range(a))
            seg = arm + "".join(rnd.choice("ACGT") for _ in range(g)) + rc(arm)
            s[i:i + len(seg)] = list(seg)
    for _ in range(12):
        i = rnd.randrange(0, 29990)
        s[i:i + rnd.randrange(1, 5)] = list("NNNN")[: 0] or ["N"]
    pal = "".join(rnd.choice("ACGT") for _ in range(6))
    s[100:112] = list(pal + rc(pal))  # one guaranteed 12-mer RC palindrome
    open(sys.argv[2], "w").write(">t\n" + "".join(s)[:30000] + "\n")
elif cmd == "inv":
    s, a, gmax = read(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    n = 0
    for i in range(len(s) - a + 1):
        L = s[i:i + a]
        if not ok(L):
            continue
        for g in range(gmax + 1):
            j = i + a + g
            Rt = s[j:j + a]
            if len(Rt) < a:
                break
            if ok(Rt) and L == rc(Rt):
                n += 1
    print(n)
elif cmd == "palin":
    s, k = read(sys.argv[2]), int(sys.argv[3])
    print(sum(1 for i in range(len(s) - k + 1) if ok(s[i:i + k]) and s[i:i + k] == rc(s[i:i + k])))
elif cmd == "mismatch":
    s, k, d, T = read(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    occ = [s[i:i + k] for i in range(len(s) - k + 1) if ok(s[i:i + k])]
    rnd = random.Random(7)
    for _ in range(T):
        t = "".join(rnd.choice("ACGT") for _ in range(k))
        code = 0
        for c in t:
            code = code * 4 + "ACGT".index(c)
        print(code, sum(1 for w in occ if sum(x != y for x, y in zip(w, t)) <= d))
