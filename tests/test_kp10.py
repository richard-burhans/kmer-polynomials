"""recon/kp10: KP-7 symmetric hashes in power-sum sketches recover a known difference exactly; sketches are linear."""

import numpy as np

from prototype import field as fq
from prototype.kmers import mix64
from prototype.recon import Side, power_sketch, reconcile
from recon.kp7hash import params, rc_codes, sym_hash


def sides(n=20000, d=300, seed=5):
    rng = np.random.default_rng(seed)
    base = np.unique(rng.integers(0, 1 << 62, n, dtype=np.uint64))
    extra = np.unique(rng.integers(0, 1 << 62, 2 * d, dtype=np.uint64))
    extra = np.setdiff1d(extra, base)
    a = np.union1d(base, extra[:d])
    b = np.union1d(base[: base.size - d // 2], extra[d:])
    p = params()
    mk = lambda c: Side.build(c, sym_hash(c, p), mix64(c))
    return a, b, mk(a), mk(b), p


def test_exact_recovery_power_sums():
    a, b, A, B, _ = sides()
    res = reconcile(A, B, 8, 16, 2, None, method="power")
    assert np.array_equal(res.a_minus_b, np.setdiff1d(a, b))
    assert np.array_equal(res.b_minus_a, np.setdiff1d(b, a))


def test_overflow_splits_and_still_exact():
    a, b, A, B, _ = sides(d=3000)
    res = reconcile(A, B, 6, 8, 2, None, method="power")  # ~70 differences per bucket against m = 8: must split
    assert res.splits > 0 and res.rounds > 1
    assert np.array_equal(res.a_minus_b, np.setdiff1d(a, b))
    assert np.array_equal(res.b_minus_a, np.setdiff1d(b, a))


def test_sketch_is_linear():
    a, _, A, _, p = sides()
    mask = np.random.default_rng(1).random(a.size) < 0.4
    bk = np.arange(1 << 8)
    full, nf = power_sketch(A, 8, bk, 6)
    s1, n1 = power_sketch(Side.build(a[mask], sym_hash(a[mask], p), mix64(a[mask])), 8, bk, 6)
    s2, n2 = power_sketch(Side.build(a[~mask], sym_hash(a[~mask], p), mix64(a[~mask])), 8, bk, 6)
    assert np.array_equal(full[0], (s1[0] + s2[0]) % fq.Q) and np.array_equal(full[1], (s1[1] + s2[1]) % fq.Q)
    assert np.array_equal(nf, n1 + n2)


def test_hash_is_strand_symmetric():
    c = np.random.default_rng(2).integers(0, 1 << 62, 5000, dtype=np.uint64)
    p = params()
    h1, h2 = sym_hash(c, p), sym_hash(rc_codes(c), p)
    assert np.array_equal(h1[0], h2[0]) and np.array_equal(h1[1], h2[1])


# ── the Eisenstein variant: F_q, q = 2^61 - 1, one 8-byte element per entry ──────────────────────────────────────

from recon import fq61, kp7e_hash, recon61  # noqa: E402


def sides_e(n=20000, d=300, seed=5):
    a, b, _, _, _ = sides(n, d, seed)
    p = kp7e_hash.params()
    mk = lambda c: recon61.SideE.build(c, kp7e_hash.sym_hash(c, p), mix64(c))
    return a, b, mk(a), mk(b), p


def test_eisenstein_exact_recovery():
    a, b, A, B, _ = sides_e()
    res = recon61.reconcile(A, B, 8, 16, 2)
    assert np.array_equal(res.a_minus_b, np.setdiff1d(a, b))
    assert np.array_equal(res.b_minus_a, np.setdiff1d(b, a))


def test_eisenstein_overflow_splits_and_still_exact():
    a, b, A, B, _ = sides_e(d=3000)
    res = recon61.reconcile(A, B, 6, 8, 2)
    assert res.splits > 0 and res.rounds > 1
    assert np.array_equal(res.a_minus_b, np.setdiff1d(a, b))
    assert np.array_equal(res.b_minus_a, np.setdiff1d(b, a))


def test_eisenstein_sketch_is_linear():
    a, _, A, _, p = sides_e()
    mask = np.random.default_rng(1).random(a.size) < 0.4
    bk = np.arange(1 << 8)
    full, nf = recon61.power_sketch(A, 8, bk, 6)
    mk = lambda c: recon61.SideE.build(c, kp7e_hash.sym_hash(c, p), mix64(c))
    s1, n1 = recon61.power_sketch(mk(a[mask]), 8, bk, 6)
    s2, n2 = recon61.power_sketch(mk(a[~mask]), 8, bk, 6)
    assert np.array_equal(full, fq61.add(s1, s2)) and np.array_equal(nf, n1 + n2)


def test_eisenstein_hash_strand_symmetric_and_matches_c():
    assert kp7e_hash.self_test(2000) == 0
