"""KP-9/KP-10 prototype tests: exact recovery of a known difference, overflow detection, Newton's identities."""

import numpy as np
import pytest

from prototype import field as fq
from prototype import kmers, recon


def sides(codesA, codesB, x0):
    hA = kmers.poly_hash(codesA, x0)
    hB = kmers.poly_hash(codesB, x0)
    return (recon.Side.build(codesA, hA, kmers.mix64(codesA)),
            recon.Side.build(codesB, hB, kmers.mix64(codesB)))


def random_sets(rng, n_common, n_a, n_b):
    pool = np.unique(rng.integers(0, 1 << 62, n_common + n_a + n_b + 1000, dtype=np.uint64))
    rng.shuffle(pool)
    common, onlyA, onlyB = pool[:n_common], pool[n_common:n_common + n_a], pool[n_common + n_a:n_common + n_a + n_b]
    return np.concatenate([common, onlyA]), np.concatenate([common, onlyB]), np.sort(onlyA), np.sort(onlyB)


@pytest.mark.parametrize("method", ["eval", "power"])
@pytest.mark.parametrize("n_a,n_b", [(0, 0), (5, 0), (0, 7), (13, 11), (40, 3)])
def test_known_difference_recovered_exactly(method, n_a, n_b):
    rng = np.random.default_rng(1 + n_a * 100 + n_b)
    A, B, onlyA, onlyB = random_sets(rng, 3000, n_a, n_b)
    x0 = fq.random_elements(rng, 1)
    x0 = (x0[0][0], x0[1][0])
    sa, sb = sides(A, B, x0)
    m, c = 32, 4
    z = fq.random_elements(rng, m + c)
    res = recon.reconcile(sa, sb, bits=2, m=m, c=c, z=z, method=method)
    assert np.array_equal(res.a_minus_b, onlyA)
    assert np.array_equal(res.b_minus_a, onlyB)


@pytest.mark.parametrize("method", ["eval", "power"])
def test_overflow_is_detected_and_split(method):
    """One bucket, a difference of 60 > m = 16: decoding must refuse it and splitting must still recover it exactly."""
    rng = np.random.default_rng(7)
    A, B, onlyA, onlyB = random_sets(rng, 500, 35, 25)
    x0 = fq.random_elements(rng, 1)
    x0 = (x0[0][0], x0[1][0])
    sa, sb = sides(A, B, x0)
    m, c = 16, 4
    z = fq.random_elements(rng, m + c)
    res = recon.reconcile(sa, sb, bits=0, m=m, c=c, z=z, method=method)
    assert res.splits >= 1 and res.rounds >= 2
    assert np.array_equal(res.a_minus_b, onlyA)
    assert np.array_equal(res.b_minus_a, onlyB)


def test_decode_refuses_overflow_directly():
    rng = np.random.default_rng(11)
    A, B, _, _ = random_sets(rng, 200, 30, 30)
    x0 = (np.int64(12345), np.int64(678))
    sa, sb = sides(A, B, x0)
    m, c = 16, 4
    z = fq.random_elements(rng, m + c)
    skA, nA = recon.eval_sketch(sa, 0, np.array([0]), z)
    skB, nB = recon.eval_sketch(sb, 0, np.array([0]), z)
    out = recon.decode_eval((skA[0][0], skA[1][0]), (skB[0][0], skB[1][0]), nA[0], nB[0], z, m, sa.h, sb.h)
    assert out is None


def test_newton_identities_give_the_elementary_symmetric_polynomials():
    """Power sums of a set, through Newton, give prod (1 - a t): compare with direct expansion."""
    rng = np.random.default_rng(3)
    a = fq.random_elements(rng, 9)
    m = 9
    S_re = np.empty(m, np.int64)
    S_im = np.empty(m, np.int64)
    x = a
    for j in range(m):
        S_re[j], S_im[j] = x[0].sum() % fq.Q, x[1].sum() % fq.Q
        x = fq.mul(x, a)
    r = recon.newton_series((S_re, S_im), m)
    # direct: multiply out prod (1 - a_i t); new[n] = old[n] - a_i * old[n-1]
    coef = (np.array([1], np.int64), np.array([0], np.int64))
    for i in range(a[0].size):
        old = (np.concatenate([coef[0], [0]]), np.concatenate([coef[1], [0]]))
        prev = (np.concatenate([[0], coef[0]]), np.concatenate([[0], coef[1]]))
        ai = (np.full(prev[0].size, a[0][i]), np.full(prev[0].size, a[1][i]))
        coef = fq.sub(old, fq.mul(ai, prev))
    assert np.array_equal(r[0], coef[0][:m + 1]) and np.array_equal(r[1], coef[1][:m + 1])


def test_canonical_kmers_are_strand_invariant():
    seq = b"ACGTTGCANNACGGTACCATGGTACGTAGCTAGCTAGGATCCAGTAC" * 3
    rc = seq[::-1].translate(bytes.maketrans(b"ACGTN", b"TGCAN"))
    k = 11
    assert np.array_equal(kmers.kmer_set(kmers.encode(seq), k), kmers.kmer_set(kmers.encode(rc), k))


def test_field_inverse():
    rng = np.random.default_rng(5)
    a = fq.random_elements(rng, 1000)
    one = fq.mul(a, fq.inv(a))
    nz = ~fq.is_zero(a)
    assert np.all(one[0][nz] == 1) and np.all(one[1][nz] == 0)
