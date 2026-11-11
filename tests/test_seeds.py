"""x and lam are drawn independently (round-1 adjudication L1 F1), in both Python ports of the C hashes."""

import pytest

from recon import kp7e_hash, kp7hash


def _xl(mod, p):
    if mod is kp7e_hash:
        return p.x, p.lam
    return (p.x.re, p.x.im), (p.lam.re, p.lam.im)


@pytest.mark.parametrize("mod", [kp7e_hash, kp7hash])
def test_x_depends_only_on_seed_x_and_lam_only_on_seed_lam(mod):
    x1, l1 = _xl(mod, mod.params(31, 11, 22))
    x2, l2 = _xl(mod, mod.params(31, 11, 23))
    x3, l3 = _xl(mod, mod.params(31, 12, 22))
    assert x1 == x2 and l1 != l2
    assert x1 != x3 and l1 == l3


@pytest.mark.parametrize("mod", [kp7e_hash, kp7hash])
def test_default_is_os_randomness_per_parameter(mod):
    a, b = _xl(mod, mod.params(31)), _xl(mod, mod.params(31))
    assert a[0] != b[0] and a[1] != b[1]
    xs, ls = _xl(mod, mod.params(31, 11, None)), _xl(mod, mod.params(31, 11, 22))
    assert xs[0] == ls[0] and xs[1] != ls[1]


@pytest.mark.parametrize("mod", [kp7e_hash, kp7hash])
def test_equal_seeds_give_unrelated_x_and_lam(mod):
    x, lam = _xl(mod, mod.params(31, 5, 5))
    assert x != lam
