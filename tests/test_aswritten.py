"""recon/aswritten: as-written S values equal S of canonical codes and are orientation-free; F values are not; k = 63
values match the C reference on strings."""

from recon import aswritten


def test_aswritten_self_test():
    assert aswritten.self_test() == 0
