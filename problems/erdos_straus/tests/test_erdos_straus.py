"""Tests for the Erdős–Straus tooling.

The checker is the ground truth: every witness any solver returns is re-verified
against exact `Fraction` arithmetic AND the independent `verify`. We also pin the
parametric identities and the bounded-range "every n is solvable" property — the
exact, falsifiable artifacts the source document calls for (not a proof).
"""

import os
import sys
from fractions import Fraction

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from erdos_straus import (
    hard_prime_residues,
    identity,
    is_hard_prime,
    residue_summary,
    solve,
    split_two,
    verify,
    verify_primes_up_to,
    verify_range,
    verify_up_to,
    witness,
)


def test_verify_accepts_a_true_certificate():
    assert verify(2, 1, 2, 2)            # 1/1? no: 1/2+1/2... check below
    assert verify(5, 2, 4, 20)           # 1/2 + 1/4 + 1/20 = 4/5


@pytest.mark.parametrize("cert", [(5, 2, 4, 21), (5, 0, 4, 20), (5, -2, 4, 20), (1, 1, 1, 1)])
def test_verify_rejects_bad_certificates(cert):
    assert not verify(*cert)


def test_solver_witnesses_are_exactly_four_over_n():
    for n in list(range(2, 2000)) + [9241, 9973, 100003]:
        x, y, z = solve(n)
        assert Fraction(1, x) + Fraction(1, y) + Fraction(1, z) == Fraction(4, n)
        assert verify(n, x, y, z)
        assert x <= y <= z


def test_every_n_in_a_range_is_solvable():
    rep = verify_range(2, 5000)
    assert rep.all_solved
    assert rep.solved == rep.checked
    assert rep.failures == []


@pytest.mark.parametrize("n", [2, 4, 100, 2 * 9973])         # even
def test_even_identity(n):
    assert identity(n) is not None and verify(n, *identity(n))


@pytest.mark.parametrize("n", [3, 9, 21, 3 * 9973])          # divisible by 3
def test_div3_identity(n):
    assert identity(n) is not None and verify(n, *identity(n))


@pytest.mark.parametrize("n", [7, 11, 19, 23, 4 * 1000 - 1])  # n ≡ 3 (mod 4)
def test_three_mod_four_identity(n):
    assert identity(n) is not None and verify(n, *identity(n))


def test_hard_residues_fall_through_to_search():
    # n ≡ 1, 5 (mod 12) have no O(1) identity here but must still be solvable.
    summ = residue_summary(12)
    assert set(summ["fall_through_to_search"]) == {1, 5}
    for n in (13, 17, 25, 29, 73, 97):       # ≡ 1 mod 4 and coprime to 3
        assert identity(n) is None
        assert verify(n, *solve(n))


@pytest.mark.parametrize("p,q", [(2, 5), (3, 7), (1, 6), (4, 9), (5, 12)])
def test_split_two_is_exact(p, q):
    res = split_two(p, q)
    if res is not None:
        y, z = res
        assert Fraction(1, y) + Fraction(1, z) == Fraction(p, q)
        assert 0 < y <= z


def test_prime_reduction_scaling_is_an_exact_identity():
    # The crux lemma: 4/m = 1/a+1/b+1/c  =>  4/(m*t) = 1/(at)+1/(bt)+1/(ct).
    for m in (5, 7, 13, 67):
        a, b, c = solve(m)
        for t in (1, 2, 3, 10, 97):
            assert verify(m * t, a * t, b * t, c * t)


def test_witness_solves_via_prime_reduction():
    for n in list(range(2, 3000)) + [25, 65, 85, 221, 999983, 1000000]:
        w = witness(n, {})
        assert w is not None and verify(n, *w)
        assert w[0] <= w[1] <= w[2]


def test_only_primes_one_mod_four_are_hard():
    assert is_hard_prime(13) and is_hard_prime(97)
    assert not is_hard_prime(7)        # ≡ 3 mod 4 -> identity
    assert not is_hard_prime(2) and not is_hard_prime(3)
    assert not is_hard_prime(25)       # composite (covered by reducing to 5)
    # every value that lacked an O(1) identity AND is prime must be ≡ 1 mod 4
    for n in range(5, 2000):
        if identity(n) is None and witness(n, {}) and _is_prime(n):
            assert n % 4 == 1


def _is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def test_verify_primes_certifies_all_n_below_bound():
    # all primes <= N solved  ==>  all n <= N solvable (by the reduction)
    rp = verify_primes_up_to(50_000)
    assert rp.all_solved
    an = verify_up_to(50_000)
    assert an.all_solved and an.solved == an.checked


def test_mordell_frontier_residues_are_the_six_squares():
    d = hard_prime_residues(100_000, 840)
    assert d["mordell_frontier_residues"] == [1, 121, 169, 289, 361, 529]
    # = {1, 11^2, 13^2, 17^2, 19^2, 23^2} mod 840
    assert sorted({r * r % 840 for r in (1, 11, 13, 17, 19, 23)}) == d["mordell_frontier_residues"]


def test_split_two_matches_a_brute_scan():
    # independent oracle: brute-force the two-term split and compare existence.
    def brute(p, q):
        for y in range(q // p + 1, 2 * q // p + 1):
            s = Fraction(p, q) - Fraction(1, y)
            if s > 0 and s.numerator == 1 and s.denominator >= y:
                return True
        return False
    for p in range(1, 6):
        for q in range(1, 40):
            assert (split_two(p, q) is not None) == brute(p, q), (p, q)
