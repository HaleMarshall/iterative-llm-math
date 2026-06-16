"""Tests for the polynomial-identity miner.

A mined identity is only accepted if (a) it expands to the zero polynomial in k
and (b) it integer-re-verifies. These tests independently confirm mined
identities hold far beyond the sample range, and — crucially — that the miner
returns NO identity for the open Mordell-frontier residues (no false positives).
"""

import os
import sys

import pytest
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from erdos_straus import verify
from identity_mining import mine

K = sp.symbols("k")
MINEABLE = [13, 17, 29, 37, 53, 61]            # coprime, ≡1 mod4, not frontier
FRONTIER = [1, 121, 169, 289, 361, 529]        # {1,11²,13²,17²,19²,23²} mod 840


@pytest.mark.parametrize("r", MINEABLE)
def test_mine_returns_a_proved_identity(r):
    m = mine(r)
    assert m is not None and m.proved
    # re-prove symbolically, independently of the miner's own check
    n = 840 * K + r
    assert sp.expand(4 * m.X * m.Y * m.Z - n * (m.Y * m.Z + m.X * m.Z + m.X * m.Y)) == 0


@pytest.mark.parametrize("r", MINEABLE)
def test_mined_identity_holds_far_beyond_the_sample(r):
    m = mine(r)
    for k in (0, 1, 7, 50, 200, 1000):          # sample range was only k<14
        x, y, z = (int(m.X.subs(K, k)), int(m.Y.subs(K, k)), int(m.Z.subs(K, k)))
        assert x > 0 and y > 0 and z > 0
        assert verify(840 * k + r, x, y, z)


@pytest.mark.parametrize("r", FRONTIER)
def test_no_false_identity_on_the_open_frontier(r):
    # The frontier classes are where no identity is known; the miner must not
    # fabricate one.
    assert mine(r) is None


def test_first_term_has_the_expected_linear_shape():
    m = mine(13)
    # x is linear in k with slope 210 (= 840/4, since x ≈ n/4)
    poly = sp.Poly(m.X, K)
    assert poly.degree() == 1
    assert poly.all_coeffs()[0] == 210
