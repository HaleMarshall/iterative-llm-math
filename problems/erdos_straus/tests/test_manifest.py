"""Tests for the theorem-manifest certificates.

These independently re-validate the per-identity proof bundle (polynomial
identity, integrality via forward differences, positivity via coefficients,
concrete witnesses) without rebuilding the whole manifest (which is slow).
"""

import os
import sys

import pytest
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from erdos_straus import verify
from identity_mining import mine
from theorem_manifest import _forward_differences, _identity_record, _why_eight_remain

K = sp.symbols("k")


def test_forward_differences_detect_integer_valued_polynomials():
    # binomial C(k,2) = k(k-1)/2 is integer-valued though it has a 1/2 coefficient
    P = K * (K - 1) / 2
    diffs = _forward_differences(P)
    assert all(isinstance(d, int) for d in diffs)
    # and it really is integer at many k
    assert all(float(P.subs(K, k)).is_integer() for k in range(20))


@pytest.mark.parametrize("r", [13, 17, 29, 53])
def test_identity_record_certificates_are_all_valid(r):
    rec = _identity_record(mine(r), 840)
    p = rec["proofs"]
    # polynomial identity
    assert p["polynomial_identity"]["is_zero_polynomial"]
    # integrality: every forward-difference list is all-integer
    for term in ("x", "y", "z"):
        assert p["integrality"]["per_term"][term]["all_integer"]
    # positivity: all coefficients >= 0 and constant term > 0
    for term in ("x", "y", "z"):
        assert p["positivity"]["per_term"][term]["all_nonnegative"]
        assert p["positivity"]["per_term"][term]["constant_term_positive"]
    # the 5 witnesses verify against the independent checker
    assert len(rec["first_5_witnesses"]) == 5
    for w in rec["first_5_witnesses"]:
        assert w["verified"] and verify(w["n"], w["x"], w["y"], w["z"])


def test_witnesses_match_the_polynomial_evaluated():
    rec = _identity_record(mine(13), 840)
    X = sp.sympify(rec["x_of_k"]); Y = sp.sympify(rec["y_of_k"]); Z = sp.sympify(rec["z_of_k"])
    for w in rec["first_5_witnesses"]:
        assert w["x"] == int(X.subs(K, w["k"]))
        assert w["y"] == int(Y.subs(K, w["k"]))
        assert w["z"] == int(Z.subs(K, w["k"]))


def test_diagnosis_of_the_eight_classes_is_recorded():
    d = _why_eight_remain()
    assert set(d) == {"interpolation_degree", "ansatz",
                      "missing_known_identities", "richer_families_required"}
    assert "NOT the cause" in d["interpolation_degree"]
    assert "THE cause" in d["ansatz"]
