"""Tests for the Don's-conjecture tooling.

The independent verifier (`verify.py`) is the ground truth: the BFS lengths from
the generator must match the lengths of explicit reaching words it produces, and
every reaching word must actually reach its subset. We pin the standardized
enumeration, complete-reachability detection, and the headline bounded result
(Don's and Zhu's bounds hold for every standardized CRA DFA up to n=6).
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import don_dfa as dd
import verify as vf


def test_image_action():
    # a maps 0,2->1 ; 1->2 on 3 states
    a = (1, 2, 1)
    assert dd.image(0b111, a) == 0b110         # {0,1,2} -> {1,2}
    assert dd.image(0b001, a) == 0b010         # {0} -> {1}


def test_a_known_small_dfa_is_completely_reachable():
    n, b = 3, dd.cyclic_b(3)
    a = (1, 2, 1)
    assert dd.is_completely_reachable(a, b, n)


def test_non_cra_example_is_detected():
    # b a permutation and a = identity-ish that cannot shrink Q -> not CRA
    n = 3
    b = dd.cyclic_b(3)
    a = (0, 1, 2)                                # identity: image never shrinks
    assert not dd.is_completely_reachable(a, b, n)


def test_bfs_lengths_match_independent_reaching_words():
    n, b = 4, dd.cyclic_b(4)
    a = (1, 2, 3, 1)
    if not dd.is_completely_reachable(a, b, n):
        pytest.skip("need a CRA example")
    dist = dd.reaching_lengths(a, b, n)
    for mask, d in dist.items():
        w = vf.shortest_reaching_word(a, b, n, mask)
        assert w is not None
        assert len(w) == d                       # generator length == verifier word length
        assert vf.reaches(a, b, n, w, mask)       # the word really reaches it


@pytest.mark.parametrize("n,expected_count", [(3, 4), (4, 16), (5, 96), (6, 552)])
def test_standardized_enumeration_counts(n, expected_count):
    assert sum(1 for _ in dd.enumerate_standardized(n)) == expected_count


def test_standardized_maps_have_the_right_structure():
    n = 6
    for a in dd.standardized_a_maps(n):
        assert 0 not in a                        # excl(a) = {0}: nothing maps to 0
        from collections import Counter
        c = Counter(a)
        dup = [s for s, cnt in c.items() if cnt >= 2]
        assert dup == [a[0]]                      # the only duplicated value is a(0)
        assert c[a[0]] == 2                       # deficiency exactly 1


@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_don_and_zhu_bounds_hold_for_all_standardized_cra(n):
    for a, b in dd.enumerate_standardized(n):
        rep = dd.check_bounds(a, b, n)
        assert rep.completely_reachable
        assert rep.zhu_holds                      # proved theorem (sanity)
        assert rep.don_holds                      # the OPEN conjecture, here verified


def test_don_bound_is_tight_somewhere():
    # the bound n(n-k) is achieved with equality (ratio exactly 1.0), not loose.
    n, b = 5, dd.cyclic_b(5)
    a = (1, 1, 2, 3, 4)
    cert = vf.certify_don_bound(a, b, n)
    assert cert["don_bound_holds"]
    assert cert["worst_case_vs_don_bound"]["ratio"] == 1.0


def test_the_n5_singleton_certificate():
    # reaching {0} needs a word of length exactly n(n-1) = 20, independently checked.
    n, b = 5, dd.cyclic_b(5)
    a = (1, 1, 2, 3, 4)
    w = vf.shortest_reaching_word(a, b, n, 0b00001)
    assert len(w) == 20 == dd.don_bound(5, 1)
    assert vf.reaches(a, b, n, w, 0b00001)


def test_verifier_rejects_a_wrong_word():
    n, b = 4, dd.cyclic_b(4)
    a = (1, 2, 3, 1)
    assert not vf.reaches(a, b, n, "a", (1 << n) - 1)   # one 'a' shrinks Q, ≠ Q
