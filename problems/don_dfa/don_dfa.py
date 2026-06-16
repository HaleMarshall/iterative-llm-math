"""Don's conjecture for standardized completely reachable DFAs — exact tooling.

Setup (Zhu, arXiv:2402.19089; Casas–Volkov, arXiv:2205.09404). A DFA has states
Z_n = {0,…,n-1} and two letters {a, b}. A word w acts on the *whole* state set by
the IMAGE action: Q·w = δ(Q, w) = { δ(q, w) : q ∈ Q }. A subset S has a "reaching
word" w if Q·w = S, and the DFA is **completely reachable** (CRA) if every
non-empty subset has a reaching word.

**Don's conjecture (open):** in every n-state CRA, every k-element subset has a
reaching word of length ≤ n(n−k).

**Zhu's theorem (proved):** for *standardized* DFAs the weaker bound
n(n−k)+n−1 always holds. Whether the stronger n(n−k) bound holds for standardized
DFAs is **open**.

**Standardized** (the normal form used here): b is the cycle q ↦ q+1 (mod n); a is
a deficiency-1 map with excl(a) = {0} (no state maps to 0) and dupl(a) = {a(0)}
(exactly one collision, at a(0)); and the automaton is completely reachable.

This module is exact (integer/set arithmetic; BFS over the subset lattice). It does
NOT prove Don's conjecture; it verifies the bounds on bounded families and searches
for counterexamples. Subsets are encoded as bitmasks for speed.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from itertools import permutations


def image(mask: int, t: tuple[int, ...]) -> int:
    """Image of the state subset `mask` under transition table `t` (a bitmask)."""
    out = 0
    q = 0
    m = mask
    while m:
        if m & 1:
            out |= 1 << t[q]
        m >>= 1
        q += 1
    return out


def reaching_lengths(a: tuple[int, ...], b: tuple[int, ...], n: int) -> dict[int, int]:
    """Shortest reaching-word length for every subset reachable from Q, by BFS
    over the image action. Keys are bitmasks; Q (full set) has length 0."""
    full = (1 << n) - 1
    dist = {full: 0}
    dq = deque([full])
    while dq:
        s = dq.popleft()
        d = dist[s]
        for t in (a, b):
            ns = image(s, t)
            if ns not in dist:
                dist[ns] = d + 1
                dq.append(ns)
    return dist


def is_completely_reachable(a: tuple[int, ...], b: tuple[int, ...], n: int) -> bool:
    """True iff every non-empty subset of states is reachable from Q."""
    return len(reaching_lengths(a, b, n)) == (1 << n) - 1


def don_bound(n: int, k: int) -> int:
    return n * (n - k)


def zhu_bound(n: int, k: int) -> int:
    return n * (n - k) + n - 1


@dataclass
class BoundReport:
    n: int
    reachable_subsets: int
    completely_reachable: bool
    don_violations: list      # (subset_bitmask, k, shortest_len) with len > n(n-k)
    zhu_violations: list      # (subset_bitmask, k, shortest_len) with len > n(n-k)+n-1
    max_len_over_don_ratio: float

    @property
    def don_holds(self) -> bool:
        return not self.don_violations

    @property
    def zhu_holds(self) -> bool:
        return not self.zhu_violations


def check_bounds(a: tuple[int, ...], b: tuple[int, ...], n: int) -> BoundReport:
    """Check Don's and Zhu's bounds on every reachable subset."""
    dist = reaching_lengths(a, b, n)
    don_bad, zhu_bad = [], []
    worst = 0.0
    for s, d in dist.items():
        k = bin(s).count("1")
        if k == n:
            continue                      # Q itself, length 0
        db = don_bound(n, k)
        if d > db:
            don_bad.append((s, k, d))
        if d > zhu_bound(n, k):
            zhu_bad.append((s, k, d))
        if db > 0:
            worst = max(worst, d / db)
    return BoundReport(n, len(dist), len(dist) == (1 << n) - 1,
                       don_bad, zhu_bad, round(worst, 4))


# --- standardized DFA enumeration -------------------------------------------
def cyclic_b(n: int) -> tuple[int, ...]:
    return tuple((q + 1) % n for q in range(n))


def standardized_a_maps(n: int):
    """Yield every transition table `a` in the standardized normal form:
    excl(a) = {0} (no state maps to 0), dupl(a) = {a(0)} (exactly one collision,
    at a(0)). Equivalently: a is a surjection Z_n → Z_n\\{0} whose only repeated
    value is a(0)."""
    nonzero = list(range(1, n))
    for v in nonzero:                              # v = a(0), the doubled value
        for partner in range(1, n):                # the state != 0 with a(partner)=v
            # remaining states (not 0, not partner) bijection onto nonzero\{v}
            rest_states = [q for q in range(1, n) if q != partner]
            targets = [t for t in nonzero if t != v]
            if len(rest_states) != len(targets):
                continue
            for perm in permutations(targets):
                a = [0] * n
                a[0] = v
                a[partner] = v
                for q, t in zip(rest_states, perm):
                    a[q] = t
                yield tuple(a)


def enumerate_standardized(n: int):
    """Yield (a, b) for every standardized DFA on n states that is, in addition,
    completely reachable (the defining membership condition)."""
    b = cyclic_b(n)
    for a in standardized_a_maps(n):
        if is_completely_reachable(a, b, n):
            yield a, b


__all__ = ["image", "reaching_lengths", "is_completely_reachable",
           "don_bound", "zhu_bound", "BoundReport", "check_bounds",
           "cyclic_b", "standardized_a_maps", "enumerate_standardized"]
