"""Erdős–Straus conjecture — exact verifier, fast solver, congruence mining.

The conjecture: for every integer n >= 2 there exist positive integers x, y, z with

        4/n = 1/x + 1/y + 1/z.

It is verified computationally past 10^17 and is widely believed, but **open**.
This module does NOT prove it. Following the source document's methodology
(exact generation → independent verification → congruence mining) it provides:

* ``verify(n, x, y, z)``  — an exact, integer-only certificate checker.
* ``solve(n)``           — a witness (x, y, z), via parametric identities with a
                           bounded exact-search fallback.
* ``verify_range(lo,hi)``— confirm every n in a range has a re-verified witness.
* ``residue_summary(m)`` — which residues mod m are handled by an O(1) identity
                           vs. fall through to search — the "easy/hard" split the
                           conjecture's congruence reductions are built on.

All arithmetic is exact integer/`Fraction` — no floating point anywhere.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt


# --- the checker (kept separate from every solver) --------------------------
def verify(n: int, x: int, y: int, z: int) -> bool:
    """True iff x, y, z are positive integers with 4/n == 1/x + 1/y + 1/z.

    Cross-multiplied so there is no rounding:
        4/n = 1/x + 1/y + 1/z  <=>  4*x*y*z == n*(y*z + x*z + x*y).
    """
    if n < 2 or x <= 0 or y <= 0 or z <= 0:
        return False
    return 4 * x * y * z == n * (y * z + x * z + x * y)


# --- O(1) parametric identities (each provably exact; checked by `verify`) ---
def identity(n: int) -> tuple[int, int, int] | None:
    """Return a witness from a classical identity when n's residue allows it.

    These cover every n that is even, divisible by 3, or ≡ 3 (mod 4). The only
    residues that fall through are n ≡ 1 (mod 4) coprime to 3 — exactly the class
    the conjecture genuinely rests on.
    """
    if n % 2 == 0:                                   # n = 2m
        m = n // 2
        return (m, m + 1, m * (m + 1))               # 1/m + 1/(m+1) + 1/(m(m+1)) = 2/m
    if n % 3 == 0:                                   # n = 3m
        m = n // 3
        return (m, 3 * m + 1, 3 * m * (3 * m + 1))   # 1/m + 1/(3m) split into 3 terms
    if n % 4 == 3:                                   # n = 4k - 1
        k = (n + 1) // 4
        nk = n * k
        return (k, nk + 1, nk * (nk + 1))            # 1/k + 1/(nk) split into 3 terms
    return None


# --- exact two-term splitter via Simon's favorite factoring -----------------
def _divisors_up_to_sqrt_of_square(q: int):
    """Yield divisors A of q*q with A <= q (i.e. A <= sqrt(q*q))."""
    # factor q, then build divisors of q^2 = prod p^(2a)
    qq = q
    factors: dict[int, int] = {}
    d = 2
    while d * d <= qq:
        while qq % d == 0:
            factors[d] = factors.get(d, 0) + 1
            qq //= d
        d += 1 if d == 2 else 2
    if qq > 1:
        factors[qq] = factors.get(qq, 0) + 1
    divs = [1]
    for p, a in factors.items():
        divs = [dd * p ** e for dd in divs for e in range(2 * a + 1)]
    limit = q
    for dd in sorted(set(divs)):
        if dd > limit:
            break
        yield dd


def split_two(p: int, q: int) -> tuple[int, int] | None:
    """Positive integers y <= z with 1/y + 1/z == p/q, or None.

    From 1/y + 1/z = p/q one gets (p*y - q)(p*z - q) = q^2, so y and z come from
    a divisor pair of q^2 — far faster than scanning y up to 2q/p.
    """
    g = gcd(p, q)
    p, q = p // g, q // g
    if p <= 0:
        return None
    qq = q * q
    for a in _divisors_up_to_sqrt_of_square(q):
        b = qq // a
        if (a + q) % p == 0 and (b + q) % p == 0:
            y, z = (a + q) // p, (b + q) // p
            if 0 < y <= z:
                return (y, z)
    return None


# --- the solver -------------------------------------------------------------
def solve(n: int, max_x: int | None = None) -> tuple[int, int, int] | None:
    """A witness (x, y, z) with x <= y <= z, or None.

    Tries the O(1) identities first; otherwise scans the smallest term
    x in (n/4, 3n/4] and splits the remainder exactly with `split_two`.
    """
    if n < 2:
        return None
    quick = identity(n)
    if quick is not None:
        x, y, z = sorted(quick)
        return (x, y, z)

    target = Fraction(4, n)
    x_lo = n // 4 + 1
    x_hi = (3 * n) // 4
    if max_x is not None:
        x_hi = min(x_hi, x_lo + max_x)
    for x in range(x_lo, x_hi + 1):
        r = target - Fraction(1, x)
        if r <= 0:
            continue
        split = split_two(r.numerator, r.denominator)
        if split is not None:
            y, z = split
            if y < x:           # keep x <= y <= z by retrying with y as anchor? simply sort
                return tuple(sorted((x, y, z)))
            return (x, y, z)
    return None


@dataclass
class RangeReport:
    lo: int
    hi: int
    checked: int
    solved: int
    by_identity: int
    failures: list

    @property
    def all_solved(self) -> bool:
        return not self.failures


def verify_range(lo: int = 2, hi: int = 10_000) -> RangeReport:
    """Find and independently re-verify a witness for every n in [lo, hi]."""
    failures, solved, by_id = [], 0, 0
    for n in range(max(2, lo), hi + 1):
        if identity(n) is not None:
            by_id += 1
        w = solve(n)
        if w is not None and verify(n, *w):
            solved += 1
        else:
            failures.append(n)
    return RangeReport(lo, hi, hi - max(2, lo) + 1, solved, by_id, failures)


def residue_summary(modulus: int = 12) -> dict:
    """Which residues mod `modulus` an O(1) identity handles, and which fall
    through to search. The fall-through classes are where the conjecture's real
    difficulty lives (n ≡ 1 mod 4, coprime to 3)."""
    handled, fell_through = [], []
    for r in range(modulus):
        n = r if r >= 2 else r + modulus
        (handled if identity(n) is not None else fell_through).append(r)
    return {"modulus": modulus, "handled_by_identity": handled,
            "fall_through_to_search": fell_through}


__all__ = ["verify", "identity", "split_two", "solve", "verify_range",
           "RangeReport", "residue_summary"]


if __name__ == "__main__":
    rep = verify_range(2, 20_000)
    print(f"[{rep.lo},{rep.hi}]: {rep.solved}/{rep.checked} solved, "
          f"{rep.by_identity} by O(1) identity, failures={len(rep.failures)}")
    print("residues mod 12:", residue_summary(12))
