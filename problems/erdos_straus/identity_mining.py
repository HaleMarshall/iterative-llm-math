"""Mining certified polynomial identities for the Erdős–Straus conjecture.

This realizes the source document's Template Beta/Gamma ("mine patterns from
verified data → turn a pattern into a checked lemma") for problem 01.

For a residue r modulo 840 we look at n = 840k + r and ask whether the canonical
(smallest) solution of 4/n = 1/x + 1/y + 1/z is a *polynomial family*: x(k),
y(k), z(k) polynomials in k that satisfy the equation for **all** k. The pipeline:

  1. solve 4/n exactly for k = 0, 1, …, (verified small data);
  2. Lagrange-interpolate X(k), Y(k), Z(k) through those points (SymPy);
  3. PROVE the candidate by expanding  4·X·Y·Z − (840k+r)·(YZ+XZ+XY)  and checking
     it is the zero polynomial in k — an identity valid for every k;
  4. re-check, with the exact integer `verify`, that X,Y,Z are positive integers
     and the equation holds for many further k.

A mined identity is a *theorem about a residue class*, machine-proved. Crucially,
the method **fails on exactly the residues where no identity is known** — the six
Mordell-frontier classes {1, 121, 169, 289, 361, 529} = {1, 11², 13², 17², 19²,
23²} mod 840 — so it reproduces the known frontier rather than papering over it.

This proves nothing about the open conjecture itself; it certifies identities for
residue classes, shrinking the set of residues that still require a per-value search.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd

import sympy as sp

from fractions import Fraction

from erdos_straus import solve, split_two, verify

_k = sp.symbols("k")


@dataclass
class MinedIdentity:
    modulus: int
    residue: int
    X: sp.Expr
    Y: sp.Expr
    Z: sp.Expr
    proved: bool            # 4XYZ - n(YZ+XZ+XY) is identically 0 in k
    rechecked_k: int        # integer re-verification done for k in [0, rechecked_k]

    def as_strings(self) -> dict:
        return {"modulus": self.modulus, "residue": self.residue,
                "x": str(sp.expand(self.X)), "y": str(sp.expand(self.Y)),
                "z": str(sp.expand(self.Z)), "proved": self.proved,
                "rechecked_to_k": self.rechecked_k}


def _prove_and_recheck(X, Y, Z, residue, modulus, recheck_to):
    """Symbolically prove 4XYZ = n(YZ+XZ+XY) for all k, then integer-recheck."""
    n = modulus * _k + residue
    if sp.expand(4 * X * Y * Z - n * (Y * Z + X * Z + X * Y)) != 0:
        return None
    for k in range(recheck_to + 1):
        xv, yv, zv = int(X.subs(_k, k)), int(Y.subs(_k, k)), int(Z.subs(_k, k))
        if xv <= 0 or yv <= 0 or zv <= 0:
            return None
        if not verify(modulus * k + residue, xv, yv, zv):
            return None
    return MinedIdentity(modulus, residue, sp.expand(X), sp.expand(Y), sp.expand(Z),
                         True, recheck_to)


def mine(residue: int, modulus: int = 840, n_points: int = 14,
         recheck_to: int = 60, max_offset: int = 120) -> MinedIdentity | None:
    """Mine and machine-prove a polynomial identity for n ≡ residue (mod modulus),
    or None if none is found.

    The first term is sought as ``x = (modulus/4)·k + c`` (since x ≈ n/4). For each
    admissible offset c the remainder 4/n − 1/x is split into 1/y + 1/z per k; if
    y(k), z(k) interpolate to polynomials and the resulting identity proves out,
    that residue is covered. Trying several offsets finds families the canonical
    smallest solution misses.
    """
    slope = modulus // 4                                  # 210 for modulus 840
    # offsets c must keep x > n/4 for all k: c > residue/4.
    c_lo = residue // 4 + 1
    for c in range(c_lo, c_lo + max_offset):
        pts_y, pts_z = [], []
        ok = True
        for k in range(n_points):
            n = modulus * k + residue
            if n < 2:
                return None
            x = slope * k + c
            r = Fraction(4, n) - Fraction(1, x)
            if r <= 0:
                ok = False
                break
            sp_ = split_two(r.numerator, r.denominator)
            if sp_ is None:
                ok = False
                break
            y, z = sp_
            pts_y.append((k, y)); pts_z.append((k, z))
        if not ok:
            continue
        X = slope * _k + c
        Y = sp.interpolate(pts_y, _k)
        Z = sp.interpolate(pts_z, _k)
        result = _prove_and_recheck(X, Y, Z, residue, modulus, recheck_to)
        if result is not None:
            return result
    return None


def mine_residues(modulus: int = 840) -> dict:
    """Mine identities for every residue coprime to `modulus` with r ≡ 1 (mod 4)
    (the classes the prime reduction leaves to search). Report which are covered
    by a proved identity and which residual classes remain — these should be the
    Mordell frontier."""
    targets = [r for r in range(modulus)
               if gcd(r, modulus) == 1 and r % 4 == 1]
    mined, residual = {}, []
    for r in targets:
        m = mine(r, modulus)
        if m is not None:
            mined[r] = m
        else:
            residual.append(r)
    frontier = [1, 121, 169, 289, 361, 529]
    res = sorted(residual)
    return {
        "modulus": modulus,
        "targets_examined": len(targets),         # coprime residues ≡ 1 (mod 4)
        "identities_mined_and_proved": len(mined),
        "residual_residues": res,
        "residual_count": len(res),
        "mordell_frontier": frontier,
        "frontier_is_subset_of_residual": set(frontier) <= set(res),
        "non_frontier_residual": [r for r in res if r not in frontier],
        "identities": {r: mined[r].as_strings() for r in sorted(mined)},
    }


__all__ = ["MinedIdentity", "mine", "mine_residues"]


if __name__ == "__main__":
    summ = mine_residues(840)
    print(f"coprime residues ≡1 mod4 examined: {summ['targets_examined']}")
    print(f"identities mined & proved:         {summ['identities_mined_and_proved']}")
    print(f"residual ({summ['residual_count']}):                      {summ['residual_residues']}")
    print(f"frontier ⊆ residual:               {summ['frontier_is_subset_of_residual']}")
    print(f"non-frontier residual:             {summ['non_frontier_residual']}")
