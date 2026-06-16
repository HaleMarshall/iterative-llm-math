# Erdős–Straus conjecture

**Open.** For every integer n ≥ 2, does `4/n = 1/x + 1/y + 1/z` have a solution in
positive integers? Believed true; verified past 10¹⁷. This is **not** a proof — it
is an exact verifier, a fast solver, and the congruence structure of the problem.

## API (`erdos_straus.py`)

```python
from erdos_straus import verify, solve, verify_range, residue_summary

verify(5, 2, 4, 20)        # True  — 1/2 + 1/4 + 1/20 = 4/5
solve(9973)                # (2494, 8290888, 103108227451928)
verify_range(2, 200000).all_solved   # True
residue_summary(12)        # which residues need search vs. an O(1) identity
```

- **`verify(n,x,y,z)`** — the checker: `4·xyz == n·(yz+xz+xy)`, integer-only.
- **`solve(n)`** — O(1) identities (n even / ÷3 / ≡ 3 mod 4) + a bounded exact
  search using the `(py−q)(pz−q) = q²` split (divisors of `q²`, not every denominator).
- **`witness(n)`** — solves via the **prime reduction**: reduce to the smallest
  prime factor p, solve `4/p`, scale by `n/p`. Search is invoked (and memoized)
  only on primes p ≡ 1 (mod 4).
- **`verify_up_to` / `verify_primes_up_to`** — bounded verification (all n, or all
  primes ⇒ all n by the reduction).
- **`hard_prime_residues`** — distribution of hard primes over residues mod 840,
  surfacing the Mordell frontier `{1, 121, 169, 289, 361, 529}`.

## The structure (what makes this "complete")

1. **Reduction to primes** (exact, unit-tested): `m | n` and `4/m = 1/a+1/b+1/c`
   ⟹ `4/n = 1/(a·n/m)+1/(b·n/m)+1/(c·n/m)`. So the conjecture collapses to primes.
2. Primes `2, 3` and `p ≡ 3 (mod 4)` have O(1) identities ⟹ only **primes
   `p ≡ 1 (mod 4)`** ever need search.
3. **Mordell frontier**: identities exist for all residues mod 840 except the six
   squares `{1, 11², 13², 17², 19², 23²}`. So the whole conjecture rests on primes
   in those six classes.

## Results

Every n ∈ [2, 10⁶] has a re-verified witness; equivalently all 78,498 primes ≤ 10⁶
are certified (39,175 hard primes searched; 2,370 in the frontier classes). See
[`results/erdos_straus_results.json`](results/erdos_straus_results.json).
*(The literature has verified the conjecture past 10¹⁷; the value here is the
tooling and reduction structure, not the record. The conjecture remains open.)*

## Tests

```bash
pytest tests/ -q     # 32 tests
```

Witnesses are re-checked against exact `Fraction` arithmetic and the independent
checker; the prime-reduction scaling is verified as an exact identity; `split_two`
is cross-checked against a brute scan; the frontier residues are pinned to the six squares.
