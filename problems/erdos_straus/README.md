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
- **`solve(n)`** — O(1) identities for n even / divisible by 3 / ≡ 3 (mod 4); a
  bounded exact search for the residual class (n ≡ 1 mod 4, coprime to 3). The
  two-term split uses `(py−q)(pz−q) = q²`, so it scans divisors of `q²` rather
  than every candidate denominator.
- **`verify_range` / `residue_summary`** — bounded verification and the easy/hard split.

## Results

All n ∈ [2, 200000] have a re-verified witness: ⅚ via an O(1) identity, ⅙ via
bounded search — matching the residue analysis (only n ≡ 1, 5 mod 12 fall
through). See [`results/erdos_straus_results.json`](results/erdos_straus_results.json).

## Tests

```bash
pytest tests/ -q     # 27 tests
```

Witnesses are re-checked against exact `Fraction` arithmetic and the independent
checker; `split_two` is cross-checked against a brute-force scan.
