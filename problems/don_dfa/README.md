# Don's conjecture for standardized completely reachable DFAs

**Open** (for standardized DFAs). In an n-state completely reachable automaton,
can every k-element subset be reached by a word of length ≤ `n(n−k)`? This is
**not** a proof — it is exact enumeration, an independent reaching-word verifier,
and a bounded counterexample search.

## Setup (Zhu, arXiv:2402.19089; Casas–Volkov, arXiv:2205.09404)

- **Image action:** `Q·w = {δ(q,w) : q ∈ Q}`; `w` reaches `S` iff `Q·w = S`.
- **Completely reachable (CRA):** every non-empty subset has a reaching word.
- **Standardized:** `b = ` cycle `q↦q+1 (mod n)`; `a` deficiency-1 with
  `excl(a)={0}`, `dupl(a)={a(0)}` (one merge); and the automaton is CRA.
- **Don's conjecture (open):** length ≤ `n(n−k)`. **Zhu (proved):** ≤ `n(n−k)+n−1`.

## API

- `don_dfa.py` — `enumerate_standardized(n)`, `reaching_lengths`,
  `is_completely_reachable`, `check_bounds` (Don + Zhu).
- `verify.py` (separate arbiter) — `shortest_reaching_word`, `reaches`,
  `certify_don_bound` (explicit, re-checked reaching words).

```python
import don_dfa as dd, verify as vf
for a, b in dd.enumerate_standardized(6):
    assert dd.check_bounds(a, b, 6).don_holds       # Don's bound holds
```

## Result

Across **all 39,164 standardized CRA DFAs for n = 3…8**: **0** Don-bound
violations and **0** Zhu-bound violations. So **Don's conjecture holds for every
standardized completely reachable DFA with n ≤ 8**, and `n(n−k)` is **tight**
(achieved with equality, never exceeded). Certificate (n=5, `a=(1,1,2,3,4)`):
reaching `{0}` needs `abbbbabbbbabbbbabbbb`, length `20 = n(n−k)`, independently
verified. See [`results/don_dfa_results.json`](results/don_dfa_results.json).

## Tests

```bash
pytest tests/ -q     # 16 tests
```

The BFS lengths are cross-checked against explicit reaching words produced and
verified by the independent `verify.py`; the standardized enumeration structure
and counts are pinned.
