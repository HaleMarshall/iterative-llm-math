# iterative-llm-math

Exact verifiers, bounded searches, and reproducible certificates for **eight open
mathematical problems**, built one iteration at a time.

🔗 **Site (GitHub Pages):** https://halemarshall.github.io/iterative-llm-math/

> **These are open problems. Nothing here claims to solve one.** Following the
> methodology of *Open Problems for Iterative Local LLM Research*, the goal is to
> produce artifacts that survive contact with a checker: an independent exact
> verifier, an exact search to a stated bound, and the structural data
> (congruence classes, infeasibility certificates, small-instance
> classifications) that real progress is built from.

## The method (every problem)

```
exact generation  →  independent verification  →  congruence / obstruction mining
```

Non-negotiable rules carried from the source document:
- the **checker is a separate concern** from any generator and is the only arbiter;
- **exact arithmetic only** — no floating point in anything proof-relevant;
- **partial, bounded, and negative results count**, and are reported as exactly that.

## Portfolio & status

| # | Problem | Status |
|---|---------|--------|
| 01 | Erdős–Straus conjecture (`4/n = 1/x+1/y+1/z`) | ✅ verifier + solver, all n ≤ 2·10⁵ verified |
| 02 | Don's conjecture for standardized completely reachable DFAs | ⏳ planned (next) |
| 03 | Odd covering systems with distinct odd moduli | ⏳ planned |
| 04 | Optimal Golomb rulers | ⏳ planned |
| 05 | Even Barker sequences | ⏳ planned |
| 06 | Ryser's circulant Hadamard conjecture | ⏳ planned |
| 07 | Perfect cuboid | ⏳ planned |
| 08 | Erdős–Moser equation | ⏳ planned |

## Layout

```
problems/<name>/        exact code: verifier + solver/search
problems/<name>/tests/  pytest, witnesses re-checked against independent arithmetic
problems/<name>/results/ machine-readable run outputs (committed)
docs/                   the GitHub Pages site (shared design system)
```

## Run the tests

```bash
cd problems/erdos_straus && pytest tests/ -q
```
