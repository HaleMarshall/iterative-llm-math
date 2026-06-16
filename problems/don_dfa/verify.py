"""Independent verifier for Don's-conjecture reaching words.

Kept deliberately separate from the BFS in `don_dfa.py` (the "generator"): this
module re-derives the image action from scratch and is the sole arbiter that a
claimed reaching word actually reaches its subset. It also produces explicit
shortest reaching *words* (not just lengths) as certificates, and checks them.
"""

from __future__ import annotations

from collections import deque


def _image(mask: int, t: tuple[int, ...]) -> int:
    out = 0
    for q in range(len(t)):
        if mask >> q & 1:
            out |= 1 << t[q]
    return out


def apply_word(a: tuple[int, ...], b: tuple[int, ...], n: int, word: str) -> int:
    """Q · word under the image action, returned as a bitmask."""
    s = (1 << n) - 1
    for c in word:
        if c == "a":
            s = _image(s, a)
        elif c == "b":
            s = _image(s, b)
        else:
            raise ValueError(f"word letters must be 'a' or 'b', got {c!r}")
    return s


def reaches(a, b, n: int, word: str, target_mask: int) -> bool:
    """True iff applying `word` to Q yields exactly `target_mask`."""
    return apply_word(a, b, n, word) == target_mask


def shortest_reaching_word(a, b, n: int, target_mask: int) -> str | None:
    """A shortest word w with Q·w == target_mask (independent BFS with parent
    pointers), or None if unreachable. Letter 'a' is tried before 'b' so the
    result is the lexicographically-least shortest word."""
    full = (1 << n) - 1
    if target_mask == full:
        return ""
    parent = {full: None}        # mask -> (prev_mask, letter)
    dq = deque([full])
    while dq:
        s = dq.popleft()
        for letter, t in (("a", a), ("b", b)):
            ns = _image(s, t)
            if ns not in parent:
                parent[ns] = (s, letter)
                if ns == target_mask:
                    word = []
                    cur = ns
                    while parent[cur] is not None:
                        prev, lt = parent[cur]
                        word.append(lt)
                        cur = prev
                    return "".join(reversed(word))
                dq.append(ns)
    return None


def certify_don_bound(a, b, n: int) -> dict:
    """Produce, for every reachable subset, a shortest reaching word and check it.
    Returns the worst case against Don's bound n(n−k), with an explicit certificate
    (the word) — all re-verified by `reaches`."""
    full = (1 << n) - 1
    # reachable masks via the same independent BFS
    reachable = {full}
    dq = deque([full])
    while dq:
        s = dq.popleft()
        for t in (a, b):
            ns = _image(s, t)
            if ns not in reachable:
                reachable.add(ns)
                dq.append(ns)

    worst = {"ratio": 0.0}
    don_violation = None
    for s in reachable:
        k = bin(s).count("1")
        if k == n:
            continue
        w = shortest_reaching_word(a, b, n, s)
        assert w is not None and reaches(a, b, n, w, s)      # certificate holds
        d, bound = len(w), n * (n - k)
        if d > bound and don_violation is None:
            don_violation = {"subset": s, "k": k, "len": d, "word": w, "bound": bound}
        if bound > 0 and d / bound > worst["ratio"]:
            worst = {"ratio": round(d / bound, 4), "subset_states": [i for i in range(n) if s >> i & 1],
                     "k": k, "shortest_len": d, "don_bound": bound, "word": w}
    return {"don_bound_holds": don_violation is None,
            "don_violation": don_violation,
            "worst_case_vs_don_bound": worst}


__all__ = ["apply_word", "reaches", "shortest_reaching_word", "certify_don_bound"]
