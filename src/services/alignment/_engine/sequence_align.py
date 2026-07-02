"""
Needleman-Wunsch global sequence alignment with affine gap penalties.

Used to match residues between two protein chains *before* structural
superposition: align the one-letter Cα sequences, then feed the matched
residue pairs to :func:`~src.services.alignment._engine.superpose.superpose`.

Ported from patinae-algos ``align/sequence_align.rs`` — the same three-matrix
(M / X / Y) affine-gap DP with traceback.  The Rust version defaults to a
BLOSUM62 substitution matrix; to stay dependency-free and compact this port
uses simple match/mismatch scoring by default (which is exact for identical or
near-identical chains — the common "same protein, different conformation" case)
and accepts an optional ``substitution`` callable for homolog alignment.
"""
from __future__ import annotations

import numpy as np

# Traceback source matrices.
_M, _X, _Y = 0, 1, 2
_NEG_INF = -1.0e18


class SequenceAlignment:
    """Result of a global sequence alignment.

    Attributes:
        matched_pairs: List of ``(source_pos, target_pos)`` for aligned (non-gap)
                       columns — indices into the input sequences.
        score:         Optimal alignment score.
        n_matched:     Number of aligned columns.
        identity:      Fraction of aligned columns with identical symbols.
    """

    __slots__ = ("matched_pairs", "score", "n_matched", "identity")

    def __init__(self, matched_pairs, score, n_matched, identity):
        self.matched_pairs = matched_pairs
        self.score = score
        self.n_matched = n_matched
        self.identity = identity


def global_align(source, target, match=1.0, mismatch=-1.0,
                 gap_open=-10.0, gap_extend=-0.5, substitution=None) -> SequenceAlignment:
    """Global (Needleman-Wunsch) alignment of two symbol sequences.

    Args:
        source, target: Sequences of single-character symbols (e.g. one-letter
                        residue codes), as strings or lists of chars.
        match:          Score for identical symbols (used if *substitution* is None).
        mismatch:       Score for differing symbols (used if *substitution* is None).
        gap_open:       Affine gap-opening penalty (negative).
        gap_extend:     Affine gap-extension penalty (negative).
        substitution:   Optional ``f(a, b) -> float`` overriding match/mismatch.

    Returns:
        :class:`SequenceAlignment`.
    """
    a = list(source)
    b = list(target)
    m, n = len(a), len(b)
    if m == 0 or n == 0:
        return SequenceAlignment([], 0.0, 0, 0.0)

    if substitution is None:
        def substitution(x, y, _m=match, _mm=mismatch):
            return _m if x == y else _mm

    rows, cols = m + 1, n + 1
    dp_m = np.full((rows, cols), _NEG_INF)
    dp_x = np.full((rows, cols), _NEG_INF)   # gap in target (consume source)
    dp_y = np.full((rows, cols), _NEG_INF)   # gap in source (consume target)
    tb_m = np.zeros((rows, cols), dtype=np.int8)
    tb_x = np.zeros((rows, cols), dtype=np.int8)
    tb_y = np.zeros((rows, cols), dtype=np.int8)

    dp_m[0, 0] = 0.0
    for i in range(1, rows):
        dp_x[i, 0] = gap_open + i * gap_extend
        tb_x[i, 0] = _X
    for j in range(1, cols):
        dp_y[0, j] = gap_open + j * gap_extend
        tb_y[0, j] = _Y

    for i in range(1, rows):
        ai = a[i - 1]
        for j in range(1, cols):
            sub = substitution(ai, b[j - 1])

            # M[i, j]: match/mismatch — best of the three diagonal predecessors.
            fm = dp_m[i - 1, j - 1]
            fx = dp_x[i - 1, j - 1]
            fy = dp_y[i - 1, j - 1]
            if fm >= fx and fm >= fy:
                dp_m[i, j] = fm + sub
                tb_m[i, j] = _M
            elif fx >= fy:
                dp_m[i, j] = fx + sub
                tb_m[i, j] = _X
            else:
                dp_m[i, j] = fy + sub
                tb_m[i, j] = _Y

            # X[i, j]: gap in target (consume source i) — from (i-1, j).
            open_x = dp_m[i - 1, j] + gap_open + gap_extend
            ext_x = dp_x[i - 1, j] + gap_extend
            if open_x >= ext_x:
                dp_x[i, j] = open_x
                tb_x[i, j] = _M
            else:
                dp_x[i, j] = ext_x
                tb_x[i, j] = _X

            # Y[i, j]: gap in source (consume target j) — from (i, j-1).
            open_y = dp_m[i, j - 1] + gap_open + gap_extend
            ext_y = dp_y[i, j - 1] + gap_extend
            if open_y >= ext_y:
                dp_y[i, j] = open_y
                tb_y[i, j] = _M
            else:
                dp_y[i, j] = ext_y
                tb_y[i, j] = _Y

    # Best terminal cell / matrix.
    sm, sx, sy = dp_m[m, n], dp_x[m, n], dp_y[m, n]
    if sm >= sx and sm >= sy:
        score, cur = sm, _M
    elif sx >= sy:
        score, cur = sx, _X
    else:
        score, cur = sy, _Y

    # Traceback.
    pairs = []
    i, j = m, n
    while i > 0 or j > 0:
        if cur == _M:
            if i == 0 or j == 0:
                break
            nxt = tb_m[i, j]
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
            cur = nxt
        elif cur == _X:
            if i == 0:
                break
            nxt = tb_x[i, j]
            i -= 1
            cur = nxt if nxt in (_M, _X) else _M
        else:  # _Y
            if j == 0:
                break
            nxt = tb_y[i, j]
            j -= 1
            cur = nxt if nxt in (_M, _Y) else _M

    pairs.reverse()

    n_matched = len(pairs)
    n_identical = sum(1 for si, ti in pairs if a[si] == b[ti])
    identity = (n_identical / n_matched) if n_matched else 0.0
    return SequenceAlignment(pairs, float(score), n_matched, identity)
