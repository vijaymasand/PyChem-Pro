"""
Iterative superposition with outlier rejection (pure NumPy).

Wraps :func:`kabsch` with repeated rejection of poorly-fitting atom pairs,
giving a robust fit when only *part* of two structures corresponds — a flexible
loop, an induced-fit side chain, or a ligand plus a handful of mismatched
atoms.  This mirrors PyMOL's ``align`` refinement cycles.

Ported from patinae-algos ``align/superpose.rs``.  Each cycle rejects any pair
whose ``distance / current_RMSD`` exceeds ``cutoff``, then re-fits on the
survivors, for up to ``cycles`` rounds (or until nothing more is rejected).
"""
from __future__ import annotations

import numpy as np

from src.services.alignment._engine.kabsch import kabsch, KabschResult


class SuperposeResult:
    """Outcome of an iterative superposition.

    Attributes:
        transform:         Final :class:`KabschResult` (rotation + translation).
        initial_rmsd:      RMSD before any outlier rejection.
        final_rmsd:        RMSD after rejection.
        cycles_performed:  Number of rejection cycles actually run.
        n_rejected:        Pairs discarded as outliers.
        n_aligned:         Pairs used in the final fit.
        pairs:             The surviving ``(src_idx, tgt_idx)`` pairs.
    """

    __slots__ = ("transform", "initial_rmsd", "final_rmsd",
                 "cycles_performed", "n_rejected", "n_aligned", "pairs")

    def __init__(self, transform, initial_rmsd, final_rmsd,
                 cycles_performed, n_rejected, n_aligned, pairs):
        self.transform = transform
        self.initial_rmsd = initial_rmsd
        self.final_rmsd = final_rmsd
        self.cycles_performed = cycles_performed
        self.n_rejected = n_rejected
        self.n_aligned = n_aligned
        self.pairs = pairs


def superpose(source_coords, target_coords, pairs, cycles=5, cutoff=2.0) -> SuperposeResult:
    """Iteratively superpose *source* onto *target* over *pairs*.

    Args:
        source_coords: ``(Ns, 3)`` array; ``pairs[k][0]`` indexes into it.
        target_coords: ``(Nt, 3)`` array; ``pairs[k][1]`` indexes into it.
        pairs:         List of ``(source_index, target_index)`` correspondences.
        cycles:        Max outlier-rejection cycles (0 → plain single Kabsch).
        cutoff:        Reject a pair when ``distance / RMSD > cutoff``.

    Returns:
        :class:`SuperposeResult`.

    Raises:
        ValueError: If fewer than 3 pairs are supplied or all are rejected.
    """
    src_all = np.asarray(source_coords, dtype=np.float64)
    tgt_all = np.asarray(target_coords, dtype=np.float64)

    active = list(pairs)
    if len(active) < 3:
        raise ValueError(f"superposition needs at least 3 pairs, got {len(active)}")

    def extract(prs):
        si = [p[0] for p in prs]
        ti = [p[1] for p in prs]
        return src_all[si], tgt_all[ti]

    src, tgt = extract(active)
    result = kabsch(src, tgt)
    initial_rmsd = result.rmsd
    cycles_performed = 0

    for cycle in range(cycles):
        if len(active) < 3:
            break

        src, tgt = extract(active)
        transformed = result.apply(src)
        rms = result.rmsd

        if rms > 1e-6:
            dist = np.linalg.norm(transformed - tgt, axis=1)
            keep = np.nonzero(dist / rms <= cutoff)[0]
        else:
            # Perfect fit — nothing to reject.
            keep = np.arange(len(active))

        # Converged: no pair rejected this cycle.
        if len(keep) == len(active):
            cycles_performed = cycle + 1
            break

        # Too few survivors to re-fit.
        if len(keep) < 3:
            if len(keep) == 0:
                raise ValueError("all pairs rejected as outliers")
            cycles_performed = cycle + 1
            break

        active = [active[i] for i in keep]
        src, tgt = extract(active)
        result = kabsch(src, tgt)
        cycles_performed = cycle + 1

    n_aligned = len(active)
    n_rejected = len(pairs) - n_aligned
    return SuperposeResult(
        transform=result,
        initial_rmsd=initial_rmsd,
        final_rmsd=result.rmsd,
        cycles_performed=cycles_performed,
        n_rejected=n_rejected,
        n_aligned=n_aligned,
        pairs=active,
    )
