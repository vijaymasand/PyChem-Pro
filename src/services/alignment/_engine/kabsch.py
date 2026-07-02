"""
Kabsch algorithm for optimal rigid-body superposition (pure NumPy).

Given two sets of *corresponding* 3-D points, finds the rotation ``R`` and
translation ``t`` that minimise the (optionally weighted) RMSD between
``R · source + t`` and ``target``.

Ported from patinae-algos ``align/kabsch.rs`` (Rust) — identical maths, using
:func:`numpy.linalg.svd` in place of the hand-rolled 3×3 SVD.  A reflection is
removed via the sign of ``det(V · Uᵀ)`` so ``R`` is always a proper rotation
(``det = +1``), never a mirror image.

Reference: Kabsch W. (1976) *Acta Crystallogr.* A32:922-923.
"""
from __future__ import annotations

import numpy as np


class KabschResult:
    """Optimal transform mapping ``source`` onto ``target``.

    Attributes:
        rotation:    (3, 3) proper-rotation matrix, applied before translation.
        translation: (3,) vector, applied after rotation.
        rmsd:        RMSD of the fit (Å), over the point pairs supplied.
        n_atoms:     Number of point pairs used.
    """

    __slots__ = ("rotation", "translation", "rmsd", "n_atoms")

    def __init__(self, rotation, translation, rmsd, n_atoms):
        self.rotation = rotation
        self.translation = translation
        self.rmsd = rmsd
        self.n_atoms = n_atoms

    def apply(self, coords):
        """Return a transformed copy of an ``(N, 3)`` coordinate array."""
        c = np.asarray(coords, dtype=np.float64)
        return c @ self.rotation.T + self.translation

    def __repr__(self):
        return (f"KabschResult(rmsd={self.rmsd:.4f}, n_atoms={self.n_atoms})")


def kabsch(source, target, weights=None) -> KabschResult:
    """Compute the optimal superposition of *source* onto *target*.

    Args:
        source:  ``(N, 3)`` coordinates to be moved.
        target:  ``(N, 3)`` reference coordinates (same length as *source*).
        weights: Optional ``(N,)`` non-negative per-pair weights (e.g. atomic
                 masses or 1/B-factor).  ``None`` weights every pair equally.

    Returns:
        :class:`KabschResult` with ``rotation``, ``translation`` and ``rmsd``.

    Raises:
        ValueError: If shapes disagree, are not ``(N, 3)``, or ``N < 3``.
    """
    src = np.asarray(source, dtype=np.float64)
    tgt = np.asarray(target, dtype=np.float64)
    if src.ndim != 2 or src.shape[1] != 3 or src.shape != tgt.shape:
        raise ValueError(
            f"source/target must be equal (N, 3) arrays, got {src.shape} and {tgt.shape}")
    n = src.shape[0]
    if n < 3:
        raise ValueError(f"Kabsch needs at least 3 point pairs, got {n}")

    # 1. (Weighted) centroids
    if weights is not None:
        w = np.asarray(weights, dtype=np.float64).reshape(-1)
        if w.shape[0] != n:
            raise ValueError("weights length must match number of points")
        w = np.clip(w, 0.0, None)
        total_w = float(w.sum())
        if total_w <= 0.0:
            raise ValueError("weights sum to zero")
        centroid_src = (src * w[:, None]).sum(axis=0) / total_w
        centroid_tgt = (tgt * w[:, None]).sum(axis=0) / total_w
    else:
        w = None
        centroid_src = src.mean(axis=0)
        centroid_tgt = tgt.mean(axis=0)

    # 2. Centre both point sets
    p = src - centroid_src
    q = tgt - centroid_tgt

    # 3. Cross-covariance matrix H = Pᵀ · (W) · Q
    h = (p * w[:, None]).T @ q if w is not None else p.T @ q

    # 4. SVD of H
    u, _s, vt = np.linalg.svd(h)
    v = vt.T

    # 5. R = V · diag(1, 1, d) · Uᵀ, where d removes any reflection
    d = 1.0 if np.linalg.det(v @ u.T) >= 0.0 else -1.0
    r = (v * np.array([1.0, 1.0, d])) @ u.T

    # 6. Translation maps the source centroid onto the target centroid
    t = centroid_tgt - r @ centroid_src

    # 7. RMSD after superposition
    rotated = p @ r.T
    diff_sq = np.sum((rotated - q) ** 2, axis=1)
    if w is not None:
        rmsd_val = float(np.sqrt(np.sum(w * diff_sq) / total_w))
    else:
        rmsd_val = float(np.sqrt(diff_sq.mean()))

    return KabschResult(r, t, rmsd_val, n)


def rmsd(coords_a, coords_b) -> float:
    """RMSD between two equal-length coordinate sets (no superposition)."""
    a = np.asarray(coords_a, dtype=np.float64)
    b = np.asarray(coords_b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"coordinate sets must have equal shape, got {a.shape} and {b.shape}")
    if a.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def apply_transform(coords, result: KabschResult):
    """Apply a :class:`KabschResult` to ``(N, 3)`` coordinates (returns a copy)."""
    return result.apply(coords)
