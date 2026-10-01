"""
Shared NumPy primitives used by every per-term calculator.

- scatter_add: fast atomic-index gradient accumulation
- batch_distance / batch_angle / batch_dihedral: vectorized geometry

Designed to be the *only* place that knows how to scatter forces back
onto atoms. Calculators pass (force_vec, atom_indices) to scatter_add
and never touch np.add.at directly.
"""
from __future__ import annotations
import numpy as np

_EPS = 1e-12


def scatter_add(values: np.ndarray, indices: np.ndarray, n_atoms: int) -> np.ndarray:
    """Accumulate per-interaction force vectors onto per-atom gradient.

    Equivalent to::
        out = np.zeros((n_atoms, 3))
        np.add.at(out, indices, values)

    but uses np.bincount per axis which is 3-5x faster on N>500 pairs
    and avoids the np.add.at GIL-holding C loop.

    Args:
        values: (Npairs, 3) float64. Force vectors per interaction.
        indices: (Npairs,) int32. Target atom index per row of values.
        n_atoms: total number of atoms (size of output's first axis).

    Returns:
        (n_atoms, 3) float64 - accumulated gradient contribution.
    """
    out = np.empty((n_atoms, 3), dtype=np.float64)
    out[:, 0] = np.bincount(indices, weights=values[:, 0], minlength=n_atoms)
    out[:, 1] = np.bincount(indices, weights=values[:, 1], minlength=n_atoms)
    out[:, 2] = np.bincount(indices, weights=values[:, 2], minlength=n_atoms)
    return out


def batch_distance(coords: np.ndarray, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """Vectorized ||coords[i] - coords[j]||_2 for paired index arrays."""
    diff = coords[i] - coords[j]
    return np.linalg.norm(diff, axis=1)


def batch_angle(coords: np.ndarray,
                i: np.ndarray, j: np.ndarray, k: np.ndarray) -> np.ndarray:
    """Vectorized angle i-j-k in radians. j is the vertex."""
    v1 = coords[i] - coords[j]
    v2 = coords[k] - coords[j]
    n1 = np.linalg.norm(v1, axis=1) + _EPS
    n2 = np.linalg.norm(v2, axis=1) + _EPS
    cos_t = np.einsum("ij,ij->i", v1, v2) / (n1 * n2)
    cos_t = np.clip(cos_t, -1.0, 1.0)
    return np.arccos(cos_t)


def batch_dihedral(coords: np.ndarray,
                   i: np.ndarray, j: np.ndarray,
                   k: np.ndarray, l: np.ndarray) -> np.ndarray:
    """Vectorized dihedral i-j-k-l in radians, atan2-based, range (-pi, pi]."""
    b1 = coords[j] - coords[i]
    b2 = coords[k] - coords[j]
    b3 = coords[l] - coords[k]

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    b2_norm = np.linalg.norm(b2, axis=1, keepdims=True) + _EPS
    b2_hat = b2 / b2_norm

    m1 = np.cross(n1, b2_hat)

    x = np.einsum("ij,ij->i", n1, n2)
    y = np.einsum("ij,ij->i", m1, n2)
    return np.arctan2(y, x)


def angle_jacobian(coords: np.ndarray,
                   i: np.ndarray, j: np.ndarray, k: np.ndarray):
    """Vectorized angle gradient: returns (g_i, g_j, g_k), each (N,3).

    Each g_x is dtheta/dr_x at angle (i, j, k), where theta is in radians.
    """
    v1 = coords[i] - coords[j]
    v2 = coords[k] - coords[j]
    r1 = np.linalg.norm(v1, axis=1) + _EPS
    r2 = np.linalg.norm(v2, axis=1) + _EPS
    v1n = v1 / r1[:, None]
    v2n = v2 / r2[:, None]
    cos_t = np.clip(np.einsum("ij,ij->i", v1n, v2n), -1.0, 1.0)
    sin_t = np.sqrt(np.maximum(1.0 - cos_t * cos_t, _EPS))
    inv_sin = 1.0 / sin_t

    g_i = (cos_t[:, None] * v1n - v2n) * (inv_sin[:, None] / r1[:, None])
    g_k = (cos_t[:, None] * v2n - v1n) * (inv_sin[:, None] / r2[:, None])
    g_j = -(g_i + g_k)
    return g_i, g_j, g_k


def dihedral_jacobian(coords: np.ndarray,
                      i: np.ndarray, j: np.ndarray,
                      k: np.ndarray, l: np.ndarray):
    """Vectorized dihedral gradient (Bekker-Berendsen-van Gunsteren 1995).

    Returns (g_i, g_j, g_k, g_l), each (N, 3). Units: rad/A.
    """
    b1 = coords[j] - coords[i]
    b2 = coords[k] - coords[j]
    b3 = coords[l] - coords[k]

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1_sq = np.einsum("ij,ij->i", n1, n1) + _EPS
    n2_sq = np.einsum("ij,ij->i", n2, n2) + _EPS
    b2_norm = np.linalg.norm(b2, axis=1) + _EPS
    b2_sq = b2_norm * b2_norm

    g_i = (b2_norm[:, None] / n1_sq[:, None]) * n1
    g_l = -(b2_norm[:, None] / n2_sq[:, None]) * n2

    f1 = np.einsum("ij,ij->i", b1, b2) / b2_sq
    f2 = np.einsum("ij,ij->i", b3, b2) / b2_sq

    g_j = -g_i * (1.0 + f1[:, None]) + g_l * f2[:, None]
    g_k = -(g_i + g_j + g_l)
    return g_i, g_j, g_k, g_l
