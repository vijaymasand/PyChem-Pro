"""
MMFF94 out-of-plane bending (Wilson angle in degrees).

For a 3-coordinated center j with neighbors {i, k, l}:
    Wilson angle chi = angle between bond j-i and the plane of j, k, l.

Energy per (i, j, k, l) entry:
    E_oop = (0.043844/2) * k_oop * chi_deg^2
          = (0.043844/2 * (180/pi)^2) * k_oop * chi_rad^2

The ArraysBuilder creates 3 entries per OOP center (one per outer atom).

Reference: Jmol Util.restorativeForceAndOutOfPlaneAngleRadians.

Note on gradient: the analytical Wilson-angle derivative (Jmol's
Util.java:300-415) is long; for Phase 1 we use a per-row finite-difference
that's correct even when multiple rows share atoms (which happens at every
3-coordinated center — three OOP entries share the center). Phase 2 will
replace with the analytical form for ~5x speedup. The per-row FD path is
O(12*N_rows) scalar Wilson-angle evaluations.
"""
from __future__ import annotations
import numpy as np

from src.services.forcefield._engine.arrays import InteractionArrays
from src.services.forcefield._engine._math import scatter_add

_EPS = 1e-12
_RAD_TO_DEG = 180.0 / np.pi
_FOOP = 0.043844 * 0.5 * _RAD_TO_DEG * _RAD_TO_DEG   # = 0.043844/2 * (180/pi)^2
_FOOP_GRAD = 0.043844 * _RAD_TO_DEG * _RAD_TO_DEG    # dE/dchi_rad = FOOP_GRAD * k_oop * chi_rad


class OopCalc:

    @staticmethod
    def energy(coords: np.ndarray, arr: InteractionArrays) -> float:
        if arr.oop_center.size == 0:
            return 0.0
        chi_rad = _wilson_angles(coords, arr.oop_i, arr.oop_center,
                                 arr.oop_j, arr.oop_k)
        return float((_FOOP * arr.oop_koop * chi_rad * chi_rad).sum())

    @staticmethod
    def gradient(coords: np.ndarray, arr: InteractionArrays) -> np.ndarray:
        return OopCalc.energy_and_gradient(coords, arr)[1]

    @staticmethod
    def energy_and_gradient(coords: np.ndarray, arr: InteractionArrays):
        n = arr.n_atoms
        if arr.oop_center.size == 0:
            return 0.0, np.zeros((n, 3), dtype=np.float64)
        chi_rad, dchi_di, dchi_dj, dchi_dk, dchi_dl = _wilson_angles_and_grad(
            coords, arr.oop_i, arr.oop_center, arr.oop_j, arr.oop_k)
        e = float((_FOOP * arr.oop_koop * chi_rad * chi_rad).sum())
        dE_dchi = _FOOP_GRAD * arr.oop_koop * chi_rad

        # i <-> arr.oop_i (outer), j <-> arr.oop_center, k <-> arr.oop_j (outer), l <-> arr.oop_k (outer)
        force_i = dE_dchi[:, None] * dchi_di
        force_j = dE_dchi[:, None] * dchi_dj
        force_k = dE_dchi[:, None] * dchi_dk
        force_l = dE_dchi[:, None] * dchi_dl

        g = (scatter_add(force_i, arr.oop_i, n)
             + scatter_add(force_j, arr.oop_center, n)
             + scatter_add(force_k, arr.oop_j, n)
             + scatter_add(force_l, arr.oop_k, n))
        return e, g


# ────────── Wilson angle helpers (vectorized) ──────────

def _wilson_angles(coords, idx_i, idx_j, idx_k, idx_l):
    """Wilson angle of bond j-i relative to plane of j, k, l (radians).

    idx_j is the central atom; idx_i, idx_k, idx_l are the three neighbors.
    chi is the angle between bond j-i and the plane spanned by j-k and j-l.
    """
    # Vectors from center j
    vi = coords[idx_i] - coords[idx_j]
    vk = coords[idx_k] - coords[idx_j]
    vl = coords[idx_l] - coords[idx_j]

    # Normal to plane (k, j, l)
    normal = np.cross(vk, vl)
    n_norm = np.linalg.norm(normal, axis=1) + _EPS
    n_hat = normal / n_norm[:, None]

    vi_norm = np.linalg.norm(vi, axis=1) + _EPS
    vi_hat = vi / vi_norm[:, None]

    # sin(chi) = dot(vi_hat, n_hat); chi is positive when i is out of the plane
    sin_chi = np.einsum("ij,ij->i", vi_hat, n_hat)
    sin_chi = np.clip(sin_chi, -1.0, 1.0)
    return np.arcsin(sin_chi)


def _wilson_from_coords(ci, cj, ck, cl):
    """Wilson angle (radians) from per-row coordinate arrays, each shape (N, 3).

    Identical maths to :func:`_wilson_angles` but takes the four atom positions
    directly rather than via index arrays, so callers can perturb a per-row copy
    without disturbing the shared global coordinates.
    """
    vi = ci - cj
    vk = ck - cj
    vl = cl - cj
    normal = np.cross(vk, vl)
    n_hat = normal / (np.linalg.norm(normal, axis=1) + _EPS)[:, None]
    vi_hat = vi / (np.linalg.norm(vi, axis=1) + _EPS)[:, None]
    sin_chi = np.clip(np.einsum("ij,ij->i", vi_hat, n_hat), -1.0, 1.0)
    return np.arcsin(sin_chi)


def _wilson_angles_and_grad(coords, idx_i, idx_j, idx_k, idx_l):
    """Return (chi_rad, dchi/dx_i, dchi/dx_j, dchi/dx_k, dchi/dx_l).

    Vectorized central finite-difference: for each of the 4 atom slots and 3
    Cartesian directions, perturb that slot's per-row coordinates for ALL rows
    at once and recompute chi. This is numerically identical to a per-row loop
    (each row owns its own slot in the (N,3) arrays, so shared atoms don't leak
    between rows) but replaces 24*N scalar Wilson-angle evaluations with 24
    vectorized ones — the dominant per-iteration cost for molecules rich in
    sp2/aromatic centers.
    """
    # Per-row copies (independent even when an atom appears in several rows).
    ci = coords[idx_i].astype(np.float64, copy=True)
    cj = coords[idx_j].astype(np.float64, copy=True)
    ck = coords[idx_k].astype(np.float64, copy=True)
    cl = coords[idx_l].astype(np.float64, copy=True)

    chi = _wilson_from_coords(ci, cj, ck, cl)
    eps = 1e-6
    inv_2eps = 1.0 / (2.0 * eps)

    grads = []
    for base in (ci, cj, ck, cl):
        g = np.zeros_like(base)
        for d in range(3):
            orig = base[:, d].copy()
            base[:, d] = orig + eps
            chi_p = _wilson_from_coords(ci, cj, ck, cl)
            base[:, d] = orig - eps
            chi_m = _wilson_from_coords(ci, cj, ck, cl)
            base[:, d] = orig
            g[:, d] = (chi_p - chi_m) * inv_2eps
        grads.append(g)

    return chi, grads[0], grads[1], grads[2], grads[3]
