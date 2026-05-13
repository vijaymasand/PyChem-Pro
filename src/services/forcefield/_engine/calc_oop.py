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


def _wilson_angle_scalar(ci, cj, ck, cl):
    """Scalar version of Wilson angle for one (i, j, k, l) tuple of coords (3,)."""
    vi = ci - cj
    vk = ck - cj
    vl = cl - cj
    normal = np.cross(vk, vl)
    n_norm = np.linalg.norm(normal) + _EPS
    n_hat = normal / n_norm
    vi_norm = np.linalg.norm(vi) + _EPS
    vi_hat = vi / vi_norm
    sin_chi = float(np.dot(vi_hat, n_hat))
    if sin_chi > 1.0:
        sin_chi = 1.0
    elif sin_chi < -1.0:
        sin_chi = -1.0
    return float(np.arcsin(sin_chi))


def _wilson_angles_and_grad(coords, idx_i, idx_j, idx_k, idx_l):
    """Return (chi_rad, dchi/dx_i, dchi/dx_j, dchi/dx_k, dchi/dx_l).

    Per-row finite-difference: for each row and each of the 4 atom positions
    and each Cartesian dimension, perturb ONLY that atom in ONLY this row's
    context, recompute chi for this row, take central difference.

    This is correct even when multiple rows share atoms.
    """
    chi = _wilson_angles(coords, idx_i, idx_j, idx_k, idx_l)
    N = chi.shape[0]
    eps = 1e-6
    inv_2eps = 1.0 / (2.0 * eps)

    g_i = np.zeros((N, 3), dtype=np.float64)
    g_j = np.zeros((N, 3), dtype=np.float64)
    g_k = np.zeros((N, 3), dtype=np.float64)
    g_l = np.zeros((N, 3), dtype=np.float64)

    for row in range(N):
        ai = idx_i[row]; aj = idx_j[row]; ak = idx_k[row]; al = idx_l[row]
        ci = coords[ai].copy()
        cj = coords[aj].copy()
        ck = coords[ak].copy()
        cl = coords[al].copy()
        for d in range(3):
            # dchi/dx_i
            ci_p = ci.copy(); ci_p[d] += eps
            ci_m = ci.copy(); ci_m[d] -= eps
            g_i[row, d] = (_wilson_angle_scalar(ci_p, cj, ck, cl)
                           - _wilson_angle_scalar(ci_m, cj, ck, cl)) * inv_2eps
            # dchi/dx_j (center)
            cj_p = cj.copy(); cj_p[d] += eps
            cj_m = cj.copy(); cj_m[d] -= eps
            g_j[row, d] = (_wilson_angle_scalar(ci, cj_p, ck, cl)
                           - _wilson_angle_scalar(ci, cj_m, ck, cl)) * inv_2eps
            # dchi/dx_k
            ck_p = ck.copy(); ck_p[d] += eps
            ck_m = ck.copy(); ck_m[d] -= eps
            g_k[row, d] = (_wilson_angle_scalar(ci, cj, ck_p, cl)
                           - _wilson_angle_scalar(ci, cj, ck_m, cl)) * inv_2eps
            # dchi/dx_l
            cl_p = cl.copy(); cl_p[d] += eps
            cl_m = cl.copy(); cl_m[d] -= eps
            g_l[row, d] = (_wilson_angle_scalar(ci, cj, ck, cl_p)
                           - _wilson_angle_scalar(ci, cj, ck, cl_m)) * inv_2eps

    return chi, g_i, g_j, g_k, g_l
