"""
MMFF94 stretch-bend coupling (new term — entirely absent in legacy code).

E_sb = 2.51210 * (kbai*(r_ij - r0_ij) + kbak*(r_jk - r0_jk)) * (theta_deg - theta0_deg)

For each angle (i, j, k):
    - r_ij is the i-j bond distance
    - r_jk is the j-k bond distance
    - theta is the i-j-k angle
"""
from __future__ import annotations
import numpy as np

from src.services.forcefield._engine.arrays import InteractionArrays
from src.services.forcefield._engine._math import (
    scatter_add, batch_distance, batch_angle, angle_jacobian,
)

_K_SB = 2.51210
_RAD_TO_DEG = 180.0 / np.pi


class StretchBendCalc:

    @staticmethod
    def energy(coords: np.ndarray, arr: InteractionArrays) -> float:
        if arr.sb_i.size == 0:
            return 0.0
        r_ij = batch_distance(coords, arr.sb_i, arr.sb_j)
        r_jk = batch_distance(coords, arr.sb_j, arr.sb_k)
        theta_rad = batch_angle(coords, arr.sb_i, arr.sb_j, arr.sb_k)
        dt_deg = theta_rad * _RAD_TO_DEG - arr.sb_theta0_deg
        d_rij = r_ij - arr.sb_r0_ij
        d_rjk = r_jk - arr.sb_r0_jk
        e = _K_SB * (arr.sb_kbai * d_rij + arr.sb_kbak * d_rjk) * dt_deg
        return float(e.sum())

    @staticmethod
    def gradient(coords: np.ndarray, arr: InteractionArrays) -> np.ndarray:
        return StretchBendCalc.energy_and_gradient(coords, arr)[1]

    @staticmethod
    def energy_and_gradient(coords: np.ndarray, arr: InteractionArrays):
        n = arr.n_atoms
        if arr.sb_i.size == 0:
            return 0.0, np.zeros((n, 3), dtype=np.float64)

        diff_ij = coords[arr.sb_i] - coords[arr.sb_j]
        diff_jk = coords[arr.sb_k] - coords[arr.sb_j]
        r_ij = np.linalg.norm(diff_ij, axis=1)
        r_jk = np.linalg.norm(diff_jk, axis=1)
        theta_rad = batch_angle(coords, arr.sb_i, arr.sb_j, arr.sb_k)
        theta_deg = theta_rad * _RAD_TO_DEG
        dt_deg = theta_deg - arr.sb_theta0_deg
        d_rij = r_ij - arr.sb_r0_ij
        d_rjk = r_jk - arr.sb_r0_jk

        sum_kr = arr.sb_kbai * d_rij + arr.sb_kbak * d_rjk
        e = float((_K_SB * sum_kr * dt_deg).sum())

        # Partials:
        #   dE/dr_ij = 2.51210 * kbai * dt_deg
        #   dE/dr_jk = 2.51210 * kbak * dt_deg
        #   dE/dtheta_rad = 2.51210 * (kbai*d_rij + kbak*d_rjk) * RAD_TO_DEG
        dE_drij = _K_SB * arr.sb_kbai * dt_deg
        dE_drjk = _K_SB * arr.sb_kbak * dt_deg
        dE_dtheta = _K_SB * sum_kr * _RAD_TO_DEG

        # Bond-ij chain rule: dr_ij/dx_i = (x_i - x_j)/r_ij, dr_ij/dx_j = -that.
        r_ij_safe = np.maximum(r_ij, 1e-12)
        r_jk_safe = np.maximum(r_jk, 1e-12)
        unit_ij = diff_ij / r_ij_safe[:, None]   # (i - j) / |i-j|
        unit_jk = diff_jk / r_jk_safe[:, None]   # (k - j) / |k-j|

        f_i_from_rij = dE_drij[:, None] * unit_ij
        f_j_from_rij = -f_i_from_rij
        f_k_from_rjk = dE_drjk[:, None] * unit_jk
        f_j_from_rjk = -f_k_from_rjk

        # Angle chain rule
        g_i, g_j, g_k = angle_jacobian(coords, arr.sb_i, arr.sb_j, arr.sb_k)
        f_i_from_t = dE_dtheta[:, None] * g_i
        f_j_from_t = dE_dtheta[:, None] * g_j
        f_k_from_t = dE_dtheta[:, None] * g_k

        force_i = f_i_from_rij + f_i_from_t
        force_j = f_j_from_rij + f_j_from_rjk + f_j_from_t
        force_k = f_k_from_rjk + f_k_from_t

        g = (scatter_add(force_i, arr.sb_i, n)
             + scatter_add(force_j, arr.sb_j, n)
             + scatter_add(force_k, arr.sb_k, n))
        return e, g
