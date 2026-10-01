"""
MMFF94 torsion (Halgren 1996) — corrected (1+cos) signs.

E_t = 0.5 * [V1*(1 + cos phi) + V2*(1 - cos 2phi) + V3*(1 + cos 3phi)]
With c = cos phi:
    E_t = 0.5 * [V1(1+c) + V2(2-2c^2) + V3(1 + c(4c^2-3))]
dE_t/dphi = 0.5 * sin phi * [-V1 + 4*V2*c + 3*V3*(1 - 4c^2)]

Legacy torsion.py used (1 - cos phi) for V1 and V3 — wrong sign.
"""
from __future__ import annotations
import numpy as np

from src.services.forcefield._engine.arrays import InteractionArrays
from src.services.forcefield._engine._math import (
    scatter_add, batch_dihedral, dihedral_jacobian,
)


class TorsionCalc:

    @staticmethod
    def energy(coords: np.ndarray, arr: InteractionArrays) -> float:
        if arr.tor_i.size == 0:
            return 0.0
        phi = batch_dihedral(coords, arr.tor_i, arr.tor_j, arr.tor_k, arr.tor_l)
        c = np.cos(phi)
        e = 0.5 * (arr.tor_v1 * (1.0 + c)
                   + arr.tor_v2 * (2.0 - 2.0 * c * c)
                   + arr.tor_v3 * (1.0 + c * (4.0 * c * c - 3.0)))
        return float(e.sum())

    @staticmethod
    def gradient(coords: np.ndarray, arr: InteractionArrays) -> np.ndarray:
        return TorsionCalc.energy_and_gradient(coords, arr)[1]

    @staticmethod
    def energy_and_gradient(coords: np.ndarray, arr: InteractionArrays):
        n = arr.n_atoms
        if arr.tor_i.size == 0:
            return 0.0, np.zeros((n, 3), dtype=np.float64)
        phi = batch_dihedral(coords, arr.tor_i, arr.tor_j, arr.tor_k, arr.tor_l)
        c = np.cos(phi); s = np.sin(phi)
        c2 = c * c
        e = float((0.5 * (arr.tor_v1 * (1.0 + c)
                          + arr.tor_v2 * (2.0 - 2.0 * c2)
                          + arr.tor_v3 * (1.0 + c * (4.0 * c2 - 3.0)))).sum())
        dEdphi = 0.5 * s * (-arr.tor_v1
                            + 4.0 * arr.tor_v2 * c
                            + 3.0 * arr.tor_v3 * (1.0 - 4.0 * c2))
        g_i, g_j, g_k, g_l = dihedral_jacobian(
            coords, arr.tor_i, arr.tor_j, arr.tor_k, arr.tor_l)
        force_i = dEdphi[:, None] * g_i
        force_j = dEdphi[:, None] * g_j
        force_k = dEdphi[:, None] * g_k
        force_l = dEdphi[:, None] * g_l
        g = (scatter_add(force_i, arr.tor_i, n)
             + scatter_add(force_j, arr.tor_j, n)
             + scatter_add(force_k, arr.tor_k, n)
             + scatter_add(force_l, arr.tor_l, n))
        return e, g
