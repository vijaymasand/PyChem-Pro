"""
MMFF94 bond stretching (Halgren 1996) — corrected cubic-quartic form.

E_b = (143.9325/2) * kb * d^2 * (1 + CS*d + (7/12)*CS^2*d^2),   CS = -2.0, d = r - r0
     = 71.96625 * kb * d^2 * (1 - 2*d + (7/3)*d^2)

dE_b/dr = 143.9325 * kb * d * (1 - 3*d + (14/3)*d^2)
"""
from __future__ import annotations
import numpy as np

from src.services.forcefield._engine.arrays import InteractionArrays
from src.services.forcefield._engine._math import scatter_add

_C0 = 71.96625      # 143.9325 / 2
_C1 = 143.9325      # dE prefactor
_7_OVER_3 = 7.0 / 3.0
_14_OVER_3 = 14.0 / 3.0


class BondCalc:
    """Stateless: only static methods."""

    @staticmethod
    def energy(coords: np.ndarray, arr: InteractionArrays) -> float:
        if arr.bond_i.size == 0:
            return 0.0
        diff = coords[arr.bond_i] - coords[arr.bond_j]
        r = np.linalg.norm(diff, axis=1)
        delta = r - arr.bond_r0
        d2 = delta * delta
        e = _C0 * arr.bond_kb * d2 * (1.0 - 2.0 * delta + _7_OVER_3 * d2)
        return float(e.sum())

    @staticmethod
    def gradient(coords: np.ndarray, arr: InteractionArrays) -> np.ndarray:
        n = arr.n_atoms
        if arr.bond_i.size == 0:
            return np.zeros((n, 3), dtype=np.float64)
        diff = coords[arr.bond_i] - coords[arr.bond_j]
        r = np.linalg.norm(diff, axis=1)
        # Avoid division by zero at coincident atoms (numerically rare)
        r_safe = np.maximum(r, 1e-12)
        delta = r - arr.bond_r0
        d2 = delta * delta
        dEdr = _C1 * arr.bond_kb * delta * (1.0 - 3.0 * delta + _14_OVER_3 * d2)
        unit = diff / r_safe[:, None]
        force = dEdr[:, None] * unit                  # (Nb, 3)
        g = scatter_add(force, arr.bond_i, n)
        g -= scatter_add(force, arr.bond_j, n)
        return g

    @staticmethod
    def energy_and_gradient(coords: np.ndarray, arr: InteractionArrays):
        """Fused: shares r and diff between energy and gradient."""
        n = arr.n_atoms
        if arr.bond_i.size == 0:
            return 0.0, np.zeros((n, 3), dtype=np.float64)
        diff = coords[arr.bond_i] - coords[arr.bond_j]
        r = np.linalg.norm(diff, axis=1)
        r_safe = np.maximum(r, 1e-12)
        delta = r - arr.bond_r0
        d2 = delta * delta
        e = float((_C0 * arr.bond_kb * d2 * (1.0 - 2.0 * delta + _7_OVER_3 * d2)).sum())
        dEdr = _C1 * arr.bond_kb * delta * (1.0 - 3.0 * delta + _14_OVER_3 * d2)
        unit = diff / r_safe[:, None]
        force = dEdr[:, None] * unit
        g = scatter_add(force, arr.bond_i, n) - scatter_add(force, arr.bond_j, n)
        return e, g
