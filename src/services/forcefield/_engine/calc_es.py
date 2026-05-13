"""
MMFF94 electrostatic — buffered Coulomb with 1-4 scaling.

E_es = (qi * qj * f) / (r + 0.05)
       where f = 249.0537 for 1-4 pairs, 332.0716 for 1-5+ pairs.
dE/dr = -E_es / (r + 0.05)

Legacy code (mmff94_service.py:450) used 332.0716 for all pairs —
missed the 0.75x scaling on 1-4 pairs that distinguishes MMFF94
from a naive Coulomb potential.
"""
from __future__ import annotations
import numpy as np

from src.services.forcefield._engine.arrays import InteractionArrays
from src.services.forcefield._engine._math import scatter_add

_BUF = 0.05


class EsCalc:

    @staticmethod
    def energy(coords: np.ndarray, arr: InteractionArrays) -> float:
        if arr.es_i.size == 0:
            return 0.0
        diff = coords[arr.es_i] - coords[arr.es_j]
        r = np.linalg.norm(diff, axis=1)
        return float((arr.es_qq * arr.es_factor / (r + _BUF)).sum())

    @staticmethod
    def gradient(coords: np.ndarray, arr: InteractionArrays) -> np.ndarray:
        return EsCalc.energy_and_gradient(coords, arr)[1]

    @staticmethod
    def energy_and_gradient(coords: np.ndarray, arr: InteractionArrays):
        n = arr.n_atoms
        if arr.es_i.size == 0:
            return 0.0, np.zeros((n, 3), dtype=np.float64)
        diff = coords[arr.es_i] - coords[arr.es_j]
        r = np.linalg.norm(diff, axis=1)
        r_safe = np.maximum(r, 1e-12)
        r_buf = r + _BUF
        e_per = arr.es_qq * arr.es_factor / r_buf
        e = float(e_per.sum())
        dEdr = -e_per / r_buf
        unit = diff / r_safe[:, None]
        force = dEdr[:, None] * unit
        g = scatter_add(force, arr.es_i, n) - scatter_add(force, arr.es_j, n)
        return e, g
