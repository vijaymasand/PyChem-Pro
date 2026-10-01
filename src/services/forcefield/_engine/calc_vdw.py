"""
MMFF94 van der Waals — buffered 7-14 potential (Halgren 1996).

x = r/rs
f1 = 1.07 / (x + 0.07)
f2 = 1.12 / (x^7 + 0.12)
E_vdw = eps * f1^7 * (f2 - 2)

dE/dr = -(7*eps*f1^7 / rs) * [(f1/1.07)(f2-2) + (f2^2 * x^6 / 1.12)]
"""
from __future__ import annotations
import numpy as np

from src.services.forcefield._engine.arrays import InteractionArrays
from src.services.forcefield._engine._math import scatter_add


class VdwCalc:

    @staticmethod
    def energy(coords: np.ndarray, arr: InteractionArrays) -> float:
        if arr.vdw_i.size == 0:
            return 0.0
        diff = coords[arr.vdw_i] - coords[arr.vdw_j]
        r = np.linalg.norm(diff, axis=1)
        rs = arr.vdw_rs
        eps = arr.vdw_eps
        x = r / np.maximum(rs, 1e-12)
        f1 = 1.07 / (x + 0.07)
        f2 = 1.12 / (x ** 7 + 0.12)
        return float((eps * (f1 ** 7) * (f2 - 2.0)).sum())

    @staticmethod
    def gradient(coords: np.ndarray, arr: InteractionArrays) -> np.ndarray:
        return VdwCalc.energy_and_gradient(coords, arr)[1]

    @staticmethod
    def energy_and_gradient(coords: np.ndarray, arr: InteractionArrays):
        n = arr.n_atoms
        if arr.vdw_i.size == 0:
            return 0.0, np.zeros((n, 3), dtype=np.float64)
        diff = coords[arr.vdw_i] - coords[arr.vdw_j]
        r = np.linalg.norm(diff, axis=1)
        r_safe = np.maximum(r, 1e-12)
        rs = arr.vdw_rs
        rs_safe = np.maximum(rs, 1e-12)
        eps = arr.vdw_eps
        x = r / rs_safe
        f1 = 1.07 / (x + 0.07)
        f2 = 1.12 / (x ** 7 + 0.12)
        f1_7 = f1 ** 7
        e_per = eps * f1_7 * (f2 - 2.0)
        e = float(e_per.sum())

        # dE/dr per Jmol MMFFVDWCalc:
        #   dE/dr = -7 eps f1^7 / rs * [ f1/1.07 * (f2 - 2) + f2^2 * x^6 / 1.12 ]
        x6 = x ** 6
        dEdr = -7.0 * eps * f1_7 / rs_safe * (
            (f1 / 1.07) * (f2 - 2.0) + f2 * f2 * x6 / 1.12
        )
        unit = diff / r_safe[:, None]
        force = dEdr[:, None] * unit
        g = scatter_add(force, arr.vdw_i, n) - scatter_add(force, arr.vdw_j, n)
        return e, g
