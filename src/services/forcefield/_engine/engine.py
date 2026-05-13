"""
MMFF94 orchestrator: sums the 7 per-term calculators into a single
energy_and_gradient call, plus a per-term decomposition for debugging
and Jmol snapshot diffing.
"""
from __future__ import annotations
from typing import Tuple, Dict
import numpy as np

from src.services.forcefield._engine.arrays import InteractionArrays
from src.services.forcefield._engine.calc_bond import BondCalc
from src.services.forcefield._engine.calc_angle import AngleCalc
from src.services.forcefield._engine.calc_stretchbend import StretchBendCalc
from src.services.forcefield._engine.calc_torsion import TorsionCalc
from src.services.forcefield._engine.calc_oop import OopCalc
from src.services.forcefield._engine.calc_vdw import VdwCalc
from src.services.forcefield._engine.calc_es import EsCalc


class MMFF94Engine:
    """Orchestrates the 7 MMFF94 term calculators.

    Construct once per optimize_geometry call (arrays are static for one
    call; rebuild between calls if molecule connectivity changes).
    """

    # Ordered tuple of (term name, calculator class)
    _TERMS = (
        ("bond", BondCalc),
        ("angle", AngleCalc),
        ("stretch_bend", StretchBendCalc),
        ("torsion", TorsionCalc),
        ("oop", OopCalc),
        ("vdw", VdwCalc),
        ("es", EsCalc),
    )

    def __init__(self, arrays: InteractionArrays):
        self.arrays = arrays

    def energy(self, coords: np.ndarray) -> float:
        return float(sum(calc.energy(coords, self.arrays) for _, calc in self._TERMS))

    def gradient(self, coords: np.ndarray) -> np.ndarray:
        g = np.zeros((self.arrays.n_atoms, 3), dtype=np.float64)
        for _, calc in self._TERMS:
            g += calc.gradient(coords, self.arrays)
        return g

    def energy_and_gradient(self, coords: np.ndarray) -> Tuple[float, np.ndarray]:
        e_total = 0.0
        g = np.zeros((self.arrays.n_atoms, 3), dtype=np.float64)
        for _, calc in self._TERMS:
            e, dg = calc.energy_and_gradient(coords, self.arrays)
            e_total += e
            g += dg
        return float(e_total), g

    def energy_decomposition(self, coords: np.ndarray) -> Dict[str, float]:
        """Per-term energies. Used for Jmol snapshot diffing — when totals
        disagree, decomposition pinpoints which term is off."""
        return {name: float(calc.energy(coords, self.arrays))
                for name, calc in self._TERMS}
