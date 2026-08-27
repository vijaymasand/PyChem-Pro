import logging
import numpy as np
from typing import Optional, Tuple, Callable, Union

from ..core.molecule import Molecule

try:
    from scipy.optimize import minimize
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

from ..forcefields.mmff94 import MMFF94
from ..semi_empirical.am1 import AM1
from ..semi_empirical.pm3 import PM3

logger = logging.getLogger(__name__)

class GeometryOptimizer:
    def __init__(self, max_iter: int = 200, gtol: float = 0.01, etol: float = 1e-6):
        self.max_iter = max_iter
        self.gtol = gtol
        self.etol = etol

    def optimize(self, molecule: Molecule, method: str = 'mmff94', **kwargs) -> Molecule:
        method = method.lower()
        if method == 'mmff94':
            calc = MMFF94(molecule)
        elif method == 'am1':
            calc = AM1(molecule)
        elif method == 'pm3':
            calc = PM3(molecule)
        else:
            raise ValueError(f"Unsupported method: {method}")

        coords = molecule.get_coordinates().flatten()

        def energy_and_grad(x):
            mol_copy = molecule.copy()
            mol_copy.set_coordinates(x.reshape(-1, 3))
            energy = calc.energy(mol_copy)
            grad = calc.gradient(mol_copy).flatten()
            return energy, grad

        if HAS_SCIPY:
            result = minimize(
                fun=lambda x: energy_and_grad(x)[0],
                jac=lambda x: energy_and_grad(x)[1],
                x0=coords,
                method='BFGS',
                options={'maxiter': self.max_iter, 'gtol': self.gtol}
            )
            optimized_coords = result.x
            final_energy = result.fun
            success = result.success
        else:
            optimized_coords, final_energy, success = self._gradient_descent(
                coords, energy_and_grad, self.max_iter, self.gtol, self.etol
            )

        if not success:
            logger.warning("Optimization did not converge")

        opt_mol = molecule.copy()
        opt_mol.set_coordinates(optimized_coords.reshape(-1, 3))
        opt_mol.energy = final_energy
        return opt_mol

    def _gradient_descent(self, x0: np.ndarray, func_grad: Callable,
                          max_iter: int, gtol: float, etol: float) -> Tuple[np.ndarray, float, bool]:
        x = x0.copy()
        f_prev, g = func_grad(x)
        for i in range(max_iter):
            g_norm = np.linalg.norm(g)
            if g_norm < gtol:
                return x, f_prev, True
            # simple line search with Armijo condition
            alpha = 1.0
            for _ in range(20):
                x_new = x - alpha * g
                f_new, g_new = func_grad(x_new)
                if f_new < f_prev - 1e-4 * alpha * np.dot(g, g):
                    break
                alpha *= 0.5
            else:
                return x, f_prev, False
            if abs(f_prev - f_new) < etol:
                return x_new, f_new, True
            x = x_new
            f_prev, g = f_new, g_new
        return x, f_prev, False
