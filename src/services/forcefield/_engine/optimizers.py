"""
Optimizers for MMFF94 (and any other force field with energy + gradient).

Both optimizers are decoupled from MMFF94 — they accept a generic
eg_fn: coords -> (E_float, grad_ndarray) callable. Useful for testing
against synthetic potentials before wiring up the engine.
"""
from __future__ import annotations
from typing import Callable, Tuple, List
import numpy as np

EgFn = Callable[[np.ndarray], Tuple[float, np.ndarray]]


def _rms(g: np.ndarray) -> float:
    return float(np.sqrt(np.mean(g * g)))


class SteepestDescent:
    """Adaptive-step steepest descent.

    Per step:
        1. Compute energy and gradient.
        2. Check RMS gradient convergence.
        3. Take a trial step of size `scale = min(step, 0.1 / max|g|)`.
        4. Accept if energy decreased; grow step by 1.2x (cap 0.1).
           Reject otherwise; shrink step by 0.5x.
        5. Stop if step shrinks below 1e-10.
    """

    def __init__(self, eg_fn: EgFn):
        self.eg = eg_fn

    def run(self, coords: np.ndarray, max_iters: int = 500,
            convergence: float = 1e-4) -> Tuple[np.ndarray, List[float], bool, int]:
        coords = coords.astype(np.float64, copy=True)
        step = 0.01
        e, g = self.eg(coords)
        traj = [float(e)]

        for k in range(max_iters):
            if _rms(g) < convergence:
                return coords, traj, True, k

            max_g = float(np.max(np.abs(g)))
            if max_g <= 0.0:
                return coords, traj, True, k

            scale = min(step, 0.1 / max_g)
            new_coords = coords - scale * g
            new_e, new_g = self.eg(new_coords)

            if new_e < e:
                coords, e, g = new_coords, new_e, new_g
                step = min(step * 1.2, 0.1)
            else:
                step *= 0.5

            traj.append(float(e))

            if step < 1e-10:
                break

        return coords, traj, False, max_iters


class LBFGS:
    """Limited-memory BFGS with strong Wolfe line search.

    Two-loop recursion for direction; backtracking with Armijo (c1=1e-4)
    and curvature (c2=0.9) conditions for step size. Memory m=10 by default.
    """

    def __init__(self, eg_fn: EgFn, memory: int = 10):
        self.eg = eg_fn
        self.m = memory

    def run(self, coords: np.ndarray, max_iters: int = 500,
            convergence: float = 1e-4) -> Tuple[np.ndarray, List[float], bool, int]:
        shape = coords.shape
        x = coords.astype(np.float64, copy=True).flatten()
        e, g_nd = self.eg(x.reshape(shape))
        g = g_nd.flatten()

        s_list: List[np.ndarray] = []
        y_list: List[np.ndarray] = []
        rho_list: List[float] = []
        traj = [float(e)]

        for k in range(max_iters):
            if _rms(g.reshape(shape)) < convergence:
                return x.reshape(shape), traj, True, k

            # Two-loop recursion
            q = g.copy()
            alphas: List[float] = []
            for i in range(len(s_list) - 1, -1, -1):
                a = rho_list[i] * float(s_list[i].dot(q))
                alphas.insert(0, a)
                q -= a * y_list[i]
            if s_list:
                yy = float(y_list[-1].dot(y_list[-1]))
                gamma = float(s_list[-1].dot(y_list[-1])) / (yy + 1e-12)
                r = gamma * q
            else:
                r = 0.01 * q
            for i in range(len(s_list)):
                beta = rho_list[i] * float(y_list[i].dot(r))
                r += s_list[i] * (alphas[i] - beta)
            direction = -r

            # Strong Wolfe line search
            slope = float(g.dot(direction))
            if slope >= 0:
                direction = -g
                slope = float(g.dot(direction))
            step_size = 1.0
            c1, c2 = 1e-4, 0.9
            new_e = e
            new_g = g
            for _ in range(20):
                new_x = x + step_size * direction
                new_e, new_g_nd = self.eg(new_x.reshape(shape))
                new_g = new_g_nd.flatten()
                # Armijo
                if new_e > e + c1 * step_size * slope:
                    step_size *= 0.5
                    continue
                # Curvature
                new_slope = float(new_g.dot(direction))
                if new_slope < c2 * slope:
                    step_size *= 1.5
                    continue
                break
            else:
                # Line search exhausted — accept whatever we have
                pass

            s = new_x - x
            y = new_g - g
            sy = float(s.dot(y))
            if sy > 1e-12:
                s_list.append(s)
                y_list.append(y)
                rho_list.append(1.0 / sy)
                if len(s_list) > self.m:
                    s_list.pop(0)
                    y_list.pop(0)
                    rho_list.pop(0)

            x, e, g = new_x, new_e, new_g
            traj.append(float(e))

        return x.reshape(shape), traj, False, max_iters
