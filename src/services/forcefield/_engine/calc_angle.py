"""
MMFF94 angle bending (Halgren 1996) — degree-based with linear-angle special case.

Non-linear: E_a = (0.043844/2) * ka * dt^2 * (1 + CB*dt)
           dt = theta_deg - theta0_deg,  CB = -0.007 per degree

Linear (theta0=180): E_a = 143.9325 * ka * (1 + cos theta_rad)
                  dE/dtheta_rad = -143.9325 * ka * sin theta_rad
"""
from __future__ import annotations
import numpy as np

from src.services.forcefield._engine.arrays import InteractionArrays
from src.services.forcefield._engine._math import (
    scatter_add, batch_angle, angle_jacobian,
)

_HALF_FA = 0.021922          # 0.043844 / 2
_FA = 0.043844
_CB = -0.007                 # per degree
_RAD_TO_DEG = 180.0 / np.pi
_FA_LIN = 143.9325


class AngleCalc:

    @staticmethod
    def energy(coords: np.ndarray, arr: InteractionArrays) -> float:
        if arr.angle_i.size == 0:
            return 0.0
        theta_rad = batch_angle(coords, arr.angle_i, arr.angle_j, arr.angle_k)
        # Non-linear part
        not_lin = ~arr.angle_is_linear
        e_total = 0.0
        if not_lin.any():
            t_deg = theta_rad[not_lin] * _RAD_TO_DEG
            dt = t_deg - arr.angle_theta0_deg[not_lin]
            ka_nl = arr.angle_ka[not_lin]
            e_total += float(np.sum(_HALF_FA * ka_nl * dt * dt * (1.0 + _CB * dt)))
        # Linear part
        if arr.angle_is_linear.any():
            t_lin = theta_rad[arr.angle_is_linear]
            ka_l = arr.angle_ka[arr.angle_is_linear]
            e_total += float(np.sum(_FA_LIN * ka_l * (1.0 + np.cos(t_lin))))
        return e_total

    @staticmethod
    def gradient(coords: np.ndarray, arr: InteractionArrays) -> np.ndarray:
        return AngleCalc.energy_and_gradient(coords, arr)[1]

    @staticmethod
    def energy_and_gradient(coords: np.ndarray, arr: InteractionArrays):
        n = arr.n_atoms
        if arr.angle_i.size == 0:
            return 0.0, np.zeros((n, 3), dtype=np.float64)

        theta_rad = batch_angle(coords, arr.angle_i, arr.angle_j, arr.angle_k)
        # dE/dtheta_rad per angle
        dEdtheta = np.zeros_like(theta_rad)
        e_total = 0.0

        not_lin = ~arr.angle_is_linear
        if not_lin.any():
            t_deg = theta_rad[not_lin] * _RAD_TO_DEG
            dt = t_deg - arr.angle_theta0_deg[not_lin]
            ka_nl = arr.angle_ka[not_lin]
            # E = 0.021922 * ka * dt^2 * (1 + CB*dt)
            e_total += float(np.sum(_HALF_FA * ka_nl * dt * dt * (1.0 + _CB * dt)))
            # dE/d(dt) = 0.021922 * ka * (2*dt + 3*CB*dt^2)
            dE_ddt = _HALF_FA * ka_nl * dt * (2.0 + 3.0 * _CB * dt)
            # dt is in degrees, theta in radians: d(dt)/d(theta_rad) = RAD_TO_DEG
            dEdtheta[not_lin] = dE_ddt * _RAD_TO_DEG

        if arr.angle_is_linear.any():
            t_lin = theta_rad[arr.angle_is_linear]
            ka_l = arr.angle_ka[arr.angle_is_linear]
            e_total += float(np.sum(_FA_LIN * ka_l * (1.0 + np.cos(t_lin))))
            dEdtheta[arr.angle_is_linear] = -_FA_LIN * ka_l * np.sin(t_lin)

        # Angle Jacobian (dtheta/dr_x)
        g_i, g_j, g_k = angle_jacobian(coords, arr.angle_i, arr.angle_j, arr.angle_k)
        force_i = dEdtheta[:, None] * g_i
        force_j = dEdtheta[:, None] * g_j
        force_k = dEdtheta[:, None] * g_k
        g = (scatter_add(force_i, arr.angle_i, n)
             + scatter_add(force_j, arr.angle_j, n)
             + scatter_add(force_k, arr.angle_k, n))
        return e_total, g
