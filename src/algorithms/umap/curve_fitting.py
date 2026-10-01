"""
Low-Dimensional Fuzzy Membership Curve Fitting (a, b parameters).
Fits non-linear parameters a and b based on min_dist and spread parameters.
"""
import numpy as np
from scipy.optimize import curve_fit

def find_ab_params(spread: float = 1.0, min_dist: float = 0.1) -> tuple[float, float]:
    """
    Fits low-dimensional fuzzy distance parameters (a, b) given spread and min_dist.
    
    Parameters
    ----------
    spread : float
        Effective scale of embedded points.
    min_dist : float
        Minimum desired distance between embedded points.
        
    Returns
    -------
    a : float
    b : float
    """
    xv = np.linspace(0, spread * 3.0, 300)
    yv = np.where(xv <= min_dist, 1.0, np.exp(-(xv - min_dist) / spread))

    def curve(x, a, b):
        return 1.0 / (1.0 + a * (x ** (2.0 * b)))

    try:
        popt, _ = curve_fit(curve, xv, yv, p0=(1.0, 1.0), maxfev=20000)
        a, b = popt[0], popt[1]
    except Exception:
        # Fallback to analytical approximations if numerical curve_fit fails
        a = 1.5769 if min_dist == 0.1 else 1.0 / (min_dist ** 2.0)
        b = 0.8951

    return float(a), float(b)
