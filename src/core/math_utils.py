# src/core/math_utils.py
import numpy as np

def calculate_distance(p1, p2, epsilon=1e-8):
    """
    Standardized Euclidean distance with safety epsilon for Force Field calculations.
    Ensures stability during energy minimization by preventing division by zero.
    """
    # Using np.subtract and np.linalg.norm is optimized for large NumPy arrays
    return np.linalg.norm(np.subtract(p1, p2)) + epsilon