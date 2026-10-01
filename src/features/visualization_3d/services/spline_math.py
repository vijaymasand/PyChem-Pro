"""
Spline Math — Faithful port of CPPCartoon's Spline.hpp

Uses uniform cubic B-spline (NOT Catmull-Rom) with forward differencing
for efficient evaluation. The key function is spline_for_planes() which
creates a separate spline for each profile vertex through 4 peptide planes.

Also includes easing functions for smooth profile transitions.
"""

import numpy as np


# ─── B-Spline Core ───────────────────────────────────────────────────────────

def bspline(vec1, vec2, vec3, vec4, n):
    """
    Uniform cubic B-spline through 4 control points using forward differencing.
    Matches CPPCartoon Spline.hpp lines 17-63 exactly.
    
    Args:
        vec1..vec4: control points (np.ndarray shape (3,))
        n: number of subdivisions (splineSteps)
    
    Returns:
        np.ndarray shape (n+1, 3) — the spline points
    """
    n1 = float(n)
    n2 = n1 * n1
    n3 = n1 * n1 * n1

    # Forward difference matrix S
    # s = [[6/n3,    0,    0, 0],
    #      [6/n3, 2/n2,    0, 0],
    #      [1/n3, 1/n2, 1/n1, 0],
    #      [   0,    0,    0, 1]]
    S = np.zeros((4, 4))
    S[0, 0] = 6.0 / n3
    S[1, 0] = 6.0 / n3;  S[1, 1] = 2.0 / n2
    S[2, 0] = 1.0 / n3;  S[2, 1] = 1.0 / n2;  S[2, 2] = 1.0 / n1
    S[3, 3] = 1.0

    # B-spline basis matrix (uniform cubic)
    oos = 1.0 / 6.0
    B = np.zeros((4, 4))
    B[0, 0] = -1*oos;  B[0, 1] =  3*oos;  B[0, 2] = -3*oos;  B[0, 3] = 1*oos
    B[1, 0] =  3*oos;  B[1, 1] = -6*oos;  B[1, 2] =  3*oos;  B[1, 3] = 0
    B[2, 0] = -3*oos;  B[2, 1] =  0;      B[2, 2] =  3*oos;  B[2, 3] = 0
    B[3, 0] =  1*oos;  B[3, 1] =  4*oos;  B[3, 2] =  1*oos;  B[3, 3] = 0

    # Geometry matrix G (4 control points as rows, 3 coords + homogeneous w=1)
    G = np.ones((4, 4))
    G[0, :3] = vec1
    G[1, :3] = vec2
    G[2, :3] = vec3
    G[3, :3] = vec4

    # Combined matrix M = S · B · G
    M = S @ B @ G

    # Generate points via forward differencing
    result = np.empty((n + 1, 3))

    # First point
    w = M[3, 3]
    if abs(w) < 1e-15:
        w = 1.0
    result[0] = M[3, :3] / w

    for k in range(n):
        M[3] += M[2]
        M[2] += M[1]
        M[1] += M[0]
        w = M[3, 3]
        if abs(w) < 1e-15:
            w = 1.0
        result[k + 1] = M[3, :3] / w

    return result


def bspline_batch(vec1, vec2, vec3, vec4, n):
    """
    Batched uniform cubic B-spline through 4 batched control points.
    Vectorized version of bspline() for arbitrary leading batch dimensions.

    Args:
        vec1..vec4: control point batches (np.ndarray shape (N, P, 3))
        n: number of subdivisions (splineSteps)

    Returns:
        np.ndarray shape (n+1, N, P, 3) — spline points for each (segment, profile) pair
    """
    n1 = float(n)
    n2 = n1 * n1
    n3 = n1 * n1 * n1

    S = np.zeros((4, 4))
    S[0, 0] = 6.0 / n3
    S[1, 0] = 6.0 / n3;  S[1, 1] = 2.0 / n2
    S[2, 0] = 1.0 / n3;  S[2, 1] = 1.0 / n2;  S[2, 2] = 1.0 / n1
    S[3, 3] = 1.0

    oos = 1.0 / 6.0
    B = np.zeros((4, 4))
    B[0, 0] = -1*oos;  B[0, 1] =  3*oos;  B[0, 2] = -3*oos;  B[0, 3] = 1*oos
    B[1, 0] =  3*oos;  B[1, 1] = -6*oos;  B[1, 2] =  3*oos;  B[1, 3] = 0
    B[2, 0] = -3*oos;  B[2, 1] =  0;      B[2, 2] =  3*oos;  B[2, 3] = 0
    B[3, 0] =  1*oos;  B[3, 1] =  4*oos;  B[3, 2] =  1*oos;  B[3, 3] = 0

    SB = S @ B  # (4, 4)

    # Stack control points: (4, N, P, 3), flatten batch dims for vectorized matmul
    G = np.stack([vec1, vec2, vec3, vec4], axis=0)
    batch_shape = G.shape[1:]          # (N, P, 3)
    G_flat = G.reshape(4, -1, 3)       # (4, N*P, 3)

    # M[i, j] = sum_k SB[i,k] * G_flat[k, j]  →  (4, N*P, 3)
    # w is always 1 for position data so the homogeneous divide is a no-op
    M = np.einsum('ij,jkl->ikl', SB, G_flat)

    n_batch = G_flat.shape[1]
    result = np.empty((n + 1, n_batch, 3))
    result[0] = M[3]

    for k in range(n):
        M[3] += M[2]
        M[2] += M[1]
        M[1] += M[0]
        result[k + 1] = M[3]

    return result.reshape((n + 1,) + batch_shape)


def spline_for_planes(pp1, pp2, pp3, pp4, n, u, v):
    """
    Create a B-spline for a specific profile point (u, v) through 4 peptide planes.
    Matches CPPCartoon Spline.hpp splineForPlanes() exactly.
    
    Each guide point = plane.Position + plane.Side * u + plane.Normal * v
    
    Args:
        pp1..pp4: PeptidePlane objects
        n: spline steps
        u, v: profile point coordinates (local x, y)
    
    Returns:
        np.ndarray shape (n+1, 3) — the spline points
    """
    g1 = pp1.position + pp1.side * u + pp1.normal * v
    g2 = pp2.position + pp2.side * u + pp2.normal * v
    g3 = pp3.position + pp3.side * u + pp3.normal * v
    g4 = pp4.position + pp4.side * u + pp4.normal * v
    return bspline(g1, g2, g3, g4, n)


# ─── Easing Functions ────────────────────────────────────────────────────────
# Matches CPPCartoon Spline.hpp lines 74-91

def linear(t):
    """Linear easing."""
    return t


def in_out_quad(t):
    """Quadratic ease-in-out. Used for smooth profile transitions."""
    if t < 0.5:
        return 2.0 * t * t
    t = 2.0 * t - 1.0
    return -0.5 * (t * (t - 2.0) - 1.0)


def out_circ(t):
    """Circular ease-out. Used for strand-to-coil transitions."""
    t = t - 1.0
    return np.sqrt(1.0 - t * t)


def in_circ(t):
    """Circular ease-in."""
    return -1.0 * (np.sqrt(1.0 - t * t) - 1.0)
