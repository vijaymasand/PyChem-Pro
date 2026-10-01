"""
Vectorized Mesh Builder — High-speed Batch Protein Mesh Generation.

This module replaces the residue-by-residue loop with a massive NumPy batch process.
It calculates all splines, profiles, and triangle indices for an entire protein chain 
in one pass, utilizing vectorized CPU operations.
"""

import numpy as np
from typing import List, Tuple, Optional
from src.features.visualization_3d.services.peptide_plane import (
    PeptidePlane, HELIX, STRAND, COIL, transition
)
from src.features.visualization_3d.services.spline_math import (
    bspline_batch, linear, in_out_quad, out_circ
)
from src.features.visualization_3d.services.profiles import (
    segment_profiles, segment_colors
)

def generate_chain_mesh_vectorized(planes: List[PeptidePlane], spline_steps: int, profile_detail: int, theme_colors: dict = None):
    """
    Generate mesh for an entire chain using optimized vectorization.
    
    Performance optimizations for large proteins:
    - Progressive mesh generation with early termination
    - Memory-efficient array operations
    - Adaptive LOD based on chain size
    
    Args:
        planes: List of PeptidePlane objects for the chain
        spline_steps: LOD subdivision steps
        profile_detail: LOD profile points
    
    Returns: (vertices, triangles, colors)
    """
    n_planes = len(planes)
    n_segments = n_planes - 3
    if n_segments <= 0:
        return None, None, None
    
    # Performance optimization: Early termination for very large proteins
    MAX_SEGMENTS = 500  # Limit segments to prevent memory explosion
    if n_segments > MAX_SEGMENTS:
        # Use simplified mesh for very large proteins
        return _generate_simplified_mesh(planes, spline_steps, profile_detail, theme_colors)

    # 1. Prepare Batch Data
    # Extract positions, sides, and normals for all 4-plane windows
    # Shape: (N_segments, 4, 3)
    pos = np.array([p.position for p in planes])
    side = np.array([p.side for p in planes])
    norm = np.array([p.normal for p in planes])
    
    # Indices for the 4 planes in each sliding window
    idx1 = np.arange(n_segments)
    idx2 = idx1 + 1
    idx3 = idx1 + 2
    idx4 = idx1 + 3
    
    # Extract profiles for all segments
    # For speed, we assume profile points (u, v) are the same for all segments of the same SS type.
    # To be perfectly faithful to CPPCartoon, we'll get them per segment.
    u1_batch = []
    v1_batch = []
    u2_batch = []
    v2_batch = []
    colors1 = []
    colors2 = []
    easings = []
    
    for i in range(n_segments):
        pp2 = planes[idx2[i]]
        pp3 = planes[idx3[i]]
        p1, p2 = segment_profiles(pp2, pp3, profile_detail)
        u1_batch.append(p1[:, 0])
        v1_batch.append(p1[:, 1])
        u2_batch.append(p2[:, 0])
        v2_batch.append(p2[:, 1])
        
        c1, c2 = segment_colors(pp2, theme_colors)
        colors1.append(c1)
        colors2.append(c2)
        
        # Select easing function
        type0 = pp2.ss1
        type1, type2 = transition(pp2)
        if type1 == STRAND and type2 != STRAND:
            easings.append(0) # linear
        elif type0 == STRAND and type1 != STRAND:
            easings.append(2) # out_circ
        else:
            easings.append(1) # in_out_quad

    u1 = np.array(u1_batch) # (N, P)
    v1 = np.array(v1_batch)
    u2 = np.array(u2_batch)
    v2 = np.array(v2_batch)
    c1 = np.array(colors1) # (N, 3)
    c2 = np.array(colors2)
    easings = np.array(easings) # (N,)

    # 2. Batch Spline Calculation
    # Guide points: g = pos + side * u + norm * v
    # pos[idx1] is (N, 3), side[idx1] is (N, 3), u1 is (N, P)
    # broadcasting: (N, 1, 3) + (N, P, 1) * (N, 1, 3) -> (N, P, 3)
    def get_guide(idx, u_in, v_in):
        p_ = pos[idx][:, np.newaxis, :]
        s_ = side[idx][:, np.newaxis, :]
        n_ = norm[idx][:, np.newaxis, :]
        uu = u_in[:, :, np.newaxis]
        vv = v_in[:, :, np.newaxis]
        return p_ + s_ * uu + n_ * vv

    g1_1 = get_guide(idx1, u1, v1)
    g1_2 = get_guide(idx2, u1, v1)
    g1_3 = get_guide(idx3, u1, v1)
    g1_4 = get_guide(idx4, u1, v1)
    
    g2_1 = get_guide(idx1, u2, v2)
    g2_2 = get_guide(idx2, u2, v2)
    g2_3 = get_guide(idx3, u2, v2)
    g2_4 = get_guide(idx4, u2, v2)

    # splines1/2 shape: (S+1, N, P, 3)
    splines1 = bspline_batch(g1_1, g1_2, g1_3, g1_4, spline_steps)
    splines2 = bspline_batch(g2_1, g2_2, g2_3, g2_4, spline_steps)
    
    S = spline_steps
    P = profile_detail
    N = n_segments

    # 3. Vectorized Blending and Mesh Construction
    # t_all: precalculate all easing values for all segments
    # t shape: (N, S+1)
    t_vals = np.linspace(0, 1, S + 1)
    t_all = np.zeros((N, S + 1))
    
    # Apply easing based on segment type
    t_all[easings == 0] = t_vals # linear
    t_all[easings == 1] = np.array([in_out_quad(t) for t in t_vals])
    t_all[easings == 2] = np.array([out_circ(t) for t in t_vals])
    
    # Reshape t for broadcasting: (S+1, N, 1, 1)
    t = t_all.T[:, :, np.newaxis, np.newaxis]
    
    # Vertices: (S+1, N, P, 3)
    v_grid = splines1 * (1.0 - t) + splines2 * t
    
    # Colors: (S+1, N, P, 3)
    # c1, c2 are (N, 3)
    # t is (S+1, N, 1, 1)
    c1_ = c1[np.newaxis, :, np.newaxis, :]
    c2_ = c2[np.newaxis, :, np.newaxis, :]
    c_grid = c1_ * (1.0 - t) + c2_ * t
    c_grid = np.broadcast_to(c_grid, (S + 1, N, P, 3))

    # Reshape to final flat buffers
    # Vertices and Colors are contiguous blocks per segment
    flat_verts = v_grid.transpose(1, 0, 2, 3).reshape(-1, 3)
    flat_colors = c_grid.transpose(1, 0, 2, 3).reshape(-1, 3)
    
    # 4. Generate Triangle Indices (Mathematically)
    # Each segment has (S+1) * P vertices.
    # Total vertices = N * (S+1) * P
    # Triangles per segment: S * P * 2 (quads)
    
    # Base indices for a single segment grid
    i_idx = np.arange(S)[:, np.newaxis]
    j_idx = np.arange(P)
    j_next = (j_idx + 1) % P
    
    row_start = i_idx * P
    next_row_start = (i_idx + 1) * P
    
    p00 = row_start + j_idx
    p01 = next_row_start + j_idx
    p10 = row_start + j_next
    p11 = next_row_start + j_next
    
    tri1_base = np.stack([p10, p11, p01], axis=-1).reshape(-1, 3)
    tri2_base = np.stack([p10, p01, p00], axis=-1).reshape(-1, 3)
    tri_base = np.vstack([tri1_base, tri2_base]) # (S*P*2, 3)
    
    # Replicate for all segments with offsets
    # Segment i starts at index i * (S+1) * P
    offsets = np.arange(N) * (S + 1) * P
    triangles = tri_base[np.newaxis, :, :] + offsets[:, np.newaxis, np.newaxis]
    triangles = triangles.reshape(-1, 3).astype(np.int32)
    
    # Return directly to ensure correctness; the renderer handles small gaps.
    return flat_verts, triangles, flat_colors


def _generate_simplified_mesh(planes: List[PeptidePlane], spline_steps: int, profile_detail: int, theme_colors: dict = None):
    """
    Generate simplified mesh for very large proteins to prevent memory explosion.
    
    Uses adaptive sampling and simplified geometry for better performance.
    """
    n_planes = len(planes)
    if n_planes <= 10:
        # Small proteins use normal vectorized approach
        return generate_chain_mesh_vectorized(planes, spline_steps, profile_detail, theme_colors)
    
    # For large proteins, use progressive sampling
    
    # Sample every N-th plane to reduce complexity
    SAMPLE_RATE = max(1, n_planes // 100)  # Sample at most 100 planes
    sampled_indices = np.arange(0, n_planes, SAMPLE_RATE)
    sampled_planes = [planes[i] for i in sampled_indices]
    
    # Generate simplified mesh
    if len(sampled_planes) < 4:
        return None, None, None
    
    return generate_chain_mesh_vectorized(sampled_planes, spline_steps, max(3, profile_detail // 2), theme_colors)
