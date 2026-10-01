# src/services/rendering/parallel_projection.py
"""
Parallel pre-render pipeline for molecular visualization.

Parallelizes projection, depth sort, and culling across CPU cores.
Only activates for molecules with 500+ atoms (overhead not worth it for smaller ones).
"""
import numpy as np
from src.core.parallel import ParallelExecutor


def _build_rotation_matrix(rot_x_deg: float, rot_y_deg: float) -> np.ndarray:
    """Build 3x3 rotation matrix from X and Y rotation angles in degrees."""
    rx = np.radians(rot_x_deg)
    ry = np.radians(rot_y_deg)
    cos_x, sin_x = np.cos(rx), np.sin(rx)
    cos_y, sin_y = np.cos(ry), np.sin(ry)
    Rx = np.array([[1, 0, 0], [0, cos_x, -sin_x], [0, sin_x, cos_x]])
    Ry = np.array([[cos_y, 0, sin_y], [0, 1, 0], [-sin_y, 0, cos_y]])
    return Ry @ Rx


# Module-level worker functions (picklable for spawn multiprocessing)

def _project_chunk(args):
    """Project a chunk of 3D coordinates to 2D screen space."""
    coords, rot_matrix, cx, cy, zoom = args
    rotated = coords @ rot_matrix.T
    screen_x = rotated[:, 0] * zoom + cx
    screen_y = -rotated[:, 1] * zoom + cy  # Y-axis flipped for screen
    depth_z = rotated[:, 2]
    return np.column_stack([screen_x, screen_y, depth_z])


def _cull_chunk(args):
    """Return boolean mask of atoms visible within viewport."""
    projected, viewport_w, viewport_h, margin = args
    visible = (
        (projected[:, 0] > -margin) &
        (projected[:, 0] < viewport_w + margin) &
        (projected[:, 1] > -margin) &
        (projected[:, 1] < viewport_h + margin)
    )
    return visible


def _sort_chunk(args):
    """Sort atom indices by depth within a chunk."""
    depths, indices = args
    order = np.argsort(depths)
    return indices[order], depths[order]


# Public API functions

def project_atoms(coords: np.ndarray, rot_x: float, rot_y: float,
                  center_x: float, center_y: float, zoom: float,
                  executor: ParallelExecutor = None) -> np.ndarray:
    """
    Project 3D atom coordinates to 2D screen space.

    Args:
        coords: Nx3 array of atom positions
        rot_x, rot_y: Rotation angles in degrees
        center_x, center_y: Screen center
        zoom: Zoom factor
        executor: ParallelExecutor (uses parallel for N >= 500)

    Returns:
        Nx3 array of (screen_x, screen_y, depth_z)
    """
    rot_matrix = _build_rotation_matrix(rot_x, rot_y)
    n = len(coords)

    if n < 500 or executor is None:
        return _project_chunk((coords, rot_matrix, center_x, center_y, zoom))

    # Split into chunks for parallel processing
    chunk_size = max(100, n // executor.num_workers)
    chunks = []
    for i in range(0, n, chunk_size):
        chunk_coords = coords[i:i+chunk_size]
        chunks.append((chunk_coords, rot_matrix, center_x, center_y, zoom))

    results = executor.map(_project_chunk, chunks)
    return np.vstack(results)


def cull_offscreen(projected: np.ndarray, viewport_w: int, viewport_h: int,
                   margin: int = 50, executor: ParallelExecutor = None) -> np.ndarray:
    """
    Return boolean mask of atoms visible within viewport + margin.
    """
    n = len(projected)

    if n < 500 or executor is None:
        return _cull_chunk((projected, viewport_w, viewport_h, margin))

    chunk_size = max(100, n // executor.num_workers)
    chunks = []
    for i in range(0, n, chunk_size):
        chunks.append((projected[i:i+chunk_size], viewport_w, viewport_h, margin))

    results = executor.map(_cull_chunk, chunks)
    return np.concatenate(results)


def depth_sort(projected: np.ndarray, executor: ParallelExecutor = None) -> np.ndarray:
    """
    Sort atom indices by depth (back to front) for painter's algorithm.
    Returns sorted indices.
    """
    n = len(projected)
    indices = np.arange(n)

    if n < 1000 or executor is None:
        return indices[np.argsort(projected[:, 2])]

    # Parallel: sort chunks then merge
    chunk_size = max(200, n // executor.num_workers)
    chunks = []
    for i in range(0, n, chunk_size):
        end = min(i + chunk_size, n)
        chunks.append((projected[i:end, 2], indices[i:end]))

    sorted_chunks = executor.map(_sort_chunk, chunks)

    # Merge sorted chunks
    all_indices = np.concatenate([c[0] for c in sorted_chunks])
    all_depths = np.concatenate([c[1] for c in sorted_chunks])
    final_order = np.argsort(all_depths)
    return all_indices[final_order]
