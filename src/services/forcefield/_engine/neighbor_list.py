"""
Cell-list neighbor pair builder for protein-scale MMFF94.

For n_atoms >= 200, full N^2 enumeration of VdW/ES pairs becomes wasteful
(4e4 pairs for a 200-atom protein where >95% are beyond any meaningful
cutoff). This module builds a 3D grid of cells of side length = cutoff,
then enumerates pairs only within and between adjacent cells. Built
once per optimize_geometry call; reused across all steps.

Cutoff defaults (MMFF94 convention):
    - 10 A for VdW
    - 15 A for ES
"""
from __future__ import annotations
import numpy as np

VDW_CUTOFF = 10.0
ES_CUTOFF = 15.0


def build_neighbor_pairs(coords: np.ndarray, cutoff: float) -> np.ndarray:
    """Return all unique pairs (i, j) with i < j and ||coords[i]-coords[j]|| <= cutoff.

    Args:
        coords: (n_atoms, 3) float64.
        cutoff: distance threshold in A.

    Returns:
        (Npairs, 2) int32. May be empty.
    """
    n = coords.shape[0]
    if n < 2:
        return np.zeros((0, 2), dtype=np.int32)

    # Cell side = cutoff. Cells indexed by (gx, gy, gz).
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    extent = maxs - mins
    nx = max(1, int(np.ceil(extent[0] / cutoff)))
    ny = max(1, int(np.ceil(extent[1] / cutoff)))
    nz = max(1, int(np.ceil(extent[2] / cutoff)))

    # Atom -> cell index
    cell_xyz = np.floor((coords - mins) / cutoff).astype(np.int32)
    np.clip(cell_xyz, 0, [nx - 1, ny - 1, nz - 1], out=cell_xyz)

    # cell_key = gx + nx*gy + nx*ny*gz
    keys = (cell_xyz[:, 0]
            + nx * cell_xyz[:, 1]
            + nx * ny * cell_xyz[:, 2])

    # Build cell -> list of atoms (using sort + segment).
    order = np.argsort(keys, kind="stable")
    sorted_keys = keys[order]
    # boundaries[k] = first index in sorted_keys with key == k
    n_cells = nx * ny * nz
    boundaries = np.searchsorted(sorted_keys, np.arange(n_cells + 1))

    # For every atom i, enumerate atoms j in its cell and 26 neighbor cells with j > i.
    pairs = []
    for i in range(n):
        gx, gy, gz = cell_xyz[i]
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    cx, cy, cz = gx + dx, gy + dy, gz + dz
                    if cx < 0 or cx >= nx or cy < 0 or cy >= ny or cz < 0 or cz >= nz:
                        continue
                    cell_key = cx + nx * cy + nx * ny * cz
                    lo = boundaries[cell_key]
                    hi = boundaries[cell_key + 1]
                    for slot in range(lo, hi):
                        j = order[slot]
                        if j <= i:
                            continue
                        diff = coords[i] - coords[j]
                        d2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]
                        if d2 <= cutoff * cutoff:
                            pairs.append((i, j))

    if not pairs:
        return np.zeros((0, 2), dtype=np.int32)
    return np.asarray(pairs, dtype=np.int32)
