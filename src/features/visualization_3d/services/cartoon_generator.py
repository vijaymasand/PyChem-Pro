"""
Cartoon Generator — Faithful port of CPPCartoon's createChainMesh.

Orchestrates the full pipeline:
  1. Extract PeptidePlanes from molecule (3-residue windows)
  2. Ensure consistent side-vector orientation
  3. Slide a 4-plane window and generate mesh segments
  4. Cache the result for interactive performance

Includes Adaptive LOD for large proteins.
"""

import numpy as np
from typing import Dict, Tuple, Optional

from src.features.visualization_3d.services.peptide_plane import compute_chain_planes
from src.features.visualization_3d.services.mesh_builder import (
    MeshAccumulator, create_segment_mesh, discontinuity
)

# ─── LOD Settings ────────────────────────────────────────────────────────────
# CPPCartoon defaults: splineSteps = 16, profileDetail = 8
# We scale based on protein size for performance on older hardware.

# LOD Settings - Optimized for professional smooth rendering
LOD_SETTINGS = {
    'ultra_low': {'spline_steps': 3, 'profile_detail': 4},      # ~85% reduction, smoother
    'low': {'spline_steps': 4, 'profile_detail': 6},           # ~80% reduction, smoother
    'medium': {'spline_steps': 6, 'profile_detail': 8},         # ~70% reduction, smoother
    'high': {'spline_steps': 8, 'profile_detail': 10},         # ~60% reduction, very smooth
    'ultra_high': {'spline_steps': 12, 'profile_detail': 16}    # ~40% reduction, ultra smooth
}


def _select_lod(n_atoms: int) -> dict:
    """Select LOD settings based on atom count for performance."""
    if n_atoms < 500:
        return LOD_SETTINGS['high']
    elif n_atoms < 1500:
        return LOD_SETTINGS['medium']
    elif n_atoms < 3000:
        return LOD_SETTINGS['low']
    else:
        return LOD_SETTINGS['ultra_low']


def generate_cartoon_mesh(molecule, spline_steps=None, profile_detail=None):
    """
    Generate the full protein cartoon mesh using Multiprocessing.
    Utilizes 50% of available CPU cores.
    """
    import os
    from concurrent.futures import ProcessPoolExecutor
    from src.features.visualization_3d.services.vectorized_mesh_builder import generate_chain_mesh_vectorized

    # Auto-select LOD if not specified
    if spline_steps is None or profile_detail is None:
        lod = _select_lod(len(molecule.atoms))
        spline_steps = spline_steps or lod['spline_steps']
        profile_detail = profile_detail or lod['profile_detail']
    
    # Performance logging
    import time
    start_time = time.perf_counter()

    # 1. Extract peptide planes per chain
    chains_data = compute_chain_planes(molecule)
    if not chains_data:
        return None, None, None

    # 2. Multiprocessing Orchestration for CPU-bound mesh generation
    # Use ProcessPoolExecutor for true parallelism on CPU-bound tasks
    # This avoids GIL limitations and provides better performance on multi-core systems
    from concurrent.futures import ProcessPoolExecutor
    
    # Sort chain IDs for consistent results
    sorted_chain_ids = sorted(chains_data.keys())
    chains_to_process = [chains_data[cid] for cid in sorted_chain_ids]

    results = []
    # Use ProcessPoolExecutor for CPU-bound mesh generation
    # Use 50% of CPU cores to avoid overwhelming the system
    num_workers = max(1, (os.cpu_count() or 4) // 2)
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(generate_chain_mesh_vectorized, chain, spline_steps, profile_detail)
            for chain in chains_to_process
        ]
        for future in futures:
            try:
                results.append(future.result())
            except Exception as e:
                import logging
                logging.error(f"Error processing chain: {e}")

    # 3. Merge Results
    final_vertices = []
    final_triangles = []
    final_colors = []
    
    for vertices, triangles, colors in results:
        if vertices is not None:
            final_vertices.append(vertices)
            final_triangles.append(triangles)
            final_colors.append(colors)
    
    if not final_vertices:
        return None, None, None
    
    # Concatenate all results
    all_vertices = np.concatenate(final_vertices, axis=0)
    all_triangles = np.concatenate(final_triangles, axis=0)
    all_colors = np.concatenate(final_colors, axis=0)
    
    # Performance timing
    end_time = time.perf_counter()
    mesh_time = end_time - start_time
    
    return all_vertices, all_triangles, all_colors


class CartoonGenerator:
    """
    Cached cartoon mesh generator.
    Generates mesh once and caches it for interactive rotation/zoom.
    """
    
    def __init__(self):
        self._cache = {}  # molecule id -> (vertices, triangles, colors)
    
    def get_mesh(self, molecule, spline_steps=None, profile_detail=None):
        """Retrieve cached mesh or generate a new one."""
        # Include LOD settings in cache key
        mol_id = (id(molecule), spline_steps, profile_detail)
        if mol_id in self._cache:
            return self._cache[mol_id]
        
        mesh = generate_cartoon_mesh(molecule, spline_steps=spline_steps, profile_detail=profile_detail)
        if mesh[0] is not None:
            self._cache[mol_id] = mesh
        
        return mesh
    
    def invalidate(self, molecule=None):
        """Clear cache for a specific molecule, or all."""
        if molecule is None:
            self._cache.clear()
        else:
            mol_id = id(molecule)
            self._cache.pop(mol_id, None)
