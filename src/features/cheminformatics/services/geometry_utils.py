"""
Surface and Volume Utilities — Robust calculation of SASA and Molecular Volume.
Uses Shrake-Rupley algorithm (Fibonacci sphere) for surface and numerical integration for volume.
"""

import numpy as np
from typing import List, Tuple

def calculate_sasa(coords: np.ndarray, symbols: List[str], probe: float = 1.4, n_points: int = 256) -> float:
    """
    Calculate Solvent Accessible Surface Area (SASA).
    
    Args:
        coords: (N, 3) array of atom coordinates
        symbols: List of element symbols
        probe: Solvent probe radius (default 1.4 for water)
        n_points: Number of points per atom for the sphere
        
    Returns:
        float: SASA in Å²
    """
    if len(coords) == 0:
        return 0.0
        
    # vdW Radii (Bondi/Mantina values)
    VDW_RADII = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
        'P': 1.80, 'S': 1.80, 'CL': 1.75, 'BR': 1.85, 'I': 1.98,
        'SI': 2.10, 'B': 2.00, 'SE': 1.90, 'TE': 2.06
    }
    
    radii = np.array([VDW_RADII.get(s.upper(), 1.7) + probe for s in symbols])
    
    # Pre-generate unit sphere points (Fibonacci)
    indices = np.arange(0, n_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/n_points)
    theta = np.pi * (1 + 5**0.5) * indices
    unit_points = np.stack([
        np.cos(theta) * np.sin(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(phi)
    ], axis=1)
    
    total_area = 0.0
    for i in range(len(coords)):
        atom_pos = coords[i]
        radius = radii[i]
        
        # Scale unit sphere to SAS radius
        sphere_points = unit_points * radius + atom_pos
        
        # Check for overlaps with all other atoms
        # We only check atoms within a reasonable distance
        is_exposed = np.ones(n_points, dtype=bool)
        
        # Distances to all other atoms
        diffs = coords - atom_pos
        dists_sq = np.sum(diffs**2, axis=1)
        
        # Potential neighbors (within R_i + R_j + 2*probe)
        # But SASA logic: point on sphere i is buried if dist(point, atom_j) < R_j
        for j in range(len(coords)):
            if i == j: continue
            
            # Optimization: only check nearby atoms
            d_ij = np.sqrt(dists_sq[j])
            if d_ij > radius + radii[j]:
                continue
            
            # Distance from sphere points to atom j
            p_diffs = sphere_points[is_exposed] - coords[j]
            p_dists_sq = np.sum(p_diffs**2, axis=1)
            
            # Mask buried points
            is_exposed[is_exposed] &= (p_dists_sq >= radii[j]**2)
            
            if not np.any(is_exposed):
                break
                
        # Area contributed by this atom: (fraction exposed) * 4 * pi * R^2
        exposed_fraction = np.sum(is_exposed) / n_points
        total_area += exposed_fraction * 4 * np.pi * radius**2
        
    return total_area

def calculate_sasa_per_atom(coords: np.ndarray, symbols: List[str], probe: float = 1.4, n_points: int = 256) -> np.ndarray:
    """
    Calculate SASA for each atom individually (accounting for overlaps).
    Returns an array of areas for each atom.
    """
    if len(coords) == 0:
        return np.array([])
        
    VDW_RADII = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
        'P': 1.80, 'S': 1.80, 'CL': 1.75, 'BR': 1.85, 'I': 1.98
    }
    
    radii = np.array([VDW_RADII.get(s.upper(), 1.7) + probe for s in symbols])
    
    indices = np.arange(0, n_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/n_points)
    theta = np.pi * (1 + 5**0.5) * indices
    unit_points = np.stack([
        np.cos(theta) * np.sin(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(phi)
    ], axis=1)
    
    atom_areas = np.zeros(len(coords))
    for i in range(len(coords)):
        atom_pos = coords[i]
        radius = radii[i]
        sphere_points = unit_points * radius + atom_pos
        is_exposed = np.ones(n_points, dtype=bool)
        
        diffs = coords - atom_pos
        dists_sq = np.sum(diffs**2, axis=1)
        
        for j in range(len(coords)):
            if i == j: continue
            d_ij = np.sqrt(dists_sq[j])
            if d_ij > radius + radii[j]:
                continue
            
            p_diffs = sphere_points[is_exposed] - coords[j]
            p_dists_sq = np.sum(p_diffs**2, axis=1)
            is_exposed[is_exposed] &= (p_dists_sq >= radii[j]**2)
            
            if not np.any(is_exposed):
                break
                
        exposed_fraction = np.sum(is_exposed) / n_points
        atom_areas[i] = exposed_fraction * 4 * np.pi * radius**2
        
    return atom_areas

def calculate_volume(coords: np.ndarray, symbols: List[str], n_points: int = 100) -> float:
    """
    Calculate Van der Waals volume using a simple but effective overlap-aware method.
    Estimates volume by integration.
    """
    if len(coords) == 0:
        return 0.0
        
    VDW_RADII = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
        'P': 1.80, 'S': 1.80, 'CL': 1.75, 'BR': 1.85, 'I': 1.98
    }
    
    radii = np.array([VDW_RADII.get(s.upper(), 1.7) for s in symbols])
    
    # Simple estimation for now: sum of spheres corrected by overlap factor
    # A better way is a grid or analytical formula (but those are complex)
    # We'll use a Monte Carlo style estimation in a bounding box for accuracy
    
    min_coords = np.min(coords - radii[:, np.newaxis], axis=0)
    max_coords = np.max(coords + radii[:, np.newaxis], axis=0)
    
    # Bounding box volume
    bbox_vol = np.prod(max_coords - min_coords)
    
    # Sample points in bbox
    n_samples = 10000
    samples = np.random.uniform(min_coords, max_coords, (n_samples, 3))
    
    # Check how many samples are inside any VdW sphere
    inside_count = 0
    for i in range(len(samples)):
        p = samples[i]
        dists_sq = np.sum((coords - p)**2, axis=1)
        if np.any(dists_sq <= radii**2):
            inside_count += 1
            
    return bbox_vol * (inside_count / n_samples)
