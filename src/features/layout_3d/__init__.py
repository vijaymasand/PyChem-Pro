"""
Native 3D Coordinate Generator — Pure Python implementation.

Generates 3D coordinates for molecules from scratch without RDKit.
Uses geometric and force-field inspired methods.
"""

import math
import random
from typing import List, Tuple, Optional
from src.core.domain.models.molecule import Molecule
from src.core.domain.models.atom import Atom
from src.core.domain.models.bond import BondType


# Covalent radii (Angstroms) for common elements
COVALENT_RADII = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39,
    'Si': 1.11, 'B': 0.84, 'Na': 1.54, 'K': 1.96, 'Ca': 1.74,
    'Mg': 1.36, 'Fe': 1.32, 'Zn': 1.22, 'Cu': 1.32, 'Au': 1.36,
    'Ag': 1.45, 'Li': 1.28, 'Al': 1.21, 'Se': 1.20, 'As': 1.21
}

# Default radius for unknown elements
DEFAULT_RADIUS = 1.0

# Ideal bond lengths (Angstroms) based on atom pairs
IDEAL_BOND_LENGTHS = {
    ('C', 'C'): 1.54,  # Single
    ('C', 'N'): 1.47,
    ('C', 'O'): 1.43,
    ('C', 'H'): 1.10,
    ('C', 'S'): 1.82,
    ('C', 'P'): 1.85,
    ('C', 'F'): 1.35,
    ('C', 'Cl'): 1.77,
    ('C', 'Br'): 1.94,
    ('C', 'I'): 2.14,
    ('N', 'N'): 1.45,
    ('N', 'O'): 1.40,
    ('N', 'H'): 1.03,
    ('O', 'O'): 1.48,
    ('O', 'H'): 0.96,
    ('O', 'S'): 1.70,
    ('O', 'P'): 1.63,
    ('S', 'S'): 2.05,
    ('S', 'H'): 1.34,
    ('P', 'H'): 1.42,
    ('P', 'P'): 2.21,
    ('F', 'H'): 0.92,
    ('Cl', 'H'): 1.27,
    ('Br', 'H'): 1.41,
    ('I', 'H'): 1.61,
}


def get_covalent_radius(symbol: str) -> float:
    """Get covalent radius for an element."""
    return COVALENT_RADII.get(symbol, DEFAULT_RADIUS)


def get_ideal_bond_length(atom1: str, atom2: str) -> float:
    """Get ideal bond length between two atoms."""
    key = (atom1, atom2)
    reverse_key = (atom2, atom1)
    
    if key in IDEAL_BOND_LENGTHS:
        return IDEAL_BOND_LENGTHS[key]
    elif reverse_key in IDEAL_BOND_LENGTHS:
        return IDEAL_BOND_LENGTHS[reverse_key]
    else:
        # Sum of covalent radii as fallback
        return get_covalent_radius(atom1) + get_covalent_radius(atom2)


class CoordinateGenerator3D:
    """
    Generates 3D coordinates from molecular topology.
    
    Algorithm:
    1. Start with first atom at origin
    2. Place bonded atoms using bond vectors
    3. Use geometric construction for branches
    4. Simple force-field relaxation
    """
    
    def __init__(self, molecule: Molecule):
        self.molecule = molecule
        self.coords = {}  # atom_idx -> (x, y, z)
        
    def generate(self, optimize: bool = True, max_steps: int = 100) -> dict:
        """
        Generate 3D coordinates for all atoms.
        
        Args:
            optimize: Whether to perform geometry optimization
            max_steps: Maximum optimization steps
            
        Returns:
            Dictionary mapping atom indices to (x, y, z) coordinates
        """
        if not self.molecule.atoms:
            return {}
            
        # Build adjacency list
        adjacency = self._build_adjacency()
        
        # Place atoms using breadth-first placement
        self._place_atoms_bfs(adjacency)
        
        # Optimize geometry if requested
        if optimize and len(self.molecule.atoms) > 1:
            self._optimize_geometry(max_steps)
        
        # Assign coordinates to atoms
        for idx, (x, y, z) in self.coords.items():
            self.molecule.atoms[idx].x = x
            self.molecule.atoms[idx].y = y
            self.molecule.atoms[idx].z = z
        
        return self.coords
    
    def _build_adjacency(self) -> dict:
        """Build adjacency list for the molecule."""
        adjacency = {i: [] for i in range(len(self.molecule.atoms))}
        for bond in self.molecule.bonds:
            adjacency[bond.begin_atom_idx].append(bond.end_atom_idx)
            adjacency[bond.end_atom_idx].append(bond.begin_atom_idx)
        return adjacency
    
    def _place_atoms_bfs(self, adjacency: dict):
        """Place atoms using breadth-first search from first atom."""
        if not self.molecule.atoms:
            return
            
        # Start with first atom at origin
        self.coords[0] = (0.0, 0.0, 0.0)
        placed = {0}
        queue = [0]
        
        while queue:
            current = queue.pop(0)
            current_atom = self.molecule.atoms[current]
            
            for neighbor in adjacency[current]:
                if neighbor in placed:
                    continue
                    
                neighbor_atom = self.molecule.atoms[neighbor]
                
                # Get ideal bond length
                bond_length = get_ideal_bond_length(
                    current_atom.symbol, 
                    neighbor_atom.symbol
                )
                
                # Calculate position for neighbor
                new_pos = self._calculate_atom_position(
                    current, neighbor, bond_length, adjacency
                )
                
                self.coords[neighbor] = new_pos
                placed.add(neighbor)
                queue.append(neighbor)
    
    def _calculate_atom_position(self, parent_idx: int, atom_idx: int, 
                                  bond_length: float, adjacency: dict) -> Tuple[float, float, float]:
        """Calculate 3D position for a new atom based on its parent."""
        px, py, pz = self.coords[parent_idx]
        parent_atom = self.molecule.atoms[parent_idx]
        atom = self.molecule.atoms[atom_idx]
        
        # Count already placed neighbors of parent (excluding current atom)
        placed_neighbors = [n for n in adjacency[parent_idx] 
                          if n in self.coords and n != atom_idx]
        
        num_placed = len(placed_neighbors)
        
        # Determine geometry based on number of neighbors
        if num_placed == 0:
            # First neighbor - place along +X axis
            return (px + bond_length, py, pz)
        
        elif num_placed == 1:
            # Second neighbor - place in XY plane with ~109.5° angle
            neighbor_idx = placed_neighbors[0]
            nx, ny, nz = self.coords[neighbor_idx]
            
            # Vector from parent to first neighbor
            v1x, v1y, v1z = nx - px, ny - py, nz - pz
            v1_len = math.sqrt(v1x**2 + v1y**2 + v1z**2)
            
            if v1_len > 0:
                # Normalize
                v1x, v1y, v1z = v1x/v1_len, v1y/v1_len, v1z/v1_len
                
                # Create perpendicular vector for tetrahedral geometry
                # Tetrahedral angle is ~109.5°, so angle between bonds is ~109.5°
                # cos(109.5°) ≈ -1/3
                
                # Find a perpendicular direction
                if abs(v1x) < 0.9:
                    perp_x, perp_y, perp_z = 1.0, 0.0, 0.0
                else:
                    perp_x, perp_y, perp_z = 0.0, 1.0, 0.0
                
                # Make perpendicular to v1
                dot = v1x * perp_x + v1y * perp_y + v1z * perp_z
                perp_x -= dot * v1x
                perp_y -= dot * v1y
                perp_z -= dot * v1z
                
                # Normalize perpendicular
                perp_len = math.sqrt(perp_x**2 + perp_y**2 + perp_z**2)
                if perp_len > 0:
                    perp_x, perp_y, perp_z = perp_x/perp_len, perp_y/perp_len, perp_z/perp_len
                
                # Tetrahedral angle: new bond at ~109.5° from existing bond
                # cos(109.5°) ≈ -1/3
                cos_theta = -1.0/3.0
                sin_theta = math.sqrt(1 - cos_theta**2)
                
                # New direction: mix of -v1 and perpendicular component
                new_x = cos_theta * (-v1x) + sin_theta * perp_x
                new_y = cos_theta * (-v1y) + sin_theta * perp_y
                new_z = cos_theta * (-v1z) + sin_theta * perp_z
                
                # Scale by bond length
                return (px + bond_length * new_x, 
                       py + bond_length * new_y, 
                       pz + bond_length * new_z)
            
            # Fallback: place at random offset
            return (px + bond_length, py, pz)
        
        else:
            # Third or more neighbors - use random directions with repulsion
            best_pos = None
            best_score = float('inf')
            
            for _ in range(50):  # Try multiple random positions
                # Random direction
                theta = random.uniform(0, 2 * math.pi)
                phi = random.uniform(0, math.pi)
                
                dx = math.sin(phi) * math.cos(theta)
                dy = math.sin(phi) * math.sin(theta)
                dz = math.cos(phi)
                
                pos = (px + bond_length * dx, 
                      py + bond_length * dy, 
                      pz + bond_length * dz)
                
                # Score based on repulsion from existing neighbors
                score = 0
                for neighbor_idx in placed_neighbors:
                    nx, ny, nz = self.coords[neighbor_idx]
                    dist_sq = (pos[0] - nx)**2 + (pos[1] - ny)**2 + (pos[2] - nz)**2
                    score += 1.0 / (dist_sq + 0.1)  # Avoid very close positions
                
                if score < best_score:
                    best_score = score
                    best_pos = pos
            
            return best_pos if best_pos else (px + bond_length, py, pz)
    
    def _optimize_geometry(self, max_steps: int = 100):
        """Simple force-field inspired geometry optimization."""
        if len(self.coords) < 2:
            return
            
        # Get bond information
        bonds_info = []
        for bond in self.molecule.bonds:
            atom1 = self.molecule.atoms[bond.begin_atom_idx]
            atom2 = self.molecule.atoms[bond.end_atom_idx]
            ideal_length = get_ideal_bond_length(atom1.symbol, atom2.symbol)
            bonds_info.append((bond.begin_atom_idx, bond.end_atom_idx, ideal_length))
        
        # Simple gradient descent
        for step in range(max_steps):
            forces = {idx: [0.0, 0.0, 0.0] for idx in self.coords}
            
            # Bond stretching forces
            for i, j, ideal_len in bonds_info:
                if i not in self.coords or j not in self.coords:
                    continue
                    
                xi, yi, zi = self.coords[i]
                xj, yj, zj = self.coords[j]
                
                dx = xj - xi
                dy = yj - yi
                dz = zj - zi
                
                dist = math.sqrt(dx**2 + dy**2 + dz**2)
                if dist < 0.01:
                    dist = 0.01
                
                # Force proportional to deviation from ideal length
                force_mag = 0.1 * (dist - ideal_len) / dist
                
                fx = force_mag * dx
                fy = force_mag * dy
                fz = force_mag * dz
                
                forces[i][0] += fx
                forces[i][1] += fy
                forces[i][2] += fz
                forces[j][0] -= fx
                forces[j][1] -= fy
                forces[j][2] -= fz
            
            # Non-bonded repulsion (simple)
            atom_indices = list(self.coords.keys())
            for i in range(len(atom_indices)):
                for j in range(i + 1, len(atom_indices)):
                    idx1, idx2 = atom_indices[i], atom_indices[j]
                    
                    # Check if bonded
                    is_bonded = False
                    for bi, bj, _ in bonds_info:
                        if (bi == idx1 and bj == idx2) or (bi == idx2 and bj == idx1):
                            is_bonded = True
                            break
                    
                    if is_bonded:
                        continue
                    
                    x1, y1, z1 = self.coords[idx1]
                    x2, y2, z2 = self.coords[idx2]
                    
                    dx = x2 - x1
                    dy = y2 - y1
                    dz = z2 - z1
                    
                    dist_sq = dx**2 + dy**2 + dz**2
                    dist = math.sqrt(dist_sq)
                    
                    # Van der Waals repulsion at short distances
                    vdw_sum = (get_covalent_radius(self.molecule.atoms[idx1].symbol) + 
                              get_covalent_radius(self.molecule.atoms[idx2].symbol))
                    
                    if dist < 2.0 * vdw_sum:
                        force_mag = 0.5 * (2.0 * vdw_sum - dist) / dist if dist > 0.01 else 0.5
                        
                        fx = force_mag * dx
                        fy = force_mag * dy
                        fz = force_mag * dz
                        
                        forces[idx1][0] -= fx
                        forces[idx1][1] -= fy
                        forces[idx1][2] -= fz
                        forces[idx2][0] += fx
                        forces[idx2][1] += fy
                        forces[idx2][2] += fz
            
            # Update positions
            max_movement = 0.0
            for idx in self.coords:
                fx, fy, fz = forces[idx]
                
                # Limit step size
                step_size = math.sqrt(fx**2 + fy**2 + fz**2)
                if step_size > 0.1:
                    scale = 0.1 / step_size
                    fx *= scale
                    fy *= scale
                    fz *= scale
                
                x, y, z = self.coords[idx]
                new_x = x + fx
                new_y = y + fy
                new_z = z + fz
                
                movement = math.sqrt(fx**2 + fy**2 + fz**2)
                max_movement = max(max_movement, movement)
                
                self.coords[idx] = (new_x, new_y, new_z)
            
            # Convergence check
            if max_movement < 0.001:
                break


def generate_3d_coordinates(molecule: Molecule, optimize: bool = True, 
                           max_opt_steps: int = 100) -> dict:
    """
    Generate 3D coordinates for a molecule from scratch.
    
    Args:
        molecule: Molecule object
        optimize: Whether to perform geometry optimization
        max_opt_steps: Maximum optimization steps
        
    Returns:
        Dictionary mapping atom indices to (x, y, z) coordinates
    """
    generator = CoordinateGenerator3D(molecule)
    return generator.generate(optimize=optimize, max_steps=max_opt_steps)
