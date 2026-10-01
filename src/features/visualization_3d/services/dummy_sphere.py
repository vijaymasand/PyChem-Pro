"""
Dummy sphere creation for molecular visualization.

Allows users to create dummy spheres at user-defined positions
such as center of mass (COM), centroid, or custom coordinates.
"""

import numpy as np
from typing import Tuple, List, Optional
from src.core.domain.models.molecule import Molecule
from src.core.domain.models.atom import Atom


class DummySphere:
    """
    Represents a dummy sphere in 3D space for visualization purposes.
    """
    
    def __init__(self, position: Tuple[float, float, float], 
                 radius: float = 0.5, color: str = '#ffff00', 
                 label: str = "Dummy", sphere_id: str = None,
                 alpha: float = 1.0):
        """
        Initialize a dummy sphere.
        
        Args:
            position: (x, y, z) coordinates
            radius: Sphere radius
            color: Hex color string
            label: Display label for the sphere
            sphere_id: Unique identifier for the sphere
            alpha: Transparency (0.0 to 1.0, where 1.0 is opaque)
        """
        self.position = np.array(position)
        self.radius = radius
        self.color = color
        self.label = label
        self.sphere_id = sphere_id or f"sphere_{id(self)}"
        self.alpha = alpha
        self.visible = True
        
    def get_position(self) -> Tuple[float, float, float]:
        """Get sphere position as tuple."""
        return tuple(self.position)
    
    def set_position(self, position: Tuple[float, float, float]):
        """Set sphere position."""
        self.position = np.array(position)
    
    def get_distance_to(self, point: Tuple[float, float, float]) -> float:
        """Calculate distance from sphere center to a point."""
        point_array = np.array(point)
        return np.linalg.norm(self.position - point_array)


class DummySphereManager:
    """
    Manages creation and manipulation of dummy spheres.
    """
    
    def __init__(self, molecule: Molecule):
        self.molecule = molecule
        if not hasattr(molecule, 'dummy_spheres'):
            molecule.dummy_spheres = []
        self.spheres = molecule.dummy_spheres
        self.next_id = len(self.spheres) + 1
        
    def create_sphere_at_position(self, position: Tuple[float, float, float], 
                              radius: float = 0.5, 
                              color: str = '#ffff00',
                              label: str = "Dummy",
                              alpha: float = 1.0) -> str:
        """
        Create a dummy sphere at specified position.
        
        Args:
            position: (x, y, z) coordinates
            radius: Sphere radius
            color: Hex color string
            label: Display label
            alpha: Transparency (0.0 to 1.0)
            
        Returns:
            Sphere ID
        """
        sphere_id = f"dummy_{self.next_id}"
        sphere = DummySphere(position, radius, color, label, sphere_id, alpha=alpha)
        self.spheres.append(sphere)
        self.next_id += 1
        
        return sphere_id
    
    def create_sphere_at_com(self, radius: float = 0.5, 
                          color: str = '#ff00ff',
                          label: str = "COM",
                          alpha: float = 1.0) -> str:
        """
        Create a dummy sphere at the center of mass.
        
        Args:
            radius: Sphere radius
            color: Hex color string
            label: Display label
            alpha: Transparency (0.0 to 1.0)
            
        Returns:
            Sphere ID
        """
        com = self.calculate_center_of_mass()
        return self.create_sphere_at_position(com, radius, color, label, alpha=alpha)
    
    def create_sphere_at_centroid(self, radius: float = 0.5,
                              color: str = '#00ff00',
                              label: str = "Centroid") -> str:
        """
        Create a dummy sphere at the geometric centroid.
        
        Args:
            radius: Sphere radius
            color: Hex color string
            label: Display label
            
        Returns:
            Sphere ID
        """
        centroid = self.calculate_geometric_centroid()
        return self.create_sphere_at_position(centroid, radius, color, label)
    
    def create_sphere_at_atom_center(self, atom_indices: List[int],
                                 radius: float = 0.3,
                                 color: str = '#00ffff',
                                 label: str = "Atom Center") -> str:
        """
        Create a dummy sphere at the center of specified atoms.
        
        Args:
            atom_indices: List of atom indices
            radius: Sphere radius
            color: Hex color string
            label: Display label
            
        Returns:
            Sphere ID
        """
        center = self.calculate_atom_center(atom_indices)
        return self.create_sphere_at_position(center, radius, color, label)
    
    def calculate_center_of_mass(self) -> Tuple[float, float, float]:
        """
        Calculate the center of mass of the molecule, including H-atoms.
        
        Returns:
            (x, y, z) coordinates of COM
        """
        if not self.molecule or not self.molecule.atoms:
            return (0.0, 0.0, 0.0)
        
        total_mass = 0.0
        weighted_sum = np.zeros(3)
        
        for atom in self.molecule.atoms:
            if not atom.has_coords:
                continue
                
            # Inclusion of H-atoms is explicitly requested. 
            # Atom.mass returns the isotope mass or standard element mass.
            mass = atom.mass
            position = np.array([atom.x, atom.y, atom.z])
            
            weighted_sum += mass * position
            total_mass += mass
        
        if total_mass > 0:
            com = weighted_sum / total_mass
        else:
            com = np.zeros(3)
            
        return tuple(com)
    
    def calculate_geometric_centroid(self) -> Tuple[float, float, float]:
        """
        Calculate the geometric centroid of all atoms (including H-atoms).
        
        Returns:
            (x, y, z) coordinates of centroid
        """
        if not self.molecule or not self.molecule.atoms:
            return (0.0, 0.0, 0.0)
        
        positions = []
        for atom in self.molecule.atoms:
            if not atom.has_coords:
                continue
            positions.append([atom.x, atom.y, atom.z])
        
        if not positions:
            return (0.0, 0.0, 0.0)
            
        positions_array = np.array(positions)
        centroid = np.mean(positions_array, axis=0)
        
        return tuple(centroid)
    
    def calculate_atom_center(self, atom_indices: List[int]) -> Tuple[float, float, float]:
        """
        Calculate the geometric center of specified atoms.
        
        Args:
            atom_indices: List of atom indices
            
        Returns:
            (x, y, z) coordinates of atom center
        """
        if not atom_indices:
            return (0.0, 0.0, 0.0)
        
        positions = []
        for idx in atom_indices:
            if 0 <= idx < len(self.molecule.atoms):
                atom = self.molecule.atoms[idx]
                positions.append([atom.x, atom.y, atom.z])
        
        if not positions:
            return (0.0, 0.0, 0.0)
        
        positions_array = np.array(positions)
        center = np.mean(positions_array, axis=0)
        
        return tuple(center)
    
    def _get_atomic_mass(self, element_symbol: str) -> float:
        """
        Get atomic mass for an element.
        
        Args:
            element_symbol: Chemical element symbol
            
        Returns:
            Atomic mass in atomic mass units
        """
        # Simplified atomic mass table
        atomic_masses = {
            'H': 1.008,
            'He': 4.003,
            'Li': 6.941,
            'Be': 9.012,
            'B': 10.811,
            'C': 12.011,
            'N': 14.007,
            'O': 15.999,
            'F': 18.998,
            'Ne': 20.180,
            'Na': 22.990,
            'Mg': 24.305,
            'Al': 26.982,
            'Si': 28.086,
            'P': 30.974,
            'S': 32.065,
            'Cl': 35.453,
            'Ar': 39.948,
            'K': 39.098,
            'Ca': 40.078,
            'Sc': 44.956,
            'Ti': 47.867,
            'V': 50.942,
            'Cr': 51.996,
            'Mn': 54.938,
            'Fe': 55.845,
            'Co': 58.933,
            'Ni': 58.693,
            'Cu': 63.546,
            'Zn': 65.380,
            'Ga': 69.723,
            'Ge': 72.630,
            'As': 74.922,
            'Se': 78.971,
            'Br': 79.904,
            'Kr': 83.798,
            'Rb': 85.468,
            'Sr': 87.620,
            'Y': 88.906,
            'Zr': 91.224,
            'Nb': 92.906,
            'Mo': 95.950,
            'Tc': 98.000,
            'Ru': 101.070,
            'Rh': 102.906,
            'Pd': 106.420,
            'Ag': 107.868,
            'Cd': 112.411,
            'In': 114.818,
            'Sn': 118.710,
            'Sb': 121.760,
            'Te': 127.600,
            'I': 126.905,
            'Xe': 131.293,
        }
        
        return atomic_masses.get(element_symbol, 12.011)  # Default to carbon mass
    
    def get_sphere(self, sphere_id: str) -> Optional[DummySphere]:
        """
        Get a sphere by ID.
        
        Args:
            sphere_id: Sphere identifier
            
        Returns:
            DummySphere object or None if not found
        """
        for sphere in self.spheres:
            if sphere.sphere_id == sphere_id:
                return sphere
        return None
    
    def remove_sphere(self, sphere_id: str) -> bool:
        """
        Remove a sphere by ID.
        
        Args:
            sphere_id: Sphere identifier
            
        Returns:
            True if sphere was removed, False if not found
        """
        for i, sphere in enumerate(self.spheres):
            if sphere.sphere_id == sphere_id:
                del self.spheres[i]
                return True
        return False
    
    def clear_all_spheres(self):
        """Remove all dummy spheres."""
        self.spheres.clear()
        self.next_id = 1
    
    def get_all_spheres(self) -> List[DummySphere]:
        """Get list of all dummy spheres."""
        return self.spheres.copy()
    
    def get_sphere_summary(self) -> dict:
        """
        Get summary of all spheres.
        
        Returns:
            Dictionary with sphere information
        """
        return {
            'total_spheres': len(self.spheres),
            'spheres': [
                {
                    'id': sphere.sphere_id,
                    'position': sphere.get_position(),
                    'radius': sphere.radius,
                    'color': sphere.color,
                    'label': sphere.label,
                    'visible': sphere.visible,
                    'alpha': sphere.alpha
                }
                for sphere in self.spheres
            ]
        }


def create_dummy_sphere_at_com(molecule: Molecule, radius: float = 0.5) -> str:
    """
    Convenience function to create a dummy sphere at center of mass.
    
    Args:
        molecule: Molecule object
        radius: Sphere radius
        
    Returns:
        Sphere ID
    """
    manager = DummySphereManager(molecule)
    return manager.create_sphere_at_com(radius)


def create_dummy_sphere_at_position(molecule: Molecule, 
                               position: Tuple[float, float, float],
                               radius: float = 0.5) -> str:
    """
    Convenience function to create a dummy sphere at custom position.
    
    Args:
        molecule: Molecule object
        position: (x, y, z) coordinates
        radius: Sphere radius
        
    Returns:
        Sphere ID
    """
    manager = DummySphereManager(molecule)
    return manager.create_sphere_at_position(position, radius)
