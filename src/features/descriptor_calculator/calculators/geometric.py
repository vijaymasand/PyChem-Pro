"""
Geometric descriptor calculations.
Includes: surface areas, volumes, radius of gyration, shape descriptors, principal moments.
"""
import numpy as np
from . import BaseCalculator


class GeometricCalculator(BaseCalculator):
    """Calculator for geometric (3D) molecular descriptors."""

    def calc_sasa(self, molecule, selection) -> float:
        """Calculate solvent accessible surface area (SASA) using Shrake-Rupley algorithm."""
        from src.features.cheminformatics.services.geometry_utils import calculate_sasa

        self.ensure_coordinates(molecule)
        atoms = [molecule.atoms[idx] for idx in selection.atom_indices if idx < len(molecule.atoms)]
        if not atoms or not all(a.has_coords for a in atoms):
            return 0.0
            
        coords = np.array([[a.x, a.y, a.z] for a in atoms])
        symbols = [a.symbol for a in atoms]
        
        return calculate_sasa(coords, symbols)

    def calc_molecular_volume(self, molecule, selection) -> float:
        """Calculate molecular volume using numerical integration."""
        from src.features.cheminformatics.services.geometry_utils import calculate_volume

        self.ensure_coordinates(molecule)
        atoms = [molecule.atoms[idx] for idx in selection.atom_indices if idx < len(molecule.atoms)]
        if not atoms or not all(a.has_coords for a in atoms):
            return 0.0
            
        coords = np.array([[a.x, a.y, a.z] for a in atoms])
        symbols = [a.symbol for a in atoms]
        
        return calculate_volume(coords, symbols)

    def calc_radius_of_gyration(self, molecule, selection) -> float:
        """Calculate radius of gyration."""
        self.ensure_coordinates(molecule)
        coords = []
        for idx in selection.atom_indices:
            if idx < len(molecule.atoms):
                atom = molecule.atoms[idx]
                if atom.x is None or atom.y is None or atom.z is None:
                    continue
                coords.append([atom.x, atom.y, atom.z])

        if not coords:
            return 0.0

        coords = np.array(coords)
        center = np.mean(coords, axis=0)
        distances = np.sum((coords - center) ** 2, axis=1)
        return np.sqrt(np.mean(distances))

    def calc_asphericity(self, molecule, selection) -> float:
        """Calculate molecular asphericity."""
        eigenvalues = self.get_gyration_eigenvalues(molecule, selection)
        if len(eigenvalues) < 3:
            return 0.0

        l1, l2, l3 = sorted(eigenvalues, reverse=True)
        sum_l = l1 + l2 + l3
        if sum_l == 0:
            return 0.0
        return 1.0 - 3.0 * (l1*l2 + l2*l3 + l1*l3) / (sum_l ** 2)

    def calc_eccentricity(self, molecule, selection) -> float:
        """Calculate molecular eccentricity."""
        eigenvalues = self.get_gyration_eigenvalues(molecule, selection)
        if len(eigenvalues) < 3:
            return 0.0

        l1, l2, l3 = sorted(eigenvalues, reverse=True)
        if l1 > 0:
            return np.sqrt(l1**2 - l3**2) / l1
        return 0.0

    def calc_principal_moment_1(self, molecule, selection) -> float:
        """Calculate first principal moment of inertia."""
        eigenvalues = self.get_gyration_eigenvalues(molecule, selection)
        return max(eigenvalues) if len(eigenvalues) > 0 else 0.0

    def calc_principal_moment_2(self, molecule, selection) -> float:
        """Calculate second principal moment of inertia."""
        eigenvalues = self.get_gyration_eigenvalues(molecule, selection)
        if len(eigenvalues) >= 2:
            return sorted(eigenvalues, reverse=True)[1]
        return 0.0

    def calc_principal_moment_3(self, molecule, selection) -> float:
        """Calculate third principal moment of inertia."""
        eigenvalues = self.get_gyration_eigenvalues(molecule, selection)
        if len(eigenvalues) >= 3:
            return sorted(eigenvalues, reverse=True)[2]
        return 0.0

    def calc_molecular_diameter(self, molecule, selection) -> float:
        """Calculate molecular diameter (max distance between atoms)."""
        coords = []
        for idx in selection.atom_indices:
            if idx < len(molecule.atoms):
                atom = molecule.atoms[idx]
                if atom.x is not None and atom.y is not None and atom.z is not None:
                    coords.append(np.array([atom.x, atom.y, atom.z]))

        if len(coords) < 2:
            return 0.0

        max_dist = 0.0
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                dist = np.linalg.norm(coords[i] - coords[j])
                max_dist = max(max_dist, dist)

        return max_dist
