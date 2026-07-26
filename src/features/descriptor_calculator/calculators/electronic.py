"""
Electronic descriptor calculations.
Includes: partial charges, dipole moment, polar surface area.
"""
import numpy as np
from . import BaseCalculator


class ElectronicCalculator(BaseCalculator):
    """Calculator for electronic molecular descriptors."""

    def _atom_charges(self, molecule, selection):
        """Partial charges of the selected atoms, computed if not present yet."""
        self.ensure_charges(molecule)
        return [getattr(molecule.atoms[idx], 'partial_charge', 0.0) or 0.0
                for idx in self.get_selected_atoms(molecule, selection)]

    def calc_total_charge(self, molecule, selection) -> float:
        """Calculate total molecular charge."""
        return float(sum(self._atom_charges(molecule, selection)))

    def calc_max_partial_charge(self, molecule, selection) -> float:
        """Calculate maximum partial charge."""
        return max(self._atom_charges(molecule, selection), default=0.0)

    def calc_min_partial_charge(self, molecule, selection) -> float:
        """Calculate minimum partial charge."""
        return min(self._atom_charges(molecule, selection), default=0.0)

    def calc_max_positive_charge(self, molecule, selection) -> float:
        """Calculate maximum positive partial charge."""
        charges = [c for c in self._atom_charges(molecule, selection) if c > 0]
        return max(charges) if charges else 0.0

    def calc_max_negative_charge(self, molecule, selection) -> float:
        """Calculate maximum negative partial charge (most negative)."""
        charges = [c for c in self._atom_charges(molecule, selection) if c < 0]
        return min(charges) if charges else 0.0

    def calc_dipole_moment(self, molecule, selection) -> float:
        """Calculate molecular dipole moment."""
        self.ensure_coordinates(molecule)
        self.ensure_charges(molecule)
        dipole = np.array([0.0, 0.0, 0.0])

        for idx in self.get_selected_atoms(molecule, selection):
            atom = molecule.atoms[idx]
            if atom.x is None or atom.y is None or atom.z is None:
                continue
            charge = getattr(atom, 'partial_charge', 0.0) or 0.0
            dipole += charge * np.array([atom.x, atom.y, atom.z])

        return float(np.linalg.norm(dipole) * 2.54)  # Convert to Debye

    def calc_polar_surface_area(self, molecule, selection) -> float:
        """Calculate polar surface area (PSA) using real SASA contributions from polar atoms."""
        from src.features.cheminformatics.services.geometry_utils import calculate_sasa
        from ...cheminformatics.services.atom_properties import AtomPropertyAnalyzer
        
        self.ensure_coordinates(molecule)
        analyzer = AtomPropertyAnalyzer(molecule)
        polar_symbols = set(analyzer.POLAR_ATOMS)

        atoms = [molecule.atoms[idx] for idx in selection.atom_indices if idx < len(molecule.atoms)]
        if not atoms or not all(a.has_coords for a in atoms):
            return 0.0
            
        coords = np.array([[a.x, a.y, a.z] for a in atoms])
        symbols = [a.symbol for a in atoms]
        
        # Calculate SASA and only sum contributions from polar atoms
        # I'll update calculate_sasa to return per-atom contributions
        from src.features.cheminformatics.services.geometry_utils import calculate_sasa_per_atom
        atom_sasa = calculate_sasa_per_atom(coords, symbols)
        
        total_psa = 0.0
        for i, symbol in enumerate(symbols):
            if symbol.upper() in polar_symbols:
                total_psa += atom_sasa[i]
                
        return total_psa

    def calc_apolar_surface_area(self, molecule, selection) -> float:
        """Calculate apolar (nonpolar) surface area."""
        from src.features.cheminformatics.services.geometry_utils import calculate_sasa
        from ...cheminformatics.services.atom_properties import AtomPropertyAnalyzer

        self.ensure_coordinates(molecule)
        atoms = [molecule.atoms[idx] for idx in selection.atom_indices if idx < len(molecule.atoms)]
        if not atoms or not all(a.has_coords for a in atoms):
            return 0.0
            
        coords = np.array([[a.x, a.y, a.z] for a in atoms])
        symbols = [a.symbol for a in atoms]
        
        # 1. Total SASA
        total_sasa = calculate_sasa(coords, symbols)
        
        # 2. Polar SASA (sum of SASA contributions from polar atoms)
        analyzer = AtomPropertyAnalyzer(molecule)
        polar_indices = set()
        for i, atom in enumerate(molecule.atoms):
             if atom.symbol in analyzer.POLAR_ATOMS:
                 polar_indices.add(i)
        
        return total_sasa - self.calc_polar_surface_area(molecule, selection)

    def calc_mean_absolute_charge(self, molecule, selection) -> float:
        """Calculate mean absolute partial charge."""
        charges = self._atom_charges(molecule, selection)
        if not charges:
            return 0.0
        return float(sum(abs(c) for c in charges) / len(charges))
