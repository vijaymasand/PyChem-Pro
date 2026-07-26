"""
Base calculator module with shared utilities for descriptor calculations.
"""
from typing import List, Set, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from ...core.domain.models.molecule import Molecule
    from ..descriptor_types import AtomSelection


class BaseCalculator:
    """Base class for descriptor calculators with shared utilities."""

    def __init__(self):
        pass

    def ensure_perception(self, molecule: "Molecule"):
        """
        Ensure the molecule has aromaticity perceived and hybridization assigned.
        This is critical for accurate constitutional and topological descriptors.
        """
        from ...smiles_parser.rules.aromaticity import perceive_aromaticity
        
        # Always run perception if it hasn't been run or to ensure sync
        # Molecule.find_rings() is already called inside perceive_aromaticity
        perceive_aromaticity(molecule)
        molecule.propagate_aromaticity()
        molecule.assign_hybridization()

    def get_selected_atoms(self, molecule: "Molecule", selection: "AtomSelection") -> List[int]:
        """Get list of selected atom indices that are valid."""
        return [idx for idx in selection.atom_indices if idx < len(molecule.atoms)]

    def get_selected_set(self, selection: "AtomSelection") -> Set[int]:
        """Get set of selected atom indices."""
        return set(selection.atom_indices)

    def count_atoms_by_symbol(self, molecule: "Molecule", selection: "AtomSelection", symbol: str) -> int:
        """Count atoms with given symbol in selection."""
        return sum(1 for idx in selection.atom_indices
                  if idx < len(molecule.atoms) and molecule.atoms[idx].symbol == symbol)

    def get_distance_matrix(self, molecule: "Molecule") -> List[List[float]]:
        """Get distance matrix for molecule using Floyd-Warshall algorithm."""
        n_atoms = len(molecule.atoms)
        distance_matrix = [[float('inf')] * n_atoms for _ in range(n_atoms)]

        # Initialize with direct bonds
        for i in range(n_atoms):
            distance_matrix[i][i] = 0.0
            for neighbor_idx in molecule.get_neighbors(i):
                distance_matrix[i][neighbor_idx] = 1.0

        # Floyd-Warshall algorithm for shortest paths
        for k in range(n_atoms):
            for i in range(n_atoms):
                for j in range(n_atoms):
                    if distance_matrix[i][k] + distance_matrix[k][j] < distance_matrix[i][j]:
                        distance_matrix[i][j] = distance_matrix[i][k] + distance_matrix[k][j]

        return distance_matrix

    def ensure_coordinates(self, molecule: "Molecule") -> bool:
        """Make sure the molecule has 3D coordinates, generating them if needed.

        A structure parsed from SMILES carries no geometry, which would leave
        every shape, surface and charge-surface descriptor at zero. The project's
        own generator is used, the same way the LogP batch tool does it.
        Returns True when coordinates are available afterwards. """
        atoms = getattr(molecule, 'atoms', None)
        if not atoms:
            return False
        if all(a.x is not None for a in atoms):
            return True
        if getattr(molecule, '_descriptor_coordgen_failed', False):
            return False
        try:
            from src.features.layout_3d import generate_3d_coordinates
            generate_3d_coordinates(molecule, optimize=True, max_opt_steps=100)
        except Exception:
            pass
        if all(a.x is not None for a in atoms):
            return True
        molecule._descriptor_coordgen_failed = True
        return False

    def ensure_charges(self, molecule: "Molecule"):
        """Compute Gasteiger partial charges once when the molecule has none.

        Without this every charge based descriptor silently reports zero for
        structures that were never passed through a charge calculation. """
        atoms = getattr(molecule, 'atoms', None)
        if not atoms:
            return
        if any(abs(getattr(a, 'partial_charge', 0.0) or 0.0) > 1e-9 for a in atoms):
            return
        if getattr(molecule, '_descriptor_charges_failed', False):
            return
        try:
            from ...cheminformatics.electrostatics.gasteiger import compute_gasteiger_charges
            compute_gasteiger_charges(molecule)
        except Exception:
            molecule._descriptor_charges_failed = True

    def get_gyration_eigenvalues(self, molecule: "Molecule", selection: "AtomSelection") -> np.ndarray:
        """Helper to compute eigenvalues of the gyration tensor."""
        self.ensure_coordinates(molecule)
        coords = []
        for idx in selection.atom_indices:
            if idx < len(molecule.atoms):
                atom = molecule.atoms[idx]
                if atom.x is None or atom.y is None or atom.z is None:
                    continue
                coords.append([atom.x, atom.y, atom.z])

        if not coords or len(coords) < 2:
            return np.array([0.0, 0.0, 0.0])

        coords = np.array(coords)
        center = np.mean(coords, axis=0)
        centered_coords = coords - center

        # Gyration tensor G = (1/N) * X^T * X
        gyration_tensor = np.dot(centered_coords.T, centered_coords) / len(coords)

        try:
            eigenvalues = np.linalg.eigvals(gyration_tensor)
            return np.real(eigenvalues)
        except np.linalg.LinAlgError:
            return np.array([0.0, 0.0, 0.0])

    def get_ring_sizes(self, molecule: "Molecule", selection: "AtomSelection") -> List[int]:
        """Helper to get ring sizes in selection."""
        rings = molecule.find_rings()
        selected_set = set(selection.atom_indices)
        sizes = []
        for ring in rings:
            if all(idx in selected_set for idx in ring):
                sizes.append(len(ring))
        return sizes

    def get_bond_count_in_selection(self, molecule: "Molecule", selected_set: Set[int]) -> int:
        """Count bonds where both atoms are in selection."""
        return sum(1 for bond in molecule.bonds
                  if bond.begin_atom_idx in selected_set and bond.end_atom_idx in selected_set)


# Export all calculators. The Extended* classes subclass the basic ones and add
# the descriptor families documented in MOLECULAR_DESCRIPTORS.md; the engine
# instantiates those, so a single object per category serves both sets.
from .constitutional import ConstitutionalCalculator
from .topological import TopologicalCalculator
from .geometric import GeometricCalculator
from .electronic import ElectronicCalculator
from .quantum import QuantumCalculator
from .fingerprints import FingerprintCalculator
from .hybrid import HybridCalculator
from .constitutional_ext import ExtendedConstitutionalCalculator
from .topological_ext import ExtendedTopologicalCalculator
from .geometric_ext import ExtendedGeometricCalculator
from .electronic_ext import ExtendedElectronicCalculator
from .quantum_ext import ExtendedQuantumCalculator
from .hybrid_ext import ExtendedHybridCalculator

__all__ = [
    'BaseCalculator',
    'ConstitutionalCalculator',
    'TopologicalCalculator',
    'GeometricCalculator',
    'ElectronicCalculator',
    'QuantumCalculator',
    'FingerprintCalculator',
    'HybridCalculator',
    'ExtendedConstitutionalCalculator',
    'ExtendedTopologicalCalculator',
    'ExtendedGeometricCalculator',
    'ExtendedElectronicCalculator',
    'ExtendedQuantumCalculator',
    'ExtendedHybridCalculator',
]
