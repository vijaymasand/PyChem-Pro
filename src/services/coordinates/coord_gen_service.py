"""
Coordinate generation service with parallel conformer search.
Implements ICoordinateGenerator protocol.
"""
import random
from src.core.parallel import ParallelExecutor
from src.core.domain.models.molecule import Molecule


# Module-level worker (picklable for spawn mode)
def _generate_single_conformer(args):
    """Generate one conformer with a random seed. Returns (coords_dict, energy)."""
    mol_data, seed, max_steps = args
    random.seed(seed)

    # Reconstruct molecule from serialized data
    from src.features.layout_3d import CoordinateGenerator3D, generate_3d_coordinates
    from src.core.domain.models.molecule import Molecule
    from src.core.domain.models.atom import Atom
    from src.core.domain.models.bond import BondType

    # Rebuild molecule (can't pickle Molecule directly across processes)
    mol = Molecule("conformer")
    for symbol, charge, impl_h in mol_data['atoms']:
        a = Atom(symbol)
        a.formal_charge = charge
        a.num_implicit_h = impl_h
        mol.add_atom(a)
    for begin, end, order in mol_data['bonds']:
        bt = {1: BondType.SINGLE, 2: BondType.DOUBLE, 3: BondType.TRIPLE}.get(
            order, BondType.SINGLE
        )
        mol.add_bond(begin, end, bt)

    # Generate coordinates
    generate_3d_coordinates(mol, optimize=True, max_opt_steps=max_steps)

    # Compute energy
    try:
        from src.services.forcefield.mmff94_service import MMFF94Service
        svc = MMFF94Service()
        energy = svc.compute_energy(mol)
    except Exception:
        energy = float('inf')

    # Extract coords
    coords = {}
    for i, atom in enumerate(mol.atoms):
        if atom.has_coords:
            coords[i] = (atom.x, atom.y, atom.z)

    return coords, energy


class CoordinateGeneratorService:
    """Coordinate generation with parallel conformer search."""

    def __init__(self, executor: ParallelExecutor = None):
        self.executor = executor or ParallelExecutor()

    def generate_3d(self, mol: Molecule, optimize: bool = True,
                    max_steps: int = 200) -> None:
        """Generate 3D coordinates in-place."""
        from src.features.layout_3d import generate_3d_coordinates
        generate_3d_coordinates(mol, optimize=optimize, max_opt_steps=max_steps)

    def generate_conformers(self, mol: Molecule, n: int = 10,
                            max_steps: int = 100) -> list:
        """
        Generate N conformers in parallel, return list of (Molecule_clone, energy).
        The first result is the lowest energy conformer.
        """
        # Serialize molecule data for cross-process transfer
        mol_data = self._serialize_mol(mol)

        # Create tasks with different random seeds
        tasks = [
            (mol_data, random.randint(0, 2**31), max_steps)
            for _ in range(n)
        ]

        # Run in parallel
        results = self.executor.map(_generate_single_conformer, tasks)

        # Sort by energy, apply best coords to original molecule
        results.sort(key=lambda x: x[1])

        if results and results[0][0]:
            best_coords = results[0][0]
            for idx, (x, y, z) in best_coords.items():
                if idx < len(mol.atoms):
                    mol.atoms[idx].x = x
                    mol.atoms[idx].y = y
                    mol.atoms[idx].z = z

        return results

    @staticmethod
    def _serialize_mol(mol):
        """Serialize molecule to picklable dict."""
        atoms = [
            (a.symbol, a.formal_charge, a.num_implicit_h)
            for a in mol.atoms
        ]
        bonds = [
            (b.begin_atom_idx, b.end_atom_idx, int(b.order))
            for b in mol.bonds
        ]
        return {'atoms': atoms, 'bonds': bonds}
