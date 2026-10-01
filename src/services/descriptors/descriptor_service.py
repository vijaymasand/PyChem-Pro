"""
Descriptor calculation service with parallel batch support.
Implements IDescriptorCalculator protocol.
"""
from src.core.parallel import ParallelExecutor
from src.core.domain.models.molecule import Molecule


# Module-level workers (picklable)
def _compute_basic_descriptors(mol_data):
    """Compute basic molecular descriptors from serialized mol data."""
    from src.core.domain.models.molecule import Molecule
    from src.core.domain.models.atom import Atom
    from src.core.domain.models.bond import BondType

    # Rebuild molecule
    mol = Molecule("desc")
    for symbol, charge, mass, impl_h in mol_data['atoms']:
        a = Atom(symbol)
        a.formal_charge = charge
        a.num_implicit_h = impl_h
        mol.add_atom(a)
    for begin, end, bond_type in mol_data['bonds']:
        bt = {1: BondType.SINGLE, 2: BondType.DOUBLE,
              3: BondType.TRIPLE, 4: BondType.AROMATIC}.get(bond_type, BondType.SINGLE)
        mol.add_bond(begin, end, bt)

    return {
        'molecular_weight': mol.molecular_weight(),
        'num_atoms': mol.num_atoms,
        'num_bonds': mol.num_bonds,
        'num_heavy_atoms': mol.count_heavy_atoms,
        'formula': mol.molecular_formula(),
        'num_rings': len(mol.find_rings()),
        'total_charge': mol.total_charge(),
    }


class DescriptorService:
    """Molecular descriptor calculator with parallel batch support."""

    AVAILABLE = [
        'molecular_weight', 'num_atoms', 'num_bonds', 'num_heavy_atoms',
        'formula', 'num_rings', 'total_charge',
    ]

    def __init__(self, executor: ParallelExecutor = None):
        self.executor = executor or ParallelExecutor()

    def available_descriptors(self) -> list[str]:
        return list(self.AVAILABLE)

    def calculate(self, mol: Molecule, descriptor_names: list[str] | None = None) -> dict:
        """Calculate descriptors for a single molecule."""
        result = {
            'molecular_weight': mol.molecular_weight(),
            'num_atoms': mol.num_atoms,
            'num_bonds': mol.num_bonds,
            'num_heavy_atoms': mol.count_heavy_atoms,
            'formula': mol.molecular_formula(),
            'num_rings': len(mol.find_rings()),
            'total_charge': mol.total_charge(),
        }
        if descriptor_names:
            result = {k: v for k, v in result.items() if k in descriptor_names}
        return result

    def calculate_batch(self, molecules: list[Molecule],
                        descriptor_names: list[str] | None = None) -> list[dict]:
        """Calculate descriptors for multiple molecules in parallel."""
        if len(molecules) <= 1:
            return [self.calculate(mol, descriptor_names) for mol in molecules]

        # Serialize molecules for cross-process transfer
        mol_data_list = [self._serialize_mol(mol) for mol in molecules]

        # Compute in parallel
        results = self.executor.map(_compute_basic_descriptors, mol_data_list)

        # Filter by requested names
        if descriptor_names:
            results = [{k: v for k, v in r.items() if k in descriptor_names} for r in results]

        return results

    @staticmethod
    def _serialize_mol(mol):
        atoms = [(a.symbol, a.formal_charge, a.mass, a.num_implicit_h) for a in mol.atoms]
        bonds = [(b.begin_atom_idx, b.end_atom_idx, b.bond_type) for b in mol.bonds]
        return {'atoms': atoms, 'bonds': bonds}
