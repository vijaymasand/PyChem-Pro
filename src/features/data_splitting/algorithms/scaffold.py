import collections
from .base import BaseSplitter

# Direct PyChem-Pro internal imports (no pychem_resolver needed)
from src.features.smiles_parser.services.parser import parse_smiles
from src.features.smiles_generator import generate_smiles


def extract_bemis_murcko_scaffold(smiles: str) -> str:
    """Extracts Bemis-Murcko scaffold using PyChem-Pro natively."""
    try:
        mol = parse_smiles(smiles)
        if not mol:
            return ""

        mol = mol.clone()
        # Find all heavy atoms
        heavy_atom_indices = {atom.index for atom in mol.atoms if atom.symbol != 'H'}
        if not heavy_atom_indices:
            return ""

        # Recursively prune heavy atoms of degree <= 1 in the heavy-atom subgraph
        while True:
            to_remove = set()
            for idx in heavy_atom_indices:
                neighbors = mol._adjacency.get(idx, [])
                active_neighbors = [n for n, _ in neighbors if n in heavy_atom_indices]
                if len(active_neighbors) <= 1:
                    to_remove.add(idx)

            if not to_remove:
                break

            heavy_atom_indices -= to_remove

        if not heavy_atom_indices:
            return ""

        # Keep only the active heavy atoms
        all_atom_indices = {atom.index for atom in mol.atoms}
        indices_to_remove = all_atom_indices - heavy_atom_indices

        if indices_to_remove:
            mol.remove_atoms(indices_to_remove)

        return generate_smiles(mol)
    except Exception:
        return ""


class ScaffoldSplitter(BaseSplitter):
    """Groups molecules by Bemis-Murcko scaffolds and distributes the
    groups using a greedy bin-packing solver.  Acyclic structures are
    treated as singletons to prevent biased grouping."""

    def split(self, df, smiles_col, target_col=None, descriptor_df=None):
        n = len(df)

        # Extract scaffolds for all molecules
        scaffolds = []
        for i, row in df.iterrows():
            smiles = str(row[smiles_col])
            scaffold = extract_bemis_murcko_scaffold(smiles)
            if not scaffold:
                # Treat acyclic molecules as singletons so they don't form one huge group
                scaffold = f"acyclic_{i}"
            scaffolds.append(scaffold)

        # Group indices by scaffold
        scaffold_groups = collections.defaultdict(list)
        for idx, scaffold in enumerate(scaffolds):
            scaffold_groups[scaffold].append(idx)

        # Sort groups by size in descending order
        sorted_groups = sorted(scaffold_groups.items(), key=lambda x: len(x[1]), reverse=True)

        train_indices = []
        test_indices = []
        train_target = int(n * self.target_ratio)

        # Greedy bin-packing to distribute scaffolds
        for scaffold, group in sorted_groups:
            if len(train_indices) + len(group) <= train_target:
                train_indices.extend(group)
            elif len(train_indices) < train_target:
                # Choose the option that brings us closer to target
                diff_if_train = abs((len(train_indices) + len(group)) - train_target)
                diff_if_test = abs(len(train_indices) - train_target)
                if diff_if_train < diff_if_test:
                    train_indices.extend(group)
                else:
                    test_indices.extend(group)
            else:
                test_indices.extend(group)

        return train_indices, test_indices
