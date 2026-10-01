"""
Non-bonded exclusion set construction for MMFF94.

Given a molecule, returns the sets of 1-2 (bonded), 1-3 (across one
angle), and 1-4 (across one torsion) atom-index pairs as `frozenset`
sets. These are used by ArraysBuilder to:
    - exclude 1-2 and 1-3 pairs from VdW and ES entirely
    - flag 1-4 pairs for ES with the 0.75 scaling factor (factor 249.0537
      instead of 332.0716)

The sets are construction-time only; they are not part of the runtime
hot path.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Set, FrozenSet

from src.core.domain.models.molecule import Molecule


@dataclass(frozen=True)
class Exclusions:
    pairs_12: Set[FrozenSet[int]]
    pairs_13: Set[FrozenSet[int]]
    pairs_14: Set[FrozenSet[int]]


def build_exclusions(mol: Molecule) -> Exclusions:
    """Return (1-2, 1-3, 1-4) exclusion sets as sets of frozensets."""
    pairs_12: Set[FrozenSet[int]] = set()
    pairs_13: Set[FrozenSet[int]] = set()
    pairs_14: Set[FrozenSet[int]] = set()

    # 1-2: every bond
    for bond in mol.bonds:
        pairs_12.add(frozenset({bond.begin_atom_idx, bond.end_atom_idx}))

    # 1-3: any atom pair connected through exactly one intermediate atom.
    # Walk every atom and pair up its neighbors.
    for atom in mol.atoms:
        neighbors = mol.get_neighbors(atom.index)
        for a in range(len(neighbors)):
            for b in range(a + 1, len(neighbors)):
                p = frozenset({neighbors[a], neighbors[b]})
                if p not in pairs_12:  # ring case: i-j-i should be 1-2 not 1-3
                    pairs_13.add(p)

    # 1-4: bond i-j extended at each end by one more bond
    for bond in mol.bonds:
        a, b = bond.begin_atom_idx, bond.end_atom_idx
        for n_a in mol.get_neighbors(a):
            if n_a == b:
                continue
            for n_b in mol.get_neighbors(b):
                if n_b == a or n_b == n_a:
                    continue
                p = frozenset({n_a, n_b})
                # Only add to 1-4 if not already 1-2 or 1-3
                if p not in pairs_12 and p not in pairs_13:
                    pairs_14.add(p)

    return Exclusions(pairs_12=pairs_12, pairs_13=pairs_13, pairs_14=pairs_14)
