"""
Per-atom ring information used by HandCodedAtomTyper.

Builds a small RingInfo record for each atom in a molecule:
    in_ring: bool
    ring_sizes: set of ring sizes the atom belongs to
    in_aromatic_ring: bool

This is a thin wrapper around mol.find_rings() / atom.is_aromatic, kept
in atom_typing to avoid bloating the Molecule class.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict

from src.core.domain.models.molecule import Molecule


@dataclass(frozen=True)
class RingInfo:
    in_ring: bool
    ring_sizes: frozenset
    in_aromatic_ring: bool


def build_ring_info(mol: Molecule) -> Dict[int, RingInfo]:
    """Return {atom_index -> RingInfo} for all atoms in mol."""
    # Force ring perception (SSSR).
    rings = mol.find_rings()  # list of lists of atom indices

    per_atom_sizes: Dict[int, set] = {a.index: set() for a in mol.atoms}
    per_atom_aromatic: Dict[int, bool] = {a.index: False for a in mol.atoms}

    for ring in rings:
        size = len(ring)
        # A ring is considered aromatic if all its atoms are aromatic.
        is_arom = all(mol.atoms[i].is_aromatic for i in ring)
        for i in ring:
            per_atom_sizes[i].add(size)
            if is_arom:
                per_atom_aromatic[i] = True

    return {
        a.index: RingInfo(
            in_ring=len(per_atom_sizes[a.index]) > 0,
            ring_sizes=frozenset(per_atom_sizes[a.index]),
            in_aromatic_ring=per_atom_aromatic[a.index],
        )
        for a in mol.atoms
    }
