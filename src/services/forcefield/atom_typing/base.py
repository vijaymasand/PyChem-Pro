"""
AtomTyper Protocol — the swap interface between Phase 1's HandCoded typer
and Phase 3's SMARTS typer. Engine and parameter lookups consume only
atom.mmff_type and atom.mmff_class; this Protocol defines who writes them.
"""
from __future__ import annotations
from typing import Protocol

from src.core.domain.models.molecule import Molecule


class AtomTyper(Protocol):
    """Assigns MMFF94 atom types and equivalence classes in-place.

    Implementations MUST be idempotent: calling type_atoms() twice with
    no other modifications between produces the same result.

    Implementations MUST write atom.mmff_type (int 1-99, 0=unmatched)
    and atom.mmff_class (int 1-50, 0=unknown) on every atom in mol.atoms.
    """

    def type_atoms(self, mol: Molecule) -> None:
        """Populate mmff_type and mmff_class on every atom of mol."""
        ...

    def coverage(self, mol: Molecule) -> float:
        """Fraction of atoms with mmff_type > 0 after type_atoms().

        Returns 1.0 for fully typed, 0.0 for nothing recognized.
        """
        ...
