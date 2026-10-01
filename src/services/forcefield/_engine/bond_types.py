"""
MMFF94 bond-type index (BTIJ) — 0 or 1.

MMFF94 distinguishes two kinds of *single* bond when it looks up bond
stretch, charge-increment and (indirectly) angle/torsion parameters:

    index 0 — an ordinary single bond;
    index 1 — a single bond that joins two atoms which can *delocalise*
              (both sp²/sp or aromatic), e.g. the central bond of
              butadiene, the aryl–vinyl bond of styrene, or the
              inter-ring bond of biphenyl.

Halgren defines this via the per-type ``arom`` / ``sbmb`` flags of
MMFFPROP.PAR (single-bond-in-a-multiple-bond system).  Rather than ship
that extra table, we express the *same* condition structurally from data
the pipeline already computes — aromaticity and hybridisation — which is
robust to atom-type coverage and needs no new parameter file:

    a plain single bond (``bond_type == SINGLE`` — aromatic ring bonds are
    stored as ``AROMATIC`` and double/triple bonds are excluded) whose two
    atoms are *both* π-capable (aromatic, sp² or sp) is bond type 1.

Callers must treat this as advisory: always attempt the type-1 parameter
first and fall back to the type-0 entry when the tables have no type-1
value, so an over-eager classification can never make results worse than
the legacy ``bond_type = 0`` behaviour.
"""
from __future__ import annotations

from src.core.domain.models.bond import BondType


_PI_HYB = frozenset({"sp2", "sp"})


def _pi_capable(atom) -> bool:
    """True when the atom carries a π system (aromatic, sp² or sp)."""
    if getattr(atom, "is_aromatic", False):
        return True
    return getattr(atom, "hybridization", None) in _PI_HYB


def mmff_bond_type_index(mol, bond) -> int:
    """Return the MMFF94 bond-type index (0 or 1) for ``bond`` in ``mol``.

    Only plain single bonds can be type 1; ``AROMATIC`` ring bonds and
    multiple bonds are always type 0.
    """
    if getattr(bond, "bond_type", BondType.SINGLE) != BondType.SINGLE:
        return 0
    ai = mol.atoms[bond.begin_atom_idx]
    aj = mol.atoms[bond.end_atom_idx]
    if _pi_capable(ai) and _pi_capable(aj):
        return 1
    return 0
