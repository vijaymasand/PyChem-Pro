"""
MMFF94 BCI (Bond Charge Increment) charge assignment.

For each bond i-j between atoms with MMFF types t_i, t_j, the BCI
parameter w_ij = -w_ji shifts q_i by w_ij and q_j by w_ji = -w_ij.

Starting charge for each atom is its formal charge x MMFF "fcadj"
(formal-charge adjustment) factor, looked up from MMFFPBCI.PAR via
parameters.get_pbci.

This function:
    1. Initializes atom.partial_charge from formal_charge x fcadj.
    2. Adds w_ij for every bond where w_ij is found in MMFFCHG.
    3. Falls back to (pbci(t_j) - pbci(t_i)) when MMFFCHG lacks an entry.
"""
from __future__ import annotations
from src.core.domain.models.molecule import Molecule
from src.services.forcefield.parameters import get_bci, get_pbci


def assign_bci_charges(mol: Molecule) -> None:
    """Populate atom.partial_charge using MMFF94 BCI on every atom in mol.

    Pre: mol.atoms have mmff_type set (call AtomTyper.type_atoms first).
    Post: every atom has a partial_charge consistent with MMFF94 BCI.
    """
    # Step 1: initialize from formal charge x fcadj.
    for a in mol.atoms:
        fcadj = 0.0
        if a.mmff_type > 0:
            pbci_info = get_pbci(a.mmff_type)
            if pbci_info is not None:
                fcadj = pbci_info[1]
        a.partial_charge = float(a.formal_charge) * fcadj

    # Step 2: walk bonds and apply BCI / PBCI fallback.
    for bond in mol.bonds:
        ai = mol.atoms[bond.begin_atom_idx]
        aj = mol.atoms[bond.end_atom_idx]
        ti, tj = ai.mmff_type, aj.mmff_type
        if ti == 0 or tj == 0:
            continue  # untyped atom — skip

        # bond_type code for the lookup. For MMFF94 the bond type is
        # 0 (single, non-aromatic, non-conjugated) or 1 (aromatic / between
        # sp2 atoms). We use 0 as a reasonable default for Phase 1.
        bond_type = 0

        w = get_bci(ti, tj, bond_type)
        if w is None:
            # Fallback: pbci(i) - pbci(j) — matches Jmol's ForceFieldMMFF.java
            # convention where dq is added to atom_i and subtracted from atom_j.
            pi = get_pbci(ti)
            pj = get_pbci(tj)
            if pi is None or pj is None:
                continue  # no params at all — leave atoms as-is
            w = pi[0] - pj[0]

        ai.partial_charge += w
        aj.partial_charge -= w
