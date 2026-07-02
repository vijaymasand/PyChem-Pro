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
from src.services.forcefield._engine.bond_types import mmff_bond_type_index


def assign_bci_charges(mol: Molecule) -> None:
    """Populate atom.partial_charge using MMFF94 BCI on every atom in mol.

    Pre: mol.atoms have mmff_type set (call AtomTyper.type_atoms first).
    Post: every atom has a partial_charge consistent with MMFF94 BCI.
    """
    # Step 1: initialize q0 from the atom's formal charge.
    #
    # The MMFF94 charge model conserves total charge: Σq must equal the net
    # molecular charge.  The legacy code used ``formal_charge * fcadj`` for
    # q0, which silently *dropped* net charge — e.g. a quaternary ammonium
    # N+ (type 34) has fcadj = 0, so its +1 vanished, and neutral dipolar
    # groups (nitro, N-oxide, sulfoxide) summed to a spurious non-zero total.
    # Using the formal charge directly conserves charge and is identical
    # (0) for the common all-neutral case, so neutral molecules are
    # unaffected while ions and zwitterions become correct.  (MMFF's
    # fractional formal-charge *sharing* over symmetry-equivalent terminal
    # atoms — e.g. −0.5/−0.5 on carboxylate O — is a further refinement not
    # attempted here; the total is correct regardless.)
    for a in mol.atoms:
        a.partial_charge = float(a.formal_charge)

    # Step 2: walk bonds and apply BCI / PBCI fallback.
    for bond in mol.bonds:
        ai = mol.atoms[bond.begin_atom_idx]
        aj = mol.atoms[bond.end_atom_idx]
        ti, tj = ai.mmff_type, aj.mmff_type
        if ti == 0 or tj == 0:
            continue  # untyped atom — skip

        # Use the true MMFF94 bond-type index (1 for delocalised single
        # bonds — biphenyl/butadiene/styrene — which carry their own charge
        # increments).  Fall back to the type-0 increment when the table has
        # no type-1 entry, so this can never do worse than the legacy code.
        bt = mmff_bond_type_index(mol, bond)
        w = get_bci(ti, tj, bt)
        if w is None and bt != 0:
            w = get_bci(ti, tj, 0)
        if w is not None:
            # MMFFCHG stores the increment for the (i,j) record such that the
            # type-i atom gains -w and the type-j atom gains +w — e.g. the
            # (bt=0, 5, 37, -0.15) H–C(aromatic) entry must yield H=+0.15,
            # C=-0.15 (benzene reference).  That is the *opposite* sign to the
            # ``ai += w`` apply used below (and to the pbci fallback, which is
            # correct), so negate the tabulated value.  Without this every
            # table-BCI bond came out reversed — aromatic carbons appeared
            # positive with their hydrogens negative.
            w = -w
        else:
            # Fallback: pbci(i) - pbci(j) — matches Jmol's ForceFieldMMFF.java
            # convention where dq is added to atom_i and subtracted from atom_j.
            pi = get_pbci(ti)
            pj = get_pbci(tj)
            if pi is None or pj is None:
                continue  # no params at all — leave atoms as-is
            w = pi[0] - pj[0]

        ai.partial_charge += w
        aj.partial_charge -= w
