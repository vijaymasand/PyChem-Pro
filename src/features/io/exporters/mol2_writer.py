"""
Tripos MOL2 Writer — Write molecules in Tripos MOL2 format.

Implements the standard MOL2 format sections:
- MOLECULE
- ATOM (with partial charges and SYBYL atom types)
- BOND
"""


def write_mol2(molecule, filepath=None):
    """
    Write a molecule in Tripos MOL2 format.

    Args:
        molecule: Molecule with 3D coordinates, SYBYL types, and partial charges
        filepath: Path to write file (if None, returns string)

    Returns:
        MOL2 file content as string
    """
    # Ensure SYBYL types are assigned
    molecule.assign_sybyl_types()

    lines = []

    # ── MOLECULE section ──
    lines.append("@<TRIPOS>MOLECULE")
    lines.append(molecule.name or "Untitled")

    num_atoms = len(molecule.atoms)
    num_bonds = len(molecule.bonds)
    num_subst = 1
    num_feat = 0
    num_sets = 0

    lines.append(f" {num_atoms} {num_bonds} {num_subst} {num_feat} {num_sets}")
    lines.append("SMALL")    # molecule type
    lines.append("GASTEIGER")  # charge type
    lines.append("")
    lines.append("")

    # ── ATOM section ──
    lines.append("@<TRIPOS>ATOM")

    for i, atom in enumerate(molecule.atoms):
        atom_id = i + 1
        atom_name = _generate_atom_name(atom, i)
        x = atom.x if atom.x is not None else 0.0
        y = atom.y if atom.y is not None else 0.0
        z = atom.z if atom.z is not None else 0.0
        sybyl_type = atom.sybyl_type or atom.symbol
        subst_id = 1
        subst_name = "MOL"
        charge = atom.partial_charge or 0.0

        # Format: atom_id atom_name x y z atom_type [subst_id subst_name charge]
        line = (
            f"{atom_id:7d} {atom_name:<8s} "
            f"{x:10.4f} {y:10.4f} {z:10.4f} "
            f"{sybyl_type:<8s} "
            f"{subst_id:4d} {subst_name:<8s} "
            f"{charge:10.4f}"
        )
        lines.append(line)

    # ── BOND section ──
    lines.append("@<TRIPOS>BOND")

    for i, bond in enumerate(molecule.bonds):
        bond_id = i + 1
        a1 = bond.begin_atom_idx + 1  # 1-indexed
        a2 = bond.end_atom_idx + 1

        # MOL2 bond types
        if bond.is_aromatic:
            bond_type = "ar"
        elif bond.order == 1.0:
            bond_type = "1"
        elif bond.order == 2.0:
            bond_type = "2"
        elif bond.order == 3.0:
            bond_type = "3"
        elif bond.order == 1.5:
            bond_type = "ar"
        else:
            bond_type = "1"

        # Check if bond is in aromatic ring
        ba = molecule.atoms[bond.begin_atom_idx]
        bb = molecule.atoms[bond.end_atom_idx]
        if ba.is_aromatic and bb.is_aromatic and bond.is_in_ring:
            bond_type = "ar"

        line = f"{bond_id:6d} {a1:4d} {a2:4d} {bond_type:>4s}"
        lines.append(line)

    # ── SUBSTRUCTURE section ──
    lines.append("@<TRIPOS>SUBSTRUCTURE")
    lines.append(f"     1 MOL         1 TEMP              0 ****  ****    0 ROOT")

    content = "\n".join(lines)

    if filepath:
        with open(filepath, 'w', newline='\n') as f:
            f.write(content)

    return content


def _generate_atom_name(atom, index):
    """Generate a unique atom name for MOL2 (e.g., C1, C2, N1, O1, H1)."""
    return f"{atom.symbol}{index + 1}"
