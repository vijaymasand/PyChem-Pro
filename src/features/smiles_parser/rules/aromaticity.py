"""
Aromaticity module — Kekulization and aromaticity perception.

Implements:
1. Kekulization: Convert aromatic bonds to alternating single/double bonds
2. Aromaticity perception: Identify aromatic rings using Hückel's rule (4n+2)
"""

from src.core.domain.models.bond import BondType


def kekulize(molecule):
    """
    Convert aromatic bonds to alternating single/double bonds (Kekulization).

    Uses a greedy matching algorithm to assign double bonds to aromatic systems
    while satisfying valence rules.

    Args:
        molecule: Molecule with aromatic bonds marked as BondType.AROMATIC
    """
    # Find all aromatic bonds
    aromatic_bonds = [b for b in molecule.bonds if b.bond_type == BondType.AROMATIC]
    if not aromatic_bonds:
        return

    # Find aromatic atoms
    aromatic_atoms = set()
    for bond in aromatic_bonds:
        aromatic_atoms.add(bond.begin_atom_idx)
        aromatic_atoms.add(bond.end_atom_idx)

    # Build aromatic subgraph adjacency
    adj = {a: [] for a in aromatic_atoms}
    for bond in aromatic_bonds:
        adj[bond.begin_atom_idx].append((bond.end_atom_idx, bond.index))
        adj[bond.end_atom_idx].append((bond.begin_atom_idx, bond.index))

    # Determine which atoms "need" a double bond
    # Carbon in aromatic ring needs one double bond for valence 4
    # Nitrogen in aromatic ring: if 2 aromatic neighbors, needs double bond (pyridine-N)
    #                            if 3 aromatic neighbors, doesn't need (pyrrole-N)
    needs_double = set()
    for atom_idx in aromatic_atoms:
        atom = molecule.get_atom(atom_idx)
        sym = atom.symbol
        num_aromatic_neighbors = len(adj[atom_idx])
        num_other_bonds = 0
        has_external_double = False

        for n_idx, bond in molecule.get_neighbor_bonds(atom_idx):
            if bond.bond_type != BondType.AROMATIC:
                num_other_bonds += bond.order
                if bond.is_double:
                    has_external_double = True

        if sym == 'C':
            if not has_external_double:
                needs_double.add(atom_idx)
        elif sym == 'N':
            # Pyridine-type N (2 aromatic bonds, no H): needs double bond
            # Pyrrole-type N (2 aromatic bonds, 1 H): doesn't need double bond
            explicit_h = atom.num_explicit_h if atom.num_explicit_h is not None else 0
            if num_aromatic_neighbors == 2 and explicit_h == 0 and atom.formal_charge == 0:
                # Check if it could be pyrrole-type (contributes 2 electrons)
                # vs pyridine-type (contributes 1 electron)
                # Heuristic: if atom has fewer total bonds, it's pyridine-type
                total_conn = num_aromatic_neighbors + num_other_bonds
                if total_conn <= 2:
                    needs_double.add(atom_idx)
            elif num_aromatic_neighbors == 3:
                # N+ in aromatic ring — needs a double bond
                if atom.formal_charge > 0:
                    needs_double.add(atom_idx)
        elif sym in ('O', 'S', 'Se', 'Te'):
            # Furan-type: doesn't need double bond (contributes lone pair)
            pass
        elif sym == 'B':
            needs_double.add(atom_idx)
        elif sym == 'P':
            if num_aromatic_neighbors == 2:
                needs_double.add(atom_idx)

    # Greedy matching: assign double bonds to maximize matching
    matched = set()
    double_bonds = set()

    # Sort atoms by number of aromatic neighbors (fewest first = most constrained)
    sorted_atoms = sorted(needs_double, key=lambda a: len(adj[a]))

    for atom_idx in sorted_atoms:
        if atom_idx in matched:
            continue

        # Find an unmatched neighbor that also needs a double bond
        for neighbor_idx, bond_idx in adj[atom_idx]:
            if neighbor_idx not in matched and neighbor_idx in needs_double:
                matched.add(atom_idx)
                matched.add(neighbor_idx)
                double_bonds.add(bond_idx)
                break

    # Second pass: try to match remaining unmatched atoms with any neighbor
    for atom_idx in sorted_atoms:
        if atom_idx in matched:
            continue
        for neighbor_idx, bond_idx in adj[atom_idx]:
            if neighbor_idx not in matched:
                matched.add(atom_idx)
                matched.add(neighbor_idx)
                double_bonds.add(bond_idx)
                break

    # Apply the kekulization
    for bond in aromatic_bonds:
        if bond.index in double_bonds:
            bond.bond_type = BondType.DOUBLE
        else:
            bond.bond_type = BondType.SINGLE


def perceive_aromaticity(molecule):
    """
    Perceive aromatic rings using Hückel's rule (4n+2 π electrons).

    Updates atom.is_aromatic and bond.bond_type for identified aromatic systems.

    Args:
        molecule: Molecule with rings already detected
    """
    rings = molecule.find_rings()

    for ring in rings:
        if len(ring) < 3 or len(ring) > 22:
            continue

        # Planarity requirement: every ring atom must be able to be sp2/sp.
        # Detect this from BONDS, not atom.hybridization — hybridization is not
        # assigned until AFTER this pass during parsing, so relying on it let an
        # sp3 ring carbon (e.g. the saturated centre of a dihydropyridine) pass
        # and falsely flagged the whole ring aromatic.  A carbon with no double
        # bond and no aromatic bond has four sigma bonds → sp3 → breaks
        # aromaticity.  (Heteroatoms may be planar via a lone pair, so they are
        # left to the π-electron count.)
        all_planar = True
        for atom_idx in ring:
            atom = molecule.get_atom(atom_idx)
            if atom.symbol == 'C':
                neighbors = molecule.get_neighbor_bonds(atom_idx)
                has_pi = any(b.is_double or b.is_aromatic for _, b in neighbors)
                if not has_pi and not atom.is_aromatic:
                    all_planar = False
                    break

        if not all_planar:
            continue

        # Count pi electrons
        pi_electrons = _count_pi_electrons(molecule, ring)

        # Hückel's rule: 4n + 2
        if pi_electrons >= 2 and (pi_electrons - 2) % 4 == 0:
            # Mark atoms and bonds as aromatic
            for atom_idx in ring:
                molecule.get_atom(atom_idx).is_aromatic = True

            for i in range(len(ring)):
                a = ring[i]
                b = ring[(i + 1) % len(ring)]
                bond = molecule.get_bond_between(a, b)
                if bond:
                    # Mark bond type as AROMATIC for fingerprinting and typing
                    bond.bond_type = BondType.AROMATIC
                    bond.is_in_ring = True


def _count_pi_electrons(molecule, ring):
    """
    Count π electrons contributed by atoms in a ring.

    Rules:
    - C with double bond in ring: 1 electron
    - C sp2 (no double bond in ring but aromatic): 1 electron
    - N with double bond in ring (pyridine): 1 electron
    - N with lone pair (pyrrole): 2 electrons
    - O, S with lone pair (furan, thiophene): 2 electrons
    - B (empty p orbital): 0 electrons
    """
    ring_set = set(ring)
    pi = 0

    for atom_idx in ring:
        atom = molecule.get_atom(atom_idx)
        sym = atom.symbol

        # Check for double bond within the ring
        has_ring_double = False
        for n_idx, bond in molecule.get_neighbor_bonds(atom_idx):
            if n_idx in ring_set and bond.is_double:
                has_ring_double = True
                break

        if sym == 'C':
            if has_ring_double or atom.is_aromatic:
                pi += 1
            else:
                # Check for external double bond (like in azulene)
                has_external_double = any(
                    b.is_double for n, b in molecule.get_neighbor_bonds(atom_idx)
                    if n not in ring_set
                )
                if has_external_double:
                    pi += 1
                else:
                    pi += 0  # sp3 carbon (no ring/external double, not aromatic)
                             # contributes no π electron

        elif sym == 'N':
            if has_ring_double:
                pi += 1  # Pyridine-type
            else:
                # Check total connectivity
                num_bonds = molecule.degree(atom_idx)
                explicit_h = atom.num_explicit_h or 0
                total_conn = num_bonds + explicit_h + atom.num_implicit_h
                if total_conn >= 3 or atom.formal_charge == 0:
                    pi += 2  # Pyrrole-type (lone pair)
                else:
                    pi += 1

        elif sym in ('O', 'S', 'Se', 'Te'):
            pi += 2  # Contribute lone pair

        elif sym == 'B':
            pi += 0  # Empty p orbital

        elif sym == 'P':
            if has_ring_double:
                pi += 1
            else:
                pi += 2

        else:
            # Default: 1 pi electron if has double bond, 2 otherwise
            pi += 1 if has_ring_double else 2

    return pi
