"""
Stereochemistry module — E/Z (cis/trans) and R/S (tetrahedral) assignment.

Processes stereo bonds (/ and \\) from SMILES to assign E/Z geometry,
and @/@@ chirality to assign R/S configuration.
"""

from src.core.domain.models.atom import Chirality


def assign_stereo(molecule, bond_stereo_info):
    """
    Process stereochemistry information from SMILES parsing.

    Args:
        molecule: Molecule object
        bond_stereo_info: list of (bond_idx, stereo_type) from parser

    Updates:
        - Double bonds: marks E/Z based on / and \\ bond patterns
        - Chiral centers: validates @/@@ chirality
    """
    # ── E/Z Stereochemistry ──
    _assign_ez_stereo(molecule, bond_stereo_info)

    # ── R/S Chirality ──
    _validate_chirality(molecule)


def _assign_ez_stereo(molecule, bond_stereo_info):
    """
    Assign E/Z stereochemistry to double bonds based on / and \\ stereo bonds.

    In SMILES: F/C=C/F means trans (E), F/C=C\\F means cis (Z)
    The / and \\ indicate the direction of the bond relative to the carbon.
    """
    if not bond_stereo_info:
        return

    # Find double bonds
    double_bonds = [b for b in molecule.bonds if b.is_double]

    for dbond in double_bonds:
        # Find stereo bonds adjacent to this double bond
        left_atom = dbond.begin_atom_idx
        right_atom = dbond.end_atom_idx

        left_stereo = None
        right_stereo = None

        # Check bonds on left carbon
        for n_idx, bond in molecule.get_neighbor_bonds(left_atom):
            if n_idx == right_atom:
                continue
            for bid, stereo in bond_stereo_info:
                if bid == bond.index:
                    left_stereo = (n_idx, stereo, bond)
                    break

        # Check bonds on right carbon
        for n_idx, bond in molecule.get_neighbor_bonds(right_atom):
            if n_idx == left_atom:
                continue
            for bid, stereo in bond_stereo_info:
                if bid == bond.index:
                    right_stereo = (n_idx, stereo, bond)
                    break

        if left_stereo and right_stereo:
            # Determine E/Z
            # Both / or both \ = trans (E)
            # Mixed / and \ = cis (Z)
            l_stereo = left_stereo[1]
            r_stereo = right_stereo[1]

            # Store E/Z info on the double bond
            # We store as a property rather than modifying bond type
            molecule.properties.setdefault('ez_stereo', {})[dbond.index] = {
                'left_atom': left_atom,
                'right_atom': right_atom,
                'left_substituent': left_stereo[0],
                'right_substituent': right_stereo[0],
                'left_stereo': l_stereo,
                'right_stereo': r_stereo,
                'is_trans': (l_stereo == r_stereo),  # same direction = trans
            }


def _validate_chirality(molecule):
    """
    Validate and process tetrahedral chirality centers.

    @  = counterclockwise viewing from first neighbor
    @@ = clockwise viewing from first neighbor
    """
    for atom in molecule.atoms:
        if atom.chirality == Chirality.NONE:
            continue

        # Chiral center must have 4 different neighbors (including implicit H)
        neighbors = molecule.get_neighbors(atom.index)
        num_h = atom.total_h

        total_neighbors = len(neighbors) + num_h
        if total_neighbors < 3:
            # Not a valid chiral center — clear chirality
            atom.chirality = Chirality.NONE
            continue

        # Store neighbor ordering for 3D coordinate generation
        # The order in SMILES determines the chirality reference
        molecule.properties.setdefault('chirality_neighbors', {})[atom.index] = {
            'neighbors': list(neighbors),
            'num_h': num_h,
            'chirality': atom.chirality,
        }


def get_cip_priority(atom, molecule, visited=None):
    """
    Calculate CIP priority for an atom (simplified version).

    Higher atomic number = higher priority.
    For tie-breaking, considers neighbors recursively.

    Args:
        atom: Atom to calculate priority for
        molecule: Molecule containing the atom
        visited: Set of already-visited atom indices (for recursion)

    Returns:
        Tuple of atomic numbers for comparison (higher = higher priority)
    """
    if visited is None:
        visited = set()

    visited.add(atom.index)
    priority = [atom.atomic_number]

    # Get neighbor priorities
    neighbor_priorities = []
    for n_idx in molecule.get_neighbors(atom.index):
        if n_idx in visited:
            continue
        n_atom = molecule.get_atom(n_idx)
        # For double bonds, count the neighbor twice
        bond = molecule.get_bond_between(atom.index, n_idx)
        count = int(bond.order) if bond else 1
        for _ in range(count):
            neighbor_priorities.append(n_atom.atomic_number)

    # Sort descending (highest priority first)
    neighbor_priorities.sort(reverse=True)
    priority.extend(neighbor_priorities)

    return tuple(priority)
