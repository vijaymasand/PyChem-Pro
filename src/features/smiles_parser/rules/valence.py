"""
Valence module — Calculate implicit hydrogens and validate valence.

Standard valence rules for organic chemistry elements.
"""

from src.core.domain.models.bond import BondType


# Default valences for common elements
# Each element maps to a tuple of allowed valences (in order of preference)
DEFAULT_VALENCES = {
    'H':  (1,),
    'He': (0,),
    'B':  (3,),
    'C':  (4,),
    'N':  (3, 5),
    'O':  (2,),
    'F':  (1,),
    'Ne': (0,),
    'Na': (1,),
    'Mg': (2,),
    'Al': (3,),
    'Si': (4,),
    'P':  (3, 5),
    'S':  (2, 4, 6),
    'Cl': (1,),
    'Ar': (0,),
    'K':  (1,),
    'Ca': (2,),
    'Se': (2, 4, 6),
    'Br': (1,),
    'Kr': (0,),
    'Te': (2, 4, 6),
    'I':  (1, 3, 5, 7),
    'Xe': (0,),
    'As': (3, 5),
    'Ge': (4,),
    'Sn': (2, 4),
    'Pb': (2, 4),
    'Bi': (3, 5),
}


def calculate_implicit_hydrogens(atom, molecule):
    """
    Calculate the number of implicit hydrogens for an atom based on
    its element's typical valence and current bond orders.

    For bracket atoms with explicit H count, the explicit count is used directly.
    For organic subset atoms, implicit H is calculated from the default valence.

    Args:
        atom: Atom object
        molecule: Molecule containing the atom

    Returns:
        Number of implicit hydrogens (int)
    """
    # Bracket atoms with explicit H count override implicit calculation
    if atom._in_bracket and atom.num_explicit_h is not None:
        return 0  # Explicit H specified, no additional implicit H

    symbol = atom.symbol

    # Get default valences for this element
    valences = DEFAULT_VALENCES.get(symbol)
    if valences is None:
        return 0  # Unknown element — no implicit H

    # Calculate current bond order sum
    bond_order_sum = 0.0
    for _, bond in molecule.get_neighbor_bonds(atom.index):
        if bond.bond_type == BondType.AROMATIC:
            bond_order_sum += 1.0  # Aromatic bonds count as 1 for implicit H calc
        else:
            bond_order_sum += bond.order

    # Add explicit H from bracket notation
    explicit_h = atom.num_explicit_h if atom.num_explicit_h is not None else 0

    # Total connections before implicit H
    total = bond_order_sum + explicit_h

    # Adjust for charge
    charge = atom.formal_charge

    # For aromatic atoms, subtract 1 from valence (they have one implicit
    # bond to the aromatic system that isn't explicit in SMILES)
    if atom.is_aromatic:
        # Aromatic atoms may contribute to pi system — adjust valence
        # For C aromatic: valence effectively 3 (one pi bond shared)
        adjusted_valences = []
        for v in valences:
            adjusted_valences.append(v)
        valences = tuple(adjusted_valences)

    # Find the smallest allowed valence that accommodates current bonds
    target_valence = None
    for v in sorted(valences):
        # Charge adjusts the target valence
        adjusted_v = v - charge  # +charge reduces available valence for H
        if adjusted_v >= total:
            target_valence = adjusted_v
            break

    if target_valence is None:
        # If no valence fits, use largest and set implicit H to 0
        return 0

    implicit_h = int(target_valence - total)
    return max(0, implicit_h)


def validate_valence(atom, molecule):
    """
    Check if an atom's valence is valid.

    Returns:
        (is_valid, message) tuple
    """
    symbol = atom.symbol
    valences = DEFAULT_VALENCES.get(symbol)
    if valences is None:
        return True, "No valence rules for this element"

    # Total valence = bond orders + implicit H + explicit H
    bond_order_sum = 0.0
    for _, bond in molecule.get_neighbor_bonds(atom.index):
        bond_order_sum += bond.order

    total_h = atom.total_h
    total_valence = bond_order_sum + total_h

    # Adjust for charge
    max_allowed = max(valences) - atom.formal_charge

    if total_valence > max_allowed + 0.5:  # 0.5 tolerance for floating point
        return False, (f"Atom {atom.index} ({symbol}): valence {total_valence} "
                       f"exceeds maximum {max_allowed}")

    return True, "OK"
