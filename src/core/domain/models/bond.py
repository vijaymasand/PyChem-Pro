"""
Bond class — Represents a chemical bond between two atoms.
"""


class BondStereo:
    """Bond stereochemistry types for E/Z isomerism."""
    NONE = 0
    UP = 1       # / in SMILES
    DOWN = 2     # \ in SMILES


class BondType:
    """Bond order constants."""
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 4  # Stored as aromatic, resolved to 1.5 for calculations
    AMIDE = 5     # Peptide/amide bond


# Bond order mapping for numerical calculations
BOND_ORDER_MAP = {
    BondType.SINGLE: 1.0,
    BondType.DOUBLE: 2.0,
    BondType.TRIPLE: 3.0,
    BondType.AROMATIC: 1.5,
    BondType.AMIDE: 1.0,
}


class Bond:
    """
    Represents a chemical bond between two atoms.

    Attributes:
        begin_atom_idx: Index of the first atom
        end_atom_idx: Index of the second atom
        bond_type: BondType constant
        stereo: BondStereo for E/Z designation
        is_in_ring: Whether this bond is part of a ring
        is_rotatable: Whether bond allows free rotation
        index: Position in molecule's bond list
    """

    __slots__ = (
        'begin_atom_idx', 'end_atom_idx', 'bond_type',
        'stereo', 'is_in_ring', 'is_rotatable', 'index'
    )

    def __init__(self, begin_atom_idx, end_atom_idx, bond_type=BondType.SINGLE,
                 stereo=BondStereo.NONE):
        self.begin_atom_idx = begin_atom_idx
        self.end_atom_idx = end_atom_idx
        self.bond_type = bond_type
        self.stereo = stereo
        self.is_in_ring = False
        self.is_rotatable = False
        self.index = -1

    @property
    def order(self):
        """Numerical bond order for calculations."""
        return BOND_ORDER_MAP.get(self.bond_type, 1.0)

    @property
    def is_single(self):
        return self.bond_type == BondType.SINGLE

    @property
    def is_double(self):
        return self.bond_type == BondType.DOUBLE

    @property
    def is_triple(self):
        return self.bond_type == BondType.TRIPLE

    @property
    def is_aromatic(self):
        return self.bond_type == BondType.AROMATIC

    @property
    def is_amide(self):
        return self.bond_type == BondType.AMIDE

    def other_atom(self, atom_idx):
        """Return the index of the atom on the other end of the bond."""
        if atom_idx == self.begin_atom_idx:
            return self.end_atom_idx
        elif atom_idx == self.end_atom_idx:
            return self.begin_atom_idx
        else:
            raise ValueError(f"Atom {atom_idx} is not part of this bond")

    def contains_atom(self, atom_idx):
        """Check if the given atom is part of this bond."""
        return atom_idx == self.begin_atom_idx or atom_idx == self.end_atom_idx

    def __repr__(self):
        type_map = {
            BondType.SINGLE: '-',
            BondType.DOUBLE: '=',
            BondType.TRIPLE: '#',
            BondType.AROMATIC: ':',
        }
        sym = type_map.get(self.bond_type, '?')
        return f"Bond({self.begin_atom_idx}{sym}{self.end_atom_idx})"
