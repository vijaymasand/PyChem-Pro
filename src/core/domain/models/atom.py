"""
Atom class — Represents a single atom in a molecular graph.
"""

from src.core.domain.models.elements import get_element, get_element_by_aromatic


class Chirality:
    """Tetrahedral chirality types."""
    NONE = 0
    COUNTERCLOCKWISE = 1  # @ in SMILES
    CLOCKWISE = 2         # @@ in SMILES


class Atom:
    """
    Represents an atom with all properties needed for molecular modeling.

    Attributes:
        index: Position in the molecule's atom list
        element: Element object from periodic table
        formal_charge: Integer formal charge
        isotope: Mass number (0 = natural abundance)
        chirality: Chirality enum value
        num_implicit_h: Number of implicit hydrogens
        num_explicit_h: Number of explicit hydrogens from SMILES bracket
        is_aromatic: Whether atom is in an aromatic ring
        x, y, z: 3D coordinates (None until generated)
        partial_charge: Gasteiger partial charge (None until computed)
        atom_class: SMILES atom class (:n)
        ring_bonds: Set of ring closure bond indices
        hybridization: Hybridization type (sp, sp2, sp3)
        sybyl_type: SYBYL atom type for MOL2 export
        pdb_name: PDB atom name (e.g., "CA", "N") for protein structures
        res_name: Residue name (e.g., "ALA", "GLY") for protein structures
        chain_id: Chain identifier for protein structures
        res_seq: Residue sequence number for protein structures
        b_factor: B-factor/temperature factor for protein structures
        is_hetatm: Whether atom is from HETATM record in PDB
        ss_type: Secondary structure type (H=helix, E=sheet, C=coil)
    """

    __slots__ = (
        'index', 'element', 'formal_charge', 'isotope', 'chirality',
        'num_implicit_h', 'num_explicit_h', 'is_aromatic',
        'x', 'y', 'z', 'x2d', 'y2d', 'z2d', 'partial_charge', 'atom_class',
        'ring_bonds', '_in_bracket', 'radical_electrons',
        'hybridization', 'sybyl_type', 'pdb_name', 'res_name',
        'chain_id', 'res_seq', 'b_factor', 'is_hetatm', 'ss_type',
        'sasa', 'sasa_points',
        'mmff_type', 'mmff_class',
    )

    def __init__(self, symbol, is_aromatic=False, formal_charge=0,
                 isotope=0, chirality=Chirality.NONE, num_explicit_h=None,
                 atom_class=0, in_bracket=False, radical_electrons=0):
        self.index = -1  # Set when added to molecule

        # Resolve element
        if is_aromatic:
            self.element = get_element_by_aromatic(symbol)
        else:
            self.element = get_element(symbol)

        self.formal_charge = formal_charge
        self.isotope = isotope
        self.chirality = chirality
        self.num_explicit_h = num_explicit_h  # None = not specified (use implicit)
        self.num_implicit_h = 0   # Calculated later
        self.is_aromatic = is_aromatic
        self.radical_electrons = radical_electrons

        # 3D coordinates
        self.x = None
        self.y = None
        self.z = None
        
        # 2D geometry coordinates
        self.x2d = None
        self.y2d = None
        self.z2d = None

        # Computed properties
        self.partial_charge = 0.0
        self.atom_class = atom_class
        self.ring_bonds = set()
        self._in_bracket = in_bracket

        # Assigned during processing
        self.hybridization = None  # sp, sp2, sp3
        self.sybyl_type = None     # For MOL2 output
        
        # PDB-specific attributes (for protein structures)
        self.pdb_name = None       # PDB atom name (e.g., "CA", "N")
        self.res_name = None       # Residue name (e.g., "ALA", "GLY")
        self.chain_id = None       # Chain identifier
        self.res_seq = None        # Residue sequence number
        self.b_factor = None       # B-factor/temperature factor
        self.is_hetatm = False     # Whether atom is from HETATM record
        self.ss_type = None        # Secondary structure type (H, E, C)
        
        self.sasa = 0.0            # Solvent Accessible Surface Area
        self.sasa_points = []      # 3D points representing the accessible surface

        # MMFF94 typing (set by AtomTyper; 0 = unknown)
        self.mmff_type = 0
        self.mmff_class = 0

    @property
    def symbol(self):
        return self.element.symbol

    @property
    def atomic_number(self):
        return self.element.atomic_number

    @property
    def mass(self):
        """Atomic mass (isotope if specified, or standard weight)."""
        iso = self.isotope if self.isotope is not None else 0
        return iso if iso > 0 else self.element.mass

    @property
    def total_h(self):
        """Total number of hydrogens (explicit + implicit)."""
        explicit = self.num_explicit_h if self.num_explicit_h is not None else 0
        return explicit + self.num_implicit_h

    @property
    def coords(self):
        """Return coordinates as tuple, or None if not set."""
        if self.x is not None:
            return (self.x, self.y, self.z)
        return None

    @coords.setter
    def coords(self, xyz):
        """Set coordinates from tuple/list."""
        self.x, self.y, self.z = xyz

    @property
    def has_coords(self):
        return self.x is not None

    def total_degree(self, bonds):
        """Total number of connections (bonds + implicit H)."""
        bond_count = sum(1 for b in bonds if b.begin_atom_idx == self.index or b.end_atom_idx == self.index)
        return bond_count + self.total_h

    def heavy_degree(self, bonds):
        """Number of non-hydrogen connections."""
        return sum(1 for b in bonds if b.begin_atom_idx == self.index or b.end_atom_idx == self.index)

    def total_valence(self, bonds):
        """Sum of bond orders + implicit H."""
        val = 0.0
        for b in bonds:
            if b.begin_atom_idx == self.index or b.end_atom_idx == self.index:
                val += b.order
        val += self.total_h
        return val

    def __repr__(self):
        charge_str = ""
        if self.formal_charge > 0:
            charge_str = f"+{self.formal_charge}"
        elif self.formal_charge < 0:
            charge_str = str(self.formal_charge)
        arom_str = " (arom)" if self.is_aromatic else ""
        return f"Atom({self.symbol}{charge_str}{arom_str}, idx={self.index})"
