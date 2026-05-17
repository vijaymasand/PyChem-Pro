"""
Periodic Table of Elements — Complete chemical data for SMILES processing.

Contains: symbol, name, atomic_number, mass, typical valences, max_valence,
covalent_radius (Å), vdw_radius (Å), electronegativity (Pauling scale),
and CPK color (hex) for 3D rendering.
"""


class Element:
    """Represents a chemical element with all properties needed for molecular modeling."""

    __slots__ = (
        'symbol', 'name', 'atomic_number', 'mass', 'typical_valences',
        'max_valence', 'covalent_radius', 'vdw_radius', 'electronegativity',
        'color', 'aromatic_symbol'
    )

    def __init__(self, symbol, name, atomic_number, mass, typical_valences,
                 max_valence, covalent_radius, vdw_radius, electronegativity,
                 color, aromatic_symbol=None):
        self.symbol = symbol
        self.name = name
        self.atomic_number = atomic_number
        self.mass = mass
        self.typical_valences = typical_valences  # tuple of common valences
        self.max_valence = max_valence
        self.covalent_radius = covalent_radius    # Angstroms
        self.vdw_radius = vdw_radius              # Angstroms
        self.electronegativity = electronegativity  # Pauling scale
        self.color = color                         # hex color for rendering
        self.aromatic_symbol = aromatic_symbol      # lowercase symbol if can be aromatic

    def __repr__(self):
        return f"Element({self.symbol}, Z={self.atomic_number})"


# ──────────────────────────────────────────────────────────────────────
#  Complete Periodic Table Data
#  Sources: IUPAC 2021, Cordero et al. (covalent), Bondi (VdW), Pauling
# ──────────────────────────────────────────────────────────────────────

_ELEMENTS_DATA = [
    # (symbol, name, Z, mass, typical_valences, max_val, cov_r, vdw_r, EN, color, aromatic)
    ("H",  "Hydrogen",    1,   1.008,   (1,),    1, 0.31, 1.20, 2.20, "#d0d0d0", None),
    ("He", "Helium",      2,   4.003,   (0,),    0, 0.28, 1.40, 0.00, "#D9FFFF", None),
    ("Li", "Lithium",     3,   6.941,   (1,),    1, 1.28, 1.82, 0.98, "#CC80FF", None),
    ("Be", "Beryllium",   4,   9.012,   (2,),    2, 0.96, 1.53, 1.57, "#C2FF00", None),
    ("B",  "Boron",       5,  10.81,    (3,),    4, 0.84, 1.92, 2.04, "#FFB5B5", "b"),
    ("C",  "Carbon",      6,  12.011,   (4,),    4, 0.76, 1.70, 2.55, "#55ff7f", "c"),
    ("N",  "Nitrogen",    7,  14.007,   (3, 5),  5, 0.71, 1.55, 3.04, "#3050F8", "n"),
    ("O",  "Oxygen",      8,  15.999,   (2,),    2, 0.66, 1.52, 3.44, "#FF0D0D", "o"),
    ("F",  "Fluorine",    9,  18.998,   (1,),    1, 0.57, 1.47, 3.98, "#90E050", None),
    ("Ne", "Neon",       10,  20.180,   (0,),    0, 0.58, 1.54, 0.00, "#B3E3F5", None),
    ("Na", "Sodium",     11,  22.990,   (1,),    1, 1.66, 2.27, 0.93, "#AB5CF2", None),
    ("Mg", "Magnesium",  12,  24.305,   (2,),    2, 1.41, 1.73, 1.31, "#8AFF00", None),
    ("Al", "Aluminum",   13,  26.982,   (3,),    3, 1.21, 1.84, 1.61, "#BFA6A6", None),
    ("Si", "Silicon",    14,  28.086,   (4,),    4, 1.11, 2.10, 1.90, "#F0C8A0", None),
    ("P",  "Phosphorus", 15,  30.974,   (3, 5),  5, 1.07, 1.80, 2.19, "#FF8000", "p"),
    ("S",  "Sulfur",     16,  32.06,    (2, 4, 6), 6, 1.05, 1.80, 2.58, "#FFFF30", "s"),
    ("Cl", "Chlorine",   17,  35.45,    (1,),    1, 1.02, 1.75, 3.16, "#1FF01F", None),
    ("Ar", "Argon",      18,  39.948,   (0,),    0, 1.06, 1.88, 0.00, "#80D1E3", None),
    ("K",  "Potassium",  19,  39.098,   (1,),    1, 2.03, 2.75, 0.82, "#8F40D4", None),
    ("Ca", "Calcium",    20,  40.078,   (2,),    2, 1.76, 2.31, 1.00, "#3DFF00", None),
    ("Sc", "Scandium",   21,  44.956,   (3,),    3, 1.70, 2.11, 1.36, "#E6E6E6", None),
    ("Ti", "Titanium",   22,  47.867,   (4,),    4, 1.60, 2.00, 1.54, "#BFC2C7", None),
    ("V",  "Vanadium",   23,  50.942,   (5,),    5, 1.53, 2.00, 1.63, "#A6A6AB", None),
    ("Cr", "Chromium",   24,  51.996,   (3, 6),  6, 1.39, 2.00, 1.66, "#8A99C7", None),
    ("Mn", "Manganese",  25,  54.938,   (2, 4, 7), 7, 1.39, 2.00, 1.55, "#9C7AC7", None),
    ("Fe", "Iron",       26,  55.845,   (2, 3),  3, 1.32, 2.00, 1.83, "#E06633", None),
    ("Co", "Cobalt",     27,  58.933,   (2, 3),  3, 1.26, 2.00, 1.88, "#F090A0", None),
    ("Ni", "Nickel",     28,  58.693,   (2,),    2, 1.24, 1.63, 1.91, "#50D050", None),
    ("Cu", "Copper",     29,  63.546,   (1, 2),  2, 1.32, 1.40, 1.90, "#C88033", None),
    ("Zn", "Zinc",       30,  65.38,    (2,),    2, 1.22, 1.39, 1.65, "#7D80B0", None),
    ("Ga", "Gallium",    31,  69.723,   (3,),    3, 1.22, 1.87, 1.81, "#C28F8F", None),
    ("Ge", "Germanium",  32,  72.63,    (4,),    4, 1.20, 2.11, 2.01, "#668F8F", None),
    ("As", "Arsenic",    33,  74.922,   (3, 5),  5, 1.19, 1.85, 2.18, "#BD80E3", None),
    ("Se", "Selenium",   34,  78.96,    (2, 4, 6), 6, 1.20, 1.90, 2.55, "#FFA100", "se"),
    ("Br", "Bromine",    35,  79.904,   (1,),    1, 1.20, 1.85, 2.96, "#A62929", None),
    ("Kr", "Krypton",    36,  83.798,   (0,),    0, 1.16, 2.02, 3.00, "#5CB8D1", None),
    ("Rb", "Rubidium",   37,  85.468,   (1,),    1, 2.20, 3.03, 0.82, "#702EB0", None),
    ("Sr", "Strontium",  38,  87.62,    (2,),    2, 1.95, 2.49, 0.95, "#00FF00", None),
    ("Y",  "Yttrium",    39,  88.906,   (3,),    3, 1.90, 2.00, 1.22, "#94FFFF", None),
    ("Zr", "Zirconium",  40,  91.224,   (4,),    4, 1.75, 2.00, 1.33, "#94E0E0", None),
    ("Nb", "Niobium",    41,  92.906,   (5,),    5, 1.64, 2.00, 1.60, "#73C2C9", None),
    ("Mo", "Molybdenum", 42,  95.96,    (6,),    6, 1.54, 2.00, 2.16, "#54B5B5", None),
    ("Tc", "Technetium", 43,  98.0,     (7,),    7, 1.47, 2.00, 1.90, "#3B9E9E", None),
    ("Ru", "Ruthenium",  44, 101.07,    (4,),    8, 1.46, 2.00, 2.20, "#248F8F", None),
    ("Rh", "Rhodium",    45, 102.91,    (3,),    6, 1.42, 2.00, 2.28, "#0A7D8C", None),
    ("Pd", "Palladium",  46, 106.42,    (2, 4),  4, 1.39, 1.63, 2.20, "#006985", None),
    ("Ag", "Silver",     47, 107.87,    (1,),    1, 1.45, 1.72, 1.93, "#C0C0C0", None),
    ("Cd", "Cadmium",    48, 112.41,    (2,),    2, 1.44, 1.58, 1.69, "#FFD98F", None),
    ("In", "Indium",     49, 114.82,    (3,),    3, 1.42, 1.93, 1.78, "#A67573", None),
    ("Sn", "Tin",        50, 118.71,    (2, 4),  4, 1.39, 2.17, 1.96, "#668080", None),
    ("Sb", "Antimony",   51, 121.76,    (3, 5),  5, 1.39, 2.06, 2.05, "#9E63B5", None),
    ("Te", "Tellurium",  52, 127.60,    (2, 4, 6), 6, 1.38, 2.06, 2.10, "#D47A00", "te"),
    ("I",  "Iodine",     53, 126.90,    (1,),    1, 1.39, 1.98, 2.66, "#940094", None),
    ("Xe", "Xenon",      54, 131.29,    (0,),    0, 1.40, 2.16, 2.60, "#429EB0", None),
]

# Build the periodic table lookup structures
PERIODIC_TABLE = {}          # symbol -> Element
_BY_ATOMIC_NUMBER = {}       # Z -> Element
_AROMATIC_MAP = {}           # aromatic_symbol -> Element

for _data in _ELEMENTS_DATA:
    _elem = Element(*_data)
    PERIODIC_TABLE[_elem.symbol] = _elem
    _BY_ATOMIC_NUMBER[_elem.atomic_number] = _elem
    if _elem.aromatic_symbol:
        _AROMATIC_MAP[_elem.aromatic_symbol] = _elem

# Add remaining elements (55-118) with basic data
_EXTENDED = [
    ("Cs", "Cesium", 55, 132.91, (1,), 1, 2.44, 3.43, 0.79, "#57178F"),
    ("Ba", "Barium", 56, 137.33, (2,), 2, 2.15, 2.68, 0.89, "#00C900"),
    ("La", "Lanthanum", 57, 138.91, (3,), 3, 2.07, 2.00, 1.10, "#70D4FF"),
    ("Ce", "Cerium", 58, 140.12, (3, 4), 4, 2.04, 2.00, 1.12, "#FFFFC7"),
    ("Pr", "Praseodymium", 59, 140.91, (3,), 3, 2.03, 2.00, 1.13, "#D9FFC7"),
    ("Nd", "Neodymium", 60, 144.24, (3,), 3, 2.01, 2.00, 1.14, "#C7FFC7"),
    ("Pm", "Promethium", 61, 145.0, (3,), 3, 1.99, 2.00, 1.13, "#A3FFC7"),
    ("Sm", "Samarium", 62, 150.36, (3,), 3, 1.98, 2.00, 1.17, "#8FFFC7"),
    ("Eu", "Europium", 63, 151.96, (2, 3), 3, 1.98, 2.00, 1.20, "#61FFC7"),
    ("Gd", "Gadolinium", 64, 157.25, (3,), 3, 1.96, 2.00, 1.20, "#45FFC7"),
    ("Tb", "Terbium", 65, 158.93, (3,), 3, 1.94, 2.00, 1.10, "#30FFC7"),
    ("Dy", "Dysprosium", 66, 162.50, (3,), 3, 1.92, 2.00, 1.22, "#1FFFC7"),
    ("Ho", "Holmium", 67, 164.93, (3,), 3, 1.92, 2.00, 1.23, "#00FF9C"),
    ("Er", "Erbium", 68, 167.26, (3,), 3, 1.89, 2.00, 1.24, "#00E675"),
    ("Tm", "Thulium", 69, 168.93, (3,), 3, 1.90, 2.00, 1.25, "#00D452"),
    ("Yb", "Ytterbium", 70, 173.05, (2, 3), 3, 1.87, 2.00, 1.10, "#00BF38"),
    ("Lu", "Lutetium", 71, 174.97, (3,), 3, 1.87, 2.00, 1.27, "#00AB24"),
    ("Hf", "Hafnium", 72, 178.49, (4,), 4, 1.75, 2.00, 1.30, "#4DC2FF"),
    ("Ta", "Tantalum", 73, 180.95, (5,), 5, 1.70, 2.00, 1.50, "#4DA6FF"),
    ("W",  "Tungsten", 74, 183.84, (6,), 6, 1.62, 2.00, 2.36, "#2194D6"),
    ("Re", "Rhenium", 75, 186.21, (7,), 7, 1.51, 2.00, 1.90, "#267DAB"),
    ("Os", "Osmium", 76, 190.23, (4,), 8, 1.44, 2.00, 2.20, "#266696"),
    ("Ir", "Iridium", 77, 192.22, (3, 4), 4, 1.41, 2.00, 2.20, "#175487"),
    ("Pt", "Platinum", 78, 195.08, (2, 4), 4, 1.36, 1.75, 2.28, "#D0D0E0"),
    ("Au", "Gold", 79, 196.97, (1, 3), 3, 1.36, 1.66, 2.54, "#FFD123"),
    ("Hg", "Mercury", 80, 200.59, (1, 2), 2, 1.32, 1.55, 2.00, "#B8B8D0"),
    ("Tl", "Thallium", 81, 204.38, (1, 3), 3, 1.45, 1.96, 1.62, "#A6544D"),
    ("Pb", "Lead", 82, 207.2, (2, 4), 4, 1.46, 2.02, 2.33, "#575961"),
    ("Bi", "Bismuth", 83, 208.98, (3, 5), 5, 1.48, 2.07, 2.02, "#9E4FB5"),
    ("Po", "Polonium", 84, 209.0, (2, 4), 4, 1.40, 1.97, 2.00, "#AB5C00"),
    ("At", "Astatine", 85, 210.0, (1,), 1, 1.50, 2.02, 2.20, "#754F45"),
    ("Rn", "Radon", 86, 222.0, (0,), 0, 1.50, 2.20, 0.00, "#428296"),
    ("Fr", "Francium", 87, 223.0, (1,), 1, 2.60, 3.48, 0.70, "#420066"),
    ("Ra", "Radium", 88, 226.0, (2,), 2, 2.21, 2.83, 0.90, "#007D00"),
    ("Ac", "Actinium", 89, 227.0, (3,), 3, 2.15, 2.00, 1.10, "#70ABFA"),
    ("Th", "Thorium", 90, 232.04, (4,), 4, 2.06, 2.00, 1.30, "#00BAFF"),
    ("Pa", "Protactinium", 91, 231.04, (5,), 5, 2.00, 2.00, 1.50, "#00A1FF"),
    ("U",  "Uranium", 92, 238.03, (6,), 6, 1.96, 1.86, 1.38, "#008FFF"),
    ("Np", "Neptunium", 93, 237.0, (5,), 5, 1.90, 2.00, 1.36, "#0080FF"),
    ("Pu", "Plutonium", 94, 244.0, (4,), 4, 1.87, 2.00, 1.28, "#006BFF"),
    ("Am", "Americium", 95, 243.0, (3,), 3, 1.80, 2.00, 1.30, "#545CF2"),
    ("Cm", "Curium", 96, 247.0, (3,), 3, 1.69, 2.00, 1.30, "#785CE3"),
    ("Bk", "Berkelium", 97, 247.0, (3,), 3, 1.60, 2.00, 1.30, "#8A4FE3"),
    ("Cf", "Californium", 98, 251.0, (3,), 3, 1.60, 2.00, 1.30, "#A136D4"),
    ("Es", "Einsteinium", 99, 252.0, (3,), 3, 1.60, 2.00, 1.30, "#B31FD4"),
    ("Fm", "Fermium", 100, 257.0, (3,), 3, 1.60, 2.00, 1.30, "#B31FBA"),
    ("Md", "Mendelevium", 101, 258.0, (3,), 3, 1.60, 2.00, 1.30, "#B30DA6"),
    ("No", "Nobelium", 102, 259.0, (2, 3), 3, 1.60, 2.00, 1.30, "#BD0D87"),
    ("Lr", "Lawrencium", 103, 266.0, (3,), 3, 1.60, 2.00, 1.30, "#C70066"),
    ("Rf", "Rutherfordium", 104, 267.0, (4,), 4, 1.50, 2.00, 0.00, "#CC0059"),
    ("Db", "Dubnium", 105, 268.0, (5,), 5, 1.50, 2.00, 0.00, "#D1004F"),
    ("Sg", "Seaborgium", 106, 269.0, (6,), 6, 1.50, 2.00, 0.00, "#D90045"),
    ("Bh", "Bohrium", 107, 270.0, (7,), 7, 1.50, 2.00, 0.00, "#E00038"),
    ("Hs", "Hassium", 108, 269.0, (8,), 8, 1.50, 2.00, 0.00, "#E6002E"),
    ("Mt", "Meitnerium", 109, 278.0, (3,), 3, 1.50, 2.00, 0.00, "#EB0026"),
    ("Ds", "Darmstadtium", 110, 281.0, (2,), 2, 1.50, 2.00, 0.00, "#000000"),
    ("Rg", "Roentgenium", 111, 282.0, (1,), 1, 1.50, 2.00, 0.00, "#000000"),
    ("Cn", "Copernicium", 112, 285.0, (2,), 2, 1.50, 2.00, 0.00, "#000000"),
    ("Nh", "Nihonium", 113, 286.0, (1,), 1, 1.50, 2.00, 0.00, "#000000"),
    ("Fl", "Flerovium", 114, 289.0, (2,), 2, 1.50, 2.00, 0.00, "#000000"),
    ("Mc", "Moscovium", 115, 290.0, (1, 3), 3, 1.50, 2.00, 0.00, "#000000"),
    ("Lv", "Livermorium", 116, 293.0, (2,), 2, 1.50, 2.00, 0.00, "#000000"),
    ("Ts", "Tennessine", 117, 294.0, (1,), 1, 1.50, 2.00, 0.00, "#000000"),
    ("Og", "Oganesson", 118, 294.0, (0,), 0, 1.50, 2.00, 0.00, "#000000"),
]

for _data in _EXTENDED:
    _elem = Element(*_data)
    PERIODIC_TABLE[_elem.symbol] = _elem
    _BY_ATOMIC_NUMBER[_elem.atomic_number] = _elem


def get_element(symbol_or_z):
    """Look up an element by symbol (str) or atomic number (int)."""
    if isinstance(symbol_or_z, int):
        elem = _BY_ATOMIC_NUMBER.get(symbol_or_z)
        if elem is None:
            raise ValueError(f"Unknown atomic number: {symbol_or_z}")
        return elem
    elif isinstance(symbol_or_z, str):
        elem = PERIODIC_TABLE.get(symbol_or_z)
        if elem is None:
            # Fallback for alternative cases, e.g. "CL" -> "Cl" in PDB
            elem = PERIODIC_TABLE.get(symbol_or_z.capitalize())
        if elem is None:
            raise ValueError(f"Unknown element symbol: {symbol_or_z}")
        return elem
    else:
        raise TypeError(f"Expected str or int, got {type(symbol_or_z)}")


def get_element_by_aromatic(aromatic_sym):
    """Look up element by aromatic SMILES symbol (e.g., 'c' -> Carbon)."""
    elem = _AROMATIC_MAP.get(aromatic_sym)
    if elem is None:
        raise ValueError(f"Unknown aromatic symbol: {aromatic_sym}")
    return elem


def bond_length(elem1_symbol, elem2_symbol, bond_order=1):
    """Estimate bond length from covalent radii and bond order."""
    e1 = get_element(elem1_symbol)
    e2 = get_element(elem2_symbol)
    r = e1.covalent_radius + e2.covalent_radius
    # Empirical correction for multiple bonds
    if bond_order == 2:
        r *= 0.87
    elif bond_order == 3:
        r *= 0.78
    elif bond_order == 1.5:  # aromatic
        r *= 0.93
    return r
