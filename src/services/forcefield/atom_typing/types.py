"""
MMFF94 atom-type metadata.

`MMFF94_TYPE_TO_CLASS` maps numeric type (1-99) to its equivalence class.
The equivalence class is used for parameter fallback: if a specific type-pair
has no entry in the parameter tables, the engine retries with the class-pair.

Source: Jmol's mmff94_atom_types.txt (HType column).
"""

# Type -> equivalence class (Halgren 1996, Jmol's HType column).
# 0 entries (e.g., undefined types) get class 0.
MMFF94_TYPE_TO_CLASS = {
    1:  1,   # CR    alkyl C
    2:  2,   # C=C   vinylic C
    3:  3,   # C=O   carbonyl C
    4:  4,   # CSP   sp C (alkyne)
    5:  5,   # HC    H on C
    6:  6,   # OR    ether/alcohol O
    7:  7,   # O=C   carbonyl O
    8:  8,   # NR    amine N (sp3)
    9:  9,   # N=C   imine N (sp2)
    10: 10,  # NC=O  amide N
    11: 11,  # F     fluorine
    12: 12,  # CL    chlorine
    13: 13,  # BR    bromine
    14: 14,  # I     iodine
    15: 15,  # S     sulfide S
    16: 16,  # S=C   thiocarbonyl S (also sulfoxide)
    17: 17,  # S=O   sulfone S
    18: 16,  # SS    disulfide / thiocarbonyl class
    20: 20,  # C-4ring  C in cyclobutyl
    21: 21,  # HO    H on -OH
    22: 22,  # C-3ring  C in cyclopropyl
    23: 23,  # H on N (amine)
    24: 23,  # H on N+ (ammonium) - same class as amine H
    25: 25,  # P=O phosphate / phosphine oxide
    26: 26,  # P phosphite (PD3)
    28: 23,  # H on N (amide) - same class as amine H
    30: 20,  # C=C in 4-ring - same class as C-4ring
    31: 31,  # H on water O (HOH)
    32: 32,  # O- anion (carboxylate, sulfonate, etc.)
    34: 8,   # NR+ ammonium - same class as amine N
    35: 6,   # alcohol -OH O - same class as ether O
    36: 23,  # H on guanidinium N
    37: 37,  # CB aromatic carbon
    38: 38,  # NPYD pyridine N
    39: 39,  # NPYL pyrrole N
    40: 8,   # aniline N - same class as amine N
    41: 3,   # CO2-/CR4R - same class as carbonyl C
    45: 45,  # NO2 nitro N
    49: 6,   # oxonium O - same class as ether O
    57: 3,   # guanidinium C - same class as carbonyl C
    67: 8,   # N-oxide N - same class as amine N
    68: 8,   # amine N+ - same class as amine N
    71: 23,  # H on -SH
}


def class_for_type(mmff_type: int) -> int:
    """Look up equivalence class for an MMFF94 atom type.

    Returns 0 for unknown types.
    """
    return MMFF94_TYPE_TO_CLASS.get(mmff_type, 0)
