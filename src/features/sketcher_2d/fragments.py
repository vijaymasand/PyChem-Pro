# -*- coding: utf-8 -*-
""" Expansion of atom labels (OH, Ph, COOH, Boc ...) into real atoms.

Typing a group label on an atom is how structures are built quickly in ChemDraw.
Instead of keeping the label as a piece of text - which the SMILES export and the
3D import could not make sense of - the label is expanded into the atoms and
bonds it stands for, so everything downstream keeps working.
"""
from math import atan2, cos, sin, pi

from .app_data import periodic_table

# labels chemists type that OASA's own abbreviation table does not carry
extra_groups = {
    "OH": "O", "HO": "O", "SH": "S", "HS": "S", "NH2": "N", "H2N": "N", "NH": "N",
    "CH3": "C", "H3C": "C", "CH2": "C", "D": "[2H]",
    "OMe": "OC", "MeO": "OC", "OEt": "OCC", "EtO": "OCC", "OiPr": "OC(C)C",
    "OtBu": "OC(C)(C)C", "OBn": "OCc1ccccc1", "OPh": "Oc1ccccc1",
    "NO2": "N(=O)=O", "O2N": "N(=O)=O", "NO": "N=O", "CN": "C#N", "NC": "C#N",
    "COOH": "C(=O)O", "CO2H": "C(=O)O", "HOOC": "C(=O)O", "CHO": "C=O", "OHC": "C=O",
    "COOMe": "C(=O)OC", "CO2Me": "C(=O)OC", "COOEt": "C(=O)OCC", "CO2Et": "C(=O)OCC",
    "COCl": "C(=O)Cl", "CONH2": "C(=O)N", "SO2NH2": "S(=O)(=O)N",
    "SO3H": "S(=O)(=O)O", "SO2Cl": "S(=O)(=O)Cl",
    "CF3": "C(F)(F)F", "F3C": "C(F)(F)F", "CCl3": "C(Cl)(Cl)Cl", "CBr3": "C(Br)(Br)Br",
    "OTf": "OS(=O)(=O)C(F)(F)F", "OTs": "OS(=O)(=O)c1ccc(C)cc1", "OMs": "OS(=O)(=O)C",
    "OAc": "OC(=O)C", "AcO": "OC(=O)C", "NHAc": "NC(=O)C", "NHBoc": "NC(=O)OC(C)(C)C",
    "NMe2": "N(C)C", "NEt2": "N(CC)CC", "SMe": "SC", "SPh": "Sc1ccccc1",
}


def known_labels():
    """ every label the sketcher can expand, sorted for a combo box """
    names = set(extra_groups)
    names.update(_oasa_groups())
    return sorted(names, key=lambda s: (len(s), s.lower()))


def _oasa_groups():
    try:
        from src.vendors.oasa.known_groups import name_to_smiles
        return dict(name_to_smiles)
    except Exception:
        return {}


def group_to_smiles(label):
    """ SMILES of a group label, None when the label is not a known group """
    label = (label or "").strip()
    if not label or label in periodic_table:
        return None
    if label in extra_groups:
        return extra_groups[label]
    groups = _oasa_groups()
    if label in groups:
        return groups[label]
    # allow a case insensitive match as a last resort (ph -> Ph)
    for name, smiles in list(extra_groups.items()) + list(groups.items()):
        if name.lower() == label.lower():
            return smiles
    return None


def is_expandable(label):
    return group_to_smiles(label) is not None


def free_direction(atom):
    """ the direction with the most room around atom (screen coordinates) """
    neighbors = atom.neighbors
    if not neighbors:
        return -pi / 2
    if len(neighbors) == 1:
        return atan2(atom.y - neighbors[0].y, atom.x - neighbors[0].x)
    angles = sorted(atan2(nb.y - atom.y, nb.x - atom.x) for nb in neighbors)
    angles.append(angles[0] + 2 * pi)
    gap, best = 0, angles[0]
    for i in range(len(angles) - 1):
        diff = angles[i + 1] - angles[i]
        if diff > gap:
            gap, best = diff, angles[i] + diff / 2
    return best


def expand_label(atom, label):
    """ replaces atom by the structure the label stands for.

    The bonds of the replaced atom are handed over to the first atom of the
    group, so a substituent simply appears where the atom used to be.
    Returns True when the label was expanded. """
    smiles = group_to_smiles(label)
    if not smiles:
        return False
    mol = atom.molecule
    if mol is None:
        return False
    try:
        from src.features.smiles_parser.services.parser import parse_smiles
        from src.core.domain.models.bond import BondType
        parsed = parse_smiles(smiles)
    except Exception:
        return False

    heavy = [a for a in parsed.atoms if a.symbol != 'H']
    if not heavy:
        return False

    length = mol.preferred_bond_length()
    head_index = heavy[0].index
    coords = {}
    for a in heavy:
        x2d = getattr(a, 'x2d', None)
        y2d = getattr(a, 'y2d', None)
        coords[a.index] = (x2d or 0.0, y2d or 0.0)
    ox, oy = coords[head_index]
    coords = {i: (x - ox, y - oy) for i, (x, y) in coords.items()}

    # turn the group so that it grows into the free space around the atom
    body = [p for i, p in coords.items() if i != head_index]
    rotation = 0.0
    if body:
        bx = sum(p[0] for p in body) / len(body)
        by = sum(p[1] for p in body) / len(body)
        if bx or by:
            rotation = free_direction(atom) - atan2(by, bx)
    cos_r, sin_r = cos(rotation), sin(rotation)

    x0, y0 = atom.x, atom.y
    new_atoms = {}
    for a in heavy:
        rx, ry = coords[a.index]
        new_atom = mol.new_atom(a.symbol)
        new_atom.set_pos(x0 + length * (rx * cos_r - ry * sin_r),
                         y0 + length * (rx * sin_r + ry * cos_r))
        if getattr(a, 'formal_charge', 0):
            new_atom.set_charge(a.formal_charge)
        new_atoms[a.index] = new_atom

    types = {BondType.DOUBLE: "double", BondType.TRIPLE: "triple"}
    for b in parsed.bonds:
        if b.begin_atom_idx not in new_atoms or b.end_atom_idx not in new_atoms:
            continue
        new_bond = mol.new_bond()
        new_bond.set_type(types.get(b.bond_type, "single"))
        new_bond.connect_atoms(new_atoms[b.begin_atom_idx], new_atoms[b.end_atom_idx])

    head = new_atoms[head_index]
    for bond in list(atom.neighbor_edges):
        bond.replace_atom(atom, head)
    mol.remove_atom(atom)
    atom.delete_from_paper()
    mol.draw()
    return True
