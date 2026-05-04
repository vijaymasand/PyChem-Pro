# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3)
from math import cos, sin, atan2
from math import pi as PI

from .drawing_parents import DrawableObject
from .atom import Atom
from .bond import Bond
from src.vendors.oasa.molecule import molecule as OasaMolecule
from . import common
from . import geometry as geo

global molecule_id_no
molecule_id_no = 1

class Molecule(OasaMolecule, DrawableObject):
    redraw_priority = 1
    meta__undo_properties = ("scale_val",)
    meta__undo_copy = ("atoms", "bonds")
    meta__undo_children_to_record = ("atoms", "bonds")
    meta__same_objects = {"vertices": "atoms", "edges": "bonds"}
    meta__scalables = ("scale_val",)
    class_name = "Molecule"

    def __init__(self):
        DrawableObject.__init__(self)
        OasaMolecule.__init__(self)
        self.name = ""
        self.atoms = self.vertices
        self.bonds = self.edges
        self.delocalizations = []
        global molecule_id_no
        self.id = 'mol' + str(molecule_id_no)
        molecule_id_no += 1
        self._last_used_atom = None
        self._sign = 1
        self.scale_val = 1.0
        # Initialize OASA molecule properties
        self._fix_valency_errors = False

    @property
    def children(self):
        return self.atoms + list(self.bonds)

    def new_atom(self, symbol="C"):
        atom = Atom(symbol)
        self.add_atom(atom)
        return atom

    def new_bond(self):
        bond = Bond()
        self.add_bond(bond)
        return bond

    def add_atom(self, atom):
        self.add_vertex(atom)
        atom.molecule = self

    def remove_atom(self, atom):
        self.delete_vertex(atom)
        atom.molecule = None

    def add_bond(self, bond):
        self.edges.add(bond)
        bond.molecule = self
        if hasattr(self, "_flush_cache"):
            self._flush_cache()

    def remove_bond(self, bond):
        if bond in self.edges:
            self.edges.remove(bond)
        bond.molecule = None
        if hasattr(self, "_flush_cache"):
            self._flush_cache()

    def delete_from_paper(self):
        for atom in self.atoms:
            atom.delete_from_paper()
        for bond in list(self.bonds):
            bond.delete_from_paper()
        super().delete_from_paper()

    def eat_molecule(self, food_mol):
        if food_mol is self: return
        for atom in food_mol.atoms[:]:
            food_mol.remove_atom(atom)
            self.add_atom(atom)
        for bond in list(food_mol.bonds):
            food_mol.remove_bond(bond)
            self.add_bond(bond)
        if food_mol.paper:
            food_mol.paper.removeObject(food_mol)

    def split_fragments(self):
        frags = list(self.get_connected_components())
        if len(frags) <= 1: return []
        new_mols = []
        for frag in frags[1:]:
            new_mol = Molecule()
            if self.paper: self.paper.addObject(new_mol)
            for atom in frag:
                self.remove_atom(atom)
                new_mol.add_atom(atom)
                for bond in list(atom.neighbor_edges):
                    if bond in self.bonds:
                        self.remove_bond(bond)
                        new_mol.add_bond(bond)
            new_mols.append(new_mol)
        return new_mols

    def find_place(self, a, distance, added_order=1):
        neighbors = a.neighbors
        if len(neighbors) == 0:
            x = a.x + cos(PI / 6) * distance
            y = a.y - sin(PI / 6) * distance
        elif len(neighbors) == 1:
            neigh = neighbors[0]
            angle = atan2(neigh.y - a.y, neigh.x - a.x)
            x = a.x + cos(angle + PI) * distance
            y = a.y + sin(angle + PI) * distance
        else:
            angles = [atan2(n.y - a.y, n.x - a.x) for n in neighbors]
            angles.sort()
            angles.append(angles[0] + 2 * PI)
            max_diff = 0
            best_angle = 0
            for i in range(len(angles) - 1):
                diff = angles[i + 1] - angles[i]
                if diff > max_diff:
                    max_diff = diff
                    best_angle = angles[i] + diff / 2
            x = a.x + cos(best_angle) * distance
            y = a.y + sin(best_angle) * distance
        return x, y

    def bounding_box(self):
        if not self.atoms: return [0, 0, 0, 0]
        bboxes = [atom.bounding_box() for atom in self.atoms]
        return common.bbox_of_bboxes(bboxes)

    def handle_overlap(self):
        to_process = self.atoms[:]
        replacement_dict = {}
        while len(to_process):
            a1 = to_process.pop(0)
            i = 0
            while i < len(to_process):
                a2 = to_process[i]
                if abs(a2.x - a1.x) <= 5 and abs(a2.y - a1.y) <= 5:
                    replacement_dict[a2] = a1
                    to_process.pop(i)
                    for bond in list(a2.neighbor_edges):
                        neighbor = bond.atom_connected_to(a2)
                        if neighbor is a1 or neighbor in a1.neighbors:
                            # If it's a bond between the merging atoms, or a duplicate bond to a common neighbor, delete it.
                            bond.disconnect_atoms()
                            self.remove_bond(bond)
                            bond.delete_from_paper()
                        else:
                            bond.replace_atom(a2, a1)
                else:
                    i += 1
        for atom in replacement_dict.keys():
            self.remove_atom(atom)
            atom.delete_from_paper()

    def scale(self, factor, center=None):
        if center is None: center = self.get_center()
        cx, cy = center
        for a in self.atoms:
            a.x = cx + (a.x - cx) * factor
            a.y = cy + (a.y - cy) * factor
        self.scale_val *= factor
        self.draw()

    def get_center(self):
        if not self.atoms: return 0, 0
        avg_x = sum(a.x for a in self.atoms) / len(self.atoms)
        avg_y = sum(a.y for a in self.atoms) / len(self.atoms)
        return avg_x, avg_y

    def bounding_box(self):
        if not self.atoms: return [0, 0, 0, 0]
        x1 = min(a.x for a in self.atoms)
        y1 = min(a.y for a in self.atoms)
        x2 = max(a.x for a in self.atoms)
        y2 = max(a.y for a in self.atoms)
        return [x1 - 10, y1 - 10, x2 + 10, y2 + 10]

    def rotate(self, angle, center=None):
        if center is None: center = self.get_center()
        cx, cy = center
        s, c = sin(angle), cos(angle)
        for a in self.atoms:
            x, y = a.x - cx, a.y - cy
            a.x = x * c - y * s + cx
            a.y = x * s + y * c + cy

    def flip_horizontal(self, center_x):
        for a in self.atoms:
            a.x = 2 * center_x - a.x
        for b in self.bonds:
            b.flip_horizontal(center_x)
        self.draw()

    def flip_vertical(self, center_y):
        for a in self.atoms:
            a.y = 2 * center_y - a.y
        for b in self.bonds:
            b.flip_vertical(center_y)
        self.draw()

    def get_bond_between_atoms(self, a1, a2):
        for bond in self.bonds:
            if bond.atom_connected_to(a1) == a2:
                return bond
        return None

    def move_by(self, dx, dy):
        # Update all atom positions first
        for atom in self.atoms:
            atom.x += dx
            atom.y += dy
        # Then redraw everything to avoid stretching bonds during the update
        if self.paper:
            self.draw()

    def clone(self):
        new_mol = Molecule()
        new_mol.scale_val = self.scale_val
        atom_map = {}
        for atom in self.atoms:
            new_atom = atom.copy()
            new_atom.molecule = new_mol
            new_atom.color = atom.color
            new_mol.add_atom(new_atom)
            atom_map[atom] = new_atom
        for bond in self.bonds:
            new_bond = Bond()
            new_bond.set_type(bond.type)
            new_bond.color = bond.color
            new_bond.connect_atoms(atom_map[bond.atoms[0]], atom_map[bond.atoms[1]])
            new_mol.add_bond(new_bond)
        return new_mol

    def set_selected(self, select):
        for atom in self.atoms:
            atom.set_selected(select)
        for bond in self.bonds:
            bond.set_selected(select)

    def draw(self):
        for atom in self.atoms:
            atom.draw()
        for bond in self.bonds:
            bond.draw()

