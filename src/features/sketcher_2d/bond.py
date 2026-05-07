# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3)
from math import degrees
import operator
from functools import reduce

from .app_data import Settings
from src.vendors.oasa.bond import bond as OasaBond
from .drawing_parents import DrawableObject, Color, PenStyle, LineCap, Font, Align
from . import geometry as geo
from . import common

global bond_id_no
bond_id_no = 1

class Bond(OasaBond, DrawableObject):
    focus_priority = 6
    redraw_priority = 3
    is_toplevel = False
    meta__undo_properties = ("molecule", "type", "color",
                "second_line_side", "auto_second_line_side", "show_delocalization")
    meta__undo_copy = ("atoms",)
    meta__same_objects = {"vertices":"atoms"}

    def __init__(self):
        DrawableObject.__init__(self)
        OasaBond.__init__(self)
        self.atoms = self.vertices
        self.molecule = None
        self.type = "single"
        global bond_id_no
        self.id = 'b' + str(bond_id_no)
        bond_id_no += 1
        self._main_items = []
        self._focus_item = None
        self._selection_item = None
        self.show_delocalization = True
        self.second_line_side = None
        self.auto_second_line_side = True
        self.line_width = Settings.bond_width
        self.line_spacing = Settings.bond_spacing
        self.coord_head_dimensions = (6 * self.line_width, 2.5 * self.line_width, 2 * self.line_width)
        self.double_length_ratio = 0.75

    def __str__(self):
        return self.id

    @property
    def atom1(self):
        return self.atoms[0]

    @property
    def atom2(self):
        return self.atoms[1]

    @property
    def order(self):
        order_dict = {"double":2, "E_or_Z":2, "triple":3, "delocalized":1.5, "partial":0.5, "hbond":0}
        return order_dict.get(self.type, 1)

    def set_type(self, bond_type):
        if bond_type == self.type: return
        self.auto_second_line_side = True
        self.second_line_side = None
        self.type = bond_type
        for atom in self.atoms:
            atom.update_occupied_valency()

    def connect_atoms(self, atom1, atom2):
        atom1.add_neighbor(atom2, self)
        atom1.on_bond_count_change()
        atom2.add_neighbor(atom1, self)
        atom2.on_bond_count_change()
        self.atoms.clear()
        self.atoms += [atom1, atom2]
        if self.molecule and hasattr(self.molecule, "_flush_cache"):
            self.molecule._flush_cache()

    def disconnect_atoms(self):
        if not self.atoms: return
        self.atoms[0].remove_neighbor(self.atoms[1])
        self.atoms[0].on_bond_count_change()
        self.atoms[1].remove_neighbor(self.atoms[0])
        self.atoms[1].on_bond_count_change()
        self.atoms.clear()
        if self.molecule and hasattr(self.molecule, "_flush_cache"):
            self.molecule._flush_cache()

    def atom_connected_to(self, atom):
        return self.atoms[0] if atom is self.atoms[1] else self.atoms[1]

    def replace_atom(self, target, substitute):
        atom1, atom2 = self.atoms
        self.disconnect_atoms()
        if atom1 is target:
            self.connect_atoms(substitute, atom2)
        else:
            self.connect_atoms(atom1, substitute)

    def set_focus(self, focus):
        if not self.paper: return
        if focus:
            self._focus_item = self.paper.addLine(self.atoms[0].pos + self.atoms[1].pos,
                                width=self.line_width + 8, color=Settings.focus_color)
            self.paper.toSelectionLayer(self._focus_item)
        else:
            if self._focus_item:
                self.paper.removeItem(self._focus_item)
                self._focus_item = None

    def set_selected(self, select):
        if not self.paper: return
        if self._selection_item:
            try: self.paper.removeItem(self._selection_item)
            except: pass
            self._selection_item = None
            
        if select:
            self._selection_item = self.paper.addLine(self.atoms[0].pos + self.atoms[1].pos,
                                    self.line_width + 4, Settings.selection_color)
            self.paper.toSelectionLayer(self._selection_item)
        else:
            if self._selection_item:
                self.paper.removeItem(self._selection_item)
                self._selection_item = None

    def clear_drawings(self):
        if not self.paper: return
        for item in self._main_items:
            self.paper.removeFocusable(item)
            self.paper.removeItem(item)
        self._main_items = []
        if self._focus_item: self.set_focus(False)
        if self._selection_item: self.set_selected(False)

    def draw(self):
        self.clear_drawings()
        self._midline = self._where_to_draw_from_and_to()
        if not self._midline: return
        self.paper = self.molecule.paper
        self._line_width = max(self.line_width * self.molecule.scale_val, 1)
        method = "_draw_%s" % self.type
        if hasattr(self, method):
            getattr(self, method)()
        
        # Add a transparent hit area to make the bond easier to select
        # Expanded to +16 to make them much easier to grab in dense fused rings
        hit_item = self.paper.addLine(self._midline, self._line_width + 16, color=Color.transparent)
        self._main_items.append(hit_item)

        for item in self._main_items:
            self.paper.toBondLayer(item)
            self.paper.addFocusable(item, self)

    def _where_to_draw_from_and_to(self):
        line = self.atoms[0].pos + self.atoms[1].pos
        bbox1 = self.atoms[0].bounding_box()
        bbox2 = self.atoms[1].bounding_box()
        if geo.rect_intersects_rect(bbox1, bbox2): return None
        x1, y1, x2, y2 = line
        if self.atoms[0].visible:
            x1, y1 = geo.circle_get_point(self.atoms[0].pos, 9, self.atoms[1].pos)
        if self.atoms[1].visible:
            x2, y2 = geo.circle_get_point(self.atoms[1].pos, 9, self.atoms[0].pos)
        return [x1, y1, x2, y2]

    def _draw_single(self):
        self._main_items = [self.paper.addLine(self._midline, self._line_width, color=self.color)]

    def _draw_double(self):
        if self.second_line_side == None:
            self.second_line_side = self._calc_second_line_side()
        d = self.second_line_side * self.line_spacing * self.molecule.scale_val
        if self.second_line_side == 0:
            d = 0.5 * self.line_spacing * self.molecule.scale_val
            line0 = calc_second_line(self, self._midline, -d)
            line1 = calc_second_line(self, self._midline, d)
            self._main_items = [self.paper.addLine(line0, self._line_width, color=self.color),
                               self.paper.addLine(line1, self._line_width, color=self.color)]
        else:
            line1 = calc_second_line(self, self._midline, d)
            self._main_items = [self.paper.addLine(self._midline, self._line_width, color=self.color),
                               self.paper.addLine(line1, self._line_width, color=self.color)]

    def _draw_triple(self):
        d = 0.75 * self.line_spacing * self.molecule.scale_val
        line1 = calc_second_line(self, self._midline, d)
        line2 = calc_second_line(self, self._midline, -d)
        self._main_items = [self.paper.addLine(self._midline, self._line_width, color=self.color),
                           self.paper.addLine(line1, self._line_width, color=self.color),
                           self.paper.addLine(line2, self._line_width, color=self.color)]

    def _draw_wedge(self):
        d = 0.5 * self.line_spacing * self.molecule.scale_val
        p1 = geo.line_get_point_at_distance(self._midline, d)
        p2 = geo.line_get_point_at_distance(self._midline, -d)
        p0 = (self._midline[0], self._midline[1])
        self._main_items = [self.paper.addPolygon([p0, p1, p2], color=self.color, fill=self.color)]

    def _draw_hashed_wedge(self):
        d = 0.5 * self.line_spacing * self.molecule.scale_val
        p1 = geo.line_get_point_at_distance(self._midline, d)
        p2 = geo.line_get_point_at_distance(self._midline, -d)
        p0 = (self._midline[0], self._midline[1])
        line_count = 8
        lines = []
        for i in range(line_count):
            t = i / (line_count - 1)
            x1 = p0[0] + (p1[0] - p0[0]) * t
            y1 = p0[1] + (p1[1] - p0[1]) * t
            x2 = p0[0] + (p2[0] - p0[0]) * t
            y2 = p0[1] + (p2[1] - p0[1]) * t
            lines.append([x1, y1, x2, y2])
        self._main_items = [self.paper.addLine(line, 1, self.color) for line in lines]

    def _calc_second_line_side(self):
        if not self.atoms or len(self.atoms) < 2:
            return 0
        p1 = self.atoms[0].pos
        p2 = self.atoms[1].pos
        line = p1 + p2
        
        # Collect all neighbors of the two atoms, excluding the atoms of the bond itself.
        neighbor_atoms = set()
        for atom in self.atoms:
            for n in atom.neighbors:
                if n not in self.atoms:
                    neighbor_atoms.add(n)
        
        if not neighbor_atoms:
            return 0
            
        sides = {}
        for n in neighbor_atoms:
            # Use a small threshold to avoid precision issues with side detection
            s = geo.line_get_side_of_point(line, n.pos, threshold=0.01)
            if s != 0:
                sides[n] = s
                
        if not sides:
            return 0
            
        left_neighbors = [n for n, s in sides.items() if s == 1]
        right_neighbors = [n for n, s in sides.items() if s == -1]
        
        if len(left_neighbors) > len(right_neighbors):
            return 1
        elif len(right_neighbors) > len(left_neighbors):
            return -1
        else:
            # Tie (e.g. fused bond or trans-alkene).
            # Use canonical tie-break based on atom IDs to ensure stability across flips.
            if not left_neighbors or not right_neighbors:
                return 0 # No neighbors on either side that aren't on the line
            
            def atom_sort_key(a):
                try: return int(a.id[1:])
                except: return a.id
                
            min_left = min(atom_sort_key(n) for n in left_neighbors)
            min_right = min(atom_sort_key(n) for n in right_neighbors)
            
            # Prefer the side with the globally "first" atom.
            # This ensures the double bond stays in the same ring after flipping.
            return 1 if min_left < min_right else -1

    def copy(self):
        new_bond = Bond()
        new_bond.type = self.type
        return new_bond

    def get_center(self):
        if not self.atoms or len(self.atoms) < 2: return 0, 0
        return (self.atoms[0].x + self.atoms[1].x) / 2, (self.atoms[0].y + self.atoms[1].y) / 2

    def flip_horizontal(self, center_x):
        if self.second_line_side is not None:
            if self.auto_second_line_side:
                self.second_line_side = None
            else:
                self.second_line_side *= -1

    def flip_vertical(self, center_y):
        if self.second_line_side is not None:
            if self.auto_second_line_side:
                self.second_line_side = None
            else:
                self.second_line_side *= -1

def calc_second_line(bond, mid_line, distance):
    x, y, x0, y0 = geo.line_get_parallel(mid_line, distance)
    dx, dy = x - x0, y - y0
    k = 0.15
    return [x - k * dx, y - k * dy, x0 + k * dx, y0 + k * dy]
