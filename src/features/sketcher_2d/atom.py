# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3)
from math import cos, sin
from math import pi as PI
import re

from .app_data import (Settings, periodic_table, auto_hydrogen_elements,
                       lone_pair_elements, electron_deficient_elements)
from src.vendors.oasa.atom import atom as OasaAtom
from .drawing_parents import DrawableObject, Color, Font, Align
from .common import list_difference, find_matching_parentheses
from . import geometry as geo
from .app_data import periodic_table as PT, Settings

global atom_id_no
atom_id_no = 1

class Atom(OasaAtom, DrawableObject):
    focus_priority = 6
    redraw_priority = 2
    is_toplevel = False
    meta__undo_properties = ("symbol", "x", "y", "z", "charge", "isotope", "explicit_hydrogens", "show_symbol", "auto_hydrogens", "molecule", "color", "lonepairs", "radical_plus", "radical_minus")
    meta__undo_copy = ("properties_",)
    meta__undo_children_to_record = []
    meta__scalables = ("x", "y", "z")

    def __init__(self, symbol='C'):
        DrawableObject.__init__(self)
        OasaAtom.__init__(self, symbol=symbol)
        self.x, self.y = None, None
        self.z = 0
        # Properties
        self.molecule = None
        self.symbol = symbol
        self.is_group = symbol not in periodic_table
        self.isotope = None
        self.oxidation_num = None
        self.charge = 0
        self.lonepairs = 0
        self.lonepair_type = 'dots'
        self.radical = 0
        self.radical_plus = 0
        self.radical_minus = 0
        self.hydrogens = 0
        self.auto_hydrogens = True
        self.show_symbol = symbol != 'C'
        self.visible = None
        self._hydrogens_text = ""
        self.hydrogen_pos = None
        self.marks_pos = []
        self.mark_types = []
        self._text = None
        self.text_layout = "Auto"
        self._alignment = None
        global atom_id_no
        self.id = 'a' + str(atom_id_no)
        atom_id_no += 1
        self.font_name = Settings.atom_font_name
        self.font_size = Settings.atom_font_size
        self.radical_size = Settings.electron_dot_size
        self.circle_charge = False
        self._main_items = []
        self._mark_items = []
        self._focusable_item = None
        self._focus_item = None
        self._selection_item = None
        self._update_valency()

    def __str__(self):
        return self.id

    # OASA refuses any symbol that is not a real element. The sketcher also has to
    # carry group labels (OH, Ph, COOH...), so the setter is relaxed here.
    def _set_symbol(self, symbol):
        if symbol in periodic_table:
            OasaAtom._set_symbol(self, symbol)
        else:
            self._clean_cache()
            self._symbol = symbol
            self.symbol_number = 0
            self.valency = 1

    symbol = property(OasaAtom._get_symbol, _set_symbol, None, "atom symbol")

    @property
    def parent(self):
        return self.molecule

    @property
    def bonds(self):
        return self.neighbor_edges

    @property
    def free_valency(self):
        return self.valency - self.occupied_valency

    @property
    def occupied_valency(self):
        val = 0
        for bond in self.neighbor_edges:
            val += bond.order
        return int(val)

    @property
    def pos(self):
        return self.x, self.y

    def set_pos(self, x, y):
        self.x, self.y = x, y

    def set_symbol(self, symbol):
        self.symbol = symbol
        self.show_symbol = symbol != "C"
        self.visible = None
        self.is_group = symbol not in periodic_table
        self._text = None
        self._alignment = None
        self.text_layout = "Auto"
        self.isotope = None
        self.auto_hydrogens = True
        self.hydrogens = 0
        self._update_valency()

    def set_hydrogens(self, count):
        if count == -1:
            self.auto_hydrogens = True
        else:
            self.hydrogens = count
            self.auto_hydrogens = False
            if count == 0:
                self.hydrogen_pos = None
        self.update_occupied_valency()
        self._update_hydrogens()

    def set_charge(self, val):
        self.charge = val
        self.marks_pos = []
        self._update_valency()

    def update_visibility(self):
        if not self.molecule or not self.molecule.paper:
            self.visible = False
            return
        if self.show_symbol or not self.neighbors or (hasattr(self.molecule.paper, 'show_carbon') and self.molecule.paper.show_carbon == "All"):
            self.visible = True
            return
        if hasattr(self.molecule.paper, 'show_carbon') and self.molecule.paper.show_carbon == "Terminal" and len(self.neighbors) == 1:
            self.visible = True
            return
        self.visible = False

    def clear_drawings(self):
        if not self.paper: return
        if self._focusable_item:
            try:
                self.paper.removeFocusable(self._focusable_item)
                self.paper.removeItem(self._focusable_item)
            except RuntimeError: pass
            except Exception: pass
            self._focusable_item = None
        for item in self._main_items + self._mark_items:
            try:
                self.paper.removeItem(item)
            except RuntimeError: pass
            except Exception: pass
        self._main_items = []
        self._mark_items = []
        if self._focus_item:
            try:
                self.paper.removeItem(self._focus_item)
            except RuntimeError: pass
            except Exception: pass
            self._focus_item = None
        if self._selection_item:
            try:
                self.paper.removeItem(self._selection_item)
            except RuntimeError: pass
            except Exception: pass
            self._selection_item = None

    def draw(self):
        if not self.molecule or not self.molecule.paper: return
        self.paper = self.molecule.paper
        self.clear_drawings()
        if self.is_group:
            self._draw_functional_group()
        else:
            if self.visible == None:
                self.update_visibility()
            if not self.visible:
                self._draw_invisible_atom()
            else:
                self._draw_visible_atom()
            self._draw_marks()
        if self._focusable_item:
            self.paper.addFocusable(self._focusable_item, self)
        
        # Redraw connected bonds to ensure they use correct trim distance (padding)
        # especially if this atom just changed visibility (e.g. C -> N)
        if self.molecule:
            for bond in self.neighbor_edges:
                bond.draw()

    def _draw_invisible_atom(self):
        r = Settings.bond_length / 3
        rect = self.x - r, self.y - r, self.x + r, self.y + r
        self._focusable_item = self.paper.addEllipse(rect, color=Color.transparent)

    def _draw_visible_atom(self):
        font = Font(self.font_name, self.font_size * self.molecule.scale_val)
        symbol_item = self.paper.addChemicalFormula(html_formula(self.symbol),
            (self.x, self.y), Align.HCenter, 0, font, color=self.color)
        self._main_items = [symbol_item]
        Sx, Sy, Sw, Sh = self.paper.itemBoundingRect(symbol_item)
        if self.hydrogens:
            self._decide_hydrogen_pos()
            H_item = self.paper.addChemicalFormula(self._hydrogens_text,
                (Sx, self.y), Align.Left, 0, font, color=self.color)
            Hx, Hy, Hw, Hh = self.paper.itemBoundingRect(H_item)
            offsets = [(Sw, 0), (0, Hh * 0.7), (-Hw, 0), (0, -Sh * 0.8)]
            self.paper.moveItemsBy([H_item], *offsets[self.hydrogen_pos])
            self._main_items.append(H_item)
        
        # Add a circular white halo behind the text to hide bond ends
        # Use smaller radius for heteroatoms (N, O, S, P, etc.)
        halo_r = 3 if self.symbol != 'C' else 11
        halo_rect = [self.x - halo_r, self.y - halo_r, self.x + halo_r, self.y + halo_r]
        halo = self.paper.addEllipse(halo_rect, color=Color.transparent, fill=Color.white)
        halo.setZValue(-1)
        self._main_items.insert(0, halo)

        focusable_rect = self.paper.itemBoundingBox(symbol_item)
        # Add padding to make it easier to hit
        x1, y1, x2, y2 = focusable_rect
        p = 4
        focusable_rect = [x1 - p, y1 - p, x2 + p, y2 + p]
        self._focusable_item = self.paper.addRect(focusable_rect, color=Color.transparent)

    def _draw_functional_group(self):
        font = Font(self.font_name, self.font_size * self.molecule.scale_val)
        if self._alignment == None: self._update_alignment()
        if self._text == None: self._update_text()
        offset = self.paper.getCharWidth(self.symbol[0], font) / 2
        self._main_items = [self.paper.addChemicalFormula(html_formula(self._text),
            (self.x, self.y), self._alignment, offset, font, color=self.color)]
        rect = self.paper.itemBoundingBox(self._main_items[0])
        # Add padding
        x1, y1, x2, y2 = rect
        p = 4
        rect = [x1 - p, y1 - p, x2 + p, y2 + p]
        self._focusable_item = self.paper.addRect(rect, color=Color.transparent)

    def _draw_marks(self):
        self._decide_marks_pos()
        ax, ay, scale = self.x, self.y, self.molecule.scale_val
        abs_pos = [(ax + dx * scale, ay + dy * scale) for dx, dy in self.marks_pos]
        font = Font(self.font_name, 0.75 * self.font_size * self.molecule.scale_val)
        
        for pos_i, (mark_type, val, angle) in enumerate(self.mark_types):
            x, y = abs_pos[pos_i]
            if mark_type == "charge":
                text = val > 0 and "+" or "−"
                count = abs(val) > 1 and ("%i" % abs(val)) or ""
                text = count + text
                item = self.paper.addHtmlText(text, (x, y), font=font, align=Align.HCenter | Align.VCenter, color=self.color)
                self._mark_items.append(item)
            elif mark_type == "lonepair":
                dot_dist = 4.0 * scale
                dot_r = 1.0 * scale
                px = -sin(angle) * dot_dist / 2
                py = cos(angle) * dot_dist / 2
                x1, y1 = x + px, y + py
                rect1 = [x1 - dot_r, y1 - dot_r, x1 + dot_r, y1 + dot_r]
                item1 = self.paper.addEllipse(rect1, color=self.color, fill=self.color)
                x2, y2 = x - px, y - py
                rect2 = [x2 - dot_r, y2 - dot_r, x2 + dot_r, y2 + dot_r]
                item2 = self.paper.addEllipse(rect2, color=self.color, fill=self.color)
                self._mark_items.extend([item1, item2])
            elif mark_type == "radical_plus":
                text = "•<sup>+</sup>"
                item = self.paper.addHtmlText(text, (x, y), font=font, align=Align.HCenter | Align.VCenter, color=self.color)
                self._mark_items.append(item)
            elif mark_type == "radical_minus":
                text = "•<sup>−</sup>"
                item = self.paper.addHtmlText(text, (x, y), font=font, align=Align.HCenter | Align.VCenter, color=self.color)
                self._mark_items.append(item)

    def bounding_box(self):
        if self._main_items and self.paper:
            bboxes = [self.paper.itemBoundingBox(item) for item in self._main_items]
            from . import common
            return common.bbox_of_bboxes(bboxes)
        return [self.x - 2, self.y - 2, self.x + 2, self.y + 2]

    def move_by(self, dx, dy):
        self.x += dx
        self.y += dy

    def flip_horizontal(self, center_x):
        self.x = 2 * center_x - self.x
        for b in self.neighbor_edges:
            b.flip_horizontal(center_x)
        self.draw()

    def flip_vertical(self, center_y):
        self.y = 2 * center_y - self.y
        for b in self.neighbor_edges:
            b.flip_vertical(center_y)
        self.draw()

    def set_focus(self, focus):
        if not self.paper: return
        if focus:
            if not self._focusable_item: return
            rect = self.paper.itemBoundingBox(self._focusable_item)
            self._focus_item = self.paper.addEllipse(rect, color=Color.black, fill=Settings.focus_color)
            self.paper.toSelectionLayer(self._focus_item)
        else:
            if self._focus_item:
                try:
                    self.paper.removeItem(self._focus_item)
                except RuntimeError:
                    pass
                except Exception:
                    pass
                self._focus_item = None

    def set_selected(self, select):
        if not self.paper: return
        if self._selection_item:
            try: self.paper.removeItem(self._selection_item)
            except RuntimeError: pass
            except Exception: pass
            self._selection_item = None
            
        if select:
            rect = self.x - 4, self.y - 4, self.x + 4, self.y + 4
            self._selection_item = self.paper.addEllipse(rect, fill=Settings.selection_color)
            self.paper.toSelectionLayer(self._selection_item)

    def eat_atom(self, victim):
        """ Merge victim atom into self. """
        if victim is self: return
        mol1 = self.molecule
        mol2 = victim.molecule
        if mol1 != mol2:
            mol1.eat_molecule(mol2)
        
        for bond in list(victim.neighbor_edges):
            if bond.atom_connected_to(victim) in self.neighbors:
                # Duplicate bond
                bond.disconnect_atoms()
                mol1.remove_bond(bond)
                bond.delete_from_paper()
            else:
                bond.replace_atom(victim, self)
        
        mol1.remove_atom(victim)
        victim.delete_from_paper()
        self.on_bond_count_change()

    def move_by(self, dx, dy):
        self.x, self.y = self.x + dx, self.y + dy

    def update_occupied_valency(self):
        val = 0 if self.auto_hydrogens else self.hydrogens
        for bond in self.neighbor_edges:
            val += bond.order
        val = int(val)
        # We don't set self.occupied_valency because it's a read-only property in OasaAtom
        self._update_valency(val)

    def _charge_corrected_valency(self, valency):
        """ a charge changes the number of bonds an atom can carry """
        if not self.charge:
            return valency
        if self.symbol in lone_pair_elements:
            # a lone pair becomes a bonding site when positive (N+ -> 4, O+ -> 3)
            return max(0, valency + self.charge)
        if self.symbol in electron_deficient_elements:
            return max(0, valency - self.charge)
        # C, Si, ... : both the cation and the anion lose one bonding site
        return max(0, valency - abs(self.charge))

    def _update_valency(self, current_occupied=None):
        if current_occupied is None:
            current_occupied = self.occupied_valency

        if not self.auto_hydrogens or self.is_group:
            self.valency = current_occupied
            return
        valencies = [self._charge_corrected_valency(v)
                     for v in periodic_table[self.symbol]["valency"]]
        # first valency (in preference order) that can hold the existing bonds
        self.valency = max(valencies)
        for val in valencies:
            if val >= current_occupied:
                self.valency = val
                break
        self._update_hydrogens()

    def raise_valency_to_senseful_value(self):
        while self.free_valency < 0:
            found = False
            for v in periodic_table[self.symbol]['valency']:
                if v > self.valency:
                    self.valency = v
                    found = True
                    break
            if not found: break

    def _update_hydrogens(self):
        if self.auto_hydrogens:
            if self.symbol in auto_hydrogen_elements:
                self.hydrogens = self.free_valency > 0 and self.free_valency or 0
            else:
                self.hydrogens = 0
        if self.hydrogens:
            self._hydrogens_text = self.hydrogens == 1 and "H" or "H<sub>%i</sub>" % self.hydrogens

    def _update_text(self):
        self._text = self.symbol
        if self.text_layout == "Auto" and self._alignment == Align.Right:
            self._text = get_reverse_formula(self._text)

    def _update_alignment(self):
        p = 0
        for atom in self.neighbors:
            if atom.x < self.x: p -= 1
            elif atom.x > self.x: p += 1
        self._alignment = p > 0 and Align.Right or Align.Left

    def on_bond_count_change(self):
        self.hydrogen_pos = None
        self._alignment = None
        if self.text_layout == "Auto": self._text = None
        self.update_occupied_valency()
        self.visible = None

    def _decide_hydrogen_pos(self):
        if self.hydrogen_pos != None: return
        if len(self.neighbor_edges) == 0:
            self.hydrogen_pos = 0 if self.hydrogens > 2 else 2
            return
        elif len(self.neighbor_edges) == 1:
            self.hydrogen_pos = self.neighbors[0].x > self.x + 1 and 2 or 0
            return
        coords = [(a.x, a.y) for a in self.neighbors]
        angles = [geo.line_get_angle_from_east([self.x, self.y, x, y]) for x, y in coords]
        angles.append(2 * PI + min(angles))
        angles.sort(reverse=True)
        diffs = list_difference(angles)
        i = diffs.index(max(diffs))
        angle = (angles[i] + angles[i + 1]) / 2
        self.hydrogen_pos = int(round(angle * 2 / PI)) % 4

    def _decide_marks_pos(self):
        self.marks_pos.clear()
        self.mark_types.clear()
        
        active_marks = []
        if self.charge != 0:
            active_marks.append(("charge", self.charge))
        if getattr(self, "lonepairs", 0) > 0:
            for _ in range(self.lonepairs):
                active_marks.append(("lonepair", None))
        if getattr(self, "radical_plus", 0):
            active_marks.append(("radical_plus", None))
        if getattr(self, "radical_minus", 0):
            active_marks.append(("radical_minus", None))
            
        if not active_marks:
            return
            
        gaps = []
        if not self.neighbors:
            gaps = [
                (2 * PI, 3 * PI / 2),
                (2 * PI, 0.0),
                (2 * PI, PI),
                (2 * PI, PI / 2)
            ]
        else:
            angles = [(a.x, a.y) for a in self.neighbors]
            angles = [geo.line_get_angle_from_east([self.x, self.y, x, y]) for x, y in angles]
            angles.append(2 * PI + min(angles))
            angles.sort(reverse=True)
            for i in range(len(angles) - 1):
                gap_size = angles[i] - angles[i+1]
                bisector = (angles[i] + angles[i+1]) / 2
                gaps.append((gap_size, bisector))
            gaps.sort(key=lambda x: x[0], reverse=True)
            
        for idx, mark in enumerate(active_marks):
            if idx < len(gaps):
                angle = gaps[idx][1]
            else:
                angle = gaps[0][1] + (idx - len(gaps) + 1) * (PI / 6)
            
            dist = 11
            if self.visible and self.symbol != 'C':
                dist = 13
                
            dx = cos(angle) * dist
            dy = sin(angle) * dist
            self.marks_pos.append((dx, dy))
            self.mark_types.append((mark[0], mark[1], angle))

    def transform(self, tr):
        self.x, self.y = tr.transform(self.x, self.y)

    def copy(self):
        new_atom = Atom(self.symbol)
        for attr in self.meta__undo_properties:
            if hasattr(self, attr):
                setattr(new_atom, attr, getattr(self, attr))
        return new_atom

def html_formula(formula):
    return re.sub(r"\d+", lambda m: "<sub>" + m.group(0) + "</sub>", formula)

def get_reverse_formula(formula):
    parts = []
    i = 0
    while i < len(formula):
        if formula[i] == "(":
            end = find_matching_parentheses(formula, i)
            part = formula[i:end + 1]
            i = end + 1
            while i < len(formula) and formula[i].isdigit():
                part += formula[i]
                i += 1
            parts.append(part)
            continue
        match = re.search(r"[A-Z][a-z]?\d*", formula[i:])
        if match:
            start, end = match.span()
            start, end = i + start, i + end
            if i != start: parts.append(formula[i:start])
            parts.append(formula[start:end])
            i = end
        else:
            parts.append(formula[i:])
            break
    return "".join(reversed(parts))
