# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3)
from .app_data import App, Settings, periodic_table
from .drawing_parents import Color, Align, PenStyle
from .tool_helpers import draw_objs_recursively
from .molecule import Molecule
from .atom import Atom
from .bond import Bond
from .arrow import Arrow
from .text_label import TextLabel
from src.shared.qt_compat import Qt, QApplication
from . import geometry as geo

toolsettings = {
    'structure': 'C',
    'bond_type': 'single',
    'mode': 'bond', # bond or atom
}

# single key shortcuts that retype the atom under the cursor
element_shortcuts = {
    'c': 'C', 'n': 'N', 'o': 'O', 's': 'S', 'p': 'P',
    'f': 'F', 'l': 'Cl', 'b': 'Br', 'i': 'I', 'h': 'H',
}

def handle_atom_keypress(key, text):
    focused = App.paper.focused_obj
    if isinstance(focused, Atom):
        import time
        current_time = time.time()
        
        # Check for Enter or Return to edit the label via dialog
        if key in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            if App.window and hasattr(App.window, '_menu_edit_label'):
                App.window._menu_edit_label(focused)
                return True
        
        last_atom = getattr(App, 'last_typed_atom', None)
        last_time = getattr(App, 'last_type_time', 0)
        buffer = getattr(App, 'typing_buffer', "")
        
        if key == Qt.Key.Key_Backspace:
            if last_atom is focused and (current_time - last_time) < 1.5 and len(buffer) > 0:
                buffer = buffer[:-1]
                if not buffer:
                    buffer = "C" # Default to Carbon if completely backspaced
                App.typing_buffer = buffer
                App.last_type_time = current_time
                focused.set_symbol(buffer)
                focused.draw()
                App.paper.save_state_to_undo_stack()
                return True
        
        char = text
        if char and (char.isalnum() or char in ('+', '-')):
            mapped_char = element_shortcuts.get(char.lower(), char)
            
            if last_atom is focused and (current_time - last_time) < 1.5:
                buffer += mapped_char
            else:
                buffer = mapped_char
            
            App.last_typed_atom = focused
            App.last_type_time = current_time
            App.typing_buffer = buffer
            
            focused.set_symbol(buffer)
            focused.draw()
            App.paper.save_state_to_undo_stack()
            App.paper.locked_focus_obj = None
            return True
    return False

# ring templates : name -> (size, {index: symbol}, aromatic, indices that never
# carry a double bond). Index 0 is the "first" vertex of the ring.
ring_templates = {
    "cyclopropane":    (3, {}, False, ()),
    "cyclobutane":     (4, {}, False, ()),
    "cyclopentane":    (5, {}, False, ()),
    "cyclohexane":     (6, {}, False, ()),
    "cycloheptane":    (7, {}, False, ()),
    "cyclooctane":     (8, {}, False, ()),
    "benzene":         (6, {}, True, ()),
    "naphthalene":     (6, {}, True, ()),   # second ring is fused on afterwards
    "pyridine":        (6, {0: "N"}, True, ()),
    "pyrrole":         (5, {0: "N"}, True, (0,)),
    "furan":           (5, {0: "O"}, True, (0,)),
    "thiophene":       (5, {0: "S"}, True, (0,)),
    "cyclopentadiene": (5, {}, True, (0,)),
}

class Tool:
    def __init__(self): pass
    def on_mouse_press(self, x, y): pass
    def on_mouse_release(self, x, y): pass
    def on_mouse_move(self, x, y): pass
    def on_mouse_double_click(self, x, y): pass
    def on_key_press(self, key, text): pass
    def on_right_click(self, x, y): pass
    def clear(self): pass
    def show_status(self, msg): pass

class StructureTool(Tool):
    """ Draws atoms and bonds.

    Click on empty paper  -> new atom
    Click on an atom      -> change its symbol, or grow a new bond from it
    Click on a bond       -> apply/cycle the bond type
    Drag from an atom     -> new bond, snapped to 15 degree steps (Shift = free),
                             landing on an existing atom joins the two structures
    """
    snap_angle = 15

    def __init__(self):
        Tool.__init__(self)
        self.reset()
        self.preview_item = None
        self.atom_with_preview_bond = None

    def reset(self):
        self.atom1 = None
        self.mouse_press_pos = None
        self.snap_target = None

    def clear(self):
        self.clear_preview()
        self.reset()

    # ------------------------------------------------------------------ helpers
    @staticmethod
    def snap_radius():
        """ how close to an atom a bond end has to be to connect to it """
        return max(10.0, Settings.bond_length * 0.4)

    def _target_atom(self, x, y):
        exclude = (self.atom1,) if self.atom1 else ()
        return App.paper.find_closest_atom(x, y, self.snap_radius(), exclude=exclude)

    def _snapped_pos(self, atom, x, y):
        """ fixed bond length and 15 degree steps, unless Shift is held """
        from math import atan2, cos, sin, pi, hypot
        if QApplication.keyboardModifiers() & Qt.KeyboardModifier.ShiftModifier:
            return x, y
        length = atom.molecule.preferred_bond_length()
        dx, dy = x - atom.x, y - atom.y
        if hypot(dx, dy) < 1e-6:
            return atom.x + length, atom.y
        step = pi * self.snap_angle / 180.0
        angle = round(atan2(dy, dx) / step) * step
        return atom.x + length * cos(angle), atom.y + length * sin(angle)

    # ------------------------------------------------------------------- events
    def on_mouse_press(self, x, y):
        self.mouse_press_pos = (x, y)
        self.atom1 = None
        focused_obj = App.paper.focused_obj
        if isinstance(focused_obj, Atom):
            self.atom1 = focused_obj

    def on_mouse_move(self, x, y):
        if not App.paper.dragging:
            focused_obj = App.paper.focused_obj
            if isinstance(focused_obj, Atom):
                if focused_obj is not self.atom_with_preview_bond:
                    self.clear_preview()
                    self.show_preview(focused_obj)
            else:
                self.clear_preview()
            return
        # Dragging to draw bond
        if not self.atom1:
            mol = Molecule()
            App.paper.addObject(mol)
            # a group label is expanded on the far end of the bond, the atom the
            # drag starts from stays a plain atom
            symbol = toolsettings['structure']
            self.atom1 = mol.new_atom(symbol if symbol in periodic_table else "C")
            self.atom1.set_pos(*self.mouse_press_pos)
            self.atom1.draw()

        target = self._target_atom(x, y)
        self.snap_target = target
        # the paper freezes the focus while the mouse is down, so the tool
        # highlights the atom the bond would snap to by itself
        App.paper.changeFocusTo(target)
        end = (target.x, target.y) if target else self._snapped_pos(self.atom1, x, y)

        self.clear_preview()
        self.preview_item = App.paper.addLine(self.atom1.pos + end, style=PenStyle.dashed)

    def on_mouse_release(self, x, y):
        if not App.paper.dragging:
            self.on_mouse_click(x, y)
            return

        self.clear_preview()
        if self.atom1:
            target = self._target_atom(x, y)
            if target:
                pos = (target.x, target.y)
            else:
                pos = self._snapped_pos(self.atom1, x, y)
                # the snapped end may still land on top of an existing atom
                target = App.paper.find_closest_atom(pos[0], pos[1], 5, exclude=(self.atom1,))
            self.make_bond(self.atom1, target, pos)
            App.paper.save_state_to_undo_stack()
        App.paper.changeFocusTo(None)
        self.reset()

    def make_bond(self, atom1, atom2, pos):
        """ bonds atom1 to atom2, creating atom2 at pos when it is None.

        Atoms of two different molecules are merged into a single molecule,
        which is what makes it possible to link two separately drawn fragments. """
        mol = atom1.molecule
        created = False
        if atom2 is None:
            atom2 = mol.new_atom("C")
            atom2.set_pos(*pos)
            created = True
        elif atom2 is atom1:
            return
        else:
            if atom2.molecule is mol:
                existing = mol.get_bond_between_atoms(atom1, atom2)
                if existing:
                    # dropping a bond onto an existing one only changes its type
                    self.change_bond(existing)
                    return
            else:
                mol.eat_molecule(atom2.molecule)

        bond = mol.new_bond()
        bond.set_type(toolsettings['bond_type'])
        bond.connect_atoms(atom1, atom2)
        atom1.draw()
        atom2.draw()
        bond.draw()
        if created:
            self.apply_structure(atom2)

    @staticmethod
    def apply_structure(atom):
        """ writes the current element or group label onto atom.

        A group label (OH, Ph, COOH ...) is expanded into real atoms, so that
        SMILES export and the 3D import keep working. """
        symbol = toolsettings['structure']
        if symbol == atom.symbol:
            return
        if symbol not in periodic_table:
            from .fragments import expand_label
            if expand_label(atom, symbol):
                return
        atom.set_symbol(symbol)
        atom.draw()

    def change_bond(self, bond):
        wanted = toolsettings['bond_type']
        if bond.type != wanted:
            bond.set_type(wanted)
        elif wanted in ("wedge", "hashed_wedge"):
            bond.reverse()  # clicking a wedge again flips its direction
        else:
            types = ['single', 'double', 'triple']
            if bond.type in types:
                bond.set_type(types[(types.index(bond.type) + 1) % len(types)])
        bond.draw()
        for a in bond.atoms:
            a.draw()

    def on_mouse_click(self, x, y):
        focused_obj = App.paper.focused_obj
        if not focused_obj:
            mol = Molecule()
            App.paper.addObject(mol)
            atom = mol.new_atom("C")
            atom.set_pos(x, y)
            atom.draw()
            self.apply_structure(atom)
        elif isinstance(focused_obj, Atom):
            # Change symbol or add branch
            if focused_obj.symbol != toolsettings['structure']:
                self.apply_structure(focused_obj)
            else:
                bond_length = focused_obj.molecule.preferred_bond_length()
                nx, ny = focused_obj.molecule.find_place(focused_obj, bond_length)
                self.make_bond(focused_obj, None, (nx, ny))
        elif isinstance(focused_obj, Bond):
            self.change_bond(focused_obj)

        App.paper.save_state_to_undo_stack()
        self.reset()

    def show_preview(self, atom):
        bond_length = atom.molecule.preferred_bond_length()
        nx, ny = atom.molecule.find_place(atom, bond_length)
        self.preview_item = App.paper.addLine(atom.pos + (nx, ny), style=PenStyle.dashed)
        self.atom_with_preview_bond = atom

    def on_key_press(self, key, text):
        return handle_atom_keypress(key, text)

    def clear_preview(self):
        if self.preview_item:
            App.paper.removeItem(self.preview_item)
            self.preview_item = None
            self.atom_with_preview_bond = None

class EraserTool(Tool):
    """ Removes atoms, bonds and other objects. Dragging erases continuously. """

    def __init__(self):
        Tool.__init__(self)
        self.erased = False

    def on_mouse_press(self, x, y):
        self.erased = False
        self._erase_at(x, y)

    def on_mouse_move(self, x, y):
        if App.paper.dragging:
            self._erase_at(x, y)

    def on_mouse_release(self, x, y):
        # the state is saved after the deletion, otherwise undo/redo is off by one
        if self.erased:
            App.paper.save_state_to_undo_stack()
        self.erased = False

    def _erase_at(self, x, y):
        # Refresh focus at the moment of click for maximum precision
        zoom = App.paper.view.transform().m11() if App.paper.view else 1.0
        d = max(2, 8 / zoom)
        focused = App.paper.find_closest_object(x, y, d)
        if focused and delete_object(focused):
            self.erased = True
            App.paper.locked_focus_obj = None
            App.paper.redraw_dirty_objects()


def delete_object(obj):
    """ removes an atom, bond or any other object from the paper.

    Atoms take their bonds with them and a structure that falls apart is split
    into separate molecules. Returns True when something was removed. """
    if obj is None:
        return False
    App.paper.unfocusObject(obj)
    if App.paper.locked_focus_obj is obj:
        App.paper.locked_focus_obj = None

    if isinstance(obj, Atom):
        mol = obj.molecule
        if mol is None:
            return False
        for bond in list(obj.bonds):
            bond.disconnect_atoms()
            mol.remove_bond(bond)
            bond.delete_from_paper()
        mol.remove_atom(obj)
        obj.delete_from_paper()
        if not mol.atoms:
            App.paper.removeObject(mol)
        else:
            # neighbours lost a bond, their hydrogens have to be redrawn
            for frag in [mol] + mol.split_fragments():
                frag.draw()
    elif isinstance(obj, Bond):
        mol = obj.molecule
        if mol is None:
            return False
        obj.disconnect_atoms()
        mol.remove_bond(obj)
        obj.delete_from_paper()
        for frag in [mol] + mol.split_fragments():
            frag.draw()
    elif isinstance(obj, Molecule):
        obj.delete_from_paper()
    elif hasattr(obj, 'delete_from_paper'):
        obj.delete_from_paper()
    else:
        return False
    return True

class RotateTool(Tool):
    def __init__(self):
        self.mouse_press_pos = None
        self.center = None
        self.molecule = None

    def on_mouse_press(self, x, y):
        self.mouse_press_pos = (x, y)
        focused = App.paper.focused_obj
        
        if focused:
            if hasattr(focused, 'molecule') and focused.molecule:
                self.molecule = focused.molecule
            else:
                self.molecule = focused
        
        if not self.molecule:
            molecules = [obj for obj in App.paper.objects if hasattr(obj, 'atoms')]
            if len(molecules) == 1:
                self.molecule = molecules[0]
            elif len(molecules) > 1:
                min_dist = float('inf')
                for m in molecules:
                    cx, cy = m.get_center()
                    dist = (x - cx)**2 + (y - cy)**2
                    if dist < min_dist:
                        min_dist = dist
                        self.molecule = m
                        
        if self.molecule:
            self.center = self.molecule.get_center()

    def on_mouse_move(self, x, y):
        if not App.paper.dragging or not self.molecule: return
        import math
        angle1 = math.atan2(self.mouse_press_pos[1] - self.center[1], self.mouse_press_pos[0] - self.center[0])
        angle2 = math.atan2(y - self.center[1], x - self.center[0])
        angle = angle2 - angle1
        self.molecule.rotate(angle, self.center)
        
        if hasattr(self.molecule, 'atoms'):
            for a in self.molecule.atoms:
                a.on_bond_count_change()
                a.draw()
        if hasattr(self.molecule, 'bonds'):
            for b in self.molecule.bonds:
                b.draw()
        elif hasattr(self.molecule, 'draw'):
            self.molecule.draw()
            
        self.mouse_press_pos = (x, y)

    def on_mouse_release(self, x, y):
        if self.molecule:
            App.paper.save_state_to_undo_stack()
        self.molecule = None

class TemplateTool(Tool):
    """ Places ring templates.

    Click on empty paper -> free standing ring (drag to spin it around)
    Click on an atom     -> ring sharing that atom (spiro / substituent ring)
    Click on a bond      -> ring fused along that bond (naphthalene, indole, ...)
    """

    def __init__(self, template_name):
        self.template_name = template_name
        self.free_ring_mol = None
        self.press_pos = None
        self.changed = False

    def spec(self, name=None):
        return ring_templates.get(name or self.template_name, ring_templates["cyclohexane"])

    # ------------------------------------------------------------------- events
    def on_mouse_press(self, x, y):
        self.press_pos = (x, y)
        self.free_ring_mol = None
        self.changed = False

        focused = App.paper.focused_obj
        if isinstance(focused, Atom):
            # an atom under the cursor wins over the bonds it belongs to,
            # otherwise every substituent would silently fuse to its ring
            self._attach_to_atom(focused, x, y)
        else:
            bond = focused if isinstance(focused, Bond) else self._bond_near(x, y)
            if bond:
                self._fuse_to_bond(bond, x, y)
            else:
                self._place_at(x, y)
        self.changed = True

    def on_mouse_move(self, x, y):
        # dragging right after dropping a free ring spins it around its centre
        if not App.paper.dragging or not self.free_ring_mol:
            return
        from math import atan2
        mol = self.free_ring_mol
        cx, cy = mol.get_center()
        angle = (atan2(y - cy, x - cx) -
                 atan2(self.press_pos[1] - cy, self.press_pos[0] - cx))
        mol.rotate(angle, (cx, cy))
        for atom in mol.atoms:
            atom.on_bond_count_change()
        self.press_pos = (x, y)
        mol.draw()

    def on_mouse_release(self, x, y):
        if self.changed:
            App.paper.save_state_to_undo_stack()
        self.free_ring_mol = None
        self.changed = False

    def clear(self):
        self.free_ring_mol = None
        self.changed = False

    def _bond_near(self, x, y):
        """ closest bond within a few pixels, used when nothing is focused """
        objs = App.paper.objectsInRect([x - 10, y - 10, x + 10, y + 10])
        bond, min_dist = None, 10.0
        for obj in objs:
            if isinstance(obj, Bond) and len(obj.atoms) == 2:
                dist = geo.dist_to_segment((x, y), obj.atoms[0].pos, obj.atoms[1].pos)
                if dist < min_dist:
                    min_dist, bond = dist, obj
        return bond

    # ---------------------------------------------------------------- placement
    def _place_at(self, x, y):
        from math import cos, sin, pi
        mol = Molecule()
        App.paper.addObject(mol)
        n, symbols, aromatic, no_double = self.spec()
        length = Settings.bond_length
        radius = length / (2 * sin(pi / n))
        # first vertex points up, which gives the familiar upright ring
        start = -pi / 2
        positions = [(x + radius * cos(start + i * 2 * pi / n),
                      y + radius * sin(start + i * 2 * pi / n)) for i in range(n)]
        atoms = [mol.new_atom("C") for _ in range(n)]
        for atom, pos in zip(atoms, positions):
            atom.set_pos(*pos)
        self._finish_ring(mol, atoms, self.template_name)
        self.free_ring_mol = mol
        mol.draw()

    def _attach_to_atom(self, atom, x, y):
        """ builds a ring that has `atom` as one of its vertices """
        from math import cos, sin, pi, atan2
        mol = atom.molecule
        n, symbols, aromatic, no_double = self.spec()
        length = mol.preferred_bond_length()
        radius = length / (2 * sin(pi / n))
        angle = self._free_direction(atom, x, y)
        cx, cy = atom.x + radius * cos(angle), atom.y + radius * sin(angle)

        start = atan2(atom.y - cy, atom.x - cx)
        step = 2 * pi / n
        atoms = [atom]
        for i in range(1, n):
            new_atom = mol.new_atom("C")
            new_atom.set_pos(cx + radius * cos(start + i * step),
                             cy + radius * sin(start + i * step))
            atoms.append(new_atom)
        # keep a hetero atom away from the atom we attached to
        self._finish_ring(mol, atoms, self.template_name, symbol_offset=n // 2)

    def _fuse_to_bond(self, bond, x, y, name=None):
        """ builds a ring that shares `bond` with the existing structure """
        from math import cos, sin, pi, tan, atan2, hypot
        name = name or self.template_name
        mol = bond.molecule
        a1, a2 = bond.atoms
        n, symbols, aromatic, no_double = self.spec(name)
        length = geo.point_distance(a1.pos, a2.pos) or mol.preferred_bond_length()

        mx, my = (a1.x + a2.x) / 2, (a1.y + a2.y) / 2
        bx, by = a2.x - a1.x, a2.y - a1.y
        blen = hypot(bx, by) or 1.0
        nx, ny = -by / blen, bx / blen           # unit normal of the bond
        # grow the ring on the free side: away from the atoms already attached
        ref = [nb for a in (a1, a2) for nb in a.neighbors if nb not in (a1, a2)]
        if ref:
            rx = sum(p.x for p in ref) / len(ref) - mx
            ry = sum(p.y for p in ref) / len(ref) - my
            if rx * nx + ry * ny > 0:
                nx, ny = -nx, -ny
        elif (x - mx) * nx + (y - my) * ny < 0:
            nx, ny = -nx, -ny

        apothem = length / (2 * tan(pi / n))
        cx, cy = mx + nx * apothem, my + ny * apothem

        radius = length / (2 * sin(pi / n))
        ang1 = atan2(a1.y - cy, a1.x - cx)
        ang2 = atan2(a2.y - cy, a2.x - cx)
        # walk around the ring in the a1 -> a2 direction
        sign = 1 if sin(ang2 - ang1) > 0 else -1
        step = sign * 2 * pi / n

        atoms = [a1, a2]
        for i in range(2, n):
            new_atom = mol.new_atom("C")
            new_atom.set_pos(cx + radius * cos(ang2 + (i - 1) * step),
                             cy + radius * sin(ang2 + (i - 1) * step))
            atoms.append(new_atom)
        return self._finish_ring(mol, atoms, name,
                                 symbol_offset=1 + n // 2, fixed=(bond,))

    # ----------------------------------------------------------------- internals
    def _finish_ring(self, mol, atoms, name, symbol_offset=0, fixed=()):
        """ creates the missing bonds, the double bond pattern and the labels """
        n, symbols, aromatic, no_double = self.spec(name)
        for i, symbol in symbols.items():
            atoms[(i + symbol_offset) % n].set_symbol(symbol)
        bonds = self._close_ring(mol, atoms)
        if aromatic:
            blocked = [atoms[(i + symbol_offset) % n] for i in no_double]
            self._kekulize(atoms, bonds, blocked=blocked, fixed=fixed)
        if name == "naphthalene" and len(bonds) > 3:
            # second ring on the bond opposite to where the first one is attached
            far = bonds[len(bonds) // 2]
            self._fuse_to_bond(far, *far.get_center(), name="benzene")
        mol.handle_overlap()
        mol.draw()
        return atoms, bonds

    @staticmethod
    def _close_ring(mol, atoms):
        """ bonds[i] joins atoms[i] and atoms[i+1], reusing bonds that exist """
        bonds = []
        n = len(atoms)
        for i in range(n):
            a1, a2 = atoms[i], atoms[(i + 1) % n]
            bond = mol.get_bond_between_atoms(a1, a2)
            if not bond:
                bond = mol.new_bond()
                bond.set_type("single")
                bond.connect_atoms(a1, a2)
            bonds.append(bond)
        return bonds

    @staticmethod
    def _kekulize(atoms, bonds, blocked=(), fixed=()):
        """ alternating double bonds, skipping atoms that are already saturated.

        Both alternations are tried and the one placing more double bonds wins,
        so a ring fused to an aromatic system keeps a sensible Kekule structure
        instead of overfilling the shared atoms. """
        n = len(bonds)
        busy = set(blocked)
        for atom in atoms:
            for b in atom.neighbor_edges:
                if b not in bonds and b.type in ("double", "triple"):
                    busy.add(atom)
        for b in fixed:
            if b.type == "double":
                busy.update(b.atoms)

        best = []
        for phase in (0, 1):
            used, plan = set(busy), []
            for i in range(n):
                bond, a1, a2 = bonds[i], atoms[i], atoms[(i + 1) % n]
                if bond in fixed or i % 2 != phase:
                    continue
                if a1 in used or a2 in used:
                    continue
                plan.append(bond)
                used.update((a1, a2))
            if len(plan) > len(best):
                best = plan

        center = (sum(a.x for a in atoms) / n, sum(a.y for a in atoms) / n)
        for bond in bonds:
            if bond in fixed:
                continue
            bond.set_type("double" if bond in best else "single")
            if bond in best:
                # the inner line of the double bond points into the ring
                v1, v2 = bond.atoms
                bond.auto_second_line_side = False
                bond.second_line_side = geo.line_get_side_of_point(
                    [v1.x, v1.y, v2.x, v2.y], center)

    @staticmethod
    def _free_direction(atom, x, y):
        """ direction with the most room around atom, the ring centre goes there """
        from math import atan2, pi
        neighbors = atom.neighbors
        if not neighbors:
            if (x, y) == (atom.x, atom.y):
                return -pi / 2
            return atan2(y - atom.y, x - atom.x)
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

class ArrowTool(Tool):
    """ Draws arrows. A plain click drops an arrow of the standard length,
    dragging sets length and direction (snapped to 15 degree steps). """

    def __init__(self, arrow_type="reaction", curvature=1.0):
        self.arrow_type = arrow_type
        self.curvature = curvature
        self.p1 = None
        self.preview_item = None

    def clear(self):
        if self.preview_item:
            App.paper.removeItem(self.preview_item)
            self.preview_item = None
        self.p1 = None

    def _end_point(self, x, y):
        from math import atan2, cos, sin, pi, hypot
        length = hypot(x - self.p1[0], y - self.p1[1])
        if length < Settings.min_arrow_length:
            # too short to be meaningful, use a standard horizontal arrow
            return self.p1[0] + 2 * Settings.min_arrow_length, self.p1[1]
        if QApplication.keyboardModifiers() & Qt.KeyboardModifier.ShiftModifier:
            return x, y
        step = pi / 12  # 15 degrees
        angle = round(atan2(y - self.p1[1], x - self.p1[0]) / step) * step
        return self.p1[0] + length * cos(angle), self.p1[1] + length * sin(angle)

    def on_mouse_press(self, x, y):
        self.p1 = (x, y)

    def on_mouse_move(self, x, y):
        if not App.paper.dragging or not self.p1: return
        if self.preview_item:
            App.paper.removeItem(self.preview_item)
        self.preview_item = App.paper.addLine(self.p1 + self._end_point(x, y),
                                              style=PenStyle.dashed)

    def on_mouse_release(self, x, y):
        if self.p1:
            if self.preview_item:
                App.paper.removeItem(self.preview_item)
                self.preview_item = None

            from .arrow import Arrow
            arrow = Arrow(self.p1, self._end_point(x, y),
                          type=self.arrow_type, curvature=self.curvature)
            App.paper.addObject(arrow)
            arrow.draw()
            App.paper.save_state_to_undo_stack()
        self.p1 = None

class TextTool(Tool):
    def __init__(self):
        self.current_text = None

    def on_mouse_press(self, x, y):
        focused = App.paper.focused_obj
        if isinstance(focused, TextLabel):
            self.current_text = focused
            App.paper.changeFocusTo(focused)
            App.paper.locked_focus_obj = focused
        else:
            # Create new text with empty string
            text = TextLabel(x, y, "")
            App.paper.addObject(text)
            text.draw()
            self.current_text = text
            App.paper.changeFocusTo(text)
            App.paper.locked_focus_obj = text
        App.paper.save_state_to_undo_stack()

    def on_key_press(self, key, text):
        return False

class SelectTool(Tool):
    def __init__(self):
        self.mouse_press_pos = None
        self.objs = []            # Selected objects (for display/delete)
        self._move_targets = []   # Objects to move when dragging
        self.dragging_curvature = None
        self.selection_rect_item = None
        self.selection_start_pos = None
        self.mode = "move"  # "move" or "select"
        self.lasso_points = []
        self.lasso_item = None
        self.rotating = False
        self.right_click_press = None
        self.is_copying = False
        self.changed = False

    def on_mouse_press(self, x, y):
        self.mouse_press_pos = (x, y)
        self.selection_start_pos = (x, y)
        self.lasso_points = [(x, y)]
        self.changed = False
        focused = App.paper.focused_obj
        if focused:
            self.mode = "move"

            # If we click on something already selected (or its parent molecule is selected),
            # move the current selection as a group.
            is_selected = focused in self.objs or (hasattr(focused, 'molecule') and focused.molecule in self.objs)

            # Narrowing will happen on mouse release if no dragging occurred.
            modifiers = QApplication.keyboardModifiers()
            has_modifiers = bool(modifiers & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier))
            
            if is_selected and not has_modifiers:
                # Keep current selection (moving as a group)
                if not self._move_targets:
                    self._move_targets = list(self.objs)
            elif is_selected:
                # Keep current selection (moving as a group)
                if not self._move_targets:
                    self._move_targets = list(self.objs)
            else:
                # If we click on something NOT selected, select the specific item visually
                self.objs = [focused]
                # BUT move the whole molecule by default for intuitive click-and-drag
                if hasattr(focused, 'molecule') and focused.molecule:
                    self._move_targets = [focused.molecule]
                else:
                    self._move_targets = [focused]
            
            # Update selection visuals
            # 1. Clear all selection visuals from top-level objects
            for o in App.paper.objects:
                if hasattr(o, 'set_selected'):
                    o.set_selected(False)
            
            # 2. Set selection for exactly what is in self.objs
            for o in self.objs:
                if hasattr(o, 'set_selected'):
                    o.set_selected(True)
            
            # Check for curvature dragging (if clicked near p2 of an arrow)
            if isinstance(focused, Arrow):
                from math import sqrt
                dist_to_tip = sqrt((x - focused.p2[0])**2 + (y - focused.p2[1])**2)
                if dist_to_tip < 15:
                    self.dragging_curvature = focused
        else:
            self.mode = "select"
            self.objs = []
            self._move_targets = []
            for o in App.paper.objects:
                if hasattr(o, 'set_selected'):
                    o.set_selected(False)

    def on_mouse_double_click(self, x, y):
        focused = App.paper.focused_obj
        if isinstance(focused, Atom):
            if App.window and hasattr(App.window, '_menu_edit_label'):
                App.window._menu_edit_label(focused)
                return
        if focused:
            # Double click selects the whole molecule or parent object
            if hasattr(focused, 'molecule') and focused.molecule:
                obj = focused.molecule
            elif hasattr(focused, 'parent') and focused.parent:
                obj = focused.parent
            else:
                obj = focused
            
            self.objs = [obj]
            self._move_targets = [obj]
            
            # Update selection visuals
            for o in App.paper.objects:
                if hasattr(o, 'set_selected'):
                    o.set_selected(False)
            for o in self.objs:
                if hasattr(o, 'set_selected'):
                    o.set_selected(True)

    def on_right_click(self, x, y):
        self.right_click_press = (x, y)
        focused = App.paper.focused_obj
        if focused:
            # Get the top-level object (molecule, arrow, text, etc.)
            if hasattr(focused, 'molecule') and focused.molecule:
                obj = focused.molecule
            elif hasattr(focused, 'parent') and focused.parent:
                obj = focused.parent
            else:
                obj = focused
            
            if obj not in self.objs:
                self.objs = [obj]
                self._move_targets = [obj]
            
            for o in App.paper.objects:
                if hasattr(o, 'set_selected'):
                    o.set_selected(o in self.objs)
            
            self.rotating = True

    def on_mouse_move(self, x, y):
        if not App.paper.dragging: return
        
        if self.rotating and self.right_click_press:
            # Rotate and Resize selected objects with right click
            if not self._move_targets: return
            obj = self._move_targets[0]
            if not hasattr(obj, 'get_center'): return
            
            center = obj.get_center()
            import math
            
            # Rotation
            angle1 = math.atan2(self.right_click_press[1] - center[1], self.right_click_press[0] - center[0])
            angle2 = math.atan2(y - center[1], x - center[0])
            angle = (angle2 - angle1) * 0.4 # Reduced sensitivity
            
            # Resizing
            dist1 = math.sqrt((self.right_click_press[0] - center[0])**2 + (self.right_click_press[1] - center[1])**2)
            dist2 = math.sqrt((x - center[0])**2 + (y - center[1])**2)
            
            scale_factor = 1.0
            if dist1 > 5:
                scale_factor = dist2 / dist1
            
            # Apply transformations
            if hasattr(obj, 'rotate'):
                obj.rotate(angle, center)
            if scale_factor != 1.0 and hasattr(obj, 'scale'):
                obj.scale(scale_factor, center)
            self.changed = True

            if hasattr(obj, 'atoms'):
                for a in obj.atoms:
                    a.on_bond_count_change()
                    a.draw()
                for b in obj.bonds:
                    b.draw()
            else:
                obj.draw()
            
            self.right_click_press = (x, y)
            return
        
        if self.dragging_curvature:
            arrow = self.dragging_curvature
            # Calculate curvature based on perpendicular distance from line p1-p2
            import math
            p1, p2 = arrow.p1, arrow.p2
            dist_12 = math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)
            if dist_12 < 1: return
            
            # Perpendicular distance using cross product
            # (p2-p1) cross (p_mouse-p1) / dist_12
            side_dist = ((p2[0] - p1[0]) * (y - p1[1]) - (p2[1] - p1[1]) * (x - p1[0])) / dist_12
            
            # Normalize curvature (dist * 0.3 * curvature = side_dist => curvature = side_dist / (dist * 0.3))
            arrow.curvature = side_dist / (dist_12 * 0.3)
            arrow.draw()
            self.changed = True
            return

        if self.mode == "select":
            # Draw selection rectangle (for now, rectangular selection)
            # Could be extended to lasso by tracking points
            if self.selection_rect_item:
                App.paper.removeItem(self.selection_rect_item)
            
            x1, y1 = self.selection_start_pos
            x2, y2 = x, y
            rect = [min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2)]
            self.selection_rect_item = App.paper.addRect(rect, color=(100, 100, 255), width=1, fill=(100, 100, 255, 50))
            return

        if not self._move_targets: return
        
        # Ctrl + Drag to copy behavior
        if not self.is_copying and (QApplication.keyboardModifiers() & Qt.KeyboardModifier.ControlModifier):
            # Create clones of current move targets
            new_targets = []
            new_objs = []
            seen_mols = set()
            
            for obj in self._move_targets:
                target_to_clone = None
                
                # Check if this object belongs to a molecule we've already decided to clone
                mol_id = None
                if hasattr(obj, 'molecule') and obj.molecule:
                    mol_id = id(obj.molecule)
                elif hasattr(obj, 'parent') and obj.parent and hasattr(obj.parent, 'atoms'):
                    mol_id = id(obj.parent)
                
                if mol_id and mol_id in seen_mols:
                    continue # Skip this object, it will be handled by the molecule clone
                
                if hasattr(obj, 'clone') and not hasattr(obj, 'molecule'): # Top-level (Molecule, Arrow, etc.)
                    target_to_clone = obj
                    if hasattr(obj, 'atoms'): # It's a molecule
                        seen_mols.add(id(obj))
                elif hasattr(obj, 'molecule') and obj.molecule:
                    target_to_clone = obj.molecule
                    seen_mols.add(id(obj.molecule))
                
                if target_to_clone:
                    clone = target_to_clone.clone()
                    App.paper.addObject(clone)
                    # We don't call clone.draw() here because move_by below will call it
                    new_targets.append(clone)
                    new_objs.append(clone)
                else:
                    # If it's a standalone object (or something we want to move but not clone)
                    # This check is slightly redundant now but safe
                    new_targets.append(obj)
            
            if new_objs:
                self._move_targets = new_targets
                self.objs = new_objs
                self.is_copying = True
                # NO: self.mouse_press_pos = (x, y) reset removed. 
                # Keeping the original press pos allows dx, dy to correctly offset 
                # the clones to the current mouse position in the first frame.
                # Select the new objects
                for o in App.paper.objects:
                    if hasattr(o, 'set_selected'):
                        o.set_selected(o in self.objs)
        
        # Move selected objects
        dx = x - self.mouse_press_pos[0]
        dy = y - self.mouse_press_pos[1]

        if dx == 0 and dy == 0: return
        self.changed = True

        for obj in self._move_targets:
            if hasattr(obj, 'move_by'):
                obj.move_by(dx, dy)
                # If it's a top-level object or a partial molecule selection (Atom), redraw
                if hasattr(obj, 'atoms'): # Molecule
                    obj.draw()
                else: # Atom, Arrow, TextLabel
                    obj.draw()
                    # For atoms, we must also redraw connected bonds
                    if hasattr(obj, 'neighbor_edges'):
                        for bond in obj.neighbor_edges:
                            bond.draw()
        
        self.mouse_press_pos = (x, y)

    def on_mouse_release(self, x, y):
        if self.rotating:
            rotated = self.changed
            self.rotating = False
            self.right_click_press = None
            if self.changed:
                App.paper.save_state_to_undo_stack()
            self.is_copying = False
            self.changed = False
            if not rotated:
                # a right click that did not turn into a drag opens the menu
                App.paper.context_menu_requested.emit(App.paper.focused_obj, x, y)
            return
        
        if self.selection_rect_item:
            App.paper.removeItem(self.selection_rect_item)
            self.selection_rect_item = None
            
            # Select individual objects within the rectangle.
            # Unlike double-click which selects the whole molecule,
            # rectangular selection picks only the atoms/bonds inside.
            x1, y1 = self.selection_start_pos
            x2, y2 = x, y
            rect = [min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2)]
            
            selected = App.paper.objectsInRect(rect)
            self.objs = []
            self._move_targets = []
            
            # Pre-collect atoms in rect for bond selection rule
            atoms_in_rect = {obj for obj in selected if isinstance(obj, Atom)}
            
            for obj in selected:
                if isinstance(obj, Bond):
                    # STRICTER RULE: Only select bond if both atoms are also in the rect
                    if obj.atoms[0] in atoms_in_rect and obj.atoms[1] in atoms_in_rect:
                        if obj not in self.objs:
                            self.objs.append(obj)
                else:
                    if obj not in self.objs:
                        self.objs.append(obj)
            
            # Identify if we've captured entire molecules or just parts
            potential_mols = {} # mol_id -> [selected_atoms_count, total_atoms_count, mol_obj]
            for obj in self.objs:
                if hasattr(obj, 'molecule') and obj.molecule:
                    mol = obj.molecule
                    if id(mol) not in potential_mols:
                        potential_mols[id(mol)] = [0, len(mol.atoms), mol]
                    if hasattr(obj, 'symbol'): # It's an atom
                        potential_mols[id(mol)][0] += 1
            
            # If an entire molecule is within the rectangle, use the molecule as move target.
            # Otherwise, use individual atoms.
            handled_mol_ids = set()
            for mol_id, (sel_count, total_count, mol) in potential_mols.items():
                if sel_count == total_count:
                    self._move_targets.append(mol)
                    handled_mol_ids.add(mol_id)
            
            for obj in self.objs:
                if isinstance(obj, Bond): continue # Bonds move via atoms
                if hasattr(obj, 'molecule') and obj.molecule:
                    if id(obj.molecule) not in handled_mol_ids:
                        if obj not in self._move_targets:
                            self._move_targets.append(obj)
                elif obj not in self._move_targets:
                    self._move_targets.append(obj)
            
            # Update selection visuals
            # 1. Clear all selection visuals from top-level objects
            for o in App.paper.objects:
                if hasattr(o, 'set_selected'):
                    o.set_selected(False)
            
            # 2. Set selection for exactly what is in self.objs
            for o in self.objs:
                if hasattr(o, 'set_selected'):
                    o.set_selected(True)
        
        # only a real modification belongs on the undo stack, a plain selection
        # click would otherwise flush the whole history
        if self.changed:
            App.paper.save_state_to_undo_stack()

        # Check if we should narrow the selection (clicked without dragging)
        if self.mode == "move" and not self.is_copying and not self.rotating and not self.selection_rect_item:
            import math
            dx = x - self.selection_start_pos[0]
            dy = y - self.selection_start_pos[1]
            if math.sqrt(dx*dx + dy*dy) < 3:
                focused = App.paper.focused_obj
                # Narrow selection if clicking a member of a group without modifiers
                modifiers = QApplication.keyboardModifiers()
                has_modifiers = bool(modifiers & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier))
                
                if focused and (focused in self.objs or (hasattr(focused, 'molecule') and focused.molecule in self.objs)) and len(self.objs) >= 1 and not has_modifiers:
                    # If the selected object is a molecule, narrowing it to the specific atom/bond
                    if len(self.objs) > 1 or (len(self.objs) == 1 and hasattr(self.objs[0], 'atoms')):
                        self.objs = [focused]
                        if isinstance(focused, Bond):
                            self._move_targets = [focused.atoms[0], focused.atoms[1]]
                        elif isinstance(focused, Atom):
                            self._move_targets = [focused]
                        else:
                            self._move_targets = [focused]
                        
                        for o in App.paper.objects:
                            if hasattr(o, 'set_selected'):
                                o.set_selected(False)
                        for o in self.objs:
                            if hasattr(o, 'set_selected'):
                                o.set_selected(True)
        
        self.dragging_curvature = None
        self.mode = "move"
        self.is_copying = False
        self.changed = False

    def on_key_press(self, key, text):
        if handle_atom_keypress(key, text):
            return True
        focused = App.paper.focused_obj
        if isinstance(focused, Arrow):
            if text in ('+', '='):
                focused.curvature += 0.2
                focused.draw()
                App.paper.save_state_to_undo_stack()
                return True
            elif text in ('-', '_'):
                focused.curvature -= 0.2
                focused.draw()
                App.paper.save_state_to_undo_stack()
                return True
        elif isinstance(focused, TextLabel):
            return False
        return False

class ChargeTool(Tool):
    """ Adds or removes a formal charge on the clicked atom. """

    def __init__(self, delta):
        Tool.__init__(self)
        self.delta = delta

    def on_mouse_press(self, x, y):
        focused = App.paper.focused_obj
        if not isinstance(focused, Atom):
            focused = App.paper.find_closest_atom(x, y, 12)
        if isinstance(focused, Atom):
            focused.set_charge(focused.charge + self.delta)
            focused.draw()
            App.paper.save_state_to_undo_stack("Charge")


class ShapeTool(Tool):
    def __init__(self, shape_class):
        super().__init__()
        self.shape_class = shape_class
        self.p1 = None
        self.current_shape = None

    def on_mouse_press(self, x, y):
        self.p1 = (x, y)
        self.current_shape = self.shape_class(x, y, x, y)
        App.paper.addObject(self.current_shape)
        self.current_shape.draw()

    def on_mouse_move(self, x, y):
        if App.paper.dragging and self.p1:
            x1, y1 = self.p1
            # Constraints (Shift for Square/Circle)
            if QApplication.keyboardModifiers() & Qt.KeyboardModifier.ShiftModifier:
                dx = x - x1
                dy = y - y1
                size = max(abs(dx), abs(dy))
                x = x1 + (size if dx > 0 else -size)
                y = y1 + (size if dy > 0 else -size)
            
            self.current_shape.x2 = x
            self.current_shape.y2 = y
            self.current_shape.draw()

    def on_mouse_release(self, x, y):
        if self.current_shape:
            bbox = self.current_shape.bounding_box()
            if (bbox[2] - bbox[0]) < 5 and (bbox[3] - bbox[1]) < 5:
                self.current_shape.delete_from_paper()
            else:
                App.paper.save_state_to_undo_stack(f"Add {self.current_shape.class_name}")
        self.p1 = None
        self.current_shape = None


class LonePairTool(Tool):
    """ Adds or removes lone pairs on the clicked atom. """
    def __init__(self):
        super().__init__()

    def on_mouse_press(self, x, y):
        focused = App.paper.focused_obj
        if not isinstance(focused, Atom):
            focused = App.paper.find_closest_atom(x, y, 12)
        if isinstance(focused, Atom):
            focused.lonepairs = (getattr(focused, "lonepairs", 0) + 1) % 4
            focused.draw()
            App.paper.save_state_to_undo_stack("Lone Pair")


class RadicalPlusTool(Tool):
    """ Toggles a radical dot with a plus sign on the clicked atom. """
    def __init__(self):
        super().__init__()

    def on_mouse_press(self, x, y):
        focused = App.paper.focused_obj
        if not isinstance(focused, Atom):
            focused = App.paper.find_closest_atom(x, y, 12)
        if isinstance(focused, Atom):
            focused.radical_plus = 1 - getattr(focused, "radical_plus", 0)
            if focused.radical_plus:
                focused.radical_minus = 0
            focused.draw()
            App.paper.save_state_to_undo_stack("Radical +")


class RadicalMinusTool(Tool):
    """ Toggles a radical dot with a minus sign on the clicked atom. """
    def __init__(self):
        super().__init__()

    def on_mouse_press(self, x, y):
        focused = App.paper.focused_obj
        if not isinstance(focused, Atom):
            focused = App.paper.find_closest_atom(x, y, 12)
        if isinstance(focused, Atom):
            focused.radical_minus = 1 - getattr(focused, "radical_minus", 0)
            if focused.radical_minus:
                focused.radical_plus = 0
            focused.draw()
            App.paper.save_state_to_undo_stack("Radical −")
