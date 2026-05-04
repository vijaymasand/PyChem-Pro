# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3)
from .app_data import App, Settings
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
    def __init__(self):
        Tool.__init__(self)
        self.reset()
        self.preview_item = None
        self.atom_with_preview_bond = None

    def reset(self):
        self.atom1 = None
        self.mouse_press_pos = None

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
            self.atom1 = mol.new_atom(toolsettings['structure'])
            self.atom1.set_pos(*self.mouse_press_pos)
            self.atom1.draw()
        
        if self.preview_item:
            App.paper.removeItem(self.preview_item)
        self.preview_item = App.paper.addLine(self.atom1.pos + (x, y), style=PenStyle.dashed)

    def on_mouse_release(self, x, y):
        if not App.paper.dragging:
            self.on_mouse_click(x, y)
            return
        
        self.clear_preview()
        if self.atom1:
            atom2 = self.atom1.molecule.new_atom(toolsettings['structure'])
            atom2.set_pos(x, y)
            bond = self.atom1.molecule.new_bond()
            bond.set_type(toolsettings['bond_type'])
            bond.connect_atoms(self.atom1, atom2)
            
            # Check for overlap
            focused_obj = App.paper.focused_obj
            if isinstance(focused_obj, Atom) and focused_obj is not self.atom1:
                focused_obj.eat_atom(atom2)
                atom2 = focused_obj
            
            self.atom1.draw()
            atom2.draw()
            bond.draw()
            App.paper.save_state_to_undo_stack()
        self.reset()

    def on_mouse_click(self, x, y):
        focused_obj = App.paper.focused_obj
        if not focused_obj:
            mol = Molecule()
            App.paper.addObject(mol)
            atom = mol.new_atom(toolsettings['structure'])
            atom.set_pos(x, y)
            atom.draw()
        elif isinstance(focused_obj, Atom):
            # Change symbol or add branch
            if focused_obj.symbol != toolsettings['structure']:
                focused_obj.set_symbol(toolsettings['structure'])
                focused_obj.draw()
            else:
                # Add branch
                bond_length = Settings.bond_length * focused_obj.molecule.scale_val
                nx, ny = focused_obj.molecule.find_place(focused_obj, bond_length)
                atom2 = focused_obj.molecule.new_atom(toolsettings['structure'])
                atom2.set_pos(nx, ny)
                bond = focused_obj.molecule.new_bond()
                bond.set_type(toolsettings['bond_type'])
                bond.connect_atoms(focused_obj, atom2)
                focused_obj.draw()
                atom2.draw()
                bond.draw()
        elif isinstance(focused_obj, Bond):
            # Cycle bond type
            types = ['single', 'double', 'triple']
            if focused_obj.type in types:
                idx = (types.index(focused_obj.type) + 1) % len(types)
                focused_obj.set_type(types[idx])
                focused_obj.draw()
                for a in focused_obj.atoms: a.draw()
        
        App.paper.save_state_to_undo_stack()
        self.reset()

    def show_preview(self, atom):
        bond_length = Settings.bond_length * atom.molecule.scale_val
        nx, ny = atom.molecule.find_place(atom, bond_length)
        self.preview_item = App.paper.addLine(atom.pos + (nx, ny), style=PenStyle.dashed)
        self.atom_with_preview_bond = atom

    def on_key_press(self, key, text):
        focused = App.paper.focused_obj
        if isinstance(focused, Atom):
            mapping = {
                'c': 'C', 'n': 'N', 'o': 'O', 's': 'S', 'p': 'P', 
                'f': 'F', 'l': 'Cl', 'b': 'Br', 'i': 'I', 'h': 'H'
            }
            symbol = mapping.get(text.lower())
            if symbol:
                focused.set_symbol(symbol)
                focused.draw()
                App.paper.save_state_to_undo_stack()
                App.paper.locked_focus_obj = None
                return True
        return False

    def clear_preview(self):
        if self.preview_item:
            App.paper.removeItem(self.preview_item)
            self.preview_item = None
            self.atom_with_preview_bond = None

class EraserTool(Tool):
    def on_mouse_press(self, x, y):
        # Refresh focus at the moment of click for maximum precision
        zoom = App.paper.view.transform().m11()
        d = max(2, 8 / zoom)
        focused = App.paper.find_closest_object(x, y, d)
        
        if focused:
            # Save undo state BEFORE deletion
            App.paper.save_state_to_undo_stack()
            
            from .atom import Atom
            from .bond import Bond
            from .arrow import Arrow
            from .text_label import TextLabel

            if isinstance(focused, Atom):
                mol = focused.molecule
                # Remove connected bonds first
                for bond in list(focused.bonds):
                    bond.disconnect_atoms()
                    mol.remove_bond(bond)
                    bond.delete_from_paper()
                mol.remove_atom(focused)
                focused.delete_from_paper()
                # Split or delete molecule if empty
                if not mol.atoms:
                    App.paper.removeObject(mol)
                else:
                    mol.split_fragments()
            elif isinstance(focused, Bond):
                mol = focused.molecule
                focused.disconnect_atoms()
                mol.remove_bond(focused)
                focused.delete_from_paper()
                mol.split_fragments()
            elif isinstance(focused, Arrow):
                focused.delete_from_paper()
            elif isinstance(focused, TextLabel):
                focused.delete_from_paper()
            
            # Draw all dirty objects
            App.paper.redraw_dirty_objects()

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
    def __init__(self, template_name):
        self.template_name = template_name

    def on_mouse_press(self, x, y):
        # Look for a bond in the vicinity to fuse to, even if an atom is focused
        objs = App.paper.objectsInRect([x - 7, y - 7, x + 7, y + 7])
        bond = None
        for obj in objs:
            if isinstance(obj, Bond):
                bond = obj
                break
        
        if bond:
            self._fuse_to_bond(bond, x, y)
        else:
            focused = App.paper.focused_obj
            if isinstance(focused, Bond):
                self._fuse_to_bond(focused, x, y)
            else:
                self._place_at(x, y)
        App.paper.save_state_to_undo_stack()

    def _place_at(self, x, y):
        mol = Molecule()
        App.paper.addObject(mol)
        self._generate_ring(mol, x, y, 30, 0)

    def _fuse_to_bond(self, bond, x, y):
        mol = bond.molecule
        a1, a2 = bond.atoms
        # Calculate angle and length of the bond
        from math import atan2, sqrt, pi, cos, sin
        angle = atan2(a2.y - a1.y, a2.x - a1.x)
        length = sqrt((a2.x - a1.x)**2 + (a2.y - a1.y)**2)
        
        n = {"benzene": 6, "cyclohexane": 6, "cyclopentane": 5, "pyridine": 6, "furan": 5, "pyrrole": 5}.get(self.template_name, 6)
        ext_angle = 2 * pi / n
        
        # Determine side based on click position
        # Cross product of (a2-a1) and (click-a1)
        # In Qt coords (Y down), a positive cross product means one side, negative the other.
        side = (a2.x - a1.x) * (y - a1.y) - (a2.y - a1.y) * (x - a1.x)
        if side > 0:
            ext_angle = -ext_angle
        
        # We start from a2 and go around. The first bond is a2->a1 (already exists)
        # So we need to add n-2 new atoms.
        curr_angle = angle + pi - ext_angle
        curr_atom = a1
        
        atoms = [a2, a1]
        for i in range(n - 2):
            nx = curr_atom.x + length * cos(curr_angle)
            ny = curr_atom.y + length * sin(curr_angle)
            
            # Check if an atom already exists at this position
            new_atom = mol.new_atom("C")
            new_atom.set_pos(nx, ny)
            atoms.append(new_atom)
            
            new_bond = mol.new_bond()
            new_bond.set_type("single")
            new_bond.connect_atoms(curr_atom, new_atom)
            
            curr_atom = new_atom
            curr_angle -= ext_angle
            
        # Close the ring
        final_bond = mol.new_bond()
        final_bond.connect_atoms(curr_atom, a2)
        
        # If it's benzene or pyridine, we need to fix double bond pattern and side
        if self.template_name in ("benzene", "pyridine"):
            # Naphthalene-like pattern: shared bond can be double or single.
            # Let's make the new ring alternate starting from a2-curr_atom as double.
            # And set second_line_side to point inwards.
            
            # Calculate center of new ring for second_line_side
            new_center_x = sum(a.x for a in atoms) / n
            new_center_y = sum(a.y for a in atoms) / n
            
            # Bonds in order: [a2-a1, a1-new1, new1-new2, new2-new3, new3-new4, new4-a2]
            shared_bond = bond
            # List of bonds for the new ring
            ring_bonds = [shared_bond]
            for i in range(2, len(atoms)):
                ring_bonds.append(mol.get_bond_between_atoms(atoms[i-1], atoms[i]))
            ring_bonds.append(final_bond)

            # Check if shared atoms already have double bonds elsewhere
            a1_has_double = any(bo.type == "double" for bo in a1.neighbor_edges if bo != bond)
            a2_has_double = any(bo.type == "double" for bo in a2.neighbor_edges if bo != bond)
            
            # Count total double bonds on shared atoms (including the shared bond)
            a1_double_count = sum(1 for bo in a1.neighbor_edges if bo.type == "double")
            a2_double_count = sum(1 for bo in a2.neighbor_edges if bo.type == "double")
            
            for i, b in enumerate(ring_bonds):
                if not b or i == 0: continue # Don't force change the shared bond
                
                # Extended logic to support up to 10 rings fusion
                # Check if adding a double bond would exceed valency (max 4 for carbon)
                # Each atom can have at most 2 double bonds (since each double bond counts as 2)
                can_add_double = True
                if b.atoms[0] in (a1, a2):
                    atom = b.atoms[0]
                    if atom == a1 and a1_double_count >= 2:
                        can_add_double = False
                    elif atom == a2 and a2_double_count >= 2:
                        can_add_double = False
                if b.atoms[1] in (a1, a2):
                    atom = b.atoms[1]
                    if atom == a1 and a1_double_count >= 2:
                        can_add_double = False
                    elif atom == a2 and a2_double_count >= 2:
                        can_add_double = False
                
                # If neither atom has a double bond (and shared isn't double), we can add 3 double bonds (i=1,3,5)
                # If they do have double bonds but still have capacity, we can still add double bonds
                # For extensive ring systems, prefer single bonds when valency is tight
                if not a1_has_double and not a2_has_double and bond.type != "double":
                    new_type = "double" if (i % 2 == 1) else "single"
                elif can_add_double:
                    # Allow double bond placement based on alternating pattern
                    new_type = "double" if (i % 2 == 0) else "single"
                else:
                    # Must use single bond to avoid exceeding valency
                    new_type = "single"
                
                b.set_type(new_type)
                
                if new_type == "double":
                    # Point double bond line towards the new ring center using native geometry logic
                    b.auto_second_line_side = False
                    v1 = b.atoms[0]
                    v2 = b.atoms[1]
                    from . import geometry as geo
                    b.second_line_side = geo.line_get_side_of_point([v1.x, v1.y, v2.x, v2.y], (new_center_x, new_center_y))
        else:
            final_bond.set_type("single")
        
        mol.handle_overlap()
        for a in mol.atoms: a.draw()
        for b in mol.bonds: b.draw()

    def _generate_ring(self, mol, x, y, radius, start_angle):
        from math import cos, sin, pi
        n = {"benzene": 6, "cyclohexane": 6, "cyclopentane": 5, "pyridine": 6, "furan": 5, "pyrrole": 5}.get(self.template_name, 6)
        points = []
        for i in range(n):
            angle = start_angle + i * (2 * pi / n)
            points.append((x + radius * cos(angle), y + radius * sin(angle)))
        
        atoms = []
        for i, p in enumerate(points):
            symbol = "C"
            if self.template_name == "pyridine" and i == 0: symbol = "N"
            elif self.template_name == "furan" and i == 0: symbol = "O"
            elif self.template_name == "pyrrole" and i == 0: symbol = "N"
            a = mol.new_atom(symbol)
            a.set_pos(*p)
            atoms.append(a)
        
        for i in range(n):
            b = mol.new_bond()
            is_double = (self.template_name in ("benzene", "pyridine") and i % 2 == 0)
            b.set_type("double" if is_double else "single")
            b.connect_atoms(atoms[i], atoms[(i+1)%n])
        
        for a in mol.atoms: a.draw()
        for b in mol.bonds: b.draw()

class ArrowTool(Tool):
    def __init__(self, arrow_type="reaction", curvature=1.0):
        self.arrow_type = arrow_type
        self.curvature = curvature
        self.p1 = None
        self.preview_item = None

    def on_mouse_press(self, x, y):
        self.p1 = (x, y)

    def on_mouse_move(self, x, y):
        if not App.paper.dragging or not self.p1: return
        if self.preview_item:
            App.paper.removeItem(self.preview_item)
        self.preview_item = App.paper.addLine(self.p1 + (x, y), style=PenStyle.dashed)

    def on_mouse_release(self, x, y):
        if self.p1:
            if self.preview_item:
                App.paper.removeItem(self.preview_item)
                self.preview_item = None
            
            from .arrow import Arrow
            arrow = Arrow(self.p1, (x, y), type=self.arrow_type, curvature=self.curvature)
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

    def on_mouse_press(self, x, y):
        self.mouse_press_pos = (x, y)
        self.selection_start_pos = (x, y)
        self.lasso_points = [(x, y)]
        focused = App.paper.focused_obj
        if focused:
            self.mode = "move"
            
            # If we click on something already selected (or its parent molecule is selected),
            # move the current selection as a group.
            is_selected = focused in self.objs or (hasattr(focused, 'molecule') and focused.molecule in self.objs)
            
            # Narrow selection if clicking a member of a group without modifiers
            modifiers = QApplication.keyboardModifiers()
            has_modifiers = bool(modifiers & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier))
            
            if is_selected and not has_modifiers:
                # Narrow to just the focused object
                self.objs = [focused]
                if isinstance(focused, Bond):
                    self._move_targets = [focused.atoms[0], focused.atoms[1]]
                elif isinstance(focused, Atom):
                    self._move_targets = [focused]
                else:
                    self._move_targets = [focused]
            elif is_selected:
                # Keep current selection (moving as a group)
                self._move_targets = list(self.objs)
                # If a top-level molecule is selected, it's already in _move_targets
                # If only individual parts are selected, they are already in _move_targets
            else:
                # If we click on something NOT selected, select it
                self.objs = [focused]
                if isinstance(focused, Bond):
                    # For bonds, we move the two connected atoms to allow stretching
                    self._move_targets = [focused.atoms[0], focused.atoms[1]]
                elif isinstance(focused, Atom):
                    # For atoms, we move only the atom itself
                    self._move_targets = [focused]
                elif hasattr(focused, 'molecule') and focused.molecule:
                    # Fallback for other molecule parts (if any)
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
                    o.set_selected(o in self.objs)

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
            self.rotating = False
            self.right_click_press = None
            if self._move_targets:
                App.paper.save_state_to_undo_stack()
            self.is_copying = False
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
        
        if self.objs or self.dragging_curvature:
            App.paper.save_state_to_undo_stack()
        self.dragging_curvature = None
        self.mode = "move"
        self.is_copying = False

    def on_key_press(self, key, text):
        focused = App.paper.focused_obj
        if isinstance(focused, Atom):
            mapping = {
                'c': 'C', 'n': 'N', 'o': 'O', 's': 'S', 'p': 'P', 
                'f': 'F', 'l': 'Cl', 'b': 'Br', 'i': 'I', 'h': 'H'
            }
            symbol = mapping.get(text.lower())
            if symbol:
                focused.set_symbol(symbol)
                focused.draw()
                App.paper.save_state_to_undo_stack()
                App.paper.locked_focus_obj = None
                return True
        elif isinstance(focused, Arrow):
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
