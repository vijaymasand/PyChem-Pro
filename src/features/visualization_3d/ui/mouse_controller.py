"""
Mouse interaction controller for MolViewer3D.

Handles rotation, pan, zoom, atom picking, rubber-band selection,
lasso selection, double-click ligand selection, measurement picks,
and auto-rotate.
"""

import math
from src.shared.qt_compat import Qt, QPointF, QRectF, QWheelEvent


class MouseController:
    """
    Processes mouse / key events on behalf of a MolViewer3D widget.

    The controller reads and mutates camera state (``rot_x``, ``rot_y``,
    ``pan_x``, ``pan_y``, ``zoom``) and selection state
    (``selected_atoms``, ``_measure_atoms``, ``_measurements``) directly
    on the *viewer* reference it is given at construction time.
    """

    def __init__(self, viewer):
        self._v = viewer  # MolViewer3D instance

    # ─── Event handlers ─────────────────────────────────────────

    def handle_mouse_press(self, event):
        """Handle mouse press: start rotation, pan, rubber-band or lasso selection."""
        v = self._v
        v._last_mouse_pos = event.position()
        v._mouse_button = event.button()
        v._has_dragged = False

        modifiers = event.modifiers()

        # Ctrl + left-click begins lasso selection
        if (event.button() == Qt.MouseButton.LeftButton
                and modifiers & Qt.KeyboardModifier.ControlModifier):
            v._is_lasso = True
            v._lasso_path = [QPointF(event.position())]
            v._is_selecting = False   # not rubber-band
            return

        # Shift + left-click begins rubber-band selection (PyMOL-style)
        if (event.button() == Qt.MouseButton.LeftButton
                and modifiers & Qt.KeyboardModifier.ShiftModifier):
            v._is_selecting = True
            v._sel_rect_origin = event.position()
            v._sel_rect_end = event.position()
            v._is_lasso = False
            return

        # Normal press: cancel any pending lasso/rubber-band state
        v._is_lasso = False
        v._is_selecting = False

    def handle_mouse_move(self, event):
        """Handle mouse move: rotate, pan, update rubber-band / lasso, or hover."""
        v = self._v
        if v._last_mouse_pos is None:
            self._detect_hover(event.position())
            return

        # ── Lasso drag ────────────────────────────────────────────
        if getattr(v, '_is_lasso', False):
            v._lasso_path.append(QPointF(event.position()))
            v._has_dragged = True
            v.update()
            return

        # ── Rubber-band selection drag ────────────────────────────
        if v._is_selecting:
            v._sel_rect_end = event.position()
            v.update()
            return

        # ── Camera drag ───────────────────────────────────────────
        dx = event.position().x() - v._last_mouse_pos.x()
        dy = event.position().y() - v._last_mouse_pos.y()

        if v._mouse_button == Qt.MouseButton.LeftButton:
            v.rot_y += dx * 0.5
            v.rot_x += dy * 0.5
            if abs(dx) > 1 or abs(dy) > 1:
                v._has_dragged = True
                v._is_interacting = True
        elif v._mouse_button == Qt.MouseButton.RightButton:
            cx = v.width() / 2.0
            cy = v.height() / 2.0
            ang1 = math.atan2(v._last_mouse_pos.y() - cy, v._last_mouse_pos.x() - cx)
            ang2 = math.atan2(event.position().y() - cy, event.position().x() - cx)
            angle_diff = ang2 - ang1
            if angle_diff > math.pi: angle_diff -= 2 * math.pi
            elif angle_diff < -math.pi: angle_diff += 2 * math.pi
            if not hasattr(v, 'rot_z'): v.rot_z = 0.0
            v.rot_z -= math.degrees(angle_diff)
            if abs(dx) > 1 or abs(dy) > 1:
                v._has_dragged = True
                v._is_interacting = True
        elif v._mouse_button == Qt.MouseButton.MiddleButton:
            v.pan_x += dx
            v.pan_y += dy
            if abs(dx) > 1 or abs(dy) > 1:
                v._has_dragged = True
                v._is_interacting = True

        v._last_mouse_pos = event.position()
        v.update()

    def handle_mouse_release(self, event):
        """Handle mouse release: commit lasso / selection rect, clicks, measurements."""
        v = self._v

        # ── Finish lasso selection ────────────────────────────────
        if getattr(v, '_is_lasso', False) and event.button() == Qt.MouseButton.LeftButton:
            self._commit_lasso_selection()
            v._is_lasso = False
            v._lasso_path = []
            v._last_mouse_pos = None
            v._mouse_button = None
            v.update()
            return

        # ── Finish rubber-band selection ──────────────────────────
        if v._is_selecting and event.button() == Qt.MouseButton.LeftButton:
            self._commit_rubber_band_selection()
            v._is_selecting = False
            v._sel_rect_origin = None
            v._sel_rect_end = None
            v._last_mouse_pos = None
            v._mouse_button = None
            v.update()
            return

        # ── Single click handling ─────────────────────────────────
        was_click = False
        if v._last_mouse_pos is not None:
            moved = (abs(event.position().x() - v._last_mouse_pos.x()) +
                     abs(event.position().y() - v._last_mouse_pos.y()))
            was_click = not getattr(v, '_has_dragged', False) and moved < 3

        btn = event.button()

        if was_click and btn == Qt.MouseButton.LeftButton:
            atom_idx = self._hit_test(event.position())
            if atom_idx != -1:
                modifiers = event.modifiers()
                if modifiers & Qt.KeyboardModifier.ShiftModifier:
                    if atom_idx in v.selected_atoms:
                        v.selected_atoms.remove(atom_idx)
                        if hasattr(v, '_measure_atoms') and atom_idx in v._measure_atoms:
                            v._measure_atoms.remove(atom_idx)
                    else:
                        v.selected_atoms.add(atom_idx)
                        if not hasattr(v, '_measure_atoms'):
                            v._measure_atoms = []
                        v._measure_atoms.append(atom_idx)
                else:
                    v.selected_atoms = {atom_idx}
                    if not hasattr(v, '_measure_atoms'):
                        v._measure_atoms = []
                    v._measure_atoms.append(atom_idx)

                    if len(v._measure_atoms) > 3:
                        v._measure_atoms.pop(0)

                if hasattr(v, '_measure_atoms'):
                    if len(v._measure_atoms) == 2:
                        self._complete_distance_measurement()
                    elif len(v._measure_atoms) == 3:
                        self._complete_angle_measurement()

                v.selection_changed.emit(set(v.selected_atoms))
                v.atom_clicked.emit(atom_idx)
                v.update()
            else:
                if v.selected_atoms:
                    v.selected_atoms.clear()
                if hasattr(v, '_measure_atoms'):
                    v._measure_atoms.clear()
                if hasattr(v, '_measurements'):
                    v._measurements.clear()
                v.selection_changed.emit(set())
                v.update()

        v._last_mouse_pos = None
        v._mouse_button = None
        if getattr(v, '_is_interacting', False):
            v._is_interacting = False
            v.update()

    def handle_double_click(self, event):
        """
        Handle double-click: select the entire ligand/fragment that the clicked atom belongs to.

        For PDB molecules: selects all HETATM atoms with the same res_name + res_seq + chain_id.
        For small molecules (SDF/MOL2/MOL) or non-protein PDB: selects all atoms in the
        same connected fragment as the clicked atom.
        """
        v = self._v
        if event.button() != Qt.MouseButton.LeftButton:
            return

        atom_idx = self._hit_test(event.position())
        if atom_idx == -1:
            return

        ligand_atoms = self._get_ligand_atoms(atom_idx)
        if ligand_atoms:
            v.selected_atoms = set(ligand_atoms)
            v.selection_changed.emit(set(v.selected_atoms))
            v.update()

    def _get_ligand_atoms(self, atom_idx):
        """
        Return the set of atom indices belonging to the same ligand as atom_idx.

        Strategy:
          1. If the molecule is a protein (is_protein property) and the atom is a HETATM,
             group by (res_name, res_seq, chain_id) — this captures multi-residue ligands
             within a PDB file.
          2. Otherwise (SDF/MOL2/MOL/non-HETATM click), return the entire connected
             fragment that contains atom_idx via BFS on the bond graph.
        """
        v = self._v
        if not v.molecule or atom_idx >= len(v.molecule.atoms):
            return set()

        atom = v.molecule.atoms[atom_idx]
        is_protein = getattr(v.molecule, 'properties', {}).get('is_protein', False)

        if is_protein and getattr(atom, 'is_hetatm', False):
            # Group by residue identity within HETATM records
            target_res = (
                getattr(atom, 'res_name', None),
                getattr(atom, 'res_seq', None),
                getattr(atom, 'chain_id', None),
            )
            # Exclude water molecules
            if target_res[0] and target_res[0].upper() in ('HOH', 'WAT', 'SOL', 'DOD'):
                return {atom_idx}

            ligand_indices = set()
            for a in v.molecule.atoms:
                if getattr(a, 'is_hetatm', False):
                    res = (
                        getattr(a, 'res_name', None),
                        getattr(a, 'res_seq', None),
                        getattr(a, 'chain_id', None),
                    )
                    if res == target_res:
                        ligand_indices.add(a.index)
            return ligand_indices

        else:
            # BFS to find all atoms in the same connected fragment
            adjacency = v.molecule._adjacency
            visited = set()
            queue = [atom_idx]
            while queue:
                current = queue.pop()
                if current in visited:
                    continue
                visited.add(current)
                for neighbor, _ in adjacency.get(current, []):
                    if neighbor not in visited:
                        queue.append(neighbor)
            return visited

    def handle_wheel(self, event: QWheelEvent):
        v = self._v
        delta = event.angleDelta().y()
        factor = 1.1 if delta > 0 else 0.9
        v.zoom *= factor
        v.zoom = max(5, min(200, v.zoom))
        v._is_interacting = True
        v.update()
        # Reset after a short delay since wheel events don't have a release event
        from PySide6.QtCore import QTimer
        QTimer.singleShot(200, lambda: self._reset_interaction(v))

    def _reset_interaction(self, v):
        if getattr(v, '_is_interacting', False):
            v._is_interacting = False
            v.update()

    def handle_key_press(self, event):
        """Handle keyboard shortcuts -- Delete key removes selected atoms."""
        v = self._v
        if event.key() == Qt.Key.Key_Delete and v.selected_atoms:
            v.delete_requested.emit(set(v.selected_atoms))
            return True  # handled
        return False  # not handled, caller should call super()

    # ─── Auto-rotate ─────────────────────────────────────────────

    def auto_rotate_step(self):
        self._v.rot_y += 0.8
        self._v.update()

    # ─── Rubber-band Helpers ─────────────────────────────────────

    def _commit_rubber_band_selection(self):
        """
        Finalise a Shift+drag selection: find all projected atoms whose
        screen positions fall inside the rubber-band rectangle and add
        them to ``selected_atoms``.  Emits ``selection_changed``.
        """
        v = self._v
        if not v._sel_rect_origin or not v._sel_rect_end:
            return

        # Build normalised rectangle
        x1 = min(v._sel_rect_origin.x(), v._sel_rect_end.x())
        y1 = min(v._sel_rect_origin.y(), v._sel_rect_end.y())
        x2 = max(v._sel_rect_origin.x(), v._sel_rect_end.x())
        y2 = max(v._sel_rect_origin.y(), v._sel_rect_end.y())

        # Tiny rectangle treated as a missed click -- clear selection
        if (x2 - x1) < 4 and (y2 - y1) < 4:
            v.selected_atoms.clear()
            v.selection_changed.emit(set())
            return

        rect = QRectF(x1, y1, x2 - x1, y2 - y1)
        projected = v._renderer._project_atoms(v)
        newly_selected = set()
        for atom_idx, sx, sy, sz, radius, color in projected:
            if rect.contains(QPointF(sx, sy)):
                newly_selected.add(atom_idx)

        if newly_selected:
            v.selected_atoms |= newly_selected
        v.selection_changed.emit(set(v.selected_atoms))

    # ─── Lasso Helpers ───────────────────────────────────────────

    def _commit_lasso_selection(self):
        """
        Finalise a Ctrl+drag lasso: find all projected atoms whose screen
        positions fall inside the freehand polygon and add them to
        ``selected_atoms``.  Emits ``selection_changed``.
        """
        v = self._v
        lasso = getattr(v, '_lasso_path', [])
        if len(lasso) < 3:
            # Too small — treat as a click
            v.selected_atoms.clear()
            v.selection_changed.emit(set())
            return

        from src.shared.qt_compat import QPainterPath, QPointF as _QPointF
        path = QPainterPath()
        path.moveTo(lasso[0])
        for pt in lasso[1:]:
            path.lineTo(pt)
        path.closeSubpath()

        projected = v._renderer._project_atoms(v)
        newly_selected = set()
        for atom_idx, sx, sy, sz, radius, color in projected:
            if path.contains(_QPointF(sx, sy)):
                newly_selected.add(atom_idx)

        if newly_selected:
            v.selected_atoms |= newly_selected
        v.selection_changed.emit(set(v.selected_atoms))

    # ─── Hit-testing / Hover ─────────────────────────────────────

    def _detect_hover(self, pos):
        v = self._v
        atom_idx = self._hit_test(pos)
        if atom_idx != v._hovered_atom:
            v._hovered_atom = atom_idx
            v.atom_hovered.emit(atom_idx)
            v.update()

    def _hit_test(self, pos):
        """
        Return the atom index closest to *pos* (screen coordinates).

        Improved accuracy:
        - Candidates must be within 1.5× their projected sphere radius of the cursor.
        - Among candidates, the one with the **smallest screen-space distance to centre**
          wins (not the one closest to the camera).  Depth is used only as a tiebreaker
          for atoms that are essentially at the same screen position.
        - This ensures ligand atoms on top of (but smaller than) cartoon ribbon atoms
          are always preferred when the cursor is visually over them.
        """
        v = self._v
        projected = v._renderer._project_atoms(v)
        if not projected:
            return -1

        mx = pos.x()
        my = pos.y()

        best_idx = -1
        best_screen_dist = float('inf')
        best_z = float('inf')

        for atom_idx, sx, sy, sz, radius, color in projected:
            dx = mx - sx
            dy = my - sy
            screen_dist = math.sqrt(dx * dx + dy * dy)

            # Hit radius: use at least 6 px so tiny hydrogens are still clickable,
            # but cap at 1.5× the visual radius for accuracy near dense regions.
            hit_radius = max(6.0, radius * 1.5)

            if screen_dist <= hit_radius:
                # Prefer the atom whose centre is closest to the cursor.
                # Break ties (within 2 px) by depth (front-most wins).
                if screen_dist < best_screen_dist - 2.0:
                    best_screen_dist = screen_dist
                    best_z = sz
                    best_idx = atom_idx
                elif screen_dist < best_screen_dist + 2.0 and sz < best_z:
                    best_z = sz
                    best_idx = atom_idx
                    best_screen_dist = screen_dist

        return best_idx

    # ─── Measurements ────────────────────────────────────────────

    def _complete_distance_measurement(self):
        """Record distance between 2 picked atoms."""
        v = self._v
        if len(v._measure_atoms) < 2:
            return
        i, j = v._measure_atoms[-2], v._measure_atoms[-1]
        a1 = v.molecule.atoms[i]
        a2 = v.molecule.atoms[j]
        if a1.has_coords and a2.has_coords:
            dx = a1.x - a2.x
            dy = a1.y - a2.y
            dz = (a1.z or 0) - (a2.z or 0)
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            v._measurements.append(('dist', i, j, dist))

    def _complete_angle_measurement(self):
        """Record angle between 3 picked atoms (vertex is the middle one)."""
        v = self._v
        if len(v._measure_atoms) < 3:
            return
        i, j, k = v._measure_atoms[-3], v._measure_atoms[-2], v._measure_atoms[-1]
        a1 = v.molecule.atoms[i]
        a2 = v.molecule.atoms[j]
        a3 = v.molecule.atoms[k]
        if a1.has_coords and a2.has_coords and a3.has_coords:
            import numpy as np
            v1 = np.array([a1.x - a2.x, a1.y - a2.y, (a1.z or 0) - (a2.z or 0)])
            v2 = np.array([a3.x - a2.x, a3.y - a2.y, (a3.z or 0) - (a2.z or 0)])
            cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
            angle = math.degrees(math.acos(max(-1, min(1, cos_a))))
            v._measurements.append(('angle', i, j, k, angle))
        # Reset picks after angle measurement
        v._measure_atoms.clear()
