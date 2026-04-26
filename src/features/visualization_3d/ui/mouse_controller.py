"""
Mouse interaction controller for MolViewer3D.

Handles rotation, pan, zoom, atom picking, rubber-band selection,
measurement picks, and auto-rotate.
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
        """Handle mouse press: start rotation, pan, or rubber-band selection."""
        v = self._v
        v._last_mouse_pos = event.position()
        v._mouse_button = event.button()
        v._has_dragged = False

        # Shift + left-click begins rubber-band selection (PyMOL-style)
        if (event.button() == Qt.MouseButton.LeftButton
                and event.modifiers() & Qt.KeyboardModifier.ShiftModifier):
            v._is_selecting = True
            v._sel_rect_origin = event.position()
            v._sel_rect_end = event.position()

    def handle_mouse_move(self, event):
        """Handle mouse move: rotate, pan, update rubber-band, or hover."""
        v = self._v
        if v._last_mouse_pos is None:
            self._detect_hover(event.position())
            return

        # Rubber-band selection drag
        if v._is_selecting:
            v._sel_rect_end = event.position()
            v.update()  # repaint to show the rectangle
            return

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
        """Handle mouse release: commit selection rect, clicks, measurements."""
        v = self._v

        # --- Finish rubber-band selection ---
        if v._is_selecting and event.button() == Qt.MouseButton.LeftButton:
            self._commit_rubber_band_selection()
            v._is_selecting = False
            v._sel_rect_origin = None
            v._sel_rect_end = None
            v._last_mouse_pos = None
            v._mouse_button = None
            v.update()
            return

        was_click = False
        if v._last_mouse_pos is not None:
            moved = (abs(event.position().x() - v._last_mouse_pos.x()) +
                     abs(event.position().y() - v._last_mouse_pos.y()))
            # Consider it a click if the mouse didn't drag at all
            was_click = not getattr(v, '_has_dragged', False) and moved < 3

        btn = event.button()

        if was_click and btn == Qt.MouseButton.LeftButton:
            atom_idx = self._hit_test(event.position())
            if atom_idx >= 0:
                v.atom_clicked.emit(atom_idx)
            else:
                # Click on empty space -> clear selection
                if v.selected_atoms:
                    v.selected_atoms.clear()
                    v.selection_changed.emit(set())
                    v.update()

        if was_click and btn == Qt.MouseButton.RightButton:
            atom_idx = self._hit_test(event.position())
            if atom_idx >= 0:
                # Toggle atom selection
                if atom_idx in v.selected_atoms:
                    v.selected_atoms.discard(atom_idx)
                else:
                    v.selected_atoms.add(atom_idx)
                v.selection_changed.emit(set(v.selected_atoms))
                # Add to measurement picks
                v._measure_atoms.append(atom_idx)
                if len(v._measure_atoms) == 2:
                    self._complete_distance_measurement()
                elif len(v._measure_atoms) == 3:
                    self._complete_angle_measurement()
                v.update()
            else:
                # Right-click on empty: clear measurements
                v._measure_atoms.clear()
                v._measurements.clear()
                v.update()

        v._last_mouse_pos = None
        v._mouse_button = None
        if getattr(v, '_is_interacting', False):
            v._is_interacting = False
            v.update()

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

    # ─── Hit-testing / Hover ─────────────────────────────────────

    def _detect_hover(self, pos):
        v = self._v
        atom_idx = self._hit_test(pos)
        if atom_idx != v._hovered_atom:
            v._hovered_atom = atom_idx
            v.atom_hovered.emit(atom_idx)
            v.update()

    def _hit_test(self, pos):
        v = self._v
        projected = v._renderer._project_atoms(v)
        if not projected:
            return -1

        mx = pos.x()
        my = pos.y()
        best_idx = -1
        best_z = float('inf')

        for atom_idx, sx, sy, sz, radius, color in projected:
            dx = mx - sx
            dy = my - sy
            dist_sq = dx * dx + dy * dy
            if dist_sq <= (radius + 3) ** 2:
                if sz < best_z:
                    best_z = sz
                    best_idx = atom_idx

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
