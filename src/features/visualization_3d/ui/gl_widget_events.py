import math
import numpy as np
from src.shared.qt_compat import Qt, QVector3D, QVector4D, QColor, QApplication, QWheelEvent, QPointF

class GLEventMixin:
    def _auto_fit(self):
        """Fit the camera so all atoms are visible."""
        if self._positions is None or len(self._positions) == 0:
            return

        coords = self._positions
        # Filter out zero-position atoms (atoms without coordinates)
        mask = np.any(coords != 0, axis=1)
        if not np.any(mask):
            return
        valid = coords[mask]

        span = np.max(valid, axis=0) - np.min(valid, axis=0)
        max_span = max(float(np.max(span)), 1.0)

        viewport_size = min(self.width(), self.height())
        if viewport_size < 10:
            viewport_size = 400  # Sensible default before first show
        self.zoom = min(100.0, max(10.0, viewport_size * 0.3 / max_span))
        self.pan_x = 0.0
        self.pan_y = 0.0

    def reset_view(self):
        """Reset camera to defaults and re-fit."""
        self.rot_x = 20.0
        self.rot_y = -30.0
        self.rot_z = 0.0
        self.pan_x = 0.0
        self.pan_y = 0.0
        if self.molecule:
            self._auto_fit()
        self.update()

    def focus_on_atoms(self, atom_indices, padding_angstroms=0.0):
        """Center the view on the given atoms and zoom in."""
        if self._positions is None or len(self._positions) == 0 or not atom_indices:
            return

        old_to_new = getattr(self.molecule, 'properties', {}).get('_old_to_new_idx', {})
        
        coords = []
        for idx in atom_indices:
            new_idx = old_to_new.get(idx, idx)
            if new_idx < len(self._positions):
                coords.append(self._positions[new_idx])

        if not coords:
            return

        coords = np.array(coords)
        centroid = np.mean(coords, axis=0)
        span = np.max(coords, axis=0) - np.min(coords, axis=0)
        max_span = max(float(np.max(span)), 1.0)
        
        max_span += 2.0 * padding_angstroms

        viewport_size = min(self.width(), self.height())
        if viewport_size < 10:
            viewport_size = 400
            
        self.zoom = min(100.0, max(15.0, viewport_size * 0.4 / max_span))

        cos_x = math.cos(math.radians(self.rot_x))
        sin_x = math.sin(math.radians(self.rot_x))
        cos_y = math.cos(math.radians(self.rot_y))
        sin_y = math.sin(math.radians(self.rot_y))
        cos_z = math.cos(math.radians(self.rot_z))
        sin_z = math.sin(math.radians(self.rot_z))

        x, y, z = centroid[0], centroid[1], centroid[2]
        x1 = x * cos_y + z * sin_y
        z1 = -x * sin_y + z * cos_y
        y1 = y * cos_x - z1 * sin_x
        
        x2 = x1 * cos_z - y1 * sin_z
        y2 = x1 * sin_z + y1 * cos_z

        self.pan_x = -x2 * self.zoom
        self.pan_y = y2 * self.zoom

        self.update()

    def mousePressEvent(self, event):
        self._last_mouse_pos = event.position()
        self._mouse_button = event.button()
        self._mouse_moved = False

        modifiers = event.modifiers()

        # Ctrl + left-click begins lasso selection
        if (event.button() == Qt.MouseButton.LeftButton
                and modifiers & Qt.KeyboardModifier.ControlModifier):
            self._is_lasso = True
            self._lasso_path = [QPointF(event.position())]
            return

        # Shift + left-click begins rubber-band selection
        if (event.button() == Qt.MouseButton.LeftButton
                and modifiers & Qt.KeyboardModifier.ShiftModifier):
            self._is_lasso = False
            self._is_selecting = True
            self._sel_rect_origin = event.position()
            self._sel_rect_end = event.position()
            return

        self._is_lasso = False

    def mouseMoveEvent(self, event):
        if self._last_mouse_pos is None:
            return

        dx = event.position().x() - self._last_mouse_pos.x()
        dy = event.position().y() - self._last_mouse_pos.y()
        
        if abs(dx) > 2 or abs(dy) > 2:
            self._mouse_moved = True

        # ── Lasso drag ────────────────────────────────────────────
        if getattr(self, '_is_lasso', False):
            self._lasso_path.append(QPointF(event.position()))
            self.update()
            return

        if self._mouse_button == Qt.MouseButton.LeftButton:
            self.rot_y += dx * 0.5
            self.rot_x += dy * 0.5
        elif self._mouse_button == Qt.MouseButton.RightButton:
            cx = self.width() / 2.0
            cy = self.height() / 2.0
            ang1 = math.atan2(self._last_mouse_pos.y() - cy, self._last_mouse_pos.x() - cx)
            ang2 = math.atan2(event.position().y() - cy, event.position().x() - cx)
            angle_diff = ang2 - ang1
            if angle_diff > math.pi: angle_diff -= 2 * math.pi
            elif angle_diff < -math.pi: angle_diff += 2 * math.pi
            self.rot_z -= math.degrees(angle_diff)
        elif self._mouse_button == Qt.MouseButton.MiddleButton:
            self.pan_x += dx
            self.pan_y += dy

        self._last_mouse_pos = event.position()
        self.update()

    def mouseReleaseEvent(self, event):
        # ── Finish lasso ──────────────────────────────────────────
        if getattr(self, '_is_lasso', False) and event.button() == Qt.MouseButton.LeftButton:
            self._commit_lasso_selection_gl()
            self._is_lasso = False
            self._lasso_path = []
            self._last_mouse_pos = None
            self._mouse_button = None
            self.update()
            return

        if not getattr(self, '_mouse_moved', False) and self._mouse_button == Qt.MouseButton.LeftButton:
            self._handle_click(event.position())
        self._last_mouse_pos = None
        self._mouse_button = None
        self._mouse_moved = False

    def mouseDoubleClickEvent(self, event):
        """Double-click: select entire ligand/fragment that the clicked atom belongs to."""
        if event.button() != Qt.MouseButton.LeftButton:
            return

        atom_idx = self._hit_test(event.position())
        if atom_idx == -1:
            return

        ligand_atoms = self._get_ligand_atoms_gl(atom_idx)
        if ligand_atoms:
            self.selected_atoms = set(ligand_atoms)
            if not hasattr(self, '_selected_atoms_ordered'):
                self._selected_atoms_ordered = []
            self._selected_atoms_ordered = list(ligand_atoms)[:3]
            self.selection_changed.emit(set(self.selected_atoms))
            self.update()

    def _get_ligand_atoms_gl(self, atom_idx):
        """
        Return the set of atom indices for the ligand/fragment containing atom_idx.

        For PDB proteins: groups HETATM atoms by (res_name, res_seq, chain_id).
        For SDF/MOL2/MOL: BFS to find the whole connected fragment.
        """
        if not self.molecule or atom_idx >= len(self.molecule.atoms):
            return set()

        atom = self.molecule.atoms[atom_idx]
        is_protein = getattr(self.molecule, 'properties', {}).get('is_protein', False)

        if is_protein and getattr(atom, 'is_hetatm', False):
            target_res = (
                getattr(atom, 'res_name', None),
                getattr(atom, 'res_seq', None),
                getattr(atom, 'chain_id', None),
            )
            if target_res[0] and target_res[0].upper() in ('HOH', 'WAT', 'SOL', 'DOD'):
                return {atom_idx}

            ligand_indices = set()
            for a in self.molecule.atoms:
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
            # BFS over molecule adjacency
            adjacency = self.molecule._adjacency
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

    def _commit_lasso_selection_gl(self):
        """
        Finalise Ctrl+drag lasso for the GL viewer: project all atoms, test each
        projected point against the freehand polygon, add enclosed atoms to selected_atoms.
        """
        lasso = getattr(self, '_lasso_path', [])
        if len(lasso) < 3:
            self.selected_atoms.clear()
            self.selection_changed.emit(set())
            return

        from src.shared.qt_compat import QPainterPath, QPointF as _QPointF

        path = QPainterPath()
        path.moveTo(lasso[0])
        for pt in lasso[1:]:
            path.lineTo(pt)
        path.closeSubpath()

        if self._positions is None or len(self._positions) == 0:
            return

        w = self.width()
        h = self.height()
        proj = self._get_projection_matrix()
        view = self._get_view_matrix()
        mvp = proj * view

        old_to_new = getattr(self.molecule, 'properties', {}).get('_old_to_new_idx', {})
        new_to_old = {v: k for k, v in old_to_new.items()}

        newly_selected = set()
        for i in range(len(self._positions)):
            vec = QVector3D(
                float(self._positions[i][0]),
                float(self._positions[i][1]),
                float(self._positions[i][2]),
            )
            clip = mvp.map(vec)
            view_pos = view.map(vec)
            if view_pos.z() > 0:
                continue

            sx = (clip.x() + 1.0) * 0.5 * w
            sy = (1.0 - clip.y()) * 0.5 * h

            if path.contains(_QPointF(sx, sy)):
                newly_selected.add(new_to_old.get(i, i))

        if newly_selected:
            self.selected_atoms |= newly_selected
        self.selection_changed.emit(set(self.selected_atoms))

    def _handle_click(self, pos):
        from src.shared.qt_compat import QApplication
        modifiers = QApplication.keyboardModifiers()
        
        atom_idx = self._hit_test(pos)
        if atom_idx != -1:
            if modifiers & Qt.KeyboardModifier.ShiftModifier:
                if atom_idx in self.selected_atoms:
                    self.selected_atoms.remove(atom_idx)
                    if hasattr(self, '_selected_atoms_ordered') and atom_idx in self._selected_atoms_ordered:
                        self._selected_atoms_ordered.remove(atom_idx)
                else:
                    self.selected_atoms.add(atom_idx)
                    if not hasattr(self, '_selected_atoms_ordered'):
                        self._selected_atoms_ordered = []
                    self._selected_atoms_ordered.append(atom_idx)
            else:
                self.selected_atoms = {atom_idx}
                if not hasattr(self, '_selected_atoms_ordered'):
                    self._selected_atoms_ordered = []
                self._selected_atoms_ordered.append(atom_idx)
                if len(self._selected_atoms_ordered) > 3:
                    self._selected_atoms_ordered.pop(0)
            self.selection_changed.emit(set(self.selected_atoms))
            self.atom_clicked.emit(atom_idx)
        else:
            self.selected_atoms.clear()
            self._selected_atoms_ordered = []
            self.selection_changed.emit(set())
            
        self.update()

    def _hit_test(self, pos):
        """
        Return the atom index closest to *pos* in screen-space.

        Improved accuracy over the original:
        - Candidates must lie within max(6px, 1.5× projected radius) of the cursor.
        - The candidate with the **smallest screen-space distance to centre** wins.
        - Depth (view-space z) is used only to break ties within 2 px of each other.
        - This prevents large cartoon-mesh atoms from hijacking clicks that visually
          land on a foreground ligand atom.
        """
        if self._positions is None or len(self._positions) == 0:
            return -1

        w = self.width()
        h = self.height()
        proj = self._get_projection_matrix()
        view = self._get_view_matrix()
        mvp = proj * view
        
        click_x = pos.x()
        click_y = pos.y()
        
        best_idx = -1
        best_screen_dist = float('inf')
        best_z = float('inf')
        
        old_to_new = getattr(self.molecule, 'properties', {}).get('_old_to_new_idx', {})
        new_to_old = {v: k for k, v in old_to_new.items()}
        
        for i in range(len(self._positions)):
            vec = QVector3D(
                float(self._positions[i][0]),
                float(self._positions[i][1]),
                float(self._positions[i][2]),
            )
            clip = mvp.map(vec)
            
            # Behind camera — skip
            view_pos = view.map(vec)
            if view_pos.z() > 0:
                continue
                
            sx = (clip.x() + 1.0) * 0.5 * w
            sy = (1.0 - clip.y()) * 0.5 * h
            
            dx = sx - click_x
            dy = sy - click_y
            screen_dist = math.sqrt(dx * dx + dy * dy)
            
            # Projected sphere radius at this depth
            if view_pos.z() < -0.1:
                scale = (h / 2.0) / (math.tan(math.radians(22.5)) * -view_pos.z())
            else:
                scale = 1.0
            radius = self._radii[i] * scale * self.sphere_scale

            # Hit tolerance: at least 6 px (to keep tiny atoms clickable),
            # but only up to 1.5× the actual sphere radius for precision.
            hit_radius = max(6.0, radius * 1.5)
            
            if screen_dist <= hit_radius:
                vz = view_pos.z()
                # Prefer closest to cursor centre; break ties (within 2 px) by depth
                if screen_dist < best_screen_dist - 2.0:
                    best_screen_dist = screen_dist
                    best_z = vz
                    best_idx = new_to_old.get(i, i)
                elif screen_dist < best_screen_dist + 2.0 and vz < best_z:
                    best_z = vz
                    best_idx = new_to_old.get(i, i)
                    best_screen_dist = screen_dist
                    
        return best_idx

    def wheelEvent(self, event: QWheelEvent):
        delta = event.angleDelta().y()
        factor = 1.1 if delta > 0 else 0.9
        self.zoom *= factor
        self.zoom = max(5.0, min(200.0, self.zoom))
        self.update()

    def contextMenuEvent(self, event):
        """Show context menu for selected atoms/residues and general viewer options."""
        from src.shared.qt_compat import QMenu
        menu = QMenu(self)

        selected_res_seqs = set()
        bs_action, wf_action, sf_action = None, None, None
        show_sc_action, hide_sc_action = None, None
        label_res_action, clear_label_action = None, None

        if self.selected_atoms:
            # Styles
            style_menu = menu.addMenu("Set Style")
            bs_action = style_menu.addAction("Ball and Stick")
            wf_action = style_menu.addAction("Wireframe")
            sf_action = style_menu.addAction("Space Fill")

            # Determine if selected contains residues
            if self.molecule:
                for idx in self.selected_atoms:
                    atom = self.molecule.atoms[idx]
                    rs = getattr(atom, 'res_seq', None)
                    if rs is not None:
                        selected_res_seqs.add(rs)

            if selected_res_seqs:
                sidechain_menu = menu.addMenu("Side Chains")
                show_sc_action = sidechain_menu.addAction("Show")
                hide_sc_action = sidechain_menu.addAction("Hide")

                # Label action
                label_res_action = menu.addAction("Label Residue Color...")
                clear_label_action = menu.addAction("Clear Residue Label")

            menu.addSeparator()

        action = menu.exec(event.globalPos())
        if not action:
            return

        if not hasattr(self, 'custom_atom_modes'):
            self.custom_atom_modes = {}
        if not hasattr(self, 'sidechain_res_vis'):
            self.sidechain_res_vis = {}

        if bs_action and action == bs_action:
            for idx in self.selected_atoms:
                self.custom_atom_modes[idx] = 'ball_and_stick'
            self.update()
        elif wf_action and action == wf_action:
            for idx in self.selected_atoms:
                self.custom_atom_modes[idx] = 'wireframe'
            self.update()
        elif sf_action and action == sf_action:
            for idx in self.selected_atoms:
                self.custom_atom_modes[idx] = 'spacefill'
            self.update()
        elif show_sc_action and action == show_sc_action:
            for rs in selected_res_seqs:
                self.sidechain_res_vis[rs] = True
            if self.gl_available: self._update_gl_buffers()
            self.update()
        elif hide_sc_action and action == hide_sc_action:
            for rs in selected_res_seqs:
                self.sidechain_res_vis[rs] = False
            if self.gl_available: self._update_gl_buffers()
            self.update()
        elif label_res_action and action == label_res_action:
            from PySide6.QtWidgets import QColorDialog
            color = QColorDialog.getColor(Qt.white, self, "Select Residue Label Color")
            if color.isValid():
                for rs in selected_res_seqs:
                    self.labeled_residues[rs] = color
                self.update()
        elif clear_label_action and action == clear_label_action:
            for rs in selected_res_seqs:
                if rs in self.labeled_residues:
                    del self.labeled_residues[rs]
            self.update()

    @property
    def show_sidechains(self): return getattr(self, '_show_sidechains', False)

    @show_sidechains.setter
    def show_sidechains(self, val):
        if getattr(self, '_show_sidechains', None) == val:
            return
        self._show_sidechains = val
        if self.gl_available:
            self._update_gl_buffers()
        self.update()

    def clear(self):
        """Remove the current molecule and clear the view."""
        self.set_molecule(None)
