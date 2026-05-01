"""
2D Molecular Viewer — ChemDraw-quality skeletal formula renderer.

Renders molecules in publication-quality 2D skeletal notation using QPainter.
Uses CoordinateGenerator2D for layout.

ChemDraw-like aesthetics:
- Clean black bonds on white/dark background
- Proper bond width proportions relative to bond length
- Double bond inner line with proper spacing
- Triple bond with even triple lines
- Wedge and hash bonds for stereochemistry
- Heteroatom labels in element-specific colors
- **Direction-aware H labels** (NH vs HN, OH vs HO) — prevents overlap
- Implicit H notation with subscript counts
- Background clearing behind labels
- Proper skeletal formula (C vertices hidden)
- Hydrogen toggle support
- Selection highlighting
- Distance / angle measurement overlay
- Atom tooltip on hover
"""

import math
from src.shared.qt_compat import QWidget, Qt, QPointF, QRectF, Signal
from src.shared.qt_compat import (
    QPainter, QColor, QPen, QBrush, QFont, QFontMetrics,
    QImage, QPainterPath, QPolygonF, QMenu, QAction
)
from src.shared.ui.theme import COLORS

from src.features.visualization_2d.ui.bond_renderer_2d import BondRenderer2D
from src.features.visualization_2d.ui.atom_renderer_2d import (
    AtomRenderer2D,
    ELEMENT_COLORS,
    ELEMENT_COLORS_LIGHT,
)

# Set to True during development to enable [DEBUG] prints
_DEBUG = False


class MolViewer2D(QWidget):
    """
    Publication-quality 2D skeletal formula viewer.
    Mimics the aesthetic of ChemDraw / RDKit depiction.

    Selection
    ---------
    **Shift + left-drag** draws a rubber-band rectangle (PyMOL-style).
    Atoms whose 2D screen positions fall inside the rectangle on
    mouse-release are added to ``selected_atoms``.  A plain left-click
    on empty space clears the selection.

    Atom Dragging
    -------------
    **Left-click + drag** on a selected atom to reposition it.
    The atom follows the mouse cursor, and connected bonds stretch
    to maintain connectivity. This allows interactive adjustment
    of the 2D layout without breaking molecular structure.

    Deletion
    --------
    Pressing the **Delete** key while atoms are selected emits
    ``delete_requested`` so the main window can remove those atoms
    from the domain model and refresh both viewers.

    Context Menu (Right-Click)
    --------------------------
    Right-clicking in the 2D view opens a context menu with:
    
    - **Fit to Screen** — Auto-fit molecule to view ('F' key)
    - **Reset View** — Reset zoom and pan to default ('R' key)

    Keyboard Shortcuts
    ------------------
    - **F**      — Fit molecule to screen
    - **R**      — Reset view (fit + center)
    - **Delete** — Delete selected atoms
    """

    # --- Signals ---
    selection_changed = Signal(object)   # emits set of selected atom indices
    delete_requested = Signal(object)    # emits set of atom indices to delete

    def __init__(self, parent=None):
        super().__init__(parent)
        self.molecule = None
        self.coords_2d = {}
        self.selected_atoms = set()
        self.setMinimumSize(400, 400)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)  # receive key events
        self.setMouseTracking(True)  # For hover tooltips

        # Display settings
        self.bg_color = QColor(COLORS['viewer_bg'])
        self._scale = 40.0
        self._offset_x = 0.0
        self._offset_y = 0.0
        self.show_hydrogens = False  # ChemDraw default: hide H
        self.show_all_labels = False
        self.show_protein_placeholder = False  # For large proteins

        # Rendering quality tuning (relative to bond length)
        self._bond_width_ratio = 0.028     # Bond line width / bond pixel length
        self._double_offset_ratio = 0.15   # Offset for inner double bond line
        self._triple_offset_ratio = 0.18   # Offset for triple bond lines
        self._label_padding = 4            # px padding around labels
        self._wedge_width_ratio = 0.12     # Half-width of wedge tip

        # Mouse
        self._last_mouse_pos = None
        self._mouse_button = None
        self._has_dragged_right = False
        self._is_dark_bg = True
        self._hovered_atom = -1            # For tooltip on hover

        # Rubber-band selection rectangle (screen coords)
        self._sel_rect_origin = None   # QPointF or None
        self._sel_rect_end = None      # QPointF or None
        self._is_selecting = False     # True while Shift+left-drag is active

        # Atom drag state for interactive repositioning
        self._is_dragging_atom = False  # True while dragging a selected atom
        self._dragged_atom_idx = None   # Index of atom being dragged
        self._drag_start_pos = None     # Starting screen position of drag

        # Explicit H atom drag state for heteroatom H repositioning
        self._is_dragging_h = False     # True while dragging an explicit H
        self._dragged_h_key = None      # (parent_idx, h_idx) tuple for dragged H
        self._dragged_h_parent_pos = None  # Screen position of parent atom

        # Export state (set during export_image calls)
        self._original_scale = None    # Pre-export scale for font sizing

        # Sub-renderers
        self._atom_renderer = AtomRenderer2D(self)
        self._bond_renderer = BondRenderer2D(self)

    def set_molecule(self, molecule, coords_2d=None):
        """
        Set the molecule and its 2D coordinates.
        If coordinates are not provided, they are generated using the
        SMILES-based OASA approach (mol -> SMILES -> oasa -> 2D coordinates).
        
        This provides a unified, high-quality 2D layout for all molecule sources
        (SMILES, MOL, SDF, MOL2) by leveraging OASA's robust SMILES parser
        which generates professional 2D coordinates during parsing.
        """
        self.molecule = molecule
        self.selected_atoms = set()
        self.coords_2d = {}

        if coords_2d:
            if _DEBUG:
                print(f"[DEBUG] Using explicitly provided 2D coordinates: {len(coords_2d)}")
            self.coords_2d = coords_2d
        elif molecule and molecule.atoms:
            # Check if coordinates already exist on atoms (e.g. from parser x2d/y2d or sketcher x/y)
            has_x2d = all(hasattr(a, 'x2d') and a.x2d is not None for a in molecule.atoms if a.symbol != 'H')
            has_xy = all(hasattr(a, 'x') and a.x is not None for a in molecule.atoms if a.symbol != 'H')

            if has_x2d:
                # Use the native OASA coordinates directly — these are the
                # high-quality coords generated by text_to_mol(calc_coords=1)
                # during parse_smiles.  Avoids the lossy domain→OASA
                # round-trip in CoordinateGenerator2DSMILES.
                if _DEBUG: print("[DEBUG] Using native OASA x2d/y2d coordinates.")
                self.coords_2d = {a.index: [a.x2d, a.y2d] for a in molecule.atoms
                                  if hasattr(a, 'x2d') and a.x2d is not None}
            elif has_xy and not has_x2d:
                # Use manual coordinates from sketcher
                if _DEBUG: print("[DEBUG] Using manual coordinates (x, y) from sketcher.")
                self.coords_2d = {a.index: [a.x, -a.y] for a in molecule.atoms}
            
            if not self.coords_2d:
                # Fallback: Use SMILES-based OASA coordinate generation
                from src.features.layout_2d.generators.coordgen2d_smiles_pure_oasa import CoordinateGenerator2DSMILES
                if _DEBUG:
                    print("[DEBUG] Fallback: Using CoordinateGenerator2DSMILES...")
                generator = CoordinateGenerator2DSMILES(molecule, force_regenerate=True)
                self.coords_2d = generator.generate()

            # Ensure ring perception is up-to-date for correct double bond placement
            molecule.find_rings()

            # Professionally Kekulize aromatic bonds using robust graph matching
            from src.features.smiles_parser.rules.aromaticity import kekulize
            kekulize(molecule)

        else:
            self.coords_2d = {}

        # Clear custom H angles when loading new molecule
        self._atom_renderer.clear_explicit_h_angles()

        self._auto_fit()
        self.update()

    def clear(self):
        self.molecule = None
        self.coords_2d = {}
        self.selected_atoms = set()
        self._atom_renderer.clear_explicit_h_angles()
        self.update()

    def set_selected(self, atom_indices):
        self.selected_atoms = set(atom_indices)
        self.update()

    # ─── Rendering ────────────────────────────────────────────────

    def paintEvent(self, event):
        painter = QPainter(self)
        self._render(painter, self.width(), self.height())
        painter.end()

    def _render(self, painter, width, height, is_export=False, export_scale=1.0):
        """Render the 2D view into *painter* at (width, height).

        When called from an export/print path, pass ``is_export=True``.
        ``export_scale`` (default 1.0) is a DPR multiplier used by
        ``file_operations.print_views`` and future high-DPI image save
        callers — when > 1.0 the widget zoom/pan are temporarily scaled
        up for the duration of the render so the molecule fills the
        larger canvas, and ``_original_scale`` is set so the font system
        (``AtomRenderer2D._get_font``) stabilises text size instead of
        bloating it linearly with zoom. State is restored in ``finally``.
        """
        # Save state we may temporarily mutate for export scaling
        _saved_scale = self._scale
        _saved_ox = self._offset_x
        _saved_oy = self._offset_y
        _saved_orig = self._original_scale
        _mutated = False
        try:
            if is_export and export_scale and export_scale != 1.0:
                self._original_scale = _saved_scale
                self._scale = _saved_scale * export_scale
                self._offset_x = _saved_ox * export_scale
                self._offset_y = _saved_oy * export_scale
                _mutated = True

            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
            painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)

            self._is_dark_bg = (self.bg_color.lightnessF() < 0.4)
            painter.fillRect(0, 0, width, height, self.bg_color)

            # Show placeholder for large proteins instead of trying to generate 2D
            if self.show_protein_placeholder and self.molecule:
                self._draw_protein_placeholder(painter, width, height)
                return

            if not self.molecule or not self.coords_2d:
                if not is_export:
                    self._draw_placeholder(painter, width, height)
                return

            visible = self._get_visible_atoms()
            has_label = self._compute_labels(visible)

            # Layer 1: Bonds
            self._bond_renderer.draw_all_bonds(painter, visible, has_label)

            # Layer 2: Selection (draw before labels so circles are visible)
            if self.selected_atoms:
                self._draw_selection(painter, visible)

            # Layer 3: Atom labels
            self._atom_renderer.draw_all_labels(painter, visible, has_label)

            # Layer 4: Distance / angle measurement overlay
            if len(self.selected_atoms) in (2, 3):
                self._draw_measurement_overlay(painter, visible)

            # Layer 5: Rubber-band selection rectangle (while Shift+dragging)
            if self._is_selecting and self._sel_rect_origin and self._sel_rect_end:
                self._draw_rubber_band(painter)

            # Layer 6: Hover tooltip
            if self._hovered_atom >= 0 and not is_export:
                self._draw_hover_tooltip(painter)

            if not is_export:
                self._draw_overlay(painter)
        finally:
            if _mutated:
                self._scale = _saved_scale
                self._offset_x = _saved_ox
                self._offset_y = _saved_oy
                self._original_scale = _saved_orig

    def _get_visible_atoms(self):
        """Return set of atom indices that should be rendered."""
        if not self.molecule:
            return set()
        visible = set()
        for atom in self.molecule.atoms:
            if atom.index not in self.coords_2d:
                continue
            if atom.symbol == 'H' and not self.show_hydrogens:
                # Always hide hydrogen atoms - they are shown as part of parent atom label (CH3, OH, NH, etc.)
                continue
            visible.add(atom.index)
        return visible

    def _compute_labels(self, visible):
        """Decide which atoms get text labels."""
        has_label = {}
        for idx in visible:
            atom = self.molecule.atoms[idx]
            if self.show_all_labels:
                has_label[idx] = True
            elif atom.symbol == 'H':
                has_label[idx] = True  # Always label visible H
            elif atom.symbol == 'C' and atom.formal_charge == 0:
                # Skeletal: hide C labels unless terminal or charged
                degree = sum(1 for n in self.molecule.get_neighbors(idx) if n in visible)
                has_label[idx] = (degree <= 1)
            else:
                has_label[idx] = True  # All heteroatoms get labels
        return has_label

    def _to_screen(self, mx, my):
        sx = mx * self._scale + self._offset_x
        sy = -my * self._scale + self._offset_y
        return sx, sy

    def _bond_color(self):
        return QColor(220, 225, 233) if self._is_dark_bg else QColor(30, 30, 35)

    def _bond_width(self):
        return max(1.5, self._scale * self._bond_width_ratio)

    # ─── Selection ────────────────────────────────────────────────

    def _draw_selection(self, painter, visible):
        for idx in self.selected_atoms:
            if idx not in visible or idx not in self.coords_2d:
                continue
            mx, my = self.coords_2d[idx]
            sx, sy = self._to_screen(mx, my)
            r = self._scale * 0.40  # Increased from 0.28 to be visible beyond label backgrounds

            # Outer glow
            pen = QPen(QColor(255, 200, 50, 140), max(2.5, self._scale * 0.045))
            painter.setPen(pen)
            painter.setBrush(QBrush(QColor(255, 200, 50, 35)))
            painter.drawEllipse(QPointF(sx, sy), r, r)

    # ─── Distance / Angle Measurement ─────────────────────────────

    def _draw_measurement_overlay(self, painter, visible):
        """Draw distance or angle measurement when 2 or 3 atoms are selected."""
        selected = sorted(self.selected_atoms)
        if len(selected) == 2:
            self._draw_distance(painter, selected[0], selected[1])
        elif len(selected) == 3:
            self._draw_angle(painter, selected[0], selected[1], selected[2])

    def _draw_distance(self, painter, idx1, idx2):
        """Draw a dashed line between two selected atoms with distance label."""
        if idx1 not in self.coords_2d or idx2 not in self.coords_2d:
            return

        mx1, my1 = self.coords_2d[idx1]
        mx2, my2 = self.coords_2d[idx2]
        sx1, sy1 = self._to_screen(mx1, my1)
        sx2, sy2 = self._to_screen(mx2, my2)

        # Dashed measurement line
        pen = QPen(QColor(255, 160, 0, 180), 1.5)
        pen.setStyle(Qt.PenStyle.DashDotLine)
        painter.setPen(pen)
        painter.drawLine(QPointF(sx1, sy1), QPointF(sx2, sy2))

        # Distance label
        dist_2d = math.hypot(mx2 - mx1, my2 - my1)
        # Try 3D distance if available
        a1 = self.molecule.atoms[idx1]
        a2 = self.molecule.atoms[idx2]
        label = ""
        if a1.has_coords and a2.has_coords:
            dist_3d = math.sqrt((a2.x - a1.x)**2 + (a2.y - a1.y)**2 + (a2.z - a1.z)**2)
            label = f"{dist_3d:.2f} Å"
        else:
            label = f"{dist_2d:.2f}"

        mid_x = (sx1 + sx2) / 2
        mid_y = (sy1 + sy2) / 2

        font = QFont('Arial', 10)
        font.setBold(True)
        painter.setFont(font)
        fm = QFontMetrics(font)
        tw = fm.horizontalAdvance(label)

        # Background pill
        bg = QRectF(mid_x - tw / 2 - 4, mid_y - 9, tw + 8, 18)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QBrush(QColor(40, 40, 60, 200)))
        painter.drawRoundedRect(bg, 4, 4)

        # Text
        painter.setPen(QColor(255, 200, 80))
        painter.drawText(QPointF(mid_x - tw / 2, mid_y + 5), label)

    def _draw_angle(self, painter, idx1, idx2, idx3):
        """Draw angle measurement arc and label for three selected atoms."""
        for idx in (idx1, idx2, idx3):
            if idx not in self.coords_2d:
                return

        # Vertex is the middle atom
        ax, ay = self.coords_2d[idx1]
        bx, by = self.coords_2d[idx2]
        cx, cy = self.coords_2d[idx3]

        # Try 3D angle
        a1 = self.molecule.atoms[idx1]
        a2 = self.molecule.atoms[idx2]
        a3 = self.molecule.atoms[idx3]
        if a1.has_coords and a2.has_coords and a3.has_coords:
            v1 = (a1.x - a2.x, a1.y - a2.y, a1.z - a2.z)
            v2 = (a3.x - a2.x, a3.y - a2.y, a3.z - a2.z)
            dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
            m1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
            m2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
            if m1 > 0 and m2 > 0:
                cos_a = max(-1, min(1, dot / (m1 * m2)))
                angle_deg = math.degrees(math.acos(cos_a))
            else:
                angle_deg = 0
        else:
            # 2D angle
            v1 = (ax - bx, ay - by)
            v2 = (cx - bx, cy - by)
            dot = v1[0]*v2[0] + v1[1]*v2[1]
            m1 = math.hypot(*v1)
            m2 = math.hypot(*v2)
            if m1 > 0 and m2 > 0:
                cos_a = max(-1, min(1, dot / (m1 * m2)))
                angle_deg = math.degrees(math.acos(cos_a))
            else:
                angle_deg = 0

        # Draw label at vertex
        sx, sy = self._to_screen(bx, by)
        label = f"{angle_deg:.1f}°"

        font = QFont('Arial', 10)
        font.setBold(True)
        painter.setFont(font)
        fm = QFontMetrics(font)
        tw = fm.horizontalAdvance(label)

        lx, ly = sx + 15, sy - 15
        bg = QRectF(lx - 4, ly - 12, tw + 8, 18)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QBrush(QColor(40, 40, 60, 200)))
        painter.drawRoundedRect(bg, 4, 4)

        painter.setPen(QColor(100, 255, 180))
        painter.drawText(QPointF(lx, ly + 2), label)

    # ─── Hover Tooltip ────────────────────────────────────────────

    def _draw_hover_tooltip(self, painter):
        """Draw a tooltip near the hovered atom showing basic info."""
        if self._hovered_atom < 0 or not self.molecule:
            return
        if self._hovered_atom >= len(self.molecule.atoms):
            return
        if self._hovered_atom not in self.coords_2d:
            return

        atom = self.molecule.atoms[self._hovered_atom]
        mx, my = self.coords_2d[self._hovered_atom]
        sx, sy = self._to_screen(mx, my)

        lines = [
            f"{atom.symbol}{self._hovered_atom + 1}",
            f"q: {atom.partial_charge:.3f}" if atom.partial_charge != 0.0 else "",
        ]
        if atom.has_coords:
            lines.append(f"({atom.x:.2f}, {atom.y:.2f}, {atom.z:.2f})")
        lines = [l for l in lines if l]

        font = QFont('Segoe UI', 9)
        painter.setFont(font)
        fm = QFontMetrics(font)
        line_h = fm.height()
        max_w = max(fm.horizontalAdvance(l) for l in lines) if lines else 40
        box_w = max_w + 12
        box_h = line_h * len(lines) + 8

        # Position tooltip above-right of atom
        tx = sx + 12
        ty = sy - box_h - 8

        # Background
        painter.setPen(QPen(QColor(100, 140, 255, 180), 1))
        painter.setBrush(QBrush(QColor(20, 24, 45, 220)))
        painter.drawRoundedRect(QRectF(tx, ty, box_w, box_h), 5, 5)

        # Text
        painter.setPen(QColor(230, 230, 255))
        y = ty + line_h
        for line in lines:
            painter.drawText(QPointF(tx + 6, y), line)
            y += line_h

    # ─── Overlay ──────────────────────────────────────────────────

    def _draw_overlay(self, painter):
        if not self.molecule:
            return
        font = QFont('Segoe UI', 10)
        painter.setFont(font)
        painter.setPen(QColor(COLORS['text_muted']))

        info = f"2D | {len(self.molecule.atoms)} atoms, {len(self.molecule.bonds)} bonds"
        if not self.show_hydrogens:
            info += " | H hidden"
        if self.selected_atoms:
            info += f" | {len(self.selected_atoms)} selected"
        painter.drawText(8, 18, info)

    def _draw_placeholder(self, painter, width, height):
        font = QFont('Segoe UI', 16)
        painter.setFont(font)
        painter.setPen(QColor(COLORS['text_muted']))
        painter.drawText(QRectF(0, 0, width, height), Qt.AlignmentFlag.AlignCenter,
                         "2D structure will appear here")

    # ─── Image Export ─────────────────────────────────────────────

    def export_image(self, filepath, dpi=300, bg_white=True):
        """Export the 2D view as a high-resolution image.

        Args:
            filepath: Output file path.
            dpi: Resolution in dots per inch.
            bg_white: If True (default), use white background for publication quality.
        """
        scale_factor = dpi / 96.0
        img_w = int(self.width() * scale_factor)
        img_h = int(self.height() * scale_factor)

        image = QImage(img_w, img_h, QImage.Format.Format_ARGB32_Premultiplied)
        image.setDotsPerMeterX(int(dpi / 0.0254))
        image.setDotsPerMeterY(int(dpi / 0.0254))

        orig_bg = self.bg_color
        orig_scale = self._scale
        orig_ox = self._offset_x
        orig_oy = self._offset_y

        # Store original scale so _get_font() uses it instead of the inflated one
        self._original_scale = orig_scale

        # Always default to white background for professional exports
        if bg_white:
            self.bg_color = QColor(255, 255, 255)

        self._scale *= scale_factor
        self._offset_x *= scale_factor
        self._offset_y *= scale_factor

        painter = QPainter(image)
        self._render(painter, img_w, img_h, is_export=True)
        painter.end()

        # Restore
        self._original_scale = None
        self.bg_color = orig_bg
        self._scale = orig_scale
        self._offset_x = orig_ox
        self._offset_y = orig_oy

        return image.save(filepath)

    # ─── Mouse Interaction ────────────────────────────────────────

    def mousePressEvent(self, event):
        """Handle mouse press: start pan, rubber-band selection, or atom drag."""
        self._last_mouse_pos = event.position()
        self._mouse_button = event.button()

        if event.button() == Qt.MouseButton.RightButton:
            self._has_dragged_right = False

        # Shift + left-click begins rubber-band selection (PyMOL-style)
        if (event.button() == Qt.MouseButton.LeftButton
                and event.modifiers() & Qt.KeyboardModifier.ShiftModifier):
            self._is_selecting = True
            self._sel_rect_origin = event.position()
            self._sel_rect_end = event.position()
            return

        # Left-click on a selected atom starts drag mode
        if event.button() == Qt.MouseButton.LeftButton:
            atom_idx = self._hit_test_2d(event.position())
            if atom_idx >= 0 and atom_idx in self.selected_atoms:
                # Start dragging the selected atom
                self._is_dragging_atom = True
                self._dragged_atom_idx = atom_idx
                self._drag_start_pos = event.position()
                if _DEBUG:
                    print(f"[DEBUG 2D] Started dragging atom {atom_idx}")
                return

            # Check for explicit H atom hit (for heteroatoms like N, O, S, P)
            h_key = self._atom_renderer.hit_test_explicit_h(event.position())
            if h_key is not None:
                # Start dragging the explicit H atom
                self._is_dragging_h = True
                self._dragged_h_key = h_key
                self._drag_start_pos = event.position()
                # Get parent position for reference
                h_data = self._atom_renderer.get_explicit_h_data(h_key)
                if h_data:
                    self._dragged_h_parent_pos = QPointF(h_data[2], h_data[3])
                if _DEBUG:
                    print(f"[DEBUG 2D] Started dragging H atom {h_key}")

    def mouseMoveEvent(self, event):
        """Handle mouse move: pan, atom drag, H atom repositioning, update rubber-band, or hover."""
        # Explicit H atom dragging - rotate H around parent atom
        if self._is_dragging_h and self._dragged_h_key is not None:
            # Update H position based on mouse position (rotates around parent)
            self._atom_renderer.update_h_position(
                self._dragged_h_key[0], self._dragged_h_key[1], event.position())
            self._last_mouse_pos = event.position()
            self.update()
            return

        # Atom dragging - update atom position in molecule coordinates
        if self._is_dragging_atom and self._dragged_atom_idx is not None:
            if self._dragged_atom_idx in self.coords_2d:
                # Convert screen delta to molecule coordinates
                dx_screen = event.position().x() - self._last_mouse_pos.x()
                dy_screen = event.position().y() - self._last_mouse_pos.y()

                # Convert screen pixels to molecule units
                dx_mol = dx_screen / self._scale
                dy_mol = -dy_screen / self._scale  # Y is inverted in screen coords

                # Update atom position
                x, y = self.coords_2d[self._dragged_atom_idx]
                self.coords_2d[self._dragged_atom_idx] = [x + dx_mol, y + dy_mol]

                # Update atom's cached coordinates
                atom = self.molecule.get_atom(self._dragged_atom_idx)
                if atom:
                    atom.x2d = x + dx_mol
                    atom.y2d = y + dy_mol

                self._last_mouse_pos = event.position()
                self.update()
            return

        # Hover detection (only when not dragging)
        if self._last_mouse_pos is None:
            old_hover = self._hovered_atom
            self._hovered_atom = self._hit_test_2d(event.position())
            if self._hovered_atom != old_hover:
                self.update()
            return

        # Rubber-band selection drag
        if self._is_selecting:
            self._sel_rect_end = event.position()
            self.update()
            return

        # Normal pan and rotate operation
        dx = event.position().x() - self._last_mouse_pos.x()
        dy = event.position().y() - self._last_mouse_pos.y()

        if self._mouse_button == Qt.MouseButton.LeftButton:
            self._offset_x += dx
            self._offset_y += dy
        elif self._mouse_button == Qt.MouseButton.RightButton:
            if abs(dx) > 0 or abs(dy) > 0:
                self._has_dragged_right = True
                self._rotate_2d_mouse(event.position(), self._last_mouse_pos)

        self._last_mouse_pos = event.position()
        self.update()

    def mouseReleaseEvent(self, event):
        """Handle mouse release: end drag, commit selection, or deselect on click."""
        # --- Finish H atom dragging ---
        if self._is_dragging_h and event.button() == Qt.MouseButton.LeftButton:
            if _DEBUG:
                print(f"[DEBUG 2D] Finished dragging H atom {self._dragged_h_key}")
            self._is_dragging_h = False
            self._dragged_h_key = None
            self._dragged_h_parent_pos = None
            self._drag_start_pos = None
            self._last_mouse_pos = None
            self._mouse_button = None
            self.update()
            return

        # --- Finish atom dragging ---
        if self._is_dragging_atom and event.button() == Qt.MouseButton.LeftButton:
            if _DEBUG:
                print(f"[DEBUG 2D] Finished dragging atom {self._dragged_atom_idx}")
            self._is_dragging_atom = False
            self._dragged_atom_idx = None
            self._drag_start_pos = None
            self._last_mouse_pos = None
            self._mouse_button = None
            self.update()
            return

        # --- Finish rubber-band selection ---
        if self._is_selecting and event.button() == Qt.MouseButton.LeftButton:
            self._commit_rubber_band_selection()
            self._is_selecting = False
            self._sel_rect_origin = None
            self._sel_rect_end = None
            self._last_mouse_pos = None
            self._mouse_button = None
            self.update()
            return

        # Detect simple click (<3px movement)
        was_click = False
        if self._last_mouse_pos is not None:
            moved = (abs(event.position().x() - self._last_mouse_pos.x()) +
                     abs(event.position().y() - self._last_mouse_pos.y()))
            was_click = moved < 3

        if was_click and event.button() == Qt.MouseButton.LeftButton:
            # Hit-test: find nearest atom under cursor
            atom_idx = self._hit_test_2d(event.position())
            if atom_idx >= 0:
                # Toggle single atom selection
                if atom_idx in self.selected_atoms:
                    self.selected_atoms.discard(atom_idx)
                else:
                    self.selected_atoms.add(atom_idx)
                self.selection_changed.emit(set(self.selected_atoms))
                self.update()
            else:
                # Click on empty space → clear selection
                if self.selected_atoms:
                    self.selected_atoms.clear()
                    self.selection_changed.emit(set())
                    self.update()

        self._last_mouse_pos = None
        self._mouse_button = None

    def contextMenuEvent(self, event):
        """
        Handle right-click context menu.
        
        Provides options for 2D view manipulation:
        - Fit to screen: Auto-fit molecule to view
        - Reset view: Reset zoom and pan to default
        """
        if self._has_dragged_right:
            return

        if not self.molecule or not self.coords_2d:
            return
        
        menu = QMenu(self)
        
        # Fit to screen action
        fit_action = QAction("Fit to Screen", self)
        fit_action.setShortcut("F")
        fit_action.triggered.connect(self._auto_fit)
        menu.addAction(fit_action)
        
        # Reset view action
        reset_action = QAction("Reset View", self)
        reset_action.setShortcut("R")
        reset_action.triggered.connect(self._reset_view)
        menu.addAction(reset_action)
        
        menu.exec(event.globalPos())
    
    def _reset_view(self):
        """Reset view to default state (fit + center)."""
        self._auto_fit()
        self.update()

    def wheelEvent(self, event):
        pos = event.position()
        delta = event.angleDelta().y()
        factor = 1.12 if delta > 0 else 1.0 / 1.12

        # Zoom towards mouse cursor
        old_mx = (pos.x() - self._offset_x) / self._scale
        old_my = -(pos.y() - self._offset_y) / self._scale

        self._scale *= factor
        self._scale = max(10, min(250, self._scale))

        self._offset_x = pos.x() - old_mx * self._scale
        self._offset_y = pos.y() + old_my * self._scale

        self.update()

    def keyPressEvent(self, event):
        """Handle keyboard shortcuts.

        F      — Fit molecule to screen
        R      — Reset view
        Delete — Delete selected atoms
        """
        key = event.key()
        if key == Qt.Key.Key_Delete and self.selected_atoms:
            self.delete_requested.emit(set(self.selected_atoms))
        elif key == Qt.Key.Key_F:
            self._auto_fit()
            self.update()
        elif key == Qt.Key.Key_R:
            self._auto_fit()
            self.update()
        else:
            super().keyPressEvent(event)

    # ─── Rubber-band Helpers ──────────────────────────────────────

    def _hit_test_2d(self, pos):
        """
        Find the atom closest to the screen position *pos*, within a
        tolerance distance.  Returns the atom index or -1 if nothing
        is close enough.
        """
        if not self.molecule or not self.coords_2d:
            return -1

        best_idx = -1
        best_dist_sq = float('inf')
        # Tolerance in screen pixels
        tolerance = max(12, self._scale * 0.35)
        tol_sq = tolerance * tolerance

        for atom in self.molecule.atoms:
            if atom.index not in self.coords_2d:
                continue
            if atom.symbol == 'H' and not self.show_hydrogens:
                continue
            mx, my = self.coords_2d[atom.index]
            sx, sy = self._to_screen(mx, my)
            dx = pos.x() - sx
            dy = pos.y() - sy
            dist_sq = dx * dx + dy * dy
            if dist_sq < tol_sq and dist_sq < best_dist_sq:
                best_dist_sq = dist_sq
                best_idx = atom.index
        return best_idx

    def _commit_rubber_band_selection(self):
        """
        Finalise a Shift+drag selection: find all atoms whose 2D screen
        positions fall inside the rubber-band rectangle and add them to
        ``selected_atoms``.  Emits ``selection_changed``.
        """
        if not self._sel_rect_origin or not self._sel_rect_end:
            return

        x1 = min(self._sel_rect_origin.x(), self._sel_rect_end.x())
        y1 = min(self._sel_rect_origin.y(), self._sel_rect_end.y())
        x2 = max(self._sel_rect_origin.x(), self._sel_rect_end.x())
        y2 = max(self._sel_rect_origin.y(), self._sel_rect_end.y())

        # Tiny rectangle → clear selection
        if (x2 - x1) < 4 and (y2 - y1) < 4:
            self.selected_atoms.clear()
            self.selection_changed.emit(set())
            return

        rect = QRectF(x1, y1, x2 - x1, y2 - y1)
        newly_selected = set()

        visible = self._get_visible_atoms() if self.molecule else set()
        for idx in visible:
            if idx not in self.coords_2d:
                continue
            mx, my = self.coords_2d[idx]
            sx, sy = self._to_screen(mx, my)
            if rect.contains(QPointF(sx, sy)):
                newly_selected.add(idx)

        if newly_selected:
            self.selected_atoms |= newly_selected
        self.selection_changed.emit(set(self.selected_atoms))

    def _draw_rubber_band(self, painter):
        """Draw the semi-transparent selection rectangle overlay."""
        x1 = min(self._sel_rect_origin.x(), self._sel_rect_end.x())
        y1 = min(self._sel_rect_origin.y(), self._sel_rect_end.y())
        x2 = max(self._sel_rect_origin.x(), self._sel_rect_end.x())
        y2 = max(self._sel_rect_origin.y(), self._sel_rect_end.y())
        rect = QRectF(x1, y1, x2 - x1, y2 - y1)

        painter.setPen(QPen(QColor(255, 200, 50, 200), 1.5,
                            Qt.PenStyle.DashLine))
        painter.setBrush(QBrush(QColor(255, 200, 50, 40)))
        painter.drawRect(rect)

    def _rotate_2d_mouse(self, current_pos, last_pos):
        """Rotate 2D coordinates circularly based on mouse drag around centroid."""
        if not self.coords_2d:
            return
            
        # Get centroid in molecule coordinates
        xs = [c[0] for c in self.coords_2d.values() if c[0] is not None]
        ys = [c[1] for c in self.coords_2d.values() if c[1] is not None]
        if not xs or not ys:
            return
        cx_mol = sum(xs) / len(xs)
        cy_mol = sum(ys) / len(ys)
        
        # Convert centroid to screen coordinates
        cx_screen, cy_screen = self._to_screen(cx_mol, cy_mol)
        
        # Calculate angle difference
        ang1 = math.atan2(last_pos.y() - cy_screen, last_pos.x() - cx_screen)
        ang2 = math.atan2(current_pos.y() - cy_screen, current_pos.x() - cx_screen)
        angle_diff = ang2 - ang1
        
        # Handle wrap-around
        if angle_diff > math.pi:
            angle_diff -= 2 * math.pi
        elif angle_diff < -math.pi:
            angle_diff += 2 * math.pi
            
        # The screen Y grows down, molecule Y grows up.
        # Rotating on screen visually corresponds to rotating molecule coordinates
        # by the negative angle difference.
        angle = -angle_diff
        
        cos_a = math.cos(angle)
        sin_a = math.sin(angle)
        
        for idx in self.coords_2d:
            x, y = self.coords_2d[idx]
            if x is not None and y is not None:
                # Translate to origin
                tx = x - cx_mol
                ty = y - cy_mol
                # Rotate
                rx = tx * cos_a - ty * sin_a
                ry = tx * sin_a + ty * cos_a
                # Translate back
                self.coords_2d[idx] = [rx + cx_mol, ry + cy_mol]
                
                # Update cached atom coordinates
                atom = self.molecule.get_atom(idx)
                if atom:
                    atom.x2d = rx + cx_mol
                    atom.y2d = ry + cy_mol

    # ─── Utils ────────────────────────────────────────────────────

    def _auto_fit(self):
        if not self.coords_2d:
            return

        # Filter to visible atoms
        visible = self._get_visible_atoms() if self.molecule else set(self.coords_2d.keys())
        if not visible:
            visible = set(self.coords_2d.keys())

        visible_coords = {k: v for k, v in self.coords_2d.items() if k in visible}
        if not visible_coords:
            return

        xs = [c[0] for c in visible_coords.values()]
        ys = [c[1] for c in visible_coords.values()]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        span_x = max_x - min_x or 1
        span_y = max_y - min_y or 1
        cx = (min_x + max_x) / 2
        cy = (min_y + max_y) / 2

        w = self.width() or 600
        h = self.height() or 500
        margin = 70

        scale_x = (w - margin * 2) / span_x
        scale_y = (h - margin * 2) / span_y
        self._scale = min(scale_x, scale_y, 120)
        self._offset_x = w / 2 - cx * self._scale
        self._offset_y = h / 2 + cy * self._scale

    def _draw_protein_placeholder(self, painter, width, height):
        """Draw placeholder for large proteins instead of generating 2D coordinates."""
        # Set up font
        font = QFont('Segoe UI', 14)
        painter.setFont(font)

        # Draw protein icon representation
        cx, cy = width / 2, height / 2

        # Draw simplified protein representation
        painter.setPen(QPen(QColor(100, 100, 100), 2))

        # Draw a simple helix representation
        helix_width = 60
        helix_height = 120
        helix_x = cx - helix_width / 2
        helix_y = cy - helix_height / 2

        # Draw helix as a spiral
        for i in range(20):
            y = helix_y + i * (helix_height / 20)
            x_offset = 15 * math.sin(i * 0.5)
            painter.drawEllipse(QPointF(helix_x + helix_width/2 + x_offset, y), 3, 3)

        # Draw sheet representation
        sheet_x = cx + 40
        sheet_y = cy - 40
        arrow_points = [
            QPointF(sheet_x, sheet_y),
            QPointF(sheet_x + 40, sheet_y + 20),
            QPointF(sheet_x + 40, sheet_y + 10),
            QPointF(sheet_x + 50, sheet_y + 15),
            QPointF(sheet_x + 40, sheet_y + 20),
            QPointF(sheet_x + 40, sheet_y + 30),
            QPointF(sheet_x, sheet_y + 10)
        ]
        painter.setBrush(QBrush(QColor(50, 150, 220, 100)))
        painter.drawPolygon(QPolygonF(arrow_points))

        # Draw text information
        painter.setPen(QColor(60, 60, 60))
        title_font = QFont('Segoe UI', 16, QFont.Weight.Bold)
        painter.setFont(title_font)
        title = "Large Protein Structure"
        title_rect = QRectF(0, cy - 180, width, 30)
        painter.drawText(title_rect, Qt.AlignmentFlag.AlignCenter, title)

        # Draw stats
        info_font = QFont('Segoe UI', 12)
        painter.setFont(info_font)
        painter.setPen(QColor(80, 80, 80))

        num_atoms = len(self.molecule.atoms) if self.molecule else 0
        num_residues = len(set((a.res_seq, a.chain_id) for a in self.molecule.atoms
                              if hasattr(a, 'res_seq') and hasattr(a, 'chain_id'))) if self.molecule else 0

        info_lines = [
            f"Atoms: {num_atoms:,}",
            f"Residues: {num_residues:,}",
            "",
            "2D view skipped for large structures (>700 atoms)",
            "for better performance. Use 3D viewer instead."
        ]

        y_offset = cy + 80
        for line in info_lines:
            if line:
                rect = QRectF(0, y_offset, width, 25)
                painter.drawText(rect, Qt.AlignmentFlag.AlignCenter, line)
            y_offset += 25

        # Draw border
        painter.setPen(QPen(QColor(200, 200, 200), 1))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawRect(10, 10, width - 20, height - 20)

    def _optimize_2d_layout(self):
        """
        Optimize the 2D molecular layout using CoordinateGenerator2DSMILES.
        
        This method regenerates 2D coordinates using the same reliable generator
        that is used when initially loading molecules. The CoordinateGenerator2DSMILES
        properly handles the SMILES roundtrip and coordinate mapping.
        
        Returns:
            bool: True if optimization succeeded, False otherwise
        """
        if not self.molecule:
            return False
        
        try:
            from src.features.layout_2d.generators.coordgen2d_smiles_pure_oasa import CoordinateGenerator2DSMILES
            
            if _DEBUG:
                print("[DEBUG 2D] Starting 2D layout optimization with CoordinateGenerator2DSMILES...")
            
            # Use CoordinateGenerator2DSMILES with force_regenerate=True
            # This does the proper SMILES roundtrip and generates fresh coordinates
            generator = CoordinateGenerator2DSMILES(self.molecule, force_regenerate=True)
            new_coords = generator.generate()
            
            if not new_coords:
                if _DEBUG:
                    print("[DEBUG 2D] Coordinate generation failed")
                return False
            
            if _DEBUG:
                print(f"[DEBUG 2D] Generated {len(new_coords)} coordinates")
            
            # Update the coordinates
            self.coords_2d = new_coords
            
            # Update atom cached coordinates
            for idx, (x, y) in new_coords.items():
                atom = self.molecule.get_atom(idx)
                if atom:
                    atom.x2d = x
                    atom.y2d = y
                    atom.z2d = 0.0
            
            # Auto-fit and update
            self._auto_fit()
            self.update()
            
            if _DEBUG:
                print("[DEBUG 2D] Layout optimization complete")
            
            return True
            
        except Exception as e:
            if _DEBUG:
                print(f"[DEBUG 2D] Optimization error: {e}")
                import traceback
                traceback.print_exc()
            return False
    
    def _place_explicit_hydrogens_after_opt(self):
        """
        Place explicit hydrogen atoms around their parent heavy atoms
        after 2D coordinate optimization.
        
        This ensures H atoms (which are skipped by domain_to_oasa_mol)
        get reasonable positions near their bonded heavy atoms.
        """
        import math
        
        placed = 0
        BOND_LENGTH_H = 0.9  # Slightly shorter than heavy atom bonds
        
        for atom in self.molecule.atoms:
            # Skip if not hydrogen
            if atom.symbol != 'H':
                continue
            # Skip if already has coordinates
            if atom.index in self.coords_2d:
                continue
            
            # Find the heavy atom this H is bonded to
            parent_idx = None
            for bond in self.molecule.bonds:
                if bond.begin_atom_idx == atom.index and bond.end_atom_idx in self.coords_2d:
                    parent_idx = bond.end_atom_idx
                    break
                elif bond.end_atom_idx == atom.index and bond.begin_atom_idx in self.coords_2d:
                    parent_idx = bond.begin_atom_idx
                    break
            
            if parent_idx is not None:
                px, py = self.coords_2d[parent_idx]
                # Place H at a deterministic angle based on atom index
                # This distributes H atoms evenly around parent
                angle = (atom.index * 2.4) % (2 * math.pi)  # Golden angle approximation
                self.coords_2d[atom.index] = [
                    px + BOND_LENGTH_H * math.cos(angle),
                    py + BOND_LENGTH_H * math.sin(angle)
                ]
                placed += 1
        
        if _DEBUG and placed > 0:
            print(f"[DEBUG 2D] Placed {placed} explicit hydrogen atoms")
    
    def _fix_overlaps_with_repulsion(self):
        """
        Fix atom overlaps using gentle repulsion forces.
        
        After OASA coordinate generation, some atoms may still overlap
        (especially in complex multi-ring systems). This method applies
        a simple repulsion-based algorithm to push overlapping atoms apart
        while preserving the overall layout topology.
        
        The algorithm:
        1. Detects pairs of atoms that are too close (distance < threshold)
        2. Applies gentle repulsion forces to push them apart
        3. Runs a few iterations to resolve all overlaps
        4. Maintains bonded atom distances to preserve structure
        """
        import math
        
        if not self.coords_2d or len(self.coords_2d) < 2:
            return
        
        # Get list of bonded pairs (to preserve bond distances)
        bonded_pairs = set()
        for bond in self.molecule.bonds:
            idx1, idx2 = bond.begin_atom_idx, bond.end_atom_idx
            if idx1 in self.coords_2d and idx2 in self.coords_2d:
                bonded_pairs.add((min(idx1, idx2), max(idx1, idx2)))
        
        # Minimum acceptable distance between non-bonded atoms
        # Increased to 1.0 to ensure atoms are clearly separated (one bond length)
        MIN_DISTANCE = 1.0  # OASA units (bond length is typically ~1.0)
        # Repulsion strength - higher to push atoms apart more aggressively
        REPULSION_STRENGTH = 0.25
        # Number of iterations
        MAX_ITERATIONS = 100
        # Early stopping if no overlaps
        CONVERGENCE_THRESHOLD = 0.005
        
        for iteration in range(MAX_ITERATIONS):
            # Calculate repulsion forces for each atom
            forces = {idx: [0.0, 0.0] for idx in self.coords_2d}
            max_displacement = 0.0
            overlaps_found = 0
            
            # Check all pairs
            indices = list(self.coords_2d.keys())
            for i in range(len(indices)):
                for j in range(i + 1, len(indices)):
                    idx1, idx2 = indices[i], indices[j]
                    x1, y1 = self.coords_2d[idx1]
                    x2, y2 = self.coords_2d[idx2]
                    
                    dx = x2 - x1
                    dy = y2 - y1
                    dist_sq = dx * dx + dy * dy
                    dist = math.sqrt(dist_sq)
                    
                    # Skip if already well-separated
                    if dist > MIN_DISTANCE:
                        continue
                    
                    # Skip bonded pairs (preserve bond geometry)
                    if (min(idx1, idx2), max(idx1, idx2)) in bonded_pairs:
                        continue
                    
                    overlaps_found += 1
                    
                    # Calculate repulsion (stronger when closer)
                    if dist < 0.001:
                        # Avoid division by zero - push in random direction
                        angle = (idx1 + idx2) * 0.7
                        dx = math.cos(angle)
                        dy = math.sin(angle)
                        dist = 0.001
                    
                    # Repulsion force inversely proportional to distance
                    # Use quadratic falloff for stronger repulsion at close range
                    overlap = MIN_DISTANCE - dist
                    force = REPULSION_STRENGTH * overlap * (1.0 + overlap / dist)
                    
                    # Normalize direction
                    fx = dx / dist * force
                    fy = dy / dist * force
                    
                    # Apply forces (opposite directions)
                    forces[idx1][0] -= fx
                    forces[idx1][1] -= fy
                    forces[idx2][0] += fx
                    forces[idx2][1] += fy
            
            # Also apply gentle attraction to bonded pairs to maintain bond lengths
            # This prevents the molecule from becoming too stretched
            TARGET_BOND_LENGTH = 1.0
            BOND_ATTRACTION_STRENGTH = 0.05
            
            for idx1, idx2 in bonded_pairs:
                if idx1 not in self.coords_2d or idx2 not in self.coords_2d:
                    continue
                x1, y1 = self.coords_2d[idx1]
                x2, y2 = self.coords_2d[idx2]
                
                dx = x2 - x1
                dy = y2 - y1
                dist = math.sqrt(dx * dx + dy * dy)
                
                if dist < 0.001:
                    continue
                
                # Attraction if bond is too long, repulsion if too short
                deviation = dist - TARGET_BOND_LENGTH
                force = BOND_ATTRACTION_STRENGTH * deviation
                
                fx = dx / dist * force
                fy = dy / dist * force
                
                forces[idx1][0] += fx
                forces[idx1][1] += fy
                forces[idx2][0] -= fx
                forces[idx2][1] -= fy
            
            # Apply forces to positions
            for idx in self.coords_2d:
                fx, fy = forces[idx]
                if abs(fx) > max_displacement:
                    max_displacement = abs(fx)
                if abs(fy) > max_displacement:
                    max_displacement = abs(fy)
                
                x, y = self.coords_2d[idx]
                self.coords_2d[idx] = [x + fx, y + fy]
            
            if _DEBUG and iteration % 10 == 0:
                print(f"[DEBUG 2D] Overlap fix iter {iteration}: {overlaps_found} overlaps, max_disp={max_displacement:.4f}")
            
            # Check convergence
            if max_displacement < CONVERGENCE_THRESHOLD:
                if _DEBUG:
                    print(f"[DEBUG 2D] Overlap fix converged at iteration {iteration}")
                break
        
        if _DEBUG:
            print(f"[DEBUG 2D] Overlap fix complete after {iteration + 1} iterations")
    
    def _center_coords_2d(self):
        """
        Center the 2D coordinates around the origin.
        
        This maintains consistent molecule positioning after optimization
        by recentering the coordinate system.
        """
        if not self.coords_2d:
            return
        
        # Calculate centroid
        xs = [c[0] for c in self.coords_2d.values() if c[0] is not None]
        ys = [c[1] for c in self.coords_2d.values() if c[1] is not None]
        
        if not xs or not ys:
            return
        
        cx = sum(xs) / len(xs)
        cy = sum(ys) / len(ys)
        
        # Center all coordinates
        for idx in self.coords_2d:
            x, y = self.coords_2d[idx]
            if x is not None and y is not None:
                self.coords_2d[idx] = [x - cx, y - cy]
