"""
3D Molecular Viewer — High-quality software-rendered molecular visualization.

Features:
- Smooth radial-gradient sphere rendering (realistic 3D look)
- Ball-and-stick / spacefill / wireframe rendering modes
- CPK coloring
- Mouse rotation, zoom, and pan
- Atom highlighting on hover
- High-DPI image export

Rendering is done directly with QPainter in this widget's paintEvent.
The PainterRenderer helper holds the gradient cache, off-screen culling,
and LOD logic for performance on large molecules.
"""

import logging
import math
import numpy as np
from src.shared.qt_compat import QWidget, Qt, QTimer, Signal, QPointF, QRectF
from src.shared.qt_compat import (
    QPainter, QColor, QPen, QBrush, QFont, QWheelEvent,
    QRadialGradient, QLinearGradient, QImage, QConicalGradient, QPainterPath,
    QMenu, QAction
)
from src.shared.ui.theme import COLORS

logger = logging.getLogger(__name__)


def _hex_to_rgb(hex_color):
    """Convert hex color string to (r, g, b) tuple."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))


# Atom radii for display (scaled for visual appeal)
DISPLAY_RADIUS = {
    'H': 0.25, 'He': 0.31, 'C': 0.40, 'N': 0.38, 'O': 0.36, 'F': 0.32,
    'P': 0.44, 'S': 0.42, 'Cl': 0.39, 'Br': 0.41, 'I': 0.44, 'B': 0.38,
    'Si': 0.44, 'Se': 0.42, 'Na': 0.50, 'K': 0.55, 'Ca': 0.48, 'Fe': 0.44,
}


class SoftwareMolViewer3D(QWidget):
    """
    Software-rendered 3D molecular viewer with mouse interaction.
    Uses QPainter with QRadialGradient for smooth, realistic sphere rendering.

    Selection
    ---------
    **Shift + left-drag** draws a rubber-band rectangle (PyMOL-style).
    Atoms whose projected screen positions fall inside the rectangle on
    mouse-release are added to ``selected_atoms``.  A plain left-click on
    empty space clears the selection.

    Deletion
    --------
    Pressing the **Delete** key while atoms are selected emits
    ``delete_requested`` so the main window can remove those atoms from
    the domain model and refresh both viewers.
    """

    # --- Signals ---
    atom_hovered = Signal(int)
    atom_clicked = Signal(int)
    selection_changed = Signal(object)   # emits set of selected atom indices
    delete_requested = Signal(object)    # emits set of atom indices to delete

    def __init__(self, parent=None):
        super().__init__(parent)
        self.molecule = None
        self.setMinimumSize(400, 400)
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)  # receive key events

        # Camera state
        self.rot_x = 20.0
        self.rot_y = -30.0
        self.rot_z = 0.0
        self.pan_x = 0.0
        self.pan_y = 0.0
        self.zoom = 40.0
        self.auto_scale = True

        # Mouse interaction
        self._last_mouse_pos = None
        self._mouse_button = None
        self._hovered_atom = -1
        self.selected_atoms = set()  # Set of atom indices to highlight

        # Rubber-band selection rectangle (screen coords, set during Shift+drag)
        self._sel_rect_origin = None   # QPointF or None
        self._sel_rect_end = None      # QPointF or None
        self._is_selecting = False     # True while Shift+left-drag is active

        # Measurement state
        self._measure_atoms = []   # List of picked atoms for distance/angle
        self._measurements = []    # List of completed measurements

        # Protein/Ligand display state
        self.show_ligands_in_cartoon = True
        self.visible_sidechains = set()  # Set of residue sequence numbers
        self.interaction_lines = []      # List of (atom1_idx, atom2_idx, type, color)
        self.labels = {}                 # atom_idx -> text

        # Rendering options
        self.render_mode = 'ball_and_stick'

        # Animation
        self._auto_rotate = False
        self._rotation_timer = QTimer(self)

        # Rendering settings
        self.show_hydrogens = True
        self.show_labels = False
        self.show_sidechains = False
        self.render_mode = 'ball_and_stick'  # 'spacefill', 'wireframe', 'cartoon', 'ribbon', 'backbone'
        self.custom_atom_modes = {}
        self.custom_atom_colors = {}  # atom_idx -> (r, g, b) tuple
        self.sidechain_res_vis = {}
        self.labeled_residues = {}  # mapping res_seq to QColor
        self.use_ssao = False  # Fake real-time ray-tracing toggle
        self.use_gouraud = False  # Gouraud normal smoothing toggle
        self.bg_color = QColor(COLORS['viewer_bg'])

        # User-adjustable radius scales (1.0 = default)
        self.sphere_scale = 0.6   # Multiplier for atom sphere radius (60% default)
        self.stick_scale = 1.0    # Multiplier for bond stick width
        self.line_scale = 1.0     # Multiplier for wireframe line width
        self.label_font_size = 9  # Fixed label font size in points
        self.label_color = QColor(255, 255, 255, 230)  # Default white labels

        # Light direction (normalized) — top-left-front
        self._light_dir = np.array([-0.4, -0.5, 1.0])
        self._light_dir = self._light_dir / np.linalg.norm(self._light_dir)

        # --- Delegates ---
        from src.features.visualization_3d.ui.painter_renderer import PainterRenderer
        from src.features.visualization_3d.ui.mouse_controller import MouseController

        self._renderer = PainterRenderer()
        self._mouse_ctrl = MouseController(self)

        self._rotation_timer.timeout.connect(self._mouse_ctrl.auto_rotate_step)

    # ─── Molecule Loading ─────────────────────────────────────────

    def set_molecule(self, molecule):
        self.molecule = molecule
        self.custom_atom_colors = {}
        self.labels = {}  # Clear custom labels when loading new molecule
        if molecule and len(molecule.atoms) > 0:
            self._auto_fit()
            # Auto-switch to cartoon for proteins, reset for small molecules
            is_protein = getattr(molecule, 'properties', {}).get('is_protein', False)
            if is_protein:
                self.render_mode = 'cartoon'
            else:
                self.render_mode = 'ball_and_stick'
        self.update()

    def clear(self):
        self.molecule = None
        self.labels = {}
        self.update()

    def clear_labels(self):
        """Universal function to clear all custom atom labels."""
        self.labels = {}
        self.show_labels = False
        self.update()

    def toggle_auto_rotate(self):
        self._auto_rotate = not self._auto_rotate
        if self._auto_rotate:
            self._rotation_timer.start(33)
        else:
            self._rotation_timer.stop()

    def reset_view(self):
        self.rot_x = 20.0
        self.rot_y = -30.0
        self.rot_z = 0.0
        self.pan_x = 0.0
        self.pan_y = 0.0
        if self.molecule:
            self._auto_fit()
        self.update()

    # ─── Rendering ─────────────────────────────────────────────────

    def paintEvent(self, event):
        painter = QPainter(self)
        self._renderer.render(self, painter, self.width(), self.height())
        painter.end()

    # Keep thin wrappers so any external code calling these still works.

    def _render(self, painter, width, height, is_export=False, export_scale=1.0):
        """Core rendering logic — used by both paintEvent and export."""
        self._renderer.render(self, painter, width, height, is_export, export_scale)

    def _project_atoms(self, vp_width=None, vp_height=None):
        """Project 3D atom coordinates to 2D screen coordinates."""
        return self._renderer._project_atoms(self, vp_width, vp_height)

    def _draw_atom_sphere(self, painter, atom_idx, sx, sy, sz, radius, rgb, alpha=1.0):
        self._renderer._draw_atom_sphere(self, painter, atom_idx, sx, sy, sz, radius, rgb, alpha)

    def _draw_selection_ring(self, painter, sx, sy, radius):
        self._renderer._draw_selection_ring(painter, sx, sy, radius)

    def _draw_sasa_surface(self, painter, width, height):
        self._renderer._draw_sasa_surface(self, painter, width, height)

    def set_selected(self, atom_indices):
        """Set which atoms are highlighted (from console select commands)."""
        self.selected_atoms = set(atom_indices)
        self.update()

    def set_atom_colors(self, colors_dict):
        """Set custom colors for specific atoms.
        
        Args:
            colors_dict: Dictionary mapping atom index to (r, g, b) tuple or None to reset.
        """
        if colors_dict is None:
            self.custom_atom_colors = {}
        else:
            self.custom_atom_colors.update(colors_dict)
        
        if hasattr(self, '_renderer'):
            self._renderer.invalidate_cache()
        self.update()

    def _draw_bonds(self, painter, projected, custom_only=False):
        self._renderer._draw_bonds(self, painter, projected, custom_only)

    def _draw_bond_line(self, painter, x1, y1, x2, y2, c1, c2, width, shade, dashed=False, is_custom=False):
        self._renderer._draw_bond_line(painter, x1, y1, x2, y2, c1, c2, width, shade, dashed, is_custom)

    def _draw_label(self, painter, atom_idx, sx, sy, radius):
        self._renderer._draw_label(self, painter, atom_idx, sx, sy, radius)

    def _draw_residue_label(self, painter, text, sx, sy, color, radius):
        self._renderer._draw_residue_label(self, painter, text, sx, sy, color, radius)

    def _draw_overlay(self, painter):
        self._renderer._draw_overlay(self, painter)

    def _draw_placeholder(self, painter, width, height):
        self._renderer._draw_placeholder(painter, width, height)

    # ─── Image Export ─────────────────────────────────────────────

    def export_image(self, filepath, dpi=300, bg_white=True):
        """
        Export the current view as a high-resolution image.

        Args:
            filepath: Output file path (.png, .jpg, .tiff, .bmp)
            dpi: Resolution in dots per inch (72, 150, 300, 600, etc.)
            bg_white: If True, use white background instead of dark

        Returns:
            True if successful
        """
        # Calculate pixel dimensions from current widget size and DPI
        scale_factor = dpi / 96.0  # 96 DPI is the default screen DPI
        img_width = int(self.width() * scale_factor)
        img_height = int(self.height() * scale_factor)

        # Create high-res image
        image = QImage(img_width, img_height, QImage.Format.Format_ARGB32_Premultiplied)
        image.setDotsPerMeterX(int(dpi / 0.0254))
        image.setDotsPerMeterY(int(dpi / 0.0254))

        # Save and modify state for export
        original_bg = self.bg_color
        original_zoom = self.zoom
        original_pan_x = self.pan_x
        original_pan_y = self.pan_y

        if bg_white:
            self.bg_color = QColor(255, 255, 255)

        # Scale zoom and pan for higher resolution
        self.zoom *= scale_factor
        self.pan_x *= scale_factor
        self.pan_y *= scale_factor

        # Render into image
        painter = QPainter(image)
        self._render(painter, img_width, img_height, is_export=True, export_scale=scale_factor)
        painter.end()

        # Restore state
        self.bg_color = original_bg
        self.zoom = original_zoom
        self.pan_x = original_pan_x
        self.pan_y = original_pan_y

        # Save
        result = image.save(filepath)
        return result

    # ─── Mouse Interaction ────────────────────────────────────────

    def mousePressEvent(self, event):
        self._mouse_ctrl.handle_mouse_press(event)

    def mouseMoveEvent(self, event):
        self._mouse_ctrl.handle_mouse_move(event)

    def mouseReleaseEvent(self, event):
        self._mouse_ctrl.handle_mouse_release(event)

    def wheelEvent(self, event: QWheelEvent):
        self._mouse_ctrl.handle_wheel(event)

    def keyPressEvent(self, event):
        if not self._mouse_ctrl.handle_key_press(event):
            super().keyPressEvent(event)

    def contextMenuEvent(self, event):
        """Show context menu for selected atoms/residues and general viewer options."""
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

        ssao_action = menu.addAction("Disable Fake Ray-Tracing" if self.use_ssao else "Enable Fake Ray-Tracing (SSAO)")
        gouraud_action = menu.addAction("Disable Smooth Shading" if self.use_gouraud else "Enable Smooth Shading (Gouraud)")

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
            self.update()
        elif hide_sc_action and action == hide_sc_action:
            for rs in selected_res_seqs:
                self.sidechain_res_vis[rs] = False
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
        elif action == ssao_action:
            self.use_ssao = not self.use_ssao
            self.update()
        elif action == gouraud_action:
            self.use_gouraud = not self.use_gouraud
            self.update()

    # ─── Rubber-band Helpers ──────────────────────────────────────

    def _commit_rubber_band_selection(self):
        self._mouse_ctrl._commit_rubber_band_selection()

    def _draw_rubber_band(self, painter):
        self._renderer._draw_rubber_band(self, painter)

    def _detect_hover(self, pos):
        self._mouse_ctrl._detect_hover(pos)

    def _hit_test(self, pos):
        return self._mouse_ctrl._hit_test(pos)

    # ─── Camera Utilities ─────────────────────────────────────────

    def _auto_fit(self):
        if not self.molecule or not self.molecule.atoms:
            return

        coords = []
        for atom in self.molecule.atoms:
            if atom.has_coords:
                coords.append([atom.x, atom.y, atom.z])

        if not coords:
            return

        coords = np.array(coords)
        span = np.max(coords, axis=0) - np.min(coords, axis=0)
        max_span = max(span) if max(span) > 0 else 1.0

        viewport_size = min(self.width(), self.height())
        self.zoom = min(100, max(10, viewport_size * 0.3 / max_span))
        self.pan_x = 0
        self.pan_y = 0

    def focus_on_atoms(self, atom_indices, padding_angstroms=0.0):
        """Center the view on the given atoms and zoom in."""
        if not self.molecule or not atom_indices:
            return

        coords = []
        for idx in atom_indices:
            if idx < len(self.molecule.atoms):
                atom = self.molecule.atoms[idx]
                if atom.has_coords:
                    coords.append([atom.x, atom.y, atom.z])

        if not coords:
            return

        coords = np.array(coords)
        centroid = np.mean(coords, axis=0)
        span = np.max(coords, axis=0) - np.min(coords, axis=0)
        max_span = max(span) if max(span) > 0 else 1.0
        
        # Add padding to the effective span to ensure surrounding area is visible
        max_span += 2.0 * padding_angstroms

        viewport_size = min(self.width(), self.height())
        # Zoom tighter than auto_fit (0.4 vs 0.3) but capped at 100
        self.zoom = min(100, max(15, viewport_size * 0.4 / max_span))

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

    def _auto_rotate_step(self):
        self._mouse_ctrl.auto_rotate_step()

    # ─── Measurements ─────────────────────────────────────────────

    def _complete_distance_measurement(self):
        self._mouse_ctrl._complete_distance_measurement()

    def _complete_angle_measurement(self):
        self._mouse_ctrl._complete_angle_measurement()

    def _draw_measurements(self, painter, projected):
        self._renderer._draw_measurements(self, painter, projected)

    # ─── Protein Rendering ─────────────────────────────────────────────

    def _draw_protein(self, painter, projected, width, height):
        self._renderer._draw_protein(self, painter, projected, width, height)

    def _group_residues(self):
        return self._renderer._group_residues(self)

    def _get_ss_color(self, ss_type):
        return self._renderer._get_ss_color(ss_type)

    def _draw_cartoon(self, painter, residues):
        self._renderer._draw_cartoon(self, painter, residues)

    def _draw_pyMOL_cartoon_chain(self, painter, points):
        self._renderer._draw_pyMOL_cartoon_chain(painter, points)

    def _draw_pyMOL_helix(self, painter, points):
        self._renderer._draw_pyMOL_helix(painter, points)

    def _draw_pyMOL_sheet(self, painter, points):
        self._renderer._draw_pyMOL_sheet(painter, points)

    def _draw_pyMOL_coil(self, painter, points):
        self._renderer._draw_pyMOL_coil(painter, points)

    def _draw_ribbon(self, painter, residues):
        self._renderer._draw_ribbon(self, painter, residues)

    def _draw_backbone(self, painter, residues):
        self._renderer._draw_backbone(self, painter, residues)

    def _draw_side_chains(self, painter, projected):
        self._renderer._draw_side_chains(self, painter, projected)

    def _draw_cylinder(self, painter, x1, y1, x2, y2, color, width):
        self._renderer._draw_cylinder(painter, x1, y1, x2, y2, color, width)

    def _draw_arrow(self, painter, x1, y1, x2, y2, color, width):
        self._renderer._draw_arrow(painter, x1, y1, x2, y2, color, width)

    def _draw_tube(self, painter, x1, y1, x2, y2, color, width):
        self._renderer._draw_tube(painter, x1, y1, x2, y2, color, width)

    def _draw_smooth_ribbon(self, painter, x1, y1, x2, y2, color, width):
        self._renderer._draw_smooth_ribbon(painter, x1, y1, x2, y2, color, width)

    def _draw_large_molecule_fast(self, painter, projected, sorted_atoms):
        self._renderer._draw_large_molecule_fast(self, painter, projected, sorted_atoms)

    def _draw_atom_simple(self, painter, sx, sy, radius, color):
        self._renderer._draw_atom_simple(painter, sx, sy, radius, color)

    def _draw_performance_indicator(self, painter, num_atoms):
        self._renderer._draw_performance_indicator(self, painter, num_atoms)

    def _draw_dummy_spheres(self, painter, width, height):
        self._renderer._draw_dummy_spheres(self, painter, width, height)


class MolViewer3D(QWidget):
    # Public Signals (forwarding events from active child)
    atom_hovered = Signal(int)
    atom_clicked = Signal(int)
    selection_changed = Signal(object)
    delete_requested = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        from src.shared.qt_compat import QStackedWidget, QVBoxLayout
        
        # Create layouts
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.stacked = QStackedWidget(self)
        layout.addWidget(self.stacked)
        
        # Instantiate children
        self.software_viewer = SoftwareMolViewer3D(self)
        
        # Instantiate GL widget (safe if GL fails)
        from src.features.visualization_3d.ui.gl_widget import GLMoleculeWidget
        self.gl_viewer = GLMoleculeWidget(self)
        
        # Add to stack
        self.stacked.addWidget(self.software_viewer)
        self.stacked.addWidget(self.gl_viewer)
        
        # Connect signals
        self.software_viewer.atom_hovered.connect(self.atom_hovered.emit)
        self.software_viewer.atom_clicked.connect(self.atom_clicked.emit)
        self.software_viewer.selection_changed.connect(self.selection_changed.emit)
        self.software_viewer.delete_requested.connect(self.delete_requested.emit)
        
        self.gl_viewer.atom_hovered.connect(self.atom_hovered.emit)
        self.gl_viewer.atom_clicked.connect(self.atom_clicked.emit)
        self.gl_viewer.selection_changed.connect(self.selection_changed.emit)
        self.gl_viewer.delete_requested.connect(self.delete_requested.emit)
        
        # Current active viewer
        self.active_viewer = self.software_viewer
        self.stacked.setCurrentWidget(self.software_viewer)
        
    @property
    def molecule(self):
        return self.active_viewer.molecule

    @molecule.setter
    def molecule(self, value):
        self.software_viewer.molecule = value
        self.gl_viewer.molecule = value

    def rebuild_mesh(self):
        """Force the GL viewer to rebuild its mesh buffers (e.g. after color changes)."""
        if hasattr(self.gl_viewer, '_update_gl_buffers') and self.gl_viewer.gl_available:
            self.gl_viewer._update_gl_buffers()
        self.update()

    def set_selected(self, atom_indices):
        self.software_viewer.set_selected(atom_indices)
        if hasattr(self.gl_viewer, 'set_selected'):
            self.gl_viewer.set_selected(atom_indices)
        else:
            self.gl_viewer.selected_atoms = set(atom_indices)
            self.gl_viewer.update()

    def set_atom_colors(self, colors_dict):
        self.software_viewer.set_atom_colors(colors_dict)
        if hasattr(self.gl_viewer, 'set_atom_colors'):
            self.gl_viewer.set_atom_colors(colors_dict)

    @property
    def selected_atoms(self):
        return self.active_viewer.selected_atoms

    @selected_atoms.setter
    def selected_atoms(self, value):
        self.software_viewer.selected_atoms = value
        # GLWidget doesn't support selected_atoms natively, but we mirror it
        if hasattr(self.gl_viewer, 'selected_atoms'):
            self.gl_viewer.selected_atoms = value

    # Mirror camera and settings
    @property
    def rot_x(self): return self.active_viewer.rot_x
    @rot_x.setter
    def rot_x(self, val):
        self.software_viewer.rot_x = val
        self.gl_viewer.rot_x = val
        
    @property
    def rot_y(self): return self.active_viewer.rot_y
    @rot_y.setter
    def rot_y(self, val):
        self.software_viewer.rot_y = val
        self.gl_viewer.rot_y = val

    @property
    def rot_z(self): return self.active_viewer.rot_z
    @rot_z.setter
    def rot_z(self, val):
        self.software_viewer.rot_z = val
        self.gl_viewer.rot_z = val

    @property
    def pan_x(self): return self.active_viewer.pan_x
    @property
    def pan_y(self): return self.active_viewer.pan_y
    @property
    def zoom(self): return self.active_viewer.zoom

    @pan_x.setter
    def pan_x(self, val):
        self.software_viewer.pan_x = val
        self.gl_viewer.pan_x = val

    @pan_y.setter
    def pan_y(self, val):
        self.software_viewer.pan_y = val
        self.gl_viewer.pan_y = val

    @zoom.setter
    def zoom(self, val):
        self.software_viewer.zoom = val
        self.gl_viewer.zoom = val

    @property
    def sphere_scale(self): return self.software_viewer.sphere_scale
    @sphere_scale.setter
    def sphere_scale(self, val):
        self.software_viewer.sphere_scale = val
        self.gl_viewer.sphere_scale = val

    @property
    def stick_scale(self): return self.software_viewer.stick_scale
    @stick_scale.setter
    def stick_scale(self, val):
        self.software_viewer.stick_scale = val
        self.gl_viewer.stick_scale = val

    @property
    def line_scale(self): return self.software_viewer.line_scale
    @line_scale.setter
    def line_scale(self, val):
        self.software_viewer.line_scale = val

    @property
    def show_hydrogens(self): return self.software_viewer.show_hydrogens
    @show_hydrogens.setter
    def show_hydrogens(self, val):
        self.software_viewer.show_hydrogens = val
        self.gl_viewer.show_hydrogens = val

    @property
    def render_mode(self): return self.software_viewer.render_mode
    @render_mode.setter
    def render_mode(self, val):
        self.software_viewer.render_mode = val
        self.gl_viewer.render_mode = val

    @property
    def use_ssao(self): return getattr(self.software_viewer, 'use_ssao', False)
    @use_ssao.setter
    def use_ssao(self, val):
        self.software_viewer.use_ssao = val

    @property
    def use_gouraud(self): return getattr(self.software_viewer, 'use_gouraud', False)
    @use_gouraud.setter
    def use_gouraud(self, val):
        self.software_viewer.use_gouraud = val

    @property
    def bg_color(self): return self.software_viewer.bg_color
    @bg_color.setter
    def bg_color(self, val):
        self.software_viewer.bg_color = val
        self.gl_viewer.bg_color = val

    # Method delegation
    def set_molecule(self, molecule):
        # Decide which viewer to use
        from src.services.rendering.renderer_factory import RendererFactory
        factory = RendererFactory()
        
        # Check if we should use GL for this molecule
        use_gl = False
        if molecule:
            # Debug logging
            is_protein = getattr(molecule, 'properties', {}).get('is_protein', False)
            n_atoms = len(molecule.atoms) if hasattr(molecule, 'atoms') else 0
            print(f"[DEBUG] set_molecule: is_protein={is_protein}, n_atoms={n_atoms}")
            
            # Check size threshold
            if factory.should_use_gl(molecule):
                # Ensure the GL widget actually succeeded in initialisation
                if factory.check_gl_available(self.gl_viewer):
                    use_gl = True
                    print(f"[DEBUG] OpenGL available and threshold met, using GL")
                else:
                    print(f"[DEBUG] OpenGL not available (gl_available=False)")
            else:
                print(f"[DEBUG] Threshold not met for GL")
        
        if use_gl:
            # Sync camera parameters from software to GL before switching
            if self.active_viewer == self.software_viewer:
                self.gl_viewer.rot_x = self.software_viewer.rot_x
                self.gl_viewer.rot_y = self.software_viewer.rot_y
                self.gl_viewer.rot_z = self.software_viewer.rot_z
                self.gl_viewer.pan_x = self.software_viewer.pan_x
                self.gl_viewer.pan_y = self.software_viewer.pan_y
                self.gl_viewer.zoom = self.software_viewer.zoom
                
            self.active_viewer = self.gl_viewer
            self.stacked.setCurrentWidget(self.gl_viewer)
            self.gl_viewer.set_molecule(molecule)
            print("[Viewer3D] Dynamic switch: Using accelerated OpenGL shaders (GPU)")
        else:
            # Sync camera parameters from GL to software before switching
            if self.active_viewer == self.gl_viewer:
                self.software_viewer.rot_x = self.gl_viewer.rot_x
                self.software_viewer.rot_y = self.gl_viewer.rot_y
                self.software_viewer.rot_z = self.gl_viewer.rot_z
                self.software_viewer.pan_x = self.gl_viewer.pan_x
                self.software_viewer.pan_y = self.gl_viewer.pan_y
                self.software_viewer.zoom = self.gl_viewer.zoom
                
            self.active_viewer = self.software_viewer
            self.stacked.setCurrentWidget(self.software_viewer)
            self.software_viewer.set_molecule(molecule)
            print("[Viewer3D] Dynamic switch: Using QPainter software engine (CPU)")

    def clear(self):
        self.software_viewer.clear()
        self.gl_viewer.clear()

    def clear_labels(self):
        self.software_viewer.clear_labels()

    def toggle_auto_rotate(self):
        self.software_viewer.toggle_auto_rotate()

    def reset_view(self):
        self.active_viewer.reset_view()

    def focus_on_atoms(self, atom_indices, padding_angstroms=0.0):
        if self.active_viewer == self.software_viewer:
            self.software_viewer.focus_on_atoms(atom_indices, padding_angstroms)
        else:
            self.gl_viewer._auto_fit()  # Fallback for GL

    def export_image(self, filepath, dpi=300, bg_white=True):
        if self.active_viewer == self.software_viewer:
            return self.software_viewer.export_image(filepath, dpi, bg_white)
        else:
            # GL fallback
            return self.gl_viewer.grabGeometry().save(filepath)

    def update(self):
        self.active_viewer.update()
        super().update()
