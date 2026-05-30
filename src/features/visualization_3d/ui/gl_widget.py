"""
OpenGL-accelerated 3D Molecular Viewer — Hardware-rendered path for large molecules.

Provides a QOpenGLWidget-based renderer that can handle 500+ atoms efficiently
using hardware acceleration. Falls back gracefully when OpenGL is unavailable
or insufficient (< 3.3).

Cross-platform notes
--------------------
- macOS: Apple deprecated OpenGL at 4.1 max. We request 3.3 Core Profile.
- Windows: Older Intel iGPUs may lack 3.3. We detect and report via gl_available.
- Linux: Mesa driver quality varies. All GL errors are caught and reported.

If anything goes wrong during GL initialisation the widget sets
``self.gl_available = False`` so the factory layer can fall back to QPainter
rendering without a crash.
"""

import math
import sys
import traceback
import ctypes

import numpy as np

from src.shared.qt_compat import (
    Qt, QColor, QPainter, QPen, QBrush, QFont, QPointF, QRectF,
    QRadialGradient, QWheelEvent, QVector3D, QVector4D, QMatrix4x4,
    QOpenGLShaderProgram, QOpenGLShader, QOpenGLBuffer, Signal
)
from src.shared.ui.theme import COLORS
from src.features.visualization_3d.ui.shaders import (
    MESH_VERTEX_SHADER, MESH_FRAGMENT_SHADER,
    SPHERE_VERTEX_SHADER, SPHERE_FRAGMENT_SHADER,
    LINE_VERTEX_SHADER, LINE_FRAGMENT_SHADER
)
import time


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

def _hex_to_rgb(hex_color: str):
    """Convert hex colour string to (r, g, b) tuple (0-255)."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))


def _hex_to_rgb_float(hex_color: str):
    """Convert hex colour string to (r, g, b) tuple in 0-1 range."""
    r, g, b = _hex_to_rgb(hex_color)
    return (r / 255.0, g / 255.0, b / 255.0)


# Atom radii for display (scaled for visual appeal) — mirrors painter_renderer
DISPLAY_RADIUS = {
    'H': 0.25, 'He': 0.31, 'C': 0.40, 'N': 0.38, 'O': 0.36, 'F': 0.32,
    'P': 0.44, 'S': 0.42, 'Cl': 0.39, 'Br': 0.41, 'I': 0.44, 'B': 0.38,
    'Si': 0.44, 'Se': 0.42, 'Na': 0.50, 'K': 0.55, 'Ca': 0.48, 'Fe': 0.44,
}


def _element_color(symbol: str) -> tuple:
    """Return (r, g, b) in 0-255 for *symbol* using the Element colour table.

    Falls back to grey if the element lookup fails.
    """
    try:
        from src.core.domain.models.elements import get_element
        elem = get_element(symbol)
        if elem and elem.color:
            return _hex_to_rgb(elem.color)
    except Exception:
        pass
    return (180, 180, 180)


def _element_color_float(symbol: str) -> tuple:
    """Return (r, g, b) in 0.0-1.0 for *symbol*."""
    r, g, b = _element_color(symbol)
    return (r / 255.0, g / 255.0, b / 255.0)


def _display_radius(symbol: str) -> float:
    """Return display radius for *symbol* (Angstroms, aesthetic scale)."""
    return DISPLAY_RADIUS.get(symbol, 0.35)


# ---------------------------------------------------------------------------
#  Try to import QOpenGLWidget — it lives in different sub-packages depending
#  on the PySide6 / PyQt6 version.
# ---------------------------------------------------------------------------

_QOpenGLWidget = None

try:
    from PySide6.QtOpenGLWidgets import QOpenGLWidget as _QOpenGLWidget
except ImportError:
    try:
        from PyQt6.QtOpenGLWidgets import QOpenGLWidget as _QOpenGLWidget
    except ImportError:
        pass

if _QOpenGLWidget is None:
    # Neither framework provides QOpenGLWidget — provide a plain QWidget stub
    # so that importing this module never crashes.
    from src.shared.qt_compat import QWidget as _QOpenGLWidget  # type: ignore[assignment]


# ---------------------------------------------------------------------------
#  GLMoleculeWidget
# ---------------------------------------------------------------------------

class GLMoleculeWidget(_QOpenGLWidget):
    """
    Hardware-accelerated 3D molecular renderer.

    Public interface is intentionally a subset of :class:`MolViewer3D` so
    the two can be swapped transparently inside a ``QStackedWidget``.

    Attributes
    ----------
    gl_available : bool
        ``True`` only after ``initializeGL`` succeeds with a usable context.
        The factory checks this to decide whether to keep using the GL path.
    molecule : Molecule | None
        Currently displayed molecule (domain model).
    rot_x, rot_y : float
        Euler-ish rotation angles (degrees).
    pan_x, pan_y : float
        Screen-space translation (pixels).
    zoom : float
        Zoom factor (pixels per Angstrom, roughly).
    """
    
    # Signals (matching MolViewer3D)
    atom_hovered = Signal(int)
    atom_clicked = Signal(int)
    selection_changed = Signal(object)
    delete_requested = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)

        # --- Public state ---
        self.molecule = None
        self.gl_available = None
        self.selected_atoms = set()
        self._selected_atoms_ordered = []

        # Camera state (matches MolViewer3D defaults)
        self.rot_x = 20.0
        self.rot_y = -30.0
        self.rot_z = 0.0
        self.pan_x = 0.0
        self.pan_y = 0.0
        self.zoom = 40.0

        # Rendering settings (subset matching MolViewer3D)
        self._show_hydrogens = False
        self.render_mode = 'ball_and_stick'
        self._sphere_scale = 0.6
        self._stick_scale = 1.0
        self._line_scale = 1.0
        self.bg_color = QColor(COLORS['viewer_bg'])
        self.custom_atom_colors = {}
        
        # Docking pose view settings
        self.show_ligands_in_cartoon = True
        self.visible_sidechains = set()
        self.interaction_lines = []
        self.custom_atom_modes = {}
        self.labels = {}
        self.labeled_residues = {}
        self.residue_label_settings = {}
        self.label_font_size = 9
        self.label_color = QColor(255, 255, 255, 230)

        # Atom data — packed numpy arrays for fast rendering
        self._positions = None   # (N, 3) float32
        self._colors = None      # (N, 3) float32  (0-1 range)
        self._radii = None       # (N,)   float32
        self._symbols = None     # list[str]  (length N)

        # Bond data
        self._bond_starts = None   # (M, 3) float32
        self._bond_ends = None     # (M, 3) float32
        self._bond_start_colors = None  # (M, 3) float32
        self._bond_end_colors = None    # (M, 3) float32

        # Mouse interaction state
        self._last_mouse_pos = None
        self._mouse_button = None

        # GL resources
        self._shader_mesh = None
        self._shader_line = None
        self._shader_sphere = None
        self._vbo_atoms = None
        self._vbo_mesh = None
        self._vbo_lines = None
        self._num_lines = 0
        self._ligand_start = 0
        self._ligand_bond_start = 0
        self._num_mesh_vertices = 0
        self._vao_atoms = None
        self._vao_mesh = None
        self._vao_lines = None

        # GL version info (populated in initializeGL)
        self._gl_version_string = ''
        self._gl_renderer_string = ''

    @property
    def show_hydrogens(self):
        return getattr(self, '_show_hydrogens', True)

    @show_hydrogens.setter
    def show_hydrogens(self, value):
        if getattr(self, '_show_hydrogens', True) != value:
            self._show_hydrogens = value
            if self.molecule:
                self.set_molecule(self.molecule)

    @property
    def sphere_scale(self):
        return self._sphere_scale

    @sphere_scale.setter
    def sphere_scale(self, value):
        if self._sphere_scale != value:
            self._sphere_scale = value
            if self.gl_available:
                self._update_gl_buffers()
            self.update()

    @property
    def stick_scale(self):
        return self._stick_scale

    @stick_scale.setter
    def stick_scale(self, value):
        if self._stick_scale != value:
            self._stick_scale = value
            self.update()

    @property
    def line_scale(self):
        return self._line_scale

    @line_scale.setter
    def line_scale(self, value):
        if self._line_scale != value:
            self._line_scale = value
            self.update()

        self.setMinimumSize(400, 400)
        self.setMouseTracking(True)

        # Request OpenGL 3.3 Core Profile (portable across all platforms)
        try:
            from PySide6.QtGui import QSurfaceFormat
            fmt = QSurfaceFormat()
            fmt.setVersion(3, 3)
            fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
            fmt.setDepthBufferSize(24)
            fmt.setSamples(4)
            self.setFormat(fmt)
        except Exception:
            try:
                from PyQt6.QtGui import QSurfaceFormat
                fmt = QSurfaceFormat()
                fmt.setVersion(3, 3)
                fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
                fmt.setDepthBufferSize(24)
                fmt.setSamples(4)
                self.setFormat(fmt)
            except Exception:
                pass  # Format request is best-effort

    # ------------------------------------------------------------------
    #  OpenGL lifecycle
    # ------------------------------------------------------------------

    def initializeGL(self):
        """Called once when the GL context is first created."""
        t_init = time.time()
        try:
            ctx = self.context()
            if ctx is None or not ctx.isValid():
                raise RuntimeError('No valid OpenGL context')

            gl = ctx.functions()
            if gl is None:
                raise RuntimeError('Could not obtain GL functions')
            gl.initializeOpenGLFunctions()

            # Query version string
            version_raw = gl.glGetString(0x1F02)  # GL_VERSION
            renderer_raw = gl.glGetString(0x1F01)  # GL_RENDERER

            self._gl_version_string = str(version_raw) if version_raw else 'unknown'
            self._gl_renderer_string = str(renderer_raw) if renderer_raw else 'unknown'

            # Parse major.minor from version string  (e.g. "4.1 INTEL-..." )
            major, minor = self._parse_gl_version(self._gl_version_string)
            if major < 3 or (major == 3 and minor < 3):
                raise RuntimeError(
                    f'OpenGL {major}.{minor} detected — 3.3+ required. '
                    f'Renderer: {self._gl_renderer_string}'
                )

            print(f'[GL] Context ready — OpenGL {major}.{minor}  '
                  f'Renderer: {self._gl_renderer_string}')

            # Basic GL state setup
            bg = self.bg_color
            gl.glClearColor(bg.redF(), bg.greenF(), bg.blueF(), 1.0)
            gl.glEnable(0x0B71)   # GL_DEPTH_TEST
            gl.glDepthFunc(0x0201)  # GL_LESS
            gl.glEnable(0x0BE2)   # GL_BLEND
            gl.glBlendFunc(0x0302, 0x0303)  # GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA
            gl.glEnable(0x8861)   # GL_MULTISAMPLE

            self._init_shaders()
            
            from src.shared.qt_compat import QOpenGLVertexArrayObject
            self._vao_atoms = QOpenGLVertexArrayObject()
            self._vao_atoms.create()
            
            self._vao_mesh = QOpenGLVertexArrayObject()
            self._vao_mesh.create()
            
            self._vao_lines = QOpenGLVertexArrayObject()
            self._vao_lines.create()

            self.gl_available = True
            print(f"[GL] Context ready — OpenGL {major}.{minor} Renderer: {self._gl_renderer_string}")
            print(f"[GL] Initialization complete in {time.time()-t_init:.3f}s")

        except Exception as exc:
            print(f'[GL] OpenGL initialisation failed: {exc}')
            traceback.print_exc()
            self.gl_available = False
        self.selected_atoms = set()

    def resizeGL(self, w, h):
        """Handle widget resize — update the GL viewport."""
        if not self.gl_available:
            return
        try:
            gl = self.context().functions()
            gl.glViewport(0, 0, w, h)
        except Exception:
            pass

    # ------------------------------------------------------------------
    #  Rendering
    # ------------------------------------------------------------------

    def paintGL(self):
        """Main render loop using VBOs and Shaders."""
        t_frame = time.time()
        if not self.gl_available:
            # Fallback to software (very slow for proteins)
            painter = QPainter(self)
            try:
                painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
                self._render_with_painter(painter, self.width(), self.height())
            finally:
                painter.end()
            return

        try:
            gl = self.context().functions()
            bg = self.bg_color
            gl.glClearColor(bg.redF(), bg.greenF(), bg.blueF(), 1.0)
            gl.glClear(0x00004100)  # COLOR | DEPTH
            
            if self.molecule is None:
                return

            # Set camera matrices
            proj = self._get_projection_matrix()
            view = self._get_view_matrix()
            
            # Use high-quality rendering states
            gl.glEnable(0x0B71) # GL_DEPTH_TEST
            gl.glEnable(0x8861) # GL_MULTISAMPLE
            
            # 1. Draw Mesh (Protein Cartoon)
            if self._vbo_mesh and getattr(self, '_num_mesh_vertices', 0) > 0:
                if self._vao_mesh: self._vao_mesh.bind()
                self._draw_mesh(gl, proj, view)
                if self._vao_mesh: self._vao_mesh.release()
                
            
            is_mesh_active = self._vbo_mesh and getattr(self, '_num_mesh_vertices', 0) > 0

            # 2. Draw Atoms (Sphere Impostors)
            if self._vbo_atoms and len(self._positions) > 0:
                if self._vao_atoms: self._vao_atoms.bind()
                
                self._draw_atoms(gl, proj, view)
                    
                if self._vao_atoms: self._vao_atoms.release()

            # 3. Draw Bonds (Sticks/Cylinders)
            if self._vbo_lines and self._num_lines > 0:
                self._shader_line.bind()
                self._shader_line.setUniformValue('projection', proj)
                self._shader_line.setUniformValue('view', view)
                
                # PySide6 lacks setUniformValue(str, float) overload, use native GL function
                loc = self._shader_line.uniformLocation('stick_scale')
                if loc != -1:
                    gl.glUniform1f(loc, float(self.stick_scale))
                if self._vao_lines: self._vao_lines.bind()
                self._vbo_lines.bind()
                self._shader_line.enableAttributeArray(0)
                self._shader_line.setAttributeBuffer(0, 0x1406, 0, 3, 56)
                self._shader_line.enableAttributeArray(1)
                self._shader_line.setAttributeBuffer(1, 0x1406, 12, 3, 56)
                self._shader_line.enableAttributeArray(2)
                self._shader_line.setAttributeBuffer(2, 0x1406, 24, 3, 56)
                self._shader_line.enableAttributeArray(3)
                self._shader_line.setAttributeBuffer(3, 0x1406, 36, 3, 56)
                self._shader_line.enableAttributeArray(4)
                self._shader_line.setAttributeBuffer(4, 0x1406, 48, 2, 56)
                
                gl.glDrawArrays(0x0004, 0, self._num_lines)
                    
                self._shader_line.disableAttributeArray(0)
                self._shader_line.disableAttributeArray(1)
                self._shader_line.disableAttributeArray(2)
                self._shader_line.disableAttributeArray(3)
                self._shader_line.disableAttributeArray(4)
                self._vbo_lines.release()
                if self._vao_lines: self._vao_lines.release()
                self._shader_line.release()
            
            # Log only if frame is very slow
            dt = time.time() - t_frame
            if dt > 0.1: # < 10 FPS
                print(f"[Performance] Slow frame: {dt:.3f}s")

        except Exception as e:
            print(f"[GL] Render error: {e}")
            traceback.print_exc()
            
        # Draw 2D overlays (Selection rings, interactions, labels) using QPainter over the GL surface
        if (self.selected_atoms or self.interaction_lines or self.labels or self.labeled_residues) and self._positions is not None:
            painter = QPainter(self)
            try:
                painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
                
                w = self.width()
                h = self.height()
                mvp = proj * view
                
                old_to_new = getattr(self.molecule, 'properties', {}).get('_old_to_new_idx', {})

                # Helper to project 3D point to 2D screen coordinate
                def project_point(pos_array):
                    vec4 = QVector4D(float(pos_array[0]), float(pos_array[1]), float(pos_array[2]), 1.0)
                    clip = mvp.map(vec4)
                    
                    if abs(clip.w()) > 1e-6:
                        ndc_x = clip.x() / clip.w()
                        ndc_y = clip.y() / clip.w()
                        ndc_z = clip.z() / clip.w()
                    else:
                        ndc_x = ndc_y = ndc_z = 0.0
                        
                    sx = (ndc_x + 1.0) * 0.5 * w
                    sy = (1.0 - ndc_y) * 0.5 * h
                    
                    vec3 = QVector3D(float(pos_array[0]), float(pos_array[1]), float(pos_array[2]))
                    view_pos = view.map(vec3)
                    return sx, sy, view_pos.z()

                # Draw Interaction Lines
                if self.interaction_lines:
                    for idx1, idx2, type_str, color_hex in self.interaction_lines:
                        n1 = old_to_new.get(idx1, idx1)
                        n2 = old_to_new.get(idx2, idx2)
                        if n1 >= len(self._positions) or n2 >= len(self._positions):
                            continue
                        
                        x1, y1, z1 = project_point(self._positions[n1])
                        x2, y2, z2 = project_point(self._positions[n2])
                        
                        # Simple culling
                        if z1 > 0 or z2 > 0: continue # Behind camera
                        
                        color = QColor(color_hex)
                        pen = QPen(color, 1.5, Qt.PenStyle.DashLine)
                        painter.setPen(pen)
                        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

                # Draw Selection Rings
                if self.selected_atoms:
                    for atom_idx in self.selected_atoms:
                        new_idx = old_to_new.get(atom_idx, atom_idx)
                        if new_idx >= len(self._positions):
                            continue
                        
                        sx, sy, z_dist = project_point(self._positions[new_idx])
                        if z_dist > 0: continue
                        
                        if z_dist < -0.1:
                            scale = (h / 2.0) / (math.tan(math.radians(22.5)) * -z_dist)
                        else:
                            scale = 1.0
                            
                        radius = self._radii[new_idx] * scale * self.sphere_scale
                        self._draw_selection_ring(painter, sx, sy, radius)

                # Draw Measurements based on ordered selection
                if hasattr(self, '_selected_atoms_ordered') and len(self._selected_atoms_ordered) >= 2:
                    font = painter.font()
                    font.setBold(True)
                    painter.setFont(font)
                    painter.setPen(QPen(Qt.GlobalColor.yellow, 2, Qt.PenStyle.DashLine))

                    if len(self._selected_atoms_ordered) == 2:
                        idx1, idx2 = self._selected_atoms_ordered
                        n1, n2 = old_to_new.get(idx1, idx1), old_to_new.get(idx2, idx2)
                        
                        if n1 < len(self._positions) and n2 < len(self._positions):
                            x1, y1, z1 = project_point(self._positions[n1])
                            x2, y2, z2 = project_point(self._positions[n2])
                            
                            if z1 <= 0 and z2 <= 0:
                                painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
                                
                                # Compute actual 3D distance using buffer positions
                                p1 = self._positions[n1]
                                p2 = self._positions[n2]
                                dx, dy, dz = p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]
                                dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                                
                                painter.drawText(int((x1+x2)/2 + 5), int((y1+y2)/2 - 5), f"{dist:.2f} Å")
                                
                    elif len(self._selected_atoms_ordered) == 3:
                        idx1, idx2, idx3 = self._selected_atoms_ordered
                        n1, n2, n3 = old_to_new.get(idx1, idx1), old_to_new.get(idx2, idx2), old_to_new.get(idx3, idx3)
                        
                        if n1 < len(self._positions) and n2 < len(self._positions) and n3 < len(self._positions):
                            x1, y1, z1 = project_point(self._positions[n1])
                            x2, y2, z2 = project_point(self._positions[n2])
                            x3, y3, z3 = project_point(self._positions[n3])
                            
                            if z1 <= 0 and z2 <= 0 and z3 <= 0:
                                painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
                                painter.drawLine(QPointF(x2, y2), QPointF(x3, y3))
                                
                                # Compute actual 3D angle using buffer positions
                                p1, p2, p3 = self._positions[n1], self._positions[n2], self._positions[n3]
                                v1 = np.array([p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]])
                                v2 = np.array([p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]])
                                norm1, norm2 = np.linalg.norm(v1), np.linalg.norm(v2)
                                if norm1 > 1e-6 and norm2 > 1e-6:
                                    cos_a = np.dot(v1, v2) / (norm1 * norm2)
                                    angle = math.degrees(math.acos(max(-1, min(1, cos_a))))
                                    painter.drawText(int(x2 + 10), int(y2 + 10), f"{angle:.1f}°")
                
                # Draw Labels
                if self.labels:
                    from src.features.visualization_3d.services.atom_rendering import draw_label
                    for atom_idx, text in self.labels.items():
                        new_idx = old_to_new.get(atom_idx, atom_idx)
                        if new_idx >= len(self._positions):
                            continue
                        
                        sx, sy, z_dist = project_point(self._positions[new_idx])
                        if z_dist > 0: continue
                        
                        if z_dist < -0.1:
                            scale = (h / 2.0) / (math.tan(math.radians(22.5)) * -z_dist)
                        else:
                            scale = 1.0
                            
                        radius = self._radii[new_idx] * scale * self.sphere_scale
                        draw_label(painter, text, sx, sy, radius, self.label_font_size, 1.0, self.label_color)

                # Draw Residue Labels
                if self.labeled_residues and self.residue_label_settings.get('show_labels', True):
                    from src.features.visualization_3d.services.atom_rendering import draw_residue_label
                    for atom in self.molecule.atoms:
                        rs = getattr(atom, 'res_seq', None)
                        if rs is not None and rs in self.labeled_residues:
                            if hasattr(atom, 'pdb_name') and atom.pdb_name.strip() == 'CA':
                                new_idx = old_to_new.get(atom.index, atom.index)
                                if new_idx >= len(self._positions):
                                    continue
                                
                                sx, sy, z_dist = project_point(self._positions[new_idx])
                                if z_dist > 0: continue
                                
                                if z_dist < -0.1:
                                    scale = (h / 2.0) / (math.tan(math.radians(22.5)) * -z_dist)
                                else:
                                    scale = 1.0
                                    
                                radius = self._radii[new_idx] * scale * self.sphere_scale
                                res_name = getattr(atom, 'res_name', 'UNK')
                                lbl_color = self.labeled_residues[rs]
                                draw_residue_label(painter, f"{res_name}{rs}", sx, sy, lbl_color, radius, self.label_font_size, 1.0, self.residue_label_settings)

            finally:
                painter.end()

    # ------------------------------------------------------------------
    #  QPainter-based rendering (used on top of the GL surface)
    # ------------------------------------------------------------------

    def _render_with_painter(self, painter: QPainter, w: int, h: int):
        """Software render pass using QPainter.

        This mirrors the essential rendering logic from ``PainterRenderer``
        but is deliberately simplified — no protein modes, no measurements,
        no atom picking overlays.  The goal is fast sphere + stick rendering
        for large molecules.
        """
        # Background (only if GL didn't already clear)
        if not self.gl_available:
            painter.fillRect(0, 0, w, h, self.bg_color)

        if self.molecule is None or self._positions is None:
            self._draw_placeholder(painter, w, h)
            return

        n = len(self._positions)
        if n == 0:
            self._draw_placeholder(painter, w, h)
            return

        # ── Project atoms ─────────────────────────────────────────────
        cx = w / 2.0 + self.pan_x
        cy = h / 2.0 + self.pan_y

        cos_x = math.cos(math.radians(self.rot_x))
        sin_x = math.sin(math.radians(self.rot_x))
        cos_y = math.cos(math.radians(self.rot_y))
        sin_y = math.sin(math.radians(self.rot_y))
        cos_z = math.cos(math.radians(self.rot_z))
        sin_z = math.sin(math.radians(self.rot_z))

        # Vectorised rotation
        pos = self._positions  # (N, 3)
        x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]

        x1 = x * cos_y + z * sin_y
        z1 = -x * sin_y + z * cos_y
        y1 = y * cos_x - z1 * sin_x
        z2 = y * sin_x + z1 * cos_x

        x2 = x1 * cos_z - y1 * sin_z
        y2 = x1 * sin_z + y1 * cos_z

        sx = cx + x2 * self.zoom
        sy = cy - y2 * self.zoom
        sz = z2

        # Display radius scaled
        display_r = self._radii * self.zoom * self.sphere_scale
        depth_factor = np.clip(1.0 + sz * 0.02, 0.5, 1.5)
        display_r = display_r * depth_factor

        # Depth sort (back to front)
        order = np.argsort(-sz)

        # ── Draw bonds ────────────────────────────────────────────────
        if self._bond_starts is not None and len(self._bond_starts) > 0:
            self._draw_bonds_painter(
                painter, cx, cy,
                cos_x, sin_x, cos_y, sin_y, cos_z, sin_z
            )

        # ── Draw atoms ────────────────────────────────────────────────
        large_molecule = n > 300
        painter.setPen(Qt.PenStyle.NoPen)

        for idx in order:
            asx = float(sx[idx])
            asy = float(sy[idx])
            asz = float(sz[idx])
            ar = float(display_r[idx])
            cr, cg, cb = self._colors[idx]

            # Convert 0-1 float colour to 0-255 int
            ri = int(cr * 255)
            gi = int(cg * 255)
            bi = int(cb * 255)

            if large_molecule and ar < 1.5:
                # Skip very small atoms in large molecules for speed
                continue

            if large_molecule:
                # Simple filled circle for speed
                color = QColor(ri, gi, bi)
                painter.setBrush(QBrush(color))
                painter.drawEllipse(QRectF(asx - ar, asy - ar, ar * 2, ar * 2))
            else:
                # Smooth radial gradient sphere
                self._draw_sphere(painter, asx, asy, ar, ri, gi, bi)

        # Performance indicator
        if large_molecule:
            self._draw_performance_label(painter, n)

    def _draw_sphere(self, painter: QPainter, sx, sy, radius, r, g, b):
        """Draw a single atom as a gradient sphere (same look as MolViewer3D)."""
        highlight_x = sx - radius * 0.30
        highlight_y = sy - radius * 0.30

        gradient = QRadialGradient(
            QPointF(sx, sy), radius, QPointF(highlight_x, highlight_y)
        )

        hl_r = min(255, r + 160)
        hl_g = min(255, g + 160)
        hl_b = min(255, b + 160)

        sh_r = max(0, int(r * 0.20))
        sh_g = max(0, int(g * 0.20))
        sh_b = max(0, int(b * 0.20))

        mid_r = max(0, min(255, int(r * 0.85)))
        mid_g = max(0, min(255, int(g * 0.85)))
        mid_b = max(0, min(255, int(b * 0.85)))

        gradient.setColorAt(0.0, QColor(hl_r, hl_g, hl_b))
        gradient.setColorAt(0.25, QColor(r, g, b))
        gradient.setColorAt(0.7, QColor(mid_r, mid_g, mid_b))
        gradient.setColorAt(1.0, QColor(sh_r, sh_g, sh_b))

        painter.setBrush(QBrush(gradient))
        painter.drawEllipse(QRectF(sx - radius, sy - radius, radius * 2, radius * 2))

    def _draw_selection_ring(self, painter: QPainter, sx: float, sy: float, radius: float):
        pen = QPen(QColor(255, 255, 0))
        pen.setWidth(2)
        painter.setPen(pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawEllipse(QRectF(sx - radius - 3, sy - radius - 3, (radius + 3) * 2, (radius + 3) * 2))

    def _draw_bonds_painter(self, painter, cx, cy, cos_x, sin_x, cos_y, sin_y, cos_z=1.0, sin_z=0.0):
        """Draw bonds as simple coloured lines (split-coloured)."""
        starts = self._bond_starts   # (M, 3)
        ends = self._bond_ends       # (M, 3)

        if starts is None or ends is None or len(starts) == 0:
            return

        # Project start points
        sx1_t = starts[:, 0] * cos_y + starts[:, 2] * sin_y
        sz1_tmp = -starts[:, 0] * sin_y + starts[:, 2] * cos_y
        sy1_t = starts[:, 1] * cos_x - sz1_tmp * sin_x
        
        sx1 = sx1_t * cos_z - sy1_t * sin_z
        sy1 = sx1_t * sin_z + sy1_t * cos_z

        psx1 = cx + sx1 * self.zoom
        psy1 = cy - sy1 * self.zoom

        # Project end points
        sx2_t = ends[:, 0] * cos_y + ends[:, 2] * sin_y
        sz2_tmp = -ends[:, 0] * sin_y + ends[:, 2] * cos_y
        sy2_t = ends[:, 1] * cos_x - sz2_tmp * sin_x
        
        sx2 = sx2_t * cos_z - sy2_t * sin_z
        sy2 = sx2_t * sin_z + sy2_t * cos_z

        psx2 = cx + sx2 * self.zoom
        psy2 = cy - sy2 * self.zoom

        bond_width = max(1.0, self.zoom * 0.06 * self.stick_scale)

        m = len(starts)
        for i in range(m):
            x1f = float(psx1[i])
            y1f = float(psy1[i])
            x2f = float(psx2[i])
            y2f = float(psy2[i])

            # Midpoint for split colouring
            mx = (x1f + x2f) / 2.0
            my = (y1f + y2f) / 2.0

            c1 = self._bond_start_colors[i]
            c2 = self._bond_end_colors[i]

            pen1 = QPen(QColor(int(c1[0] * 255), int(c1[1] * 255), int(c1[2] * 255)),
                        bond_width)
            pen1.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(pen1)
            painter.drawLine(QPointF(x1f, y1f), QPointF(mx, my))

            pen2 = QPen(QColor(int(c2[0] * 255), int(c2[1] * 255), int(c2[2] * 255)),
                        bond_width)
            pen2.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(pen2)
            painter.drawLine(QPointF(mx, my), QPointF(x2f, y2f))

    def _draw_placeholder(self, painter: QPainter, w, h):
        """Draw a placeholder message when no molecule is loaded."""
        if not self.gl_available:
            painter.fillRect(0, 0, w, h, self.bg_color)
        painter.setPen(QColor(150, 150, 150))
        font = QFont('Segoe UI', 12)
        painter.setFont(font)

        label = 'No molecule loaded'
        if self.gl_available:
            label += '  (GL accelerated)'
        else:
            label += '  (software fallback)'

        painter.drawText(QRectF(0, 0, w, h), Qt.AlignmentFlag.AlignCenter, label)

    def _draw_performance_label(self, painter: QPainter, n_atoms: int):
        """Show a small GL / atom-count badge in the corner."""
        painter.setPen(QColor(120, 120, 180, 180))
        font = QFont('Segoe UI', 8)
        painter.setFont(font)
        tag = 'GL' if self.gl_available else 'SW'
        text = f'{tag} | {n_atoms} atoms'
        painter.drawText(8, self.height() - 8, text)

    # ------------------------------------------------------------------
    #  Molecule data
    # ------------------------------------------------------------------

    def set_atom_colors(self, colors_dict):
        """Set custom colors for specific atoms and rebuild GPU buffers."""
        if colors_dict is None:
            self.custom_atom_colors = {}
        else:
            self.custom_atom_colors.update(colors_dict)
            
        if self.molecule and self._colors is not None:
            old_to_new = getattr(self.molecule, 'properties', {}).get('_old_to_new_idx', {})
            for atom in self.molecule.atoms:
                idx = atom.index
                new_idx = old_to_new.get(idx, idx)
                if new_idx >= len(self._colors):
                    continue
                if idx in self.custom_atom_colors:
                    c = self.custom_atom_colors[idx]
                    if isinstance(c, tuple) and len(c) == 3:
                        if isinstance(c[0], int):
                            c = [x / 255.0 for x in c]
                    elif hasattr(c, 'redF'):
                        c = [c.redF(), c.greenF(), c.blueF()]
                    self._colors[new_idx] = c
                else:
                    self._colors[new_idx] = _element_color_float(atom.symbol)
            
            # Re-upload bond colors
            self._pack_bonds(self.molecule, self.molecule.atoms, self._colors) # The atoms array here is ignored for unpacking since we remap it inside _pack_bonds correctly now
            
            if self.gl_available:
                self._update_gl_buffers()
        self.update()

    def set_molecule(self, mol):
        """Load a molecule for rendering.

        Packs atom positions, colours, and radii into flat numpy arrays.
        """
        self.molecule = mol

        if mol is None or not hasattr(mol, 'atoms') or len(mol.atoms) == 0:
            self._positions = None
            self._colors = None
            self._radii = None
            self._symbols = None
            self._bond_starts = None
            self._bond_ends = None
            self._bond_start_colors = None
            self._bond_end_colors = None
            self.update()
            return

        # Split atoms into protein vs ligand to render ligand bonds/atoms over/during mesh rendering
        protein_atoms = []
        ligand_atoms = []
        for a in mol.atoms:
            if a.symbol == 'H' and not self.show_hydrogens:
                continue
            if getattr(a, 'is_hetatm', False):
                ligand_atoms.append(a)
            else:
                protein_atoms.append(a)
                
        # Remap atom indices for bonds
        old_to_new = {}
        ordered_atoms = protein_atoms + ligand_atoms
        for i, a in enumerate(ordered_atoms):
            old_to_new[a.index] = i
        mol.properties['_old_to_new_idx'] = old_to_new
        
        self._ligand_start = len(protein_atoms)
        n = len(ordered_atoms)

        positions = np.zeros((n, 3), dtype=np.float32)
        colors = np.zeros((n, 3), dtype=np.float32)
        radii = np.zeros(n, dtype=np.float32)
        symbols = []
        
        self._atom_is_backbone = np.zeros(n, dtype=bool)
        self._atom_is_sidechain = np.zeros(n, dtype=bool)
        self._atom_res_seqs = np.zeros(n, dtype=int)
        
        _AMINO_ACIDS = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
        _BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'OXT'}

        for i, atom in enumerate(ordered_atoms):
            if atom.has_coords:
                positions[i] = [atom.x, atom.y, atom.z if atom.z is not None else 0.0]
            c = _element_color_float(atom.symbol)
            colors[i] = c
            radii[i] = _display_radius(atom.symbol)
            symbols.append(atom.symbol)
            
            rs = getattr(atom, 'res_seq', -1)
            if rs is not None:
                self._atom_res_seqs[i] = rs
                
            if i < self._ligand_start:
                if atom.res_name in _AMINO_ACIDS and atom.pdb_name.strip() not in _BACKBONE_ATOMS:
                    self._atom_is_sidechain[i] = True
                elif atom.res_name in _AMINO_ACIDS:
                    self._atom_is_backbone[i] = True
            
        # Center to origin
        if n > 0:
            self._centroid = np.mean(positions, axis=0)
            positions -= self._centroid
        else:
            self._centroid = np.zeros(3, dtype=np.float32)

        self._positions = positions
        self._colors = colors
        self._radii = radii
        self._symbols = symbols

        # Pack bond data
        self._pack_bonds(mol, ordered_atoms, colors)

        # Auto-fit camera
        self._auto_fit()

        if self.gl_available:
            self._update_gl_buffers()

        self.update()

    def _init_shaders(self):
        """Initialize GLSL shader programs."""
        try:
            from src.shared.qt_compat import QOpenGLShaderProgram, QOpenGLShader
            
            self._shader_mesh = QOpenGLShaderProgram()
            self._shader_mesh.addShaderFromSourceCode(QOpenGLShader.ShaderTypeBit.Vertex, MESH_VERTEX_SHADER)
            self._shader_mesh.addShaderFromSourceCode(QOpenGLShader.ShaderTypeBit.Fragment, MESH_FRAGMENT_SHADER)
            self._shader_mesh.link()
            
            self._shader_sphere = QOpenGLShaderProgram()
            self._shader_sphere.addShaderFromSourceCode(QOpenGLShader.ShaderTypeBit.Vertex, SPHERE_VERTEX_SHADER)
            self._shader_sphere.addShaderFromSourceCode(QOpenGLShader.ShaderTypeBit.Fragment, SPHERE_FRAGMENT_SHADER)
            self._shader_sphere.link()
            
            self._shader_line = QOpenGLShaderProgram()
            self._shader_line.addShaderFromSourceCode(QOpenGLShader.ShaderTypeBit.Vertex, LINE_VERTEX_SHADER)
            self._shader_line.addShaderFromSourceCode(QOpenGLShader.ShaderTypeBit.Fragment, LINE_FRAGMENT_SHADER)
            self._shader_line.link()
        except Exception as e:
            print(f"[GL] Shader init error: {e}")
            self.gl_available = False
        self.selected_atoms = set()

    def _update_gl_buffers(self):
        """Upload molecule data to GPU buffers."""
        if not self.gl_available: return
        import time
        t_upd = time.time()
        try:
            self.makeCurrent()
            gl = self.context().functions()
            
            # Calculate Atom Visibility Masks
            active_radii = self._radii * self.sphere_scale
            is_mesh_active = self.molecule and self.molecule.properties.get('is_protein')
            show_all_sidechains = getattr(self, 'show_sidechains', False)
            sidechain_res_vis = getattr(self, 'sidechain_res_vis', {})
            visible_sidechains = getattr(self, 'visible_sidechains', set())
            
            if is_mesh_active and hasattr(self, '_atom_is_backbone'):
                # Mask out all backbone atoms
                active_radii[self._atom_is_backbone] = 0.0
                
                # Mask out sidechains if not shown
                if not show_all_sidechains:
                    for i in range(len(active_radii)):
                        if self._atom_is_sidechain[i]:
                            rs = self._atom_res_seqs[i]
                            if not sidechain_res_vis.get(rs, False) and rs not in visible_sidechains:
                                active_radii[i] = 0.0

            # 1. Atoms Buffer (Center, Color, Radius, Offset)
            if self._vbo_atoms: self._vbo_atoms.destroy()
            n = len(self._positions)
            if n > 0:
                import numpy as np
                data = np.zeros(n * 6 * 9, dtype=np.float32)
                for i in range(6):
                    data[i*9+0::54] = self._positions[:, 0]
                    data[i*9+1::54] = self._positions[:, 1]
                    data[i*9+2::54] = self._positions[:, 2]
                    data[i*9+3::54] = self._colors[:, 0]
                    data[i*9+4::54] = self._colors[:, 1]
                    data[i*9+5::54] = self._colors[:, 2]
                    data[i*9+6::54] = active_radii
                offsets = np.array([[-1, -1], [1, -1], [1, 1], [-1, -1], [1, 1], [-1, 1]], dtype=np.float32)
                for i in range(6):
                    data[i*9+7::54] = offsets[i, 0]
                    data[i*9+8::54] = offsets[i, 1]
                from src.shared.qt_compat import QOpenGLBuffer
                self._vbo_atoms = QOpenGLBuffer(QOpenGLBuffer.Type.VertexBuffer)
                self._vbo_atoms.create()
                self._vbo_atoms.bind()
                self._vbo_atoms.allocate(data.tobytes(), len(data.tobytes()))
            
            # 2. Protein Mesh Buffer
            if self.molecule and self.molecule.properties.get('is_protein'):
                from src.features.visualization_3d.services.cartoon_generator import generate_cartoon_mesh
                from src.shared.ui.theme import COLORS
                t_mesh = time.time()
                # Always use high quality mesh for GPU renderer (it can easily handle it)
                v, t, c = generate_cartoon_mesh(self.molecule, spline_steps=24, profile_detail=16, theme_colors=COLORS)
                if v is not None:
                    if self._vbo_mesh: self._vbo_mesh.destroy()
                    indices = t.flatten()
                    v_flat = v[indices] - self._centroid
                    c_flat = c[indices]
                    self._num_mesh_vertices = len(v_flat)
                    mdata = np.zeros(self._num_mesh_vertices * 9, dtype=np.float32)
                    mdata[0::9] = v_flat[:, 0]
                    mdata[1::9] = v_flat[:, 1]
                    mdata[2::9] = v_flat[:, 2]
                    mdata[3::9] = c_flat[:, 0]
                    mdata[4::9] = c_flat[:, 1]
                    mdata[5::9] = c_flat[:, 2]
                    # Vectorized normal computation (replaces slow per-triangle Python loop)
                    n_tris = self._num_mesh_vertices // 3
                    tri_verts = v_flat[:n_tris * 3].reshape(n_tris, 3, 3)  # (T, 3_verts, 3_xyz)
                    edge1 = tri_verts[:, 1, :] - tri_verts[:, 0, :]
                    edge2 = tri_verts[:, 2, :] - tri_verts[:, 0, :]
                    normals = np.cross(edge1, edge2)  # (T, 3)
                    mags = np.linalg.norm(normals, axis=1, keepdims=True)
                    mags[mags < 1e-6] = 1.0
                    normals /= mags
                    # Broadcast same normal to all 3 vertices of each triangle
                    normals_expanded = np.repeat(normals, 3, axis=0)  # (T*3, 3)
                    mdata[6::9] = normals_expanded[:, 0]
                    mdata[7::9] = normals_expanded[:, 1]
                    mdata[8::9] = normals_expanded[:, 2]
                    from src.shared.qt_compat import QOpenGLBuffer
                    self._vbo_mesh = QOpenGLBuffer(QOpenGLBuffer.Type.VertexBuffer)
                    self._vbo_mesh.create()
                    self._vbo_mesh.bind()
                    self._vbo_mesh.allocate(mdata.tobytes(), len(mdata.tobytes()))
                    self._ibo_mesh = None
                    print(f"[GL] Mesh buffer updated in {time.time()-t_mesh:.3f}s")
            else:
                self._vbo_mesh = None
                
            # Calculate Bond Visibility Mask
            if self._bond_starts is not None and len(self._bond_starts) > 0:
                N_bonds = len(self._bond_starts)
                self._num_lines = N_bonds * 6
                
                line_data = np.empty((N_bonds, 6, 14), dtype=np.float32)
                
                starts = self._bond_starts.copy()
                ends = self._bond_ends.copy()
                
                if is_mesh_active and hasattr(self, '_bond_is_backbone'):
                    # Create boolean mask
                    bond_mask = np.ones(N_bonds, dtype=bool)
                    # Hide backbone bonds
                    if len(self._bond_is_backbone) == N_bonds:
                        bond_mask[self._bond_is_backbone] = False
                    
                    if not show_all_sidechains:
                        if not sidechain_res_vis and not visible_sidechains:
                            bond_mask[:self._ligand_bond_start // 6] = False
                        else:
                            # Use vectorized masking for explicitly visible sidechains
                            rs_array = self._bond_res_seq
                            mask_to_hide = np.ones(N_bonds, dtype=bool)
                            
                            for rs in visible_sidechains:
                                mask_to_hide[rs_array == rs] = False
                            for rs, vis in sidechain_res_vis.items():
                                if vis:
                                    mask_to_hide[rs_array == rs] = False
                            
                            # Hide sidechain bonds whose residue is NOT visible
                            # (we only apply this to protein bonds, i.e., before _ligand_bond_start)
                            is_protein_bond = np.arange(N_bonds) < (self._ligand_bond_start // 6)
                            bond_mask[mask_to_hide & (rs_array != -1) & is_protein_bond] = False
                    
                    starts[~bond_mask] = 0.0
                    ends[~bond_mask] = 0.0
                
                colors1 = self._bond_start_colors
                colors2 = self._bond_end_colors
                
                line_data[:, :, 0:3] = starts[:, np.newaxis, :]
                line_data[:, :, 3:6] = ends[:, np.newaxis, :]
                line_data[:, :, 6:9] = colors1[:, np.newaxis, :]
                line_data[:, :, 9:12] = colors2[:, np.newaxis, :]
                
                corners = np.array([
                    [-1.0, 0.0], [ 1.0, 0.0], [ 1.0, 1.0],
                    [-1.0, 0.0], [ 1.0, 1.0], [-1.0, 1.0]
                ], dtype=np.float32)
                line_data[:, :, 12:14] = corners
                
                line_data = line_data.reshape(-1, 14)
                
                from src.shared.qt_compat import QOpenGLBuffer
                self._vbo_lines = QOpenGLBuffer(QOpenGLBuffer.Type.VertexBuffer)
                self._vbo_lines.create()
                self._vbo_lines.bind()
                self._vbo_lines.allocate(line_data.tobytes(), len(line_data.tobytes()))
            else:
                self._vbo_lines = None
                self._num_lines = 0

            print(f"[Performance] Total GL Buffer update took {time.time()-t_upd:.3f}s")

        except Exception as e:
            print(f"[GL] Buffer update error: {e}")

    def _get_projection_matrix(self):
        proj = QMatrix4x4()
        aspect = self.width() / self.height() if self.height() > 0 else 1.0
        proj.perspective(45.0, aspect, 0.1, 1000.0)
        return proj

    def _get_view_matrix(self):
        view = QMatrix4x4()
        view.translate(self.pan_x / self.zoom, -self.pan_y / self.zoom, -100.0 + self.zoom)
        view.rotate(self.rot_x, 1, 0, 0)
        view.rotate(self.rot_y, 0, 1, 0)
        view.rotate(self.rot_z, 0, 0, 1)
        return view

    def _draw_atoms(self, gl, proj, view, start_idx=0, count=None):
        self._shader_sphere.bind()
        self._shader_sphere.setUniformValue("projection", proj)
        self._shader_sphere.setUniformValue("view", view)
        
        self._vbo_atoms.bind()
        self._shader_sphere.enableAttributeArray(0)
        self._shader_sphere.setAttributeBuffer(0, 0x1406, 0, 3, 36)
        self._shader_sphere.enableAttributeArray(1)
        self._shader_sphere.setAttributeBuffer(1, 0x1406, 12, 3, 36)
        self._shader_sphere.enableAttributeArray(2)
        self._shader_sphere.setAttributeBuffer(2, 0x1406, 24, 1, 36)
        self._shader_sphere.enableAttributeArray(3)
        self._shader_sphere.setAttributeBuffer(3, 0x1406, 28, 2, 36)
        
        if count is None:
            count = len(self._positions) - start_idx
        
        gl.glDrawArrays(0x0004, start_idx * 6, count * 6)
        
    def _draw_mesh(self, gl, proj, view):
        self._shader_mesh.bind()
        self._shader_mesh.setUniformValue("projection", proj)
        self._shader_mesh.setUniformValue("view", view)
        from src.shared.qt_compat import QMatrix4x4, QVector3D
        self._shader_mesh.setUniformValue("model", QMatrix4x4())
        self._shader_mesh.setUniformValue("lightPos", QVector3D(50, 50, 100))
        self._shader_mesh.setUniformValue("viewPos", QVector3D(0, 0, 100))
        self._shader_mesh.setUniformValue("lightColor", QVector3D(1, 1, 1))
        
        self._vbo_mesh.bind()
        self._shader_mesh.enableAttributeArray(0)
        self._shader_mesh.setAttributeBuffer(0, 0x1406, 0, 3, 36)
        self._shader_mesh.enableAttributeArray(1)
        self._shader_mesh.setAttributeBuffer(1, 0x1406, 12, 3, 36)
        self._shader_mesh.enableAttributeArray(2)
        self._shader_mesh.setAttributeBuffer(2, 0x1406, 24, 3, 36)
        
        gl.glDrawArrays(0x0004, 0, self._num_mesh_vertices)

    def _pack_bonds(self, mol, atoms, atom_colors):
        """Pack bond endpoint data into flat numpy arrays."""
        bonds = getattr(mol, 'bonds', None)
        if not bonds:
            self._bond_starts = None
            self._bond_ends = None
            self._bond_start_colors = None
            self._bond_end_colors = None
            self._ligand_bond_start = 0
            return

        n_atoms = len(atoms)
        
        prot_starts, prot_ends, prot_start_colors, prot_end_colors = [], [], [], []
        lig_starts, lig_ends, lig_start_colors, lig_end_colors = [], [], [], []

        old_to_new = getattr(mol, 'properties', {}).get('_old_to_new_idx', {})
        
        # Precompute rings for bond offsets (Skip for large structures to avoid 25s freeze)
        rings = []
        if hasattr(mol, 'find_rings') and len(atoms) < 1500:
            rings = mol.find_rings()
        bond_to_ring = {}
        for ring in rings:
            for i in range(len(ring)):
                a, b = ring[i], ring[(i+1) % len(ring)]
                bond_to_ring[frozenset([a, b])] = ring
                
        self._bond_is_backbone = []
        self._bond_res_seq = []
        
        for bond in bonds:
            bi = bond.begin_atom_idx
            ei = bond.end_atom_idx

            nbi = old_to_new.get(bi)
            nei = old_to_new.get(ei)

            if nbi is None or nei is None:
                continue
            
            # Map original indices to ordered indices
            new_bi = old_to_new.get(bi, bi)
            new_ei = old_to_new.get(ei, ei)
            
            if new_bi >= n_atoms or new_ei >= n_atoms:
                continue
            
            # Fetch original atoms from the molecule for property checks
            a1 = mol.atoms[bi]
            a2 = mol.atoms[ei]
            
            if not a1.has_coords or not a2.has_coords:
                continue
            # Skip hydrogen bonds if hydrogens hidden
            if not self.show_hydrogens and (a1.symbol == 'H' or a2.symbol == 'H'):
                continue
            
            is_backbone_bond = self._atom_is_backbone[new_bi] and self._atom_is_backbone[new_ei]

            p1 = self._positions[new_bi]
            p2 = self._positions[new_ei]
            
            # Check if this bond belongs to the ligand
            is_ligand = getattr(a1, 'is_hetatm', False) or getattr(a2, 'is_hetatm', False)
            
            bond_lines = []
            if is_ligand and (bond.is_double or bond.is_aromatic or bond.order >= 1.5):
                # Calculate a 3D offset perpendicular to the bond
                d = p2 - p1
                length = np.linalg.norm(d)
                if length > 1e-4:
                    d = d / length
                    offset_dir = None
                    # Try to use ring center for offset direction
                    ring = bond_to_ring.get(frozenset([bi, ei]))
                    if ring:
                        center = np.zeros(3, dtype=np.float32)
                        valid_atoms = 0
                        for idx in ring:
                            r_new_idx = old_to_new.get(idx)
                            if r_new_idx is not None and r_new_idx < n_atoms:
                                center += self._positions[r_new_idx]
                                valid_atoms += 1
                        if valid_atoms > 0:
                            center /= valid_atoms
                            mid = (p1 + p2) / 2.0
                            to_center = center - mid
                            to_center -= np.dot(to_center, d) * d
                            norm = np.linalg.norm(to_center)
                            if norm > 1e-4:
                                offset_dir = to_center / norm
                    
                    if offset_dir is None:
                        # Fallback orthogonal vector
                        up = np.array([0.0, 1.0, 0.0], dtype=np.float32)
                        if abs(np.dot(d, up)) > 0.9:
                            up = np.array([1.0, 0.0, 0.0], dtype=np.float32)
                        offset_dir = np.cross(d, up)
                        offset_dir /= np.linalg.norm(offset_dir)
                    
                    offset = offset_dir * 0.15 * self.stick_scale
                    
                    if bond.is_aromatic or bond.order == 1.5:
                        bond_lines.append((p1, p2))
                        bond_lines.append((p1 + offset, p2 + offset))
                    else:
                        bond_lines.append((p1 + offset*0.7, p2 + offset*0.7))
                        bond_lines.append((p1 - offset*0.7, p2 - offset*0.7))
                else:
                    bond_lines.append((p1, p2))
            else:
                bond_lines.append((p1, p2))
            
            for start_p, end_p in bond_lines:
                if is_ligand:
                    lig_starts.append(start_p)
                    lig_ends.append(end_p)
                    lig_start_colors.append(atom_colors[new_bi])
                    lig_end_colors.append(atom_colors[new_ei])
                else:
                    prot_starts.append(start_p)
                    prot_ends.append(end_p)
                    prot_start_colors.append(atom_colors[new_bi])
                    prot_end_colors.append(atom_colors[new_ei])
                    self._bond_is_backbone.append(is_backbone_bond)
                    
                    rs = getattr(a1, 'res_seq', -1)
                    if rs == -1:
                        rs = getattr(a2, 'res_seq', -1)
                    self._bond_res_seq.append(rs)

        self._ligand_bond_start = len(prot_starts) * 6
        
        self._bond_is_backbone.extend([False] * len(lig_starts))
        self._bond_is_backbone = np.array(self._bond_is_backbone, dtype=bool)
        
        self._bond_res_seq.extend([-1] * len(lig_starts))
        self._bond_res_seq = np.array(self._bond_res_seq, dtype=int)

        starts = prot_starts + lig_starts
        ends = prot_ends + lig_ends
        start_colors = prot_start_colors + lig_start_colors
        end_colors = prot_end_colors + lig_end_colors

        if starts:
            self._bond_starts = np.array(starts, dtype=np.float32)
            self._bond_ends = np.array(ends, dtype=np.float32)
            self._bond_start_colors = np.array(start_colors, dtype=np.float32)
            self._bond_end_colors = np.array(end_colors, dtype=np.float32)
        else:
            self._bond_starts = None
            self._bond_ends = None
            self._bond_start_colors = None
            self._bond_end_colors = None

    # ------------------------------------------------------------------
    #  Camera
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    #  Mouse interaction
    # ------------------------------------------------------------------

    def mousePressEvent(self, event):
        self._last_mouse_pos = event.position()
        self._mouse_button = event.button()
        self._mouse_moved = False

    def mouseMoveEvent(self, event):
        if self._last_mouse_pos is None:
            return

        dx = event.position().x() - self._last_mouse_pos.x()
        dy = event.position().y() - self._last_mouse_pos.y()
        
        if abs(dx) > 2 or abs(dy) > 2:
            self._mouse_moved = True

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
        if not getattr(self, '_mouse_moved', False) and self._mouse_button == Qt.MouseButton.LeftButton:
            self._handle_click(event.position())
        self._last_mouse_pos = None
        self._mouse_button = None
        self._mouse_moved = False

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
        if self._positions is None or len(self._positions) == 0:
            return -1

        w = self.width()
        h = self.height()
        proj = self._get_projection_matrix()
        view = self._get_view_matrix()
        mvp = proj * view
        
        click_x = pos.x()
        click_y = pos.y()
        
        closest_idx = -1
        min_dist_sq = float('inf')
        closest_z = float('inf')
        
        old_to_new = getattr(self.molecule, 'properties', {}).get('_old_to_new_idx', {})
        new_to_old = {v: k for k, v in old_to_new.items()}
        
        for i in range(len(self._positions)):
            vec = QVector3D(float(self._positions[i][0]), float(self._positions[i][1]), float(self._positions[i][2]))
            clip = mvp.map(vec)
            
            # Behind camera
            view_pos = view.map(vec)
            if view_pos.z() > 0:
                continue
                
            sx = (clip.x() + 1.0) * 0.5 * w
            sy = (1.0 - clip.y()) * 0.5 * h
            
            dx = sx - click_x
            dy = sy - click_y
            dist_sq = dx*dx + dy*dy
            
            # Determine radius to see if it's a hit
            if view_pos.z() < -0.1:
                scale = (h / 2.0) / (math.tan(math.radians(22.5)) * -view_pos.z())
            else:
                scale = 1.0
            radius = self._radii[i] * scale * self.sphere_scale
            hit_radius_sq = (max(10, radius * 1.5)) ** 2
            
            if dist_sq <= hit_radius_sq:
                if view_pos.z() < closest_z:
                    closest_z = view_pos.z()
                    closest_idx = new_to_old.get(i, i)
                    min_dist_sq = dist_sq
                    
        return closest_idx

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

    # ------------------------------------------------------------------
    #  Utilities
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_gl_version(version_str: str) -> tuple:
        """Extract (major, minor) ints from a GL_VERSION string.

        Examples:
            '4.1 INTEL-...'        -> (4, 1)
            '3.3.0 NVIDIA 535...'  -> (3, 3)
            'OpenGL ES 3.2 ...'    -> (3, 2)
            'garbage'              -> (0, 0)
        """
        import re
        m = re.search(r'(\d+)\.(\d+)', version_str)
        if m:
            return int(m.group(1)), int(m.group(2))
        return (0, 0)

    def clear(self):
        """Remove the current molecule and clear the view."""
        self.set_molecule(None)
