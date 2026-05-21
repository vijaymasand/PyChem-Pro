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
    QRadialGradient, QWheelEvent, QVector3D, QMatrix4x4,
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

        # Camera state (matches MolViewer3D defaults)
        self.rot_x = 20.0
        self.rot_y = -30.0
        self.rot_z = 0.0
        self.pan_x = 0.0
        self.pan_y = 0.0
        self.zoom = 40.0

        # Rendering settings (subset matching MolViewer3D)
        self.show_hydrogens = True
        self.render_mode = 'ball_and_stick'
        self.sphere_scale = 0.6
        self.stick_scale = 1.0
        self.bg_color = QColor(COLORS['viewer_bg'])

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
                
                # If protein mesh is active, only draw ligand spheres!
                if is_mesh_active and self._ligand_start < len(self._positions):
                    num_ligand_atoms = len(self._positions) - self._ligand_start
                    self._draw_atoms(gl, proj, view, start_idx=self._ligand_start, count=num_ligand_atoms)
                elif not is_mesh_active:
                    self._draw_atoms(gl, proj, view)
                    
                if self._vao_atoms: self._vao_atoms.release()

            # 3. Draw Bonds (Lines)
            if self._vbo_lines and self._num_lines > 0:
                self._shader_line.bind()
                self._shader_line.setUniformValue('projection', proj)
                self._shader_line.setUniformValue('view', view)
                
                gl.glLineWidth(2.0)
                if self._vao_lines: self._vao_lines.bind()
                self._vbo_lines.bind()
                self._shader_line.enableAttributeArray(0)
                self._shader_line.setAttributeBuffer(0, 0x1406, 0, 3, 24)
                self._shader_line.enableAttributeArray(1)
                self._shader_line.setAttributeBuffer(1, 0x1406, 12, 3, 24)
                
                if is_mesh_active:
                    num_ligand_lines = self._num_lines - self._ligand_bond_start
                    if num_ligand_lines > 0:
                        gl.glDrawArrays(0x0001, self._ligand_bond_start, num_ligand_lines)
                else:
                    gl.glDrawArrays(0x0001, 0, self._num_lines)
                    
                self._shader_line.disableAttributeArray(0)
                self._shader_line.disableAttributeArray(1)
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

        for i, atom in enumerate(ordered_atoms):
            if atom.has_coords:
                positions[i] = [atom.x, atom.y, atom.z if atom.z is not None else 0.0]
            c = _element_color_float(atom.symbol)
            colors[i] = c
            radii[i] = _display_radius(atom.symbol)
            symbols.append(atom.symbol)
            
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
                    data[i*9+6::54] = self._radii * self.sphere_scale
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
                t_mesh = time.time()
                # Always use high quality mesh for GPU renderer (it can easily handle it)
                v, t, c = generate_cartoon_mesh(self.molecule, spline_steps=24, profile_detail=16)
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
                
            # 3. Create Line VBO
            if self._bond_starts is not None and len(self._bond_starts) > 0:
                self._num_lines = len(self._bond_starts) * 2
                line_data = np.empty((len(self._bond_starts), 2, 6), dtype=np.float32)
                line_data[:, 0, 0:3] = self._bond_starts
                line_data[:, 0, 3:6] = self._bond_start_colors
                line_data[:, 1, 0:3] = self._bond_ends
                line_data[:, 1, 3:6] = self._bond_end_colors
                line_data = line_data.reshape(-1, 6)
                
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
        for bond in bonds:
            bi = bond.begin_atom_idx
            ei = bond.end_atom_idx
            
            # Map original indices to ordered indices
            new_bi = old_to_new.get(bi, bi)
            new_ei = old_to_new.get(ei, ei)
            
            if new_bi >= n_atoms or new_ei >= n_atoms:
                continue
            a1 = atoms[new_bi]
            a2 = atoms[new_ei]
            if not a1.has_coords or not a2.has_coords:
                continue
            # Skip hydrogen bonds if hydrogens hidden
            if not self.show_hydrogens and (a1.symbol == 'H' or a2.symbol == 'H'):
                continue

            p1 = self._positions[new_bi]
            p2 = self._positions[new_ei]
            
            # Check if this bond belongs to the ligand
            is_ligand = getattr(a1, 'is_hetatm', False) or getattr(a2, 'is_hetatm', False)
            
            if is_ligand:
                lig_starts.append(p1)
                lig_ends.append(p2)
                lig_start_colors.append(atom_colors[new_bi])
                lig_end_colors.append(atom_colors[new_ei])
            else:
                prot_starts.append(p1)
                prot_ends.append(p2)
                prot_start_colors.append(atom_colors[new_bi])
                prot_end_colors.append(atom_colors[new_ei])

        self._ligand_bond_start = len(prot_starts) * 2

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

    # ------------------------------------------------------------------
    #  Mouse interaction
    # ------------------------------------------------------------------

    def mousePressEvent(self, event):
        self._last_mouse_pos = event.position()
        self._mouse_button = event.button()

    def mouseMoveEvent(self, event):
        if self._last_mouse_pos is None:
            return

        dx = event.position().x() - self._last_mouse_pos.x()
        dy = event.position().y() - self._last_mouse_pos.y()

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
        self._last_mouse_pos = None
        self._mouse_button = None

    def wheelEvent(self, event: QWheelEvent):
        delta = event.angleDelta().y()
        factor = 1.1 if delta > 0 else 0.9
        self.zoom *= factor
        self.zoom = max(5.0, min(200.0, self.zoom))
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
