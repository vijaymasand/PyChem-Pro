"""
OpenGL-accelerated 3D Molecular Viewer — Hardware-rendered path for large molecules.

Cross-platform notes
--------------------
- macOS: Apple deprecated OpenGL at 4.1 max. We request 3.3 Core Profile.
- Windows: Older Intel iGPUs may lack 3.3. We detect and report via gl_available.
- Linux: Mesa driver quality varies. All GL errors are caught and reported.
"""

import math
import sys
import traceback
import ctypes
import time
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

from src.features.visualization_3d.ui.gl_widget_helpers import _element_color_float, _display_radius
from src.features.visualization_3d.ui.gl_widget_data import GLDataMixin
from src.features.visualization_3d.ui.gl_widget_events import GLEventMixin

_QOpenGLWidget = None

try:
    from PySide6.QtOpenGLWidgets import QOpenGLWidget as _QOpenGLWidget
except ImportError:
    try:
        from PyQt6.QtOpenGLWidgets import QOpenGLWidget as _QOpenGLWidget
    except ImportError:
        pass

if _QOpenGLWidget is None:
    from src.shared.qt_compat import QWidget as _QOpenGLWidget  # type: ignore[assignment]


class GLMoleculeWidget(GLDataMixin, GLEventMixin, _QOpenGLWidget):
    """Hardware-accelerated 3D molecular renderer.

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
    Zoom factor (pixels per Angstrom, roughly)."""
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
        self._mouse_moved = False

        # Rubber-band selection state
        self._is_selecting = False
        self._sel_rect_origin = None
        self._sel_rect_end = None

        # Lasso selection state
        self._is_lasso = False
        self._lasso_path = []

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
            fmt.setSamples(8)  # Increased multi-sampling for smoother edges
            self.setFormat(fmt)
        except Exception:
            try:
                from PyQt6.QtGui import QSurfaceFormat
                fmt = QSurfaceFormat()
                fmt.setVersion(3, 3)
                fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
                fmt.setDepthBufferSize(24)
                fmt.setSamples(8)
                self.setFormat(fmt)
            except Exception:
                pass  # Format request is best-effort

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

            # print(f'[GL] Context ready — OpenGL {major}.{minor}  '
            #       f'Renderer: {self._gl_renderer_string}')

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
            # print(f"[GL] Context ready — OpenGL {major}.{minor} Renderer: {self._gl_renderer_string}")
            # print(f"[GL] Initialization complete in {time.time()-t_init:.3f}s")

        except Exception as exc:
            # print(f'[GL] OpenGL initialisation failed: {exc}')
            # traceback.print_exc()
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
            gl.glEnable(0x809E) # GL_SAMPLE_ALPHA_TO_COVERAGE
            
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
                # print(f"[Performance] Slow frame: {dt:.3f}s")
                pass

        except Exception as e:
            # print(f"[GL] Render error: {e}")
            # traceback.print_exc()
            pass
            
        # Draw 2D overlays (Selection rings, interactions, labels, dummy spheres) using QPainter over the GL surface
        has_dummy_spheres = hasattr(self.molecule, 'dummy_spheres') and bool(self.molecule.dummy_spheres)
        if (self.selected_atoms or self.interaction_lines or self.labels or self.labeled_residues or has_dummy_spheres) and (self._positions is not None or has_dummy_spheres):
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

                # Draw Dummy Spheres (COM, Centroid, Custom)
                if has_dummy_spheres:
                    self._draw_dummy_spheres_gl(painter, w, h, proj, view, mvp)

                # Draw Interaction Lines
                if self.interaction_lines and self._positions is not None:
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
                if self.selected_atoms and self._positions is not None:
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
                if hasattr(self, '_selected_atoms_ordered') and len(self._selected_atoms_ordered) >= 2 and self._positions is not None:
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
                if self.labels and self._positions is not None:
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
                if self.labeled_residues and self.residue_label_settings.get('show_labels', True) and self._positions is not None:
                    from src.features.visualization_3d.services.atom_rendering import draw_residue_label, get_formatted_residue_label
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
                                label_text = get_formatted_residue_label(res_name, rs, self.residue_label_settings)
                                draw_residue_label(painter, label_text, sx, sy, lbl_color, radius, self.label_font_size, 1.0, self.residue_label_settings)

            finally:
                painter.end()

    def _draw_dummy_spheres_gl(self, painter, w, h, proj, view, mvp):
        """Draw dummy spheres (COM, centroid, custom) over the OpenGL viewport using QPainter."""
        if not hasattr(self.molecule, 'dummy_spheres') or not self.molecule.dummy_spheres:
            return

        dummy_spheres = self.molecule.dummy_spheres
        centroid = getattr(self, '_centroid', None)
        if centroid is None:
            centroid = np.zeros(3, dtype=np.float32)

        spheres_with_depth = []
        for sphere in dummy_spheres:
            if not getattr(sphere, 'visible', True):
                continue

            pos = getattr(sphere, 'position', (0.0, 0.0, 0.0))
            x = float(pos[0]) - float(centroid[0])
            y = float(pos[1]) - float(centroid[1])
            z = float(pos[2]) - float(centroid[2])

            vec4 = QVector4D(x, y, z, 1.0)
            clip = mvp.map(vec4)

            if abs(clip.w()) > 1e-6:
                ndc_x = clip.x() / clip.w()
                ndc_y = clip.y() / clip.w()
                ndc_z = clip.z() / clip.w()
            else:
                ndc_x = ndc_y = ndc_z = 0.0

            sx = (ndc_x + 1.0) * 0.5 * w
            sy = (1.0 - ndc_y) * 0.5 * h

            vec3 = QVector3D(x, y, z)
            view_pos = view.map(vec3)
            z_dist = view_pos.z()

            if z_dist > 0:  # Behind camera
                continue

            spheres_with_depth.append((sphere, sx, sy, z_dist))

        sorted_spheres = sorted(spheres_with_depth, key=lambda s: (s[3], -getattr(s[0], 'radius', 0.5)))

        tan_half_fov = math.tan(math.radians(22.5))
        for sphere, sx, sy, z_dist in sorted_spheres:
            radius = getattr(sphere, 'radius', 0.5)
            alpha = getattr(sphere, 'alpha', 1.0)

            if z_dist < -0.001:
                scale = (h / 2.0) / (tan_half_fov * -z_dist)
            else:
                scale = 1.0

            label = getattr(sphere, 'label', '')
            if label in ['COM', 'Centroid']:
                display_r = radius * scale * getattr(self, 'sphere_scale', 1.0)
            else:
                display_r = radius * scale

            color_hex = getattr(sphere, 'color', '#ffff00')
            from src.shared.ui.theme import COLORS
            if label == 'COM':
                color_hex = COLORS.get('sphere_com', color_hex)
            elif label == 'Centroid':
                color_hex = COLORS.get('sphere_centroid', color_hex)

            color_q = QColor(color_hex)
            alpha_val = max(0.0, min(1.0, alpha))

            grad = QRadialGradient(QPointF(sx - display_r * 0.3, sy - display_r * 0.3), max(1.0, display_r * 1.3))
            highlight = QColor(255, 255, 255, int(alpha_val * 220))
            mid_color = QColor(color_q.red(), color_q.green(), color_q.blue(), int(alpha_val * 230))
            dark = QColor(max(0, color_q.red() - 70), max(0, color_q.green() - 70), max(0, color_q.blue() - 70), int(alpha_val * 255))

            grad.setColorAt(0.0, highlight)
            grad.setColorAt(0.35, mid_color)
            grad.setColorAt(1.0, dark)

            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(grad))
            r_val = max(2.0, display_r)
            painter.drawEllipse(QRectF(sx - r_val, sy - r_val, r_val * 2.0, r_val * 2.0))

            is_custom = label and label.lower() == 'custom'
            if label and alpha_val > 0.2 and not is_custom:
                painter.setPen(QColor(255, 255, 255, int(alpha_val * 255)))
                font = painter.font()
                font.setPointSize(8)
                font.setBold(False)
                painter.setFont(font)
                label_rect = QRectF(sx + r_val + 5, sy - 10, 150, 20)
                painter.drawText(label_rect, Qt.AlignmentFlag.AlignLeft, label)

        # Draw lasso polygon overlay (while Ctrl+dragging)
        if getattr(self, '_is_lasso', False) and getattr(self, '_lasso_path', []):
            from src.shared.qt_compat import QPainterPath
            lasso = self._lasso_path
            if len(lasso) >= 2:
                path = QPainterPath()
                path.moveTo(lasso[0])
                for pt in lasso[1:]:
                    path.lineTo(pt)
                path.closeSubpath()

                lasso_painter = QPainter(self)
                try:
                    lasso_painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
                    pen = QPen(QColor(0, 220, 255, 220), 1.8, Qt.PenStyle.DashLine)
                    pen.setDashPattern([6, 3])
                    lasso_painter.setPen(pen)
                    lasso_painter.setBrush(QBrush(QColor(0, 200, 255, 30)))
                    lasso_painter.drawPath(path)
                finally:
                    lasso_painter.end()

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

