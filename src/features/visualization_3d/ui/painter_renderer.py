"""
QPainter-based rendering engine for MolViewer3D.

Handles 3D-to-2D projection, depth sorting, atom sphere rendering,
bond drawing, labels, protein representations, measurements overlay,
and dummy sphere rendering.
"""

import math
import time
from src.shared.qt_compat import (
    Qt, QPointF, QRectF,
    QPainter, QColor, QPen, QBrush, QFont, QPainterPath,
    QPolygonF, QRadialGradient,
)
from src.shared.ui.theme import COLORS


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


class PainterRenderer:
    """
    Rendering helper for :class:`MolViewer3D`.

    Every public method receives the viewer (or its state) explicitly.
    Lightweight caches (gradient, spatial hash) are maintained for
    performance but are automatically invalidated as needed.
    """

    def __init__(self):
        # Gradient cache: (element_symbol, radius_bucket, is_hovered, use_ssao) -> gradient color stops
        # We cache the *color stops* (not the positioned QRadialGradient) because
        # the gradient centre changes every frame.  Rebuilding a QRadialGradient
        # from cached stops is cheap; computing the colour arithmetic is not.
        self._gradient_cache = {}

    # ─── Cache management ────────────────────────────────────────

    def invalidate_cache(self):
        """Clear all rendering caches (call when colour scheme changes)."""
        self._gradient_cache.clear()

    def _get_cached_gradient(self, element_symbol, radius, sx, sy, rgb,
                             is_hovered=False, use_ssao=False, sz=0.0,
                             alpha=1.0):
        """Return a positioned QRadialGradient, reusing cached colour stops.

        The cache key is ``(element_symbol, radius_bucket, is_hovered,
        use_ssao)`` so atoms of the same element at the same visual size
        share the same computed colour stops.  The gradient is then
        repositioned to *(sx, sy)* which is a trivial operation.
        """
        bucket = round(radius * 2)  # half-pixel buckets
        # Quantise depth for shading so nearby atoms share the same entry
        depth_bucket = round(sz * 4)
        key = (element_symbol, bucket, is_hovered, use_ssao, depth_bucket,
               int(alpha * 255))

        if key not in self._gradient_cache:
            r, g, b = rgb

            # Depth-based shading (mirrors atom_rendering.draw_atom_sphere)
            if use_ssao:
                depth_shade = max(0.2, min(1.0, 1.0 - sz * 0.05))
            else:
                depth_shade = max(0.5, min(1.0, 1.0 - sz * 0.025))
            r = int(r * depth_shade)
            g = int(g * depth_shade)
            b = int(b * depth_shade)

            if is_hovered:
                r = min(255, r + 50)
                g = min(255, g + 50)
                b = min(255, b + 50)

            r = max(0, min(255, r))
            g = max(0, min(255, g))
            b = max(0, min(255, b))

            alpha_int = int(alpha * 255)
            highlight = QColor(min(255, r + 160), min(255, g + 160),
                               min(255, b + 160), alpha_int)
            base = QColor(r, g, b, alpha_int)
            mid = QColor(max(0, min(255, int(r * 0.85))),
                         max(0, min(255, int(g * 0.85))),
                         max(0, min(255, int(b * 0.85))), alpha_int)
            shadow = QColor(max(0, int(r * 0.20)),
                            max(0, int(g * 0.20)),
                            max(0, int(b * 0.20)), alpha_int)

            # Also cache the flat colour for the simple-dot LOD path
            flat_color = QColor(r, g, b, alpha_int)

            # Rim lighting color
            rim_color = QColor(min(255, r + 30), min(255, g + 30), min(255, b + 30), alpha_int)

            self._gradient_cache[key] = (highlight, base, mid, shadow, rim_color,
                                         flat_color, radius > 4 and alpha > 0.5)

        stops = self._gradient_cache[key]
        highlight, base_c, mid_c, shadow_c, rim_c, _flat, do_specular = stops

        # Build a positioned gradient (cheap — no colour math)
        highlight_x = sx - radius * 0.35
        highlight_y = sy - radius * 0.35
        gradient = QRadialGradient(QPointF(sx, sy), radius,
                                   QPointF(highlight_x, highlight_y))
        gradient.setColorAt(0.0, highlight)
        gradient.setColorAt(0.08, highlight)
        gradient.setColorAt(0.35, base_c)
        gradient.setColorAt(0.75, mid_c)
        gradient.setColorAt(0.92, shadow_c)
        gradient.setColorAt(1.0, rim_c)

        return gradient, do_specular

    # ─── Spatial hash for bond lookup ────────────────────────────

    @staticmethod
    def _build_spatial_hash(projected, cell_size=50):
        """Build a grid mapping cell coordinates to projected-atom indices.

        Used by ``_draw_bonds`` to avoid a full scan of the projected list
        when locating atom endpoints.  Reduces bond lookup from O(N) dict
        build to an O(1)-amortised grid probe.
        """
        grid = {}
        proj_map = {}
        for entry in projected:
            atom_idx = entry[0]
            sx, sy = entry[1], entry[2]
            proj_map[atom_idx] = entry
            gx = int(sx // cell_size)
            gy = int(sy // cell_size)
            grid.setdefault((gx, gy), []).append(atom_idx)
        return grid, proj_map

    # ─── Projection ──────────────────────────────────────────────

    def _project_atoms(self, v, vp_width=None, vp_height=None):
        """Project 3D atom coordinates to 2D screen coordinates.

        *v* is the MolViewer3D instance whose camera / molecule state is
        used for the projection.
        """
        if not v.molecule:
            return []

        w = vp_width or v.width()
        h = vp_height or v.height()
        cx = w / 2 + v.pan_x
        cy = h / 2 + v.pan_y

        cos_x = math.cos(math.radians(v.rot_x))
        sin_x = math.sin(math.radians(v.rot_x))
        cos_y = math.cos(math.radians(v.rot_y))
        sin_y = math.sin(math.radians(v.rot_y))
        cos_z = math.cos(math.radians(getattr(v, 'rot_z', 0.0)))
        sin_z = math.sin(math.radians(getattr(v, 'rot_z', 0.0)))

        projected = []
        for atom in v.molecule.atoms:
            if not atom.has_coords:
                continue
            if atom.symbol == 'H' and not v.show_hydrogens:
                continue

            x, y, z = atom.x, atom.y, atom.z

            # Rotation around Y then X
            x1 = x * cos_y + z * sin_y
            z1 = -x * sin_y + z * cos_y
            y1 = y * cos_x - z1 * sin_x
            z2 = y * sin_x + z1 * cos_x
            
            x2 = x1 * cos_z - y1 * sin_z
            y2 = x1 * sin_z + y1 * cos_z

            sx = cx + x2 * v.zoom
            sy = cy - y2 * v.zoom
            sz = z2

            # Display radius with user scale
            base_r = DISPLAY_RADIUS.get(atom.symbol, 0.35)
            atom_render_mode = getattr(v, 'custom_atom_modes', {}).get(atom.index, v.render_mode)
            if atom_render_mode == 'spacefill':
                base_r = atom.element.vdw_radius * 0.5
            elif atom_render_mode == 'wireframe':
                base_r = 0.1
            display_r = base_r * v.zoom * v.sphere_scale

            # Depth-based size attenuation
            depth_factor = 1.0 + sz * 0.02
            display_r *= max(0.5, min(1.5, depth_factor))

            # Color from element or custom override
            custom_colors = getattr(v, 'custom_atom_colors', {})
            if atom.index in custom_colors:
                color = custom_colors[atom.index]
            else:
                color = _hex_to_rgb(atom.element.color)

            projected.append((atom.index, sx, sy, sz, display_r, color))

        return projected

    # ─── Main render entry-point ─────────────────────────────────

    def render(self, v, painter, width, height, is_export=False, export_scale=1.0):
        """Core rendering logic -- used by both paintEvent and export."""
        # Performance optimization: Adaptive antialiasing based on molecule size
        num_atoms = len(v.molecule.atoms) if v.molecule else 0
        if num_atoms > 8000:
            # Disable antialiasing for very large molecules to improve performance
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        elif num_atoms > 5000:
            # Use high-quality antialiasing for medium molecules
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
        elif num_atoms > 2000:
            # Use standard antialiasing for small-medium molecules
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
        else:
            # Use full antialiasing for small molecules
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
            # Note: HighQualityAntialiasing not available in PySide6
            
        painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)

        v._export_scale = export_scale

        # Background
        painter.fillRect(0, 0, width, height, v.bg_color)

        if not v.molecule or len(v.molecule.atoms) == 0:
            if not is_export:
                self._draw_placeholder(painter, width, height)
            return

        # Performance optimization: Limit gradient cache size and add mesh caching
        if len(self._gradient_cache) > 2048:  # Reduced cache size
            self._gradient_cache.clear()
        
        # Initialize mesh cache for performance
        if not hasattr(self, '_mesh_cache'):
            self._mesh_cache = {}  # molecule_id -> (vertices, triangles, colors, timestamp)

        # Project atoms
        projected = self._project_atoms(v, width, height)
        if not projected:
            return

        # Sort by depth (furthest first = painter's algorithm)
        sorted_atoms = sorted(projected, key=lambda x: x[3])

        # Off-screen culling margin (pixels)
        _cull_margin = 50

        # --- Protein modes ---
        if v.render_mode in ('cartoon', 'ribbon', 'backbone'):
            self._draw_protein(v, painter, projected, width, height)

            # Draw any custom-styled atoms over the protein backbone
            cam = getattr(v, 'custom_atom_modes', {})
            
            # AUTOMATIC HETEROATOM/SIDECHAIN LOGIC
            atoms_to_draw = set()
            for atom_idx in cam:
                atoms_to_draw.add(atom_idx)
            
            # Add HETATMs (ligands) automatically
            if getattr(v, 'show_ligands_in_cartoon', True):
                for atom in v.molecule.atoms:
                    if getattr(atom, 'is_hetatm', False):
                        # Filter out water by default
                        res_name = getattr(atom, 'res_name', '').upper()
                        if res_name not in ('HOH', 'WAT', 'SOL', 'DOD'):
                            atoms_to_draw.add(atom.index)
            
            # Add specific residues (sidechains)
            visible_sidechains = getattr(v, 'visible_sidechains', set())
            if visible_sidechains:
                for atom in v.molecule.atoms:
                    if atom.res_seq in visible_sidechains:
                        # Draw sidechain atoms (not backbone CA, C, N, O)
                        pdb_name = getattr(atom, 'pdb_name', '').strip()
                        if pdb_name not in ('CA', 'C', 'N', 'O'):
                            atoms_to_draw.add(atom.index)
            
            # Ensure any manually styled atoms (like docking ligands) are visible
            if hasattr(v, 'custom_atom_modes'):
                atoms_to_draw.update(v.custom_atom_modes.keys())
            
            if atoms_to_draw:
                self._draw_bonds(v, painter, projected, custom_only=True, 
                                 width=width, height=height, force_indices=atoms_to_draw)
                
                # DRAW INTERACTION LINES
                if getattr(v, 'interaction_lines', []):
                    self._draw_interaction_lines(v, painter, projected, width, height)

                for atom_idx, sx, sy, sz, radius, color in sorted_atoms:
                    if atom_idx in atoms_to_draw:
                        # Off-screen culling
                        if (sx < -_cull_margin or sx > width + _cull_margin or
                                sy < -_cull_margin or sy > height + _cull_margin):
                            continue
                        self._draw_atom_sphere(v, painter, atom_idx, sx, sy, sz, radius, color)
            

        else:
            # Optimize rendering for VERY large molecules only.
            # Below 8000 atoms we still use gradient spheres so PDB proteins
            # look like real 3D spheres. The gradient cache + off-screen
            # culling + LOD keep the frame rate acceptable.
            num_atoms = len(sorted_atoms)
            use_simple_rendering = num_atoms > 4000  # Reduced threshold for better performance

            if use_simple_rendering and v.render_mode == 'ball_and_stick':
                self._draw_large_molecule_fast(v, painter, projected, sorted_atoms,
                                               width, height)
            else:
                # Draw bonds first (behind atoms)
                self._draw_bonds(v, painter, projected, width=width, height=height)

                # Draw atoms with smooth gradient spheres
                for atom_idx, sx, sy, sz, radius, color in sorted_atoms:
                    # Off-screen culling
                    if (sx < -_cull_margin or sx > width + _cull_margin or
                            sy < -_cull_margin or sy > height + _cull_margin):
                        continue

                    if use_simple_rendering:
                        self._draw_atom_simple(painter, sx, sy, radius, color)
                    elif radius < 2:
                        # LOD: tiny atom — draw as simple filled circle (no gradient)
                        self._draw_atom_dot(painter, sx, sy, radius, color)
                    else:
                        self._draw_atom_sphere(v, painter, atom_idx, sx, sy, sz, radius, color)

        # DRAW CUSTOM LABELS (e.g. lipophilic contributions)
        if getattr(v, 'labels', {}):
            self._draw_labels(v, painter, projected)

        # Draw selection highlights
        if v.selected_atoms:
            for atom_idx, sx, sy, sz, radius, color in sorted_atoms:
                if atom_idx in v.selected_atoms:
                    self._draw_selection_ring(painter, sx, sy, radius)

        # Draw rubber-band selection rectangle (while Shift+dragging)
        if v._is_selecting and v._sel_rect_origin and v._sel_rect_end:
            self._draw_rubber_band(v, painter)

        # Draw SASA point-cloud surface
        if getattr(v, 'show_sasa_surface', False):
            self._draw_sasa_surface(v, painter, width, height)

        # Draw labels if enabled
        if v.show_labels:
            for atom_idx, sx, sy, sz, radius, color in sorted_atoms:
                self._draw_label(v, painter, atom_idx, sx, sy, radius)

        # Draw specific residue labels dynamically
        if hasattr(v, 'labeled_residues') and v.labeled_residues:
            settings = getattr(v, 'residue_label_settings', {})
            if settings.get('show_labels', True):
                for atom_idx, sx, sy, sz, radius, color in sorted_atoms:
                    atom = v.molecule.atoms[atom_idx]
                    rs = getattr(atom, 'res_seq', None)
                    if rs is not None and rs in v.labeled_residues:
                        if hasattr(atom, 'pdb_name') and atom.pdb_name.strip() == 'CA':
                            res_name = getattr(atom, 'res_name', 'UNK')
                            lbl_color = v.labeled_residues[rs]
                            self._draw_residue_label(v, painter, f"{res_name}{rs}", sx, sy, lbl_color, radius, settings)

        # Draw dummy spheres (COM, centroid, custom)
        self._draw_dummy_spheres(v, painter, width, height)

        # Draw measurements
        if v._measure_atoms or v._measurements:
            self._draw_measurements(v, painter, projected)

        # Draw info overlay (not on exports)
        if not is_export:
            self._draw_overlay(v, painter)

            # Show performance indicator for large molecules
            if len(projected) > 500:
                self._draw_performance_indicator(v, painter, len(projected))

    # ─── Atom drawing ────────────────────────────────────────────

    def _draw_atom_dot(self, painter, sx, sy, radius, rgb):
        """LOD fast-path: draw a tiny atom as a plain filled circle (no gradient)."""
        if isinstance(rgb, tuple):
            r, g, b = rgb
            color = QColor(r, g, b)
        else:
            color = rgb
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QBrush(color))
        dot_r = max(1.5, radius)
        painter.drawEllipse(QPointF(sx, sy), dot_r, dot_r)

    def _draw_atom_sphere(self, v, painter, atom_idx, sx, sy, sz, radius, rgb, alpha=1.0):
        """Draw an atom sphere with performance optimizations for large proteins.

        Performance optimizations:
        - Simplified gradients for large molecules
        - Reduced specular calculations
        - Early termination for distant atoms
        """
        is_hovered = (atom_idx == v._hovered_atom)
        use_ssao = getattr(v, 'use_ssao', False)

        # Resolve element symbol for cache key
        if atom_idx >= 0 and v.molecule and atom_idx < len(v.molecule.atoms):
            atom = v.molecule.atoms[atom_idx]
            elem_sym = atom.symbol
            # Check if this atom is part of a custom render set (ligand/sidechain)
            is_custom = atom_idx in getattr(v, 'custom_atom_modes', {})
        else:
            elem_sym = '_dummy'
            is_custom = False

        # Performance optimization: Skip complex rendering for very large molecules
        # but ALWAYS draw ligands/sidechains (custom modes) in high quality.
        num_atoms = len(v.molecule.atoms) if v.molecule else 0
        use_simple_rendering = num_atoms > 10000 and not is_hovered and not is_custom

        if use_simple_rendering:
            # Simple filled circle for very large proteins
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(QColor(*rgb)))
            painter.drawEllipse(QPointF(sx, sy), max(2, radius), max(2, radius))
            return

        gradient, do_specular = self._get_cached_gradient(
            elem_sym, radius, sx, sy, rgb,
            is_hovered=is_hovered, use_ssao=use_ssao, sz=sz, alpha=alpha
        )
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QBrush(gradient))
        painter.drawEllipse(QRectF(sx - radius, sy - radius,
                                    radius * 2, radius * 2))

        # Specular for performance
        if do_specular:
            spec_radius = radius * 0.12
            spec_x = sx - radius * 0.3
            spec_y = sy - radius * 0.3
            spec_grad = QRadialGradient(QPointF(spec_x, spec_y), spec_radius,
                                         QPointF(spec_x, spec_y))
            
            opacity = 180 if is_hovered else 80
            spec_grad.setColorAt(0.0, QColor(255, 255, 255, int(opacity * alpha)))
            spec_grad.setColorAt(1.0, QColor(255, 255, 255, 0))
            painter.setBrush(QBrush(spec_grad))
            painter.drawEllipse(QRectF(spec_x - spec_radius,
                                        spec_y - spec_radius,
                                        spec_radius * 2,
                                        spec_radius * 2))

    def _draw_selection_ring(self, painter, sx, sy, radius):
        from src.features.visualization_3d.services.atom_rendering import draw_selection_ring
        draw_selection_ring(painter, sx, sy, radius)

    def _draw_sasa_surface(self, v, painter, width, height):
        from src.features.visualization_3d.services.atom_rendering import draw_sasa_surface
        draw_sasa_surface(painter, v.molecule, width, height,
                          v.pan_x, v.pan_y, v.rot_x, v.rot_y, v.zoom,
                          v.selected_atoms, getattr(v, 'show_sasa_selected_only', False))

    def _draw_atom_simple(self, painter, sx, sy, radius, color):
        """Draw atom as simple filled circle (fast for large molecules)."""
        painter.setPen(Qt.PenStyle.NoPen)
        # Convert RGB tuple to QColor if needed
        if isinstance(color, (tuple, list)) and len(color) >= 3:
            from PySide6.QtGui import QColor
            color = QColor(int(color[0]), int(color[1]), int(color[2]))
        painter.setBrush(QBrush(color))
        painter.drawEllipse(QPointF(sx, sy), max(2, radius), max(2, radius))

    # ─── Bond drawing ────────────────────────────────────────────

    def _draw_bonds(self, v, painter, projected, custom_only=False,
                    width=None, height=None, force_indices=None):
        """Draw bonds as lines/cylinders between atoms."""
        if not v.molecule:
            return

        proj_map = {p[0]: p for p in projected}
        cam = getattr(v, 'custom_atom_modes', {})

        # Pre-calculate ring information for 'inside-ring' rendering
        rings = v.molecule.find_rings()
        bond_to_ring = {}
        for ring in rings:
            for i in range(len(ring)):
                a, b = ring[i], ring[(i+1) % len(ring)]
                bond_to_ring[frozenset([a, b])] = ring

        # Culling bounds (with generous margin for long bonds)
        _bond_cull = width is not None and height is not None
        _bond_margin = 100

        for bond in v.molecule.bonds:
            i = bond.begin_atom_idx
            j = bond.end_atom_idx

            if i not in proj_map or j not in proj_map:
                continue

            if custom_only:
                is_forced = force_indices and (i in force_indices or j in force_indices)
                if not is_forced and i not in cam and j not in cam:
                    continue

            _, x1, y1, z1, r1, c1 = proj_map[i]
            _, x2, y2, z2, r2, c2 = proj_map[j]

            # Off-screen culling: skip if both endpoints are outside viewport
            if _bond_cull:
                if ((x1 < -_bond_margin and x2 < -_bond_margin) or
                        (x1 > width + _bond_margin and x2 > width + _bond_margin) or
                        (y1 < -_bond_margin and y2 < -_bond_margin) or
                        (y1 > height + _bond_margin and y2 > height + _bond_margin)):
                    continue

            avg_z = (z1 + z2) / 2
            depth_shade = max(0.35, min(1.0, 1.0 - avg_z * 0.025))

            bond_render_mode = v.render_mode
            cam = getattr(v, 'custom_atom_modes', {})

            is_custom = False
            if i in cam and j in cam:
                bond_render_mode = cam[i]
                is_custom = True
            elif i in cam:
                bond_render_mode = cam[i]
                is_custom = True
            elif j in cam:
                bond_render_mode = cam[j]
                is_custom = True

            base_width = max(2, v.zoom * 0.07 * v.stick_scale)
            if bond_render_mode == 'wireframe':
                base_width = max(1.5, v.zoom * 0.04 * v.line_scale)
            elif bond_render_mode == 'spacefill':
                base_width = max(1, v.zoom * 0.04 * v.stick_scale)
            elif bond_render_mode == 'ball_and_stick':
                base_width = max(3, v.zoom * 0.1 * v.stick_scale)

            if bond.is_double or bond.order == 2.0:
                dx = y2 - y1
                dy = -(x2 - x1)
                length = math.sqrt(dx * dx + dy * dy) + 1e-10
                dx /= length
                dy /= length
                offset = base_width * 0.9

                # Inside-ring heuristic
                ring = bond_to_ring.get(frozenset([i, j]))
                if ring:
                    # Projected ring center
                    rsx = sum(proj_map[idx][1] for idx in ring) / len(ring)
                    rsy = sum(proj_map[idx][2] for idx in ring) / len(ring)
                    # Midpoint of bond
                    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
                    # Dot product of offset vector (dx, dy) and vector to center (rsx-mx, rsy-my)
                    if dx * (rsx - mx) + dy * (rsy - my) < 0:
                        dx, dy = -dx, -dy

                for sign in (-1, 1):
                    ox = dx * offset * sign
                    oy = dy * offset * sign
                    self._draw_bond_line(painter, x1 + ox, y1 + oy, x2 + ox, y2 + oy,
                                         c1, c2, base_width * 0.55, depth_shade, is_custom=is_custom)

            elif bond.is_triple or bond.order == 3.0:
                dx = y2 - y1
                dy = -(x2 - x1)
                length = math.sqrt(dx * dx + dy * dy) + 1e-10
                dx /= length
                dy /= length
                offset = base_width * 1.3

                self._draw_bond_line(painter, x1, y1, x2, y2,
                                     c1, c2, base_width * 0.5, depth_shade, is_custom=is_custom)
                for sign in (-1, 1):
                    ox = dx * offset * sign
                    oy = dy * offset * sign
                    self._draw_bond_line(painter, x1 + ox, y1 + oy, x2 + ox, y2 + oy,
                                         c1, c2, base_width * 0.45, depth_shade, is_custom=is_custom)
            elif bond.is_aromatic or bond.order == 1.5:
                dx = y2 - y1
                dy = -(x2 - x1)
                length = math.sqrt(dx * dx + dy * dy) + 1e-10
                dx /= length
                dy /= length
                offset = base_width * 0.7

                # Inside-ring heuristic
                ring = bond_to_ring.get(frozenset([i, j]))
                if ring:
                    rsx = sum(proj_map[idx][1] for idx in ring) / len(ring)
                    rsy = sum(proj_map[idx][2] for idx in ring) / len(ring)
                    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
                    if dx * (rsx - mx) + dy * (rsy - my) < 0:
                        dx, dy = -dx, -dy

                self._draw_bond_line(painter, x1, y1, x2, y2,
                                     c1, c2, base_width * 0.7, depth_shade, is_custom=is_custom)

                ox = dx * offset
                oy = dy * offset
                self._draw_bond_line(painter, x1 + ox, y1 + oy, x2 + ox, y2 + oy,
                                     c1, c2, base_width * 0.4, depth_shade, dashed=True, is_custom=is_custom)
            else:
                self._draw_bond_line(painter, x1, y1, x2, y2,
                                     c1, c2, base_width, depth_shade, is_custom=is_custom)

    def _draw_bond_line(self, painter, x1, y1, x2, y2, c1, c2, width, shade, dashed=False, is_custom=False):
        from src.features.visualization_3d.services.atom_rendering import draw_bond_line
        draw_bond_line(painter, x1, y1, x2, y2, c1, c2, width, shade, dashed, is_custom)

    # ─── Labels ──────────────────────────────────────────────────

    def _draw_label(self, v, painter, atom_idx, sx, sy, radius):
        # Skip symbol if a custom label (e.g. lipophilicity) is already assigned to this atom
        if atom_idx in getattr(v, 'labels', {}):
            return
            
        from src.features.visualization_3d.services.atom_rendering import draw_label
        atom = v.molecule.atoms[atom_idx]
        draw_label(painter, atom.symbol, sx, sy, radius, v.label_font_size, getattr(v, '_export_scale', 1.0), v.label_color)

    def _draw_residue_label(self, v, painter, text, sx, sy, color, radius, settings=None):
        from src.features.visualization_3d.services.atom_rendering import draw_residue_label
        settings = settings or getattr(v, 'residue_label_settings', {})
        draw_residue_label(painter, text, sx, sy, color, radius, v.label_font_size, getattr(v, '_export_scale', 1.0), settings)

    def _draw_overlay(self, v, painter):
        from src.features.visualization_3d.services.overlay_rendering import draw_overlay
        draw_overlay(painter, v.molecule, v._hovered_atom)

    def _draw_placeholder(self, painter, width, height):
        from src.features.visualization_3d.services.overlay_rendering import draw_placeholder
        draw_placeholder(painter, width, height)

    # ─── Rubber-band rectangle ───────────────────────────────────

    def _draw_rubber_band(self, v, painter):
        """Draw the semi-transparent selection rectangle overlay."""
        x1 = min(v._sel_rect_origin.x(), v._sel_rect_end.x())
        y1 = min(v._sel_rect_origin.y(), v._sel_rect_end.y())
        x2 = max(v._sel_rect_origin.x(), v._sel_rect_end.x())
        y2 = max(v._sel_rect_origin.y(), v._sel_rect_end.y())
        rect = QRectF(x1, y1, x2 - x1, y2 - y1)

        # Semi-transparent yellow fill
        painter.setPen(QPen(QColor(255, 200, 50, 200), 1.5,
                            Qt.PenStyle.DashLine))
        painter.setBrush(QBrush(QColor(255, 200, 50, 40)))
        painter.drawRect(rect)

    # ─── Large molecule fast path ────────────────────────────────

    def _draw_large_molecule_fast(self, v, painter, projected, sorted_atoms,
                                   width=None, height=None):
        """Fast rendering for large molecules (>500 atoms)."""
        proj_map = {i: (atom_idx, sx, sy, sz, radius, color)
                    for atom_idx, sx, sy, sz, radius, color in projected
                    for i in [atom_idx]}

        _cull = width is not None and height is not None
        _margin = 100

        # Draw bonds as simple lines
        pen = QPen(QColor(100, 100, 100), 1)
        painter.setPen(pen)
        for bond in v.molecule.bonds:
            a1, a2 = bond.begin_atom_idx, bond.end_atom_idx
            if a1 in proj_map and a2 in proj_map:
                _, x1, y1, *_ = proj_map[a1]
                _, x2, y2, *_ = proj_map[a2]
                # Off-screen culling for bonds
                if _cull:
                    if ((x1 < -_margin and x2 < -_margin) or
                            (x1 > width + _margin and x2 > width + _margin) or
                            (y1 < -_margin and y2 < -_margin) or
                            (y1 > height + _margin and y2 > height + _margin)):
                        continue
                painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

        # Draw atoms as simple circles (no gradients), with culling
        _atom_margin = 50
        for atom_idx, sx, sy, sz, radius, color in sorted_atoms:
            if _cull and (sx < -_atom_margin or sx > width + _atom_margin or
                          sy < -_atom_margin or sy > height + _atom_margin):
                continue
            self._draw_atom_simple(painter, sx, sy, radius, color)

    def _draw_performance_indicator(self, v, painter, num_atoms):
        """Draw performance indicator for large molecules."""
        font = QFont('Segoe UI', 10)
        painter.setFont(font)

        indicator_text = f"Large molecule ({num_atoms} atoms) - Fast rendering active"
        text_rect = QRectF(10, v.height() - 30, 400, 25)

        painter.fillRect(text_rect, QColor(0, 0, 0, 120))

        painter.setPen(QColor(255, 255, 100))
        painter.drawText(text_rect, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                         indicator_text)

    def _draw_interaction_lines(self, v, painter, projected, width, height):
        """Draw dashed lines for molecular interactions."""
        proj_map = {p[0]: p for p in projected}
        
        for idx1, idx2, type_str, color_hex in v.interaction_lines:
            if idx1 not in proj_map or idx2 not in proj_map:
                continue
                
            _, x1, y1, z1, _, _ = proj_map[idx1]
            _, x2, y2, z2, _, _ = proj_map[idx2]
            
            # Off-screen culling
            _margin = 50
            if ((x1 < -_margin and x2 < -_margin) or
                (x1 > width + _margin and x2 > width + _margin) or
                (y1 < -_margin and y2 < -_margin) or
                (y1 > height + _margin and y2 > height + _margin)):
                continue
                
            color = QColor(color_hex)
            pen = QPen(color, 1.5, Qt.PenStyle.DashLine)
            painter.setPen(pen)
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

    def _draw_labels(self, v, painter, projected):
        """Draw text labels for specific atoms (e.g. lipophilic contributions)."""
        from src.features.visualization_3d.services.atom_rendering import draw_label
        proj_map = {p[0]: p for p in projected}
        
        for idx, text in v.labels.items():
            if idx in proj_map:
                _, sx, sy, sz, radius, _ = proj_map[idx]
                draw_label(painter, text, sx, sy, radius, v.label_font_size, getattr(v, '_export_scale', 1.0), v.label_color)

    # ─── Measurements ────────────────────────────────────────────

    def _draw_measurements(self, v, painter, projected):
        """Draw distance/angle measurement overlays."""
        proj_map = {p[0]: p for p in projected}

        # Draw dotted lines for pending picks
        if len(v._measure_atoms) >= 1:
            pen = QPen(QColor(255, 200, 50, 200), 2, Qt.PenStyle.DashLine)
            painter.setPen(pen)
            for m in range(len(v._measure_atoms) - 1):
                i = v._measure_atoms[m]
                j = v._measure_atoms[m + 1]
                if i in proj_map and j in proj_map:
                    _, x1, y1, *_ = proj_map[i]
                    _, x2, y2, *_ = proj_map[j]
                    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

        font = QFont('Segoe UI', 11)
        font.setBold(True)
        painter.setFont(font)

        for meas in v._measurements:
            if meas[0] == 'dist':
                _, i, j, dist = meas
                if i in proj_map and j in proj_map:
                    _, x1, y1, *_ = proj_map[i]
                    _, x2, y2, *_ = proj_map[j]
                    pen = QPen(QColor(50, 220, 255, 200), 2, Qt.PenStyle.DashDotLine)
                    painter.setPen(pen)
                    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
                    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
                    label = f"{dist:.2f} A"
                    painter.setPen(Qt.PenStyle.NoPen)
                    painter.setBrush(QBrush(QColor(0, 0, 0, 180)))
                    painter.drawRoundedRect(QRectF(mx - 30, my - 12, 60, 20), 4, 4)
                    painter.setPen(QColor(50, 220, 255))
                    painter.drawText(QRectF(mx - 30, my - 12, 60, 20),
                                     Qt.AlignmentFlag.AlignCenter, label)

            elif meas[0] == 'angle':
                _, i, j, k, angle = meas
                if j in proj_map:
                    _, vx, vy, *_ = proj_map[j]
                    if i in proj_map:
                        _, x1, y1, *_ = proj_map[i]
                        pen = QPen(QColor(255, 150, 50, 200), 2, Qt.PenStyle.DashDotLine)
                        painter.setPen(pen)
                        painter.drawLine(QPointF(vx, vy), QPointF(x1, y1))
                    if k in proj_map:
                        _, x3, y3, *_ = proj_map[k]
                        painter.drawLine(QPointF(vx, vy), QPointF(x3, y3))

                    label = f"{angle:.1f} deg"
                    painter.setPen(Qt.PenStyle.NoPen)
                    painter.setBrush(QBrush(QColor(0, 0, 0, 180)))
                    painter.drawRoundedRect(QRectF(vx - 35, vy - 30, 70, 20), 4, 4)
                    painter.setPen(QColor(255, 150, 50))
                    painter.drawText(QRectF(vx - 35, vy - 30, 70, 20),
                                     Qt.AlignmentFlag.AlignCenter, label)

    # ─── Protein Rendering ───────────────────────────────────────

    def _draw_protein(self, v, painter, projected, width, height):
        """Draw protein structures with cartoon/ribbon/backbone representations."""
        if not v.molecule:
            return

        # Performance optimization: Check mesh cache first
        mol_id = id(v.molecule) if v.molecule else 0
        current_time = time.perf_counter()
        
        # Check if we have a cached mesh that's still valid
        if hasattr(self, '_mesh_cache'):
            cached_data = self._mesh_cache.get(mol_id)
            if cached_data:
                vertices, triangles, colors, cache_time = cached_data
                # Use cache if it's less than 5 seconds old
                if current_time - cache_time < 5.0:
                    print(f"[Performance] Using cached mesh for {len(vertices)} vertices")
                    # Draw cached mesh directly
                    from src.features.visualization_3d.services.protein_rendering import draw_cached_mesh
                    draw_cached_mesh(painter, vertices, triangles, colors, v.zoom)
                    return

        try:
            from src.features.visualization_3d.services.protein_rendering import (
                render_protein_cartoon, render_protein_ribbon
            )

            if v.render_mode == 'cartoon':
                mesh_data = render_protein_cartoon(
                    painter=painter,
                    molecule=v.molecule,
                    width=width,
                    height=height,
                    rot_x=v.rot_x,
                    rot_y=v.rot_y,
                    pan_x=v.pan_x,
                    pan_y=v.pan_y,
                    zoom=v.zoom,
                    rot_z=getattr(v, 'rot_z', 0.0),
                    color_scheme="secondary_structure",
                    use_ssao=getattr(v, 'use_ssao', False),
                    use_gouraud=getattr(v, 'use_gouraud', False),
                    is_interacting=getattr(v, '_is_interacting', False)
                )
            elif v.render_mode == 'ribbon':
                mesh_data = render_protein_ribbon(
                    painter=painter,
                    molecule=v.molecule,
                    width=width,
                    height=height,
                    rot_x=v.rot_x,
                    rot_y=v.rot_y,
                    pan_x=v.pan_x,
                    pan_y=v.pan_y,
                    zoom=v.zoom,
                    rot_z=getattr(v, 'rot_z', 0.0),
                    color_scheme="rainbow"
                )
            elif v.render_mode == 'backbone':
                residues = self._group_residues(v)
                self._draw_backbone(v, painter, residues)
                return  # Backbone doesn't use mesh caching
            
            # Cache the mesh data for future use
            if mesh_data and hasattr(self, '_mesh_cache'):
                vertices, triangles, colors = mesh_data
                self._mesh_cache[mol_id] = (vertices, triangles, colors, current_time)
                print(f"[Performance] Cached mesh with {len(vertices)} vertices")

        except Exception as e:
            print(f"[Performance] Error in protein rendering: {e}")
            import traceback
            traceback.print_exc()

            # Optionally draw side chains
            if getattr(v, 'show_sidechains', False) or getattr(v, 'sidechain_res_vis', {}):
                self._draw_side_chains(v, painter, projected)

        except Exception as e:
            import traceback
            traceback.print_exc()
            residues = self._group_residues(v)
            if v.render_mode == 'cartoon':
                self._draw_cartoon(v, painter, residues)
            elif v.render_mode == 'ribbon':
                self._draw_ribbon(v, painter, residues)
            elif v.render_mode == 'backbone':
                self._draw_backbone(v, painter, residues)

    def _group_residues(self, v):
        """Group atoms into residues for protein rendering."""
        residues = {}

        for i, atom in enumerate(v.molecule.atoms):
            if not hasattr(atom, 'chain_id') or not hasattr(atom, 'res_seq'):
                continue

            chain_id = atom.chain_id or 'A'
            res_seq = atom.res_seq

            key = (chain_id, res_seq)
            if key not in residues:
                residues[key] = {
                    'chain_id': chain_id,
                    'res_seq': res_seq,
                    'res_name': getattr(atom, 'res_name', 'UNK'),
                    'atoms': [],
                    'ss_type': getattr(atom, 'ss_type', 'C')
                }

            residues[key]['atoms'].append(i)

            if residues[key]['ss_type'] == 'C' and hasattr(atom, 'ss_type'):
                residues[key]['ss_type'] = atom.ss_type

        return residues

    def _get_ss_color(self, ss_type):
        """Get color for secondary structure type."""
        colors = {
            'H': QColor(220, 50, 50),
            'E': QColor(50, 150, 220),
            'C': QColor(180, 180, 180),
        }
        return colors.get(ss_type, QColor(180, 180, 180))

    def _draw_cartoon(self, v, painter, residues):
        """Draw PyMOL-style cartoon representation."""
        if not residues or len(residues) == 0:
            return

        sorted_residues = sorted(residues.items(),
                                 key=lambda x: (x[0][0], x[0][1]))

        chains = {}
        for (chain_id, res_seq), residue in sorted_residues:
            if chain_id not in chains:
                chains[chain_id] = []

            ca_idx = None
            for atom_idx in residue['atoms']:
                atom = v.molecule.atoms[atom_idx]
                if hasattr(atom, 'pdb_name') and atom.pdb_name == 'CA':
                    ca_idx = atom_idx
                    break

            if ca_idx is None:
                continue

            ca_atom = v.molecule.atoms[ca_idx]
            if not ca_atom.has_coords:
                continue

            w, h = v.width(), v.height()
            cx = w / 2 + v.pan_x
            cy = h / 2 + v.pan_y

            cos_x = math.cos(math.radians(v.rot_x))
            sin_x = math.sin(math.radians(v.rot_x))
            cos_y = math.cos(math.radians(v.rot_y))
            sin_y = math.sin(math.radians(v.rot_y))
            cos_z = math.cos(math.radians(getattr(v, 'rot_z', 0.0)))
            sin_z = math.sin(math.radians(getattr(v, 'rot_z', 0.0)))

            x, y, z = ca_atom.x, ca_atom.y, ca_atom.z
            x1 = x * cos_y + z * sin_y
            z1 = -x * sin_y + z * cos_y
            y1 = y * cos_x - z1 * sin_x
            z2 = y * sin_x + z1 * cos_x

            x2 = x1 * cos_z - y1 * sin_z
            y2 = x1 * sin_z + y1 * cos_z

            sx = cx + x2 * v.zoom
            sy = cy - y2 * v.zoom

            chains[chain_id].append({
                'res_seq': res_seq,
                'ss_type': residue['ss_type'],
                'x': sx,
                'y': sy,
                'z': z2
            })

        for chain_id, points in chains.items():
            if len(points) < 2:
                continue
            self._draw_pyMOL_cartoon_chain(painter, points)

    def _draw_pyMOL_cartoon_chain(self, painter, points):
        """Draw PyMOL-style cartoon chain with proper secondary structure."""
        if not points or len(points) < 2:
            return

        segments = []
        current_segment = []
        current_ss = None

        for point in points:
            ss_type = point['ss_type']
            if current_ss is None:
                current_ss = ss_type
                current_segment = [point]
            elif ss_type == current_ss:
                current_segment.append(point)
            else:
                if len(current_segment) >= 1:
                    segments.append((current_ss, current_segment))
                current_ss = ss_type
                current_segment = [point]

        if len(current_segment) >= 1:
            segments.append((current_ss, current_segment))

        for ss_type, segment in segments:
            if len(segment) < 2:
                continue

            if ss_type == 'H':
                self._draw_pyMOL_helix(painter, segment)
            elif ss_type == 'E':
                self._draw_pyMOL_sheet(painter, segment)
            else:
                self._draw_pyMOL_coil(painter, segment)

    def _draw_pyMOL_helix(self, painter, points):
        """Draw PyMOL-style alpha helix with configurable color."""
        if len(points) < 2:
            return

        from src.shared.ui.theme import COLORS
        color = QColor(COLORS.get('ss_helix', '#dc3232'))

        path = QPainterPath()
        path.moveTo(points[0]['x'], points[0]['y'])

        for i in range(1, len(points)):
            path.lineTo(points[i]['x'], points[i]['y'])

        pen = QPen(color, 15)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        painter.setPen(pen)
        painter.drawPath(path)

        highlight_pen = QPen(QColor(255, 255, 255, 150), 6)
        highlight_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(highlight_pen)
        painter.drawPath(path)

    def _draw_pyMOL_sheet(self, painter, points):
        """Draw PyMOL-style beta sheet with configurable color."""
        if len(points) < 2:
            return

        from src.shared.ui.theme import COLORS
        color = QColor(COLORS.get('ss_sheet', '#3296dc'))

        path = QPainterPath()
        path.moveTo(points[0]['x'], points[0]['y'])

        for i in range(1, len(points)):
            path.lineTo(points[i]['x'], points[i]['y'])

        pen = QPen(color, 12)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        painter.setPen(pen)
        painter.drawPath(path)

        if len(points) >= 2:
            last_point = points[-1]
            second_last = points[-2]

            dx = last_point['x'] - second_last['x']
            dy = last_point['y'] - second_last['y']
            angle = math.atan2(dy, dx)

            arrow_length = 20
            arrow_width = 12

            base_x = last_point['x'] - arrow_length * math.cos(angle)
            base_y = last_point['y'] - arrow_length * math.sin(angle)

            perp_angle = angle + math.pi / 2
            offset_x = arrow_width / 2 * math.cos(perp_angle)
            offset_y = arrow_width / 2 * math.sin(perp_angle)

            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(color))

            rect_points = [
                QPointF(base_x - offset_x, base_y - offset_y),
                QPointF(base_x + offset_x, base_y + offset_y),
                QPointF(last_point['x'] + offset_x, last_point['y'] + offset_y),
                QPointF(last_point['x'] - offset_x, last_point['y'] - offset_y)
            ]

            polygon = QPolygonF(rect_points)
            painter.drawPolygon(polygon)

    def _draw_pyMOL_coil(self, painter, points):
        """Draw PyMOL-style coil with configurable color."""
        if len(points) < 2:
            return

        from src.shared.ui.theme import COLORS
        color = QColor(COLORS.get('ss_coil', '#b4b4b4'))

        path = QPainterPath()
        path.moveTo(points[0]['x'], points[0]['y'])

        for i in range(1, len(points)):
            path.lineTo(points[i]['x'], points[i]['y'])

        pen = QPen(color, 8)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        painter.setPen(pen)
        painter.drawPath(path)

    def _draw_ribbon(self, v, painter, residues):
        """Draw ribbon representation with smooth curves."""
        sorted_residues = sorted(residues.items(),
                                 key=lambda x: (x[0][0], x[0][1]))

        points_by_chain = {}

        for (chain_id, res_seq), residue in sorted_residues:
            ca_idx = None
            for atom_idx in residue['atoms']:
                atom = v.molecule.atoms[atom_idx]
                if hasattr(atom, 'pdb_name') and atom.pdb_name == 'CA':
                    ca_idx = atom_idx
                    break

            if ca_idx is None:
                continue

            ca_atom = v.molecule.atoms[ca_idx]
            if not ca_atom.has_coords:
                continue

            w, h = v.width(), v.height()
            cx = w / 2 + v.pan_x
            cy = h / 2 + v.pan_y

            cos_x = math.cos(math.radians(v.rot_x))
            sin_x = math.sin(math.radians(v.rot_x))
            cos_y = math.cos(math.radians(v.rot_y))
            sin_y = math.sin(math.radians(v.rot_y))

            x, y, z = ca_atom.x, ca_atom.y, ca_atom.z
            x1 = x * cos_y + z * sin_y
            z1 = -x * sin_y + z * cos_y
            y1 = y * cos_x - z1 * sin_x
            z2 = y * sin_x + z1 * cos_x

            sx = cx + x1 * v.zoom
            sy = cy - y1 * v.zoom

            if chain_id not in points_by_chain:
                points_by_chain[chain_id] = []

            points_by_chain[chain_id].append((sx, sy, residue['ss_type']))

        for chain_id, points in points_by_chain.items():
            if len(points) < 2:
                continue

            for i in range(len(points) - 1):
                x1, y1, ss1 = points[i]
                x2, y2, ss2 = points[i + 1]
                color = self._get_ss_color(ss1)
                self._draw_smooth_ribbon(painter, x1, y1, x2, y2, color, 6)

    def _draw_backbone(self, v, painter, residues):
        """Draw simple backbone trace."""
        sorted_residues = sorted(residues.items(),
                                 key=lambda x: (x[0][0], x[0][1]))

        prev_point = None

        for (chain_id, res_seq), residue in sorted_residues:
            ca_idx = None
            for atom_idx in residue['atoms']:
                atom = v.molecule.atoms[atom_idx]
                if hasattr(atom, 'pdb_name') and atom.pdb_name == 'CA':
                    ca_idx = atom_idx
                    break

            if ca_idx is None:
                continue

            ca_atom = v.molecule.atoms[ca_idx]
            if not ca_atom.has_coords:
                continue

            w, h = v.width(), v.height()
            cx = w / 2 + v.pan_x
            cy = h / 2 + v.pan_y

            cos_x = math.cos(math.radians(v.rot_x))
            sin_x = math.sin(math.radians(v.rot_x))
            cos_y = math.cos(math.radians(v.rot_y))
            sin_y = math.sin(math.radians(v.rot_y))

            x, y, z = ca_atom.x, ca_atom.y, ca_atom.z
            x1 = x * cos_y + z * sin_y
            z1 = -x * sin_y + z * cos_y
            y1 = y * cos_x - z1 * sin_x
            z2 = y * sin_x + z1 * cos_x

            sx = cx + x1 * v.zoom
            sy = cy - y1 * v.zoom

            if prev_point:
                color = self._get_ss_color(residue['ss_type'])
                pen = QPen(color, 3)
                pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                painter.setPen(pen)
                painter.drawLine(QPointF(*prev_point), QPointF(sx, sy))

            color = self._get_ss_color(residue['ss_type'])
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(color))
            painter.drawEllipse(QPointF(sx, sy), 3, 3)

            prev_point = (sx, sy)

    def _draw_side_chains(self, v, painter, projected):
        """Draw side chain atoms as sticks."""
        show_all = getattr(v, 'show_sidechains', False)
        vis_map = getattr(v, 'sidechain_res_vis', {})
        cam = getattr(v, 'custom_atom_modes', {})

        proj_map = {i: (atom_idx, sx, sy, sz, radius, color)
                    for atom_idx, sx, sy, sz, radius, color in projected
                    for i in [atom_idx]}

        for bond in v.molecule.bonds:
            a1, a2 = bond.begin_atom_idx, bond.end_atom_idx

            if a1 in cam or a2 in cam:
                continue

            atom1, atom2 = v.molecule.atoms[a1], v.molecule.atoms[a2]

            res_seq1 = getattr(atom1, 'res_seq', None)
            res_seq2 = getattr(atom2, 'res_seq', None)

            if not show_all:
                if not (vis_map.get(res_seq1, False) or vis_map.get(res_seq2, False)):
                    continue

            if (hasattr(atom1, 'pdb_name') and atom1.pdb_name in ['N', 'CA', 'C', 'O'] and
                    hasattr(atom2, 'pdb_name') and atom2.pdb_name in ['N', 'CA', 'C', 'O']):
                continue

            if a1 in proj_map and a2 in proj_map:
                _, x1, y1, *_ = proj_map[a1]
                _, x2, y2, *_ = proj_map[a2]
                pen = QPen(QColor(150, 150, 150), 2)
                pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                painter.setPen(pen)
                painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

        for atom_idx, sx, sy, sz, radius, color in projected:
            if atom_idx in cam:
                continue

            atom = v.molecule.atoms[atom_idx]

            res_seq = getattr(atom, 'res_seq', None)
            if not show_all and not vis_map.get(res_seq, False):
                continue

            if hasattr(atom, 'pdb_name') and atom.pdb_name in ['N', 'CA', 'C', 'O']:
                continue

            painter.setPen(Qt.PenStyle.NoPen)
            if isinstance(color, tuple):
                color = QColor(*color)
            painter.setBrush(QBrush(color))
            painter.drawEllipse(QPointF(sx, sy), 2, 2)

    # ─── Primitive helpers ───────────────────────────────────────

    def _draw_cylinder(self, painter, x1, y1, x2, y2, color, width):
        """Draw a cylinder between two points."""
        pen = QPen(color, width)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

    def _draw_arrow(self, painter, x1, y1, x2, y2, color, width):
        """Draw an arrow between two points."""
        pen = QPen(color, width)
        pen.setCapStyle(Qt.PenCapStyle.SquareCap)
        painter.setPen(pen)

        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

        angle = math.atan2(y2 - y1, x2 - x1)
        arrow_length = width * 2
        arrow_angle = 0.5

        ax1 = x2 - arrow_length * math.cos(angle - arrow_angle)
        ay1 = y2 - arrow_length * math.sin(angle - arrow_angle)
        ax2 = x2 - arrow_length * math.cos(angle + arrow_angle)
        ay2 = y2 - arrow_length * math.sin(angle + arrow_angle)

        painter.drawLine(QPointF(x2, y2), QPointF(ax1, ay1))
        painter.drawLine(QPointF(x2, y2), QPointF(ax2, ay2))

    def _draw_tube(self, painter, x1, y1, x2, y2, color, width):
        """Draw a smooth tube between two points."""
        pen = QPen(color, width)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

    def _draw_smooth_ribbon(self, painter, x1, y1, x2, y2, color, width):
        """Draw a smooth ribbon segment."""
        pen = QPen(color, width)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

    # ─── Dummy spheres ───────────────────────────────────────────

    def _draw_dummy_spheres(self, v, painter, width, height):
        """Draw dummy spheres (COM, centroid, custom) with alpha support and shell sorting."""
        if not v.molecule:
            return

        if not hasattr(v.molecule, 'dummy_spheres'):
            return

        dummy_spheres = v.molecule.dummy_spheres
        if not dummy_spheres:
            return

        cx = width / 2 + v.pan_x
        cy = height / 2 + v.pan_y

        cos_x = math.cos(math.radians(v.rot_x))
        sin_x = math.sin(math.radians(v.rot_x))
        cos_y = math.cos(math.radians(v.rot_y))
        sin_y = math.sin(math.radians(v.rot_y))

        spheres_with_depth = []
        for sphere in dummy_spheres:
            if not sphere.visible:
                continue

            x, y, z = sphere.position
            x1 = x * cos_y + z * sin_y
            z1 = -x * sin_y + z * cos_y
            y1 = y * cos_x - z1 * sin_x
            z2 = y * sin_x + z1 * cos_x

            spheres_with_depth.append((sphere, z2))

        sorted_spheres = sorted(spheres_with_depth, key=lambda s: (s[1], -s[0].radius))

        for sphere, sz in sorted_spheres:
            x, y, z = sphere.position
            radius = sphere.radius
            alpha = getattr(sphere, 'alpha', 1.0)

            from src.shared.ui.theme import COLORS
            color_hex = sphere.color
            if sphere.label == 'COM':
                color_hex = COLORS.get('sphere_com', color_hex)
            elif sphere.label == 'Centroid':
                color_hex = COLORS.get('sphere_centroid', color_hex)

            x1 = x * cos_y + z * sin_y
            z1 = -x * sin_y + z * cos_y
            y1 = y * cos_x - z1 * sin_x

            sx = cx + x1 * v.zoom
            sy = cy - y1 * v.zoom

            if getattr(sphere, 'label', '') in ['COM', 'Centroid']:
                display_r = radius * v.zoom * v.sphere_scale
            else:
                display_r = radius * v.zoom
            depth_factor = 1.0 + sz * 0.02
            display_r *= max(0.5, min(1.5, depth_factor))

            color_rgb = _hex_to_rgb(color_hex)

            self._draw_atom_sphere(v, painter, -1, sx, sy, sz, display_r, color_rgb, alpha=alpha)

            is_custom = hasattr(sphere, 'label') and sphere.label.lower() == 'custom'
            if hasattr(sphere, 'label') and sphere.label and alpha > 0.3 and not is_custom:
                painter.setPen(QColor(255, 255, 255, int(alpha * 255)))
                painter.setFont(QFont('Segoe UI', 8))
                label_rect = QRectF(sx + display_r + 5, sy - 10, 150, 20)
                painter.drawText(label_rect, Qt.AlignmentFlag.AlignLeft, sphere.label)
