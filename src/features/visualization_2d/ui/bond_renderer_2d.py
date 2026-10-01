"""
2D Bond Renderer — Robust bond drawing with stereochemistry and label clearance.

Handles single, double, triple, aromatic, and stereo (wedge/hash) bonds.
Calculates label clearance (shrink) to prevent bonds from overlapping atom labels.
"""

import math
from src.shared.qt_compat import Qt, QPointF, QPainterPath, QPen, QBrush, QPainter


class BondRenderer2D:
    """Stateless helper for drawing molecular bonds on a QPainter."""

    def __init__(self, viewer):
        self._v = viewer  # MolViewer2D instance

    def draw_all_bonds(self, painter, visible, has_label):
        """Iterate through all bonds and draw them if visible."""
        v = self._v
        bond_color = v._bond_color()
        bw = v._bond_width()

        # Group indices by bond set to avoid double drawing
        drawn_bonds = set()

        for i in visible:
            if i not in v.coords_2d:
                continue

            x1, y1 = v._to_screen(*v.coords_2d[i])
            neighbors = v.molecule.get_neighbors(i)

            for j in neighbors:
                if j not in visible or j not in v.coords_2d:
                    continue

                bond_idx = tuple(sorted((i, j)))
                if bond_idx in drawn_bonds:
                    continue
                drawn_bonds.add(bond_idx)

                bond = v.molecule.get_bond_between(i, j)
                x2, y2 = v._to_screen(*v.coords_2d[j])

                dx, dy = x2 - x1, y2 - y1
                dist = math.hypot(dx, dy)
                if dist < 2:
                    continue
                nx, ny = dx / dist, dy / dist

                # Shrink amount based on font metrics and bond direction
                shrink1 = self._label_shrink(i, has_label, dx, dy, painter)
                shrink2 = self._label_shrink(j, has_label, -dx, -dy, painter)

                if dist > (shrink1 + shrink2 + 4):
                    bx1, by1 = x1 + nx * shrink1, y1 + ny * shrink1
                    bx2, by2 = x2 - nx * shrink2, y2 - ny * shrink2
                else:
                    # Not enough room to draw bond
                    continue

                order = bond.order if bond else 1
                if order == 2:
                    self._draw_double(painter, bx1, by1, bx2, by2, bond_color, bw, bond)
                elif order == 3:
                    self._draw_triple(painter, bx1, by1, bx2, by2, bond_color, bw)
                elif order == 1.5:
                    self._draw_aromatic(painter, bx1, by1, bx2, by2, bond_color, bw, bond)
                else:
                    self._draw_single(painter, bx1, by1, bx2, by2, bond_color, bw, bond)

    # ─── Label shrink helper ──────────────────────────────────────

    def _label_shrink(self, idx, has_label, dx, dy, painter):
        """Compute how much to shrink bonds near labeled atoms.
        
        Uses RenderingConfig for consistent spacing and direction-aware logic
        to ensure bonds point to the anchor atom center.
        """
        v = self._v
        if not has_label.get(idx, False):
            return 0
        atom = v.molecule.atoms[idx]
        
        parts = v._atom_renderer._label_parts(atom)
        font = v._atom_renderer._get_font()
        sub_font = v._atom_renderer._get_subscript_font()

        painter.setFont(font)
        fm = painter.fontMetrics()
        painter.setFont(sub_font)
        fm_sub = painter.fontMetrics()
        
        if not parts:
            return v._label_padding
            
        first_text, _ = parts[0]
        first_w = fm.horizontalAdvance(first_text) if not parts[0][1] else fm_sub.horizontalAdvance(first_text)
        
        from .rendering_config import RenderingConfig
        char_gap, sub_gap, export_factor = RenderingConfig.get_gaps(v)

        total_w = 0
        for i, (text, is_sub) in enumerate(parts):
            m = fm_sub if is_sub else fm
            total_w += m.horizontalAdvance(text)
            if i > 0:
                total_w += sub_gap if is_sub else (char_gap if not parts[i-1][1] else 0)
        
        scaled_padding = int(v._label_padding * export_factor)
        h_offset = RenderingConfig.get_h_offset(export_factor)

        if dx > 0:
            # Bond leaves to the right
            clear_x = (total_w - first_w / 2 + scaled_padding) - h_offset
        else:
            # Bond leaves to the left
            clear_x = (first_w / 2 + scaled_padding) + h_offset

        # Approximate vertical clearance
        clear_y = fm.height() / 2.0 + scaled_padding

        dist = math.hypot(dx, dy)
        if dist < 0.001:
            return 0

        if abs(dx) > 0.001:
            t_x = clear_x * dist / abs(dx)
        else:
            t_x = float('inf')

        if abs(dy) > 0.001:
            t_y = clear_y * dist / abs(dy)
        else:
            t_y = float('inf')

        return min(t_x, t_y)

    # ─── Single bond ──────────────────────────────────────────────

    def _draw_single(self, painter, x1, y1, x2, y2, color, width, bond=None):
        if bond and hasattr(bond, 'stereo'):
            if bond.stereo == 'up':
                self._draw_wedge(painter, x1, y1, x2, y2, color, width)
                return
            elif bond.stereo == 'down':
                self._draw_hash(painter, x1, y1, x2, y2, color, width)
                return

        pen = QPen(color, width, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin)
        painter.setPen(pen)
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

    # ─── Stereo bonds ─────────────────────────────────────────────

    def _draw_wedge(self, painter, x1, y1, x2, y2, color, width):
        dx, dy = x2 - x1, y2 - y1
        dist = math.hypot(dx, dy)
        if dist < 2: return
        nx, ny = -dy / dist, dx / dist
        wedge_width = self._v._scale * self._v._wedge_width_ratio
        p1, p2, p3 = QPointF(x1, y1), QPointF(x2 - nx*wedge_width, y2 - ny*wedge_width), QPointF(x2 + nx*wedge_width, y2 + ny*wedge_width)
        path = QPainterPath(); path.moveTo(p1); path.lineTo(p2); path.lineTo(p3); path.closeSubpath()
        painter.setPen(Qt.PenStyle.NoPen); painter.setBrush(QBrush(color)); painter.drawPath(path)

    def _draw_hash(self, painter, x1, y1, x2, y2, color, width):
        dx, dy = x2 - x1, y2 - y1
        dist = math.hypot(dx, dy)
        if dist < 2: return
        nx, ny = -dy / dist, dx / dist
        hash_w = self._v._scale * self._v._wedge_width_ratio
        pen = QPen(color, width * 0.5, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        for i in range(1, 7):
            t = i / 6.0
            hx, hy = x1 + dx * t, y1 + dy * t
            hw = hash_w * t
            painter.drawLine(QPointF(hx - nx*hw, hy - ny*hw), QPointF(hx + nx*hw, hy + ny*hw))

    # ─── Double/Triple/Aromatic ───────────────────────────────────

    def _draw_double(self, painter, x1, y1, x2, y2, color, width, bond=None):
        dx, dy = x2 - x1, y2 - y1
        dist = math.hypot(dx, dy)
        if dist < 2: return
        ux, uy = dx/dist, dy/dist
        nx, ny = -uy, ux
        offset = self._v._scale * self._v._double_offset_ratio
        side = self._double_bond_side(bond, x1, y1, x2, y2, nx, ny)
        onx, ony = nx * offset * side, ny * offset * side

        painter.setPen(QPen(color, width, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap))
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

        shrink = 0.18
        painter.setPen(QPen(color, width * 0.85, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap))
        painter.drawLine(QPointF(x1 + dx*shrink + onx, y1 + dy*shrink + ony),
                         QPointF(x2 - dx*shrink + onx, y2 - dy*shrink + ony))

    def _draw_triple(self, painter, x1, y1, x2, y2, color, width):
        dx, dy = x2 - x1, y2 - y1
        dist = math.hypot(dx, dy)
        if dist < 2: return
        ux, uy = dx/dist, dy/dist
        nx, ny = -uy, ux
        offset = self._v._scale * self._v._triple_offset_ratio
        painter.setPen(QPen(color, width * 0.8, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap))
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
        for s in (1, -1):
            ox, oy = nx * offset * s, ny * offset * s
            painter.drawLine(QPointF(x1 + ox, y1 + oy), QPointF(x2 + ox, y2 + oy))

    def _draw_aromatic(self, painter, x1, y1, x2, y2, color, width, bond=None):
        dx, dy = x2 - x1, y2 - y1
        dist = math.hypot(dx, dy)
        if dist < 2: return
        ux, uy = dx/dist, dy/dist
        nx, ny = -uy, ux
        offset = self._v._scale * self._v._double_offset_ratio
        side = self._double_bond_side(bond, x1, y1, x2, y2, nx, ny)
        onx, ony = nx * offset * side, ny * offset * side
        painter.setPen(QPen(color, width, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap))
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
        shrink = 0.18
        pen_dash = QPen(color, max(1.0, width * 0.65), Qt.PenStyle.DashLine)
        pen_dash.setDashPattern([3.5, 2.5]); pen_dash.setCapStyle(Qt.PenCapStyle.FlatCap)
        painter.setPen(pen_dash)
        painter.drawLine(QPointF(x1 + dx*shrink + onx, y1 + dy*shrink + ony),
                         QPointF(x2 - dx*shrink + onx, y2 - dy*shrink + ony))

    def _double_bond_side(self, bond, x1, y1, x2, y2, nx, ny):
        v = self._v
        if not bond or not (hasattr(v.molecule, '_rings') and v.molecule._rings): return 1
        for ring in v.molecule._rings:
            if bond.begin_atom_idx in ring and bond.end_atom_idx in ring:
                cx_r = sum(v.coords_2d[r][0] for r in ring) / len(ring)
                cy_r = sum(v.coords_2d[r][1] for r in ring) / len(ring)
                csx, csy = v._to_screen(cx_r, cy_r)
                dot = (csx - (x1+x2)/2) * nx + (csy - (y1+y2)/2) * ny
                return 1 if dot > 0 else -1
        return 1
