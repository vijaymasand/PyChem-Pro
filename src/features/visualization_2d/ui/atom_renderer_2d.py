"""
2D Atom Renderer — ChemDraw-quality atom labels, charges, and coloring.

Handles atom symbols, implicit hydrogen display, element-specific coloring,
formal charge notation, and direction-aware H placement.
Extracted from mol_viewer_2d.py as a pure refactor (no behavior change).
"""

import math
from src.shared.qt_compat import Qt, QPointF, QRectF
from src.shared.qt_compat import (
    QPainter, QColor, QPen, QBrush, QFont, QFontMetrics,
)


# ── ChemDraw-style element colors ────────────────────────────────
# Against dark bg: brighter; against white bg: standard ChemDraw colors
ELEMENT_COLORS = {
    'C':  QColor(220, 220, 220),   # Near-white (hidden in skeletal)
    'H':  QColor(170, 170, 180),   # Light gray
    'N':  QColor(48, 128, 255),    # Vivid blue
    'O':  QColor(255, 50, 50),     # Vivid red
    'S':  QColor(230, 195, 50),    # Gold-yellow
    'P':  QColor(255, 140, 30),    # Warm orange
    'F':  QColor(80, 210, 80),     # Bright green
    'Cl': QColor(80, 210, 80),     # Green
    'Br': QColor(180, 80, 40),     # Dark orange-brown
    'I':  QColor(170, 50, 170),    # Purple
    'B':  QColor(255, 180, 160),   # Salmon
    'Se': QColor(255, 160, 0),     # Orange
    'Si': QColor(200, 180, 140),   # Tan
}

# Against white (export) background
ELEMENT_COLORS_LIGHT = {
    'C':  QColor(30, 30, 30),
    'H':  QColor(80, 80, 80),
    'N':  QColor(10, 50, 220),
    'O':  QColor(220, 20, 20),
    'S':  QColor(180, 150, 0),
    'P':  QColor(220, 100, 0),
    'F':  QColor(0, 160, 0),
    'Cl': QColor(0, 160, 0),
    'Br': QColor(150, 60, 20),
    'I':  QColor(130, 0, 130),
    'B':  QColor(200, 130, 100),
    'Se': QColor(200, 120, 0),
    'Si': QColor(150, 130, 90),
}


class AtomRenderer2D:
    """Stateless helper that draws atom labels on a QPainter."""

    def __init__(self, viewer):
        self._v = viewer  # MolViewer2D instance
        self._explicit_h_positions = {}
        self._explicit_h_angles = {}

    def _get_font(self):
        """Get the main atom label font, scaled properly."""
        from .rendering_config import RenderingConfig
        v = self._v
        # For exports: 
        # - when export_scale > 1 (high-DPI export): use original scale to prevent bloated fonts
        # - when export_scale < 1 (print preview): use current scale so fonts scale with molecule
        if v._original_scale is not None:
            export_factor = v._scale / v._original_scale if v._original_scale > 0 else 1.0
            if export_factor > 1.0:
                # High-DPI export: use original scale to keep font stable
                base_scale = v._original_scale
            else:
                # Print preview: use current scale so font scales with molecule
                base_scale = v._scale
        else:
            base_scale = v._scale
        size = int(base_scale * RenderingConfig.FONT_SIZE_RATIO)
        size = max(RenderingConfig.MIN_FONT_SIZE, min(size, RenderingConfig.MAX_FONT_SIZE))
        font = QFont('Arial', size)
        font.setWeight(QFont.Weight.DemiBold)
        return font

    def _get_subscript_font(self):
        """Get the subscript font scaled properly."""
        v = self._v
        if v._original_scale is not None:
            export_factor = v._scale / v._original_scale if v._original_scale > 0 else 1.0
            if export_factor > 1.0:
                base_scale = v._original_scale
            else:
                base_scale = v._scale
        else:
            base_scale = v._scale
        size = max(6, int(base_scale * 0.20))
        if v._original_scale is not None:
            size = min(size, 10)
        font = QFont('Arial', size)
        font.setWeight(QFont.Weight.DemiBold)
        return font

    def _get_h_placement_side(self, atom_idx, visible):
        """Determine H placement side (right/left)."""
        v = self._v
        if atom_idx not in v.coords_2d: return 'right'
        ax, ay = v.coords_2d[atom_idx]
        neighbors = v.molecule.get_neighbors(atom_idx)
        v_neighbors = [n for n in neighbors if n in visible and n in v.coords_2d]
        if not v_neighbors: return 'right'
        sum_dx = sum((v.coords_2d[n][0] - ax) for n in v_neighbors)
        return 'left' if sum_dx > 0.05 else 'right'

    def draw_all_labels(self, painter, visible, has_label):
        v = self._v
        font = self._get_font()
        sub_font = self._get_subscript_font()
        
        painter.setFont(font)
        fm = painter.fontMetrics()
        
        painter.setFont(sub_font)
        fm_sub = painter.fontMetrics()

        from .rendering_config import RenderingConfig
        char_gap, sub_gap, export_factor = RenderingConfig.get_gaps(v)

        for idx in visible:
            if not has_label.get(idx, False):
                continue
            atom = v.molecule.atoms[idx]
            mx, my = v.coords_2d[idx]
            sx, sy = v._to_screen(mx, my)

            parts = self._label_parts(atom)
            color = self._element_color(atom.symbol)

            total_w = 0
            for i, (text, is_sub) in enumerate(parts):
                m = fm_sub if is_sub else fm
                total_w += m.horizontalAdvance(text)
                if i > 0:
                    if is_sub:
                        total_w += sub_gap
                    elif not parts[i-1][1]:
                        total_w += char_gap

            h = fm.height()
            scaled_padding = int(v._label_padding * export_factor)
            pad = scaled_padding + (1 if export_factor <= 1.0 else int(export_factor))

            # Stable horizontal centering using font's natural advance
            first_text, first_is_sub = parts[0]
            m_first = fm_sub if first_is_sub else fm
            first_w = m_first.horizontalAdvance(first_text)
            
            # Start cx so that the center of the heavy atom symbol is at sx (with a slight shift for export balance)
            cx = sx - first_w / 2 + RenderingConfig.get_h_offset(export_factor)
            
            # Draw background rectangle aligned with the actual text layout
            # Use tighter padding to prevent covering bonds (especially terminal CH3)
            bg_pad = max(1, int(1.5 * export_factor))
            bg_rect = QRectF(cx - bg_pad, sy - h / 2,
                             total_w + bg_pad * 2, h + 1)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(v.bg_color))
            painter.drawRect(bg_rect)

            # Use centralized vertical offset for consistent heteroatom placement
            baseline_y = sy + RenderingConfig.get_v_offset(fm, export_factor)

            for i, (text, is_sub) in enumerate(parts):
                f = sub_font if is_sub else font
                m = fm_sub if is_sub else fm
                tw = m.horizontalAdvance(text)

                if i > 0:
                    if is_sub:
                        cx += sub_gap
                    elif not parts[i-1][1]:
                        cx += char_gap

                painter.setFont(f)
                painter.setPen(color)

                if is_sub:
                    painter.drawText(
                        QPointF(cx, baseline_y + fm.descent() + (2 * export_factor)), text)
                else:
                    painter.drawText(QPointF(cx, baseline_y), text)
                cx += tw

        self._draw_explicit_h_atoms(painter, visible, has_label, font, fm, export_factor > 1.0, export_factor)

    def _draw_explicit_h_atoms(self, painter, visible, has_label, font, fm, is_export, export_factor):
        v = self._v
        h_color = self._element_color('H')
        bond_color = v._bond_color()
        bw = v._bond_width()
        h = fm.height()
        from .rendering_config import RenderingConfig
        baseline_offset = RenderingConfig.get_v_offset(fm, export_factor)
        h_bond_len = v._scale * 0.70

        self._explicit_h_positions.clear()

        for idx in visible:
            if not has_label.get(idx, False): continue
            atom = v.molecule.atoms[idx]
            if atom.symbol in ('C', 'H'): continue
            h_count = self._get_total_h_count(atom)
            if h_count == 0: continue

            mx, my = v.coords_2d[idx]
            sx, sy = v._to_screen(mx, my)
            base_symbol_w = fm.horizontalAdvance(atom.symbol) / 2
            scaled_padding = v._label_padding * export_factor
            label_half_w = base_symbol_w + scaled_padding
            h_positions = self._compute_h_positions(idx, visible, h_count, h_bond_len)

            for h_idx, (h_pos_x, h_pos_y) in enumerate(h_positions):
                h_key = (idx, h_idx)
                self._explicit_h_positions[h_key] = (h_pos_x, h_pos_y, sx, sy)
                dx, dy = h_pos_x - sx, h_pos_y - sy
                dist = math.hypot(dx, dy)
                if dist < 1: continue

                nx, ny = dx/dist, dy/dist
                bond_start_x, bond_start_y = sx + nx * label_half_w, sy + ny * label_half_w
                h_tw = fm.horizontalAdvance('H')
                h_label_half = h_tw / 2 + scaled_padding
                bond_end_x, bond_end_y = h_pos_x - nx * h_label_half, h_pos_y - ny * h_label_half

                pen = QPen(bond_color, bw * 0.8, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap)
                painter.setPen(pen)
                painter.drawLine(QPointF(bond_start_x, bond_start_y), QPointF(bond_end_x, bond_end_y))

                tw = fm.horizontalAdvance('H')
                bg_rect = QRectF(h_pos_x - tw / 2 - scaled_padding, h_pos_y - h / 2 - 2, tw + scaled_padding * 2, h + 4)
                painter.setPen(Qt.PenStyle.NoPen)
                painter.setBrush(QBrush(v.bg_color))
                painter.drawRect(bg_rect)

                painter.setFont(font); painter.setPen(h_color)
                painter.drawText(QPointF(h_pos_x - tw / 2, h_pos_y + baseline_offset), 'H')

    def _compute_h_positions(self, atom_idx, visible, h_count, h_bond_len):
        v = self._v
        mx, my = v.coords_2d[atom_idx]
        sx, sy = v._to_screen(mx, my)
        
        custom_angles = [self._explicit_h_angles.get((atom_idx, i)) for i in range(h_count)]
        if all(a is not None for a in custom_angles):
            return [(sx + h_bond_len * math.cos(a), sy + h_bond_len * math.sin(a)) for a in custom_angles]

        neighbors = v.molecule.get_neighbors(atom_idx)
        v_neighbors = [n for n in neighbors if n in visible and n in v.coords_2d]
        bond_angles = []
        for n_idx in v_neighbors:
            nx, ny = v.coords_2d[n_idx]
            nsx, nsy = v._to_screen(nx, ny)
            bond_angles.append(math.atan2(nsy - sy, nsx - sx))

        positions = []
        if not bond_angles:
            for i in range(h_count):
                a = math.pi / 2 + (i - (h_count - 1) / 2) * (math.pi / 6)
                positions.append((sx + h_bond_len * math.cos(a), sy + h_bond_len * math.sin(a)))
        elif len(bond_angles) == 1:
            base = bond_angles[0] + math.pi
            for i in range(h_count):
                # Spread hydrogens by 120 degrees (2*pi/3) for better aesthetics
                spread = (math.pi * 2 / 3) if h_count > 1 else (math.pi / 3)
                a = base + (i - (h_count - 1) / 2) * spread
                positions.append((sx + h_bond_len * math.cos(a), sy + h_bond_len * math.sin(a)))
        else:
            bond_angles.sort()
            gaps = []
            for i in range(len(bond_angles)):
                a1, a2 = bond_angles[i], bond_angles[(i + 1) % len(bond_angles)]
                gaps.append(((a2 - a1) % (2 * math.pi), a1))
            gaps.sort(key=lambda g: g[0], reverse=True)
            for i in range(min(h_count, len(gaps))):
                a = gaps[i][1] + gaps[i][0] / 2
                positions.append((sx + h_bond_len * math.cos(a), sy + h_bond_len * math.sin(a)))

        return positions

    def _get_total_h_count(self, atom):
        v = self._v
        h_count = atom.num_implicit_h
        for n_idx in v.molecule.get_neighbors(atom.index):
            if v.molecule.atoms[n_idx].symbol == 'H': h_count += 1
        return h_count

    def _label_parts(self, atom):
        h_count = self._get_total_h_count(atom)
        if atom.symbol != 'C' and atom.symbol != 'H':
            parts = [(atom.symbol, False)]
            if atom.formal_charge != 0:
                ch = ('+' if atom.formal_charge > 0 else '') + str(atom.formal_charge)
                if atom.formal_charge == 1: ch = '+'
                if atom.formal_charge == -1: ch = '-'
                parts.append((ch, True))
            return parts

        parts = [(atom.symbol, False)]
        if h_count == 1: parts.append(('H', False))
        elif h_count > 1: parts.extend([('H', False), (str(h_count), True)])
        if atom.formal_charge != 0:
            ch = ('+' if atom.formal_charge > 0 else '') + str(atom.formal_charge)
            if atom.formal_charge == 1: ch = '+'
            if atom.formal_charge == -1: ch = '-'
            parts.append((ch, True))
        return parts

    def _build_label(self, atom):
        if atom.symbol not in ('C', 'H'): return atom.symbol
        label = atom.symbol
        h_count = self._get_total_h_count(atom)
        if h_count == 1: label += 'H'
        elif h_count > 1: label += f'H{h_count}'
        if atom.formal_charge != 0:
            ch = ('+' if atom.formal_charge > 0 else '') + str(atom.formal_charge)
            if atom.formal_charge == 1: ch = '+'
            if atom.formal_charge == -1: ch = '-'
            label += ch
        return label

    def _element_color(self, symbol):
        v = self._v
        palette = ELEMENT_COLORS if v._is_dark_bg else ELEMENT_COLORS_LIGHT
        return palette.get(symbol, QColor(200, 200, 200) if v._is_dark_bg else QColor(30, 30, 30))

    def hit_test_explicit_h(self, pos, tolerance=12):
        tol_sq = tolerance**2
        for h_key, (hx, hy, _, _) in self._explicit_h_positions.items():
            if (pos.x()-hx)**2 + (pos.y()-hy)**2 < tol_sq: return h_key
        return None

    def update_h_position(self, parent_idx, h_idx, new_pos):
        v = self._v
        if parent_idx not in v.coords_2d: return
        mx, my = v.coords_2d[parent_idx]
        sx, sy = v._to_screen(mx, my)
        self._explicit_h_angles[(parent_idx, h_idx)] = math.atan2(new_pos.y() - sy, new_pos.x() - sx)

    def get_explicit_h_data(self, h_key): return self._explicit_h_positions.get(h_key)
    def clear_explicit_h_angles(self): self._explicit_h_angles.clear()
