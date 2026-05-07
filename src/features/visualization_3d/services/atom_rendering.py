"""
Module for 3D rendering of atoms, bonds, and molecular surfaces in PyChem.
"""
import math
import numpy as np

from src.shared.qt_compat import Qt, QPointF, QRectF
from src.shared.qt_compat import QPainter, QColor, QPen, QBrush, QFont, QRadialGradient
from src.shared.ui.theme import COLORS

def _hex_to_rgb(hex_str: str):
    """Convert hex string (e.g. #FF0000) to (R, G, B) tuple."""
    hex_str = hex_str.lstrip('#')
    return tuple(int(hex_str[i:i+2], 16) for i in (0, 2, 4))

def draw_atom_sphere(painter: QPainter, sx: float, sy: float, sz: float, radius: float, rgb: tuple, 
                     is_hovered: bool = False, use_ssao: bool = False, alpha: float = 1.0):
    """
    Draw a single atom as a smooth, realistic sphere using QRadialGradient.
    """
    r, g, b = rgb

    # Depth-based shading
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

    highlight_x = sx - radius * 0.30
    highlight_y = sy - radius * 0.30

    gradient = QRadialGradient(QPointF(sx, sy), radius, QPointF(highlight_x, highlight_y))

    highlight_r = min(255, r + 160)
    highlight_g = min(255, g + 160)
    highlight_b = min(255, b + 160)

    shadow_r = max(0, int(r * 0.20))
    shadow_g = max(0, int(g * 0.20))
    shadow_b = max(0, int(b * 0.20))

    mid_r = max(0, min(255, int(r * 0.85)))
    mid_g = max(0, min(255, int(g * 0.85)))
    mid_b = max(0, min(255, int(b * 0.85)))

    alpha_int = int(alpha * 255)
    gradient.setColorAt(0.0, QColor(highlight_r, highlight_g, highlight_b, alpha_int))
    gradient.setColorAt(0.25, QColor(r, g, b, alpha_int))
    gradient.setColorAt(0.7, QColor(mid_r, mid_g, mid_b, alpha_int))
    gradient.setColorAt(1.0, QColor(shadow_r, shadow_g, shadow_b, alpha_int))

    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(QBrush(gradient))
    painter.drawEllipse(QRectF(sx - radius, sy - radius, radius * 2, radius * 2))

    # Extra specular dot
    if radius > 6 and alpha > 0.5:
        spec_radius = radius * 0.18
        spec_x = sx - radius * 0.28
        spec_y = sy - radius * 0.28

        spec_grad = QRadialGradient(QPointF(spec_x, spec_y), spec_radius, QPointF(spec_x, spec_y))
        spec_grad.setColorAt(0.0, QColor(255, 255, 255, int(180 * alpha)))
        spec_grad.setColorAt(0.5, QColor(255, 255, 255, int(60 * alpha)))
        spec_grad.setColorAt(1.0, QColor(255, 255, 255, 0))

        painter.setBrush(QBrush(spec_grad))
        painter.drawEllipse(QRectF(spec_x - spec_radius, spec_y - spec_radius, spec_radius * 2, spec_radius * 2))

def draw_selection_ring(painter: QPainter, sx: float, sy: float, radius: float):
    """Draw a glowing yellow ring around a selected atom."""
    ring_r = radius + max(3, radius * 0.25)
    pen = QPen(QColor(255, 200, 50, 200), max(2, radius * 0.12))
    painter.setPen(pen)
    painter.setBrush(QBrush(QColor(255, 200, 50, 30)))
    painter.drawEllipse(QRectF(sx - ring_r, sy - ring_r, ring_r * 2, ring_r * 2))


def draw_label(painter: QPainter, label: str, sx: float, sy: float, radius: float, base_font_size: int, export_scale: float = 1.0):
    font_size = base_font_size
    if export_scale > 1.0:
        scale_factor = min(export_scale * 0.3, 0.8)
        font_size = max(6, int(font_size * scale_factor))

    font = QFont('Segoe UI', font_size)
    font.setBold(True)
    painter.setFont(font)

    offset_x = int(radius * 0.5 + 2)
    offset_y = int(-radius * 0.3)

    painter.setPen(QColor(0, 0, 0, 180))
    painter.drawText(int(sx + offset_x + 1), int(sy + offset_y + 1), label)
    painter.setPen(QColor(255, 255, 255, 230))
    painter.drawText(int(sx + offset_x), int(sy + offset_y), label)

def draw_residue_label(painter: QPainter, text: str, sx: float, sy: float, color: QColor, radius: float, base_font_size: int, export_scale: float = 1.0, settings: dict = None):
    settings = settings or {}
    font_family = settings.get('font_family', 'Segoe UI')
    font_size = settings.get('font_size', base_font_size + 2)
    is_bold = settings.get('bold', True)
    is_italic = settings.get('italic', False)
    draw_bg = settings.get('background', False)

    if export_scale > 1.0:
        scale_factor = min(export_scale * 0.3, 0.8)
        font_size = max(8, int(font_size * scale_factor))

    font = QFont(font_family, font_size)
    font.setBold(is_bold)
    font.setItalic(is_italic)
    painter.setFont(font)

    offset_x = int(radius * 0.5 + 4)
    offset_y = int(-radius * 0.3 - 4)

    metrics = painter.fontMetrics()
    text_rect = metrics.boundingRect(text)
    text_rect.moveTo(int(sx + offset_x), int(sy + offset_y) - metrics.ascent())

    if draw_bg:
        bg_rect = text_rect.adjusted(-2, -2, 2, 2)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(255, 255, 255, 200))
        painter.drawRoundedRect(bg_rect, 3, 3)

    painter.setPen(QColor(0, 0, 0, 200))
    painter.drawText(int(sx + offset_x + 1), int(sy + offset_y + 1), text)

    painter.setPen(QColor(color))
    painter.drawText(int(sx + offset_x), int(sy + offset_y), text)

def draw_sasa_surface(painter: QPainter, molecule, width: float, height: float, pan_x: float, pan_y: float, 
                      rot_x: float, rot_y: float, zoom: float, selected_atoms: set, show_selected_only: bool):
    """Projects and renders the 3D analytical point-cloud surface for solvent accessibility."""
    if not molecule: return
        
    cx = width / 2 + pan_x
    cy = height / 2 + pan_y

    cos_x = math.cos(math.radians(rot_x))
    sin_x = math.sin(math.radians(rot_x))
    cos_y = math.cos(math.radians(rot_y))
    sin_y = math.sin(math.radians(rot_y))
    
    dot_r = max(1, min(4, zoom * 0.065))
    
    painter.setPen(Qt.PenStyle.NoPen)
    for atom in molecule.atoms:
        if not getattr(atom, 'sasa_points', None):
            continue
            
        if show_selected_only and selected_atoms:
            if atom.index not in selected_atoms:
                continue
            
        color = _hex_to_rgb(atom.element.color)
        brush = QBrush(QColor(color[0], color[1], color[2], 120))
        painter.setBrush(brush)
        
        for (x, y, z) in atom.sasa_points:
            x1 = x * cos_y + z * sin_y
            z1 = -x * sin_y + z * cos_y
            y1 = y * cos_x - z1 * sin_x
            z2 = y * sin_x + z1 * cos_x

            sx = cx + x1 * zoom
            sy = cy - y1 * zoom
            
            painter.drawRect(QRectF(sx - dot_r/2, sy - dot_r/2, dot_r, dot_r))

def draw_bond_line(painter: QPainter, x1: float, y1: float, x2: float, y2: float, 
                   c1: tuple, c2: tuple, width: float, shade: float, 
                   dashed: bool = False, is_custom: bool = False):
    """Draw a split-colored bond line with rounded caps."""
    mx = (x1 + x2) / 2
    my = (y1 + y2) / 2

    use_theme_stick_colors = ('stick_default' in COLORS) and not is_custom
    
    if use_theme_stick_colors:
        stick_color_hex = COLORS.get('stick_default', '#808080')
        r, g, b = _hex_to_rgb(stick_color_hex)
        
        color = QColor(max(0, min(255, int(r*shade))), max(0, min(255, int(g*shade))), max(0, min(255, int(b*shade))))
        pen = QPen(color, max(1, width))
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        if dashed:
            pen.setStyle(Qt.PenStyle.DashLine)
        painter.setPen(pen)
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
    else:
        r1, g1, b1 = c1
        color1 = QColor(max(0, min(255, int(r1*shade))), max(0, min(255, int(g1*shade))), max(0, min(255, int(b1*shade))))
        pen1 = QPen(color1, max(1, width))
        pen1.setCapStyle(Qt.PenCapStyle.RoundCap)
        if dashed: pen1.setStyle(Qt.PenStyle.DashLine)
        painter.setPen(pen1)
        painter.drawLine(QPointF(x1, y1), QPointF(mx, my))

        r2, g2, b2 = c2
        color2 = QColor(max(0, min(255, int(r2*shade))), max(0, min(255, int(g2*shade))), max(0, min(255, int(b2*shade))))
        pen2 = QPen(color2, max(1, width))
        pen2.setCapStyle(Qt.PenCapStyle.RoundCap)
        if dashed: pen2.setStyle(Qt.PenStyle.DashLine)
        painter.setPen(pen2)
        painter.drawLine(QPointF(mx, my), QPointF(x2, y2))
