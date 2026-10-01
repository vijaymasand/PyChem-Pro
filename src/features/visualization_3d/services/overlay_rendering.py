"""
Module for 2D UI overlays and measurements in PyChem's 3D viewer.
"""
from src.shared.qt_compat import Qt, QPointF, QRectF
from src.shared.qt_compat import QPainter, QColor, QPen, QBrush, QFont
from src.shared.ui.theme import COLORS

def draw_overlay(painter: QPainter, molecule, hovered_atom: int):
    """Draw text HUD overlay."""
    if not molecule:
        return

    font = QFont('Segoe UI', 11)
    painter.setFont(font)
    painter.setPen(QColor(COLORS['text_secondary']))

    y = 20
    texts = [
        f"Atoms: {len(molecule.atoms)}",
        f"Bonds: {len(molecule.bonds)}",
    ]

    if 0 <= hovered_atom < len(molecule.atoms):
        atom = molecule.atoms[hovered_atom]
        texts.append("-------------")
        texts.append(f"Atom: {atom.symbol}{hovered_atom + 1}")
        if atom.has_coords:
            texts.append(f"Pos: ({atom.x:.2f}, {atom.y:.2f}, {atom.z:.2f})")
        texts.append(f"Charge: {atom.partial_charge:.4f}")

    for text in texts:
        painter.drawText(10, y, text)
        y += 18

def draw_placeholder(painter: QPainter, width: int, height: int):
    """Draw placeholder text when no molecule is loaded."""
    font = QFont('Segoe UI', 16)
    painter.setFont(font)
    painter.setPen(QColor(COLORS['text_muted']))
    rect = QRectF(0, 0, width, height)
    painter.drawText(rect, Qt.AlignmentFlag.AlignCenter,
                     "Enter a SMILES string and click\n'Convert to 3D' to visualize")

def draw_measurements(painter: QPainter, measurements: list, measure_atoms: list, proj_map: dict):
    """Draw active measurements (distances/angles) as UI overlays."""
    # Draw dotted lines for pending picks
    if len(measure_atoms) >= 1:
        pen = QPen(QColor(255, 200, 50, 200), 2, Qt.PenStyle.DashLine)
        painter.setPen(pen)
        for m in range(len(measure_atoms) - 1):
            i = measure_atoms[m]
            j = measure_atoms[m + 1]
            if i in proj_map and j in proj_map:
                _, x1, y1, *_ = proj_map[i]
                _, x2, y2, *_ = proj_map[j]
                painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

    font = QFont('Segoe UI', 11)
    font.setBold(True)
    painter.setFont(font)

    for meas in measurements:
        if meas[0] == 'dist':
            _, i, j, dist = meas
            if i in proj_map and j in proj_map:
                _, x1, y1, *_ = proj_map[i]
                _, x2, y2, *_ = proj_map[j]
                # Dashed yellow line
                pen = QPen(QColor(50, 220, 255, 200), 2, Qt.PenStyle.DashDotLine)
                painter.setPen(pen)
                painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
                # Label at midpoint
                mx, my = (x1 + x2) / 2, (y1 + y2) / 2
                label = f"{dist:.2f} A"
                # Background
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
                # Draw from j
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
