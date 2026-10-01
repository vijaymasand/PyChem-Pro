# -*- coding: utf-8 -*-
from .drawing_parents import DrawableObject, Color
from . import geometry as geo
from math import pi, cos, sin, atan2, sqrt

class Arrow(DrawableObject):
    meta__undo_properties = ("p1", "p2", "type", "curvature", "color")
    def __init__(self, p1, p2, type="reaction", curvature=1.0):
        DrawableObject.__init__(self)
        self.p1 = p1
        self.p2 = p2
        self.type = type
        self.curvature = curvature
        self._items = []

    def draw(self):
        self.clear_drawings()
        if not self.paper: return
        
        if self.type == "reaction":
            self._draw_reaction_arrow()
        elif self.type == "curly":
            self._draw_curly_arrow()
        elif self.type == "equilibrium":
            self._draw_parallel_arrow(heads="half")
        elif self.type == "reversible":
            self._draw_parallel_arrow(heads="full")
        elif self.type == "fish_up":
            self._draw_curly_arrow(head_type="fish_up")
        elif self.type == "fish_down":
            self._draw_curly_arrow(head_type="fish_down")

        for item in self._items:
            if hasattr(self.paper, 'addFocusable'):
                self.paper.addFocusable(item, self)

    def _draw_reaction_arrow(self):
        # Main line
        line = self.p1 + self.p2
        self._items.append(self.paper.addLine(line, color=self.color, width=2))
        
        # Arrow head (Open style - two lines)
        angle = atan2(self.p2[1] - self.p1[1], self.p2[0] - self.p1[0])
        head_len = 10
        head_angle = pi / 6 
        h1 = (self.p2[0] - head_len * cos(angle - head_angle),
              self.p2[1] - head_len * sin(angle - head_angle))
        h2 = (self.p2[0] - head_len * cos(angle + head_angle),
              self.p2[1] - head_len * sin(angle + head_angle))
        
        self._items.append(self.paper.addLine(self.p2 + h1, color=self.color, width=2))
        self._items.append(self.paper.addLine(self.p2 + h2, color=self.color, width=2))

    def _draw_curly_arrow(self, head_type="full"):
        # Quadratic Bezier curve
        mid_x = (self.p1[0] + self.p2[0]) / 2
        mid_y = (self.p1[1] + self.p2[1]) / 2
        
        angle = atan2(self.p2[1] - self.p1[1], self.p2[0] - self.p1[0])
        dist = sqrt((self.p2[0] - self.p1[0])**2 + (self.p2[1] - self.p1[1])**2)
        offset = dist * 0.3
        
        cp_x = mid_x - offset * sin(angle) * self.curvature
        cp_y = mid_y + offset * cos(angle) * self.curvature
        
        from src.shared.qt_compat import QPainterPath, QPointF, QGraphicsPathItem, QPen, QColor
        path = QPainterPath(QPointF(*self.p1))
        path.quadTo(QPointF(cp_x, cp_y), QPointF(*self.p2))
        
        item = QGraphicsPathItem(path)
        item.setPen(QPen(QColor(*self.color), 1.5))
        self.paper.addItem(item)
        self._items.append(item)
        
        # Head
        tangent_angle = atan2(self.p2[1] - cp_y, self.p2[0] - cp_x)
        head_len = 10
        if head_type in ("full", "fish_up"):
            h1 = (self.p2[0] - head_len * cos(tangent_angle - pi/6),
                  self.p2[1] - head_len * sin(tangent_angle - pi/6))
            self._items.append(self.paper.addLine(self.p2 + h1, color=self.color, width=1.5))
        if head_type in ("full", "fish_down"):
            h2 = (self.p2[0] - head_len * cos(tangent_angle + pi/6),
                  self.p2[1] - head_len * sin(tangent_angle + pi/6))
            self._items.append(self.paper.addLine(self.p2 + h2, color=self.color, width=1.5))

    def _draw_parallel_arrow(self, heads="half"):
        from src.shared.qt_compat import QPainterPath, QPointF, QGraphicsPathItem, QPen, QColor
        
        angle = atan2(self.p2[1] - self.p1[1], self.p2[0] - self.p1[0])
        dist = sqrt((self.p2[0] - self.p1[0])**2 + (self.p2[1] - self.p1[1])**2)
        offset_dist = dist * 0.3 * self.curvature
        
        dx = 3 * sin(angle)
        dy = 3 * cos(angle)
        
        # We draw two parallel Bezier curves
        mid_x = (self.p1[0] + self.p2[0]) / 2
        mid_y = (self.p1[1] + self.p2[1]) / 2
        
        cp_x_base = mid_x - offset_dist * sin(angle)
        cp_y_base = mid_y + offset_dist * cos(angle)

        # Line 1: Forward
        p1a = (self.p1[0] - dx, self.p1[1] + dy)
        p2a = (self.p2[0] - dx, self.p2[1] + dy)
        cp_a = (cp_x_base - dx, cp_y_base + dy)
        
        path1 = QPainterPath(QPointF(*p1a))
        path1.quadTo(QPointF(*cp_a), QPointF(*p2a))
        item1 = QGraphicsPathItem(path1)
        item1.setPen(QPen(QColor(*self.color), 1.5))
        self.paper.addItem(item1)
        self._items.append(item1)
        
        # Head 1
        tangent1 = atan2(p2a[1] - cp_a[1], p2a[0] - cp_a[0])
        head_len = 8
        h1 = (p2a[0] - head_len * cos(tangent1 - pi/6), p2a[1] - head_len * sin(tangent1 - pi/6))
        self._items.append(self.paper.addLine(p2a + h1, color=self.color, width=1.5))
        if heads == "full":
            h1b = (p2a[0] - head_len * cos(tangent1 + pi/6), p2a[1] - head_len * sin(tangent1 + pi/6))
            self._items.append(self.paper.addLine(p2a + h1b, color=self.color, width=1.5))

        # Line 2: Backward
        p1b = (self.p1[0] + dx, self.p1[1] - dy)
        p2b = (self.p2[0] + dx, self.p2[1] - dy)
        cp_b = (cp_x_base + dx, cp_y_base - dy)
        
        path2 = QPainterPath(QPointF(*p2b))
        path2.quadTo(QPointF(*cp_b), QPointF(*p1b))
        item2 = QGraphicsPathItem(path2)
        item2.setPen(QPen(QColor(*self.color), 1.5))
        self.paper.addItem(item2)
        self._items.append(item2)
        
        # Head 2
        tangent2 = atan2(p1b[1] - cp_b[1], p1b[0] - cp_b[0])
        h2 = (p1b[0] - head_len * cos(tangent2 - pi/6), p1b[1] - head_len * sin(tangent2 - pi/6))
        self._items.append(self.paper.addLine(p1b + h2, color=self.color, width=1.5))
        if heads == "full":
            h2b = (p1b[0] - head_len * cos(tangent2 + pi/6), p1b[1] - head_len * sin(tangent2 + pi/6))
            self._items.append(self.paper.addLine(p1b + h2b, color=self.color, width=1.5))

    def clear_drawings(self):
        if not self.paper: return
        for item in self._items:
            try:
                if hasattr(self.paper, 'removeFocusable'):
                    self.paper.removeFocusable(item)
                self.paper.removeItem(item)
            except (RuntimeError, AttributeError):
                pass
        self._items = []
        if hasattr(self, '_focus_item') and self._focus_item:
            try: self.paper.removeItem(self._focus_item)
            except: pass
            self._focus_item = None
        if hasattr(self, '_selection_item') and self._selection_item:
            try: self.paper.removeItem(self._selection_item)
            except: pass
            self._selection_item = None

    def set_focus(self, focus):
        if not self.paper: return
        if focus:
            bbox = self.bounding_box()
            from .app_data import Settings
            self._focus_item = self.paper.addRect(bbox, color=(200, 200, 255), width=1)
            self.paper.toSelectionLayer(self._focus_item)
        else:
            if hasattr(self, '_focus_item') and self._focus_item:
                try: self.paper.removeItem(self._focus_item)
                except: pass
                self._focus_item = None

    def set_selected(self, select):
        if not self.paper: return
        if select:
            bbox = self.bounding_box()
            from .app_data import Settings
            self._selection_item = self.paper.addRect(bbox, color=Settings.selection_color, width=1)
            self.paper.toSelectionLayer(self._selection_item)
        else:
            if hasattr(self, '_selection_item') and self._selection_item:
                try: self.paper.removeItem(self._selection_item)
                except: pass
                self._selection_item = None

    def bounding_box(self):
        x1, y1 = self.p1
        x2, y2 = self.p2
        return [min(x1, x2)-10, min(y1, y2)-10, max(x1, x2)+10, max(y1, y2)+10]

    def get_center(self):
        return (self.p1[0] + self.p2[0]) / 2, (self.p1[1] + self.p2[1]) / 2

    def move_by(self, dx, dy):
        self.p1 = (self.p1[0] + dx, self.p1[1] + dy)
        self.p2 = (self.p2[0] + dx, self.p2[1] + dy)
        self.draw()

    def flip_horizontal(self, center_x):
        self.p1 = (2 * center_x - self.p1[0], self.p1[1])
        self.p2 = (2 * center_x - self.p2[0], self.p2[1])
        self.draw()

    def flip_vertical(self, center_y):
        self.p1 = (self.p1[0], 2 * center_y - self.p1[1])
        self.p2 = (self.p2[0], 2 * center_y - self.p2[1])
        self.curvature *= -1
        self.draw()

    def rotate(self, angle, center):
        from math import sin, cos
        cx, cy = center
        s, c = sin(angle), cos(angle)
        
        # Rotate p1
        x, y = self.p1[0] - cx, self.p1[1] - cy
        self.p1 = (x * c - y * s + cx, x * s + y * c + cy)
        
        # Rotate p2
        x, y = self.p2[0] - cx, self.p2[1] - cy
        self.p2 = (x * c - y * s + cx, x * s + y * c + cy)
        
        self.draw()

    def scale(self, factor, center=None):
        if center is None: center = self.get_center()
        cx, cy = center
        self.p1 = (cx + (self.p1[0] - cx) * factor, cy + (self.p1[1] - cy) * factor)
        self.p2 = (cx + (self.p2[0] - cx) * factor, cy + (self.p2[1] - cy) * factor)
        self.draw()

    def clone(self):
        new_arrow = Arrow(self.p1, self.p2, type=self.type, curvature=self.curvature)
        new_arrow.color = self.color
        return new_arrow
