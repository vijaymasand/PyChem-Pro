# -*- coding: utf-8 -*-
from .drawing_parents import DrawableObject, Color, PenStyle, LineCap
from .app_data import Settings
from . import common

class Shape(DrawableObject):
    def __init__(self, x1, y1, x2, y2):
        super().__init__()
        self.x1, self.y1, self.x2, self.y2 = x1, y1, x2, y2
        self.line_width = 1.2
        self.fill_color = None # Transparent by default
        self._main_items = []
        self._selection_item = None
        self._focus_item = None
        self.meta__undo_properties = ("x1", "y1", "x2", "y2", "color", "fill_color", "line_width")

    def move_by(self, dx, dy):
        self.x1 += dx
        self.y1 += dy
        self.x2 += dx
        self.y2 += dy
        self.draw()

    def bounding_box(self):
        return [min(self.x1, self.x2), min(self.y1, self.y2), max(self.x1, self.x2), max(self.y1, self.y2)]

    def clear_drawings(self):
        if not self.paper: return
        for item in self._main_items:
            self.paper.removeFocusable(item)
            self.paper.removeItem(item)
        self._main_items = []
        if self._focus_item:
            self.paper.removeItem(self._focus_item)
            self._focus_item = None
        if self._selection_item:
            self.paper.removeItem(self._selection_item)
            self._selection_item = None

    def set_focus(self, focus):
        if not self.paper: return
        if focus:
            bbox = self.bounding_box()
            self._focus_item = self.paper.addRect([bbox[0]-2, bbox[1]-2, bbox[2]+2, bbox[3]+2], color=Settings.focus_color, width=1)
        else:
            if self._focus_item:
                self.paper.removeItem(self._focus_item)
                self._focus_item = None

    def set_selected(self, select):
        if not self.paper: return
        if self._selection_item:
            self.paper.removeItem(self._selection_item)
            self._selection_item = None
        if select:
            bbox = self.bounding_box()
            self._selection_item = self.paper.addRect([bbox[0]-1, bbox[1]-1, bbox[2]+1, bbox[3]+1], color=Settings.selection_color, width=2)
        else:
            if self._selection_item:
                self.paper.removeItem(self._selection_item)
                self._selection_item = None

class Rectangle(Shape):
    def clone(self):
        c = Rectangle(self.x1, self.y1, self.x2, self.y2)
        c.color = self.color
        c.fill_color = self.fill_color
        c.line_width = self.line_width
        return c

    def draw(self):
        self.clear_drawings()
        if not self.paper: return
        rect = [min(self.x1, self.x2), min(self.y1, self.y2), max(self.x1, self.x2), max(self.y1, self.y2)]
        item = self.paper.addRect(rect, color=self.color, width=self.line_width, fill=self.fill_color)
        self._main_items = [item]
        self.paper.addFocusable(item, self)

class Ellipse(Shape):
    def clone(self):
        c = Ellipse(self.x1, self.y1, self.x2, self.y2)
        c.color = self.color
        c.fill_color = self.fill_color
        c.line_width = self.line_width
        return c

    def draw(self):
        self.clear_drawings()
        if not self.paper: return
        rect = [min(self.x1, self.x2), min(self.y1, self.y2), max(self.x1, self.x2), max(self.y1, self.y2)]
        item = self.paper.addEllipse(rect, color=self.color, width=self.line_width, fill=self.fill_color)
        self._main_items = [item]
        self.paper.addFocusable(item, self)
