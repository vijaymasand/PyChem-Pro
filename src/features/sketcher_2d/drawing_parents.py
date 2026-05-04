# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3)

class Color:
    black = (0, 0, 0)
    white = (255, 255, 255)
    darkGray = (128, 128, 128)
    gray = (160, 160, 160)
    lightGray = (192, 192, 192)
    red = (255, 0, 0)
    green = (0, 255, 0)
    blue = (0, 0, 255)
    cyan = (0, 255, 255)
    magenta = (255, 0, 255)
    yellow = (255, 255, 0)
    transparent = (0, 0, 0, 0)

def hex_color(color):
    clr = '#'
    for x in color[:3]:
        clr += "%.2x" % x
    return clr

class PenStyle:
    no_line = 0
    solid = 1
    dashed = 2
    dotted = 3

class LineCap:
    butt = 0x00
    square = 0x10
    round = 0x20

class Font:
    def __init__(self, name="Sans Serif", size=10):
        self.name = name
        self.size = size
        self.bold = False
        self.italic = False

class DrawableObject:
    focus_priority = 10
    redraw_priority = 10
    is_toplevel = True
    meta__undo_properties = ()
    meta__undo_copy = ()
    meta__undo_children_to_record = ()
    meta__same_objects = {}
    meta__scalables = ()

    def __init__(self):
        self.paper = None
        self.color = (0, 0, 0)

    @property
    def class_name(self):
        return self.__class__.__name__

    @property
    def parent(self):
        return None

    @property
    def children(self):
        return []

    def draw(self):
        pass

    def clear_drawings(self):
        pass

    def bounding_box(self):
        return None

    def scale(self, factor, center=None):
        pass

    def rotate(self, angle, center=None):
        pass

    def transform(self, tr):
        pass

    def move_by(self, dx, dy):
        pass

    def set_focus(self, focus):
        pass

    def set_selected(self, selected):
        pass

    def mark_dirty(self):
        if self.paper:
            self.paper.dirty_objects.add(self)

    def delete_from_paper(self):
        if not self.paper:
            return
        if self in self.paper.dirty_objects:
            self.paper.dirty_objects.remove(self)
        self.paper.unfocusObject(self)
        self.paper.deselectObject(self)
        self.clear_drawings()
        if self.is_toplevel:
            self.paper.removeObject(self)

class Align:
    Left = 0x01
    Right = 0x02
    HCenter = 0x04
    Top = 0x20
    Bottom = 0x40
    VCenter = 0x80
    Baseline = 0x100
