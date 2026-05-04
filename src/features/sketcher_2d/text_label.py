# -*- coding: utf-8 -*-
from .drawing_parents import DrawableObject, Color, Font, Align
from src.shared.qt_compat import QGraphicsTextItem, QFont, QColor, QGraphicsItem

global text_id_no
text_id_no = 1

class TextLabel(DrawableObject):
    meta__undo_properties = ("x", "y", "text", "font_name", "font_size", "color")
    def __init__(self, x, y, text="Text"):
        DrawableObject.__init__(self)
        self.x = x
        self.y = y
        self.text = text
        self.font_name = "Arial"
        self.font_size = 14
        self._text_item = None
        self._bg_item = None
        self._focus_item = None
        self._selection_item = None
        global text_id_no
        self.id = 'text' + str(text_id_no)
        text_id_no += 1

    @property
    def pos(self):
        return self.x, self.y

    def set_pos(self, x, y):
        self.x, self.y = x, y

    def draw(self):
        self.clear_drawings()
        if not self.paper: return
        
        font = Font(self.font_name, self.font_size)
        
        self._text_item = self.paper.addHtmlText(self.text, (self.x, self.y), font=font, 
                                                 align=Align.HCenter | Align.VCenter, color=self.color)
        if hasattr(self._text_item, 'setFlag'):
            self._text_item.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsFocusable, True)
            
        self.paper.addFocusable(self._text_item, self)

    def clear_drawings(self):
        if not self.paper: return
        if self._text_item:
            try:
                self.paper.removeFocusable(self._text_item)
                self.paper.removeItem(self._text_item)
            except:
                pass
            self._text_item = None
        if self._bg_item:
            try:
                self.paper.removeItem(self._bg_item)
            except:
                pass
            self._bg_item = None
        if self._focus_item:
            try:
                self.paper.removeItem(self._focus_item)
            except:
                pass
            self._focus_item = None
        if self._selection_item:
            try:
                self.paper.removeItem(self._selection_item)
            except:
                pass
            self._selection_item = None

    def set_focus(self, focus):
        if not self.paper: return
        
        # Import Qt here to use for flags
        from src.shared.qt_compat import Qt
        
        if focus:
            if hasattr(self._text_item, 'setTextInteractionFlags'):
                self._text_item.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
                self._text_item.setFocus()
            bbox = self.bounding_box()
            from .app_data import Settings
            self._focus_item = self.paper.addRect(bbox, color=(200, 200, 255), width=1)
            self.paper.toSelectionLayer(self._focus_item)
        else:
            if hasattr(self._text_item, 'setTextInteractionFlags'):
                self._text_item.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
                self._text_item.clearFocus()
                new_text = self._text_item.toPlainText()
                if new_text != self.text:
                    self.text = new_text
                    self.draw()
                if self.paper:
                    self.paper.text_editing_finished.emit()
            if self._focus_item:
                try:
                    self.paper.removeItem(self._focus_item)
                except:
                    pass
                self._focus_item = None

    def set_selected(self, select):
        if not self.paper: return
        if select:
            bbox = self.bounding_box()
            from .app_data import Settings
            self._selection_item = self.paper.addRect(bbox, color=Settings.selection_color, width=1)
            self.paper.toSelectionLayer(self._selection_item)
        else:
            if self._selection_item:
                try:
                    self.paper.removeItem(self._selection_item)
                except:
                    pass
                self._selection_item = None

    def bounding_box(self):
        if self._text_item and self.paper:
            return self.paper.itemBoundingBox(self._text_item)
        return [self.x - 20, self.y - 10, self.x + 20, self.y + 10]

    def get_center(self):
        return self.x, self.y

    def move_by(self, dx, dy):
        self.x += dx
        self.y += dy
        self.draw()

    def flip_horizontal(self, center_x):
        self.x = 2 * center_x - self.x
        self.draw()

    def flip_vertical(self, center_y):
        self.y = 2 * center_y - self.y
        self.draw()

    def rotate(self, angle, center):
        # Text doesn't rotate in this implementation
        pass

    def scale(self, factor, center=None):
        if center is None: center = self.get_center()
        cx, cy = center
        self.x = cx + (self.x - cx) * factor
        self.y = cy + (self.y - cy) * factor
        self.font_size *= factor
        self.draw()

    def clone(self):
        new_text = TextLabel(self.x, self.y, text=self.text)
        new_text.color = self.color
        new_text.font_name = self.font_name
        new_text.font_size = self.font_size
        return new_text

