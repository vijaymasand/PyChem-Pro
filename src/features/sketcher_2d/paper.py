# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3)
from .app_data import App, Settings
from .drawing_parents import Color, Font, Align, PenStyle, LineCap
from .common import bbox_of_bboxes
from .tool_helpers import draw_objs_recursively
from .undo_manager import UndoManager

from src.shared.qt_compat import (QGraphicsScene, QGraphicsItem, QGraphicsTextItem,
    QRectF, QPointF, Qt, QColor, QPen, QBrush, QPolygonF, QPainterPath,
    QFontMetricsF, QFont, QTransform, Signal)

class Paper(QGraphicsScene):
    text_editing_finished = Signal()

    def __init__(self, view=None):
        QGraphicsScene.__init__(self, view)
        self.view = view
        if view:
            view.setScene(self)
        self.paper_item = None
        self.objects = []
        self.dirty_objects = set()
        self.mouse_pressed = False
        self.dragging = False
        self.modifier_keys = set()
        self.focusable_items = set()
        self.do_not_focus = set()
        self.focused_obj = None
        self.locked_focus_obj = None
        self.selected_objs = []
        item = self.addText("")
        self.textitem_margin = item.boundingRect().width() / 2
        self.removeItem(item)
        self.undo_manager = UndoManager(self)

    def setSize(self, w, h):
        self.setSceneRect(0, 0, w, h)
        if self.paper_item:
            self.removeItem(self.paper_item)
        self.paper_item = self.addRect([0, 0, w, h], style=PenStyle.no_line, fill=(255, 255, 255))
        self.paper_item.setZValue(-10)

    def addObject(self, obj):
        if obj not in self.objects:
            self.objects.append(obj)
            obj.paper = self

    def removeObject(self, obj):
        if obj in self.objects:
            obj.paper = None
            self.objects.remove(obj)

    def objectsInRect(self, rect):
        x1, y1, x2, y2 = rect
        gfx_items = self.items(QRectF(x1, y1, x2 - x1, y2 - y1))
        # self.items() returns items in Z-order descending (topmost first)
        result = []
        seen = set()
        for itm in gfx_items:
            if itm in self.focusable_items:
                obj = itm.object
                if obj not in seen:
                    result.append(obj)
                    seen.add(obj)
        return result

    def redraw_dirty_objects(self):
        draw_objs_recursively(self.dirty_objects)
        self.dirty_objects.clear()

    def addLine(self, line, width=1, color=Color.black, style=PenStyle.solid, cap=LineCap.square):
        pen = QPen(QColor(*color), width, Qt.PenStyle(style), Qt.PenCapStyle(cap))
        return QGraphicsScene.addLine(self, *line, pen)

    def addRect(self, rect, width=1, color=Color.black, style=PenStyle.solid, fill=None):
        x1, y1, x2, y2 = rect
        pen = QPen(QColor(*color), width, Qt.PenStyle(style))
        brush = QBrush(QColor(*fill)) if fill else QBrush()
        return QGraphicsScene.addRect(self, x1, y1, x2 - x1, y2 - y1, pen, brush)

    def addPolygon(self, points, width=1, color=Color.black, style=PenStyle.solid, fill=None):
        polygon = QPolygonF([QPointF(*p) for p in points])
        pen = QPen(QColor(*color), width, Qt.PenStyle(style))
        brush = QBrush(QColor(*fill)) if fill else QBrush()
        return QGraphicsScene.addPolygon(self, polygon, pen, brush)

    def addPolyline(self, points, width=1, color=Color.black, style=PenStyle.solid):
        shape = QPainterPath(QPointF(*points[0]))
        for pt in points[1:]: shape.lineTo(QPointF(*pt))
        pen = QPen(QColor(*color), width, Qt.PenStyle(style))
        return QGraphicsScene.addPath(self, shape, pen)

    def addEllipse(self, rect, width=1, color=Color.black, fill=None):
        x1, y1, x2, y2 = rect
        pen = QPen(QColor(*color), width)
        brush = QBrush(QColor(*fill)) if fill else QBrush()
        return QGraphicsScene.addEllipse(self, x1, y1, x2 - x1, y2 - y1, pen, brush)

    def addHtmlText(self, text, pos, font=None, align=Align.Left | Align.Baseline, color=(0, 0, 0)):
        item = QGraphicsTextItem()
        item.setDefaultTextColor(QColor(*color))
        if font:
            qf = QFont(font.name)
            sz = int(round(font.size))
            if sz > 0:
                qf.setPixelSize(sz)
            item.setFont(qf)
        item.setHtml(text)
        self.addItem(item)
        fm = QFontMetricsF(item.font())
        item_w = item.boundingRect().width()
        x, y = pos[0] - self.textitem_margin, pos[1] - self.textitem_margin
        w, h = item_w - 2 * self.textitem_margin, fm.height()
        if align & Align.HCenter: x -= w / 2
        elif align & Align.Right: x -= w
        if align & Align.Baseline: y -= fm.ascent()
        elif align & Align.Bottom: y -= h
        elif align & Align.VCenter: y -= h / 2
        item.setPos(x, y)
        return item

    def addChemicalFormula(self, text, pos, align, offset, font, color=(0, 0, 0)):
        item = QGraphicsTextItem()
        item.setDefaultTextColor(QColor(*color))
        if font:
            qf = QFont(font.name)
            sz = int(round(font.size))
            if sz > 0:
                qf.setPixelSize(sz)
            item.setFont(qf)
        item.setHtml(text)
        self.addItem(item)
        rect = item.boundingRect()
        w, h = rect.width(), rect.height()
        x, y, w = pos[0] - self.textitem_margin, pos[1] - h / 2, w - 2 * self.textitem_margin
        if align == Align.HCenter: x -= w / 2
        elif align == Align.Left: x -= offset
        elif align == Align.Right: x -= w - offset
        item.setPos(x, y)
        return item

    def itemBoundingRect(self, item):
        if item is None: return 0, 0, 0, 0
        rect = item.sceneBoundingRect()
        x, y, w, h = rect.x(), rect.y(), rect.width(), rect.height()
        if isinstance(item, QGraphicsTextItem):
            return x + self.textitem_margin, y + self.textitem_margin, w - 2 * self.textitem_margin, h - 2 * self.textitem_margin
        return x, y, w, h

    def itemBoundingBox(self, item):
        if item is None: return [0, 0, 0, 0]
        rect = item.sceneBoundingRect()
        x1, y1, x2, y2 = rect.left(), rect.top(), rect.right(), rect.bottom()
        if isinstance(item, QGraphicsTextItem):
            return [x1 + self.textitem_margin, y1 + self.textitem_margin, x2 - self.textitem_margin, y2 - self.textitem_margin]
        return [x1, y1, x2, y2]

    def allObjectsBoundingBox(self):
        bboxes = []
        for obj in self.objects:
            bboxes.append(obj.bounding_box())
        return bbox_of_bboxes(bboxes)

    def toBondLayer(self, item): item.setZValue(-2)
    def toSelectionLayer(self, item): item.setZValue(-3)

    def moveItemsBy(self, items, dx, dy):
        for item in items: item.moveBy(dx, dy)

    def getCharWidth(self, char, font):
        qf = QFont(font.name)
        sz = int(round(font.size))
        if sz > 0:
            qf.setPixelSize(sz)
        return QFontMetricsF(qf).horizontalAdvance(char)

    def addFocusable(self, graphics_item, obj):
        graphics_item.object = obj
        self.focusable_items.add(graphics_item)

    def removeFocusable(self, graphics_item):
        if graphics_item in self.focusable_items:
            self.focusable_items.remove(graphics_item)
            graphics_item.object = None

    def changeFocusTo(self, focused_obj):
        if self.focused_obj is focused_obj: return
        if self.focused_obj: self.focused_obj.set_focus(False)
        if focused_obj: focused_obj.set_focus(True)
        self.focused_obj = focused_obj

    def unfocusObject(self, obj):
        if obj is self.focused_obj:
            obj.set_focus(False)
            self.focused_obj = None

    def selectObject(self, obj):
        if obj not in self.selected_objs:
            obj.set_selected(True)
            self.selected_objs.append(obj)

    def deselectObject(self, obj):
        if obj in self.selected_objs:
            obj.set_selected(False)
            self.selected_objs.remove(obj)

    def deselectAll(self):
        for obj in self.selected_objs: obj.set_selected(False)
        self.selected_objs = []

    def mousePressEvent(self, ev):
        if ev.button() not in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton): return
        self.mouse_pressed = True
        self.mouse_button = ev.button()
        self.dragging = False
        self._mouse_press_pos = (ev.scenePos().x(), ev.scenePos().y())
        x, y = self._mouse_press_pos
        
        if ev.button() == Qt.MouseButton.RightButton:
            if App.tool and hasattr(App.tool, 'on_right_click'):
                App.tool.on_right_click(x, y)
            return

        objs = self.objectsInRect([x - 5, y - 5, x + 5, y + 5])
        if objs:
            self.locked_focus_obj = objs[0]
            self.changeFocusTo(objs[0])
        else:
            self.locked_focus_obj = None
            self.changeFocusTo(None)

        if App.tool: App.tool.on_mouse_press(*self._mouse_press_pos)
        QGraphicsScene.mousePressEvent(self, ev)

    def mouseMoveEvent(self, ev):
        x, y = ev.scenePos().x(), ev.scenePos().y()
        if self.mouse_pressed and not self.dragging:
            self.dragging = True
        # Only change focus when not pressing mouse (to preserve focus during typing)
        if not self.mouse_pressed:
            if not self.locked_focus_obj:
                objs = self.objectsInRect([x - 5, y - 5, x + 5, y + 5])
                focused_obj = objs[0] if objs else None
                self.changeFocusTo(focused_obj)
        if App.tool: App.tool.on_mouse_move(x, y)
        QGraphicsScene.mouseMoveEvent(self, ev)

    def mouseReleaseEvent(self, ev):
        if not self.mouse_pressed or getattr(self, 'mouse_button', None) != ev.button():
            return QGraphicsScene.mouseReleaseEvent(self, ev)
            
        self.mouse_pressed = False
        if App.tool and hasattr(App.tool, 'on_mouse_release'): 
            App.tool.on_mouse_release(ev.scenePos().x(), ev.scenePos().y())
        if self.dragging:
            self.locked_focus_obj = None
        self.dragging = False
        QGraphicsScene.mouseReleaseEvent(self, ev)
        
    def mouseDoubleClickEvent(self, ev):
        if ev.button() != Qt.MouseButton.LeftButton: return
        x, y = ev.scenePos().x(), ev.scenePos().y()
        if App.tool and hasattr(App.tool, 'on_mouse_double_click'):
            App.tool.on_mouse_double_click(x, y)
        QGraphicsScene.mouseDoubleClickEvent(self, ev)

    def save_state_to_undo_stack(self, name=''):
        self.undo_manager.save_current_state(name)

    def undo(self):
        self.changeFocusTo(None)
        self.focusable_items.clear()
        self.undo_manager.undo()

    def redo(self):
        self.changeFocusTo(None)
        self.focusable_items.clear()
        self.undo_manager.redo()
