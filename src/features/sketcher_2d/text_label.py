# -*- coding: utf-8 -*-
import re
from .drawing_parents import DrawableObject, Color, Font, Align
from src.shared.qt_compat import QGraphicsTextItem, QFont, QColor, QGraphicsItem

global text_id_no
text_id_no = 1

def format_text_to_html(text):
    if not text:
        return ""
    
    # 1. Replace explicit LaTeX/markdown-style markup
    # Superscripts: ^{...}, ^(...) or ^c (single character/digit/sign)
    text = re.sub(r'\^\{([^}]+)\}', r'<sup>\1</sup>', text)
    text = re.sub(r'\^\(([^)]+)\)', r'<sup>\1</sup>', text)
    text = re.sub(r'\^([a-zA-Z0-9+\-])', r'<sup>\1</sup>', text)
    
    # Subscripts: _{...}, _(...) or _c (single character/digit/sign)
    text = re.sub(r'\_\{([^}]+)\}', r'<sub>\1</sub>', text)
    text = re.sub(r'\_\(([^)]+)\)', r'<sub>\1</sub>', text)
    text = re.sub(r'\_([a-zA-Z0-9+\-])', r'<sub>\1</sub>', text)
    
    # 2. Split by HTML tags to avoid formatting content inside tags
    parts = re.split(r'(<[^>]+>)', text)
    
    for i in range(len(parts)):
        # Even indices are plain text, odd indices are HTML tags
        if i % 2 == 0 and parts[i]:
            # Split by whitespace to process word by word
            tokens = re.split(r'(\s+)', parts[i])
            for j in range(len(tokens)):
                token = tokens[j]
                if not token or token.isspace():
                    continue
                
                # Strip trailing punctuation for chemical formula check
                match = re.match(r'^(.*?)([\.,;:!\?]*)$', token)
                if match:
                    clean_word, suffix = match.groups()
                else:
                    clean_word, suffix = token, ""
                
                # Format the clean word
                formatted_word = clean_word
                
                # Rule A: Hybridization (e.g. sp3, dsp2, sp3d2)
                # Word consists only of lowercase s, p, d, f and digits, has at least one of s,p,d,f, and has at least one digit
                if re.match(r'^[spdf\d]+$', clean_word) and any(c in 'spdf' for c in clean_word) and any(c.isdigit() for c in clean_word):
                    formatted_word = re.sub(r'(\d+)', r'<sup>\1</sup>', clean_word)
                else:
                    # Rule B: Check for charge at the end (e.g. NH4+, SO42-, Ca2+, H+)
                    # Requires at least one chemical formula character before the charge (non-greedy, max 1 digit for charge magnitude)
                    charge_match = re.match(r'^([a-zA-Z0-9()\[\]\-\*·.]+?)(\d?[+-])$', clean_word)
                    if charge_match:
                        base, charge = charge_match.groups()
                        charge_html = f"<sup>{charge}</sup>"
                    else:
                        base = clean_word
                        charge_html = ""
                    
                    # Rule C: Isotopes (e.g. 13C, 14N)
                    # Leading digits followed by an uppercase letter
                    base = re.sub(r'^(\d+)([A-Z][a-z]?)', r'<sup>\1</sup>\2', base)
                    
                    # Rule D: Subscript digits preceded by a letter, ), or ]
                    base = re.sub(r'([a-zA-Z)\]])(\d+)', r'\1<sub>\2</sub>', base)
                    
                    formatted_word = base + charge_html
                
                tokens[j] = formatted_word + suffix
            
            parts[i] = "".join(tokens)
            
    return "".join(parts)

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
        
        html_text = format_text_to_html(self.text)
        self._text_item = self.paper.addHtmlText(html_text, (self.x, self.y), font=font, 
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
                if not self._text_item.toPlainText().strip():
                    self.text = ""
                else:
                    new_text = self._text_item.toHtml()
                    if new_text != self.text:
                        self.text = new_text
                self.draw()  # Always redraw to restore HTML rendering!
                if self.paper:
                    # the widget drops labels that were left empty
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

