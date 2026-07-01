"""
Object Panel — PyMOL-style list of loaded structures.

Sits at the bottom of the left control column.  Each loaded structure
(:class:`~src.app.molecule_scene.MoleculeObject`) gets a row with:

* a **visibility** checkbox (show/hide it in the 3D overlay),
* its **name** (click to make it the active object — chemistry, 2D and export
  target the active object),
* a small **menu** button for per-object representation, colour and removal.

A master **all** row on top toggles every object at once.  The panel is purely
a view: it emits signals and the controller mutates the
:class:`~src.app.molecule_scene.MoleculeScene`, then calls :meth:`set_objects`
to refresh.
"""

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QCheckBox, QPushButton,
    QScrollArea, QFrame, QMenu, QColorDialog, QSizePolicy,
    Qt, Signal, QColor,
)
from src.shared.ui.theme import COLORS, theme_signals


_REPRESENTATIONS = [
    ("Default", None),
    ("Ball & Stick", "ball_and_stick"),
    ("Space Fill", "spacefill"),
    ("Wireframe", "wireframe"),
]


class _ObjectRow(QWidget):
    """A single object row: [✓] name … menu."""

    def __init__(self, obj, panel):
        super().__init__()
        self.obj_id = obj.id
        self._panel = panel

        row = QHBoxLayout(self)
        row.setContentsMargins(2, 1, 2, 1)
        row.setSpacing(4)

        self.check = QCheckBox()
        self.check.setChecked(obj.visible)
        self.check.setToolTip("Show / hide this object")
        self.check.toggled.connect(
            lambda v: panel.visibility_toggled.emit(self.obj_id, v))
        row.addWidget(self.check)

        self.name_btn = QPushButton(obj.name)
        self.name_btn.setFlat(True)
        self.name_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.name_btn.setToolTip(
            f"{obj.name} — {obj.num_atoms} atoms\nClick to make active")
        self.name_btn.setSizePolicy(QSizePolicy.Policy.Expanding,
                                    QSizePolicy.Policy.Fixed)
        self.name_btn.clicked.connect(
            lambda: panel.active_changed.emit(self.obj_id))
        row.addWidget(self.name_btn, 1)

        self.menu_btn = QPushButton("⋯")  # ⋯
        self.menu_btn.setFixedWidth(24)
        self.menu_btn.setToolTip("Object options")
        self.menu_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.menu_btn.clicked.connect(self._show_menu)
        row.addWidget(self.menu_btn)

    def set_active(self, active):
        """Highlight this row when it is the active object."""
        if active:
            self.name_btn.setStyleSheet(
                f"text-align:left; font-weight:600;"
                f"color:{COLORS['accent']};")
        else:
            self.name_btn.setStyleSheet(
                f"text-align:left; color:{COLORS['text_primary']};")

    def _show_menu(self):
        panel = self._panel
        oid = self.obj_id
        menu = QMenu(self)

        rep_menu = menu.addMenu("Representation")
        for label, rep in _REPRESENTATIONS:
            act = rep_menu.addAction(label)
            act.triggered.connect(
                lambda checked=False, r=rep: panel.representation_changed.emit(oid, r))

        color_menu = menu.addMenu("Color")
        act_cpk = color_menu.addAction("By Element (CPK)")
        act_cpk.triggered.connect(
            lambda checked=False: panel.color_changed.emit(oid, None))
        act_auto = color_menu.addAction("Auto (distinct)")
        act_auto.triggered.connect(
            lambda checked=False: panel.color_changed.emit(oid, 'auto'))
        act_custom = color_menu.addAction("Custom…")
        act_custom.triggered.connect(lambda checked=False: self._pick_color())

        menu.addSeparator()
        act_remove = menu.addAction("Remove")
        act_remove.triggered.connect(
            lambda checked=False: panel.remove_requested.emit(oid))

        self.menu_btn.menu_ref = menu  # keep alive
        menu.exec(self.menu_btn.mapToGlobal(
            self.menu_btn.rect().bottomLeft()))

    def _pick_color(self):
        color = QColorDialog.getColor(QColor(*COLORS_accent_rgb()), self,
                                      "Choose Object Color")
        if color.isValid():
            self._panel.color_changed.emit(
                self.obj_id, (color.red(), color.green(), color.blue()))


def COLORS_accent_rgb():
    c = QColor(COLORS['accent'])
    return (c.red(), c.green(), c.blue())


class ObjectPanel(QWidget):
    """Scrollable list of loaded objects with show/hide + per-object options."""

    visibility_toggled = Signal(int, bool)
    all_visibility_toggled = Signal(bool)
    active_changed = Signal(int)
    representation_changed = Signal(int, object)
    color_changed = Signal(int, object)
    remove_requested = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("objectPanel")
        self._rows = []
        self._active_id = -1
        self._init_ui()
        self._apply_theme()
        try:
            theme_signals().theme_changed.connect(self._apply_theme)
        except Exception:
            pass

    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(10, 8, 10, 8)
        outer.setSpacing(6)

        header = QLabel("Objects")
        header.setObjectName("labelSection")
        outer.addWidget(header)

        # Master "all" row.
        all_row = QHBoxLayout()
        all_row.setContentsMargins(2, 0, 2, 0)
        all_row.setSpacing(4)
        self.all_check = QCheckBox()
        self.all_check.setChecked(True)
        self.all_check.setToolTip("Show / hide all objects")
        self.all_check.toggled.connect(self.all_visibility_toggled.emit)
        all_row.addWidget(self.all_check)
        all_label = QLabel("all")
        all_label.setObjectName("labelMuted")
        all_row.addWidget(all_label, 1)
        outer.addLayout(all_row)

        # Scrollable list of object rows.
        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        self._list_container = QWidget()
        self._list_layout = QVBoxLayout(self._list_container)
        self._list_layout.setContentsMargins(0, 0, 0, 0)
        self._list_layout.setSpacing(1)
        self._list_layout.addStretch()
        self._scroll.setWidget(self._list_container)
        outer.addWidget(self._scroll, 1)

        self._empty_label = QLabel("No structures loaded")
        self._empty_label.setObjectName("labelHint")
        self._empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        outer.addWidget(self._empty_label)

    def set_objects(self, objects, active_id):
        """Rebuild the row list from the scene's objects."""
        # Clear existing rows.
        for row in self._rows:
            row.setParent(None)
            row.deleteLater()
        self._rows = []
        self._active_id = active_id

        for obj in objects:
            row = _ObjectRow(obj, self)
            row.set_active(obj.id == active_id)
            # Insert before the trailing stretch.
            self._list_layout.insertWidget(self._list_layout.count() - 1, row)
            self._rows.append(row)

        has_objs = bool(objects)
        self._empty_label.setVisible(not has_objs)
        self._scroll.setVisible(has_objs)

        if has_objs:
            all_visible = all(o.visible for o in objects)
            self.all_check.blockSignals(True)
            self.all_check.setChecked(all_visible)
            self.all_check.blockSignals(False)

    def _apply_theme(self):
        bg = COLORS['bg_secondary']
        self.setStyleSheet(
            f"QWidget#objectPanel {{ background-color: {bg};"
            f" border-top: 1px solid {COLORS['border']}; }}")
        for row in self._rows:
            row.set_active(row.obj_id == self._active_id)
