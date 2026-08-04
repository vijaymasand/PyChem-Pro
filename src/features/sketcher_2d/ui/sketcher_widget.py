# -*- coding: utf-8 -*-
from src.shared.qt_compat import *
from src.shared.qt_compat import QImage, QPainter, QFileDialog, QRectF

from ..paper import Paper
from ..tools import (StructureTool, EraserTool, ShapeTool, toolsettings,
                     ring_templates)
from ..fragments import known_labels
from ..shapes import Rectangle, Ellipse
from ..app_data import App, Settings, common_elements, periodic_table

# order the templates the way a chemist looks for them
ring_names = ["benzene", "cyclohexane", "cyclopentane", "cyclopropane",
              "cyclobutane", "cycloheptane", "cyclooctane", "naphthalene",
              "pyridine", "pyrrole", "furan", "thiophene", "cyclopentadiene"]
ring_names += [n for n in ring_templates if n not in ring_names]


class SketcherWidget(QWidget):
    molecule_imported = Signal(str) # SMILES string

    def __init__(self, parent=None):
        super().__init__(parent)
        self._init_ui()
        App.paper = self.paper # Shared App object

    def _add_tool_action(self, action, shortcut=None):
        """ registers a checkable tool action on the toolbar and its shortcut """
        if shortcut:
            action.setShortcut(QKeySequence(shortcut))
            action.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            self.addAction(action)
        self.tool_group.addAction(action)
        self.toolbar.addAction(action)
        return action

    def _init_ui(self):
        # Main layout: Horizontal (Toolbar on left, Canvas on right)
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(5, 5, 5, 5)
        main_layout.setSpacing(5)

        # Toolbar (Vertical)
        self.toolbar = QToolBar()
        self.toolbar.setOrientation(Qt.Orientation.Vertical)
        self.toolbar.setIconSize(QSize(24, 24))
        self.toolbar.setMovable(False)
        self.toolbar.setStyleSheet("""
            QToolBar {
                background-color: transparent;
                border-right: 1px solid palette(mid);
                padding: 5px;
            }
            QToolButton {
                margin: 2px;
                padding: 5px;
                border-radius: 4px;
            }
            QToolButton:checked {
                background-color: palette(highlight);
                color: palette(highlighted-text);
            }
        """)
        main_layout.addWidget(self.toolbar)

        # Tool group for exclusive selection
        self.tool_group = QActionGroup(self)

        # Basic Tools. The single key shortcuts avoid c/n/o/s/p/f/l/b/i/h,
        # those retype the atom under the cursor.
        self.action_select = QAction("Select", self)
        self.action_select.setToolTip("Select & Move (V)")
        self.action_select.setCheckable(True)
        self.action_select.setChecked(True)
        self._add_tool_action(self.action_select, "V")

        self.action_bond = QAction("Bond", self)
        self.action_bond.setToolTip("Draw bonds (D)\n"
                                    "Drag from an atom, bonds snap to 15° steps\n"
                                    "Hold Shift while dragging for a free bond")
        self.action_bond.setCheckable(True)
        self._add_tool_action(self.action_bond, "D")

        self.action_eraser = QAction("Eraser", self)
        self.action_eraser.setToolTip("Remove atoms or bonds, drag to keep erasing (E)")
        self.action_eraser.setCheckable(True)
        self._add_tool_action(self.action_eraser, "E")

        self.action_text = QAction("Text", self)
        self.action_text.setToolTip("Add text (T)")
        self.action_text.setCheckable(True)
        self._add_tool_action(self.action_text, "T")

        self.action_rotate = QAction("Rotate", self)
        self.action_rotate.setToolTip("Rotate molecule (R)")
        self.action_rotate.setCheckable(True)
        self._add_tool_action(self.action_rotate, "R")

        self.toolbar.addSeparator()

        # Bond Types
        self.toolbar.addWidget(QLabel("Bond"))
        self.bond_combo = QComboBox()
        self.bond_combo.addItems(["single", "double", "triple", "wedge", "hashed_wedge"])
        self.bond_combo.setToolTip("Type used for new bonds.\n"
                                   "Clicking an existing bond applies it, "
                                   "clicking again cycles the order / flips the wedge.")
        self.bond_combo.currentTextChanged.connect(self._on_bond_type_changed)
        self.toolbar.addWidget(self.bond_combo)

        self.toolbar.addSeparator()

        # Atom / group label
        self.toolbar.addWidget(QLabel("Atom"))
        self.element_combo = QComboBox()
        self.element_combo.setEditable(True)
        self.element_combo.addItems(common_elements)
        self.element_combo.addItems([g for g in known_labels() if g not in common_elements])
        self.element_combo.setCurrentText("C")
        self.element_combo.setToolTip(
            "Element or group used when drawing.\n"
            "Group labels (OH, Ph, COOH, NO2, Boc ...) are expanded into atoms.\n"
            "You can also hover an atom and press c/n/o/s/p/f/l/b/i/h.")
        self.element_combo.currentTextChanged.connect(self._on_element_changed)
        self.toolbar.addWidget(self.element_combo)

        self.action_charge_plus = QAction("Charge +", self)
        self.action_charge_plus.setToolTip("Click an atom to raise its charge")
        self.action_charge_plus.setCheckable(True)
        self._add_tool_action(self.action_charge_plus)

        self.action_charge_minus = QAction("Charge −", self)
        self.action_charge_minus.setToolTip("Click an atom to lower its charge")
        self.action_charge_minus.setCheckable(True)
        self._add_tool_action(self.action_charge_minus)

        self.action_lone_pair = QAction("Lone Pair", self)
        self.action_lone_pair.setToolTip("Click an atom to add/remove lone pairs")
        self.action_lone_pair.setCheckable(True)
        self._add_tool_action(self.action_lone_pair)

        self.action_radical_plus = QAction("Radical +", self)
        self.action_radical_plus.setToolTip("Click an atom to toggle radical cation (dot with plus)")
        self.action_radical_plus.setCheckable(True)
        self._add_tool_action(self.action_radical_plus)

        self.action_radical_minus = QAction("Radical −", self)
        self.action_radical_minus.setToolTip("Click an atom to toggle radical anion (dot with minus)")
        self.action_radical_minus.setCheckable(True)
        self._add_tool_action(self.action_radical_minus)

        self.toolbar.addSeparator()

        # Rings
        self.toolbar.addWidget(QLabel("Rings"))
        self.ring_combo = QComboBox()
        self.ring_combo.addItems(ring_names)
        self.ring_combo.setCurrentText("benzene")
        self.ring_combo.setToolTip("Ring template.\n"
                                   "Click empty paper for a free ring (drag to spin it),\n"
                                   "an atom to hang the ring on it, "
                                   "or a bond to fuse the ring onto it.")
        self.ring_combo.currentTextChanged.connect(self._on_ring_changed)
        self.toolbar.addWidget(self.ring_combo)

        self.action_ring = QAction("Ring", self)
        self.action_ring.setToolTip("Place the selected ring (G)")
        self.action_ring.setCheckable(True)
        self._add_tool_action(self.action_ring, "G")

        self.toolbar.addSeparator()
        self.toolbar.addWidget(QLabel("Arrows"))
        self.action_arrow = QAction("Reaction", self)
        self.action_arrow.setCheckable(True)
        self.tool_group.addAction(self.action_arrow)
        self.toolbar.addAction(self.action_arrow)

        self.action_equilibrium = QAction("Equil.", self)
        self.action_equilibrium.setCheckable(True)
        self.tool_group.addAction(self.action_equilibrium)
        self.toolbar.addAction(self.action_equilibrium)

        self.action_reversible = QAction("Revers.", self)
        self.action_reversible.setCheckable(True)
        self.tool_group.addAction(self.action_reversible)
        self.toolbar.addAction(self.action_reversible)

        self.action_curly = QAction("Curly CCW", self)
        self.action_curly.setCheckable(True)
        self.tool_group.addAction(self.action_curly)
        self.toolbar.addAction(self.action_curly)

        self.action_curly_cw = QAction("Curly CW", self)
        self.action_curly_cw.setCheckable(True)
        self.tool_group.addAction(self.action_curly_cw)
        self.toolbar.addAction(self.action_curly_cw)

        self.action_fish_up = QAction("Fish Up", self)
        self.action_fish_up.setCheckable(True)
        self.tool_group.addAction(self.action_fish_up)
        self.toolbar.addAction(self.action_fish_up)

        self.action_fish_down = QAction("Fish Down", self)
        self.action_fish_down.setCheckable(True)
        self.tool_group.addAction(self.action_fish_down)
        self.toolbar.addAction(self.action_fish_down)

        self.toolbar.addSeparator()

        self.action_rectangle = QAction("Rect", self)
        self.action_rectangle.setCheckable(True)
        self.action_rectangle.setToolTip("Rectangle (Shift for Square)")
        self.toolbar.addAction(self.action_rectangle)
        self.tool_group.addAction(self.action_rectangle)

        self.action_ellipse = QAction("Ellipse", self)
        self.action_ellipse.setCheckable(True)
        self.action_ellipse.setToolTip("Ellipse (Shift for Circle)")
        self.toolbar.addAction(self.action_ellipse)
        self.tool_group.addAction(self.action_ellipse)

        self.toolbar.addSeparator()

        # back to the select tool, the way Esc works in every drawing program
        self.action_escape = QAction("Esc", self)
        self.action_escape.setShortcut(QKeySequence(Qt.Key.Key_Escape))
        self.action_escape.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        self.action_escape.triggered.connect(self._activate_select_tool)
        self.addAction(self.action_escape)

        # Right side: Canvas and Right Toolbar
        self.view = QGraphicsView()
        self.view.setRenderHint(QPainter.Antialiasing)
        self.view.setRenderHint(QPainter.SmoothPixmapTransform)
        self.view.setBackgroundBrush(QBrush(QColor(255, 255, 255)))
        self.view.setStyleSheet("QGraphicsView { border: 1px solid palette(mid); border-radius: 4px; }")
        main_layout.addWidget(self.view, 1)

        # Right Toolbar (Vertical)
        self.right_toolbar = QToolBar()
        self.right_toolbar.setOrientation(Qt.Orientation.Vertical)
        self.right_toolbar.setIconSize(QSize(24, 24))
        self.right_toolbar.setMovable(False)
        self.right_toolbar.setStyleSheet(self.toolbar.styleSheet())
        main_layout.addWidget(self.right_toolbar)

        # Zoom Controls (Right Toolbar)
        self.action_zoom_in = QAction("Zoom In", self)
        self.action_zoom_in.triggered.connect(lambda: self._zoom(1.2))
        self.right_toolbar.addAction(self.action_zoom_in)

        self.action_zoom_out = QAction("Zoom Out", self)
        self.action_zoom_out.triggered.connect(lambda: self._zoom(0.8))
        self.right_toolbar.addAction(self.action_zoom_out)

        self.action_zoom_reset = QAction("Reset", self)
        self.action_zoom_reset.triggered.connect(self._zoom_reset)
        self.right_toolbar.addAction(self.action_zoom_reset)

        self.right_toolbar.addSeparator()

        # Undo/Redo/Clear (Right Toolbar)
        self.action_undo = QAction("Undo", self)
        self.action_undo.setShortcut(QKeySequence.Undo)
        self.action_undo.triggered.connect(self._on_undo)
        self.right_toolbar.addAction(self.action_undo)

        self.action_redo = QAction("Redo", self)
        self.action_redo.setShortcut(QKeySequence.Redo)
        self.action_redo.triggered.connect(self._on_redo)
        self.right_toolbar.addAction(self.action_redo)

        self.right_toolbar.addSeparator()

        self.action_flip_h = QAction("Flip H", self)
        self.action_flip_h.setToolTip("Flip Horizontal")
        self.action_flip_h.triggered.connect(self._on_flip_h)
        self.right_toolbar.addAction(self.action_flip_h)

        self.action_flip_v = QAction("Flip V", self)
        self.action_flip_v.setToolTip("Flip Vertical")
        self.action_flip_v.triggered.connect(self._on_flip_v)
        self.right_toolbar.addAction(self.action_flip_v)

        self.action_cleanup = QAction("Clean Up", self)
        self.action_cleanup.setToolTip("Lay the structure out again with even "
                                       "bond lengths and angles")
        self.action_cleanup.triggered.connect(self._on_cleanup)
        self.right_toolbar.addAction(self.action_cleanup)

        self.right_toolbar.addSeparator()

        # Copy/Paste (Right Toolbar)
        self.action_copy = QAction("Copy", self)
        self.action_copy.setShortcut(QKeySequence.Copy)
        self.action_copy.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        self.action_copy.triggered.connect(self._on_copy)
        self.right_toolbar.addAction(self.action_copy)
        self.addAction(self.action_copy)

        self.action_paste = QAction("Paste", self)
        self.action_paste.setShortcut(QKeySequence.Paste)
        self.action_paste.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        self.action_paste.triggered.connect(self._on_paste)
        self.right_toolbar.addAction(self.action_paste)
        self.addAction(self.action_paste)

        self.right_toolbar.addSeparator()

        self.action_export = QAction("Export", self)
        self.action_export.setToolTip("Export Selected as Image")
        self.action_export.triggered.connect(self._on_export)
        self.right_toolbar.addAction(self.action_export)

        self.action_save_mol = QAction("Save", self)
        self.action_save_mol.setToolTip("Save the structure as MOL / SDF (Ctrl+S)")
        self.action_save_mol.setShortcut(QKeySequence.Save)
        self.action_save_mol.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        self.action_save_mol.triggered.connect(self._on_save_mol)
        self.right_toolbar.addAction(self.action_save_mol)
        self.addAction(self.action_save_mol)

        self.action_open_mol = QAction("Open", self)
        self.action_open_mol.setToolTip("Open a MOL / SDF file (Ctrl+O)")
        self.action_open_mol.setShortcut(QKeySequence.Open)
        self.action_open_mol.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        self.action_open_mol.triggered.connect(self._on_open_mol)
        self.right_toolbar.addAction(self.action_open_mol)
        self.addAction(self.action_open_mol)

        self.right_toolbar.addSeparator()

        self.action_smiles = QAction("SMILES", self)
        self.action_smiles.setToolTip("SMILES to 2D")
        self.action_smiles.triggered.connect(self._on_smiles_to_2d)
        self.right_toolbar.addAction(self.action_smiles)

        self.action_copy_smiles = QAction("Get SMILES", self)
        self.action_copy_smiles.setToolTip("Copy the SMILES of the structure to the clipboard")
        self.action_copy_smiles.triggered.connect(self._on_copy_smiles)
        self.right_toolbar.addAction(self.action_copy_smiles)

        self.right_toolbar.addSeparator()

        self.action_clear = QAction("Clear", self)
        self.action_clear.setToolTip("Remove everything from the page")
        self.action_clear.triggered.connect(self._on_clear)
        self.right_toolbar.addAction(self.action_clear)

        self.right_toolbar.addSeparator()

        self.action_select_all = QAction("Select All", self)
        self.action_select_all.setShortcut(QKeySequence.SelectAll)
        self.action_select_all.triggered.connect(self._on_select_all)
        self.addAction(self.action_select_all)
        
        self.action_page_setup = QAction("Page", self)
        self.action_page_setup.setToolTip("Page Setup")
        self.action_page_setup.triggered.connect(self._on_page_setup)
        self.right_toolbar.addAction(self.action_page_setup)
        
        self.right_toolbar.addSeparator()
        
        self.action_color = QAction("Color", self)
        self.action_color.setToolTip("Change Color of Selected")
        self.action_color.triggered.connect(self._on_color)
        self.right_toolbar.addAction(self.action_color)

        self.right_toolbar.addSeparator()
        self.right_toolbar.addWidget(QLabel("Carbons"))
        self.carbon_combo = QComboBox()
        self.carbon_combo.addItems(["Never", "Terminal", "All"])
        self.carbon_combo.setToolTip("Show the C label on carbon atoms")
        self.carbon_combo.currentTextChanged.connect(self._on_show_carbon_changed)
        self.right_toolbar.addWidget(self.carbon_combo)

        # Spacer for Right Toolbar
        right_spacer = QWidget()
        right_spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.right_toolbar.addWidget(right_spacer)
        
        self.action_import = QAction("Import", self)
        self.action_import.setToolTip("Import to 3D View (Ctrl+I)")
        self.action_import.setShortcut(QKeySequence("Ctrl+I"))
        self.action_import.triggered.connect(self._on_import)
        self.right_toolbar.addAction(self.action_import)
        
        # Style the button specifically to make it visible
        import_btn = self.right_toolbar.widgetForAction(self.action_import)
        if import_btn:
            import_btn.setStyleSheet("font-weight: bold; color: palette(highlight);")

        self.paper = Paper(self.view)
        self.view.setScene(self.paper)
        self.paper.setSize(2000, 2000)
        self.view.centerOn(0, 0)
        
        from ..tools import (TemplateTool, ArrowTool, SelectTool, TextTool,
                             RotateTool, ChargeTool, LonePairTool, RadicalPlusTool,
                             RadicalMinusTool)
        self.tools = {
            "bond": StructureTool(),
            "eraser": EraserTool(),
            "text": TextTool(),
            "select": SelectTool(),
            "rotate": RotateTool(),
            "charge_plus": ChargeTool(1),
            "charge_minus": ChargeTool(-1),
            "lone_pair": LonePairTool(),
            "radical_plus": RadicalPlusTool(),
            "radical_minus": RadicalMinusTool(),
            "arrow": ArrowTool("reaction", curvature=0.0),
            "equilibrium": ArrowTool("equilibrium", curvature=0.0),
            "reversible": ArrowTool("reversible", curvature=0.0),
            "curly": ArrowTool("curly", curvature=1.0),
            "curly_cw": ArrowTool("curly", curvature=-1.0),
            "fish_up": ArrowTool("fish_up", curvature=1.0),
            "fish_down": ArrowTool("fish_down", curvature=-1.0)
        }
        for name in ring_templates:
            self.tools[name] = TemplateTool(name)
        self.current_tool = self.tools["select"]
        App.tool = self.current_tool

        self.tool_group.triggered.connect(self._on_tool_changed)
        self.paper.text_editing_finished.connect(self._on_text_editing_finished)
        self.paper.context_menu_requested.connect(self._on_context_menu)

    # ------------------------------------------------------------ context menu
    def _on_context_menu(self, obj, x, y):
        menu = self._build_context_menu(obj)
        point = self.view.viewport().mapToGlobal(self.view.mapFromScene(QPointF(x, y)))
        menu.exec(point)

    def _build_context_menu(self, obj):
        from ..atom import Atom
        from ..bond import Bond

        menu = QMenu(self)
        if isinstance(obj, Atom):
            elements = menu.addMenu("Element")
            for symbol in common_elements:
                action = elements.addAction(symbol)
                action.triggered.connect(
                    lambda checked=False, s=symbol, a=obj: self._menu_set_symbol(a, s))
            charge = menu.addMenu("Charge")
            for label, value in (("Increase (+)", 1), ("Decrease (−)", -1), ("Neutral", 0)):
                action = charge.addAction(label)
                action.triggered.connect(
                    lambda checked=False, v=value, a=obj: self._menu_set_charge(a, v))
            action = menu.addAction("Edit label…")
            action.triggered.connect(lambda checked=False, a=obj: self._menu_edit_label(a))
            menu.addSeparator()
            action = menu.addAction("Delete atom")
            action.triggered.connect(lambda checked=False, a=obj: self._menu_delete(a))
        elif isinstance(obj, Bond):
            types = menu.addMenu("Bond type")
            for name in ("single", "double", "triple", "wedge", "hashed_wedge"):
                action = types.addAction(name.replace("_", " "))
                action.triggered.connect(
                    lambda checked=False, t=name, b=obj: self._menu_set_bond_type(b, t))
            action = menu.addAction("Flip wedge / swap ends")
            action.triggered.connect(lambda checked=False, b=obj: self._menu_flip_bond(b))
            menu.addSeparator()
            action = menu.addAction("Delete bond")
            action.triggered.connect(lambda checked=False, b=obj: self._menu_delete(b))
        elif obj is not None:
            action = menu.addAction("Delete")
            action.triggered.connect(lambda checked=False, o=obj: self._menu_delete(o))
        else:
            menu.addAction(self.action_select_all)
            menu.addAction(self.action_paste)
            menu.addAction(self.action_clear)

        if isinstance(obj, (Atom, Bond)):
            menu.addSeparator()
            action = menu.addAction("Clean up structure")
            action.triggered.connect(
                lambda checked=False, m=obj.molecule: self._cleanup_molecule(m))
        return menu

    def _menu_set_symbol(self, atom, symbol):
        atom.set_symbol(symbol)
        atom.draw()
        self.paper.save_state_to_undo_stack("Set element")

    def _menu_set_charge(self, atom, value):
        atom.set_charge(atom.charge + value if value else 0)
        atom.draw()
        self.paper.save_state_to_undo_stack("Set charge")

    def _menu_edit_label(self, atom):
        text, ok = QInputDialog.getText(self, "Atom label",
                                        "Element or group (OH, Ph, COOH, NO2 ...):",
                                        text=atom.symbol)
        if not ok or not text.strip():
            return
        text = text.strip()
        from ..fragments import expand_label
        if text not in periodic_table and expand_label(atom, text):
            self.paper.save_state_to_undo_stack("Expand group")
            return
        atom.set_symbol(text)
        atom.draw()
        self.paper.save_state_to_undo_stack("Set label")

    def _menu_set_bond_type(self, bond, bond_type):
        bond.set_type(bond_type)
        bond.draw()
        for atom in bond.atoms:
            atom.draw()
        self.paper.save_state_to_undo_stack("Set bond type")

    def _menu_flip_bond(self, bond):
        bond.reverse()
        bond.draw()
        self.paper.save_state_to_undo_stack("Flip bond")

    def _menu_delete(self, obj):
        from ..tools import delete_object
        if delete_object(obj):
            self.paper.redraw_dirty_objects()
            self.paper.save_state_to_undo_stack("Delete")

    def _cleanup_molecule(self, mol):
        """ regenerates a tidy 2D layout, ChemDraw's "Clean Up Structure" """
        if mol is None or len(mol.atoms) < 2:
            return
        from src.vendors.oasa.coords_generator import coords_generator
        center = mol.get_center()
        try:
            coords_generator().calculate_coords(mol, bond_length=1, force=1)
        except Exception as e:
            QMessageBox.warning(self, "Clean Up", "Could not lay out the structure: %s" % e)
            return
        scale = Settings.bond_length
        cx = sum(a.x for a in mol.atoms) / len(mol.atoms)
        cy = sum(a.y for a in mol.atoms) / len(mol.atoms)
        for atom in mol.atoms:
            atom.x = center[0] + (atom.x - cx) * scale
            atom.y = center[1] + (atom.y - cy) * scale
            atom.on_bond_count_change()
        mol.scale_val = 1.0
        for bond in mol.bonds:
            if bond.auto_second_line_side:
                bond.second_line_side = None
        mol.draw()
        self.paper.save_state_to_undo_stack("Clean up structure")

    def _on_cleanup(self):
        mol = self._target_molecule()
        if not mol:
            QMessageBox.warning(self, "Clean Up", "There is no structure to lay out.")
            return
        self._cleanup_molecule(mol)

    def _on_text_editing_finished(self):
        from ..tools import TextTool
        from ..text_label import TextLabel
        # a label that was left empty is invisible and can never be clicked
        # again, so it is dropped instead of littering the page
        for obj in list(self.paper.objects):
            if isinstance(obj, TextLabel) and not (obj.text or "").strip():
                if self.paper.focused_obj is obj:
                    self.paper.focused_obj = None
                if self.paper.locked_focus_obj is obj:
                    self.paper.locked_focus_obj = None
                obj.delete_from_paper()
        # ChemDraw style: after text is finished, switch back to select tool
        if isinstance(self.current_tool, TextTool):
            self._activate_select_tool()

    def wheelEvent(self, event):
        angle = event.angleDelta().y()
        factor = 1.1 if angle > 0 else 0.9
        self._zoom(factor, at_cursor=True)
        super().wheelEvent(event)

    def keyPressEvent(self, event):
        if App.tool and hasattr(App.tool, 'on_key_press'):
            if App.tool.on_key_press(event.key(), event.text()):
                return

        if event.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            from ..tools import SelectTool
            from ..atom import Atom
            from ..bond import Bond
            from ..molecule import Molecule
            if isinstance(self.current_tool, SelectTool) and self.current_tool.objs:
                for obj in self.current_tool.objs[:]:
                    if isinstance(obj, Molecule):
                        # Delete entire molecule
                        obj.delete_from_paper()
                    elif isinstance(obj, Atom):
                        # Delete single atom: disconnect bonds, remove from molecule
                        mol = obj.molecule
                        if mol:
                            for bond in list(obj.neighbor_edges):
                                bond.clear_drawings()
                                bond.disconnect_atoms()
                                mol.remove_bond(bond)
                            mol.remove_atom(obj)
                            obj.clear_drawings()
                            obj.paper = None
                            # Split disconnected fragments
                            mol.split_fragments()
                            # Remove molecule if empty
                            if len(mol.atoms) == 0:
                                mol.delete_from_paper()
                            else:
                                mol.draw()
                    elif isinstance(obj, Bond):
                        # Delete single bond: disconnect, split fragments
                        mol = obj.molecule
                        if mol:
                            obj.clear_drawings()
                            obj.disconnect_atoms()
                            mol.remove_bond(obj)
                            mol.split_fragments()
                            mol.draw()
                    else:
                        # TextLabel, Arrow, or other drawable object
                        obj.delete_from_paper()
                
                self.current_tool.objs = []
                self.current_tool._move_targets = []
                for o in self.paper.objects:
                    if hasattr(o, 'set_selected'):
                        o.set_selected(False)
                self.paper.changeFocusTo(None)
                self.paper.locked_focus_obj = None
                # the state is saved after the deletion, otherwise the first
                # redo would bring the deleted objects back
                self.paper.save_state_to_undo_stack("Delete")
        super().keyPressEvent(event)

    def _zoom(self, factor, at_cursor=False):
        anchor = (QGraphicsView.ViewportAnchor.AnchorUnderMouse if at_cursor
                  else QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self.view.setTransformationAnchor(anchor)
        self.view.scale(factor, factor)
        self.view.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)

    def _zoom_reset(self):
        self.view.resetTransform()

    def showEvent(self, event):
        App.paper = self.paper
        App.tool = self.current_tool
        super().showEvent(event)

    def _activate_select_tool(self):
        self.action_select.setChecked(True)
        self._on_tool_changed(self.action_select)

    def _on_tool_changed(self, action):
        # let the tool that is going away drop its previews and half done work
        previous = getattr(self, 'current_tool', None)
        if previous and hasattr(previous, 'clear'):
            previous.clear()

        if action == self.action_bond:
            self.current_tool = self.tools["bond"]
        elif action == self.action_eraser:
            self.current_tool = self.tools["eraser"]
        elif action == self.action_text:
            self.current_tool = self.tools["text"]
        elif action == self.action_select:
            self.current_tool = self.tools["select"]
        elif action == self.action_rotate:
            self.current_tool = self.tools["rotate"]
        elif action == self.action_charge_plus:
            self.current_tool = self.tools["charge_plus"]
        elif action == self.action_charge_minus:
            self.current_tool = self.tools["charge_minus"]
        elif action == self.action_lone_pair:
            self.current_tool = self.tools["lone_pair"]
        elif action == self.action_radical_plus:
            self.current_tool = self.tools["radical_plus"]
        elif action == self.action_radical_minus:
            self.current_tool = self.tools["radical_minus"]
        elif action == self.action_ring:
            self.current_tool = self.tools[self.ring_combo.currentText()]
        elif action == self.action_arrow:
            self.current_tool = self.tools["arrow"]
        elif action == self.action_equilibrium:
            self.current_tool = self.tools["equilibrium"]
        elif action == self.action_reversible:
            self.current_tool = self.tools["reversible"]
        elif action == self.action_curly:
            self.current_tool = self.tools["curly"]
        elif action == self.action_curly_cw:
            self.current_tool = self.tools["curly_cw"]
        elif action == self.action_fish_up:
            self.current_tool = self.tools["fish_up"]
        elif action == self.action_fish_down:
            self.current_tool = self.tools["fish_down"]
        elif action == self.action_rectangle:
            self.current_tool = ShapeTool(Rectangle)
        elif action == self.action_ellipse:
            self.current_tool = ShapeTool(Ellipse)
        App.tool = self.current_tool

    def _on_bond_type_changed(self, text):
        toolsettings['bond_type'] = text
        # picking a bond type is only useful with the bond tool in hand
        if not self.action_bond.isChecked():
            self.action_bond.setChecked(True)
            self._on_tool_changed(self.action_bond)

    def _on_element_changed(self, text):
        text = (text or "").strip()
        if not text:
            return
        toolsettings['structure'] = text
        if not self.action_bond.isChecked():
            self.action_bond.setChecked(True)
            self._on_tool_changed(self.action_bond)

    def _on_ring_changed(self, text):
        if text not in self.tools:
            return
        self.action_ring.setChecked(True)
        self._on_tool_changed(self.action_ring)

    def _on_show_carbon_changed(self, text):
        self.paper.show_carbon = text
        for obj in self.paper.objects:
            if hasattr(obj, 'atoms'):
                for atom in obj.atoms:
                    atom.visible = None
                obj.draw()

    def _on_clear(self):
        if self.paper.objects:
            answer = QMessageBox.question(
                self, "Clear", "Remove everything from the page?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes:
                return
        self._on_clear_confirmed()

    def _on_clear_confirmed(self):
        # clear_all() also drops the python side references, a plain
        # QGraphicsScene.clear() would leave them pointing at deleted items
        self.paper.clear_all()
        for tool in self.tools.values():
            if hasattr(tool, 'clear'):
                tool.clear()
        if hasattr(self.current_tool, 'objs'):
            self.current_tool.objs = []
            self.current_tool._move_targets = []

    def _on_import(self):
        target_mol = self._target_molecule()
        if not target_mol:
            QMessageBox.warning(self, "Import", "There is no structure to import.")
            return
        try:
            smiles = self._molecule_smiles(target_mol)
            if smiles:
                self.molecule_imported.emit(smiles)
        except Exception as e:
            QMessageBox.critical(
                self, "Import",
                "Could not convert the structure to SMILES: %s" % e)

    def _on_undo(self):
        self.paper.undo()

    def _on_redo(self):
        self.paper.redo()

    # ------------------------------------------------------------ file / SMILES
    def _target_molecule(self):
        """ the molecule the export actions work on: selected one, else the last """
        from ..tools import SelectTool
        if isinstance(self.current_tool, SelectTool) and self.current_tool.objs:
            for obj in self.current_tool.objs:
                if obj.class_name == "Molecule":
                    return obj
                if getattr(obj, 'molecule', None):
                    return obj.molecule
                parent = getattr(obj, 'parent', None)
                if parent is not None and getattr(parent, 'class_name', None) == "Molecule":
                    return parent
        mols = [obj for obj in self.paper.objects if obj.class_name == "Molecule"]
        return mols[-1] if mols else None

    def _molecule_smiles(self, mol):
        from ..fileformat_smiles import Smiles
        groups = [a.symbol for a in mol.atoms if a.symbol not in periodic_table]
        if groups:
            raise ValueError("the structure contains unexpanded labels: %s"
                             % ", ".join(sorted(set(groups))))
        return Smiles().generate(mol)

    def _on_copy_smiles(self):
        mol = self._target_molecule()
        if not mol:
            QMessageBox.warning(self, "SMILES", "There is no structure to convert.")
            return
        try:
            smiles = self._molecule_smiles(mol)
        except Exception as e:
            QMessageBox.critical(self, "SMILES", "Could not generate SMILES: %s" % e)
            return
        QApplication.clipboard().setText(smiles or "")
        QMessageBox.information(self, "SMILES", "Copied to clipboard:\n\n%s" % smiles)

    def _on_save_mol(self):
        mol = self._target_molecule()
        if not mol:
            QMessageBox.warning(self, "Save", "There is no structure to save.")
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Structure", "", "MDL Molfile (*.mol);;SD File (*.sdf)")
        if not file_path:
            return
        try:
            text = self._to_molfile(mol)
            with open(file_path, "w") as f:
                f.write(text)
                if file_path.lower().endswith(".sdf"):
                    f.write("$$$$\n")
        except Exception as e:
            QMessageBox.critical(self, "Save", "Could not save the structure: %s" % e)

    def _to_molfile(self, mol):
        """ MDL molfile of a sketched molecule (y is flipped, the paper grows down) """
        orders = {"single": 1, "double": 2, "triple": 3, "wedge": 1, "hashed_wedge": 1}
        stereo = {"wedge": 1, "hashed_wedge": 6}
        atoms = list(mol.atoms)
        index = {a: i + 1 for i, a in enumerate(atoms)}
        scale = mol.preferred_bond_length() or 1.0
        lines = ["", "  PyChem 2D Sketcher", "",
                 "%3d%3d  0  0  0  0  0  0  0  0999 V2000" % (len(atoms), len(mol.bonds))]
        for a in atoms:
            symbol = a.symbol if a.symbol in periodic_table else "C"
            lines.append("%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0"
                         % (a.x / scale, -a.y / scale, 0.0, symbol))
        for b in mol.bonds:
            if len(b.atoms) != 2:
                continue
            lines.append("%3d%3d%3d%3d  0  0  0" % (index[b.atoms[0]], index[b.atoms[1]],
                                                    orders.get(b.type, 1),
                                                    stereo.get(b.type, 0)))
        charged = [(index[a], a.charge) for a in atoms if a.charge]
        for i in range(0, len(charged), 8):
            chunk = charged[i:i + 8]
            lines.append("M  CHG%3d%s" % (len(chunk),
                         "".join("%4d%4d" % c for c in chunk)))
        lines.append("M  END")
        return "\n".join(lines) + "\n"

    def _on_open_mol(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Structure", "", "Molfiles (*.mol *.sdf *.mdl);;All files (*)")
        if not file_path:
            return
        try:
            with open(file_path) as f:
                text = f.read()
            self._from_molfile(text)
        except Exception as e:
            import traceback
            traceback.print_exc()
            QMessageBox.critical(self, "Open", "Could not read the structure: %s" % e)

    def _from_molfile(self, text):
        from ..molecule import Molecule
        lines = text.splitlines()
        if len(lines) < 4:
            raise ValueError("not a molfile")
        counts = lines[3]
        n_atoms, n_bonds = int(counts[0:3]), int(counts[3:6])
        view_rect = self.view.mapToScene(self.view.viewport().rect()).boundingRect()
        cx, cy = view_rect.center().x(), view_rect.center().y()
        scale = Settings.bond_length

        mol = Molecule()
        self.paper.addObject(mol)
        atoms = []
        for line in lines[4:4 + n_atoms]:
            x, y = float(line[0:10]), float(line[10:20])
            symbol = line[31:34].strip() or "C"
            atom = mol.new_atom(symbol if symbol in periodic_table else "C")
            atom.set_pos(cx + x * scale, cy - y * scale)
            atoms.append(atom)
        types = {1: "single", 2: "double", 3: "triple", 4: "double"}
        for line in lines[4 + n_atoms:4 + n_atoms + n_bonds]:
            a1, a2, order = int(line[0:3]), int(line[3:6]), int(line[6:9])
            stereo = int(line[9:12]) if len(line) >= 12 and line[9:12].strip() else 0
            bond = mol.new_bond()
            bond.set_type({1: "wedge", 6: "hashed_wedge"}.get(stereo, types.get(order, "single")))
            bond.connect_atoms(atoms[a1 - 1], atoms[a2 - 1])
        for line in lines[4 + n_atoms + n_bonds:]:
            if line.startswith("M  CHG"):
                values = line[6:].split()
                for i in range(1, len(values) - 1, 2):
                    atoms[int(values[i]) - 1].set_charge(int(values[i + 1]))
            elif line.startswith("M  END") or line.startswith("$$$$"):
                break
        mol.draw()
        self.paper.save_state_to_undo_stack("Open molfile")

    def _on_flip_h(self):
        from ..tools import SelectTool
        if isinstance(self.current_tool, SelectTool) and self.current_tool._move_targets:
            self.paper.save_state_to_undo_stack("Flip Horizontal")
            # Flip selected objects horizontally around their combined center
            bboxes = []
            for obj in self.current_tool._move_targets:
                if hasattr(obj, 'bounding_box'):
                    bb = obj.bounding_box()
                    if bb: bboxes.append(bb)
            
            if not bboxes: return
            from ..common import bbox_of_bboxes
            whole_bbox = bbox_of_bboxes(bboxes)
            center_x = (whole_bbox[0] + whole_bbox[2]) / 2
            
            for obj in self.current_tool._move_targets:
                if hasattr(obj, 'flip_horizontal'):
                    obj.flip_horizontal(center_x)
                # Ensure redraw
                if hasattr(obj, 'draw'):
                    obj.draw()
                # If it's an atom, redraw connected bonds
                if hasattr(obj, 'neighbor_edges'):
                    for b in obj.neighbor_edges:
                        b.draw()
            self.paper.update()

    def _on_flip_v(self):
        from ..tools import SelectTool
        if isinstance(self.current_tool, SelectTool) and self.current_tool._move_targets:
            self.paper.save_state_to_undo_stack("Flip Vertical")
            # Flip selected objects vertically around their combined center
            bboxes = []
            for obj in self.current_tool._move_targets:
                if hasattr(obj, 'bounding_box'):
                    bb = obj.bounding_box()
                    if bb: bboxes.append(bb)
            
            if not bboxes: return
            from ..common import bbox_of_bboxes
            whole_bbox = bbox_of_bboxes(bboxes)
            center_y = (whole_bbox[1] + whole_bbox[3]) / 2
            
            for obj in self.current_tool._move_targets:
                if hasattr(obj, 'flip_vertical'):
                    obj.flip_vertical(center_y)
                # Ensure redraw
                if hasattr(obj, 'draw'):
                    obj.draw()
                # If it's an atom, redraw connected bonds
                if hasattr(obj, 'neighbor_edges'):
                    for b in obj.neighbor_edges:
                        b.draw()
            self.paper.update()

    def _on_select_all(self):
        from ..tools import SelectTool
        if not isinstance(self.current_tool, SelectTool):
            self.action_select.setChecked(True)
            self._on_tool_changed(self.action_select)
        
        # Select everything top-level
        self.current_tool.objs = list(self.paper.objects)
        self.current_tool._move_targets = list(self.paper.objects)
        for o in self.paper.objects:
            if hasattr(o, 'set_selected'):
                o.set_selected(True)
        self.paper.update()

    def _on_page_setup(self):
        from src.shared.qt_compat import QDialog, QFormLayout, QSpinBox, QDialogButtonBox, QVBoxLayout, QComboBox, QRadioButton, QHBoxLayout, QLabel
        dialog = QDialog(self)
        dialog.setWindowTitle("Page Setup")
        layout = QVBoxLayout(dialog)
        
        form = QFormLayout()
        
        # Presets
        preset_combo = QComboBox()
        preset_combo.addItems(["Custom", "A4", "Letter", "Legal"])
        form.addRow("Paper Size:", preset_combo)
        
        # Units
        unit_combo = QComboBox()
        unit_combo.addItems(["Pixels", "Inches"])
        form.addRow("Units:", unit_combo)
        
        # Orientation
        orient_layout = QHBoxLayout()
        portrait_radio = QRadioButton("Portrait")
        landscape_radio = QRadioButton("Landscape")
        portrait_radio.setChecked(True)
        orient_layout.addWidget(portrait_radio)
        orient_layout.addWidget(landscape_radio)
        form.addRow("Orientation:", orient_layout)
        
        w_spin = QSpinBox()
        w_spin.setRange(100, 10000)
        w_spin.setValue(int(self.paper.width()))
        
        h_spin = QSpinBox()
        h_spin.setRange(100, 10000)
        h_spin.setValue(int(self.paper.height()))
        
        form.addRow("Width:", w_spin)
        form.addRow("Height:", h_spin)
        layout.addLayout(form)
        
        # Logic for presets and units
        def update_values():
            unit = unit_combo.currentText()
            preset = preset_combo.currentText()
            is_portrait = portrait_radio.isChecked()
            
            # DPI = 100 for conversion
            dpi = 100.0
            
            if preset != "Custom":
                if preset == "A4":
                    w, h = 8.27, 11.69
                elif preset == "Letter":
                    w, h = 8.5, 11.0
                elif preset == "Legal":
                    w, h = 8.5, 14.0
                
                if not is_portrait:
                    w, h = h, w
                
                if unit == "Pixels":
                    w_spin.setValue(int(w * dpi))
                    h_spin.setValue(int(h * dpi))
                else:
                    w_spin.setValue(int(w))
                    h_spin.setValue(int(h))
            
            if unit == "Inches":
                w_spin.setRange(1, 100)
                h_spin.setRange(1, 100)
            else:
                w_spin.setRange(100, 10000)
                h_spin.setRange(100, 10000)

        preset_combo.currentIndexChanged.connect(update_values)
        unit_combo.currentIndexChanged.connect(update_values)
        portrait_radio.toggled.connect(update_values)
        
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)
        
        if dialog.exec():
            w = w_spin.value()
            h = h_spin.value()
            if unit_combo.currentText() == "Inches":
                w *= 100
                h *= 100
            self.paper.setSize(w, h)
            self.paper.save_state_to_undo_stack("Page Setup")

    def _on_color(self):
        from src.shared.qt_compat import QColorDialog, QColor
        
        # Get current color from first selected object if possible
        initial_color = QColor(0, 0, 0)
        from ..tools import SelectTool
        if isinstance(self.current_tool, SelectTool) and self.current_tool.objs:
            obj = self.current_tool.objs[0]
            if hasattr(obj, 'color'):
                c = obj.color
                initial_color = QColor(c[0], c[1], c[2])
        
        color = QColorDialog.getColor(initial_color, self, "Select Color")
        if color.isValid():
            new_color = (color.red(), color.green(), color.blue())
            
            # Apply to selected
            if isinstance(self.current_tool, SelectTool):
                targets = self.current_tool.objs
                if not targets: return
                
                for obj in targets:
                    self._apply_color_recursive(obj, new_color)
                
                self.paper.save_state_to_undo_stack("Change Color")
                self.paper.update()

    def _apply_color_recursive(self, obj, color):
        # Apply to the object itself
        if hasattr(obj, 'color'):
            obj.color = color
        
        # If it's a molecule, apply to atoms and bonds
        # Only do this if the object is explicitly a Molecule
        from ..molecule import Molecule
        if isinstance(obj, Molecule):
            if hasattr(obj, 'atoms'):
                for a in obj.atoms:
                    a.color = color
                    a.draw()
            if hasattr(obj, 'bonds'):
                for b in obj.bonds:
                    b.color = color
                    b.draw()
        
        if hasattr(obj, 'draw'):
            obj.draw()

    def _on_copy(self):
        from ..tools import SelectTool
        if isinstance(self.current_tool, SelectTool):
            App.clipboard = []
            seen_mols = set()
            # Use both objs and _move_targets to ensure we capture what the user intended
            targets = self.current_tool.objs
            if not targets and hasattr(self.current_tool, '_move_targets'):
                targets = self.current_tool._move_targets
                
            for obj in targets:
                if hasattr(obj, 'clone'):
                    App.clipboard.append(obj.clone())
                elif hasattr(obj, 'molecule') and obj.molecule:
                    mol = obj.molecule
                    if id(mol) not in seen_mols:
                        App.clipboard.append(mol.clone())
                        seen_mols.add(id(mol))
                elif hasattr(obj, 'parent') and obj.parent and hasattr(obj.parent, 'clone'):
                    if id(obj.parent) not in seen_mols:
                        App.clipboard.append(obj.parent.clone())
                        seen_mols.add(id(obj.parent))
            print(f"Copied {len(App.clipboard)} objects to clipboard")

    def _on_paste(self):
        if not App.clipboard:
            return
        
        new_objs = []
        offset = 20
        
        # Deselect current selection before pasting
        for o in self.paper.objects:
            if hasattr(o, 'set_selected'):
                o.set_selected(False)
        
        for obj in App.clipboard:
            if not hasattr(obj, 'clone'): continue
            new_obj = obj.clone()
            # Add to paper BEFORE moving or drawing so that paper property is available
            self.paper.addObject(new_obj)
            
            if hasattr(new_obj, 'move_by'):
                new_obj.move_by(offset, offset)
            elif hasattr(new_obj, 'atoms'):
                for a in new_obj.atoms:
                    a.move_by(offset, offset)
            
            new_obj.draw()
            new_objs.append(new_obj)
        
        # Select the newly pasted objects
        from ..tools import SelectTool
        if isinstance(self.current_tool, SelectTool):
            self.current_tool.objs = new_objs
            self.current_tool._move_targets = list(new_objs)
            for o in self.paper.objects:
                if hasattr(o, 'set_selected'):
                    o.set_selected(o in new_objs)
        
        self.paper.save_state_to_undo_stack("Paste")
        
        # Shift clipboard for consecutive pastes
        for obj in App.clipboard:
            if hasattr(obj, 'move_by'):
                obj.move_by(offset, offset)
            elif hasattr(obj, 'atoms'):
                # Handle molecules that don't have move_by (though they should)
                for a in obj.atoms:
                    a.move_by(offset, offset)

    def _on_export(self):
        from ..tools import SelectTool
        from ..common import bbox_of_bboxes
        from ..paper import Paper
        
        # 1. Identify selected objects (top-level)
        selected_objs = []
        if isinstance(self.current_tool, SelectTool):
            targets = []
            if hasattr(self.current_tool, '_move_targets') and self.current_tool._move_targets:
                targets = self.current_tool._move_targets
            else:
                targets = self.current_tool.objs
                
            for obj in targets:
                top = obj
                if hasattr(obj, 'molecule') and obj.molecule:
                    top = obj.molecule
                elif hasattr(obj, 'parent') and obj.parent:
                    top = obj.parent
                
                if top not in selected_objs:
                    selected_objs.append(top)
        
        if not selected_objs:
            QMessageBox.warning(self, "Export", "Please select the objects you want to export.")
            return

        # 2. Get combined bounding box
        bboxes = []
        for obj in selected_objs:
            if hasattr(obj, 'bounding_box'):
                bb = obj.bounding_box()
                if bb: bboxes.append(bb)
        
        if not bboxes:
            QMessageBox.warning(self, "Export", "Could not determine the bounding box of selected objects.")
            return
            
        whole_bbox = bbox_of_bboxes(bboxes)
        margin = 10
        x1, y1, x2, y2 = whole_bbox
        x1 -= margin; y1 -= margin; x2 += margin; y2 += margin
        width = x2 - x1
        height = y2 - y1

        # 3. Get Export Settings
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export as Image", "", 
            "PNG Image (*.png);;JPEG Image (*.jpg *.jpeg);;BMP Image (*.bmp)"
        )
        if not file_path: return
        
        dpi, ok = QInputDialog.getInt(self, "Export DPI", "Enter DPI:", value=300, minValue=72, maxValue=1200)
        if not ok: return

        # 4. Render to high-DPI image
        scale = dpi / 100.0
        img = QImage(int(width * scale), int(height * scale), QImage.Format.Format_ARGB32)
        img.fill(Qt.GlobalColor.white)
        
        painter = QPainter(img)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        
        # Use a temporary scene to render only selected objects
        temp_paper = Paper()
        for obj in selected_objs:
            if hasattr(obj, 'clone'):
                clone = obj.clone()
                temp_paper.addObject(clone)
                clone.draw()
        
        # Render the specific source rectangle into the image
        temp_paper.render(painter, QRectF(0, 0, width * scale, height * scale), QRectF(x1, y1, width, height))
        painter.end()
        
        if img.save(file_path):
            QMessageBox.information(self, "Export", f"Successfully exported to {file_path}")
        else:
            QMessageBox.critical(self, "Export", "Failed to save image.")

    def _on_smiles_to_2d(self):
        text, ok = QInputDialog.getText(self, "SMILES to 2D", "Enter SMILES string:")
        if not ok or not text.strip():
            return
            
        try:
            from src.features.smiles_parser.services.parser import parse_smiles
            from ..molecule import Molecule
            from ..atom import Atom
            from ..bond import Bond
            from src.core.domain.models.bond import BondType
            
            # parse_smiles uses OASA's text_to_mol(calc_coords=1,
            # localize_aromatic_bonds=True), which already:
            #   1. Generates high-quality 2D coordinates (x2d/y2d)
            #   2. Kekulizes aromatic bonds (alternating single/double)
            # No need to call kekulize or CoordinateGenerator2DSMILES.
            mol_domain = parse_smiles(text.strip())
            
            # Place at center of view
            view_rect = self.view.mapToScene(self.view.viewport().rect()).boundingRect()
            cx, cy = view_rect.center().x(), view_rect.center().y()
            
            new_mol = Molecule()
            self.paper.addObject(new_mol)
            
            atom_map = {}
            scale = 40  # Scale for better visibility
            
            for atom in mol_domain.atoms:
                if atom.symbol == 'H':
                    continue  # Skip explicit H atoms in sketcher
                new_atom = Atom(atom.symbol)
                new_atom.charge = atom.formal_charge
                if hasattr(atom, 'x2d') and atom.x2d is not None:
                    new_atom.x = cx + atom.x2d * scale
                    new_atom.y = cy + atom.y2d * scale
                else:
                    new_atom.x, new_atom.y = cx, cy
                new_mol.add_atom(new_atom)
                atom_map[atom.index] = new_atom
            
            for bond in mol_domain.bonds:
                # Skip bonds involving H atoms
                if bond.begin_atom_idx not in atom_map or bond.end_atom_idx not in atom_map:
                    continue
                new_bond = Bond()
                btype = "single"
                if bond.bond_type == BondType.DOUBLE: btype = "double"
                elif bond.bond_type == BondType.TRIPLE: btype = "triple"
                # OASA already kekulized aromatic bonds to single/double,
                # so BondType.AROMATIC should be rare here. Treat as single.
                
                new_bond.set_type(btype)
                new_bond.connect_atoms(atom_map[bond.begin_atom_idx], atom_map[bond.end_atom_idx])
                new_mol.add_bond(new_bond)
            
            new_mol.draw()
            self.paper.save_state_to_undo_stack("Import SMILES")
            self.view.centerOn(cx, cy)
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            QMessageBox.critical(self, "Error", f"Could not convert SMILES: {str(e)}")
