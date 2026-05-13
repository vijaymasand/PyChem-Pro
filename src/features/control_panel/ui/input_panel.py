"""
Input Panel — SMILES input, conversion controls, and molecule info display.
"""

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QGroupBox, QProgressBar, QTextEdit,
    QComboBox, QCheckBox, QFrame, QSizePolicy, QSpinBox,
    QScrollArea, QSlider, Qt, Signal, QFont, QPalette, QColor
)
from src.shared.ui.theme import COLORS, theme_signals


class InputPanel(QWidget):
    """
    Left panel with SMILES input, controls, and molecule information.
    """

    convert_requested = Signal(str)
    export_sdf_requested = Signal()
    export_mol2_requested = Signal()
    export_image_requested = Signal(int, bool)  # (dpi, white_bg) — white_bg always True now
    sphere_scale_changed = Signal(float)
    stick_scale_changed = Signal(float)
    line_scale_changed = Signal(float)
    show_sidechains_changed = Signal(bool)
    show_sasa_changed = Signal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        # objectName 'leftPanel' is kept for the global stylesheet
        # hairline-border rule. The background fill is painted via
        # an explicit inline stylesheet applied in _apply_theme()
        # because Qt's QScrollArea palette caching (on macOS) refuses
        # to follow a live theme swap through the app-level stylesheet
        # cascade.  Inline stylesheets win and give us reliable swaps.
        self.setObjectName("leftPanel")
        self.setFixedWidth(300)
        self._init_ui()
        self._apply_theme()
        try:
            theme_signals().theme_changed.connect(self._apply_theme)
        except Exception:
            pass

    def _apply_theme(self):
        """
        Paint the panel, container, and scrollbar from the CURRENT
        COLORS values. Two mechanisms are used in parallel:

        * An inline Qt stylesheet — governs the rounded shapes,
          padding, hover/focus states.
        * A QPalette assignment — governs the native macOS QScrollBar
          and QScrollArea paint, which reads from palette roles and
          NOT from the stylesheet cascade.

        Both are needed because neither alone covers every Qt paint
        path on every platform.
        """
        bg = COLORS['bg_secondary']
        border = COLORS['border']
        scroll_bg = COLORS['scrollbar_bg']
        scroll_handle = COLORS['scrollbar_handle']
        accent = COLORS['accent']
        text = COLORS['text_primary']

        bg_color = QColor(bg)
        text_color = QColor(text)
        handle_color = QColor(scroll_handle)

        # ── The panel itself (stylesheet) ────────────────────────
        self.setStyleSheet(
            "QWidget#leftPanel {"
            f"  background-color: {bg};"
            f"  border-right: 1px solid {border};"
            "}"
        )
        # Panel palette (for native child widgets that read it).
        pal = self.palette()
        pal.setColor(QPalette.ColorRole.Window, bg_color)
        pal.setColor(QPalette.ColorRole.Base, bg_color)
        pal.setColor(QPalette.ColorRole.WindowText, text_color)
        self.setPalette(pal)

        # ── Inner scroll-area container ─────────────────────────
        if hasattr(self, '_scroll_container') and self._scroll_container is not None:
            self._scroll_container.setStyleSheet(
                f"background-color: {bg};"
            )
            cpal = self._scroll_container.palette()
            cpal.setColor(QPalette.ColorRole.Window, bg_color)
            cpal.setColor(QPalette.ColorRole.Base, bg_color)
            self._scroll_container.setPalette(cpal)

        # ── Scroll area + its QScrollBar ────────────────────────
        if hasattr(self, '_scroll') and self._scroll is not None:
            self._scroll.setStyleSheet(
                f"QScrollArea {{ background-color: {bg}; border: none; }}"
                f"QScrollBar:vertical {{"
                f"  background: {scroll_bg};"
                f"  width: 10px;"
                f"  margin: 0;"
                f"  border: none;"
                f"}}"
                f"QScrollBar::handle:vertical {{"
                f"  background: {scroll_handle};"
                f"  border-radius: 3px;"
                f"  min-height: 24px;"
                f"  margin: 2px 3px;"
                f"}}"
                f"QScrollBar::handle:vertical:hover {{"
                f"  background: {accent};"
                f"}}"
                f"QScrollBar::add-line:vertical,"
                f"QScrollBar::sub-line:vertical {{"
                f"  height: 0; background: none;"
                f"}}"
                f"QScrollBar::add-page:vertical,"
                f"QScrollBar::sub-page:vertical {{"
                f"  background: {scroll_bg};"
                f"}}"
            )
            # Palette on the scroll area AND its viewport — this is
            # what macOS native QScrollBar paint actually reads.
            spal = self._scroll.palette()
            spal.setColor(QPalette.ColorRole.Window, bg_color)
            spal.setColor(QPalette.ColorRole.Base, bg_color)
            spal.setColor(QPalette.ColorRole.Button, handle_color)
            self._scroll.setPalette(spal)
            vp = self._scroll.viewport()
            if vp is not None:
                vpal = vp.palette()
                vpal.setColor(QPalette.ColorRole.Window, bg_color)
                vpal.setColor(QPalette.ColorRole.Base, bg_color)
                vp.setPalette(vpal)
                vp.setAutoFillBackground(True)
            # The scrollbar itself.
            vbar = self._scroll.verticalScrollBar()
            if vbar is not None:
                bpal = vbar.palette()
                bpal.setColor(QPalette.ColorRole.Window, QColor(scroll_bg))
                bpal.setColor(QPalette.ColorRole.Base, QColor(scroll_bg))
                bpal.setColor(QPalette.ColorRole.Button, handle_color)
                bpal.setColor(QPalette.ColorRole.ButtonText, handle_color)
                vbar.setPalette(bpal)
                vbar.setAutoFillBackground(True)

        # Belt and braces: unpolish + polish every descendant so any
        # cached palettes are discarded.
        for child in self.findChildren(QWidget):
            try:
                st = child.style()
                st.unpolish(child)
                st.polish(child)
                child.update()
            except Exception:
                pass
        self.update()

    # ─── Spacing scale ──────────────────────────────────────────
    SPACE_XS = 4
    SPACE_SM = 6
    SPACE_MD = 10
    SPACE_LG = 14
    SPACE_XL = 20

    def _init_ui(self):
        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.setSpacing(0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        # The scroll area's viewport inherits the palette from its
        # parent but we explicitly fill the inner container in
        # _apply_theme() so nothing shows through with a stale colour.
        scroll.setAutoFillBackground(True)
        self._scroll = scroll

        container = QWidget()
        container.setAutoFillBackground(True)
        self._scroll_container = container
        layout = QVBoxLayout(container)
        # Tight outer padding. The old 16 px all-around felt airy.
        layout.setSpacing(self.SPACE_MD)
        layout.setContentsMargins(self.SPACE_MD, self.SPACE_MD,
                                  self.SPACE_MD, self.SPACE_MD)

        # ── Brand header ────────────────────────────────────────
        title = QLabel("PyChem")
        title.setObjectName("labelTitle")
        title.setAlignment(Qt.AlignmentFlag.AlignLeft)
        layout.addWidget(title)

        subtitle = QLabel("Molecular workbench")
        subtitle.setObjectName("labelSubtitle")
        subtitle.setAlignment(Qt.AlignmentFlag.AlignLeft)
        layout.addWidget(subtitle)

        # ── SMILES Input ────────────────────────────────────────
        # IMPORTANT: inner layout uses ZERO content margins so the
        # QGroupBox's own padding is the ONLY inset. Otherwise the
        # two paddings compound and the input field floats in a sea
        # of dead space.
        input_group = QGroupBox("Input")
        input_layout = QVBoxLayout(input_group)
        input_layout.setSpacing(self.SPACE_XS)
        input_layout.setContentsMargins(0, 0, 0, 0)

        self.smiles_input = QLineEdit()
        self.smiles_input.setPlaceholderText("Enter SMILES...")
        self.smiles_input.setClearButtonEnabled(True)
        self.smiles_input.returnPressed.connect(self._on_convert)
        input_layout.addWidget(self.smiles_input)

        hint = QLabel("e.g.  c1ccccc1   ·   CCO   ·   CC(=O)O")
        hint.setObjectName("labelHint")
        input_layout.addWidget(hint)

        layout.addWidget(input_group)

        # ── Primary action ──────────────────────────────────────
        self.convert_btn = QPushButton("Convert to 3D")
        self.convert_btn.setFixedHeight(32)
        self.convert_btn.setObjectName("btnSuccess")
        self.convert_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.convert_btn.clicked.connect(self._on_convert)
        layout.addWidget(self.convert_btn)

        # ── Progress ────────────────────────────────────────────
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(False)
        self.progress_bar.setTextVisible(False)
        layout.addWidget(self.progress_bar)

        # ── View Options ────────────────────────────────────────
        view_group = QGroupBox("View")
        view_layout = QVBoxLayout(view_group)
        view_layout.setSpacing(self.SPACE_XS)
        view_layout.setContentsMargins(0, 0, 0, 0)

        self.show_h_check = QCheckBox("Show hydrogens")
        self.show_h_check.setChecked(False)
        view_layout.addWidget(self.show_h_check)

        self.show_labels_check = QCheckBox("Show atom labels")
        self.show_labels_check.setChecked(False)
        view_layout.addWidget(self.show_labels_check)

        label_btn_row = QHBoxLayout()
        label_btn_row.setSpacing(4)
        label_btn_row.setContentsMargins(20, 0, 0, 0) # Indent under checkbox
        self.label_color_btn = QPushButton("Color")
        self.label_color_btn.setObjectName("btnTertiary")
        self.label_color_btn.setFixedWidth(50)
        label_btn_row.addWidget(self.label_color_btn)
        
        self.clear_labels_btn = QPushButton("Clear")
        self.clear_labels_btn.setObjectName("btnTertiary")
        self.clear_labels_btn.setFixedWidth(50)
        label_btn_row.addWidget(self.clear_labels_btn)
        label_btn_row.addStretch()
        view_layout.addLayout(label_btn_row)

        self.show_sidechains_check = QCheckBox("Show side chains")
        self.show_sidechains_check.setChecked(False)
        view_layout.addWidget(self.show_sidechains_check)

        self.show_sasa_check = QCheckBox("SASA surface")
        self.show_sasa_check.setChecked(False)
        view_layout.addWidget(self.show_sasa_check)

        self.sasa_selected_only_check = QCheckBox("SASA — selection only")
        self.sasa_selected_only_check.setChecked(False)
        view_layout.addWidget(self.sasa_selected_only_check)

        view_layout.addSpacing(self.SPACE_XS)

        render_row = QHBoxLayout()
        render_row.setSpacing(self.SPACE_SM)
        render_label = QLabel("Style")
        render_label.setObjectName("labelMuted")
        render_label.setFixedWidth(40)
        render_row.addWidget(render_label)
        self.render_combo = QComboBox()
        self.render_combo.addItems([
            "Ball & Stick", "Space Fill", "Wireframe",
            "Cartoon", "Ribbon", "Backbone",
        ])
        render_row.addWidget(self.render_combo, 1)
        view_layout.addLayout(render_row)

        view_layout.addSpacing(self.SPACE_XS)

        self.sphere_slider, self.sphere_label = self._make_slider(
            view_layout, "Sphere", 10, 300, 100, None,
            lambda v: self._on_scale_changed('sphere', v)
        )
        self.stick_slider, self.stick_label = self._make_slider(
            view_layout, "Stick", 10, 300, 100, None,
            lambda v: self._on_scale_changed('stick', v)
        )
        self.line_slider, self.line_label = self._make_slider(
            view_layout, "Line", 10, 300, 100, None,
            lambda v: self._on_scale_changed('line', v)
        )

        layout.addWidget(view_group)

        # ── Export ──────────────────────────────────────────────
        export_group = QGroupBox("Export")
        export_layout = QVBoxLayout(export_group)
        export_layout.setSpacing(self.SPACE_XS)
        export_layout.setContentsMargins(0, 0, 0, 0)

        struct_row = QHBoxLayout()
        struct_row.setSpacing(self.SPACE_SM)
        self.sdf_btn = QPushButton("Save SDF")
        self.sdf_btn.setObjectName("btnSecondary")
        self.sdf_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.sdf_btn.clicked.connect(self.export_sdf_requested.emit)
        self.sdf_btn.setEnabled(False)
        struct_row.addWidget(self.sdf_btn)

        self.mol2_btn = QPushButton("Save MOL2")
        self.mol2_btn.setObjectName("btnSecondary")
        self.mol2_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.mol2_btn.clicked.connect(self.export_mol2_requested.emit)
        self.mol2_btn.setEnabled(False)
        struct_row.addWidget(self.mol2_btn)
        export_layout.addLayout(struct_row)

        dpi_row = QHBoxLayout()
        dpi_row.setSpacing(self.SPACE_SM)
        dpi_label = QLabel("DPI")
        dpi_label.setObjectName("labelMuted")
        dpi_label.setFixedWidth(40)
        dpi_row.addWidget(dpi_label)

        self.dpi_spin = QSpinBox()
        self.dpi_spin.setRange(72, 1200)
        self.dpi_spin.setValue(300)
        self.dpi_spin.setSingleStep(50)
        self.dpi_spin.setFixedWidth(96)
        dpi_row.addWidget(self.dpi_spin)
        dpi_row.addStretch()
        export_layout.addLayout(dpi_row)

        self.export_img_btn = QPushButton("Export image")
        self.export_img_btn.setObjectName("btnSecondary")
        self.export_img_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.export_img_btn.clicked.connect(self._on_export_image)
        self.export_img_btn.setEnabled(False)
        export_layout.addWidget(self.export_img_btn)

        layout.addWidget(export_group)

        # ── Molecule Info ───────────────────────────────────────
        info_group = QGroupBox("Molecule")
        info_layout = QVBoxLayout(info_group)
        info_layout.setContentsMargins(0, 0, 0, 0)
        info_layout.setSpacing(0)

        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setMinimumHeight(96)
        self.info_text.setMaximumHeight(128)
        self.info_text.setFrameShape(QFrame.Shape.NoFrame)
        self.info_text.setStyleSheet(
            # Inside a QGroupBox already, so drop the border and
            # blend with the card fill to remove visual clutter.
            "QTextEdit { border: none; background: transparent; padding: 0; }"
        )
        self.info_text.setPlaceholderText(
            "No molecule loaded."
        )
        info_layout.addWidget(self.info_text)

        layout.addWidget(info_group)

        layout.addStretch()

        # ── Footer ──────────────────────────────────────────────
        footer = QLabel("PyChem  ·  v1.0")
        footer.setObjectName("labelMuted")
        footer.setAlignment(Qt.AlignmentFlag.AlignLeft)
        layout.addWidget(footer)

        scroll.setWidget(container)
        outer_layout.addWidget(scroll)

    def _on_convert(self):
        # QLineEdit exposes .text(); fall back to .toPlainText() if
        # the widget was ever swapped back to a multi-line editor.
        if hasattr(self.smiles_input, 'text'):
            smiles = self.smiles_input.text().strip()
        else:
            smiles = self.smiles_input.toPlainText().strip()
        if smiles:
            self.convert_requested.emit(smiles)

    def set_progress(self, value):
        if value <= 0:
            self.progress_bar.setVisible(False)
        else:
            self.progress_bar.setVisible(True)
            self.progress_bar.setValue(value)

    def update_molecule_info(self, molecule):
        if not molecule:
            self.info_text.clear()
            return

        def _row(key, value):
            return f"{key:<9} {value}"

        info_lines = [
            _row("Formula", molecule.molecular_formula()),
            _row("Weight",  f"{molecule.molecular_weight():.2f} Da"),
            _row("Atoms",   str(len(molecule.atoms))),
            _row("Bonds",   str(len(molecule.bonds))),
            _row("Charge",  str(molecule.total_charge())),
            _row("Rings",   str(len(molecule.find_rings()))),
        ]

        charges = [a.partial_charge for a in molecule.atoms
                   if a.partial_charge != 0]
        if charges:
            info_lines.append(
                _row("Q range",
                     f"{min(charges):+.3f}  to  {max(charges):+.3f}")
            )

        self.info_text.setPlainText("\n".join(info_lines))

    def enable_tools(self, enabled=True):
        self.sdf_btn.setEnabled(enabled)
        self.mol2_btn.setEnabled(enabled)
        self.export_img_btn.setEnabled(enabled)

    def _on_export_image(self):
        dpi = self.dpi_spin.value()
        self.export_image_requested.emit(dpi, True)  # Always white background

    def _make_slider(self, parent_layout, label_text, min_val, max_val,
                     default, style, callback):
        """
        Create a labelled slider row:   LABEL  [slider]  VALUE

        The slider and its label / value rely on the global theme —
        no inline stylesheets are applied so live theme switching
        picks up colour changes automatically.
        """
        row = QHBoxLayout()
        row.setSpacing(10)
        row.setContentsMargins(0, 0, 0, 0)

        label = QLabel(label_text)
        label.setObjectName("labelMuted")
        label.setFixedWidth(44)
        row.addWidget(label)

        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(min_val, max_val)
        slider.setValue(default)
        row.addWidget(slider, 1)

        val_label = QLabel(f"{default}%")
        val_label.setObjectName("labelData")
        val_label.setFixedWidth(38)
        val_label.setAlignment(Qt.AlignmentFlag.AlignRight |
                               Qt.AlignmentFlag.AlignVCenter)
        row.addWidget(val_label)

        slider.valueChanged.connect(lambda v: val_label.setText(f"{v}%"))
        slider.valueChanged.connect(callback)

        parent_layout.addLayout(row)
        return slider, val_label

    def _on_scale_changed(self, which, value):
        """Emit scale changed signal."""
        scale = value / 100.0
        if which == 'sphere':
            self.sphere_scale_changed.emit(scale)
        elif which == 'stick':
            self.stick_scale_changed.emit(scale)
        elif which == 'line':
            self.line_scale_changed.emit(scale)
