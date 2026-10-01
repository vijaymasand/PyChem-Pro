"""
Custom Sphere Dialog — Unified UI for adding multi-layered translucent spheres.

Features:
- Position entry with COM defaulting
- Table-based layer management (radius, color, transparency)
- Layer-wise transparency (alpha) support
"""

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
    QPushButton, QTableWidget, QTableWidgetItem,
    QHeaderView, QDoubleSpinBox, QColorDialog,
    QGroupBox, QFormLayout, QFrame
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor

from src.features.visualization_3d.services.dummy_sphere import DummySphereManager


class CustomSphereDialog(QDialog):
    """
    Dialog for defining a multi-layered sphere at specific coordinates.
    """
    
    def __init__(self, molecule, parent=None):
        super().__init__(parent)
        self.molecule = molecule
        self.manager = DummySphereManager(molecule)
        self.com = self.manager.calculate_center_of_mass()
        
        self.setWindowTitle("Add Custom Sphere Shells")
        self.setMinimumWidth(450)
        self.setModal(True)
        
        self._init_ui()
        self._reset_to_com()
        
    def _init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # --- Position Section ---
        pos_group = QGroupBox("Sphere Center Position (Å)")
        pos_layout = QFormLayout(pos_group)
        
        self.spin_x = QDoubleSpinBox()
        self.spin_y = QDoubleSpinBox()
        self.spin_z = QDoubleSpinBox()
        
        for s in [self.spin_x, self.spin_y, self.spin_z]:
            s.setRange(-999.9, 999.9)
            s.setDecimals(3)
            s.setSingleStep(0.1)
            
        pos_layout.addRow("X Coordinate:", self.spin_x)
        pos_layout.addRow("Y Coordinate:", self.spin_y)
        pos_layout.addRow("Z Coordinate:", self.spin_z)
        
        btn_layout = QHBoxLayout()
        reset_com_btn = QPushButton("Reset to COM")
        reset_com_btn.clicked.connect(self._reset_to_com)
        btn_layout.addStretch()
        btn_layout.addWidget(reset_com_btn)
        pos_layout.addRow(btn_layout)
        
        main_layout.addWidget(pos_group)
        
        # --- Layers Section ---
        layer_group = QGroupBox("Sphere Shells/Layers")
        layer_layout = QVBoxLayout(layer_group)
        
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["Radius (Å)", "Color", "Opacity (%)"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setAlternatingRowColors(True)
        
        layer_layout.addWidget(self.table)
        
        btn_row = QHBoxLayout()
        add_btn = QPushButton("Add Layer")
        add_btn.clicked.connect(self._add_layer)
        remove_btn = QPushButton("Remove Selected")
        remove_btn.clicked.connect(self._remove_layer)
        
        btn_row.addWidget(add_btn)
        btn_row.addWidget(remove_btn)
        btn_row.addStretch()
        
        layer_layout.addLayout(btn_row)
        main_layout.addWidget(layer_group)
        
        # Help text
        hint = QLabel("<i>Tip: Outer layers with lower opacity create a shell effect.</i>")
        hint.setStyleSheet("color: #888; font-size: 11px;")
        main_layout.addWidget(hint)
        
        # --- Dialog Buttons ---
        diag_btns = QHBoxLayout()
        ok_btn = QPushButton("OK")
        ok_btn.setDefault(True)
        ok_btn.clicked.connect(self.accept)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        
        diag_btns.addStretch()
        diag_btns.addWidget(ok_btn)
        diag_btns.addWidget(cancel_btn)
        main_layout.addLayout(diag_btns)
        
        # Add default layer
        self._add_layer(radius=1.0, color="#ffff00", alpha=0.5)

    def _reset_to_com(self):
        """Set coordinates to Center of Mass."""
        self.spin_x.setValue(self.com[0])
        self.spin_y.setValue(self.com[1])
        self.spin_z.setValue(self.com[2])
        
    def _add_layer(self, radius=None, color=None, alpha=None):
        """Add a new sphere layer to the table."""
        row = self.table.rowCount()
        self.table.insertRow(row)
        
        # Radius
        radius_spin = QDoubleSpinBox()
        radius_spin.setRange(0.1, 50.0)
        radius_spin.setValue(radius if radius is not None else (row + 1) * 0.5)
        radius_spin.setDecimals(2)
        radius_spin.setSuffix(" Å")
        self.table.setCellWidget(row, 0, radius_spin)
        
        # Color button
        color_val = color if color else "#ffff00"
        color_btn = QPushButton()
        color_btn.setFixedSize(60, 24)
        color_btn.setStyleSheet(f"background-color: {color_val}; border: 1px solid #666;")
        color_btn.setProperty("hex_color", color_val)
        color_btn.clicked.connect(lambda: self._choose_color(color_btn))
        self.table.setCellWidget(row, 1, color_btn)
        
        # Opacity (Alpha as %)
        opacity_spin = QDoubleSpinBox()
        opacity_spin.setRange(0.0, 100.0)
        opacity_spin.setValue((alpha * 100) if alpha is not None else 50.0)
        opacity_spin.setSuffix(" %")
        self.table.setCellWidget(row, 2, opacity_spin)
        
    def _remove_layer(self):
        """Remove selected layer."""
        idx = self.table.currentRow()
        if idx >= 0:
            self.table.removeRow(idx)
        elif self.table.rowCount() > 0:
            self.table.removeRow(self.table.rowCount() - 1)
            
    def _choose_color(self, btn):
        """Pick a color for the specific layer."""
        current = QColor(btn.property("hex_color"))
        color = QColorDialog.getColor(current, self, "Pick Shell Color")
        if color.isValid():
            hex_c = color.name().lower()
            btn.setStyleSheet(f"background-color: {hex_c}; border: 1px solid #666;")
            btn.setProperty("hex_color", hex_c)
            
    def get_result(self):
        """
        Get the defined sphere parameters.
        
        Returns:
            Tuple: (x, y, z), List of dicts [(radius, color, alpha), ...]
        """
        pos = (self.spin_x.value(), self.spin_y.value(), self.spin_z.value())
        layers = []
        
        for row in range(self.table.rowCount()):
            radius = self.table.cellWidget(row, 0).value()
            color = self.table.cellWidget(row, 1).property("hex_color")
            alpha = self.table.cellWidget(row, 2).value() / 100.0
            layers.append({
                'radius': radius,
                'color': color,
                'alpha': alpha
            })
            
        return pos, layers
