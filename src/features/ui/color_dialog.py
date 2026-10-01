"""
Color Customization Dialog for molecular visualization.

Allows users to customize atom colors and selected atom colors.
"""

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
    QPushButton, QColorDialog, QGroupBox, QGridLayout,
    QComboBox, QSpinBox, QCheckBox, QTabWidget
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor

from src.shared.ui.theme import COLORS


class AtomColorDialog(QDialog):
    """
    Dialog for customizing atom colors in molecular visualization.
    """
    
    # Signal emitted when colors are changed
    colors_changed = Signal(dict)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Atom Color Customization")
        self.setModal(True)
        self.resize(500, 600)
        
        # Current color settings
        self.atom_colors = COLORS.copy()
        
        self._init_ui()
        
    def _init_ui(self):
        """Initialize the user interface."""
        layout = QVBoxLayout(self)
        
        # Create tab widget for different color categories
        tab_widget = QTabWidget()
        
        # Atom colors tab
        atom_tab = self._create_atom_colors_tab()
        tab_widget.addTab(atom_tab, "Atom Colors")
        
        # Selection colors tab
        selection_tab = self._create_selection_colors_tab()
        tab_widget.addTab(selection_tab, "Selection Colors")
        
        # Charge colors tab
        charge_tab = self._create_charge_colors_tab()
        tab_widget.addTab(charge_tab, "Charge Colors")
        
        layout.addWidget(tab_widget)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        reset_btn = QPushButton("Reset to Default")
        reset_btn.clicked.connect(self._reset_to_default)
        
        apply_btn = QPushButton("Apply")
        apply_btn.clicked.connect(self._apply_changes)
        
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        
        button_layout.addWidget(reset_btn)
        button_layout.addStretch()
        button_layout.addWidget(apply_btn)
        button_layout.addWidget(cancel_btn)
        
        layout.addLayout(button_layout)
        
    def _create_atom_colors_tab(self):
        """Create tab for element-specific atom colors."""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        # Group for common elements
        common_group = QGroupBox("Common Elements")
        common_layout = QGridLayout(common_group)
        
        # Element color mappings
        elements = [
            ('Hydrogen', 'atom_h', '#d0d0d0'),
            ('Carbon', 'atom_c', '#55ff7f'),
            ('Nitrogen', 'atom_n', '#3050f8'),
            ('Oxygen', 'atom_o', '#ff0d0d'),
            ('Fluorine', 'atom_f', '#90e050'),
            ('Phosphorus', 'atom_p', '#ff8000'),
            ('Sulfur', 'atom_s', '#ffff30'),
        ]
        
        self.element_color_buttons = {}
        
        for i, (name, color_key, default_color) in enumerate(elements):
            row = i // 2
            col = (i % 2) * 3
            
            # Label
            label = QLabel(f"{name}:")
            label.setStyleSheet("color: #e8e8f0; font-weight: 500;")
            common_layout.addWidget(label, row, col)
            
            # Color button
            color_btn = QPushButton()
            color_btn.setFixedSize(40, 25)
            current_color = self.atom_colors.get(color_key, default_color)
            color_btn.setStyleSheet(f"background-color: {current_color}; border: 1px solid #666;")
            self.element_color_buttons[color_key] = (color_btn, label)
            color_btn.clicked.connect(lambda checked, k=color_key: self._choose_color(k))
            common_layout.addWidget(color_btn, row, col + 1)
            
            # Hex value label
            hex_label = QLabel(current_color.upper())
            hex_label.setStyleSheet("color: #9898b0; font-family: monospace; font-size: 11px;")
            common_layout.addWidget(hex_label, row, col + 2)
        
        layout.addWidget(common_group)
        
        # Group for halogens
        halogen_group = QGroupBox("Halogens")
        halogen_layout = QGridLayout(halogen_group)
        
        halogens = [
            ('Chlorine', 'atom_cl', '#1ff01f'),
            ('Bromine', 'atom_br', '#a62929'),
            ('Iodine', 'atom_i', '#940094'),
        ]
        
        for i, (name, color_key, default_color) in enumerate(halogens):
            row = i // 2
            col = (i % 2) * 3
            
            # Label
            label = QLabel(f"{name}:")
            label.setStyleSheet("color: #e8e8f0; font-weight: 500;")
            halogen_layout.addWidget(label, row, col)
            
            # Color button
            color_btn = QPushButton()
            color_btn.setFixedSize(40, 25)
            current_color = self.atom_colors.get(color_key, default_color)
            color_btn.setStyleSheet(f"background-color: {current_color}; border: 1px solid #666;")
            self.element_color_buttons[color_key] = (color_btn, label)
            color_btn.clicked.connect(lambda checked, k=color_key: self._choose_color(k))
            halogen_layout.addWidget(color_btn, row, col + 1)
            
            # Hex value label
            hex_label = QLabel(current_color.upper())
            hex_label.setStyleSheet("color: #9898b0; font-family: monospace; font-size: 11px;")
            halogen_layout.addWidget(hex_label, row, col + 2)
        
        layout.addWidget(halogen_group)
        layout.addStretch()
        
        return widget
        
    def _create_selection_colors_tab(self):
        """Create tab for selection and highlighting colors."""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        # Selection colors group
        selection_group = QGroupBox("Selection Colors")
        selection_layout = QGridLayout(selection_group)
        
        selection_colors = [
            ('Selected Atoms', 'atom_selected', '#ff00ff'),
            ('Highlighted Atoms', 'atom_highlight', '#ffff00'),
        ]
        
        self.selection_color_buttons = {}
        
        for i, (name, color_key, default_color) in enumerate(selection_colors):
            # Label
            label = QLabel(f"{name}:")
            label.setStyleSheet("color: #e8e8f0; font-weight: 500;")
            selection_layout.addWidget(label, i, 0)
            
            # Color button
            color_btn = QPushButton()
            color_btn.setFixedSize(40, 25)
            current_color = self.atom_colors.get(color_key, default_color)
            color_btn.setStyleSheet(f"background-color: {current_color}; border: 1px solid #666;")
            self.selection_color_buttons[color_key] = (color_btn, label)
            color_btn.clicked.connect(lambda checked, k=color_key: self._choose_color(k))
            selection_layout.addWidget(color_btn, i, 1)
            
            # Hex value label
            hex_label = QLabel(current_color.upper())
            hex_label.setStyleSheet("color: #9898b0; font-family: monospace; font-size: 11px;")
            selection_layout.addWidget(hex_label, i, 2)
        
        layout.addWidget(selection_group)
        layout.addStretch()
        
        return widget
        
    def _create_charge_colors_tab(self):
        """Create tab for charge-based coloring."""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        # Charge colors group
        charge_group = QGroupBox("Charge-Based Colors")
        charge_layout = QGridLayout(charge_group)
        
        charge_colors = [
            ('Positive Charges', 'atom_positive', '#0000ff'),
            ('Negative Charges', 'atom_negative', '#ff0000'),
        ]
        
        self.charge_color_buttons = {}
        
        for i, (name, color_key, default_color) in enumerate(charge_colors):
            # Label
            label = QLabel(f"{name}:")
            label.setStyleSheet("color: #e8e8f0; font-weight: 500;")
            charge_layout.addWidget(label, i, 0)
            
            # Color button
            color_btn = QPushButton()
            color_btn.setFixedSize(40, 25)
            current_color = self.atom_colors.get(color_key, default_color)
            color_btn.setStyleSheet(f"background-color: {current_color}; border: 1px solid #666;")
            self.charge_color_buttons[color_key] = (color_btn, label)
            color_btn.clicked.connect(lambda checked, k=color_key: self._choose_color(k))
            charge_layout.addWidget(color_btn, i, 1)
            
            # Hex value label
            hex_label = QLabel(current_color.upper())
            hex_label.setStyleSheet("color: #9898b0; font-family: monospace; font-size: 11px;")
            charge_layout.addWidget(hex_label, i, 2)
        
        layout.addWidget(charge_group)
        
        # Options for charge-based coloring
        options_group = QGroupBox("Charge Coloring Options")
        options_layout = QVBoxLayout(options_group)
        
        self.enable_charge_coloring = QCheckBox("Enable charge-based coloring")
        self.enable_charge_coloring.setChecked(False)
        options_layout.addWidget(self.enable_charge_coloring)
        
        # Charge threshold
        threshold_layout = QHBoxLayout()
        threshold_layout.addWidget(QLabel("Charge threshold:"))
        
        self.charge_threshold = QSpinBox()
        self.charge_threshold.setRange(1, 100)
        self.charge_threshold.setValue(10)
        self.charge_threshold.setSuffix(" %")
        threshold_layout.addWidget(self.charge_threshold)
        threshold_layout.addStretch()
        
        options_layout.addLayout(threshold_layout)
        layout.addWidget(options_group)
        layout.addStretch()
        
        return widget
        
    def _choose_color(self, color_key):
        """Open color dialog for choosing a color."""
        current_color = self.atom_colors.get(color_key, '#ffffff')
        color = QColorDialog.getColor(QColor(current_color), self, f"Choose Color for {color_key}")
        
        if color.isValid():
            hex_color = color.name().lower()
            self.atom_colors[color_key] = hex_color
            
            # Update button and label
            if color_key in self.element_color_buttons:
                btn, label = self.element_color_buttons[color_key]
            elif color_key in self.selection_color_buttons:
                btn, label = self.selection_color_buttons[color_key]
            elif color_key in self.charge_color_buttons:
                btn, label = self.charge_color_buttons[color_key]
            else:
                return
                
            btn.setStyleSheet(f"background-color: {hex_color}; border: 1px solid #666;")
            
            # Find and update hex label
            for i in range(btn.parent().layout().count()):
                widget = btn.parent().layout().itemAt(i).widget()
                if isinstance(widget, QLabel) and widget.text().startswith('#'):
                    widget.setText(hex_color.upper())
                    break
    
    def _reset_to_default(self):
        """Reset all colors to default values."""
        default_colors = {
            'atom_h': '#d0d0d0', 'atom_c': '#55ff7f', 'atom_n': '#3050f8',
            'atom_o': '#ff0d0d', 'atom_f': '#90e050', 'atom_p': '#ff8000',
            'atom_s': '#ffff30', 'atom_cl': '#1ff01f', 'atom_br': '#a62929',
            'atom_i': '#940094', 'atom_selected': '#ff00ff', 
            'atom_highlight': '#ffff00', 'atom_positive': '#0000ff',
            'atom_negative': '#ff0000'
        }
        
        self.atom_colors.update(default_colors)
        
        # Update all buttons and labels
        for color_key, hex_color in default_colors.items():
            if color_key in self.element_color_buttons:
                btn, label = self.element_color_buttons[color_key]
            elif color_key in self.selection_color_buttons:
                btn, label = self.selection_color_buttons[color_key]
            elif color_key in self.charge_color_buttons:
                btn, label = self.charge_color_buttons[color_key]
            else:
                continue
                
            btn.setStyleSheet(f"background-color: {hex_color}; border: 1px solid #666;")
            
            # Update hex label
            for i in range(btn.parent().layout().count()):
                widget = btn.parent().layout().itemAt(i).widget()
                if isinstance(widget, QLabel) and widget.text().startswith('#'):
                    widget.setText(hex_color.upper())
                    break
    
    def _apply_changes(self):
        """Apply the color changes and emit signal."""
        # Update the global COLORS dictionary
        COLORS.update(self.atom_colors)
        
        # Emit signal with new colors
        self.colors_changed.emit(self.atom_colors.copy())
        
        # Close dialog
        self.accept()
    
    def get_colors(self):
        """Get current color settings."""
        return self.atom_colors.copy()
    
    def set_colors(self, colors):
        """Set color settings from dictionary."""
        self.atom_colors.update(colors)
        
        # Update UI to reflect new colors
        for color_key, hex_color in colors.items():
            if color_key in self.element_color_buttons:
                btn, label = self.element_color_buttons[color_key]
            elif color_key in self.selection_color_buttons:
                btn, label = self.selection_color_buttons[color_key]
            elif color_key in self.charge_color_buttons:
                btn, label = self.charge_color_buttons[color_key]
            else:
                continue
                
            btn.setStyleSheet(f"background-color: {hex_color}; border: 1px solid #666;")
            
            # Update hex label
            for i in range(btn.parent().layout().count()):
                widget = btn.parent().layout().itemAt(i).widget()
                if isinstance(widget, QLabel) and widget.text().startswith('#'):
                    widget.setText(hex_color.upper())
                    break


def get_atom_color(element_symbol, custom_colors=None):
    """
    Get color for an atom element.
    
    Args:
        element_symbol: Chemical element symbol
        custom_colors: Optional custom color dictionary
        
    Returns:
        Hex color string
    """
    colors = custom_colors if custom_colors else COLORS
    
    # Try element-specific color first
    color_key = f'atom_{element_symbol.lower()}'
    if color_key in colors:
        return colors[color_key]
    
    # Fall back to default color
    return colors.get('atom_default', '#808080')


def get_charge_color(charge, threshold=0.1, custom_colors=None):
    """
    Get color based on partial charge.
    
    Args:
        charge: Partial charge value
        threshold: Charge magnitude threshold for coloring
        custom_colors: Optional custom color dictionary
        
    Returns:
        Hex color string
    """
    colors = custom_colors if custom_colors else COLORS
    
    if abs(charge) < threshold:
        return None  # Don't color small charges
    
    if charge > 0:
        return colors.get('atom_positive', '#0000ff')
    else:
        return colors.get('atom_negative', '#ff0000')
