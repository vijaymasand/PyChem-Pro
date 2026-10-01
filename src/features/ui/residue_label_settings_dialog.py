"""
Residue Label Settings Dialog — GUI for customizing residue label appearance.

Provides options for customizing residue label color, font size, and other visual properties
used in molecular visualization.
"""

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
    QPushButton, QColorDialog, QSpinBox, QGroupBox,
    QFormLayout, QCheckBox, QComboBox
)
from PySide6.QtGui import QFont, QFontDatabase
from PySide6.QtCore import Qt, Signal


class ResidueLabelSettingsDialog(QDialog):
    """Dialog for customizing residue label appearance."""
    
    settings_changed = Signal(dict)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Residue Label Settings")
        self.setMinimumWidth(400)
        
        self._init_ui()
        self._load_current_settings()
    
    def _init_ui(self):
        """Initialize the user interface."""
        layout = QVBoxLayout(self)
        
        # Color selection
        color_group = QGroupBox("Label Color")
        color_layout = QHBoxLayout()
        
        self.color_label = QLabel("Color:")
        self.color_button = QPushButton("    ")
        self.color_button.setFixedWidth(60)
        self.color_button.clicked.connect(self._choose_color)
        
        color_layout.addWidget(self.color_label)
        color_layout.addWidget(self.color_button)
        color_layout.addStretch()
        color_group.setLayout(color_layout)
        
        # Font settings
        font_group = QGroupBox("Font Settings")
        font_layout = QFormLayout()
        
        # Font family
        self.font_combo = QComboBox()
        font_db = QFontDatabase()
        font_families = font_db.families()
        self.font_combo.addItems(font_families)
        font_layout.addRow("Font Family:", self.font_combo)
        
        # Font size
        self.font_size_spin = QSpinBox()
        self.font_size_spin.setRange(6, 48)
        self.font_size_spin.setValue(12)
        font_layout.addRow("Font Size:", self.font_size_spin)
        
        # Font style
        self.bold_check = QCheckBox("Bold")
        self.italic_check = QCheckBox("Italic")
        font_layout.addRow("Style:", self.bold_check)
        font_layout.addRow("", self.italic_check)
        
        font_group.setLayout(font_layout)
        
        # Label options
        options_group = QGroupBox("Label Options")
        options_layout = QFormLayout()
        
        self.show_labels_check = QCheckBox("Show Residue Labels")
        self.show_labels_check.setChecked(True)
        options_layout.addRow("", self.show_labels_check)
        
        self.background_check = QCheckBox("Label Background")
        self.background_check.setChecked(False)
        options_layout.addRow("", self.background_check)
        
        options_group.setLayout(options_layout)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        self.reset_btn = QPushButton("Reset to Default")
        self.reset_btn.clicked.connect(self._reset_to_default)
        
        self.apply_btn = QPushButton("Apply")
        self.apply_btn.clicked.connect(self._apply_settings)
        self.apply_btn.setDefault(True)
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.reject)
        
        button_layout.addWidget(self.reset_btn)
        button_layout.addStretch()
        button_layout.addWidget(self.apply_btn)
        button_layout.addWidget(self.cancel_btn)
        
        # Add all groups to main layout
        layout.addWidget(color_group)
        layout.addWidget(font_group)
        layout.addWidget(options_group)
        layout.addSpacing(10)
        layout.addLayout(button_layout)
    
    def _load_current_settings(self):
        """Load current settings from viewer."""
        if hasattr(self.parent(), 'viewer_3d'):
            viewer = self.parent().viewer_3d
            
            # Get current label color
            if hasattr(viewer, 'labeled_residues') and viewer.labeled_residues:
                # Use first available color as current
                first_color = next(iter(viewer.labeled_residues.values()), None)
                if first_color:
                    self.current_color = first_color
                    self.color_button.setStyleSheet(f"background-color: {first_color.name()};")
            
            # Load font and display settings if they exist
            if hasattr(viewer, 'residue_label_settings'):
                settings = viewer.residue_label_settings
                if 'font_family' in settings:
                    self.font_combo.setCurrentText(settings['font_family'])
                else:
                    self.font_combo.setCurrentText("Arial")
                if 'font_size' in settings:
                    self.font_size_spin.setValue(settings['font_size'])
                else:
                    self.font_size_spin.setValue(12)
                if 'bold' in settings:
                    self.bold_check.setChecked(settings['bold'])
                if 'italic' in settings:
                    self.italic_check.setChecked(settings['italic'])
                if 'show_labels' in settings:
                    self.show_labels_check.setChecked(settings['show_labels'])
                if 'background' in settings:
                    self.background_check.setChecked(settings['background'])
            else:
                # Set default font values
                self.font_combo.setCurrentText("Arial")
                self.font_size_spin.setValue(12)
                self.bold_check.setChecked(False)
                self.italic_check.setChecked(False)
    
    def _choose_color(self):
        """Open color picker dialog."""
        color = QColorDialog.getColor(self.current_color if hasattr(self, 'current_color') else Qt.black, 
                                     self, "Select Residue Label Color")
        if color.isValid():
            self.current_color = color
            self.color_button.setStyleSheet(f"background-color: {color.name()};")
    
    def _reset_to_default(self):
        """Reset settings to default values."""
        self.current_color = Qt.black
        self.color_button.setStyleSheet("background-color: black;")
        self.font_combo.setCurrentText("Arial")
        self.font_size_spin.setValue(12)
        self.bold_check.setChecked(False)
        self.italic_check.setChecked(False)
        self.show_labels_check.setChecked(True)
        self.background_check.setChecked(False)
    
    def _apply_settings(self):
        """Apply settings and emit signal."""
        settings = {
            'color': self.current_color if hasattr(self, 'current_color') else Qt.black,
            'font_family': self.font_combo.currentText(),
            'font_size': self.font_size_spin.value(),
            'bold': self.bold_check.isChecked(),
            'italic': self.italic_check.isChecked(),
            'show_labels': self.show_labels_check.isChecked(),
            'background': self.background_check.isChecked()
        }
        
        self.settings_changed.emit(settings)
        self.accept()
    
    def get_settings(self):
        """Get current settings."""
        return {
            'color': self.current_color if hasattr(self, 'current_color') else Qt.black,
            'font_family': self.font_combo.currentText(),
            'font_size': self.font_size_spin.value(),
            'bold': self.bold_check.isChecked(),
            'italic': self.italic_check.isChecked(),
            'show_labels': self.show_labels_check.isChecked(),
            'background': self.background_check.isChecked()
        }


def show_residue_label_settings_dialog(parent=None) -> dict:
    """Show residue label settings dialog and return selected settings."""
    dialog = ResidueLabelSettingsDialog(parent)
    
    if dialog.exec() == QDialog.Accepted:
        return dialog.get_settings()
    
    return None
