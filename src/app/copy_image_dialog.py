"""
Copy as Image Dialog — Dialog for copying molecular views to clipboard.

Provides options for format (PNG, JPEG, BMP), DPI, and whether to
render selected atoms only or the entire molecule.
"""

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QSpinBox, QCheckBox, QPushButton, QDialogButtonBox,
    QGroupBox, QRadioButton, QWidget
)
from PySide6.QtCore import Qt


class CopyImageDialog(QDialog):
    """
    Dialog for configuring image copy options.
    
    Allows user to select:
    - Source: 2D view, 3D view, or Both
    - Format: PNG, JPEG, BMP
    - DPI: Resolution for the output image
    - Selection only: Copy only selected atoms or full molecule
    """
    
    # Format mapping to file extensions
    FORMATS = {
        "PNG": "png",
        "JPEG": "jpg",
        "BMP": "bmp"
    }
    
    def __init__(self, parent=None, has_selection_2d=False, has_selection_3d=False):
        super().__init__(parent)
        self.setWindowTitle("Copy as Image")
        self.setMinimumWidth(350)
        
        self._has_selection_2d = has_selection_2d
        self._has_selection_3d = has_selection_3d
        
        self._build_ui()
        
    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(12)
        
        # Source selection
        source_group = QGroupBox("Source")
        source_layout = QVBoxLayout(source_group)
        
        self._2d_radio = QRadioButton("2D View")
        self._2d_radio.setChecked(True)
        source_layout.addWidget(self._2d_radio)
        
        self._3d_radio = QRadioButton("3D View")
        source_layout.addWidget(self._3d_radio)
        
        self._both_radio = QRadioButton("Both Views (side by side)")
        source_layout.addWidget(self._both_radio)
        
        layout.addWidget(source_group)
        
        # Format selection
        format_layout = QHBoxLayout()
        format_layout.addWidget(QLabel("Format:"))
        self._format_combo = QComboBox()
        self._format_combo.addItems(list(self.FORMATS.keys()))
        self._format_combo.setCurrentText("PNG")
        format_layout.addWidget(self._format_combo)
        layout.addLayout(format_layout)
        
        # DPI selection
        dpi_layout = QHBoxLayout()
        dpi_layout.addWidget(QLabel("DPI (Resolution):"))
        self._dpi_spin = QSpinBox()
        self._dpi_spin.setRange(72, 1200)
        self._dpi_spin.setValue(300)
        self._dpi_spin.setSingleStep(50)
        dpi_layout.addWidget(self._dpi_spin)
        layout.addLayout(dpi_layout)
        
        # Selection only option
        self._selection_only_check = QCheckBox("Copy selection only (if nothing selected, copies full molecule)")
        self._selection_only_check.setChecked(True)
        self._selection_only_check.setEnabled(self._has_selection_2d or self._has_selection_3d)
        layout.addWidget(self._selection_only_check)
        
        # High Quality / Ray Tracing option
        self._hq_check = QCheckBox("Ray (slow) - High Anti-Aliasing (3D view only)")
        self._hq_check.setChecked(True)
        self._hq_check.setEnabled(self._has_selection_3d or True) # Always enable since they can select 3D view
        layout.addWidget(self._hq_check)
        
        # Status label
        self._status_label = QLabel()
        self._status_label.setStyleSheet("color: gray; font-size: 11px;")
        layout.addWidget(self._status_label)
        self._update_status()
        
        # Connect source changes to update status
        self._2d_radio.toggled.connect(self._update_status)
        self._3d_radio.toggled.connect(self._update_status)
        self._both_radio.toggled.connect(self._update_status)
        
        # Buttons
        button_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        
    def _update_status(self):
        """Update the status label based on current selection."""
        if self._2d_radio.isChecked():
            if self._has_selection_2d:
                self._status_label.setText("2D View: Selection available")
            else:
                self._status_label.setText("2D View: No selection (will copy full molecule)")
        elif self._3d_radio.isChecked():
            if self._has_selection_3d:
                self._status_label.setText("3D View: Selection available")
            else:
                self._status_label.setText("3D View: No selection (will copy full molecule)")
        else:
            self._status_label.setText("Both views will be combined side by side")
    
    def get_settings(self):
        """
        Get the current dialog settings.
        
        Returns:
            dict with keys: source, format, dpi, selection_only
        """
        if self._2d_radio.isChecked():
            source = "2d"
        elif self._3d_radio.isChecked():
            source = "3d"
        else:
            source = "both"
            
        return {
            "source": source,
            "format": self.FORMATS[self._format_combo.currentText()],
            "dpi": self._dpi_spin.value(),
            "selection_only": self._selection_only_check.isChecked(),
            "high_quality": self._hq_check.isChecked()
        }
