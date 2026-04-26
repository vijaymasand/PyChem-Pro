"""
Protein Color Dialog — GUI for customizing protein secondary structure colors.

Provides PySide6-native color pickers for customizing helix, sheet, coil, and turn colors
used in protein cartoon and ribbon representations.
"""

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QColorDialog, QGroupBox, QGridLayout, QFrame
)
from PySide6.QtGui import QColor
from PySide6.QtCore import Qt

from src.shared.ui.theme import COLORS


class ProteinColorDialog(QDialog):
    """
    Dialog for customizing protein secondary structure colors.
    
    Allows users to customize colors for:
    - Alpha helices
    - Beta sheets  
    - Coils/loops
    - Turns
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Protein Color Settings")
        self.setMinimumWidth(400)
        self.setMinimumHeight(350)
        
        # Store current colors
        self.current_colors = {
            'ss_helix': COLORS.get('ss_helix', '#dc3232'),
            'ss_sheet': COLORS.get('ss_sheet', '#3296dc'),
            'ss_coil': COLORS.get('ss_coil', '#b4b4b4'),
            'ss_turn': COLORS.get('ss_turn', '#00d4aa'),
        }
        
        # Color preview widgets
        self.color_previews = {}
        
        self._init_ui()
        self._apply_styles()
    
    def _init_ui(self):
        """Initialize the dialog UI."""
        layout = QVBoxLayout(self)
        layout.setSpacing(16)
        layout.setContentsMargins(20, 20, 20, 20)
        
        # Title
        title = QLabel("Protein Secondary Structure Colors")
        title.setStyleSheet("font-size: 16px; font-weight: bold; color: #e8e8f0;")
        layout.addWidget(title)
        
        # Subtitle
        subtitle = QLabel("Customize colors for cartoon and ribbon representations")
        subtitle.setStyleSheet("font-size: 12px; color: #9898b0; margin-bottom: 10px;")
        layout.addWidget(subtitle)
        
        # Color settings group
        color_group = QGroupBox("Secondary Structure Colors")
        color_layout = QGridLayout(color_group)
        color_layout.setSpacing(12)
        color_layout.setContentsMargins(16, 20, 16, 16)
        
        # Helix color
        self._add_color_row(color_layout, 0, 'ss_helix', "Alpha Helix", 
                           "α-helices (H)", self.current_colors['ss_helix'])
        
        # Sheet color
        self._add_color_row(color_layout, 1, 'ss_sheet', "Beta Sheet", 
                           "β-sheets (E)", self.current_colors['ss_sheet'])
        
        # Coil color
        self._add_color_row(color_layout, 2, 'ss_coil', "Coil/Loop", 
                           "Coils and loops (C)", self.current_colors['ss_coil'])
        
        # Turn color
        self._add_color_row(color_layout, 3, 'ss_turn', "Turn", 
                           "Turns (T)", self.current_colors['ss_turn'])
        
        layout.addWidget(color_group)
        
        # Presets section
        presets_group = QGroupBox("Color Presets")
        presets_layout = QHBoxLayout(presets_group)
        presets_layout.setSpacing(8)
        
        # PyMOL preset
        pymol_btn = QPushButton("PyMOL Default")
        pymol_btn.setToolTip("Red helices, blue sheets (classic PyMOL)")
        pymol_btn.clicked.connect(self._apply_pymol_preset)
        presets_layout.addWidget(pymol_btn)
        
        # Jmol preset
        jmol_btn = QPushButton("Jmol Default")
        jmol_btn.setToolTip("Cyan helices, yellow sheets (classic Jmol)")
        jmol_btn.clicked.connect(self._apply_jmol_preset)
        presets_layout.addWidget(jmol_btn)
        
        # Pastel preset
        pastel_btn = QPushButton("Pastel")
        pastel_btn.setToolTip("Soft pastel colors")
        pastel_btn.clicked.connect(self._apply_pastel_preset)
        presets_layout.addWidget(pastel_btn)
        
        layout.addWidget(presets_group)
        
        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setStyleSheet("background-color: #2a3a5c; margin: 10px 0;")
        layout.addWidget(separator)
        
        # Buttons
        button_layout = QHBoxLayout()
        button_layout.setSpacing(10)
        button_layout.addStretch()
        
        self.apply_btn = QPushButton("Apply")
        self.apply_btn.setObjectName("btnSuccess")
        self.apply_btn.setStyleSheet(self._get_success_button_style())
        self.apply_btn.clicked.connect(self._apply_colors)
        button_layout.addWidget(self.apply_btn)
        
        close_btn = QPushButton("Close")
        close_btn.setObjectName("btnSecondary")
        close_btn.clicked.connect(self.reject)
        button_layout.addWidget(close_btn)
        
        layout.addLayout(button_layout)
    
    def _add_color_row(self, layout, row, color_key, label, description, initial_color):
        """Add a color selection row to the layout."""
        # Label
        label_widget = QLabel(label)
        label_widget.setStyleSheet("font-weight: 600; color: #e8e8f0;")
        layout.addWidget(label_widget, row, 0)
        
        # Description
        desc_widget = QLabel(description)
        desc_widget.setStyleSheet("color: #9898b0; font-size: 11px;")
        layout.addWidget(desc_widget, row, 1)
        
        # Color preview button
        preview_btn = QPushButton()
        preview_btn.setFixedSize(40, 28)
        preview_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {initial_color};
                border: 2px solid #2a3a5c;
                border-radius: 4px;
            }}
            QPushButton:hover {{
                border-color: #6c63ff;
            }}
        """)
        preview_btn.clicked.connect(lambda: self._pick_color(color_key, preview_btn))
        layout.addWidget(preview_btn, row, 2)
        
        # Change button
        change_btn = QPushButton("Change")
        change_btn.setStyleSheet(self._get_secondary_button_style())
        change_btn.clicked.connect(lambda: self._pick_color(color_key, preview_btn))
        layout.addWidget(change_btn, row, 3)
        
        # Store reference
        self.color_previews[color_key] = preview_btn
    
    def _pick_color(self, color_key, preview_btn):
        """Open color picker and update color."""
        current_color = QColor(self.current_colors[color_key])
        
        color = QColorDialog.getColor(
            current_color,
            self,
            f"Select {color_key.replace('ss_', '').title()} Color",
            QColorDialog.ColorDialogOption.ShowAlphaChannel
        )
        
        if color.isValid():
            hex_color = color.name()
            self.current_colors[color_key] = hex_color
            
            # Update preview button
            preview_btn.setStyleSheet(f"""
                QPushButton {{
                    background-color: {hex_color};
                    border: 2px solid #2a3a5c;
                    border-radius: 4px;
                }}
                QPushButton:hover {{
                    border-color: #6c63ff;
                }}
            """)
    
    def _apply_pymol_preset(self):
        """Apply PyMOL default colors."""
        self.current_colors = {
            'ss_helix': '#dc3232',  # Red
            'ss_sheet': '#3296dc',  # Blue
            'ss_coil': '#b4b4b4',   # Gray
            'ss_turn': '#00d4aa',   # Teal
        }
        self._update_previews()
    
    def _apply_jmol_preset(self):
        """Apply Jmol default colors."""
        self.current_colors = {
            'ss_helix': '#00ffff',  # Cyan
            'ss_sheet': '#ffff00',  # Yellow
            'ss_coil': '#b4b4b4',   # Gray
            'ss_turn': '#ff69b4',   # Pink
        }
        self._update_previews()
    
    def _apply_pastel_preset(self):
        """Apply pastel colors."""
        self.current_colors = {
            'ss_helix': '#ff9999',  # Light red
            'ss_sheet': '#99ccff',  # Light blue
            'ss_coil': '#cccccc',   # Light gray
            'ss_turn': '#99ffcc',   # Light teal
        }
        self._update_previews()
    
    def _update_previews(self):
        """Update all color preview buttons."""
        for color_key, preview_btn in self.color_previews.items():
            hex_color = self.current_colors[color_key]
            preview_btn.setStyleSheet(f"""
                QPushButton {{
                    background-color: {hex_color};
                    border: 2px solid #2a3a5c;
                    border-radius: 4px;
                }}
                QPushButton:hover {{
                    border-color: #6c63ff;
                }}
            """)
    def _apply_colors(self):
        """Apply selected colors to theme but keep dialog open."""
        # Update theme colors - ensure we update the actual theme COLORS
        from src.shared.ui.theme import COLORS as THEME_COLORS
        for key, value in self.current_colors.items():
            THEME_COLORS[key] = value
            
        # CRITICAL: Invalidate both the mesh cache and image cache so the new colors take effect
        try:
            from src.features.visualization_3d.services.protein_rendering import _cartoon_gen, render_protein_cartoon
            _cartoon_gen.invalidate()
            if hasattr(render_protein_cartoon, '_img_cache'):
                render_protein_cartoon._img_cache.clear()
        except ImportError:
            pass
        
        # Determine if we have a parent and trigger repaint
        parent = self.parent()
        if parent:
            if hasattr(parent, 'viewer_3d'):
                if hasattr(parent.viewer_3d, 'repaint'):
                    parent.viewer_3d.repaint()
                if hasattr(parent.viewer_3d, 'update'):
                    parent.viewer_3d.update()
            
            if hasattr(parent, 'status_bar'):
                parent.status_bar.showMessage("Applied protein colors.")
                
        # Force Qt event processing so it repaints immediately even when dialog is open
        from PySide6.QtCore import QCoreApplication
        QCoreApplication.processEvents()
        # No self.accept() here! We just want to apply without closing.
        return self.current_colors.copy()
        
    def get_colors(self):
        """Return the current colors dict."""
        return self.current_colors.copy()
    
    def _apply_styles(self):
        """Apply dark theme styles to dialog."""
        self.setStyleSheet(f"""
            QDialog {{
                background-color: {COLORS['bg_primary']};
                color: {COLORS['text_primary']};
            }}
            QGroupBox {{
                background-color: {COLORS['bg_tertiary']};
                border: 1px solid {COLORS['border']};
                border-radius: 8px;
                margin-top: 12px;
                padding: 12px 16px 16px 16px;
                font-weight: 600;
                font-size: 13px;
                color: {COLORS['accent2']};
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 6px;
            }}
            QLabel {{
                color: {COLORS['text_primary']};
            }}
        """)
    
    def _get_success_button_style(self):
        """Get style for success/apply button."""
        return f"""
            QPushButton {{
                background-color: {COLORS['success']};
                color: white;
                border: none;
                padding: 10px 24px;
                border-radius: 6px;
                font-weight: 600;
                font-size: 13px;
            }}
            QPushButton:hover {{
                background-color: #00e8a0;
            }}
        """
    
    def _get_secondary_button_style(self):
        """Get style for secondary buttons."""
        return f"""
            QPushButton {{
                background-color: {COLORS['bg_widget']};
                color: {COLORS['text_primary']};
                border: 1px solid {COLORS['border']};
                padding: 8px 16px;
                border-radius: 6px;
                font-size: 12px;
            }}
            QPushButton:hover {{
                background-color: {COLORS['bg_hover']};
                border-color: {COLORS['accent']};
            }}
        """


def show_protein_color_dialog(parent=None):
    """
    Show the protein color dialog and return selected colors.
    
    Args:
        parent: Parent widget for the dialog
        
    Returns:
        dict: Selected colors if accepted, None if cancelled
    """
    dialog = ProteinColorDialog(parent)
    
    if dialog.exec() == QDialog.DialogCode.Accepted:
        return dialog.get_colors()
    
    return None
