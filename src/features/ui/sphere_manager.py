"""
Sphere Management Panel for managing dummy spheres.

Provides a panel to list, remove, and manage all created spheres.
"""

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QListWidget, 
    QPushButton, QLabel, QGroupBox, QCheckBox, QSpinBox,
    QColorDialog, QInputDialog, QMessageBox, QSplitter
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor

from src.shared.ui.theme import COLORS


class SphereManagerPanel(QWidget):
    """
    Panel for managing dummy spheres in the 3D viewer.
    """
    
    # Signals
    sphere_selected = Signal(str)  # Emitted when a sphere is selected
    sphere_removed = Signal(str)    # Emitted when a sphere is removed
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.sphere_manager = None
        self.current_molecule = None
        self.sphere_list = []
        
        self._init_ui()
        
    def _init_ui(self):
        """Initialize the user interface."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)
        
        # Title
        title = QLabel("Sphere Manager")
        title.setStyleSheet(f"font-size: 14px; font-weight: bold; color: {COLORS['text_primary']};")
        layout.addWidget(title)
        
        # Sphere list
        list_group = QGroupBox("Created Spheres")
        list_layout = QVBoxLayout(list_group)
        
        self.sphere_list_widget = QListWidget()
        self.sphere_list_widget.setStyleSheet(f"""
            QListWidget {{
                background-color: {COLORS['bg_widget']};
                color: {COLORS['text_primary']};
                border: 1px solid {COLORS['border']};
                border-radius: 4px;
                padding: 4px;
            }}
            QListWidget::item {{
                padding: 4px;
                margin: 2px;
                border-radius: 3px;
            }}
            QListWidget::item:selected {{
                background-color: {COLORS['accent']};
                color: white;
            }}
            QListWidget::item:hover {{
                background-color: {COLORS['bg_hover']};
            }}
        """)
        self.sphere_list_widget.itemClicked.connect(self._on_sphere_selected)
        
        list_layout.addWidget(self.sphere_list_widget)
        layout.addWidget(list_group)
        
        # Sphere actions
        actions_group = QGroupBox("Sphere Actions")
        actions_layout = QVBoxLayout(actions_group)
        
        # Remove selected sphere
        self.remove_btn = QPushButton("Remove Selected")
        self.remove_btn.clicked.connect(self._remove_selected_sphere)
        self.remove_btn.setEnabled(False)
        actions_layout.addWidget(self.remove_btn)
        
        # Remove all spheres
        self.clear_all_btn = QPushButton("Clear All Spheres")
        self.clear_all_btn.clicked.connect(self._clear_all_spheres)
        self.clear_all_btn.setEnabled(False)
        actions_layout.addWidget(self.clear_all_btn)
        
        layout.addWidget(actions_group)
        
        # Sphere properties
        props_group = QGroupBox("Sphere Properties")
        props_layout = QVBoxLayout(props_group)
        
        # Visibility checkbox
        self.visible_checkbox = QCheckBox("Show Spheres")
        self.visible_checkbox.setChecked(True)
        self.visible_checkbox.stateChanged.connect(self._toggle_visibility)
        props_layout.addWidget(self.visible_checkbox)
        
        # Size adjustment
        size_layout = QHBoxLayout()
        size_layout.addWidget(QLabel("Default Size:"))
        self.size_spinbox = QSpinBox()
        self.size_spinbox.setRange(1, 50)
        self.size_spinbox.setValue(30)  # 0.3 units
        self.size_spinbox.setSuffix(" × 0.1")
        self.size_spinbox.valueChanged.connect(self._update_default_size)
        size_layout.addWidget(self.size_spinbox)
        size_layout.addStretch()
        props_layout.addLayout(size_layout)
        
        layout.addWidget(props_group)
        
        # Statistics
        stats_group = QGroupBox("Statistics")
        stats_layout = QVBoxLayout(stats_group)
        
        self.stats_label = QLabel("No spheres created")
        self.stats_label.setStyleSheet(f"color: {COLORS['text_secondary']};")
        stats_layout.addWidget(self.stats_label)
        
        layout.addWidget(stats_group)
        
        layout.addStretch()
        
    def set_molecule(self, molecule):
        """Set the current molecule and update sphere list."""
        self.current_molecule = molecule
        
        if molecule:
            try:
                from src.features.visualization_3d.services.dummy_sphere import DummySphereManager
                self.sphere_manager = DummySphereManager(molecule)
                self._update_sphere_list()
            except Exception as e:
                print(f"Error setting molecule: {e}")
        else:
            self.sphere_manager = None
            self._clear_sphere_list()
    
    def _update_sphere_list(self):
        """Update the sphere list display."""
        if not self.sphere_manager:
            return
            
        self.sphere_list_widget.clear()
        self.sphere_list.clear()
        
        spheres = self.sphere_manager.get_all_spheres()
        
        for sphere in spheres:
            # Create list item with sphere info
            item_text = f"{sphere.label} - {sphere.get_position()[0]:.2f}, {sphere.get_position()[1]:.2f}, {sphere.get_position()[2]:.2f}"
            
            # Add color indicator
            color_indicator = f"● {sphere.color}"
            full_text = f"{color_indicator} {item_text}"
            
            self.sphere_list_widget.addItem(full_text)
            self.sphere_list.append(sphere)
        
        self._update_statistics()
        self._update_button_states()
    
    def _clear_sphere_list(self):
        """Clear the sphere list display."""
        self.sphere_list_widget.clear()
        self.sphere_list.clear()
        self._update_statistics()
        self._update_button_states()
    
    def _on_sphere_selected(self, item):
        """Handle sphere selection in the list."""
        index = self.sphere_list_widget.row(item)
        if 0 <= index < len(self.sphere_list):
            sphere = self.sphere_list[index]
            self.sphere_selected.emit(sphere.sphere_id)
            self.remove_btn.setEnabled(True)
    
    def _remove_selected_sphere(self):
        """Remove the selected sphere."""
        current_row = self.sphere_list_widget.currentRow()
        if 0 <= current_row < len(self.sphere_list):
            sphere = self.sphere_list[current_row]
            
            reply = QMessageBox.question(
                self, "Remove Sphere", 
                f"Remove sphere '{sphere.label}'?",
                QMessageBox.Yes | QMessageBox.No
            )
            
            if reply == QMessageBox.Yes:
                if self.sphere_manager.remove_sphere(sphere.sphere_id):
                    self._update_sphere_list()
                    self.sphere_removed.emit(sphere.sphere_id)
                    self.remove_btn.setEnabled(False)
    
    def _clear_all_spheres(self):
        """Clear all spheres."""
        reply = QMessageBox.question(
            self, "Clear All Spheres",
            "Remove all spheres? This action cannot be undone.",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            if self.sphere_manager:
                self.sphere_manager.clear_all_spheres()
                self._update_sphere_list()
    
    def _toggle_visibility(self, state):
        """Toggle sphere visibility."""
        # This would need to be implemented in the 3D viewer
        # For now, just update the checkbox state
        pass
    
    def _update_default_size(self, value):
        """Update default sphere size."""
        # This would affect future sphere creation
        pass
    
    def _update_statistics(self):
        """Update statistics display."""
        if self.sphere_manager:
            summary = self.sphere_manager.get_sphere_summary()
            total = summary['total_spheres']
            
            if total == 0:
                self.stats_label.setText("No spheres created")
            else:
                self.stats_label.setText(f"Total spheres: {total}")
        else:
            self.stats_label.setText("No molecule loaded")
    
    def _update_button_states(self):
        """Update button enabled states."""
        has_spheres = len(self.sphere_list) > 0
        has_selection = self.sphere_list_widget.currentItem() is not None
        
        self.remove_btn.setEnabled(has_selection)
        self.clear_all_btn.setEnabled(has_spheres)
    
    def add_sphere(self, sphere_id, position, radius, color, label):
        """Add a new sphere to the list."""
        if self.sphere_manager:
            self._update_sphere_list()
    
    def refresh(self):
        """Refresh the sphere list."""
        self._update_sphere_list()


class SphereDockWidget(QWidget):
    """
    Dockable widget for sphere management.
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Sphere Manager")
        self.resize(300, 400)
        
        layout = QVBoxLayout(self)
        
        # Create sphere manager panel
        self.sphere_panel = SphereManagerPanel()
        layout.addWidget(self.sphere_panel)
        
        # Quick actions
        quick_actions = QGroupBox("Quick Actions")
        quick_layout = QVBoxLayout(quick_actions)
        
        # Add COM sphere
        com_btn = QPushButton("Add COM Sphere")
        com_btn.clicked.connect(self._add_com_sphere)
        quick_layout.addWidget(com_btn)
        
        # Add centroid sphere
        centroid_btn = QPushButton("Add Centroid")
        centroid_btn.clicked.connect(self._add_centroid_sphere)
        quick_layout.addWidget(centroid_btn)
        
        layout.addWidget(quick_actions)
    
    def set_molecule(self, molecule):
        """Set the current molecule."""
        self.sphere_panel.set_molecule(molecule)
    
    def _add_com_sphere(self):
        """Add COM sphere."""
        # This would call the main window's COM sphere method
        if self.parent():
            getattr(self.parent(), '_add_com_sphere', lambda: None)()
    
    def _add_centroid_sphere(self):
        """Add centroid sphere."""
        # This would call the main window's centroid sphere method
        if self.parent():
            getattr(self.parent(), '_add_centroid_sphere', lambda: None)()
