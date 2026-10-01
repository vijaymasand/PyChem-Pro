"""
Installed Plugins Dialog - PySide6 GUI for selecting and managing plugins.

Provides a dialog window to view, select, and activate/deactivate plugins.
Similar to PyMOL's plugin management interface.
"""

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QListWidget, QListWidgetItem,
    QPushButton, QLabel, QCheckBox, QGroupBox, QTextEdit, QSplitter,
    QMessageBox, QMenu, QFrame, QWidget, QSizePolicy
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QAction

from src.shared.ui.theme import COLORS


class InstalledPluginsDialog(QDialog):
    """
    Dialog for managing installed plugins.
    
    Allows users to:
    - View all installed plugins
    - Enable/disable individual plugins
    - Select active plugins for the toolbar
    - View plugin information
    """
    
    # Signal emitted when plugin selection changes
    plugins_changed = Signal(list)  # List of active plugin names
    
    def __init__(self, plugin_interface=None, parent=None):
        super().__init__(parent)
        self.plugin_interface = plugin_interface
        self.setWindowTitle("Installed Plugins")
        self.setMinimumSize(600, 450)
        self.resize(700, 500)
        
        # Store plugin items
        self.plugin_items = {}
        self.active_plugins = []
        
        self._init_ui()
        self._apply_styles()
        self._load_plugins()
    
    def _init_ui(self):
        """Initialize the dialog UI."""
        layout = QVBoxLayout(self)
        layout.setSpacing(12)
        layout.setContentsMargins(16, 16, 16, 16)
        
        # Title
        title_label = QLabel("Installed Plugins")
        title_label.setStyleSheet(f"""
            font-size: 18px;
            font-weight: bold;
            color: {COLORS['text_primary']};
            margin-bottom: 8px;
        """)
        layout.addWidget(title_label)
        
        # Description
        desc_label = QLabel("Select plugins to activate in the toolbar. Enabled plugins will add buttons to the main toolbar for quick access.")
        desc_label.setStyleSheet(f"color: {COLORS['text_secondary']}; font-size: 12px;")
        desc_label.setWordWrap(True)
        layout.addWidget(desc_label)
        
        layout.addSpacing(8)
        
        # Main splitter
        splitter = QSplitter(Qt.Orientation.Horizontal)
        
        # Left side: Plugin list
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(8)
        
        # Plugin list header
        list_header = QLabel("Available Plugins")
        list_header.setStyleSheet(f"""
            font-size: 13px;
            font-weight: 600;
            color: {COLORS['accent2']};
        """)
        left_layout.addWidget(list_header)
        
        # Plugin list widget
        self.plugin_list = QListWidget()
        self.plugin_list.setSelectionMode(QListWidget.SelectionMode.SingleSelection)
        self.plugin_list.currentItemChanged.connect(self._on_plugin_selected)
        left_layout.addWidget(self.plugin_list)
        
        # Select all/none buttons
        btn_layout = QHBoxLayout()
        
        self.select_all_btn = QPushButton("Select All")
        self.select_all_btn.clicked.connect(self._select_all_plugins)
        btn_layout.addWidget(self.select_all_btn)
        
        self.select_none_btn = QPushButton("Select None")
        self.select_none_btn.clicked.connect(self._select_no_plugins)
        btn_layout.addWidget(self.select_none_btn)
        
        btn_layout.addStretch()
        left_layout.addLayout(btn_layout)
        
        splitter.addWidget(left_widget)
        
        # Right side: Plugin details
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(8)
        
        # Details header
        details_header = QLabel("Plugin Details")
        details_header.setStyleSheet(f"""
            font-size: 13px;
            font-weight: 600;
            color: {COLORS['accent2']};
        """)
        right_layout.addWidget(details_header)
        
        # Details group box
        details_group = QGroupBox()
        details_group.setStyleSheet(f"""
            QGroupBox {{
                background-color: {COLORS['bg_tertiary']};
                border: 1px solid {COLORS['border']};
                border-radius: 6px;
                padding: 12px;
            }}
        """)
        details_layout = QVBoxLayout(details_group)
        
        # Plugin name
        self.name_label = QLabel("No plugin selected")
        self.name_label.setStyleSheet(f"""
            font-size: 14px;
            font-weight: bold;
            color: {COLORS['text_primary']};
        """)
        details_layout.addWidget(self.name_label)
        
        # Separator
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet(f"background-color: {COLORS['border']};")
        details_layout.addWidget(sep)
        
        # Enable checkbox
        self.enable_check = QCheckBox("Enable plugin in toolbar")
        self.enable_check.setStyleSheet(f"color: {COLORS['text_primary']};")
        self.enable_check.stateChanged.connect(self._on_enable_changed)
        details_layout.addWidget(self.enable_check)
        
        details_layout.addSpacing(8)
        
        # Description
        desc_title = QLabel("Description:")
        desc_title.setStyleSheet(f"font-weight: 600; color: {COLORS['text_secondary']};")
        details_layout.addWidget(desc_title)
        
        self.description_text = QTextEdit()
        self.description_text.setReadOnly(True)
        self.description_text.setMaximumHeight(80)
        self.description_text.setStyleSheet(f"""
            QTextEdit {{
                background-color: {COLORS['bg_secondary']};
                color: {COLORS['text_primary']};
                border: 1px solid {COLORS['border']};
                border-radius: 4px;
                padding: 4px;
            }}
        """)
        details_layout.addWidget(self.description_text)
        
        # Info fields
        self.version_label = QLabel("Version: -")
        self.version_label.setStyleSheet(f"color: {COLORS['text_secondary']}; font-size: 11px;")
        details_layout.addWidget(self.version_label)
        
        self.author_label = QLabel("Author: -")
        self.author_label.setStyleSheet(f"color: {COLORS['text_secondary']}; font-size: 11px;")
        details_layout.addWidget(self.author_label)
        
        self.status_label = QLabel("Status: Not loaded")
        self.status_label.setStyleSheet(f"color: {COLORS['text_secondary']}; font-size: 11px;")
        details_layout.addWidget(self.status_label)
        
        details_layout.addStretch()
        right_layout.addWidget(details_group)
        
        # Action buttons
        action_layout = QHBoxLayout()
        
        self.load_btn = QPushButton("Load Plugin")
        self.load_btn.clicked.connect(self._load_selected_plugin)
        action_layout.addWidget(self.load_btn)
        
        self.unload_btn = QPushButton("Unload Plugin")
        self.unload_btn.clicked.connect(self._unload_selected_plugin)
        action_layout.addWidget(self.unload_btn)
        
        action_layout.addStretch()
        right_layout.addLayout(action_layout)
        
        splitter.addWidget(right_widget)
        
        # Set splitter sizes (40% left, 60% right)
        splitter.setSizes([280, 420])
        layout.addWidget(splitter, 1)
        
        # Bottom buttons
        bottom_layout = QHBoxLayout()
        bottom_layout.addStretch()
        
        self.apply_btn = QPushButton("Apply")
        self.apply_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {COLORS['success']};
                color: white;
                border: none;
                padding: 8px 24px;
                border-radius: 4px;
                font-weight: 600;
            }}
            QPushButton:hover {{
                background-color: #00e8a0;
            }}
        """)
        self.apply_btn.clicked.connect(self._apply_changes)
        bottom_layout.addWidget(self.apply_btn)
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.reject)
        bottom_layout.addWidget(self.cancel_btn)
        
        layout.addLayout(bottom_layout)
    
    def _apply_styles(self):
        """Apply dark theme styles."""
        self.setStyleSheet(f"""
            QDialog {{
                background-color: {COLORS['bg_primary']};
                color: {COLORS['text_primary']};
            }}
            QListWidget {{
                background-color: {COLORS['bg_secondary']};
                color: {COLORS['text_primary']};
                border: 1px solid {COLORS['border']};
                border-radius: 4px;
                padding: 4px;
            }}
            QListWidget::item {{
                padding: 6px;
                border-radius: 3px;
            }}
            QListWidget::item:selected {{
                background-color: {COLORS['accent']};
                color: white;
            }}
            QListWidget::item:hover {{
                background-color: {COLORS['bg_hover']};
            }}
            QPushButton {{
                background-color: {COLORS['bg_widget']};
                color: {COLORS['text_primary']};
                border: 1px solid {COLORS['border']};
                padding: 6px 16px;
                border-radius: 4px;
            }}
            QPushButton:hover {{
                background-color: {COLORS['bg_hover']};
                border-color: {COLORS['accent']};
            }}
            QCheckBox {{
                color: {COLORS['text_primary']};
            }}
            QCheckBox::indicator {{
                width: 16px;
                height: 16px;
            }}
        """)
    
    def _load_plugins(self):
        """Load plugins from plugin interface."""
        self.plugin_list.clear()
        self.plugin_items = {}
        
        if not self.plugin_interface or not self.plugin_interface.plugin_manager:
            # Add some example plugins for demonstration
            example_plugins = [
                ("Molecular Descriptors", "Calculate molecular descriptors", "1.0", "System"),
                ("Protein Analysis", "Analyze protein structures", "1.2", "System"),
                ("Ligand Tools", "Ligand preparation tools", "0.9", "System"),
                ("Visualization", "Advanced visualization options", "1.1", "System"),
            ]
            for name, desc, version, author in example_plugins:
                self._add_plugin_item(name, desc, version, author, False)
            return
        
        # Load actual plugins
        try:
            plugins = self.plugin_interface.plugin_manager.discover_plugins()
            for plugin_name, plugin_info in plugins.items():
                is_loaded = self.plugin_interface.plugin_manager.is_plugin_loaded(plugin_name)
                self._add_plugin_item(
                    plugin_name,
                    plugin_info.get('description', 'No description'),
                    plugin_info.get('version', 'Unknown'),
                    plugin_info.get('author', 'Unknown'),
                    is_loaded
                )
        except Exception as e:
            print(f"Error loading plugins: {e}")
    
    def _add_plugin_item(self, name, description, version, author, is_loaded):
        """Add a plugin item to the list."""
        item = QListWidgetItem(name)
        item.setData(Qt.ItemDataRole.UserRole, {
            'name': name,
            'description': description,
            'version': version,
            'author': author,
            'loaded': is_loaded
        })
        
        if is_loaded:
            item.setCheckState(Qt.CheckState.Checked)
        else:
            item.setCheckState(Qt.CheckState.Unchecked)
        
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
        self.plugin_list.addItem(item)
        self.plugin_items[name] = item
    
    def _on_plugin_selected(self, current, previous):
        """Handle plugin selection change."""
        if not current:
            self.name_label.setText("No plugin selected")
            self.description_text.clear()
            self.version_label.setText("Version: -")
            self.author_label.setText("Author: -")
            self.status_label.setText("Status: Not loaded")
            self.enable_check.setEnabled(False)
            return
        
        data = current.data(Qt.ItemDataRole.UserRole)
        if data:
            self.name_label.setText(data['name'])
            self.description_text.setText(data['description'])
            self.version_label.setText(f"Version: {data['version']}")
            self.author_label.setText(f"Author: {data['author']}")
            
            if data['loaded']:
                self.status_label.setText("Status: Loaded and active")
                self.status_label.setStyleSheet(f"color: {COLORS['success']}; font-size: 11px;")
            else:
                self.status_label.setText("Status: Not loaded")
                self.status_label.setStyleSheet(f"color: {COLORS['text_secondary']}; font-size: 11px;")
            
            self.enable_check.setEnabled(True)
            self.enable_check.setChecked(current.checkState() == Qt.CheckState.Checked)
    
    def _on_enable_changed(self, state):
        """Handle enable checkbox change."""
        current_item = self.plugin_list.currentItem()
        if current_item:
            if state == Qt.CheckState.Checked:
                current_item.setCheckState(Qt.CheckState.Checked)
            else:
                current_item.setCheckState(Qt.CheckState.Unchecked)
    
    def _select_all_plugins(self):
        """Select all plugins."""
        for i in range(self.plugin_list.count()):
            item = self.plugin_list.item(i)
            item.setCheckState(Qt.CheckState.Checked)
    
    def _select_no_plugins(self):
        """Deselect all plugins."""
        for i in range(self.plugin_list.count()):
            item = self.plugin_list.item(i)
            item.setCheckState(Qt.CheckState.Unchecked)
    
    def _load_selected_plugin(self):
        """Load the selected plugin."""
        current = self.plugin_list.currentItem()
        if not current:
            QMessageBox.warning(self, "No Selection", "Please select a plugin to load.")
            return
        
        data = current.data(Qt.ItemDataRole.UserRole)
        if self.plugin_interface and self.plugin_interface.plugin_manager:
            try:
                if self.plugin_interface.plugin_manager.load_plugin(data['name']):
                    QMessageBox.information(self, "Success", f"Plugin '{data['name']}' loaded successfully.")
                    self._load_plugins()  # Refresh list
                else:
                    QMessageBox.warning(self, "Error", f"Failed to load plugin '{data['name']}'.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error loading plugin: {str(e)}")
    
    def _unload_selected_plugin(self):
        """Unload the selected plugin."""
        current = self.plugin_list.currentItem()
        if not current:
            QMessageBox.warning(self, "No Selection", "Please select a plugin to unload.")
            return
        
        data = current.data(Qt.ItemDataRole.UserRole)
        if self.plugin_interface and self.plugin_interface.plugin_manager:
            try:
                if self.plugin_interface.plugin_manager.unload_plugin(data['name']):
                    QMessageBox.information(self, "Success", f"Plugin '{data['name']}' unloaded successfully.")
                    self._load_plugins()  # Refresh list
                else:
                    QMessageBox.warning(self, "Error", f"Failed to unload plugin '{data['name']}'.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error unloading plugin: {str(e)}")
    
    def _apply_changes(self):
        """Apply plugin selection changes."""
        active_plugins = []
        for i in range(self.plugin_list.count()):
            item = self.plugin_list.item(i)
            if item.checkState() == Qt.CheckState.Checked:
                data = item.data(Qt.ItemDataRole.UserRole)
                active_plugins.append(data['name'])
        
        self.active_plugins = active_plugins
        self.plugins_changed.emit(active_plugins)
        self.accept()
    
    def get_active_plugins(self):
        """Return the list of active plugins."""
        return self.active_plugins


def show_installed_plugins_dialog(plugin_interface=None, parent=None):
    """
    Show the installed plugins dialog.
    
    Args:
        plugin_interface: The plugin interface instance
        parent: Parent widget for the dialog
        
    Returns:
        tuple: (dialog_result, list of active plugins)
    """
    dialog = InstalledPluginsDialog(plugin_interface, parent)
    result = dialog.exec()
    
    if result == QDialog.DialogCode.Accepted:
        return True, dialog.get_active_plugins()
    else:
        return False, []
