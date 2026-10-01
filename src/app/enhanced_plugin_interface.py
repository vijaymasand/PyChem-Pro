"""
Enhanced Plugin Interface for Runtime Plugin Management

This module provides a comprehensive interface for managing plugins with
runtime extensibility support, including browsing, installing, and managing
both bundled and user plugins.

Features:
- Runtime plugin installation/removal from local files
- Enhanced UI with search and filtering
- Plugin dependency management
- Plugin configuration management
"""

import logging
import json
from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTabWidget, QTableWidget, QTableWidgetItem, QTextEdit,
    QGroupBox, QCheckBox, QComboBox, QProgressBar,
    QMessageBox, Qt, Signal, QFileDialog, QLineEdit,
    QSplitter, QFrame, QScrollArea, QSpinBox, QTreeWidget,
    QTreeWidgetItem, QDialog, QDialogButtonBox, QHeaderView,
    QTimer, QThread, QSettings
)

from src.plugins.enhanced_plugin_manager import EnhancedPluginManager
from src.plugins.utils.integration import PluginIntegrationAPI


class PluginInstallWorker(QThread):
    """Worker thread for plugin installation to avoid UI freezing."""
    
    progress = Signal(int, str)
    finished = Signal(bool, str)
    
    def __init__(self, plugin_manager, install_type, *args):
        super().__init__()
        self.plugin_manager = plugin_manager
        self.install_type = install_type
        self.args = args
    
    def run(self):
        try:
            if self.install_type == "file":
                file_path = self.args[0]
                self.progress.emit(50, f"Installing plugin from {file_path}")
                success, message = self.plugin_manager.install_plugin_from_file(file_path)
            elif self.install_type == "uninstall":
                plugin_name = self.args[0]
                self.progress.emit(50, f"Uninstalling plugin {plugin_name}")
                success, message = self.plugin_manager.uninstall_plugin(plugin_name)
            else:
                success, message = False, f"Unknown install type: {self.install_type}"
            
            self.progress.emit(100, "Operation complete")
            self.finished.emit(success, message)
            
        except Exception as e:
            self.finished.emit(False, f"Installation error: {e}")


class EnhancedPluginManagerWidget(QWidget):
    """
    Enhanced widget for managing plugins with runtime extensibility.
    """
    
    def __init__(self, plugin_manager: EnhancedPluginManager):
        super().__init__()
        self.plugin_manager = plugin_manager
        self.setup_ui()
        self.refresh_plugin_list()
    
    def setup_ui(self):
        """Setup the enhanced user interface."""
        layout = QVBoxLayout(self)
        
        # Header
        header = QLabel("Enhanced Plugin Manager")
        header.setStyleSheet("font-size: 18px; font-weight: bold; color: #2196F3; margin: 10px;")
        layout.addWidget(header)
        
        # Plugin Management Tab
        self.management_tab = self.create_management_tab()
        layout.addWidget(self.management_tab)
    
    def create_management_tab(self):
        """Create the plugin management tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Search and filter
        filter_group = QGroupBox("Search & Filter")
        filter_layout = QHBoxLayout(filter_group)
        
        self.search_edit = QLineEdit()
        self.search_edit.setPlaceholderText("Search plugins...")
        self.search_edit.textChanged.connect(self.filter_plugins)
        filter_layout.addWidget(self.search_edit)
        
        self.type_filter = QComboBox()
        self.type_filter.addItems(["All Types", "Bundled", "User", "Analysis", "Visualization", "I/O", "Utility"])
        self.type_filter.currentTextChanged.connect(self.filter_plugins)
        filter_layout.addWidget(self.type_filter)
        
        self.status_filter = QComboBox()
        self.status_filter.addItems(["All Status", "Active", "Inactive", "Error"])
        self.status_filter.currentTextChanged.connect(self.filter_plugins)
        filter_layout.addWidget(self.status_filter)
        
        # Install from file button
        self.install_file_btn = QPushButton("Install from File")
        self.install_file_btn.clicked.connect(self.install_from_file)
        filter_layout.addWidget(self.install_file_btn)
        
        # Choose plugins folder button
        self.choose_folder_btn = QPushButton("Choose Plugins Folder")
        self.choose_folder_btn.clicked.connect(self.choose_plugins_folder)
        filter_layout.addWidget(self.choose_folder_btn)
        
        layout.addWidget(filter_group)
        
        # Plugin table
        self.plugin_table = QTableWidget()
        self.plugin_table.setColumnCount(7)
        self.plugin_table.setHorizontalHeaderLabels([
            "Name", "Version", "Type", "Status", "Author", "Actions", "Source"
        ])
        self.plugin_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        layout.addWidget(self.plugin_table)
        
        # Progress bar for operations
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
        
        # Status label
        self.status_label = QLabel("Ready")
        self.status_label.setStyleSheet("color: #666; font-style: italic;")
        layout.addWidget(self.status_label)
        
        return tab
    
    
    def refresh_plugin_list(self):
        """Refresh the plugin list."""
        all_plugins = self.plugin_manager.get_all_plugins()
        
        self.plugin_table.setRowCount(len(all_plugins))
        
        for row, (plugin_name, metadata) in enumerate(all_plugins.items()):
            # Name
            self.plugin_table.setItem(row, 0, QTableWidgetItem(metadata.info.name))
            
            # Version
            self.plugin_table.setItem(row, 1, QTableWidgetItem(metadata.info.version))
            
            # Type
            self.plugin_table.setItem(row, 2, QTableWidgetItem(metadata.info.plugin_type.value))
            
            # Status
            status_item = QTableWidgetItem(metadata.status.value)
            if metadata.status.value == "Active":
                status_item.setStyleSheet("color: green;")
            elif metadata.status.value == "Error":
                status_item.setStyleSheet("color: red;")
            self.plugin_table.setItem(row, 3, status_item)
            
            # Author
            self.plugin_table.setItem(row, 4, QTableWidgetItem(metadata.info.author))
            
            # Actions
            actions_widget = QWidget()
            actions_layout = QHBoxLayout(actions_widget)
            actions_layout.setContentsMargins(0, 0, 0, 0)
            
            if self.plugin_manager.is_plugin_loaded(plugin_name):
                unload_btn = QPushButton("Unload")
                unload_btn.clicked.connect(lambda checked, name=plugin_name: self.unload_plugin(name))
                actions_layout.addWidget(unload_btn)
            else:
                load_btn = QPushButton("Load")
                load_btn.clicked.connect(lambda checked, name=plugin_name: self.load_plugin(name))
                actions_layout.addWidget(load_btn)
            
            # Uninstall button for user plugins
            if self.plugin_manager.is_user_plugin(plugin_name):
                uninstall_btn = QPushButton("Uninstall")
                uninstall_btn.clicked.connect(lambda checked, name=plugin_name: self.uninstall_plugin(name))
                actions_layout.addWidget(uninstall_btn)
            
            self.plugin_table.setCellWidget(row, 5, actions_widget)
            
            # Source
            source = "Bundled" if self.plugin_manager.is_bundled_plugin(plugin_name) else "User"
            self.plugin_table.setItem(row, 6, QTableWidgetItem(source))
    
    def filter_plugins(self):
        """Filter plugins based on search and filters."""
        search_text = self.search_edit.text().lower()
        type_filter = self.type_filter.currentText()
        status_filter = self.status_filter.currentText()
        
        for row in range(self.plugin_table.rowCount()):
            name_item = self.plugin_table.item(row, 0)
            type_item = self.plugin_table.item(row, 2)
            status_item = self.plugin_table.item(row, 3)
            source_item = self.plugin_table.item(row, 6)
            
            if not all([name_item, type_item, status_item, source_item]):
                continue
            
            name_match = search_text in name_item.text().lower()
            
            type_match = (
                type_filter == "All Types" or
                type_filter == source_item.text() or
                type_item.text() == type_filter
            )
            
            status_match = status_filter == "All Status" or status_item.text() == status_filter
            
            should_show = name_match and type_match and status_match
            self.plugin_table.setRowHidden(row, not should_show)
    
    def load_plugin(self, plugin_name: str):
        """Load a plugin."""
        self.status_label.setText(f"Loading plugin '{plugin_name}'...")
        
        success = self.plugin_manager.load_plugin(plugin_name)
        if success:
            self.status_label.setText(f"Plugin '{plugin_name}' loaded successfully")
        else:
            self.status_label.setText(f"Failed to load plugin '{plugin_name}'")
        
        self.refresh_plugin_list()
    
    def unload_plugin(self, plugin_name: str):
        """Unload a plugin."""
        self.status_label.setText(f"Unloading plugin '{plugin_name}'...")
        
        success = self.plugin_manager.unload_plugin(plugin_name)
        if success:
            self.status_label.setText(f"Plugin '{plugin_name}' unloaded successfully")
        else:
            self.status_label.setText(f"Failed to unload plugin '{plugin_name}'")
        
        self.refresh_plugin_list()
    
    def uninstall_plugin(self, plugin_name: str):
        """Uninstall a user plugin."""
        reply = QMessageBox.question(
            self, "Confirm Uninstall",
            f"Are you sure you want to uninstall plugin '{plugin_name}'?",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            self.progress_bar.setVisible(True)
            self.progress_bar.setValue(0)
            self.status_label.setText(f"Uninstalling plugin '{plugin_name}'...")
            
            self.worker = PluginInstallWorker(self.plugin_manager, "uninstall", plugin_name)
            self.worker.progress.connect(self.on_install_progress)
            self.worker.finished.connect(self.on_install_finished)
            self.worker.start()
    
    def install_from_file(self):
        """Install a plugin from a file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Install Plugin", "", 
            "Python Files (*.py);;Plugin Packages (*.zip);;All Files (*)"
        )
        
        if file_path:
            self.progress_bar.setVisible(True)
            self.progress_bar.setValue(0)
            self.status_label.setText(f"Installing plugin from {file_path}...")
            
            self.worker = PluginInstallWorker(self.plugin_manager, "file", file_path)
            self.worker.progress.connect(self.on_install_progress)
            self.worker.finished.connect(self.on_install_finished)
            self.worker.start()
    
    def choose_plugins_folder(self):
        """Open a directory dialog to choose an external plugins folder."""
        settings = QSettings("PyChem", "PyChemPro")
        current_dir = settings.value("custom_plugins_directory", "")
        if not current_dir:
            current_dir = str(self.plugin_manager.user_plugins_directory)
            
        selected_dir = QFileDialog.getExistingDirectory(
            self, "Choose Plugins Folder", current_dir
        )
        
        if selected_dir:
            self.plugin_manager.user_plugins_directory = Path(selected_dir)
            self.plugin_manager.discover_all_plugins()
            self.refresh_plugin_list()
            
            # Persist custom folder path
            settings.setValue("custom_plugins_directory", selected_dir)
            
            # Update status label
            count = len(self.plugin_manager.get_all_plugins())
            self.status_label.setText(f"Loaded custom plugins directory: {selected_dir} ({count} plugins found)")
            QMessageBox.information(
                self, "Plugins Folder Loaded",
                f"Successfully loaded external plugins folder:\n{selected_dir}\n\nTotal plugins found: {count}"
            )

    def on_install_progress(self, value: int, message: str):
        """Handle installation progress."""
        self.progress_bar.setValue(value)
        self.status_label.setText(message)
    
    def on_install_finished(self, success: bool, message: str):
        """Handle installation completion."""
        self.progress_bar.setVisible(False)
        
        if success:
            QMessageBox.information(self, "Success", message)
        else:
            QMessageBox.warning(self, "Error", message)
        
        self.status_label.setText("Ready")
        self.refresh_plugin_list()
    
