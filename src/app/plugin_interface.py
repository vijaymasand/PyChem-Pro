"""
Enhanced Plugin Interface for Main Application

This module provides an improved interface between the main application
and the plugin system, with better UI, browsing capabilities, and
comprehensive plugin management.

The PluginBrowserWidget has been extracted to plugin_card.py and
install/template utilities to plugin_installer.py.
"""

import logging
import json
from typing import Dict, List, Optional, Any
from pathlib import Path

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTabWidget, QTableWidget, QTableWidgetItem, QTextEdit,
    QGroupBox, QCheckBox, QComboBox, QProgressBar,
    QMessageBox, Qt, Signal, QFileDialog, QLineEdit,
    QSplitter, QFrame, QScrollArea, QSpinBox
)
from src.plugins import PluginManager, PLUGIN_API_VERSION
from src.plugins.utils.integration import PluginIntegrationAPI

# Import from extracted modules
from src.app.plugin_card import PluginBrowserWidget
from src.app.plugin_installer import (
    install_plugin_from_file,
    save_template_to_file,
    get_analysis_template,
    get_visualization_template,
    get_io_template,
)


class EnhancedPluginManagerWidget(QWidget):
    """
    Enhanced widget for managing plugins.

    Provides comprehensive plugin management with browsing,
    installation, and detailed information.
    """

    def __init__(self, plugin_manager: PluginManager):
        super().__init__()
        self.plugin_manager = plugin_manager
        self.setup_ui()
        self.refresh_plugin_list()

    def setup_ui(self):
        """Setup the enhanced user interface."""
        layout = QVBoxLayout(self)

        # Header
        header = QLabel("Plugin Manager")
        header.setStyleSheet("font-size: 18px; font-weight: bold; color: #2196F3; margin: 10px;")
        layout.addWidget(header)

        # Create tab widget
        self.tab_widget = QTabWidget()
        layout.addWidget(self.tab_widget)

        # Plugin Management Tab
        self.management_tab = self.create_management_tab()
        self.tab_widget.addTab(self.management_tab, "Manage Plugins")

        # Plugin Browser Tab
        self.browser_tab = PluginBrowserWidget(self.plugin_manager)
        self.browser_tab.plugin_installed.connect(self.on_plugin_installed)
        self.tab_widget.addTab(self.browser_tab, "Browse & Install")

        # Plugin Documentation Tab
        self.docs_tab = self.create_documentation_tab()
        self.tab_widget.addTab(self.docs_tab, "Documentation")

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
        self.type_filter.addItems(["All Types", "Analysis", "Visualization", "I/O", "Utility"])
        self.type_filter.currentTextChanged.connect(self.filter_plugins)
        filter_layout.addWidget(self.type_filter)

        self.status_filter = QComboBox()
        self.status_filter.addItems(["All Status", "Active", "Inactive", "Error"])
        self.status_filter.currentTextChanged.connect(self.filter_plugins)
        filter_layout.addWidget(self.status_filter)

        layout.addWidget(filter_group)

        # Plugin table
        self.plugin_table = QTableWidget()
        self.plugin_table.setColumnCount(7)
        self.plugin_table.setHorizontalHeaderLabels([
            "Name", "Version", "Type", "Status", "Author", "Description", "Actions"
        ])
        self.plugin_table.horizontalHeader().setStretchLastSection(True)
        self.plugin_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.plugin_table.itemSelectionChanged.connect(self.on_plugin_selected)
        layout.addWidget(self.plugin_table)

        # Control buttons
        button_layout = QHBoxLayout()

        self.refresh_btn = QPushButton("Refresh")
        self.refresh_btn.clicked.connect(self.refresh_plugin_list)
        self.refresh_btn.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """)
        button_layout.addWidget(self.refresh_btn)

        self.load_all_btn = QPushButton("Load All")
        self.load_all_btn.clicked.connect(self.load_all_plugins)
        button_layout.addWidget(self.load_all_btn)

        self.unload_all_btn = QPushButton("Unload All")
        self.unload_all_btn.clicked.connect(self.unload_all_plugins)
        button_layout.addWidget(self.unload_all_btn)

        self.reload_btn = QPushButton("Reload Selected")
        self.reload_btn.clicked.connect(self.reload_selected_plugin)
        button_layout.addWidget(self.reload_btn)

        layout.addLayout(button_layout)

        # Plugin details
        details_group = QGroupBox("Plugin Details")
        details_layout = QVBoxLayout(details_group)

        self.details_text = QTextEdit()
        self.details_text.setReadOnly(True)
        self.details_text.setMaximumHeight(150)
        details_layout.addWidget(self.details_text)

        layout.addWidget(details_group)

        return tab

    def create_documentation_tab(self):
        """Create the documentation tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)

        # Documentation content
        docs_text = QTextEdit()
        docs_text.setReadOnly(True)
        docs_text.setHtml(self.get_plugin_documentation())
        layout.addWidget(docs_text)

        return tab

    def get_plugin_documentation(self):
        """Get plugin documentation HTML."""
        return """
        <h2>Plugin System Documentation</h2>

        <h3>Overview</h3>
        <p>The SMILES application supports a flexible plugin system that allows you to extend functionality with custom analysis, visualization, and I/O plugins.</p>

        <h3>Plugin Types</h3>
        <ul>
        <li><strong>Analysis Plugins:</strong> Calculate molecular properties, perform statistical analysis, generate reports</li>
        <li><strong>Visualization Plugins:</strong> Create custom molecular visualizations, charts, and interactive displays</li>
        <li><strong>I/O Plugins:</strong> Handle file import/export for different formats, connect to external databases</li>
        <li><strong>Utility Plugins:</strong> Provide helper functions, tools, and utilities</li>
        </ul>

        <h3>Creating Plugins</h3>
        <p>Use the "Browse & Install" tab to create plugins from templates. Templates include:</p>
        <ul>
        <li>Basic plugin structure</li>
        <li>UI components</li>
        <li>Error handling</li>
        <li>Integration with the main application</li>
        </ul>

        <h3>Plugin Structure</h3>
        <p>Each plugin must include:</p>
        <ul>
        <li>A class inheriting from <code>BasePlugin</code></li>
        <li>A widget class inheriting from <code>PluginWidget</code></li>
        <li>Plugin information via <code>PluginInfo</code></li>
        <li>Proper initialization and cleanup methods</li>
        </ul>

        <h3>Installation</h3>
        <p>To install plugins:</p>
        <ol>
        <li>Go to the "Browse & Install" tab</li>
        <li>Click "Browse..." to select a plugin file</li>
        <li>Click "Install Plugin"</li>
        <li>The plugin will be automatically detected and available</li>
        </ol>

        <h3>API Reference</h3>
        <p>Plugins have access to:</p>
        <ul>
        <li>Current molecule data</li>
        <li>Qt widgets for UI creation</li>
        <li>Logging system</li>
        <li>Plugin integration API</li>
        <li>Event system for molecule changes</li>
        </ul>

        <h3>Best Practices</h3>
        <ul>
        <li>Follow the plugin template structure</li>
        <li>Include proper error handling</li>
        <li>Use the logging system for debugging</li>
        <li>Test plugins thoroughly before distribution</li>
        <li>Document plugin functionality clearly</li>
        </ul>

        <h3>Troubleshooting</h3>
        <p>Common issues:</p>
        <ul>
        <li><strong>Plugin not loading:</strong> Check for syntax errors and missing dependencies</li>
        <li><strong>UI not appearing:</strong> Ensure the widget is properly created and returned</li>
        <li><strong>Errors on molecule change:</strong> Implement proper null checks in <code>on_molecule_changed</code></li>
        </ul>
        """

    def refresh_plugin_list(self):
        """Refresh the plugin list with enhanced information."""
        try:
            # Discover plugins
            plugins = self.plugin_manager.discover_plugins()

            # Store all plugins for filtering
            self.all_plugins = plugins

            # Update table
            self.update_plugin_table(plugins)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error refreshing plugin list: {e}")

    def update_plugin_table(self, plugins):
        """Update the plugin table with filtered plugins."""
        self.plugin_table.setRowCount(len(plugins))

        for row, (name, metadata) in enumerate(plugins.items()):
            # Name
            name_item = QTableWidgetItem(name)
            name_item.setToolTip(f"Plugin: {name}")
            self.plugin_table.setItem(row, 0, name_item)

            # Version
            version_item = QTableWidgetItem(metadata.info.version)
            self.plugin_table.setItem(row, 1, version_item)

            # Type
            type_item = QTableWidgetItem(metadata.info.plugin_type.value)
            self.plugin_table.setItem(row, 2, type_item)

            # Status
            status = metadata.status.value if metadata.status else "inactive"
            status_item = QTableWidgetItem(status)
            if metadata.status and metadata.status.value == "active":
                status_item.setStyleSheet("color: green; font-weight: bold;")
            elif metadata.status and metadata.status.value == "error":
                status_item.setStyleSheet("color: red; font-weight: bold;")
            self.plugin_table.setItem(row, 3, status_item)

            # Author
            author_item = QTableWidgetItem(metadata.info.author)
            self.plugin_table.setItem(row, 4, author_item)

            # Description
            desc_item = QTableWidgetItem(metadata.info.description[:50] + "..." if len(metadata.info.description) > 50 else metadata.info.description)
            desc_item.setToolTip(metadata.info.description)
            self.plugin_table.setItem(row, 5, desc_item)

            # Actions
            actions_widget = QWidget()
            actions_layout = QHBoxLayout(actions_widget)
            actions_layout.setContentsMargins(0, 0, 0, 0)

            if self.plugin_manager.is_plugin_loaded(name):
                unload_btn = QPushButton("Unload")
                unload_btn.clicked.connect(lambda checked, n=name: self.unload_plugin(n))
                unload_btn.setStyleSheet("background-color: #f44336; color: white; padding: 4px 8px;")
                actions_layout.addWidget(unload_btn)
            else:
                load_btn = QPushButton("Load")
                load_btn.clicked.connect(lambda checked, n=name: self.load_plugin(n))
                load_btn.setStyleSheet("background-color: #4CAF50; color: white; padding: 4px 8px;")
                actions_layout.addWidget(load_btn)

            self.plugin_table.setCellWidget(row, 6, actions_widget)

        self.plugin_table.resizeColumnsToContents()

    def filter_plugins(self):
        """Filter plugins based on search and filters."""
        if not hasattr(self, 'all_plugins'):
            return

        search_text = self.search_edit.text().lower()
        type_filter = self.type_filter.currentText()
        status_filter = self.status_filter.currentText()

        filtered = {}

        for name, metadata in self.all_plugins.items():
            # Search filter
            if search_text and search_text not in name.lower() and search_text not in metadata.info.description.lower():
                continue

            # Type filter
            if type_filter != "All Types" and metadata.info.plugin_type.value != type_filter.lower():
                continue

            # Status filter
            status = metadata.status.value if metadata.status else "inactive"
            if status_filter != "All Status" and status != status_filter.lower():
                continue

            filtered[name] = metadata

        self.update_plugin_table(filtered)

    def on_plugin_selected(self):
        """Handle plugin selection."""
        current_row = self.plugin_table.currentRow()
        if current_row >= 0:
            name_item = self.plugin_table.item(current_row, 0)
            if name_item:
                plugin_name = name_item.text()
                self.show_plugin_details(plugin_name)

    def show_plugin_details(self, plugin_name):
        """Show detailed information about a plugin."""
        try:
            plugins = self.plugin_manager.discover_plugins()
            if plugin_name in plugins:
                metadata = plugins[plugin_name]

                details = f"""
<strong>Plugin Details:</strong>

Name: {metadata.info.name}
Version: {metadata.info.version}
Type: {metadata.info.plugin_type.value}
Author: {metadata.info.author}
Status: {metadata.status.value if metadata.status else 'inactive'}

Description:
{metadata.info.description}

Dependencies: {', '.join(metadata.info.dependencies) if metadata.info.dependencies else 'None'}

File: {metadata.module_path}
Min API Version: {metadata.info.min_api_version}
"""

                self.details_text.setHtml(details)

        except Exception as e:
            self.details_text.setPlainText(f"Error loading details: {str(e)}")

    def load_plugin(self, plugin_name: str):
        """Load a specific plugin."""
        try:
            if self.plugin_manager.load_plugin(plugin_name):
                QMessageBox.information(self, "Success", f"Plugin {plugin_name} loaded successfully")
                self.refresh_plugin_list()
            else:
                QMessageBox.critical(self, "Error", f"Failed to load plugin {plugin_name}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading plugin: {e}")

    def unload_plugin(self, plugin_name: str):
        """Unload a specific plugin."""
        try:
            if self.plugin_manager.unload_plugin(plugin_name):
                QMessageBox.information(self, "Success", f"Plugin {plugin_name} unloaded successfully")
                self.refresh_plugin_list()
            else:
                QMessageBox.critical(self, "Error", f"Failed to unload plugin {plugin_name}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error unloading plugin: {e}")

    def load_all_plugins(self):
        """Load all available plugins."""
        try:
            plugins = self.plugin_manager.discover_plugins()
            loaded_count = 0

            for plugin_name in plugins.keys():
                if not self.plugin_manager.is_plugin_loaded(plugin_name):
                    if self.plugin_manager.load_plugin(plugin_name):
                        loaded_count += 1

            QMessageBox.information(self, "Success", f"Loaded {loaded_count} plugins")
            self.refresh_plugin_list()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading plugins: {e}")

    def unload_all_plugins(self):
        """Unload all loaded plugins."""
        try:
            plugins = self.plugin_manager.discover_plugins()
            unloaded_count = 0

            for plugin_name in plugins.keys():
                if self.plugin_manager.is_plugin_loaded(plugin_name):
                    if self.plugin_manager.unload_plugin(plugin_name):
                        unloaded_count += 1

            QMessageBox.information(self, "Success", f"Unloaded {unloaded_count} plugins")
            self.refresh_plugin_list()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error unloading plugins: {e}")

    def reload_selected_plugin(self):
        """Reload the selected plugin."""
        current_row = self.plugin_table.currentRow()
        if current_row >= 0:
            name_item = self.plugin_table.item(current_row, 0)
            if name_item:
                plugin_name = name_item.text()

                # Unload first
                if self.plugin_manager.is_plugin_loaded(plugin_name):
                    self.plugin_manager.unload_plugin(plugin_name)

                # Then load
                if self.plugin_manager.load_plugin(plugin_name):
                    QMessageBox.information(self, "Success", f"Plugin {plugin_name} reloaded successfully")
                    self.refresh_plugin_list()
                else:
                    QMessageBox.critical(self, "Error", f"Failed to reload plugin {plugin_name}")

    def on_plugin_installed(self, plugin_name):
        """Handle plugin installation."""
        self.refresh_plugin_list()


class PluginInterface:
    """
    Enhanced interface between the main application and plugin system.

    Provides comprehensive plugin management with browsing,
    installation, and integration capabilities.
    """

    def __init__(self, parent=None):
        super().__init__()
        self.parent = parent
        self.plugin_manager = None
        self.plugin_widget = None
        self.plugin_tab_widget = None  # Add this missing attribute
        self.plugin_tabs = {}  # Add this missing attribute
        self.logger = logging.getLogger("plugin.interface")

    def initialize_plugin_system(self):
        """Initialize the plugin system."""
        try:
            from src.plugins import PluginManager

            # Get plugins directory
            plugins_dir = Path("plugins")
            if not plugins_dir.exists():
                plugins_dir.mkdir(exist_ok=True)

            # Create plugin manager
            self.plugin_manager = PluginManager(str(plugins_dir))

            # Create and set plugin API
            from src.plugins.utils.integration import PluginIntegrationAPI
            self.plugin_api = PluginIntegrationAPI(self.parent)
            self.plugin_manager.set_api(self.plugin_api)

            # Discover plugins
            self.plugin_manager.discover_plugins()

            self.logger.info("Plugin system initialized successfully with API")
            return True

        except Exception as e:
            self.logger.error(f"Failed to initialize plugin system: {e}")
            return False

    def setup_event_handlers(self):
        """Setup event handlers for plugin system."""
        if self.plugin_manager:
            # Handle molecule changes
            self.plugin_manager.add_event_handler('molecule_changed', self.on_molecule_changed)

    def create_plugin_manager_widget(self) -> EnhancedPluginManagerWidget:
        """Create the enhanced plugin manager widget."""
        if self.plugin_manager:
            return EnhancedPluginManagerWidget(self.plugin_manager)
        else:
            return QWidget()

    def create_plugin_tabs(self, tab_widget: QTabWidget):
        """Create tabs for all loaded plugins."""
        if not self.plugin_manager:
            return

        try:
            plugins = self.plugin_manager.discover_plugins()
            loaded_count = 0

            for plugin_name, metadata in plugins.items():
                if self.plugin_manager.is_plugin_loaded(plugin_name):
                    try:
                        plugin = self.plugin_manager.get_plugin(plugin_name)
                        if plugin and hasattr(plugin, 'create_widget'):
                            widget = plugin.create_widget()
                            if widget:
                                # Handle both PluginWidget wrappers and direct QWidgets
                                real_widget = widget.widget if hasattr(widget, 'widget') and widget.widget else widget
                                tab_widget.addTab(real_widget, plugin_name)
                                loaded_count += 1

                    except Exception as e:
                        self.logger.error(f"Error creating tab for plugin {plugin_name}: {e}")

            self.logger.info(f"Created {loaded_count} plugin tabs")

        except Exception as e:
            self.logger.error(f"Error creating plugin tabs: {e}")

    def on_molecule_changed(self, event):
        """
        Handle molecule changes from plugin system.

        Args:
            event: Molecule changed event
        """
        molecule = event.data.get('molecule')

        # Notify all plugins
        if self.plugin_manager:
            self.plugin_manager.on_molecule_changed(molecule)

    def get_plugin_manager(self) -> PluginManager:
        """Get the plugin manager."""
        return self.plugin_manager

    def get_plugin_api(self) -> PluginIntegrationAPI:
        """Get the plugin API."""
        if hasattr(self, 'plugin_api'):
            return self.plugin_api
        if self.plugin_manager:
            return PluginIntegrationAPI(self.plugin_manager)
        return None

    def show_plugin_manager(self):
        """Show the plugin manager dialog."""
        if not self.plugin_manager:
            QMessageBox.warning(self.parent, "Plugin System",
                              "Plugin system not initialized")
            return

        # Create and show plugin manager widget
        manager_widget = self.create_plugin_manager_widget()
        manager_widget.setWindowTitle("Plugin Manager")
        manager_widget.resize(800, 600)
        manager_widget.show()

        # Store reference if needed
        self.plugin_widget = manager_widget

    def refresh_plugin_tabs(self):
        """Refresh plugin tabs."""
        if not self.plugin_tab_widget:
            return

        # Clear existing plugin tabs (except manager tab)
        for i in range(self.plugin_tab_widget.count() - 1, 0, -1):
            tab_text = self.plugin_tab_widget.tabText(i)
            if tab_text != "Plugin Manager":
                self.plugin_tab_widget.removeTab(i)

        # Add tabs for loaded plugins
        loaded_plugins = self.plugin_manager.get_loaded_plugins()
        for plugin_name, plugin in loaded_plugins.items():
            try:
                widget = plugin.create_widget()
                if widget:
                    self.plugin_tab_widget.addTab(widget, plugin_name)
                    self.plugin_tabs[plugin_name] = widget
                    self.logger.info(f"Added tab for plugin: {plugin_name}")
            except Exception as e:
                self.logger.error(f"Error creating tab for plugin {plugin_name}: {e}")

    def show_plugin_tab(self, plugin_name: str):
        """
        Show a specific plugin tab.

        Args:
            plugin_name: Name of the plugin to show
        """
        if not self.plugin_tab_widget:
            return

        # Find the tab index for the plugin
        for i in range(self.plugin_tab_widget.count()):
            if self.plugin_tab_widget.tabText(i) == plugin_name:
                self.plugin_tab_widget.setCurrentIndex(i)
                self.plugin_manager.on_plugin_activated(plugin_name)
                return

        # If tab not found, try to create it
        plugin = self.plugin_manager.get_plugin(plugin_name)
        if plugin:
            try:
                widget = plugin.create_widget()
                if widget:
                    index = self.plugin_tab_widget.addTab(widget, plugin_name)
                    self.plugin_tab_widget.setCurrentIndex(index)
                    self.plugin_manager.on_plugin_activated(plugin_name)
            except Exception as e:
                self.logger.error(f"Error creating tab for plugin {plugin_name}: {e}")

    def cleanup(self):
        """Clean up the plugin interface."""
        try:
            # Unload all plugins
            loaded_plugins = list(self.plugin_manager.get_loaded_plugins().keys())
            for plugin_name in loaded_plugins:
                self.plugin_manager.unload_plugin(plugin_name)

            self.logger.info("Plugin interface cleaned up")

        except Exception as e:
            self.logger.error(f"Error cleaning up plugin interface: {e}")

    def get_plugin_statistics(self) -> Dict[str, Any]:
        """Get plugin system statistics."""
        stats = self.plugin_manager.get_statistics()
        stats.update({
            'api_version': PLUGIN_API_VERSION,
            'loaded_plugins': len(self.plugin_manager.get_loaded_plugins()),
            'available_plugins': len(self.plugin_manager.get_all_plugins()),
            'active_tabs': len(self.plugin_tabs)
        })
        return stats
