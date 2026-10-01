"""
Plugin Installation and Template Utilities

This module provides install/uninstall/update logic and file management
for the plugin system, as well as plugin template generation.

Extracted from plugin_interface.py — pure refactor, no behavior change.
"""

import logging
from pathlib import Path

from src.shared.qt_compat import QMessageBox

logger = logging.getLogger("plugin.installer")


def install_plugin_from_file(plugin_manager, file_path: str, parent_widget=None):
    """
    Install a plugin from a file path.

    Copies the plugin file to the plugins directory, validates it,
    discovers it, and attempts to auto-load it.

    Args:
        plugin_manager: The PluginManager instance.
        file_path: Path to the .py plugin file.
        parent_widget: Optional parent widget for message boxes.

    Returns:
        tuple: (success: bool, messages: list[str])
    """
    messages = []
    try:
        messages.append(f"Installing plugin from: {file_path}")

        # Copy plugin to plugins directory
        plugins_dir = Path(plugin_manager.plugins_directory)
        plugin_file = Path(file_path)

        if not plugin_file.exists():
            raise FileNotFoundError(f"Plugin file not found: {file_path}")

        # Read plugin file
        with open(plugin_file, 'r', encoding='utf-8') as f:
            content = f.read()

        # Validate plugin
        from src.plugins.utils.validation import validate_plugin_file
        validation_result = validate_plugin_file(str(plugin_file))

        if not validation_result.get('valid', False):
            errors = validation_result.get('errors', [])
            raise ValueError(f"Plugin validation failed: {errors}")

        # Copy to plugins directory
        dest_path = plugins_dir / plugin_file.name
        with open(dest_path, 'w', encoding='utf-8') as f:
            f.write(content)

        messages.append(f"\u2713 Plugin copied to: {dest_path}")
        messages.append("\u2713 Refreshing plugin list...")

        # Refresh plugin list
        plugin_manager.discover_plugins()

        # Auto-load the installed plugin
        try:
            plugin_name = plugin_file.stem
            if plugin_manager.load_plugin(plugin_name):
                messages.append(f"\u2713 Plugin '{plugin_name}' loaded successfully!")
            else:
                messages.append(f"\u26a0 Plugin installed but could not be loaded automatically")
        except Exception as load_error:
            messages.append(f"\u26a0 Plugin installed but failed to load: {str(load_error)}")

        messages.append("\u2713 Plugin installed successfully!")
        return True, messages

    except Exception as e:
        messages.append(f"\u2717 Installation failed: {str(e)}")
        return False, messages


def save_template_to_file(plugin_type: str, file_path: str, parent_widget=None):
    """
    Save a plugin template to a file.

    Args:
        plugin_type: One of "analysis", "visualization", "io".
        file_path: Destination file path.
        parent_widget: Optional parent widget for message boxes.

    Returns:
        tuple: (success: bool, message: str)
    """
    templates = {
        "analysis": get_analysis_template(),
        "visualization": get_visualization_template(),
        "io": get_io_template()
    }

    template_content = templates.get(plugin_type)
    if not template_content:
        return False, f"Unknown template type: {plugin_type}"

    try:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(template_content)
        return True, f"\u2713 Template saved to: {file_path}"
    except Exception as e:
        return False, f"\u2717 Template creation failed: {str(e)}"


def get_analysis_template():
    """Get analysis plugin template."""
    return '''"""
Analysis Plugin Template

This is a template for creating analysis plugins that calculate
molecular properties and display results.
"""

from typing import Dict, Any
import logging

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QTextEdit, Qt
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType


class AnalysisWidget(PluginWidget):
    """Widget for the analysis plugin."""

    def __init__(self, plugin: 'AnalysisPlugin'):
        super().__init__(plugin)
        self.setup_ui()

    def setup_ui(self):
        """Setup the user interface."""
        self.widget = QWidget()
        layout = QVBoxLayout(self.widget)

        # Title
        title = QLabel("Analysis Plugin")
        title.setStyleSheet("font-size: 16px; font-weight: bold; margin: 10px;")
        layout.addWidget(title)

        # Calculate button
        self.calculate_btn = QPushButton("Calculate Properties")
        self.calculate_btn.clicked.connect(self.calculate_properties)
        layout.addWidget(self.calculate_btn)

        # Results table
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(2)
        self.results_table.setHorizontalHeaderLabels(["Property", "Value"])
        layout.addWidget(self.results_table)

        # Status text
        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(100)
        self.status_text.setReadOnly(True)
        layout.addWidget(self.status_text)

    def calculate_properties(self):
        """Calculate molecular properties."""
        if not self.plugin.current_molecule:
            self.status_text.append("No molecule loaded")
            return

        try:
            mol = self.plugin.current_molecule
            self.status_text.append("Calculating properties...")

            # Add your property calculations here
            # Example:
            # weight = calculate_molecular_weight(mol)
            # self.results_table.setItem(0, 0, QTableWidgetItem("Molecular Weight"))
            # self.results_table.setItem(0, 1, QTableWidgetItem(str(weight)))

            self.status_text.append("Calculation completed!")

        except Exception as e:
            self.status_text.append(f"Error: {str(e)}")

    def on_molecule_changed(self, molecule):
        """Handle molecule changes."""
        self.plugin.current_molecule = molecule
        if molecule:
            self.calculate_btn.setEnabled(True)
            self.status_text.append(f"Molecule loaded: {len(molecule.atoms)} atoms")
        else:
            self.calculate_btn.setEnabled(False)
            self.status_text.clear()


class AnalysisPlugin(BasePlugin):
    """Analysis plugin template."""

    def __init__(self):
        super().__init__()
        self.current_molecule = None

    def get_info(self) -> PluginInfo:
        """Get plugin information."""
        return PluginInfo(
            name="Analysis Plugin",
            version="1.0.0",
            description="Template for analysis plugins",
            author="Your Name",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        )

    def create_widget(self) -> 'AnalysisWidget':
        """Create the plugin widget."""
        return AnalysisWidget(self)

    def initialize(self):
        """Initialize the plugin."""
        self.logger.info("Analysis plugin initialized")
        return True

    def cleanup(self):
        """Clean up plugin resources."""
        self.logger.info("Analysis plugin cleaned up")

    def on_molecule_changed(self, molecule):
        """Handle molecule changes."""
        self.current_molecule = molecule
        if hasattr(self, 'widget') and self.widget:
            self.widget.on_molecule_changed(molecule)
'''


def get_visualization_template():
    """Get visualization plugin template."""
    return '''"""
Visualization Plugin Template

This is a template for creating visualization plugins that
display molecules in different ways.
"""

from typing import Dict, Any
import logging

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTextEdit, Qt, QPainter, QColor
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType


class VisualizationWidget(PluginWidget):
    """Widget for the visualization plugin."""

    def __init__(self, plugin: 'VisualizationPlugin'):
        super().__init__(plugin)
        self.setup_ui()

    def setup_ui(self):
        """Setup the user interface."""
        self.widget = QWidget()
        layout = QVBoxLayout(self.widget)

        # Title
        title = QLabel("Visualization Plugin")
        title.setStyleSheet("font-size: 16px; font-weight: bold; margin: 10px;")
        layout.addWidget(title)

        # Visualization area
        self.viz_widget = QWidget()
        self.viz_widget.setMinimumHeight(300)
        self.viz_widget.setStyleSheet("border: 1px solid #ccc; background-color: white;")
        layout.addWidget(self.viz_widget)

        # Controls
        controls_layout = QHBoxLayout()

        self.render_btn = QPushButton("Render")
        self.render_btn.clicked.connect(self.render_molecule)
        controls_layout.addWidget(self.render_btn)

        layout.addLayout(controls_layout)

        # Status text
        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(100)
        self.status_text.setReadOnly(True)
        layout.addWidget(self.status_text)

    def render_molecule(self):
        """Render the molecule."""
        if not self.plugin.current_molecule:
            self.status_text.append("No molecule loaded")
            return

        try:
            self.status_text.append("Rendering molecule...")

            # Add your visualization code here
            # Example: Custom drawing, charts, etc.

            self.status_text.append("Rendering completed!")

        except Exception as e:
            self.status_text.append(f"Error: {str(e)}")

    def on_molecule_changed(self, molecule):
        """Handle molecule changes."""
        self.plugin.current_molecule = molecule
        if molecule:
            self.render_btn.setEnabled(True)
            self.status_text.append(f"Molecule loaded: {len(molecule.atoms)} atoms")
        else:
            self.render_btn.setEnabled(False)
            self.status_text.clear()


class VisualizationPlugin(BasePlugin):
    """Visualization plugin template."""

    def __init__(self):
        super().__init__()
        self.current_molecule = None

    def get_info(self) -> PluginInfo:
        """Get plugin information."""
        return PluginInfo(
            name="Visualization Plugin",
            version="1.0.0",
            description="Template for visualization plugins",
            author="Your Name",
            plugin_type=PluginType.VISUALIZATION,
            dependencies=[]
        )

    def create_widget(self) -> 'VisualizationWidget':
        """Create the plugin widget."""
        return VisualizationWidget(self)

    def initialize(self):
        """Initialize the plugin."""
        self.logger.info("Visualization plugin initialized")
        return True

    def cleanup(self):
        """Clean up plugin resources."""
        self.logger.info("Visualization plugin cleaned up")

    def on_molecule_changed(self, molecule):
        """Handle molecule changes."""
        self.current_molecule = molecule
        if hasattr(self, 'widget') and self.widget:
            self.widget.on_molecule_changed(molecule)
'''


def get_io_template():
    """Get I/O plugin template."""
    return '''"""
I/O Plugin Template

This is a template for creating I/O plugins that handle
file import/export operations.
"""

from typing import Dict, Any
import logging

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTextEdit, Qt, QFileDialog
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType


class IOWidget(PluginWidget):
    """Widget for the I/O plugin."""

    def __init__(self, plugin: 'IOPlugin'):
        super().__init__(plugin)
        self.setup_ui()

    def setup_ui(self):
        """Setup the user interface."""
        self.widget = QWidget()
        layout = QVBoxLayout(self.widget)

        # Title
        title = QLabel("I/O Plugin")
        title.setStyleSheet("font-size: 16px; font-weight: bold; margin: 10px;")
        layout.addWidget(title)

        # File operations
        file_layout = QHBoxLayout()

        self.import_btn = QPushButton("Import File")
        self.import_btn.clicked.connect(self.import_file)
        file_layout.addWidget(self.import_btn)

        self.export_btn = QPushButton("Export File")
        self.export_btn.clicked.connect(self.export_file)
        file_layout.addWidget(self.export_btn)

        layout.addLayout(file_layout)

        # Status text
        self.status_text = QTextEdit()
        self.status_text.setReadOnly(True)
        layout.addWidget(self.status_text)

    def import_file(self):
        """Import a file."""
        try:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Import File", "", "All Files (*)"
            )

            if file_path:
                self.status_text.append(f"Importing: {file_path}")

                # Add your import logic here
                # Example: Read file, parse data, create molecule

                self.status_text.append("Import completed!")

        except Exception as e:
            self.status_text.append(f"Import error: {str(e)}")

    def export_file(self):
        """Export a file."""
        if not self.plugin.current_molecule:
            self.status_text.append("No molecule to export")
            return

        try:
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Export File", "", "All Files (*)"
            )

            if file_path:
                self.status_text.append(f"Exporting: {file_path}")

                # Add your export logic here
                # Example: Write molecule data to file

                self.status_text.append("Export completed!")

        except Exception as e:
            self.status_text.append(f"Export error: {str(e)}")

    def on_molecule_changed(self, molecule):
        """Handle molecule changes."""
        self.plugin.current_molecule = molecule
        if molecule:
            self.status_text.append(f"Molecule loaded: {len(molecule.atoms)} atoms")
        else:
            self.status_text.clear()


class IOPlugin(BasePlugin):
    """I/O plugin template."""

    def __init__(self):
        super().__init__()
        self.current_molecule = None

    def get_info(self) -> PluginInfo:
        """Get plugin information."""
        return PluginInfo(
            name="I/O Plugin",
            version="1.0.0",
            description="Template for I/O plugins",
            author="Your Name",
            plugin_type=PluginType.IO,
            dependencies=[]
        )

    def create_widget(self) -> 'IOWidget':
        """Create the plugin widget."""
        return IOWidget(self)

    def initialize(self):
        """Initialize the plugin."""
        self.logger.info("I/O plugin initialized")
        return True

    def cleanup(self):
        """Clean up plugin resources."""
        self.logger.info("I/O plugin cleaned up")

    def on_molecule_changed(self, molecule):
        """Handle molecule changes."""
        self.current_molecule = molecule
        if hasattr(self, 'widget') and self.widget:
            self.widget.on_molecule_changed(molecule)
'''
