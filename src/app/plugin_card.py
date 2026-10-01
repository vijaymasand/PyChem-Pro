"""
Plugin Browser Widget

This module provides the PluginBrowserWidget for browsing,
installing, and creating plugin templates.

Extracted from plugin_interface.py — pure refactor, no behavior change.
"""

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QTextEdit, QLineEdit, QFileDialog, QMessageBox,
    Signal
)
from src.plugins import PluginManager
from src.app.plugin_installer import (
    install_plugin_from_file,
    save_template_to_file,
)


class PluginBrowserWidget(QWidget):
    """
    Widget for browsing and installing new plugins.
    """

    plugin_installed = Signal(str)

    def __init__(self, plugin_manager: PluginManager):
        super().__init__()
        self.plugin_manager = plugin_manager
        self.setup_ui()

    def setup_ui(self):
        """Setup the plugin browser interface."""
        layout = QVBoxLayout(self)

        # Title and description
        title = QLabel("Plugin Browser")
        title.setStyleSheet("font-size: 16px; font-weight: bold; color: #2196F3;")
        layout.addWidget(title)

        desc = QLabel("Browse and install new plugins from local files or repositories")
        desc.setStyleSheet("color: #666; margin-bottom: 10px;")
        layout.addWidget(desc)

        # Installation section
        install_group = QGroupBox("Install Plugin")
        install_layout = QVBoxLayout(install_group)

        # File selection
        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select plugin file (.py)...")
        file_layout.addWidget(self.file_path_edit)

        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self.browse_plugin_file)
        browse_btn.setStyleSheet("""
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
        file_layout.addWidget(browse_btn)
        install_layout.addLayout(file_layout)

        # Install button
        self.install_btn = QPushButton("Install Plugin")
        self.install_btn.clicked.connect(self.install_plugin)
        self.install_btn.setEnabled(False)
        self.install_btn.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                border: none;
                padding: 10px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #1976D2;
            }
            QPushButton:disabled {
                background-color: #ccc;
            }
        """)
        install_layout.addWidget(self.install_btn)

        layout.addWidget(install_group)

        # Plugin templates
        templates_group = QGroupBox("Plugin Templates")
        templates_layout = QVBoxLayout(templates_group)

        template_desc = QLabel("Create new plugins from templates:")
        templates_layout.addWidget(template_desc)

        template_buttons = QHBoxLayout()

        analysis_btn = QPushButton("Analysis Plugin")
        analysis_btn.clicked.connect(lambda: self.create_template("analysis"))
        template_buttons.addWidget(analysis_btn)

        viz_btn = QPushButton("Visualization Plugin")
        viz_btn.clicked.connect(lambda: self.create_template("visualization"))
        template_buttons.addWidget(viz_btn)

        io_btn = QPushButton("I/O Plugin")
        io_btn.clicked.connect(lambda: self.create_template("io"))
        template_buttons.addWidget(io_btn)

        templates_layout.addLayout(template_buttons)
        layout.addWidget(templates_group)

        # Status
        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(100)
        self.status_text.setReadOnly(True)
        layout.addWidget(self.status_text)

        # Connect file path edit
        self.file_path_edit.textChanged.connect(self.on_file_path_changed)

    def browse_plugin_file(self):
        """Browse for plugin file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Plugin File", "", "Python Files (*.py);;All Files (*)"
        )
        if file_path:
            self.file_path_edit.setText(file_path)

    def on_file_path_changed(self):
        """Handle file path change."""
        self.install_btn.setEnabled(bool(self.file_path_edit.text().strip()))

    def install_plugin(self):
        """Install plugin from file."""
        file_path = self.file_path_edit.text().strip()
        if not file_path:
            return

        self.status_text.clear()

        success, messages = install_plugin_from_file(
            self.plugin_manager, file_path, parent_widget=self
        )

        for msg in messages:
            self.status_text.append(msg)

        if success:
            from pathlib import Path
            plugin_name = Path(file_path).stem
            self.plugin_installed.emit(plugin_name)
            self.file_path_edit.clear()
            QMessageBox.information(self, "Success",
                                   f"Plugin '{plugin_name}' installed successfully!")
        else:
            QMessageBox.critical(self, "Installation Error",
                               f"Failed to install plugin:\n{messages[-1] if messages else 'Unknown error'}")

    def create_template(self, plugin_type):
        """Create a new plugin from template."""
        try:
            # Save template file
            file_path, _ = QFileDialog.getSaveFileName(
                self, f"Save {plugin_type.title()} Plugin Template",
                f"{plugin_type}_plugin.py", "Python Files (*.py)"
            )

            if file_path:
                success, message = save_template_to_file(plugin_type, file_path, parent_widget=self)
                self.status_text.append(message)

                if success:
                    QMessageBox.information(self, "Template Created",
                                           f"Plugin template saved to:\n{file_path}")
                else:
                    QMessageBox.critical(self, "Template Error",
                                       f"Failed to create template:\n{message}")

        except Exception as e:
            self.status_text.append(f"\u2717 Template creation failed: {str(e)}")
            QMessageBox.critical(self, "Template Error",
                               f"Failed to create template:\n{str(e)}")
