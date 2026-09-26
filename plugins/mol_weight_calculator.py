"""
Molecular Weight Calculator Plugin

This plugin demonstrates how to create a simple analysis plugin that
calculates and displays molecular weight information.

Features:
- Molecular weight calculation
- Element composition analysis
- Export results to CSV
- Real-time updates when molecule changes
"""

from typing import Dict, Any, List
import logging
import csv
from io import StringIO

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QTextEdit, QProgressBar,
    Qt, Signal, QFileDialog, QMessageBox
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType


class MolWeightWidget(PluginWidget):
    """
    Widget for the molecular weight calculator plugin.
    """

    def __init__(self, plugin: 'MolWeightCalculatorPlugin'):
        """
        Initialize the molecular weight widget.

        Args:
            plugin: The plugin instance
        """
        super().__init__(plugin)
        self.setup_ui()

    def setup_ui(self):
        """Setup the user interface."""
        self.widget = QWidget()
        _base_layout = QVBoxLayout(self.widget)
        
        _top_bar = QHBoxLayout()
        _top_bar.addStretch()
        self.btn_maximize = QPushButton("Maximize")
        self.btn_maximize.clicked.connect(self.toggle_maximize)
        _top_bar.addWidget(self.btn_maximize)
        self.btn_fullscreen = QPushButton("Full Screen")
        self.btn_fullscreen.clicked.connect(lambda: self._open_popout(fullscreen=True))
        _top_bar.addWidget(self.btn_fullscreen)
        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(self.widget.close)
        _top_bar.addWidget(self.btn_close)
        _base_layout.addLayout(_top_bar)
        
        self.main_content_widget = QWidget()
        _base_layout.addWidget(self.main_content_widget)
        
        layout = QVBoxLayout(self.main_content_widget)

        # Title
        title = QLabel("Molecular Weight Calculator")
        title.setStyleSheet("font-size: 16px; font-weight: bold; margin: 10px;")
        layout.addWidget(title)

        # Description
        desc = QLabel("Calculate molecular weight and element composition")
        desc.setStyleSheet("color: #666; margin-bottom: 10px;")
        layout.addWidget(desc)

        # Calculate button
        self.calculate_btn = QPushButton("Calculate Molecular Weight")
        self.calculate_btn.clicked.connect(self.calculate_properties)
        self.calculate_btn.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """)
        layout.addWidget(self.calculate_btn)

        # Results table
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(2)
        self.results_table.setHorizontalHeaderLabels(["Property", "Value"])
        self.results_table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.results_table)

        # Element composition table
        self.element_table = QTableWidget()
        self.element_table.setColumnCount(3)
        self.element_table.setHorizontalHeaderLabels(["Element", "Count", "Mass Contribution"])
        self.element_table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.element_table)

        # Export button
        self.export_btn = QPushButton("Export Results")
        self.export_btn.clicked.connect(self.export_results)
        self.export_btn.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #1976D2;
            }
        """)
        layout.addWidget(self.export_btn)

        # Status text
        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(100)
        self.status_text.setReadOnly(True)
        layout.addWidget(self.status_text)

        # Initially disable buttons
        self.update_ui_state(False)

    def calculate_properties(self):
        """Calculate molecular weight and composition."""
        if not self.plugin.current_molecule:
            QMessageBox.warning(self.widget, "No Molecule",
                              "Please load a molecule first.")
            return

        try:
            self.status_text.clear()
            self.status_text.append("Calculating molecular properties...")

            mol = self.plugin.current_molecule

            # Calculate molecular weight
            total_weight = 0.0
            element_counts = {}
            element_masses = {}

            # Atomic weights (simplified periodic table)
            atomic_weights = {
                'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
                'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'Br': 79.904,
                'I': 126.904, 'F': 18.998, 'Na': 22.990, 'K': 39.098,
                'Mg': 24.305, 'Ca': 40.078, 'Fe': 55.845, 'Zn': 65.38,
                'Cu': 63.546, 'Mn': 54.938, 'Co': 58.933, 'Ni': 58.693
            }

            for atom in mol.atoms:
                element = atom.element.symbol
                count = element_counts.get(element, 0) + 1
                element_counts[element] = count

                weight = atomic_weights.get(element, 0.0)
                mass_contrib = element_masses.get(element, 0.0) + weight
                element_masses[element] = mass_contrib

                total_weight += weight

            # Update main results table
            self.results_table.setRowCount(4)
            self.results_table.setItem(0, 0, QTableWidgetItem("Molecular Weight"))
            self.results_table.setItem(0, 1, QTableWidgetItem(f"{total_weight:.3f} g/mol"))

            self.results_table.setItem(1, 0, QTableWidgetItem("Number of Atoms"))
            self.results_table.setItem(1, 1, QTableWidgetItem(str(len(mol.atoms))))

            self.results_table.setItem(2, 0, QTableWidgetItem("Number of Bonds"))
            self.results_table.setItem(2, 1, QTableWidgetItem(str(len(mol.bonds))))

            self.results_table.setItem(3, 0, QTableWidgetItem("Molecule Name"))
            self.results_table.setItem(3, 1, QTableWidgetItem(mol.name or "Unnamed"))

            # Update element composition table
            self.element_table.setRowCount(len(element_counts))
            for i, (element, count) in enumerate(element_counts.items()):
                self.element_table.setItem(i, 0, QTableWidgetItem(element))
                self.element_table.setItem(i, 1, QTableWidgetItem(str(count)))
                self.element_table.setItem(i, 2, QTableWidgetItem(f"{element_masses[element]:.3f}"))

            self.status_text.append(f"✓ Molecular weight: {total_weight:.3f} g/mol")
            self.status_text.append(f"✓ {len(element_counts)} different elements")
            self.status_text.append("✓ Calculation completed successfully!")

            self.update_ui_state(True)

        except Exception as e:
            self.status_text.append(f"✗ Error: {str(e)}")
            QMessageBox.critical(self.widget, "Calculation Error",
                               f"Failed to calculate properties:\n{str(e)}")

    def export_results(self):
        """Export results to CSV file."""
        try:
            filename, _ = QFileDialog.getSaveFileName(
                self.widget, "Export Results", "", "CSV Files (*.csv)"
            )

            if not filename:
                return

            output = StringIO()
            writer = csv.writer(output)

            # Write molecular weight results
            writer.writerow(["Molecular Weight Results"])
            writer.writerow(["Property", "Value"])
            for row in range(self.results_table.rowCount()):
                prop = self.results_table.item(row, 0).text()
                value = self.results_table.item(row, 1).text()
                writer.writerow([prop, value])

            writer.writerow([])  # Empty row

            # Write element composition
            writer.writerow(["Element Composition"])
            writer.writerow(["Element", "Count", "Mass Contribution"])
            for row in range(self.element_table.rowCount()):
                element = self.element_table.item(row, 0).text()
                count = self.element_table.item(row, 1).text()
                mass = self.element_table.item(row, 2).text()
                writer.writerow([element, count, mass])

            # Save to file
            with open(filename, 'w', newline='') as f:
                f.write(output.getvalue())

            self.status_text.append(f"✓ Results exported to: {filename}")
            QMessageBox.information(self.widget, "Export Successful",
                                   f"Results saved to:\n{filename}")

        except Exception as e:
            self.status_text.append(f"✗ Export error: {str(e)}")
            QMessageBox.critical(self.widget, "Export Error",
                               f"Failed to export results:\n{str(e)}")

    def update_ui_state(self, has_results):
        """Update UI button states."""
        self.export_btn.setEnabled(has_results)

    def on_molecule_changed(self, molecule):
        """Handle molecule changes."""
        self.plugin.current_molecule = molecule
        if molecule:
            self.calculate_btn.setEnabled(True)
            self.status_text.clear()
            self.status_text.append(f"Molecule loaded: {molecule.name or 'Unnamed'}")
            self.status_text.append(f"  - {len(molecule.atoms)} atoms")
            self.status_text.append(f"  - {len(molecule.bonds)} bonds")
        else:
            self.calculate_btn.setEnabled(False)
            self.status_text.clear()
            self.status_text.append("No molecule loaded")
            # Clear tables
            self.results_table.setRowCount(0)
            self.element_table.setRowCount(0)
            self.update_ui_state(False)


    def _open_popout(self, fullscreen: bool = False) -> None:
        if getattr(self, '_active_popout', None) is not None:
            self._active_popout.raise_()
            self._active_popout.activateWindow()
            return
        
        try:
            from src.shared.qt_compat import QDialog, QVBoxLayout, Qt
        except:
            from PySide6.QtWidgets import QDialog, QVBoxLayout
            from PySide6.QtCore import Qt

        dialog = QDialog()
        dialog.setWindowTitle('Plugin')
        dialog.setWindowFlags(Qt.Window)
        vbox = QVBoxLayout(dialog)
        vbox.setContentsMargins(4, 4, 4, 4)

        self.main_content_widget.setParent(dialog)
        vbox.addWidget(self.main_content_widget)

        dialog.finished.connect(self._restore_from_popout)
        self._active_popout = dialog

        if fullscreen:
            dialog.showFullScreen()
        else:
            dialog.showMaximized()

    def _restore_from_popout(self) -> None:
        w = getattr(self, 'widget', self)
        self.main_content_widget.setParent(w)
        w.layout().addWidget(self.main_content_widget)
        self._active_popout = None

    def toggle_maximize(self):
        self._open_popout(fullscreen=False)

class MolWeightCalculatorPlugin(BasePlugin):
    """
    Molecular Weight Calculator Plugin.

    A simple plugin that calculates molecular weight and element composition.
    """

    def __init__(self):
        """Initialize the plugin."""
        super().__init__(PluginInfo(
            name="Molecular Weight Calculator",
            version="1.0.0",
            description="Calculate molecular weight and element composition",
            author="SMILES Team",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        ))
        self.current_molecule = None
        self.widget = None

    def get_info(self) -> PluginInfo:
        """Get plugin information."""
        return PluginInfo(
            name="Molecular Weight Calculator",
            version="1.0.0",
            description="Calculate molecular weight and element composition",
            author="SMILES Team",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        )

    def create_widget(self) -> 'MolWeightWidget':
        """Create the plugin widget."""
        if self.widget is None:
            self.widget = MolWeightWidget(self)
        return self.widget

    def initialize(self):
        """Initialize the plugin."""
        self.logger.info("Molecular Weight Calculator plugin initialized")
        return True

    def cleanup(self):
        """Clean up plugin resources."""
        if self.widget:
            self.widget.widget.deleteLater()
            self.widget = None
        self.logger.info("Molecular Weight Calculator plugin cleaned up")

    def on_molecule_changed(self, molecule):
        """Handle molecule changes."""
        self.current_molecule = molecule
        if self.widget:
            self.widget.on_molecule_changed(molecule)
