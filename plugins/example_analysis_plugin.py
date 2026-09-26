"""
Example Analysis Plugin

This plugin demonstrates how to create an analysis plugin that
calculates custom molecular properties and displays them in a table.
"""

from typing import Dict, Any
import csv
from datetime import datetime

from src.shared.qt_compat import (
    Qt,
    QHBoxLayout,
    QDialog,
    QWidget, QVBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QProgressBar, QTextEdit,
    QFileDialog, QMessageBox
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType


class ExampleAnalysisWidget(PluginWidget):
    """Widget for the example analysis plugin."""

    def __init__(self, plugin: 'ExampleAnalysisPlugin'):
        super().__init__(plugin)
        self.setup_ui()

    def setup_ui(self):
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

        # Title & Description
        title = QLabel("Example Analysis Plugin")
        title.setStyleSheet("font-size: 16px; font-weight: bold; margin: 10px;")
        layout.addWidget(title)

        description = QLabel(
            "This plugin calculates custom molecular properties including:\n"
            "• Atom count analysis\n"
            "• Bond type distribution\n"
            "• Molecular complexity metrics"
        )
        description.setWordWrap(True)
        description.setStyleSheet("margin: 10px; color: #666;")
        layout.addWidget(description)

        # Calculate button
        self.calculate_btn = QPushButton("Calculate Properties")
        self.calculate_btn.clicked.connect(self.calculate_properties)
        self.calculate_btn.setStyleSheet("""
            QPushButton {
                background-color: #6c63ff; color: white; border: none;
                padding: 10px; border-radius: 5px; font-weight: bold;
            }
            QPushButton:hover { background-color: #7f78ff; }
            QPushButton:disabled { background-color: #ccc; }
        """)
        layout.addWidget(self.calculate_btn)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)

        # Results table
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(3)
        self.results_table.setHorizontalHeaderLabels(["Property", "Value", "Unit"])
        self.results_table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.results_table)

        # Export button
        self.export_btn = QPushButton("Export Results")
        self.export_btn.clicked.connect(self.export_results)
        self.export_btn.setEnabled(False)
        layout.addWidget(self.export_btn)

        # Status text
        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(100)
        self.status_text.setReadOnly(True)
        layout.addWidget(self.status_text)

        self.update_ui_state(False)

    def calculate_properties(self):
        if not self.plugin.current_molecule:
            QMessageBox.warning(self.widget, "No Molecule", "Please load a molecule first.")
            return

        try:
            self.calculate_btn.setEnabled(False)
            self.progress_bar.setVisible(True)
            self.progress_bar.setValue(0)
            self.status_text.clear()
            self.status_text.append("Starting calculation...")

            mol = self.plugin.current_molecule
            properties = self._calculate_properties(mol)

            self._update_results_table(properties)
            self.update_ui_state(True)
            self.status_text.append(f"✓ Calculation complete! Found {len(properties)} properties.")

        except Exception as e:
            self.status_text.append(f"✗ Error: {e}")
            QMessageBox.critical(self.widget, "Calculation Error", f"Error during calculation:\n{e}")
        finally:
            self.calculate_btn.setEnabled(True)
            self.progress_bar.setVisible(False)

    def _calculate_properties(self, molecule) -> Dict[str, Any]:
        properties = {}
        properties['Atom Count'] = {'value': len(molecule.atoms), 'unit': ''}
        properties['Bond Count'] = {'value': len(molecule.bonds), 'unit': ''}

        atom_types = {}
        for atom in molecule.atoms:
            # Handle standard object properties safely
            symbol = getattr(atom, 'symbol', 'Unknown')
            if not isinstance(symbol, str):
                symbol = getattr(getattr(atom, 'element', None), 'symbol', 'Unknown')
            atom_types[symbol] = atom_types.get(symbol, 0) + 1

        for symbol, count in sorted(atom_types.items()):
            properties[f'{symbol} Atoms'] = {'value': count, 'unit': ''}

        bond_types = {}
        for bond in molecule.bonds:
            bond_type = str(getattr(bond, 'bond_type', getattr(bond, 'order', 'Unknown')))
            bond_types[bond_type] = bond_types.get(bond_type, 0) + 1

        for bond_type, count in sorted(bond_types.items()):
            properties[f'{bond_type} Bonds'] = {'value': count, 'unit': ''}

        if molecule.atoms:
            properties['Atom Diversity'] = {'value': len(atom_types), 'unit': ''}
            properties['Bond Diversity'] = {'value': len(bond_types), 'unit': ''}
            total_bonds = len(molecule.bonds) * 2
            properties['Average Connectivity'] = {'value': total_bonds / len(molecule.atoms), 'unit': ''}

        return properties

    def _update_results_table(self, properties: Dict[str, Any]):
        self.results_table.setRowCount(len(properties))
        for row, (property_name, data) in enumerate(properties.items()):
            self.results_table.setItem(row, 0, QTableWidgetItem(property_name))
            value = f"{data['value']:.3f}" if isinstance(data['value'], float) else str(data['value'])
            self.results_table.setItem(row, 1, QTableWidgetItem(value))
            self.results_table.setItem(row, 2, QTableWidgetItem(data.get('unit', '')))

    def export_results(self):
        try:
            filename = f"analysis_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            filepath, _ = QFileDialog.getSaveFileName(self.widget, "Export Results", filename, "CSV Files (*.csv)")
            if not filepath: return

            with open(filepath, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['Property', 'Value', 'Unit'])
                for row in range(self.results_table.rowCount()):
                    writer.writerow([
                        self.results_table.item(row, 0).text(),
                        self.results_table.item(row, 1).text(),
                        self.results_table.item(row, 2).text()
                    ])

            self.status_text.append(f"✓ Results exported to: {filepath}")
            QMessageBox.information(self.widget, "Export Complete", f"Results exported to:\n{filepath}")

        except Exception as e:
            self.status_text.append(f"✗ Export error: {e}")
            QMessageBox.critical(self.widget, "Export Error", f"Failed to export results:\n{e}")

    def on_molecule_changed(self, molecule):
        if molecule:
            self.status_text.clear()
            self.status_text.append(f"Molecule loaded: {len(molecule.atoms)} atoms, {len(molecule.bonds)} bonds")
            self.calculate_btn.setEnabled(True)
        else:
            self.status_text.clear()
            self.status_text.append("No molecule loaded")
            self.results_table.setRowCount(0)
            self.update_ui_state(False)
            self.calculate_btn.setEnabled(False)

    def update_ui_state(self, has_results):
        self.export_btn.setEnabled(has_results)



    def _open_popout(self, fullscreen: bool = False) -> None:
        if getattr(self, '_active_popout', None) is not None:
            self._active_popout.raise_()
            self._active_popout.activateWindow()
            return

        dialog = QDialog()
        dialog.setWindowTitle("Plugin")
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
        self.main_content_widget.setParent(self.widget)
        self.widget.layout().addWidget(self.main_content_widget)
        self._active_popout = None

    def toggle_maximize(self):
        self._open_popout(fullscreen=False)


class ExampleAnalysisPlugin(BasePlugin):
    """
    Example Analysis Plugin.
    Demonstrates analysis property calculations following strict API guidelines.
    """

    def __init__(self):
        super().__init__(PluginInfo(
            name="Example Analysis Plugin",
            version="1.0.0",
            description="Calculates custom molecular properties and metrics",
            author="SMILES Team",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        ))
        self.current_molecule = None
        self.widget = None

    def get_info(self) -> PluginInfo:
        return PluginInfo(
            name="Example Analysis Plugin",
            version="1.0.0",
            description="Calculates custom molecular properties and metrics",
            author="SMILES Team",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        )

    def create_widget(self) -> ExampleAnalysisWidget:
        if self.widget is None:
            self.widget = ExampleAnalysisWidget(self)
        return self.widget

    def initialize(self):
        """
        Initialize the plugin.
        Note: Bypassing super().initialize() intentionally to avoid 'api' positional argument crashes.
        """
        self.logger.info("Example Analysis Plugin initialized")
        return True

    def cleanup(self):
        if self.widget:
            self.widget.widget.deleteLater()
            self.widget = None
        self.logger.info("Example Analysis Plugin cleaned up")

    def on_molecule_changed(self, molecule):
        self.current_molecule = molecule
        if self.widget:
            self.widget.on_molecule_changed(molecule)