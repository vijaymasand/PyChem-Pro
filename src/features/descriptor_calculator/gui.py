"""
Molecular Descriptor Calculator GUI

Comprehensive GUI for molecular descriptor calculation with selection support,
progress tracking, and export capabilities.
"""

import sys
import csv
import os
from typing import Dict, List, Any, Optional
from datetime import datetime

from src.shared.qt_compat import *

# Only define GUI classes if Qt framework is available
if QT_FRAMEWORK is not None:

    from ..cheminformatics.services.atom_properties import AtomPropertyAnalyzer
    from .descriptor_engine import DescriptorEngine
    from .pydes.engine import PyDesEngine
    from .descriptor_types import (
        DescriptorCategory, DescriptorInfo, DescriptorResult,
        CalculationProgress, AtomSelection, SelectionType
    )
    from ...shared.ui.theme import COLORS

    class DescriptorCalculationThread(QThread):
        """Thread for descriptor calculations to avoid GUI freezing."""
        
        progress_updated = Signal(object)
        calculation_finished = Signal(dict)
        error_occurred = Signal(str)
        
        def __init__(self, molecules, n_jobs=-1):
            super().__init__()
            self.molecules = molecules if isinstance(molecules, list) else [molecules]
            self.n_jobs = n_jobs
        
        def run(self):
            try:
                # Use the new PyDesEngine for batch calculation
                results_list = PyDesEngine.calculate_batch(self.molecules, n_jobs=self.n_jobs)
                # For UI display, we'll just emit the first result if there's only one, or a summary
                self.calculation_finished.emit({'batch_results': results_list})
                
            except Exception as e:
                self.error_occurred.emit(str(e))

    class SelectionBuilder(QWidget):
        """Widget for building atom selections."""
        
        def __init__(self, molecule, parent=None):
            super().__init__(parent)
            self.molecule = molecule
            self.selection_history = []
            self.init_ui()
        
        def init_ui(self):
            layout = QVBoxLayout(self)
            
            # Selection type group
            type_group = QGroupBox("Selection Type")
            type_layout = QVBoxLayout(type_group)
            
            self.selection_type_group = QButtonGroup()
            
            self.all_radio = QRadioButton("All Atoms")
            self.all_radio.setChecked(True)
            self.selection_type_group.addButton(self.all_radio, 0)
            
            self.custom_radio = QRadioButton("Custom Selection")
            self.selection_type_group.addButton(self.custom_radio, 1)
            
            type_layout.addWidget(self.all_radio)
            type_layout.addWidget(self.custom_radio)
            
            layout.addWidget(type_group)
            
            # Custom selection input
            custom_group = QGroupBox("Custom Selection")
            custom_layout = QVBoxLayout(custom_group)
            
            self.custom_input = QLineEdit()
            self.custom_input.setPlaceholderText("Enter atom indices (e.g., 1,2,3-5,7)")
            custom_layout.addWidget(QLabel("Atom Indices:"))
            custom_layout.addWidget(self.custom_input)
            
            layout.addWidget(custom_group)
            
            # Current selection display
            current_group = QGroupBox("Current Selection")
            current_layout = QVBoxLayout(current_group)
            
            self.current_selection = QTextEdit()
            self.current_selection.setMaximumHeight(100)
            self.current_selection.setReadOnly(True)
            current_layout.addWidget(QLabel("Selected Atoms:"))
            current_layout.addWidget(self.current_selection)
            
            layout.addWidget(current_group)
            
            # Selection buttons
            button_layout = QHBoxLayout()
            
            self.preview_btn = QPushButton("Preview Selection")
            self.preview_btn.clicked.connect(self.preview_selection)
            button_layout.addWidget(self.preview_btn)
            
            self.clear_btn = QPushButton("Clear Selection")
            self.clear_btn.clicked.connect(self.clear_selection)
            button_layout.addWidget(self.clear_btn)
            
            layout.addLayout(button_layout)
            
            # Update selection display
            self.update_selection_display()
        
        def get_selection(self) -> AtomSelection:
            """Get current atom selection."""
            selected_type = self.selection_type_group.checkedId()
            
            if selected_type == 0:  # All atoms
                return AtomSelection(
                    SelectionType.ALL,
                    list(range(len(self.molecule.atoms))),
                    "All atoms"
                )
            elif selected_type == 1:  # Custom
                indices = self.parse_custom_selection()
                return AtomSelection(
                    SelectionType.CUSTOM,
                    indices,
                    f"Custom: {self.custom_input.text()}"
                )
            
            return AtomSelection(SelectionType.ALL, [], "Empty selection")
        
        def parse_custom_selection(self) -> List[int]:
            """Parse custom selection input."""
            text = self.custom_input.text().strip()
            if not text:
                return []
            
            indices = []
            parts = text.split(',')
            
            for part in parts:
                part = part.strip()
                if '-' in part:
                    # Range like 3-7
                    start, end = map(int, part.split('-'))
                    indices.extend(range(start, end + 1))
                else:
                    # Single index
                    indices.append(int(part))
            
            return sorted(set(indices))
        
        def preview_selection(self):
            """Preview current selection."""
            selection = self.get_selection()
            self.update_selection_display()
        
        def clear_selection(self):
            """Clear current selection."""
            self.custom_input.clear()
            self.all_radio.setChecked(True)
            self.update_selection_display()
        
        def update_selection_display(self):
            """Update the selection display."""
            selection = self.get_selection()
            
            if len(selection.atom_indices) <= 20:
                indices_str = ", ".join(map(str, selection.atom_indices))
            else:
                indices_str = f"{len(selection.atom_indices)} atoms: {selection.atom_indices[:10]}..."
            
            display_text = f"Type: {selection.selection_type.value}\n"
            display_text += f"Description: {selection.description}\n"
            display_text += f"Indices: {indices_str}"
            
            self.current_selection.setText(display_text)

    class DescriptorConfigWidget(QWidget):
        """Widget for configuring descriptor categories."""
        
        def __init__(self, parent=None):
            super().__init__(parent)
            self.engine = DescriptorEngine()
            self.init_ui()
        
        def init_ui(self):
            layout = QVBoxLayout(self)
            
            # Category selection
            category_group = QGroupBox("Descriptor Categories")
            category_layout = QVBoxLayout(category_group)
            
            self.category_checkboxes = {}
            for category in DescriptorCategory:
                checkbox = QCheckBox(f"{category.value} ({len(self.engine.descriptors[category])} descriptors)")
                checkbox.setChecked(True)
                self.category_checkboxes[category] = checkbox
                category_layout.addWidget(checkbox)
            
            layout.addWidget(category_group)
            
            # Quick select buttons
            button_layout = QHBoxLayout()
            
            self.select_all_btn = QPushButton("Select All")
            self.select_all_btn.clicked.connect(self.select_all_categories)
            button_layout.addWidget(self.select_all_btn)
            
            self.select_none_btn = QPushButton("Select None")
            self.select_none_btn.clicked.connect(self.select_none_categories)
            button_layout.addWidget(self.select_none_btn)
            
            layout.addLayout(button_layout)
            
            # Update preview
            self.update_descriptor_preview()
            
            # Connect checkbox signals
            for checkbox in self.category_checkboxes.values():
                checkbox.toggled.connect(self.update_descriptor_preview)
        
        def get_selected_categories(self) -> List[DescriptorCategory]:
            """Get selected descriptor categories."""
            selected = []
            for category, checkbox in self.category_checkboxes.items():
                if checkbox.isChecked():
                    selected.append(category)
            return selected
        
        def select_all_categories(self):
            """Select all categories."""
            for checkbox in self.category_checkboxes.values():
                checkbox.setChecked(True)
        
        def select_none_categories(self):
            """Select no categories."""
            for checkbox in self.category_checkboxes.values():
                checkbox.setChecked(False)
        
        def update_descriptor_preview(self):
            """Update descriptor preview."""
            selected_categories = self.get_selected_categories()
            
            preview_text = f"Selected {len(selected_categories)} categories:\n"
            total_descriptors = 0
            
            for category in selected_categories:
                descriptors = self.engine.descriptors[category]
                total_descriptors += len(descriptors)
                preview_text += f"  • {category.value}: {len(descriptors)} descriptors\n"
            
            preview_text += f"\nTotal: {total_descriptors} descriptors"
            
            # Note: This would need a QTextEdit widget to display
            # For now, just print to console
            print(preview_text)

    class DescriptorCalculatorDialog(QMainWindow):
        """Main descriptor calculator dialog."""
        
        def __init__(self, molecule=None, parent=None):
            super().__init__(parent)
            self.molecules = [molecule] if molecule else []
            self.current_results = {}
            self.batch_results = []
            self.calculation_thread = None
            
            self.setWindowTitle("Molecular Descriptor Calculator (PyDes)")
            self.setGeometry(100, 100, 1200, 800)
            
            self.init_ui()
            self.apply_styles()
        
        def init_ui(self):
            """Initialize the user interface."""
            central_widget = QWidget()
            self.setCentralWidget(central_widget)
            
            # Main layout
            main_layout = QHBoxLayout(central_widget)
            
            # Left panel (tabs)
            left_panel = QWidget()
            left_layout = QVBoxLayout(left_panel)
            
            self.tab_widget = QTabWidget()
            
            # Selection tab (If single molecule)
            if self.molecules and len(self.molecules) == 1:
                self.selection_widget = SelectionBuilder(self.molecules[0])
                self.tab_widget.addTab(self.selection_widget, "Selection")
            else:
                # Placeholder for multi-molecule
                lbl = QLabel(f"Loaded {len(self.molecules)} molecules for batch processing.")
                lbl.setAlignment(Qt.AlignCenter)
                self.tab_widget.addTab(lbl, "Selection")
            
            # Configuration tab
            self.config_widget = DescriptorConfigWidget()
            self.tab_widget.addTab(self.config_widget, "Configuration")
            
            left_layout.addWidget(self.tab_widget)
            
            # Right panel (controls and results)
            right_panel = QWidget()
            right_layout = QVBoxLayout(right_panel)
            
            # Results display
            results_group = QGroupBox("Results")
            results_layout = QVBoxLayout(results_group)
            
            self.results_table = QTableWidget()
            self.results_table.setColumnCount(3)
            self.results_table.setHorizontalHeaderLabels(["Descriptor", "Value", "Category"])
            self.results_table.horizontalHeader().setStretchLastSection(True)
            results_layout.addWidget(self.results_table)
            
            right_layout.addWidget(results_group)
            
            # Export controls
            export_group = QGroupBox("Export")
            export_layout = QVBoxLayout(export_group)
            
            # Export buttons
            export_buttons_layout = QHBoxLayout()
            
            self.export_results_btn = QPushButton("Export Results")
            self.export_results_btn.clicked.connect(self.export_results)
            self.export_results_btn.setEnabled(False)
            export_buttons_layout.addWidget(self.export_results_btn)
            
            self.export_documentation_btn = QPushButton("Export Documentation")
            self.export_documentation_btn.clicked.connect(self.export_documentation)
            self.export_documentation_btn.setEnabled(False)
            export_buttons_layout.addWidget(self.export_documentation_btn)
            
            export_layout.addLayout(export_buttons_layout)
            right_layout.addWidget(export_group)
            
            # Calculation controls
            controls_group = QGroupBox("Calculation Controls")
            controls_layout = QVBoxLayout(controls_group)
            
            # Progress bar
            self.progress_bar = QProgressBar()
            self.progress_bar.setVisible(False)
            controls_layout.addWidget(self.progress_bar)
            
            self.progress_label = QLabel("Ready")
            controls_layout.addWidget(self.progress_label)
            
            # Control buttons
            button_layout = QHBoxLayout()
            
            self.load_files_btn = QPushButton("Load Files...")
            self.load_files_btn.clicked.connect(self.load_files)
            button_layout.addWidget(self.load_files_btn)
            
            self.calculate_btn = QPushButton("Calculate PyDes Descriptors")
            self.calculate_btn.clicked.connect(self.calculate_descriptors)
            button_layout.addWidget(self.calculate_btn)
            
            self.stop_btn = QPushButton("Stop")
            self.stop_btn.clicked.connect(self.stop_calculation)
            self.stop_btn.setEnabled(False)
            button_layout.addWidget(self.stop_btn)
            
            controls_layout.addLayout(button_layout)
            right_layout.addWidget(controls_group)
            
            # Add panels to main layout
            splitter = QSplitter(Qt.Horizontal)
            splitter.addWidget(left_panel)
            splitter.addWidget(right_panel)
            splitter.setSizes([400, 800])
            
            main_layout.addWidget(splitter)
        
        def apply_styles(self):
            """Apply custom styles."""
            self.setStyleSheet(f"""
                QMainWindow {{
                    background-color: {COLORS['bg_primary']};
                    color: {COLORS['text_primary']};
                }}
                
                QWidget {{
                    background-color: {COLORS['bg_secondary']};
                    color: {COLORS['text_primary']};
                    font-family: 'Segoe UI', 'SF Pro Display', 'Roboto', sans-serif;
                    font-size: 13px;
                }}
                
                QGroupBox {{
                    background-color: {COLORS['bg_tertiary']};
                    border: 1px solid {COLORS['border']};
                    border-radius: 8px;
                    margin-top: 12px;
                    padding-top: 8px;
                    font-weight: bold;
                }}
                
                QTabWidget::pane {{
                    border: 1px solid {COLORS['border']};
                    background-color: {COLORS['bg_secondary']};
                }}
                
                QTabBar::tab {{
                    background-color: {COLORS['bg_widget']};
                    color: {COLORS['text_primary']};
                    padding: 8px 15px;
                    margin-right: 2px;
                    border: 1px solid {COLORS['border']};
                    border-bottom: none;
                }}
                
                QTabBar::tab:selected {{
                    background-color: {COLORS['accent']};
                    color: white;
                }}
                
                QPushButton {{
                    background-color: {COLORS['bg_widget']};
                    color: {COLORS['text_primary']};
                    border: 1px solid {COLORS['border']};
                    padding: 8px 15px;
                    border-radius: 4px;
                    font-weight: 500;
                }}
                
                QPushButton:hover {{
                    background-color: {COLORS['accent']};
                    color: white;
                }}
                
                QPushButton:disabled {{
                    background-color: {COLORS['border']};
                    color: {COLORS['text_secondary']};
                }}
                
                QLineEdit {{
                    background-color: {COLORS['bg_widget']};
                    color: {COLORS['text_primary']};
                    border: 1px solid {COLORS['border']};
                    padding: 6px 10px;
                    border-radius: 4px;
                }}
                
                QTextEdit {{
                    background-color: {COLORS['bg_widget']};
                    color: {COLORS['text_primary']};
                    border: 1px solid {COLORS['border']};
                    padding: 8px;
                    border-radius: 4px;
                }}
                
                QTableWidget {{
                    background-color: {COLORS['bg_widget']};
                    color: {COLORS['text_primary']};
                    border: 1px solid {COLORS['border']};
                    gridline-color: {COLORS['border']};
                    selection-background-color: {COLORS['accent']};
                }}
                
                QHeaderView::section {{
                    background-color: {COLORS['bg_tertiary']};
                    color: {COLORS['text_primary']};
                    border: 1px solid {COLORS['border']};
                    padding: 6px;
                    font-weight: bold;
                }}
                
                QProgressBar {{
                    border: 1px solid {COLORS['border']};
                    border-radius: 4px;
                    text-align: center;
                }}
                
                QProgressBar::chunk {{
                    background-color: {COLORS['accent']};
                    border-radius: 3px;
                }}
                
                QCheckBox {{
                    color: {COLORS['text_primary']};
                }}
                
                QRadioButton {{
                    color: {COLORS['text_primary']};
                }}
                
                QLabel {{
                    color: {COLORS['text_primary']};
                }}
            """)
        
        def load_files(self):
            """Load multiple molecules for batch processing."""
            from src.shared.qt_compat import QFileDialog
            from pychem.api import load
            filepaths, _ = QFileDialog.getOpenFileNames(
                self, "Select Molecule Files", "", "Molecule Files (*.mol *.sdf *.mol2 *.pdb);;All Files (*)"
            )
            if filepaths:
                self.molecules = []
                for fp in filepaths:
                    try:
                        mol = load(fp, parallel=False)
                        if mol:
                            if not getattr(mol, 'name', None):
                                mol.name = os.path.basename(fp)
                            self.molecules.append(mol)
                    except Exception as e:
                        print(f"Error loading {fp}: {e}")
                
                self.progress_label.setText(f"Loaded {len(self.molecules)} molecules.")
                # Update the tab placeholder
                if len(self.molecules) > 0:
                    lbl = QLabel(f"Loaded {len(self.molecules)} molecules for batch processing.")
                    lbl.setAlignment(Qt.AlignCenter)
                    self.tab_widget.removeTab(0)
                    self.tab_widget.insertTab(0, lbl, "Selection")
                    self.tab_widget.setCurrentIndex(0)

        def calculate_descriptors(self):
            """Start descriptor calculation."""
            if not self.molecules:
                QMessageBox.warning(self, "No Molecules", "Please load at least one molecule.")
                return
            
            # Start calculation thread
            self.calculate_btn.setEnabled(False)
            self.load_files_btn.setEnabled(False)
            self.stop_btn.setEnabled(True)
            self.progress_bar.setVisible(True)
            self.progress_bar.setRange(0, 0) # Indeterminate progress for batch
            self.progress_label.setText(f"Starting calculation for {len(self.molecules)} molecules...")
            
            self.calculation_thread = DescriptorCalculationThread(self.molecules, n_jobs=-1)
            
            self.calculation_thread.calculation_finished.connect(self.calculation_completed)
            self.calculation_thread.error_occurred.connect(self.calculation_error)
            
            self.calculation_thread.start()
        
        def stop_calculation(self):
            """Stop descriptor calculation."""
            if self.calculation_thread and self.calculation_thread.isRunning():
                self.calculation_thread.terminate()
                self.calculation_thread.wait()
                
                self.calculate_btn.setEnabled(True)
                self.load_files_btn.setEnabled(True)
                self.stop_btn.setEnabled(False)
                self.progress_bar.setVisible(False)
                self.progress_bar.setRange(0, 100)
                self.progress_label.setText("Calculation stopped")
        
        def update_progress(self, progress: CalculationProgress):
            """Update calculation progress."""
            self.progress_bar.setValue(int(progress.percentage))
            self.progress_label.setText(f"Calculating {progress.current_descriptor} "
                                  f"({progress.completed}/{progress.total})")
        
        def calculation_completed(self, results: Dict[str, Any]):
            """Handle calculation completion."""
            self.batch_results = results.get('batch_results', [])
            if self.batch_results:
                self.current_results = self.batch_results[0]
            
            self.calculate_btn.setEnabled(True)
            self.load_files_btn.setEnabled(True)
            self.stop_btn.setEnabled(False)
            self.progress_bar.setVisible(False)
            self.progress_bar.setRange(0, 100)
            num_desc = len(self.current_results) if self.current_results else 0
            self.progress_label.setText(f"Calculation completed: {num_desc} descriptors for {len(self.batch_results)} molecules")
            
            # Populate results table (showing first molecule as preview)
            if self.current_results:
                self.populate_results_table(self.current_results)
            
            # Enable export buttons
            self.export_results_btn.setEnabled(True)
            
            QMessageBox.information(self, "Calculation Complete", 
                               f"Successfully calculated {num_desc} descriptors for {len(self.batch_results)} molecules.\n"
                               f"Use the Export buttons to save results to CSV.")
        
        def populate_results_table(self, results: Dict[str, Any]):
            """Populate the results table with calculated descriptors for preview."""
            self.results_table.setRowCount(len(results))
            
            for row, (name, value) in enumerate(results.items()):
                # Descriptor name
                self.results_table.setItem(row, 0, QTableWidgetItem(name))
                
                # Value
                value_str = f"{value:.6f}" if isinstance(value, float) else str(value)
                self.results_table.setItem(row, 1, QTableWidgetItem(value_str))
                
                # Category 
                self.results_table.setItem(row, 2, QTableWidgetItem("PyDes"))
            
            # Resize columns to fit content
            self.results_table.resizeColumnsToContents()
        
        def get_descriptor_category(self, descriptor_name: str) -> str:
            """Get the category for a descriptor name."""
            from src.features.descriptor_calculator.descriptor_engine import DescriptorEngine
            
            engine = DescriptorEngine()
            
            for category, descriptors in engine.descriptors.items():
                for desc in descriptors:
                    if desc.name == descriptor_name:
                        return category.value
            
            return "Unknown"
        
        def calculation_error(self, error_message: str):
            """Handle calculation error."""
            self.calculate_btn.setEnabled(True)
            self.load_files_btn.setEnabled(True)
            self.stop_btn.setEnabled(False)
            self.progress_bar.setVisible(False)
            self.progress_bar.setRange(0, 100)
            self.progress_label.setText("Calculation failed")
            
            QMessageBox.critical(self, "Calculation Error", 
                              f"Error during calculation:\n{error_message}")
        
        def export_results(self):
            """Export results to CSV with user-defined location and name."""
            if not self.current_results:
                QMessageBox.warning(self, "No Results", 
                                "No calculation results to export.")
                return
            
            # Get save location from user
            from src.shared.qt_compat import QFileDialog
            default_name = f"molecular_descriptors_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            filepath, _ = QFileDialog.getSaveFileName(
                self, 
                "Export Descriptor Results",
                default_name,
                "CSV Files (*.csv);;All Files (*)"
            )
            
            if not filepath:
                return  # User cancelled
            
            try:
                self.export_results_to_csv(filepath)
                QMessageBox.information(self, "Export Complete", 
                                    f"Results exported successfully to:\n{filepath}")
            except Exception as e:
                QMessageBox.critical(self, "Export Error", 
                                  f"Error exporting results:\n{str(e)}")
        
        def export_documentation(self):
            """Export descriptor documentation to CSV with user-defined location and name."""
            # Get save location from user
            from src.shared.qt_compat import QFileDialog
            default_name = f"descriptor_documentation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            filepath, _ = QFileDialog.getSaveFileName(
                self, 
                "Export Descriptor Documentation",
                default_name,
                "CSV Files (*.csv);;All Files (*)"
            )
            
            if not filepath:
                return  # User cancelled
            
            try:
                self.export_documentation_to_csv(filepath)
                QMessageBox.information(self, "Export Complete", 
                                    f"Documentation exported successfully to:\n{filepath}")
            except Exception as e:
                QMessageBox.critical(self, "Export Error", 
                                  f"Error exporting documentation:\n{str(e)}")
        
        def export_results_to_csv(self, filepath: str):
            """Export calculation results to CSV file in transposed batch format."""
            if not self.batch_results:
                return
                
            import pandas as pd
            df = pd.DataFrame(self.batch_results)
            
            # Ensure Molecule_Name is the first column
            if 'Molecule_Name' in df.columns:
                cols = ['Molecule_Name'] + [col for col in df.columns if col != 'Molecule_Name']
                df = df[cols]
                
            df.to_csv(filepath, index=False)
        
        def get_descriptor_info(self, descriptor_name: str) -> tuple[str, str]:
            """Get description and unit for a descriptor name."""
            from src.features.descriptor_calculator.descriptor_engine import DescriptorEngine
            
            engine = DescriptorEngine()
            
            for category, descriptors in engine.descriptors.items():
                for desc in descriptors:
                    if desc.name == descriptor_name:
                        return (desc.description or '', desc.unit or '')
            
            return ('', '')
        
        def export_documentation_to_csv(self, filepath: str):
            """Export descriptor documentation to CSV file."""
            from src.features.descriptor_calculator.descriptor_engine import DescriptorEngine
            
            engine = DescriptorEngine()
            
            with open(filepath, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                
                # Write header
                writer.writerow(['Category', 'Descriptor', 'Description', 'Formula',
                                 'Unit', 'Typical Range'])

                # Write documentation for all available descriptors
                for category, descriptors in engine.descriptors.items():
                    for desc in descriptors:
                        writer.writerow([
                            category.value,
                            desc.name,
                            desc.description or '',
                            desc.formula or '',
                            desc.unit or '',
                            getattr(desc, 'range', '') or ''
                        ])

    def show_descriptor_calculator(molecule=None, parent=None):
        """Show descriptor calculator dialog."""
        try:
            dialog = DescriptorCalculatorDialog(molecule, parent)
            dialog.show()
            return dialog
        except Exception as e:
            print(f"[DEBUG GUI] Error creating dialog: {e}")
            import traceback
            traceback.print_exc()
            raise

# End of Qt availability check
else:
    # Define dummy functions when Qt is not available
    def show_descriptor_calculator(molecule, parent=None):
        """Dummy function when Qt is not available."""
        print("GUI not available: No Qt framework installed")
        print("Install PySide6 or PyQt6 to use the GUI")
        print("  pip install PySide6")
        print("  or")
        print("  pip install PyQt6")
        return None
