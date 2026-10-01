"""
Molecular Descriptor Calculator GUI

Comprehensive GUI for molecular descriptor calculation with:
- Multiprocessing via imap_unordered (per-molecule parallelism)
- Streaming CSV writes: each molecule's row is flushed to disk immediately
  after its descriptors are computed — no end-of-batch write
- Live progress: determinate progress bar, molecule counter, elapsed time, ETA
- Live preview table: updates with each completed molecule
- Clean Stop: partial CSV is preserved
"""

import sys
import csv
import os
import time
import threading
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

    # -----------------------------------------------------------------------
    # Background calculation thread
    # -----------------------------------------------------------------------

    class DescriptorCalculationThread(QThread):
        """
        Background thread that:
          1. Opens the output CSV and writes the header row.
          2. Dispatches molecules to a multiprocessing pool via imap_unordered.
          3. For every completed molecule: writes its CSV row (flushed immediately)
             and emits molecule_done so the GUI can update live.
          4. Emits calculation_finished when all molecules are processed
             (or calculation_stopped if the user cancelled mid-run).
        """

        # (result_dict, completed_count, total_count) — fired per molecule
        molecule_done  = Signal(dict, int, int)
        # Fired once at the end — carries (total_written, output_csv_path, was_stopped)
        calculation_finished = Signal(int, str, bool)
        error_occurred = Signal(str)

        def __init__(self, molecules, output_path: str, n_jobs: int = -1, parent=None):
            super().__init__(parent)
            self.molecules = molecules if isinstance(molecules, list) else [molecules]
            self.output_path = output_path
            self.n_jobs = n_jobs
            # Threading event used to signal an early stop from the main thread
            self._stop_event = threading.Event()

        def request_stop(self):
            """Ask the thread to stop after the current molecule finishes."""
            self._stop_event.set()

        def run(self):
            total = len(self.molecules)
            if total == 0:
                self.calculation_finished.emit(0, self.output_path, False)
                return

            csv_file = None
            csv_writer = None
            header_written = False
            completed = 0

            try:
                # Open the CSV once for the whole run
                csv_file = open(self.output_path, 'w', newline='', encoding='utf-8')

                def _on_molecule_result(result: dict, mol_idx: int, mol_total: int):
                    nonlocal csv_writer, header_written, completed

                    # Write header on first result (we now know all column names)
                    if not header_written:
                        fieldnames = list(result.keys())
                        # Ensure Molecule_Name is first
                        if 'Molecule_Name' in fieldnames:
                            fieldnames = ['Molecule_Name'] + [
                                f for f in fieldnames if f != 'Molecule_Name'
                            ]
                        csv_writer = csv.DictWriter(
                            csv_file, fieldnames=fieldnames, extrasaction='ignore'
                        )
                        csv_writer.writeheader()
                        csv_file.flush()
                        header_written = True

                    # Write this molecule's row and flush immediately
                    csv_writer.writerow(result)
                    csv_file.flush()
                    completed += 1

                    # Notify GUI (thread-safe via Qt signal)
                    self.molecule_done.emit(result, completed, total)

                # Run the streaming calculation
                PyDesEngine.calculate_batch_streaming(
                    molecules=self.molecules,
                    result_callback=_on_molecule_result,
                    n_jobs=self.n_jobs,
                    stop_event=self._stop_event,
                )

                was_stopped = self._stop_event.is_set()
                self.calculation_finished.emit(completed, self.output_path, was_stopped)

            except Exception as e:
                import traceback as _tb
                self.error_occurred.emit(f"{e}\n\n{_tb.format_exc()}")
            finally:
                if csv_file:
                    csv_file.close()

    # -----------------------------------------------------------------------
    # Selection builder (unchanged from original)
    # -----------------------------------------------------------------------

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
            self.update_selection_display()

        def get_selection(self) -> AtomSelection:
            selected_type = self.selection_type_group.checkedId()
            if selected_type == 0:
                return AtomSelection(
                    SelectionType.ALL,
                    list(range(len(self.molecule.atoms))),
                    "All atoms"
                )
            elif selected_type == 1:
                indices = self.parse_custom_selection()
                return AtomSelection(
                    SelectionType.CUSTOM,
                    indices,
                    f"Custom: {self.custom_input.text()}"
                )
            return AtomSelection(SelectionType.ALL, [], "Empty selection")

        def parse_custom_selection(self) -> List[int]:
            text = self.custom_input.text().strip()
            if not text:
                return []
            indices = []
            for part in text.split(','):
                part = part.strip()
                if '-' in part:
                    start, end = map(int, part.split('-'))
                    indices.extend(range(start, end + 1))
                else:
                    indices.append(int(part))
            return sorted(set(indices))

        def preview_selection(self):
            self.update_selection_display()

        def clear_selection(self):
            self.custom_input.clear()
            self.all_radio.setChecked(True)
            self.update_selection_display()

        def update_selection_display(self):
            selection = self.get_selection()
            if len(selection.atom_indices) <= 20:
                indices_str = ", ".join(map(str, selection.atom_indices))
            else:
                indices_str = f"{len(selection.atom_indices)} atoms: {selection.atom_indices[:10]}..."
            display_text = (
                f"Type: {selection.selection_type.value}\n"
                f"Description: {selection.description}\n"
                f"Indices: {indices_str}"
            )
            self.current_selection.setText(display_text)

    # -----------------------------------------------------------------------
    # Descriptor category config widget (unchanged)
    # -----------------------------------------------------------------------

    class DescriptorConfigWidget(QWidget):
        """Widget for configuring descriptor categories."""

        def __init__(self, parent=None):
            super().__init__(parent)
            self.engine = DescriptorEngine()
            self.init_ui()

        def init_ui(self):
            layout = QVBoxLayout(self)

            category_group = QGroupBox("Descriptor Categories")
            category_layout = QVBoxLayout(category_group)

            self.category_checkboxes = {}
            for category in DescriptorCategory:
                checkbox = QCheckBox(
                    f"{category.value} ({len(self.engine.descriptors[category])} descriptors)"
                )
                checkbox.setChecked(True)
                self.category_checkboxes[category] = checkbox
                category_layout.addWidget(checkbox)

            layout.addWidget(category_group)

            button_layout = QHBoxLayout()
            self.select_all_btn = QPushButton("Select All")
            self.select_all_btn.clicked.connect(self.select_all_categories)
            button_layout.addWidget(self.select_all_btn)

            self.select_none_btn = QPushButton("Select None")
            self.select_none_btn.clicked.connect(self.select_none_categories)
            button_layout.addWidget(self.select_none_btn)

            layout.addLayout(button_layout)
            self.update_descriptor_preview()

            for checkbox in self.category_checkboxes.values():
                checkbox.toggled.connect(self.update_descriptor_preview)

        def get_selected_categories(self) -> List[DescriptorCategory]:
            return [cat for cat, cb in self.category_checkboxes.items() if cb.isChecked()]

        def select_all_categories(self):
            for cb in self.category_checkboxes.values():
                cb.setChecked(True)

        def select_none_categories(self):
            for cb in self.category_checkboxes.values():
                cb.setChecked(False)

        def update_descriptor_preview(self):
            selected = self.get_selected_categories()
            total = sum(len(self.engine.descriptors[c]) for c in selected)
            print(
                f"Selected {len(selected)} categories — {total} descriptors total"
            )

    # -----------------------------------------------------------------------
    # Main dialog
    # -----------------------------------------------------------------------

    class DescriptorCalculatorDialog(QMainWindow):
        """
        Main descriptor calculator dialog.

        Workflow:
          1. User loads molecules (or a single molecule is passed in).
          2. User picks an output CSV path via the "Browse…" button.
          3. User clicks "Calculate" — the thread starts immediately.
          4. Progress bar and counters update molecule-by-molecule.
          5. The preview table shows the most recently finished molecule.
          6. "Stop" cleanly terminates the pool; partial CSV is preserved.
          7. "Open Output File" opens the completed (or partial) CSV.
        """

        def __init__(self, molecule=None, parent=None):
            super().__init__(parent)
            self.molecules = [molecule] if molecule else []
            self.calculation_thread = None
            self._calc_start_time = None
            self._last_result = {}          # most recently completed molecule
            self._completed_count = 0

            self.setWindowTitle("Molecular Descriptor Calculator (PyDes)")
            self.setGeometry(100, 100, 1280, 820)

            self.init_ui()
            self.apply_styles()

        # ------------------------------------------------------------------
        # UI construction
        # ------------------------------------------------------------------

        def init_ui(self):
            central_widget = QWidget()
            self.setCentralWidget(central_widget)
            main_layout = QHBoxLayout(central_widget)

            # ---- Left panel ------------------------------------------------
            left_panel = QWidget()
            left_layout = QVBoxLayout(left_panel)

            self.tab_widget = QTabWidget()

            # Molecule info tab
            if self.molecules and len(self.molecules) == 1:
                self.selection_widget = SelectionBuilder(self.molecules[0])
                self.tab_widget.addTab(self.selection_widget, "Selection")
            else:
                lbl = QLabel(f"Loaded {len(self.molecules)} molecules for batch processing.")
                lbl.setAlignment(Qt.AlignCenter)
                self.tab_widget.addTab(lbl, "Selection")

            # Configuration tab
            self.config_widget = DescriptorConfigWidget()
            self.tab_widget.addTab(self.config_widget, "Configuration")

            left_layout.addWidget(self.tab_widget)

            # ---- Right panel -----------------------------------------------
            right_panel = QWidget()
            right_layout = QVBoxLayout(right_panel)

            # Output path row
            output_group = QGroupBox("Output File")
            output_layout = QHBoxLayout(output_group)
            self.output_path_edit = QLineEdit()
            default_name = (
                f"molecular_descriptors_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            )
            self.output_path_edit.setText(
                os.path.join(os.path.expanduser("~"), "Documents", default_name)
            )
            self.output_path_edit.setPlaceholderText("Output CSV file path…")
            output_layout.addWidget(self.output_path_edit)

            self.browse_btn = QPushButton("Browse…")
            self.browse_btn.clicked.connect(self.browse_output_path)
            output_layout.addWidget(self.browse_btn)

            right_layout.addWidget(output_group)

            # Stats bar (counters + time)
            stats_group = QGroupBox("Run Statistics")
            stats_layout = QHBoxLayout(stats_group)

            self.counter_label = QLabel("Molecules: 0 / 0")
            self.counter_label.setObjectName("stats_label")
            stats_layout.addWidget(self.counter_label)

            stats_layout.addStretch()

            self.elapsed_label = QLabel("Elapsed: —")
            self.elapsed_label.setObjectName("stats_label")
            stats_layout.addWidget(self.elapsed_label)

            stats_layout.addStretch()

            self.eta_label = QLabel("ETA: —")
            self.eta_label.setObjectName("stats_label")
            stats_layout.addWidget(self.eta_label)

            right_layout.addWidget(stats_group)

            # Timer for clock updates
            self._timer = QTimer(self)
            self._timer.setInterval(1000)
            self._timer.timeout.connect(self._update_time_labels)

            # Progress bar
            self.progress_bar = QProgressBar()
            self.progress_bar.setVisible(False)
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)
            right_layout.addWidget(self.progress_bar)

            self.progress_label = QLabel("Ready — choose an output file and click Calculate")
            right_layout.addWidget(self.progress_label)

            # Live preview table
            results_group = QGroupBox("Live Preview (most recent molecule)")
            results_layout = QVBoxLayout(results_group)

            self.results_table = QTableWidget()
            self.results_table.setColumnCount(3)
            self.results_table.setHorizontalHeaderLabels(["Descriptor", "Value", "Category"])
            self.results_table.horizontalHeader().setStretchLastSection(True)
            results_layout.addWidget(self.results_table)

            right_layout.addWidget(results_group)

            # Controls row
            controls_group = QGroupBox("Calculation Controls")
            controls_layout = QVBoxLayout(controls_group)

            button_layout = QHBoxLayout()

            self.load_files_btn = QPushButton("Load Files…")
            self.load_files_btn.clicked.connect(self.load_files)
            button_layout.addWidget(self.load_files_btn)

            self.calculate_btn = QPushButton("⚡  Calculate Descriptors")
            self.calculate_btn.clicked.connect(self.calculate_descriptors)
            button_layout.addWidget(self.calculate_btn)

            self.stop_btn = QPushButton("⏹  Stop")
            self.stop_btn.clicked.connect(self.stop_calculation)
            self.stop_btn.setEnabled(False)
            button_layout.addWidget(self.stop_btn)

            self.open_output_btn = QPushButton("📂  Open Output File")
            self.open_output_btn.clicked.connect(self.open_output_file)
            self.open_output_btn.setEnabled(False)
            button_layout.addWidget(self.open_output_btn)

            self.export_documentation_btn = QPushButton("📄  Export Documentation")
            self.export_documentation_btn.clicked.connect(self.export_documentation)
            button_layout.addWidget(self.export_documentation_btn)

            controls_layout.addLayout(button_layout)
            right_layout.addWidget(controls_group)

            # Assemble splitter
            splitter = QSplitter(Qt.Horizontal)
            splitter.addWidget(left_panel)
            splitter.addWidget(right_panel)
            splitter.setSizes([380, 900])

            main_layout.addWidget(splitter)

        # ------------------------------------------------------------------
        # Styles
        # ------------------------------------------------------------------

        def apply_styles(self):
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
                QGroupBox::title {{
                    subcontrol-origin: margin;
                    left: 10px;
                    padding: 0 4px;
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
                    min-height: 18px;
                }}
                QProgressBar::chunk {{
                    background-color: {COLORS['accent']};
                    border-radius: 3px;
                }}
                QCheckBox {{ color: {COLORS['text_primary']}; }}
                QRadioButton {{ color: {COLORS['text_primary']}; }}
                QLabel {{ color: {COLORS['text_primary']}; }}
                QLabel#stats_label {{
                    font-size: 12px;
                    color: {COLORS['text_secondary']};
                    padding: 2px 6px;
                }}
            """)

        # ------------------------------------------------------------------
        # File operations
        # ------------------------------------------------------------------

        def browse_output_path(self):
            """Let the user pick the output CSV location."""
            from src.shared.qt_compat import QFileDialog
            default_name = (
                f"molecular_descriptors_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            )
            filepath, _ = QFileDialog.getSaveFileName(
                self,
                "Choose Output CSV File",
                os.path.join(os.path.expanduser("~"), "Documents", default_name),
                "CSV Files (*.csv);;All Files (*)"
            )
            if filepath:
                self.output_path_edit.setText(filepath)

        def load_files(self):
            """Load multiple molecules for batch processing."""
            from src.shared.qt_compat import QFileDialog
            from pychem.api import load
            filepaths, _ = QFileDialog.getOpenFileNames(
                self, "Select Molecule Files", "",
                "Molecule Files (*.mol *.sdf *.mol2 *.pdb);;All Files (*)"
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

                self.progress_label.setText(
                    f"Loaded {len(self.molecules)} molecules. Ready to calculate."
                )
                self.counter_label.setText(f"Molecules: 0 / {len(self.molecules)}")

                # Refresh the Selection tab
                if self.molecules:
                    lbl = QLabel(
                        f"Loaded {len(self.molecules)} molecules for batch processing."
                    )
                    lbl.setAlignment(Qt.AlignCenter)
                    self.tab_widget.removeTab(0)
                    self.tab_widget.insertTab(0, lbl, "Selection")
                    self.tab_widget.setCurrentIndex(0)

        def open_output_file(self):
            """Open the output CSV in the OS default application."""
            path = self.output_path_edit.text().strip()
            if path and os.path.exists(path):
                import subprocess
                try:
                    if sys.platform == 'win32':
                        os.startfile(path)
                    elif sys.platform == 'darwin':
                        subprocess.Popen(['open', path])
                    else:
                        subprocess.Popen(['xdg-open', path])
                except Exception as e:
                    QMessageBox.warning(self, "Cannot Open File",
                                        f"Could not open the file:\n{e}")
            else:
                QMessageBox.warning(self, "File Not Found",
                                    "Output file does not exist yet.")

        # ------------------------------------------------------------------
        # Calculation lifecycle
        # ------------------------------------------------------------------

        def calculate_descriptors(self):
            """Validate inputs and start the streaming calculation thread."""
            if not self.molecules:
                QMessageBox.warning(self, "No Molecules",
                                    "Please load at least one molecule.")
                return

            output_path = self.output_path_edit.text().strip()
            if not output_path:
                QMessageBox.warning(self, "No Output File",
                                    "Please choose an output CSV file path first.")
                return

            # Ensure the parent directory exists
            out_dir = os.path.dirname(output_path)
            if out_dir and not os.path.exists(out_dir):
                try:
                    os.makedirs(out_dir, exist_ok=True)
                except OSError as e:
                    QMessageBox.critical(self, "Cannot Create Directory",
                                         f"Cannot create output directory:\n{e}")
                    return

            # Reset state
            self._completed_count = 0
            self._last_result = {}
            self._calc_start_time = time.monotonic()

            # Update UI
            self.calculate_btn.setEnabled(False)
            self.load_files_btn.setEnabled(False)
            self.browse_btn.setEnabled(False)
            self.stop_btn.setEnabled(True)
            self.open_output_btn.setEnabled(False)
            self.progress_bar.setVisible(True)
            self.progress_bar.setRange(0, len(self.molecules))
            self.progress_bar.setValue(0)
            self.counter_label.setText(f"Molecules: 0 / {len(self.molecules)}")
            self.elapsed_label.setText("Elapsed: 0s")
            self.eta_label.setText("ETA: calculating…")
            self.progress_label.setText(
                f"Starting calculation for {len(self.molecules)} molecules…"
            )
            self.results_table.setRowCount(0)

            self._timer.start()

            # Build and start the thread
            self.calculation_thread = DescriptorCalculationThread(
                molecules=self.molecules,
                output_path=output_path,
                n_jobs=-1,
            )
            self.calculation_thread.molecule_done.connect(self._on_molecule_done)
            self.calculation_thread.calculation_finished.connect(self._on_calculation_finished)
            self.calculation_thread.error_occurred.connect(self._on_calculation_error)
            self.calculation_thread.start()

        def stop_calculation(self):
            """Request a clean stop — partial CSV is preserved."""
            if self.calculation_thread and self.calculation_thread.isRunning():
                self.calculation_thread.request_stop()
                self.stop_btn.setEnabled(False)
                self.progress_label.setText("Stop requested — finishing current molecule…")

        # ------------------------------------------------------------------
        # Slots: per-molecule updates
        # ------------------------------------------------------------------

        def _on_molecule_done(self, result: dict, completed: int, total: int):
            """Called after each molecule completes — updates all live widgets."""
            self._completed_count = completed
            self._last_result = result

            # Progress bar
            self.progress_bar.setValue(completed)

            # Counter
            self.counter_label.setText(f"Molecules: {completed} / {total}")

            # Molecule name
            mol_name = result.get('Molecule_Name', f'Molecule {completed}')
            self.progress_label.setText(
                f"Completed: {mol_name}  ({completed} / {total})"
            )

            # Live preview table — show this molecule's descriptors
            self._populate_results_table(result)

        def _update_time_labels(self):
            """Fired every second by QTimer to refresh Elapsed / ETA."""
            if self._calc_start_time is None:
                return
            elapsed = time.monotonic() - self._calc_start_time
            total = len(self.molecules)
            done = self._completed_count

            self.elapsed_label.setText(f"Elapsed: {self._fmt_duration(elapsed)}")

            if done > 0 and total > 0:
                rate = elapsed / done          # seconds per molecule
                remaining = rate * (total - done)
                self.eta_label.setText(f"ETA: {self._fmt_duration(remaining)}")
            else:
                self.eta_label.setText("ETA: calculating…")

        @staticmethod
        def _fmt_duration(seconds: float) -> str:
            seconds = int(seconds)
            if seconds < 60:
                return f"{seconds}s"
            m, s = divmod(seconds, 60)
            if m < 60:
                return f"{m}m {s}s"
            h, m = divmod(m, 60)
            return f"{h}h {m}m"

        # ------------------------------------------------------------------
        # Slots: completion / error
        # ------------------------------------------------------------------

        def _on_calculation_finished(self, written: int, csv_path: str, was_stopped: bool):
            """Fired when the thread exits (completed or stopped)."""
            self._timer.stop()
            self._reset_controls()

            elapsed = time.monotonic() - self._calc_start_time if self._calc_start_time else 0
            self.elapsed_label.setText(f"Elapsed: {self._fmt_duration(elapsed)}")
            self.eta_label.setText("—")

            self.progress_bar.setValue(written)

            if was_stopped:
                self.progress_label.setText(
                    f"Stopped — {written} molecules written to CSV (partial results preserved)."
                )
                QMessageBox.information(
                    self, "Calculation Stopped",
                    f"Calculation was stopped.\n"
                    f"{written} molecules were written to:\n{csv_path}\n\n"
                    f"The partial CSV is complete and usable."
                )
            else:
                self.progress_label.setText(
                    f"Done — {written} molecules written in {self._fmt_duration(elapsed)}."
                )
                QMessageBox.information(
                    self, "Calculation Complete",
                    f"Successfully calculated descriptors for {written} molecules.\n\n"
                    f"Results saved to:\n{csv_path}"
                )

            if written > 0 and os.path.exists(csv_path):
                self.open_output_btn.setEnabled(True)

        def _on_calculation_error(self, error_message: str):
            """Handle calculation error."""
            self._timer.stop()
            self._reset_controls()
            self.progress_label.setText("Calculation failed — see error details.")
            QMessageBox.critical(self, "Calculation Error",
                                  f"Error during calculation:\n\n{error_message}")

        def _reset_controls(self):
            """Re-enable UI controls after a run ends."""
            self.calculate_btn.setEnabled(True)
            self.load_files_btn.setEnabled(True)
            self.browse_btn.setEnabled(True)
            self.stop_btn.setEnabled(False)
            self.progress_bar.setVisible(False)

        # ------------------------------------------------------------------
        # Results table
        # ------------------------------------------------------------------

        def _populate_results_table(self, results: Dict[str, Any]):
            """Fill the preview table with the most recently completed molecule."""
            # Exclude Molecule_Name from the descriptor rows
            rows = [(k, v) for k, v in results.items() if k != 'Molecule_Name']
            self.results_table.setRowCount(len(rows))

            for row, (name, value) in enumerate(rows):
                self.results_table.setItem(row, 0, QTableWidgetItem(name))
                value_str = f"{value:.6f}" if isinstance(value, float) else str(value)
                self.results_table.setItem(row, 1, QTableWidgetItem(value_str))
                self.results_table.setItem(row, 2, QTableWidgetItem("PyDes"))

            self.results_table.resizeColumnsToContents()

        # ------------------------------------------------------------------
        # Documentation export (unchanged in purpose)
        # ------------------------------------------------------------------

        def export_documentation(self):
            """Export descriptor documentation to CSV."""
            from src.shared.qt_compat import QFileDialog
            default_name = (
                f"descriptor_documentation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            )
            filepath, _ = QFileDialog.getSaveFileName(
                self,
                "Export Descriptor Documentation",
                default_name,
                "CSV Files (*.csv);;All Files (*)"
            )
            if not filepath:
                return
            try:
                self._write_documentation_csv(filepath)
                QMessageBox.information(self, "Export Complete",
                                        f"Documentation exported to:\n{filepath}")
            except Exception as e:
                QMessageBox.critical(self, "Export Error",
                                     f"Error exporting documentation:\n{str(e)}")

        def _write_documentation_csv(self, filepath: str):
            engine = DescriptorEngine()
            with open(filepath, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(['Category', 'Descriptor', 'Description',
                                  'Formula', 'Unit', 'Typical Range'])
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

        # ------------------------------------------------------------------
        # Legacy helpers (kept for any existing callers)
        # ------------------------------------------------------------------

        def get_descriptor_category(self, descriptor_name: str) -> str:
            engine = DescriptorEngine()
            for category, descriptors in engine.descriptors.items():
                for desc in descriptors:
                    if desc.name == descriptor_name:
                        return category.value
            return "Unknown"

        def get_descriptor_info(self, descriptor_name: str) -> tuple:
            engine = DescriptorEngine()
            for category, descriptors in engine.descriptors.items():
                for desc in descriptors:
                    if desc.name == descriptor_name:
                        return (desc.description or '', desc.unit or '')
            return ('', '')

    # -----------------------------------------------------------------------
    # Public factory
    # -----------------------------------------------------------------------

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
    def show_descriptor_calculator(molecule=None, parent=None):
        """Dummy function when Qt is not available."""
        print("GUI not available: No Qt framework installed")
        print("Install PySide6 or PyQt6 to use the GUI")
        print("  pip install PySide6")
        print("  or")
        print("  pip install PyQt6")
        return None
