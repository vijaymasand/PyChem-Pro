"""
Radius of Gyration (Rg) & Radial Distribution Analyzer Plugin for PyChem-Pro

This plugin calculates:
- Center of Mass (CoM)
- Radius of Gyration (Rg)
- Normalized radial distribution of elements (C, N, O, etc.) in molecular files.

Features:
- Clean importable functions for Jupyter Notebooks or other Python scripts.
- Pure Python direct molecule parser for .mol2, .sdf, and .pdb files (no PyMOL dependency).
- Parallel processing for directories containing hundreds of structures.
- Modern PySide6 GUI with progress tracking, results table, and element distribution plots.
"""

import os
import sys
import csv
import math
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

# =========================================================================
# CORE COMPUTATIONAL LOGIC & DIRECT PARSERS (Importable & independent of GUI/PyMOL)
# =========================================================================

def get_atomic_mass(element):
    """Returns atomic mass for common elements. Fallback is Carbon mass (12.011)."""
    atomic_weights = {
        'H': 1.008, 'HE': 4.0026, 'LI': 6.94, 'BE': 9.0122, 'B': 10.81, 'C': 12.011, 
        'N': 14.007, 'O': 15.999, 'F': 18.998, 'NE': 20.180, 'NA': 22.990, 'MG': 24.305, 
        'AL': 26.982, 'SI': 28.085, 'P': 30.974, 'S': 32.06, 'CL': 35.45, 'AR': 39.948, 
        'K': 39.098, 'CA': 40.078, 'FE': 55.845, 'ZN': 65.38, 'CU': 63.546, 'BR': 79.904, 'I': 126.904
    }
    return atomic_weights.get(element.upper(), 12.011)


def parse_mol2(filepath):
    """Parses atom elements and coordinates from a .mol2 file."""
    atoms = []
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        in_atom_section = False
        for line in lines:
            if line.startswith('@<TRIPOS>ATOM'):
                in_atom_section = True
                continue
            elif line.startswith('@<TRIPOS>'):
                in_atom_section = False
                continue
                
            if in_atom_section:
                parts = line.split()
                if len(parts) >= 6:
                    atom_type = parts[5].split('.')[0]
                    element = ''.join([c for c in atom_type if c.isalpha()]).upper()
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    atoms.append((element, x, y, z))
    except Exception as e:
        print(f"Error parsing MOL2 file {filepath}: {e}")
    return atoms


def parse_pdb(filepath):
    """Parses atom elements and coordinates from a .pdb file."""
    atoms = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    element = line[76:78].strip()
                    if not element:
                        # Fallback to atom name
                        element = ''.join([c for c in line[12:16].strip() if c.isalpha()])
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atoms.append((element.upper(), x, y, z))
    except Exception as e:
        print(f"Error parsing PDB file {filepath}: {e}")
    return atoms


def parse_sdf(filepath):
    """Parses atom elements and coordinates from a .sdf file."""
    atoms = []
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        if len(lines) >= 4:
            counts_line = lines[3]
            num_atoms = int(counts_line[:3].strip())
            for idx in range(4, 4 + num_atoms):
                if idx < len(lines):
                    parts = lines[idx].split()
                    if len(parts) >= 4:
                        x = float(parts[0])
                        y = float(parts[1])
                        z = float(parts[2])
                        element = parts[3]
                        atoms.append((element.upper(), x, y, z))
    except Exception as e:
        print(f"Error parsing SDF file {filepath}: {e}")
    return atoms


def parse_molecule_file(filepath):
    """Routes file parsing to the correct handler based on extension."""
    ext = os.path.splitext(filepath)[1].lower()
    if ext == '.mol2':
        return parse_mol2(filepath)
    elif ext == '.pdb':
        return parse_pdb(filepath)
    elif ext in ('.sdf', '.mol'):
        return parse_sdf(filepath)
    return []


def calculate_com_and_rg(atoms):
    """Calculates CoM (Center of Mass) and Rg (Radius of Gyration) from atom tuples."""
    if not atoms:
        return None, None
        
    coords = np.array([[a[1], a[2], a[3]] for a in atoms])
    masses = np.array([get_atomic_mass(a[0]) for a in atoms])
    
    total_mass = np.sum(masses)
    if total_mass == 0:
        return None, None
        
    # Center of Mass
    com = np.sum(coords * masses[:, np.newaxis], axis=0) / total_mass
    
    # Radius of Gyration
    sq_dists = np.sum((coords - com)**2, axis=1)
    rg = math.sqrt(np.sum(masses * sq_dists) / total_mass)
    
    return com, rg


def analyze_molecule(filepath, bin_width=0.1):
    """
    Performs full Rg and normalized radial element analysis on a single file.
    """
    atoms = parse_molecule_file(filepath)
    com, rg = calculate_com_and_rg(atoms)
    if rg is None or rg == 0:
        return None
        
    # Element bin counts
    mol_bins = {}
    for atom in atoms:
        elem, x, y, z = atom
        dist = math.sqrt((x-com[0])**2 + (y-com[1])**2 + (z-com[2])**2)
        norm_dist = dist / rg
        bin_idx = int(math.floor(norm_dist / bin_width))
        
        if bin_idx not in mol_bins:
            mol_bins[bin_idx] = {}
        mol_bins[bin_idx][elem] = mol_bins[bin_idx].get(elem, 0) + 1
        
    return {
        'filename': os.path.basename(filepath),
        'rg': rg,
        'com': com.tolist(),
        'bins': mol_bins
    }


def _process_single_file_worker(args):
    """Worker function for multiprocessing."""
    filepath, bin_width = args
    return analyze_molecule(filepath, bin_width)


def process_folder_parallel(input_dir, extensions=('.mol2', '.sdf', '.pdb'), bin_width=0.1, num_processors=1):
    """
    Processes a directory of molecular structures in parallel.
    """
    files = [
        os.path.join(input_dir, f) for f in os.listdir(input_dir)
        if f.lower().endswith(extensions)
    ]
    if not files:
        return []
        
    tasks = [(f, bin_width) for f in files]
    
    if num_processors > 1 and len(tasks) > 1:
        try:
            with ProcessPoolExecutor(max_workers=num_processors) as executor:
                results = list(executor.map(_process_single_file_worker, tasks))
        except Exception:
            results = [_process_single_file_worker(t) for t in tasks]
    else:
        results = [_process_single_file_worker(t) for t in tasks]
        
    return [r for r in results if r is not None]

# =========================================================================
# GUI LAYER (Dynamically loaded based on Qt availability)
# =========================================================================

try:
    from src.shared.qt_compat import (
        QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
        QProgressBar, QTextEdit, QFileDialog, QMessageBox, QComboBox,
        QSpinBox, QTableWidget, QTableWidgetItem, QAbstractItemView, Qt
    )
    from src.plugins.base_plugin import BasePlugin, PluginWidget
    from src.plugins.plugin_types import PluginInfo, PluginType
    HAS_PYCHEM_BASE = True
except ImportError:
    HAS_PYCHEM_BASE = False
    try:
        from PySide6.QtWidgets import (
            QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
            QProgressBar, QTextEdit, QFileDialog, QMessageBox, QComboBox,
            QSpinBox, QTableWidget, QTableWidgetItem, QAbstractItemView, QMainWindow, QApplication
        )
        from PySide6.QtCore import Qt
    except ImportError:
        try:
            from PyQt6.QtWidgets import (
                QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
                QProgressBar, QTextEdit, QFileDialog, QMessageBox, QComboBox,
                QSpinBox, QTableWidget, QTableWidgetItem, QAbstractItemView, QMainWindow, QApplication
            )
            from PyQt6.QtCore import Qt
        except ImportError:
            class QWidget: pass
            class BasePlugin: pass
            class PluginWidget: pass

try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
except ImportError:
    try:
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
    except ImportError:
        FigureCanvas = None
        NavigationToolbar = None

from matplotlib.figure import Figure

class RadiusOfGyrationWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.input_dir = None
        self.results = []
        self.setup_ui()
        
    def setup_ui(self):
        main_layout = QHBoxLayout(self)
        
        # Left Panel (Controls)
        left_panel = QVBoxLayout()
        main_layout.addLayout(left_panel, 1)
        
        title = QLabel("Radius of Gyration (Rg) Analyzer")
        title.setStyleSheet("font-size: 16px; font-weight: bold; color: #333;")
        left_panel.addWidget(title)
        
        self.btn_select_dir = QPushButton("Select Folder with Molecules")
        self.btn_select_dir.clicked.connect(self.select_directory)
        left_panel.addWidget(self.btn_select_dir)
        
        self.lbl_dir = QLabel("No directory selected.")
        self.lbl_dir.setStyleSheet("color: gray; font-style: italic;")
        left_panel.addWidget(self.lbl_dir)
        
        # File Type Selection
        left_panel.addWidget(QLabel("File Format:"))
        self.combo_format = QComboBox()
        self.combo_format.addItems([".mol2", ".sdf", ".pdb", "All Formats"])
        left_panel.addWidget(self.combo_format)
        
        # Processor Slider
        settings_layout = QHBoxLayout()
        proc_vbox = QVBoxLayout()
        proc_vbox.addWidget(QLabel("Processors:"))
        self.spin_processors = QSpinBox()
        self.spin_processors.setRange(1, os.cpu_count() or 4)
        self.spin_processors.setValue(1)
        proc_vbox.addWidget(self.spin_processors)
        settings_layout.addLayout(proc_vbox)
        left_panel.addLayout(settings_layout)
        
        self.btn_run = QPushButton("Process Directory")
        self.btn_run.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold; padding: 6px;")
        self.btn_run.clicked.connect(self.run_processing)
        left_panel.addWidget(self.btn_run)
        
        self.btn_save = QPushButton("Save Results (CSV)")
        self.btn_save.setEnabled(False)
        self.btn_save.clicked.connect(self.save_results)
        left_panel.addWidget(self.btn_save)
        
        # Table of results
        self.table_results = QTableWidget()
        self.table_results.setColumnCount(3)
        self.table_results.setHorizontalHeaderLabels(["Filename", "Rg (Å)", "Center of Mass"])
        left_panel.addWidget(self.table_results)
        
        # Right Panel (Plot Canvas)
        if FigureCanvas is not None:
            self.plot_layout = QVBoxLayout()
            self.figure = Figure(figsize=(5, 4))
            self.canvas = FigureCanvas(self.figure)
            self.toolbar = NavigationToolbar(self.canvas, self)
            self.plot_layout.addWidget(self.toolbar)
            self.plot_layout.addWidget(self.canvas)
            main_layout.addLayout(self.plot_layout, 1)
        else:
            main_layout.addWidget(QLabel("Matplotlib canvas could not be loaded."), 1)

    def select_directory(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Input Directory")
        if dir_path:
            self.input_dir = dir_path
            self.lbl_dir.setText(f"Folder: {os.path.basename(dir_path)}")

    def run_processing(self):
        if not self.input_dir:
            QMessageBox.warning(self, "Warning", "Please select an input directory.")
            return
            
        fmt = self.combo_format.currentText()
        exts = ('.mol2', '.sdf', '.pdb') if fmt == "All Formats" else (fmt,)
        cores = self.spin_processors.value()
        
        self.results = process_folder_parallel(self.input_dir, exts, 0.1, cores)
        
        if not self.results:
            QMessageBox.warning(self, "No Files", "No molecules were processed.")
            return
            
        # Update Table
        self.table_results.setRowCount(len(self.results))
        for r_idx, res in enumerate(self.results):
            self.table_results.setItem(r_idx, 0, QTableWidgetItem(res['filename']))
            self.table_results.setItem(r_idx, 1, QTableWidgetItem(f"{res['rg']:.3f}"))
            com_str = ", ".join([f"{c:.2f}" for c in res['com']])
            self.table_results.setItem(r_idx, 2, QTableWidgetItem(f"[{com_str}]"))
            
        self.btn_save.setEnabled(True)
        self.update_distribution_plot()

    def update_distribution_plot(self):
        if FigureCanvas is None or not self.results:
            return
            
        # Aggregate distributions
        global_bins = {}
        all_elements = set()
        total_mols = len(self.results)
        
        for res in self.results:
            for b_idx, elems in res['bins'].items():
                if b_idx not in global_bins:
                    global_bins[b_idx] = {}
                for elem, count in elems.items():
                    global_bins[b_idx][elem] = global_bins[b_idx].get(elem, 0) + count
                    all_elements.add(elem)
                    
        if not global_bins:
            return
            
        max_bin = max(global_bins.keys())
        sorted_elements = sorted(list(all_elements))
        
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        
        bin_width = 0.1
        x_vals = np.arange(max_bin + 1) * bin_width
        
        for elem in sorted_elements:
            y_vals = []
            for b_idx in range(max_bin + 1):
                total_cnt = global_bins.get(b_idx, {}).get(elem, 0)
                y_vals.append(total_cnt / total_mols)
            ax.plot(x_vals, y_vals, marker='o', label=f'Avg {elem}')
            
        ax.set_title("Averaged Radial Distribution of Atoms")
        ax.set_xlabel("Normalized Distance (r/Rg)")
        ax.set_ylabel("Average Count per Molecule")
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.5)
        self.canvas.draw()

    def save_results(self):
        if not self.results:
            return
            
        base_path, _ = QFileDialog.getSaveFileName(self, "Save Output Base Filename", "", "CSV Files (*.csv)")
        if not base_path:
            return
            
        base_name = os.path.splitext(base_path)[0]
        rg_csv = f"{base_name}_rg.csv"
        dist_csv = f"{base_name}_distribution.csv"
        
        # 1. Save Rg & CoM
        with open(rg_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Filename", "Rg", "CoM_X", "CoM_Y", "CoM_Z"])
            for res in self.results:
                writer.writerow([res['filename'], f"{res['rg']:.3f}"] + res['com'])
                
        # 2. Save Distribution
        global_bins = {}
        all_elements = set()
        total_mols = len(self.results)
        
        for res in self.results:
            for b_idx, elems in res['bins'].items():
                if b_idx not in global_bins:
                    global_bins[b_idx] = {}
                for elem, count in elems.items():
                    global_bins[b_idx][elem] = global_bins[b_idx].get(elem, 0) + count
                    all_elements.add(elem)
                    
        max_bin = max(global_bins.keys()) if global_bins else 0
        sorted_elements = sorted(list(all_elements))
        
        with open(dist_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Normalized_Bin_Start", "Normalized_Bin_End"] + [f"Avg_{e}" for e in sorted_elements])
            for i in range(max_bin + 1):
                start = i * 0.1
                end = (i + 1) * 0.1
                row = [f"{start:.1f}", f"{end:.1f}"]
                for elem in sorted_elements:
                    avg_val = global_bins.get(i, {}).get(elem, 0) / total_mols
                    row.append(f"{avg_val:.4f}")
                writer.writerow(row)
                
        QMessageBox.information(self, "Success", f"Saved:\n1. {rg_csv}\n2. {dist_csv}")


# =========================================================================
# PLUGIN INTERFACE (Integrates into PyChem-Pro if imported)
# =========================================================================

class RadiusOfGyrationPlugin(BasePlugin):
    def __init__(self):
        super().__init__(PluginInfo(
            name="Radius of Gyration Analyzer",
            version="1.1.0",
            description="Calculates Radius of Gyration and Atom Distributions in parallel.",
            author="DeepMind Antigravity",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        ))
        self.widget = None

    def create_widget(self):
        if self.widget is None:
            self.widget = RadiusOfGyrationWidget()
        return self.widget

    def initialize(self):
        return True

    def cleanup(self):
        if self.widget:
            self.widget.deleteLater()
            self.widget = None


# Standalone runner
if __name__ == "__main__":
    try:
        app = QApplication(sys.argv)
        window = QMainWindow()
        window.setWindowTitle("Radius of Gyration Analyzer")
        window.setCentralWidget(RadiusOfGyrationWidget())
        window.resize(1100, 700)
        window.show()
        sys.exit(app.exec())
    except Exception as e:
        print(f"Failed to start GUI: {e}")