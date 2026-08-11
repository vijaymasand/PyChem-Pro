"""
Chemical Space Statistical Comparator Plugin for PyChem-Pro

This plugin compares multiple datasets against a reference dataset using:
- Kolmogorov-Smirnov (KS) Test
- Kullback-Leibler (KL) Divergence
- Significance scoring

Features:
- Clean importable functions for Jupyter Notebooks or other Python scripts.
- Parallelized comparison testing.
- Sleek PySide6 GUI with interactive results table and overlay distribution plotting.
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
from concurrent.futures import ProcessPoolExecutor

# =========================================================================
# CORE COMPUTATIONAL LOGIC (Importable & independent of GUI)
# =========================================================================

def calculate_ks_kl(ref_data, comp_data, bins=50):
    """
    Calculates KS statistic, p-value, and KL divergence between two datasets.
    """
    ref_data = ref_data[~np.isnan(ref_data)]
    comp_data = comp_data[~np.isnan(comp_data)]
    
    if len(ref_data) < 2 or len(comp_data) < 2:
        return np.nan, np.nan, np.nan
        
    # 1. KS Test
    ks_stat, p_value = stats.ks_2samp(ref_data, comp_data)
    
    # 2. KL Divergence
    min_val = min(np.min(ref_data), np.min(comp_data))
    max_val = max(np.max(ref_data), np.max(comp_data))
    
    if min_val == max_val:
        return ks_stat, p_value, 0.0
        
    bin_edges = np.linspace(min_val, max_val, bins)
    
    hist_ref, _ = np.histogram(ref_data, bins=bin_edges, density=True)
    hist_comp, _ = np.histogram(comp_data, bins=bin_edges, density=True)
    
    # Laplace smoothing
    epsilon = 1e-10
    hist_ref += epsilon
    hist_comp += epsilon
    
    hist_ref /= np.sum(hist_ref)
    hist_comp /= np.sum(hist_comp)
    
    kl_div = stats.entropy(pk=hist_ref, qk=hist_comp)
    
    return ks_stat, p_value, kl_div


def _compare_dataset_worker(args):
    """Worker function for parallel dataset comparisons."""
    comp_name, comp_data, ref_data, col_name = args
    ks, p, kl = calculate_ks_kl(ref_data, comp_data)
    
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    
    return {
        'Comparison': comp_name,
        'Feature': col_name,
        'Sample Sizes': f"n={len(ref_data)} vs n={len(comp_data)}",
        'KS Statistic': ks,
        'P-Value': p,
        'Significance': sig,
        'KL Divergence': kl
    }


def compare_datasets_parallel(reference_df, comparison_dfs, column_name, num_processors=1):
    """
    Compares reference dataset with multiple target datasets on a selected column.
    
    Args:
        reference_df (pd.DataFrame): Reference dataset.
        comparison_dfs (dict): Dict of dataset name -> DataFrame.
        column_name (str): Column to run the statistical comparison on.
        num_processors (int): Number of CPU cores to utilize.
        
    Returns:
        pd.DataFrame: Table containing statistical comparison results.
    """
    if column_name not in reference_df.columns:
        return pd.DataFrame()
        
    ref_data = reference_df[column_name].dropna().values
    if len(ref_data) == 0:
        return pd.DataFrame()
        
    tasks = []
    for name, df in comparison_dfs.items():
        if column_name in df.columns:
            comp_data = df[column_name].dropna().values
            if len(comp_data) > 0:
                tasks.append((name, comp_data, ref_data, column_name))
                
    if not tasks:
        return pd.DataFrame()
        
    if num_processors > 1 and len(tasks) > 1:
        try:
            with ProcessPoolExecutor(max_workers=num_processors) as executor:
                results = list(executor.map(_compare_dataset_worker, tasks))
        except Exception:
            results = [_compare_dataset_worker(t) for t in tasks]
    else:
        results = [_compare_dataset_worker(t) for t in tasks]
        
    return pd.DataFrame(results)

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
            # Stubs
            class QWidget: pass
            class BasePlugin: pass
            class PluginWidget: pass

# Matplotlib Qt Integration
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

class StatComparisonWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.loaded_datasets = {} # {filename: dataframe}
        self.common_columns = []
        self.setup_ui()
        
    def setup_ui(self):
        main_layout = QHBoxLayout(self)
        
        # Left Panel (Controls)
        left_panel = QVBoxLayout()
        main_layout.addLayout(left_panel, 1)
        
        title = QLabel("Chemical Space Statistical Comparator")
        title.setStyleSheet("font-size: 16px; font-weight: bold; color: #333;")
        left_panel.addWidget(title)
        
        # Load button
        self.btn_load = QPushButton("Load CSV Files")
        self.btn_load.clicked.connect(self.load_files)
        left_panel.addWidget(self.btn_load)
        
        self.lbl_status = QLabel("No files loaded.")
        self.lbl_status.setStyleSheet("color: gray; font-style: italic;")
        left_panel.addWidget(self.lbl_status)
        
        # Column selection dropdown
        left_panel.addWidget(QLabel("Analysis Column:"))
        self.combo_columns = QComboBox()
        left_panel.addWidget(self.combo_columns)
        
        # Reference dataset dropdown
        left_panel.addWidget(QLabel("Reference Dataset:"))
        self.combo_ref = QComboBox()
        left_panel.addWidget(self.combo_ref)
        
        # Processors and settings
        settings_layout = QHBoxLayout()
        proc_vbox = QVBoxLayout()
        proc_vbox.addWidget(QLabel("Processors:"))
        self.spin_processors = QSpinBox()
        self.spin_processors.setRange(1, os.cpu_count() or 4)
        self.spin_processors.setValue(1)
        proc_vbox.addWidget(self.spin_processors)
        settings_layout.addLayout(proc_vbox)
        left_panel.addLayout(settings_layout)
        
        # Run button
        self.btn_run = QPushButton("Run Comparison")
        self.btn_run.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold; padding: 6px;")
        self.btn_run.clicked.connect(self.run_comparison)
        left_panel.addWidget(self.btn_run)
        
        # Save results button
        self.btn_save = QPushButton("Save Results to CSV")
        self.btn_save.setEnabled(False)
        self.btn_save.clicked.connect(self.save_results)
        left_panel.addWidget(self.btn_save)
        
        # Treeview (Table)
        self.table_results = QTableWidget()
        self.table_results.setColumnCount(7)
        self.table_results.setHorizontalHeaderLabels([
            "Comparison", "Feature", "Sample Sizes", "KS Statistic", "P-Value", "Significance", "KL Divergence"
        ])
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

    def load_files(self):
        file_paths, _ = QFileDialog.getOpenFileNames(self, "Open CSV Files", "", "CSV Files (*.csv)")
        if not file_paths:
            return
            
        self.loaded_datasets = {}
        first_cols = None
        
        for path in file_paths:
            filename = os.path.basename(path)
            try:
                df = pd.read_csv(path)
                cleaned_data = {}
                
                for col in df.columns:
                    numeric_series = pd.to_numeric(df[col], errors='coerce')
                    if numeric_series.notna().sum() > 10:
                        cleaned_data[col] = numeric_series
                        
                numeric_df = pd.DataFrame(cleaned_data)
                if not numeric_df.empty:
                    self.loaded_datasets[filename] = numeric_df
                    if first_cols is None:
                        first_cols = set(numeric_df.columns)
                    else:
                        first_cols = first_cols.intersection(set(numeric_df.columns))
            except Exception as e:
                print(f"Failed to clean {filename}: {e}")
                
        if not self.loaded_datasets:
            QMessageBox.critical(self, "Error", "No valid numerical data found.")
            return
            
        self.common_columns = sorted(list(first_cols)) if first_cols else []
        self.combo_columns.clear()
        self.combo_columns.addItems(self.common_columns)
        
        self.combo_ref.clear()
        self.combo_ref.addItems(list(self.loaded_datasets.keys()))
        
        self.lbl_status.setText(f"Loaded {len(self.loaded_datasets)} file(s). Shared columns: {len(self.common_columns)}")

    def run_comparison(self):
        if not self.loaded_datasets:
            QMessageBox.warning(self, "Warning", "Please load CSV files first.")
            return
            
        ref_name = self.combo_ref.currentText()
        target_col = self.combo_columns.currentText()
        cores = self.spin_processors.value()
        
        if not ref_name or not target_col:
            QMessageBox.warning(self, "Warning", "Select reference and target column.")
            return
            
        ref_df = self.loaded_datasets[ref_name]
        comparison_dfs = {name: df for name, df in self.loaded_datasets.items() if name != ref_name}
        
        # Parallel stats computation
        results_df = compare_datasets_parallel(ref_df, comparison_dfs, target_col, cores)
        
        # Populate Table
        self.table_results.setRowCount(len(results_df))
        for r_idx, row in results_df.iterrows():
            self.table_results.setItem(r_idx, 0, QTableWidgetItem(str(row['Comparison'])))
            self.table_results.setItem(r_idx, 1, QTableWidgetItem(str(row['Feature'])))
            self.table_results.setItem(r_idx, 2, QTableWidgetItem(str(row['Sample Sizes'])))
            self.table_results.setItem(r_idx, 3, QTableWidgetItem(f"{row['KS Statistic']:.4f}"))
            self.table_results.setItem(r_idx, 4, QTableWidgetItem(f"{row['P-Value']:.2e}"))
            self.table_results.setItem(r_idx, 5, QTableWidgetItem(str(row['Significance'])))
            self.table_results.setItem(r_idx, 6, QTableWidgetItem(f"{row['KL Divergence']:.4f}"))
            
        self.btn_save.setEnabled(True)
        self.results_df = results_df
        
        # Update Plot
        if FigureCanvas is not None:
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            
            ref_vals = ref_df[target_col].dropna().values
            ax.hist(ref_vals, bins=30, alpha=0.5, label=f'Ref: {ref_name}', density=True)
            
            for name, df in comparison_dfs.items():
                if target_col in df.columns:
                    ax.hist(df[target_col].dropna().values, bins=30, alpha=0.5, label=f'Comp: {name}', density=True)
                    
            ax.set_title(f'Overlay Distribution of {target_col}')
            ax.set_xlabel(target_col)
            ax.set_ylabel('Density')
            ax.legend()
            self.canvas.draw()

    def save_results(self):
        if not hasattr(self, 'results_df') or self.results_df.empty:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save Comparison Table", "", "CSV Files (*.csv)")
        if path:
            self.results_df.to_csv(path, index=False)
            QMessageBox.information(self, "Saved", f"Results exported to:\n{path}")


# =========================================================================
# PLUGIN INTERFACE (Integrates into PyChem-Pro if imported)
# =========================================================================

class StatComparisonPlugin(BasePlugin):
    def __init__(self):
        super().__init__(PluginInfo(
            name="Chemical Space Statistical Comparator",
            version="1.1.0",
            description="Compares chemical datasets using KS-Test & KL Divergence in parallel.",
            author="DeepMind Antigravity",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        ))
        self.widget = None

    def create_widget(self):
        if self.widget is None:
            self.widget = StatComparisonWidget()
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
        window.setWindowTitle("Statistical Space Comparator")
        window.setCentralWidget(StatComparisonWidget())
        window.resize(1100, 700)
        window.show()
        sys.exit(app.exec())
    except Exception as e:
        print(f"Failed to start GUI: {e}")