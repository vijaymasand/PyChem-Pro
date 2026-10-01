"""
PyChem-Pro Fingerprint Acquisition & QSAR Prep Plugin
Designed to calculate molecular fingerprints for QSAR model building and chemical space analysis.

Supports:
  - Morgan (ECFP4/ECFP6-like, customizable radius=1..4 and bit length=1024..4096)
  - Topological (path-based, Daylight-like, bit length=1024..4096)
  - MACCS Keys (166 structural keys)
  - Combined (Morgan + Topological, doubled feature space)

QSAR Features & Optimizations:
  - Near-Zero Variance Bit Pruning (filters constant/uninformative bits)
  - Min/Max Bit Frequency Filtering (filters rare or ubiquitous features)
  - Redundant Bit Deduplication
  - Automatic Preservation of Activity Target Columns & Train/Test Split Columns (1=Train, 2=Test)
  - Fast Pairwise Tanimoto Similarity Matrix & Dataset Diversity Index Calculation
  - Full Integration into PyChem-Pro GUI Plugin Architecture (Dockable/Popout Widget + Active Molecule Inspector)

Usage (CLI):
  python calculate_fingerprints.py -i input.csv -o output_fp.csv -t morgan --variance-thresh 0.001 --tanimoto-matrix

Usage (GUI / Plugin):
  Discovered automatically by PyChem-Pro PluginManager or executed directly via:
  python calculate_fingerprints.py
"""

import os
import sys
import time
import argparse
import multiprocessing
from typing import Dict, Any, List, Optional, Tuple

import numpy as np
import pandas as pd

# Ensure the project root and plugins directory are in python path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.insert(0, project_root)
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from src.features.smiles_parser.services.parser import parse_smiles
from src.features.cheminformatics.topology.fingerprints import (
    get_morgan_fingerprint,
    get_topological_fingerprint,
    get_maccs_keys,
    tanimoto_similarity
)

# PyChem-Pro Plugin Framework
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType

# Qt Imports from PyChem-Pro compatibility layer
from src.shared.qt_compat import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QFileDialog, QMessageBox, QComboBox, QPushButton, QLabel,
    QLineEdit, QProgressBar, QThread, Signal, QFont, Qt, QTextEdit,
    QGroupBox, QFormLayout, QSpinBox, QDoubleSpinBox, QCheckBox,
    QTableWidget, QTableWidgetItem, QTabWidget, QDialog
)

# Fingerprint type mappings
FP_TYPES = {
    'morgan': 'Morgan (ECFP-like)',
    'topological': 'Topological (Path-based)',
    'maccs': 'MACCS Keys (166 bits)',
    'combined': 'Combined (Morgan + Topological)'
}


# =============================================================================
# Pure NumPy / Algorithmic Helpers for Fingerprints & QSAR Feature Processing
# =============================================================================

def _compute_single_smiles(args):
    """
    Worker function to process a single SMILES string in parallel pool.
    args tuple: (idx, smiles, name, fp_type, radius, n_bits)
    """
    idx, smiles, name, fp_type, radius, n_bits = args
    try:
        mol = parse_smiles(smiles)
        if mol is None:
            return idx, None, "Parsing returned None"
        
        if fp_type == 'morgan':
            vec = get_morgan_fingerprint(mol, radius=radius, n_bits=n_bits)
        elif fp_type == 'topological':
            vec = get_topological_fingerprint(mol, min_path=1, max_path=7, n_bits=n_bits)
        elif fp_type == 'maccs':
            vec = get_maccs_keys(mol)
            if len(vec) == 167:
                vec = vec[1:]  # Drop dummy 0-index if present for standard 166 bits
        elif fp_type == 'combined':
            half_bits = max(512, n_bits // 2)
            morgan = get_morgan_fingerprint(mol, radius=radius, n_bits=half_bits)
            topo = get_topological_fingerprint(mol, min_path=1, max_path=7, n_bits=half_bits)
            vec = morgan + topo
        else:
            vec = get_morgan_fingerprint(mol, radius=radius, n_bits=n_bits)
            
        return idx, vec, None
    except Exception as e:
        return idx, None, str(e)


def prune_fingerprint_matrix(
    raw_matrix: np.ndarray,
    variance_thresh: float = 0.0,
    min_freq: float = 0.0,
    max_freq: float = 1.0,
    remove_duplicates: bool = False
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Prunes binary fingerprint matrix for QSAR model optimization.
    
    Returns:
        pruned_matrix (N x B_retained), retained_bit_indices, stats_dict
    """
    X = np.array(raw_matrix, dtype=np.uint8)
    n_samples, orig_nbits = X.shape

    if n_samples == 0 or orig_nbits == 0:
        return X, np.arange(orig_nbits), {
            'original_bits': orig_nbits,
            'retained_bits': orig_nbits,
            'removed_bits': 0,
            'avg_active_bits': 0.0,
            'density_pct': 0.0
        }

    # 1. Variance Filter
    variances = np.var(X, axis=0) if n_samples > 1 else np.ones(orig_nbits)
    var_mask = variances >= variance_thresh

    # 2. Frequency Filter
    frequencies = np.mean(X, axis=0)
    freq_mask = (frequencies >= min_freq) & (frequencies <= max_freq)

    valid_mask = var_mask & freq_mask
    valid_indices = np.where(valid_mask)[0]

    # Fallback safety: Keep at least non-zero variance bits if filters are too strict
    if len(valid_indices) == 0:
        valid_indices = np.where(variances > 0)[0]
        if len(valid_indices) == 0:
            valid_indices = np.arange(orig_nbits)

    X_sub = X[:, valid_indices]

    # 3. Deduplication of identical bit columns
    if remove_duplicates and X_sub.shape[1] > 1:
        # Transpose bit columns and find unique patterns across compounds
        _, unique_col_idx = np.unique(X_sub, axis=1, return_index=True)
        unique_col_idx.sort()
        valid_indices = valid_indices[unique_col_idx]
        X_sub = X[:, valid_indices]

    stats = {
        'original_bits': orig_nbits,
        'retained_bits': len(valid_indices),
        'removed_bits': orig_nbits - len(valid_indices),
        'avg_active_bits': float(np.mean(np.sum(X_sub, axis=1))) if len(valid_indices) > 0 else 0.0,
        'density_pct': float(np.mean(X_sub) * 100) if len(valid_indices) > 0 else 0.0
    }

    return X_sub, valid_indices, stats


def compute_tanimoto_matrix(X: np.ndarray) -> np.ndarray:
    """
    Computes NxN Tanimoto similarity matrix using vector dot-products.
    """
    X_f = X.astype(np.float32)
    intersection = np.dot(X_f, X_f.T)
    bit_sums = np.sum(X_f, axis=1)
    union = bit_sums[:, None] + bit_sums[None, :] - intersection
    with np.errstate(divide='ignore', invalid='ignore'):
        sim = np.where(union > 0, intersection / union, 1.0)
    return sim


# =============================================================================
# QThread Background Worker for Parallel Processing
# =============================================================================

class FingerprintWorker(QThread):
    progress = Signal(int, int)
    finished = Signal(pd.DataFrame, object, int, int, dict, list)
    error = Signal(str)

    def __init__(self, smiles_list: List[str], ids_list: List[str], id_col_name: str,
                 fp_type: str, radius: int, n_bits: int, n_jobs: int,
                 variance_thresh: float = 0.0, min_freq: float = 0.0, max_freq: float = 1.0,
                 remove_duplicates: bool = False, compute_tanimoto: bool = False,
                 extra_df: Optional[pd.DataFrame] = None, target_cols: Optional[List[str]] = None,
                 split_col: Optional[str] = None):
        super().__init__()
        self.smiles_list = smiles_list
        self.ids_list = ids_list
        self.id_col_name = id_col_name
        self.fp_type = fp_type
        self.radius = radius
        self.n_bits = n_bits
        self.n_jobs = n_jobs
        self.variance_thresh = variance_thresh
        self.min_freq = min_freq
        self.max_freq = max_freq
        self.remove_duplicates = remove_duplicates
        self.compute_tanimoto = compute_tanimoto
        self.extra_df = extra_df if extra_df is not None else pd.DataFrame()
        self.target_cols = target_cols if target_cols else []
        self.split_col = split_col

    def run(self):
        try:
            worker_args = [
                (idx, smiles, name, self.fp_type, self.radius, self.n_bits)
                for idx, (smiles, name) in enumerate(zip(self.smiles_list, self.ids_list))
            ]
            
            results = [None] * len(worker_args)
            errors = []
            processed_count = 0
            total = len(worker_args)
            chunk_size = max(1, total // (self.n_jobs * 20))
            
            with multiprocessing.Pool(processes=self.n_jobs) as pool:
                for idx, vec, err in pool.imap_unordered(_compute_single_smiles, worker_args, chunksize=chunk_size):
                    results[idx] = vec
                    if err:
                        errors.append(f"Row {idx+1} | ID: {self.ids_list[idx]} | Error: {err}")
                    processed_count += 1
                    self.progress.emit(processed_count, total)

            success_rows = []
            success_ids = []
            success_smiles = []
            success_indices = []
            failed_count = 0

            for idx, vec in enumerate(results):
                if vec is not None:
                    success_rows.append(vec)
                    success_ids.append(self.ids_list[idx])
                    success_smiles.append(self.smiles_list[idx])
                    success_indices.append(idx)
                else:
                    failed_count += 1

            if len(success_rows) == 0:
                self.error.emit("No molecules were successfully processed.")
                return

            # Prune bit matrix for QSAR modeling
            raw_matrix = np.array(success_rows, dtype=np.uint8)
            pruned_matrix, retained_bits, stats = prune_fingerprint_matrix(
                raw_matrix,
                variance_thresh=self.variance_thresh,
                min_freq=self.min_freq,
                max_freq=self.max_freq,
                remove_duplicates=self.remove_duplicates
            )

            col_names = [f"FP_bit{b}" for b in retained_bits]
            output_df = pd.DataFrame(pruned_matrix, columns=col_names)

            # Insert SMILES and ID
            output_df.insert(0, 'SMILES', success_smiles)
            output_df.insert(0, self.id_col_name, success_ids)

            # Preserve QSAR Split and Activity Target Columns
            sub_extra = self.extra_df.iloc[success_indices] if len(self.extra_df) == len(self.smiles_list) else None
            if sub_extra is not None:
                if self.split_col and self.split_col in sub_extra.columns:
                    output_df.insert(2, self.split_col, sub_extra[self.split_col].values)

                for tcol in self.target_cols:
                    if tcol in sub_extra.columns and tcol not in output_df.columns:
                        insert_pos = 3 if (self.split_col and self.split_col in output_df.columns) else 2
                        output_df.insert(insert_pos, tcol, sub_extra[tcol].values)

            # Calculate Tanimoto similarity matrix if enabled
            tanimoto_df = None
            if self.compute_tanimoto:
                tan_matrix = compute_tanimoto_matrix(pruned_matrix)
                tanimoto_df = pd.DataFrame(tan_matrix, index=success_ids, columns=success_ids)
                tanimoto_df.insert(0, self.id_col_name, success_ids)

            self.finished.emit(output_df, tanimoto_df, len(success_rows), failed_count, stats, errors)
        except Exception as e:
            self.error.emit(str(e))


# =============================================================================
# Plugin Widget Interface
# =============================================================================

class CalculateFingerprintsWidget(PluginWidget):
    """
    Main Plugin Widget for Fingerprint Acquisition & QSAR Prep.
    Includes docking popout controls, 3 main feature tabs, and real-time inspector.
    """

    def __init__(self, plugin: Optional['CalculateFingerprintsPlugin'] = None):
        super().__init__(plugin)
        self.input_df = None
        self.worker = None
        self.latest_output_df = None
        self.latest_tanimoto_df = None
        self.setup_ui()

    def setup_ui(self):
        self.widget = QWidget()
        _base_layout = QVBoxLayout(self.widget)
        _base_layout.setContentsMargins(5, 5, 5, 5)

        # Top Control Bar (Popout, Fullscreen, Dock Restore)
        _top_bar = QHBoxLayout()
        _top_bar.addStretch()
        self.btn_maximize = QPushButton("Maximize Window")
        self.btn_maximize.clicked.connect(self.toggle_maximize)
        _top_bar.addWidget(self.btn_maximize)

        self.btn_fullscreen = QPushButton("Full Screen")
        self.btn_fullscreen.clicked.connect(lambda: self._open_popout(fullscreen=True))
        _top_bar.addWidget(self.btn_fullscreen)

        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(self.widget.close)
        _top_bar.addWidget(self.btn_close)
        _base_layout.addLayout(_top_bar)

        # Main Container
        self.main_content_widget = QWidget()
        _base_layout.addWidget(self.main_content_widget)
        main_layout = QVBoxLayout(self.main_content_widget)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)

        # Main Styling
        self.widget.setStyleSheet("""
            QWidget {
                background-color: #1e1e24;
                color: #ffffff;
                font-family: 'Segoe UI', Arial, sans-serif;
            }
            QGroupBox {
                background-color: #2a2a35;
                border: 1px solid #3e3e4f;
                border-radius: 8px;
                margin-top: 10px;
                padding-top: 12px;
                font-weight: bold;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 5px;
                color: #00bcd4;
            }
            QLabel {
                font-size: 13px;
                color: #e0e0e0;
            }
            QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {
                background-color: #16161d;
                border: 1px solid #3e3e4f;
                border-radius: 4px;
                padding: 6px;
                color: #ffffff;
            }
            QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {
                border: 1px solid #00bcd4;
            }
            QPushButton {
                background-color: #3f51b5;
                color: white;
                border: none;
                border-radius: 4px;
                padding: 8px 16px;
                font-weight: bold;
                font-size: 13px;
            }
            QPushButton:hover {
                background-color: #4f51c5;
            }
            QPushButton#btn_run {
                background-color: #008b94;
                font-size: 14px;
                padding: 10px;
            }
            QPushButton#btn_run:hover {
                background-color: #00a4b4;
            }
            QProgressBar {
                border: 1px solid #3e3e4f;
                border-radius: 4px;
                text-align: center;
                background-color: #16161d;
            }
            QProgressBar::chunk {
                background-color: #00bcd4;
            }
            QTextEdit {
                background-color: #16161d;
                border: 1px solid #3e3e4f;
                border-radius: 6px;
                font-family: 'Consolas', 'Courier New', monospace;
                font-size: 12px;
                color: #a0ffd0;
            }
            QTabWidget::pane {
                border: 1px solid #3e3e4f;
                background-color: #1e1e24;
                border-radius: 6px;
            }
            QTabBar::tab {
                background-color: #2a2a35;
                color: #cccccc;
                padding: 8px 16px;
                border-top-left-radius: 4px;
                border-top-right-radius: 4px;
                margin-right: 2px;
            }
            QTabBar::tab:selected {
                background-color: #00bcd4;
                color: #ffffff;
                font-weight: bold;
            }
        """)

        # Header Title
        header = QLabel("Fingerprint Acquisition & QSAR Prep Engine")
        header.setFont(QFont("Segoe UI", 15))
        header.setStyleSheet("color: #00bcd4; font-weight: bold; margin-bottom: 5px;")
        main_layout.addWidget(header)

        # Tab Widget
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)

        # Build Tabs
        self.tab_batch = QWidget()
        self.tab_inspector = QWidget()
        self.tab_matrix = QWidget()

        self.tabs.addTab(self.tab_batch, "1. Batch Dataset & QSAR Prep")
        self.tabs.addTab(self.tab_inspector, "2. Active Molecule Inspector")
        self.tabs.addTab(self.tab_matrix, "3. Tanimoto & Dataset Diversity")

        self._build_batch_tab()
        self._build_inspector_tab()
        self._build_matrix_tab()

    # -------------------------------------------------------------------------
    # Tab 1: Batch Dataset Processing & QSAR Pruning
    # -------------------------------------------------------------------------
    def _build_batch_tab(self):
        layout = QVBoxLayout(self.tab_batch)

        # 1. Dataset Load Box
        grp_file = QGroupBox("1. Dataset Import")
        lay_file = QHBoxLayout(grp_file)
        self.txt_file = QLineEdit()
        self.txt_file.setReadOnly(True)
        self.txt_file.setPlaceholderText("Select CSV, TSV, or Excel dataset file...")
        lay_file.addWidget(self.txt_file)

        btn_browse = QPushButton("Browse...")
        btn_browse.clicked.connect(self.browse_file)
        lay_file.addWidget(btn_browse)
        layout.addWidget(grp_file)

        # 2. Configuration & QSAR Preservation Box
        grp_config = QGroupBox("2. Feature & QSAR Column Setup")
        lay_config = QFormLayout(grp_config)
        lay_config.setContentsMargins(15, 12, 15, 12)

        self.cmb_smiles = QComboBox()
        self.cmb_smiles.setEnabled(False)
        lay_config.addRow("SMILES Column:", self.cmb_smiles)

        self.cmb_id = QComboBox()
        self.cmb_id.setEnabled(False)
        lay_config.addRow("ID / Name Column:", self.cmb_id)

        self.cmb_target = QComboBox()
        self.cmb_target.setEnabled(False)
        lay_config.addRow("Activity Target Column (Optional):", self.cmb_target)

        self.cmb_split = QComboBox()
        self.cmb_split.setEnabled(False)
        lay_config.addRow("Train/Test Split Column (1=Train, 2=Test):", self.cmb_split)

        self.cmb_fp_type = QComboBox()
        for key, val in FP_TYPES.items():
            self.cmb_fp_type.addItem(val, key)
        self.cmb_fp_type.currentIndexChanged.connect(self._on_fp_type_changed)
        lay_config.addRow("Fingerprint Type:", self.cmb_fp_type)

        # Parameters row
        param_layout = QHBoxLayout()
        self.spn_radius = QSpinBox()
        self.spn_radius.setRange(1, 4)
        self.spn_radius.setValue(2)
        param_layout.addWidget(QLabel("Morgan Radius:"))
        param_layout.addWidget(self.spn_radius)

        self.cmb_nbits = QComboBox()
        for b in [1024, 2048, 4096]:
            self.cmb_nbits.addItem(f"{b} bits", b)
        self.cmb_nbits.setCurrentIndex(1)
        param_layout.addWidget(QLabel("Bit Vector Size:"))
        param_layout.addWidget(self.cmb_nbits)
        lay_config.addRow("Algorithm Options:", param_layout)

        layout.addWidget(grp_config)

        # 3. QSAR Bit Pruning & Similarity Options Box
        grp_qsar = QGroupBox("3. QSAR Feature Pruning & Matrix Options")
        lay_qsar = QFormLayout(grp_qsar)

        self.spn_var = QDoubleSpinBox()
        self.spn_var.setRange(0.0, 0.25)
        self.spn_var.setSingleStep(0.001)
        self.spn_var.setDecimals(4)
        self.spn_var.setValue(0.0000)
        self.spn_var.setToolTip("Removes constant/near-zero variance bits with Var < threshold")
        lay_qsar.addRow("Variance Threshold (Var < T removed):", self.spn_var)

        freq_layout = QHBoxLayout()
        self.spn_min_freq = QDoubleSpinBox()
        self.spn_min_freq.setRange(0.0, 50.0)
        self.spn_min_freq.setValue(0.0)
        self.spn_min_freq.setSuffix("%")

        self.spn_max_freq = QDoubleSpinBox()
        self.spn_max_freq.setRange(50.0, 100.0)
        self.spn_max_freq.setValue(100.0)
        self.spn_max_freq.setSuffix("%")

        freq_layout.addWidget(QLabel("Min:"))
        freq_layout.addWidget(self.spn_min_freq)
        freq_layout.addWidget(QLabel("Max:"))
        freq_layout.addWidget(self.spn_max_freq)
        lay_qsar.addRow("Bit Occurrence Frequency %:", freq_layout)

        self.chk_dedup = QCheckBox("Remove Identical Duplicate Bit Columns")
        lay_qsar.addRow("Deduplication:", self.chk_dedup)

        self.chk_tanimoto = QCheckBox("Calculate & Export Tanimoto Similarity Matrix (N x N)")
        lay_qsar.addRow("Distance Kernel:", self.chk_tanimoto)

        layout.addWidget(grp_qsar)

        # 4. Action & Progress
        self.btn_run = QPushButton("Run Acquisition && Save QSAR Dataset...")
        self.btn_run.setObjectName("btn_run")
        self.btn_run.setEnabled(False)
        self.btn_run.clicked.connect(self.start_calculation)
        layout.addWidget(self.btn_run)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)

        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setPlaceholderText("Calculation status and feature stats will be reported here...")
        layout.addWidget(self.log_output)

        self.log("System Ready. Please select an input compound dataset.")

    def _on_fp_type_changed(self, index):
        fp_type = self.cmb_fp_type.currentData()
        self.spn_radius.setEnabled(fp_type in ['morgan', 'combined'])

    # -------------------------------------------------------------------------
    # Tab 2: Active Molecule Real-time Inspector
    # -------------------------------------------------------------------------
    def _build_inspector_tab(self):
        layout = QVBoxLayout(self.tab_inspector)

        grp_mol = QGroupBox("Single Molecule Fingerprint Real-Time Inspector")
        lay_mol = QFormLayout(grp_mol)

        self.txt_inspect_smiles = QLineEdit()
        self.txt_inspect_smiles.setPlaceholderText("Enter SMILES string (or click load active molecule below)...")
        lay_mol.addRow("SMILES String:", self.txt_inspect_smiles)

        btn_load_active = QPushButton("Fetch Active Molecule from PyChem-Pro Workspace")
        btn_load_active.clicked.connect(self.fetch_active_molecule)
        lay_mol.addRow("", btn_load_active)

        btn_inspect = QPushButton("Calculate & Inspect Fingerprints")
        btn_inspect.clicked.connect(self.inspect_single_smiles)
        lay_mol.addRow("", btn_inspect)

        layout.addWidget(grp_mol)

        # Summary box
        grp_stats = QGroupBox("Fingerprint Vector Summary")
        lay_stats = QVBoxLayout(grp_stats)
        self.txt_inspect_result = QTextEdit()
        self.txt_inspect_result.setReadOnly(True)
        lay_stats.addWidget(self.txt_inspect_result)
        layout.addWidget(grp_stats)

        # Similarity comparator box
        grp_comp = QGroupBox("Tanimoto Similarity Lead Comparison")
        lay_comp = QFormLayout(grp_comp)

        self.txt_ref_smiles = QLineEdit()
        self.txt_ref_smiles.setPlaceholderText("Enter reference/lead molecule SMILES...")
        lay_comp.addRow("Reference SMILES:", self.txt_ref_smiles)

        btn_compare = QPushButton("Calculate Pairwise Tanimoto Similarity")
        btn_compare.clicked.connect(self.compare_smiles_tanimoto)
        lay_comp.addRow("", btn_compare)

        self.lbl_similarity_res = QLabel("Similarity Score: N/A")
        self.lbl_similarity_res.setFont(QFont("Segoe UI", 12))
        self.lbl_similarity_res.setStyleSheet("color: #00bcd4; font-weight: bold;")
        lay_comp.addRow("", self.lbl_similarity_res)

        layout.addWidget(grp_comp)

    # -------------------------------------------------------------------------
    # Tab 3: Tanimoto & Dataset Diversity Matrix Tool
    # -------------------------------------------------------------------------
    def _build_matrix_tab(self):
        layout = QVBoxLayout(self.tab_matrix)

        grp_mat = QGroupBox("Dataset Tanimoto Matrix & Chemical Diversity Metric")
        lay_mat = QVBoxLayout(grp_mat)

        self.lbl_diversity_stats = QLabel("Load or calculate dataset fingerprints in Tab 1 to inspect matrix stats.")
        self.lbl_diversity_stats.setWordWrap(True)
        lay_mat.addWidget(self.lbl_diversity_stats)

        self.tbl_tanimoto = QTableWidget()
        self.tbl_tanimoto.setColumnCount(0)
        self.tbl_tanimoto.setRowCount(0)
        lay_mat.addWidget(self.tbl_tanimoto)

        btn_export_matrix = QPushButton("Export Tanimoto Matrix CSV...")
        btn_export_matrix.clicked.connect(self.export_tanimoto_matrix)
        lay_mat.addWidget(btn_export_matrix)

        layout.addWidget(grp_mat)

    # -------------------------------------------------------------------------
    # Dataset Loading & Callbacks
    # -------------------------------------------------------------------------
    def log(self, msg: str):
        self.log_output.append(f"[{time.strftime('%H:%M:%S')}] {msg}")

    def browse_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self.widget, "Open SMILES Dataset", "", "Data Files (*.csv *.tsv *.xlsx *.xls)"
        )
        if not file_path:
            return

        self.txt_file.setText(file_path)
        self.log(f"Loading file: {file_path}")

        try:
            if file_path.endswith(('.xlsx', '.xls')):
                self.input_df = pd.read_excel(file_path, nrows=5)
            else:
                with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                    first_line = f.readline()
                    sep = '\t' if '\t' in first_line else ','
                self.input_df = pd.read_csv(file_path, sep=sep, nrows=5)

            cols = list(self.input_df.columns)
            self.cmb_smiles.clear()
            self.cmb_id.clear()
            self.cmb_target.clear()
            self.cmb_split.clear()

            self.cmb_smiles.addItems(cols)
            self.cmb_id.addItems(cols)
            
            self.cmb_target.addItem("(None)", "")
            self.cmb_target.addItems(cols)

            self.cmb_split.addItem("(None)", "")
            self.cmb_split.addItems(cols)

            cols_lower = [c.lower() for c in cols]

            # Auto-detect SMILES column
            if 'smiles' in cols_lower:
                self.cmb_smiles.setCurrentIndex(cols_lower.index('smiles'))

            # Auto-detect ID column
            id_targets = ['immpat id', 'id', 'name', 'pubchem id', 'mol_id', 'compound_id']
            id_idx = [cols_lower.index(x) for x in id_targets if x in cols_lower]
            if id_idx:
                self.cmb_id.setCurrentIndex(id_idx[0])

            # Auto-detect Target Activity column
            target_matches = [cols_lower.index(x) for x in ['pic50', 'ic50', 'activity', 'logp', 'exp_val', 'value', 'label'] if x in cols_lower]
            if target_matches:
                self.cmb_target.setCurrentIndex(target_matches[0] + 1)

            # Auto-detect Split column
            split_matches = [cols_lower.index(x) for x in ['split', 'set', 'subset'] if x in cols_lower]
            if split_matches:
                self.cmb_split.setCurrentIndex(split_matches[0] + 1)

            self.cmb_smiles.setEnabled(True)
            self.cmb_id.setEnabled(True)
            self.cmb_target.setEnabled(True)
            self.cmb_split.setEnabled(True)
            self.btn_run.setEnabled(True)
            self.log(f"Header loaded successfully ({len(cols)} columns detected).")

        except Exception as e:
            QMessageBox.critical(self.widget, "Error Loading File", f"Could not read dataset headers:\n{str(e)}")
            self.log(f"Error: {str(e)}")

    def start_calculation(self):
        input_path = self.txt_file.text()
        smiles_col = self.cmb_smiles.currentText()
        id_col = self.cmb_id.currentText()
        target_col = self.cmb_target.currentText() if self.cmb_target.currentIndex() > 0 else None
        split_col = self.cmb_split.currentText() if self.cmb_split.currentIndex() > 0 else None

        fp_type = self.cmb_fp_type.currentData()
        radius = self.spn_radius.value()
        n_bits = self.cmb_nbits.currentData()

        var_thresh = self.spn_var.value()
        min_freq = self.spn_min_freq.value() / 100.0
        max_freq = self.spn_max_freq.value() / 100.0
        remove_dedup = self.chk_dedup.isChecked()
        compute_tanimoto = self.chk_tanimoto.isChecked()

        output_path, _ = QFileDialog.getSaveFileName(
            self.widget, "Save QSAR Fingerprint Dataset", "fingerprints_qsar.csv", "CSV Files (*.csv);;TSV Files (*.tsv)"
        )
        if not output_path:
            return

        self.log("Loading dataset records into memory...")
        try:
            if input_path.endswith(('.xlsx', '.xls')):
                full_df = pd.read_excel(input_path)
            else:
                with open(input_path, 'r', encoding='utf-8', errors='ignore') as f:
                    first_line = f.readline()
                    sep = '\t' if '\t' in first_line else ','
                full_df = pd.read_csv(input_path, sep=sep)
        except Exception as e:
            QMessageBox.critical(self.widget, "Dataset Load Error", f"Failed to load dataset:\n{str(e)}")
            return

        initial_len = len(full_df)
        full_df = full_df.dropna(subset=[smiles_col])
        full_df = full_df[full_df[smiles_col].astype(str).str.strip() != '']
        valid_len = len(full_df)

        if initial_len != valid_len:
            self.log(f"Filtered out {initial_len - valid_len} empty/invalid SMILES rows.")

        if valid_len == 0:
            QMessageBox.warning(self.widget, "Warning", "No valid SMILES entries found.")
            return

        smiles_list = full_df[smiles_col].astype(str).tolist()
        ids_list = full_df[id_col].astype(str).tolist()

        target_cols = [target_col] if target_col else []

        n_jobs = max(1, multiprocessing.cpu_count() // 2)

        self.log(f"Computing {FP_TYPES[fp_type]} ({n_bits} bits) for {valid_len} compounds on {n_jobs} cores...")
        self.btn_run.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)

        self.worker = FingerprintWorker(
            smiles_list=smiles_list,
            ids_list=ids_list,
            id_col_name=id_col,
            fp_type=fp_type,
            radius=radius,
            n_bits=n_bits,
            n_jobs=n_jobs,
            variance_thresh=var_thresh,
            min_freq=min_freq,
            max_freq=max_freq,
            remove_duplicates=remove_dedup,
            compute_tanimoto=compute_tanimoto,
            extra_df=full_df,
            target_cols=target_cols,
            split_col=split_col
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.finished.connect(
            lambda out_df, tan_df, s_cnt, f_cnt, st, errs: self.on_finished(
                out_df, tan_df, s_cnt, f_cnt, st, errs, output_path
            )
        )
        self.worker.error.connect(self.on_error)
        self.worker.start()

    def on_progress(self, current: int, total: int):
        pct = int(100 * current / total)
        self.progress_bar.setValue(pct)
        if current % 200 == 0 or current == total:
            self.log(f"Processed {current}/{total} molecules ({pct}%)...")

    def on_finished(self, output_df: pd.DataFrame, tanimoto_df: Optional[pd.DataFrame],
                    success_count: int, fail_count: int, stats: dict, errors: list, output_path: str):
        self.progress_bar.setVisible(False)
        self.btn_run.setEnabled(True)

        self.latest_output_df = output_df
        self.latest_tanimoto_df = tanimoto_df

        self.log("Saving processed QSAR dataset to disk...")
        try:
            sep = '\t' if output_path.endswith('.tsv') else ','
            output_df.to_csv(output_path, index=False, sep=sep)
            self.log(f"✓ Saved QSAR dataset: {output_path}")

            # Save Tanimoto matrix if present
            if tanimoto_df is not None:
                base_name, ext = os.path.splitext(output_path)
                tan_path = f"{base_name}_tanimoto_matrix{ext}"
                tanimoto_df.to_csv(tan_path, index=False, sep=sep)
                self.log(f"✓ Saved Tanimoto Matrix: {tan_path}")
                self._update_matrix_tab_display(tanimoto_df)

            summary = (
                f"Fingerprint Acquisition & QSAR Prep Complete!\n\n"
                f"Processed Molecules: {success_count} success, {fail_count} failed\n"
                f"Original Bit Features: {stats['original_bits']}\n"
                f"Retained QSAR Bits:    {stats['retained_bits']} (Removed {stats['removed_bits']} bits)\n"
                f"Avg Active Bits/Mol:   {stats['avg_active_bits']:.2f} ({stats['density_pct']:.2f}% bit density)\n"
                f"Dataset Shape:         {output_df.shape}\n\n"
                f"Output Path: {output_path}"
            )
            QMessageBox.information(self.widget, "QSAR Fingerprint Success", summary)
            self.log(summary.replace('\n\n', ' | ').replace('\n', ' | '))

            if errors:
                self.log(f"--- Encounted {len(errors)} warnings/errors ---")
                for err in errors[:5]:
                    self.log(err)

        except Exception as e:
            QMessageBox.critical(self.widget, "Export Error", f"Could not write output files:\n{str(e)}")
            self.log(f"Save error: {str(e)}")

    def on_error(self, err_msg: str):
        self.progress_bar.setVisible(False)
        self.btn_run.setEnabled(True)
        QMessageBox.critical(self.widget, "Worker Error", f"Error during execution:\n{err_msg}")
        self.log(f"Fatal error: {err_msg}")

    # -------------------------------------------------------------------------
    # Inspector & Matrix Helpers
    # -------------------------------------------------------------------------
    def fetch_active_molecule(self):
        if self.plugin and hasattr(self.plugin, 'current_molecule') and self.plugin.current_molecule:
            mol = self.plugin.current_molecule
            smiles = getattr(mol, 'smiles', '')
            if smiles:
                self.txt_inspect_smiles.setText(smiles)
                self.inspect_single_smiles()
            else:
                QMessageBox.information(self.widget, "Molecule Active", "Active molecule found in workspace without SMILES string.")
        else:
            QMessageBox.information(self.widget, "No Active Molecule", "No active molecule currently loaded in PyChem-Pro window.")

    def inspect_single_smiles(self):
        smiles = self.txt_inspect_smiles.text().strip()
        if not smiles:
            QMessageBox.warning(self.widget, "Missing SMILES", "Please enter a valid SMILES string.")
            return

        try:
            mol = parse_smiles(smiles)
            if mol is None:
                QMessageBox.warning(self.widget, "Parse Error", "Failed to parse SMILES string.")
                return

            morg = get_morgan_fingerprint(mol, radius=2, n_bits=2048)
            topo = get_topological_fingerprint(mol, min_path=1, max_path=7, n_bits=2048)
            maccs = get_maccs_keys(mol)

            active_morg = [i for i, b in enumerate(morg) if b == 1]
            active_topo = [i for i, b in enumerate(topo) if b == 1]
            active_maccs = [i for i, b in enumerate(maccs) if b == 1]

            res = [
                f"Compound Analysis: {smiles}",
                f"Atoms: {len(mol.atoms)} | Heavy Atoms: {sum(1 for a in mol.atoms if a.atomic_number > 1)}",
                "-" * 60,
                f"Morgan (ECFP4 2048):       {len(active_morg)} active bits ({len(active_morg)/2048*100:.2f}% density)",
                f"  Indices (first 15):     {active_morg[:15]}...",
                f"Topological (Path 2048):   {len(active_topo)} active bits ({len(active_topo)/2048*100:.2f}% density)",
                f"  Indices (first 15):     {active_topo[:15]}...",
                f"MACCS Keys (166 bits):     {len(active_maccs)} active keys",
                f"  Active Key Indices:     {active_maccs}",
            ]
            self.txt_inspect_result.setText("\n".join(res))
        except Exception as e:
            QMessageBox.critical(self.widget, "Error", f"Failed to compute single molecule fingerprints:\n{e}")

    def compare_smiles_tanimoto(self):
        smiles1 = self.txt_inspect_smiles.text().strip()
        smiles2 = self.txt_ref_smiles.text().strip()

        if not smiles1 or not smiles2:
            QMessageBox.warning(self.widget, "Missing SMILES", "Please specify both SMILES strings.")
            return

        try:
            mol1 = parse_smiles(smiles1)
            mol2 = parse_smiles(smiles2)
            if not mol1 or not mol2:
                QMessageBox.warning(self.widget, "Parse Error", "Failed to parse one of the SMILES strings.")
                return

            fp1 = get_morgan_fingerprint(mol1, radius=2, n_bits=2048)
            fp2 = get_morgan_fingerprint(mol2, radius=2, n_bits=2048)

            sim = tanimoto_similarity(fp1, fp2)
            self.lbl_similarity_res.setText(f"Morgan Tanimoto Similarity: {sim:.4f}")
        except Exception as e:
            QMessageBox.critical(self.widget, "Error", f"Similarity calculation error:\n{e}")

    def _update_matrix_tab_display(self, tanimoto_df: pd.DataFrame):
        if tanimoto_df is None or len(tanimoto_df) == 0:
            return

        matrix_vals = tanimoto_df.iloc[:, 1:].values.astype(float)
        np.fill_diagonal(matrix_vals, np.nan)
        mean_sim = float(np.nanmean(matrix_vals)) if not np.isnan(matrix_vals).all() else 1.0
        min_sim = float(np.nanmin(matrix_vals)) if not np.isnan(matrix_vals).all() else 1.0
        max_sim = float(np.nanmax(matrix_vals)) if not np.isnan(matrix_vals).all() else 1.0
        diversity = 1.0 - mean_sim

        self.lbl_diversity_stats.setText(
            f"<b>Dataset Size:</b> {len(tanimoto_df)} compounds | "
            f"<b>Mean Tanimoto:</b> {mean_sim:.4f} | "
            f"<b>Min:</b> {min_sim:.4f} | <b>Max:</b> {max_sim:.4f} | "
            f"<b>Dataset Diversity Index (1 - Mean):</b> <font color='#00bcd4'>{diversity:.4f}</font>"
        )

        n_rows = min(30, len(tanimoto_df))
        n_cols = min(30, tanimoto_df.shape[1])
        self.tbl_tanimoto.setRowCount(n_rows)
        self.tbl_tanimoto.setColumnCount(n_cols)
        self.tbl_tanimoto.setHorizontalHeaderLabels(list(tanimoto_df.columns[:n_cols]))

        for r in range(n_rows):
            for c in range(n_cols):
                val = tanimoto_df.iloc[r, c]
                txt = f"{val:.3f}" if isinstance(val, (float, np.floating)) else str(val)
                self.tbl_tanimoto.setItem(r, c, QTableWidgetItem(txt))

    def export_tanimoto_matrix(self):
        if self.latest_tanimoto_df is None:
            QMessageBox.warning(self.widget, "No Matrix", "No Tanimoto matrix currently loaded. Check 'Calculate Tanimoto Matrix' in Tab 1 first.")
            return

        path, _ = QFileDialog.getSaveFileName(self.widget, "Export Tanimoto Matrix", "tanimoto_matrix.csv", "CSV Files (*.csv)")
        if path:
            self.latest_tanimoto_df.to_csv(path, index=False)
            QMessageBox.information(self.widget, "Success", f"Matrix exported to:\n{path}")

    def on_molecule_changed(self, molecule):
        if molecule:
            smiles = getattr(molecule, 'smiles', '')
            if smiles:
                self.txt_inspect_smiles.setText(smiles)
                self.inspect_single_smiles()

    # -------------------------------------------------------------------------
    # Docking / Window Popout Controls
    # -------------------------------------------------------------------------
    def _open_popout(self, fullscreen: bool = False) -> None:
        if getattr(self, '_active_popout', None) is not None:
            self._active_popout.raise_()
            self._active_popout.activateWindow()
            return

        dialog = QDialog()
        dialog.setWindowTitle("PyChem-Pro Plugin - Fingerprint Acquisition & QSAR Prep")
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


# =============================================================================
# Plugin Core Class
# =============================================================================

class CalculateFingerprintsPlugin(BasePlugin):
    """
    Fingerprint Acquisition & QSAR Prep Plugin for PyChem-Pro.
    Subclasses BasePlugin for discovery by PyChem-Pro PluginManager.
    """

    def __init__(self):
        info = PluginInfo(
            name="Fingerprint Acquisition & QSAR Prep",
            version="2.0.0",
            description="Calculates molecular fingerprints (Morgan, Topological, MACCS, Combined) with QSAR bit pruning, Tanimoto similarity matrices, and activity column preservation.",
            author="PyChem-Pro Team",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        )
        super().__init__(info)
        self.current_molecule = None
        self.widget = None

    def get_info(self) -> PluginInfo:
        return self.info

    def create_widget(self) -> CalculateFingerprintsWidget:
        if self.widget is None:
            self.widget = CalculateFingerprintsWidget(self)
        return self.widget

    def initialize(self, main_window=None, api=None):
        self._main_window = main_window
        self._api = api
        self._is_initialized = True
        self.logger.info("Fingerprint Acquisition & QSAR Prep Plugin initialized successfully.")
        return True

    def cleanup(self):
        if self.widget:
            if hasattr(self.widget, 'widget') and self.widget.widget:
                self.widget.widget.deleteLater()
            self.widget = None
        self.logger.info("Fingerprint Acquisition & QSAR Prep Plugin cleaned up.")

    def on_molecule_changed(self, molecule):
        self.current_molecule = molecule
        if self.widget:
            self.widget.on_molecule_changed(molecule)


# =============================================================================
# CLI Execution Mode
# =============================================================================

def run_cli(args):
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' does not exist.")
        sys.exit(1)

    print("=" * 70)
    print("      PyChem-Pro Fingerprint Acquisition & QSAR Prep Utility")
    print("=" * 70)
    print(f"Input file:          {args.input}")
    print(f"Fingerprint type:    {FP_TYPES.get(args.type, args.type)}")
    print(f"Variance threshold:  {args.variance_thresh}")
    print(f"Frequency filter:    [{args.min_freq*100:.1f}%, {args.max_freq*100:.1f}%]")

    try:
        if args.input.endswith(('.xlsx', '.xls')):
            df = pd.read_excel(args.input)
        else:
            with open(args.input, 'r', encoding='utf-8', errors='ignore') as f:
                first_line = f.readline()
                sep = '\t' if '\t' in first_line else ','
            df = pd.read_csv(args.input, sep=sep)
    except Exception as e:
        print(f"Error loading file: {e}")
        sys.exit(1)

    print(f"Loaded dataset containing {len(df)} rows.")

    # Match SMILES column
    smiles_col = args.smiles_col
    if smiles_col not in df.columns:
        cols_lower = [c.lower() for c in df.columns]
        if 'smiles' in cols_lower:
            smiles_col = df.columns[cols_lower.index('smiles')]
        else:
            print(f"Error: SMILES column '{args.smiles_col}' not found.")
            sys.exit(1)

    # Match ID column
    id_col = args.id_col
    if id_col not in df.columns:
        cols_lower = [c.lower() for c in df.columns]
        matched = False
        for target in ['immpat id', 'id', 'name', 'pubchem id', 'mol_id']:
            if target in cols_lower:
                id_col = df.columns[cols_lower.index(target)]
                matched = True
                break
        if not matched:
            df['Temp_ID'] = [f"MOL_{i+1}" for i in range(len(df))]
            id_col = 'Temp_ID'

    df = df.dropna(subset=[smiles_col])
    df = df[df[smiles_col].astype(str).str.strip() != '']
    valid_len = len(df)

    smiles_list = df[smiles_col].astype(str).tolist()
    ids_list = df[id_col].astype(str).tolist()

    n_jobs = args.jobs if args.jobs and args.jobs > 0 else max(1, multiprocessing.cpu_count() // 2)
    print(f"Processing {valid_len} molecules using {n_jobs} parallel processes...")

    worker_args = [(idx, s, name, args.type, args.radius, args.bits) for idx, (s, name) in enumerate(zip(smiles_list, ids_list))]
    results = [None] * len(worker_args)

    start_time = time.time()
    with multiprocessing.Pool(processes=n_jobs) as pool:
        for idx, vec, err in pool.imap_unordered(_compute_single_smiles, worker_args):
            results[idx] = vec

    elapsed = time.time() - start_time
    success_rows = [v for v in results if v is not None]

    if len(success_rows) == 0:
        print("Error: No molecules were successfully processed.")
        sys.exit(1)

    raw_matrix = np.array(success_rows, dtype=np.uint8)
    pruned_matrix, retained_bits, stats = prune_fingerprint_matrix(
        raw_matrix,
        variance_thresh=args.variance_thresh,
        min_freq=args.min_freq,
        max_freq=args.max_freq,
        remove_duplicates=args.remove_duplicates
    )

    col_names = [f"FP_bit{b}" for b in retained_bits]
    output_df = pd.DataFrame(pruned_matrix, columns=col_names)
    output_df.insert(0, 'SMILES', smiles_list[:len(success_rows)])
    output_df.insert(0, id_col, ids_list[:len(success_rows)])

    if args.split_col and args.split_col in df.columns:
        output_df.insert(2, args.split_col, df[args.split_col].values[:len(success_rows)])

    if args.target_col and args.target_col in df.columns:
        output_df.insert(3 if (args.split_col and args.split_col in output_df.columns) else 2, args.target_col, df[args.target_col].values[:len(success_rows)])

    sep = '\t' if args.output.endswith('.tsv') else ','
    output_df.to_csv(args.output, index=False, sep=sep)

    print("-" * 70)
    print(f"Elapsed Time:      {elapsed:.2f} sec")
    print(f"Original Bits:     {stats['original_bits']}")
    print(f"Retained QSAR Bits:{stats['retained_bits']} (Pruned {stats['removed_bits']} bits)")
    print(f"Output Saved To:   {args.output}")

    if args.tanimoto_matrix:
        base, ext = os.path.splitext(args.output)
        tan_path = f"{base}_tanimoto_matrix{ext}"
        tan_mat = compute_tanimoto_matrix(pruned_matrix)
        tan_df = pd.DataFrame(tan_mat, index=ids_list[:len(success_rows)], columns=ids_list[:len(success_rows)])
        tan_df.insert(0, id_col, ids_list[:len(success_rows)])
        tan_df.to_csv(tan_path, index=False, sep=sep)
        print(f"Tanimoto Matrix:   {tan_path}")
    print("=" * 70)


# =============================================================================
# Entry Point
# =============================================================================

def main():
    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser(description="PyChem-Pro Fingerprint Acquisition & QSAR Prep CLI")
        parser.add_argument('-i', '--input', required=True, help="Input CSV/TSV/XLSX file")
        parser.add_argument('-o', '--output', required=True, help="Output fingerprint CSV file")
        parser.add_argument('-t', '--type', default='morgan', choices=list(FP_TYPES.keys()))
        parser.add_argument('-r', '--radius', type=int, default=2, help="Morgan radius")
        parser.add_argument('-b', '--bits', type=int, default=2048, help="Bit vector size")
        parser.add_argument('-s', '--smiles-col', default='SMILES', help="SMILES column name")
        parser.add_argument('-id', '--id-col', default='IMMPAT ID', help="ID column name")
        parser.add_argument('--target-col', default=None, help="Activity target column to preserve")
        parser.add_argument('--split-col', default=None, help="Train/Test Split column to preserve")
        parser.add_argument('--variance-thresh', type=float, default=0.0, help="Variance filter threshold")
        parser.add_argument('--min-freq', type=float, default=0.0, help="Min bit frequency (0.0 to 1.0)")
        parser.add_argument('--max-freq', type=float, default=1.0, help="Max bit frequency (0.0 to 1.0)")
        parser.add_argument('--remove-duplicates', action='store_true', help="Remove identical duplicate bit columns")
        parser.add_argument('--tanimoto-matrix', action='store_true', help="Export Tanimoto similarity matrix")
        parser.add_argument('-j', '--jobs', type=int, default=None, help="Parallel CPU cores")

        args = parser.parse_args()
        run_cli(args)
    else:
        app = QApplication(sys.argv)
        gui = QMainWindow()
        gui.setWindowTitle("PyChem-Pro Fingerprint Acquisition & QSAR Prep Engine")
        w = CalculateFingerprintsWidget()
        gui.setCentralWidget(w.widget)
        gui.resize(900, 650)
        gui.show()
        sys.exit(app.exec())


if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()
