"""Ramachandran Plot Plugin with KDE density regions."""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from itertools import cycle
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from enum import Enum

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QFileDialog, QMessageBox, Qt, QListWidget, QListWidgetItem,
    QGroupBox, QCheckBox, QTableWidget, QTableWidgetItem, QHeaderView,
    QSpinBox, QGridLayout, QFrame, QTabWidget, QSplitter, QTextEdit,
    QScrollArea
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType


class RegionType(Enum):
    """Classification based on KDE density."""
    FAVORED = "favored"
    ALLOWED = "allowed"
    GENEROUS = "generous"
    OUTLIER = "outlier"
    NOT_CALCULATED = "not_calculated"


DEFAULT_DENSITY_THRESHOLDS = {
    RegionType.FAVORED: -1.0,
    RegionType.ALLOWED: -2.0,
    RegionType.GENEROUS: -3.5,
}


@dataclass
class ResidueInfo:
    """Container for residue information with phi/psi angles."""
    chain: str
    name: str
    number: int
    phi: Optional[float]
    psi: Optional[float]
    region: RegionType = RegionType.NOT_CALCULATED
    
    @property
    def res_id(self) -> str:
        return f"{self.chain}:{self.name}{self.number}"
    
    @property
    def has_complete_angles(self) -> bool:
        return self.phi is not None and self.psi is not None


@dataclass
class StructureData:
    """Container for parsed PDB structure data."""
    filepath: str
    residues: List[ResidueInfo]
    color: Tuple[float, float, float, float]
    
    @property
    def phi_psi_pairs(self) -> List[Tuple[float, float]]:
        return [(r.phi, r.psi) for r in self.residues if r.has_complete_angles]
    
    @property
    def region_counts(self) -> Dict[RegionType, int]:
        counts = {rt: 0 for rt in RegionType}
        for r in self.residues:
            if r.has_complete_angles:
                counts[r.region] += 1
        return counts


def calculate_dihedral(p0: np.ndarray, p1: np.ndarray,
                       p2: np.ndarray, p3: np.ndarray) -> float:
    """Calculate dihedral angle in degrees for four points."""
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    
    b1_norm = np.linalg.norm(b1)
    if b1_norm < 1e-10:
        return 0.0
    b1 = b1 / b1_norm
    
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    
    v_norm = np.linalg.norm(v)
    w_norm = np.linalg.norm(w)
    if v_norm < 1e-10 or w_norm < 1e-10:
        return 0.0
    
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    
    return float(np.degrees(np.arctan2(y, x)))


class PDBParser:
    """Robust PDB parser for backbone atom extraction."""
    
    def __init__(self):
        self.warnings: List[str] = []
    
    def parse(self, filepath: str) -> Dict[str, Dict]:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"PDB file not found: {filepath}")
        
        residues: Dict[str, Dict] = {}
        self.warnings = []
        
        with open(filepath, "r") as f:
            for line_num, line in enumerate(f, 1):
                if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
                    continue
                try:
                    atom_name = line[12:16].strip()
                    if atom_name not in ["N", "CA", "C"]:
                        continue
                    
                    res_name = line[17:20].strip()
                    chain_id = line[21].strip() or "A"
                    res_num_str = line[22:26].strip()
                    
                    try:
                        res_num = int(res_num_str)
                    except ValueError:
                        res_num = int("".join(c for c in res_num_str if c.isdigit()))
                    
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    res_id = f"{chain_id}:{res_name}{res_num}"
                    
                    if res_id not in residues:
                        residues[res_id] = {
                            "atoms": {},
                            "name": res_name,
                            "num": res_num,
                            "chain": chain_id
                        }
                    residues[res_id]["atoms"][atom_name] = np.array([x, y, z])
                except (ValueError, IndexError) as e:
                    self.warnings.append(f"Line {line_num}: {e}")
                    continue
        
        if not residues:
            raise ValueError("No valid backbone atoms found")
        return residues


def extract_phi_psi(residues: Dict[str, Dict]) -> List[ResidueInfo]:
    """Calculate phi/psi angles from parsed residue data."""
    sorted_keys = sorted(
        residues.keys(),
        key=lambda x: (residues[x]["chain"], residues[x]["num"])
    )
    
    result: List[ResidueInfo] = []
    
    for i, key in enumerate(sorted_keys):
        current = residues[key]
        atoms = current["atoms"]
        
        if not all(k in atoms for k in ["N", "CA", "C"]):
            continue
        
        phi = None
        psi = None
        
        if i > 0:
            prev_key = sorted_keys[i - 1]
            prev = residues[prev_key]
            if (prev["chain"] == current["chain"] and
                prev["num"] == current["num"] - 1 and
                "C" in prev["atoms"]):
                phi = calculate_dihedral(
                    prev["atoms"]["C"],
                    atoms["N"],
                    atoms["CA"],
                    atoms["C"]
                )
        
        if i < len(sorted_keys) - 1:
            next_key = sorted_keys[i + 1]
            next_res = residues[next_key]
            if (next_res["chain"] == current["chain"] and
                next_res["num"] == current["num"] + 1 and
                "N" in next_res["atoms"]):
                psi = calculate_dihedral(
                    atoms["N"],
                    atoms["CA"],
                    atoms["C"],
                    next_res["atoms"]["N"]
                )
        
        result.append(ResidueInfo(
            chain=current["chain"],
            name=current["name"],
            number=current["num"],
            phi=phi,
            psi=psi
        ))
    
    return result


class KdeDensityMap:
    """Handler for KDE density map from Top8000 dataset."""
    
    def __init__(self, kde_path: Optional[str] = None):
        if kde_path is None:
            base_dir = os.path.dirname(os.path.dirname(__file__))
            kde_path = os.path.join(base_dir, "src", "data", "kde.dat")
        
        self.kde_path = kde_path
        self.thresholds = DEFAULT_DENSITY_THRESHOLDS.copy()
        self._density_grid: Optional[np.ndarray] = None
        self._log_density: Optional[np.ndarray] = None
        self._load_kde()
    
    def _load_kde(self) -> None:
        try:
            with open(self.kde_path, "rb") as f:
                raw_data = f.read()
            
            z = np.frombuffer(raw_data, dtype=np.float64)
            if z.size != 10000:
                z = np.frombuffer(raw_data, dtype=np.float32)
            
            # KDE data is typically stored in (phi, psi) order.
            # Transpose to (psi, phi) so indexing matches [psi_idx, phi_idx]
            # which maps to [row, col] -> (Y, X).
            self._density_grid = np.reshape(z, (100, 100)).T
            self._log_density = np.log10(self._density_grid + 1e-10)
            
            # Auto-adjust thresholds if data range is unexpected
            max_log = np.max(self._log_density)
            if max_log < -1.0:
                # If peak density is lower than -1.0, traditional thresholds will fail.
                # Shift relative to the peak: Fav=peak-0.5, Allw=peak-1.5, Gen=peak-3.0
                self.thresholds[RegionType.FAVORED] = max_log - 1.0
                self.thresholds[RegionType.ALLOWED] = max_log - 2.5
                self.thresholds[RegionType.GENEROUS] = max_log - 4.5
        except Exception as e:
            raise RuntimeError(f"Failed to load KDE: {e}")
    
    @property
    def density_grid(self) -> np.ndarray:
        if self._density_grid is None:
            raise RuntimeError("KDE not loaded")
        return self._density_grid
    
    @property
    def log_density(self) -> np.ndarray:
        if self._log_density is None:
            raise RuntimeError("KDE not loaded")
        return self._log_density
    
    def get_density_at(self, phi: float, psi: float) -> float:
        phi_idx = int((phi + 180) / 3.6)
        psi_idx = int((psi + 180) / 3.6)
        phi_idx = max(0, min(99, phi_idx))
        psi_idx = max(0, min(99, psi_idx))
        return float(self._log_density[psi_idx, phi_idx])
    
    def classify_point(self, phi: float, psi: float) -> RegionType:
        density = self.get_density_at(phi, psi)
        
        # Debug first few classifications
        if not hasattr(self, '_debug_count'): self._debug_count = 0
        if self._debug_count < 5:
            print(f"DEBUG: point ({phi:.1f}, {psi:.1f}) density={density:.3f} thresholds={self.thresholds}")
            self._debug_count += 1
            
        if density > self.thresholds[RegionType.FAVORED]:
            return RegionType.FAVORED
        elif density > self.thresholds[RegionType.ALLOWED]:
            return RegionType.ALLOWED
        elif density > self.thresholds[RegionType.GENEROUS]:
            return RegionType.GENEROUS
        return RegionType.OUTLIER


class RamachandranWidget(PluginWidget):
    """Enhanced widget for Ramachandran Plot analysis."""
    
    def __init__(self, plugin):
        super().__init__(plugin)
        self.pdb_files: List[str] = []
        self.structure_data: Dict[str, StructureData] = {}
        self.colors = cycle(plt.cm.tab10.colors)
        
        try:
            self.kde_map = KdeDensityMap()
        except RuntimeError as e:
            QMessageBox.warning(None, "KDE Error", str(e))
            self.kde_map = None
        
        self.setup_ui()
    
    def setup_ui(self):
        self.widget = QWidget()
        self.widget.setStyleSheet(self._get_stylesheet())
        
        main_layout = QHBoxLayout(self.widget)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(15)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setWidget(self._create_control_panel())
        
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(scroll)
        splitter.addWidget(self._create_plot_panel())
        splitter.setSizes([300, 900])
        main_layout.addWidget(splitter)
        
        self.update_plot()
    
    def cleanup(self):
        """Explicitly close matplotlib figure to release resources."""
        try:
            if hasattr(self, 'figure'):
                plt.close(self.figure)
            self.structure_data.clear()
            if hasattr(self, 'plugin') and self.plugin and hasattr(self.plugin, 'logger'):
                self.plugin.logger.info("RamachandranWidget cleaned up and figure closed.")
            else:
                print("RamachandranWidget cleaned up and figure closed.")
        except Exception as e:
            if hasattr(self, 'plugin') and self.plugin and hasattr(self.plugin, 'logger'):
                self.plugin.logger.error(f"Error during RamachandranWidget cleanup: {e}")
            else:
                print(f"Error during RamachandranWidget cleanup: {e}")

    def _get_stylesheet(self):
        return """
            QWidget { background-color: #ffffff; color: #1e293b; font-family: 'Segoe UI'; font-size: 11px; }
            QGroupBox { border: 1px solid #e2e8f0; border-radius: 6px; margin-top: 10px;
                       padding-top: 8px; font-weight: bold; color: #4f46e5; }
            QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 4px; }
            QPushButton { background-color: #4f46e5; border: none; border-radius: 4px;
                         color: white; padding: 4px 8px; font-weight: bold; min-height: 24px; }
            QPushButton:hover { background-color: #6366f1; }
            QPushButton#danger_btn { background-color: #ef4444; }
            QPushButton#danger_btn:hover { background-color: #f87171; }
            QPushButton#secondary_btn { background-color: #64748b; }
            QListWidget { background-color: #f8fafc; border: 1px solid #e2e8f0;
                         border-radius: 4px; padding: 2px; }
            QCheckBox::indicator { width: 14px; height: 14px; border-radius: 3px;
                                  border: 1px solid #cbd5e1; }
            QCheckBox::indicator:checked { background-color: #4f46e5; }
            QTableWidget { background-color: #ffffff; gridline-color: #e2e8f0;
                         border-radius: 4px; font-size: 9px; }
            QHeaderView::section { background-color: #f1f5f9; font-weight: bold; padding: 2px; }
            QTabWidget::pane { border: 1px solid #e2e8f0; border-radius: 4px; }
            QTabBar::tab { background: #f1f5f9; border: 1px solid #e2e8f0;
                          padding: 4px 8px; border-top-left-radius: 4px;
                          border-top-right-radius: 4px; }
            QTabBar::tab:selected { background: #4f46e5; color: white; }
            QTextEdit { background-color: #f8fafc; border: 1px solid #e2e8f0;
                       border-radius: 4px; font-family: Consolas; font-size: 9px; }
        """
    
    def _create_control_panel(self):
        panel = QFrame()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(10, 10, 10, 10)
        panel.setFixedWidth(320)
        
        title = QLabel("RAMACHANDRAN ANALYSIS")
        title.setStyleSheet("font-size: 18px; font-weight: 800; color: #1e293b;")
        layout.addWidget(title)
        
        file_group = QGroupBox("STRUCTURE MANAGEMENT")
        file_layout = QVBoxLayout(file_group)
        
        self.file_list = QListWidget()
        self.file_list.setMaximumHeight(120)
        file_layout.addWidget(self.file_list)
        
        btn_layout = QHBoxLayout()
        self.add_btn = QPushButton("Add PDB")
        self.add_btn.clicked.connect(self.on_add_pdb)
        btn_layout.addWidget(self.add_btn)
        
        self.remove_btn = QPushButton("Remove")
        self.remove_btn.setObjectName("danger_btn")
        self.remove_btn.clicked.connect(self.on_remove_pdb)
        btn_layout.addWidget(self.remove_btn)
        file_layout.addLayout(btn_layout)
        
        self.clear_btn = QPushButton("Clear All")
        self.clear_btn.setObjectName("secondary_btn")
        self.clear_btn.clicked.connect(self.on_clear_all)
        file_layout.addWidget(self.clear_btn)
        
        layout.addWidget(file_group)
        
        options_group = QGroupBox("DISPLAY OPTIONS")
        options_grid = QGridLayout(options_group)
        
        self.show_kde_check = QCheckBox("Show KDE Density")
        self.show_kde_check.setChecked(True)
        self.show_kde_check.stateChanged.connect(self.update_plot)
        options_grid.addWidget(self.show_kde_check, 0, 0)
        
        self.show_contours_check = QCheckBox("Show Contours")
        self.show_contours_check.setChecked(True)
        self.show_contours_check.stateChanged.connect(self.update_plot)
        options_grid.addWidget(self.show_contours_check, 0, 1)
        
        self.show_points_check = QCheckBox("Show Points")
        self.show_points_check.setChecked(True)
        self.show_points_check.stateChanged.connect(self.update_plot)
        options_grid.addWidget(self.show_points_check, 1, 0)
        
        self.color_by_region_check = QCheckBox("Color by Region")
        self.color_by_region_check.setChecked(False)
        self.color_by_region_check.stateChanged.connect(self.update_plot)
        options_grid.addWidget(self.color_by_region_check, 1, 1)
        
        self.highlight_outliers_check = QCheckBox("Highlight Outliers")
        self.highlight_outliers_check.setChecked(True)
        self.highlight_outliers_check.setToolTip("Mark outliers in bright red DIAMONDS for visibility.")
        self.highlight_outliers_check.stateChanged.connect(self.update_plot)
        options_grid.addWidget(self.highlight_outliers_check, 2, 0)
        
        layout.addWidget(options_group)
        
        stats_group = QGroupBox("REGION STATISTICS (%)")
        stats_layout = QVBoxLayout(stats_group)
        
        self.stats_table = QTableWidget(0, 5)
        self.stats_table.setHorizontalHeaderLabels(["Structure", "Fav", "Allw", "Gen", "Out"])
        self.stats_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.stats_table.setMaximumHeight(150)
        self.stats_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.stats_table.setShowGrid(False)
        self.stats_table.verticalHeader().setVisible(False)
        stats_layout.addWidget(self.stats_table)
        
        layout.addWidget(stats_group)
        
        export_group = QGroupBox("EXPORT")
        export_layout = QGridLayout(export_group)
        
        self.dpi_spin = QSpinBox()
        self.dpi_spin.setRange(100, 2400)
        self.dpi_spin.setValue(600)
        self.dpi_spin.setSingleStep(100)
        export_layout.addWidget(QLabel("DPI:"), 0, 0)
        export_layout.addWidget(self.dpi_spin, 0, 1)
        
        self.export_plot_btn = QPushButton("Export Plot")
        self.export_plot_btn.clicked.connect(self.on_export_plot)
        export_layout.addWidget(self.export_plot_btn, 1, 0)
        
        self.export_csv_btn = QPushButton("Export CSV")
        self.export_csv_btn.clicked.connect(self.on_export_csv)
        export_layout.addWidget(self.export_csv_btn, 1, 1)
        
        self.export_report_btn = QPushButton("Report")
        self.export_report_btn.setObjectName("secondary_btn")
        self.export_report_btn.clicked.connect(self.on_export_report)
        export_layout.addWidget(self.export_report_btn, 2, 0, 1, 2)
        
        layout.addWidget(export_group)
        
        help_text = QLabel("Double-click plot to label residues. Right-click to clear.")
        help_text.setStyleSheet("color: #64748b; font-size: 10px;")
        help_text.setWordWrap(True)
        layout.addWidget(help_text)
        
        layout.addStretch()
        return panel
    
    def _create_plot_panel(self):
        panel = QFrame()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.tabs = QTabWidget()
        
        plot_tab = QWidget()
        plot_layout = QVBoxLayout(plot_tab)
        plot_layout.setContentsMargins(10, 10, 10, 10)
        
        self.figure = Figure(figsize=(8, 8), facecolor='#ffffff')
        self.figure.set_dpi(100)
        
        import matplotlib as mpl
        mpl.rcParams['path.simplify'] = False
        
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setStyleSheet("background-color: transparent;")
        self.ax = self.figure.add_subplot(111)
        self.canvas.mpl_connect('button_press_event', self.on_plot_click)
        
        plot_layout.addWidget(self.canvas)
        self.tabs.addTab(plot_tab, "Plot")
        
        details_tab = QWidget()
        details_layout = QVBoxLayout(details_tab)
        self.details_text = QTextEdit()
        self.details_text.setReadOnly(True)
        details_layout.addWidget(self.details_text)
        self.tabs.addTab(details_tab, "Details")
        
        layout.addWidget(self.tabs)
        return panel
    
    def on_add_pdb(self):
        files, _ = QFileDialog.getOpenFileNames(
            self.widget, "Open PDB Files", "",
            "PDB Files (*.pdb *.ent);;All Files (*.*)"
        )
        
        if not files:
            return
        
        parser = PDBParser()
        
        for filepath in files:
            if filepath in self.pdb_files:
                continue
            
            try:
                residues = parser.parse(filepath)
                residue_infos = extract_phi_psi(residues)
                
                if not residue_infos:
                    QMessageBox.warning(self.widget, "Warning",
                        f"No valid residues in {os.path.basename(filepath)}")
                    continue
                
                if self.kde_map:
                    for r in residue_infos:
                        if r.has_complete_angles:
                            r.region = self.kde_map.classify_point(r.phi, r.psi)
                
                color = next(self.colors)
                self.structure_data[filepath] = StructureData(
                    filepath=filepath,
                    residues=residue_infos,
                    color=color
                )
                
                self.pdb_files.append(filepath)
                self.file_list.addItem(os.path.basename(filepath))
                
            except Exception as e:
                QMessageBox.critical(self.widget, "Error",
                    f"Failed to parse {os.path.basename(filepath)}:\n{str(e)}")
        
        self.update_plot()
        self.update_details()
    
    def on_remove_pdb(self):
        for item in self.file_list.selectedItems():
            basename = item.text()
            for filepath in list(self.pdb_files):
                if os.path.basename(filepath) == basename:
                    self.pdb_files.remove(filepath)
                    del self.structure_data[filepath]
                    break
            self.file_list.takeItem(self.file_list.row(item))
        self.update_plot()
        self.update_details()
    
    def on_clear_all(self):
        self.pdb_files.clear()
        self.structure_data.clear()
        self.file_list.clear()
        self.update_plot()
        self.update_details()
    
    def on_export_plot(self):
        file_path, _ = QFileDialog.getSaveFileName(
            self.widget, "Save Plot", "ramachandran_plot.png",
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)"
        )
        if not file_path:
            return
        
        try:
            self.figure.savefig(file_path, dpi=self.dpi_spin.value(),
                              bbox_inches='tight', facecolor='white')
            QMessageBox.information(self.widget, "Success",
                f"Plot saved to:\n{file_path}")
        except Exception as e:
            QMessageBox.critical(self.widget, "Error", str(e))
    
    def on_export_csv(self):
        if not self.structure_data:
            QMessageBox.information(self.widget, "Info", "No data to export")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self.widget, "Export CSV", "ramachandran_data.csv",
            "CSV (*.csv)"
        )
        if not file_path:
            return
        
        try:
            with open(file_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Structure', 'Chain', 'Residue', 'Number',
                                'Phi', 'Psi', 'Region'])
                
                for filepath, data in self.structure_data.items():
                    name = os.path.basename(filepath)
                    for r in data.residues:
                        writer.writerow([
                            name, r.chain, r.name, r.number,
                            f"{r.phi:.2f}" if r.phi else "",
                            f"{r.psi:.2f}" if r.psi else "",
                            r.region.value
                        ])
            
            QMessageBox.information(self.widget, "Success",
                f"Data exported to:\n{file_path}")
        except Exception as e:
            QMessageBox.critical(self.widget, "Error", str(e))
    
    def on_export_report(self):
        if not self.structure_data:
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self.widget, "Save Report", "ramachandran_report.txt",
            "Text (*.txt)"
        )
        if not file_path:
            return
        
        try:
            with open(file_path, 'w') as f:
                f.write("=" * 60 + "\n")
                f.write("RAMACHANDRAN ANALYSIS REPORT\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"Structures: {len(self.structure_data)}\n")
                f.write("KDE Map: Top8000 dataset\n\n")
                
                for filepath, data in self.structure_data.items():
                    f.write("-" * 60 + "\n")
                    f.write(f"Structure: {os.path.basename(filepath)}\n")
                    f.write("-" * 60 + "\n\n")
                    
                    stats = data.region_counts
                    total = sum(stats.values()) or 1
                    
                    f.write(f"Total: {total}\n")
                    for rt in RegionType:
                        pct = stats[rt] / total * 100
                        f.write(f"  {rt.value.capitalize():10s}: {stats[rt]:4d} ({pct:5.1f}%)\n")
                    f.write("\n")
            
            QMessageBox.information(self.widget, "Success",
                f"Report saved to:\n{file_path}")
        except Exception as e:
            QMessageBox.critical(self.widget, "Error", str(e))
    
    def on_plot_click(self, event):
        if not event.dblclick:
            if event.button == 3:
                self.update_plot()
            return
        
        if event.xdata is None or event.ydata is None:
            return
        
        best_dist = float('inf')
        best_residue = None
        
        for data in self.structure_data.values():
            for r in data.residues:
                if not r.has_complete_angles:
                    continue
                dist = (r.phi - event.xdata)**2 + (r.psi - event.ydata)**2
                if dist < best_dist:
                    best_dist = dist
                    best_residue = r
        
        if best_dist < 100 and best_residue:
            self.ax.annotate(
                best_residue.res_id,
                xy=(best_residue.phi, best_residue.psi),
                xytext=(10, 10), textcoords='offset points',
                fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.7),
                arrowprops=dict(arrowstyle='->')
            )
            self.canvas.draw()
            self._show_residue_details(best_residue)
    
    def _plot_kde_density(self):
        if not self.kde_map or not self.show_kde_check.isChecked():
            return
        
        density = self.kde_map.density_grid
        log_density = self.kde_map.log_density
        
        display_data = log_density
        self.ax.imshow(
            display_data, cmap=plt.cm.YlOrRd,
            extent=(-180, 180, -180, 180),
            alpha=0.3, aspect='auto', origin='lower',
            vmin=-7, vmax=0
        )
        
        if self.show_contours_check.isChecked():
            # Create coordinate meshgrid for contours
            phi_range = np.linspace(-180, 180, 100)
            psi_range = np.linspace(-180, 180, 100)
            X, Y = np.meshgrid(phi_range, psi_range)
            
            contour_data = density
            levels = [10**i for i in range(-7, 0)]
            self.ax.contour(
                X, Y, contour_data, 
                colors='black', linewidths=0.5,
                levels=levels, alpha=0.4
            )
    
    def _plot_points(self):
        if not self.show_points_check.isChecked():
            return
        
        self.stats_table.setRowCount(0)
        
        region_colors = {
            RegionType.FAVORED: '#22c55e',
            RegionType.ALLOWED: '#3b82f6',
            RegionType.GENEROUS: '#f59e0b',
            RegionType.OUTLIER: '#ef4444',
        }
        
        region_markers = {
            RegionType.FAVORED: 'o',      # Circle
            RegionType.ALLOWED: 's',      # Square
            RegionType.GENEROUS: '^',     # Triangle Up
            RegionType.OUTLIER: 'D',      # Diamond
        }
        
        for filepath in self.pdb_files:
            data = self.structure_data[filepath]
            base_color = data.color
            
            by_region = {rt: [] for rt in RegionType}
            for r in data.residues:
                if r.has_complete_angles:
                    by_region[r.region].append(r)
            
            for region in RegionType:
                residues = by_region[region]
                if not residues:
                    continue
                
                phis = [r.phi for r in residues]
                psis = [r.psi for r in residues]
                
                if self.color_by_region_check.isChecked():
                    color = region_colors[region]
                    alpha = 0.8
                else:
                    color = base_color
                    alpha = 0.7 if region != RegionType.OUTLIER else 0.5
                
                # Apply special highlighting for outliers if requested
                if region == RegionType.OUTLIER and self.highlight_outliers_check.isChecked():
                    color = '#ef4444'
                    alpha = 1.0
                    marker_size = 80  # Prominent size
                    linewidth = 1.2
                    z_order = 15
                else:
                    marker_size = 40 if region == RegionType.OUTLIER else 30
                    linewidth = 0.6
                    z_order = 10 if region != RegionType.OUTLIER else 11
                
                marker = region_markers[region]
                
                self.ax.scatter(
                    phis, psis, s=marker_size, c=[color], alpha=alpha,
                    marker=marker, edgecolors='white', linewidths=linewidth,
                    zorder=z_order
                )
            
            stats = data.region_counts
            total = sum(stats.values()) or 1
            
            row = self.stats_table.rowCount()
            self.stats_table.insertRow(row)
            self.stats_table.setItem(row, 0, QTableWidgetItem(os.path.basename(filepath)[:12]))
            for i, rt in enumerate([RegionType.FAVORED, RegionType.ALLOWED,
                                   RegionType.GENEROUS, RegionType.OUTLIER], 1):
                self.stats_table.setItem(row, i, QTableWidgetItem(f"{stats[rt]/total*100:.1f}"))
        
        if len(self.pdb_files) > 1 and not self.color_by_region_check.isChecked():
            legend_elements = [
                Line2D([0], [0], marker='o', color='w',
                       label=os.path.basename(fp)[:12],
                       markerfacecolor=self.structure_data[fp].color, markersize=10)
                for fp in self.pdb_files
            ]
            self.ax.legend(handles=legend_elements, loc='upper right',
                          fontsize=9, framealpha=0.95)
        elif self.color_by_region_check.isChecked():
            region_legend = [
                Line2D([0], [0], marker='o', color='w', label='Favored',
                       markerfacecolor=region_colors[RegionType.FAVORED], markersize=8),
                Line2D([0], [0], marker='s', color='w', label='Allowed',
                       markerfacecolor=region_colors[RegionType.ALLOWED], markersize=8),
                Line2D([0], [0], marker='^', color='w', label='Generous',
                       markerfacecolor=region_colors[RegionType.GENEROUS], markersize=8),
                Line2D([0], [0], marker='D', color='w', label='Outlier',
                       markerfacecolor=region_colors[RegionType.OUTLIER], markersize=8),
            ]
            self.ax.legend(handles=region_legend, loc='upper right',
                          fontsize=9, framealpha=0.95, title='Regions')
    
    def update_plot(self):
        self.ax.clear()
        self.ax.set_facecolor('#ffffff')
        self.figure.set_facecolor('#ffffff')
        
        self._plot_kde_density()
        
        ticks = [-180, -135, -90, -45, 0, 45, 90, 135, 180]
        self.ax.set_xlim(-180, 180)
        self.ax.set_ylim(-180, 180)
        self.ax.set_xticks(ticks)
        self.ax.set_yticks(ticks)
        
        self.ax.tick_params(colors='#334155', labelsize=10)
        self.ax.set_xlabel(r'$\phi$ (degrees)', fontsize=12, fontweight='bold')
        self.ax.set_ylabel(r'$\psi$ (degrees)', fontsize=12, fontweight='bold')
        
        for spine in self.ax.spines.values():
            spine.set_edgecolor('#64748b')
            spine.set_linewidth(1.0)
        
        self.ax.axhline(0, color='#94a3b8', lw=0.8, alpha=0.6)
        self.ax.axvline(0, color='#94a3b8', lw=0.8, alpha=0.6)
        self.ax.grid(color='#cbd5e1', alpha=0.3, ls='--', lw=0.6)
        
        self._plot_points()
        
        if self.structure_data:
            title = f"Ramachandran Plot ({len(self.structure_data)} structure(s))"
        else:
            title = "Ramachandran Plot - Load PDB files to begin"
        
        self.ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
        self.figure.tight_layout()
        self.canvas.draw()
    
    def update_details(self):
        if not self.structure_data:
            self.details_text.setPlainText("No structures loaded.\n\nClick 'Add PDB' to begin.")
            return
        
        text_parts = []
        for filepath, data in self.structure_data.items():
            text_parts.append(f"Structure: {os.path.basename(filepath)}")
            text_parts.append("=" * 50)
            
            stats = data.region_counts
            total = sum(stats.values())
            
            text_parts.append(f"\nTotal residues: {total}")
            for rt in RegionType:
                pct = stats[rt] / total * 100 if total else 0
                text_parts.append(f"  {rt.value.capitalize():10s}: {stats[rt]:4d} ({pct:5.1f}%)")
            
            outliers = [r for r in data.residues if r.region == RegionType.OUTLIER and r.has_complete_angles]
            if outliers:
                text_parts.append(f"\nOutliers ({len(outliers)}):")
                for r in outliers[:15]:
                    phi_val = f"{r.phi:.1f}" if r.phi is not None else "None"
                    psi_val = f"{r.psi:.1f}" if r.psi is not None else "None"
                    text_parts.append(f"  {r.res_id}: phi={phi_val}, psi={psi_val}")
                if len(outliers) > 15:
                    text_parts.append(f"  ... and {len(outliers) - 15} more")
            
            text_parts.append("\n" + "-" * 50 + "\n")
        
        self.details_text.setPlainText("\n".join(text_parts))
    
    def _show_residue_details(self, residue):
        density = self.kde_map.get_density_at(residue.phi, residue.psi) if self.kde_map else 0
        phi_val = f"{residue.phi:.2f}" if residue.phi is not None else "N/A"
        psi_val = f"{residue.psi:.2f}" if residue.psi is not None else "N/A"
        
        text = f"""
Selected Residue
{'=' * 40}
Residue: {residue.res_id}

Angles:
  Phi: {phi_val}
  Psi: {psi_val}

Classification: {residue.region.value.upper().replace('_', ' ')}
Density: {density:.3f} (log10)
"""
        self.details_text.setPlainText(text)
        self.tabs.setCurrentIndex(1)


class RamachandranPlugin(BasePlugin):
    """Enhanced Ramachandran plot plugin with KDE density regions."""
    
    def __init__(self):
        super().__init__(PluginInfo(
            name="Ramachandran Plotter",
            version="1.0.0",
            description="KDE-based Ramachandran analysis using Top8000 dataset.",
            author="PyChem Team",
            plugin_type=PluginType.ANALYSIS,
            keywords=["protein", "backbone", "dihedral", "ramachandran", "kde"]
        ))
        self._widget = None
    
    def create_widget(self):
        if self._widget is None:
            self._widget = RamachandranWidget(self)
        return self._widget
    
    def initialize(self, main_window=None, api=None):
        try:
            if main_window is not None and api is not None:
                super().initialize(main_window, api)
            else:
                self._main_window = main_window
                self._api = api
                self._is_initialized = True
            self.logger.info("Ramachandran Plotter v2.0.0 initialized")
            return True
        except Exception as e:
            self.logger.error(f"Init failed: {e}")
            return False
    
    def cleanup(self):
        if self._widget:
            try:
                self._widget.cleanup()
            except:
                pass
            self._widget = None
        self.logger.info("Ramachandran Plotter cleaned up")
