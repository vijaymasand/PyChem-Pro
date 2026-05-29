"""
Chemical Space Visualization - UI Engine
PyQt5 / PySide6 User Interface with Matplotlib integration.
"""
import sys
import numpy as np
import pandas as pd
from typing import Optional

from src.shared.qt_compat import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QComboBox, QCheckBox, QSlider, QGroupBox, QFileDialog, QSpinBox, QDoubleSpinBox,
    Qt, QThread, Signal, QMessageBox, QProgressDialog, QColorDialog, QSplitter, QTableWidget, QTableWidgetItem,
    QDialog, QFormLayout
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.colors as mcolors
import matplotlib.cm as cm

from .data_engine import DataEngine
from .featurization_engine import FeaturizationEngine
from .analysis_engine import AnalysisEngine

class WorkerThread(QThread):
    progress = Signal(int, str)
    finished = Signal(bool, str, object) # success, message, result

    def __init__(self, func, *args, **kwargs):
        super().__init__()
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def run(self):
        try:
            # We can pass self.progress.emit as a callback if the function supports it
            if 'progress_callback' in self.func.__code__.co_varnames:
                self.kwargs['progress_callback'] = self.progress.emit
                
            res = self.func(*self.args, **self.kwargs)
            self.finished.emit(True, "Complete", res)
        except Exception as e:
            self.finished.emit(False, str(e), None)


class ExportImageDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Export High-Res Image")
        self.layout = QFormLayout(self)
        
        self.dpi_spin = QSpinBox()
        self.dpi_spin.setRange(72, 1200)
        self.dpi_spin.setValue(300)
        self.layout.addRow("DPI:", self.dpi_spin)
        
        self.font_spin = QSpinBox()
        self.font_spin.setRange(6, 48)
        self.font_spin.setValue(12)
        self.layout.addRow("Font Size:", self.font_spin)
        
        self.point_spin = QSpinBox()
        self.point_spin.setRange(10, 500)
        self.point_spin.setValue(50)
        self.layout.addRow("Point Size:", self.point_spin)
        
        self.alpha_spin = QDoubleSpinBox()
        self.alpha_spin.setRange(0.1, 1.0)
        self.alpha_spin.setSingleStep(0.1)
        self.alpha_spin.setValue(0.7)
        self.layout.addRow("Point Transparency:", self.alpha_spin)
        
        self.width_spin = QDoubleSpinBox()
        self.width_spin.setRange(2.0, 30.0)
        self.width_spin.setValue(10.0)
        self.layout.addRow("Width (inches):", self.width_spin)
        
        self.height_spin = QDoubleSpinBox()
        self.height_spin.setRange(2.0, 30.0)
        self.height_spin.setValue(8.0)
        self.layout.addRow("Height (inches):", self.height_spin)
        
        self.transparent_chk = QCheckBox("Transparent Background")
        self.layout.addRow("", self.transparent_chk)
        
        self.btn_box = QHBoxLayout()
        self.btn_save = QPushButton("Save Image")
        self.btn_cancel = QPushButton("Cancel")
        self.btn_box.addWidget(self.btn_save)
        self.btn_box.addWidget(self.btn_cancel)
        self.layout.addRow(self.btn_box)
        
        self.btn_save.clicked.connect(self.accept)
        self.btn_cancel.clicked.connect(self.reject)


class ChemicalSpaceWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Chemical Space Visualization")
        self.resize(1200, 800)

        # Engines
        self.data_engine = DataEngine()
        self.feat_engine = FeaturizationEngine()
        self.analysis_engine = AnalysisEngine()

        # State
        self.X = None
        self.df_valid = None
        self.annot_labels = []
        self.cbar = None

        self.setup_ui()

    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        # Splitter
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)

        # Left Panel (Controls)
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setAlignment(Qt.AlignTop)

        # 1. Data Loading
        data_group = QGroupBox("1. Data Loading")
        data_layout = QVBoxLayout(data_group)
        
        self.btn_load_csv = QPushButton("Load CSV")
        self.btn_load_csv.clicked.connect(self.load_csv)
        data_layout.addWidget(self.btn_load_csv)
        
        self.lbl_stats = QLabel("No data loaded")
        data_layout.addWidget(self.lbl_stats)

        smiles_layout = QHBoxLayout()
        smiles_layout.addWidget(QLabel("SMILES Column:"))
        self.combo_smiles = QComboBox()
        smiles_layout.addWidget(self.combo_smiles)
        data_layout.addLayout(smiles_layout)

        left_layout.addWidget(data_group)

        # 2. Featurization
        feat_group = QGroupBox("2. Featurization")
        feat_layout = QVBoxLayout(feat_group)
        
        self.chk_fingerprints = QCheckBox("Use Fingerprints")
        self.chk_fingerprints.setChecked(True)
        feat_layout.addWidget(self.chk_fingerprints)
        
        self.chk_descriptors = QCheckBox("Use Descriptors")
        self.chk_descriptors.setChecked(True)
        feat_layout.addWidget(self.chk_descriptors)
        
        self.chk_normalize = QCheckBox("Normalize (StandardScaler)")
        self.chk_normalize.setChecked(True)
        feat_layout.addWidget(self.chk_normalize)

        left_layout.addWidget(feat_group)

        # 3. Dimensionality Reduction
        dim_group = QGroupBox("3. Embedding (PCA/t-SNE)")
        dim_layout = QVBoxLayout(dim_group)
        
        method_layout = QHBoxLayout()
        method_layout.addWidget(QLabel("Method:"))
        self.combo_method = QComboBox()
        self.combo_method.addItems(["t-SNE", "PCA"])
        method_layout.addWidget(self.combo_method)
        dim_layout.addLayout(method_layout)

        tsne_layout = QHBoxLayout()
        tsne_layout.addWidget(QLabel("t-SNE Perplexity:"))
        self.spin_perplexity = QSpinBox()
        self.spin_perplexity.setRange(5, 100)
        self.spin_perplexity.setValue(30)
        tsne_layout.addWidget(self.spin_perplexity)
        dim_layout.addLayout(tsne_layout)

        left_layout.addWidget(dim_group)

        # 4. Clustering & Overlays
        clus_group = QGroupBox("4. Clustering")
        clus_layout = QVBoxLayout(clus_group)
        
        clus_method_layout = QHBoxLayout()
        clus_method_layout.addWidget(QLabel("Algorithm:"))
        self.combo_cluster = QComboBox()
        self.combo_cluster.addItems(["None", "K-Means", "DBSCAN"])
        clus_method_layout.addWidget(self.combo_cluster)
        clus_layout.addLayout(clus_method_layout)
        
        self.chk_outliers = QCheckBox("Highlight Outliers (IsoForest)")
        clus_layout.addWidget(self.chk_outliers)

        left_layout.addWidget(clus_group)

        # 5. Appearance
        app_group = QGroupBox("5. Appearance")
        app_layout = QVBoxLayout(app_group)
        
        color_layout = QHBoxLayout()
        color_layout.addWidget(QLabel("Color by:"))
        self.combo_color = QComboBox()
        self.combo_color.addItem("Default / Cluster")
        color_layout.addWidget(self.combo_color)
        
        color_layout.addWidget(QLabel("Palette:"))
        self.combo_cmap = QComboBox()
        self.combo_cmap.addItems(["viridis", "plasma", "inferno", "magma", "cividis", "coolwarm", "Spectral", 'twilight', 'twilight_shifted', 'turbo', 'berlin', 'managua', 'vanimo', 'Blues', 'BrBG', 'BuGn', 'BuPu', 'CMRmap', 'GnBu', 'Greens', 'Greys', 'OrRd', 'Oranges', 'PRGn', 'PiYG', 'PuBu', 'PuBuGn', 'PuOr', 'PuRd', 'Purples', 'RdBu', 'RdGy', 'RdPu', 'RdYlBu', 'RdYlGn', 'Reds'])
        color_layout.addWidget(self.combo_cmap)
        app_layout.addLayout(color_layout)

        label_layout = QHBoxLayout()
        label_layout.addWidget(QLabel("Label by:"))
        self.combo_label = QComboBox()
        label_layout.addWidget(self.combo_label)
        app_layout.addLayout(label_layout)

        left_layout.addWidget(app_group)

        # Action Buttons
        self.btn_run = QPushButton("🚀 Run Pipeline")
        self.btn_run.setStyleSheet("background-color: #2b5c8f; color: white; font-weight: bold; padding: 10px;")
        self.btn_run.clicked.connect(self.run_pipeline)
        left_layout.addWidget(self.btn_run)

        self.btn_export = QPushButton("Export Results (CSV)")
        self.btn_export.clicked.connect(self.export_csv)
        left_layout.addWidget(self.btn_export)

        self.btn_export_img = QPushButton("Export High-Res Image")
        self.btn_export_img.clicked.connect(self.export_image)
        left_layout.addWidget(self.btn_export_img)

        left_layout.addStretch()

        # Right Panel (Matplotlib Canvas)
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        right_layout.addWidget(self.toolbar)
        right_layout.addWidget(self.canvas)

        # Hover annotation
        self.ax = self.figure.add_subplot(111)
        self.ax.set_title("Chemical Space")
        self.annot = self.ax.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="w", alpha=0.9),
                                      arrowprops=dict(arrowstyle="->"))
        self.annot.set_visible(False)
        self.canvas.mpl_connect("motion_notify_event", self.hover)
        self.canvas.mpl_connect("button_press_event", self.on_click)

        # Add to splitter
        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        splitter.setSizes([300, 900])

    def load_csv(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Select CSV", "", "CSV Files (*.csv)")
        if file_name:
            success, msg = self.data_engine.load_csv(file_name)
            if success:
                cols = self.data_engine.get_columns()
                self.combo_smiles.clear()
                self.combo_smiles.addItems(cols)
                
                # Auto-select SMILES column if found
                for i, col in enumerate(cols):
                    if 'smiles' in col.lower():
                        self.combo_smiles.setCurrentIndex(i)
                        break

                self.combo_label.clear()
                self.combo_label.addItems(cols)
                # Auto-select same column for labels
                self.combo_label.setCurrentIndex(self.combo_smiles.currentIndex())

                num_cols = self.data_engine.get_numeric_columns()
                self.combo_color.clear()
                self.combo_color.addItem("Default / Cluster")
                self.combo_color.addItems(num_cols)

                self.lbl_stats.setText(f"Total Rows: {self.data_engine.total_rows}")
            else:
                QMessageBox.warning(self, "Load Error", msg)

    def run_pipeline(self):
        if self.data_engine.df is None:
            QMessageBox.warning(self, "Warning", "Please load a CSV first.")
            return

        smiles_col = self.combo_smiles.currentText()
        if not smiles_col:
            return

        # 1. Cleanse Data
        self.data_engine.cleanse_data(smiles_col)
        self.lbl_stats.setText(f"Valid Compounds: {self.data_engine.valid_compounds}\nSkipped: {self.data_engine.skipped_rows}")

        # 2. Featurization (Async)
        self.btn_run.setEnabled(False)
        self.progress_dialog = QProgressDialog("Featurizing molecules...", "Cancel", 0, 100, self)
        self.progress_dialog.setWindowTitle("Processing Data")
        self.progress_dialog.setWindowModality(Qt.WindowModal)

        use_fp = self.chk_fingerprints.isChecked()
        use_desc = self.chk_descriptors.isChecked()
        norm = self.chk_normalize.isChecked()

        # Run feat engine
        self.worker = WorkerThread(self.feat_engine.featurize_dataset, 
                                   self.data_engine.df.copy(), 
                                   smiles_col, 
                                   use_fp, use_desc, norm)
        self.worker.progress.connect(self.progress_dialog.setValue)
        self.worker.finished.connect(self.on_featurization_complete)
        self.worker.start()

    def on_featurization_complete(self, success, msg, result):
        self.progress_dialog.close()
        if not success:
            QMessageBox.critical(self, "Error", f"Featurization failed: {msg}")
            self.btn_run.setEnabled(True)
            return

        X, df_valid = result
        if X is None or len(X) == 0:
            QMessageBox.critical(self, "Error", "No valid features generated.")
            self.btn_run.setEnabled(True)
            return

        self.X = X
        self.df_valid = df_valid

        # 3. Embedding (Async)
        self.progress_dialog = QProgressDialog("Running dimensionality reduction...", "Cancel", 0, 100, self)
        self.progress_dialog.setWindowTitle("Processing Data")
        self.progress_dialog.setRange(0, 0) # Indeterminate
        
        method = self.combo_method.currentText()
        perp = self.spin_perplexity.value()

        if method == "PCA":
            self.worker = WorkerThread(self.analysis_engine.run_pca, self.X, 2)
        else:
            self.worker = WorkerThread(self.analysis_engine.run_tsne, self.X, perp)

        self.worker.finished.connect(self.on_embedding_complete)
        self.worker.start()

    def on_embedding_complete(self, success, msg, _):
        self.progress_dialog.close()
        if not success:
            QMessageBox.critical(self, "Error", f"Embedding failed: {msg}")
            self.btn_run.setEnabled(True)
            return

        # 4. Clustering (Sync, usually fast)
        clus_method = self.combo_cluster.currentText()
        if clus_method == "K-Means":
            self.analysis_engine.run_kmeans()
        elif clus_method == "DBSCAN":
            self.analysis_engine.run_dbscan()
        else:
            self.analysis_engine.clusters = None

        if self.chk_outliers.isChecked():
            self.analysis_engine.detect_outliers(self.X)
        else:
            self.analysis_engine.outliers = None

        self.update_plot()
        self.btn_run.setEnabled(True)

    def update_plot(self):
        # Clear the entire figure to prevent grid space stacking (squishing bug)
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.cbar = None
        self.annot_labels = []

        if self.analysis_engine.embedding is None:
            self.canvas.draw()
            return

        emb = self.analysis_engine.embedding
        x, y = emb[:, 0], emb[:, 1]

        # Determine colors
        color_col = self.combo_color.currentText()
        colors = None
        cmap = None

        if color_col != "Default / Cluster" and color_col in self.df_valid.columns:
            # Continuous color
            colors = self.df_valid[color_col].values
            cmap = self.combo_cmap.currentText()
        elif self.analysis_engine.clusters is not None:
            # Discrete cluster color
            colors = self.analysis_engine.clusters
            cmap = 'tab10'
        else:
            colors = 'blue'

        # Scatter
        self.scatter = self.ax.scatter(x, y, c=colors, cmap=cmap, alpha=0.7, edgecolors='w', s=50)

        if cmap != 'tab10' and color_col != "Default / Cluster":
            self.cbar = self.figure.colorbar(self.scatter, ax=self.ax, fraction=0.03, pad=0.02, aspect=30)
            self.cbar.set_label(color_col)

        # Outliers
        if self.analysis_engine.outliers is not None:
            outlier_idx = self.analysis_engine.outliers
            self.ax.scatter(x[outlier_idx], y[outlier_idx], facecolors='none', edgecolors='r', s=100, label='Outliers')
            self.ax.legend()

        self.ax.set_title(f"Chemical Space ({self.combo_method.currentText()})")
        self.ax.set_xlabel("Dimension 1")
        self.ax.set_ylabel("Dimension 2")

        # Re-add annotation object
        self.annot = self.ax.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="w", alpha=0.9),
                                      arrowprops=dict(arrowstyle="->"))
        self.annot.set_visible(False)

        self.figure.tight_layout()
        self.canvas.draw()

    def update_annot(self, ind):
        pos = self.scatter.get_offsets()[ind["ind"][0]]
        self.annot.xy = pos
        
        # Adjust xytext to prevent overlapping colorbar on the right
        x = pos[0]
        xlim = self.ax.get_xlim()
        x_range = xlim[1] - xlim[0]
        if x > xlim[0] + 0.6 * x_range:
            self.annot.set_position((-20, 20))
            self.annot.set_horizontalalignment('right')
        else:
            self.annot.set_position((20, 20))
            self.annot.set_horizontalalignment('left')
        
        idx = ind["ind"][0]
        label_col = self.combo_label.currentText()
        if not label_col:
            label_col = self.combo_smiles.currentText()
            
        label_val = self.df_valid.iloc[idx][label_col]
        
        text = f"Idx: {idx}\n{label_col}: {label_val}"
        
        if self.analysis_engine.clusters is not None:
            text += f"\nCluster: {self.analysis_engine.clusters[idx]}"
            
        color_col = self.combo_color.currentText()
        if color_col != "Default / Cluster" and color_col in self.df_valid.columns:
            val = self.df_valid.iloc[idx][color_col]
            text += f"\n{color_col}: {val}"

        self.annot.set_text(text)
        self.annot.get_bbox_patch().set_alpha(0.9)

    def hover(self, event):
        vis = self.annot.get_visible()
        if event.inaxes == self.ax and hasattr(self, 'scatter'):
            cont, ind = self.scatter.contains(event)
            if cont:
                self.update_annot(ind)
                self.annot.set_visible(True)
                self.canvas.draw_idle()
            else:
                if vis:
                    self.annot.set_visible(False)
                    self.canvas.draw_idle()

    def on_click(self, event):
        if event.inaxes == self.ax and hasattr(self, 'scatter'):
            cont, ind = self.scatter.contains(event)
            if cont:
                idx = ind["ind"][0]
                pos = self.scatter.get_offsets()[idx]
                
                label_col = self.combo_label.currentText()
                if not label_col:
                    label_col = self.combo_smiles.currentText()
                    
                label_val = str(self.df_valid.iloc[idx][label_col])
                
                # Add persistent, draggable annotation with a connecting line
                annot = self.ax.annotate(
                    label_val,
                    xy=pos,
                    xytext=(20, 20),
                    textcoords="offset points",
                    fontsize=9,
                    fontweight='bold',
                    color='black',
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9),
                    arrowprops=dict(arrowstyle="-", color="gray", alpha=0.8, connectionstyle="arc3")
                )
                
                # Try setting it as draggable (works dynamically in matplotlib 3+)
                if hasattr(annot, 'set_draggable'):
                    annot.set_draggable(True)
                elif hasattr(annot, 'draggable'):
                    annot.draggable(True)
                    
                self.annot_labels.append(annot)
                self.canvas.draw_idle()

    def export_csv(self):
        if self.df_valid is None or self.analysis_engine.embedding is None:
            QMessageBox.warning(self, "Warning", "No data to export.")
            return

        file_name, _ = QFileDialog.getSaveFileName(self, "Export CSV", "chemical_space.csv", "CSV Files (*.csv)")
        if file_name:
            out_df = self.df_valid.copy()
            out_df['X_Coord'] = self.analysis_engine.embedding[:, 0]
            out_df['Y_Coord'] = self.analysis_engine.embedding[:, 1]
            
            if self.analysis_engine.clusters is not None:
                out_df['Cluster'] = self.analysis_engine.clusters
                
            if self.analysis_engine.outliers is not None:
                out_df['Is_Outlier'] = self.analysis_engine.outliers
                
            out_df.to_csv(file_name, index=False)
            QMessageBox.information(self, "Success", f"Data exported to {file_name}")

    def export_image(self):
        if self.analysis_engine.embedding is None:
            QMessageBox.warning(self, "Warning", "No plot to export.")
            return
            
        dialog = ExportImageDialog(self)
        if dialog.exec() == QDialog.Accepted:
            file_name, _ = QFileDialog.getSaveFileName(self, "Save Image", "chemical_space.png", "PNG Images (*.png);;PDF Files (*.pdf);;SVG Files (*.svg)")
            if not file_name:
                return
                
            dpi = dialog.dpi_spin.value()
            font_size = dialog.font_spin.value()
            point_size = dialog.point_spin.value()
            alpha = dialog.alpha_spin.value()
            w = dialog.width_spin.value()
            h = dialog.height_spin.value()
            transparent = dialog.transparent_chk.isChecked()
            
            orig_size = self.figure.get_size_inches()
            orig_scatter_sizes = self.scatter.get_sizes() if hasattr(self, 'scatter') else [50]
            orig_scatter_alpha = self.scatter.get_alpha() if hasattr(self, 'scatter') else 0.7
            
            try:
                self.figure.set_size_inches(w, h)
                for text in self.annot_labels:
                    text.set_fontsize(font_size)
                self.ax.title.set_fontsize(font_size + 4)
                self.ax.xaxis.label.set_fontsize(font_size + 2)
                self.ax.yaxis.label.set_fontsize(font_size + 2)
                self.ax.tick_params(axis='both', which='major', labelsize=font_size)
                
                if hasattr(self, 'scatter'):
                    self.scatter.set_sizes([point_size] * len(self.scatter.get_offsets()))
                    self.scatter.set_alpha(alpha)
                    
                if self.cbar is not None:
                    self.cbar.ax.tick_params(labelsize=font_size)
                    self.cbar.set_label(self.cbar.ax.get_ylabel(), size=font_size+2)
                    
                self.figure.savefig(file_name, dpi=dpi, transparent=transparent, bbox_inches='tight')
                QMessageBox.information(self, "Success", f"Image successfully exported to\n{file_name}")
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export image: {e}")
            finally:
                self.figure.set_size_inches(orig_size)
                for text in self.annot_labels:
                    text.set_fontsize(9) 
                self.ax.title.set_fontsize(12)
                self.ax.xaxis.label.set_fontsize(10)
                self.ax.yaxis.label.set_fontsize(10)
                self.ax.tick_params(axis='both', which='major', labelsize=10)
                
                if hasattr(self, 'scatter'):
                    self.scatter.set_sizes(orig_scatter_sizes)
                    self.scatter.set_alpha(orig_scatter_alpha)
                    
                if self.cbar is not None:
                    self.cbar.ax.tick_params(labelsize=10)
                    self.cbar.set_label(self.cbar.ax.get_ylabel(), size=10)
                    
                self.canvas.draw()
