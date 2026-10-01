import os
import multiprocessing
import pandas as pd

from src.shared.qt_compat import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QFileDialog, QComboBox, QTextEdit, QMessageBox,
    QDoubleSpinBox, QSpinBox, QGroupBox, QFormLayout, QSplitter,
    QButtonGroup, QRadioButton, QListWidget, QListWidgetItem, QAbstractItemView,
    QDialog, QDialogButtonBox, QScrollArea, QSizePolicy, QThread, Signal, Qt, QPixmap
)

from .split_engine import DataSplitEngine
from .algorithms import ALGOS_NEED_FEATURES
from .utils.descriptors import DescriptorMode
from .utils.fingerprints import ALL_FP_TYPES

class SlowCalcWarningDialog(QDialog):
    """Consent dialog shown before slow descriptor calculation."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("⚠ Time Warning")
        self.setFixedWidth(440)
        lay = QVBoxLayout(self)
        msg = QLabel(
            "<b>Descriptor calculation with PyDes may take a very long time</b> "
            "for large datasets (minutes to hours).<br><br>"
            "Make sure you have selected the correct SMILES column and that "
            "molecules are valid before proceeding.<br><br>"
            "Do you want to continue?"
        )
        msg.setWordWrap(True)
        lay.addWidget(msg)
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Yes | QDialogButtonBox.StandardButton.No)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        lay.addWidget(btns)

class SplitWorker(QThread):
    log_signal      = Signal(str)
    finished_signal = Signal(bool, str, dict)
    plot_signal     = Signal(str)

    def __init__(self, engine, input_file, output_file, smiles_col, name_col,
                 target_col, ratio, algorithm, n_jobs,
                 desc_mode, csv_desc_cols, fp_type, metric="euclidean"):
        super().__init__()
        self.engine        = engine
        self.input_file    = input_file
        self.output_file   = output_file
        self.smiles_col    = smiles_col
        self.name_col      = name_col if name_col != "(None)" else None
        self.target_col    = target_col if target_col != "(None)" else None
        self.ratio         = ratio
        self.algorithm     = algorithm
        self.n_jobs        = n_jobs
        self.desc_mode     = desc_mode          # DescriptorMode enum
        self.csv_desc_cols = csv_desc_cols      # list[str] or []
        self.fp_type       = fp_type            # fingerprint type string
        self.metric        = metric

    def run(self):
        try:
            self.log_signal.emit(f"Loading: {os.path.basename(self.input_file)}")
            self.engine.load_csv(self.input_file)
            self.log_signal.emit(f"Loaded {len(self.engine.df)} records.")

            needs_features = self.algorithm in ALGOS_NEED_FEATURES
            if needs_features:
                mode_label = {
                    DescriptorMode.USE_CSV:      "CSV descriptors",
                    DescriptorMode.FINGERPRINTS:  f"fingerprints ({self.fp_type})",
                    DescriptorMode.CALCULATE:     "PyDes full descriptors",
                }[self.desc_mode]
                self.log_signal.emit(f"Feature source: {mode_label}  |  Workers: {self.n_jobs}")

            self.log_signal.emit(f"Running {self.algorithm} splitter…")
            result = self.engine.split(
                algorithm=self.algorithm,
                target_ratio=self.ratio,
                smiles_col=self.smiles_col,
                id_col=self.name_col,
                target_col=self.target_col,
                desc_mode=self.desc_mode,
                csv_desc_cols=self.csv_desc_cols,
                fp_type=self.fp_type,
                metric=self.metric,
                n_jobs=self.n_jobs
            )

            self.log_signal.emit(f"Saving → {os.path.basename(self.output_file)}")
            result.annotated_df.to_csv(self.output_file, index=False)

            if needs_features and self.engine._desc_df is not None and not self.engine._desc_df.empty:
                self.log_signal.emit("Generating PCA plot…")
                plot_path = os.path.join(os.path.dirname(self.output_file), "chemical_space_split.png")
                if self.engine.plot_pca(result.train_indices, result.test_indices, plot_path):
                    self.plot_signal.emit(plot_path)

            self.finished_signal.emit(True, "Dataset successfully split and saved!", result.stats)

        except Exception as e:
            import traceback
            self.log_signal.emit(f"ERROR: {e}\n{traceback.format_exc()}")
            self.finished_signal.emit(False, str(e), {})


class DataSplittingWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("QSAR Data Splitting Tool")
        self.engine = DataSplitEngine()
        
        # Don't apply hardcoded dark theme; PyChem-Pro handles app-wide themes
        
        self.input_path   = ""
        self.output_path  = ""

        self._setup_ui()

    def _setup_ui(self):
        cw = QWidget()
        self.setCentralWidget(cw)
        main = QHBoxLayout(cw)
        main.setContentsMargins(5, 5, 5, 5)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)
        main.addWidget(splitter)

        # ── LEFT PANEL ─────────────────────────────
        left_content = QWidget()
        left_content.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        ll = QVBoxLayout(left_content)
        ll.setContentsMargins(10, 10, 5, 10)
        ll.setSpacing(6)

        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setWidget(left_content)

        # 1. File import
        g1 = QGroupBox("1. Dataset Import")
        g1l = QVBoxLayout(g1)
        self.btn_input = QPushButton("Select Input CSV File")
        self.btn_input.clicked.connect(self._select_input)
        self.lbl_input = QLabel("No file selected.")
        self.lbl_input.setWordWrap(True)
        g1l.addWidget(self.btn_input)
        g1l.addWidget(self.lbl_input)
        ll.addWidget(g1)

        # 2. Column mapping
        g2 = QGroupBox("2. Variable Mapping")
        g2l = QFormLayout(g2)
        self.combo_smiles  = QComboBox()
        self.combo_name    = QComboBox()
        self.combo_target  = QComboBox()
        g2l.addRow("SMILES Column:",          self.combo_smiles)
        g2l.addRow("Molecule ID Column:",      self.combo_name)
        g2l.addRow("Target/Activity Column:",  self.combo_target)
        ll.addWidget(g2)

        # 4. Feature source
        g3 = QGroupBox("4. Feature Source (for distance-based algorithms)")
        g3l = QVBoxLayout(g3)

        self._feat_grp = QButtonGroup(self)
        self.rb_csv   = QRadioButton("Use descriptor columns from CSV")
        self.rb_fp    = QRadioButton("Calculate Fingerprints (fast, no 3-D needed)")
        self.rb_calc  = QRadioButton("Calculate Full Descriptors via PyDes  ⚠ slow")
        self.rb_fp.setChecked(True)
        for rb in (self.rb_csv, self.rb_fp, self.rb_calc):
            self._feat_grp.addButton(rb)
            g3l.addWidget(rb)

        self.rb_csv.toggled.connect(self._on_feat_mode_change)
        self.rb_fp.toggled.connect(self._on_feat_mode_change)

        self.lbl_csv_cols  = QLabel("Select descriptor columns (Ctrl+click for multi):")
        self.lst_csv_cols  = QListWidget()
        self.lst_csv_cols.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.lst_csv_cols.setMaximumHeight(100)
        self.lbl_csv_cols.setVisible(False)
        self.lst_csv_cols.setVisible(False)
        g3l.addWidget(self.lbl_csv_cols)
        g3l.addWidget(self.lst_csv_cols)

        self.lbl_desc_detected = QLabel("")
        self.lbl_desc_detected.setWordWrap(True)
        g3l.addWidget(self.lbl_desc_detected)

        fp_row = QHBoxLayout()
        fp_row.addWidget(QLabel("Fingerprint type:"))
        self.combo_fp = QComboBox()
        self.combo_fp.addItems(ALL_FP_TYPES)
        fp_row.addWidget(self.combo_fp)
        self.fp_row_widget = QWidget()
        self.fp_row_widget.setLayout(fp_row)
        g3l.addWidget(self.fp_row_widget)

        # 3. Splitting config
        g4 = QGroupBox("3. Splitting Configuration")
        g4l = QFormLayout(g4)
        self.spin_ratio = QDoubleSpinBox()
        self.spin_ratio.setRange(0.05, 0.95)
        self.spin_ratio.setValue(0.80)
        self.spin_ratio.setSingleStep(0.05)
        self.combo_algo = QComboBox()
        self.combo_algo.addItems(self.engine.available_algorithms())
        self.combo_algo.currentTextChanged.connect(self._on_algo_change)
        
        self.lbl_metric = QLabel("Distance Metric:")
        self.combo_metric = QComboBox()
        self.combo_metric.addItems(["Euclidean", "Manhattan", "Cosine", "Tanimoto / Jaccard"])
        
        self.spin_cores = QSpinBox()
        mc = multiprocessing.cpu_count()
        self.spin_cores.setRange(1, mc)
        self.spin_cores.setValue(max(1, mc // 2))
        g4l.addRow("Train Ratio:",        self.spin_ratio)
        g4l.addRow("Algorithm:",          self.combo_algo)
        g4l.addRow(self.lbl_metric,       self.combo_metric)
        g4l.addRow("Processor Count:",    self.spin_cores)
        ll.addWidget(g4)
        ll.addWidget(g3)

        # 5. Output
        g5 = QGroupBox("5. Export Destination")
        g5l = QVBoxLayout(g5)
        self.btn_output = QPushButton("Set Destination Path")
        self.btn_output.clicked.connect(self._select_output)
        self.lbl_output = QLabel("No export path set.")
        self.lbl_output.setWordWrap(True)
        g5l.addWidget(self.btn_output)
        g5l.addWidget(self.lbl_output)
        ll.addWidget(g5)

        btn_layout = QHBoxLayout()
        self.btn_run = QPushButton("▶  Run Splitting")
        self.btn_run.clicked.connect(self._run)
        
        self.btn_stop = QPushButton("⏹  Stop Splitting")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self._stop)
        
        btn_layout.addWidget(self.btn_run)
        btn_layout.addWidget(self.btn_stop)
        ll.addLayout(btn_layout)

        ll.addWidget(QLabel("Execution Log:"))
        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)
        self.txt_log.setMinimumHeight(120)
        self.txt_log.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        ll.addWidget(self.txt_log)
        ll.addStretch(0)

        splitter.addWidget(left_scroll)

        # ── RIGHT PANEL ───────────────────────────────────────────────────
        right = QWidget()
        rl    = QVBoxLayout(right)
        rl.setContentsMargins(5, 10, 10, 10)
        rl.setSpacing(8)

        gs = QGroupBox("Split Statistics")
        gsl = QFormLayout(gs)
        self.lbl_total = QLabel("-")
        self.lbl_train = QLabel("-")
        self.lbl_test  = QLabel("-")
        self.lbl_ks    = QLabel("-")
        gsl.addRow("Total Compounds:", self.lbl_total)
        gsl.addRow("Train Set:",        self.lbl_train)
        gsl.addRow("Test Set:",         self.lbl_test)
        gsl.addRow("KS Distribution:",  self.lbl_ks)
        rl.addWidget(gs)

        gp = QGroupBox("Chemical Space Projection (PCA)")
        gpl = QVBoxLayout(gp)
        self.lbl_plot = QLabel("Visualization will appear here after splitting.")
        self.lbl_plot.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.lbl_plot.setMinimumSize(400, 320)
        gpl.addWidget(self.lbl_plot)
        rl.addWidget(gp)

        splitter.addWidget(right)
        splitter.setStretchFactor(0, 40)
        splitter.setStretchFactor(1, 60)

        self._on_feat_mode_change()
        self._on_algo_change(self.combo_algo.currentText())

    def _on_feat_mode_change(self):
        csv_active = self.rb_csv.isChecked()
        fp_active  = self.rb_fp.isChecked()
        self.lbl_csv_cols.setVisible(csv_active)
        self.lst_csv_cols.setVisible(csv_active)
        self.fp_row_widget.setVisible(fp_active)
        
        if fp_active:
            self.combo_metric.setCurrentText("Tanimoto / Jaccard")
        else:
            self.combo_metric.setCurrentText("Euclidean")

    def _on_algo_change(self, algo):
        needs = algo in ALGOS_NEED_FEATURES
        self.rb_csv.setEnabled(needs)
        self.rb_fp.setEnabled(needs)
        self.rb_calc.setEnabled(needs)
        self.fp_row_widget.setEnabled(needs)
        self.lst_csv_cols.setEnabled(needs)
        
        is_distance_algo = algo in ["CADEX (Kennard-Stone)", "Sphere Exclusion", "Boruta-based", "Duplex"]
        self.lbl_metric.setVisible(is_distance_algo)
        self.combo_metric.setVisible(is_distance_algo)
        
        if algo == "Boruta-based":
            self.txt_log.append("[Info] Boruta requires a target/activity column.")

    def _select_input(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open CSV", "", "CSV Files (*.csv)")
        if not path:
            return
        
        if self.engine.load_csv(path):
            self.input_path = path
            self.lbl_input.setText(path)
            
            cols = self.engine.original_columns
            num_cols = self.engine.get_numeric_columns()
            
            for cb in (self.combo_smiles, self.combo_name, self.combo_target):
                cb.clear()
            self.combo_name.addItem("(None)")
            self.combo_target.addItem("(None)")
            self.combo_smiles.addItems(cols)
            self.combo_name.addItems(cols)
            self.combo_target.addItems(cols)
            
            guessed = self.engine.guess_columns()
            if guessed["smiles"]: self.combo_smiles.setCurrentText(guessed["smiles"])
            if guessed["id"]: self.combo_name.setCurrentText(guessed["id"])
            if guessed["target"]: self.combo_target.setCurrentText(guessed["target"])

            self.lst_csv_cols.clear()
            for c in num_cols:
                self.lst_csv_cols.addItem(QListWidgetItem(c))

            non_meta = [c for c in num_cols if not any(k in c.lower() for k in ("smiles","id","name","title","activity","target","pic50","ic50","pki"))]
            if len(non_meta) >= 10:
                self.lbl_desc_detected.setText(f"✔ {len(non_meta)} numeric columns detected — 'Use CSV descriptors' is recommended.")
                self.rb_csv.setChecked(True)
                for i in range(self.lst_csv_cols.count()):
                    item = self.lst_csv_cols.item(i)
                    if item.text() in non_meta:
                        item.setSelected(True)
            else:
                self.lbl_desc_detected.setText("")
                self.rb_fp.setChecked(True)

            self.txt_log.append(f"[Info] Loaded {len(cols)} columns from {os.path.basename(path)}")
            base = os.path.splitext(os.path.basename(path))[0]
            self.output_path = os.path.join(os.path.dirname(path), f"{base}_split.csv")
            self.lbl_output.setText(self.output_path)
            self._on_feat_mode_change()

    def _select_output(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", self.output_path or "split.csv", "CSV Files (*.csv)")
        if path:
            self.output_path = path
            self.lbl_output.setText(path)

    def _run(self):
        if not self.input_path:
            QMessageBox.warning(self, "Missing File", "Please select an input CSV."); return
        if not self.output_path:
            QMessageBox.warning(self, "Missing Path", "Please set a destination path."); return

        algo   = self.combo_algo.currentText()
        target = self.combo_target.currentText()

        if algo == "Boruta-based" and target == "(None)":
            QMessageBox.warning(self, "Target Required", "Boruta requires a target/activity column."); return

        if self.rb_csv.isChecked():
            mode = DescriptorMode.USE_CSV
            csv_cols = [self.lst_csv_cols.item(i).text() for i in range(self.lst_csv_cols.count()) if self.lst_csv_cols.item(i).isSelected()]
            if algo in ALGOS_NEED_FEATURES and not csv_cols:
                QMessageBox.warning(self, "No Columns Selected", "Please select at least one descriptor column from the list."); return
        elif self.rb_fp.isChecked():
            mode     = DescriptorMode.FINGERPRINTS
            csv_cols = []
        else:
            dlg = SlowCalcWarningDialog(self)
            if dlg.exec() != QDialog.DialogCode.Accepted:
                return
            mode     = DescriptorMode.CALCULATE
            csv_cols = []

        fp_type = self.combo_fp.currentText()

        for w in (self.btn_run, self.btn_input, self.btn_output):
            w.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.txt_log.clear()
        for lbl in (self.lbl_total, self.lbl_train, self.lbl_test, self.lbl_ks):
            lbl.setText("Calculating…")
        self.lbl_plot.clear()
        self.lbl_plot.setText("Generating projection…")

        metric_map = {"Euclidean": "euclidean", "Manhattan": "cityblock", "Cosine": "cosine", "Tanimoto / Jaccard": "jaccard"}
        selected_metric = metric_map.get(self.combo_metric.currentText(), "euclidean")

        self.worker = SplitWorker(
            engine=self.engine,
            input_file=self.input_path,
            output_file=self.output_path,
            smiles_col=self.combo_smiles.currentText(),
            name_col=self.combo_name.currentText(),
            target_col=target,
            ratio=self.spin_ratio.value(),
            algorithm=algo,
            n_jobs=self.spin_cores.value(),
            desc_mode=mode,
            csv_desc_cols=csv_cols,
            fp_type=fp_type,
            metric=selected_metric,
        )
        self.worker.log_signal.connect(self.txt_log.append)
        self.worker.plot_signal.connect(self._show_plot)
        self.worker.finished_signal.connect(self._on_finished)
        self.worker.start()

    def _show_plot(self, path):
        px = QPixmap(path)
        self.lbl_plot.setPixmap(px.scaled(self.lbl_plot.size(), Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation))

    def _on_finished(self, success, message, stats):
        for w in (self.btn_run, self.btn_input, self.btn_output):
            w.setEnabled(True)
        self.btn_stop.setEnabled(False)

        if success:
            dupes = stats.get("duplicates_dropped", 0)
            if dupes > 0:
                message += f"\n\nNote: {dupes} duplicate SMILES were removed from the dataset."
                self.txt_log.append(f"[Info] {dupes} duplicate SMILES were removed.")
                
            QMessageBox.information(self, "Done", message)
            self.lbl_total.setText(str(stats.get("total_count", "-")))
            self.lbl_train.setText(f"{stats.get('train_count',0)} ({stats.get('train_pct',0):.1f}%)")
            self.lbl_test.setText(f"{stats.get('test_count',0)} ({stats.get('test_pct',0):.1f}%)")
            if "ks_pvalue" in stats:
                p = stats["ks_pvalue"]
                verdict = "Similar ✔" if p > 0.05 else "Differ ✘"
                self.lbl_ks.setText(f"p={p:.4f}  {verdict}")
            else:
                self.lbl_ks.setText("N/A (no numeric target)")
            self.txt_log.append("\n[Done] Splitting completed successfully.")
        else:
            QMessageBox.critical(self, "Error", f"Failed:\n{message}")
            for lbl in (self.lbl_total, self.lbl_train, self.lbl_test, self.lbl_ks):
                lbl.setText("Failed")
            self.lbl_plot.setText("Visualization unavailable.")


    def _stop(self):
        if hasattr(self, 'worker') and self.worker.isRunning():
            self.txt_log.append("\n[Info] Stopping split process, please wait...")
            
            from .utils.descriptors import ACTIVE_POOLS as DESCRIPTOR_POOLS
            from .utils.fingerprints import ACTIVE_POOLS as FINGERPRINT_POOLS
            
            for pool in list(DESCRIPTOR_POOLS):
                try: pool.terminate()
                except Exception: pass
            DESCRIPTOR_POOLS.clear()
            
            for pool in list(FINGERPRINT_POOLS):
                try: pool.terminate()
                except Exception: pass
            FINGERPRINT_POOLS.clear()
            
            self.worker.terminate()
            self.worker.wait()
            
            self.txt_log.append("[Warning] Splitting was stopped by the user.")
            
            for w in (self.btn_run, self.btn_input, self.btn_output):
                w.setEnabled(True)
            self.btn_stop.setEnabled(False)
            self.lbl_plot.clear()
            self.lbl_plot.setText("Splitting stopped.")
