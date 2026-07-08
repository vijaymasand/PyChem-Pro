"""
QSAR Dataset Curator Plugin (Universal)

Works with any CSV format (ChEMBL, BindingDB, PubChem, custom, etc.).
User selects SMILES column, activity column, extra columns to keep,
and configures all curation steps via the GUI.
"""

import math
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QFileDialog, QComboBox, QTextEdit, QMessageBox, QCheckBox,
    QGroupBox, QScrollArea, QListWidget, QListWidgetItem,
    QTabWidget, QTableWidget, QTableWidgetItem, QHeaderView,
    QSplitter, QLineEdit, QDoubleSpinBox, QFormLayout,
    QThread, Signal, Qt
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType


# ---------------------------------------------------------------------------
# Background worker
# ---------------------------------------------------------------------------
class CurationWorker(QThread):
    log_signal      = Signal(str)
    progress_signal = Signal(int, int)   # current, total
    finished_signal = Signal(bool, str, object)  # success, msg, dataframe

    def __init__(self, config: dict):
        super().__init__()
        self.cfg = config
        self.lfc = rdMolStandardize.LargestFragmentChooser()
        self.uc  = rdMolStandardize.Uncharger()
        self.allowed_atoms = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}

    # ------------------------------------------------------------------
    def _read(self, path):
        try:
            return pd.read_csv(path)
        except Exception:
            return pd.read_csv(path, encoding='latin1')

    def _curate_mol(self, smi):
        try:
            mol = Chem.MolFromSmiles(str(smi))
            if mol is None:
                return None
            if self.cfg.get('filter_organometallics', True):
                if any(a.GetAtomicNum() not in self.allowed_atoms for a in mol.GetAtoms()):
                    return None
            if self.cfg.get('remove_salts', True):
                mol = self.lfc.choose(mol)
            if self.cfg.get('neutralize', True):
                mol = self.uc.uncharge(mol)
            return Chem.MolToSmiles(mol)
        except Exception:
            return None

    # ------------------------------------------------------------------
    def run(self):
        try:
            cfg = self.cfg
            smiles_col   = cfg['smiles_col']
            activity_col = cfg['activity_col']
            extra_cols   = cfg.get('extra_cols', [])
            strategy     = cfg.get('duplicate_strategy', 'Mean (Average)')
            out_path     = cfg['output_file']

            self.log_signal.emit("📂  Loading dataset…")
            df = self._read(cfg['source_file'])
            self.log_signal.emit(f"   Loaded {len(df):,} rows, {len(df.columns)} columns.")

            # ── Validate required columns ──────────────────────────────
            for col in [smiles_col, activity_col]:
                if col not in df.columns:
                    self.finished_signal.emit(False, f"Column '{col}' not found in file.", None)
                    return

            # ── Optional ChEMBL-style type filter ─────────────────────
            type_filter = cfg.get('activity_type_filter', '').strip()
            if type_filter and 'Standard Type' in df.columns:
                types = [t.strip() for t in type_filter.split(',')]
                df = df[df['Standard Type'].isin(types)]
                self.log_signal.emit(f"   After activity-type filter ({type_filter}): {len(df):,} rows.")

            # ── Optional unit filter ───────────────────────────────────
            unit_filter = cfg.get('unit_filter', '').strip()
            if unit_filter and 'Standard Units' in df.columns:
                df = df[df['Standard Units'] == unit_filter]
                self.log_signal.emit(f"   After unit filter ({unit_filter}): {len(df):,} rows.")

            # ── Optional relation filter ───────────────────────────────
            if cfg.get('filter_relations', False) and 'Standard Relation' in df.columns:
                df['Standard Relation'] = df['Standard Relation'].astype(str).str.replace("'", "").str.strip()
                df = df[df['Standard Relation'] == '=']
                self.log_signal.emit(f"   After relation filter (=): {len(df):,} rows.")

            # ── Drop rows missing SMILES or activity ──────────────────
            df = df.dropna(subset=[smiles_col, activity_col])
            df[activity_col] = pd.to_numeric(df[activity_col], errors='coerce')
            df = df.dropna(subset=[activity_col])
            self.log_signal.emit(f"   After dropping missing values: {len(df):,} rows.")

            # ── Optional activity value range filter ──────────────────
            act_min = cfg.get('activity_min')
            act_max = cfg.get('activity_max')
            if act_min is not None:
                df = df[df[activity_col] >= act_min]
            if act_max is not None:
                df = df[df[activity_col] <= act_max]
            if act_min is not None or act_max is not None:
                self.log_signal.emit(f"   After activity range filter: {len(df):,} rows.")

            # ── RDKit structural curation ──────────────────────────────
            self.log_signal.emit("🔬  Applying RDKit structural curation…")
            total = len(df)
            canonical = []
            for i, smi in enumerate(df[smiles_col]):
                canonical.append(self._curate_mol(smi))
                if i % 200 == 0:
                    self.progress_signal.emit(i, total)
            self.progress_signal.emit(total, total)
            df = df.copy()
            df['Canonical_SMILES'] = canonical
            df = df.dropna(subset=['Canonical_SMILES'])
            self.log_signal.emit(f"   After structural curation: {len(df):,} rows.")

            # ── pIC50 / pKi conversion ─────────────────────────────────
            if cfg.get('convert_pic50', False):
                unit = cfg.get('act_unit', 'nM')
                factor = {'nM': 1e-9, 'uM': 1e-6, 'mM': 1e-3, 'M': 1.0}.get(unit, 1e-9)
                df['pActivity'] = df[activity_col].apply(
                    lambda v: -math.log10(v * factor) if v > 0 else None
                )
                df = df.dropna(subset=['pActivity'])
                activity_col_out = 'pActivity'
                self.log_signal.emit(f"   Converted to pActivity ({unit} -> pIC50/pKi): {len(df):,} rows.")
            else:
                activity_col_out = activity_col

            # ── Duplicate handling ────────────────────────────────────
            self.log_signal.emit(f"🔁  Handling duplicates via '{strategy}'…")
            agg_map = {'Mean (Average)': 'mean', 'Median': 'median', 'Lowest (Most Potent)': 'min', 'Highest': 'max'}
            agg_func = agg_map.get(strategy, 'mean')

            agg_dict = {activity_col_out: agg_func}
            valid_extra = [c for c in extra_cols if c in df.columns and c != 'Canonical_SMILES']
            for c in valid_extra:
                agg_dict[c] = 'first'

            curated_df = df.groupby('Canonical_SMILES', as_index=False).agg(agg_dict)
            self.log_signal.emit(f"   Final dataset: {len(curated_df):,} unique molecules.")

            # ── Save ──────────────────────────────────────────────────
            curated_df.to_csv(out_path, index=False)
            self.log_signal.emit(f"✅  Saved → {out_path}")
            self.finished_signal.emit(True, "Curation completed successfully!", curated_df)

        except Exception as e:
            self.finished_signal.emit(False, str(e), None)


# ---------------------------------------------------------------------------
# Main Widget
# ---------------------------------------------------------------------------
class QsarCuratorWidget(PluginWidget):
    def __init__(self, plugin: 'QsarCuratorPlugin'):
        super().__init__(plugin)
        self.source_path = ""
        self.output_path = ""
        self.worker      = None
        self._df_columns = []
        self.setup_ui()

    # ------------------------------------------------------------------
    def setup_ui(self):
        self.widget = QWidget()
        root = QVBoxLayout(self.widget)
        root.setSpacing(6)

        # ── Header ────────────────────────────────────────────────────
        title = QLabel("⚗️  QSAR Dataset Curator  —  Universal")
        title.setStyleSheet("font-size:15px;font-weight:bold;color:#1a73e8;margin-bottom:2px;")
        sub = QLabel("Supports ChEMBL · BindingDB · PubChem · Any CSV format")
        sub.setStyleSheet("color:#666;font-size:11px;margin-bottom:6px;")
        root.addWidget(title)
        root.addWidget(sub)

        # ── Tabs ──────────────────────────────────────────────────────
        tabs = QTabWidget()
        tabs.addTab(self._build_data_tab(),    "📁 Data")
        tabs.addTab(self._build_columns_tab(), "🗂 Columns")
        tabs.addTab(self._build_filters_tab(), "🔧 Filters")
        tabs.addTab(self._build_options_tab(), "⚙️ Options")
        tabs.addTab(self._build_preview_tab(), "📊 Preview")
        root.addWidget(tabs)

        # ── Run + progress area ───────────────────────────────────────
        self.btn_run = QPushButton("🚀  Run Curation")
        self.btn_run.setStyleSheet(
            "background:#1a73e8;color:white;font-weight:bold;padding:9px;"
            "border-radius:5px;font-size:13px;"
        )
        self.btn_run.clicked.connect(self.run_curation)
        root.addWidget(self.btn_run)

        root.addWidget(QLabel("Execution Log:"))
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setStyleSheet("font-family:monospace;font-size:11px;background:#1e1e1e;color:#d4d4d4;")
        self.log_output.setMinimumHeight(110)
        root.addWidget(self.log_output)

    # ------------------------------------------------------------------
    def _build_data_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)

        grp = QGroupBox("File Selection")
        g = QVBoxLayout(grp)

        self.btn_source = QPushButton("📂  Browse Source CSV…")
        self.btn_source.clicked.connect(self.select_source)
        self.lbl_source = QLabel("No file selected.")
        self.lbl_source.setStyleSheet("color:gray;font-size:11px;")

        self.btn_output = QPushButton("💾  Browse Output CSV…")
        self.btn_output.clicked.connect(self.select_output)
        self.lbl_output = QLabel("No output file selected.")
        self.lbl_output.setStyleSheet("color:gray;font-size:11px;")

        g.addWidget(self.btn_source)
        g.addWidget(self.lbl_source)
        g.addWidget(self.btn_output)
        g.addWidget(self.lbl_output)
        lay.addWidget(grp)
        lay.addStretch()
        return w

    # ------------------------------------------------------------------
    def _build_columns_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)

        info = QLabel("Load a source CSV first, then assign columns.")
        info.setStyleSheet("color:#888;font-size:11px;")
        lay.addWidget(info)

        form = QFormLayout()
        form.setSpacing(8)

        self.combo_smiles = QComboBox()
        self.combo_smiles.setToolTip("Column containing SMILES strings")
        form.addRow("SMILES column:", self.combo_smiles)

        self.combo_activity = QComboBox()
        self.combo_activity.setToolTip("Column containing numeric activity values")
        form.addRow("Activity column:", self.combo_activity)

        lay.addLayout(form)

        grp2 = QGroupBox("Extra columns to keep in output")
        g2 = QVBoxLayout(grp2)
        hint = QLabel("Check any additional columns to retain (IDs, names, etc.)")
        hint.setStyleSheet("color:#888;font-size:10px;")
        g2.addWidget(hint)
        self.list_extra = QListWidget()
        self.list_extra.setMinimumHeight(120)
        g2.addWidget(self.list_extra)
        lay.addWidget(grp2)
        lay.addStretch()
        return w

    # ------------------------------------------------------------------
    def _build_filters_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)

        grp1 = QGroupBox("ChEMBL-Style Filters  (ignored when columns absent)")
        g1 = QFormLayout(grp1)
        g1.setSpacing(6)

        self.edit_type_filter = QLineEdit()
        self.edit_type_filter.setPlaceholderText("e.g. IC50,Ki  (comma-separated, blank = skip)")
        g1.addRow("Activity Type (Standard Type):", self.edit_type_filter)

        self.edit_unit_filter = QLineEdit()
        self.edit_unit_filter.setPlaceholderText("e.g. nM   (blank = skip)")
        g1.addRow("Unit Filter (Standard Units):", self.edit_unit_filter)

        self.chk_relation = QCheckBox("Remove ambiguous relations (keep '=' only)")
        g1.addRow("", self.chk_relation)
        lay.addWidget(grp1)

        grp2 = QGroupBox("Activity Value Range Filter")
        g2 = QFormLayout(grp2)
        g2.setSpacing(6)

        self.chk_act_range = QCheckBox("Enable range filter")
        g2.addRow("", self.chk_act_range)

        range_row = QHBoxLayout()
        self.spin_act_min = QDoubleSpinBox()
        self.spin_act_min.setRange(0, 1e12)
        self.spin_act_min.setValue(0.1)
        self.spin_act_min.setDecimals(4)
        self.spin_act_max = QDoubleSpinBox()
        self.spin_act_max.setRange(0, 1e12)
        self.spin_act_max.setValue(100000)
        self.spin_act_max.setDecimals(4)
        range_row.addWidget(QLabel("Min:"))
        range_row.addWidget(self.spin_act_min)
        range_row.addWidget(QLabel("Max:"))
        range_row.addWidget(self.spin_act_max)
        g2.addRow("Value range:", range_row)
        lay.addWidget(grp2)
        lay.addStretch()
        return w

    # ------------------------------------------------------------------
    def _build_options_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)

        grp1 = QGroupBox("Structural Curation Steps")
        g1 = QVBoxLayout(grp1)
        self.chk_org   = QCheckBox("Remove organometallics / inorganic atoms")
        self.chk_salt  = QCheckBox("Remove salts / retain largest fragment")
        self.chk_neut  = QCheckBox("Neutralize charges")
        self.chk_org.setChecked(True)
        self.chk_salt.setChecked(True)
        self.chk_neut.setChecked(True)
        for c in [self.chk_org, self.chk_salt, self.chk_neut]:
            g1.addWidget(c)
        lay.addWidget(grp1)

        grp2 = QGroupBox("Duplicate Handling")
        g2 = QFormLayout(grp2)
        self.combo_strategy = QComboBox()
        self.combo_strategy.addItems(['Mean (Average)', 'Median', 'Lowest (Most Potent)', 'Highest'])
        g2.addRow("Strategy:", self.combo_strategy)
        lay.addWidget(grp2)

        grp3 = QGroupBox("pIC50 / pKi Conversion  (−log₁₀)")
        g3 = QFormLayout(grp3)
        self.chk_pic50 = QCheckBox("Convert activity values to pActivity")
        self.combo_unit_pic50 = QComboBox()
        self.combo_unit_pic50.addItems(['nM', 'uM', 'mM', 'M'])
        g3.addRow("", self.chk_pic50)
        g3.addRow("Source unit:", self.combo_unit_pic50)
        lay.addWidget(grp3)

        lay.addStretch()
        return w

    # ------------------------------------------------------------------
    def _build_preview_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lbl = QLabel("Preview of curated output (first 50 rows):")
        lbl.setStyleSheet("font-size:11px;color:#555;")
        lay.addWidget(lbl)
        self.preview_table = QTableWidget()
        self.preview_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        lay.addWidget(self.preview_table)
        return w

    # ------------------------------------------------------------------
    def select_source(self):
        fn, _ = QFileDialog.getOpenFileName(self.widget, "Select Source CSV", "", "CSV Files (*.csv);;All Files (*)")
        if not fn:
            return
        self.source_path = fn
        short = fn.split('/')[-1].split('\\')[-1]
        self.lbl_source.setText(f"✔  {short}")
        self.lbl_source.setStyleSheet("color:#1a73e8;font-size:11px;")
        self._load_columns(fn)

    def select_output(self):
        fn, _ = QFileDialog.getSaveFileName(self.widget, "Select Output CSV", "curated_dataset.csv", "CSV Files (*.csv)")
        if fn:
            self.output_path = fn
            short = fn.split('/')[-1].split('\\')[-1]
            self.lbl_output.setText(f"✔  {short}")
            self.lbl_output.setStyleSheet("color:green;font-size:11px;")

    def _load_columns(self, path):
        """Populate column selectors after reading CSV header."""
        try:
            try:
                df_head = pd.read_csv(path, nrows=0)
            except Exception:
                df_head = pd.read_csv(path, nrows=0, encoding='latin1')
            cols = list(df_head.columns)
            self._df_columns = cols

            for combo in [self.combo_smiles, self.combo_activity]:
                combo.clear()
                combo.addItems(cols)

            # Smart defaults
            for i, c in enumerate(cols):
                cl = c.lower()
                if 'smiles' in cl or cl == 'smi':
                    self.combo_smiles.setCurrentIndex(i)
                    break
            for i, c in enumerate(cols):
                cl = c.lower()
                if any(k in cl for k in ['value', 'activity', 'ic50', 'ki', 'affinity', 'potency']):
                    self.combo_activity.setCurrentIndex(i)
                    break

            # Extra columns list
            self.list_extra.clear()
            for c in cols:
                item = QListWidgetItem(c)
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(Qt.CheckState.Unchecked)
                self.list_extra.addItem(item)

            self.log_message(f"📋  Detected {len(cols)} columns: {', '.join(cols[:8])}{'…' if len(cols) > 8 else ''}")
        except Exception as e:
            self.log_message(f"⚠️  Could not read columns: {e}")

    def log_message(self, msg):
        self.log_output.append(msg)

    def _get_extra_cols(self):
        extra = []
        for i in range(self.list_extra.count()):
            item = self.list_extra.item(i)
            if item.checkState() == Qt.CheckState.Checked:
                extra.append(item.text())
        return extra

    def run_curation(self):
        if not self.source_path or not self.output_path:
            QMessageBox.warning(self.widget, "Missing Files", "Please select both source and output CSV files.")
            return
        if not self.combo_smiles.currentText() or not self.combo_activity.currentText():
            QMessageBox.warning(self.widget, "Missing Columns", "Please assign SMILES and Activity columns.")
            return

        cfg = {
            'source_file'            : self.source_path,
            'output_file'            : self.output_path,
            'smiles_col'             : self.combo_smiles.currentText(),
            'activity_col'           : self.combo_activity.currentText(),
            'extra_cols'             : self._get_extra_cols(),
            'activity_type_filter'   : self.edit_type_filter.text(),
            'unit_filter'            : self.edit_unit_filter.text(),
            'filter_relations'       : self.chk_relation.isChecked(),
            'activity_min'           : self.spin_act_min.value() if self.chk_act_range.isChecked() else None,
            'activity_max'           : self.spin_act_max.value() if self.chk_act_range.isChecked() else None,
            'filter_organometallics' : self.chk_org.isChecked(),
            'remove_salts'           : self.chk_salt.isChecked(),
            'neutralize'             : self.chk_neut.isChecked(),
            'duplicate_strategy'     : self.combo_strategy.currentText(),
            'convert_pic50'          : self.chk_pic50.isChecked(),
            'act_unit'               : self.combo_unit_pic50.currentText(),
        }

        self.btn_run.setEnabled(False)
        self.log_output.clear()
        self.log_message("▶  Starting curation…")

        self.worker = CurationWorker(cfg)
        self.worker.log_signal.connect(self.log_message)
        self.worker.finished_signal.connect(self.curation_finished)
        self.worker.start()

    def curation_finished(self, success, message, result_df):
        self.btn_run.setEnabled(True)
        if success:
            self._populate_preview(result_df)
            QMessageBox.information(self.widget, "Success ✅", message)
        else:
            QMessageBox.critical(self.widget, "Error ❌", f"Curation failed:\n{message}")

    def _populate_preview(self, df):
        if df is None or df.empty:
            return
        preview = df.head(50)
        self.preview_table.clear()
        self.preview_table.setRowCount(len(preview))
        self.preview_table.setColumnCount(len(preview.columns))
        self.preview_table.setHorizontalHeaderLabels(list(preview.columns))
        for r, (_, row) in enumerate(preview.iterrows()):
            for c, val in enumerate(row):
                cell = QTableWidgetItem(str(val))
                cell.setFlags(cell.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.preview_table.setItem(r, c, cell)


# ---------------------------------------------------------------------------
# Plugin class
# ---------------------------------------------------------------------------
class QsarCuratorPlugin(BasePlugin):
    def __init__(self):
        super().__init__(PluginInfo(
            name="QSAR Dataset Curator",
            version="2.0.0",
            description="Universal QSAR curator — works with ChEMBL, BindingDB, PubChem, or any CSV",
            author="SMILES Team",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        ))
        self.widget = None

    def get_info(self) -> PluginInfo:
        return self.info

    def create_widget(self) -> 'QsarCuratorWidget':
        if self.widget is None:
            self.widget = QsarCuratorWidget(self)
        self.logger.info("QSAR Curator widget created")
        return self.widget

    def initialize(self):
        self.logger.info("QSAR Curator plugin initialized")
        return True

    def cleanup(self):
        if self.widget:
            self.widget.widget.deleteLater()
            self.widget = None
        self.logger.info("QSAR Curator plugin cleaned up")