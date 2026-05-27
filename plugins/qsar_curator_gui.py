"""
QSAR Dataset Curator Plugin

A data curation tool that standardizes SMILES, neutralizes charges,
removes salts/organometallics, and handles duplicates from ChEMBL CSVs.
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger

# Disable RDKit warnings to keep the console/logs clean
RDLogger.DisableLog('rdApp.*')

# Strictly using ONLY the imports allowed by the host application's qt_compat
from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QFileDialog, QComboBox, QTextEdit, QMessageBox,
    QThread, Signal, Qt
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType


class CurationWorker(QThread):
    """Background thread to process the dataset without freezing the GUI."""
    log_signal = Signal(str)
    finished_signal = Signal(bool, str)

    def __init__(self, source_file, output_file, duplicate_strategy):
        super().__init__()
        self.source_file = source_file
        self.output_file = output_file
        self.duplicate_strategy = duplicate_strategy

        # Initialize RDKit standardizers
        self.lfc = rdMolStandardize.LargestFragmentChooser()
        self.uc = rdMolStandardize.Uncharger()

        # Allowed atomic numbers (Organic subset: H, B, C, N, O, F, Si, P, S, Cl, Br, I)
        self.allowed_atoms = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}

    def curate_molecule(self, smiles):
        """Applies QSAR standardization rules to a single SMILES string."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            # 1. Filter out organometallics / inorganics
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() not in self.allowed_atoms:
                    return None

            # 2. Retain largest fragment (removes salts/solvents)
            mol = self.lfc.choose(mol)

            # 3. Neutralize charges
            mol = self.uc.uncharge(mol)

            # 4. Generate canonical SMILES
            return Chem.MolToSmiles(mol)
        except Exception:
            return None

    def run(self):
        try:
            self.log_signal.emit("Loading dataset...")
            def read_csv_safe(path, **kwargs):
                try:
                    return pd.read_csv(path, **kwargs)
                except Exception as e:
                    err_name = type(e).__name__
                    if err_name.startswith('Uni') and err_name.endswith('Error'):
                        return pd.read_csv(path, encoding='latin1', **kwargs)
                    raise

            # Load CSV. Pandas handles standard comma separation natively.
            df = read_csv_safe(self.source_file)
            initial_count = len(df)
            self.log_signal.emit(f"Loaded {initial_count} rows.")

            # Filter 1: Keep only IC50 or Ki
            if 'Standard Type' in df.columns:
                df = df[df['Standard Type'].isin(['IC50', 'Ki'])]
                self.log_signal.emit(f"Rows after filtering Standard Type (IC50/Ki): {len(df)}")

            # Filter 2: Keep only nanomolar (nM) units
            if 'Standard Units' in df.columns:
                df = df[df['Standard Units'] == 'nM']
                self.log_signal.emit(f"Rows after filtering Standard Units (nM): {len(df)}")

            # Filter 3: Remove ambiguous relations (Keep only '=')
            if 'Standard Relation' in df.columns:
                df['Standard Relation'] = df['Standard Relation'].astype(str).str.replace("'", "").str.strip()
                df = df[df['Standard Relation'] == '=']
                self.log_signal.emit(f"Rows after removing ambiguous relations: {len(df)}")

            # Drop missing values in critical columns and ensure numerics
            df = df.dropna(subset=['Smiles', 'Standard Value'])
            df['Standard Value'] = pd.to_numeric(df['Standard Value'], errors='coerce')
            df = df.dropna(subset=['Standard Value'])

            # Apply RDKit Curation
            self.log_signal.emit("Applying RDKit structural curation (salts, neutralization, organometallics)...")
            df['Canonical_Smiles'] = df['Smiles'].apply(self.curate_molecule)

            # Drop structures that failed curation (e.g., organometallics)
            df = df.dropna(subset=['Canonical_Smiles'])
            self.log_signal.emit(f"Rows after structural curation: {len(df)}")

            # Handle Duplicates
            self.log_signal.emit(f"Handling duplicates using strategy: {self.duplicate_strategy}...")
            if self.duplicate_strategy == 'Mean (Average)':
                agg_func = 'mean'
            elif self.duplicate_strategy == 'Median':
                agg_func = 'median'
            else:  # Lowest (Most potent)
                agg_func = 'min'

            # Group by the canonical SMILES and aggregate the activity value
            agg_dict = {'Standard Value': agg_func}
            if 'Molecule ChEMBL ID' in df.columns: agg_dict['Molecule ChEMBL ID'] = 'first'
            if 'Standard Type' in df.columns: agg_dict['Standard Type'] = 'first'
            if 'Standard Units' in df.columns: agg_dict['Standard Units'] = 'first'

            curated_df = df.groupby('Canonical_Smiles', as_index=False).agg(agg_dict)

            final_count = len(curated_df)
            self.log_signal.emit(f"Final curated dataset contains {final_count} unique molecules.")

            # Save to CSV
            curated_df.to_csv(self.output_file, index=False)
            self.log_signal.emit(f"Successfully saved to {self.output_file}")

            self.finished_signal.emit(True, "Curation completed successfully!")

        except Exception as e:
            self.finished_signal.emit(False, str(e))


class QsarCuratorWidget(PluginWidget):
    """
    Widget for the QSAR Curator plugin.
    """
    def __init__(self, plugin: 'QsarCuratorPlugin'):
        super().__init__(plugin)
        self.source_path = ""
        self.output_path = ""
        self.worker = None
        self.setup_ui()

    def setup_ui(self):
        self.widget = QWidget()
        layout = QVBoxLayout(self.widget)

        # Title & Info
        title = QLabel("QSAR Dataset Curator")
        title.setStyleSheet("font-size: 16px; font-weight: bold; margin-bottom: 5px;")
        layout.addWidget(title)

        info = QLabel("Standardizes SMILES, neutralizes charges, removes organometallics, and filters ChEMBL data.")
        info.setStyleSheet("color: #555; margin-bottom: 10px;")
        layout.addWidget(info)

        # File Selection Area
        file_layout = QVBoxLayout()

        self.btn_source = QPushButton("Select Source ChEMBL CSV")
        self.btn_source.clicked.connect(self.select_source)
        self.lbl_source = QLabel("No source file selected.")
        self.lbl_source.setStyleSheet("color: gray;")

        self.btn_output = QPushButton("Select Output Destination")
        self.btn_output.clicked.connect(self.select_output)
        self.lbl_output = QLabel("No output file selected.")
        self.lbl_output.setStyleSheet("color: gray;")

        file_layout.addWidget(self.btn_source)
        file_layout.addWidget(self.lbl_source)
        file_layout.addWidget(self.btn_output)
        file_layout.addWidget(self.lbl_output)
        layout.addLayout(file_layout)

        # Duplicate Handling Options
        strat_layout = QHBoxLayout()
        strat_layout.addWidget(QLabel("Duplicate Handling Strategy:"))
        self.combo_strategy = QComboBox()
        self.combo_strategy.addItems(['Mean (Average)', 'Median', 'Lowest (Most Potent)'])
        strat_layout.addWidget(self.combo_strategy)
        layout.addLayout(strat_layout)

        # Run Button
        self.btn_run = QPushButton("🚀 Run Curation")
        self.btn_run.setStyleSheet("background-color: #2b5c8f; color: white; font-weight: bold; padding: 8px;")
        self.btn_run.clicked.connect(self.run_curation)
        layout.addWidget(self.btn_run)

        # Log Output
        layout.addWidget(QLabel("Execution Log:"))
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        layout.addWidget(self.log_output)

    def select_source(self):
        file_name, _ = QFileDialog.getOpenFileName(self.widget, "Select Source CSV", "", "CSV Files (*.csv)")
        if file_name:
            self.source_path = file_name
            self.lbl_source.setText(self.source_path)
            self.lbl_source.setStyleSheet("color: black;")

    def select_output(self):
        file_name, _ = QFileDialog.getSaveFileName(self.widget, "Select Output CSV", "curated_dataset.csv", "CSV Files (*.csv)")
        if file_name:
            self.output_path = file_name
            self.lbl_output.setText(self.output_path)
            self.lbl_output.setStyleSheet("color: black;")

    def log_message(self, message):
        self.log_output.append(message)

    def run_curation(self):
        if not self.source_path or not self.output_path:
            QMessageBox.warning(self.widget, "Missing Information", "Please select both source and output files.")
            return

        self.btn_run.setEnabled(False)
        self.log_output.clear()
        self.log_message("Starting curation process...")

        strategy = self.combo_strategy.currentText()

        # Initialize and start the background thread
        self.worker = CurationWorker(self.source_path, self.output_path, strategy)
        self.worker.log_signal.connect(self.log_message)
        self.worker.finished_signal.connect(self.curation_finished)
        self.worker.start()

    def curation_finished(self, success, message):
        self.btn_run.setEnabled(True)
        if success:
            QMessageBox.information(self.widget, "Success", message)
        else:
            QMessageBox.critical(self.widget, "Error", f"An error occurred during curation:\n{message}")


class QsarCuratorPlugin(BasePlugin):
    """
    QSAR Curator Plugin.
    """
    def __init__(self):
        super().__init__(PluginInfo(
            name="QSAR Dataset Curator",
            version="1.0.0",
            description="Standardizes SMILES and cleans ChEMBL datasets",
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
        """
        Bypass super().initialize() to avoid application argument crash.
        """
        self.logger.info("QSAR Curator plugin initialized successfully")
        return True

    def cleanup(self):
        if self.widget:
            self.widget.widget.deleteLater()
            self.widget = None
        self.logger.info("QSAR Curator plugin cleaned up")