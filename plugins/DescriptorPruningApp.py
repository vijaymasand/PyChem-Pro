"""
Descriptor Pruning Plugin

Auto-discovers the Dependent Variable and prunes highly constant
and correlated descriptors across multiple CSV files.
"""

import os
import pandas as pd
import numpy as np
from sklearn.feature_selection import VarianceThreshold

# Strictly using ONLY the imports allowed by your qt_compat file
from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QTextEdit, QProgressBar,
    Qt, Signal, QFileDialog, QMessageBox
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType

class DescriptorPruningWidget(PluginWidget):
    """Qt Widget for the Descriptor Pruning Plugin."""

    def __init__(self, plugin: 'DescriptorPruningPlugin'):
        super().__init__(plugin)
        self.desc_files = []
        self.setup_ui()

    def setup_ui(self):
        self.widget = QWidget()
        layout = QVBoxLayout(self.widget)

        # 1. File Selection
        layout.addWidget(QLabel("<b>1. Select CSV Files</b>"))

        file_layout = QHBoxLayout()
        self.file_table = QTableWidget()
        self.file_table.setColumnCount(1)
        self.file_table.setHorizontalHeaderLabels(["Selected Files"])
        self.file_table.horizontalHeader().setStretchLastSection(True)
        self.file_table.setMaximumHeight(100)
        file_layout.addWidget(self.file_table)

        self.browse_btn = QPushButton("Browse Files")
        self.browse_btn.clicked.connect(self.select_files)
        file_layout.addWidget(self.browse_btn)
        layout.addLayout(file_layout)

        # 2. Settings & Thresholds (Built using an editable QTableWidget)
        layout.addWidget(QLabel("<b>2. Settings & Thresholds</b>"))

        self.settings_table = QTableWidget()
        self.settings_table.setColumnCount(2)
        self.settings_table.setHorizontalHeaderLabels(["Setting", "Value"])
        self.settings_table.horizontalHeader().setStretchLastSection(True)

        # Define our default settings
        settings_data = [
            ("Dependent Variable (DV) Name", "pIC50"),
            ("Threshold limit for constant", "0.90"),
            ("Minimum Variance", "0.001"),
            ("Min Correlation with DV", "0.05"),
            ("Max Inter-Descriptor Corr", "0.90"),
            ("Log Top 'N' Influential", "10")
        ]

        self.settings_table.setRowCount(len(settings_data))
        for i, (name, val) in enumerate(settings_data):
            name_item = QTableWidgetItem(name)
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable) # Make name read-only
            val_item = QTableWidgetItem(val)

            self.settings_table.setItem(i, 0, name_item)
            self.settings_table.setItem(i, 1, val_item)

        layout.addWidget(self.settings_table)

        # 3. Execution UI
        layout.addWidget(QLabel("<b>3. Execute</b>"))

        self.progress = QProgressBar()
        self.progress.setValue(0)
        layout.addWidget(self.progress)

        self.run_btn = QPushButton("Run Auto-Discovery & Pruning")
        self.run_btn.setStyleSheet("background-color: #2980b9; color: white; font-weight: bold; padding: 8px;")
        self.run_btn.clicked.connect(self.run_pipeline)
        layout.addWidget(self.run_btn)

        self.status_text = QTextEdit()
        self.status_text.setReadOnly(True)
        layout.addWidget(self.status_text)

    def select_files(self):
        files, _ = QFileDialog.getOpenFileNames(self.widget, "Select CSV files", "", "CSV Files (*.csv)")
        if files:
            self.desc_files = files
            self.file_table.setRowCount(len(files))
            for i, f in enumerate(files):
                self.file_table.setItem(i, 0, QTableWidgetItem(os.path.basename(f)))
            self.status_text.append(f"Loaded {len(files)} files.")

    def run_pipeline(self):
        if not self.desc_files:
            QMessageBox.warning(self.widget, "Error", "No files selected.")
            return

        # Extract values from the settings table safely
        try:
            target_dv = self.settings_table.item(0, 1).text().strip()
            p_const = float(self.settings_table.item(1, 1).text())
            var_thresh = float(self.settings_table.item(2, 1).text())
            corr_dv_min = float(self.settings_table.item(3, 1).text())
            corr_inter = float(self.settings_table.item(4, 1).text())
            top_n = int(self.settings_table.item(5, 1).text())
        except ValueError:
            QMessageBox.warning(self.widget, "Input Error", "Please ensure all numeric settings contain valid numbers.")
            return

        save_path, _ = QFileDialog.getSaveFileName(self.widget, "Save Master CSV", "Final_Master_Descriptors.csv", "CSV Files (*.csv)")
        if not save_path:
            return

        try:
            self.run_btn.setEnabled(False)
            self.status_text.append("Starting pipeline...")
            y_data = None
            ids = None

            def read_csv_safe(path, **kwargs):
                try:
                    return pd.read_csv(path, **kwargs)
                except Exception as e:
                    err_name = type(e).__name__
                    if err_name.startswith('Uni') and err_name.endswith('Error'):
                        return pd.read_csv(path, encoding='latin1', **kwargs)
                    raise

            # --- PHASE 1: AUTO-DISCOVERY ---
            for f in self.desc_files:
                df_head = read_csv_safe(f, nrows=2)
                if target_dv in df_head.columns:
                    full_df = read_csv_safe(f)
                    y_data = full_df[target_dv]
                    ids = full_df.iloc[:, 0]
                    self.status_text.append(f"Found DV '{target_dv}' in {os.path.basename(f)}")
                    break

            if y_data is None:
                QMessageBox.critical(self.widget, "Error", f"Could not find '{target_dv}' in any file.")
                self.run_btn.setEnabled(True)
                return

            # --- PHASE 2: PROCESSING ---
            all_processed_dfs = []
            for i, f_path in enumerate(self.desc_files):
                self.status_text.append(f"Processing {os.path.basename(f_path)}...")
                df = read_csv_safe(f_path)
                X = df.select_dtypes(include=[np.number])

                if target_dv in X.columns:
                    X = X.drop(columns=[target_dv])

                # 1. Constant/Low Variance Filter
                vt = VarianceThreshold(threshold=(p_const * (1 - p_const)))
                X = X.loc[:, vt.fit(X).get_support()]
                X = X.loc[:, X.var() > var_thresh]

                # 2. Correlation with DV Filter
                corrs_with_y = X.apply(lambda col: col.corr(y_data)).abs()
                X = X.loc[:, corrs_with_y >= corr_dv_min]

                # 3. Inter-correlation Tournament
                corr_matrix = X.corr().abs()
                upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))

                to_drop = []
                for col in upper.columns:
                    high_corr = upper.index[upper[col] > corr_inter]
                    for feature in high_corr:
                        if corrs_with_y[col] >= corrs_with_y[feature]:
                            to_drop.append(feature)
                        else:
                            to_drop.append(col)

                all_processed_dfs.append(X.drop(columns=list(set(to_drop))))
                self.progress.setValue(int(((i + 1) / len(self.desc_files)) * 100))

            # --- PHASE 3: FINAL MERGE & OUTPUT ---
            master = pd.concat(all_processed_dfs, axis=1)
            final_corrs = master.apply(lambda col: col.corr(y_data))

            # Log generation
            log_filename = save_path.replace('.csv', '_Verification_Log.txt')
            top_features = final_corrs.abs().sort_values(ascending=False).head(top_n)

            with open(log_filename, 'w') as f_log:
                f_log.write("CSV Pruning Log\n")
                f_log.write(f"Retained Descriptors: {master.shape[1]}\n")
                f_log.write("-" * 40 + "\n")
                for feat, val in top_features.items():
                    f_log.write(f"{feat}: r = {final_corrs[feat]:.4f}\n")

            master.insert(0, 'Molecule_ID', ids)
            master[target_dv] = y_data.values

            master.to_csv(save_path, index=False)

            self.status_text.append(f"<b>Success!</b> Retained {master.shape[1]-2} descriptors.")
            QMessageBox.information(self.widget, "Success", f"Master File and Log Created successfully at:\n{save_path}")

        except Exception as e:
            self.status_text.append(f"<span style='color:red'>Error: {str(e)}</span>")
            QMessageBox.critical(self.widget, "Processing Error", str(e))
        finally:
            self.run_btn.setEnabled(True)


class DescriptorPruningPlugin(BasePlugin):
    """
    Plugin wrapper for the Descriptor Pruning Tool.
    """
    def __init__(self):
        super().__init__(PluginInfo(
            name="Descriptor Pruning (Auto-Discovery)",
            version="1.0.0",
            description="Removes highly constant and correlated descriptors across CSV datasets.",
            author="SMILES Team",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        ))
        self.widget = None

    def get_info(self) -> PluginInfo:
        return self.info

    def create_widget(self) -> 'DescriptorPruningWidget':
        if self.widget is None:
            self.widget = DescriptorPruningWidget(self)
        return self.widget

    # FIX: Added *args, **kwargs or api=None so the software can pass the api object without crashing
    def initialize(self, api=None):
        self.logger.info("Descriptor Pruning plugin initialized.")
        return True

    def cleanup(self):
        if self.widget:
            self.widget.widget.deleteLater()
            self.widget = None
        self.logger.info("Descriptor Pruning plugin cleaned up.")