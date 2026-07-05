import os
import csv
from pathlib import Path
import traceback

from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
    QLineEdit, QFileDialog, QComboBox, QProgressBar, QTextEdit, 
    QThread, Signal, Qt, QSpinBox
)

from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType
from src.features.smiles_parser.services.parser import parse_smiles
from src.services.forcefield.mmff94_service import MMFF94Service
from src.features.io.exporters.sdf_writer import write_sdf
from src.features.io.exporters.mol2_writer import write_mol2
import concurrent.futures
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing

def process_molecule(mol_id, smiles, output_dir, output_format):
    # Imports inside function for multiprocessing compatibility
    from src.features.smiles_parser.services.parser import parse_smiles
    from src.services.forcefield.mmff94_service import MMFF94Service
    from src.features.io.exporters.sdf_writer import write_sdf
    from src.features.io.exporters.mol2_writer import write_mol2
    from pathlib import Path
    
    try:
        mol = parse_smiles(smiles, use_bkchem_tokenizer=True)
        if not mol or len(mol.atoms) == 0:
            raise ValueError(f"Failed to parse SMILES: {smiles}")
            
        mol.name = mol_id
        
        # Pre-seed 3D coordinates to 0.0 to prevent NoneType errors during H addition
        for a in mol.atoms:
            if getattr(a, 'x', None) is None: a.x = 0.0
            if getattr(a, 'y', None) is None: a.y = 0.0
            if getattr(a, 'z', None) is None: a.z = 0.0
        
        mmff94 = MMFF94Service()
        mmff94.add_hydrogens(mol)
        mmff94.optimize_geometry(mol, max_iters=500, convergence=1e-4)
        
        out_path = Path(output_dir)
        if output_format == "SDF":
            filepath = str(out_path / f"{mol_id}.sdf")
            write_sdf(mol, filepath)
        elif output_format == "MOL2":
            filepath = str(out_path / f"{mol_id}.mol2")
            write_mol2(mol, filepath)
            
        return mol_id, True, filepath, None
    except Exception as e:
        return mol_id, False, None, str(e)

class BatchConversionWorker(QThread):
    progress_updated = Signal(int, int) # current, total
    log_message = Signal(str)
    finished_conversion = Signal(int, int) # success, failed
    
    def __init__(self, csv_path, output_dir, output_format, num_threads):
        super().__init__()
        self.csv_path = csv_path
        self.output_dir = output_dir
        self.output_format = output_format
        self.num_threads = num_threads
        self.is_running = True
        
    def run(self):
        try:
            self.log_message.emit("Starting batch conversion...")
            rows = []
            with open(self.csv_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                is_first = True
                for r in reader:
                    if len(r) >= 2:
                        mol_id_val = r[0].strip()
                        smiles_val = r[1].strip()
                        
                        # Skip header row if present
                        if is_first:
                            is_first = False
                            if "smiles" in smiles_val.lower() or "chembl" in mol_id_val.lower():
                                continue
                                
                        rows.append((mol_id_val, smiles_val))
                    is_first = False
            
            total = len(rows)
            if total == 0:
                self.log_message.emit("CSV file is empty or missing columns.")
                self.finished_conversion.emit(0, 0)
                return
                
            self.log_message.emit(f"Found {total} rows to process.")
            
            success_count = 0
            failed_count = 0
            
            num_workers = self.num_threads
            self.log_message.emit(f"Using {num_workers} parallel workers (Thread Pool).")
            
            with ThreadPool(processes=num_workers) as pool:
                # Submit all tasks
                results = []
                for mol_id, smiles in rows:
                    res = pool.apply_async(process_molecule, (mol_id, smiles, self.output_dir, self.output_format))
                    results.append(res)
                
                completed = 0
                for res in results:
                    if not self.is_running:
                        self.log_message.emit("Conversion stopped by user. Cancelling remaining tasks...")
                        pool.terminate()
                        break
                        
                    mol_id, success, filepath, error_msg = res.get()
                    completed += 1
                    
                    self.log_message.emit(f"Processing ({completed}/{total}): {mol_id}")
                    if success:
                        success_count += 1
                        self.log_message.emit(f"  Success -> {filepath}")
                    else:
                        failed_count += 1
                        self.log_message.emit(f"  Error for {mol_id}: {error_msg}")
                        
                    self.progress_updated.emit(completed, total)
                
            self.log_message.emit("Batch conversion complete.")
            self.finished_conversion.emit(success_count, failed_count)
            
        except Exception as e:
            self.log_message.emit(f"Fatal error during batch conversion: {str(e)}\n{traceback.format_exc()}")
            self.finished_conversion.emit(0, 0)
            
    def stop(self):
        self.is_running = False

class Batch3DConverterWidget(QWidget):
    def __init__(self, plugin):
        super().__init__()
        self.plugin = plugin
        self.worker = None
        self._setup_ui()
        
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Title
        title = QLabel("Batch SMILES to 3D Converter")
        title.setStyleSheet("font-size: 16px; font-weight: bold; margin-bottom: 10px;")
        layout.addWidget(title)
        
        # Instructions
        layout.addWidget(QLabel("Select a CSV file containing Mol ID in the first column and SMILES in the second."))
        
        # Source File Selection
        h1 = QHBoxLayout()
        h1.addWidget(QLabel("CSV File:"))
        self.csv_field = QLineEdit()
        self.csv_field.setReadOnly(True)
        h1.addWidget(self.csv_field)
        btn_browse_csv = QPushButton("Browse...")
        btn_browse_csv.clicked.connect(self._browse_csv)
        h1.addWidget(btn_browse_csv)
        layout.addLayout(h1)
        
        # Output Directory Selection
        h2 = QHBoxLayout()
        h2.addWidget(QLabel("Output Dir:"))
        self.out_field = QLineEdit()
        self.out_field.setReadOnly(True)
        h2.addWidget(self.out_field)
        btn_browse_out = QPushButton("Browse...")
        btn_browse_out.clicked.connect(self._browse_out)
        h2.addWidget(btn_browse_out)
        layout.addLayout(h2)
        
        # Format Selection and Threads
        h3 = QHBoxLayout()
        h3.addWidget(QLabel("Output Format:"))
        self.format_combo = QComboBox()
        self.format_combo.addItems(["SDF", "MOL2"])
        h3.addWidget(self.format_combo)
        
        h3.addSpacing(20)
        h3.addWidget(QLabel("Threads:"))
        self.thread_spin = QSpinBox()
        total_cores = multiprocessing.cpu_count()
        self.thread_spin.setRange(1, total_cores * 2)
        self.thread_spin.setValue(max(1, total_cores // 2))
        h3.addWidget(self.thread_spin)
        
        h3.addStretch()
        layout.addLayout(h3)
        
        # Action Buttons
        h_btn = QHBoxLayout()
        self.btn_start = QPushButton("Start Conversion")
        self.btn_start.clicked.connect(self._start_conversion)
        self.btn_start.setStyleSheet("background-color: #2e7d32; color: white;")
        
        self.btn_stop = QPushButton("Stop")
        self.btn_stop.clicked.connect(self._stop_conversion)
        self.btn_stop.setEnabled(False)
        self.btn_stop.setStyleSheet("background-color: #c62828; color: white;")
        
        h_btn.addWidget(self.btn_start)
        h_btn.addWidget(self.btn_stop)
        layout.addLayout(h_btn)
        
        # Progress Bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        layout.addWidget(self.progress_bar)
        
        # Log Output
        self.log_area = QTextEdit()
        self.log_area.setReadOnly(True)
        layout.addWidget(self.log_area)
        
    def _browse_csv(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select CSV", "", "CSV Files (*.csv)")
        if file:
            self.csv_field.setText(file)
            
    def _browse_out(self):
        dir_ = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_:
            self.out_field.setText(dir_)
            
    def _start_conversion(self):
        csv_path = self.csv_field.text()
        out_dir = self.out_field.text()
        
        if not csv_path or not out_dir:
            self.log_area.append("Please select both a CSV file and an output directory.")
            return
            
        fmt = self.format_combo.currentText()
        threads = self.thread_spin.value()
        
        self.btn_start.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.progress_bar.setValue(0)
        self.log_area.clear()
        
        self.worker = BatchConversionWorker(csv_path, out_dir, fmt, threads)
        self.worker.progress_updated.connect(self._on_progress)
        self.worker.log_message.connect(self._on_log)
        self.worker.finished_conversion.connect(self._on_finished)
        self.worker.start()
        
    def _stop_conversion(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.btn_stop.setEnabled(False)
            
    def _on_progress(self, current, total):
        if total > 0:
            val = int((current / total) * 100)
            self.progress_bar.setValue(val)
            
    def _on_log(self, msg):
        self.log_area.append(msg)
        # scroll to bottom
        scrollbar = self.log_area.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())
        
    def _on_finished(self, success, failed):
        self.btn_start.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self.log_area.append(f"\n--- Conversion Summary ---")
        self.log_area.append(f"Successfully converted: {success}")
        self.log_area.append(f"Failed: {failed}")
        
    def cleanup(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.worker.wait(2000)

class Batch3DConverterPlugin(BasePlugin):
    def __init__(self):
        info = PluginInfo(
            name="Batch 3D Converter",
            version="1.0.0",
            description="Batch converts SMILES from a CSV file into 3D SDF/Mol2 files using MMFF94.",
            author="DeepMind Agent",
            plugin_type=PluginType.UTILITY
        )
        super().__init__(info)
        self._widget = None
        
    def create_widget(self):
        if self._widget is None:
            self._widget = Batch3DConverterWidget(self)
        return self._widget
        
    def cleanup(self):
        if self._widget:
            self._widget.cleanup()
        super().cleanup()
