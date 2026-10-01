from src.shared.qt_compat import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                                 QDoubleSpinBox, QCheckBox, QPushButton, QComboBox, 
                                 Qt, QFileDialog, Signal)
import os

class DockingPoseDialog(QDialog):
    pick_requested = Signal() # Request picking from viewer

    def __init__(self, ligands, parent=None):
        super().__init__(parent)
        self.setWindowTitle("3D Molecular Docking Pose")
        self.setMinimumWidth(350)
        
        self.ligands = ligands # List of Set[int]
        self.selected_ligand_idx = 0
        
        self._init_ui()

    def _init_ui(self):
        layout = QVBoxLayout(self)
        
        # Ligand Selection
        if len(self.ligands) > 1:
            layout.addWidget(QLabel("Select target ligand fragment:"))
            self.ligand_combo = QComboBox()
            for i, frag in enumerate(self.ligands):
                self.ligand_combo.addItem(f"Fragment {i+1} ({len(frag)} atoms)")
            layout.addWidget(self.ligand_combo)
        else:
            layout.addWidget(QLabel("Ligand detected automatically."))
            self.ligand_combo = QComboBox()
            if self.ligands:
                self.ligand_combo.addItem(f"Detected Fragment ({len(self.ligands[0])} atoms)")
            self.ligand_combo.setEnabled(False)
            layout.addWidget(self.ligand_combo)

        # Pick Button
        self.pick_btn = QPushButton("Pick Ligand from 3D View")
        self.pick_btn.clicked.connect(self._on_pick_clicked)
        layout.addWidget(self.pick_btn)

        # Distance Limit
        dist_layout = QHBoxLayout()
        dist_layout.addWidget(QLabel("Nearby residue distance threshold (Å):"))
        self.dist_spin = QDoubleSpinBox()
        self.dist_spin.setRange(2.0, 15.0)
        self.dist_spin.setValue(5.0)
        dist_layout.addWidget(self.dist_spin)
        layout.addLayout(dist_layout)
        
        # Interaction Toggles
        layout.addWidget(QLabel("Visual Interactions:"))
        self.hb_check = QCheckBox("Hydrogen Bonds")
        self.hb_check.setChecked(True)
        layout.addWidget(self.hb_check)
        
        self.salt_check = QCheckBox("Salt Bridges / Ionic")
        self.salt_check.setChecked(True)
        layout.addWidget(self.salt_check)
        
        self.hydro_check = QCheckBox("Hydrophobic Contacts")
        self.hydro_check.setChecked(False)
        layout.addWidget(self.hydro_check)
        
        layout.addSpacing(10)
        
        # Report Button
        self.report_btn = QPushButton("Save Interaction Report (CSV)...")
        self.report_btn.clicked.connect(self._on_save_report)
        layout.addWidget(self.report_btn)
        
        layout.addSpacing(10)
        
        # Quick Action Buttons
        layout.addWidget(QLabel("Quick Actions:"))
        action_layout = QHBoxLayout()
        
        self.label_btn = QPushButton("Label Nearby Residues")
        self.label_btn.clicked.connect(self._on_label_nearby)
        self.label_btn.setToolTip("Label residues within 5.0 Å of ligand")
        action_layout.addWidget(self.label_btn)
        
        self.clear_labels_btn = QPushButton("Clear Labels")
        self.clear_labels_btn.clicked.connect(self._on_clear_labels)
        self.clear_labels_btn.setToolTip("Clear all residue labels")
        action_layout.addWidget(self.clear_labels_btn)
        
        self.zoom_btn = QPushButton("Zoom to Ligand")
        self.zoom_btn.clicked.connect(self._on_zoom_to_ligand)
        self.zoom_btn.setToolTip("Zoom to ligand with 7 Å surrounding area")
        action_layout.addWidget(self.zoom_btn)
        
        layout.addLayout(action_layout)
        layout.addSpacing(10)
        
        # Buttons
        btn_layout = QHBoxLayout()
        self.apply_btn = QPushButton("Apply View")
        self.apply_btn.setDefault(True)
        self.apply_btn.clicked.connect(self.accept)
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.reject)
        
        btn_layout.addWidget(self.cancel_btn)
        btn_layout.addWidget(self.apply_btn)
        layout.addLayout(btn_layout)

    def _on_pick_clicked(self):
        """Request the main window to pick an atom."""
        self.hide()
        # Note: Signal emission is handled by the caller who created the dialog
        # but here we just use a callback or a flag if signal is tricky in this context
        self.pick_requested_flag = True
        self.done(100) # Custom return code for picking

    def select_fragment_containing(self, atom_idx):
        """Select the fragment that contains the given atom index."""
        for i, frag in enumerate(self.ligands):
            if atom_idx in frag:
                if self.ligand_combo:
                    self.ligand_combo.setCurrentIndex(i)
                return True
        return False

    def _on_save_report(self):
        self.save_report_requested = True
        self.accept()

    def _on_label_nearby(self):
        """Label residues within 5.0 Å of ligand."""
        self.label_nearby_requested = True
        self.accept()

    def _on_clear_labels(self):
        """Clear all residue labels."""
        self.clear_labels_requested = True
        self.accept()

    def _on_zoom_to_ligand(self):
        """Zoom to ligand with 7 Å surrounding area."""
        self.zoom_to_ligand_requested = True
        self.accept()

    def get_config(self):
        l_idx = self.ligand_combo.currentIndex() if self.ligand_combo else 0
        return {
            'ligand_indices': self.ligands[l_idx],
            'distance': self.dist_spin.value(),
            'show_hbonds': self.hb_check.isChecked(),
            'show_salt': self.salt_check.isChecked(),
            'show_hydro': self.hydro_check.isChecked(),
            'save_report': getattr(self, 'save_report_requested', False),
            'label_nearby': getattr(self, 'label_nearby_requested', False),
            'clear_labels': getattr(self, 'clear_labels_requested', False),
            'zoom_to_ligand': getattr(self, 'zoom_to_ligand_requested', False)
        }
