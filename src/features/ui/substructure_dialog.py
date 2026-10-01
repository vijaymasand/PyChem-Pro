"""
Substructure Matcher GUI Dialog.
"""
from src.shared.qt_compat import *
from src.shared.ui.theme import COLORS

class SubstructureDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("SMILES/SMARTS Substructure Matcher")
        self.setMinimumWidth(450)
        
        # Apply theme
        self.setStyleSheet(f"""
            QDialog {{ background-color: {COLORS['bg_primary']}; color: {COLORS['text_primary']}; }}
            QLabel {{ color: {COLORS['text_primary']}; }}
            QLineEdit {{
                background-color: {COLORS['bg_widget']};
                color: {COLORS['text_primary']};
                border: 1px solid {COLORS['border']};
                border-radius: 4px;
                padding: 6px;
            }}
            QPushButton {{
                background-color: {COLORS['accent']};
                color: {COLORS['bg_primary']};
                border: none;
                border-radius: 4px;
                padding: 6px 12px;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: #5555ff;
            }}
        """)
        
        layout = QVBoxLayout(self)
        
        info_label = QLabel("Enter a SMILES/SMARTS string to find inside the loaded molecule:")
        layout.addWidget(info_label)
        
        search_layout = QHBoxLayout()
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("e.g. c1ccccc1")
        search_layout.addWidget(self.search_input)
        
        self.search_btn = QPushButton("Find && Highlight")
        search_layout.addWidget(self.search_btn)
        
        layout.addLayout(search_layout)
        
        self.status_label = QLabel("")
        self.status_label.setStyleSheet(f"color: {COLORS['text_secondary']}; font-style: italic;")
        layout.addWidget(self.status_label)
        
        self.search_btn.clicked.connect(self._on_search)
        self.search_input.returnPressed.connect(self._on_search)

    def _on_search(self):
        query = self.search_input.text().strip()
        if not query:
            self.status_label.setText("Please enter a valid query.")
            return
            
        main_window = self.parent()
        if not main_window or not hasattr(main_window, 'molecule') or not main_window.molecule:
            self.status_label.setText("No molecule loaded in the Main Window.")
            return
            
        try:
            from src.features.cheminformatics.services.substructure_matcher import find_substructure_matches
            matches = find_substructure_matches(main_window.molecule, query)
            
            if not matches:
                self.status_label.setText("No matches found.")
                main_window.viewer_3d.set_selected(set())
                main_window.viewer_2d.set_selected(set())
            else:
                self.status_label.setText(f"Found {len(matches)} match(es)! (Selecting all regions).")
                
                # Flatten the matches into a set of atom indices
                flat_indices = set()
                for m in matches:
                    flat_indices.update(m)
                
                main_window.viewer_3d.set_selected(flat_indices)
                main_window.viewer_2d.set_selected(flat_indices)
                
                # Zoom / Focus on the selected atoms in the 3D viewer
                if hasattr(main_window.viewer_3d, 'focus_on_atoms'):
                    main_window.viewer_3d.focus_on_atoms(flat_indices)
                    
        except Exception as e:
            self.status_label.setText(f"Error parsing or querying: {e}")
