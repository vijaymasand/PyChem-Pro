# -*- coding: utf-8 -*-
from src.shared.qt_compat import *
from .sketcher_widget import SketcherWidget

class SketcherWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("2D Sketcher")
        self.resize(1000, 700)
        
        self.sketcher = SketcherWidget(self)
        self.setCentralWidget(self.sketcher)
        
        # Connect the imported signal to the main window if possible
        # This will be handled in MainWindow when opening this window.
        
    def closeEvent(self, event):
        # We might want to warn about unsaved changes if we had a proper save system
        super().closeEvent(event)
