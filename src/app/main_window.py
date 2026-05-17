"""
Main Window — Central application window coordinating all components.
"""

import traceback
import time
from src.shared.qt_compat import *

# Performance optimization imports
from src.core.performance import ParallelFileLoader, get_profiler, profile_operation

_DEBUG = False

from src.shared.ui.theme import (
    get_stylesheet, COLORS,
    ThemeMode, apply_initial_theme, set_theme as _apply_theme,
    save_theme_preference, theme_signals, current_mode,
)
from src.features.control_panel.ui.input_panel import InputPanel
from src.app.enhanced_plugin_interface import EnhancedPluginManagerWidget
from src.plugins.enhanced_plugin_manager import EnhancedPluginManager
from src.plugins.plugin_config import PluginConfigManager
from src.core.domain.models.bond import BondType
from src.features.sketcher_2d.ui.sketcher_widget import SketcherWidget
from src.features.visualization_3d.services.docking_pose_service import DockingPoseService
from src.features.ui.docking_pose_dialog import DockingPoseDialog

# ── Extracted modules ────────────────────────────────────────────
from src.app.conversion_worker import ConversionWorker
from src.app import menu_bar as _menu_bar
from src.app import file_operations as _file_ops
from src.app import chemistry_actions as _chem
from src.app import viewer_coordinator as _viewer
from src.app import molecule_controller as _mol_ctrl


class ThemedTabWidget(QTabWidget):
    """
    QTabWidget that paints the empty region beside the tab bar.

    Background
    ----------
    On macOS the QTabBar geometry only covers the actual tab labels
    (e.g. 250 px wide for two short tabs) even when the parent
    QTabWidget is 1000 px wide. Everything from ``tab_bar.right`` to
    ``self.right`` within the tab bar's vertical range is **not**
    painted by any widget -- QTabBar's geometry ends at its own
    width, ``QTabWidget::pane`` starts below the tab bar, and
    QMacStyle ignores ``QTabWidget { background-color }`` stylesheet
    rules for that parent region.

    The result is a stale dark strip that survives every palette,
    stylesheet, autoFillBackground, and polish combination.

    Fix
    ---
    Override ``paintEvent`` and explicitly ``fillRect`` the orphan
    region with the current theme's ``bg_tertiary`` colour using a
    raw QPainter, which bypasses QMacStyle entirely.  A single-pixel
    hairline is drawn along the bottom of that region so the border
    under the tab bar is continuous with the border under each tab.
    """

    def paintEvent(self, event):
        super().paintEvent(event)
        tab_bar = self.tabBar()
        if tab_bar is None or not tab_bar.isVisible():
            return
        x = tab_bar.x() + tab_bar.width()
        w = self.width() - x
        h = tab_bar.height()
        if w <= 0 or h <= 0:
            return
        try:
            from src.shared.ui.theme import COLORS
            fill = QColor(COLORS.get('bg_tertiary', '#ECECEC'))
            line = QColor(COLORS.get('border',      '#D5D5D5'))
            p = QPainter(self)
            p.fillRect(x, tab_bar.y(), w, h, fill)
            p.setPen(QPen(line, 1))
            p.drawLine(x, tab_bar.y() + h - 1,
                       self.width(), tab_bar.y() + h - 1)
            p.end()
        except Exception:
            pass


class MainWindow(QMainWindow):
    """
    Main application window.

    Layout:
    +---------------------------------------------+
    | Menu Bar                                     |
    +---------------------------------------------+
    | Python Console (collapsible)                 |
    +-----------+---------------------------------+
    |           |  [3D View]  [2D View]  (tabs)   |
    | Input     |                                  |
    | Panel     |    Molecular Viewer              |
    |           |                                  |
    +-----------+---------------------------------+
    | Status Bar                                   |
    +---------------------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.molecule = None
        self._worker = None
        self._thread = None
        self._undo_stack = []       # Molecule.clone() snapshots for Ctrl+Z (max 10)
        self._UNDO_LIMIT = 10
        self._com_radius = 0.5
        self._centroid_radius = 0.4
        self.sketcher_window = None

        self.setWindowTitle("PyChem -- Molecular Viewer and Cheminformatics Software")
        self.setMinimumSize(1100, 700)
        self.resize(1280, 800)
        self.setAcceptDrops(True)

        # Enable smooth resizing
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent, False)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, False)

        # Apply the saved theme preference (defaults to System). This
        # installs the stylesheet on QApplication so every child widget
        # inherits it. Must run BEFORE the menu/central widget so they
        # pick up the right colours on first paint.
        apply_initial_theme()

        # Keep menu check marks and left-panel styles in sync when the
        # OS colour scheme changes (System mode only — triggered by
        # theme_signals()) or when the user picks Light/Dark manually.
        theme_signals().theme_changed.connect(self._on_theme_changed)

        # Initialize plugin system first
        self._init_plugin_system()

        self._init_menu_bar()
        self._init_central_widget()
        self._init_status_bar()
        self._connect_signals()

    # ── Menu bar (delegated) ─────────────────────────────────────

    def _init_menu_bar(self):
        """Create menu bar."""
        _menu_bar.build_menu_bar(self)

    # ── Plugin system ────────────────────────────────────────────

    def _init_plugin_system(self):
        """Initialize the enhanced plugin system."""
        try:
            # Initialize plugin configuration manager
            self.plugin_config_manager = PluginConfigManager()
            
            # Initialize enhanced plugin manager
            self.plugin_manager = EnhancedPluginManager()
            
            # Set up plugin API integration
            from src.plugins.utils.integration import PluginIntegrationAPI
            plugin_api = PluginIntegrationAPI(self)
            self.plugin_manager.set_api(plugin_api)
            
            # Connect to plugin events to refresh UI automatically
            self.plugin_manager.add_event_handler('plugin_loaded', lambda e: self._refresh_plugin_tabs())
            self.plugin_manager.add_event_handler('plugin_unloaded', lambda e: self._refresh_plugin_tabs())
            
            # Initialize enhanced plugin interface
            self.enhanced_plugin_interface = EnhancedPluginManagerWidget(self.plugin_manager)
            
            # Use QTimer from qt_compat
            QTimer.singleShot(100, self._delayed_init_plugin_system)
        except Exception as e:
            traceback.print_exc()
            print(f"Error preparing enhanced plugin system: {e}")
            self.plugin_manager = None
            self.enhanced_plugin_interface = None

    def _delayed_init_plugin_system(self):
        """Perform the heavy plugin discovery after the main event loop has started."""
        try:
            if self.plugin_manager:
                # Discover all plugins (bundled + user)
                self.plugin_manager.discover_all_plugins()
                
                # Auto-load enabled plugins
                self._auto_load_plugins()
                
                # Refresh UI
                self._refresh_plugin_list()
                self._refresh_plugin_tabs()
        except Exception as e:
            print(f"Error initializing enhanced plugin system: {e}")
    
    def _auto_load_plugins(self):
        """Auto-load plugins that are configured to auto-load."""
        try:
            all_plugins = self.plugin_manager.get_all_plugins()
            plugin_configs = self.plugin_config_manager.get_all_plugin_configs()
            
            # Get load order based on dependencies
            plugin_names = list(all_plugins.keys())
            load_order = self.plugin_config_manager.get_load_order(plugin_names)
            
            loaded_count = 0
            for plugin_name in load_order:
                # Check if plugin should auto-load
                if self.plugin_config_manager.should_auto_load(plugin_name):
                    # Check if plugin is enabled
                    if self.plugin_config_manager.is_plugin_enabled(plugin_name):
                        if not self.plugin_manager.is_plugin_loaded(plugin_name):
                            if self.plugin_manager.load_plugin(plugin_name):
                                loaded_count += 1
            
            if loaded_count > 0:
                self.status_bar.showMessage(f"Auto-loaded {loaded_count} plugins", 3000)
                
        except Exception as e:
            print(f"Error auto-loading plugins: {e}")

    # ── Central widget ───────────────────────────────────────────

    def _init_central_widget(self):
        """Create central widget with console on top, tabbed viewer below."""
        central = QWidget()
        self.setCentralWidget(central)

        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Vertical splitter: console row on top, content below
        v_splitter = QSplitter(Qt.Orientation.Vertical)
        v_splitter.setHandleWidth(3)

        # Console row: 75% console + 25% right button panel (PyMOL-style)
        console_row = QSplitter(Qt.Orientation.Horizontal)
        console_row.setHandleWidth(2)

        # Python console (left 75% normally, but now takes full row since right panel removed)
        self.console = PythonConsole()
        console_row.addWidget(self.console)

        # The python console now takes all horizontal space in this splitter.
        console_row.setStretchFactor(0, 1)

        v_splitter.addWidget(console_row)

        # Horizontal splitter: input panel + tabbed viewer (bottom)
        h_splitter = QSplitter(Qt.Orientation.Horizontal)
        h_splitter.setHandleWidth(2)

        # Input panel (left). Its objectName is 'leftPanel' and the
        # #leftPanel rule in the global stylesheet handles the fill
        # and the hairline right-hand border — no inline override.
        #
        # NOTE: the panel has setFixedWidth(300) inside its __init__.
        # We must ensure the h_splitter initial sizes match that so
        # the first paint does not clip the sidebar before the fixed
        # width is honoured on the second layout pass.
        self.input_panel = InputPanel()
        h_splitter.addWidget(self.input_panel)

        # Tabbed viewer area (right).  Use ThemedTabWidget which
        # overrides paintEvent to fill the "orphan" strip beside the
        # tab bar with the theme's bg_tertiary colour -- see the
        # ThemedTabWidget docstring for the full explanation of why
        # QTabWidget/QTabBar stylesheet rules cannot cover that region
        # on macOS.
        self.viewer_tabs = ThemedTabWidget()
        self.viewer_tabs.setTabPosition(QTabWidget.TabPosition.North)
        self.viewer_tabs.setElideMode(Qt.TextElideMode.ElideNone)
        self.viewer_tabs.setUsesScrollButtons(False)
        self._viewer_tab_bar = self.viewer_tabs.tabBar()
        if self._viewer_tab_bar is not None:
            self._viewer_tab_bar.setExpanding(False)
            self._viewer_tab_bar.setUsesScrollButtons(False)
            # Fill the empty area behind the tabs so nothing native
            # bleeds through on theme swap.
            self._viewer_tab_bar.setAutoFillBackground(True)
            self._viewer_tab_bar.setDrawBase(False)

        self.viewer_3d = MolViewer3D()
        self.viewer_2d = MolViewer2D()

        self.viewer_tabs.addTab(self.viewer_3d, "3D View")
        self.viewer_tabs.addTab(self.viewer_2d, "2D View")

        # Apply theme-specific inline stylesheet to the viewer tabs so
        # the empty strip beside the tabs and the pane track theme
        # swaps without relying on the global cascade (which Qt's
        # native QTabBar paint on macOS ignores).
        self._apply_viewer_tabs_theme()
        try:
            theme_signals().theme_changed.connect(self._apply_viewer_tabs_theme)
        except Exception:
            pass

        # Plugin dock widget (right side) instead of tab bar
        self._init_plugin_dock()

        h_splitter.addWidget(self.viewer_tabs)
        h_splitter.setStretchFactor(0, 0)
        h_splitter.setStretchFactor(1, 1)
        # Match the initial sizes to InputPanel's fixed 300 px so the
        # sidebar renders at its final width on the very first paint.
        h_splitter.setSizes([300, 980])

        v_splitter.addWidget(h_splitter)
        v_splitter.setStretchFactor(0, 0)
        v_splitter.setStretchFactor(1, 1)

        main_layout.addWidget(v_splitter)

        # Set viewer references for console selection commands
        self.console.set_viewer(self.viewer_3d, self.viewer_2d)

    # ── Plugin dock ──────────────────────────────────────────────

    def _init_plugin_dock(self):
        """
        Create a dockable sidebar for plugin management.

        Styling comes from the global theme — QListWidget, QPushButton
        and QTextEdit are already styled there. Only the dock header
        label and container background need inline styling, and those
        are reapplied on theme change via ``_apply_plugin_dock_theme``.
        """
        self.plugin_dock = QDockWidget("Plugins", self)
        self.plugin_dock.setAllowedAreas(
            Qt.DockWidgetArea.RightDockWidgetArea
            | Qt.DockWidgetArea.LeftDockWidgetArea
        )
        self.plugin_dock.setMinimumWidth(220)

        self._plugin_dock_container = QWidget()
        dock_layout = QVBoxLayout(self._plugin_dock_container)
        dock_layout.setContentsMargins(12, 14, 12, 14)
        dock_layout.setSpacing(10)

        self._plugin_dock_header = QLabel("Installed plugins")
        self._plugin_dock_header.setObjectName("labelSection")
        dock_layout.addWidget(self._plugin_dock_header)

        # Plugin list — global QListWidget rule styles it.
        self.plugin_list = QListWidget()
        dock_layout.addWidget(self.plugin_list)

        # Plugin action buttons — #btnSecondary uses global rules.
        btn_row = QHBoxLayout()
        btn_row.setSpacing(8)

        load_btn = QPushButton("Load")
        load_btn.setObjectName("btnSecondary")
        load_btn.setToolTip("Load selected plugin")
        load_btn.clicked.connect(self._load_selected_plugin)
        btn_row.addWidget(load_btn)

        unload_btn = QPushButton("Unload")
        unload_btn.setObjectName("btnSecondary")
        unload_btn.setToolTip("Unload selected plugin")
        unload_btn.clicked.connect(self._unload_selected_plugin)
        btn_row.addWidget(unload_btn)

        dock_layout.addLayout(btn_row)

        details_label = QLabel("Details")
        details_label.setObjectName("labelSection")
        dock_layout.addWidget(details_label)

        self.plugin_detail = QTextEdit()
        self.plugin_detail.setReadOnly(True)
        self.plugin_detail.setMaximumHeight(140)
        self.plugin_detail.setPlaceholderText(
            "Select a plugin to see details."
        )
        dock_layout.addWidget(self.plugin_detail)

        self._apply_plugin_dock_theme()
        try:
            theme_signals().theme_changed.connect(self._apply_plugin_dock_theme)
        except Exception:
            pass

        self.plugin_dock.setWidget(self._plugin_dock_container)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.plugin_dock)
        self.plugin_dock.hide()  # Hidden by default

        # Populate plugin list
        self._refresh_plugin_list()

    def _refresh_plugin_list(self):
        """Refresh plugin list in dock widget."""
        self.plugin_list.clear()
        if self.plugin_manager:
            try:
                plugins = self.plugin_manager.get_all_plugins()
                loaded = self.plugin_manager.get_loaded_plugins()
                for name, meta in plugins.items():
                    item = QListWidgetItem(name)
                    if name in loaded:
                        item.setForeground(QColor("green"))
                        # Add badge for user plugins
                        if self.plugin_manager.is_user_plugin(name):
                            item.setText(f"{name} (User)")
                    self.plugin_list.addItem(item)
            except Exception as e:
                print(f"Error refreshing plugin list: {e}")

    def _on_plugin_selected(self, row):
        """Show details for selected plugin."""
        if row < 0 or not self.plugin_manager:
            self.plugin_detail.clear()
            return
        try:
            plugins = self.plugin_manager.get_all_plugins()
            names = list(plugins.keys())
            if row < len(names):
                name = names[row]
                meta = plugins[name]
                info = meta.info
                loaded = self.plugin_manager.get_loaded_plugins()
                details = f"Name: {name}\n"
                details += f"Status: {'Loaded' if name in loaded else 'Not loaded'}\n"
                details += f"Version: {info.version}\n"
                details += f"Type: {info.plugin_type.value if hasattr(info.plugin_type, 'value') else info.plugin_type}\n"
                details += f"Description: {info.description}\n"
                details += f"Author: {info.author}\n"
                if self.plugin_manager.is_user_plugin(name):
                    details += f"Source: User Plugin\n"
                else:
                    details += f"Source: Bundled Plugin\n"
                self.plugin_detail.setPlainText(details)
        except Exception:
            pass

    def _load_selected_plugin(self):
        """Load currently selected plugin."""
        row = self.plugin_list.currentRow()
        if row < 0 or not self.plugin_manager:
            return
        try:
            plugins = self.plugin_manager.get_all_plugins()
            names = list(plugins.keys())
            if row < len(names):
                plugin_name = names[row].replace(" (User)", "")
                if self.plugin_manager.load_plugin(plugin_name):
                    self._refresh_plugin_list()
                    self._refresh_plugin_tabs()
                    self.status_bar.showMessage(f"Plugin '{plugin_name}' loaded")
        except Exception as e:
            self.status_bar.showMessage(f"Error loading plugin: {e}")

    def _unload_selected_plugin(self):
        """Unload currently selected plugin."""
        row = self.plugin_list.currentRow()
        if row < 0 or not self.plugin_manager:
            return
        try:
            plugins = self.plugin_manager.get_all_plugins()
            names = list(plugins.keys())
            if row < len(names):
                plugin_name = names[row].replace(" (User)", "")
                if self.plugin_manager.unload_plugin(plugin_name):
                    self._refresh_plugin_list()
                    self._refresh_plugin_tabs()
                    self._on_plugin_selected(row) # Refresh details text immediately
                    self.status_bar.showMessage(f"Plugin '{plugin_name}' unloaded")
        except Exception as e:
            self.status_bar.showMessage(f"Error unloading plugin: {e}")

    def _toggle_plugin_dock(self):
        """Toggle visibility of the plugin dock widget."""
        if self.plugin_dock.isVisible():
            self.plugin_dock.hide()
        else:
            self._refresh_plugin_list()
            self.plugin_dock.show()

    def _refresh_plugins(self):
        """Discover and refresh the plugin list."""
        if self.plugin_manager:
            self.plugin_manager.discover_all_plugins()
            self._refresh_plugin_list()
            self.status_bar.showMessage("Plugins refreshed")

    def _load_all_plugins(self):
        """Load all discovered plugins."""
        if self.plugin_manager:
            plugins = self.plugin_manager.get_all_plugins()
            loaded_count = 0
            for name in plugins:
                if self.plugin_manager.load_plugin(name):
                    loaded_count += 1
            self._refresh_plugin_list()
            self._refresh_plugin_tabs()
            self.status_bar.showMessage(f"Loaded {loaded_count} plugins")

    def _unload_all_plugins(self):
        """Unload all loaded plugins."""
        if self.plugin_manager:
            loaded = self.plugin_manager.get_loaded_plugins()
            names = list(loaded.keys())
            unloaded_count = 0
            for name in names:
                if self.plugin_manager.unload_plugin(name):
                    unloaded_count += 1
            self._refresh_plugin_list()
            self._refresh_plugin_tabs()
            self.status_bar.showMessage(f"Unloaded {unloaded_count} plugins")

    def _refresh_plugin_tabs(self):
        """Update the plugin tabs to match loaded plugins."""
        # Remove all tabs from index 2 onwards (plugin tabs)
        # 0: 3D View, 1: 2D View
        while self.viewer_tabs.count() > 2:
            self.viewer_tabs.removeTab(2)
        
        # Recreate tabs for all loaded plugins
        self._create_plugin_tabs()

    # ── Status bar ───────────────────────────────────────────────

    def _init_status_bar(self):
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready -- Enter a SMILES string to begin")

    # ── Signal wiring (delegated) ────────────────────────────────

    def _connect_signals(self):
        """Connect all widget signals."""
        _mol_ctrl.connect_signals(self)

    # ── Molecule controller delegates ────────────────────────────

    def _convert_smiles(self, smiles):
        _mol_ctrl.convert_smiles(self, smiles)

    def _on_sketcher_import(self, smiles):
        """Handle SMILES imported from the 2D sketcher."""
        if not smiles:
            return
        # Set the SMILES in the input panel for visibility
        self.input_panel.smiles_input.setText(smiles)
        # Start conversion
        self._convert_smiles(smiles)
        # Switch to 3D View to show the result
        self.viewer_tabs.setCurrentWidget(self.viewer_3d)

    def _open_sketcher(self):
        """Open the 2D Sketcher in an independent window."""
        if not self.sketcher_window:
            from src.features.sketcher_2d.ui.sketcher_window import SketcherWindow
            self.sketcher_window = SketcherWindow(self)
            self.sketcher_window.sketcher.molecule_imported.connect(self._on_sketcher_import)
        
        self.sketcher_window.show()
        self.sketcher_window.raise_()
        self.sketcher_window.activateWindow()

    def _on_conversion_done(self, molecule, error):
        _mol_ctrl.on_conversion_done(self, molecule, error)

    def _set_molecule(self, molecule):
        _mol_ctrl.set_molecule(self, molecule)

    def _delete_selected_atoms(self, atom_indices=None):
        _mol_ctrl.delete_selected_atoms(self, atom_indices)

    def _undo_delete(self):
        _mol_ctrl.undo_delete(self)

    def _close_molecule(self):
        _mol_ctrl.close_molecule(self)

    def _deselect_all(self):
        _mol_ctrl.deselect_all(self)

    def _copy_as_image(self):
        """Copy current view as image to clipboard."""
        _file_ops.copy_as_image(self)

    def _show_substructure_dialog(self):
        _mol_ctrl.show_substructure_dialog(self)

    def _notify_plugins_molecule_changed(self):
        """Notify all loaded plugins about molecule changes."""
        if self.plugin_manager and self.molecule:
            self.plugin_manager.set_current_molecule(self.molecule)

    # ── File operations delegates ────────────────────────────────

    def _open_smiles_file(self):
        _file_ops.open_smiles_file(self)

    def _import_structure_file(self, filepath=None):
        _file_ops.import_structure_file(self, filepath)

    def _export_sdf(self):
        _file_ops.export_sdf(self)

    def _export_mol2(self):
        _file_ops.export_mol2(self)

    def _export_image(self, dpi, white_bg):
        _file_ops.export_image(self, dpi, white_bg)

    def _print_views(self):
        _file_ops.print_views(self)

    def _print_preview(self):
        _file_ops.print_preview(self)

    # ── Theme ────────────────────────────────────────────────────

    def _apply_plugin_dock_theme(self):
        """Re-apply theme-dependent inline styling on the plugin dock."""
        if hasattr(self, '_plugin_dock_container'):
            self._plugin_dock_container.setStyleSheet(
                f"background-color: {COLORS['bg_secondary']};"
            )

    def _apply_viewer_tabs_theme(self):
        """
        Paint the viewer QTabWidget's tab bar, pane, and empty strip
        beside the tabs from the current COLORS so theme swaps cannot
        leave a dark native strip on the right.

        macOS's native QTabBar paint reads the widget's QPalette, NOT
        the stylesheet cascade. Inline stylesheets alone are enough on
        programmatic grab() but on a live window the native paint
        still wins. We therefore set BOTH the palette and the inline
        stylesheet on QTabWidget and QTabBar, plus the backing
        autoFillBackground flag, so every painter path picks up the
        right colour no matter which one Qt chooses to use.
        """
        if not hasattr(self, 'viewer_tabs'):
            return
        bg_pane = COLORS['bg_primary']
        bg_bar = COLORS['bg_tertiary']
        bg_tab = COLORS['bg_tertiary']
        bg_tab_sel = COLORS['bg_primary']
        border = COLORS['border']
        text = COLORS['text_primary']
        text_dim = COLORS['text_secondary']
        accent = COLORS['accent']
        hover = COLORS['bg_hover']

        # 1. Palette — this is what macOS native QTabBar actually reads.
        bar_color = QColor(bg_bar)
        text_color = QColor(text)
        pal = self.viewer_tabs.palette()
        pal.setColor(QPalette.ColorRole.Window, bar_color)
        pal.setColor(QPalette.ColorRole.Base, bar_color)
        pal.setColor(QPalette.ColorRole.Button, bar_color)
        pal.setColor(QPalette.ColorRole.WindowText, text_color)
        pal.setColor(QPalette.ColorRole.ButtonText, text_color)
        self.viewer_tabs.setPalette(pal)
        self.viewer_tabs.setAutoFillBackground(True)

        if self._viewer_tab_bar is not None:
            pal2 = self._viewer_tab_bar.palette()
            pal2.setColor(QPalette.ColorRole.Window, bar_color)
            pal2.setColor(QPalette.ColorRole.Base, bar_color)
            pal2.setColor(QPalette.ColorRole.Button, bar_color)
            pal2.setColor(QPalette.ColorRole.WindowText, text_color)
            pal2.setColor(QPalette.ColorRole.ButtonText, text_color)
            self._viewer_tab_bar.setPalette(pal2)
            self._viewer_tab_bar.setAutoFillBackground(True)

        # 2. Stylesheet — this is what applies to the tab shapes and
        #    pseudo-states (hover / selected). The palette covers the
        #    empty strip the stylesheet can't reach reliably.
        self.viewer_tabs.setStyleSheet(
            f"QTabWidget {{ background-color: {bg_bar}; }}"
            f"QTabWidget::pane {{"
            f"  border: 1px solid {border};"
            f"  background-color: {bg_pane};"
            f"  top: 0;"
            f"}}"
            f"QTabBar {{"
            f"  background-color: {bg_bar};"
            f"  border: none;"
            f"}}"
            f"QTabBar::tab {{"
            f"  background-color: {bg_tab};"
            f"  color: {text_dim};"
            f"  padding: 9px 22px;"
            f"  min-width: 80px;"
            f"  border: none;"
            f"  border-right: 1px solid {border};"
            f"  border-bottom: 1px solid {border};"
            f"  font-size: 12px;"
            f"  font-weight: 500;"
            f"}}"
            f"QTabBar::tab:selected {{"
            f"  background-color: {bg_tab_sel};"
            f"  color: {text};"
            f"  border-bottom: 2px solid {accent};"
            f"}}"
            f"QTabBar::tab:hover:!selected {{"
            f"  color: {text};"
            f"  background-color: {hover};"
            f"}}"
        )

        # 3. Force a repaint.
        if self._viewer_tab_bar is not None:
            self._viewer_tab_bar.style().unpolish(self._viewer_tab_bar)
            self._viewer_tab_bar.style().polish(self._viewer_tab_bar)
            self._viewer_tab_bar.update()
        self.viewer_tabs.style().unpolish(self.viewer_tabs)
        self.viewer_tabs.style().polish(self.viewer_tabs)
        self.viewer_tabs.update()

    def _set_theme(self, mode):
        """User picked a theme from View → Theme menu."""
        _apply_theme(mode)
        save_theme_preference(mode)

    def _on_theme_changed(self):
        """
        Called after any theme swap (user action or OS change while
        SYSTEM is selected). Keeps the menu radio group in sync with
        the active mode.
        """
        actions = getattr(self, '_theme_actions', None)
        if not actions:
            return
        active = current_mode()
        for mode, action in actions.items():
            try:
                action.setChecked(mode == active)
            except Exception:
                pass

    def _update_recent_files_menu(self):
        _file_ops.update_recent_files_menu(self)

    def _open_recent_file(self, filepath):
        _file_ops.open_recent_file(self, filepath)

    def _add_recent_file(self, filepath):
        _file_ops.add_recent_file(self, filepath)

    def dragEnterEvent(self, event):
        _file_ops.handle_drag_enter(self, event)

    def dropEvent(self, event):
        _file_ops.handle_drop(self, event)

    # ── Chemistry action delegates ───────────────────────────────

    def _optimize_geometry(self, checked=False, method=None):
        _chem.optimize_geometry(self, checked, method)

    def _compute_charges(self, checked=False, method=None):
        _chem.compute_charges(self, checked, method)

    def _open_descriptor_calculator(self):
        _chem.open_descriptor_calculator(self)

    def _perceive_aromaticity_action(self):
        _chem.perceive_aromaticity_action(self)

    def _generate_smiles(self):
        _chem.generate_smiles(self)

    def _copy_smiles_to_clipboard(self):
        _chem.copy_smiles_to_clipboard(self)

    # ── Lipophilicity Actions ─────────────────────────────────────

    def _calculate_logp_action(self):
        _chem.calculate_logp_action(self)

    def _show_lipophilic_contributions(self):
        _chem.show_lipophilic_contributions(self)

    def _color_by_lipophilicity(self):
        _chem.color_by_lipophilicity(self)

    # ── Viewer coordinator delegates ─────────────────────────────

    def _toggle_hydrogens(self):
        _viewer.toggle_hydrogens(self)

    def _toggle_labels(self):
        _viewer.toggle_labels(self)

    def _toggle_sidechains(self):
        _viewer.toggle_sidechains(self)

    def _toggle_sasa(self, checked):
        _viewer.toggle_sasa(self, checked)

    def _toggle_sasa_selected_only(self, checked):
        _viewer.toggle_sasa_selected_only(self, checked)

    def _change_render_mode(self, index):
        _viewer.change_render_mode(self, index)

    def _change_3d_bg_color(self):
        _viewer.change_3d_bg_color(self)

    def _change_2d_bg_color(self):
        _viewer.change_2d_bg_color(self)

    def _toggle_com_sphere(self):
        _viewer.toggle_com_sphere(self)

    def _toggle_centroid_sphere(self):
        _viewer.toggle_centroid_sphere(self)

    def _show_color_dialog(self):
        _viewer.show_color_dialog(self)

    def _basic_color_toggle(self):
        _viewer.basic_color_toggle(self)

    def _show_protein_color_dialog(self):
        _viewer.show_protein_color_dialog(self)

    def _change_residue_colors(self):
        """Change residue label settings for residues around ligand."""
        from src.features.ui.residue_label_settings_dialog import show_residue_label_settings_dialog
        
        # Show comprehensive residue label settings dialog
        settings = show_residue_label_settings_dialog(self)
        
        if not settings:
            return
            
        self.viewer_3d.residue_label_settings = settings
        
        if not hasattr(self.viewer_3d, 'labeled_residues') or not self.viewer_3d.labeled_residues:
            self.viewer_3d.update()
            return

        color = settings.get('color', Qt.black)
        count = 0
        for rs in self.viewer_3d.labeled_residues:
            self.viewer_3d.labeled_residues[rs] = color
            count += 1
            
        self.viewer_3d.update()
        self.status_bar.showMessage(f"Updated settings for {count} residue labels", 3000)

    def _update_atom_colors(self, colors):
        _viewer.update_atom_colors(self, colors)

    def _add_com_sphere(self):
        _viewer.add_com_sphere(self)

    def _add_centroid_sphere(self):
        _viewer.add_centroid_sphere(self)

    def _add_custom_sphere(self):
        _viewer.add_custom_sphere(self)

    def _create_dummy_atom_sphere(self, position, radius, color, label):
        _viewer.create_dummy_atom_sphere(self, position, radius, color, label)

    def _clear_all_spheres(self):
        _viewer.clear_all_spheres(self)

    # ── Plugin management ────────────────────────────────────────

    def _show_plugin_manager(self):
        """Show enhanced plugin manager dialog."""
        if self.enhanced_plugin_interface:
            self.enhanced_plugin_interface.show()
        else:
            QMessageBox.warning(self, "Plugins", "Enhanced plugin system not available")

    def _show_installed_plugins_dialog(self):
        """Show the installed plugins selection dialog."""
        try:
            from src.app.installed_plugins_dialog import show_installed_plugins_dialog

            result, active_plugins = show_installed_plugins_dialog(
                self.plugin_manager,
                self
            )

            if result:
                print(f"Plugins selected: {active_plugins}")
                self.status_bar.showMessage(f"Selected {len(active_plugins)} plugins for toolbar")

        except Exception as e:
            print(f"Error showing installed plugins dialog: {e}")
            import traceback
            traceback.print_exc()
            QMessageBox.warning(self, "Error", f"Could not open installed plugins dialog: {str(e)}")

    def _refresh_plugins(self):
        """Refresh enhanced plugin system."""
        if self.plugin_manager:
            try:
                self.plugin_manager.discover_all_plugins()
                self._auto_load_plugins()
                self._refresh_plugin_tabs()
                QMessageBox.information(self, "Plugins", "Enhanced plugin system refreshed successfully")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to refresh plugins:\n{str(e)}")
        else:
            QMessageBox.warning(self, "Plugins", "Enhanced plugin system not available")

    def _load_all_plugins(self):
        """Load all available plugins."""
        if self.plugin_manager:
            try:
                plugins = self.plugin_manager.get_all_plugins()
                loaded_count = 0

                for plugin_name in plugins.keys():
                    if not self.plugin_manager.is_plugin_loaded(plugin_name):
                        if self.plugin_manager.load_plugin(plugin_name):
                            loaded_count += 1

                self._refresh_plugin_tabs()
                QMessageBox.information(self, "Plugins", f"Loaded {loaded_count} plugins")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load plugins:\n{str(e)}")
        else:
            QMessageBox.warning(self, "Plugins", "Enhanced plugin system not available")

    def _unload_all_plugins(self):
        """Unload all loaded plugins."""
        if self.plugin_manager:
            try:
                plugins = self.plugin_manager.get_all_plugins()
                unloaded_count = 0

                for plugin_name in plugins.keys():
                    if self.plugin_manager.is_plugin_loaded(plugin_name):
                        if self.plugin_manager.unload_plugin(plugin_name):
                            unloaded_count += 1

                self._refresh_plugin_tabs()
                QMessageBox.information(self, "Plugins", f"Unloaded {unloaded_count} plugins")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to unload plugins:\n{str(e)}")
        else:
            QMessageBox.warning(self, "Plugins", "Enhanced plugin system not available")

    def _refresh_plugin_tabs(self):
        """Refresh plugin tabs, preserving the core viewer tabs."""
        if self.plugin_manager and hasattr(self, 'viewer_tabs'):
            # Preserve the first 2 core tabs: 3D View, 2D View (Sketcher moved to separate window)
            while self.viewer_tabs.count() > 2:
                self.viewer_tabs.removeTab(2)
            self._create_plugin_tabs()

    def _create_plugin_tabs(self):
        """Create tabs for all loaded plugins."""
        if not self.plugin_manager:
            return
            
        loaded_plugins = self.plugin_manager.get_loaded_plugins()
        for plugin_name, plugin in loaded_plugins.items():
            try:
                widget = self.plugin_manager.get_plugin_widget(plugin_name)
                if hasattr(widget, 'get_widget'):
                    # If it has get_widget, use it (standard PluginWidget)
                    tab_widget = widget.get_widget()
                elif hasattr(widget, 'widget') and widget.widget is not widget and isinstance(widget.widget, QWidget):
                    # If it has a .widget attribute that is a different object and is a QWidget
                    tab_widget = widget.widget
                else:
                    # Otherwise assume the object itself is the widget
                    tab_widget = widget
                
                self.viewer_tabs.addTab(tab_widget, plugin_name)
            except Exception as e:
                print(f"Error creating tab for plugin {plugin_name}: {e}")


    # ── About ────────────────────────────────────────────────────

    def _show_about(self):
        QMessageBox.about(
            self, "About PyChem",
            "<h2>PyChem</h2>"
            "<p>Version 1.0.0</p>"
            "<p>A molecular viewer and cheminformatics software with:</p>"
            "<ul>"
            "<li>Full-spec SMILES parser</li>"
            "<li>2D and 3D structure visualization</li>"
            "<li>Distance geometry 3D generation</li>"
            "<li>Geometry optimization</li>"
            "<li>Gasteiger-Marsili partial charges</li>"
            "<li>MOL/SDF/MOL2 import and export</li>"
            "<li>Embedded Python console</li>"
            "</ul>"
            "<p>Built entirely without external cheminformatics libraries.</p>"
            "<p>Developed by Dr. Vijay Masand, Mr. Gaurav Masand, and Mr. Krish Masand.</p>"
            "<p>Support us by giving it a Star on GitHub to help others find it.</p>"
            "<p>We need financial support, please consider donating at <a href='https://buymeacoffee.com/vijaymasand'>https://buymeacoffee.com/vijaymasand</a>.</p>"
        )

    # ── Window events ────────────────────────────────────────────

    def resizeEvent(self, event):
        """Handle window resize events safely."""
        try:
            super().resizeEvent(event)

            if hasattr(self, 'viewer_3d') and self.viewer_3d:
                self.viewer_3d.update()
            if hasattr(self, 'viewer_2d') and self.viewer_2d:
                self.viewer_2d.update()

        except Exception as e:
            print(f"Resize error: {e}")
            pass
    def _open_docking_pose_view(self):
        """Analyze and display molecular docking pose."""
        mol = self.molecule
        if not mol or not mol.atoms:
            QMessageBox.warning(self, "Docking Pose", "No molecule loaded.")
            return

        service = DockingPoseService(mol)
        ligands = service.find_ligands()
        
        if not ligands:
            QMessageBox.warning(self, "Docking Pose", "No suitable ligand fragments detected.")
            return
            
        dialog = DockingPoseDialog(ligands, self)
        
        while True:
            result = dialog.exec_()
            if result == 100: # Picking requested
                self.status_bar.showMessage("Click on an atom in the 3D View to select its fragment...")
                # Create a local event loop to wait for atom_clicked
                from src.shared.qt_compat import QEventLoop
                loop = QEventLoop()
                picked_atom = [-1]
                
                def on_clicked(idx):
                    picked_atom[0] = idx
                    loop.quit()
                
                self.viewer_3d.atom_clicked.connect(on_clicked)
                loop.exec_()
                self.viewer_3d.atom_clicked.disconnect(on_clicked)
                
                if picked_atom[0] != -1:
                    found = dialog.select_fragment_containing(picked_atom[0])
                    if not found:
                        self.status_bar.showMessage("Selected atom is not part of any detected fragment.")
                    else:
                        self.status_bar.showMessage(f"Selected fragment containing atom {picked_atom[0]}")
                dialog.show()
                continue
            elif result == 0: # Cancel
                return
            else: # Apply
                break

        config = dialog.get_config()
        ligand_indices = config['ligand_indices']
        dist = config['distance']
        
        # Handle quick action buttons
        if config['label_nearby']:
            # Label residues within 5.0 Å of ligand
            if hasattr(self, 'console') and self.console:
                self.console._label_residues('within 5.0 ligand')
            return
        elif config['clear_labels']:
            # Clear all residue labels
            if hasattr(self, 'console') and self.console:
                self.console._clear_residue_labels()
            return
        elif config['zoom_to_ligand']:
            # Zoom to ligand with 7 Å surrounding area
            self._zoom_to_ligand_area(ligand_indices, 7.0)
            return
        
        nearby_res_seqs = service.find_nearby_residues(ligand_indices, dist)
        interactions = service.detect_interactions(ligand_indices, nearby_res_seqs)
        
        # Update viewer state
        self.viewer_3d.show_ligands_in_cartoon = True  # Ensure ligands are visible
        self.viewer_3d.visible_sidechains = nearby_res_seqs
        
        # Interaction colors
        interaction_colors = {
            "Hydrogen Bond": "#00FF00",
            "Salt Bridge": "#FF00FF",
            "Hydrophobic": "#FFFF00",
            "Pi-Stacking": "#00FFFF"
        }
        
        lines = []
        for inter in interactions:
            if inter.type == "Hydrogen Bond" and not config['show_hbonds']: continue
            if inter.type == "Salt Bridge" and not config['show_salt']: continue
            if inter.type == "Hydrophobic" and not config['show_hydro']: continue
            
            color = interaction_colors.get(inter.type, "#FFFFFF")
            lines.append((inter.atom1_idx, inter.atom2_idx, inter.type, color))
        
        self.viewer_3d.interaction_lines = lines
        
        # Custom modes and labels
        cam = {}
        labels = {}
        processed_residues = set()
        
        # Ligand atoms as ball-and-stick
        for idx in ligand_indices:
            cam[idx] = 'ball_and_stick'
        
        # Residue atoms as sticks and labels
        for atom in mol.atoms:
            if atom.res_seq in nearby_res_seqs:
                pdb_name = getattr(atom, 'pdb_name', '').strip()
                if pdb_name not in ('CA', 'C', 'N', 'O'):
                    cam[atom.index] = 'stick'
                
                # Disable automatic residue labeling to reduce visual clutter
                # if pdb_name == 'CA' and atom.res_seq not in processed_residues:
                #     res_label = f"{getattr(atom, 'res_name', 'UNK')}{atom.res_seq}"
                #     labels[atom.index] = res_label
                #     processed_residues.add(atom.res_seq)
        
        self.viewer_3d.custom_atom_modes = cam
        self.viewer_3d.labels = labels
        self.viewer_3d.update()
        
        if config['save_report']:
            self._save_docking_report(service, interactions)

    def _zoom_to_ligand_area(self, ligand_indices, radius_angstroms=7.0):
        """Zoom to ligand with specified radius in angstroms."""
        if not self.viewer_3d.molecule or not ligand_indices:
            return
        
        # Use the existing focus_on_atoms method which handles centering correctly, 
        # passing the radius to include the surrounding area
        self.viewer_3d.focus_on_atoms(ligand_indices, padding_angstroms=radius_angstroms)
        self.status_bar.showMessage(f"Zoomed to ligand with {radius_angstroms} Å surrounding area", 3000)

    def _save_docking_report(self, service, interactions):
        """Save interaction report to CSV."""
        path, _ = QFileDialog.getSaveFileName(self, "Save Docking Report", "docking_report.csv", "CSV Files (*.csv)")
        if path:
            try:
                csv_content = service.generate_report_csv(interactions)
                with open(path, 'w') as f:
                    f.write(csv_content)
                QMessageBox.information(self, "Report Saved", f"Report saved to {path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save report: {e}")
