"""
Molecule Controller — Mediates services <-> viewers, undo/redo,
selection sync, conversion orchestration.

All functions receive the MainWindow instance as ``window``.
"""

import time

from src.shared.qt_compat import (
    QApplication, QMessageBox, QThread,
)
from src.core.performance import get_profiler, profile_operation
from src.features.ui.substructure_dialog import SubstructureDialog
from src.app import viewer_coordinator as _viewer

_DEBUG = False


# ── Signal wiring ────────────────────────────────────────────────

def connect_signals(window):
    """Connect all widget signals."""
    window.input_panel.convert_requested.connect(window._convert_smiles)
    window.input_panel.export_sdf_requested.connect(window._export_sdf)
    window.input_panel.export_mol2_requested.connect(window._export_mol2)
    window.input_panel.export_image_requested.connect(window._export_image)

    # 2D Sketcher integration
    if hasattr(window, 'sketcher_2d'):
        window.sketcher_2d.molecule_imported.connect(window._on_sketcher_import)

    # Radius sliders
    window.input_panel.sphere_scale_changed.connect(
        lambda v: setattr(window.viewer_3d, 'sphere_scale', v) or window.viewer_3d.update())
    window.input_panel.stick_scale_changed.connect(
        lambda v: setattr(window.viewer_3d, 'stick_scale', v) or window.viewer_3d.update())
    window.input_panel.line_scale_changed.connect(
        lambda v: setattr(window.viewer_3d, 'line_scale', v) or window.viewer_3d.update())

    # View options
    window.input_panel.show_h_check.toggled.connect(window._toggle_hydrogens)
    window.input_panel.show_labels_check.toggled.connect(window._toggle_labels)
    window.input_panel.label_color_btn.clicked.connect(lambda: _viewer.change_label_color(window))
    window.input_panel.clear_labels_btn.clicked.connect(lambda: _viewer.clear_labels(window))
    window.input_panel.show_sidechains_check.toggled.connect(window._toggle_sidechains)
    window.input_panel.show_sasa_check.toggled.connect(window._toggle_sasa)
    window.input_panel.sasa_selected_only_check.toggled.connect(window._toggle_sasa_selected_only)
    window.input_panel.render_combo.currentIndexChanged.connect(window._change_render_mode)

    # ---------- Selection synchronization (3D <-> 2D) ----------
    window._syncing_selection = False

    def _sync_from_3d(selected):
        """Forward 3D selection to 2D viewer."""
        if window._syncing_selection:
            return
        window._syncing_selection = True
        window.viewer_2d.selected_atoms = set(selected)
        window.viewer_2d.update()
        window._syncing_selection = False

    def _sync_from_2d(selected):
        """Forward 2D selection to 3D viewer."""
        if window._syncing_selection:
            return
        window._syncing_selection = True
        window.viewer_3d.selected_atoms = set(selected)
        window.viewer_3d.update()
        window._syncing_selection = False

    window.viewer_3d.selection_changed.connect(_sync_from_3d)
    window.viewer_2d.selection_changed.connect(_sync_from_2d)

    # ---------- Delete requests from either viewer ----------
    window.viewer_3d.delete_requested.connect(window._delete_selected_atoms)
    window.viewer_2d.delete_requested.connect(window._delete_selected_atoms)


# ── SMILES conversion ───────────────────────────────────────────

def convert_smiles(window, smiles):
    """Convert SMILES to 3D structure."""
    from src.app.conversion_worker import ConversionWorker

    window.status_bar.showMessage(f"Converting: {smiles}")
    window.input_panel.set_progress(5)
    window.input_panel.convert_btn.setEnabled(False)

    window._thread = QThread()
    window._worker = ConversionWorker(smiles)
    window._worker.moveToThread(window._thread)

    window._thread.started.connect(window._worker.run)
    window._worker.progress.connect(window.input_panel.set_progress)
    window._worker.finished.connect(window._on_conversion_done)
    window._worker.finished.connect(window._thread.quit)
    window._worker.finished.connect(window._worker.deleteLater)
    window._thread.finished.connect(window._thread.deleteLater)

    window._thread.start()


def on_conversion_done(window, molecule, error):
    """Handle conversion result."""
    window.input_panel.convert_btn.setEnabled(True)
    window.input_panel.set_progress(0)

    if error:
        window.status_bar.showMessage(f"Error: {error}")
        QMessageBox.critical(window, "Conversion Error",
                            f"Failed to convert SMILES:\n\n{error}")
        return

    set_molecule(window, molecule)
    window.status_bar.showMessage(
        f"Converted: {molecule.molecular_formula()} -- "
        f"{len(molecule.atoms)} atoms, {len(molecule.bonds)} bonds")


# ── Set molecule ─────────────────────────────────────────────────

def set_molecule(window, molecule):
    """Set molecule across all viewers and panels."""

    if _DEBUG:
        print("[DEBUG] _set_molecule called: {} atoms, is_protein={}".format(
            len(molecule.atoms), molecule.properties.get('is_protein', False)))

    set_start_time = time.time()

    # Always compute Center of Mass (fast O(N))
    with get_profiler().time_operation("center_of_mass"):
        try:
            from src.features.cheminformatics.services.spatial_properties import compute_center_of_mass
            compute_center_of_mass(molecule)
        except Exception as e:
            if _DEBUG:
                print(f"[DEBUG] Skipping COM calculation: {e}")

    window.molecule = molecule
    if _DEBUG:
        print("[DEBUG] Setting 3D viewer molecule...")
    with get_profiler().time_operation("3d_viewer_update"):
        window.viewer_3d.set_molecule(molecule)
    if _DEBUG:
        print("[DEBUG] 3D viewer updated")

    # Skip 2D generation for large proteins and PDB files
    source_format = molecule.properties.get('source_format', 'unknown')
    is_pdb = source_format in ('pdb', 'ent')

    should_skip_2d = False
    is_protein = molecule.properties.get('is_protein', False)
    num_atoms = len(molecule.atoms)

    if is_pdb and num_atoms >= 700:
        should_skip_2d = True
    elif is_protein and num_atoms > 1000:
        should_skip_2d = True

    if should_skip_2d:
        if _DEBUG:
            print(f"[DEBUG] Skipping 2D layout: is_pdb={is_pdb}, num_atoms={num_atoms}")
        window.viewer_2d.clear()
        window.viewer_2d.show_protein_placeholder = True
        window.viewer_2d.update()
    else:
        if _DEBUG:
            print("[DEBUG] Generating 2D coordinates...")
        with get_profiler().time_operation("2d_viewer_update"):
            window.viewer_2d.set_molecule(molecule)
        window.viewer_2d.show_protein_placeholder = False

    if _DEBUG:
        print("[DEBUG] Updating console and input panel...")
    with get_profiler().time_operation("console_panel_update"):
        window.console.set_molecule(molecule)
        window.input_panel.update_molecule_info(molecule)
        window.input_panel.enable_tools(True)

    # Notify plugins of molecule change
    if _DEBUG:
        print("[DEBUG] Notifying plugins...")
    with get_profiler().time_operation("plugin_notification"):
        notify_plugins_molecule_changed(window)

    for action in window.opt_menu_actions:
        action.setEnabled(True)
    for action in window.chg_menu_actions:
        action.setEnabled(True)

    # Enable new feature buttons in menu
    window.color_action.setEnabled(True)
    window.protein_color_action.setEnabled(True)
    window.aromaticity_action.setEnabled(True)
    if hasattr(window, 'docking_pose_action'):
        window.docking_pose_action.setEnabled(True)
    for action in window.tools_submenu_actions:
        action.setEnabled(True)
    
    # Enable Lipophilicity actions
    if hasattr(window, '_lipo_menu_actions'):
        for action in window._lipo_menu_actions:
            action.setEnabled(True)
    
    # Enable Copy as Image action
    if hasattr(window, '_copy_image_action'):
        window._copy_image_action.setEnabled(True)


# ── Delete / Undo / Close ────────────────────────────────────────

def delete_selected_atoms(window, atom_indices=None):
    """Delete the given atom indices from the molecule."""
    if not window.molecule:
        return

    indices = set(atom_indices) if atom_indices else set(window.viewer_3d.selected_atoms)
    if not indices:
        window.status_bar.showMessage("Nothing selected to delete")
        return

    # Snapshot for undo
    window._undo_stack.append(window.molecule.clone())
    if len(window._undo_stack) > window._UNDO_LIMIT:
        window._undo_stack.pop(0)
    window._undo_action.setEnabled(True)

    # Remove atoms
    removed = window.molecule.remove_atoms(indices)
    window.status_bar.showMessage(
        f"Deleted {removed} atom(s) \u2014 {len(window.molecule.atoms)} remaining  "
        f"(Ctrl+Z to undo)")

    # Refresh viewers
    if len(window.molecule.atoms) == 0:
        close_molecule(window)
        return

    window.viewer_3d.selected_atoms.clear()
    window.viewer_2d.selected_atoms.clear()
    window.viewer_3d.set_molecule(window.molecule)
    window.viewer_2d.set_molecule(window.molecule)
    window.input_panel.update_molecule_info(window.molecule)


def undo_delete(window):
    """Restore the molecule to its state before the last delete."""
    if not window._undo_stack:
        window.status_bar.showMessage("Nothing to undo")
        return

    window.molecule = window._undo_stack.pop()

    window.viewer_3d.selected_atoms.clear()
    window.viewer_2d.selected_atoms.clear()
    window.viewer_3d.set_molecule(window.molecule)
    window.viewer_2d.set_molecule(window.molecule)
    window.console.set_molecule(window.molecule)
    window.input_panel.update_molecule_info(window.molecule)
    if not window._undo_stack:
        window._undo_action.setEnabled(False)

    window.status_bar.showMessage(
        f"Undo: restored {len(window.molecule.atoms)} atoms, "
        f"{len(window.molecule.bonds)} bonds")


def close_molecule(window):
    """Clear the current molecule without closing the application window."""
    window.molecule = None
    window._undo_stack.clear()
    window._undo_action.setEnabled(False)

    window.viewer_3d.selected_atoms.clear()
    window.viewer_2d.selected_atoms.clear()
    window.viewer_3d.clear()
    window.viewer_2d.clear()
    
    # Disable Copy as Image action
    if hasattr(window, '_copy_image_action'):
        window._copy_image_action.setEnabled(False)

    if hasattr(window, 'console') and window.console:
        window.console.set_molecule(None)

    if hasattr(window, 'input_panel') and window.input_panel:
        window.input_panel.enable_tools(False)

    window.status_bar.showMessage("Ready \u2014 Enter a SMILES string to begin")


def deselect_all(window):
    """Clear atom selection in both viewers."""
    window.viewer_3d.selected_atoms.clear()
    window.viewer_2d.selected_atoms.clear()
    window.viewer_3d.update()
    window.viewer_2d.update()
    window.status_bar.showMessage("Selection cleared")


def show_substructure_dialog(window):
    """Show the SMILES/SMARTS substructure matcher dialog."""
    if not hasattr(window, '_substructure_dialog') or not window._substructure_dialog:
        window._substructure_dialog = SubstructureDialog(window)
    window._substructure_dialog.show()
    window._substructure_dialog.raise_()
    window._substructure_dialog.activateWindow()


# ── Plugins ──────────────────────────────────────────────────────

def notify_plugins_molecule_changed(window):
    """Notify all plugins that the molecule has changed.

    Goes directly to ``window.plugin_manager`` (set in MainWindow.__init__).
    The previous code referenced ``window.plugin_interface.plugin_manager``,
    but that attribute name was never set on MainWindow — the plugin manager
    lives on ``window.plugin_manager`` and the UI widget on
    ``window.enhanced_plugin_interface``.
    """
    plugin_manager = getattr(window, 'plugin_manager', None)
    if plugin_manager is not None:
        plugin_manager.set_current_molecule(window.molecule)
