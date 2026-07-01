"""
Molecule Controller — Mediates services <-> viewers, undo/redo,
selection sync, conversion orchestration, and the multi-object scene.

All functions receive the MainWindow instance as ``window``.

Multi-object model
------------------
``window.scene`` (:class:`~src.app.molecule_scene.MoleculeScene`) holds every
loaded structure.  The 3D viewer is fed a single *merged* display molecule
(built by :meth:`MoleculeScene.build_overlay`) so the existing render pipeline
is reused for the PyMOL-style overlay.  The **active** object is merged *first*,
so its local atom indices equal the merged indices — that keeps ``window.molecule``
(the active object) index-compatible with the 2D viewer and every existing
per-atom feature.  ``window.molecule`` always mirrors the active object.
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

    # ---------- Object panel (multi-object scene) ----------
    op = getattr(window, 'object_panel', None)
    if op is not None:
        op.visibility_toggled.connect(
            lambda oid, v: on_object_visibility(window, oid, v))
        op.all_visibility_toggled.connect(
            lambda v: on_all_visibility(window, v))
        op.active_changed.connect(
            lambda oid: on_object_active(window, oid))
        op.representation_changed.connect(
            lambda oid, r: on_object_representation(window, oid, r))
        op.color_changed.connect(
            lambda oid, c: on_object_color(window, oid, c))
        op.remove_requested.connect(
            lambda oid: on_object_remove(window, oid))

    # ---------- Selection synchronization (3D <-> 2D) ----------
    window._syncing_selection = False

    def _sync_from_3d(selected):
        """Forward 3D selection to 2D viewer (active object's atoms only).

        The active object is merged first, so its merged indices are
        ``0 .. n_active-1`` — identical to the 2D viewer's local indices.
        Atoms belonging to other overlaid objects are simply not shown in 2D.
        """
        if window._syncing_selection:
            return
        window._syncing_selection = True
        active = window.scene.active_object()
        n = len(active.molecule.atoms) if active and active.molecule else 0
        window.viewer_2d.selected_atoms = {i for i in selected if i < n}
        window.viewer_2d.update()
        window._syncing_selection = False

    def _sync_from_2d(selected):
        """Forward 2D selection to 3D viewer (identity for the active object)."""
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
    window._pending_smiles = smiles

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

    name = molecule.name or getattr(window, '_pending_smiles', None) \
        or molecule.molecular_formula() or "molecule"
    add_molecule_object(window, molecule, name=name, source_format='smiles')
    window.status_bar.showMessage(
        f"Converted: {molecule.molecular_formula()} -- "
        f"{len(molecule.atoms)} atoms, {len(molecule.bonds)} bonds")


# ── Scene: add / refresh / active ────────────────────────────────

def add_molecule_object(window, molecule, name=None, source_path=None,
                        source_format='unknown', make_active=True):
    """Add a structure to the scene as a new object and refresh the overlay."""
    with get_profiler().time_operation("center_of_mass"):
        try:
            from src.features.cheminformatics.services.spatial_properties import (
                compute_center_of_mass)
            compute_center_of_mass(molecule)
        except Exception:
            pass

    name = name or molecule.name or molecule.molecular_formula() or "molecule"
    window.scene.add_object(molecule, name, source_path, source_format, make_active)
    refresh_scene(window, refit=True, reset_selection=True)


def set_molecule(window, molecule):
    """Backward-compatible entry: add the molecule as a new active object."""
    add_molecule_object(window, molecule)


def refresh_scene(window, refit=False, reset_selection=False):
    """Rebuild the merged 3D overlay and re-point the 2D/chemistry context.

    Args:
        refit: re-fit the 3D camera to the whole scene (use on load).
        reset_selection: clear the atom selection (use when the visible-object
            set or merge order changes, since merged indices then shift).
    """
    scene = window.scene
    overlay = scene.build_overlay()
    window._scene_atom_origin = overlay.atom_origin
    window._scene_object_ranges = overlay.object_ranges

    active = scene.active_object()
    window.molecule = active.molecule if active else None

    _apply_overlay_3d(window, overlay, refit, reset_selection)
    _apply_active_secondary(window, window.molecule)

    if hasattr(window, 'object_panel'):
        window.object_panel.set_objects(scene.objects, scene.active_id)


def _apply_overlay_3d(window, overlay, refit, reset_selection):
    """Push the merged display molecule + per-object styling into the 3D viewer.

    ``viewer_3d`` is the software/GL wrapper, so the molecule swap goes through
    its ``set_molecule`` (which handles software/GL selection, auto-fit and
    centroid).  On a light refresh the camera is saved and restored so toggling
    visibility / style / colour does not jump the view.
    """
    v3 = window.viewer_3d

    if reset_selection:
        v3.selected_atoms = set()
        window.viewer_2d.selected_atoms = set()

    if overlay.molecule is None:
        v3.clear()
        return

    cam = None
    if not refit:
        cam = (v3.rot_x, v3.rot_y, v3.rot_z, v3.pan_x, v3.pan_y, v3.zoom)

    v3.set_molecule(overlay.molecule)

    if cam is not None:
        v3.rot_x, v3.rot_y, v3.rot_z = cam[0], cam[1], cam[2]
        v3.pan_x, v3.pan_y = cam[3], cam[4]
        v3.zoom = cam[5]

    # Re-apply overlay styling — set_molecule clears colours/labels and chooses
    # a render mode from is_protein.
    v3.labels = {}
    v3.custom_atom_modes = overlay.custom_modes
    v3.software_viewer.custom_atom_colors = overlay.custom_colors
    if hasattr(v3.gl_viewer, 'custom_atom_colors'):
        v3.gl_viewer.custom_atom_colors = overlay.custom_colors
    v3.render_mode = _compute_global_mode(window.scene, overlay)

    if hasattr(v3.software_viewer, '_renderer'):
        v3.software_viewer._renderer.invalidate_cache()

    # The GL renderer (used for large structures) bakes per-object colours and
    # representations into its mesh, so rebuild it when such styling is active.
    # rebuild_mesh() degrades to a plain update() when GL is not in use.
    if (overlay.custom_colors or overlay.custom_modes) and \
            getattr(v3, 'active_viewer', None) is getattr(v3, 'gl_viewer', None):
        v3.rebuild_mesh()
    else:
        v3.update()


def _compute_global_mode(scene, overlay):
    """Global render mode for atoms without a per-object representation.

    Cartoon / ribbon / backbone are whole-chain modes — they apply globally
    (driven by the active object's choice, else any protein present).
    """
    active = scene.active_object()
    if active is not None and active.representation in ('cartoon', 'ribbon', 'backbone'):
        return active.representation
    return 'cartoon' if overlay.any_protein else 'ball_and_stick'


def _apply_active_secondary(window, molecule):
    """Apply the active object to the 2D view, console, panels and plugins."""
    if molecule is None:
        window.viewer_2d.clear()
        if hasattr(window, 'console') and window.console:
            window.console.set_molecule(None)
        if hasattr(window, 'input_panel') and window.input_panel:
            window.input_panel.enable_tools(False)
        enable_molecule_actions(window, False)
        notify_plugins_molecule_changed(window)
        return

    # 2D view — skip layout for large proteins / PDB structures.
    source_format = molecule.properties.get('source_format', 'unknown')
    is_pdb = source_format in ('pdb', 'ent')
    is_protein = molecule.properties.get('is_protein', False)
    num_atoms = len(molecule.atoms)

    should_skip_2d = (is_pdb and num_atoms >= 700) or (is_protein and num_atoms > 1000)

    if should_skip_2d:
        window.viewer_2d.clear()
        window.viewer_2d.show_protein_placeholder = True
        window.viewer_2d.update()
    else:
        with get_profiler().time_operation("2d_viewer_update"):
            window.viewer_2d.set_molecule(molecule)
        window.viewer_2d.show_protein_placeholder = False

    with get_profiler().time_operation("console_panel_update"):
        window.console.set_molecule(molecule)
        window.input_panel.update_molecule_info(molecule)
        window.input_panel.enable_tools(True)

    notify_plugins_molecule_changed(window)
    enable_molecule_actions(window, True)


def enable_molecule_actions(window, enabled):
    """Enable/disable menu actions that require a loaded molecule."""
    for attr in ('opt_menu_actions', 'chg_menu_actions',
                 'tools_submenu_actions', '_lipo_menu_actions'):
        for action in getattr(window, attr, []) or []:
            try:
                action.setEnabled(enabled)
            except Exception:
                pass
    for attr in ('color_action', 'protein_color_action', 'aromaticity_action',
                 'docking_pose_action', '_copy_image_action'):
        action = getattr(window, attr, None)
        if action is not None:
            try:
                action.setEnabled(enabled)
            except Exception:
                pass


# ── Object-panel handlers ────────────────────────────────────────

def on_object_visibility(window, obj_id, visible):
    window.scene.set_visible(obj_id, visible)
    refresh_scene(window, reset_selection=True)


def on_all_visibility(window, visible):
    window.scene.set_all_visible(visible)
    refresh_scene(window, reset_selection=True)


def on_object_active(window, obj_id):
    window.scene.set_active(obj_id)
    # Active object is merged first, so indices shift — reset selection.
    refresh_scene(window, reset_selection=True)


def on_object_representation(window, obj_id, representation):
    window.scene.set_representation(obj_id, representation)
    refresh_scene(window)


def on_object_color(window, obj_id, color):
    window.scene.set_color(obj_id, color)
    refresh_scene(window)


def on_object_remove(window, obj_id):
    window.scene.remove_object(obj_id)
    if window.scene.is_empty():
        _clear_all_views(window)
        window._undo_stack.clear()
        window._undo_action.setEnabled(False)
    else:
        refresh_scene(window, reset_selection=True)
    window.status_bar.showMessage("Object removed")


def set_active_object(window, obj_id):
    on_object_active(window, obj_id)


# ── Delete / Undo / Close ────────────────────────────────────────

def delete_selected_atoms(window, atom_indices=None):
    """Delete selected atoms, routing them back to their source objects."""
    scene = window.scene
    if scene.is_empty():
        return

    indices = set(atom_indices) if atom_indices else set(window.viewer_3d.selected_atoms)
    if not indices:
        window.status_bar.showMessage("Nothing selected to delete")
        return

    origin = window._scene_atom_origin
    per_obj = {}
    for midx in indices:
        if midx in origin:
            oid, local_idx = origin[midx]
            per_obj.setdefault(oid, set()).add(local_idx)

    if not per_obj:
        return

    # Snapshot affected objects for undo.
    snapshot = []
    for oid in per_obj:
        obj = scene.get(oid)
        if obj is not None and obj.molecule is not None:
            snapshot.append((obj, obj.molecule.clone()))
    window._undo_stack.append(snapshot)
    if len(window._undo_stack) > window._UNDO_LIMIT:
        window._undo_stack.pop(0)
    window._undo_action.setEnabled(True)

    total_removed = 0
    for oid, local_idxs in per_obj.items():
        obj = scene.get(oid)
        if obj is not None and obj.molecule is not None:
            total_removed += obj.molecule.remove_atoms(local_idxs)
    # Drop any object that lost all its atoms.
    for obj, _clone in list(snapshot):
        if obj.molecule is not None and len(obj.molecule.atoms) == 0:
            scene.remove_object(obj.id)

    if scene.is_empty():
        _clear_all_views(window)
        window.status_bar.showMessage(
            f"Deleted {total_removed} atom(s) — scene cleared (Ctrl+Z to undo)")
        return

    refresh_scene(window, reset_selection=True)
    window.status_bar.showMessage(
        f"Deleted {total_removed} atom(s) (Ctrl+Z to undo)")


def undo_delete(window):
    """Restore the objects affected by the last delete."""
    if not window._undo_stack:
        window.status_bar.showMessage("Nothing to undo")
        return

    snapshot = window._undo_stack.pop()
    scene = window.scene
    for obj, clone in snapshot:
        obj.molecule = clone
        if obj not in scene.objects:
            scene.objects.append(obj)        # re-insert an emptied object
    if scene.active_id == -1 and scene.objects:
        scene.active_id = scene.objects[0].id

    refresh_scene(window, reset_selection=True)

    if not window._undo_stack:
        window._undo_action.setEnabled(False)
    window.status_bar.showMessage("Undo: restored deleted atoms")


def close_molecule(window):
    """Close the active object (or clear the scene when it is the last one)."""
    scene = window.scene
    active = scene.active_object()
    if active is not None:
        scene.remove_object(active.id)

    window._undo_stack.clear()
    window._undo_action.setEnabled(False)

    if scene.is_empty():
        _clear_all_views(window)
    else:
        refresh_scene(window, refit=False, reset_selection=True)


def _clear_all_views(window):
    """Reset every viewer/panel to the empty state."""
    window.molecule = None
    window._scene_atom_origin = {}
    window._scene_object_ranges = {}

    window.viewer_3d.selected_atoms.clear()
    window.viewer_2d.selected_atoms.clear()
    window.viewer_3d.clear()
    window.viewer_2d.clear()

    if hasattr(window, 'console') and window.console:
        window.console.set_molecule(None)
    if hasattr(window, 'input_panel') and window.input_panel:
        window.input_panel.enable_tools(False)
    if hasattr(window, 'object_panel'):
        window.object_panel.set_objects(window.scene.objects, window.scene.active_id)

    enable_molecule_actions(window, False)
    window.status_bar.showMessage("Ready — Enter a SMILES string to begin")


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
    """Notify all plugins that the active molecule has changed."""
    plugin_manager = getattr(window, 'plugin_manager', None)
    if plugin_manager is not None:
        plugin_manager.set_current_molecule(window.molecule)
