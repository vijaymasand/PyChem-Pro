"""
Alignment Actions — superimpose the open structures in the 3D scene.

The scene (see :mod:`src.app.molecule_scene`) can hold several structures at
once; the *active* object is the reference frame.  "Align Visible → Active"
rigid-body fits every other visible object onto the active one (Kabsch fit with
outlier rejection) and moves those objects' coordinates in place, so the overlay
redraws superimposed.

All functions receive the MainWindow instance as ``window``.  The heavy lifting
lives in the GUI-agnostic :class:`~src.services.alignment.align_service.AlignmentService`
(also exposed as ``pychem.align`` for notebooks/scripts).
"""

from src.shared.qt_compat import QApplication, QMessageBox


def align_visible_to_active(window):
    """Superimpose every other visible scene object onto the active object."""
    scene = getattr(window, 'scene', None)
    if scene is None or scene.is_empty():
        QMessageBox.information(window, "Align Structures",
                                "Load at least two structures to align.")
        return

    active = scene.active_object()
    if active is None:
        QMessageBox.information(window, "Align Structures",
                                "No active object to align onto. Click a row in "
                                "the Objects panel to make it the reference.")
        return

    mobiles = [o for o in scene.visible_objects() if o.id != active.id]
    if not mobiles:
        QMessageBox.information(
            window, "Align Structures",
            "Need at least two visible structures.\n\n"
            f"'{active.name}' is the reference — make one or more other objects "
            "visible to align them onto it.")
        return

    window.status_bar.showMessage(f"Aligning {len(mobiles)} structure(s) onto '{active.name}'...")
    QApplication.processEvents()

    from pychem._bridge import get_registry
    service = get_registry().alignment

    successes = []
    failures = []
    for obj in mobiles:
        try:
            # Moves obj.molecule's coordinates in place onto the reference.
            result = service.align(obj.molecule, active.molecule)
            successes.append((obj.name, result))
        except Exception as exc:  # correspondence too small, no coords, etc.
            failures.append((obj.name, str(exc)))

    # Rebuild the overlay from the (now-moved) molecules and re-centre the view.
    window._refresh_scene(refit=True)

    _report(window, active.name, successes, failures)


def _report(window, ref_name, successes, failures):
    """Show a per-object RMSD summary and update the status bar."""
    lines = [f"Reference: {ref_name}", ""]
    for name, res in successes:
        lines.append(
            f"✓ {name}: RMSD {res.rmsd:.3f} Å  "
            f"({res.n_aligned} atoms, method '{res.method}'"
            + (f", {res.n_rejected} rejected" if res.n_rejected else "")
            + ")")
    for name, err in failures:
        lines.append(f"✗ {name}: {err}")

    if successes:
        best = min(res.rmsd for _, res in successes)
        window.status_bar.showMessage(
            f"Aligned {len(successes)} structure(s) onto '{ref_name}' "
            f"(best RMSD {best:.3f} Å)"
            + (f"; {len(failures)} failed" if failures else ""))
        icon = QMessageBox.Icon.Information
    else:
        window.status_bar.showMessage("Alignment failed for all structures.")
        icon = QMessageBox.Icon.Warning

    box = QMessageBox(window)
    box.setIcon(icon)
    box.setWindowTitle("Alignment Result")
    box.setText("\n".join(lines))
    box.setStandardButtons(QMessageBox.StandardButton.Ok)
    box.exec()
