"""
Viewer Coordinator — Synchronises 2D/3D viewer state, view options,
color dialogs, sphere management.

All functions receive the MainWindow instance as ``window``.
"""

from src.shared.qt_compat import (
    QApplication, QColorDialog, QMessageBox, QInputDialog, QPaintEvent,
)

_DEBUG = False


# ── View toggles ─────────────────────────────────────────────────

def toggle_hydrogens(window):
    checked = window.input_panel.show_h_check.isChecked()
    window.viewer_3d.show_hydrogens = checked
    window.viewer_3d.update()


def toggle_labels(window):
    checked = window.input_panel.show_labels_check.isChecked()
    window.viewer_3d.show_labels = checked
    window.viewer_3d.update()


def clear_labels(window):
    """Clear all custom labels and uncheck the show labels box."""
    window.viewer_3d.clear_labels()
    window.input_panel.show_labels_check.setChecked(False)
    window.status_bar.showMessage("Labels cleared")


def change_label_color(window):
    """Open color picker to change label color."""
    from PySide6.QtWidgets import QColorDialog
    current = window.viewer_3d.label_color
    color = QColorDialog.getColor(current, window, "Choose Label Color")
    if color.isValid():
        window.viewer_3d.label_color = color
        window.viewer_3d.update()
        window.status_bar.showMessage(f"Label color changed to {color.name()}")


def toggle_sidechains(window):
    window.viewer_3d.show_sidechains = window.input_panel.show_sidechains_check.isChecked()
    window.viewer_3d.update()


def toggle_sasa(window, checked):
    window.viewer_3d.show_sasa_surface = checked
    if checked and window.molecule:
        has_points = any(
            getattr(a, 'sasa_points', None)
            for a in window.molecule.atoms
        )
        if not has_points:
            window.status_bar.showMessage("Computing SASA surface...")
            QApplication.processEvents()
            try:
                from src.features.cheminformatics.services.spatial_properties import compute_sasa
                compute_sasa(window.molecule)
                window.status_bar.showMessage("SASA surface computed")
            except Exception as e:
                window.status_bar.showMessage(f"SASA computation failed: {e}")
    window.viewer_3d.update()


def toggle_sasa_selected_only(window, checked):
    window.viewer_3d.show_sasa_selected_only = checked
    window.viewer_3d.update()


def change_render_mode(window, index):
    modes = ['ball_and_stick', 'spacefill', 'wireframe', 'cartoon', 'ribbon', 'backbone']
    if 0 <= index < len(modes):
        window.viewer_3d.render_mode = modes[index]
        window.viewer_3d.update()


# ── Background color ─────────────────────────────────────────────

def change_3d_bg_color(window):
    """Open color picker to change 3D viewer background."""
    from PySide6.QtGui import QColor
    current = window.viewer_3d.bg_color
    color = QColorDialog.getColor(
        current, window, "Choose 3D Background Color",
        QColorDialog.ColorDialogOption.ShowAlphaChannel)
    if color.isValid():
        window.viewer_3d.bg_color = color
        window.viewer_3d.update()
        window.status_bar.showMessage(
            f"3D Background: {color.name()}")


def change_2d_bg_color(window):
    """Open color picker to change 2D viewer background."""
    from PySide6.QtGui import QColor
    current = window.viewer_2d.bg_color
    color = QColorDialog.getColor(
        current, window, "Choose 2D Background Color",
        QColorDialog.ColorDialogOption.ShowAlphaChannel)
    if color.isValid():
        window.viewer_2d.bg_color = color
        window.viewer_2d.update()
        window.status_bar.showMessage(
            f"2D Background: {color.name()}")


# ── COM / Centroid toggles (View menu checkboxes) ────────────────

def toggle_com_sphere(window):
    """Toggle COM dummy sphere visibility."""
    if not window.molecule:
        return

    from src.features.visualization_3d.services.dummy_sphere import DummySphereManager
    manager = DummySphereManager(window.molecule)

    for s in manager.get_all_spheres():
        if s.label == 'COM':
            manager.remove_sphere(s.sphere_id)

    if hasattr(window, '_toggle_com_action') and window._toggle_com_action.isChecked():
        manager.create_sphere_at_com(radius=window._com_radius)

    window.viewer_3d.update()


def toggle_centroid_sphere(window):
    """Toggle Centroid dummy sphere visibility."""
    if not window.molecule:
        return

    from src.features.visualization_3d.services.dummy_sphere import DummySphereManager
    manager = DummySphereManager(window.molecule)

    for s in manager.get_all_spheres():
        if s.label == 'Centroid':
            manager.remove_sphere(s.sphere_id)

    if hasattr(window, '_toggle_centroid_action') and window._toggle_centroid_action.isChecked():
        manager.create_sphere_at_centroid(radius=window._centroid_radius)

    window.viewer_3d.update()


# ── Sphere tools (Applications > Tools) ──────────────────────────

def add_com_sphere(window):
    """Add dummy sphere at center of mass with radius selection."""
    if not window.molecule:
        window.status_bar.showMessage("No molecule loaded for COM sphere")
        return

    radius, ok = QInputDialog.getDouble(
        window, "COM Sphere Radius",
        "Enter radius (\u00c5):", window._com_radius, 0.1, 50.0, 2
    )
    if not ok:
        return

    window._com_radius = radius

    try:
        from src.features.visualization_3d.services.dummy_sphere import DummySphereManager
        manager = DummySphereManager(window.molecule)

        for s in manager.get_all_spheres():
            if s.label == 'COM':
                manager.remove_sphere(s.sphere_id)

        sphere_id = manager.create_sphere_at_com(radius=window._com_radius)
        com_pos = manager.calculate_center_of_mass()

        window.viewer_3d.update()
        window.status_bar.showMessage(
            f"COM sphere added at ({com_pos[0]:.2f}, {com_pos[1]:.2f}, {com_pos[2]:.2f}) "
            f"with radius {window._com_radius} \u00c5")

    except Exception as e:
        window.status_bar.showMessage(f"COM sphere error: {e}")


def add_centroid_sphere(window):
    """Add dummy sphere at geometric centroid with radius selection."""
    if not window.molecule:
        window.status_bar.showMessage("No molecule loaded for centroid sphere")
        return

    radius, ok = QInputDialog.getDouble(
        window, "Centroid Sphere Radius",
        "Enter radius (\u00c5):", window._centroid_radius, 0.1, 50.0, 2
    )
    if not ok:
        return

    window._centroid_radius = radius

    try:
        from src.features.visualization_3d.services.dummy_sphere import DummySphereManager
        manager = DummySphereManager(window.molecule)

        for s in manager.get_all_spheres():
            if s.label == 'Centroid':
                manager.remove_sphere(s.sphere_id)

        sphere_id = manager.create_sphere_at_centroid(radius=window._centroid_radius)
        centroid_pos = manager.calculate_geometric_centroid()

        window.viewer_3d.update()
        window.status_bar.showMessage(
            f"Centroid sphere added at ({centroid_pos[0]:.2f}, {centroid_pos[1]:.2f}, "
            f"{centroid_pos[2]:.2f}) with radius {window._centroid_radius} \u00c5")

    except Exception as e:
        window.status_bar.showMessage(f"Centroid sphere error: {e}")


def add_custom_sphere(window):
    """Add multi-layered dummy spheres at custom position using unified dialog."""
    if not window.molecule:
        window.status_bar.showMessage("No molecule loaded for custom sphere")
        return

    try:
        from src.features.ui.custom_sphere_dialog import CustomSphereDialog
        dialog = CustomSphereDialog(window.molecule, window)

        if dialog.exec():
            pos, layers = dialog.get_result()

            from src.features.visualization_3d.services.dummy_sphere import DummySphereManager
            manager = DummySphereManager(window.molecule)

            created_count = 0
            for layer in layers:
                sphere_id = manager.create_sphere_at_position(
                    pos,
                    radius=layer['radius'],
                    color=layer['color'],
                    alpha=layer['alpha'],
                    label='Custom'
                )
                created_count += 1

            if hasattr(window.viewer_3d, 'update'):
                window.viewer_3d.update()

            window.status_bar.showMessage(
                f"Added {created_count} custom sphere shell(s) at "
                f"({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f})")

    except Exception as e:
        window.status_bar.showMessage(f"Custom sphere error: {e}")
        import traceback
        traceback.print_exc()


def create_dummy_atom_sphere(window, position, radius, color, label):
    """Create a dummy atom to represent a sphere in the 3D viewer."""
    try:
        from src.core.domain.models.atom import Atom

        dummy_atom = Atom('Xe')
        dummy_atom.coords = position
        dummy_atom.partial_charge = 0.0

        original_count = len(window.molecule.atoms)
        window.molecule.add_atom(dummy_atom)

        if hasattr(window.viewer_3d, 'set_molecule'):
            window.viewer_3d.set_molecule(window.molecule)
        else:
            window.viewer_3d.update()

        if not hasattr(window, '_dummy_atoms'):
            window._dummy_atoms = []
        window._dummy_atoms.append(dummy_atom)

    except Exception as e:
        print(f"Dummy atom creation error: {e}")


def clear_all_spheres(window):
    """Clear all dummy spheres."""
    if not window.molecule:
        window.status_bar.showMessage("No molecule loaded")
        return

    try:
        from src.features.visualization_3d.services.dummy_sphere import DummySphereManager

        manager = DummySphereManager(window.molecule)
        original_count = len(manager.get_all_spheres())

        if original_count == 0:
            window.status_bar.showMessage("No spheres to clear")
            return

        manager.clear_all_spheres()

        if hasattr(window, '_dummy_atoms'):
            for dummy_atom in window._dummy_atoms:
                if dummy_atom in window.molecule.atoms:
                    window.molecule.atoms.remove(dummy_atom)
            window._dummy_atoms.clear()

        if hasattr(window.viewer_3d, 'set_molecule'):
            window.viewer_3d.set_molecule(window.molecule)
        else:
            window.viewer_3d.update()

        window.status_bar.showMessage(f"Cleared {original_count} sphere(s)")

    except Exception as e:
        window.status_bar.showMessage(f"Clear spheres error: {e}")


# ── Color dialogs ────────────────────────────────────────────────

def show_color_dialog(window):
    """Show PySide6-native color dialog to avoid threading issues."""
    if not window.molecule:
        window.status_bar.showMessage("No molecule loaded for color customization")
        return

    try:
        from src.features.ui.pyside6_color_dialog import show_pyside6_color_dialog, apply_pyside6_colors

        selected_colors = show_pyside6_color_dialog(window)

        if selected_colors:
            apply_pyside6_colors(selected_colors)
            update_atom_colors(window, selected_colors)

            color_count = len(selected_colors)
        else:
            window.status_bar.showMessage("Color selection cancelled")

    except Exception as e:
        import traceback
        traceback.print_exc()
        window.status_bar.showMessage(f"Color dialog error: {str(e)}")


def basic_color_toggle(window):
    """Basic color toggle as last resort."""
    try:
        from src.shared.ui.theme import COLORS

        if not hasattr(window, '_color_toggle_index'):
            window._color_toggle_index = 0

        basic_schemes = [
            {'atom_c': '#55ff7f', 'atom_h': '#d0d0d0', 'atom_o': '#ff0d0d', 'atom_n': '#3050f8'},
            {'atom_c': '#ff6b6b', 'atom_o': '#4ecdc4', 'atom_n': '#45b7d1'},
            {'atom_c': '#2ecc71', 'atom_o': '#e74c3c', 'atom_n': '#3498db'},
        ]

        current_scheme = basic_schemes[window._color_toggle_index % len(basic_schemes)]
        COLORS.update(current_scheme)

        update_atom_colors(window, current_scheme)
        window._color_toggle_index += 1

        scheme_names = ["Default", "Pastel", "Vibrant"]
        current_name = scheme_names[(window._color_toggle_index - 1) % len(scheme_names)]

        window.status_bar.showMessage(f"Applied {current_name} color scheme")

    except Exception as e:
        window.status_bar.showMessage(f"Color customization unavailable: {e}")


def show_protein_color_dialog(window):
    """Show protein color customization dialog."""
    if not window.molecule:
        window.status_bar.showMessage("No molecule loaded for protein color customization")
        return

    try:
        from src.features.visualization_3d.ui.protein_color_dialog import show_protein_color_dialog as _show

        selected_colors = _show(window)

        if selected_colors:
            if hasattr(window.viewer_3d, 'repaint'):
                window.viewer_3d.repaint()

            if hasattr(window.viewer_3d, 'update'):
                window.viewer_3d.update()

            try:
                from src.features.visualization_3d.services.protein_rendering import _cartoon_gen
                _cartoon_gen.invalidate()
            except ImportError:
                pass

            from PySide6.QtCore import QCoreApplication
            QCoreApplication.processEvents()

            window.status_bar.showMessage(
                f"Applied protein colors: Helix={selected_colors.get('ss_helix', 'N/A')}, "
                f"Sheet={selected_colors.get('ss_sheet', 'N/A')}")
        else:
            window.status_bar.showMessage("Protein color selection cancelled")

    except Exception as e:
        import traceback
        traceback.print_exc()
        window.status_bar.showMessage(f"Protein color dialog error: {str(e)}")


def update_atom_colors(window, colors):
    """Update atom colors in 3D viewer."""
    try:
        if not window.molecule:
            return

        for atom in window.molecule.atoms:
            element_symbol = atom.element.symbol.lower()
            color_key = f'atom_{element_symbol}'

            if color_key in colors:
                new_color = colors[color_key]
                atom.element.color = new_color

                from src.shared.ui.theme import COLORS
                COLORS[color_key] = new_color

        from src.shared.ui.theme import COLORS

        sphere_keys = ['sphere_default', 'sphere_com', 'sphere_centroid', 'sphere_custom']
        for key in sphere_keys:
            if key in colors:
                COLORS[key] = colors[key]

        sphere_colors_selected = any(key in colors for key in sphere_keys)

        if window.molecule and sphere_colors_selected:
            try:
                from src.features.visualization_3d.services.dummy_sphere import DummySphere
                import numpy as np

                if not hasattr(window.molecule, 'dummy_spheres'):
                    window.molecule.dummy_spheres = []

                if hasattr(window.molecule, 'atoms') and len(window.molecule.atoms) > 0:
                    positions = [
                        [atom.x, atom.y, atom.z]
                        for atom in window.molecule.atoms if hasattr(atom, 'x')
                    ]
                    if positions:
                        positions = np.array(positions)
                        center = positions.mean(axis=0)

                        if 'sphere_default' in colors:
                            default_sphere = DummySphere(
                                position=tuple(center + np.array([0.5, 0.0, 0.0])),
                                radius=0.5,
                                color=colors['sphere_default'],
                                label="Default"
                            )
                            window.molecule.dummy_spheres.append(default_sphere)

                        if 'sphere_com' in colors:
                            com_sphere = DummySphere(
                                position=tuple(center + np.array([0.0, 0.5, 0.0])),
                                radius=0.5,
                                color=colors['sphere_com'],
                                label="COM"
                            )
                            window.molecule.dummy_spheres.append(com_sphere)

                        if 'sphere_centroid' in colors:
                            centroid_sphere = DummySphere(
                                position=tuple(center + np.array([0.0, -0.5, 0.0])),
                                radius=0.5,
                                color=colors['sphere_centroid'],
                                label="Centroid"
                            )
                            window.molecule.dummy_spheres.append(centroid_sphere)

                        if hasattr(window.viewer_3d, 'set_molecule'):
                            window.viewer_3d.set_molecule(window.molecule)
            except Exception as e:
                pass

        stick_keys = ['stick_default', 'stick_single', 'stick_double', 'stick_triple', 'stick_selected', 'stick_highlight']
        for key in stick_keys:
            if key in colors:
                COLORS[key] = colors[key]

        if hasattr(window.viewer_3d, 'sphere_scale'):
            for key, value in colors.items():
                if hasattr(window.viewer_3d, key):
                    setattr(window.viewer_3d, key, value)

        if hasattr(window.viewer_3d, 'update_atom_colors'):
            window.viewer_3d.update_atom_colors(colors)

        elif hasattr(window.viewer_3d, 'set_colors'):
            window.viewer_3d.set_colors(colors)

        elif hasattr(window.viewer_3d, 'color_atoms'):
            window.viewer_3d.color_atoms(colors)

        elif hasattr(window.viewer_3d, 'update_colors'):
            window.viewer_3d.update_colors(colors)

        else:
            if hasattr(window.viewer_3d, 'COLORS'):
                from src.shared.ui.theme import COLORS as _COLORS
                window.viewer_3d.COLORS.update(_COLORS)

            if hasattr(window.viewer_3d, 'set_molecule'):
                window.viewer_3d.set_molecule(window.molecule)

            if hasattr(window.viewer_3d, '_renderer') and hasattr(window.viewer_3d._renderer, 'invalidate_cache'):
                window.viewer_3d._renderer.invalidate_cache()

            if hasattr(window.viewer_3d, 'update'):
                window.viewer_3d.update()

            if hasattr(window.viewer_3d, 'repaint'):
                window.viewer_3d.repaint()

            if hasattr(window.viewer_3d, 'refresh'):
                window.viewer_3d.refresh()

            if hasattr(window.viewer_3d, 'redraw'):
                window.viewer_3d.redraw()

            if hasattr(window.viewer_3d, 'paintEvent'):
                try:
                    window.viewer_3d.paintEvent(QPaintEvent())
                except Exception:
                    pass

            if hasattr(window.viewer_3d, 'updateGeometry'):
                window.viewer_3d.updateGeometry()

    except Exception as e:
        import traceback
        traceback.print_exc()
