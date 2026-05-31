"""
File Operations — Open, save, export, drag-and-drop, and recent-file handling.

All functions receive the MainWindow instance as ``window``.
"""

import os
import time

from src.shared.qt_compat import (
    QApplication, QFileDialog, QMessageBox, QAction, QSettings, Qt, QDialog,
)
from src.core.performance import ParallelFileLoader, get_profiler

_DEBUG = False


# ── Open / Import ─────────────────────────────────────────────────

def open_smiles_file(window):
    filepath, _ = QFileDialog.getOpenFileName(
        window, "Open SMILES File", "",
        "Text Files (*.txt *.smi);;All Files (*)")
    if filepath:
        try:
            with open(filepath, 'r') as f:
                content = f.read().strip()
            smiles = content.split()[0]
            window.input_panel.smiles_input.setText(smiles)
            window._convert_smiles(smiles)
            add_recent_file(window, filepath)
        except Exception as e:
            QMessageBox.critical(window, "Error", f"Failed to read file:\n{e}")


def import_structure_file(window, filepath=None):
    """Import a molecule from MOL, SDF, MOL2, or PDB file."""
    if not filepath:
        filepath, _ = QFileDialog.getOpenFileName(
            window, "Import Structure File", "",
            "All Structure Files (*.mol *.sdf *.mol2 *.pdb *.ent);;"
            "PDB Files (*.pdb *.ent);;"
            "MOL Files (*.mol);;"
            "SDF Files (*.sdf);;"
            "MOL2 Files (*.mol2);;"
            "All Files (*)")
    if not filepath:
        # print("[DEBUG] No file selected")  # Commented out for reduced verbosity
        return

    # print(f"[DEBUG] Selected file: {filepath}")  # Commented out for reduced verbosity

    file_size = os.path.getsize(filepath) / 1024  # KB
    # print(f"[DEBUG] File size: {file_size:.1f} KB")  # Commented out for reduced verbosity

    if file_size > 100:
        window.status_bar.showMessage(f"Importing large file ({file_size:.1f}KB): {filepath}")
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
    else:
        window.status_bar.showMessage(f"Importing: {filepath}")

    QApplication.processEvents()

    load_start_time = time.time()

    try:
        loader = ParallelFileLoader(parallel_threshold_kb=100.0, use_parallel=False)

        ext = filepath.lower().rsplit('.', 1)[-1] if '.' in filepath else ''
        # print(f"[DEBUG] File extension: {ext}")  # Commented out for reduced verbosity

        with get_profiler().time_operation("file_load"):
            mol = loader.load_file(filepath)

        load_time = time.time() - load_start_time
        # print(f"[DEBUG] File loaded in {load_time:.3f}s")  # Commented out for reduced verbosity

        if mol is None:
            raise ValueError("Failed to parse molecule - no atoms loaded")

        mol.properties['source_format'] = ext if ext in ('pdb', 'ent') else 'unknown'

        # print(f"[DEBUG] Molecule loaded: {len(mol.atoms)} atoms, {len(mol.bonds)} bonds")  # Commented out for reduced verbosity
        # print(f"[DEBUG] Molecule name: {mol.name}")  # Commented out for reduced verbosity
        # print(f"[DEBUG] Properties: {mol.properties}")  # Commented out for reduced verbosity

        num_atoms = len(mol.atoms)
        if num_atoms <= 150:
            # print("[DEBUG] Assigning hybridization...")  # Commented out for reduced verbosity
            try:
                mol.assign_hybridization()
                mol.assign_sybyl_types()

                from src.features.smiles_parser.rules.aromaticity import perceive_aromaticity
                perceive_aromaticity(mol)
            except Exception as e:
                # print(f"[DEBUG] Molecule typing/aromaticity failed: {e}")  # Commented out for reduced verbosity
                pass

            # print("[DEBUG] Computing charges...")  # Commented out for reduced verbosity
            has_charges = any(a.partial_charge != 0 for a in mol.atoms)
            if not has_charges:
                try:
                    from src.features.cheminformatics.electrostatics.gasteiger import compute_gasteiger_charges
                    compute_gasteiger_charges(mol)
                except Exception as e:
                    # print(f"[DEBUG] charge computation failed: {e}")  # Commented out for reduced verbosity
                    pass
        else:
            # print("[DEBUG] Skipping hybridization/charges for large structure to ensure fast load speed")  # Commented out for reduced verbosity
            pass

        # print("[DEBUG] Calling _set_molecule...")  # Commented out for reduced verbosity
        window._set_molecule(mol)

        is_protein = mol.properties.get('is_protein', False)
        info = f"Imported: {mol.name or filepath} -- {len(mol.atoms)} atoms, {len(mol.bonds)} bonds"
        if is_protein:
            info += " [PROTEIN]"
        window.status_bar.showMessage(info)
        add_recent_file(window, filepath)

    except Exception as e:
        window.status_bar.showMessage(f"Import error: {e}")
        QMessageBox.critical(window, "Import Error",
                            f"Failed to import file:\n\n{e}")
    finally:
        QApplication.restoreOverrideCursor()


# ── Export ────────────────────────────────────────────────────────

def export_sdf(window):
    if not window.molecule:
        return
    filepath, _ = QFileDialog.getSaveFileName(
        window, "Save SDF File", "",
        "SDF Files (*.sdf);;MOL Files (*.mol);;All Files (*)")
    if filepath:
        try:
            from src.features.io.exporters.sdf_writer import write_sdf
            write_sdf(window.molecule, filepath)
            window.status_bar.showMessage(f"Saved: {filepath}")
        except Exception as e:
            QMessageBox.critical(window, "Export Error", str(e))


def export_mol2(window):
    if not window.molecule:
        return
    filepath, _ = QFileDialog.getSaveFileName(
        window, "Save MOL2 File", "",
        "MOL2 Files (*.mol2);;All Files (*)")
    if filepath:
        try:
            from src.features.io.exporters.mol2_writer import write_mol2
            write_mol2(window.molecule, filepath)
            window.status_bar.showMessage(f"Saved: {filepath}")
        except Exception as e:
            QMessageBox.critical(window, "Export Error", str(e))


def export_image(window, dpi, white_bg):
    if not window.molecule:
        return

    current_viewer = window.viewer_tabs.currentWidget()

    filepath, _ = QFileDialog.getSaveFileName(
        window, f"Export Image ({dpi} DPI)", "",
        "PNG Image (*.png);;TIFF Image (*.tiff *.tif);;JPEG Image (*.jpg *.jpeg);;BMP Image (*.bmp);;All Files (*)")
    if filepath:
        try:
            success = current_viewer.export_image(filepath, dpi=dpi, bg_white=white_bg)
            if success:
                window.status_bar.showMessage(f"Image saved: {filepath} ({dpi} DPI)")
            else:
                QMessageBox.warning(window, "Export Error", "Failed to save image.")
        except Exception as e:
            QMessageBox.critical(window, "Export Error", str(e))


# ── Print ────────────────────────────────────────────────────────

def _render_print_content(window, painter, page_rect, is_preview=False):
    """
    Render the print content (2D and 3D views) to a QPainter.
    
    Args:
        window: MainWindow instance
        painter: QPainter to render to
        page_rect: QRectF of the printable area
        is_preview: If True, add a watermark indicating preview mode
    """
    from PySide6.QtGui import QImage, QColor, QFont, QPainter
    from PySide6.QtCore import Qt, QRectF
    
    page_w = int(page_rect.width())
    page_h = int(page_rect.height())
    margin = int(min(page_w, page_h) * 0.03)

    # Split page into top (2D) and bottom (3D) halves
    half_h = (page_h - 3 * margin) // 2
    top_rect = QRectF(margin, margin, page_w - 2 * margin, half_h)
    bot_rect = QRectF(margin, margin * 2 + half_h,
                      page_w - 2 * margin, half_h)

    # Render each viewer into its own QImage
    def _render_viewer(viewer, width, height, viewer_name="viewer"):
        # Check if 2D or 3D viewer based on attributes
        # 3D viewer has 'zoom', 'pan_x', 'pan_y'
        # 2D viewer has '_scale', '_offset_x', '_offset_y'
        is_3d = hasattr(viewer, 'zoom') and hasattr(viewer, 'pan_x')
        
        if is_3d:
            # 3D viewer - use the working scaling approach
            # Calculate export scale based on target dimensions vs widget size
            widget_width = viewer.width()
            widget_height = viewer.height()
            scale_x = width / widget_width if widget_width > 0 else 1.0
            scale_y = height / widget_height if widget_height > 0 else 1.0
            export_scale = min(scale_x, scale_y)
            
            img = QImage(int(width), int(height),
                         QImage.Format.Format_ARGB32_Premultiplied)
            img.fill(QColor(255, 255, 255))
            p = QPainter(img)
            try:
                old_bg = getattr(viewer, 'bg_color', None)
                if old_bg is not None:
                    viewer.bg_color = QColor(255, 255, 255)
                
                # Scale zoom and pan by export_scale to maintain the same view
                _saved_zoom = viewer.zoom
                _saved_pan_x = viewer.pan_x
                _saved_pan_y = viewer.pan_y
                try:
                    viewer.zoom *= export_scale
                    viewer.pan_x *= export_scale
                    viewer.pan_y *= export_scale
                    viewer._render(p, int(width), int(height),
                                   is_export=True, export_scale=export_scale)
                finally:
                    viewer.zoom = _saved_zoom
                    viewer.pan_x = _saved_pan_x
                    viewer.pan_y = _saved_pan_y
                
                if old_bg is not None:
                    viewer.bg_color = old_bg
            except Exception as e:
                print(f"[PRINT RENDER ERROR] {viewer_name}: {e}")
                import traceback
                traceback.print_exc()
                p.setPen(QColor(255, 0, 0))
                font = QFont("Arial", 10)
                p.setFont(font)
                p.drawText(10, 20, f"Render error: {str(e)[:50]}")
            finally:
                p.end()
            return img
        else:
            # 2D viewer - render at widget's original size and let Qt scale to fit
            # This preserves font readability better than scaling down the molecule
            widget_width = viewer.width()
            widget_height = viewer.height()
            
            img = QImage(int(widget_width), int(widget_height),
                         QImage.Format.Format_ARGB32_Premultiplied)
            img.fill(QColor(255, 255, 255))
            p = QPainter(img)
            try:
                old_bg = getattr(viewer, 'bg_color', None)
                if old_bg is not None:
                    viewer.bg_color = QColor(255, 255, 255)
                
                if hasattr(viewer, 'coords_2d') and viewer.coords_2d:
                    viewer._render(p, int(widget_width), int(widget_height), is_export=True)
                else:
                    # Draw placeholder for empty 2D view
                    p.setPen(QColor(150, 150, 150))
                    font = QFont("Arial", 12)
                    p.setFont(font)
                    p.drawText(10, int(widget_height/2), "2D view not available")
                
                if old_bg is not None:
                    viewer.bg_color = old_bg
            except Exception as e:
                print(f"[PRINT RENDER ERROR] {viewer_name}: {e}")
                import traceback
                traceback.print_exc()
                p.setPen(QColor(255, 0, 0))
                font = QFont("Arial", 10)
                p.setFont(font)
                p.drawText(10, 20, f"Render error: {str(e)[:50]}")
            finally:
                p.end()
            return img

    # Render 2D
    v2d = window.viewer_2d
    img_2d = _render_viewer(v2d, top_rect.width(), top_rect.height(), "2D")

    # Render 3D
    v3d = window.viewer_3d
    img_3d = _render_viewer(v3d, bot_rect.width(), bot_rect.height(), "3D")

    # Draw headers
    title_font = QFont("Helvetica", 9)
    title_font.setBold(True)
    painter.setFont(title_font)
    painter.setPen(QColor(0, 0, 0))

    mol_name = getattr(window.molecule, 'name', None) or 'Molecule'
    try:
        mol_formula = window.molecule.molecular_formula()
    except:
        mol_formula = "?"
    num_atoms = getattr(window.molecule, 'num_atoms', 0)
    num_bonds = getattr(window.molecule, 'num_bonds', 0)
    
    header_text = (
        f"PyChem -- {mol_name}  |  "
        f"{mol_formula}  |  "
        f"{num_atoms} atoms, "
        f"{num_bonds} bonds"
    )
    painter.drawText(QRectF(margin, margin * 0.25,
                             page_w - 2 * margin, margin * 0.7),
                      int(Qt.AlignmentFlag.AlignLeft |
                          Qt.AlignmentFlag.AlignVCenter),
                      header_text)

    # Draw view labels (2D View / 3D View)
    label_font = QFont("Helvetica", 9)
    label_font.setBold(True)
    painter.setFont(label_font)
    painter.drawText(QRectF(margin, top_rect.y() - margin * 0.6,
                             page_w - 2 * margin, margin * 0.6),
                      int(Qt.AlignmentFlag.AlignLeft),
                      "2D View")
    painter.drawText(QRectF(margin, bot_rect.y() - margin * 0.6,
                             page_w - 2 * margin, margin * 0.6),
                      int(Qt.AlignmentFlag.AlignLeft),
                      "3D View")

    # Draw the two images into their respective halves
    painter.drawImage(top_rect, img_2d)
    painter.drawImage(bot_rect, img_3d)
    
    # Add preview watermark if in preview mode
    if is_preview:
        watermark_font = QFont("Helvetica", 48)
        watermark_font.setBold(True)
        painter.setFont(watermark_font)
        painter.setPen(QColor(200, 200, 200, 100))  # Semi-transparent gray
        painter.drawText(page_rect, int(Qt.AlignmentFlag.AlignCenter), "PREVIEW")


def print_preview(window):
    """
    Show a print preview dialog for the 2D and 3D views.
    
    This allows users to see how the print will look before sending to printer.
    """
    if not window.molecule:
        window.status_bar.showMessage("No molecule to preview")
        return

    try:
        from PySide6.QtPrintSupport import QPrinter, QPrintPreviewDialog
        from PySide6.QtGui import QPainter
    except ImportError:
        QMessageBox.critical(
            window, "Print Preview Error",
            "Qt print support is not available. Please install PySide6 with "
            "print support (PySide6 >= 6.5 normally includes it)."
        )
        return

    # Configure printer with sensible defaults
    printer = QPrinter(QPrinter.PrinterMode.HighResolution)
    printer.setDocName(f"PyChem - {window.molecule.name or 'molecule'}")

    # Create preview dialog
    preview_dialog = QPrintPreviewDialog(printer, window)
    preview_dialog.setWindowTitle("Print Preview - 2D + 3D Views")
    
    # Connect the paint requested signal
    def handle_paint_request(printer_obj):
        from PySide6.QtPrintSupport import QPrinter
        from PySide6.QtGui import QFont, QColor, QPainter
        painter = QPainter(printer_obj)
        try:
            page_rect = printer_obj.pageRect(QPrinter.Unit.DevicePixel)
            _render_print_content(window, painter, page_rect, is_preview=True)
        except Exception as e:
            print(f"[PRINT PREVIEW ERROR] {e}")
            import traceback
            traceback.print_exc()
            # Draw error message on page
            painter.setPen(QColor(255, 0, 0))
            painter.setFont(QFont("Arial", 12))
            painter.drawText(100, 100, f"Preview Error: {str(e)}")
        finally:
            painter.end()
    
    preview_dialog.paintRequested.connect(handle_paint_request)
    preview_dialog.exec()


def print_views(window):
    """
    Print both 2D and 3D views of the current molecule to a single page.

    Opens a native QPrintDialog, then renders:
      * the 2D viewer into the top half of the printable area
      * the 3D viewer into the bottom half of the printable area
    Both views are fitted and centered inside their half, preserving
    aspect ratio, so the molecule is fully visible on the page.
    """
    if not window.molecule:
        window.status_bar.showMessage("No molecule to print")
        return

    try:
        from PySide6.QtPrintSupport import QPrinter, QPrintDialog
        from PySide6.QtGui import QPainter
    except ImportError:
        QMessageBox.critical(
            window, "Print Error",
            "Qt print support is not available. Please install PySide6 with "
            "print support (PySide6 >= 6.5 normally includes it)."
        )
        return

    # Configure printer with sensible defaults
    printer = QPrinter(QPrinter.PrinterMode.HighResolution)
    printer.setDocName(f"PyChem - {window.molecule.name or 'molecule'}")

    dialog = QPrintDialog(printer, window)
    dialog.setWindowTitle("Print 2D + 3D Views")
    if dialog.exec() != QPrintDialog.DialogCode.Accepted:
        return

    try:
        # Get printable area in device pixels
        page_rect = printer.pageRect(QPrinter.Unit.DevicePixel)
        
        # Paint onto the printer
        page_painter = QPainter(printer)
        try:
            _render_print_content(window, page_painter, page_rect, is_preview=False)
        finally:
            page_painter.end()

        window.status_bar.showMessage("Print job sent")
    except Exception as e:
        QMessageBox.critical(window, "Print Error",
                             f"Failed to print:\n{e}")
        window.status_bar.showMessage(f"Print error: {e}")


# ── Drag & Drop ──────────────────────────────────────────────────

def handle_drag_enter(window, event):
    if event.mimeData().hasUrls():
        urls = event.mimeData().urls()
        if urls and urls[0].isLocalFile():
            ext = urls[0].toLocalFile().lower()
            if ext.endswith(('.smi', '.mol', '.sdf', '.mol2', '.pdb', '.ent')):
                event.acceptProposedAction()
                return
    event.ignore()


def handle_drop(window, event):
    urls = event.mimeData().urls()
    if urls and urls[0].isLocalFile():
        filepath = urls[0].toLocalFile()
        ext = filepath.lower()
        if ext.endswith('.smi'):
            try:
                with open(filepath, 'r') as f:
                    smiles = f.read().strip().split()[0]
                window.input_panel.smiles_input.setText(smiles)
                window._convert_smiles(smiles)
                add_recent_file(window, filepath)
            except Exception as e:
                window.status_bar.showMessage(f"Error reading SMILES: {e}")
        else:
            window._import_structure_file(filepath)
        event.acceptProposedAction()


# ── Recent Files ─────────────────────────────────────────────────

def update_recent_files_menu(window):
    """Update recent files context menu."""
    window.recent_menu.clear()
    settings = QSettings("PyChem", "Viewer")
    files = settings.value("recent_files", [])
    if not isinstance(files, list):
        files = []
    for f in files:
        action = QAction(f, window)
        action.triggered.connect(lambda checked=False, path=f: open_recent_file(window, path))
        window.recent_menu.addAction(action)
    if not files:
        action = QAction("No Recent Files", window)
        action.setEnabled(False)
        window.recent_menu.addAction(action)


def open_recent_file(window, filepath):
    if not os.path.exists(filepath):
        QMessageBox.warning(window, "Error", f"File not found:\n{filepath}")
        return
    ext = filepath.lower()
    if ext.endswith('.smi'):
        try:
            with open(filepath, 'r') as f:
                smiles = f.read().strip().split()[0]
            window.input_panel.smiles_input.setText(smiles)
            window._convert_smiles(smiles)
        except Exception as e:
            window.status_bar.showMessage(f"Error reading SMILES: {e}")
    else:
        window._import_structure_file(filepath)


def add_recent_file(window, filepath):
    settings = QSettings("PyChem", "Viewer")
    files = settings.value("recent_files", [])
    if not isinstance(files, list):
        files = []
    if filepath in files:
        files.remove(filepath)
    files.insert(0, filepath)
    if len(files) > 5:
        files = files[:5]
    settings.setValue("recent_files", files)
    update_recent_files_menu(window)


# ── Copy as Image ─────────────────────────────────────────────────

def copy_as_image(window):
    """
    Copy the current molecular view as an image to the clipboard.
    
    Opens a dialog to select:
    - Source: 2D view, 3D view, or Both
    - Format: PNG, JPEG, BMP
    - DPI: Resolution (72-1200)
    - Selection only: Copy only selected atoms (if any)
    
    If nothing is selected, the full molecule is copied automatically.
    The image is placed on the clipboard for pasting into MS Word/PowerPoint.
    """
    if not window.molecule:
        QMessageBox.information(window, "Copy as Image", "No molecule loaded.")
        return
    
    from PySide6.QtGui import QImage, QPainter, QClipboard
    from PySide6.QtCore import Qt, QRect, QSize
    from src.app.copy_image_dialog import CopyImageDialog
    
    # Check for selections in both viewers
    has_sel_2d = bool(window.viewer_2d.selected_atoms) if hasattr(window, 'viewer_2d') else False
    has_sel_3d = bool(window.viewer_3d.selected_atoms) if hasattr(window, 'viewer_3d') else False
    
    # Show dialog
    dialog = CopyImageDialog(window, has_sel_2d, has_sel_3d)
    if dialog.exec() != QDialog.DialogCode.Accepted:
        return
    
    settings = dialog.get_settings()
    source = settings["source"]
    fmt = settings["format"]
    dpi = settings["dpi"]
    selection_only = settings["selection_only"]
    high_quality = settings.get("high_quality", False)
    
    try:
        # Calculate scale factor based on DPI (96 is standard screen DPI)
        scale_factor = dpi / 96.0
        
        if source == "2d":
            image = _render_viewer_to_image(
                window.viewer_2d, scale_factor, selection_only, fmt, high_quality=False
            )
        elif source == "3d":
            image = _render_viewer_to_image(
                window.viewer_3d, scale_factor, selection_only, fmt, high_quality=high_quality
            )
        else:  # both
            image = _render_both_viewers(window, scale_factor, selection_only, fmt, high_quality=high_quality)
        
        if image is None:
            QMessageBox.warning(window, "Copy as Image", "Failed to render image.")
            return
        
        # Set DPI metadata on the image
        image.setDotsPerMeterX(int(dpi / 0.0254))
        image.setDotsPerMeterY(int(dpi / 0.0254))
        
        # Copy to clipboard
        clipboard = QApplication.clipboard()
        clipboard.setImage(image)
        
        # Show success message
        sel_text = " (selection)" if selection_only and (has_sel_2d or has_sel_3d) else ""
        window.status_bar.showMessage(
            f"Copied {source.upper()} view{sel_text} to clipboard ({dpi} DPI, {fmt.upper()})"
        )
        
    except Exception as e:
        QMessageBox.critical(window, "Copy as Image Error", str(e))


def _render_viewer_to_image(viewer, scale_factor, selection_only, fmt, high_quality=False):
    """
    Render a single viewer (2D or 3D) to a QImage.
    
    Uses the viewer's built-in export_image method for correct rendering,
    then crops to the molecule bounding box.
    """
    from PySide6.QtGui import QImage, QColor
    from PySide6.QtCore import Qt
    import io
    
    is_3d = not hasattr(viewer, 'coords_2d')
    SCREEN_MARGIN = 50  # Pixels to add around molecule
    
    # Calculate bounding box of molecule in original screen coordinates
    mol = viewer.molecule if hasattr(viewer, 'molecule') else None
    min_x = min_y = float('inf')
    max_x = max_y = float('-inf')
    
    if not is_3d and viewer.coords_2d:
        # 2D: get atom positions in screen coords using _to_screen transformation
        atom_indices = viewer.selected_atoms if (selection_only and viewer.selected_atoms) else set(viewer.coords_2d.keys())
        for idx in atom_indices:
            if idx in viewer.coords_2d:
                x, y = viewer.coords_2d[idx]
                # Match _to_screen: sx = mx * _scale + _offset_x, sy = -my * _scale + _offset_y
                sx = x * viewer._scale + viewer._offset_x
                sy = -y * viewer._scale + viewer._offset_y
                min_x, max_x = min(min_x, sx), max(max_x, sx)
                min_y, max_y = min(min_y, sy), max(max_y, sy)
                
    elif is_3d and mol and mol.atoms:
        # 3D: project atoms to screen
        atoms = mol.atoms if not (selection_only and hasattr(viewer, 'selected_atoms') and viewer.selected_atoms) else \
                [a for a in mol.atoms if a.index in viewer.selected_atoms]
        cx, cy = viewer.width() / 2, viewer.height() / 2
        zoom = getattr(viewer, 'zoom', 1.0)
        pan_x = getattr(viewer, 'pan_x', 0) or 0
        pan_y = getattr(viewer, 'pan_y', 0) or 0
        for atom in atoms:
            x3d, y3d = getattr(atom, 'x3d', None), getattr(atom, 'y3d', None)
            if x3d is None or y3d is None:
                continue
            sx = cx + (x3d * zoom) + pan_x
            sy = cy - (y3d * zoom) + pan_y
            min_x, max_x = min(min_x, sx), max(max_x, sx)
            min_y, max_y = min(min_y, sy), max(max_y, sy)
    
    # Use full viewer if no atoms found
    if min_x == float('inf'):
        bbox_x, bbox_y, bbox_w, bbox_h = 0, 0, viewer.width(), viewer.height()
    else:
        bbox_x = max(0, min_x - SCREEN_MARGIN)
        bbox_y = max(0, min_y - SCREEN_MARGIN)
        bbox_w = min(viewer.width(), max_x - min_x + 2 * SCREEN_MARGIN)
        bbox_h = min(viewer.height(), max_y - min_y + 2 * SCREEN_MARGIN)
    
    # Create a temp file for export
    import tempfile
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
        temp_path = f.name
    
    try:
        # Use viewer's export to get full image
        dpi = int(96 * scale_factor)
        # print(f"[DEBUG CopyImage] Exporting to {temp_path} at dpi={dpi}")  # Commented out for reduced verbosity
        if hasattr(viewer, 'export_image'):
            import inspect
            sig = inspect.signature(viewer.export_image)
            if 'high_quality' in sig.parameters:
                success = viewer.export_image(temp_path, dpi=dpi, bg_white=True, high_quality=high_quality)
            else:
                success = viewer.export_image(temp_path, dpi=dpi, bg_white=True)
        else:
            success = False
        # print(f"[DEBUG CopyImage] Export success: {success}")  # Commented out for reduced verbosity
        if not success:
            return None
            
        # Load the exported image
        full_image = QImage(temp_path)
        # print(f"[DEBUG CopyImage] Image loaded: {not full_image.isNull()}, size: {full_image.width()}x{full_image.height()}")  # Commented out for reduced verbosity
        if full_image.isNull():
            return None
        
        # Calculate crop region (scaled to exported image coords)
        export_scale = full_image.width() / viewer.width()
        crop_x = int(bbox_x * export_scale)
        crop_y = int(bbox_y * export_scale)
        crop_w = int(bbox_w * export_scale)
        crop_h = int(bbox_h * export_scale)
        
        # print(f"[DEBUG CopyImage] Bounding box: ({bbox_x:.1f}, {bbox_y:.1f}, {bbox_w:.1f}, {bbox_h:.1f})")  # Commented out for reduced verbosity
        # print(f"[DEBUG CopyImage] Export scale: {export_scale:.2f}")  # Commented out for reduced verbosity
        # print(f"[DEBUG CopyImage] Initial crop: ({crop_x}, {crop_y}, {crop_w}, {crop_h})")  # Commented out for reduced verbosity
        
        # Clamp bounds
        crop_x = max(0, min(crop_x, full_image.width() - 1))
        crop_y = max(0, min(crop_y, full_image.height() - 1))
        crop_w = max(10, min(crop_w, full_image.width() - crop_x))
        crop_h = max(10, min(crop_h, full_image.height() - crop_y))
        
        # print(f"[DEBUG CopyImage] Final crop: ({crop_x}, {crop_y}, {crop_w}, {crop_h})")  # Commented out for reduced verbosity

        # For 2D, return full image (cropping is causing issues)
        if not is_3d:
            # print(f"[DEBUG CopyImage] Returning full 2D image")  # Commented out for reduced verbosity
            return full_image
        else:
            return full_image.copy(crop_x, crop_y, crop_w, crop_h)
        
    finally:
        # Cleanup temp file
        try:
            import os
            os.unlink(temp_path)
        except:
            pass


def _render_both_viewers(window, scale_factor, selection_only, fmt, high_quality=False):
    """
    Render both 2D and 3D viewers side by side.
    
    Returns:
        QImage containing both views
    """
    from PySide6.QtGui import QImage, QPainter, QColor
    from PySide6.QtCore import Qt
    
    # Get dimensions
    v2d = window.viewer_2d
    v3d = window.viewer_3d
    
    # Use larger dimensions scaled up
    w = max(v2d.width(), v3d.width())
    h = max(v2d.height(), v3d.height())
    
    img_w = int(w * scale_factor * 2.1)  # 2 viewers side by side with gap
    img_h = int(h * scale_factor)
    
    # Create combined image
    image = QImage(img_w, img_h, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(QColor(255, 255, 255))
    
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    
    try:
        # Render 2D on left
        img_2d = _render_viewer_to_image(v2d, scale_factor, selection_only, fmt, high_quality=False)
        if img_2d:
            painter.drawImage(0, 0, img_2d)
        
        # Render 3D on right
        img_3d = _render_viewer_to_image(v3d, scale_factor, selection_only, fmt, high_quality=high_quality)
        if img_3d:
            x_offset = int(img_w * 0.5)
            painter.drawImage(x_offset, 0, img_3d)
        
    finally:
        painter.end()
    
    return image
