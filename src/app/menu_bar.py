"""
Menu Bar — Construction of all menus and their actions.

All functions receive the MainWindow instance as ``window``.
"""

from src.shared.qt_compat import QAction, QKeySequence


def build_menu_bar(window):
    """Create the full menu bar on *window*."""
    menu_bar = window.menuBar()

    _build_file_menu(window, menu_bar)
    _build_edit_menu(window, menu_bar)
    _build_tools_menu(window, menu_bar)
    _build_plugins_menu(window, menu_bar)
    _build_view_menu(window, menu_bar)
    _build_applications_menu(window, menu_bar)
    _build_support_menu(window, menu_bar)
    _build_help_menu(window, menu_bar)


# ── File ──────────────────────────────────────────────────────────

def _build_file_menu(window, menu_bar):
    file_menu = menu_bar.addMenu("&File")

    open_smiles = QAction("Open &SMILES File...", window)
    open_smiles.setShortcut(QKeySequence.StandardKey.Open)
    open_smiles.triggered.connect(window._open_smiles_file)
    file_menu.addAction(open_smiles)

    import_mol = QAction("&Import MOL/SDF/MOL2...", window)
    import_mol.setShortcut(QKeySequence("Ctrl+I"))
    import_mol.triggered.connect(lambda: window._import_structure_file(None))
    file_menu.addAction(import_mol)

    window.recent_menu = file_menu.addMenu("Open &Recent")
    window._update_recent_files_menu()

    file_menu.addSeparator()

    save_sdf = QAction("Save as &SDF...", window)
    save_sdf.setShortcut(QKeySequence("Ctrl+S"))
    save_sdf.triggered.connect(window._export_sdf)
    file_menu.addAction(save_sdf)

    save_mol2 = QAction("Save as &MOL2...", window)
    save_mol2.setShortcut(QKeySequence("Ctrl+Shift+S"))
    save_mol2.triggered.connect(window._export_mol2)
    file_menu.addAction(save_mol2)

    save_img = QAction("Export &Image...", window)
    save_img.setShortcut(QKeySequence("Ctrl+E"))
    save_img.triggered.connect(lambda: window._export_image(300, True))
    file_menu.addAction(save_img)

    file_menu.addSeparator()

    close_action = QAction("&Close Molecule", window)
    close_action.setShortcut(QKeySequence("Ctrl+W"))
    close_action.triggered.connect(window._close_molecule)
    file_menu.addAction(close_action)
    window._close_action = close_action

    print_preview_action = QAction("Print Pre&view...", window)
    print_preview_action.setShortcut(QKeySequence("Ctrl+Shift+P"))
    print_preview_action.triggered.connect(window._print_preview)
    file_menu.addAction(print_preview_action)
    window._print_preview_action = print_preview_action

    print_action = QAction("&Print...", window)
    print_action.setShortcut(QKeySequence("Ctrl+P"))
    print_action.triggered.connect(window._print_views)
    file_menu.addAction(print_action)
    window._print_action = print_action

    file_menu.addSeparator()

    exit_action = QAction("E&xit", window)
    exit_action.setShortcut(QKeySequence("Ctrl+Q"))
    exit_action.triggered.connect(window.close)
    file_menu.addAction(exit_action)


# ── Edit ──────────────────────────────────────────────────────────

def _build_edit_menu(window, menu_bar):
    edit_menu = menu_bar.addMenu("&Edit")

    window._undo_action = QAction("&Undo Delete", window)
    window._undo_action.setShortcut(QKeySequence("Ctrl+Z"))
    window._undo_action.triggered.connect(window._undo_delete)
    window._undo_action.setEnabled(False)
    edit_menu.addAction(window._undo_action)

    delete_sel_action = QAction("&Delete Selected", window)
    delete_sel_action.setShortcut(QKeySequence("Delete"))
    delete_sel_action.triggered.connect(
        lambda: window._delete_selected_atoms(
            set(window.viewer_3d.selected_atoms) if hasattr(window, 'viewer_3d') else set()))
    edit_menu.addAction(delete_sel_action)

    edit_menu.addSeparator()

    copy_image_action = QAction("&Copy as Image...", window)
    copy_image_action.setShortcut(QKeySequence("Ctrl+Shift+I"))
    copy_image_action.triggered.connect(window._copy_as_image)
    edit_menu.addAction(copy_image_action)
    window._copy_image_action = copy_image_action

    edit_menu.addSeparator()

    deselect_action = QAction("Deselect &All", window)
    deselect_action.setShortcut(QKeySequence("Escape"))
    deselect_action.triggered.connect(window._deselect_all)
    edit_menu.addAction(deselect_action)


# ── Tools ─────────────────────────────────────────────────────────

def _build_tools_menu(window, menu_bar):
    tools_menu = menu_bar.addMenu("&Tools")

    tools_menu.addSeparator()

    auto_rotate = QAction("Toggle &Auto-Rotate", window)
    auto_rotate.setShortcut(QKeySequence("Ctrl+R"))
    auto_rotate.triggered.connect(lambda: window.viewer_3d.toggle_auto_rotate())
    tools_menu.addAction(auto_rotate)

    reset_view = QAction("&Reset View", window)
    reset_view.setShortcut(QKeySequence("Ctrl+0"))
    reset_view.triggered.connect(lambda: window.viewer_3d.reset_view())
    tools_menu.addAction(reset_view)

    tools_menu.addSeparator()
    sketcher_action = QAction("2D-&Sketcher", window)
    sketcher_action.triggered.connect(window._open_sketcher)
    tools_menu.addAction(sketcher_action)
    
    tools_menu.addSeparator()
    
    # --- Lipophilicity Menu ---
    lipo_menu = tools_menu.addMenu("&Lipophilicity")
    
    logp_action = QAction("LogP", window)
    logp_action.triggered.connect(window._calculate_logp_action)
    lipo_menu.addAction(logp_action)
    
    lipo_contrib_action = QAction("Lipophilic Contribution", window)
    lipo_contrib_action.triggered.connect(window._show_lipophilic_contributions)
    lipo_menu.addAction(lipo_contrib_action)
    
    lipo_color_action = QAction("Color 3D-view (Gradient)", window)
    lipo_color_action.triggered.connect(window._color_by_lipophilicity)
    lipo_menu.addAction(lipo_color_action)
    
    window._lipo_menu_actions = [logp_action, lipo_contrib_action, lipo_color_action]


# ── Plugins ───────────────────────────────────────────────────────

def _build_plugins_menu(window, menu_bar):
    plugins_menu = menu_bar.addMenu("&Plugins")

    installed_plugins_action = QAction("&Installed Plugins...", window)
    installed_plugins_action.setShortcut(QKeySequence("Ctrl+Shift+I"))
    installed_plugins_action.triggered.connect(window._toggle_plugin_dock)
    plugins_menu.addAction(installed_plugins_action)

    plugins_menu.addSeparator()

    refresh_plugins_action = QAction("&Refresh Plugins", window)
    refresh_plugins_action.triggered.connect(window._refresh_plugins)
    plugins_menu.addAction(refresh_plugins_action)

    load_all_plugins_action = QAction("&Load All Plugins", window)
    load_all_plugins_action.triggered.connect(window._load_all_plugins)
    plugins_menu.addAction(load_all_plugins_action)

    unload_all_plugins_action = QAction("&Unload All Plugins", window)
    unload_all_plugins_action.triggered.connect(window._unload_all_plugins)
    plugins_menu.addAction(unload_all_plugins_action)


# ── View ──────────────────────────────────────────────────────────

def _build_view_menu(window, menu_bar):
    from src.shared.ui.theme import COLORS  # noqa: F401 – only needed at call time
    from src.shared.qt_compat import QActionGroup

    view_menu = menu_bar.addMenu("&View")

    view_menu.addSeparator()

    toggle_com = QAction("Toggle &COM Sphere", window)
    toggle_com.setCheckable(True)
    toggle_com.setChecked(False)
    toggle_com.triggered.connect(window._toggle_com_sphere)
    window._toggle_com_action = toggle_com
    view_menu.addAction(toggle_com)

    toggle_centroid = QAction("Toggle &Centroid Sphere", window)
    toggle_centroid.setCheckable(True)
    toggle_centroid.setChecked(False)
    toggle_centroid.triggered.connect(window._toggle_centroid_sphere)
    window._toggle_centroid_action = toggle_centroid
    view_menu.addAction(toggle_centroid)

    view_menu.addSeparator()

    bg_color_3d_action = QAction("3D Background Color...", window)
    bg_color_3d_action.setShortcut(QKeySequence("Ctrl+Shift+B"))
    bg_color_3d_action.triggered.connect(window._change_3d_bg_color)
    view_menu.addAction(bg_color_3d_action)

    bg_color_2d_action = QAction("2D Background Color...", window)
    bg_color_2d_action.setShortcut(QKeySequence("Ctrl+Alt+B"))
    bg_color_2d_action.triggered.connect(window._change_2d_bg_color)
    view_menu.addAction(bg_color_2d_action)

    residue_color_action = QAction("Residue Label Settings...", window)
    residue_color_action.setShortcut(QKeySequence("Ctrl+Shift+R"))
    residue_color_action.triggered.connect(window._change_residue_colors)
    view_menu.addAction(residue_color_action)

    view_menu.addSeparator()

    # ── Theme submenu (Light / Dark / System) ─────────────────
    from src.shared.ui.theme import ThemeMode, current_mode
    theme_menu = view_menu.addMenu("&Theme")

    theme_group = QActionGroup(window)
    theme_group.setExclusive(True)

    active_mode = current_mode()

    theme_system_action = QAction("&System (follow OS)", window)
    theme_system_action.setCheckable(True)
    theme_system_action.setChecked(active_mode == ThemeMode.SYSTEM)
    theme_system_action.triggered.connect(
        lambda: window._set_theme(ThemeMode.SYSTEM))
    theme_group.addAction(theme_system_action)
    theme_menu.addAction(theme_system_action)

    theme_menu.addSeparator()

    theme_light_action = QAction("&Light", window)
    theme_light_action.setCheckable(True)
    theme_light_action.setChecked(active_mode == ThemeMode.LIGHT)
    theme_light_action.triggered.connect(
        lambda: window._set_theme(ThemeMode.LIGHT))
    theme_group.addAction(theme_light_action)
    theme_menu.addAction(theme_light_action)

    theme_dark_action = QAction("&Dark", window)
    theme_dark_action.setCheckable(True)
    theme_dark_action.setChecked(active_mode == ThemeMode.DARK)
    theme_dark_action.triggered.connect(
        lambda: window._set_theme(ThemeMode.DARK))
    theme_group.addAction(theme_dark_action)
    theme_menu.addAction(theme_dark_action)

    # Expose the actions on the window so _set_theme can keep the
    # check marks in sync when the OS scheme changes under SYSTEM.
    window._theme_actions = {
        ThemeMode.SYSTEM: theme_system_action,
        ThemeMode.LIGHT: theme_light_action,
        ThemeMode.DARK: theme_dark_action,
    }


# ── Applications ──────────────────────────────────────────────────

def _build_applications_menu(window, menu_bar):
    applications_menu = menu_bar.addMenu("&Applications")

    descriptor_calc = QAction("Molecular &Descriptor Calculator", window)
    descriptor_calc.setShortcut(QKeySequence("Ctrl+D"))
    descriptor_calc.triggered.connect(window._open_descriptor_calculator)
    applications_menu.addAction(descriptor_calc)

    applications_menu.addSeparator()


    # Optimize geometry sub-menu
    opt_menu = applications_menu.addMenu("&Optimize Geometry")

    opt_mmff94 = QAction("MMFF94", window)
    opt_mmff94.triggered.connect(lambda checked: window._optimize_geometry(method="MMFF94"))
    opt_menu.addAction(opt_mmff94)

    opt_am1 = QAction("AM1", window)
    opt_am1.triggered.connect(lambda checked: window._optimize_geometry(method="AM1"))
    opt_menu.addAction(opt_am1)

    window.opt_menu_actions = [opt_mmff94, opt_am1]

    # Compute charges sub-menu
    chg_menu = applications_menu.addMenu("&Compute Charges")

    chg_gasteiger = QAction("Gasteiger", window)
    chg_gasteiger.triggered.connect(lambda checked: window._compute_charges(method="Gasteiger"))
    chg_menu.addAction(chg_gasteiger)

    chg_mmff94 = QAction("MMFF94", window)
    chg_mmff94.triggered.connect(lambda checked: window._compute_charges(method="MMFF94"))
    chg_menu.addAction(chg_mmff94)

    chg_am1 = QAction("AM1", window)
    chg_am1.triggered.connect(lambda checked: window._compute_charges(method="AM1"))
    chg_menu.addAction(chg_am1)

    chg_pm3 = QAction("PM3", window)
    chg_pm3.triggered.connect(lambda checked: window._compute_charges(method="PM3"))
    chg_menu.addAction(chg_pm3)

    window.chg_menu_actions = [chg_gasteiger, chg_mmff94, chg_am1, chg_pm3]

    applications_menu.addSeparator()

    window.aromaticity_action = QAction("Perceive Aromaticity", window)
    window.aromaticity_action.triggered.connect(window._perceive_aromaticity_action)
    applications_menu.addAction(window.aromaticity_action)

    window.docking_pose_action = QAction("Molecular Docking Pose in 3D", window)
    window.docking_pose_action.triggered.connect(window._open_docking_pose_view)
    applications_menu.addAction(window.docking_pose_action)

    applications_menu.addSeparator()

    find_substructure_action = QAction("Find &Substructure...", window)
    find_substructure_action.setShortcut(QKeySequence("Ctrl+F"))
    find_substructure_action.triggered.connect(window._show_substructure_dialog)
    applications_menu.addAction(find_substructure_action)

    copy_smiles_action = QAction("Copy S&MILES", window)
    copy_smiles_action.setShortcut(QKeySequence("Ctrl+Shift+C"))
    copy_smiles_action.triggered.connect(window._copy_smiles_to_clipboard)
    applications_menu.addAction(copy_smiles_action)

    applications_menu.addSeparator()

    smiles_action = QAction("&SMILES Generation", window)
    smiles_action.setShortcut(QKeySequence("Ctrl+Shift+S"))
    smiles_action.triggered.connect(window._generate_smiles)
    applications_menu.addAction(smiles_action)

    applications_menu.addSeparator()

    window.color_action = QAction("Atom &Colors...", window)
    window.color_action.triggered.connect(window._show_color_dialog)
    applications_menu.addAction(window.color_action)

    window.protein_color_action = QAction("&Protein Colors...", window)
    window.protein_color_action.triggered.connect(window._show_protein_color_dialog)
    applications_menu.addAction(window.protein_color_action)

    applications_menu.addSeparator()

    tools_submenu = applications_menu.addMenu("&Sphere Tools")

    window.com_sphere_action = QAction("COM Sphere", window)
    window.com_sphere_action.triggered.connect(window._add_com_sphere)
    tools_submenu.addAction(window.com_sphere_action)

    window.centroid_action_tool = QAction("Centroid", window)
    window.centroid_action_tool.triggered.connect(window._add_centroid_sphere)
    tools_submenu.addAction(window.centroid_action_tool)

    window.custom_sphere_action = QAction("Custom Sphere", window)
    window.custom_sphere_action.triggered.connect(window._add_custom_sphere)
    tools_submenu.addAction(window.custom_sphere_action)

    tools_submenu.addSeparator()

    window.clear_spheres_action = QAction("Clear Spheres", window)
    window.clear_spheres_action.triggered.connect(window._clear_all_spheres)
    tools_submenu.addAction(window.clear_spheres_action)

    window.tools_submenu_actions = [
        window.com_sphere_action, window.centroid_action_tool,
        window.custom_sphere_action, window.clear_spheres_action,
    ]

    # Disable molecule-dependent actions initially
    for action in (window.opt_menu_actions + window.chg_menu_actions +
                   [window.color_action, window.protein_color_action,
                    window.aromaticity_action, window.docking_pose_action] + 
                    window.tools_submenu_actions + window._lipo_menu_actions):
        action.setEnabled(False)


# ── Support ───────────────────────────────────────────────────────

def _build_support_menu(window, menu_bar):
    support_menu = menu_bar.addMenu("&Support")

    coffee_action = QAction("☕ Buy Me a Coffee", window)
    coffee_action.setStatusTip("Support PyChem development")
    coffee_action.triggered.connect(lambda: _open_url("https://buymeacoffee.com/vijaymasand"))
    support_menu.addAction(coffee_action)


def _open_url(url):
    try:
        from PySide6.QtGui import QDesktopServices
        from PySide6.QtCore import QUrl
        QDesktopServices.openUrl(QUrl(url))
    except Exception:
        import webbrowser
        webbrowser.open(url)


# ── Help ──────────────────────────────────────────────────────────

def _build_help_menu(window, menu_bar):
    help_menu = menu_bar.addMenu("&Help")

    about_action = QAction("&About", window)
    about_action.triggered.connect(window._show_about)
    help_menu.addAction(about_action)
