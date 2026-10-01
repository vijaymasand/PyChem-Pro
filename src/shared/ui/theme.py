"""
PyChem Theme — Refined research-workbench aesthetic.

Inspired by the Google Cloud Console: soft low-contrast surfaces,
a single muted blue accent, small rounded corners, system-native
typography. Ships two palettes (light / dark) sharing the same
semantic token names so the stylesheet template is reused.

Architecture
------------
* ``COLORS``            — the currently active palette, mutated in
                          place when the theme is switched so every
                          ``from src.shared.ui.theme import COLORS``
                          reference picks up the new values.
* ``set_theme(mode)``   — switch to ``ThemeMode.LIGHT``,
                          ``DARK`` or ``SYSTEM``; emits
                          ``theme_signals().theme_changed``.
* ``current_mode()``    — the mode the user chose (may be SYSTEM).
* ``effective_mode()``  — the resolved mode actually in use
                          (either LIGHT or DARK).
* ``get_stylesheet()``  — the Qt stylesheet for the active palette.
* ``theme_signals()``   — ``QObject`` exposing ``theme_changed``
                          signal so viewers can invalidate caches.

The module does NOT touch any canvas / renderer colour keys.  The
``atom_*``, ``ss_*`` and ``viewer_bg`` tokens are identical in both
palettes — the molecular renderer, protein cartoon service and
atom / bond CPK maps continue to read the same values they always
have.
"""

from __future__ import annotations

from enum import Enum
from typing import Optional


# ─── Canvas / renderer colours ─────────────────────────────────────
#
# These keys are read by the 3D renderer, the protein cartoon
# service and the atom/bond colour tables. They are identical in
# both light and dark themes — we never touch canvas colours.

_CANVAS_COLORS = {
    'viewer_bg':       '#ffffff',

    # CPK atom colours
    'atom_default':    '#808080',
    'atom_h':          '#d0d0d0',
    'atom_c':          '#55ff7f',
    'atom_n':          '#3050f8',
    'atom_o':          '#ff0d0d',
    'atom_f':          '#90e050',
    'atom_p':          '#ff8000',
    'atom_s':          '#ffff30',
    'atom_cl':         '#1ff01f',
    'atom_br':         '#a62929',
    'atom_i':          '#940094',
    'atom_selected':   '#ff00ff',
    'atom_highlight':  '#ffff00',
    'atom_positive':   '#0000ff',
    'atom_negative':   '#ff0000',

    # Protein secondary-structure colours
    'ss_helix':        '#dc3232',
    'ss_sheet':        '#3296dc',
    'ss_coil':         '#b4b4b4',
    'ss_turn':         '#00d4aa',
    'ss_default':      '#808080',
}


# ─── Light palette (JetBrains-style, muted steel) ─────────────────
#
# Editor surfaces from VS Code Light+ (clean whites) paired with a
# more restrained steel-blue accent in the spirit of IntelliJ /
# Darcula's light variant. The punchy VS Code `#007ACC` status-bar
# blue reads "Microsoft product" — this palette leans professional
# and neutral instead.

_LIGHT_CHROME = {
    # Surfaces
    'bg_primary':      '#FFFFFF',
    'bg_secondary':    '#F4F4F4',
    'bg_tertiary':     '#ECECEC',
    'bg_widget':       '#FFFFFF',
    'bg_hover':        '#E6E6E6',
    'bg_active':       '#D6E4F5',

    # Lines
    'border':          '#D5D5D5',
    'border_focus':    '#2B5B8F',

    # Type
    'text_primary':    '#222222',
    'text_secondary':  '#555555',
    'text_muted':      '#8A8A8A',

    # Accent — muted professional steel (not a shout)
    'accent':          '#2B5B8F',   # deep steel blue on white
    'accent_hover':    '#1F4373',
    'accent2':         '#2E7F7E',   # muted teal
    'accent3':         '#79602B',   # muted olive

    # Semantic
    'success':         '#2E7F7E',
    'warning':         '#A07A16',
    'error':           '#B9423A',

    # Scrollbars
    'scrollbar_bg':    '#F4F4F4',
    'scrollbar_handle':'#BFBFBF',
}


# ─── Dark palette (JetBrains Darcula-inspired, muted steel) ──────
#
# Warm neutral greys modelled on the JetBrains Darcula family
# rather than pure VS Code Dark+. The accent is a muted steel blue
# `#5585B5` — clearly visible on dark surfaces but never neon.

_DARK_CHROME = {
    # Surfaces
    'bg_primary':      '#1E1F22',
    'bg_secondary':    '#2B2D30',
    'bg_tertiary':     '#313438',
    'bg_widget':       '#2B2D30',
    'bg_hover':        '#3C3F41',
    'bg_active':       '#3D5A7A',

    # Lines
    'border':          '#3C3F41',
    'border_focus':    '#5585B5',

    # Type
    'text_primary':    '#DFE1E5',
    'text_secondary':  '#A1A4AB',
    'text_muted':      '#6F7480',

    # Accent — muted steel
    'accent':          '#5585B5',   # never eye-melting, always readable
    'accent_hover':    '#6D9CCC',
    'accent2':         '#6AAB9C',
    'accent3':         '#BFA46F',

    # Semantic
    'success':         '#6AAB9C',
    'warning':         '#BFA46F',
    'error':           '#CF6C6C',

    # Scrollbars
    'scrollbar_bg':    '#1E1F22',
    'scrollbar_handle':'#494C52',
}


# ─── Active palette (mutated in place on switch) ──────────────────
#
# ``COLORS`` starts out populated with the light palette so that any
# module importing it at import time gets valid values. The actual
# theme chosen by the user is applied at application startup via
# ``set_theme()``.

COLORS: dict = {**_LIGHT_CHROME, **_CANVAS_COLORS}


# ─── Typography ───────────────────────────────────────────────────
#
# System-native fonts on every platform (this is exactly what GCP
# does on desktop). Distinctive without bundling.

FONT_UI = (
    "'SF Pro Text', 'Segoe UI Variable', 'Segoe UI', "
    "'Helvetica Neue', 'Noto Sans', sans-serif"
)

FONT_MONO = (
    "'SF Mono', 'JetBrains Mono', 'Cascadia Code', "
    "'Menlo', 'Consolas', monospace"
)


# ─── Theme mode management ────────────────────────────────────────

class ThemeMode(Enum):
    """User-selected theme preference."""
    SYSTEM = "system"
    LIGHT = "light"
    DARK = "dark"


_current_mode: ThemeMode = ThemeMode.SYSTEM
_effective_mode: ThemeMode = ThemeMode.LIGHT


def current_mode() -> ThemeMode:
    """Return the user's chosen theme mode (may be SYSTEM)."""
    return _current_mode


def effective_mode() -> ThemeMode:
    """Return the resolved mode actually in use (LIGHT or DARK)."""
    return _effective_mode


def detect_system_mode() -> ThemeMode:
    """
    Query the OS for its current colour scheme preference.

    Uses ``QStyleHints.colorScheme()`` on Qt 6.5+. Falls back to
    LIGHT if the application instance or the API is unavailable
    (e.g. headless environments or older Qt builds).
    """
    try:
        from PySide6.QtCore import Qt
        from PySide6.QtGui import QGuiApplication
        app = QGuiApplication.instance()
        if app is None:
            return ThemeMode.LIGHT
        hints = app.styleHints()
        scheme = getattr(hints, 'colorScheme', None)
        if scheme is None:
            return ThemeMode.LIGHT
        value = scheme()
        if value == Qt.ColorScheme.Dark:
            return ThemeMode.DARK
        return ThemeMode.LIGHT
    except Exception:
        return ThemeMode.LIGHT


# ─── Theme signal plumbing ────────────────────────────────────────

_signals_instance = None


def theme_signals():
    """
    Return a shared ``QObject`` whose ``theme_changed`` signal fires
    after every theme swap. Instantiated lazily so the module can be
    imported before the ``QApplication`` exists.
    """
    global _signals_instance
    if _signals_instance is None:
        from PySide6.QtCore import QObject, Signal

        class _ThemeSignals(QObject):
            theme_changed = Signal()

        _signals_instance = _ThemeSignals()
    return _signals_instance


# ─── Theme application ───────────────────────────────────────────

def set_theme(mode: ThemeMode) -> None:
    """
    Switch to ``mode``, mutate ``COLORS`` in place, and re-apply the
    Qt stylesheet to the running ``QApplication``. Emits
    ``theme_signals().theme_changed`` so listeners can refresh.

    Safe to call before a ``QApplication`` exists — in that case the
    stylesheet re-application is skipped and will happen on the next
    ``set_theme()`` call made from ``main.py``.
    """
    global _current_mode, _effective_mode
    _current_mode = mode

    if mode == ThemeMode.SYSTEM:
        _effective_mode = detect_system_mode()
    else:
        _effective_mode = mode

    chrome = _DARK_CHROME if _effective_mode == ThemeMode.DARK else _LIGHT_CHROME
    new_palette = {**chrome, **_CANVAS_COLORS}

    # Mutate in place so every existing reference to COLORS sees the
    # new values. Do NOT replace the dict object.
    COLORS.clear()
    COLORS.update(new_palette)

    # Re-apply the stylesheet if a QApplication is running.
    # Qt defers stylesheet repolish to the next paint event, and some
    # child widgets (notably QScrollArea viewports) cache their own
    # background from the palette rather than the stylesheet — which
    # means a plain ``app.setStyleSheet()`` does NOT update every
    # widget in the tree on a live theme swap. We therefore:
    #   1. Clear + reapply the stylesheet.
    #   2. Walk EVERY widget in the application (not just top-level)
    #      and force an unpolish + polish cycle so the new rules are
    #      evaluated against each one individually.
    #   3. Call update() so the next paint redraws with the new
    #      background.
    try:
        from PySide6.QtWidgets import QApplication
        app = QApplication.instance()
        if app is not None:
            new_ss = get_stylesheet()
            app.setStyleSheet("")
            app.setStyleSheet(new_ss)
            for w in app.allWidgets():
                try:
                    st = w.style()
                    st.unpolish(w)
                    st.polish(w)
                    w.update()
                except Exception:
                    pass
    except Exception:
        pass

    # Notify listeners (viewer caches, widgets with inline stylesheets)
    try:
        theme_signals().theme_changed.emit()
    except Exception:
        pass


def apply_initial_theme() -> ThemeMode:
    """
    Read the saved theme preference from ``QSettings`` and apply it.
    Intended to be called once after the ``QApplication`` is created
    but before the main window is shown. Returns the effective mode.
    """
    try:
        from PySide6.QtCore import QSettings
        settings = QSettings("PyChem", "Viewer")
        saved = settings.value("theme/mode", "system")
        mode_map = {
            "system": ThemeMode.SYSTEM,
            "light":  ThemeMode.LIGHT,
            "dark":   ThemeMode.DARK,
        }
        mode = mode_map.get(str(saved).lower(), ThemeMode.SYSTEM)
    except Exception:
        mode = ThemeMode.SYSTEM

    set_theme(mode)

    # Auto-swap when the OS colour scheme changes while SYSTEM is
    # selected. Wire the connection once here.
    _connect_system_watcher()

    return _effective_mode


def save_theme_preference(mode: ThemeMode) -> None:
    """Persist the user's choice to ``QSettings``."""
    try:
        from PySide6.QtCore import QSettings
        settings = QSettings("PyChem", "Viewer")
        settings.setValue("theme/mode", mode.value)
    except Exception:
        pass


_system_watcher_connected = False


def _connect_system_watcher() -> None:
    """
    When the user has SYSTEM selected, listen for OS-level colour
    scheme changes and automatically re-apply the theme so dark/light
    follows the OS without a restart.
    """
    global _system_watcher_connected
    if _system_watcher_connected:
        return
    try:
        from PySide6.QtGui import QGuiApplication
        app = QGuiApplication.instance()
        if app is None:
            return
        hints = app.styleHints()
        if hints is None:
            return
        sig = getattr(hints, 'colorSchemeChanged', None)
        if sig is None:
            return

        def _on_os_scheme_changed(_new_scheme):
            if _current_mode == ThemeMode.SYSTEM:
                set_theme(ThemeMode.SYSTEM)

        sig.connect(_on_os_scheme_changed)
        _system_watcher_connected = True
    except Exception:
        pass


# ─── Stylesheet template ──────────────────────────────────────────

def get_stylesheet() -> str:
    """Generate the complete Qt stylesheet for the current palette."""
    c = COLORS
    return f"""
    /* ── Global ────────────────────────────────────────────── */
    QMainWindow {{
        background-color: {c['bg_primary']};
        color: {c['text_primary']};
    }}
    QWidget {{
        background-color: transparent;
        color: {c['text_primary']};
        font-family: {FONT_UI};
        font-size: 12px;
    }}
    /* The left panel's background is painted explicitly by the
       InputPanel widget's own inline stylesheet from _apply_theme()
       — Qt's QScrollArea palette caching on macOS prevents the
       cascade from updating reliably on live theme swaps. We still
       keep a minimal rule here for the hairline right-hand border so
       it updates via the global re-polish. */
    QWidget#leftPanel QScrollArea {{
        border: none;
    }}

    /* ── Menu bar ──────────────────────────────────────────── */
    QMenuBar {{
        background-color: {c['bg_secondary']};
        color: {c['text_primary']};
        border-bottom: 1px solid {c['border']};
        padding: 2px 4px;
        font-size: 13px;
    }}
    QMenuBar::item {{
        background-color: transparent;
        padding: 6px 12px;
        border-radius: 4px;
        margin: 1px;
    }}
    QMenuBar::item:selected {{
        background-color: {c['bg_hover']};
        color: {c['accent']};
    }}
    QMenuBar::item:pressed {{
        background-color: {c['bg_active']};
    }}

    QMenu {{
        background-color: {c['bg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 6px;
        padding: 4px;
        font-size: 13px;
    }}
    QMenu::item {{
        padding: 7px 26px 7px 16px;
        border-radius: 4px;
        color: {c['text_primary']};
    }}
    QMenu::item:selected {{
        background-color: {c['bg_hover']};
        color: {c['accent']};
    }}
    QMenu::item:disabled {{
        color: {c['text_muted']};
    }}
    QMenu::separator {{
        height: 1px;
        background: {c['border']};
        margin: 4px 8px;
    }}
    QMenu::indicator {{
        width: 14px;
        height: 14px;
        left: 6px;
    }}

    /* ── Buttons ───────────────────────────────────────────── */
    /* Default button: quiet surface with subtle hairline border. */
    QPushButton {{
        background-color: {c['bg_tertiary']};
        color: {c['text_primary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 6px 14px;
        font-size: 12px;
        font-weight: 500;
        min-height: 18px;
    }}
    QPushButton:hover {{
        background-color: {c['bg_hover']};
        border-color: {c['accent']};
        color: {c['accent']};
    }}
    QPushButton:pressed {{
        background-color: {c['bg_active']};
    }}
    QPushButton:disabled {{
        background-color: {c['bg_secondary']};
        color: {c['text_muted']};
        border-color: {c['border']};
    }}
    QPushButton:focus {{
        outline: none;
    }}

    /* Primary CTA — the single blue hero button */
    QPushButton#btnSuccess {{
        background-color: {c['accent']};
        color: #FFFFFF;
        border: 1px solid {c['accent']};
        font-weight: 600;
    }}
    QPushButton#btnSuccess:hover {{
        background-color: {c['accent_hover']};
        border-color: {c['accent_hover']};
        color: #FFFFFF;
    }}
    QPushButton#btnSuccess:pressed {{
        background-color: {c['accent_hover']};
    }}
    QPushButton#btnSuccess:disabled {{
        background-color: {c['bg_tertiary']};
        color: {c['text_muted']};
        border-color: {c['border']};
    }}

    QPushButton#btnSecondary {{
        background-color: {c['bg_tertiary']};
        color: {c['text_primary']};
        border: 1px solid {c['border']};
    }}
    QPushButton#btnSecondary:hover {{
        background-color: {c['bg_hover']};
        color: {c['accent']};
        border-color: {c['accent']};
    }}

    QPushButton#btnTertiary {{
        background-color: transparent;
        color: {c['text_secondary']};
        border: 1px solid {c['border']};
        padding: 2px 6px;
        font-size: 10px;
        font-weight: 400;
        min-height: 14px;
    }}
    QPushButton#btnTertiary:hover {{
        background-color: {c['bg_hover']};
        color: {c['accent']};
        border-color: {c['accent']};
    }}

    QPushButton#btnDanger {{
        background-color: {c['bg_tertiary']};
        color: {c['error']};
        border: 1px solid {c['error']};
    }}
    QPushButton#btnDanger:hover {{
        background-color: {c['error']};
        color: #FFFFFF;
    }}

    /* ── Text inputs ───────────────────────────────────────── */
    QLineEdit, QTextEdit, QPlainTextEdit {{
        background-color: {c['bg_widget']};
        color: {c['text_primary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 6px 10px;
        font-family: {FONT_MONO};
        font-size: 12px;
        selection-background-color: {c['accent']};
        selection-color: #FFFFFF;
    }}
    QLineEdit:hover, QTextEdit:hover, QPlainTextEdit:hover {{
        border-color: {c['text_muted']};
    }}
    QLineEdit:focus, QTextEdit:focus, QPlainTextEdit:focus {{
        border-color: {c['accent']};
    }}
    QLineEdit:disabled, QTextEdit:disabled, QPlainTextEdit:disabled {{
        color: {c['text_muted']};
        background-color: {c['bg_secondary']};
    }}

    /* ── Labels ────────────────────────────────────────────── */
    QLabel {{
        color: {c['text_primary']};
        padding: 0;
        background: transparent;
    }}
    QLabel#labelTitle {{
        color: {c['text_primary']};
        font-family: {FONT_UI};
        font-size: 17px;
        font-weight: 700;
        padding: 0;
    }}
    QLabel#labelSubtitle {{
        color: {c['text_secondary']};
        font-family: {FONT_UI};
        font-size: 11px;
        font-weight: 400;
        padding: 0;
    }}
    QLabel#labelSection {{
        color: {c['text_secondary']};
        font-family: {FONT_UI};
        font-size: 10px;
        font-weight: 600;
        letter-spacing: 0.4px;
        padding: 1px 0;
    }}
    QLabel#labelMuted {{
        color: {c['text_muted']};
        font-family: {FONT_UI};
        font-size: 10px;
        font-weight: 400;
    }}
    QLabel#labelHint {{
        color: {c['text_muted']};
        font-family: {FONT_MONO};
        font-size: 10px;
        font-weight: 400;
    }}
    QLabel#labelData {{
        color: {c['text_primary']};
        font-family: {FONT_MONO};
        font-size: 11px;
    }}

    /* ── Group box ─────────────────────────────────────────── */
    /*
       Tight layout. The title sits in the 18 px margin strip ABOVE
       the frame so it never overlaps the border line. QGroupBox
       padding is kept small -- the inner layout should set its own
       content margins to zero because the QGroupBox padding already
       provides the inset.  This is the single biggest cause of the
       "floating in empty space" look that plagued earlier iterations.
    */
    QGroupBox {{
        background-color: {c['bg_tertiary']};
        border: 1px solid {c['border']};
        border-radius: 5px;
        margin-top: 18px;
        padding: 10px 10px 10px 10px;
        font-family: {FONT_UI};
        font-size: 10px;
        font-weight: 600;
        color: {c['text_secondary']};
    }}
    QGroupBox::title {{
        subcontrol-origin: margin;
        subcontrol-position: top left;
        left: 2px;
        top: 0;
        padding: 0;
        background: transparent;
        color: {c['text_secondary']};
    }}

    /* ── Progress bar ──────────────────────────────────────── */
    QProgressBar {{
        background-color: {c['bg_widget']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        text-align: center;
        color: {c['text_secondary']};
        font-size: 10px;
        min-height: 8px;
        max-height: 8px;
    }}
    QProgressBar::chunk {{
        background-color: {c['accent']};
        border-radius: 3px;
    }}

    /* ── Status bar ────────────────────────────────────────── */
    QStatusBar {{
        background-color: {c['bg_secondary']};
        color: {c['text_secondary']};
        border-top: 1px solid {c['border']};
        padding: 3px 10px;
        font-size: 12px;
    }}
    QStatusBar::item {{ border: none; }}

    /* ── Scrollbars ────────────────────────────────────────── */
    QScrollBar:vertical {{
        background: transparent;
        width: 10px;
        margin: 0;
        border: none;
    }}
    QScrollBar::handle:vertical {{
        background: {c['scrollbar_handle']};
        border-radius: 3px;
        min-height: 24px;
        margin: 2px 3px;
    }}
    QScrollBar::handle:vertical:hover {{
        background: {c['accent']};
    }}
    QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
        height: 0; background: none;
    }}
    QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {{
        background: transparent;
    }}
    QScrollBar:horizontal {{
        background: transparent;
        height: 10px;
        margin: 0;
        border: none;
    }}
    QScrollBar::handle:horizontal {{
        background: {c['scrollbar_handle']};
        border-radius: 3px;
        min-width: 24px;
        margin: 3px 2px;
    }}
    QScrollBar::handle:horizontal:hover {{
        background: {c['accent']};
    }}
    QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
        width: 0; background: none;
    }}
    QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal {{
        background: transparent;
    }}

    /* ── Splitters ─────────────────────────────────────────── */
    QSplitter::handle {{
        background-color: {c['border']};
    }}
    QSplitter::handle:horizontal {{ width: 1px; }}
    QSplitter::handle:vertical {{ height: 1px; }}
    QSplitter::handle:hover {{
        background-color: {c['accent']};
    }}

    /* ── Tooltip ───────────────────────────────────────────── */
    QToolTip {{
        background-color: {c['bg_tertiary']};
        color: {c['text_primary']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 6px 10px;
        font-size: 12px;
    }}

    /* ── Tab widget (VS Code editor-tab style) ─────────────── */
    /*
       All three sub-regions of QTabWidget need an explicit
       background so Qt never falls through to the native dark
       strip behaviour in light mode.
    */
    QTabWidget {{
        background-color: {c['bg_primary']};
    }}
    QTabWidget::pane {{
        border: 1px solid {c['border']};
        background-color: {c['bg_primary']};
        top: 0;
    }}
    QTabWidget::tab-bar {{
        left: 0;
        alignment: left;
    }}
    QTabBar {{
        qproperty-drawBase: 0;
        background-color: {c['bg_tertiary']};
        border-bottom: 1px solid {c['border']};
    }}
    QTabBar::tab {{
        background-color: {c['bg_tertiary']};
        color: {c['text_secondary']};
        padding: 9px 22px;
        min-width: 80px;
        border: none;
        border-right: 1px solid {c['border']};
        border-bottom: 1px solid {c['border']};
        font-size: 12px;
        font-weight: 500;
    }}
    QTabBar::tab:selected {{
        background-color: {c['bg_primary']};
        color: {c['text_primary']};
        border-bottom: 2px solid {c['accent']};
    }}
    QTabBar::tab:hover:!selected {{
        color: {c['text_primary']};
        background-color: {c['bg_hover']};
    }}
    QTabBar::scroller {{
        width: 0;
    }}

    /* ── File dialog ───────────────────────────────────────── */
    QFileDialog {{
        background-color: {c['bg_primary']};
    }}
    QFileDialog QListView, QFileDialog QTreeView {{
        background-color: {c['bg_widget']};
        color: {c['text_primary']};
        border: 1px solid {c['border']};
    }}

    /* ── Combo box ─────────────────────────────────────────── */
    QComboBox {{
        background-color: {c['bg_widget']};
        border: 1px solid {c['border']};
        border-radius: 4px;
        padding: 4px 8px;
        color: {c['text_primary']};
        font-size: 12px;
        min-height: 18px;
    }}
    QComboBox:hover {{
        border-color: {c['text_muted']};
    }}
    QComboBox:focus, QComboBox:on {{
        border-color: {c['accent']};
    }}
    QComboBox::drop-down {{
        border: none;
        width: 22px;
        background: transparent;
    }}
    QComboBox QAbstractItemView {{
        background-color: {c['bg_secondary']};
        border: 1px solid {c['border']};
        border-radius: 6px;
        color: {c['text_primary']};
        selection-background-color: {c['bg_hover']};
        selection-color: {c['accent']};
        outline: none;
        padding: 2px;
    }}

    /* ── Check box ─────────────────────────────────────────── */
    QCheckBox {{
        color: {c['text_primary']};
        spacing: 7px;
        font-size: 12px;
        padding: 1px 0;
    }}
    QCheckBox::indicator {{
        width: 13px;
        height: 13px;
        border: 1px solid {c['border']};
        border-radius: 2px;
        background-color: {c['bg_widget']};
    }}
    QCheckBox::indicator:hover {{
        border-color: {c['accent']};
    }}
    QCheckBox::indicator:checked {{
        background-color: {c['accent']};
        border-color: {c['accent']};
    }}
    QCheckBox::indicator:disabled {{
        background-color: {c['bg_secondary']};
        border-color: {c['border']};
    }}

    /* ── Radio button ──────────────────────────────────────── */
    QRadioButton {{
        color: {c['text_primary']};
        spacing: 8px;
        font-size: 13px;
    }}
    QRadioButton::indicator {{
        width: 14px;
        height: 14px;
        border: 1px solid {c['border']};
        border-radius: 7px;
        background-color: {c['bg_widget']};
    }}
    QRadioButton::indicator:checked {{
        background-color: {c['accent']};
        border-color: {c['accent']};
    }}

    /* ── Spin box ──────────────────────────────────────────── */
    QSpinBox, QDoubleSpinBox {{
        background-color: {c['bg_widget']};
        border: 1px solid {c['border']};
        border-radius: 6px;
        padding: 5px 8px;
        color: {c['text_primary']};
        font-family: {FONT_MONO};
        font-size: 12px;
        selection-background-color: {c['accent']};
        selection-color: #FFFFFF;
    }}
    QSpinBox:focus, QDoubleSpinBox:focus {{
        border-color: {c['accent']};
    }}
    QSpinBox::up-button, QDoubleSpinBox::up-button {{
        subcontrol-origin: border;
        subcontrol-position: top right;
        width: 16px;
        background: {c['bg_tertiary']};
        border-left: 1px solid {c['border']};
        border-top-right-radius: 5px;
    }}
    QSpinBox::down-button, QDoubleSpinBox::down-button {{
        subcontrol-origin: border;
        subcontrol-position: bottom right;
        width: 16px;
        background: {c['bg_tertiary']};
        border-left: 1px solid {c['border']};
        border-top: 1px solid {c['border']};
        border-bottom-right-radius: 5px;
    }}
    QSpinBox::up-button:hover, QSpinBox::down-button:hover,
    QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover {{
        background: {c['bg_hover']};
    }}

    /* ── Slider ────────────────────────────────────────────── */
    QSlider::groove:horizontal {{
        background: {c['border']};
        height: 3px;
        border-radius: 2px;
    }}
    QSlider::sub-page:horizontal {{
        background: {c['accent']};
        height: 3px;
        border-radius: 2px;
    }}
    QSlider::add-page:horizontal {{
        background: {c['border']};
        height: 3px;
        border-radius: 2px;
    }}
    QSlider::handle:horizontal {{
        background: {c['bg_secondary']};
        width: 14px;
        height: 14px;
        margin: -6px 0;
        border: 2px solid {c['accent']};
        border-radius: 8px;
    }}
    QSlider::handle:horizontal:hover {{
        background: {c['accent']};
    }}
    QSlider::handle:horizontal:pressed {{
        background: {c['accent_hover']};
        border-color: {c['accent_hover']};
    }}

    /* ── List widget ───────────────────────────────────────── */
    QListWidget {{
        background-color: {c['bg_widget']};
        color: {c['text_primary']};
        border: 1px solid {c['border']};
        border-radius: 6px;
        font-size: 12px;
        padding: 2px;
        outline: none;
    }}
    QListWidget::item {{
        padding: 7px 10px;
        border: none;
        border-radius: 4px;
    }}
    QListWidget::item:hover {{
        background-color: {c['bg_hover']};
    }}
    QListWidget::item:selected {{
        background-color: {c['bg_hover']};
        color: {c['accent']};
    }}

    /* ── Dock widget ───────────────────────────────────────── */
    QDockWidget {{
        background-color: {c['bg_secondary']};
        color: {c['text_primary']};
        font-size: 12px;
        font-weight: 500;
    }}
    QDockWidget::title {{
        background: {c['bg_tertiary']};
        color: {c['text_primary']};
        padding: 8px 12px;
        border-bottom: 1px solid {c['border']};
        text-align: left;
    }}

    /* ── Header view (tables) ──────────────────────────────── */
    QHeaderView::section {{
        background-color: {c['bg_tertiary']};
        color: {c['text_secondary']};
        padding: 8px 10px;
        border: none;
        border-right: 1px solid {c['border']};
        border-bottom: 1px solid {c['border']};
        font-size: 11px;
        font-weight: 600;
    }}
    QTableWidget, QTableView {{
        background-color: {c['bg_widget']};
        color: {c['text_primary']};
        gridline-color: {c['border']};
        border: 1px solid {c['border']};
        selection-background-color: {c['bg_hover']};
        selection-color: {c['accent']};
        font-family: {FONT_MONO};
        font-size: 12px;
        alternate-background-color: {c['bg_secondary']};
    }}

    /* ── Dialog ────────────────────────────────────────────── */
    QDialog {{
        background-color: {c['bg_primary']};
        color: {c['text_primary']};
    }}
    QMessageBox {{
        background-color: {c['bg_primary']};
    }}
    """
