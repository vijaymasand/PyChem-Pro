"""
SMILES to 3D Molecular Structure Converter
═══════════════════════════════════════════

Main entry point. Initializes the application, checks license,
and launches the GUI.

Performance optimizations enabled:
- Parallel file loading for large molecules
- Lazy imports for heavy modules
- Startup profiling
"""

import sys
import os
import time

# Startup timing
_startup_start = time.perf_counter()

def main():
    """Application entry point with performance optimizations."""
    # Ensure the src package is importable
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    # Import and initialize startup optimizer early
    try:
        from src.core.performance import StartupOptimizer, get_profiler
        optimizer = StartupOptimizer()
        optimizer.optimize_imports()  # Disable bytecode writing for speed
        _use_optimizer = True
    except ImportError:
        _use_optimizer = False
        print("[PERF] Performance optimizations not available")

    # Profile Qt import (can be slow)
    qt_import_start = time.perf_counter()
    from src.shared.qt_compat import QApplication, QFont, QIcon, Qt
    qt_import_time = time.perf_counter() - qt_import_start
    # print(f"[PERF] Qt imports: {qt_import_time:.3f}s")  # Commented out for reduced verbosity

    # High DPI support
    QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)

    app = QApplication(sys.argv)
    app.setApplicationName("SMILES to 3D Converter")
    app.setApplicationVersion("1.0.0")
    app.setOrganizationName("SMILES3D")

    icon_path = os.path.join(script_dir, "assets", "icon.svg")
    if os.path.exists(icon_path):
        app.setWindowIcon(QIcon(icon_path))

    # Set default font
    font = QFont("Segoe UI", 10)
    app.setFont(font)

    # License check
    from src.core.security.license import LicenseManager
    license_mgr = LicenseManager()
    is_valid, message, features = license_mgr.validate_license()

    if not is_valid:
        # Auto-generate a development license (for first run)
        # In production, this would show a license activation dialog
        license_mgr.generate_and_save(days_valid=3650)
        is_valid, message, features = license_mgr.validate_license()

    # Launch main window with timing
    main_window_start = time.perf_counter()
    from src.app.main_window import MainWindow
    window = MainWindow()
    window.show()
    main_window_time = time.perf_counter() - main_window_start
    # print(f"[PERF] MainWindow creation: {main_window_time:.3f}s")  # Commented out for reduced verbosity
    
    # Total startup time
    total_startup = time.perf_counter() - _startup_start
    # print(f"[PERF] Total startup time: {total_startup:.3f}s")  # Commented out for reduced verbosity
    
    # Print performance report if available
    if _use_optimizer:
        try:
            from src.core.performance import print_performance_report
            # Schedule report print on app exit or provide menu option
        except ImportError:
            pass

    sys.exit(app.exec())


if __name__ == "__main__":
    # Required for Windows multiprocessing support
    # This must be guarded for multiprocessing to work correctly
    main()
