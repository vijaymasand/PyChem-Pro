"""
Build script — compile PyChem into a native standalone binary with Nuitka.

Produces:
  macOS   → dist/PyChem.app     (wrap with `python package_dmg.py` to get a .dmg)
  Windows → dist/PyChem.exe     (single-file)
  Linux   → dist/PyChem         (single-file)

Usage:
    python build.py                 # onefile/app-bundle
    python build.py --dir           # standalone directory (faster iteration)
"""

import os
import sys
import platform
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
APP_NAME = "PyChem"
VERSION = "1.0.0"
COMPANY = "PyChem"
DESCRIPTION = "PyChem - Molecular Visualization Software"


def build(onefile: bool = True) -> None:
    system = platform.system()
    icon_icns = ROOT / "assets" / "icon.icns"
    icon_ico = ROOT / "assets" / "icon.ico"

    cmd = [
        sys.executable, "-m", "nuitka",
        "--standalone",
        "--enable-plugin=pyside6",
        "--remove-output",
        "--assume-yes-for-downloads",
        "--output-dir=dist",
        f"--output-filename={APP_NAME}",
        f"--company-name={COMPANY}",
        f"--product-name={APP_NAME}",
        f"--file-version={VERSION}",
        f"--product-version={VERSION}",
        f"--file-description={DESCRIPTION}",
        "--include-package=src",
        "--include-package=pychem",
        "--include-data-files=plugins/*.py=plugins/",
        "--include-data-dir=assets=assets",
        "--include-qt-plugins=sensible,iconengines,imageformats",
        "--include-module=PySide6.QtSvg",
    ]

    if onefile:
        if system == "Darwin":
            # macOS app bundles are self-contained; onefile breaks .app structure
            pass
        else:
            cmd.append("--onefile")

    if system == "Windows":
        cmd.append("--windows-console-mode=disable")
        if icon_ico.exists():
            cmd.append(f"--windows-icon-from-ico={icon_ico}")
    elif system == "Darwin":
        cmd.append("--macos-create-app-bundle")
        cmd.append(f"--macos-app-name={APP_NAME}")
        cmd.append(f"--macos-app-version={VERSION}")
        if icon_icns.exists():
            cmd.append(f"--macos-app-icon={icon_icns}")
    elif system == "Linux":
        if icon_ico.exists():
            cmd.append(f"--linux-icon={icon_ico}")

    cmd.append("main.py")

    print("=" * 60)
    print(f"Building {APP_NAME} v{VERSION}")
    print(f"Platform: {system} ({platform.machine()})")
    print(f"Mode: {'onefile / app bundle' if onefile else 'standalone directory'}")
    print("=" * 60)
    print("Command:", " ".join(str(c) for c in cmd), "\n")

    result = subprocess.run(cmd, cwd=ROOT)
    if result.returncode != 0:
        print(f"\nBuild failed (exit {result.returncode})")
        sys.exit(result.returncode)
    print("\nBuild complete -> dist/")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--dir", action="store_true", help="standalone dir (no onefile)")
    args = p.parse_args()
    build(onefile=not args.dir)
