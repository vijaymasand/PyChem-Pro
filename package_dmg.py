"""
Wrap the compiled dist/PyChem.app into a signed-layout DMG.

Assumes `python build.py` has already produced dist/PyChem.app on macOS.
Output: dist/PyChem-<version>-macos.dmg
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
APP = ROOT / "dist" / "PyChem.app"
VERSION = "1.0.0"
DMG_PATH = ROOT / "dist" / f"PyChem-{VERSION}-macos.dmg"
SETTINGS = ROOT / "dmg_settings.py"


def main() -> None:
    if not APP.exists():
        sys.exit(f"Missing {APP}. Run `python build.py` first.")
    if DMG_PATH.exists():
        DMG_PATH.unlink()

    cmd = [
        sys.executable, "-m", "dmgbuild",
        "-s", str(SETTINGS),
        "-D", f"app={APP}",
        "PyChem",
        str(DMG_PATH),
    ]
    print("Building DMG:", " ".join(cmd))
    r = subprocess.run(cmd, cwd=ROOT)
    if r.returncode != 0:
        sys.exit(r.returncode)
    print(f"\nDMG ready -> {DMG_PATH}")
    print(f"Size: {DMG_PATH.stat().st_size / (1024*1024):.1f} MB")


if __name__ == "__main__":
    main()
