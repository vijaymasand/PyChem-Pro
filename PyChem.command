#!/bin/bash
# PyChem launcher for macOS — double-click this file in Finder to run PyChem.
# Automatically installs Python and dependencies on first run.

# Move to the folder containing this script so relative paths work
cd "$(dirname "$0")"

# ── Helper: pretty section headers ──────────────────────────────────────────
info()  { echo ""; echo "▶  $*"; }
ok()    { echo "✓  $*"; }
fail()  { echo "✗  $*"; }

# ── Step 1: Find Python 3.10+ ────────────────────────────────────────────────
info "Checking for Python 3.10+..."

PYTHON=""
for candidate in python3.13 python3.12 python3.11 python3.10 python3; do
    # Also search Homebrew prefixes directly in case PATH is not set yet
    for prefix in "" "/opt/homebrew/bin/" "/usr/local/bin/"; do
        full="${prefix}${candidate}"
        if command -v "$full" &>/dev/null; then
            version=$("$full" -c "import sys; print(sys.version_info >= (3,10))" 2>/dev/null)
            if [ "$version" = "True" ]; then
                PYTHON="$full"
                break 2
            fi
        fi
    done
done

# ── Step 2: Auto-install Python via Homebrew if not found ────────────────────
if [ -z "$PYTHON" ]; then
    fail "Python 3.10+ not found."
    info "Setting up Homebrew to install Python automatically..."

    # Activate Homebrew if already installed but not on PATH
    if [ -f "/opt/homebrew/bin/brew" ]; then
        eval "$(/opt/homebrew/bin/brew shellenv)"
    elif [ -f "/usr/local/bin/brew" ]; then
        eval "$(/usr/local/bin/brew shellenv)"
    fi

    # Install Homebrew if still not available
    if ! command -v brew &>/dev/null; then
        info "Installing Homebrew (you will be asked for your Mac password)..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

        # Activate newly installed Homebrew
        if [ -f "/opt/homebrew/bin/brew" ]; then
            eval "$(/opt/homebrew/bin/brew shellenv)"
        elif [ -f "/usr/local/bin/brew" ]; then
            eval "$(/usr/local/bin/brew shellenv)"
        fi
    fi

    if command -v brew &>/dev/null; then
        info "Installing Python 3.12 via Homebrew..."
        brew install python@3.12
        # Point to the freshly installed binary directly
        PYTHON="$(brew --prefix python@3.12)/bin/python3"
    fi
fi

# ── Step 3: Final check — if still nothing, show alert ──────────────────────
if [ -z "$PYTHON" ] || ! command -v "$PYTHON" &>/dev/null; then
    fail "Could not install Python automatically."
    osascript -e '
        set btn to button returned of (display alert "Python Installation Failed" message "PyChem could not install Python automatically." & return & return & "Please install Python 3.10+ manually from python.org, then open PyChem again." as critical buttons {"Quit", "Open python.org"} default button "Open python.org")
        if btn is "Open python.org" then
            open location "https://www.python.org/downloads/"
        end if
    '
    exit 1
fi

ok "Using Python: $("$PYTHON" --version)"

# ── Step 4: Set Finder icon (once) ──────────────────────────────────────────
if [ ! -f ".launcher_icon_set" ]; then
    ICON="$(pwd)/assets/icon.icns"
    TARGET="$(pwd)/PyChem.command"
    osascript -l JavaScript -e "
        ObjC.import('AppKit');
        var img = \$.NSImage.alloc.initWithContentsOfFile('$ICON');
        \$.NSWorkspace.sharedWorkspace.setIconForFileOptions(img, '$TARGET', 0);
    " 2>/dev/null && touch .launcher_icon_set
fi

# ── Step 5: Create venv and install dependencies (first run only) ────────────
if [ ! -d "venv" ]; then
    info "First run — installing dependencies (this takes ~2 minutes)..."
    "$PYTHON" -m venv venv
    venv/bin/pip install --upgrade pip -q
    venv/bin/pip install -r requirements.txt
    ok "Dependencies installed."
fi

# ── Step 6: Launch PyChem ────────────────────────────────────────────────────
info "Launching PyChem..."
venv/bin/python main.py
