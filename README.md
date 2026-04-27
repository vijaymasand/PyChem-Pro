<table border="0" cellpadding="0" cellspacing="0">
  <tr>
    <td width="96" valign="middle"><img src="assets/icon.svg" alt="PyChem-Pro logo" width="80" height="80"></td>
    <td valign="middle"><h1>PyChem-Pro</h1></td>
  </tr>
</table>

> A pure-Python desktop application and library for chemistry and cheminformatics — molecular visualization, SMILES parsing, PDB loading, MMFF94 geometry optimization, descriptors, and a plugin ecosystem. No RDKit. No OpenBabel. Everything implemented from scratch.

**Repository:** https://github.com/vijaymasand/PyChem-Pro

[![Buy Me A Coffee](https://www.buymeacoffee.com/assets/img/custom_images/yellow_img.png)](https://buymeacoffee.com/vijaymasand)

---

## Table of Contents

- [About PyChem-Pro](#about-pychem-pro)
- [Key Features](#key-features)
- [Architecture Overview](#architecture-overview)
- [Tech Stack](#tech-stack)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Running the Application](#running-the-application)
  - [Recommended: One-Click Launchers](#recommended-one-click-launchers)
  - [Manual Setup (all platforms)](#manual-setup-all-platforms)
- [Quick Start — Python API](#quick-start--python-api)
- [Project Structure](#project-structure)
- [Services Layer (Public Protocols)](#services-layer-public-protocols)
- [MMFF94 Force Field](#mmff94-force-field)
- [Rendering Pipeline](#rendering-pipeline)
- [Multiprocessing](#multiprocessing)
- [Plugins](#plugins)
- [Development Workflow](#development-workflow)
- [Contributing](#contributing)
- [Cross-Platform Notes](#cross-platform-notes)
- [Performance Targets](#performance-targets)
- [Roadmap](#roadmap)
- [License](#license)
- [Citation](#citation)
- [Acknowledgments](#acknowledgments)

---

## About PyChem-Pro

PyChem-Pro is a desktop chemistry application and Python library that combines molecular visualization (like PyMOL) with cheminformatics primitives (like RDKit). Unlike most tools in the space, PyChem-Pro is **pure Python with NumPy**. There is no C++ extension, no dependency on any cheminformatics library. Every feature — SMILES parsing, 3D coordinate generation, force field optimization, descriptor calculation, Shrake-Rupley SASA, ring perception, protein cartoon rendering — is implemented from scratch and readable end-to-end.

This makes PyChem-Pro:

- **Educational.** Graduate students can open any file and see how a real MMFF94 optimizer or a Catmull-Rom ribbon renderer actually works.
- **Portable.** No compilers. No system libraries beyond Python and Qt. Runs on Windows, macOS, and Linux identically.
- **Extensible.** A service-oriented architecture with Protocol interfaces lets researchers swap the force field, the renderer, or the file loader without touching the rest of the codebase.
- **Academic-friendly.** Intended for inclusion in university cheminformatics curricula and as an open-source reference implementation.

## Key Features

### Chemistry
- SMILES parser with Huckel aromaticity perception, stereochemistry, bracketed atoms, radicals
- PDB / MOL / MOL2 / SDF file readers with automatic bond detection and CONECT parsing
- 3D coordinate generation from molecular topology (BFS + force-field relaxation)
- **MMFF94 force field** with hydrogen addition, bond stretching, angle bending, torsion, Van der Waals, and BCI partial charge assignment
- L-BFGS and steepest descent optimizers
- Gasteiger partial charges, AM1 and PM3 semi-empirical charge models
- Molecular descriptor calculator (constitutional, topological, electronic, geometric, hybrid, quantum, fingerprints)
- Morgan (ECFP) and topological fingerprints
- Center of Mass and Shrake-Rupley SASA

### Visualization
- Hardware-aware software 3D rendering (QPainter with gradient sphere shading)
- Ball-and-stick, space-fill, wireframe, cartoon, ribbon, and backbone modes
- PyMOL-style protein cartoons with alpha-helix / beta-sheet / coil detection
- Color-by secondary structure, rainbow, chain, B-factor
- 2D chemical structure viewer with single/double/triple/aromatic bond rendering, wedge/hash stereo bonds, atom labels, formal charges
- Natural-language atom selection (`sele('organic')`, `sele('within 5.0 COM')`, `sele('chain A and helix')`)
- High-DPI image export (72 / 150 / 300 / 600 DPI)
- Offline ray-tracer for publication-quality images
- Print (Ctrl+P) — prints 2D and 3D views on a single page

### Architecture
- Service-oriented with `typing.Protocol` interfaces
- EventBus for decoupled pub/sub communication
- ParallelExecutor at 50% CPU cores across file loading, force field, rendering, descriptors, ray-tracing, and conformer search
- Public `pychem` Python package usable in Jupyter notebooks without PySide6

### Plugin System
- `BasePlugin` / `PluginWidget` API for custom analysis, visualization, I/O, and utility plugins
- Built-in templates for analysis, visualization, and I/O plugin types
- QSAR modeling, Ramachandran, descriptor pruning, docking pose visualization, molecular weight calculator — all shipped as plugins

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                     PUBLIC API LAYER (pychem/)                       │
│  No Qt dependency. Jupyter-friendly.                                 │
│  import pychem; pychem.parse_smiles("CCO"); pychem.optimize(mol)     │
└──────────────────────────────┬──────────────────────────────────────┘
                               │
┌──────────────────────────────▼──────────────────────────────────────┐
│                     SERVICE REGISTRY (src/core/registry.py)          │
│  Plain Python factory, no DI framework.                              │
│  registry.forcefield, registry.loader, registry.descriptors, ...     │
└───┬──────────┬──────────┬──────────┬──────────┬──────────┬─────────┘
    │          │          │          │          │          │
┌───▼───┐ ┌───▼───┐ ┌───▼───┐ ┌───▼───┐ ┌───▼───┐ ┌───▼────┐
│Force  │ │Render │ │Loader │ │Coord  │ │Descr  │ │Plugin  │
│Field  │ │Service│ │Service│ │Gen    │ │Calc   │ │Manager │
└───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘ └───┬────┘
    │         │         │         │         │         │
┌───▼─────────▼─────────▼─────────▼─────────▼─────────▼──────────────┐
│                     CORE DOMAIN (src/core/domain/)                   │
│  Molecule, Atom, Bond, Element — zero external dependencies          │
└─────────────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────────────┐
│                     INFRASTRUCTURE (src/core/)                       │
│  EventBus, ParallelExecutor, Protocols, Performance Profiler         │
└─────────────────────────────────────────────────────────────────────┘
```

**Dependency rules (enforced by convention):**

1. Core domain imports nothing from services, features, app, or Qt.
2. Infrastructure imports only from core domain.
3. Services import from core domain + protocols. Never from other service implementations.
4. Features import from services (via protocols) + core domain + Qt.
5. App imports from features + services + core domain + Qt.
6. Public API imports from services + core domain. Never from Qt.
7. Plugins import from public API + Qt compat layer. Never from service internals.

Every subsystem implements a `typing.Protocol` interface (`IForceField`, `IRenderer`, `ILoader`, `ICoordinateGenerator`, `IDescriptorCalculator`) so the implementation can be swapped without touching consumers.

---

## Tech Stack

| Layer | Technology |
|-------|------------|
| Language | Python 3.10+ (3.13/3.14 supported) |
| GUI Framework | PySide6 (Qt 6.5+) |
| Numerical | NumPy 1.24+ |
| Multiprocessing | Python `concurrent.futures.ProcessPoolExecutor` with forced `spawn` context for cross-platform safety |
| Bundler (optional) | Nuitka 2.0+ for producing standalone binaries |
| Plugin stack (optional) | pandas, scipy, scikit-learn, matplotlib, packaging — only needed if you run plugins that require them |

**No chemistry dependencies.** PyChem-Pro does not use any cheminformatics or external chemistry library.

### Vendored

- **OASA** (Open Structure Access) — a subset is vendored under `src/vendors/oasa/` and used for SMILES/InChI utilities and 2D coordinate generation fallbacks. It is treated as frozen third-party code and not modified.

---

## System Requirements

- **Operating System:** Windows 10/11, macOS 12+, or any Linux distribution with Qt 6 support (Ubuntu 22.04+, Fedora 37+, Arch, etc.)
- **Python:** 3.10 or newer
- **RAM:** 4 GB minimum, 8 GB recommended for large proteins (>5000 atoms)
- **Display:** OpenGL 3.3+ optional (used by the future hardware-accelerated renderer path). The default QPainter renderer works on any display.
- **CPU:** Multi-core recommended. The `ParallelExecutor` uses 50% of available cores automatically.

### Cross-Platform Notes

PyChem-Pro is designed to work identically on Windows, macOS, and Linux. Two design decisions enforce this:

1. **Multiprocessing:** `ProcessPoolExecutor` is always created with the `spawn` start method (`multiprocessing.get_context('spawn')`). This avoids `fork()`-related crashes on macOS (Qt + CoreFoundation assertion) and matches the only option available on Windows.
2. **File paths:** All path handling uses `os.path` or `pathlib` — no hard-coded separators.
3. **Fonts:** The UI uses a best-available fallback sequence and will warn but not fail if a platform-specific font is missing.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/vijaymasand/PyChem-Pro.git
cd PyChem-Pro
```

### 2. Create a virtual environment

```bash
# macOS / Linux
python3 -m venv venv
source venv/bin/activate

# Windows (cmd)
python -m venv venv
venv\Scripts\activate.bat

# Windows (PowerShell)
python -m venv venv
venv\Scripts\Activate.ps1
```

### 3. Install dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

This installs the runtime dependencies:
- `PySide6 >= 6.5` (Qt GUI framework)
- `numpy >= 1.24`, `matplotlib >= 3.7`, `pandas >= 2.0`, `scipy >= 1.10`
- `scikit-learn >= 1.3`, `rdkit >= 2023.3` (used by QSAR and descriptor plugins)
- `psutil`, `pillow`

Build tools (Nuitka, dmgbuild) are in `requirements-build.txt` and only needed if you want to compile a standalone binary — not required for running from source.

### 4. Optional — install plugin dependencies

The built-in plugins in `plugins/` have optional dependencies. Install them only if you plan to use those plugins:

```bash
pip install packaging pandas scipy scikit-learn matplotlib rdkit
```

Without these, the core application runs fine; the affected plugins are simply skipped at startup with a warning.

---

## Running the Application

### Recommended: One-Click Launchers

The easiest way to run PyChem-Pro is to double-click the launcher for your platform. No terminal, no Python knowledge, no setup steps required.

| Platform | File to double-click |
|----------|----------------------|
| macOS | `PyChem.command` (in Finder) |
| Windows | `PyChem.bat`, then use the `PyChem - Launch.lnk` shortcut created on first run |

Both launchers handle everything automatically and follow the same steps:

| Step | What happens |
|------|-------------|
| 1 | Checks whether Python 3.10+ is already installed |
| 2 | **If Python is missing** — installs it automatically: macOS uses Homebrew (asks for your password once); Windows uses `winget` (built into Windows 10/11), falling back to a silent download from python.org |
| 3 | Creates a PyChem-Pro icon shortcut for easy access (first run only) |
| 4 | Creates an isolated Python virtual environment inside the project folder (`venv/`) so no system packages are touched (first run only) |
| 5 | Installs all required packages from `requirements.txt` into that venv (first run only) |
| 6 | Launches PyChem |

> **First run** takes 2–5 minutes while packages download and install. Every run after that starts in seconds — the venv is reused as-is.

> The launchers are provided for your convenience and are the recommended way to get started. If you prefer to manage your Python environment yourself, or if a launcher does not work on your machine, use the manual setup below — it is straightforward and gives you full control.

---

### Manual Setup (all platforms)

#### 1. Install Python 3.10+

- **macOS:** Download from [python.org/downloads](https://www.python.org/downloads/) or run `brew install python@3.12`
- **Linux:** `sudo apt install python3` (Ubuntu/Debian) or `sudo dnf install python3` (Fedora)
- **Windows:** Download from [python.org/downloads](https://www.python.org/downloads/) — check "Add Python to PATH" during install

Verify:
```bash
python3 --version   # should print 3.10 or newer
```

#### 2. Clone the repository

```bash
git clone https://github.com/vijaymasand/PyChem-Pro.git
cd PyChem-Pro
```

#### 3. Create a virtual environment

```bash
# macOS / Linux
python3 -m venv venv
source venv/bin/activate

# Windows (cmd)
python -m venv venv
venv\Scripts\activate.bat
```

#### 4. Install dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

#### 5. Run PyChem-Pro

```bash
python main.py
```

The first launch generates a 10-year development license automatically. No network access required.

### Expected startup output

```
Using PySide6 framework
Plugin system initialized successfully
```

---

## Quick Start — Python API

PyChem-Pro provides a public package at the top level (`pychem/`) that does not require PySide6 to be imported. You can use it in Jupyter notebooks or scripts.

```python
import pychem

# Parse SMILES
mol = pychem.parse_smiles("CCO")  # ethanol
print(mol.molecular_formula())    # C2H6O

# Generate 3D coordinates
pychem.generate_3d(mol)

# Full MMFF94 optimization (adds hydrogens, assigns BCI charges, minimizes)
result = pychem.optimize(mol, max_iters=500)
print(f"Converged: {result.converged}, E = {result.final_energy:.2f} kcal/mol")
print(f"Steps: {result.num_steps}")

# Partial charges are now assigned
for atom in mol.atoms:
    print(f"  {atom.symbol}{atom.index}: q = {atom.partial_charge:+.4f}")

# Molecular descriptors
desc = pychem.descriptors(mol)
print(desc)
# {'molecular_weight': 46.07, 'num_atoms': 9, 'num_bonds': 8,
#  'num_heavy_atoms': 3, 'formula': 'C2H6O', 'num_rings': 0, 'total_charge': 0}

# Batch descriptor calculation (uses ParallelExecutor)
mols = [pychem.parse_smiles(s) for s in ["CCO", "c1ccccc1", "CC(=O)O", "CCCCCC"]]
batch = pychem.descriptors_batch(mols)

# Load a protein from PDB
protein = pychem.load("1AKE.pdb")
print(f"{protein.num_atoms} atoms, {protein.num_bonds} bonds")

# Add explicit hydrogens only (no optimization)
h_count = pychem.add_hydrogens(mol)

# Compute MMFF94 BCI charges only
pychem.compute_charges(mol)
```

### Available API functions

```python
pychem.parse_smiles(smiles: str) -> Molecule
pychem.load(path: str) -> Molecule
pychem.generate_3d(mol: Molecule, optimize: bool = True, max_steps: int = 200) -> None
pychem.optimize(mol: Molecule, max_iters: int = 500, method: str = 'lbfgs') -> OptimizationResult
pychem.add_hydrogens(mol: Molecule) -> int
pychem.compute_charges(mol: Molecule) -> None
pychem.descriptors(mol: Molecule, names: list[str] | None = None) -> dict
pychem.descriptors_batch(mols: list[Molecule], names: list[str] | None = None) -> list[dict]
```

---

## Project Structure

```
PyChem/
├── pychem/                        # Public API (no Qt dependency)
│   ├── __init__.py                #   import pychem
│   ├── api.py                     #   public facade functions
│   └── _bridge.py                 #   ServiceRegistry singleton
│
├── src/
│   ├── core/                      # Infrastructure and core domain
│   │   ├── domain/models/         # Molecule, Atom, Bond, Element
│   │   ├── protocols/             # IForceField, IRenderer, ILoader, ...
│   │   ├── events.py              # EventBus + event dataclasses
│   │   ├── parallel.py            # ParallelExecutor (50% CPU cores)
│   │   ├── registry.py            # ServiceRegistry
│   │   ├── performance/           # Profiler, parallel loader
│   │   └── security/              # License manager
│   │
│   ├── services/                  # Service implementations
│   │   ├── forcefield/            # MMFF94Service, HydrogenAdder,
│   │   │                          # AngleBending, Torsion, parameters
│   │   ├── rendering/             # RendererFactory, parallel_projection
│   │   ├── loading/               # LoaderService
│   │   ├── coordinates/           # CoordinateGeneratorService
│   │   └── descriptors/           # DescriptorService
│   │
│   ├── app/                       # GUI application (thin shell)
│   │   ├── main_window.py         # QMainWindow wiring
│   │   ├── menu_bar.py            # Menu construction
│   │   ├── toolbar.py             # Toolbar construction
│   │   ├── file_operations.py     # Open / Save / Export / Print
│   │   ├── chemistry_actions.py   # Optimize / Charges / Descriptors
│   │   ├── viewer_coordinator.py  # View toggles, COM/centroid spheres
│   │   ├── molecule_controller.py # Signal wiring, undo/redo, selection
│   │   ├── conversion_worker.py   # QThread SMILES->3D worker
│   │   ├── plugin_interface.py    # Plugin browser UI
│   │   ├── plugin_card.py
│   │   └── plugin_installer.py
│   │
│   ├── features/                  # Feature modules (split by domain)
│   │   ├── visualization_3d/      # MolViewer3D + painter_renderer +
│   │   │                          # mouse_controller + protein_rendering
│   │   ├── visualization_2d/      # MolViewer2D + bond/atom renderers
│   │   ├── cheminformatics/       # AM1, PM3, Gasteiger, MMFF94 legacy
│   │   ├── layout_3d/             # 3D coordinate generator
│   │   ├── layout_2d/             # 2D layout generators
│   │   ├── smiles_parser/         # OpenSMILES parser
│   │   ├── smiles_generator/      # SMILES writer
│   │   ├── io/                    # File readers/writers
│   │   ├── descriptor_calculator/ # Descriptor GUI and engines
│   │   ├── scripting_console/     # Python REPL + atom selection
│   │   ├── control_panel/         # Input panel
│   │   └── ui/                    # Dialogs (colors, spheres, etc.)
│   │
│   ├── plugins/                   # Plugin manager infrastructure
│   ├── shared/                    # Qt compat, theme
│   └── vendors/oasa/              # Vendored OASA library (frozen)
│
├── plugins/                       # Built-in and user plugins
│   ├── docking_pose_visualizer.py
│   ├── ramachandran_plugin.py
│   ├── qsar_modeler_plugin.py
│   ├── mol_weight_calculator.py
│   └── ...
│
├── tests/                         # Unit tests
├── testing/                       # Development / debugging scripts
├── docs/                          # Architecture specs, user guides
│   └── superpowers/
│       ├── specs/                 # Design documents
│       └── plans/                 # Implementation plans
├── main.py                        # Application entry point
├── build.py                       # Nuitka bundler script
├── requirements.txt
└── README.md                      # You are here
```

---

## Services Layer (Public Protocols)

All major subsystems are exposed through `typing.Protocol` interfaces so alternative implementations can be dropped in without touching consumers.

```python
from src.core.protocols import (
    IForceField, IRenderer, ILoader,
    ICoordinateGenerator, IDescriptorCalculator,
    OptimizationResult, Camera, RenderMode,
)
```

### Example: registering a custom force field

```python
from src.core.protocols.forcefield import IForceField, OptimizationResult
from src.core.domain.models.molecule import Molecule

class MyForceField:
    """Custom force field implementing the IForceField protocol."""

    def add_hydrogens(self, mol: Molecule) -> int:
        ...

    def assign_atom_types(self, mol: Molecule) -> None:
        ...

    def assign_charges(self, mol: Molecule) -> None:
        ...

    def optimize_geometry(self, mol, max_iters=500, convergence=1e-4,
                          method='lbfgs') -> OptimizationResult:
        ...

    def compute_energy(self, mol: Molecule) -> float:
        ...

# Install into the registry at startup
from src.core.registry import ServiceRegistry
registry = ServiceRegistry()
registry.forcefield = MyForceField()
```

The class does not need to inherit from anything — `typing.Protocol` uses structural subtyping.

---

## MMFF94 Force Field

PyChem-Pro ships a pure-Python implementation of the Merck Molecular Force Field (MMFF94). At present, it is not perfect. It is simplified relative to the full 75K+ rule set but covers standard organic chemistry:

| Term | Implementation |
|------|----------------|
| Bond stretching | Hookean with MMFF94 `r0` and `kb` parameters |
| Angle bending | Cubic-corrected `0.021914 ka (theta - theta0)^2 [1 + cb (theta - theta0)]` with analytical gradient |
| Torsion (dihedral) | Fourier `V1 cos(phi) + V2 cos(2 phi) + V3 cos(3 phi)` with analytical gradient |
| Van der Waals | Lennard-Jones 12-6 with 1-2 and 1-3 exclusions, parallelized for >200 pairs |
| Partial charges | Bond Charge Increment (BCI) method |
| Hydrogen addition | Hybridization-aware 3D placement (tetrahedral / trigonal / linear) |

### Optimization pipeline

1. Assign hybridization (sp / sp2 / sp3)
2. Add explicit hydrogens with ideal 3D positions
3. Assign BCI partial charges
4. Build interaction lists (bonds, angles, torsions, VdW pairs)
5. Minimize energy via steepest descent or L-BFGS with Armijo line search

Verified analytical gradients against numerical finite differences within `1e-4` for angles and `1e-3` for torsions.

### Files

```
src/services/forcefield/
├── mmff94_service.py       # Unified IForceField implementation
├── hydrogen.py             # HydrogenAdder with hybridization-aware placement
├── angle_bending.py        # AngleBendingCalculator
├── torsion.py              # TorsionCalculator
└── parameters.py           # Consolidated bond/angle/torsion/VdW/BCI tables
```

---

## Rendering Pipeline

The 3D viewer (`src/features/visualization_3d/ui/mol_viewer_3d.py`) paints directly via QPainter. Rendering logic is extracted into `painter_renderer.py` (the engine) and `mouse_controller.py` (interaction).

### Optimizations

- **Gradient cache.** `QRadialGradient` color stops are cached by `(element, radius_bucket, is_hovered, use_ssao, depth_bucket)`. Rebuilding a positioned gradient from cached stops is free; the expensive color arithmetic runs once per unique atom type.
- **Off-screen culling.** Atoms and bonds outside the viewport plus a 50-100 pixel margin are skipped before any draw call.
- **LOD (Level of Detail).** Atoms projected to less than 2 px radius are drawn as plain filled circles instead of gradient spheres.
- **Parallel pre-render.** For large molecules, projection, depth sort, and visibility culling can be split across `ParallelExecutor` workers (in `src/services/rendering/parallel_projection.py`). Draw calls themselves remain on the main thread per Qt's threading model.

### Protein cartoon rendering

`src/features/visualization_3d/services/protein_rendering.py` implements:

- Simplified DSSP-style secondary structure detection (~85% accuracy on standard tests)
- Catmull-Rom spline smoothing for ribbon paths
- PyMOL-style cartoon tubes for helices, flat ribbons with arrow heads for sheets, thin coils for loops
- Color schemes: secondary structure, rainbow, by chain, by B-factor

---

## Multiprocessing

A single `ParallelExecutor` at `src/core/parallel.py` wraps `concurrent.futures.ProcessPoolExecutor` and is shared across all services. It uses **50% of available CPU cores** by default.

### Forced spawn start method

```python
import multiprocessing as mp
_mp_context = mp.get_context('spawn')
```

Cross-platform rationale:
- **macOS:** `fork()` + Qt causes CoreFoundation assertions. `spawn` avoids it.
- **Windows:** `fork()` is not available. `spawn` is the only option.
- **Linux:** `fork()` works, but `spawn` is safer with Qt and matches other platforms.

### Where multiprocessing is used

| Service | What is parallelized |
|---------|----------------------|
| File loading | PDB / MOL2 atom record parsing split into N chunks |
| MMFF94 force field | VdW pairwise computation split into N chunks (only when >200 pairs) |
| Coordinate generation | N conformers with independent random seeds, return lowest energy |
| Descriptor calculation | Batch molecules distributed across workers |
| 3D rendering | Projection, depth sort, culling split across workers (large molecules only) |
| Ray-tracing export | 64x64 pixel tiles processed in parallel |

All worker functions are **module-level** (required by the `spawn` start method for pickling).

---

## Plugins

PyChem-Pro has a first-class plugin system. Plugins appear as dockable panels in the GUI and can subscribe to molecule-changed events.

### Minimal plugin

```python
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType
from src.shared.qt_compat import QVBoxLayout, QPushButton, QLabel

class HelloPluginWidget(PluginWidget):
    def setup_ui(self):
        layout = QVBoxLayout(self.widget)
        self.label = QLabel("No molecule loaded")
        layout.addWidget(self.label)

    def on_molecule_changed(self, molecule):
        if molecule:
            self.label.setText(f"{molecule.num_atoms} atoms")

class HelloPlugin(BasePlugin):
    def get_info(self) -> PluginInfo:
        return PluginInfo(
            name="Hello Plugin",
            version="1.0.0",
            description="Shows atom count",
            author="Your Name",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[],
        )

    def create_widget(self):
        return HelloPluginWidget(self)

    def initialize(self):
        return True

    def cleanup(self):
        pass

    def on_molecule_changed(self, molecule):
        if hasattr(self, 'widget') and self.widget:
            self.widget.on_molecule_changed(molecule)
```

Drop the file into `plugins/` and restart PyChem, or use **Plugins → Installed Plugins...** to load it manually.

### Built-in plugins

| Plugin | Purpose |
|--------|---------|
| `mol_weight_calculator.py` | Molecular weight and element composition |
| `ramachandran_plugin.py` | Ramachandran plot for protein backbones |
| `qsar_modeler_plugin.py` | QSAR model building and validation |
| `qsar_rfa_mars_plugin.py` | Random Forest + MARS QSAR |
| `docking_pose_visualizer.py` | Docking pose comparison |
| `DescriptorPruningApp.py` | Descriptor correlation pruning |
| `example_analysis_plugin.py` | Template for analysis plugins |
| `example_visualization_plugin.py` | Template for visualization plugins |

Optional plugins declare their external dependencies in `PluginInfo.dependencies`. If a required package is missing, the plugin is skipped with a warning at startup.

---

## Development Workflow

### Coding standards

- **Style.** Follow PEP 8, 4-space indentation, max 100 chars per line.
- **Typing.** Use `typing.Protocol` for service interfaces, dataclasses for DTOs, and type hints everywhere reasonable.
- **No chemistry dependencies.** NumPy is the only permitted numerical dependency.
- **Module-level functions for multiprocessing.** Any function passed to `ParallelExecutor.map` must be defined at module scope (required by the `spawn` start method).
- **Cross-platform paths.** Use `os.path` or `pathlib`. Never hard-code `/` or `\`.
- **Large files.** Keep individual Python files under ~800 lines. If a file grows beyond that, split it into focused modules.

### Commit conventions

- Use imperative present-tense commit subjects: `feat: add X`, `fix: correct Y`, `refactor: split Z`, `perf: cache W`, `docs: explain V`.
- First line 72 chars max. Body wrapped at 72 chars.
- Do not add `Co-Authored-By: Claude` or any AI attribution.
- Co-authorship for human collaborators is welcome via the standard `Co-Authored-By: Name <email>` trailer.

### Running the application during development

```bash
source venv/bin/activate
python main.py
```

Startup profiling is enabled by default and prints timing information for Qt imports, MainWindow creation, and total startup time.

### Building a standalone binary (optional)

```bash
python build.py
```

This uses Nuitka to produce a compiled executable in `build/`. Windows and macOS builds are tested; Linux builds work but are not packaged for distribution.

---

## Contributing

Contributions are welcome — bug reports, feature requests, pull requests, new plugins, documentation improvements, and test coverage.

### Before submitting a pull request

1. Run the full unit test suite and make sure nothing regresses:
   ```bash
   for t in tests/test_*.py; do python3 "$t"; done
   ```
2. Launch the GUI and verify the affected workflow manually.
3. Update relevant documentation in `docs/` if you changed an interface.
4. Follow the commit conventions above.

### Reporting bugs

Please include:
- Operating system and version
- Python version (`python --version`)
- PySide6 version (`pip show PySide6`)
- Full traceback if there is an error
- Minimal reproducer (SMILES string, PDB id, or attached file)

---

## Performance Targets

Approximate targets measured on a 4-core Intel machine with 8 GB RAM. Actual numbers depend on the molecule.

| Metric | Target |
|--------|--------|
| MMFF94 ethanol optimization | ~0.3 s (full force field with H addition) |
| MMFF94 insulin (51 residues) | <30 s |
| 3D render FPS (100 atoms) | ~60 FPS |
| 3D render FPS (5000 atoms) | ~30 FPS |
| 3D render FPS (10K atoms) | ~60 FPS on GL path (future) |
| PDB load (10K atoms) | ~1 s (parallel) |
| Descriptor calculation (100 mols) | ~50% faster than sequential via ParallelExecutor |
| Startup time | ~1.5 s |

---

## Roadmap

The architecture spec at `docs/superpowers/specs/2026-04-11-architecture-redesign-design.md` captures the long-term plan. Short-term items:

- Hardware-accelerated OpenGL renderer with protein cartoon support
- Stereochemistry enforcement in 3D coordinate generation (R/S, E/Z)
- Ring template library for faster and more accurate cyclic structure generation
- Conformer ensemble analysis (RMSD clustering, Boltzmann-weighted properties)
- In-app 2D chemical structure editor / sketcher
- Extended MMFF94 parameter tables (halogens, phosphorus, metals)
- Docking interface (scoring function only, no compiled dependencies)

See `docs/` for detailed specs and implementation plans.

---

## License

PyChem-Pro is licensed under the **[Polyform Noncommercial License 1.0.0](LICENSE)**.

Copyright (c) 2026 Vijay Masand and Gaurav Masand

- **Personal and educational use** — free to use, modify, and distribute.
- **Commercial use** — strictly prohibited. This includes using PyChem as part of a paid service, distributing it for a fee, or relying on it to generate revenue.

See the `LICENSE` file for the full legal text.

Vendored OASA code retains its original BSD-style license.

---

## Citation

If you use PyChem-Pro in academic work, please cite:

```
@software{pychem,
  title   = {PyChem-Pro: A Pure-Python Cheminformatics and Molecular Visualization Toolkit},
  author  = {Masand, Gaurav and Masand, Vijay},
  year    = {2026},
  url     = {https://github.com/vijaymasand/PyChem},
  note    = {Accessed: YYYY-MM-DD}
}
```

A proper DOI will be issued once the first tagged release is published.

---

## Acknowledgments

- The **MMFF94 force field** is the work of Thomas A. Halgren (Merck, 1996). Our implementation is a simplified pure-Python version of the published specification.
- The **OASA** library (Open Structure Access, Beda Kosata) is vendored under `src/vendors/oasa/` and provides SMILES/InChI utilities and 2D layout primitives.
- **PyMOL** and **Jmol** inspired the protein cartoon representations.
- **PySide6 / Qt** from the Qt Project provides the GUI framework.

---

## Contact

- GitHub: https://github.com/vijaymasand/PyChem-Pro
- Issues: https://github.com/vijaymasand/PyChem-Pro/issues

Questions, feedback, and pull requests are welcome.
