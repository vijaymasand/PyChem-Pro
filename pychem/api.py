"""
Public facade functions for PyChem.

These thin wrappers delegate to the ServiceRegistry.
They are the stable public API surface.

mmCIF/PDBx support was added in PyChem-Pro v1.0 (2026-06-21).
Use ``read_mmcif``, ``write_mmcif``, ``read_mmcif_models``, and the
generic ``export()`` dispatcher for PDBx/mmCIF I/O.
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Optional, List

if TYPE_CHECKING:
    from src.core.domain.models.molecule import Molecule


def load(path: str, parallel: bool = True) -> "Molecule":
    """Load a molecular file (PDB, MOL, MOL2, SDF, CIF, mmCIF).

    Uses parallel chunk parsing for large files by default.

    Args:
        path:     Path to the molecular file.  Supported extensions:
                  ``.pdb``, ``.ent``, ``.mol``, ``.mol2``, ``.sdf``, ``.sd``,
                  ``.cif``, ``.mmcif``.
        parallel: Allow parallel processing for large files.

    Returns:
        Molecule instance.
    """
    from pychem._bridge import get_registry
    return get_registry().loader.load(path, parallel=parallel)


def parse_smiles(smiles: str):
    """Parse a SMILES string into a Molecule."""
    from src.features.smiles_parser.services.parser import parse_smiles as _parse
    return _parse(smiles)


def generate_3d(mol, optimize: bool = True, max_steps: int = 200) -> None:
    """Generate 3D coordinates in-place on the molecule."""
    from src.features.layout_3d import generate_3d_coordinates
    generate_3d_coordinates(mol, optimize=optimize, max_opt_steps=max_steps)


def optimize(mol, max_iters: int = 500, convergence: float = 1e-4, method: str = 'lbfgs'):
    """Optimize molecular geometry using MMFF94. Returns OptimizationResult."""
    from pychem._bridge import get_registry
    return get_registry().forcefield.optimize_geometry(
        mol, max_iters=max_iters, convergence=convergence, method=method
    )


def compute_charges(mol) -> None:
    """Assign MMFF94 partial charges (BCI method) in-place on the molecule."""
    from pychem._bridge import get_registry
    get_registry().forcefield.assign_charges(mol)


def add_hydrogens(mol) -> int:
    """Add explicit hydrogens with 3D positions. Returns count added."""
    from pychem._bridge import get_registry
    return get_registry().forcefield.add_hydrogens(mol)


def descriptors(mol, names=None) -> dict:
    """Calculate molecular descriptors."""
    from pychem._bridge import get_registry
    return get_registry().descriptors.calculate(mol, descriptor_names=names)


def descriptors_batch(molecules, names=None) -> list:
    """Calculate molecular descriptors for multiple molecules in parallel."""
    from pychem._bridge import get_registry
    return get_registry().descriptors.calculate_batch(molecules, descriptor_names=names)


# ── mmCIF / PDBx specific functions ──────────────────────────────────────────

def read_mmcif(path: str, model: int = 1) -> "Molecule":
    """Read a PDBx/mmCIF file and return a Molecule.

    Args:
        path:   Path to a ``.cif`` or ``.mmcif`` file.
        model:  Model number to load (default=1; relevant for NMR ensembles).

    Returns:
        Molecule instance with atoms, bonds, secondary structure, and
        metadata stored in ``molecule.properties``.
    """
    from pdbx.mmcif_molecule import read_mmcif as _read
    return _read(path, model=model)


def read_mmcif_string(cif_text: str, model: int = 1) -> "Molecule":
    """Parse CIF text (string) and return a Molecule.

    Args:
        cif_text: Raw PDBx/mmCIF file content as a string.
        model:    Model number to load (default=1).

    Returns:
        Molecule instance.
    """
    from pdbx.mmcif_molecule import read_mmcif as _read
    return _read(cif_text, from_string=True, model=model)


def read_mmcif_models(path: str) -> List["Molecule"]:
    """Read all models from a multi-model mmCIF file (e.g. NMR ensembles).

    Args:
        path: Path to a ``.cif`` or ``.mmcif`` file.

    Returns:
        List of Molecule objects, one per model number.
    """
    from pdbx.mmcif_molecule import read_mmcif_models as _read_all
    return _read_all(path)


def write_mmcif(mol: "Molecule", filepath: Optional[str] = None,
                entry_id: Optional[str] = None) -> str:
    """Write a Molecule to PDBx/mmCIF format.

    Args:
        mol:        The Molecule to serialise.
        filepath:   If given, write to this path (UTF-8).
        entry_id:   PDB-style 4-character entry ID (e.g. ``"1ABC"``).

    Returns:
        CIF content as a string (also writes to *filepath* if given).
    """
    from pdbx.mmcif_molecule import write_mmcif as _write
    return _write(mol, filepath=filepath, entry_id=entry_id)


def export(mol: "Molecule", filepath: str, **kwargs) -> Optional[str]:
    """Export a Molecule to a file, dispatching by extension.

    Supported formats: ``.mol``, ``.sdf``, ``.mol2``, ``.cif``, ``.mmcif``.

    Args:
        mol:      The Molecule to export.
        filepath: Output file path — the extension determines the format.
        **kwargs: Additional keyword arguments forwarded to the format writer.

    Returns:
        File content as a string for text formats.

    Raises:
        ValueError: If the file extension is not supported.
    """
    import os
    ext = os.path.splitext(filepath)[1].lower()

    if ext == ".mol":
        from src.features.io.exporters.sdf_writer import write_mol
        return write_mol(mol, filepath=filepath)
    elif ext in (".sdf", ".sd"):
        from src.features.io.exporters.sdf_writer import write_sdf
        return write_sdf(mol, filepath=filepath, **kwargs)
    elif ext == ".mol2":
        from src.features.io.exporters.mol2_writer import write_mol2
        return write_mol2(mol, filepath=filepath)
    elif ext in (".cif", ".mmcif"):
        from src.features.io.exporters.mmcif_exporter import write_mmcif_molecule
        return write_mmcif_molecule(mol, filepath=filepath, **kwargs)
    else:
        raise ValueError(f"Unsupported export format: {ext}")
