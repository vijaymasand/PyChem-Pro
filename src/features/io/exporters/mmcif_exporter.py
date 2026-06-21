"""
mmCIF / PDBx Exporter — Write molecules in PDBx/mmCIF format.

Follows the same pattern as ``mol2_writer.py`` and ``sdf_writer.py``:
a thin wrapper around the core bridge in ``pdbx.mmcif_molecule``.

Usage::

    from src.features.io.exporters.mmcif_exporter import write_mmcif_molecule

    cif_text = write_mmcif_molecule(mol)                    # -> str
    write_mmcif_molecule(mol, filepath="output.cif")        # -> str (also writes file)
    write_mmcif_molecule(mol, filepath="out.cif", entry_id="4XYZ")
"""

from __future__ import annotations

from typing import Optional


def write_mmcif_molecule(molecule, filepath: Optional[str] = None,
                         entry_id: Optional[str] = None) -> str:
    """Write *molecule* to PDBx/mmCIF format.

    Args:
        molecule:  A PyChem-Pro ``Molecule`` instance.
        filepath:  Optional output path.  If given, the file is written in
                   UTF-8.  The returned string always contains the CIF text.
        entry_id:  PDB-style 4-character entry ID.  Defaults to
                   ``molecule.name`` (upper-cased, trimmed to 4 chars).

    Returns:
        CIF content as a string.
    """
    from pdbx.mmcif_molecule import write_mmcif as _write
    return _write(molecule, filepath=filepath, entry_id=entry_id)
