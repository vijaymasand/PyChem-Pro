# PyChem-Pro — pdbx/__init__.py
# Migrated from Python 2 to Python 3 by the PyChem-Pro project (2026-06-21).
#
"""
pdbx — PDBx/mmCIF reader, parser, and writer for PyChem-Pro.

Public API
----------
::

    from pdbx import read_mmcif, write_mmcif, read_mmcif_models

    # Read a .cif / .mmcif file into a PyChem-Pro Molecule
    mol = read_mmcif("1abc.cif")

    # Write a Molecule back to mmCIF
    write_mmcif(mol, "out.cif")

    # Read all NMR models from a multi-model CIF
    models = read_mmcif_models("nmr_ensemble.cif")

Lower-level container API (for advanced use):
::

    from pdbx import PdbxReader, PdbxWriter, DataContainer, DataCategory
"""

from pdbx.mmcif_molecule import read_mmcif, write_mmcif, read_mmcif_models
from pdbx.reader.PdbxReader import PdbxReader
from pdbx.writer.PdbxWriter import PdbxWriter
from pdbx.reader.PdbxContainers import (
    DataContainer,
    DataCategory,
    DefinitionContainer,
    CifName,
)

__all__ = [
    # High-level Molecule bridge
    "read_mmcif",
    "write_mmcif",
    "read_mmcif_models",
    # Low-level container API
    "PdbxReader",
    "PdbxWriter",
    "DataContainer",
    "DataCategory",
    "DefinitionContainer",
    "CifName",
]
