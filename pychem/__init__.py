"""
PyChem — Pure Python cheminformatics toolkit.

Public API for use in Jupyter notebooks and scripts.
No PySide6 dependency required for this package.
"""
from pychem.api import load, parse_smiles, generate_3d, optimize, descriptors, compute_charges, add_hydrogens
from src.core.domain.models.molecule import Molecule
from src.core.domain.models.atom import Atom
from src.core.domain.models.bond import Bond
from src.features.descriptor_calculator.pydes.api import pydes

__version__ = '1.0.0'

__all__ = [
    'Molecule', 'Atom', 'Bond',
    'load', 'parse_smiles', 'generate_3d', 'optimize', 'descriptors',
    'compute_charges', 'add_hydrogens', 'pydes'
]
