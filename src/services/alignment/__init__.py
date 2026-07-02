"""
Structural alignment service (pure NumPy).

Superimpose molecules with the Kabsch algorithm and iterative outlier
rejection; build atom correspondences by index, element, or protein Cα
sequence.  Ported from the patinae-algos (Rust) ``align`` crate, adapted to
PyChem's :class:`Molecule` model with no new dependencies.

Public surface::

    from src.services.alignment import AlignmentService, AlignmentResult

or, more conveniently, via the ``pychem`` facade::

    import pychem
    pychem.align(mobile, reference)
    pychem.rmsd(mol_a, mol_b)
"""
from src.services.alignment.align_service import AlignmentService, AlignmentResult

__all__ = ["AlignmentService", "AlignmentResult"]
