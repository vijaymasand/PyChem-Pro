# src/services/rendering/renderer_factory.py
"""
Renderer factory -- picks QPainter or OpenGL based on molecule size.

Threshold: 500 atoms
- Below: QPainter (high quality gradients, no GPU needed)
- Above: OpenGL (hardware accelerated, handles 10K+ atoms)
- Fallback: if OpenGL fails, always use QPainter

Cross-platform: works on Windows, macOS, Linux.
macOS: Apple deprecated OpenGL but 3.3 still works via Core Profile.
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)

ATOM_THRESHOLD = 500

# Protein molecules should always use OpenGL for better performance
PROTEIN_ATOM_THRESHOLD = 100  # Lower threshold for proteins


class RendererFactory:
    """Decides which renderer to use for a given molecule.

    Usage::

        factory = RendererFactory()

        # During set_molecule:
        if factory.should_use_gl(molecule):
            if factory.check_gl_available(gl_widget):
                # switch to GL widget
            else:
                # stay on QPainter (GL not available)
        else:
            # use QPainter (small molecule)
    """

    def __init__(self):
        self._gl_available: Optional[bool] = None  # Lazy check; None = not yet tested

    # ------------------------------------------------------------------
    #  Decision helpers
    # ------------------------------------------------------------------

    def should_use_gl(self, mol) -> bool:
        """Return True if the molecule is large enough to benefit from GL.

        Parameters
        ----------
        mol : Molecule | None
            The molecule to check.  Must have a ``num_atoms`` property
            or an ``atoms`` list.

        Returns
        -------
        bool
        """
        if mol is None:
            return False

        # Support both .num_atoms property and len(.atoms)
        n = getattr(mol, 'num_atoms', None)
        if n is None:
            atoms = getattr(mol, 'atoms', None)
            n = len(atoms) if atoms is not None else 0

        # Check if this is a protein molecule
        is_protein = False
        if hasattr(mol, 'properties'):
            is_protein = mol.properties.get('is_protein', False)

        # Use lower threshold for proteins to ensure OpenGL acceleration
        threshold = PROTEIN_ATOM_THRESHOLD if is_protein else ATOM_THRESHOLD

        print(f"[DEBUG] should_use_gl: n={n}, is_protein={is_protein}, threshold={threshold}, result={n >= threshold}")
        return n >= threshold

    def check_gl_available(self, gl_widget) -> bool:
        """Check whether OpenGL is actually usable.

        Call this after the GL widget has been shown at least once so that
        ``initializeGL`` has had a chance to run.

        Parameters
        ----------
        gl_widget : GLMoleculeWidget | None
            The GL widget instance (may be None if creation failed).

        Returns
        -------
        bool
        """
        if self._gl_available is not None:
            return self._gl_available

        if gl_widget is not None and hasattr(gl_widget, 'gl_available'):
            # If gl_available is None, it hasn't initialized yet, but we should assume it's available
            # so that it gets shown and initializeGL() runs.
            self._gl_available = gl_widget.gl_available is not False
            print(f"[DEBUG] check_gl_available: gl_widget.gl_available={gl_widget.gl_available}")
        else:
            self._gl_available = False
            print(f"[DEBUG] check_gl_available: gl_widget is None or has no gl_available attribute")

        if self._gl_available:
            logger.info("OpenGL renderer available -- will use for large molecules")
        else:
            logger.info("OpenGL renderer NOT available -- QPainter fallback for all sizes")

        return self._gl_available

    def reset(self):
        """Reset the cached GL availability flag (e.g. after context recreation)."""
        self._gl_available = None

    @property
    def gl_available(self) -> Optional[bool]:
        """Current GL availability state (None = not yet tested)."""
        return self._gl_available
