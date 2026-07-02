"""
Alignment service — superimpose one molecule onto another and report RMSD.

This is the public, GUI-agnostic entry point for structural alignment.  It
chooses an atom-correspondence strategy (see :mod:`correspondence`), fits an
optimal rigid-body transform with iterative outlier rejection (see
:mod:`_engine.superpose`), and — by default — applies that transform *in place*
to every atom of the mobile molecule so overlays and exports pick it up.

Typical use (Jupyter / scripts)::

    import pychem
    ref = pychem.load("ref.pdb")
    mob = pychem.load("model2.pdb")
    result = pychem.align(mob, ref)          # moves `mob` onto `ref`
    print(result.rmsd, result.n_aligned)

The Rust reference (patinae-algos ``align``) additionally offers Combinatorial
Extension (CE) for sequence-independent fold matching; that is intentionally not
ported here — the ``sequence`` method covers same/homologous proteins, which is
the case that arises when several structures are opened together.
"""
from __future__ import annotations

import numpy as np

from src.services.alignment import correspondence as corr
from src.services.alignment._engine.kabsch import kabsch
from src.services.alignment._engine.superpose import superpose


class AlignmentResult:
    """Outcome of an :meth:`AlignmentService.align` call.

    Attributes:
        rmsd:           RMSD (Å) over the aligned atom pairs after fitting.
        initial_rmsd:   RMSD before outlier rejection.
        n_aligned:      Atom pairs used in the final fit.
        n_rejected:     Pairs discarded as outliers.
        rotation:       (3, 3) rotation applied to the mobile molecule.
        translation:    (3,) translation applied after rotation.
        method:         Correspondence strategy actually used.
        selection:      Atom selection actually used.
        transformed:    Whether the mobile molecule's coordinates were moved.
        pairs:          Final ``(mobile_idx, reference_idx)`` atom pairs.
    """

    __slots__ = ("rmsd", "initial_rmsd", "n_aligned", "n_rejected",
                 "rotation", "translation", "method", "selection",
                 "transformed", "pairs")

    def __init__(self, rmsd, initial_rmsd, n_aligned, n_rejected, rotation,
                 translation, method, selection, transformed, pairs):
        self.rmsd = rmsd
        self.initial_rmsd = initial_rmsd
        self.n_aligned = n_aligned
        self.n_rejected = n_rejected
        self.rotation = rotation
        self.translation = translation
        self.method = method
        self.selection = selection
        self.transformed = transformed
        self.pairs = pairs

    def __repr__(self):
        return (f"AlignmentResult(rmsd={self.rmsd:.4f} Å, n_aligned={self.n_aligned}, "
                f"n_rejected={self.n_rejected}, method='{self.method}', "
                f"selection='{self.selection}')")


class AlignmentService:
    """Superimpose molecules via Kabsch fitting with outlier rejection.

    Args:
        loader: Optional object with a ``load(path)`` method, used when
                :meth:`align` / :meth:`rmsd` are given file paths instead of
                :class:`Molecule` instances.  Falls back to the default
                :class:`LoaderService` when omitted.
    """

    def __init__(self, loader=None):
        self._loader = loader

    # ── Public API ───────────────────────────────────────────────────

    def align(self, mobile, reference, method='auto', selection='auto',
              cycles=5, cutoff=2.0, weights=None, transform=True) -> AlignmentResult:
        """Superimpose *mobile* onto *reference*.

        Args:
            mobile:    Molecule (or file path) to be moved.
            reference: Molecule (or file path) held fixed.
            method:    ``'auto'`` | ``'index'`` | ``'element'`` | ``'sequence'``.
                       ``'auto'`` picks ``'sequence'`` for two proteins, else
                       ``'index'``.
            selection: ``'auto'`` | ``'all'`` | ``'heavy'`` | ``'ca'`` |
                       ``'backbone'``.  ``'auto'`` uses ``'ca'`` with the
                       sequence method, otherwise ``'heavy'``.
            cycles:    Outlier-rejection cycles (0 → plain single Kabsch).
            cutoff:    Reject a pair when ``distance / RMSD > cutoff``.
            weights:   Optional per-pair weights; forces a single weighted
                       Kabsch (no iterative rejection).
            transform: If True, move every atom of *mobile* in place.

        Returns:
            :class:`AlignmentResult`.

        Raises:
            ValueError: If fewer than 3 atom correspondences can be built.
        """
        mob = self._as_molecule(mobile)
        ref = self._as_molecule(reference)

        method, selection = self._resolve(mob, ref, method, selection)
        pairs = self._build_pairs(mob, ref, method, selection)
        if len(pairs) < 3:
            raise ValueError(
                f"only {len(pairs)} atom correspondence(s) found via "
                f"method='{method}', selection='{selection}' — need at least 3. "
                f"Try a different method/selection or check the structures share atoms.")

        src_all = self._all_coords(mob)
        tgt_all = self._all_coords(ref)

        if weights is not None or cycles <= 0:
            # Single (optionally weighted) Kabsch — no outlier rejection.
            si = [p[0] for p in pairs]
            ti = [p[1] for p in pairs]
            kres = kabsch(src_all[si], tgt_all[ti], weights=weights)
            rotation, translation = kres.rotation, kres.translation
            rmsd_val = initial = kres.rmsd
            n_aligned, n_rejected = kres.n_atoms, 0
            final_pairs = pairs
        else:
            sres = superpose(src_all, tgt_all, pairs, cycles=cycles, cutoff=cutoff)
            rotation = sres.transform.rotation
            translation = sres.transform.translation
            rmsd_val = sres.final_rmsd
            initial = sres.initial_rmsd
            n_aligned, n_rejected = sres.n_aligned, sres.n_rejected
            final_pairs = sres.pairs

        if transform:
            self._apply_in_place(mob, rotation, translation)

        return AlignmentResult(
            rmsd=rmsd_val, initial_rmsd=initial, n_aligned=n_aligned,
            n_rejected=n_rejected, rotation=rotation, translation=translation,
            method=method, selection=selection, transformed=bool(transform),
            pairs=final_pairs)

    def align_many(self, mobiles, reference, **kwargs):
        """Align each molecule in *mobiles* onto *reference*.

        Returns a list of :class:`AlignmentResult` in the same order.
        """
        ref = self._as_molecule(reference)
        return [self.align(m, ref, **kwargs) for m in mobiles]

    def rmsd(self, mol_a, mol_b, method='auto', selection='auto',
             superpose=False, cycles=5, cutoff=2.0) -> float:
        """RMSD between two molecules over their atom correspondence.

        Args:
            superpose: If False (default), RMSD is computed on the current
                       coordinates *as they are* (no fitting).  If True, an
                       optimal fit is computed first (molecules are not moved).
        """
        mob = self._as_molecule(mol_a)
        ref = self._as_molecule(mol_b)
        method, selection = self._resolve(mob, ref, method, selection)
        pairs = self._build_pairs(mob, ref, method, selection)
        if len(pairs) < 1:
            raise ValueError("no atom correspondences found for RMSD")

        src_all = self._all_coords(mob)
        tgt_all = self._all_coords(ref)
        si = [p[0] for p in pairs]
        ti = [p[1] for p in pairs]
        src, tgt = src_all[si], tgt_all[ti]

        if superpose:
            if len(pairs) < 3:
                raise ValueError("need at least 3 pairs to superpose before RMSD")
            return float(kabsch(src, tgt).rmsd)
        diff = src - tgt
        return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))

    # ── Internals ────────────────────────────────────────────────────

    def _resolve(self, mob, ref, method, selection):
        """Fill in ``'auto'`` method/selection from molecule content."""
        method = (method or 'auto').lower()
        selection = (selection or 'auto').lower()
        both_protein = corr.is_protein(mob) and corr.is_protein(ref)

        if method == 'auto':
            method = 'sequence' if both_protein else 'index'

        if selection == 'auto':
            selection = 'ca' if method == 'sequence' else 'heavy'
        return method, selection

    def _build_pairs(self, mob, ref, method, selection):
        if method == 'index':
            return corr.pairs_by_index(mob, ref, selection)
        if method == 'element':
            return corr.pairs_by_element(mob, ref, selection)
        if method == 'sequence':
            return corr.pairs_by_sequence(mob, ref)
        raise ValueError(f"unknown alignment method: {method!r}")

    @staticmethod
    def _all_coords(mol):
        """``(N, 3)`` array of all atom coordinates (NaN where unset)."""
        arr = np.full((len(mol.atoms), 3), np.nan, dtype=np.float64)
        for i, a in enumerate(mol.atoms):
            if a.x is not None:
                arr[i, 0] = a.x
                arr[i, 1] = a.y
                arr[i, 2] = a.z
        return arr

    @staticmethod
    def _apply_in_place(mol, rotation, translation):
        """Rotate+translate every atom of *mol* that carries coordinates."""
        idx = [i for i, a in enumerate(mol.atoms) if a.x is not None]
        if not idx:
            return
        coords = np.array([[mol.atoms[i].x, mol.atoms[i].y, mol.atoms[i].z]
                           for i in idx], dtype=np.float64)
        moved = coords @ rotation.T + translation
        for k, i in enumerate(idx):
            atom = mol.atoms[i]
            atom.x = float(moved[k, 0])
            atom.y = float(moved[k, 1])
            atom.z = float(moved[k, 2])

    def _as_molecule(self, obj):
        """Accept a Molecule or a path string; load the latter."""
        if isinstance(obj, str):
            return self._load(obj)
        if obj is None or not hasattr(obj, 'atoms'):
            raise TypeError(f"expected a Molecule or file path, got {type(obj)!r}")
        return obj

    def _load(self, path):
        if self._loader is not None:
            return self._loader.load(path)
        from src.services.loading.loader_service import LoaderService
        return LoaderService().load(path)
