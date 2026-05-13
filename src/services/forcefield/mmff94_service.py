# src/services/forcefield/mmff94_service.py
"""
MMFF94 force field service — thin adapter over the vectorized
_engine package (Phase 1 of the Jmol port).

Public surface unchanged from the legacy file:
    - MMFF94Service class
    - add_hydrogens, assign_atom_types, assign_charges, compute_energy,
      optimize_geometry, optimize_geometry_batch

Internals delegate to:
    - HandCodedAtomTyper (atom_typing/) — sets atom.mmff_type, mmff_class
    - assign_bci_charges (_engine/charges.py) — sets atom.partial_charge
    - ArraysBuilder (_engine/arrays.py) — builds NumPy interaction arrays
    - MMFF94Engine (_engine/engine.py) — computes energy + gradient
    - SteepestDescent / LBFGS (_engine/optimizers.py) — runs the loop

See docs/superpowers/specs/2026-05-11-mmff94-jmol-port-design.md.
"""
from __future__ import annotations
import numpy as np

from src.core.protocols.forcefield import OptimizationResult
from src.core.domain.models.molecule import Molecule
from src.services.forcefield.hydrogen import HydrogenAdder
from src.services.forcefield.atom_typing import AtomTyper, default_typer
from src.services.forcefield._engine.charges import assign_bci_charges
from src.services.forcefield._engine.arrays import ArraysBuilder
from src.services.forcefield._engine.engine import MMFF94Engine
from src.services.forcefield._engine.optimizers import SteepestDescent, LBFGS


# ─── Module-level worker for batch multiprocessing ─────────────────

def _optimize_worker(args):
    """Picklable worker for optimize_geometry_batch.

    Lives at module scope so multiprocessing.spawn can pickle it.
    """
    mol, max_iters, convergence, method = args
    service = MMFF94Service(executor=None)  # No nested pools
    result = service.optimize_geometry(mol, max_iters, convergence, method)
    return mol, result


class MMFF94Service:
    """Vectorized MMFF94 force field service.

    Phases 1-2: hand-coded ~43-type atom typing, all 7 term calculators
    with corrected formulas, vectorized NumPy. Phase 3 will swap
    HandCodedAtomTyper for SmartsAtomTyper (full SMARTS) behind the
    same AtomTyper interface — engine and parameters unchanged.
    """

    def __init__(self, executor=None, atom_typer: AtomTyper | None = None):
        self.executor = executor
        self.atom_typer = atom_typer if atom_typer is not None else default_typer()
        self.h_adder = HydrogenAdder()

    # ─── IForceField protocol ────────────────────────────────────

    def add_hydrogens(self, mol: Molecule) -> int:
        return self.h_adder.add_hydrogens(mol)

    def assign_atom_types(self, mol: Molecule) -> None:
        self.atom_typer.type_atoms(mol)

    def assign_charges(self, mol: Molecule) -> None:
        assign_bci_charges(mol)

    def compute_energy(self, mol: Molecule) -> float:
        self._setup(mol)
        coords = self._coords(mol)
        engine = MMFF94Engine(ArraysBuilder.build(mol))
        return float(engine.energy(coords))

    def optimize_geometry(self, mol: Molecule, max_iters: int = 500,
                          convergence: float = 1e-4,
                          method: str = "lbfgs") -> OptimizationResult:
        # 0. Seed 3D coordinates from 2D/SMILES if missing (existing logic).
        self._seed_coords_from_2d(mol)

        # 1. Set up: hybridization → aromaticity → H's → types → charges.
        self._setup(mol)

        # 2. Build interaction arrays for this molecule.
        arrays = ArraysBuilder.build(mol)
        engine = MMFF94Engine(arrays)

        # 3. Extract coordinates as a (N,3) float64 NumPy array.
        coords = self._coords(mol)

        # 4. Pick optimizer.
        opt_cls = {"lbfgs": LBFGS, "steepest_descent": SteepestDescent}.get(method, LBFGS)
        coords_out, traj, converged, steps = opt_cls(engine.energy_and_gradient).run(
            coords, max_iters=max_iters, convergence=convergence
        )

        # 5. Write coords back to atoms.
        for i, atom in enumerate(mol.atoms):
            atom.x = float(coords_out[i, 0])
            atom.y = float(coords_out[i, 1])
            atom.z = float(coords_out[i, 2])

        # 6. Final RMS gradient.
        final_e, final_g = engine.energy_and_gradient(coords_out)
        rms_grad = float(np.sqrt(np.mean(final_g * final_g)))

        return OptimizationResult(
            converged=converged,
            final_energy=float(traj[-1] if traj else final_e),
            num_steps=steps,
            energy_trajectory=[float(x) for x in traj],
            rms_gradient=rms_grad,
        )

    def optimize_geometry_batch(self, molecules, max_iters: int = 500,
                                convergence: float = 1e-4,
                                method: str = "lbfgs"):
        if not molecules:
            return []
        if self.executor is None or len(molecules) == 1:
            return [self.optimize_geometry(m, max_iters, convergence, method)
                    for m in molecules]

        args = [(m, max_iters, convergence, method) for m in molecules]
        out = self.executor.map(_optimize_worker, args)
        results = []
        for orig, (opt_mol, res) in zip(molecules, out):
            for j, a in enumerate(opt_mol.atoms):
                orig.atoms[j].x = a.x
                orig.atoms[j].y = a.y
                orig.atoms[j].z = a.z
            results.append(res)
        return results

    # ─── Internals ───────────────────────────────────────────────

    def _setup(self, mol: Molecule) -> None:
        mol.assign_hybridization()
        mol.perceive_aromaticity()
        self.h_adder.add_hydrogens(mol)
        self.atom_typer.type_atoms(mol)
        assign_bci_charges(mol)

    def _coords(self, mol: Molecule) -> np.ndarray:
        c = np.array([[a.x, a.y, a.z] for a in mol.atoms], dtype=np.float64)
        if np.allclose(c[:, 2], 0.0, atol=0.01):
            # Break perfect Z=0 planarity to avoid getting stuck in the
            # 2D-seeded local minimum.
            c[:, 2] += np.random.RandomState(31).uniform(-0.1, 0.1, len(c))
        return c

    def _seed_coords_from_2d(self, mol: Molecule) -> None:
        """Ensure every atom has 3D coords before optimization.

        Priority: keep existing 3D → else use OASA 2D (x2d/y2d) → else
        generate 2D layout now → last resort, origin. Z is seeded with
        a small per-atom random offset to avoid flat starting points.
        This is the unchanged v1 behavior, kept verbatim because it's
        a separate concern from MMFF94 correctness.
        """
        needs_seed = any(
            a.x is None or a.y is None or a.z is None for a in mol.atoms
        )
        if not needs_seed:
            return

        have_2d = all(
            getattr(a, "x2d", None) is not None
            and getattr(a, "y2d", None) is not None
            for a in mol.atoms
        )
        if not have_2d:
            try:
                from src.features.layout_2d.generators.coordgen2d_smiles_pure_oasa import (
                    CoordinateGenerator2DSMILES,
                )
                CoordinateGenerator2DSMILES(mol, force_regenerate=True).generate()
            except Exception:
                try:
                    from src.features.layout_2d.generators.coordgen2d import (
                        CoordinateGenerator2D,
                    )
                    CoordinateGenerator2D(mol, force_regenerate=False).generate()
                except Exception:
                    pass

        mol.assign_hybridization()

        rng = np.random.RandomState(17)
        for a in mol.atoms:
            if a.x is None:
                a.x = float(a.x2d) if getattr(a, "x2d", None) is not None else 0.0
            if a.y is None:
                a.y = float(a.y2d) if getattr(a, "y2d", None) is not None else 0.0
            if a.z is None:
                hyb = getattr(a, "hybridization", "sp3") or "sp3"
                if hyb == "sp":
                    a.z = float(rng.uniform(-0.1, 0.1))
                elif hyb == "sp2":
                    a.z = float(rng.uniform(-0.2, 0.2))
                else:
                    a.z = float(rng.uniform(-0.5, 0.5))
