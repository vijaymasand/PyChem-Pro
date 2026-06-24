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
from src.core.math_utils import calculate_distance


# ─── Module-level worker for batch multiprocessing ─────────────────

def _optimize_worker(args):
    """Picklable worker for optimize_geometry_batch."""
    mol, max_iters, convergence, method = args
    service = MMFF94Service(executor=None)
    result = service.optimize_geometry(mol, max_iters, convergence, method)
    return mol, result


class MMFF94Service:
    """Vectorized MMFF94 force field service."""

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
        # 1. We MUST seed coordinates first, otherwise HydrogenAdder crashes 
        # because the 'parent' atoms have NoneType coordinates.
        self._seed_coords_from_2d(mol)
        
        # 2. Now run the setup (aromaticity, types, charges, and H-addition)
        self._setup(mol)
        
        # 3. Proceed with calculation
        coords = self._coords(mol)
        engine = MMFF94Engine(ArraysBuilder.build(mol))
        return float(engine.energy(coords))

    def optimize_geometry(self, mol: Molecule, max_iters: int = 500,
                          convergence: float = 1e-4,
                          method: str = "lbfgs") -> OptimizationResult:
        self._seed_coords_from_2d(mol)
        self._setup(mol)
        arrays = ArraysBuilder.build(mol)
        engine = MMFF94Engine(arrays)
        coords = self._coords(mol)

        opt_cls = {"lbfgs": LBFGS, "steepest_descent": SteepestDescent}.get(method, LBFGS)
        coords_out, traj, converged, steps = opt_cls(engine.energy_and_gradient).run(
            coords, max_iters=max_iters, convergence=convergence
        )

        for i, atom in enumerate(mol.atoms):
            atom.x, atom.y, atom.z = map(float, coords_out[i])

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
        if not molecules: return []
        if self.executor is None or len(molecules) == 1:
            return [self.optimize_geometry(m, max_iters, convergence, method) for m in molecules]

        args = [(m, max_iters, convergence, method) for m in molecules]
        out = self.executor.map(_optimize_worker, args)
        results = []
        for orig, (opt_mol, res) in zip(molecules, out):
            for j, a in enumerate(opt_mol.atoms):
                orig.atoms[j].x, orig.atoms[j].y, orig.atoms[j].z = a.x, a.y, a.z
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
            c[:, 2] += np.random.RandomState(31).uniform(-0.1, 0.1, len(c))
        return c

    def _seed_coords_from_2d(self, mol: Molecule) -> None:
        """Ensure every atom has 3D coords before optimization."""
        needs_seed = any(a.x is None or a.y is None or a.z is None for a in mol.atoms)
        if not needs_seed: return

        have_2d = all(getattr(a, "x2d", None) is not None for a in mol.atoms)
        if not have_2d:
            try:
                from src.features.layout_2d.generators.coordgen2d_smiles_pure_oasa import CoordinateGenerator2DSMILES
                CoordinateGenerator2DSMILES(mol, force_regenerate=True).generate()
            except Exception:
                try:
                    from src.features.layout_2d.generators.coordgen2d import CoordinateGenerator2D
                    CoordinateGenerator2D(mol, force_regenerate=False).generate()
                except Exception: pass

        mol.assign_hybridization()
        for a in mol.atoms:
            if a.x is None: a.x = float(a.x2d) if getattr(a, "x2d", None) is not None else 0.0
            if a.y is None: a.y = float(a.y2d) if getattr(a, "y2d", None) is not None else 0.0
            if a.z is None:
                a.z = self._compute_initial_z(a, mol)

    def _compute_initial_z(self, atom, mol) -> float:
        """Detailed initial z-seeding based on hybridization and atom type."""
        import random
        random.seed(atom.index)
        hyb = getattr(atom, 'hybridization', 'sp3') or 'sp3'
        
        if hyb == 'sp': z_base = random.uniform(-0.1, 0.1)
        elif hyb == 'sp2': z_base = random.uniform(-0.2, 0.2)
        else: z_base = random.uniform(-0.5, 0.5)
        
        scales = {'H': 0.3, 'C': 1.0, 'N': 1.0, 'O': 1.0, 'S': 1.2, 'P': 1.2}
        z_scale = scales.get(atom.symbol, 0.8)
        return z_base * z_scale
            