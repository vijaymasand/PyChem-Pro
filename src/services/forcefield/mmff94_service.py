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


def _segment_crosses_triangle(p1, p2, v0, v1, v2) -> bool:
    """Möller–Trumbore intersection of segment p1→p2 with triangle v0,v1,v2."""
    d = p2 - p1
    e1 = v1 - v0
    e2 = v2 - v0
    h = np.cross(d, e2)
    a = float(e1 @ h)
    if -1e-9 < a < 1e-9:
        return False  # segment parallel to the triangle plane
    f = 1.0 / a
    s = p1 - v0
    u = f * float(s @ h)
    if u < 0.0 or u > 1.0:
        return False
    q = np.cross(s, e1)
    v = f * float(d @ q)
    if v < 0.0 or u + v > 1.0:
        return False
    t = f * float(e2 @ q)
    return 1e-6 < t < 1.0 - 1e-6  # strictly inside the segment


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
        """Assign MMFF94 BCI partial charges in place.

        MMFF94 BCI charges are defined on the *hydrogen-complete, MMFF-typed*
        molecule, so the molecule must be fully prepared first — otherwise a
        freshly parsed/loaded structure (no explicit H, ``mmff_type == 0`` on
        every atom) has every bond skipped and comes out with all-zero charges.

        This runs the same preparation as :meth:`optimize_geometry`: seed 3-D
        coordinates if missing (needed only so hydrogens can be *placed*; BCI
        charges themselves are geometry-independent), then add explicit
        hydrogens, perceive aromaticity, and assign atom types before the BCI
        pass.  ``add_hydrogens`` is idempotent, so calling this after an
        optimization simply re-derives the (identical) charges without adding
        duplicate atoms; and for an already-loaded 3-D structure the coordinate
        seeding is a no-op that leaves the geometry untouched.
        """
        self._seed_coords_from_2d(mol)
        self._setup(mol)

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

        # 1. Set up: hybridization → aromaticity → add H → types → charges.
        #    Hydrogens are added FIRST (a bare heavy-atom skeleton has the
        #    wrong valence/geometry) and charges are assigned here, BEFORE the
        #    minimisation, because the MMFF94 energy includes an electrostatic
        #    term (es_qq = q_i·q_j) that the optimizer needs.  BCI charges are
        #    geometry-independent, so assigning them on the seeded start
        #    geometry gives exactly the charges of the final optimised one.
        self._setup(mol)

        # 2. Build interaction arrays for this molecule.
        arrays = ArraysBuilder.build(mol)
        engine = MMFF94Engine(arrays)

        # 3. Extract coordinates as a (N,3) float64 NumPy array.
        coords = self._coords(mol)

        # 4. Optimize.  A near-flat 2-D-seeded start can trap the local
        #    optimizer in a *tangled* minimum — e.g. a substituent chain
        #    threaded through a ring, seen as a bond passing through the ring
        #    face.  Such a geometry carries a non-bonded clash and a high
        #    (strained) energy, so we detect it and, only then, retry from
        #    randomised 3-D starts and keep the lowest-energy untangled result.
        #    Well-behaved molecules pass on the first try and pay nothing extra.
        opt_cls = {"lbfgs": LBFGS, "steepest_descent": SteepestDescent}.get(method, LBFGS)

        def _run(c0):
            return opt_cls(engine.energy_and_gradient).run(
                c0, max_iters=max_iters, convergence=convergence)

        coords_out, traj, converged, steps = _run(coords)
        # Only a *tangled* result warrants the extra restarts — a clash-free
        # geometry that merely hit the iteration cap (common for floppy
        # molecules relaxing soft modes) is fine and pays nothing.
        if self._bad_contacts(mol, coords_out) > 0:
            coords_out, traj, converged, steps = self._reoptimize_untangled(
                mol, coords, _run, (coords_out, traj, converged, steps))

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

    # ─── Tangle detection & recovery ─────────────────────────────

    _CLASH_DIST = 2.0   # non-bonded heavy-atom separation (Å) below which a
                        # contact is a genuine clash (a threaded/overlapping
                        # geometry), not a normal van-der-Waals contact.

    def _excluded_pairs(self, mol):
        """1-2 (bonded) and 1-3 (geminal) atom pairs — legitimately close."""
        excl = set()
        for b in mol.bonds:
            i, j = b.begin_atom_idx, b.end_atom_idx
            excl.add((i, j) if i < j else (j, i))
        for i in range(len(mol.atoms)):
            neigh = mol.get_neighbors(i)
            for a in range(len(neigh)):
                for b in range(a + 1, len(neigh)):
                    x, y = neigh[a], neigh[b]
                    excl.add((x, y) if x < y else (y, x))
        return excl

    def _bad_contacts(self, mol, coords) -> int:
        """Number of non-bonded heavy-atom clashes plus bonds threading a ring.

        A threaded ring (chain passing through the ring face) always produces
        at least one such close contact, so this is a reliable, cheap signal
        for "the optimizer landed in a tangled minimum".
        """
        sym = [a.symbol for a in mol.atoms]
        heavy = [i for i in range(len(sym)) if sym[i] != 'H']
        excl = self._excluded_pairs(mol)
        thr2 = self._CLASH_DIST * self._CLASH_DIST
        clashes = 0
        for x in range(len(heavy)):
            i = heavy[x]
            ci = coords[i]
            for y in range(x + 1, len(heavy)):
                j = heavy[y]
                if (i, j) in excl:
                    continue
                d = ci - coords[j]
                if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < thr2:
                    clashes += 1
        return clashes + self._ring_threads(mol, coords)

    @staticmethod
    def _ring_threads(mol, coords) -> int:
        """Count bonds whose segment passes through a ring face.

        Each ring is triangulated as a fan from its centroid and every
        non-ring bond is tested against those triangles with a Möller–Trumbore
        segment/triangle intersection.  This catches a chain threaded through a
        ring at *any* angle (unlike a plane-crossing test, which misses a bond
        lying nearly in the ring plane) and does not false-trigger on ordinary
        substituents that merely sit next to the ring.
        """
        rings = mol.find_rings()
        if not rings:
            return 0
        bonds = [(b.begin_atom_idx, b.end_atom_idx) for b in mol.bonds]
        threads = 0
        for ring in rings:
            idx = list(ring)
            pts = coords[idx]
            c = pts.mean(axis=0)
            rs = set(idx)
            tris = [(c, pts[k], pts[(k + 1) % len(idx)]) for k in range(len(idx))]
            for i, j in bonds:
                if i in rs or j in rs:
                    continue  # ring bonds and ring-attached bonds are legitimate
                p1, p2 = coords[i], coords[j]
                if any(_segment_crosses_triangle(p1, p2, v0, v1, v2)
                       for (v0, v1, v2) in tris):
                    threads += 1
        return threads

    def _reoptimize_untangled(self, mol, seed_coords, run_fn, best):
        """Retry the minimisation from randomised 3-D starts, keeping the
        lowest-energy, least-tangled result.

        The kicks use fixed RNG seeds so the outcome is fully deterministic.
        Selection is lexicographic — fewest bad contacts first, then lowest
        energy — so an untangled geometry always beats a tangled one even if
        the tangled one happens to score a slightly lower raw energy.
        """
        b_coords, b_traj, b_conv, b_steps = best
        b_bad = self._bad_contacts(mol, b_coords)
        b_e = b_traj[-1] if b_traj else float("inf")

        for k in range(6):
            if b_bad == 0:
                break  # untangled — stop, even if soft modes are still relaxing
            sigma = 0.8 + 0.4 * k  # escalate the 3-D spread each attempt
            kick = np.random.RandomState(1009 + k).normal(0.0, sigma, seed_coords.shape)
            co, tr, cv, st = run_fn(seed_coords + kick)
            bad = self._bad_contacts(mol, co)
            e = tr[-1] if tr else float("inf")
            if (bad, e) < (b_bad, b_e):
                b_coords, b_traj, b_conv, b_steps = co, tr, cv, st
                b_bad, b_e = bad, e
        return b_coords, b_traj, b_conv, b_steps

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
