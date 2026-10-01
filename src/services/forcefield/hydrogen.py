"""
HydrogenAdder — Adds explicit hydrogens with 3D coordinate placement.

Fixes MMFF94 force field geometry by ensuring all implicit hydrogens are
converted to explicit atoms with hybridization-aware 3D positions
(tetrahedral for sp3, trigonal planar for sp2, linear for sp).
"""

import numpy as np

from src.core.domain.models.atom import Atom
from src.core.domain.models.bond import BondType


# Ideal X-H bond lengths in Angstroms
_H_BOND_LENGTHS = {
    'C': 1.09,
    'N': 1.01,
    'O': 0.96,
    'S': 1.34,
}
_DEFAULT_H_BOND_LENGTH = 1.00

# Standard valences used for PDB inference (first typical valence)
_STANDARD_VALENCE = {
    'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1,
    'P': 3, 'S': 2, 'Cl': 1, 'Br': 1, 'I': 1,
    'B': 3, 'Si': 4,
}

# Tetrahedral reference directions (unit vectors pointing to vertices)
_TETRAHEDRAL_DIRS = np.array([
    [1.0, 1.0, 1.0],
    [1.0, -1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
]) / np.sqrt(3.0)

# Trigonal planar reference directions (in the XY plane)
_TRIGONAL_DIRS = np.array([
    [1.0, 0.0, 0.0],
    [-0.5, np.sqrt(3.0) / 2.0, 0.0],
    [-0.5, -np.sqrt(3.0) / 2.0, 0.0],
])

# Linear reference directions
_LINEAR_DIRS = np.array([
    [1.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0],
])


class HydrogenAdder:
    """
    Adds explicit hydrogen atoms to a molecular graph with 3D placement.

    Pipeline:
        1. Assign hybridization to all atoms
        2. Infer missing hydrogen counts (for PDB-loaded molecules)
        3. Add explicit H atoms and bonds to the molecular graph
        4. Place H atoms in 3D using hybridization-aware geometry
    """

    def add_hydrogens(self, mol):
        """
        Add explicit hydrogens to the molecule with 3D coordinates.

        Args:
            mol: Molecule instance with heavy atoms already placed in 3D.

        Returns:
            int: Number of hydrogen atoms added.
        """
        # Step 1: Assign hybridization so we know geometry
        mol.assign_hybridization()

        # Step 2: Infer missing H counts for PDB atoms that have no H info
        self._infer_missing_h_counts(mol)

        # Step 3: Record which atoms need H, then add explicit H atoms
        heavy_count = mol.num_atoms
        count = mol.add_explicit_hydrogens()

        # Step 4: Place each new H atom in 3D
        self._place_hydrogens_3d(mol, heavy_count)

        return count

    def _infer_missing_h_counts(self, mol):
        """
        For atoms with no hydrogen information (typical of PDB loading),
        infer the number of missing hydrogens from standard valence minus
        the sum of bond orders.
        """
        for atom in mol.atoms:
            if atom.symbol == 'H':
                continue

            # If the atom already has implicit or explicit H info, skip
            if atom.num_implicit_h > 0:
                continue
            if atom.num_explicit_h is not None and atom.num_explicit_h > 0:
                continue

            # Get standard valence for this element
            std_valence = _STANDARD_VALENCE.get(atom.symbol)
            if std_valence is None:
                # Try from element's typical_valences
                if hasattr(atom, 'element') and atom.element.typical_valences:
                    std_valence = atom.element.typical_valences[0]
                else:
                    continue

            # Sum of bond orders to existing neighbors
            neighbor_bonds = mol.get_neighbor_bonds(atom.index)
            bond_order_sum = sum(b.order for _, b in neighbor_bonds)

            # Account for formal charge
            missing = int(std_valence - bond_order_sum - atom.formal_charge)

            # Only infer if the atom already has a mostly-complete heavy-atom
            # skeleton (avoids over-hydrogenating isolated terminal atoms that
            # were not part of a PDB-style backbone).  Heuristic: the number
            # of missing H should not exceed the number of existing bonds.
            n_bonds = len(neighbor_bonds)
            if missing > 0 and missing <= n_bonds:
                atom.num_implicit_h = missing

    def _place_hydrogens_3d(self, mol, heavy_count):
        """
        Place newly added hydrogen atoms in 3D space.

        Args:
            mol: Molecule with H atoms already added (indices >= heavy_count).
            heavy_count: Number of atoms before H addition (all heavy atoms).
        """
        for h_idx in range(heavy_count, mol.num_atoms):
            h_atom = mol.atoms[h_idx]
            if h_atom.symbol != 'H':
                continue

            # Find the parent atom (the heavy atom bonded to this H)
            neighbors = mol.get_neighbors(h_idx)
            if not neighbors:
                continue
            parent_idx = neighbors[0]
            parent = mol.atoms[parent_idx]

            self._place_single_h(mol, h_atom, parent, heavy_count)

    def _place_single_h(self, mol, h_atom, parent, heavy_count):
        """
        Place a single hydrogen atom relative to its parent using
        hybridization-aware geometry.
        """
        # Determine bond length
        bond_length = _H_BOND_LENGTHS.get(parent.symbol, _DEFAULT_H_BOND_LENGTH)

        # Gather existing neighbor directions (from parent to each neighbor,
        # excluding the current H atom being placed)
        existing_dirs = []
        parent_pos = np.array([parent.x, parent.y, parent.z])

        for n_idx in mol.get_neighbors(parent.index):
            if n_idx == h_atom.index:
                continue
            n_atom = mol.atoms[n_idx]
            if not n_atom.has_coords:
                continue
            n_pos = np.array([n_atom.x, n_atom.y, n_atom.z])
            diff = n_pos - parent_pos
            norm = np.linalg.norm(diff)
            if norm > 1e-8:
                existing_dirs.append(diff / norm)

        # Choose direction for this H
        hybridization = parent.hybridization or 'sp3'
        direction = self._ideal_h_direction(hybridization, existing_dirs)

        # Place the H atom
        h_pos = parent_pos + direction * bond_length
        h_atom.x, h_atom.y, h_atom.z = float(h_pos[0]), float(h_pos[1]), float(h_pos[2])

    def _ideal_h_direction(self, hybridization, existing_dirs):
        """
        Choose an ideal direction for a new hydrogen bond given the
        hybridization type and directions to existing neighbors.

        Returns a unit vector pointing from parent toward the new H position.
        """
        n_existing = len(existing_dirs)

        if hybridization == 'sp3':
            return self._pick_from_template(_TETRAHEDRAL_DIRS, existing_dirs)
        elif hybridization == 'sp2':
            return self._pick_from_template(_TRIGONAL_DIRS, existing_dirs)
        elif hybridization == 'sp':
            return self._pick_from_template(_LINEAR_DIRS, existing_dirs)
        else:
            # Default: sp3
            return self._pick_from_template(_TETRAHEDRAL_DIRS, existing_dirs)

    def _pick_from_template(self, template_dirs, existing_dirs):
        """
        Given a set of template directions (e.g. tetrahedral vertices) and
        a list of existing neighbor directions, pick the template direction
        that is most "away" from all existing neighbors.

        If there are no existing neighbors, return the first template direction.
        If existing neighbors exist, rotate the template to best match existing
        neighbors, then return the best unused slot.
        """
        if not existing_dirs:
            return template_dirs[0].copy()

        n_existing = len(existing_dirs)
        n_template = len(template_dirs)

        if n_existing >= n_template:
            # All template slots used; place opposite to the mean of existing
            mean_dir = np.mean(existing_dirs, axis=0)
            norm = np.linalg.norm(mean_dir)
            if norm > 1e-8:
                return -mean_dir / norm
            else:
                return self._find_perpendicular(existing_dirs[0])

        # Strategy: find the direction that maximizes minimum angle from
        # all existing neighbors. We do this by:
        # 1. Computing the "anti-mean" direction (opposite to average of existing)
        # 2. If there's only one existing dir, pick a nice tetrahedral/trigonal angle

        if n_existing == 1:
            # One existing neighbor: place at the ideal angle
            existing = np.array(existing_dirs[0])
            if n_template == 4:
                # Tetrahedral: ~109.5 degrees from existing
                return self._place_at_angle(existing, np.arccos(-1.0 / 3.0))
            elif n_template == 3:
                # Trigonal: 120 degrees from existing
                return self._place_at_angle(existing, 2.0 * np.pi / 3.0)
            elif n_template == 2:
                # Linear: 180 degrees
                return -existing
            else:
                return self._place_at_angle(existing, np.arccos(-1.0 / 3.0))

        elif n_existing == 2:
            # Two existing neighbors
            d1 = np.array(existing_dirs[0])
            d2 = np.array(existing_dirs[1])

            if n_template == 4:
                # Tetrahedral: place in the plane bisecting d1, d2 but opposite
                bisector = d1 + d2
                bisector_norm = np.linalg.norm(bisector)
                if bisector_norm > 1e-8:
                    bisector = bisector / bisector_norm
                    # Cross product gives the perpendicular to the d1-d2 plane
                    cross = np.cross(d1, d2)
                    cross_norm = np.linalg.norm(cross)
                    if cross_norm > 1e-8:
                        cross = cross / cross_norm
                        # Tetrahedral angle from bisector
                        anti = -bisector
                        # Tilt away from bisector plane
                        angle = np.arccos(-1.0 / 3.0) / 2.0
                        return (anti * np.cos(angle) + cross * np.sin(angle))
                    else:
                        return -bisector
                else:
                    # d1 and d2 are opposite: pick perpendicular
                    return self._find_perpendicular(d1)

            elif n_template == 3:
                # Trigonal planar: the third direction is opposite the bisector
                bisector = d1 + d2
                bisector_norm = np.linalg.norm(bisector)
                if bisector_norm > 1e-8:
                    return -bisector / bisector_norm
                else:
                    return self._find_perpendicular(d1)

            elif n_template == 2:
                # Linear: shouldn't reach here with 2 existing, but fallback
                mean = d1 + d2
                norm = np.linalg.norm(mean)
                if norm > 1e-8:
                    return -mean / norm
                return self._find_perpendicular(d1)

        elif n_existing == 3:
            # Three existing neighbors (only for sp3 with one H left)
            mean_dir = np.mean(existing_dirs, axis=0)
            norm = np.linalg.norm(mean_dir)
            if norm > 1e-8:
                return -mean_dir / norm
            else:
                return self._find_perpendicular(existing_dirs[0])

        # Fallback
        mean_dir = np.mean(existing_dirs, axis=0)
        norm = np.linalg.norm(mean_dir)
        if norm > 1e-8:
            return -mean_dir / norm
        return np.array([0.0, 0.0, 1.0])

    def _place_at_angle(self, existing_dir, angle):
        """
        Place a new direction at a given angle from an existing direction.
        Rotates around a perpendicular axis.
        """
        perp = self._find_perpendicular(existing_dir)
        # Rodrigues' rotation: rotate existing_dir around perp by angle
        d = np.array(existing_dir)
        result = (d * np.cos(angle) +
                  np.cross(perp, d) * np.sin(angle) +
                  perp * np.dot(perp, d) * (1.0 - np.cos(angle)))
        norm = np.linalg.norm(result)
        if norm > 1e-8:
            return result / norm
        return perp

    def _find_perpendicular(self, v):
        """
        Find a unit vector perpendicular to v using Gram-Schmidt.
        """
        v = np.array(v, dtype=float)
        norm = np.linalg.norm(v)
        if norm < 1e-8:
            return np.array([1.0, 0.0, 0.0])
        v = v / norm

        # Choose a candidate that's not parallel to v
        if abs(v[0]) < 0.9:
            candidate = np.array([1.0, 0.0, 0.0])
        else:
            candidate = np.array([0.0, 1.0, 0.0])

        # Gram-Schmidt: remove component along v
        perp = candidate - np.dot(candidate, v) * v
        perp_norm = np.linalg.norm(perp)
        if perp_norm < 1e-8:
            return np.array([0.0, 0.0, 1.0])
        return perp / perp_norm
