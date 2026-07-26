"""
Extended geometric (3D) descriptors.

Principal moments of inertia and the shape ratios built on them, size along the
principal axes, planarity, and density/globularity style ratios. Every one of
these needs 3D coordinates; without them they return 0.0 rather than raising,
so a descriptor run on a molecule parsed from SMILES still completes.
"""
import numpy as np

from .geometric import GeometricCalculator


class ExtendedGeometricCalculator(GeometricCalculator):
    """Geometric descriptors, basic set plus shape and size families."""

    # ------------------------------------------------------------ helpers
    def _coords_and_masses(self, molecule, selection):
        self.ensure_coordinates(molecule)
        coords, masses = [], []
        for idx in selection.atom_indices:
            if idx >= len(molecule.atoms):
                continue
            atom = molecule.atoms[idx]
            if not getattr(atom, 'has_coords', False) or atom.x is None:
                continue
            coords.append([atom.x, atom.y, atom.z])
            masses.append(getattr(atom, 'mass', 12.011))
        if not coords:
            return None, None
        return np.array(coords, dtype=float), np.array(masses, dtype=float)

    def _inertia_eigenvalues(self, molecule, selection):
        """Principal moments of inertia in amu*A^2, ascending."""
        coords, masses = self._coords_and_masses(molecule, selection)
        if coords is None or len(coords) < 2:
            return None
        center = np.average(coords, axis=0, weights=masses)
        centered = coords - center
        tensor = np.zeros((3, 3))
        for (x, y, z), m in zip(centered, masses):
            tensor[0, 0] += m * (y * y + z * z)
            tensor[1, 1] += m * (x * x + z * z)
            tensor[2, 2] += m * (x * x + y * y)
            tensor[0, 1] -= m * x * y
            tensor[0, 2] -= m * x * z
            tensor[1, 2] -= m * y * z
        tensor[1, 0], tensor[2, 0], tensor[2, 1] = tensor[0, 1], tensor[0, 2], tensor[1, 2]
        return np.sort(np.linalg.eigvalsh(tensor))

    def _principal_axes_extent(self, molecule, selection):
        """Size of the molecule along its three principal axes."""
        coords, masses = self._coords_and_masses(molecule, selection)
        if coords is None or len(coords) < 2:
            return None
        centered = coords - coords.mean(axis=0)
        _, _, axes = np.linalg.svd(centered, full_matrices=False)
        projections = centered.dot(axes.T)
        return np.sort(projections.max(axis=0) - projections.min(axis=0))[::-1]

    # ------------------------------------------------- moments of inertia
    def calc_inertia_a(self, molecule, selection) -> float:
        moments = self._inertia_eigenvalues(molecule, selection)
        return float(moments[0]) if moments is not None else 0.0

    def calc_inertia_b(self, molecule, selection) -> float:
        moments = self._inertia_eigenvalues(molecule, selection)
        return float(moments[1]) if moments is not None else 0.0

    def calc_inertia_c(self, molecule, selection) -> float:
        moments = self._inertia_eigenvalues(molecule, selection)
        return float(moments[2]) if moments is not None else 0.0

    def calc_npr1(self, molecule, selection) -> float:
        """Normalised principal moment ratio I1/I3, 0.5 for a rod, 1 for a disc."""
        moments = self._inertia_eigenvalues(molecule, selection)
        if moments is None or moments[2] <= 0:
            return 0.0
        return float(moments[0] / moments[2])

    def calc_npr2(self, molecule, selection) -> float:
        """Normalised principal moment ratio I2/I3, 1 for both rods and discs."""
        moments = self._inertia_eigenvalues(molecule, selection)
        if moments is None or moments[2] <= 0:
            return 0.0
        return float(moments[1] / moments[2])

    def calc_inertial_shape_factor(self, molecule, selection) -> float:
        """I2 / (I1 * I3), large for elongated molecules."""
        moments = self._inertia_eigenvalues(molecule, selection)
        if moments is None or moments[0] <= 0 or moments[2] <= 0:
            return 0.0
        return float(moments[1] / (moments[0] * moments[2]))

    def calc_spherocity_index(self, molecule, selection) -> float:
        """3 * smallest gyration eigenvalue / trace, 1 for a sphere."""
        eigenvalues = np.sort(self.get_gyration_eigenvalues(molecule, selection))
        total = float(np.sum(eigenvalues))
        if total <= 0:
            return 0.0
        return float(3.0 * eigenvalues[0] / total)

    # ------------------------------------------------------------- extent
    def calc_molecular_length(self, molecule, selection) -> float:
        extent = self._principal_axes_extent(molecule, selection)
        return float(extent[0]) if extent is not None else 0.0

    def calc_molecular_width(self, molecule, selection) -> float:
        extent = self._principal_axes_extent(molecule, selection)
        return float(extent[1]) if extent is not None else 0.0

    def calc_molecular_thickness(self, molecule, selection) -> float:
        extent = self._principal_axes_extent(molecule, selection)
        return float(extent[2]) if extent is not None else 0.0

    def calc_length_to_width_ratio(self, molecule, selection) -> float:
        extent = self._principal_axes_extent(molecule, selection)
        if extent is None or extent[1] <= 0:
            return 0.0
        return float(extent[0] / extent[1])

    def calc_bounding_box_volume(self, molecule, selection) -> float:
        extent = self._principal_axes_extent(molecule, selection)
        if extent is None:
            return 0.0
        return float(np.prod(extent))

    def calc_span(self, molecule, selection) -> float:
        """Largest distance from the centroid to any atom."""
        coords, _ = self._coords_and_masses(molecule, selection)
        if coords is None:
            return 0.0
        center = coords.mean(axis=0)
        return float(np.max(np.linalg.norm(coords - center, axis=1)))

    def calc_plane_of_best_fit(self, molecule, selection) -> float:
        """Mean distance of the atoms from their best fit plane, 0 when planar."""
        coords, _ = self._coords_and_masses(molecule, selection)
        if coords is None or len(coords) < 4:
            return 0.0
        centered = coords - coords.mean(axis=0)
        _, _, axes = np.linalg.svd(centered, full_matrices=False)
        normal = axes[2]
        return float(np.mean(np.abs(centered.dot(normal))))

    def calc_mass_weighted_radius_of_gyration(self, molecule, selection) -> float:
        coords, masses = self._coords_and_masses(molecule, selection)
        if coords is None:
            return 0.0
        total_mass = float(np.sum(masses))
        if total_mass <= 0:
            return 0.0
        center = np.average(coords, axis=0, weights=masses)
        squared = np.sum(masses * np.sum((coords - center) ** 2, axis=1))
        return float(np.sqrt(squared / total_mass))

    def calc_mean_atomic_distance(self, molecule, selection) -> float:
        """Mean through space distance between atom pairs (3D Wiener / pairs)."""
        coords, _ = self._coords_and_masses(molecule, selection)
        if coords is None or len(coords) < 2:
            return 0.0
        diff = coords[:, None, :] - coords[None, :, :]
        distances = np.linalg.norm(diff, axis=-1)
        n = len(coords)
        return float(np.sum(np.triu(distances, 1)) / (n * (n - 1) / 2))

    def calc_geometric_wiener_index(self, molecule, selection) -> float:
        """Sum of all through space interatomic distances."""
        coords, _ = self._coords_and_masses(molecule, selection)
        if coords is None or len(coords) < 2:
            return 0.0
        diff = coords[:, None, :] - coords[None, :, :]
        return float(np.sum(np.triu(np.linalg.norm(diff, axis=-1), 1)))

    # ------------------------------------------------ gravitational indices
    def calc_gravitational_index(self, molecule, selection) -> float:
        """Sum of m_i*m_j / r_ij^2 over all atom pairs (Katritzky)."""
        coords, masses = self._coords_and_masses(molecule, selection)
        if coords is None or len(coords) < 2:
            return 0.0
        diff = coords[:, None, :] - coords[None, :, :]
        distances = np.linalg.norm(diff, axis=-1)
        with np.errstate(divide='ignore', invalid='ignore'):
            weights = np.outer(masses, masses) / distances ** 2
        weights[~np.isfinite(weights)] = 0.0
        return float(np.sum(np.triu(weights, 1)))

    def calc_gravitational_index_bonded(self, molecule, selection) -> float:
        """Same, restricted to bonded pairs."""
        chosen = set(selection.atom_indices)
        total = 0.0
        for bond in molecule.bonds:
            i, j = bond.begin_atom_idx, bond.end_atom_idx
            if i not in chosen or j not in chosen:
                continue
            a, b = molecule.atoms[i], molecule.atoms[j]
            if a.x is None or b.x is None:
                continue
            distance = np.linalg.norm(np.array([a.x, a.y, a.z]) - np.array([b.x, b.y, b.z]))
            if distance > 0:
                total += getattr(a, 'mass', 12.0) * getattr(b, 'mass', 12.0) / distance ** 2
        return float(total)

    # -------------------------------------------------------- area/volume
    def calc_surface_to_volume_ratio(self, molecule, selection) -> float:
        volume = self.calc_molecular_volume(molecule, selection)
        if volume <= 0:
            return 0.0
        return float(self.calc_sasa(molecule, selection) / volume)

    def calc_globularity(self, molecule, selection) -> float:
        """Surface of the equivalent sphere over the real surface, 1 for a sphere."""
        volume = self.calc_molecular_volume(molecule, selection)
        surface = self.calc_sasa(molecule, selection)
        if volume <= 0 or surface <= 0:
            return 0.0
        equivalent = (36.0 * np.pi * volume ** 2) ** (1.0 / 3.0)
        return float(equivalent / surface)

    def calc_molecular_density(self, molecule, selection) -> float:
        """Molecular weight per unit volume (g/mol per A^3)."""
        volume = self.calc_molecular_volume(molecule, selection)
        if volume <= 0:
            return 0.0
        mass = sum(getattr(molecule.atoms[i], 'mass', 0.0)
                   for i in selection.atom_indices if i < len(molecule.atoms))
        return float(mass / volume)
