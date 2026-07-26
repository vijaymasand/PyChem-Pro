"""
Extended topological (2D) descriptors.

Adds the descriptor families that QSAR work usually expects on top of the
basic indices: Kier-Hall connectivity and shape indices, the modern degree
based indices, distance/detour style indices, graph spectral descriptors and
2D autocorrelations. Everything is computed on the hydrogen suppressed graph,
as the original definitions require.
"""
import numpy as np

from .topological import TopologicalCalculator
from .graph_utils import get_subgraph


class ExtendedTopologicalCalculator(TopologicalCalculator):
    """Topological descriptors, basic set plus the extended families."""

    def graph(self, molecule, selection):
        return get_subgraph(molecule, selection, heavy_only=True)

    # ------------------------------------------------- connectivity (chi)
    def _chi_path(self, graph, order, deltas):
        """Kier-Hall path connectivity index of the given order."""
        safe = np.where(deltas > 0, deltas, np.nan)
        if order == 0:
            values = 1.0 / np.sqrt(safe)
            return float(np.nansum(values))
        total = 0.0
        for path in graph.simple_paths(order):
            product = np.prod([safe[i] for i in path])
            if np.isfinite(product) and product > 0:
                total += 1.0 / np.sqrt(product)
        return float(total)

    def _chi_cluster(self, graph, size, deltas):
        """Cluster connectivity: a central atom with `size` neighbours."""
        from itertools import combinations
        safe = np.where(deltas > 0, deltas, np.nan)
        total = 0.0
        for center in range(graph.n):
            neighbors = graph.adj[center]
            if len(neighbors) < size:
                continue
            for group in combinations(neighbors, size):
                product = safe[center] * np.prod([safe[i] for i in group])
                if np.isfinite(product) and product > 0:
                    total += 1.0 / np.sqrt(product)
        return float(total)

    def calc_chi0n(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_path(graph, 0, graph.simple_delta())

    def calc_chi1n(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_path(graph, 1, graph.simple_delta())

    def calc_chi2n(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_path(graph, 2, graph.simple_delta())

    def calc_chi3n(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_path(graph, 3, graph.simple_delta())

    def calc_chi4n(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_path(graph, 4, graph.simple_delta())

    def calc_chi0v(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_path(graph, 0, graph.valence_delta())

    def calc_chi1v(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_path(graph, 1, graph.valence_delta())

    def calc_chi2v(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_path(graph, 2, graph.valence_delta())

    def calc_chi3v(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_path(graph, 3, graph.valence_delta())

    def calc_chi4v(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_path(graph, 4, graph.valence_delta())

    def calc_chi3_cluster(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_cluster(graph, 3, graph.simple_delta())

    def calc_chi4_cluster(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return self._chi_cluster(graph, 4, graph.simple_delta())

    # ------------------------------------------------------ Kier shape
    # Kier alpha values, the covalent radius of an atom relative to sp3 carbon
    # (Hall & Kier). Anything not listed falls back to the radius ratio.
    KIER_ALPHA = {
        ('C', 'sp3'): 0.00, ('C', 'sp2'): -0.13, ('C', 'sp'): -0.22,
        ('N', 'sp3'): -0.04, ('N', 'sp2'): -0.20, ('N', 'sp'): -0.29,
        ('O', 'sp3'): -0.04, ('O', 'sp2'): -0.20,
        ('P', 'sp3'): 0.43, ('P', 'sp2'): 0.30,
        ('S', 'sp3'): 0.35, ('S', 'sp2'): 0.22,
        ('F', None): -0.07, ('Cl', None): 0.29,
        ('Br', None): 0.48, ('I', None): 0.73,
    }

    def _kier_alpha(self, graph):
        """Kier alpha correction: covalent radius relative to sp3 carbon."""
        try:
            graph.molecule.assign_hybridization()
        except Exception:
            pass
        alpha = 0.0
        for atom in graph.atoms:
            symbol = atom.symbol
            hybridization = getattr(atom, 'hybridization', None)
            value = self.KIER_ALPHA.get((symbol, hybridization))
            if value is None:
                value = self.KIER_ALPHA.get((symbol, None))
            if value is None:
                radius = getattr(atom.element, 'covalent_radius', 0.77) or 0.77
                value = radius / 0.77 - 1.0
            alpha += value
        return alpha

    def _kappa(self, graph, order):
        a = self._kier_alpha(graph)
        A = graph.n
        if A < order + 1:
            return 0.0
        paths = graph.path_count(order)
        denominator = (paths + a) ** 2
        if denominator <= 0:
            return 0.0
        if order == 1:
            numerator = (A + a) * (A + a - 1) ** 2
        elif order == 2:
            numerator = (A + a - 1) * (A + a - 2) ** 2
        else:
            if A % 2:
                numerator = (A + a - 1) * (A + a - 3) ** 2
            else:
                numerator = (A + a - 3) * (A + a - 2) ** 2
        return float(max(numerator / denominator, 0.0))

    def calc_kappa1_alpha(self, molecule, selection) -> float:
        return self._kappa(self.graph(molecule, selection), 1)

    def calc_kappa2_alpha(self, molecule, selection) -> float:
        return self._kappa(self.graph(molecule, selection), 2)

    def calc_kappa3_alpha(self, molecule, selection) -> float:
        return self._kappa(self.graph(molecule, selection), 3)

    def calc_kier_flexibility(self, molecule, selection) -> float:
        """Kier molecular flexibility index Phi = kappa1 * kappa2 / A."""
        graph = self.graph(molecule, selection)
        if graph.n == 0:
            return 0.0
        return self._kappa(graph, 1) * self._kappa(graph, 2) / graph.n

    # --------------------------------------------- degree based indices
    def calc_abc_index(self, molecule, selection) -> float:
        """Atom-bond connectivity index."""
        graph = self.graph(molecule, selection)
        degrees = graph.degrees
        total = 0.0
        for u, v, _ in graph.edges:
            du, dv = degrees[u], degrees[v]
            if du > 0 and dv > 0 and du + dv > 2:
                total += np.sqrt((du + dv - 2.0) / (du * dv))
        return float(total)

    def calc_augmented_zagreb_index(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        degrees = graph.degrees
        total = 0.0
        for u, v, _ in graph.edges:
            du, dv = degrees[u], degrees[v]
            if du + dv > 2:
                total += (du * dv / (du + dv - 2.0)) ** 3
        return float(total)

    def calc_forgotten_index(self, molecule, selection) -> float:
        """Forgotten index F = sum of cubed degrees."""
        graph = self.graph(molecule, selection)
        return float(np.sum(graph.degrees ** 3))

    def calc_modified_zagreb_index(self, molecule, selection) -> float:
        """Modified first Zagreb index = sum 1/degree^2."""
        graph = self.graph(molecule, selection)
        degrees = graph.degrees
        return float(sum(1.0 / (d * d) for d in degrees if d > 0))

    def calc_sum_connectivity_index(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        degrees = graph.degrees
        return float(sum(1.0 / np.sqrt(degrees[u] + degrees[v])
                         for u, v, _ in graph.edges
                         if degrees[u] + degrees[v] > 0))

    def calc_geometric_arithmetic_index(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        degrees = graph.degrees
        total = 0.0
        for u, v, _ in graph.edges:
            du, dv = degrees[u], degrees[v]
            if du + dv > 0:
                total += 2.0 * np.sqrt(du * dv) / (du + dv)
        return float(total)

    def calc_harmonic_index(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        degrees = graph.degrees
        return float(sum(2.0 / (degrees[u] + degrees[v])
                         for u, v, _ in graph.edges
                         if degrees[u] + degrees[v] > 0))

    def calc_narumi_simple_index(self, molecule, selection) -> float:
        """Log of the product of vertex degrees (the plain product overflows)."""
        graph = self.graph(molecule, selection)
        degrees = [d for d in graph.degrees if d > 0]
        if not degrees:
            return 0.0
        return float(np.sum(np.log(degrees)))

    def calc_narumi_harmonic_index(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        degrees = [d for d in graph.degrees if d > 0]
        if not degrees:
            return 0.0
        return float(len(degrees) / sum(1.0 / d for d in degrees))

    def calc_narumi_geometric_index(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        degrees = [d for d in graph.degrees if d > 0]
        if not degrees:
            return 0.0
        return float(np.exp(np.mean(np.log(degrees))))

    # ------------------------------------------- distance based indices
    def calc_schultz_index(self, molecule, selection) -> float:
        """Schultz molecular topological index, sum of degree*(A+D) rows."""
        graph = self.graph(molecule, selection)
        if graph.n == 0:
            return 0.0
        dist = np.where(np.isfinite(graph.distance), graph.distance, 0.0)
        combined = graph.adjacency + dist
        return float(np.sum(graph.degrees[:, None] * combined))

    def calc_gutman_index(self, molecule, selection) -> float:
        """Gutman (Schultz of the second kind) index."""
        graph = self.graph(molecule, selection)
        if graph.n == 0:
            return 0.0
        dist = np.where(np.isfinite(graph.distance), graph.distance, 0.0)
        degrees = graph.degrees
        weights = np.outer(degrees, degrees) * dist
        return float(np.sum(np.triu(weights, 1)))

    def calc_szeged_index(self, molecule, selection) -> float:
        """Szeged index, the distance based generalisation of Wiener."""
        graph = self.graph(molecule, selection)
        dist = graph.distance
        total = 0.0
        for u, v, _ in graph.edges:
            closer_u = np.sum(dist[u] < dist[v])
            closer_v = np.sum(dist[v] < dist[u])
            total += closer_u * closer_v
        return float(total)

    def calc_pi_index(self, molecule, selection) -> float:
        """Padmakar-Ivan index (vertex version)."""
        graph = self.graph(molecule, selection)
        dist = graph.distance
        total = 0.0
        for u, v, _ in graph.edges:
            total += np.sum(dist[u] < dist[v]) + np.sum(dist[v] < dist[u])
        return float(total)

    def calc_mean_wiener_index(self, molecule, selection) -> float:
        """Wiener index per vertex pair."""
        graph = self.graph(molecule, selection)
        mask = graph.finite_pairs()
        pairs = int(np.sum(mask))
        if pairs == 0:
            return 0.0
        return float(np.sum(graph.distance[mask]) / pairs)

    def calc_topological_radius(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        if graph.n == 0:
            return 0.0
        return float(np.min(graph.eccentricity))

    def calc_average_eccentricity(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        if graph.n == 0:
            return 0.0
        return float(np.mean(graph.eccentricity))

    def calc_petitjean_index(self, molecule, selection) -> float:
        """Topological shape index (diameter - radius) / radius."""
        graph = self.graph(molecule, selection)
        if graph.n == 0:
            return 0.0
        radius = float(np.min(graph.eccentricity))
        diameter = float(np.max(graph.eccentricity))
        return (diameter - radius) / radius if radius > 0 else 0.0

    def calc_eccentric_distance_sum(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        dist = np.where(np.isfinite(graph.distance), graph.distance, 0.0)
        return float(np.sum(graph.eccentricity * dist.sum(axis=1)))

    def calc_kirchhoff_index(self, molecule, selection) -> float:
        """Resistance distance (quasi-Wiener) index from Laplacian eigenvalues."""
        graph = self.graph(molecule, selection)
        if graph.n < 2 or not graph.is_connected:
            return 0.0
        eigenvalues = np.sort(np.linalg.eigvalsh(graph.laplacian))[1:]
        positive = eigenvalues[eigenvalues > 1e-9]
        if positive.size == 0:
            return 0.0
        return float(graph.n * np.sum(1.0 / positive))

    # --------------------------------------------------------- spectral
    def _adjacency_eigenvalues(self, graph):
        if graph.n == 0:
            return np.array([0.0])
        return np.linalg.eigvalsh(graph.adjacency)

    def calc_graph_energy(self, molecule, selection) -> float:
        """Graph energy, the sum of absolute adjacency eigenvalues."""
        graph = self.graph(molecule, selection)
        return float(np.sum(np.abs(self._adjacency_eigenvalues(graph))))

    def calc_spectral_radius(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return float(np.max(np.abs(self._adjacency_eigenvalues(graph))))

    def calc_estrada_index(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        eigenvalues = np.clip(self._adjacency_eigenvalues(graph), -50, 50)
        return float(np.sum(np.exp(eigenvalues)))

    def calc_laplacian_spectral_radius(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        if graph.n == 0:
            return 0.0
        return float(np.max(np.linalg.eigvalsh(graph.laplacian)))

    def calc_algebraic_connectivity(self, molecule, selection) -> float:
        """Fiedler value, zero when the selection falls apart."""
        graph = self.graph(molecule, selection)
        if graph.n < 2:
            return 0.0
        eigenvalues = np.sort(np.linalg.eigvalsh(graph.laplacian))
        return float(max(eigenvalues[1], 0.0))

    def calc_log_spanning_tree_count(self, molecule, selection) -> float:
        """Log10 of the number of spanning trees (Kirchhoff matrix tree theorem)."""
        graph = self.graph(molecule, selection)
        if graph.n < 2 or not graph.is_connected:
            return 0.0
        eigenvalues = np.sort(np.linalg.eigvalsh(graph.laplacian))[1:]
        positive = eigenvalues[eigenvalues > 1e-9]
        if positive.size == 0:
            return 0.0
        return float(np.sum(np.log10(positive)) - np.log10(graph.n))

    def _walk_count(self, graph, order):
        if graph.n == 0:
            return 0.0
        power = np.linalg.matrix_power(graph.adjacency, order)
        return float(np.sum(power))

    def calc_walk_count_2(self, molecule, selection) -> float:
        return self._walk_count(self.graph(molecule, selection), 2)

    def calc_walk_count_3(self, molecule, selection) -> float:
        return self._walk_count(self.graph(molecule, selection), 3)

    def calc_walk_count_4(self, molecule, selection) -> float:
        return self._walk_count(self.graph(molecule, selection), 4)

    def calc_self_returning_walk_3(self, molecule, selection) -> float:
        """Trace of A^3, six times the number of triangles."""
        graph = self.graph(molecule, selection)
        if graph.n == 0:
            return 0.0
        return float(np.trace(np.linalg.matrix_power(graph.adjacency, 3)))

    # -------------------------------------------------- autocorrelation
    def _moreau_broto(self, graph, prop, lag):
        values = graph.relative_property(prop)
        dist = graph.distance
        mask = np.triu(dist == lag)
        if not mask.any():
            return 0.0
        products = np.outer(values, values)
        return float(np.sum(products[mask]))

    def _moran(self, graph, prop, lag):
        values = graph.property_vector(prop)
        n = graph.n
        if n < 2:
            return 0.0
        mean = values.mean()
        deviation = values - mean
        variance = np.sum(deviation ** 2) / n
        if variance <= 1e-12:
            return 0.0
        mask = np.triu(graph.distance == lag)
        pairs = int(np.sum(mask))
        if pairs == 0:
            return 0.0
        covariance = np.sum(np.outer(deviation, deviation)[mask]) / pairs
        return float(covariance / variance)

    def _geary(self, graph, prop, lag):
        values = graph.property_vector(prop)
        n = graph.n
        if n < 2:
            return 0.0
        mean = values.mean()
        variance = np.sum((values - mean) ** 2) / (n - 1)
        if variance <= 1e-12:
            return 0.0
        mask = np.triu(graph.distance == lag)
        pairs = int(np.sum(mask))
        if pairs == 0:
            return 0.0
        differences = (values[:, None] - values[None, :]) ** 2
        return float(np.sum(differences[mask]) / (2 * pairs) / variance)

    def calc_ats_en_1(self, molecule, selection) -> float:
        return self._moreau_broto(self.graph(molecule, selection), 'electronegativity', 1)

    def calc_ats_en_2(self, molecule, selection) -> float:
        return self._moreau_broto(self.graph(molecule, selection), 'electronegativity', 2)

    def calc_ats_en_3(self, molecule, selection) -> float:
        return self._moreau_broto(self.graph(molecule, selection), 'electronegativity', 3)

    def calc_ats_polarizability_1(self, molecule, selection) -> float:
        return self._moreau_broto(self.graph(molecule, selection), 'polarizability', 1)

    def calc_ats_polarizability_2(self, molecule, selection) -> float:
        return self._moreau_broto(self.graph(molecule, selection), 'polarizability', 2)

    def calc_ats_polarizability_3(self, molecule, selection) -> float:
        return self._moreau_broto(self.graph(molecule, selection), 'polarizability', 3)

    def calc_ats_volume_1(self, molecule, selection) -> float:
        return self._moreau_broto(self.graph(molecule, selection), 'vdw_volume', 1)

    def calc_ats_volume_2(self, molecule, selection) -> float:
        return self._moreau_broto(self.graph(molecule, selection), 'vdw_volume', 2)

    def calc_ats_volume_3(self, molecule, selection) -> float:
        return self._moreau_broto(self.graph(molecule, selection), 'vdw_volume', 3)

    def calc_moran_en_1(self, molecule, selection) -> float:
        return self._moran(self.graph(molecule, selection), 'electronegativity', 1)

    def calc_moran_en_2(self, molecule, selection) -> float:
        return self._moran(self.graph(molecule, selection), 'electronegativity', 2)

    def calc_moran_mass_1(self, molecule, selection) -> float:
        return self._moran(self.graph(molecule, selection), 'mass', 1)

    def calc_geary_en_1(self, molecule, selection) -> float:
        return self._geary(self.graph(molecule, selection), 'electronegativity', 1)

    def calc_geary_en_2(self, molecule, selection) -> float:
        return self._geary(self.graph(molecule, selection), 'electronegativity', 2)

    def calc_geary_mass_1(self, molecule, selection) -> float:
        return self._geary(self.graph(molecule, selection), 'mass', 1)

    # ------------------------------------------- Galvez topological charge
    def _charge_transfer_matrix(self, graph):
        dist = graph.distance
        with np.errstate(divide='ignore', invalid='ignore'):
            reciprocal = np.where((dist > 0) & np.isfinite(dist), 1.0 / dist ** 2, 0.0)
        product = graph.adjacency.dot(reciprocal)
        return product - product.T

    def _galvez_index(self, graph, lag):
        if graph.n < 2:
            return 0.0
        transfer = self._charge_transfer_matrix(graph)
        mask = np.triu(graph.distance == lag)
        if not mask.any():
            return 0.0
        return float(np.sum(np.abs(transfer[mask])))

    def calc_topological_charge_1(self, molecule, selection) -> float:
        return self._galvez_index(self.graph(molecule, selection), 1)

    def calc_topological_charge_2(self, molecule, selection) -> float:
        return self._galvez_index(self.graph(molecule, selection), 2)

    def calc_topological_charge_3(self, molecule, selection) -> float:
        return self._galvez_index(self.graph(molecule, selection), 3)

    def calc_mean_topological_charge(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        if graph.n < 2:
            return 0.0
        total = sum(self._galvez_index(graph, lag) for lag in (1, 2, 3))
        return total / (graph.n - 1)

    # ------------------------------------------------ information content
    def _shannon_entropy(self, labels):
        if not labels:
            return 0.0
        total = len(labels)
        counts = {}
        for label in labels:
            counts[label] = counts.get(label, 0) + 1
        entropy = -sum((c / total) * np.log2(c / total) for c in counts.values())
        return max(entropy, 0.0)  # a single class gives -0.0 otherwise

    def calc_atom_type_information(self, molecule, selection) -> float:
        """Shannon entropy of the element distribution, per atom."""
        graph = self.graph(molecule, selection)
        return float(self._shannon_entropy([a.symbol for a in graph.atoms]))

    def calc_total_atom_information(self, molecule, selection) -> float:
        """Total information content on the element distribution."""
        graph = self.graph(molecule, selection)
        return float(graph.n * self._shannon_entropy([a.symbol for a in graph.atoms]))

    def calc_bond_type_information(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        labels = []
        for _, _, bond in graph.edges:
            if bond.is_aromatic:
                labels.append('ar')
            elif bond.is_triple:
                labels.append('3')
            elif bond.is_double:
                labels.append('2')
            else:
                labels.append('1')
        return float(self._shannon_entropy(labels))

    def calc_distance_information(self, molecule, selection) -> float:
        """Bonchev-Trinajstic information index on the distance distribution."""
        graph = self.graph(molecule, selection)
        mask = graph.finite_pairs()
        distances = graph.distance[mask]
        if distances.size == 0:
            return 0.0
        labels = [int(d) for d in distances]
        return float(len(labels) * self._shannon_entropy(labels))

    def calc_vertex_degree_information(self, molecule, selection) -> float:
        graph = self.graph(molecule, selection)
        return float(self._shannon_entropy([int(d) for d in graph.degrees]))
