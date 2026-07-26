"""
Shared graph utilities for the descriptor calculators.

The classical topological indices (Wiener, chi, kappa, Balaban, ...) are defined
on the hydrogen suppressed molecular graph, so :class:`SubGraph` drops hydrogens
by default. Everything is computed once per graph and cached, which matters
because a full descriptor run asks for the distance matrix dozens of times.
"""
import numpy as np

# valence electron count, used by the Kier-Hall valence connectivity indices
VALENCE_ELECTRONS = {
    'H': 1, 'Li': 1, 'Be': 2, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7,
    'Na': 1, 'Mg': 2, 'Al': 3, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7,
    'K': 1, 'Ca': 2, 'As': 5, 'Se': 6, 'Br': 7, 'Te': 6, 'I': 7,
}

# atomic polarizability (A^3), Miller & Savchik
ATOMIC_POLARIZABILITY = {
    'H': 0.667, 'C': 1.76, 'N': 1.10, 'O': 0.802, 'F': 0.557, 'Si': 5.38,
    'P': 3.63, 'S': 2.90, 'Cl': 2.18, 'Br': 3.05, 'I': 5.35, 'B': 3.03,
    'Se': 3.77, 'As': 4.31,
}

# van der Waals volume (A^3) derived from Bondi radii
VDW_VOLUME = {
    'H': 7.24, 'C': 20.58, 'N': 15.60, 'O': 14.71, 'F': 13.31, 'Si': 38.79,
    'P': 24.43, 'S': 24.43, 'Cl': 22.45, 'Br': 26.52, 'I': 32.52, 'B': 29.66,
    'Se': 28.73, 'As': 26.52,
}


def valence_electrons(symbol):
    """Number of valence electrons of an element (0 when unknown)."""
    return VALENCE_ELECTRONS.get(symbol, 4)


_GRAPH_CACHE = {}
_GRAPH_CACHE_LIMIT = 8


def get_subgraph(molecule, selection=None, heavy_only=True):
    """Cached SubGraph. A full descriptor run asks for the same graph ~50 times.

    The cached graph keeps a reference to its molecule, so an id() collision with
    a later molecule cannot happen while the entry is alive. """
    if selection is None:
        selection_key = None
    else:
        selection_key = (len(selection.atom_indices), tuple(selection.atom_indices))
    key = (id(molecule), len(molecule.atoms), len(molecule.bonds),
           heavy_only, selection_key)
    graph = _GRAPH_CACHE.get(key)
    if graph is None:
        graph = SubGraph(molecule, selection, heavy_only=heavy_only)
        if len(_GRAPH_CACHE) >= _GRAPH_CACHE_LIMIT:
            _GRAPH_CACHE.clear()
        _GRAPH_CACHE[key] = graph
    return graph


class SubGraph:
    """The selected part of a molecule seen as a plain graph.

    Attributes:
        indices   original atom indices, in ascending order
        atoms     the corresponding Atom objects
        adj       adjacency list over local indices
        edges     list of (local_u, local_v, bond)
    """

    def __init__(self, molecule, selection=None, heavy_only=True):
        n_all = len(molecule.atoms)
        if selection is None:
            chosen = set(range(n_all))
        else:
            chosen = {i for i in selection.atom_indices if 0 <= i < n_all}

        indices = sorted(chosen)
        if heavy_only:
            heavy = [i for i in indices if molecule.atoms[i].symbol != 'H']
            if heavy:  # a pure hydrogen selection keeps its atoms
                indices = heavy

        self.molecule = molecule
        self.indices = indices
        self.n = len(indices)
        self.position = {a: k for k, a in enumerate(indices)}
        self.atoms = [molecule.atoms[i] for i in indices]

        self.adj = [[] for _ in range(self.n)]
        self.edges = []
        for bond in molecule.bonds:
            u = self.position.get(bond.begin_atom_idx)
            v = self.position.get(bond.end_atom_idx)
            if u is None or v is None:
                continue
            self.adj[u].append(v)
            self.adj[v].append(u)
            self.edges.append((u, v, bond))

        self._distance = None
        self._degrees = None
        self._h_counts = None
        self._eccentricity = None
        self._adjacency = None

    # ------------------------------------------------------------ basics
    @property
    def n_bonds(self):
        return len(self.edges)

    @property
    def degrees(self):
        """Number of heavy neighbours of every vertex."""
        if self._degrees is None:
            self._degrees = np.array([len(a) for a in self.adj], dtype=float)
        return self._degrees

    @property
    def h_counts(self):
        """Hydrogens on every vertex, implicit ones and suppressed ones alike."""
        if self._h_counts is None:
            counts = []
            for original in self.indices:
                atom = self.molecule.atoms[original]
                explicit = sum(1 for nb in self.molecule.get_neighbors(original)
                               if self.molecule.atoms[nb].symbol == 'H')
                counts.append(getattr(atom, 'total_h', 0) + explicit)
            self._h_counts = np.array(counts, dtype=float)
        return self._h_counts

    @property
    def adjacency(self):
        if self._adjacency is None:
            matrix = np.zeros((self.n, self.n))
            for u, v, _ in self.edges:
                matrix[u, v] = 1.0
                matrix[v, u] = 1.0
            self._adjacency = matrix
        return self._adjacency

    @property
    def laplacian(self):
        return np.diag(self.degrees) - self.adjacency

    @property
    def distance(self):
        """All pairs shortest path lengths, np.inf between disconnected parts."""
        if self._distance is None:
            self._distance = self._bfs_distances()
        return self._distance

    def _bfs_distances(self):
        n = self.n
        dist = np.full((n, n), np.inf)
        for start in range(n):
            dist[start, start] = 0.0
            frontier = [start]
            depth = 0
            seen = {start}
            while frontier:
                depth += 1
                nxt = []
                for node in frontier:
                    for nb in self.adj[node]:
                        if nb not in seen:
                            seen.add(nb)
                            dist[start, nb] = depth
                            nxt.append(nb)
                frontier = nxt
        return dist

    @property
    def eccentricity(self):
        """Longest shortest path from every vertex (within its component)."""
        if self._eccentricity is None:
            dist = self.distance
            ecc = np.zeros(self.n)
            for i in range(self.n):
                finite = dist[i][np.isfinite(dist[i])]
                ecc[i] = finite.max() if finite.size else 0.0
            self._eccentricity = ecc
        return self._eccentricity

    @property
    def is_connected(self):
        return self.n > 0 and np.isfinite(self.distance).all()

    def finite_pairs(self):
        """Upper triangle mask of vertex pairs that are connected."""
        mask = np.isfinite(self.distance) & ~np.eye(self.n, dtype=bool)
        return np.triu(mask)

    # -------------------------------------------------------- properties
    def property_vector(self, name):
        """Per-vertex property vector used by the autocorrelation descriptors."""
        values = []
        for atom in self.atoms:
            symbol = atom.symbol
            if name == 'mass':
                values.append(getattr(atom, 'mass', 12.011))
            elif name == 'charge':
                values.append(getattr(atom, 'partial_charge', 0.0) or 0.0)
            elif name == 'electronegativity':
                values.append(getattr(atom.element, 'electronegativity', 2.55) or 0.0)
            elif name == 'polarizability':
                values.append(ATOMIC_POLARIZABILITY.get(symbol, 1.76))
            elif name == 'vdw_volume':
                values.append(VDW_VOLUME.get(symbol, 20.58))
            elif name == 'atomic_number':
                values.append(getattr(atom, 'atomic_number', 6))
            else:
                values.append(0.0)
        return np.array(values, dtype=float)

    def relative_property(self, name):
        """Property divided by the carbon value, the usual autocorrelation scaling."""
        carbon = {'mass': 12.011, 'electronegativity': 2.55,
                  'polarizability': 1.76, 'vdw_volume': 20.58,
                  'atomic_number': 6.0}.get(name, 1.0)
        return self.property_vector(name) / carbon if carbon else self.property_vector(name)

    def simple_delta(self):
        """Kier-Hall simple delta: the heavy atom degree."""
        return self.degrees

    def valence_delta(self):
        """Kier-Hall valence delta.

        delta_v = Zv - h for the first two periods and (Zv - h)/(Z - Zv - 1)
        beyond, which is what makes sulfur and the halogens behave. """
        deltas = np.zeros(self.n)
        h_counts = self.h_counts
        for i, atom in enumerate(self.atoms):
            z = getattr(atom, 'atomic_number', 6)
            zv = valence_electrons(atom.symbol)
            h = h_counts[i]
            if z <= 10:
                deltas[i] = max(zv - h, 0.0)
            else:
                denominator = z - zv - 1
                deltas[i] = max((zv - h) / denominator, 0.0) if denominator > 0 else 0.0
        return deltas

    # ------------------------------------------------------------- paths
    def simple_paths(self, length):
        """All simple paths with `length` bonds, each reported once.

        Molecular graphs are small and sparse, so the plain depth first walk is
        fast enough for the path orders the connectivity indices need. """
        if length <= 0:
            return [[i] for i in range(self.n)]
        paths = []

        def walk(path, visited):
            if len(path) == length + 1:
                if path[0] < path[-1]:
                    paths.append(list(path))
                return
            for nb in self.adj[path[-1]]:
                if nb in visited:
                    continue
                visited.add(nb)
                path.append(nb)
                walk(path, visited)
                path.pop()
                visited.remove(nb)

        for start in range(self.n):
            walk([start], {start})
        return paths

    def path_count(self, length):
        return len(self.simple_paths(length))

    def rings(self):
        """SSSR rings of the molecule restricted to the selected atoms."""
        try:
            rings = self.molecule.find_rings()
        except Exception:
            return []
        chosen = set(self.indices)
        return [r for r in rings if all(i in chosen for i in r)]
