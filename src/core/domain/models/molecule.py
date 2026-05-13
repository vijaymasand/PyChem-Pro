"""
Molecule class — Molecular graph with atoms, bonds, ring perception, and utilities.
"""

from collections import deque
from src.core.domain.models.atom import Atom
from src.core.domain.models.bond import Bond, BondType


class Molecule:
    """
    Represents a molecular graph with atoms and bonds.

    Provides graph operations: adjacency queries, ring perception (SSSR),
    fragment detection, molecular formula computation, and hybridization assignment.
    """

    def __init__(self, name=""):
        self.name = name
        self.atoms = []          # list of Atom
        self.bonds = []          # list of Bond
        self._adjacency = {}     # atom_idx -> list of (neighbor_idx, bond_idx)
        self._rings = None       # cached SSSR
        self.properties = {}     # arbitrary key-value properties
        self._suppress_ring_invalidation = False  # bulk-load flag

    # ─── Graph Building ─────────────────────────────────────────────

    def add_atom(self, atom):
        """Add an atom to the molecule. Returns the atom index."""
        atom.index = len(self.atoms)
        self.atoms.append(atom)
        self._adjacency[atom.index] = []
        if not self._suppress_ring_invalidation:
            self._rings = None
        return atom.index

    def add_bond(self, begin_idx, end_idx, bond_type=BondType.SINGLE, stereo=0):
        """Add a bond between two atoms. Returns the bond index."""
        if begin_idx == end_idx:
            raise ValueError("Cannot bond atom to itself")
        # Check if bond already exists
        for neighbor_idx, bond_idx in self._adjacency.get(begin_idx, []):
            if neighbor_idx == end_idx:
                return bond_idx  # Bond already exists
        bond = Bond(begin_idx, end_idx, bond_type, stereo)
        bond.index = len(self.bonds)
        self.bonds.append(bond)
        self._adjacency[begin_idx].append((end_idx, bond.index))
        self._adjacency[end_idx].append((begin_idx, bond.index))
        if not self._suppress_ring_invalidation:
            self._rings = None
        return bond.index

    def begin_bulk_load(self):
        """Suppress ring-cache invalidation during bulk atom/bond insertion."""
        self._suppress_ring_invalidation = True

    def end_bulk_load(self):
        """Re-enable ring-cache invalidation and clear cached rings."""
        self._suppress_ring_invalidation = False
        self._rings = None

    # ─── Graph Editing (Remove) ─────────────────────────────────────

    def remove_atoms(self, indices_to_remove):
        """
        Remove a set of atoms (and all bonds touching them) from the molecule.

        After removal, the remaining atoms and bonds are re-indexed so that
        ``self.atoms[i].index == i`` and ``self.bonds[j].index == j`` hold
        for every surviving element.  The adjacency map is rebuilt from
        scratch.

        Args:
            indices_to_remove: Iterable of atom indices to delete.

        Returns:
            int: Number of atoms actually removed.
        """
        remove_set = set(indices_to_remove)
        if not remove_set:
            return 0

        # --- 1. Identify surviving atoms and build old→new index map ---
        old_to_new = {}          # old_atom_idx → new_atom_idx
        new_atoms = []
        for old_idx, atom in enumerate(self.atoms):
            if old_idx not in remove_set:
                new_idx = len(new_atoms)
                old_to_new[old_idx] = new_idx
                atom.index = new_idx
                new_atoms.append(atom)

        # --- 2. Keep only bonds whose *both* endpoints survive ---
        new_bonds = []
        for bond in self.bonds:
            if bond.begin_atom_idx in remove_set or bond.end_atom_idx in remove_set:
                continue  # discard — touches a deleted atom
            bond.begin_atom_idx = old_to_new[bond.begin_atom_idx]
            bond.end_atom_idx = old_to_new[bond.end_atom_idx]
            bond.index = len(new_bonds)
            new_bonds.append(bond)

        # --- 3. Rebuild adjacency map ---
        new_adjacency = {i: [] for i in range(len(new_atoms))}
        for bond in new_bonds:
            new_adjacency[bond.begin_atom_idx].append(
                (bond.end_atom_idx, bond.index))
            new_adjacency[bond.end_atom_idx].append(
                (bond.begin_atom_idx, bond.index))

        # --- 4. Commit ---
        removed_count = len(self.atoms) - len(new_atoms)
        self.atoms = new_atoms
        self.bonds = new_bonds
        self._adjacency = new_adjacency
        self._rings = None  # invalidate ring cache
        return removed_count

    def remove_bonds(self, bond_indices_to_remove):
        """
        Remove a set of bonds from the molecule without removing any atoms.

        Surviving bonds are re-indexed and the adjacency map is rebuilt.

        Args:
            bond_indices_to_remove: Iterable of bond indices to delete.

        Returns:
            int: Number of bonds actually removed.
        """
        remove_set = set(bond_indices_to_remove)
        if not remove_set:
            return 0

        new_bonds = []
        for bond in self.bonds:
            if bond.index in remove_set:
                continue
            bond.index = len(new_bonds)
            new_bonds.append(bond)

        # Rebuild adjacency
        new_adjacency = {i: [] for i in range(len(self.atoms))}
        for bond in new_bonds:
            new_adjacency[bond.begin_atom_idx].append(
                (bond.end_atom_idx, bond.index))
            new_adjacency[bond.end_atom_idx].append(
                (bond.begin_atom_idx, bond.index))

        removed_count = len(self.bonds) - len(new_bonds)
        self.bonds = new_bonds
        self._adjacency = new_adjacency
        self._rings = None
        return removed_count

    # ─── Graph Queries ──────────────────────────────────────────────

    def get_atom(self, idx):
        """Get atom by index."""
        return self.atoms[idx]

    def get_bond(self, idx):
        """Get bond by index."""
        return self.bonds[idx]

    def get_bond_between(self, atom1_idx, atom2_idx):
        """Get the bond between two atoms, or None."""
        for neighbor_idx, bond_idx in self._adjacency.get(atom1_idx, []):
            if neighbor_idx == atom2_idx:
                return self.bonds[bond_idx]
        return None

    def get_neighbors(self, atom_idx):
        """Get list of neighbor atom indices."""
        return [n for n, _ in self._adjacency.get(atom_idx, [])]

    def get_neighbor_bonds(self, atom_idx):
        """Get list of (neighbor_idx, bond) tuples for an atom."""
        result = []
        for neighbor_idx, bond_idx in self._adjacency.get(atom_idx, []):
            result.append((neighbor_idx, self.bonds[bond_idx]))
        return result

    def degree(self, atom_idx):
        """Number of explicit bonds to an atom."""
        return len(self._adjacency.get(atom_idx, []))

    @property
    def num_atoms(self):
        return len(self.atoms)

    @property
    def num_bonds(self):
        return len(self.bonds)

    @property
    def count_heavy_atoms(self):
        """Count non-hydrogen atoms."""
        return sum(1 for a in self.atoms if a.symbol != 'H')

    def propagate_aromaticity(self):
        """
        Mark atoms as aromatic if they are part of an aromatic bond.
        This is useful for formats like MOL2 and SDF that mark aromatic bonds explicitly.
        """
        from src.core.domain.models.bond import BondType
        for bond in self.bonds:
            if bond.bond_type == BondType.AROMATIC:
                self.atoms[bond.begin_atom_idx].is_aromatic = True
                self.atoms[bond.end_atom_idx].is_aromatic = True

    # ─── Ring Perception (SSSR) ─────────────────────────────────────

    def find_rings(self, max_size=20):
        """
        Find the Smallest Set of Smallest Rings (SSSR) using a BFS-based approach.
        Returns a list of rings, each ring is a list of atom indices.
        """
        if self._rings is not None:
            return self._rings

        n = len(self.atoms)
        if n == 0:
            self._rings = []
            return self._rings

        # Expected number of rings from Euler's formula
        # For a connected graph: num_rings = num_bonds - num_atoms + 1
        # For multiple fragments: num_rings = num_bonds - num_atoms + num_fragments
        num_fragments = len(self.get_fragments())
        expected_rings = len(self.bonds) - len(self.atoms) + num_fragments

        if expected_rings <= 0:
            self._rings = []
            return self._rings

        rings = []
        # BFS from each atom to find cycles
        for start in range(n):
            # BFS to find shortest paths that form rings
            visited = {start: None}  # atom -> parent
            queue = deque([(start, -1)])  # (atom, parent_bond_idx)

            while queue and len(rings) < expected_rings * 3:
                current, parent_bond = queue.popleft()

                for neighbor, bond_idx in self._adjacency.get(current, []):
                    if bond_idx == parent_bond:
                        continue

                    if neighbor in visited:
                        # Found a ring — trace back from both sides
                        ring = self._trace_ring(visited, current, neighbor, start)
                        if ring and len(ring) <= max_size:
                            ring_set = frozenset(ring)
                            # Avoid duplicates
                            if not any(frozenset(r) == ring_set for r in rings):
                                rings.append(ring)
                    else:
                        visited[neighbor] = current
                        queue.append((neighbor, bond_idx))

        # Keep only the SSSR (smallest independent set)
        rings.sort(key=len)
        sssr = []
        used_edges = set()

        for ring in rings:
            edges = set()
            for i in range(len(ring)):
                a, b = ring[i], ring[(i + 1) % len(ring)]
                edge = (min(a, b), max(a, b))
                edges.add(edge)

            # Check if this ring contributes a new edge
            if not edges.issubset(used_edges):
                sssr.append(ring)
                used_edges.update(edges)
                if len(sssr) >= expected_rings:
                    break

        # Mark ring bonds and atoms
        for ring in sssr:
            for i in range(len(ring)):
                a_idx = ring[i]
                b_idx = ring[(i + 1) % len(ring)]
                bond = self.get_bond_between(a_idx, b_idx)
                if bond:
                    bond.is_in_ring = True
                self.atoms[a_idx].ring_bonds.add(
                    bond.index if bond else -1
                )

        self._rings = sssr
        return self._rings

    def _trace_ring(self, visited, node_a, node_b, start):
        """Trace a ring from BFS tree when a cycle is detected."""
        path_a = []
        current = node_a
        while current is not None:
            path_a.append(current)
            current = visited.get(current)

        path_b = []
        current = node_b
        while current is not None:
            path_b.append(current)
            current = visited.get(current)

        # Find common ancestor
        set_a = set(path_a)
        common = None
        for node in path_b:
            if node in set_a:
                common = node
                break

        if common is None:
            return None

        # Build ring: path from common to node_a + reverse path from common to node_b
        ring = []
        for node in path_a:
            ring.append(node)
            if node == common:
                break

        rev_path = []
        for node in path_b:
            if node == common:
                break
            rev_path.append(node)

        ring.extend(reversed(rev_path))
        return ring if len(ring) >= 3 else None

    def is_in_ring(self, atom_idx):
        """Check if an atom is in any ring."""
        rings = self.find_rings()
        return any(atom_idx in ring for ring in rings)

    def get_ring_size(self, atom_idx):
        """Get the smallest ring size containing this atom."""
        rings = self.find_rings()
        sizes = [len(ring) for ring in rings if atom_idx in ring]
        return min(sizes) if sizes else 0

    # ─── Fragment Detection ─────────────────────────────────────────

    def get_fragments(self):
        """
        Get connected components as lists of atom indices.
        Uses BFS traversal.
        """
        visited = set()
        fragments = []

        for start in range(len(self.atoms)):
            if start in visited:
                continue
            fragment = []
            queue = deque([start])
            while queue:
                node = queue.popleft()
                if node in visited:
                    continue
                visited.add(node)
                fragment.append(node)
                for neighbor in self.get_neighbors(node):
                    if neighbor not in visited:
                        queue.append(neighbor)
            fragments.append(fragment)

        return fragments

    # ─── Hybridization ──────────────────────────────────────────────

    def assign_hybridization(self):
        """
        Assign hybridization (sp, sp2, sp3) for each atom based on bond orders.
        """
        for atom in self.atoms:
            if atom.symbol == 'H':
                atom.hybridization = 's'
                continue

            neighbors = self.get_neighbor_bonds(atom.index)
            total_bond_order = sum(b.order for _, b in neighbors)
            num_neighbors = len(neighbors)

            # Check for triple bond
            has_triple = any(b.is_triple for _, b in neighbors)
            has_double = any(b.is_double for _, b in neighbors)
            num_doubles = sum(1 for _, b in neighbors if b.is_double)

            if has_triple or (num_doubles >= 2 and num_neighbors == 2):
                atom.hybridization = 'sp'
            elif has_double or getattr(atom, 'is_aromatic', False):
                atom.hybridization = 'sp2'
            else:
                atom.hybridization = 'sp3'

    def perceive_aromaticity(self):
        """
        Identify aromatic rings and atoms using Hückel's 4n+2 rule.
        Marks atoms and bonds as aromatic.
        """
        rings = self.find_rings()
        for ring in rings:
            if len(ring) not in (5, 6):
                continue
            
            # Simple Hückel-like check
            potential_aromatic = True
            pi_electrons = 0
            
            for atom_idx in ring:
                atom = self.atoms[atom_idx]
                neighbors = self.get_neighbor_bonds(atom_idx)
                num_doubles = sum(1 for _, b in neighbors if b.is_double)
                
                if atom.symbol == 'C':
                    if num_doubles >= 1 or self.is_in_ring(atom_idx):
                        pi_electrons += 1
                    else:
                        potential_aromatic = False
                elif atom.symbol == 'N':
                    # Pyrrole-type (2e) or Pyridine-type (1e)
                    if num_doubles >= 1:
                        pi_electrons += 1
                    else:
                        pi_electrons += 2
                elif atom.symbol in ('O', 'S'):
                    pi_electrons += 2
                else:
                    potential_aromatic = False
                
                if not potential_aromatic:
                    break
            
            if potential_aromatic and (pi_electrons % 4 == 2):
                # Mark as aromatic
                for atom_idx in ring:
                    self.atoms[atom_idx].is_aromatic = True
                for i in range(len(ring)):
                    a, b = ring[i], ring[(i+1)%len(ring)]
                    bond = self.get_bond_between(a, b)
                    if bond:
                        from src.core.domain.models.bond import BondType
                        bond.bond_type = BondType.AROMATIC

    # ─── SYBYL Atom Types ──────────────────────────────────────────

    def assign_sybyl_types(self):
        """
        Assign SYBYL atom types for MOL2 output.
        """
        self.assign_hybridization()

        for atom in self.atoms:
            sym = atom.symbol

            if sym == 'H':
                atom.sybyl_type = 'H'
            elif sym == 'C':
                if atom.is_aromatic:
                    atom.sybyl_type = 'C.ar'
                elif atom.hybridization == 'sp':
                    atom.sybyl_type = 'C.1'
                elif atom.hybridization == 'sp2':
                    atom.sybyl_type = 'C.2'
                else:
                    atom.sybyl_type = 'C.3'
            elif sym == 'N':
                if atom.is_aromatic:
                    atom.sybyl_type = 'N.ar'
                elif atom.hybridization == 'sp':
                    atom.sybyl_type = 'N.1'
                elif atom.hybridization == 'sp2':
                    # Check if amide
                    neighbors = self.get_neighbor_bonds(atom.index)
                    is_amide = False
                    for n_idx, bond in neighbors:
                        if self.atoms[n_idx].symbol == 'C':
                            # Check if this carbon is a carbonyl (has C=O)
                            for nn_idx, nb in self.get_neighbor_bonds(n_idx):
                                if self.atoms[nn_idx].symbol == 'O' and nb.is_double:
                                    is_amide = True
                                    # Set bond type to AMIDE for fingerprinting logic
                                    from src.core.domain.models.bond import BondType
                                    bond.bond_type = BondType.AMIDE
                    atom.sybyl_type = 'N.am' if is_amide else 'N.2'
                elif atom.formal_charge == 1:
                    atom.sybyl_type = 'N.4'
                else:
                    atom.sybyl_type = 'N.3'
            elif sym == 'O':
                if atom.is_aromatic:
                    atom.sybyl_type = 'O.ar' if atom.is_aromatic else 'O.3'
                elif atom.hybridization == 'sp2':
                    atom.sybyl_type = 'O.2'
                elif atom.formal_charge == -1:
                    atom.sybyl_type = 'O.co2'
                else:
                    atom.sybyl_type = 'O.3'
            elif sym == 'S':
                if atom.is_aromatic:
                    atom.sybyl_type = 'S.ar'
                elif atom.hybridization == 'sp2':
                    atom.sybyl_type = 'S.2'
                else:
                    atom.sybyl_type = 'S.3'
            elif sym == 'P':
                atom.sybyl_type = 'P.3'
            elif sym in ('F', 'Cl', 'Br', 'I'):
                atom.sybyl_type = sym
            else:
                atom.sybyl_type = sym

    # ─── Molecular Formula ──────────────────────────────────────────

    def molecular_formula(self):
        """
        Generate molecular formula string (Hill system: C first, H second, then alphabetical).
        """
        counts = {}
        for atom in self.atoms:
            sym = atom.symbol
            counts[sym] = counts.get(sym, 0) + 1
            # Add implicit hydrogens
            if atom.total_h > 0:
                counts['H'] = counts.get('H', 0) + atom.total_h

        parts = []
        # Hill system: C first, then H, then alphabetical
        if 'C' in counts:
            c_count = counts.pop('C')
            parts.append(f"C{c_count}" if c_count > 1 else 'C')
            if 'H' in counts:
                h_count = counts.pop('H')
                parts.append(f"H{h_count}" if h_count > 1 else 'H')

        for sym in sorted(counts.keys()):
            c = counts[sym]
            parts.append(f"{sym}{c}" if c > 1 else sym)

        return ''.join(parts)

    def molecular_weight(self):
        """Calculate molecular weight including implicit hydrogens."""
        mw = 0.0
        for atom in self.atoms:
            mw += atom.mass
            mw += atom.total_h * 1.008  # hydrogen mass
        return mw

    def total_charge(self):
        """Sum of all formal charges."""
        return sum(a.formal_charge for a in self.atoms)

    # ─── Explicit Hydrogen Management ───────────────────────────────

    def add_explicit_hydrogens(self):
        """
        Convert implicit hydrogens to explicit atoms and bonds.
        Returns the number of H atoms added.
        """
        added = 0
        n_atoms = len(self.atoms)
        atoms_to_process = list(range(n_atoms))

        for atom_idx in atoms_to_process:
            atom = self.atoms[atom_idx]
            num_h = atom.total_h

            for _ in range(num_h):
                h = Atom('H')
                h_idx = self.add_atom(h)
                self.add_bond(atom_idx, h_idx, BondType.SINGLE)
                added += 1

            atom.num_implicit_h = 0
            atom.num_explicit_h = 0

        return added

    # ─── Utilities ──────────────────────────────────────────────────

    def clone(self):
        """Create a deep copy of the molecule."""
        import copy
        mol = Molecule(self.name)
        mol.atoms = [copy.copy(a) for a in self.atoms]
        mol.bonds = [copy.copy(b) for b in self.bonds]
        mol._adjacency = {k: list(v) for k, v in self._adjacency.items()}
        mol.properties = dict(self.properties)
        return mol

    def __repr__(self):
        return (f"Molecule(name='{self.name}', atoms={len(self.atoms)}, "
                f"bonds={len(self.bonds)})")

    def __str__(self):
        return f"{self.name or 'Molecule'}: {self.molecular_formula()} ({len(self.atoms)} atoms)"
